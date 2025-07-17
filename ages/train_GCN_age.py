from skbio import TreeNode
import pandas as pd
from Bio import SeqIO
import networkx as nx
import pickle
import copy
##### pytorche
import torch
from collections import defaultdict
from typing import Any, Dict, Iterable, List, Literal, Optional, Tuple, Union
import os

os.environ['CUDA_LAUNCH_BLOCKING']="1"
os.environ['TORCH_USE_CUDA_DSA'] = "1"
os.environ['TORCH'] = torch.__version__
print(torch.__version__)
from torch import Tensor
from torch.utils.dlpack import from_dlpack, to_dlpack

import torch_geometric
from torch_geometric.utils.num_nodes import maybe_num_nodes
from torch_geometric.utils.convert import from_networkx
from torch.nn import Linear, Sequential, ReLU, Softmax, Embedding, Sigmoid, LeakyReLU
from torch_geometric.nn.norm import BatchNorm
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool, global_add_pool, global_max_pool, SAGPooling, RGCNConv, GINEConv
import torch.nn as nn

from torch_geometric.loader import DataLoader
import matplotlib.pyplot as plt
from tqdm import tqdm
import math
import json

from torch_geometric.data import Data, Dataset
from sklearn.preprocessing import LabelEncoder
import torch
from sklearn.model_selection import train_test_split



class EdgeCentricRGCN(torch.nn.Module):
    def __init__(self, num_edge_types, hidden_dim=256, output_dim=5):
        super().__init__()
        # Encodage des types d'arêtes
        self.hidden_dim = hidden_dim
        # Projection des caractéristiques initiales des noeuds (tous à 1)
        #self.node_encoder = Linear(1, hidden_dim)

        # Prendre en compte les arretes dans les deux sens
        self.rgcn_in = RGCNConv(hidden_dim, hidden_dim, num_relations=num_edge_types)
        self.rgcn_out = RGCNConv(hidden_dim, hidden_dim, num_relations=num_edge_types)

        # Couches de convolution GINE avec traitement des arêtes
        self.conv1 = GINEConv(
            nn=Sequential(
                Linear(hidden_dim, hidden_dim),
                ReLU(),
                BatchNorm(hidden_dim),
                Linear(hidden_dim, hidden_dim )
                
            ))

        self.conv2 = GINEConv(
            nn=Sequential(
                Linear(hidden_dim, hidden_dim),
                ReLU(),
                BatchNorm(hidden_dim),
                Linear(hidden_dim, hidden_dim )
                

            ))

        self.mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            ReLU(),
            BatchNorm(hidden_dim),
            nn.Linear(hidden_dim, output_dim),  # 1 seule valeur prédite
            nn.LogSoftmax()
        )

    def forward(self, data):
        edge_index = data.edge_index

        edge_weight = data.edge_weight.to(data.edge_index.device)
        node_weight = data.node_weight.to(data.edge_index.device)

        node_encoder = Linear(1,self.hidden_dim).to(data.edge_index.device)

        num_nodes = len(node_weight)
        num_edges = len(edge_weight)
        node_weight= torch.unsqueeze(node_weight,0)
        node_weight = torch.transpose(node_weight,0,1)

        x=node_encoder(node_weight).to(data.edge_index.device)

        edge_weight = torch.unsqueeze(edge_weight,0).to(data.edge_index.device)
        edge_weight = torch.transpose(edge_weight,0,1).to(data.edge_index.device)

        edge_encoder = Linear(1, self.hidden_dim).to(data.edge_index.device)
        edge_weight = edge_encoder(edge_weight).to(data.edge_index.device)

        x = F.relu(self.conv1(x, data.edge_index, edge_weight)).to(data.edge_index.device)
        x = F.relu(self.conv2(x, data.edge_index, edge_weight)).to(data.edge_index.device)

        x = global_add_pool(x, data.batch).to(data.edge_index.device)
        return self.mlp(x).squeeze(-1)

EDGE_TYPES = [
    "phylo_dist",
    "seq_dist"
]

class AGPGraphDataset(Dataset):
    def __init__(self, graphs, labels, normalise=True):
        super().__init__()
        self.graphs = graphs
        self.labels = labels
        #self.normalise = normalise

        # Création d'une correspondance ordonnée
        self.ordered_keys = list(self.labels)

        # Créer un mapping des edge types vers des entiers
        self.edge_type_encoder = LabelEncoder()

        # LabelEncoder pour les types d'arêtes
        self.edge_type_encoder.fit(EDGE_TYPES)


    def len(self):
        return len(self.graphs)

    def get(self, idx):
        ind = idx = int(idx) #index of the graph
        graph_id = idx  # Clé réelle du JSON
        this_graph = self.graphs[idx]
        this_label = self.labels[idx]
        
        
        # Conversion des noeuds
        node_features = [n[1]["int_id"] for n in this_graph.nodes(data=True)]        
        unnormalized_node_weights = [n[1]["weight"] for n in this_graph.nodes(data=True)]
        three_quarters_weight, max_weight = np.quantile(unnormalized_node_weights, [0, 0.25,0.5,0.75,1])[3:5] 
        node_weight = [x/max_weight if x<three_quarters_weight else 1.0 for x in unnormalized_node_weights]
        x = torch.tensor(node_features, dtype=torch.float32)
        node_weight = torch.tensor(node_weight, dtype=torch.float32)

        edge_coo = [[e[2]["source_int_id"] for e in self.graphs[ind].edges(data=True)],[e[2]["dest_int_id"] for e in self.graphs[ind].edges(data=True)]]
        edge_index = torch.tensor(edge_coo, dtype=torch.long)

        edge_type = [e[2]["type"] for e in self.graphs[ind].edges(data=True)]
        edge_type = self.edge_type_encoder.transform(edge_type)
        edge_type = torch.tensor(edge_type, dtype=torch.long)
        edge_three_quarter_weight = 396.536581
        edge_max_weight = 865
        edge_weight = [e[2]["phylo_weight"]/edge_max_weight if e[2]["phylo_weight"]<edge_three_quarter_weight else 1.0 for e in self.graphs[ind].edges(data=True)]
        #edge_weight = [1-e for e in edge_weight] 
        edge_weight = torch.tensor(edge_weight, dtype=torch.float32)
        return Data(x=x, edge_index=edge_index, edge_type=edge_type, edge_weight = edge_weight, y=torch.tensor(this_label, dtype=torch.long), node_weight=node_weight)

    def denormalize_prediction(self, y_norm):
        """Convertit la sortie normalisée vers la vraie valeur cible (exp du log inverse)."""
        return y_norm  * self.target_max

import numpy as np
def train(train_graphs, test_graphs, train_labels,test_labels):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print("device", device)
    train_dataset = AGPGraphDataset(train_graphs, train_labels, normalise=False)
    test_dataset = AGPGraphDataset(test_graphs, test_labels, normalise=False)
    #print(train_dataset)

    sample = train_dataset.get(0)
    print("Sample x:", sample.node_weight)
    print("len x", len(sample.x))
    print("Sample y:", sample.y)
    print("Sample edge_index:", sample.edge_index.size())
    print("Sample edge_type:", sample.edge_type.size())
    print("sample edge weight", sample.edge_weight)
    print("sample node weight", sample.node_weight.size())
    
    bs = 256
    # DataLoaders
    train_loader = DataLoader(train_dataset, batch_size=bs, shuffle=True,
                              num_workers=54, pin_memory=True)
    test_loader = DataLoader(test_dataset, batch_size=bs, shuffle=True,
                             num_workers=54, pin_memory=True)

    model = EdgeCentricRGCN(num_edge_types=len(EDGE_TYPES)).to(device)
    model = model.to(torch.float32)
    optimizer = torch.optim.Adam(model.parameters(), lr=1e-2)
    #optimizer = torch.optim.SGD(model.parameters(), lr = 1e-3)
    loss_fn = nn.NLLLoss()

    # Boucle d'entraînement
    for epoch in range(20):
        model.train()
        total_loss = 0
        trainAccuracy = 0
        totalTrainAccuracy = 0
        valCorrect = 0
        for batch in train_loader:
            #print(batch)
            batch = batch.to(device)
            optimizer.zero_grad()
            
            with torch.cuda.amp.autocast(dtype=torch.float32):
                pred = model(batch)
                pred_real = pred
                y_real = batch.y
                
                #trainAccuracyo = sum([1 if y_real[x] == pred_real[x].argmax().cpu().numpy()[0] else 0 for x in range(len(pred_real))])
                #print("output",[x.argmax().numpy()[0] for x in pred_real])
                #print("train accuracy result",trainAccuracyo)
                #print("denominator", y_real.size(0))
                #trainAccuracy = trainAccuracyo / y_real.size(0)
                #totalTrainAccuracy += trainAccuracy
            loss = loss_fn(pred, batch.y)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            total_loss += loss.item()
        #print("loss", loss.item())
        print("total loss", loss)
        model.eval()
        #avgAccuracy = totalTrainAccuracy / len(train_loader)
        #print('train accuracy', avgAccuracy)
        #total_mae = 0
        #total = 0
        #all_preds = []
        #all_trues = []

        #with torch.no_grad():
        #    testAccuracy = 0
        #    totalTestAccuracy = 0
        #    correct = 0
        #    for batch in test_loader:
        #        batch = batch.to(device)
        #        pred = model(batch)

                #pred_real = pred
                #y_real = batch.y
                #print("y predictions", pred)
                #print("y real", y_real)
                #all_preds.append(pred_real.cpu())
                #all_trues.append(y_real.cpu())
                #this_round_loss =  loss_fn(pred, batch.y).tolist()
                #total_mae += this_round_loss
                #print("epoch loss and total loss:",this_round_loss, total_mae)
                #correct +=  sum([1 if y_real[x] == pred_real[x].argmax().cpu().numpy()[0] else 0 for x in range(len(pred_real))])
                #total += batch.num_graphs



        #all_preds = torch.cat(all_preds)
        #all_trues = torch.cat(all_trues)
        #t_acc = correct/total
        #print('test accuracy', t_acc)
        #print(f"Epoch {epoch} | Train Loss (CEL): {total_loss / (bs * len(train_loader)):.4f} | "
        #      f"Test CEL: {total_mae / total:.4f} | " f"Train accuracy: {avgAccuracy:.4f}| " f"Test accuracy: {t_acc:.4f}")
        #print(predictions)
        # Plot prédictions vs vraies valeurs (chaque 10 epochs)
        #if epoch % 10 == 0:
        #    plt.figure(figsize=(6, 6))
        #    plt.scatter(all_trues, all_preds, alpha=0.6, edgecolors='k')
        #    plt.plot([min(all_trues), max(all_trues)], [min(all_trues), max(all_trues)], 'r--')  # y = x
        #    plt.xlabel("Vraies valeurs")
        #    plt.ylabel("Prédictions")
        #    plt.title(f"Époque {epoch} - Prédictions vs Réalité")
        #    plt.grid(True)
        #    plt.tight_layout()
        #    plt.show()
        #print("mean of predictions", np.mean(all_preds.numpy()), "mean of labels", np.mean(all_trues.numpy()))

    return model



with open('marius_split_AGES_with_int.pickle', 'rb') as f:
    T2D_dataset = pickle.load(f)

GRAPHS = []
LABELS = []
for i in T2D_dataset:
    print(i)
    for j in T2D_dataset[i]:
        print(j)
        GRAPHS.append(j[0])
        LABELS.append(j[1])



#GRAPHS,LABELS = pickle.load(open("T2D_BIG_DATASET_FULL","rb"))

X_train, X_test, y_train, y_test = train_test_split(GRAPHS, LABELS, test_size=0.3, random_state=133, shuffle=True)

model = train(X_train, X_test, y_train, y_test)
checkpoint = {
    "model_state": model.state_dict()
}
torch.save(checkpoint,"./best_model_age.pth")