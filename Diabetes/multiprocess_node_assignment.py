from collections import defaultdict
import pickle
from multiprocessing import Pool


import numpy as np
from Bio import SeqIO
import pandas as pd
import networkx as nx

aln_dict = SeqIO.to_dict(SeqIO.parse("/media/data/fileshare/T2D_seqs_CMaligned_with_Rfam.sto", "stockholm")) #nos séquences spécifiques à notre dataset

seq_to_ids = {str(v.seq).replace("-","").replace("U","T"): k for k, v in aln_dict.items()}

GG = pickle.load(open('BIG_GRAPH_with_edge_weights_and_node_int_ids.pickle','rb'))
print('loaded GG')



def build_graph_dataset(inputte):
    has_diabetes, start, end = inputte
    #G = copy.deepcopy(GG)
    abundance_seqs_T2D = list(df_T2D["#OTU ID"])
    sind = 0
    to_return = []
    for sample in list(df_T2D.columns)[start:end]:
        sind+=1
        abundance_dict = {}
        abundance_values = list(df_T2D[sample])
        for row in range(len(abundance_seqs_T2D)):
            if abundance_seqs_T2D[row] in seq_to_ids:
                abundance_dict[seq_to_ids[abundance_seqs_T2D[row]]] = abundance_values[row]
        max_ab = max(abundance_dict.values())
        if max_ab==0:
            continue
        #print("maximum abundance")
        for i in abundance_dict:
            abundance_dict[i] = 1000 * abundance_dict[i]/max_ab

        nodes_to_keep = [node for node, value in abundance_dict.items() if value > 0]
        G = GG.subgraph(nodes_to_keep).copy()
        node_index = {}
        for nind,node in enumerate(G.nodes(data=True)):
            node[1]['weight'] = abundance_dict[node[0]]
            node[1]["int_id"] = nind
            node_index[node[0]] = nind
            
        for e in G.edges(data=True):
            e[2]["source_int_id"] = node_index[e[0]]
            e[2]["dest_int_id"] = node_index[e[1]]
        to_return.append((G, has_diabetes))
    return start, to_return

import numpy as np
df_T2D = pd.read_csv('FILTERED_T2D_BIOM/AGP.data.biom.filtered.T2D.tsv', sep='\t')
steps = np.linspace(1,len(df_T2D.columns),21)
starts = [int(x) for x in steps[0:-1]]
ends = [int(x) for x in steps[1:]]
job_inputs = []
dataset = defaultdict(list)

for ind in range(len(starts)):
    job_inputs.append((True,starts[ind],ends[ind]))


with Pool(20) as p:
    for start, to_return in (p.imap_unordered(build_graph_dataset, job_inputs)):
        print(start)
        dataset[start].extend(to_return)

with open('marius_split_T2Dtrue_with_int.pickle', 'wb') as f:
    pickle.dump(dataset, f)
#===================================

df_T2D = pd.read_csv('FILTERED_noT2D_BIOM/AGP.data.biom.filtered.noT2D.tsv', sep='\t')
steps = np.linspace(1,len(df_T2D.columns),21)
starts = [int(x) for x in steps[0:-1]]
ends = [int(x) for x in steps[1:]]
job_inputs = []
dataset = defaultdict(list)

for ind in range(len(starts)):
    job_inputs.append((False,starts[ind],ends[ind]))


with Pool(20) as p:
    for start, to_return in (p.imap_unordered(build_graph_dataset, job_inputs)):
        print(start)
        dataset[start].extend(to_return)

with open('marius_split_T2Dfalse_with_int.pickle', 'wb') as f:
    pickle.dump(dataset, f)