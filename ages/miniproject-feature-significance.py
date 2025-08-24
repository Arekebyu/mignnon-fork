import numpy as np
import os
import pandas as pd
import pickle
from sklearn.cluster import AgglomerativeClustering
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.inspection import permutation_importance
import time

if os.path.exists("./ages/feature_importance.pkl"):
    with open("./ages/feature_importance.pkl", "rb") as file:
        importance = pickle.load(file)
    for age in importance:
        print("importance of predicting features for age", age, "in descending order is : \n", importance)
    quit()
# constants
LABEL_OF_INTEREST = 'age_cat'
PATH_TO_METADATA = "/media/data/BIAS/microbiome_GCN/AGES/AGP_ages.metadata.txt"
PATH_TO_OTU_DATA = "/media/data/BIAS/microbiome_GCN/AGES/AGP.data.biom.filtered.ages.tsv"
PATH_TO_DIST_DICT = "/media/data/BIAS/microbiome_GCN/AGES/MATRICES_AGES.pickle"

# loading the necessary files
print("Loading metadata")
df_metadata = pd.read_csv(PATH_TO_METADATA, sep='\t', low_memory=False).set_index("#SampleID")
print("Loading otu_count")
df_otu_count = pd.read_csv(PATH_TO_OTU_DATA,sep='\t').set_index("#OTU ID")
print("Loading distance dict")
with open(PATH_TO_DIST_DICT, "rb") as f:
    dist_dict = pickle.load(f)['TreeDistMatrix']
print("Done loading files")

# some more variables that will be used later on
features = df_otu_count.index
labels = df_metadata[LABEL_OF_INTEREST]
sample_ids = labels.keys()
num_features = len(features)

# convert distance into numpy matrix (why is it a map?)
dist_mat = np.zeros(shape=(num_features, num_features))
for i, list_of_dists in enumerate(dist_dict.values()):
    for j, dists in enumerate(list_of_dists.values()):
        dist_mat[i][j] = dists

# Apply hierarchical clustering onto features to group similar ones together 
# The metric is the distance matrix loaded in the previous step
# The number of clusters was determined (by eye) to be around 16
num_nodes = 16
hierarchical_clustering = AgglomerativeClustering(n_clusters=num_nodes, metric='precomputed', linkage='average')
hierarchical_clustering.fit(dist_mat)
hierarchical_clustering.labels_

clusters = [[] for _ in range(num_nodes)]

# send the ith term to their respective cluster
for index, cluster in enumerate(hierarchical_clustering.labels_):
    clusters[cluster].append(features[index])
new_order = []
# group clusters back together
for cluster in clusters:
    new_order = new_order + cluster

# reorder the features and take transpose
ordered_ages_df = df_otu_count.loc[new_order].T

# removes 15 or so data points corresponding to zeroes for all features
ordered_ages_df = ordered_ages_df.loc[~((ordered_ages_df == 0).all(axis='columns'))]
labels = labels[ordered_ages_df.index]

# normalizing 
sum = ordered_ages_df.sum(axis='columns')
ages_percentage_df = (ordered_ages_df.div(sum, axis='index'))

# convert into ids
ENCODINGS = {x : i for i, x in enumerate(sorted(ages_percentage_df.keys()))}
temp = ages_percentage_df.rename(mapper=ENCODINGS, axis='columns')

# find most significant predictors
selected_features = set()
for i, age in enumerate(labels.unique()):
    ''' Idea behind the criterion below is that we want to find features where 
        there is a large difference between an age's average value for that 
        feature and other ages' features.
        Standard deviation is used because mean may yield small values despite 
        having large differences.
    '''
    own_deviation = temp.loc[labels[labels == age].index].std()
    other_deviation = temp.loc[labels[labels != age].index].std()
    evaluation = (own_deviation - other_deviation) ** 2
    selected_features |= set(list(evaluation.sort_values(ascending=False).iloc[0:700].index))

# selects those best features according to our values
selected_features = list(selected_features)
# selected_features = list(set(temp.columns).difference(selected_features))
'''Note that if you want to remove the bottom worst features, set ascending
    to True on line 77, comment out line 80 and uncomment line 81'''

selected_ages_df = temp.loc[:, selected_features]

# turns age into numerical value 
AGE_ENCODINGS = {x : i for i, x in enumerate(sorted(labels.unique()))}
y = labels.map(AGE_ENCODINGS).reindex(selected_ages_df.index)

# prepare data to train random forest
X = selected_ages_df.to_numpy()
y = y.to_numpy()

# test/train split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.10 , random_state=11)

print("Training model...")
random_forest = RandomForestClassifier(random_state=11, n_estimators=200)
random_forest.fit(X_train, y_train)
print("Done!")

# begin analysis
results = {}
print("Beginning permutation importance...")
for cls in AGE_ENCODINGS.values():
    begin = time.time()
    mask = y_train == cls
    imp = permutation_importance(random_forest, X_train[mask], y_train[mask], n_repeats=1)
    results[cls] = imp.importances_mean
    print(cls, "finished in ", time.time() - begin)

inv_encodings = {v: k for k, v in ENCODINGS.items()}
importance = {}
for age in results:
    temp = results[age]
    temp = {inv_encodings[i]: imp for i, imp in enumerate(temp)}
    temp = dict(sorted(temp.items(), key=lambda item: item[1], reverse=True))
    importance[sorted(labels.unique())[age]] = temp
    print("importance of predicting features for age", age, "in descending order is : \n", temp)
with open("./ages/featureImportance.pkl", "wb") as file:
    pickle.dump(importance, file, pickle.HIGHEST_PROTOCOL)
