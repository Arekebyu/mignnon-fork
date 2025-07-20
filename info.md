In this file, you will find the data for the microbiome status prediction project, as well as comments on the code.
Pour le français, vous pouvez référer aux diapos présentées au meeting de labo.

To build a model for label prediction from microbiome profile, you need:
- a .biom file, which contains the abundance data (each row is one sequence, each column is a sample)
- a .tsv metadata file, which contains a typical csv file, where each row is one sample, each column is a variable
- a .nwk reference phylogeny file, which gives us evolutionary information about our data

option: a taxonomy file, estimated from the phylogeny, that gives us a list of species.
option: a sequences file (often fasta) that contains the sequences in the BIOM file.

# CODE
The code is in python 3.11
So the code here is all over the place. I think I have met each of you to explain what I am doing with the code so that should help make it straightforward.

There are two folders with similar code: Diabetes and ages, corresponding to two different projects.
The code is divided in basically 3 steps, each with a notebook.

1. extracting data from the american gut project database. Here we grab the data we need, download it, subset it and clean it according to our needs.
2. associating the data with phylogeny and obtaining distances between the sequences so that we can build the graph
3. building graphs and using them to train models.

The first step is not actually done through a notebook, it's done in command line and requires QIIME2. It's already done for you since it's what we use to get the metadata and the BIOM files we need.
The second step is the data engineering, it's ultimately very arbitrary, I did it one way but many other ways are possible.All that matters is that we take as input a BIOM table and a dirty metadata file, and output a pandas dataframe, clear labels, and a distance matrix of some kind.
The last step is the actual preparation of the graph and training on the model. The notebook shows the overall process, but I implemented the two key steps as python files to allow for multiprocessing and execution with nohup.

# DATA 
### Project Diabetes - small dataset (a few hundred sequences)
Here are the files for the diabetes project. All the paths are on frayere

**diabetes positive**
BIOM file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_T2D_BIOM/AGP.data.biom.filtered.T2D.tsv
metadata file : /media/data/BIAS/microbiome_GCN/AGP/AGP_T2D.metadata.txt
sequences file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_T2D_BIOM/AGP_seqs.fa 

**diabetes negative**
BIOM file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_noT2D_BIOM/AGP.data.biom.filtered.noT2D.tsv
metadata file : /media/data/BIAS/microbiome_GCN/AGP/AGP_T2Dcontrol.metadata.txt
sequences file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_noT2D_BIOM/AGP_seqs.fa 

Reference phylogeny: /media/data/BIAS/microbiome_GCN/2024.09.phylogeny.asv.nwk 
Reference taxonomy: /media/data/BIAS/microbiome_GCN/2024.09.taxonomy.asv.tsv

### Project Diabetes - BIG dataset (as many sequences as possible)
**diabetes positive**
BIOM file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_T2D_BIOM/AGP.data.biom.filtered.T2D.tsv
metadata file : /media/data/BIAS/microbiome_GCN/AGP/FILTERED_T2D_BIOM/AGP_T2D.metadata.txt
sequences file: /media/data/BIAS/microbiome_GCN/AGP/FILTERED_T2D_BIOM/AGP_seqs.fa 

**diabetes negative**
BIOM file: /media/data/BIAS/microbiome_GCN/BIGDIABETES/FILTERED_noT2D_BIOM/AGP.data.biom.filtered.noT2D.tsv
metadata file : /media/data/BIAS/microbiome_GCN/BIGDIABETES/FILTERED_noT2D_BIOM/AGP_T2Dcontrol.metadata.txt
sequences file: /media/data/BIAS/microbiome_GCN/BIGDIABETES/FILTERED_noT2D_BIOM/AGP_seqs.fa 

(these are same as previous)
Reference phylogeny: /media/data/BIAS/microbiome_GCN/2024.09.phylogeny.asv.nwk 
Reference taxonomy: /media/data/BIAS/microbiome_GCN/2024.09.taxonomy.asv.tsv

### Project AGE OF POOP dataset
BIOM file: /media/data/BIAS/microbiome_GCN/AGES/AGP.data.biom.filtered.ages.tsv
metadata file : /media/data/BIAS/microbiome_GCN/AGES/AGP_ages.metadata.txt
sequences file: /media/data/BIAS/microbiome_GCN/AGES/AGP_ages_seqs.fa

Precomputed phylogeny distances file: /media/data/BIAS/microbiome_GCN/AGES/MATRICES_ages.pickle

(these are same as previous)
Reference phylogeny: /media/data/BIAS/microbiome_GCN/2024.09.phylogeny.asv.nwk 
Reference taxonomy: /media/data/BIAS/microbiome_GCN/2024.09.taxonomy.asv.tsv
