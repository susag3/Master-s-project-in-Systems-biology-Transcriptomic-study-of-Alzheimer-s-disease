#This code performs community detection of the CSD network using the Louvain method, and calculated the modularity score.

#Code was written in Python 3.7.4. 
#First need to pip3-install python_louvain.

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import community
from community import community_louvain

#Read in the whole network as edgelist, separated by whitespace
CSD = nx.read_edgelist('CSD_edgelist2.txt', create_using=nx.Graph(), nodetype = str)

#Compute the best community partition, returns dictionary with gene as key and module number as value
partition = community_louvain.best_partition(CSD)

#Calculate modularity score
mod_score = community_louvain.modularity(partition, CSD)
print(mod_score)

#Write results to file; gene symbol as first column and corresponding module number as second column
with open ('Louvain_new.txt', "w") as f:
    for key, value in partition.items():
        #print(key + '\t' + '\t'.join(str(value)), file=f)
        f.write(key + '\t')
        f.write(str(value))
        f.write('\n')
    





