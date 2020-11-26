#first in terminal: pip3 install python_louvain
#this script only works in ubuntu terminal

import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import community
from community import community_louvain
#plt.figure(figsize=(5,2.5)) #size of actual image

#read in the whole network as edgelist, separated by whitespace
CSD = nx.read_edgelist('CSD_edgelist2.txt', create_using=nx.Graph(), nodetype = str)
#CSD_edgelist.txt is tab-separated, didnt work

#First compute the best partition, returns dictionary with gene as key and module number as value
partition = community_louvain.best_partition(CSD)

mod_score = community_louvain.modularity(partition, CSD)
print(mod_score)

with open ('Louvain_new.txt', "w") as f:
    for key, value in partition.items():
        #print(key + '\t' + '\t'.join(str(value)), file=f)
        f.write(key + '\t')
        f.write(str(value))
        f.write('\n')
    





