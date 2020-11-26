#This code calculates node homogeneity (H) for each node in the CSD network and generates a box plot of H as a function of node degree.

from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

#Create preliminary dictionary and list for further use
homo = defaultdict(lambda:{'C': 0, 'S': 0, 'D': 0})
homolist = [[] for _ in range(70)] #range must be equal to or larger than the maximum node degree in network

#Read in CSD network edge list, including link scores, correlation values and last column is link type (C,S or D)
with open('CSDSelectionDetailedold.txt', newline='') as f:
    for line in f:
        line = line.rstrip().split('\t')
        homo[line[0]][line[-1]] += 1 #counter of number of links of a specific type

for gene, links in homo.items():
    degree = sum(links.values()) #degree is total number of links
    homoval = round(sum(link**2 for link in links.values())/(degree**2), 2) #Homogeneity equation, round up to 2 decimals
    homo[gene].update({'Homo': homoval, 'Degree': degree})
    homolist[degree] += [homoval]

#Print gene symbol and its values in dictionary (# of C-,S-,D- links, homogeneity score and node degree)
for gene, vals in homo.items():
    print(gene, vals)

#Box plot generation
fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
ax.set_xlabel('Degree')
ax.set_ylabel('Homogeneity')
ax.boxplot(homolist, showmeans=True) #show mean values
plt.show()

