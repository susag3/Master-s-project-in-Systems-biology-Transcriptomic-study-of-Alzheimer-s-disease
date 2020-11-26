from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

homo = defaultdict(lambda:{'C': 0, 'S': 0, 'D': 0})
homolist = [[] for _ in range(70)]
with open('CSDSelectionDetailedold.txt', newline='') as f:
    for line in f:
        line = line.rstrip().split('\t')
        homo[line[0]][line[-1]] += 1

for gene, links in homo.items():
    degree = sum(links.values()) #degree is total number of links
    homoval = round(sum(link**2 for link in links.values())/(degree**2), 2) #Homogeneity equation, round up to 2 decimals
    homo[gene].update({'Homo': homoval, 'Degree': degree})
    homolist[degree] += [homoval]

for gene, vals in homo.items():
    if gene == 'ENPP2' or gene == 'CADPS' or gene == 'MDH1' or gene == 'VSNL1' or gene == 'LCAT':
    #gene == 'AQR' or gene == 'AL158206.1' or gene == 'HPRT1' or gene == 'GTF2I' or gene == 'YWHAH' or gene == 'TOM1L2' or gene == 'TMEM178A' or gene == 'LCAT' or gene == 'PLTP' or gene == 'NAPB' or gene == 'GOT1':
        print(gene, vals)

fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)
ax.set_xlabel('Degree')
ax.set_ylabel('Homogeneity')
ax.boxplot(homolist, showmeans=True)
#plt.show()

