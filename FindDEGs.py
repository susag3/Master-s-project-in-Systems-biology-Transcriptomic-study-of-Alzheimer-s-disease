#This script finds the significantly upregulated and downregulated genes from the DEA, and writes to individual files.

#Opens the outputfile from DEA.py, including all expression changes 
with open('diffstats_allregionsold_volcano.txt') as f:
    upDEGs = []
    downDEGs = []
    firstline = f.readline()
    for line in f:
        splitline = line.rstrip().split('\t')
        if float(splitline[6])<0.05: #genes with FDR < 0.05 are significant
            if float(splitline[3])>=0.2: #log2FC > 0.2 are upregulated
                upDEGs.append(splitline[0]) #add gene symbol
            elif float(splitline[3])<=(-0.2): #log2FC < -0.2 are downregulated
                downDEGs.append(splitline[0]) #add gene symbol
print(len(upDEGs), upDEGs) #Number of up-regulated genes and the list of gene symbols
print(len(downDEGs), downDEGs) #Number of down-regulated genes and the list of gene symbols

#Create a file where the gene symbol of up-DEGs are on invidual lines
with open('up-DEGs.txt', 'w') as file1:
    for i, gene in enumerate(upDEGs):
        file1.write(gene)
        file1.write('\n')

#Create a file where the gene symbol of down-DEGs are on invidual lines
with open('down-DEGs.txt', 'w') as file2:
    for i, gene in enumerate(downDEGs):
        file2.write(gene)
        file2.write('\n')

