#This script finds the genes in common for two input text files and prints it to a new file. 
#Used to find genes previously associated with AD in the CSD network.

import pandas as pd 

#List of AD-affiliated genes for AD from MalaCards
df1 = pd.read_csv("totaldiseasegenes_sorted.txt", header = None) 
df1.columns = ['Gene']

#list of all genes in the CSD network
df2 = pd.read_csv("CSD_genelist.txt", header = None)
df2.columns = ['Gene']

#Merge dataframes to compare for unique genes
df = df1.merge(df2, on=["Gene"])
common_genes=sorted(df["Gene"].drop_duplicates(keep="first").values.tolist())
print(common_genes[:10], len(common_genes)) #prints first 10 genes and the number of common genes

#Write the list of AD-related genes in CSD network to a new file
with open ("ADgenesCSD.txt", "w") as outfile: 
    for gene in common_genes:
        outfile.write(gene) #write gene to file
        outfile.write("\n") #newline puts the genes in one column
