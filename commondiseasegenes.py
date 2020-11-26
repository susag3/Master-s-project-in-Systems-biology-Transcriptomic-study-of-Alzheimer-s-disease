#Finds the genes in common for two text files and print it to a new file
import pandas as pd 

df1 = pd.read_csv("totaldiseasegenes_sorted.txt", header = None) #input file
df1.columns = ['Gene']
print(df1)
df2 = pd.read_csv("CSD_genelist.txt", header = None)
df2.columns = ['Gene']
print(df2)
df = df1.merge(df2, on=["Gene"])
common_genes=sorted(df["Gene"].drop_duplicates(keep="first").values.tolist())
print(common_genes[:10], len(common_genes)) #prints first 10 genes and the number of common genes

#common_list =  list(set(list1).intersection(list2)) #intersection (common genes) of list1 and list2
#common_list.sort() #alphabetical order
#print(len(common_list)) #number of genes in common
#diff = set(list1).symmetric_difference(set(list2)) #genes in either, but not both
#print (diff, len(diff))

with open ("ADgenesCSD.txt", "w") as outfile: 
    for gene in common_genes:
        outfile.write(gene) #write gene to file
        outfile.write("\n") #newline puts the genes in one column