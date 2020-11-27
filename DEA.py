'''
DEA: Differential Expression Analysis
This script calculates mean expression of all genes in control and case patients separately.
Then calculates log2 fold change (difference in mean), raw p-value, T-statistic
and adjusted p-value (FDR by Benjamini-Hochberg)
User must change input and output files (lines 15, 20, 53), corresponding to desired tissue calculation
'''
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind #T-test
import matplotlib.pyplot as plt
from math import log10

#Expression for only control
data2 = pd.read_csv("ADcontrolold.txt", sep = "\t", header = None, skiprows=[0], index_col=0) #expression data control, skip first row
df2 = pd.DataFrame(data2)
df2['mean'] = df2.mean(axis=1) #mean expression of each row (gene)

#Expression for only AD patients
data3 = pd.read_csv("ADexp.txt", sep = "\t", header = None, skiprows=[0], index_col=0) #expression data sick, skip first row
df3 = pd.DataFrame(data3)
df3['mean'] = df3.mean(axis=1) #mean expression of each row (gene)

#Calculate fold change from control to AD 
log2FC = (df3['mean']-df2['mean']) #fold change is just difference in mean, because values already log2
newdf = pd.DataFrame(index=df2.index.copy()) #new dataframe with gene name as index (first column)
newdf.index.names = ['Gene name'] #Gene name as column name
newdf['Mean control exp.'] = df2['mean'] #add mean expression for control as new column, round up to 2 decimals
newdf['Mean AD exp.'] = df3['mean'] #add mean expression for sick as new column
newdf['log2FC'] = log2FC #add log2FC as new column to dataframe

##Perform multiple t-test, row-by-row (for every gene in both dataframes)
df_m = pd.merge(df2, df3, left_index=True, right_index=True)
T_stat, p_vals = ttest_ind(df_m.iloc[:, df2.shape[1]:-1], df_m.iloc[:, :df2.shape[1]-1], axis=1) #compare y with x
print(f'T-statistic= {T_stat}, Raw p-value = {p_vals}')
newdf['Raw p-value'] = p_vals #add p-values to new column in dataframe
newdf['T-statistic'] = np.round(T_stat, decimals = 2) #add T-statistic to new column in dataframe

def fdr(p_vals): #function that adjusts p-vals, returns Benjamini-Hochberg adjusted P-value
    from scipy.stats import rankdata
    ranked_p_values = rankdata(p_vals)
    fdr = p_vals * (len(p_vals) / ranked_p_values)
    fdr[fdr > 1] = 1

    return fdr

newdf['FDR'] = fdr(p_vals) #perform function and add adjusted P-value
newdf.sort_values(by=['log2FC'], ascending = False, inplace = True) #sort from highest to lowest log2FC

#Write dataframe to file
newdf.to_csv('diffstats_allregionsold.txt', index=True, sep = '\t') 
