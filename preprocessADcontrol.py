import sys
import os
import statistics


with open("AffyAD2.txt") as f:
    #AffyAD2 has gene name as first column and AFFY U133 Plus 2.0 probe as second column

    #Make dictionary of probeID with gene name
    probeToGene = {}
    badprobes = [] #probes that link to more than one gene, want to filter them out (done in line 19)

    f.readline()
    #make the gene dictionary, probeID first column (key) and gene name in second (value)
    for line in f:
        geneName, probeNumber = line.rstrip().split(",")
        if (geneName and probeNumber) and (probeNumber not in badprobes) and geneName!= "0": #if both contains something, not bad probe, not zero-gene
            if probeNumber in probeToGene and probeToGene[probeNumber] != geneName: #probenumber exists and 
                badprobes.append(probeNumber)
                del probeToGene[probeNumber] #delete the probe
            else:
                probeToGene[probeNumber] = geneName #add the key-value (probe-gene) pair to the dict
    
    probesLookup = {} #dictionary where key is gene and value(s) are probe(s)
    for probeNumber in probeToGene:
        if probeToGene[probeNumber] in probesLookup:
            probesLookup[probeToGene[probeNumber]].append(probeNumber) #add probe to the probe list for gene if it exists in gene dictionary
        else:
            probesLookup[probeToGene[probeNumber]] = [probeNumber] #if not exists, add key-value pair to dict
    #print(list((gene, len(probes), probes) for gene, probes in probesLookup.items() if len(probes)>1)) #print the genes which have more (good) probes

##Sample extraction to columns## 
with open("E-GEOD-48350-experiment-design", "r", encoding="utf8") as d:
    firstline = d.readline()
    old_ind = [] #all >= 60 years
    HC = [] #hippocampus individuals
    EC = [] #entorhinal cortex individuals
    SFG = [] #superior frontal gyrus individuals
    PCG = [] #postcentral gyrus individuals
    for line in d: #line-by-line in file
        splitLine = line.rstrip().split("\t")
        if splitLine[2].startswith(("6", "7", "8", "9")): #if column three in file starts with either 6,7,8,9 (it means age>=60)
            old_ind.append(splitLine[0]) #add samplename (column 1) to list of old individuals
            #print(splitLine[2], splitLine[0]) #prints age column and samplename
            if splitLine[12] == 'hippocampus': #hippo samples over 60
                HC.append(splitLine[0])
            elif "entorhinal cortex" in splitLine[12]: #EC samples over 60
                EC.append(splitLine[0])
            elif "superior frontal gyrus" in splitLine[12]: #SFG samples over 60
                SFG.append(splitLine[0])
            elif "postcentral gyrus" in splitLine[12]: #PCG samples over 60
                PCG.append(splitLine[0])
    AD = [s for s in old_ind if "GSM117" in s] #Alzheimers over 60 (all of the AD samples)
    Contr = [s for s in old_ind if ("GSM300" in s) or ("GSM318840" in s) or ("GSM350078" in s)] #control over 60


with open("E-GEOD-48350", "r", encoding="utf8") as f: #read from this file, Alzheimers disease set
    #Find columns based on the lists above, to extract correct expression columns
    line = f.readline() #use the header
    splitLine = line.rstrip().split("\t")
    SampleCol = [] 
    for i, col in enumerate(splitLine): #go through columns in header
        if (col in HC) and (col in Contr): #hippo over 60 and case/control
            SampleCol.append(i) #column number(s) with desired individual(s)

    probeExpTable = {} #make dict of probes with expression level for each patient (equal to input-file)
    for line in f:
        splitLine = line.rstrip().split("\t")
        if splitLine[2] in probeToGene: #probes that are in AffyAD2.txt, third column
            probeExpTable[splitLine[2]] = splitLine
            if splitLine[1] != probeToGene[splitLine[2]]: #check
                print (splitLine[1:2], probeToGene[splitLine[2]])

    geneExpTable = {}        
    for gene in probesLookup:
        geneExpTable[gene] = [] #list of gene expressions per gene
        for i in SampleCol: #calculate for each sample(patient with diagnose) in dataset
            allProbesExpVals = [] #exp vals of all the probes for a gene
            for probe in probesLookup[gene]:
                allProbesExpVals.append(float(probeExpTable[probe][i])) #add gene expressions for all probes of a gene
                #print(probeExpTable[probe])           
            #print(allProbesExpVals)
            geneExpTable[gene].append(str(statistics.mean(allProbesExpVals))) #average of expression values

with open("HC_Contr.txt", "w", encoding="utf8") as of: #write to this file
    print('Gene names Expression', file=of) #create header of the file
    for gene, expValues in geneExpTable.items():
        print(gene + '\t' + '\t'.join(expValues), file=of) #prints gene and expression values to the file with '\n'

print(len(SampleCol), SampleCol, len(geneExpTable)) #number of columns and genes
print (os.access('HC_Contr.txt', os.R_OK)) #True if file exists


   