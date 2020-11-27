#Takes input from 4 files, 'UseableCValues.txt', 'UseableSValues.txt', 'UseableDValues.txt', 'AllValues.txt' (output from findCSD.py) and generates 4 networks - one for each of the C/S/D-type interactions, as well as an aggregate network containing all 3, and a network file including all the rho values (from AllValues.txt)

#Outputs selected pairs, along with metric type and value, to SelectedSNodes.txt, SelectedCNodes.txt, SelectedDNodes.txt, and CSDSelection.txt and CSDSelectionDetailed.txt


##############################################
#Parameters to vary are selSize and noSels

#selSize indicates approximate proportion of selected nodes, being equal to the 1/(desired p-value). Increasing selSize selects fewer edges. 
#noSels may be increased as desireable, at the expense of longer running time. May also be decreased, at the expense of accuracy


import math
import numpy
import random

selSize = 200000
noSels = 10000


cValueFile = 'UseableCValues.txt'
sValueFile = 'UseableSValues.txt'
dValueFile = 'UseableDValues.txt'

f = open(cValueFile)

valueList = []

for line in f:
    valueList.append(float(line))


cutoffAtSel = []
for i in range(noSels):

    selection = [valueList[i] for i in random.sample(xrange(0, len(valueList)), selSize)]
    cutoffAtSel.append(max(selection))

cCutoff = numpy.mean(cutoffAtSel)

f.close()



f = open(sValueFile)

valueList = []

for line in f:
    valueList.append(float(line))


cutoffAtSel = []
for i in range(noSels):

    selection = [valueList[i] for i in random.sample(xrange(0, len(valueList)), selSize)]
    cutoffAtSel.append(max(selection))

sCutoff = numpy.mean(cutoffAtSel)
f.close()


f = open(dValueFile)

valueList = []

for line in f:
    valueList.append(float(line))


cutoffAtSel = []
for i in range(noSels):

    selection = [valueList[i] for i in random.sample(xrange(0, len(valueList)), selSize)]
    cutoffAtSel.append(max(selection))

dCutoff = numpy.mean(cutoffAtSel)
f.close()

cNetF = open('CNetworkold.txt', 'w')
sNetF = open('SNetworkold.txt', 'w')
dNetF = open('DNetworkold.txt', 'w')
csdNetF = open('CSDSelectionold.txt', 'w')
detailedF = open('CSDSelectionDetailedold.txt', 'w')

f = open('AllValues.txt')

f.readline()

for line in f:

    splitLine = line.rstrip().split('\t')
    if float(splitLine[4]) > cCutoff:
        print >> cNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[4])+'\t'+'C'
        print >> csdNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[4])+'\t'+'C'
        print >> detailedF, line.rstrip()+'\t'+'C'

    if float(splitLine[5]) > sCutoff:
        print >> sNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[5])+'\t'+'S'
        print >> csdNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[5])+'\t'+'S'
        print >> detailedF, line.rstrip()+'\t'+'S' 

    if float(splitLine[6]) > dCutoff:
        print >> dNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[6])+'\t'+'D'
        print >> csdNetF, str(splitLine[0])+'\t'+str(splitLine[1])+'\t'+str(splitLine[6])+'\t'+'D'
        print >> detailedF, line.rstrip()+'\t'+'D' 



cNetF.close()
sNetF.close()
dNetF.close()
detailedF.close()
f.close()
