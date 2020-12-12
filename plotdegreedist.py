import matplotlib.pyplot as plt
import numpy as np

X, Y = [], []
#open file from Cytoscape where first/second column is degree/number of nodes
for line in open('degreedist_data.txt', 'r'): 
    values = [float(s) for s in line.split()]
    X.append(values[0])
    Y.append(values[1])

fig, ax = plt.subplots()
ax.loglog(X,Y, 'o', color='black') #log-log scale
ax.axis ([1,70, 1, 1000])
plt.xlabel('Degree')
plt.ylabel('Number of nodes')

#Fit power law, equation calcuated in Cytoscape
b = np.linspace(1,69,69)
y=781.94*b**(-1.876)
plt.loglog(b,y,'r--') #red dotted line

plt.show()
