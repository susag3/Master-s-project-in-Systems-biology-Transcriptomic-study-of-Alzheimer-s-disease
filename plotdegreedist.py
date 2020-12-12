import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import ScalarFormatter

X, Y = [], []
for line in open('degreedist_data.txt', 'r'):
    values = [float(s) for s in line.split()]
    X.append(values[0])
    Y.append(values[1])

fig, ax = plt.subplots()
ax.loglog(X,Y, 'o', color='black')
ax.axis ([1,70, 1, 1000])
plt.xlabel('Degree')
plt.ylabel('Number of nodes')

from matplotlib.ticker import ScalarFormatter
for axis in [ax.xaxis, ax.yaxis]:
    axis.set_major_formatter(ScalarFormatter())

b = np.linspace(1,69,69)
y=781.94*b**(-1.876)
plt.loglog(b,y,'r--')
plt.show()
