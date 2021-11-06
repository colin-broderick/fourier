import matplotlib.pyplot as plt
import numpy as np


data = np.genfromtxt("data.txt")
plt.scatter(data[:,0], data[:,1])
plt.axis("equal")
plt.show()
