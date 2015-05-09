from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np


dimensions = np.loadtxt(open('dimensions.dat', 'r'))
N_x = int(dimensions[0])
N_y = int(dimensions[1])

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(0, N_x, 1)
Y = np.arange(0, N_y, 1)
X, Y = np.meshgrid(X, Y)
psi = np.loadtxt(open('psi.dat', 'r'))

surf = ax.plot_surface(X, Y, psi, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
ax.set_zlim(0, 1.01)

ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()