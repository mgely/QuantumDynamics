from pylab import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.close('figure 1')

# test data obtained from simulation
to_plot = np.loadtxt(open('data.dat', 'r'))

plot(to_plot.T)

show()