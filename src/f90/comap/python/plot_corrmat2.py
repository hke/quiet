from setup_matplotlib import *
import matplotlib.cm as cm
import numpy as N
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable


# Create the plot
width=8.8
fig = plt.figure(figsize=(1.4*cm2inch(width), cm2inch(width)))

# this should be changed for making a panel of multiple figures
ax = fig.add_subplot(111)



filename = 'corrmat.dat'
map = N.loadtxt(filename)
map = map.T

im0 = imshow(map, interpolation='nearest', vmin=-0.01, vmax=0.01)

cbar = colorbar(im0, ax=ax, format="%.3f")


plt.savefig("corrmat.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

