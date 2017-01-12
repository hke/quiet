from setup_matplotlib import *
import matplotlib.cm as cm
import numpy as N
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

directory="."
#prefix='testbed1046'
#prefix=sys.argv[1]
#print prefix
#prefix='test'

# Create the plot
width=300
fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)))

# this should be changed for making a panel of multiple figures
#f, axarr = plt.subplots(1,1)
fig = plt.figure(figsize=(cm2inch(width), cm2inch(width)))

# this should be changed for making a panel of multiple figures
ax = fig.add_subplot(111)


filename = 'corrmat.dat'
map = N.loadtxt(filename)
map = map.T


cmap = matplotlib.cm.jet
cmap.set_bad('white',1.)
ax.set_title('Correlation matrix')
im0 = ax.imshow(map, interpolation='nearest', cmap=cmap, vmin=-1, vmax=1)
divider0 = make_axes_locatable(ax)
# Append axes to the right of ax3, with 5% width of ax3
cax0 = divider0.append_axes("right", size="5%", pad=0.15)
# Create colorbar in the appended axes
# Tick locations can be set with the kwarg `ticks`
# and the format of the ticklabels with kwarg `format`
cbar0 = plt.colorbar(im0, cax=cax0, format="%.0f")
# Remove xticks from ax3
#axarr[0].xaxis.set_visible(False)
#axarr[0].yaxis.set_visible(False)
# Manually set ticklocations
#axarr[0].set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

#cbar0 = f.colorbar(cax0, orientation='horizontal')




plt.savefig("corrmat.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)
#fig.show()
