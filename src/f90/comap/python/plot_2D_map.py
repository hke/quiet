from setup_matplotlib import *
import matplotlib.cm as cm
import numpy as N
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

directory="."
#prefix='testbed1046'
prefix=sys.argv[1]
#print prefix
#prefix='test'

# Create the plot
width=300
fig = plt.figure(figsize=(cm2inch(width), 1.5*cm2inch(width)))

# this should be changed for making a panel of multiple figures
f, axarr = plt.subplots(4,1)

filename = prefix+'_freq0001_map.dat'
map = N.loadtxt(filename)
map = map.T

masked_array = N.ma.array (map, mask=(map==0.))
cmap = matplotlib.cm.jet
cmap.set_bad('white',1.)
axarr[0].set_title('Map [uK]')
im0 = axarr[0].imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=-300, vmax=300)
divider0 = make_axes_locatable(axarr[0])
# Append axes to the right of ax3, with 5% width of ax3
cax0 = divider0.append_axes("right", size="5%", pad=0.15)
# Create colorbar in the appended axes
# Tick locations can be set with the kwarg `ticks`
# and the format of the ticklabels with kwarg `format`
cbar0 = plt.colorbar(im0, cax=cax0, format="%.0f")
# Remove xticks from ax3
axarr[0].xaxis.set_visible(False)
axarr[0].yaxis.set_visible(False)
# Manually set ticklocations
#axarr[0].set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

#cbar0 = f.colorbar(cax0, orientation='horizontal')


filename = prefix+'_freq0001_rms.dat'
rms = N.loadtxt(filename)
rms = rms.T

masked_array = N.ma.array (map/rms, mask=(rms==0.))
cmap = matplotlib.cm.jet
cmap.set_bad('white',1.)
axarr[1].set_title('Map/RMS [sigma]')
im1 = axarr[1].imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=-3, vmax=3)
divider1 = make_axes_locatable(axarr[1])
# Append axes to the right of ax3, with 5% width of ax3
cax1 = divider1.append_axes("right", size="5%", pad=0.15)
# Create colorbar in the appended axes
# Tick locations can be set with the kwarg `ticks`
# and the format of the ticklabels with kwarg `format`
cbar1 = plt.colorbar(im1, cax=cax1, format="%.1f")
# Remove xticks from ax3
#axarr[0].xaxis.set_visible(False)
axarr[1].xaxis.set_visible(False)
axarr[1].yaxis.set_visible(False)
# Manually set ticklocations
#axarr[0].set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

#cbar0 = f.colorbar(cax0, orientation='horizontal')


masked_array = N.ma.array (rms, mask=(rms==0.))
cmap = matplotlib.cm.jet
cmap.set_bad('white',1.)
axarr[2].set_title('RMS [uK]')
im2 = axarr[2].imshow(masked_array, interpolation='nearest', cmap=cmap, vmin=0, vmax=300)
divider2 = make_axes_locatable(axarr[2])
# Append axes to the right of ax3, with 5% width of ax3
cax2 = divider2.append_axes("right", size="5%", pad=0.15)
# Create colorbar in the appended axes
# Tick locations can be set with the kwarg `ticks`
# and the format of the ticklabels with kwarg `format`
cbar2 = plt.colorbar(im2, cax=cax2, format="%.0f")
# Remove xticks from ax3
#axarr[0].xaxis.set_visible(False)
axarr[2].xaxis.set_visible(False)
axarr[2].yaxis.set_visible(False)
# Manually set ticklocations
#axarr[0].set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

#cbar0 = f.colorbar(cax0, orientation='horizontal')


filename = prefix+'_freq0001_nhit.dat'
map = N.loadtxt(filename)
map = map.T

masked_array = N.ma.array (map, mask=(map==0.))
cmap = matplotlib.cm.jet
cmap.set_bad('white',1.)
axarr[3].set_title('Hit count')
im3 = axarr[3].imshow(masked_array, interpolation='nearest', cmap=cmap)
divider3 = make_axes_locatable(axarr[3])
# Append axes to the right of ax3, with 5% width of ax3
cax3 = divider3.append_axes("right", size="5%", pad=0.15)
# Create colorbar in the appended axes
# Tick locations can be set with the kwarg `ticks`
# and the format of the ticklabels with kwarg `format`
cbar3 = plt.colorbar(im3, cax=cax3, format="%.0f")
# Remove xticks from ax3
axarr[3].xaxis.set_visible(False)
axarr[3].yaxis.set_visible(False)
#axarr[0].xaxis.set_visible(False)
# Manually set ticklocations
#axarr[0].set_yticks([0.0, 2.5, 3.14, 4.0, 5.2, 7.0])

#cbar0 = f.colorbar(cax0, orientation='horizontal')



plt.savefig(prefix+".pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

