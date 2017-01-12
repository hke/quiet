from setup_matplotlib import *
import matplotlib.cm as cm
import numpy as N
from pylab import *
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable

directory="."
#prefix='testbed1046'
#prefix='test'
prefix=sys.argv[1]


# Create the plot
width=300
fig = plt.figure(figsize=(cm2inch(width), 1.5*cm2inch(width)))

# this should be changed for making a panel of multiple figures
f, axarr = plt.subplots(2,1)

filename = prefix+'_tod.dat'
data = N.loadtxt(filename)

data = data[0:300000,:]

axarr[0].plot(data[:,0], data[:,2])
axarr[0].plot(data[:,0], data[:,2]-data[:,1],color='red')
axarr[0].plot([0,max(data[:,0])], [0,0], linestyle='--')
axarr[0].xaxis.set_visible(False)
axarr[0].set_xlim(0,max(data[:,0]))
axarr[0].set_ylim(0.9*min(data[:,2]),1.1*max(data[:,2]))
axarr[0].set_ylabel('TOD [spec unit]')

axarr[1].plot(data[:,0], data[:,1])
axarr[1].plot([0,max(data[:,0])], [0,0], linestyle='--')
axarr[1].set_xlim([0,max(data[:,0])])
axarr[1].set_xlabel('Time [sec]')
axarr[1].set_ylabel('TOD [spec unit]')


#axarr[0].ylim([-ymax, ymax]);
#axarr[0].xlim([0, xmax]);

plt.savefig(prefix+"_tod.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

