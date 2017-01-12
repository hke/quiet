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

data = np.loadtxt(prefix+'_freq0001_Pa.dat') # l, C_l

# Create the plot
width=8.8
fig = plt.figure(figsize=(1.4*cm2inch(width), cm2inch(width)))

# this should be changed for making a panel of multiple figures
ax = fig.add_subplot(111)

# power spectrum
plt.plot(data[:,0], data[:,1], "k", label=prefix)

xmin = -30
xmax = 60
ymax = 1.1*max(data[:,1])

plt.plot([0,0], [0,ymax], "k", linestyle=':')

mu    = sum(data[:,0]*data[:,1]) * (data[1,0]-data[0,0])
sigma = N.sqrt(sum((data[:,0]-mu)**2*data[:,1]) * (data[1,0]-data[0,0]))

plt.text(xmin+0.05*(xmax-xmin), 0.93*ymax, r"$A = {mu} \pm {sigma}$".format(mu=mu, sigma=sigma))

# labels
plt.xlabel(r"$A \quad [\mu\mathrm{K}]$"); plt.ylabel(r"$P(A)$")
ax.yaxis.labelpad = 10*width/17.; ax.xaxis.labelpad = 10*width/17. # distance of axis label to tick labels

# reduce ticks for small figures
if width < 10:
    ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
    
# grid
#plt.grid(True, which="major", axis="both")

# legend
leg = plt.legend(frameon=True)
# remove box around legend
leg.get_frame().set_edgecolor("white")
leg.get_frame().set_alpha(.8)

# axes limits
plt.ylim([0, ymax]);
plt.xlim([xmin, xmax]);
#plt.ylim([-2, 2]); plt.xlim([0, 18]);

# reduce white space around figure
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)

# set vertical y axis ticklables
for ticklabel in ax.yaxis.get_ticklabels():
    ticklabel.set_rotation("vertical")


plt.savefig(prefix+"_freq0001_Pa.pdf", bbox_inches='tight', bbox_extra_artists=[],pad_inches=0.03)

