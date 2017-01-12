import numpy as np
import matplotlib.pyplot as plt

nfiles = 22
filenames = [] 
for i in range(nfiles):
    filenames.append("timestream"+str(i)+".txt")

Data = []
for name in filenames:
    Data.append(np.loadtxt(name))
data2 = []
data2 = np.vstack(Data)

time = data2[:,0]
ra = data2[:,1]; dec = data2[:,2]
az = data2[:,3]; el = data2[:,4]

time[:] = (time[:]-time[0])*24*3600

# Flags, what to plot
hitmap = True
ra_dec_all = False
ra_dec_some = False
subplot_one = False
subplot_three = False

if ra_dec_all:
    plt.plot(ra,dec)
    plt.ylabel('Declination (deg)')
    plt.xlabel('Right Ascension (deg)')

if hitmap:
    plt.figure()
    plt.hist2d(ra,dec,bins=50)
    plt.ylabel('Declination (deg)')
    plt.xlabel('Right Ascension (deg)')
    plt.colorbar()

if ra_dec_some:
    # Right ascension vs declination for different elevations
    plt.figure()
    plt.plot(ra[np.argwhere(el==el[0])],dec[np.argwhere(el==el[0])],'r',label=el[0])
    plt.plot(ra[np.argwhere(el==el[9000])],dec[np.argwhere(el==el[9000])],'b',label=el[9000])
    plt.plot(ra[np.argwhere(el==el[20000])],dec[np.argwhere(el==el[20000])],'k',label=el[20000])
    plt.plot(ra[np.argwhere(el==el[-1])],dec[np.argwhere(el==el[-1])],'g',label=el[-1])
    plt.ylabel('Declination (deg)')
    plt.xlabel('Right Ascension (deg)')
    plt.legend()

if subplot_one:
    plt.figure()
    plt.subplot(4,1,1)
    plt.plot(time,az)
    plt.ylabel('Azimuth (deg)')
    plt.subplot(4,1,2)
    plt.plot(time,el)
    plt.ylabel('Elevation (deg)')
    plt.subplot(4,1,3)
    plt.plot(time,ra)
    plt.ylabel('Right ascension (deg)')
    plt.subplot(4,1,4)
    plt.plot(time,dec)
    plt.ylabel('Declination (deg)')
    plt.xlabel('time (s)')
    #plt.subplot_tool()

if subplot_three:
    # This part plots az, al, ra, and dec as a function of time
    plt.figure()
    plt.subplot(4,3,1)
    plt.plot(time[np.argwhere(el==el[0])],az[np.argwhere(el==el[0])],'r')
    plt.ylabel('Azimuth (deg)')
    plt.subplot(4,3,2)
    plt.plot(time[np.argwhere(el==el[7000])],az[np.argwhere(el==el[7000])],'b')
    plt.subplot(4,3,3)
    plt.plot(time[np.argwhere(el==el[-1])],az[np.argwhere(el==el[-1])],'g')
    plt.subplot(4,3,4)
    plt.plot(time[np.argwhere(el==el[0])],el[np.argwhere(el==el[0])],'r')
    plt.ylabel('Elevation (deg)')
    plt.subplot(4,3,5)
    plt.plot(time[np.argwhere(el==el[7000])],el[np.argwhere(el==el[7000])],'b')
    plt.subplot(4,3,6)
    plt.plot(time[np.argwhere(el==el[-1])],el[np.argwhere(el==el[-1])],'g')
    plt.subplot(4,3,7)
    plt.plot(time[np.argwhere(el==el[0])],ra[np.argwhere(el==el[0])],'r')
    plt.ylabel('Right ascension (deg)')
    plt.subplot(4,3,8)
    plt.plot(time[np.argwhere(el==el[7000])],ra[np.argwhere(el==el[7000])],'b')
    plt.subplot(4,3,9)
    plt.plot(time[np.argwhere(el==el[-1])],ra[np.argwhere(el==el[-1])],'g')
    plt.subplot(4,3,10)
    plt.plot(time[np.argwhere(el==el[0])],dec[np.argwhere(el==el[0])],'r')
    plt.ylabel('Declination (deg)')
    plt.xlabel('time (s)')
    plt.subplot(4,3,11)
    plt.plot(time[np.argwhere(el==el[7000])],dec[np.argwhere(el==el[7000])],'b')
    plt.xlabel('time (s)')
    plt.subplot(4,3,12)
    plt.plot(time[np.argwhere(el==el[-1])],dec[np.argwhere(el==el[-1])],'g')
    plt.xlabel('time (s)')
    #plt.subplot_tool()


plt.show()

