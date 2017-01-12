#!/bin/bash

# this runs the ppcl code on any diff pixels created in the maps directory, currently it must be run
# from the powspec directory
# it uses cosine apodization (not ideal)

   files=`ls ../maps/*diff_pixels.txt`
   for i in $files; do
     echo $i 
     ba=`echo $i | cut -d'/' -f3` 
     echo $ba
     ma=`echo $i | cut -d'_' -f1-3`_map.fits 
     echo 'commencing ppcl'
     # note that i is the pixels, ma is the map

     # get the map, pixels into the right units:
     sed 's/uk/mk/g' $ma > temp.fits
     mv temp.fits $ma
     grep n $i > header.dat
     grep -v n $i > pix.dat
     python /usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/scalethis.py pix.dat temp.dat
     cat header.dat temp.dat > $i
     rm header.dat
     rm temp.dat
     rm pix.dat
     
     # copy the null spectrum (0 for each el bin): nullspectrum.dat
     # copy the definition of the el bins, and the statement of which window
     # functions to use: est_pure0.dat
     cp /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/nullspectrum.dat .
     cp /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/est_pure0.dat .

     # I. Defining weight functions
     #   ppcl-cos-apodization:  cosine apodization on a spherical cap shaped region
     #   ppcl-fkp:              algorithmically generated weight function using FKP ansatz (fast)
     #   ppcl-optimize:         algorithmically generated weight function using method from astro-ph/0612571 (slow)
     #   ppcl-fkp2:             equivalent to ppcl-fkp, except that it can handle T fluctuation by option -t.
     #   FOR COS APODIZATION: theta = degrees south from the north pole = 90 - latitude
     #                        phi   = RA = longitude
     #   patch_coords = {patch_2:( 185, -42), patch_2a: (181, -39), 
     #     patch_4: (75, -40), patch_4a: (78, -39), patch_6: (10, -45), patch_6a: (12, -48),
     #     patch_7: (330, -30), patch_7a: (332, -36), patch_a: (148, 6),
     #     patch_b: (21, -48), patch_gc:(-94.5, -25)


     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-cos-apodization -m  -t 137.3 -p 12.5 $i 0.0 7.5 > tophat.dat
     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-cos-apodization -t 137.3 -p 12.5 $i 1.73 7.5 > cosine_1.73.dat
     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-cos-apodization -t 137.3 -p 12.5 $i 4.62 7.5 > cosine_4.62.dat


     sed -e 's/tophat_13deg.dat/tophat.dat/g' -e 's/cos_8deg_13deg.dat/cosine_1.73.dat/g' -e 's/cos_3deg_13deg.dat/cosine_4.62.dat/g' est_pure0.dat > est_input.dat   

     # II. Computing transfer matrix
     #   ppcl-compute-xfer:     computes transfer matrix given a set of weight functions and bandpowers
     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-compute-xfer -o est_output.dat est_input.dat 23 

     # try this, can compute the estimator (w/o noise bias subtracted)
     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-estimate est_output.dat $ma > $ba.clest.dat

     # III. Generating simulations
     #   ppcl-simulate:         generates simulation with signal + uncorrelated noise
     #   ppcl-estimate:         computes power spectrum of externally generated simulation, WITHOUT noise bias subtracted
     #   ppcl-mc:               equivalent to Monte Carlo iterations of ppcl-simulate -> ppcl-estimate, but noise bias IS subtracted
     #   ppcl-estimate2:        equivalent to ppcl-estimate, except it cal calculate TT, TE, TB, and EB in addition to EE and BB.
     /usit/titan/u1/newburgh/repository/ppcl/trunk/ppcl/ppcl-mc -n 500 nullspectrum.dat $i 23 est_output.dat > $ba.cls.dat 

   done
