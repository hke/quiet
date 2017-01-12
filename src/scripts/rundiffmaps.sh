#!/bin/bash

# this creates difference maps, and differenced noise pixels files for any set of .a,.b files in the 
# directory you run it in

for i in *a_map.fits; do
   
   filea=$i
   bam=`echo $i | cut -d'.' -f1`
   fileb=$bam.b_map.fits
   if [ -f $fileb ]; then
      diffmap=$bam.diff_map.fits
      echo 'differencing maps: ', $filea, $fileb, 'putting into : ' $diffmap
      /usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/makeDifferenceMap.py $filea $fileb $diffmap
   fi 
   done   

for j in *a_pixels.txt; do
   filea=$j
   bap=`echo $j | cut -d'.' -f1`
   fileb=$bap.b_pixels.txt
   if [ -f $fileb ]; then
     diffpix=$bap.diff_pixels.txt
     echo 'differencing pixels: ' $filea, $fileb, 'putting into: ', $diffpix
     /usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/differenceNoiseMap.py $filea $fileb > $diffpix
   fi
   done
