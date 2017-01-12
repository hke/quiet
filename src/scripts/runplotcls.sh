#!/bin/bash

# this runs the ppcl code on any diff pixels created in the maps directory, currently it must be run
# from the powspec directory
# it uses cosine apodization (not ideal)

   files=`ls *.cls.dat`
   for i in $files; do
     echo $i
     erfile=$i
     clfile=`echo $i | cut -d'.' -f1-3`.clest.dat
     grep '#' $clfile > temp.dat
     mv temp.dat $clfile

     python /usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/plotPowerSpec.py $erfile $clfile

     done


