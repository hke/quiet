#!/bin/bash

# this creates .gif files for any fits files in the directory you run it in
# it also creates .png files for any pixel files in the directory you run it in

for i in *.fits; do patch2img patch_6a $i -P -min -0.1 -max 0.1; echo 'converting map: ' $i; done

for j in *pixels.txt; do python /usit/titan/u1/newburgh/repository/quiet_svn/oslo/src/python/pix2img.py $j $j.png -p patch_6a -m -0.1 -M 0.1 ; echo 'converting pixel file: ' $j ; done 
