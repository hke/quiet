#!/bin/bash

# this removes the info at the top of accept files, when ran from the directory where the accept lists are,
# it will do it for all accept files in the directory

for i in *_accept.txt; do
   grep -v '#' $i > temp.txt
   mv temp.txt $i
   done 
