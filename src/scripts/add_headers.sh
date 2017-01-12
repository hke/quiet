#!/bin/bash

# this adds the header to the  accept files, when ran from the directory where the accept lists are,
# it will do it for all accept files in the directory

for i in *_accept.txt; do
    echo 'adding header to: ' $i
    baso=`echo $i | cut -d'_' -f1-3` 
    bast=`echo $baso | cut -d'-' -f4`
    bad=`pwd`
    basdi=`echo $bad | cut -d'/' -f1-6`
    run=$basdi/cuts/$bast
    basr=`echo $bast | cut -d'_' -f3`

    echo '# header_start' > temp.dat
    grep '#' $run >> temp.dat
    echo '# info: weather1, static, dead' >> temp.dat
    echo '# jack-knife: ' $basr >> temp.dat
    echo '# header_end' >> temp.dat
    cat temp.dat $i > temp2.dat
    mv temp2.dat $i
    rm temp.dat
    
   
    mv temp.txt $i
    done 
