#!/bin/bash

echo 'patch (2a,4a,6a,7b,noise), binsize (5,10,30,60,120)'

rm cuts/weather_sigma_*$2_*_$1*.txt
rm cuts/howmany*$2*$1*.txt
rm cuts/weather_*$2_*$1*.txt


#python plot_allmods_hists.py 5 -1 $1
python plot_allmods_hists.py $2 -1 $1 # 1 sigma threshold

cat cuts/weather_sigma_acceptlist_$2_patch_$1.txt | sort -k 1n,1 -k 2n,2 -k 3n,3  | uniq > cuts/weather_acceptlist_$2_patch_$1.txt_sorted

python make_accept_from_list.py cuts/weather_acceptlist_$2_patch_$1.txt_sorted
cat holder.txt | awk '{print $1, $2, $NF, $0}' | awk '{$4=""; print}' | awk '{$4=""; print}' | awk '{NF--;print}'> cuts/weather_acceptlist_asrun_$2_patch_$1.txt
rm holder.txt
