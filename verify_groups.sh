#!/bin/bash

xloc_names=($(cut -f1 ../order-data-files/groups.tsv | sort -u))

for xloc in "${xloc_names[@]}"
do
    grep -P "$xloc\t" ../order-data-files/groups.tsv > /tmp/$xloc.1.tmp
    grep -P "$xloc\t" ../order-data-files/groups2.tsv | awk '{gsub("Arabidopsis\t", ""); print}' > /tmp/$xloc.2.tmp
#     echo "#######"
#     echo "#######"
#     echo "#######"
#     echo $xloc
#     echo "#######"
#     echo "#######"
#     echo "#######"
    diff /tmp/$xloc.1.tmp /tmp/$xloc.2.tmp
    rm /tmp/$xloc.*
done
