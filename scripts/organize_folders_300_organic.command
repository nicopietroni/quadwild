#!/bin/bash

for i in ../../test/300/Organic/*.obj; do


    complete_file="${i}"
    only_file_name=$(basename $complete_file)
    dir_name="${complete_file%.*}"
    dir_name+="/"

    mkdir -p "$dir_name" | cut -f 1 -d '.'

    cp $i $dir_name

 done
