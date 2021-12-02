#!/bin/bash

for i in ../../test/300/Organic/*.obj; do


    complete_file="${i}"
    only_file_name=$(basename $complete_file)
    dir_name="${complete_file%.*}"
    dir_name+="/"


    new_filename="${dir_name}""${only_file_name}"
    new_filename2=$(echo $new_filename | sed 's/\(.*\)\..*/\1/')
     new_filename2="$new_filename2""_rem.obj"
    echo "$new_filename2"

    ./field_tracing "${new_filename2}" batch
done
