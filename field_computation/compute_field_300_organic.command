#!/bin/bash

for i in ../../test/300/Organic/*.obj; do


    complete_file="${i}"
    only_file_name=$(basename $complete_file)
    dir_name="${complete_file%.*}"
    dir_name+="/"


    new_filename="${dir_name}""${only_file_name}"

    ./field_computation "${new_filename}" basic_setup_organic.txt batch
done
