#!/bin/bash

echo "Processing all OBJs in $1"

for i in $1/*.obj; do

    complete_file="${i}"

    echo "Processing File ${i}"

    only_file_name=$(basename $complete_file)
    dir_name="${complete_file%.*}"
    dir_name+="/"

    mkdir -p "$dir_name" | cut -f 1 -d '.'

    cp $i $dir_name

    new_filename="${dir_name}""${only_file_name}"

    ./quadwild  "${new_filename}" basic_setup_mechanical.txt
done
