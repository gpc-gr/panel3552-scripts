#!/bin/bash

set -eu

source_root=/u4$(cd $1; pwd)
destination_root=/u4$(mkdir -p $2; cd $2; pwd)

source_id=$(basename $source_root | sed 's/\.wrong_name//g')
destination_id=$(basename $destination_root)

sub_dir_names=("01-chunk" "02-merge" "03-vqsr")
for sub_dir_name in "${sub_dir_names[@]}"; do
    mkdir -p $destination_root/$sub_dir_name
    for source in $(find $source_root/$sub_dir_name ! -type d -mindepth 1 -maxdepth 1); do
        destination=$destination_root/$sub_dir_name/$(basename $source | sed "s/${source_id}/${destination_id}/g")
        echo "* $source"
        echo "    => $destination"
        ln -s $source $destination
    done
done

