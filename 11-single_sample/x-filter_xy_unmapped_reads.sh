#!/bin/bash

set -eu
trap 'exit 100' ERR
trap 'exit 100' SIGXCPU

source=$1
output=$2


/home/gpc-gr/panel-build37/work/11-single_sample/scripts/x-filter_xy_unmapped_reads2.py $source --threads $NSLOTS\
    | gzip -c\
    > $output

