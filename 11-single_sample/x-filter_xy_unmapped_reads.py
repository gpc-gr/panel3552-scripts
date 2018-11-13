#!/usr/bin/env pypy

from __future__ import print_function

import argparse
import io
import os
import subprocess
import sys


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', type=os.path.abspath)
    parser.add_argument('--threads', type=int, default=1)
    args = parser.parse_args()

    #
    records = {}

    for line in _read_bam(args.bam, args.threads):
        cols = line.strip().split('\t')
        id = cols[0]
        chromosome1 = cols[2]
        chromosome2 = cols[6]
        sequence = cols[9]

        if len(sequence) not in (162, 259):
            # TODO: consider clipped reads
            #print('length = {}'.format(len(sequence)), file=sys.stderr)
            continue

        if ((chromosome1 == '*') and (chromosome2 == '*')) or (chromosome1 in ('X', 'Y')) or (chromosome2 in ('X', 'Y')):
            if id in records:
                _write_fastq(records.pop(id))
                _write_fastq(cols)
            else:
                records[id] = cols
   
    #
    print('[WARN] {} records remailing'.format(len(records)), file=sys.stderr)


def _read_bam(path, num_threads):
    process = subprocess.Popen("""
        export MODULEPATH=/u4/share1/home/gpc-gr/local/modulefiles
        module load samtools/1.6
        
        samtools view --threads {} -f1 {}
    """.format(num_threads, path), shell=True, stdout=subprocess.PIPE, bufsize=-1)

    with io.open(process.stdout.fileno(), closefd=False) as fin:
        for line in fin:
            yield line

    process.wait()


def _write_fastq(cols):
    print('@' + cols[0])
    print(cols[9])
    print('+')
    print(cols[10])


if __name__ == '__main__':
    main()

