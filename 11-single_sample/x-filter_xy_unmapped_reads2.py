#!/usr/bin/env pypy

from __future__ import print_function

import argparse
import os
import sys


def main():
    #
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    #
    records = {}

    for line in sys.stdin:
        #
        if line.startswith('@'):
            sys.stdout.write(line)
            continue

        #
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
                sys.stdout.write(records.pop(id))
                sys.stdout.write(line)
            else:
                records[id] = line
   
    #
    print('[WARN] {} records remailing'.format(len(records)), file=sys.stderr)


if __name__ == '__main__':
    main()

