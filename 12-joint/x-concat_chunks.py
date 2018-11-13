#!/u4/share1/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import csv
import gzip
import glob 
import os
import sys


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('--disconcordant-variant-list')
    parser.add_argument('sources', nargs='+', type=os.path.abspath)
    args = parser.parse_args()

    #
    source_paths = list(sorted(args.sources))

    #
    blacklisted_variants = set()
    if args.disconcordant_variant_list:
        with open(args.disconcordant_variant_list) as fin:
            for record in csv.DictReader(fin, delimiter='\t'):
                key = '{chromosome}:{chromosomal_position}'.format(**record)
                blacklisted_variants.add(key)

    #
    processed_variants = set()
    previous_position = 0

    for index, path in enumerate(source_paths):
        with gzip.open(path, 'rt') as fin:
            for line in fin:
                if line.startswith('#'):
                    if index == 0:
                        sys.stdout.write(line)

                else:
                    cols = line.split('\t')
                    key = '{}:{}'.format(cols[0], cols[1])
                    if (key in blacklisted_variants) or (key in processed_variants):
                        continue
                    
                    current_position = int(cols[1])
                    if previous_position >= current_position:
                        print('[WARN] disconcordant variant found: {}'.format(key), file=sys.stderr)
                        continue

                    processed_variants.add(key)
                    previous_position = current_position
                    sys.stdout.write(line)


if __name__ == '__main__':
    main()

