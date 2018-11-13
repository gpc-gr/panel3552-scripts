#!/u4/share1/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import collections
import gzip
import os
import sys


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('vcf', type=os.path.abspath)
    args = parser.parse_args()

    #
    histogram = collections.defaultdict(int)

    with _open(args.vcf) as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            cols = line.split('\t')
            reference = cols[2]
            alternatives = cols[3].split(',')

            for alternative in alternatives:
                if alternative != '*':
                    histogram[len(alternative) - len(reference)] += 1

    #
    print('length\tfrequency')
    for length, frequency in sorted(histogram.items()):
        print('{}\t{}'.format(length, frequency))


def _open(path):
    if path.endswith('.vcf.gz'):
        return gzip.open(path, 'rt')
    elif path.endswith('.vcf'):
        return open(path, 'r')
    elif path in ('/dev/stdin', '-'):
        return sys.stdin
    elif path in ('/dev/stderr',):
        return sys.stderr
    else:
        raise Exception


if __name__ == '__main__':
    main()

