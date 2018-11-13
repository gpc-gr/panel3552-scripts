#!/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import csv
import os
import re
import subprocess
from configparser import SafeConfigParser


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--idtable',
        type=os.path.abspath,
        default='/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample/idtable-current.tsv'
    )
    parser.add_argument(
        '--dry-run',
        action='store_true'
    )
    parser.add_argument(
        'batch_root',
        type=os.path.abspath,
    )
    args = parser.parse_args()

    #
    #if args.batch_root.startswith('/u4/share1/home/gpc-gr/'):
    #    batch_root = args.batch_root
    #elif args.batch_root.startswith('/share1/home/gpc-gr/'):
    #    batch_root = '/u4' + args.batch_root
    #elif args.batch_root.startswith('/home/gpc-gr/'):
    #    batch_root = '/u4/share1' + args.batch_root
    #else:
    #    raise Exception(args.batch_root)
    batch_root = args.batch_root

    #
    env_name = os.path.basename(batch_root) + '-env.ini'
    env_path = os.path.join(batch_root, env_name)

    config = SafeConfigParser()
    with open(env_path) as fin:
        config.readfp(fin)
    
    #expected_grid_cell = config.get('grid', 'cell')
    #actual_grid_cell = os.environ.get('SGE_CELL')
    #if actual_grid_cell != expected_grid_cell:
    #    raise Exception('Unexpected SGE_CELL: expected: {}, actual: {}'.format(expected_grid_cell, actual_grid_cell))

    #
    id_list_name = os.path.basename(batch_root) + '-sample_ids.tsv'
    id_list_path = os.path.join(batch_root, id_list_name)

    with open(id_list_path) as fin:
        target_tmmids = [r['tmmid'] for r in csv.DictReader(fin, delimiter='\t')]
      
    #
    reseq_makefile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'reseq.mk')

    with open(args.idtable) as fin:
        for record in csv.DictReader(fin, delimiter='\t'):
            #
            if record['tmmid'] not in target_tmmids:
                continue

            #
            sample_root = os.path.join(batch_root, record['tmmid'])
            fastqs = [n for i, n in enumerate(sorted(set(record['fastqs'].split(',')))) if i % 2 == 0]
            assert len(fastqs) == 2

            options = [
                'ROOT={}'.format(sample_root),
                'TMMID={}'.format(record['tmmid']),
                'FASTQ1={}'.format(re.sub('_R1(_(\d+))?.fastq.gz', '', fastqs[0])),
                'FASTQ2={}'.format(re.sub('_R1(_(\d+))?.fastq.gz', '', fastqs[1])),
                'SEX={}'.format(record['sex']),
                'INCLUDE_NORMAL_BAM_METRICS=true',
                'INCLUDE_NORMAL_GVCF=true'
            ]
            if args.dry_run:
                options.append('DRY_RUN=true')
           
            print('[INFO] Processing {}'.format(record['tmmid']))
            os.chdir(sample_root)
            subprocess.check_call(['make', '-f', reseq_makefile] + options + ['reseq'])


if __name__ == '__main__':
    main()

