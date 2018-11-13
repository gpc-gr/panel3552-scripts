#!/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import csv
import glob
import os


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--root',
        type=os.path.abspath,
        default='/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample'
    )
    parser.add_argument(
        '--idtable',
        type=os.path.abspath,
        default='/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample/idtable-current.tsv'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=100
    )
    parser.add_argument(
        '--force',
        action='store_true'
    )
    parser.add_argument(
        'batch_name'
    )
    parser.add_argument(
        'grid_cell'
    )
    args = parser.parse_args()

    #
    already_processed_tmmids = set()
    for path in sorted(glob.glob(os.path.join(args.root, 'batch_*/batch_*-sample_ids.tsv'))):
        print('[INFO] Reading {}'.format(path))

        with open(path) as fin:
            tmmids = set(r['tmmid'] for r in csv.DictReader(fin, delimiter='\t'))
            already_processed_tmmids.update(tmmids)

    #
    target_records = []
    with open(args.idtable) as fin:
        for record in csv.DictReader(fin, delimiter='\t'):
            #
            if record['tmmid'] in already_processed_tmmids:
                continue

            if record['fastqs'] == '-':
                print('[WARN] FASTQ missing: {}'.format(record['tmmid']))
                continue

            #
            target_records.append(record)
            if len(target_records) >= args.batch_size:
                break

    if not target_records:
        print('[INFO] No batch left')
        return

    #
    batch_root = os.path.join(args.root, args.batch_name)
    assert args.force or (not os.path.exists(batch_root))
    os.makedirs(batch_root, exist_ok=True)
    
    batch_sample_ids_path = os.path.join(batch_root, '{}-sample_ids.tsv'.format(args.batch_name))
    with open(batch_sample_ids_path, 'w') as fout:
        print('tmmid', file=fout)
        print('\n'.join(r['tmmid'] for r in target_records), file=fout)

    for record in target_records:
        #
        sample_root = os.path.join(batch_root, record['tmmid'])
        if not os.path.exists(sample_root):
            os.makedirs(sample_root)

        sample_log_root = os.path.join(sample_root, 'logs')
        if not os.path.exists(sample_log_root):
            os.makedirs(sample_log_root)

        #
        for fastq_name in record['fastqs'].split(','):
            #
            source = '/u4/share1/home/gpc-gr/panel-build37/work/fastqs/' + fastq_name
            if not os.path.exists(source):
                source = '/u4/share1/primary_analysis_data/illumina/' + fastq_name

            assert os.path.exists(source)

            #
            destination = os.path.join(sample_root, fastq_name)
            os.symlink(source, destination)

    batch_env_path = os.path.join(batch_root, '{}-env.ini'.format(args.batch_name))
    with open(batch_env_path, 'w') as fout:
        print('[grid]', file=fout)
        print('cell = {}'.format(args.grid_cell), file=fout)
        print(file=fout)


if __name__ == '__main__':
    main()

