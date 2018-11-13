import csv
import glob
import itertools
import os
import shutil


target_tmmids = set()
with open('idtable-current.tsv') as fin:
    for record in csv.DictReader(fin, delimiter='\t'):
        if record['sex'] == '1':
            target_tmmids.add(record['tmmid'])


for batch_name in os.listdir('.'):
    if (not batch_name.startswith('batch_')) or batch_name.endswith('_test'):
        continue

    with open(os.path.join(batch_name, batch_name + '-sample_ids.tsv')) as fin:
        batch_tmmids = set(r['tmmid'] for r in csv.DictReader(fin, delimiter='\t'))

    for tmmid in target_tmmids.intersection(batch_tmmids):
        sample_root = os.path.join(batch_name, tmmid)
        files = itertools.chain(
            glob.glob(os.path.join(sample_root, 'logs', '*PAR2*')),
            glob.glob(os.path.join(sample_root, '*PAR2*'))
        )

        for path in files:
            if path.endswith('.orig'):
                continue

            source = path
            destination = path + '.orig'

            print(source)
            shutil.move(source, destination)

