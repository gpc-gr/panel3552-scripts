#!/u4/share1/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import csv
import glob
import os
import subprocess

import jinja2


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reference-fasta',
        dest='reference_fasta',
        type=os.path.abspath,
        default='/u4/share1/home/gpc-gr/panel-build37/work/00-reference/hs37d5.fa'
    )
    parser.add_argument(
        '--batch-root',
        dest='batch_root',
        type=os.path.abspath,
        default='/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample'
    )
    parser.add_argument(
        '--chromosome',
        dest='chromosome',
        default=None
    )
    parser.add_argument(
        'panel'
    )
    parser.add_argument(
        'output_root',
        type=os.path.abspath
    )
    args = parser.parse_args()

    #
    reference_id = os.path.splitext(os.path.basename(args.reference_fasta))[0]
    reference_fasta = args.reference_fasta
    reference_region_bed = os.path.join(os.path.dirname(args.reference_fasta), reference_id + '-window_3m.bed')

    #
    source_gvcfs = _collect_sample_gvcfs(args.batch_root, args.panel)
    source_gvcf_groups = [source_gvcfs[i:i+500] for i in range(0, len(source_gvcfs), 500)]

    #
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))
    template = env.get_template('s01-joint_genotyping-panel-v2.sh.j2')
    
    with open(reference_region_bed) as fin:
        for record in csv.reader(fin, delimiter='\t'):
            #
            chromosome = record[0]
            start = int(record[1])
            end = int(record[2])
            output = os.path.join(
                args.output_root,
                '01-chunk',
                os.path.basename(args.output_root) + '-{}_{:09d}_{:09d}.vcf.gz'.format(chromosome, start, end))
            
            if chromosome not in set(map(str, range(1, 22+1))):
                continue
            if args.chromosome and (args.chromosome != chromosome):
                continue

            if os.path.exists(output):
                continue

            #
            script_name = '{}-{}-s51-joint_genotyping-{}_{:09d}_{:09d}.sh'.format(
                'gpcgr_joint',
                os.path.basename(args.output_root),
                chromosome,
                start,
                end)
            script = template.render({
                'script_name': script_name,
                'module_path': os.environ['MODULEPATH'],
                'reference_fasta': reference_fasta,
                'batch': os.path.basename(args.output_root),
                'root': os.path.join(args.output_root, '01-chunk'),
                'region': '{}:{}-{}'.format(chromosome, start, end),
                'padding': 1 * 1000,
                'source_gvcf_groups': source_gvcf_groups,
                'output': output 
            })

            with open(os.path.join(args.output_root, '01-chunk', 'logs', 'orig', script_name), 'w') as fout:
                fout.write(script)

    #subprocess.check_call(['qsub', script_name])


def _collect_sample_gvcfs(batch_root, panel):
    #
    target_tmmids = {}
    with open(os.path.join(batch_root, 'idtable-current.tsv')) as fin:
        for record in csv.DictReader(fin, delimiter='\t'):
            if record['panel'] <= panel:
                target_tmmids[record['tmmid']] = record['panel'], record['igid']

    #
    records = []
    for path in glob.glob(os.path.join(batch_root, 'batch_*/batch_*-sample_ids.tsv')):
        #
        batch_id = os.path.basename(os.path.dirname(path))
        
        #
        with open(path) as fin:
            for record in csv.DictReader(fin, delimiter='\t'):
                tmmid = record['tmmid']
                if tmmid not in target_tmmids:
                    continue

                panel, igid = target_tmmids[tmmid]
                gvcf = os.path.join(batch_root, batch_id, tmmid, tmmid + '.bwamem.hc3.g.vcf.gz')
                records.append((panel, tmmid, igid, gvcf))

    #
    records = sorted(records, key=lambda r: (r[0], r[2]))
    return [r[3] for r in records]


if __name__ == '__main__':
    main()

