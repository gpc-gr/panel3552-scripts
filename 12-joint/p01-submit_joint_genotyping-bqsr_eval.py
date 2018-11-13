#!/u4/share1/home/gpc-gr/panel-build37/work/env3/bin/python3

from __future__ import print_function

import argparse
import csv
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
        '--source-root',
        dest='source_roots',
        required=True,
        nargs='+',
        type=os.path.abspath
    )
    parser.add_argument(
        '--output-root',
        dest='output_root',
        required=True,
        type=os.path.abspath
    )
    parser.add_argument(
        '--use-bqsred-gvcf',
        action='store_true'
    )
    args = parser.parse_args()

    #
    reference_id = os.path.splitext(os.path.basename(args.reference_fasta))[0]
    reference_fasta = args.reference_fasta
    reference_region_bed = os.path.join(os.path.dirname(args.reference_fasta), reference_id + '-window_3m.bed')

    #
    source_gvcfs = _collect_sample_gvcfs(args.source_roots, args.use_bqsred_gvcf)

    #
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))
    template = env.get_template('s01-joint_genotyping.sh.j2')
    
    with open(reference_region_bed) as fin:
        for record in csv.reader(fin, delimiter='\t'):
            #
            chromosome = record[0]
            start = int(record[1])
            end = int(record[2])

            if chromosome != '11':
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
                'gvcfs': source_gvcfs,
                'output': os.path.basename(args.output_root) + '_{}_{:09d}_{:09d}.vcf.gz'.format(chromosome, start, end)
            })

            with open(os.path.join(args.output_root, '01-chunk', 'logs', script_name), 'w') as fout:
                fout.write(script)

    #subprocess.check_call(['qsub', script_name])


def _collect_sample_gvcfs(source_roots, use_bqsred_gvcf):
    result = []
    for source_root in source_roots:
        for tmmid in _load_batch_sample_ids(source_root):
            sample_gvcf_name = tmmid + '.bwamem' + ('.bqsr' if use_bqsred_gvcf else '') + '.hc3.g.vcf.gz'
            sample_gvcf_path = os.path.join(source_root, tmmid, sample_gvcf_name)
            assert os.path.exists(sample_gvcf_path)

            result.append(sample_gvcf_path)

    return result


def _load_batch_sample_ids(batch_root):
    path = os.path.join(batch_root, os.path.basename(batch_root) + '-sample_ids.tsv')
    with open(path) as fin:
        return [r['tmmid'] for r in csv.DictReader(fin, delimiter='\t')]


if __name__ == '__main__':
    main()

