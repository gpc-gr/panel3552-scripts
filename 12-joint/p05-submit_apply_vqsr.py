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
    parser.add_argument('data_root', type=os.path.abspath)
    parser.add_argument('--job-prefix', default='gpcgr_joint')
    args = parser.parse_args()

    #
    reference_id = os.path.splitext(os.path.basename(args.reference_fasta))[0]
    reference_fasta = args.reference_fasta

    batch = os.path.basename(args.data_root)

    #
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))

    template = env.get_template('s05-apply_vqsr.sh.j2')
    script_name = '{}-{}-s05-apply_vqsr.sh'.format(args.job_prefix, batch)
    script_path = os.path.join(args.data_root, '03-vqsr', 'logs', script_name)

    with open(script_path, 'w') as fout:
        fout.write(template.render({
            'module_path': '/u4/share1/home/gpc-gr/local/modulefiles',
            'job_prefix': args.job_prefix,
            'batch': batch,
            'reference_fasta': args.reference_fasta,
            'region': '11',
            'source_vcf': os.path.join(args.data_root, '02-merge', batch + '.vcf.gz'),
            'source_snv_recal': os.path.join(args.data_root, '03-vqsr', batch + '.snv.recall'),
            'source_snv_tranches': os.path.join(args.data_root, '03-vqsr', batch + '.snv.tranches'),
            'snv_filter_level': '99.5',
            'source_indel_recal': os.path.join(args.data_root, '03-vqsr', batch + '.indel.recall'),
            'source_indel_tranches': os.path.join(args.data_root, '03-vqsr', batch + '.snv.tranches'),
            'indel_filter_level': '97.0',
            'output_vcf_prefix': os.path.join(args.data_root, '03-vqsr', batch)
        }))


if __name__ == '__main__':
    main()

