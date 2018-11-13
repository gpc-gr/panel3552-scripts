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
    source = os.path.join(args.data_root, '02-merge', batch + '.vcf.gz')

    #
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))

    #
    template = env.get_template('s04-vqsr_snv.sh.j2')
    script_name = '{}-{}-s04-vqsr_snv.sh'.format(args.job_prefix, batch)
    script_path = os.path.join(args.data_root, '03-vqsr', 'logs', script_name)

    with open(script_path, 'w') as fout:
        fout.write(template.render({
            'job_prefix': args.job_prefix,
            'batch': batch,
            'root': os.path.join(args.data_root, '03-vqsr'),
            'reference_fasta': args.reference_fasta,
            'resource_prefix': '/u4/share1/home/gpc-gr/panel-build37/work/03-gatk_resource_bundle/vqsr',
            'module_path': '/u4/share1/home/gpc-gr/local/modulefiles',
            'output_prefix': os.path.join(args.data_root, '03-vqsr', batch),
            'source': source
        }))

    #
    template = env.get_template('s04-vqsr_indel.sh.j2')
    script_name = '{}-{}-s04-vqsr_indel.sh'.format(args.job_prefix, batch)
    script_path = os.path.join(args.data_root, '03-vqsr', 'logs', script_name)

    with open(script_path, 'w') as fout:
        fout.write(template.render({
            'job_prefix': args.job_prefix,
            'batch': batch,
            'root': os.path.join(args.data_root, '03-vqsr'),
            'reference_fasta': args.reference_fasta,
            'resource_prefix': '/u4/share1/home/gpc-gr/panel-build37/work/03-gatk_resource_bundle/vqsr',
            'module_path': '/u4/share1/home/gpc-gr/local/modulefiles',
            'output_prefix': os.path.join(args.data_root, '03-vqsr', batch),
            'source': source
        }))


if __name__ == '__main__':
    main()

