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
        '--output-root',
        dest='output_root',
        required=True,
        type=os.path.abspath
    )
    parser.add_argument(
        '--include-162pe',
        action='store_true'
    )
    parser.add_argument(
        '--include-259pe',
        action='store_true'
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
    source_gvcfs = []
    if args.include_162pe:
        template = '/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample/batch_9999/{id}_162PE/{id}_162PE.bwamem.hc3.g.vcf.gz'
        source_gvcfs.extend(template.format(id=id) for id in ['NA12878', 'NA12891', 'NA12892'])
    if args.include_259pe:
        template = '/u4/share1/home/gpc-gr/panel-build37/work/11-single_sample/batch_9999/{id}_259PE/{id}_259PE.bwamem.hc3.g.vcf.gz'
        source_gvcfs.extend(template.format(id=id) for id in ['NA12878', 'NA12891', 'NA12892'])

    if args.use_bqsred_gvcf:
        source_gvcfs = [p.replace('.bwamem.hc3.g.vcf.gz', '.bwamem.bqsr.hc3.g.vcf.gz') for p in source_gvcfs]
    
    #
    template_dir = os.path.join(os.path.dirname(__file__), 'templates')
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir))
    template = env.get_template('s01-joint_genotyping-bqsr_eval.sh.j2')
    
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
                'output': os.path.basename(args.output_root) + '-{}_{:09d}_{:09d}.vcf.gz'.format(chromosome, start, end)
            })

            script_path = os.path.join(args.output_root, '01-chunk', 'logs', script_name)
            with open(script_path, 'w') as fout:
                fout.write(script)

            subprocess.check_call(['qsub', script_path])


if __name__ == '__main__':
    main()

