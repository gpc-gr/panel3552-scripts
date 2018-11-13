#!/usr/bin/env python
# -*- coding: utf-8 -*-


import gzip
import subprocess
import sys
import os.path
import vcf_converter


JAVA     = '/usr/local/pkg/java/jdk1.8.0_144/bin/java -XX:+UseSerialGC'
GATK_JAR = '/home/gpc-gr/local/packages/gatk-3.7-0-gcfedb67/GenomeAnalysisTK.jar'
#GATK_JAR = 'jar_files/GenomeAnalysisTK.jar'


class Setting:

    def __init__(self, rate, vcf_suffix, vcf_file_prefix_list, output_prefix):
        self.rate = rate
        self.vcf_suffix = vcf_suffix
        self.output_prefix = output_prefix
        self.vcf_file_prefix_list = vcf_file_prefix_list


def system(command):
    subprocess.call(command, shell=True)


def mk_dir(dir_name):
    system('mkdir -p ' + dir_name)


def write_lines(lines, output_file):
    mk_dir(os.path.dirname(output_file))
    with open(output_file, 'w') as fp:
        for line in lines:
            fp.write(line)
            fp.write('\n')


def get_contig_size(reference, contig_name):
    contig_size = 0
    with open(reference, 'r') as fp:
        in_chrM_flag = False
        for line in fp:
            if line.startswith('>'):
                if line.startswith('>' + contig_name):
                    in_chrM_flag = True
                else:
                    in_chrM_flag = False
                continue
            if in_chrM_flag:
                contig_size += len(line.strip())
    return contig_size


def convert_to_org_position(position, base_shift, mt_size):
    position += base_shift
    if position > mt_size:
        position -= mt_size
    return position


def get_head_tail_records(vcf_file, base_shift, mt_size, margin):
    head_records = []
    tail_records = []

    with open(vcf_file) as fp:
        for line in fp:
            if line.startswith('#CHROM'):
                break
        variant_dict = {}
        for line in fp:
            items = line.strip().split()
            position = int(items[1])
            org_position = convert_to_org_position(position, base_shift, mt_size)
            items[1] = str(org_position)
            if org_position < margin:
                head_records.append('\t'.join(items))
            elif org_position >= mt_size - margin:
                tail_records.append('\t'.join(items))
    return head_records, tail_records


def merge_vcf(vcf_file1, vcf_file2, base_shift, mt_size, margin, output_file):
    head_records, tail_records = get_head_tail_records(vcf_file2, base_shift, mt_size, margin)

    vcf_lines = []
    with open(vcf_file1) as fp:
        for line in fp:
            vcf_lines.append(line.strip())
            if line.startswith('#CHROM'):
                break
        vcf_lines.extend(head_records)
        for line in fp:
            items = line.strip().split()
            position = int(items[1])
            if position < margin:
                continue
            elif position >= mt_size - margin:
                continue
            vcf_lines.append(line.strip())
    vcf_lines.extend(tail_records)
    write_lines(vcf_lines, output_file)


def get_tmm_id_list(vcf_file):
    tmm_id_list  = []
    root, ext = os.path.splitext(vcf_file)
    if ext == '.gz':
        fp = gzip.open(vcf_file)
    else:
        fp = open(vcf_file)
    for line in fp:
        if line.startswith('#CHROM'):
            items = line.strip().split()
            for i in xrange(9, len(items)):
                tmm_id = items[i]
                tmm_id_list.append(tmm_id)
            break
    return tmm_id_list


def get_vcf_file_prefix_list(tmm_id_list):
    return map(
        lambda x : '/home/gpc-gr/panel-build37/work/13-mt/mt_pipeline/debug_analysis/' + x + '/vcf/' + x,
        tmm_id_list)


def joint_genotyping(vcf_files, reference_fasta, options, output_prefix):
    mk_dir(os.path.dirname(output_prefix))
    write_lines(vcf_files, output_prefix + '.list')

    command  = JAVA + ' -jar ' + GATK_JAR
    command += ' -T GenotypeGVCFs'
    command += ' --max_alternate_alleles 100'
    command += ' -R ' + reference_fasta
    command += ' -V ' + output_prefix + '.list'
    if options is not None:
        command += ' ' + ' '.join(options)
    command += ' -o ' + output_prefix + '.vcf;'
    command += 'rm ' + output_prefix + '.list;'
    system(command)


def analysis(vcf_file_prefix_list, reference_fasta_prefix, shift_size, mtDNA_size, margin, setting):
    output_prefix = setting.output_prefix
    vcf_suffix = setting.vcf_suffix

    options = [
        '--heterozygosity '       + str(0.001   * setting.rate),
        '--indel_heterozygosity ' + str(1.25E-4 * setting.rate)]

    reference_fasta = reference_fasta_prefix + '_with_mtDNA.fa'
    joint_genotyping(
        map(lambda x : x + '.' + vcf_suffix, vcf_file_prefix_list),
        reference_fasta,
        options,
        output_prefix)


    reference_fasta = reference_fasta_prefix + '_with_shifted_mtDNA.fa'
    joint_genotyping(
        map(lambda x : x + '_shifted.' + vcf_suffix, vcf_file_prefix_list),
        reference_fasta,
        options,
        output_prefix + '.shifted')

    merge_vcf(
        output_prefix + '.vcf',
        output_prefix + '.shifted.vcf',
        shift_size,
        mtDNA_size,
        margin,
        output_prefix + '.merged.vcf')

    vcf_converter.main(output_prefix + '.merged.vcf', output_prefix)


def main(argv):
    shift_size = 10000
    margin     = 2500

    org_mtDNA_fasta = 'org_data/rCRS.fa'
    mtDNA_size = get_contig_size(org_mtDNA_fasta, 'MT')

    reference_fasta_prefix = 'data/fasta/hs37d5'

    panel_types = []
    vcf_file_dict = {}

    panel_type = '4KJPN'
    vcf_file = '/home/gpc-gr/panel-build37/work/13-mt/02-merge/p4007_4kjpn_alt_MT/p4007_4kjpn_alt_MT.vcf.gz'
    vcf_file_dict.update({panel_type : vcf_file})
    panel_types.append(panel_type)

    panel_type = '3.5KJPN'
    vcf_file = '/home/gpc-gr/panel-build37/work/13-mt/02-merge/p3552_3.5kjpn_alt_MT/p3552_3.5kjpn_alt_MT.vcf.gz'
    vcf_file_dict.update({panel_type : vcf_file})
    panel_types.append(panel_type)

    setting_list = []
    for panel_type in panel_types:
        vcf_file = vcf_file_dict[panel_type]
        tmm_id_list = get_tmm_id_list(vcf_file)
        vcf_file_prefix_list = get_vcf_file_prefix_list(tmm_id_list)
        output_prefix = 'analysis/joint_genotyping/mt_panel_' + panel_type + '.het2'
        setting = Setting(2,   'het2.vcf',   vcf_file_prefix_list, output_prefix)
        setting_list.append(setting)
        output_prefix = 'analysis/joint_genotyping/mt_panel_' + panel_type + '.het2.5'
        setting = Setting(2.5, 'het2.5.vcf', vcf_file_prefix_list, output_prefix)
        setting_list.append(setting)

    for setting in setting_list:
        vcf_file_prefix_list = setting.vcf_file_prefix_list
        analysis(vcf_file_prefix_list, reference_fasta_prefix, shift_size, mtDNA_size, margin, setting)


if __name__ == '__main__':
    argv = sys.argv
    main(argv)
