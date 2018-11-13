#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys

from optparse import OptionParser

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
            position = (int(items[1]) + base_shift) % mt_size
            items[1] = str(position)
            if position < margin:
                head_records.append('\t'.join(items))
            elif position >= mt_size - margin:
                tail_records.append('\t'.join(items))
    return head_records, tail_records

def merge_individual_vcf(vcf_file1, vcf_file2, mt_size, base_shift, margin, output_file):
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
            if position >= margin and position < mt_size - margin:
                vcf_lines.append(line.strip())
    vcf_lines.extend(tail_records)
    write_lines(vcf_lines, output_file)

def write_lines(lines, output_file):
    with open(output_file, 'w') as fp:
        for line in lines:
            fp.write(line)
            fp.write('\n')

def main(argv=None):

    parser = OptionParser()
    parser.add_option("--vcf-file1", dest="vcf_file1", default=None, help="vcf file 1")
    parser.add_option("--vcf-file2", dest="vcf_file2", default=None, help="vcf file 2")
    parser.add_option("--mt-size", dest="mt_size", default=None, help="MT size")
    parser.add_option("--base-shift", dest="base_shift", default=None, help="MT size")
    parser.add_option("--margin", dest="margin", default=None, help="margin")
    parser.add_option("--output-file", dest="output_file", default=None, help="output_file")

    options, args = parser.parse_args(args=sys.argv)

    merge_individual_vcf(
        options.vcf_file1,
        options.vcf_file2,
        int(options.mt_size),
        int(options.base_shift),
        int(options.margin),
        options.output_file)

if __name__ == '__main__':
    main(argv = sys.argv)
