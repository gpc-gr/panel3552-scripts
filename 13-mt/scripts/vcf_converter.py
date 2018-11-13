#!/usr/bin/env python
# -*- coding: utf-8 -*-


import gzip
import subprocess
import sys
import os.path


NO_REWARD_POSITIONS = set({513, 16179, 16193})
DEL = '<DEL>'


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


class Variants:

    def __init__(self, position, ref, alts, genotypes):
        self.position      = position
        self.ref           = ref
        self.org_genotypes = genotypes
        self.genotypes     = list(genotypes)
        self.ref_fragments = self.convert_to_fragments(ref)
        self.fragments_dict = {ref : self.ref_fragments}
        for alt in alts:
            REWARD = 0.1
            if position in NO_REWARD_POSITIONS:
                REWARD = 0
            fragments = get_fragments(ref, alt, REWARD)
            self.fragments_dict.update({alt : fragments})
        for i, genotype in enumerate(genotypes):
            if genotype is None or genotype == '.':
                continue
            if genotype == '*':
                self.genotypes[i] = self.ref[0]
            else:
                self.genotypes[i] = self.fragments_dict[genotype][0]

    def convert_to_fragments(self, fragment):
        fragments = []
        for c in fragment:
            fragments.append(c)
        return fragments


def get_fragments(ref, alt, REWARD):
    if len(ref) == 1 and len(alt) == 1:
        return [alt]
    aligned_ref_sequence, aligned_alt_sequence = Smith_Waterman(ref, alt, REWARD)
    fragments = []
    fragment = None
    for i in xrange(len(aligned_ref_sequence)):
        if aligned_ref_sequence[i] == '_':
            fragment.append(aligned_alt_sequence[i])
        else:
            if fragment is not None:
                fragments.append(''.join(fragment))
            fragment = []
            if aligned_alt_sequence[i] != '_':
                fragment.append(aligned_alt_sequence[i])
    if fragment is not None:
        fragments.append(''.join(fragment))
    for i, fragment in enumerate(fragments):
        if fragment == '':
            fragments[i] = DEL
    return fragments


class SW_CELL:

    def __init__(self, cost, trace):
        self.cost  = cost
        self.trace = trace


def Smith_Waterman(sequence1, sequence2, REWARD):
    num_rows    = len(sequence1)
    num_columns = len(sequence2)

    SW_matrix = [None] * (num_rows + 1)
    for r in xrange(num_rows + 1):
        SW_matrix[r] = [None] * (num_columns + 1)
        for c in xrange(num_columns + 1):
            SW_matrix[r][c] = SW_CELL(None, None)

    INDEL_COST    = 6
    MISMATCH_COST = 4
    HORIZONTAL = 'H'
    VERTICAL   = 'V'
    DIAGONAL   = 'D'

    INF_COST = float('inf')
    SW_matrix[0][0].cost = 0
    for r in xrange(1, num_rows + 1):
        SW_matrix[r][0].cost  = SW_matrix[r-1][0].cost + INF_COST
        SW_matrix[r][0].trace = VERTICAL
    for c in xrange(1, num_columns + 1):
        SW_matrix[0][c].cost  = SW_matrix[0][c-1].cost + INF_COST
        SW_matrix[0][c].trace = HORIZONTAL
    for r in xrange(1, num_rows + 1):
        for c in xrange(1, num_columns + 1):
            best_trace = DIAGONAL
            best_cost = SW_matrix[r-1][c-1].cost
            if sequence1[r-1] != sequence2[c-1]:
                best_cost += MISMATCH_COST
            cost = SW_matrix[r-1][c].cost
            cost += INDEL_COST
            if c == num_columns:
                cost -= REWARD
            if best_cost > cost:
                best_cost  = cost
                best_trace = VERTICAL
            cost = SW_matrix[r][c-1].cost
            cost += INDEL_COST
            if best_cost > cost:
                best_cost  = cost
                best_trace = HORIZONTAL
            SW_matrix[r][c].cost  = best_cost
            SW_matrix[r][c].trace = best_trace

    aligned_sequence1 = []
    aligned_sequence2 = []

    r = num_rows
    c = num_columns

    s1 = num_rows    - 1
    s2 = num_columns - 1

    while True:
        trace = SW_matrix[r][c].trace
        if trace == VERTICAL:
            aligned_sequence1.append(sequence1[s1])
            aligned_sequence2.append('_')
            r  -= 1
            s1 -= 1
        elif trace == HORIZONTAL:
            aligned_sequence1.append('_')
            aligned_sequence2.append(sequence2[s2])
            c  -= 1
            s2 -= 1
        else:
            aligned_sequence1.append(sequence1[s1])
            aligned_sequence2.append(sequence2[s2])
            r  -= 1
            c  -= 1
            s1 -= 1
            s2 -= 1
        if r == 0 and c == 0:
            break
    aligned_sequence1.reverse()
    aligned_sequence2.reverse()
    return aligned_sequence1, aligned_sequence2


def get_vcf_header_lines(vcf_file):
    header_lines = []
    root, ext = os.path.splitext(vcf_file)
    if ext == '.gz':
        fp = gzip.open(vcf_file)
    else:
        fp = open(vcf_file)
    for line in fp:
        header_lines.append(line.strip())
        if line.startswith('#CHROM'):
            break
    fp.close()
    return header_lines


def get_variants_list(vcf_file):
    variants_list = []
    root, ext = os.path.splitext(vcf_file)
    if ext == '.gz':
        fp = gzip.open(vcf_file)
    else:
        fp = open(vcf_file)
    for line in fp:
        if line.startswith('#CHROM'):
            items = line.strip().split()
            samples = items[9:]
            break
    for line in fp:
        variants = get_variants(line.strip())
        variants_list.append(variants)
    fp.close()
    return variants_list, samples


def get_variants(record):
    items = record.split()
    position = int(items[1])
    ref  = items[3]
    alts = items[4].split(',')
    alleles = [ref]
    alleles.extend(alts)

    genotypes = []
    for i in xrange(9, len(items)):
        if items[i].split(':')[0] == '.':
            genotypes.append('.')
            continue
        genotype_num = int(items[i].split(':')[0])
        genotype = alleles[genotype_num]
        genotypes.append(genotype)
    return Variants(position, ref, alts, genotypes)


def get_variants_dict(vcf_file):
    variants_list, samples = get_variants_list(vcf_file)
    sample_size = len(samples)
    variants_dict = {}
    for variants in variants_list:
        variants_dict.update({variants.position : variants})
    for variants in variants_list:
        position = variants.position
        ref = variants.ref
        ref_fragments = variants.ref_fragments
        ref_size = len(ref)
        for i, genotype in enumerate(variants.org_genotypes):
            if genotype == '.' \
                or genotype == '*' \
                or genotype == ref:
                continue
            for step in xrange(1, ref_size):
                downstream_position = position + step
                if downstream_position not in variants_dict:
                    downstream_variants = Variants(downstream_position, ref[step], [], [None] * sample_size)
                    variants_dict.update({downstream_position : downstream_variants})
                downstream_allele   = variants.fragments_dict[genotype][step]
                downstream_variants = variants_dict[downstream_position]
                downstream_genotype = downstream_variants.genotypes[i]

                if downstream_genotype is None \
                    or downstream_genotype == '.' \
                    or downstream_genotype == downstream_variants.ref[0]:
                    downstream_variants.genotypes[i] = downstream_allele
        variants.ref = ref[0]
    return variants_dict, samples


def get_vcf_record(variants, samples, snv_flag):
    alt_allele_count_dict = {}
    ref = variants.ref
    genotypes = variants.genotypes
    for i, genotype in enumerate(genotypes):
        if genotype is None:
            genotypes[i] = ref
            continue
        if genotype == '.':
            continue
        if snv_flag:
            if genotype == DEL:
                continue
            else:
                genotype = genotype[0]
        if genotype == ref:
            continue
        if genotype not in alt_allele_count_dict:
            alt_allele_count_dict.update({genotype : 0})
        alt_allele_count_dict[genotype] += 1
    alts = alt_allele_count_dict.keys()
    if len(alts) == 0:
        return None
    alts.sort(key = lambda x : -alt_allele_count_dict[x])
    allele_num_dict = {ref : str(0)}
    for i, alt in enumerate(alts):
        allele_num_dict.update({alt : str(i + 1)})
    position = variants.position
    genotype_nums = []
    for genotype in genotypes:
        if snv_flag:
            if genotype == DEL:
                genotype = '.'
            else:
                genotype = genotype[0]
        if genotype == '.':
            genotype_nums.append(genotype)
        else:
            genotype_num = allele_num_dict[genotype]
            genotype_nums.append(genotype_num)
    line  = '\t'.join(['MT', str(position), '.', ref, ','.join(alts), '.', '.', '.', 'GT'])
    line += '\t' + '\t'.join(genotype_nums)
    return line


def main(vcf_file, output_prefix):
    variants_dict, samples = get_variants_dict(vcf_file)
    positions = variants_dict.keys()
    positions.sort()

    lines = get_vcf_header_lines(vcf_file)
    snv_flag = False
    for position in positions:
        vcf_record = get_vcf_record(variants_dict[position], samples, snv_flag)
        if vcf_record is not None:
            lines.append(vcf_record)
    write_lines(lines, output_prefix + '.analysis.vcf')

    lines = get_vcf_header_lines(vcf_file)
    snv_flag = True
    for position in positions:
        vcf_record = get_vcf_record(variants_dict[position], samples, snv_flag)
        if vcf_record is not None:
            lines.append(vcf_record)
    write_lines(lines, output_prefix + '.snv.vcf')


if __name__ == '__main__':
    vcf_file = 'analysis/joint_genotyping/mt_panel_4KJPN.merged.vcf'
    output_prefix = 'analysis/joint_genotyping/mt_panel_4KJPN'
    main(vcf_file, output_prefix)
