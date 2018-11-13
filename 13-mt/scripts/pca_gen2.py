#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
import subprocess
import sys


R     = '/usr/local/pkg/R/3.4.2/bin/R'
PLINK = '/usr/local/pkg/plink/1.90b5.1/plink'


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


def get_EAS_keep_file(pop_file, output_file):
    lines = []
    with open(pop_file) as fp:
        for line in fp:
            break
        for line in fp:
            items = line.strip().split()
            if items[2] == 'EAS':
                lines.append(items[0] + ' ' + items[0] + ' 0 0 0 -9')
    write_lines(lines, output_file)


def make_bed(vcf_files, keep_file, output_prefix):
    mk_dir(os.path.dirname(output_prefix))
    system('rm -f ' + output_prefix + '.tped')
    for i, vcf_file in enumerate(vcf_files):
        i += 1
        command  = ''
        command += PLINK
        command += ' --vcf ' + vcf_file
        command += ' --maf 0.01'
        command += ' --hwe 0.00001'
        if keep_file is not None:
            command += ' --keep ' + keep_file
        command += ' --recode transpose'
        command += ' --out ' + output_prefix + '_' + str(i)
        command += '; '
        command += 'cat ' + output_prefix + '_' + str(i) + '.tped'
        command += ' >> ' + output_prefix + '.tped'
        command += '; '
        command += 'rm -f ' + output_prefix + '_' + str(i) + '.tped'
        command += '; '

        if i == 1:
            command += 'mv ' + output_prefix + '_' + str(i) + '.tfam'
            command += ' '   + output_prefix + '.tfam'
            command += ';'
        else:
            command += 'rm -f ' + output_prefix + '_' + str(i) + '.tfam'
            command += ';'
        system(command)

    command  = ''
    command += PLINK
    command += ' --tped ' + output_prefix + '.tped'
    command += ' --tfam ' + output_prefix + '.tfam'
    command += ' --make-bed'
    command += ' --out ' + output_prefix
    system(command)


def prepare_renamed_SNP_bim(bim_file, output_file):
    lines = []
    name_count_dict = {}
    with open(bim_file) as fp:
        for line in fp:
            items = line.strip().split()
            SNP_name = items[1]
            if SNP_name == '.':
                continue
            if SNP_name not in name_count_dict:
                name_count_dict.update({SNP_name : 0})
    with open(bim_file) as fp:
        for line in fp:
            items = line.strip().split()
            SNP_name = items[1]
            if SNP_name == '.':
                SNP_name = 'chr' + items[0] + ':' +  items[3] + ':' + items[4] + ':' + items[5]
            if SNP_name in name_count_dict:
                name_count_dict[SNP_name] += 1
                if name_count_dict[SNP_name] > 1:
                    SNP_name += ':' + str(name_count_dict[SNP_name])
            else:
                name_count_dict.update({SNP_name : 1})
            items[1] = SNP_name
            lines.append(' '.join(items))
    write_lines(lines, output_file)


def is_same_allele_pair(alleles1, alleles2):
    if alleles1[0] == alleles2[0] and alleles1[1] == alleles2[1]:
        return True
    if alleles1[0] == alleles2[1] and alleles1[1] == alleles2[0]:
        return True
    return False


def get_common_bim(bim_file1, bim_file2, output_prefix):
    bim_dict = {}
    with open(bim_file1) as fp:
        for line in fp:
            items = line.strip().split()
            key = items[0] + ':' + items[3]
            alleles = [items[4], items[5]]
            exclude_flag = False
            for allele in alleles:
                if allele == '*' or allele == 0 or len(allele) > 1:
                    exclude_flag = True
                    break
            if exclude_flag:
                continue
            if key not in bim_dict:
                bim_dict.update({key : [alleles, line.strip(), False]})
    lines1 = []
    lines2 = []
    with open(bim_file2) as fp:
        for line in fp:
            items = line.strip().split()
            key = items[0] + ':' + items[3]
            alleles = [items[4], items[5]]
            if key not in bim_dict:
                continue
            if bim_dict[key][2]:
                continue
            if is_same_allele_pair(alleles, bim_dict[key][0]):
                lines1.append(bim_dict[key][1])
                lines2.append(line.strip())
                bim_dict[key][2] = True
    write_lines(lines1, output_prefix + '_1.common_bim')
    write_lines(lines2, output_prefix + '_2.common_bim')


def merge_bed(bed_prefix1, bed_prefix2, output_prefix):
    prepare_renamed_SNP_bim(bed_prefix1 + '.bim', bed_prefix1 + '.bim.renamed')
    prepare_renamed_SNP_bim(bed_prefix2 + '.bim', bed_prefix2 + '.bim.renamed')
    get_common_bim(bed_prefix1 + '.bim.renamed', bed_prefix2 + '.bim.renamed', output_prefix)

    system('rm -f ' + output_prefix + '.ped')
    for i, bed_prefix in enumerate([bed_prefix1, bed_prefix2]):
        i += 1
        command  = ''
        command += PLINK
        command += ' --bed ' + bed_prefix + '.bed'
        command += ' --fam ' + bed_prefix + '.fam'
        command += ' --bim ' + bed_prefix + '.bim.renamed'
        command += ' --extract ' + output_prefix + '_' + str(i) + '.common_bim'
        command += ' --make-bed'
        command += ' --out ' + output_prefix + '_' + str(i)
        command += '; '

        command += PLINK
        command += ' --bfile ' + output_prefix + '_' + str(i)
        command += ' --out '   + output_prefix + '_' + str(i)
        command += ' --recode'
        command += '; '

        command += 'cat ' + output_prefix + '_' + str(i) + '.ped'
        command += ' >> ' + output_prefix + '.ped'
        command += ';'
        system(command)

    command  = ''
    command += PLINK
    command += ' --ped ' + output_prefix + '.ped'
    command += ' --map ' + output_prefix + '_2.map'
    command += ' --make-bed'
    command += ' --out ' + output_prefix
    command += ';'
    system(command)


def pca(bed_prefix, pop_file, output_prefix):
    mk_dir(os.path.dirname(output_prefix))

    command  = PLINK
    command += ' --bfile ' + bed_prefix
    command += ' --maf  0.05'
    command += ' --geno 0.01'
    command += ' --hwe  0.05'
    command += ' --make-bed'
    command += ' --out ' + output_prefix
    system(command)

    command  = PLINK
    command += ' --bfile ' + output_prefix
    command += ' --autosome'
    command += ' --indep-pairwise 200 4 0.1'
    command += ' --out ' + output_prefix
    system(command)

    command  = PLINK
    command += ' --bfile ' + output_prefix
    command += ' --extract ' + output_prefix + '.prune.in'
    command += ' --out ' + output_prefix
    command += ' --pca header'
    system(command)

    info_dict = {}
    color_dict = {'ToMMo': 'magenta', 'JPT': 'blue', 'CHB' : 'red', 'CHS': 'yellow', 'KHV' : 'green', 'CDX': 'brown'}
    with open(output_prefix + '.eigenvec') as fp:
        for line in fp:
            break
        for line in fp:
            items = line.strip().split()
            sample_key = items[0] + ' ' + items[1]
            info_dict.update({sample_key : [items[2], items[3], 'magenta']})
    with open(pop_file) as fp:
        for line in fp:
            break
        for line in fp:
            items = line.strip().split()
            sample_key = items[0] + ' ' + items[0]
            pop = items[1]
            if sample_key in info_dict:
                info_dict[sample_key][2] = color_dict[pop]

    lines = []
    for sample_key in info_dict.keys():
        if info_dict[sample_key][2] == 'magenta':
            lines.append(' '.join(info_dict[sample_key]))
    for sample_key in info_dict.keys():
        if info_dict[sample_key][2] != 'magenta':
            lines.append(' '.join(info_dict[sample_key]))
    write_lines(lines, output_prefix + '.dat')

    labels = ['ToMMo', 'JPT', 'CHB', 'CHS', 'CDX', 'KHV']

    r_command  = 'data.table <- read.table("' + output_prefix + '.dat");'
    r_command += 'output_file <- "' + output_prefix + '.pdf";'
    r_command += 'pdf(output_file);'
    r_command += 'plot(x=as.numeric(data.table[,1]), y=as.numeric(data.table[,2]), col=as.character(data.table[,3]), xlab="PC1", ylab="PC2", pch=20);'
    r_command += 'legend("bottomright", legend=c("' + '","'.join(labels) + '"), col=c("' + '","'.join([color_dict[label] for label in labels]) + '"), pch=20, lty="blank");'
    r_command += 'dev.off();'
    system(R + ' -q -e \'' + r_command + '\'')


def main():
    onekgp_vcf_file = '/home/gpc-gr/local/downloads/1000genomes_phase3/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
    pop_file        = '/share2/home/gpc-gr/personal/tadaka/20180423-mt_panel_comparison/onekgp/integrated_call_samples_v3.20130502.ALL.panel'

    onekgp_vcf_files = ['/home/gpc-gr/local/downloads/1000genomes_phase3/ALL.chr' + str(chr_num) + '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' for chr_num in range(1, 23)]

    ToMMo_vcf_files  = ['/home/gpc-gr/panel-build37/release/panel/p3552_3.5kjpn_alt/vqsr_99.5_99.0/p3552_3.5kjpn_alt-passed_99.5_99.0-' + str(chr_num) + '.vcf.gz' for chr_num in range(1, 23)]

    onekgp_EAS_bed_prefix = 'pca/bed/onekgp_EAS'
    ToMMo_bed_prefix      = 'pca/bed/ToMMo_3.5kjpn'

    EAS_keep_file = onekgp_EAS_bed_prefix + '.keep'
    get_EAS_keep_file(pop_file, EAS_keep_file)

    #make_bed(onekgp_vcf_files, EAS_keep_file, onekgp_EAS_bed_prefix)
    #make_bed(ToMMo_vcf_files,  None,          ToMMo_bed_prefix)

    merged_bed_prefix = 'pca/bed/ToMMo_3.5kjpn_with_EAS'
    #merge_bed(ToMMo_bed_prefix, onekgp_EAS_bed_prefix, merged_bed_prefix)

    pca(merged_bed_prefix, pop_file, 'pca/images2/ToMMo_3.5kjpn_with_EAS')


if __name__ == '__main__':
    main()
