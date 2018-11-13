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
    command += ' --indep-pairwise 200 100 0.1'
    command += ' --out ' + output_prefix
    system(command)

    command  = PLINK
    command += ' --bfile ' + output_prefix
    command += ' --extract ' + output_prefix + '.prune.in'
    command += ' --out ' + output_prefix
    command += ' --pca header'
    system(command)

    info_dict = {}
    color_dict = {'ToMMo': 'magenta', 'OUTLIER': 'black'}
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

    r_command  = 'data.table <- read.table("' + output_prefix + '.dat");'
    r_command += 'output_file <- "' + output_prefix + '.pdf";'
    r_command += 'pdf(output_file);'
    r_command += 'plot(x=as.numeric(data.table[,1]), y=as.numeric(data.table[,2]), col=as.character(data.table[,3]), xlab="PC1", ylab="PC2", pch=20);'
    r_command += 'dev.off();'
    system(R + ' -q -e \'' + r_command + '\'')


def main():
    ToMMo_bed_prefix  = 'pca/bed/ToMMo_3.5kjpn'
    merged_bed_prefix = 'pca/bed/ToMMo_3.5kjpn_only'

    system('ln -s ' + os.path.basename(ToMMo_bed_prefix) + '.bed '         + merged_bed_prefix +'.bed')
    system('ln -s ' + os.path.basename(ToMMo_bed_prefix) + '.bim.renamed ' + merged_bed_prefix +'.bim')
    system('ln -s ' + os.path.basename(ToMMo_bed_prefix) + '.fam '         + merged_bed_prefix +'.fam')

    outliers = ['tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******',
                'tmm*******']

    outlier_file = 'pca/data/outliners.txt'
    write_lines(['sample pop'] + [outlier + ' OUTLIER' for outlier in outliers], outlier_file)

    pca(merged_bed_prefix, outlier_file, 'pca/images/' + os.path.basename(merged_bed_prefix))


if __name__ == '__main__':
    main()
