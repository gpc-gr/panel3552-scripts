#!/usr/bin/env python
# -*- coding: utf-8 -*-


import gzip
import subprocess
import sys
import os.path
import vcf_converter


def system(command):
    subprocess.call(command, shell=True)


def mk_dir(dir_name):
    system('mkdir -p ' + dir_name)


def set_tabix(vcf_file):

    command  = 'bgzip -f ' + vcf_file + ' > ' + vcf_file + '.gz;'
    command += 'tabix -p vcf ' + vcf_file + '.gz;'
    system(command)


def main(argv):

    source_vcf = 'analysis/joint_genotyping/mt_panel_3.5KJPN.het2.snv.vcf'
    target_vcf = 'analysis/release/p3552_3.5kjpn_alt_MT/full/p3552_3.5kjpn_alt-full-MT.vcf'
    mk_dir(os.path.dirname(target_vcf))
    system('cp ' + source_vcf + ' ' + target_vcf)
    set_tabix(target_vcf)

    source_vcf = 'analysis/joint_genotyping/mt_panel_4KJPN.het2.snv.vcf'
    target_vcf = 'analysis/release/p4007_4kjpn_alt_MT/full/p4007_4kjpn_alt-full-MT.vcf'
    mk_dir(os.path.dirname(target_vcf))
    system('cp ' + source_vcf + ' ' + target_vcf)
    set_tabix(target_vcf)


if __name__ == '__main__':
    argv = sys.argv
    main(argv)
