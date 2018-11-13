#!/u4/share1/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/python3

from __future__ import print_function

import argparse
import gzip
import glob
import os
import re
import subprocess
import sys

import jinja2


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('source_root', type=os.path.abspath)
    args = parser.parse_args()

    #
    file_region_map = {}
    for path in glob.glob(os.path.join(args.source_root, '*.vcf.gz')):
        if '.g.vcf.gz' in path:
            continue
   
        try:
            chromosome, start, end, padding = _get_chunk_region(path)
        except:
            continue
        if not chromosome.isdigit():
            continue

        file_region_map[path] = chromosome, start, end, padding

    #
    print('type\tchromosome\tchromosomal_position\treference\talternative\tvcf1\tvcf2')

    for path1, (chromosome1, start1, end1, padding1) in file_region_map.items():
        for path2, (chromosome2, start2, end2, padding2) in file_region_map.items():
            #
            if (path1 >= path2) or (chromosome1 != chromosome2) or (end1 != start2):
                continue

            #
            comp_chromosome = chromosome1
            comp_start = end1 - padding2
            comp_end = end1 + padding1
            comp_start2 = start2 - padding2
            comp_end2 = start2 + padding1
            comp_region = '{}:{}-{}'.format(comp_chromosome, comp_start, comp_end)
            assert (comp_start == comp_start2) and (comp_end == comp_end2)
          
            print('[INFO] chunk1={}, chunk2={}, start={}, end={}'.format(
                os.path.basename(path1), os.path.basename(path2), comp_start, comp_end
            ), file=sys.stderr)

            #
            variants1 = _read_vcf(path1, comp_region)
            variants2 = _read_vcf(path2, comp_region)

            for key in sorted(set(variants1) - set(variants2)):
                row = list(key) + [os.path.basename(path1), os.path.basename(path2), comp_start, comp_end]
                print('VCF1_ONLY\t' + '\t'.join(map(str, row)))

            for key in sorted(set(variants2) - set(variants1)):
                row = list(key) + [os.path.basename(path1), os.path.basename(path2), comp_start, comp_end]
                print('VCF2_ONLY\t' + '\t'.join(map(str, row)))

            for key in sorted(set(variants1).intersection(set(variants2))):
                if variants1[key] == variants2[key]:
                    continue

                row = list(key) + [os.path.basename(path1), os.path.basename(path2), comp_start, comp_end]
                print('COMMON\t' + '\t'.join(map(str, row)))


def _get_chunk_region(path):
    with gzip.open(path, 'rt') as fin:
        for line in fin:
            if line.startswith('##GATKCommandLine.GenotypeGVCFs'):
                return _extract_region(line)
            
    raise Exception


def _extract_region(line):
    chromosome = None
    start = None
    end = None
    padding = None

    for col in line.split():
        if col.startswith('intervals'):
            match = re.search('\[([a-zA-Z0-9_\.]+):(\d+)-(\d+)\]', col)
            chromosome = match.group(1)
            start = int(match.group(2))
            end = int(match.group(3))
    
        elif col.startswith('interval_padding'):
            padding = int(col.split('=')[1])

    for value in (chromosome, start, end, padding):
        assert value is not None

    return chromosome, start, end, padding


def _read_vcf(path, region):
    command = """
        export MODULEPATH=/u4/share1/home/gpc-gr/local/modulefiles
        module load bcftools/1.6

        bcftools view --regions {region} {path}
    """.format(region=region, path=path)
    output = subprocess.check_output(command, shell=True).decode('utf-8')

    variants = {}
    for line in output.splitlines():
        if line.startswith('#'):
            continue

        cols = line.strip().split('\t')
        chromosome = cols[0]
        chromosomal_position = cols[1]
        ref = cols[3]
        alts = cols[4]
        genotypes = tuple(s.split(':')[0] for s in cols[9:])

        variants[chromosome, chromosomal_position, ref, alts] = genotypes

    return variants


if __name__ == '__main__':
    main()

