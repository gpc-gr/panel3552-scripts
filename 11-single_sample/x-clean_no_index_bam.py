#!/home/gpc-gr/panel-build37/work/env3/bin/python3

from __future__ import print_function

import argparse
import glob
import os
import subprocess


def main():
    #
    parser = argparse.ArgumentParser()
    parser.add_argument('root', type=os.path.abspath, default=os.getcwd())
    args = parser.parse_args()

    #
    for name in os.listdir(args.root):
        if not (name.startswith('T') or name.startswith('t')):
            continue

        sample_id = name
        sample_root = os.path.join(args.root, sample_id)
        sample_merged_bam = os.path.join(sample_root, sample_id + '.bwamem.bam')

        if not os.path.exists(sample_merged_bam):
            continue

        sample_noindex_bams = glob.glob(os.path.join(sample_root, '*_L00[12].bwamem.bam'))
        if not sample_noindex_bams:
            continue

        if not all(_check_no_index_bam(p) for p in sample_noindex_bams):
            continue

        print('\n'.join(sample_noindex_bams))


def _check_no_index_bam(path):
    script = """
        JAVA=/usr/local/pkg/java/jdk1.8.0_144/bin/java
        PICARD_JAR=/u4/share1/home/gpc-gr/local/packages/picard-2.10.6/picard-2.10.6.jar

        $JAVA -jar $PICARD_JAR CheckTerminatorBlock I={bam} 2>/dev/null
        echo $?
    """.format(bam=path)

    output = subprocess.check_output(script, shell=True).decode('utf-8')
    return output.strip() == '0'


if __name__ == '__main__':
    main()

