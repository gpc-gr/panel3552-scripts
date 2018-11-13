#!/usr/bin/env python
# -*- coding: utf-8 -*-

import subprocess
import sys
import time
import os.path

BWA_MEM    = '/home/gpc-gr/local/packages/bwa-0.7.12/bin/bwa'
SAMTOOLS   = '/home/gpc-gr/local/packages/samtools-1.6/samtools.sh'
JAVA       = '/usr/local/pkg/java/jdk1.8.0_144/bin/java -XX:+UseSerialGC'
GATK_JAR   = '/home/gpc-gr/local/packages/gatk-3.7-0-gcfedb67/GenomeAnalysisTK.jar'
PICARD_JAR = '/home/gpc-gr/local/packages/picard-2.10.6/picard-2.10.6.jar'

class Fasta_setting:

    def __init__(self, org_reference, org_mtDNA_fasta, reference_prefix):
        self.org_reference    = org_reference
        self.org_mtDNA_fasta  = org_mtDNA_fasta
        self.reference_prefix = reference_prefix

class Data_item:
    def __init__(self, bam_file, project_id):
        self.bam_file   = bam_file
        self.project_id = project_id

def execute(command):
    subprocess.call(command, shell=True)

def my_subprocess(command):
    p = subprocess.Popen(
            command,
            stdin  = subprocess.PIPE,
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE,
            shell  = True)
    return p

def mk_dir(dir_name):
    execute('mkdir -p ' + dir_name)

def load_fasta(fasta_file):
    sequence = []
    with open(fasta_file, 'r') as fp:
        for line in fp:
            if line.startswith('>MT'):
                continue
            sequence.extend(line.strip())
    return sequence

def prepare_reference_wo_chrM(reference, output_file):
    mk_dir(os.path.dirname(output_file))
    with open(output_file, 'w') as w_fp:
        with open(reference, 'r') as fp:
            in_chrM_flag = False
            for line in fp:
                if line.startswith('>'):
                    if line.startswith('>MT'):
                        in_chrM_flag = True
                    else:
                        in_chrM_flag = False
                if in_chrM_flag == False:
                    w_fp.write(line)

def merge_fasta(fasta_files, output_file):
    mk_dir(os.path.dirname(output_file))
    for i in xrange(len(fasta_files)):
        fasta_file = fasta_files[i]
        if i == 0:
            execute('cat ' + fasta_file + ' > '   + output_file)
        else:
            execute('cat ' + fasta_file + ' >> '  + output_file)

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

def prepare_shifted_fasta(reference, contig_name, shift_size, output_file):
    mk_dir(os.path.dirname(output_file))

    genome = ''
    one_line_base_count = 0
    with open(reference, 'r') as fp:
        for line in fp:
            if line.startswith('>' + contig_name):
                continue
            if one_line_base_count == 0:
                one_line_base_count = len(line.strip())
            genome += line.strip()

    shifted_genome  = genome[shift_size:]
    shifted_genome += genome[:shift_size]

    shifted_genome_size = len(shifted_genome)

    mk_dir(os.path.dirname(output_file))
    with open(output_file, 'w') as fp:
        fp.write('>' + contig_name)
        fp.write('\n')
        start_position = 0
        while start_position < shifted_genome_size:
            end_position = start_position + one_line_base_count
            if end_position >= shifted_genome_size:
                end_position = shifted_genome_size
            fp.write(shifted_genome[start_position:end_position])
            fp.write('\n')
            start_position = end_position

def prepare_dict(reference):
    output_file = os.path.dirname(reference) + '/' + '.'.join(os.path.basename(reference).split('.')[:-1]) + '.dict'
    command  = ''
    command += 'rm -f ' + output_file + ';'
    command += JAVA
    command5+= ' -jar jar_files/CreateSequenceDictionary.jar'
    command += ' R=' + reference + ' O=' + output_file + ';'
    command += SAMTOOLS + ' faidx ' + reference + ';'
    execute(command)

def prepare_bwa_index(fasta_file):
    execute(BWA_MEM + ' index ' + fasta_file)

def set_index(fasta_file):
    prepare_dict(fasta_file)
    prepare_bwa_index(fasta_file)

def reference_preparation(fasta_setting):
    reference_prefix        = fasta_setting.reference_prefix
    mk_dir(os.path.dirname(reference_prefix))
    org_reference           = fasta_setting.org_reference
    org_mtDNA_fasta         = fasta_setting.org_mtDNA_fasta
    reference_shifted_mtDNA = reference_prefix + '_with_shifted_mtDNA.fa'
    reference               = reference_prefix + '_with_mtDNA.fa'
    reference_wo_mtDNA      = reference_prefix + '_wo_mtDNA.fa'
    shifted_mtDNA_fasta     = reference_prefix + '_shifted_mtDNA_only.fa'

    prepare_shifted_fasta(org_mtDNA_fasta, 'MT', fasta_setting.shift_size, shifted_mtDNA_fasta)

    prepare_reference_wo_chrM(org_reference, reference_wo_mtDNA)

    merge_fasta([reference_wo_mtDNA, org_mtDNA_fasta],     reference)
    set_index(reference)

    merge_fasta([reference_wo_mtDNA, shifted_mtDNA_fasta], reference_shifted_mtDNA)
    set_index(reference_shifted_mtDNA)

def extract_mt_and_unmapped(bam_file, sample_id, output_prefix, qsub_log_dir):
    mk_dir(os.path.dirname(output_prefix))
    execute('ls -lh ' + bam_file)
    command  = ''
    for i in xrange(10):
        command += 'ls -lh ' + bam_file + ';'
    command += SAMTOOLS + ' view -u -f  4 -F264 ' + bam_file + ' -o ' + output_prefix + '_unmapped1.bam;'
    command += SAMTOOLS + ' view -u -f  8 -F260 ' + bam_file + ' -o ' + output_prefix + '_unmapped2.bam;'
    command += SAMTOOLS + ' view -u -f 12 -F256 ' + bam_file + ' -o ' + output_prefix + '_unmapped3.bam;'
    command += SAMTOOLS + ' view -u -F 268 ' + bam_file + ' MT -o ' + output_prefix + '_chrM.bam;'
    tmp_bam_files = []
    tmp_bam_files.append(output_prefix + '_unmapped1.bam')
    tmp_bam_files.append(output_prefix + '_unmapped2.bam')
    tmp_bam_files.append(output_prefix + '_unmapped3.bam')
    tmp_bam_files.append(output_prefix + '_chrM.bam')

    merged_bam_prefix = output_prefix + '_chrM_and_unmapped'
    command += SAMTOOLS + ' merge -f ' + merged_bam_prefix + '.bam ' + ' '.join(tmp_bam_files) + ';'
    command += SAMTOOLS + ' sort -n '  + merged_bam_prefix + '.bam -o ' + merged_bam_prefix + '.sorted.bam;'
    command += JAVA + ' -Xmx4g -Xms4g'
    command += ' -jar jar_files/DuplicatedReadRemover.jar'
    command += ' ' + merged_bam_prefix + '.sorted.bam'
    command += ' | ' + SAMTOOLS + ' view -Sb - -o ' + merged_bam_prefix + '.nodup.bam;'

    for tmp_bam_file in tmp_bam_files:
        command += ' rm -f '  + tmp_bam_file + ';'
    command += ' rm -f '  + merged_bam_prefix + '.bam;'
    command += ' rm -f '  + merged_bam_prefix + '.sorted.bam;'

    job_name_base = 'extract_mt_and_unmapped'
    job_name = job_name_base + '_' + sample_id
    qsub_dir = qsub_log_dir + '/' + job_name_base
    memory_size = '10G'
    qsub(command, qsub_dir, job_name, memory_size, None)

    print('Job Name:   ' + job_name)
    print('Input:      ' + bam_file)
    print('Output:     ' + merged_bam_prefix + '.nodup.bam')

def get_fastq(bam_file, sample_id, output_prefix, qsub_log_dir):
    mk_dir(os.path.dirname(output_prefix))
    command  = JAVA + ' -Xmx8G -Xms8G'
    command += ' -jar ' + PICARD_JAR
    command += ' SamToFastq'
    command += ' I='  + bam_file
    command += ' F='  + output_prefix + '_read1.fastq'
    command += ' F2=' + output_prefix + '_read2.fastq'
    command += ' NON_PF=true RE_REVERSE=true || true;'

    job_name_base = 'extract_fastq'
    job_name      = job_name_base + '_' + sample_id
    hold_job_name = 'extract_mt_and_unmapped_' + sample_id
    qsub_dir = qsub_log_dir + '/' + job_name_base
    memory_size = '12G'
    options = ['-hold_jid ' + hold_job_name]
    qsub(command, qsub_dir, job_name, memory_size, options)

    print('Job Name:   ' + job_name)
    print('H-job Name: ' + hold_job_name)
    print('Input:      ' + bam_file)
    print('Output:     ' + output_prefix + '_read[1/2].fastq')

def bwa_mem_mapping(fastq_files, bwa_index, sample_id, output_prefix, qsub_log_dir):
    mk_dir(os.path.dirname(output_prefix))
    num_threads = 10
    command  = ''
    command += BWA_MEM + ' mem -t ' + str(num_threads - 1)
    command += ' -R "@RG\tID:1\tSM:' + sample_id + '\tPL:ILLUMINA"'
    command += ' -K 10000000'
    command += ' ' + bwa_index
    command += ' ' + ' '.join(fastq_files)
    command += ' | ' + SAMTOOLS + ' view -Sb - > ' + output_prefix + '.bam;'
    command += SAMTOOLS + ' sort --threads ' + str(num_threads) + ' -o ' + output_prefix + '.sorted.bam' + ' ' + output_prefix + '.bam;'
    command += SAMTOOLS + ' index ' + output_prefix + '.sorted.bam;'

    command += ' rm -f ' + output_prefix + '.bam;'

    job_name_base = 'bwa_mem'
    job_name      = job_name_base + '_' + sample_id
    hold_job_name = 'extract_fastq_' + sample_id
    qsub_dir = qsub_log_dir + '/' + job_name_base
    memory_size = '4G'
    options = ['-pe def_slot ' + str(num_threads), '-hold_jid ' + hold_job_name]
    qsub(command, qsub_dir, job_name, memory_size, options)

    print('Job Name:   ' + job_name)
    print('H-job Name: ' + hold_job_name)
    print('Input:      ' + ', '.join(fastq_files))
    print('Output:     ' + output_prefix + '.sorted.bam')

def variant_call(bam_file, reference, sample_id, output_file, qsub_log_dir):
    mk_dir(os.path.dirname(output_file))
    command  = JAVA + ' -Xmx8G -Xms8G'
    command += ' -jar ' + GATK_JAR
    command += ' -R ' + reference
    command += ' -T HaplotypeCaller'
    command += ' -I ' + bam_file
    command += ' --emitRefConfidence GVCF'
    command += ' -variant_index_type LINEAR'
    command += ' -variant_index_parameter 128000'
    command += ' -ploidy 1'
    command += ' -L org_data/targets.interval_list'
    command += ' -o ' + output_file
    command += ' --heterozygosity 0.0025'
    command += ' --indel_heterozygosity 0.0003125'
    job_name_base = 'variant_call'
    job_name      = job_name_base + '_' + sample_id
    hold_job_name = 'bwa_mem_' + sample_id
    qsub_dir = qsub_log_dir + '/' + job_name_base
    memory_size = '12G'
    options = ['-hold_jid ' + hold_job_name]
    qsub(command, qsub_dir, job_name_base + '_' + sample_id, memory_size, options)

    print('Job Name:   ' + job_name)
    print('H-job Name: ' + hold_job_name)
    print('Input:      ' + bam_file)
    print('Output:     ' + output_file)

def merge_individual_vcf(vcf_file1, vcf_file2, sample_id, mt_size, shift_size, output_file, qsub_log_dir):
    mk_dir(os.path.dirname(output_file))
    margin = 2500

    command  = 'python scripts/merge_vcf_files.py'
    command += ' --vcf-file1 ' + vcf_file1
    command += ' --vcf-file2 ' + vcf_file2
    command += ' --mt-size ' + str(mt_size)
    command += ' --base-shift ' + str(shift_size)
    command += ' --margin ' + str(margin)
    command += ' --output-file ' + output_file + ';'

    command += 'rm -f ' + output_file + '.gz;'
    command += 'rm -f ' + output_file + '.gz.tbi;'
    command += 'bgzip ' + output_file + ';'
    command += 'tabix -p vcf ' + output_file + '.gz;'

    job_name_base = 'merge_vcf'
    job_name      = job_name_base + '_' + sample_id
    hold_job_name = 'variant_call_' + sample_id
    qsub_dir = qsub_log_dir + '/' + job_name_base
    options = ['-hold_jid ' + hold_job_name]
    memory_size = '2G'
    qsub(command, qsub_dir, job_name_base + '_' + sample_id, memory_size, options)

    print('Job Name:   ' + job_name)
    print('H-job Name: ' + hold_job_name)
    print('Input:      ' + ', '.join([vcf_file1, vcf_file2]))
    print('Output:     ' + output_file)

def qsub(command, qsub_dir, job_name, memory_size, additional_options):
    mk_dir(qsub_dir + '/scripts')
    mk_dir(qsub_dir + '/log')
    options = ['-V', '-cwd',
        '-N ' + job_name, '-o ' + qsub_dir + '/log', '-e ' + qsub_dir + '/log',
        '-l s_vmem=' + str(memory_size) + ',mem_req=' + str(memory_size)]
    if additional_options is not None:
        options.extend(additional_options)
    options.append('-soft -l sjob,lmem')

    job_id_file = './log/job_id.txt'
    mk_dir(os.path.dirname(job_id_file))

    print('output dir: ' + qsub_dir + '/scripts')
    print('options: ' + ' '.join(options))

    script_file  = qsub_dir + '/scripts/' + job_name
    script_file += '_' + str(int(time.time()))
    script_file += '_' + str(os.getpid())
    script_file += '.sh'

    print('create ' + script_file)
    gen_qsub_script_file(command, script_file)

    print('qsub ' + ' '.join(options) + ' ' + script_file)

    command = 'qsub -terse ' + ' '.join(options) + ' ' + script_file
    print('command: ' + command)

    qsub_result = my_subprocess(command)
    if qsub_result.wait() == 0:
        with open(job_id_file, 'a') as fp:
            for line in qsub_result.stdout:
                fp.write(line + '\n')
                print(line)
    else:
        sys.stderr.write('Failed to apply command ' + command + '\n')

def gen_qsub_script_file(command, script_file):
    mk_dir(os.path.dirname(script_file))

    with open(script_file, 'w') as fp:
        fp.write('#!/bin/bash' + '\n')
        fp.write('#$ -S /bin/bash' + '\n')
        fp.write('\n')
        fp.write('set -euxo pipefail\n')
        fp.write('trap \'exit 100\' ERR SIGXCPU\n')
        fp.write('\n')
        fp.write(command + '\n')
    execute('chmod 0755 ' + script_file)

def main_analysis(bam_file, sample_id, fasta_setting, output_dir, qsub_log_dir):

    extracted_bam_prefix = '/'.join([output_dir, sample_id, 'extracted', sample_id])

    #extract_mt_and_unmapped(bam_file, sample_id, extracted_bam_prefix, qsub_log_dir)

    #fastq_prefix = '/'.join([output_dir, sample_id, 'extracted', sample_id])
    #get_fastq(extracted_bam_prefix + '_chrM_and_unmapped.nodup.bam', sample_id, fastq_prefix, qsub_log_dir)

    #fastq_files = []
    #fastq_files.append(fastq_prefix + '_read1.fastq')
    #fastq_files.append(fastq_prefix + '_read2.fastq')

    reference = fasta_setting.reference_prefix + '_with_mtDNA.fa'

    bam_prefix = '/'.join([output_dir, sample_id, 'bwa_mem', sample_id])
    #bwa_mem_mapping(fastq_files, reference, sample_id, bam_prefix, qsub_log_dir)
    vcf_file = '/'.join([output_dir, sample_id, 'vcf', sample_id + '.het2.5.vcf'])
    variant_call(bam_prefix + '.sorted.bam', reference, sample_id, vcf_file, qsub_log_dir)

    reference = fasta_setting.reference_prefix + '_with_shifted_mtDNA.fa'
    bam_prefix_with_shifted_mtDNA = '/'.join([output_dir, sample_id, 'bwa_mem', sample_id + '_shifted'])
    #bwa_mem_mapping(fastq_files, reference, sample_id, bam_prefix_with_shifted_mtDNA, qsub_log_dir)
    vcf_file_with_shifted_mtDNA   = '/'.join([output_dir, sample_id, 'vcf', sample_id + '_shifted.het2.5.vcf'])
    variant_call(bam_prefix_with_shifted_mtDNA + '.sorted.bam', reference, sample_id, vcf_file_with_shifted_mtDNA, qsub_log_dir)

    merged_vcf_file = '/'.join([output_dir, sample_id, 'vcf', sample_id + '_merged.het2.5.vcf'])
    merge_individual_vcf(vcf_file, vcf_file_with_shifted_mtDNA, sample_id, fasta_setting.mtDNA_size, fasta_setting.shift_size, merged_vcf_file, qsub_log_dir)

def get_data_item_list(sample_data_file):
    data_item_list = []
    with open(sample_data_file) as fp:
        for line in fp:
            break
        for line in fp:
            items = line.strip().split()
            project_id = items[3]   # <- tmmid
            bam_file = '/u4/share1/home/gpc-gr/panel-build37/release/individual/{tmmid}/{tmmid}.bwamem.bam'.format(tmmid=project_id)

            data_item = Data_item(bam_file, project_id)
            data_item_list.append(data_item)
    return data_item_list

def main(argv):
    org_reference    = 'org_data/hs37d5.fa'
    org_mtDNA_fasta  = 'org_data/rCRS.fa'
    reference_prefix = 'data/fasta/hs37d5'
    shift_size       = 10000

    fasta_setting = Fasta_setting(org_reference, org_mtDNA_fasta, reference_prefix)
    mtDNA_size = get_contig_size(org_mtDNA_fasta, 'MT')
    fasta_setting.mtDNA_size = mtDNA_size
    fasta_setting.shift_size = shift_size

    #reference_preparation(fasta_setting)

    sample_data_file = '/home/gpc-gr/panel-build37/work/11-single_sample/idtable-current.tsv'
    data_item_list = get_data_item_list(sample_data_file)

    output_dir = 'debug_analysis'
    qsub_dir   = 'qsub_debug'

    for i in xrange(len(data_item_list)):
        data_item = data_item_list[i]
        project_id = data_item.project_id
        bam_file   = data_item.bam_file

        merged_vcf_file = 'debug_analysis/{tmmid}/vcf/{tmmid}_merged.het2.5.vcf.gz'.format(tmmid=project_id)

        if (0 <= i < 4100) and (not os.path.exists(merged_vcf_file)):
            main_analysis(bam_file, project_id, fasta_setting, output_dir, qsub_dir)
            import time; time.sleep(1)


if __name__ == '__main__':
    argv = sys.argv
    main(argv)
