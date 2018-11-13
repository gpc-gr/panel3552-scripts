#
#
#

# ================================================================================
#
# ================================================================================

ROOT				:= $(shell pwd)
TMMID				=
FASTQ1				=
FASTQ2				=
SEX				=
FASTQ1_R1			= $(wildcard $(FASTQ1)_R1*.fastq.gz)
FASTQ1_R2			= $(wildcard $(FASTQ1)_R2*.fastq.gz)
FASTQ2_R1			= $(wildcard $(FASTQ2)_R1*.fastq.gz)
FASTQ2_R2			= $(wildcard $(FASTQ2)_R2*.fastq.gz)
INCLUDE_NORMAL_BAM_METRICS	:= false
INCLUDE_NORMAL_GVCF		:= false
INCLUDE_NORMAL_GVCF_METRICS	:= false
INCLUDE_NORMAL_VCF		:= false
INCLUDE_NORMAL_VCF_METRICS	:= false
INCLUDE_BQSR_BAM_METRICS	:= false
INCLUDE_BQSR_GVCF		:= false
INCLUDE_BQSR_GVCF_METRICS	:= false
INCLUDE_BQSR_VCF		:= false
INCLUDE_BQSR_VCF_METRICS	:= false
INCLUDE_BQSR_PLOT		:= false
DRY_RUN				:= false
SYNC				:= false

JOB_PREFIX			:= gpcgr

GPC_GR_HOME			:= /share2/home/gpc-gr
WORK_ROOT			:= $(GPC_GR_HOME)/panel-build37/work
MODULE_PATH			:= $(GPC_GR_HOME)/local/modulefiles
REFERENCE_FASTA			:= $(WORK_ROOT)/00-reference/hs37d5.fa
REFERENCE_FASTA_PAR3		:= $(WORK_ROOT)/00-reference/hs37d5_PAR3.fa
DBSNP_VCF			:= $(WORK_ROOT)/03-gatk_resource_bundle/bqsr/dbsnp_138.b37.vcf.gz
MILLS_INDEL_VCF			:= $(WORK_ROOT)/03-gatk_resource_bundle/bqsr/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
ONEK_INDEL_VCF			:= $(WORK_ROOT)/03-gatk_resource_bundle/bqsr/1000G_phase1.indels.b37.vcf.gz

TQSUB				:= /u4/share1/home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/tqsub
TDIR				:= $(realpath $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))/templates
TCOMMON				= '{"job_prefix": "$(JOB_PREFIX)", "module_path": "$(MODULE_PATH)"}'
TCOMMON				+= '{"reference_fasta": "$(REFERENCE_FASTA)", "reference_fasta_par3": "$(REFERENCE_FASTA_PAR3)"}'
TCOMMON				+= '{"root": "$(ROOT)", "tmmid": "$(TMMID)", "sex": "$(SEX)"}'

ifeq ($(DRY_RUN),true)
TQSUB				+= --dry-run
endif
ifeq ($(SYNC),true)
TQSUB				+= --sync
endif


# ================================================================================
# preparation
# ================================================================================

prepare:\
	$(ROOT)/logs


$(ROOT)/logs:
	@echo "[INFO] Preparing $(ROOT)/logs"
	@mkdir -p $(ROOT)/logs


# ================================================================================
# reseq
# ================================================================================

PRODUCTS = $(TMMID).bwamem.bam

ifeq ($(INCLUDE_NORMAL_BAM_METRICS),true)
PRODUCTS += $(TMMID).bwamem.alignment_summary_metrics
PRODUCTS += $(TMMID).bwamem.idxstats
PRODUCTS += $(TMMID).bwamem.wgs_metrics_autosome
endif
ifeq ($(INCLUDE_NORMAL_GVCF),true)
PRODUCTS += $(TMMID).bwamem.hc3.g.vcf.gz
PRODUCTS += $(TMMID).bwamem.hc3.chrXY_PAR2.g.vcf.gz
PRODUCTS += $(TMMID).bwamem.chrXY_PAR3.bam
PRODUCTS += $(TMMID).bwamem.hc3.chrXY_PAR3.g.vcf.gz
endif
#ifeq ($(INCLUDE_NORMAL_VCF),true)
#PRODUCTS += $(TMMID).bwamem.hc3.vcf.gz
#endif
ifeq ($(INCLUDE_BQSR_BAM_METRICS),true)
PRODUCTS += $(TMMID).bwamem.bqsr.alignment_summary_metrics
endif
ifeq ($(INCLUDE_BQSR_GVCF),true)
PRODUCTS += $(TMMID).bwamem.bqsr.hc3.g.vcf.gz
endif
#ifeq ($(INCLUDE_BQSR_VCF),true)
#PRODUCTS += $(TMMID).bwamem.bqsr.hc3.vcf.gz
#endif
ifeq ($(INCLUDE_BQSR_PLOT),true)
PRODUCTS += $(TMMID).bwamem.bqsr.pdf
endif


.SECONDARY: $(FASTQ1).bwamem.bam $(FASTQ2).bwamem.bam


.DEFAULT: reseq
reseq: $(PRODUCTS)


# ================================================================================
# mapping (s1x)
# ================================================================================

$(FASTQ1).bwamem.bam: $(FASTQ1_R1) $(FASTQ1_R2)
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ1).sh"
	@$(TQSUB)\
		$(TDIR)/s11-bwa_mem_align.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ1).sh\
		$(TCOMMON)\
		'{"fastq_prefix": "$(FASTQ1)"}'\
		'{"fastq_r1": "$(FASTQ1_R1)"}'\
		'{"fastq_r2": "$(FASTQ1_R2)"}'\
		'{"output": "$(FASTQ1).bwamem.bam"}'


$(FASTQ2).bwamem.bam: $(FASTQ2_R1) $(FASTQ2_R2)
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ2).sh"
	@$(TQSUB)\
		$(TDIR)/s11-bwa_mem_align.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ2).sh\
		$(TCOMMON)\
		'{"fastq_prefix": "$(FASTQ2)"}'\
		'{"fastq_r1": "$(FASTQ2_R1)"}'\
		'{"fastq_r2": "$(FASTQ2_R2)"}'\
		'{"output": "$(FASTQ2).bwamem.bam"}'


$(TMMID).bwamem.bam: $(FASTQ1).bwamem.bam $(FASTQ2).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh"
	@$(TQSUB)\
		$(TDIR)/s12-rmdup.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh\
		$(TCOMMON)\
		'{"prev_script_names": ["$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ1).sh", "$(JOB_PREFIX)-$(TMMID)-s11-bwa_mem_align-$(FASTQ2).sh"]}'\
		'{"fastq1_prefix": "$(FASTQ1)"}'\
		'{"fastq2_prefix": "$(FASTQ2)"}'\
		'{"source1": "$(FASTQ1).bwamem.bam"}'\
		'{"source2": "$(FASTQ2).bwamem.bam"}'\
		'{"output": "$(TMMID).bwamem.bam"}'\
		'{"metrics_output": "$(TMMID).bwamem.rmdup_metrics.txt"}'


$(TMMID).bwamem.alignment_summary_metrics: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s13-bam_metrics.sh"
	@$(TQSUB)\
		$(TDIR)/s13-bam_metrics.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s13-bam_metrics.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh"}'\
		'{"source": "$(TMMID).bwamem.bam"}'\
		'{"output_prefix": "$(TMMID).bwamem"}'


$(TMMID).bwamem.idxstats: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s14-bam_samtools_metrics.sh"
	@$(TQSUB)\
		$(TDIR)/s14-bam_samtools_metrics.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s14-bam_samtools_metrics.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s13-bam_metrics.sh"}'\
		'{"source": "$(TMMID).bwamem.bam"}'\
		'{"output_prefix": "$(TMMID).bwamem"}'


$(TMMID).bwamem.wgs_metrics_autosome: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s17-picard_wgs_metrics.sh"
	@$(TQSUB)\
		$(TDIR)/s17-picard_wgs_metrics.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s17-picard_wgs_metrics.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh"}'\
		'{"source": "$(TMMID).bwamem.bam"}'\
		'{"output_prefix": "$(TMMID).bwamem.wgs_metrics"}'


$(TMMID).chrXY.interleaved.fastq.gz: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s15-extract_xy_unmapped.sh"
	@$(TQSUB)\
		$(TDIR)/s15-extract_xy_unmapped.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s15-extract_xy_unmapped.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh"}'\
		'{"source": "$(TMMID).bwamem.bam"}'\
		'{"output": "$(TMMID).chrXY.interleaved.fastq.gz"}'


$(TMMID).bwamem.chrXY_PAR3.bam: $(TMMID).chrXY.interleaved.fastq.gz
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s16-bwa_mem_align_PAR3.sh"
	@$(TQSUB)\
		$(TDIR)/s16-bwa_mem_align_PAR3.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s16-bwa_mem_align_PAR3.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s15-extract_xy_unmapped.sh"}'\
		'{"fastq_interleaved": "$(TMMID).chrXY.interleaved.fastq.gz"}'\
		'{"output": "$(TMMID).bwamem.chrXY_PAR3.bam"}'\
		'{"metrics_output": "$(TMMID).bwamem.chrXY_PAR3.rmdup_metrics.txt"}'


# ================================================================================
# BQSR (s2x)
# ================================================================================

T_BQSR_SITES = '{"dbsnp": "$(DBSNP_VCF)", "mills_indel": "$(MILLS_INDEL_VCF)", "onek_indel": "$(ONEK_INDEL_VCF)"}'
T_BQSR_TABLES = '{"before_table": "$(TMMID).bwamem.bqsr.table", "after_table": "$(TMMID).bwamem.bqsr_after.table"}'


$(TMMID).bwamem.bqsr.table: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s21-bqsr_table.sh"
	@$(TQSUB)\
		$(TDIR)/s21-bqsr_table.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s21-bqsr_table.sh\
		$(TCOMMON)\
		$(T_BQSR_SITES)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh"}'\
		'{"source_bam": "$(TMMID).bwamem.bam"}'\
		'{"output_table": "$(TMMID).bwamem.bqsr.table"}'


$(TMMID).bwamem.bqsr_after.table: $(TMMID).bwamem.bam $(TMMID).bwamem.bqsr.table
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s22-bqsr_after_table.sh"
	@$(TQSUB)\
		$(TDIR)/s22-bqsr_after_table.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s22-bqsr_after_table.sh\
		$(TCOMMON)\
		$(T_BQSR_SITES)\
		$(T_BQSR_TABLES)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s21-bqsr_table.sh"}'\
		'{"bam": "$(TMMID).bwamem.bam"}'


$(TMMID).bwamem.bqsr.pdf: $(TMMID).bwamem.bqsr.table $(TMMID).bwamem.bqsr_after.table
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s23-bqsr_plot.sh"
	@$(TQSUB)\
		$(TDIR)/s23-bqsr_plot.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s23-bqsr_plot.sh\
		$(TCOMMON)\
		$(T_BQSR_TABLES)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s22-bqsr_after_table.sh"}'\
		'{"plot": "$(TMMID).bwamem.bqsr.pdf"}'


$(TMMID).bwamem.bqsr.bam: $(TMMID).bwamem.bam $(TMMID).bwamem.bqsr.table
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s24-bqsr.sh"
	@$(TQSUB)\
		$(TDIR)/s24-bqsr.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s24-bqsr.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s21-bqsr_table.sh"}'\
		'{"source_bam": "$(TMMID).bwamem.bam"}'\
		'{"output_bam": "$(TMMID).bwamem.bqsr.bam"}'\
		'{"bqsr_table": "$(TMMID).bwamem.bqsr.table"}'


$(TMMID).bwamem.bqsr.alignment_summary_metrics: $(TMMID).bqsr.bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s25-bqsr_bam_metrics.sh"
	@$(TQSUB)\
		$(TDIR)/s25-bqsr_bam_metrics.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s25-bqsr_bam_metrics.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s24-bqsr.sh"}'\
		'{"source": "$(TMMID).bqsr.bwamem.bam"}'\
		'{"output_prefix": "$(TMMID).bwamem.bqsr"}'


# ================================================================================
# variant call (s3x, s4x)
# ================================================================================

$(TMMID).bwamem.hc3.g.vcf.gz: $(TMMID).bwamem.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s31-gvcf.sh"
	@$(TQSUB)\
		$(TDIR)/s31-gvcf.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s31-gvcf.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s14-bam_samtools_metrics.sh"}'\
		'{"source": "$(TMMID).bwamem.bam"}'\
		'{"output": "$(TMMID).bwamem.hc3.g.vcf.gz"}'


$(TMMID).bwamem.hc3.chrXY_PAR2.g.vcf.gz: $(TMMID).bwamem.bam $(TMMID).bwamem.hc3.g.vcf.gz
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s32-gvcf_chrXY_PAR2.sh"
	@$(TQSUB)\
		$(TDIR)/s32-gvcf_chrXY_PAR2.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s32-gvcf_chrXY_PAR2.sh\
		$(TCOMMON)\
		'{"prev_script_names": ["$(JOB_PREFIX)-$(TMMID)-s12-rmdup.sh", "$(JOB_PREFIX)-$(TMMID)-s31-gvcf.sh"]}'\
		'{"source_bam": "$(TMMID).bwamem.bam"}'\
		'{"source_gvcf": "$(TMMID).bwamem.hc3.g.vcf.gz"}'\
		'{"output": "$(TMMID).bwamem.hc3.chrXY_PAR2.g.vcf.gz"}'


$(TMMID).bwamem.hc3.chrXY_PAR3.g.vcf.gz: $(TMMID).bwamem.chrXY_PAR3.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s33-gvcf_chrXY_PAR3.sh"
	@$(TQSUB)\
		$(TDIR)/s33-gvcf_chrXY_PAR3.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s33-gvcf_chrXY_PAR3.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s16-bwa_mem_align_PAR3.sh"}'\
		'{"source": "$(TMMID).bwamem.chrXY_PAR3.bam"}'\
		'{"output": "$(TMMID).bwamem.hc3.chrXY_PAR3.g.vcf.gz"}'


$(TMMID).bwamem.bqsr.hc3.g.vcf.gz: $(TMMID).bwamem.bqsr.bam
	@echo "[INFO] Preparing $(ROOT)/logs/$(TMMID)-s41-bqsr_gvcf.sh"
	@$(TQSUB)\
		$(TDIR)/s41-bqsr_gvcf.sh.j2\
		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s41-bqsr_gvcf.sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(TMMID)-s24-bqsr.sh"}'\
		'{"source": "$(TMMID).bwamem.bqsr.bam"}'\
		'{"output": "$(TMMID).bwamem.bqsr.hc3.g.vcf.gz"}'


#$(TMMID).bwamem.bqsr.hc3.vcf.gz: $(TMMID).bwamem.bqsr.bam
#	@echo "[INFO] Preparing $(ROOT)/logs/$(TMMID)-s42-bqsr_vcf.sh"
#	@$(TQSUB)\
#		$(TDIR)/s42-bqsr_vcf.sh.j2\
#		$(ROOT)/logs/$(JOB_PREFIX)-$(TMMID)-s42-bqsr_vcf.sh\
#		$(TCOMMON)\
#		'{"source": "$(TMMID).bwamem.bqsr.bam", "output": "$(TMMID).bwamem.bqsr.hc3.vcf.gz"}'

# ================================================================================
# clean
# ================================================================================

clean:
	rm $(FASTQ1).bam $(FASTQ2).bam

