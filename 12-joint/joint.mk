#
#
#

# ================================================================================
#
# ================================================================================

ROOT				= $(shell pwd)
PANEL_ID			:= $(shell basename $(shell pwd))
REFERENCE_FASTA			:= /u4/share1/home/gpc-gr/panel-build37/work/00-reference/hs37d5.fa
GATK_RESOURCE_BUNDLE_PREFIX	:= /u4/share1/home/gpc-gr/panel-build37/work/03-gatk_resource_bundle/vqsr
CHROMOSOMES			:= 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
VQSR_SNV_THRESHOLD		:= 99.0
VQSR_INDEL_THRESHOLD		:= 97.0

JOB_PREFIX			:= gpcgr_joint_merge

GPC_GR_HOME			:= /u4/share1/home/gpc-gr
MODULE_PATH			:= $(GPC_GR_HOME)/local/modulefiles

TQSUB				:= /home/gpc-gr/panel-build37/work/progenv/pyenv-3.6/bin/tqsub --ignore-if-already-queued
TDIR				:= $(realpath $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))/templates
TCOMMON				:= '{"job_prefix": "$(JOB_PREFIX)", "module_path": "$(MODULE_PATH)", "reference_fasta": "$(REFERENCE_FASTA)", "resource_prefix": "$(GATK_RESOURCE_BUNDLE_PREFIX)", "root": "$(ROOT)"}'


# ================================================================================
#
# ================================================================================

RAW_VCFS		= $(addsuffix .vcf.gz,$(addprefix 02-merge/$(PANEL_ID)-full-,$(CHROMOSOMES)))
VQSR_FULL_VCFS		= $(addsuffix -vqsr_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD).vcf.gz,$(addprefix 03-vqsr/$(PANEL_ID)-full-,$(CHROMOSOMES)))
VQSR_PASS_VCFS		= $(addsuffix .vcf.gz,$(addprefix 03-vqsr/$(PANEL_ID)-passed_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)-,$(CHROMOSOMES)))

#RAW_VCF_ANNOVARS	= $(addsuffix .annovar.gz,$(RAW_VCFS))
RAW_VCF_ANNOVARS	= $(addsuffix .annovar.all.gz,$(RAW_VCFS))

PRODUCTS		= $(RAW_VCFS) $(VQSR_FULL_VCFS) $(VQSR_PASS_VCFS) $(RAW_VCF_ANNOVARS)


.PHONY: all
all: $(PRODUCTS)


# ================================================================================
# 01-joint_genotyping
# ================================================================================

# run manually


# ================================================================================
# 02-merge
# ================================================================================

02-merge/$(PANEL_ID)-disconcordant_variants.tsv:
	@$(TQSUB)\
		$(TDIR)/s02-check_chunk_concordance.sh.j2\
		$(ROOT)/02-merge/logs/$(JOB_PREFIX)-$(PANEL_ID)-s02-check_chunk_concordance.sh\
		$(TCOMMON)\
		'{"chunk_root": "$(ROOT)/01-chunk"}'\
		'{"output": "$(ROOT)/02-merge/$(PANEL_ID)-disconcordant_variants.tsv"}'


define CHROMOSOME_VCF
02-merge/$(PANEL_ID)-full-$(1).vcf.gz: 02-merge/$(PANEL_ID)-disconcordant_variants.tsv
	@$(TQSUB)\
		$(TDIR)/s03-merge_chunk.sh.j2\
		$(ROOT)/02-merge/logs/$(JOB_PREFIX)-$(PANEL_ID)-s03-merge_chunk-$(1).sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(PANEL_ID)-s02-check_chunk_concordance.sh"}'\
		'{"sources": "$(ROOT)/01-chunk/$(PANEL_ID)-$(1)_*.vcf.gz"}'\
		'{"output": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz"}'\
		'{"region": "$(1)"}'

02-merge/$(PANEL_ID)-full-$(1).vcf.gz.annovar.gz: 02-merge/$(PANEL_ID)-full-$(1).vcf.gz
	@$(TQSUB)\
		$(TDIR)/s11-annovar_raw.sh.j2\
		$(ROOT)/02-merge/logs/$(JOB_PREFIX)-$(PANEL_ID)-s11-annovar_raw-$(1).sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(PANEL_ID)-s03-merge_chunk-$(1).sh"}'\
		'{"source": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz"}'\
		'{"output": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz.annovar.gz"}'

02-merge/$(PANEL_ID)-full-$(1).vcf.gz.annovar.all.gz: 02-merge/$(PANEL_ID)-full-$(1).vcf.gz
	@$(TQSUB)\
		$(TDIR)/s12-annovar_raw_all_samples.sh.j2\
		$(ROOT)/02-merge/logs/$(JOB_PREFIX)-$(PANEL_ID)-s12-annovar_raw_all_samples-$(1).sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(PANEL_ID)-s03-merge_chunk-$(1).sh"}'\
		'{"source": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz"}'\
		'{"output": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz.annovar.all.gz"}'
endef

$(foreach chromosome,$(CHROMOSOMES),$(eval $(call CHROMOSOME_VCF,$(chromosome))))


# ================================================================================
# 03-vqsr
# ================================================================================

FULL_CHROM_VCF_MERGE_SCRIPTS = $(addsuffix .sh,$(addprefix $(JOB_PREFIX)-$(PANEL_ID)-s03-merge_chunk-,$(CHROMOSOMES)))

03-vqsr/$(PANEL_ID)-full-vqsr_snv.recal: $(RAW_VCFS)
	@$(TQSUB)\
		$(TDIR)/s04-vqsr_snv.sh.j2\
		$(ROOT)/03-vqsr/logs/$(JOB_PREFIX)-$(PANEL_ID)-s04-vqsr_snv.sh\
		$(TCOMMON)\
		'{"prev_script_names": "$(FULL_CHROM_VCF_MERGE_SCRIPTS)"}'\
		'{"sources": "$(addprefix $(ROOT)/,$(RAW_VCFS))"}'\
		'{"output_prefix": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_snv"}'


03-vqsr/$(PANEL_ID)-full-vqsr_indel.recal: $(RAW_VCFS)
	@$(TQSUB)\
		$(TDIR)/s04-vqsr_indel.sh.j2\
		$(ROOT)/03-vqsr/logs/$(JOB_PREFIX)-$(PANEL_ID)-s04-vqsr_indel.sh\
		$(TCOMMON)\
		'{"prev_script_names": "$(FULL_CHROM_VCF_MERGE_SCRIPTS)"}'\
		'{"sources": "$(addprefix $(ROOT)/,$(RAW_VCFS))"}'\
		'{"output_prefix": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_indel"}'


define CHROMOSOME_VQSR
03-vqsr/$(PANEL_ID)-full-$(1)-vqsr_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD).vcf.gz:\
		02-merge/$(PANEL_ID)-full-$(1).vcf.gz\
		03-vqsr/$(PANEL_ID)-full-vqsr_snv.recal\
		03-vqsr/$(PANEL_ID)-full-vqsr_indel.recal

	@$(TQSUB)\
		$(TDIR)/s05-apply_vqsr.sh.j2\
		$(ROOT)/03-vqsr/logs/$(JOB_PREFIX)-$(PANEL_ID)-s05-apply_vqsr-$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)_$(1).sh\
		$(TCOMMON)\
		'{"prev_script_names": ["$(JOB_PREFIX)-$(PANEL_ID)-s04-vqsr_snv.sh", "$(JOB_PREFIX)-$(PANEL_ID)-s04-vqsr_indel.sh"]}'\
		'{"region": "$(1)"}'\
		'{"source_vcf": "$(ROOT)/02-merge/$(PANEL_ID)-full-$(1).vcf.gz"}'\
		'{"snv_recal": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_snv.recal"}'\
		'{"snv_tranches": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_snv.tranches"}'\
		'{"snv_filter_level": "$(VQSR_SNV_THRESHOLD)"}'\
		'{"indel_recal": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_indel.recal"}'\
		'{"indel_tranches": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-vqsr_indel.tranches"}'\
		'{"indel_filter_level": "$(VQSR_INDEL_THRESHOLD)"}'\
		'{"output_vcf": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-$(1)-vqsr_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD).vcf.gz"}'


03-vqsr/$(PANEL_ID)-passed_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)-$(1).vcf.gz:\
		03-vqsr/$(PANEL_ID)-full-$(1)-vqsr_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD).vcf.gz

	@$(TQSUB)\
		$(TDIR)/s06-filter_passed.sh.j2\
		$(ROOT)/03-vqsr/logs/$(JOB_PREFIX)-$(PANEL_ID)-s06-filter_passed-$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)_$(1).sh\
		$(TCOMMON)\
		'{"prev_script_name": "$(JOB_PREFIX)-$(PANEL_ID)-s05-apply_vqsr-$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)_$(1).sh"}'\
		'{"region": "$(1)"}'\
		'{"source": "$(ROOT)/03-vqsr/$(PANEL_ID)-full-$(1)-vqsr_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD).vcf.gz"}'\
		'{"output": "$(ROOT)/03-vqsr/$(PANEL_ID)-passed_$(VQSR_SNV_THRESHOLD)_$(VQSR_INDEL_THRESHOLD)-$(1).vcf.gz"}'
endef

$(foreach chromosome,$(CHROMOSOMES),$(eval $(call CHROMOSOME_VQSR,$(chromosome))))

