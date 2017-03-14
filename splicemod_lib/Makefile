# @author : Kimberly Insigne
# 
# This is a pipeline for the splicemod sequencing data in the SMN1 intron
# backbone done in 06/2016

# paths
SHELL := /bin/bash
SOURCE_DATA := /data/projects/splicing/raw-data/201606_smn1_intron
DATA := raw_data/201606_smn1_intron

# programs
BBMERGE := bbmerge.sh

# string processing to making dependencies easier
READS := $(wildcard $(DATA)/*.fastq*)
NO_GZ := $(READS:.gz=)
BASE := $(subst _R1,, $(filter-out %_R2, $(notdir $(NO_GZ:.fastq=))))


all : produced_data/all_alignments.csv


.PHONY: setup


setup :
	@# Actions:
	@#	Create symbolic links from original data folder to local data folder, 
	@#	simplify names for later dependecies. 
	@# Dependencies:
	@#	None. This rule will run every time, overwriting symlinks if they exist.
	@ echo "Creating data links ... (will run everytime)"
	@ mkdir -p $(DATA)
	@ ../scripts/symlink_reformat_fastq_names.sh $(SOURCE_DATA) $(DATA)


produced_data/%_merged.fastq : $(DATA)/%_R1.fastq $(DATA)/%_R2.fastq
	@# Actions: 
	@#	Merge uncompressed paired read FASTQ files with bbmerge
	@# Dependencies:
	@#	uncompressed fastq files in data directory
	$(BBMERGE) in1=$(word 1, $^) in2=$(word 2, $^) out=$@ \
		outu1=produced_data/$(subst .fastq,_unmerged.fastq,$(notdir $(word 1, $^))) \
		outu2=produced_data/$(subst .fastq,_unmerged.fastq,$(notdir $(word 2, $^)))


produced_data/%_perfect_alignment.txt : produced_data/%_merged.fastq ../scripts/perfect_alignment.py
	@# Actions:
	@#	Align sequences to reference FASTA file
	@# Dependencies:
	@#	merged FASTQ files from previous step.
	@ echo "Aligning" $(subst _perfect_alignment.txt,, $(notdir $@))
	@ python ../scripts/perfect_alignment.py -i $< -r ref/splicemod_ref.fa -o $@


produced_data/all_alignments.csv : $(patsubst %,produced_data/%_perfect_alignment.txt,$(BASE)) ../scripts/aggregate_alignments.py
	@# Actions:
	@#	Aggregate all perfect alignments into one file
	@# Dependencies:
	@#	perfect alignment files in produced_data/ from the previous step
	@ echo "Aggregating alignments..."
	@ ls produced_data/*perfect_alignment.txt > produced_data/alignment_list.txt
	@ python ../scripts/aggregate_alignments.py produced_data/alignment_list.txt $@
	@ rm produced_data/alignment_list.txt


print-%  : ; @echo $* = $($*)