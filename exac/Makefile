# @author Kimberly Insigne
# kiminsigne@gmail.com
#
# This Makefile merges paired read files, aligns each merged file to reference,
# and aggregates the alignments.
# 
# This is meant to be a general purpose template and the sample list should be
# specified and changed for the specific use. This uses the 09/2016 ExAC data as 
# an example.

# paths
SHELL := /bin/bash
SOURCE_DATA := /data/projects/splicing/raw_data/201609_exac
DATA := raw_data/201609_exac
SCRIPTS := ../scripts # path to master splicing scripts

# programs
BBMERGE := bbmerge.sh

# string processing to making dependencies easier
READS := $(wildcard $(DATA)/*.fastq*)
NO_GZ := $(READS:.gz=)
BASE := $(subst _R1,, $(filter-out %_R2, $(notdir $(NO_GZ:.fastq=))))

# -----------------------------------------------------------------------------
# Recipes

all: setup align

align : produced_data/all_alignments.csv

clean :
	@ rm -f produced_data/*

.PHONY : clean

# keep all intermediate files
.SECONDARY : 

setup :
	@# Actions:
	@#	Create symbolic links from original data folder to local data folder, 
	@#	simplify names for later dependecies. 
	@# Dependencies:
	@#	None. This rule will run every time, overwriting symlinks if they exist.
	@ echo "Creating data links... (will run everytime)"
	@ ../scripts/symlink_reformat_fastq_names.sh $(SOURCE_DATA) $(DATA)

# -----------------------------------------------------------------------------
# Primary processing pipeline

produced_data/%_merged.fastq : $(DATA)/%_R1.fastq.gz $(DATA)/%_R2.fastq.gz
	@# Actions: 
	@#	Merge compressed paired read FASTQ files with bbmerge
	@# Dependencies:
	@#	compressed fastq files in data directory
	$(BBMERGE) in1=$(word 1, $^) in2=$(word 2, $^) out=$@ \
		outu1=produced_data/$(subst .fastq.gz,_unmerged.fastq,$(notdir $(word 1, $^))) \
		outu2=produced_data/$(subst .fastq.gz,_unmerged.fastq,$(notdir $(word 2, $^)))


produced_data/%_merged.fastq : $(DATA)/%_R1.fastq $(DATA)/%_R2.fastq
	@# Actions: 
	@#	Merge uncompressed paired read FASTQ files with bbmerge
	@# Dependencies:
	@#	uncompressed fastq files in data directory
	$(BBMERGE) in1=$(word 1, $^) in2=$(word 2, $^) out=$@ \
		outu1=produced_data/$(subst .fastq,_unmerged.fastq,$(notdir $(word 1, $^))) \
		outu2=produced_data/$(subst .fastq,_unmerged.fastq,$(notdir $(word 2, $^)))

# merge : $(patsubst %,produced_data/%_merged.fastq,$(BASE))
# 	@# Actions:
# 	@#	Dummy rule that triggers merging of paired reads, use this as placeholder
# 	@#	for end of pipleine
# 	@# Dependencies:
# 	@#	merged fastq files

produced_data/%_perfect_alignment.txt : produced_data/%_merged.fastq $(SCRIPTS)/perfect_alignment.py
	@# Actions:
	@#	Align sequences to reference FASTA file
	@# Dependencies:
	@#	merged FASTQ files from previous step.
	@ echo "Aligning" $(subst _perfect_alignment.txt,, $(notdir $@))
	@ python $(SCRIPTS)/perfect_alignment.py -i $< -r ref/exac_ref_v2_all.fa -o $@


produced_data/all_alignments.csv : $(patsubst %,produced_data/%_perfect_alignment.txt,$(BASE)) $(SCRIPTS)/aggregate_alignments.py
	@# Actions:
	@#	Aggregate all perfect alignments into one file
	@# Dependencies:
	@#	perfect alignment files in produced_data/ from the previous step
	@ echo "Aggregating alignments..."
	@ ls produced_data/*perfect_alignment.txt > produced_data/alignment_list.txt
	@ python $(SCRIPTS)/aggregate_alignments.py produced_data/alignment_list.txt $@
	@ rm produced_data/alignment_list.txt

