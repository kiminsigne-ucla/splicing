'''
@author Kimberly Insigne
March 28, 2017

This scripts update the sequence information in the splicemod reference file 
containing all sequences, concatenated from the original mutant and natural 
reference files. It uses the genomic coordinates to pull the sequence from the
genome and reverse complement if on the negative strand. The sequence currently
in the reference file may refer to the sequence synthesized on the chip, in which
case it is the version with the lowest A's (easier to synthesize). It could also
be the sequence without strand consideration. To be sure, let's just pull from 
the genome.
'''

import MySQLdb
# use cogent to extract sequences from the genome
from cogent.db.ensembl import HostAccount, Genome
import pandas as pd

# set up connection
# host, user, password, port
pycog = HostAccount('ensembldb.ensembl.org', 'anonymous', '', 3306)
hs37 = Genome('human', Release=78, account=pycog)

# read in reference
data = pd.read_table('../produced_data/splicemod_data_clean.txt', sep = '\t')
# reformat so intron/exon length columns are int and not float, convert NA to 0 so we can use int
data[['intron1_len', 'exon_len', 'intron2_len']] = data[['intron1_len', 'exon_len', 'intron2_len']].fillna(0.0).astype(int)

# grab sequences with chr, start and end information
lib = data[data.chr.notnull() & data.start.notnull() & data.end.notnull()]


def grab_region(chr, start, end, strand, genome):

	start = int(start)
	end = int(end)

	cog_reg = genome.getRegion(CoordName=chr, Start=start, End=end)
	cog_seq = cog_reg.Seq

	if strand == -1:
		cog_seq = cog_seq.rc()
	seq = cog_seq._seq

	return seq

lib['updated_seq'] = [grab_region(lib.chr.iloc[i],
								  lib.start.iloc[i],
								  lib.end.iloc[i],
								  lib.strand.iloc[i],
								  hs37) for i in range(len(lib))]

lib.to_csv('splicemod_data_clean_updated_ref.txt', sep='\t', index=False)


