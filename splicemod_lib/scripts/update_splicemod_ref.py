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

