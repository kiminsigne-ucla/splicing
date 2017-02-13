"""
This script returns the average conservation (phastCons) score for a 
given interval region, as defined in the sorted input bed file.

"""

import argparse
import pandas as pd
import numpy as np


# wig functions courtesty of Chris Hartl
class WigReader(object):
    def __init__(self, file_):
        self._handle = gzip.open(file_) if file_.endswith('gz') else open(file_)
        self._header = self._handle.readline()
        self._header_info = dict([t.split('=') for t in self._header.strip().split(' ')[1:]])
        self._header_info = {k: (v if k == 'chrom' else int(v)) for k, v in self._header_info.items()}
        self._offset = 0

    @property
    def chrom(self):
        return self._header_info['chrom']

    def start(self):
        return self._header_info['start']

    def step(self):
        return self._header_info['step']

    def __iter__(self):
        for line in self._handle:
            if line.startswith('fixedStep'):
                self._header = line
                self._header_info = dict([t.split('=') for t in self._header.strip().split()[1:]])
                self._header_info = {k: (v if k == 'chrom' else int(v)) for k, v in self._header_info.items()}
                self._offset = 0
                continue
            data = (self.start() + self._offset * self.step(), float(line.strip()))
            self._offset += 1
            yield data


def in_interval(pos, interval):
    return interval[0] <= pos <= interval[1]


def get_interval(itr, interval):
    e = next(itr)
    while e[0] < interval[0]:
        e = next(itr)
    while in_interval(e[0], interval):
        yield e
        e = next(itr)
    # note: skips first base outside of interval, but ours are
    # not overlapping


def bed_iter(filename):
	"""
	Iterate through a bed file, yielding all intervals for one chromosome
	"""
	with open(filename) as bed:
		# initialize 
		chr_info = []
		interval = next(bed).strip().split('\t')
		# format: chr, start, end, name
		curr_chrom = interval[0]
		chr_info.append(interval)
		for line in bed:
			interval = line.strip().split('\t')
			chrom = interval[0]
			if chrom != curr_chrom:
				yield chr_info
				chr_info = []
				curr_chrom = chrom
			else:
				chr_info.append(interval)



def main(bed, cons_folder, output):

	outfile = open(output, 'w')

	bed_chroms = bed_iter(bed)
	cons_suffix = '.phastCons100way.wigFix'

	for chrom_region in bed_chroms:
		chrom = chrom_region[0][0]
		print "Working on chromosome ", str(chrom)
		# open appropriate wig file
		wig = WigReader(cons_folder + chrom + cons_suffix)
		for region in chrom_region:
			# outputs tuple, (position, score)
			cons_scores = get_interval(iter(wig), (int(region[1]), int(region[2])))
			# get average
			cons_mean = np.mean([x[1] for x in cons_scores])
			# add conservation mean to region information, output to file
			region.append(str(cons_mean))
			outfile.write('\t'.join(region)+'\n')

	outfile.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser('wrapper for accessing phastCons files.\
		Returns average conservation score for regions defined in BED file')
	parser.add_argument('bed', help='BED file of interval regions.')
	parser.add_argument('cons_folder', help='path to folder with phastCons files')
	parser.add_argument('output', help='name of output file')

	args = parser.parse_args()
	main(args.bed, args.cons_folder, args.output)
		



