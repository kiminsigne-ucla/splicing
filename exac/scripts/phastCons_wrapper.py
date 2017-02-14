"""
This script returns the average conservation (phastCons) score for a 
given interval region, as defined in the sorted input bed file.

"""

import argparse
import pandas as pd
import numpy as np
import gzip
import multiprocessing
from functools import partial


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


# def get_interval(itr, interval):
#     e = next(itr)
#     while e[0] < interval[0]:
#         e = next(itr)
#     while in_interval(e[0], interval):
#         yield e
#         e = next(itr)
#     # note: skips first base outside of interval, but ours are
#     # not overlapping


def get_interval(itr, interval, curr_base=None):
	if not curr_base:
		# initialize
		e = next(itr)
	else:
		# start at current base, do not iterate through an extra base
		e = curr_base
	while e[0] < interval[0]:
		e = next(itr)
	while in_interval(e[0], interval):
		yield e 
		e = next(itr)
	yield e      
    # note: will always yield one base outside of region so that no base
    # is ever skipped


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


def chrom_extract_cons(chrom_region, cons_folder):
	"""
	This function takes as input a chromosome (e.g. chr9) and regions in that
	chromosome as parsed by bed_iter(). It opens the appropriate wig file and
	extracts the mean conservation scores for each region.

	chrom_region: chromosome regions, list of lists. Each element in the list is a
	region, which is itself a list with four fields (chromosome, start, end, name).
	This is the format returned from bed_iter()
	cons_folder: directory where conservation files are located

	return: chrom_region with an extra field for conservation added to each region
	"""
	# suffix common to all files
	cons_suffix = '.phastCons100way.wigFix.gz'
	# get chromosome from chrom_region
	chrom = chrom_region[0][0]
	# open appropriate wig file
	wig = WigReader(cons_folder + chrom + cons_suffix)

	# initialize iter with first region
	# region format: [chrN, start, end, name]
	# outputs tuple, (position, score)
	first_region = chrom_region.pop(0)
	cons_scores = list(get_interval(iter(wig), (int(first_region[1]), int(first_region[2]))))
	# pop first base outside interval from end of list
	curr_base = cons_scores.pop()
	# get average
	cons_mean = np.mean([x[1] for x in cons_scores])
	first_region.append(str(cons_mean))


	for region in chrom_region:
		# initialize with current base so that no base is ever skipped
		cons_scores = list(get_interval(iter(wig), (int(region[1]), int(region[2])),
			curr_base=curr_base))
		if len(cons_scores) >= 1:
			# pop first base outside interval from end of list
			curr_base = cons_scores.pop()
		# get average
		cons_mean = np.mean([x[1] for x in cons_scores])
		region.append(str(cons_mean))
	
	# finally, add the first region we popped in the beginning
	chrom_region.insert(1, first_region)
	return chrom_region


def main(bed, cons_folder, output, num_processes):

	outfile = open(output, 'w')

	bed_chroms = bed_iter(bed)
	
	pool = multiprocessing.Pool(processes=num_processes)

	# assign constant parameter to function
	func = partial(chrom_extract_cons, cons_folder=cons_folder)
	# map each chromosome in bed_chroms to the function and run in parallel
	cons_scores = pool.map(func, bed_chroms)
    # write to file
	with open(output, 'w') as outfile:
		for chrom_region in cons_scores:
			for region in chrom_region:
				outfile.write('\t'.join(region)+'\n')


if __name__ == '__main__':
	parser = argparse.ArgumentParser('wrapper for accessing phastCons files.\
		Returns average conservation score for regions defined in BED file')
	parser.add_argument('bed', help='BED file of interval regions.')
	parser.add_argument('cons_folder', help='path to folder with compressed phastCons files')
	parser.add_argument('output', help='name of output file')
	parser.add_argument('num_processes', type=int, default=3,
		help='number of parallel processes')

	args = parser.parse_args()
	main(args.bed, args.cons_folder, args.output, args.num_processes)
		



