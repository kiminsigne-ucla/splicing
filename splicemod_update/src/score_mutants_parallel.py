"""
This script takes as input an "imperfect" reference fasta file. It contains one entry for each mutant, with the header
describing the closest perfect reference sequence and alignment details. It will score each imperfect sequence similarly
to score_exons.py, scoring the splice donor/acceptor and splicing motifs.
"""

import argparse
from itertools import islice
from functools import partial
import multiprocessing
import Bio
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

import score
import feature


def parse_ref(filename):
    lib = {}
    # get intron/exon/intron boundaries of sequences needed for sequence features of imperfect sequences
    with open(filename) as infile:
        while True:
            line = list(islice(infile, 2))
            if not line:
                break
            header, seq = line
            # replace '= ' with '=' to make parsing easier
            header = header.replace('= ', '=')
            name = header.split(' ')[0][1:]
            fields = header.strip().split(' ')[2:]
            fields = {x.split('=')[0]: x.split('=')[1] for x in fields if '=' in x}
            if fields['strand'] != '':
                if fields['strand'] == '-':
                    fields['strand'] = -1
                elif fields['strand'] == '+':
                    fields['strand'] = 1
                if fields['strand'] == '.':
                    fields['strand'] = None
                else:
                    fields['strand'] = int(fields['strand'])
            if fields['len'] != '':
                fields['len'] = map(int, fields['len'].split('.'))
                if fields['strand'] == -1:
                    # switch order of introns so that upstream intron size is always first
                    fields['len'] = [fields['len'][2], fields['len'][1], fields['len'][0]]
            fields['seq'] = seq.strip().upper()
            lib[name] = fields
    return lib


def seq_iter(filename):
    """
    Iterate through a tab-delimited sorted sequence file of imperfect sequences, first column header, second column
    sequence. This function returns all mutants of the same sequence in chunks, so the intron/exon features only need to
    be created once and the original sequence scored once.

    :param filename: Name of sorted tab-delimited sequence file
    :return: a dictionary containing all mutants of the same perfect sequence
    """

    with open(filename) as infile:
        # initialize
        mutants = []
        header, seq = next(infile).strip().split('\t')
        name = header.split('.')[0] # remove trailing mutant number e.g. ENSE00000332835_000.0010
        mutants.append((header, seq))
        curr_name = name
        for line in infile:
            header, seq = line.strip().split('\t')
            name = header.split('.')[0]
            if name != curr_name:
                yield mutants
                mutants = []
                curr_name = name
            else:
                mutants.append((header, seq))


def score_mutants(mut_list, ref, min_length):
    """
    This function takes a dictionary of mutants (all of the same reference sequence) and
    - filters out those that are too short
    - creates intron/exon features for all mutants based on the reference sequence
    - scores each mutant and reference sequence

    :param mut_list: list of tuples of mutants of the same reference sequence
    :param ref: dictionary of reference sequences
    :param min_length: Minimum length of mutant
    :return: tab-delimited result string, containing mutant header, sequence, mutant score, and original sequence score
    """

    scores = []

    # convert mutant list to dictionary
    mutants = {x[0]: x[1] for x in mut_list}
    if len(mutants) == 0:
        return ['']
    # get corresponding reference sequence
    name = mutants.keys()[0].split()[0].split('.')[0]
    ref_features = ref.get(name, None)
    if not ref_features:
        for header, seq in mutants.items():
            scores.append('\t'.join([header, seq, 'None', 'None']) + '\n')
        return scores

    # create SeqRecord for original feature
    orig_record = SeqRecord(Seq(ref_features['seq']), id=name)
    orig_record.annotations = ref[name]

    # initialize with intron/exon features
    # get corresponding boundaries from reference
    upstr_intron_size, exon_size, downstr_intron_size = ref_features['len']

    # calculate the feature start/ends
    upstr_loc = (0, upstr_intron_size)
    exon_loc = (upstr_loc[1], upstr_loc[1] + exon_size)
    downstr_loc = (exon_loc[1], exon_loc[1] + downstr_intron_size)

    # create feature for mutant and for original sequence
    upstr_intron_feature = SeqFeature(
        FeatureLocation(*upstr_loc),
        type="intron",
        strand=orig_record.annotations["strand"])

    exon_feature = SeqFeature(
        FeatureLocation(*exon_loc),
        type="exon",
        strand=orig_record.annotations["strand"])

    downstr_intron_feature = SeqFeature(
        FeatureLocation(*downstr_loc),
        type="intron",
        strand=orig_record.annotations["strand"])

    orig_record.features.extend([upstr_intron_feature, exon_feature, downstr_intron_feature])

    orig_record, orig_score = score_sequence(orig_record)

    for header, seq in mutants.items():
        if len(seq) < min_length:
            continue

        # create SeqRecord object
        record = SeqRecord(Seq(seq.strip()), id=header.split()[0][1:])
        annot_list = ['unknown', 'cigar', 'md', 'alignment']
        for k, v in zip(annot_list, header.strip().split()[1:]):
            record.annotations[k] = v

        # initialize strand
        record.annotations['strand'] = ref_features['strand']

        record.features.extend([upstr_intron_feature, exon_feature, downstr_intron_feature])

        record, mut_score = score_sequence(record)

        scores.append('\t'.join([header, seq, str(mut_score), str(orig_score)]) + '\n')

    return scores


def score_sequence(record):
    """
    :param SeqRecord: SeqRecord object, containing features to score
    :return: sequence score
    """

    # find motifs
    record.populate_attribs()
    # make sure correct motifs are present
    record.get_correct_motifs()

    seq_score = score.score_seq_record(record)

    return record, seq_score


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='Name of input tab-separated .txt file of imperfect sequences')
    parser.add_argument('output', help='Name of output file')
    parser.add_argument('ref', help='Reference fasta file for library')
    parser.add_argument('min_length', type=int, help='Minimum length of mutant')
    parser.add_argument('num_processes', type=int, help='Number of parallel processes')

    args = parser.parse_args()

    ref = parse_ref(args.ref)

    mut_chunks = seq_iter(args.input)

    # assign constant parameter to function
    func = partial(score_mutants, ref=ref, min_length=args.min_length)
    pool = multiprocessing.Pool(processes=args.num_processes)
    mut_scores = pool.map(func, mut_chunks)

    with open(args.output, 'w') as outfile:
        for scores in mut_scores:
            for x in scores:
                outfile.write(x)