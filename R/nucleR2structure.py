#!/usr/bin/env python3

import argparse
import itertools
import re
import sys


# simple useful functions #####################################################

def firsts(xs):
    """
    Given a list of tuples, return a lists of their first elements
    """
    return (x for x, *_ in xs)


def seconds(xs):
    """
    Given a list of tuples, return a lists of their second elements
    """
    return (x for _, x, *_ in xs)


def mapcat(f, *x):
    """
    Map a function and concatenate the results
    (a -> [b]) -> [a] -> [b]
    """
    return list(itertools.chain.from_iterable(map(f, *x)))


# functions to parse inputs ###################################################

def get_args(margin=4):
    """
    Parse the command line arguments
    """
    parser = argparse.ArgumentParser()
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument(
        "--calls",
        help="gff file with the nucleR calls",
        required=True
    )
    required_named.add_argument(
        "--genome_file",
        help="directory containing the reference genome in fasta format",
        required=True
    )
    required_named.add_argument(
        "--range",
        help="genomic range to process (ex: chrI:2000-2600)",
        required=True
    )
    required_named.add_argument(
        "--seq_output",
        help="output file for the sequence",
        required=True
    )
    required_named.add_argument(
        "--nucs_output",
        help="output file for the nucleosomes",
        required=True
    )
    parser.add_argument(
        "--margin",
        help="number of bases around selected nucleosomes",
        type=int,
        default=margin
    )
    args = parser.parse_args()
    return (args.calls,
            args.genome_file,
            args.range,
            args.margin,
            args.nucs_output,
            args.seq_output)


def parse_range(x, pattern=r"(?P<chr>.+):(?P<start>\d+)\.\.(?P<end>\d+)"):
    """
    Parse the raw range string (p.e.: chrI:300..4000) into a tuple containing
    the chromosome, the start and the end
    """
    match = re.compile(pattern).match(x)
    chr = match.group("chr")
    start = int(match.group("start"))
    end = int(match.group("end"))
    return chr, start, end


def parse_line(line):
    """
    Parse a row in the gff into a tuple containing the chromosome, the start,
    the end and the score
    """
    line_chr, _, _, line_start, line_end, line_score, *_ = line.split("\t")
    return line_chr, int(line_start), int(line_end), float(line_score)


def nuc_in_range(nuc, start, end):
    """
    Is a nucleosome inside a given range?
    """
    s, e, _ = nuc
    return s >= start and e <= end


def filter_lines(x, chr, start, end):
    """
    Return wether an entry in the gff is in the proper chromosome and range
    """
    c, *nuc = x
    return chr == c and nuc_in_range(nuc, start, end)


def get_nucs(fh, chr, start, end):
    """
    Given a filehandler for a gff file and a range, return the nucleosomes
    that are in that range (as tuples containing start, end and score)
    """
    for line in fh:
        parsed = parse_line(line)
        if filter_lines(parsed, chr, start, end):
            _, s, e, score = parsed
            yield s, e, score


def read_fasta(f):
    """
    Read a fasta file and return it as a dictionary with sequence ids as keys
    and the sequences as values
    """
    def iter(xs):
        for x in xs:
            lines = x.split()
            chr = lines[0]
            seq = "".join(lines[1:])
            yield chr, seq
    with open(f) as fh:
        txt = fh.read()
    res = {k: v for k, v in iter(txt.split(">")[1:])}
    return res


# functions to remove overlapping nucleosomes #################################

def nuc_grouper(nucs):
    """
    Given a list of nucleosomes, return a list of lists of nucleosomes,
    grouping the ones that overlap
    """
    mrans = merge_rans(nucs)
    res = [[x for x in nucs if nuc_in_range(x, start, end)]
            for start, end in mrans]
    return res


def merge_rans(ranges):
    """
    Given a lists of nucleosomes, return a list with their merged ranges
    """
    ranges = list(ranges)
    starts = list(firsts(ranges))
    ends = list(seconds(ranges))

    ii = [1] * len(starts)
    jj = [-1] * len(ends)

    stends = list(zip(starts+ends, ii+jj))
    stends.sort(key=lambda x: x[0])

    xs = list(itertools.accumulate(seconds(stends)))

    si = (i for i, _ in enumerate(xs) if xs[i] != 0 and xs[i-1] == 0)
    ei = (i for i, x in enumerate(xs) if x == 0)

    cstarts = (stends[i][0] for i in si)
    cends = (stends[i][0] for i in ei)

    return list(zip(cstarts, cends))


def rm_least(xs):
    """
    From a list of nucleosomes, remove the one with the lowest score
    """
    i, _ = min(enumerate(x for _, _, x in xs), key=lambda a: a[1])
    xs.pop(i)
    return xs


def unoverlapper(nucs):
    """
    Given a list of potentially overlapped nucleosomes, return a new list of
    nucleosomes with no overlaps. Discard overlapping nucleosomes with lowest
    scores.
    Group the list of nucleosomes into sublists of overlapping nucleosomes
    (a list with only one nucleosome, means it's not overlapped). For each
    subgroup, check if there's an overlap (more than one nucleosome in that
    sublist). If there's an overlap, remove the nucleosome with the lowest
    score, and recursively try again until no overlaps are found. In the end,
    concantenate all the non-overlapped nucleosomes into one single list.
    """
    def go(x):
        if len(x) == 1:  # no overlaps
            return x
        else:  # some overlaps; remove one, group the results and try again
            return mapcat(go, nuc_grouper(rm_least(x)))
    return mapcat(go, nuc_grouper(nucs))


# other helper functions ######################################################

def get_range(nucs):
    """
    Given a list of nucleosomes, return the whole range they occupy
    """
    start = min(firsts(nucs))
    end = max(seconds(nucs))
    return start, end


def shift(x, by=0):
    """
    Given a nucleosome (tuple with start, end and score), shift upstream it
    by a given amount
    """
    start, end, score = x
    return start-by, end-by, score


def expand_range(x, by=0):
    """
    Given a range, expand it (upstream and downstream) by a given amount
    """
    start, end = x
    return start-by, end+by


def cut_seq(seq, ranges):
    """
    Given a sequence and a list of ranges, return the sub-sequences for those
    ranges
    """
    ii = [0] + list(seconds(ranges))
    jj = list(firsts(ranges)) + [len(seq)]
    subs = list(seq[i:j] for i, j in zip(ii, jj))
    return subs


def get_cut_positions(nucs):
    """
    Get the relative positions of the start of the nucleosomes after cutting
    out the nucleosome range. The resulting positions should be spaced only by
    linker fragments
    """
    diff_sums = itertools.accumulate(y-x for x, y, _ in nucs)
    cuttings = [0] + list(diff_sums)[:-1]
    return [x-c for (x, _, _), c in zip(nucs, cuttings)]


###############################################################################


def main():

    # get user supplied inputs and arguments #############################

    calls_file, genome_file, range_str, margin, nucs_f, seq_f = get_args()
    chr, start, end = parse_range(range_str)

    # read and parse input files #########################################

    with open(calls_file) as fh:
        nucs = list(get_nucs(fh, chr, start, end))

    fasta = read_fasta(genome_file)
    seq = fasta[chr]

    # do all the processing ##############################################

    unoverlapped = unoverlapper(nucs)

    sel_start, sel_end = expand_range(get_range(nucs), margin)
    shifted_nucs = [shift(x, sel_start) for x in unoverlapped]
    sub_seq = seq[sel_start:sel_end]

    nucs_pos = get_cut_positions(shifted_nucs)
    sub_seqs = cut_seq(sub_seq, shifted_nucs)

    # save the output files ##############################################

    out_nucs = ' '.join(map(str, nucs_pos)) + "\n"
    out_seq = ''.join(sub_seqs) + "\n"

    with open(nucs_f, 'w') as fh:
        fh.write(out_nucs)
    with open(seq_f, 'w') as fh:
        fh.write(out_seq)

    return 0


###############################################################################

if __name__ == "__main__":
    sys.exit(main())

###############################################################################
