import itertools
import os
import re

import dateutil.rrule as rrule
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

COMPLEMENTS_DICT = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'}


def get_filename(name):
    if name is None:
        return name
    return os.path.basename(name).split('.')[0]


# TODO not needed?
# def to_int(data):
#     data = data.apply(lambda x: x.astype(int) if x is not None else math.nan)
#     return data


def unique_path(f):
    """
    Function to check if given filepath exists, and if yes, append it with next non-taken number.
    """
    fnew = f
    root, ext = os.path.splitext(f)
    i = 0
    while os.path.exists(fnew):
        i += 1
        fnew = '%s_%i%s' % (root, i, ext)
    return fnew


def is_or_create_dir(dir):
    """
    Function to create directory if it does not exists.
    :param dir: directory to be created
    """
    if not os.path.isdir(dir):
        os.makedirs(dir)


def get_bin_number(text):
    """
    Get appendix from file number representing nanopore bin.
    """
    return list(map(int, re.findall(r'\d+', text)))[-1]


def get_n_percent(df, n, tail=False):
    """
    Function to extract top (tail=False) or bottom (tail=True) n percent of data sorted by strand bias.
    :param df: dataframe from which we want to extract part of data
    :param n: number of percent of data we want to extract
    :param tail: if false, data are taken from top of the dataset, else from tail
    :return: DataFrame containing top or bottom n% of input DF

    """
    if tail:
        return df.tail(int(len(df) * (n / 100)))
    else:
        return df.head(int(len(df) * (n / 100)))


def hours_aligned(start, end, interval=1):
    """
    Function that finds closest whole hour and splits the duration into 'interval' hours long chunk until it
    reaches end timestamp.
    :param start: oldest timestamp of dataset
    :param end: newest timestamp of dataset
    :param interval: desired length of one chunk
    :return: list of timestamps
    """
    chunks = []
    rule = rrule.rrule(rrule.HOURLY, byminute=0, bysecond=0, dtstart=start, interval=interval)
    for x in rule.between(start, end, inc=False):
        chunks.append(x)
    chunks.append(end)
    return chunks


def get_ratio(x, y):
    """
    Function to get ratio of two numbers.
    :param x: x
    :param y: y
    """
    if x < y:
        return round(x / y, 6)
    return round(y / x, 6)


def get_strand_bias_percentage(ratio):
    """
    Function to count Strand Bias percentage from ratio of two numbers.
    :param ratio: ratio of kmer and its reverse complement
    :return: deviation from 1 in %
    """
    return 100 - (ratio * 100)


def gc_percentage(string):
    """
    Function to compute percentage of Cs and Gs in string.
    :param string: string to compute percentage in
    :return: percentage of CG content
    """
    cg = (string.count("C") + string.count("G")) / len(string) * 100
    return round(cg, 3)


def split_forwards_and_backwards(k):
    """
    Function to divide all possible kmers of length k into two groups, where one contains complements of another
    """
    # generate all possible combinations of length k
    all_kmers = [''.join(i) for i in list(itertools.product("ACGT", repeat=k))]

    # use sets for faster lookup
    forwards = set()
    complements = set()

    for kmer in all_kmers:
        if kmer in complements:
            continue
        else:
            kmers = [kmer, get_reverse_complement(kmer)]
            kmers.sort()
            forwards.add(kmers[0])
            complements.add(kmers[1])

    return list(forwards), list(complements)


def check_if_nanopore(path):
    """
    Function to check if file contains timestamp in description (specific for nanopore data)
    :param path: path to file
    :return: bool
    """
    for rec in SeqIO.parse(path, path.split('.')[-1]):
        if "start_time=" in rec.description:
            return True
        return False


def get_reverse_complement(seq):
    """
    Function to produce reverse complement to given sequence.
    :param seq: string sequence from which reverse complement is wanted
    :return: string containing reverse complement
    """
    reverse_complement = ""
    for n in reversed(seq):
        reverse_complement += COMPLEMENTS_DICT[n.capitalize()]
    return reverse_complement


def parse_fasta(path):
    """
    Function to parse Jellyfish fasta output into list of sequences and list of their counts
    :param path: path to Jellyfish output
    :return: list of sequences, list of their counts
    """
    seq = []
    seq_count = []
    with open(path) as fasta_file:
        c = 0
        for count, sequence in SimpleFastaParser(fasta_file):
            if sequence == get_reverse_complement(sequence):
                continue
            c += 1
            seq.append(sequence)
            seq_count.append(float(count))
    return seq, seq_count
