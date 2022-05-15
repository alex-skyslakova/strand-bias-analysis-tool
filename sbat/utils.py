import itertools
import os
import re
import sys

import dateutil.rrule as rrule
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

# dict of nucleotides
# duplicities for small letters are in order to better up performance
# (direct dict lookup is faster than converting char to uppercase and then lookup)
COMPLEMENTS_DICT = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A',
                    'a': 'T',
                    'c': 'G',
                    't': 'A',
                    'g': 'C'}

# basic ISO suffices
UNITS = {"K": 10 ** 3,
         "M": 10 ** 6,
         "G": 10 ** 9,
         "T": 10 ** 12}


def get_filename(filepath):
    """
    Function to convert path to file into a filename without filetype suffix.
    :param filepath: path to be transformed
    :return: name of the file
    """
    if filepath is None:
        return filepath
    return os.path.basename(filepath).split('.')[0]


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
    return round(100 - (ratio * 100), 4)


def parse_iso_size(size):
    """
    Function to split number in shortened format (i.e. 500K) into classic representation (500000
    :param size: number to be converted
    :return: long representation of number
    """
    if not size.isalnum():
        sys.exit("SIZE: expecting positive number, optionally in ISO format")
    if not size.isnumeric():
        _, number, unit = re.split(r'(\d+)', size)
        if unit in UNITS.keys():
            result = int(float(number) * UNITS[unit])
        else:
            sys.exit("unsupported suffix: {}, use one of K M G T".format(unit))
    else:
        result = size
    try:
        return int(result)
    except:
        sys.exit("--size parameter must be integer")


def gc_percentage(string):
    """
    Function to compute percentage of Cs and Gs in string.
    :param string: string to compute percentage in
    :return: percentage of GC content
    """
    gc = (string.count("C") + string.count("G")) / len(string) * 100
    return round(gc, 3)


def split_forwards_and_backwards(all_kmers):
    """
    Function to divide all possible kmers of length k into two groups, where one contains complements of another
    """

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
        reverse_complement += COMPLEMENTS_DICT[n]
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


def select_more_frequent(row, seq=False):
    """
    Fuction to return count (sequence if seq=True) of sequence or its rev. complement based on which is more frequent.
    :param row: one row of DataFrame with SB statistics
    :param seq: if true, sequence is returned instead of its count
    :return: count/sequence of more frequent out of seq and its rev. compl.
    """
    if row["seq_count"] > row["rev_complement_count"]:
        if seq:
            return row["seq"]
        return row["seq_count"]
    else:
        if seq:
            return row["rev_complement"]
        return row["rev_complement_count"]
