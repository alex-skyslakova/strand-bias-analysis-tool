import itertools
import os
import subprocess

import math
import time
from shutil import which
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import utils
from datetime import datetime

COMPLEMENTS_DICT = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'}

__version__ = '0.1.0'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        nargs='+',
                        type=str,
                        help='fasta or fastq file to count and analyze strand bias on')
    parser.add_argument('-v', '--version',
                        action="store_true",
                        default=False,
                        help='current version of the application'
                        )
    parser.add_argument('-n', '--no-jellyfish',
                        default=False,
                        help='skip k-mer counting. Requires input in fasta file where id=count, seq=k-mer')
    parser.add_argument('-o', '--output',
                        # action='output',
                        nargs=1,
                        default="sbat/",
                        help='output directory')
    parser.add_argument('-m', '--mer',
                        nargs=1,
                        default=0,
                        help='k-mer to count and analyze bias for. When set to 0, bias is analyzed for all k-mer '
                             'in 5-10 range. MER must be >= 3 for analysis')
    parser.add_argument('-s', '--size',
                        type=int,
                        nargs=1,
                        help='size of hash table for jellyfish')
    parser.add_argument('-t', '--threads',
                        nargs=1,
                        default=1,
                        help='number of threads jellyfish shall use for computations')
    parser.add_argument('-c', '--keep-computations',
                        action="store_true",
                        default=False,
                        help='keep jellyfish outputs and computations produced as partial results')
    args = parser.parse_args()

    if args.version:
        version()

    for input_file in args.input:
        if not os.path.isfile(input_file):
            print("no such file: " + input_file)
            parser.print_usage()

    if args.output is not None and not os.path.isdir(args.output[0]):
        print("no such directory: " + args.output[0])
        parser.print_usage()

    if not args.no_jellyfish:
        if args.threads < 1:
            print("number of threads must be a positive integer")
            return
        if args.mer < 3:
            print("MER must be a positive number higher or equal to 3")
        for file in args.input:
            print(file)
            # run_jellyfish(input_file=args.input, k=args.mer, t=args.threads)
    else:
        print("would run only analysis")


def version():
    """
    Prints current version of the tool.
    """
    print("StrandBiasAnalysisTool v" + __version__)


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


def run_jellyfish(input_file, output_dir, k=7, t=1, s='500M'):
    dump_file = os.path.join(output_dir, "mer_counts.jf")
    calculate = "jellyfish count -m " + str(k) + " -s " + s + " -t " + str(t) + " -o " + dump_file + " " + input_file
    print(calculate)
    subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
    dump = "jellyfish dump " + dump_file
    output_file = os.path.join(output_dir, "output_" + str(k) + "_" + os.path.basename(input_file))
    with open(output_file, "w") as outfile:
        subprocess.run(dump.split(" "), stdout=outfile)


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
            # print(count, sequence)
            if sequence == get_reverse_complement(sequence):
                continue
            c += 1
            seq.append(sequence)
            seq_count.append(float(count))
    return seq, seq_count


def jellyfish_to_dataframe(path, k):
    """
    Function to create dataframe with statistics based on Jellyfish output
    :param path: path to Jellyfish output
    :param k: length of kmers in given file
    """
    seq, seq_count = parse_fasta(path)

    # create dataframe with k-mers and their counts
    jellyfish_data = pd.DataFrame(
        data={'seq': seq,
              'seq_count': seq_count},
        index=seq,
        columns=['seq', 'seq_count'])

    # add column with reverse complements
    jellyfish_data['rev_complement'] = jellyfish_data.apply(lambda row: get_reverse_complement(row["seq"]),
                                                            axis=1)
    # split sequence set into forward and backward sequences (so that k-mer and its reverse complement
    # are not together in group)
    fwd_kmers, bwd_kmers = split_forwards_and_backwards(k)

    # remove backward group from DataFrame, as it is already represented as reverse complement to some
    # other k-mer in the DataFrame
    jellyfish_forward = jellyfish_data.drop(bwd_kmers, errors="ignore").set_index("rev_complement")

    # join forward DF with original one on index (equal to forward sequence) in order to connect info about
    # forward and backward datasets
    jellyfish_data = jellyfish_data.reset_index().join(jellyfish_forward, on="index", rsuffix="_", lsuffix="").drop(
        columns=["seq_", "index"], axis=1).dropna()
    jellyfish_data.rename(columns={"seq_count_": "rev_complement_count"}, inplace=True)

    # calculate ratio of forward and backward k-mer frequencies
    jellyfish_data["ratio"] = round(jellyfish_data["seq_count"] / jellyfish_data["rev_complement_count"], 5)
    # calculate deviation from 100% accuracy
    jellyfish_data["strand_bias_%"] = jellyfish_data.apply(lambda row: get_strand_bias_percentage(row["ratio"]), axis=1)
    # calculate CG content percentage
    jellyfish_data["CG_%"] = jellyfish_data.apply(lambda row: gc_percentage(row["seq"]), axis=1)
    # sort data by bias in descending order
    jellyfish_data = jellyfish_data.sort_values(by=["strand_bias_%"], ascending=False)
    # save DataFrame into csv
    jellyfish_data.to_csv("df_" + os.path.basename(path.split(".")[0]) + ".csv", index=False)


def get_strand_bias_percentage(ratio):
    """
    Function to count difference between kmers ratio and expected value (1, which means no bias) in percents.
    :param ratio: ratio of kmer and its reverse complement
    :return: deviation from 1 in %
    """
    dev = ratio - 1 if ratio >= 1 else 1 - ratio
    return round(dev * 100, 3)


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
            forwards.add(kmer)
            complements.add(get_reverse_complement(kmer))

    return list(forwards), list(complements)

def to_string(my_list):
    return ''.join(my_list)


def convert(seq, forward):
    if seq in forward:
        return get_reverse_complement(seq)
    return seq


for k in range(2, 11):
    print(datetime.now())
    print(k)
    run_jellyfish("pacbio_m54238_180628_014238.Q20.fastq", k)
    jellyfish_to_dataframe("output_" + str(k) + "_pacbio_m54238_180628_014238.Q20.fastq", k)
    print(datetime.now())
