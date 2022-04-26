import itertools
import os
import subprocess

import math
import sys
import time
from shutil import which
import argparse

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
import utils
from datetime import datetime
import analysis

COMPLEMENTS_DICT = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'}

__version__ = '0.1.0'

JELLYFISH_TO_DF_FORMAT = "output_{0}_{1}.fasta"
JELLYFISH_TO_DF_BATCHES = "output_{0}_nanopore_{1}_batch_{2}.fasta"

def arg_parser():
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
    parser.add_argument('-j', '--no-jellyfish',
                        action='store_true',
                        default=False,
                        help='skip k-mer counting. Requires input in fasta file where id=count, seq=k-mer')
    parser.add_argument('-o', '--output',
                        # action='output',
                        nargs=1,
                        default=["sbat/"],
                        help='output directory')
    parser.add_argument('-m', '--mer',
                        nargs=1,
                        default=[0],
                        type=int,
                        help='k-mer to count and analyze bias for. When set to 0, bias is analyzed for all k-mer '
                             'in 5-10 range. MER must be >= 3 for analysis')
    parser.add_argument('-s', '--size',
                        nargs=1,
                        help='size of hash table for jellyfish')
    parser.add_argument('-t', '--threads',
                        nargs=1,
                        type=int,
                        default=[1],
                        help='number of threads jellyfish shall use for computations')
    parser.add_argument('-c', '--keep-computations',
                        action="store_true",
                        default=False,
                        help='keep jellyfish outputs and computations produced as partial results')
    parser.add_argument('-n', '--detect-nanopore',
                        action="store_true",
                        default=False,
                        help='identify nanopore datasets among inputs and provide advanced analysis')
    args = parser.parse_args()
    return args, parser


def main():
    args, parser = arg_parser()

    if args.version:
        version()
        return 1

    for input_file in args.input:
        if not os.path.isfile(input_file):
            print("no such file: " + input_file)
            parser.print_usage()
            return 1

    if not os.path.isdir(args.output[0]):
        os.mkdir(args.output[0])

    if args.mer[0] < 3 and args.mer[0] != 0:
        print("MER must be a positive number higher or equal to 3")
        return 1

    if args.mer[0] == 0:
        start_k = 5
        end_k = 10
    else:  # run with specified k
        start_k = args.mer[0]
        end_k = args.mer[0]

    if args.threads[0] < 1:
        print("number of threads must be a positive integer")
        return 1

    jellyfish_outdir = os.path.join(args.output[0], 'jellyfish', '')
    if not os.path.isdir(jellyfish_outdir):
        os.mkdir(jellyfish_outdir)

    for file in args.input:
        print(args.output[0], str(time.time()).split('.')[0])
        sb_analysis = analysis.init_analysis(args.output[0], str(time.time()).split('.')[0])
        print("input: " + file)

        for k in range(start_k, end_k + 1):
            if not args.no_jellyfish:
                print("running computation and analysis for K=" + str(k))
                jf_output = run_jellyfish(input_file=file, output_dir=jellyfish_outdir, k=k, t=args.threads[0])

                if args.detect_nanopore and check_if_nanopore(file):
                    print('would analyze nanopore')
                    nanopore_analysis(file, args.output[0], start_k=start_k, end_k=end_k,
                                      threads=args.threads[0])  # TODO hash
            else:
                print("jellyfish disabled, running only analysis...")
                jf_output = file
            jellyfish_to_dataframe(jf_output, k, sb_analysis=sb_analysis, output="out/sbat/dump")


def version():
    """
    Prints current version of the tool.
    """
    print("StrandBiasAnalysisTool v" + __version__)


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


def run_jellyfish(input_file, output_dir, k=7, t=1, s='500M'):
    dump_file = os.path.join(output_dir, "mer_counts.jf")
    calculate = "jellyfish count -m " + str(k) + " -s " + s + " -t " + str(t) + " -o " + dump_file + " " + input_file
    print(calculate)
    subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
    dump = "jellyfish dump " + dump_file
    output_file = os.path.join(output_dir, "output_" + str(k) + "_" + os.path.basename(input_file))
    with open(output_file, "w") as outfile:
        subprocess.run(dump.split(" "), stdout=outfile)
    return output_file


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


def jellyfish_to_dataframe(path, k, output="sbat/dump", save=False, sb_analysis=None, batch=None):
    """
    Function to create dataframe with statistics based on Jellyfish output
    :param path: path to Jellyfish output
    :param k: length of kmers in given file
    """
    seq, seq_count = parse_fasta(path)
    if len(seq) == 0:
        sys.exit("no data parsed from {}".format(path))

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

    filename = utils.unique_path("df_{}.csv".format(os.path.basename(path.split(".")[0])))
    if sb_analysis:
        analysis.fill_sb_analysis_from_df("df_" + os.path.basename(path.split(".")[0]) + ".csv", jellyfish_data,
                                          k, batch, sb_analysis)
    if save:
        filename = utils.unique_path(os.path.join(output, filename))
        jellyfish_data.to_csv(filename, index=False)

    return jellyfish_data


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




# input = complete name and path
def nanopore_analysis(input, output_dir, start_k=5, end_k=10, threads=1, hash_size="300M", bin_interval=1):
    batch_files = analysis.bin_nanopore(input, bin_interval)
    #batch_files = ["nanopore/subsamples_50M/nanopore_nanopore_GM24385_11_batch_" + str(i) + ".fasta" for i in range(49)]
    if len(batch_files) < 2:
        print("data duration is too short for analysis of hour-long batches, aborting...")
        return
    sb_analysis = analysis.init_analysis(os.path.dirname(input), utils.get_filename(input))
    dataframe = pd.DataFrame(
        data={},
        columns=['seq', 'seq_count', 'rev_complement', 'rev_complement_count', 'ratio', 'strand_bias_%', 'CG_%'])
    for k in range(start_k, end_k + 1):
        batch_dfs = []
        for index, file in enumerate(batch_files):
            run_jellyfish(file, output_dir, k, threads, hash_size)
            df_file = os.path.join(output_dir, JELLYFISH_TO_DF_BATCHES.format(k, utils.get_filename(input), index))
            current_df = jellyfish_to_dataframe(df_file, k, sb_analysis=sb_analysis, batch=index)
            dataframe = pd.concat([dataframe, current_df])
            batch_dfs.append(current_df)
            #analysis.fill_sb_analysis_from_df(df_file, dataframe, k, index, sb_analysis)
        analysis.draw_conf_interval_graph(batch_dfs, output_dir, utils.get_filename(input), k)
        analysis.draw_basic_stats_lineplot(output_dir, utils.get_filename(input), sb_analysis, k, x_axis="batch")


if __name__ == '__main__':
    main()

# nanopore_analysis("nanotest.fasta", "", )
