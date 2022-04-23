import itertools
import os
import subprocess

import math
from shutil import which

from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from datetime import datetime

COMPLEMENTS_DICT = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'}


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


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None


def run_jellyfish(input_file, k=7):
    calculate = "jellyfish count -m " + str(k) + " -s 500M " + input_file
    subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
    dump = "jellyfish dump mer_counts.jf"
    with open("output_" + str(k) + "_" + os.path.basename(input_file), "w") as outfile:
        subprocess.run(dump.split(" "), stdout=outfile)


def parse_fasta(path):
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


def get_reverse_complement_count(act_row, data):
    try:
        count = data.loc[act_row['rev_complement']]['seq_count']
        return float(count)
    except KeyError:
        return math.nan


def to_int(data):
    data = data.apply(lambda x: x.astype(int) if x is not None else math.nan)
    return data


def get_strand_bias_percentage(ratio):
    dev = ratio - 1 if ratio >= 1 else 1 - ratio
    return round(dev * 100, 3)


def gc_percentage(string):
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
