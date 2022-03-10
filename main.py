import os
import subprocess

from Bio import SeqIO
import math

complements_dict = {'A': 'T',
                    'C': 'G',
                    'G': 'C',
                    'T': 'A'}


class KMerStatistics:
    def __init__(self, k):
        self.k = k
        self.kmers = []
        self.input_file = ""

    def add_kmer(self, mer, forward_count, backward_count):
        self.kmers.append(KMer(mer, forward_count, backward_count))

    def get_sorted_kmers(self):
        self.kmers.sort(key=lambda x: x.ratio, reverse=True)

    '''def print_by_ratio(self):
        self.get_sorted_kmers()
        for i in self.kmers:
            print("***")
            print("kmer: " + i.forward + " (" + str(i.forward_count) + ")")
            print("kmer: " + i.backward + " (" + str(i.backward_count) + ")")
            print("ratio: " + str(i.ratio))

    def print_by_ratio_csv(self):
        if self.k <= 7:
            self.get_sorted_kmers()
            kmers_to_csv(self.k, self.kmers, self.input_file)

        else:
            percentage = round(len(self.kmers) * 0.02)
            result = get_n_most_significant(self.kmers, percentage)
            kmers_to_csv(self.k, result, self.input_file)'''


class KMer:
    def __init__(self, mer, forward_count, backward_count):
        self.forward = mer
        self.forward_count = forward_count
        self.backward = get_reverse_complement(mer)
        self.backward_count = backward_count
        self.ratio = round(forward_count / backward_count if backward_count != 0 else math.nan, 3)


def read_fasta(path, output_type='list'):
    seq_list = [rec.seq for rec in SeqIO.parse(path, "fasta")]
    seq_dict = {rec.seq: int(rec.id) for rec in SeqIO.parse(path, "fasta")}  # dictionary of pairs 'sequence: count'
    if output_type == "list":
        return seq_list
    if output_type == "dict":
        return seq_dict


def get_reverse_complement(seq):
    """
    Function to produce reverse complement to given sequence.
    :param seq: string sequence from which reverse complement is wanted
    :return: string containing reverse complement
    """
    reverse_complement = ""
    for n in reversed(seq):
        reverse_complement += complements_dict[n.capitalize()]
    return reverse_complement


def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    from shutil import which

    return which(name) is not None


'''def get_kmer_ratios(path, k):
    jellyfish_result_dict = read_fasta(path, output_type="dict")
    stats = KMerStatistics(k)
    kmers = list(jellyfish_result_dict.keys())
    while kmers:

        kmer = kmers.pop()
        rev_complement = get_reverse_complement(kmer)
        if kmer == rev_complement:
            stats.add_kmer(kmer, jellyfish_result_dict[kmer] // 2, jellyfish_result_dict[kmer] // 2)
        else:
            stats.add_kmer(kmer, jellyfish_result_dict[kmer], jellyfish_result_dict[rev_complement])
            kmers.remove(rev_complement)
    return stats'''


def run_jellyfish(input_file):
    k = 11
    calculate = "jellyfish count -m " + str(k) + " -s 500M " + input_file
    subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
    dump = "jellyfish dump mer_counts.jf"
    with open("output_" + str(k) + "_" + os.path.basename(input_file), "w") as outfile:
        subprocess.run(dump.split(" "), stdout=outfile)
    #stats = get_kmer_ratios("output_" + str(k) + "_" + os.path.basename(input_file), k)
    #stats.input_file = os.path.basename(input_file)


def kmers_to_csv(k, kmers, output_file_name):
    with open("_ordered_output_" + str(k) + "_" + os.path.basename(output_file_name.split(".")[0] + ".csv"),
              'w') as output:
        output.write("forward_kmer,forward_count,backward_kmer,backward_count,ratio\n")
        for i in kmers:
            string = i.forward + "," + str(i.forward_count) + "," + i.backward + "," + str(
                i.backward_count) + "," + str(i.ratio) + "\n"
            output.write(str(string))


'''
def get_n_most_significant(kmers, n):
    best_n = []
    worst_n = []
    used_indexes_best = []
    used_indexes_worst = []
    for i in range(n):
        current_worst = None
        current_best = None
        for j in range(len(kmers)):
            if j in used_indexes_best or j in used_indexes_worst:
                continue
            kmer = kmers[j]
            if  current_worst is None or kmer.ratio > current_worst.ratio:
                current_worst = kmer
                used_indexes_worst.append(j)
                continue
            if current_best is None or kmer.ratio < current_best.ratio:
                current_best = kmer
                used_indexes_best.append(j)
                continue
        best_n.append(current_best)
        worst_n.append(current_worst)
    best_n.reverse()
    return worst_n + best_n

input_file = "GRCh38_latest_genomic.fna"
k = 9

stats = get_kmer_ratios("output_" + str(k) + "_" + os.path.basename(input_file), k)
print("Here")
stats.input_file = input_file
#stats.print_by_ratio()
stats.print_by_ratio_csv()'''


# //////////////////////////////////////////////////


def run_jellyfish2(input_file, k=7):
    calculate = "jellyfish count -m " + str(k) + " -s 500M " + input_file
    subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
    dump = "jellyfish dump mer_counts.jf"
    with open("output_" + str(k) + "_" + os.path.basename(input_file), "w") as outfile:
        subprocess.run(dump.split(" "), stdout=outfile)
    #stats = get_kmer_ratios("output_" + str(k) + "_" + os.path.basename(input_file), k)
    #stats.input_file = os.path.basename(input_file)


from Bio.SeqIO.FastaIO import SimpleFastaParser
import pandas as pd
from datetime import datetime


def jellyfish_to_DataFrame(path, k):
    seq = []
    seq_count = []
    #print(path)
    with open(path) as fasta_file:
        c= 0
        for count, sequence in SimpleFastaParser(fasta_file):
            #print(count, sequence)
            c += 1
            seq.append(sequence)
            seq_count.append(float(count))
    #print(c)
    jellyfish_data = pd.DataFrame(
        data={'seq': seq,
              'seq_count': seq_count},
        index=seq,
        columns=['seq', 'seq_count'])

    jellyfish_data['rev_complement'] = jellyfish_data.apply(lambda row: get_reverse_complement(row["seq"]),
                                                            axis=1)  # add column with
    fwd_kmers, bwd_kmers = get_kmers_to_remove(k)
    jellyfish_forward = jellyfish_data.drop(bwd_kmers, errors="ignore").set_index("rev_complement")
    jellyfish_data = jellyfish_data.reset_index().join(jellyfish_forward, on="index", rsuffix="_", lsuffix="").drop(
        columns=["seq_", "index"], axis=1).dropna()
    jellyfish_data.rename(columns={"seq_count_": "rev_complement_count"}, inplace=True)

    print(jellyfish_data.head())


    # jellyfish_data['rc_count'] = jellyfish_data.apply(lambda row: get_reverse_complement_count(row, jellyfish_data), axis=1)
    print(datetime.now())
    # for i in range(len(jellyfish_data)//2):
    #    rc = jellyfish_data.iloc[[i]]["seq"][0]
    #    jellyfish_data.drop(get_reverse_complement(rc), inplace=True)

    jellyfish_data["ratio"] = round(jellyfish_data["seq_count"] / jellyfish_data["rev_complement_count"], 5)
    #print(datetime.now())
    #print(jellyfish_data.head())

    jellyfish_data["deviation_%"] = jellyfish_data.apply(lambda row: get_deviation_percentage(row["ratio"]), axis=1)
    jellyfish_data["CG_%"] = jellyfish_data.apply(lambda row: gc_percentage(row["seq"]), axis=1)
    #print(datetime.now())
    jellyfish_data = jellyfish_data.sort_values(by=["deviation_%"], ascending=False)
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


def get_deviation_percentage(ratio):
    dev = ratio - 1 if ratio >= 1 else 1 - ratio
    return round(dev * 100, 3)


def gc_percentage(string):
    cg = (string.count("C") + string.count("G")) / len(string) * 100
    return round(cg, 3)


def kmers(n):
    return [p for p in itertools.product("ACGT", repeat=n)]


def repetitions(n):
    return list(itertools.combinations_with_replacement("ACGT", n))


def get_kmers_to_remove(k):
    all_kmers = [''.join(i) for i in  list(itertools.product("ACGT", repeat=k))]
    forwards = all_kmers[:len(all_kmers)//2]
    complements = []

    for i in forwards:#range(len(a) - 1, (len(a))//2 - 1, -1):
        complements.append(get_reverse_complement(i))
    return forwards, complements

def to_string(my_list):
    return ''.join(my_list)



import itertools
x = "ACGT"

def convert(seq, forward):
    if seq in forward:
        return get_reverse_complement(seq)
    return seq

for k in range(2, 11):

    print(datetime.now())
    print(k)
    run_jellyfish2("pacbio_m54238_180628_014238.Q20.fastq", k)
    jellyfish_to_DataFrame("output_" + str(k) + "_pacbio_m54238_180628_014238.Q20.fastq", k)
    print(datetime.now())

