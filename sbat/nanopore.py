import datetime
import math
import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from pytz import utc
from dateutil.parser import parse as dparse

from sbat.utils import get_filename, hours_aligned, unique_path
import bisect

NANOPORE_BIN_FORMAT = 'nanopore_{0}_bin_{1}.fasta'
NANOPORE_BIN_DF_FORMAT = 'output_{0}_nanopore_{1}_bin_{2}.csv'


class Nanopore:
    """
    Class taking care of Nanopore analysis.
    """
    def __init__(self, interval=1):
        self.common = None  # instance of primary Analysis
        self.jf = None  # instance of Jellyfish
        self.subs_reads = math.inf  # max number of reads in bin
        self.subs_bases = math.inf  # max number of nucleotides in bin
        self.bin_interval = interval  # time interval of one bin

    def init_common(self, analysis, jf):
        """
        Method for assigning instances of Jellyfish and Analysis.
        """
        self.common = analysis
        self.jf = jf

    def nanopore_analysis(self, input):
        """
        Method for analysing split bins, calls method for splitting dataset and then iterates over newly
        created datasets and performs analysis.
        """
        bin_files = self.bin_nanopore(input, self.bin_interval)
        if len(bin_files) < 2:
            print("data duration is too short for analysis of {}-hour-long bins, aborting...".format(self.bin_interval))
            return
        self.common.np_sb_analysis_file = self.common.init_analysis(nanopore=True)
        dataframe = pd.DataFrame(
            data={},
            columns=['seq', 'seq_count', 'rev_complement', 'rev_complement_count', 'ratio', 'strand_bias_%', 'CG_%'])

        # nanopore analysis runs only for the lowest k in the range due to more accurate results from small ks
        # it can be easily upgraded to run for all ks in range, but at the expense of speed
        bin_dfs = []
        for index, file in enumerate(bin_files):
            if file is None:
                bin_dfs.append(file)
                continue
            df_file = self.jf.run_jellyfish(file, self.common.start_k)
            current_df = self.common.jellyfish_to_dataframe(df_file, self.common.start_k, bin=index)
            dataframe = pd.concat([dataframe, current_df])
            bin_dfs.append(current_df)
        self.common.plot_conf_interval_graph(bin_dfs, self.common.start_k, start_index=self.common.start_k, nanopore=True)
        self.common.plot_basic_stats_lineplot(
            get_filename(input), self.common.np_sb_analysis_file, self.common.start_k, nanopore=True)
        self.common.track_most_common_kmer_change_freq(bin_dfs, self.common.start_k)

    def bin_nanopore(self, fastq, interval=1):
        """
        Method for splitting Nanopore dataset into smaller datasets based on time of sequencing. Each bin will
        contain data of 'interval' hours. Used for Nanopore bias comparation per hours.
        :param fastq: input dataset path
        :param interval: number of hours that should fall into one bin
        :return: list of newly created filenames with None on indexes, where bin was empty
        """
        print("splitting dataset...")
        subsampling = self.subs_reads != math.inf or self.subs_bases != math.inf
        file_type = 'fastq' if fastq.split('.')[-1] == 'fastq' else 'fasta'

        start = utc.localize(datetime.datetime.now())
        end = utc.localize(datetime.datetime(1970, 1, 1, 0, 0, 0))

        # go through dataset and find oldest and newest timestamp of sequencing
        for record in SeqIO.parse(fastq, file_type):
            record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
            if record_time < start:
                start = record_time
            if record_time > end:
                end = record_time
                
        # find closest whole hour to start and yield timestamps for each 'interval' hours until reaching end
        bins = hours_aligned(start, end, interval)
        reads_per_bin = [0 for _ in range(len(bins))]  # keep info about reads
        bases_per_bin = [0 for _ in range(len(bins))]  # keep info about bases
        bin_files = ["" for _ in range(len(bins))]  # filenames of newly created bins

        for record in SeqIO.parse(fastq, file_type):
            record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
            # find position of current record_time in list of timestamps -- equal to bin positioning
            bin = bisect.bisect_left(bins, record_time)
            if subsampling and (
                    bases_per_bin[bin] >= self.subs_bases or reads_per_bin[bin] >= self.subs_reads):
                continue
            reads_per_bin[bin] += 1
            bases_per_bin[bin] += len(record.seq)

            filename = os.path.join(self.common.dump_dir, NANOPORE_BIN_FORMAT.format(get_filename(fastq), bin))
            bin_files[bin] = filename
            f = open(filename, 'a')
            f.write(record.format('fasta'))
            f.close()

        self.plot_bin_distribution(reads_per_bin, get_filename(fastq), "Reads")
        self.plot_bin_distribution(bases_per_bin, get_filename(fastq), "Nucleotides")

        # return list of filenames for non-empty bins and None for the empty ones
        return [str(i) or None for i in bin_files]

    def plot_bin_distribution(self, counts_per_bin, filename, what_of="Reads"):
        """
        Method for plotting distribution of nucleotides or reads in bins.
        :param counts_per_bin: list of counts of nucleotides or reads in given file, where index of list = bin number
        :param filename: name of input dataset, used for determining name of output plot file
        :param what_of: "Reads" or "Nucleotides" based on what is to be analysed (just for plot name)
        :return: None
        """
        if counts_per_bin is None or counts_per_bin == []:
            return
        bins = [x for x in range(len(counts_per_bin))]

        fig, ax = plt.subplots(figsize=(18, 12))

        ax.set_ylabel('{} counts'.format(what_of), size=24)
        ax.set_title('{} Counts per Bin'.format(what_of), size=30)
        ax.set_xlabel('Bins', size=24)
        x = np.arange(len(bins)) * 2
        ax.set_xticks(x)
        plt.setp(ax.get_yticklabels(), fontsize=20)
        ax.set_xticklabels(bins, size=20)
        freq = (len(bins)//20) + 1
        for i, label in enumerate(ax.xaxis.get_ticklabels()):
            if i % freq == 0:
                continue
            label.set_visible(False)
        bar_width = 1
        cols = ax.bar(x - bar_width / 2, counts_per_bin, bar_width, label='{} counts per time slot (bin)'.format(what_of))
        plt.legend(fontsize=18)
        for col in cols:
            height = col.get_height()
            ax.annotate('{}'.format(height),
                        xy=(col.get_x() + col.get_width() / 2, height),
                        xytext=(0, 4),
                        textcoords="offset points",
                        ha='center', va='bottom', rotation=90, fontsize=16, clip_on=True)

        fig_name = unique_path(os.path.join(self.common.fig_dir,
                                            "fig_{}_per_bins_{}.png".format(what_of.lower(), filename)))
        plt.savefig(fig_name)
