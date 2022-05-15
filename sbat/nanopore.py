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

NANOPORE_BIN_FORMAT = 'nanopore_{0}_batch_{1}.fasta'
NANOPORE_BIN_DF_FORMAT = 'output_{0}_nanopore_{1}_batch_{2}.csv'


class Nanopore:
    def __init__(self, interval=1):
        self.common = None
        self.jf = None
        self.subs_reads = math.inf
        self.subs_bases = math.inf
        self.bin_interval = interval

    def init_common(self, analysis, jf):
        self.common = analysis
        self.jf = jf

    # input = complete name and path
    def nanopore_analysis(self, input):
        batch_files = self.bin_nanopore(input, self.bin_interval)
        #batch_files_ = ["nano_full/sbat_out/dump/nanopore_{}_batch_".format(get_filename(input)) + str(i) + ".fasta" for i in range(71)]
        #batch_files = []
        #for i in batch_files_:
        #    if os.path.isfile(i):
        #        batch_files.append(i)
        #print(batch_files)
        if len(batch_files) < 2:
            print("data duration is too short for analysis of hour-long batches, aborting...")
            return
        self.common.np_sb_analysis_file = self.common.init_analysis(nanopore=True)
        dataframe = pd.DataFrame(
            data={},
            columns=['seq', 'seq_count', 'rev_complement', 'rev_complement_count', 'ratio', 'strand_bias_%', 'CG_%'])

        for k in range(self.common.start_k, self.common.end_k + 1):
            batch_dfs = []
            for index, file in enumerate(batch_files):
                if file is None:
                    batch_dfs.append(file)
                    continue
                df_file = self.jf.run_jellyfish(file, k)
                current_df = self.common.jellyfish_to_dataframe(df_file, k, batch=index)
                dataframe = pd.concat([dataframe, current_df])
                batch_dfs.append(current_df)
            self.common.plot_conf_interval_graph(batch_dfs, k, start_index=self.common.start_k, nanopore=True)
            self.common.draw_basic_stats_lineplot(
                get_filename(input), self.common.np_sb_analysis_file, k, nanopore=True)
            self.common.track_most_common_kmer_change_freq(batch_dfs, k)

    # Nanopore bias comparation per hours
    def bin_nanopore(self, fastq, interval=1):
        print("in bin")
        subsampling = self.subs_reads != math.inf or self.subs_bases != math.inf
        file_type = 'fastq' if fastq.split('.')[-1] == 'fastq' else 'fasta'

        start = utc.localize(datetime.datetime.now())
        end = utc.localize(datetime.datetime(1970, 1, 1, 0, 0, 0))

        for record in SeqIO.parse(fastq, file_type):
            record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
            if record_time < start:
                start = record_time
            if record_time > end:
                end = record_time
        batches = hours_aligned(start, end, interval)
        reads_per_batch = [0 for _ in range(len(batches))]
        bases_per_batch = [0 for _ in range(len(batches))]
        batch_files = ["" for _ in range(len(batches))]

        for record in SeqIO.parse(fastq, file_type):
            print("here")
            record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
            batch = bisect.bisect_left(batches, record_time)
            if subsampling and (
                    bases_per_batch[batch] >= self.subs_bases or reads_per_batch[batch] >= self.subs_reads):
                continue
            reads_per_batch[batch] += 1
            bases_per_batch[batch] += len(record.seq)
            filename = os.path.join(self.common.dump_dir, NANOPORE_BIN_FORMAT.format(get_filename(fastq), batch))
            batch_files[batch] = filename
            #if filename not in batchfiles:
            #    batchfiles.append(filename)
            f = open(filename, 'a')
            f.write(record.format('fasta'))
            f.close()

        self.plot_bin_distribution(reads_per_batch, get_filename(fastq), "Reads")
        self.plot_bin_distribution(bases_per_batch, get_filename(fastq), "Nucleotides")
        #batch_files.sort(key=utils.get_bin_number)
        return [str(i) or None for i in batch_files] #batch_files

    def plot_bin_distribution(self, counts_per_bin, filename, what_of="Reads"):
        bins = [x for x in range(len(counts_per_bin))]

        fig, ax = plt.subplots(figsize=(18, 12))

        ax.set_ylabel('{} counts'.format(what_of), size=20)
        ax.set_title('{} Counts per Bin'.format(what_of), size=20)
        ax.set_xlabel('Bins', size=20)
        x = np.arange(len(bins)) * 2
        ax.set_xticks(x)
        plt.setp(ax.get_yticklabels(), fontsize=15)
        ax.set_xticklabels(bins, size=15)
        bar_width = 1
        pps = ax.bar(x - bar_width / 2, counts_per_bin, bar_width, label='{} counts'.format(what_of))
        for p in pps:
            height = p.get_height()
            ax.annotate('{}'.format(height),
                        xy=(p.get_x() + p.get_width() / 2, height),
                        xytext=(0, 4),  # 4 points vertical offset
                        textcoords="offset points",
                        ha='center', va='bottom', rotation=90, fontsize=15)

        fig_name = unique_path(os.path.join(self.common.fig_dir, "fig_{}_per_bins_{}.png".format(what_of.lower(), filename)))
        plt.savefig(fig_name)
