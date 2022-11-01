import math
import os
import sys

import numpy as np
import pandas as pd
import statistics
from math import sqrt
from matplotlib import pyplot as plt
from sbat import utils


class Analysis:
    """
    Class representing primary analysis.
    """
    def __init__(self, file=None, output_dir='', start_k=5, end_k=10, threads=1, margin=5, interactive=False):
        self.start_k = start_k
        self.end_k = end_k
        self.filepath = file
        self.filename = utils.get_filename(file)  # string used in further naming of the created files and figures
        self.out_dir = os.path.join(output_dir)
        self.fig_dir = os.path.join(output_dir, 'figures')
        # temp files, dfs; will be removed in the end unless --kep-computations
        self.dump_dir = os.path.join(output_dir, 'dump')
        self.sb_analysis_file = None
        self.np_sb_analysis_file = None
        self.margin = margin
        self.threads = threads
        self.keep_computations = False
        self.interactive = interactive
        self.plotter = None

    def set_file(self, file):
        """
        Method to init input file and file identificator names in Analysis.
        :param file: path passed as input
        """
        self.filepath = file
        self.filename = utils.get_filename(file)

    def set_output(self, dir):
        """
        Init creation of directories based on setting of output dir.
        :param dir: output dir specified by user or default 'sbat_out'
        """
        self.out_dir = os.path.join(dir, self.out_dir)
        self.fig_dir = os.path.join(dir, self.fig_dir)
        self.dump_dir = os.path.join(dir, self.dump_dir)
        utils.is_or_create_dir(self.out_dir)
        utils.is_or_create_dir(self.fig_dir)
        utils.is_or_create_dir(self.dump_dir)

    def jellyfish_to_dataframe(self, path, k, bin=None):
        """
        Function to create dataframe with statistics based on Jellyfish output
        :param path: path to Jellyfish output
        :param k: length of kmers in given file
        :param bin: bin number
        """
        seq, seq_count = utils.parse_fasta(path)
        if len(seq) == 0:
            sys.exit("no data parsed from {}".format(path))

        # create dataframe with k-mers and their counts
        jellyfish_data = pd.DataFrame(
            data={'seq': seq,
                  'seq_count': seq_count},
            index=seq,
            columns=['seq', 'seq_count'])

        # add column with reverse complements
        jellyfish_data['rev_complement'] = jellyfish_data.apply(lambda row: utils.get_reverse_complement(row["seq"]),
                                                                axis=1)
        # split sequence set into forward and backward sequences (so that k-mer and its reverse complement
        # are not together in group)
        fwd_kmers, bwd_kmers = utils.split_forwards_and_backwards(seq)

        # remove forward group from DataFrame, to gain backward group that will be mapped to mix of k-mers
        jellyfish_backward = jellyfish_data.drop(fwd_kmers, errors="ignore").set_index("rev_complement")

        # join forward DF with original one on index (equal to forward sequence) in order to connect info about
        # forward and backward datasets
        jellyfish_data = jellyfish_data.reset_index().join(jellyfish_backward, on="index", rsuffix="_", lsuffix="").drop(
            columns=["seq_", "index"], axis=1).dropna()

        if len(jellyfish_data.index) != 0:
            jellyfish_data.rename(columns={"seq_count_": "rev_complement_count"}, inplace=True)
            # calculate ratio of forward and backward k-mer frequencies
            jellyfish_data["ratio"] = jellyfish_data.apply(
                lambda row: utils.get_ratio(row["seq_count"], row["rev_complement_count"]), axis=1)
            # calculate deviation from 100% accuracy
            jellyfish_data["strand_bias_%"] = jellyfish_data.apply(
                lambda row: utils.get_strand_bias_percentage(row["ratio"]),
                axis=1)
            # calculate GC content percentage
            jellyfish_data["GC_%"] = jellyfish_data.apply(lambda row: utils.gc_percentage(row["seq"]), axis=1)
            # sort data by bias in descending order
            jellyfish_data = jellyfish_data.sort_values(by=["strand_bias_%"], ascending=False)

        filename = utils.unique_path("df_{}.csv".format(os.path.basename(path.split(".")[0])))
        # if sb_analysis:
        self.fill_sb_analysis_from_df(jellyfish_data, k, bin)
        if self.keep_computations:
            filename = utils.unique_path(os.path.join(self.dump_dir, filename))
            jellyfish_data.to_csv(filename, index=False)

        return jellyfish_data

    def init_analysis(self, nanopore=False):
        analysis = pd.DataFrame(
            data={},
            index=None,
            columns=['file', 'k', 'bin', 'bias_mean', 'bias_median', 'bias_modus', 'percentile_5', 'percentile_95'])

        if nanopore:
            analysis_name = utils.unique_path(os.path.join(self.out_dir, 'np_sb_analysis_' + self.filename + '.csv'))
        else:
            analysis_name = utils.unique_path(os.path.join(self.out_dir, 'sb_analysis_' + self.filename + '.csv'))
            self.sb_analysis_file = analysis_name
        print("analysis stored in: {}".format(analysis_name))
        analysis.to_csv(analysis_name, index=False)
        return analysis_name


    def calculate_gc_plot_data(self, dfs):
        if all(x is None for x in dfs):
            return None
        data = CalculatedGCData()

        for i, df in enumerate(dfs):
            if df is None or len(utils.get_n_percent(df, self.margin).index) == 0:
                # skip DF if it's None or has too little values for retrieving N percent
                continue
            data.kmers.append(i + self.start_k)
            df_head = utils.get_n_percent(df, self.margin)  # get N percent with the highest bias
            data.upper_gc.append(None if len(df_head.index) == 0 else round(df_head["GC_%"].mean(), 2))
            data.upper_biases.append(None if len(df_head.index) == 0 else round(df_head["strand_bias_%"].mean(), 2))

            df_tail = utils.get_n_percent(df, self.margin, True)  # get N percent with the lowest bias
            data.lower_gc.append(None if len(df_head.index) == 0 else round(df_tail["GC_%"].mean(), 2))
            data.lower_biases.append(None if len(df_head.index) == 0 else round(df_tail["strand_bias_%"].mean(), 2))

        if not data.kmers:  # no dataframe is big enough to provide data
            return None
        return data

    def plot_gc_from_dataframe(self, data):
        if data is None:
            return

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18, 25))
        x_label = 'Mean GC content [%]'
        y_label = 'Mean Strand bias [%]'

        for ax in [ax1, ax2, ax3]:
            ax.set_xlabel(x_label, fontsize=27)
            ax.set_ylabel(y_label, fontsize=27)
            ax.xaxis.set_tick_params(labelsize=25)
            ax.yaxis.set_tick_params(labelsize=25)

        ax1.set_title("GC Content vs Strand Bias in Top " + str(self.margin) + "% of SB Score", fontsize=30)
        ax1.scatter(data.upper_gc, data.upper_biases, marker="^", color="red", s=300,
                    label="SB in bottom {}% for given k".format(self.margin))

        ax2.set_title("GC content vs Strand Bias in Bottom " + str(self.margin) + "% of SB Score", fontsize=30)
        ax2.scatter(data.lower_gc, data.lower_biases, marker="v", color="green", s=300,
                    label="SB in bottom {}% for given k".format(self.margin))

        ax3.set_title("GC content vs Strand Bias in Bottom and Top " + str(self.margin) + "% of SB Score", fontsize=30)
        ax3.scatter(data.lower_gc, data.lower_biases, marker="v", color="green", s=300,
                    label="SB in bottom {}% for given k".format(self.margin))
        ax3.scatter(data.upper_gc, data.upper_biases, marker="^", color="red", s=300,
                    label="SB in bottom {}% for given k".format(self.margin))
        plt.legend()
        for i, txt in enumerate(data.kmers):
            ax1.annotate(" " + str(txt), (data.upper_gc[i], data.upper_biases[i]), fontsize=24)
            ax2.annotate(" " + str(txt), (data.lower_gc[i], data.lower_biases[i]), fontsize=24)
        fig_path = os.path.join(self.fig_dir, "fig_gc_{0}%_{1}.png".format(str(self.margin), self.filename))
        fig_path = utils.unique_path(fig_path)
        plt.tight_layout(pad=1.5)
        plt.savefig(fig_path)

    def plot_conf_interval_graph(self, dataframes, k='', start_index=0, nanopore=False):
        if not nanopore:  # CI plot for different k sizes
            x = [x for x in range(min(self.start_k, 5), max(self.end_k + 1, 11))]
        else:  # CI plot for bins
            x = [x for x in range(len(dataframes))]

        plt.figure(figsize=(12, 8))
        plt.xticks(x, x, fontsize=20)
        plt.yticks(fontsize=20)
        plt.ylabel("Strand bias [%]", fontsize=22)

        if not nanopore:
            plt.title('Confidence Interval among Different Sizes of K', fontsize=25)
            plt.xlabel("K", fontsize=22)
        else:
            plt.title('Confidence Interval among Bins for K={}'.format(k), fontsize=25)
            plt.xlabel("Bins", fontsize=22)

        y = []
        for index, df in enumerate(dataframes):
            legend = True if index == 0 else False
            if df is None or df.shape[0] < 3:
                continue
            if k is None or k == "":
                index = start_index + index

            mean, ci = plot_confidence_interval(index, df['strand_bias_%'], legend=legend)
            y.append(round(mean, 3))

        plt.xlim(x[0] - 0.5, x[-1] + 0.5)
        plt.ylim(0, min(100, max(y) + 5))

        if len(x) == len(y) and len(x) > 1:  # there cannot be gaps in data (bins)
            try:
                z = np.polyfit(x, y, 3)  # polynomial fit
                p = np.poly1d(z)
                plt.plot(x, p(x), 'r--', label="polynomial fit of 3rd degree", linewidth=2)
                plt.legend(fontsize=18)
            except Exception as e:  # in case of a problem, just skip the polynomial fit and continue
                print("Error occurred during polynomial fitting: {}\nskipping...".format(e))
        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_ci_{0}_{1}.png'.format(self.filename, k)))
        plt.savefig(fig_name)
        plt.close()

    def plot_basic_stats_lineplot(self, name, statfile, k=None, nanopore=False):
        df = pd.read_csv(statfile)
        if k is not None:
            df = df.loc[df['k'] == k].dropna()

        if df.shape[0] <= 1:
            # no point in drawing this plot for 0 or 1 record
            return

        # Plot a simple line chart
        fig = plt.figure(figsize=(9, 6))
        plt.ylabel("Strand bias", fontsize=16)
        if nanopore:  # plotting SB x Bins for specific size of K (nanopore analysis)
            x_axis = "bin"
            plt.title('Mean and Median of Strand Bias for K={}'.format(k), fontsize=18)
            plt.xlabel("Bins", fontsize=16)
            fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}_k{1}.png'.format(name, k)))
            bin_count = df[x_axis].max() + 1
            # print label of each bin, if there is too much of them, print only every count // 10
            tick_labels = range(0, df[x_axis].max() + 1)[::((bin_count//20) + 1)]
            plt.xticks(tick_labels, fontsize=14)

        else:  # plotting SB x sizes of K (primary analysis)
            x_axis = "k"
            plt.title('Mean and Median of Strand Bias', fontsize=20)
            plt.xlabel("K", fontsize=16)
            fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}.png'.format(name)))
            plt.xticks(df[x_axis], fontsize=14)
        plt.yticks(fontsize=14)
        plt.plot(df[x_axis], df['bias_mean'], '-bo', label='Mean value of strand bias')
        plt.plot(df[x_axis], df['bias_median'], '-go', label='Median value of strand bias')

        plt.legend(fontsize=14)
        plt.savefig(fig_name)
        plt.close()

    def fill_sb_analysis_from_df(self, df, k, bin):
        if df is None or df.empty:
            bias_mean, bias_median, bias_modus, percentile_5, percentile_95 = [math.nan for _ in range(5)]
        else:
            bias_mean = round(df['strand_bias_%'].mean(), 2)
            bias_median = round(df['strand_bias_%'].median(), 2)
            bias_modus = round(df['strand_bias_%'].mode().iloc[0], 2)
            percentile_5 = round(df['strand_bias_%'].quantile(0.05), 2)
            percentile_95 = round(df['strand_bias_%'].quantile(0.95), 2)

        import csv
        stat = [self.filename, k, bin, bias_mean, bias_median, bias_modus, percentile_5, percentile_95]
        if bin is None:
            sb = self.sb_analysis_file
        else:
            sb = self.np_sb_analysis_file
        with open(sb, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(stat)

    def plot_kmers_vs_bias(self, df, k):
        # get more frequent out of k-mer and its rev. complement
        df["more_freq_count"] = df.apply(lambda row: utils.select_more_frequent(row), axis=1)
        df = df.sort_values(by=['more_freq_count'], ascending=False)
        kmers = df["seq"]
        bias = df["strand_bias_%"]

        plt.figure(figsize=(26, 10))
        plt.margins(x=0.01)
        plt.title("K-mers of Length {} vs Strand Bias".format(k), fontsize=28)
        plt.xlabel('K-mers', fontsize=25)
        plt.ylabel('Strand Bias [%]', fontsize=25)
        plt.yticks(fontsize=18)
        ax = plt.scatter(kmers, bias, marker="o", color="green", s=6, label="specific k-mer")
        plt.legend(fontsize=20)
        if k > 5:  # no reason in printing k-mers, there is too much of them to fit for k > 5
            ax.axes.xaxis.set_ticks([])
        else:
            plt.xticks(range(len(kmers)), kmers, rotation=90, fontsize=3.5)

        fig_name = utils.unique_path(
            os.path.join(self.fig_dir, 'fig_kmer_vs_bias_{0}_k{1}.png'.format(self.filename, k)))
        plt.savefig(fig_name, dpi=250)
        plt.close()

    def track_most_common_kmer_change_freq(self, dfs, k):
        if dfs is None or len(dfs) < 1:
            return
        fwds, bwds = utils.split_forwards_and_backwards(dfs[0]["seq"])
        fwds.sort()
        kmer_changes = pd.DataFrame(
            data={'seq': fwds},
            columns=['seq'])
        last_bin = len(dfs) - 1
        valid_bins = len(dfs)

        for bin, df in enumerate(dfs):
            if df is None or df.empty:
                valid_bins -= 1
                continue
            df["freq_count"] = df.apply(lambda row: utils.select_more_frequent(row), axis=1)
            df = df.sort_values(by=['seq'])[['seq', 'freq_count', 'GC_%', 'strand_bias_%',]]
            if bin != 0:
                df = df[["seq", "strand_bias_%"]]
            kmer_changes = kmer_changes.merge(df, how='right', on='seq', suffixes=("", "_bin_{}".format(bin))).dropna()
            last_bin = bin

        if valid_bins < 1:
            return
        if valid_bins == 1:
            filename = utils.unique_path(os.path.join(self.out_dir, "kmer_freq_k{}_{}.csv".format(k, self.filename)))
            kmer_changes = kmer_changes.sort_values(by=['freq_count'], ascending=False)
            kmer_changes.to_csv(filename, index=None)
            return

        kmer_changes.rename(columns={'strand_bias_%': 'strand_bias_%_bin_0'}, inplace=True)
        kmer_changes = kmer_changes.sort_values(by=['freq_count'], ascending=False)
        kmer_changes['diff'] = abs(kmer_changes['strand_bias_%_bin_{}'.format(last_bin)] - kmer_changes['strand_bias_%_bin_0'])
        kmer_changes = kmer_changes.round(3)
        filename = utils.unique_path(os.path.join(self.out_dir, "kmer_diff_k{}_{}.csv".format(k, self.filename)))
        kmer_changes.to_csv(filename, index=None)


def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.2, legend=False):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / sqrt(len(values))
    left = x - horizontal_line_width / 2
    top = mean - confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean + confidence_interval
    plt.plot([x, x], [top, bottom], color=color, linewidth=2, label="Confidence Level 95%" if legend else "")
    plt.plot([left, right], [top, top], color=color, linewidth=2)
    plt.plot([left, right], [bottom, bottom], color=color, linewidth=2)
    plt.plot(x, mean, 'o', color='#f44336', linewidth=2, label="Mean Strand Bias" if legend else "")

    return mean, confidence_interval


class CalculatedGCData:
    def __init__(self):
        self.upper_gc = []
        self.upper_biases = []
        self.lower_gc = []
        self.lower_biases = []
        self.kmers = []