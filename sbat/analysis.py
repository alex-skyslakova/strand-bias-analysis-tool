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
    def __init__(self, file=None, output_dir='', start_k=5, end_k=10, threads=1, whisker=5):
        self.start_k = start_k
        self.end_k = end_k
        self.filepath = file
        self.filename = utils.get_filename(file)
        self.out_dir = os.path.join(output_dir)
        self.fig_dir = os.path.join(output_dir, 'figures')
        # 'output' files, dfs; will be removed in the end unless --kep-computations
        self.dump_dir = os.path.join(output_dir, 'dump')
        self.sb_analysis_file = None
        self.np_sb_analysis_file = None
        self.whisker = whisker
        self.threads = threads
        self.keep_computations = False

    def set_file(self, file):
        self.filepath = file
        self.filename = utils.get_filename(file)

    def set_output(self, dir):
        self.out_dir = os.path.join(dir, self.out_dir)
        self.fig_dir = os.path.join(dir, self.fig_dir)
        self.dump_dir = os.path.join(dir, self.dump_dir)
        utils.is_or_create_dir(self.out_dir)
        utils.is_or_create_dir(self.fig_dir)
        utils.is_or_create_dir(self.dump_dir)

    def jellyfish_to_dataframe(self, path, k, batch=None):
        """
        Function to create dataframe with statistics based on Jellyfish output
        :param path: path to Jellyfish output
        :param k: length of kmers in given file
        :param batch: bin number
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
        self.fill_sb_analysis_from_df(jellyfish_data, k, batch)
        if self.keep_computations:
            filename = utils.unique_path(os.path.join(self.dump_dir, filename))
            jellyfish_data.to_csv(filename, index=False)

        return jellyfish_data

    def init_analysis(self, nanopore=False):
        analysis = pd.DataFrame(
            data={},
            index=None,
            columns=['file', 'k', 'batch', 'bias_mean', 'bias_median', 'bias_modus', 'percentile_5', 'percentile_95'])

        if nanopore:
            analysis_name = utils.unique_path(os.path.join(self.out_dir, 'np_sb_analysis_' + self.filename + '.csv'))
        else:
            analysis_name = utils.unique_path(os.path.join(self.out_dir, 'sb_analysis_' + self.filename + '.csv'))
            self.sb_analysis_file = analysis_name
        print("analysis stored in: {}".format(analysis_name))
        analysis.to_csv(analysis_name, index=False)
        return analysis_name

    def plot_gc_from_dataframe(self, dfs):
        if all(x is None for x in dfs):
            return
        upper_gc = []
        upper_biases = []
        lower_gc = []
        lower_biases = []
        kmers = []
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(18, 22))

        for i, df in enumerate(dfs):
            if df is None or len(utils.get_n_percent(df, self.whisker).index) == 0:
                # skip DF if it's None or has too little values for retrieving N percent
                continue
            kmers.append(i + self.start_k)

            df_head = utils.get_n_percent(df, self.whisker)  # get N percent with the highest bias
            upper_gc.append(df_head["GC_%"].mean().round(2))
            upper_biases.append(df_head["strand_bias_%"].mean().round(2))

            df_tail = utils.get_n_percent(df, self.whisker, True)  # get N percent with the lowest bias
            lower_gc.append(None if len(df_head.index) == 0 else df_tail["GC_%"].mean().round(2))
            lower_biases.append(None if len(df_head.index) == 0 else df_tail["strand_bias_%"].mean().round(2))

        if not kmers:  # no dataframe is big enough to provide data
            return

        x_label = 'Mean GC content [%]'
        y_label = 'Mean Strand bias [%]'

        for ax in [ax1, ax2, ax3]:
            ax.set_xlabel(x_label, fontsize=18)
            ax.set_ylabel(y_label, fontsize=18)
            ax.xaxis.set_tick_params(labelsize=18)
            ax.yaxis.set_tick_params(labelsize=18)

        ax1.set_title("CG Content vs Strand Bias in Top " + str(self.whisker) + "% of SB Score", fontsize=20)
        ax1.scatter(upper_cg, upper_biases, marker="^", color="red")

        ax2.set_title("CG content vs Strand Bias in Bottom " + str(self.whisker) + "% of SB Score", fontsize=20)
        ax2.scatter(lower_cg, lower_biases, marker="v", color="green")

        ax3.set_title("CG content vs Strand Bias in Bottom and Top " + str(self.whisker) + "% of SB Score", fontsize=20)
        ax3.scatter(lower_cg, lower_biases, marker="v", color="green")
        ax3.scatter(upper_cg, upper_biases, marker="^", color="red")

        for i, txt in enumerate(kmers):
            ax1.annotate(" " + str(txt), (upper_cg[i], upper_biases[i]), fontsize=15)
            ax2.annotate(" " + str(txt), (lower_cg[i], lower_biases[i]), fontsize=15)
        fig_path = os.path.join(self.fig_dir, "fig_cg_{0}%_{1}.png".format(str(self.whisker), self.filename))
        fig_path = unique_path(fig_path)
        plt.savefig(fig_path)

    def plot_conf_interval_graph(self, dataframes, k='', start_index=0, nanopore=False):  # TODO fix args
        if not nanopore:
            x = [x for x in range(min(self.start_k, 5), max(self.end_k + 1, 11))]
        else:
            x = [x for x in range(len(dataframes))]

        plt.figure(figsize=(15, 7))
        plt.xticks(x, x, fontsize=12)
        plt.yticks(fontsize=12)

        if not nanopore:
            plt.title('Confidence Interval among Different Sizes of K', fontsize=18)
        else:
            plt.title('Confidence Interval among Bins for K={}'.format(k), fontsize=18)
        y = []

        plt.ylabel("Strand bias [%]", fontsize=18)
        if not nanopore:
            plt.xlabel("K", fontsize=18)
        else:
            plt.xlabel("Bins", fontsize=18)

        for index, df in enumerate(dataframes):
            if df is None or df.shape[0] < 3:
                continue

            if k is None or k == "":
                index = start_index + index

            mean, ci = plot_confidence_interval(index, df['strand_bias_%'])
            y.append(round(mean, 3))

        plt.xlim(x[0] - 0.5, x[-1] + 0.5)
        plt.ylim(0, min(100, max(y) + 5))

        if len(x) > 1 and len(y) > 1:
            try:
                z = np.polyfit(x, y, 3)  # polynomial fit
                p = np.poly1d(z)
                plt.plot(x, p(x), 'r--', label="polynomial fit of 3rd degree")
                plt.legend()
            except Exception as e:
                print("Error occurred during polynomial fitting: {}\nskipping...".format(e))
        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_ci_{0}_{1}.png'.format(self.filename, k)))
        plt.savefig(fig_name)
        plt.close()

    def draw_basic_stats_lineplot(self, name, statfile, k=None, nanopore=False):
        df = pd.read_csv(statfile)
        if k is not None:
            df = df.loc[df['k'] == k].dropna()

        if df.shape[0] <= 1:
            # no point in drawing this plot for 0 or 1 record
            return

        # Plot a simple line chart
        plt.figure()
        plt.ylabel("Strand bias", fontsize=14)
        if nanopore:  # plotting SB x Bins for specific size of K (nanopore analysis)
            x_axis = "batch"
            plt.title('Mean and Median of Strand Bias for K={}'.format(k))
            plt.xlabel("Bins", fontsize=14)
            fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}_k{1}.png'.format(name, k)))
            plt.xticks(df[x_axis], fontsize=12)
        else:  # plotting SB x sizes of K (primary analysis)
            x_axis = "k"
            plt.title('Mean and Median of Strand Bias')
            plt.xlabel("K", fontsize=14)
            fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}.png'.format(name)))
            plt.xticks(df[x_axis], fontsize=12)
        plt.yticks(fontsize=12)
        plt.plot(df[x_axis], df['bias_mean'], '-bo', label='Mean value of strand bias')
        plt.plot(df[x_axis], df['bias_median'], '-go', label='Median value of strand bias')

        plt.legend()
        plt.savefig(fig_name)
        plt.close()

    def fill_sb_analysis_from_df(self, df, k, batch):
        if df is None or df.empty:
            bias_mean, bias_median, bias_modus, percentile_5, percentile_95 = [math.nan for _ in range(5)]
        else:
            bias_mean = df['strand_bias_%'].mean().round(2)
            bias_median = df['strand_bias_%'].median().round(2)
            bias_modus = df['strand_bias_%'].mode().iloc[0].round(2)
            percentile_5 = round(df['strand_bias_%'].quantile(0.05), 2)
            percentile_95 = round(df['strand_bias_%'].quantile(0.95), 2)

        import csv
        stat = [self.filename, k, batch, bias_mean, bias_median, bias_modus, percentile_5, percentile_95]
        if batch is None:
            sb = self.sb_analysis_file
        else:
            sb = self.np_sb_analysis_file
        with open(sb, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(stat)

    def plot_kmers_vs_bias(self, df, k):
        df["more_freq_count"] = df.apply(lambda row: utils.select_more_frequent(row), axis=1)
        df = df.sort_values(by=['more_freq_count'], ascending=False)
        kmers = df["seq"]
        bias = df["strand_bias_%"]
        freq = df["more_freq_count"]

        plt.figure(figsize=(35, 10))

        plt.title("K-mers of Length {} vs Strand Bias".format(k), fontsize=25)
        plt.xlabel('K-mers', fontsize=25)
        plt.ylabel('Strand bias [%]', fontsize=25)
        ax = plt.scatter(kmers, bias, marker="o", color="green", s=6)
        if k > 5:
            ax.axes.xaxis.set_ticks([])
        else:
            plt.xticks(range(len(freq)), freq, rotation=90, fontsize=3.5)

        fig_name = utils.unique_path(
            os.path.join(self.fig_dir, 'fig_kmer_vs_bias_{0}_k{1}.png'.format(self.filename, k)))
        plt.savefig(fig_name, dpi=400)
        plt.close()

    def track_most_common_kmer_change_freq(self, dfs, k):
        if dfs == None or len(dfs) < 2:
            return
        fwds, bwds = utils.split_forwards_and_backwards(dfs[0]["seq"])
        fwds.sort()
        kmer_changes = pd.DataFrame(
            data={'seq': fwds},
            columns=['seq'])
        last_bin = len(dfs) - 1
        valid_bins = len(dfs)

        for batch, df in enumerate(dfs):
            if df is None or df.empty:
                valid_bins -= 1
                continue
            df["freq_count"] = df.apply(lambda row: utils.select_more_frequent(row), axis=1)
            df = df.sort_values(by=['seq'])[['seq', 'freq_count', 'GC_%', 'strand_bias_%',]]
            if batch != 0:
                df = df[["seq", "strand_bias_%"]]
            kmer_changes = kmer_changes.merge(df, how='right', on='seq', suffixes=("", "_batch_{}".format(batch))).dropna()
            last_bin = batch

        if valid_bins <= 1:
            return

        kmer_changes.rename(columns={'strand_bias_%': 'strand_bias_%_batch_0'}, inplace=True)
        kmer_changes = kmer_changes.sort_values(by=['freq_count'], ascending=False)
        kmer_changes['diff'] = abs(kmer_changes['strand_bias_%_batch_{}'.format(last_bin)] - kmer_changes['strand_bias_%_batch_0'])
        kmer_changes = kmer_changes.round(3)
        filename = utils.unique_path(os.path.join(self.out_dir, "kmer_diff_k{}_{}.csv".format(k, self.filename)))
        kmer_changes.to_csv(filename, index=None)


def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.25):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / sqrt(len(values))
    left = x - horizontal_line_width / 2
    top = mean - confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean + confidence_interval
    plt.plot([x, x], [top, bottom], color=color)
    plt.plot([left, right], [top, top], color=color)
    plt.plot([left, right], [bottom, bottom], color=color)
    plt.plot(x, mean, 'o', color='#f44336')

    return mean, confidence_interval
