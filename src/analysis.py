import os
import sys
import numpy as np
import pandas as pd
import statistics
from math import sqrt
from matplotlib import pyplot as plt
import utils


# CG content

# TODO not needed?
def select_cg_range(value):
    for border in range(10, 100 + 1, 10):
        if value <= border:
            return border


# TODO not needed?
def select_label(value):
    if value == 0 or value == 10:
        return '0-10% CG'
    return '{}-{}% CG'.format(str(value - 9), str(value))


def strand_bias_bin(bias, bin_size, max=100):
    bin = bias // bin_size
    if bin > max:
        return max
    return bin


class Analysis:
    def __init__(self, file=None, output_dir='', start_k=5, end_k=10, threads=1, whisker=5):
        self.start_k = start_k
        self.end_k = end_k
        self.filepath = file
        self.filename = utils.get_filename(file)
        self.out_dir = os.path.join(output_dir, 'sbat')
        self.fig_dir = os.path.join(output_dir, 'sbat', 'figures')
        self.dump_dir = os.path.join(output_dir, 'sbat', 'dump') # output files, mer.jf, dfs, removed in the end unless --kep-computations
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
        fwd_kmers, bwd_kmers = utils.split_forwards_and_backwards(k)

        # remove backward group from DataFrame, as it is already represented as reverse complement to some
        # other k-mer in the DataFrame
        jellyfish_forward = jellyfish_data.drop(bwd_kmers, errors="ignore").set_index("rev_complement")

        # join forward DF with original one on index (equal to forward sequence) in order to connect info about
        # forward and backward datasets
        jellyfish_data = jellyfish_data.reset_index().join(jellyfish_forward, on="index", rsuffix="_", lsuffix="").drop(
            columns=["seq_", "index"], axis=1).dropna()

        if len(jellyfish_data.index) == 0:
            return None

        jellyfish_data.rename(columns={"seq_count_": "rev_complement_count"}, inplace=True)

        # calculate ratio of forward and backward k-mer frequencies
        jellyfish_data["ratio"] = jellyfish_data.apply(
            lambda row: utils.get_ratio(row["seq_count"], row["rev_complement_count"]), axis=1)
        # calculate deviation from 100% accuracy
        jellyfish_data["strand_bias_%"] = jellyfish_data.apply(lambda row: utils.get_strand_bias_percentage(row["ratio"]),
                                                               axis=1)
        # calculate CG content percentage
        jellyfish_data["CG_%"] = jellyfish_data.apply(lambda row: utils.gc_percentage(row["seq"]), axis=1)
        # sort data by bias in descending order
        jellyfish_data = jellyfish_data.sort_values(by=["strand_bias_%"], ascending=False)

        filename = utils.unique_path("df_{}.csv".format(os.path.basename(path.split(".")[0])))
        #if sb_analysis:
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

    def plot_cg_from_dataframe(self, dfs):
        if all(x is None for x in dfs):
            return
        upper_cg = []
        upper_biases = []
        lower_cg = []
        lower_biases = []
        kmers = []
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(15, 15))

        for i, df in enumerate(dfs):
            if df is None or len(utils.get_n_percent(df, self.whisker).index) == 0:
                continue
            kmers.append(i + self.start_k)
            df_head = utils.get_n_percent(df, self.whisker)

            upper_cg.append(df_head["CG_%"].mean().round(2))
            upper_biases.append(df_head["strand_bias_%"].mean().round(2))

            df_tail = utils.get_n_percent(df, self.whisker, True)
            lower_cg.append(None if len(df_head.index) == 0 else df_tail["CG_%"].mean().round(2))
            lower_biases.append(None if len(df_head.index) == 0 else df_tail["strand_bias_%"].mean().round(2))

        x_label = 'Mean CG content [%]'
        y_label = 'Mean Strand bias [%]'
        for ax in [ax1, ax2, ax3]:
            ax.set_xlabel(x_label)
            ax.set_ylabel(y_label)

        ax1.set_title("CG content vs Strand Bias in top " + str(self.whisker) + "% of SB score")
        ax1.scatter(upper_cg, upper_biases, marker="^", color="red")

        ax2.set_title("CG content vs Strand Bias in bottom " + str(self.whisker) + "% of SB score")
        ax2.scatter(lower_cg, lower_biases, marker="v", color="green")

        ax3.set_title("CG content vs Strand Bias in bottom and top " + str(self.whisker) + "% of SB score")
        ax3.scatter(lower_cg, lower_biases, marker="v", color="green")
        ax3.scatter(upper_cg, upper_biases, marker="^", color="red")

        for i, txt in enumerate(kmers):
            if dfs[i] is None:
                continue
            ax1.annotate(" " + str(txt), (upper_cg[i], upper_biases[i]), fontsize=15)
            ax2.annotate(" " + str(txt), (lower_cg[i], lower_biases[i]), fontsize=15)
        fig_path = os.path.join(self.fig_dir, "fig_cg_{0}%_{1}.png".format(str(self.whisker), self.filename))
        fig_path = utils.unique_path(fig_path)
        plt.savefig(fig_path)

    def plot_conf_interval_graph(self, dataframes, k='', start_index=0): # TODO fix args
        if k is None or k=="":
            x = [x + start_index for x in range(len(dataframes))]
        else:
            x = [x for x in range(len(dataframes))]
        plt.figure(figsize=(15, 7))
        plt.xticks(x, x)
        if k is None or k=="":
            plt.title('Confidence Interval among Bins')
        else:
            plt.title('Confidence Interval for K={}'.format(k))
        y = []

        plt.ylabel("Strand bias [%]")
        if k is None or k == '':
            plt.xlabel("K")
        else:
            plt.xlabel("Bins")

        for index, df in enumerate(dataframes):
            if df is None or df.shape[0] < 3:
                continue

            if k is None or k == "":
                index = start_index + index

            mean, ci = plot_confidence_interval(index, df['strand_bias_%'])
            y.append(mean)

        if len(x) > 1 and len(y) > 1:
            try:
                z = np.polyfit(x, y, 3)  # polynomial fit
                p = np.poly1d(z)
                plt.plot(x, p(x), 'r--')
            except Exception as e:
                print("Error occurred during fitting linear regression: {}\nskipping...".format(e))
        print(k)
        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_ci_{0}_{1}.png'.format(self.filename, k)))
        print(fig_name)
        plt.savefig(fig_name)
        plt.close()

    def draw_basic_stats_lineplot(self, name, statfile, k=None, x_axis='k'):
        df = pd.read_csv(statfile)

        # if k is set and x_axis=batch, plotting nanopore split data
        if k is not None:
            df = df.loc[df['k'] == k]

        # Plot a simple line chart
        plt.figure()
        plt.title('Mean and Median of Strand Bias in Nanopore Data')
        plt.ylabel("Strand bias")
        if x_axis == 'k':
            plt.xlabel("K")
        else:
            plt.xlabel("Bins")
        plt.plot(df[x_axis], df['bias_mean'], color='b', label='Mean value of strand bias')
        plt.plot(df[x_axis], df['bias_median'], color='g', label='Median value of strand bias')
        #plt.plot(df[x_axis], df['bias_modus'], color='r', label='Mode value of strand bias')

        plt.legend()
        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}_{1}.png'.format(name, k)))
        plt.savefig(fig_name)
        plt.close()

    def fill_sb_analysis_from_df(self, df, k, batch):
        filename = self.filename
        bias_mean = df['strand_bias_%'].mean().round(2)
        bias_median = df['strand_bias_%'].median().round(2)
        bias_modus = df['strand_bias_%'].mode().iloc[0].round(2)
        percentile_5 = round(df['strand_bias_%'].quantile(0.05), 2)
        percentile_95 = round(df['strand_bias_%'].quantile(0.95), 2)

        import csv
        stat = [filename, k, batch, bias_mean, bias_median, bias_modus, percentile_5, percentile_95]
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

        plt.figure(figsize=(35, 10))
        #plt.show()
        plt.title("K-mers of length {} vs strand bias".format(k))
        plt.xlabel('K-mers')
        plt.ylabel('Strand bias [%]')
        ax = plt.scatter(kmers, bias, marker="o", color="green", s=6)
        if k > 5:
            ax.axes.xaxis.set_ticks([])
        else:
            plt.xticks(range(len(df["seq"])), kmers, rotation=90, fontsize=3.5)

        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_kmer_vs_bias_{0}_k{1}.png'.format(self.filename, k)))
        print(fig_name)
        plt.savefig(fig_name, dpi=400)
        plt.close()

#
# def track_most_common_kmer_change(dfs, k):
#     fwds, bwds = utils.split_forwards_and_backwards(k)
#     fwds = fwds.sort()
#     kmer_changes = pd.DataFrame(
#         data={'seq': fwds},
#         columns=['seq'])
#     for batch, df in enumerate(dfs):
#         df = df.sort_values(by=['seq'])[['seq', 'strand_bias_%']]
#         kmer_changes = kmer_changes.merge(df, how='right', on='seq', suffixes=("", "_batch_{}".format(batch)))
#
#     df.rename({'strand_bias_%': 'strand_bias_%_batch_0'}, inplace=True)
#     return kmer_changes
#
#
# a = pd.DataFrame(
#     data={'seq': ['a', 'c', 'g', 't']},
#           columns=['seq'])
# b = pd.DataFrame(
#     data={'seq': ['a', 'c', 'g', 't'], 'sb': [8,2,3,4]},
#           columns=['seq', 'sb'])
#
# c = pd.DataFrame(
#     data={'seq': ['a', 'c', 'g', 't'], 'sb': [8,25,3,4]},
#           columns=['seq', 'sb'])
#
#
# dfs = []
# for i in range(0, 24):
#     df = pd.read_csv('D:\Alex\School\sbapr\strand-bias-analysis-tool\out\sbat\dump\df_output_6_nanopore_41ec7eae4495de82ef09c416dbd6d5c983ff8c4e_batch_{}_1.csv'.format(i))
#     dfs.append(df)
# #df7 = pd.read_csv('D:\Alex\School\sbapr\strand-bias-analysis-tool\out\sbat\dump\df_output_5_nanopore_41ec7eae4495de82ef09c416dbd6d5c983ff8c4e_batch_2.csv')
#
#
# r = track_most_common_kmer_change(dfs, 5)
# r['diff'] = abs(r['strand_bias_%_batch_23'] - r['strand_bias_%'])
#
# r.sort_values(by=['diff'], ascending=False, inplace=True)
#
# print(r)
# filtered = r[['seq', 'diff', 'strand_bias_%_batch_23', 'strand_bias_%']]
# print(filtered)
# filtered.to_csv("filtered")

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

def select_more_frequent(row, seq=False):
    if row["seq_count"] > row["rev_complement_count"]:
        if seq:
            return row["seq"]
        return row["seq_count"]
    else:
        if seq:
            return row["rev_complement"]
        return row["rev_complement_count"]
