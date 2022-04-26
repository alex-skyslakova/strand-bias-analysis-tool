import os
import Bio.SeqIO as SeqIO
import numpy as np
from dateutil.parser import parse as dparse
import pandas as pd
from pytz import utc

import statistics
from math import sqrt
from matplotlib import pyplot as plt
from collections import Counter
import bisect
import dateutil.rrule as rrule
import datetime

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
        self.sb_analysis_file = None #os.path.join(output_dir, 'sbat', 'sb_analysis_' + utils.get_filename(file) + ".csv")
        self.whisker = whisker
        self.threads = threads

    def set_output(self, dir):
        self.out_dir = os.path.join(dir, self.out_dir)
        self.fig_dir = os.path.join(dir, self.fig_dir)
        self.dump_dir = os.path.join(dir, self.dump_dir)
        utils.is_or_create_dir(self.out_dir)
        utils.is_or_create_dir(self.fig_dir)
        utils.is_or_create_dir(self.dump_dir)

    def init_analysis(self):
        analysis = pd.DataFrame(
            data={},
            index=None,
            columns=['file', 'k', 'batch', 'bias_mean', 'bias_median', 'bias_modus', 'percentile_5', 'percentile_95'])
        analysis_name = utils.unique_path(os.path.join(self.out_dir, 'sb_analysis_' + self.filename + '.csv'))
        print("analysis stored in: {}".format(analysis_name))
        analysis.to_csv(analysis_name, index=False)
        self.sb_analysis_file = analysis_name
        return analysis_name

    def plot_cg_from_dataframe(self):
        upper_cg = []
        upper_biases = []
        lower_cg = []
        lower_biases = []
        kmers = [x for x in range(self.start_k, self.end_k + 1)]
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))

        for k in range(self.start_k, self.end_k + 1):
            df = pd.read_csv(self.dump_dir + "df_output_{}_{}.csv".format(k, self.filename))

            df_head = utils.get_n_percent(df, self.whisker)
            upper_cg.append(df_head["CG_%"].mean().round(2))
            upper_biases.append(df_head["strand_bias_%"].mean().round(2))

            df_tail = utils.get_n_percent(df, self.whisker, True)
            lower_cg.append(df_tail["CG_%"].mean().round(2))
            lower_biases.append(df_tail["strand_bias_%"].mean().round(2))

        ax1.set_title("CG content vs Strand Bias in top " + str(self.whisker) + "% of SB score")
        ax1.set_xlabel('Mean CG content [%]')
        ax1.set_ylabel('Mean Strand bias [deviation from normal value in %]')
        ax1.scatter(upper_cg, upper_biases, marker="^", color="red")

        ax2.set_title("CG content vs Strand Bias in bottom " + str(self.whisker) + "% of SB score")
        ax2.set_xlabel('Mean CG content [%]')
        ax2.set_ylabel('Mean Strand bias [deviation from normal value in %]')
        ax2.scatter(lower_cg, lower_biases, marker="v", color="green")
        for i, txt in enumerate(kmers):
            ax1.annotate(" " + str(txt), (upper_cg[i], upper_biases[i]), fontsize=15)
            ax2.annotate(" " + str(txt), (lower_cg[i], lower_biases[i]), fontsize=15)
        fig_path = os.path.join(self.fig_dir, "fig_cg_{0}%_{1}.png".format(str(self.whisker), self.filename))
        fig_path = utils.unique_path(fig_path)
        plt.savefig(fig_path)

    def plot_conf_interval_graph(self, dataframes, k='', start_index=0): # TODO fix args
        x = [x + start_index for x in range(len(dataframes))]
        plt.figure(figsize=(15, 7))
        plt.xticks(x, x)
        plt.title('Confidence Interval')
        y = []

        plt.ylabel("Strand bias")
        if k is None or k == '':
            plt.xlabel("K")
        else:
            plt.xlabel("Bins")

        for index, df in enumerate(dataframes):
            index = start_index + index
            if df.shape[0] < 3:
                plt.close()
                return
            mean, ci = plot_confidence_interval(index, df['strand_bias_%'])
            y.append(mean)

        if len(x) > 1 and len(y) > 1:
            try:
                z = np.polyfit(x, y, 1)  # polynomial fit
                p = np.poly1d(z)
                plt.plot(x, p(x), 'r--')
            except Exception as e:
                print("Error occurred during fitting linear regression: {}\nskipping...".format(e))

        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_ci_{0}_{1}.png'.format(self.filename, k)))
        plt.savefig(fig_name)
        plt.close()

    def draw_basic_stats_lineplot(self, name, statfile, k, x_axis='k'):
        df = pd.read_csv(statfile)

        # if k is set and x_axis=batch, plotting nanopore split data
        if k is not None:
            df = df.loc[df['k'] == k]

        # Plot a simple line chart
        plt.figure()
        plt.title('Mean, Mode and Median of Strand Bias in Nanopore Data')
        plt.ylabel("Strand bias")
        if x_axis == 'k':
            plt.xlabel("K")
        else:
            plt.xlabel("Bins")
        plt.plot(df[x_axis], df['bias_mean'], color='b', label='Mean value of strand bias')
        plt.plot(df[x_axis], df['bias_median'], color='g', label='Median value of strand bias')
        plt.plot(df[x_axis], df['bias_modus'], color='r', label='Mode value of strand bias')

        plt.legend()
        fig_name = utils.unique_path(os.path.join(self.fig_dir, 'fig_lineplot_{0}_{1}.png'.format(name, k)))
        plt.savefig(fig_name)
        plt.close()

    def fill_sb_analysis_from_df(self, df, k, batch):
        filename = self.filename
        bias_mean = df['strand_bias_%'].mean().round(2)
        bias_median = df['strand_bias_%'].median().round(2)
        bias_modus = df['strand_bias_%'].mode().iloc[0]
        percentile_5 = round(df['strand_bias_%'].quantile(0.05), 2)
        percentile_95 = round(df['strand_bias_%'].quantile(0.95), 2)

        import csv
        stat = [filename, k, batch, bias_mean, bias_median, bias_modus, percentile_5, percentile_95]
        with open(self.sb_analysis_file, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(stat)





# TODO possibly remove
'''
def fill_sb_analysis(output_csv, k, batch, stat_file):
    df = pd.read_csv(output_csv)
    fill_sb_analysis_from_df(output_csv, df, k, batch, stat_file)
'''



def get_binned_sb_info(df):
    pass#df.apply(lambda row: get_strand_bias_percentage(row["ratio"]), axis=1





def plot_confidence_interval(x, values, z=1.96, color='#2187bb', horizontal_line_width=0.25):
    mean = statistics.mean(values)
    stdev = statistics.stdev(values)
    confidence_interval = z * stdev / sqrt(len(values))
    print(mean, stdev)
    left = x - horizontal_line_width / 2
    top = mean - confidence_interval
    right = x + horizontal_line_width / 2
    bottom = mean + confidence_interval
    plt.plot([x, x], [top, bottom], color=color)
    plt.plot([left, right], [top, top], color=color)
    plt.plot([left, right], [bottom, bottom], color=color)
    plt.plot(x, mean, 'o', color='#f44336')

    return mean, confidence_interval

'''
def analysis(analysis_full_path, output_dir, name, start_kmer=5, end_kmer=10, batch_count=49):
    #sb_analysis = init_analysis(analysis_full_path, name.split('.')[0])
    #df_analysis = pd.read_csv(sb_analysis)
    for kmer in range(start_kmer, end_kmer + 1):
        draw_conf_interval_graph(analysis_full_path, output_dir, utils.get_filename(name), kmer, batch_count)
        draw_basic_stats_lineplot(output_dir, utils.get_filename(name),
                                  analysis_full_path, kmer)

'''




'''
x = []
for i in range(49):
    x.append(pd.read_csv(r"Users\alexa\Documents\School\sbapr\strand-bias-analysis-tool\nanopore\output\df_output_8_nanopore_nanopore_GM24385_11_batch_{0}.csv".format(i)))

plot_conf_interval_graph(x, '', 'test_analysis_k8')'''
#analysis("nanopore/batches/", "nanopore/batches/", "nanopore_nanopore_GM24385_11")
