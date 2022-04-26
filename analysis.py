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

NANOPORE_BIN_FORMAT = 'nanopore_{0}_batch_{1}.fasta'
NANOPORE_BIN_DF_FORMAT = 'output_{0}_nanopore_{1}_batch_{2}.csv'

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


def plot_cg_from_dataframe(input_dir, output_dir, filename, percent=5, start_kmer=5, end_kmer=10):
    upper_cg = []
    upper_biases = []
    lower_cg = []
    lower_biases = []
    kmers = [x for x in range(start_kmer, end_kmer + 1)]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))  # plt.figure(figsize=(10,7))

    for k in range(start_kmer, end_kmer + 1):
        df = pd.read_csv(input_dir + "df_output_" + str(k) + "_" + filename + ".csv")

        df_head = get_n_percent(df, percent)
        upper_cg.append(df_head["CG_%"].mean().round(2))
        upper_biases.append(df_head["strand_bias_%"].mean().round(2))

        df_tail = get_n_percent(df, percent, True)
        lower_cg.append(df_tail["CG_%"].mean().round(2))
        lower_biases.append(df_tail["strand_bias_%"].mean().round(2))

    ax1.set_title("CG content vs Strand Bias in top " + str(percent) + "% of SB score")
    ax1.set_xlabel('Mean CG content [%]')
    ax1.set_ylabel('Mean Strand bias [deviation from normal value in %]')
    ax1.scatter(upper_cg, upper_biases, marker="^", color="red")

    ax2.set_title("CG content vs Strand Bias in bottom " + str(percent) + "% of SB score")
    ax2.set_xlabel('Mean CG content [%]')
    ax2.set_ylabel('Mean Strand bias [deviation from normal value in %]')
    ax2.scatter(lower_cg, lower_biases, marker="v", color="green")
    for i, txt in enumerate(kmers):
        ax1.annotate(" " + str(txt), (upper_cg[i], upper_biases[i]), fontsize=15)
        ax2.annotate(" " + str(txt), (lower_cg[i], lower_biases[i]), fontsize=15)
    plt.savefig(output_dir + "fig_cg_" + str(percent) + "%_" + filename + ".png")


def get_n_percent(df, n, tail=False):
    if tail:
        return df.tail(int(len(df) * (n / 100)))
    else:
        return df.head(int(len(df) * (n / 100)))


# Nanopore bias comparation per hours

def bin_nanopore(fastq, interval=1):
    file_type = 'fastq' if fastq.split('.')[-1] == 'fastq' else 'fasta'

    batchfiles = []
    start = utc.localize(datetime.datetime.now())
    end = utc.localize(datetime.datetime(1970, 1, 1, 0, 0, 0))

    for record in SeqIO.parse(fastq, file_type):
        record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
        if record_time < start:
            start = record_time
        if record_time > end:
            end = record_time
    batches = hours_aligned(start, end, interval)
    sequences_per_batch = [0 for _ in range(len(batches))]

    for record in SeqIO.parse(fastq, file_type):
        record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
        batch = bisect.bisect_left(batches, record_time)
        sequences_per_batch[batch] += 1
        filename = utils.unique_path(os.path.join(os.path.dirname(fastq), NANOPORE_BIN_FORMAT.format(utils.get_filename(fastq), batch)))
        if filename not in batchfiles:
            batchfiles.append(filename)
        f = open(filename, 'a')
        f.write(record.format('fasta'))
        f.close()

    plot_bin_distribution(sequences_per_batch, utils.get_filename(fastq))
    return batchfiles


def plot_bin_distribution(sequences_per_bin, filename):
    bins = [x for x in range(len(sequences_per_bin))]

    fig, ax = plt.subplots(figsize=(18, 12))

    ax.set_ylabel('Sequence counts', size=15)
    ax.set_title('Sequence counts per bin - ' + filename, size=15)
    ax.set_xlabel('Bins', size=15)
    x = np.arange(len(bins)) * 2
    ax.set_xticks(x)
    ax.set_xticklabels(bins, size=10)
    bar_width = 1
    pps = ax.bar(x - bar_width / 2, sequences_per_bin, bar_width, label='Sequence counts')
    for p in pps:
        height = p.get_height()
        ax.annotate('{}'.format(height),
                    xy=(p.get_x() + p.get_width() / 2, height),
                    xytext=(0, 4),  # 4 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom', rotation=90)

    fig_name = utils.unique_path("fig_seq_per_bins_{}.png".format(filename))
    plt.savefig(fig_name)


def filter_time(descr, time_from, time_to):
    # extract start time from description
    time = dparse([i for i in descr.split() if i.startswith('start_time')][0].split('=')[1])
    if dparse(time_from) <= time <= dparse(time_to):
        return True
    else:
        return False


def hours_aligned(start, end, interval=1):
    chunks = []
    rule = rrule.rrule(rrule.HOURLY, byminute=0, bysecond=0, dtstart=start, interval=interval)
    for x in rule.between(start, end, inc=False):
        chunks.append(x)
    chunks.append(end)
    return chunks

# TODO possibly remove
def fill_sb_analysis(output_csv, k, batch, stat_file):
    df = pd.read_csv(output_csv)
    fill_sb_analysis_from_df(output_csv, df, k, batch, stat_file)


def fill_sb_analysis_from_df(df_name, df, k, batch, stat_file):
    filename = df_name
    bias_mean = df['strand_bias_%'].mean().round(2)
    bias_median = df['strand_bias_%'].median().round(2)
    bias_modus = df['strand_bias_%'].mode().iloc[0]
    percentile_5 = round(df['strand_bias_%'].quantile(0.05), 2)
    percentile_95 = round(df['strand_bias_%'].quantile(0.95), 2)

    import csv
    stat = [filename, k, batch, bias_mean, bias_median, bias_modus, percentile_5, percentile_95]
    with open(stat_file, 'a', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(stat)

def get_binned_sb_info(df):
    pass#df.apply(lambda row: get_strand_bias_percentage(row["ratio"]), axis=1


def init_analysis(out_dir, filename):
    analysis = pd.DataFrame(
        data={},
        index=None,
        columns=['file', 'k', 'batch', 'bias_mean', 'bias_median', 'bias_modus', 'percentile_5', 'percentile_95'])
    analysis_name = utils.unique_path(os.path.join(out_dir, 'sb_analysis_' + filename + '.csv'))
    print("analysis stored in: {}".format(analysis_name))
    analysis.to_csv(analysis_name, index=False)
    return analysis_name


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


def analysis(analysis_full_path, output_dir, name, start_kmer=5, end_kmer=10, batch_count=49):
    #sb_analysis = init_analysis(analysis_full_path, name.split('.')[0])
    #df_analysis = pd.read_csv(sb_analysis)
    for kmer in range(start_kmer, end_kmer + 1):
        draw_conf_interval_graph(analysis_full_path, output_dir, utils.get_filename(name), kmer, batch_count)
        draw_basic_stats_lineplot(output_dir, utils.get_filename(name),
                                  analysis_full_path, kmer)


def draw_conf_interval_graph(dataframes, output_dir, name, k='', start_index=0):
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
        print('in ' + str(index))
         #df = df_batch.loc[df_batch["seq"].str.len() == k] #pd.read_csv(input_dir + 'df_output_' + str(k) + '_' + name + '_batch_' + str(batch) + '.csv')
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

    fig_name = utils.unique_path(os.path.join(output_dir, 'fig_ci_' + name + '_' + str(k) + '.png'))
    plt.savefig(fig_name)
    plt.close()



def draw_basic_stats_lineplot(output_dir, name, statfile, k, x_axis='k'):
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
    fig_name = utils.unique_path(os.path.join(output_dir, 'fig_lineplot_' + name + '_' + str(k) + '.png'))
    plt.savefig(fig_name)
    plt.close()
