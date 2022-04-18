import os
import Bio.SeqIO as SeqIO
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


# CG content

def select_cg_range(value):
    for border in range(10, 100 + 1, 10):
        if value <= border:
            return border


def select_label(value):
    if value == 0 or value == 10:
        return "0-10% CG"
    return "{}-{}% CG".format(str(value - 9), str(value))


def plot_cg_content(df, output_name, order):
    ranged_cg_df = df.apply(lambda row: select_cg_range(row["CG_%"]), axis=1)
    c = Counter(ranged_cg_df)
    print(c.keys())
    labels = [select_label(key) for key in c.keys()]  # str(j*100/i) + "% CG" for j in range(1,i+1)

    print(labels)
    data = list(c.values())
    print(data)
    # Creating plot
    plt.figure(figsize=(10, 7))
    plt.pie(data, labels=labels)
    plt.savefig(output_name + ".png")


def plot_cg_from_dataframe(df_name):
    index = 0
    df = pd.read_csv(df_name)
    plt.figure(figsize=(10, 10))
    fig, axs = plt.subplots(3, 3)
    plot_cg_content(df, "cg_whole_file_" + os.path.basename(df_name), 0)

    for percent in [1, 3, 5, 10]:
        index = index + 1
        df_head = df.head(int(len(df) * (percent / 100)))
        df_tail = df.tail(int(len(df) * (percent / 100)))
        plot_cg_content(df_head, "cg_" + str(percent) + "%_highest_bias_" + os.path.basename(df_name), index + 1)
        plot_cg_content(df_tail, "cg_" + str(percent) + "%_lowest_bias_" + os.path.basename(df_name), index + 2)
        index = index + 2


# Nanopore bias comparation per hours

def bin_nanopore(fastq):
    file_type = "fastq" if fastq.split(".")[-1] == "fastq" else "fasta"

    batchfiles = []
    start = utc.localize(datetime.datetime.now())
    end = utc.localize(datetime.datetime(1970, 1, 1, 0, 0, 0))
    for record in SeqIO.parse(fastq, file_type):
        record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
        if record_time < start:
            start = record_time
        if record_time > end:
            end = record_time
    batches = hours_aligned(start, end)

    for record in SeqIO.parse(fastq, file_type):
        record_time = dparse([i for i in record.description.split() if i.startswith('start_time')][0].split('=')[1])
        batch = bisect.bisect_left(batches, record_time)
        filename = "nanopore_{0}_batch_{1}.fasta".format(get_filename(fastq), batch, file_type)
        if filename not in batchfiles:
            batchfiles.append(filename)
        f = open(filename, 'a')
        f.write(record.format("fasta"))
        f.close()


def filter_time(descr, tfrom, tto):
    # extract start time from description
    time = dparse([i for i in descr.split() if i.startswith('start_time')][0].split('=')[1])
    if dparse(tfrom) <= time <= dparse(tto):
        return True
    else:
        return False


def hours_aligned(start, end):
    chunks = []
    rule = rrule.rrule(rrule.HOURLY, byminute=0, bysecond=0, dtstart=start)
    for x in rule.between(start, end, inc=False):
        chunks.append(x)
    chunks.append(end)
    return chunks


def fill_sb_analysis(df_output, batch, stat_file):
    df = pd.read_csv(df_output)
    filename = df_output
    bias_mean = df["strand_bias_%"].mean().round(2)
    bias_median = df["strand_bias_%"].median().round(2)
    bias_modus = df["strand_bias_%"].round(2).mode().iloc[0]
    bias_5_top = df.head(int(len(df) * (5 / 100)))["strand_bias_%"].median().round(2)
    bias_5_bottom = df.tail(int(len(df) * (5 / 100)))["strand_bias_%"].median().round(2)

    import csv
    stat = [filename, batch, bias_mean, bias_median, bias_modus, bias_5_top, bias_5_bottom]
    with open(stat_file, 'a') as f:
        writer = csv.writer(f)
        writer.writerow(stat)


def get_filename(name):
    return os.path.basename(name).split(".")[0]


def init_analysis(out_dir, filename):
    analysis = pd.DataFrame(
        data={},
        index=None,
        columns=['file', 'batch', 'bias_mean', 'bias_median', 'bias_modus', 'bias_top_5', 'bias_bottom_5'])
    analysis.to_csv(out_dir + "sb_analysis_" + filename + ".csv", index=False)


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


def analysis(name, dir, start_kmer=5, end_kmer=10):
    init_analysis(dir, name.split(".")[0])
    for i in range(0, 49):
        kmer = 10
        fill_sb_analysis("df_output_" + str(kmer) + "_" + name.split("."[0]) + "_batch_" + str(i) + ".csv", i,
                         "sb_analysis_" + str(name.split("."[0])) + ".csv")
    for k in range(start_kmer, end_kmer + 1):
        draw_conf_interval_graph(name)


def draw_conf_interval_graph(name, k, batch_count):
    plt.xticks([x for x in range(49)], [str(x) for x in range(49)])
    plt.title('Confidence Interval')

    for batch in range(batch_count):
        print("in " + str(batch))
        df = pd.read_csv("df_output_" + str(k) + "_" + name + "_batch_" + str(batch) + ".csv")
        plot_confidence_interval(batch, df["strand_bias_%"])
    plt.savefig("fig_ci_" + name + "_" + str(k) + ".png")
    plt.show()


def draw_basic_stats_lineplot(name, k):
    df = pd.read_csv("strand_bias_analysis_6_2022-04-18.csv")

    # Plot a simple line chart
    plt.plot(df['batch'], df['bias_mean'], color='b', label='Mean value of strand bias')
    plt.plot(df['batch'], df['bias_median'], color='g', label='Median value of strand bias')
    plt.plot(df['batch'], df['bias_modus'], color='r', label='Mode value of strand bias')

    plt.legend()
    plt.savefig("fig_lineplot_" + name + "_" + str(k) + ".png")
    plt.show()
