import datetime
import os
import shutil
import sys
import argparse

from sbat.jellyfish import Jellyfish
from sbat.nanopore import Nanopore
from sbat.analysis import Analysis
from sbat.utils import check_if_nanopore, parse_iso_size
import re as re


__version__ = '0.1.0'


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
                        nargs='+',
                        type=str,
                        help='fasta or fastq file to count and analyze strand bias on')
    parser.add_argument('-v', '--version',
                        action="store_true",
                        default=False,
                        help='current version of the application'
                        )
    parser.add_argument('-j', '--no-jellyfish',
                        action='store_true',
                        default=False,
                        help='skip k-mer counting. Requires input in fasta file where id=count, seq=k-mer')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=["sbat_out/"],
                        help='output directory')
    parser.add_argument('-m', '--mer',
                        nargs=1,
                        default=[0],
                        type=int,
                        help='k-mer to count and analyze bias for. When set to 0, bias is analyzed for all k-mer '
                             'in 5-10 range. MER must be >= 3 for analysis')
    parser.add_argument('-s', '--size',
                        nargs=1,
                        default="100M",
                        help='size of hash table for jellyfish')
    parser.add_argument('-t', '--threads',
                        nargs=1,
                        type=int,
                        default=[1],
                        help='number of threads jellyfish shall use for computations')
    parser.add_argument('-r', '--subsample-reads',
                        nargs=1,
                        type=int, # TODO convert to format of hash?
                        help='select number of reads you want to use in analysis per one bin, default all')
    parser.add_argument('-b', '--subsample-bases',
                        nargs=1,
                        type=int, # TODO convert to format of hash?
                        help='select number of nucleotides you want to use in analysis per one bin, default all')
    parser.add_argument('-i', '--bin-interval',
                        nargs=1,
                        type=int,
                        default=[1],
                        help='number of hours that would fall into one bin when analysing nanopore')
    parser.add_argument('-w', '--whisker',
                        nargs=1,
                        type=int,
                        default=[5],
                        help='number of % taken into account when comparing highest N% and lowest N% of SB levels')
    parser.add_argument('-c', '--keep-computations',
                        action="store_true",
                        default=False,
                        help='keep jellyfish outputs and computations produced as partial results')
    parser.add_argument('-n', '--detect-nanopore',
                        action="store_true",
                        default=False,
                        help='identify nanopore datasets among inputs and provide advanced analysis')
    args = parser.parse_args()
    analysis_args = args_checker(args)
    return analysis_args


def args_checker(args):
    jf = None
    nano = None
    a_args = Analysis()

    if args.version:
        version()
        sys.exit(0)

    for input_file in args.input:
        if not os.path.isfile(input_file):
            sys.exit("no such file: {}".format(input_file))

    if not os.path.isdir(args.output[0]):
        os.mkdir(args.output[0])

    a_args.set_output(args.output[0])
    a_args.keep_computations = args.keep_computations
    if args.whisker[0] <= 0 or args.whisker[0] > 100:
        sys.exit("whisker must be number from interval (0, 100]")
    else:
        a_args.whisker = args.whisker[0]
    if args.mer[0] < 3 and args.mer[0] != 0:
        sys.exit("MER must be a positive number higher or equal to 3")

    if args.no_jellyfish and args.detect_nanopore:
        sys.exit("cannot detect nanopore when jellyfish off - nanopore requires jellyfish for analyses")

    if not args.no_jellyfish:
        jf = Jellyfish()
        jf.set_outdir(args.output[0])
        if args.threads[0] < 1:
            sys.exit("number of threads must be a positive integer")
        else:
            a_args.threads = args.threads[0]
            jf.threads = args.threads[0]
            print(args.size)
        jf.hash_size = parse_iso_size(args.size)

    if args.mer[0] == 0:
        # default boundaries to iterate upon
        a_args.start_k = 5
        a_args.end_k = 10
    else:
        # run with specified k
        a_args.start_k = args.mer[0]
        a_args.end_k = args.mer[0]

    if args.detect_nanopore:
        nano = Nanopore()
        nano.init_common(a_args, jf)

        if args.subsample_reads is not None:
            nano.subs_reads = parse_iso_size(args.subsample_reads[0])
        if args.subsample_bases is not None:
            nano.subs_bases = parse_iso_size(args.subsample_bases[0])
        if args.bin_interval is not None:
            nano.bin_interval = args.bin_interval[0]

    return a_args, jf, nano, args.input


def main():
    analysis, jf, nano, input_files = arg_parser()
    for file in input_files:
        analysis.set_file(file)
        analysis.init_analysis()
        print("input: " + file)
        dfs = []
        run_nano = nano is not None and check_if_nanopore(file)
        for k in range(analysis.start_k, analysis.end_k + 1):
            if jf is not None:
                print("running computation and analysis for K=" + str(k))
                jf_output = jf.run_jellyfish(file, k)

                if run_nano:
                    print('running nanopore analysis...')
                    nano.nanopore_analysis(file)
                    run_nano = False  # run analysis only once for each input file
            else:
                print("jellyfish disabled, running only analysis...")
                jf_output = file
            df = analysis.jellyfish_to_dataframe(jf_output, k)  # convert jellyfish results to DataFrame
            dfs.append(df)
            if df is not None:
                analysis.plot_kmers_vs_bias(df, k)
        analysis.draw_basic_stats_lineplot(analysis.filename, analysis.sb_analysis_file)
        analysis.plot_conf_interval_graph(dfs, start_index=analysis.start_k)
        analysis.plot_gc_from_dataframe(dfs)
    if not analysis.keep_computations:
        shutil.rmtree(analysis.dump_dir)
        if jf is not None:
            shutil.rmtree(jf.jf_dir)


def version():
    """
    Prints current version of the tool.
    """
    print("StrandBiasAnalysisTool v" + __version__)


if __name__ == '__main__':
    main()
