import datetime
import os
import shutil
import sys
import argparse

from sbat.analysis import Analysis
from sbat.jellyfish import Jellyfish
from sbat.nanopore import Nanopore
from sbat.interactive_plots import IPlotter, analyze_bias, run_server

from sbat.utils import check_if_nanopore, parse_iso_size
import re as re


__version__ = '0.0.7'
RUN_SERVER = False


class ParseSBAT(argparse.ArgumentParser):
    """
    Extension of common ArgumentParser class in order to allow pass only one or two arguments to 'mer' parameter.
    """
    def _match_argument(self, action, arg_strings_pattern):
        if action.dest == 'mer':
            narg_pattern = '(-*A{1,2})'
            match = re.match(narg_pattern, arg_strings_pattern)
            if match:
                return len(match.group(1))
            else:
                raise argparse.ArgumentError(action, "expected {} or {} arguments".format(1, 2))
        else:
            return super()._match_argument(action, arg_strings_pattern)


def arg_parser():
    """
    Function for parsing arguments of the application.
    """
    parser = ParseSBAT()
    parser.add_argument('input',
                        nargs='*',
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
                        nargs=2,
                        default=[0],
                        type=int,
                        metavar=("START_K", "[END_K]"),
                        help='k-mer size to count and analyze bias for. When only START_K is set, sbat computes '
                             'for only this k. If also END_K is set, range START_K-END_K is used. Default is range '
                             '5-10. MER must be >= 3 for analysis')
    parser.add_argument('-s', '--size',
                        nargs=1,
                        default=["100M"],
                        help='size of hash table for jellyfish, default 100M')
    parser.add_argument('-t', '--threads',
                        nargs=1,
                        type=int,
                        default=[1],
                        help='number of threads jellyfish shall use for computations')
    parser.add_argument('-r', '--subsample-reads',
                        nargs=1,
                        help='select number of reads you want to use in analysis per one bin, default all')
    parser.add_argument('-b', '--subsample-bases',
                        nargs=1,
                        help='select number of nucleotides you want to use in analysis per one bin, default all')
    parser.add_argument('-i', '--bin-interval',
                        nargs=1,
                        type=int,
                        default=[1],
                        help='number of hours that would fall into one bin when analysing nanopore')
    parser.add_argument('-g', '--margin',
                        nargs=1,
                        type=int,
                        default=[5],
                        help='number of %% taken into account when comparing highest N%% and lowest N%% of SB levels')
    parser.add_argument('-c', '--keep-computations',
                        action="store_true",
                        default=False,
                        help='keep jellyfish outputs and computations produced as partial results')
    parser.add_argument('-n', '--detect-nanopore',
                        action="store_true",
                        default=False,
                        help='identify nanopore datasets among inputs and provide advanced analysis')
    parser.add_argument('-p', '--interactive-plots',
                        action="store_true",
                        default=False,
                        help='run bokeh server with interactive version of plots')
    args = parser.parse_args()
    analysis_args = args_checker(args)
    return analysis_args


def args_checker(args):
    """
    Function for validating arguments of the application.
    :param args: output of ArgumentParser class containing all the arguments passed to the application
    :return: triplet of Analysis, Jellyfish and Nanopore objects initialized with app's input
    """
    jf = None
    nano = None
    a_args = Analysis()

    if args.version:
        version()
        sys.exit(0)

    if args.input is None or args.input == []:
        sys.exit("error: the following argument is required: input")

    for input_file in args.input:
        if not os.path.isfile(input_file):
            sys.exit("no such file: {}".format(input_file))

    if not os.path.isdir(args.output[0]):
        os.mkdir(args.output[0])

    a_args.set_output(args.output[0])
    a_args.keep_computations = args.keep_computations
    if args.margin[0] <= 0 or args.margin[0] > 100:
        sys.exit("margin must be number from interval (0, 100]")
    else:
        a_args.margin = args.margin[0]
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
        jf.hash_size = parse_iso_size(args.size[0])

    if len(args.mer) == 1:
        if args.mer[0] == 0:
            # MER not set, default boundaries to iterate upon
            a_args.start_k = 5
            a_args.end_k = 10
        elif args.mer[0] < 3:
            sys.exit("MER must be a positive number higher or equal to 3")
        else:
            # set only first boundary, run only for this k
            a_args.start_k = args.mer[0]
            a_args.end_k = args.mer[0]
    else:
        # run with specified boundaries
        if args.mer[1] < args.mer[0]:
            sys.exit("END_K must be bigger or equal to START_K")
        a_args.start_k = args.mer[0]
        a_args.end_k = args.mer[1]

    if args.interactive_plots:
        a_args.interactive = True
        global RUN_SERVER
        RUN_SERVER = True

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


def version():
    """
    Prints current version of the tool.
    """
    print("StrandBiasAnalysisTool v" + __version__)


if __name__ == '__main__':
    analysis, jf, nano, input_files = arg_parser()
    analyze_bias(analysis, jf, nano, input_files)
    if RUN_SERVER:
        run_server()
