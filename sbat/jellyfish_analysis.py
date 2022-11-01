import os
import subprocess
import sys

from sbat.utils import is_or_create_dir


class Jellyfish:
    """
    Class for representation of Jellyfish.
    """
    def __init__(self, threads=1, hash_size="100M"):
        self.threads = threads
        self.hash_size = hash_size
        self.jf_dir = None

    def set_outdir(self, output_dir):
        self.jf_dir = os.path.join(output_dir, 'jellyfish')
        is_or_create_dir(self.jf_dir)

    def run_jellyfish(self, input_file, k):
        """
        Method to execute jellyfish with correct arguments.
        :param input_file: input file
        :param k: size of k for which should jellyfish detect
        :return: file containing computed k-mers
        """
        if input_file is None:
            return None
        dump_file = os.path.join(self.jf_dir, "mer_counts.jf")

        try:
            calculate = "jellyfish count -m {0} -s {1} -t {2} -o {3} {4}".format(k, self.hash_size, self.threads,
                                                                                 dump_file, input_file)
            print(calculate)
            subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
            dump = "jellyfish dump {0}".format(dump_file)
            output_file = os.path.join(self.jf_dir, "output_" + str(k) + "_" + os.path.basename(input_file))
            with open(output_file, "w") as outfile:
                subprocess.run(dump.split(" "), stdout=outfile)
            return output_file
        except Exception as e:
            print(e)
            sys.exit("problem running jellyfish (possibly not installed), exiting...")
