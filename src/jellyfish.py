import os
import subprocess

import utils


class Jellyfish:
    def __init__(self, threads=1, hash_size="500M", output_dir=''):
        self.threads = threads
        self.hash_size = hash_size
        self.jf_dir = self.set_outdir(output_dir)


    def set_outdir(self, output_dir):
        self.jf_dir = os.path.join(output_dir, 'sbat', 'jellyfish')
        utils.is_or_create_dir(self.jf_dir)

    def run_jellyfish(self, input_file, k):
        dump_file = os.path.join(self.jf_dir, "mer_counts.jf")
        calculate = "jellyfish count -m {0} -s {1} -t {2} -o {3} {4}".format(k, self.hash_size, self.threads,
                                                                             dump_file, input_file)
        print(calculate)
        subprocess.run(calculate.split(" "), stdout=subprocess.PIPE)
        dump = "jellyfish dump {0}".format(dump_file)
        output_file = os.path.join(self.jf_dir, "output_" + str(k) + "_" + os.path.basename(input_file))
        with open(output_file, "w") as outfile:
            subprocess.run(dump.split(" "), stdout=outfile)
        return output_file