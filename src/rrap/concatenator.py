import subprocess
import os
from pathlib import Path

class Concatenator:
    def __init__(self, args):
        self.args = args
        self.cat_file_name = os.path.join(self.args.o, 'allgenomes_cat_{0}.fna'.format(self.args.n))

    def concatenate(self):
        # replace contigs with file name

        # delete cat file if it already exists
        if os.path.exists(self.cat_file_name):
            if self.args.verbosity:
                print("deleting preexisting concatenated file: {0}".format(self.cat_file_name))
            subprocess.run('rm -f {0}'.format(self.cat_file_name), shell=True)

        if self.args.verbosity:
            print("concatenating genomes from the dir {0}".format(self.args.rg))

        if self.args.contig_merge:
            # get list of all ref genome file paths
            cmd = subprocess.run('ls {0}'.format(os.path.join(self.args.rg, "*")),
                                    shell=True, capture_output=True, text=True)
            rf_list = cmd.stdout.split()

            # read into each rf genome and change headers to group contigs under organism header (file name)
            for rf_path in rf_list:
                self.combine_contigs(rf_path)
        else:
            # concatenate all ref genomes into one file
            subprocess.run("cat {0} > {1}".format(os.path.join(self.args.rg, "*"), self.cat_file_name), shell=True)
        return self.cat_file_name

    def combine_contigs(self, rf_path):
        # replace initial header line with line that contains just organism e.g. (>Candidatus Fonsibacter)
        # remove all subsequent header lines
        # this line remains in place to deal with files that contain many contigs
        rf_name = Path(rf_path).stem

        # read rf genome and add initial header/remove contig headers
        with open(rf_path) as file:
            lines = file.readlines()
            lines = [line for line in lines[1:] if not line.startswith(">")]
            lines.insert(0, ">" + rf_name + "\n")

        # write modified lines into concatenated fle
        with open(self.cat_file_name, 'a') as f:
            for line in lines:
                f.write("%s" % line)
