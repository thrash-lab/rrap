import subprocess
import os


class Indexer:
    def __init__(self, args, output_path, cat_file_path):
        self.args = args
        self.output_path = output_path
        self.cat_file_path = cat_file_path

    def index(self):
        if self.args.rg:
            threads_addon = ""
            if self.args.threads:
                threads_addon = "--threads {}".format(self.args.threads)
        
            quiet_addon = ""
            if self.args.verbosity:
                print("Running the following command:")
                print("bowtie2-build {0} {1}\n".format(self.cat_file_path, self.output_path, threads_addon))
            else:
                quiet_addon = "--quiet"

            subprocess.run("bowtie2-build {3} {2} {0} {1}".format(self.cat_file_path, self.output_path, 
                                                                  threads_addon, quiet_addon),
                           shell=True)
        else:
            raise Exception("reference genome dir not specified")
