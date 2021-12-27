import subprocess


class Indexer:
    def __init__(self, args, index_dir_path, cat_file_path):
        self.args = args
        self.index_dir_path = index_dir_path
        self.cat_file_path = cat_file_path

    def index(self):
        if self.args.rg:
            threads_addon = ""
            if self.args.threads:
                threads_addon = "--threads {}".format(self.args.threads)
            print("Attempting to run command:")
            print("bowtie2-build {0} {1}".format(self.cat_file_path, self.index_dir_path, threads_addon))
            subprocess.run("bowtie2-build {0} {1}".format(self.cat_file_path, self.index_dir_path, threads_addon),
                           shell=True)
        else:
            raise Exception("reference genome dir not specified")
