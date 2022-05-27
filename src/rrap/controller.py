"""
 Copyright 2020, Conner Kojima <cykojima@usc.edu>

 This file is the main driver of rrap, a bionformatics pipeline that
 runs read recruitment and rpkm visualization given reference 
 genome and metagenome directories.

"""

import os
import argparse
import subprocess
from pathlib import Path
from rrap import concatenator
from rrap import indexer
from rrap import read_recruiter
from rrap import visualizer


def main():
    c = Controller()
    c.run()


class Controller:
    def __init__(self):
        self.p = argparse.ArgumentParser(prog="RRAP",
                                         description="Run read recruitment on a set of cleaned fna files")

        # argument groups
        self.inputs = None
        self.outputs = None
        self.optional = None
        self.subcommands = None
        self.arg_groups = []

        # args
        self.args = None

        # pipes
        self.concatenator = None
        self.indexer = None
        self.read_recruiter = None
        self.visualizer = None

        # intermediate products
        self.cat_file_path = None
        self.index_dir_path = None

        self.rpkm_heater_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rpkm_heater.py")
        self.stats_dir_path = None

    def run(self):
        self.add_arguments()

        print("---------making output dir if needed-------------")
        self.set_output_dir()
        self.cat_file_name = os.path.join(self.args.o, 'allgenomes_cat_{0}.fna'.format(self.args.n))

        if not self.args.index_pass:
            if self.args.crg:
                self.cat_file_path = self.args.crg
            else:
                print("---------concatenating reference genomes-------------")
                self.concatenator = concatenator.Concatenator(self.args)
                self.cat_file_path = self.concatenator.concatenate()

            print("---------indexing reference genomes-------------")
            self.indexer = indexer.Indexer(self.args, self.index_dir_path, self.cat_file_path)
            self.indexer.index()
        else:
            print("---------skipped indexing reference genomes-------------")

        if not self.args.rr_pass:
            print("---------read recruitment and data transform-------------")
            self.read_recruiter = read_recruiter.ReadRecruiter(self.args, self.index_dir_path,
                                                               self.cat_file_path, self.stats_dir_path,
                                                               self.bam_dir_path)
            self.read_recruiter.read_recruit()
        else:
            print("---------skipped read recruitment and data transform-------------")

        if not self.args.vis_pass:
            print("---------visualization-------------")
            self.visualizer = visualizer.Visualizer(self.args, self.rpkm_heater_path, self.stats_dir_path)
            self.visualizer.calculate_rpkm()
            if self.args.extra_vis:
                self.visualizer.plot_heatmaps()
        else:
            print("---------skipped visualization-------------")


    def add_arguments(self):
        # TODO specify argument groups
        self.inputs = self.p.add_argument_group("## input arguments")
        self.outputs = self.p.add_argument_group("## output arguments")
        self.options = self.p.add_argument_group("## options")

        self.arg_groups.extend([self.inputs, self.outputs, self.optional])

        # Add the arguments
        self.inputs.add_argument('-i', help='text file of all dir paths that contain cleaned metaG fna files', 
                                 required=True)
        self.inputs.add_argument('-crg', help=' path for concatenated reference genome fa file', 
                                 required=False)
        self.inputs.add_argument('-rg', help='input directory for reference genomes', 
                                 required=True)
        self.outputs.add_argument('-o', help='output directory path', 
                                  required=True)
        self.inputs.add_argument('-n', help='name of the project', required=True, metavar='project name')

        self.inputs.add_argument('-sort_gen', help='txt file of sorted genomes (if --extra-vis flag is used)', required=False)
        self.inputs.add_argument('-sort_samples', help='txt file of sorted samples (if --extra-vis flag is used)', required=False)
        self.inputs.add_argument("--threads", help='number of available threads', required=False)
        self.inputs.add_argument("-suffix", default="_pass_1.fastq", 
                                  help="everything in metaG file name that is after the acc for the forward (R1) read files \n"
                                  "e.g. (-QUALITY_PASSED_R1.fastq for <sample_acc>-QUALITY_PASSED_R1.fastq) \n" 
                                  "Otherwise, RRAP assumes that the forward pass file name is formatted as <acc>_pass_1.fastq"
                                  "NOTE: suffixes that contain a dash must specify '--' as an escape character e.g. '-suffix \"-- -QUALITY_PASSED_R1.fastq\"'")

        # specify optional args
        self.options.add_argument("--merge-contigs", default=False, dest='contig_merge',
                                  action='store_true', help="Concatenate contigs under individual organisms")
        self.options.add_argument("--skip-indexing", default=False, dest='index_pass',
                                  action='store_true', 
                                  help='Specify if the indexing step has already been completed and can be skipped. \
                                        If this flag is used, please check that bowtie2 index files exist e.g. \
                                        <output_dir_path>/index_dir/<project_name>.x.bt2')
        self.options.add_argument("--skip-rr", default=False, dest='rr_pass',
                                  action='store_true',
                                  help='Specify if the read recruitment step has already been completed and can be' \
                                       'skipped. If this flag is used, please check that read recruitment files exist' \
                                       'e.g. <output_dir_path>/')
        self.options.add_argument("--skip-vis", default=False, dest='vis_pass',
                                  action='store_true',
                                  help='Specify if the visualization step can be skipped')
        self.options.add_argument("--extra-vis", default=False, dest='extra_vis',
                                  action='store_true', 
                                  help='create csv with normalized RPKM values (log10) and plot heatmap using normalized values')

        self.args = self.p.parse_args()
        

    def set_output_dir(self):
        if self.args.o:
            # create output dir and create inner stats and index dir
            self.stats_dir_path = os.path.join(self.args.o, "stats", self.args.n)
            self.index_dir_path = os.path.join(self.args.o, "index", self.args.n)
            self.bam_dir_path = os.path.join(self.args.o, "bam", self.args.n)

            # make output (if it does not exist)
            if not os.path.isdir(self.args.o):
                subprocess.run("mkdir " + self.args.o, shell=True)

            # make stats, index, and bam dir
            if not os.path.isdir(os.path.join(self.args.o, "stats")):
                subprocess.run("mkdir " + os.path.join(self.args.o, "stats"), shell=True)
            if not os.path.isdir(os.path.join(self.args.o, "index")):
                subprocess.run("mkdir " + os.path.join(self.args.o, "index"), shell=True)
            if not os.path.isdir(os.path.join(self.args.o, "bam")):
                subprocess.run("mkdir " + os.path.join(self.args.o, "bam"), shell=True)

            # make project name specific dirs
            if not os.path.isdir(self.stats_dir_path):
                subprocess.run("mkdir " + self.stats_dir_path, shell=True)
            if not os.path.isdir(self.index_dir_path):
                subprocess.run("mkdir " + self.index_dir_path, shell=True)
            if not os.path.isdir(self.bam_dir_path):
                subprocess.run("mkdir " + self.bam_dir_path, shell=True)


if __name__ == "__main__":
    main()
