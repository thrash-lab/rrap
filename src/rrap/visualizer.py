import subprocess
import os
import pandas as pd
import numpy as np
import sys


class Visualizer:
    def __init__(self, args, rpkm_heater_path, stats_dir_path):
        self.args = args
        self.rpkm_heater_path = rpkm_heater_path
        self.stats_dir_path = stats_dir_path

    def visualize(self):
        self.plot_heatmaps()


    def calculate_rpkm(self):
        if self.args.verbosity:
            print("calculating RPKM values with and without log normalization:")
        # make appropriate dir
        rpkm_output_dir = self.create_rpkm_output_dir()

        # df holds rpkm values with genome acc as the row names and metaG acc as the headers
        df = pd.DataFrame()
        tot_reads_df = pd.read_csv(os.path.join(self.args.o, "total_reads_{0}.csv".format(self.args.n)))

        # loop through files with bam.stats suffix in generated stats dir
        for file in os.listdir(self.stats_dir_path):
            if file.endswith(".bam.stats"):

                # rpkm dict will hold rpkm values for single metaG
                rpkm = {}
                # read csv
                entry = pd.read_csv(os.path.join(self.stats_dir_path, file), delimiter="\t", header=None)
                # add headers to pd
                entry.columns = ["genome", "gen_length", "r_mapped", "r_unmapped"]
                # calculate total mapped reads
                tot_reads = entry['r_mapped'].sum()

                # specify metaG accession in dict e.g. <acc>.bam.stats --> <acc>
                rpkm['ACC'] = file[:-10]

                # find tot_sample_reads
                total_sample_reads = int(tot_reads_df.loc[tot_reads_df["ACC"] == rpkm['ACC']]['total reads'])

                # specify rpkm values for each genome for this specific metaG
                for i in range(len(entry['genome'])):
                    if entry['genome'][i] != "*":
                        rpkm[entry['genome'][i]] = \
                            [entry['r_mapped'][i]/((entry['gen_length'][i]/1000)*(total_sample_reads/1000000))]

                # convert dict to data frame and transpose
                rpkm_df = pd.DataFrame(rpkm)
                rpkm_df.set_index('ACC', inplace=True)
                rpkm_df = rpkm_df.transpose()

                # add rpkm data to overall dataframe
                df = pd.concat([df, rpkm_df], axis=1, join='outer')
        
        # make sure column and index orderings are consistent across runs
        df = df.sort_values("ACC", axis=1)
        df = df.sort_index()

        # ignore divide by 0 warnings to by mumpy
        if not sys.warnoptions:
            import warnings
            warnings.simplefilter("ignore")

        # create log file
        df_log10 = np.log10(df)

        if self.args.verbosity:
            print("RPKM table")
            print(df, "\n")
            print("RPKM table log10 normalized")
            print(df_log10, "\n")

        df.to_csv(os.path.join(rpkm_output_dir, self.args.n + "_rpkm_noLog.csv"), index_label='ACC')
        df_log10.to_csv(os.path.join(rpkm_output_dir, self.args.n + "_rpkm_log10.csv"), index_label='ACC')
        

    def create_rpkm_output_dir(self):
        # make appropriate dir
        rpkm_base_dir = os.path.join(self.args.o, "rpkm")
        rpkm_output_dir = os.path.join(rpkm_base_dir, self.args.n)

        if not os.path.isdir(rpkm_base_dir):
            subprocess.run("mkdir {}".format(rpkm_base_dir), shell=True)
        if not os.path.isdir(rpkm_output_dir):  
            subprocess.run("mkdir {}".format(rpkm_output_dir), shell=True)
        
        return rpkm_output_dir

