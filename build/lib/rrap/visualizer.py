import subprocess
import os
from pathlib import Path


class Visualizer:
    def __init__(self, args, rpkm_heater_path, stats_dir_path):
        self.args = args
        self.rpkm_heater_path = rpkm_heater_path
        self.stats_dir_path = stats_dir_path

    def visualize(self):
        self.plot_heatmaps()
        self.calculate_coverage_and_depth()

    def plot_heatmaps(self):
        # TODO
        rpkm_output_dir = os.path.join(self.args.o, "rpkm")
        subprocess.run("mkdir {}".format(rpkm_output_dir), shell=True)
        sort_gen_addon = ""
        sort_samples_addon = ""

        # incorporate optional sort_gen and sort_sample options for rpkm_heater
        if self.args.sort_gen:
            sort_gen_addon = "-sort_gen {0}".format(self.args.sort_gen)
        if self.args.sort_samples:
            sort_samples_addon = "-sort_samples {0}".format(self.args.sort_samples)

        cmd = "rpkm_heater -map -i {1} -o {2} -project {3} {4} {5}".format(self.rpkm_heater_path,
                                                                           self.stats_dir_path,
                                                                           rpkm_output_dir,
                                                                           self.args.n,
                                                                           sort_gen_addon,
                                                                           sort_samples_addon)

        print("running: " + cmd)
        subprocess.run(cmd, shell=True)

    def calculate_coverage_and_depth(self):
        # TODO
        pass
