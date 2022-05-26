import subprocess
import os
from pathlib import Path


class ReadRecruiter:
    def __init__(self, args, index_dir_path, cat_file_path, stats_dir_path):
        self.args = args
        self.index_dir_path = index_dir_path
        self.cat_file_path = cat_file_path
        self.stats_dir_path = stats_dir_path

        if self.args.suffix[0:2] == "--":
            self.args.suffix = self.args.suffix[3:]

    def read_recruit(self):
        if self.args.i:

            print("*** get list of acc ***")
            # retrieve individual clean dir paths from -i flag
            with open(self.args.i) as file:
                clean_dir_paths = file.readlines()
                clean_dir_paths = [line.rstrip() for line in clean_dir_paths]

            for clean_dir_path in clean_dir_paths:
                # create list of tuples with each tuple containing the forward and reverse file name for a sample
                tuple_list = self.find_acc(clean_dir_path)

                print("Running read recruitment for dir: " + clean_dir_path)
                self.align_reads(tuple_list, clean_dir_path, self.index_dir_path)

                output_dir_path = os.path.join(clean_dir_path, "rrap_output_dir_" + self.args.n)
                print("copying .bam.stats files to stats dir: " + clean_dir_path)
                subprocess.run('cp {0} {1}'.format(os.path.join(output_dir_path, "*.bam.stats"),
                                                   self.stats_dir_path), shell=True)
                print("read recruitment complete: " + clean_dir_path)
        else:
            pass

    def find_acc(self, path):
        # get list of all fastq files (r1 and r2)
        # self.args.suffix.split(".")[-1] denotes the fastq extension e.g.fastq or fq
        cmd = subprocess.run('ls {0}'.format(os.path.join(path, "*" + self.args.suffix.split(".")[-1])),
                             shell=True, capture_output=True, text=True)
        acc_list = cmd.stdout.split()
        acc_list = [line.rstrip() for line in acc_list]

        # make sure that half the fastq files contain the forward pass (r1) suffix
        print(acc_list)
        print(self.args.suffix)

        contains_suffix = [sample for sample in acc_list if self.args.suffix in sample]
        if len(contains_suffix) == 0:
            raise IOError("Incorrect suffix specified for forward pass fastq files in path {0} Check the -suffix argument".format(path))
        elif len(contains_suffix) * 2 != len(acc_list):
            raise IOError("Incorrect suffix specified for forward pass fastq files in path {0}. Check that all fastq files have the same"
                          "suffix. Check the -suffix argument".format(path))

        # sort list alphanumerically
        acc_list.sort()

        # create list of tuples with each tuple containing an even/odd pair e.g. (0, 1), (2, 3), (4, 5)
        tuple_list = []
        for i in range(int(len(acc_list)/2)):
            # general format of file name is <acc><suffix>
            acc = acc_list[i*2].split(self.args.suffix)[0]

            # acc_list[i*2] should denote r1 file, acc_list[i*2] should denote r2 file 
            tuple_list.append((acc, acc_list[i*2], acc_list[(i*2)+1]))

        return tuple_list


    def align_reads(self, tuple_list, clean_dir_path, index_path):
        # create an output dir in each fastq dir specified with -i
        output_dir_path = os.path.join(clean_dir_path, "rrap_output_dir_" + self.args.n)
        subprocess.run("mkdir {0}".format(output_dir_path), shell=True)

        # for each acc
        for sample in tuple_list:
            print("sample: ", sample)
            # create helpful variables
            acc = sample[0].split("/")[-1]
            r1_path = sample[1]
            r2_path = sample[2]
            acc_path = sample[0]

            print("working on sample:", acc, "\n")

            # only run bowtie2 if .bam.stats file doesn't exist
            if not os.path.exists(os.path.expanduser(os.path.join(output_dir_path, "{0}.bam.stats".format(os.path.basename(acc))))):
                # run read recruitment
                print("previous file does not exist: running read recruitment")

                threads_addon = ""
                if self.args.threads:
                    threads_addon = "--threads {}".format(self.args.threads)

                command = 'bowtie2 {3} -x "{0}" -1 "{4}" -2 "{5}" ' \
                          '--no-unal -S "{2}.sam"'.format(index_path, acc_path, acc, threads_addon, r1_path, r2_path)
                subprocess.run(command, shell=True)

                self.generate_stats_file(acc, output_dir_path)

    def generate_stats_file(self, acc, output_dir_path):
        # convert output SAM fle into a BAM file
        subprocess.run('samtools view -F 4 -bS "{0}.sam" > "{0}-RAW.bam"'.format(acc), shell=True)

        # sort and index the BAM file
        subprocess.run('samtools sort "{0}-RAW.bam" -o "{0}.bam"'.format(acc), shell=True)
        subprocess.run('samtools index "{0}.bam"'.format(acc), shell=True)
        subprocess.run('samtools idxstats "{0}.bam" > "{0}.bam.stats"'.format(acc), shell=True)

        # copy bam.stats files into output dir
        subprocess.run('mv {0}.bam.stats {1}'.format(acc, output_dir_path), shell=True)

        # TODO remove temporary files
        subprocess.run('rm "{0}.sam" "{0}-RAW.bam"'.format(acc), shell=True)