import subprocess
import os
from pathlib import Path

class ReadRecruiter:
    def __init__(self, args, index_dir_path, cat_file_path, stats_dir_path):
        print("got here")
        self.args = args
        self.index_dir_path = index_dir_path
        self.cat_file_path = cat_file_path
        self.stats_dir_path = stats_dir_path


    def read_recruit(self):
        if self.args.i:

            print("*** get list of acc ***")
            # retrieve individual clean dir paths from -i flag
            with open(self.args.i) as file:
                clean_dir_paths = file.readlines()
                clean_dir_paths = [line.rstrip() for line in clean_dir_paths]

            # create an intermediate acc_list for each clean_dir_path
            for clean_dir_path in clean_dir_paths:
                # get list of cleaned fastq files
                acc_list = self.list_files_with_suffix(clean_dir_path, 'QUALITY_PASSED_R2.fastq')

                print("Running read recruitment for dir: " + clean_dir_path)
                self.align_reads(acc_list, clean_dir_path, self.index_dir_path)

                output_dir_path = os.path.join(clean_dir_path, "rrap_output_dir_" + self.args.n)
                print("copying .bam.stats files to stats dir: " + clean_dir_path)
                subprocess.run('cp {0} {1}'.format(os.path.join(output_dir_path, "*.bam.stats"),
                                                   self.stats_dir_path), shell=True)
                print("read recruitment complete: " + clean_dir_path)
        else:
            pass

    def list_files_with_suffix(self, path, suffix):
        cmd = subprocess.run('ls {0}'.format(os.path.join(path, "*" + suffix)),
                             shell=True, capture_output=True, text=True)
        acc_list = cmd.stdout.split()
        acc_list = [line.rstrip() for line in acc_list]
        acc_list = [line[:-24] for line in acc_list]
        return acc_list

    def align_reads(self, acc_list, clean_dir_path, index_path):
        # access every sample acc in given acc list

        output_dir_path = os.path.join(clean_dir_path, "rrap_output_dir_" + self.args.n)
        subprocess.run("mkdir {0}".format(output_dir_path), shell=True)

        # for each acc
        for acc in acc_list:
            print("working on sample:", acc, "\n")

            # create helpful variables
            acc_path = os.path.join(output_dir_path, os.path.basename(acc))

            print("checking if file exists: ", os.path.expanduser(os.path.join(output_dir_path, "{0}.bam.stats".format(os.path.basename(acc)))))
            # only run bowtie2 if .bam.stats file doesn't exist
            if not os.path.exists(os.path.expanduser(os.path.join(output_dir_path, "{0}.bam.stats".format(os.path.basename(acc))))):
                # run read recruitment
                print("previous file does not exist: running read recruitment")

                threads_addon = ""
                if self.args.threads:
                    threads_addon = "--threads {}".format(self.args.threads)

                command = 'bowtie2 {3} -x "{0}" -1 "{2}-QUALITY_PASSED_R1.fastq" -2 ' \
                          '"{2}-QUALITY_PASSED_R2.fastq" ' \
                          '--no-unal -S "{2}.sam"'.format(index_path, acc_path, acc, threads_addon)
                print("aligning reads: ", command)
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
