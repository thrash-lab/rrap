## Overview
*RRAP* is a bionformatics pipeline that aligns reads from metagenomes to reference
genomes. *RRAP* uses the following programs: bowtie2, samtools, and rpkm_heater. It
essentially acts a wrapper that handles the logistics of read recruitment and walks
the user through the process of concatenating genomes, creating an index, aligning
reads, and extracting RPKM values from BAM files. 

## Obtaining RRAP

*RRAP* is available as a conda package. The recommended method of 
installation is with the following command `conda install -c kojiconner rrap`.

## Usage

The most common way to use *RRAP* is with some variation of the following command:
`rrap -i  <path_to_txt_file_that_contains_metagenome_path_dirs> -rg 
<path_to_reference_genome_dir> -o <path_to_output_dir> -n <project_name> 
--extra-vis -suffix <suffix>`.

## Quickstart
The quickest way to get started with *RRAP* is to run the command `rrap_test`
through the terminal. The terminal will automatically run the command: 
`rrap -i  <path_to_test_data_dir>/metaG_paths.txt -rg <path_to_test_data_dir>/reference
 -o <path_to_test_data_dir>/output -n test --merge-contigs --extra-vis -suffix 
_toy_R1.fastq` which will walk the user through processing a sample dataset. In
addition, the program will run some tests to ensure that *RRAP* is running as
intended. 
