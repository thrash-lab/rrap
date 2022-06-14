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
intended. More specifically, it will use the following reference genomes: HIMB59 and 
HTCC1062 (1, 2). The test command will align reads from the metagenomes with the 
following accessions: ERR864073, ERR864077, SRR11803378 (3, 4). After mapping reads, 
RRAP will report RPKM values and confirm that the program is working as intended.

NOTE: Because metagenomic files require a lot of storage, RRAP uses a randomized subset
of 5000 reads for each metagenome. The subsetting was performed with seqkit (5).

### References (for dataset)
1. Grote, J., Thrash, J. C., Huggett, M. J., Landry, Z. C., Carini, P., Giovannoni, S. J., & Rappé, M. S. (2012). Streamlining and core genome conservation among highly divergent members of the SAR11 clade. MBio, 3(5), e00252-12.
2. Rappé, M. S., Connon, S. A., Vergin, K. L., & Giovannoni, S. J. (2002). Cultivation of the ubiquitous SAR11 marine bacterioplankton clade. Nature, 418(6898), 630-633.
3. Fortunato, C. S., & Crump, B. C. (2015). Microbial gene abundance and expression patterns across a river to ocean salinity gradient. PLoS One, 10(11), e0140578.
4. Sakowski, E. G., Arora-Williams, K., Tian, F., Zayed, A. A., Zablocki, O., Sullivan, M. B., & Preheim, S. P. (2021). Interaction dynamics and virus–host range for estuarine actinophages captured by epicPCR. Nature Microbiology, 6(5), 630-642.
5. Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), e0163962.

## Links to Other Repositories
* https://github.com/samtools/samtools
* https://github.com/BenLangmead/bowtie2
* https://github.com/thrash-lab/rpkm_heater
