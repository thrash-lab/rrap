## Overview
*RRAP* is a bionformatics pipeline that aligns reads from metagenomes to reference
genomes. *RRAP* uses the following programs: bowtie2 and SAMTools. It
essentially acts a wrapper that handles the logistics of read recruitment and walks
the user through the process of concatenating genomes, creating an index, aligning
reads, and extracting RPKM values from BAM files. 

## Obtaining RRAP
*RRAP* is available as a conda package. Create and activate a new conda environment 
with the following commands: `conda create -n rrap_env` and `conda activate rrap_env`.
Next, RRAP can be installed with `conda install -c kojiconner -c bioconda rrap`.
The RRAP environment can be exited with the `conda deactivate` command. Just remember
to use `conda activate rrap_env` to startup RRAP again.

## Usage
The most common way to use *RRAP* is with some variation of the following command:
`rrap -i  <path_to_txt_file_that_contains_metagenome_dir_paths> -rg 
<reference_genome_dir_path> -o <output_dir_path> -n <project_name> 
--heat-map -suffix <suffix>`. The txt file should contain a metagenome dir path
on each line. 

## Quickstart
The quickest way to get started with *RRAP* is to run the command `rrap_test`
through the terminal. The terminal will automatically run the command: 
`rrap -i  <path_to_test_data_dir>/metaG_paths.txt -rg <path_to_test_data_dir>/reference
 -o <path_to_test_data_dir>/output -n test --merge-contigs --extra-vis -suffix 
_toy_R1.fastq` which will walk the user through processing a sample dataset. In
addition, the program will run some tests to ensure that *RRAP* is running as
intended. The sample dataset contains the following reference genomes: HIMB59 and 
HTCC1062 (1, 2). The test command will align reads from the metagenomes with the 
following accessions: ERR864073, ERR864077, SRR11803378 (3, 4). After mapping reads, 
RRAP will report RPKM values and confirm that the program is working as intended.

NOTE: Because metagenomic files require a lot of storage, RRAP uses a randomized subset
of 10,000 reads from the metagenomes ERR864073 and SRR11803378. Likewise, 20,000 reads 
were subset from the metagenome ERR864077. Subsetting was performed with seqkit (5).

### References (for dataset)
1. Grote, J., Thrash, J. C., Huggett, M. J., Landry, Z. C., Carini, P., Giovannoni, S. J., & Rappé, M. S. (2012). Streamlining and core genome conservation among highly divergent members of the SAR11 clade. MBio, 3(5), e00252-12.
2. Giovannoni, S. J., Tripp, H. J., Givan, S., Podar, M., Vergin, K. L., Baptista, D., Bibbs, L., Eads, J., Richardson, T. H., Noordewier, M., Rappé, M. S., Short, J. M., Carrington, J. C., & Mathur, E. J. (2005). Genome streamlining in a cosmopolitan oceanic bacterium. Science (New York, N.Y.), 309(5738), 1242–1245. https://doi.org/10.1126/science.1114057
3. Fortunato, C. S., & Crump, B. C. (2015). Microbial gene abundance and expression patterns across a river to ocean salinity gradient. PLoS One, 10(11), e0140578.
4. Sakowski, E. G., Arora-Williams, K., Tian, F., Zayed, A. A., Zablocki, O., Sullivan, M. B., & Preheim, S. P. (2021). Interaction dynamics and virus–host range for estuarine actinophages captured by epicPCR. Nature Microbiology, 6(5), 630-642.
5. Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), e0163962.

## Links to Other Repositories
* https://github.com/samtools/samtools
* https://github.com/BenLangmead/bowtie2
