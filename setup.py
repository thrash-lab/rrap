from setuptools import setup

setup(name='rrap',
      version='1.3.2',
      description='A metagenomic read recruitment data pipeline',
      url='https://github.com/Kojiconner/rrap_metag_pkg/',
      author='Conner Kojima',
      author_email='cykojima@usc.edu',
      license='MIT',
      packages=['rrap'],
      scripts=['bin/rrap', 'bin/rrap_test'],
      install_requires=[
        "joblib==0.16.0",
        "scikit-learn==0.23.1",
        "scipy>=1.5.2",
        "threadpoolctl==2.1.0"
      ],
      python_requires=">=3.7.3",
      package_dir={"": "src"},
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
      ],
      zip_safe=False,
      data_files=[("tests/reference", ["tests/reference/HTCC1062.fasta", "tests/reference/HIMB59.fasta"]),
                  ("tests/metaGs/dir_01", ["tests/metaGs/dir_01/ERR864073_toy_R1.fastq", "tests/metaGs/dir_01/ERR864073_toy_R2.fastq",
                                           "tests/metaGs/dir_01/ERR864077_toy_R1.fastq", "tests/metaGs/dir_01/ERR864077_toy_R2.fastq"]),
                  ("tests/metaGs/dir_02", ["tests/metaGs/dir_02/SRR11803378_toy_R1.fastq", "tests/metaGs/dir_02/SRR11803378_toy_R2.fastq"]),
                  ("tests/output", ["tests/output/allgenomes_cat_val.fna", "tests/output/total_reads_val.csv"]),
                  ("tests/output/bam/val", ["tests/output/bam/val/ERR864077.bam", "tests/output/bam/val/ERR864077.bam.bai",
                                             "tests/output/bam/val/SRR11803378.bam", "tests/output/bam/val/SRR11803378.bam.bai",
                                             "tests/output/bam/val/ERR864073.bam", "tests/output/bam/val/ERR864073.bam.bai"]),
                  ("tests/output/index/val", ["tests/output/index/val/val.1.bt2", "tests/output/index/val/val.2.bt2",
                                               "tests/output/index/val/val.3.bt2", "tests/output/index/val/val.4.bt2",
                                               "tests/output/index/val/val.rev.1.bt2", "tests/output/index/val/val.rev.2.bt2"]),
                  ("tests/output/rpkm/val", ["tests/output/rpkm/val/val_rpkm_log10.csv", "tests/output/rpkm/val/val_rpkm_noLog.csv"]),
                  ("tests/output/stats/val", ["tests/output/stats/val/ERR864077.bam.stats", "tests/output/stats/val/SRR11803378.bam.stats",
                                               "tests/output/stats/val/ERR864073.bam.stats"])])
