# ECTyper (an easy typer)
**ecyper** wraps a standalone serotyping module for _Escherichia coli_.
Supports _fasta_ and _fastq_ file formats.

[![Build Status](https://travis-ci.org/phac-nml/ecoli_serotyping.svg?branch=master)](https://travis-ci.org/phac-nml/ecoli_serotyping)

# Dependencies:
- python 3.6.3.*
- pandas 0.21.0.*
- samtools 1.5.*
- bowtie2 2.3.0.*
- mash 1.1.*
- bcftools 1.6.*
- biopython 1.69.*
- blast 2.2.31 .*
- seqtk 1.2.*

# Installation
1. Get `miniconda` if you do not already have `miniconda` or `anaconda`:
    1. `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh`
    1. `bash miniconda.sh -b -p $HOME/miniconda`
    1. `export PATH="$HOME/miniconda/bin:$PATH"`
2. Install ectyper  
    * Directly via `conda`
    	1. `conda install -c bioconda ectyper`
    * Through `github`
    	1. Install dependencies
          `conda install pandas samtools bowtie2 mash bcftools biopython nose blast seqtk tqdm python=3.6`
    	1. Download git repository then unzip
          `wget https://github.com/phac-nml/ecoli_serotyping/archive/master.zip`
    	1. Install ectyper inside unzipped directory
          `python setup.py install`

# Basic Usage
1. Put all of your fasta/fastq files in one folder (concatenate paired files if you want the result to be considered a single entity)
1. `ectyper -i [file path]`
1. View the results on the console or in `output/[datetime]/output.csv`

# Example Usage
* `ectyper -i ecoliA.fasta`  for a single file
* `ectyper -i ecoliA.fasta,ecoliB.fastq,ecoliC.fna`	for multiple files  
* `ectyper -i ecoli_folder`	for a folder

# Advanced Usage
```
usage: ectyper [-h] -i INPUT [-d PERCENTIDENTITY] [-l PERCENTLENGTH]
               [--verify] [-s] [--detailed] [-o OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Location of E. coli genome file(s). Can be a single
                        file or a directory
  -d PERCENTIDENTITY, --percentIdentity PERCENTIDENTITY
                        Percent identity required for an allele match [default 90]
  -l PERCENTLENGTH, --percentLength PERCENTLENGTH
                        Percent length required for an allele match [default 50]
  --verify              Enable E. coli species verification
  -s, --species         Enable species identification when a non-E. coli
                        genome is found Note: refseq downloading is required
                        when running this option for the first time
  --detailed            Enable detailed program output
  -o OUTPUT, --output OUTPUT
                        Directory location of output files.
```
* The first time species identification is enabled you will need to wait for **ectyper** to download the reference sequences.

# Building the conda package
Python 2.7
(requires a custom version of process32 from the channel kevinkle)

`conda build -c kevinkle -c bioconda recipe/ --python=2.7`

Python 3.6

`conda build -c bioconda recipe/ --python=3.6`
