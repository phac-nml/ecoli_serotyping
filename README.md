# ECTyper (an easy typer)
**ECTyper** wraps a standalone serotyping module.  
Support _fasta_ and _fastq_ format

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
    	1. `conda install -c anaconda -c bioconda -c ubcsamsung ectyper`  
    * Through `github`  
    	1. Install dependencies
          `conda install pandas samtools bowtie2 mash bcftools biopython nose blast seqtk tqdm python=3.6`
    	1. Download git repository then unzip
          `wget https://github.com/phac-nml/ecoli_serotyping/archive/master.zip`
    	1. Install ectyper inside unzipped directory
          `python setup.py install`

# Basic Usage
1. Put all you fasta/fastq file in one folder (concatenate paired files if you want to result to be considered as single entity)
1. `ectyper -i [file path]`
1. View result on console or in `output/output.csv`
* If you want to enable species identification,
    1. Download refseq from https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
    2. Put refseq into ectyper/Data/
    3. Specify `-s` as argument when executing ectyper

# Example Usage
* `ectyper -i ecoliA.fasta`  for single file
* `ectyper -i ecoliA.fasta,ecoliB.fastq,ecoliC.fna`	for multiple file  
* `ectyper -i ecoli_folder`	for folder
