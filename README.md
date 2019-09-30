[![Master branch build status](https://api.travis-ci.org/phac-nml/ecoli_serotyping.svg?branch=master "Master Build Status")](https://travis-ci.org/phac-nml/ecoli_serotyping)

# ECTyper (an easy typer)
**ectyper** is a standalone serotyping module for _Escherichia coli_. It supports _fasta_ and _fastq_ file formats.

# Dependencies:
- python >=3.5
- pytest >=3.5
- pandas 0.23.1
- samtools 1.8
- bowtie2 2.3.4.1
- mash 2.0
- bcftools 1.8
- biopython 1.70
- blast 2.7.1
- seqtk 1.2
- requests >=2.0

# Installation
1. Get `miniconda` if you do not already have `miniconda` or `anaconda`:
    1. `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh`
    1. `bash miniconda.sh -b -p $HOME/miniconda`
    1. `echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc`
    1.  `source ~/.bashrc`
2. Install ectyper  
    * Directly via `conda` 
    	1. `conda install -c bioconda ectyper` 
    * Through `github`
    	1. Install dependencies
          `pandas samtools bowtie2 mash bcftools biopython pytest blast seqtk tqdm python=3.6`
    	1. Download git repository then unzip
          `wget https://github.com/phac-nml/ecoli_serotyping/archive/master.zip`
    	1. Install ectyper inside unzipped directory
          `python setup.py install`

# Basic Usage
1. Put the fasta/fastq files for serotyping analyses in one folder (concatenate paired files if you would like them to be considered a single entity)
1. `ectyper -i [file path]`
1. View the results on the console or in `ectyper_[datetime]/output.csv`

# Example Usage
* `ectyper -i ecoliA.fasta`  for a single file
* `ectyper -i ecoliA.fasta -o output_dir` for a single file, results stored in `output_dir`
* `ectyper -i ecoliA.fasta,ecoliB.fastq,ecoliC.fna`	for multiple files  
* `ectyper -i ecoli_folder`	for a folder

# Advanced Usage
```
usage: ectyper [-h] [-V] -i INPUT [-c CORES] [-d PERCENTIDENTITY]
               [-l PERCENTLENGTH] [--verify] [-o OUTPUT] [-r REFSEQ] [-s]

ectyper v0.8.1 Prediction of Escherichia coli serotype from raw reads or
assembled genome sequences

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        Location of E. coli genome file(s). Can be a single
                        file, a comma-separated list of files, or a directory
  -c CORES, --cores CORES
                        The number of cores to run ectyper with
  -d PERCENTIDENTITY, --percentIdentity PERCENTIDENTITY
                        Percent identity required for an allele match [default
                        90]
  -l PERCENTLENGTH, --percentLength PERCENTLENGTH
                        Percent length required for an allele match [default
                        50]
  --verify              Enable E. coli species verification
  -o OUTPUT, --output OUTPUT
                        Directory location of output files
  -r REFSEQ, --refseq REFSEQ
                        Location of pre-computed MASH RefSeq sketch. If
                        provided, genomes identified as non-E. coli will have
                        their species identified using MASH. For best results
                        the pre-sketched RefSeq archive https://gembox.cbcb.um
                        d.edu/mash/refseq.genomes.k21s1000.msh is recommended
  -s, --sequence        Prints the allele sequences if enabled as the final
                        columns of the output
```

