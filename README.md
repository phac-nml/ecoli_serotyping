# ECTyper (an easy typer)
The `ectyper` wraps a standalone serotyping module.  
Support FASTA and FASTQ format

# Dependencies:
* python (v3.6.2)
* bowtie2 (v2.2.6)
* SAMtools (v1.6)
* bcftools (v1.2)
* mash (v2.0)

# Installation
1. Get `miniconda` if you do not already have `miniconda` or `anaconda`:
    1. `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh`
    1. `bash miniconda.sh -b -p $HOME/miniconda`
    1. `export PATH="$HOME/miniconda/bin:$PATH"`
1. Add `bioconda` to conda channel
    1. `conda config --add channels bioconda`
1. Install dependencies
    1. `conda install samtools bowtie2 mash bcftools biopython nose blast tqdm python=3.6`
1. Install ectyper
    1. `python setup.py install`

# Basic Usage
1. Put all you fasta/fastq file in one folder
1. `ectyper -i [dir_path or file_path]`
1. View result on console or in `output.json`
* If you want to enable species identification,
    1. Download refseq from https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
    2. Put refseq into ectyper/Data/
    3. Specify `-s` as argument when executing ectyper

# Benchmark
## Genbank Database (452 genomes/836 metadata)
* concordance to provided metadata: 776 (92.82%)
* discrepancies with provided metadata: 7 (0.84%)
## Enterobase Database (5355 genomes/7775 metadata)
* concordance to provided metadata: 6022 (77.45%)
* discrepancies with provided metadata: 263 (3.38%)