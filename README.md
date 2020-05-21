[![Master branch build status](https://api.travis-ci.org/phac-nml/ecoli_serotyping.svg?branch=master "Master Build Status")](https://travis-ci.org/phac-nml/ecoli_serotyping)

# ECTyper (an easy typer)
**ectyper** is a standalone versatile serotyping module for _Escherichia coli_. It supports both _fasta_ (assembled) and _fastq_ (raw reads) file formats.
The tool provides convenient species identification coupled to quality control module giving a complete, transparent and reference laboratories suitable report on E.coli serotpying.


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
1. If you do not have conda environment, get and install `miniconda` or `anaconda`:
    1. `wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh`
    1. `bash miniconda.sh -b -p $HOME/miniconda`
    1. `echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc`
    1. `source ~/.bashrc`
2. Install ectyper either  
    * Directly via `conda` (simplest option)
    	1. `conda install -c bioconda ectyper` 
    * As a source code through `github` repository
    	1. Install dependencies
          `pandas samtools bowtie2 mash bcftools biopython pytest blast seqtk tqdm python=3.6`
    	1. Clone the repository (most recent build) or checkout a particular release
          `git clone https://github.com/galaxyproject/galaxy.git`
          `git checkout v1.0.0`
    	1. Install ectyper
          `python setup.py install`

# Basic Usage
1. Put the fasta/fastq files for serotyping analyses in one folder (concatenate paired raw reads files if you would like them to be considered a single entity)
1. `ectyper -i [file path] -o [output_dir]`
1. View the results on the console or in `cat [output folder]/output.csv`

# Example Usage
* `ectyper -i ecoliA.fasta`  for a single file
* `ectyper -i ecoliA.fasta -o output_dir` for a single file, results stored in `output_dir`
* `ectyper -i ecoliA.fasta,ecoliB.fastq,ecoliC.fna`	for multiple files  (comma-delimited)
* `ectyper -i ecoli_folder`	for a folder (all files in the folder will be checked by the tool)

# Advanced Usage
```
usage: ectyper [-h] [-V] -i INPUT [-c CORES] [-opid PERCENTIDENTITYOTYPE]
               [-hpid PERCENTIDENTITYHTYPE] [-oplen PERCENTLENGTHOTYPE]
               [-hplen PERCENTLENGTHHTYPE] [--verify] [-o OUTPUT] [-r REFSEQ] [-s] [--debug]
               [--dbpath DBPATH]

ectyper v1.0.0 database v1.0 Prediction of Escherichia coli serotype from raw reads or assembled
genome sequences. The default settings are recommended.

optional arguments:
  -h, --help            show this help message and exit
  -V, --version         show program's version number and exit
  -i INPUT, --input INPUT
                        Location of E. coli genome file(s). Can be a single file, a comma-
                        separated list of files, or a directory
  -c CORES, --cores CORES
                        The number of cores to run ectyper with
  -opid PERCENTIDENTITYOTYPE, --percentIdentityOtype PERCENTIDENTITYOTYPE
                        Percent identity required for an O antigen allele match [default 90]
  -hpid PERCENTIDENTITYHTYPE, --percentIdentityHtype PERCENTIDENTITYHTYPE
                        Percent identity required for an H antigen allele match [default 95]
  -oplen PERCENTLENGTHOTYPE, --percentLengthOtype PERCENTLENGTHOTYPE
                        Percent length required for an O antigen allele match [default 95]
  -hplen PERCENTLENGTHHTYPE, --percentLengthHtype PERCENTLENGTHHTYPE
                        Percent length required for an H antigen allele match [default 50]
  --verify              Enable E. coli species verification
  -o OUTPUT, --output OUTPUT
                        Directory location of output files
  -r REFSEQ, --refseq REFSEQ
                        Location of pre-computed MASH RefSeq sketch. If provided, genomes
                        identified as non-E. coli will have their species identified using MASH.
                        For best results the pre-sketched RefSeq archive
                        https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh is
                        recommended
  -s, --sequence        Prints the allele sequences if enabled as the final columns of the
                        output
  --debug               Print more detailed log including debug messages
  --dbpath DBPATH       Path to a custom database of O and H antigen alleles in JSON format. Check
                        Data/ectyper_database.json for more information
```

# Fine tunning parameters
ECTyper requires minimum options to run (`-i` and `-o`) but allows for extensive configuration to accomodate wide variaty of typing scenarios

| Parameter|      Explanation                                                 | Usage scenario                                                                    |
|----------|:----------------------------------------------------------------:|:----------------------------------------------------------------------------------:
| `-opid`  |  Specify minimum `%identity` threshold just for O antigen match| Poor coverage of O antigen genes or for exploratory work (recommended value is 90)  |
| `-opcov` |  Minimum `%covereage` threshold for a valid match against reference O antigen alleles | Poor coverage of O antigen genes and a user wants to get O antigen call regardless (recommend value is 95)|
| `-hpid`  |  Specify minimum `%identity` threshold just for H antigen match| Poor coverage of O antigen genes or for exploratory work (recommend value is 95)    |
| `-hpcov` |  Minimum `%covereage` threshold for a valid match against reference H antigen alleles | Poor coverage of O antigen genes and a user wants to get O antigen call regardless (recommend value is 95)|
|`--verify`|  Verify species of the input and run QC module providing information on the reliability of the result and any typing issues | User not sure if sample is E.coli and wants to obtain if serotype prediction is of sufficient quality for reporting purposes|
| `-r`     |  Specify custom MASH sketch of reference genomes that will be used for species inference | User has a new assembled genome that is not available in NCBI RefSeq database. Make sure to add metadata to `assembly_summary_refseq.txt` and provide custom accession number that start with `GCF_` prefix|
|`--dbpath`|  Provide custom appended database of O and H antigen reference alleles in JSON format following structure and field names as default database `ectyper_alleles_db.json` | User wants to add new alleles to the alleles database to improve typing performance |


# Quality Control (QC) module
To provide an easier interpretation of the results and typing metrics, following QC codes were developed. 
These codes allow to quickly filter "reportable" and "non-reportable" samples. The QC module is tightly linked to ECTyper allele database, specifically, `MinPident` and `MinPcov` fields.
For each reference allele minimum `%identity` and `%coverage` values were determined as a function of potential "cross-talk" between antigens (i.e. multiple potential antigen calls at a given setting).
The QC module covers the following serotyping scenarios. More scenarios might be added in future versions depending on user needs.

| QC flag          |      Explanation                                                 | 
|------------------|:-----------------------------------------------------------------|
|PASS (REPORTABLE) |Both O and H antigen alleles meet min `%identity` or `%coverage` thresholds (ensuring no antigen cross-talk) and single antigen predicted for O and H|
|FAIL (-:- TYPING) |Sample is E.coli and O and H antigens are not typed. Serotype:  -:- |
|WARNING MIXED O-TYPE|A mixed O antigen call is predicted requiring wet-lab confirmation |
|WARNING (WRONG SPECIES)| Sample is non-E.coli (e.g. E.albertii, Shigella, etc.) based on RefSeq assemblies|
|WARNING (-:H TYPING)| Sample is E.coli and O antigen is not predicted (e.g. -:H18)|
|WARNING (O:- TYPING)| Sample is E.coli and O antigen is not predicted (e.g. O17:-)|
|WARNING (O NON-REPORT)|O antigen alleles do not meet min %identity or %coverage thresholds|
|WARNING (H NON-REPORT)|H antigen alleles do not meet min %id or %cov thresholds|
|WARNING (O and H NON-REPORT)| Both O and H antigen alleles do not meet min %identity or %coverage thresholds|


 


# Availability
1. Terminal: Source code from this repository
1. Terminal: Conda package available from BioConda channel at [https://anaconda.org/bioconda/ectyper](https://anaconda.org/bioconda/ectyper)
1. Terminal: Docker/Singularity images (coming soon) 
1. Web-based: Galaxy wrapper available for installation from [Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/nml/ectyper/)
1. Web-based: IRIDA plug-in available from [https://github.com/phac-nml/irida-plugin-ectyper](https://github.com/phac-nml/irida-plugin-ectyper)
