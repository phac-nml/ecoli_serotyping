[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=ectyper)
![GitHub release (latest SemVer)](https://img.shields.io/github/v/release/phac-nml/ecoli_serotyping)
![Conda](https://img.shields.io/conda/dn/bioconda/ectyper)
[![PyPI version](https://badge.fury.io/py/ectyper.svg)](https://pypi.org/project/ectyper/)
![GitHub issues](https://img.shields.io/github/issues/phac-nml/ecoli_serotyping)
![Docker Pulls](https://img.shields.io/docker/pulls/kbessonov/ectyper)

# ECTyper (an easy typer)
`ECTyper` is a standalone versatile serotyping module for _Escherichia coli_. It supports both _fasta_ (assembled) and _fastq_ (raw reads) file formats.
The tool provides convenient species identification coupled to quality control module giving a complete, transparent and reference laboratories suitable report on E.coli serotyping.

# Introduction
*Escherichia coli* is a priority foodborne pathogen of public health concern and popular model organism. Phenotypic characterization such as serotyping, toxin typing and pathotyping provide critical information for surveillance and outbreak detection activities and research including source attribution, outbreak cluster assignment, pathogenicy potential, risk assessement and others. 

`ECTyper` uses whole-genome sequencing (WGS) for E.coli characterizion including species identification, *in silico* serotyping covering O and H antigens, Shiga toxin typing and DEC pathotyping. It is a versatile, scallable, easy to use tool allowing to obtain key information on E.coli accepting both raw and assembled inputs.

As WGS becomes standard within public health and research laboratories, it is important to harness the high thourghput and resolution potential of this technology providing accurate and rapid at scale typing of E.coli both in public health, clinical and research contexts.

## Citation
Bessonov, Kyrylo, Chad Laing, James Robertson, Irene Yong, Kim Ziebell, Victor PJ Gannon, Anil Nichani, Gitanjali Arya, John HE Nash, and Sara Christianson. "ECTyper: in silico Escherichia coli serotype and species prediction from raw and assembled whole-genome sequence data." Microbial genomics 7, no. 12 (2021): 000728. [https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000728](https://www.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000728)

## Contact
For any questions, issues or comments please make a Github issue or reach out to [Kyrylo Bessonov](kyrylo.bessonov@phac-aspc.gc.ca).

# Installation
Multiple installation options are available depending on the user context and needs. The most convinient installation is as a `conda` package as it will install all required dependencies. 

### Images
Docker and Singularity images are also available from [https://biocontainers.pro/tools/ectyper](https://biocontainers.pro/tools/ectyper) that could be useful for NextFlow or hassle-free deployment

### Databases
ECTyper uses multiple databases 
  - the species identification database is available from [https://zenodo.org/records/10211569](https://zenodo.org/records/10211569)
  - the O and H antigen allele sequences are stored in [ectyper_alleles_db.json](ectyper/Data/ectyper_alleles_db.json)
  - the toxin and pathotype signature marker sequences are stored in [ectyper_patho_stx_toxin_typing_database.json](ectyper/Data/ectyper_patho_stx_toxin_typing_database.json)

## Option 1: As a conda package
Optionally if you do not have a conda environment, get and install `miniconda` or `anaconda`:

  ```
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
  bash miniconda.sh -b -p $HOME/miniconda
  echo ". $HOME/miniconda/etc/profile.d/conda.sh" >> ~/.bashrc
  source ~/.bashrc
  ```
    
Install the latest `ectyper` conda package from a `bioconda` channel 
	
  ```
  conda install -c bioconda ectyper
  ```

## Option 2: Install using pip
Install using `pip3` utility including python but missing on [non-python dependencies](#dependencies)
  ```
  pip3 install ectyper
  ```
## Option 3: From source code
Second option is to install from the source allowing to excercise maximum control over installation process.

Install dependencies. On Ubuntu distro run
  ```
  apt install samtools bowtie2 mash bcftools ncbi-blast+ seqtk
  ```

Install python dependencies via `pip`:
  ```
  pip3 install pandas biopython
  ```
Clone the repository or checkout a particular release (e.g `v1.0.0`, `v2.0.0` etc.):
  ```
  git clone https://github.com/phac-nml/ecoli_serotyping.git
  git checkout v1.0.0 #optionally checkout a specific release version
  ```

Finally, install ectyper  
```
python3 setup.py install # option 1
pip3 install .   # option 2
```
## Compatibility
### Dependencies:
- python >= 3.5
- bcftools >= 1.8
- blast == 2.7.1
- seqtk >= 1.2
- samtools >= 1.8
- bowtie2 >= 2.3.4.1
- mash >= 2.0

### Python packages:
- biopython >= 1.70
- pandas >= 0.23.1
- requests >= 2.0

# Getting started
## Basic Usage
1. Put the fasta/fastq files for serotyping analyses in one folder (concatenate paired raw reads files if you would like them to be considered a single entity)
1. `ectyper -i [file path] -o [output_dir]`
1. View the results on the console or in `cat [output folder]/output.csv`

## Example Input Scenarios
* `ectyper -i ecoliA.fasta`  for a single file (the output folder will be named using `ectyper_<date>_<time>` pattern)
* `ectyper -i ecoliA.fasta -o output_dir` for a single file, results stored in `output_dir` folder
* `ectyper -i ecoliA.fasta ecoliB.fastq ecoli_folder/`	for multiple files and directory separated by space
* `ectyper -i ecoliA.fasta ecoliB.fastq,ecoliC.fna`
* `ectyper -i ecoli_folder`	scan for input files in a folder and subdirectories (all files in the folder will be checked by the tool)
* `ectyper -i ecoli_folder/*.fasta` scan for FASTA input files in a folder and subdirectories

## Advanced Usage
```
usage: ectyper [-h] [-V] -i INPUT [-c CORES] [-opid PERCENTIDENTITYOTYPE]
               [-hpid PERCENTIDENTITYHTYPE] [-oplen PERCENTLENGTHOTYPE]
               [-hplen PERCENTLENGTHHTYPE] [--verify] [-o OUTPUT] [-r REFSEQ] [-s] [--debug]
               [--dbpath DBPATH]

ectyper v2.0.0 database v1.0 Prediction of Escherichia coli serotype from raw reads or assembled
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


## Configuration and fine-tunning parameters
`ECTyper` requires minimum options to run (`-i` and `-o`) but allows for extensive configuration to accomodate wide variaty of typing scenarios

| Parameter|      Explanation                                                 | Usage scenario                                                                    |
|----------|:----------------------------------------------------------------:|:----------------------------------------------------------------------------------:
| `-opid`  |  Specify minimum `%identity` threshold just for O antigen match| Poor coverage of O antigen genes or for exploratory work (recommended value is 90)  |
| `-opcov` |  Minimum `%covereage` threshold for a valid match against reference O antigen alleles | Poor coverage of O antigen genes and a user wants to get O antigen call regardless (recommend value is 95)|
| `-hpid`  |  Specify minimum `%identity` threshold just for H antigen match| Poor coverage of O antigen genes or for exploratory work (recommend value is 95)    |
| `-hpcov` |  Minimum `%covereage` threshold for a valid match against reference H antigen alleles | Poor coverage of O antigen genes and a user wants to get O antigen call regardless (recommend value is 95)|
|`--verify`|  Verify species of the input and run QC module providing information on the reliability of the result and any typing issues | User not sure if sample is E.coli and wants to obtain if serotype prediction is of sufficient quality for reporting purposes|
| `-r`     |  Specify custom MASH sketch of reference genomes that will be used for species inference | User has a new assembled genome that is not available in NCBI RefSeq database. Make sure to add metadata to `assembly_summary_refseq.txt` and provide custom accession number that start with `GCF_` prefix|
|`--dbpath`|  Provide custom appended database of O and H antigen reference alleles in JSON format following structure and field names as default database `ectyper_alleles_db.json` | User wants to add new alleles to the alleles database to improve typing performance |

# Data Input
Both raw and assembled reads are accepted in FASTA and FASTQ formats from any sequencing platform. The tool was designed for single sample inputs, but was shown to work on multi-taxa metagenomic raw reads FASTQ inputs.

# Data Output
The output of the tool is stored in text files with the main report stored in `output.tsv` tab-delimited text file.

The BLASTN hits of the O and H antigen database are stored in `blastn_output_alleles.txt` tab-delimited file.

The log messages are stored in `ectyper.log` text file
```
{out folder name}
├── blastn_output_alleles.txt
├── ectyper.log
└── output.tsv
```

## Quality Control (QC) module
To provide an easier interpretation of the results and typing metrics, following QC codes were developed. 
These codes allow to quickly filter "reportable" and "non-reportable" samples. The QC module is tightly linked to ECTyper allele database, specifically, `MinPident` and `MinPcov` fields.
For each reference allele minimum `%identity` and `%coverage` values were determined as a function of potential "cross-talk" between antigens (i.e. multiple potential antigen calls at a given setting).
The QC module covers the following serotyping scenarios. More scenarios might be added in future versions depending on user needs.

| QC flag          |      Explanation                                          | 
|------------------|:-----------------------------------------------------------------|
|PASS (REPORTABLE) |Both O and H antigen alleles meet min `%identity` or `%coverage` thresholds (ensuring no antigen cross-talk) and single antigen predicted for O and H|
|FAIL (-:- TYPING) |Sample is E.coli and O and H antigens are not typed. Serotype:  -:- |
|WARNING MIXED O-TYPE|A mixed O antigen call is predicted requiring wet-lab confirmation |
|WARNING (WRONG SPECIES)| A sample is non-E.coli (e.g. E.albertii, Shigella, etc.) based on RefSeq assemblies|
|WARNING (-:H TYPING)| A sample is E.coli and O antigen is not predicted (e.g. -:H18)|
|WARNING (O:- TYPING)| A sample is E.coli and O antigen is not predicted (e.g. O17:-)|
|WARNING (O NON-REPORT)|O antigen alleles do not meet min %identity or %coverage thresholds|
|WARNING (H NON-REPORT)|H antigen alleles do not meet min %id or %cov thresholds|
|WARNING (O and H NON-REPORT)| Both O and H antigen alleles do not meet min %identity or %coverage thresholds|
 
## Report format
`ECTyper` capitalizes on a concise minimum output coupled to easy results interpretation and reporting. `ECTyper v1.0` serotyping results are available in a tab-delimited `output.tsv` file consisting of the 16 columns listed below:

1. **Name**: Sample name (usually a unique identifier) 
2. **Species**: the species column provides valuable species identification information in case of inadvertent sample contamination or mislabelling events
3. **O-type**: O antigen
4. **H-type**: H antigen
5. **Serotype**: Predicted O and H antigen(s)
6. **QC**: The Quality Control value summarizing the overall quality of prediction
7. **Evidence**: How many alleles in total used to both call O and H antigens
8. **GeneScores**: ECTyper O and H antigen gene scores in 0 to 1 range
9. **AllelesKeys**: Best matching `ECTyper` database allele keys used to call the serotype  
10. **GeneIdentities(%)**: `%identity` values of the query alleles
11. **GeneCoverages(%)**: `%coverage` values of the query alleles
12. **GeneContigNames**: the contig names where the query alleles were found
13. **GeneRanges**: genomic coordinates of the query alleles
14. **GeneLengths**: allele lengths of the query alleles
15. **Database**: database release version and date
16. **Warnings**: any additional warnings linked to the quality control status or any other error message(s). 


Selected columns from the `ECTyper` typical report are shown below. 

|Name|Species|Serotype|Evidence|QC|GeneScores|AlleleKeys|GeneIdentities(%)   |    GeneCoverages(%)   |     GeneContigNames| GeneRanges   |   GeneLengths  |   Database |        Warnings|
|------|:------|:-------|:--|:---------|:------------|:-----|:-----|:----|:----|:----|:----|:---|:--|
|15-520|Escherichia coli|O174:H21|Based on 3 allele(s)|PASS (REPORTABLE)|wzx:1; wzy:1; fliC:1;|O104-5-wzx-origin;O104-13-wzy;H7-6-fliC-origin;|100;100;100;|    100;100;100;|contig00049;contig00001;contig00019;|   22302-23492;178-1290;6507-8264;| 1191;1113;1758;| v1.0 (2020-05-07)  |     -  |
EC20151709|Escherichia coli|O157:H43|Based on 3 allele(s)|PASS (REPORTABLE)|wzx:1;wzy:0.999;fliC:1|O157-5-wzx-origin;O157-9-wzy-origin;H43-1-fliC-origin;|100;99.916;99.934;   |   100;100;100;  |  contig00002;contig00002;contig00003; |   62558-63949;64651-65835;59962-61467;   | 1392;1185;1506; |v1.0 (2020-05-07)  |     - |


FAQs

## FAQ

**Does ECTyper can be run on multiple samples in a directory?**

ECTyper proves flexible ways to specify inputs located in different locations. One can provide multiple paths to several directories separated by space. In addition, one can specify file type to look for in a given diretory(ies). Note that paths that contain a star `*` symbol would only look for files in specified directory and would not look in subdirectories. For example,

- Process all files in `folder1` and `folder2` directories and file `sample.fasta` located in `folder3` 

    `ectyper -i folder1/ folder2/ folder3/sample.fasta -o ectyper_results` 
- Process all fasta files in `folder1` and all fastq files in `folder2`. All sub-directories in those 2 folders will be ignored. To process those sub-folders either specify path to them or provide paths to directories without the `*` wildcard symbol. 
  
  `ectyper -i folder1/*.fasta folder2/*.fastq` 

**Why ECTyper sometimes provides serotype results separated by forward  slash / for O-antigen**

Some O-antigens display very high degree of homology and are very hard to discern even using wet-lab agglutination assays. Even using both `wzx` and `wzy` genes it is not possible to reliably resolve those O-antigens. The 16 high similarity groups were identified by [Joensen, Katrine G., et al.](https://journals.asm.org/doi/full/10.1128/jcm.00008-15). Thus, if a given O-antigen is a member of any of those high similarity groups, all potential O-antigens are reported separated by `/` such as group 9 reporeted as `O17/O44/O73/O77/O106`.


# Availability

|Resource|Description|Type|
|--------|:----------|:---|
|[PyPI](https://pypi.org/project/ectyper/)| PyPI pacakge that could be installed via `pip` utility|Terminal|
|[Conda](https://anaconda.org/bioconda/ectyper)   | Conda package available from BioConda channel|Terminal|
|[Docker](https://hub.docker.com/r/kbessonov/ectyper)| Images containing completely initialized ECTyper with all dependencies |Terminal|
|[Singluarity](https://biocontainers.pro/tools/ectyper) | Images containing completely initialized ECTyper with all dependencies |Terminal|
|[GitHub](https://github.com/phac-nml/ecoli_serotyping) | Install source code as any Python package|Terminal|
|[Galaxy ToolShed](https://toolshed.g2.bx.psu.edu/view/nml/ectyper/)| Galaxy wrapper available for installation on a private/public instance|Web-based|
|[Galaxy Europe](https://usegalaxy.eu/root?tool_id=ectyper)| Galaxy public server to execute your analysis from anywhere|Web-based| 
|[IRIDA plugin](https://github.com/phac-nml/irida-plugin-ectyper)| IRIDA instances could easily install additional pipeline|Web-based|

# Legal and Compliance Information

Copyright Government of Canada 2024

Written by: National Microbiology Laboratory, Public Health Agency of Canada

Licensed under the Apache License, Version 2.0 (the "License"); you may not use this work except in compliance with the License. You may obtain a copy of the License at:

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.