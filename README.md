# ECTyper (an easy typer)
The `ectyper` wraps a standalone serotyping module.  
Support FASTA and FASTQ format.gi

# Dependencies:
* python (v3.6.2)
* bowtie2 (v2.2.6)
* SAMtools (v1.6)
* bcftools (v1.2)
* mash (v2.0)

# Basic usage
1. Run `pip [path_to_ectyper/]`
1. Put all you fasta/fastq file in one folder
1. `ectyper -i [dir_path or file_path]`
1. View result on console  

* If you want to enable species identification,
    1. Download refseq from https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh
    2. Put refseq into ectyper/Data/
    3. Specify `-s` as argument when executing ectyper