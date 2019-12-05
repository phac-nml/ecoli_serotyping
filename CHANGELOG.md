**v0.8.1**
* supports E.coli species detection via 10 short 1000 bp sequences based on E.coli core genomes
* contains signatures for the 178 O and 53 H antigen types
* non-E.coli isolates are identified using mash screen function against entire RefSeq database of all genomes


**v0.9.0**
* improved O-antigen serotyping coverage of complex samples that lack some O-antigen signatures. 
* fall back to single O-antigen signature detection when both signature pairs are not found and species is E.coli
* improved O-antigen identification favoring presence of both alleles (e.g. wzx and wzy) to support the final call. 
The sum of scores for both alleles of the same antigen is used in ranking now
* verification for E.albertii species against RefSeq genomes
* addition of Quality Control flags in the output (as an extra column in the `results.tsv`) for ease of results interpretation
* species reporting to better resolve E.coli vs Shigella vs E.albertii cases and contamination
* serotype prediction for all Escherichia genus to ease classification of cryptic spqeices
* automatic download of the RefSeq sketch and associated meta data every 6 months for improved species identification
* improved species identification for the FASTQ files. All raw reads are used for species identification
* better complex cases handling and error recovery in cases of poor reference allele coverage. 
* serotype reporting based on multiple evidences including reference allele database and closest reference genome
* query length coverage default threshold lowered from 50% to 10% for truncated alleles. This greatly improved sensitivity.
Did not found any noticeable negative effect on specificity and accuracy based on EnteroBase, PNC2018-2019 and internal datasets
* users is warned about potential false postive result in case of non-E.coli species or assembly errors
* wrote additional unit tests


**v0.9.1**
* Implemented better species identification error handling for cases when accession number
is not found in the assembly stats. This is especially important for custom
MASH sketches that might have accession numbers not found in `assembly_summary.txt`  stats file
* Added `--assemblystats` command line flag for cases when user will supply its own MASH sketch for species 
identification
* Replaced the `assembly_stats.txt` file with the GenBank equivalent as
some of RefSeq sequences found in the MASH sketch were removed and GenBank is a superset 
that covers all accessions found in sketch
* added `--assemblername` flag to allow to assemble raw reads via 
`Shovill`, `Spades`, `Canu`, `Unicycler` assembler.