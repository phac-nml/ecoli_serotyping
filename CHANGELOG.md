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
* Corrected issues #76 and 77. Improved behaviour in case of no `--verify` switch is 
specified and MASH distance to RefSeq genomes fails
* Added additional test cases for non-E.coli genomes

**v1.0.0**
* Updated database of the reference alleles removing truncated and overlapping alleles for greater speed and accuracy
* Simplified report output removing confidence score categories
* Changed `--verify` flag logic. 
  * Now species verification via MASH distance to RefSeq genomes and E.coli specific alleles is done only if `--verify`
parameter is specified. 
  * If `--verify` is not specified, all input genomes are treated as E.coli without doing any species verification

**v2.0.0**
* Updated species identification module now based on GTDB + custom Escherichia and Shigella sketch covering all known bacterial species
* Implemented pathotyping covering 7 DEC *Escherichia coli* pathotypes (`DAEC`, `EAEC`, `EHEC`, `EIEC`, `EPEC`, `ETEC` and `STEC`) supporting simultaneous presence of multiple signatures (e.g. `ETEC/STEC`). Note that `EHEC` is reported as `EHEC-STEC` as this is a more severe subtype of `STEC`. 
* Implemented Shiga 1 and 2 toxin typing supporting multiple toxin signatures present in a single sample.
  * A total of 4 *stx1* subtypes are supported: `stx1a`, `stx1c`, `stx1d` and `stx1e`.
  * A total of 15 *stx2* subtypes are supported: `stx2a`, `stx2b`, `stx2c`, `stx2d`, `stx2e`, `stx2f`, `stx2g` ,`stx2h`, `stx2i`, `stx2j`, `stx2k`, `stx2l`, `stx2m`, `stx2n`, `stx2o`.
* new database of pathotypes and toxins in JSON clear transparent format composed of the key virulence factors based on both BioNumerics and literature sources  
* support for gzip compressed inputs `fastq.gz` and `fasta.gz` saving storage and increasing versatility
* other toxin typing covering enterohemolysin A (`ehxA`), hemolysin E (`hlyE`),  hemolysin A (`hlyA`)

  