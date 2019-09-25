**v0.8.1**
* supports E.coli species detection via 10 short 1000 bp sequences based on E.coli core genomes
* contains signatures for the 178 O and 53 H antigen types
* non-E.coli isolates are identified using mash screen function against entire RefSeq database of all genomes


**v0.9.0**
* improves O-antigen serotyping coverage of complex samples that lack some O-antigen signatures (poor coverage,large fragmentation). 
* fall back to single O-antigen signature detection when both signature pairs are not found and species is E.coli
* verification for E.albertii species against RefSeq genomes
* addition of Quality Control flags in the output (as an extra column in the `results.tsv`) for ease of results interpretation
* species reporting to better resolve E.coli vs Shigella cases and contamination
* better complex cases handling and error recovery in cases of poor reference allele coverage. 
* serotype reporting based on multiple evidences including reference allele database and closest reference genome
* query length coverage deafult threshold lowered from 50% to 30% for truncated alleles. User is warned about potential false postive result


