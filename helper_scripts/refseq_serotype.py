#!/usr/bin/env python

"""
Given the `assembly_summary_refseq.txt` file, create a list of Escherichia
coli genomes that need to be downloaded, and then download them to the
specified directory.

Columns in the file:
--------------------
[0] assembly_accession
[1] bioproject
[2] biosample
[3] wgs_master
[4] refseq_category
[5] taxid
[6] species_taxid
[7] organism_name
[8] infraspecific_name
[9] isolate
[10] version_status
[11] assembly_level
[12] release_type
[13] genome_rep
[14] seq_rel_date
[15] asm_name
[16] submitter
[17] gbrs_paired_asm
[18] paired_asm_comp
[19] ftp_path
[20] excluded_from_refseq
[21] relation_to_type_material



"""

import sys
import re
import csv

results = {}
with open(sys.argv[1], 'r', encoding='utf-8') as csvfh:
    csvfh.readline()
    csvfh.readline()
    csvreader = csv.reader(csvfh, delimiter='\t')
    for row in csvreader:
        m = re.search(r"Escherichia coli", row[7])
        if m:
            ftp_link = row[19]
            mftp = re.search(r"(GCF_.*$)", ftp_link)
            if mftp:
                ftp_link = ftp_link + '/' + mftp.group(1) + "_genomic.fna.gz"

            msero = re.search(r"O\d+", row[7])
            if msero:
                results[row[0]]={
                    'O':msero.group(0),
                    'ftp':ftp_link
                }

            mflag = re.search(r":(H\d+)", row[7])
            if mflag:
                if row[0] not in results:
                    results[row[0]]={}

                results[row[0]]['H'] = mflag.group(1)
                results[row[0]]['ftp'] = ftp_link

for k,v in results.items():
    outline=[k]

    for kk, vv in reversed(sorted(v.items())):
        outline.append(vv)
    print("\t".join(outline))
