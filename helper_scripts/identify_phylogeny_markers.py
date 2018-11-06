#!/usr/bin/env python

import sys
import pandas as pd
import numpy as np
from Bio import SeqIO
import re

binary_file = sys.argv[1]
fasta_file = sys.argv[2]
output_file = sys.argv[3]

FRAGMENT_CUTOFF=206

fragments_to_keep = []
with open(binary_file, 'r') as bfh:

    #histogram selection of "majority core"
    df = pd.read_csv(bfh, delimiter='\t', header=0, index_col=0)
    rowsums = df.sum(axis=1)
    print(len(rowsums))

    h = np.histogram(rowsums, 100)
    #print(h)

    #keep only the fragments with a rowsum > FRAGMENT_CUTOFF
    final_fragments_df = df.loc[rowsums >= FRAGMENT_CUTOFF]
    print(final_fragments_df.index)


    print(len(final_fragments_df))
    rowsum_check = final_fragments_df.sum(axis=1)
    #print(rowsum_check)

    #open the fasta file, keeping only the segments that were >= FRAGMENT_CUTOFF

    with open(fasta_file, "r") as ffh:
        counter = 1
        with open(output_file, 'w') as ofh:
            for record in SeqIO.parse(ffh, "fasta"):
                desc = record.description
                m = re.match(r"lcl\|(\d+)", desc)
                if m:
                    # the dataframe index is of int64
                    # need to cast the match, or it won't match the index
                    fragment_name = int(m.group(1))

                    if fragment_name in final_fragments_df.index:
                        record.id = "dominant_core_" + str(counter)
                        record.description = ""
                        SeqIO.write(record, ofh, "fasta")
                        counter += 1



