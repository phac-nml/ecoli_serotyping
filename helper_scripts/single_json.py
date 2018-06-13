#!/usr/env/bin python

import os
import json
import sys
from Bio import SeqIO

"""
Add the fasta sequence to the JSON to have a single source of truth
"""

fasta_file = sys.argv[1]
json_file = sys.argv[2]
out_file = sys.argv[3]

# Gather the fasta sequences
fasta_dict = {}
with open(fasta_file, 'r') as fastafh:
    for record in SeqIO.parse(fastafh, "fasta"):
        fasta_dict[str(record.id)] = str(record.seq)


# Combine the fasta sequences with the other data and print it out
with open(json_file, 'r') as jsonfh:
    json_data = json.load(jsonfh)

    for k,v in fasta_dict.items():
        if k in json_data["O"]:
            json_data["O"][k]["seq"] = v

        elif k in json_data["H"]:
            json_data["H"][k]["seq"] = v

        else:
            print("{} does not exist".format(k), file=sys.stderr)

    with open(out_file, 'w') as jsonofh:
        jsonofh.write(json.dumps(json_data))


