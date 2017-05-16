#!/usr/bin/env python

"""
Generate a JSON dictionary of all known headers and their allele output.
ectyper will read this in and store as a dictionary for exact lookup.
Output format in JSON:
[
    "O":{'8__wzx__wzx-O1__112 GU299791.1;O antigen flippase;O1' : 'O1'}
        ...
    "H":{'1__fliC__fliC-H1__1 AB028471.1;flagellin;H1' : 'H1'}
    ...
]
"""

import os
import json
import re
import collections
import definitions

# The input file of O and H antigens stored in the DATA directory
input_file = os.path.join(definitions.ROOT_DIR,  'Data/EcOH.fasta')

file_handle = open(input_file, 'r')

# We will iterate through every line, looking only at the headers
# and storing them verbatim under either O or H, with allele as data.
o_dictionary = collections.defaultdict(str)
h_dictionary = collections.defaultdict(str)

for line in file_handle:
    clean_line = line.strip()

    # look for fasta headers only
    regex = re.compile('^>(.+(_|;)((O|H)\d+)$)')
    m = regex.search(clean_line)

    if m:
        #now check whether O/H and store alleles
        if m.group(4) == 'O':
            o_dictionary[m.group(1)]=m.group(3)
        elif m.group(4) == 'H':
            h_dictionary[m.group(1)]=m.group(3)
        else:
            print("ERROR, no O/H for " + clean_line)

#create a single dictionary for JSON output
allele_dictionary = {'O':o_dictionary,
                     'H':h_dictionary}

#output the JSON to the Data/allele_types.json file
out_file = os.path.join(definitions.ROOT_DIR, 'Data/allele_types.json')
out_handle = open(out_file, 'w')
out_handle.write(json.dumps(allele_dictionary))
