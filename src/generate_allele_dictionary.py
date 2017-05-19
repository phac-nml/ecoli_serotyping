#!/usr/bin/env python

"""
Generate a JSON dictionary of all known headers and their allele output.
ectyper will read this in and store as a dictionary for exact lookup.
Blast stops at first space, so to include entire header we will substitute
all non-word characters with '_'
Output format in JSON:
[
    "O":{'8__wzx__wzx-O1__112 GU299791.1;O antigen flippase;O1' : 'O1'}
        ...
    "H":{'1__fliC__fliC-H1__1 AB028471.1;flagellin;H1' : 'H1'}
    ...
]

Output file is always to Data/allele_types.json
"""

import os
import json
import re
import collections
import definitions

# The input file of O and H antigens stored in the DATA directory
input_file = os.path.join(definitions.ROOT_DIR, 'Data/EcOH.fasta')
file_handle = open(input_file, 'r')

# We will iterate through every line, looking only at the headers
# and storing them verbatim under either O or H, with allele as data.
o_dictionary = collections.defaultdict(str)
h_dictionary = collections.defaultdict(str)

for line in file_handle:
    clean_line = line.strip()

    # single_line = re.sub('[\s]', '_', clean_line)

    # look for fasta headers only
    # match gene name, O or H, and allele type
    # note that this does not match those designated as "Onovel" or "Hnovel"
    # those will have to be manually curated / update
    regex = re.compile('^((>\d+__|>)([a-zA-Z]+).+([_;])(([OH])\d+)$)')
    m = regex.search(clean_line)

    if m:
        # we need to match what blast will match for comparing
        blast_match_regex = re.compile('^>(\S+)')
        bm = blast_match_regex.search(clean_line)
        # bm.group(1) is the key we will store by

        # in the match, the groups are as follows
        # m.group(1) = entire header
        # m.group(2) = beginning matches, >\d+__|>
        # m.group(3) = gene name
        # m.group(4) = _|;
        # m.group(5) = allele name
        # m.group(6) = O|H

        # the regex includes the '>' carat, so we need to use a slice that
        # excludes it
        # now check whether O/H and store alleles and gene names
        if m.group(6) == 'O':
            o_dictionary[bm.group(1)] = {'gene': m.group(3),
                                         'allele': m.group(5)}
        elif m.group(6) == 'H':
            h_dictionary[bm.group(1)] = {'gene': m.group(3),
                                         'allele': m.group(5)}
        else:
            print("ERROR, no O/H for " + clean_line)

# create a single dictionary for JSON output
allele_dictionary = {'O': o_dictionary,
                     'H': h_dictionary}

# output the JSON to the Data/allele_types.json file
out_file = os.path.join(definitions.ROOT_DIR, 'Data/allele_types.json')
out_handle = open(out_file, 'w')
out_handle.write(json.dumps(allele_dictionary))
