#!/usr/bin/env python

"""
Given the `assembly_summary_refseq.txt` file, create a list
"""


import sys
import re

with open(sys.argv[1], "r") as fh:
    for line in fh:
        la = line.split()
