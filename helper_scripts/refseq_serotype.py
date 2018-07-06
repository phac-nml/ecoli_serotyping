#!/usr/bin/env python

import sys
import re

with open(sys.argv[1], "r") as fh:
    for line in fh:
        la = line.split()
