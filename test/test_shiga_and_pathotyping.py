import sys
import pytest
import tempfile
import os
from ectyper import ectyper, definitions
import subprocess
import pandas as pd
import logging
import re

TEST_ROOT = os.path.dirname(__file__)
LOG=logging.getLogger("TEST")
LOG.setLevel(logging.INFO)

def set_input(input,
              percent_iden=None,
              verify=True,
              output=tempfile.mkdtemp(),
              cores=1,
              debug=False,
              pathotype = False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input,
            '-c', str(cores),
            ]

    if percent_iden:
        args += ['-d', str(percent_iden)]
    if verify:
        args += ['--verify']
    if output:
        args += ['-o', output]
    if debug:
        args+=['--debug']    
    if pathotype:
        args+=['--pathotype']    

    sys.argv[1:] = args


def test_single_stx2_subtyping(caplog):   
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,'Data/EscherichiaO28H5.fasta') 
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, verify=True, debug=False, output=tmpdir, pathotype=True)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
        rows = outfp.readlines()
      
    columnslist = rows[0].split('\t')   
    secondrow = rows[1]
    secondrow_values=secondrow.split('\t')
    assert len(secondrow_values) == len(columnslist), f"Unequal number of columns {len(columnslist)} and values {len(secondrow_values)}"
    assert "STEC" in secondrow   
    assert "stx2a" in secondrow

def test_stx1_stx2_subtyping_pathotyping(caplog): 
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,'Data/Escherichia.fna')  
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, verify=True, debug=True, output=tmpdir, pathotype=True)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
        secondrow = outfp.readlines()[1]
    assert "EHEC" in secondrow   
    assert "stx1a;stx2a" in secondrow
    assert "AP010958.1;AP010958.1" in secondrow


def test_multi_stx_non_overlap_ranges(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,'Data/CP041431_STEC316.fasta.gz')  
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, verify=True, debug=False, output=tmpdir, pathotype=True)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
        secondrow = outfp.readlines()[1]
    assert "STEC" in secondrow   
    assert "stx2e" in secondrow
    assert "stx2k" in secondrow    

def test_multi_stx_non_overlap_different_contigs(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,'Data/SRR7947260.fasta.gz')  
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, verify=True, debug=False, output=tmpdir, pathotype=True)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
        secondrow = outfp.readlines()[1]
    assert "ETEC/STEC" in secondrow   
    assert "stx2a" in secondrow
    assert "stx2g" in secondrow   
    assert "contig00064;contig00074" in secondrow


def test_multi_stx_overlap_same_contig(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,'Data/SRR7612273.fasta.gz')  
    tmpdir = tempfile.mkdtemp()
    set_input(input=file, cores=4, verify=True, debug=True, output=tmpdir, pathotype=True)
    ectyper.run_program()
    with open(os.path.join(tmpdir,"output.tsv")) as outfp:
        secondrow = outfp.readlines()[1]
    assert "STEC" in secondrow   
    assert "stx2a" in secondrow
    assert "stx2d" in secondrow   
    assert "contig00078;contig00078" in secondrow

