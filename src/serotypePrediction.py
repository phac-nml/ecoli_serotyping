#!/usr/bin/env python

import re
import src.genomeFunctions

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_record, args):
    """
    Entry point for serotype prediction

    :param blast_record: = {'qseqid':la[0],
                        'qlen':la[1],
                        'sseqid':la[2],
                        'length':la[3],
                        'pident':la[4]
                        }
    :param args: argparse commandline options
    :return: dict = {
        otype = O
        htype = H
    }
    """

    # Initially check that the result passes the length / identity filters
    if not record_passes_cutoffs(blast_record, args):
        return {}


    # serotype header formats
    # >1__fliC__fliC-H1__1 AB028471.1;flagellin;H1
    # >2__flkA__flkA-H3__77 AB128916;flagellin;H3
    # >3__fllA__fllA-H44__82 #VALUE!;flagellin;H44
    # >4__flmA__flmA-H54__84 AB128918;flagellin;H54
    # >5__flnA__flnA-H17__85 CP002291;flagellin;H17
    # >6__wzm__wzm-O101-Gp15__86 AB812046.1;O-antigen ABC transporter permease;O101
    # >7__wzt__wzt-O101-Gp15__99 AB812046.1;ATP-binding protein;O101
    # >8__wzx__wzx-O1__112 GU299791.1;O antigen flippase;O1
    # >9__wzy__wzy-O1__314 GU299791.1;O antigen polyermase;O1
    # >8__wzx__wzx-O116var1__513 ERS085345;O antigen flippase;O116-like
    # >6__wzm__wzm-Onovel32__523 SRR2544773;O-antigen ABC transporter permease;Onovel32
    # >7__wzt__wzt-Onovel32__524 SRR2544773;ATP-binding protein;Onovel32
    # >fliC_54_AJ605766_H17
    # >wzx_36_DQ462205-FJ539194_O28ac
    # >gnd|1_O26
    #

    serotype_results = {
        # 'otype':None,
        # 'htype':None,
        # 'flic':None,
        # 'flka':None,
        # 'flla':None,
        # 'flma':None,
        # 'wzx':None,
        # 'wzy':None,
        # 'wzm':None,
        # 'wzt':None,
        # 'gnd':None
    }

    h_re_patterns = (
        # Look for any flagella match
        re.compile('(fl\w\w)-(H\d+)__'),
        re.compile('(fl\w\w)_.+_(H\d+)$')
    )

    o_re_patterns = (
        # look for any somatic match
        re.compile('(wz\w)-(O\d+)'),
        re.compile('^(wz\w)_.+_(O\d+)\w*$'),
        re.compile('(gnd)\|\d+_(O\d+)')
    )

    all_re_patters = (h_re_patterns, o_re_patterns)

    # Once a match has been found, we don't need to check any others
    # Checks for H-antigens first
    for pattern in all_re_patters:
        for rep in pattern:
            m = rep.search(blast_record['qseqid'])

            if m:
                serotype_results[m.group(1).lower()] = m.group(2)
                return serotype_results


def record_passes_cutoffs(blast_record, args):
    """
    For serotype prediction, need to ensure the blast hit is equal to or 
    greater than the cutoffs supplied
    
    :param blast_record: = {'qseqid':la[0],
                        'qlen':la[1],
                        'sseqid':la[2],
                        'length':la[3],
                        'pident':la[4]
                        }
    :param args: argparse commandline options
    :return: {True|False}
    """

    if (float(blast_record['qlen']) / float(blast_record[
        'length']) * 100) >= float(args.percentLength) and (
        float(blast_record['pident']) >= float(args.percentIdentity)):
        return True
    else:
        return False
