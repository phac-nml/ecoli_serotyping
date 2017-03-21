#!/usr/bin/env python

import re
import src.genomeFunctions

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_record, args, results_dict):
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

    genome_name = src.genomeFunctions.get_genome_name(blast_record['sseqid'])


    serotype_results = {
        'otype':None,
        'htype':None,
        'flic':None,
        'flka':None,
        'flla':None,
        'flma':None,
        'wzx':None,
        'wzy':None,
        'wzm':None,
        'wzt':None
    }

    h_re_patterns = (
        # Look for any flagella match
        # We will sub-match for the gene if this hits
        re.compile('fl\w\w-(H\d+)__'),
        re.compile('fl\w\w_.+_(H\d+)$')
    )

    o_re_patterns = (
        # look for any somatic match
        # We will sub-match if it hits
        re.compile('wz\w-(O\d+)'),
        re.compile('>wz\w_.+_(O\d+)\w*$'),
        re.compile('gnd\|\d+_(O\d+)')
    )

    for rep in h_re_patterns:
        m = rep.search(blast_record['qseqid'])

        if m:
            serotype_results['otype']=m.group(1)


    return serotype_results