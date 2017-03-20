#!/usr/bin/env python

"""
    Virulence prediction for E. coli
"""


def predict_virulence_factors(blast_record):
    """
    Entry point for virulence factor prediction
    :return: dict(
        otype = O
        htype = H
    )
    """
    for alignment in blast_record.alignments:
        return alignment.title.strip()
