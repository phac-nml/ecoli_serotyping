#!/usr/bin/env python

import src.blastFunctions

"""
    Virulence prediction for E. coli
"""


def predict_virulence_factors(blast_record, args):
    """
    Entry point for virulence factor prediction
    :return: dict(
        otype = O
        htype = H
    )
    """

    # Initially check that the result passes the length / identity filters
    if not src.blastFunctions.record_passes_cutoffs(blast_record, args):
        return {}

    return {}