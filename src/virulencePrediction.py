#!/usr/bin/env python

"""
    Virulence prediction for E. coli
"""


def predict_virulence_factors(blast_record, args, results_dict):
    """
    Entry point for virulence factor prediction
    :return: dict(
        otype = O
        htype = H
    )
    """

