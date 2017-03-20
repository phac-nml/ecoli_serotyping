#!/usr/bin/env python

"""
    Serotype prediction for E. coli
"""


def predict_serotype(blast_record):
    """
    Entry point for serotype prediction
    :return: dict(
        otype = O
        htype = H
    )
    """

    return blast_record.query.strip()
