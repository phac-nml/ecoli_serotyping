#!/usr/bin/env python

from Bio.Blast import NCBIXML


def filterPredictions(predictionDict, percent_identity, percent_length):

    newDict = {}
    for title, hsp in predictionDict.iteritem():
        hsp_length = hsp.align_length/len(hsp.query)
        hsp_identity = hsp.positives/hsp.align_length

        percent_identity = percent_identity/100
        percent_length = percent_length/100

        if hsp_length >= percent_length:
            if hsp_identity >= percent_identity:
                newDict[title] = hsp

    print newDict
    return newDict
