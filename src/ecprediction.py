#!/usr/bin/env python

def findPerfectMatches(alignmentsDict):
    """
    Identifying the identical matches and the ones that will need a prediction.
    The method stores them in two dictionaries: one for identical matches and one that will need a prediction.

    :param alignmentsDict:
    :return identicalDict, predictionDict:
    """
    identicalDict = {}
    predictionDict = {}

    for title,hsp in alignmentsDict.iteritems():
        if len(hsp.query) == hsp.positives:
            identicalDict[title] = hsp
        else:
            predictionDict[title] = hsp

    return identicalDict, predictionDict


def filterPredictions(predictionDict, percent_identity, percent_length):
    """
    Filtering the given dictionary with a loop to remove the matches that matched lower than the expected % length and/or the
    % identity. Instead of removing them from the initial dictionary, the method creates a new one and assigns the
    good matches.

    :param predictionDict:
    :param percent_identity:
    :param percent_length:
    :return newDict:
    """
    newDict = {}

    print "Filtering..."

    percent_identity = float(percent_identity)/100
    percent_length = float(percent_length)/100

    for title, hsp in predictionDict.iteritems():

        hsp_length = float(hsp.align_length)/len(hsp.query)
        hsp_identity = float(hsp.positives)/hsp.align_length

        if hsp_length >= percent_length:
            if hsp_identity >= percent_identity:
                newDict[title] = hsp

    return newDict

