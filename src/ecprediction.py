#!/usr/bin/env python

def findPerfectMatches(genomesDict):
    """
    Identifying the identical matches and the ones that will need a prediction.
    The method stores them in two dictionaries: one for identical matches and one that will need a prediction.

    :param genomesDict:
    :return identicalDict, predictionDict:
    """
    identicalDict = {}
    predictionDict = {}

    for genome_name, value in genomesDict.iteritems():
        tempPDict = {}
        tempIDict = {}
        if isinstance(value, dict):
            for title, hsp in value.iteritems():
                if len(hsp.query) == hsp.positives:
                    tempIDict[title] = hsp
                    identicalDict[genome_name] = tempIDict
                else:
                    tempPDict[title] = hsp
                    predictionDict[genome_name]= tempPDict

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

    for genome_name, value in predictionDict.iteritems():
        for title, hsp in value.iteritems():
         hsp_length = float(hsp.align_length)/len(hsp.query)
         hsp_identity = float(hsp.positives)/hsp.align_length

         if hsp_length >= percent_length:
           if hsp_identity >= percent_identity:
               newDict[title] = hsp

    return newDict

