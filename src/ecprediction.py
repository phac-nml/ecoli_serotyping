#!/usr/bin/env python
import re


# def findPerfectMatches(genomesDict):
#     """
#     Identifying the identical matches and the ones that will need a prediction.
#     The method stores them in two dictionaries: one for identical matches and one that will need a prediction.
#
#     :param genomesDict:
#     :return identicalDict, predictionDict:
#     """
#     identicalDict = {}
#     predictionDict = {}
#
#     for genome_name, value in genomesDict.iteritems():
#         tempPDict = {}
#         tempIDict = {}
#         if isinstance(value, dict):
#             for title, hsp in value.iteritems():
#                 if len(hsp.query) == hsp.positives:
#                     print title
#                     tempIDict[title] = hsp
#                     identicalDict[genome_name] = tempIDict
#                 else:
#                     tempPDict[title] = hsp
#                     predictionDict[genome_name]= tempPDict
#
#     return identicalDict, predictionDict


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
        tempDict = {}
        if isinstance(value, dict):
            for title, hsp in value.iteritems():
             hsp_length = float(hsp.align_length)/len(hsp.query)
             hsp_identity = float(hsp.positives)/hsp.align_length
             if hsp_length >= percent_length:
                if hsp_identity >= percent_identity:
                   tempDict[title] = [hsp, hsp_length*hsp_identity]
                   newDict[genome_name] = tempDict

    return newDict


def findTopMatch(predictionDict):
 topMatchDict = {}

 for genome_name, value in predictionDict.iteritems():
     tempDict = {}
     flic_perc = 0
     flnaa_perc = 0
     flla_perc = 0
     fmla_perc = 0
     flka_perc = 0
     gnd_perc = 0
     wzx_perc = 0
     wzy_perc = 0
     wzm_perc = 0
     wzt_perc = 0

     for title, hsp in value.iteritems():

        if re.search('(fliC-H\d)', title):
            if hsp[1] > flic_perc:
                tempDict['fliC'] = {title: hsp[0]}

        elif re.search('(wzx-O\d)', title):
            if hsp[1] > wzx_perc:
                tempDict['wzx'] = {title:hsp[0]}

        elif re.search('(wzy-O\d)', title):
            if hsp[1] > wzy_perc:
                tempDict['wzy'] = {title:hsp[0]}

        elif re.search('(wzt-O\d)', title):
            if hsp[1] > wzt_perc:
                tempDict['wzt'] = {title:hsp[0]}

        elif re.search('(wzm-O\d)', title):
            if hsp[1] > wzm_perc:
                tempDict['wzm'] = {title:hsp[0]}

        elif re.search('(flnaA-H\d)', title):
            if hsp[1] > flnaa_perc:
                tempDict['flnaA'] = {title:hsp[0]}

        elif re.search('(fllA-H\d)', title):
            if hsp[1] > flla_perc:
                tempDict['fllA'] = {title:hsp[0]}

        elif re.search('(fmlA-H\d)', title):
            if hsp[1] > fmla_perc:
                tempDict['fmlA'] = {title:hsp[0]}

        elif re.search('(flkA-H\d)', title):
            if hsp[1] > flka_perc:
                tempDict['flkA'] = {title:hsp[0]}

        elif re.search('(gnd-O\d)', title):
            if hsp[1] > gnd_perc:
                tempDict['gnd'] = {title:hsp[0]}

     topMatchDict[genome_name] = tempDict

 print topMatchDict
 return topMatchDict



