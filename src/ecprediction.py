#!/usr/bin/env python
import re


def getProductPercentage(length, hsp):
    hsp_length = float(hsp.positives)/length
    hsp_identity = float(hsp.identities)/hsp.align_length
    return hsp_identity*hsp_length


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

    for genome_name, alignment in predictionDict.iteritems():
        tempDict = {}
        if isinstance(alignment, dict):
            for title, info in alignment.iteritems():
                al_length = info.keys()[0]
                hsp = info[al_length]
                hsp_length = float(hsp.positives)/int(al_length)
                hsp_identity = float(hsp.identities)/hsp.align_length
                if hsp_length >= percent_length:
                   if hsp_identity >= percent_identity:
                      tempDict[title] = [hsp, al_length]
                      newDict[genome_name] = tempDict
                   else:
                      tempDict[title] = [hsp, 'NA']
                      newDict[genome_name] = tempDict
                else:
                   tempDict[title] = [hsp, 'NA']
                   newDict[genome_name] = tempDict

        elif not alignment:
            newDict[genome_name] = 'NA'

    return newDict


def sortMatches(predictionDict):
     """
     Going through the prediction dictionary and finding the top possibles matches (can be perfect match or not).
     The return dictionary stores multiple types: fliC, wzx, wzy, wzm, wzt, fllA, flnaA, flkA, fmlA, gnd.

     :param predictionDict:
     :return topMatchDict:
     """
     print "Sorting matches..."
     sortedMatchDict = {}

     for genome_name, value in predictionDict.iteritems():
         tempDict = {'fliC':[],'flnaA':[], 'fllA':[], 'fmlA': [], 'flkA': [], 'gnd': [], 'wzx': [], 'wzy':[], 'wzm':[], 'wzt': []}

         if value == 'NA':
             sortedMatchDict[genome_name] = 'NA'

         else:
             for title, hsp in value.iteritems():
                 productPercentage = 'NA'
                 if hsp[1] != 'NA':
                     productPercentage = getProductPercentage(hsp[1], hsp[0])


                 if re.search('(fliC-H\d)', title):
                    tempDict['fliC'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(wzx-O\d)', title):
                    tempDict['wzx'].append({'title':title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(wzy-O\d)', title):
                    tempDict['wzy'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(wzt-O\d)', title):
                    tempDict['wzt'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(wzm-O\d)', title):
                    tempDict['wzm'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(flnaA-H\d)', title):
                    tempDict['flnaA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(fllA-H\d)', title):
                    tempDict['fllA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(fmlA-H\d)', title):
                    tempDict['fmlA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(flkA-H\d)', title):
                    tempDict['flkA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})

                 elif re.search('(gnd-O\d)', title):
                    tempDict['gnd'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'perc': productPercentage})



             for type in tempDict.keys():
                 if not tempDict[type]:
                     emptyList = tempDict.pop(type, None)
                 else:
                     sortedList = sorted(tempDict[type], key=lambda k: k['perc'], reverse = True)
                     tempDict[type] = sortedList

             sortedMatchDict[genome_name] = tempDict

     return sortedMatchDict


def searchType(title, type):
    match = re.search('(%s-[a-zA-Z]\d+)' % type, title)
    match = str(match.group())
    match = match.split('-')[1]
    return match


def findType(topMatch, matchDict, topMatchList):


    for type, alignmentList in matchDict.iteritems():
        perc = alignmentList[0].get('perc')
        topPerc = topMatch.get('perc')

        if topPerc == 'NA':
            topMatch = alignmentList[0]
            topMatchList = alignmentList

        elif perc == 'NA':
            continue
        elif perc > topPerc:
            topMatch = alignmentList[0]
            topMatchList = alignmentList
        elif perc == topPerc and str(alignmentList[0].get('title')) != str(topMatch.get('title')):
            item1 = 0
            item2 = 0
            match1 = searchType(alignmentList[0].get('title'), type)
            topMatch1 = searchType(topMatch.get('title'),'[a-zA-Z]')

            for alignment in alignmentList:
                match2 = searchType(alignment.get('title'), type)
                if match2 == topMatch1:
                    item1 = alignment.get('perc')
                    break

            for alignment in topMatchList:
                topMatch2 = searchType(alignment.get('title'), '[a-zA-Z]')
                if topMatch2 == match1:
                    item2 == alignment.get('perc')
                    break

            if item1>item2:
                topMatch = alignmentList[0]
                topMatchList = alignmentList
        else:
            continue


    return topMatch




def findTopMatch(matchDict):

    topMatchDict={}

    print "Finding O and H types..."

    for genome_name, typeDict in matchDict.iteritems():
        tempDict = {}

        if typeDict == 'NA':
            topMatchDict[genome_name] = {'otype': 'NA', 'htype': 'NA', 'predictionstrength': 'tbd'}

        else:

            oTypeDict = {
                k: typeDict[k] for k in typeDict.keys() if
                searchType(str(typeDict[k][0].get('title')), k).startswith('O')
            }

            hTypeDict = {
                k: typeDict[k] for k in typeDict.keys() if
                searchType(str(typeDict[k][0].get('title')), k).startswith('H')
            }



            oKey = oTypeDict.keys()[0]
            oTypeList = oTypeDict[oKey]
            oTypeMatch = findType(oTypeList[0], oTypeDict, oTypeList)

            if oTypeMatch['perc'] == 'NA':
                tempDict['otype'] = 'NA'
            else:
                tempDict['otype'] = oTypeMatch

            hKey = hTypeDict.keys()[0]
            hTypeList = hTypeDict[hKey]
            hTypeMatch = findType(hTypeList[0], hTypeDict, hTypeList)

            if hTypeMatch['perc'] == 'NA':
                tempDict['htype'] = 'NA'
            else:
                tempDict['htype'] = hTypeMatch

            topMatchDict[genome_name] = tempDict

    return topMatchDict