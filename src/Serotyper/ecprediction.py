#!/usr/bin/env python
import re
import logging


def getProductPercentage(length, hsp):
        """
        Calculating the product of percentages (identity and length), to be used for comparisons when sorting the matches
        later on. If the product is higher than 1, it returns 0, as it's not permitted.

        :param length: Length of the alignment
        :param hsp: HSP
        :return product: Returns 0 if the product is greater than 1.
        """
        logging.info('Getting product percentage of \n' + str(hsp))
        if 0>= length:
            return 0

        hsp_length = abs(1-(1-float(hsp.positives)/length))
        hsp_identity = float(hsp.identities)/hsp.align_length

        product = hsp_identity*hsp_length

        if product>1:
            logging.warning("Product was higher than 1. Now set to 0.")
            return 0
        else:
            return product



def filterPredictions(predictionDict, percent_identity, percent_length):
    """
    Filtering the given dictionary with a loop to remove the matches that matched lower than the expected % length and/or the
    % identity. Instead of removing them from the initial dictionary, the method creates a new one and assigns the
    good matches.

    :param predictionDict: Unfiltered dictionary containing all the information about possible matches.
    :param percent_identity: Percentage of identity desired by the user.
    :param percent_length: Percentage of length desired by the user.
    :return newDict: Dictionary with matches that fill the percentage requirements.
    """

    logging.info("Filtering predictions with " + str(percent_identity) + "% identity and "+ str(percent_length) + "% length.")

    newDict = {}

    percent_identity = float(percent_identity)/100
    percent_length = float(percent_length)/100

    for genome_name, alignment in predictionDict.iteritems():
        logging.info(str(genome_name))
        tempDict = {}

        if isinstance(alignment, dict):
            for title, info in alignment.iteritems():
                al_length = info.keys()[0]
                hsp = info[al_length]
                hsp_length = float(hsp.positives)/int(al_length)
                hsp_identity = float(hsp.identities)/hsp.align_length

                #Comparing with the cutoffs
                if hsp_length >= percent_length:
                   if hsp_identity >= percent_identity:
                      tempDict[title] = [hsp, al_length]
                      newDict[genome_name] = tempDict
                   else:
                      tempDict[title] = [hsp, 'NA']
                      newDict[genome_name] = tempDict
                      logging.info("Percentage of identity was not high enough. \nAlignment: " + str(title))
                else:
                   tempDict[title] = [hsp, 'NA']
                   newDict[genome_name] = tempDict
                   logging.info("Percentage of length was not high enough. \nAlignment: " + str(title))

        else:
            newDict[genome_name] = 'NA'
            logging.info("There are no matches for genome " + str(genome_name))

    return newDict


def sortMatches(predictionDict):
     """
     Going through the prediction dictionary and finding the top possibles matches (can be perfect match or not).
     The return dictionary stores multiple types: fliC, wzx, wzy, wzm, wzt, fllA, flnaA, flkA, fmlA, gnd.

     :param predictionDict: Unsorted dictionary with all the possible matches.
     :return topMatchDict: Sorted dictionary with all the matches ordered by percentage.
     """


     logging.info("Sorting matches from " + str(predictionDict))
     sortedMatchDict = {}

     for genome_name, value in predictionDict.iteritems():
         tempDict = {'fliC':[],'flnaA':[], 'fllA':[], 'fmlA': [], 'flkA': [], 'gnd': [], 'wzx': [], 'wzy':[], 'wzm':[], 'wzt': []}

         if value == 'NA':
             sortedMatchDict[genome_name] = 'NA'

         else:
             for title, hsp in value.iteritems():
                 productPercentage = 0
                 if hsp[1] != 'NA':
                     productPercentage = getProductPercentage(hsp[1], hsp[0])


                 if re.search('fliC', title):
                    tempDict['fliC'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('wzx', title):
                    tempDict['wzx'].append({'title':title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('wzy', title):
                    tempDict['wzy'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('wzt', title):
                    tempDict['wzt'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('wzm', title):
                    tempDict['wzm'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('flnaA', title):
                    tempDict['flnaA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('fllA', title):
                    tempDict['fllA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('fmlA', title):
                    tempDict['fmlA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('flkA', title):
                    tempDict['flkA'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})

                 elif re.search('gnd', title):
                    tempDict['gnd'].append({'title': title, 'hsp': hsp[0], 'length': hsp[1], 'percentage': productPercentage})


            #Sorting the lists in the resulting dictionary
             for type in tempDict.keys():
                 if not tempDict[type]:
                     emptyList = tempDict.pop(type, None)
                 else:
                     sortedList = sorted(tempDict[type], key=lambda k: k['percentage'], reverse = True)
                     tempDict[type] = sortedList

             sortedMatchDict[genome_name] = tempDict

     return sortedMatchDict


def searchType(title, type):
    """
    Searching the name of the type in the title. If there is no match, the method returns None.

    :param title: String that will be searched to find the type.
    :param type: Usually O or H, to simplify the search.
    :return match: Result found when searching the string. If it is None, then the method returns the string 'none'.
    """
    match = re.search('(%s\d+$|-%s\d+|_%s\d+)' % (type, type, type), title)

    if match != None:
        match = str(match.group())
        if match.startswith('-'):
            match = match.split('-')[1]
        if match.startswith('_'):
            match = match.split('_')[1]
    else:
        match = re.search('OR$', title)

        if match != None:
            match = str(match.group())
        else:
            logging.warning("There are no matches in " + str(title) + " with the type " + str(type))
            match = 'none'

    return match


def findTopMatch(topMatch, matchDict, topMatchList, ohType):
    """
    Filtering the type dictionaries to find the best match to the user's query(ies). Help method to the findTopMatch method.
    Starts with an assigned top match, compares it to the other matches one by one(with percentage):
    1- if the top match's percentage is 0, then assign a new different match as the top match and go to the next match
    2- if a match's percentage is 0, go to the next match
    3- if the top match has a smaller percentage than a match, assign said match with the greater percentage as top match and go to
    the next match
    4- if the top match and another match have the same percentage but are not the same types, filter through their match
    lists to find similar matches, compare the percentages, and then step 3 is applied


    :param topMatch: Current top match to be compared with other matches.
    :param matchDict: Dictionary containing all the matches for a genome.
    :param topMatchList: List containing the current top match.
    :return topMatch: Top match for the genome.
    """


    for type, alignmentList in matchDict.iteritems():

        perc = alignmentList[0].get('percentage')
        topPerc = topMatch.get('percentage')

        if topPerc == 0:
            #Switch to next top match
            topMatch = alignmentList[0]
            topMatchList = alignmentList

            logging.info("Current top match (" + str(topMatch.get('title')) +  ") has percentage of 0. Setting " + str(alignmentList[0].get('title')) + " as the new top match.")
            continue

        elif perc == 0:
            #Stay with current match
            logging.info("The percentage of " + str(alignmentList[0].get('title')) + " is 0. Continuing with current top match (" + str(topMatch.get('title')) +".")
            continue

        elif perc > topPerc:
            #Switch to next top match with a higher percentage
            topMatch = alignmentList[0]
            topMatchList = alignmentList

            logging.info("Setting " + str(alignmentList[0].get('title')) + " as new top match.")
            continue

        elif perc == topPerc and str(alignmentList[0].get('title')) != str(topMatch.get('title')):
            #Search through the top match dictionary lists and the alignment list to obtain the best recurrences
            # for either the top match or the compared match

            item1 = 0
            item2 = 0
            match1 = searchType(alignmentList[0].get('title'), ohType)
            topMatch1 = searchType(topMatch.get('title'),ohType)

            for alignment in alignmentList:
                match2 = searchType(alignment.get('title'), ohType)
                if match2 == topMatch1:
                    item1 = alignment.get('percentage')
                    break

            for alignment in topMatchList:
                topMatch2 = searchType(alignment.get('title'), ohType)
                if topMatch2 == match1:
                    item2 == alignment.get('percentage')
                    break

            if item1>item2:
                #Switch to next top match with higher percentage from the best recurrence
                topMatch = alignmentList[0]
                topMatchList = alignmentList
                logging.info("Setting " + str(alignmentList[0].get('title')) + " as new top match.")

            continue

        else:
            continue

    return topMatch




def findTopMatches(matchDict):
    """
    Dividing the matches dictionary in two, O-types and H-types, to then find the best match.
    If the O-type or H-type dictionary is empty, its O-type or H-type will be NA.
    The top match dictionary contains the genomes names, their O-type, their H-type and their prediction strength.

    :param matchDict: Dictionary containing all matches for every genome.
    :return topMatchDict: Dictionary containing only the top matches for every genome.
    """

    topMatchDict={}

    for genome_name, typeDict in matchDict.iteritems():

        logging.info("Finding top matches for genome " + str(genome_name) + ".")
        tempDict = {}

        if typeDict == 'NA':
            topMatchDict[genome_name] = {'O type': 'NA', 'H type': 'NM', 'prediction_strength': 'NA'}

        else:

            oTypeDict = {
                k: typeDict[k] for k in typeDict.keys() if
                searchType(str(typeDict[k][0].get('title')), 'O').startswith('O')
            }

            if not oTypeDict:
                tempDict['O type'] = 'NA'
            else:
                #Find O type top match
                oKey = oTypeDict.keys()[0]
                oTypeList = oTypeDict[oKey]
                oTypeMatch = findTopMatch(oTypeList[0], oTypeDict, oTypeList, 'O')

                if oTypeMatch['length'] == 'NA':
                    tempDict['O type'] = 'NA'
                else:
                    tempDict['O type'] = oTypeMatch

            hTypeDict = {
                k: typeDict[k] for k in typeDict.keys() if
                searchType(str(typeDict[k][0].get('title')), 'H').startswith('H')
            }
            if not hTypeDict:
                tempDict['H type'] = 'NM'
            else:
                #Find H type top match
                hKey = hTypeDict.keys()[0]
                hTypeList = hTypeDict[hKey]
                hTypeMatch = findTopMatch(hTypeList[0], hTypeDict, hTypeList, 'H')
                if hTypeMatch['length'] == 'NA':
                    tempDict['H type'] = 'NM'
                else:
                    tempDict['H type'] = hTypeMatch

            tempDict['prediction_strength'] = 'Top match'
            topMatchDict[genome_name] = tempDict
    return topMatchDict