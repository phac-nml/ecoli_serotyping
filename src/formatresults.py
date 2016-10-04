#!/usr/bin/env python
import csv
import os.path

from flask import jsonify
from ecprediction import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def toResultDict(topMatches, verbose):
    resultDict = {}
    if verbose == 'true':
        for genome_name, serotype in topMatches.iteritems():
            if serotype == 'NA':
                resultDict[genome_name] = 'No matches were found for this genome, thus no prediction could be made.'
            else:
                tempDict1 = {}
                for type, info in serotype.iteritems():
                    if info == 'NA' and type != 'prediction_strength':
                        tempDict1[type] = 'No prediction could be made for O types for this genome.'
                    elif info == 'NM':
                        tempDict1[type] = 'No prediction could be made for H types for this genome (non-motile).'
                    elif type != 'prediction_strength':
                        tempDict2 = {}
                        for key, value in info.iteritems():
                            if key == 'title':
                                otype = searchType(value, 'O')
                                htype = searchType(value, 'H')

                                if otype != 'none':
                                    tempDict2['RESULT'] =  str(otype)
                                else:
                                    tempDict2['RESULT'] = str(htype)
                            elif key != 'hsp':
                                tempDict2[key] = str(value)
                            else:
                                tempDict2['score'] = str(value.score)
                                tempDict2['bits'] = str(value.bits)
                                tempDict2['e value'] = str(value.expect)
                        tempDict1[type] = tempDict2
                    else:
                        tempDict1['prediction strength'] = str(info)
                    resultDict[genome_name] = tempDict1
    else:
        for genome_name, serotype in topMatches.iteritems():
            if serotype == 'NA':
                resultDict[genome_name] = 'No matches were found for this genome, thus no prediction could be made.'
            else:
                tempDict = {}
                if topMatches[genome_name]['otype'] != 'NA':
                    tempDict['O type'] = searchType(topMatches[genome_name]['otype']['title'], 'O')
                else:
                    tempDict['O type'] = 'No prediction could be made for O types for this genome.'

                if topMatches[genome_name]['htype'] != 'NM':
                    tempDict['H type'] = searchType(topMatches[genome_name]['htype']['title'], 'H')
                else:
                    tempDict['H type'] = 'No prediction could be made for H types for this genome (non-motile).'

                resultDict[genome_name] = tempDict

    return resultDict


def toCSV(genomes_parsed, verbose):
    header = ['Genome', 'O Type', 'H Type']

    serotype_file = open(SCRIPT_DIRECTORY + '../temp/Results/Serotype_Results.csv', 'wb')
    csvwriter = csv.DictWriter(serotype_file, header)
    csvwriter.writeheader()

    for genome_name, type in genomes_parsed.iteritems():
        row = {'Genome': genome_name}

        if verbose == 'false':
            if isinstance(type, dict):
                keys = sorted(type.keys())
                row.update({'H Type': genomes_parsed[genome_name][keys[0]], 'O Type': genomes_parsed[genome_name][keys[1]]})
            else:
                row.update({'H Type': genomes_parsed[genome_name], 'O Type': genomes_parsed[genome_name]})
        elif isinstance(type, dict):
            oTempStr =''
            hTempStr = ''

            if isinstance(genomes_parsed[genome_name]['htype'], dict):
                for title, info in genomes_parsed[genome_name]['htype'].iteritems():
                    hTempStr+= str(title) + ": " + str(info) + "\n"
                row.update({'H Type': hTempStr})
            else:
                row.update({'H Type': genomes_parsed[genome_name]['htype']})

            if isinstance(genomes_parsed[genome_name]['otype'], dict):
                for title, info in genomes_parsed[genome_name]['otype'].iteritems():
                    oTempStr += str(title) + ": " + str(info) + "\n"
                row.update({'O Type': oTempStr})
            else:
                row.update({'O Type': genomes_parsed[genome_name]['otype']})

        else:
            row.update({'H Type': genomes_parsed[genome_name], 'O Type': genomes_parsed[genome_name]})
        csvwriter.writerow(row)

    serotype_file.close()



def formatResults(topMatches, verbose):

    json_data =  toResultDict(topMatches,verbose)
    toCSV(json_data, verbose)

    return jsonify(json_data)
