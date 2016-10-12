#!/usr/bin/env python
import csv
import os.path

from flask import url_for
from ecprediction import *

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"


def toHTML():
    with open(SCRIPT_DIRECTORY + '../temp/Results/Serotype_Results.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        stylesheet_url = url_for('static',filename='style.css')
        html_info = "<!DOCTYPE html> \
                     <html lang='en'> \
                     <head> \
                        <meta charset='UTF-8'> \
                        <title>Results</title> \
                        <link rel='stylesheet' type='text/css' href=" + str(stylesheet_url) + ">  \
                     </head>\
                     <body>"
        result_table = "<div id='results-div'><table class=results><caption class='results-caption'>SEROTYPE PREDICTION RESULTS" \
                       "</caption>" \
                       "<thead class='results-head'><tr>" \
                            "<th>Genome</th>" \
                            "<th>O Type</th>" \
                            "<th>H Type</th>" \
                       "</tr></thead><tbody>"

        for line in reader:
            result_genome = str(line['Genome']).replace(";", '<br>')
            result_otype = str(line['O Type']).replace(";", '<br>')
            result_htype = str(line['H Type']).replace(";", '<br>')

            result_table+= "<tr><td id='result-genome'>" + result_genome + "</td><td>" + result_otype + "</td><td>" + \
                           result_htype + "</td></tr>"

        result_table+= "</tbody></table></div>"

        return_button = "<div id='return-button'>" \
                        "<input type='button' class='button return-button' onclick='" + \
                        'location.href="http://127.0.0.1:5000/ectyper/upload";' + \
                        "' value='Return to main page'/></div>"

        return html_info + result_table + return_button + "</body></html>"

def toResultDict(topMatches, verbose):
    """
    Filtering through the top matches dictionary to only keep the necessary result information (may be changed later
    on).

    :param topMatches: Top matches dictionary containing all the unfiltered results.
    :param verbose: Boolean stating whether the user wants full information or not.
    :return resultDict: Top matches dictionary containing the filtered results.
    """

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
    """
    Writing final results in a CSV file.

    :param genomes_parsed: Top matches dictionary with the final results.
    :param verbose: Boolean stating whether the user wants full information or not.
    """

    header = ['Genome', 'O Type', 'H Type']

    with open(SCRIPT_DIRECTORY + '../temp/Results/Serotype_Results.csv', 'wb') as csvfile:
        csvwriter = csv.DictWriter(csvfile, header)
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
                        if title == 'RESULT':
                            hTempStr= str(title) + ": " + str(info) + "; \n" + hTempStr
                        else:
                            hTempStr+= str(title) + ": " + str(info) + "; \n"
                    row.update({'H Type': hTempStr})
                else:
                    row.update({'H Type': genomes_parsed[genome_name]['htype']})

                if isinstance(genomes_parsed[genome_name]['otype'], dict):
                    for title, info in genomes_parsed[genome_name]['otype'].iteritems():
                        if title == 'RESULT':
                            oTempStr = str(title) + ": " + str(info) + "; \n" + oTempStr
                        else:
                            oTempStr += str(title) + ": " + str(info) + "; \n"
                    row.update({'O Type': oTempStr})
                else:
                    row.update({'O Type': genomes_parsed[genome_name]['otype']})

            else:
                row.update({'H Type': genomes_parsed[genome_name], 'O Type': genomes_parsed[genome_name]})
            csvwriter.writerow(row)



def formatResults(topMatches, verbose):

    json_data =  toResultDict(topMatches,verbose)
    toCSV(json_data, verbose)

    return json_data
