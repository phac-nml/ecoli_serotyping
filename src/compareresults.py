#!/usr/bin/env python

from ecvalidatingfiles import *
from ecprediction import *

import csv

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def initializeDanishDB():
    """
    Generating the Serotype Finder database to later query it.

    :return 0 if the database was created successfully, 1 otherwise:
    """

    REL_DIR = SCRIPT_DIRECTORY + '../temp/'

    return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/serotypefinder/O_H_type.fsa ", "-dbtype", "nucl", "-title", "SerotypeFinderDB", "-out", REL_DIR + "SerotypeFinderDB"])



def writeFile(topMatches, results):
    """
    Writing a table containing all the results (top matches) from both the Serotype Finder database and the EcOH database.
    The method compares them, and add comments if they have different results.

    :param topMatches:
    :param results:
    """

    with open(SCRIPT_DIRECTORY + "../temp/results.csv", 'wb') as csvfile:
        header = ['Genomes','Type', 'Serotype Finder', 'E.C. Typer', 'Notes']
        w = csv.DictWriter(csvfile, header)
        w.writeheader()

        for genome_name, serotype in sorted(topMatches.items()):
            keys = sorted(serotype.keys())

            row1 = {'Genomes': genome_name}
            row1.update({'Type':keys[0]})

            if serotype[keys[0]] == 'NM':
                row1.update({'Serotype Finder': 'No prediction could be made (non-motile).'})
            if results[genome_name][keys[0]] == 'NM':
                row1.update({'E.C. Typer': 'No prediction could be made (non-motile).'})

            else:
                tempDict1 = {}
                tempDict2 = {}
                for k in serotype[keys[0]]:
                    if k != 'hsp':
                        if k == 'title':
                            match1 = searchType(serotype[keys[0]][k], 'H')
                            match2 = searchType(results[genome_name][keys[0]][k], 'H')
                            tempDict1['result'] = (match1)
                            tempDict2['result'] = (match2)

                            if match1 != match2:
                                row1.update({'Notes': "The Serotype Finder and the E.C Typer do not share the same results for the " + str(keys[0]) + " in this genome."})
                        else:
                            tempDict1[k] = serotype[keys[0]][k]
                            tempDict2[k] = results[genome_name][keys[0]][k]
                            if tempDict1[k] != tempDict2[k]:
                                row1.update({'Notes': "The Serotype Finder and the E.C. Typer do not share the same " + str(k) +" for the " + str(keys[0]) + " in this genome."})
                row1.update({'Serotype Finder': tempDict1, 'E.C. Typer': tempDict2})

            w.writerow(row1)

            row2 = {'Genomes': ''}
            row2.update({'Type': keys[1]})

            if serotype[keys[1]] == 'NA':
                row2.update({'Serotype Finder': 'No prediction could be made'})
            if results[genome_name][keys[1]] == 'NA':
                row2.update({'E.C. Typer': 'No prediction could be made'})
            else:
                tempDict1 = {}
                tempDict2 = {}
                for k in serotype[keys[1]]:
                    if k != 'hsp':
                        if k == 'title':
                            match1 = searchType(serotype[keys[1]][k], 'O')
                            match2 = searchType(results[genome_name][keys[1]][k], 'O')
                            tempDict1['result'] = (match1)
                            tempDict2['result'] = (match2)

                            if match1 != match2:
                                row2.update({'Notes': "The Serotype Finder and the E.C Typer do not share the same results for the " + str(keys[1]) + " in this genome."})
                        else:
                            tempDict1[k] = serotype[keys[1]][k]
                            tempDict2[k] = results[genome_name][keys[1]][k]
                            if tempDict1[k] != tempDict2[k]:
                                row2.update({'Notes': "The Serotype Finder and the E.C. Typer do not share the same " + str(k) +" for the " + str(keys[1]) + " in this genome."})
                row2.update({'Serotype Finder': tempDict1, 'E.C. Typer': tempDict2})
            w.writerow(row2)

            row3 = {'Genomes': ''}
            row3.update({'Type': keys[2], 'Serotype Finder': serotype[keys[2]], 'E.C. Typer': results[genome_name][keys[2]] })
            w.writerow(row3)


def compareResults(filesList, results):

    if initializeDanishDB() == 0:
        resultsList = runBlastQuery(filesList, 'SerotypeFinderDB')
        genomesDict = parseResults(resultsList)
        predictionDict = filterPredictions(genomesDict, 85, 60)
        matchDict = sortMatches(predictionDict)
        topMatches = findTopMatches(matchDict)
        writeFile(topMatches, results)






