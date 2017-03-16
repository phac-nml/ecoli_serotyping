#!/usr/bin/env python

"""
    Predictive genomics for _E. coli_ from the command line.
    Currently includes serotyping and VF finding.
"""

import logging.config
import definitions

logging.config.fileConfig(definitions.LOGGER_CONFIG)
log = logging.getLogger(__name__)

# TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
# sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))
#
# from sharedmethods import *
#
# SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
# TEST= False
#  logging.basicConfig(filename=SCRIPT_DIRECTORY + 'ectyper.log',level=logging.INFO)
#  logging.info('Starting E. Coli Serotyper')


def runProgram():
    """
    Wrapper for both the serotyping and virulence finder
    :return: success or failure
    """

    # createDirs()
    #args = parseCommandLine()

    # if args.input == None:
    #     logging.info('No inputs were given.')
    #     print 'Error'
    #
    # else:
    #     roughGenomesList = getFilesList(args.input)
    #     genomesList = checkFiles(roughGenomesList)
    #
    #     if isinstance(genomesList, list):
    #         setGlobalDicts()
    #         clearGlobalDicts()
    #
    #         if initializeDB() == 0:
    #             results_file = runBlastQuery(genomesList, 'ECTyperDB')
    #             genomesDict = parseResults(results_file)
    #             predictionDict = filterPredictions(genomesDict,
    #                                                args.percentIdentity,
    #                                                args.percentLength)
    #             matchDict = sortMatches(predictionDict)
    #             topMatches = findTopMatches(matchDict)
    #
    #             logging.info('Top matches are ' + str(topMatches))
    #             json_results = formatResults(topMatches, args.verbose)
    #
    #             if TEST == True:
    #                 compareResults(roughGenomesList, topMatches)
    #
    #             if args.galaxy == 1:
    #                 if args.outputFormat == 'json':
    #                     json.dump(json_results, args.output)
    #                 else:
    #                     toGalaxyCSV(json_results, args.verbose, args.output)
    #             else:
    #                 if args.csv == 1:
    #                     toCSV(json_results, args.verbose)
    #
    #                 print json_results
    #                 logging.info('Successfully ended the program.')
    #
    #
    #         else:
    #             logging.error(
    #                 'There was an error while generating the database.')
    #             print 'Error'
    #     else:
    #         print genomesList
    log.info("Done")