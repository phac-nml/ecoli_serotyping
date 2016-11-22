#!/usr/bin/env python

from formatresults import *
from compareresults import *

TEMP_SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
sys.path.append(os.path.abspath(TEMP_SCRIPT_DIRECTORY + '../'))

from createdirs import createDirs

TEST= False

def runProgram():

  logging.basicConfig(filename=SCRIPT_DIRECTORY + 'ectyper.log',level=logging.INFO)

  createDirs()

  args = parseCommandLine()

  if args.input == None:
    print 'Error'

  else:
    roughGenomesList = getFilesList(args.input)
    genomesList = checkFiles(roughGenomesList)

    if isinstance(genomesList, list):
      if initializeDB() == 0:
        results_file = runBlastQuery(genomesList, 'ECTyperDB')
        genomesDict = parseResults(results_file)
        predictionDict = filterPredictions(genomesDict, args.percentIdentity, args.percentLength)
        matchDict = sortMatches(predictionDict)
        topMatches = findTopMatches(matchDict)

        logging.info('Top matches are ' + str(topMatches))
        json_results = formatResults(topMatches, args.verbose)

        if TEST == True:
          compareResults(roughGenomesList, topMatches)

        if args.csv == 1:
          toCSV(json_results, args.verbose)

        print json_results


      else:
        logging.error('There was an error while generating the database.')
    else:
      print genomesList



if __name__=='__main__':
  runProgram()