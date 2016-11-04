#!/usr/bin/env python

from formatresults import *
from compareresults import *

TEST= False
SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def runProgram():

  logging.basicConfig(filename='ectyper.log',level=logging.INFO)
  if not os.path.isdir(SCRIPT_DIRECTORY + '../temp'):
    os.mkdir(SCRIPT_DIRECTORY + '../temp')
  args = parseCommandLine()

  roughGenomesList = getFilesList(args.input)
  genomesList = checkFiles(roughGenomesList)

  if isinstance(genomesList, list):
    if initializeDB() == 0:
      resultsList = runBlastQuery(genomesList, 'ECTyperDB')
      genomesDict = parseResults(resultsList)
      predictionDict = filterPredictions(genomesDict, args.percentIdentity, args.percentLength)
      matchDict = sortMatches(predictionDict)
      topMatches = findTopMatches(matchDict)
      logging.info('Top matches are ' + str(topMatches))

      if TEST == True:
        compareResults(roughGenomesList, topMatches)

      if args.csv == 'true':
        toCSV(topMatches, args.verbose)
      print formatResults(topMatches, args.verbose)

    else:
      logging.error('There was an error while generating the database.')
  else:
    print genomesList



if __name__=='__main__':
  runProgram()