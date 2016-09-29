#!/usr/bin/env python

from compareresults import *

TEST= True

if __name__=='__main__':
  logging.basicConfig(filename='ectyper.log',level=logging.INFO)

  args = parseCommandLine()
  roughGenomesList = getFilesList(args.input)
  genomesList = checkFiles(roughGenomesList)

  if initializeDB() == 0:
    resultsList = runBlastQuery(genomesList, 'ECTyperDB')
    genomesDict = parseResults(resultsList)
    predictionDict = filterPredictions(genomesDict, args.pi, args.pl)
    matchDict = sortMatches(predictionDict)
    topMatches = findTopMatches(matchDict)
    logging.info('Top matches are ' + str(topMatches))

    if TEST == True:
      compareResults(roughGenomesList, topMatches)

  else:
    logging.error('There was an error while generating the database.')