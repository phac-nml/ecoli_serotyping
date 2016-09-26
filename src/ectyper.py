#!/usr/bin/env python

import logging

from ecprediction import *
from ecvalidatingfiles import *


if __name__=='__main__':
  logging.basicConfig(filename='ectyper.log',level=logging.INFO)

  args = parseCommandLine()
  roughGenomesList = getFilesList(args.input)
  genomesList = checkFiles(roughGenomesList)

  if initializeDB() == 0:
    resultsList = runBlastQuery(genomesList)
    genomesDict = parseResults(resultsList)
    predictionDict = filterPredictions(genomesDict, args.pi, args.pl)
    matchDict = sortMatches(predictionDict)
    topMatches = findTopMatches(matchDict)
    logging.info('Top matches are ' + str(topMatches))

  else:
    logging.error('There was an error while generating the database.')