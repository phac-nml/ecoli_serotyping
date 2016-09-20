#!/usr/bin/env python

from ecprediction import *
from ecvalidatingfiles import *


if __name__=='__main__':
  args = parseCommandLine()
  roughGenomesList = getListGenomes(args.input)
  genomesList = checkFiles(roughGenomesList)

  if initializeDB() == 0:
    resultsList = runBlastQuery(genomesList)
    genomesDict = parseResults(resultsList)
    #identicalDict, predictionDict = findPerfectMatches(genomesDict)
    predictionDict = filterPredictions(genomesDict, args.pi, args.pl)
    topMatchDict = findTopMatch(predictionDict)

  else:
    print("Oops something happened")