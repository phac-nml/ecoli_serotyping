#!/usr/bin/env python

from ecprediction import *
from ecvalidatingfiles import *


if __name__=='__main__':
  args = parseCommandLine()
  roughGenomesList = getListGenomes(args.input)
  genomesList = checkFiles(roughGenomesList)

  if initializeDB() == 0:
    resultsList = runBlastQuery(genomesList)
    alignmentsDict = parseResults(resultsList)
    identicalDict, predictionDict = findPerfectMatches(alignmentsDict)
    predictionDict = filterPredictions(predictionDict, args.pi, args.pl)

  else:
    print("Oops something happened")