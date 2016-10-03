#!/usr/bin/env python

from flask import Flask
from flask import jsonify
from compareresults import *

TEST= False


app = Flask(__name__)



def formatResults(topMatches, verbose):

  resultDict = {}
  if verbose == 'true':
    for genome_name, serotype in topMatches.iteritems():
      if serotype == 'NA':
        resultDict[genome_name] = 'No matches were found for this genome, thus no prediction could be made.'
      else:
        tempDict1 = {}
        for type, info in serotype.iteritems():
          if info == 'NA' and type != 'prediction_strength':
            tempDict1[type] = 'No matches were found for O types for this genome, thus no prediction could be made.'
          elif info == 'NM':
            tempDict1[type] = 'No matches were found for H types for this genome, thus no prediction could be made (non-motile).'
          elif type != 'prediction_strength':
            tempDict2 = {}
            for key, value in info.iteritems():
              if key == 'title':
                otype = searchType(value, 'O')
                htype = searchType(value, 'H')

                if otype != 'none':
                  tempDict2['O type'] =  str(otype)
                else:
                  tempDict2['H type'] = str(htype)
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
          tempDict['O type'] = 'No matches were found for O types for this genome, thus no prediction could be made.'

        if topMatches[genome_name]['htype'] != 'NM':
          tempDict['H type'] = searchType(topMatches[genome_name]['htype']['title'], 'H')
        else:
          tempDict['H type'] = 'No matches were found for H types for this genome, thus no prediction could be made (non-motile).'

        resultDict[genome_name] = tempDict

  return jsonify(resultDict)


@app.route("/")
def runProgram():
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

    return formatResults(topMatches, args.v)

  else:
    logging.error('There was an error while generating the database.')

if __name__=='__main__':
 app.run(debug=True)