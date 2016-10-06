#!/usr/bin/env python

from flask import *
from formatresults import *
from compareresults import *

TEST= False
GENOMES = {}

app = Flask(__name__)

@app.route('/ectyper/upload', methods =['POST'])
def uploadFiles():
  if 'file' not in request.files:
    return redirect(request.url)
  if len(request.files) == 0:
    abort(404)
  else:
    search =  subprocess.call([SCRIPT_DIRECTORY + "../src/ectyper.py", request.files])
    if search == 0:
      return redirect(url_for('runProgram'))


@app.route("/ectyper/results", methods=['GET'])
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
    GENOMES = findTopMatches(matchDict)
    logging.info('Top matches are ' + str(GENOMES))

    if TEST == True:
      compareResults(roughGenomesList, GENOMES)

    return formatResults(GENOMES, args.v)

  else:
    logging.error('There was an error while generating the database.')

@app.errorhandler(404)
def not_found(error):
  return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
 app.run(debug=True)