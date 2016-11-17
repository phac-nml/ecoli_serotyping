#!/usr/bin/env python

from ectyper import *
from flask_uploads import *
from flask import *

import shutil
import Tkinter
import tkMessageBox
import subprocess
import ast
import formatresults

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
OUTPUT = ''
RESULTS = False
VERBOSITY = 'false'
PERC_ID = 90
PERC_LEN = 90
I = 0


app = Flask(__name__)
app.config['UPLOADED_FASTAFILES_DEST'] = SCRIPT_DIRECTORY + '../../temp/Uploads'
app.config['UPLOADED_FASTAFILES_ALLOW'] = set(['fasta', 'fsa_nt', 'fsa'])
files = UploadSet('fastafiles')
configure_uploads(app, (files,))

logging.basicConfig(filename=SCRIPT_DIRECTORY + 'ectyperdebug.log',level=logging.INFO)

@app.route('/upload', methods =['POST', 'GET'])
def uploadFiles():
    """
    Uploading files method. When the request is POST, the program runs the ectyper through the command line with the
    information provided by the user on the upload file form (files, percentages, result format).
    When the request is GET, the program returns the upload file form.

    :return: Results in table or JSON format or the upload a file form.
    """

    if request.method == 'POST':
       global OUTPUT, I, RESULTS, VERBOSITY, PERC_ID, PERC_LEN, IS_CURL


       PERC_ID = request.form['perc-id']
       PERC_LEN = request.form['perc-len']
       resultFiles = request.files.getlist('file')
       VERBOSITY = request.form['verbosity']

       logging.info('In uploadFiles, method == POST')

       if 'table-checkbox' in request.form:
            RESULTS = True
       else:
            RESULTS = False

       if len(resultFiles) == 1:
           filename = resultFiles[0].filename
           if not filename:
              root = Tkinter.Tk()
              root.withdraw()
              tkMessageBox.showwarning('Oops!','No files were uploaded. Please try again.')

              return render_template('uploadfile.html')

           else:
            files.save(resultFiles[0])
            logging.info('Single file')
            OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "--input", files.path(filename),
                                                "-pl", PERC_LEN, "-pi", PERC_ID, '-v', VERBOSITY, '-csv', 'false'])
       else:
           if os.path.isdir(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I)):
               shutil.rmtree(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))

           os.makedirs(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))
           logging.info('Directory')
           for file in resultFiles:
               files.save(file,'temp_dir'+ str(I))
           OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "--input", SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I),
                                                "-pl", PERC_LEN, "-pi", PERC_ID, '-v', VERBOSITY, '-csv', 'false'])
       I +=1
       return redirect(url_for('getResults'))
    return render_template('uploadfile.html')


@app.route('/curl-upload', methods=['POST'])
def curl_uploadFiles():
    """
    Uploading files method for the command line. When the request is POST, the program runs the ectyper through
    the command line with the information provided by the user in the curl command line.
    When the request is GET, the program returns a string.

    :return: Results in JSON format.
    """
    if request.method == 'POST':
        global OUTPUT, I, RESULTS, VERBOSITY, PERC_ID, PERC_LEN, IS_CURL

        resultFiles = request.files.getlist('files[]')
        if len(resultFiles) == 1:
            filename = resultFiles[0].filename
            if not filename:
                return 'No files were uploaded.'
            else:
                files.save(resultFiles[0])
                OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "-input", files.path(filename),
                                                  "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-v', VERBOSITY, '-csv', 'false'])
        else:
            if os.path.isdir(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I)):
                shutil.rmtree(SCRIPT_DIRECTORY + '../temp/Uploads/temp_dir' + str(I))

            os.makedirs(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))
            for file in resultFiles:
                files.save(file,'temp_dir'+ str(I))
            OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "-input",
                                              SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I),
                                              "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-v', VERBOSITY, '-csv', 'false'])
        I +=1
        return jsonify(ast.literal_eval(OUTPUT))
    return 'No HTTP requests were made.'


@app.route('/results', methods=['GET'])
def getResults():
    """
    Results formatting method. Used by the uploadFiles() method.

    :return: Results in the format desired by the user (table or JSON).
    """
    logging.info('In getResults')

    if 'Error' in OUTPUT:
        logging.info('No valid files uploaded')
        root = Tkinter.Tk()
        root.withdraw()
        tkMessageBox.showwarning('Oops!','No valid files were uploaded. Valid files are: .fasta, .fsa_nt.')
        return render_template('uploadfile.html')

    elif RESULTS:
        logging.info('To HTML table')
        return formatresults.toHTML(ast.literal_eval(OUTPUT), VERBOSITY)
    else:
        logging.info('To JSON format')
        return jsonify(ast.literal_eval(OUTPUT))


@app.errorhandler(404)
def not_found(error):
  return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
    app.run(debug=True)