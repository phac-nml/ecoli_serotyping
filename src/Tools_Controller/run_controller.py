#!/usr/bin/env python

from flask_uploads import *
from tool_formatresults import *

import shutil
import Tkinter
import tkMessageBox
import subprocess
import ast
import logging

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
OUTPUT = ''
RESULTS = False
VERBOSITY = 0
PERC_ID = 90
PERC_LEN = 90
I = 0


app = Flask(__name__)
app.config['UPLOADED_FASTAFILES_DEST'] = SCRIPT_DIRECTORY + '../../temp/Uploads'
app.config['UPLOADED_FASTAFILES_ALLOW'] = set(['fasta', 'fsa_nt', 'fsa'])
files = UploadSet('fastafiles')
configure_uploads(app, (files,))



@app.route('/superphy/controller', methods =['POST', 'GET'])
def uploadFiles():

    if request.method == 'POST':
        global OUTPUT, I, RESULTS, VERBOSITY, PERC_ID, PERC_LEN


        PERC_ID = request.form['perc-id']
        PERC_LEN = request.form['perc-len']
        resultFiles = request.files.getlist('file')
        VERBOSITY = request.form['verbosity']
        RESULTS = request.form['table-radiobutton']
        threshold = request.form['threshold']

        serotyper_out = 0
        vf_out = 0
        amr_out = 0
        all_vfs = 0
        all_amr = 0

        filename = resultFiles[0].filename
        if len(resultFiles) == 1 and not filename:
            root = Tkinter.Tk()
            root.withdraw()
            tkMessageBox.showwarning('Oops!','No files were uploaded. Please try again.')

            return render_template('controller.html')


        if 'serotyper-checkbox' in request.form:
            serotyper_out = 1

        if 'vf-checkbox' in request.form:
            vf_out = 1

        if 'amr-checkbox' in request.form:
            amr_out = 1

        if 'all-vfs-checkbox' in request.form:
            all_vfs = 1

        if 'all-amr-checkbox' in request.form:
            all_amr = 1

        logging.info('In uploadFiles, method == POST')

        if os.path.isdir(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I)):
            shutil.rmtree(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))

        os.makedirs(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))
        logging.info('Directory')
        for file in resultFiles:
            files.save(file,'temp_dir'+ str(I))

        OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "tools_controller.py", "--input", SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I),
                                          "-s", str(serotyper_out), "-vf", str(vf_out), "-amr", str(amr_out), "-min", str(threshold), "-p", str(all_amr), "-avf", str(all_vfs),
                                              "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-sv', str(VERBOSITY), '-csv', '0'])
        I +=1
        return redirect(url_for('getResults'))
    return render_template('controller.html')



# @app.route('/curl-upload', methods=['POST'])
# def curl_uploadFiles():
#
#
#
#
@app.route('/superphy/controller/results', methods=['GET'])
def getResults():
    logging.info('In getResults')

    if 'Error' in OUTPUT:
        logging.info('No valid files uploaded')
        root = Tkinter.Tk()
        root.withdraw()
        tkMessageBox.showwarning('Oops!','No valid files were uploaded. Valid files are: .fasta, .fsa_nt.')
        return render_template('uploadfile.html')

    # elif RESULTS:
    #     logging.info('To HTML table')
    #     return ectyper_formatting.toHTML(ast.literal_eval(OUTPUT), VERBOSITY)
    else:
        logging.info('To JSON format')
        return toJSON(ast.literal_eval(OUTPUT))




@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
    app.run(debug=True)