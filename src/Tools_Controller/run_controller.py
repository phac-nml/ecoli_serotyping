#!/usr/bin/env python

from flask_uploads import *
from flask import *
from tools_formatresults import *

import shutil
import Tkinter
import tkMessageBox
import subprocess
import ast
import logging

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
OUTPUT = ''
FORMAT = 0
VERBOSITY = 0
PERC_ID = 90
PERC_LEN = 90
I = 0


app = Flask(__name__)
app.config['UPLOADED_FASTAFILES_DEST'] = SCRIPT_DIRECTORY + '../../temp/Uploads'
app.config['UPLOADED_FASTAFILES_ALLOW'] = set(['fasta', 'fsa_nt', 'fsa', 'fna'])
files = UploadSet('fastafiles')
configure_uploads(app, (files,))

logging.basicConfig(filename=SCRIPT_DIRECTORY + 'controllerdebug.log',level=logging.INFO)

@app.route('/ectyper', methods =['POST', 'GET'])
def uploadFiles():

    if request.method == 'POST':
        global OUTPUT, I, FORMAT, VERBOSITY, PERC_ID, PERC_LEN


        PERC_ID = request.form['perc-id']
        PERC_LEN = request.form['perc-len']
        resultFiles = request.files.getlist('files[]')
        VERBOSITY = request.form['verbosity']
        FORMAT = request.form['table-radiobutton']
        threshold = request.form['threshold']

        serotyper_out = 0
        vf_out = 0
        amr_out = 0
        all_vfs = 0
        perf_amr = 1

        filename = resultFiles[0].filename
        if len(resultFiles) == 1 and not filename:
            root = Tkinter.Tk()
            root.withdraw()
            tkMessageBox.showwarning('Oops!','No files were uploaded. Please try again.')

            return render_template('controller.html')

        if 'serotyper-checkbox' not in request.form and\
            'vf-checkbox' not in request.form and \
            'amr-checkbox' not in request.form:

            root = Tkinter.Tk()
            root.withdraw()
            tkMessageBox.showwarning('Oops!','No tools from Serotyper, Virulence Factors or AMR were selected. Please try again.')

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
            perf_amr = 0


        logging.info('In uploadFiles, method == POST')

        if os.path.isdir(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I)):
            shutil.rmtree(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))

        os.makedirs(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))
        logging.info('Starting controller')
        for file in resultFiles:
            files.save(file,'temp_dir'+ str(I))

        OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "tools_controller.py", "--input", SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I),
                                          "-s", str(serotyper_out), "-vf", str(vf_out), "-amr", str(amr_out), "-min", str(threshold), "-p", str(perf_amr), "-avf", str(all_vfs),
                                              "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-sv', str(VERBOSITY), '-csv', '0'])
        OUTPUT = OUTPUT.split('\n')[0]

        I +=1

        if all_vfs == 1 or perf_amr == 0:
            return redirect(url_for('straightDownload'))

        return redirect(url_for('getResults'))
    return render_template('controller.html')



@app.route('/ectyper/curl-upload', methods=['POST'])
def curl_uploadFiles():
    """
    Uploading files method for the command line. When the request is POST, the program runs the ectyper through
    the command line with the information provided by the user in the curl command line.
    When the request is GET, the program returns a string.

    :return: Results in JSON format.
    """
    if request.method == 'POST':
        global OUTPUT, I, VERBOSITY, PERC_ID, PERC_LEN
        s = 0
        vf = 0
        amr = 0
        avf = 0
        perf_amr = 1
        csv = 0
        threshold = 0

        resultFiles = request.files.getlist('files[]')
        resultForm = ast.literal_eval(request.form['user_data'])

        if 'serotype' in resultForm:
            s = resultForm['serotype']
        if 'virulence_factors' in resultForm:
            vf = resultForm['virulence_factors']
        if 'amr' in resultForm:
            amr = resultForm['amr']
        if 'percent_identity' in resultForm:
            PERC_ID = resultForm['percent_identity']
        if 'percent_length' in resultForm:
            PERC_LEN = resultForm['percent_length']
        if 'verbose' in resultForm:
            VERBOSITY = resultForm['verbose']
        if 'all_vfs' in resultForm:
            avf = resultForm['all_vfs']
        if 'minimum' in resultForm:
            threshold = resultForm['minimum']
        if 'perfect_amr' in resultForm:
            perf_amr = resultForm['perfect_amr']
        if 'csv' in resultForm:
            csv = resultForm['csv']

        if len(resultFiles) == 1:
            filename = resultFiles[0].filename
            if not filename:
                return 'No files were uploaded.'
            else:
                files.save(resultFiles[0])
                OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "tools_controller.py", "--input", files.path(filename),
                                                  "-s", str(s), "-vf", str(vf), "-amr", str(amr), "-min", str(threshold), "-p", str(perf_amr), "-avf", str(avf),
                                                  "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-sv', str(VERBOSITY), '-csv', str(csv)])
        else:
            if os.path.isdir(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I)):
                shutil.rmtree(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))

            os.makedirs(SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I))

            for file in resultFiles:
                files.save(file,'temp_dir'+ str(I))
            OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "tools_controller.py", "--input", SCRIPT_DIRECTORY + '../../temp/Uploads/temp_dir' + str(I),
                                              "-s", str(s), "-vf", str(vf), "-amr", str(amr), "-min", str(threshold), "-p", str(perf_amr), "-avf", str(avf),
                                              "-pl", str(PERC_LEN), "-pi", str(PERC_ID), '-sv', str(VERBOSITY), '-csv', str(csv)])
        I +=1
        return jsonify(ast.literal_eval(OUTPUT))
    return 'No HTTP requests were made.'



@app.route('/ectyper/results', methods=['GET'])
def getResults():
    global FORMAT
    logging.info('In getResults')

    if 'Error' in OUTPUT:
        logging.info('No valid files uploaded')
        root = Tkinter.Tk()
        root.withdraw()
        tkMessageBox.showwarning('Oops!','No valid files were uploaded. Valid files are: .fasta, .fsa_nt.')
        return render_template('controller.html')

    elif FORMAT == '1':
        logging.info('To HTML table')
        totable_list = toTableList(ast.literal_eval(OUTPUT), VERBOSITY)
        return render_template('table_results.html', result=totable_list)
    else:
        logging.info('To JSON format')
        return render_template('JSON_results.html', result=json.dumps(ast.literal_eval(OUTPUT)))


@app.route('/ectyper/download', methods=['GET'])
def straightDownload():
    global OUTPUT
    if FORMAT == '0':
        json_result = json.dumps(ast.literal_eval(OUTPUT))
        return Response(
            json_result,
            mimetype='application/json',
            headers={"Content-disposition": "attachment; filename=Results_Summary.json"}
        )
    else:
        csv_result = getCSV(ast.literal_eval(OUTPUT), VERBOSITY)
        return Response(
            csv_result,
            mimetype='text/csv',
            headers={"Content-disposition": "attachment; filename=Results_Summary.csv"}
        )




@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
    app.run(debug=True)