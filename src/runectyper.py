#!/usr/bin/env python

from ectyper import *
from flask_uploads import *
from flask import *
import subprocess
import ast
import formatresults

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
OUTPUT = {}
RESULTS = False
I = 0


app = Flask(__name__)
app.config['UPLOADED_FASTAFILES_DEST'] = SCRIPT_DIRECTORY + '../temp/Uploads'
app.config['UPLOADED_FASTAFILES_ALLOW'] = set(['fasta', 'fsa_nt'])
files = UploadSet('fastafiles')
configure_uploads(app, (files,))

@app.route('/ectyper/upload', methods =['POST', 'GET'])
def uploadFiles():
    if request.method == 'POST':
       global OUTPUT
       global I
       global RESULTS

       perc_id = request.form['perc-id']
       perc_len = request.form['perc-len']
       resultFiles = request.files.getlist('file')
       verbosity = request.form['verbosity']
       if 'table-checkbox' in request.form:
            RESULTS = True
       else:
            RESULTS = False

       if len(resultFiles) == 1:
           filename = resultFiles[0].filename
           if not filename:
            return 'No files were uploaded.'
           else:
            files.save(resultFiles[0])
            OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "-input", files.path(filename),
                                                "-pl", perc_len, "-pi", perc_id, '-v', verbosity])
       else:
           os.makedirs(SCRIPT_DIRECTORY + '../temp/Uploads/temp_dir' + str(I))
           for file in resultFiles:
               files.save(file,'temp_dir'+ str(I))
           OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "-input", SCRIPT_DIRECTORY + '../temp/Uploads/temp_dir' + str(I),
                                                "-pl", perc_len, "-pi", perc_id, '-v', verbosity])
       I +=1
       return redirect(url_for('getResults'))
    return render_template('uploadfile.html')

@app.route('/ectyper/results', methods=['GET'])
def getResults():
    if RESULTS:
        return formatresults.toHTML()
    else:
       return jsonify(ast.literal_eval(OUTPUT))


@app.errorhandler(404)
def not_found(error):
  return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
    app.run(debug=True)