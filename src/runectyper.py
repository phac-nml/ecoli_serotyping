#!/usr/bin/env python

from ectyper import *
from flask_uploads import *
from flask import *
import subprocess

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
OUTPUT = {}


app = Flask(__name__)
app.config['UPLOADED_FASTAFILES_DEST'] = SCRIPT_DIRECTORY + '../temp/Uploads'
app.config['UPLOADED_FASTAFILES_ALLOW'] = set(['fasta', 'fsa_nt'])
files = UploadSet('fastafiles')
configure_uploads(app, (files,))

@app.route('/ectyper/upload', methods =['POST', 'GET'])
def uploadFiles():
    if request.method == 'POST':
       file = request.files['file']
       if not file:
            return 'No files were uploaded.'
       global OUTPUT
       filename = file.filename
       files.save(request.files['file'])
       OUTPUT = subprocess.check_output([SCRIPT_DIRECTORY + "ectyper.py", "-input", files.path(filename)])
       return redirect(url_for('getResults'))
    return render_template('uploadfile.html')

@app.route('/ectyper/results', methods=['GET'])
def getResults():
    return jsonify({'GENOMES' :OUTPUT})

@app.errorhandler(404)
def not_found(error):
  return make_response(jsonify({'Error': 'Not found.'}), 404)

if __name__=='__main__':
    app.run(debug=True)