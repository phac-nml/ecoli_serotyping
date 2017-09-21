from subprocess import Popen, PIPE
from Bio.SeqIO import SeqRecord, write
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import definitions
import os

def samToFasta(sam_file):
    '''
    Return file path of fasta output
    Convert a .SAM file to .fasta file
    '''
    fasta_file = os.path.splitext(sam_file)[0]+'.fasta'
    file_name_no_ext = os.path.splitext(os.path.basename(sam_file))[0]
    cmd = ['samtools', 'view', '-S', sam_file]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(
            "%r failed\n\tstatus code %s\n\tstderr %r" % (
            cmd, process.returncode, process.stderr))
    fasta_data = []
    decoded_stdout = stdout.decode('ascii')
    splited_stdout = decoded_stdout.split('\n')
    for line in splited_stdout:
        if line == '':
            continue
        columns = line.split("\t")
        if len(columns) < 10:
            return
        new_record = SeqRecord(
            Seq(columns[9],IUPAC.IUPACUnambiguousDNA),
            id=file_name_no_ext,
            description="")
        fasta_data.append(new_record)

    write(fasta_data, fasta_file, 'fasta')
    return fasta_file

def build_bowtie_index(reference_in):
    '''
    Return the path to the director of bowtie index
    '''
    index_name = 'bowtie_index'
    if os.path.isdir(reference_in):
        # concatenate all files if input is a directory
        referece_files = os.listdir(reference_in)
        reference_in = referece_files.join(',')
    else:
        # set index_name as input name if input is a file
        index_name = os.path.splitext(os.path.split(reference_in)[1])[0]
    reference_in_dir = os.path.split(reference_in)[0]+'/'
    index_out_dir = reference_in_dir+ 'bowtie_index/'
    index_out = index_out_dir + index_name

    # make directory if it does not already exist
    if not os.path.exists(index_out_dir):
        os.makedirs(index_out_dir)
    cmd = ['bowtie2-build', reference_in, index_out]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout_data, stderr_data = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(
            "%r failed\n\tstatus code %s\n\tstdout %r\n\tstderr %r" % (
            cmd, process.returncode, stdout_data.decode('ascii'), stderr_data.decode('ascii')))
    return index_out_dir

def align_reads(reads_in):
    '''
    Return the path of the .sam file generated
    '''
    reference_index = definitions.REFERENCE_INDEX
    sam_out = os.path.split(reads_in)[0] + '/' + os.path.splitext(os.path.basename(reads_in))[0] + '.sam'
    
    cmd = ['bowtie2', '-x', reference_index, '-U', reads_in, '-S', sam_out]
    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout_data, stderr_data = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(
            "%r failed\n\tstatus code %s\n\tstdout %r\n\tstderr %r" % (
            cmd, process.returncode, stdout_data.decode('ascii'), stderr_data.decode('ascii')))

    return sam_out
    # bowtie2 -x Data/test_reads/bowtie_index/EcOH -U Data/test_reads/reads_1.fq -S Data/test_reads/reads_1.sam