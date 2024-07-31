import sys
import pytest
import tempfile
import os, json
from ectyper import ectyper



TEST_ROOT = os.path.dirname(__file__)
ROOT_DIR = os.path.abspath(os.path.join(TEST_ROOT, '..'))
MASHINFILE = os.path.join(ROOT_DIR, 'ectyper/Data/EnteroRef_GTDBSketch_20231003_V2.msh')

def set_input(input,
              percent_iden=None,
              verify=True,
              output=tempfile.mkdtemp(),
              cores=1,
              print_sequence=False,
              refseqmash = False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    print(input)
    args = ['-i', input,
            '-c', str(cores),
            ]

    if percent_iden:
        args += ['-d', str(percent_iden)]
    if verify:
        args += ['--verify']
    if output:
        args += ['-o', output]
    if print_sequence:
        args += ['--sequence']
    if refseqmash:
        args += ['-r', MASHINFILE]

    sys.argv[1:] = args

def test_emtpy_BLAST_antigen_hits(tmpdir):
    """
    fastq2fasta_missing_OH_BLAST_results.fasta is a minimal input file generated from FASTQ file that has almost no O- and H- antigens resulting in empty BLAST results
    """
    fastafile=os.path.join(TEST_ROOT, 'Data/fastq2fasta_missing_OH_BLAST_results.fasta')
    set_input(input=fastafile, verify=False, refseqmash=True)
    args = ectyper.commandLineOptions.parse_command_line()
    with open(file=ectyper.definitions.SEROTYPE_ALLELE_JSON) as fp:
            ectyperdb_dict = json.load(fp)
    alleles_fasta_file = ectyper.create_alleles_fasta_file(tmpdir, ectyperdb_dict)
    db_prediction_dict = ectyper.genome_group_prediction(args.input, alleles_fasta_file, args, tmpdir, ectyperdb_dict)
    assert db_prediction_dict == {}

def test_invalid_fasta():
     fastafile=os.path.join(TEST_ROOT, 'Data/invalid_fasta.fasta')
     set_input(input=fastafile, verify=True, refseqmash=True)
     args = ectyper.commandLineOptions.parse_command_line()
     ecoli_files_dict, other_files_dict,filesnotfound_dict = ectyper.speciesIdentification.verify_ecoli_and_inputs(fasta_fastq_files_dict={fastafile:''}, ofiles=[], filesnotfound=[], args=args)
     assert ecoli_files_dict == {}
     assert other_files_dict['invalid_fasta']['error'] != ''
     assert 'is invalid/empty' in other_files_dict['invalid_fasta']['error'], f'"is invalid/empty" not found in "{other_files_dict['invalid_fasta']['error']}"' 
     assert filesnotfound_dict == {}

     
    
