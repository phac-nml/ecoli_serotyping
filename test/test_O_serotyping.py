import sys
import os
import logging
from ectyper import ectyper
from ectyper.predictionFunctions import quality_control_results
LOG=logging.getLogger("TEST")
TEST_ROOT = os.path.dirname(__file__)

def set_input(input,
              percent_iden=None,
              #output=tempfile.mkdtemp(),
              output="tmp",
              cores=1,
              print_sequence=False,
              sketchpath=False,
              debug=False):
    """
    Create the sys.argv[] without need for commandline input.
    :param input: Input file given by testing function
    :param percent_iden: Percent identity for comparison
    :param output: Location of output
    :return: None
    """
    args = ['-i', input,
            '--verify',
            '-c', str(cores)
            ]
    if sketchpath:
        args+=['-r', os.path.join(TEST_ROOT, 'Data/refseqsketch/refseq.genomes.k21s1000.msh')]
    if percent_iden:
        args += ['-d', str(percent_iden)]
    if output:
        args += ['-o', output]
    if print_sequence:
        args += ['--sequence']
    if debug:
        args+=['--debug']
    sys.argv[1:] = args

def test_Otyping(caplog):
    """
    Giving E.coli fasta genomes with truncated wzx and wzy genes with reference coverage <50 predict O and H antigens
    :return: None
    """
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT, 'Data/Escherichia_O26H11.fasta')#+","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(input=file,cores=4,print_sequence=True)
    ectyper.run_program()


    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         secondrow = outfp.readlines()[1].split("\t")
         Otype = secondrow[2]
         Htype = secondrow[3]

    assert Otype == "O26", "Expected O26 but reported O-type:" + Otype
    assert Htype == "H11", "Expected H11 but reported H-type:" + Htype


def test_closeOalles_O42_O28(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EscherichiaO28H5.fasta')  # +","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(input=file, cores=4, print_sequence=True, debug=False)
    ectyper.run_program()
    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         secondrow = outfp.readlines()[1]
    print(secondrow)
    assert "EscherichiaO28H5\tEscherichia coli\tO28\tH25\tO28:H25\tPASS\tHIGH" in secondrow

def test_QCmodule(): #QCmodule
    final_dict={'Sample1': {'O': 'O26',
                                       'H': 'H11',
                                       'H11': {'fliC': 1.0,
                                               '≈fliC': 'AAC'},
                                       'O26': {'≈wzx' : 'ATG',
                                               'wzx'  : 0.46,
                                               '≈wzy' : 'CTT',
                                               'wzy'  : 0.08
                                               },
                                        'species': 'Escherichia coli'
                                       }
                }

    assert  quality_control_results("Sample1",final_dict) == {'AlleleNames': ['wzx', 'wzy', 'fliC'], 'NumberOfAlleles': 3, 'QCflag': 'PASS', 'ConfidenceLevel': 'LOW'}


    final_dict = {'Sample2': {'O': '-',
                                         'H': 'H11',
                                         'H11': {'fliC': 1.0,
                                                 '≈fliC': 'AAC'
                                                },
                                         'species': 'Escherichia coli'
                                         }
                  }


    assert quality_control_results("Sample2", final_dict) == {'AlleleNames': ['fliC'], 'NumberOfAlleles': 1, 'QCflag': 'FAIL','ConfidenceLevel': '-'}

    final_dict = {'Sample3': {'O': 'O26',
                                         'O26': {'≈wzx' : 'ATG','wzx'  : 0.56,
                                                },
                                         'H': '-',
                                          'species': 'Escherichia coli'
                                         }
                  }

    assert quality_control_results("Sample3", final_dict) == {'AlleleNames': ['wzx'], 'NumberOfAlleles': 1, 'QCflag': 'PASS','ConfidenceLevel': 'LOW'}

    final_dict = {'Sample4': {'O': 'O26',
                              'O26': {'≈wzx': 'ATG','wzx': 0.95,
                                      '≈wzy': 'ATG','wzy': 0.87
                                      },
                              'H': '-',
                              'species': 'Escherichia coli'
                              }
                  }

    assert quality_control_results("Sample4", final_dict) == {'AlleleNames': ['wzx','wzy'], 'NumberOfAlleles': 2,
                                                              'QCflag': 'PASS','ConfidenceLevel': 'MEDIUM'}

    final_dict = {'Sample5': {'O': 'O26',
                              'O26': {'≈wzx': 'ATG', 'wzx': 0.99
                                      },
                              'H': '-',
                              'species': 'Escherichia coli'
                              }
                  }

    assert quality_control_results("Sample5", final_dict) == {'AlleleNames': ['wzx'], 'NumberOfAlleles': 1,
                                                              'QCflag': 'PASS', 'ConfidenceLevel': 'HIGH'}


def test_Shigella_typing(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/DRR015915_Shigella_boydii.fasta')  # +","+os.path.join(TEST_ROOT, 'Data/Escherichia.fna')
    set_input(input=file, cores=4, print_sequence=True)
    ectyper.run_program()
    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         secondrow = outfp.readlines()[1].split("\t")
         species = secondrow[1]
    assert species == "Shigella boydii"

def test_download_refseq_mash(caplog):
    caplog.set_level(logging.DEBUG)
    response = ectyper.speciesIdentification.get_refseq_mash()
    assert response == True,"Something went wrong with the RefSeq sketch download. Check internet connection ..."


def test_mixofspecies(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/Campylobacter.fasta') +","+os.path.join(TEST_ROOT, 'Data/Salmonella.fasta')+","\
                         + os.path.join(TEST_ROOT, 'Data/Escherichia.fastq')
    set_input(input=file, cores=4, print_sequence=True, sketchpath=False)
    ectyper.run_program()

    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         rows = outfp.readlines()
    rows=rows[1:] #remove header line

    serovars=[]; genomenames=[]; QCflag=[]; confidence=[]
    for row in rows:
        rowlist = row.split("\t")
        print(rowlist)
        serovars.append(rowlist[4])
        genomenames.append(rowlist[0])
        QCflag.append(rowlist[5])
        confidence.append(rowlist[6])

    assert serovars == ['-:-', 'O22:H8', '-:-']
    assert genomenames == ["Campylobacter","Escherichia","Salmonella"]
    assert QCflag == ["-","PASS","-"]
    assert confidence == ["-","HIGH","-"]

def test_Ealbertii_1(caplog): #error
    LOG.info("Starting 1 of 3 test on EnteroBase on sample ESC_HA8355AA_AS: Escherichia albertii O65:H5")
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/ESC_HA8355AA_AS_Ealberii_O65H5.fasta')
    set_input(input=file, cores=4, print_sequence=True, sketchpath=False)
    ectyper.run_program()
    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #remove header line
    print(secondrow)
    assert "ESC_HA8355AA_AS_Ealberii_O65H5\tEscherichia albertii\tO65\tH5\tO65:H5\tNA" in secondrow  #suggest to be H52 and not H5

def test_Ealbertii_2(): #error
    LOG.info("Starting 2 of 3 test on EnteroBase on sample on ESC_HA8509AA_AS: Escherichia albertii O5:H5")

    file = os.path.join(TEST_ROOT,
                        'Data/ESC_HA8509AA_AS_EalbertiiO5H5.fasta')
    set_input(input=file, cores=4, print_sequence=True, sketchpath=False)
    ectyper.run_program()

    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #remove header line
    print(secondrow)
    assert "ESC_HA8509AA_AS_EalbertiiO5H5\tEscherichia albertii\t-\tH5\t-:H5\tNA" in secondrow

def test_Ealbertii_3(caplog):
    LOG.info("Starting 3 of 3 test Escherichia albertii O49:NM") #can not type O49 due to poor sequence quality of uncertainty of wet-lab O49 typing
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/Ealbertii_O49NM.fasta')
    set_input(input=file, cores=4, print_sequence=True, sketchpath=False)
    ectyper.run_program()

    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #remove header line
    print(secondrow)
    assert "Ealbertii_O49NM\tEscherichia albertii\t-\tH5\t-:H5\tNA" in secondrow


def test_Ecoli_O17H18(caplog):
    caplog.set_level(logging.DEBUG)
    file = os.path.join(TEST_ROOT,
                        'Data/EscherichiaO17H18.fasta')
    set_input(input=file, cores=4, print_sequence=False, sketchpath=False,debug=True)
    ectyper.run_program()

    with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
         rows = outfp.readlines()
    secondrow=rows[1:][0] #remove header line
    print(secondrow)
    assert "Escherichia coli\tO17\tH18\tO17:H18\tPASS\tHIGH" in secondrow


# def test_EC20151996():
#     file = os.path.join(TEST_ROOT,
#                         'Data/EC20151996rawreads.fastq')
#     set_input(input=file, cores=4, print_sequence=False, sketchpath=False, debug=True)
#     ectyper.run_program()
#     with open(os.path.join(TEST_ROOT,"tmp/output.tsv")) as outfp:
#          rows = outfp.readlines()
#     secondrow=rows[1:][0] #remove header line
#     print(secondrow)
#     assert "EC20151996rawreads\tEscherichia albertii\tO3\tH17\tO3:H17\tNA\t-\tBased on 3 allele(s)" in secondrow