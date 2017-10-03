#!usr/bin/env python

import src.blastFunctions
import src.genomeFunctions
import src.serotypePrediction


def test_parse_serotype():
    """

    """
    dict_of_qseqid = {'8__wzx__wzx-O1__112 GU299791.1;O antigen '
                      'flippase;O1': 'O1',
                      '1__fliC__fliC-H1__1 AB028471.1;flagellin;H1': 'H1',
                      '1__fliC__fliC-H19__19 AY250002.1;flagellin;H19': 'H19',
                      '2__flkA__flkA-H3__77 AB128916;flagellin;H3': 'H3',
                      '3__fllA__fllA-H55__83 AB269771;flagellin;H55': 'H55',
                      '4__flmA__flmA-H54__84 AB128918;flagellin;H54': 'H54',
                      '5__flnA__flnA-H17__85 CP002291;flagellin;H17': 'H17',
                      '6__wzm__wzm-O101-Gp15__86 AB812046.1;O-antigen ABC '
                      'transporter permease;O101': 'O101',
                      '7__wzt__wzt-O60__102 AB812022.1;ATP-binding '
                      'protein;O60': 'O60',
                      '8__wzx__wzx-O10__113 AB811599.1;O antigen '
                      'flippase;O10': 'O10',
                      '9__wzy__wzy-O173__409 GU068046.1;O antigen '
                      'polyermase;O173': 'O173',
                      '9__wzy__wzy-O116var1__512 ERS085345;O antigen '
                      'polyermase;O116-like': 'O116',
                      '9__wzy__wzy-O108var1__522 ERS150902;O antigen '
                      'polyermase;O108-like': 'O108',
                      '8__wzx__wzx-Onovel5__548 SRR2970275;O antigen '
                      'flippase;Onovel5':None,
                      'fliC_54_AJ605766_H17':'H17',
                      'fliC_111197_JH965342_H29':'H29',
                      'fllA_1_AB269771_H55':'H55',
                      'wzx_36_DQ462205-FJ539194_O28ac':'O28',
                      'wzx_103_KB021482_O104':'O104',
                      'wzy_152_AIGX01000028_O45':'O45',
                      'wzx_406_56-54Cigleris_O128ab':'O128',
                      'wzt_3_AB010293_O9':'O9',
                      'wzm_3_AB010293_O9':'O9',
                      'gnd|14_O103':'O103'
                      }

    blast_record_type1 = {
        'qseqid': '1__fliC__fliC-H1__1 AB028471.1;flagellin;H1',
        'qlen': 900,
        'sseqid': 'genomeA',
        'length': 900,
        'pident': 100,
        'sstart': 2,
        'send': 901,
        'sframe': 1
    }

    assert src.genomeFunctions.get_genome_name(
        blast_record_type1['sseqid']) == 'genomeA'

    serotype_parser = src.genomeFunctions.get_parsing_dict('serotype')
    # type1Result = src.serotypePrediction.parse_serotype(blast_record_type1,
    #                                                     serotype_parser[
    #                                                         'data'])

    # assert type1Result == {'fliC': {'antigen': 'H1',
    #                                 'blast_record': blast_record_type1}}

    # for qseqid in dict_of_qseqid.keys():
    #     ps = src.serotypePrediction.parse_serotype(
    #         {'qseqid': qseqid}, serotype_parser['data']
    #     )

        # for name in ps.keys():
        #     assert ps[name]['antigen'] == dict_of_qseqid[qseqid]