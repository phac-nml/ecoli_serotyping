#!usr/bin/env python

import src.blastFunctions
import src.genomeFunctions
import src.serotypePrediction

def test_parse_serotype():
    """

    """

    blast_record_type1 = {
        'qseqid': '1__fliC__fliC-H1__1 AB028471.1;flagellin;H1',
        'qlen': 900,
        'sseqid': 'genomeA',
        'length': 900,
        'pident': 100,
        'sstart':2,
        'send':901,
        'sframe':1
    }

    assert src.genomeFunctions.get_genome_name(
        blast_record_type1['sseqid']) == 'genomeA'

    serotype_parser = src.genomeFunctions.get_parsing_dict('serotype')
    type1Result = src.serotypePrediction.parse_serotype(blast_record_type1,
                                                            serotype_parser[
                                                             'data'])

    assert type1Result == {'fliC':{'antigen':'H1',
                                   'blast_record':blast_record_type1}}