
import src.virulencePrediction
import src.serotypePrediction
import src.genomeFunctions
import src.blastFunctions
from argparse import Namespace
import argparse
import pandas as pd
import re
import shutil
import os.path
from pathlib import Path
import definitions
import pickle

CONST_CSV = False
CONST_INPUT = 'filename'
CONST_MINGENOMES = 1
CONST_PERIDENT = 90
CONST_PERLEN = 90
CONST_SER = True
CONST_VIR = True





def test_sero_pred() :
    testing_dict = {"JHJE01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHJE01000023.1", "qlen": "643", "send": "23634", "pident": "95.17", "sstart": "24275", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JASU01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JASU01000004.1", "qlen": "1497", "send": "106075", "pident": "99.53", "sstart": "104587", "qseqid": "1__fliC__fliC-H2__20", "length": "1489", "sframe": "1"}, "antigen": "H2"}, "gnd": {"blast_record": {"sseqid": "JASU01000004.1", "qlen": "643", "send": "4950", "pident": "95.61", "sstart": "4313", "qseqid": "gnd|3_O26", "length": "638", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AUQB01": {"otype": {"ant_number": "O104", "strength": "99.88"}, "serotype": {"wzx": {"blast_record": {"sseqid": "AUQB01000007.1", "qlen": "867", "send": "142231", "pident": "99.88", "sstart": "143098", "qseqid": "wzx_106_AFOB02000091_O104", "length": "868", "sframe": "1"}, "antigen": "O104"}, "gnd": {"blast_record": {"sseqid": "AUQB01000007.1", "qlen": "643", "send": "140285", "pident": "93.76", "sstart": "140924", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHJY01": {"otype": {"ant_number": "O157", "strength": "99.41"}, "serotype": {"wzx": {"blast_record": {"sseqid": "JHJY01000064.1", "qlen": "1356", "send": "13868", "pident": "99.41", "sstart": "12512", "qseqid": "wzx_200_AKLI01000048_O157", "length": "1357", "sframe": "1"}, "antigen": "O157"}, "fliC": {"blast_record": {"sseqid": "JHJY01000002.1", "qlen": "1578", "send": "31633", "pident": "99.81", "sstart": "30064", "qseqid": "1__fliC__fliC-H16__13", "length": "1570", "sframe": "1"}, "antigen": "H16"}, "gnd": {"blast_record": {"sseqid": "JHJY01000015.1", "qlen": "643", "send": "21997", "pident": "95.32", "sstart": "22637", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "JHJY01000064.1", "qlen": "1185", "send": "11775", "pident": "99.41", "sstart": "10590", "qseqid": "9__wzy__wzy-O157__388", "length": "1186", "sframe": "1"}, "antigen": "O157"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHLA01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHLA01000026.1", "qlen": "1467", "send": "29170", "pident": "99.93", "sstart": "27712", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "JHLA01000102.1", "qlen": "643", "send": "21307", "pident": "94.54", "sstart": "21947", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AHZD02": {"otype": {"ant_number": "O157", "strength": "99.71"}, "serotype": {"wzx": {"blast_record": {"sseqid": "AHZD02000001.1", "qlen": "1356", "send": "597864", "pident": "99.71", "sstart": "596508", "qseqid": "wzx_200_AKLI01000048_O157", "length": "1357", "sframe": "1"}, "antigen": "O157"}, "gnd": {"blast_record": {"sseqid": "AHZD02000001.1", "qlen": "643", "send": "608771", "pident": "95.48", "sstart": "608131", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "AHZD02000001.1", "qlen": "741", "send": "595529", "pident": "99.60", "sstart": "594788", "qseqid": "wzy_199_AKLQ01000042_O157", "length": "742", "sframe": "1"}, "antigen": "O157"}}, "htype": {"ant_number": "-", "strength": 0}}, "BBUR01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "BBUR01000041.1", "qlen": "643", "send": "23954", "pident": "93.30", "sstart": "23313", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGSG01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "AGSG01000154.1", "qlen": "643", "send": "34526", "pident": "96.26", "sstart": "35166", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHJX01": {"otype": {"ant_number": "O86", "strength": "99.33"}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHJX01000037.1", "qlen": "643", "send": "9214", "pident": "95.79", "sstart": "8574", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "JHJX01000037.1", "qlen": "1341", "send": "6396", "pident": "99.33", "sstart": "5053", "qseqid": "9__wzy__wzy-O86__504", "length": "1344", "sframe": "1"}, "antigen": "O86"}}, "htype": {"ant_number": "-", "strength": 0}}, "AJQW01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "AJQW01000001.1", "qlen": "643", "send": "259080", "pident": "93.60", "sstart": "259720", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHMM01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHMM01000032.1", "qlen": "1332", "send": "30790", "pident": "96.89", "sstart": "29471", "qseqid": "fliC_311_AGSG01000116_H25", "length": "1320", "sframe": "1"}, "antigen": "H25"}, "gnd": {"blast_record": {"sseqid": "JHMM01000025.1", "qlen": "643", "send": "9695", "pident": "94.23", "sstart": "9055", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHOC01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHOC01000031.1", "qlen": "1479", "send": "8751", "pident": "94.19", "sstart": "7273", "qseqid": "1__fliC__fliC-H40__46", "length": "1480", "sframe": "1"}, "antigen": "H40"}, "gnd": {"blast_record": {"sseqid": "JHOC01000023.1", "qlen": "643", "send": "20091", "pident": "95.76", "sstart": "20727", "qseqid": "gnd|1_O26", "length": "637", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JFBE01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JFBE01000053.1", "qlen": "1788", "send": "9402", "pident": "99.33", "sstart": "7616", "qseqid": "1__fliC__fliC-H12__9", "length": "1787", "sframe": "1"}, "antigen": "H12"}, "gnd": {"blast_record": {"sseqid": "JFBE01000001.1", "qlen": "643", "send": "357789", "pident": "94.39", "sstart": "357148", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHHC01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHHC01000004.1", "qlen": "1758", "send": "151501", "pident": "99.72", "sstart": "153257", "qseqid": "1__fliC__fliC-H7__70", "length": "1757", "sframe": "1"}, "antigen": "H7"}, "gnd": {"blast_record": {"sseqid": "JHHC01000011.1", "qlen": "643", "send": "92154", "pident": "93.93", "sstart": "91513", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JASR01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JASR01000005.1", "qlen": "643", "send": "126708", "pident": "93.46", "sstart": "127349", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "FM180568.1": {"otype": {"ant_number": "O90", "strength": "99.23"}, "serotype": {"gnd": {"blast_record": {"sseqid": "FM180568.1", "qlen": "643", "send": "2218784", "pident": "96.10", "sstart": "2219424", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "FM180568.1", "qlen": "1158", "send": "2221571", "pident": "99.23", "sstart": "2222737", "qseqid": "9__wzy__wzy-O90-Gp4__507", "length": "1167", "sframe": "1"}, "antigen": "O90"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHLP01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHLP01000018.1", "qlen": "643", "send": "44845", "pident": "96.10", "sstart": "44205", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGTH01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "AGTH01000188.1", "qlen": "643", "send": "64859", "pident": "94.23", "sstart": "64219", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGTJ01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AGTJ01000023.1", "qlen": "1833", "send": "74760", "pident": "98.36", "sstart": "76591", "qseqid": "1__fliC__fliC-H19__18", "length": "1832", "sframe": "1"}, "antigen": "H19"}, "gnd": {"blast_record": {"sseqid": "AGTJ01000128.1", "qlen": "643", "send": "61427", "pident": "93.93", "sstart": "62068", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHJJ01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHJJ01000027.1", "qlen": "643", "send": "4317", "pident": "94.24", "sstart": "3676", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AUZS01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AUZS01000098.1", "qlen": "1578", "send": "66600", "pident": "99.94", "sstart": "68169", "qseqid": "1__fliC__fliC-H16__13", "length": "1570", "sframe": "1"}, "antigen": "H16"}, "gnd": {"blast_record": {"sseqid": "AUZS01000181.1", "qlen": "643", "send": "39416", "pident": "95.16", "sstart": "40056", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AE005174.2": {"otype": {"ant_number": "O157", "strength": "99.93"}, "serotype": {"wzx": {"blast_record": {"sseqid": "AE005174.2", "qlen": "1356", "send": "2853432", "pident": "99.93", "sstart": "2854788", "qseqid": "wzx_200_AKLI01000048_O157", "length": "1357", "sframe": "1"}, "antigen": "O157"}, "fliC": {"blast_record": {"sseqid": "AE005174.2", "qlen": "1758", "send": "2699592", "pident": "99.94", "sstart": "2701348", "qseqid": "1__fliC__fliC-H7__70", "length": "1757", "sframe": "1"}, "antigen": "H7"}, "gnd": {"blast_record": {"sseqid": "AE005174.2", "qlen": "643", "send": "2842523", "pident": "93.45", "sstart": "2843163", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "AE005174.2", "qlen": "741", "send": "2855767", "pident": "99.87", "sstart": "2856508", "qseqid": "wzy_199_AKLQ01000042_O157", "length": "742", "sframe": "1"}, "antigen": "O157"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHKX01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHKX01000022.1", "qlen": "1578", "send": "142546", "pident": "99.94", "sstart": "144115", "qseqid": "1__fliC__fliC-H16__13", "length": "1570", "sframe": "1"}, "antigen": "H16"}, "gnd": {"blast_record": {"sseqid": "JHKX01000028.1", "qlen": "643", "send": "22791", "pident": "93.46", "sstart": "23432", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AP010953.1": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AP010953.1", "qlen": "1467", "send": "2715723", "pident": "99.93", "sstart": "2717181", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "AP010953.1", "qlen": "643", "send": "2854070", "pident": "95.48", "sstart": "2854710", "qseqid": "gnd|3_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "CP000800.1": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "CP000800.1", "qlen": "643", "send": "2282905", "pident": "95.17", "sstart": "2283546", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHFY01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHFY01000043.1", "qlen": "1467", "send": "29069", "pident": "99.93", "sstart": "27611", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "JHFY01000107.1", "qlen": "643", "send": "24140", "pident": "96.26", "sstart": "24780", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGTI01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "AGTI01000214.1", "qlen": "643", "send": "19836", "pident": "95.48", "sstart": "20477", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHNM01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHNM01000013.1", "qlen": "1758", "send": "64322", "pident": "99.72", "sstart": "66078", "qseqid": "1__fliC__fliC-H7__70", "length": "1757", "sframe": "1"}, "antigen": "H7"}, "gnd": {"blast_record": {"sseqid": "JHNM01000013.1", "qlen": "643", "send": "149493", "pident": "93.76", "sstart": "150133", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHNV01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHNV01000079.1", "qlen": "1050", "send": "80528", "pident": "99.13", "sstart": "81565", "qseqid": "1__fliC__fliC-H4__42", "length": "1038", "sframe": "1"}, "antigen": "H4"}, "gnd": {"blast_record": {"sseqid": "JHNV01000080.1", "qlen": "643", "send": "50713", "pident": "93.30", "sstart": "50072", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "CP001855.1": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "CP001855.1", "qlen": "1788", "send": "1970831", "pident": "97.59", "sstart": "1972617", "qseqid": "1__fliC__fliC-H12__9", "length": "1787", "sframe": "1"}, "antigen": "H12"}, "gnd": {"blast_record": {"sseqid": "CP001855.1", "qlen": "643", "send": "2128359", "pident": "94.85", "sstart": "2128999", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AJVS01": {"otype": {"ant_number": "O111", "strength": "99.81"}, "serotype": {"fliC": {"blast_record": {"sseqid": "AJVS01000516.1", "qlen": "1467", "send": "100018", "pident": "99.93", "sstart": "101476", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "AJVS01000357.1", "qlen": "643", "send": "24352", "pident": "100.00", "sstart": "24993", "qseqid": "gnd|62_O111", "length": "642", "sframe": "1"}, "antigen": "O111"}, "wzy": {"blast_record": {"sseqid": "AJVS01000357.1", "qlen": "1053", "send": "27545", "pident": "99.81", "sstart": "28598", "qseqid": "9__wzy__wzy-O111__331", "length": "1054", "sframe": "1"}, "antigen": "O111"}}, "htype": {"ant_number": "-", "strength": 0}}, "AFOB02": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AFOB02000091.1", "qlen": "1050", "send": "10204", "pident": "99.23", "sstart": "11241", "qseqid": "1__fliC__fliC-H4__42", "length": "1038", "sframe": "1"}, "antigen": "H4"}, "gnd": {"blast_record": {"sseqid": "AFOB02000091.1", "qlen": "643", "send": "153567", "pident": "93.92", "sstart": "154207", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AQGP01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AQGP01000031.1", "qlen": "1497", "send": "96129", "pident": "99.60", "sstart": "94641", "qseqid": "1__fliC__fliC-H2__20", "length": "1489", "sframe": "1"}, "antigen": "H2"}, "gnd": {"blast_record": {"sseqid": "AQGP01000004.1", "qlen": "643", "send": "586322", "pident": "95.76", "sstart": "585686", "qseqid": "gnd|1_O26", "length": "637", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHNW01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHNW01000003.1", "qlen": "643", "send": "573492", "pident": "95.76", "sstart": "572856", "qseqid": "gnd|1_O26", "length": "637", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHGM01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHGM01000068.1", "qlen": "1467", "send": "27763", "pident": "99.93", "sstart": "29221", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "JHGM01000016.1", "qlen": "643", "send": "24199", "pident": "95.48", "sstart": "24839", "qseqid": "gnd|3_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AHAW01": {"otype": {"ant_number": "O90", "strength": "99.23"}, "serotype": {"fliC": {"blast_record": {"sseqid": "AHAW01000030.1", "qlen": "1461", "send": "347977", "pident": "98.09", "sstart": "349437", "qseqid": "fliC_109_AM231154_H27", "length": "1465", "sframe": "1"}, "antigen": "H27"}, "gnd": {"blast_record": {"sseqid": "AHAW01000039.1", "qlen": "643", "send": "23047", "pident": "96.10", "sstart": "23687", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "AHAW01000039.1", "qlen": "1158", "send": "25834", "pident": "99.23", "sstart": "27000", "qseqid": "9__wzy__wzy-O90-Gp4__507", "length": "1167", "sframe": "1"}, "antigen": "O90"}}, "htype": {"ant_number": "-", "strength": 0}}, "CP001846.1": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "CP001846.1", "qlen": "1758", "send": "2394221", "pident": "100.00", "sstart": "2395977", "qseqid": "1__fliC__fliC-H7__70", "length": "1757", "sframe": "1"}, "antigen": "H7"}}, "htype": {"ant_number": "-", "strength": 0}}, "AP010958.1": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AP010958.1", "qlen": "1497", "send": "2243220", "pident": "99.60", "sstart": "2244708", "qseqid": "1__fliC__fliC-H2__20", "length": "1489", "sframe": "1"}, "antigen": "H2"}, "gnd": {"blast_record": {"sseqid": "AP010958.1", "qlen": "643", "send": "2510848", "pident": "96.26", "sstart": "2511488", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHMH01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHMH01000073.1", "qlen": "1497", "send": "11321", "pident": "99.60", "sstart": "12809", "qseqid": "1__fliC__fliC-H2__20", "length": "1489", "sframe": "1"}, "antigen": "H2"}, "gnd": {"blast_record": {"sseqid": "JHMH01000098.1", "qlen": "643", "send": "23365", "pident": "96.10", "sstart": "24005", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JASO01": {"otype": {"ant_number": "O145", "strength": "99.75"}, "serotype": {"fliC": {"blast_record": {"sseqid": "JASO01000012.1", "qlen": "1332", "send": "541217", "pident": "96.97", "sstart": "539898", "qseqid": "fliC_311_AGSG01000116_H25", "length": "1320", "sframe": "1"}, "antigen": "H25"}, "gnd": {"blast_record": {"sseqid": "JASO01000012.1", "qlen": "643", "send": "423158", "pident": "95.32", "sstart": "422518", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "JASO01000012.1", "qlen": "1188", "send": "413346", "pident": "99.75", "sstart": "412158", "qseqid": "9__wzy__wzy-O145__374", "length": "1189", "sframe": "1"}, "antigen": "O145"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGTL01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AGTL01000001.1", "qlen": "1740", "send": "42390", "pident": "92.29", "sstart": "40648", "qseqid": "1__fliC__fliC-H28__29", "length": "1752", "sframe": "1"}, "antigen": "H28"}, "gnd": {"blast_record": {"sseqid": "AGTL01000031.1", "qlen": "643", "send": "21956", "pident": "95.32", "sstart": "22596", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHNJ01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"gnd": {"blast_record": {"sseqid": "JHNJ01000059.1", "qlen": "643", "send": "1699", "pident": "93.93", "sstart": "1058", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}, "AGTK01": {"otype": {"ant_number": "O45", "strength": "99.85"}, "serotype": {"fliC": {"blast_record": {"sseqid": "AGTK01000036.1", "qlen": "1497", "send": "86453", "pident": "99.60", "sstart": "87941", "qseqid": "1__fliC__fliC-H2__20", "length": "1489", "sframe": "1"}, "antigen": "H2"}, "gnd": {"blast_record": {"sseqid": "AGTK01000379.1", "qlen": "643", "send": "55280", "pident": "94.70", "sstart": "54639", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "AGTK01000379.1", "qlen": "687", "send": "50328", "pident": "99.85", "sstart": "49641", "qseqid": "wzy_152_AIGX01000028_O45", "length": "688", "sframe": "1"}, "antigen": "O45"}}, "htype": {"ant_number": "-", "strength": 0}}, "AJVU01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AJVU01000606.1", "qlen": "1479", "send": "16465", "pident": "94.19", "sstart": "17943", "qseqid": "1__fliC__fliC-H40__46", "length": "1480", "sframe": "1"}, "antigen": "H40"}, "gnd": {"blast_record": {"sseqid": "AJVU01000301.1", "qlen": "643", "send": "23327", "pident": "100.00", "sstart": "22686", "qseqid": "gnd|62_O111", "length": "642", "sframe": "1"}, "antigen": "O111"}}, "htype": {"ant_number": "-", "strength": 0}}, "JHHB01": {"otype": {"ant_number": "O123", "strength": "99.92"}, "serotype": {"fliC": {"blast_record": {"sseqid": "JHHB01000020.1", "qlen": "1467", "send": "21571", "pident": "99.93", "sstart": "23029", "qseqid": "1__fliC__fliC-H11__5", "length": "1459", "sframe": "1"}, "antigen": "H11"}, "gnd": {"blast_record": {"sseqid": "JHHB01000003.1", "qlen": "643", "send": "23142", "pident": "93.60", "sstart": "23782", "qseqid": "gnd|1_O26", "length": "641", "sframe": "1"}, "antigen": "O26"}, "wzy": {"blast_record": {"sseqid": "JHHB01000003.1", "qlen": "1251", "send": "30180", "pident": "99.92", "sstart": "31431", "qseqid": "9__wzy__wzy-O123-O186-Gp5__349", "length": "1252", "sframe": "1"}, "antigen": "O123"}}, "htype": {"ant_number": "-", "strength": 0}}, "AQGO01": {"otype": {"ant_number": "-", "strength": 0}, "serotype": {"fliC": {"blast_record": {"sseqid": "AQGO01000013.1", "qlen": "1713", "send": "42746", "pident": "97.45", "sstart": "41019", "qseqid": "1__fliC__fliC-H46__52", "length": "1728", "sframe": "1"}, "antigen": "H46"}, "gnd": {"blast_record": {"sseqid": "AQGO01000048.1", "qlen": "643", "send": "4318", "pident": "95.17", "sstart": "3677", "qseqid": "gnd|1_O26", "length": "642", "sframe": "1"}, "antigen": "O26"}}, "htype": {"ant_number": "-", "strength": 0}}}





    #print(open('/home/james/Downloads/genbank_ecoli.txt').read())
    #print("PANDAS DATA")
    dataframe = pd.read_csv('/home/james/Downloads/genbank_ecoli.txt', sep=None)
    df = dataframe[['ftp_path', 'organism_name', '# assembly_accession', 'asm_name']]

    #drop all rows without serotype
    df2 = df[df['organism_name'].str.contains(('O\d+:H\d+'))]
    #print(df2)

    m = re.compile('O\d+')

    dict_list = []

    #set all organism names to serotype H_:O_
    for index, row in df2.iterrows():
        m = re.search('O\d+:H\d+', row['organism_name'])
        df2.set_value(index, 'organism_name', m.group(0))
        m_o = re.search('(?<=O)\d+', row['organism_name'])
        m_h = re.search('(?<=H)\d+', row['organism_name'])
        #dict = dict_list.append({'path': row['ftp_path'], 'H': m_h.group(0), 'O': m_o.group(0)})




    #drop duplicates in organism_name column

    df2 = df2.drop_duplicates('organism_name')

    for index, row in df2.iterrows():
        m = re.search('O\d+:H\d+', row['organism_name'])
        m_f = re.search('(?<=/)GCA_.+', row['ftp_path'])
        df2.set_value(index, 'organism_name', m.group(0))
        #df2.set_value(index, 'ftp_path', row['ftp_path'] + '/' + row['# assembly_accession'] + '_' + row['asm_name'] + '_genomic.fna.gz')
        m_o = re.search('O\d+', row['organism_name'])
        m_h = re.search('H\d+', row['organism_name'])
        dict_list.append({'H': m_h.group(0), 'O': m_o.group(0), 'file' : m_f.group(0) + '_genomic.fna'})

    #prints the ftp_paths to a text file
    #df2['ftp_path'].to_csv('/home/james/Downloads/ftp_sero', sep='\t', index=False)


    dir = '/home/james/Downloads/ecoli_sero_tests/'


    j = 0
    for entry in dict_list:
        file = dict_list[j]['file']
        genome_header = src.genomeFunctions.get_fasta_header_from_file(dir + file)
        genome_name = src.genomeFunctions.get_genome_name(genome_header)
        dict_list[j]['name'] = genome_name
        dict_list[j]['output'] = testing_dict[genome_name]
        j = j + 1





    print(dict_list)
    print(df2)
    print("DONE")


    test_data = dict_list




    args = Namespace(csv = CONST_CSV, input = CONST_INPUT, minimumGenomes = CONST_MINGENOMES, percentIdentity =CONST_PERIDENT, percentLength = CONST_PERLEN, serotyper = CONST_SER,  virulenceFactors = CONST_VIR)
    for entry in test_data:
       # f = open('delete_these', 'a')
        print("ENTRY")
        print(entry)
        print(entry['name'])
        print({entry['name'] :entry['output']})
        output = entry['output']['serotype']
        print('-------------------------------------------------------------------------------------')
        if entry['output']['otype']['strength'] == 0:
           if 'wzy' in output :
                print(src.serotypePrediction.parse_serotype(output['wzy']['blast_record'], args))
                result = src.serotypePrediction.parse_serotype(output['wzy']['blast_record'], args)['serotype']['wzy']['antigen']
                comparison = entry['O']
                print (result)
                print (comparison)
               # if result != comparison :
                   # f.write(entry['file']+ '\n')
           elif 'wzx' in output :
             print(src.serotypePrediction.parse_serotype(output['wzx']['blast_record'], args))
             result = src.serotypePrediction.parse_serotype(output['wzx']['blast_record'], args)['serotype']['wzx']['antigen']
             comparison = entry['O']
             print (result)
             print (comparison)
            # if result != comparison:
                # f.write(entry['file'] + '\n')
        else :
            print('printing')
            result = src.serotypePrediction.predict_serotype({entry['name'] : entry['output']})[entry['name']]['otype']['ant_number']
            comparison = entry['O']
            print('results =')
            print(result)
            print(comparison)
           # if result != comparison:
              #  f.write(entry['file'] + '\n')

        if entry['output']['htype']['strength'] == 0:
            if 'fliC' in output:
                print(src.serotypePrediction.parse_serotype(output['fliC']['blast_record'], args))
                result = src.serotypePrediction.parse_serotype(output['fliC']['blast_record'], args)['serotype']['fliC']['antigen']
                comparison = entry['H']
                print (result)
                print (comparison)
               # if result != comparison :
                    #f.write(entry['file']+ '\n')
        else :
            print('printing')
            result = src.serotypePrediction.predict_serotype({entry['name'] : entry['output']})[entry['name']]['htype']['ant_number']
            comparison = entry['H']
            print(result)
            print(comparison)
            #if result != comparison:
               # f.write(entry['file'] + '\n')
        #f.close()

test_sero_pred()

{'JASO01': {'serotype': {'fliC': {'antigen': 'H25', 'blast_record': {'qseqid': 'fliC_311_AGSG01000116_H25', 'length': '1320', 'send': '541217', 'sframe': '1', 'sseqid': 'JASO01000012.1', 'qlen': '1332', 'pident': '96.97', 'sstart': '539898'}}, 'gnd': {'antigen': 'O26', 'blast_record': {'qseqid': 'gnd|1_O26', 'length': '641', 'send': '423158', 'sframe': '1', 'sseqid': 'JASO01000012.1', 'qlen': '643', 'pident': '95.32', 'sstart': '422518'}}, 'wzy': {'antigen': 'O145', 'blast_record': {'qseqid': '9__wzy__wzy-O145__374', 'length': '1189', 'send': '413346', 'sframe': '1', 'sseqid': 'JASO01000012.1', 'qlen': '1188', 'pident': '99.75', 'sstart': '412158'}}}}, 'default_factory': {}}