#!/usr/bin/env python

import argparse
import os
import re
from Bio import SeqIO
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
GENOMES = {}

def parseCommandLine():
  """
  Initalizing the two main commands of the command line for the project.
  - input: refers to the location of the file(s) that will be processed
  - database: refers to the location of the BLAST database that will be used

  :return parser.parse_args():
  """

  parser = argparse.ArgumentParser()

  parser.add_argument("input", help="Location of new file(s). Can be a single file or a directory")
  parser.add_argument("-out", "-output", help="", default="STDOUT")

  return parser.parse_args()


def getListGenomes(data):
  """
  Creating a list out of the files entered (where each file name is its absolute path). This creates a uniform
  format that works for both single files and directories.
  :param data:
  :return filesList:
  """

  filesList = []

  if os.path.isdir(data):
    print("Using files from " + data)

    for root, dirs, files in os.walk(data):
     for filename in files:
       filesList.append(os.path.join(root,filename))

  else:
    print("Using file " + data)
    filesList.append(os.path.abspath(data))

  return sorted(filesList)


def getGenomeName(recordID, filename):
  """
  Getting the name of the genome by hierarchy to store or search in the GENOMES dictionary in the later called methods.

  :param recordID:
  :param filename:
  :return genomeName:
  """

  recordID = str(recordID)

  if re.search('lcl\|([\w-]*)', recordID):
    match = re.search('lcl\|([\w-]*)', recordID)
    match = str(match.group())
    genome_name = match.split('|')[1]

  elif re.search('(^[a-zA-Z][a-zA-Z]\w{6}\.\d)',recordID):
    match = re.search('(\w{8}\.\d)',recordID)
    genome_name = str(match.group())

  elif re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID):
    match = re.search('(ref\|\w{2}_\w{6}|gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID)
    match = str(match.group())
    genome_name = match.split('|')[1]

  elif re.search('gi\|\d{8}', recordID):
    match = re.search('gi\|\d{8}', recordID)
    match = str(match.group())
    genome_name = match.split('|')[1]

  else:
    genome_name = filename

  return genome_name


def checkFiles(genomesList):
  """
  Creating a list containing only valid .fasta files out of the previously formed list (from the getListGenomes method).
  If the newly created list is empty, the program exits with a warning.
  This filters out the invalid files.

  :param: genomesList
  :return: newGenomesList
  """
  newGenomesList = []

  for file in genomesList:
    flag = 0
    for record in SeqIO.parse(file, "fasta"):
      match = re.search('(^[a-zA-Z]+)', str(record.seq))
      if not match:
        break
      else:
        flag = 1

    if flag>0:
      newGenomesList.append(file)

      filename = os.path.basename(file)
      filename = os.path.splitext(filename)
      genome_name = getGenomeName(list(SeqIO.parse(file,"fasta"))[0], filename)

      GENOMES[genome_name] = ''
    else:
      print("File " + file+ " is in invalid format")

  if not newGenomesList:
    print("No valid fasta files \n Exiting")
    exit(1)

  else:
    return sorted(newGenomesList)


def initializeDB():
  """
  Generating the database if it was not already generated. The database can be found in the tmp folder.

  :return int 0 or 1: 0 being that the database was created successfully (or already existed).
  """

  REL_DIR = SCRIPT_DIRECTORY + '../temp/'

  if os.path.isfile(REL_DIR + 'ECTyperDB.nin'):
    print "Database already exists"
    return 0
  else:
    return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"])


def runBlastQuery(genomesList):
  """
  Generating the .xml files containing the genomes that were queried by using the command line.

  :param genomesList:
  :return resultList:
  """

  REL_DIR = SCRIPT_DIRECTORY + '../temp/'
  resultsList = []

  print "Searching the database..."

  for file in genomesList:
    filename = os.path.basename(file)
    filename = os.path.splitext(filename)

    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=file, db= REL_DIR + "ECTyperDB", outfmt=5, out= REL_DIR + filename[0] +".xml")

    stdout, stderr = blastn_cline()
    resultsList.append(REL_DIR + filename[0] + ".xml")

  print("Generated " + str(len(resultsList)) + " .xml file(s)")
  return sorted(resultsList)


def parseResults(resultsList):
  """
  Searching the result list to find the hsps necessary to identify the sequences.
  The method stores them in a dictionary that will be stored in GENOMES dictionary.

  :param resultsList:
  :return alignmentsDict:
  """
  alignmentsDict = {}

  for result in resultsList:
    result_handle = open(result)
    blast_records = NCBIXML.parse(result_handle)
    filename = os.path.basename(result)
    filename = os.path.splitext(filename)

    for blast_record in blast_records:
      genome_name = getGenomeName(blast_record.query, filename)
      for alignment in blast_record.alignments:
        alignmentsDict[alignment.title] = alignment.hsps[0]
      GENOMES[genome_name] = alignmentsDict

  print("Found " + str(len(alignmentsDict)) + " objects")
  return alignmentsDict


def findPerfectMatches(alignmentsDict):
  """
  Identifying the identical matches and the ones that will need a prediction.
  The method stores them in different dictionaries.

  :param alignmentsDict:
  :return identicalDict, predictionDict:
  """
  identicalDict = {}
  predictionDict = {}

  for title,hsp in alignmentsDict.iteritems():
    if len(hsp.query) == hsp.positives:
      identicalDict[title] = hsp
      #print hsp
    else:
      predictionDict[title] = hsp

  print identicalDict
  return identicalDict, predictionDict



if __name__=='__main__':
  args = parseCommandLine()
  roughGenomesList = getListGenomes(args.input)
  genomesList = checkFiles(roughGenomesList)

  if initializeDB() == 0:
    resultsList = runBlastQuery(genomesList)
    alignmentsDict = parseResults(resultsList)
    identicalDict, predictionDict = findPerfectMatches(alignmentsDict)

  else:
    print("Oops something happened")