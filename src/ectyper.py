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
  :param input_data:
  :return files:
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
    genomeName = match.split('|')[1]

  elif re.search('(^[a-zA-Z][a-zA-Z]\w{6}\.\d)',recordID):
    match = re.search('(\w{8}\.\d)',recordID)
    genomeName = str(match.group())

  elif re.search('ref\|\w{2}_\w{8}|(gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID):
    match = re.search('ref\|\w{2}_\w{8}|(gb\|\w{8}|emb\|\w{8}|dbj\|\w{8})',recordID)
    match = str(match.group())
    genomeName = match.split('|')[1]

  elif re.search('gi\|\d{8}', recordID):
    match = re.search('gi\|\d{8}', recordID)
    match = str(match.group())
    genomeName = match.split('|')[1]

  else:
    genomeName = filename

  return genomeName


def checkFiles(listGenomes):
  """
  Creating a list containig only valid .fasta files out of the previously formed list (from the getListGenomes method).
  If the newly created list is empty, the program exits with a warning.
  This filters out the invalid files.

  :param: listGenomes
  :return: newListGenomes
  """
  newListGenomes = []

  for file in listGenomes:
    flag = 0
    for record in SeqIO.parse(file, "fasta"):
      match = re.search('(^[a-zA-Z]+)', str(record.seq))
      if not match:
        break
      else:
        flag = 1

    if flag>0:
      newListGenomes.append(file)

      fileName = os.path.basename(file)
      fileName = os.path.splitext(fileName)
      genomeName = getGenomeName(list(SeqIO.parse(file,"fasta"))[0], fileName)

      GENOMES[genomeName] = ''
    else:
      print("File " + file+ " is in invalid format")

  if not newListGenomes:
    print("No valid fasta files \n Exiting")
    exit(1)

  else:
    return sorted(newListGenomes)


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


def runBlast(listGenomes):
  """
  Generating the .xml files containing the genomes that were queried by using the command line.

  :param listGenomes:
  :return resultList:
  """

  REL_DIR = SCRIPT_DIRECTORY + '../temp/'
  resultList = []

  print "Searching the database..."

  for file in listGenomes:
    fileName = os.path.basename(file)
    fileName = os.path.splitext(fileName)

    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=file, db= REL_DIR + "ECTyperDB", outfmt=5, out= REL_DIR + fileName[0] +".xml")

    stdout, stderr = blastn_cline()
    resultList.append(REL_DIR+fileName[0] + ".xml")

  print("Generated " + str(len(resultList)) + " .xml file(s)")
  return resultList


def parseResults(results):

  alignments = {}

  for result in results:
    result_handle = open(result)
    blast_records = NCBIXML.parse(result_handle)
    fileName = os.path.basename(result)
    fileName = os.path.splitext(fileName)

    for blast_record in blast_records:
      genomeName = getGenomeName(blast_record.query, fileName)
      for alignment in blast_record.alignments:
        #hspList.append(alignment.hsps[0])
         #print("HSP: \n" + str(alignment.hsps[0]))
        alignments[alignment.title] = alignment.hsps[0]
      GENOMES[genomeName] = alignments

  print("Found " + str(len(alignments)) + " objects")
  return alignments

def findPerfectMatches(alignments):

  identicalDict = {}
  predictionDict = {}

  for title,hsp in alignments.iteritems():
    if len(hsp.query) == hsp.positives:
      identicalDict[title] = hsp
      #print hsp
    else:
      predictionDict[title] = hsp

  print identicalDict
  return identicalDict, predictionDict



if __name__=='__main__':
  args = parseCommandLine()
  roughListGenomes = getListGenomes(args.input)
  listGenomes = checkFiles(roughListGenomes)

  if initializeDB() == 0:
    results = runBlast(listGenomes)
    alignments = parseResults(results)
    identicalDict, predictionDict = findPerfectMatches(alignments)

  else:
    print("Oops something happened")