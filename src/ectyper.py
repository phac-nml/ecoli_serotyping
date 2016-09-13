#!/usr/bin/env python

import argparse
import os
import re
from Bio import SeqIO
import subprocess
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"


def parseCommandLine():
  """
  Initalizing the two main commands of the command line for the project.
  - input: refers to the location of the file(s) that will be processed
  - database: refers to the location of the BLAST database that will be used
  :return parser.parse_args():
  """

  parser = argparse.ArgumentParser()

  parser.add_argument("input", help="Location of new file(s). Can be a single file or a directory")
 # parser.add_argument("database", help="Location of the BLAST database.")
  parser.add_argument("-out", "-output", help="", default="STDOUT")

  return parser.parse_args()


def getListGenomes(input_data):
  """
  Creating a list out of the files entered (where each file name is its absolute path). This creates a uniform
  format that works for both single files and directories.
  :param input_data:
  :return files:
  """

  filesList = []

  if os.path.isdir(input_data):
    print("Using files from " + input_data)

    for root, dirs, files in os.walk(input_data):
     for filename in files:
       filesList.append(os.path.join(root,filename))

  else:
    print("Using file " + input_data)
    filesList.append(os.path.abspath(input_data))

  return sorted(filesList)

def getGenomeName(file):

  record = SeqIO.parse(file,"fasta")[0]
  recordId = str(record.id)
  genomeName = ''

  if re.search('lcl\|([\w-]*)', recordId):
    genomeName =
  elif re.search('',recordId):
    genomeName =
  elif re.search('',recordId):
    genomeName =
  elif re.search('',recordId):
    genomeName =
  elif re.search('', recordId):
    genomeName =
  elif re.search('', recordId):
    genomeName =
  elif re.search('', recordId):
    genomeName =
  else:
    genomeName =



def checkFiles(listGenomes):
  """
  Creating a list containig only valid .fasta files out of the previously formed list (from the getListGenomes method).
  If the newly created list is empty, the program exits with a warning.
  This filters out the invalid files.
  :param listGenomes:
  :return newListGenomes:
  """
  newListGenomes = []

  for filename in listGenomes:
    flag = 0
    for record in SeqIO.parse(filename, "fasta"):
      match = re.search('(^[a-zA-Z]+)', str(record.seq))
      if not match:
        break
      else:
        flag = 1

    if flag>0:
      newListGenomes.append(filename)
      genomeName = getGenomeName(filename)
    else:
      print("File " + filename+ " is in invalid format")

  if not newListGenomes:
    print("No valid fasta files \n Exiting")
    exit(1)

  else:
    return sorted(newListGenomes)


def initializeDB():

  REL_DIR = SCRIPT_DIRECTORY + '../temp/'

  if os.path.isfile(REL_DIR + 'ECTyperDB.nin'):
    print "Database already exists"
    return 0
  else:
    return subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", REL_DIR + "ECTyperDB"])


def runBlast(listGenomes):

  REL_DIR = SCRIPT_DIRECTORY + '../temp/'
  resultList = []

  print "Searching the database..."

  for file in listGenomes:
    filename = os.path.basename(file)
    filename = os.path.splitext(filename)

    blastn_cline = NcbiblastnCommandline(cmd="blastn", query=file, db= REL_DIR + "ECTyperDB", outfmt=5, out= REL_DIR + filename[0] +".xml")

    stdout, stderr = blastn_cline()
    resultList.append(REL_DIR+filename[0] + ".xml")

  print("Generated " + str(len(resultList)) + " .xml file(s)")
  return resultList


def parseResults(results):

  hspList = []

  for result in results:
    result_handle = open(result)
    blast_records = NCBIXML.parse(result_handle)

    for blast_record in blast_records:
      for alignment in blast_record.alignments:
         hspList.append(alignment.hsps[0])
         #print("HSP: \n" + str(alignment.hsps[0]))

  print("Found " + str(len(hspList)) + " objects")
  return hspList


def findPerfectMatches(hspList):

  identicalList = []
  predictionList = []

  for hsp in hspList:
    if len(hsp.query) == hsp.positives:
      identicalList.append(hsp)
      print hsp
    else:
      predictionList.append(hsp)

  return identicalList, predictionList



if __name__=='__main__':
  args = parseCommandLine()
  roughListGenomes = getListGenomes(args.input)
  listGenomes = checkFiles(roughListGenomes)

  if initializeDB() == 0:
    results = runBlast(listGenomes)
    hspList = parseResults(results)
    identicalList, predictionList = findPerfectMatches(hspList)

  else:
    print("Oops something happened")