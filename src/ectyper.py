#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO



def parseCommandLine():
  """
  Initalizing the two main commands of the command line for the project.
  - input: refers to the location of the file(s) that will be processed
  - database: refers to the location of the BLAST database that will be used
  :return parser.parse_args():
  """

  parser = argparse.ArgumentParser()

  parser.add_argument("input", help="Location of new file(s).")
  parser.add_argument("database", help="Location of the BLAST database.")

  return parser.parse_args()


def getListGenomes(input_data):
  """
  Creating a list out of the files entered (where each file name is its absolute path). This creates a uniform
  format that works for both single files and directories.
  :param input_data:
  :return files:
  """

  files = []

  if os.path.isdir(input_data):
    #print("Using files from " + input_data)
    tmp = os.listdir(input_data)

    fullPath = os.path.normpath(input_data)

    for filename in tmp:
      files.append(fullPath + "/" + filename)

    print files

  else:
    # print("Using file " + input_data)
    files.append(os.path.abspath(input_data))
    print files

  return files


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
    counter=0
    for seq_record in SeqIO.parse(filename, "fasta"):
      counter+=1
    if counter>0:
      newListGenomes.append(filename)
    else:
      print("File is in invalid format")

  if not newListGenomes:
    print("No valid fasta files")
    exit(1)

  else:
    return newListGenomes



if __name__=='__main__':
  args = parseCommandLine()
  roughListGenomes = getListGenomes(args.input)
  listGenomes = checkFiles(roughListGenomes)
  print(listGenomes)