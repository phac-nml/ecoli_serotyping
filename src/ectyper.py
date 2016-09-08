#!/usr/bin/env python

import argparse
import os
from Bio import SeqIO



def parseCommandLine():

  parser = argparse.ArgumentParser()

  parser.add_argument("input", help="Location of new file(s).")
  parser.add_argument("database", help="Location of the BLAST database.")

  return parser.parse_args()


def getListGenomes(input_data):

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
  #call functions
  args = parseCommandLine()
  roughListGenomes = getListGenomes(args.input)
  listGenomes = checkFiles(roughListGenomes)
  print(listGenomes)