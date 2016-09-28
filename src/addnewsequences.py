#!/usr/bin/env python

from ecvalidatingfiles import *
from Bio import SeqIO

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
TITLES = {}


def parseCommandLine():
    """
    Initalizing the two main commands of the command line for the project.
    - names: location of the file(s) containing the desired names to format
    - sequences: location of the files containing the sequences to add to the database

    :return parser.parse_args():
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("names", help="Location of the file(s) of the names to format. Can be a single file or a directory.")
    parser.add_argument("sequences", help="Location of the file(s) of the sequences to add to the database. Can be a single file or a directory.")

    return parser.parse_args()


def formatNames(filesList):
    """
    Formatting the sequence names so that they are separated by type. If a sequence shares multiple types, a title will
    be assigned to each.

    :param filesList:
    """

    filesDict = {}

    for file in filesList:
        filesDict[file] = [line.rstrip('\n') for line in open(file)]


    for file, lines in filesDict.iteritems():
        for line in lines:
            matches = re.findall('(O\d+\w*|H\d+\w*)', line)
            if matches != None:
                oldLine = line
                TITLES[oldLine] = []
                for match in matches:
                    line  = re.sub(r'(O\d+\w*|H\d+\w*|[/]*)', '', line)
                    match = str(match)
                    TITLES[oldLine].append(line+match)




def addSeqToDB(filesList):
    """
    Filtering through the sequences to find the ones with titles that match those of the names list.
    If matching, the sequence is added to the database.

    :param filesList:
    """


    with open(SCRIPT_DIRECTORY + "../Data/EcOH.fasta", 'a+') as handle:
        for file in filesList:
            for record in SeqIO.parse(file, "fasta"):

                if record.id in TITLES:
                    for title in TITLES[record.id]:
                        record.id = title
                        record.description = title
                        record.name = title
                        SeqIO.write(record, handle, "fasta")

    subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", SCRIPT_DIRECTORY+ "../temp/ECTyperDB"])


if __name__=='__main__':

    args = parseCommandLine()
    namesList = getFilesList(args.names)
    seqsList = getFilesList(args.sequences)
    formatNames(namesList)
    addSeqToDB(seqsList)



