#!/usr/bin/env python

from ecvalidatingfiles import *
from Bio import SeqIO

SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"
TITLES = {}

### This module is subject to change, depending on the input format of the sequence names ###

def parseCommandLine():
    """
    Initalizing the two main commands of the command line for the project.
    - names: location of the file(s) containing the desired names to format
    - sequences: location of the files containing the sequences to add to the database

    :return parser.parse_args():
    """

    parser = argparse.ArgumentParser()
    parser.add_argument("sequences", help="Location of the file(s) of the sequences to add to the database. Can be a single file or a directory.")

    return parser.parse_args()


def formatName(title, count):
    """
    Formatting the sequence name. If a sequence shares multiple types, a title will be assigned to each.

    :param title: String to be searched
    :param count: Integer that will be modified if the sequence is a gnd sequence.
    :return count:
    """
    old_title = title

    if title not in TITLES:
        TITLES[old_title] = []

    #Already normally formatted
    match = re.search(r'(-O\d+|-H\d+)', title)
    if match != None:
        TITLES[old_title].append(title)
        return count

    #gnd format
    match = re.search('serogroup=O\w+', title)
    if match != None:
        match = str(match.group())
        match = match.split('=')[1]
        if match != 'OUT':
            new_title = 'gnd|' + str(count) +  "_"  + str(match)
        else:
            new_title = 'gnd|' + str(count) +  "_OR"
        TITLES[old_title].append(new_title)
        count +=1
    else:
        #Danish database format
        matches = re.findall(r'(O\d+\w*|H\d+\w*)', title)
        if matches != None:
            for match in matches:
                new_title  = re.sub(r'(O\d+\w*|H\d+\w*|[/]*)', '', title)
                match = str(match)
                TITLES[old_title].append(new_title+match)

    return count



def addSeqToDB(filesList):
    """
    Filtering through the sequences to format their names and then add them to the database.

    :param filesList:
    """

    count =1
    for record in SeqIO.parse(SCRIPT_DIRECTORY + "../Data/gnd_sequences.fasta", "fasta"):
        count +=1

    with open(SCRIPT_DIRECTORY + "../Data/EcOH.fasta", 'a+') as handle:
        for file in filesList:
            for record in SeqIO.parse(file, "fasta"):
                temp_count = formatName(record.description, count)

                if record.description in TITLES:
                    for title in TITLES[record.description]:
                        record.id = title
                        record.description = title
                        record.name = title
                        SeqIO.write(record, handle, "fasta")

                        if temp_count > count:
                            SeqIO.write(record, open(SCRIPT_DIRECTORY + "../Data/gnd_sequences.fasta", 'a+'), "fasta")
                            count = temp_count

    subprocess.call(["/usr/bin/makeblastdb", "-in", SCRIPT_DIRECTORY + "../Data/EcOH.fasta ", "-dbtype", "nucl", "-title", "ECTyperDB", "-out", SCRIPT_DIRECTORY+ "../temp/ECTyperDB"])


if __name__=='__main__':

    args = parseCommandLine()
    seqsList = getFilesList(args.sequences)
    addSeqToDB(seqsList)




