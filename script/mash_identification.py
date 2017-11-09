import os
import sys
import tempfile

import pandas as pd
from Bio import SeqIO
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
from ectyper import speciesIdentification, loggingFunctions
import subprocess
from tqdm import tqdm
import logging

LOG = logging.getLogger()

def main():
    loggingFunctions.initialize_logging('mash_identification.log')
    genome_dir = '/home/sam/Projects/MoreSerotype/temp/genomes'
    genomes = os.listdir(genome_dir)
    species_results = []
    for genome in tqdm(genomes):
        genome_file = os.path.join(genome_dir, genome)
        if os.path.getsize(genome_file) < 1000:
            LOG.info("%s is too small.", genome_file)
            continue
        with tempfile.TemporaryDirectory() as temp_dir:
            new_file = os.path.join(temp_dir, genome)
            new_fh = open(new_file,'w')
            header = '> %s\n' %genome
            new_fh.write(header)
            for record in SeqIO.parse(genome_file, 'fasta'):
                new_fh.write(str(record.seq))
                new_fh.write('nnnnnnnnnnnnnnnnnnnn')
            new_fh.close()
            genome_name = os.path.basename(genome)
            try:
                species = speciesIdentification.get_species(new_file)
            except subprocess.CalledProcessError as err:
                LOG.debug(err.stderr, err.stdout)
                species = ''
            new_entry = {
                'genome name':genome_name,
                'species':species
            }
            species_results.append(new_entry)
    df = pd.DataFrame(species_results)
    df.to_csv('species_results.csv')

if __name__ == '__main__':
    main()