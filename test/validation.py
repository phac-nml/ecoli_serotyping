'''
Validation script to run ectyper on all enterobase genome
'''
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from ectyper import ectyper

if __name__=='__main__':
    genome_dir = '/home/sam/Projects/MoreSerotype/temp/genomes'
    genomes_chunk = []
    for genome_file in os.listdir(genome_dir):
        if len(genomes_chunk) < 10:
            genomes_chunk.append(genome_file)
        else:
            # 
    args = ['-i', input]
    sys.argv[1:] = args
    ectyper.run_program()
