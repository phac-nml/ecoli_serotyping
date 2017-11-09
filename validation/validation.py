'''
Validation script to run ectyper on all enterobase genome
'''
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
import tempfile
import shutil
import random

from ectyper import ectyper

if __name__=='__main__':
    # temp_dir = os.path.join(tempfile.gettempdir(), 'sample')
    # os.mkdir(temp_dir)
    # print(temp_dir)
    genome_dir = '/home/sam/genomes'
    # samples = random.sample(os.listdir(genome_dir), 10)
    # for basename in samples:
    #     src_file = os.path.join(genome_dir, basename)
    #     dst_file = os.path.join(temp_dir, basename)
    #     try:
    #         shutil.copy(src_file, dst_file)
    #     # eg. src and dest are the same file
    #     except shutil.Error as e:
    #         print('Error: %s' % e)
    #     # eg. source or destination doesn't exist
    #     except IOError as e:
    #         print('Error: %s' % e.strerror)
    args = ['-i', genome_dir, '-out', 'new.json']
    sys.argv[1:] = args
    ectyper.run_program()
