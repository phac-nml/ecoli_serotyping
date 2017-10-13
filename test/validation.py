'''
Validation script to run ectyper on all enterobase genome
'''
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from ectyper import ectyper

if __name__=='__main__':
    input = '/home/sam/Projects/MoreSerotype/temp/genomes'
    args = ['-i', input]
    sys.argv[1:] = args
    ectyper.run_program()
