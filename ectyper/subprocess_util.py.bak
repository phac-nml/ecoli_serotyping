'''
Utilities
'''
import logging
import subprocess
import timeit

LOG = logging.getLogger(__name__)

def run_subprocess(cmd, is_shell=False):
    '''
    Run cmd command on subprocess

    Args:
        cmd (str): cmd command

    Returns:
        stdout (str): The stdout of cmd
    '''
    start_time = timeit.default_timer()
    LOG.debug("Running: {0}".format(cmd))
    comp_proc = subprocess.run(
        cmd,
        shell=is_shell,
        universal_newlines=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stderr = comp_proc.stderr
    stdout = comp_proc.stdout
    elapsed_time = timeit.default_timer() - start_time
    LOG.debug("Subprocess finished successfully in {:0.3f} sec.".format(elapsed_time))
    return stdout
