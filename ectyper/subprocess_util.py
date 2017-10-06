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
    :param
        cmd (str): cmd command
    :return
        stdout(str)
    '''
    start_time = timeit.default_timer()
    LOG.info("Running: %s", cmd)
    comp_proc = subprocess.run(
        cmd,
        shell=is_shell,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stderr = comp_proc.stderr
    stdout = comp_proc.stdout
    elapsed_time = timeit.default_timer() - start_time
    if comp_proc.returncode == 0:
        LOG.info("Subprocess finish successfully in %0.3f sec.", elapsed_time)
        LOG.debug(stderr)
        return stdout
    LOG.fatal("Subprocess terminated with error in %d0.3f sec", elapsed_time)
    LOG.debug(stderr)
    exit(1)
