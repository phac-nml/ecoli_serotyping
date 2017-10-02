'''
Utilities
'''
import logging, subprocess

LOG = logging.getLogger(__name__)

def run_subprocess(cmd, shell=False):
    '''
    Run cmd command on subprocess
    :param
        cmd (str): cmd command
    :return
        output(str)
    '''
    LOG.info("Running: %s", cmd)
    comp_proc = subprocess.run(
        cmd,
        shell=shell,
        universal_newlines=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    output = '/n'.join([
        str(comp_proc.returncode),
        comp_proc.stdout,
        comp_proc.stderr
    ])
    if comp_proc.returncode == 0:
        LOG.info("Subprocess finish successfully")
        LOG.debug(output)
        return output
    LOG.fatal("Subprocess error")
    LOG.debug(output)
    exit(1)