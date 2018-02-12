'''
Utilities
'''
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
import logging
import timeit
import os
import sys

if os.name == 'posix' and sys.version_info[0] < 3:
    import subprocess32 as subprocess
else:
    import subprocess

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
