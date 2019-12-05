#!/usr/bin/env python

import logging
import subprocess
import timeit

LOG = logging.getLogger(__name__)


def run_subprocess(cmd, input_data=None, un=False, ignorereturncode=False):
    """
    Run cmd command on subprocess

    Args:
        cmd (str): cmd command

    Returns:
        stdout (str): The stdout of cmd
    """

    start_time = timeit.default_timer()
    comp_proc = subprocess.run(
        cmd,
        shell=False,
        input = input_data,
        check=False,
        universal_newlines=un,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )


    if comp_proc.returncode == 0 or ignorereturncode == True:
        elapsed_time = timeit.default_timer() - start_time
        LOG.debug("Subprocess {} finished successfully in {:0.3f} sec.".format(cmd, elapsed_time))
        return comp_proc
    else:
        LOG.error("Error in subprocess. The following command failed: {}".format(cmd))
        LOG.error("Subprocess failed with error: \"{}\"".format(comp_proc.stderr.decode("utf-8")))
        LOG.critical("ectyper has stopped")
        exit("subprocess failure")
