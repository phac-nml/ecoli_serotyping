#!/usr/bin/env python

"""
    Set up the logging
"""

import logging
import tempfile
import os

LOG = logging.getLogger(__name__)

def initialize_logging(log_file='default.log'):
    '''
    setup logging to stream INFO log to console
    and ALL log to log_file
    :param
        log_file (str): log filename [optional]
    :return
        Log filename
    '''
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s.%(msecs)02d %(name)-20s %(levelname)-8s %(message)s',
                        datefmt='%y-%m-%d %H:%M:%S',
                        filename=log_file,
                        filemode='w')
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    # set a format which is simpler for console use
    formatter = logging.Formatter(
        fmt='%(name)-20s: %(levelname)-8s %(message)s',
        datefmt='%H:%M:%S'
    )
    # tell the handler to use this format
    console.setFormatter(formatter)
    # add the handler to the root logger
    logging.getLogger('').addHandler(console)

    # Now, we can log to the root logger, or any other logger. First the root...
    LOG.info('Logger initialized.')
    return log_file