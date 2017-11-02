#!/usr/bin/env python

"""
    Set up the logging
"""

import logging
import os
from ectyper import definitions

LOG = logging.getLogger(__name__)

def initialize_logging(log_file=None):
    '''
    setup logging to stream INFO log to console
    and ALL log to log_file
    :param
        log_file (str): log filename [optional]
    :return
        Log filename
    '''
    if log_file is None:
        log_file = os.path.join(definitions.WORKPLACE_DIR, 'default.log')


    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s.%(msecs)02d %(name)-25s %(levelname)-8s %(message)s',
                        datefmt='%y-%m-%d %H:%M:%S',
                        filename=log_file,
                        filemode='w')
    # If logger already exist, return the existing logger
    if len(logging.getLogger().handlers)>1:
        return log_file
    # add the handler to the root logger
    ch_handler = logging.StreamHandler()
    ch_handler.setLevel('INFO')
    ch_handler.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logging.getLogger('').addHandler(ch_handler)

    # Now, we can log to the root logger, or any other logger. First the root...
    LOG.info('Logger initialized. Debug log stored at %s', log_file)
    return log_file