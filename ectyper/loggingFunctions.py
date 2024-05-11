import logging, os
from ectyper import commandLineOptions

def create_logger():
    """
    Create the logger for ectyper

    :return: The root logger for the program
    """

    log = logging.getLogger('ectyper')
    formatter = logging.Formatter(
        '%(asctime)s %(name)s:%(lineno)d %(levelname)-5s %(message)s')
    
    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    log.addHandler(console)

    return log


