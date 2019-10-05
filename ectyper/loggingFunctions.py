import logging

def create_logger():
    """
    Create the logger for ectyper

    :return: The root logger for the program
    """

    log = logging.getLogger('ectyper')
    formatter = logging.Formatter(
        '%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
    log.setLevel(logging.DEBUG)

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setFormatter(formatter)


    #console.setLevel(logging.INFO)
    log.addHandler(console)

    return log


