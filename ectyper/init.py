from ectyper import (definitions, speciesIdentification)
import logging
logging.basicConfig(level=logging.DEBUG)




def main():
    logging.info("Initializing Species ID database ...")
    speciesIdentification.get_species_mash(definitions.SPECIES_ID_SKETCH)
    logging.info("Done")