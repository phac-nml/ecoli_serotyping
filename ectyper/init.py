import argparse
from ectyper import (definitions, speciesIdentification)
import logging
logging.basicConfig(level=logging.DEBUG)





def main():
    logging.info("Initializing Species ID database ...")
    parser = argparse.ArgumentParser(
        description='Species ID database ectyper initializer'
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force species ID database initialization"
    )
    args = parser.parse_args()

    speciesIdentification.get_species_mash(definitions.SPECIES_ID_SKETCH, args.force)
    logging.info("Done")