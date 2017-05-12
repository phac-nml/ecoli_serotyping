#!/usr/bin/env python

"""
Convert the output from ectyper to table format
"""

import logging
import collections

log = logging.getLogger(__name__)

def results_dict_to_table(list_of_files, list_of_genomes, results_dict):
    """
    Take the results dict, and converts it to output in table, structure.
    Builds up an array of output.
    
    :param results_dict: 
    :return: results_as_table
    """

    output_aoa = [["filename"
                  ,"genome"
                  ,"contig ID"
                  ,"analysis"
                  ,"hit"
                  ,"description"
                  ,"orientation"
                  ,"start"
                  ,"stop"]]

    for f, g in zip(list_of_files, list_of_genomes):

        if "otype" in results_dict[g]:
            serotype = results_dict[g]["otype"]["ant_number"] + ":" +  results_dict[g]["htype"]["ant_number"]
            serotype_line = [f, g, "NA","Serotype", serotype, "NA", "NA", "NA", "NA"]
            output_aoa.append(serotype_line)

        if "vf" not in results_dict[g]:
            log.warning("No VFs for genome: " + g)
            continue

        for vf in results_dict[g]["vf"]:
            analysis = "Virulence Factor"
            br = results_dict[g]["vf"][vf]["blast_record"]
            contig_id = br["sseqid"]
            hit = vf
            description = results_dict[g]["vf"][vf]["description"]
            orientation = br["sframe"]
            startbp = br["sstart"]
            stopbp = br["send"]

            next_line = [f
                         ,g
                         ,contig_id
                         ,analysis
                         ,hit
                         ,description
                         ,orientation
                         ,startbp
                         ,stopbp]

            output_aoa.append(next_line)

    return(output_aoa)

