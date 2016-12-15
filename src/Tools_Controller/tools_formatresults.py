#!/usr/bin/env python
from collections import OrderedDict


def tosimpleString(serotype_dict):

    end_string = ''

    for key, value in serotype_dict.iteritems():
        if key == 'RESULT':
            end_string = str(key) + ": " + str(value) + "; " + end_string
        else:
            end_string+= str(key) + ": " + str(value) + "; "

    return end_string


def toTableList(data, verbose):
    """
    Obtaining the necessary information to build the HTML table containing the results.
    First, going through the result dictionary to gather all the genes/antigen for each tool and storing them in
    another dictionary. Second going through the result dictionary to sort it so that the HTML table headers
    have corresponding HTML table data.

    :param data: Results dictionary.
    :param verbose: amount of information about the serotype desired (1= full, 0= serotype only)
    :return: list of the main headers, the sorted result dictionary and the subheaders dictionary.
    """

    main_headers = {}
    temp_dict = {}
    result_Dict = {}
    for genome_name, results_info in data.iteritems():
        result_Dict[genome_name] = {}

        for result_name, genes_info in results_info.iteritems():
            result_Dict[genome_name][result_name] = []

            if result_name not in main_headers:
                main_headers[result_name] = 0
                temp_dict[result_name] = []

            for gene_name in genes_info.keys():
                if gene_name not in temp_dict[result_name]:
                    main_headers[result_name] += 1
                    temp_dict[result_name].append(str(gene_name))
                if result_name == 'Serotype' and isinstance(genes_info[gene_name], dict) and verbose=='1':
                    temp_str = tosimpleString(data[genome_name][result_name][gene_name])
                    data[genome_name][result_name][gene_name] = ''
                    data[genome_name][result_name][gene_name] = temp_str

            temp_dict[result_name] = sorted(temp_dict[result_name], key=str.lower)

    for genome_name, results_info in data.iteritems():
        for result_name, genes_info in results_info.iteritems():
            for gene_name in temp_dict[result_name]:
                if gene_name not in genes_info:
                    result_Dict[genome_name][result_name].append('0')
                else:
                    result_Dict[genome_name][result_name].append(str(genes_info[gene_name]).replace(';','\n'))

    return [main_headers, result_Dict, temp_dict]


def getCSV(data, verbose):
    """
    Generating a CSV string from the results dictionary to be downloaded afterwards.

    :param data: Results dictionary.
    :return: string containing the headers, the subheaders and the data from the results dictionary.
    """

    table_list = toTableList(data, verbose)

    content = table_list[1]
    main_headers = table_list[0]
    headers = table_list[2]

    headers_str = ""

    for result_name, appearances in main_headers.iteritems():
        if result_name == 'Serotype':
            temp_str = ',' + result_name
            for x in range(1, appearances):
                temp_str += ","
            headers_str = temp_str + headers_str
        else:
            headers_str += ',' + result_name
            for x in range(1, appearances):
                headers_str += ","

    headers_str += '\r\n'

    subheaders_str = ''
    for  result_name, genes in headers.iteritems():
        if result_name == 'Serotype':
            temp_str = ''
            for gene_name in genes:
                temp_str += ',' + gene_name
            subheaders_str = temp_str + subheaders_str
        else:
            for gene_name in genes:
                subheaders_str += ',' + gene_name


    subheaders_str = "Genomes" + subheaders_str + '\r\n'

    content_str = ''
    for genome_name, result_info in content.iteritems():
        for result_name, gene_info in result_info.iteritems():
            if result_name == 'Serotype':
                temp_str = ''
                for gene_name in gene_info:
                    gene_str = gene_name.replace('\n', ';')
                    temp_str +=',' + gene_str
                content_str = temp_str + content_str
            else:
                for gene_name in gene_info:
                    gene_str = gene_name.replace('\n', ';')
                    content_str +=',' + gene_str
        content_str = genome_name + content_str + '\r\n'

    return headers_str + subheaders_str + content_str












