#!/usr/bin/env python

from flask import *



def toJSON(data):

    stylesheet_url = url_for('static',filename='css/controller.css')
    js_url = url_for('static', filename='js/downloadresults.js')
    html_info = "<!DOCTYPE html> \
                 <html lang='en'> \
                 <head> \
                    <meta charset='UTF-8'> \
                    <title>Results</title> \
                    <link rel='stylesheet' type='text/css' href=" + str(stylesheet_url) + ">" \
                    "<script type='text/javascript' src='" + js_url + "'></script>\
                 </head>\
                 <body>"\
                "<div id='result-button'>" \
                    "<input type='button' class='button result-button' onclick='" + \
                    'location.href="/superphy/controller";' + \
                    "' value='Return to main page'/>" \
                    "<button type='button' id='download-button-json' class='button result-button'>Download the results" \
                    "</button></div>" \
                    "<div><h2 class='results-caption'>RESULTS</h2><p>"


    return html_info + json.dumps(data) + "</p></div></body></html>"


def toTable(data):
    """
    Rendering the results dictionary in an HTML table.

    :param data: Results dictionary.
    :param verbose: Boolean defining the type of information the user desires.
    :return: String containing the HTML page containing the table.
    """

    stylesheet_url = url_for('static',filename='css/controller.css')
    js_url = url_for('static', filename='js/results.js')
    html_info = "<!DOCTYPE html> \
                 <html lang='en'> \
                 <head> \
                    <meta charset='UTF-8'> \
                    <title>Results</title> \
                    <link rel='stylesheet' type='text/css' href=" + str(stylesheet_url) + ">" \
                    "<script type='text/javascript' src='" + js_url + "'></script>  \
                 </head>\
                 <body>"
    result_table = "<div id='results-div'><table class=results><caption class='results-caption'>RESULTS" \
                   "</caption>" \
                   "<thead class='results-head'><tr>" \
                   "<th>Genome</th>" \
                   "<th>O Type</th>" \
                   "<th>H Type</th>" \
                   "</tr></thead><tbody>"

    for genome_name in data.keys():
        result_genome = str(genome_name)
        result_otype = str(data[genome_name]['otype']).replace(";", '<br>')
        result_htype = str(data[genome_name]['htype']).replace(";", '<br>')

        result_table+= "<tr><td id='result-genome'>" + result_genome + "</td><td>" + result_otype + "</td><td>" + \
                       result_htype + "</td></tr>"

    result_table+= "</tbody></table></div>"


    return_button = "<div id='result-button'>" \
                    "<input type='button' class='button result-button' onclick='" + \
                    'location.href="/upload";' + \
                    "' value='Return to main page'/>" \
                    "<button type='button' id='download-button' class='button result-button'>Download the results" \
                    "</button></div>"

    return html_info + result_table + return_button + "</body></html>"