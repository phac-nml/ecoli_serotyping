#!/usr/bin/env python

import os
SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__)) + "/"

def createDirs():

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/Serotyping_Database/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/Serotyping_Database/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/databases/VF_Database/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/databases/VF_Database/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Results/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/Results/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/xml/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/xml/')

    if not os.path.isdir(SCRIPT_DIRECTORY + '../temp/Uploads/'):
        os.mkdir(SCRIPT_DIRECTORY + '../temp/Uploads/')
