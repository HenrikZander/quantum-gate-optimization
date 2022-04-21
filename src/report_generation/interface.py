######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: interface.py

# Author(s): Henrik Zander

# Date created: 8 March 2022

# Copyright 2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports
import json
import os

# Local Python imports
import report_generation.report_two_qubit_gate as rprt2

######################################################################################################################################################################
# Functions to handle JSON objects

def getFromjson(path):
    with open(path, 'r') as file:
        dict = json.load(file)
    return dict


######################################################################################################################################################################
# Interface functions

def generateReport2(resultPath):
    if rprt2.numOfAvailableGraphs() != 9:
        raise Exception("Incorrect number of graphs.")
    
    if not os.path.isfile(resultPath):
        raise Exception("Parameter 'resultPath' is not a file!")
    
    result = getFromjson(resultPath)
    rprt2.generateReport(result)


######################################################################################################################################################################
# path = "C:\\Users\\Henrik Zander\\Desktop\\WACQT\\quantum-gate-optimization\\Results\\Qubit Pair 01\\Solutions\\210728_iSWAP_1.json"
# generateReport2(path)