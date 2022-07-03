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

# Date created: 21 April 2022

# Copyright 2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports.
import json
import os

# Local Python imports.
import graph_generation.fidelity_plot as fidelityGraph

######################################################################################################################################################################
# Functions to handle JSON objects

def getFromjson(path):
    with open(path, 'r') as file:
        dict = json.load(file)
    return dict


######################################################################################################################################################################
# Interface functions

def plotFidelity(solutionPath, saveLocationPath, saveGraph=True):

    if not os.path.isfile(solutionPath):
        raise Exception("Provided path does not point at a file!")

    solution = getFromjson(solutionPath)

    try:
        times = solution["times"]
        fidelities = solution["fidelitiesAtTimes"]
        modulationTime = solution["modulationTime"]
    except (Exception):
        print("Data corrupted or not available.")
        return

    try:
        fidelityGraph.generate(times, fidelities, modulationTime, saveLocationPath, saveGraph=saveGraph)
    except (Exception):
        print("Error occured during graph generation.")
        return


def plotPopulation(solutionPath, saveLocationPath, saveGraph=True):
    
    if not os.path.isfile(solutionPath):
        raise Exception("Provided path does not point at a file!")

    solution = getFromjson(solutionPath)

######################################################################################################################################################################
