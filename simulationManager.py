######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: simulationManager.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 21 July 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import model
from dataManager import *
import numpy as np
import matplotlib.pyplot as plt


######################################################################################################################################################################

def simulate():
    pass

def plotFidelity(solutionPath, useSavedPlot=False, saveToFile=False):
    solutionDict = getFromjson(fileName=solutionPath)
    x = (solutionDict['theta'], solutionDict['delta'], solutionDict['omegaPhi'], solutionDict['modulationTime'])
    
    circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
    circuitData = {}
    for key in circuitDataKeys:
        circuitData[key] = [2*np.pi*item for item in solutionDict[key]]

    if solutionDict['gateType'] == 'iSWAP':
        iSWAP = True
    else:
        iSWAP = False
    if solutionDict['gateType'] == 'CZ':
        CZ = True
    else:
        CZ = False
    
    if useSavedPlot:
        F = solutionDict['fidelitiesAtTimes']
        times = solutionDict['times']
    else:
        indices = np.linspace(-116, -1, 116).astype(int)
        F, times = model.getGateFidelity(x, N=4, iSWAP=iSWAP, CZ=CZ, tIndices=indices, circuitData=circuitData)

    plt.figure(figsize=(8, 7))
    plt.plot(times, F)
    plt.plot([x[-1], x[-1]], [0, 1], 'r--')
    plt.grid()
    plt.ylim([0.99, 1])
    plt.xlim([times[0], times[-1]])
    plt.legend(["Fidelitet", "$t_{MOD}$"], fontsize=19, loc="lower right")
    #plt.title("Grindfidelitet kring $t_{MOD}$", fontsize=17)
    plt.xlabel("Tid efter grindstart [ns]", fontsize=26)
    plt.ylabel("Fidelitet", fontsize=26)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    plt.show()

    if (saveToFile == True):
        solutionDict['times'] = times
        solutionDict['fidelitiesAtTimes'] = F
        solutionDict['gateFidelity'] = F[-76]
        dumpTojson(solutionDict,solutionPath)