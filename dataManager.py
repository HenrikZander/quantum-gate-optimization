######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: dataManager.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 19 July 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import os
import json
import numpy as np
from datetime import *
import copy

######################################################################################################################################################################
# Functions that handle data loading and dumping to json-files.


def getFromjson(fileName):
    with open(fileName, 'r') as jsonFile:
        jsonDict = json.load(jsonFile)
    return jsonDict


def dumpTojson(data, fileName):
    with open(fileName, 'w') as jsonFile:
        json.dump(data, jsonFile, ensure_ascii=False, indent=4)


def getCircuitData(solutionDict):
    circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
    circuitData = {}
    for key in circuitDataKeys:
        circuitData[key] = [2*np.pi*item for item in solutionDict[key]]
    return circuitData


def createSolName(ymd, gateType, solNumber):
    return ymd + "_" + gateType + "_" + str(solNumber)


def addNewSolution(x, gateType, N, solNumber=1, creationTime=datetime.today(), folder='Demo Circuit', circuitFile=None, circuitData=None, riseTime=25):
    ymd = creationTime.strftime('%Y%m%d')[2:]
    creationTime = creationTime.strftime('%Y-%m-%d %H:%M:%S')

    solName = createSolName(ymd, gateType, solNumber)

    # print(solName)

    filePath = folder + "\\" + solName + ".json"

    if circuitData is not None:
        solDict = {}
        circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
        for key in circuitDataKeys:
            solDict[key] = [item/(2*np.pi) for item in circuitData[key]]
    else:
        if circuitFile is None:
            circuitFile = folder + '\circuit.json'
        solDict = getFromjson(circuitFile)

    while (solNumber < 1000):
        if not os.path.isfile(filePath):
            newInfoDict = {
                "creationTime": creationTime,
                "gateType": gateType,
                "nOptimizationLvls": N,
                "riseTime": riseTime,
                "theta": x[0],
                "delta": x[1],
                "omegaPhi": x[2],
                "modulationTime": x[3],
                'gateFidelity': None,
                'times': None,
                'fidelitiesAtTimes': None,
                # Listor med avvikelser från korrekt Theta, delta, omegaPhi, opTime
                'deviations': [[], [], [], []],
                # Listor med fideliteten evaluerad vid dessa avvikelser
                'fidelities1D_Dev': [[], [], [], []],
                'fidelities2D_Theta_delta': None,
                'fidelities2D_Theta_omegaPhi': None,
                'fidelities2D_Theta_opTime': None,
                'fidelities2D_delta_omegaPhi': None,
                'fidelities2D_delta_opTime': None,
                'fidelities2D_omegaPhi_opTime': None
            }
            solDict.update(newInfoDict)

            dumpTojson(solDict, filePath)
            return
        else:
            try:
                existingSolDict = getFromjson(filePath)
                if (existingSolDict['creationTime'] == creationTime):
                    print("Can't add solution: Solution already exists!")
                    return
            except FileNotFoundError:
                pass

            solNumber += 1
            solName = createSolName(ymd, gateType, solNumber)
            filePath = folder + "\\" + solName + ".json"


def saveSolutionsTojson(results, gateType, N, folder, circuitFile=None, circuitData=None, dateAndTime=datetime.today()):
    for i in range(len(results)):
        x = results[i].x.tolist()
        addNewSolution(x, gateType, N, folder=folder, circuitFile=circuitFile, circuitData=circuitData, creationTime=dateAndTime)


def saveResToFile(result, algorithmName, iterations, runtime, algorithmDE=False, algorithmSHG=False, fileName="result.txt", dateAndTime=datetime.today()):
    resultFile = open(fileName, "a")
    
    dateStr = f'Result (timeout) from: {dateAndTime}\n'
    iterStr = f'Total iterations performed (or minima found): {iterations}\n'
    runtimeStr = f'Total runtime for algorithm: {runtime} seconds.\n'
    dividerStr = "##################################################\n"
    strList = [dividerStr, dateStr, iterStr, runtimeStr, dividerStr, "\n"]
    resultFile.writelines(strList)
    
    for res in result:
        x = res[0]
        fun = res[1]
        
        if algorithmDE:
            resultStr = f'The {algorithmName} algorithm gave a minimum at the point {x} with a convergence of {fun}.\n'
        elif algorithmSHG:
            resultStr = f'The {algorithmName} algorithm was examining the point {x}.\n'
        else: 
            resultStr = f'The {algorithmName} algorithm gave a minimum of {fun} at the point {x}.\n'
        
        resultFile.writelines([resultStr, "\n"])
    resultFile.write(dividerStr)
    resultFile.close()


def evaluateResult(x, fun, resultList, N=5):
    """
    This function is used to save the N best minima from an optimizer. The 
    function evaluates if the new minima fun is smaller than any of the 
    ones in resultList. If the length of resultList is shorter than N, the
    new minima will just be appended to resultList. 
    """
    if len(resultList[1]) < N:
        resultList[0].append(x)
        resultList[1].append(fun)
    else:
        if np.amax(resultList[1]) > fun:
            indexOfMax = np.argmax(resultList[1])
            resultList[0][indexOfMax] = x
            resultList[1][indexOfMax] = fun
    return resultList

######################################################################################################################################################################