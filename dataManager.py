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
from pathlib import Path
import math

######################################################################################################################################################################
# Helper function for the callback functions


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
# Functions that handle data loading and dumping to json-files.


def getFromjson(filePath):
    with open(filePath, 'r') as jsonFile:
        jsonDict = json.load(jsonFile)
    return jsonDict


def dumpTojson(data, filePath):
    with open(filePath, 'w') as jsonFile:
        json.dump(data, jsonFile, ensure_ascii=False, indent=4)


def createSolName(ymd, gateType, solNumber):
    return ymd + "_" + gateType + "_" + str(solNumber)


def addNewSolution(x, gateType, N, solNumber=1, creationTime=datetime.today(), folder='Demo Circuit', circuitFile=None, circuitData=None, riseTime=25, signalType=None):
    ymd = creationTime.strftime('%Y%m%d')[2:]
    creationTime = creationTime.strftime('%Y-%m-%d %H:%M:%S')

    solName = createSolName(ymd, gateType, solNumber)

    # print(solName)

    filePath = Path(folder, solName + ".json") #folder + "/" + solName + ".json"

    if circuitData is not None:
        solutionDict = {}
        circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
        for key in circuitDataKeys:
            solutionDict[key] = [item/(2*np.pi) for item in circuitData[key]]
    else:
        if circuitFile is None:
            circuitFile = folder + '/circuit.json'
        solutionDict = getFromjson(circuitFile)

    if signalType == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
    elif signalType == 'cos':
        x0name = 'theta'
        x1name = 'delta'
    
    while (solNumber < 1000):
        if not os.path.isfile(filePath):
            newInfoDict = {
                "creationTime": creationTime,
                "comboSignal": False,
                "signalType": signalType,
                "gateType": gateType,
                "nOptimizationLvls": N,
                "riseTime": riseTime,
                x0name: x[0],
                x1name: x[1],
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
            
            solutionDict.update(newInfoDict)

            dumpTojson(solutionDict, filePath)
            return
        else:
            try:
                existingSolutionDict = getFromjson(filePath)
                if (existingSolutionDict['creationTime'] == creationTime):
                    print("Can't add solution: Solution already exists!")
                    return
            except FileNotFoundError:
                pass

            solNumber += 1
            solName = createSolName(ymd, gateType, solNumber)
            filePath = Path(folder, solName + ".json") # folder + "/" + solName + ".json"


def saveSolutionsTojson(results, gateType, N, folder, circuitFile=None, circuitData=None, dateAndTime=datetime.today(), signalType=None):
    for i in range(len(results)):
        x = results[i].x.tolist()
        addNewSolution(x, gateType, N, folder=folder, circuitFile=circuitFile, circuitData=circuitData, creationTime=dateAndTime, signalType=signalType)


def generateComboSolutionFile(solPath1, solPath2, comboGateType, comboSolPath=None):
    solDict1 = getFromjson(solPath1)
    solDict2 = getFromjson(solPath2)

    if solDict1['signalType'] != solDict2['signalType']:
        print("The two provided solutions must have the same signal type!")
        return
    elif not math.isclose(solDict1['riseTime'], solDict2['riseTime']):
        print("Currently, superposing two signals with different rise times isn't supported.")
        return

    circuitData1 = getCircuitData(solDict1, convertToGradPerSecond=False)
    circuitData2 = getCircuitData(solDict2, convertToGradPerSecond=False)

    if not areCircuitDatasIdentical(circuitData1, circuitData2):
        print("The two provided solutions must use the same circuit parameters!")
        return

    if solDict1['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
    elif solDict1['signalType'] == 'cos':
        x0name = 'theta'
        x1name = 'delta'

    if not math.isclose(solDict1[x0name], solDict2[x0name]):
        print("The two superposed signals need to have the same DC value (Theta or A)")
        return

    creationTime = datetime.today()

    comboSolDict = {
        "creationTime": creationTime.strftime('%Y-%m-%d %H:%M:%S'),
        "comboSignal": True,
        "signalType": solDict1['signalType'],
        "gateType": comboGateType,
        "nOptimizationLvls": max(solDict1['nOptimizationLvls'], solDict2['nOptimizationLvls']),
        "riseTime": solDict1['riseTime'],
        x0name: solDict1[x0name],
        x1name+"s": [solDict1[x1name], solDict2[x1name]],
        "omegaPhis": [solDict1['omegaPhi'], solDict2['omegaPhi']],
        "modulationTimes": [solDict1['modulationTime'], solDict2['modulationTime']],
        'gateFidelity': None,
        'times': None,
        'fidelitiesAtTimes': None
    }
    circuitData1.update(comboSolDict)
    
    if comboSolPath is None:
        folder = str(solPath1.parent)

        solNumber = 1
        while solNumber < 1000:
            comboSolName = creationTime.strftime('%Y%m%d')[2:] + '_' + 'combo' + '_' + comboGateType + '_' + str(solNumber)
            comboSolPath = Path(folder, comboSolName + ".json")

            if not os.path.isfile(comboSolPath):
                break

            solNumber += 1
    dumpTojson(circuitData1, comboSolPath)


def areCircuitDatasIdentical(circuitData1, circuitData2):
    for key in circuitData2:
            for index, value in enumerate(circuitData2[key]):
                if not math.isclose(value, circuitData1[key][index]):
                    return False
    return True

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

def getCircuitData(solutionDict, convertToGradPerSecond=True):
    if convertToGradPerSecond:
        koeff = 2*np.pi
    else:
        koeff = 1
    
    circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
    circuitData = {}
    for key in circuitDataKeys:
        circuitData[key] = [koeff*item for item in solutionDict[key]]
    return circuitData


######################################################################################################################################################################