from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
from variables import *
from functions import *
from plotting import *
from dataManager import *

solsDict = getFromjson('solutions.json')
for solName in solsDict:
    testSolutionNew = getFromjson('circuit files/qubitPair01.json')

    for key in solsDict[solName]:
        if key != 'x':
            testSolutionNew[key] = solsDict[solName][key]
        else:
            testSolutionNew['riseTime'] = 25
            testSolutionNew['theta'] = solsDict[solName]['x'][0]
            testSolutionNew['delta'] = solsDict[solName]['x'][1]
            testSolutionNew['omegaPhi'] = solsDict[solName]['x'][2] / (2 * np.pi)
            testSolutionNew['modulationTime'] = solsDict[solName]['x'][3]
    
    dumpTojson(testSolutionNew, "solutions/" + solName[2:-5] + ".json")

'''
def getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex):
    labelsList = []
    eigenIndices = []
    for label in eigenEnergyDict:
        knownThetas = eigenEnergyDict[label][0]
        knownEigInds = eigenEnergyDict[label][3]
        for i, kTh in enumerate(knownThetas):
            if (kTh > theta):
                iAfter = i
                iBefore = i-1
                break
        if (knownEigInds[iBefore] == knownEigInds[iAfter]):
            eigenIndices.append(knownEigInds[iAfter])
        else:
            eigenIndices.append(-knownEigInds[iAfter])
        labelsList.append(label)
    usedLabels = ["" for _ in range(maxUsedIndex+1)]
    #print(eigenIndices)
    for i, ei in enumerate(eigenIndices):
        if abs(ei) > maxUsedIndex:
            pass
        elif ei >= 0:
            usedLabels[ei] = labelsList[i]
        else:
            usedLabels[abs(ei)] = "|???>"
    return usedLabels

maxProjIndex = 16
theta = -0.49
eigenEnergyDict = getFromjson('eigenenergies.json')

labelsUsed = getEigenstateLabels(eigenEnergyDict, theta, maxProjIndex)

print(labelsUsed)
'''

