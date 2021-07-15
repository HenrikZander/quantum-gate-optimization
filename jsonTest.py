from os import name
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
from variables import *
from functions import *
from plotting import *

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

