from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
from variables import *
from functions import *
from plotting import *
from dataManager import *
import os

bounds = ([0., 1.],[-2.0, 3],[3, 4])
print(np.asarray(bounds)[:,1])

'''
dirPath = "Results/Qubit Pair 03/Solutions"
directory = os.fsencode(dirPath)

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename.endswith(".json"):
        filePath = dirPath + "/" + filename
        solDict = getFromjson(filePath)

        if solDict['gateType'] == 'CZ':
            print(filename)
            print(solDict['gateFidelity'])
'''
'''
dirPath = "Results/Qubit Pair 03"
directory = os.fsencode(dirPath)

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename == "eigenenergies.json":
        filePath = dirPath + "/" + filename
        eigDict = getFromjson(filePath)

        for key in eigDict:
            if (key != "circuit") and (key != "Thetas"):
                #del eigDict[key][0]
                eigDict[key][0], eigDict[key][1] = eigDict[key][1], eigDict[key][0]
        
        dumpTojson(eigDict, filePath)
'''

