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

dirPath = "Results/Qubit Pair 03/Solutions"
directory = os.fsencode(dirPath)

for file in os.listdir(directory):
    filename = os.fsdecode(file)

    if filename.endswith(".json"):
        filePath = dirPath + "/" + filename
        solDict = getFromjson(filePath)

        print(filename)
        print(solDict['gateFidelity'])


