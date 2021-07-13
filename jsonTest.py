from os import name
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime
from variables import *
from functions import *
from plotting import *

N=4
xCurr = [-3.66144215e-01, 7.76542110e-02, 5.31979461e+00, 1.35483199e+02]
fidelitiesAtTimes = [(0.999580975343943+0j), (0.9995711538296317+0j), (0.9995456227712491+0j), (0.9995303291269018+0j)]

#creationTime = datetime.today().strftime('%Y-%m-%d %H:%M:%S')
ymd = "210707"
creationTime = "2021-07-07 00:05:03"
gateType = "CZ"
usedNlvls = str(N) + "lvl"
solNumber = str(1)
solName = "x_" + ymd + "_" + gateType + "_" + solNumber + "_" + usedNlvls

print(solName)

solutionsDict = getjsonDict('solutions.json')
solutionsDict[solName] = {
    "creationTime": creationTime,
    "gateType": gateType,
    "nOptimizationLvls": N,
    'x': xCurr,
    'gateFidelity': fidelitiesAtTimes[0].real,
    'times': None,
    'fidelitiesAtTimes': None,
    'fidelities1D_Dev': [[], [], [], []], # Listor med fideliten evaluerad runt korrekt Theta, delta, omegaPhi, opTime
    'fidelities2D_Theta_delta': None,
    'fidelities2D_Theta_omegaPhi': None,
    'fidelities2D_Theta_opTime': None,
    'fidelities2D_delta_omegaPhi': None,
    'fidelities2D_delta_opTime': None,
    'fidelities2D_omegaPhi_opTime': None
}

dumpjsonDict(solutionsDict, 'solutions.json')

'''
solsDict = getjsonDict('solutions.json')
nameList = [xName for xName, _ in solsDict.items()]
print(nameList)
'''
'''
fidsTheta = [0.993, 0.998, 0.9996, 0.999, 0.995]

with open('solutions.json', 'r') as solsFile:
    solsDict = json.load(solsFile)
    for s in solsDict['solutions']:
        if(s['name'] == "x_2021-07-07_00:05:03_CZ_4lvl"):
            s['fidelities1D_Dev'][0] = fidsTheta
    
with open('solutions.json', 'w') as solsFile:
    json.dump(solsDict, solsFile, ensure_ascii=False, indent=4)
'''

'''
data = {}
data['people'] = []
data['people'].append({
    'name': 'Scott',
    'website': 'stackabuse.com',
    'from': 'Nebraska'
})
data['people'].append({
    'name': 'Larry',
    'website': 'google.com',
    'from': 'Michigan'
})
data['people'].append({
    'name': 'Tim',
    'website': 'apple.com',
    'from': 'Alabama'
})

with open('data.json', 'w', encoding='utf-8') as outfile:
    json.dump(data, outfile, ensure_ascii=False, indent=4)

with open('data.json') as json_file:
    data = json.load(json_file)
    for p in data['people']:
        print('Name: ' + p['name'])
        print('Website: ' + p['website'])
        print('From: ' + p['from'])
        print('')
'''


