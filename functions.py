######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: functions.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 27 February 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import os
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time
import random
import scipy
from numba import njit

from plotting import *
from model import *

import json
from datetime import datetime
import math
from dataManager import *
import copy

######################################################################################################################################################################
# Global variables

i = 0
global maxRuntime
maxRuntime = 40000

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
# json helper functions

'''
def getSolutionNameList():
    solsDict = getFromjson('solutions.json')
    return [xName for xName in solsDict]
'''

'''
def createSolNameOld(ymd, gateType, solNumberStr, usedNlvls):
    return "x_" + ymd + "_" + gateType + "_" + solNumberStr + "_" + usedNlvls
'''

def createSolName(ymd, gateType, solNumberStr):
    return ymd + "_" + gateType + "_" + solNumberStr

'''
def addNewSolutionOld(x, gateType, N, solNumber=1, creationTime=datetime.today(), fileName='solutions.json'):
    ymd = creationTime.strftime('%Y%m%d')[2:]
    creationTime = creationTime.strftime('%Y-%m-%d %H:%M:%S')
    usedNlvls = str(N) + "lvl"
    solName = createSolName(ymd, gateType, str(solNumber), usedNlvls)

    # print(solName)

    solsDict = getFromjson(fileName)
    while (solNumber < 1000):
        if solName not in solsDict:
            solsDict[solName] = {
                "creationTime": creationTime,
                "gateType": gateType,
                "nOptimizationLvls": N,
                'x': x,
                'gateFidelity': None,
                'times': None,
                'fidelitiesAtTimes': None,
                # Listor med avvikelser från korrekt Theta, delta, omegaPhi, opTime
                'deviations': [[], [], [], []],
                # Listor med fideliten evaluerad vid dessa avvikelser
                'fidelities1D_Dev': [[], [], [], []],
                'fidelities2D_Theta_delta': None,
                'fidelities2D_Theta_omegaPhi': None,
                'fidelities2D_Theta_opTime': None,
                'fidelities2D_delta_omegaPhi': None,
                'fidelities2D_delta_opTime': None,
                'fidelities2D_omegaPhi_opTime': None
            }

            dumpTojson(solsDict, fileName)
            return None
        else:
            if (solsDict[solName]['creationTime'] == creationTime):
                print("Can't add solution: Solution already exists in solutions.json")
                return None
            solNumber += 1
            solName = createSolName(ymd, gateType, str(solNumber), usedNlvls)
'''

def addNewSolution(x, gateType, N, solNumber=1, creationTime=datetime.today(), folder='solutions_qubitPair01', circuitFile='circuit files/qubitPair01.json', circuitData=None, riseTime=25):
    ymd = creationTime.strftime('%Y%m%d')[2:]
    creationTime = creationTime.strftime('%Y-%m-%d %H:%M:%S')

    solName = createSolName(ymd, gateType, str(solNumber))

    # print(solName)

    filePath = "./" + folder + "/" + solName + ".json"

    if circuitData is None:
        solDict = getFromjson(circuitFile)
    else:
        solDict = copy.deepcopy(circuitData)

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
                # Listor med fideliten evaluerad vid dessa avvikelser
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
                    print("Can't add solution: Solution already exists in solutions.json")
                    return
            except FileNotFoundError:
                pass

            solNumber += 1
            solName = createSolName(ymd, gateType, str(solNumber))
            filePath = "./" + folder + "/" + solName + ".json"


def saveSolutionsTojson(results, gateType, N, folder, circuitFile=None, circuitData=None, dateAndTime=datetime.today()):
    for i in range(len(results)):
        x = results[i].x.tolist()
        addNewSolution(x, gateType, N, folder=folder, circuitFile=circuitFile, circuitData=circuitData, creationTime=dateAndTime)


######################################################################################################################################################################
# Simulation functions


def getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex):
    theta = -abs(theta)
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


'''
def simulateHamiltonian(x=None, xName=None, sinStepHamiltonian=True, rotatingFrame=False, initialStateIndex=1, highestProjectionIndex=8, N=4, circuitData=None):
    """
    This function simulates the population transfers between 
    different eigenstates of the 4-level hamiltonian from 
    McKay eq. 11. It also calculates the gate fidelity for both
    the iSWAP and CZ gate using the 4-level system.
    ---------------------------------------------------------
    INPUT:
            x (array(float)) {Optional}: The parameter set that the simulation is be performed for. If x is not provided, xName needs to be.
            sinStepHamiltonian (boolean) {Optional}: If True the amplitude of AC part of the flux signal that is applied to the tunable bus will be sinusodially modulated.
            rotatingFrame (boolean) {Optional}: If True the states will be transformed into the rotating frame after the time evolution has been completed.
            initialStateIndex (int) {Optional}: This parameter decides which eigenstate that will be the initial state in the time evolution. If initialStateIndex=0 the eigenstate with the lowest associated energy will be the inital state.
            highestProjectionIndex (int) {Optional}: The eigenstates between, and including, the one with the lowest energy up to the (highestProjectionIndex)-lowest eigenstate will be projected onto.
            N (int) {Optional}: How many energy levels that should be accounted for in the simulations.
            xName (string) {Optional}: The name of (i.e. key to) the solution in solutions.json that one wants to simulate. If xName is not provided, x needs to be.
    ---------------------------------------------------------
    OUTPUT:
            gateFidelity_iSWAP, gateFidelity_CZ (float, float): The gate fidelity for both the iSWAP and CZ gate.
    ---------------------------------------------------------
    """

    if (xName is None):
        if (x is None):
            print("Must supply either the name of an x stored in a json file or a list x manually (or both)")
            return None
    else:
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getFromjson('solutions.json')
        x = solsDict[xName]['x']

    # Calculate the dimension of the tensor states and set the simulation time.
    D = N**3
    simulationTime = int(x[-1]) + 10

    # Calculate the eigenstates and eigenenergies of the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x, N=N, getBBHamiltonianComps=True)

    # Calculate the tunable bus frequency when only the DC part of the flux is active.
    omegaTBDC = coeffomegaTB(omegas[2], x[0])

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux only has it's DC part.
    eigenStatesAndEnergies = getThetaEigenstates(x, hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)

    # Get eigenindices.
    eigIndices = getIndices(N, eigenStatesAndEnergies[1])

    # Calculate the unitary for transforming the hamiltonian to the eigen basis.
    eigenBasisUnitary = getEBUnitary(x, eigenStatesAndEnergies, N, D)

    # Get the hamiltonian that has a sinusodially modulated AC flux and also is in the eigen basis.
    hamiltonian = getHamiltonian(x, N=N, eigEs=eigenStatesAndEnergies[0], U_e=eigenBasisUnitary, sinStepHamiltonian=sinStepHamiltonian)

    # Change the simulation settings and create the timestamps for where the evolution is to be evaluated.
    options = solver.Options()
    options.nsteps = 10000
    timeStamps = np.linspace(0, simulationTime, simulationTime*5)

    # Create the initial state.
    initialState = Qobj(basis(D, initialStateIndex), dims=[[N, N, N], [1, 1, 1]])

    # Time evolve the initial state.
    result = sesolve(hamiltonian, initialState, timeStamps, e_ops=[], options=options, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': omegas[2], 'operationTime': x[3], 'omegaTBTh': omegaTBDC})
    states = result.states

    # Transform into the rotating frame.
    Hrot = np.array(hamiltonian[0])
    Hrot[eigIndices[-1], eigIndices[-1]] = Hrot[eigIndices[-2], eigIndices[-2]] + Hrot[eigIndices[-3], eigIndices[-3]] - Hrot[eigIndices[0], eigIndices[0]]
    Hrot = Qobj(Hrot, dims=[[N, N, N], [N, N, N]])

    if rotatingFrame:
        for i, t in enumerate(timeStamps):
            unitaryRotatingFrame = getRFUnitary(Hrot, t)
            states[i] = unitaryRotatingFrame*states[i]

    # Calculate the expectation values a projection operator.
    projectionOperators = []
    #pIndexes = [0,1,2,5,highestProjectionIndex]

    for index in range(highestProjectionIndex+1):
        operator = Qobj(basis(D, index), dims=[[N, N, N], [1, 1, 1]]) * Qobj(basis(D, index), dims=[[N, N, N], [1, 1, 1]]).dag()
        projectionOperators.append(operator)

    expectationValues = expect(projectionOperators, states)

    # Calculate gate fidelity for both iSWAP and CZ.
    gateFidelity_iSWAP, timesiSWAP = getGateFidelity(x, N=N, iSWAP=True, circuitData=circuitData)
    gateFidelity_CZ, timesCZ = getGateFidelity(x, N=N, CZ=True, circuitData=circuitData)

    # Print fidelity
    print(f'################################################\n\nGate fidelity for iSWAP: {gateFidelity_iSWAP}.\n\nGate fidelity for CZ: {gateFidelity_CZ}.\n\n################################################')

    # Plot the expectation values.
    plt.figure(figsize=(8, 7))
    #labels = ["|000>", "|010>", "|100>", "|001>", "|020>", "|110>", "|011>", "|200>", "|101>"]
    # labels = []
    #labels = ['|000>','|010>','|100>','|001>','|020>','|110>','|200>','|011>','|101>','|002>','|003>','|120>']
    eigenEnergyDict = getFromjson('eigenenergies.json')
    theta = x[0]
    maxUsedIndex = highestProjectionIndex
    labels = getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex)

    linestyle = ['-', '--', '-.', ':']
    for index, values in enumerate(expectationValues):
        plt.plot(timeStamps, values, ls=linestyle[math.floor(index / 10) % 4])
        # eigenOrder = (findEigenIndex(x, eigenStateIndex=index))[1]
        # labels.append(f'|{eigenOrder[0]}{eigenOrder[1]}{eigenOrder[2]}>')

    plt.grid()
    plt.ylim([0, 1.1])
    plt.xlim([0, timeStamps[-1]])
    leg = plt.legend(labels, fontsize=19, loc='center right')
    plt.xlabel("Tid efter grindstart [ns]", fontsize=26)
    plt.ylabel("Population", fontsize=26)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(5.0)
    plt.show()

    return gateFidelity_iSWAP, gateFidelity_CZ
'''

def getCircuitData(solutionDict):
    circuitDataKeys = ['frequencies', 'anharmonicities', 'couplings']
    circuitData = {}
    for key in circuitDataKeys:
        circuitData[key] = [2*np.pi*item for item in solutionDict[key]]
    return circuitData


def simulatePopTransfer(solutionPath, sinStepHamiltonian=True, rotatingFrame=True, initialStateIndex=1, highestProjectionIndex=12, N=4):
    """
    This function simulates the population transfers between 
    different eigenstates of the 4-level hamiltonian from 
    McKay eq. 11.
    ---------------------------------------------------------
    INPUT:
            solutionPath (string): The path to the solution file that the simulation is be performed for.
            sinStepHamiltonian (boolean) {Optional}: If True the amplitude of AC part of the flux signal that is applied to the tunable bus will be sinusodially modulated.
            rotatingFrame (boolean) {Optional}: If True the states will be transformed into the rotating frame after the time evolution has been completed.
            initialStateIndex (int) {Optional}: This parameter determines which eigenstate will be the initial state in the time evolution. If initialStateIndex=0,
                the eigenstate with the lowest associated energy will be the initial state.
            highestProjectionIndex (int) {Optional}: Determines how many eigenstates to project the evolving state onto. If highestProjectionIndex = n, the population of
                the n+1 lowest eigenstates are shown (since the lowest possible eigenIndex is zero).
            N (int) {Optional}: How many energy levels that should be accounted for in the simulations.
    ---------------------------------------------------------
    """

    # Load a solution from a file and create a list x defining the control signal, converting omegaPhi from GHz to Grad/s.
    solutionDict = getFromjson(fileName=solutionPath)
    x = [solutionDict['theta'], solutionDict['delta'], 2*np.pi*solutionDict['omegaPhi'], solutionDict['modulationTime']]

    # Load circuit data into a dictionary and change units from GHz to Grad/s.
    circuitData = getCircuitData(solutionDict)

    # Calculate the dimension of the tensor states and set the simulation time.
    D = N**3
    simulationTime = int(x[-1]) + 10

    # Calculate the eigenstates and eigenenergies of the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData)

    # Calculate the tunable bus frequency when only the DC part of the flux is active.
    omegaTBDC = coeffomegaTB(circuitData['frequencies'][2], x[0])

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux only has its DC part.
    eigenStatesAndEnergies = getThetaEigenstates(x, hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)

    # Get eigenindices.
    eigIndices = getIndices(N, eigenStatesAndEnergies[1])

    # Calculate the unitary for transforming the hamiltonian to the eigen basis.
    eigenBasisUnitary = getEBUnitary(x, eigenStatesAndEnergies, N, D)

    # Get the hamiltonian that has a sinusodially modulated AC flux and also is in the eigen basis.
    hamiltonian = getHamiltonian(x, N=N, eigEs=eigenStatesAndEnergies[0], U_e=eigenBasisUnitary, sinStepHamiltonian=sinStepHamiltonian, circuitData=circuitData)

    # Change the simulation settings and create the timestamps for where the evolution is to be evaluated.
    options = solver.Options()
    options.nsteps = 10000
    timeStamps = np.linspace(0, simulationTime, simulationTime*5)

    # Create the initial state.
    initialState = Qobj(basis(D, initialStateIndex), dims=[[N, N, N], [1, 1, 1]])

    # Time evolve the initial state.
    result = sesolve(hamiltonian, initialState, timeStamps, [], options=options, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': circuitData['frequencies'][2], 'operationTime': x[3], 'omegaTBTh': omegaTBDC, 'riseTime': solutionDict['riseTime']})
    states = result.states

    # Transform into the rotating frame.
    Hrot = np.array(hamiltonian[0])
    Hrot[eigIndices[-1], eigIndices[-1]] = Hrot[eigIndices[-2], eigIndices[-2]] + Hrot[eigIndices[-3], eigIndices[-3]] - Hrot[eigIndices[0], eigIndices[0]]
    Hrot = Qobj(Hrot, dims=[[N, N, N], [N, N, N]])

    if rotatingFrame:
        for i, t in enumerate(timeStamps):
            unitaryRotatingFrame = getRFUnitary(Hrot, t)
            states[i] = unitaryRotatingFrame*states[i]

    # Calculate the expectation values a projection operator.
    projectionOperators = []
    #pIndexes = [0,1,2,5,highestProjectionIndex]

    for index in range(highestProjectionIndex+1):
        operator = Qobj(basis(D, index), dims=[[N, N, N], [1, 1, 1]]) * Qobj(basis(D, index), dims=[[N, N, N], [1, 1, 1]]).dag()
        projectionOperators.append(operator)

    expectationValues = expect(projectionOperators, states)

    # Plot the expectation values.
    plt.figure(figsize=(8, 7))
    #labels = ["|000>", "|010>", "|100>", "|001>", "|020>", "|110>", "|011>", "|200>", "|101>"]
    # labels = []
    #labels = ['|000>','|010>','|100>','|001>','|020>','|110>','|200>','|011>','|101>','|002>','|003>','|120>']
    eigenEnergyDict = getFromjson('eigenenergies.json')
    theta = x[0]
    maxUsedIndex = highestProjectionIndex
    labels = getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex)

    linestyle = ['-', '--', '-.', ':']
    for index, values in enumerate(expectationValues):
        plt.plot(timeStamps, values, ls=linestyle[math.floor(index / 10) % 4])

    plt.grid()
    plt.ylim([0, 1.1])
    plt.xlim([0, timeStamps[-1]])
    leg = plt.legend(labels, fontsize=19, loc='center right')
    plt.xlabel("Tid efter grindstart [ns]", fontsize=26)
    plt.ylabel("Population", fontsize=26)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.tight_layout()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(5.0)
    plt.show()

'''
def plotFidelityOld(x=None, iSWAP=False, CZ=False, xName=None, useSavedPlot=False, saveTojson=False, circuitData=None):
    if (xName is not None):
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getFromjson('solutions.json')
        x = solsDict[xName]['x']
        x[2] = x[2] / (2*np.pi)
        if useSavedPlot:
            F = solsDict[xName]['fidelitiesAtTimes']
            times = solsDict[xName]['times']
        else:
            indices = np.linspace(-116, -1, 116).astype(int)
            F, times = getGateFidelity(x, N=4, iSWAP=iSWAP, CZ=CZ, tIndices=indices, circuitData=circuitData)
    else:
        if (x is None):
            print("Must supply either the name of an x stored in a json file or a list x manually (or both)")
            return None
        else:
            indices = np.linspace(-116, -1, 116).astype(int)
            F, times = getGateFidelity(x, N=4, iSWAP=iSWAP, CZ=CZ, tIndices=indices, circuitData=circuitData)

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

    if (saveTojson == True):
        if (xName is None):
            print("Warning: Couldn't save results. To save results, you need to specify a solution listed in solutions.json")
        else:
            solsDict[xName]['times'] = times
            solsDict[xName]['fidelitiesAtTimes'] = F
            solsDict[xName]['gateFidelity'] = F[-76]
            dumpTojson(solsDict,'solutions.json')
'''

def plotFidelity(solutionPath, useSavedPlot=False, saveToFile=False):
    solutionDict = getFromjson(fileName=solutionPath)
    x = (solutionDict['theta'], solutionDict['delta'], solutionDict['omegaPhi'], solutionDict['modulationTime'])
    
    circuitData = getCircuitData(solutionDict)

    if solutionDict['gateType'] == 'iSWAP':
        iSWAP = True
    else:
        iSWAP = False
    if solutionDict['gateType'] == 'SWAP':
        SWAP = True
    else:
        SWAP = False
    if solutionDict['gateType'] == 'CZ':
        CZ = True
    else:
        CZ = False
    
    if useSavedPlot:
        F = solutionDict['fidelitiesAtTimes']
        times = solutionDict['times']
    else:
        indices = np.linspace(-116, -1, 116).astype(int)
        F, times = getGateFidelity(x, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=indices, circuitData=circuitData)

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


def deltaPulsePlot():
    operationTime = 75
    x = np.linspace(0, operationTime, 500)
    y = []

    for time in x:
        y.append(sinBox(time, operationTime))

    plt.plot(x, y)
    plt.plot([25, 25], [0, 1.1], 'r--', [operationTime-25, operationTime-25], [0, 1.1], 'r--')
    plt.xlabel('Tid [ns]', fontsize=16)
    plt.ylabel('Del av full amplitud, $\delta(t) / \delta_0$', fontsize=16)
    plt.xlim([0, operationTime])
    plt.ylim([0, 1.1])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.annotate('', xy=(0, 0.5), xytext=(25, 0.5), arrowprops=dict(arrowstyle='<->'))
    plt.annotate("$t_{Stig}$", xy=(6, 0.53), fontsize=18)
    plt.annotate('', xy=(operationTime-25, 0.5), xytext=(operationTime, 0.5), arrowprops=dict(arrowstyle='<->'))
    plt.annotate("$t_{Fall}$", xy=(operationTime-10, 0.53), fontsize=18)
    plt.show()

'''
def getRobustnessPlotOld(x=None, xName=None, iSWAP=False, CZ=False, checkTheta=False, checkDelta=False, checkOmegaPhi=False, checkOpTime=False, nPointsList=[9], maxDevs=[5e-4, 1e-3, 2*np.pi * 2e-3, 4e0], useSavedPlot=False, saveTojson=False, circuitData=None):
    if (xName is None):
        if (x is None):
            print("Must supply either the name of an x stored in a json file or a list x manually (or both)")
            return None
    else:
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getFromjson('solutions.json')
        x = solsDict[xName]['x']
    
    if (useSavedPlot and saveTojson):
        print("Re-saving something in a json file that's already saved there is rather pointless, don't ya think?")
        saveTojson = False
    
    checkList = [checkTheta, checkDelta, checkOmegaPhi, checkOpTime]
    if (sum(checkList) < 1):
        print("You need to specify at least one parameter to check")
    elif (sum(checkList) > 2):
        print("Too many parameters! There is currently only support for checking one or two parameters simultaneously")
    elif (len(nPointsList) > 2):
        print("Can't plot in more than two dimensions!")
    elif (len(nPointsList) < 1):
        print("Can't plot in less than one dimension!")
    else:
        xDev = [xi for xi in x]
        plt.figure(figsize=(11, 7))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        axlabels = ["Avvikelse från hittat $\Theta$ [$\Phi_0$]", "Avvikelse från hittat $\delta$ [$\Phi_0$]", "Avvikelse från hittat $\omega_\Phi$ [MHz]", "Avvikelse från hittat $t_{MOD}$ [ns]"]
        legendStrs = ["$\Theta = %.4f$" %x[0], "$\delta = %.4f$" %x[1], "$\omega_\Phi = %.1f$" %(x[2]/(2*np.pi)*1000), "$t_{MOD} = %.1f$" %x[3]]
        # maxDevs = [3e-4, 1e-3, 2*np.pi * 8e-3, 40e0] # Testing extended bounds, looking for chevrony stuff

        xIndices = []
        if (checkTheta):
            xIndices.append(0)
        if (checkDelta):
            xIndices.append(1)
        if (checkOmegaPhi):
            xIndices.append(2)
        if (checkOpTime):
            xIndices.append(3)

        if (sum(checkList) == 1):
            xIndex = xIndices[0]
            if (useSavedPlot):
                deviations = np.array(solsDict[xName]['deviations'][xIndex])
                fidelities = solsDict[xName]['fidelities1D_Dev'][xIndex]
            else:
                deviations = np.linspace(-maxDevs[xIndex], maxDevs[xIndex], nPointsList[0])
                fidelities = []
                for i, d in enumerate(deviations):
                    xDev[xIndex] = x[xIndex] + d
                    fidelity, _ = getGateFidelity(xDev, N=4, iSWAP=iSWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
                    fidelities.append(fidelity[0])
                    statusBar((i+1)*100/nPointsList[0])

            if (xIndex == 2):
                # Re-scale x-axis to MHz
                nTicks = 9
                locs = np.linspace(deviations[0], deviations[-1], nTicks)
                newDevMax = deviations[-1]/(2*np.pi)*1e3
                newTicks = np.linspace(-newDevMax,newDevMax,nTicks)
                plt.xticks(locs,newTicks)
            
            plt.plot(deviations, fidelities, 'b-')
            plt.plot([0, 0], [0, 1], 'r--')
            plt.legend(["Fidelitet", legendStrs[xIndex]], fontsize=19, loc="upper right")
            plt.grid()
            plt.xlabel(axlabels[xIndex], fontsize=26)
            plt.ylabel("Fidelitet", fontsize=26)
            plt.xlim([deviations[0], deviations[-1]])
            plt.ylim([0.99, 1])
            plt.tight_layout()
            plt.show()

            if (saveTojson):
                solsDict[xName]['deviations'][xIndex] = deviations.tolist()
                solsDict[xName]['fidelities1D_Dev'][xIndex] = fidelities
                dumpTojson(solsDict, 'solutions.json')
        elif (sum(checkList) == 2):
            if (len(nPointsList) == 1):
                nPointsList.append(nPointsList[0])
            
            plt.xlabel(axlabels[xIndices[0]], fontsize=26)
            plt.ylabel(axlabels[xIndices[1]], fontsize=26)
            iLegendStr = legendStrs[xIndices[0]]
            jLegendStr = legendStrs[xIndices[1]]

            if (xIndices == [0,1]):
                listName = "fidelities2D_Theta_delta"
            elif (xIndices == [0,2]):
                listName = "fidelities2D_Theta_omegaPhi"
            elif (xIndices == [0,3]):
                listName = "fidelities2D_Theta_opTime"
            elif (xIndices == [1,2]):
                listName = "fidelities2D_delta_omegaPhi"
            elif (xIndices == [1,3]):
                listName = "fidelities2D_delta_opTime"
            elif (xIndices == [2,3]):
                listName = "fidelities2D_omegaPhi_opTime"
            else:
                print("Error: If you see this, something went wrong")
                return None

            if (useSavedPlot):
                fidelities2D = np.array(solsDict[xName][listName][0])
                iDeviations = np.array(solsDict[xName][listName][1])
                jDeviations = np.array(solsDict[xName][listName][2])
            else:
                iDeviations = np.linspace(-maxDevs[xIndices[0]], maxDevs[xIndices[0]], nPointsList[0])
                jDeviations = np.linspace(-maxDevs[xIndices[1]], maxDevs[xIndices[1]], nPointsList[1])
                
                fidelities = []
                for j, jDev in enumerate(jDeviations):
                    xDev[xIndices[1]] = x[xIndices[1]] + jDev
                    for i, iDev in enumerate(iDeviations):
                        xDev[xIndices[0]] = x[xIndices[0]] + iDev
                        fidelity, _ = getGateFidelity(xDev, N=4, iSWAP=iSWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
                        fidelities.append(fidelity[0])
                        statusBar((j*nPointsList[0] + (i+1))*100/(nPointsList[0]*nPointsList[1]))
                fidelities2D = np.array(fidelities).reshape(nPointsList[1], nPointsList[0])

            nyTicks = nxTicks = 9

            xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
            xticks = np.linspace(iDeviations[0], iDeviations[-1], nxTicks)
            plt.xticks(xlocs,xticks)
            ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
            yticks = np.linspace(jDeviations[0], jDeviations[-1], nyTicks)
            plt.yticks(ylocs,yticks)

            if (xIndices[0] == 2):
                # Re-scale x-axis to MHz
                xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
                newDevMax = iDeviations[-1]/(2*np.pi)*1e3
                xticks = np.linspace(-newDevMax, newDevMax, nxTicks)
                plt.xticks(xlocs, xticks)
            elif (xIndices[1] == 2):
                # Re-scale y-axis to MHz
                ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
                newDevMax = jDeviations[-1]/(2*np.pi)*1e3
                yticks = np.linspace(-newDevMax, newDevMax, nyTicks)
                plt.yticks(ylocs, yticks)

            plt.imshow(fidelities2D, interpolation="bilinear", origin="lower")
            plt.colorbar(label="Fidelitet", orientation="vertical")

            plt.axvline(x=(nPointsList[0]-1)/2, color='k')
            plt.axhline(y=(nPointsList[1]-1)/2, color='r')
            plt.legend([iLegendStr, jLegendStr], fontsize=19, loc="upper right")
            plt.show()

            if (saveTojson):
                solsDict[xName][listName] = [[],[],[]]
                solsDict[xName][listName][0] = fidelities2D.tolist()
                solsDict[xName][listName][1] = iDeviations.tolist()
                solsDict[xName][listName][2] = jDeviations.tolist()
                dumpTojson(solsDict, 'solutions.json')
'''

def getRobustnessPlot(solutionPath, checkTheta=False, checkDelta=False, checkOmegaPhi=False, checkOpTime=False, nPointsList=[9], maxDevs=[5e-4, 1e-3, 2e-3, 4e0], useSavedPlot=False, saveToFile=False):
    solDict = getFromjson(solutionPath)
    x = [solDict['theta'], solDict['delta'], solDict['omegaPhi'], solDict['modulationTime']]
    circuitData = getCircuitData(solDict)

    if solDict['gateType'] == 'iSWAP':
        iSWAP = True
    else:
        iSWAP = False
    if solDict['gateType'] == 'SWAP':
        SWAP = True
    else:
        SWAP = False
    if solDict['gateType'] == 'CZ':
        CZ = True
    else:
        CZ = False

    if (useSavedPlot and saveToFile):
        print("Re-saving data in a file where that data is already saved is rather pointless, don't ya think?")
        saveToFile = False

    checkList = [checkTheta, checkDelta, checkOmegaPhi, checkOpTime]
    if (sum(checkList) < 1):
        print("You need to specify at least one parameter to check")
    elif (sum(checkList) > 2):
        print("Too many parameters! There is currently only support for checking one or two parameters simultaneously")
    elif (len(nPointsList) > 2):
        print("Can't plot in more than two dimensions!")
    elif (len(nPointsList) < 1):
        print("Can't plot in less than one dimension!")
    else:
        xDev = [xi for xi in x]
        plt.figure(figsize=(11, 7))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        axlabels = ["Avvikelse från hittat $\Theta$ [$\Phi_0$]", "Avvikelse från hittat $\delta$ [$\Phi_0$]", "Avvikelse från hittat $\omega_\Phi$ [MHz]", "Avvikelse från hittat $t_{MOD}$ [ns]"]
        legendStrs = ["$\Theta = %.4f$" %x[0], "$\delta = %.4f$" %x[1], "$\omega_\Phi = %.1f$" %(x[2]*1000), "$t_{MOD} = %.1f$" %x[3]]
        # maxDevs = [3e-4, 1e-3, 2*np.pi * 8e-3, 40e0] # Testing extended bounds, looking for chevrony stuff

        xIndices = []
        if (checkTheta):
            xIndices.append(0)
        if (checkDelta):
            xIndices.append(1)
        if (checkOmegaPhi):
            xIndices.append(2)
        if (checkOpTime):
            xIndices.append(3)

        if (sum(checkList) == 1):
            xIndex = xIndices[0]
            if (useSavedPlot):
                deviations = np.array(solDict['deviations'][xIndex])
                fidelities = solDict['fidelities1D_Dev'][xIndex]
            else:
                deviations = np.linspace(-maxDevs[xIndex], maxDevs[xIndex], nPointsList[0])
                fidelities = []
                for i, d in enumerate(deviations):
                    xDev[xIndex] = x[xIndex] + d
                    fidelity, _ = getGateFidelity(xDev, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
                    fidelities.append(fidelity[0])
                    statusBar((i+1)*100/nPointsList[0])

            if (xIndex == 2):
                # Re-scale x-axis to MHz
                nTicks = 9
                locs = np.linspace(deviations[0], deviations[-1], nTicks)
                newDevMax = deviations[-1]*1e3
                newTicks = np.linspace(-newDevMax,newDevMax,nTicks)
                plt.xticks(locs,newTicks)
            
            plt.plot(deviations, fidelities, 'b-')
            plt.plot([0, 0], [0, 1], 'r--')
            plt.legend(["Fidelitet", legendStrs[xIndex]], fontsize=19, loc="upper right")
            plt.grid()
            plt.xlabel(axlabels[xIndex], fontsize=26)
            plt.ylabel("Fidelitet", fontsize=26)
            plt.xlim([deviations[0], deviations[-1]])
            plt.ylim([0.99, 1])
            plt.tight_layout()
            plt.show()

            if (saveToFile):
                solDict['deviations'][xIndex] = deviations.tolist()
                solDict['fidelities1D_Dev'][xIndex] = fidelities
                dumpTojson(solDict, solutionPath)
        elif (sum(checkList) == 2):
            if (len(nPointsList) == 1):
                nPointsList.append(nPointsList[0])
            
            plt.xlabel(axlabels[xIndices[0]], fontsize=26)
            plt.ylabel(axlabels[xIndices[1]], fontsize=26)
            iLegendStr = legendStrs[xIndices[0]]
            jLegendStr = legendStrs[xIndices[1]]

            if (xIndices == [0,1]):
                listName = "fidelities2D_Theta_delta"
            elif (xIndices == [0,2]):
                listName = "fidelities2D_Theta_omegaPhi"
            elif (xIndices == [0,3]):
                listName = "fidelities2D_Theta_opTime"
            elif (xIndices == [1,2]):
                listName = "fidelities2D_delta_omegaPhi"
            elif (xIndices == [1,3]):
                listName = "fidelities2D_delta_opTime"
            elif (xIndices == [2,3]):
                listName = "fidelities2D_omegaPhi_opTime"
            else:
                print("Error: If you see this, something went wrong")
                return

            if (useSavedPlot):
                fidelities2D = np.array(solDict[listName][0])
                iDeviations = np.array(solDict[listName][1])
                jDeviations = np.array(solDict[listName][2])
            else:
                iDeviations = np.linspace(-maxDevs[xIndices[0]], maxDevs[xIndices[0]], nPointsList[0])
                jDeviations = np.linspace(-maxDevs[xIndices[1]], maxDevs[xIndices[1]], nPointsList[1])
                
                fidelities = []
                for j, jDev in enumerate(jDeviations):
                    xDev[xIndices[1]] = x[xIndices[1]] + jDev
                    for i, iDev in enumerate(iDeviations):
                        xDev[xIndices[0]] = x[xIndices[0]] + iDev
                        fidelity, _ = getGateFidelity(xDev, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
                        fidelities.append(fidelity[0])
                        statusBar((j*nPointsList[0] + (i+1))*100/(nPointsList[0]*nPointsList[1]))
                fidelities2D = np.array(fidelities).reshape(nPointsList[1], nPointsList[0])

            nyTicks = nxTicks = 9

            xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
            xticks = np.linspace(iDeviations[0], iDeviations[-1], nxTicks)
            plt.xticks(xlocs,xticks)
            ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
            yticks = np.linspace(jDeviations[0], jDeviations[-1], nyTicks)
            plt.yticks(ylocs,yticks)

            if (xIndices[0] == 2):
                # Re-scale x-axis to MHz
                xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
                newDevMax = iDeviations[-1]*1e3
                xticks = np.linspace(-newDevMax, newDevMax, nxTicks)
                plt.xticks(xlocs, xticks)
            elif (xIndices[1] == 2):
                # Re-scale y-axis to MHz
                ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
                newDevMax = jDeviations[-1]*1e3
                yticks = np.linspace(-newDevMax, newDevMax, nyTicks)
                plt.yticks(ylocs, yticks)

            plt.imshow(fidelities2D, interpolation="bilinear", origin="lower")
            plt.colorbar(label="Fidelitet", orientation="vertical")

            plt.axvline(x=(nPointsList[0]-1)/2, color='k')
            plt.axhline(y=(nPointsList[1]-1)/2, color='r')
            plt.legend([iLegendStr, jLegendStr], fontsize=19, loc="upper right")
            plt.show()

            if (saveToFile):
                solDict[listName] = [[],[],[]]
                solDict[listName][0] = fidelities2D.tolist()
                solDict[listName][1] = iDeviations.tolist()
                solDict[listName][2] = jDeviations.tolist()
                dumpTojson(solDict, solutionPath)

'''
def plotEigenenergies(x=None, xName=None, N=3, simPoints=200, numOfEnergyLevels=None, useSavedPlot=False, saveTojson=False):
    if (xName is None):
        if (x is None):
            print("Must supply either the name of an x stored in a json file or a list x manually (or both)")
            return None
    else:
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getFromjson('solutions.json')
        x = solsDict[xName]['x']

    x[0] = - abs(x[0])

    if (useSavedPlot and saveTojson):
        print("Re-saving something in a json file that's already saved there is rather pointless, don't ya think?")
        saveTojson = False

    if (useSavedPlot):
        eigenEnergyDict = getFromjson('eigenenergies.json')
        energyOfEigenstate = [eigenEnergyDict[key] for key in eigenEnergyDict]
    else:
        if numOfEnergyLevels is None:
            numOfEnergyLevels = N**3
        
        energyOfEigenstate = [ [ [], [], (x, x, x), [] ] for x in range(N**3)]

        i = 0
        for q1 in range(N):
            for q2 in range(N):
                for qTB in range(N):
                    energyOfEigenstate[i][2] = (q1,q2,qTB)
                    i = i + 1

        HBareBasisComponents = getHamiltonian(x, N=N, getBBHamiltonianComps=True)
        thetas = np.linspace(-0.5, 0, simPoints)
        
        for i, theta in enumerate(thetas):
            omegaTBTh = coeffomegaTB(omegas[2], theta)
            eigenStatesAndEnergiesBareBasis = getThetaEigenstates(x, HBareBasisComponents[0]+HBareBasisComponents[1], HBareBasisComponents[2], omegaTBTh)
            order = eigenstateOrder(eigenStatesAndEnergiesBareBasis[1][0:numOfEnergyLevels], N) # eigenStatesAndEnergiesBareBasis[0][0:numOfEnergyLevels],

            for entry in order:
                _, state, energyIndex = entry
                energy = eigenStatesAndEnergiesBareBasis[0][energyIndex]

                saveEnergyAndFlux(energyOfEigenstate, state, theta, energy, energyIndex)

            statusBar((i+1)*100/simPoints)


    ############################
    # Plot energies!
    print("Plotting!")
    linestyle = ['-', '--', '-.', ':']
    labels = []
    plt.figure(figsize=(8, 7))

    for index, item in enumerate(energyOfEigenstate):
        flux, energy, state, _ = item
        if not (len(flux) == 0):
            plt.plot(flux, energy, ls=linestyle[math.floor(index / 10) % 4])
            labels.append(indexToString(state))

    plt.plot([x[0], x[0]], [-200, 200], 'r--')
    plt.xlabel('Magnetic Flux [$\Phi$]', fontsize=16)
    plt.ylabel('Energi', fontsize=16)
    plt.xlim([-0.5, 0])
    plt.ylim([-1, 100])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    leg = plt.legend(labels, bbox_to_anchor=(1.1, 1), fontsize=10, loc="upper right")
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.show()
    ############################
    if (saveTojson):
        eigenEnergyDict = {labels[i]: item for i, item in enumerate(energyOfEigenstate)}
        dumpTojson(eigenEnergyDict, 'eigenenergies.json')
'''

def plotEigenenergies(solutionPath=None, eigenenergiesPath='eigenenergies.json', N=3, simPoints=200, numOfEnergyLevels=None, useSavedPlot=False, saveToFile=False):
    solDict = getFromjson(solutionPath)
    x = [solDict['theta'], solDict['delta'], solDict['omegaPhi'], solDict['modulationTime']]
    circuitData = getCircuitData(solDict)

    x[0] = - abs(x[0])

    if (useSavedPlot and saveToFile):
        print("Re-saving something in a json file that's already saved there is rather pointless, don't ya think?")
        saveToFile = False

    if (useSavedPlot):
        eigenEnergyDict = getFromjson(eigenenergiesPath)
        energyOfEigenstate = [eigenEnergyDict[key] for key in eigenEnergyDict]
    else:
        if numOfEnergyLevels is None:
            numOfEnergyLevels = N**3
        
        energyOfEigenstate = [ [ [], [], (x, x, x), [] ] for x in range(N**3)]

        i = 0
        for q1 in range(N):
            for q2 in range(N):
                for qTB in range(N):
                    energyOfEigenstate[i][2] = (q1,q2,qTB)
                    i = i + 1

        HBareBasisComponents = getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData)
        thetas = np.linspace(-0.5, 0, simPoints)
        
        for i, theta in enumerate(thetas):
            omegaTBTh = coeffomegaTB(circuitData['frequencies'][2], theta)
            eigenStatesAndEnergiesBareBasis = getThetaEigenstates(x, HBareBasisComponents[0]+HBareBasisComponents[1], HBareBasisComponents[2], omegaTBTh)
            order = eigenstateOrder(eigenStatesAndEnergiesBareBasis[1][0:numOfEnergyLevels], N) # eigenStatesAndEnergiesBareBasis[0][0:numOfEnergyLevels],

            for entry in order:
                _, state, energyIndex = entry
                energy = eigenStatesAndEnergiesBareBasis[0][energyIndex]

                saveEnergyAndFlux(energyOfEigenstate, state, theta, energy, energyIndex)

            statusBar((i+1)*100/simPoints)


    ############################
    # Plot energies!
    print("Plotting!")
    linestyle = ['-', '--', '-.', ':']
    labels = []
    plt.figure(figsize=(8, 7))

    for index, item in enumerate(energyOfEigenstate):
        flux, energy, state, _ = item
        if not (len(flux) == 0):
            plt.plot(flux, energy, ls=linestyle[math.floor(index / 10) % 4])
            labels.append(indexToString(state))

    plt.plot([x[0], x[0]], [-200, 200], 'r--')
    plt.xlabel('Magnetic Flux [$\Phi$]', fontsize=16)
    plt.ylabel('Energi', fontsize=16)
    plt.xlim([-0.5, 0])
    plt.ylim([-1, 100])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    leg = plt.legend(labels, bbox_to_anchor=(1.1, 1), fontsize=10, loc="upper right")
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.show()
    ############################
    if (saveToFile):
        eigenEnergyDict = {labels[i]: item for i, item in enumerate(energyOfEigenstate)}
        dumpTojson(eigenEnergyDict, eigenenergiesPath)

'''
def findEigenIndex(x, eigenStateIndex=0, N=4, printResult=False):

    # Get eigenindices and dimension.
    eigIndices = getIndicesOld(N)
    D = N**3

    # Get the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x, N=N, getBBHamiltonianComps=True)

    # Calculate the tunable bus frequency when Phi=0.
    omegaTBDC = coeffomegaTB(omegas[2], 0)

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux is zero.
    eigenStatesAndEnergies = getThetaEigenstates(x, hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)
    eigenState = eigenStatesAndEnergies[1][eigenStateIndex]

    cleanEigenState = [[0] for _ in range(D)]
    for i, val in enumerate(cleanEigenState):
        cleanEigenState[i][0] = np.abs(eigenState[i].item(0))

    eigenState = Qobj(cleanEigenState, dims=[[N, N, N], [1, 1, 1]], shape=(D, 1))
    result = [100, [-1, -1, -1]]
    for q1 in range(N):
        for q2 in range(N):
            for qTB in range(N):
                diff = (Qobj(tensor(basis(N, q1), basis(N, q2), basis(N, qTB)), dims=[[N, N, N], [1, 1, 1]]) - eigenState).norm()
                if diff < result[0]:
                    result[0] = diff
                    result[1] = [q1, q2, qTB]
    if printResult:
        print(f'The eigenstate with eigen index {eigenStateIndex} is the |{result[1][0]}{result[1][1]}{result[1][2]}> state. The norm difference is {result[0]}.')

    return result
'''

def indexToString(indexTuple):
    return f'|{indexTuple[0]}{indexTuple[1]}{indexTuple[2]}>'


def saveEnergyAndFlux(itemList, state, flux, energy, energyIndex):
    for item in itemList:
        if item[2] == state:
            item[0].append(flux)
            item[1].append(energy)
            item[3].append(energyIndex)
            break

######################################################################################################################################################################