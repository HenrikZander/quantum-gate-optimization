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

from qutip import *
import model
from dataManager import *
import numpy as np
import matplotlib.pyplot as plt
import math

######################################################################################################################################################################

def simulate():
    pass


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


def plotFidelity(solutionPath, useSavedPlot=False, saveToFile=False):
    solutionDict = getFromjson(fileName=solutionPath)
    x = (solutionDict['theta'], solutionDict['delta'], solutionDict['omegaPhi'], solutionDict['modulationTime'])
    
    circuitData = getCircuitData(solutionDict)

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
    hamiltonianBareBasis = model.getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData)

    # Calculate the tunable bus frequency when only the DC part of the flux is active.
    omegaTBDC = model.coeffomegaTB(circuitData['frequencies'][2], x[0])

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux only has its DC part.
    eigenStatesAndEnergies = model.getThetaEigenstates(x, hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)

    # Get eigenindices.
    eigIndices = model.getIndices(N, eigenStatesAndEnergies[1])

    # Calculate the unitary for transforming the hamiltonian to the eigen basis.
    eigenBasisUnitary = model.getEBUnitary(x, eigenStatesAndEnergies, N, D)

    # Get the hamiltonian that has a sinusodially modulated AC flux and also is in the eigen basis.
    hamiltonian = model.getHamiltonian(x, N=N, eigEs=eigenStatesAndEnergies[0], U_e=eigenBasisUnitary, sinStepHamiltonian=sinStepHamiltonian, circuitData=circuitData)

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
            unitaryRotatingFrame = model.getRFUnitary(Hrot, t)
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