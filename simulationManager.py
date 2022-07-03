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
from plotting import *
from pathlib import Path
from threading import Thread
import simulate2QubitFrame as gui

######################################################################################################################################################################


def simulate(solutionData, guiData):
    ###################### Assign the given data to variables ######################

    energyLevels = guiData["energy-levels"]
    solutionPath = Path(guiData["solution-path"])

    eigenenergiesPath = solutionPath.parents[1]
    eigenenergiesPath = eigenenergiesPath.joinpath(Path("eigenenergies.json"))

    eigenstateIndex = guiData["eigenstateIndex"]

    ############################## Generate the plots ##############################

    if guiData["populationTransferPlot"]:
        wantEigEnPlot = True
        if not eigenenergiesPath.is_file():
            # ToDo: Make it possible to run plotEigenenergies via gui properly.
            print("No eigenenergies.json found. Generating.")
            plotEigenenergies(solutionPath, eigenenergiesPath, N=energyLevels, simPoints=500, numOfEnergyLevels=None, useSavedPlot=False, saveToFile=True)
        elif wantEigEnPlot:
            plotEigenenergies(solutionPath, eigenenergiesPath, N=energyLevels, simPoints=500, numOfEnergyLevels=None, useSavedPlot=True, saveToFile=False)
        simulatePopTransfer(solutionPath, eigenenergiesPath, initialStateIndex=eigenstateIndex, highestProjectionIndex=max(12, eigenstateIndex+5), N=energyLevels)
    elif guiData["fidelityPlot"]:
        try:
            plotFidelity(solutionPath, useSavedPlot=True, saveToFile=False)
        except:
            process = Thread(target=generateFidelityData, args=((solutionPath,)))
            process.start()
    elif guiData["stabilityPlot"]:
        # ToDo: Make sure this works and change to pdf generation later on.
        getRobustnessPlot(solutionPath, checkX0=False, checkX1=True, checkOmegaPhi=True, checkOpTime=False, nPointsList=[9], maxDevs=[5e-4, 5e-3, 5e-3, 4e0], useSavedPlot=False, saveToFile=True) # 1e-2?/1e-3 standard for Bmax/deltaMax
    
    ################################################################################


def generateFidelityData(solutionPath):
    gui.disableStartSimulationButton()
    gui.writeStatus("Fidelity data not available. Generating fidelity data!")

    process = Thread(target=plotFidelity, args=((solutionPath, False, True, False, 5)))
    process.start()
    process.join()

    gui.enableStartSimulationButton()
    gui.writeStatus("Fidelity data generated! Restart simulation!")


def getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex):
    theta = -abs(theta)
    labelsList = []
    eigenIndices = []
    knownThetas = eigenEnergyDict['Thetas']

    for label in eigenEnergyDict:
        if (label != "circuit") and (label != "Thetas"):
            knownEigInds = eigenEnergyDict[label][2]
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
        #print(usedLabels)
    return usedLabels


def plotOmegaTBvariations(solutionPath):
    solutionDict = getFromjson(solutionPath)
    circuitData = getCircuitData(solutionDict)

    if solutionDict['signalType'] == 'cos':
        x0name = 'theta'
        x1name = 'delta'
    elif solutionDict['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
    
    if solutionDict['comboSignal']:
        x = [solutionDict[x0name], solutionDict[x1name+'s'][0], solutionDict['omegaPhis'][0], solutionDict['modulationTimes'][0]]
        otherx = [solutionDict[x0name], solutionDict[x1name+'s'][1], solutionDict['omegaPhis'][1], solutionDict['modulationTimes'][1]]
        simulationTime = int(max(x[-1], otherx[-1])) + 10
        omegaTBDC = model.coeffomegaTB(circuitData['frequencies'][2], x[0], signalType=solutionDict['signalType'])
        args = {x0name: x[0], x1name+'s': [x[1], otherx[1]], 'omegaPhis': [2*np.pi*x[2], 2*np.pi*otherx[2]], 'omegatb0': circuitData['frequencies'][2], 'operationTimes': [x[3], otherx[3]], 'omegaTBTh': omegaTBDC, 'riseTime': solutionDict['riseTime']}
        if solutionDict['signalType'] == 'cos':
            omegaTBfun = model.omegaTBComboCosSinBox
        elif solutionDict['signalType'] == 'arccos':
            omegaTBfun = model.omegaTBComboArccosSinBox
    else:
        x = [solutionDict[x0name], solutionDict[x1name], solutionDict['omegaPhi'], solutionDict['modulationTime']]
        otherx = None
        simulationTime = int(x[-1]) + 10
        omegaTBDC = model.coeffomegaTB(circuitData['frequencies'][2], x[0], signalType=solutionDict['signalType'])
        args = {x0name: x[0], x1name: x[1], 'omegaPhi': 2*np.pi*x[2], 'omegatb0': circuitData['frequencies'][2], 'operationTime': x[3], 'omegaTBTh': omegaTBDC, 'riseTime': solutionDict['riseTime']}
        if solutionDict['signalType'] == 'cos':
            omegaTBfun = model.omegaTBSinBox
        elif solutionDict['signalType'] == 'arccos':
            omegaTBfun = model.omegaTBArccosSinBox

    timeStamps = np.linspace(0, simulationTime, simulationTime*20)

    oTBs = []
    for t in timeStamps:
        oTBs.append(omegaTBfun(t, args))
    
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot()

    ax.plot(timeStamps, oTBs)
    ax.plot([x[-1], x[-1]], [-50, 50], 'r--')
    if solutionDict['comboSignal']:
        ax.plot([otherx[-1], otherx[-1]], [-50, 50], 'g--')
    ax.grid()
    ax.set_ylim([-20, 20])
    ax.set_xlim([timeStamps[0], timeStamps[-1]])
    if solutionDict['comboSignal']:
        ax.legend(["Deviation in $\omega_{TB}$ [Grad/s]", "$t_{MOD}$ for gate component 1", "$t_{MOD}$ for gate component 2"], fontsize=19, loc="lower right")
    else:
        ax.legend(["Deviation in $\omega_{TB}$ [Grad/s]", "$t_{MOD}$"], fontsize=19, loc="lower right")
    #ax.title("Grindfidelitet kring $t_{MOD}$", fontsize=17)
    ax.set_xlabel("Time elapsed since gate start [ns]", fontsize=26)
    ax.set_ylabel("$\omega_{TB}$", fontsize=26)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    plt.tight_layout()
    plt.show()


def simulatePopTransfer(solutionPath, eigenenergiesPath=None, sinBoxHamiltonian=True, rotatingFrame=True, initialStateIndex=1, highestProjectionIndex=12, N=4):
    """
    This function simulates the population transfers between 
    different eigenstates of the 4-level hamiltonian from 
    McKay eq. 11.
    ---------------------------------------------------------
    INPUT:
            solutionPath (string): The path to the solution file that the simulation is be performed for.
            sinBoxHamiltonian (boolean) {Optional}: If True the amplitude of AC part of the flux signal that is applied to the tunable bus will be sinusodially modulated.
            rotatingFrame (boolean) {Optional}: If True the states will be transformed into the rotating frame after the time evolution has been completed.
            initialStateIndex (int) {Optional}: This parameter determines which eigenstate will be the initial state in the time evolution. If initialStateIndex=0,
                the eigenstate with the lowest associated energy will be the initial state.
            highestProjectionIndex (int) {Optional}: Determines how many eigenstates to project the evolving state onto. If highestProjectionIndex = n, the population of
                the n+1 lowest eigenstates are shown (since the lowest possible eigenIndex is zero).
            N (int) {Optional}: How many energy levels that should be accounted for in the simulations.
    ---------------------------------------------------------
    """
    if eigenenergiesPath is None:
        print("You have to provide a json file with the energies of all relevant eigenstates for all relevant values of Phi!")
        print("Such a file can be generated by running the function plotEigenenergies")
        return
    else:
        eigenEnergyDict = getFromjson(eigenenergiesPath)

    # Load a solution from a file and create a list x defining the control signal, converting omegaPhi from GHz to Grad/s.
    solutionDict = getFromjson(solutionPath)

    if solutionDict['signalType'] == 'cos':
        x0name = 'theta'
        x1name = 'delta'
    elif solutionDict['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
    
    if solutionDict['comboSignal']:
        x = [solutionDict[x0name], solutionDict[x1name+'s'][0], solutionDict['omegaPhis'][0], solutionDict['modulationTimes'][0]]
        otherx = [solutionDict[x0name], solutionDict[x1name+'s'][1], solutionDict['omegaPhis'][1], solutionDict['modulationTimes'][1]]
    else:
        x = [solutionDict[x0name], solutionDict[x1name], solutionDict['omegaPhi'], solutionDict['modulationTime']]
        otherx = None
    
    # Load circuit data into a dictionary and change units from GHz to Grad/s.
    circuitData = getCircuitData(solutionDict)
    eigenCircuitData = getCircuitData(eigenEnergyDict['circuit'])

    if not areCircuitDatasIdentical(circuitData, eigenCircuitData):
            print("The solution and the eigenenergy file must use the same circuit parameters!")
            return
    
    # Calculate the dimension of the tensor states and set the simulation time.
    D = N**3
    if solutionDict['comboSignal']:
        simulationTime = int(max(x[-1], otherx[-1])) + 10
    else:
        simulationTime = int(x[-1]) + 10

    # Calculate the eigenstates and eigenenergies of the bare basis hamiltonian.
    hamiltonianBareBasis = model.getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData, signalType=solutionDict['signalType'], otherx=otherx)

    # Calculate the tunable bus frequency when only the DC part of the flux is active.
    omegaTBDC = model.coeffomegaTB(circuitData['frequencies'][2], x[0], signalType=solutionDict['signalType'])

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux only has its DC part.
    eigenStatesAndEnergies = model.getThetaEigenstates(hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)

    # Get eigenindices.
    eigIndices = model.getIndices(N, eigenStatesAndEnergies[1])

    # Calculate the unitary for transforming the hamiltonian to the eigenbasis.
    eigenBasisUnitary = model.getEBUnitary(eigenStatesAndEnergies, N, D)

    # Get the hamiltonian that has a sinusodially modulated AC flux and also is in the eigenbasis.
    hamiltonian = model.getHamiltonian(x, N=N, eigEs=eigenStatesAndEnergies[0], U_e=eigenBasisUnitary, sinBoxHamiltonian=sinBoxHamiltonian, circuitData=circuitData, signalType=solutionDict['signalType'], otherx=otherx)

    # Change the simulation settings and create the timestamps for where the evolution is to be evaluated.
    options = solver.Options()
    options.nsteps = 10000
    timeStamps = np.linspace(0, simulationTime, simulationTime*5)

    # Create the initial state.
    initialState = Qobj(basis(D, initialStateIndex), dims=[[N, N, N], [1, 1, 1]])

    # Time evolve the initial state.
    if solutionDict['comboSignal']:
        args = {x0name: x[0], x1name+'s': [x[1], otherx[1]], 'omegaPhis': [2*np.pi*x[2], 2*np.pi*otherx[2]], 'omegatb0': circuitData['frequencies'][2], 'operationTimes': [x[3], otherx[3]], 'omegaTBTh': omegaTBDC, 'riseTime': solutionDict['riseTime']}
    else:
        args = {x0name: x[0], x1name: x[1], 'omegaPhi': 2*np.pi*x[2], 'omegatb0': circuitData['frequencies'][2], 'operationTime': x[3], 'omegaTBTh': omegaTBDC, 'riseTime': solutionDict['riseTime']}
    
    result = sesolve(hamiltonian, initialState, timeStamps, [], options=options, args=args)
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
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot()
    
    if solutionDict['signalType'] == 'arccos':
        theta = np.arccos(x[0] ** 2) / np.pi
    elif solutionDict['signalType'] == 'cos':
        theta = x[0]
    
    labels = getEigenstateLabels(eigenEnergyDict, theta, maxUsedIndex=highestProjectionIndex)

    linestyle = ['-', '--', '-.', ':']
    for index, values in enumerate(expectationValues):
        ax.plot(timeStamps, values, ls=linestyle[math.floor(index / 10) % 4])

    ax.grid()
    ax.set_ylim([0, 1.1])
    ax.set_xlim([0, timeStamps[-1]])
    leg = ax.legend(labels, fontsize=19, loc='center right')
    ax.set_xlabel("Time elapsed since gate start [ns]", fontsize=26)
    ax.set_ylabel("Population", fontsize=26)
    ax.tick_params(axis='both', which='major', labelsize=16)
    ax.tick_params(axis='both', which='minor', labelsize=16)
    plt.tight_layout()
    for legobj in leg.legendHandles:
        legobj.set_linewidth(5.0)
    plt.show()


def plotFidelity(solutionPath, useSavedPlot=False, saveToFile=False, plot=True, energyLevels=4):
    solutionDict = getFromjson(solutionPath)

    if solutionDict['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
    elif solutionDict['signalType'] == 'cos':
        x0name = 'theta'
        x1name = 'delta'
    
    if solutionDict['comboSignal']:
        x = [solutionDict[x0name], solutionDict[x1name+'s'][0], solutionDict['omegaPhis'][0], solutionDict['modulationTimes'][0]]
        otherx = [solutionDict[x0name], solutionDict[x1name+'s'][1], solutionDict['omegaPhis'][1], solutionDict['modulationTimes'][1]]
    else:
        x = [solutionDict[x0name], solutionDict[x1name], solutionDict['omegaPhi'], solutionDict['modulationTime']]
        otherx = None
    
    circuitData = getCircuitData(solutionDict)

    if useSavedPlot:
        F = solutionDict['fidelitiesAtTimes']
        times = solutionDict['times']

        if len(F) == 0:
            raise Exception("Fidelities not previously generated.")
    else:
        indices = np.linspace(-116, -1, 116).astype(int)
        gateType = solutionDict['gateType']
        
        F, times = model.getGateFidelity(x, gateType=gateType, N=energyLevels, tIndices=indices, circuitData=circuitData, signalType=solutionDict['signalType'], otherx=otherx)

    if plot:
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot()

        ax.plot(times, F)
        ax.plot([x[-1], x[-1]], [0, 1], 'r--')
        if solutionDict['comboSignal']:
            ax.plot([otherx[-1], otherx[-1]], [0, 1], 'g--')
        ax.grid()
        ax.set_ylim([0.99, 1])
        ax.set_xlim([times[0], times[-1]])
        if solutionDict['comboSignal']:
            ax.legend(["Fidelity", "$t_{MOD}$ for gate component 1", "$t_{MOD}$ for gate component 2"], fontsize=19, loc="lower right")
        else:
            ax.legend(["Fidelity", "$t_{MOD}$"], fontsize=19, loc="lower right")
        #ax.title("Grindfidelitet kring $t_{MOD}$", fontsize=17)
        ax.set_xlabel("Time elapsed since gate start [ns]", fontsize=26)
        ax.set_ylabel("Fidelity", fontsize=26)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)
        plt.tight_layout()
        plt.show()

    if (saveToFile == True):
        solutionDict['times'] = times
        solutionDict['fidelitiesAtTimes'] = F
        solutionDict['gateFidelity'] = F[-76]
        dumpTojson(solutionDict,solutionPath)


def getRobustnessPlot(solutionPath, checkX0=False, checkX1=False, checkOmegaPhi=False, checkOpTime=False, nPointsList=[9], maxDevs=[5e-4, 1e-3, 2e-3, 4e0], useSavedPlot=False, saveToFile=False):
    solutionDict = getFromjson(solutionPath)
    
    if solutionDict['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        x1name = 'acAmplitude'
        axlabel_x0 = "Deviation from optimal $A$ [1]"
        axlabel_x1 = "Deviation from optimal $B$ [1]"
    elif solutionDict['signalType'] == 'cos':
        x0name = 'theta'
        x1name = 'delta'
        axlabel_x0 = "Deviation from optimal $\Theta$ [$\Phi_0$]"
        axlabel_x1 = "Deviation from optimal $\delta$ [$\Phi_0$]"
    
    x = [solutionDict[x0name], solutionDict[x1name], solutionDict['omegaPhi'], solutionDict['modulationTime']]
    circuitData = getCircuitData(solutionDict)

    if solutionDict['signalType'] == 'arccos':
        legendStr_x0 = "$A = %.4f$" %x[0]
        legendStr_x1 = "$B = %.4f$" %x[1]
    elif solutionDict['signalType'] == 'cos':
        legendStr_x0 = "$\Theta = %.4f$" %x[0]
        legendStr_x1 = "$\delta = %.4f$" %x[1]

    if (useSavedPlot and saveToFile):
        print("Re-saving data in a file where that data is already saved is rather pointless, don't ya think?")
        saveToFile = False

    checkList = [checkX0, checkX1, checkOmegaPhi, checkOpTime]
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

        axlabels = [axlabel_x0, axlabel_x1, "Deviation from optimal $\omega_\Phi$ [MHz]", "Deviation from optimal $t_{MOD}$ [ns]"]
        legendStrs = [legendStr_x0, legendStr_x1, "$\omega_\Phi = %.1f$" %(x[2]*1000), "$t_{MOD} = %.1f$" %x[3]]
        # maxDevs = [3e-4, 1e-3, 2*np.pi * 8e-3, 40e0] # Testing extended bounds, looking for chevrony stuff

        xIndices = []
        if (checkX0):
            xIndices.append(0)
        if (checkX1):
            xIndices.append(1)
        if (checkOmegaPhi):
            xIndices.append(2)
        if (checkOpTime):
            xIndices.append(3)

        if (sum(checkList) == 1):
            xIndex = xIndices[0]
            if (useSavedPlot):
                deviations = np.array(solutionDict['deviations'][xIndex])
                fidelities = solutionDict['fidelities1D_Dev'][xIndex]
            else:
                deviations = np.linspace(-maxDevs[xIndex], maxDevs[xIndex], nPointsList[0])
                fidelities = []
                for i, d in enumerate(deviations):
                    xDev[xIndex] = x[xIndex] + d
                    fidelity, _ = model.getGateFidelity(xDev, gateType=solutionDict['gateType'], N=4, tIndices=[-76], circuitData=circuitData, signalType=solutionDict['signalType'])
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
            plt.legend(["Fidelity", legendStrs[xIndex]], fontsize=19, loc="upper right")
            plt.grid()
            plt.xlabel(axlabels[xIndex], fontsize=26)
            plt.ylabel("Fidelity", fontsize=26)
            plt.xlim([deviations[0], deviations[-1]])
            plt.ylim([0.99, 1])
            plt.tight_layout()
            plt.show()

            if (saveToFile):
                solutionDict['deviations'][xIndex] = deviations.tolist()
                solutionDict['fidelities1D_Dev'][xIndex] = fidelities
                dumpTojson(solutionDict, solutionPath)
        elif (sum(checkList) == 2):
            if (len(nPointsList) == 1):
                nPointsList.append(nPointsList[0])
            
            plt.xlabel(axlabels[xIndices[0]], fontsize=26)
            plt.ylabel(axlabels[xIndices[1]], fontsize=26)
            iLegendStr = legendStrs[xIndices[0]]
            jLegendStr = legendStrs[xIndices[1]]

            if solutionDict['signalType'] == 'arccos':
                xParts = ['dcAmplitude', 'acAmplitude', 'omegaPhi', 'opTime']
            elif solutionDict['signalType'] == 'cos':
                xParts = ['Theta', 'delta', 'omegaPhi', 'opTime']
            listName = "fidelities2D_" + xParts[xIndices[0]] + "_" + xParts[xIndices[1]]

            if (useSavedPlot):
                fidelities2D = np.array(solutionDict[listName][0])
                iDeviations = np.array(solutionDict[listName][1])
                jDeviations = np.array(solutionDict[listName][2])
            else:
                iDeviations = np.linspace(-maxDevs[xIndices[0]], maxDevs[xIndices[0]], nPointsList[0])
                jDeviations = np.linspace(-maxDevs[xIndices[1]], maxDevs[xIndices[1]], nPointsList[1])
                
                fidelities = []
                for j, jDev in enumerate(jDeviations):
                    xDev[xIndices[1]] = x[xIndices[1]] + jDev
                    for i, iDev in enumerate(iDeviations):
                        xDev[xIndices[0]] = x[xIndices[0]] + iDev
                        fidelity, _ = model.getGateFidelity(xDev, gateType=solutionDict['gateType'], N=4, tIndices=[-76], circuitData=circuitData, signalType=solutionDict['signalType'])
                        fidelities.append(fidelity[0])
                        statusBar((j*nPointsList[0] + (i+1))*100/(nPointsList[0]*nPointsList[1]))
                fidelities2D = np.array(fidelities).reshape(nPointsList[1], nPointsList[0])

            nyTicks = nxTicks = 5

            xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
            xticks = np.around(np.linspace(iDeviations[0], iDeviations[-1], nxTicks), 6) # This is suboptimal
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
                solutionDict[listName] = [[],[],[]]
                solutionDict[listName][0] = fidelities2D.tolist()
                solutionDict[listName][1] = iDeviations.tolist()
                solutionDict[listName][2] = jDeviations.tolist()
                dumpTojson(solutionDict, solutionPath)


def plotEigenenergies(solutionPath, eigenenergiesPath, N=3, simPoints=200, numOfEnergyLevels=None, useSavedPlot=False, saveToFile=False):
    solutionDict = getFromjson(solutionPath)

    if solutionDict['signalType'] == 'arccos':
        x0name = 'dcAmplitude'
        #x1name = 'acAmplitude'
    elif solutionDict['signalType'] == 'cos':
        x0name = 'theta'
        #x1name = 'delta'
    
    #x = [solutionDict[x0name], solutionDict[x1name], solutionDict['omegaPhi'], solutionDict['modulationTime']]
    circuitData = getCircuitData(solutionDict)

    if solutionDict['signalType'] == 'arccos':
        th = - np.arccos(solutionDict[x0name] ** 2) / np.pi
    elif solutionDict['signalType'] == 'cos':
        th = - abs(solutionDict[x0name])

    if (useSavedPlot and saveToFile):
        print("Re-saving something in a json file that's already saved there is rather pointless, don't ya think?")
        saveToFile = False


    if (useSavedPlot):
        eigenEnergyDict = getFromjson(eigenenergiesPath)
        eigenCircuitData = getCircuitData(eigenEnergyDict['circuit'])
        if not areCircuitDatasIdentical(circuitData, eigenCircuitData):
            print("The solution and the eigenenergy file must use the same circuit parameters!")
            return
        energyOfEigenstate = [eigenEnergyDict[key] for key in eigenEnergyDict if (key != "Thetas") and (key != "circuit")]
        thetas = eigenEnergyDict['Thetas']
    else:
        if numOfEnergyLevels is None:
            numOfEnergyLevels = N**3
        
        energyOfEigenstate = [ [ (k, k, k), [], [] ] for k in range(N**3)]

        i = 0
        for q1 in range(N):
            for q2 in range(N):
                for qTB in range(N):
                    energyOfEigenstate[i][0] = (q1,q2,qTB)
                    i = i + 1

        HBareBasisComponents = model.getHamiltonian([th], N=N, getBBHamiltonianComps=True, circuitData=circuitData, signalType=solutionDict['signalType'])
        thetas = np.linspace(-0.5, 0, simPoints)
        
        for i, theta in enumerate(thetas):
            omegaTBTh = model.coeffomegaTB(circuitData['frequencies'][2], theta, signalType=solutionDict['signalType'])
            eigenStatesAndEnergiesBareBasis = model.getThetaEigenstates(HBareBasisComponents[0]+HBareBasisComponents[1], HBareBasisComponents[2], omegaTBTh)
            order = model.eigenstateOrder(eigenStatesAndEnergiesBareBasis[1][0:numOfEnergyLevels], N) # eigenStatesAndEnergiesBareBasis[0][0:numOfEnergyLevels],

            for entry in order:
                _, state, energyIndex = entry
                energy = eigenStatesAndEnergiesBareBasis[0][energyIndex]

                saveEnergyAndFlux(energyOfEigenstate, state, energy, energyIndex)

            statusBar((i+1)*100/simPoints)


    ############################
    # Plot energies!
    print("Plotting!")
    linestyle = ['-', '--', '-.', ':']
    labelsOrdered = []
    plt.figure(figsize=(8, 7))

    # Find the index in thetas corresponding to (just above) th
    thetaIndex = 0
    for i, theta in enumerate(thetas):
        if theta > th:
            thetaIndex = i
            break
    
    # Reorder energyOfEigenstate according to rising energy at Phi = th
    energyOfEigenstateOrdered = sorted(energyOfEigenstate, key=lambda stateEnergiesIndices: stateEnergiesIndices[2][thetaIndex])

    for index, item in enumerate(energyOfEigenstateOrdered):
        state, energy, _ = item
        plt.plot(thetas, [e/(2*np.pi) for e in energy], ls=linestyle[math.floor(index / 10) % 4])
        labelsOrdered.append(indexToString(state))

    plt.plot([th, th], [-30, 30], 'r--')
    plt.xlabel('Magnetic Flux [$\Phi$]', fontsize=16)
    plt.ylabel('Energy [GHz]', fontsize=16)
    plt.xlim([-0.5, 0])
    plt.ylim([-1, 15])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    leg = plt.legend(labelsOrdered, bbox_to_anchor=(1.1, 1), fontsize=10, loc="upper right")
    for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
    plt.show()
    ############################
    if (saveToFile):
        eeDict = {}
        eeDict['circuit'] = getCircuitData(solutionDict, convertToGradPerSecond=False)
        eeDict['Thetas'] = thetas.tolist()
        updateDict = {indexToString(item[0]): item for item in energyOfEigenstate}
        eeDict.update(updateDict)
        dumpTojson(eeDict, eigenenergiesPath)


def indexToString(indexTuple):
    return f'|{indexTuple[0]}{indexTuple[1]}{indexTuple[2]}>'


def saveEnergyAndFlux(itemList, state, energy, energyIndex):
    for item in itemList:
        if item[0] == state:
            item[1].append(energy)
            item[2].append(energyIndex)
            break