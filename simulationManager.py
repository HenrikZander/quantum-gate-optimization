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
        F, times = model.getGateFidelity(x, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=indices, circuitData=circuitData)

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
                    fidelity, _ = model.getGateFidelity(xDev, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
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
                        fidelity, _ = model.getGateFidelity(xDev, N=4, iSWAP=iSWAP, SWAP=SWAP, CZ=CZ, tIndices=[-76], circuitData=circuitData)
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

        HBareBasisComponents = model.getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData)
        thetas = np.linspace(-0.5, 0, simPoints)
        
        for i, theta in enumerate(thetas):
            omegaTBTh = model.coeffomegaTB(circuitData['frequencies'][2], theta)
            eigenStatesAndEnergiesBareBasis = model.getThetaEigenstates(x, HBareBasisComponents[0]+HBareBasisComponents[1], HBareBasisComponents[2], omegaTBTh)
            order = model.eigenstateOrder(eigenStatesAndEnergiesBareBasis[1][0:numOfEnergyLevels], N) # eigenStatesAndEnergiesBareBasis[0][0:numOfEnergyLevels],

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


def indexToString(indexTuple):
    return f'|{indexTuple[0]}{indexTuple[1]}{indexTuple[2]}>'


def saveEnergyAndFlux(itemList, state, flux, energy, energyIndex):
    for item in itemList:
        if item[2] == state:
            item[0].append(flux)
            item[1].append(energy)
            item[3].append(energyIndex)
            break