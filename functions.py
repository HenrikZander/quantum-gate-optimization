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

# Last modified: 22 May 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy
from numba import njit

from plotting import *
from model import *
from variables import *

import json
from datetime import datetime

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
# Callback functions

def callbackDE(x,convergence=None):
    """
    Callback function for the Differential Evolution optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    print(f'Num of iterations: {i+1}')
    print(f'The currently best minimum that the Differential Evolution algorithm has found has a convergence of {convergence} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = [(x,convergence)]
        saveResToFile(result, "Differential Evolution", i, passedTime, algorithmDE=True)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def callbackDA(x, fun, context):
    """
    Callback function for the Dual Anneling optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    bestResults = evaluateResult(x, fun, bestResults)
    print(f'Num of minima found: {i+1}')
    print(f'The Dual Anneling algorithm gave a minimum of {fun} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = []
        for i in range(len(bestResults[1])):
            result.append((bestResults[0][i],bestResults[1][i]))
        saveResToFile(result, "Dual Anneling", i, passedTime)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def callbackSHG(x):
    """
    Callback function for the Simplical Homology Global optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    print(f'Num of iterations: {i+1}')
    print(f'The current point of examination by the Simplicial Homology Global algorithm is {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False

######################################################################################################################################################################

def findMinimum(costFunction, bounds, runSHG=True, runDA=True, runDE=True):
    """
    This function aims to find the global minimum of the specified cost 
    function, within the parameter bounds that has been specified. A 
    number of different algorithms are implemented and they can be used 
    consecutively during one function call.
    ---------------------------------------------------------
    INPUT:
            costFunction (function pointer): The cost function that is to be minimized. The input parameters to the cost function should be specified as an array x.
            bounds (array(tuple(int))): The parameter bounds that the optimizer will do it's best to not leave. The length of this array determines the expected dimension of the cost function input x.
            runSHG (boolean) {Optional}: If True the function will use the Simplicial Homology Global algorithm to minimize the cost function.
            runDA (boolean) {Optional}: If True the function will use the Dual Anneling algorithm to minimize the cost function.
            runDE (boolean) {Optional}: If True the function will use the Differential Evolution algorithm to minimize the cost function.
    ---------------------------------------------------------
    OUTPUT:
            result (array(scipy.optimize.OptimizeResult)): An array of Scipy result structures. See Scipy documentation for more information.
    ---------------------------------------------------------
    """
    
    global startTime 
    
    #Initial formatting tasks for the result that will be printed in the terminal.
    print("")
    message = "##################################################\n"
    result = []
    algorithmsUsed = []
    runtime = []

    #Optimization using the Simplicial Homology Global algorithm.
    if runSHG:
        startTime = time.time()
        resSHG = scipy.optimize.shgo(costFunction, bounds, iters=4, callback=callbackSHG)
        timeSHG = time.time() - startTime
        message += f'The optimizaton using the \"Simplicial Homology Global\"-algorithm took {round(timeSHG,2)}s to execute and ended on a minimum of {resSHG.fun} at the point {resSHG.x}.\n'
        message += f'Function evaluations performed: {resSHG.nfev}\n'
        result.append(resSHG)
        algorithmsUsed.append("Simplicial Homology Global")
        runtime.append(timeSHG)
    
    #Optimization using the Dual Annealing algorithm. 
    if runDA:
        startTime = time.time()
        resDA = scipy.optimize.dual_annealing(costFunction, bounds, callback=callbackDA)
        timeDA = time.time() - startTime
        message += f'The optimizaton using the \"Dual Annealing\"-algorithm took {round(timeDA,2)}s to execute and ended on a minimum of {resDA.fun} at the point {resDA.x}.\n'
        message += f'Function evaluations performed: {resDA.nfev}\n'
        result.append(resDA)
        algorithmsUsed.append("Dual Annealing")
        runtime.append(timeDA)
    
    #Optimization using the Differential Evolution algorithm.
    if runDE:
        startTime = time.time()
        resDE = scipy.optimize.differential_evolution(costFunction, bounds, callback=callbackDE, workers=-1, updating='deferred', maxiter=3000)
        timeDE = time.time() - startTime
        message += f'The optimizaton using the \"Differential Evolution\"-algorithm took {round(timeDE,2)}s to execute and ended on a minimum of {resDE.fun} at the point {resDE.x}.\n'
        message += f'Function evaluations performed: {resDE.nfev}\n'
        result.append(resDE)
        algorithmsUsed.append("Differential Evolution")
        runtime.append(timeDE)
    
    print("")        
    print(message + "##################################################")
    saveAllFinalResults(result, algorithmsUsed, runtime)
    return result

######################################################################################################################################################################
# json helper functions

def getjsonDict(fileName):
    with open(fileName, 'r') as jsonFile:
        jsonDict = json.load(jsonFile)
    return jsonDict

def dumpjsonDict(dict, fileName):
    with open(fileName, 'w') as jsonFile:
        json.dump(dict, jsonFile, ensure_ascii=False, indent=4)

def getSolutionNameList():
    solsDict = getjsonDict('solutions.json')
    return [xName for xName in solsDict]

def createSolName(ymd,gateType,solNumberStr,usedNlvls):
    return "x_" + ymd + "_" + gateType + "_" + solNumberStr + "_" + usedNlvls

def addNewSolution(x, gateType, N, solNumber=1, creationTime=datetime.today()):
    ymd = creationTime.strftime('%Y%m%d')[2:]
    creationTime = creationTime.strftime('%Y-%m-%d %H:%M:%S')
    usedNlvls = str(N) + "lvl"
    solName = createSolName(ymd,gateType,str(solNumber),usedNlvls)

    # print(solName)

    solsDict = getjsonDict('solutions.json')
    while True:
        if solName not in solsDict:
            solsDict[solName] = {
                "creationTime": creationTime,
                "gateType": gateType,
                "nOptimizationLvls": N,
                'x': x,
                'gateFidelity': None,
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

            dumpjsonDict(solsDict, 'solutions.json')
            return None
        else:
            if (solsDict[solName]['creationTime'] == creationTime):
                print("Can't add solution: Solution already exists in solutions.json")
                return None
            solNumber += 1
            solName = createSolName(ymd,gateType,str(solNumber),usedNlvls)
    

######################################################################################################################################################################
# Simulation functions

def simulateHamiltonian(x=None, sinStepHamiltonian=True, rotatingFrame=False, initialStateIndex=1, highestProjectionIndex=8, N=4, xName=None, saveFidelity=False):
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
    
    if (xName is not None):
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getjsonDict('solutions.json')
        x = solsDict[xName]['x']
    elif (x is None):
        print("At least one of the x and xName input fields needs to not be None!")
        return None


    # Calculate the dimension of the tensor states and set the simulation time.
    D = N**3
    simulationTime = int(x[-1]) + 10

    # Calculate the eigenstates and eigenenergies of the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x,N=N,getBBHamiltonianComps=True)

    # Calculate the tunable bus frequency when only the DC part of the flux is active.
    omegaTBDC = coeffomegaTB(omegas[2],x[0])

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
    timeStamps = np.linspace(0,simulationTime,simulationTime*5)
    
    # Create the initial state.
    initialState = Qobj(basis(D,initialStateIndex),dims=[[N,N,N],[1,1,1]])

    # Time evolve the initial state.
    result = sesolve(hamiltonian, initialState, timeStamps, [], options=options, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': omegas[2], 'operationTime': x[3], 'omegaTBTh': omegaTBDC})
    states = result.states

    # Transform into the rotating frame.
    Hrot = np.array(hamiltonian[0])
    Hrot[eigIndices[-1],eigIndices[-1]] = Hrot[eigIndices[-2],eigIndices[-2]] + Hrot[eigIndices[-3],eigIndices[-3]] - Hrot[eigIndices[0],eigIndices[0]]
    Hrot = Qobj(Hrot,dims=[[N,N,N],[N,N,N]])

    if rotatingFrame:
        for i, t in enumerate(timeStamps):
            unitaryRotatingFrame = getRFUnitary(Hrot, t)
            states[i] = unitaryRotatingFrame*states[i]
    
    # Calculate the expectation values a projection operator.
    projectionOperators = []
    #pIndexes = [0,1,2,5,highestProjectionIndex]

    for index in range(highestProjectionIndex+1):
        operator = Qobj(basis(D,index),dims=[[N,N,N],[1,1,1]]) * Qobj(basis(D,index),dims=[[N,N,N],[1,1,1]]).dag()
        projectionOperators.append(operator)
    
    expectationValues = expect(projectionOperators, states)
    
    #Calculate gate fidelity for both iSWAP and CZ.
    gateFidelity_iSWAP, timesiSWAP = getGateFidelity(x,N=N,wantiSWAP=True)
    gateFidelity_CZ, timesCZ = getGateFidelity(x,N=N,wantCZ=True)
    
    # Print fidelity
    print(f'################################################\n\nGate fidelity for iSWAP: {gateFidelity_iSWAP}.\n\nGate fidelity for CZ: {gateFidelity_CZ}.\n\n################################################')
    
    # Plot the expectation values.
    plt.figure(figsize=(8,7))
    #labels = ["|000>", "|010>", "|100>", "|001>", "|020>", "|110>", "|011>", "|200>", "|101>"]
    # labels = []
    labels = ['|000>','|010>','|100>','|001>','|020>','|110>','|200>','|011>','|101>','|002>','|003>','|120>']
    
    for index, values in enumerate(expectationValues):
        plt.plot(timeStamps, values)
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


def plotFidelity(x=None, wantiSWAP=False, wantCZ=False, xName=None, saveResults=False):
    if (xName is not None):
        if (x is not None):
            print("Warning: Ignoring manually supplied x")
        solsDict = getjsonDict('solutions.json')
        F = solsDict[xName]['fidelitiesAtTimes']
        times = solsDict[xName]['times']
        x = solsDict[xName]['x']
        if (F is None) or (times is None):
            indices = np.linspace(-116, -1, 116).astype(int)
            F, times = getGateFidelity(x, N=4, wantiSWAP=wantiSWAP, wantCZ=wantCZ, tIndices=indices)
    else:
        if (x is None):
            print("Must supply either the name of an x stored in a json file or a list x manually (or both)")
            return None
        else:
            indices = np.linspace(-116, -1, 116).astype(int)
            F, times = getGateFidelity(x, N=4, wantiSWAP=wantiSWAP, wantCZ=wantCZ, tIndices=indices)
            
    plt.figure(figsize=(8,7))
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

    if (saveResults == True):
        if (xName is None):
            print("Warning: Couldn't save results. To save results, you need to specify a solution listed in solutions.json")
        else:
            solsDict[xName]['times'] = times
            solsDict[xName]['fidelitiesAtTimes'] = F
            solsDict[xName]['gateFidelity'] = F[-76]
            dumpjsonDict(solsDict,'solutions.json')


def deltaPulsePlot():
    operationTime = 75
    x = np.linspace(0,operationTime, 500)
    y = []

    for time in x:
        y.append(sinBox(time,operationTime))

    plt.plot(x, y)
    plt.plot([25, 25], [0, 1.1], 'r--', [operationTime-25, operationTime-25], [0, 1.1], 'r--')
    plt.xlabel('Tid [ns]', fontsize=16)
    plt.ylabel('Del av full amplitud, $\delta(t) / \delta_0$', fontsize=16)
    plt.xlim([0, operationTime])
    plt.ylim([0, 1.1])
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.annotate('', xy=(0,0.5), xytext=(25,0.5), arrowprops=dict(arrowstyle='<->'))
    plt.annotate("$t_{Stig}$", xy=(6,0.53), fontsize=18)
    plt.annotate('', xy=(operationTime-25,0.5), xytext=(operationTime,0.5), arrowprops=dict(arrowstyle='<->'))
    plt.annotate("$t_{Fall}$", xy=(operationTime-10,0.53), fontsize=18)
    plt.show()


def getRobustnessPlot(x, wantiSWAP=False, wantCZ=False, checkTheta=False, checkDelta=False, checkOmegaPhi=False, checkOpTime=False, nPointsList=[9], maxDevs=[3e-4, 1e-3, 2*np.pi * 2e-3, 4e0], saveResults=False):
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
        plt.figure(figsize=(11,7))
        plt.xticks(fontsize=16)
        plt.yticks(fontsize=16)

        axlabels = ["Avvikelse från hittat $\Theta$ [$\Phi_0$]", "Avvikelse från hittat $\delta$ [$\Phi_0$]", "Avvikelse från hittat $\omega_\Phi$ [MHz]", "Avvikelse från hittat $t_{MOD}$ [ns]"]
        legendStrs = ["$\Theta = %.4f$" %x[0], "$\delta = %.4f$" %x[1], "$\omega_\Phi = %.1f$" %(x[2]/(2*np.pi)*1000), "$t_{MOD} = %.1f$" %x[3]]
        # maxDevs = [3e-4, 1e-3, 2*np.pi * 8e-3, 40e0] # Testing extended bounds, looking for chevrons

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
            deviations = np.linspace(-maxDevs[xIndex], maxDevs[xIndex], nPointsList[0])

            if (xIndex == 2):
                # Re-scale x-axis to MHz
                nTicks = 9
                locs = np.linspace(-maxDevs[2], maxDevs[2], nTicks)
                newDevMax = maxDevs[2]/(2*np.pi)*1e3
                newTicks = np.linspace(-newDevMax,newDevMax,nTicks)
                plt.xticks(locs,newTicks)

            fidelities = []
            for i, d in enumerate(deviations):
                xDev[xIndex] = x[xIndex] + d
                fidelity, _ = getGateFidelity(xDev, N=4, wantiSWAP=wantiSWAP, wantCZ=wantCZ, tIndices=[-76])
                fidelities.append(fidelity)
                statusBar((i+1)*100/nPointsList[0])
            
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
        elif (sum(checkList) == 2):
            if (len(nPointsList) == 1):
                nPointsList.append(nPointsList[0])
            
            iDeviations = np.linspace(-maxDevs[xIndices[0]], maxDevs[xIndices[0]], nPointsList[0])
            jDeviations = np.linspace(-maxDevs[xIndices[1]], maxDevs[xIndices[1]], nPointsList[1])
            plt.xlabel(axlabels[xIndices[0]], fontsize=26)
            plt.ylabel(axlabels[xIndices[1]], fontsize=26)
            iLegendStr = legendStrs[xIndices[0]]
            jLegendStr = legendStrs[xIndices[1]]

            fidelities = []
            for j, jDev in enumerate(jDeviations):
                xDev[xIndices[1]] = x[xIndices[1]] + jDev
                for i, iDev in enumerate(iDeviations):
                    xDev[xIndices[0]] = x[xIndices[0]] + iDev
                    fidelity, _ = getGateFidelity(xDev, N=4, wantiSWAP=wantiSWAP, wantCZ=wantCZ, tIndices=[-76])
                    fidelities.append(fidelity[0])
                    statusBar((j*nPointsList[0] + (i+1))*100/(nPointsList[0]*nPointsList[1]))
            fidelities2D = np.array(fidelities).reshape(nPointsList[1], nPointsList[0])

            nyTicks = nxTicks = 9

            xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
            xticks = np.linspace(-maxDevs[xIndices[0]], maxDevs[xIndices[0]], nxTicks)
            plt.xticks(xlocs,xticks)
            ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
            yticks = np.linspace(-maxDevs[xIndices[1]], maxDevs[xIndices[1]], nyTicks)
            plt.yticks(ylocs,yticks)

            if (xIndices[0] == 2):
                # Re-scale x-axis to MHz
                xlocs = np.linspace(0,nPointsList[0]-1,nxTicks)
                newDevMax = maxDevs[2]/(2*np.pi)*1e3
                xticks = np.linspace(-newDevMax, newDevMax, nxTicks)
                plt.xticks(xlocs,xticks)
            elif (xIndices[1] == 2):
                # Re-scale y-axis to MHz
                ylocs = np.linspace(0,nPointsList[1]-1,nyTicks)
                newDevMax = maxDevs[2]/(2*np.pi)*1e3
                yticks = np.linspace(-newDevMax, newDevMax, nyTicks)
                plt.yticks(ylocs,yticks)

            plt.imshow(fidelities2D, interpolation="bilinear", origin="lower")
            plt.colorbar(label="Fidelitet", orientation="vertical")

            plt.axvline(x=(nPointsList[0]-1)/2, color='k')
            plt.axhline(y=(nPointsList[1]-1)/2, color='r')
            plt.legend([iLegendStr, jLegendStr], fontsize=19, loc="upper right")
            plt.show()


def indexToString(indexTuple):
    return f'|{indexTuple[0]}{indexTuple[1]}{indexTuple[2]}>'


def saveEnergyAndFlux(itemList, state, flux, energy):
    for item in itemList:
        if item[2] == state:
            item[0].append(flux)
            item[1].append(energy)
            break


def plotEigenenergies(x, N=3, simPoints=200, numOfEnergyLevels=None):
    if numOfEnergyLevels is None:
        numOfEnergyLevels = N**3
    
    energyOfEigenstate = [ [ [],[],(x, x, x) ] for x in range(N**3)]

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

            saveEnergyAndFlux(energyOfEigenstate, state, theta, energy)

        statusBar((i+1)*100/simPoints)

    ############################
    #Plot energies!
    print("Plotting!")
    linestyle = [':', '--', '-.']
    labels = []
    plt.figure(figsize=(8,7))

    for index, item in enumerate(energyOfEigenstate):
        flux, energy, state = item
        if not (len(flux) == 0):
            plt.plot(flux, energy, ls=linestyle[index % 3])
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


def findEigenIndex(x, eigenStateIndex=0, N=4, printResult=False):

    # Get eigenindices and dimension.
    eigIndices = getIndicesOld(N)
    D = N**3

    # Get the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x,N=N,getBBHamiltonianComps=True)

    # Calculate the tunable bus frequency when Phi=0.
    omegaTBDC = coeffomegaTB(omegas[2],0)

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux is zero.
    eigenStatesAndEnergies = getThetaEigenstates(x, hamiltonianBareBasis[0]+hamiltonianBareBasis[1], hamiltonianBareBasis[2], omegaTBDC)
    eigenState = eigenStatesAndEnergies[1][eigenStateIndex]

    cleanEigenState = [[0] for _ in range(D)]
    for i, val in enumerate(cleanEigenState):
        cleanEigenState[i][0] = np.abs(eigenState[i].item(0))

    eigenState = Qobj(cleanEigenState,dims=[[N,N,N],[1,1,1]], shape=(D,1))
    result = [100,[-1,-1,-1]]
    for q1 in range(N):
        for q2 in range(N):
            for qTB in range(N):
                diff = ( Qobj(tensor(basis(N,q1),basis(N,q2),basis(N,qTB)),dims=[[N,N,N],[1,1,1]]) - eigenState ).norm()
                if diff < result[0]:
                    result[0] = diff
                    result[1] = [q1,q2,qTB]
    if printResult:
        print(f'The eigenstate with eigen index {eigenStateIndex} is the |{result[1][0]}{result[1][1]}{result[1][2]}> state. The norm difference is {result[0]}.')

    return result


######################################################################################################################################################################
# Cost function definitions


def averageFidelity(F, gateTimeWeight=2):
    F[0] = gateTimeWeight*F[0]
    return np.sum(F)/(len(F)+gateTimeWeight-1)


def costiSWAP2(x):
    F, _ = getGateFidelity(x,N=2,wantiSWAP=True)
    return -averageFidelity(F)


def costiSWAP3(x):
    F, _ = getGateFidelity(x,N=3,wantiSWAP=True)
    return -averageFidelity(F)


def costiSWAP4(x):
    F, _ = getGateFidelity(x,N=4,wantiSWAP=True)
    return -averageFidelity(F)


def costCZ2(x):
    F, _ = getGateFidelity(x,N=2,wantCZ=True)
    return -averageFidelity(F)


def costCZ3(x):
    F, _ = getGateFidelity(x,N=3,wantCZ=True)
    return -averageFidelity(F)


def costCZ4(x):
    F, _ = getGateFidelity(x,N=4,wantCZ=True)
    return -averageFidelity(F)


######################################################################################################################################################################
# Optimize gate function


def optimizeGate(iSWAP=False,CZ=False,energyLevels=2, timeoutAlgorithm=55000, maxAllowedGateTime=240, runSHG=False, runDA=False, runDE=False, wantTradGate=False, wantCZ_20=False):
    """
    The function tries to optimize the choosen gate for the
    choosen parameters, using the optimization algorithms 
    that the user decide would best fit the problem.
    ---------------------------------------------------------
    INPUT:
            energyLevels (int) {Optional}: How many energy levels that should be accounted for in the simulations.
            timeoutAlgorithm (int) {Optional}: How many seconds each optimization algorithm can run, before a forced termination is initiated.
            maxAllowedGateTime (int) {Optional}: The longest gate time we will allow a solution to have.
            runSHG (boolean) {Optional}: If True, the function will use the Simplicial Homology Global algorithm to optimize the gate.
            runDA (boolean) {Optional}: If True, the function will use the Dual Anneling algorithm to optimize the gate.
            runDE (boolean) {Optional}: If True, the function will use the Differential Evolution algorithm to optimize the gate.
            wantTradGate (boolean) {Optional}: If True, the function will restrict the parameter space to focus on ''traditional'', 
                well-understood iSWAP/CZ gates.
            wantCZ_20 (boolean) {Optional}: This boolean is only relevant if the optimizer is looking for a traditional CZ gate.
                If True, the function will try to find a CZ gate using the 11<->20 transition. If False, the 11<->02 transition will be used instead.

        Set only ONE of these to True!:
            iSWAP (boolean) {Optional}: Optimize the iSWAP quantum gate.
            CZ (boolean) {Optional}: Optimize the CZ quantum gate.
    ---------------------------------------------------------
    OUTPUT:
            No output.
    ---------------------------------------------------------
    """

    #Set the longest allowed optimizing time for each algorithm.
    global maxRuntime
    maxRuntime = timeoutAlgorithm

    # Get the parameter bounds.
    wantTradCZ = (wantTradGate and CZ)
    wantTradiSWAP = (wantTradGate and iSWAP)
    parameterBounds = getParameterBounds(maxAllowedGateTime=maxAllowedGateTime, wantTradCZ=wantTradCZ, wantTradiSWAP=wantTradiSWAP, wantCZ_20=wantCZ_20)
    
    # Check the gate specifiers.
    if not (iSWAP or CZ):
        print("Specify which gate you would like to optimize the fidelity for!")
    elif (iSWAP and CZ):
        print("You can only optimize the fidelity for one gate at a time!")
    else:
        # Decide which cost function to use.
        if iSWAP:
            if energyLevels==2:
                cost = costiSWAP2
            elif energyLevels==3:
                cost = costiSWAP3
            else:
                cost = costiSWAP4
        else:
            if energyLevels==2:
                cost = costCZ2
            elif energyLevels==3:
                cost = costCZ3
            else:
                cost = costCZ4
        
        # Optimize the fidelity for the choosen gate. 
        findMinimum(cost, parameterBounds, runSHG=runSHG, runDA=runDA, runDE=runDE)


######################################################################################################################################################################
# Scrap functions

"""
def simulateEigenEnergies(x, numOfEnergyLevels=2, pointResolution=500):
    hamiltonian = hamiltonianModule.getHamiltonian(x, getEigenEnergies=True)
    
    eigenEnergies = [[] for _ in range(numOfEnergyLevels)] #[[]]*numOfEnergyLevels
    bounds = hamiltonianModule.getParameterBounds()
    start, end = bounds[0][0]+bounds[1][0], bounds[0][1]+bounds[1][1]
    listOfPhis = np.linspace(start, end, pointResolution)
    
    H = hamiltonian(0)
    eigenList = H.eigenstates()   
    # print(eigenList[1][1:4])
    
    fig, ax = plt.subplots()
    labels = ["|010><010|"]
    for phi in listOfPhis:
        H = hamiltonian(phi)
        eigenList = H.eigenstates()
        
        for i, energyLevel in enumerate(eigenEnergies):
            energyLevel.append(eigenList[0][i])
    
    for index, energy in enumerate(eigenEnergies, start=1):
        ax.plot(listOfPhis, energy)
        #labels.append("Eigenenergy " + str(index))
    ax.set_xlabel('Flux')
    ax.set_ylabel('Energy')
    ax.legend(labels)
    plt.show()
"""

######################################################################################################################################################################