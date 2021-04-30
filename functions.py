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

# Last modified: 29 April 2021

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
        resDE = scipy.optimize.differential_evolution(costFunction, bounds, callback=callbackDE, workers=-1, updating='deferred')
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
# Simulation functions


def simulateHamiltonian(x0, sinStepHamiltonian=True, rotatingFrame=False, initialStateIndex=1, highestProjectionIndex=5):
    """
    This function simulates the population transfers between 
    different eigenstates of the 4-level hamiltonian from 
    McKay eq. 11. It also calculates the gate fidelity for both
    the iSWAP and CZ gate using the 4-level system.
    ---------------------------------------------------------
    INPUT:
            x0 (array(float)): The parameter set that the simulation is be performed for.
            sinStepHamiltonian (boolean) {Optional}: If True the amplitude of AC part of the flux signal that is applied to the tunnable bus will be sinusodially modulated.
            rotatingFrame (boolean) {Optional}: If True the states will be transformed into the rotating frame after the time evolution has been completed.
            initialStateIndex (int) {Optional}: This parameter decides which eigenstate that will be the initial state in the time evolution. If initialStateIndex=0 the eigenstate with the lowest associated energy will be the inital state.
            highestProjectionIndex (int) {Optional}: The eigenstates between, and including, the one with the lowest energy up to the (highestProjectionIndex)-lowest eigenstate will be projected onto.
    ---------------------------------------------------------
    OUTPUT:
            gateFidelity_iSWAP, gateFidelity_CZ (float, float): The gate fidelity for both the iSWAP and CZ gate.
    ---------------------------------------------------------
    """

    # Calculate the dimension of the tensor states and set the simulation time.
    N = 4
    D = N**3
    simulationTime = int(x0[-1]) + 25

    # Calculate the eigenstates and eigenenergies of the bare basis hamiltonian.
    hamiltonianBareBasis = getHamiltonian(x0,N=N,getBBHamiltonianComps=True)

    # Calculate the tunnable bus frequency when only the DC part of the flux is active.
    omegaTBDC = coeffomegaTB(x0[3],x0[0])

    # Calculate eigenstates and eigenenergies of the hamiltonian in the bare basis when the flux only has it's DC part.
    eigenStatesAndEnergies = getThetaEigenstates(x0, hamiltonianBareBasis[0], hamiltonianBareBasis[1], omegaTBDC)

    # Calculate the unitary for transforming the hamiltonian to the eigen basis.
    eigenBasisUnitary = getEBUnitary(x0, eigenStatesAndEnergies, N, D)

    # Get the hamiltonian that has a sinusodially modulated AC flux and also is in the eigen basis.
    hamiltonian = getHamiltonian(x0, N=N, eigEs=eigenStatesAndEnergies[0], U_e=eigenBasisUnitary, sinStepHamiltonian=sinStepHamiltonian)
    
    # Change the simulation settings and create the timestamps for where the evolution is to be evaluated. 
    options = solver.Options()
    options.nsteps = 10000
    timeStamps = np.linspace(0,simulationTime,simulationTime*5)
    
    # Create the initial state.
    initialState = Qobj(basis(D,initialStateIndex),dims=[[N,N,N],[1,1,1]])

    # Time evolve the initial state.
    result = sesolve(hamiltonian, initialState, timeStamps, [], options=options, args={'theta': x0[0], 'delta': x0[1], 'omegaphi': x0[2], 'omegatb0': x0[3], 'operationTime': x0[4], 'omegaTBTh': omegaTBDC})
    states = result.states

    # Transform into the rotating frame.
    if rotatingFrame:
        for i, t in enumerate(timeStamps):
            unitaryRotatingFrame = getRFUnitary(hamiltonian[0], t)
            states[i] = unitaryRotatingFrame*states[i]
    
    # Calculate the expectation values a projection operator.
    projectionOperators = []

    for index in range(highestProjectionIndex+1):
        operator = Qobj(basis(D,index),dims=[[N,N,N],[1,1,1]]) * Qobj(basis(D,index),dims=[[N,N,N],[1,1,1]]).dag()
        projectionOperators.append(operator)
    
    expectationValues = expect(projectionOperators, states)

    #Calculate gate fidelity for both iSWAP and CZ.
    gateFidelity_iSWAP = getGateFidelity(x0,N=4,wantiSWAP=True)
    gateFidelity_CZ = getGateFidelity(x0,N=4,wantCZ=True)

    # Print fidelity
    print(f'################################################\nGate fidelity for iSWAP: {gateFidelity_iSWAP}.\n\nGate fidelity for CZ: {gateFidelity_CZ}.\n################################################')

    # Plot the expectation values.
    plotExpect(timeStamps, expectationValues)

    return gateFidelity_iSWAP, gateFidelity_CZ


######################################################################################################################################################################
# Cost function definitions


def costiSWAP2(x):
    return -getGateFidelity(x,N=2,wantiSWAP=True)


def costiSWAP3(x):
    return -getGateFidelity(x,N=3,wantiSWAP=True)


def costiSWAP4(x):
    return -getGateFidelity(x,N=4,wantiSWAP=True)


def costCZ2(x):
    return -getGateFidelity(x,N=2,wantCZ=True)


def costCZ3(x):
    return -getGateFidelity(x,N=3,wantCZ=True)


def costCZ4(x):
    return -getGateFidelity(x,N=4,wantCZ=True)


######################################################################################################################################################################
# Optimize gate function


def optimizeGate(iSWAP=False,CZ=False,energyLevels=2, timeoutAlgorithm=40000, maxAllowedGateTime=240, runSHG=False, runDA=False, runDE=False):
    """
    The function tries to optimize the choosen gate for the
    choosen parameters, using the optimization algorithms 
    that the user decide would best fit.
    ---------------------------------------------------------
    INPUT:
            energyLevels (int) {Optional}: How many energy levels that should be accounted for in the simulations.
            timeoutAlgorithm (int) {Optional}: How many seconds each optimization algorithm can run, before a forced termination is initiated.
            maxAllowedGateTime (int) {Optional}: The longest gate time we will allow a solution to have.
            runSHG (boolean) {Optional}: If True the function will use the Simplicial Homology Global algorithm to optimize the gate.
            runDA (boolean) {Optional}: If True the function will use the Dual Anneling algorithm to optimize the gate.
            runDE (boolean) {Optional}: If True the function will use the Differential Evolution algorithm to optimize the gate.

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
    parameterBounds = getParameterBounds(maxAllowedGateTime=maxAllowedGateTime)
    
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


# Not sure if this one should be kept.
"""
def simulateEigenEnergies(hamiltonianModule, x, numOfEnergyLevels=4, pointResolution=500):
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