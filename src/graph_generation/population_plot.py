######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: population_plot.py

# Author(s): Henrik Zander

# Date created: 9 January 2022

# Copyright 2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports.
from fpdf import FPDF
from pathlib import Path
import json

######################################################################################################################################################################

# ToDo: Slpit up this function according to the new file structure.
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

######################################################################################################################################################################