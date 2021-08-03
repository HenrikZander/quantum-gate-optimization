######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: main.py

# Author(s): Henrik Zander

# Date created: 27 February 2021

# Copyright 2021, Henrik Zander, All rights reserved.

######################################################################################################################################################################

from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import multiprocessing

from variables import *
from functions import *
from plotting import *
import model
from datetime import datetime
import time as t
from dataManager import *

######################################################################################################################################################################
# Interesting solutions in the format x = [Theta, delta, omegaPhi, operationTime].

# Poor solutions:

# x_0425_2_iSWAP = [-4.31076591e-01,  2.75072184e-02,  3.52975522e+00,  3.67717108e+01, 8.39876656e+01]
# x_0502_iSWAP_1_3lvl = [-3.28050858e-01,  7.83176577e-02,  2.02152681e+00,  3.06207420e+01, 1.03221706e+02]
# x_0502_iSWAP_3_3lvl = [3.96689737e-01, 6.81485794e-02, 4.23032165e+00, 2.82614635e+01, 1.27222763e+02]
# x_0503_iSWAP_1_4lvl = [-3.98464392e-01,  4.89012371e-02,  3.55650561e+00,  3.64347107e+01, 9.98881018e+01] # Fideliteten stiger snabbt efter grinden stängts av. Men dålig fidelitet.
# x_0509_iSWAP_2_3lvl = [-3.19991864e-01,  9.38707220e-02,  2.52146782e+00,  3.87216564e+01, 9.45620521e+01] # Inte tillräckligt hög fidelitet.

# x_0426_CZ = [4.44287328e-01, 2.05934595e-02, 2.04012274e+00, 3.07227234e+01, 1.78093842e+02] # Oscillerande fidelitet
# x_0501_CZ_3lvl = [4.74835237e-01, 3.32890870e-02, 1.92795763e+00, 3.11288724e+01, 1.52910662e+02]

# Good solutions:

x_0425_1_iSWAP = [0.241082420, 0.108883457, 2.50702612, 39.0692027, 218.127752]                                                                                # F = [(0.9987825701889115+0j), (0.9987615809690565+0j), (0.9987404349878605+0j), (0.9987187305959936+0j)]
x_0501_iSWAP_3lvl = [2.87022047e-01, 1.07442615e-01, 2.50566092e+00, 4.51975847e+01, 1.72672744e+02]                                                           # F = [(0.9970235907731843+0j), (0.9969987243890179+0j), (0.996973740723554+0j), (0.9969486009884246+0j)]
x_0504_iSWAP_1_3lvl = [4.56199639e-01, 3.21813167e-02, 3.68027336e+00, 4.63027865e+01, 5.16555865e+01] # Väldigt hög fidelitet och väldigt snabb.              # F = [(0.9996386695754572+0j), (0.998520052448225+0j), (0.9958511334297339+0j), (0.9916456385963169+0j)]
x_0502_iSWAP_2_3lvl = [-3.23785448e-01,  3.74318501e-02,  2.45144003e+00, 9.51244239e+01] # ometaTB0 = 4.02247024e+01 # Fidelitet stiger efter att griden stängts av.     # F = [(0.9945143458523513+0j), (0.9946846702621196+0j), (0.9948477877462383+0j), (0.9950036555493578+0j)]
x_0508_iSWAP_1_3lvl = [-3.77602967e-01,  7.79255423e-02,  3.55268857e+00,  2.83153193e+01, 9.99260963e+01] # Hög fidelitet okej operationstid.                 # F = [(0.9996315297932361+0j), (0.9995749397182396+0j), (0.9977996490730219+0j), (0.9943159415431349+0j)]
x_0509_iSWAP_1_3lvl = [4.67882816e-01, 4.10434975e-02, 3.57317326e+00, 4.45333793e+01, 6.83961701e+01]                                                         # F = [(0.9963800344710967+0j), (0.9948909833881396+0j), (0.9930313677770233+0j), (0.9908035658030836+0j)]

x_0508_CZ_1_4lvl = [-3.81665839e-01,  5.24394496e-02,  3.96158842e+00,  2.93721722e+01, 1.57298060e+02] # Väldigt hög fidelitet                                # F = [(0.9980325063437839+0j), (0.9998233507400702+0j), (0.9995032108558896+0j), (0.9970728189390012+0j)]
x_0509_CZ_1_4lvl = [4.20660096e-01, 6.85315342e-02, 4.86571561e-02, 3.30879395e+01, 9.29363646e+01] # Ännu högre fidelitet och mer stabil efteråt.             # F = [(0.9997256108907016+0j), (0.9999223226253042+0j), (0.9990775163025432+0j), (0.9971938931736171+0j)]

# New solution with a correct identification of the eigenstates:
x_0522_iSWAP_1_3lvl = [-0.21598281,  0.18286676,  2.49402702, 98.72330723]
x_0523_iSWAP_1_3lvl_crashed = [-4.26842714e-01, 2.92423430e-02, 2.48443602e+00, 9.99094553e+01]

x_0522_CZ_1_4lvl = [-3.42740346e-01,  5.14288612e-02,  2.97036730e+00,  1.13806249e+02]
x_0523_CZ_1_4lvl = [4.59175954e-01, 5.89394560e-03, 3.51909227e-01, 6.33353067e+01]
x_0525_CZ_1_4lvl = [-2.53065028e-01,  1.49440908e-01,  2.12682870e+00,  1.75050178e+02]
x_0525_CZ_2_4lvl = [-3.42944957e-01,  4.54371850e-02,  3.03927755e+00,  1.75860451e+02]
x_0525_CZ_3_4lvl = [-3.33260202e-01,  9.12536869e-02,  3.03332840e+00,  1.53983390e+02]
x_0525_CZ_4_4lvl = [3.31034106e-01, 7.56628827e-02, 2.03924480e+00, 1.16525277e+02]

# Presented solutions in the report:
x_0502_iSWAP_2_3lvl = [-3.23785448e-01,  3.74318501e-02,  2.45144003e+00, 9.51244239e+01] # omegaTB0 = 4.02247024e+01 # Fidelitet stiger efter att griden stängts av.     # F = [(0.9945143458523513+0j), (0.9946846702621196+0j), (0.9948477877462383+0j), (0.9950036555493578+0j)]
x_0522_CZ_1_4lvl = [-3.42740346e-01,  5.14288612e-02,  2.97036730e+00,  1.13806249e+02]

# Only solutions below this point are to be added to solutions.json:

# Solutions for qubit pair 01 at Chalmers:
x_210706_iSWAP_4lvl = [-0.24300583, 0.18855054, 3.52505319, 139.68066173]
x_210707_CZ_1_4lvl = [-3.66144215e-01, 7.76542110e-02, 5.31979461e+00, 1.35483199e+02]
x_210707_CZ_2_4lvl = [-3.97197613e-01, 4.50017315e-02, 5.04744811e+00, 1.20619159e+02]
x_210708_CZ_1_4lvl = [-3.87809473e-01, 4.70332021e-02, 2.01877517e+00, 1.37128800e+02]
x_210710_iSWAP_1_4lvl = [3.91437575e-01, 6.30717089e-02, 3.54292316e+00, 9.99427537e+01]
x_210713_iSWAP_1_4lvl = [4.03166181e-01, 2.16061627e-02, 3.50015579e+00, 7.98055976e+01] # Bad: fidelity only 0.98
x_210713_iSWAP_2_4lvl = [-0.2452259, 0.19856197, 3.52316351, 62.77634052] # Terrible: fidelity only 0.85
x_210713_iSWAP_3_4lvl = [3.77667126e-01, 5.17617928e-02, 3.51608253e+00, 8.83255578e+01] # Bad: fidelity only 0.987
x_210715_CZ_1_4lvl = [4.28265258e-01, 2.28050495e-02, 4.28583819e+00, 1.15226156e+02]
x_210716_CZ_2_4lvl = [4.03140410e-01, 2.95892756e-02, 4.14970423e+00, 1.06035386e+02]
# Solution to use in simulations:
solName = "210801_SWAP_2"
#xUsed = x_210713_iSWAP_1_4lvl
xName = 'x_' + solName + '_4lvl' # Preferred if solution is listed in solutions.json
#dateStr = "2021-07-13 15:12:18.004557"
solPath = './solutions_qubitPair01/' + solName + '.json' # './circuit files/qubitPair01.json'
solDict = getFromjson(solPath)

######################################################################################################################################################################
# The main function that auto-runs on compilation.


def main():
    start = t.time()

    # for index in range(64):
    #     findEigenIndex(xUsed, eigenStateIndex=index, N=4, printResult=True)
    # optimizeGate(CZ=True, energyLevels=4, maxAllowedGateTime=120, runDE=True, wantTradGate=False, wantCZ_20=False)
    # optimizeGate(CZ=True, energyLevels=4, maxAllowedGateTime=240, runDE=True)
    # optimizeGate(CZ=True, energyLevels=4, maxAllowedGateTime=200, runDE=True)
    # optimizeGate(CZ=True, energyLevels=4, maxAllowedGateTime=170, runDE=True)
    # optimizeGate(CZ=True, energyLevels=4, maxAllowedGateTime=140, runDE=True)
    # simulateHamiltonian(xName=xName, sinStepHamiltonian=True, rotatingFrame=True, initialStateIndex=3, N=4, highestProjectionIndex=12, circuitData=getFromjson('./circuit files/qubitPair01.json'))
    # simulatePopTransfer(solutionPath=solPath, initialStateIndex=1, highestProjectionIndex=12)
    # plotFidelityOld(xName=xName, iSWAP=True, useSavedPlot=False, saveTojson=False, circuitData=solDict)
    plotFidelity(solutionPath=solPath, useSavedPlot=False, saveToFile=True)
    # deltaPulsePlot()
    # testPlotStates()
    # testEigenstateOrder()
    # plotEigenenergiesOld(xName=xName, N=4, simPoints=500, numOfEnergyLevels=None, useSavedPlot=True, saveTojson=False)
    # plotEigenenergies(solutionPath=solPath, eigenenergiesPath='eigenenergies.json', N=4, simPoints=200, numOfEnergyLevels=None, useSavedPlot=True, saveToFile=False)
    # test()
    # getRobustnessPlot(solutionPath=solPath, useSavedPlot=False, saveToFile=True, checkOpTime=True, nPointsList=[33], maxDevs=[0.005, 0.005, 0.002, 5])
    # print(getFromjson('solutions.json'))
    # print(getSolutionNameList())
    # addNewSolution(xUsed, gateType="iSWAP", N=4, creationTime=datetime.fromisoformat(dateStr))
    # solsDict = getFromjson('solutions.json')
    # print(solsDict[xName]['fidelitiesAtTimes'][-45])

    print(f'Total running time: {t.time() - start} seconds.')


######################################################################################################################################################################
# Functions used in testing.

'''
circuitDict = getFromjson('./circuit files/qubitPair01.json')
for key in circuitDict:
    circuitDict[key] = (2*np.pi*np.array(circuitDict[key])).tolist()
dumpTojson(circuitDict, './circuit files/qubitPair01.json')
'''

def testEigenstateOrder():
    plotEigenenergies(xUsed, N=4, simPoints=200, numOfEnergyLevels=None)
    # HBBComps = model.getHamiltonian(xUsed, N=N, getBBHamiltonianComps=True)
    # omegaTBTh = model.coeffomegaTB(omegas[2], xUsed[0])
    # eigStsBB = model.getThetaEigenstates(xUsed, HBBComps[0]+HBBComps[1], HBBComps[2], omegaTBTh)
    # print(model.eigenstateOrder(eigStsBB[0], eigStsBB[1], N))


def testSpeedOfModifiedGateFidelity():
    temp = 0
    numTimes = 10
    x = x_0503_iSWAP_1_4lvl

    # Unmodified gate fidelity
    start = time.time()
    for i in range(numTimes):
        temp = model.getGateFidelity(x, N=4, iSWAP=True)
    print(f'Total running time for UNMODIFIED gate fidelity: {time.time() - start} seconds.')
    print(f'F_avg = {temp}.\n')
    
    # Modified gate fidelity
    start = time.time()
    for i in range(numTimes):
        temp = model_modified.getGateFidelity(x, N=4, iSWAP=True)
    print(f'Total running time for MODIFIED gate fidelity: {time.time() - start} seconds.')
    print(f'F_avg = {temp}.\n')


def testEvaluateResult():
    a = [[[1, 1], [2, 2], [3, 3]], [1, 2, 3]]
    print(a)
    print(evaluateResult([10, 10], 4, a))


def testGenerateCostFunction():
    #x0 = [Theta, delta, omegaPhi, omegas[2]]
    # x0 = [0.27704506, -0.20963235, 2.30482495, 46.81909043]
    x0 = [-1.11740019,  0.90085377,  2.52085547, 49.17815508]

    cost = generateCostFunction(McKay1.getHamiltonian, McKay1.getEigenStates)
    print(cost(x0))


def testFun(x):
    """A simple test funtion to be used by the testGateOptimizer-function"""
    return 3*(1-x[0])**2*np.exp(-(x[0]**2) - (x[1]+1)**2) - 10*(x[0]/5 - x[0]**3 - x[1]**5)*np.exp(-x[0]**2-x[1]**2) - 1/3*np.exp(-(x[0]+1)**2 - x[1]**2)


def testFindMinimum():
    """Used to test the optimizeGate-function in functions.py"""
    res = findMinimum(testFun, [(-4, 4), (-4, 4)])
    # print(res[0])


def testStateAnimations():
    """Used to test the animatesStates-function in plotting.py"""
    H = sx
    state = basis(2, 0)
    times = np.linspace(0, 5, 120)

    result = sesolve(H, state, times, [])
    animateStates(result, "temp")


def testPlotStates():
    """Used to test the plotStates-function in plotting.py"""
    H = sx
    state = basis(2, 0)
    times = np.linspace(0, np.pi/2, 20)

    result = sesolve(H, state, times, [])
    plotStates(result)


######################################################################################################################################################################
# Code for testing Numba and how it can be used to speed up the code.
# x = [Theta, delta, omegaPhi, operationTime]


def testNumbaSpeedup(hamiltonianModule):
    fun = hamiltonianModule.omegaTB
    temp = 0
    start = time.time()
    for i in range(10000000):
        temp = fun(1,{'theta': 0.2467117, 'delta': 0.0433039, 'omegaphi':0.46700076, 'omegatb0': 27.07391683, 'operationTime': 150})
    print(f'Total running time for none-sinstep omegaTB: {time.time() - start} seconds.')
    print("Temp=", temp)
    
    fun = hamiltonianModule.omegaTBSinStep
    temp = 0
    start = time.time()
    for i in range(10000000):
        temp = fun(1,{'theta': 0.2467117, 'delta': 0.0433039, 'omegaphi':0.46700076, 'omegatb0': 27.07391683, 'operationTime': 150})
    print(f'Total running time for sinstep omegaTB: {time.time() - start} seconds.')
    print("Temp=", temp)


######################################################################################################################################################################
# Code for testing the possibility of multiprocessing the costfunction.


def testMultiprocessing(hamiltonianModule, x):
    numOfTimes = 10

    start = time.time()
    processes = []
    
    for _ in range(numOfTimes):
        p = multiprocessing.Process(target=costParallel, args=[x])
        p.start()
        processes.append(p)
    
    for process in processes:
        process.join()
    print(f'Total running time for cost evaluation in parallel: {time.time() - start} seconds.')
    
    start = time.time()
    for _ in range(numOfTimes):
        costParallel(x)
    print(f'Total running time for cost evaluation in a series: {time.time() - start} seconds.')


def testPrint():
    global cancelOptimization
    print(f'Test {cancelOptimization}')

######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.


if __name__ == "__main__":
    main()


######################################################################################################################################################################
