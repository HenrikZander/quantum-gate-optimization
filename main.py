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

# Last modified: 29 April 2021

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
import time

######################################################################################################################################################################
# Interesting solutions in the format x = [Theta, delta, omegaPhi, omegaTB0, operationTime].

# Using H_rot = H_Th, including the interacting part of the hamiltonian:
x_0425_1_iSWAP = [0.241082420, 0.108883457, 2.50702612, 39.0692027, 218.127752]
x_0425_2_iSWAP = [-4.31076591e-01,  2.75072184e-02,  3.52975522e+00,  3.67717108e+01, 8.39876656e+01]
x_0426_CZ = [4.44287328e-01, 2.05934595e-02, 2.04012274e+00, 3.07227234e+01, 1.78093842e+02]

# Using H_rot = H_0^Th, excluding the interacting part of the hamiltonian:
x_0501_CZ_3lvl = [4.74835237e-01, 3.32890870e-02, 1.92795763e+00, 3.11288724e+01, 1.52910662e+02]

######################################################################################################################################################################
# The main function that auto-runs on compilation.


def main():
    start = time.time()
    # testMultiprocessing(currentHamiltonianModule, pointOfInterest)
    # testNumbaSpeedup(currentHamiltonianModule)
    # testFindMinimum()
    # testGateOptimizer()
    # testGenerateCostFunction()
    # optimizeGate(CZ=True, energyLevels=3, runDE=True)
    # optimizeGateParallel(currentHamiltonianModule, runDE=True)
    simulateHamiltonian(x_0501_CZ_3lvl, sinStepHamiltonian=True, rotatingFrame=True, initialStateIndex=5, N=4)
    # simulateEigenEnergies(currentHamiltonianModule, pointOfInterest)
    print(f'Total running time: {time.time() - start} seconds.')


######################################################################################################################################################################
# Functions used in testing.


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
    """Used to test the animatesStates-function in plotting.py.py"""
    H = testHamiltonian()
    state = basis(2, 0)
    times = np.linspace(0, 5, num=120)

    result = sesolve(H, state, times, [])
    animateStates(result, "temp")


######################################################################################################################################################################
# Code for testing Numba and how it can be used to speed up the code.
# x = [Theta, delta, omegaPhi, omegaTB0, operationTime]


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


######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.


if __name__ == "__main__":
    main()


######################################################################################################################################################################