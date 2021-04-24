from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import multiprocessing

from variables import *
from functions import *
from plotting import *
import McKay11EB_3lvl
import time

# print(tensor(gSt,eSt,gSt))
# print(tensor(eSt,gSt,gSt))
# print(tensor(gSt,gSt,eSt))


def main():
    start = time.time()
    pointOfInterest = [0.38416891,   0.18051406,   3.84487738,  30.23212348, 101.4918881] #x = [Theta, delta, omegaPhi, omegaTB0, operationTime]
    currentHamiltonianModule = McKay11EB_3lvl
    testMultiprocessing(currentHamiltonianModule, pointOfInterest)
    # testNumbaSpeedup(currentHamiltonianModule)
    # testFindMinimum()
    # testGateOptimizer()
    # testGenerateCostFunction()
    # optimizeGate(currentHamiltonianModule, runDE=True, sinStepHamiltonian=True)
    # optimizeGateParallell(currentHamiltonianModule, runDE=True)
    # simulateHamiltonian(currentHamiltonianModule, pointOfInterest, simulationTime=200, sinStepHamiltonian=True)
    # simulateEigenEnergies(currentHamiltonianModule, pointOfInterest)
    print(f'Total running time: {time.time() - start} seconds.')


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

#######################################################################################
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


#######################################################################################
# Code for testing the possibility of multiprocessing the costfunction.


def testMultiprocessing(hamiltonianModule, x):
    numOfTimes = 10

    start = time.time()
    processes = []
    
    for _ in range(numOfTimes):
        p = multiprocessing.Process(target=costParallell, args=[x])
        p.start()
        processes.append(p)
    
    for process in processes:
        process.join()
    print(f'Total running time for cost evaluation in parallell: {time.time() - start} seconds.')
    
    start = time.time()
    for _ in range(numOfTimes):
        costParallell(x)
    print(f'Total running time for cost evaluation in a series: {time.time() - start} seconds.')

#######################################################################################


if __name__ == "__main__":
    main()
