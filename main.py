from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *
import McKay1
import time

# print(tensor(gSt,eSt,gSt))
# print(tensor(eSt,gSt,gSt))
# print(tensor(gSt,gSt,eSt))


def main():
    start = time.time()
    pointOfInterest = [ 0.2467117 ,  0.0433039 ,  0.46700076, 27.07391683]
    currentHamiltonianModule = McKay1
    # testFindMinimum()
    # testGateOptimizer()
    # testGenerateCostFunction()
    # optimizeGate(currentHamiltonianModule, runDE=True)
    # simulateHamiltonian(currentHamiltonianModule, pointOfInterest, simulationTime=100)
    # simulateEigenEnergies(currentHamiltonianModule, pointOfInterest)
    print(f'Total running time: {time.time() - start} seconds.')


def testEvaluateResult():
    a = [[[1,1],[2,2],[3,3]],[1,2,3]]
    print(a)
    print(evaluateResult([10,10],4,a))


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
    res = findMinimum(testFun, [(-4, 4),(-4,4)])
    #print(res[0])


def testStateAnimations():
    """Used to test the animatesStates-function in plotting.py.py"""
    H = testHamiltonian()
    state = basis(2,0)
    times = np.linspace(0,5,num=120)
    
    result = sesolve(H, state, times, [])
    animateStates(result, "temp")


if __name__ == "__main__":
    main()