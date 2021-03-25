from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *
import McKay1
import time


def main():
    start = time.time()
    #testFindMinimum()
    #testGateOptimizer()
    testGenerateCostFunction()
    print(f'Total running time: {time.time() - start} seconds.')

def testGateOptimizer():
    hamiltonian = McKay1.getHamiltonian1
    #x = [Theta, delta, omegaPhi, omegaTB0]
    parameterBounds = [(-2,2),(0,2),(20,60),(20,60)]
    initialGuess = [Theta, delta, omegaPhi, omegas[2]]
    initialState = tensor(excitedState,groundState,groundState)
    
    projectionOperators = McKay1.getProjectionOperators()
    optimizeGate(hamiltonian, parameterBounds, projectionOperators, initialGuess=initialGuess, initialState=initialState, runBayesian=True)
    

def testGenerateCostFunction():
    #x0 = [Theta, delta, omegaPhi, omegas[2]]
    x0 = [0.27704506, -0.20963235, 2.30482495, 46.81909043]
    initialState = tensor(excitedState,groundState,groundState)
    
    pSt1 = tensor(eSt,gSt,gSt) # 100
    pOp1 = pSt1 * pSt1.dag()
    pSt2 = tensor(gSt,eSt,gSt) # 010
    pOp2 = pSt2 * pSt2.dag()
    pStTB = tensor(gSt,gSt,eSt) # 001
    pOpTB = pStTB * pStTB.dag()
    projectionOperators = McKay1.getProjectionOperators()
    
    hamiltonian = McKay1.getHamiltonian1
    cost = generateCostFunction(hamiltonian, projectionOperators, initialState=initialState)
    print(cost(x0))
    
    
def testFun(x):
    """A simple test funtion to be used by the testGateOptimizer-function"""
    return 3*(1-x[0])**2*np.exp(-(x[0]**2) - (x[1]+1)**2) - 10*(x[0]/5 - x[0]**3 - x[1]**5)*np.exp(-x[0]**2-x[1]**2) - 1/3*np.exp(-(x[0]+1)**2 - x[1]**2) 


def testFindMinimum():
    """Used to test the optimizeGate-function in functions.py"""
    res = findMinimum(testFun, [(-4, 4),(-4,4)], runBayesian=True)
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