from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *


def main():
    testGateOptimizer()
    return 0


def testFun(x):
    """A simple test funtion to be used by the testGateOptimizer-function"""
    return 3*(1-x[0])**2*np.exp(-(x[0]**2) - (x[1]+1)**2) - 10*(x[0]/5 - x[0]**3 - x[1]**5)*np.exp(-x[0]**2-x[1]**2) - 1/3*np.exp(-(x[0]+1)**2 - x[1]**2) 


def testGateOptimizer():
    """Used to test the optimizeGate-function in functions.py"""
    res = findMinimum(testFun, [(-4, 4),(-4,4)], runBayesianWithBH=True)
    #print(res[0])


def testStateAnimations():
    """Used to test the animatesStates-function in plotting.py.py"""
    H = testHamiltonian()
    state = basis(2,0)
    times = np.linspace(0,5,num=120)
    
    result = sesolve(H, state, times, [])
    animateStates(result, "temp")

    
main()