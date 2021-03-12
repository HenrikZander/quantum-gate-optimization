from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *


def main():
    pass


def testFun(x):
    """A simple test funtion to be used by the testGateOptimizer-function"""
    return x[0]**2 + 7*np.sin(np.pi*x[0]) + 2


def testGateOptimizer():
    """Used to test the optimizeGate-function in functions.py"""
    x = np.linspace(-10,10,500)
    y=[]
    for tempx in x:
        y.append(testFun([tempx]))
    plotVectors(x,y)
    res = optimizeGate(testFun, [(-7, 7)])
    print(res.x[0], res.fun)


def testStateAnimations():
    """Used to test the animatesStates-function in plotting.py.py"""
    H = testHamiltonian()
    state = basis(2,0)
    times = np.linspace(0,5,num=120)
    
    result = sesolve(H, state, times, [])
    animateStates(result, "temp")

    
main()