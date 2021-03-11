from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *


def main():
    """H = hamiltonian(omega,g)
    state = tensor(basis(2,1),basis(2,0),basis(2,0))"""
    """H = testHamiltonian()
    state = basis(2,0)
    times = np.linspace(0,5,num=120)
    
    result = sesolve(H, state, times, [])
    animateStates(result, "temp")"""
    x = np.linspace(-10,10,500)
    y = testFun(x)
    plotVectors(x,y)
    res = optimizeGate(testFun, [(-7, 7)])
    print(res.x[0])


def testFun(x):
    return np.power(x,2) + 7*np.sin(np.pi*x)


main()