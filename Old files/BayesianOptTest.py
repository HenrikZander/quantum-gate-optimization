from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *

# Test av den Bayesianska optimeringen

def testFun(x):
    """A simple test funtion"""
    return 3*(1-x[0])**2*np.exp(-(x[0]**2) - (x[1]+1)**2) - 10*(x[0]/5 - x[0]**3 - x[1]**5)*np.exp(-x[0]**2-x[1]**2) - 1/3*np.exp(-(x[0]+1)**2 - x[1]**2)

res = gp_minimize(testFun, [(-4, 4),(-4,4)])
print(res.x[0], res.x[1], res.fun)