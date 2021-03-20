from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *
from functions import *
from plotting import *

# Test av den Bayesianska optimeringen

def testFun(x):
    """A simple test funtion"""
    return 7*x[0]**2 + 10*np.sin(5*np.pi*x[0]) + 2

res = gp_minimize(testFun, [(-1.5,0.5)])
print(res.x[0], res.fun)