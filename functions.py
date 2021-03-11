from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from skopt import gp_minimize

from variables import *


def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result

#def hamiltonian(omega, g):
#    return HBAR*( (-omega[0]/2)*sz1 + (-omega[1]/2)*sz2 + g[0]*(sp1*smTB + sm1*spTB) + g[1]*(sp2*smTB + sm2*spTB) + (-omega[2]/2)*szTB )


def testHamiltonian():
    return 2 * np.pi * 0.1 * sigmax()


def Phi(t):
    return Theta + delta*np.cos(omegaPhi*t)


def omegaTB(t, args):
    return omegas[2]*np.sqrt(np.abs(np.cos(PI*Phi(t))))


def optimizeGate(hamiltonian, bounds, maxFunEvals=500, numOfStartsPoints=50):
    res = gp_minimize(hamiltonian, bounds, acq_func="EI", n_calls=maxFunEvals, n_initial_points=numOfStartsPoints, random_state=1234)