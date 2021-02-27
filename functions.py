from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *


def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)

    result = 0
    for n in range(0, N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n

    result *= x ** (N + 1)

    return result

def hamiltonian(omega, g):
    return HBAR*( (-omega[0]/2)*sz1 + (-omega[1]/2)*sz2 + g[0]*(sp1*smTB + sm1*spTB) + g[1]*(sp2*smTB + sm2*spTB) + (-omega[2]/2)*szTB )


def testHamiltonian():
    return 2 * np.pi * 0.1 * sigmax()


def omegaPhi(omega):
    pass


def omegaTilde(omega):
    pass


def delta(time):
    pass