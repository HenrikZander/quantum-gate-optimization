from qutip import *
import numpy as np
import matplotlib.pyplot as plt

from variables import *


def fidelityNew(state, targetState):
    return np.abs(np.sqrt(targetState.overlap(state)))
    

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