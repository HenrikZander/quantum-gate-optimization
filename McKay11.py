from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *


def getProjectionOperators():
    pSt1 = tensor(eSt3,gSt3,gSt3) # 100
    pOp1 = pSt1 * pSt1.dag()
    pSt2 = tensor(gSt3,eSt3,gSt3) # 010
    pOp2 = pSt2 * pSt2.dag()
    return [pOp2]

def getAllProjectionOperators():
    pSt1 = tensor(eSt3,gSt3,gSt3) # 100
    pOp1 = pSt1 * pSt1.dag()
    pSt2 = tensor(gSt3,eSt3,gSt3) # 010
    pOp2 = pSt2 * pSt2.dag()
    pStTB = tensor(gSt3,gSt3,eSt3) # 001
    pOpTB = pStTB * pStTB.dag()
    return [pOp1,pOp2,pOpTB]


def getHamiltonian(x):
    #The format of x is the following: x = [Theta, delta, omegaPhi, omegaTB0]
    H0 = omegas[0]*ad3_1*a3_1 - (alphas[0]/2.0)*(1-ad3_1*a3_1)*ad3_1*a3_1 + omegas[1]*ad3_2*a3_2 - (alphas[1]/2.0)*(1-ad3_2*a3_2)*ad3_2*a3_2 - (alphas[2]/2.0)*(1-ad3_TB*a3_TB)*ad3_TB*a3_TB + gs[0]*(ad3_1 + a3_1)*(ad3_TB + a3_TB) + gs[1]*(ad3_2 + a3_2)*(ad3_TB + a3_TB)
    H1 = ad3_TB*a3_TB
    def Phi(t):
        return x[0] + x[1]*np.cos(x[2]*t)
    def omegaTB(t, args):
        return x[3]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
    return [H0, [H1, omegaTB]]


def getInitialState():
    return tensor(excitedState3,groundState3,groundState3) #|100> ket

def getInitialGuess():
    return [Theta, delta, omegaPhi, omegas[2]]

def getParameterBounds():
    #Format of x: x = [Theta, delta, omegaPhi, omegaTB0]
    return [(-0.5,0.5),(0,0.25),(0,5),(20,60)]

# I did not make a version of McKay1's 'timeEvolutionH1()', because it seems to use omegaTB without defining it

def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result