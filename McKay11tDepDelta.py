from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay11 import *

x1 = [1.0-0.56682599, 0.2387175, 2.52545336, 47.62454201] # Theta, delta, omegaPhi, omegaTB0
x2 = [2.0-1.42465745, 0.57298549, 2.52355871, 49.01101584] # Theta, delta, omegaPhi, omegaTB0
x3 = [-0.380725997,  0.0937065073, 0.0159535441,  39.0881746] # Theta, delta, omegaPhi, omegaTB0
x4 = [-0.33206001,  0.47870109, 1.26839156, 26.54997568] # Theta, delta, omegaPhi, omegaTB0

x0401DE = [-0.380725997,  0.0937065073,  0.0159535441,  39.0881746] # Theta, delta, omegaPhi, omegaTB0 # tOp = 112.5 / 295.5
x0401DA = [-0.35236044,  0.10932302,  0.84559989, 31.59465913] # Theta, delta, omegaPhi, omegaTB0 # tOp = 50.0 / 50.28

x0416DE = [0.2467117, 0.0433039, 0.46700076, 27.07391683]

xTest = [-0.25, 0.1, 2*PI*0.4, 2*PI*7.5] # funkar hyfsat -0.45, 0.0326, 0.18, 2*PI*7.5

xUsed = x0416DE

tMax = 76.16+5.0
ts = np.linspace(0,tMax,500)
psi0 = getInitialState()

pOps = getAllProjectionOperators()

eigpsi0 = getInitialEigenState(xUsed)
eigpOps = getAllThetaEigenProjectionOperators(xUsed)

opTime = tMax-5.0
tRise = 25.0
tWait = opTime - 2*tRise
def sinBox(t):
    return sinstep(t, 0, tRise) - sinstep(t, tWait + tRise, tWait + 2*tRise)
def Phi(x,t):
    return x[0] + sinBox(t) * x[1]*np.cos(x[2]*t)

Phis = Phi(xUsed,ts)

H = getSinStepHamiltonian(xUsed,opTime,tRise)

output = sesolve(H, eigpsi0, ts) # Verkar funka # , e_ops=eigpOps

plotEigenExpect(output, Phis, xUsed)