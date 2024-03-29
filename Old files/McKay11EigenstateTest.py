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

x0401DE = [-0.380725997,  0.0937065073,  0.0159535441,  39.0881746] # Theta, delta, omegaPhi, omegaTB0 # tOp = 112.5
x0401DA = [-0.35236044,  0.10932302,  0.84559989, 31.59465913] # Theta, delta, omegaPhi, omegaTB0 # tOp = 50

H0 = omegas[0]*ad3_1*a3_1 - (alphas[0]/2.0)*(1-ad3_1*a3_1)*ad3_1*a3_1 + omegas[1]*ad3_2*a3_2 - (alphas[1]/2.0)*(1-ad3_2*a3_2)*ad3_2*a3_2 - (alphas[2]/2.0)*(1-ad3_TB*a3_TB)*ad3_TB*a3_TB + gs[0]*(ad3_1 + a3_1)*(ad3_TB + a3_TB) + gs[1]*(ad3_2 + a3_2)*(ad3_TB + a3_TB)
H1 = ad3_TB*a3_TB

xUsed = x0401DA

# eigSts = getZeroEigenstates(xUsed, H_const=H0, H_omegaTB=H1)
# print(eigSts[1])

Phis = [0, 0.05]
pOps = []
eigenStateIndices = [2, 1, 3]
for k in range(len(eigenStateIndices)):
    pOps.append([])
    for Phi in Phis:
        pOps[k].append(getEigenProjectionOperator(xUsed,Phi,eigenStateIndices[k]))



print(pOps[0])
