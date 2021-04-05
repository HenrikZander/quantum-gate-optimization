from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay11 import *

ts = np.linspace(0,200,500)
psi0 = getInitialState()

pOps = getAllProjectionOperators()

x1 = [1.0-0.56682599, 0.2387175, 2.52545336, 47.62454201] # Theta, delta, omegaPhi, omegaTB0
x2 = [2.0-1.42465745, 0.57298549, 2.52355871, 49.01101584] # Theta, delta, omegaPhi, omegaTB0
x3 = [-0.380725997,  0.0937065073, 0.0159535441,  39.0881746] # Theta, delta, omegaPhi, omegaTB0
x4 = [-0.33206001,  0.47870109, -1.26839156, 26.54997568] # Theta, delta, omegaPhi, omegaTB0

x0401DE = [-0.380725997,  0.0937065073,  0.0159535441,  39.0881746] # Theta, delta, omegaPhi, omegaTB0
x0401DA = [-0.35236044,  0.10932302,  0.84559989, 31.59465913] # Theta, delta, omegaPhi, omegaTB0

xUsed = x0401DE

eigpOps = getAllEigenProjectionOperators(xUsed)
H = getHamiltonian(xUsed)

output = sesolve(H, psi0, ts, e_ops=eigpOps) # Verkar funka
plotExpect(output)