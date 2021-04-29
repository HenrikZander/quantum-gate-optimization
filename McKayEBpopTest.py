from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay11EB_4lvl import *

x0417DE = [0.38416891, 0.18051406, 3.84487738, 30.23212348, 101.4918881]
x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]
x0418DEb = [-0.422849290, 0.0583265568, 4.40324141, 50.0444593, 131.164936]

x0424iSWAP_3lvl = [0.241082420, 0.108883457, 2.50702612, 39.0692027, 218.127752]
x0425iSWAP_3lvl = [-0.431076591, 0.0275072184, 3.52975522, 36.7717108, 83.9876656]

x0426CZ_4lvl = [4.44287328e-01, 2.05934595e-02, 2.04012274e+00, 3.07227234e+01, 1.78093842e+02]
x0426CZ_3lvl = [-2.58309069e-01, 5.46369915e-02, 4.42420464e+00, 2.83414251e+01, 1.79648639e+02]

x = x0426CZ_4lvl

opTime = x[4]

tMax = opTime + 5.0
ts = np.linspace(0,tMax,500)

psi0 = getInitialState(5)
pOps = getCZProjectionOperators()

# Get both parts of the hamiltonian in the bare basis.
HBBComps = getHamiltonian(x,getBBHamiltonianComps=True)

# Get the number of considered energy levels for each qubit and the resulting dimension of the comnined tensor state.
# From here on, unless otherwise stated, every state is considered a tensor state.
n = getnLevels()
D = getD()

# Calculate omegaTB at Phi = Theta
omegaTBTh = coeffomegaTB(x[3],x[0])
# Calculate eigenstates and eigenenergies in the bare basis at Phi = Theta
eigStsBB = getThetaEigenstates(x, HBBComps[0], HBBComps[1], omegaTBTh)
# Get unitary for transformation into eigenbasis
U_e = getEBUnitary(x, eigStsBB, n, D)
# NB: r and U_e are ordered based on eigenenergies

# Calculate the eigenbasis hamiltonian
HEB = getHamiltonian(x, eigEs=eigStsBB[0], U_e=U_e, sinStepHamiltonian=True)
output = sesolve(HEB, psi0, ts, e_ops=pOps, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': x[3], 'operationTime': x[4], 'omegaTBTh': omegaTBTh})

plotExpect(output)