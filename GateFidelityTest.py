from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1 import *
# import scipy.integrate as integrate

H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
H1 = (-1/2)*szTB

x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]

xUsed = x0418DE

# Calculate eigenstates at Phi = Theta and store them in a list r
eigpsi0 = getThetaEigenstates(xUsed, H0, H1)
r = eigpsi0[1]

# We are especially interested in |000>, |010>, |100> and |110>
eigIndices = [0, 1, 2, 4]

# Simulate evolution of eigenstates
opTime = xUsed[4]
ts = np.linspace(0,opTime,500)
H = getSinStepHamiltonian(xUsed, operationTime=opTime)

# Initialise a list c of the time-evolved eigenstates
c = r
# Calculate final states and store them in c
for i in range(len(c)):
    output = sesolve(H, c[i], ts)
    c[i] = output.states[-1]

# Since H0 is time-independent (right?), we can calculate rotating frame unitary like
U_rf = (1j*H0*ts[-1]).expm()

# Transform c into the rotating frame
c_rf = U_rf * c

# Calculate M-matrix such that M_ij = <r_i|c_j>_rf:
# Initialize as a 4x4 zero list
M = [[0]*4]*4
# Assign values to elements M_ij
for i in range(4):
    for j in range(4):
        M[i][j] = r[eigIndices[i]].overlap(c_rf[eigIndices[j]])

# Calculate phases:
phi = np.angle(M[0][0])
theta1 = np.angle(M[0][1]) + PI/2 - phi
theta2 = np.angle(M[1][0]) + PI/2 - phi

# Change type of M to matrix to simplify calculating fidelity
M = np.matrix(M)

# Ideal iSWAP gate matrix (with phases):
U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(-PI/2 + theta1 + phi)), 0], [0, np.exp(1j*(-PI/2 + theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + 2*phi))]])

# Calculate gate fidelity
N = 2
F_avg = (np.absolute(np.trace(M*U.H))**2 + np.trace(M.H*M)) / (N*(N+1))

print(F_avg)