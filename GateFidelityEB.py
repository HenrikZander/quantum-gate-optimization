from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1 import *
import scipy.integrate as integrate

H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
H1 = (-1/2)*szTB

x0417DE = [0.38416891, 0.18051406, 3.84487738, 30.23212348, 101.4918881]
x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]
x0418DEb = [-0.422849290, 0.0583265568, 4.40324141, 50.0444593, 131.164936]

xUsed = x0418DE

nLevels = 2
D = nLevels**3

# Calculate eigenstates in the bare basis at Phi = Theta and use them to construct U_e
eigStsBB = getThetaEigenstates(xUsed, H0, H1)

# Define a list r of eigenstates in the eigenbasis
r = []
U_e = Qobj()
for i in range(D):
    r.append(Qobj(basis(D,i),dims=[[2,2,2],[1,1,1]]))
    U_e += r[i] * eigStsBB[1][i].dag()
# NB: r and U_e are ordered based on eigenenergies

# We are especially interested in |000>, |010>, |100> and |110>
eigIndices = [0, 1, 2, 4]

# Simulate evolution of eigenstates
# Calculate the eigenbasis hamiltonian
opTime = xUsed[4]
ts = np.linspace(0,opTime,500)
HBB = getSinStepHamiltonian(xUsed, operationTime=opTime)

HEB = [U_e * HBB[0] * U_e.dag(), [U_e * HBB[1][0] * U_e.dag(), HBB[1][1]]]

# Initialise a list c of the time-evolved eigenstates
c = r
# Calculate final states and store them in c
for i in range(len(c)):
    output = sesolve(HEB, c[i], ts)
    c[i] = output.states[-1]
# NB: c is ordered based on eigenenergies

# Calculate U_rf:
HBB_Th = H0 + xUsed[3]*np.sqrt(np.abs(np.cos(PI*xUsed[0]))) * H1
HEB_Th = U_e * HBB_Th * U_e.dag()

U_rf = (1j*HEB_Th*ts[-1]).expm()

# Transform c into the rotating frame
c_rf = U_rf * c

# Calculate M-matrix such that M_ij = <r_i|c_j>_rf:
# Initialize as a 4x4 zero nested list
M = [[0]*4]*4
# Assign values to elements M_ij
for i in range(4):
    for j in range(4):
        M[i][j] = r[eigIndices[i]].overlap(c_rf[eigIndices[j]])

# Calculate phases (iSWAP):
phi = np.angle(M[0][0])
theta1 = np.angle(M[0][1]) + PI/2 - phi
theta2 = np.angle(M[1][0]) + PI/2 - phi
# Calculate phases (CZ):
# phi = np.angle(M[0][0])
# theta1 = np.angle(M[0][0]) - phi
# theta2 = np.angle(M[1][1]) - phi

# Change M's type to matrix to simplify calculating fidelity
M = np.matrix(M)

# Ideal iSWAP gate matrix (with phases):
U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(-PI/2 + theta1 + phi)), 0], [0, np.exp(1j*(-PI/2 + theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + 2*phi))]])
# Ideal CZ gate matrix (with phases):
# U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, np.exp(1j*(theta1 + phi)), 0, 0], [0, 0, np.exp(1j*(theta2 + phi)), 0], [0, 0, 0, np.exp(1j*(PI + theta1 + theta2 + 2*phi))]])

# Calculate gate fidelity
N = 2
F_avg = (np.absolute(np.trace(M*U.H))**2 + np.trace(M.H*M)) / (N*(N+1))

print(F_avg)