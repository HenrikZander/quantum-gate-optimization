from qutip import *
import numpy as np
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

eigIndices = [0, 1, 2, 4]
eigpsi0 = getThetaEigenstates(xUsed, H0, H1)

r = eigpsi0[1]
c = r

opTime = xUsed[4]
tMax = opTime + 5.0
ts = np.linspace(0,opTime,500)
H = getSinStepHamiltonian(xUsed, operationTime=opTime)

for i in range(len(c)):
    output = sesolve(H, c[i], ts)
    c[i] = output.states[-1]

tempOp = 1j*H0*ts[-1]
U_rf = tempOp.expm()

c_rf = U_rf * c

M = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

for i in range(4):
    for j in range(4):
        M[i][j] = r[eigIndices[i]].overlap(c_rf[eigIndices[j]])

# Phases:
phi = np.angle(M[0][0])
theta1 = np.angle(M[0][1]) + PI/2 - phi
theta2 = np.angle(M[1][0]) + PI/2 - phi

U = [[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(-PI/2 + theta1 + phi)), 0], [0, np.exp(1j*(-PI/2 + theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + 2*phi))]]

Mmat = np.matrix(M)
Umat = np.matrix(U)

N = 2
F_avg = (np.absolute(np.trace(Mmat*Umat.H))**2 + np.trace(Mmat.H*Mmat)) / (N*(N+1))

print(F_avg)