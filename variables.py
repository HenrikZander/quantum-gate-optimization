from qutip import *
import numpy as np

# We reduce all constans with units which means that hbar = 1 and Phi0 = 1 for example. 
# We also use nanoseconds as our time unit, which means our frequenzy unit becomes GHz and Gigarad/s.

PI = np.pi
# HBAR  = 1.054571817e-34
# Q_e = 1.602176634e-19
# PHI0 = HBAR * PI / Q_e

# These are in Gigarad/s:
omegas = 2*PI * np.array([4.2 ,3.8, 7.5, ])
gs = 2*PI * np.array([0.050, 0.050])
alpha = -2*PI * 0.150
omegaPhi = np.abs(omegas[0] - omegas[1])


Theta = -0.1 # DC flux on the tunable bus
delta = 0.21 # Amplitude of AC flux on the tunable bus


groundState = basis(2,0) #Ground state
excitedState = basis(2,1) #Excited state
gSt = groundState
eSt = excitedState

# Definition of basic operators 
sp = sigmap() #Raising operator for a two-level system
sm = sigmam() #Lowering operator for a two-level system
sz = sigmaz() #Pauli sigma-z matrice
sx = sigmax() #Pauli sigma-x matrice
sy = sigmay() #Pauli sigma-y matrice
I = qeye(2) #Indentity operator for a two-level system

# The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
# Upgraded states are defined as qubit one, qubit two and tunable bus in that order
sx1 = tensor(sx,I,I)
sx2 = tensor(I,sx,I)
sy1 = tensor(sy,I,I)
sy2 = tensor(I,sy,I)
sz1 = tensor(sz,I,I)
sz2 = tensor(I,sz,I)
szTB = tensor(I,I,sz)
sp1 = tensor(sp,I,I)
sp2 = tensor(I,sp,I)
spTB = tensor(I,I,sp)
sm1 = tensor(sm,I,I)
sm2 = tensor(I,sm,I)
smTB = tensor(I,I,sm)