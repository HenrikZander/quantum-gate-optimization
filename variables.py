######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: variables.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 27 February 2021

# Last modified: 22 May 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

from qutip import *
import numpy as np

######################################################################################################################################################################

# We reduce all constans with units which means that hbar = 1 and Phi0 = 1, for example. 
# We also use nanoseconds as our time unit, which means that our frequency unit becomes GHz and Gigarad/s.

PI = np.pi
# HBAR  = 1.054571817e-34
# Q_e = 1.602176634e-19
# PHI0 = HBAR * PI / Q_e

# These are in Gigarad/s:
omegas = np.array([2*PI*4.2 ,2*PI*3.8, 2*PI*7.0]) # 4.02247024e+01])
gs = 2*PI * np.array([0.05, 0.05])
alphas = -2*PI * np.array([0.150, 0.150, 0.150])
omegaPhi = np.abs(omegas[0] - omegas[1])


Theta = -0.1 # DC flux applied to the tunable bus
delta = 0.21 # Amplitude of AC flux applied to the tunable bus

# Ground and excited states for a two-level system
groundState = basis(2,0) # Ground state
excitedState = basis(2,1) # Excited state
gSt = groundState
eSt = excitedState

# States for a three-level system
groundState3 = basis(3,0) # Ground state
excitedState3 = basis(3,1) # First excited state
upperState3 = basis(3,2) # Second excited state
gSt3 = groundState3
eSt3 = excitedState3
uSt3 = upperState3

# Definition of basic operators for a two-level system
sp = sigmap() # Raising operator for a two-level system
sm = sigmam() # Lowering operator for a two-level system
sx = sigmax() # Pauli sigma-x operator
sy = sigmay() # Pauli sigma-y operator
sz = sigmaz() # Pauli sigma-z operator
I = qeye(2) # Identity operator for a two-level system

# The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
# Upgraded states are defined as qubit one, qubit two and tunable bus in that order
sx1 = tensor(sx,I,I)
sx2 = tensor(I,sx,I)
sxTB = tensor(I,I,sx)

sy1 = tensor(sy,I,I)
sy2 = tensor(I,sy,I)
syTB = tensor(I,I,sy)

sz1 = tensor(sz,I,I)
sz2 = tensor(I,sz,I)
szTB = tensor(I,I,sz)

sp1 = tensor(sp,I,I)
sp2 = tensor(I,sp,I)
spTB = tensor(I,I,sp)

sm1 = tensor(sm,I,I)
sm2 = tensor(I,sm,I)
smTB = tensor(I,I,sm)

# Definition of basic operators for a three-level system
a3 = destroy(3) # Annihilation operator for a three-level system
ad3 = create(3) # Creation operator for a three-level system
I3 = qeye(3) # Identity operator for a three-level system

# The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
# Upgraded states are defined as qubit one, qubit two and tunable bus in that order

a3_1 = tensor(a3,I3,I3)
a3_2 = tensor(I3,a3,I3)
a3_TB = tensor(I3,I3,a3)

ad3_1 = tensor(ad3,I3,I3)
ad3_2 = tensor(I3,ad3,I3)
ad3_TB = tensor(I3,I3,ad3)

# Definition of basic operators for a four-level system
a4 = destroy(4) # Annihilation operator for a four-level system
ad4 = create(4) # Creation operator for a four-level system
I4 = qeye(4) # Identity operator for a four-level system

# The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
# Upgraded states are defined as qubit one, qubit two and tunable bus in that order

a4_1 = tensor(a4,I4,I4)
a4_2 = tensor(I4,a4,I4)
a4_TB = tensor(I4,I4,a4)

ad4_1 = tensor(ad4,I4,I4)
ad4_2 = tensor(I4,ad4,I4)
ad4_TB = tensor(I4,I4,ad4)

######################################################################################################################################################################