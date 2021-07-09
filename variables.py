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

# PI = np.pi
# HBAR  = 1.054571817e-34
# Q_e = 1.602176634e-19
# PHI0 = HBAR * np.pi / Q_e

# These are in Gigarad/s:
'''
# Our original testing values
omegas = np.array([2*np.pi*4.2, 2*np.pi*3.8, 2*np.pi*7.0]) # 4.02247024e+01]) # omega1 omega2 omegaTB0
gs = 2*np.pi * np.array([0.05, 0.05]) # g1 g2
alphas = -2*np.pi * np.array([0.150, 0.150, 0.150]) # alpha1 alpha2 alphaTB
'''
# Qubits 0 and 1 from 5-qubit design, April 2021
omegas = 2*np.pi * np.array([4.7381, 4.1762, 8.960]) # omega1 omega2 omegaTB0
gs = 2*np.pi * np.array([0.0545, 0.05552]) # g1 g2
alphas = -2*np.pi * np.array([0.221, 0.230, 0.08]) # alpha1 alpha2 alphaTB
'''
# Qubits 0 and 3 from 5-qubit design, April 2021
omegas = 2*np.pi * np.array([4.7381, 4.2936, 8.692]) # omega1 omega2 omegaTB0
gs = 2*np.pi * np.array([0.055, 0.05574]) # g1 g2
alphas = -2*np.pi * np.array([0.221, 0.226, 0.08]) # alpha1 alpha2 alphaTB
'''

omegaPhi_iSWAP = np.abs(omegas[0] - omegas[1])
omegaPhi_CZ_20 = np.abs(omegas[0] + alphas[0] - omegas[1])
omegaPhi_CZ_02 = np.abs(omegas[0] - (omegas[1] + alphas[1]))
phi_crossing = np.arccos((np.maximum(omegas[0], omegas[1])/omegas[2])**2)/np.pi

######################################################################################################################################################################