from qutip import *
import numpy as np
import matplotlib.pyplot as plt

#from variables import *
#from functions import *
#from plotting import *
# comment

ISWAP=iswap(N=None,targets=[0,1])

avg_gate_fidelity=(abs((M*U.dag()).tr())^2)/(ISWAP.shape[1]*(ISWAP.shape[1]+1))

phi=np.angle(M11)
theta1=np.angle(M12)+(np.pi/2)-phi
theta2=np.angle(M21)+(np.pi/2)-phi

r1=tensor(gSt,gSt) # |00>
r2=tensor(gSt,eSt) # |01>
r3=tensor(eSt,gSt) # |10>
r4=tensor(eSt,eSt) # |11>

U_iswap=[[np.exp((1j*phi)), 0, 0, 0],[0, 0, np.exp((1j*(-(np.pi/2)+theta1+phi))), 0],[0, np.exp((1j*(-(np.pi/2)+theta2+phi))), 0, 0],[0, 0, 0, np.exp((theta1+theta2+2*phi))]]

