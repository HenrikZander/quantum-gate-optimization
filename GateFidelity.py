from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1 import *
import scipy.integrate as integrate

#from variables import *
#from functions import *
#from plotting import *
# comment


# Define the initial states
r1=tensor(gSt,gSt) # |00>
r2=tensor(gSt,eSt) # |01>
r3=tensor(eSt,gSt) # |10>
r4=tensor(eSt,eSt) # |11>
r=[r1, r2, r3, r4]
pOp00 = r1 * r1.dag()
pOp01 = r2 * r2.dag()
pOp10 = r3 * r3.dag()
pOp11 = r4 * r4.dag()

ts=np.linspace(0,500,500):

M = np.zeros((4,4))

for t in ts:

    # Define the optimized hamiltioninans parameters
    x = [1.0-0.56682599, 0.2387175, 2.52545336, 47.62454201] # Theta, delta, omegaPhi, omegaTB0

    def Phi(t):
        return x[0] + x[1]*np.cos(x[2]*t)
    def omegaTB(t, args):
        return x[3]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
    
    H_offdiag = gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H_diag =(-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + (-1/2)*omegaTB*szTB
    I=np.matrix(integrate.quad(lambda t: H_diag, 0, t))
    U_I=1j*exp(*I)
    H=U_I*(H_diag+H_offdiag)*U_I.dag()-1j*U_I*U_I.dag() # den sista termen med U_I.dag() ska även differentiera, men osäker på tisberoendet av H_diag

    # Calculate c_j in to the roatingframe 
    c1=sesolve(H, r1, t)
    c2=sesolve(H, r2, t)
    c3=sesolve(H, r3, t)
    c4=sesolve(H, r4, t)
    c=[c1, c2, c3, c4]

       

    # Calculating the 4x4 M matrix
    for i in range(4):

        for j in range(4):

            M[i+1][j+1]=r[i].dag*c[j]


    
    if ISWAP:
        gate=iswap(N=None,targets=[0,1])
        U=[[np.exp((1j*phi)), 0, 0, 0],[0, 0, np.exp((1j*(-(np.pi/2)+theta1+phi))), 0],[0, np.exp((1j*(-(np.pi/2)+theta2+phi))), 0, 0],[0, 0, 0, np.exp((theta1+theta2+2*phi))]]
        phi=np.angle(M[1][1])
        theta1=np.angle(M[1][2])+(np.pi/2)-phi
        theta2=np.angle(M[2][1])+(np.pi/2)-phi

    if CZ:
        gate=cphase(np.pi, N=2, control=0, target=1)
        # Måste sätta in U_CZ

    avg_gt_fd=(abs((M*U.dag()).tr())^2)/(ISWAP.shape[1]*(ISWAP.shape[1]+1))