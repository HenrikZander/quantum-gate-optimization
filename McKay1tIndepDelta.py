from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *

ts = np.linspace(0,2e-6,1000)
psi0 = tensor(eSt,gSt,gSt)

pSt1 = tensor(eSt,gSt,gSt) # 100
pOp1 = pSt1 * pSt1.dag()
pSt2 = tensor(gSt,eSt,gSt) # 010
pOp2 = pSt2 * pSt2.dag()
pStTB = tensor(gSt,gSt,eSt) # 001
pOpTB = pStTB * pStTB.dag()


H0 = -omegas[0] / 2.0 * sz1 + gs[0] * (sp1*smTB + sm1*spTB) - omegas[1] / 2.0 * sz2 + gs[1] * (sp2*smTB + sm2*spTB)

H1 = -omegas[2] / 2.0 * szTB

# Rimliga startvärden?
Theta = -0.108 * PHI0 # från McKay
delta = 0.155 * PHI0 # från McKay
omegaPhi = 2*PI * 850.6e6 # från McKay

def H1_coeff(t, args):
    return np.sqrt(np.abs(np.cos(PI * Phi(t, Theta, delta, omegaPhi) / PHI0)))

H = [H0, [H1, H1_coeff]]

output = mcsolve(H, psi0, ts, e_ops=[pOp1, pOp2, pOpTB]) # detta funkar inte?
plotExpect(output)