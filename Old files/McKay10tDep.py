from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1 import *

ts = np.linspace(0,1,1000)
psi0 = tensor(eSt,gSt,gSt) #100

pSt1 = tensor(eSt,gSt,gSt) # 100
pOp1 = pSt1 * pSt1.dag()
pSt2 = tensor(gSt,eSt,gSt) # 010
pOp2 = pSt2 * pSt2.dag()
pStTB = tensor(gSt,gSt,eSt) # 001
pOpTB = pStTB * pStTB.dag()


k0 = 10.0 # k0 = deltamax * dJ/dPhi
tRise = 0.1
tWait = PI / (2.0*k0) - tRise

H1 = 0.5 * (sx1 * sx2 + sy1 * sy2)

smoothness = 1
def H1_coeff(t, args):
    return k0 * ( smoothstep(t, 0, tRise, smoothness) - smoothstep(t, tWait + tRise, tWait + 2*tRise, smoothness) )

H = [H1, H1_coeff]

output = sesolve(H, psi0, ts, e_ops=[pOp1, pOp2, pOpTB])
plotExpect(output)

Xtarg = tensor(gSt,eSt,gSt)
XofT = mesolve(H, psi0, ts)

print(XofT.times[-80])

Xfinal = XofT.states[-80]
print(Xfinal)

print(abs((Xtarg.dag() * Xfinal).tr()))
print(abs(Xfinal.overlap(Xtarg)))
print(fidelity(Xfinal, Xtarg))