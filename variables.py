from qutip import *
import numpy as np

# Vi vill ha enheter s.a. hbar = 1, Phi0 = 1, och använder ns som tidsenhet (=> GHz, Gigarad/s som frekvensenheter)

PI = np.pi
# HBAR  = 1.054571817e-34
# Q_e = 1.602176634e-19
# PHI0 = HBAR * PI / Q_e

# These are in Gigarad/s:
omegas = 2*PI * np.array([4.2 ,3.8, 7.5])
gs = 2*PI * np.array([0.050, 0.050])
alpha = -2*PI * 0.150

gSt = basis(2,0)
eSt = basis(2,1)

sp = sigmap()
sm = sigmam()
sz = sigmaz()
sx = sigmax()
sy = sigmay()
I = qeye(2)

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