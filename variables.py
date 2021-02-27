from qutip import *
import numpy as np

HBAR  = 1.054571817e-34
PI = np.pi

# These are in rad/s:
omegas = 2*PI * np.array([4.2e9 ,3.8e9, 7.5e9])
gs = 2*PI * np.array([50e6, 50e6])
alpha = -2*PI * 150e6

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