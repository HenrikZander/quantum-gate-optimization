from qutip import *

sp = sigmap()
sm = sigmam()
sz = sigmaz()
sx = sigmax()
sy = sigmay()
I = identity(2)

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

omega = [5e9 ,7e9, 7.5e9]
g = [0.25, 0.25]
delta = 0.1
HBAR  = 1.054571817e-34