from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1EB import *

x0417DE = [0.38416891, 0.18051406, 3.84487738, 30.23212348, 101.4918881]
x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]
x0418DEb = [-0.422849290, 0.0583265568, 4.40324141, 50.0444593, 131.164936]

xUsed = x0418DE

# HBBComps = getHamiltonian(xUsed,getBBHamiltonianComps=True)
# eigSts = getThetaEigenstates(xUsed,HBBComps[0],HBBComps[1])

# print(eigSts[0])

start=time.time()
for _ in range(10):
    gfid = getGateFidelity(xUsed,wantiSWAP=True)
print(f'Time: {time.time()-start}')
print(gfid)