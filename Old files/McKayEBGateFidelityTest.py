from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from model import *

x0417DE = [0.38416891, 0.18051406, 3.84487738, 30.23212348, 101.4918881]
x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]
x0418DEb = [-0.422849290, 0.0583265568, 4.40324141, 50.0444593, 131.164936]

x0424DE = [0.241082420, 0.108883457, 2.50702612, 39.0692027, 218.127752]
x0425DE = [-0.431076591, 0.0275072184, 3.52975522, 36.7717108, 83.9876656]

x0426CZ_4lvl = [4.44287328e-01, 2.05934595e-02, 2.04012274e+00, 3.07227234e+01, 1.78093842e+02]
x0426CZ_3lvl = [-2.58309069e-01, 5.46369915e-02, 4.42420464e+00, 2.83414251e+01, 1.79648639e+02]

xUsed = x0425DE
'''
start=time.time()
for _ in range(10):
    gfid = getGateFidelity(xUsed,wantiSWAP=True)
print(f'Time: {time.time()-start}')
'''
gfid = getGateFidelity(xUsed,N=4,wantiSWAP=True)
print(gfid)