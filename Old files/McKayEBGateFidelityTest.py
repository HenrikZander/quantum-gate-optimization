from qutip import *
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay11EB_3lvl import *

x0417DE = [0.38416891, 0.18051406, 3.84487738, 30.23212348, 101.4918881]
x0418DE = [0.35779783, 0.16539164, 1.57643895, 28.25386814, 95.16742431]
x0418DEb = [-0.422849290, 0.0583265568, 4.40324141, 50.0444593, 131.164936]

x0424DE = [0.241082420, 0.108883457, 2.50702612, 39.0692027, 218.127752]
x0425DE = [-0.431076591, 0.0275072184, 3.52975522, 36.7717108, 83.9876656]

xUsed = x0425DE

start=time.time()
for _ in range(10):
    gfid = getGateFidelity(xUsed,wantiSWAP=True)
print(f'Time: {time.time()-start}')

#gfid = getGateFidelity(xUsed,wantiSWAP=True)
print(gfid)