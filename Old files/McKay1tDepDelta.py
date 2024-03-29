from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from variables import *
from functions import *
from plotting import *
from McKay1 import *

ts = np.linspace(0,500,2000)
psi0 = tensor(eSt,gSt,gSt)

pOps = getAllProjectionOperators()


x1 = [1.0-0.56682599, 0.2387175, 2.52545336, 47.62454201] # Theta, delta, omegaPhi, omegaTB0
x2 = [2.0-1.42465745, 0.57298549, 2.52355871, 49.01101584] # Theta, delta, omegaPhi, omegaTB0

H = getSStepHamiltonian(x1,265.6)

output = sesolve(H, psi0, ts, e_ops=pOps) # Verkar göra det den ska
plotExpect(output)