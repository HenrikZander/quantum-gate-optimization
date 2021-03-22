from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from skopt import gp_minimize
import time
import scipy

from variables import *


def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result

#def hamiltonian(omega, g):
#    return HBAR*( (-omega[0]/2)*sz1 + (-omega[1]/2)*sz2 + g[0]*(sp1*smTB + sm1*spTB) + g[1]*(sp2*smTB + sm2*spTB) + (-omega[2]/2)*szTB )


def testHamiltonian():
    return 2 * np.pi * 0.1 * sigmax()


def Phi(t):
    return Theta + delta*np.cos(omegaPhi*t)


def omegaTB(t, args):
    return omegas[2]*np.sqrt(np.abs(np.cos(PI*Phi(t))))


def optimizeGate(costFunction, bounds, runBayesian=False, runSHG=True, runDA=True, runDE=True, runBH=True):
    print("")
    message = "##################################################\n"
    result = []
    
    #Optimization using the Bayesian optimization algoritm
    if(runBayesian):
        startTime = time.time()
        resBayesian = gp_minimize(costFunction, bounds)
        timeBayesian = time.time() - startTime
        message += f'The optimizaton using \"gp_minimize()\" took {round(timeBayesian,2)}s to execute and ended on a minimum of {resBayesian.fun} at the point {resBayesian.x}.\n'
        result.append(resBayesian)
    
    if(runSHG):
        startTime = time.time()
        resSHG = scipy.optimize.shgo(costFunction, bounds)
        timeSHG = time.time() - startTime
        message += f'The optimizaton using the \"Simplicial Homology Global\"-algorithm took {round(timeSHG,2)}s to execute and ended on a minimum of {resSHG.fun} at the point {resSHG.x}.\n'
        result.append(resSHG)
    
    if(runDA):
        startTime = time.time()
        resDA = scipy.optimize.dual_annealing(costFunction, bounds)
        timeDA = time.time() - startTime
        message += f'The optimizaton using the \"Dual Annealing\"-algorithm took {round(timeDA,2)}s to execute and ended on a minimum of {resDA.fun} at the point {resDA.x}.\n'
        result.append(resDA)
        
    if(runDE):
        startTime = time.time()
        resDE = scipy.optimize.differential_evolution(costFunction, bounds)
        timeDE = time.time() - startTime
        message += f'The optimizaton using the \"Differential Evolution\"-algorithm took {round(timeDE,2)}s to execute and ended on a minimum of {resDE.fun} at the point {resDE.x}.\n'
        result.append(resDE)
    
    if(runBH):
        startTime = time.time()
        resBH = scipy.optimize.basinhopping(costFunction, (0,0), niter=500)
        timeBH = time.time() - startTime
        message += f'The optimizaton using the \"Basin-hopping\"-algorithm took {round(timeBH,2)}s to execute and ended on a minimum of {resBH.fun} at the point {resBH.x}.\n'
        result.append(resBH)
    
    print(message + "##################################################")
    return result