from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from skopt import gp_minimize
import time
import scipy

from variables import *
from plotting import *

i=0

def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result

"""
def hamiltonian(omega, g):
    return HBAR*( (-omega[0]/2)*sz1 + (-omega[1]/2)*sz2 + g[0]*(sp1*smTB + sm1*spTB) + g[1]*(sp2*smTB + sm2*spTB) + (-omega[2]/2)*szTB )
"""

def testHamiltonian():
    return 2 * np.pi * 0.1 * sigmax()

"""
def Phi(t):
    return Theta + delta*np.cos(omegaPhi*t)


def omegaTB(t, args):
    return omegas[2]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
"""

def findMinimum(costFunction, bounds, x0=None, runBayesian=False, runSHG=True, runDA=True, runDE=True, runBH=True, runBayesianWithBH=False, numOfIterBH = 500):
    """
    This function finds the minimum of the specified cost function, either within 
    the given bounds or close to the inital guess 'x0' depending on the algorithm 
    in use. A number of different algorithms can be used simultaneously.
    """
    
    #Initial formatting tasks.
    print("")
    message = "##################################################\n"
    setInitialGuess = False
    coordinatesOrigin = [0]*len(bounds)
    result = []
    
    #If the initial guess for the Basin-hopping algorithm hasen't been specified set the guess to the origin.
    if x0 is None:
        x0 = [0]*len(bounds)
    
    #If we want to use the result from the Bayesian algorithm as our starting point for the Basin-hopping algorithm, 
    #run the Bayesian algorithm, set the result to the initial guess and then run the Basin-hopping algoritm.
    #Also reduce the amount of iterations to be performed by the Basin-hopping algorithm.
    if runBayesianWithBH:
        setInitialGuess = True
        runBayesian = True
        runBH = True
        numOfIterBH = 100
    
    #Optimization using the Bayesian optimization algoritm.
    if runBayesian:
        startTime = time.time()
        resBayesian = gp_minimize(costFunction, bounds)
        timeBayesian = time.time() - startTime
        message += f'The optimizaton using \"gp_minimize()\" took {round(timeBayesian,2)}s to execute and ended on a minimum of {resBayesian.fun} at the point {resBayesian.x}.\n'
        result.append(resBayesian)
    
    #Set the initial guess for the Basin-hopping algorithm to the result of the Bayesian algorithm.
    if setInitialGuess:
        x0 = resBayesian.x
    
    #Optimization using the Simplicial Homology Global algorithm.
    if runSHG:
        startTime = time.time()
        resSHG = scipy.optimize.shgo(costFunction, bounds, iters=4)
        timeSHG = time.time() - startTime
        message += f'The optimizaton using the \"Simplicial Homology Global\"-algorithm took {round(timeSHG,2)}s to execute and ended on a minimum of {resSHG.fun} at the point {resSHG.x}.\n'
        message += f'Function evaluations performed: {resSHG.nfev}\n'
        result.append(resSHG)
    
    #Optimization using the Dual Annealing algorithm. 
    if runDA:
        startTime = time.time()
        resDA = scipy.optimize.dual_annealing(costFunction, bounds)
        timeDA = time.time() - startTime
        message += f'The optimizaton using the \"Dual Annealing\"-algorithm took {round(timeDA,2)}s to execute and ended on a minimum of {resDA.fun} at the point {resDA.x}.\n'
        message += f'Function evaluations performed: {resDA.nfev}\n'
        result.append(resDA)
    
    #Optimization using the Differential Evolution algorithm.
    if runDE:
        startTime = time.time()
        resDE = scipy.optimize.differential_evolution(costFunction, bounds)
        timeDE = time.time() - startTime
        message += f'The optimizaton using the \"Differential Evolution\"-algorithm took {round(timeDE,2)}s to execute and ended on a minimum of {resDE.fun} at the point {resDE.x}.\n'
        message += f'Function evaluations performed: {resDE.nfev}\n'
        result.append(resDE)
    
    #Optimization using the Basin-hopping algorithm.
    if runBH:
        startTime = time.time()
        resBH = scipy.optimize.basinhopping(costFunction, x0, niter=numOfIterBH)
        timeBH = time.time() - startTime
        message += f'The optimizaton using the \"Basin-hopping\"-algorithm took {round(timeBH,2)}s to execute and ended on a minimum of {resBH.fun} at the point {resBH.x}.\n'
        message += f'Function evaluations performed: {resBH.nfev}\n'
        result.append(resBH)
        
    if runBayesianWithBH:
        message += "##################################################\n"
        message += f'The optimizaton using the \"Basin-hopping\"-algorithm together with an initial guess from the \"gp_minimize()\" algorithm took {round(timeBH + timeBayesian,2)}s in total to execute and ended on a minimum of {resBH.fun} at the point {resBH.x}.\n'
            
    print(message + "##################################################")
    return result


def generateCostFunction(hamiltonian, projectionOperators, useGateFidelity=False, initialState=None, timeStamps=None):
    if timeStamps is None:
        timeStamps = np.linspace(0,1700,850)
    
    def costFunction(x):
        H = hamiltonian(x)
        print(x)
        options = solver.Options()
        options.nsteps = 3000
        if useGateFidelity:
            #Generate a hamiltonian based on the parameters in the list x.
            #Calculate the time evolution of the Hamiltonian
            #Calculate the maximum gate fidelity
            return 0#gateFidelity
        elif not initialState is None:
            result = sesolve(H, initialState, timeStamps, projectionOperators, options=options)
            allExpectedValues = result.expect
            #plotExpect(result)
            expectValue = np.amin(allExpectedValues[0])
            print(expectValue)
            return expectValue
        else:
            print("Calculation failed! Check format of initial quantum state if you have set useGateFidelity equal to True.")
            return None
    
    return costFunction


def optimizeGate(hamiltonian, parameterBounds, projectionOperators, initialGuess=None, initialState=None, runBayesian=False, runSHG=False, runDA=False, runDE=False, runBH=False, runBayesianWithBH=False):
    costFunction = generateCostFunction(hamiltonian, projectionOperators, initialState=initialState)
    result = findMinimum(costFunction, parameterBounds, initialGuess, runBayesian=runBayesian, runSHG=runSHG, runDA=runDA, runDE=runDE, runBH=runBH, runBayesianWithBH=runBayesianWithBH)
    
    
    
    
    
    
    
    
    
