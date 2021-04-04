from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from skopt import gp_minimize
import time
import scipy
from plotting import *


i = 0
global maxRuntime
maxRuntime = 10800

def smoothstep(x, x_min=0, x_max=1, N=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0
    for n in range(N + 1):
         result += comb(N + n, n) * comb(2 * N + 1, N - n) * (-x) ** n
    result *= x ** (N + 1)
    return result

def sinstep(x, x_min=0, x_max=1):
    x = np.clip((x - x_min) / (x_max - x_min), 0, 1)
    result = 0.5 - 0.5*np.cos(np.pi*x)
    return result

def getEigenstates(x, H_const, H_omegaTB):
    H = H_const + x[3]*np.sqrt(np.abs(np.cos(np.pi*x[0]))) * H_omegaTB
    return H.eigenstates()


def evaluateResult(x, fun, resultList, N=5):
    if len(resultList[1]) < N:
        resultList[0].append(x)
        resultList[1].append(fun)
    else:
        if np.amax(resultList[1]) > fun:
            indexOfMax = np.argmax(resultList[1])
            resultList[0][indexOfMax] = x
            resultList[1][indexOfMax] = fun
    return resultList


def bayesianCallback(result):
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    fun = result.fun
    x = result.x
    
    bestResults = evaluateResult(x, fun, bestResults)
    print(f'Num of found minima: {i+1}')
    print(f'The Bayesian Optimization algorithm gave a minimum of {fun} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = []
        for i in range(len(bestResults[1])):
            result.append((bestResults[0][i],bestResults[1][i]))
        saveResToFile(result, "Bayesian Optimization", i, passedTime)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False

    
def callbackBH(x, fun, accept):
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    bestResults = evaluateResult(x, fun, bestResults)
    print(f'Num of minima found: {i+1}')
    print(f'The Basin-Hopping algorithm gave a minimum of {fun} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = []
        for i in range(len(bestResults[1])):
            result.append((bestResults[0][i],bestResults[1][i]))
        saveResToFile(result, "Basin-Hopping", i, passedTime)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def callbackDE(x,convergence=None):
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    print(f'Num of iterations: {i+1}')
    print(f'The currently best minimum that the Differential Evolution algorithm has found has a convergence of {convergence} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = [(x,convergence)]
        saveResToFile(result, "Differential Evolution", i, passedTime, algorithmDE=True)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def callbackDA(x, fun, context):
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    bestResults = evaluateResult(x, fun, bestResults)
    print(f'Num of minima found: {i+1}')
    print(f'The Dual Anneling algorithm gave a minimum of {fun} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = []
        for i in range(len(bestResults[1])):
            result.append((bestResults[0][i],bestResults[1][i]))
        saveResToFile(result, "Dual Anneling", i, passedTime)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def callbackSHG(x):
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[],[]]
    
    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime
    
    print(f'Num of iterations: {i+1}')
    print(f'The current point of examination by the Simplicial Homology Global algorithm is {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        i = 0
        print("Optimization timeout triggered!")
        return True
    else: 
        return False


def findMinimum(costFunction, bounds, x0=None, runBayesian=False, runSHG=True, runDA=True, runDE=True, runBH=True, runBayesianWithBH=False, numOfMinimaToFindBH = 500):
    """
    This function finds the minimum of the specified cost function, either within 
    the given bounds or close to the inital guess 'x0' depending on the algorithm 
    in use. A number of different algorithms can be used simultaneously.
    """
    
    global startTime 
    
    #Initial formatting tasks.
    print("")
    message = "##################################################\n"
    setInitialGuess = False
    coordinatesOrigin = [0]*len(bounds)
    result = []
    algorithmsUsed = []
    runtime = []
    
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
        resBayesian = gp_minimize(costFunction, bounds, callback=bayesianCallback, n_jobs=-1, n_calls=4000, n_initial_points=100)
        timeBayesian = time.time() - startTime
        message += f'The optimizaton using \"gp_minimize()\" took {round(timeBayesian,2)}s to execute and ended on a minimum of {resBayesian.fun} at the point {resBayesian.x}.\n'
        result.append(resBayesian)
        algorithmsUsed.append("Bayesian Optimization")
        runtime.append(timeBayesian)
    
    #Set the initial guess for the Basin-hopping algorithm to the result of the Bayesian algorithm.
    if setInitialGuess:
        x0 = resBayesian.x
    
    #Optimization using the Simplicial Homology Global algorithm.
    if runSHG:
        startTime = time.time()
        resSHG = scipy.optimize.shgo(costFunction, bounds, iters=4, callback=callbackSHG)
        timeSHG = time.time() - startTime
        message += f'The optimizaton using the \"Simplicial Homology Global\"-algorithm took {round(timeSHG,2)}s to execute and ended on a minimum of {resSHG.fun} at the point {resSHG.x}.\n'
        message += f'Function evaluations performed: {resSHG.nfev}\n'
        result.append(resSHG)
        algorithmsUsed.append("Simplicial Homology Global")
        runtime.append(timeSHG)
    
    #Optimization using the Dual Annealing algorithm. 
    if runDA:
        startTime = time.time()
        resDA = scipy.optimize.dual_annealing(costFunction, bounds, callback=callbackDA)
        timeDA = time.time() - startTime
        message += f'The optimizaton using the \"Dual Annealing\"-algorithm took {round(timeDA,2)}s to execute and ended on a minimum of {resDA.fun} at the point {resDA.x}.\n'
        message += f'Function evaluations performed: {resDA.nfev}\n'
        result.append(resDA)
        algorithmsUsed.append("Dual Annealing")
        runtime.append(timeDA)
    
    #Optimization using the Differential Evolution algorithm.
    if runDE:
        startTime = time.time()
        resDE = scipy.optimize.differential_evolution(costFunction, bounds, callback=callbackDE)
        timeDE = time.time() - startTime
        message += f'The optimizaton using the \"Differential Evolution\"-algorithm took {round(timeDE,2)}s to execute and ended on a minimum of {resDE.fun} at the point {resDE.x}.\n'
        message += f'Function evaluations performed: {resDE.nfev}\n'
        result.append(resDE)
        algorithmsUsed.append("Differential Evolution")
        runtime.append(timeDE)
    
    #Optimization using the Basin-hopping algorithm.
    if runBH:
        startTime = time.time()
        resBH = scipy.optimize.basinhopping(costFunction, x0, niter=numOfMinimaToFindBH, callback=callbackBH)
        timeBH = time.time() - startTime
        message += f'The optimizaton using the \"Basin-Hopping\"-algorithm took {round(timeBH,2)}s to execute and ended on a minimum of {resBH.fun} at the point {resBH.x}.\n'
        message += f'Function evaluations performed: {resBH.nfev}\n'
        result.append(resBH)
        algorithmsUsed.append("Basin-Hopping")
        runtime.append(timeBH)
        
    if runBayesianWithBH:
        message += "##################################################\n"
        message += f'The optimizaton using the \"Basin-hopping\"-algorithm together with an initial guess from the \"gp_minimize()\" algorithm took {round(timeBH + timeBayesian,2)}s in total to execute and ended on a minimum of {resBH.fun} at the point {resBH.x}.\n'
    
    print("")        
    print(message + "##################################################")
    saveAllFinalResults(result, algorithmsUsed, runtime)
    return result


def generateCostFunction(hamiltonian, projectionOperators, useGateFidelity=False, initialState=None, timeStamps=None, maximumGateTime=250):
    if timeStamps is None:
        timeStamps = np.linspace(0,maximumGateTime,maximumGateTime*3)
    
    def costFunction(x):
        H = hamiltonian(x)
        options = solver.Options()
        options.nsteps = 10000
        if useGateFidelity:
            #Generate a hamiltonian based on the parameters in the list x.
            #Calculate the time evolution of the Hamiltonian
            #Calculate the maximum gate fidelity
            return 0#gateFidelity
        elif not initialState is None:
            result = sesolve(H, initialState, timeStamps, projectionOperators, options=options)
            allExpectedValues = result.expect
            expectValue = -np.amax(allExpectedValues[0])
            return expectValue
        else:
            print("Calculation failed! Check format of initial quantum state if you have set useGateFidelity equal to True.")
            return None
    
    return costFunction


def optimizeGate(hamiltonianModule, maximumGateTime=250, runBayesian=False, runSHG=False, runDA=False, runDE=False, runBH=False, runBayesianWithBH=False):
    hamiltonian = hamiltonianModule.getHamiltonian
    parameterBounds = hamiltonianModule.getParameterBounds()
    projectionOperators = hamiltonianModule.getProjectionOperators()
    initialState = hamiltonianModule.getInitialState()
    initialGuess = hamiltonianModule.getInitialGuess()
    
    costFunction = generateCostFunction(hamiltonian, projectionOperators, initialState=initialState, maximumGateTime=maximumGateTime)
    findMinimum(costFunction, parameterBounds, initialGuess, runBayesian=runBayesian, runSHG=runSHG, runDA=runDA, runDE=runDE, runBH=runBH, runBayesianWithBH=runBayesianWithBH)

    
def simulateHamiltonian(hamiltonianModule, x0, simulationTime=500):
    hamiltonian = hamiltonianModule.getHamiltonian(x0)
    initialState = hamiltonianModule.getInitialState()
    projectionOperators = hamiltonianModule.getAllProjectionOperators()
    
    options = solver.Options()
    options.nsteps = 10000
    timeStamps = np.linspace(0,simulationTime,simulationTime*3)
    
    result = sesolve(hamiltonian, initialState, timeStamps, projectionOperators, options=options)
    plotExpect(result)
    