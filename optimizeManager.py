######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: optimizeManager.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 19 July 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import time
import optimize2QubitFrame as gui
import dataManager
import model
import numpy as np
import scipy
from datetime import datetime

######################################################################################################################################################################
# Global variables

i = 0
global maxRuntime
maxRuntime = 40000

######################################################################################################################################################################
# Cost function definition


def averageFidelity(F, gateTimeWeight=2):
    F[0] = gateTimeWeight*F[0]
    return np.sum(F)/(len(F)+gateTimeWeight-1)


def cost(x, gateType, N, circuitData, signalType):
    F, _ = model.getGateFidelity(x, gateType=gateType, N=N, circuitData=circuitData, printResults=True, signalType=signalType)
    return -averageFidelity(F)


######################################################################################################################################################################
# Optimize gate function


def optimize2QubitGate(userData):
    """
    The function tries to optimize the selected gate for the
    supplied parameters, using the optimization algorithms 
    that the user decides would best fit the problem.
    ---------------------------------------------------------
    INPUT:
            energyLevels (int) {Optional}: How many energy levels that should be accounted for in the simulations.
            timeoutAlgorithm (int) {Optional}: How many seconds each optimization algorithm can run, before a forced termination is initiated.
            maxAllowedGateTime (int) {Optional}: The longest gate time we will allow a solution to have.
            runSHG (boolean) {Optional}: If True, the function will use the Simplicial Homology Global algorithm to optimize the gate.
            runDA (boolean) {Optional}: If True, the function will use the Dual Anneling algorithm to optimize the gate.
            runDE (boolean) {Optional}: If True, the function will use the Differential Evolution algorithm to optimize the gate.
            wantTradGate (boolean) {Optional}: If True, the function will restrict the parameter space to focus on ''traditional'', 
                well-understood iSWAP/CZ gates.
            wantCZ_20 (boolean) {Optional}: This boolean is only relevant if the optimizer is looking for a traditional CZ gate.
                If True, the function will try to find a CZ gate using the 11<->20 transition. If False, the 11<->02 transition will be used instead.

        Set only ONE of these to True!:
            iSWAP (boolean) {Optional}: Optimize the iSWAP quantum gate.
            CZ (boolean) {Optional}: Optimize the CZ quantum gate.
    ---------------------------------------------------------
    OUTPUT:
            No output.
    ---------------------------------------------------------
    """
    #
    # Divide data from user into variables
    #
    signalType = userData['signalType']

    if signalType == 'arccos':
        x0name, x1name = "dcAmplitude", "acAmplitude"
    elif signalType == 'cos':
        x0name, x1name = "theta", "delta"

    gateType = userData['gateType']
    parameterBounds = (userData[x0name], userData[x1name], userData["omegaPhi"], userData["modulationTime"])
    energyLevels = userData["energy-levels"]
    runSHG=userData["runSHG"]
    runDA=userData["runDA"]
    runDE=userData["runDE"]

    # Dessa ska tas in mha gui:t
    solutionsFolder = userData["save-folder"] # "Qubit Pair 03/Solutions"

    findMinimum(cost, parameterBounds, argumentsToOptimizer=(gateType, energyLevels, userData, signalType), runSHG=runSHG, runDA=runDA, runDE=runDE, solutionsFolder=solutionsFolder)
    gui.enableStopButton()
    if gui.getRunOptimizer():
        gui.processFinished()
    else:
        gui.processStopped()


######################################################################################################################################################################
# Callback functions


def callbackDE(x, convergence=None):
    """
    Callback function for the Differential Evolution optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[], []]

    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime

    print('##################################################')
    print(f'Num of iterations: {i+1}')
    print(f'The best gate that the Differential Evolution algorithm has found thus far lies at the point {x}, and the current convergence is {convergence}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    print('##################################################')
    i += 1

    gui.setProgressValue(convergence)

    if convergence > 1:
        gui.setStatus("Status: Finished! Polishing solution.")
        gui.disableStartStopButtons()
    
    if gui.getRunOptimizer():
        return False
    else:
        return True


def callbackDA(x, fun, context):
    """
    Callback function for the Dual Annealing optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[], []]

    currentTime = time.time()
    passedTime = currentTime - startTime
    iterTime = currentTime - lastTime
    lastTime = currentTime

    bestResults = dataManager.evaluateResult(x, fun, bestResults)
    print(f'Num of minima found: {i+1}')
    print(f'The Dual Annealing algorithm gave a minimum of {fun} at the point {x}.')
    print(f'Total time passed: {passedTime} seconds.')
    print(f'Iteration time: {iterTime} seconds.\n')
    i += 1
    if passedTime > maxRuntime:
        result = []
        for i in range(len(bestResults[1])):
            result.append((bestResults[0][i], bestResults[1][i]))
        dataManager.saveResToFile(result, "Dual Annealing", i, passedTime)
        i = 0
        print("Optimization timeout triggered!")
        return True
    else:
        return False


def callbackSHG(x):
    """
    Callback function for the Simplical Homology Global optimizer.
    """
    global i
    global lastTime
    global bestResults
    if i == 0:
        lastTime = startTime
        bestResults = [[], []]

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


######################################################################################################################################################################


def findMinimum(costFunction, bounds, argumentsToOptimizer, runSHG=False, runDA=False, runDE=True, solutionsFolder=None):
    """
    This function aims to find the global minimum of the specified cost 
    function, within the parameter bounds that has been specified. A 
    number of different algorithms are implemented and they can be used 
    consecutively during one function call.
    ---------------------------------------------------------
    INPUT:
            costFunction (function pointer): The cost function that is to be minimized. The input parameters to the cost function should be specified as an array x.
            bounds (array(tuple(int))): The parameter bounds that the optimizer will do it's best to not leave. The length of this array determines the expected dimension of the cost function input x.
            runSHG (boolean) {Optional}: If True the function will use the Simplicial Homology Global algorithm to minimize the cost function.
            runDA (boolean) {Optional}: If True the function will use the Dual Anneling algorithm to minimize the cost function.
            runDE (boolean) {Optional}: If True the function will use the Differential Evolution algorithm to minimize the cost function.
    ---------------------------------------------------------
    OUTPUT:
            result (array(scipy.optimize.OptimizeResult)): An array of Scipy result structures. See Scipy documentation for more information.
    ---------------------------------------------------------
    """

    global startTime
    global i
    i = 0
    # Initial formatting tasks for the result that will be printed in the terminal.
    print("")
    message = "##################################################\n"
    result = []
    algorithmsUsed = []
    runtime = []

    '''
    if argumentsToOptimizer[5]:
        # Constraining A + B to be \in [0,1] when an arccos signal is used.
        def constraintFunction(x):
            return np.array(x[0] + x[1])
        nlc = scipy.optimize.NonlinearConstraint(constraintFunction, 0, 1)
    '''
    
    # Optimization using the Simplicial Homology Global algorithm.
    if runSHG:
        startTime = time.time()
        resSHG = scipy.optimize.shgo(costFunction, bounds, iters=4, callback=callbackSHG)
        timeSHG = time.time() - startTime
        message += f'The optimizaton using the \"Simplicial Homology Global\"-algorithm took {round(timeSHG,2)}s to execute and ended on a minimum of {resSHG.fun} at the point {resSHG.x}.\n'
        message += f'Function evaluations performed: {resSHG.nfev}\n'
        result.append(resSHG)
        algorithmsUsed.append("Simplicial Homology Global")
        runtime.append(timeSHG)

    i = 0

    # Optimization using the Dual Annealing algorithm.
    if runDA:
        startTime = time.time()
        resDA = scipy.optimize.dual_annealing(costFunction, bounds, callback=callbackDA)
        timeDA = time.time() - startTime
        message += f'The optimizaton using the \"Dual Annealing\"-algorithm took {round(timeDA,2)}s to execute and ended on a minimum of {resDA.fun} at the point {resDA.x}.\n'
        message += f'Function evaluations performed: {resDA.nfev}\n'
        result.append(resDA)
        algorithmsUsed.append("Dual Annealing")
        runtime.append(timeDA)

    i = 0

    # Optimization using the Differential Evolution algorithm.
    if runDE: #  constraints=(nlc),
        startTime = time.time()
        resDE = scipy.optimize.differential_evolution(costFunction, bounds, callback=callbackDE, workers=-1, updating='deferred', maxiter=100000, args=argumentsToOptimizer, polish=True)
        timeDE = time.time() - startTime
        message += f'The optimizaton using the \"Differential Evolution\"-algorithm took {round(timeDE,2)}s to execute and ended on a minimum of {resDE.fun} at the point {resDE.x}.\n'
        message += f'Function evaluations performed: {resDE.nfev}\n'

        ## Start of modified code snippet from the polishing done by scipy.differential_evolution
        polish_method = 'L-BFGS-B'

        # We aren't currently using constraints
        '''
        if self._wrapped_constraints:
            polish_method = 'trust-constr'

            constr_violation = self._constraint_violation_fn(DE_result.x)
            if np.any(constr_violation > 0.):
                warnings.warn("differential evolution didn't find a"
                                " solution satisfying the constraints,"
                                " attempting to polish from the least"
                                " infeasible solution", UserWarning)
        '''
        # Manual polish, no constraints
        polStart = time.time()
        resPol = scipy.optimize.minimize(costFunction,
                            np.copy(resDE.x),
                            method=polish_method,
                            bounds=bounds, #self.limits.T,
                            args=argumentsToOptimizer)
                            #callback=callbackDE) # Perhaps we want callback?
                            # constraints=self.constraints)
        timePol = time.time() - polStart
        #print('Polishing done!\n')
        message += f'The polishing took {round(timePol,2)}s to execute and ended on a minimum of {resPol.fun} at the point {resPol.x}.\n'
        message += f'Function evaluations performed: {resPol.nfev}\n'
        totTime = time.time()-startTime
        message += f'Total time elapsed: {round(totTime,2)}s\n'
        message += f'Total number of function evaluations: {resDE.nfev + resPol.nfev}\n'
        # Polished solution is only accepted if there is an improvement in
        # cost function, the polishing was successful and the solution lies
        # within the bounds.
        if (resPol.fun < resDE.fun and
                resPol.success and
                np.all(resPol.x <= np.asarray(bounds)[:,1]) and
                np.all(np.asarray(bounds)[:,0] <= resPol.x)):
            saveRes = resPol
        else:
            print("\nDidn't save polished solution!")
            print([resPol.fun < resDE.fun, resPol.success, np.all(resPol.x <= np.asarray(bounds)[:,1]), np.all(np.asarray(bounds)[:,0] <= resPol.x)])
            saveRes = resDE
        ## End of modified code snippet from the polishing done by scipy.differential_evolution

        result.append(saveRes)
        algorithmsUsed.append("Differential Evolution")
        runtime.append(totTime)
    
    print("")
    print(message + "##################################################")
    dateAndTime = datetime.today()

    ############ Saving to result.txt (this should be removed later) ############
    saveAllFinalResults(result, algorithmsUsed, runtime, dateAndTime=dateAndTime)
    #############################################################################
    
    if solutionsFolder is not None:
        gateType, N, circuitData, signalType = argumentsToOptimizer
        dataManager.saveSolutionsTojson(result, gateType, N, solutionsFolder, circuitData=circuitData, dateAndTime=dateAndTime, signalType=signalType)
    return result, dateAndTime


def saveAllFinalResults(results, algorithm, runtime, fileName="result.txt", dateAndTime=datetime.today()):
    resultFile = open(fileName, "a")
    for i in range(len(results)):
        dateStr = f'Result (full) from: {dateAndTime}\n'
        algorithmStr = f'Algorithm: {algorithm[i]}\n'
        runtimeStr = f'Total runtime: {runtime[i]} seconds\n'
        dividerStr = "##################################################\n"
        resultStr = str(results[i])
        
        strList = [dividerStr, dateStr, algorithmStr, runtimeStr, dividerStr, resultStr, "\n", dividerStr]
        resultFile.writelines(strList)
    resultFile.close()

######################################################################################################################################################################