######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: plotting.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 27 February 2021

# Last modified: 29 April 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

from qutip import *
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

######################################################################################################################################################################
# Functions for saving optimization results to a .txt file.

def saveAllFinalResults(results, algorithm, runtime, fileName="result.txt"):
    resultFile = open(fileName, "a")
    todaysDate = datetime.now()
    for i in range(len(results)):
        dateStr = f'Result (full) from: {todaysDate}\n'
        algorithmStr = f'Algorithm: {algorithm[i]}\n'
        runtimeStr = f'Total runtime: {runtime[i]} seconds\n'
        dividerStr = "##################################################\n"
        resultStr = str(results[i])
        
        strList = [dividerStr, dateStr, algorithmStr, runtimeStr, dividerStr, resultStr, "\n", dividerStr]
        resultFile.writelines(strList)
    resultFile.close()
    
    
def saveResToFile(result, algorithmName, iterations, runtime, algorithmDE=False, algorithmSHG=False, fileName="result.txt"):
    resultFile = open(fileName, "a")
    todaysDate = datetime.now()
    
    dateStr = f'Result (timeout) from: {todaysDate}\n'
    iterStr = f'Total iterations performed (or minima found): {iterations}\n'
    runtimeStr = f'Total runtime for algorithm: {runtime} seconds.\n'
    dividerStr = "##################################################\n"
    strList = [dividerStr, dateStr, iterStr, runtimeStr, dividerStr, "\n"]
    resultFile.writelines(strList)
    
    for res in result:
        x = res[0]
        fun = res[1]
        
        if algorithmDE:
            resultStr = f'The {algorithmName} algorithm gave a minimum at the point {x} with a convergence of {fun}.\n'
        elif algorithmSHG:
            resultStr = f'The {algorithmName} algorithm was examining the point {x}.\n'
        else: 
            resultStr = f'The {algorithmName} algorithm gave a minimum of {fun} at the point {x}.\n'
        
        resultFile.writelines([resultStr, "\n"])
    resultFile.write(dividerStr)
    resultFile.close()


######################################################################################################################################################################
# Functions that handle plotting of states, expectation value and similar.


def plotStates(result):
    """Plots all the states in the "result"-structure on the same Bloch-sphere."""
    states = result.states
    fig = Bloch()
    fig.add_states(states)
    fig.make_sphere()
    plt.show()


def plotExpect(timeStamps, expectationValues):
    """Plots the expectation values for an arbitrary number of projection operators."""
    fig, ax = plt.subplots()
    labels = []
    for i, e in enumerate(expectationValues):
        ax.plot(timeStamps, e)
        labels.append("Operator " + str(i))
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')
    ax.legend(labels)
    plt.show()

# plotEigenExpect nedan är egentligen inte särskilt nödvändig, tror jag, men jag låter den ligga kvar här tillsvidare.

# Ideally, Phis should be calculated from x and result.times, 
# but this requires specifying whether Phi is smoothstepped or not,
# so I included it as an input variable instead:
def plotEigenExpect(result, Phis, x, eigenStateIndices=[2, 1, 3]):
    """Plots the expectation values for an arbitrary number of eigenstate projection operators."""
    stateEvolution = result.states
    times = result.times

    expvalsList = []

    # The nested for loops below make a nested list expvalsList,
    # where expvalsList[k] consists of the expectation values of the (time dependent) projection operator
    # projecting onto the eigenstate corresponding to eigenStateIndices[k],
    # i.e. expvalsList[k](t) = |<Psi(t)|eigpOp_eigenStateIndices[k](t)>|^2

    for k in range(len(eigenStateIndices)):
        expvalsList.append([])
        for tIndex in range(len(times)):
            eigpOp = getEigenProjectionOperator(x,Phis[tIndex],eigenStateIndices[k])
            expvalsList[k].append(expect(eigpOp,stateEvolution[tIndex]))

    fig, ax = plt.subplots()
    labels = ["Qubit 1", "Qubit 2", "Coupler"] #[]
    #i = 1
    for e in expvalsList:
        ax.plot(times, e)
        #labels.append("Z Projection " + str(i))
        #i = i + 1
    ax.set_xlabel('Time')
    ax.set_ylabel('Expectation values')
    ax.legend(labels)
    plt.show()


def createFilename(path,index):
    """Creates a filename for a .png-file using an index and the path to where the file should be placed."""
    path = path + "\\"
    numStr = str(index)
    if index < 10:
        path = path + "00" + numStr
    elif index < 100:
        path = path + "0" + numStr
    else:
        path = path + numStr
    return path + ".png"


def createAnimationVideo(states, dirc, framerate=60):
    """Creates the frames for the animation and stiches them together using the ffmpeg software."""
    b = Bloch()
    b.vector_color = ['r']
    b.view = [-40,30]
    currentPath = os.getcwd()
    path = currentPath + "\\" + dirc
    
    for i in range(len(states)):
        b.clear()
        b.add_states(states[i])
        filename = createFilename(path, i)
        b.save(name=filename)
    os.system('ffmpeg -r ' + str(framerate) + ' -i ' + '.\\' + dirc + '\\' + '%03d.png -vf scale=-2:1080 -pix_fmt yuv420p bloch.mp4')
    #Idea: create a new folder if this one already exsists.
    #ffmpeg -r 20 -i %03d.png -vf scale=-2:1080 -pix_fmt yuv420p bloch.mp4
    

def animateStates(result, dirc):
    """Animates two-level state vectors on the Bloch sphere."""
    states = result.states
    createAnimationVideo(states, dirc)
    
def plotVectors(x,y):
    """Plots the values in the lists x and y."""
    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()


######################################################################################################################################################################