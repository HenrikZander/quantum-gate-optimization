from qutip import *
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from variables import *

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


def plotStates(result):
    """Plots all the states in the "result"-structure on the same Bloch-sphere."""
    states = result.states
    fig = Bloch()
    fig.add_states(states)
    fig.show()


def plotExpect(result):
    """Plots the expectation values for an arbitraty amount of collapse operators."""
    lists = result.expect
    fig, ax = plt.subplots()
    labels = ["Qubit 1", "Qubit 2", "Coupler"] #[]
    #i = 1
    for e in lists:
        ax.plot(result.times, e)
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