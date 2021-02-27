from qutip import *
import os
import numpy as np
import matplotlib.pyplot as plt

from variables import *


def plotStates(result):
    states = result.states
    fig = Bloch()
    fig.add_states(states)
    fig.show()


def plotExpect(result):
    lists = result.expect
    fig, ax = plt.subplots()
    for e in lists:
        ax.plot(result.times, e)
    ax.set_xlabel('Time')
    ax.set_ylabel('Expextation values')
    ax.legend(('Sigma1-z','Sigma2_z'))
    plt.show()


def createFilename(path,index):
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
    """Animates two-level state vectors on the Bloch sphere"""
    states = result.states
    createAnimationVideo(states, dirc)