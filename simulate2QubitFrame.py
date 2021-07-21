######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: simluate2QubitWindow.py

# Author(s): Henrik Zander

# Date created: 20 July 2021

# Copyright 2021, Henrik Zander, All rights reserved.

######################################################################################################################################################################

from tkinter import *
from tkinter import filedialog
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time
from threading import Thread
import numpy as np

import dataManager
#import simulateManager

######################################################################################################################################################################
# Global variables and function to initiate global variables

relativeHeight = 1
relativeWidth = 1

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
subscriptZero = "0".translate(SUB)


def initiateGlobalVariables(root, givenHeight):
    global height
    global width
    height = givenHeight
    width = int(1.62*height)

    global frequencyQ1
    global frequencyQ2
    global frequencyCoupler
    frequencyQ1 = DoubleVar(root)
    frequencyQ2 = DoubleVar(root)
    frequencyCoupler = DoubleVar(root)

    global anharmonicityQ1
    global anharmonicityQ2
    global anharmonicityCoupler
    anharmonicityQ1 = DoubleVar(root)
    anharmonicityQ2 = DoubleVar(root)
    anharmonicityCoupler = DoubleVar(root)

    global couplingQ1
    global couplingQ2
    couplingQ1 = DoubleVar(root)
    couplingQ2 = DoubleVar(root)

    global solutionPath
    global resultFolderPath
    solutionPath = StringVar(root)
    resultFolderPath = StringVar(root)
    resultFolderPath.set("C:/Gateside/Results")

    global energyLevels
    energyLevels = IntVar(root)


######################################################################################################################################################################
# High level functions.

def simulateSolution():
    pass

######################################################################################################################################################################
# Functions to interact with the global variables.


def getAllVariables():
    data = {}

    data["frequencies"] = (2*np.pi*frequencyQ1.get(), 2 * np.pi*frequencyQ2.get(), 2*np.pi*frequencyCoupler.get())

    data["anharmonicities"] = (2*np.pi*anharmonicityQ1.get(), 2 * np.pi*anharmonicityQ2.get(), 2*np.pi*anharmonicityCoupler.get())

    data["couplings"] = (2*np.pi*couplingQ1.get(), 2*np.pi*couplingQ2.get())

    data["save-folder"] = resultFolderPath.get()

    data["runDE"] = runDifferentialEvolution.get()
    data["runSHG"] = runSimplicalHomologyGlobal.get()
    data["runDA"] = runDualAnneling.get()

    data[selectedGate.get()] = True
    data["energy-levels"] = energyLevels.get()

    data["rise-time"] = riseTime.get()

    data["theta"] = (thetaLower.get(), thetaUpper.get())
    data["delta"] = (deltaLower.get(), deltaUpper.get())
    data["omega-phi"] = (2*np.pi*omegaPhiLower.get(), 2*np.pi*omegaPhiUpper.get())
    data["modulation-time"] = (modulationTimeLower.get(), modulationTimeUpper.get())

    return data


def setCircuitVariables(data):
    frequencyQ1.set(data["frequencies"][0])
    frequencyQ2.set(data["frequencies"][1])
    frequencyCoupler.set(data["frequencies"][2])

    anharmonicityQ1.set(data["anharmonicities"][0])
    anharmonicityQ2.set(data["anharmonicities"][1])
    anharmonicityCoupler.set(data["anharmonicities"][2])

    couplingQ1.set(data["couplings"][0])
    couplingQ2.set(data["couplings"][1])


def getCircuitVariables():
    data = {}

    data["frequencies"] = [frequencyQ1.get(), frequencyQ2.get(), frequencyCoupler.get()]
    data["anharmonicities"] = [anharmonicityQ1.get(), anharmonicityQ2.get(), anharmonicityCoupler.get()]
    data["couplings"] = [couplingQ1.get(), couplingQ2.get()]

    return data


def setDefaultCircuitValuesFromVariables(configData, newCircuitData):
    configData["frequencies"] = newCircuitData["frequencies"]
    configData["anharmonicities"] = newCircuitData["anharmonicities"]
    configData["couplings"] = newCircuitData["couplings"]

    return configData


######################################################################################################################################################################
# Button callbacks


def loadSolution():
    solutionFile = filedialog.askopenfilename(title="Select solution", defaultextension='.json', filetypes=[("JSON files (.json)", '*.json')])
    if solutionFile:
        solutionData = dataManager.getFromjson(solutionFile)
        setCircuitVariables(solutionData)


def selectSaveFolder():
    folderPath = filedialog.askdirectory(title="Choose folder to save results in", initialdir="C:/")
    if folderPath:
        resultFolderPath.set(folderPath)


def loadCircuit():
    circuitFile = filedialog.askopenfilename(title="Select circuit", defaultextension='.json', filetypes=[("JSON files (.json)", '*.json')])
    if circuitFile:
        circuitData = dataManager.getFromjson(circuitFile)
        setCircuitVariables(circuitData)


def changeDefaultCircuit():
    configData = dataManager.getFromjson("config.json")
    newCircuitData = getCircuitVariables()
    setDefaultCircuitValuesFromVariables(configData, newCircuitData)
    dataManager.dumpTojson(configData, "config.json")


def useDefaultCircuit():
    configData = dataManager.getFromjson("config.json")
    setCircuitVariables(configData)


def exportCircuit():
    circuitData = getCircuitVariables()
    path = filedialog.asksaveasfilename(title="Save circuit", defaultextension='.json', filetypes=[("JSON files (.json)", '*.json')])
    if path:
        dataManager.dumpTojson(circuitData, path)

######################################################################################################################################################################
# Functions that generate the window for simulation control of two qubit gates.

##################################################################################
# Functions that generate the circuit parameter frames and result saving functionality.


def generateSaveFolderInput(circuitFrame, entryCharacterWidth):
    separator = ttk.Separator(circuitFrame, orient='horizontal')
    separator.grid(row=5, column=0, columnspan=3, sticky="nesw", pady=(10,0))

    saveFolderTitleFrame = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="red")
    saveFolderTitleFrame.grid(row=6, column=0, columnspan=3)

    saveFolderTitle = Label(saveFolderTitleFrame, text='Configure result output', font=('Helvetica', 12))
    titleFont = tkFont.Font(saveFolderTitle, saveFolderTitle.cget("font"))
    titleFont.configure(underline=True)
    saveFolderTitle.configure(font=titleFont)
    saveFolderTitle.place(anchor='center', relx=0.5, rely=0.5)

    saveFolderFrameOuter = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="orange")
    saveFolderFrameOuter.grid(row=7, column=0, columnspan=3)

    saveFolderFrameInner = Frame(saveFolderFrameOuter)
    saveFolderFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    saveFolderEntry = Entry(saveFolderFrameInner, state="readonly", width=45, readonlybackground="white", textvariable=resultFolderPath)
    saveFolderEntry.pack(side=LEFT)

    saveFolderButton = Button(saveFolderFrameInner, text="Select Folder", command=selectSaveFolder, background="#21e4d7")
    saveFolderButton.pack(side=LEFT)


def generateCouplingInputs(circuitFrame, entryCharacterWidth):
    couplingInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    couplingInputFrameOuter.grid(row=3, column=0, columnspan=3)

    couplingInputFrameInner = LabelFrame(couplingInputFrameOuter, text="Coupling strength to coupler [GHz]: ")
    couplingInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    couplingLabelQ1 = Label(couplingInputFrameInner, text="Q1: ")
    couplingLabelQ1.pack(side=LEFT, padx=(0, 0))

    couplingEntryQ1 = Entry(couplingInputFrameInner, width=entryCharacterWidth, textvariable=couplingQ1, state="readonly")
    couplingEntryQ1.pack(side=LEFT)

    couplingLabelQ2 = Label(couplingInputFrameInner, text="Q2: ")
    couplingLabelQ2.pack(side=LEFT, padx=(10, 0), pady=5)

    couplingEntryQ2 = Entry(couplingInputFrameInner, width=entryCharacterWidth, textvariable=couplingQ2, state="readonly")
    couplingEntryQ2.pack(side=LEFT, padx=(0, 5))


def generateAnharmonicityInputs(circuitFrame, entryCharacterWidth):
    anharmonicityInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    anharmonicityInputFrameOuter.grid(row=2, column=0, columnspan=3)

    anharmonicityInputFrameInner = LabelFrame(anharmonicityInputFrameOuter, text="Anharmonicites [GHz]: ")
    anharmonicityInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    anharmonicityLabelQ1 = Label(anharmonicityInputFrameInner, text="Q1: ")
    anharmonicityLabelQ1.pack(side=LEFT, padx=(0, 0))

    anharmonicityEntryQ1 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityQ1, state="readonly")
    anharmonicityEntryQ1.pack(side=LEFT)

    anharmonicityLabelQ2 = Label(anharmonicityInputFrameInner, text="Q2: ")
    anharmonicityLabelQ2.pack(side=LEFT, padx=(10, 0))

    anharmonicityEntryQ2 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityQ2, state="readonly")
    anharmonicityEntryQ2.pack(side=LEFT)

    anharmonicityLabelCoupler = Label(anharmonicityInputFrameInner, text="Coupler: ")
    anharmonicityLabelCoupler.pack(side=LEFT, padx=(10, 0), pady=5)

    anharmonicityEntryCoupler = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityCoupler, state="readonly")
    anharmonicityEntryCoupler.pack(side=LEFT, padx=(0, 5))


def generateFrequencyInputs(circuitFrame, entryCharacterWidth):
    frequencyInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    frequencyInputFrameOuter.grid(row=1, column=0, columnspan=3)

    frequencyInputFrameInner = LabelFrame(frequencyInputFrameOuter, text="Frequencies [GHz]: ")
    frequencyInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    frequencyLabelQ1 = Label(frequencyInputFrameInner, text="Q1: ")
    frequencyLabelQ1.pack(side=LEFT, padx=(0, 0))

    frequencyEntryQ1 = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyQ1, state="readonly")
    frequencyEntryQ1.pack(side=LEFT)

    frequencyLabelQ2 = Label(frequencyInputFrameInner, text="Q2: ")
    frequencyLabelQ2.pack(side=LEFT, padx=(10, 0))

    frequencyEntryQ2 = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyQ2, state="readonly")
    frequencyEntryQ2.pack(side=LEFT)

    frequencyLabelCoupler = Label(frequencyInputFrameInner, text="Coupler: ")
    frequencyLabelCoupler.pack(side=LEFT, padx=(10, 0), pady=5)

    frequencyEntryCoupler = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyCoupler, state="readonly")
    frequencyEntryCoupler.pack(side=LEFT, padx=(0, 5))


def generateCircuitInputControls(circuitFrame):
    controlsInputFrameOuter = Frame(circuitFrame, height=80, width=relativeWidth*width*0.40)  # , background="orange")
    controlsInputFrameOuter.grid(row=4, column=0, columnspan=3)

    controlsInputFrameInner = Frame(controlsInputFrameOuter)
    controlsInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    loadCircuitButton = Button(controlsInputFrameInner, text="Load Circuit", command=loadCircuit, padx=3, pady=3, background="#21e4d7")
    loadCircuitButton.grid(row=0, column=0, padx=4, pady=4)

    changeDefaultCircuitButton = Button(controlsInputFrameInner, text="Set to Default", command=changeDefaultCircuit, padx=3, pady=3, background="#21e4d7")
    changeDefaultCircuitButton.grid(row=0, column=1, padx=4, pady=4)

    useDefaultCircuitButton = Button(controlsInputFrameInner, text="Use Default", command=useDefaultCircuit, padx=3, pady=3, background="#21e4d7")
    useDefaultCircuitButton.grid(row=0, column=2, padx=4, pady=4)

    generateCircuitButton = Button(controlsInputFrameInner, text="Export", command=exportCircuit, padx=3, pady=3, background="#21e4d7")
    generateCircuitButton.grid(row=0, column=3, padx=4, pady=4)


def generateCircuitInputs(circuitFrame):
    entryCharacterWidth = 10

    circuitTitleFrame = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="red")
    circuitTitleFrame.grid(row=0, column=0, columnspan=3)

    circuitTitle = Label(circuitTitleFrame, text='Circuit parameters', font=('Helvetica', 12))
    titleFont = tkFont.Font(circuitTitle, circuitTitle.cget("font"))
    titleFont.configure(underline=True)
    circuitTitle.configure(font=titleFont)
    circuitTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateFrequencyInputs(circuitFrame, entryCharacterWidth)
    generateAnharmonicityInputs(circuitFrame, entryCharacterWidth)
    generateCouplingInputs(circuitFrame, entryCharacterWidth)
    generateSaveFolderInput(circuitFrame, entryCharacterWidth)


def simulateCircuitParameterFrame(topFrame):
    circuitFrame = Frame(topFrame, height=relativeHeight * height, width=relativeWidth*width*0.40)
    circuitFrame.grid(column=0, row=0)
    circuitFrame.grid_propagate(0)

    generateCircuitInputs(circuitFrame)
    #generateCircuitInputControls(circuitFrame)

##################################################################################
# Functions that generate the controls for simulation configuration and widgets that show of a selected solution.


def generateSelectEnergyLevels(settingsFrameRight):
    selectEnergyLevelsFrameOuter = Frame(settingsFrameRight, height=55, width=relativeWidth*width*0.30)  # , background="blue")
    selectEnergyLevelsFrameOuter.grid(row=1, column=0)

    selectEnergyLevelsFrameInner = Frame(selectEnergyLevelsFrameOuter)  # , background="orange")
    selectEnergyLevelsFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectEnergyLevelsTitle = Label(selectEnergyLevelsFrameInner, text="Number of energy levels per qubit:")
    selectEnergyLevelsTitle.pack(side=LEFT, padx=(0, 5))

    selectEnergyLevels = ttk.Spinbox(selectEnergyLevelsFrameInner, from_=2, to=7, textvariable=energyLevels, width=4, state="readonly")
    selectEnergyLevels.set(3)
    selectEnergyLevels.pack(side=LEFT)


def generateRightSettings(settingsFrameRight):
    generateSelectEnergyLevels(settingsFrameRight)


def generateSelectSolutionInput(solutionFrame):
    selectSolutionFrameOuter = Frame(solutionFrame, height=50, width=relativeWidth*width*0.60, background="green")
    selectSolutionFrameOuter.grid(column=0, row=1, columnspan=2)
    selectSolutionFrameOuter.grid_propagate(0)

    selectSolutionFrameInner = Frame(selectSolutionFrameOuter)
    selectSolutionFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    selectSolutionEntry = Entry(selectSolutionFrameInner, state="readonly", width=50, readonlybackground="white", textvariable=solutionPath)
    selectSolutionEntry.pack(side=LEFT)

    selectSolutionButton = Button(selectSolutionFrameInner, text="Select Solution", command=loadSolution, background="#21e4d7")
    selectSolutionButton.pack(side=LEFT)


def solutionSelectionFrame(topFrame):
    selectSolutionTitleFrameHeight = 50
    selectSolutionFrameHeight = 110

    solutionFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.60, background="yellow")
    solutionFrame.grid(row=0, column=2)
    solutionFrame.grid_propagate(0)

    selectSolutionTitleFrame = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.60, background="orange")
    selectSolutionTitleFrame.grid(row=0, column=0, columnspan=2)

    settingsTitle = Label(selectSolutionTitleFrame, text='Simulation Configuration', font=('Helvetica', 12))
    titleSettingsFont = tkFont.Font(settingsTitle, settingsTitle.cget("font"))
    titleSettingsFont.configure(underline=True)
    settingsTitle.configure(font=titleSettingsFont)
    settingsTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateSelectSolutionInput(solutionFrame)

    solutionSettingsFrameLeft = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.30, background="blue")
    solutionSettingsFrameLeft.grid(column=0, row=2)
    solutionSettingsFrameLeft.grid_propagate(0)

    solutionSettingsFrameRight = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.30, background="red")
    solutionSettingsFrameRight.grid(column=1, row=2)
    solutionSettingsFrameRight.grid_propagate(0)

    separator = ttk.Separator(solutionFrame, orient='horizontal')
    separator.grid(row=3, column=0, columnspan=2, sticky="nesw")

    solutionPreviewTitleFrame = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.60, background="orange")
    solutionPreviewTitleFrame.grid(row=4, column=0, columnspan=2)

    solutionPreviewTitle = Label(solutionPreviewTitleFrame, text='Parameters of selected solution', font=('Helvetica', 12))
    solutionPreviewTitle.configure(font=titleSettingsFont)
    solutionPreviewTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateRightSettings(solutionSettingsFrameRight)


##################################################################################

######################################################################################################################################################################
# Functions that generate different windows


def simulateControlFrame(root, parentWidget, givenHeight, givenWidth):
    initiateGlobalVariables(root, givenHeight)

    global relativeHeight
    global relativeWidth
    global height
    global width

    height = givenHeight
    width = givenWidth

    windowFrame = Frame(parentWidget, height=height, width=width)
    windowFrame.place(anchor="center", relx=0.5, rely=0.5)
    windowFrame.grid_propagate(0)

    topFrame = Frame(windowFrame, height=relativeHeight*height, width=width)
    topFrame.grid(row=0, column=0)
    topFrame.grid_propagate(0)

    separator = ttk.Separator(topFrame, orient='vertical')
    separator.grid(row=0, column=1, sticky="nesw")

    simulateCircuitParameterFrame(topFrame)
    solutionSelectionFrame(topFrame)

    return windowFrame


######################################################################################################################################################################