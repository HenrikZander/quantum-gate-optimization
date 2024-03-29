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
import functions
import simulationManager

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
    resultFolderPath.set("/Users/emilingelsten/Documents/KA-kod/quantum-gate-optimization/Results/5Q4C_minZZ_approx/Solutions") # C:/Gateside/Results
    
    global selectedSimulationOutput
    selectedSimulationOutput = StringVar(root)

    global energyLevels
    global eigenstateIndex
    energyLevels = IntVar(root)
    eigenstateIndex = IntVar(root)


######################################################################################################################################################################
# High level functions.


def simulateSolution(solutionData):
    writeStatus("Simulation started! Please wait.")
    popTransfer, fidelityPlot, stabilityPlot = identifySimulationPlot()
    guiData = getAllVariables()

    guiData["populationTransferPlot"] = popTransfer
    guiData["fidelityPlot"] = fidelityPlot
    guiData["stabilityPlot"] = stabilityPlot

    simulationManager.simulate(solutionData,guiData)
    writeStatus("Simulation complete.")


def identifySimulationPlot():
    if selectedSimulationOutput.get() == "popTrans":
        return (True, False, False)
    elif selectedSimulationOutput.get() == "fidelity":
        return (False, True, False)
    elif selectedSimulationOutput.get() == "stability":
        return (False, False, True)
    else:
        return (False, False, False)


######################################################################################################################################################################
# Functions to interact with the global variables.


def writeStatus(statusString):
    solutionPreviewField.config(state="normal")
    solutionPreviewField.delete("8.0", "end")
    solutionPreviewField.insert("8.0", "\nStatus: " + statusString)
    solutionPreviewField.config(state="disabled")


def invalidSolutionStatus():
    solutionPreviewField.config(state="normal")
    solutionPreviewField.delete("1.0", "end")
    solutionPreviewField.insert("1.0", "Status: Invalid solution! Try again!")
    solutionPreviewField.config(state="disabled")


def writeSolutionParameterStatus(solutionData):
    if solutionData["signalType"] == "arccos":
        solutionPreviewField.config(state="normal")
        solutionPreviewField.delete("1.0", "end")
        solutionPreviewField.insert("1.0", f'Resulting gate: {solutionData["gateType"]}\n')
        solutionPreviewField.insert("2.0", f'Signal type: \u03A6 = 1 / \u03C0 \u2022 arccos[ (A + B(t) \u2022 cos(\u03C9t))\u00B2 ]\n\n')
        solutionPreviewField.insert("4.0", f'Strength of DC-part A [1]: {solutionData["dcAmplitude"]}\n')
        solutionPreviewField.insert("5.0", f'Amplitude of B(t) [1]: {solutionData["acAmplitude"]}\n')
        solutionPreviewField.insert("6.0", f'Frequency \u03C9 of AC-flux [GHz]: {solutionData["omegaPhi"]}\n')
        solutionPreviewField.insert("7.0", f'Total modulation time of AC-flux [ns]: {solutionData["modulationTime"]}\n\n')
        solutionPreviewField.insert("9.0", "Status: Solution format valid! Ready to simulate.")
        solutionPreviewField.config(state="disabled")
    else:
        solutionPreviewField.config(state="normal")
        solutionPreviewField.delete("1.0", "end")
        solutionPreviewField.insert("1.0", f'Resulting gate: {solutionData["gateType"]}\n')
        solutionPreviewField.insert("2.0", f'Signal type: \u03A6 = \u0398 + \u03B4(t) \u2022 cos(\u03C9t)\n\n')
        solutionPreviewField.insert("4.0", f'Strength of DC-flux \u0398 [\u03A6{subscriptZero}]: {solutionData["theta"]}\n')
        solutionPreviewField.insert("5.0", f'Amplitude of \u03B4(t) [\u03A6{subscriptZero}]: {solutionData["delta"]}\n')
        solutionPreviewField.insert("6.0", f'Frequency \u03C9 of AC-flux [GHz]: {solutionData["omegaPhi"]}\n')
        solutionPreviewField.insert("7.0", f'Total modulation time of AC-flux [ns]: {solutionData["modulationTime"]}\n\n')
        solutionPreviewField.insert("9.0", "Status: Solution format valid! Ready to simulate.")
        solutionPreviewField.config(state="disabled")


def getAllVariables():
    data = {}

    data["frequencies"] = (2*np.pi*frequencyQ1.get(), 2*np.pi*frequencyQ2.get(), 2*np.pi*frequencyCoupler.get())
    data["anharmonicities"] = (2*np.pi*anharmonicityQ1.get(), 2*np.pi*anharmonicityQ2.get(), 2*np.pi*anharmonicityCoupler.get())
    data["couplings"] = (2*np.pi*couplingQ1.get(), 2*np.pi*couplingQ2.get())

    data["energy-levels"] = energyLevels.get()
    data["solution-path"] = solutionPath.get()
    data["eigenstateIndex"] = int(eigenstateIndex.get())

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


def enableStartSimulationButton():
    startSimulationButton.config(background="green", command=startSimulation)


def disableStartSimulationButton():
    startSimulationButton.config(background="gray", command=NONE)


######################################################################################################################################################################
# Button callbacks


def toggleInitialStateSelection():
    if selectedSimulationOutput.get() == "popTrans":
        selectEigenstateIndex.config(state=NORMAL)
        selectEnergyLevels.config(state="readonly")

    else:
        selectEigenstateIndex.config(state=DISABLED)
        selectEnergyLevels.config(state=DISABLED)


def startSimulation():
    solutionData = dataManager.getFromjson(solutionPath.get())
    simulateSolution(solutionData)


def loadSolution():
    solutionFilePath = filedialog.askopenfilename(title="Select solution", defaultextension='.json', filetypes=[("JSON files (.json)", '*.json')], initialdir="./Results")
    if solutionFilePath:
        try:
            solutionData = dataManager.getFromjson(solutionFilePath)
            setCircuitVariables(solutionData)
            writeSolutionParameterStatus(solutionData)
            enableStartSimulationButton()
        except:
            print("Exception!")
            disableStartSimulationButton()
            invalidSolutionStatus()
        solutionPath.set(solutionFilePath)


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
    separator.grid(row=8, column=0, columnspan=3, sticky="nesw", pady=(10,0))

    saveFolderTitleFrame = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="red")
    saveFolderTitleFrame.grid(row=9, column=0, columnspan=3)

    saveFolderTitle = Label(saveFolderTitleFrame, text='Configure result output', font=('Helvetica', 12))
    titleFont = tkFont.Font(saveFolderTitle, saveFolderTitle.cget("font"))
    titleFont.configure(underline=True)
    saveFolderTitle.configure(font=titleFont)
    saveFolderTitle.place(anchor='center', relx=0.5, rely=0.5)

    saveFolderFrameOuter = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="orange")
    saveFolderFrameOuter.grid(row=10, column=0, columnspan=3)

    saveFolderFrameInner = Frame(saveFolderFrameOuter)
    saveFolderFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    saveFolderEntry = Entry(saveFolderFrameInner, state="readonly", width=45, readonlybackground="white", textvariable=resultFolderPath)
    saveFolderEntry.pack(side=LEFT)

    saveFolderButton = Button(saveFolderFrameInner, text="Select Folder", command=selectSaveFolder, background="#21e4d7")
    saveFolderButton.pack(side=LEFT)


def generateSimulationOutputSelection(parent):
    separator = ttk.Separator(parent, orient='horizontal')
    separator.grid(row=4, column=0, columnspan=3, sticky="nesw", pady=(10,0))

    selectSimulationOutputTitleFrame = Frame(parent, height=50, width=relativeWidth*width*0.40)  # , background="red")
    selectSimulationOutputTitleFrame.grid(row=5, column=0, columnspan=3)

    selectSimulationOutputTitle = Label(selectSimulationOutputTitleFrame, text='Simulation Output', font=('Helvetica', 12))
    titleFont = tkFont.Font(selectSimulationOutputTitle, selectSimulationOutputTitle.cget("font"))
    titleFont.configure(underline=True)
    selectSimulationOutputTitle.configure(font=titleFont)
    selectSimulationOutputTitle.place(anchor='center', relx=0.5, rely=0.5)

    selectSimulationOutputFrameOuter = Frame(parent, height=55, width=relativeWidth*width*0.40)  # , background="orange")
    selectSimulationOutputFrameOuter.grid(row=6, column=0, columnspan=3)
    selectSimulationOutputFrameOuter.pack_propagate(0)

    selectSimulationOutputFrameInner = LabelFrame(selectSimulationOutputFrameOuter, text="Simulation output: ")
    selectSimulationOutputFrameInner.place(anchor='n', relx=0.5, rely=0)

    selectPopulationTransfer = Radiobutton(selectSimulationOutputFrameInner, text="Population Transfer", value="popTrans", variable=selectedSimulationOutput, command=toggleInitialStateSelection)
    selectPopulationTransfer.select()
    selectPopulationTransfer.pack(side=LEFT)

    selectFidelityPlot = Radiobutton(selectSimulationOutputFrameInner, text="Fidelity Plot", value="fidelity", variable=selectedSimulationOutput, command=toggleInitialStateSelection)
    selectFidelityPlot.pack(side=LEFT)

    selectStabilityPlot = Radiobutton(selectSimulationOutputFrameInner, text="Stability Plot", value="stability", variable=selectedSimulationOutput, command=toggleInitialStateSelection)
    selectStabilityPlot.pack(side=LEFT)

    ######################### Select eigenstate index #########################
    
    selectEigenstateIndexFrameOuter = Frame(parent, height=55, width=relativeWidth*width*0.40)  # , background="blue")
    selectEigenstateIndexFrameOuter.grid(row=7, column=0, columnspan=3)

    selectEigenstateIndexFrameInner = Frame(selectEigenstateIndexFrameOuter)  # , background="orange")
    selectEigenstateIndexFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectEigenstateIndexTitle = Label(selectEigenstateIndexFrameInner, text="Initial eigenstate index (0 = lowest energy state):")
    selectEigenstateIndexTitle.pack(side=LEFT, padx=(0, 5))

    global selectEigenstateIndex
    selectEigenstateIndex = ttk.Spinbox(selectEigenstateIndexFrameInner, from_=0, to=100000, textvariable=eigenstateIndex, width=6)
    selectEigenstateIndex.set(0)
    selectEigenstateIndex.pack(side=LEFT)

    ###########################################################################


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
    generateSimulationOutputSelection(circuitFrame)
    generateSaveFolderInput(circuitFrame, entryCharacterWidth)


def simulateCircuitParameterFrame(topFrame):
    circuitFrame = Frame(topFrame, height=relativeHeight * height, width=relativeWidth*width*0.40)
    circuitFrame.grid(column=0, row=0)
    circuitFrame.grid_propagate(0)

    generateCircuitInputs(circuitFrame)


##################################################################################
# Functions that generate the controls for simulation configuration and widgets that show of a selected solution.


def generateSolutionPreviewField(solutionFrame):
    solutionPreviewFrameOuter = Frame(solutionFrame, height=relativeHeight*height-210, width=relativeWidth*width*0.60)# , background="green")
    solutionPreviewFrameOuter.grid(row=5, column=0, columnspan=2)
    solutionPreviewFrameOuter.grid_propagate(0)

    solutionPreviewFrameInner = Frame(solutionPreviewFrameOuter)
    solutionPreviewFrameInner.place(anchor='n', relx=0.5, rely=0)

    global solutionPreviewField
    solutionPreviewField = Text(solutionPreviewFrameInner, height=17, width=65, state="disabled", background="lightgray", font=("Helvetica", 11), spacing1=5)
    solutionPreviewField.pack(anchor="center")


def generateStartSimulationFrame(settingsFrame):
    startSimulationFrameOuter = Frame(settingsFrame, height=55, width=relativeWidth*width*0.30)# , background="blue")
    startSimulationFrameOuter.grid(row=1, column=0)

    startSimulationFrameInner = Frame(startSimulationFrameOuter)# , background="orange")
    startSimulationFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    global startSimulationButton
    startSimulationButton = Button(startSimulationFrameInner, text="Start Simulation", command=NONE, background="gray")
    startSimulationButton.pack(side=LEFT)


def generateSelectEnergyLevels(settingsFrame):
    selectEnergyLevelsFrameOuter = Frame(settingsFrame, height=55, width=relativeWidth*width*0.30)  # , background="blue")
    selectEnergyLevelsFrameOuter.grid(row=1, column=0)

    selectEnergyLevelsFrameInner = Frame(selectEnergyLevelsFrameOuter)  # , background="orange")
    selectEnergyLevelsFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectEnergyLevelsTitle = Label(selectEnergyLevelsFrameInner, text="Number of energy levels per qubit:")
    selectEnergyLevelsTitle.pack(side=LEFT, padx=(0, 5))

    global selectEnergyLevels
    selectEnergyLevels = ttk.Spinbox(selectEnergyLevelsFrameInner, from_=2, to=7, textvariable=energyLevels, width=4, state="readonly")
    selectEnergyLevels.set(3)
    selectEnergyLevels.pack(side=LEFT)


def generateLeftSettings(settingsFrameLeft):
    generateSelectEnergyLevels(settingsFrameLeft)


def generateRightSettings(settingsFrameRight):
    generateStartSimulationFrame(settingsFrameRight)


def generateSelectSolutionInput(solutionFrame):
    selectSolutionFrameOuter = Frame(solutionFrame, height=50, width=relativeWidth*width*0.60)# , background="green")
    selectSolutionFrameOuter.grid(column=0, row=1, columnspan=2)
    selectSolutionFrameOuter.grid_propagate(0)

    selectSolutionFrameInner = Frame(selectSolutionFrameOuter)
    selectSolutionFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    selectSolutionEntry = Entry(selectSolutionFrameInner, state="readonly", width=65, readonlybackground="white", textvariable=solutionPath)
    selectSolutionEntry.pack(side=LEFT)

    selectSolutionButton = Button(selectSolutionFrameInner, text="Select Solution", command=loadSolution, background="#21e4d7")
    selectSolutionButton.pack(side=LEFT)


def solutionSelectionFrame(topFrame):
    selectSolutionTitleFrameHeight = 50
    selectSolutionFrameHeight = 110

    solutionFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.60)# , background="yellow")
    solutionFrame.grid(row=0, column=2)
    solutionFrame.grid_propagate(0)

    selectSolutionTitleFrame = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.60)# , background="orange")
    selectSolutionTitleFrame.grid(row=0, column=0, columnspan=2)

    settingsTitle = Label(selectSolutionTitleFrame, text='Simulation Configuration', font=('Helvetica', 12))
    titleSettingsFont = tkFont.Font(settingsTitle, settingsTitle.cget("font"))
    titleSettingsFont.configure(underline=True)
    settingsTitle.configure(font=titleSettingsFont)
    settingsTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateSelectSolutionInput(solutionFrame)

    solutionSettingsFrameLeft = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.30) #, background="blue")
    solutionSettingsFrameLeft.grid(column=0, row=2)
    solutionSettingsFrameLeft.grid_propagate(0)

    solutionSettingsFrameRight = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.30) #, background="red")
    solutionSettingsFrameRight.grid(column=1, row=2)
    solutionSettingsFrameRight.grid_propagate(0)

    separator = ttk.Separator(solutionFrame, orient='horizontal')
    separator.grid(row=3, column=0, columnspan=2, sticky="nesw")

    solutionPreviewTitleFrame = Frame(solutionFrame, height=selectSolutionTitleFrameHeight, width=relativeWidth*width*0.60) #, background="orange")
    solutionPreviewTitleFrame.grid(row=4, column=0, columnspan=2)

    solutionPreviewTitle = Label(solutionPreviewTitleFrame, text='Current Simulation Status', font=('Helvetica', 12))
    solutionPreviewTitle.configure(font=titleSettingsFont)
    solutionPreviewTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateLeftSettings(solutionSettingsFrameLeft)
    generateRightSettings(solutionSettingsFrameRight)
    generateSolutionPreviewField(solutionFrame)


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