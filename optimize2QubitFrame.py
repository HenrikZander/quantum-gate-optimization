######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: optimize2QubitFrame.py

# Author(s): Henrik Zander

# Date created: 12 July 2021

# Copyright 2021, Henrik Zander, All rights reserved.

######################################################################################################################################################################

from tkinter import *
from tkinter import filedialog
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time
from threading import Thread
import numpy as np
import os

import dataManager
import optimizeManager
import boundaryConditionPopUpWindow as boundaryWindow

######################################################################################################################################################################
# Global variables and function to initiate global variables

relativeHeight = 0.85 #0.8
relativeWidth = 1

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
subscriptZero = "0".translate(SUB)

runOptimizer = False


def initiateGlobalVariables(rootWindow, givenHeight):
    global root
    global height
    global width
    height = givenHeight
    width = int(1.62*height)
    root = rootWindow

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

    global resultFolderPath
    global numberOfConsecutiveSessions
    resultFolderPath = StringVar(root)
    resultFolderPath.set(os.path.abspath("") + "/Results")
    numberOfConsecutiveSessions = IntVar(root)

    global runDifferentialEvolution
    global runSimplicalHomologyGlobal
    global runDualAnneling
    runDifferentialEvolution = BooleanVar(root)
    runSimplicalHomologyGlobal = BooleanVar(root)
    runDualAnneling = BooleanVar(root)

    global selectedGate
    global energyLevels
    selectedGate = StringVar(root)
    energyLevels = IntVar(root)

    global signalShape
    global riseTime
    signalShape = StringVar(root)
    riseTime = IntVar(root)

    global x0Lower
    global x0Upper
    x0Lower = DoubleVar(root)
    x0Upper = DoubleVar(root)

    global x1Lower
    global x1Upper
    x1Lower = DoubleVar(root)
    x1Upper = DoubleVar(root)

    global omegaPhiLower
    global omegaPhiUpper
    omegaPhiLower = DoubleVar(root)
    omegaPhiUpper = DoubleVar(root)

    global modulationTimeLower
    global modulationTimeUpper
    modulationTimeLower = DoubleVar(root)
    modulationTimeUpper = DoubleVar(root)

    global progressValue
    global statusString
    progressValue = DoubleVar(root)
    statusString = StringVar(root)

    global x0LabelVar
    global x1LabelVar
    x0LabelVar = StringVar(root)
    x1LabelVar = StringVar(root)


######################################################################################################################################################################
# High level functions.


def enableStopButton():
    startOptimizeButton.config(background="grey", command=NONE)
    stopOptimizeButton.config(background="#FF0000", command=stopOptimizing)


def resetStartStopButtons():
    startOptimizeButton.config(background="#00FF00", command=startOptimizing)
    stopOptimizeButton.config(background="grey", command=NONE)


def disableStartStopButtons():
    startOptimizeButton.config(background="grey", command=NONE)
    stopOptimizeButton.config(background="grey", command=NONE)


def processStopped():
    statusString.set("Status: Optimization stopped. Ready to optimize.")
    progressValue.set(0)
    resetStartStopButtons()


def processFinished():
    statusString.set("Status: Optimization finished. Ready to optimize.")
    progressValue.set(0)

'''
def identifyGate(userData):
    try:
        iSWAP = userData["iSWAP"]
        SWAP = False
        CZ = False
    except:
        iSWAP = False
        try:
            SWAP = userData["SWAP"]
            CZ = False
        except:
            SWAP = False
            CZ = True

    return (iSWAP, SWAP, CZ)
'''

def callOptimizeGate():
    dataFromUser = getAllVariablesForTheOptimizer()

    optimizeManager.optimize2QubitGate(userData=dataFromUser)


def scheduleOptimizingSessions():
    global process
    global progressValue
    global runOptimizer
    global numberOfConsecutiveSessions
    sessions = numberOfConsecutiveSessions.get()

    for i in range(sessions):
        if not runOptimizer:
            break
        progressValue.set(0)
        statusString.set(f'Status: Optimization ({i+1}/{sessions}) running, please wait.')
        ##########################################
        # This decides which function to run
        process = Thread(target=callOptimizeGate)
        process.start()
        process.join()
        ##########################################
    runOptimizer = False
    resetStartStopButtons()


######################################################################################################################################################################
# Functions to interact with the global variables.


def getX0X1names():
    if signalShape.get() == '\u03A6 = 1 / \u03C0 \u2022 arccos[ (A + B(t) \u2022 cos(\u03C9t))\u00B2 ]':
        x0name = "dcAmplitude"
        x1name = "acAmplitude"
    elif signalShape.get() == '\u03A6 = \u0398 + \u03B4(t) \u2022 cos(\u03C9t)':
        x0name = "theta"
        x1name = "delta"
    
    return x0name, x1name


def getSignalType():
    if signalShape.get() == '\u03A6 = 1 / \u03C0 \u2022 arccos[ (A + B(t) \u2022 cos(\u03C9t))\u00B2 ]':
        signalType = 'arccos'
    elif signalShape.get() == '\u03A6 = \u0398 + \u03B4(t) \u2022 cos(\u03C9t)':
        signalType = 'cos'
    
    return signalType


def setBoundaryValues(boundaryValues):
    x0Lower.set(boundaryValues[2][0])
    x0Upper.set(boundaryValues[2][1])

    x1Lower.set(boundaryValues[3][0])
    x1Upper.set(boundaryValues[3][1])

    omegaPhiLower.set(boundaryValues[4][0])
    omegaPhiUpper.set(boundaryValues[4][1])

    modulationTimeLower.set(boundaryValues[5][0])
    modulationTimeUpper.set(boundaryValues[5][1])


def getAllVariablesForTheOptimizer():
    data = {}

    data["frequencies"] = (2*np.pi*frequencyQ1.get(), 2*np.pi*frequencyQ2.get(), 2*np.pi*frequencyCoupler.get())

    data["anharmonicities"] = (2*np.pi*anharmonicityQ1.get(), 2*np.pi*anharmonicityQ2.get(), 2*np.pi*anharmonicityCoupler.get())

    data["couplings"] = (2*np.pi*couplingQ1.get(), 2*np.pi*couplingQ2.get())

    data["save-folder"] = resultFolderPath.get()

    data["runDE"] = runDifferentialEvolution.get()
    data["runSHG"] = runSimplicalHomologyGlobal.get()
    data["runDA"] = runDualAnneling.get()

    data['signalType'] = getSignalType()
    data["gateType"] = selectedGate.get()
    data["energy-levels"] = energyLevels.get()

    data["rise-time"] = riseTime.get()

    x0name, x1name = getX0X1names()

    data[x0name] = (x0Lower.get(), x0Upper.get())
    data[x1name] = (x1Lower.get(), x1Upper.get())
    data["omegaPhi"] = (omegaPhiLower.get(), omegaPhiUpper.get())
    data["modulationTime"] = (modulationTimeLower.get(), modulationTimeUpper.get())

    return data


def getRunOptimizer():
    global runOptimizer
    return runOptimizer


def setProgressValue(value):
    progressValue.set(value)


def setStatus(currentStatus):
    global statusString
    statusString.set(currentStatus)

    
def useDefaultBoundaryValues():
    data = boundaryWindow.getPreset(0)
    setBoundaryValues(data)


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


def setDefaultBoundaryValues():
    newValues = [[x0Lower.get(), x0Upper.get()], [x1Lower.get(), x1Upper.get()], [omegaPhiLower.get(), omegaPhiUpper.get()], [modulationTimeLower.get(), modulationTimeUpper.get()]]
    boundaryWindow.changePreset(0, newValues)


def updateX0X1Labels(event):
    signalType = getSignalType()
    if signalType == 'arccos':
        x0LabelVar.set("DC signal strength A [1]:")
        x1LabelVar.set("Amplitude of B(t) [1]:")
    elif signalType == 'cos':
        x0LabelVar.set("Strength of DC-flux \u0398 [\u03A6"+subscriptZero+"]:")
        x1LabelVar.set("Amplitude of \u03B4(t) [\u03A6"+subscriptZero+"]:")

######################################################################################################################################################################
# Button callbacks


def presetBoundaryConditions():
    boundaryWindow.updateStaticPresets()
    boundaryWindow.selectPresetWindow(root, selectedGate.get())


def setBoundaryDefault():
    setDefaultBoundaryValues()


def useBoundaryDefault():
    useDefaultBoundaryValues()


def selectSaveFolder():
    folderPath = filedialog.askdirectory(title="Choose folder to save results in", initialdir=resultFolderPath.get())
    if folderPath:
        resultFolderPath.set(folderPath)


def startOptimizing():
    global runOptimizer
    runOptimizer = True

    enableStopButton()
    statusString.set("Status: Optimization starting, please wait.")

    optimizerHandler = Thread(target=scheduleOptimizingSessions)
    optimizerHandler.start()


def stopOptimizing():
    global runOptimizer
    runOptimizer = False

    disableStartStopButtons()
    statusString.set("Status: Optimization stopping, please wait.")


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
# Functions that generate the window for optimization control of two qubit gates.

##################################################################################
# Functions that generate the status frame.


def optimizeStatusFrame(bottomFrame):
    global startOptimizeButton
    startOptimizeButton = Button(bottomFrame, text="Start Optimizing", command=startOptimizing, padx=3, pady=3, background="#00FF00", relief=FLAT)
    startOptimizeButton.grid(column=2, row=0, sticky='s')

    global stopOptimizeButton
    stopOptimizeButton = Button(bottomFrame, text="Stop Optimizing", padx=3, pady=3, background="grey", relief=FLAT)
    stopOptimizeButton.grid(column=2, row=1, sticky='n')

    progressFrame = Frame(bottomFrame, height=(1-relativeHeight)*height, width=width*0.80)
    progressFrame.grid(column=0, row=0, rowspan=2, columnspan=2)
    progressFrame.grid_propagate(0)

    progressLabel = Label(bottomFrame, font=("Helvetica", 11), textvariable=statusString)
    statusString.set("Status: Ready to optimize.")
    progressLabel.place(relx=0.05, rely=0.5, anchor='sw')

    global progress
    progress = ttk.Progressbar(bottomFrame, length=width*0.67, orient=HORIZONTAL, mode='determinate', value=0, variable=progressValue, maximum=1)
    progress.place(relx=0.05, rely=0.5, anchor='nw')


##################################################################################
# Functions that generate the circuit parameter frames and result saving functionality.


def generateSaveFolderInput(circuitFrame, entryCharacterWidth):
    separator = ttk.Separator(circuitFrame, orient='horizontal')
    separator.grid(row=5, column=0, columnspan=3, sticky="nesw")

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

    saveFolderButton = Button(saveFolderFrameInner, text="Select Folder", command=selectSaveFolder, background="#21e4d7", relief=FLAT)
    saveFolderButton.pack(side=LEFT)

    numOfSessionsFrameOuter = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)  # , background="red")
    numOfSessionsFrameOuter.grid(row=8, column=0, columnspan=3)

    numOfSessionsFrameInner = Frame(numOfSessionsFrameOuter)
    numOfSessionsFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    numOfSessionsLabel = Label(numOfSessionsFrameInner, text="Number of consecutive optimizations to run: ")
    numOfSessionsLabel.pack(side=LEFT)

    numOfSessionsInputBox = ttk.Spinbox(numOfSessionsFrameInner, from_=1, to=10, state="readonly", width=4, textvariable=numberOfConsecutiveSessions)
    numOfSessionsInputBox.set(1)
    numOfSessionsInputBox.pack(side=LEFT)


def generateCouplingInputs(circuitFrame, entryCharacterWidth):
    couplingInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    couplingInputFrameOuter.grid(row=3, column=0, columnspan=3)

    couplingInputFrameInner = LabelFrame(couplingInputFrameOuter, text="Coupling strength to coupler [GHz]: ")
    couplingInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    couplingLabelQ1 = Label(couplingInputFrameInner, text="Q1: ")
    couplingLabelQ1.pack(side=LEFT, padx=(0, 0))

    couplingEntryQ1 = Entry(couplingInputFrameInner, width=entryCharacterWidth, textvariable=couplingQ1)
    couplingEntryQ1.pack(side=LEFT)

    couplingLabelQ2 = Label(couplingInputFrameInner, text="Q2: ")
    couplingLabelQ2.pack(side=LEFT, padx=(10, 0), pady=5)

    couplingEntryQ2 = Entry(couplingInputFrameInner, width=entryCharacterWidth, textvariable=couplingQ2)
    couplingEntryQ2.pack(side=LEFT, padx=(0, 5))


def generateAnharmonicityInputs(circuitFrame, entryCharacterWidth):
    anharmonicityInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    anharmonicityInputFrameOuter.grid(row=2, column=0, columnspan=3)

    anharmonicityInputFrameInner = LabelFrame(anharmonicityInputFrameOuter, text="Anharmonicites [GHz]: ")
    anharmonicityInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    anharmonicityLabelQ1 = Label(anharmonicityInputFrameInner, text="Q1: ")
    anharmonicityLabelQ1.pack(side=LEFT, padx=(0, 0))

    anharmonicityEntryQ1 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityQ1)
    anharmonicityEntryQ1.pack(side=LEFT)

    anharmonicityLabelQ2 = Label(anharmonicityInputFrameInner, text="Q2: ")
    anharmonicityLabelQ2.pack(side=LEFT, padx=(10, 0))

    anharmonicityEntryQ2 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityQ2)
    anharmonicityEntryQ2.pack(side=LEFT)

    anharmonicityLabelCoupler = Label(anharmonicityInputFrameInner, text="Coupler: ")
    anharmonicityLabelCoupler.pack(side=LEFT, padx=(10, 0), pady=5)

    anharmonicityEntryCoupler = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth, textvariable=anharmonicityCoupler)
    anharmonicityEntryCoupler.pack(side=LEFT, padx=(0, 5))


def generateFrequencyInputs(circuitFrame, entryCharacterWidth):
    frequencyInputFrameOuter = Frame(circuitFrame, height=60, width=relativeWidth*width*0.40)  # , background="orange")
    frequencyInputFrameOuter.grid(row=1, column=0, columnspan=3)

    frequencyInputFrameInner = LabelFrame(frequencyInputFrameOuter, text="Frequencies [GHz]: ")
    frequencyInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    frequencyLabelQ1 = Label(frequencyInputFrameInner, text="Q1: ")
    frequencyLabelQ1.pack(side=LEFT, padx=(0, 0))

    frequencyEntryQ1 = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyQ1)
    frequencyEntryQ1.pack(side=LEFT)

    frequencyLabelQ2 = Label(frequencyInputFrameInner, text="Q2: ")
    frequencyLabelQ2.pack(side=LEFT, padx=(10, 0))

    frequencyEntryQ2 = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyQ2)
    frequencyEntryQ2.pack(side=LEFT)

    frequencyLabelCoupler = Label(frequencyInputFrameInner, text="Coupler: ")
    frequencyLabelCoupler.pack(side=LEFT, padx=(10, 0), pady=5)

    frequencyEntryCoupler = Entry(frequencyInputFrameInner, width=entryCharacterWidth, textvariable=frequencyCoupler)
    frequencyEntryCoupler.pack(side=LEFT, padx=(0, 5))


def generateCircuitInputControls(circuitFrame):
    controlsInputFrameOuter = Frame(circuitFrame, height=80, width=relativeWidth*width*0.40)  # , background="orange")
    controlsInputFrameOuter.grid(row=4, column=0, columnspan=3)

    controlsInputFrameInner = Frame(controlsInputFrameOuter)
    controlsInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    loadCircuitButton = Button(controlsInputFrameInner, text="Load Circuit", command=loadCircuit, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    loadCircuitButton.grid(row=0, column=0, padx=4, pady=4)

    changeDefaultCircuitButton = Button(controlsInputFrameInner, text="Save as Default", command=changeDefaultCircuit, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    changeDefaultCircuitButton.grid(row=0, column=1, padx=4, pady=4)

    useDefaultCircuitButton = Button(controlsInputFrameInner, text="Use Default", command=useDefaultCircuit, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    useDefaultCircuitButton.grid(row=0, column=2, padx=4, pady=4)

    generateCircuitButton = Button(controlsInputFrameInner, text="Export", command=exportCircuit, padx=3, pady=3, background="#21e4d7", relief=FLAT)
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


def optimizeCircuitParameterFrame(topFrame):
    # , background='blue')#, highlightbackground="black", highlightthickness=1)
    circuitFrame = Frame(topFrame, height=relativeHeight * height, width=relativeWidth*width*0.40)
    circuitFrame.grid(column=0, row=0)
    circuitFrame.grid_propagate(0)

    generateCircuitInputs(circuitFrame)
    generateCircuitInputControls(circuitFrame)


##################################################################################
# Functions that generate the controls for optimizer configuration, boundary search limits and signal shape.


def generateModulationTimeBoundaryInput(inputBoundaryFrame, entryCharacterWidth):
    modulationTimeInputFrameOuter = Frame(inputBoundaryFrame, height=35, width=relativeWidth*width*0.60)  # , background="green")
    modulationTimeInputFrameOuter.grid(row=4, column=0, columnspan=3)

    modulationTimeInputFrameInner = Frame(modulationTimeInputFrameOuter)  # , background="blue")
    modulationTimeInputFrameInner.place(anchor="e", relx=0.82, rely=0.5)

    modulationTimeLabel = Label(modulationTimeInputFrameInner, text="Total modulation time of AC-flux [ns]:")
    modulationTimeLabel.pack(side=LEFT, padx=(0, 5))

    lowerModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Lower limit:")
    lowerModulationTimeLabel.pack(side=LEFT)

    lowerModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeLower)
    lowerModulationTimeEntry.pack(side=LEFT, padx=(0, 5))

    upperModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Upper limit:")
    upperModulationTimeLabel.pack(side=LEFT)

    upperModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeUpper)
    upperModulationTimeEntry.pack(side=LEFT)


def generateOmegaPhiBoundaryInput(inputBoundaryFrame, entryCharacterWidth):
    omegaPhiInputFrameOuter = Frame(inputBoundaryFrame, height=35, width=relativeWidth*width*0.60)  # , background="yellow")
    omegaPhiInputFrameOuter.grid(row=3, column=0, columnspan=3)

    omegaPhiInputFrameInner = Frame(omegaPhiInputFrameOuter)  # , background="blue")
    omegaPhiInputFrameInner.place(anchor="e", relx=0.82, rely=0.5)

    omegaPhiLabel = Label(omegaPhiInputFrameInner, text="Frequency \u03C9 of AC-flux [GHz]:")
    omegaPhiLabel.pack(side=LEFT, padx=(0, 5))

    lowerOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Lower limit:")
    lowerOmegaPhiLabel.pack(side=LEFT)

    lowerOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiLower)
    lowerOmegaPhiEntry.pack(side=LEFT, padx=(0, 5))

    upperOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Upper limit:")
    upperOmegaPhiLabel.pack(side=LEFT)

    upperOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiUpper)
    upperOmegaPhiEntry.pack(side=LEFT)


def generateX1BoundaryInput(inputBoundaryFrame, entryCharacterWidth):
    x1InputFrameOuter = Frame(inputBoundaryFrame, height=35, width=relativeWidth*width*0.60)  # , background="green")
    x1InputFrameOuter.grid(row=2, column=0, columnspan=3)

    x1InputFrameInner = Frame(x1InputFrameOuter)  # , background="blue")
    x1InputFrameInner.place(anchor="e", relx=0.82, rely=0.5)

    x1LabelVar.set("Amplitude of \u03B4(t) [\u03A6"+subscriptZero+"]:")
    x1Label = Label(x1InputFrameInner, textvariable=x1LabelVar)
    x1Label.pack(side=LEFT, padx=(0, 5))

    lowerX1Label = Label(x1InputFrameInner, text="Lower limit:")
    lowerX1Label.pack(side=LEFT)

    lowerX1Entry = Entry(x1InputFrameInner, width=entryCharacterWidth, textvariable=x1Lower)
    lowerX1Entry.pack(side=LEFT, padx=(0, 5))

    upperX1Label = Label(x1InputFrameInner, text="Upper limit:")
    upperX1Label.pack(side=LEFT)

    upperX1Entry = Entry(x1InputFrameInner, width=entryCharacterWidth, textvariable=x1Upper)
    upperX1Entry.pack(side=LEFT)


def generateX0BoundaryInput(inputBoundaryFrame, entryCharacterWidth):
    x0InputFrameOuter = Frame(inputBoundaryFrame, height=35, width=relativeWidth*width*0.60)  # , background="yellow")
    x0InputFrameOuter.grid(row=1, column=0, columnspan=3)

    x0InputFrameInner = Frame(x0InputFrameOuter)  # , background="blue")
    x0InputFrameInner.place(anchor="e", relx=0.82, rely=0.5)

    x0LabelVar.set("Strength of DC-flux \u0398 [\u03A6"+subscriptZero+"]:")
    x0Label = Label(x0InputFrameInner, textvariable=x0LabelVar)
    x0Label.pack(side=LEFT, padx=(0, 5))

    lowerX0Label = Label(x0InputFrameInner, text="Lower limit:")
    lowerX0Label.pack(side=LEFT)

    lowerX0Entry = Entry(x0InputFrameInner, width=entryCharacterWidth, textvariable=x0Lower)
    lowerX0Entry.pack(side=LEFT, padx=(0, 5))

    upperX0Label = Label(x0InputFrameInner, text="Upper limit:")
    upperX0Label.pack(side=LEFT)

    upperX0Entry = Entry(x0InputFrameInner, width=entryCharacterWidth, textvariable=x0Upper)
    upperX0Entry.pack(side=LEFT)


def generateBoundaryInput(inputBoundaryFrame):
    entryCharacterWidth = 6

    riseTimeFrameOuter = Frame(inputBoundaryFrame, height=35, width=relativeWidth*width*0.60)  # , background="red")
    riseTimeFrameOuter.grid(row=0, column=0, columnspan=3)

    riseTimeFrame = Frame(riseTimeFrameOuter)  # , background="blue")
    riseTimeFrame.place(anchor="center", relx=0.5, rely=0.5)

    riseTimeLabel = Label(riseTimeFrame, text="Rise time of \u03B4(t) (from 0% to 100%) [ns]:")
    riseTimeLabel.pack(side=LEFT, padx=(0, 5))

    riseTimeEntry = ttk.Spinbox(riseTimeFrame, from_=15, to=100, textvariable=riseTime, width=4, state="readonly")
    riseTimeEntry.set(25)
    riseTimeEntry.pack(side=LEFT)

    generateX0BoundaryInput(inputBoundaryFrame, entryCharacterWidth)
    generateX1BoundaryInput(inputBoundaryFrame, entryCharacterWidth)
    generateOmegaPhiBoundaryInput(inputBoundaryFrame, entryCharacterWidth)
    generateModulationTimeBoundaryInput(inputBoundaryFrame, entryCharacterWidth)


def generateSelectEnergyLevels(settingsFrameRight):
    selectEnergyLevelsFrameOuter = Frame(settingsFrameRight, height=55, width=relativeWidth*width*0.30)  # , background="blue")
    selectEnergyLevelsFrameOuter.grid(row=1, column=0)

    selectEnergyLevelsFrameInner = Frame(selectEnergyLevelsFrameOuter)  # , background="orange")
    selectEnergyLevelsFrameInner.place(anchor="n", relx=0.5, rely=0.1)

    selectEnergyLevelsTitle = Label(selectEnergyLevelsFrameInner, text="Number of energy levels per qubit:")
    selectEnergyLevelsTitle.pack(side=LEFT, padx=(0, 5))

    selectEnergyLevels = ttk.Spinbox(selectEnergyLevelsFrameInner, from_=2, to=7, textvariable=energyLevels, width=4, state="readonly")
    selectEnergyLevels.set(4)
    selectEnergyLevels.pack(side=LEFT)


def generateSelectGate(settingsFrameRight):
    selectGateFrameOuter = Frame(settingsFrameRight, height=55, width=relativeWidth*width*0.30)  # , background="orange")
    selectGateFrameOuter.grid(row=0, column=0)
    selectGateFrameOuter.pack_propagate(0)

    selectGateFrameInner = LabelFrame(selectGateFrameOuter, text="Gate to optimize for: ")
    selectGateFrameInner.place(anchor='n', relx=0.5, rely=0)

    selectISWAP = Radiobutton(selectGateFrameInner, text="iSWAP", value="iSWAP", variable=selectedGate)
    selectISWAP.select()
    selectISWAP.pack(side=LEFT)

    selectSWAP = Radiobutton(selectGateFrameInner, text="SWAP", value="SWAP", variable=selectedGate)
    selectSWAP.pack(side=LEFT)

    selectCZ = Radiobutton(selectGateFrameInner, text="CZ", value="CZ", variable=selectedGate)
    selectCZ.pack(side=LEFT)


def generateSelectAlgorithm(settingsFrameLeft):
    selectAlgorithmFrameOuter = Frame(settingsFrameLeft, height=110, width=relativeWidth*width*0.30)  # , background="orange")
    selectAlgorithmFrameOuter.grid(row=0, column=0)

    selectAlgorithmFrameInner = LabelFrame(selectAlgorithmFrameOuter, text="Algorithms to use: ")
    selectAlgorithmFrameInner.place(anchor='n', relx=0.5, rely=0)

    selectDifferentialEvolutionBox = Checkbutton(selectAlgorithmFrameInner, text="Differential Evolution (recommended)", variable=runDifferentialEvolution)
    selectDifferentialEvolutionBox.select()
    selectDifferentialEvolutionBox.grid(row=0, column=0, sticky='w')

    selectSimplicalHomologyBox = Checkbutton(selectAlgorithmFrameInner, text="Simplical Homology Global (Disabled)", variable=runSimplicalHomologyGlobal, state=DISABLED)
    selectSimplicalHomologyBox.grid(row=1, column=0, sticky='w')

    selectDualAnnelingBox = Checkbutton(selectAlgorithmFrameInner, text="Dual Anneling (Disabled)", variable=runDualAnneling, state=DISABLED)
    selectDualAnnelingBox.grid(row=2, column=0, sticky='w')


def generateBoundarySettings(settingsBoundaryFrame):
    selectSignalFrameOuter = Frame(settingsBoundaryFrame, height=35, width=relativeWidth*width*0.60)# , background="blue")
    selectSignalFrameOuter.grid(row=0, column=0)

    selectSignalFrameInner = Frame(selectSignalFrameOuter)
    selectSignalFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectSignalTitle = Label(selectSignalFrameInner, text="Shape of magnetic flux signal (\u03A6):")
    selectSignalTitle.pack(side=LEFT, padx=(0, 5))

    selectSignal = ttk.Combobox(selectSignalFrameInner, state="readonly", values=('\u03A6 = \u0398 + \u03B4(t) \u2022 cos(\u03C9t)', '\u03A6 = 1 / \u03C0 \u2022 arccos[ (A + B(t) \u2022 cos(\u03C9t))\u00B2 ]'), textvariable=signalShape, width=38)
    selectSignal.bind("<<ComboboxSelected>>", updateX0X1Labels)
    selectSignal.current(0)
    selectSignal.pack(side=LEFT)

    inputBoundaryFrame = Frame(settingsBoundaryFrame, height=175, width=relativeWidth*width*0.60)
    inputBoundaryFrame.grid(row=1, column=0)
    inputBoundaryFrame.grid_propagate(0)

    generateBoundaryInput(inputBoundaryFrame)


def generateBoundaryControls(settingsBoundaryFrame):
    boundaryControlsFrameOuter = Frame(settingsBoundaryFrame, height=relativeHeight*height-2*210, width=relativeWidth*width*0.60)# , background="green")
    boundaryControlsFrameOuter.grid(row=2, column=0)

    boundaryControlsFrameInner = Frame(boundaryControlsFrameOuter)
    boundaryControlsFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    setDefaultBoundaryConditionButton = Button(boundaryControlsFrameInner, text="Save as Default", command=setBoundaryDefault, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    setDefaultBoundaryConditionButton.grid(row=0, column=0, padx=4, pady=4)

    useDefaultBoundaryConditionButton = Button(boundaryControlsFrameInner, text="Use Default", command=useBoundaryDefault, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    useDefaultBoundaryConditionButton.grid(row=0, column=1, padx=4, pady=4)

    presetBoundaryConditionsButton = Button(boundaryControlsFrameInner, text="Preset Boundary Conditions", command=presetBoundaryConditions, padx=3, pady=3, background="#21e4d7", relief=FLAT)
    presetBoundaryConditionsButton.grid(row=0, column=2, padx=4, pady=4)


def generateRightSettings(settingsFrameRight):
    generateSelectGate(settingsFrameRight)
    generateSelectEnergyLevels(settingsFrameRight)


def generateLeftSettings(settingsFrameLeft):
    generateSelectAlgorithm(settingsFrameLeft)


def optimizerSettingsFrame(topFrame):
    settingsTitleFrameHeight = 50
    optimizeSettingsFrameHeight = 110

    separator = ttk.Separator(topFrame, orient='vertical')
    separator.grid(row=0, column=1, sticky="nesw")

    settingsFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.60)# , background="yellow")
    settingsFrame.grid(row=0, column=2)
    settingsFrame.grid_propagate(0)

    settingsTitleFrame = Frame(settingsFrame, height=settingsTitleFrameHeight, width=relativeWidth*width*0.60)  # , background="orange")
    settingsTitleFrame.grid(row=0, column=0, columnspan=2)

    settingsTitle = Label(settingsTitleFrame, text='Optimizer configuration', font=('Helvetica', 12))
    titleSettingsFont = tkFont.Font(settingsTitle, settingsTitle.cget("font"))
    titleSettingsFont.configure(underline=True)
    settingsTitle.configure(font=titleSettingsFont)
    settingsTitle.place(anchor='center', relx=0.5, rely=0.5)

    settingsFrameLeft = Frame(settingsFrame, height=optimizeSettingsFrameHeight, width=relativeWidth*width*0.30)  # , background="green")
    settingsFrameLeft.grid(column=0, row=1)
    settingsFrameLeft.grid_propagate(0)

    settingsFrameRight = Frame(settingsFrame, height=optimizeSettingsFrameHeight, width=relativeWidth*width*0.30)  # , background="red")
    settingsFrameRight.grid(column=1, row=1)
    settingsFrameRight.grid_propagate(0)

    separator = ttk.Separator(settingsFrame, orient='horizontal')
    separator.grid(row=2, column=0, columnspan=2, sticky="nesw")

    settingsBoundaryTitleFrame = Frame(settingsFrame, height=settingsTitleFrameHeight, width=relativeWidth*width*0.60)  # , background="orange")
    settingsBoundaryTitleFrame.grid(row=3, column=0, columnspan=2)

    settingsBoundaryTitle = Label(settingsBoundaryTitleFrame, text='Magnetic flux signal and search boundaries', font=('Helvetica', 12))
    settingsBoundaryTitle.configure(font=titleSettingsFont)
    settingsBoundaryTitle.place(anchor='center', relx=0.5, rely=0.5)

    settingsBoundaryFrame = Frame(settingsFrame, height=relativeHeight*height-optimizeSettingsFrameHeight - 2*settingsTitleFrameHeight, width=relativeWidth*width*0.60)#, background="orange")
    settingsBoundaryFrame.grid(row=4, column=0, columnspan=2)
    settingsBoundaryFrame.grid_propagate(0)

    generateLeftSettings(settingsFrameLeft)
    generateRightSettings(settingsFrameRight)
    generateBoundarySettings(settingsBoundaryFrame)
    generateBoundaryControls(settingsBoundaryFrame)


##################################################################################

######################################################################################################################################################################
# Functions that generate different windows


def optimizeControlFrame(root, parentWidget, givenHeight, givenWidth):
    initiateGlobalVariables(root, givenHeight)

    global relativeHeight
    global relativeWidth
    global height
    global width

    height = givenHeight
    width = givenWidth

    windowFrame = Frame(parentWidget, height=height, width=width)# , background="green")
    windowFrame.place(anchor="center", relx=0.5, rely=0.5)
    windowFrame.grid_propagate(0)

    topFrame = Frame(windowFrame, height=relativeHeight*height, width=width)# , background="green")
    topFrame.grid(row=0, column=0)
    topFrame.grid_propagate(0)

    separator = ttk.Separator(windowFrame, orient='horizontal')
    separator.grid(row=1, column=0, sticky="nesw")

    bottomFrame = Frame(windowFrame, height=(1-relativeHeight)*height, width=width)
    bottomFrame.grid(row=2, column=0)
    bottomFrame.grid_propagate(0)

    optimizeStatusFrame(bottomFrame)
    optimizeCircuitParameterFrame(topFrame)
    optimizerSettingsFrame(topFrame)

    useDefaultCircuit()
    useDefaultBoundaryValues()

    return windowFrame


######################################################################################################################################################################