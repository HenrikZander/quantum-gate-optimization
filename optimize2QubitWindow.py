from tkinter import *
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time

from main import *

######################################################################################################################################################################
# Global variables and function to generate global variables

height = 600
width = int(1.62*height)
cancelOptimization = False
relativeHeight = 0.85
relativeWidth = 1

def generateGlobalVariables(root):
    global runDifferentialEvolution
    global runSimplicalHomologyGlobal
    global runDualAnneling
    global selectedGate
    global energyLevels

    runDifferentialEvolution = BooleanVar(root)
    runSimplicalHomologyGlobal = BooleanVar(root)
    runDualAnneling = BooleanVar(root)
    selectedGate = StringVar(root)
    energyLevels = IntVar(root)

######################################################################################################################################################################
# Button callbacks

def startOptimizing():
    progress.start(5)
    # print("Start")
    startOptimizeButton.config(background="grey", command=NONE)
    stopOptimizeButton.config(background="red", command=stopOptimizing)
    progressLabel.config(text="Status: Optimization running, please wait.")


def stopOptimizing():
    # print("Stop")
    startOptimizeButton.config(background="green", command=startOptimizing)
    stopOptimizeButton.config(background="grey", command=NONE)
    progressLabel.config(text="Status: Optimization stopped. Ready to optimize.")
    progress.stop()
    progress["value"] = 0


def loadCircuit():
    print(runDifferentialEvolution.get())
    print("Circuit loaded.")


def setDefaultCircuit():
    print("Default Circuit set.")


def useDefaultCircuit():
    print("Reset to default circuit.")


def exportCircuit():
    print("Circuit exported to JSON.")

######################################################################################################################################################################
# Functions that generate the window for optimization control of two qubit gates.

def optimizeStatusFrame(bottomFrame):
    global startOptimizeButton 
    startOptimizeButton = Button(bottomFrame, text="Start Optimizing", command=startOptimizing, padx=3, pady=3, background="green")
    startOptimizeButton.grid(column=2, row=0)

    global stopOptimizeButton 
    stopOptimizeButton = Button(bottomFrame, text="Stop Optimizing", padx=3, pady=3, background="grey")
    stopOptimizeButton.grid(column=2, row=1)

    progressFrame = Frame(bottomFrame, height=(1-relativeHeight)*height, width=relativeWidth*width*0.80)
    progressFrame.grid(column=0, row=0, rowspan=2, columnspan=2)
    progressFrame.grid_propagate(0)

    global progressLabel
    progressLabel = Label(bottomFrame, text="Status: Ready to optimize.", font=("Helvetica", 11))
    progressLabel.place(relx=0.05, rely=0.5, anchor='sw')

    global progress
    progress = ttk.Progressbar(bottomFrame, length=width*0.67, orient=HORIZONTAL, mode='determinate', value=0)
    progress.place(relx=0.05, rely=0.5, anchor='nw')


def generateCircuitInputControls(circuitFrame):
    controlsInputFrameOuter = Frame(circuitFrame, height=120, width=relativeWidth*width*0.40)#, background="orange")
    controlsInputFrameOuter.grid(row=7, column=0, columnspan=3)

    controlsInputFrameInner = Frame(controlsInputFrameOuter)
    controlsInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    global loadCircuitButton
    loadCircuitButton = Button(controlsInputFrameInner, text="Load Circuit", command=loadCircuit, padx=3, pady=3, background="#21e4d7")
    loadCircuitButton.grid(row=0, column=0, padx=4, pady=4)

    global setDefaultCircuitButton
    setDefaultCircuitButton = Button(controlsInputFrameInner, text="Set to Default", command=setDefaultCircuit, padx=3, pady=3, background="#21e4d7")
    setDefaultCircuitButton.grid(row=0, column=1, padx=4, pady=4)

    global useDefaultCircuitButton
    useDefaultCircuitButton = Button(controlsInputFrameInner, text="Use Default", command=useDefaultCircuit, padx=3, pady=3, background="#21e4d7")
    useDefaultCircuitButton.grid(row=0, column=2, padx=4, pady=4)

    global generateCircuitButton
    generateCircuitButton = Button(controlsInputFrameInner, text="Export", command=exportCircuit, padx=3, pady=3, background="#21e4d7")
    generateCircuitButton.grid(row=0, column=3, padx=4, pady=4)


def generateCouplingInputs(circuitFrame, entryCharacterWidth):
    couplingTitleFrame = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="blue")
    couplingTitleFrame.grid(row=5, column=0, columnspan=3)
    
    couplingTitle = Label(couplingTitleFrame, text="Coupling strength to coupler [GHz]:")
    couplingTitle.place(anchor='center', relx=0.5, rely=0.5)

    couplingInputFrameOuter = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="orange")
    couplingInputFrameOuter.grid(row=6, column=0, columnspan=3)

    couplingInputFrameInner = Frame(couplingInputFrameOuter)
    couplingInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    couplingLabelQ1 = Label(couplingInputFrameInner, text="Q1: ")
    couplingLabelQ1.pack(side=LEFT, padx=(0,0))

    global couplingEntryQ1
    couplingEntryQ1 = Entry(couplingInputFrameInner, width=entryCharacterWidth)
    couplingEntryQ1.pack(side=LEFT)

    couplingLabelQ2 = Label(couplingInputFrameInner, text="Q2: ")
    couplingLabelQ2.pack(side=LEFT, padx=(10,0))

    global couplingEntryQ2
    couplingEntryQ2 = Entry(couplingInputFrameInner, width=entryCharacterWidth)
    couplingEntryQ2.pack(side=LEFT)

    couplingLabelCoupler = Label(couplingInputFrameInner, text="Coupler: ")
    couplingLabelCoupler.pack(side=LEFT, padx=(10,0))

    global couplingEntryCoupler
    couplingEntryCoupler = Entry(couplingInputFrameInner, width=entryCharacterWidth)
    couplingEntryCoupler.pack(side=LEFT)


def generateAnharmonicityInputs(circuitFrame, entryCharacterWidth):
    anharmonicityTitleFrame = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="blue")
    anharmonicityTitleFrame.grid(row=3, column=0, columnspan=3)
    
    anharmonicityTitle = Label(anharmonicityTitleFrame, text="Anharmonicites [GHz]:")
    anharmonicityTitle.place(anchor='center', relx=0.5, rely=0.5)

    anharmonicityInputFrameOuter = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="orange")
    anharmonicityInputFrameOuter.grid(row=4, column=0, columnspan=3)

    anharmonicityInputFrameInner = Frame(anharmonicityInputFrameOuter)
    anharmonicityInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    anharmonicityLabelQ1 = Label(anharmonicityInputFrameInner, text="Q1: ")
    anharmonicityLabelQ1.pack(side=LEFT, padx=(0,0))

    global anharmonicityEntryQ1
    anharmonicityEntryQ1 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth)
    anharmonicityEntryQ1.pack(side=LEFT)

    anharmonicityLabelQ2 = Label(anharmonicityInputFrameInner, text="Q2: ")
    anharmonicityLabelQ2.pack(side=LEFT, padx=(10,0))

    global anharmonicityEntryQ2
    anharmonicityEntryQ2 = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth)
    anharmonicityEntryQ2.pack(side=LEFT)

    anharmonicityLabelCoupler = Label(anharmonicityInputFrameInner, text="Coupler: ")
    anharmonicityLabelCoupler.pack(side=LEFT, padx=(10,0))

    global anharmonicityEntryCoupler
    anharmonicityEntryCoupler = Entry(anharmonicityInputFrameInner, width=entryCharacterWidth)
    anharmonicityEntryCoupler.pack(side=LEFT)


def generateFrequencyInputs(circuitFrame, entryCharacterWidth):
    frequencyTitleFrame = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="blue")
    frequencyTitleFrame.grid(row=1, column=0, columnspan=3)
    
    frequencyTitle = Label(frequencyTitleFrame, text="Frequencies [GHz]:")
    frequencyTitle.place(anchor='center', relx=0.5, rely=0.5)

    frequencyInputFrameOuter = Frame(circuitFrame, height=30, width=relativeWidth*width*0.40)#, background="orange")
    frequencyInputFrameOuter.grid(row=2, column=0, columnspan=3)

    frequencyInputFrameInner = Frame(frequencyInputFrameOuter)
    frequencyInputFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    frequencyLabelQ1 = Label(frequencyInputFrameInner, text="Q1: ")
    frequencyLabelQ1.pack(side=LEFT, padx=(0,0))

    global frequencyEntryQ1
    frequencyEntryQ1 = Entry(frequencyInputFrameInner, width=entryCharacterWidth)
    frequencyEntryQ1.pack(side=LEFT)

    frequencyLabelQ2 = Label(frequencyInputFrameInner, text="Q2: ")
    frequencyLabelQ2.pack(side=LEFT, padx=(10,0))

    global frequencyEntryQ2
    frequencyEntryQ2 = Entry(frequencyInputFrameInner, width=entryCharacterWidth)
    frequencyEntryQ2.pack(side=LEFT)

    frequencyLabelCoupler = Label(frequencyInputFrameInner, text="Coupler: ")
    frequencyLabelCoupler.pack(side=LEFT, padx=(10,0))

    global frequencyEntryCoupler
    frequencyEntryCoupler = Entry(frequencyInputFrameInner, width=entryCharacterWidth)
    frequencyEntryCoupler.pack(side=LEFT)


def generateCircuitInputs(circuitFrame):
    entryCharacterWidth = 10

    circuitTitleFrame = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)#, background="red")
    circuitTitleFrame.grid(row=0, column=0, columnspan=3)

    circuitTitle = Label(circuitTitleFrame, text='Circuit parameters', font=('Helvetica',12))
    titleFont = tkFont.Font(circuitTitle, circuitTitle.cget("font"))
    titleFont.configure(underline = True)
    circuitTitle.configure(font=titleFont)
    circuitTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateFrequencyInputs(circuitFrame, entryCharacterWidth)
    generateAnharmonicityInputs(circuitFrame, entryCharacterWidth)
    generateCouplingInputs(circuitFrame, entryCharacterWidth)


def optimizeCircuitParameterFrame(topFrame):
    circuitFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.40)#, highlightbackground="black", highlightthickness=1, background='blue')
    circuitFrame.grid(column=0, row=0)
    circuitFrame.grid_propagate(0)

    generateCircuitInputs(circuitFrame)
    generateCircuitInputControls(circuitFrame)


def generateSelectAlgorithm(settingsFrameLeft):
    selectAlgorithmTitleFrame = Frame(settingsFrameLeft, height=30, width=relativeWidth*width*0.30, background="blue")
    selectAlgorithmTitleFrame.grid(row=0, column=0)
    
    selectAlgorithmTitle = Label(selectAlgorithmTitleFrame, text="Algorithms to use:")
    selectAlgorithmTitle.place(anchor='center', relx=0.5, rely=0.5)

    selectAlgorithmFrameOuter = Frame(settingsFrameLeft, height=80, width=relativeWidth*width*0.30, background="orange")
    selectAlgorithmFrameOuter.grid(row=1, column=0)

    selectAlgorithmFrameInner = Frame(selectAlgorithmFrameOuter)
    selectAlgorithmFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    selectDifferentialEvolutionBox = Checkbutton(selectAlgorithmFrameInner, text="Differential Evolution (recommended)", variable=runDifferentialEvolution)
    selectDifferentialEvolutionBox.select()
    selectDifferentialEvolutionBox.grid(row=0, column=0, sticky='w')

    selectSimplicalHomologyBox = Checkbutton(selectAlgorithmFrameInner, text="Simplical Homology Global", variable=runSimplicalHomologyGlobal)
    selectSimplicalHomologyBox.grid(row=1, column=0, sticky='w')

    selectDualAnnelingBox = Checkbutton(selectAlgorithmFrameInner, text="Dual Anneling", variable=runDualAnneling)
    selectDualAnnelingBox.grid(row=2, column=0, sticky='w')


def generateSelectGate(settingsFrameRight):
    selectGateTitleFrame = Frame(settingsFrameRight, height=30, width=relativeWidth*width*0.30, background="blue")
    selectGateTitleFrame.grid(row=0, column=0)
    
    selectGateTitle = Label(selectGateTitleFrame, text="Gate to optimize for:")
    selectGateTitle.place(anchor='center', relx=0.5, rely=0.5)

    selectGateFrameOuter = Frame(settingsFrameRight, height=45, width=relativeWidth*width*0.30, background="orange")
    selectGateFrameOuter.grid(row=1, column=0)

    selectGateFrameInner = Frame(selectGateFrameOuter)
    selectGateFrameInner.place(anchor='center', relx=0.5, rely=0.5)

    selectISWAP = Radiobutton(selectGateFrameInner, text="iSWAP", value="iSWAP", variable=selectedGate)
    selectISWAP.select()
    selectISWAP.pack(side=LEFT)

    selectSWAP = Radiobutton(selectGateFrameInner, text="SWAP", value="SWAP", variable=selectedGate)
    selectSWAP.pack(side=LEFT)

    selectCZ = Radiobutton(selectGateFrameInner, text="CZ", value="CZ", variable=selectedGate)
    selectCZ.pack(side=LEFT)


def generateSelectEnergyLevels(settingsFrameRight):
    selectEnergyLevelsFrameOuter = Frame(settingsFrameRight, height=35, width=relativeWidth*width*0.30, background="blue")
    selectEnergyLevelsFrameOuter.grid(row=2, column=0)

    selectEnergyLevelsFrameInner = Frame(selectEnergyLevelsFrameOuter, background="orange")
    selectEnergyLevelsFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectEnergyLevelsTitle = Label(selectEnergyLevelsFrameInner, text="Number of energy levels per qubit:")
    selectEnergyLevelsTitle.pack(side=LEFT, padx=(0,5))

    selectEnergyLevels = ttk.Spinbox(selectEnergyLevelsFrameInner, from_=2, to=7, textvariable=energyLevels, width=4, state="readonly")
    selectEnergyLevels.set(3)
    selectEnergyLevels.pack(side=LEFT)


def generateLeftSettings(settingsFrameLeft):
    generateSelectAlgorithm(settingsFrameLeft)


def generateRightSettings(settingsFrameRight):
    generateSelectGate(settingsFrameRight)
    generateSelectEnergyLevels(settingsFrameRight)


def generateBoundarySettings(settingsBoundaryFrame):
    selectSignalFrameOuter = Frame(settingsBoundaryFrame, height=35, width=relativeWidth*width*0.60, background="blue")
    selectSignalFrameOuter.grid(row=0, column=0)

    selectSignalFrameInner = Frame(selectSignalFrameOuter, background="orange")
    selectSignalFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectSignalTitle = Label(selectSignalFrameInner, text="Shape of magnetic flux signal (\u03A6):")
    selectSignalTitle.pack(side=LEFT, padx=(0,5))

    selectSignal = ttk.Combobox(selectSignalFrameInner, state="readonly", values=('\u03A6 = \u0398 + \u03B4 \u2022 cos(\u03C9t)','Another signal'))
    selectSignal.current(0)
    selectSignal.pack(side=LEFT)


def optimizerSettingsFrame(topFrame):
    settingsTitleFrameHeight = 50
    settingsFrameHeight = 110

    settingsFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.60, background="yellow")
    settingsFrame.grid(column=1, row=0)
    settingsFrame.grid_propagate(0)

    settingsTitleFrame = Frame(settingsFrame, height=settingsTitleFrameHeight, width=relativeWidth*width*0.60, background="orange")
    settingsTitleFrame.grid(row=0, column=0, columnspan=2)

    settingsTitle = Label(settingsTitleFrame, text='Optimizer configuration', font=('Helvetica',12))
    titleSettingsFont = tkFont.Font(settingsTitle, settingsTitle.cget("font"))
    titleSettingsFont.configure(underline = True)
    settingsTitle.configure(font=titleSettingsFont)
    settingsTitle.place(anchor='center', relx=0.5, rely=0.5)

    settingsFrameLeft = Frame(settingsFrame, height=settingsFrameHeight, width=relativeWidth*width*0.30, background="green")
    settingsFrameLeft.grid(column=0, row=1)
    settingsFrameLeft.grid_propagate(0)

    settingsFrameRight = Frame(settingsFrame, height=settingsFrameHeight, width=relativeWidth*width*0.30, background="red")
    settingsFrameRight.grid(column=1, row=1)
    settingsFrameRight.grid_propagate(0)

    settingsBoundaryTitleFrame = Frame(settingsFrame, height=settingsTitleFrameHeight, width=relativeWidth*width*0.60, background="orange")
    settingsBoundaryTitleFrame.grid(row=2, column=0, columnspan=2)

    settingsBoundaryTitle = Label(settingsBoundaryTitleFrame, text='Boundary search limits', font=('Helvetica',12))
    settingsBoundaryTitle.configure(font=titleSettingsFont)
    settingsBoundaryTitle.place(anchor='center', relx=0.5, rely=0.5)

    settingsBoundaryFrame = Frame(settingsFrame, height=relativeHeight*height-settingsFrameHeight, width=relativeWidth*width*0.60, background="green")
    settingsBoundaryFrame.grid(row=3, column=0, columnspan=2)
    settingsBoundaryFrame.grid_propagate(0)

    generateLeftSettings(settingsFrameLeft)
    generateRightSettings(settingsFrameRight)
    generateBoundarySettings(settingsBoundaryFrame)

######################################################################################################################################################################