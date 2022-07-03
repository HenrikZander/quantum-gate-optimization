######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: boundaryConditionPopUpWindow.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 29 July 2021

# Copyright 2022, Henrik Zander & Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import json
from tkinter import *
import math
import numpy as np
from numpy.core.numeric import cross

import dataManager
import optimize2QubitFrame as gui

######################################################################################################################################################################
# Global variables

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
subscriptZero = "0".translate(SUB)


def initiateGlobalVariables(root):
    global currentPresetIndex
    currentPresetIndex = IntVar(root)

    global presetName
    presetName = StringVar(root)

    global gate
    gate = StringVar(root)

    global x0LowerNew
    global x0UpperNew
    x0LowerNew = DoubleVar(root)
    x0UpperNew = DoubleVar(root)

    global x1LowerNew
    global x1UpperNew
    x1LowerNew = DoubleVar(root)
    x1UpperNew = DoubleVar(root)

    global omegaPhiLowerNew
    global omegaPhiUpperNew
    omegaPhiLowerNew = DoubleVar(root)
    omegaPhiUpperNew = DoubleVar(root)

    global modulationTimeLowerNew
    global modulationTimeUpperNew
    modulationTimeLowerNew = DoubleVar(root)
    modulationTimeUpperNew = DoubleVar(root)


def setVariablesToZero():
    presetName.set("")

    x0LowerNew.set(0)
    x0UpperNew.set(0)

    x1LowerNew.set(0)
    x1UpperNew.set(0)

    omegaPhiLowerNew.set(0)
    omegaPhiUpperNew.set(0)

    modulationTimeLowerNew.set(0)
    modulationTimeUpperNew.set(0)


def setVariablesToPreset(index, copy=False):
    parameterValues = getPreset(index)

    if copy:
        presetName.set(parameterValues[0] + " (copy)")
    else:
        presetName.set(parameterValues[0])

    signalType = gui.getSignalType()
    #print(gui.signalShape.get())
    #print(signalType)
    if signalType == 'arccos':
        if parameterValues[2][0] * parameterValues[2][1] > 0:
            dcAmplitudeLower = getAfromTheta(max(parameterValues[2], key=abs))
            dcAmplitudeUpper = getAfromTheta(min(parameterValues[2], key=abs))
        else:
            dcAmplitudeLower = getAfromTheta(max(parameterValues[2], key=abs))
            dcAmplitudeUpper = 1

        # OBS: Detta är egentligen ass, men det blir typ det rimligaste iom att
        # det inte går att begränsa B dynamiskt beroende på A (vilket egentligen 
        # är vad man vill göra):
        acAmplitudeLower = parameterValues[3][0]
        acAmplitudeUpper = parameterValues[3][1]
        
        x0LowerNew.set(dcAmplitudeLower)
        x0UpperNew.set(dcAmplitudeUpper)

        x1LowerNew.set(acAmplitudeLower)
        x1UpperNew.set(acAmplitudeUpper)
    elif signalType == 'cos':
        x0LowerNew.set(parameterValues[2][0])
        x0UpperNew.set(parameterValues[2][1])

        x1LowerNew.set(parameterValues[3][0])
        x1UpperNew.set(parameterValues[3][1])
    else:
        print("Signaltype not supported! Must be either arccos or cos!")

    omegaPhiLowerNew.set(parameterValues[4][0])
    omegaPhiUpperNew.set(parameterValues[4][1])

    modulationTimeLowerNew.set(parameterValues[5][0])
    modulationTimeUpperNew.set(parameterValues[5][1])


def getPreviewVariableValues():
    return presetName.get(), [[x0LowerNew.get(), x0UpperNew.get()], [x1LowerNew.get(), x1UpperNew.get()], [omegaPhiLowerNew.get(), omegaPhiUpperNew.get()], [modulationTimeLowerNew.get(), modulationTimeUpperNew.get()]]


def getAfromTheta(theta):
    dcAmplitude = np.sqrt(np.abs(np.cos(np.pi * theta)))
    return dcAmplitude


def getBfromThetaDelta(theta,delta):
    '''
    This function assumes that theta ± delta lies within [-1,1] Phi_0, and that delta is positive.
    '''
    if theta + delta > 0.5:
        acAmplitude = -1/2 * np.sqrt(np.abs(np.cos(np.pi * (theta - delta))))
    elif theta - delta < -0.5:
        acAmplitude = 1/2 * np.sqrt(np.abs(np.cos(np.pi * (theta + delta))))
    else:
        if abs(theta) + abs(delta) > 0.5:
            print("Weird input in getBfromThetaDelta.")
        acAmplitude = 1/2 * ( np.sqrt(np.abs(np.cos(np.pi * (theta + delta)))) - np.sqrt(np.abs(np.cos(np.pi * (theta - delta)))) )
    return acAmplitude


def loadSelectionPreview(event):
    selection = event.widget.curselection()
    if selection:
        index = selection[0]
        setVariablesToPreset(index)


######################################################################################################################################################################
# Functions for dynamically updating certain default presets


def getCrossings():
    freq1 = gui.frequencyQ1.get()
    freq2 = gui.frequencyQ2.get()
    freqTB0 = gui.frequencyCoupler.get()
    anh1 = gui.anharmonicityQ1.get()
    anh2 = gui.anharmonicityQ2.get()

    if freq1 > freq2:
        highFrequency = freq1
        highAnh = anh1
        lowFrequency = freq2
        lowAnh = anh2

        upperCZ = 'CZ_20'
        lowerCZ = 'CZ_02'
    else:
        highFrequency = freq2
        highAnh = anh2
        lowFrequency = freq1
        lowAnh = anh1

        upperCZ = 'CZ_02'
        lowerCZ = 'CZ_20'

    nearCrossing_iSWAP = np.arccos((highFrequency/freqTB0)**2) / np.pi
    farCrossing_iSWAP = np.arccos((lowFrequency/freqTB0)**2) / np.pi

    nearCrossing_upperCZ = np.arccos(((highFrequency + highAnh)/freqTB0)**2) / np.pi

    farCrossing_lowerCZ = np.arccos(((lowFrequency + lowAnh)/freqTB0)**2) / np.pi

    crossingsDict = {}
    crossingsDict['iSWAP'] = [nearCrossing_iSWAP, farCrossing_iSWAP]
    crossingsDict[upperCZ] = [nearCrossing_upperCZ, farCrossing_iSWAP]
    crossingsDict[lowerCZ] = [nearCrossing_iSWAP, farCrossing_lowerCZ]
    return crossingsDict


def updateStaticPresets():
    configData = dataManager.getFromjson("config.json")
    presetList = configData["boundaryPresets"]
    crossingsDict = getCrossings()

    maxTheta = 0.48
    minThetaSWAP = 0.3
    tradFreqDeviation = 0.04
    tradFreq_iSWAP = abs(gui.frequencyQ1.get() - gui.frequencyQ2.get())
    tradFreq_CZ_20 = abs(gui.frequencyQ1.get() + gui.anharmonicityQ1.get() - gui.frequencyQ2.get())
    tradFreq_CZ_02 = abs(gui.frequencyQ1.get() - gui.frequencyQ2.get() - gui.anharmonicityQ2.get())

    for p in presetList:
        if p[0] == "Traditional iSWAP":
            p[2] = [-crossingsDict['iSWAP'][0], 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_iSWAP - 2*tradFreqDeviation, tradFreq_iSWAP + tradFreqDeviation]
            p[5] = [50, 100]
        elif p[0] == "Traditional CZ (|11> <-> |20>)":
            p[2] = [-crossingsDict['CZ_20'][0], 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_CZ_20 - 2*tradFreqDeviation, tradFreq_CZ_20 + tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Traditional CZ (|11> <-> |02>)":
            p[2] = [-crossingsDict['CZ_02'][0], 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_CZ_02 - 2*tradFreqDeviation, tradFreq_CZ_02 + tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Intercrossing iSWAP":
            p[2] = [-crossingsDict['iSWAP'][1], -crossingsDict['iSWAP'][0]]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_iSWAP - 0.5*tradFreqDeviation, tradFreq_iSWAP + 2*tradFreqDeviation]
            p[5] = [50, 100]
        elif p[0] == "Intercrossing CZ (|11> <-> |20>)":
            p[2] = [-crossingsDict['CZ_20'][1], -crossingsDict['CZ_20'][0]]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_CZ_20 - 0.5*tradFreqDeviation, tradFreq_CZ_20 + 2*tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Intercrossing CZ (|11> <-> |02>)":
            p[2] = [-crossingsDict['CZ_02'][1], -crossingsDict['CZ_02'][0]]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_CZ_02 - 0.5*tradFreqDeviation, tradFreq_CZ_02 + 2*tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Constrained SWAP":
            p[2] = [-maxTheta, -minThetaSWAP]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_iSWAP - 2*tradFreqDeviation, tradFreq_iSWAP + 2*tradFreqDeviation]
            p[5] = [50, 150]
    
    dataManager.dumpTojson(configData,'config.json')


######################################################################################################################################################################
# The pop up window that lets the user select a preset boundary condition.


def selectPresetWindow(root, selectedGate):
    global suppliedRoot
    global pop
    suppliedRoot = root

    pop = Toplevel(root)
    pop.title("Boundary condition manager")

    height = 500
    width = 935

    pop.geometry(str(width)+"x"+str(height))
    pop.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    pop.iconphoto(False, programIcon)
    
    initiateGlobalVariables(suppliedRoot)
    gate.set(selectedGate)

    generateBoundaryManager(pop, height, width)
    loadAllPresets()

    pop.grab_set()


######################################################################################################################################################################
# The pop up window that is used to preview boundary presets and create new ones.


def interactWithPresetWindow(index, addNew=False):
    global suppliedRoot
    global popPreview
    popPreview = Toplevel(suppliedRoot)
    popPreview.protocol("WM_DELETE_WINDOW", closePreviewWindow)

    height = 300
    width = 575

    currentPresetIndex.set(index)

    popPreview.geometry(str(width)+"x"+str(height))
    popPreview.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    popPreview.iconphoto(False, programIcon)

    if addNew:
        popPreview.title("Add new boundary preset")
        generateBoundaryPresetPreview(popPreview, height, width)
        generatePresetPreviewControls(popPreview, height, width)

        if index is None:
            setVariablesToZero()
        else:
            setVariablesToPreset(index, copy=True)
    else:
        popPreview.title("Edit boundary preset")
        generateBoundaryPresetPreview(popPreview, height, width, edit=True)
        generatePresetPreviewControls(popPreview, height, width, edit=True)

        setVariablesToPreset(index)

    popPreview.grab_set()


######################################################################################################################################################################
# Function that generates the widgets in the windows.


def generateBoundaryManager(pop, height, width):
    leftFrame = Frame(pop, height=height, width=width*0.4)# , background="green")
    leftFrame.pack_propagate(0)
    leftFrame.grid(row=0, column=0)

    rightFrame = Frame(pop, height=height, width=width*0.6)# , background="red")
    rightFrame.pack_propagate(0)
    rightFrame.grid(row=0, column=1)

    windowTitleFrame = Frame(leftFrame, width=width*0.4, height=height*0.05)
    windowTitleFrame.pack()

    windowTitle = Label(windowTitleFrame, text="Available presets:", font=('Helvetica', 10))
    windowTitle.place(anchor="center", relx=0.5, rely=0.5)

    listboxFrame = Frame(leftFrame)
    listboxFrame.pack()

    global presetsBox
    listboxScrollbar = Scrollbar(listboxFrame, orient=VERTICAL)
    presetsBox = Listbox(listboxFrame, width=50, height=25, relief=FLAT, selectmode=SINGLE, highlightthickness=0, bd=0, background="#E0E0E0", yscrollcommand=listboxScrollbar.set)
    listboxScrollbar.config(command=presetsBox.yview)
    listboxScrollbar.pack(side=RIGHT, fill=Y)
    presetsBox.pack()

    presetsBox.bind("<<ListboxSelect>>", loadSelectionPreview)

    buttonFrameOuter = Frame(leftFrame, width=width*0.4, height=height*0.1)
    buttonFrameOuter.pack()

    buttonFrameInner = Frame(buttonFrameOuter)
    buttonFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectPresetButton = Button(buttonFrameInner, text="Select", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=selectPreset)
    selectPresetButton.pack(side=LEFT, padx=4, pady=4)

    addNewPresetButton = Button(buttonFrameInner, text="Add New/Copy Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=addNewPreset)
    addNewPresetButton.pack(side=LEFT, padx=4, pady=4)

    editPresetButton = Button(buttonFrameInner, text="Edit Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=editPreset)
    editPresetButton.pack(side=LEFT, padx=4, pady=4)

    selectionPreview(rightFrame, height, width)


def selectionPreview(parentWidget, height, width):
    entryCharacterWidth = 6
    width = width*0.58

    parameterFrame = LabelFrame(parentWidget, text="Parameters of selected preset: ", width=width, height=height*0.32)
    parameterFrame.grid_propagate(0)
    parameterFrame.place(anchor="nw", relx=0.02, rely=0.05)

    width = width*0.95
    generateX0InputWidgets(parameterFrame, height, width, entryCharacterWidth, readOnlyState=True)
    generateX1InputWidgets(parameterFrame, height, width, entryCharacterWidth, readOnlyState=True)
    generateOmegaInputWidgets(parameterFrame, height, width, entryCharacterWidth, readOnlyState=True)
    generateModulationTimeInputWidgets(parameterFrame, height, width, entryCharacterWidth, readOnlyState=True)


def generateBoundaryPresetPreview(parentWidget, height, width, edit=False):
    entryCharacterWidth = 6
    width = width*0.9

    nameFrame = LabelFrame(parentWidget, text="Preset name: ")
    nameFrame.pack(pady=(8,0))

    nameEntry = Entry(nameFrame, width=55, textvariable=presetName)
    nameEntry.pack(padx=5, pady=5)

    if edit:
        nameEntry.config(state="readonly")

    parameterFrame = LabelFrame(parentWidget, text="Preset parameters: ")
    parameterFrame.pack(pady=(8,0))

    generateX0InputWidgets(parameterFrame, height, width, entryCharacterWidth)
    generateX1InputWidgets(parameterFrame, height, width, entryCharacterWidth)
    generateOmegaInputWidgets(parameterFrame, height, width, entryCharacterWidth)
    generateModulationTimeInputWidgets(parameterFrame, height, width, entryCharacterWidth)


def generatePresetPreviewControls(parentWidget, height, width, edit=False):
    previewControlFrame = Frame(parentWidget, height=height*0.2, width=width)
    previewControlFrame.pack()

    if edit:
        editorControlFrame = Frame(previewControlFrame)
        editorControlFrame.place(anchor='center', relx=0.5, rely=0.5)

        editorSaveButton = Button(editorControlFrame, text="Save and close", command=saveEditedPreset, padx=3, pady=3, background="#21e4d7", relief=FLAT)
        editorSaveButton.grid(row=0, column=0, padx=4, pady=4)

        editorDeleteButton = Button(editorControlFrame, text="Delete", command=deleteEditedPreset, padx=3, pady=3, background="#21e4d7", relief=FLAT)
        editorDeleteButton.grid(row=0, column=1, padx=4, pady=4)
    else:
        creatorControlFrame = Frame(previewControlFrame)
        creatorControlFrame.place(anchor='center', relx=0.5, rely=0.5)

        creatorSaveButton = Button(creatorControlFrame, text="Save and Close", command=saveNewPreset, padx=3, pady=3, background="#21e4d7", relief=FLAT)
        creatorSaveButton.grid(row=0, column=0, padx=4, pady=4)

        creatorCancelButton = Button(creatorControlFrame, text="Cancel", command=cancelNewPreset, padx=3, pady=3, background="#21e4d7", relief=FLAT)
        creatorCancelButton.grid(row=0, column=1, padx=4, pady=4)


######################################################################################################################################################################
# Button callback functions


def closePreviewWindow():
    global popPreview
    popPreview.grab_release()
    popPreview.destroy()

    global pop
    pop.grab_set()


def deleteEditedPreset():
    deletePreset(currentPresetIndex.get())
    closePreviewWindow()


def saveEditedPreset():
    _, presetParameterValues = getPreviewVariableValues()
    changePreset(currentPresetIndex.get(), presetParameterValues)

    closePreviewWindow()


def saveNewPreset():
    name, presetParameterValues = getPreviewVariableValues()
    addPreset(name, presetParameterValues)

    closePreviewWindow()


def cancelNewPreset():
    closePreviewWindow()


def selectPreset():
    global presetsBox
    selectedPresetIndex = presetsBox.curselection()
    if len(selectedPresetIndex) > 0:
        selectedPresetIndex = selectedPresetIndex[0]

        boundaryValues = getPreset(selectedPresetIndex)

        if gui.getSignalType() == 'arccos':
            for i in range(2):
                boundaryValues[2][i] = getAfromTheta(boundaryValues[2][i])
            # OBS: Jag laddar här in deltaBounds som om de vore BBounds.
        
        gui.setBoundaryValues(boundaryValues)
        # print(boundaryValues)

        global pop
        pop.grab_release()
        pop.destroy()
    else:
        print("No preset selected!")


def addNewPreset():
    global presetsBox
    selectedPresetIndexTuple = presetsBox.curselection()

    if len(selectedPresetIndexTuple) > 0:
        selectedPresetIndex = selectedPresetIndexTuple[0]
        interactWithPresetWindow(selectedPresetIndex, addNew=True)
    else:
        interactWithPresetWindow(None, addNew=True)


def editPreset():
    global presetsBox
    selectedPresetIndexTuple = presetsBox.curselection()

    if len(selectedPresetIndexTuple) > 0:
        selectedPresetIndex = selectedPresetIndexTuple[0]
        interactWithPresetWindow(selectedPresetIndex)
    else:
        print("No preset selected!")


######################################################################################################################################################################
# Functions for handling all the preset items


def loadAllPresets():
    global presetsBox
    presetsBox.delete(0,END)
    
    selectedGate = gate.get()
    configData = dataManager.getFromjson("config.json")
    presetList = configData["boundaryPresets"]

    for item in presetList:
        if item[6] in [selectedGate, "common"]:
            presetsBox.insert(END, item[0])


def addPreset(name, boundaryValues):
    if gui.getSignalType() == 'arccos':
        for i in range(2):
            boundaryValues[0][i] = - 1/np.pi * np.arccos(boundaryValues[0][i] ** 2)
        # OBS: Jag sparar här BBounds som om de vore deltaBounds, eftersom jag ovan valde att använda ett preset:s deltaBounds som BBounds rakt av.
    configData = dataManager.getFromjson("config.json")
    configData["boundaryPresets"].append([name, True, *boundaryValues, gate.get()])
    dataManager.dumpTojson(configData, "config.json")

    loadAllPresets()


def getPreset(index):
    # Get list of presets from config.json
    configData = dataManager.getFromjson("config.json")
    presetList = configData["boundaryPresets"]
    # Only care about common presets and presets specific to the currently selected gate
    try:
        selectedGate = gate.get()
    except: # One could perhaps argue for choosing "common" instead of "iSWAP" as the default here, but I don't think it matters that much /Emil
        selectedGate = "iSWAP"
    presetList[:] = [preset for preset in presetList if preset[6] in [selectedGate, "common"]]
    #print(presetList)
    return presetList[index]


def changePreset(index, newBoundaryValues):
    if gui.getSignalType() == 'arccos':
        for i in range(2):
            newBoundaryValues[0][i] = - 1/np.pi * np.arccos(newBoundaryValues[0][i] ** 2)
        # OBS: Jag sparar här BBounds som om de vore deltaBounds, eftersom jag ovan valde att använda ett preset:s deltaBounds som BBounds rakt av.
    configData = dataManager.getFromjson("config.json")
    oldBoundaryValues = configData["boundaryPresets"][index]
    configData["boundaryPresets"][index] = [oldBoundaryValues[0], oldBoundaryValues[1], *newBoundaryValues]
    dataManager.dumpTojson(configData, "config.json")

    loadAllPresets()


def deletePreset(index):
    configData = dataManager.getFromjson("config.json")
    allPreset = configData["boundaryPresets"]
    del allPreset[index]
    configData["boundaryPresets"] = allPreset
    dataManager.dumpTojson(configData, "config.json")

    loadAllPresets()


######################################################################################################################################################################
# Functions that help generate the widgets in the creator/editor window.


def generateX0InputWidgets(parentWidget, height, width, entryCharacterWidth, readOnlyState=False):
    if readOnlyState:
        state = "readonly"
    else:
        state = NORMAL

    x0InputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="yellow")
    x0InputFrameOuter.grid(row=0, column=0)

    x0InputFrameInner = Frame(x0InputFrameOuter)# , background="blue")
    x0InputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    x0Label = Label(x0InputFrameInner, textvariable=gui.x0LabelVar)
    x0Label.pack(side=LEFT, padx=(0, 5))

    lowerX0Label = Label(x0InputFrameInner, text="Lower limit:")
    lowerX0Label.pack(side=LEFT)

    lowerX0Entry = Entry(x0InputFrameInner, width=entryCharacterWidth, textvariable=x0LowerNew, state=state)
    lowerX0Entry.pack(side=LEFT, padx=(0, 5))

    upperX0Label = Label(x0InputFrameInner, text="Upper limit:")
    upperX0Label.pack(side=LEFT)

    upperX0Entry = Entry(x0InputFrameInner, width=entryCharacterWidth, textvariable=x0UpperNew, state=state)
    upperX0Entry.pack(side=LEFT)


def generateX1InputWidgets(parentWidget, height, width, entryCharacterWidth, readOnlyState=False):
    if readOnlyState:
        state = "readonly"
    else:
        state = NORMAL
    
    x1InputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="green")
    x1InputFrameOuter.grid(row=1, column=0)

    x1InputFrameInner = Frame(x1InputFrameOuter)  # , background="blue")
    x1InputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    x1Label = Label(x1InputFrameInner, textvariable=gui.x1LabelVar)
    x1Label.pack(side=LEFT, padx=(0, 5))

    lowerX1Label = Label(x1InputFrameInner, text="Lower limit:")
    lowerX1Label.pack(side=LEFT)

    lowerX1Entry = Entry(x1InputFrameInner, width=entryCharacterWidth, textvariable=x1LowerNew, state=state)
    lowerX1Entry.pack(side=LEFT, padx=(0, 5))

    upperX1Label = Label(x1InputFrameInner, text="Upper limit:")
    upperX1Label.pack(side=LEFT)

    upperX1Entry = Entry(x1InputFrameInner, width=entryCharacterWidth, textvariable=x1UpperNew, state=state)
    upperX1Entry.pack(side=LEFT)


def generateOmegaInputWidgets(parentWidget, height, width, entryCharacterWidth, readOnlyState=False):
    if readOnlyState:
        state = "readonly"
    else:
        state = NORMAL
    
    omegaPhiInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="yellow")
    omegaPhiInputFrameOuter.grid(row=2, column=0)

    omegaPhiInputFrameInner = Frame(omegaPhiInputFrameOuter)  # , background="blue")
    omegaPhiInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    omegaPhiLabel = Label(omegaPhiInputFrameInner, text="Frequency \u03C9 of AC-flux [GHz]:")
    omegaPhiLabel.pack(side=LEFT, padx=(0, 5))

    lowerOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Lower limit:")
    lowerOmegaPhiLabel.pack(side=LEFT)

    lowerOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiLowerNew, state=state)
    lowerOmegaPhiEntry.pack(side=LEFT, padx=(0, 5))

    upperOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Upper limit:")
    upperOmegaPhiLabel.pack(side=LEFT)

    upperOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiUpperNew, state=state)
    upperOmegaPhiEntry.pack(side=LEFT)


def generateModulationTimeInputWidgets(parentWidget, height, width, entryCharacterWidth, readOnlyState=False):
    if readOnlyState:
        state = "readonly"
    else:
        state = NORMAL
    
    modulationTimeInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="green")
    modulationTimeInputFrameOuter.grid(row=4, column=0, columnspan=3)

    modulationTimeInputFrameInner = Frame(modulationTimeInputFrameOuter)  # , background="blue")
    modulationTimeInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    modulationTimeLabel = Label(modulationTimeInputFrameInner, text="Total modulation time of AC-flux [ns]:")
    modulationTimeLabel.pack(side=LEFT, padx=(0, 5))

    lowerModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Lower limit:")
    lowerModulationTimeLabel.pack(side=LEFT)

    lowerModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeLowerNew, state=state)
    lowerModulationTimeEntry.pack(side=LEFT, padx=(0, 5))

    upperModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Upper limit:")
    upperModulationTimeLabel.pack(side=LEFT)

    upperModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeUpperNew, state=state)
    upperModulationTimeEntry.pack(side=LEFT)


######################################################################################################################################################################