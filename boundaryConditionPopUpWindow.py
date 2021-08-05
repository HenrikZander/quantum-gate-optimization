######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: boundaryConditionPopUpWindow.py

# Author(s): Henrik Zander

# Date created: 29 July 2021

# Copyright 2021, Henrik Zander, All rights reserved.

######################################################################################################################################################################

import json
from tkinter import *
import math
import numpy as np

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

    global thetaLowerNew
    global thetaUpperNew
    thetaLowerNew = DoubleVar(root)
    thetaUpperNew = DoubleVar(root)

    global deltaLowerNew
    global deltaUpperNew
    deltaLowerNew = DoubleVar(root)
    deltaUpperNew = DoubleVar(root)

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

    thetaLowerNew.set(0)
    thetaUpperNew.set(0)

    deltaLowerNew.set(0)
    deltaUpperNew.set(0)

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

    thetaLowerNew.set(parameterValues[2][0])
    thetaUpperNew.set(parameterValues[2][1])

    deltaLowerNew.set(parameterValues[3][0])
    deltaUpperNew.set(parameterValues[3][1])

    omegaPhiLowerNew.set(parameterValues[4][0])
    omegaPhiUpperNew.set(parameterValues[4][1])

    modulationTimeLowerNew.set(parameterValues[5][0])
    modulationTimeUpperNew.set(parameterValues[5][1])


def getPreviewVariableValues():
    return presetName.get(), [[thetaLowerNew.get(), thetaUpperNew.get()], [deltaLowerNew.get(), deltaUpperNew.get()], [omegaPhiLowerNew.get(), omegaPhiUpperNew.get()], [modulationTimeLowerNew.get(), modulationTimeUpperNew.get()]]


######################################################################################################################################################################
# Functions for dynamically updating certain default presets


def getCrossings():
    lowFrequency = np.minimum(gui.frequencyQ1.get(), gui.frequencyQ2.get())
    highFrequency = np.maximum(gui.frequencyQ1.get(), gui.frequencyQ2.get())

    nearCrossing = np.arccos((highFrequency/gui.frequencyCoupler.get())**2) / np.pi
    farCrossing = np.arccos((lowFrequency/gui.frequencyCoupler.get())**2) / np.pi
    return [nearCrossing, farCrossing]


def updateStaticPresets():
    configData = dataManager.getFromjson("config.json")
    presetList = configData["boundaryPresets"]
    crossings = getCrossings()

    maxTheta = 0.48
    minThetaSWAP = 0.3
    tradFreqDeviation = 0.02
    tradFreq_iSWAP = abs(gui.frequencyQ1.get() - gui.frequencyQ2.get())
    tradFreq_CZ_20 = abs(gui.frequencyQ1.get() + gui.anharmonicityQ1.get() - gui.frequencyQ2.get())
    tradFreq_CZ_02 = abs(gui.frequencyQ1.get() - gui.frequencyQ2.get() - gui.anharmonicityQ2.get())

    for p in presetList:
        if p[0] == "Traditional iSWAP":
            p[2] = [-crossings[0], 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_iSWAP - tradFreqDeviation, tradFreq_iSWAP + tradFreqDeviation]
            p[5] = [50, 100]
        elif p[0] == "Traditional CZ (|11> <-> |20>)":
            p[2] = [-maxTheta, 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_CZ_20 - tradFreqDeviation, tradFreq_CZ_20 + tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Traditional CZ (|11> <-> |02>)":
            p[2] = [-maxTheta, 0]
            p[3] = [0, 0.25]
            p[4] = [tradFreq_CZ_02 - tradFreqDeviation, tradFreq_CZ_02 + tradFreqDeviation]
            p[5] = [50, 120]
        elif p[0] == "Intercrossing iSWAP":
            p[2] = [-crossings[1], -crossings[0]]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_iSWAP - tradFreqDeviation, tradFreq_iSWAP + tradFreqDeviation]
            p[5] = [50, 100]
        elif p[0] == "Constrained SWAP":
            p[2] = [-maxTheta, -minThetaSWAP]
            p[3] = [0, 0.15]
            p[4] = [tradFreq_iSWAP - 2*tradFreqDeviation, tradFreq_iSWAP + 2*tradFreqDeviation]
            p[5] = [50, 150]
    
    dataManager.dumpTojson(configData,'config.json')


######################################################################################################################################################################
# The pop up window that lets the user select a preset boundary condition.


def selectPresetWindow(root):
    global suppliedRoot
    global pop
    suppliedRoot = root

    pop = Toplevel(root)
    pop.title("Boundary condition manager")

    height = 500
    width = 450

    pop.geometry(str(width)+"x"+str(height))
    pop.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    pop.iconphoto(False, programIcon)
    
    initiateGlobalVariables(suppliedRoot)
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
    windowTitleFrame = Frame(pop, width=width*0.85, height=height*0.05)
    windowTitleFrame.pack()

    windowTitle = Label(windowTitleFrame, text="Available presets:", font=('Helvetica', 10))
    windowTitle.place(anchor="w", relx=0, rely=0.5)

    listboxFrame = Frame(pop)
    listboxFrame.pack()

    global presetsBox
    listboxScrollbar = Scrollbar(listboxFrame, orient=VERTICAL)
    presetsBox = Listbox(listboxFrame, width=50, height=25, relief=FLAT, selectmode=SINGLE, highlightthickness=0, bd=0, background="#E0E0E0", yscrollcommand=listboxScrollbar.set)
    listboxScrollbar.config(command=presetsBox.yview)
    listboxScrollbar.pack(side=RIGHT, fill=Y)
    presetsBox.pack()

    buttonFrameOuter = Frame(pop, width=width*0.85, height=height*0.1)
    buttonFrameOuter.pack()

    buttonFrameInner = Frame(buttonFrameOuter)
    buttonFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectPresetButton = Button(buttonFrameInner, text="Select", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=selectPreset)
    selectPresetButton.pack(side=LEFT, padx=4, pady=4)

    addNewPresetButton = Button(buttonFrameInner, text="Add New/Copy Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=addNewPreset)
    addNewPresetButton.pack(side=LEFT, padx=4, pady=4)

    editPresetButton = Button(buttonFrameInner, text="Edit Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=editPreset)
    editPresetButton.pack(side=LEFT, padx=4, pady=4)


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

    generateThetaInputWidgets(parameterFrame, height, width, entryCharacterWidth)
    generateDeltaInputWidgets(parameterFrame, height, width, entryCharacterWidth)
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
    
    configData = dataManager.getFromjson("config.json")
    presetList = configData["boundaryPresets"]

    for item in presetList:
        presetsBox.insert(END, item[0])


def addPreset(name, boundaryValues):
    configData = dataManager.getFromjson("config.json")
    configData["boundaryPresets"].append([name, True, *boundaryValues])
    dataManager.dumpTojson(configData, "config.json")

    loadAllPresets()


def getPreset(index):
    configData = dataManager.getFromjson("config.json")
    return configData["boundaryPresets"][index]


def changePreset(index, newBoundaryValues):
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


def generateThetaInputWidgets(parentWidget, height, width, entryCharacterWidth):
    thetaInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="yellow")
    thetaInputFrameOuter.grid(row=0, column=0)

    thetaInputFrameInner = Frame(thetaInputFrameOuter)  # , background="blue")
    thetaInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    thetaLabel = Label(thetaInputFrameInner, text="Strength of DC-flux \u0398 [\u03A6"+subscriptZero+"]:")
    thetaLabel.pack(side=LEFT, padx=(0, 5))

    lowerThetaLabel = Label(thetaInputFrameInner, text="Lower limit:")
    lowerThetaLabel.pack(side=LEFT)

    lowerThetaEntry = Entry(thetaInputFrameInner, width=entryCharacterWidth, textvariable=thetaLowerNew)
    lowerThetaEntry.pack(side=LEFT, padx=(0, 5))

    upperThetaLabel = Label(thetaInputFrameInner, text="Upper limit:")
    upperThetaLabel.pack(side=LEFT)

    upperThetaEntry = Entry(thetaInputFrameInner, width=entryCharacterWidth, textvariable=thetaUpperNew)
    upperThetaEntry.pack(side=LEFT)


def generateDeltaInputWidgets(parentWidget, height, width, entryCharacterWidth):
    deltaInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="green")
    deltaInputFrameOuter.grid(row=1, column=0)

    deltaInputFrameInner = Frame(deltaInputFrameOuter)  # , background="blue")
    deltaInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    deltaLabel = Label(deltaInputFrameInner, text="Amplitude of \u03B4(t) [\u03A6"+subscriptZero+"]:")
    deltaLabel.pack(side=LEFT, padx=(0, 5))

    lowerDeltaLabel = Label(deltaInputFrameInner, text="Lower limit:")
    lowerDeltaLabel.pack(side=LEFT)

    lowerDeltaEntry = Entry(deltaInputFrameInner, width=entryCharacterWidth, textvariable=deltaLowerNew)
    lowerDeltaEntry.pack(side=LEFT, padx=(0, 5))

    upperDeltaLabel = Label(deltaInputFrameInner, text="Upper limit:")
    upperDeltaLabel.pack(side=LEFT)

    upperDeltaEntry = Entry(deltaInputFrameInner, width=entryCharacterWidth, textvariable=deltaUpperNew)
    upperDeltaEntry.pack(side=LEFT)


def generateOmegaInputWidgets(parentWidget, height, width, entryCharacterWidth):
    omegaPhiInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="yellow")
    omegaPhiInputFrameOuter.grid(row=2, column=0)

    omegaPhiInputFrameInner = Frame(omegaPhiInputFrameOuter)  # , background="blue")
    omegaPhiInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    omegaPhiLabel = Label(omegaPhiInputFrameInner, text="Frequency \u03C9 of AC-flux [GHz]:")
    omegaPhiLabel.pack(side=LEFT, padx=(0, 5))

    lowerOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Lower limit:")
    lowerOmegaPhiLabel.pack(side=LEFT)

    lowerOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiLowerNew)
    lowerOmegaPhiEntry.pack(side=LEFT, padx=(0, 5))

    upperOmegaPhiLabel = Label(omegaPhiInputFrameInner, text="Upper limit:")
    upperOmegaPhiLabel.pack(side=LEFT)

    upperOmegaPhiEntry = Entry(omegaPhiInputFrameInner, width=entryCharacterWidth, textvariable=omegaPhiUpperNew)
    upperOmegaPhiEntry.pack(side=LEFT)


def generateModulationTimeInputWidgets(parentWidget, height, width, entryCharacterWidth):
    modulationTimeInputFrameOuter = Frame(parentWidget, height=35, width=width)# , background="green")
    modulationTimeInputFrameOuter.grid(row=4, column=0, columnspan=3)

    modulationTimeInputFrameInner = Frame(modulationTimeInputFrameOuter)  # , background="blue")
    modulationTimeInputFrameInner.place(anchor="e", relx=0.9, rely=0.5)

    modulationTimeLabel = Label(modulationTimeInputFrameInner, text="Total modulation time of AC-flux [ns]:")
    modulationTimeLabel.pack(side=LEFT, padx=(0, 5))

    lowerModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Lower limit:")
    lowerModulationTimeLabel.pack(side=LEFT)

    lowerModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeLowerNew)
    lowerModulationTimeEntry.pack(side=LEFT, padx=(0, 5))

    upperModulationTimeLabel = Label(modulationTimeInputFrameInner, text="Upper limit:")
    upperModulationTimeLabel.pack(side=LEFT)

    upperModulationTimeEntry = Entry(modulationTimeInputFrameInner, width=entryCharacterWidth, textvariable=modulationTimeUpperNew)
    upperModulationTimeEntry.pack(side=LEFT)


######################################################################################################################################################################