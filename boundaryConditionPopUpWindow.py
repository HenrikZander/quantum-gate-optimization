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

from tkinter import *
import math

import dataManager
import optimize2QubitFrame as gui

######################################################################################################################################################################
# The pop up window that let's the user select a preset boundary condition.


def selectPresetWindow(root):
    global suppliedRoot
    global pop
    suppliedRoot = root

    pop = Toplevel(root)
    pop.title("Boundary condition manager")

    height = 500
    width = 350

    pop.geometry(str(width)+"x"+str(height))
    pop.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    pop.iconphoto(False, programIcon)
    
    generateBoundaryManager(pop, height, width)
    loadAllPresets()

    pop.grab_set()


######################################################################################################################################################################
# The pop up window that is used to preview boundary presets and create new ones.


def interactWithPresetWindow(addNew=False):
    global suppliedRoot
    global popPreview
    popPreview = Toplevel(suppliedRoot)

    height = 500
    width = 400

    popPreview.geometry(str(width)+"x"+str(height))
    popPreview.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    popPreview.iconphoto(False, programIcon)

    if addNew:
        popPreview.title("Add new boundary preset")
        generateBoundaryPresetCreator(popPreview, height, width)
    else:
        popPreview.title("Edit boundary preset")
        generateBoundaryPresetEditor(popPreview, height, width)

    popPreview.grab_set()


######################################################################################################################################################################
# Function that generates the widgets in the windows.


def generateBoundaryManager(pop, height, width):
    windowTitleFrame = Frame(pop, width=width*0.85, height=height*0.05)
    windowTitleFrame.pack()

    windowTitle = Label(windowTitleFrame, text="Available presets:", font=('Helvetica', 10))
    windowTitle.place(anchor="w", relx=0, rely=0.5)

    global presetsBox
    presetsBox = Listbox(pop, width=50, height=25, relief=FLAT, selectmode=SINGLE, highlightthickness=0, bd=0, background="#E0E0E0")
    presetsBox.pack()

    buttonFrameOuter = Frame(pop, width=width*0.85, height=height*0.1)
    buttonFrameOuter.pack()

    buttonFrameInner = Frame(buttonFrameOuter)
    buttonFrameInner.place(anchor="center", relx=0.5, rely=0.5)

    selectPresetButton = Button(buttonFrameInner, text="Select", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=selectPreset)
    selectPresetButton.pack(side=LEFT, padx=4, pady=4)

    addNewPresetButton = Button(buttonFrameInner, text="Add Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=addNewPreset)
    addNewPresetButton.pack(side=LEFT, padx=4, pady=4)

    editPresetButton = Button(buttonFrameInner, text="Edit Preset", padx=3, pady=3, background="#21e4d7", relief=FLAT, command=editPreset)
    editPresetButton.pack(side=LEFT, padx=4, pady=4)


def generateBoundaryPresetCreator(popPreview, height, width):
    pass


def generateBoundaryPresetEditor(popPreview, height, width):
    pass


######################################################################################################################################################################
# Button callback functions


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
    interactWithPresetWindow(addNew=True)


def editPreset():
    interactWithPresetWindow()


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


######################################################################################################################################################################