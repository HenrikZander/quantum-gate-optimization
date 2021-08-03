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

######################################################################################################################################################################
# The pop up window that let's the user select a preset boundary condition.


def selectPresetWindow(root):
    pop = Toplevel(root)
    pop.title("Boundary condition manager")

    height = 500
    width = 350

    pop.geometry(str(width)+"x"+str(height))
    pop.resizable(width=False, height=False)

    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    pop.iconphoto(False, programIcon)
    
    generateBoundaryManager(pop, height, width)

    pop.grab_set()


######################################################################################################################################################################
# Function that generates the widgets in the window.


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


######################################################################################################################################################################
# Button callback functions


def selectPreset():
    pass


def addNewPreset():
    pass


def editPreset():
    pass


######################################################################################################################################################################