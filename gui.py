######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: gui.py

# Author(s): Henrik Zander

# Date created: 5 July 2021

# Copyright 2021, Henrik Zander, All rights reserved.

######################################################################################################################################################################

from tkinter import *
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time

from main import *
import optimize2QubitWindow

######################################################################################################################################################################
# Global variables

height = 600
width = int(1.62*height)

######################################################################################################################################################################
# Main 

def main():
    global root

    root = Tk()
    root.title("GateSide")
    root.geometry(str(width)+"x"+str(height))
    root.resizable(width=False, height=False)
    programIcon = PhotoImage(file = "../GUI/Gateside_Logomark.png")
    root.iconphoto(False, programIcon)

    optimizeControlWindow(root, height)

    root.mainloop()

######################################################################################################################################################################
# Functions that is used to manage the config.JSON-file.



######################################################################################################################################################################
# Functions that generate different windows

def optimizeControlWindow(root, height):
    optimize2QubitWindow.initiateGlobalVariables(root, height)
    relativeHeight = optimize2QubitWindow.relativeHeight
    relativeWidth = optimize2QubitWindow.relativeWidth

    topFrame = Frame(root, height=relativeHeight*height, width=relativeWidth*width)#, background="green")
    topFrame.grid(row=0, column=0)
    topFrame.grid_propagate(0)

    separator = ttk.Separator(root, orient='horizontal')
    separator.grid(row=1, column=0, sticky="nesw")

    bottomFrame = Frame(root, height=(1-relativeHeight-0.02)*height, width=relativeWidth*width)
    bottomFrame.grid(row=2, column=0)
    bottomFrame.grid_propagate(0)

    optimize2QubitWindow.optimizeStatusFrame(bottomFrame)
    optimize2QubitWindow.optimizeCircuitParameterFrame(topFrame)
    optimize2QubitWindow.optimizerSettingsFrame(topFrame)

    optimize2QubitWindow.useDefaultCircuit()
    optimize2QubitWindow.setDefaultBoundaryValues()

######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################