from tkinter import *
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time

from main import *
from optimize2QubitWindow import *

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
    initiateGlobalVariables(root, height)

    topFrame = Frame(root, height=relativeHeight*height, width=relativeWidth*width)#, background="green")
    topFrame.grid(row=0, column=0)
    topFrame.grid_propagate(0)

    separator = ttk.Separator(root, orient='horizontal')
    separator.grid(row=1, column=0, sticky="nesw")

    bottomFrame = Frame(root, height=(1-relativeHeight-0.02)*height, width=relativeWidth*width)
    bottomFrame.grid(row=2, column=0)
    bottomFrame.grid_propagate(0)

    optimizeStatusFrame(bottomFrame)
    optimizeCircuitParameterFrame(topFrame)
    optimizerSettingsFrame(topFrame)

######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################