from tkinter import *
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time

from main import *
from optimize2QubitWindow import *

######################################################################################################################################################################
# Global variables and function to generate global variables

height = 600
width = int(1.62*height)
cancelOptimization = False
relativeHeight = 0.85
relativeWidth = 1

######################################################################################################################################################################
# Main 

def main():
    global root

    root = Tk()
    generateGlobalVariables(root)
    root.title("GateSide")
    root.geometry(str(width)+"x"+str(height))
    root.resizable(width=False, height=False)
    programIcon = PhotoImage(file = "../GUI/Gateside_Logomark.png")
    root.iconphoto(False, programIcon)

    optimizeControlWindow(root)

    root.mainloop()

######################################################################################################################################################################
# Generate different windows

def optimizeControlWindow(root):
    topFrame = Frame(root, height=relativeHeight*height, width=relativeWidth*width, background="green")
    topFrame.grid(column=0, row=0)
    topFrame.grid_propagate(0)

    bottomFrame = Frame(root, height=(1-relativeHeight)*height, width=relativeWidth*width)
    bottomFrame.grid(column=0, row=1)
    bottomFrame.grid_propagate(0)

    optimizeStatusFrame(bottomFrame)
    optimizeCircuitParameterFrame(topFrame)
    optimizerSettingsFrame(topFrame)

######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################