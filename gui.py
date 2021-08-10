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
import math
import platform
import optimize2QubitFrame
import simulate2QubitFrame

######################################################################################################################################################################
# Global variables

height = 0
width = 0

######################################################################################################################################################################
# Main 


def main():
    global root
    global height
    global width

    root = Tk()
    height, width = calculateWindowSize(root)
    root.title("GateSide")
    root.geometry(str(width)+"x"+str(math.ceil(height*1.04)))
    root.resizable(width=False, height=False)
    programIcon = PhotoImage(file = "./assets/Gateside_Logomark.png")
    root.iconphoto(False, programIcon)

    twoQubit(root, height, width)

    root.mainloop()


######################################################################################################################################################################
# Function that calculate the height of the window.


def calculateWindowSize(root):
    screen_height = root.winfo_screenheight()
    screen_width = root.winfo_screenwidth()

    operatingSystem = platform.system()
    print(operatingSystem)
    
    if operatingSystem == "Windows":
        height = math.ceil(screen_height*0.6944444)
        width = math.ceil(screen_width*0.6328125)
    elif operatingSystem == "Darwin":
        height = math.ceil(screen_height*0.9)
        width = math.ceil(screen_width*0.78)
    else:
        height = math.ceil(screen_height*0.8)
        width = math.ceil(screen_width*0.73)

    return height, width


######################################################################################################################################################################
# Functions that are used to generate different windows.


def twoQubit(root, height, width):
    twoQubitNotepad = ttk.Notebook(root, height=height, width=width)
    twoQubitNotepad.pack()

    mainMenuFrame = Frame(twoQubitNotepad, height=height, width=width , background="green")
    mainMenuFrame.place(anchor="center", relx=0.5, rely=0.5)

    optimizeFrame = optimize2QubitFrame.optimizeControlFrame(root, twoQubitNotepad, height, width)

    simulateFrame = simulate2QubitFrame.simulateControlFrame(root, twoQubitNotepad, height, width)
    
    twoQubitNotepad.add(mainMenuFrame, text="Main Menu")
    twoQubitNotepad.add(optimizeFrame, text="Optimize")
    twoQubitNotepad.add(simulateFrame, text="Simulate")


######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################