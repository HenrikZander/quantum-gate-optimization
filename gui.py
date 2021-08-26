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
from PIL import ImageTk,Image

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

    global programNotepad
    programNotepad = ttk.Notebook(root, height=height, width=width)
    programNotepad.pack()

    mainMenu(programNotepad, height=height, width=width)
    # twoQubit(root, programNotepad, height, width)

    root.mainloop()


######################################################################################################################################################################
# Button callback functions.


def openTwoQubitWindow():
    closeAllTabs()
    twoQubit(root, programNotepad, height, width)


def openThreeQubitWindow():
    closeAllTabs()


######################################################################################################################################################################
# Helper functions for the window generation.


def calculateWindowSize(root):
    screen_height = root.winfo_screenheight()
    screen_width = root.winfo_screenwidth()

    operatingSystem = platform.system()

    if operatingSystem == "Windows":
        height = math.ceil(screen_height*0.6944444)
        width = math.ceil(screen_width*0.6328125)
    elif operatingSystem == "Darwin":
        height = math.ceil(screen_height*0.75)
        width = math.ceil(screen_width*0.75)
    else:
        height = math.ceil(screen_height*0.8)
        width = math.ceil(screen_width*0.73)

    return height, width


def closeAllTabs():
    global programNotepad
    tabs = programNotepad.tabs()

    for tab in tabs[1:]:
        programNotepad.forget(tab)


######################################################################################################################################################################
# Functions that are used to generate different windows.


def mainMenu(parentWidget, height=height, width=width):
    mainMenuFrame = Frame(parentWidget, height=height, width=width)
    mainMenuFrame.place(anchor="center", relx=0.5, rely=0.5)

    parentWidget.add(mainMenuFrame, text="Main Menu")

    ####################################################################################
    # Generation of actual main menu widgets.
    ####################################################################################

    mainMenuLogo = ImageTk.PhotoImage(Image.open("./assets/gateside_logo.png").resize((400,100)))
    mainMenuLogoLabel = Label(mainMenuFrame, image=mainMenuLogo)
    mainMenuLogoLabel.image = mainMenuLogo
    mainMenuLogoLabel.pack(pady=30)

    twoQubitButton = Button(mainMenuFrame, text="Two Qubit Gates", command=openTwoQubitWindow, padx=3, pady=3, background="#21e4d7", relief=FLAT, font=('Helvetica', 12))
    twoQubitButton.pack(pady=10)

    threeQubitButton = Button(mainMenuFrame, text="Three Qubit Gates", command=openThreeQubitWindow, padx=3, pady=3, background="#21e4d7", relief=FLAT, font=('Helvetica', 12))
    threeQubitButton.pack(pady=10)

    copyrightLabel = Label(mainMenuFrame, text="Copyright © Henrik Zander, Emil Ingelsten 2021; All rights reserved.")
    copyrightLabel.pack(side=BOTTOM)

    ####################################################################################
    return mainMenuFrame


def twoQubit(root, notepad, height, width):
    optimizeFrame = optimize2QubitFrame.optimizeControlFrame(root, notepad, height, width)
    simulateFrame = simulate2QubitFrame.simulateControlFrame(root, notepad, height, width)
    
    notepad.add(optimizeFrame, text="Optimize")
    notepad.add(simulateFrame, text="Simulate")


######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################