from tkinter import *
import tkinter.ttk as ttk
import tkinter.font as tkFont
import time

######################################################################################################################################################################
# Global variables

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
    root.title("GateSide")
    root.geometry(str(width)+"x"+str(height))
    root.resizable(width=False, height=False)
    programIcon = PhotoImage(file = "../GUI/gateside_icon.png")
    root.iconphoto(False, programIcon)

    optimizeControlWindow(root)

    root.mainloop()

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
    print("Circuit loaded.")


def setDefaultCircuit():
    print("Default Circuit set.")


def useDefaultCircuit():
    print("Reset to default circuit.")


def exportCircuit():
    print("Circuit exported to JSON.")

######################################################################################################################################################################
# Window for optimization control of two qubit gates.

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
    setDefaultCircuitButton = Button(controlsInputFrameInner, text="Set Default", command=setDefaultCircuit, padx=3, pady=3, background="#21e4d7")
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

    gateTitleFrame = Frame(circuitFrame, height=50, width=relativeWidth*width*0.40)#, background="orange")
    gateTitleFrame.grid(row=0, column=0, columnspan=3)

    gateTitle = Label(gateTitleFrame, text='Circuit parameters', font=('Helvetica',11))
    titleFont = tkFont.Font(gateTitle, gateTitle.cget("font"))
    titleFont.configure(underline = True, size=12)
    gateTitle.configure(font=titleFont)
    gateTitle.place(anchor='center', relx=0.5, rely=0.5)

    generateFrequencyInputs(circuitFrame, entryCharacterWidth)
    generateAnharmonicityInputs(circuitFrame, entryCharacterWidth)
    generateCouplingInputs(circuitFrame, entryCharacterWidth)


def optimizeCircuitParameterFrame(topFrame):
    circuitFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.40)#, highlightbackground="black", highlightthickness=1, background='blue')
    circuitFrame.grid(column=0, row=0)
    circuitFrame.grid_propagate(0)

    generateCircuitInputs(circuitFrame)
    generateCircuitInputControls(circuitFrame)


def optimizerSettingsFrame(topFrame):
    settingsFrame = Frame(topFrame, height=relativeHeight*height, width=relativeWidth*width*0.60)
    settingsFrame.grid(column=1, row=0)
    settingsFrame.grid_propagate(0)


def optimizeControlWindow(root):
    topFrame = Frame(root, height=relativeHeight*height, width=relativeWidth*width, background="green")
    topFrame.grid(column=0, row=0)
    topFrame.grid_propagate(0)

    bottomFrame = Frame(root, height=(1-relativeHeight)*height, width=relativeWidth*width)
    bottomFrame.grid(column=0, row=1)
    bottomFrame.grid_propagate(0)

    optimizeStatusFrame(bottomFrame)
    optimizeCircuitParameterFrame(topFrame)
    #optimizerSettingsFrame(topFrame)

######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################