from tkinter import *
from tkinter.ttk import *
import time

height = 500
width = int(1.62*height)
cancelOptimization = False

def main():
    root = Tk()
    root.title("GateSide by WACQT")
    root.geometry(str(width)+"x"+str(height))

    optimizeControlWindow(root)

    root.mainloop()
######################################################################################################################################################################
# Different windows

def startOptimizing():
    global cancelOptimization

    while not cancelOptimization:
        print("Doing stuff!")
        time.sleep(2)


def stopOptimizing():
    global cancelOptimization
    cancelOptimization = True



def optimizeControlWindow(root, relativeHeight=0.75, relativeWidth=1):
    
    style = Style()
    style.configure('w1.TFrame',background="#00F")

    topFrame = Frame(root, style='w1.TFrame', height=relativeHeight*height, width=relativeWidth*width)
    topFrame.grid(column=0, row=0)
    topFrame.grid_propagate(0)

    style = Style()
    style.configure('w2.TFrame',background="#0F0")

    bottomFrame = Frame(root, style='w2.TFrame', height=(1-relativeHeight)*height, width=relativeWidth*width)
    bottomFrame.grid(column=0, row=1)
    bottomFrame.grid_propagate(0)

    global startOptimizeButton 
    startOptimizeButton = Button(bottomFrame, text="Start Optimizing", command=startOptimizing)
    startOptimizeButton.grid()

    global stopOptimizeButton 
    stopOptimizeButton = Button(bottomFrame, text="Stop Optimizing", command=stopOptimizing)
    stopOptimizeButton.grid(column=0, row=1)








######################################################################################################################################################################
# Auto-run the main()-function if main.py is compiled.

if __name__ == "__main__":
    main()

######################################################################################################################################################################