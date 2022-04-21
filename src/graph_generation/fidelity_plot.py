######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: fidelity_plot.py

# Author(s): Henrik Zander

# Date created: 21 April 2022

# Copyright 2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports.
import matplotlib.pyplot as plt

# Local Python imports.

######################################################################################################################################################################
# Function for creating a fidelity plot of a certain gate.

def generate(times, fidelities, modulationTime, saveLocationPath, saveGraph=True):
        minimumFidelity = max(fidelities)

        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot()

        ax.plot(times, fidelities)
        ax.plot([modulationTime, modulationTime], [0, 1], 'r--')

        ax.grid()
        ax.set_ylim([minimumFidelity*0.996, 1])
        ax.set_xlim([times[0], times[-1]])
        leg = ax.legend(["Fidelity", "$t_{MOD}$"], fontsize=20, loc="lower right")

        ax.set_xlabel("Time since gate start [ns]", fontsize=23)
        ax.set_ylabel("Fidelity", fontsize=23)
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(axis='both', which='minor', labelsize=16)

        for legobj in leg.legendHandles:
            legobj.set_linewidth(2.0)

        # plt.title("Gate fidelity around $t_{MOD}$", fontsize=18)
        plt.tight_layout()

        if saveGraph:
            plt.savefig(saveLocationPath)
        else:
            plt.show()

######################################################################################################################################################################
