######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: GateSide.py

# Author(s): Henrik Zander

# Date created: 9 January 2022

# Copyright 2022, Henrik Zander, All rights reserved.

######################################################################################################################################################################

# Global Python imports
from pathlib import Path

# Local Python imports
import report_generation.interface as report_generator
import graph_generation.interface as graph_generator

######################################################################################################################################################################

path = "C:\\Users\\henri\\Desktop\\WACQT\\quantum-gate-optimization\\Results\\Qubit Pair 01\\Solutions\\210728_iSWAP_1.json"
# 210728_iSWAP_1.json

graph_generator.plotFidelity(path, 'temp\\fig1.jpeg')

report_generator.generateReport2(path)