######################################################################################################################################################################

#      .oooooo.   oooo                  oooo
#     d8P'  `Y8b  `888                  `888
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P'


# File name: dataManager.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 19 July 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

import json

######################################################################################################################################################################
# Functions that handle data loading and dumping to json-files.


def getFromjson(fileName):
    with open(fileName, 'r') as jsonFile:
        jsonDict = json.load(jsonFile)
    return jsonDict


def dumpTojson(data, fileName):
    with open(fileName, 'w') as jsonFile:
        json.dump(data, jsonFile, ensure_ascii=False, indent=4)


######################################################################################################################################################################