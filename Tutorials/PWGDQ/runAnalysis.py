#!/usr/bin/env python3

from ast import parse
import sys
import json
import os
import argparse

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='config.json', help='config file name')
parser.add_argument('-runData', help="Run over data", action="store_true")
parser.add_argument('-runMC', help="Run over MC", action="store_true")
parser.add_argument('--arg', help='Configuration argument')
parser.add_argument("--aod-writer-json", help = "Name of the json configuration file", action = "store", type = str)

extrargs = parser.parse_args()

# Make some checks on provided arguments
if len(sys.argv) < 3:
  print("ERROR: Invalid syntax! The command line should look like this:")
  print("  ./runAnalysis.py <yourConfig.json> <runData|runMC> [task:param:value] ...")
  sys.exit()

# Load the configuration file provided as the first parameter
config = {}
with open(extrargs.cfgFileName) as configFile:
  config = json.load(configFile)

# Check whether we run over data or MC
if not (extrargs.runMC or extrargs.runData):
  print("ERROR: You have to specify either runMC or runData !")
  sys.exit()

runOverMC = False
if (extrargs.runMC):
  runOverMC = True

print("runOverMC ",runOverMC)

if extrargs.arg != "" and extrargs.arg is not None:
  args = [line.split(':') for line in extrargs.arg.split(',') if line]
  for threeIndex in args:
    if len(threeIndex) != 3:
      print("ERROR: Wrong parameter syntax for --arg: ", threeIndex, " in ", extrargs.arg)
      print("Correct syntax: task:param:value,task:param:value ... ")
      print("Example: --arg analysis-same-event-pairing:processDecayToEESkimmed:true")
      sys.exit()
  for arg in args:
    config[arg[0]][arg[1]] = arg[2]

taskNameInConfig = "analysis-event-selection"
taskNameInCommandLine = "o2-analysis-dq-table-reader"
if runOverMC == True:
  taskNameInCommandLine = "o2-analysis-dq-efficiency"

if not taskNameInConfig in config:
  print("ERROR: Task to be run not found in the configuration file!")
  sys.exit()

# Write the updated configuration file into a temporary file
updatedConfigFileName = "tempConfig.json"
with open(updatedConfigFileName,'w') as outputFile:
  json.dump(config, outputFile, indent = 2)
  
commandToRun = (taskNameInCommandLine + " --configuration json://" + updatedConfigFileName + " -b")
if extrargs.aod_writer_json:
    commandToRun = (taskNameInCommandLine + " --configuration json://" + updatedConfigFileName + " --aod-writer-json " + extrargs.aod_writer_json + " -b")
    
print("====================================================================================================================")
print("Command to run:")
print(commandToRun)
print("====================================================================================================================")
os.system(commandToRun)
