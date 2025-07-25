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
parser.add_argument('--add_mc_conv', help="Add the converter from mcparticle to mcparticle+001", action="store_true")
parser.add_argument('--add_fdd_conv', help="Add the fdd converter", action="store_true")
parser.add_argument('--add_track_prop', help="Add track propagation to the innermost layer (TPC or ITS)", action="store_true")
parser.add_argument("--add_weakdecay_ind", help = "Add Converts V0 and cascade version 000 to 001", action = "store_true")
parser.add_argument("--add_col_conv", help = "Add the converter from collision to collision+001", action = "store_true")
extrargs = parser.parse_args()

commonDeps = ["o2-analysis-timestamp", "o2-analysis-event-selection", "o2-analysis-multiplicity-table"]
barrelDeps = ["o2-analysis-trackselection", "o2-analysis-trackextension","o2-analysis-pid-tof-base", "o2-analysis-pid-tof", "o2-analysis-pid-tof-full", "o2-analysis-pid-tof-beta", "o2-analysis-pid-tpc-full"]
#barrelDeps = ["o2-analysis-trackselection","o2-analysis-pid-tof-base", "o2-analysis-pid-tof", "o2-analysis-pid-tof-full", "o2-analysis-pid-tof-beta", "o2-analysis-pid-tpc-full"]
muonDeps = ["o2-analysis-fwdtrackextension"]
specificDeps = {
  "processFull" : [],
  "processFullWithCov" : [],
  "processFullWithCovAndEventFilter" : ["o2-analysis-dq-filter-pp","o2-analysis-fwdtrackextension"],
  "processFullWithCent" : ["o2-analysis-centrality-table"],
  "processFullWithCentAndMults" : ["o2-analysis-centrality-table"],
  "processBarrelOnly" : [],
  "processBarrelOnlyWithCov" : [],
  "processBarrelOnlyWithV0Bits" : ["o2-analysis-dq-v0-selector"],
  "processBarrelOnlyWithDalitzBits" : ["o2-analysis-dq-dalitz-selection"],
  "processBarrelOnlyWithEventFilter" : ["o2-analysis-dq-filter-pp","o2-analysis-fwdtrackextension"],
  "processBarrelOnlyWithMults" : [],
  "processBarrelOnlyWithCovAndEventFilter" : ["o2-analysis-dq-filter-pp","o2-analysis-fwdtrackextension"],
  "processBarrelOnlyWithQvector" : ["o2-analysis-centrality-table", "o2-analysis-dq-flow"],
  "processBarrelOnlyWithCent" : ["o2-analysis-centrality-table"],
  "processBarrelOnlyWithCentAndMults" : ["o2-analysis-centrality-table"],
  "processMuonOnly" : [],
  "processMuonOnlyWithCov" : [],
  "processMuonOnlyWithCent" : ["o2-analysis-centrality-table"],
  "processMuonOnlyWithMults" : [],
  "processMuonOnlyWithCentAndMults" : ["o2-analysis-centrality-table"],
  "processMuonOnlyWithCovAndCent" : ["o2-analysis-centrality-table"],
  "processMuonOnlyWithQvector" : ["o2-analysis-centrality-table", "o2-analysis-dq-flow"],
  "processMuonOnlyWithFilter" : ["o2-analysis-dq-filter-pp"],
  "processAmbiguousMuonOnly" : [],
  "processAmbiguousMuonOnlyWithCov" : [],
  "processAmbiguousBarrelOnly" : []
  # "processFullWithCentWithV0Bits": ["o2-analysis-centrality-table","o2-analysis-dq-v0-selector", "o2-analysis-weak-decay-indices"],
  # "processFullWithEventFilterWithV0Bits": ["o2-analysis-dq-filter-pp","o2-analysis-dq-v0-selector", "o2-analysis-weak-decay-indices"],
}

# Definition of all the tables we may write
tables = {
  "ReducedEvents" : {"table": "AOD/REDUCEDEVENT/0"},
  "ReducedEventsExtended" : {"table": "AOD/REEXTENDED/0"},
  "ReducedEventsVtxCov" : {"table": "AOD/REVTXCOV/0"},
  "ReducedEventsQvector" : {"table": "AOD/REQVECTOR/0"},
  "ReducedMCEventLabels" : {"table": "AOD/REMCCOLLBL/0"},
  "ReducedMCEvents" : {"table": "AOD/REDUCEDMCEVENT/0"},
  "ReducedTracks" : {"table": "AOD/REDUCEDTRACK/0"},
  "ReducedTracksBarrel" : {"table": "AOD/RTBARREL/0"},
  "ReducedTracksBarrelCov" : {"table": "AOD/RTBARRELCOV/0"},
  "ReducedTracksBarrelPID" : {"table": "AOD/RTBARRELPID/0"},
  "ReducedTracksBarrelLabels" : {"table": "AOD/RTBARRELLABELS/0"},
  "ReducedMCTracks" : {"table": "AOD/REDUCEDMCTRACK/0"},
  "ReducedMuons" : {"table": "AOD/RTMUON/0"},
  "ReducedMuonsExtra" : {"table": "AOD/RTMUONEXTRA/0"},
  "ReducedMuonsCov" : {"table": "AOD/RTMUONCOV/0"},
  "ReducedMuonsLabels" : {"table": "AOD/RTMUONSLABELS/0"},
  "AmbiguousTracksMid" : {"table": "AOD/AMBIGUOUSTRACK/0"},
  "AmbiguousTracksFwd" : {"table": "AOD/AMBIGUOUSFWDTR/0"},
  "DalitzBits" : {"table": "AOD/DALITZBITS/0"}
}
# Tables to be written, per process function
commonTables = ["ReducedEvents", "ReducedEventsExtended", "ReducedEventsVtxCov"]
barrelCommonTables = ["ReducedTracks","ReducedTracksBarrel","ReducedTracksBarrelPID"]
muonCommonTables = ["ReducedMuons", "ReducedMuonsExtra"]
specificTables = {
  "processFull": [],
  "processFullWithCov": ["ReducedTracksBarrelCov", "ReducedMuonsCov"],
  "processFullWithCovAndEventFilter": ["ReducedTracksBarrelCov", "ReducedMuonsCov"],
  "processFullWithCent": [],
  "processFullWithCentAndMults": [],
  "processBarrelOnly": [],
  "processBarrelOnlyWithCov": ["ReducedTracksBarrelCov"],
  "processBarrelOnlyWithV0Bits": [],
  "processBarrelOnlyWithDalitzBits": ["DalitzBits"],
  "processBarrelOnlyWithQvector": ["ReducedEventsQvector"],
  "processBarrelOnlyWithEventFilter": [],
  "processBarrelOnlyWithMults": [],
  "processBarrelOnlyWithCovAndEventFilter": ["ReducedTracksBarrelCov"],
  "processBarrelOnlyWithCent": [],
  "processBarrelOnlyWithCentAndMults": [],
  "processMuonOnly": [],
  "processMuonOnlyWithCov": ["ReducedMuonsCov"],
  "processMuonOnlyWithCent": [],
  "processMuonOnlyWithMults": [],
  "processMuonOnlyWithCentAndMults": [],
  "processMuonOnlyWithCovAndCent": ["ReducedMuonsCov"],
  "processMuonOnlyWithQvector": ["ReducedEventsQvector"],
  "processMuonOnlyWithFilter": [],
  "processAmbiguousMuonOnly": ["AmbiguousTracksFwd"],
  "processAmbiguousMuonOnlyWithCov": ["AmbiguousTracksFwd", "ReducedMuonsCov"],
  "processAmbiguousBarrelOnly": ["AmbiguousTracksMid"]
}

# Make some checks on provided arguments
if len(sys.argv) < 3:
  print("ERROR: Invalid syntax! The command line should look like this:")
  print("  ./runTableMaker.py <yourConfig.json> <runData|runMC> [task:param:value] ...")
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

# Delete trackextension dependency if track-propagation dependency provided (for compatibility)
if extrargs.add_track_prop:
  barrelDeps.remove("o2-analysis-trackextension")

if extrargs.arg != "" and extrargs.arg is not None:
  args = [line.split(':') for line in extrargs.arg.split(',') if line]
  for threeIndex in args:
    if len(threeIndex) != 3:
      print("ERROR: Wrong parameter syntax for --arg: ", threeIndex, " in ", extrargs.arg)
      print("Correct syntax: task:param:value,task:param:value ... ")
      print("Example: --arg table-maker:processBarrelOnly:true")
      sys.exit()
  for arg in args:
    config[arg[0]][arg[1]] = arg[2]

taskNameInConfig = "table-maker"
taskNameInCommandLine = "o2-analysis-dq-table-maker"
if runOverMC == True:
  taskNameInConfig = "table-maker-m-c"
  taskNameInCommandLine = "o2-analysis-dq-table-maker-mc"

if not taskNameInConfig in config:
  print("ERROR: Task to be run not found in the configuration file!")
  sys.exit()

# Write the updated configuration file into a temporary file
updatedConfigFileName = "tempConfig.json"
with open(updatedConfigFileName,'w') as outputFile:
  json.dump(config, outputFile, indent = 2)

# Check which dependencies need to be run
depsToRun = {}
for dep in commonDeps:
  depsToRun[dep] = 1

for processFunc in specificDeps.keys():
  if not processFunc in config[taskNameInConfig].keys():
    continue
  if config[taskNameInConfig][processFunc] == "true":
    if "processFull" in processFunc or "processBarrel" in processFunc or "processAmbiguousBarrel" in processFunc:
      for dep in barrelDeps:
        depsToRun[dep] = 1
    if "processFull" in processFunc or "processMuon" in processFunc or "processAmbiguousMuon" in processFunc:
      for dep in muonDeps:
        depsToRun[dep] = 1
    for dep in specificDeps[processFunc]:
      depsToRun[dep] = 1

# Check which tables are required in the output
tablesToProduce = {}
for table in commonTables:
  tablesToProduce[table] = 1

if runOverMC == True:
  tablesToProduce["ReducedMCEvents"] = 1
  tablesToProduce["ReducedMCEventLabels"] = 1

for processFunc in specificDeps.keys():
  if not processFunc in config[taskNameInConfig].keys():
    continue
  if config[taskNameInConfig][processFunc] == "true":
    print("processFunc ========")
    print(processFunc)
    if "processFull" in processFunc or "processBarrel" in processFunc:
      print("common barrel tables==========")
      for table in barrelCommonTables:
        print(table)
        tablesToProduce[table] = 1
      if runOverMC == True:
        tablesToProduce["ReducedTracksBarrelLabels"] = 1
    if "processFull" in processFunc or "processMuon" in processFunc:
      print("common muon tables==========")
      for table in muonCommonTables:
        print(table)
        tablesToProduce[table] = 1
      if runOverMC == True:
        tablesToProduce["ReducedMuonsLabels"] = 1
    if runOverMC == True:
      tablesToProduce["ReducedMCTracks"] = 1
    print("specific tables==========")
    for table in specificTables[processFunc]:
      print(table)
      tablesToProduce[table] = 1

# Generate the aod-writer output descriptor json file
writerConfig = {}
writerConfig["OutputDirector"] = {
  "debugmode": True,
  "resfile": "reducedAod",
  "resfilemode": "RECREATE",
  "ntfmerge": 1,
  "OutputDescriptors": []
}
iTable = 0
for table in tablesToProduce.keys():
  writerConfig["OutputDirector"]["OutputDescriptors"].insert(iTable, tables[table])
  iTable += 1

writerConfigFileName = "aodWriterTempConfig.json"
with open(writerConfigFileName,'w') as writerConfigFile:
  json.dump(writerConfig, writerConfigFile, indent = 2)

print(writerConfig)
#sys.exit()

commandToRun = taskNameInCommandLine + " --configuration json://" + updatedConfigFileName + " --severity error --shm-segment-size 12000000000 --aod-writer-json " + writerConfigFileName + " -b"
for dep in depsToRun.keys():
  commandToRun += " | " + dep + " --configuration json://" + updatedConfigFileName + " -b"

if extrargs.add_mc_conv:
  commandToRun += " | o2-analysis-mc-converter --configuration json://" + updatedConfigFileName + " -b"

if extrargs.add_fdd_conv:
  commandToRun += " | o2-analysis-fdd-converter --configuration json://" + updatedConfigFileName + " -b"

if extrargs.add_track_prop:
  commandToRun += " | o2-analysis-track-propagation --configuration json://" + updatedConfigFileName + " -b"

if extrargs.add_weakdecay_ind:
  commandToRun += " | o2-analysis-weak-decay-indices --configuration json://" + updatedConfigFileName + " -b"

if extrargs.add_col_conv:
  commandToRun += " | o2-analysis-collision-converter --configuration json://" + updatedConfigFileName + " -b"

print("====================================================================================================================")
print("Command to run:")
print(commandToRun)
print("====================================================================================================================")
print("Tables to produce:")
print(tablesToProduce.keys())
print("====================================================================================================================")
#sys.exit()

os.system(commandToRun)

