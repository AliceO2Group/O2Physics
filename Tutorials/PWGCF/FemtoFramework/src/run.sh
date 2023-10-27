#!/bin/bash

# in case the central server is not reachable, set SE yourself
# https://alimonitor.cern.ch/stats?page=SE/table
# export alien_CLOSE_SE=ALICE::UPB::EOS

# make sure you can connect to central server
# alien.py

# number of task you want to run, run task 1 by default
Task="${1:-1}"

# config file to use, use config.json by default
ConfigFile="${2:-config.json}"

# options passed to each workflow
Options=("-b" "--configuration" "json://${ConfigFile}")

# command to be executed
Command="o2-analysistutorial-cf-femtodream-tutorial-${Task} ${Options[*]}"

# print comand before executing it
echo "$Command"
eval "$Command"

exit 0
