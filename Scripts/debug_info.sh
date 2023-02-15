#!/bin/bash

# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# @brief Extract useful debugging information from an O2 log file.
#
# This script extracts the following information from an O2 log file:
# - devices (omitting internal-dpl...)
# - input files
# - warnings
# - errors
#
# @author Vít Kučera <vit.kucera@cern.ch>, Inha University
# @date 2023-02-15

LOG=""
DEVICES=0
INPUT=0
WARNINGS=0
ERRORS=0

# Print devices.
PrintDevices() {
  echo -e "\nDevices:"
  grep "\\[INFO\\] Starting " "$LOG" | grep -v "internal-dpl" | cut -d" " -f3 | sort
}

# Print input files.
PrintInputFiles() {
  echo -e "\nInput files:"
  grep "\\[INFO\\] Read info: " "$LOG" | cut -d= -f2 | cut -d, -f1 | sort
}

# Print warnings.
PrintWarnings() {
  echo -e "\nWarnings:"
  grep \
  -e "\\[WARN\\]" \
  -e "Warning in " \
  "$LOG" | sort -u
}

# Print errors.
PrintErrors() {
  echo -e "\nErrors:"
  grep \
  -e "\\[ERROR\\]" \
  -e "\\[FATAL\\]" \
  -e "segmentation" \
  -e "Segmentation" \
  -e "command not found" \
  -e "Error:" \
  -e "Error in " \
  "$LOG" | sort -u
}

# Print out a help message.
Help() {
  echo "Usage: $(basename "$0") -f LOGFILE [-d] [-i] [-w] [-e] [-h] "
  echo "-d  Show devices."
  echo "-i  Show input files."
  echo "-w  Show warnings."
  echo "-e  Show errors."
}

####################################################################################################

# Parse command line options.
while getopts ":hf:diwe" opt; do
  case ${opt} in
    h)
      Help; exit 0;;
    f)
      LOG="$OPTARG";;
    d)
      DEVICES=1;;
    i)
      INPUT=1;;
    w)
      WARNINGS=1;;
    e)
      ERRORS=1;;
    \?)
      echo "Error: Invalid option: $OPTARG" 1>&2; Help; exit 1;;
    :)
      echo "Error: Invalid option: $OPTARG requires an argument." 1>&2; Help; exit 1;;
  esac
done

# Check a valid input file.
[ "$LOG" ] || { echo "Error: Provide a log file." 1>&2; exit 1; }
[ -f "$LOG" ] || { echo "Error: File $LOG does not exist." 1>&2; exit 1; }

[ $DEVICES -eq 1 ] && PrintDevices
[ $INPUT -eq 1 ] && PrintInputFiles
[ $WARNINGS -eq 1 ] && PrintWarnings
[ $ERRORS -eq 1 ] && PrintErrors

exit 0
