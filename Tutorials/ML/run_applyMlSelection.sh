#!/bin/bash

# Copyright 2019-2023 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

# @brief Bash script to execute the ML selection task on Run 3 MC input
#
# The input AO2D.root files is expected in the working directory.
#
# @author Vít Kučera <vit.kucera@cern.ch>, Inha University
# @author Fabio Catalano <fabio.catalano@cern.ch>, CERN
# @date 2023-10-29

# log file where the terminal output will be saved
LOGFILE="stdout.log"

# directory of this script
DIR_THIS="$(dirname "$(realpath "$0")")"

# O2 configuration file (in the same directory)
JSON="$DIR_THIS/dpl-config_applyMlSelection.json"

# command line options of O2 workflows
OPTIONS="-b --configuration json://$JSON --aod-memory-rate-limit 2000000000 --shm-segment-size 16000000000 --resources-monitoring 2 --min-failure-level error"

# execute the mini task workflow and its dependencies
# shellcheck disable=SC2086 # Ignore unquoted options.
o2-analysistutorial-apply-ml-selection $OPTIONS | \
o2-analysis-hf-candidate-creator-3prong $OPTIONS | \
o2-analysis-hf-track-index-skim-creator $OPTIONS | \
o2-analysis-track-to-collision-associator $OPTIONS | \
o2-analysis-timestamp $OPTIONS | \
o2-analysis-track-propagation $OPTIONS | \
o2-analysis-event-selection $OPTIONS | \
o2-analysis-trackselection $OPTIONS \
> "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi
