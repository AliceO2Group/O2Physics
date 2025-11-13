#!/bin/bash

OPTION="-b --configuration json://myConfigExtractor.json"

# Specify the executable binary directly (adjust the path as needed)
EXECUTABLE=~/alice/sw/BUILD/O2Physics-latest/O2Physics/stage/bin/o2-analysistutorial-mm-my-example-task-pid-feature-extractor
CONFIG_PATH="$(pwd)/myConfigExtractor.json"
OPTION="-b --configuration json://$CONFIG_PATH"


echo "Starting O2Physics PID Feature Extraction Workflow..."
echo "Using configuration: myConfigExtractor.json"

o2-analysis-timestamp ${OPTION} | \
o2-analysis-event-selection ${OPTION} | \
o2-analysis-tracks-extra-v002-converter ${OPTION} | \
o2-analysis-track-propagation ${OPTION} | \
o2-analysis-pid-tpc-base ${OPTION} | \
o2-analysis-pid-tpc ${OPTION} | \
o2-analysis-pid-tof-base ${OPTION} | \
o2-analysis-pid-tof ${OPTION} | \
o2-analysis-pid-tof-beta ${OPTION} | \
o2-analysis-multiplicity-table ${OPTION} | \
o2-analysis-mccollision-converter ${OPTION} | \
${EXECUTABLE} ${OPTION}

echo "Check output files (e.g. pid_features.csv, pid_features.root)"
