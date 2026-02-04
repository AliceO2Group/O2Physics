#!/bin/bash
# log file where the terminal output will be saved
STEP="skeleton"
LOGFILE="log-${STEP}.txt"

#directory of this script
DIR_THIS=$PWD

OPTION="-b --configuration json://configuration_skeleton.json"

o2-analysistutorial-lf-strangeness-derived-skeleton ${OPTION} --aod-file @input_data.txt > "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
    echo "No problems!"
    mkdir -p "${DIR_THIS}/results/${STEP}"
    mv AnalysisResults.root "${DIR_THIS}/results/${STEP}/AnalysisResults.root"
    mv dpl-config.json "${DIR_THIS}/results/${STEP}/${STEP}.json"
else
    echo "Error: Exit code ${rc}"
    echo "Check the log file ${LOGFILE}"
    exit ${rc}
fi