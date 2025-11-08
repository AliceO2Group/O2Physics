#!/bin/bash
# log file where the terminal output will be saved
STEP="deriveddata"
LOGFILE="log-${STEP}.txt"

#directory of this script
DIR_THIS=$PWD

OPTION="-b --configuration json://configuration.json"

o2-analysis-trackselection ${OPTION} |
o2-analysis-ft0-corrected-table ${OPTION} |
o2-analysis-multcenttable ${OPTION} |
o2-analysis-event-selection-service ${OPTION} |
o2-analysis-pid-tpc-service ${OPTION} |
o2-analysis-pid-tof-merge ${OPTION} |
o2-analysis-propagationservice ${OPTION} |
o2-analysis-lf-strangederivedbuilder ${OPTION} --aod-file @input_data.txt --aod-writer-json OutputDirector.json >"$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
    echo "No problems!"
    mkdir -p "${DIR_THIS}/results/${STEP}"
    mv AnalysisResults.root "${DIR_THIS}/results/${STEP}/AnalysisResults.root"
    mv AO2D.root "${DIR_THIS}/results/${STEP}/AO2D.root"
    mv dpl-config.json "${DIR_THIS}/results/${STEP}/${STEP}.json"
else
    echo "Error: Exit code ${rc}"
    echo "Check the log file ${LOGFILE}"
    exit ${rc}
fi
