#!/bin/bash
# log file where the terminal output will be saved
STEP="2"
LOGFILE="log-STEP${STEP}.txt"

#directory of this script
DIR_THIS=$PWD

OPTION="-b --configuration json://configuration_step2.json"

o2-analysis-trackselection ${OPTION} |
o2-analysis-ft0-corrected-table ${OPTION} |
o2-analysis-multcenttable ${OPTION} |
o2-analysis-event-selection-service ${OPTION} |
o2-analysis-pid-tpc-service ${OPTION} |
o2-analysis-propagationservice ${OPTION} |
o2-analysistutorial-lf-strangeness-step2 ${OPTION} --aod-file @input_data.txt > "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
    echo "No problems!"
    mkdir -p "${DIR_THIS}/results/step${STEP}"
    mv AnalysisResults.root "${DIR_THIS}/results/step${STEP}/AnalysisResults.root"
    mv dpl-config.json "${DIR_THIS}/results/step${STEP}/step${STEP}.json"
else
    echo "Error: Exit code ${rc}"
    echo "Check the log file ${LOGFILE}"
    exit ${rc}
fi