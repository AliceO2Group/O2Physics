#!/bin/bash
# log file where the terminal output will be saved
STEP="deriveddata"
LOGFILE="log-${STEP}.txt"

#directory of this script
DIR_THIS=$PWD

OPTION="-b --configuration json://configurationMC.json"

o2-analysis-pid-tof-base "${OPTION}" |
    o2-analysis-mccollisionextra "${OPTION}" |
    o2-analysis-lf-lambdakzerobuilder "${OPTION}" |
    o2-analysis-lf-cascadebuilder "${OPTION}" |
    o2-analysis-lf-cascademcbuilder "${OPTION}" |
    o2-analysis-centrality-table "${OPTION}" |
    o2-analysis-lf-lambdakzeromcbuilder "${OPTION}" |
    o2-analysis-mccollision-converter "${OPTION}" |
    o2-analysis-ud-sgcand-producer "${OPTION}" |
    o2-analysis-timestamp "${OPTION}" |
    o2-analysis-ft0-corrected-table "${OPTION}" |
    o2-analysis-track-propagation "${OPTION}" |
    o2-analysis-pid-tpc-base "${OPTION}" |
    o2-analysis-pid-tpc "${OPTION}" |
    o2-analysis-multiplicity-table "${OPTION}" |
    o2-analysis-trackselection "${OPTION}" |
    o2-analysis-pid-tof-full "${OPTION}" |
    o2-analysis-pid-tof-beta "${OPTION}" |
    o2-analysis-event-selection "${OPTION}" |
    o2-analysis-lf-strangederivedbuilder "${OPTION}" --aod-file @input_dataMC.txt --aod-writer-json OutputDirectorMC.json >"$LOGFILE" 2>&1

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
