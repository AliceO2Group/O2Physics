#!/bin/bash
# log file where the terminal output will be saved
STEP="skeleton"
LOGFILE="log-${STEP}.txt"

# directory of this script
DIR_THIS="$(dirname "$(realpath "$0")")"

# run from the folder of this script, so that the relative paths in input_data.txt work
cd "$DIR_THIS" || exit 1

# input: text file with the list of AO2D files to analyse (one path per line)
# default: real data, input_data.txt -> data/AO2D.root
# to run this step on the MC instead:  INPUTFILE=input_MC.txt ./runSkeleton.sh
INPUTFILE="${INPUTFILE:-input_data.txt}"

# Task to run. Default: the checkpoint of this step.
# To run YOUR code (the skeleton you are modifying) with the configuration of this step:
#   TASK=o2-analysistutorial-lf-strangeness-derived-v0-skeleton ./runSkeleton.sh
TASK="${TASK:-o2-analysistutorial-lf-strangeness-derived-v0-skeleton}"

OPTION="-b --configuration json://$DIR_THIS/jsonConfigs/$STEP.json"

# If checkTables.C reports old table versions, prepend the suggested converters, e.g.:
#   CONVERTERS='o2-analysis-lf-straevselsconverter6 ${OPTION} | ' ./runSkeleton.sh
CONVERTERS="${CONVERTERS:-}"
WORKFLOW="${TASK} ${OPTION}"

eval "${CONVERTERS}${WORKFLOW} --aod-file @${INPUTFILE}" > "$LOGFILE" 2>&1

# report status
rc=$?
if [ $rc -eq 0 ]; then
  echo "No problems!"
  mkdir -p "$DIR_THIS/results/$STEP"
  mv "AnalysisResults.root" "$DIR_THIS/results/$STEP/AnalysisResults.root"
  mv "dpl-config.json" "$DIR_THIS/results/$STEP/dpl-config-$STEP.json"
else
  echo "Error: Exit code $rc"
  echo "Check the log file $LOGFILE"
  exit $rc
fi
