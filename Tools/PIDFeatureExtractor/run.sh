#!/bin/bash

# PID Feature Extractor Workflow (simplified, table-based)
# Detectors: TPC, TOF, TRD, ITS, EMCal, HMPID + centrality (FT0C) - the 7
# detector groups / 34-feature contract used by the FSE PID model.
# Event level: Centrality FT0C (Pb-Pb Run 3)
#
# Output is two AOD-joinable tables (PidFeaturesData / PidFeaturesMc)
# written via the framework's own Produces<> mechanism; set the usual
# --aod-writer-json if you want to control where the output
# AnalysisResults-style file goes. CSV export alongside the table is
# available via exportCsv (off by default) - see pidFeatureExtractorConfig.json.
#
# DPG cuts (eta/pT/DCA/TPC-cluster) and Bayesian PID (TPC alone is enough;
# TOF folded in if present; priors configurable) are both optional, set in
# the same config block - see README.md for details.
#
# Mode is a JSON choice: pidFeatureExtractorConfig.json -> "pid-feature-extractor"
#   "processData": 1, "processMc": 0   -> real/raw data (no MC truth)
#   "processData": 0, "processMc": 1   -> MC (reconstructed + truth)
# This task does not itself abort if you leave both/neither on - the
# framework will just run whichever you've enabled or do nothing productive
# if neither is. Double check your config.

CONFIG="$(pwd)/pidFeatureExtractorConfig.json"
OPTION="-b --configuration json://${CONFIG}"

# CRITICAL: Add shared memory flag (~half your available system RAM)
SHM_SIZE="--shm-segment-size 4000000000"

EXTRACTOR=~/alice/sw/BUILD/O2Physics-latest/O2Physics/stage/bin/o2-analysis-pid-feature-extractor
INFERENCE=~/alice/sw/BUILD/O2Physics-latest/O2Physics/stage/bin/o2-analysis-pid-onnx-inference

echo "Starting O2Physics PID Feature Extraction + Inference Workflow..."
echo "Using configuration: ${CONFIG}"
echo "Shared memory segment size: ${SHM_SIZE}"
echo ""

# Pipeline:
#   timestamp → event selection → track propagation → track selection
#   (needed for requireGlobalTrackInFilter())
#   → TPC PID → TOF PID → TOF beta
#   → multiplicity → centrality (needed for CentFT0Cs, both MC and data)
#   → feature extractor → ONNX inference
#
# TRD, ITS, EMCal: already in TracksExtra — no extra task needed
# HMPID:           already in AO2D O2hmpid_001 — no extra task needed
#
# The inference task consumes PidFeaturesData/PidFeaturesMc directly (DPL
# wires the table dependency automatically since both tasks run in one
# workflow here) - drop the last pipe stage if you only want the features
# table and don't want to run inference.

o2-analysis-timestamp ${OPTION} ${SHM_SIZE} | \
o2-analysis-event-selection ${OPTION} ${SHM_SIZE} | \
o2-analysis-track-propagation ${OPTION} ${SHM_SIZE} | \
o2-analysis-trackselection ${OPTION} ${SHM_SIZE} | \
o2-analysis-pid-tpc-base ${OPTION} ${SHM_SIZE} | \
o2-analysis-pid-tpc ${OPTION} ${SHM_SIZE} | \
o2-analysis-pid-tof-base ${OPTION} ${SHM_SIZE} | \
o2-analysis-pid-tof ${OPTION} ${SHM_SIZE} | \
o2-analysis-pid-tof-beta ${OPTION} ${SHM_SIZE} | \
o2-analysis-multiplicity-table ${OPTION} ${SHM_SIZE} | \
o2-analysis-centrality-table ${OPTION} ${SHM_SIZE} | \
${EXTRACTOR} ${OPTION} ${SHM_SIZE} | \
${INFERENCE} ${OPTION} ${SHM_SIZE}

EXIT_CODE=$?

echo ""
if [ $EXIT_CODE -eq 0 ]; then
  echo "✓ Workflow completed successfully!"
else
  echo "✗ Workflow failed with exit code: $EXIT_CODE"
fi

exit $EXIT_CODE
