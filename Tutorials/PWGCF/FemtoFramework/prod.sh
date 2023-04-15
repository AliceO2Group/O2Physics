#!/bin/bash

# in case the central server is not reachable, set SE yourself
# https://alimonitor.cern.ch/stats?page=SE/table
# export alien_CLOSE_SE=ALICE::UPB::EOS

# make sure you can connect to central server
# alien-token-init

# pass config file as argument or use the default
ConfigFile="${1:-prod-config.json}"

# pass the same options to all workflows
Options=("-b" "--configuration" "json://${ConfigFile}")

o2-analysis-timestamp "${Options[@]}" |
        o2-analysis-event-selection "${Options[@]}" |
        o2-analysis-multiplicity-table "${Options[@]}" |
        o2-analysis-track-propagation "${Options[@]}" |
        o2-analysis-pid-tpc-base "${Options[@]}" |
        o2-analysis-pid-tpc "${Options[@]}" |
        o2-analysis-pid-tof-base "${Options[@]}" |
        o2-analysis-pid-tof "${Options[@]}" |
        o2-analysis-lf-lambdakzerobuilder "${Options[@]}" |
        o2-analysis-cf-femtodream-producer "${Options[@]}" --aod-writer-resfile FemtoAO2D --aod-writer-keep AOD/FEMTODREAMPARTS/0,AOD/FEMTODREAMCOLS/0,AOD/FEMTODEBUGPARTS/0

exit 0
