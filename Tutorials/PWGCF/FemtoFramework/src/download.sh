#!/bin/bash

# available periods and runs can be found here: https://alimonitor.cern.ch/hyperloop/datasets
# For development, try to find files with reasonable size (<1GB) for faster execution
File="/alice/data/2022/LHC22m/523559/apass3/1800/o2_ctf_run00523559_orbit0438833895_tf0000935394_epn151/001/AO2D.root"

alien_cp "alien://${File}" "file://AO2D.root"

[ ! -f Input.txt ] && echo "AO2D.root" >>Input.txt

exit 0
