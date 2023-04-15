#!/bin/bash
## See https://alice-o2.web.cern.ch/copyright for details of the copyright
## holders. All rights not expressly granted are reserved.
##
## This software is distributed under the terms of the GNU General Public
## License v3 (GPL Version 3), copied verbatim in the file "COPYING".
##
## In applying this license CERN does not waive the privileges and immunities
## granted to it by virtue of its status as an Intergovernmental Organization
## or submit itself to any jurisdiction.

### \author Anton Riedel

# in case the central server is not reachable, set SE yourself
# https://alimonitor.cern.ch/stats?page=SE/table
# export alien_CLOSE_SE=ALICE::UPB::EOS

# make sure you can connect to central server
# alien-token-init

# pass config file as argument or use the default
ConfigFile="${1:-run-config.json}"

# pass the same options to all workflows
Options=("-b" "--configuration" "json://${ConfigFile}")

o2-analysis-timestamp "${Options[@]}" |
        o2-analysis-event-selection "${Options[@]}" |
        o2-analysis-multiplicity-table "${Options[@]}" |
        o2-analysis-track-propagation "${Options[@]}" |
        o2-analysis-pid-tpc-base "${Options[@]}" |
        o2-analysis-pid-tpc "${Options[@]}" |
        o2-analysistutorial-cf-femtodream-tutorial-0 "${Options[@]}"
