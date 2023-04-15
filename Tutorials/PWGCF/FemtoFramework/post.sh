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

# pass config file as argument or use the default
ConfigFile="${1:-post-config.json}"

# pass the same options to all workflows
Options="-b --configuration json://${ConfigFile}"

o2-analysis-cf-femtodream-debug-track ${Options} |
        o2-analysis-cf-femtodream-pair-track-track ${Options}

exit 0
