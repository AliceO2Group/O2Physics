#!/bin/bash

# pass config file as argument or use the default
ConfigFile="${1:-post-config.json}"

# pass the same options to all workflows
Options=("-b" "--configuration" "json://${ConfigFile}")

o2-analysis-cf-femtodream-debug-track "${Options[@]}" |
        o2-analysis-cf-femtodream-pair-track-track "${Options[@]}"

exit 0
