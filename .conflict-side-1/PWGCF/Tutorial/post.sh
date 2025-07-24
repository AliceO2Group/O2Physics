#!/bin/bash

config_file=$1

o2-analysis-cf-femtodream-debug-track -b --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-cf-femtodream-pair-track-track --configuration json://"$config_file" --aod-memory-rate-limit 600000000
