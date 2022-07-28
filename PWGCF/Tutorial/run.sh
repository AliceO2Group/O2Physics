#!/bin/bash

config_file=$1

o2-analysis-timestamp --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-event-selection --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-multiplicity-table --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-track-propagation --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-pid-tpc --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-cf-cf-tutorial-5 --configuration json://"$config_file" --aod-memory-rate-limit 600000000
