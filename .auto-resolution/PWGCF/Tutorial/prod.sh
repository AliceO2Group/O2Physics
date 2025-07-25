#!/bin/bash

config_file=$1

o2-analysis-trackselection --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-timestamp --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-multiplicity-table --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-event-selection --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-track-propagation --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-pid-tpc --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-pid-tof-base --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-pid-tof --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-lf-lambdakzerobuilder --configuration json://"$config_file" --aod-memory-rate-limit 600000000 |
  o2-analysis-cf-femtodream-producer --configuration json://"$config_file" --aod-writer-resfile FemtoAO2D --aod-writer-keep AOD/FEMTODREAMPARTS/0,AOD/FEMTODREAMCOLS/0,AOD/FEMTODEBUGPARTS/0 --aod-memory-rate-limit 600000000
