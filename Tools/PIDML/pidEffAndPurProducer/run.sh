#!/bin/bash

config_file="my-dpl-config.json"
FILENAME=AO2D_18G4_run285064_023

o2-analysis-tracks-extra-converter --configuration json://$config_file -b |
o2-analysis-timestamp --configuration json://$config_file -b |
o2-analysis-trackextension --configuration json://$config_file -b |
o2-analysis-trackselection --configuration json://$config_file -b |
o2-analysis-multiplicity-table --configuration json://$config_file -b |
o2-analysis-bc-converter --configuration json://$config_file -b |
o2-analysis-collision-converter --configuration json://$config_file -b |
o2-analysis-zdc-converter --configuration json://$config_file -b |
o2-analysis-pid-tof-base --configuration json://$config_file -b |
o2-analysis-pid-tof-beta --configuration json://$config_file -b |
o2-analysis-pid-tof-full --configuration json://$config_file -b |
o2-analysis-pid-tpc-full --configuration json://$config_file -b |
o2-analysis-pid-tpc-base --configuration json://$config_file -b |
o2-analysis-pid-eff-and-pur-producer --aod-file $AO2DS_DIR/$FILENAME.root --configuration json://$config_file -b \
  --aod-writer-keep AOD/PIDEFFANDPURRES/0:::PID_RESULTS_$FILENAME
mv AnalysisResults.root MC_HISTS_$FILENAME.root
