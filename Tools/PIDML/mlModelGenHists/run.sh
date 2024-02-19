#!/bin/bash

config_file="my-dpl-config.json"

  o2-analysis-tracks-extra-converter --configuration json://$config_file -b |
  #o2-analysis-track-propagation --configuration json://$config_file -b |
o2-analysis-timestamp --configuration json://$config_file -b |
  o2-analysis-trackextension --configuration json://$config_file -b |
  o2-analysis-trackselection --configuration json://$config_file -b |
  o2-analysis-multiplicity-table --configuration json://$config_file -b |
  #o2-analysis-fdd-converter --configuration json://$config_file -b |
  o2-analysis-bc-converter --configuration json://$config_file -b |
  o2-analysis-collision-converter --configuration json://$config_file -b |
  o2-analysis-zdc-converter --configuration json://$config_file -b |
  o2-analysis-pid-tof-base --configuration json://$config_file -b |
  o2-analysis-pid-tof-beta --configuration json://$config_file -b |
  o2-analysis-pid-tof-full --configuration json://$config_file -b |
  o2-analysis-pid-tpc-full --configuration json://$config_file -b |
  o2-analysis-pid-tpc-base --configuration json://$config_file -b |
  o2-analysis-ml-model-gen-hists --configuration json://$config_file -b
