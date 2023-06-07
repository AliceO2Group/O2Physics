export OPTIONS="-b --configuration json://config.json --shm-segment-size 10000000000"
o2-analysis-timestamp ${OPTIONS} \
| o2-analysis-event-selection ${OPTIONS} \
| o2-analysis-multiplicity-table ${OPTIONS} \
| o2-analysis-trackselection ${OPTIONS} \
| o2-analysis-track-propagation ${OPTIONS} \
| o2-analysis-je-emcal-correction-task --clusterDefinition "kV3Default" ${OPTIONS} \
| o2-analysistutorial-em-emc-cluster-step0 ${OPTIONS} \
| o2-analysis-zdc-converter
