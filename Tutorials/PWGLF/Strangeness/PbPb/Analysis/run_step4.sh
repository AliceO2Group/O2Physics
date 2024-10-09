OPTION="-b --configuration json://configuration_step4.json"

o2-analysis-lf-cascadepid ${OPTION} | o2-analysis-lf-cascadespawner ${OPTION} | o2-analysistutorial-lf-strangeness-pbpb-step4 ${OPTION} --aod-file @input_data.txt