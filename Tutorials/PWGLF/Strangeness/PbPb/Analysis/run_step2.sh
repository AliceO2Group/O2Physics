OPTION="-b --configuration json://configuration_step2.json"

o2-analysis-lf-cascadepid ${OPTION} | o2-analysis-lf-cascadespawner ${OPTION} | o2-analysistutorial-lf-strangeness-pbpb-step2 ${OPTION} --aod-file @input_data.txt