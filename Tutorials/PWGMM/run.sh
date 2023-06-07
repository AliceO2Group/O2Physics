OPTION="-b --configuration json://myConfig.json"

o2-analysistutorial-mm-my-example-task ${OPTION} | o2-analysis-track-propagation ${OPTION} | o2-analysis-timestamp ${OPTION}
