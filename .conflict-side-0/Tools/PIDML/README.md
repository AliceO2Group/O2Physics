# PID ML in O2

Particle identification is essential in most of the analyzes.
The PID ML interface will help you to make use of the machine learning models to improve purity and efficiency of particle kinds for your analysis.
A single model is tailored to a specific particle kind, e.g., pions with PID 211. For each track, the model returns a float value in [0, 1] which measures the ''certainty'' of the model that this track is of given kind.

## PidONNXModel

This class represents a single ML model from an ONNX file. It requires the following parameters:
- path to local top directory with PID ML models. If `useCCDB` is set to `true`, this is the location of the downloaded files
- CCDB path to the top directory with PID ML models
- boolean flag: whether CCDB should be used
- CCDB Api instance created in an analysis task
- timestamp of the input analysis data -- neded to choose appropriate model
- PID to be checked
- minimum certainty for accepting a track to be of given PID
- *p* limits array - specifiying p limits for each detector configuration (TPC, TPC+TOF, TPC+TOF+TRD)

Let's assume your `PidONNXModel` instance is named `pidModel`.
Then, inside your analysis task `process()` function, you can iterate over tracks and call: `pidModel.applyModel(track);` to get the certainty of the model.
You can also use `pidModel.applyModelBoolean(track);` to receive a true/false answer, whether the track can be accepted based on the minimum certainty provided to the `PidONNXModel` constructor.

You can check [a simple analysis task example](https://github.com/AliceO2Group/O2Physics/blob/master/Tools/PIDML/simpleApplyPidOnnxModel.cxx).
It uses configurable parameters and shows how to calculate the data timestamp. Note that the calculation of the timestamp requires subscribing to `aod::Collisions` and `aod::BCsWithTimestamps`.
For Hyperloop tests, you can set `cfgUseFixedTimestamp` to true with `cfgTimestamp` set to the default value.

On the other hand, it is possible to use locally stored models, and then the timestamp is not used, so it can be a dummy value. `processTracksOnly` presents how to analyze on local-only PID ML models.

## PidONNXInterface

This is a wrapper around PidONNXModel that contains several models. It has the possibility of automatically selecting best fitting model for the data. You can also configure manually model selection. Some obligatory class parameters are the same as for `PidONNXModel`:
- path to local top directory with PID ML models. If `useCCDB` is set to `true`, this is the location of the downloaded files
- CCDB path to the top directory with PID ML models
- boolean flag: whether CCDB should be used
- CCDB Api instance created in an analysis task
- timestamp of the input analysis data -- neded to choose appropriate model

Then, obligatory parameters for the interface:
- a vector of int output PIDs
- a 2-dimensional LabeledArray of *p* limits for each PID, for each detector configuration. It describes the minimum *p* values at which each next detector should be included for predicting given PID
- a vector of minimum certainties for each PID for accepting a track to be of this PID
- boolean flag: whether to switch on auto mode. If true, then *p*T limits and minimum certainties can be passed as an empty array and an empty vector, and the interface will fill them with default configuration:
  - *p* limits: same values for all PIDs: 0.0 (TPC), 0.5 (TPC + TOF), 0.8 (TPC + TOF + TRD)
  - minimum certainties: 0.5 for all PIDs

You can use the interface in the same way as the model, by calling `applyModel(track)` or `applyModelBoolean(track)`. The interface will then call the respective method of the model selected with the aforementioned interface parameters.

In the future, the interface will be extended with a more sophisticated model selection strategy. Moreover, it will also allow for using a backup model in the case the best fit model doesn't exist.

There is again [a simple analysis task example](https://github.com/AliceO2Group/O2Physics/blob/master/Tools/PIDML/simpleApplyPidOnnxInterface.cxx) for using `PidONNXInterface`. It is analogous to the `PidONNXModel` example.

## Notes on running the current task

Currently, only models for run 285064 (timestamp interval: 1524176895000 - 1524212953000) are uploaded to CCDB, so you can use hardcoded timestamp 1524176895000 for tests.

Both model and interface analysis examples can be run with a script:

### Script for Run2 Converted to Run3 data
```bash
#!/bin/bash

config_file="my-config.json"

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
    o2-analysis-simple-apply-pid-onnx-model --configuration json://$config_file -b
```
Remember to set every setting, which states that helper task should process Run2 data to `true`.

### Script for Run3 data
```bash
#!/bin/bash

config_file="my-config.json"

o2-analysis-timestamp --configuration json://$config_file -b |
    o2-analysis-event-selection --configuration json://$config_file -b |
    o2-analysis-trackselection --configuration json://$config_file -b |
    o2-analysis-multiplicity-table --configuration json://$config_file -b |
    o2-analysis-track-propagation --configuration json://$config_file -b |
    o2-analysis-pid-tof-base --configuration json://$config_file -b |
    o2-analysis-pid-tof-beta --configuration json://$config_file -b |
    o2-analysis-pid-tof-full --configuration json://$config_file -b |
    o2-analysis-pid-tpc-full --configuration json://$config_file -b |
    o2-analysis-pid-tpc-base --configuration json://$config_file -b |
    o2-analysis-simple-apply-pid-onnx-model --configuration json://$config_file -b
```
Remember to set every setting, which states that helper task should process Run3 data to `true`.


Replace "model" with "interface" in the last line if you want to run the interface workflow.
