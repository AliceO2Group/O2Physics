# PID ML in O2

Particle identification is essential in most of the analyzes. The following tutorial will help you to make use of the machine learning models to improve purity and efficiency of particle kinds for your analysis. A single model is tailored to a specific particle kind, e.g., pions with PID 211. For each track, the model returns a float value in [0, 1] which measures the ''certainty'' of the model that this track is of given kind.

## PidONNXModel

This class represents a single ML model from an ONNX file. It requires the following parameters:
- path to local top directory with PID ML models. If `useCCDB` is set to `true`, this is the location of the downloaded files
- CCDB path to the top directory with PID ML models
- boolean flag: whether CCDB should be used
- CCDB Api instance created in an analysis task
- timestamp of the input analysis data -- neded to choose appropriate model
- PID to be checked
- detector setup: what detectors should be used for identification. It is described by enum PidMLDetector. Currently available setups: TPC, TPC+TOF, TPC+TOF+TRD
- minimum certainty for accepting a track to be of given PID

Let's assume your `PidONNXModel` instance is named `pidModel`. Then, inside your analysis task `process()` function, you can iterate over tracks and call: `pidModel.applyModel(track);` to get the certainty of the model. You can also use `pidModel.applyModelBoolean(track);` to receive a true/false answer, whether the track can be accepted based on the minimum certainty provided to the `PidONNXModel` constructor.

You can check a [simple analysis task example](https://github.com/AliceO2Group/O2Physics/blob/master/Tools/PIDML/simpleApplyPidOnnxModel.cxx). It uses configurable parameters and shows how to calculate the data timestamp. Note that, however, calculation of the timestamp requires subscribing to `aod::Collisions` and `aod::BCsWithTimestamps`. On the other hand, it is possible to use locally stored models, and then the timestamp is not used, so it can be a dummy value. `processTracksOnly` presents how to analyze on local-only PID ML models.

## PidONNXInterface

This is a wrapper around PidONNXModel that contains several models. 
