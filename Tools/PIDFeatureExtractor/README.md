# PID Feature Extractor + ONNX Inference

This provides particle identification for ALICE Run 3 Pb-Pb analyses using
a trained ML model (a detector-aware attention model conditioned on which
detectors each track actually has hits in - TPC, TOF, TRD, ITS, EMCal,
HMPID, plus event centrality). Two tasks:

- **`pidFeatureExtractor.cxx`** reads AO2D data and writes out the model's
  input features - kinematics, per-detector PID signals, and detector
  presence flags - to a ROOT file (and optionally CSV).
- **`pidOnnxInference.cxx`** takes that file, runs the trained ONNX model
  over it, and writes back a probability for each particle species
  (pion / kaon / proton / electron) per track.

You run the extractor first, then inference on its output. Both are
regular AOD-subscribing tasks; `PidOnnxInference` does its real work by
reading the extractor's output file directly rather than the AOD data it's
subscribed to (see "Running" below for what that means in practice).

## PidFeatureExtractor

An ordinary AOD-subscribing analysis task. It reads track and collision
data and, for each track passing the (optional, off by default) quality
cuts, writes one row containing:

- kinematics (momentum, eta, phi, DCA)
- per-detector signals for TPC, TOF, TRD, ITS, EMCal, and HMPID, each with
  a flag saying whether that detector actually has a hit on this track
- event centrality
- a Bayesian PID posterior, for comparison against the ML model
- for MC only: the true particle ID and whether it's a physical primary

Mode is a runtime switch - enable `processData` for real data or
`processMc` for MC (reconstructed + truth), not both.

### PidFeatureExtractor options

| Option                 | Default          | What it does                                                     |
|------------------------|------------------|------------------------------------------------------------------|
| `outputPath`           | `pid_features`   | Output file base name                                            |
| `exportROOT`           | `true`           | Write a ROOT file                                                |
| `exportCsv`            | `false`          | Also write CSV                                                   |
| `etaMin` / `etaMax`    | `-99` / `99`     | Eta cut - wide open by default (no cut)                          |
| `ptMin` / `ptMax`      | `0` / `9999`     | pT cut, GeV/c - wide open by default                             |
| `dcaXYMax` / `dcaZMax` | `9999` / `9999`  | DCA cuts, cm - wide open by default                              |
| `itsMinClusters`       | `0`              | Minimum ITS clusters - `0` = no cut                              |
| `tpcMinClusters`       | `0`              | Minimum TPC clusters - `0` = no cut                              |
| `computeBayesianPid`   | `true`           | Compute the comparison Bayesian posterior                        |
| `bayesianPriors`       | flat (`1,1,1,1`) | Per-species priors `[pi, ka, pr, el]` for the Bayesian posterior |

All the cuts default to "off" - tighten them in your config if you want
quality selection applied here rather than downstream.

## PidOnnxInference

Takes the file `PidFeatureExtractor` wrote and runs the trained ONNX model
over it, row by row, in `init()` - not per-collision. The model can be
loaded either from CCDB or from a local file, which is handled by
`o2::analysis::MlResponse` (`Tools/ML/MlResponse.h`).

This is still a normal AOD-subscribing task, so it needs a valid AO2D
file to run at all, the same as any other task in this repository - but
it doesn't actually use that data; `process()` is intentionally empty.
Point it at any valid AO2D (the same one you ran the extractor against is
the obvious choice) purely to satisfy the pipeline.

By default it assumes every detector group is present and usable, exactly
as the input data says. If you want to see how the model behaves with a
detector deliberately left out - for testing, or to match a specific
detector configuration - each group can be switched off independently;
turning one off overrides the data for that group, the same way a genuine
detector miss would look.

### PidOnnxInference options

| Option              | Default                     | What it does                                                                 |
|---------------------|-----------------------------|------------------------------------------------------------------------------|
| `inputRootFile`     | `pid_features_data.root`    | File written by `PidFeatureExtractor`                                        |
| `inputTreeName`     | `pid_features`              | Tree name inside it                                                          |
| `outputPath`        | `pid_predictions`           | Output file base name                                                        |
| `exportCsv`         | `false`                     | Also write CSV                                                               |
| `loadModelFromCcdb` | `true`                      | Load the model from CCDB; set `false` to use a local file instead            |
| `ccdbUrl`           | `http://alice-ccdb.cern.ch` |                                                                              |
| `modelPathsCcdb`    | *(placeholder)*             | CCDB path to your model - set this to a real path before running             |
| `timestampCcdb`     | `-1`                        | `-1` = latest                                                                |
| `onnxFileNames`     | `pid_feature_model.onnx`    | Local model file, used when `loadModelFromCcdb` is `false`                   |
| `useTPC`            | `true`                      | Include TPC. Set `false` to exclude it from inference regardless of the data |
| `useTOF`            | `true`                      | Include TOF                                                                  |
| `useTRD`            | `true`                      | Include TRD                                                                  |
| `useITS`            | `true`                      | Include ITS                                                                  |
| `useEMCal`          | `true`                      | Include EMCal                                                                |
| `useHMPID`          | `true`                      | Include HMPID                                                                |
| `useCentrality`     | `true`                      | Include event centrality                                                     |

Output columns are `mlProbPi`, `mlProbKa`, `mlProbPr`, `mlProbEl` (one
probability per species) and `mlPredictedClass` (the most likely species,
as an index: `0`=pion, `1`=kaon, `2`=proton, `3`=electron).

## Running

Both use the usual `--configuration json://your-config.json` mechanism,
and both are AOD-subscribing tasks - `PidOnnxInference` just doesn't use
the AOD data it's given, it reads `PidFeatureExtractor`'s output file
instead. Run the extractor first:

```bash
#!/bin/bash

config_file="my-config.json"

o2-analysis-timestamp --configuration json://$config_file -b |
    o2-analysis-event-selection --configuration json://$config_file -b |
    o2-analysis-track-propagation --configuration json://$config_file -b |
    o2-analysis-trackselection --configuration json://$config_file -b |
    o2-analysis-pid-tpc-base --configuration json://$config_file -b |
    o2-analysis-pid-tpc --configuration json://$config_file -b |
    o2-analysis-pid-tof-base --configuration json://$config_file -b |
    o2-analysis-pid-tof --configuration json://$config_file -b |
    o2-analysis-pid-tof-beta --configuration json://$config_file -b |
    o2-analysis-multiplicity-table --configuration json://$config_file -b |
    o2-analysis-centrality-table --configuration json://$config_file -b |
    o2-analysis-pid-feature-extractor --configuration json://$config_file -b
```

Then run inference, once the extractor has finished and its output file
exists. Any valid AO2D works as input here, since its content is unused -
reusing the same one is the simplest choice:

```bash
#!/bin/bash

config_file="my-config.json"

o2-analysis-timestamp --configuration json://$config_file -b |
    o2-analysis-pid-onnx-inference --configuration json://$config_file -b
```
