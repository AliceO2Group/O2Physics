# Heavy-flavour analysis tutorial

## Introduction

Welcome to the heavy-flavour analysis tutorial!

This directory contains the code and configuration for running a minimalistic version of the D<sup>0</sup> analysis on Run 3 data.

See the [PWG-HF O<sup>2</sup> documentation](https://aliceo2group.github.io/analysis-framework/docs/advanced-specifics/pwghf.html) for the overview of the full heavy-flavour analysis framework.

## Setup

1. Create and enter a dedicated working directory on your device.
  This is the directory, where you will put your input files, run the code, and produce the output files, so it is a good idea to have it outside the O2Physics repository.
1. Put the input file(s) `AO2D.root` (and `AnalysisResults_trees.root`) in the working directory (or in a dedicated directory for input data).
1. Copy the [`input_skim.txt`](input_skim.txt) and [`input_task.txt`](input_task.txt) files into the working directory and adjust the paths inside if needed.
1. Load the O2Physics environment.

## Skim production

The mini skim creator workflow `o2-analysistutorial-hf-skim-creator-mini` (implemented in [`skimCreatorMini.cxx`](skimCreatorMini.cxx)) is a simplified version of the `o2-analysis-hf-track-index-skim-creator` workflow (implemented in [`trackIndexSkimCreator.cxx`](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/TableProducer/trackIndexSkimCreator.cxx)) which is used for the central production of linked derived data for all HF analyses and which performs:

- the HF event selection,
- the HF secondary-track selection,
- the HF secondary-vertex reconstruction and loose selection of found HF decay candidates.

The mini skim creator processes the `AO2D.root` file(s) and produces a derived table [`HfT2Prongs`](DataModelMini.h) with track index skims of 2-prong decay candidates.
This table contains paired track indices which point to tracks in the track table in the parent `AO2D.root` file.

These derived skims are produced by executing the [`run_skim.sh`](run_skim.sh) bash script in the working directory:

```bash
~/alice/O2Physics/Tutorials/PWGHF/run_skim.sh
```

It processes files specified in `./input_skim.txt`, uses the configuration from [`dpl-config_skim.json`](dpl-config_skim.json), and produces the `./AnalysisResults_trees.root` file with the derived table and the `./AnalysisResults.root` file with control histograms.

## Mini task

The mini task workflow `o2-analysistutorial-hf-task-mini` (implemented in [`taskMini.cxx`](taskMini.cxx)) is a simplified version of the D<sup>0</sup> analysis chain part which includes:

- the 2-prong candidate creator (for the full candidate reconstruction),
- the D<sup>0</sup> candidate selector (for the candidate selection),
- the D<sup>0</sup> analysis task (for the analysis of selected candidates and filling of output histograms).

The first step (candidate creator) consumes the track index skim table and therefore needs the derived `AnalysisResults_trees.root` file as input.

Processing the derived file requires access to the parent `AO2D.root` file.
The absolute path to the parent file is stored in the derived file but it can be overridden with the command line parameter `--aod-parent-base-path-replacement "old-path-to-parent;new-path-to-parent"`.
(If the parent and the derived files are both in the same directory, `new-path-to-parent` can be empty.)

Run the mini task by executing the [`run_task.sh`](run_task.sh) bash script in the working directory:

```bash
~/alice/O2Physics/Tutorials/PWGHF/run_task.sh
```

It processes files specified in `./input_task.txt`, uses the configuration from [`dpl-config_task.json`](dpl-config_task.json), and produces the `./AnalysisResults.root` file with histograms.

## Exercise tips

Organise your working environment so that you can easily switch between running the code in the working directory and modifying the code in the O2Physics repository.

When you execute the bash script, the terminal output is saved in the `./stdout.log` log file.
If an error occurs, the script will report the non-zero exit code and ask you to check the log file to find the problem.
If you get errors that have to be tolerated, you need to remove the `--min-failure-level error` option from the command.

The full O<sup>2</sup> configuration is dumped at the end of processing into the `./dpl-config.json` file.
It is useful to compare it with the input configuration file to spot mismatches.

See the [rebuilding instructions](https://aliceo2group.github.io/analysis-framework/docs/gettingstarted/installing.html#building-partially-for-development-using-ninja) in the analysis framework documentation for advice on compiling your code changes.

See the [Troubleshooting](https://aliceo2group.github.io/analysis-framework/docs/troubleshooting/) section of the analysis framework documentation for debugging advice.

Consider using the [Shell rc utilities](https://aliceo2group.github.io/analysis-framework/docs/tools/#shell-rc-utilities) for easier recompilation and debugging.
