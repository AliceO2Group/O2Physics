# Heavy-flavour analysis tutorial

## Introduction

Welcome to the heavy-flavour analysis tutorial!

This directory contains the code and configuration for running a minimalistic version of the D<sup>0</sup> analysis on Run 3 real data.

See the [official PWG-HF O<sup>2</sup> documentation](https://aliceo2group.github.io/analysis-framework/docs/advanced-specifics/pwghf.html) for the overview of the full heavy-flavour analysis framework.

## Setup

1. Create a dedicated working directory on your device.
  This is the directory, where you will put your input files, run the code and produce the output files, so it is a good idea to have it outside the O2Physics repository.
2. Put the input file(s) `AO2D.root` (and `AnalysisResults_trees.root`) in the working directory.
3. Load the O2Physics environment.

## Skim production

The file `AnalysisResults_trees.root` contains a derived table with track index skims of 2-prong decay candidates.
This table contains paired track indices which point to the track table in the parent `AO2D.root` file.
It is produced from the `AO2D.root` file by a dedicated workflow `o2-analysistutorial-hf-skim-creator-mini` (implemented in [`skimCreatorMini.cxx`](skimCreatorMini.cxx)).

If you need to produce these derived skims, you can do that by executing the [`run_skim.sh`](run_skim.sh) bash script in the working directory:

```bash
bash ~/alice/O2Physics/Tutorials/PWGHF/run_skim.sh
```

It will use the configuration from [`dpl-config_skim.json`](dpl-config_skim.json).

## Mini task

The mini task workflow `o2-analysistutorial-hf-task-mini` (implemented in [`taskMini.cxx`](taskMini.cxx)) is a simplified version of the D<sup>0</sup> analysis chain part which includes:

- the 2-prong candidate creator,
- the D<sup>0</sup> candidate selector,
- the D<sup>0</sup> analysis task.

The first step (candidate creator) consumes the track index skim table and therefore needs the derived `AnalysisResults_trees.root` file as input.

Processing the derived file requires access to the parent `AO2D.root` file.
The absolute path to the parent file is stored in the derived file but it can be overridden with the parameter `aod-parent-base-path-replacement` in the JSON configuration,
where one has to provide a replacement mask in the format `"old-path-to-parent;new-path-to-parent"`.
(If the parent and the derived files are both in the same directory, `new-path-to-parent` can be empty.)

Run the mini task by executing the [`run_task.sh`](run_task.sh) bash script in the working directory:

```bash
bash ~/alice/O2Physics/Tutorials/PWGHF/run_task.sh
```

It will use the configuration from [`dpl-config_task.json`](dpl-config_task.json) and produce the output file `AnalysisResults.root` with histograms in the working directory.

### Exercise tips

Organise your working environment so that you can easily switch between running the code in the working directory and modifying the code in the O2Physics repository.

When you execute the bash script, the terminal output is saved in the `stdout.log` log file in the working directory.
If an error occurs, the script will report the non-zero exit code and ask you to check the log file to find the problem.

The full O<sup>2</sup> configuration is dumped at the end of processing into the `dpl-config.json` file in the working directory.

See the [rebuilding instructions](https://aliceo2group.github.io/analysis-framework/docs/gettingstarted/installing.html#building-partially-for-development-using-ninja) in the official O<sup>2</sup> analysis framework documentation for advice on compiling your code changes.

See the [Troubleshooting](https://aliceo2group.github.io/analysis-framework/docs/troubleshooting/) section of the official O<sup>2</sup> analysis framework documentation for debugging advice.
