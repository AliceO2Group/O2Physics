# Strangeness hands-on — V0s in pp @ √s = 5.36 TeV (derived data)

In this hands-on session you will write an analysis task for **K0S, Λ and Λbar** on strangeness derived
data of pp collisions at 5.36 TeV, in 5 steps, starting from a skeleton.

You will learn how to:
- use the strangeness derived data and structure your analysis task;
- reconstruct K0S, Λ and Λbar;
- apply topological, kinematic and PID (TPC and TOF) selections;
- access the MC information of the V0s to compute the efficiency.

## Preparation

Update O2Physics and enter the environment:

```bash
alienv enter O2Physics/latest
which o2-analysistutorial-lf-strangeness-derived-v0-skeleton    # must exist
which o2-analysistutorial-lf-strangeness-derived-v0-step0       # checkpoint
```

Go to the folder of the hands-on in your O2Physics source directory:

```bash
cd ~/alice/O2Physics/Tutorials/PWGLF/Strangeness/Derived/AnalysisV0
```

Download the two input files (strangeness derived data, ~<200> MB each):

```bash
mkdir -p data MC
# real data -> data/AO2D.root   (used by skeleton and steps 0-3)
rsync -r -u --progress USERNAME@lxplus.cern.ch:/eos/user/<u>/<user>/tutorialData/data/AO2D.root data
# Monte Carlo -> MC/AO2D.root   (used by step 4, which needs the MC truth)
rsync -r -u --progress USERNAME@lxplus.cern.ch:/eos/user/<u>/<user>/tutorialData/MC/AO2D.root MC
```

Alternatively, download them from CERNBox links.

The run scripts do not contain the file paths: they read them from two text files (one path per line,
relative to this folder or absolute), which you have to create (the paths are also given in the slides):

```bash
echo "data/AO2D.root" > input_data.txt   # used by runSkeleton.sh, runStep0-3.sh
echo "MC/AO2D.root" > input_MC.txt       # used by runStep4.sh
```

If you save the files elsewhere, edit these two text files. To run steps 0-3 on the MC (e.g. to compare
with data): `INPUTFILE=input_MC.txt ./runStep2.sh`.

## How to work

The code you will modify is the skeleton in this folder:
`Tutorials/PWGLF/Strangeness/Derived/AnalysisV0/strangeness_derived_v0_skeleton.cxx`.
The checkpoints of each step are in the same folder (`strangeness_derived_v0_stepX.cxx`).

```bash
./runSkeleton.sh                                   # "No problems!" must be printed
root -l results/skeleton/AnalysisResults.root
```

After modifying and recompiling the skeleton, run it with the configuration of a given step:

```bash
TASK=o2-analysistutorial-lf-strangeness-derived-v0-skeleton ./runStep1.sh
```

Without `TASK=...` the checkpoint of the step is run: useful to compare with your results.

| Step     | What is added                                                                           | Tables                                          | What to look at                                 |
|----------|-----------------------------------------------------------------------------------------|-------------------------------------------------|-------------------------------------------------|
| skeleton | loop over events, sel8, \|z\| < 10 cm                                                   | `StraCollisions`, `StraEvSels`                  | `eventSelection/hVertexZRec`                    |
| 0        | loop over V0s, K0S, Λ, Λbar invariant masses                                            | `V0Cores`                                       | `kzeroShort/hMassK0Short`, `lambda/hMassLambda` |
| 1        | topological selections (`Filter` + dynamic columns), rapidity, competing mass rejection | `V0Cores`                                       | `hMass*Selected`, `hV0CosPA`, `hV0Radius`       |
| 2        | TPC PID of the daughters (`posTrackExtra_as<>`)                                         | `V0Extras`, `DauTrackExtras`, `DauTrackTPCPIDs` | `hNSigmaPos/NegPionFromK0s`, `hMass*Selected`   |
| 3        | TOF PID optimised for secondary particles (`strangenesstofpid` task)                    | `V0TOFNSigmas`                                  | `hMass*SelectedWithTOF`                         |
| 4        | MC truth: true reconstructed and generated V0s                                          | `V0CoreMCLabels`, `V0MCCores`, `V0MCCollRefs`   | `hPt*TrueRec`, `genParticles/hPt*Gen`           |

Step 1 selections (in `jsonConfigs/step1.json`):

| Selection                    | Value       | Getter                         | Column                       |
|------------------------------|-------------|--------------------------------|------------------------------|
| DCA pos./neg. daughter to PV | > 0.06 cm   | `dcapostopv()`, `dcanegtopv()` | standard                     |
| DCA between daughters        | < 0.5 cm    | `dcaV0daughters()`             | standard                     |
| cosine of pointing angle     | > 0.995     | `v0cosPA()`                    | **standard** in derived data |
| decay radius                 | > 0.9 cm    | `v0radius()`                   | dynamic                      |
| rapidity                     | \|y\| < 0.5 | `yK0Short()`, `yLambda()`      | dynamic                      |

Hints for the steps:
- **Step 0**: look at the `aod::V0Cores` table in `PWGLF/DataModel/LFStrangenessTables.h` and at the getters
  `.mK0Short()`, `.mLambda()`, `.mAntiLambda()`.
- **Step 1**: keep only standard V0s (`v0Type() == 1`); a `Filter` cannot cut on dynamic columns.
- **Step 2**: K0S → π⁺π⁻, Λ → pπ⁻, Λbar → p̄π⁺. You need `LFStrangenessPIDTables.h`.
- **Step 3**: look at `.positiveHasTOF()`, `.negativeHasTOF()` and `.tofNSigmaK0PiPlus()`, `.tofNSigmaLaPr()`, … in `LFStrangenessPIDTables.h`.
- **Step 4**: `has_v0MCCore()`, `v0MCCore_as<aod::V0MCCores>()`, `.pdgCode()` with `PDG_t::kK0Short`,
  `PDG_t::kLambda0`, `PDG_t::kLambda0Bar`; select only primary particles (`.isPhysicalPrimary()`).

## Checkpoints

- **Step 0**: peaks at ~0.497 GeV/c² (K0S) and ~1.116 GeV/c² (Λ, Λbar), already visible thanks to the derived data pre-selections.
- **Step 1**: the background is reduced.
- **Step 2**: the Nσ bands of the pions from K0S are centred at 0; the background under the Λ and Λbar peaks is further reduced.
- **Step 3**: TOF further suppresses the combinatorial background.
- **Step 4**: `hPt*TrueRec / hPt*Gen` gives the efficiency (× acceptance) as a function of pT.

## If something does not run

```bash
root -l -b -q "checkTables.C(\"$(cat input_data.txt)\")"   # data: step 4 "CANNOT run" is expected
root -l -b -q "checkTables.C(\"$(cat input_MC.txt)\")"    # MC: all steps must be ok
```

The macro tells you which steps can run with the file and, if a table has an old version, which converter
to use:

```bash
CONVERTERS='o2-analysis-lf-straevselsconverter6 ${OPTION} | ' ./runStep0.sh
```
