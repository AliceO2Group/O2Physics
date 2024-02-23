// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_

// ...
#include <vector>
#include <string>

// Task configuration:
Configurable<string> cfTaskName{
  "cfTaskName", "Default task name",
  "set task name - use eventually to determine weights for this task"};
Configurable<bool> cfDryRun{
  "cfDryRun", false,
  "book all histos and run without storing and calculating anything"};
Configurable<bool> cfVerbose{
  "cfVerbose", false,
  "run or not in verbose mode (but not for function calls per particle)"};
Configurable<bool> cfVerboseForEachParticle{
  "cfVerboseForEachParticle", false,
  "run or not in verbose mode (also for function calls per particle)"};
Configurable<bool> cfDoAdditionalInsanityChecks{
  "cfDoAdditionalInsanityChecks", false,
  "do additional insanity checks at run time (this leads to small loss of performance)"};
Configurable<bool> cfUseCCDB{
  "cfUseCCDB", true,
  "access personal files from CCDB or from home dir in AliEn"};
Configurable<string> cfWhatToProcess{
  "cfWhatToProcess", "Rec",
  "Rec = process only reconstructed, Sim = process only simulated, RecSim = process both reconstructed and simulated. For converted Run 1 and Run 2 data, append still either _Run2 or _Run1"};
Configurable<unsigned int> cfRandomSeed{
  "cfRandomSeed", 0,
  "0 = random seed is guaranteed to be unique in space and time"};
Configurable<bool> cfUseFisherYates{
  "cfUseFisherYates", false,
  "use or not Fisher-Yates algorithm to randomize particle indices"};
Configurable<int> cfFixedNumberOfRandomlySelectedTracks{
  "cfFixedNumberOfRandomlySelectedTracks", -1,
  "Set to some integer > 0, to apply and use. Set to <=0, to ignore."};
Configurable<bool> cfUseStopwatch{
  "cfUseStopwatch", false,
  "If true, some basic info on time execution is printed, here and there."};

// Default booking:
Configurable<vector<int>> cfBookEventHistograms{"cfBookEventHistograms", {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) event histogram, ordering is the same as in enum eEventHistograms"};
Configurable<vector<int>> cfBookParticleHistograms{"cfBookParticleHistograms", {1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) particle histogram, ordering is the same as in enum eParticleHistograms"};

// Correlations:
Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false,
                                           "calculate or not correlations"};

// Test0:
Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false,
                                    "calculate or not Test0"};
Configurable<bool> cfCalculateTest0AsFunctionOfIntegrated{"cfCalculateTest0AsFunctionOfIntegrated", true,
                                                          "calculate or not Test0 as a function of integrated"};
Configurable<bool> cfCalculateTest0AsFunctionOfMultiplicity{"cfCalculateTest0AsFunctionOfMultiplicity", true,
                                                            "calculate or not Test0 as a function of multiplicity"};
Configurable<bool> cfCalculateTest0AsFunctionOfCentrality{"cfCalculateTest0AsFunctionOfCentrality", true,
                                                          "calculate or not Test0 as a function of centrality"};
Configurable<bool> cfCalculateTest0AsFunctionOfPt{"cfCalculateTest0AsFunctionOfPt", false,
                                                  "calculate or not Test0 as a function of pt"};
Configurable<bool> cfCalculateTest0AsFunctionOfEta{"cfCalculateTest0AsFunctionOfEta", false,
                                                   "calculate or not Test0 as a function of eta"};
Configurable<string> cfFileWithLabels{"cfFileWithLabels",
                                      "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Particle weights:
Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false,
                                   "use or not phi weights"};
Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false,
                                  "use or not pt weights"};
Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false,
                                   "use or not eta weights"};
Configurable<bool> cfUseDiffPhiPtWeights{"cfUseDiffPhiPtWeights", false,
                                         "use or not differential phi(pt) weights"};
Configurable<bool> cfUseDiffPhiEtaWeights{"cfUseDiffPhiEtaWeights", false,
                                          "use or not differential phi(eta) weights"};
Configurable<string> cfFileWithWeights{"cfFileWithWeights",
                                       "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Nested loops:
Configurable<bool> cfCalculateNestedLoops{"cfCalculateNestedLoops", false,
                                          "cross-check for all events all correlations with nested loops"};
Configurable<bool> cfCalculateCustomNestedLoops{"cfCalculateCustomNestedLoops", false,
                                                "cross-check e-b-e all correlations with custom nested loops"};

// Event cuts:
Configurable<int> cNumberOfEvents_min{
  "cNumberOfEvents_min", -1,
  "minimum number of events to process (set to -1 to ignore)"};
Configurable<int> cNumberOfEvents_max{"cNumberOfEvents_max", 1000000000,
                                      "maximum number of events to process (use only in the local analysis, because after that BailOut() is forced)"};
Configurable<int> cTotalMultiplicity_min{
  "cTotalMultiplicity_min", -1,
  "minimum total multiplicity of an event to be processed (set to -1 to "
  "ignore)"};
Configurable<int> cTotalMultiplicity_max{
  "cTotalMultiplicity_max", 1000000000,
  "maximum total multiplicity of an event to be processed"};
Configurable<int> cSelectedTracks_min{
  "cSelectedTracks_min", -1,
  "minimum number of selected tracks (set to -1 to ignore)"};
Configurable<int> cSelectedTracks_max{"cSelectedTracks_max", 1000000000,
                                      "maximum number of selected tracks"};
Configurable<float> cCentrality_min{"cCentrality_min", 0.,
                                    "minimum centrality"};
Configurable<float> cCentrality_max{"cCentrality_max", 100.,
                                    "maximum centrality"};
Configurable<float> cVertex_x_min{"cVertex_x_min", -10.0,
                                  "minimum vertex x range [cm]"};
Configurable<float> cVertex_x_max{"cVertex_x_max", 10.0,
                                  "maximum vertex x range [cm]"};
Configurable<float> cVertex_y_min{"cVertex_y_min", -10.0,
                                  "minimum vertex y range [cm]"};
Configurable<float> cVertex_y_max{"cVertex_y_max", 10.0,
                                  "maximum vertex y range [cm]"};
Configurable<float> cVertex_z_min{"cVertex_z_min", -10.0,
                                  "minimum vertex z range [cm]"};
Configurable<float> cVertex_z_max{"cVertex_z_max", 10.0,
                                  "maximum vertex z range [cm]"};
Configurable<int> cNContributors_min{
  "cNContributors_min", -1,
  "minimum number of vertex contributors (set to -1 to ignore)"};
Configurable<int> cNContributors_max{"cNContributors_max", 1000000000,
                                     "maximum number of vertex contributors"};
Configurable<int> cImpactParameter_min{
  "cImpactParameter_min", -1,
  "minimum value of impact parameter (can be used only for sim) (set to -1 to ignore)"};
Configurable<int> cImpactParameter_max{"cImpactParameter_max", 1000000000,
                                       "maximum value of impact parameter (can be used osnly for sim)"};

// Particle cuts:
Configurable<float> pt_min{"pt_min", 0.2, "minimum track pt value [GeV/c]"};
Configurable<float> pt_max{"pt_max", 5.0, "maximum track pt value [GeV/c]"};

// Results histograms:
Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false,
                                           "save or not results histograms"};

// Fixed-length binning (default): TBI 20240113 I do it via string + tokenize + Atof(), until I figure out how to pass array directly via configurable
Configurable<string> cFixedLength_pt_bins{"cFixedLength_pt_bins", "1.0,2.0,5.0", "fixed-length pt bins"};

// Variable-length binning: TBI 20240113 I do it via string + tokenize + Atof(), until I figure out how to pass array directly via configurable
Configurable<bool> cUseVariableLength_mult_bins{"cUseVariableLength_mult_bins", false, "use or not variable-length multiplicity bins"};
Configurable<string> cVariableLength_mult_bins{"cVariableLength_mult_bins", "0.,100.,250.,1000.", "variable-length multiplicity bins"};
Configurable<bool> cUseVariableLength_cent_bins{"cUseVariableLength_cent_bins", false, "use or not variable-length centrality bins"};
Configurable<string> cVariableLength_cent_bins{"cVariableLength_cent_bins", "0.,10.,50.,100.", "variable-length centrality bins"};
Configurable<bool> cUseVariableLength_pt_bins{"cUseVariableLength_pt_bins", false, "use or not variable-length pt bins"};
Configurable<string> cVariableLength_pt_bins{"cVariableLength_pt_bins", "1.0,2.0,5.0", "variable-length pt bins"};
Configurable<bool> cUseVariableLength_eta_bins{"cUseVariableLength_eta_bins", false, "use or not variable-length eta bins"};
Configurable<string> cVariableLength_eta_bins{"cVariableLength_eta_bins", "-0.8,-0.4,0.0,0.4,0.8", "variable-length eta bins"};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
