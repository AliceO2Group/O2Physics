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
Configurable<string> cfTaskName{"cfTaskName", "Default task name", "set task name - use eventually to determine weights for this task"};
Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without storing and calculating anything"};
Configurable<bool> cfVerbose{"cfVerbose", false, "run or not in verbose mode (but not for function calls per particle)"};
Configurable<bool> cfVerboseForEachParticle{"cfVerboseForEachParticle", false, "run or not in verbose mode (also for function calls per particle)"};
Configurable<bool> cfDoAdditionalInsanityChecks{"cfDoAdditionalInsanityChecks", false, "do additional insanity checks at run time (this leads to small loss of performance)"};
Configurable<bool> cfInsanityCheckForEachParticle{"cfInsanityCheckForEachParticle", false, "do insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging."};
Configurable<bool> cfUseCCDB{"cfUseCCDB", true, "access personal files from CCDB or from home dir in AliEn"};
Configurable<unsigned int> cfRandomSeed{"cfRandomSeed", 0, "0 = random seed is guaranteed to be unique in space and time"};
Configurable<bool> cfUseFisherYates{"cfUseFisherYates", false, "use or not Fisher-Yates algorithm to randomize particle indices"};
Configurable<int> cfFixedNumberOfRandomlySelectedTracks{"cfFixedNumberOfRandomlySelectedTracks", -1, "Set to some integer > 0, to apply and use. Set to <=0, to ignore."};
Configurable<bool> cfUseStopwatch{"cfUseStopwatch", false, "If true, some basic info on time execution is printed, here and there."};

// Default booking:
Configurable<vector<int>> cfBookEventHistograms{"cfBookEventHistograms", {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) event histogram, ordering is the same as in enum eEventHistograms"};
Configurable<vector<int>> cfBookParticleHistograms{"cfBookParticleHistograms", {1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) particle histogram, ordering is the same as in enum eParticleHistograms"};
Configurable<vector<int>> cfBookParticleHistograms2D{"cfBookParticleHistograms2D", {1, 1}, "Book (1) or do not book (0) particle histogram, ordering is the same as in enum eParticleHistograms2D"};

// Q-vector:
Configurable<bool> cfCalculateQvectors{"cfCalculateQvectors", true, "calculate or not Q-vectors (all, also diff. ones). If I want only to fill control histograms, then I can set here false"};
// Correlations:
Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false, "calculate or not correlations"};

// Test0:
Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false, "calculate or not Test0"};
Configurable<bool> cfCalculateTest0AsFunctionOfIntegrated{"cfCalculateTest0AsFunctionOfIntegrated", true, "calculate or not Test0 as a function of integrated"};
Configurable<bool> cfCalculateTest0AsFunctionOfMultiplicity{"cfCalculateTest0AsFunctionOfMultiplicity", true, "calculate or not Test0 as a function of multiplicity"};
Configurable<bool> cfCalculateTest0AsFunctionOfCentrality{"cfCalculateTest0AsFunctionOfCentrality", true, "calculate or not Test0 as a function of centrality"};
Configurable<bool> cfCalculateTest0AsFunctionOfPt{"cfCalculateTest0AsFunctionOfPt", false, "calculate or not Test0 as a function of pt"};
Configurable<bool> cfCalculateTest0AsFunctionOfEta{"cfCalculateTest0AsFunctionOfEta", false, "calculate or not Test0 as a function of eta"};
Configurable<string> cfFileWithLabels{"cfFileWithLabels", "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Particle weights:
Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use or not phi weights"};
Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use or not pt weights"};
Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use or not eta weights"};
Configurable<bool> cfUseDiffPhiPtWeights{"cfUseDiffPhiPtWeights", false, "use or not differential phi(pt) weights"};
Configurable<bool> cfUseDiffPhiEtaWeights{"cfUseDiffPhiEtaWeights", false, "use or not differential phi(eta) weights"};
Configurable<string> cfFileWithWeights{"cfFileWithWeights", "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"

// Nested loops:
Configurable<bool> cfCalculateNestedLoops{"cfCalculateNestedLoops", false, "cross-check for all events all correlations with nested loops"};
Configurable<bool> cfCalculateCustomNestedLoops{"cfCalculateCustomNestedLoops", false, "cross-check e-b-e all correlations with custom nested loops"};
Configurable<bool> cfCalculateKineCustomNestedLoops{"cfCalculateKineCustomNestedLoops", false, "cross-check e-b-e all differential (vs. pt, eta, etc.) correlations with custom nested loops"};
Configurable<int> cfMaxNestedLoop{"cfMaxNestedLoop", -1, "if set to e.g. 4, all nested loops beyond that, e.g. 6-p and 8-p, are NOT calculated"};

// Internal validation:
Configurable<bool> cfUseInternalValidation{"cfUseInternalValidation", false, "perform internal validation using flow analysis on-the-fly"};
Configurable<unsigned int> cfnEventsInternalValidation{"cfnEventsInternalValidation", 0, "number of events simulated on-the-fly for internal validation"};
Configurable<string> cfHarmonicsOptionInternalValidation{"cfHarmonicsOptionInternalValidation", "constant", "for internal validation, set whether flow amplitudes are \"constant\" or \"correlared\""};
Configurable<bool> cfRescaleWithTheoreticalInput{"cfRescaleWithTheoreticalInput", false, "if kTRUE, all correlators are rescaled with theoretical input, so that all results in profiles are 1"};
Configurable<vector<float>> cfInternalValidationAmplitudes{"cfInternalValidationAmplitudes", {0.0, 0.1, 0.2}, "{v1, v2, v3, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of vn's is gMaxHarmonic."};
Configurable<vector<float>> cfInternalValidationPlanes{"cfInternalValidationPlanes", {0.0, 0.0, 0.0}, "{Psi1, Psi2, Psi3, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of Psin's is gMaxHarmonic."};
Configurable<vector<int>> cfMultRangeInternalValidation{"cfMultRangeInternalValidation", {1000, 1001}, "{min, max}, with convention: min <= M < max"};

// Event cuts:
Configurable<string> cfTrigger{
  "cfTrigger", "some supported trigger",
  "set here some supported trigger (kINT7, ...) "};
Configurable<bool> cfUseSel7{"cfUseSel7", false, "use for Run 2 data and MC (see official doc)"};
Configurable<bool> cfUseSel8{"cfUseSel8", false, "use for Run 3 data and MC (see official doc)"};
Configurable<string> cfCentralityEstimator{"cfCentralityEstimator", "some supported centrality estimator", "set here some supported centrality estimator (CentFT0M, CentFV0A, CentNTPV, ... for Run 3, CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2) "};
Configurable<vector<int>> cNumberOfEvents{"cNumberOfEvents", {-1, 1000000000}, "Number of events to process: {min, max}, with convention: min <= N < max"};
Configurable<vector<int>> cTotalMultiplicity{"cTotalMultiplicity", {-1, 1000000000}, "Total multiplicity range: {min, max}, with convention: min <= M < max"};
Configurable<vector<int>> cSelectedTracks{"cSelectedTracks", {-1, 1000000000}, "Selected tracks range: {min, max}, with convention: min <= M < max"};
Configurable<vector<float>> cCentrality{"cCentrality", {-10., 110.}, "Centrality range: {min, max}, with convention: min <= cent < max"};
Configurable<vector<float>> cVertex_x{"cVertex_x", {-10., 10.}, "Vertex x position range: {min, max}[cm], with convention: min <= Vx < max"};
Configurable<vector<float>> cVertex_y{"cVertex_y", {-10., 10.}, "Vertex y position range: {min, max}[cm], with convention: min <= Vy < max"};
Configurable<vector<float>> cVertex_z{"cVertex_z", {-10., 10.}, "Vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
Configurable<vector<int>> cNContributors{"cNContributors", {-1, 1000000000}, "Number of vertex contributors: {min, max}, with convention: min <= IP < max"};
Configurable<vector<float>> cImpactParameter{"cImpactParameter", {-1, 1000000000}, "Impact parameter range (can be used osnly for sim): {min, max}, with convention: min <= IP < max"};

// Particle cuts:
Configurable<float> pt_min{"pt_min", 0.2, "minimum track pt value [GeV/c]"};
Configurable<float> pt_max{"pt_max", 5.0, "maximum track pt value [GeV/c]"};

// Results histograms:
Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false, "save or not results histograms"};

// Fixed-length binning (default):
Configurable<vector<float>> cFixedLength_mult_bins{"cFixedLength_mult_bins", {2000, 0., 20000.}, "nMultBins, multMin, multMax"};
Configurable<vector<float>> cFixedLength_cent_bins{"cFixedLength_cent_bins", {110, 0., 110.}, "nCentBins, centMin, centMax"};
Configurable<vector<float>> cFixedLength_pt_bins{"cFixedLength_pt_bins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
Configurable<vector<float>> cFixedLength_eta_bins{"cFixedLength_eta_bins", {1000, -2., 2.}, "nEtaBins, etaMin, etaMax"};

// Variable-length binning: TBI 20240113 I do it via string + tokenize + Atof(), use arrays eventually as for FixedLength case above.
Configurable<bool> cUseVariableLength_mult_bins{"cUseVariableLength_mult_bins", false, "use or not variable-length multiplicity bins"};
Configurable<string> cVariableLength_mult_bins{"cVariableLength_mult_bins", "0.,100.,250.,1000.", "variable-length multiplicity bins"};
Configurable<bool> cUseVariableLength_cent_bins{"cUseVariableLength_cent_bins", false, "use or not variable-length centrality bins"};
Configurable<string> cVariableLength_cent_bins{"cVariableLength_cent_bins", "0.,10.,50.,100.", "variable-length centrality bins"};
Configurable<bool> cUseVariableLength_pt_bins{"cUseVariableLength_pt_bins", false, "use or not variable-length pt bins"};
Configurable<string> cVariableLength_pt_bins{"cVariableLength_pt_bins", "1.0,2.0,5.0", "variable-length pt bins"};
Configurable<bool> cUseVariableLength_eta_bins{"cUseVariableLength_eta_bins", false, "use or not variable-length eta bins"};
Configurable<string> cVariableLength_eta_bins{"cVariableLength_eta_bins", "-0.8,-0.4,0.0,0.4,0.8", "variable-length eta bins"};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
