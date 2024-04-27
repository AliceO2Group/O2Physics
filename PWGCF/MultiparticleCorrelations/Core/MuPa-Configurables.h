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

// *) Task configuration:
struct : ConfigurableGroup {
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
} cf_tc;

// *) QA:
//    TBI 20240426 add here configurables for QA histograms

// *) Event histograms:
struct : ConfigurableGroup {
  Configurable<vector<int>> cfBookEventHistograms{"cfBookEventHistograms", {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) event histogram, ordering is the same as in enum eEventHistograms"};
  // TBI 20240426 add here some more configurables for event histogtams
} cf_eh;

// *) Event cuts:
struct : ConfigurableGroup {
  Configurable<vector<int>> cfUseEventCuts{"cfUseEventCuts", {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0}, "Use (1) or do not use (0) event cuts, ordering is the same as in enums eEventHistograms + eEventCuts"};
  Configurable<string> cfTrigger{"cfTrigger", "some supported trigger", "set here some supported trigger (kINT7, ...) "};
  Configurable<bool> cfUseSel7{"cfUseSel7", false, "use for Run 2 data and MC (see official doc)"};
  Configurable<bool> cfUseSel8{"cfUseSel8", false, "use for Run 3 data and MC (see official doc)"};
  Configurable<string> cfCentralityEstimator{"cfCentralityEstimator", "some supported centrality estimator", "set here some supported centrality estimator (CentFT0M, CentFV0A, CentNTPV, ... for Run 3, CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2) "};
  Configurable<vector<int>> cfNumberOfEvents{"cfNumberOfEvents", {-1, 1000000000}, "Number of events to process: {min, max}, with convention: min <= N < max"};
  Configurable<vector<int>> cfTotalMultiplicity{"cfTotalMultiplicity", {-1, 1000000000}, "Total multiplicity range: {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfSelectedTracks{"cfSelectedTracks", {-1, 1000000000}, "Selected tracks range: {min, max}, with convention: min <= M < max"};
  Configurable<vector<float>> cfCentrality{"cfCentrality", {-10., 110.}, "Centrality range: {min, max}, with convention: min <= cent < max"};
  Configurable<vector<float>> cfVertex_x{"cfVertex_x", {-10., 10.}, "Vertex x position range: {min, max}[cm], with convention: min <= Vx < max"};
  Configurable<vector<float>> cfVertex_y{"cfVertex_y", {-10., 10.}, "Vertex y position range: {min, max}[cm], with convention: min <= Vy < max"};
  Configurable<vector<float>> cfVertex_z{"cfVertex_z", {-10., 10.}, "Vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<vector<int>> cfNContributors{"cfNContributors", {-1, 1000000000}, "Number of vertex contributors: {min, max}, with convention: min <= IP < max"};
  Configurable<vector<float>> cfImpactParameter{"cfImpactParameter", {-1, 1000000000}, "Impact parameter range (can be used osnly for sim): {min, max}, with convention: min <= IP < max"};
  // TBI 20240426 do I need to add separate support for booleans to use each specific cut?
} cf_ec;

// *) Particle histograms:
struct : ConfigurableGroup {
  Configurable<vector<int>> cfBookParticleHistograms{"cfBookParticleHistograms", {1, 1, 1, 1, 1, 1, 1}, "Book (1) or do not book (0) particle histogram, ordering is the same as in enum eParticleHistograms"};
  Configurable<vector<int>> cfBookParticleHistograms2D{"cfBookParticleHistograms2D", {1, 1}, "Book (1) or do not book (0) particle histogram, ordering is the same as in enum eParticleHistograms2D"};
} cf_ph;

// *) Particle cuts:
struct : ConfigurableGroup {
  Configurable<vector<int>> cfUseParticleCuts{"cfUseParticleCuts", {0, 1, 1, 1, 1, 1, 1, 1}, "Use (1) or do not use (0) event cuts, ordering is the same as in enums eEventHistograms + eEventCuts"};
  Configurable<vector<float>> cfPhi{"cfPhi", {0.0, TMath::TwoPi()}, "phi range: {min, max}[rad], with convention: min <= phi < max"};
  Configurable<vector<float>> cfPt{"cfPt", {0.2, 5.0}, "pt range: {min, max}[GeV], with convention: min <= pt < max"};
  Configurable<vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};
  // TBI 20240426 add suport for etpcNClsCrossedRows, eDCA_xy, eDCA_z, eDPG, ...
  // TBI 20240426 do I need to add separate support for booleans to use each specific cut?
} cf_pc;

// Q-vector:
Configurable<bool> cfCalculateQvectors{"cfCalculateQvectors", true, "calculate or not Q-vectors (all, also diff. ones). If I want only to fill control histograms, then I can set here false"};
// Correlations:
Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false, "calculate or not correlations"};

// *) Test0:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false, "calculate or not Test0"};
  Configurable<bool> cfCalculateTest0AsFunctionOfIntegrated{"cfCalculateTest0AsFunctionOfIntegrated", true, "calculate or not Test0 as a function of integrated"};
  Configurable<bool> cfCalculateTest0AsFunctionOfMultiplicity{"cfCalculateTest0AsFunctionOfMultiplicity", true, "calculate or not Test0 as a function of multiplicity"};
  Configurable<bool> cfCalculateTest0AsFunctionOfCentrality{"cfCalculateTest0AsFunctionOfCentrality", true, "calculate or not Test0 as a function of centrality"};
  Configurable<bool> cfCalculateTest0AsFunctionOfPt{"cfCalculateTest0AsFunctionOfPt", false, "calculate or not Test0 as a function of pt"};
  Configurable<bool> cfCalculateTest0AsFunctionOfEta{"cfCalculateTest0AsFunctionOfEta", false, "calculate or not Test0 as a function of eta"};
  Configurable<string> cfFileWithLabels{"cfFileWithLabels", "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_t0;

// *) Particle weights:
struct : ConfigurableGroup {
  Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use or not phi weights"};
  Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use or not pt weights"};
  Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use or not eta weights"};
  Configurable<bool> cfUseDiffPhiPtWeights{"cfUseDiffPhiPtWeights", false, "use or not differential phi(pt) weights"};
  Configurable<bool> cfUseDiffPhiEtaWeights{"cfUseDiffPhiEtaWeights", false, "use or not differential phi(eta) weights"};
  Configurable<string> cfFileWithWeights{"cfFileWithWeights", "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_pw;

// *) Nested loops:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateNestedLoops{"cfCalculateNestedLoops", false, "cross-check for all events all correlations with nested loops"};
  Configurable<bool> cfCalculateCustomNestedLoops{"cfCalculateCustomNestedLoops", false, "cross-check e-b-e all correlations with custom nested loops"};
  Configurable<bool> cfCalculateKineCustomNestedLoops{"cfCalculateKineCustomNestedLoops", false, "cross-check e-b-e all differential (vs. pt, eta, etc.) correlations with custom nested loops"};
  Configurable<int> cfMaxNestedLoop{"cfMaxNestedLoop", -1, "if set to e.g. 4, all nested loops beyond that, e.g. 6-p and 8-p, are NOT calculated"};
} cf_nl;

// *) Internal validation:
struct : ConfigurableGroup {
  Configurable<bool> cfUseInternalValidation{"cfUseInternalValidation", false, "perform internal validation using flow analysis on-the-fly"};
  Configurable<bool> cfInternalValidationForceBailout{"cfInternalValidationForceBailout", false, "force bailout (use only locally, since there is no graceful exit (yet))"};
  Configurable<unsigned int> cfnEventsInternalValidation{"cfnEventsInternalValidation", 0, "number of events simulated on-the-fly for internal validation"};
  Configurable<string> cfHarmonicsOptionInternalValidation{"cfHarmonicsOptionInternalValidation", "constant", "for internal validation, set whether flow amplitudes are \"constant\" or \"correlared\""};
  Configurable<bool> cfRescaleWithTheoreticalInput{"cfRescaleWithTheoreticalInput", false, "if kTRUE, all correlators are rescaled with theoretical input, so that all results in profiles are 1"};
  Configurable<vector<float>> cfInternalValidationAmplitudes{"cfInternalValidationAmplitudes", {0.01, 0.02, 0.03, 0.04}, "{v1, v2, v3, v4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of vn's is gMaxHarmonic."};
  Configurable<vector<float>> cfInternalValidationPlanes{"cfInternalValidationPlanes", {0.0, 0.0, 0.0, 0.0}, "{Psi1, Psi2, Psi3, Psi4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of Psin's is gMaxHarmonic."};
  Configurable<vector<int>> cfMultRangeInternalValidation{"cfMultRangeInternalValidation", {1000, 1001}, "{min, max}, with convention: min <= M < max"};
} cf_iv;

// Results histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false, "save or not results histograms"};
  // Fixed-length binning (default):
  Configurable<vector<float>> cfFixedLength_mult_bins{"cfFixedLength_mult_bins", {2000, 0., 20000.}, "nMultBins, multMin, multMax"};
  Configurable<vector<float>> cfFixedLength_cent_bins{"cfFixedLength_cent_bins", {110, 0., 110.}, "nCentBins, centMin, centMax"};
  Configurable<vector<float>> cfFixedLength_pt_bins{"cfFixedLength_pt_bins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
  Configurable<vector<float>> cfFixedLength_eta_bins{"cfFixedLength_eta_bins", {1000, -2., 2.}, "nEtaBins, etaMin, etaMax"};
  // Variable-length binning: TBI 20240113 I do it via string + tokenize + Atof(), use arrays eventually as for FixedLength case above.
  Configurable<bool> cfUseVariableLength_mult_bins{"cfUseVariableLength_mult_bins", false, "use or not variable-length multiplicity bins"};
  Configurable<string> cfVariableLength_mult_bins{"cfVariableLength_mult_bins", "0.,100.,250.,1000.", "variable-length multiplicity bins"};
  Configurable<bool> cfUseVariableLength_cent_bins{"cfUseVariableLength_cent_bins", false, "use or not variable-length centrality bins"};
  Configurable<string> cfVariableLength_cent_bins{"cfVariableLength_cent_bins", "0.,10.,50.,100.", "variable-length centrality bins"};
  Configurable<bool> cfUseVariableLength_pt_bins{"cfUseVariableLength_pt_bins", false, "use or not variable-length pt bins"};
  Configurable<string> cfVariableLength_pt_bins{"cfVariableLength_pt_bins", "1.0,2.0,5.0", "variable-length pt bins"};
  Configurable<bool> cfUseVariableLength_eta_bins{"cfUseVariableLength_eta_bins", false, "use or not variable-length eta bins"};
  Configurable<string> cfVariableLength_eta_bins{"cfVariableLength_eta_bins", "-0.8,-0.4,0.0,0.4,0.8", "variable-length eta bins"};
} cf_res;

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
