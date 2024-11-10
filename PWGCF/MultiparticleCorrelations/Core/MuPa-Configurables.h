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
  // std::string prefix = "Task configuration"; // AA: now these configurables also appear grouped on hyperloop => TBI 20240522 check if this work, and if further modifications in init are needed
  Configurable<string> cfTaskName{"cfTaskName", "Default task name", "set task name - use eventually to determine weights for this task"};
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without storing and calculating anything"};
  Configurable<bool> cfVerbose{"cfVerbose", false, "run or not in verbose mode (but not for simple utility functions or function calls per particle)"};
  Configurable<bool> cfVerboseUtility{"cfVerboseUtility", false, "run or not in verbose mode, also for simple utility functions (but not for function calls per particle)"};
  Configurable<bool> cfVerboseForEachParticle{"cfVerboseForEachParticle", false, "run or not in verbose mode (also for function calls per particle)"};
  Configurable<bool> cfDoAdditionalInsanityChecks{"cfDoAdditionalInsanityChecks", false, "do additional insanity checks at run time (this leads to small loss of performance)"};
  Configurable<bool> cfInsanityCheckForEachParticle{"cfInsanityCheckForEachParticle", false, "do insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging."};
  Configurable<bool> cfUseCCDB{"cfUseCCDB", true, "if requested, access personal files from CCDB (true) or from home dir in AliEn (false)"};
  Configurable<unsigned int> cfRandomSeed{"cfRandomSeed", 0, "0 = random seed is guaranteed to be unique in space and time"};
  Configurable<bool> cfUseFisherYates{"cfUseFisherYates", false, "use or not Fisher-Yates algorithm to randomize particle indices"};
  Configurable<int> cfFixedNumberOfRandomlySelectedTracks{"cfFixedNumberOfRandomlySelectedTracks", -1, "set to some integer > 0, to apply and use. Set to <=0, to ignore."};
  Configurable<bool> cfUseStopwatch{"cfUseStopwatch", false, "if true, some basic info on time execution is printed, here and there. Very loosely, this can be used for execution time profiling."};
  Configurable<float> cfFloatingPointPrecision{"cfFloatingPointPrecision", 0.000001, "two floats are the same if TMath::Abs(f1 - f2) < fFloatingPointPrecision"};
} cf_tc;

// *) QA:
struct : ConfigurableGroup {
  Configurable<bool> cfCheckUnderflowAndOverflow{"cfCheckUnderflowAndOverflow", false, "check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins (use only locally)"};
  Configurable<bool> cfFillQAEventHistograms2D{"cfFillQAEventHistograms2D", false, "if false, all QA 2D event histograms are not filled. if true, only the ones for which fBookQAEventHistograms2D[...] is true, are filled"};
  Configurable<vector<string>> cfBookQAEventHistograms2D{"cfBookQAEventHistograms2D", {"MultTPC_vs_NContributors-1", "Vertex_z_vs_MultTPC-1", "Vertex_z_vs_NContributors-1", "CentFT0M_vs_CentNTPV-1", "CentRun2V0M_vs_CentRun2SPDTracklets-1", "CentRun2V0M_vs_NContributors-1", "TrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange-1"}, "book (1) or do not book (0) this QA 2D event histogram"};
  Configurable<bool> cfFillQAParticleHistograms2D{"cfFillQAParticleHistograms2D", false, "if false, all QA 2D particle histograms are not filled. if true, only the ones for which fBookQAParticleHistograms2D[...] is true, are filled"};
  Configurable<vector<string>> cfBookQAParticleHistograms2D{"cfBookQAParticleHistograms2D", {"Pt_vs_dcaXY-1"}, "book (1) or do not book (0) this QA 2D particle histogram"};
} cf_qa;

// *) Event histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillEventHistograms{"cfFillEventHistograms", true, "if false, all event histograms are not filled. if true, only the ones for which fBookEventHistograms[...] is true, are filled"};
  Configurable<vector<string>> cfBookEventHistograms{"cfBookEventHistograms", {"NumberOfEvents-1", "TotalMultiplicity-1", "SelectedTracks-1", "MultFV0M-1", "MultFT0M-1", "MultTPC-1", "MultNTracksPV-1", "MultTracklets-1", "Centrality-1", "Vertex_x-1", "Vertex_y-1", "Vertex_z-1", "NContributors-1", "ImpactParameter-1", "Occupancy-1"}, "Book (1) or do not book (0) event histogram"};
} cf_eh;

// *) Event cuts:
struct : ConfigurableGroup {
  Configurable<vector<string>> cfUseEventCuts{"cfUseEventCuts", {"NumberOfEvents-1", "TotalMultiplicity-1", "SelectedTracks-1", "MultFV0M-1", "MultFT0M-1", "MultTPC-1", "MultNTracksPV-1", "MultTracklets-1", "Centrality-1", "Vertex_x-1", "Vertex_y-1", "Vertex_z-1", "NContributors-1", "ImpactParameter-1", "Occupancy-1", "Trigger-0", "Sel7-1", "Sel8-1", "CentralityEstimator-1", "SelectedEvents-1", "NoSameBunchPileup-1", "IsGoodZvtxFT0vsPV-1", "IsVertexITSTPC-1", "IsVertexTOFmatched-1", "IsVertexTRDmatched-1", "OccupancyEstimator-1"}, "use (1) or do not use (0) event cuts"};
  Configurable<bool> cfUseEventCutCounterAbsolute{"cfUseEventCutCounterAbsolute", false, "profile and save how many times each event cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseEventCutCounterSequential{"cfUseEventCutCounterSequential", false, "profile and save how many times each event cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<bool> cfPrintCutCounterContent{"cfPrintCutCounterContent", false, "if true, prints on the screen after each event the content of fEventCutCounterHist[*][*] (all which were booked)"};
  // Remark: Preserve below the same ordering as in enum's eEventHistograms + eEventCuts. In hyperloop, in any case this ordering is lost, because there it's alphabetical TBI 20240521 check this, after I added now std::string prefix thingie
  Configurable<vector<int>> cfNumberOfEvents{"cfNumberOfEvents", {-1, 1000000000}, "total number of events to process (whether or not they survive event cuts): {min, max}, with convention: min <= N < max"};
  Configurable<vector<int>> cfTotalMultiplicity{"cfTotalMultiplicity", {-1, 1000000000}, "total multiplicity range: {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfSelectedTracks{"cfSelectedTracks", {-1, 1000000000}, "selected tracks range: {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfMultFV0M{"cfMultFV0M", {-1, 1000000000}, "MultFV0M range {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfMultFT0M{"cfMultFT0M", {-1, 1000000000}, "MultFT0M range {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfMultTPC{"cfMultTPC", {-1, 1000000000}, "MultTPC range {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfMultNTracksPV{"cfMultNTracksPV", {-1, 1000000000}, "MultNTracksPV range {min, max}, with convention: min <= M < max"};
  Configurable<vector<int>> cfMultTracklets{"cfMultTracklets", {-1, 1000000000}, "MultTracklets range {min, max}, with convention: min <= M < max"};
  Configurable<vector<float>> cfCentrality{"cfCentrality", {-10., 110.}, "centrality range: {min, max}, with convention: min <= cent < max"};
  Configurable<vector<float>> cfVertex_x{"cfVertex_x", {-10., 10.}, "vertex x position range: {min, max}[cm], with convention: min <= Vx < max"};
  Configurable<vector<float>> cfVertex_y{"cfVertex_y", {-10., 10.}, "vertex y position range: {min, max}[cm], with convention: min <= Vy < max"};
  Configurable<vector<float>> cfVertex_z{"cfVertex_z", {-10., 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<vector<int>> cfNContributors{"cfNContributors", {-1, 1000000000}, "Number of vertex contributors: {min, max}, with convention: min <= N < max"};
  Configurable<vector<float>> cfImpactParameter{"cfImpactParameter", {-1, 1000000000}, "Impact parameter range (can be used only for sim): {min, max}, with convention: min <= IP < max"};
  Configurable<vector<float>> cfOccupancy{"cfOccupancy", {-2, 1000000000}, "Range for occupancy (use cfOccupancyEstimator to set specific estimator): {min, max}, with convention: min <= X < max"};
  Configurable<string> cfTrigger{"cfTrigger", "some supported trigger", "set here some supported trigger (kINT7, ...) "};
  Configurable<bool> cfUseSel7{"cfUseSel7", false, "use for Run 1 and 2 data and MC (see official doc)"};
  Configurable<bool> cfUseSel8{"cfUseSel8", false, "use for Run 3 data and MC (see official doc)"};
  Configurable<string> cfCentralityEstimator{"cfCentralityEstimator", "some supported centrality estimator", "set here some supported centrality estimator (CentFT0M, CentFV0A, CentNTPV, ... for Run 3, and CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2 and 1) "};
  Configurable<vector<int>> cfSelectedEvents{"cfSelectedEvents", {-1, 1000000000}, "Selected number of events to process (i.e. only events which survive event cuts): {min, max}, with convention: min <= N < max"};
  Configurable<bool> cfUseNoSameBunchPileup{"cfUseNoSameBunchPileup", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsGoodZvtxFT0vsPV{"cfUseIsGoodZvtxFT0vsPV", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexITSTPC{"cfUseIsVertexITSTPC", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexTOFmatched{"cfUseIsVertexTOFmatched", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexTRDmatched{"cfUseIsVertexTRDmatched", false, "TBI 20240521 explanation"};
  Configurable<string> cfOccupancyEstimator{"cfOccupancyEstimator", "some supported occupancy estimator", "set here some supported occupancy estimator (TrackOccupancyInTimeRange, FT0COccupancyInTimeRange, ..."};
} cf_ec;

// *) Particle histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillParticleHistograms{"cfFillParticleHistograms", true, "if false, all 1D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms[...] is kTRUE, are filled"};
  Configurable<vector<string>> cfBookParticleHistograms{"cfBookParticleHistograms", {"Phi-1", "Pt-1", "Eta-1", "Charge-1", "tpcNClsFindable-1", "tpcNClsShared-1", "tpcNClsFound-1", "tpcNClsCrossedRows-1", "itsNCls-1", "itsNClsInnerBarrel-1", "tpcCrossedRowsOverFindableCls-1", "tpcFoundOverFindableCls-1", "tpcFractionSharedCls-1", "dcaXY-1", "dcaZ-1", "PDG-1"}, "Book (1) or do not book (0) particle histogram"};
  Configurable<bool> cfFillParticleHistograms2D{"cfFillParticleHistograms2D", true, "if false, all 2D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms2D[...] is kTRUE, are filled"};
  Configurable<vector<string>> cfBookParticleHistograms2D{"cfBookParticleHistograms2D", {"Phi_vs_Pt-1", "Phi_vs_Eta-1"}, "Book (1) or do not book (0) this 2D particle histogram"};
} cf_ph;

// *) Particle cuts:
struct : ConfigurableGroup {
  Configurable<vector<string>> cfUseParticleCuts{"cfUseParticleCuts", {"Phi-1", "Pt-1", "Eta-1", "Charge-1", "tpcNClsFindable-1", "tpcNClsShared-1", "tpcNClsFound-1", "tpcNClsCrossedRows-1", "itsNCls-1", "itsNClsInnerBarrel-1", "tpcCrossedRowsOverFindableCls-1", "tpcFoundOverFindableCls-1", "tpcFractionSharedCls-1", "dcaXY-1", "dcaZ-1", "PDG-1", "trackCutFlagFb1-0", "trackCutFlagFb2-0", "isQualityTrack-0", "isPrimaryTrack-0", "isInAcceptanceTrack-0", "isGlobalTrack-0", "PtDependentDCAxyParameterization-0"}, "Use (1) or do not use (0) particle cuts"};
  Configurable<bool> cfUseParticleCutCounterAbsolute{"cfUseParticleCutCounterAbsolute", false, "profile and save how many times each particle cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseParticleCutCounterSequential{"cfUseParticleCutCounterSequential", false, "profile and save how many times each particle cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<vector<float>> cfPhi{"cfPhi", {0.0, TMath::TwoPi()}, "phi range: {min, max}[rad], with convention: min <= phi < max"};
  Configurable<vector<float>> cfPt{"cfPt", {0.2, 5.0}, "pt range: {min, max}[GeV], with convention: min <= pt < max"};
  Configurable<vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cfCharge{"cfCharge", {-1.5, 1.5}, "particle charge. {-1.5,0} = only negative, {0,1.5} = only positive"};
  Configurable<vector<float>> cftpcNClsFindable{"cftpcNClsFindable", {-1000., 1000.}, "tpcNClsFindable range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcNClsShared{"cftpcNClsShared", {-1000., 1000.}, "tpcNClsShared range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcNClsFound{"cftpcNClsFound", {-1000., 1000.}, "tpcNClsFound range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcNClsCrossedRows{"cftpcNClsCrossedRows", {-1000., 1000.}, "tpcNClsCrossedRows range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cfitsNCls{"cfitsNCls", {-1000., 1000.}, "itsNCls range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cfitsNClsInnerBarrel{"cfitsNClsInnerBarrel", {-1000., 1000.}, "itsNClsInnerBarrel range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcCrossedRowsOverFindableCls{"cftpcCrossedRowsOverFindableCls", {-1000., 1000.}, "tpcCrossedRowsOverFindableCls range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcFoundOverFindableCls{"cftpcFoundOverFindableCls", {-1000., 1000.}, "tpcFoundOverFindableCls range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cftpcFractionSharedCls{"cftpcFractionSharedCls", {-1000., 1000.}, "tpcFractionSharedCls range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cfdcaXY{"cfdcaXY", {-1000., 1000.}, "dcaXY range: {min, max}, with convention: min <= dcaXY < max (yes, DCA can be negative!)"};
  Configurable<vector<float>> cfdcaZ{"cfdcaZ", {-1000., 1000.}, "dcaZ range: {min, max}, with convention: min <= dcaZ < max (yes, DCA can be negative!)"};
  Configurable<vector<float>> cfPDG{"cfPDG", {-5000., 5000.}, "PDG code"};
  Configurable<bool> cftrackCutFlagFb1{"cftrackCutFlagFb1", false, "TBI 20240510 add description"};
  Configurable<bool> cftrackCutFlagFb2{"cftrackCutFlagFb2", false, "TBI 20240510 add description"};
  Configurable<bool> cfisQualityTrack{"cfisQualityTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisPrimaryTrack{"cfisPrimaryTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisInAcceptanceTrack{"cfisInAcceptanceTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisGlobalTrack{"cfisGlobalTrack", false, "TBI 20240510 add description"};
  Configurable<string> cfPtDependentDCAxyParameterization{"cfPtDependentDCAxyParameterization", "some formula TBI add some default formula, e.g. 0.0105+0.0350/x^1.1", "set here formula for pt-dependence DCAxy cut, in the following example format 0.0105+0.0350/x^1.1"};
  // TBI 20240426 do I need to add separate support for booleans to use each specific cut?
} cf_pc;

// *) Q-vector:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateQvectors{"cfCalculateQvectors", false, "calculate or not Q-vectors (all, also diff. ones). If I want only to fill control histograms, then set here false"};
} cf_qv;

// *) Multiparticle correlations:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false, "calculate or not multiparticle correlations"};
} cf_mupa;

// *) Test0:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false, "calculate or not Test0"};
  Configurable<bool> cfCalculateTest0AsFunctionOfIntegrated{"cfCalculateTest0AsFunctionOfIntegrated", false, "calculate or not Test0 as a function of integrated"};
  Configurable<bool> cfCalculateTest0AsFunctionOfMultiplicity{"cfCalculateTest0AsFunctionOfMultiplicity", false, "calculate or not Test0 as a function of multiplicity"};
  Configurable<bool> cfCalculateTest0AsFunctionOfCentrality{"cfCalculateTest0AsFunctionOfCentrality", false, "calculate or not Test0 as a function of centrality"};
  Configurable<bool> cfCalculateTest0AsFunctionOfPt{"cfCalculateTest0AsFunctionOfPt", false, "calculate or not Test0 as a function of pt"};
  Configurable<bool> cfCalculateTest0AsFunctionOfEta{"cfCalculateTest0AsFunctionOfEta", false, "calculate or not Test0 as a function of eta"};
  Configurable<bool> cfCalculateTest0AsFunctionOfOccupancy{"cfCalculateTest0AsFunctionOfOccupancy", false, "calculate or not Test0 as a function of occupancy"};
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

// *) Toy NUA:
struct : ConfigurableGroup {
  Configurable<vector<int>> cfApplyNUAPDF{"cfApplyNUAPDF", {0, 0, 0}, "Apply (1) or do not apply (0) NUA on variable, ordering is the same as in enum eNUAPDF (phi, pt, eta)"};
  Configurable<vector<int>> cfUseDefaultNUAPDF{"cfUseDefaultNUAPDF", {1, 1, 1}, "Use (1) or do not use (0) default NUA profile, ordering is the same as in enum eNUAPDF (phi, pt, eta)"};
  Configurable<vector<string>> cfCustomNUAPDFHistNames{"cfCustomNUAPDFHistNames", {"a", "bb", "ccc"}, "the names of histograms holding custom NUA in an external file."};
  Configurable<string> cfFileWithCustomNUA{"cfFileWithCustomNUA", "/home/abilandz/DatasetsO2/customNUA.root", "path to external ROOT file which holds all histograms with custom NUA"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_nua;

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

// *) Results histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false, "save or not results histograms"};

  // Fixed-length binning (default):
  Configurable<vector<float>> cfFixedLength_mult_bins{"cfFixedLength_mult_bins", {2000, 0., 20000.}, "nMultBins, multMin, multMax"};
  Configurable<vector<float>> cfFixedLength_cent_bins{"cfFixedLength_cent_bins", {110, 0., 110.}, "nCentBins, centMin, centMax"};
  Configurable<vector<float>> cfFixedLength_pt_bins{"cfFixedLength_pt_bins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};
  Configurable<vector<float>> cfFixedLength_eta_bins{"cfFixedLength_eta_bins", {1000, -2., 2.}, "nEtaBins, etaMin, etaMax"};
  Configurable<vector<float>> cfFixedLength_occu_bins{"cfFixedLength_occu_bins", {400, 0., 4000.}, "nOccuBins, occuMin, occuMax"};

  // Variable-length binning (per request):
  Configurable<bool> cfUseVariableLength_mult_bins{"cfUseVariableLength_mult_bins", false, "use or not variable-length multiplicity bins"};
  Configurable<vector<float>> cfVariableLength_mult_bins{"cfVariableLength_mult_bins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length multiplicity bins"};
  Configurable<bool> cfUseVariableLength_cent_bins{"cfUseVariableLength_cent_bins", false, "use or not variable-length centrality bins"};
  Configurable<vector<float>> cfVariableLength_cent_bins{"cfVariableLength_cent_bins", {0., 10., 50., 100.}, "variable-length centrality bins"};
  Configurable<bool> cfUseVariableLength_pt_bins{"cfUseVariableLength_pt_bins", false, "use or not variable-length pt bins"};
  Configurable<vector<float>> cfVariableLength_pt_bins{"cfVariableLength_pt_bins", {0.20, 0.30, 0.40, 0.65, 1.00, 2.00, 5.00}, "variable-length pt bins"};
  Configurable<bool> cfUseVariableLength_eta_bins{"cfUseVariableLength_eta_bins", false, "use or not variable-length eta bins"};
  Configurable<vector<float>> cfVariableLength_eta_bins{"cfVariableLength_eta_bins", {3.0, -1.0, -0.4, 0.0, 0.4, 1.0, 3.0}, "variable-length eta bins"};
  Configurable<bool> cfUseVariableLength_occu_bins{"cfUseVariableLength_occu_bins", false, "use or not variable-length occupancy bins"};
  Configurable<vector<float>> cfVariableLength_occu_bins{"cfVariableLength_occu_bins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length occupancy bins"};

} cf_res;

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
