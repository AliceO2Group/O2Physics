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
  Configurable<string> cfTaskIsConfiguredFromJson{"cfTaskIsConfiguredFromJson", "no", "always set manaully to \"yes\" via JSON, merely to ensure that settings are not ignored silently"};
  Configurable<string> cfTaskName{"cfTaskName", "Default task name", "set task name - use eventually to determine weights for this task"};
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without storing and calculating anything"};
  Configurable<bool> cfVerbose{"cfVerbose", false, "run or not in verbose mode (but not for simple utility functions or function calls per particle)"};
  Configurable<bool> cfVerboseUtility{"cfVerboseUtility", false, "run or not in verbose mode, also for simple utility functions (but not for function calls per particle)"};
  Configurable<bool> cfVerboseForEachParticle{"cfVerboseForEachParticle", false, "run or not in verbose mode (also for function calls per particle)"};
  Configurable<bool> cfVerboseEventCounter{"cfVerboseEventCounter", false, "print or not only event counter"};
  Configurable<bool> cfVerboseEventCut{"cfVerboseEventCut", false, "print or not which event cut didn't survive"};
  Configurable<bool> cfPlainPrintout{"cfPlainPrintout", false, "print in color or in plain (use the latter in HL)"};
  Configurable<bool> cfDoAdditionalInsanityChecks{"cfDoAdditionalInsanityChecks", false, "do additional insanity checks at run time (this leads to small loss of performance)"};
  Configurable<bool> cfInsanityCheckForEachParticle{"cfInsanityCheckForEachParticle", false, "do insanity checks at run time for each particle, at the expense of losing a lot of performance. Use only during debugging."};
  Configurable<unsigned int> cfRandomSeed{"cfRandomSeed", 0, "0 = random seed is guaranteed to be unique in space and time"};
  Configurable<bool> cfUseFisherYates{"cfUseFisherYates", false, "use or not Fisher-Yates algorithm to randomize particle indices"};
  Configurable<int> cfFixedNumberOfRandomlySelectedTracks{"cfFixedNumberOfRandomlySelectedTracks", -1, "set to some integer > 0, to apply and use. Set to <=0, to ignore."};
  Configurable<bool> cfUseStopwatch{"cfUseStopwatch", false, "if true, some basic info on time execution is printed, here and there. Very loosely, this can be used for execution time profiling."};
  Configurable<float> cfFloatingPointPrecision{"cfFloatingPointPrecision", 0.000001, "two floats are the same if TMath::Abs(f1 - f2) < fFloatingPointPrecision"};
  Configurable<int> cfSequentialBailout{"cfSequentialBailout", 0, "if fSequentialBailout > 0, then each fSequentialBailout events the function BailOut() is called. Can be used for real analysis and for IV"};
  Configurable<bool> cfUseSpecificCuts{"cfUseSpecificCuts", false, "if true, analysis-specific cuts set via configurable cfWhichSpecificCuts are applied after DefaultCuts(). "};
  Configurable<string> cfWhichSpecificCuts{"cfWhichSpecificCuts", "some supported set of analysis-specific cuts (e.g. LHC23zzh, ...)", "determine which set of analysis-specific cuts will be applied after DefaultCuts(). Use in combination with tc.fUseSpecificCuts"};

} cf_tc;

// *) QA:
struct : ConfigurableGroup {
  Configurable<bool> cfCheckUnderflowAndOverflow{"cfCheckUnderflowAndOverflow", false, "check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins (use only locally)"};
  Configurable<int> cfRebin{"cfRebin", 1, "number of bins of selected heavy 2D histograms are devided with this number"};
  Configurable<bool> cfFillQAEventHistograms2D{"cfFillQAEventHistograms2D", false, "if false, all QA 2D event histograms are not filled. if true, only the ones for which fBookQAEventHistograms2D[...] is true, are filled"};
  Configurable<vector<string>> cfBookQAEventHistograms2D{"cfBookQAEventHistograms2D", {"1-Multiplicity_vs_ReferenceMultiplicity", "1-Multiplicity_vs_NContributors", "1-Multiplicity_vs_Centrality", "1-Multiplicity_vs_Vertex_z", "1-Multiplicity_vs_Occupancy", "1-ReferenceMultiplicity_vs_NContributors", "1-ReferenceMultiplicity_vs_Centrality", "1-ReferenceMultiplicity_vs_Vertex_z", "1-ReferenceMultiplicity_vs_Occupancy", "1-NContributors_vs_Centrality", "1-NContributors_vs_Vertex_z", "1-NContributors_vs_Occupancy", "1-Centrality_vs_Vertex_z", "1-Centrality_vs_Occupancy", "0-Centrality_vs_ImpactParameter", "1-Vertex_z_vs_Occupancy", "0-MultNTracksPV_vs_MultNTracksGlobal", "1-CentFT0C_vs_CentFT0CVariant1", "1-eCentFT0C_vs_CentFT0M", "1-eCentFT0C_vs_CentFV0A", "0-CentFT0C_vs_CentNTPV", "0-CentFT0C_vs_CentNGlobal", "0-CentFT0M_vs_CentNTPV", "0-CentRun2V0M_vs_CentRun2SPDTracklets", "1-TrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange", "1-CurrentRunDuration_vs_InteractionRate"}, "book (1) or do not book (0) this QA 2D event histogram"};
  Configurable<bool> cfFillQAParticleHistograms2D{"cfFillQAParticleHistograms2D", false, "if false, all QA 2D particle histograms are not filled. if true, only the ones for which fBookQAParticleHistograms2D[...] is true, are filled"};
  Configurable<vector<string>> cfBookQAParticleHistograms2D{"cfBookQAParticleHistograms2D", {"1-Pt_vs_dcaXY"}, "book (1) or do not book (0) this QA 2D particle histogram"};
  Configurable<bool> cfFillQAParticleEventHistograms2D{"cfFillQAParticleEventHistograms2D", false, "if false, all QA 2D particle event histograms are not filled. if true, only the ones for which fBookQAParticleEventHistograms2D[...] is true, are filled"};
  Configurable<vector<string>> cfBookQAParticleEventHistograms2D{"cfBookQAParticleEventHistograms2D", {"1-CurrentRunDuration_vs_itsNCls", "1-CurrentRunDuration_vs_itsNClsNegEtaEbyE", "1-CurrentRunDuration_vs_itsNClsPosEtaEbyE", "1-CurrentRunDuration_vs_Eta0804EbyE", "1-CurrentRunDuration_vs_Eta0400EbyE", "1-CurrentRunDuration_vs_Eta0004EbyE", "1-CurrentRunDuration_vs_Eta0408EbyE", "1-CurrentRunDuration_vs_Pt0005EbyE", "1-CurrentRunDuration_vs_Pt0510EbyE", "1-CurrentRunDuration_vs_Pt1050EbyE"}, "book (1) or do not book (0) this QA 2D particle event histogram"};

} cf_qa;

// *) Event histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillEventHistograms{"cfFillEventHistograms", true, "if false, all event histograms are not filled. if true, only the ones for which fBookEventHistograms[...] is true, are filled"};
  Configurable<vector<string>> cfBookEventHistograms{"cfBookEventHistograms", {"1-NumberOfEvents", "1-TotalMultiplicity", "1-Multiplicity", "1-ReferenceMultiplicity", "1-Centrality", "1-Vertex_x", "1-Vertex_y", "1-Vertex_z", "1-NContributors", "0-ImpactParameter", "0-EventPlaneAngle", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "0-MultMCNParticlesEta08"}, "Book (1) or do not book (0) event histogram"};
} cf_eh;

// *) Event cuts:
struct : ConfigurableGroup {
  Configurable<vector<string>> cfUseEventCuts{"cfUseEventCuts", {"1-NumberOfEvents", "1-TotalMultiplicity", "1-Multiplicity", "1-ReferenceMultiplicity", "1-Centrality", "1-Vertex_x", "1-Vertex_y", "1-Vertex_z", "1-NContributors", "1-ImpactParameter", "0-EventPlaneAngle", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "0-MultMCNParticlesEta08", "0-Trigger", "0-Sel7", "1-Sel8", "1-MultiplicityEstimator", "1-ReferenceMultiplicityEstimator", "1-CentralityEstimator", "1-SelectedEvents", "1-NoSameBunchPileup", "1-IsGoodZvtxFT0vsPV", "1-IsVertexITSTPC", "1-IsVertexTOFmatched", "1-IsVertexTRDmatched", "0-NoCollInTimeRangeStrict", "0-NoCollInTimeRangeStandard", "0-NoCollInRofStrict", "0-NoCollInRofStandard", "0-NoHighMultCollInPrevRof", "0-IsGoodITSLayer3", "0-IsGoodITSLayer0123", "0-IsGoodITSLayersAll", "1-OccupancyEstimator", "1-MinVertexDistanceFromIP", "1-CentralityWeights"}, "use (1) or do not use (0) event cuts"};
  Configurable<bool> cfUseEventCutCounterAbsolute{"cfUseEventCutCounterAbsolute", false, "profile and save how many times each event cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseEventCutCounterSequential{"cfUseEventCutCounterSequential", false, "profile and save how many times each event cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<bool> cfPrintCutCounterContent{"cfPrintCutCounterContent", false, "if true, prints on the screen after each event the content of fEventCutCounterHist[*][*] (all which were booked)"};
  // Remark: Preserve below the same ordering as in enum's eEventHistograms + eEventCuts. In hyperloop, in any case this ordering is lost, because there it's alphabetical TBI 20240521 check this, after I added now std::string prefix thingie
  Configurable<vector<int>> cfNumberOfEvents{"cfNumberOfEvents", {-1, 1000000000}, "total number of events to process (whether or not they survive event cuts): {min, max}, with convention: min <= N < max"};
  Configurable<vector<int>> cfTotalMultiplicity{"cfTotalMultiplicity", {-1, 1000000000}, "total multiplicity range: {min, max}, with convention: min <= M < max"};
  Configurable<vector<float>> cfMultiplicity{"cfMultiplicity", {-1., 1000000000.}, "multiplicity (defined via cfMultiplicityEstimator) range {min, max}, with convention: min <= M < max"};
  Configurable<vector<float>> cfReferenceMultiplicity{"cfReferenceMultiplicity", {-1., 1000000000.}, "reference multiplicity (defined via cfReferenceMultiplicityEstimator)  range {min, max}, with convention: min <= M < max"};
  Configurable<vector<float>> cfCentrality{"cfCentrality", {-10., 110.}, "centrality range: {min, max}, with convention: min <= cent < max"};
  Configurable<vector<float>> cfVertex_x{"cfVertex_x", {-10., 10.}, "vertex x position range: {min, max}[cm], with convention: min <= Vx < max"};
  Configurable<vector<float>> cfVertex_y{"cfVertex_y", {-10., 10.}, "vertex y position range: {min, max}[cm], with convention: min <= Vy < max"};
  Configurable<vector<float>> cfVertex_z{"cfVertex_z", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<float> cfMinVertexDistanceFromIP{"cfMinVertexDistanceFromIP", {0.000001}, "if sqrt(vx^2+vy^2+vz^2) < cfMinVertexDistanceFromIP [cm], the event is reject. IP = nominal Interaction Point."};
  Configurable<vector<int>> cfNContributors{"cfNContributors", {2, 1000000000}, "Number of vertex contributors: {min, max}, with convention: min <= N < max"};
  Configurable<vector<float>> cfImpactParameter{"cfImpactParameter", {-1, 1000000000}, "Impact parameter range (can be used only for sim): {min, max}, with convention: min <= IP < max"};
  Configurable<vector<float>> cfEventPlaneAngle{"cfEventPlaneAngle", {-o2::constants::math::PI, o2::constants::math::TwoPI}, "Event Plane Angle range (can be used only for sim): {min, max}, with convention: min <= EP < max"};
  Configurable<vector<float>> cfOccupancy{"cfOccupancy", {-0.0001, 1000000000}, "Range for occupancy (use cfOccupancyEstimator to set specific estimator): {min, max}, with convention: min <= X < max. Important: remember that 0. has to be included, therefore I set -0.0001 by default for low edge"};
  Configurable<vector<float>> cfInteractionRate{"cfInteractionRate", {0., 1000000000.}, "Range for interaction rate: {min, max}, with convention: min <= X < max"};
  Configurable<vector<float>> cfCurrentRunDuration{"cfCurrentRunDuration", {-2, 1000000000}, "Range for current run duration (i.e. seconds since start of run) in seconds: {min, max}, with convention: min <= X < max. Only collisions taken in this range (measured from SOR) are taken for analysis"};
  Configurable<vector<float>> cfMultMCNParticlesEta08{"cfMultMCNParticlesEta08", {-1, 1000000000}, "Range for MultMCNParticlesEta08 : {min, max}, with convention: min <= X < max"};
  Configurable<string> cfTrigger{"cfTrigger", "some supported trigger (e.g. INT7 for Run 2, TVXinTRD for Run 3, etc...)", "set here some supported trigger"};
  Configurable<bool> cfUseSel7{"cfUseSel7", false, "use for Run 1 and 2 data and MC (see official doc)"};
  Configurable<bool> cfUseSel8{"cfUseSel8", false, "use for Run 3 data and MC (see official doc)"};
  Configurable<string> cfMultiplicityEstimator{"cfMultiplicityEstimator", "SelectedTracks", "all results vs. mult are calculated against this multiplicity. Can be set to SelectedTracks (calculated internally), ReferenceMultiplicity (calculated outside of my code), etc."};
  //  Configurable<string> cfReferenceMultiplicityEstimator{"cfReferenceMultiplicityEstimator", "some supported option for ref. mult. (MultFT0C, MultFV0M, MultTPC, etc.)", "Reference multiplicity, calculated outside of my code. Can be MultFT0C, MultFV0M, MultTPC, etc."};
  Configurable<string> cfReferenceMultiplicityEstimator{"cfReferenceMultiplicityEstimator", "MultFT0C", "Reference multiplicity, calculated outside of my code. Can be MultFT0C, MultFV0M, MultTPC, etc."};
  //  Configurable<string> cfCentralityEstimator{"cfCentralityEstimator", "some supported centrality estimator (e.g. CentFT0C, ...)", "set here some supported centrality estimator (CentFT0C, CentFT0M, CentFV0A, CentNTPV, ... for Run 3, and CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2 and 1) "};
  Configurable<string> cfCentralityEstimator{"cfCentralityEstimator", "CentFT0C", "set here some supported centrality estimator (CentFT0C, CentFT0M, CentFV0A, CentNTPV, ... for Run 3, and CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2 and 1) "};
  Configurable<vector<int>> cfSelectedEvents{"cfSelectedEvents", {-1, 1000000000}, "Selected number of events to process (i.e. only events which survive event cuts): {min, max}, with convention: min <= N < max"};
  Configurable<bool> cfUseNoSameBunchPileup{"cfUseNoSameBunchPileup", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsGoodZvtxFT0vsPV{"cfUseIsGoodZvtxFT0vsPV", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexITSTPC{"cfUseIsVertexITSTPC", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexTOFmatched{"cfUseIsVertexTOFmatched", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsVertexTRDmatched{"cfUseIsVertexTRDmatched", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseNoCollInTimeRangeStrict{"cfUseNoCollInTimeRangeStrict", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseNoCollInTimeRangeStandard{"cfUseNoCollInTimeRangeStandard", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseNoCollInRofStrict{"cfUseNoCollInRofStrict", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseNoCollInRofStandard{"cfUseNoCollInRofStandard", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseNoHighMultCollInPrevRof{"cfUseNoHighMultCollInPrevRof", false, "TBI 20240521 explanation"};
  Configurable<bool> cfUseIsGoodITSLayer3{"cfUseIsGoodITSLayer3", false, "TBI 20241220 explanation (or see enum)"};
  Configurable<bool> cfUseIsGoodITSLayer0123{"cfUseIsGoodITSLayer0123", false, "TBI 20241220 explanation (or see enum)"};
  Configurable<bool> cfUseIsGoodITSLayersAll{"cfUseIsGoodITSLayersAll", false, "TBI 20241220 explanation (or see enum)"};
  Configurable<string> cfOccupancyEstimator{"cfOccupancyEstimator", "FT0COccupancyInTimeRange", "set here some supported occupancy estimator (TrackOccupancyInTimeRange, FT0COccupancyInTimeRange, ..."};
} cf_ec;

// *) Particle histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillParticleHistograms{"cfFillParticleHistograms", true, "if false, all 1D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms[...] is kTRUE, are filled"};
  Configurable<vector<string>> cfBookParticleHistograms{"cfBookParticleHistograms", {"1-Phi", "1-Pt", "1-Eta", "1-Charge", "1-tpcNClsFindable", "1-tpcNClsShared", "1-itsChi2NCl", "1-tpcNClsFound", "1-tpcNClsCrossedRows", "1-itsNCls", "1-itsNClsInnerBarrel", "1-tpcCrossedRowsOverFindableCls", "1-tpcFoundOverFindableCls", "1-tpcFractionSharedCls", "1-tpcChi2NCl", "1-dcaXY", "1-dcaZ", "0-PDG"}, "Book (1) or do not book (0) particle histogram"};
  Configurable<bool> cfFillParticleHistograms2D{"cfFillParticleHistograms2D", false, "if false, all 2D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms2D[...] is kTRUE, are filled"};
  Configurable<vector<string>> cfBookParticleHistograms2D{"cfBookParticleHistograms2D", {"1-Phi_vs_Pt", "1-Phi_vs_Eta"}, "Book (1) or do not book (0) 2D particle histograms"};
} cf_ph;

// *) Particle cuts:
struct : ConfigurableGroup {
  Configurable<vector<string>> cfUseParticleCuts{"cfUseParticleCuts", {"1-Phi", "1-Pt", "1-Eta", "1-Charge", "1-tpcNClsFindable", "1-tpcNClsShared", "1-itsChi2NCl", "1-tpcNClsFound", "1-tpcNClsCrossedRows", "1-itsNCls", "1-itsNClsInnerBarrel", "1-tpcCrossedRowsOverFindableCls", "1-tpcFoundOverFindableCls", "1-tpcFractionSharedCls", "1-tpcChi2NCl", "1-dcaXY", "1-dcaZ", "1-PDG", "0-trackCutFlag", "0-trackCutFlagFb1", "0-trackCutFlagFb2", "0-isQualityTrack", "0-isPrimaryTrack", "0-isInAcceptanceTrack", "0-isGlobalTrack", "1-isPVContributor", "0-PtDependentDCAxyParameterization"}, "Use (1) or do not use (0) particle cuts"};
  Configurable<bool> cfUseParticleCutCounterAbsolute{"cfUseParticleCutCounterAbsolute", false, "profile and save how many times each particle cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseParticleCutCounterSequential{"cfUseParticleCutCounterSequential", false, "profile and save how many times each particle cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<vector<float>> cfPhi{"cfPhi", {0.0, o2::constants::math::TwoPI}, "phi range: {min, max}[rad], with convention: min <= phi < max"};
  Configurable<vector<float>> cfPt{"cfPt", {0.2, 5.0}, "pt range: {min, max}[GeV], with convention: min <= pt < max"};
  Configurable<vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};
  Configurable<vector<float>> cfCharge{"cfCharge", {-1.5, 1.5}, "particle charge. {-1.5,0} = only negative, {0,1.5} = only positive"};
  Configurable<vector<float>> cftpcNClsFindable{"cftpcNClsFindable", {-1000., 1000.}, "tpcNClsFindable range: {min, max}, with convention: min <= cftpcNClsFindable < max"};
  Configurable<vector<float>> cftpcNClsShared{"cftpcNClsShared", {-1000., 1000.}, "tpcNClsShared range: {min, max}, with convention: min <= cftpcNClsShared < max"};
  Configurable<vector<float>> cfitsChi2NCl{"cfitsChi2NCl", {-1000., 1000.}, "itsChi2NCl range: {min, max}, with convention: min <= cfitsChi2NCl < max"};
  Configurable<vector<float>> cftpcNClsFound{"cftpcNClsFound", {70., 1000.}, "tpcNClsFound range: {min, max}, with convention: min <= cftpcNClsFound < max"};
  Configurable<vector<float>> cftpcNClsCrossedRows{"cftpcNClsCrossedRows", {70., 1000.}, "tpcNClsCrossedRows range: {min, max}, with convention: min <= tpcNClsCrossedRows < max"};
  Configurable<vector<float>> cfitsNCls{"cfitsNCls", {5., 1000.}, "itsNCls range: {min, max}, with convention: min <= itsNCls < max"};
  Configurable<vector<float>> cfitsNClsInnerBarrel{"cfitsNClsInnerBarrel", {-1000., 1000.}, "itsNClsInnerBarrel range: {min, max}, with convention: min <= cfitsNClsInnerBarrel < max"};
  Configurable<vector<float>> cftpcCrossedRowsOverFindableCls{"cftpcCrossedRowsOverFindableCls", {0.8, 1000.}, "tpcCrossedRowsOverFindableCls range: {min, max}, with convention: min <= cftpcCrossedRowsOverFindableCls < max"};
  Configurable<vector<float>> cftpcFoundOverFindableCls{"cftpcFoundOverFindableCls", {0.8, 1000.}, "tpcFoundOverFindableCls range: {min, max}, with convention: min <= cftpcFoundOverFindableCls < max"};
  Configurable<vector<float>> cftpcFractionSharedCls{"cftpcFractionSharedCls", {-1000., 0.4}, "tpcFractionSharedCls range: {min, max}, with convention: min <= cftpcFractionSharedCls < max"};
  Configurable<vector<float>> cftpcChi2NCl{"cftpcChi2NCl", {-1000., 4.}, "tpcChi2NCl range: {min, max}, with convention: min <= cftpcChi2NCl < max"};
  Configurable<vector<float>> cfdcaXY{"cfdcaXY", {-2.4, 2.4}, "dcaXY range: {min, max}, with convention: min <= dcaXY < max (yes, DCA can be negative!)"};
  Configurable<vector<float>> cfdcaZ{"cfdcaZ", {-3.2, 3.2}, "dcaZ range: {min, max}, with convention: min <= dcaZ < max (yes, DCA can be negative!)"};
  Configurable<vector<float>> cfPDG{"cfPDG", {-5000., 5000.}, "PDG code"};
  Configurable<bool> cftrackCutFlag{"cftrackCutFlag", false, "general selection, particle cuts are tuned centrally. Use only in Run 3."};
  Configurable<bool> cftrackCutFlagFb1{"cftrackCutFlagFb1", false, "general selection + 1 point in ITS IB. Use only in Run 3."};
  Configurable<bool> cftrackCutFlagFb2{"cftrackCutFlagFb2", false, "general selection + 2 point in ITS IB. Use only in Run 3."};
  Configurable<bool> cfisQualityTrack{"cfisQualityTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisPrimaryTrack{"cfisPrimaryTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisInAcceptanceTrack{"cfisInAcceptanceTrack", false, "TBI 20250113 obsolete - see enum, to be removed"};
  Configurable<bool> cfisGlobalTrack{"cfisGlobalTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisPVContributor{"cfisPVContributor", false, "Has this track contributed to the collision vertex fit"};
  Configurable<string> cfPtDependentDCAxyParameterization{"cfPtDependentDCAxyParameterization", "some formula TBI add some default formula, e.g. 0.0105+0.0350/x^1.1", "set here formula for pt-dependence DCAxy cut, in the following example format 0.0105+0.0350/x^1.1"};
  // TBI 20240426 do I need to add separate support for booleans to use each specific cut?
} cf_pc;

// *) Q-vector:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateQvectors{"cfCalculateQvectors", true, "calculate or not Q-vectors (all, also diff. ones). If I want only to fill control histograms, then set here false"};
} cf_qv;

// *) Multiparticle correlations:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false, "calculate or not multiparticle correlations"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfIntegrated{"cfCalculateCorrelationsAsFunctionOfIntegrated", true, "calculate or not correlations as a function of integrated"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfMultiplicity{"cfCalculateCorrelationsAsFunctionOfMultiplicity", true, "calculate or not correlations as a function of multiplicity"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfCentrality{"cfCalculateCorrelationsAsFunctionOfCentrality", true, "calculate or not correlations as a function of centrality"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfPt{"cfCalculateCorrelationsAsFunctionOfPt", false, "calculate or not correlations as a function of pt"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfEta{"cfCalculateCorrelationsAsFunctionOfEta", false, "calculate or not correlations as a function of eta"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfOccupancy{"cfCalculateCorrelationsAsFunctionOfOccupancy", true, "calculate or not correlations as a function of occupancy"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfInteractionRate{"cfCalculateCorrelationsAsFunctionOfInteractionRate", true, "calculate or not correlations as a function of interaction rate"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfCurrentRunDuration{"cfCalculateCorrelationsAsFunctionOfCurrentRunDuration", true, "calculate or not correlations as a function of current run duration (i.e. vs. seconds since start of run)"};
  Configurable<bool> cfCalculateCorrelationsAsFunctionOfVz{"cfCalculateCorrelationsAsFunctionOfVz", true, "calculate or not correlations as a function of vertex z position"};
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
  Configurable<bool> cfCalculateTest0AsFunctionOfInteractionRate{"cfCalculateTest0AsFunctionOfInteractionRate", false, "calculate or not Test0 as a function of interaction rate"};
  Configurable<bool> cfCalculateTest0AsFunctionOfCurrentRunDuration{"cfCalculateTest0AsFunctionOfCurrentRunDuration", false, "calculate or not Test0 as a function of current run duration (i.e. vs. seconds since start of run)"};
  Configurable<bool> cfCalculateTest0AsFunctionOfVz{"cfCalculateTest0AsFunctionOfVz", false, "calculate or not Test0 as a function of vertex z position"};
  Configurable<string> cfFileWithLabels{"cfFileWithLabels", "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
  Configurable<bool> cfUseDefaultLabels{"cfUseDefaultLabels", false, "use default internally hardwired labels, only for testing purposes"};
  Configurable<string> cfWhichDefaultLabels{"cfWhichDefaultLabels", "standard", "only for testing purposes, select one set of default labels, see GetDefaultObjArrayWithLabels for supported options"};
} cf_t0;

// *) Eta separation:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateEtaSeparations{"cfCalculateEtaSeparations", false, "calculate or not 2p corr. vs. eta separations"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfIntegrated{"cfCalculateEtaSeparationsAsFunctionOfIntegrated", false, "calculate or not 2p corr. vs. eta separations ..."};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfMultiplicity{"cfCalculateEtaSeparationsAsFunctionOfMultiplicity", false, "calculate or not 2p corr. vs. eta separations as a function of multiplicity"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfCentrality{"cfCalculateEtaSeparationsAsFunctionOfCentrality", false, "calculate or not 2p corr. vs. eta separations as a function of centrality"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfPt{"cfCalculateEtaSeparationsAsFunctionOfPt", false, "calculate or not 2p corr. vs. eta separations as a function of pt"};
  // Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfEta{"cfCalculateEtaSeparationsAsFunctionOfEta", false, "this one doesn't make sense in this context"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfOccupancy{"cfCalculateEtaSeparationsAsFunctionOfOccupancy", false, "calculate or not 2p corr. vs. eta separations as a function of occupancy"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfInteractionRate{"cfCalculateEtaSeparationsAsFunctionOfInteractionRate", false, "calculate or not 2p corr. vs. eta separations as a function of interaction rate"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfCurrentRunDuration{"cfCalculateEtaSeparationsAsFunctionOfCurrentRunDuration", false, "calculate or not 2p corr. vs. eta separations as a function of current run duration (i.e. vs. seconds since start of run)"};
  Configurable<bool> cfCalculateEtaSeparationsAsFunctionOfVz{"cfCalculateEtaSeparationsAsFunctionOfVz", false, "calculate or not 2p corr. vs. eta separations as a function of vertex z position"};
  Configurable<vector<float>> cfEtaSeparationsValues{"cfEtaSeparationsValues", {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}, "Eta separation between interval A (-eta) and B (+eta)"};
  Configurable<vector<string>> cfEtaSeparationsSkipHarmonics{"cfEtaSeparationsSkipHarmonics", {"0-v1", "0-v2", "0-v3", "0-v4", "1-v5", "1-v6", "1-v7", "1-v8", "1-v9"}, "For calculation of 2p correlation with eta separation these harmonics will be skipped (if first flag = \"0-v1\", v1 will be NOT be skipped in the calculus of 2p correlations with eta separations, etc.)"};
} cf_es;

// *) Particle weights:
struct : ConfigurableGroup {
  Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use or not phi weights"};
  Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use or not pt weights"};
  Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use or not eta weights"};
  Configurable<bool> cfUseDiffPhiPtWeights{"cfUseDiffPhiPtWeights", false, "use or not differential phi(pt) weights"};
  Configurable<bool> cfUseDiffPhiEtaWeights{"cfUseDiffPhiEtaWeights", false, "use or not differential phi(eta) weights"};
  Configurable<string> cfFileWithWeights{"cfFileWithWeights", "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_pw;

// *) Centrality weights:
struct : ConfigurableGroup {
  Configurable<bool> cfUseCentralityWeights{"cfUseCentralityWeights", false, "use or not centrality weights"};
  Configurable<string> cfFileWithCentralityWeights{"cfFileWithCentralityWeights", "/home/abilandz/DatasetsO2/centralityWeights.root", "path to external ROOT file which holds centrality weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_cw;

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
  Configurable<string> cfHarmonicsOptionInternalValidation{"cfHarmonicsOptionInternalValidation", "constant", "for internal validation, supported options are \"constant\", \"correlated\" and \"persistent\""};
  Configurable<bool> cfRescaleWithTheoreticalInput{"cfRescaleWithTheoreticalInput", false, "if kTRUE, all correlators are rescaled with theoretical input, so that all results in profiles are 1"};
  Configurable<vector<float>> cfInternalValidationAmplitudes{"cfInternalValidationAmplitudes", {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09}, "{v1, v2, v3, v4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of vn's is gMaxHarmonic."};
  Configurable<vector<float>> cfInternalValidationPlanes{"cfInternalValidationPlanes", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, "{Psi1, Psi2, Psi3, Psi4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of Psin's is gMaxHarmonic."};
  Configurable<vector<int>> cfMultRangeInternalValidation{"cfMultRangeInternalValidation", {1000, 1001}, "{min, max}, with convention: min <= M < max"};
} cf_iv;

// *) Results histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false, "save or not results histograms"};

  // Fixed-length binning (default):
  Configurable<vector<float>> cfFixedLength_mult_bins{"cfFixedLength_mult_bins", {2000, 0., 20000.}, "nMultBins, multMin, multMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_cent_bins{"cfFixedLength_cent_bins", {110, 0., 110.}, "nCentBins, centMin, centMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_pt_bins{"cfFixedLength_pt_bins", {1000, 0., 10.}, "nPtBins, ptMin, ptMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_eta_bins{"cfFixedLength_eta_bins", {80, -2., 2.}, "nEtaBins, etaMin, etaMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_occu_bins{"cfFixedLength_occu_bins", {200, 0., 60000.}, "nOccuBins, occuMin, occuMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_ir_bins{"cfFixedLength_ir_bins", {1000, 0., 100.}, "nirBins, irMin, irMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_crd_bins{"cfFixedLength_crd_bins", {100000, 0., 100000.}, "ncrdBins, crdMin, crdMax (only for results histograms)"};
  Configurable<vector<float>> cfFixedLength_vz_bins{"cfFixedLength_vz_bins", {400, -20., 20.}, "nvzBins, vzMin, vzMax (only for results histograms)"};

  // Variable-length binning (per request):
  Configurable<bool> cfUseVariableLength_mult_bins{"cfUseVariableLength_mult_bins", false, "use or not variable-length multiplicity bins"};
  Configurable<vector<float>> cfVariableLength_mult_bins{"cfVariableLength_mult_bins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length multiplicity bins"};
  Configurable<bool> cfUseVariableLength_cent_bins{"cfUseVariableLength_cent_bins", false, "use or not variable-length centrality bins"};
  Configurable<vector<float>> cfVariableLength_cent_bins{"cfVariableLength_cent_bins", {0., 10., 50., 100.}, "variable-length centrality bins"};
  Configurable<bool> cfUseVariableLength_pt_bins{"cfUseVariableLength_pt_bins", true, "use or not variable-length pt bins"};
  Configurable<vector<float>> cfVariableLength_pt_bins{"cfVariableLength_pt_bins", {0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 4.00, 5.00}, "variable-length pt bins"};
  Configurable<bool> cfUseVariableLength_eta_bins{"cfUseVariableLength_eta_bins", true, "use or not variable-length eta bins"};
  Configurable<vector<float>> cfVariableLength_eta_bins{"cfVariableLength_eta_bins", {-0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8}, "variable-length eta bins"};
  Configurable<bool> cfUseVariableLength_occu_bins{"cfUseVariableLength_occu_bins", false, "use or not variable-length occupancy bins"};
  Configurable<vector<float>> cfVariableLength_occu_bins{"cfVariableLength_occu_bins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length occupancy bins"};
  Configurable<bool> cfUseVariableLength_ir_bins{"cfUseVariableLength_ir_bins", false, "use or not variable-length interaction rate bins"};
  Configurable<vector<float>> cfVariableLength_ir_bins{"cfVariableLength_ir_bins", {0., 5., 10., 50., 100., 200.}, "variable-length ineraction rate bins"};
  Configurable<bool> cfUseVariableLength_crd_bins{"cfUseVariableLength_crd_bins", false, "use or not variable-length current run duration bins"};
  Configurable<vector<float>> cfVariableLength_crd_bins{"cfVariableLength_crd_bins", {0., 5., 10., 50., 100., 500.}, "variable-length current run duration bins"};
  Configurable<bool> cfUseVariableLength_vz_bins{"cfUseVariableLength_vz_bins", false, "use or not variable-length vertex z bins"};
  Configurable<vector<float>> cfVariableLength_vz_bins{"cfVariableLength_vz_bins", {-10., -8., -6., -4, -2., -1., 0., 1., 2., 4., 6., 8., 10.}, "variable-length vertex z bins"};

} cf_res;

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
