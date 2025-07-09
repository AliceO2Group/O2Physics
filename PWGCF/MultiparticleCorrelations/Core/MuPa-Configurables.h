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

/// \file MuPa-Configurables.h
/// \brief ... TBI 20250425
/// \author Ante.Bilandzic@cern.ch

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_

// ...
#include <string>
#include <vector>

// *) Task configuration:
struct : ConfigurableGroup {
  // std::string prefix = "Task configuration"; // AA: now these configurables also appear grouped on hyperloop => TBI 20240522 check if this work, and if further modifications in init are needed
  Configurable<std::string> cfTaskIsConfiguredFromJson{"cfTaskIsConfiguredFromJson", "no", "always set manaully to \"yes\" via JSON, merely to ensure that settings are not ignored silently"};
  Configurable<std::string> cfTaskName{"cfTaskName", "Default task name", "set task name - use eventually to determine weights for this task"};
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
  Configurable<float> cfFloatingPointPrecision{"cfFloatingPointPrecision", 0.000001, "two floats are the same if abs(f1 - f2) < fFloatingPointPrecision"};
  Configurable<int> cfSequentialBailout{"cfSequentialBailout", 0, "if fSequentialBailout > 0, then each fSequentialBailout events the function BailOut() is called. Can be used for real analysis and for IV"};
  Configurable<bool> cfUseSpecificCuts{"cfUseSpecificCuts", false, "if true, analysis-specific cuts set via configurable cfWhichSpecificCuts are applied after DefaultCuts(). "};
  Configurable<std::string> cfWhichSpecificCuts{"cfWhichSpecificCuts", "some supported set of analysis-specific cuts (e.g. LHC23zzh, ...)", "determine which set of analysis-specific cuts will be applied after DefaultCuts(). Use in combination with tc.fUseSpecificCuts"};
  Configurable<std::string> cfSkipTheseRuns{"cfSkipTheseRuns", "", "Set here via comma-separated list which runs will be skipped during hl analysis (a.k.a. \"bad runs\"). Leave empty to ignore. Example format and list for LHC23zzh: \"544116,544091\""};
  Configurable<bool> cfUseSetBinLabel{"cfUseSetBinLabel", false, "until hist->SetBinLabel(...) large memory consumption is resolved, for each histogram dump all that info in the y-axis title. See also local executable PostprocessLabels.C, where I do the final bin labeling offline"};
  Configurable<bool> cfUseClone{"cfUseClone", false, "until hist->Clone(...) large memory consumption is resolved, do not use cloning. See ROOT Forum thread."};
  Configurable<bool> cfUseFormula{"cfUseFormula", false, "until TFormula large memory consumption is resolved, do not use this class. See ROOT Forum thread."};
  Configurable<bool> cfUseDatabasePDG{"cfUseDatabasePDG", false, "When enabled, there is a standard memory blow-up."};
} cf_tc;

// *) QA:
struct : ConfigurableGroup {
  Configurable<bool> cfCheckUnderflowAndOverflow{"cfCheckUnderflowAndOverflow", false, "check and bail out if in event and particle histograms there are entries which went to underflow or overflow bins (use only locally)"};
  Configurable<int> cfRebin{"cfRebin", 10, "number of bins of selected heavy 2D histograms are devided with this number"};
  Configurable<bool> cfFillQAEventHistograms2D{"cfFillQAEventHistograms2D", false, "if false, all QA 2D event histograms are not filled. if true, only the ones for which fBookQAEventHistograms2D[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookQAEventHistograms2D{"cfBookQAEventHistograms2D", {"1-Multiplicity_vs_ReferenceMultiplicity", "1-Multiplicity_vs_NContributors", "1-Multiplicity_vs_Centrality", "1-Multiplicity_vs_VertexZ", "1-Multiplicity_vs_Occupancy", "1-Multiplicity_vs_InteractionRate", "1-ReferenceMultiplicity_vs_NContributors", "1-ReferenceMultiplicity_vs_Centrality", "1-ReferenceMultiplicity_vs_VertexZ", "1-ReferenceMultiplicity_vs_Occupancy", "1-ReferenceMultiplicity_vs_InteractionRate", "1-NContributors_vs_Centrality", "1-NContributors_vs_VertexZ", "1-NContributors_vs_Occupancy", "1-NContributors_vs_InteractionRate", "1-Centrality_vs_VertexZ", "1-Centrality_vs_Occupancy", "0-Centrality_vs_ImpactParameter", "1-Centrality_vs_InteractionRate", "1-VertexZ_vs_Occupancy", "1-VertexZ_vs_InteractionRate", "0-MultNTracksPV_vs_MultNTracksGlobal", "1-CentFT0C_vs_CentFT0CVariant1", "1-CentFT0C_vs_CentFT0M", "1-CentFT0C_vs_CentFV0A", "0-CentFT0C_vs_CentNTPV", "0-CentFT0C_vs_CentNGlobal", "0-CentFT0M_vs_CentNTPV", "0-CentRun2V0M_vs_CentRun2SPDTracklets", "1-TrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange", "1-CurrentRunDuration_vs_InteractionRate", "1-Multiplicity_vs_FT0CAmplitudeOnFoundBC", "1-CentFT0C_vs_FT0CAmplitudeOnFoundBC", "1-Centrality_vs_CentralitySim"}, "book (1) or do not book (0) this QA 2D event histogram"};
  Configurable<bool> cfFillQAParticleHistograms2D{"cfFillQAParticleHistograms2D", false, "if false, all QA 2D particle histograms are not filled. if true, only the ones for which fBookQAParticleHistograms2D[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookQAParticleHistograms2D{"cfBookQAParticleHistograms2D", {"1-Pt_vs_dcaXY"}, "book (1) or do not book (0) this QA 2D particle histogram"};
  Configurable<bool> cfFillQAParticleEventHistograms2D{"cfFillQAParticleEventHistograms2D", false, "if false, all QA 2D particle event histograms are not filled. if true, only the ones for which fBookQAParticleEventHistograms2D[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookQAParticleEventHistograms2D{"cfBookQAParticleEventHistograms2D", {"1-CurrentRunDuration_vs_itsNCls", "1-CurrentRunDuration_vs_itsNClsNegEtaEbyE", "1-CurrentRunDuration_vs_itsNClsPosEtaEbyE", "1-CurrentRunDuration_vs_Eta0804EbyE", "1-CurrentRunDuration_vs_Eta0400EbyE", "1-CurrentRunDuration_vs_Eta0004EbyE", "1-CurrentRunDuration_vs_Eta0408EbyE", "1-CurrentRunDuration_vs_Pt0005EbyE", "1-CurrentRunDuration_vs_Pt0510EbyE", "1-CurrentRunDuration_vs_Pt1050EbyE"}, "book (1) or do not book (0) this QA 2D particle event histogram"};

  Configurable<bool> cfFillQACorrelationsVsHistograms2D{"cfFillQACorrelationsVsHistograms2D", false, "if false, all QA 2D histograms of this category are not filled. if true, only the ones for which fBookQACorrelationsVsHistograms2D[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookQACorrelationsVsHistograms2D{"cfBookQACorrelationsVsHistograms2D", {"1-Correlations_vs_Multiplicity", "1-Correlations_vs_ReferenceMultiplicity", "1-Correlations_vs_Centrality", "1-Correlations_vs_Phi", "1-Correlations_vs_Pt", "1-Correlations_vs_Eta", "1-Correlations_vs_Charge", "1-Correlations_vs_tpcNClsFindable", "1-Correlations_vs_tpcNClsShared", "1-Correlations_vs_itsChi2NCl", "1-Correlations_vs_tpcNClsFound", "1-Correlations_vs_tpcNClsCrossedRows", "1-Correlations_vs_itsNCls", "1-Correlations_vs_itsNClsInnerBarrel", "1-Correlations_vs_tpcCrossedRowsOverFindableCls", "1-Correlations_vs_tpcFoundOverFindableCls", "1-Correlations_vs_tpcFractionSharedCls", "1-Correlations_vs_tpcChi2NCl", "1-Correlations_vs_dcaXY", "1-Correlations_vs_dcaZ"}, "book (1) or do not book (0) this QA 2D histogram"};
  Configurable<std::vector<int>> cfQACorrelationsVsHistogramsMinMaxHarmonic{"cfQACorrelationsVsHistogramsMinMaxHarmonic", {1, 5}, "harmonics are filled for min <= harmonic < max"};

  Configurable<bool> cfFillQACorrelationsVsInteractionRateVsProfiles2D{"cfFillQACorrelationsVsInteractionRateVsProfiles2D", false, "if false, all QA 2D profiles of this category are not filled. if true, only the ones for which fBookQACorrelationsVsInteractionRateVsProfiles2D[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookQACorrelationsVsInteractionRateVsProfiles2D{"cfBookQACorrelationsVsInteractionRateVsProfiles2D", {"1-CorrVsIR_vs_CurrentRunDuration", "1-CorrVsIR_vs_Multiplicity", "1-CorrVsIR_vs_ReferenceMultiplicity", "1-CorrVsIR_vs_Centrality", "1-CorrVsIR_vs_MeanPhi", "1-CorrVsIR_vs_SigmaMeanPhi", "1-CorrVsIR_vs_MeanPt", "1-CorrVsIR_vs_SigmaMeanPt", "1-CorrVsIR_vs_MeanEta", "1-CorrVsIR_vs_SigmaMeanEta"}, "book (1) or do not book (0) this QA 2D profile"};
  Configurable<std::vector<int>> cfQACorrelationsVsInteractionRateVsProfilesMinMaxHarmonic{"cfQACorrelationsVsInteractionRateVsProfilesMinMaxHarmonic", {1, 5}, "harmonics are filled for min <= harmonic < max"};

} cf_qa;

// *) Event histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillEventHistograms{"cfFillEventHistograms", true, "if false, all event histograms are not filled. if true, only the ones for which fBookEventHistograms[...] is true, are filled"};
  Configurable<std::vector<std::string>> cfBookEventHistograms{"cfBookEventHistograms", {"1-NumberOfEvents", "1-TotalMultiplicity", "1-Multiplicity", "1-ReferenceMultiplicity", "1-Centrality", "1-VertexX", "1-VertexY", "1-VertexZ", "1-NContributors", "0-ImpactParameter", "0-EventPlaneAngle", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "0-MultMCNParticlesEta08"}, "Book (1) or do not book (0) event histogram"};
} cf_eh;

// *) Event cuts:
struct : ConfigurableGroup {
  Configurable<std::vector<std::string>> cfUseEventCuts{"cfUseEventCuts", {"1-NumberOfEvents", "1-TotalMultiplicity", "1-Multiplicity", "1-ReferenceMultiplicity", "1-Centrality", "1-VertexX", "1-VertexY", "1-VertexZ", "1-NContributors", "1-ImpactParameter", "0-EventPlaneAngle", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "0-MultMCNParticlesEta08", "0-Trigger", "0-Sel7", "1-Sel8", "1-MultiplicityEstimator", "1-ReferenceMultiplicityEstimator", "1-CentralityEstimator", "1-SelectedEvents", "1-NoSameBunchPileup", "1-IsGoodZvtxFT0vsPV", "1-IsVertexITSTPC", "1-IsVertexTOFmatched", "1-IsVertexTRDmatched", "0-NoCollInTimeRangeStrict", "0-NoCollInTimeRangeStandard", "0-NoCollInRofStrict", "0-NoCollInRofStandard", "0-NoHighMultCollInPrevRof", "0-IsGoodITSLayer3", "0-IsGoodITSLayer0123", "0-IsGoodITSLayersAll", "1-OccupancyEstimator", "1-MinVertexDistanceFromIP", "0-NoPileupTPC", "0-NoPileupFromSPD", "0-NoSPDOnVsOfPileup", "1-RefMultVsNContrUp", "1-RefMultVsNContrLow", "1-CentralityCorrelationsCut", "0-FT0Bad", "0-ITSBad", "0-ITSLimAccMCRepr", "0-TPCBadTracking", "0-TPCLimAccMCRepr", "0-TPCBadPID", "1-CentralityWeights"}, "use (1) or do not use (0) event cuts"};
  Configurable<bool> cfUseEventCutCounterAbsolute{"cfUseEventCutCounterAbsolute", false, "profile and save how many times each event cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseEventCutCounterSequential{"cfUseEventCutCounterSequential", false, "profile and save how many times each event cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<bool> cfPrintCutCounterContent{"cfPrintCutCounterContent", false, "if true, prints on the screen after each event the content of fEventCutCounterHist[*][*] (all which were booked)"};
  // Remark: Preserve below the same ordering as in enum's eEventHistograms + eEventCuts. In hyperloop, in any case this ordering is lost, because there it's alphabetical TBI 20240521 check this, after I added now std::string prefix thingie
  Configurable<std::vector<int>> cfNumberOfEvents{"cfNumberOfEvents", {-1, 1000000000}, "total number of events to process (whether or not they survive event cuts): {min, max}, with convention: min <= N < max"};
  Configurable<std::vector<int>> cfTotalMultiplicity{"cfTotalMultiplicity", {-1, 1000000000}, "total multiplicity range: {min, max}, with convention: min <= M < max"};
  Configurable<std::vector<float>> cfMultiplicity{"cfMultiplicity", {-1., 1000000000.}, "multiplicity (defined via cfMultiplicityEstimator) range {min, max}, with convention: min <= M < max"};
  Configurable<std::vector<float>> cfReferenceMultiplicity{"cfReferenceMultiplicity", {-1., 1000000000.}, "reference multiplicity (defined via cfReferenceMultiplicityEstimator)  range {min, max}, with convention: min <= M < max"};
  Configurable<std::vector<float>> cfCentrality{"cfCentrality", {-10., 110.}, "centrality range: {min, max}, with convention: min <= cent < max"};
  Configurable<std::vector<float>> cfVertexX{"cfVertexX", {-10., 10.}, "vertex x position range: {min, max}[cm], with convention: min <= Vx < max"};
  Configurable<std::vector<float>> cfVertexY{"cfVertexY", {-10., 10.}, "vertex y position range: {min, max}[cm], with convention: min <= Vy < max"};
  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<float> cfMinVertexDistanceFromIP{"cfMinVertexDistanceFromIP", {0.000001}, "if sqrt(vx^2+vy^2+vz^2) < cfMinVertexDistanceFromIP [cm], the event is reject. IP = nominal Interaction Point."};
  Configurable<std::vector<int>> cfNContributors{"cfNContributors", {2, 1000000000}, "Number of vertex contributors: {min, max}, with convention: min <= N < max"};
  Configurable<std::vector<float>> cfImpactParameter{"cfImpactParameter", {-1, 1000000000}, "Impact parameter range (can be used only for sim): {min, max}, with convention: min <= IP < max"};
  Configurable<std::vector<float>> cfEventPlaneAngle{"cfEventPlaneAngle", {-o2::constants::math::PI, o2::constants::math::TwoPI}, "Event Plane Angle range (can be used only for sim): {min, max}, with convention: min <= EP < max"};
  Configurable<std::vector<float>> cfOccupancy{"cfOccupancy", {-0.0001, 1000000000}, "Range for occupancy (use cfOccupancyEstimator to set specific estimator): {min, max}, with convention: min <= X < max. Important: remember that 0. has to be included, therefore I set -0.0001 by default for low edge"};
  Configurable<std::vector<float>> cfInteractionRate{"cfInteractionRate", {0.1, 1000000000.}, "Range for interaction rate (in kHz): {min, max}, with convention: min <= X < max"};
  Configurable<std::vector<float>> cfCurrentRunDuration{"cfCurrentRunDuration", {-2, 1000000000}, "Range for current run duration (i.e. seconds since start of run) in seconds: {min, max}, with convention: min <= X < max. Only collisions taken in this range (measured from SOR) are taken for analysis"};
  Configurable<std::vector<float>> cfMultMCNParticlesEta08{"cfMultMCNParticlesEta08", {-1, 1000000000}, "Range for MultMCNParticlesEta08 : {min, max}, with convention: min <= X < max"};
  Configurable<std::string> cfTrigger{"cfTrigger", "some supported trigger (e.g. \"kINT7\" for Run 2, \"kTVXinTRD\" for Run 3, etc...)", "set here some supported trigger"};
  Configurable<bool> cfUseSel7{"cfUseSel7", false, "use for Run 1 and 2 data and MC (see official doc)"};
  Configurable<bool> cfUseSel8{"cfUseSel8", false, "use for Run 3 data and MC (see official doc)"};
  Configurable<std::string> cfMultiplicityEstimator{"cfMultiplicityEstimator", "SelectedTracks", "all results vs. mult are calculated against this multiplicity. Can be set to SelectedTracks (calculated internally), ReferenceMultiplicity (calculated outside of my code), etc."};
  Configurable<std::string> cfReferenceMultiplicityEstimator{"cfReferenceMultiplicityEstimator", "MultFT0C", "Reference multiplicity, calculated outside of my code. Can be MultFT0C, MultFV0M, MultTPC, etc."};
  //  Configurable<std::string> cfCentralityEstimator{"cfCentralityEstimator", "some supported centrality estimator (e.g. CentFT0C, ...)", "set here some supported centrality estimator (CentFT0C, CentFT0M, CentFV0A, CentNTPV, ... for Run 3, and CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2 and 1) "};
  Configurable<std::string> cfCentralityEstimator{"cfCentralityEstimator", "CentFT0C", "set here some supported centrality estimator (CentFT0C, CentFT0M, CentFV0A, CentNTPV, ... for Run 3, and CentRun2V0M, CentRun2SPDTracklets, ..., for Run 2 and 1) "};
  Configurable<std::vector<int>> cfSelectedEvents{"cfSelectedEvents", {-1, 1000000000}, "Selected number of events to process (i.e. only events which survive event cuts): {min, max}, with convention: min <= N < max"};
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
  Configurable<bool> cfUseNoPileupTPC{"cfUseNoPileupTPC", false, "TBI 20250318 explanation"};
  Configurable<bool> cfUseNoPileupFromSPD{"cfUseNoPileupFromSPD", false, "TBI 20250318 explanation"};
  Configurable<bool> cfUseNoSPDOnVsOfPileup{"cfUseNoSPDOnVsOfPileup", false, "TBI 20250318 explanation"};
  Configurable<std::string> cfOccupancyEstimator{"cfOccupancyEstimator", "FT0COccupancyInTimeRange", "set here some supported occupancy estimator (TrackOccupancyInTimeRange, FT0COccupancyInTimeRange, ..."};
  Configurable<std::string> cfRefMultVsNContrUp{"cfRefMultVsNContrUp", "1200. + 0.20*x", "set here some formula in the mandatory format \"p0 + p1*x\" for the upper boundary cut in RefMult_vs_NContr correlation"};
  Configurable<std::string> cfRefMultVsNContrLow{"cfRefMultVsNContrLow", "-650. + 0.08*x", "set here some in the mandatory format \"p0 + p1*x\" for the lower boundary cut in RefMult_vs_NContr correlation"};
  Configurable<std::string> cfCentralityCorrelationsCut{"cfCentralityCorrelationsCut", "CentFT0C_CentFT0M", "Indicate two centrality estimators for the calculation of centrality correlation cut"};
  Configurable<float> cfCentralityCorrelationsCutTreshold{"cfCentralityCorrelationsCutTreshold", 10.0, "set the treshold for centrality correlation cut"};
  Configurable<std::string> cfCentralityCorrelationsCutVersion{"cfCentralityCorrelationsCutVersion", "Absolute", "set the version of centrality correlation cut. Supported: \"Relative\" and \"Absolute\""};
  Configurable<bool> cfUseFT0Bad{"cfUseFT0Bad", false, "TBI 20250516 explanation (or see enum)"};
  Configurable<bool> cfUseITSBad{"cfUseITSBad", false, "TBI 20250516 explanation (or see enum)"};
  Configurable<bool> cfUseITSLimAccMCRepr{"cfUseITSLimAccMCRepr", false, "TBI 20250516 explanation (or see enum)"};
  Configurable<bool> cfUseTPCBadTracking{"cfUseTPCBadTracking", false, "TBI 20250516 explanation (or see enum)"};
  Configurable<bool> cfUseTPCLimAccMCRepr{"cfUseTPCLimAccMCRepr", false, "TBI 20250516 explanation (or see enum)"};
  Configurable<bool> cfUseTPCBadPID{"cfUseTPCBadPID", false, "TBI 20250516 explanation (or see enum)"};
} cf_ec;

// *) Particle histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfFillParticleHistograms{"cfFillParticleHistograms", true, "if false, all 1D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms[...] is kTRUE, are filled"};
  Configurable<std::vector<std::string>> cfBookParticleHistograms{"cfBookParticleHistograms", {"1-Phi", "1-Pt", "1-Eta", "1-Charge", "1-tpcNClsFindable", "1-tpcNClsShared", "1-itsChi2NCl", "1-tpcNClsFound", "1-tpcNClsCrossedRows", "1-itsNCls", "1-itsNClsInnerBarrel", "1-tpcCrossedRowsOverFindableCls", "1-tpcFoundOverFindableCls", "1-tpcFractionSharedCls", "1-tpcChi2NCl", "1-dcaXY", "1-dcaZ", "0-PDG"}, "Book (1) or do not book (0) particle histogram"};
  Configurable<bool> cfFillParticleHistograms2D{"cfFillParticleHistograms2D", false, "if false, all 2D particle histograms are not filled. if kTRUE, the ones for which fBookParticleHistograms2D[...] is kTRUE, are filled"};
  Configurable<std::vector<std::string>> cfBookParticleHistograms2D{"cfBookParticleHistograms2D", {"1-Phi_vs_Pt", "1-Phi_vs_Eta"}, "Book (1) or do not book (0) 2D particle histograms"};
  Configurable<int> cfRebinSparse{"cfRebinSparse", 1, "used only for all fixed-length bins which are implemented directly for sparse histograms (i.e. not inherited from results histograms)"};
  Configurable<std::vector<std::string>> cfBookParticleSparseHistograms{"cfBookParticleSparseHistograms", {"0-DWPhi", "0-DWPt", "0-DWEta"}, "Book (1) or do not book (0) particular category of sparse histograms"};
  // TBI 20250223 add eventually configurable for FillParticleSparseHistogramsDimension
} cf_ph;

// *) Particle cuts:
struct : ConfigurableGroup {
  Configurable<std::vector<std::string>> cfUseParticleCuts{"cfUseParticleCuts", {"1-Phi", "1-Pt", "1-Eta", "1-Charge", "1-tpcNClsFindable", "1-tpcNClsShared", "1-itsChi2NCl", "1-tpcNClsFound", "1-tpcNClsCrossedRows", "1-itsNCls", "1-itsNClsInnerBarrel", "1-tpcCrossedRowsOverFindableCls", "1-tpcFoundOverFindableCls", "1-tpcFractionSharedCls", "1-tpcChi2NCl", "1-dcaXY", "1-dcaZ", "1-PDG", "0-trackCutFlag", "0-trackCutFlagFb1", "0-trackCutFlagFb2", "0-isQualityTrack", "1-isPrimaryTrack", "0-isInAcceptanceTrack", "0-isGlobalTrack", "1-isPVContributor", "0-PtDependentDCAxyParameterization"}, "Use (1) or do not use (0) particle cuts"};
  Configurable<bool> cfUseParticleCutCounterAbsolute{"cfUseParticleCutCounterAbsolute", false, "profile and save how many times each particle cut counter triggered (absolute). Use with care, as this is computationally heavy"};
  Configurable<bool> cfUseParticleCutCounterSequential{"cfUseParticleCutCounterSequential", false, "profile and save how many times each particle cut counter triggered (sequential). Use with care, as this is computationally heavy"};
  Configurable<std::vector<float>> cfPhi{"cfPhi", {0.0, o2::constants::math::TwoPI}, "phi range: {min, max}[rad], with convention: min <= phi < max"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "pt range: {min, max}[GeV], with convention: min <= pt < max"};
  Configurable<std::vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range: {min, max}, with convention: min <= eta < max"};
  Configurable<std::vector<float>> cfCharge{"cfCharge", {-1.5, 1.5}, "particle charge. {-1.5,0} = only negative, {0,1.5} = only positive"};
  Configurable<std::vector<float>> cftpcNClsFindable{"cftpcNClsFindable", {-1000., 1000.}, "tpcNClsFindable range: {min, max}, with convention: min <= cftpcNClsFindable < max"};
  Configurable<std::vector<float>> cftpcNClsShared{"cftpcNClsShared", {-1000., 1000.}, "tpcNClsShared range: {min, max}, with convention: min <= cftpcNClsShared < max"};
  Configurable<std::vector<float>> cfitsChi2NCl{"cfitsChi2NCl", {-1000., 36.}, "itsChi2NCl range: {min, max}, with convention: min <= cfitsChi2NCl < max"};
  Configurable<std::vector<float>> cftpcNClsFound{"cftpcNClsFound", {70., 1000.}, "tpcNClsFound range: {min, max}, with convention: min <= cftpcNClsFound < max"};
  Configurable<std::vector<float>> cftpcNClsCrossedRows{"cftpcNClsCrossedRows", {70., 1000.}, "tpcNClsCrossedRows range: {min, max}, with convention: min <= tpcNClsCrossedRows < max"};
  Configurable<std::vector<float>> cfitsNCls{"cfitsNCls", {5., 1000.}, "itsNCls range: {min, max}, with convention: min <= itsNCls < max"};
  Configurable<std::vector<float>> cfitsNClsInnerBarrel{"cfitsNClsInnerBarrel", {-1000., 1000.}, "itsNClsInnerBarrel range: {min, max}, with convention: min <= cfitsNClsInnerBarrel < max"};
  Configurable<std::vector<float>> cftpcCrossedRowsOverFindableCls{"cftpcCrossedRowsOverFindableCls", {0.8, 1000.}, "tpcCrossedRowsOverFindableCls range: {min, max}, with convention: min <= cftpcCrossedRowsOverFindableCls < max"};
  Configurable<std::vector<float>> cftpcFoundOverFindableCls{"cftpcFoundOverFindableCls", {0.8, 1000.}, "tpcFoundOverFindableCls range: {min, max}, with convention: min <= cftpcFoundOverFindableCls < max"};
  Configurable<std::vector<float>> cftpcFractionSharedCls{"cftpcFractionSharedCls", {-1000., 0.4}, "tpcFractionSharedCls range: {min, max}, with convention: min <= cftpcFractionSharedCls < max"};
  Configurable<std::vector<float>> cftpcChi2NCl{"cftpcChi2NCl", {-1000., 4.}, "tpcChi2NCl range: {min, max}, with convention: min <= cftpcChi2NCl < max"};
  Configurable<std::vector<float>> cfdcaXY{"cfdcaXY", {-2.4, 2.4}, "dcaXY range: {min, max}, with convention: min <= dcaXY < max (yes, DCA can be negative!)"};
  Configurable<std::vector<float>> cfdcaZ{"cfdcaZ", {-3.2, 3.2}, "dcaZ range: {min, max}, with convention: min <= dcaZ < max (yes, DCA can be negative!)"};
  Configurable<std::vector<float>> cfPDG{"cfPDG", {-5000., 5000.}, "PDG code"};
  Configurable<bool> cftrackCutFlag{"cftrackCutFlag", false, "general selection, particle cuts are tuned centrally. Use only in Run 3."};
  Configurable<bool> cftrackCutFlagFb1{"cftrackCutFlagFb1", false, "general selection + 1 point in ITS IB. Use only in Run 3."};
  Configurable<bool> cftrackCutFlagFb2{"cftrackCutFlagFb2", false, "general selection + 2 point in ITS IB. Use only in Run 3."};
  Configurable<bool> cfisQualityTrack{"cfisQualityTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisPrimaryTrack{"cfisPrimaryTrack", true, "Set to true by default both in Run 3 and Run 2 TBI 20250319 validate still for Run 1"};
  Configurable<bool> cfisInAcceptanceTrack{"cfisInAcceptanceTrack", false, "TBI 20250113 obsolete - see enum, to be removed"};
  Configurable<bool> cfisGlobalTrack{"cfisGlobalTrack", false, "TBI 20240510 add description"};
  Configurable<bool> cfisPVContributor{"cfisPVContributor", false, "Has this track contributed to the collision vertex fit"};
  Configurable<std::string> cfPtDependentDCAxyParameterization{"cfPtDependentDCAxyParameterization", "some formula TBI add some default formula, e.g. 0.0105+0.0350/x^1.1", "set here formula for pt-dependence DCAxy cut, in the following example format 0.0105+0.0350/x^1.1"};
  // TBI 20240426 do I need to add separate support for booleans to use each specific cut?
} cf_pc;

// *) Q-vector:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateQvectors{"cfCalculateQvectors", true, "calculate or not Q-vectors (all, also diff. ones). If I want only to fill control histograms, then set here false"};
} cf_qv;

// *) Multiparticle correlations:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateCorrelations{"cfCalculateCorrelations", false, "calculate or not multiparticle correlations"};
  Configurable<std::vector<std::string>> cfCalculateCorrelationsAsFunctionOf{"cfCalculateCorrelationsAsFunctionOf", {"1-Integrated", "1-Multiplicity", "1-Centrality", "0-Pt", "0-Eta", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "1-Vz", "0-Charge"}, "calculate or not correlations as a function of specified variable"};
} cf_mupa;

// *) Test0:
struct : ConfigurableGroup {

  // 1D:
  Configurable<bool> cfCalculateTest0{"cfCalculateTest0", false, "calculate or not Test0"};
  Configurable<std::vector<std::string>> cfCalculateTest0AsFunctionOf{"cfCalculateTest0AsFunctionOf", {"1-Integrated", "1-Multiplicity", "1-Centrality", "1-Pt", "1-Eta", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "1-Vz", "1-Charge"}, "calculate or not correlations as a function of specified variable"};

  // 2D:
  Configurable<bool> cfCalculate2DTest0{"cfCalculate2DTest0", false, "calculate or not 2D Test0 using TProfile2D"};
  Configurable<std::vector<std::string>> cfCalculate2DTest0AsFunctionOf{"cfCalculate2DTest0AsFunctionOf", {"1-Centrality_Pt", "1-Centrality_Eta", "1-Centrality_Charge", "1-Centrality_Vz", "1-Pt_Eta", "1-Pt_Charge", "1-Eta_Charge"}, "calculate or not correlations in 2D as a function of two specified variables."};

  // 3D:
  Configurable<bool> cfCalculate3DTest0{"cfCalculate3DTest0", false, "calculate or not 3D Test0 using TProfile3D"};
  Configurable<std::vector<std::string>> cfCalculate3DTest0AsFunctionOf{"cfCalculate3DTest0AsFunctionOf", {"1-Centrality_Pt_Eta", "1-Centrality_Pt_Charge", "1-Centrality_Pt_Vz", "1-Centrality_Eta_Vz", "1-Centrality_Eta_Charge", "1-Centrality_Vz_Charge", "1-Pt_Eta_Charge"}, "calculate or not correlations in 3D as a function of three specified variables."};

  Configurable<std::string> cfFileWithLabels{"cfFileWithLabels", "/home/abilandz/DatasetsO2/labels.root", "path to external ROOT file which specifies all labels"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
  Configurable<bool> cfUseDefaultLabels{"cfUseDefaultLabels", false, "use default internally hardwired labels, only for testing purposes"};
  Configurable<std::string> cfWhichDefaultLabels{"cfWhichDefaultLabels", "standard", "only for testing purposes, select one set of default labels, see GetDefaultObjArrayWithLabels for supported options"};
} cf_t0;

// *) Eta separation:
struct : ConfigurableGroup {
  Configurable<bool> cfCalculateEtaSeparations{"cfCalculateEtaSeparations", false, "calculate or not 2p corr. vs. eta separations"};
  Configurable<std::vector<std::string>> cfCalculateEtaSeparationsAsFunctionOf{"cfCalculateEtaSeparationsAsFunctionOf", {"1-Integrated", "1-Multiplicity", "1-Centrality", "1-Pt", "0-Eta", "1-Occupancy", "1-InteractionRate", "1-CurrentRunDuration", "1-Vz", "1-Charge"}, "calculate or not correlations as a function of specified variable"};
  Configurable<std::vector<float>> cfEtaSeparationsValues{"cfEtaSeparationsValues", {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}, "Eta separation between interval A (-eta) and B (+eta)"};
  Configurable<std::vector<std::string>> cfEtaSeparationsSkipHarmonics{"cfEtaSeparationsSkipHarmonics", {"0-v1", "0-v2", "0-v3", "0-v4", "1-v5", "1-v6", "1-v7", "1-v8", "1-v9"}, "For calculation of 2p correlation with eta separation these harmonics will be skipped (if first flag = \"0-v1\", v1 will be NOT be skipped in the calculus of 2p correlations with eta separations, etc.)"};
} cf_es;

// *) Particle weights:
struct : ConfigurableGroup {
  Configurable<bool> cfUsePhiWeights{"cfUsePhiWeights", false, "use or not phi weights"};
  Configurable<bool> cfUsePtWeights{"cfUsePtWeights", false, "use or not pt weights"};
  Configurable<bool> cfUseEtaWeights{"cfUseEtaWeights", false, "use or not eta weights"};
  Configurable<bool> cfUseDiffPhiPtWeights{"cfUseDiffPhiPtWeights", false, "use or not differential phi(pt) weights"};
  Configurable<bool> cfUseDiffPhiEtaWeights{"cfUseDiffPhiEtaWeights", false, "use or not differential phi(eta) weights"};
  Configurable<std::vector<std::string>> cfWhichDiffPhiWeights{"cfWhichDiffPhiWeights", {"1-wPhi", "1-wPt", "1-wEta", "1-wCharge", "1-wCentrality", "1-wVertexZ"}, "use (1) or do not use (0) differential phi weight for particular dimension. If only phi is set to 1, integrated phi weights are used. If phi is set to 0, ALL dimensions are switched off (yes!)"};
  Configurable<std::vector<std::string>> cfWhichDiffPtWeights{"cfWhichDiffPtWeights", {"0-wPt", "0-wCharge", "0-wCentrality"}, "use (1) or do not use (0) differential pt weight for particular dimension. If only pt is set to 1, integrated pt weights are used. If pt is set to 0, ALL dimensions are switched off (yes!)"};
  Configurable<std::vector<std::string>> cfWhichDiffEtaWeights{"cfWhichDiffEtaWeights", {"0-wEta", "0-wCharge", "0-wCentrality"}, "use (1) or do not use (0) differential eta weight for particular dimension. If only eta is set to 1, integrated eta weights are used. If eta is set to 0, ALL dimensions are switched off (yes!)"};
  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/home/abilandz/DatasetsO2/weights.root", "path to external ROOT file which holds all particle weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_pw;

// *) Centrality weights:
struct : ConfigurableGroup {
  Configurable<bool> cfUseCentralityWeights{"cfUseCentralityWeights", false, "use or not centrality weights"};
  Configurable<std::string> cfFileWithCentralityWeights{"cfFileWithCentralityWeights", "/home/abilandz/DatasetsO2/centralityWeights.root", "path to external ROOT file which holds centrality weights in O2 format"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
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
  Configurable<std::vector<int>> cfApplyNUAPDF{"cfApplyNUAPDF", {0, 0, 0}, "Apply (1) or do not apply (0) NUA on variable, ordering is the same as in enum eNUAPDF (phi, pt, eta)"};
  Configurable<std::vector<int>> cfUseDefaultNUAPDF{"cfUseDefaultNUAPDF", {1, 1, 1}, "Use (1) or do not use (0) default NUA profile, ordering is the same as in enum eNUAPDF (phi, pt, eta)"};
  Configurable<std::vector<std::string>> cfCustomNUAPDFHistNames{"cfCustomNUAPDFHistNames", {"a", "bb", "ccc"}, "the names of histograms holding custom NUA in an external file."};
  Configurable<std::string> cfFileWithCustomNUA{"cfFileWithCustomNUA", "/home/abilandz/DatasetsO2/customNUA.root", "path to external ROOT file which holds all histograms with custom NUA"}; // for AliEn file prepend "/alice/cern.ch/", for CCDB prepend "/alice-ccdb.cern.ch"
} cf_nua;

// *) Internal validation:
struct : ConfigurableGroup {
  Configurable<bool> cfUseInternalValidation{"cfUseInternalValidation", false, "perform internal validation using flow analysis on-the-fly"};
  Configurable<bool> cfInternalValidationForceBailout{"cfInternalValidationForceBailout", false, "force bailout (use only locally, since there is no graceful exit (yet))"};
  Configurable<unsigned int> cfnEventsInternalValidation{"cfnEventsInternalValidation", 0, "number of events simulated on-the-fly for internal validation"};
  Configurable<std::string> cfHarmonicsOptionInternalValidation{"cfHarmonicsOptionInternalValidation", "constant", "for internal validation, supported options are \"constant\", \"correlated\" and \"persistent\""};
  Configurable<bool> cfRescaleWithTheoreticalInput{"cfRescaleWithTheoreticalInput", false, "if kTRUE, all correlators are rescaled with theoretical input, so that all results in profiles are 1"};
  Configurable<bool> cfRandomizeReactionPlane{"cfRandomizeReactionPlane", true, "set to false only when validating against theoretical value the non-isotropic correlators"};
  Configurable<std::vector<float>> cfInternalValidationAmplitudes{"cfInternalValidationAmplitudes", {0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09}, "{v1, v2, v3, v4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of vn's is gMaxHarmonic."};
  Configurable<std::vector<float>> cfInternalValidationPlanes{"cfInternalValidationPlanes", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, "{Psi1, Psi2, Psi3, Psi4, ...} + has an effect only in combination with cfHarmonicsOptionInternalValidation = \"constant\". Max number of Psin's is gMaxHarmonic."};
  Configurable<std::vector<int>> cfMultRangeInternalValidation{"cfMultRangeInternalValidation", {1000, 1001}, "{min, max}, with convention: min <= M < max"};
} cf_iv;

// *) Results histograms:
struct : ConfigurableGroup {
  Configurable<bool> cfSaveResultsHistograms{"cfSaveResultsHistograms", false, "save or not results histograms"};

  // Fixed-length binning (default):
  Configurable<std::vector<double>> cfFixedLengthMultBins{"cfFixedLengthMultBins", {2000, 0., 20000.}, "nMultBins, multMin, multMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthCentBins{"cfFixedLengthCentBins", {110, 0., 110.}, "nCentBins, centMin, centMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthPtBins{"cfFixedLengthPtBins", {1000, 0., 10.}, "nPtBins, ptMin, ptMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthEtaBins{"cfFixedLengthEtaBins", {80, -2., 2.}, "nEtaBins, etaMin, etaMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthOccuBins{"cfFixedLengthOccuBins", {200, 0., 60000.}, "nOccuBins, occuMin, occuMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthIRBins{"cfFixedLengthIRBins", {1000, 0., 100.}, "nirBins, irMin, irMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthCRDBins{"cfFixedLengthCRDBins", {100000, 0., 100000.}, "ncrdBins, crdMin, crdMax (only for results histograms)"};
  Configurable<std::vector<double>> cfFixedLengthVzBins{"cfFixedLengthVzBins", {400, -20., 20.}, "nvzBins, vzMin, vzMax (only for results histograms)"};

  // Variable-length binning (per request):
  Configurable<bool> cfUseVariableLengthMultBins{"cfUseVariableLengthMultBins", false, "use or not variable-length multiplicity bins"};
  Configurable<std::vector<double>> cfVariableLengthMultBins{"cfVariableLengthMultBins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length multiplicity bins"};
  Configurable<bool> cfUseVariableLengthCentBins{"cfUseVariableLengthCentBins", false, "use or not variable-length centrality bins"};
  Configurable<std::vector<double>> cfVariableLengthCentBins{"cfVariableLengthCentBins", {0., 10., 50., 100.}, "variable-length centrality bins"};
  Configurable<bool> cfUseVariableLengthPtBins{"cfUseVariableLengthPtBins", true, "use or not variable-length pt bins"};
  Configurable<std::vector<double>> cfVariableLengthPtBins{"cfVariableLengthPtBins", {0.20, 0.25, 0.30, 0.35, 0.40, 0.50, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.50, 3.00, 4.00, 5.00}, "variable-length pt bins"};
  Configurable<bool> cfUseVariableLengthEtaBins{"cfUseVariableLengthEtaBins", true, "use or not variable-length eta bins"};
  Configurable<std::vector<double>> cfVariableLengthEtaBins{"cfVariableLengthEtaBins", {-0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8}, "variable-length eta bins"};
  Configurable<bool> cfUseVariableLengthOccuBins{"cfUseVariableLengthOccuBins", false, "use or not variable-length occupancy bins"};
  Configurable<std::vector<double>> cfVariableLengthOccuBins{"cfVariableLengthOccuBins", {0., 5., 6., 7., 8., 9., 100., 200., 500., 1000., 10000.}, "variable-length occupancy bins"};
  Configurable<bool> cfUseVariableLengthIRBins{"cfUseVariableLengthIRBins", false, "use or not variable-length interaction rate bins"};
  Configurable<std::vector<double>> cfVariableLengthIRBins{"cfVariableLengthIRBins", {0., 5., 10., 50., 100., 200.}, "variable-length ineraction rate bins"};
  Configurable<bool> cfUseVariableLengthCRDBins{"cfUseVariableLengthCRDBins", false, "use or not variable-length current run duration bins"};
  Configurable<std::vector<double>> cfVariableLengthCRDBins{"cfVariableLengthCRDBins", {0., 5., 10., 50., 100., 500.}, "variable-length current run duration bins"};
  Configurable<bool> cfUseVariableLengthVzBins{"cfUseVariableLengthVzBins", false, "use or not variable-length vertex z bins"};
  Configurable<std::vector<double>> cfVariableLengthVzBins{"cfVariableLengthVzBins", {-10., -8., -6., -4, -2., -1., 0., 1., 2., 4., 6., 8., 10.}, "variable-length vertex z bins"};

} cf_res;

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_CONFIGURABLES_H_
