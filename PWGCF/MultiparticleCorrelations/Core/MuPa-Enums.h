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

#ifndef PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
#define PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_

enum eConfiguration {
  eTaskName = 1, // here I start from 1 exceptionally, because these enums are used as bin contents, and ROOT starts counting bins from 1
  eRunNumber,
  eDryRun,
  eVerbose,
  eVerboseUtility,
  eVerboseForEachParticle,
  eVerboseEventCounter,
  ePlainPrintout,
  eDoAdditionalInsanityChecks,
  eInsanityCheckForEachParticle,
  eWhichProcess,
  eRandomSeed,
  eUseFisherYates,
  eFixedNumberOfRandomlySelectedTracks,
  eUseStopwatch,
  eFloatingPointPrecision,
  eSequentialBailout,
  eUseSpecificCuts,
  eWhichSpecificCuts,
  eConfiguration_N
};

enum eProcess {
  eProcessRec = 0,     // Run 3, only reconstructed
  eProcessRecSim,      // Run 3, both reconstructed and simulated
  eProcessSim,         // Run 3, only simulated
  eProcessRec_Run2,    // Run 2, only reconstructed
  eProcessRecSim_Run2, // Run 2, both reconstructed and simulated
  eProcessSim_Run2,    // Run 2, only simulated
  eProcessRec_Run1,    // Run 1, only reconstructed
  eProcessRecSim_Run1, // Run 1, both reconstructed and simulated
  eProcessSim_Run1,    // Run 1, only simulated
  eProcessTest,        // minimum subscription to the tables, for testing purposes
  // Generic flags, calculated and set from individual flags above in DefaultConfiguration(), AFTER process switch was taken into account:
  eGenericRec,    // generic "Rec" case, eTest is treated for the time being as "Rec"
  eGenericRecSim, // generic "RecSim" case
  eGenericSim,    // generic "Sim" case
  eProcess_N
};

enum eRecSim { eRec = 0,
               eSim,
               eRecAndSim,
               eRec_Run2, // converted Run 2 data
               eSim_Run2,
               eRecAndSim_Run2,
               eRec_Run1, // converted Run 1 data
               eSim_Run1,
               eRecAndSim_Run1,
               eTest };

enum eBeforeAfter { eBefore = 0, // use this one for cuts
                    eAfter = 1 };

enum eMinMax { eMin = 0,
               eMax = 1 };

enum eXYZ { eX = 0,
            eY = 1,
            eZ = 2 };

enum eDefaultColors { eColor = kBlack,
                      eFillColor = kGray };

enum eWeights { wPHI = 0,
                wPT = 1,
                wETA = 2,
                eWeights_N };

enum eDiffWeights {
  wPHIPT = 0,
  wPHIETA,
  eDiffWeights_N
};

enum eVnPsin { eVn = 0,
               ePsin = 1 };

enum eEventHistograms {
  eNumberOfEvents = 0,    // Total events = eNumberOfEvents + eBefore, Selected events = eNumberOfEvents + eAfter
  eTotalMultiplicity,     // TBI 20241123 I define it as tracks.size(), but most likely this I do not need this
  eMultiplicity,          // see documentation below for ebye.fMultiplicity
  eReferenceMultiplicity, // see documentation below for ebye.fReferenceMultiplicity
  eCentrality,            // default centrality estimator
  eVertex_x,
  eVertex_y,
  eVertex_z,
  eNContributors, // number of tracks used for the vertex
  eImpactParameter,
  eOccupancy,             // from helper task o2-analysis-event-selection, see also IA's presentation in https://indico.cern.ch/event/1464946, slide 38. Use specific occupancy estimator via eOccupancyEstimator
  eMultMCNParticlesEta08, // from helper task table o2::aod::MultMCExtras
  eEventHistograms_N
};

enum eEventCuts {
  // a) For available event selection bits, check https://github.com/AliceO2Group/O2Physics/blob/master/Common/CCDB/EventSelectionParams.cxx
  // b) Some settings are configurable, check: https://github.com/AliceO2Group/O2Physics/blob/master/Common/TableProducer/eventSelection.cxx
  eTrigger = eEventHistograms_N,   // Do NOT use eTrigger for Run 3. Validated only for Run 2, and it has to be "kINT7" . TBI 20240522 investigate for Run 1
  eSel7,                           // See def. of sel7 in Ref. b) above. Event selection decision based on V0A & V0C => use only in Run 2 and Run 1. TBI 20240522 I stil need to validate this one over MC
  eSel8,                           // See def. of sel7 in Ref. b) above. Event selection decision based on TVX => use only in Run 3, both for data and MC
                                   // *) As of 20240410, kNoITSROFrameBorder (only in MC) and kNoTimeFrameBorder event selection cuts are part of Sel8
                                   //    See also email from EK from 2024041
  eMultiplicityEstimator,          // see documentation below for ebye.fMultiplicity
  eReferenceMultiplicityEstimator, // see documentation below for ebye.fReferenceMultiplicity
  eCentralityEstimator,            // the default centrality estimator, set via configurable. All supported centrality estimators, for QA, etc, are in enum eCentralityEstimators
  eSelectedEvents,                 // selected events = eNumberOfEvents + eAfter => therefore I do not need a special histogram for it
  eNoSameBunchPileup,              // reject collisions in case of pileup with another collision in the same foundBC (emails from IA on 20240404 and EK on 20240410)
  eIsGoodZvtxFT0vsPV,              // small difference between z-vertex from PV and from FT0 (emails from IA on 20240404 and EK on 20240410)
  eIsVertexITSTPC,                 // at least one ITS-TPC track (reject vertices built from ITS-only tracks) (emails from IA on 20240404 and EK on 20240410
  eIsVertexTOFmatched,             // at least one of vertex contributors is matched to TOF
  eIsVertexTRDmatched,             // at least one of vertex contributors is matched to TRD
  eOccupancyEstimator,             // the default Occupancy estimator, set via configurable. All supported centrality estimators, for QA, etc, are in enum eOccupancyEstimators
  eEventCuts_N
};

enum eParticleHistograms {

  // from o2::aod::Tracks (Track parameters at their point closest to the collision vertex)
  ePhi = 0,
  ePt,
  eEta,
  eCharge, // Charge: positive: 1, negative: -1

  // from o2::aod::TracksExtra_001
  etpcNClsFindable,
  etpcNClsShared,
  etpcNClsFound,
  etpcNClsCrossedRows,
  eitsNCls,
  eitsNClsInnerBarrel,
  etpcCrossedRowsOverFindableCls,
  etpcFoundOverFindableCls,
  etpcFractionSharedCls,

  // from o2::aod::TracksDCA
  edcaXY,
  edcaZ,

  // the rest:
  ePDG,

  // counter:
  eParticleHistograms_N
};

enum eParticleHistograms2D { // All 2D histograms are first implemented in eQAParticleHistograms2D, the ones which I need regularly, are then promoted to this category.
  ePhiPt = 0,
  ePhiEta,
  eParticleHistograms2D_N
};

enum eParticleCuts {

  // from o2::aod::TrackSelection
  etrackCutFlagFb1 = eParticleHistograms_N, // do not use in Run 2 and 1
  etrackCutFlagFb2,                         // do not use in Run 2 and 1
  eisQualityTrack,                          // not validated in Run 3, but it can be used in Run 2 and Run 1 (for the latter, it yields to large NUA)
  eisPrimaryTrack,
  eisInAcceptanceTrack, // TBI 20240516 check and document how acceptance window is defined
  eisGlobalTrack,       // not validated in Run 3, but it can be used in Run 2 and Run 1 (for the latter, it yields to real holes in NUA)

  // special treatment:
  ePtDependentDCAxyParameterization,

  eParticleCuts_N
};

enum eAsFunctionOf {
  AFO_INTEGRATED = 0,
  AFO_MULTIPLICITY = 1, // vs. default multiplicity, which is (at the moment) fSelectedTracks, i.e. number of tracks in Q-vector
  AFO_CENTRALITY = 2,   // vs. default centrality estimator, see how it's calculated in DetermineCentrality(...)
  AFO_PT = 3,
  AFO_ETA = 4,
  AFO_OCCUPANCY = 5, // vs. default "occupancy" variable which is (at the moment) "TrackOccupancyInTimeRange" (alternative is "FT0COccupancyInTimeRange")
  eAsFunctionOf_N
}; // prefix is needed, to avoid conflict with enum eKinematics

enum eNUAPDF {
  ePhiNUAPDF = 0,
  ePtNUAPDF,
  eEtaNUAPDF,
  eNUAPDF_N
};

enum eqvectorKine { // Here "kine" originally meant "kinematic", i.e. vs. pt or vs. eta, now it's general.
  PTq = 0,
  ETAq,
  eqvectorKine_N
};

enum eTimer {
  eGlobal = 0,
  eLocal,
  eTimer_N
};

enum eEventCounterForDryRun {
  eFill = 0,
  ePrint
};

enum eCutModus {
  eCut = 0,           // standard, i.e. no cut counters are used
  eCutCounterBinning, // dry call to EventCuts and ParticleCuts, just to establish order of binning in CutCountets, which resembles order of cut implementation
  eCutCounterAbsolute,
  eCutCounterSequential
};

enum eCutCounter {
  eAbsolute = 0,
  eSequential,
  eCutCounter_N
};

enum eQAEventHistograms2D {
  // General (estimators can be chosen via configurables):
  eMultiplicity_vs_ReferenceMultiplicity = 0,
  eMultiplicity_vs_NContributors,
  eMultiplicity_vs_Centrality,
  eMultiplicity_vs_Vertex_z,
  eMultiplicity_vs_Occupancy,
  eReferenceMultiplicity_vs_NContributors,
  eReferenceMultiplicity_vs_Centrality,
  eReferenceMultiplicity_vs_Vertex_z,
  eReferenceMultiplicity_vs_Occupancy,
  eNContributors_vs_Centrality,
  eNContributors_vs_Vertex_z,
  eNContributors_vs_Occupancy,
  eCentrality_vs_Vertex_z,
  eCentrality_vs_Occupancy,
  eVertex_z_vs_Occupancy,
  // ...
  // Specific (estimators are hardwired):
  eCentFT0C_vs_CentNTPV,                // Run 3 centrality
  eCentFT0M_vs_CentNTPV,                // Run 3 centrality
  eCentRun2V0M_vs_CentRun2SPDTracklets, // Run 2 centrality (do not use in Run 1 converted, because there is no centrality information)
  eTrackOccupancyInTimeRange_vs_FT0COccupancyInTimeRange,
  // ...
  eQAEventHistograms2D_N
};

enum eQAParticleHistograms2D {
  ePt_vs_dcaXY,
  eQAParticleHistograms2D_N
};

enum eReferenceMultiplicityEstimators {
  // Run 3:
  eMultTPC = 0,
  eMultFV0M,      // ref. mult from helper task o2-analysis-multiplicity-table
  eMultFT0C,      // ref. mult from helper task o2-analysis-multiplicity-table
  eMultFT0M,      // ref. mult from helper task o2-analysis-multiplicity-table
  eMultNTracksPV, // ref. mult from helper task o2-analysis-multiplicity-table
  // Run 2:
  eMultTracklets, // ref. mult from helper task o2-analysis-multiplicity-table, use only for Run 2
  eReferenceMultiplicityEstimators_N
};

enum eCentralityEstimators {
  // Run 3:
  eCentFT0C = 0,
  eCentFT0M,
  eCentFV0A,
  eCentNTPV,
  // Run 2:
  eCentRun2V0M,
  eCentRun2SPDTracklets,
  eCentralityEstimators_N
};

enum eOccupancyEstimators {
  eTrackOccupancyInTimeRange, // from helper task o2-analysis-event-selection, see also IA's presentation in https://indico.cern.ch/event/1464946, slide 38
  eFT0COccupancyInTimeRange,  // from helper task o2-analysis-event-selection
  eOccupancyEstimators_N
};

enum eEventCounter {
  eTotal,     // total number of events, before any cuts are applied
  eProcessed, // number of processed events, i.e. number of events which survived cuts and on which analysis have been performed
  eEventCounter_N
};

enum eSpecificCuts {
  eLHC23zzh,
  eSpecificCuts_N
};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
