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
  eVerboseForEachParticle,
  eDoAdditionalInsanityChecks,
  eInsanityCheckForEachParticle,
  eUseCCDB,
  eWhichProcess,
  eRandomSeed,
  eUseFisherYates,
  eFixedNumberOfRandomlySelectedTracks,
  eUseStopwatch,
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

enum eBeforeAfter { eBefore = 0,
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
  eNumberOfEvents = 0, // Total events = eNumberOfEvents + eBefore, Selected events = eNumberOfEvents + eAfter
  eTotalMultiplicity,
  eSelectedTracks,
  eMultFV0M,      // ref. mult from helper task o2-analysis-multiplicity-table
  eMultFT0M,      // ref. mult from helper task o2-analysis-multiplicity-table
  eMultTPC,       // ref. mult from helper task o2-analysis-multiplicity-table
  eMultNTracksPV, // ref. mult from helper task o2-analysis-multiplicity-table
  eCentrality,    // default centrality estimator
  eVertex_x,
  eVertex_y,
  eVertex_z,
  eNContributors, // number of tracks used for the vertex
  eImpactParameter,
  eEventHistograms_N
};

enum eEventHistograms2D {
  eVertex_z_vs_MultTPC = 0,
  eVertex_z_vs_NContributors,
  eEventHistograms2D_N
};

enum eEventCuts {
  eTrigger = eEventHistograms_N, // yes, because I do not want to duplicate the same enums from eEventHistograms here
  eSel7,                         // Event selection decision based on V0A & V0C => use only in Run 2 and Run 1
  eSel8,                         // Event selection decision based on TVX => use only in Run 3
  eCentralityEstimator,
  eSelectedEvents, // selected events = eNumberOfEvents + eAfter => therefore I do not need a special histogram for it
  eEventCuts_N
};

enum eParticleHistograms {

  // from o2::aod::Tracks: (Track parameters at collision vertex)
  ePhi = 0,
  ePt,
  eEta,
  eCharge, // Charge: positive: 1, negative: -1

  // from o2::aod::TracksExtra_001:
  etpcNClsFindable,
  etpcNClsShared,
  etpcNClsFound,
  etpcNClsCrossedRows,
  eitsNCls,
  eitsNClsInnerBarrel,
  etpcCrossedRowsOverFindableCls,
  etpcFoundOverFindableCls,
  etpcFractionSharedCls,

  // from o2::aod::TracksDCA:
  edcaXY,
  edcaZ,

  // the rest:
  ePDG,

  // counter:
  eParticleHistograms_N
};

enum eParticleHistograms2D {
  ePhiPt = 0,
  ePhiEta,
  eParticleHistograms2D_N
};

enum eParticleCuts {

  // from o2::aod::TrackSelection:
  etrackCutFlagFb1 = eParticleHistograms_N, // yes, because I do not want to duplicate the same enums from eParticleHistograms here
  etrackCutFlagFb2,
  eisQualityTrack,
  eisPrimaryTrack,
  eisInAcceptanceTrack,
  eisGlobalTrack,
  eParticleCuts_N
};

enum eAsFunctionOf {
  AFO_INTEGRATED = 0,
  AFO_MULTIPLICITY = 1, // vs. default multiplicity, which is (at the moment) fSelectedTracks, i.e. number of tracks in Q-vector
  AFO_CENTRALITY = 2,   // vs. default centrality estimator, see how it's calculated in DetermineCentrality(...)
  AFO_PT = 3,
  AFO_ETA = 4,
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

enum eEventCounter {
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

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
