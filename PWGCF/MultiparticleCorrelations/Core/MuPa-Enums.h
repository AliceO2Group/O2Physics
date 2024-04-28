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
  eNumberOfEvents = 0,
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

enum eEventCuts {
  eTrigger = eEventHistograms_N, // yes, because I do not want to duplicate the same enum's from eEventHistograms here
  eSel7,
  eSel8,
  eCentralityEstimator,
  eEventCuts_N
};

enum eParticleHistograms {
  ePhi = 0,
  ePt,
  eEta,
  etpcNClsCrossedRows, // from aod::TracksExtra
  eDCA_xy,             // from aod::TracksDCA
  eDCA_z,
  ePDG,
  eParticleHistograms_N
};

enum eParticleHistograms2D {
  ePhiPt = 0,
  ePhiEta,
  eParticleHistograms2D_N
};

enum eParticleCuts {
  eTBI = eParticleHistograms_N, // yes, because I do not want to duplicate the same enum's from eParticleHistograms here
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
#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
