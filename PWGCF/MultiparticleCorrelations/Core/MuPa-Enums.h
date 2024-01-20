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
  eTaskName = 1, // here I start from 1 exceptionally, because this enum's are used as bin contents, and ROOT starts counting bins from 1
  eRunNumber,
  eVerbose,
  eVerboseForEachParticle,
  eDoAdditionalInsanityChecks,
  eUseCCDB,
  eProcessRemainingEvents,
  eWhatToProcess,
  eRandomSeed,
  eUseFisherYates,
  eFixedNumberOfRandomlySelectedTracks,
  eConfiguration_N
};

enum eRecSim { eRec = 0,
               eSim = 1,
               eRecAndSim = 2 }; // TBI 20231021 find a better name?

enum eBeforeAfter { eBefore = 0,
                    eAfter = 1 };

enum eMinMax { eMin = 0,
               eMax = 1 };

enum eDefaultColors { eColor = kBlack,
                      eFillColor = kGray };

enum eWeights { wPHI = 0,
                wPT = 1,
                wETA = 2,
                eWeights_N };

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
  eNContributors, // Number of tracks used for the vertex
  eImpactParameter,
  eEventHistograms_N
};

enum eParticleHistograms {
  ePhi = 0,
  ePt,
  eEta,
  etpcNClsCrossedRows, // from aod::TracksExtra
  eDCA_xy,             // from aod::TracksDCA
  eDCA_z,
  eParticleHistograms_N
};

enum eAsFunctionOf {
  AFO_INTEGRATED = 0,
  AFO_MULTIPLICITY = 1,
  AFO_CENTRALITY = 2,
  AFO_PT = 3,
  AFO_ETA = 4,
  eAsFunctionOf_N
}; // prefix is needed, to avoid conflict with enum eKinematics

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
