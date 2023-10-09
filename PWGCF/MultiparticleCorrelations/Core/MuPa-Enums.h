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
  eTaskName = 1,
  eRunNumber = 2,
  eVerbose = 3,
  eVerboseForEachParticle = 4,
  eUseCCDB = 5,
  eConfiguration_N
};

enum eRecoSim { eRec = 0,
                eSim = 1 };

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
  eNumberOfEvents,
  eTotalMultiplicity,
  eSelectedTracks,
  eCentrality,
  eVertex_x,
  eVertex_y,
  eVertex_z,
  eNContributors, // Number of tracks used for the vertex
  eEventHistograms_N
};

enum eParticleHistograms {
  ePhi,
  ePt,
  eEta,
  etpcNClsCrossedRows, // from aod::TracksExtra
  eDCA_xy,             // from aod::TracksDCA
  eDCA_z,
  eParticleHistograms_N
};

#endif // PWGCF_MULTIPARTICLECORRELATIONS_CORE_MUPA_ENUMS_H_
