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
//
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Header file for definition of constants in the DQ framework
/// \file DQConstants.h

#ifndef PWGDQ_CORE_DQCONSTANTS_H_
#define PWGDQ_CORE_DQCONSTANTS_H_

#include <cstdint>

namespace o2::dqconstants
{
// Numerical conventions
inline constexpr float kInvalidValue = -9999.f;
inline constexpr int kInvalidNum = -9999;
inline constexpr int kInvalidIndex = -1;

inline constexpr std::uint8_t  kFlag8  = 0u;
inline constexpr std::uint32_t kFlag32 = 0u;
inline constexpr std::uint64_t kFlag64 = 0ull;

inline constexpr float kInvalidFloat = -1.0f;
inline constexpr float kEmptyVariable = 0.0f;

inline constexpr double kMillisecondsPerMinute = 60000.0;

inline constexpr int kTwoProng = 2;
inline constexpr int kThreeProng = 3;
inline constexpr int kFourProng = 4;

// Checks
inline constexpr auto isCoordinateValid = [](float value) constexpr -> bool {
  return value >= 0.f;
};

inline constexpr auto isDCAInvalid = [](float value) constexpr -> bool {
  return value > 999;
};

inline constexpr auto isProcCodeInvalid = [](int procCode) constexpr -> bool {
  return procCode == 0;
};

inline constexpr auto checkCalibrationTypeInvalid = [](int type) constexpr -> bool {
  return type < 0 || type > 2;
};

inline constexpr auto getNumDegreesOfFreedom = [](int Ntrk) constexpr -> int {
  return 2*Ntrk - 3;
};

// Shared histogram class names
namespace histograms
{
// Event histograms names
inline constexpr char kEventBeforeCuts[] = "Event_BeforeCuts";
inline constexpr char kEventAfterCuts[] = "Event_AfterCuts";
inline constexpr char kEventMCTruth[] = "Event_MCTruth";

// Track histogram names
inline constexpr char kTrackBarrelBeforeCuts[] = "TrackBarrel_BeforeCuts";
inline constexpr char kTrackBarrelAfterCuts[] = "TrackBarrel_";

inline constexpr char kPairsBarrelSEPM[] = "PairsBarrelSEPM";
inline constexpr char kPairsBarrelSEPP[] = "PairsBarrelSEPP";
inline constexpr char kPairsBarrelSEMM[] = "PairsBarrelSEMM";

inline constexpr char kPairsBarrelMEPM[] = "PairsBarrelMEPM";
inline constexpr char kPairsBarrelMEPP[] = "PairsBarrelMEPP";
inline constexpr char kPairsBarrelMEMM[] = "PairsBarrelMEMM";

// Muon histogram names
inline constexpr char kMuonBeforeCuts[] = "Muons_BeforeCuts";
inline constexpr char kTrackMuonBeforeCuts[] = "TrackMuon_BeforeCuts";
inline constexpr char kTrackMuonAfterCuts[] = "TrackMuon_";

inline constexpr char kPairsMuonSEPM[] = "PairsMuonSEPM";
inline constexpr char kPairsMuonSEPP[] = "PairsMuonSEPP";
inline constexpr char kPairsMuonSEMM[] = "PairsMuonSEMM";

inline constexpr char kPairsMuonMEPM[] = "PairsMuonMEPM";
inline constexpr char kPairsMuonMEPP[] = "PairsMuonMEPP";
inline constexpr char kPairsMuonMEMM[] = "PairsMuonMEMM";
} // namespace histograms

namespace devices {
inline constexpr char kEventSelection[] = "analysis-event-selection";
inline constexpr char kTrackSelection[] = "analysis-track-selection";
inline constexpr char kMuonSelection[] = "analysis-muon-selection";
inline constexpr char kPrefilterSelection[] = "analysis-prefilter-selection";
inline constexpr char kSameEventPairing[] = "analysis-same-event-pairing";
inline constexpr char kAsymmetricPairing[] = "analysis-asymmetric-pairing";
inline constexpr char kDileptonTrack[] = "analysis-dilepton-track";
inline constexpr char kDileptonTrackTrack[] = "analysis-dilepton-track-track";
} // namespace devices

// Defaults
namespace defaults
{
// Event selection defaults
inline constexpr int kEventTagStartBit = 56;
inline constexpr float kMaxAbsVertexZCm = 10.f;
inline constexpr char kDefaultRCTLabel[] = "CBT";

inline constexpr int kNumDalitzSelections = 8;

inline constexpr int kNumTPCPostCalibrationTypes = 4;

// Track defaults
inline constexpr int kNumDCAPars = 2;
inline constexpr int kIndexDCAxy = 0;
inline constexpr int kIndexDCAz = 0;
inline constexpr double kInitialDCA = 1e10;

// Bimodality defaults
inline constexpr float kBimodalityInvalidBinWidth = -1.0;
inline constexpr float kBimodalityBinWidth = 0.005;

// KFParticle defaults
inline constexpr int kKFParticleMassConstrainedConstruction = 2;
inline constexpr int kKFParticleLegacyConstructionMethod = 3;

inline constexpr int kTrkCovEntries = 21;
inline constexpr int kNumTrkCovPars = 6;

inline constexpr int kTrkCovEntriesNoX = 15;
inline constexpr int kNumTrkCovParsNoX = 5;

inline constexpr int kNumTrkComponents = 3;

// Vertexing defaults
inline constexpr float kVertexingLxyzTolerance = 1.0e-8f;
inline constexpr float kInitialShift = 0.0f;

inline constexpr float kGJMinMagnitude = 1.e-6f;

// Detector defaults
inline constexpr int kPostCalibElectron = 0;
inline constexpr int kPostCalibPion = 1;
inline constexpr int kPostCalibKaon = 2;
inline constexpr int kPostCalibProton = 3;

inline constexpr float kMuonAbsorberExitZ = -505.; // cm

inline constexpr double kInvalidZnEnergy = -1.e12;
inline constexpr int kZDCNeutronNTowers = 4;

// Hadron defaults
inline constexpr float kInitialHadronMass = 0.0f;

// Lifetimes
inline constexpr float kLifetimeBPlus = 1.638e-12f; // s
inline constexpr float kLifetimeB0 = 1.517e-12f;    // s

} // namespace defaults

} // namespace o2::aod::dqconstants

#endif //PWGDQ_CORE_DQCONSTANTS_H_