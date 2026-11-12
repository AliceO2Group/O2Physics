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
inline constexpr float InvalidValueF = -9999.f;
inline constexpr int InvalidValueI = -9999;
inline constexpr int InvalidIndex = -1;

inline constexpr std::uint8_t Flag8 = 0u;
inline constexpr std::uint32_t Flag32 = 0u;
inline constexpr std::uint64_t Flag64 = 0ull;

inline constexpr float EmptyVariable = 0.0f;

inline constexpr double MillisecondsPerMinute = 60000.0;

inline constexpr int TwoProng = 2;
inline constexpr int ThreeProng = 3;
inline constexpr int FourProng = 4;

inline constexpr int TPCCalibTypeStart = 0;
inline constexpr int TPCCalibTypeEnd = 2;

inline constexpr float InvalidDCA = 998.0f;
inline constexpr float InvalidStdDev = -1.0;

// Checks
inline constexpr auto IsCoordinateInvalid = [](float value) constexpr -> bool {
  return value < 0.f;
};

inline constexpr auto IsDCAInvalid = [](float value) constexpr -> bool {
  return value > InvalidDCA;
};

inline constexpr auto IsStdDevValid = [](float value) constexpr -> bool {
  return value < 0.f;
};

inline constexpr auto IsProcCodeInvalid = [](int procCode) constexpr -> bool {
  return procCode == 0;
};

inline constexpr auto CheckCalibrationTypeInvalid = [](int type) constexpr -> bool {
  return type < TPCCalibTypeStart || type > TPCCalibTypeEnd;
};

inline constexpr auto GetNumDegreesOfFreedom = [](int Ntrk) constexpr -> int {
  return 2 * Ntrk - 3;
};

// Shared histogram class names
namespace histograms
{
// Event histograms names
inline constexpr char EventBeforeCuts[] = "Event_BeforeCuts";
inline constexpr char EventAfterCuts[] = "Event_AfterCuts";
inline constexpr char EventMCTruth[] = "Event_MCTruth";

// Track histogram names
inline constexpr char TrackBarrelBeforeCuts[] = "TrackBarrel_BeforeCuts";
inline constexpr char TrackBarrelAfterCuts[] = "TrackBarrel_";

inline constexpr char PairsBarrelSEPM[] = "PairsBarrelSEPM";
inline constexpr char PairsBarrelSEPP[] = "PairsBarrelSEPP";
inline constexpr char PairsBarrelSEMM[] = "PairsBarrelSEMM";

inline constexpr char PairsBarrelMEPM[] = "PairsBarrelMEPM";
inline constexpr char PairsBarrelMEPP[] = "PairsBarrelMEPP";
inline constexpr char PairsBarrelMEMM[] = "PairsBarrelMEMM";

// Muon histogram names
inline constexpr char MuonBeforeCuts[] = "Muons_BeforeCuts";
inline constexpr char TrackMuonBeforeCuts[] = "TrackMuon_BeforeCuts";
inline constexpr char TrackMuonAfterCuts[] = "TrackMuon_";

inline constexpr char PairsMuonSEPM[] = "PairsMuonSEPM";
inline constexpr char PairsMuonSEPP[] = "PairsMuonSEPP";
inline constexpr char PairsMuonSEMM[] = "PairsMuonSEMM";

inline constexpr char PairsMuonMEPM[] = "PairsMuonMEPM";
inline constexpr char PairsMuonMEPP[] = "PairsMuonMEPP";
inline constexpr char PairsMuonMEMM[] = "PairsMuonMEMM";
} // namespace histograms

namespace devices
{
inline constexpr char EventSelection[] = "analysis-event-selection";
inline constexpr char TrackSelection[] = "analysis-track-selection";
inline constexpr char MuonSelection[] = "analysis-muon-selection";
inline constexpr char PrefilterSelection[] = "analysis-prefilter-selection";
inline constexpr char SameEventPairing[] = "analysis-same-event-pairing";
inline constexpr char AsymmetricPairing[] = "analysis-asymmetric-pairing";
inline constexpr char DileptonTrack[] = "analysis-dilepton-track";
inline constexpr char DileptonTrackTrack[] = "analysis-dilepton-track-track";
} // namespace devices

// Defaults
namespace defaults
{
// Event selection defaults
inline constexpr int EventTagStartBit = 56;
inline constexpr float MaxAbsVertexZCm = 10.f;
inline constexpr char DefaultRCTLabel[] = "CBT";

inline constexpr int NumDalitzSelections = 8;

inline constexpr int NumTPCPostCalibrationTypes = 4;

// Track defaults
inline constexpr int NumDCAPars = 2;
inline constexpr int IndexDCAxy = 0;
inline constexpr int IndexDCAz = 1;
inline constexpr double InitialDCA = 1e10;

inline constexpr int RotationType1 = 1;
inline constexpr int RotationType2 = 2;
inline constexpr int RotationType3 = 3;

inline constexpr int BdtSingleClassOutput = 1;
inline constexpr int BdtMultiClassOutput = 3;

// Bimodality defaults
inline constexpr int BimodalMinDataSize = 3;
inline constexpr float BimodalityInvalidBinWidth = -1.0;
inline constexpr float BimodalityBinWidth = 0.005;

inline constexpr float BimodalityDCAFrac001 = 0.01;
inline constexpr float BimodalityDCAFrac002 = 0.02;
inline constexpr float BimodalityDCAFrac005 = 0.05;
inline constexpr float BimodalityDCAFrac01 = 0.1;
inline constexpr float BimodalityDCAFrac02 = 0.2;
inline constexpr float BimodalityDCAFrac05 = 0.5;
inline constexpr float BimodalityDCAFrac1 = 1.00;

// KFParticle defaults
inline constexpr int KFParticleMassConstrainedConstruction = 2;
inline constexpr int KFParticleLegacyConstructionMethod = 3;

inline constexpr int TrkCovEntries = 21;
inline constexpr int NumTrkCovPars = 6;

inline constexpr int TrkCovEntriesNoX = 15;
inline constexpr int NumTrkCovParsNoX = 5;

inline constexpr int NumTrkComponents = 3;

// Vertexing defaults
inline constexpr float InvalidVertexingError = 1.e8f;
inline constexpr float VertexingLxyzTolerance = 1.0e-8f;
inline constexpr float InitialShift = 0.0f;

inline constexpr float GJMinMagnitude = 1.e-6f;

// Detector defaults
inline constexpr int PostCalibElectron = 0;
inline constexpr int PostCalibPion = 1;
inline constexpr int PostCalibKaon = 2;
inline constexpr int PostCalibProton = 3;

inline constexpr int TPCPIDCalibrationMeanSigma = 1;
inline constexpr int TPCPIDCalibrationMeanSigmaStatus = 2;

inline constexpr int NumITSLayers = 7;

inline constexpr int NumMFTLayers = 10;
inline constexpr float MuonAbsorberExitZ = -505.; // cm

inline constexpr double InvalidZnEnergy = -1.e12;
inline constexpr int ZDCNeutronNTowers = 4;

inline constexpr float CurvatureNormalizationFactor = 2.466;

// Hadron defaults
inline constexpr float InitialHadronMass = 0.0f;

// Lifetimes
inline constexpr float LifetimeBPlus = 1.638e-12f; // s
inline constexpr float LifetimeB0 = 1.517e-12f;    // s

} // namespace defaults

} // namespace o2::dqconstants

#endif // PWGDQ_CORE_DQCONSTANTS_H_
