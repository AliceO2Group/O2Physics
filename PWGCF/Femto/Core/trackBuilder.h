// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file trackBuilder.h
/// \brief track builder
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_TRACKBUILDER_H_
#define PWGCF_FEMTO_CORE_TRACKBUILDER_H_

#include "PWGCF/Femto/Core/baseSelection.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/modes.h"
#include "PWGCF/Femto/Core/selectionContainer.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/Configurable.h"

#include "fairlogger/Logger.h"

#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femto
{
namespace trackbuilder
{

struct ConfTrackFilters : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackFilters");
  // kinematic cuts for filtering tracks
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
};

struct ConfTrackBits : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackBits");
  // track quality cuts
  o2::framework::Configurable<std::vector<float>> tpcClustersMin{"tpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};
  o2::framework::Configurable<std::vector<float>> tpcCrossedRowsMin{"tpcCrossedRowsMin", {80.f}, "Minimum number of crossed rows in TPC"};
  o2::framework::Configurable<std::vector<float>> tpcClustersOverCrossedRows{"tpcClustersOverCrossedRows", {0.83f}, "Minimum fraction of clusters over crossed rows in TPC"};
  o2::framework::Configurable<std::vector<float>> tpcSharedClustersMax{"tpcSharedClustersMax", {160.f}, "Maximum number of shared clusters in TPC"};
  o2::framework::Configurable<std::vector<float>> tpcSharedClusterFractionMax{"tpcSharedClusterFractionMax", {1.f}, "Maximum fraction of shared clusters in TPC"};
  o2::framework::Configurable<std::vector<float>> itsClustersMin{"itsClustersMin", {5.f}, "Minimum number of clusters in ITS"};
  o2::framework::Configurable<std::vector<float>> itsIbClustersMin{"itsIbClustersMin", {3.f}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
  o2::framework::Configurable<std::vector<std::string>> dcaxyMax{"dcaxyMax", {"0.004 + 0.013*pow(x, -1)"}, "Maximum |dca_xy| as a function of pT. Has to be a valid TForumal, where x=pt"};
  o2::framework::Configurable<std::vector<std::string>> dcazMax{"dcazMax", {"0.004 + 0.013*pow(x, -1)"}, "Maximum |dca_z| as a function of pT. Has to be a valid TForumal, where x=pt"};

  // Electron PID cuts
  o2::framework::Configurable<bool> requirePidElectron{"requirePidElectron", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofElectron{"minMomTofElectron", 0.3, "Minimum momentum to required TOF PID for Electron"};
  o2::framework::Configurable<std::vector<std::string>> itsElectron{"itsElectron", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Electron PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcElectron{"tpcElectron", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Electron PID"};
  o2::framework::Configurable<std::vector<std::string>> tofElectron{"tofElectron", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsElectron{"tpcitsElectron", {}, "Maximum nsigma_TPCITS for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofElectron{"tpctofElectron", {}, "Maximum nsigma_TPCTOF for Electron PID"};

  // Pion PID cuts
  o2::framework::Configurable<bool> requirePidPion{"requirePidPion", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofPion{"minMomTofPion", 0.5, "Minimum momentum to required TOF PID for Pion"};
  o2::framework::Configurable<std::vector<std::string>> itsPion{"itsPion", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Pion PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcPion{"tpcPion", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Pion PID"};
  o2::framework::Configurable<std::vector<std::string>> tofPion{"tofPion", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsPion{"tpcitsPion", {}, "Maximum nsigma_TPCITS for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tpctofPion{"tpctofPion", {}, "Maximum nsigma_TPCTOF for Pion PID"};

  // Kaon PID cuts
  o2::framework::Configurable<bool> requirePidKaon{"requirePidKaon", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofKaon{"minMomTofKaon", 0.4, "Minimum momentum to required TOF PID for Kaon"};
  o2::framework::Configurable<std::vector<std::string>> itsKaon{"itsKaon", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Kaon PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcKaon{"tpcKaon", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Kaon PID"};
  o2::framework::Configurable<std::vector<std::string>> tofKaon{"tofKaon", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsKaon{"tpcitsKaon", {}, "Maximum nsigma_TPCITS for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpctofKaon{"tpctofKaon", {}, "Maximum nsigma_TPCTOF for Kaon PID"};

  // Proton PID cuts
  o2::framework::Configurable<bool> requirePidProton{"requirePidProton", true, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofProton{"minMomTofProton", 0.75, "Minimum momentum to required TOF PID for Proton"};
  o2::framework::Configurable<std::vector<std::string>> itsProton{"itsProton", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Proton PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcProton{"tpcProton", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Proton PID"};
  o2::framework::Configurable<std::vector<std::string>> tofProton{"tofProton", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsProton{"tpcitsProton", {}, "Maximum nsigma_TPCITS for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tpctofProton{"tpctofProton", {}, "Maximum nsigma_TPCTOF for Proton PID"};

  // Deuteron PID cuts
  o2::framework::Configurable<bool> requirePidDeuteron{"requirePidDeuteron", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofDeuteron{"minMomTofDeuteron", 1.2, "Minimum momentum to required TOF PID for Deuteron"};
  o2::framework::Configurable<std::vector<std::string>> itsDeuteron{"itsDeuteron", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Deuteron PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcDeuteron{"tpcDeuteron", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Deuteron PID"};
  o2::framework::Configurable<std::vector<std::string>> tofDeuteron{"tofDeuteron", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsDeuteron{"tpcitsDeuteron", {}, "Maximum nsigma_TPCITS for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofDeuteron{"tpctofDeuteron", {}, "Maximum nsigma_TPCTOF for Deuteron PID"};

  // Triton PID cuts
  o2::framework::Configurable<bool> requirePidTriton{"requirePidTriton", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofTriton{"minMomTofTriton", 1.4, "Minimum momentum to required TOF PID for Triton"};
  o2::framework::Configurable<std::vector<std::string>> itsTriton{"itsTriton", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Triton PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcTriton{"tpcTriton", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Triton PID"};
  o2::framework::Configurable<std::vector<std::string>> tofTriton{"tofTriton", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsTriton{"tpcitsTriton", {}, "Maximum nsigma_TPCITS for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tpctofTriton{"tpctofTriton", {}, "Maximum nsigma_TPCTOF for Triton PID"};

  // Helium PID cuts
  o2::framework::Configurable<bool> requirePidHelium{"requirePidHelium", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofHelium{"minMomTofHelium", 1.6, "Minimum momentum to required TOF PID for Helium"};
  o2::framework::Configurable<std::vector<std::string>> itsHelium{"itsHelium", {}, "Ranges LowerLimit;UpperLimit for nsigma_ITS for Helium PID"};
  o2::framework::Configurable<std::vector<std::string>> tpcHelium{"tpcHelium", {}, "Ranges LowerLimit;UpperLimit for nsigma_TPC for Helium PID"};
  o2::framework::Configurable<std::vector<std::string>> tofHelium{"tofHelium", {}, "Ranges LowerLimit;UpperLimit for nsigma_TOF for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsHelium{"tpcitsHelium", {}, "Maximum nsigma_TPCITS for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tpctofHelium{"tpctofHelium", {}, "Maximum nsigma_TPCTOF for Helium PID"};
};

// define the template structure for TrackSelection
template <const char* Prefix>
struct ConfTrackSelection : public o2::framework::ConfigurableGroup {
  std::string prefix = Prefix; // Unique prefix based on the template argument
  // configuration parameters
  o2::framework::Configurable<int> pdgCodeAbs{"pdgCodeAbs", 2212, "Absolute value of PDG code. Set sign of charge to -1 for antiparticle."};
  o2::framework::Configurable<int> chargeAbs{"chargeAbs", 1, "Absolute value of charge (e.g. 1 for most tracks, 2 for He3). Set sign of charge to -1 for antiparticle"};
  o2::framework::Configurable<int> chargeSign{"chargeSign", 1, "Track charge sign: +1 for positive, -1 for negative, 0 for both"};
  // filters for kinematics
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT (GeV/c)"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT (GeV/c)"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // track selection masks
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskLowMomentum{"maskLowMomentum", 1ul, "Bitmask for selections below momentum threshold"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskHighMomentum{"maskHighMomentum", 2ul, "Bitmask for selections above momentum threshold"};
  // track rejection masks
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> rejectionMaskLowMomentum{"rejectionMaskLowMomentum", 0ul, "Bitmask for rejections below momentum threshold"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> rejectionMaskHighMomentum{"rejectionMaskHighMomentum", 0ul, "Bitmask for rejections above momentum threshold"};
  // momentum threshold for PID usage
  o2::framework::Configurable<float> pidThres{"pidThres", 1.2f, "Momentum threshold for using TPCTOF/TOF pid for tracks with large momentum (GeV/c)"};
};

// Define unique prefixes as constexpr string literals
constexpr const char PrefixTrackSelection1[] = "TrackSelection1";
constexpr const char PrefixTrackSelection2[] = "TrackSelection2";
constexpr const char PrefixTrackSelection3[] = "TrackSelection3";

// Instantiate different instances with unique prefixes
using ConfTrackSelection1 = ConfTrackSelection<PrefixTrackSelection1>;
using ConfTrackSelection2 = ConfTrackSelection<PrefixTrackSelection2>;
using ConfTrackSelection3 = ConfTrackSelection<PrefixTrackSelection3>;

/// enum for all track selections
enum TrackSels {
  // track quality cuts
  kTPCnClsMin,          ///< Min. number of TPC clusters
  kTPCcRowsMin,         ///< Min. number of crossed TPC rows
  kTPCnClsOvercRowsMin, ///< Min. fraction of TPC clusters of TPC crossed rows
  kTPCsClsMax,          ///< Max. number of shared TPC clusters
  kTPCsClsFracMax,      ///< Max. fractions of shared TPC clusters
  kITSnClsMin,          ///< Min. number of ITS clusters
  kITSnClsIbMin,        ///< Min. number of ITS clusters in the inner barrel
  kDCAxyMax,            ///< Max. |DCA_xy| (cm) as a function of pT
  kDCAzMax,             ///< Max. |DCA_z| (cm) as a function of pT

  /// track pid cuts
  kItsElectron, ///< ITS Electon PID
  kItsPion,     ///< ITS Pion PID
  kItsKaon,     ///< ITS Kaon PID
  kItsProton,   ///< ITS Proton PID
  kItsDeuteron, ///< ITS Deuteron PID
  kItsTriton,   ///< ITS Triton PID
  kItsHelium,   ///< ITS He3 PID

  kTpcElectron, ///< TPC Electon PID
  kTpcPion,     ///< TPC Pion PID
  kTpcKaon,     ///< TPC Kaon PID
  kTpcProton,   ///< TPC Proton PID
  kTpcDeuteron, ///< TPC Deuteron PID
  kTpcTriton,   ///< TPC Triton PID
  kTpcHelium,   ///< TPC He3 PID

  kTofElectron, ///< TOF Electon PID
  kTofPion,     ///< TOF Pion PID
  kTofKaon,     ///< TOF Kaon PID
  kTofProton,   ///< TOF Proton PID
  kTofDeuteron, ///< TOF Deuteron PID
  kTofTriton,   ///< TOF Triton PID
  kTofHelium,   ///< TOF He3 PID

  kTpcitsElectron, ///< TPC+ITS Electon PID
  kTpcitsPion,     ///< TPC+ITS Pion PID
  kTpcitsKaon,     ///< TPC+ITS Kaon PID
  kTpcitsProton,   ///< TPC+ITS Proton PID
  kTpcitsDeuteron, ///< TPC+ITS Deuteron PID
  kTpcitsTriton,   ///< TPC+ITS Triton PID
  kTpcitsHelium,   ///< TPC+ITS He3 PID

  kTpctofElectron, ///< TPC+TOF Electon PID
  kTpctofPion,     ///< TPC+TOF Pion PID
  kTpctofKaon,     ///< TPC+TOF Kaon PID
  kTpctofProton,   ///< TPC+TOF Proton PID
  kTpctofDeuteron, ///< TPC+TOF Deuteron PID
  kTpctofTriton,   ///< TPC+TOF Triton PID
  kTpctofHelium,   ///< TPC+TOF He3 PID

  kTrackSelsMax
};

constexpr char TrackSelHistName[] = "hTrackSelection";
constexpr char TrackSelsName[] = "Track Selection Object";
const std::unordered_map<TrackSels, std::string> trackSelectionNames = {
  {kTPCnClsMin, "Min. number of TPC clusters"},
  {kTPCcRowsMin, "Min. number of crossed TPC rows"},
  {kTPCnClsOvercRowsMin, "Min. fraction of TPC clusters over TPC crossed rows"},
  {kTPCsClsMax, "Max. number of shared TPC clusters"},
  {kTPCsClsMax, "Max. number of shared TPC clusters"},
  {kTPCsClsFracMax, "Max. fractions of shared TPC clusters"},
  {kITSnClsMin, "Min. number of ITS clusters"},
  {kITSnClsIbMin, "Min. number of ITS clusters in the inner barrel"},
  {kDCAxyMax, "Max. |DCA_xy| (cm) as a function of pT"},
  {kDCAzMax, "Max. |DCA_z| (cm) as a function of pT"},

  {kItsElectron, "ITS Electron PID"},
  {kItsPion, "ITS Pion PID"},
  {kItsKaon, "ITS Kaon PID"},
  {kItsProton, "ITS Proton PID"},
  {kItsDeuteron, "ITS Deuteron PID"},
  {kItsTriton, "ITS Triton PID"},
  {kItsHelium, "ITS He3 PID"},

  {kTpcElectron, "TPC Electron PID"},
  {kTpcPion, "TPC Pion PID"},
  {kTpcKaon, "TPC Kaon PID"},
  {kTpcProton, "TPC Proton PID"},
  {kTpcDeuteron, "TPC Deuteron PID"},
  {kTpcTriton, "TPC Triton PID"},
  {kTpcHelium, "TPC He3 PID"},

  {kTofElectron, "TOF Electron PID"},
  {kTofPion, "TOF Pion PID"},
  {kTofKaon, "TOF Kaon PID"},
  {kTofProton, "TOF Proton PID"},
  {kTofDeuteron, "TOF Deuteron PID"},
  {kTofTriton, "TOF Triton PID"},
  {kTofHelium, "TOF He3 PID"},

  {kTpcitsElectron, "TPC+ITS Electron PID"},
  {kTpcitsPion, "TPC+ITS Pion PID"},
  {kTpcitsKaon, "TPC+ITS Kaon PID"},
  {kTpcitsProton, "TPC+ITS Proton PID"},
  {kTpcitsDeuteron, "TPC+ITS Deuteron PID"},
  {kTpcitsTriton, "TPC+ITS Triton PID"},
  {kTpcitsHelium, "TPC+ITS He PID"},

  {kTpctofElectron, "TPC+TOF Electron PID"},
  {kTpctofPion, "TPC+TOF Pion PID"},
  {kTpctofKaon, "TPC+TOF Kaon PID"},
  {kTpctofProton, "TPC+TOF Proton PID"},
  {kTpctofDeuteron, "TPC+TOF Deuteron PID"},
  {kTpctofTriton, "TPC+TOF Triton PID"},
  {kTpctofHelium, "TPC+TOF He3 PID"}};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
template <const char* HistName>
class TrackSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackMaskType, kTrackSelsMax>
{
 public:
  TrackSelection() = default;
  ~TrackSelection() = default;

  template <typename T1, typename T2>
  void configure(o2::framework::HistogramRegistry* registry, T1& config, T2& filter)
  {
    mPtMin = filter.ptMin.value;
    mPtMax = filter.ptMax.value;
    mEtaMin = filter.etaMin.value;
    mEtaMax = filter.etaMax.value;
    mPhiMin = filter.phiMin.value;
    mPhiMax = filter.phiMax.value;

    // add selections for track quality
    this->addSelection(kTPCnClsMin, trackSelectionNames.at(kTPCnClsMin), config.tpcClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCcRowsMin, trackSelectionNames.at(kTPCcRowsMin), config.tpcCrossedRowsMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCnClsOvercRowsMin, trackSelectionNames.at(kTPCnClsOvercRowsMin), config.tpcClustersOverCrossedRows.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCsClsMax, trackSelectionNames.at(kTPCsClsMax), config.tpcSharedClustersMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kTPCsClsFracMax, trackSelectionNames.at(kTPCsClsFracMax), config.tpcSharedClusterFractionMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kITSnClsMin, trackSelectionNames.at(kITSnClsMin), config.itsClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kITSnClsIbMin, trackSelectionNames.at(kITSnClsIbMin), config.itsIbClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDCAxyMax, trackSelectionNames.at(kDCAxyMax), filter.ptMin.value, filter.ptMax.value, config.dcaxyMax.value, limits::kAbsUpperFunctionLimit, true, true, false);
    this->addSelection(kDCAzMax, trackSelectionNames.at(kDCAzMax), filter.ptMin.value, filter.ptMax.value, config.dcazMax.value, limits::kAbsUpperFunctionLimit, true, true, false);

    // add selections for Electron pid
    this->addSelection(kItsElectron, trackSelectionNames.at(kItsElectron), config.itsElectron.value, false, false, config.requirePidElectron);
    this->addSelection(kTpcElectron, trackSelectionNames.at(kTpcElectron), config.tpcElectron.value, false, false, config.requirePidElectron);
    this->addSelection(kTofElectron, trackSelectionNames.at(kTofElectron), config.tofElectron.value, false, false, config.requirePidElectron);
    this->addSelection(kTpcitsElectron, trackSelectionNames.at(kTpcitsElectron), config.tpcitsElectron.value, limits::kUpperLimit, false, false, config.requirePidElectron);
    this->addSelection(kTpctofElectron, trackSelectionNames.at(kTpctofElectron), config.tpctofElectron.value, limits::kUpperLimit, false, false, config.requirePidElectron);
    mElectronTofThres = config.minMomTofElectron.value;

    // add selections for Pion pid
    this->addSelection(kItsPion, trackSelectionNames.at(kItsPion), config.itsPion.value, false, false, config.requirePidPion);
    this->addSelection(kTpcPion, trackSelectionNames.at(kTpcPion), config.tpcPion.value, false, false, config.requirePidPion);
    this->addSelection(kTofPion, trackSelectionNames.at(kTofPion), config.tofPion.value, false, false, config.requirePidPion);
    this->addSelection(kTpcitsPion, trackSelectionNames.at(kTpcitsPion), config.tpcitsPion.value, limits::kUpperLimit, false, false, config.requirePidPion);
    this->addSelection(kTpctofPion, trackSelectionNames.at(kTpctofPion), config.tpctofPion.value, limits::kUpperLimit, false, false, config.requirePidPion);
    mPionTofThres = config.minMomTofPion.value;

    // add selections for Kaon pid
    this->addSelection(kItsKaon, trackSelectionNames.at(kItsKaon), config.itsKaon.value, false, false, config.requirePidKaon);
    this->addSelection(kTpcKaon, trackSelectionNames.at(kTpcKaon), config.tpcKaon.value, false, false, config.requirePidKaon);
    this->addSelection(kTofKaon, trackSelectionNames.at(kTofKaon), config.tofKaon.value, false, false, config.requirePidKaon);
    this->addSelection(kTpcitsKaon, trackSelectionNames.at(kTpcitsKaon), config.tpcitsKaon.value, limits::kUpperLimit, false, false, config.requirePidKaon);
    this->addSelection(kTpctofKaon, trackSelectionNames.at(kTpctofKaon), config.tpctofKaon.value, limits::kUpperLimit, false, false, config.requirePidKaon);
    mKaonTofThres = config.minMomTofKaon.value;

    // add selections for Proton pid
    this->addSelection(kItsProton, trackSelectionNames.at(kItsProton), config.itsProton.value, false, false, config.requirePidProton);
    this->addSelection(kTpcProton, trackSelectionNames.at(kTpcProton), config.tpcProton.value, false, false, config.requirePidProton);
    this->addSelection(kTofProton, trackSelectionNames.at(kTofProton), config.tofProton.value, false, false, config.requirePidProton);
    this->addSelection(kTpcitsProton, trackSelectionNames.at(kTpcitsProton), config.tpcitsProton.value, limits::kUpperLimit, false, false, config.requirePidProton);
    this->addSelection(kTpctofProton, trackSelectionNames.at(kTpctofProton), config.tpctofProton.value, limits::kUpperLimit, false, false, config.requirePidProton);
    mProtonTofThres = config.minMomTofProton.value;

    // add selections for Deuteron pid
    this->addSelection(kItsDeuteron, trackSelectionNames.at(kItsDeuteron), config.itsDeuteron.value, false, false, config.requirePidDeuteron);
    this->addSelection(kTpcDeuteron, trackSelectionNames.at(kTpcDeuteron), config.tpcDeuteron.value, false, false, config.requirePidDeuteron);
    this->addSelection(kTofDeuteron, trackSelectionNames.at(kTofDeuteron), config.tofDeuteron.value, false, false, config.requirePidDeuteron);
    this->addSelection(kTpcitsDeuteron, trackSelectionNames.at(kTpcitsDeuteron), config.tpcitsDeuteron.value, limits::kUpperLimit, false, false, config.requirePidDeuteron);
    this->addSelection(kTpctofDeuteron, trackSelectionNames.at(kTpctofDeuteron), config.tpctofDeuteron.value, limits::kUpperLimit, false, false, config.requirePidDeuteron);
    mDeuteronTofThres = config.minMomTofDeuteron.value;

    // add selections for Triton pid
    this->addSelection(kItsTriton, trackSelectionNames.at(kItsTriton), config.itsTriton.value, false, false, config.requirePidTriton);
    this->addSelection(kTpcTriton, trackSelectionNames.at(kTpcTriton), config.tpcTriton.value, false, false, config.requirePidTriton);
    this->addSelection(kTofTriton, trackSelectionNames.at(kTofTriton), config.tofTriton.value, false, false, config.requirePidTriton);
    this->addSelection(kTpcitsTriton, trackSelectionNames.at(kTpcitsTriton), config.tpcitsTriton.value, limits::kUpperLimit, false, false, config.requirePidTriton);
    this->addSelection(kTpctofTriton, trackSelectionNames.at(kTpctofTriton), config.tpctofTriton.value, limits::kUpperLimit, false, false, config.requirePidTriton);
    mTritonTofThres = config.minMomTofTriton.value;

    // add selections for Helium pid
    this->addSelection(kItsHelium, trackSelectionNames.at(kItsHelium), config.itsHelium.value, false, false, config.requirePidHelium);
    this->addSelection(kTpcHelium, trackSelectionNames.at(kTpcHelium), config.tpcHelium.value, false, false, config.requirePidHelium);
    this->addSelection(kTofHelium, trackSelectionNames.at(kTofHelium), config.tofHelium.value, false, false, config.requirePidHelium);
    this->addSelection(kTpcitsHelium, trackSelectionNames.at(kTpcitsHelium), config.tpcitsHelium.value, limits::kUpperLimit, false, false, config.requirePidHelium);
    this->addSelection(kTpctofHelium, trackSelectionNames.at(kTpctofHelium), config.tpctofHelium.value, limits::kUpperLimit, false, false, config.requirePidHelium);
    mHeliumTofThres = config.minMomTofHelium.value;

    this->setupContainers<HistName>(registry);
  }

  template <typename T>
  bool checkFilters(T const& track) const
  {
    return ((track.pt() > mPtMin && track.pt() < mPtMax) &&
            (track.eta() > mEtaMin && track.eta() < mEtaMax) &&
            (track.phi() > mPhiMin && track.phi() < mPhiMax));
  }

  template <typename T1>
  void evaluatePid(T1 const& Track,
                   float tofThreshold,
                   float nsigmaIts,
                   float nsigmaTpc,
                   float nsigmaTof,
                   TrackSels its,
                   TrackSels tpc,
                   TrackSels tof,
                   TrackSels tpcits,
                   TrackSels tpctof,
                   bool ignoreThreshold = false)
  {
    // above threshold without TOF: skip entirely unless forced
    // forced evaluation is used in the second pass to populate rejection bits
    if (!ignoreThreshold && Track.p() >= tofThreshold && !Track.hasTOF()) {
      return;
    }
    this->evaluateObservable(its, nsigmaIts);
    this->evaluateObservable(tpc, nsigmaTpc);
    this->evaluateObservable(tpcits, std::hypot(nsigmaTpc, nsigmaIts));
    if (Track.hasTOF()) {
      this->evaluateObservable(tof, nsigmaTof);
      this->evaluateObservable(tpctof, std::hypot(nsigmaTpc, nsigmaTof));
    }
  }

  template <typename T>
  void applySelections(T const& Track)
  {
    this->reset();
    this->evaluateObservable(kTPCnClsMin, Track.tpcNClsFound());
    this->evaluateObservable(kTPCcRowsMin, Track.tpcNClsCrossedRows());
    this->evaluateObservable(kTPCnClsOvercRowsMin, static_cast<float>(Track.tpcNClsFound()) / static_cast<float>(Track.tpcNClsCrossedRows()));
    this->evaluateObservable(kTPCsClsMax, Track.tpcNClsShared());
    this->evaluateObservable(kTPCsClsFracMax, static_cast<float>(Track.tpcNClsShared()) / static_cast<float>(Track.tpcNClsFound()));
    this->evaluateObservable(kITSnClsMin, Track.itsNCls());
    this->evaluateObservable(kITSnClsIbMin, Track.itsNClsInnerBarrel());
    this->updateLimits(kDCAxyMax, Track.pt());
    this->evaluateObservable(kDCAxyMax, Track.dcaXY());
    this->updateLimits(kDCAzMax, Track.pt());
    this->evaluateObservable(kDCAzMax, Track.dcaZ());

    // first pass: threshold-aware PID evaluation
    // determines if the track passes any optional selection and should be stored
    this->evaluatePid(Track, mElectronTofThres, Track.itsNSigmaEl(), Track.tpcNSigmaEl(), Track.tofNSigmaEl(), kItsElectron, kTpcElectron, kTofElectron, kTpcitsElectron, kTpctofElectron);
    this->evaluatePid(Track, mPionTofThres, Track.itsNSigmaPi(), Track.tpcNSigmaPi(), Track.tofNSigmaPi(), kItsPion, kTpcPion, kTofPion, kTpcitsPion, kTpctofPion);
    this->evaluatePid(Track, mKaonTofThres, Track.itsNSigmaKa(), Track.tpcNSigmaKa(), Track.tofNSigmaKa(), kItsKaon, kTpcKaon, kTofKaon, kTpcitsKaon, kTpctofKaon);
    this->evaluatePid(Track, mProtonTofThres, Track.itsNSigmaPr(), Track.tpcNSigmaPr(), Track.tofNSigmaPr(), kItsProton, kTpcProton, kTofProton, kTpcitsProton, kTpctofProton);
    this->evaluatePid(Track, mDeuteronTofThres, Track.itsNSigmaDe(), Track.tpcNSigmaDe(), Track.tofNSigmaDe(), kItsDeuteron, kTpcDeuteron, kTofDeuteron, kTpcitsDeuteron, kTpctofDeuteron);
    this->evaluatePid(Track, mTritonTofThres, Track.itsNSigmaTr(), Track.tpcNSigmaTr(), Track.tofNSigmaTr(), kItsTriton, kTpcTriton, kTofTriton, kTpcitsTriton, kTpctofTriton);
    this->evaluatePid(Track, mHeliumTofThres, Track.itsNSigmaHe(), Track.tpcNSigmaHe(), Track.tofNSigmaHe(), kItsHelium, kTpcHelium, kTofHelium, kTpcitsHelium, kTpctofHelium);

    // second pass: if the track passes minimal and any optional selection,
    // re-evaluate all species ignoring thresholds so rejection bits are fully
    // populated for all competing hypotheses. evaluate() resets each container
    // before writing, so no explicit reset is needed before this pass.
    if (this->passesAllRequiredSelections()) {
      this->evaluatePid(Track, mElectronTofThres, Track.itsNSigmaEl(), Track.tpcNSigmaEl(), Track.tofNSigmaEl(), kItsElectron, kTpcElectron, kTofElectron, kTpcitsElectron, kTpctofElectron, true);
      this->evaluatePid(Track, mPionTofThres, Track.itsNSigmaPi(), Track.tpcNSigmaPi(), Track.tofNSigmaPi(), kItsPion, kTpcPion, kTofPion, kTpcitsPion, kTpctofPion, true);
      this->evaluatePid(Track, mKaonTofThres, Track.itsNSigmaKa(), Track.tpcNSigmaKa(), Track.tofNSigmaKa(), kItsKaon, kTpcKaon, kTofKaon, kTpcitsKaon, kTpctofKaon, true);
      this->evaluatePid(Track, mProtonTofThres, Track.itsNSigmaPr(), Track.tpcNSigmaPr(), Track.tofNSigmaPr(), kItsProton, kTpcProton, kTofProton, kTpcitsProton, kTpctofProton, true);
      this->evaluatePid(Track, mDeuteronTofThres, Track.itsNSigmaDe(), Track.tpcNSigmaDe(), Track.tofNSigmaDe(), kItsDeuteron, kTpcDeuteron, kTofDeuteron, kTpcitsDeuteron, kTpctofDeuteron, true);
      this->evaluatePid(Track, mTritonTofThres, Track.itsNSigmaTr(), Track.tpcNSigmaTr(), Track.tofNSigmaTr(), kItsTriton, kTpcTriton, kTofTriton, kTpcitsTriton, kTpctofTriton, true);
      this->evaluatePid(Track, mHeliumTofThres, Track.itsNSigmaHe(), Track.tpcNSigmaHe(), Track.tofNSigmaHe(), kItsHelium, kTpcHelium, kTofHelium, kTpcitsHelium, kTpctofHelium, true);
    }

    this->assembleBitmask<HistName>();
  }

 protected:
  float mElectronTofThres = 99.f;
  float mPionTofThres = 99.f;
  float mKaonTofThres = 99.f;
  float mProtonTofThres = 99.f;
  float mDeuteronTofThres = 99.f;
  float mTritonTofThres = 99.f;
  float mHeliumTofThres = 99.f;

  float mPtMin = 0.f;
  float mPtMax = 99.f;
  float mEtaMin = -0.9;
  float mEtaMax = 0.9;
  float mPhiMin = 0;
  float mPhiMax = o2::constants::math::TwoPI;
};

struct TrackBuilderProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::FTracks> producedTracks;
  o2::framework::Produces<o2::aod::FTrackMass> producedTrackMass;
  o2::framework::Produces<o2::aod::FTrackMasks> producedTrackMasks;
  o2::framework::Produces<o2::aod::FTrackDcas> producedTrackDcas;
  o2::framework::Produces<o2::aod::FTrackExtras> producedTrackExtras;
  o2::framework::Produces<o2::aod::FElectronPids> producedElectronPids;
  o2::framework::Produces<o2::aod::FPionPids> producedPionPids;
  o2::framework::Produces<o2::aod::FKaonPids> producedKaonPids;
  o2::framework::Produces<o2::aod::FProtonPids> producedProtonPids;
  o2::framework::Produces<o2::aod::FDeuteronPids> producedDeuteronPids;
  o2::framework::Produces<o2::aod::FTritonPids> producedTritonPids;
  o2::framework::Produces<o2::aod::FHeliumPids> producedHeliumPids;
};

struct ConfTrackTables : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackTables");
  o2::framework::Configurable<int> produceTracks{"produceTracks", -1, "Produce Tracks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceTrackMasks{"produceTrackMasks", -1, "Produce TrackMasks (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceTrackMass{"produceTrackMass", -1, "Produce TrackMass (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceTrackDcas{"produceTrackDcas", -1, "Produce TrackDcas (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceTrackExtras{"produceTrackExtras", -1, "Produce TrackExtras (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceElectronPids{"produceElectronPids", -1, "Produce ElectronPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> producePionPids{"producePionPids", -1, "Produce PionPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceKaonPids{"produceKaonPids", -1, "Produce KaonPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceProtonPids{"produceProtonPids", -1, "Produce ProtonPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceDeuteronPids{"produceDeuteronPids", -1, "Produce DeuteronPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceTritonPids{"produceTritonPids", -1, "Produce TritonPids (-1: auto; 0 off; 1 on)"};
  o2::framework::Configurable<int> produceHeliumPids{"produceHeliumPids", -1, "Produce HeliumPids (-1: auto; 0 off; 1 on)"};
};

template <const char* HistName>
class TrackBuilder
{
 public:
  TrackBuilder() = default;
  ~TrackBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(o2::framework::HistogramRegistry* registry, T1& config, T2& filter, T3& table, T4& initContext)
  {
    LOG(info) << "Initialize femto track builder...";

    mProduceTracks = utils::enableTable("FTracks_001", table.produceTracks.value, initContext);
    mProduceTrackMasks = utils::enableTable("FTrackMasks_001", table.produceTrackMasks.value, initContext);
    mProduceTrackMass = utils::enableTable("FTrackMass_001", table.produceTrackMass.value, initContext);
    mProduceTrackDcas = utils::enableTable("FTrackDcas_001", table.produceTrackDcas.value, initContext);
    mProduceTrackExtras = utils::enableTable("FTrackExtras_001", table.produceTrackExtras.value, initContext);
    mProduceElectronPids = utils::enableTable("FElectronPids_001", table.produceElectronPids.value, initContext);
    mProducePionPids = utils::enableTable("FPionPids_001", table.producePionPids.value, initContext);
    mProduceKaonPids = utils::enableTable("FKaonPids_001", table.produceKaonPids.value, initContext);
    mProduceProtonPids = utils::enableTable("FProtonPids_001", table.produceProtonPids.value, initContext);
    mProduceDeuteronPids = utils::enableTable("FDeuteronPids_001", table.produceDeuteronPids.value, initContext);
    mProduceTritonPids = utils::enableTable("FTritonPids_001", table.produceTritonPids.value, initContext);
    mProduceHeliumPids = utils::enableTable("FHeliumPids_001", table.produceHeliumPids.value, initContext);

    if (mProduceTracks || mProduceTrackMasks || mProduceTrackMass || mProduceTrackDcas || mProduceTrackExtras || mProduceElectronPids || mProducePionPids || mProduceKaonPids || mProduceProtonPids || mProduceDeuteronPids || mProduceTritonPids || mProduceHeliumPids) {
      mFillAnyTable = true;
    } else {
      LOG(info) << "No tables configured, Selection object will not be configured...";
      LOG(info) << "Initialization done...";
      return;
    }
    mTrackSelection.configure(registry, config, filter);
    mTrackSelection.printSelections(TrackSelsName);
    LOG(info) << "Initialization done...";
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fillTracks(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& tracks, T5& trackProducts)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& track : tracks) {
      if (!mTrackSelection.checkFilters(track)) {
        continue;
      }
      mTrackSelection.applySelections(track);
      if (!mTrackSelection.passesAllRequiredSelections()) {
        continue;
      }

      collisionBuilder.template fillCollision<system>(collisionProducts, col);
      this->fillTrack<modes::Track::kTrack>(track, trackProducts, collisionProducts);
    }
  }

  template <modes::Track type, typename T1, typename T2, typename T3>
  void fillTrack(T1 const& track, T2& trackProducts, T3& collisionProducts)
  {
    if (mProduceTracks) {
      trackProducts.producedTracks(collisionProducts.producedCollision.lastIndex(),
                                   track.pt() * track.sign(),
                                   track.eta(),
                                   track.phi());
      indexMap.emplace(track.globalIndex(), trackProducts.producedTracks.lastIndex());
    }
    if (mProduceTrackMasks) {
      if constexpr (type == modes::Track::kTrack) {
        trackProducts.producedTrackMasks(mTrackSelection.getBitmask());
      } else {
        trackProducts.producedTrackMasks(static_cast<o2::aod::femtodatatypes::TrackMaskType>(0u));
      }
    }
    if (mProduceTrackMass) {
      trackProducts.producedTrackMass(track.mass());
    }
    if (mProduceTrackDcas) {
      trackProducts.producedTrackDcas(track.dcaXY(), track.dcaZ());
    }
    if (mProduceTrackExtras) {
      trackProducts.producedTrackExtras(track.isPVContributor(),
                                        track.itsNCls(),
                                        track.itsNClsInnerBarrel(),
                                        track.itsChi2NCl(),
                                        track.itsClusterSizes(),
                                        track.tpcSignal(),
                                        track.tpcInnerParam(),
                                        track.tpcNClsFound(),
                                        track.tpcNClsCrossedRows(),
                                        track.tpcNClsShared(),
                                        track.beta());
    }
    if (mProduceElectronPids) {
      trackProducts.producedElectronPids(track.itsNSigmaEl(), track.tpcNSigmaEl(), track.tofNSigmaEl());
    }
    if (mProducePionPids) {
      trackProducts.producedPionPids(track.itsNSigmaPi(), track.tpcNSigmaPi(), track.tofNSigmaPi());
    }
    if (mProduceKaonPids) {
      trackProducts.producedKaonPids(track.itsNSigmaKa(), track.tpcNSigmaKa(), track.tofNSigmaKa());
    }
    if (mProduceProtonPids) {
      trackProducts.producedProtonPids(track.itsNSigmaPr(), track.tpcNSigmaPr(), track.tofNSigmaPr());
    }
    if (mProduceDeuteronPids) {
      trackProducts.producedDeuteronPids(track.itsNSigmaDe(), track.tpcNSigmaDe(), track.tofNSigmaDe());
    }
    if (mProduceTritonPids) {
      trackProducts.producedTritonPids(track.itsNSigmaTr(), track.tpcNSigmaTr(), track.tofNSigmaTr());
    }
    if (mProduceHeliumPids) {
      trackProducts.producedHeliumPids(track.itsNSigmaHe(), track.tpcNSigmaHe(), track.tofNSigmaHe());
    }
  }

  template <modes::System system, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  void fillMcTracks(T1 const& col, T2& collisionBuilder, T3& collisionProducts, T4 const& mcCols, T5 const& tracks, T6 const& tracksWithItsPid, T7& trackProducts, T8 const& mcParticles, T9& mcBuilder, T10& mcProducts)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& trackWithItsPid : tracksWithItsPid) {
      if (!mTrackSelection.checkFilters(trackWithItsPid)) {
        continue;
      }
      mTrackSelection.applySelections(trackWithItsPid);
      if (!mTrackSelection.passesAllRequiredSelections()) {
        continue;
      }
      collisionBuilder.template fillMcCollision<system>(collisionProducts, col, mcCols, mcProducts, mcBuilder);
      // get track from the track table so we can dereference mc particle properly
      auto track = tracks.iteratorAt(trackWithItsPid.index());
      this->template fillMcTrack<system, modes::Track::kTrack>(col, collisionProducts, mcCols, track, trackWithItsPid, trackProducts, mcParticles, mcBuilder, mcProducts);
    }
  }

  template <modes::System system, modes::Track trackType, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  void fillMcTrack(T1 const& col, T2& collisionProducts, T3 const& mcCols, T4 const& track, T5 const& trackWithItsPid, T6& trackProducts, T7 const& mcParticles, T8& mcBuilder, T9& mcProducts)
  {
    if (!mProduceTracks) {
      return;
    }
    this->template fillTrack<trackType>(trackWithItsPid, trackProducts, collisionProducts);
    mcBuilder.template fillMcTrackWithLabel<system>(col, mcCols, track, mcParticles, mcProducts);
  }

  template <modes::Track type, typename T1, typename T2, typename T3>
  int64_t getDaughterIndex(const T1& daughter, T2& trackProducts, T3& collisionProducts)
  {
    auto result = utils::getIndex(daughter.globalIndex(), indexMap);
    if (result) {
      return result.value();
    } else {
      this->fillTrack<type>(daughter, trackProducts, collisionProducts);
      int64_t idx = trackProducts.producedTracks.lastIndex();
      return idx;
    }
  }

  template <modes::System system, modes::Track type, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  int64_t getDaughterIndex(const T1& col, T2& collisionProducts, T3 const& mcCols, const T4& daughter, T5& daughterWithItsPid, T6& trackProducts, T7 const& mcParticles, T8& mcBuilder, T9& mcProducts)
  {
    auto result = utils::getIndex(daughter.globalIndex(), indexMap);
    if (result) {
      // daugher already in track table
      return result.value();
    } else {
      this->fillMcTrack<system, type>(col, collisionProducts, mcCols, daughter, daughterWithItsPid, trackProducts, mcParticles, mcBuilder, mcProducts);
      // daughter is last track which was added added
      return trackProducts.producedTracks.lastIndex();
    }
  }

  template <typename T>
  void reset(T const& tracks)
  {
    indexMap.clear();
    indexMap.reserve(tracks.size());
  };

 private:
  TrackSelection<HistName> mTrackSelection;
  bool mFillAnyTable = false;
  bool mProduceTracks = false;
  bool mProduceTrackMasks = false;
  bool mProduceTrackMass = false;
  bool mProduceTrackDcas = false;
  bool mProduceTrackExtras = false;
  bool mProduceElectronPids = false;
  bool mProducePionPids = false;
  bool mProduceKaonPids = false;
  bool mProduceProtonPids = false;
  bool mProduceDeuteronPids = false;
  bool mProduceTritonPids = false;
  bool mProduceHeliumPids = false;

  std::unordered_map<int64_t, int64_t> indexMap; // for mapping tracks to daughers of lambdas, cascades and resonances ...
};

struct TrackBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFTracks> producedTracks;
  o2::framework::Produces<o2::aod::StoredFTrackMasks> producedTrackMasks;
};

struct ConfTrackTablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackTables");
  o2::framework::Configurable<int> limitTrack1{"limitTrack1", 1, "At least this many tracks of type 1 need to be in the collision. Ignored if set to 0."};
  o2::framework::Configurable<int> limitTrack2{"limitTrack2", 0, "At least this many tracks of type 2 need to be in the collision. Ignored if set to 0."};
};

class TrackBuilderDerivedToDerived
{
 public:
  TrackBuilderDerivedToDerived() = default;
  ~TrackBuilderDerivedToDerived() = default;

  template <typename T>
  void init(T& config)
  {
    mLimitTrack1 = config.limitTrack1.value;
    mLimitTrack2 = config.limitTrack2.value;

    if (mLimitTrack1 == 0 && mLimitTrack2 == 0) {
      LOG(fatal) << "Both track limits are 0. Breaking...";
    }
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool collisionHasTooFewTracks(T1& col, T2& /*trackTable*/, T3& partitionTrack1, T4& partitionTrack2, T5& cache)
  {
    auto trackSlice1 = partitionTrack1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto trackSlice2 = partitionTrack2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    return trackSlice1.size() < mLimitTrack1 || trackSlice2.size() < mLimitTrack2;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  void processTracks(T1& col, T2& /*trackTable*/, T3& partitionTrack1, T4& partitionTrack2, T5& cache, T6& newTrackTable, T7& newCollisionTable)
  {
    indexMap.clear();

    if (mLimitTrack1 > 0) {
      auto trackSlice1 = partitionTrack1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      for (auto const& track : trackSlice1) {
        this->fillTrack(track, newTrackTable, newCollisionTable);
      }
    }

    if (mLimitTrack2 > 0) {
      auto trackSlice2 = partitionTrack2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
      for (auto const& track : trackSlice2) {
        this->fillTrack(track, newTrackTable, newCollisionTable);
      }
    }
  }

  template <typename T1, typename T2, typename T3>
  void fillTrack(T1 const& track, T2& trackProducts, T3& collisionProducts)
  {
    if (indexMap.find(track.globalIndex()) == indexMap.end()) { // protect against double filling
      trackProducts.producedTracks(collisionProducts.producedCollision.lastIndex(),
                                   track.signedPt(),
                                   track.eta(),
                                   track.phi());
      trackProducts.producedTrackMasks(track.mask());
      indexMap.emplace(track.globalIndex(), trackProducts.producedTracks.lastIndex());
    }
  }

  template <typename T1, typename T2, typename T3>
  int64_t getDaughterIndex(const T1& daughter, T2& trackProducts, T3& collisionProducts)
  {
    auto result = utils::getIndex(daughter.globalIndex(), indexMap);
    if (result) {
      return result.value();
    } else {
      this->fillTrack(daughter, trackProducts, collisionProducts);
      int64_t idx = trackProducts.producedTracks.lastIndex();
      return idx;
    }
  }

 private:
  int mLimitTrack1 = 0;
  int mLimitTrack2 = 0;

  std::unordered_map<int64_t, int64_t> indexMap; // for mapping tracks to daughers of lambdas, cascades and resonances ...
};

} // namespace trackbuilder
//
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRACKBUILDER_H_
