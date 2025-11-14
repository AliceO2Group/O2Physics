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
  o2::framework::Configurable<bool> pidIsOptionalElectron{"pidIsOptionalElectron", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofElectron{"minMomTofElectron", 0.3, "Minimum momentum to required TOF PID for Electron"};
  o2::framework::Configurable<std::vector<float>> itsElectron{"itsElectron", {}, "Maximum |nsigma| for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tpcElectron{"tpcElectron", {}, "Maximum |nsigma| for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tofElectron{"tofElectron", {}, "Maximum |nsigma| for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsElectron{"tpcitsElectron", {}, "Maximum |nsigma| for Electron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofElectron{"tpctofElectron", {}, "Maximum |nsigma| for Electron PID"};

  // Pion PID cuts
  o2::framework::Configurable<bool> pidIsOptionalPion{"pidIsOptionalPion", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofPion{"minMomTofPion", 0.5, "Minimum momentum to required TOF PID for Pion"};
  o2::framework::Configurable<std::vector<float>> itsPion{"itsPion", {}, "Maximum |nsigma| for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tpcPion{"tpcPion", {}, "Maximum |nsigma| for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tofPion{"tofPion", {}, "Maximum |nsigma| for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsPion{"tpcitsPion", {}, "Maximum |nsigma| for Pion PID"};
  o2::framework::Configurable<std::vector<float>> tpctofPion{"tpctofPion", {}, "Maximum |nsigma| for Pion PID"};

  // Kaon PID cuts
  o2::framework::Configurable<bool> pidIsOptionalKaon{"pidIsOptionalKaon", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofKaon{"minMomTofKaon", 0.4, "Minimum momentum to required TOF PID for Kaon"};
  o2::framework::Configurable<std::vector<float>> itsKaon{"itsKaon", {}, "Maximum |nsigma| for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpcKaon{"tpcKaon", {}, "Maximum |nsigma| for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tofKaon{"tofKaon", {}, "Maximum |nsigma| for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsKaon{"tpcitsKaon", {}, "Maximum |nsigma| for Kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpctofKaon{"tpctofKaon", {}, "Maximum |nsigma| for Kaon PID"};

  // Proton PID cuts
  o2::framework::Configurable<bool> pidIsOptionalProton{"pidIsOptionalProton", true, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofProton{"minMomTofProton", 0.75, "Minimum momentum to required TOF PID for Proton"};
  o2::framework::Configurable<std::vector<float>> itsProton{"itsProton", {}, "Maximum |nsigma| for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tpcProton{"tpcProton", {}, "Maximum |nsigma| for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tofProton{"tofProton", {}, "Maximum |nsigma| for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsProton{"tpcitsProton", {3.f}, "Maximum |nsigma| for Proton PID"};
  o2::framework::Configurable<std::vector<float>> tpctofProton{"tpctofProton", {3.f}, "Maximum |nsigma| for Proton PID"};

  // Deuteron PID cuts
  o2::framework::Configurable<bool> pidIsOptionalDeuteron{"pidIsOptionalDeuteron", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofDeuteron{"minMomTofDeuteron", 1.2, "Minimum momentum to required TOF PID for Deuteron"};
  o2::framework::Configurable<std::vector<float>> itsDeuteron{"itsDeuteron", {}, "Maximum |nsigma| for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpcDeuteron{"tpcDeuteron", {}, "Maximum |nsigma| for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tofDeuteron{"tofDeuteron", {}, "Maximum |nsigma| for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsDeuteron{"tpcitsDeuteron", {}, "Maximum |nsigma| for Deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofDeuteron{"tpctofDeuteron", {}, "Maximum |nsigma| for Deuteron PID"};

  // Triton PID cuts
  o2::framework::Configurable<bool> pidIsOptionalTriton{"pidIsOptionalTriton", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofTriton{"minMomTofTriton", 1.4, "Minimum momentum to required TOF PID for Triton"};
  o2::framework::Configurable<std::vector<float>> itsTriton{"itsTriton", {}, "Maximum |nsigma| for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tpcTriton{"tpcTriton", {}, "Maximum |nsigma| for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tofTriton{"tofTriton", {}, "Maximum |nsigma| for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsTriton{"tpcitsTriton", {}, "Maximum |nsigma| for Triton PID"};
  o2::framework::Configurable<std::vector<float>> tpctofTriton{"tpctofTriton", {}, "Maximum |nsigma| for Triton PID"};

  // Helium PID cuts
  o2::framework::Configurable<bool> pidIsOptionalHelium{"pidIsOptionalHelium", false, "Make election PID optional"};
  o2::framework::Configurable<float> minMomTofHelium{"minMomTofHelium", 1.6, "Minimum momentum to required TOF PID for Helium"};
  o2::framework::Configurable<std::vector<float>> itsHelium{"itsHelium", {}, "Maximum |nsigma| for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tpcHelium{"tpcHelium", {}, "Maximum |nsigma| for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tofHelium{"tofHelium", {}, "Maximum |nsigma| for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsHelium{"tpcitsHelium", {}, "Maximum |nsigma| for Helium PID"};
  o2::framework::Configurable<std::vector<float>> tpctofHelium{"tpctofHelium", {}, "Maximum |nsigma| for Helium PID"};
};

// define the template structure for TrackSelection
template <const char* Prefix>
struct ConfTrackSelection : public o2::framework::ConfigurableGroup {
  std::string prefix = Prefix; // Unique prefix based on the template argument
  // configuration parameters
  o2::framework::Configurable<int> pdgCode{"pdgCode", 2212, "Track PDG code"};
  o2::framework::Configurable<int> chargeAbs{"chargeAbs", 1, "Absolute value of charge (e.g. 1 for most tracks, 2 for He3)"};
  o2::framework::Configurable<int> chargeSign{"chargeSign", 1, "Track charge sign: +1 for positive, -1 for negative, 0 for both"};
  // filters for kinematics
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT (GeV/c)"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT (GeV/c)"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // track selection masks
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskLowMomentum{"maskLowMomentum", 0x2u, "Bitmask for selections below momentum threshold"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskHighMomentum{"maskHighMomentum", 0x1u, "Bitmask for selections above momentum threshold"};
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
    mPtMin = filter.ptMin;
    mPtMax = filter.ptMax;
    mEtaMin = filter.etaMin;
    mEtaMax = filter.etaMax;
    mPhiMin = filter.phiMin;
    mPhiMax = filter.phiMax;

    // add selections for track quality
    this->addSelection(kTPCnClsMin, trackSelectionNames.at(kTPCnClsMin), config.tpcClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCcRowsMin, trackSelectionNames.at(kTPCcRowsMin), config.tpcCrossedRowsMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCnClsOvercRowsMin, trackSelectionNames.at(kTPCnClsOvercRowsMin), config.tpcClustersOverCrossedRows.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kTPCsClsMax, trackSelectionNames.at(kTPCsClsMax), config.tpcSharedClustersMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kTPCsClsFracMax, trackSelectionNames.at(kTPCsClsFracMax), config.tpcSharedClusterFractionMax.value, limits::kUpperLimit, true, true, false);
    this->addSelection(kITSnClsMin, trackSelectionNames.at(kITSnClsMin), config.itsClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kITSnClsIbMin, trackSelectionNames.at(kITSnClsIbMin), config.itsIbClustersMin.value, limits::kLowerLimit, true, true, false);
    this->addSelection(kDCAxyMax, trackSelectionNames.at(kDCAxyMax), filter.ptMin, filter.ptMax.value, config.dcaxyMax.value, limits::kAbsUpperFunctionLimit, true, true, false);
    this->addSelection(kDCAzMax, trackSelectionNames.at(kDCAzMax), filter.ptMin.value, filter.ptMax.value, config.dcazMax.value, limits::kAbsUpperFunctionLimit, true, true, false);

    // add selections for Electron pid
    this->addSelection(kItsElectron, trackSelectionNames.at(kItsElectron), config.itsElectron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalElectron);
    this->addSelection(kTpcElectron, trackSelectionNames.at(kTpcElectron), config.tpcElectron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalElectron);
    this->addSelection(kTofElectron, trackSelectionNames.at(kTofElectron), config.tofElectron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalElectron);
    this->addSelection(kTpcitsElectron, trackSelectionNames.at(kTpcitsElectron), config.tpcitsElectron.value, limits::kUpperLimit, false, false, config.pidIsOptionalElectron);
    this->addSelection(kTpctofElectron, trackSelectionNames.at(kTpctofElectron), config.tpctofElectron.value, limits::kUpperLimit, false, false, config.pidIsOptionalElectron);
    mElectronTofThres = config.minMomTofElectron.value;

    // add selections for Pion pid
    this->addSelection(kItsPion, trackSelectionNames.at(kItsPion), config.itsPion.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalPion);
    this->addSelection(kTpcPion, trackSelectionNames.at(kTpcPion), config.tpcPion.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalPion);
    this->addSelection(kTofPion, trackSelectionNames.at(kTofPion), config.tofPion.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalPion);
    this->addSelection(kTpcitsPion, trackSelectionNames.at(kTpcitsPion), config.tpcitsPion.value, limits::kUpperLimit, false, false, config.pidIsOptionalPion);
    this->addSelection(kTpctofPion, trackSelectionNames.at(kTpctofPion), config.tpctofPion.value, limits::kUpperLimit, false, false, config.pidIsOptionalPion);
    mPionTofThres = config.minMomTofPion.value;

    // add selections for Kaon pid
    this->addSelection(kItsKaon, trackSelectionNames.at(kItsKaon), config.itsKaon.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalKaon);
    this->addSelection(kTpcKaon, trackSelectionNames.at(kTpcKaon), config.tpcKaon.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalKaon);
    this->addSelection(kTofKaon, trackSelectionNames.at(kTofKaon), config.tofKaon.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalKaon);
    this->addSelection(kTpcitsKaon, trackSelectionNames.at(kTpcitsKaon), config.tpcitsKaon.value, limits::kUpperLimit, false, false, config.pidIsOptionalKaon);
    this->addSelection(kTpctofKaon, trackSelectionNames.at(kTpctofKaon), config.tpctofKaon.value, limits::kUpperLimit, false, false, config.pidIsOptionalKaon);
    mKaonTofThres = config.minMomTofKaon.value;

    // add selections for Proton pid
    this->addSelection(kItsProton, trackSelectionNames.at(kItsProton), config.itsProton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalProton);
    this->addSelection(kTpcProton, trackSelectionNames.at(kTpcProton), config.tpcProton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalProton);
    this->addSelection(kTofProton, trackSelectionNames.at(kTofProton), config.tofProton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalProton);
    this->addSelection(kTpcitsProton, trackSelectionNames.at(kTpcitsProton), config.tpcitsProton.value, limits::kUpperLimit, false, false, config.pidIsOptionalProton);
    this->addSelection(kTpctofProton, trackSelectionNames.at(kTpctofProton), config.tpctofProton.value, limits::kUpperLimit, false, false, config.pidIsOptionalProton);
    mProtonTofThres = config.minMomTofProton.value;

    // add selections for Deuteron pid
    this->addSelection(kItsDeuteron, trackSelectionNames.at(kItsDeuteron), config.itsDeuteron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalDeuteron);
    this->addSelection(kTpcDeuteron, trackSelectionNames.at(kTpcDeuteron), config.tpcDeuteron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalDeuteron);
    this->addSelection(kTofDeuteron, trackSelectionNames.at(kTofDeuteron), config.tofDeuteron.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalDeuteron);
    this->addSelection(kTpcitsDeuteron, trackSelectionNames.at(kTpcitsDeuteron), config.tpcitsDeuteron.value, limits::kUpperLimit, false, false, config.pidIsOptionalDeuteron);
    this->addSelection(kTpctofDeuteron, trackSelectionNames.at(kTpctofDeuteron), config.tpctofDeuteron.value, limits::kUpperLimit, false, false, config.pidIsOptionalDeuteron);
    mDeuteronTofThres = config.minMomTofDeuteron.value;

    // add selections for Triton pid
    this->addSelection(kItsTriton, trackSelectionNames.at(kItsTriton), config.itsTriton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalTriton);
    this->addSelection(kTpcTriton, trackSelectionNames.at(kTpcTriton), config.tpcTriton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalTriton);
    this->addSelection(kTofTriton, trackSelectionNames.at(kTofTriton), config.tofTriton.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalTriton);
    this->addSelection(kTpcitsTriton, trackSelectionNames.at(kTpcitsTriton), config.tpcitsTriton.value, limits::kUpperLimit, false, false, config.pidIsOptionalTriton);
    this->addSelection(kTpctofTriton, trackSelectionNames.at(kTpctofTriton), config.tpctofTriton.value, limits::kUpperLimit, false, false, config.pidIsOptionalTriton);
    mTritonTofThres = config.minMomTofTriton.value;

    // add selections for Helium pid
    this->addSelection(kItsHelium, trackSelectionNames.at(kItsHelium), config.itsHelium.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalHelium);
    this->addSelection(kTpcHelium, trackSelectionNames.at(kTpcHelium), config.tpcHelium.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalHelium);
    this->addSelection(kTofHelium, trackSelectionNames.at(kTofHelium), config.tofHelium.value, limits::kAbsUpperLimit, false, false, config.pidIsOptionalHelium);
    this->addSelection(kTpcitsHelium, trackSelectionNames.at(kTpcitsHelium), config.tpcitsHelium.value, limits::kUpperLimit, false, false, config.pidIsOptionalHelium);
    this->addSelection(kTpctofHelium, trackSelectionNames.at(kTpctofHelium), config.tpctofHelium.value, limits::kUpperLimit, false, false, config.pidIsOptionalHelium);
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
                   TrackSels tpctof)
  {
    // if track is below threshold, just check every PID
    if (Track.p() < tofThreshold) {
      this->evaluateObservable(its, nsigmaIts);
      this->evaluateObservable(tpc, nsigmaTpc);
      this->evaluateObservable(tpcits, std::hypot(nsigmaTpc, nsigmaIts));
      this->evaluateObservable(tof, nsigmaTof);
      this->evaluateObservable(tpctof, std::hypot(nsigmaTpc, nsigmaTof));
      return;
    }
    // if track is above threshold, check if TOF PID is available
    // if not, we dont check any selection and they stay at reseted values, i.e. the cut fails
    if (Track.hasTOF()) {
      // if tof inforamtion is available, check them first
      this->evaluateObservable(tof, nsigmaTof);
      this->evaluateObservable(tpctof, std::hypot(nsigmaTpc, nsigmaTof));
      // if both failed, the bitmask will be 0 and there is no need to check tpc and its information since we do not want to have this track
      // so if we just bail out here, the PID for this particle type will failed for its, tpc and tof
      if (this->passesOptionalSelection(tof) || this->passesOptionalSelection(tpctof)) {
        this->evaluateObservable(its, nsigmaIts);
        this->evaluateObservable(tpc, nsigmaTpc);
        this->evaluateObservable(tpcits, std::hypot(nsigmaTpc, nsigmaIts));
      }
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

    // evalue bitmask for pt dependent dca cuts
    this->updateLimits(kDCAxyMax, Track.pt());
    this->evaluateObservable(kDCAxyMax, Track.dcaXY());

    this->updateLimits(kDCAzMax, Track.pt());
    this->evaluateObservable(kDCAzMax, Track.dcaZ());

    this->evaluatePid(Track,
                      mElectronTofThres,
                      Track.itsNSigmaEl(),
                      Track.tpcNSigmaEl(),
                      Track.tofNSigmaEl(),
                      kItsElectron,
                      kTpcElectron,
                      kTofElectron,
                      kTpcitsElectron,
                      kTpctofElectron);

    this->evaluatePid(Track,
                      mPionTofThres,
                      Track.itsNSigmaPi(),
                      Track.tpcNSigmaPi(),
                      Track.tofNSigmaPi(),
                      kItsPion,
                      kTpcPion,
                      kTofPion,
                      kTpcitsPion,
                      kTpctofPion);

    this->evaluatePid(Track,
                      mKaonTofThres,
                      Track.itsNSigmaKa(),
                      Track.tpcNSigmaKa(),
                      Track.tofNSigmaKa(),
                      kItsKaon,
                      kTpcKaon,
                      kTofKaon,
                      kTpcitsKaon,
                      kTpctofKaon);

    this->evaluatePid(Track,
                      mProtonTofThres,
                      Track.itsNSigmaPr(),
                      Track.tpcNSigmaPr(),
                      Track.tofNSigmaPr(),
                      kItsProton,
                      kTpcProton,
                      kTofProton,
                      kTpcitsProton,
                      kTpctofProton);

    this->evaluatePid(Track,
                      mDeuteronTofThres,
                      Track.itsNSigmaDe(),
                      Track.tpcNSigmaDe(),
                      Track.tofNSigmaDe(),
                      kItsDeuteron,
                      kTpcDeuteron,
                      kTofDeuteron,
                      kTpcitsDeuteron,
                      kTpctofDeuteron);

    this->evaluatePid(Track,
                      mTritonTofThres,
                      Track.itsNSigmaTr(),
                      Track.tpcNSigmaTr(),
                      Track.tofNSigmaTr(),
                      kItsTriton,
                      kTpcTriton,
                      kTofTriton,
                      kTpcitsTriton,
                      kTpctofTriton);

    this->evaluatePid(Track,
                      mHeliumTofThres,
                      Track.itsNSigmaHe(),
                      Track.tpcNSigmaHe(),
                      Track.tofNSigmaHe(),
                      kItsHelium,
                      kTpcHelium,
                      kTofHelium,
                      kTpcitsHelium,
                      kTpctofHelium);

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
    mProduceTrackDcas = utils::enableTable("FTrackDcas_001", table.produceTrackDcas.value, initContext);
    mProduceTrackExtras = utils::enableTable("FTrackExtras_001", table.produceTrackExtras.value, initContext);
    mProduceElectronPids = utils::enableTable("FElectronPids_001", table.produceElectronPids.value, initContext);
    mProducePionPids = utils::enableTable("FPionPids_001", table.producePionPids.value, initContext);
    mProduceKaonPids = utils::enableTable("FKaonPids_001", table.produceKaonPids.value, initContext);
    mProduceProtonPids = utils::enableTable("FProtonPids_001", table.produceProtonPids.value, initContext);
    mProduceDeuteronPids = utils::enableTable("FDeuteronPids_001", table.produceDeuteronPids.value, initContext);
    mProduceTritonPids = utils::enableTable("FTritonPids_001", table.produceTritonPids.value, initContext);
    mProduceHeliumPids = utils::enableTable("FHeliumPids_001", table.produceHeliumPids.value, initContext);

    if (mProduceTracks || mProduceTrackMasks || mProduceTrackDcas || mProduceTrackExtras || mProduceElectronPids || mProducePionPids || mProduceKaonPids || mProduceProtonPids || mProduceDeuteronPids || mProduceTritonPids || mProduceHeliumPids) {
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

  template <typename T1, typename T2, typename T3, typename T4>
  void fillTracks(T1 const& tracks, T2& trackProducts, T3& collisionProducts, T4& indexMap)
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
      this->fillTrack<modes::Track::kPrimaryTrack>(track, trackProducts, collisionProducts, indexMap);
    }
  }

  template <modes::Track type, typename T1, typename T2, typename T3, typename T4>
  void fillTrack(T1 const& track, T2& trackProducts, T3& collisionProducts, T4& indexMap)
  {
    if (mProduceTracks) {
      trackProducts.producedTracks(collisionProducts.producedCollision.lastIndex(),
                                   track.pt() * track.sign(),
                                   track.eta(),
                                   track.phi());
      indexMap.emplace(track.globalIndex(), trackProducts.producedTracks.lastIndex());
    }
    if (mProduceTrackMasks) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedTrackMasks(mTrackSelection.getBitmask());
      } else {
        trackProducts.producedTrackMasks(static_cast<o2::aod::femtodatatypes::TrackMaskType>(0u));
      }
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
                                        track.beta(),
                                        track.mass());
    }
    if (mProduceElectronPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedElectronPids(track.itsNSigmaEl(), track.tpcNSigmaEl(), track.tofNSigmaEl());
      } else {
        trackProducts.producedElectronPids(0, track.tpcNSigmaEl(), track.tofNSigmaEl());
      }
    }
    if (mProducePionPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedPionPids(track.itsNSigmaPi(), track.tpcNSigmaPi(), track.tofNSigmaPi());
      } else {
        trackProducts.producedPionPids(0, track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
    }
    if (mProduceKaonPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedKaonPids(track.itsNSigmaKa(), track.tpcNSigmaKa(), track.tofNSigmaKa());
      } else {
        trackProducts.producedKaonPids(0, track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
    }
    if (mProduceProtonPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedProtonPids(track.itsNSigmaPr(), track.tpcNSigmaPr(), track.tofNSigmaPr());
      } else {
        trackProducts.producedProtonPids(0, track.tpcNSigmaPr(), track.tofNSigmaPr());
      }
    }
    if (mProduceDeuteronPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedDeuteronPids(track.itsNSigmaDe(), track.tpcNSigmaDe(), track.tofNSigmaDe());
      } else {
        trackProducts.producedDeuteronPids(0, track.tpcNSigmaDe(), track.tofNSigmaDe());
      }
    }
    if (mProduceTritonPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedTritonPids(track.itsNSigmaTr(), track.tpcNSigmaTr(), track.tofNSigmaTr());
      } else {
        trackProducts.producedTritonPids(0, track.tpcNSigmaTr(), track.tofNSigmaTr());
      }
    }
    if (mProduceHeliumPids) {
      if constexpr (type == modes::Track::kPrimaryTrack) {
        trackProducts.producedHeliumPids(track.itsNSigmaHe(), track.tpcNSigmaHe(), track.tofNSigmaHe());
      } else {
        trackProducts.producedHeliumPids(0, track.tpcNSigmaHe(), track.tofNSigmaHe());
      }
    }
  }

  template <modes::Track type, typename T1, typename T2, typename T3, typename T4>
  int64_t getDaughterIndex(const T1& daughter, T2& trackProducts, T3& collisionProducts, T4& indexMap)
  {
    auto result = utils::getIndex(daughter.globalIndex(), indexMap);
    if (result) {
      return result.value();
    } else {
      this->fillTrack<type>(daughter, trackProducts, collisionProducts, indexMap);
      int64_t idx = trackProducts.producedTracks.lastIndex();
      indexMap.emplace(daughter.globalIndex(), idx);
      return idx;
    }
  }

 private:
  TrackSelection<HistName> mTrackSelection;
  bool mFillAnyTable = false;
  bool mProduceTracks = false;
  bool mProduceTrackMasks = false;
  bool mProduceTrackDcas = false;
  bool mProduceTrackExtras = false;
  bool mProduceElectronPids = false;
  bool mProducePionPids = false;
  bool mProduceKaonPids = false;
  bool mProduceProtonPids = false;
  bool mProduceDeuteronPids = false;
  bool mProduceTritonPids = false;
  bool mProduceHeliumPids = false;
};

struct TrackBuilderDerivedToDerivedProducts : o2::framework::ProducesGroup {
  o2::framework::Produces<o2::aod::StoredFTracks> producedTracks;
  o2::framework::Produces<o2::aod::StoredFTrackMasks> producedTrackMasks;
};

struct ConfTrackTablesDerivedToDerived : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("TrackTables");
  o2::framework::Configurable<int> limitTrack1{"limitTrack1", 1, "At least this many tracks of type 1 need to be in the collision"};
  o2::framework::Configurable<int> limitTrack2{"limitTrack2", 0, "At least this many tracks of type 2 need to be in the collision"};
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
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  bool collisionHasTooFewTracks(T1& col, T2& /*trackTable*/, T3& partitionTrack1, T4& partitionTrack2, T5& cache)
  {
    auto trackSlice1 = partitionTrack1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto trackSlice2 = partitionTrack2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    if (trackSlice1.size() >= mLimitTrack1 && trackSlice2.size() >= mLimitTrack2) {
      return false;
    }
    return true;
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void processTracks(T1& col, T2& /*trackTable*/, T3& partitionTrack1, T4& partitionTrack2, T5& indexMap, T6& cache, T7& newTrackTable, T8& newCollisionTable)
  {
    auto trackSlice1 = partitionTrack1->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);
    auto trackSlice2 = partitionTrack2->sliceByCached(o2::aod::femtobase::stored::fColId, col.globalIndex(), cache);

    for (auto const& track : trackSlice1) {
      this->fillTrack(track, newTrackTable, newCollisionTable, indexMap);
    }
    for (auto const& track : trackSlice2) {
      this->fillTrack(track, newTrackTable, newCollisionTable, indexMap);
    }
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillTrack(T1 const& track, T2& trackProducts, T3& collisionProducts, T4& indexMap)
  {
    trackProducts.producedTracks(collisionProducts.producedCollision.lastIndex(),
                                 track.signedPt(),
                                 track.eta(),
                                 track.phi());
    trackProducts.producedTrackMasks(track.mask());
    indexMap.emplace(track.globalIndex(), trackProducts.producedTracks.lastIndex());
  }

  template <typename T1, typename T2, typename T3, typename T4>
  int64_t getDaughterIndex(const T1& daughter, T2& trackProducts, T3& collisionProducts, T4& indexMap)
  {
    auto result = utils::getIndex(daughter.globalIndex(), indexMap);
    if (result) {
      return result.value();
    } else {
      this->fillTrack(daughter, trackProducts, collisionProducts, indexMap);
      int64_t idx = trackProducts.producedTracks.lastIndex();
      indexMap.emplace(daughter.globalIndex(), idx);
      return idx;
    }
  }

 private:
  int mLimitTrack1 = 0;
  int mLimitTrack2 = 0;
};

} // namespace trackbuilder
//
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRACKBUILDER_H_
