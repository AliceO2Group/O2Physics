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

  o2::framework::Configurable<float> minMomentumForTof{"minMomentumForTof", 2.0f, "Minimum momentum to required TOF PID (all species)"};

  // track its pid cuts
  o2::framework::Configurable<std::vector<float>> itsElectron{"itsElectron", {}, "Maximum |nsigma| for electron PID"};
  o2::framework::Configurable<std::vector<float>> itsPion{"itsPion", {}, "Maximum |nsigma| for pion PID"};
  o2::framework::Configurable<std::vector<float>> itsKaon{"itsKaon", {}, "Maximum |nsigma| for kaon PID"};
  o2::framework::Configurable<std::vector<float>> itsProton{"itsProton", {}, "Maximum |nsigma| for proton PID"};
  o2::framework::Configurable<std::vector<float>> itsDeuteron{"itsDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
  o2::framework::Configurable<std::vector<float>> itsTriton{"itsTriton", {}, "Maximum |nsigma| for trition PID"};
  o2::framework::Configurable<std::vector<float>> itsHelium{"itsHelium", {}, "Maximum |nsigma| for helium PID"};

  // track tpc pid cuts
  o2::framework::Configurable<std::vector<float>> tpcElectron{"tpcElectron", {}, "Maximum |nsigma| for electron PID"};
  o2::framework::Configurable<std::vector<float>> tpcPion{"tpcPion", {}, "Maximum |nsigma| for pion PID"};
  o2::framework::Configurable<std::vector<float>> tpcKaon{"tpcKaon", {}, "Maximum |nsigma| for kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpcProton{"tpcProton", {}, "Maximum |nsigma| for proton PID"};
  o2::framework::Configurable<std::vector<float>> tpcDeuteron{"tpcDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpcTriton{"tpcTriton", {}, "Maximum |nsigma| for trition PID"};
  o2::framework::Configurable<std::vector<float>> tpcHelium{"tpcHelium", {}, "Maximum |nsigma| for helium PID"};

  // track tof pid cuts
  o2::framework::Configurable<std::vector<float>> tofElectron{"tofElectron", {}, "Maximum |nsigma| for electron PID"};
  o2::framework::Configurable<std::vector<float>> tofPion{"tofPion", {}, "Maximum |nsigma| for pion PID"};
  o2::framework::Configurable<std::vector<float>> tofKaon{"tofKaon", {}, "Maximum |nsigma| for kaon PID"};
  o2::framework::Configurable<std::vector<float>> tofProton{"tofProton", {}, "Maximum |nsigma| for proton PID"};
  o2::framework::Configurable<std::vector<float>> tofDeuteron{"tofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tofTriton{"tofTriton", {}, "Maximum |nsigma| for trition PID"};
  o2::framework::Configurable<std::vector<float>> tofHelium{"tofHelium", {}, "Maximum |nsigma| for helium PID"};

  // track tpcits pid cuts
  o2::framework::Configurable<std::vector<float>> tpcitsElectron{"tpcitsElectron", {}, "Maximum |nsigma| for electron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsPion{"tpcitsPion", {}, "Maximum |nsigma| for pion PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsKaon{"tpcitsKaon", {}, "Maximum |nsigma| for kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsProton{"tpcitsProton", {3.f}, "Maximum |nsigma| for proton PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsDeuteron{"tpcitsDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsTriton{"tpcitsTriton", {}, "Maximum |nsigma| for trition PID"};
  o2::framework::Configurable<std::vector<float>> tpcitsHelium{"tpcitsHelium", {}, "Maximum |nsigma| for helium PID"};

  // track tpctof pid cuts
  o2::framework::Configurable<std::vector<float>> tpctofElectron{"tpctofElectron", {}, "Maximum |nsigma| for electron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofPion{"tpctofPion", {}, "Maximum |nsigma| for pion PID"};
  o2::framework::Configurable<std::vector<float>> tpctofKaon{"tpctofKaon", {}, "Maximum |nsigma| for kaon PID"};
  o2::framework::Configurable<std::vector<float>> tpctofProton{"tpctofProton", {3.f}, "Maximum |nsigma| for proton PID"};
  o2::framework::Configurable<std::vector<float>> tpctofDeuteron{"tpctofDeuteron", {}, "Maximum |nsigma| for deuteron PID"};
  o2::framework::Configurable<std::vector<float>> tpctofTriton{"tpctofTriton", {}, "Maximum |nsigma| for trition PID"};
  o2::framework::Configurable<std::vector<float>> tpctofHelium{"tpctofHelium", {}, "Maximum |nsigma| for helium PID"};
};

// define the template structure for TrackSelection
template <const char* Prefix>
struct ConfTrackSelection : public o2::framework::ConfigurableGroup {
  std::string prefix = Prefix; // Unique prefix based on the template argument
  // configuration parameters
  o2::framework::Configurable<int> pdgCode{"pdgCode", 2212, "Track PDG code"};
  o2::framework::Configurable<int> charge{"charge", 1, "Charge of the track (use +/-1 for positive/negative tracks, except He3 needs +/-2)"};
  // filters for kinematics
  o2::framework::Configurable<float> ptMin{"ptMin", 0.2f, "Minimum pT (GeV/c)"};
  o2::framework::Configurable<float> ptMax{"ptMax", 6.f, "Maximum pT (GeV/c)"};
  o2::framework::Configurable<float> etaMin{"etaMin", -0.9f, "Minimum eta"};
  o2::framework::Configurable<float> etaMax{"etaMax", 0.9f, "Maximum eta"};
  o2::framework::Configurable<float> phiMin{"phiMin", 0.f, "Minimum phi"};
  o2::framework::Configurable<float> phiMax{"phiMax", 1.f * o2::constants::math::TwoPI, "Maximum phi"};
  // track selection masks
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskLowMomentum{"maskLowMomentum", 2u, "Bitmask for selections below momentum threshold"};
  o2::framework::Configurable<o2::aod::femtodatatypes::TrackMaskType> maskHighMomentum{"maskHighMomentum", 1u, "Bitmask for selections above momentum threshold"};
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

const char trackSelsName[] = "Track Selection Object";
const std::unordered_map<TrackSels, std::string> trackSelsToString = {
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
class TrackSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackMaskType, kTrackSelsMax>
{
 public:
  TrackSelection() {}
  virtual ~TrackSelection() = default;

  template <typename T1, typename T2>
  void configure(T1& config, T2& filter)
  {
    mPtMin = filter.ptMin;
    mPtMax = filter.ptMax;
    mEtaMin = filter.etaMin;
    mEtaMax = filter.etaMax;
    mPhiMin = filter.phiMin;
    mPhiMax = filter.phiMax;
    mMinimalMomentumForTof = config.minMomentumForTof.value;

    // add selections for track quality
    this->addSelection(config.tpcClustersMin.value, kTPCnClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.tpcCrossedRowsMin.value, kTPCcRowsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.tpcClustersOverCrossedRows.value, kTPCnClsOvercRowsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.tpcSharedClustersMax.value, kTPCsClsMax, limits::kUpperLimit, true, true);
    this->addSelection(config.tpcSharedClusterFractionMax.value, kTPCsClsFracMax, limits::kUpperLimit, true, true);
    this->addSelection(config.itsClustersMin.value, kITSnClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.itsIbClustersMin.value, kITSnClsIbMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dcaxyMax.name, filter.ptMin.value, filter.ptMax.value, config.dcaxyMax.value, kDCAxyMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.dcazMax.name, filter.ptMin.value, filter.ptMax.value, config.dcazMax.value, kDCAzMax, limits::kAbsUpperFunctionLimit, true, true);

    // add selections for its pid
    this->addSelection(config.itsElectron.value, kItsElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsPion.value, kItsPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsKaon.value, kItsKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsProton.value, kItsProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsDeuteron.value, kItsDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsTriton.value, kItsTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsHelium.value, kItsHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpc pid
    this->addSelection(config.tpcElectron.value, kTpcElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcPion.value, kTpcPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcKaon.value, kTpcKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcProton.value, kTpcProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcDeuteron.value, kTpcDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcTriton.value, kTpcTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcHelium.value, kTpcHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tof pid
    this->addSelection(config.tofElectron.value, kTofElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofPion.value, kTofPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofKaon.value, kTofKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofProton.value, kTofProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofDeuteron.value, kTofDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofTriton.value, kTofTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofHelium.value, kTofHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpcits pid
    this->addSelection(config.tpcitsElectron.value, kTpcitsElectron, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsPion.value, kTpcitsPion, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsKaon.value, kTpcitsKaon, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsProton.value, kTpcitsProton, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsDeuteron.value, kTpcitsDeuteron, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsTriton.value, kTpcitsTriton, limits::kUpperLimit, false, false);
    this->addSelection(config.tpcitsHelium.value, kTpcitsHelium, limits::kUpperLimit, false, false);
    // add selections for tpctof pid
    this->addSelection(config.tpctofElectron.value, kTpctofElectron, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofPion.value, kTpctofPion, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofKaon.value, kTpctofKaon, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofProton.value, kTpctofProton, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofDeuteron.value, kTpctofDeuteron, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofTriton.value, kTpctofTriton, limits::kUpperLimit, false, false);
    this->addSelection(config.tpctofHelium.value, kTpctofHelium, limits::kUpperLimit, false, false);
  }

  template <typename T>
  bool hasTofAboveThreshold(T const& track) const
  {
    // If track momentum exceeds threshold, we require valid TOF info
    return !(track.p() > mMinimalMomentumForTof && !track.hasTOF());
  }

  template <typename T>
  bool checkFilters(T const& track) const
  {
    return ((track.pt() > mPtMin && track.pt() < mPtMax) &&
            (track.eta() > mEtaMin && track.eta() < mEtaMax) &&
            (track.phi() > mPhiMin && track.phi() < mPhiMax));
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

    // its pid
    this->evaluateObservable(kItsElectron, Track.itsNSigmaEl());
    this->evaluateObservable(kItsPion, Track.itsNSigmaPi());
    this->evaluateObservable(kItsKaon, Track.itsNSigmaKa());
    this->evaluateObservable(kItsProton, Track.itsNSigmaPr());
    this->evaluateObservable(kItsDeuteron, Track.itsNSigmaDe());
    this->evaluateObservable(kItsTriton, Track.itsNSigmaTr());
    this->evaluateObservable(kItsHelium, Track.itsNSigmaHe());

    // tpc pid
    this->evaluateObservable(kTpcElectron, Track.tpcNSigmaEl());
    this->evaluateObservable(kTpcPion, Track.tpcNSigmaPi());
    this->evaluateObservable(kTpcKaon, Track.tpcNSigmaKa());
    this->evaluateObservable(kTpcProton, Track.tpcNSigmaPr());
    this->evaluateObservable(kTpcDeuteron, Track.tpcNSigmaDe());
    this->evaluateObservable(kTpctofTriton, Track.tpcNSigmaTr());
    this->evaluateObservable(kTpcHelium, Track.tpcNSigmaHe());

    // tof pid
    this->evaluateObservable(kTofElectron, Track.tofNSigmaEl());
    this->evaluateObservable(kTofPion, Track.tofNSigmaPi());
    this->evaluateObservable(kTofKaon, Track.tofNSigmaKa());
    this->evaluateObservable(kTofProton, Track.tofNSigmaPr());
    this->evaluateObservable(kTofDeuteron, Track.tofNSigmaDe());
    this->evaluateObservable(kTofTriton, Track.tofNSigmaTr());
    this->evaluateObservable(kTofHelium, Track.tofNSigmaHe());

    // combined tpc + its pid
    this->evaluateObservable(kTpcitsElectron, std::hypot(Track.tpcNSigmaEl(), Track.itsNSigmaEl()));
    this->evaluateObservable(kTpcitsPion, std::hypot(Track.tpcNSigmaPi(), Track.itsNSigmaPi()));
    this->evaluateObservable(kTpcitsKaon, std::hypot(Track.tpcNSigmaKa(), Track.itsNSigmaKa()));
    this->evaluateObservable(kTpcitsProton, std::hypot(Track.tpcNSigmaPr(), Track.itsNSigmaPr()));
    this->evaluateObservable(kTpcitsDeuteron, std::hypot(Track.tpcNSigmaDe(), Track.itsNSigmaDe()));
    this->evaluateObservable(kTpcitsTriton, std::hypot(Track.tpcNSigmaTr(), Track.itsNSigmaTr()));
    this->evaluateObservable(kTpcitsHelium, std::hypot(Track.tpcNSigmaHe(), Track.itsNSigmaHe()));

    // combined tpc + tof pid
    this->evaluateObservable(kTpctofElectron, std::hypot(Track.tpcNSigmaEl(), Track.tofNSigmaEl()));
    this->evaluateObservable(kTpctofPion, std::hypot(Track.tpcNSigmaPi(), Track.tofNSigmaPi()));
    this->evaluateObservable(kTpctofKaon, std::hypot(Track.tpcNSigmaKa(), Track.tofNSigmaKa()));
    this->evaluateObservable(kTpctofProton, std::hypot(Track.tpcNSigmaPr(), Track.tofNSigmaPr()));
    this->evaluateObservable(kTpctofDeuteron, std::hypot(Track.tpcNSigmaDe(), Track.tofNSigmaDe()));
    this->evaluateObservable(kTpctofTriton, std::hypot(Track.tpcNSigmaTr(), Track.tofNSigmaTr()));
    this->evaluateObservable(kTpctofHelium, std::hypot(Track.tpcNSigmaHe(), Track.tofNSigmaHe()));

    this->assembleBitmask();
  };

 protected:
  float mMinimalMomentumForTof = 2.f;
  float mPtMin = 0.f;
  float mPtMax = 0.f;
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

class TrackBuilder
{
 public:
  TrackBuilder() = default;
  virtual ~TrackBuilder() = default;

  template <typename T1, typename T2, typename T3, typename T4>
  void init(T1& config, T2& filter, T3& table, T4& initContext)
  {
    mTrackSelection.configure(config, filter);
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
      mTrackSelection.printSelections(trackSelsName, trackSelsToString);
    } else {
      LOG(info) << "No tables configured";
    }
    LOG(info) << "Initialization done...";
  }

  template <typename T1, typename T2, typename T3, typename T4>
  void fillTracks(T1 const& tracks, T2& trackProducts, T3& collisionProducts, T4& indexMap)
  {
    if (!mFillAnyTable) {
      return;
    }
    for (const auto& track : tracks) {
      if (!mTrackSelection.checkFilters(track) || !mTrackSelection.hasTofAboveThreshold(track)) {
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
    indexMap.emplace(track.globalIndex(), trackProducts.producedTracks.lastIndex());
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
  TrackSelection mTrackSelection;
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

} // namespace trackbuilder
} // namespace o2::analysis::femto

#endif // PWGCF_FEMTO_CORE_TRACKBUILDER_H_
