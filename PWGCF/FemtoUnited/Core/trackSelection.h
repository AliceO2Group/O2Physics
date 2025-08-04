// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file trackSelection.h
/// \brief Definition of track selections
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_

#include "PWGCF/FemtoUnited/Core/baseSelection.h"
#include "PWGCF/FemtoUnited/Core/dataTypes.h"

#include "Framework/Configurable.h"

#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>

namespace o2::analysis::femtounited
{
namespace trackselection
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
  o2::framework::Configurable<std::vector<float>> tpcSharedClustersMax{"tpcSharedClustersMax", {160.f}, "Maximum number of shared clusters in TPC"};
  o2::framework::Configurable<std::vector<float>> tpcSharedClusterFractionMax{"tpcSharedClusterFractionMax", {1.f}, "Maximum fraction of shared clusters in TPC"};
  o2::framework::Configurable<std::vector<float>> itsClustersMin{"itsClustersMin", {5.f}, "Minimum number of clusters in ITS"};
  o2::framework::Configurable<std::vector<float>> itsIbClustersMin{"itsIbClustersMin", {3.f}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
  o2::framework::Configurable<std::vector<std::string>> dcaxyMax{"dcaxyMax", {"0.004 + 0.013*TMath::Power(x, -1)"}, "Maximum |dca_xy| as a function of pT. Has to be a valid TForumal, where x=pt"};
  o2::framework::Configurable<std::vector<std::string>> dcazMax{"dcazMax", {"0.004 + 0.013*TMath::Power(x, -1)"}, "Maximum |dca_z| as a function of pT. Has to be a valid TForumal, where x=pt"};

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
  o2::framework::Configurable<std::vector<float>> tpcProton{"tpcProton", {3.f}, "Maximum |nsigma| for proton PID"};
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
  o2::framework::Configurable<int> sign{"sign", 1, "Sign of the track (1 for positive tracks and -1 for negative tracks)"};
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
  kTPCnClsMin,     ///< Min. number of TPC clusters
  kTPCcRowsMin,    ///< Min. number of crossed TPC rows
  kTPCsClsMax,     ///< Max. number of shared TPC clusters
  kTPCsClsFracMax, ///< Max. fractions of shared TPC clusters
  kITSnClsMin,     ///< Min. number of ITS clusters
  kITSnClsIbMin,   ///< Min. number of ITS clusters in the inner barrel
  kDCAxyMax,       ///< Max. |DCA_xy| (cm) as a function of pT
  kDCAzMax,        ///< Max. |DCA_z| (cm) as a function of pT

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
  void configure(T1 const& config, T2 const& filter)
  {
    mMinimalMomentumForTof = config.minMomentumForTof.value;
    // add selections for track quality
    this->addSelection(config.tpcClustersMin.value, trackselection::kTPCnClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.tpcCrossedRowsMin.value, trackselection::kTPCcRowsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.tpcSharedClustersMax.value, trackselection::kTPCsClsMax, limits::kUpperLimit, true, true);
    this->addSelection(config.tpcSharedClusterFractionMax.value, trackselection::kTPCsClsFracMax, limits::kUpperLimit, true, true);
    this->addSelection(config.itsClustersMin.value, trackselection::kITSnClsMin, limits::kLowerLimit, true, true);
    this->addSelection(config.itsIbClustersMin.value, trackselection::kITSnClsIbMin, limits::kLowerLimit, true, true);
    this->addSelection(config.dcaxyMax.name, filter.ptMin.value, filter.ptMax.value, config.dcaxyMax.value, trackselection::kDCAxyMax, limits::kAbsUpperFunctionLimit, true, true);
    this->addSelection(config.dcazMax.name, filter.ptMin.value, filter.ptMax.value, config.dcazMax.value, trackselection::kDCAzMax, limits::kAbsUpperFunctionLimit, true, true);
    // add selections for its pid
    this->addSelection(config.itsElectron.value, trackselection::kItsElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsPion.value, trackselection::kItsPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsKaon.value, trackselection::kItsKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsProton.value, trackselection::kItsProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsDeuteron.value, trackselection::kItsDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsTriton.value, trackselection::kItsTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.itsHelium.value, trackselection::kItsHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpc pid
    this->addSelection(config.tpcElectron.value, trackselection::kTpcElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcPion.value, trackselection::kTpcPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcKaon.value, trackselection::kTpcKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcProton.value, trackselection::kTpcProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcDeuteron.value, trackselection::kTpcDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcTriton.value, trackselection::kTpcTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpcHelium.value, trackselection::kTpcHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tof pid
    this->addSelection(config.tofElectron.value, trackselection::kTofElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofPion.value, trackselection::kTofPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofKaon.value, trackselection::kTofKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofProton.value, trackselection::kTofProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofDeuteron.value, trackselection::kTofDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofTriton.value, trackselection::kTofTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tofHelium.value, trackselection::kTofHelium, limits::kAbsUpperLimit, false, false);
    // add selections for tpctof pid
    this->addSelection(config.tpctofElectron.value, trackselection::kTpctofElectron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofPion.value, trackselection::kTpctofPion, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofKaon.value, trackselection::kTpctofKaon, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofProton.value, trackselection::kTpctofProton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofDeuteron.value, trackselection::kTpctofDeuteron, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofTriton.value, trackselection::kTpctofTriton, limits::kAbsUpperLimit, false, false);
    this->addSelection(config.tpctofHelium.value, trackselection::kTpctofHelium, limits::kAbsUpperLimit, false, false);
  }

  template <typename T>
  void applySelections(T const& Track)
  {

    this->reset();
    this->evaluateObservable(kTPCnClsMin, Track.tpcNClsFound());
    this->evaluateObservable(kTPCcRowsMin, Track.tpcNClsCrossedRows());
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

  template <typename T>
  bool hasTofAboveThreshold(T const& track) const
  {
    // If track momentum exceeds threshold, we require valid TOF info
    return !(track.p() > mMinimalMomentumForTof && !track.hasTOF());
  }

 protected:
  float mMinimalMomentumForTof = 2.f;
};
}; // namespace trackselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKSELECTION_H_
