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

/// \file PIDSelectionFilterAndAnalysis.h
/// \brief PIDSelectionFilterAndAnalysis class for PID filtering tracks
/// for derived data production
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TWOPARTICLECORRELATIONS_CORE_PIDSELECTIONFILTERANDANALYSIS_H_
#define PWGCF_TWOPARTICLECORRELATIONS_CORE_PIDSELECTIONFILTERANDANALYSIS_H_

#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

#include <TString.h>

#include <Rtypes.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#undef INCORPORATEBAYESIANPID

namespace o2
{
namespace analysis
{
namespace pwgcf
{
/* forward declaration */
class PIDSelectionFilterAndAnalysis;

///\brief Convenience class for configuration access
class PIDSelectionConfigurable
{
  friend class PIDSelectionFilterAndAnalysis;

 public:
#ifdef INCORPORATEBAYESIANPID
  PIDSelectionConfigurable(std::string pidtpcel = "", std::string pidtpcmu = "", std::string pidtpcpi = "", std::string pidtpcka = "", std::string pidtpcpr = "",
                           std::string pidtofel = "", std::string pidtofmu = "", std::string pidtofpi = "", std::string pidtofka = "", std::string pidtofpr = "",
                           std::string pidbayel = "", std::string pidbaymu = "", std::string pidbaypi = "", std::string pidbayka = "", std::string pidbaypr = "")
    : mPidTpcSelEl{pidtpcel}, mPidTpcSelMu{pidtpcmu}, mPidTpcSelPi{pidtpcpi}, mPidTpcSelKa{pidtpcka}, mPidTpcSelPr{pidtpcpr}, mPidTofSelEl{pidtofel}, mPidTofSelMu{pidtofmu}, mPidTofSelPi{pidtofpi}, mPidTofSelKa{pidtofka}, mPidTofSelPr{pidtofpr}, mPidBaySelEl{pidbayel}, mPidBaySelMu{pidbaymu}, mPidBaySelPi{pidbaypi}, mPidBaySelKa{pidbayka}, mPidBaySelPr{pidbaypr}
  {
  }
  PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                           std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr,
                           std::vector<std::string> pidbayel, std::vector<std::string> pidbaymu, std::vector<std::string> pidbaypi, std::vector<std::string> pidbayka, std::vector<std::string> pidbaypr);
#else
  PIDSelectionConfigurable(std::string pidtpcel = "", std::string pidtpcmu = "", std::string pidtpcpi = "", std::string pidtpcka = "", std::string pidtpcpr = "",
                           std::string pidtofel = "", std::string pidtofmu = "", std::string pidtofpi = "", std::string pidtofka = "", std::string pidtofpr = "")
    : mPidTpcSelEl{std::move(pidtpcel)},
      mPidTpcSelMu{std::move(pidtpcmu)},
      mPidTpcSelPi{std::move(pidtpcpi)},
      mPidTpcSelKa{std::move(pidtpcka)},
      mPidTpcSelPr{std::move(pidtpcpr)},
      mPidTofSelEl{std::move(pidtofel)},
      mPidTofSelMu{std::move(pidtofmu)},
      mPidTofSelPi{std::move(pidtofpi)},
      mPidTofSelKa{std::move(pidtofka)},
      mPidTofSelPr{std::move(pidtofpr)}
  {
  }
  PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                           std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr);
#endif

 private:
  std::string mPidTpcSelEl = "";
  std::string mPidTpcSelMu = "";
  std::string mPidTpcSelPi = "";
  std::string mPidTpcSelKa = "";
  std::string mPidTpcSelPr = "";
  std::string mPidTofSelEl = "";
  std::string mPidTofSelMu = "";
  std::string mPidTofSelPi = "";
  std::string mPidTofSelKa = "";
  std::string mPidTofSelPr = "";
  std::string mPidBaySelEl = "";
  std::string mPidBaySelMu = "";
  std::string mPidBaySelPi = "";
  std::string mPidBaySelKa = "";
  std::string mPidBaySelPr = "";

 private:
  ClassDefNV(PIDSelectionConfigurable, 1);
};

/// \brief Filter of tracks based on PID and track selection once filetered
class PIDSelectionFilterAndAnalysis : public SelectionFilterAndAnalysis
{
 public:
  PIDSelectionFilterAndAnalysis();
  PIDSelectionFilterAndAnalysis(const TString&, Selmodes);
  PIDSelectionFilterAndAnalysis(const PIDSelectionConfigurable& pidsel, Selmodes mode);
  ~PIDSelectionFilterAndAnalysis() override;

  void setPTOF(float ptof) { mPTOF = ptof; }
  void setRequireTOF(bool requiretof = false) { mRequireTOF = requiretof; }
  void setEllipticTPCTOF(bool elliptic = false) { mEllipticTPCTOF = elliptic; }

  template <typename TrackToFilter>
  uint64_t filter(TrackToFilter const& track);

  enum PIDSpecies {
    kElectron = 0, ///< electron
    kMuon,         ///< muon
    kPion,         ///< pion
    kKaon,         ///< kaon
    kProton,       ///< proton
    kNoOfSpecies,  ///< the number of considered species
    kWrongSpecies = -1
  };

  enum class PIDCuts : int {
    kITS = 0,
    kTPC,
    kTOF,
    kTPCTOF,
    kNEARSIGMA,
    kAWAYSIGMA,
    kTOFREQ,
    kNCuts
  };

  static const std::string mCutNames[static_cast<int>(PIDCuts::kNCuts)];
  static const std::vector<std::string> mSpeciesNames;

 private:
  void constructCutFromString(const TString&);
  int calculateMaskLength() override;
  void storeArmedMask() override;

  float mPTOF = 0.8f;           ///< the p threshold for cheking TOF information
  bool mRequireTOF = false;     ///< is TOF required
  bool mEllipticTPCTOF = false; ///< 2D nsigmas elliptic TPC+TOF

  PIDSpecies mSpeciesOfInterest = kWrongSpecies;
  std::vector<CutBrick<float>*> mCloseNsigmasTPC;
  std::vector<CutBrick<float>*> mCloseNsigmasTOF;
  std::vector<CutBrick<float>*> mBayesProbability;

  ClassDefOverride(PIDSelectionFilterAndAnalysis, 1);
};

/// \brief Fills the filter cuts mask
template <typename TrackToFilter>
inline uint64_t PIDSelectionFilterAndAnalysis::filter(TrackToFilter const& track)
{
  uint64_t selectedMask = 0UL;
  mSelectedMask = 0UL;
  int bit = 0;

  auto filterBrickValue = [&](auto brick, auto value) {
    if (brick != nullptr) {
      std::vector<bool> res = brick->filter(value);
      for (const auto& b : res) {
        if (b) {
          SETBIT(selectedMask, bit);
        }
        bit++;
      }
    }
  };
  filterBrickValue(mCloseNsigmasTPC[kElectron], track.tpcNSigmaEl());
  filterBrickValue(mCloseNsigmasTPC[kMuon], track.tpcNSigmaMu());
  filterBrickValue(mCloseNsigmasTPC[kPion], track.tpcNSigmaPi());
  filterBrickValue(mCloseNsigmasTPC[kKaon], track.tpcNSigmaKa());
  filterBrickValue(mCloseNsigmasTPC[kProton], track.tpcNSigmaPr());
  filterBrickValue(mCloseNsigmasTOF[kElectron], track.tofNSigmaEl());
  filterBrickValue(mCloseNsigmasTOF[kMuon], track.tofNSigmaMu());
  filterBrickValue(mCloseNsigmasTOF[kPion], track.tofNSigmaPi());
  filterBrickValue(mCloseNsigmasTOF[kKaon], track.tofNSigmaKa());
  filterBrickValue(mCloseNsigmasTOF[kProton], track.tofNSigmaPr());
#ifdef INCORPORATEBAYESIANPID
  filterBrickValue(mBayesProbability[kElectron], track.bayesEl());
  filterBrickValue(mBayesProbability[kMuon], track.bayesMu());
  filterBrickValue(mBayesProbability[kPion], track.bayesPi());
  filterBrickValue(mBayesProbability[kKaon], track.bayesKa());
  filterBrickValue(mBayesProbability[kProton], track.bayesPr());
#endif

  mSelectedMask = selectedMask;
  return mSelectedMask;
}

} // namespace pwgcf
} // namespace analysis
} // namespace o2

#endif // PWGCF_TWOPARTICLECORRELATIONS_CORE_PIDSELECTIONFILTERANDANALYSIS_H_
