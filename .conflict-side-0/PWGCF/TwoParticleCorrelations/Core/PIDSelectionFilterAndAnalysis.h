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
#ifndef PIDSELECTIONFILTERANDANALYSIS_H
#define PIDSELECTIONFILTERANDANALYSIS_H

#include "SelectionFilterAndAnalysis.h"
#include "SkimmingConfigurableCuts.h"

#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>

#undef INCORPORATEBAYESIANPID

namespace o2
{
namespace analysis
{
namespace PWGCF
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
    : mPidTpcSel_el{pidtpcel}, mPidTpcSel_mu{pidtpcmu}, mPidTpcSel_pi{pidtpcpi}, mPidTpcSel_ka{pidtpcka}, mPidTpcSel_pr{pidtpcpr}, mPidTofSel_el{pidtofel}, mPidTofSel_mu{pidtofmu}, mPidTofSel_pi{pidtofpi}, mPidTofSel_ka{pidtofka}, mPidTofSel_pr{pidtofpr}, mPidBaySel_el{pidbayel}, mPidBaySel_mu{pidbaymu}, mPidBaySel_pi{pidbaypi}, mPidBaySel_ka{pidbayka}, mPidBaySel_pr{pidbaypr}
  {
  }
  PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                           std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr,
                           std::vector<std::string> pidbayel, std::vector<std::string> pidbaymu, std::vector<std::string> pidbaypi, std::vector<std::string> pidbayka, std::vector<std::string> pidbaypr);
#else
  PIDSelectionConfigurable(std::string pidtpcel = "", std::string pidtpcmu = "", std::string pidtpcpi = "", std::string pidtpcka = "", std::string pidtpcpr = "",
                           std::string pidtofel = "", std::string pidtofmu = "", std::string pidtofpi = "", std::string pidtofka = "", std::string pidtofpr = "")
    : mPidTpcSel_el{pidtpcel},
      mPidTpcSel_mu{pidtpcmu},
      mPidTpcSel_pi{pidtpcpi},
      mPidTpcSel_ka{pidtpcka},
      mPidTpcSel_pr{pidtpcpr},
      mPidTofSel_el{pidtofel},
      mPidTofSel_mu{pidtofmu},
      mPidTofSel_pi{pidtofpi},
      mPidTofSel_ka{pidtofka},
      mPidTofSel_pr{pidtofpr}
  {
  }
  PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                           std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr);
#endif

 private:
  std::string mPidTpcSel_el = "";
  std::string mPidTpcSel_mu = "";
  std::string mPidTpcSel_pi = "";
  std::string mPidTpcSel_ka = "";
  std::string mPidTpcSel_pr = "";
  std::string mPidTofSel_el = "";
  std::string mPidTofSel_mu = "";
  std::string mPidTofSel_pi = "";
  std::string mPidTofSel_ka = "";
  std::string mPidTofSel_pr = "";
  std::string mPidBaySel_el = "";
  std::string mPidBaySel_mu = "";
  std::string mPidBaySel_pi = "";
  std::string mPidBaySel_ka = "";
  std::string mPidBaySel_pr = "";

 private:
  ClassDefNV(PIDSelectionConfigurable, 1);
};

/// \brief Filter of tracks based on PID and track selection once filetered
class PIDSelectionFilterAndAnalysis : public SelectionFilterAndAnalysis
{
 public:
  PIDSelectionFilterAndAnalysis();
  PIDSelectionFilterAndAnalysis(const TString&, selmodes);
  PIDSelectionFilterAndAnalysis(const PIDSelectionConfigurable& pidsel, selmodes mode);
  virtual ~PIDSelectionFilterAndAnalysis() override;

  void SetPTOF(float ptof) { mPTOF = ptof; }
  void SetRequireTOF(bool requiretof = false) { mRequireTOF = requiretof; }
  void SetEllipticTPCTOF(bool elliptic = false) { mEllipticTPCTOF = elliptic; }

  template <typename TrackToFilter>
  uint64_t Filter(TrackToFilter const& track);

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
  void ConstructCutFromString(const TString&);
  virtual int CalculateMaskLength() override;
  virtual void StoreArmedMask() override;

  float mPTOF = 0.8f;           ///< the p threshold for cheking TOF information
  bool mRequireTOF = false;     ///< is TOF required
  bool mEllipticTPCTOF = false; ///< 2D nsigmas elliptic TPC+TOF

  PIDSpecies mSpeciesOfInterest = kWrongSpecies;
  std::vector<CutBrick<float>*> mCloseNsigmasTPC;
  std::vector<CutBrick<float>*> mCloseNsigmasTOF;
  std::vector<CutBrick<float>*> mBayesProbability;

  ClassDef(PIDSelectionFilterAndAnalysis, 1);
};

/// \brief Fills the filter cuts mask
template <typename TrackToFilter>
inline uint64_t PIDSelectionFilterAndAnalysis::Filter(TrackToFilter const& track)
{
  uint64_t selectedMask = 0UL;
  mSelectedMask = 0UL;
  int bit = 0;

  auto filterBrickValue = [&](auto brick, auto value) {
    if (brick != nullptr) {
      std::vector<bool> res = brick->Filter(value);
      for (auto b : res) {
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

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // PIDSELECTIONFILTERANDANALYSIS_H
