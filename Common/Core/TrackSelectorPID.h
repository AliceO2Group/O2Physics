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

/// \file TrackSelectorPID.h
/// \brief PID track selector class
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#ifndef COMMON_CORE_TRACKSELECTORPID_H_
#define COMMON_CORE_TRACKSELECTORPID_H_

#include <TPDGCode.h>

#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"

/// Class for track selection using PID detectors

struct TrackSelectorPID {
  /// Selection status
  enum Status {
    NotApplicable = 0,
    Rejected,
    Conditional,
    Accepted
  };
};

template <uint64_t pdg = kPiPlus>
class TrackSelectorPidBase
{
 public:
  /// Default constructor
  TrackSelectorPidBase() = default;

  /// Conversion operator
  template <uint64_t pdgNew>
  operator TrackSelectorPidBase<pdgNew>() const
  {
    TrackSelectorPidBase<pdgNew> objNew;
    // TPC
    objNew.setRangePtTpc(mPtTpcMin, mPtTpcMax);
    objNew.setRangeNSigmaTpc(mNSigmaTpcMin, mNSigmaTpcMax);
    objNew.setRangeNSigmaTpcCondTof(mNSigmaTpcMinCondTof, mNSigmaTpcMaxCondTof);
    // TOF
    objNew.setRangePtTof(mPtTofMin, mPtTofMax);
    objNew.setRangeNSigmaTof(mNSigmaTofMin, mNSigmaTofMax);
    objNew.setRangeNSigmaTofCondTpc(mNSigmaTofMinCondTpc, mNSigmaTofMaxCondTpc);
    // RICH
    objNew.setRangePtRich(mPtRichMin, mPtRichMax);
    objNew.setRangeNSigmaRich(mNSigmaRichMin, mNSigmaRichMax);
    objNew.setRangeNSigmaRichCondTof(mNSigmaRichMinCondTof, mNSigmaRichMaxCondTof);
    // Bayesian
    objNew.setRangePtBayes(mPtBayesMin, mPtBayesMax);
    objNew.setProbBayesMin(mProbBayesMin);
    return objNew;
  }

  /// Default destructor
  ~TrackSelectorPidBase() = default;

  // TPC

  /// Set pT range where TPC PID is applicable.
  void setRangePtTpc(float ptMin, float ptMax)
  {
    mPtTpcMin = ptMin;
    mPtTpcMax = ptMax;
  }

  /// Set TPC nσ range in which a track should be accepted.
  void setRangeNSigmaTpc(float nsMin, float nsMax)
  {
    mNSigmaTpcMin = nsMin;
    mNSigmaTpcMax = nsMax;
  }

  /// Set the custom value of TPC nσ vs. pion hypothesis to be used for selections
  void setCustomNSigmaTpcPi(float nSigma)
  {
    mUseCustomNSigmaTpcPi = true;
    mCustomNSigmaTpcPi = nSigma;
  }
  /// Set the custom value of TPC nσ vs. kaon hypothesis to be used for selections
  void setCustomNSigmaTpcKa(float nSigma)
  {
    mUseCustomNSigmaTpcKa = true;
    mCustomNSigmaTpcKa = nSigma;
  }
  /// Set the custom value of TPC nσ vs. proton hypothesis to be used for selections
  void setCustomNSigmaTpcPr(float nSigma)
  {
    mUseCustomNSigmaTpcPr = true;
    mCustomNSigmaTpcPr = nSigma;
  }
  /// Set the custom value of TPC nσ vs. muon hypothesis to be used for selections
  void setCustomNSigmaTpcMu(float nSigma)
  {
    mUseCustomNSigmaTpcMu = true;
    mCustomNSigmaTpcMu = nSigma;
  }
  /// Set the custom value of TPC nσ vs. electron hypothesis to be used for selections
  void setCustomNSigmaTpcEl(float nSigma)
  {
    mUseCustomNSigmaTpcEl = true;
    mCustomNSigmaTpcEl = nSigma;
  }

  /// Set TPC nσ range in which a track should be conditionally accepted if combined with TOF. Set to 0 to disable.
  void setRangeNSigmaTpcCondTof(float nsMin, float nsMax)
  {
    mNSigmaTpcMinCondTof = nsMin;
    mNSigmaTpcMaxCondTof = nsMax;
  }

  /// Checks if track is OK for TPC PID.
  /// \param track  track
  /// \return true if track is OK for TPC PID
  template <typename T>
  bool isValidForTpc(const T& track)
  {
    auto pt = track.pt();
    return mPtTpcMin <= pt && pt <= mPtTpcMax;
  }

  /// Checks if track is compatible with given particle species hypothesis within given TPC nσ range.
  /// \param track  track
  /// \param conditionalTof  variable to store the result of selection with looser cuts for conditional accepting of track if combined with TOF
  /// \return true if track satisfies TPC PID hypothesis for given TPC nσ range
  template <typename T>
  bool isSelectedByTpc(const T& track, bool& conditionalTof)
  {
    // Accept if selection is disabled via large values.
    if (mNSigmaTpcMin < -999. && mNSigmaTpcMax > 999.) {
      return true;
    }

    // Get nσ for a given particle hypothesis.
    double nSigma = 100.;
    if constexpr (pdg == kElectron) {
      nSigma = mUseCustomNSigmaTpcEl ? mCustomNSigmaTpcEl : track.tpcNSigmaEl();
    } else if constexpr (pdg == kMuonMinus) {
      nSigma = mUseCustomNSigmaTpcMu ? mCustomNSigmaTpcMu : track.tpcNSigmaMu();
    } else if constexpr (pdg == kPiPlus) {
      nSigma = mUseCustomNSigmaTpcPi ? mCustomNSigmaTpcPi : track.tpcNSigmaPi();
    } else if constexpr (pdg == kKPlus) {
      nSigma = mUseCustomNSigmaTpcKa ? mCustomNSigmaTpcKa : track.tpcNSigmaKa();
    } else if constexpr (pdg == kProton) {
      nSigma = mUseCustomNSigmaTpcPr ? mCustomNSigmaTpcPr : track.tpcNSigmaPr();
    } else {
      errorPdg();
    }

    if (mNSigmaTpcMinCondTof < -999. && mNSigmaTpcMaxCondTof > 999.) {
      conditionalTof = true;
    } else {
      conditionalTof = mNSigmaTpcMinCondTof <= nSigma && nSigma <= mNSigmaTpcMaxCondTof;
    }
    return mNSigmaTpcMin <= nSigma && nSigma <= mNSigmaTpcMax;
  }

  /// Returns status of TPC PID selection for a given track.
  /// \param track  track
  /// \return TPC selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusTpc(const T& track)
  {
    if (!isValidForTpc(track)) {
      return TrackSelectorPID::NotApplicable;
    }
    bool condTof = false;
    if (isSelectedByTpc(track, condTof)) {
      return TrackSelectorPID::Accepted;
    } else if (condTof) {
      return TrackSelectorPID::Conditional; // potential to be accepted if combined with TOF
    } else {
      return TrackSelectorPID::Rejected;
    }
  }

  // TOF

  /// Set pT range where TOF PID is applicable.
  void setRangePtTof(float ptMin, float ptMax)
  {
    mPtTofMin = ptMin;
    mPtTofMax = ptMax;
  }

  /// Set TOF nσ range in which a track should be accepted.
  void setRangeNSigmaTof(float nsMin, float nsMax)
  {
    mNSigmaTofMin = nsMin;
    mNSigmaTofMax = nsMax;
  }

  /// Set the custom value of TOF nσ vs. pion hypothesis to be used for selections
  void setCustomNSigmaTofPi(float nSigma)
  {
    mUseCustomNSigmaTofPi = true;
    mCustomNSigmaTofPi = nSigma;
  }
  /// Set the custom value of TOF nσ vs. kaon hypothesis to be used for selections
  void setCustomNSigmaTofKa(float nSigma)
  {
    mUseCustomNSigmaTofKa = true;
    mCustomNSigmaTofKa = nSigma;
  }
  /// Set the custom value of TOF nσ vs. proton hypothesis to be used for selections
  void setCustomNSigmaTofPr(float nSigma)
  {
    mUseCustomNSigmaTofPr = true;
    mCustomNSigmaTofPr = nSigma;
  }
  /// Set the custom value of TOF nσ vs. muon hypothesis to be used for selections
  void setCustomNSigmaTofMu(float nSigma)
  {
    mUseCustomNSigmaTofMu = true;
    mCustomNSigmaTofMu = nSigma;
  }
  /// Set the custom value of TOF nσ vs. electron hypothesis to be used for selections
  void setCustomNSigmaTofEl(float nSigma)
  {
    mUseCustomNSigmaTofEl = true;
    mCustomNSigmaTofEl = nSigma;
  }

  /// Set TOF nσ range in which a track should be conditionally accepted if combined with TPC. Set to 0 to disable.
  void setRangeNSigmaTofCondTpc(float nsMin, float nsMax)
  {
    mNSigmaTofMinCondTpc = nsMin;
    mNSigmaTofMaxCondTpc = nsMax;
  }

  /// Checks if track is OK for TOF PID.
  /// \param track  track
  /// \return true if track is OK for TOF PID
  template <typename T>
  bool isValidForTof(const T& track)
  {
    auto pt = track.pt();
    return mPtTofMin <= pt && pt <= mPtTofMax;
  }

  /// Checks if track is compatible with given particle species hypothesis within given TOF nσ range.
  /// \param track  track
  /// \param conditionalTpc  variable to store the result of selection with looser cuts for conditional accepting of track if combined with TPC
  /// \return true if track satisfies TOF PID hypothesis for given TOF nσ range
  template <typename T>
  bool isSelectedByTof(const T& track, bool& conditionalTpc)
  {
    // Accept if selection is disabled via large values.
    if (mNSigmaTofMin < -999. && mNSigmaTofMax > 999.) {
      return true;
    }

    // Get nσ for a given particle hypothesis.
    double nSigma = 100.;
    if constexpr (pdg == kElectron) {
      nSigma = mUseCustomNSigmaTofEl ? mCustomNSigmaTofEl : track.tofNSigmaEl();
    } else if constexpr (pdg == kMuonMinus) {
      nSigma = mUseCustomNSigmaTofMu ? mCustomNSigmaTofMu : track.tofNSigmaMu();
    } else if constexpr (pdg == kPiPlus) {
      nSigma = mUseCustomNSigmaTofPi ? mCustomNSigmaTofPi : track.tofNSigmaPi();
    } else if constexpr (pdg == kKPlus) {
      nSigma = mUseCustomNSigmaTofKa ? mCustomNSigmaTofKa : track.tofNSigmaKa();
    } else if constexpr (pdg == kProton) {
      nSigma = mUseCustomNSigmaTofPr ? mCustomNSigmaTofPr : track.tofNSigmaPr();
    } else {
      errorPdg();
    }

    if (mNSigmaTofMinCondTpc < -999. && mNSigmaTofMaxCondTpc > 999.) {
      conditionalTpc = true;
    } else {
      conditionalTpc = mNSigmaTofMinCondTpc <= nSigma && nSigma <= mNSigmaTofMaxCondTpc;
    }
    return mNSigmaTofMin <= nSigma && nSigma <= mNSigmaTofMax;
  }

  /// Returns status of TOF PID selection for a given track.
  /// \param track  track
  /// \return TOF selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusTof(const T& track)
  {
    if (!isValidForTof(track)) {
      return TrackSelectorPID::NotApplicable;
    }
    bool condTpc = false;
    if (isSelectedByTof(track, condTpc)) {
      return TrackSelectorPID::Accepted;
    } else if (condTpc) {
      return TrackSelectorPID::Conditional; // potential to be accepted if combined with TPC
    } else {
      return TrackSelectorPID::Rejected;
    }
  }

  // RICH

  /// Set pT range where RICH PID is applicable.
  void setRangePtRich(float ptMin, float ptMax)
  {
    mPtRichMin = ptMin;
    mPtRichMax = ptMax;
  }

  /// Set RICH nσ range in which a track should be accepted.
  void setRangeNSigmaRich(float nsMin, float nsMax)
  {
    mNSigmaRichMin = nsMin;
    mNSigmaRichMax = nsMax;
  }

  /// Set RICH nσ range in which a track should be conditionally accepted if combined with TOF.
  void setRangeNSigmaRichCondTof(float nsMin, float nsMax)
  {
    mNSigmaRichMinCondTof = nsMin;
    mNSigmaRichMaxCondTof = nsMax;
  }

  /// Checks if track is OK for RICH PID.
  /// \param track  track
  /// \return true if track is OK for RICH PID
  template <typename T>
  bool isValidForRich(const T& track)
  {
    if (track.richId() < 0) {
      return false;
    }
    auto pt = track.pt();
    return mPtRichMin <= pt && pt <= mPtRichMax;
  }

  /// Checks if track is compatible with given particle species hypothesis within given RICH nσ range.
  /// \param track  track
  /// \param conditionalTof  variable to store the result of selection with looser cuts for conditional accepting of track if combined with TOF
  /// \return true if track satisfies RICH PID hypothesis for given RICH nσ range
  template <typename T>
  bool isSelectedByRich(const T& track, bool& conditionalTof)
  {
    // Accept if selection is disabled via large values.
    if (mNSigmaRichMin < -999. && mNSigmaRichMax > 999.) {
      return true;
    }

    // Get nσ for a given particle hypothesis.
    double nSigma = 100.;
    if constexpr (pdg == kElectron) {
      nSigma = track.rich().richNsigmaEl();
    } else if constexpr (pdg == kMuonMinus) {
      nSigma = track.rich().richNsigmaMu();
    } else if constexpr (pdg == kPiPlus) {
      nSigma = track.rich().richNsigmaPi();
    } else if constexpr (pdg == kKPlus) {
      nSigma = track.rich().richNsigmaKa();
    } else if constexpr (pdg == kProton) {
      nSigma = track.rich().richNsigmaPr();
    } else {
      errorPdg();
    }

    if (mNSigmaRichMinCondTof < -999. && mNSigmaRichMaxCondTof > 999.) {
      conditionalTof = true;
    } else {
      conditionalTof = mNSigmaRichMinCondTof <= nSigma && nSigma <= mNSigmaRichMaxCondTof;
    }
    return mNSigmaRichMin <= nSigma && nSigma <= mNSigmaRichMax;
  }

  /// Returns status of RICH PID selection for a given track.
  /// \param track  track
  /// \return RICH selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusRich(const T& track)
  {
    if (!isValidForRich(track)) {
      return TrackSelectorPID::NotApplicable;
    }
    bool condTof = false;
    if (isSelectedByRich(track, condTof)) {
      return TrackSelectorPID::Accepted;
    } else if (condTof) {
      return TrackSelectorPID::Conditional; // potential to be accepted if combined with TOF
    } else {
      return TrackSelectorPID::Rejected;
    }
  }

  // MID

  /// Checks if track is OK for MID PID.
  /// \param track  track
  /// \return true if track is OK for MID PID
  template <typename T>
  bool isValidForMid(const T& track)
  {
    if constexpr (pdg == kMuonMinus) {
      return track.midId() > -1;
    } else {
      errorPdg();
      return false;
    }
  }

  /// Checks if track is compatible with muon hypothesis in the MID detector.
  /// \param track  track
  /// \return true if track has been identified as muon by the MID detector
  template <typename T>
  bool isSelectedByMid(const T& track)
  {
    if constexpr (pdg == kMuonMinus) {
      return track.mid().midIsMuon() == 1; // FIXME: change to return track.midIsMuon() once the column is bool.
    } else {
      errorPdg();
      return false;
    }
  }

  /// Returns status of MID PID selection for a given track.
  /// \param track  track
  /// \return MID selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusMid(const T& track)
  {
    if constexpr (pdg == kMuonMinus) {
      if (!isValidForMid(track)) {
        return TrackSelectorPID::NotApplicable;
      }
      if (isSelectedByMid(track)) {
        return TrackSelectorPID::Accepted;
      } else {
        return TrackSelectorPID::Rejected;
      }
    } else {
      errorPdg();
      return TrackSelectorPID::Rejected;
    }
  }

  // Combined selection (TPC + TOF)

  /// Returns status of combined PID (TPC or TOF) selection for a given track.
  /// \param track  track
  /// \return status of combined PID (TPC or TOF) (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusTpcOrTof(const T& track)
  {
    int pidTpc = statusTpc(track);
    int pidTof = statusTof(track);

    if (pidTpc == TrackSelectorPID::Accepted || pidTof == TrackSelectorPID::Accepted) {
      return TrackSelectorPID::Accepted;
    }
    if (pidTpc == TrackSelectorPID::Conditional && pidTof == TrackSelectorPID::Conditional) {
      return TrackSelectorPID::Accepted;
    }
    if (pidTpc == TrackSelectorPID::Rejected || pidTof == TrackSelectorPID::Rejected) {
      return TrackSelectorPID::Rejected;
    }
    return TrackSelectorPID::NotApplicable; // (NotApplicable for one detector) and (NotApplicable or Conditional for the other)
  }

  /// Returns status of combined PID (TPC and TOF) selection for a given track when both detectors are applicable. Returns status of single PID otherwise.
  /// \param track  track
  /// \return status of combined PID (TPC and TOF) (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusTpcAndTof(const T& track)
  {
    int pidTpc = TrackSelectorPID::NotApplicable;
    if (track.hasTPC()) {
      pidTpc = statusTpc(track);
    }
    int pidTof = TrackSelectorPID::NotApplicable;
    if (track.hasTOF()) {
      pidTof = statusTof(track);
    }

    if (pidTpc == TrackSelectorPID::Accepted && pidTof == TrackSelectorPID::Accepted) {
      return TrackSelectorPID::Accepted;
    }
    if (pidTpc == TrackSelectorPID::Accepted && (pidTof == TrackSelectorPID::NotApplicable || pidTof == TrackSelectorPID::Conditional)) {
      return TrackSelectorPID::Accepted;
    }
    if ((pidTpc == TrackSelectorPID::NotApplicable || pidTpc == TrackSelectorPID::Conditional) && pidTof == TrackSelectorPID::Accepted) {
      return TrackSelectorPID::Accepted;
    }
    if (pidTpc == TrackSelectorPID::Conditional && pidTof == TrackSelectorPID::Conditional) {
      return TrackSelectorPID::Accepted;
    }
    if (pidTpc == TrackSelectorPID::Rejected || pidTof == TrackSelectorPID::Rejected) {
      return TrackSelectorPID::Rejected;
    }
    return TrackSelectorPID::NotApplicable; // (NotApplicable for one detector) and (NotApplicable or Conditional for the other)
  }

  /// Checks whether a track is identified as electron and rejected as pion by TOF or RICH.
  /// \param track  track
  /// \param useTof  switch to use TOF
  /// \param useRich  switch to use RICH
  /// \return true if track is selected by TOF or RICH
  /// \note Ported from https://github.com/feisenhu/ALICE3-LoI-LMee/blob/main/efficiency/macros/anaEEstudy.cxx
  template <typename T>
  bool isElectronAndNotPion(const T& track, bool useTof = true, bool useRich = true)
  {
    bool isSelTof = false;
    bool isSelRich = false;
    bool hasRich = track.richId() > -1;
    bool hasTof = isValidForTof(track);
    auto nSigmaTofEl = mUseCustomNSigmaTofEl ? mCustomNSigmaTofEl : track.tofNSigmaEl();
    auto nSigmaTofPi = mUseCustomNSigmaTofPi ? mCustomNSigmaTofPi : track.tofNSigmaPi();
    auto nSigmaRichEl = hasRich ? track.rich().richNsigmaEl() : -1000.;
    auto nSigmaRichPi = hasRich ? track.rich().richNsigmaPi() : -1000.;
    auto p = track.p();

    // TOF
    if (useTof && hasTof && (p < 0.6)) {
      if (p > 0.4 && hasRich) {
        if ((std::abs(nSigmaTofEl) < mNSigmaTofMax) && (std::abs(nSigmaRichEl) < mNSigmaRichMax)) {
          isSelTof = true; // is selected as electron by TOF and RICH
        }
      } else if (p <= 0.4) {
        if (std::abs(nSigmaTofEl) < mNSigmaTofMax) {
          isSelTof = true; // is selected as electron by TOF
        }
      } else {
        isSelTof = false; // This is rejecting all the heavier particles which do not have a RICH signal in the p area of 0.4-0.6 GeV/c
      }
      if (std::abs(nSigmaTofPi) < mNSigmaTofMax) {
        isSelTof = false; // is selected as pion by TOF
      }
    } else {
      isSelTof = false;
    }

    // RICH
    if (useRich && hasRich) {
      if (std::abs(nSigmaRichEl) < mNSigmaRichMax) {
        isSelRich = true; // is selected as electron by RICH
      }
      if ((std::abs(nSigmaRichPi) < mNSigmaRichMax) && (p > 1.0) && (p < 2.0)) {
        isSelRich = false; // is selected as pion by RICH
      }
    } else {
      isSelRich = false;
    }

    return isSelRich || isSelTof;
  }

  // Bayesian

  /// Set pT range where Bayes PID is applicable.
  void setRangePtBayes(float ptMin, float ptMax)
  {
    mPtBayesMin = ptMin;
    mPtBayesMax = ptMax;
  }

  /// Set minimum Bayesian probability above which a track should be accepted.
  void setProbBayesMin(float cut)
  {
    mProbBayesMin = cut;
  }

  /// Checks if track is OK for Bayesian PID.
  /// \param track  track
  /// \return true if track is OK for Bayesian PID
  template <typename T>
  bool isValidForBayes(const T& track)
  {
    auto pt = track.pt();
    return (mPtBayesMin <= pt && pt <= mPtBayesMax);
  }

  /// Bayesian maximum probability algorithm.
  /// \param track  track
  /// \return true if selected species has the highest Bayesian probability
  template <typename T>
  bool isSelectedByBayes(const T& track)
  {
    // Get index of the most probable species for a given track.
    if constexpr (pdg == kElectron) {
      return track.bayesID() == o2::track::PID::Electron;
    } else if constexpr (pdg == kMuonMinus) {
      return track.bayesID() == o2::track::PID::Muon;
    } else if constexpr (pdg == kPiPlus) {
      return track.bayesID() == o2::track::PID::Pion;
    } else if constexpr (pdg == kKPlus) {
      return track.bayesID() == o2::track::PID::Kaon;
    } else if constexpr (pdg == kProton) {
      return track.bayesID() == o2::track::PID::Proton;
    } else {
      errorPdg();
      return false;
    }
  }

  /// Checks if track is compatible with given particle species hypothesis within given Bayesian probability range.
  /// \param track  track
  /// \return true if track satisfies PID hypothesis for given Bayesian probability range
  template <typename T>
  bool isSelectedByBayesProb(const T& track)
  {
    if (mProbBayesMin < 0.) { // switch off with negative values
      return true;
    }

    // Get probability for a given particle hypothesis.
    double prob = 0.;
    if constexpr (pdg == kElectron) {
      prob = track.bayesEl();
    } else if constexpr (pdg == kMuonMinus) {
      prob = track.bayesMu();
    } else if constexpr (pdg == kPiPlus) {
      prob = track.bayesPi();
    } else if constexpr (pdg == kKPlus) {
      prob = track.bayesKa();
    } else if constexpr (pdg == kProton) {
      prob = track.bayesPr();
    } else {
      errorPdg();
    }

    return mProbBayesMin <= prob;
  }

  /// Returns status of Bayesian PID selection for a given track, based on the most probable particle species.
  /// \param track  track
  /// \return Bayesian selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusBayes(const T& track)
  {
    if (!isValidForBayes(track)) {
      return TrackSelectorPID::NotApplicable;
    }
    if (isSelectedByBayes(track)) {
      return TrackSelectorPID::Accepted;
    } else {
      return TrackSelectorPID::Rejected;
    }
  }

  /// Returns status of Bayesian PID selection for a given track, based on the probability for a given particle species.
  /// \param track  track
  /// \return Bayesian selection status (see TrackSelectorPID::Status)
  template <typename T>
  TrackSelectorPID::Status statusBayesProb(const T& track)
  {
    if (!isValidForBayes(track)) {
      return TrackSelectorPID::NotApplicable;
    }
    if (isSelectedByBayesProb(track)) {
      return TrackSelectorPID::Accepted;
    } else {
      return TrackSelectorPID::Rejected;
    }
  }

 private:
  // TPC
  float mPtTpcMin = 0.;                ///< minimum pT for TPC PID [GeV/c]
  float mPtTpcMax = 100.;              ///< maximum pT for TPC PID [GeV/c]
  float mNSigmaTpcMin = -3.;           ///< minimum number of TPC σ
  float mNSigmaTpcMax = 3.;            ///< maximum number of TPC σ
  float mNSigmaTpcMinCondTof = 0.;     ///< minimum number of TPC σ if combined with TOF
  float mNSigmaTpcMaxCondTof = 0.;     ///< maximum number of TPC σ if combined with TOF
  bool mUseCustomNSigmaTpcPi = false;  ///< enable the usage of a custom value of TPC nσ vs. pion hypothesis
  bool mUseCustomNSigmaTpcKa = false;  ///< enable the usage of a custom value of TPC nσ vs. kaon hypothesis
  bool mUseCustomNSigmaTpcPr = false;  ///< enable the usage of a custom value of TPC nσ vs. proton hypothesis
  bool mUseCustomNSigmaTpcMu = false;  ///< enable the usage of a custom value of TPC nσ vs. muon hypothesis
  bool mUseCustomNSigmaTpcEl = false;  ///< enable the usage of a custom value of TPC nσ vs. electron hypothesis
  float mCustomNSigmaTpcPi = -99999.f; ///< custom value of TPC nσ vs. pion hypothesis to be used for selections
  float mCustomNSigmaTpcKa = -99999.f; ///< custom value of TPC nσ vs. kaon hypothesis to be used for selections
  float mCustomNSigmaTpcPr = -99999.f; ///< custom value of TPC nσ vs. proton hypothesis to be used for selections
  float mCustomNSigmaTpcMu = -99999.f; ///< custom value of TPC nσ vs. muon hypothesis to be used for selections
  float mCustomNSigmaTpcEl = -99999.f; ///< custom value of TPC nσ vs. electron hypothesis to be used for selections

  // TOF
  float mPtTofMin = 0.;                ///< minimum pT for TOF PID [GeV/c]
  float mPtTofMax = 100.;              ///< maximum pT for TOF PID [GeV/c]
  float mNSigmaTofMin = -3.;           ///< minimum number of TOF σ
  float mNSigmaTofMax = 3.;            ///< maximum number of TOF σ
  float mNSigmaTofMinCondTpc = 0.;     ///< minimum number of TOF σ if combined with TPC
  float mNSigmaTofMaxCondTpc = 0.;     ///< maximum number of TOF σ if combined with TPC
  bool mUseCustomNSigmaTofPi = false;  ///< enable the usage of a custom value of TOF nσ vs. pion hypothesis
  bool mUseCustomNSigmaTofKa = false;  ///< enable the usage of a custom value of TOF nσ vs. kaon hypothesis
  bool mUseCustomNSigmaTofPr = false;  ///< enable the usage of a custom value of TOF nσ vs. proton hypothesis
  bool mUseCustomNSigmaTofMu = false;  ///< enable the usage of a custom value of TOF nσ vs. muon hypothesis
  bool mUseCustomNSigmaTofEl = false;  ///< enable the usage of a custom value of TOF nσ vs. electron hypothesis
  float mCustomNSigmaTofPi = -99999.f; ///< custom value of TOF nσ vs. pion hypothesis to be used for selections
  float mCustomNSigmaTofKa = -99999.f; ///< custom value of TOF nσ vs. kaon hypothesis to be used for selections
  float mCustomNSigmaTofPr = -99999.f; ///< custom value of TOF nσ vs. proton hypothesis to be used for selections
  float mCustomNSigmaTofMu = -99999.f; ///< custom value of TOF nσ vs. muon hypothesis to be used for selections
  float mCustomNSigmaTofEl = -99999.f; ///< custom value of TOF nσ vs. electron hypothesis to be used for selections

  // RICH
  float mPtRichMin = 0.;            ///< minimum pT for RICH PID [GeV/c]
  float mPtRichMax = 100.;          ///< maximum pT for RICH PID [GeV/c]
  float mNSigmaRichMin = -3.;       ///< minimum number of RICH σ
  float mNSigmaRichMax = 3.;        ///< maximum number of RICH σ
  float mNSigmaRichMinCondTof = 0.; ///< minimum number of RICH σ if combined with TOF
  float mNSigmaRichMaxCondTof = 0.; ///< maximum number of RICH σ if combined with TOF

  // Bayesian
  float mPtBayesMin = 0.;    ///< minimum pT for Bayesian PID [GeV/c]
  float mPtBayesMax = 100.;  ///< maximum pT for Bayesian PID [GeV/c]
  float mProbBayesMin = -1.; ///< minimum Bayesian probability [%]

  /// Throw fatal for unsupported PDG values.
  void errorPdg()
  {
    LOGF(fatal, "Species with PDG code %d not supported", pdg);
  }
};

// Predefined types
using TrackSelectorEl = TrackSelectorPidBase<kElectron>;  // El
using TrackSelectorMu = TrackSelectorPidBase<kMuonMinus>; // Mu
using TrackSelectorPi = TrackSelectorPidBase<kPiPlus>;    // Pi
using TrackSelectorKa = TrackSelectorPidBase<kKPlus>;     // Ka
using TrackSelectorPr = TrackSelectorPidBase<kProton>;    // Pr

#endif // COMMON_CORE_TRACKSELECTORPID_H_
