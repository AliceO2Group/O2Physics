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
// Class for track selection
//

#ifndef PWGEM_DILEPTON_CORE_EMTRACKCUT_H_
#define PWGEM_DILEPTON_CORE_EMTRACKCUT_H_

#include "Framework/Logger.h"

#include "TNamed.h"

class EMTrackCut : public TNamed
{
 public:
  EMTrackCut() = default;
  EMTrackCut(const char* name, const char* title) : TNamed(name, title) {}
  ~EMTrackCut() {}

  enum class EMTrackCuts : int {
    // track cut
    kTrackPtRange,
    kTrackEtaRange,
    kTrackPhiRange,
    kTrackBit,
    kNCuts
  };

  template <typename TTrack>
  bool IsSelected(TTrack const& track) const
  {
    if (!IsSelectedTrack(track, EMTrackCuts::kTrackPtRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, EMTrackCuts::kTrackEtaRange)) {
      return false;
    }
    if (!IsSelectedTrack(track, EMTrackCuts::kTrackPhiRange)) {
      return false;
    }

    if (!IsSelectedTrack(track, EMTrackCuts::kTrackBit)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool IsSelectedTrack(T const& track, const EMTrackCuts& cut) const
  {
    switch (cut) {
      case EMTrackCuts::kTrackPtRange:
        return track.pt() > mMinTrackPt && track.pt() < mMaxTrackPt;

      case EMTrackCuts::kTrackEtaRange:
        return track.eta() > mMinTrackEta && track.eta() < mMaxTrackEta;

      case EMTrackCuts::kTrackPhiRange:
        return track.phi() > mMinTrackPhi && track.phi() < mMaxTrackPhi;

      case EMTrackCuts::kTrackBit: {
        // for (int i = 0; i < 10; i++) {
        //   if ((mTrackBit & (1 << i)) > 0 && !((track.trackBit() & (1 << i)) > 0)) {
        //     return false;
        //   }
        // }
        // return true;
        return (track.trackBit() & mTrackBit) >= mTrackBit;
      }

      default:
        return false;
    }
  }

  // Setters
  void SetTrackPtRange(float minPt = 0.f, float maxPt = 1e10f);
  void SetTrackEtaRange(float minEta = -1e10f, float maxEta = 1e10f);
  void SetTrackPhiRange(float minPhi = 0.f, float maxPhi = 6.3f);
  void SetTrackBit(uint16_t bits);

 private:
  // kinematic cuts
  float mMinTrackPt{0.f}, mMaxTrackPt{1e10f};      // range in pT
  float mMinTrackEta{-1e10f}, mMaxTrackEta{1e10f}; // range in eta
  float mMinTrackPhi{0.f}, mMaxTrackPhi{6.3};      // range in phi

  // track quality cuts
  uint16_t mTrackBit{0};
  // std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  ClassDef(EMTrackCut, 1);
};

#endif // PWGEM_DILEPTON_CORE_EMTRACKCUT_H_
