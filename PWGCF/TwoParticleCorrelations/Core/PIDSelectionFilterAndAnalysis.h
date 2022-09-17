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

#include <Rtypes.h>
#include <TString.h>
#include <TObject.h>
#include <TNamed.h>
#include <TList.h>

#include "SkimmingConfigurableCuts.h"
#include "SelectionFilterAndAnalysis.h"

namespace o2
{
namespace analysis
{
namespace PWGCF
{

/// \brief Filter of tracks based on PID and track selection once filetered
class PIDSelectionFilterAndAnalysis : public SelectionFilterAndAnalysis
{
 public:
  PIDSelectionFilterAndAnalysis();
  PIDSelectionFilterAndAnalysis(const TString&, selmodes);

  void SetPTOF(float ptof) { mPTOF = ptof; }
  void SetRequireTOF(bool requiretof = false) { mRequireTOF = requiretof; }
  void SetEllipticTPCTOF(bool elliptic = false) { mEllipticTPCTOF = elliptic; }

  template <typename TrackToFilter>
  uint64_t Filter(TrackToFilter const& track);

 private:
  void ConstructCutFromString(const TString&);
  int CalculateMaskLength();

  float mPTOF = 0.8f;                                 ///< the p threshold for cheking TOF information
  bool mRequireTOF = false;                           ///< is TOF required
  bool mEllipticTPCTOF = false;                       ///< 2D nsigmas elliptic TPC+TOF
  std::vector<PIDSelectionBrick*> mInclusiveTrackPID; ///< the list of species wanted to be detected
  std::vector<PIDSelectionBrick*> mExclusiveTrackPID; ///< the list of species wanted to be rejected

  ClassDef(PIDSelectionFilterAndAnalysis, 1)
};

/// \brief Fills the filter cuts mask
template <typename TrackToFilter>
uint64_t PIDSelectionFilterAndAnalysis::Filter(TrackToFilter const& track)
{
  uint64_t selectedMask = 0UL;
  int bit = 0;

  auto filterBrickTrack = [&](auto brick, auto track) {
    std::vector<bool> res = brick->Filter(track);
    for (auto b : res) {
      if (b) {
        SETBIT(selectedMask, bit++);
      }
    }
  };

  for (auto bricksp : mInclusiveTrackPID) {
    filterBrickTrack(bricksp, track);
  }
  for (auto bricksp : mExclusiveTrackPID) {
    filterBrickTrack(bricksp, track);
  }
  return mSelectedMask = selectedMask;
}

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // PIDSELECTIONFILTERANDANALYSIS_H
