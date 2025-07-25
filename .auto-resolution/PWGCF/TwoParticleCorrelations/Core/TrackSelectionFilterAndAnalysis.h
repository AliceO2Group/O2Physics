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
#ifndef TRACKSELECTIONFILTERANDANALYSIS_H
#define TRACKSELECTIONFILTERANDANALYSIS_H

#include "SelectionFilterAndAnalysis.h"
#include "SkimmingConfigurableCuts.h"

#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>

#include <fairlogger/Logger.h>

#include <Rtypes.h>

namespace o2
{
namespace analysis
{
namespace PWGCF
{
/* forward de declaration */
class TrackSelectionFilterAndAnalysis;

///\brief Convenience class for configurable access
class TrackSelectionConfigurable
{
  friend class TrackSelectionFilterAndAnalysis;

 public:
  TrackSelectionConfigurable(std::string ttype = "",
                             std::string nclstpc = "",
                             std::string nxrtpc = "",
                             std::string nclsits = "",
                             std::string chi2clustpc = "",
                             std::string chi2clusits = "",
                             std::string xrofctpc = "",
                             std::string dcaxy = "",
                             std::string dcaz = "",
                             std::string ptrange = "",
                             std::string etarange = "")
    : mTrackTypes{ttype},
      mNClustersTPC{nclstpc},
      mNCrossedRowsTPC{nxrtpc},
      mNClustersITS{nclsits},
      mMaxChi2PerClusterTPC{chi2clustpc},
      mMaxChi2PerClusterITS{chi2clusits},
      mMinNCrossedRowsOverFindableClustersTPC{xrofctpc},
      mMaxDcaXY{dcaxy},
      mMaxDcaZ{dcaz},
      mPtRange{ptrange},
      mEtaRange{etarange}
  {
  }
  TrackSelectionConfigurable(std::vector<std::string> ttype,
                             std::vector<std::string> nclstpc,
                             std::vector<std::string> nxrtpc,
                             std::vector<std::string> nclsits,
                             std::vector<std::string> chi2clustpc,
                             std::vector<std::string> chi2clusits,
                             std::vector<std::string> xrofctpc,
                             std::vector<std::string> dcaxy,
                             std::vector<std::string> dcaz,
                             std::vector<std::string> ptrange,
                             std::vector<std::string> etarange);

 private:
  std::string mTrackTypes = "";                             /// the track types to select list
  std::string mNClustersTPC = "";                           //! the number of TPC clusters cuts
  std::string mNCrossedRowsTPC = "";                        //! the number of TPC crossed rows cuts
  std::string mNClustersITS = "";                           //! the number of ITS clusters cuts
  std::string mMaxChi2PerClusterTPC = "";                   //! the max Chi2 per TPC cluster cuts
  std::string mMaxChi2PerClusterITS = "";                   //! the max Chi2 per ITS cluster cuts
  std::string mMinNCrossedRowsOverFindableClustersTPC = ""; //! the min ration crossed TPC rows over findable TPC clusters cuts
  std::string mMaxDcaXY = "";                               //! the DCAxy cuts
  std::string mMaxDcaZ = "";                                //! the DCAz cuts
  std::string mPtRange = "";                                //! the pT range cuts
  std::string mEtaRange = "";                               //! the eta range cuts

 private:
  ClassDefNV(TrackSelectionConfigurable, 1);
};

/// \brief Filter of tracks and track selection once filetered
class TrackSelectionFilterAndAnalysis : public SelectionFilterAndAnalysis
{
 public:
  TrackSelectionFilterAndAnalysis();
  TrackSelectionFilterAndAnalysis(const TString&, selmodes);
  TrackSelectionFilterAndAnalysis(const TrackSelectionConfigurable&, selmodes);

  void SetPtRange(const TString&);
  void SetEtaRange(const TString&);

  template <typename TrackToFilter>
  uint64_t Filter(TrackToFilter const& track);

 private:
  void ConstructCutFromString(const TString&);
  int CalculateMaskLength();
  void StoreArmedMask();

  TList mTrackSign;                                         /// the track charge sign list
  TList mTrackTypes;                                        /// the track types to select list
  CutBrick<int>* mNClustersTPC;                             //! the number of TPC clusters cuts
  CutBrick<int>* mNCrossedRowsTPC;                          //! the number of TPC crossed rows cuts
  CutBrick<int>* mNClustersITS;                             //! the number of ITS clusters cuts
  CutBrick<float>* mMaxChi2PerClusterTPC;                   //! the max Chi2 per TPC cluster cuts
  CutBrick<float>* mMaxChi2PerClusterITS;                   //! the max Chi2 per ITS cluster cuts
  CutBrick<float>* mMinNCrossedRowsOverFindableClustersTPC; //! the min ration crossed TPC rows over findable TPC clusters cuts
  CutBrick<float>* mMaxDcaXY;                               //! the DCAxy cuts
  CutBrick<float>* mMaxDcaZ;                                //! the DCAz cuts
  CutBrick<float>* mPtRange;                                //! the pT range cuts
  CutBrick<float>* mEtaRange;                               //! the eta range cuts

  ClassDef(TrackSelectionFilterAndAnalysis, 1)
};

/// \brief Fills the filter cuts mask
template <typename TrackToFilter>
uint64_t TrackSelectionFilterAndAnalysis::Filter(TrackToFilter const& track)
{
  uint64_t selectedMask = 0UL;
  int bit = 0;

  auto filterTrackType = [&](TrackSelectionBrick* ttype, auto trk) {
    if (ttype->Filter(trk)) {
      SETBIT(selectedMask, bit);
    }
    bit++;
  };

  auto filterBrickValue = [&](auto brick, auto value) {
    std::vector<bool> res = brick->Filter(value);
    for (auto b : res) {
      if (b) {
        SETBIT(selectedMask, bit);
      }
      bit++;
    }
  };

  auto filterBrickValueNoMask = [](auto brick, auto value) {
    std::vector<bool> res = brick->Filter(value);
    for (auto b : res) {
      if (b) {
        return true;
      }
    }
    return false;
  };

  for (int i = 0; i < mTrackSign.GetEntries(); ++i) {
    filterBrickValue((CutBrick<float>*)mTrackSign.At(i), track.sign());
  }
  for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
    filterTrackType((TrackSelectionBrick*)mTrackTypes.At(i), track);
  }
  if (mNClustersTPC != nullptr) {
    filterBrickValue(mNClustersTPC, track.tpcNClsFound());
  }
  if (mNCrossedRowsTPC != nullptr) {
    filterBrickValue(mNCrossedRowsTPC, track.tpcNClsCrossedRows());
  }
  if (mNClustersITS != nullptr) {
    filterBrickValue(mNClustersITS, track.itsNCls());
  }
  if (mMaxChi2PerClusterTPC != nullptr) {
    filterBrickValue(mMaxChi2PerClusterTPC, track.tpcChi2NCl());
  }
  if (mMaxChi2PerClusterITS != nullptr) {
    filterBrickValue(mMaxChi2PerClusterITS, track.itsChi2NCl());
  }
  if (mMinNCrossedRowsOverFindableClustersTPC != nullptr) {
    filterBrickValue(mMinNCrossedRowsOverFindableClustersTPC, track.tpcCrossedRowsOverFindableCls());
  }
  if (mMaxDcaXY != nullptr) {
    filterBrickValue(mMaxDcaXY, track.dcaXY());
  }
  if (mMaxDcaZ != nullptr) {
    filterBrickValue(mMaxDcaZ, track.dcaZ());
  }
  if (mPtRange != nullptr) {
    if (not filterBrickValueNoMask(mPtRange, track.pt())) {
      selectedMask = 0UL;
    }
  }
  if (mEtaRange != nullptr) {
    if (not filterBrickValueNoMask(mEtaRange, track.pt())) {
      selectedMask = 0UL;
    }
  }
  return mSelectedMask = selectedMask;
}

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // TRACKSELECTIONFILTERANDANALYSIS_H
