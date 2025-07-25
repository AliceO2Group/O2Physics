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
#ifndef EVENTSELECTIONFILTERANDANALYSIS_H
#define EVENTSELECTIONFILTERANDANALYSIS_H

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
/* forward declaration */
class EventSelectionFilterAndAnalysis;

///\brief Convenience class for configurable access
class EventSelectionConfigurable
{
  friend class EventSelectionFilterAndAnalysis;

 public:
  EventSelectionConfigurable(std::string bfieldsel = "",
                             std::string multsel = "",
                             std::string trigsel = "",
                             std::string zvtxsel = "",
                             std::string pileuprej = "")
    : mBFiledSel(bfieldsel), mMultSel{multsel}, mTriggerSel{trigsel}, mZVertexSel{zvtxsel}, mPileUpRejection{pileuprej}
  {
  }
  EventSelectionConfigurable(std::vector<std::string> bfieldsel,
                             std::vector<std::string> multsel,
                             std::vector<std::string> trigsel,
                             std::vector<std::string> zvtxsel,
                             std::vector<std::string> pileuprej);

 private:
  std::string mBFiledSel = "";       //! the magnetic field selection cuts
  std::string mMultSel = "";         //! the multiplicity selection cuts
  std::string mTriggerSel = "";      //! the trigger selection cuts
  std::string mZVertexSel = "";      //! the z vertex selection cuts
  std::string mPileUpRejection = ""; //! the pile-up rejection criteria

 private:
  ClassDefNV(EventSelectionConfigurable, 1);
};

/// \brief Filter of tracks and track selection once filetered
class EventSelectionFilterAndAnalysis : public SelectionFilterAndAnalysis
{
 public:
  EventSelectionFilterAndAnalysis();
  EventSelectionFilterAndAnalysis(const TString&, selmodes);
  EventSelectionFilterAndAnalysis(const EventSelectionConfigurable&, selmodes);

  template <typename CollisionToFilter, typename AssociatedTracks>
  uint64_t Filter(CollisionToFilter const& col, AssociatedTracks const& trks, int bfield);
  std::vector<float> GetMultiplicities() { return (mMultiplicityClasses != nullptr) ? mMultiplicityClasses->GetMultiplicities() : std::vector<float>{}; }
  /// \brief Gets the index of the active multiplicity value within the multiplicities array
  int getMultiplicityIndex() { return (mMultiplicityClasses != nullptr) ? mMultiplicityClasses->GetArmedIndex() : -1; }

 private:
  struct ComplexBrickHelper {
    CutBrick<float>* mBrick = nullptr;
    int mDefaultEstimatorIndex = -1;
    std::vector<int> mAlternateEstimatorIndex = std::vector<int>{};
    std::vector<float> mValues = std::vector<float>{};
    int Length() { return mBrick->Length(); }
    virtual void initialize() = 0;
    virtual bool Filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index) = 0;
    bool Filter(uint64_t& mask, int& bit);
    void armedBrick(uint64_t&, uint64_t&, uint64_t&, int&);
  };
  struct MultiplicityBrick : public ComplexBrickHelper {
    virtual void initialize();
    std::vector<float> GetMultiplicities();
    int GetArmedIndex() { return mBrick->getArmedIndex(); }
    virtual bool Filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index);
  };
  struct PileUpRejBrick : public ComplexBrickHelper {
    virtual void initialize();
    std::vector<float> mIndepVar = std::vector<float>{};
    virtual bool Filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index);
  };

  static bool filterBrickValue(uint64_t& mask, int& bit, CutBrick<float>* brick, float value);
  void ConstructCutFromString(const TString&);
  template <typename CollisionToFilter, typename AssociatedTracks>
  void StoreMultiplicities(CollisionToFilter const&, AssociatedTracks const&);
  int CalculateMaskLength() override;
  virtual void StoreArmedMask() override;

  std::vector<CutBrick<float>*> mBFieldSelection; //! the magnetic field selection cuts
  MultiplicityBrick* mMultiplicityClasses;        //! the multiplicity classes cuts
  CutBrick<int>* mTriggerSelection;               //! the trigger selection cuts
  CutBrick<float>* mZVertex;                      //! the z vertex selection cuts
  PileUpRejBrick* mPileUpRejection;               //! the pile-up rejection criteria

  ClassDef(EventSelectionFilterAndAnalysis, 1)
};

/// \brief Stores the different multiplicities needed for proper collision filtering
/// TODO: perhaps in a next iteration check if some storage can be removed
template <typename CollisionToFilter, typename AssociatedTracks>
inline void EventSelectionFilterAndAnalysis::StoreMultiplicities(CollisionToFilter const& col, AssociatedTracks const& tracks)
{
  /* store the collision multiplicities for the different estimators */
  /* TODO: we need to adapt this to the Run 3 scenario */
  /* this is hard coded with the estimator index declaration */
  mMultiplicityClasses->mValues[0] = col.centRun2V0M();
  mMultiplicityClasses->mValues[1] = col.centRun2CL0();
  mMultiplicityClasses->mValues[2] = col.centRun2CL1();

  int ntracksTPCout = 0;
  int nTPCClusters = 0; // no TPC clusters available at collision level for the time being

  /* number of tracks with TPCout */
  /* we require tracks with pT > 0.15 and |eta| < 0.8 and additionally the TPCout flag */
  /* for the time being there is not TPCout flag, we just consider hasTPC              */
  for (auto track : tracks) {
    if (track.pt() > 0.15 and abs(track.eta()) < 0.8 and track.hasTPC()) {
      ntracksTPCout++;
    }
  }
  /* store the collision multiplicities for the different pile-up correlator estimators */
  /* TODO: we need to adapt this to the Run 3 scenario */
  /* this is hard coded with the pile-up correlator estimator index declaration */
  mPileUpRejection->mValues[0] = col.centRun2V0M();
  mPileUpRejection->mIndepVar[0] = ntracksTPCout;
  mPileUpRejection->mValues[1] = col.centRun2V0M();
  mPileUpRejection->mIndepVar[1] = col.multTracklets();
  mPileUpRejection->mValues[2] = col.centRun2V0M();
  mPileUpRejection->mIndepVar[0] = nTPCClusters;
}

inline bool EventSelectionFilterAndAnalysis::filterBrickValue(uint64_t& mask, int& bit, CutBrick<float>* brick, float value)
{
  bool atleastone = false;
  std::vector<bool> res = brick->Filter(value);
  for (auto b : res) {
    if (b) {
      atleastone = true;
      SETBIT(mask, bit);
    }
    bit++;
  }
  return atleastone;
};

inline bool EventSelectionFilterAndAnalysis::ComplexBrickHelper::Filter(uint64_t& mask, int& bit)
{
  bool acceptcollision = true;
  if (mAlternateEstimatorIndex.size() > 0) {
    bool atleastonealternative = false;
    /*   if (mBrick->IsA() != CutWithVariations<float>::Class()) {
          LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter() expected class with variations cut but it is not there");
        } */
    /* first the default */
    TList& deflst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getDefaultBricks();
    bool acc = Filter(mask, bit, (CutBrick<float>*)deflst.At(0), mDefaultEstimatorIndex);
    atleastonealternative = atleastonealternative || acc;
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int i = 0; i < varlst.GetEntries(); ++i) {
      bool acc = Filter(mask, bit, (CutBrick<float>*)varlst.At(i), mAlternateEstimatorIndex[i]);
      atleastonealternative = atleastonealternative || acc;
    }
    acceptcollision = acceptcollision && atleastonealternative;
  } else {
    /* no alternative estimators, just the default */
    bool acc = Filter(mask, bit, mBrick, mDefaultEstimatorIndex);
    acceptcollision = acceptcollision && acc;
  }
  return acceptcollision;
}

inline bool EventSelectionFilterAndAnalysis::MultiplicityBrick::Filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index)
{
  return filterBrickValue(mask, bit, brick, mValues[index]);
}

inline bool EventSelectionFilterAndAnalysis::PileUpRejBrick::Filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index)
{
  brick->setIndependentFnVar(mIndepVar[index]);
  return filterBrickValue(mask, bit, brick, mValues[index]);
}

/// \brief Fills the filter cuts mask
template <typename CollisionToFilter, typename AssociatedTracks>
inline uint64_t EventSelectionFilterAndAnalysis::Filter(CollisionToFilter const& col, const AssociatedTracks& trks, int bfield)
{
  StoreMultiplicities(col, trks);
  uint64_t selectedMask = 0UL;
  mSelectedMask = 0UL;
  int bit = 0;

  /* we require the collision be accepted by the whole set of filters */
  bool acceptcollision = true;
  if (mBFieldSelection.size() > 0) {
    bool oneatleast = false;
    for (auto brick : mBFieldSelection) {
      bool acc = filterBrickValue(selectedMask, bit, brick, bfield);
      oneatleast = oneatleast || acc;
    }
    acceptcollision = acceptcollision && oneatleast;
  }
  if (mMultiplicityClasses != nullptr) {
    bool acc = mMultiplicityClasses->ComplexBrickHelper::Filter(selectedMask, bit);
    acceptcollision = acceptcollision && acc;
  }
  if (mTriggerSelection != nullptr) {
  }
  if (mZVertex != nullptr) {
    bool acc = filterBrickValue(selectedMask, bit, mZVertex, col.posZ());
    acceptcollision = acceptcollision && acc;
  }
  if (mPileUpRejection != nullptr) {
    bool acc = mPileUpRejection->ComplexBrickHelper::Filter(selectedMask, bit);
    acceptcollision = acceptcollision && acc;
  }
  LOGF(debug, "EventSelectionFilterAndAnalysis::Filter(), %s collision", acceptcollision ? "accepted" : "rejected");
  if (acceptcollision) {
    mSelectedMask = selectedMask;
  }
  return mSelectedMask;
}

/// \brief Fills the filter cuts mask
inline std::vector<float> EventSelectionFilterAndAnalysis::MultiplicityBrick::GetMultiplicities()
{
  std::vector<float> res;

  if (mAlternateEstimatorIndex.size() > 0) {
    /* first the default */
    res.push_back(mValues[mDefaultEstimatorIndex]);
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int i = 0; i < varlst.GetEntries(); ++i) {
      res.push_back(mValues[mAlternateEstimatorIndex[i]]);
    }
  } else {
    /* no alternative estimators, just the default */
    res.push_back(mValues[mDefaultEstimatorIndex]);
  }
  return res;
}

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // EVENTSELECTIONFILTERANDANALYSIS_H
