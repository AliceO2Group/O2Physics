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

/// \file EventSelectionFilterAndAnalysis.h
/// \brief EventSelectionFilterAndAnalysis class for filtering collisions for derived data production
/// \author victor.gonzalez.sebastian@gmail.com

#ifndef PWGCF_TWOPARTICLECORRELATIONS_CORE_EVENTSELECTIONFILTERANDANALYSIS_H_
#define PWGCF_TWOPARTICLECORRELATIONS_CORE_EVENTSELECTIONFILTERANDANALYSIS_H_

#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

#include <Framework/Logger.h>

#include <TList.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace o2
{
namespace analysis
{
namespace pwgcf
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
    : mBFiledSel(std::move(bfieldsel)), mMultSel{std::move(multsel)}, mTriggerSel{std::move(trigsel)}, mZVertexSel{std::move(zvtxsel)}, mPileUpRejection{std::move(pileuprej)}
  {
  }
  EventSelectionConfigurable(const std::vector<std::string>& bfieldsel,
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
  EventSelectionFilterAndAnalysis(const TString&, Selmodes);
  EventSelectionFilterAndAnalysis(const EventSelectionConfigurable&, Selmodes);

  template <typename CollisionToFilter, typename AssociatedTracks>
  uint64_t filter(CollisionToFilter const& col, AssociatedTracks const& trks, int bfield);
  std::vector<float> getMultiplicities() { return (mMultiplicityClasses != nullptr) ? mMultiplicityClasses->getMultiplicities() : std::vector<float>{}; }
  /// \brief Gets the index of the active multiplicity value within the multiplicities array
  int getMultiplicityIndex() { return (mMultiplicityClasses != nullptr) ? mMultiplicityClasses->getArmedIndex() : -1; }

 private:
  struct ComplexBrickHelper {
    virtual ~ComplexBrickHelper() = default;
    CutBrick<float>* mBrick = nullptr;
    int mDefaultEstimatorIndex = -1;
    std::vector<int> mAlternateEstimatorIndex = std::vector<int>{};
    std::vector<float> mValues = std::vector<float>{};
    int length() { return mBrick->length(); }
    virtual void initialize() = 0;
    virtual bool filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index) = 0;
    bool filter(uint64_t& mask, int& bit);
    void armedBrick(uint64_t&, uint64_t&, uint64_t&, int&);
  };
  struct MultiplicityBrick : public ComplexBrickHelper {
    void initialize() override;
    std::vector<float> getMultiplicities();
    int getArmedIndex() { return mBrick->getArmedIndex(); }
    bool filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index) override;
  };
  struct PileUpRejBrick : public ComplexBrickHelper {
    void initialize() override;
    std::vector<float> mIndepVar = std::vector<float>{};
    bool filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index) override;
  };

  static bool filterBrickValue(uint64_t& mask, int& bit, CutBrick<float>* brick, float value);
  void constructCutFromString(const TString&);
  template <typename CollisionToFilter, typename AssociatedTracks>
  void storeMultiplicities(CollisionToFilter const&, AssociatedTracks const&);
  int calculateMaskLength() override;
  void storeArmedMask() override;

  std::vector<CutBrick<float>*> mBFieldSelection; //! the magnetic field selection cuts
  MultiplicityBrick* mMultiplicityClasses;        //! the multiplicity classes cuts
  CutBrick<int>* mTriggerSelection;               //! the trigger selection cuts
  CutBrick<float>* mZVertex;                      //! the z vertex selection cuts
  PileUpRejBrick* mPileUpRejection;               //! the pile-up rejection criteria

  ClassDefOverride(EventSelectionFilterAndAnalysis, 1)
};

/// \brief Stores the different multiplicities needed for proper collision filtering
/// TODO: perhaps in a next iteration check if some storage can be removed
template <typename CollisionToFilter, typename AssociatedTracks>
inline void EventSelectionFilterAndAnalysis::storeMultiplicities(CollisionToFilter const& col, AssociatedTracks const& tracks)
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
  constexpr float StandardMinPt = 0.15;
  constexpr float StandardMaxEta = 0.8;
  for (const auto& track : tracks) {
    if (track.pt() > StandardMinPt && std::abs(track.eta()) < StandardMaxEta && track.hasTPC()) {
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
  std::vector<bool> res = brick->filter(value);
  for (const auto& b : res) {
    if (b) {
      atleastone = true;
      SETBIT(mask, bit);
    }
    bit++;
  }
  return atleastone;
};

inline bool EventSelectionFilterAndAnalysis::ComplexBrickHelper::filter(uint64_t& mask, int& bit)
{
  bool acceptcollision = true;
  if (mAlternateEstimatorIndex.size() > 0) {
    bool atleastonealternative = false;
    /*   if (mBrick->IsA() != CutWithVariations<float>::Class()) {
          LOGF(fatal, "EventSelectionFilterAndAnalysis::filter() expected class with variations cut but it is not there");
        } */
    /* first the default */
    TList& deflst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getDefaultBricks();
    bool acc = filter(mask, bit, (CutBrick<float>*)deflst.At(0), mDefaultEstimatorIndex);
    atleastonealternative = atleastonealternative || acc;
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int i = 0; i < varlst.GetEntries(); ++i) {
      bool acc = filter(mask, bit, (CutBrick<float>*)varlst.At(i), mAlternateEstimatorIndex[i]);
      atleastonealternative = atleastonealternative || acc;
    }
    acceptcollision = acceptcollision && atleastonealternative;
  } else {
    /* no alternative estimators, just the default */
    bool acc = filter(mask, bit, mBrick, mDefaultEstimatorIndex);
    acceptcollision = acceptcollision && acc;
  }
  return acceptcollision;
}

inline bool EventSelectionFilterAndAnalysis::MultiplicityBrick::filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index)
{
  return filterBrickValue(mask, bit, brick, mValues[index]);
}

inline bool EventSelectionFilterAndAnalysis::PileUpRejBrick::filter(uint64_t& mask, int& bit, CutBrick<float>* brick, int index)
{
  brick->setIndependentFnVar(mIndepVar[index]);
  return filterBrickValue(mask, bit, brick, mValues[index]);
}

/// \brief Fills the filter cuts mask
template <typename CollisionToFilter, typename AssociatedTracks>
inline uint64_t EventSelectionFilterAndAnalysis::filter(CollisionToFilter const& col, const AssociatedTracks& trks, int bfield)
{
  storeMultiplicities(col, trks);
  uint64_t selectedMask = 0UL;
  mSelectedMask = 0UL;
  int bit = 0;

  /* we require the collision be accepted by the whole set of filters */
  bool acceptcollision = true;
  if (mBFieldSelection.size() > 0) {
    bool oneatleast = false;
    for (const auto& brick : mBFieldSelection) {
      bool acc = filterBrickValue(selectedMask, bit, brick, bfield);
      oneatleast = oneatleast || acc;
    }
    acceptcollision = acceptcollision && oneatleast;
  }
  if (mMultiplicityClasses != nullptr) {
    bool acc = mMultiplicityClasses->ComplexBrickHelper::filter(selectedMask, bit);
    acceptcollision = acceptcollision && acc;
  }
  if (mZVertex != nullptr) {
    bool acc = filterBrickValue(selectedMask, bit, mZVertex, col.posZ());
    acceptcollision = acceptcollision && acc;
  }
  if (mPileUpRejection != nullptr) {
    bool acc = mPileUpRejection->ComplexBrickHelper::filter(selectedMask, bit);
    acceptcollision = acceptcollision && acc;
  }
  LOGF(debug, "EventSelectionFilterAndAnalysis::filter(), %s collision", acceptcollision ? "accepted" : "rejected");
  if (acceptcollision) {
    mSelectedMask = selectedMask;
  }
  return mSelectedMask;
}

/// \brief Fills the filter cuts mask
inline std::vector<float> EventSelectionFilterAndAnalysis::MultiplicityBrick::getMultiplicities()
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

} // namespace pwgcf
} // namespace analysis
} // namespace o2

#endif // PWGCF_TWOPARTICLECORRELATIONS_CORE_EVENTSELECTIONFILTERANDANALYSIS_H_
