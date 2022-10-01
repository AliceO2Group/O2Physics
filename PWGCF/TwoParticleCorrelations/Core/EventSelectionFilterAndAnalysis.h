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
/* forward declaration */
class EventSelectionFilterAndAnalysis;

///\brief Convenince class for configurable access
class EventSelectionConfigurable
{
  friend class EventSelectionFilterAndAnalysis;

 public:
  EventSelectionConfigurable(std::string multsel = "",
                             std::string trigsel = "",
                             std::string zvtxsel = "",
                             std::string pileuprej = "")
    : mMultSel{multsel}, mTriggerSel{trigsel}, mZVertexSel{zvtxsel}, mPileUpRejection{pileuprej}
  {
  }
  EventSelectionConfigurable(std::vector<std::string> multsel,
                             std::vector<std::string> trigsel,
                             std::vector<std::string> zvtxsel,
                             std::vector<std::string> pileuprej);

 private:
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

  template <typename CollisionToFilter>
  uint64_t Filter(CollisionToFilter const& col);
  std::vector<float> GetMultiplicities();
  /// \brief Gets the index of the active multiplicity value within the multiplicities array
  int getMultiplicityIndex() { return mMultiplicityClasses->getArmedIndex(); }

 private:
  void ConstructCutFromString(const TString&);
  void InitializeMultiplicityFilter();
  int CalculateMaskLength() override;
  virtual void StoreArmedMask() override;

  CutBrick<float>* mMultiplicityClasses;                 //! the multiplicity default classes cuts
  CutBrick<int>* mTriggerSelection;                      //! the trigger selection cuts
  CutBrick<float>* mZVertex;                             //! the z vertex selection cuts
  CutBrick<int>* mPileUpRejection;                       //! the pile-up rejection criteria
  std::vector<float> mMultiplicities;                    // the collision multiplicities from the different implemented estimators
  int mDefaultMultiplicityEstimatorIndex;                // the default estimator index on the collision multiplicities array
  std::vector<int> mAlternateMultiplicityEstimatorIndex; // the vector of alternate estimators index on the collision multiplicities array

  ClassDef(EventSelectionFilterAndAnalysis, 1)
};

/// \brief Fills the filter cuts mask
template <typename CollisionToFilter>
inline uint64_t EventSelectionFilterAndAnalysis::Filter(CollisionToFilter const& col)
{
  /* store the collision multiplicities for the different estimators */
  /* TODO: we need to adapt this to the Run 3 scenario */
  /* this is hard coded with the estimator index declaration */
  mMultiplicities[0] = col.centRun2V0M();
  mMultiplicities[1] = col.centRun2CL0();
  mMultiplicities[2] = col.centRun2CL1();

  uint64_t selectedMask = 0UL;
  mSelectedMask = 0UL;
  int bit = 0;

  auto filterBrickValue = [&](auto brick, auto value) {
    bool atleastone = false;
    std::vector<bool> res = brick->Filter(value);
    for (auto b : res) {
      if (b) {
        atleastone = true;
        SETBIT(selectedMask, bit);
      }
      bit++;
    }
    return atleastone;
  };

  /* we require the collision be accepted by the whole set of filters */
  bool acceptcollision = true;
  if (mMultiplicityClasses != nullptr) {
    if (mAlternateMultiplicityEstimatorIndex.size() > 0) {
      bool atleastonealternative = false;
      /* we have alternative estimators so our brick is of kind cwv */
      if (mMultiplicityClasses->IsA() != CutWithVariations<float>::Class()) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter() expected class with variations cut but it is not there");
      }
      /* first the default */
      TList& deflst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getDefaultBricks();
      if (deflst.GetEntries() > 1) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter() expected only one default multiplicity class estimator");
      }
      bool acc = filterBrickValue((CutBrick<float>*)deflst.At(0), mMultiplicities[mDefaultMultiplicityEstimatorIndex]);
      atleastonealternative = atleastonealternative || acc;
      /* and now the variants */
      TList& varlst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getVariantBricks();
      for (int i = 0; i < varlst.GetEntries(); ++i) {
        if (varlst.At(i)->IsA() != CutBrickSelectorMultipleRanges<float>::Class()) {
          LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter, expected a multirange selector");
        }
        bool acc = filterBrickValue((CutBrick<float>*)varlst.At(i), mMultiplicities[mAlternateMultiplicityEstimatorIndex[i]]);
        atleastonealternative = atleastonealternative || acc;
      }
      acceptcollision = acceptcollision && atleastonealternative;
    } else {
      /* no alternative estimators, just the default */
      bool acc = filterBrickValue(mMultiplicityClasses, mMultiplicities[mDefaultMultiplicityEstimatorIndex]);
      acceptcollision = acceptcollision && acc;
    }
  }
  if (mTriggerSelection != nullptr) {
  }
  if (mZVertex != nullptr) {
    bool acc = filterBrickValue(mZVertex, col.posZ());
    acceptcollision = acceptcollision && acc;
  }
  if (mPileUpRejection != nullptr) {
  }
  LOGF(debug, "EventSelectionFilterAndAnalysis::Filter(), %s collision", acceptcollision ? "accepted" : "rejected");
  if (acceptcollision) {
    mSelectedMask = selectedMask;
  }
  return mSelectedMask;
}

/// \brief Fills the filter cuts mask
inline std::vector<float> EventSelectionFilterAndAnalysis::GetMultiplicities()
{
  std::vector<float> res;

  if (mMultiplicityClasses != nullptr) {
    if (mAlternateMultiplicityEstimatorIndex.size() > 0) {
      /* first the default */
      res.push_back(mMultiplicities[mDefaultMultiplicityEstimatorIndex]);
      /* and now the variants */
      TList& varlst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getVariantBricks();
      for (int i = 0; i < varlst.GetEntries(); ++i) {
        res.push_back(mMultiplicities[mAlternateMultiplicityEstimatorIndex[i]]);
      }
    } else {
      /* no alternative estimators, just the default */
      res.push_back(mMultiplicities[mDefaultMultiplicityEstimatorIndex]);
    }
  }
  return res;
}

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // EVENTSELECTIONFILTERANDANALYSIS_H
