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

#include <boost/regex.hpp>
#include <TObjArray.h>

#include "Framework/AnalysisTask.h"
#include "EventSelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::analysis::PWGCF;
using namespace boost;

EventSelectionConfigurable::EventSelectionConfigurable(std::vector<std::string> multsel,
                                                       std::vector<std::string> trigsel,
                                                       std::vector<std::string> zvtxsel,
                                                       std::vector<std::string> pileuprej)
  : mMultSel(""),
    mTriggerSel(""),
    mZVertexSel(""),
    mPileUpRejection("")
{
  auto storeCutString = [](auto& selvector, std::string selname) {
    if (selvector.size() != 0) {
      if (selvector.size() == 1) {
        if (selvector[0].size() != 0) {
          return TString::Format("%s{%s}", selname.c_str(), selvector[0].c_str());
        } else {
          return TString("");
        }
      } else {
        TString scut = selname + "{cwv{";
        bool def = true;
        bool firstvar = true;
        for (auto cut : selvector) {
          if (def) {
            scut += cut + ':';
            def = false;
          } else {
            if (not firstvar) {
              scut += ',';
            }
            scut += cut;
            firstvar = false;
          }
        }
        scut += "}}";
        return scut;
      }
    }
    return TString("");
  };
  mMultSel = storeCutString(multsel, "centmult");
  mTriggerSel = storeCutString(trigsel, "trigger");
  mZVertexSel = storeCutString(zvtxsel, "zvtx");
  mPileUpRejection = storeCutString(pileuprej, "pileup");
}

ClassImp(EventSelectionFilterAndAnalysis);

const int nNoOfMultiplicityEstimators = 3;
const char* multiplicityEstimators[nNoOfMultiplicityEstimators] = {"V0M", "CL0", "CL1"};

/// \brief Default constructor
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis()
  : SelectionFilterAndAnalysis(),
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr),
    mMultiplicities{},
    mDefaultMultiplicityEstimatorIndex(-1),
    mAlternateMultiplicityEstimatorIndex{}
{
}

/// \brief Constructor from regular expression
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis(const TString& cutstr, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr),
    mMultiplicities{},
    mDefaultMultiplicityEstimatorIndex(-1),
    mAlternateMultiplicityEstimatorIndex{}
{
  ConstructCutFromString(cutstr);
}

/// \brief Constructor from the event selection configurable
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis(const EventSelectionConfigurable& evtsel, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr),
    mMultiplicities{},
    mDefaultMultiplicityEstimatorIndex(-1),
    mAlternateMultiplicityEstimatorIndex{}
{
  TString cutString = "eventsel{";
  bool first = true;

  auto appendCut = [&cutString, &first](std::string str) {
    if (str.size() > 0) {
      if (first) {
        cutString += str;
        first = false;
      } else {
        cutString += "," + str;
      }
    }
  };
  appendCut(evtsel.mMultSel);
  appendCut(evtsel.mTriggerSel);
  appendCut(evtsel.mZVertexSel);
  appendCut(evtsel.mPileUpRejection);

  cutString += "}";
  ConstructCutFromString(cutString);
  if (mMaskLength > 64) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just 64 available for the time being", mMaskLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == kAnalysis) {
    StoreArmedMask();
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int EventSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
  auto addLength = [&length](auto& brick) {
    if (brick != nullptr) {
      length += brick->Length();
      LOGF(info, "EventSelectionFilterAndAnalysis::CalculateMaskLength, cumulated length: %d", length);
    }
  };

  addLength(mMultiplicityClasses);
  addLength(mTriggerSelection);
  addLength(mZVertex);
  addLength(mPileUpRejection);

  return length;
}

void EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter()
{
  /* initialize the storage of multiplicities from the differente estimators */
  for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
    mMultiplicities.push_back(-1);
  }
  /* now initialize the filter machinery */
  if (mMultiplicityClasses != nullptr) {
    if (mMultiplicityClasses->IsA() != CutWithVariations<float>::Class()) {
      /* no alternative estimators, just the default */
      LOGF(info, "Searching for default multiplicity estimator %s", mMultiplicityClasses->GetName());
      for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
        LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
        if (strcmp(multiplicityEstimators[i], mMultiplicityClasses->GetName()) == 0) {
          mDefaultMultiplicityEstimatorIndex = i;
          break;
        }
      }
      if (mDefaultMultiplicityEstimatorIndex < 0) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default multiplicity class estimator not implemented");
      }
    } else {
      /* we have alternative estimators our brick is of kind cwv */
      /* first the default */
      TList& deflst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getDefaultBricks();
      if (deflst.GetEntries() > 1) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() expected only one default multiplicity class estimator");
      }
      LOGF(info, "Searching for default multiplicity estimator %s", deflst.At(0)->GetName());
      for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
        LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
        if (strcmp(multiplicityEstimators[i], deflst.At(0)->GetName()) == 0) {
          mDefaultMultiplicityEstimatorIndex = i;
          break;
        }
      }
      if (mDefaultMultiplicityEstimatorIndex < 0) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default multiplicity class estimator not implemented");
      }
      /* and now the variants */
      TList& varlst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getVariantBricks();
      for (int ivar = 0; ivar < varlst.GetEntries(); ++ivar) {
        LOGF(info, "Searching for variant multiplicity estimator %s", varlst.At(ivar)->GetName());
        for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
          LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
          if (strcmp(multiplicityEstimators[i], varlst.At(ivar)->GetName()) == 0) {
            mAlternateMultiplicityEstimatorIndex.push_back(i);
            break;
          }
        }
      }
      if (mAlternateMultiplicityEstimatorIndex.size() != (uint)(varlst.GetEntries())) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() variant multiplicity class estimators not all implemented");
      }
    }
  } else {
    LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter expected a multiplicity filter but it is not there");
  }
}

void EventSelectionFilterAndAnalysis::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr);
  /* let's catch the first level */
  regex cutregex("^eventsel\\{?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() > 2)) {
    Fatal("EventSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};ncls{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
  }
  this->SetName("EventSelectionFilterAndAnalysisCuts");
  this->SetTitle(cutstr.Data());

  /* let's now handle the event characteristics */
  {
    LOGF(info, "Captured %s", m[1].str().c_str());
    TString lev2str = m[1].str();
    while (lev2str.Length() != 0) {
      std::set<std::string> allowed = {"lim", "th", "rg", "xrg", "cwv", "mrg"};
      lev2str.Remove(TString::kLeading, ' ');
      lev2str.Remove(TString::kLeading, ',');
      lev2str.Remove(TString::kLeading, ' ');
      if (lev2str.Length() == 0) {
        break;
      }
      regex cutregex("^(\\w+)\\{((?:[^{}]++|\\{(?2)\\})*)\\}");
      std::string in(lev2str.Data());
      smatch m;
      bool res = regex_search(in, m, cutregex);

      if (not res or m.empty() or (m.size() != 3)) {
        Fatal("EventSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};nclstpc{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
      }
      LOGF(info, "Captured %s", m[1].str().c_str());
      auto storeIntCut = [&m, &allowed](CutBrick<int>*& brickvar) {
        if (brickvar != nullptr) {
          delete brickvar;
        }
        brickvar = CutBrick<int>::constructBrick(m[1].str().c_str(), m[2].str().c_str(), allowed);
      };
      auto storeFloatCut = [&m, &allowed](CutBrick<float>*& brickvar) {
        if (brickvar != nullptr) {
          delete brickvar;
        }
        brickvar = CutBrick<float>::constructBrick(m[1].str().c_str(), m[2].str().c_str(), allowed);
      };
      if (m[1].str() == "centmult") {
        storeFloatCut(mMultiplicityClasses);
        InitializeMultiplicityFilter();
      } else if (m[1].str() == "mtrigg") {
        storeIntCut(mTriggerSelection);
      } else if (m[1].str() == "zvtx") {
        storeFloatCut(mZVertex);
      } else if (m[1].str() == "pileup") {
        storeIntCut(mPileUpRejection);
      }
      /* removing the already handled event characteristics cut */
      lev2str.Remove(0, m[0].length());
    }
  }
  mMaskLength = CalculateMaskLength();
}

/// \brief Fills the filter cuts mask
void EventSelectionFilterAndAnalysis::StoreArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask = 0UL;
  mForcedArmedMask = 0UL;
  int bit = 0;

  auto armedBrick = [&](auto brick, bool opt = false) {
    std::vector<bool> res = brick->IsArmed();
    for (auto b : res) {
      if (b) {
        SETBIT(armedMask, bit);
        if (opt) {
          SETBIT(optMask, bit);
        } else {
          SETBIT(forcedMask, bit);
        }
      }
      bit++;
    }
  };

  if (mMultiplicityClasses != nullptr) {
    if (mAlternateMultiplicityEstimatorIndex.size() > 0) {
      /* we have alternative estimators so our brick is of kind cwv */
      if (mMultiplicityClasses->IsA() != CutWithVariations<float>::Class()) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter() expected class with variations cut but it is not there");
      }
      /* first the default */
      TList& deflst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getDefaultBricks();
      if (deflst.GetEntries() > 1) {
        LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter() expected only one default multiplicity class estimator");
      }
      armedBrick((CutBrick<float>*)deflst.At(0), true);
      /* and now the variants */
      TList& varlst = dynamic_cast<CutWithVariations<float>*>(mMultiplicityClasses)->getVariantBricks();
      for (int i = 0; i < varlst.GetEntries(); ++i) {
        if (varlst.At(i)->IsA() != CutBrickSelectorMultipleRanges<float>::Class()) {
          LOGF(fatal, "EventSelectionFilterAndAnalysis::Filter, expected a multirange selector");
        }
        armedBrick((CutBrick<float>*)varlst.At(i), true);
      }
    } else {
      /* no alternative estimators, just the default */
      armedBrick(mMultiplicityClasses, true);
    }
  }
  if (mTriggerSelection != nullptr) {
  }
  if (mZVertex != nullptr) {
    armedBrick(mZVertex);
  }
  if (mPileUpRejection != nullptr) {
  }
  LOGF(info, "EventSelectionFilterAndAnalysis::StoreArmedMask(), masks 0x%08lx, 0x%08lx, 0x%08lx", armedMask, optMask, forcedMask);
  mArmedMask = armedMask;
  mOptArmedMask = optMask;
  mForcedArmedMask = forcedMask;
}
