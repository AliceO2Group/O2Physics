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
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "EventSelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;
using namespace boost;

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

/// \brief Constructor from the track selection configurable
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
  if (mMaskLength > 64) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just 64 available for the time being", mMaskLength);
  }
}

