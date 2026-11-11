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

/// \file EventSelectionFilterAndAnalysis.cxx
/// \brief Implementation of the EventSelectionFilterAndAnalysis class for
/// filtering collisions for derived data production
/// \author victor.gonzalez.sebastian@gmail.com

#include "EventSelectionFilterAndAnalysis.h"

#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

#include <Framework/Logger.h>

#include <TList.h>
#include <TObjArray.h>
#include <TString.h>

#include <boost/regex/v5/regbase.hpp>
#include <boost/regex/v5/regex.hpp>
#include <boost/regex/v5/regex_fwd.hpp>
#include <boost/regex/v5/regex_search.hpp>
#include <sys/types.h>

#include <Rtypes.h>

#include <cstdint>
#include <cstring>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::pwgcf;
using namespace boost;

EventSelectionConfigurable::EventSelectionConfigurable(const std::vector<std::string>& bfieldsel,
                                                       std::vector<std::string> multsel,
                                                       std::vector<std::string> trigsel,
                                                       std::vector<std::string> zvtxsel,
                                                       std::vector<std::string> pileuprej)
  : mBFiledSel(""),
    mMultSel(""),
    mTriggerSel(""),
    mZVertexSel(""),
    mPileUpRejection("")
{
  auto storeCutString = [](auto& selvector, const std::string& selname) {
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
        for (const auto& cut : selvector) {
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
  /* b-field is a bit specific */
  {
    TString bfield = "bfield{";
    bool first = true;
    for (const auto& pol : bfieldsel) {
      if (not first) {
        bfield += ',';
      }
      bfield += pol;
      first = false;
    }
    bfield += '}';
    mBFiledSel = bfield;
  }
  mMultSel = storeCutString(multsel, "centmult");
  mTriggerSel = storeCutString(trigsel, "trigger");
  mZVertexSel = storeCutString(zvtxsel, "zvtx");
  mPileUpRejection = storeCutString(pileuprej, "pileuprej");
}

ClassImp(EventSelectionFilterAndAnalysis);

const int nNoOfMultiplicityEstimators = 3;
const char* multiplicityEstimators[nNoOfMultiplicityEstimators] = {"V0M", "CL0", "CL1"};
const int nNoOfPileUpCorrelators = 3;
const char* multiplicityCorrelators[nNoOfPileUpCorrelators] = {"V0M_TPCOUT", "V0M_TRKLETS", "V0M_TPCCLSTS"};

/// \brief Default constructor
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis()
  : SelectionFilterAndAnalysis(),
    mBFieldSelection{},
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr)
{
}

/// \brief Constructor from regular expression
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis(const TString& cutstr, Selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mBFieldSelection{},
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr)
{
  constructCutFromString(cutstr);
}

/// \brief Constructor from the event selection configurable
EventSelectionFilterAndAnalysis::EventSelectionFilterAndAnalysis(const EventSelectionConfigurable& evtsel, Selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mBFieldSelection{},
    mMultiplicityClasses(nullptr),
    mTriggerSelection(nullptr),
    mZVertex(nullptr),
    mPileUpRejection(nullptr)
{
  TString cutString = "eventsel{";
  bool first = true;
  constexpr int MaskMaxLength = 64;

  auto appendCut = [&cutString, &first](const std::string& str) {
    if (str.size() > 0) {
      if (first) {
        cutString += str;
        first = false;
      } else {
        cutString += "," + str;
      }
    }
  };
  appendCut(evtsel.mBFiledSel);
  appendCut(evtsel.mMultSel);
  appendCut(evtsel.mTriggerSel);
  appendCut(evtsel.mZVertexSel);
  appendCut(evtsel.mPileUpRejection);

  cutString += "}";
  constructCutFromString(cutString);
  if (mMaskLength > MaskMaxLength) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just %d available for the time being", mMaskLength, MaskMaxLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == KAnalysis) {
    storeArmedMask();
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int EventSelectionFilterAndAnalysis::calculateMaskLength()
{
  int length = 0;
  auto addLength = [&length](auto& brick) {
    if (brick != nullptr) {
      length += brick->length();
      LOGF(info, "EventSelectionFilterAndAnalysis::calculateMaskLength, cumulated length: %d", length);
    }
  };

  for (const auto& brick : mBFieldSelection) {
    length += brick->length();
  }
  addLength(mMultiplicityClasses);
  addLength(mTriggerSelection);
  addLength(mZVertex);
  addLength(mPileUpRejection);

  return length;
}

void EventSelectionFilterAndAnalysis::MultiplicityBrick::initialize()
{
  /* initialize the storage of multiplicities from the different estimators */
  for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
    mValues.push_back(-1);
  }
  /* now initialize the filter machinery */
  if (mBrick->IsA() != CutWithVariations<float>::Class()) {
    /* no alternative estimators, just the default */
    LOGF(info, "Searching for default multiplicity estimator %s", mBrick->GetName());
    for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
      LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
      if (strcmp(multiplicityEstimators[i], mBrick->GetName()) == 0) {
        mDefaultEstimatorIndex = i;
        break;
      }
    }
    if (mDefaultEstimatorIndex < 0) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default multiplicity class estimator not implemented");
    }
  } else {
    /* we have alternative estimators our brick is of kind cwv */
    /* first the default */
    TList& deflst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getDefaultBricks();
    if (deflst.GetEntries() > 1) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() expected only one default multiplicity class estimator");
    }
    LOGF(info, "Searching for default multiplicity estimator %s", deflst.At(0)->GetName());
    for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
      LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
      if (strcmp(multiplicityEstimators[i], deflst.At(0)->GetName()) == 0) {
        mDefaultEstimatorIndex = i;
        break;
      }
    }
    if (mDefaultEstimatorIndex < 0) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default multiplicity class estimator not implemented");
    }
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int ivar = 0; ivar < varlst.GetEntries(); ++ivar) {
      LOGF(info, "Searching for variant multiplicity estimator %s", varlst.At(ivar)->GetName());
      for (int i = 0; i < nNoOfMultiplicityEstimators; ++i) {
        LOGF(info, "Checking multiplicity estimator %s", multiplicityEstimators[i]);
        if (strcmp(multiplicityEstimators[i], varlst.At(ivar)->GetName()) == 0) {
          mAlternateEstimatorIndex.push_back(i);
          break;
        }
      }
    }
    if (mAlternateEstimatorIndex.size() != (uint)(varlst.GetEntries())) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() variant multiplicity class estimators not all implemented");
    }
  }
}

void EventSelectionFilterAndAnalysis::PileUpRejBrick::initialize()
{
  /* initialize the storage of multiplicities from the different pile-up correlator estimators */
  for (int i = 0; i < nNoOfPileUpCorrelators; ++i) {
    mValues.push_back(-1);
    mIndepVar.push_back(-1);
  }
  /* now initialize the filter machinery */
  if (mBrick->IsA() != CutWithVariations<float>::Class()) {
    /* no alternative estimators, just the default */
    LOGF(info, "Searching for default pile-up correlator estimator %s", mBrick->GetName());
    for (int i = 0; i < nNoOfPileUpCorrelators; ++i) {
      LOGF(info, "Checking pile-up correlator estimator %s", multiplicityCorrelators[i]);
      if (strcmp(multiplicityCorrelators[i], mBrick->GetName()) == 0) {
        mDefaultEstimatorIndex = i;
        break;
      }
    }
    if (mDefaultEstimatorIndex < 0) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default pile-up correlator estimator not implemented");
    }
  } else {
    /* we have alternative estimators our brick is of kind cwv */
    /* first the default */
    TList& deflst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getDefaultBricks();
    if (deflst.GetEntries() > 1) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() expected only one default multiplicity class estimator");
    }
    LOGF(info, "Searching for default pile-up correlator estimator %s", deflst.At(0)->GetName());
    for (int i = 0; i < nNoOfPileUpCorrelators; ++i) {
      LOGF(info, "Checking pile-up correlator estimator %s", multiplicityCorrelators[i]);
      if (strcmp(multiplicityCorrelators[i], deflst.At(0)->GetName()) == 0) {
        mDefaultEstimatorIndex = i;
        break;
      }
    }
    if (mDefaultEstimatorIndex < 0) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() default pile-up correlator estimator not implemented");
    }
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int ivar = 0; ivar < varlst.GetEntries(); ++ivar) {
      LOGF(info, "Searching for variant pile-up correlator estimator %s", varlst.At(ivar)->GetName());
      for (int i = 0; i < nNoOfPileUpCorrelators; ++i) {
        LOGF(info, "Checking pile-up correlator estimator %s", multiplicityCorrelators[i]);
        if (strcmp(multiplicityCorrelators[i], varlst.At(ivar)->GetName()) == 0) {
          mAlternateEstimatorIndex.push_back(i);
          break;
        }
      }
    }
    if (mAlternateEstimatorIndex.size() != (uint)(varlst.GetEntries())) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::InitializeMultiplicityFilter() variant pile-up correlator estimators not all implemented");
    }
  }
}

void EventSelectionFilterAndAnalysis::constructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());
  /* let's catch the first level */
  regex cutregex("^eventsel\\{?(\\w+\\{[\\w\\d.,:{}=\\-\\+\\*\\/]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  constexpr int ExpectedMaxLength = 2;
  if (not res || m.empty() || (m.size() > ExpectedMaxLength)) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis::constructCutFromString", "Wrong RE: %s, try eventsel{bfield{positive-yes,negative-yes},zvtx{cwv{rg{-7.0,7.0}-yes:rg{-10.0,10.0}-no,rg{-3.0,3.0}-no}}} for instance", cutstr.Data());
  }
  this->SetName("EventSelectionFilterAndAnalysisCuts");
  this->SetTitle(cutstr.Data());

  /* let's now handle the event characteristics */
  {
    LOGF(info, "Captured %s", m[1].str().c_str());
    TString lev2str = m[1].str();
    while (lev2str.Length() != 0) {
      std::set<std::string> allowed = {"lim", "th", "rg", "xrg", "fnlim", "fnth", "fnrg", "fnxrg", "cwv", "mrg"};
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
      constexpr int ExpectedComponents = 3;

      if (not res or m.empty() or (m.size() != ExpectedComponents)) {
        Fatal("EventSelectionFilterAndAnalysis::::constructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};nclstpc{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
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
      /* b-field requires speciat treatment */
      auto storeBFieldCut = [&m](auto& brickvector) {
        LOGF(info, "Captured %s with %s", m[1].str().c_str(), m[2].str().c_str());
        auto storeBFPolarity = [&brickvector](auto regex) {
          auto addBFieldBrick = [&brickvector](auto name, auto regex, bool arm) {
            std::set<std::string> allowed = {"lim", "th"};
            CutBrick<float>* brick = CutBrick<float>::constructBrick(name, regex, allowed);
            brick->arm(arm);
            brickvector.push_back(brick);
          };
          if (TString(regex).Contains("positive")) {
            addBFieldBrick("bfp", "th{0}", TString(regex).Contains("-yes"));
          } else if (TString(regex).Contains("negative")) {
            addBFieldBrick("bfn", "lim{0}", TString(regex).Contains("-yes"));
          } else {
            LOGF(fatal, "storeBFPolarity(). Unknown polarity %s", regex);
          }
        };
        TObjArray* toks = TString(m[2].str()).Tokenize(",");
        for (int i = 0; i < toks->GetEntries(); ++i) {
          storeBFPolarity(toks->At(i)->GetName());
        }
        delete toks;
      };
      if (m[1].str() == "bfield") {
        storeBFieldCut(mBFieldSelection);
      } else if (m[1].str() == "centmult") {
        mMultiplicityClasses = new MultiplicityBrick();
        storeFloatCut(mMultiplicityClasses->mBrick);
        mMultiplicityClasses->initialize();
      } else if (m[1].str() == "mtrigg") {
        storeIntCut(mTriggerSelection);
      } else if (m[1].str() == "zvtx") {
        storeFloatCut(mZVertex);
      } else if (m[1].str() == "pileuprej") {
        mPileUpRejection = new PileUpRejBrick();
        storeFloatCut(mPileUpRejection->mBrick);
        mPileUpRejection->initialize();
      }
      /* removing the already handled event characteristics cut */
      lev2str.Remove(0, m[0].length());
    }
  }
  mMaskLength = calculateMaskLength();
}

///
void EventSelectionFilterAndAnalysis::ComplexBrickHelper::armedBrick(uint64_t& armedmask, uint64_t& optmask, uint64_t& forcedmask, int& bit)
{
  auto armedBrick = [&](auto brick, bool opt = false) {
    std::vector<bool> res = brick->isArmed();
    for (const auto& b : res) {
      if (b) {
        SETBIT(armedmask, bit);
        if (opt) {
          SETBIT(optmask, bit);
        } else {
          SETBIT(forcedmask, bit);
        }
      }
      bit++;
    }
  };

  if (mAlternateEstimatorIndex.size() > 0) {
    /* we have alternative estimators so our brick is of kind cwv */
    if (mBrick->IsA() != CutWithVariations<float>::Class()) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::filter() expected class with variations cut but it is not there");
    }
    /* first the default */
    TList& deflst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getDefaultBricks();
    if (deflst.GetEntries() > 1) {
      LOGF(fatal, "EventSelectionFilterAndAnalysis::filter() expected only one default multiplicity class estimator");
    }
    armedBrick((CutBrick<float>*)deflst.At(0), true);
    /* and now the variants */
    TList& varlst = dynamic_cast<CutWithVariations<float>*>(mBrick)->getVariantBricks();
    for (int i = 0; i < varlst.GetEntries(); ++i) {
      armedBrick((CutBrick<float>*)varlst.At(i), true);
    }
  } else {
    /* no alternative estimators, just the default */
    armedBrick(mBrick, true);
  }
}

/// \brief Fills the filter cuts mask
void EventSelectionFilterAndAnalysis::storeArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask.clear();
  mForcedArmedMask = 0UL;
  int bit = 0;

  auto armedBrick = [&](auto brick, bool opt = false) {
    std::vector<bool> res = brick->isArmed();
    for (const auto& b : res) {
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

  optMask = 0UL;
  for (const auto& brick : mBFieldSelection) {
    armedBrick(brick, true);
  }
  mOptArmedMask.push_back(optMask);
  optMask = 0UL;
  if (mMultiplicityClasses != nullptr) {
    mMultiplicityClasses->armedBrick(armedMask, optMask, forcedMask, bit);
  }
  mOptArmedMask.push_back(optMask);
  optMask = 0UL;
  if (mZVertex != nullptr) {
    armedBrick(mZVertex);
  }
  if (mPileUpRejection != nullptr) {
    mPileUpRejection->armedBrick(armedMask, optMask, forcedMask, bit);
  }
  LOGF(info, "EventSelectionFilterAndAnalysis::storeArmedMask(), masks 0x%016lx, %s, 0x%016lx", armedMask, printOptionalMasks().Data(), forcedMask);
  mArmedMask = armedMask;
  mForcedArmedMask = forcedMask;
}
