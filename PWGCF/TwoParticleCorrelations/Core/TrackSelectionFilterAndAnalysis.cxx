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

/// \file TrackSelectionFilterAndAnalysis.cxx
/// \brief Implementation of the TracksSelectionFilterAndAnalysis class
/// for track filter and selection once filtered
/// \author victor.gonzalez.sebastian@gmail.com

#include "TrackSelectionFilterAndAnalysis.h"

#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

#include <Framework/Logger.h>

#include <TObjArray.h>
#include <TString.h>

#include <boost/regex/v5/regbase.hpp>
#include <boost/regex/v5/regex.hpp>
#include <boost/regex/v5/regex_fwd.hpp>
#include <boost/regex/v5/regex_search.hpp>

#include <Rtypes.h>

#include <cstdint>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::pwgcf;
using namespace boost;

/// \brief Constructor adapted to hyperloop
TrackSelectionConfigurable::TrackSelectionConfigurable(const std::vector<std::string>& ttype,
                                                       std::vector<std::string> nclstpc,
                                                       std::vector<std::string> nxrtpc,
                                                       std::vector<std::string> nclsits,
                                                       std::vector<std::string> chi2clustpc,
                                                       std::vector<std::string> chi2clusits,
                                                       std::vector<std::string> xrofctpc,
                                                       std::vector<std::string> dcaxy,
                                                       std::vector<std::string> dcaz,
                                                       std::vector<std::string> ptrange,
                                                       std::vector<std::string> etarange)

{
  auto storeCutString = [](auto& selvector, const std::string& selname) {
    if (selvector.size() != 0) {
      if (selvector.size() == 1) {
        if (selvector[0].size() != 0) {
          return TString::Format("%s{%s}", selname.c_str(), selvector[0].c_str());
        }
        return TString("");
      }
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
    return TString("");
  };
  /* track types are a bit specific */
  {
    TString ttypes = "ttype{";
    bool first = true;
    for (const auto& type : ttype) {
      if (not first) {
        ttypes += ',';
      }
      ttypes += type;
      first = false;
    }
    ttypes += "}";
    mTrackTypes = ttypes;
  }
  mNClustersTPC = storeCutString(nclstpc, "nclstpc");
  mNCrossedRowsTPC = storeCutString(nxrtpc, "nxrtpc");
  mNClustersITS = storeCutString(nclsits, "nclsits");
  mMaxChi2PerClusterTPC = storeCutString(chi2clustpc, "chi2clustpc");
  mMaxChi2PerClusterITS = storeCutString(chi2clusits, "chi2clusits");
  mMinNCrossedRowsOverFindableClustersTPC = storeCutString(xrofctpc, "xrofctpc");
  mMaxDcaXY = storeCutString(dcaxy, "dcaxy");
  mMaxDcaZ = storeCutString(dcaz, "dcaz");
  mPtRange = storeCutString(ptrange, "pT");
  mEtaRange = storeCutString(etarange, "eta");
}

ClassImp(TrackSelectionFilterAndAnalysis);

/// \brief Default constructor
TrackSelectionFilterAndAnalysis::TrackSelectionFilterAndAnalysis()
  : mNClustersTPC(nullptr),
    mNCrossedRowsTPC(nullptr),
    mNClustersITS(nullptr),
    mMaxChi2PerClusterTPC(nullptr),
    mMaxChi2PerClusterITS(nullptr),
    mMinNCrossedRowsOverFindableClustersTPC(nullptr),
    mMaxDcaXY(nullptr),
    mMaxDcaZ(nullptr)
{
  /* we own the track sign and types cuts objects */
  mTrackSign.SetOwner(true);
  mTrackTypes.SetOwner(true);
  /* at least we initialize by default pT and eta cuts */
  mPtRange = CutBrick<float>::constructBrick("pT", "rg{0.2,10}", std::set<std::string>{"rg"});
  mEtaRange = CutBrick<float>::constructBrick("eta", "rg{-0.8,0.8}", std::set<std::string>{"rg"});
}

/// \brief Constructor from regular expression
TrackSelectionFilterAndAnalysis::TrackSelectionFilterAndAnalysis(const TString& cutstr, Selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mNClustersTPC(nullptr),
    mNCrossedRowsTPC(nullptr),
    mNClustersITS(nullptr),
    mMaxChi2PerClusterTPC(nullptr),
    mMaxChi2PerClusterITS(nullptr),
    mMinNCrossedRowsOverFindableClustersTPC(nullptr),
    mMaxDcaXY(nullptr),
    mMaxDcaZ(nullptr)
{
  /* we own the track sign and types cuts objects */
  mTrackSign.SetOwner(true);
  mTrackTypes.SetOwner(true);
  /* at least we initialize by default pT and eta cuts */
  mPtRange = CutBrick<float>::constructBrick("pT", "rg{0.1,50}", std::set<std::string>{"rg"});
  mEtaRange = CutBrick<float>::constructBrick("eta", "rg{-1.0,1.0}", std::set<std::string>{"rg"});

  constructCutFromString(cutstr);
}

/// \brief Constructor from the track selection configurable
TrackSelectionFilterAndAnalysis::TrackSelectionFilterAndAnalysis(const TrackSelectionConfigurable& trcksel, Selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mNClustersTPC(nullptr),
    mNCrossedRowsTPC(nullptr),
    mNClustersITS(nullptr),
    mMaxChi2PerClusterTPC(nullptr),
    mMaxChi2PerClusterITS(nullptr),
    mMinNCrossedRowsOverFindableClustersTPC(nullptr),
    mMaxDcaXY(nullptr),
    mMaxDcaZ(nullptr),
    mPtRange(nullptr),
    mEtaRange(nullptr)
{
  constexpr int MaxMaskLength = 64;
  /* we own the track sign and types cuts objects */
  mTrackSign.SetOwner(true);
  mTrackTypes.SetOwner(true);

  TString cutString = "tracksel{" + trcksel.mTrackTypes;
  if (!trcksel.mNClustersTPC.empty() or
      !trcksel.mNCrossedRowsTPC.empty() or
      !trcksel.mNClustersITS.empty() or
      !trcksel.mMaxChi2PerClusterTPC.empty() or
      !trcksel.mMaxChi2PerClusterITS.empty() or
      !trcksel.mMinNCrossedRowsOverFindableClustersTPC.empty() or
      !trcksel.mMaxDcaXY.empty() or
      !trcksel.mMaxDcaZ.empty() or
      !trcksel.mPtRange.empty() or
      !trcksel.mEtaRange.empty()) {
    cutString += ";";
    bool first = true;

    auto appendCut = [&cutString, &first](const std::string& str) {
      if (!str.empty()) {
        if (first) {
          cutString += str;
          first = false;
        } else {
          cutString += "," + str;
        }
      }
    };
    appendCut(trcksel.mNClustersTPC);
    appendCut(trcksel.mNCrossedRowsTPC);
    appendCut(trcksel.mNClustersITS);
    appendCut(trcksel.mMaxChi2PerClusterTPC);
    appendCut(trcksel.mMaxChi2PerClusterITS);
    appendCut(trcksel.mMinNCrossedRowsOverFindableClustersTPC);
    appendCut(trcksel.mMaxDcaXY);
    appendCut(trcksel.mMaxDcaZ);
    appendCut(trcksel.mPtRange);
    appendCut(trcksel.mEtaRange);
  }
  cutString += "}";
  constructCutFromString(cutString);
  if (mMaskLength > MaxMaskLength) {
    LOGF(fatal, "TrackSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just %d available for the time being", mMaskLength, MaxMaskLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == KAnalysis) {
    /* TODO: check cuts consistency. In principle it should be there except for the ttype valid combinations */
    storeArmedMask();
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int TrackSelectionFilterAndAnalysis::calculateMaskLength()
{
  int length = 0;
  for (int i = 0; i < mTrackSign.GetEntries(); ++i) {
    length += (dynamic_cast<CutBrick<float>*>(mTrackSign.At(i)))->length();
  }
  for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
    length += (dynamic_cast<SpecialCutBrick*>(mTrackTypes.At(i)))->length();
  }
  if (mNClustersTPC != nullptr) {
    length += mNClustersTPC->length();
  }
  if (mNCrossedRowsTPC != nullptr) {
    length += mNCrossedRowsTPC->length();
  }
  if (mNClustersITS != nullptr) {
    length += mNClustersITS->length();
  }
  if (mMaxChi2PerClusterTPC != nullptr) {
    length += mMaxChi2PerClusterTPC->length();
  }
  if (mMaxChi2PerClusterITS != nullptr) {
    length += mMaxChi2PerClusterITS->length();
  }
  if (mMinNCrossedRowsOverFindableClustersTPC != nullptr) {
    length += mMinNCrossedRowsOverFindableClustersTPC->length();
  }
  if (mMaxDcaXY != nullptr) {
    length += mMaxDcaXY->length();
  }
  if (mMaxDcaZ != nullptr) {
    length += mMaxDcaZ->length();
  }
  if (mPtRange != nullptr) {
    length += mPtRange->length();
  }
  if (mEtaRange != nullptr) {
    length += mEtaRange->length();
  }
  return length;
}

void TrackSelectionFilterAndAnalysis::setPtRange(const TString& regex)
{

  delete mPtRange;

  mPtRange = CutBrick<float>::constructBrick("pT", regex.Data(), std::set<std::string>{"rg", "th", "lim", "xrg"});
  mMaskLength = calculateMaskLength();
}

void TrackSelectionFilterAndAnalysis::setEtaRange(const TString& regex)
{

  delete mEtaRange;

  mEtaRange = CutBrick<float>::constructBrick("eta", regex.Data(), std::set<std::string>{"rg", "th", "lim", "xrg"});
  mMaskLength = calculateMaskLength();
}

void TrackSelectionFilterAndAnalysis::constructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());
  /* let's catch the first level */
  regex cutregex("^tracksel\\{ttype\\{([\\w\\d,-]+)};?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  constexpr int MaxExpectedComponents = 3;
  if (not res or m.empty() or (m.size() > MaxExpectedComponents)) {
    LOGF(fatal, "TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};ncls{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
  }
  SetName("TrackSelectionFilterAndAnalysisCuts");
  SetTitle(cutstr.Data());

  /* let's introduce the charge sign bricks into the chain */
  {
    auto addChargeBrick = [this](auto name, auto regex) {
      std::set<std::string> allowed = {"lim", "th"};
      CutBrick<float>* brick = CutBrick<float>::constructBrick(name, regex, allowed);
      brick->arm(true);
      mTrackSign.Add(brick);
    };
    addChargeBrick("chp", "th{0}");
    addChargeBrick("chn", "lim{0}");
  }

  /* let's split the handling of track types and of its characteristics */
  /* let's handle the track types */
  {
    LOGF(info, "Captured %s", m[1].str().c_str());
    TObjArray* lev2toks = TString(m[1].str()).Tokenize(",");
    for (int i = 0; i < lev2toks->GetEntries(); ++i) {
      mTrackTypes.Add(new TrackSelectionBrick(lev2toks->At(i)->GetName()));
    }
    delete lev2toks;
  }
  /* let's now handle the reco track characteristics */
  {
    LOGF(info, "Captured %s", m[2].str().c_str());
    TString lev2str = m[2].str();
    while (lev2str.Length() != 0) {
      std::set<std::string> allowed = {"lim", "th", "rg", "xrg", "cwv"};
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
        Fatal("TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};nclstpc{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
      }
      LOGF(info, "Captured %s", m[1].str().c_str());
      auto storeIntCut = [&m, &allowed](CutBrick<int>*& brickvar) {
        delete brickvar;

        brickvar = CutBrick<int>::constructBrick(m[1].str().c_str(), m[2].str().c_str(), allowed);
      };
      auto storeFloatCut = [&m, &allowed](CutBrick<float>*& brickvar) {
        delete brickvar;

        brickvar = CutBrick<float>::constructBrick(m[1].str().c_str(), m[2].str().c_str(), allowed);
      };
      if (m[1].str() == "nclstpc") {
        storeIntCut(mNClustersTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableNClustersTPCCheck();
        }
      } else if (m[1].str() == "nclsits") {
        storeIntCut(mNClustersITS);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableNClustersITSCheck();
        }
      } else if (m[1].str() == "nxrtpc") {
        storeIntCut(mNCrossedRowsTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableNCrossedRowsTPCCheck();
        }
      } else if (m[1].str() == "chi2clustpc") {
        storeFloatCut(mMaxChi2PerClusterTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableMaxChi2PerClusterTPCCheck();
        }
      } else if (m[1].str() == "chi2clusits") {
        storeFloatCut(mMaxChi2PerClusterITS);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableMaxChi2PerClusterITSCheck();
        }
      } else if (m[1].str() == "xrofctpc") {
        storeFloatCut(mMinNCrossedRowsOverFindableClustersTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableMinNCrossedRowsOverFindableClustersTPCCheck();
        }
      } else if (m[1].str() == "dcaxy") {
        storeFloatCut(mMaxDcaXY);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableMaxDcaXYCheck();
        }
      } else if (m[1].str() == "dcaz") {
        storeFloatCut(mMaxDcaZ);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          (dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)))->disableMaxDcaZCheck();
        }
      } else if (m[1].str() == "pT") {
        storeFloatCut(mPtRange);
      } else if (m[1].str() == "eta") {
        storeFloatCut(mEtaRange);
      } else {
        Fatal("TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, cut on variable %s not implemented", cutstr.Data(), m[1].str().c_str());
      }
      /* removing the already handled track characteristics cut */
      lev2str.Remove(0, m[0].length());
    }
  }
  mMaskLength = calculateMaskLength();
}

/// \brief Fills the filter cuts mask
void TrackSelectionFilterAndAnalysis::storeArmedMask()
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
  for (int i = 0; i < mTrackSign.GetEntries(); ++i) {
    armedBrick(dynamic_cast<CutBrick<float>*>(mTrackSign.At(i)), true);
  }
  mOptArmedMask.push_back(optMask);
  optMask = 0UL;
  for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
    armedBrick(dynamic_cast<TrackSelectionBrick*>(mTrackTypes.At(i)), true);
  }
  mOptArmedMask.push_back(optMask);
  optMask = 0UL;
  if (mNClustersTPC != nullptr) {
    armedBrick(mNClustersTPC);
  }
  if (mNCrossedRowsTPC != nullptr) {
    armedBrick(mNCrossedRowsTPC);
  }
  if (mNClustersITS != nullptr) {
    armedBrick(mNClustersITS);
  }
  if (mMaxChi2PerClusterTPC != nullptr) {
    armedBrick(mMaxChi2PerClusterTPC);
  }
  if (mMaxChi2PerClusterITS != nullptr) {
    armedBrick(mMaxChi2PerClusterITS);
  }
  if (mMinNCrossedRowsOverFindableClustersTPC != nullptr) {
    armedBrick(mMinNCrossedRowsOverFindableClustersTPC);
  }
  if (mMaxDcaXY != nullptr) {
    armedBrick(mMaxDcaXY);
  }
  if (mMaxDcaZ != nullptr) {
    armedBrick(mMaxDcaZ);
  }
  /* for the time being the pT and eta bricks dont go to the mask */
  LOGF(info, "TrackSelectionFilterAndAnalysis::StoreArmedMask(), masks 0x%016lx, %s, 0x%016lx", armedMask, printOptionalMasks().Data(), forcedMask);
  mArmedMask = armedMask;
  mForcedArmedMask = forcedMask;
}
