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

#include <fairlogger/Logger.h>
#include "TrackSelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::analysis::PWGCF;
using namespace boost;

/// \brief Constructor adapted to hyperloop
TrackSelectionConfigurable::TrackSelectionConfigurable(std::vector<std::string> ttype,
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
  : mTrackTypes{},
    mNClustersTPC{},
    mNCrossedRowsTPC{},
    mNClustersITS{},
    mMaxChi2PerClusterTPC{},
    mMaxChi2PerClusterITS{},
    mMinNCrossedRowsOverFindableClustersTPC{},
    mMaxDcaXY{},
    mMaxDcaZ{},
    mPtRange{},
    mEtaRange{}
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
  /* track types are a bit specific */
  {
    TString ttypes = "ttype{";
    bool first = true;
    for (auto& type : ttype) {
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
  : SelectionFilterAndAnalysis(),
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
  mPtRange = CutBrick<float>::constructBrick("pT", "rg{0.2,10}", std::set<std::string>{"rg"});
  mEtaRange = CutBrick<float>::constructBrick("eta", "rg{-0.8,0.8}", std::set<std::string>{"rg"});
}

/// \brief Constructor from regular expression
TrackSelectionFilterAndAnalysis::TrackSelectionFilterAndAnalysis(const TString& cutstr, selmodes mode)
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

  ConstructCutFromString(cutstr);
}

/// \brief Constructor from the track selection configurable
TrackSelectionFilterAndAnalysis::TrackSelectionFilterAndAnalysis(const TrackSelectionConfigurable& trcksel, selmodes mode)
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
  /* we own the track sign and types cuts objects */
  mTrackSign.SetOwner(true);
  mTrackTypes.SetOwner(true);

  TString cutString = "tracksel{" + trcksel.mTrackTypes;
  if (trcksel.mNClustersTPC.size() > 0 or
      trcksel.mNCrossedRowsTPC.size() > 0 or
      trcksel.mNClustersITS.size() > 0 or
      trcksel.mMaxChi2PerClusterTPC.size() > 0 or
      trcksel.mMaxChi2PerClusterITS.size() > 0 or
      trcksel.mMinNCrossedRowsOverFindableClustersTPC.size() > 0 or
      trcksel.mMaxDcaXY.size() > 0 or
      trcksel.mMaxDcaZ.size() > 0 or
      trcksel.mPtRange.size() > 0 or
      trcksel.mEtaRange.size() > 0) {
    cutString += ";";
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
  ConstructCutFromString(cutString);
  if (mMaskLength > 64) {
    LOGF(fatal, "TrackSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just 64 available for the time being", mMaskLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == kAnalysis) {
    /* TODO: check cuts consistency. In principle it should be there except for the ttype valid combinations */
    StoreArmedMask();
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int TrackSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
  for (int i = 0; i < mTrackSign.GetEntries(); ++i) {
    length += ((CutBrick<float>*)mTrackSign.At(i))->Length();
  }
  for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
    length += ((SpecialCutBrick*)mTrackTypes.At(i))->Length();
  }
  if (mNClustersTPC != nullptr) {
    length += mNClustersTPC->Length();
  }
  if (mNCrossedRowsTPC != nullptr) {
    length += mNCrossedRowsTPC->Length();
  }
  if (mNClustersITS != nullptr) {
    length += mNClustersITS->Length();
  }
  if (mMaxChi2PerClusterTPC != nullptr) {
    length += mMaxChi2PerClusterTPC->Length();
  }
  if (mMaxChi2PerClusterITS != nullptr) {
    length += mMaxChi2PerClusterITS->Length();
  }
  if (mMinNCrossedRowsOverFindableClustersTPC != nullptr) {
    length += mMinNCrossedRowsOverFindableClustersTPC->Length();
  }
  if (mMaxDcaXY != nullptr) {
    length += mMaxDcaXY->Length();
  }
  if (mMaxDcaZ != nullptr) {
    length += mMaxDcaZ->Length();
  }
  if (mPtRange != nullptr) {
    length += mPtRange->Length();
  }
  if (mEtaRange != nullptr) {
    length += mEtaRange->Length();
  }
  return length;
}

void TrackSelectionFilterAndAnalysis::SetPtRange(const TString& regex)
{
  if (mPtRange != nullptr) {
    delete mPtRange;
  }
  mPtRange = CutBrick<float>::constructBrick("pT", regex.Data(), std::set<std::string>{"rg", "th", "lim", "xrg"});
  mMaskLength = CalculateMaskLength();
}

void TrackSelectionFilterAndAnalysis::SetEtaRange(const TString& regex)
{
  if (mEtaRange != nullptr) {
    delete mEtaRange;
  }
  mEtaRange = CutBrick<float>::constructBrick("eta", regex.Data(), std::set<std::string>{"rg", "th", "lim", "xrg"});
  mMaskLength = CalculateMaskLength();
}

void TrackSelectionFilterAndAnalysis::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());
  /* let's catch the first level */
  regex cutregex("^tracksel\\{ttype\\{([\\w\\d,-]+)};?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() > 3)) {
    LOGF(fatal, "TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};ncls{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
  }
  SetName("TrackSelectionFilterAndAnalysisCuts");
  SetTitle(cutstr.Data());

  /* let's introduce the charge sign bricks into the chain */
  {
    auto addChargeBrick = [this](auto name, auto regex) {
      std::set<std::string> allowed = {"lim", "th"};
      CutBrick<float>* brick = CutBrick<float>::constructBrick(name, regex, allowed);
      brick->Arm(true);
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

      if (not res or m.empty() or (m.size() != 3)) {
        Fatal("TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};nclstpc{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
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
      if (m[1].str() == "nclstpc") {
        storeIntCut(mNClustersTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableNClustersTPCCheck();
        }
      } else if (m[1].str() == "nclsits") {
        storeIntCut(mNClustersITS);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableNClustersITSCheck();
        }
      } else if (m[1].str() == "nxrtpc") {
        storeIntCut(mNCrossedRowsTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableNCrossedRowsTPCCheck();
        }
      } else if (m[1].str() == "chi2clustpc") {
        storeFloatCut(mMaxChi2PerClusterTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableMaxChi2PerClusterTPCCheck();
        }
      } else if (m[1].str() == "chi2clusits") {
        storeFloatCut(mMaxChi2PerClusterITS);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableMaxChi2PerClusterITSCheck();
        }
      } else if (m[1].str() == "xrofctpc") {
        storeFloatCut(mMinNCrossedRowsOverFindableClustersTPC);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableMinNCrossedRowsOverFindableClustersTPCCheck();
        }
      } else if (m[1].str() == "dcaxy") {
        storeFloatCut(mMaxDcaXY);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableMaxDcaXYCheck();
        }
      } else if (m[1].str() == "dcaz") {
        storeFloatCut(mMaxDcaZ);
        for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
          ((TrackSelectionBrick*)mTrackTypes.At(i))->DisableMaxDcaZCheck();
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
  mMaskLength = CalculateMaskLength();
}

/// \brief Fills the filter cuts mask
void TrackSelectionFilterAndAnalysis::StoreArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask.clear();
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

  optMask = 0UL;
  for (int i = 0; i < mTrackSign.GetEntries(); ++i) {
    armedBrick((CutBrick<float>*)mTrackSign.At(i), true);
  }
  mOptArmedMask.push_back(optMask);
  optMask = 0UL;
  for (int i = 0; i < mTrackTypes.GetEntries(); ++i) {
    armedBrick((TrackSelectionBrick*)mTrackTypes.At(i), true);
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
