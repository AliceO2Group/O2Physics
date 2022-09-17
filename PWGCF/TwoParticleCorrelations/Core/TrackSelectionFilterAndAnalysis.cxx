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
#include "TrackSelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;
using namespace boost;

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
  /* we own the track types cuts objects */
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
  /* we own the track types cuts objects */
  mTrackTypes.SetOwner(true);
  /* at least we initialize by default pT and eta cuts */
  mPtRange = CutBrick<float>::constructBrick("pT", "rg{0.2,10}", std::set<std::string>{"rg"});
  mEtaRange = CutBrick<float>::constructBrick("eta", "rg{-0.8,0.8}", std::set<std::string>{"rg"});

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
  /* we own the track types cuts objects */
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
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int TrackSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
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
  regex cutregex("^tracksel\\{ttype\\{([\\w\\d,]+)};?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() > 3)) {
    Fatal("TrackSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try tracksel{ttype{FB32,FB96};ncls{th{70}},nxr{cwv{th{70},th{80}}}} for instance", cutstr.Data());
  }
  SetName("TrackSelectionFilterAndAnalysisCuts");
  SetTitle(cutstr.Data());

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

