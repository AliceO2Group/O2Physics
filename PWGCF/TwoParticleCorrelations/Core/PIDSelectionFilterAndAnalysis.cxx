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

#include "PIDSelectionFilterAndAnalysis.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

#include <TObjArray.h>

#include <boost/regex.hpp>

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;
using namespace boost;

#ifdef INCORPORATEBAYESIANPID
PIDSelectionConfigurable::PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                                                   std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr,
                                                   std::vector<std::string> pidbayel, std::vector<std::string> pidbaymu, std::vector<std::string> pidbaypi, std::vector<std::string> pidbayka, std::vector<std::string> pidbaypr)
#else
PIDSelectionConfigurable::PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                                                   std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr)
#endif
  : mPidTpcSel_el{},
    mPidTpcSel_mu{},
    mPidTpcSel_pi{},
    mPidTpcSel_ka{},
    mPidTpcSel_pr{},
    mPidTofSel_el{},
    mPidTofSel_mu{},
    mPidTofSel_pi{},
    mPidTofSel_ka{},
    mPidTofSel_pr{},
    mPidBaySel_el{},
    mPidBaySel_mu{},
    mPidBaySel_pi{},
    mPidBaySel_ka{},
    mPidBaySel_pr{}
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
  mPidTpcSel_el = storeCutString(pidtpcel, "tpcel");
  mPidTpcSel_mu = storeCutString(pidtpcmu, "tpcmu");
  mPidTpcSel_pi = storeCutString(pidtpcpi, "tpcpi");
  mPidTpcSel_ka = storeCutString(pidtpcka, "tpcka");
  mPidTpcSel_pr = storeCutString(pidtpcpr, "tpcpr");
  mPidTofSel_el = storeCutString(pidtofel, "tofel");
  mPidTofSel_mu = storeCutString(pidtofmu, "tofmu");
  mPidTofSel_pi = storeCutString(pidtofpi, "tofpi");
  mPidTofSel_ka = storeCutString(pidtofka, "tofka");
  mPidTofSel_pr = storeCutString(pidtofpr, "tofpr");
#ifdef INCORPORATEBAYESIANPID
  mPidBaySel_el = storeCutString(pidbayel, "bayel");
  mPidBaySel_mu = storeCutString(pidbaymu, "baymu");
  mPidBaySel_pi = storeCutString(pidbaypi, "baypi");
  mPidBaySel_ka = storeCutString(pidbayka, "bayka");
  mPidBaySel_pr = storeCutString(pidbaypr, "baypr");
#endif
}

ClassImp(PIDSelectionFilterAndAnalysis);

/// \brief species supported
const std::vector<std::string> PIDSelectionFilterAndAnalysis::mSpeciesNames = {"el", "mu", "pi", "ka", "pr"};

/// \brief Default constructor
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis()
  : SelectionFilterAndAnalysis(),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr),
    mBayesProbability(kNoOfSpecies, nullptr)
{
}

/// \brief Constructor from regular expression
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const TString& cutstr, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr),
    mBayesProbability(kNoOfSpecies, nullptr)
{
  ConstructCutFromString(cutstr);
}

/// \brief Constructor from the PID selection configurable
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const PIDSelectionConfigurable& pidsel, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr),
    mBayesProbability(kNoOfSpecies, nullptr)
{
  TString cutString = "pidsel{";
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
  cutString += "tpcsel{";
  first = true;
  appendCut(pidsel.mPidTpcSel_el);
  appendCut(pidsel.mPidTpcSel_mu);
  appendCut(pidsel.mPidTpcSel_pi);
  appendCut(pidsel.mPidTpcSel_ka);
  appendCut(pidsel.mPidTpcSel_pr);
  cutString += "},tofsel{";
  first = true;
  appendCut(pidsel.mPidTofSel_el);
  appendCut(pidsel.mPidTofSel_mu);
  appendCut(pidsel.mPidTofSel_pi);
  appendCut(pidsel.mPidTofSel_ka);
  appendCut(pidsel.mPidTofSel_pr);
  cutString += "},baysel{";
  first = true;
  appendCut(pidsel.mPidBaySel_el);
  appendCut(pidsel.mPidBaySel_mu);
  appendCut(pidsel.mPidBaySel_pi);
  appendCut(pidsel.mPidBaySel_ka);
  appendCut(pidsel.mPidBaySel_pr);

  cutString += "}}";
  ConstructCutFromString(cutString);
  if (mMaskLength > 64) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just 64 available for the time being", mMaskLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == kAnalysis) {
    StoreArmedMask();
  }
}

/// \brief Destructor
/// Releases the taken memory
PIDSelectionFilterAndAnalysis::~PIDSelectionFilterAndAnalysis()
{
  for (auto brick : mCloseNsigmasTPC) {
    delete brick;
  }
  for (auto brick : mCloseNsigmasTOF) {
    delete brick;
  }
  for (auto brick : mBayesProbability) {
    delete brick;
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int PIDSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
  auto addLength = [&](auto bricklst) {
    for (auto brick : bricklst) {
      if (brick != nullptr) {
        length += brick->Length();
      }
    }
  };
  addLength(mCloseNsigmasTPC);
  addLength(mCloseNsigmasTOF);
  addLength(mBayesProbability);
  return length;
}

void PIDSelectionFilterAndAnalysis::ConstructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());
  /* let's catch the first level */
  regex cutregex("^pidsel\\{?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  if (not res or m.empty() or (m.size() > 2)) {
    Fatal("PIDSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE: %s, try pidsel{tpcsel{tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}}} for instance", cutstr.Data());
  }
  this->SetName("EventSelectionFilterAndAnalysisCuts");
  this->SetTitle(cutstr.Data());

  /* let's now handle each of the detectors requirements */
  {
    LOGF(info, "Captured %s", m[1].str().c_str());
    TString lev2str = m[1].str();
    while (lev2str.Length() != 0) {
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
        LOGF(fatal, "PIDSelectionFilterAndAnalysis::::ConstructCutFromString. Wrong RE: %s, try pidsel{tpcsel{tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}}} for instance", cutstr.Data());
      }
      LOGF(info, "Captured %s", m[1].str().c_str());
      auto handleDetectorLevel = [](std::string detector, auto& bricklst, std::string cut) {
        std::set<std::string> allowed = {"lim", "th", "rg", "xrg", "cwv", "mrg"};
        TString lev3str = cut;
        while (lev3str.Length() != 0) {
          regex cutregex3l("^(\\w+)\\{((?:[^{}]++|\\{(?2)\\})*)\\}");
          smatch m3l;
          auto storeCut = [&m3l, &allowed](CutBrick<float>*& brickvar) {
            if (brickvar != nullptr) {
              delete brickvar;
            }
            brickvar = CutBrick<float>::constructBrick(m3l[1].str().c_str(), m3l[2].str().c_str(), allowed);
          };
          lev3str.Remove(TString::kLeading, ' ');
          lev3str.Remove(TString::kLeading, ',');
          lev3str.Remove(TString::kLeading, ' ');
          if (lev3str.Length() == 0) {
            break;
          }
          std::string in3l(lev3str.Data());
          bool res3l = regex_search(in3l, m3l, cutregex3l);
          if (not res3l or m3l.empty() or (m3l.size() != 3)) {
            LOGF(fatal, "PIDSelectionFilterAndAnalysis::::ConstructCutFromString. Wrong RE for detector %s: %s, try tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}} for instance", detector.c_str(), cut.c_str());
          }
          LOGF(info, "Captured %s", m3l[1].str().c_str());
          /* let's search for the involved species index */
          int spindex = -1;
          for (uint i = 0; i < mSpeciesNames.size(); ++i) {
            TString spname = detector + mSpeciesNames[i];
            LOGF(info, "Checking %s vs %s", spname.Data(), m3l[1].str().c_str());
            if (spname.EqualTo(m3l[1].str().c_str())) {
              spindex = i;
              break;
            }
          }
          if (spindex < 0) {
            LOGF(fatal, "PIDSelectionFilterAndAnalysis::::ConstructCutFromString. Wrong species %s for detector %s: %s, try tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}} for instance", m3l[1].str().c_str(), detector.c_str(), cut.c_str());
          } else {
            storeCut(bricklst[spindex]);
          }
          /* removing the already handled species requirements */
          lev3str.Remove(0, m3l[0].length());
        }
      };
      if (m[1].str() == "tpcsel") {
        handleDetectorLevel("tpc", mCloseNsigmasTPC, m[2].str());
      } else if (m[1].str() == "tofsel") {
        handleDetectorLevel("tof", mCloseNsigmasTOF, m[2].str());
      } else if (m[1].str() == "baysel") {
        handleDetectorLevel("bay", mBayesProbability, m[2].str());
      } else {
        Fatal("PIDSelectionFilterAndAnalysis::::ConstructCutFromString", "Wrong RE detector %s, use tpcsel, or tofsel only", m[1].str().c_str());
      }
      /* removing the already handled detector requirements */
      lev2str.Remove(0, m[0].length());
    }
  }
  mMaskLength = CalculateMaskLength();
}

void PIDSelectionFilterAndAnalysis::StoreArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask.clear();
  mForcedArmedMask = 0UL;
  int bit = 0;
  auto armedList = [&](auto bricklst) {
    auto armedBrick = [&](auto brick, bool opt = false) {
      if (brick != nullptr) {
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
      }
    };
    for (auto brick : bricklst) {
      armedBrick(brick, false);
    }
  };
  armedList(mCloseNsigmasTPC);
  armedList(mCloseNsigmasTOF);
  armedList(mBayesProbability);

  LOGF(info, "PIDSelectionFilterAndAnalysis::StoreArmedMask(), masks 0x%016lx, %s, 0x%016lx", armedMask, printOptionalMasks().Data(), forcedMask);
  mArmedMask = armedMask;
  mForcedArmedMask = forcedMask;
}
