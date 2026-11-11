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

/// \file PIDSelectionFilterAndAnalysis.cxx
/// \brief Implementation of the PIDSelectionFilterAndAnalysis class for
/// PID filtering tracks for derived data production
/// \author victor.gonzalez.sebastian@gmail.com

#include "PIDSelectionFilterAndAnalysis.h"

#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SkimmingConfigurableCuts.h"

#include <Framework/Logger.h>

#include <TString.h>

#include <boost/regex/v5/regbase.hpp>
#include <boost/regex/v5/regex.hpp>
#include <boost/regex/v5/regex_fwd.hpp>
#include <boost/regex/v5/regex_search.hpp>
#include <sys/types.h>

#include <Rtypes.h>

#include <cstdint>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::analysis::pwgcf;
using namespace boost;

#ifdef INCORPORATEBAYESIANPID
PIDSelectionConfigurable::PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                                                   std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr,
                                                   std::vector<std::string> pidbayel, std::vector<std::string> pidbaymu, std::vector<std::string> pidbaypi, std::vector<std::string> pidbayka, std::vector<std::string> pidbaypr)
#else
PIDSelectionConfigurable::PIDSelectionConfigurable(std::vector<std::string> pidtpcel, std::vector<std::string> pidtpcmu, std::vector<std::string> pidtpcpi, std::vector<std::string> pidtpcka, std::vector<std::string> pidtpcpr,
                                                   std::vector<std::string> pidtofel, std::vector<std::string> pidtofmu, std::vector<std::string> pidtofpi, std::vector<std::string> pidtofka, std::vector<std::string> pidtofpr)
#endif
  : mPidTpcSelEl{},
    mPidTpcSelMu{},
    mPidTpcSelPi{},
    mPidTpcSelKa{},
    mPidTpcSelPr{},
    mPidTofSelEl{},
    mPidTofSelMu{},
    mPidTofSelPi{},
    mPidTofSelKa{},
    mPidTofSelPr{},
    mPidBaySelEl{},
    mPidBaySelMu{},
    mPidBaySelPi{},
    mPidBaySelKa{},
    mPidBaySelPr{}
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
  mPidTpcSelEl = storeCutString(pidtpcel, "tpcel");
  mPidTpcSelMu = storeCutString(pidtpcmu, "tpcmu");
  mPidTpcSelPi = storeCutString(pidtpcpi, "tpcpi");
  mPidTpcSelKa = storeCutString(pidtpcka, "tpcka");
  mPidTpcSelPr = storeCutString(pidtpcpr, "tpcpr");
  mPidTofSelEl = storeCutString(pidtofel, "tofel");
  mPidTofSelMu = storeCutString(pidtofmu, "tofmu");
  mPidTofSelPi = storeCutString(pidtofpi, "tofpi");
  mPidTofSelKa = storeCutString(pidtofka, "tofka");
  mPidTofSelPr = storeCutString(pidtofpr, "tofpr");
#ifdef INCORPORATEBAYESIANPID
  mPidBaySelEl = storeCutString(pidbayel, "bayel");
  mPidBaySelMu = storeCutString(pidbaymu, "baymu");
  mPidBaySelPi = storeCutString(pidbaypi, "baypi");
  mPidBaySelKa = storeCutString(pidbayka, "bayka");
  mPidBaySelPr = storeCutString(pidbaypr, "baypr");
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
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const TString& cutstr, Selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr),
    mBayesProbability(kNoOfSpecies, nullptr)
{
  constructCutFromString(cutstr);
}

/// \brief Constructor from the PID selection configurable
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const PIDSelectionConfigurable& pidsel, Selmodes mode)
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
  constexpr int MaxLength = 64;

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
  cutString += "tpcsel{";
  first = true;
  appendCut(pidsel.mPidTpcSelEl);
  appendCut(pidsel.mPidTpcSelMu);
  appendCut(pidsel.mPidTpcSelPi);
  appendCut(pidsel.mPidTpcSelKa);
  appendCut(pidsel.mPidTpcSelPr);
  cutString += "},tofsel{";
  first = true;
  appendCut(pidsel.mPidTofSelEl);
  appendCut(pidsel.mPidTofSelMu);
  appendCut(pidsel.mPidTofSelPi);
  appendCut(pidsel.mPidTofSelKa);
  appendCut(pidsel.mPidTofSelPr);
  cutString += "},baysel{";
  first = true;
  appendCut(pidsel.mPidBaySelEl);
  appendCut(pidsel.mPidBaySelMu);
  appendCut(pidsel.mPidBaySelPi);
  appendCut(pidsel.mPidBaySelKa);
  appendCut(pidsel.mPidBaySelPr);

  cutString += "}}";
  constructCutFromString(cutString);
  if (mMaskLength > MaxLength) {
    LOGF(fatal, "EventSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just %d available for the time being", mMaskLength, MaxLength);
  }
  mCutStringSignature = cutString.ReplaceAll("-yes", "").ReplaceAll("-no", "");
  if (mode == KAnalysis) {
    storeArmedMask();
  }
}

/// \brief Destructor
/// Releases the taken memory
PIDSelectionFilterAndAnalysis::~PIDSelectionFilterAndAnalysis()
{
  for (const auto& brick : mCloseNsigmasTPC) {
    delete brick;
  }
  for (const auto& brick : mCloseNsigmasTOF) {
    delete brick;
  }
  for (const auto& brick : mBayesProbability) {
    delete brick;
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int PIDSelectionFilterAndAnalysis::calculateMaskLength()
{
  int length = 0;
  auto addLength = [&](const auto& bricklst) {
    for (const auto& brick : bricklst) {
      if (brick != nullptr) {
        length += brick->length();
      }
    }
  };
  addLength(mCloseNsigmasTPC);
  addLength(mCloseNsigmasTOF);
  addLength(mBayesProbability);
  return length;
}

void PIDSelectionFilterAndAnalysis::constructCutFromString(const TString& cutstr)
{
  LOGF(info, "Cut string: %s", cutstr.Data());
  /* let's catch the first level */
  regex cutregex("^pidsel\\{?(\\w+\\{[\\w\\d.,:{}-]+})*}$", regex_constants::ECMAScript | regex_constants::icase);
  std::string in(cutstr.Data());
  smatch m;

  bool res = regex_search(in, m, cutregex);
  constexpr int MaxComponents = 2;
  if (not res or m.empty() or (m.size() > MaxComponents)) {
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
      constexpr int ExpectedComponents = 3;
      if (not res or m.empty() or (m.size() != ExpectedComponents)) {
        LOGF(fatal, "PIDSelectionFilterAndAnalysis::::ConstructCutFromString. Wrong RE: %s, try pidsel{tpcsel{tpcpi{cwv{rg{-3.0,3.0}:rg{-2.0,2.0},rg{-3.0,5.0}}}}} for instance", cutstr.Data());
      }
      LOGF(info, "Captured %s", m[1].str().c_str());
      auto handleDetectorLevel = [](const std::string& detector, auto& bricklst, const std::string& cut) {
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
          if (not res3l or m3l.empty() or (m3l.size() != ExpectedComponents)) {
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
  mMaskLength = calculateMaskLength();
}

void PIDSelectionFilterAndAnalysis::storeArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask.clear();
  mForcedArmedMask = 0UL;
  int bit = 0;
  auto armedList = [&](const auto& bricklst) {
    auto armedBrick = [&](auto brick, bool opt = false) {
      if (brick != nullptr) {
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
      }
    };
    for (const auto& brick : bricklst) {
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
