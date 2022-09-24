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
#include "PIDSelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;
using namespace boost;

ClassImp(PIDSelectionFilterAndAnalysis);

/// \brief Default constructor
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis()
  : SelectionFilterAndAnalysis(),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasITS(kNoOfSpecies, nullptr),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr)
{
}

/// \brief Constructor from regular expression
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const TString& cutstr, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasITS(kNoOfSpecies, nullptr),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr)
{
  ConstructCutFromString(cutstr);
}

/// \brief Constructor from the PID selection configurable
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const PIDSelectionConfigurable& pidsel, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mCloseNsigmasITS(kNoOfSpecies, nullptr),
    mCloseNsigmasTPC(kNoOfSpecies, nullptr),
    mCloseNsigmasTOF(kNoOfSpecies, nullptr)
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
  /* TODO: interlace the different detectors */
  appendCut(pidsel.mPidItsSel_el);
  appendCut(pidsel.mPidItsSel_mu);
  appendCut(pidsel.mPidItsSel_pi);
  appendCut(pidsel.mPidItsSel_ka);
  appendCut(pidsel.mPidItsSel_pr);
  appendCut(pidsel.mPidTpcSel_el);
  appendCut(pidsel.mPidTpcSel_mu);
  appendCut(pidsel.mPidTpcSel_pi);
  appendCut(pidsel.mPidTpcSel_ka);
  appendCut(pidsel.mPidTpcSel_pr);
  appendCut(pidsel.mPidTofSel_el);
  appendCut(pidsel.mPidTofSel_mu);
  appendCut(pidsel.mPidTofSel_pi);
  appendCut(pidsel.mPidTofSel_ka);
  appendCut(pidsel.mPidTofSel_pr);

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

/// \brief Destructor
/// Releases the taken memory
PIDSelectionFilterAndAnalysis::~PIDSelectionFilterAndAnalysis()
{
  for (auto brick : mCloseNsigmasITS) {
    delete brick;
  }
  for (auto brick : mCloseNsigmasTPC) {
    delete brick;
  }
  for (auto brick : mCloseNsigmasTOF) {
    delete brick;
  }
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int PIDSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
  auto addLength = [&](auto bricklst) {
    for (auto brick : bricklst) {
      length += brick->Length();
    }
  };
  addLength(mCloseNsigmasITS);
  addLength(mCloseNsigmasTPC);
  addLength(mCloseNsigmasTOF);
  return length;
}

void PIDSelectionFilterAndAnalysis::ConstructCutFromString(const TString& cutstr)
{
  mMaskLength = CalculateMaskLength();
}

void PIDSelectionFilterAndAnalysis::StoreArmedMask()
{
  uint64_t armedMask = 0UL;
  uint64_t optMask = 0UL;
  uint64_t forcedMask = 0UL;
  mArmedMask = 0UL;
  mOptArmedMask = 0UL;
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
      armedBrick(brick, true);
    }
  };
  armedList(mCloseNsigmasITS);
  armedList(mCloseNsigmasTPC);
  armedList(mCloseNsigmasTOF);

  LOGF(info, "PIDSelectionFilterAndAnalysis::StoreArmedMask(), masks 0x%08lx, 0x%08lx, 0x%08lx", armedMask, optMask, forcedMask);
  mArmedMask = armedMask;
  mOptArmedMask = optMask;
  mForcedArmedMask = forcedMask;
}
