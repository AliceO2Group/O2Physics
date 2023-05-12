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
#include "SelectionFilterAndAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;

ClassImp(SelectionFilterAndAnalysis);

/// \brief Default constructor
SelectionFilterAndAnalysis::SelectionFilterAndAnalysis()
  : TNamed(),
    mMode(kFilter),
    mCutStringSignature(),
    mMaskLength(0),
    mSelectedMask(0UL),
    mArmedMask(0UL),
    mOptArmedMask{},
    mForcedArmedMask(0UL)
{
}

/// \brief Named constructor for a concrete operating mode
SelectionFilterAndAnalysis::SelectionFilterAndAnalysis(const TString& name, selmodes mode)
  : TNamed(name, name),
    mMode(mode),
    mCutStringSignature(),
    mMaskLength(0),
    mSelectedMask(0UL),
    mArmedMask(0UL),
    mOptArmedMask{},
    mForcedArmedMask(0UL)
{
}

TString SelectionFilterAndAnalysis::printOptionalMasks() const
{
  TString str = "(";
  bool first = true;
  for (auto option : mOptArmedMask) {
    if (not first) {
      str += ", ";
    }
    first = false;
    str += TString::Format("0x%016lx", u_long(option));
  }
  str += ")";
  return str;
}
