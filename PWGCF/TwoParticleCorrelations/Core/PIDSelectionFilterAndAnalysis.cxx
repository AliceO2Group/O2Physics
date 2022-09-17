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
    mInclusiveTrackPID(0),
    mExclusiveTrackPID(0)
{
}

/// \brief Constructor from regular expression
PIDSelectionFilterAndAnalysis::PIDSelectionFilterAndAnalysis(const TString& cutstr, selmodes mode)
  : SelectionFilterAndAnalysis("", mode),
    mPTOF(0.8f),
    mRequireTOF(false),
    mEllipticTPCTOF(false),
    mInclusiveTrackPID(0),
    mExclusiveTrackPID(0)
{
  ConstructCutFromString(cutstr);
}

/// \brief Calculates the length of the mask needed to store the selection cuts
int PIDSelectionFilterAndAnalysis::CalculateMaskLength()
{
  int length = 0;
  for (auto brick : mInclusiveTrackPID) {
    length += brick->Length();
  }
  for (auto brick : mExclusiveTrackPID) {
    length += brick->Length();
  }
  return length;
}

void PIDSelectionFilterAndAnalysis::ConstructCutFromString(const TString& cutstr)
{
  mMaskLength = CalculateMaskLength();
  if (mMaskLength > 64) {
    LOGF(fatal, "PIDSelectionFilterAndAnalysis not ready for filter mask of %d bits. Just 64 available for the time being", mMaskLength);
  }
}
