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
#include "FilterAndAnalysisFramework.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;
using namespace o2::framework::expressions;
using namespace o2::analysis::PWGCF;

ClassImp(FilterAndAnalysisFramework);

/// \brief Default constructor
FilterAndAnalysisFramework::FilterAndAnalysisFramework()
  : TNamed(),
    ccdb(nullptr),
    fTrackFilter(nullptr),
    fEventFilter(nullptr),
    fPIDFilter(nullptr)
{
}

/// \brief Named constructor for a concrete operating mode
FilterAndAnalysisFramework::FilterAndAnalysisFramework(EventSelectionConfigurable& evtf, TrackSelectionConfigurable& trkf, PIDSelectionConfigurable& pidf, SelectionFilterAndAnalysis::selmodes mode)
  : TNamed(),
    ccdb(nullptr),
    fTrackFilter(nullptr),
    fEventFilter(nullptr),
    fPIDFilter(nullptr)
{
  fTrackFilter = new PWGCF::TrackSelectionFilterAndAnalysis(trkf, mode);
  fEventFilter = new PWGCF::EventSelectionFilterAndAnalysis(evtf, mode);
  fPIDFilter = new PWGCF::PIDSelectionFilterAndAnalysis(pidf, mode);
}
