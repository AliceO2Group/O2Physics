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

/// \file qvecConverter.cxx
/// \brief Analysis task for neutral pion flow with EMCal
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Converts EMEventsQvec_000 into EMEventsQvec_001
struct QvecConverter {
  Produces<aod::EMEventsQvec_001> qvec001;

  void process(aod::EMEventsQvec_000 const& emEventsQVec)
  {
    for (const auto& qvec : emEventsQVec) {
      constexpr float EmptyV0 = -999.f;
      qvec001(qvec.q2xft0m(), qvec.q2yft0m(), qvec.q2xft0a(), qvec.q2yft0a(), qvec.q2xft0c(), qvec.q2yft0c(), EmptyV0, EmptyV0, qvec.q2xbpos(), qvec.q2ybpos(), qvec.q2xbneg(), qvec.q2ybneg(), qvec.q2xbtot(), qvec.q2ybtot(), qvec.q3xft0m(), qvec.q3yft0m(), qvec.q3xft0a(), qvec.q3yft0a(), qvec.q3xft0c(), qvec.q3yft0c(), EmptyV0, EmptyV0, qvec.q3xbpos(), qvec.q3ybpos(), qvec.q3xbneg(), qvec.q3ybneg(), qvec.q3xbtot(), qvec.q3ybtot());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QvecConverter>(cfgc),
  };
}
