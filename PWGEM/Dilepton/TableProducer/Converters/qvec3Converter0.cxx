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
//
// ========================
//
// This code converts q vector table into qvec3 table.
//    Please write to: daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/DataModel/dileptonTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct qvec3Converter0 {
  Produces<aod::EMEventsQvec3_000> qvec3_000;

  void process(aod::EMEventsQvec_001 const& collisions)
  {
    for (const auto& collision : collisions) {
      qvec3_000(
          collision.q3xft0m(), collision.q3yft0m(),
          collision.q3xft0a(), collision.q3yft0a(),
          collision.q3xft0c(), collision.q3yft0c(),
          collision.q3xfv0a(), collision.q3yfv0a(),
          collision.q3xbpos(), collision.q3ybpos(),
          collision.q3xbneg(), collision.q3ybneg(),
          collision.q3xbtot(), collision.q3ybtot());
    } // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qvec3Converter0>(cfgc, TaskName{"qvec3-converter0"})};
}
