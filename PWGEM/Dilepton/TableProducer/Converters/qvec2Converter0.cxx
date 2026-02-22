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
// This code converts q vector table into qvec2 table.
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

struct qvec2Converter0 {
  Produces<aod::EMEventsQvec2_000> qvec2_000;

  void process(aod::EMEventsQvec_001 const& collisions)
  {
    for (const auto& collision : collisions) {
      qvec2_000(
        collision.q2xft0m(), collision.q2yft0m(),
        collision.q2xft0a(), collision.q2yft0a(),
        collision.q2xft0c(), collision.q2yft0c(),
        collision.q2xfv0a(), collision.q2yfv0a(),
        collision.q2xbpos(), collision.q2ybpos(),
        collision.q2xbneg(), collision.q2ybneg(),
        collision.q2xbtot(), collision.q2ybtot());
    } // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qvec2Converter0>(cfgc, TaskName{"qvec2-converter0"})};
}
