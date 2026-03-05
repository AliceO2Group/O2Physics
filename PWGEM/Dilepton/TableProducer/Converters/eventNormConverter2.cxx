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
// This code runs loop over ULS ee pars for virtual photon QC.
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

struct eventNormConverter2 {
  Produces<aod::EMEventNormInfos_002> event_002;

  void process(aod::EMEventNormInfos_001 const& collisions)
  {
    for (const auto& collision : collisions) {

      int8_t posZint8 = static_cast<int8_t>(collision.posZ() * 2.f);
      if (posZint8 == 0) {
        if (collision.posZ() < 0.f) {
          posZint8 = -1;
        } else {
          posZint8 = +1;
        }
      }

      event_002(
        o2::aod::emevsel::reduceSelectionBit(collision),
        collision.rct_raw(),
        posZint8,
        105 + 110,
        collision.centFT0C() < 1.f ? static_cast<uint8_t>(collision.centFT0C() * 100.f) : static_cast<uint8_t>(collision.centFT0C() + 110.f),
        105 + 110);
    } // end of collision loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<eventNormConverter2>(cfgc, TaskName{"event-norm-converter2"})};
}
