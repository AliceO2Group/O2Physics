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
///
/// \brief A tutorial task to retrieve objects from CCDB given a time stamp.
/// \author Daiki Sekihata
/// \since 2026-03-01

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <chrono>

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

namespace o2::aod
{
namespace testccdb
{
DECLARE_SOA_CCDB_COLUMN(GRPMagField, grpMagField, o2::parameters::GRPMagField, "GLO/Config/GRPMagField");   //!
DECLARE_SOA_CCDB_COLUMN(MeanVertex, meanVertex, o2::dataformats::MeanVertexObject, "GLO/Calib/MeanVertex"); //!
} // namespace testccdb

DECLARE_SOA_TIMESTAMPED_TABLE(MyCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "MYCCDBOBJ", //!
                              testccdb::GRPMagField, testccdb::MeanVertex);
} // namespace o2::aod

struct TestCCDBTable {
  void init(o2::framework::InitContext&) {}

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::MyCCDBObjects>;

  void process(MyBCs const& bcs)
  {
    int i = 0;
    for (const auto& bc : bcs) {
      if (i >= 5) {
        return;
      }
      float l3current = bc.grpMagField().getL3Current();
      float zvtx = bc.meanVertex().getZ();
      LOGF(info, "bc.globalBC() = %llu, bc.timestamp() = %llu, L3 current = %f A, mean zvtx = %f cm", bc.globalBC(), bc.timestamp(), l3current, zvtx);
      i++;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TestCCDBTable>(cfgc),
  };
}
