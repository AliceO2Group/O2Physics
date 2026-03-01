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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CommonDataFormat/InteractionRecord.h"

#include <chrono>

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

namespace o2::aod
{
namespace testccdb
{
DECLARE_SOA_CCDB_COLUMN(GRPObject, grpObject, o2::parameters::GRPObject, "GLO/GRP/GRP"); //!
DECLARE_SOA_CCDB_COLUMN(GRPMagField, grpMagField, o2::parameters::GRPMagField, "GLO/Config/GRPMagField"); //!
} // namespace testccdb

DECLARE_SOA_TIMESTAMPED_TABLE(MyCCDBObjects, aod::Timestamps, o2::aod::timestamp::Timestamp, 1, "MYCCDBOBJ", //!
                              testccdb::GRPObject, testccdb::GRPMagField);
} // namespace o2::aod

struct TestCCDBTable {

  void init(o2::framework::InitContext&) {}

  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::MyCCDBObjects>;
  // using MyBCs = soa::Join<aod::BCs, aod::MyCCDBObjects>;

  void process(MyBCs const& bcs)
  {
    for (const auto& bc: bcs) {
      float bz = bc.grpObject().getNominalL3Field();
      float l3current = bc.grpMagField().getL3Current();
      LOGF(info, "bc.timestamp() = %lld, bz = %f kG, %f A", bc.timestamp(), bz, l3current);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TestCCDBTable>(cfgc),
  };
}
