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
/// \file   timestamp.cxx
/// \author Nicolò Jacazio
/// \since  2020-06-22
/// \brief  A task to fill the timestamp table from run number.
///         Uses headers from CCDB
///
#include "MetadataHelper.h"

#include "Common/Tools/timestampModule.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DetectorsRaw/HBFUtils.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <vector>

using namespace o2::framework;
using namespace o2::header;
using namespace o2;

o2::common::core::MetadataHelper metadataInfo; // Metadata helper

struct TimestampTask {
  Produces<aod::Timestamps> timestampTable; /// Table with SOR timestamps produced by the task
  Service<o2::ccdb::BasicCCDBManager> ccdb; /// CCDB manager to access orbit-reset timestamp
  o2::ccdb::CcdbApi ccdb_api;               /// API to access CCDB headers

  Configurable<std::string> ccdb_url{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB database"};

  o2::common::timestamp::timestampConfigurables timestampConfigurables;
  o2::common::timestamp::TimestampModule timestampMod;

  std::vector<int64_t> timestampBuffer;

  void init(o2::framework::InitContext&)
  {
    // CCDB initialization
    ccdb->setURL(ccdb_url.value);
    ccdb_api.init(ccdb_url.value);

    // timestamp configuration + init
    timestampMod.init(timestampConfigurables, metadataInfo);
  }

  void process(aod::BCs const& bcs)
  {
    timestampMod.process(bcs, ccdb, timestampBuffer, timestampTable);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{adaptAnalysisTask<TimestampTask>(cfgc)};
}
