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
/// \file   qaPIDTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Implementation for QA tasks of the TOF PID quantities
///

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the TOF QA plots
struct tofSimple {
  Service<o2::pid::tof::TOFResponse> tofResponse;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(o2::framework::InitContext& initContext)
  {
    LOG(info) << "tofSimple >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ";
    tofResponse->initSetup(ccdb, initContext);
  }
  void process(aod::BCs const&) {}
};

struct tofSimpleASD {
  Service<o2::pid::tof::TOFResponse> tofResponse;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  void init(o2::framework::InitContext& initContext)
  {
    LOG(info) << "tofSimpleASD >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ";
    tofResponse->initSetup(ccdb, initContext);
    if (!tofResponse->isInit()) {
      LOG(fatal) << "TOF response not initialized";
    } else {
      LOG(info) << "tofSimpleASD OK >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ";
    }
  }
  void process(aod::BCs const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  o2::pid::tof::TOFResponseImpl::metadataInfo.initMetadata(cfgc);
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSimple>(cfgc)};
  workflow.push_back(adaptAnalysisTask<tofSimpleASD>(cfgc));
  return workflow;
}
