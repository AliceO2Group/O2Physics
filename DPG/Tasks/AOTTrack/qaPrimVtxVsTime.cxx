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
/// \author Mattia Faggin <mattia.faggin@cern.ch>, Padova University and INFN

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "CCDB/BasicCCDBManager.h"
#include "Framework/HistogramRegistry.h"
#include "CommonDataFormat/BunchFilling.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2::framework;
using namespace o2;

struct QaPrimVtxVsTime {

  int lastRunNumber = -1;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  ConfigurableAxis binsVertexPosZ{"binsVertexPosZ", {100, -20., 20.}, ""};

  AxisSpec axisVertexPosZ{binsVertexPosZ, "Primary vertex Z (cm)"};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// @brief init function
  /// @param
  void init(InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  using BCsWithTimeStamp = soa::Join<aod::BCs, aod::Timestamps>;
  /// @brief process function
  void process(aod::Collisions const& collisions, BCsWithTimeStamp const& bcs)
  {
    int runNumber = bcs.rawIteratorAt(0).runNumber();
    if (lastRunNumber != runNumber) {
      /// execute the code in this scope only once, i.e. when the current run is considered for the first time in this DF
      lastRunNumber = runNumber;
      int64_t tsSOR = 0;
      int64_t tsEOR = 0;

      /// reject AO2Ds for which no CCDB access is possible
      if (runNumber < 500000) {
        LOG(warning) << ">>> run number " << runNumber << " < 500000. access to CCDB not possible. Exiting";
        return;
      }

      /// If we are here, the current run was never considered before.
      /// Let's add the TH2 that we need for the monitoring
      /// Let's define the x-axis according to the start-of-run (SOR) and end-of-run (EOR) times
      o2::ccdb::CcdbApi ccdb_api;
      ccdb_api.init(ccdburl);

      auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdb_api, runNumber);
      tsSOR = soreor.first;
      tsEOR = soreor.second;

      double minSec = floor(tsSOR / 1000.); /// round tsSOR to the highest integer lower than tsSOR
      double maxSec = ceil(tsEOR / 1000.);  /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>(maxSec - minSec), minSec, maxSec, "seconds (from January 1st, 1970 at UTC)"};
      histos.add("hPosZvsTime", "", kTH2F, {axisSeconds, axisVertexPosZ});
    }

    /// The rest of the code is always run
    /// Here we fill the correlation TH2 posZ vs. time
    for (auto& collision : collisions) {
      /// skip collisions w/o an assigned bunch crossing
      if (!collision.has_bc()) {
        continue;
      }

      const auto timestamp = collision.bc_as<BCsWithTimeStamp>().timestamp(); /// NB: in ms
      histos.fill(HIST("hPosZvsTime"), timestamp / 1000., collision.posZ());
    }
  }
};
//__________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<QaPrimVtxVsTime>(cfgc)};
}
