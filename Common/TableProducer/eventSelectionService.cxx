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

/// \file eventSelectionTester.cxx
/// \brief unified, self-configuring event selection task
/// \author ALICE

//===============================================================
//
// Unified, self-configuring event selection task
//
//===============================================================

#include "MetadataHelper.h"

#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/EventSelectionTools.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonUtils/NameConf.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::framework;

MetadataHelper metadataInfo; // Metadata helper

using BCsWithRun2InfosTimestampsAndMatches = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::Run2MatchedToBCSparse>;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;

struct eventselectionRun2 {
  o2::common::eventselection::bcselConfigurables bcselOpts;
  o2::common::eventselection::BcSelectionModule bcselmodule;

  o2::common::eventselection::evselConfigurables evselOpts;
  o2::common::eventselection::EventSelectionModule evselmodule;

  Produces<aod::BcSels> bcsel;
  Produces<aod::EvSels> evsel;

  // for slicing
  SliceCache cache;

  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // the best: have readable cursors
  // this: a stopgap solution to avoid spawning yet another device
  std::vector<o2::common::eventselection::bcselEntry> bcselsbuffer;

  // auxiliary
  Partition<FullTracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Preslice<FullTracks> perCollision = aod::track::collisionId;

  void init(o2::framework::InitContext& context)
  {
    // CCDB boilerplate init
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);

    // task-specific
    bcselmodule.init(context, bcselOpts, histos);
    evselmodule.init(context, evselOpts, histos, metadataInfo);
  }

  void process(BCsWithRun2InfosTimestampsAndMatches const& bcs,
               aod::Collisions const& collisions,
               aod::Zdcs const&,
               aod::FV0As const&,
               aod::FV0Cs const&,
               aod::FT0s const&,
               aod::FDDs const&,
               FullTracks const&)
  {
    bcselmodule.processRun2(ccdb, bcs, bcselsbuffer, bcsel);
    evselmodule.processRun2(ccdb, histos, collisions, tracklets, cache, bcselsbuffer, evsel);
  }
};

struct eventselectionRun3 {
  o2::common::eventselection::bcselConfigurables bcselOpts;
  o2::common::eventselection::BcSelectionModule bcselmodule;

  o2::common::eventselection::evselConfigurables evselOpts;
  o2::common::eventselection::EventSelectionModule evselmodule;

  o2::common::eventselection::lumiConfigurables lumiOpts;
  o2::common::eventselection::LumiModule lumimodule;

  Produces<aod::BcSels> bcsel;
  Produces<aod::EvSels> evsel;

  // for slicing
  SliceCache cache;

  // CCDB boilerplate declarations
  o2::framework::Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // the best: have readable cursors
  // this: a stopgap solution to avoid spawning yet another device
  std::vector<o2::common::eventselection::bcselEntry> bcselsbuffer;

  // auxiliary
  Partition<FullTracksIU> pvTracks = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  Preslice<FullTracksIU> perCollisionIU = aod::track::collisionId;

  void init(o2::framework::InitContext& context)
  {
    // CCDB boilerplate init
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setURL(ccdburl.value);

    // task-specific
    bcselmodule.init(context, bcselOpts, histos);
    evselmodule.init(context, evselOpts, histos, metadataInfo);
    lumimodule.init(context, lumiOpts, histos);
  }

  void process(aod::Collisions const& collisions,
               BCsWithRun3Matchings const& bcs,
               aod::Zdcs const&,
               aod::FV0As const&,
               aod::FT0s const& ft0s, // to resolve iterator
               aod::FDDs const&,
               FullTracksIU const&)
  {
    bcselmodule.processRun3(ccdb, histos, bcs, bcselsbuffer, bcsel);
    evselmodule.processRun3(ccdb, histos, bcs, collisions, pvTracks, ft0s, cache, bcselsbuffer, evsel);
    lumimodule.process(ccdb, histos, bcs, bcselsbuffer);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata for later too
  metadataInfo.initMetadata(cfgc);

  bool isRun3 = true, hasRunInfo = false;
  if (cfgc.options().hasOption("aod-metadata-Run") == true) {
    hasRunInfo = true;
    if (cfgc.options().get<std::string>("aod-metadata-Run") == "2") {
      isRun3 = false;
    }
  }

  LOGF(info, "Event selection autoconfiguring from metadata. Availability of info for Run 2/3 is %i", hasRunInfo);
  if (!hasRunInfo) {
    LOGF(info, "Metadata info missing or incomplete. Make sure --aod-file is provided at the end of the last workflow and that the AO2D has metadata stored.");
    LOGF(info, "Initializing with Run 3 data as default. Please note you will not be able to change settings manually.");
    LOGF(info, "You should instead make sure the metadata is read in correctly.");
    return WorkflowSpec{adaptAnalysisTask<eventselectionRun3>(cfgc)};
  } else {
    LOGF(info, "Metadata successfully read in. Is this Run 3? %i - will self-configure.", isRun3);
    if (isRun3) {
      return WorkflowSpec{adaptAnalysisTask<eventselectionRun3>(cfgc)};
    } else {
      return WorkflowSpec{adaptAnalysisTask<eventselectionRun2>(cfgc)};
    }
  }
  throw std::runtime_error("Unsupported run type / problem when configuring event selection!");
}
