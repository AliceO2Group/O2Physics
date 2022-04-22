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
#include "Framework/ConfigParamSpec.h"

using namespace o2;
using namespace o2::framework;

// custom configurable for switching between run2 and run3 selection types
void customize(std::vector<ConfigParamSpec>& workflowOptions)
{
  workflowOptions.push_back(ConfigParamSpec{"selection-run", VariantType::Int, 2, {"selection type: 2 - run 2, 3 - run 3"}});
}

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "iostream"

struct MultiplicityTableTaskIndexed {
  Produces<aod::Mults> mult;
  
  //For vertex-Z corrections in calibration
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<soa::Join<aod::Tracks, aod::TracksExtra>> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  
  TList fCalibObjects; //for storing calibration objects
  
  void init(InitContext& context)
  {
    ccdb->setURL("Users/v/victor/Centrality/Calibration"); //temporary - to be tuned  shortly
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }
  
  void processRun2(aod::Run2MatchedSparse::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracksExtra, aod::BCs const&, aod::Zdcs const&, aod::FV0As const& fv0as, aod::FV0Cs const& fv0cs, aod::FT0s const& ft0s)
  {
    float multV0A = 0.f;
    float multV0C = 0.f;
    float multT0A = 0.f;
    float multT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = 0.f;
    float multZNC = 0.f;
    
    float multZeqV0A = 0.f;
    float multZeqT0A = 0.f;
    float multZeqT0C = 0.f;
    float multZeqFDDA = 0.f;
    float multZeqFDDC = 0.f;
    float multZeqNContribs = 0.f;
    
    auto trackletsGrouped = run2tracklets->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex());
    int multTracklets = trackletsGrouped.size();
    int multTPC = tracksGrouped.size();
    int multNContribs = 0;

    if (collision.has_fv0a()) {
      for (auto amplitude : collision.fv0a().amplitude()) {
        multV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (auto amplitude : collision.fv0c().amplitude()) {
        multV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      multZNA = zdc.energyCommonZNA();
      multZNC = zdc.energyCommonZNC();
    }

    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multT0A=%5.0f multT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multV0A, multV0C, multT0A, multT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC);
    mult(multV0A, multV0C, multT0A, multT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC, multNContribs,
         multZeqV0A, multZeqT0A, multZeqT0A, multZeqFDDA, multZeqFDDC, multZeqNContribs);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun2, "Produce Run 2 multiplicity tables", true);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra> const& tracksExtra, aod::BCs const& bcs, aod::Zdcs const& zdcs, aod::FV0As const& fv0as, aod::FT0s const& ft0s, aod::FDDs const& fdds)
  {
    float multV0A = 0.f;
    float multV0C = 0.f;
    float multT0A = 0.f;
    float multT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = 0.f;
    float multZNC = 0.f;
    int multTracklets = 0;

    float multZeqV0A = 0.f;
    float multZeqT0A = 0.f;
    float multZeqT0C = 0.f;
    float multZeqFDDA = 0.f;
    float multZeqFDDC = 0.f;
    float multZeqNContribs = 0.f;
    
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto pvContribsGrouped = pvContribTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    int multTPC = tracksGrouped.size();
    int multNContribs = tracksGrouped.size();
    
    TList *lCalibObjects = ccdb->getForTimeStamp<TList>("Users/v/victor/Centrality/Calibration", 1635634560883); //temporary
    
    TProfile *hVtxZV0A     = (TProfile*) lCalibObjects->FindObject("hVtxZV0A");
    TProfile *hVtxZT0A     = (TProfile*) lCalibObjects->FindObject("hVtxZT0A");
    TProfile *hVtxZT0C     = (TProfile*) lCalibObjects->FindObject("hVtxZT0C");
    TProfile *hVtxZFDDA    = (TProfile*) lCalibObjects->FindObject("hVtxZFDDA");
    TProfile *hVtxZFDDC    = (TProfile*) lCalibObjects->FindObject("hVtxZFDDC");
    TProfile *hVtxZNTracks = (TProfile*) lCalibObjects->FindObject("hVtxZNTracks");
    
    if(fabs(collision.posZ())<15.0f){
      multZeqV0A       = multV0A/hVtxZV0A->Interpolate( collision.posZ() );
      multZeqT0A       = multT0A/hVtxZT0A->Interpolate( collision.posZ() );
      multZeqT0C       = multT0C/hVtxZT0C->Interpolate( collision.posZ() );
      multZeqFDDA      = multFDDA/hVtxZFDDA->Interpolate( collision.posZ() );
      multZeqFDDC      = multFDDC/hVtxZFDDC->Interpolate( collision.posZ() );
      multZeqNContribs = multNContribs/hVtxZNTracks->Interpolate( collision.posZ() );
    }
    
    // using FT0 row index from event selection task
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (auto amplitude : ft0.amplitudeA()) {
        multT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multT0C += amplitude;
      }
    }
    // using FDD row index from event selection task
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      for (auto amplitude : fdd.amplitudeA()) {
        multFDDA += amplitude;
      }
      for (auto amplitude : fdd.amplitudeC()) {
        multFDDC += amplitude;
      }
    }
    // using FV0 row index from event selection task
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      for (auto amplitude : fv0.amplitude()) {
        multV0A += amplitude;
      }
    }
    LOGF(debug, "multV0A=%5.0f multV0C=%5.0f multT0A=%5.0f multT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multV0A, multV0C, multT0A, multT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC);
    mult(multV0A, multV0C, multT0A, multT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC, multNContribs,
         multZeqV0A, multZeqT0A, multZeqT0A, multZeqFDDA, multZeqFDDC, multZeqNContribs);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun3, "Produce Run 3 multiplicity tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityTableTaskIndexed>(cfgc, TaskName{"multiplicity-table"})};
}
