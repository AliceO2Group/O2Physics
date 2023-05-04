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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include <CCDB/BasicCCDBManager.h>
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "iostream"

struct MultiplicityTableTaskIndexed {
  SliceCache cache;
  Produces<aod::Mults> mult;
  Produces<aod::MultZeqs> multzeq;

  // For vertex-Z corrections in calibration
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using Run2Tracks = soa::Join<aod::Tracks, aod::TracksExtra>;
  Partition<Run2Tracks> run2tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));
  Partition<Run2Tracks> tracksWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run2Tracks> pvContribTracks = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run2Tracks> pvContribTracksEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  // Configurable
  Configurable<int> doVertexZeq{"doVertexZeq", 1, "if 1: do vertex Z eq mult table"};

  int mRunNumber;
  bool lCalibLoaded;
  TList* lCalibObjects;
  TProfile* hVtxZFV0A;
  TProfile* hVtxZFT0A;
  TProfile* hVtxZFT0C;
  TProfile* hVtxZFDDA;
  TProfile* hVtxZFDDC;
  TProfile* hVtxZNTracks;

  void init(InitContext& context)
  {
    if (doprocessRun2 == false && doprocessRun3 == false) {
      LOGF(fatal, "Neither processRun2 nor processRun3 enabled. Please choose one.");
    }
    if (doprocessRun2 == true && doprocessRun3 == true) {
      LOGF(fatal, "Cannot enable processRun2 and processRun3 at the same time. Please choose one.");
    }

    mRunNumber = 0;
    lCalibLoaded = false;
    lCalibObjects = nullptr;
    hVtxZFV0A = nullptr;
    hVtxZFT0A = nullptr;
    hVtxZFT0C = nullptr;
    hVtxZFDDA = nullptr;
    hVtxZFDDC = nullptr;
    hVtxZNTracks = nullptr;

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false); // don't fatal, please - exception is caught explicitly (as it should)
  }

  void processRun2(aod::Run2MatchedSparse::iterator const& collision,
                   Run2Tracks const& tracksExtra,
                   aod::BCs const&,
                   aod::Zdcs const&,
                   aod::FV0As const& fv0as,
                   aod::FV0Cs const& fv0cs,
                   aod::FT0s const& ft0s)
  {
    float multFV0A = 0.f;
    float multFV0C = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = 0.f;
    float multZNC = 0.f;

    auto trackletsGrouped = run2tracklets->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksGrouped = tracksWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int multTracklets = trackletsGrouped.size();
    int multTPC = tracksGrouped.size();
    int multNContribs = 0;
    int multNContribsEta1 = 0;

    if (collision.has_fv0a()) {
      for (auto amplitude : collision.fv0a().amplitude()) {
        multFV0A += amplitude;
      }
    }
    if (collision.has_fv0c()) {
      for (auto amplitude : collision.fv0c().amplitude()) {
        multFV0C += amplitude;
      }
    }
    if (collision.has_ft0()) {
      auto ft0 = collision.ft0();
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    }
    if (collision.has_zdc()) {
      auto zdc = collision.zdc();
      multZNA = zdc.energyCommonZNA();
      multZNC = zdc.energyCommonZNC();
    }

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC);
    mult(multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC, multNContribs, multNContribsEta1);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun2, "Produce Run 2 multiplicity tables", true);

  using Run3Tracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  Partition<Run3Tracks> tracksIUWithTPC = (aod::track::tpcNClsFindable > (uint8_t)0);
  Partition<Run3Tracks> pvContribTracksIU = (nabs(aod::track::eta) < 0.8f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  Partition<Run3Tracks> pvContribTracksIUEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);
  void processRun3(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                   Run3Tracks const& tracksExtra,
                   soa::Join<aod::BCs, aod::Timestamps> const& bcs,
                   aod::Zdcs const& zdcs,
                   aod::FV0As const& fv0as,
                   aod::FT0s const& ft0s,
                   aod::FDDs const& fdds)
  {
    float multFV0A = 0.f;
    float multFV0C = 0.f;
    float multFT0A = 0.f;
    float multFT0C = 0.f;
    float multFDDA = 0.f;
    float multFDDC = 0.f;
    float multZNA = 0.f;
    float multZNC = 0.f;
    int multTracklets = 0;

    float multZeqFV0A = 0.f;
    float multZeqFT0A = 0.f;
    float multZeqFT0C = 0.f;
    float multZeqFDDA = 0.f;
    float multZeqFDDC = 0.f;
    float multZeqNContribs = 0.f;

    auto tracksGrouped = tracksIUWithTPC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto pvContribsGrouped = pvContribTracksIU->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto pvContribsEta1Grouped = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    int multTPC = tracksGrouped.size();
    int multNContribs = pvContribsGrouped.size();
    int multNContribsEta1 = pvContribsEta1Grouped.size();

    /* check the previous run number */
    auto bc = collision.bc_as<soa::Join<aod::BCs, aod::Timestamps>>();
    if (doVertexZeq > 0) {
      if (bc.runNumber() != mRunNumber) {
        mRunNumber = bc.runNumber(); // mark this run as at least tried
        lCalibObjects = ccdb->getForTimeStamp<TList>("Centrality/Calibration", bc.timestamp());
        if (lCalibObjects) {
          hVtxZFV0A = (TProfile*)lCalibObjects->FindObject("hVtxZFV0A");
          hVtxZFT0A = (TProfile*)lCalibObjects->FindObject("hVtxZFT0A");
          hVtxZFT0C = (TProfile*)lCalibObjects->FindObject("hVtxZFT0C");
          hVtxZFDDA = (TProfile*)lCalibObjects->FindObject("hVtxZFDDA");
          hVtxZFDDC = (TProfile*)lCalibObjects->FindObject("hVtxZFDDC");
          hVtxZNTracks = (TProfile*)lCalibObjects->FindObject("hVtxZNTracksPV");
          lCalibLoaded = true;
          // Capture error
          if (!hVtxZFV0A || !hVtxZFT0A || !hVtxZFT0C || !hVtxZFDDA || !hVtxZFDDC || !hVtxZNTracks) {
            LOGF(error, "Problem loading CCDB objects! Please check");
            lCalibLoaded = false;
          }
        } else {
          LOGF(error, "Problem loading CCDB object! Please check");
          lCalibLoaded = false;
        }
      }
    }
    // using FT0 row index from event selection task
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (auto amplitude : ft0.amplitudeA()) {
        multFT0A += amplitude;
      }
      for (auto amplitude : ft0.amplitudeC()) {
        multFT0C += amplitude;
      }
    }
    // using FDD row index from event selection task
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      for (auto amplitude : fdd.chargeA()) {
        multFDDA += amplitude;
      }
      for (auto amplitude : fdd.chargeC()) {
        multFDDC += amplitude;
      }
    }
    // using FV0 row index from event selection task
    if (collision.has_foundFV0()) {
      auto fv0 = collision.foundFV0();
      for (auto amplitude : fv0.amplitude()) {
        multFV0A += amplitude;
      }
    }
    if (fabs(collision.posZ()) < 15.0f && lCalibLoaded) {
      multZeqFV0A = hVtxZFV0A->Interpolate(0.0) * multFV0A / hVtxZFV0A->Interpolate(collision.posZ());
      multZeqFT0A = hVtxZFT0A->Interpolate(0.0) * multFT0A / hVtxZFT0A->Interpolate(collision.posZ());
      multZeqFT0C = hVtxZFT0C->Interpolate(0.0) * multFT0C / hVtxZFT0C->Interpolate(collision.posZ());
      multZeqFDDA = hVtxZFDDA->Interpolate(0.0) * multFDDA / hVtxZFDDA->Interpolate(collision.posZ());
      multZeqFDDC = hVtxZFDDC->Interpolate(0.0) * multFDDC / hVtxZFDDC->Interpolate(collision.posZ());
      multZeqNContribs = hVtxZNTracks->Interpolate(0.0) * multNContribs / hVtxZNTracks->Interpolate(collision.posZ());
    }

    LOGF(debug, "multFV0A=%5.0f multFV0C=%5.0f multFT0A=%5.0f multFT0C=%5.0f multFDDA=%5.0f multFDDC=%5.0f multZNA=%6.0f multZNC=%6.0f multTracklets=%i multTPC=%i", multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC);
    mult(multFV0A, multFV0C, multFT0A, multFT0C, multFDDA, multFDDC, multZNA, multZNC, multTracklets, multTPC, multNContribs, multNContribsEta1);
    multzeq(multZeqFV0A, multZeqFT0A, multZeqFT0C, multZeqFDDA, multZeqFDDC, multZeqNContribs);
  }
  PROCESS_SWITCH(MultiplicityTableTaskIndexed, processRun3, "Produce Run 3 multiplicity tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityTableTaskIndexed>(cfgc, TaskName{"multiplicity-table"})};
}
