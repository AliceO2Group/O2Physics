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

/// \file testPV.cxx
/// \brief a task to check PV quality in MC
/// \author daiki.sekihata@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <string>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct testPV {

  // Configurables
  // Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  struct : ConfigurableGroup {
    std::string prefix = "eventCut";
    Configurable<int> cfgEventGeneratorId{"cfgEventGeneratorId", -1, "event generator index. e.g. select gap/signal events"};
    Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for AA"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventCut;

  HistogramRegistry fRegistry{"fRegistry"};

  void init(o2::framework::InitContext&)
  {
    // ccdb->setURL(ccdburl);
    // ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // ccdb->setFatalWhenNull(false);

    mRunNumber = 0;

    addHistograms();
  }

  int mRunNumber{0};
  // Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    mRunNumber = bc.runNumber();
  }

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("hCollisionCounter", "collision counter", kTH1D, {{2, 0.5f, 2.5f}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    fRegistry.add("Vertex/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Vertex/hNContrib", "Number of PV contributors;N_{contrib}", kTH1F, {{101, -0.5, 100.5}}, false);
    fRegistry.add("Vertex/hChi2", "vertex chi2;#chi^{2}/N_{contrib}", kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Vertex/hChi2vsNContrib", "vertex #chi^{2}/N_{contrib} vs. N_{contrib};N_{contrib};#chi^{2}/N_{contrib}", kTH2F, {{101, -0.5, 100.5}, {100, 0, 10}}, false);

    fRegistry.add("Vertex/hSigmaX", "vertex #sigma_{X} vs. N_{contrib};N_{contrib};#sigma_{X} (#mum)", kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/hSigmaY", "vertex #sigma_{Y} vs. N_{contrib};N_{contrib};#sigma_{Y} (#mum)", kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/hSigmaZ", "vertex #sigma_{Z} vs. N_{contrib};N_{contrib};#sigma_{Z} (#mum)", kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);

    fRegistry.add("Vertex/hDeltaX", "vertex #DeltaX vs. N_{contrib};N_{contrib};#DeltaX = (X_{rec} #minus X_{gen})/#sigma_{X}", kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);
    fRegistry.add("Vertex/hDeltaY", "vertex #DeltaY vs. N_{contrib};N_{contrib};#DeltaY = (Y_{rec} #minus Y_{gen})/#sigma_{Y}", kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);
    fRegistry.add("Vertex/hDeltaZ", "vertex #DeltaZ vs. N_{contrib};N_{contrib};#DeltaZ = (Z_{rec} #minus Z_{gen})/#sigma_{Z}", kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);

    fRegistry.add("Vertex/hCollisionTime", "vertex time;N_{contrib};collision time (ns)", kTH2F, {{101, -0.5, 100.5}, {500, -25, 25}}, false);
    fRegistry.add("Vertex/hCollisionTimeRes", "vertex time resolution;N_{contrib};collision time resolution (ns)", kTH2F, {{101, -0.5, 100.5}, {250, 0, 25}}, false);
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (!(eventCut.cfgZvtxMin < collision.posZ() && collision.posZ() < eventCut.cfgZvtxMax)) {
      return false;
    }
    if (eventCut.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventCut.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (eventCut.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (eventCut.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventCut.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (eventCut.cfgRequireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (eventCut.cfgRequireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    return true;
  }

  template <typename TCollision, typename TMCCollision>
  void fillVertexHistograms(TCollision const& collision, TMCCollision const& mcCollision)
  {
    fRegistry.fill(HIST("Vertex/hZvtx"), collision.posZ());
    fRegistry.fill(HIST("Vertex/hNContrib"), collision.numContrib());
    fRegistry.fill(HIST("Vertex/hChi2"), collision.chi2() / collision.numContrib());
    fRegistry.fill(HIST("Vertex/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());

    fRegistry.fill(HIST("Vertex/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Vertex/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Vertex/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

    fRegistry.fill(HIST("Vertex/hDeltaX"), collision.numContrib(), (collision.posX() - mcCollision.posX()) / std::sqrt(collision.covXX()));
    fRegistry.fill(HIST("Vertex/hDeltaY"), collision.numContrib(), (collision.posY() - mcCollision.posY()) / std::sqrt(collision.covYY()));
    fRegistry.fill(HIST("Vertex/hDeltaZ"), collision.numContrib(), (collision.posZ() - mcCollision.posZ()) / std::sqrt(collision.covZZ()));

    fRegistry.fill(HIST("Vertex/hCollisionTime"), collision.numContrib(), collision.collisionTime());
    fRegistry.fill(HIST("Vertex/hCollisionTimeRes"), collision.numContrib(), collision.collisionTimeRes());
  }

  template <typename TBCs, typename TCollisions, typename TMCCollisions, typename TMCParticles>
  void run(TBCs const&, TCollisions const& collisions, TMCCollisions const&, TMCParticles const&)
  {
    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto bc = collision.template bc_as<TBCs>();
      initCCDB(bc);

      auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (eventCut.cfgEventGeneratorId > -1 && mcCollision.getSubGeneratorId() != eventCut.cfgEventGeneratorId) {
        continue;
      }

      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (!isSelectedCollision(collision)) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 2);

      fillVertexHistograms(collision, mcCollision);

    } // end of collision loop
  }

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using MyBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;

  Filter collisionFilter_evsel = eventCut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventCut.cfgZvtxMax;
  using FilteredMyCollisions = soa::Filtered<MyCollisions>;

  void processMC(FilteredMyCollisions const& collisions, MyBCs const& bcs, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    run(bcs, collisions, mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(testPV, processMC, "processMC", true);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(testPV, processDummy, "processDummy", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<testPV>(cfgc, TaskName{"test-pv"})};
}
