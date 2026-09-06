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
#include "Common/Core/Zorro.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <string>

struct testPV {

  // Configurables
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventCut";
    o2::framework::Configurable<int> cfgEventGeneratorId{"cfgEventGeneratorId", -1, "event generator index. e.g. select gap/signal events"};
    o2::framework::Configurable<float> cfgChi2PerNcontribMax{"cfgChi2PerNcontribMax", 1e+10, "max. chi2/Ncontrib of PV"};
    o2::framework::Configurable<float> cfgZvtxMin{"cfgZvtxMin", -10.f, "min. Zvtx"};
    o2::framework::Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    o2::framework::Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    o2::framework::Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    o2::framework::Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    o2::framework::Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.

    // for RCT
    o2::framework::Configurable<bool> cfgRequireGoodRCT{"cfgRequireGoodRCT", true, "require good detector flag in run condtion table"};
    o2::framework::Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "CBT", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
    o2::framework::Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "set ZDC flag for AA"};
    o2::framework::Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } eventCut;

  // for zorro
  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "zorroGroup";
    o2::framework::Configurable<std::string> cfgTriggerName{"cfgTriggerName", "fGlobalDimuon", "desired software trigger name"};
    o2::framework::Configurable<std::string> ccdbPathSoftwareTrigger{"ccdbPathSoftwareTrigger", "EventFiltering/Zorro/", "ccdb path for ZORRO objects"};
    o2::framework::Configurable<uint64_t> bcMarginForSoftwareTrigger{"bcMarginForSoftwareTrigger", 100, "Number of BCs of margin for software triggers"};
  } zorroGroup;

  o2::framework::HistogramRegistry fRegistry{"fRegistry"};
  Zorro zorro;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mRunNumber = 0;

    addHistograms();
  }

  int mRunNumber{0};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  template <bool isTriggerAnalysis, typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    if constexpr (isTriggerAnalysis) {
      zorro.setCCDBpath(zorroGroup.ccdbPathSoftwareTrigger);
      zorro.setBCtolerance(zorroGroup.bcMarginForSoftwareTrigger); // this does nothing.
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), zorroGroup.cfgTriggerName.value);
    }

    mRunNumber = bc.runNumber();
  }

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("hCollisionCounter", "collision counter", o2::framework::HistType::kTH1D, {{2, 0.5f, 2.5f}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    fRegistry.add("Vertex/hZvtx", "vertex z; Z_{vtx} (cm)", o2::framework::HistType::kTH1F, {{100, -50, +50}}, false);
    fRegistry.add("Vertex/hNContrib", "Number of PV contributors;N_{contrib}", o2::framework::HistType::kTH1F, {{101, -0.5, 100.5}}, false);
    fRegistry.add("Vertex/hChi2", "vertex chi2;#chi^{2}/N_{contrib}", o2::framework::HistType::kTH1F, {{100, 0, 10}}, false);
    fRegistry.add("Vertex/hChi2vsNContrib", "vertex #chi^{2}/N_{contrib} vs. N_{contrib};N_{contrib};#chi^{2}/N_{contrib}", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {100, 0, 10}}, false);

    fRegistry.add("Vertex/hSigmaX", "vertex #sigma_{X} vs. N_{contrib};N_{contrib};#sigma_{X} (#mum)", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/hSigmaY", "vertex #sigma_{Y} vs. N_{contrib};N_{contrib};#sigma_{Y} (#mum)", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/hSigmaZ", "vertex #sigma_{Z} vs. N_{contrib};N_{contrib};#sigma_{Z} (#mum)", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {1000, 0, 100}}, false);

    fRegistry.add("Vertex/hCollisionTime", "vertex time;N_{contrib};collision time (ns)", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {500, -25, 25}}, false);
    fRegistry.add("Vertex/hCollisionTimeRes", "vertex time resolution;N_{contrib};collision time resolution (ns)", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {250, 0, 25}}, false);

    if (doprocessMC) {
      fRegistry.add("Vertex/hDeltaX", "vertex #DeltaX vs. N_{contrib};N_{contrib};#DeltaX = (X_{rec} #minus X_{gen})/#sigma_{X}", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);
      fRegistry.add("Vertex/hDeltaY", "vertex #DeltaY vs. N_{contrib};N_{contrib};#DeltaY = (Y_{rec} #minus Y_{gen})/#sigma_{Y}", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);
      fRegistry.add("Vertex/hDeltaZ", "vertex #DeltaZ vs. N_{contrib};N_{contrib};#DeltaZ = (Z_{rec} #minus Z_{gen})/#sigma_{Z}", o2::framework::HistType::kTH2F, {{101, -0.5, 100.5}, {200, -10, 10}}, false);
    }
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
    if (collision.chi2() / collision.numContrib() > eventCut.cfgChi2PerNcontribMax) {
      return false;
    }
    return true;
  }

  template <typename TCollision>
  void fillVertexHistograms(TCollision const& collision)
  {
    fRegistry.fill(HIST("Vertex/hZvtx"), collision.posZ());
    fRegistry.fill(HIST("Vertex/hNContrib"), collision.numContrib());
    fRegistry.fill(HIST("Vertex/hChi2"), collision.chi2() / collision.numContrib());
    fRegistry.fill(HIST("Vertex/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());

    fRegistry.fill(HIST("Vertex/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Vertex/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
    fRegistry.fill(HIST("Vertex/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

    fRegistry.fill(HIST("Vertex/hCollisionTime"), collision.numContrib(), collision.collisionTime());
    fRegistry.fill(HIST("Vertex/hCollisionTimeRes"), collision.numContrib(), collision.collisionTimeRes());
  }

  template <bool isMC, bool isTriggerAnalysis, typename TBCs, typename TCollisions, typename TMCCollisions, typename TMCParticles>
  void run(TBCs const&, TCollisions const& collisions, TMCCollisions const&, TMCParticles const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<TBCs>();
      initCCDB<isTriggerAnalysis>(bc);

      if constexpr (isMC) {
        if (!collision.has_mcCollision()) {
          continue;
        }
        auto mcCollision = collision.template mcCollision_as<TMCCollisions>();
        if (eventCut.cfgEventGeneratorId > -1 && mcCollision.getSubGeneratorId() != eventCut.cfgEventGeneratorId) {
          continue;
        }
      }

      if constexpr (isTriggerAnalysis) {
        if (!zorro.isSelected(bc.globalBC(), zorroGroup.bcMarginForSoftwareTrigger)) { // triggered event
          continue;
        }
      }

      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (!isSelectedCollision(collision)) {
        continue;
      }

      fRegistry.fill(HIST("hCollisionCounter"), 2);
      fillVertexHistograms(collision);
      if constexpr (isMC) {
        auto mcCollision = collision.template mcCollision_as<TMCCollisions>();
        fRegistry.fill(HIST("Vertex/hDeltaX"), collision.numContrib(), (collision.posX() - mcCollision.posX()) / std::sqrt(collision.covXX()));
        fRegistry.fill(HIST("Vertex/hDeltaY"), collision.numContrib(), (collision.posY() - mcCollision.posY()) / std::sqrt(collision.covYY()));
        fRegistry.fill(HIST("Vertex/hDeltaZ"), collision.numContrib(), (collision.posZ() - mcCollision.posZ()) / std::sqrt(collision.covZZ()));
      }

    } // end of collision loop
  }

  using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using MyCollisionsMC = o2::soa::Join<MyCollisions, o2::aod::McCollisionLabels>;
  using MyBCs = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::BcSels>;

  o2::framework::expressions::Filter collisionFilter_evsel = eventCut.cfgZvtxMin < o2::aod::collision::posZ && o2::aod::collision::posZ < eventCut.cfgZvtxMax;
  using FilteredMyCollisions = o2::soa::Filtered<MyCollisions>;
  using FilteredMyCollisionsMC = o2::soa::Filtered<MyCollisionsMC>;

  void processData(FilteredMyCollisions const& collisions, MyBCs const& bcs)
  {
    run<false, false>(bcs, collisions, nullptr, nullptr);
  }
  PROCESS_SWITCH(testPV, processData, "processData", true);

  void processTriggeredData(FilteredMyCollisions const& collisions, MyBCs const& bcs)
  {
    run<false, true>(bcs, collisions, nullptr, nullptr);
  }
  PROCESS_SWITCH(testPV, processTriggeredData, "processTriggeredData", false);

  void processMC(FilteredMyCollisionsMC const& collisions, MyBCs const& bcs, o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles)
  {
    run<true, false>(bcs, collisions, mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(testPV, processMC, "processMC", false);

  void processDummy(MyCollisions const&) {}
  PROCESS_SWITCH(testPV, processDummy, "processDummy", false);
};
o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{adaptAnalysisTask<testPV>(cfgc, o2::framework::TaskName{"test-pv"})};
}
