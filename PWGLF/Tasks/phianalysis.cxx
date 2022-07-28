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

/// \file phianalysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct phianalysis {
  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  // Configurables

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 5.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  /// Partition for firstTrack
  Partition<aod::ResoDaughters> parts1 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter() && (nabs(aod::resodaughter::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(aod::resodaughter::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(aod::resodaughter::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
  Partition<aod::ResoDaughters> parts2 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter() && (nabs(aod::resodaughter::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(aod::resodaughter::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(aod::resodaughter::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    long now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec vtxZAxis = {100, -20, 20};

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Mass QA (quick check)
    histos.add("phiinvmass", "Invariant mass of Phi", kTH1F, {{500, 0.8, 1.3, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("phiinvmassME", "Invariant mass of Phi mixed event", kTH1F, {{500, 0.8, 1.3, "Invariant Mass (GeV/#it{c}^2)"}});

    // 3d histogram
    histos.add("h3phiinvmass", "Invariant mass of Phi", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 0.8, 1.3}});
    histos.add("h3phiinvmassME", "Invariant mass of Phi mixed event", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 0.8, 1.3}});
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();

  void process(aod::ResoCollision& collision,
               aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    // LOGF(info, "event id: %d", collision.bcId());
    auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
    auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());

    for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {
      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() > 0)
        continue;

      auto arrMom = array{
        array{trk1.px(), trk1.py(), trk1.pz()},
        array{trk2.px(), trk2.py(), trk2.pz()}};
      auto resoPt = RecoDecay::sqrtSumOfSquares(trk1.px() + trk2.px(), trk1.py() + trk2.py());
      auto arrMass = array{massKa, massKa};
      auto mass = RecoDecay::m(arrMom, arrMass);
      histos.fill(HIST("phiinvmass"), mass);
      histos.fill(HIST("h3phiinvmass"), collision.multV0M(), resoPt, mass);
    }
  }

  // Processing Event Mixing
  void processME(o2::aod::ResoCollisions& collision,
                 o2::aod::BCsWithTimestamps const&, aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    BinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 5, -1, collision, collision)) {

      auto magFieldTesla1 = 0.0;
      static o2::parameters::GRPObject* grpo1 = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", collision1.timestamp());
      if (!grpo1) {
        magFieldTesla1 = 0;
      } else {
        // LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", collision1.timestamp(), grpo1->getNominalL3Field());
        magFieldTesla1 = 0.1 * (grpo1->getNominalL3Field());
      }

      auto magFieldTesla2 = 0.0;
      static o2::parameters::GRPObject* grpo2 = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", collision2.timestamp());
      if (!grpo2) {
        magFieldTesla2 = 0;
      } else {
        // LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", collision2.timestamp(), grpo2->getNominalL3Field());
        magFieldTesla2 = 0.1 * (grpo1->getNominalL3Field());
      }

      auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex());
      auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex());

      if (magFieldTesla1 != magFieldTesla2) {
        continue;
      }
      for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {
        // Un-like sign pair only
        if (trk1.sign() * trk2.sign() > 0)
          continue;

        auto arrMom = array{
          array{trk1.px(), trk1.py(), trk1.pz()},
          array{trk2.px(), trk2.py(), trk2.pz()}};
        auto resoPt = RecoDecay::sqrtSumOfSquares(trk1.px() + trk2.px(), trk1.py() + trk2.py());
        auto arrMass = array{massKa, massKa};
        auto mass = RecoDecay::m(arrMom, arrMass);
        histos.fill(HIST("phiinvmassME"), mass);
        histos.fill(HIST("h3phiinvmassME"), collision1.multV0M(), resoPt, mass);
      }
    }
  };
  PROCESS_SWITCH(phianalysis, processME, "Process EventMixing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phianalysis>(cfgc, TaskName{"lf-phianalysis"})};
}