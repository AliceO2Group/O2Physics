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
#include <TLorentzVector.h>

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
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  /// Partition for firstTrack
  Partition<aod::ResoDaughters> parts1 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter() && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
  Partition<aod::ResoDaughters> parts2 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter() && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
  // Partition<aod::ResoDaughters> parts1 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter(); // w/o Basic DCA cuts
  // Partition<aod::ResoDaughters> parts2 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter(); // w/o Basic DCA cuts
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
    histos.add("phiinvmass", "Invariant mass of Phi", kTH1F, {{700, 0.8, 1.5, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("phiinvmassME", "Invariant mass of Phi mixed event", kTH1F, {{700, 0.8, 1.5, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("trk1pT", "pT distribution of track1", kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("trk2pT", "pT distribution of track1", kTH1F, {{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("TOF_TPC_Map1", "TOF + TPC Combined PID for Kaons;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_Nsigma1", "TOF NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TPC_Nsigma1", "TPC NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_TPC_Map2", "TOF + TPC Combined PID for Kaons;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_Nsigma2", "TOF NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TPC_Nsigma2", "TPC NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});

    // 3d histogram
    histos.add("h3phiinvmass", "Invariant mass of Phi", kTH3F, {{300, 0, 3000}, {100, 0.0f, 10.0f}, {700, 0.8, 1.5}});
    histos.add("h3phiinvmassME", "Invariant mass of Phi mixed event", kTH3F, {{300, 0, 3000}, {100, 0.0f, 10.0f}, {700, 0.8, 1.5}});
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();

  void process(aod::ResoCollision& collision,
               aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    // LOGF(info, "event id: %d", collision.bcId());
    auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
    auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {
      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() > 0)
        continue;
      if (trk1.sign() > 0) {
        //  --- PID QA Kaons +
        histos.fill(HIST("TOF_Nsigma1"), trk1.pt(), trk1.tofNSigmaKa());
        histos.fill(HIST("TPC_Nsigma1"), trk1.pt(), trk1.tpcNSigmaKa());
        histos.fill(HIST("TOF_TPC_Map1"), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());
        //  --- PID QA Kaons -
        histos.fill(HIST("TOF_Nsigma2"), trk2.pt(), trk2.tofNSigmaKa());
        histos.fill(HIST("TPC_Nsigma2"), trk2.pt(), trk2.tpcNSigmaKa());
        histos.fill(HIST("TOF_TPC_Map2"), trk2.tofNSigmaKa(), trk2.tpcNSigmaKa());
      } else {
        //  --- PID QA Kaons +
        histos.fill(HIST("TOF_Nsigma1"), trk2.pt(), trk2.tofNSigmaKa());
        histos.fill(HIST("TPC_Nsigma1"), trk2.pt(), trk2.tpcNSigmaKa());
        histos.fill(HIST("TOF_TPC_Map1"), trk2.tofNSigmaKa(), trk2.tpcNSigmaKa());
        //  --- PID QA Kaons -
        histos.fill(HIST("TOF_Nsigma2"), trk1.pt(), trk1.tofNSigmaKa());
        histos.fill(HIST("TPC_Nsigma2"), trk1.pt(), trk1.tpcNSigmaKa());
        histos.fill(HIST("TOF_TPC_Map2"), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());
      }

      histos.fill(HIST("trk1pT"), trk1.pt());
      histos.fill(HIST("trk2pT"), trk2.pt());

      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;

      histos.fill(HIST("phiinvmass"), lResonance.M());
      histos.fill(HIST("h3phiinvmass"), collision.multV0M(), lResonance.Pt(), lResonance.M());
    }
  }

  // Processing Event Mixing
  void processME(o2::aod::ResoCollisions& collision,
                 o2::aod::BCsWithTimestamps const&, aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultTPCtemp> colBinning{{CfgVtxBins, CfgMultBins}, true};

    for (auto& [collision1, collision2] : soa::selfCombinations(colBinning, 10, -1, collision, collision)) {

      auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision1.globalIndex());
      auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision2.globalIndex());

      TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
      for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {
        // Un-like sign pair only
        if (trk1.sign() * trk2.sign() > 0)
          continue;

        lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
        lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
        lResonance = lDecayDaughter1 + lDecayDaughter2;

        if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
          continue;

        histos.fill(HIST("phiinvmassME"), lResonance.M());
        histos.fill(HIST("h3phiinvmassME"), collision1.multV0M(), lResonance.Pt(), lResonance.M());
      }
    }
  };
  PROCESS_SWITCH(phianalysis, processME, "Process EventMixing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phianalysis>(cfgc, TaskName{"lf-phianalysis"})};
}
