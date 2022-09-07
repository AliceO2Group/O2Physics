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

// Inspired by phianalysis.cxx

/// \file lambda1520analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
/// \author Hirak Kumar Koley <hirak.kumar.koley@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include <CCDB/BasicCCDBManager.h>
#include "DataFormatsParameters/GRPObject.h"
#include <math.h>
#include "TMath.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace std;
using std::array;

enum EventSelection { kNoSelection = 0,
                      kVertexCut = 1 };

struct lambda1520analysis {
  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 100.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  // Configurables

  /// Event cuts
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  // Configurable<float>cfgCutPtmin{"cfgCutPtmin", 0.1f, "Minimal pT for tracks"};
  // Configurable<float>cfgCutPtmax{"cfgCutPtmax", 10.0f, "Minimal pT for tracks"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 5.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  /// Partition for first Track selecting as Kaon
  Partition<aod::ResoDaughters> parts1 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDKaonCutInFilter() && requireTOFPIDKaonCutInFilter() && (nabs(aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(aod::track::dcaXY) < static_cast<float_t>((cMaxDCArToPVcut))); // Basic DCA cuts
  // && (nabs(aod::track::dcaXY) <= (0.0105f + 0.0350f / npow(aod::track::pt, 1.1f))) // strict DCAz cut
  // && (nabs(aod::track::eta) < static_cast<float_t>(cfgCutEta))
  // && (aod::track::pt < static_cast<float_t>(cfgCutPtmax))
  // && (aod::track::pt > static_cast<float_t>(cfgCutPtmin));
  /// Partition for second Track selecting as proton
  Partition<aod::ResoDaughters> parts2 = (aod::resodaughter::partType == uint8_t(aod::resodaughter::DaughterType::kTrack)) && requireTPCPIDProtonCutInFilter() && requireTOFPIDProtonCutInFilter() && (nabs(aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(aod::track::dcaXY) < static_cast<float_t>((cMaxDCArToPVcut))); // Basic DCA cuts
  // && (nabs(aod::track::dcaXY) <= (0.0105f + 0.0350f / npow(aod::track::pt, 1.1f))) // strict DCAz cut
  // && (nabs(aod::track::eta) < static_cast<float_t>(cfgCutEta))
  // && (aod::track::pt < static_cast<float_t>(cfgCutPtmax))
  // && (aod::track::pt > static_cast<float_t>(cfgCutPtmin));

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

    //  Event QA
    auto h = histos.add<TH1>("QA/Event/EnumEvents",
                             "Event selection;"
                             ";"
                             "Selected Events",
                             kTH1F, {{10, 0, 10}});
    h->GetXaxis()->SetBinLabel(1, "Events read and Ev. sel passed");
    h->GetXaxis()->SetBinLabel(2, "posZ passed");
    histos.add("QA/Event/VertexZ",
               "Event selection;"
               "Vertex Z (cm);"
               "Selected Events",
               kTH1F, {{4000, -20, 20}});
    histos.add("QA/Event/Selected/VertexZ",
               "Event selection;"
               "Vertex Z (cm);"
               "Selected Events",
               kTH1F, {{4000, -20, 20}});

    //  PID QA
    //  --- Kaon
    histos.add("QA/Kaon/TOF_TPC_Map",
               "TOF + TPC Combined PID for Kaons;"
               "#sigma_{TOF}^{Kaon};"
               "#sigma_{TPC}^{Kaon}",
               kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    histos.add("QA/Kaon/TOF_Nsigma",
               "TOF NSigma for Kaons;"
               "#it{p}_{T} (GeV/#it{c});"
               "#sigma_{TOF}^{Kaon};",
               kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    histos.add("QA/Kaon/TPC_Nsigma",
               "TPC NSigma for Kaons;"
               "#it{p}_{T} (GeV/#it{c});"
               "#sigma_{TPC}^{Kaon};",
               kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    //  --- Proton
    histos.add("QA/Proton/TOF_TPC_Map",
               "TOF + TPC Combined PID for Protons;"
               "#sigma_{TOF}^{Proton};"
               "#sigma_{TPC}^{Proton}",
               kTH2F, {{1000, -10, 10}, {1000, -10, 10}});
    histos.add("QA/Proton/TOF_Nsigma",
               "TOF NSigma for Protons;"
               "#it{p}_{T} (GeV/#it{c});"
               "#sigma_{TOF}^{Proton};",
               kTH2F, {{1000, 0, 20}, {1000, -10, 10}});
    histos.add("QA/Proton/TPC_Nsigma",
               "TPC NSigma for Protons;"
               "#it{p}_{T} (GeV/#it{c});"
               "#sigma_{TPC}^{Proton};",
               kTH2F, {{1000, 0, 20}, {1000, -10, 10}});

    // Mass QA (quick check)
    histos.add("Analysis/lambda1520invmass", "Invariant mass of Lambda1520 K^{#pm}p^{#mp}", kTH1F, {{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("Analysis/lambda1520invmassLS", "Invariant mass of Lambda1520 Like Sign Method K^{#pm}p^{#pm}", kTH1F, {{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("Analysis/lambda1520invmassLSkp", "Invariant mass of Lambda1520 Like Sign Method K^{#plus}p^{#plus}", kTH1F, {{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});         // K+ + Pr
    histos.add("Analysis/lambda1520invmassLSkbarpbar", "Invariant mass of Lambda1520 Like Sign Method K^{#minus}p^{#minus}", kTH1F, {{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}}); // K- + anti-Pr
    histos.add("Analysis/lambda1520invmassME", "Invariant mass of Lambda1520 mixed event K^{#pm}p^{#mp}", kTH1F, {{500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});

    // 3d histogram
    histos.add("Analysis/h3lambda1520invmass", "Invariant mass of Lambda1520 K^{#pm}p^{#mp}", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 1.3, 1.8}});
    histos.add("Analysis/h3lambda1520invmassLS", "Invariant mass of Lambda1520 Like Sign Method K^{#pm}p^{#pm}", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 1.3, 1.8}});
    histos.add("Analysis/h3lambda1520invmassLSkp", "Invariant mass of Lambda1520 Like Sign Method K^{#plus}p^{#plus}", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 1.3, 1.8}});         // K+ + Pr
    histos.add("Analysis/h3lambda1520invmassLSkbarpbar", "Invariant mass of Lambda1520 Like Sign Method K^{#minus}p^{#minus}", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 1.3, 1.8}}); // K- + anti-Pr
    histos.add("Analysis/h3lambda1520invmassME", "Invariant mass of Lambda1520 mixed event K^{#pm}p^{#mp}", kTH3F, {{100, 0.0f, 100.0f}, {100, 0.0f, 10.0f}, {500, 1.3, 1.8}});
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();

  //  template <typename aod::ResoCollision>
  void process(aod::ResoCollision& collision,
               aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    //  Collision QA
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);

    //  --- Event selection: Vertex position
    histos.fill(HIST("QA/Event/VertexZ"), collision.posZ());
    if (!(fabs(collision.posZ()) < ConfEvtZvtx))
      return;
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    histos.fill(HIST("QA/Event/Selected/VertexZ"), collision.posZ());

    // LOGF(info, "event id: %d", collision.bcId());
    auto group1 = parts1->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
    auto group2 = parts2->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());

    for (auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(group1, group2))) {

      //  PID Selection
      //  --- PID QA Kaons
      histos.fill(HIST("QA/Kaon/TOF_Nsigma"), trk1.pt(), trk1.tofNSigmaKa());
      histos.fill(HIST("QA/Kaon/TPC_Nsigma"), trk1.pt(), trk1.tpcNSigmaKa());
      histos.fill(HIST("QA/Kaon/TOF_TPC_Map"), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());
      //  --- PID QA Pions
      histos.fill(HIST("QA/Proton/TOF_Nsigma"), trk2.pt(), trk2.tofNSigmaPi());
      histos.fill(HIST("QA/Proton/TPC_Nsigma"), trk2.pt(), trk2.tpcNSigmaPi());
      histos.fill(HIST("QA/Proton/TOF_TPC_Map"), trk2.tofNSigmaPi(), trk2.tpcNSigmaPi());

      auto arrMom = array{
        array{trk1.px(), trk1.py(), trk1.pz()},
        array{trk2.px(), trk2.py(), trk2.pz()}};
      auto resoPt = RecoDecay::sqrtSumOfSquares(trk1.px() + trk2.px(), trk1.py() + trk2.py());
      auto arrMass = array{massKa, massPr};
      auto mass = RecoDecay::m(arrMom, arrMass);
      // auto arrMom3 = array{trk1.px() + trk2.px(), trk1.py() + trk2.py(), trk1.pz() + trk2.pz()};
      // auto resorap = RecoDecay::y(arrMom3, mass);
      // if (abs(resorap) > 0.5) continue;

      // like sign pair only
      if (trk1.sign() * trk2.sign() > 0) {
        histos.fill(HIST("Analysis/lambda1520invmassLS"), mass);
        histos.fill(HIST("Analysis/h3lambda1520invmassLS"), collision.multV0M(), resoPt, mass);

        if (trk1.sign() > 0 && trk2.sign() > 0) { // K+ + Pr
          histos.fill(HIST("Analysis/lambda1520invmassLSkp"), mass);
          histos.fill(HIST("Analysis/h3lambda1520invmassLSkp"), collision.multV0M(), resoPt, mass);
        } else { // K- + anti-Pr
          histos.fill(HIST("Analysis/lambda1520invmassLSkbarpbar"), mass);
          histos.fill(HIST("Analysis/h3lambda1520invmassLSkbarpbar"), collision.multV0M(), resoPt, mass);
        }
      } else { // Un-like sign pair only
        histos.fill(HIST("Analysis/lambda1520invmass"), mass);
        histos.fill(HIST("Analysis/h3lambda1520invmass"), collision.multV0M(), resoPt, mass);
      }
    }
  }

  // Processing Event Mixing
  void processME(o2::aod::ResoCollisions& collision,
                 o2::aod::BCsWithTimestamps const&, aod::ResoDaughters const&, aod::Reso2TracksPIDExt const&)
  {
    ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M> colBinning{{CfgVtxBins, CfgMultBins}, true};

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
        auto arrMass = array{massKa, massPr};
        auto mass = RecoDecay::m(arrMom, arrMass);
        // auto arrMom3 = array{trk1.px() + trk2.px(), trk1.py() + trk2.py(), trk1.pz() + trk2.pz()};
        // auto resorap = RecoDecay::y(arrMom3, mass);
        // if (abs(resorap) > 0.5) continue;
        histos.fill(HIST("Analysis/lambda1520invmassME"), mass);
        histos.fill(HIST("Analysis/h3lambda1520invmassME"), collision1.multV0M(), resoPt, mass);
      }
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processME, "Process EventMixing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambda1520analysis>(cfgc, TaskName{"lf-lambda1520analysis"})};
}
