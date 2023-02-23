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

#include <CCDB/BasicCCDBManager.h>
#include <TLorentzVector.h>
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "DataFormatsParameters/GRPObject.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace std;

enum EventSelection { kNoSelection = 0,
                      kVertexCut = 1 };

struct lambda1520analysis {
  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};

  // Configurables
  // Event cuts
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cMinPtcut{"cMinPtcut", 0.1f, "Minimal pT for tracks"};

  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 5.0f, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0f, "Track DCAz cut to PV Minimum"};
  // PID selection cuts
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 3.0f, "TPC and TOF PID cut (loose, improve performance)"};
  Configurable<float> pidnSigmaTPCCuthasTOF{"pidnSigmaTPCCuthasTOF", 5.0f, "nSigma cut for TPC PID if track has TOF Information"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 10, "Number of mixed events per event"};

  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    AxisSpec vtxZAxis = {100, -20, 20};

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 10., 20.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    //  Event QA
    auto h = histos.add<TH1>("QA/Event/EnumEvents", "Event selection; ; Event Count", {HistType::kTH1F, {{10, 0, 10}}});
    h->GetXaxis()->SetBinLabel(1, "Events read and Ev. sel. passed");
    h->GetXaxis()->SetBinLabel(2, "posZ passed");
    histos.add("QA/Event/VertexZ", "Event selection; Vertex Z (cm); Selected Events", {HistType::kTH1F, {{100, -20, 20}}});
    histos.add("QA/Event/Selected/VertexZ", "Event selection; Vertex Z (cm); Selected Events", {HistType::kTH1F, {{100, -20, 20}}});

    //  PID QA
    //  --- Kaon
    histos.add("QA/Kaon/TOF_TPC_Map", "TOF + TPC Combined PID for Kaons; #sigma_{TOF}^{Kaon}; #sigma_{TPC}^{Kaon}; Counts;", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("QA/Kaon/TOF_Nsigma", "TOF NSigma for Kaons; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Kaon};", {HistType::kTH2F, {{100, 0, 10}, {1000, -10, 10}}});
    histos.add("QA/Kaon/TPC_Nsigma", "TPC NSigma for Kaons; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Kaon};", {HistType::kTH2F, {{100, 0, 10}, {1000, -10, 10}}});
    //  --- Proton
    histos.add("QA/Proton/TOF_TPC_Map", "TOF + TPC Combined PID for Protons; #sigma_{TOF}^{Proton}; #sigma_{TPC}^{Proton}; Counts;", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("QA/Proton/TOF_Nsigma", "TOF NSigma for Protons; #it{p}_{T} (GeV/#it{c}); #sigma_{TOF}^{Proton};", {HistType::kTH2F, {{100, 0, 10}, {1000, -10, 10}}});
    histos.add("QA/Proton/TPC_Nsigma", "TPC NSigma for Protons; #it{p}_{T} (GeV/#it{c}); #sigma_{TPC}^{Proton};", {HistType::kTH2F, {{100, 0, 10}, {1000, -10, 10}}});

    histos.add("QA/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {{100, 0, 10}}});
    histos.add("QA/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {{100, 0, 10}}});

    // Mass QA (quick check)
    histos.add("Analysis/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}});
    histos.add("Analysis/lambda1520invmassLS", "Invariant mass of #Lambda(1520) Like Sign Method K^{#pm}p^{#pm}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}});
    histos.add("Analysis/lambda1520invmassLSkp", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}});         // K+ + Pr
    histos.add("Analysis/lambda1520invmassLSkbarpbar", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}}); // K- + anti-Pr
    histos.add("Analysis/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}});

    // 3d histogram
    histos.add("Analysis/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("Analysis/h3lambda1520invmassLS", "Invariant mass of #Lambda(1520) Like Sign Method K^{#pm}p^{#pm}", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("Analysis/h3lambda1520invmassLSkp", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});         // K+ + Pr
    histos.add("Analysis/h3lambda1520invmassLSkbarpbar", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}}); // K- + anti-Pr
    histos.add("Analysis/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});

    if (doprocessMC) {
      histos.add("MC/h3recolambda1520invmass", "Invariant mass of Reconstructed MC #Lambda(1520)", HistType::kTH3F, {{200, 0.0f, 2000.0f, "mult_{TPC}"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"}, {500, 1.3, 1.8, "Invariant Mass (GeV/#it{c}^2)"}});
      histos.add("MC/truelambda1520pt", "pT distribution of True MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {{100, 0, 10}}});
      histos.add("MC/reconlambda1520pt", "pT distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {{100, 0, 10}}});
      histos.add("MC/reconlambda1520invmass", "Inv mass distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {{500, 1.3, 1.8}}});
    }
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    //  Collision QA
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kNoSelection);

    //  --- Event selection: Vertex position
    histos.fill(HIST("QA/Event/VertexZ"), collision.posZ());
    if (!(fabs(collision.posZ()) < ConfEvtZvtx))
      return;
    histos.fill(HIST("QA/Event/EnumEvents"), EventSelection::kVertexCut);
    histos.fill(HIST("QA/Event/Selected/VertexZ"), collision.posZ());

    // LOGF(info, "event id: %d, event multiplicity: %d", collision.bcId(), collision.multTPCtemp());

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsUpperIndexPolicy(dTracks, dTracks))) {

      // Trk1: Kaon, Trk2: Proton
      // pT-dependent PID cut
      // TPC PID cut
      if ((trk1.pt() < 0.2) && (std::abs(trk1.tpcNSigmaKa()) > 4.0))
        continue;
      if ((trk1.pt() >= 0.2) && (trk1.pt() < 0.4) && (std::abs(trk1.tpcNSigmaKa()) > 3.0))
        continue;
      if ((trk1.pt() >= 0.4) && (trk1.pt() < 0.45) && (std::abs(trk1.tpcNSigmaKa()) > 2.0))
        continue;

      if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
        if ((trk1.pt() >= 0.45) && (trk1.pt() < 0.6) && (std::abs(trk1.tpcNSigmaKa()) > 2.0))
          continue;
        if (trk1.pt() > 0.6)
          continue;
      }

      if ((trk2.pt() < 0.25) && (std::abs(trk2.tpcNSigmaPr()) > 4.0))
        continue;
      if ((trk2.pt() >= 0.25) && (trk2.pt() < 0.7) && (std::abs(trk2.tpcNSigmaPr()) > 3.0))
        continue;
      if ((trk2.pt() >= 0.7) && (trk2.pt() < 0.8) && (std::abs(trk2.tpcNSigmaPr()) > 2.0))
        continue;

      if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
        if ((trk2.pt() >= 0.8) && (trk2.pt() < 1.1) && (std::abs(trk2.tpcNSigmaPr()) > 2.0))
          continue;
        if (trk2.pt() > 1.1)
          continue;
      }

      // TOF PID cut
      if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
        if ((trk1.pt() >= 0.1) && (std::abs(trk1.tofNSigmaKa()) > static_cast<float_t>(pidnSigmaPreSelectionCut)))
          continue;
        if ((trk1.pt() >= 0.1) && (std::abs(trk1.tpcNSigmaKa()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF)))
          continue;
        // LOGF(info, "Inside the TOF PID: value for Track 1: %f, PIDflag: %d, check: %d", trk1.tofNSigmaKa(), trk1.tofPIDselectionFlag(), trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF);
      }

      if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
        if ((trk2.pt() >= 0.1) && (std::abs(trk2.tofNSigmaPr()) > static_cast<float_t>(pidnSigmaPreSelectionCut)))
          continue;
        if ((trk2.pt() >= 0.1) && (std::abs(trk2.tpcNSigmaPr()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF)))
          continue;
        // LOGF(info, "Inside the TOF PID: value for Track 2: %f, PIDflag: %d, check: %d", trk2.tofNSigmaPr(), trk2.tofPIDselectionFlag(), trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF);
      }

      //  PID Selection
      //  --- PID QA Kaons
      // if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
      histos.fill(HIST("QA/Kaon/TOF_Nsigma"), trk1.pt(), trk1.tofNSigmaKa());
      histos.fill(HIST("QA/Kaon/TPC_Nsigma"), trk1.pt(), trk1.tpcNSigmaKa());
      histos.fill(HIST("QA/Kaon/TOF_TPC_Map"), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());
      //}
      //  --- PID QA Protons
      // if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
      histos.fill(HIST("QA/Proton/TOF_Nsigma"), trk2.pt(), trk2.tofNSigmaPr());
      histos.fill(HIST("QA/Proton/TPC_Nsigma"), trk2.pt(), trk2.tpcNSigmaPr());
      histos.fill(HIST("QA/Proton/TOF_TPC_Map"), trk2.tofNSigmaPr(), trk2.tpcNSigmaPr());
      // }

      histos.fill(HIST("QA/Kaon/pT"), trk1.pt());
      histos.fill(HIST("QA/Proton/pT"), trk2.pt());

      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPr);
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;

      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        // LOGF(info, "track 1 (K) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk1.globalIndex(), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());  //trk1=Kaon
        // LOGF(info, "track 2 (p) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk2.globalIndex(), trk2.tofNSigmaPr(), trk2.tpcNSigmaPr());  //trk2=Proton
        // LOGF(info, "ResoMass: %f", lResonance.M());
        histos.fill(HIST("Analysis/lambda1520invmass"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmass"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
      } else {
        if (trk1.globalIndex() == trk2.globalIndex())
          continue;
        // LOGF(info, "track 1 (K) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk1.globalIndex(), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());  //trk1=Kaon
        // LOGF(info, "track 2 (p) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk2.globalIndex(), trk2.tofNSigmaPr(), trk2.tpcNSigmaPr());  //trk2=Proton
        // LOGF(info, "ResoMass: %f", lResonance.M());
        histos.fill(HIST("Analysis/lambda1520invmassLS"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassLS"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
      }

      // Like sign pair ++
      if (trk1.sign() > 0 && trk2.sign() > 0) {

        histos.fill(HIST("Analysis/lambda1520invmassLSkp"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassLSkp"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
      }

      // Like sign pair --
      if (trk1.sign() < 0 && trk2.sign() < 0) {
        histos.fill(HIST("Analysis/lambda1520invmassLSkbarpbar"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassLSkbarpbar"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
      }

      if constexpr (IsMC) {
        if (abs(trk1.pdgCode()) != kKPlus || abs(trk2.pdgCode()) != kProton) // check if the tracks are kaons and Protons
          continue;
        auto mother1 = trk1.motherId();
        auto mother2 = trk2.motherId();
        if (mother1 == mother2) {         // Same mother
          if (trk1.motherPDG() == 3124) { // lambda1520
            histos.fill(HIST("MC/reconlambda1520invmass"), lResonance.M());
            histos.fill(HIST("MC/reconlambda1520pt"), lResonance.Pt());
            histos.fill(HIST("MC/h3recolambda1520invmass"), collision.multTPCtemp(), lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processData(aod::ResoCollisions& collisions,
                   aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "[DATA] Processing %d collisions", collisions.size());
    for (auto& collision : collisions) {
      Partition<aod::ResoTracks> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
      fillHistograms<false>(collision, colTracks);
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processData, "Process Event for data", true);

  void processMC(aod::ResoCollisions& collisions,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks, aod::McParticles const& mcParticles)
  {
    LOGF(debug, "[MC] MC events: %d", collisions.size());
    for (auto& collision : collisions) {
      Partition<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex());
      fillHistograms<true>(collision, colTracks);
    }

    // Not related to the real collisions
    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 3124) {         // Lambda(1520)
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          // LOGF(info, "[Rapidity cut] Lambda(1520): %d, y: %f", part.pdgCode(), part.y());
          continue;
        }
        bool pass1 = false;
        bool pass2 = false;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (abs(dau.pdgCode()) == kKPlus) { // Decay to Kaons
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kProton) { // Decay to Protons
            pass1 = true;
          }
        }
        if (!pass1 || !pass2)
          continue;
        histos.fill(HIST("MC/truelambda1520pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processMC, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVetZTPCtemp = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultTPCtemp>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVetZTPCtemp colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVetZTPCtemp> pairs{colBinning, cfgNoMixedEvents, -1, collisions, tracksTuple}; // -1 is the number of the bin to skip

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      Partition<aod::ResoTracks> selectedTracks1 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks1.bindTable(tracks1);
      Partition<aod::ResoTracks> selectedTracks2 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks2.bindTable(tracks2);

      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(selectedTracks1, selectedTracks2))) {
        // Un-like sign pair only
        if (trk1.sign() * trk2.sign() > 0)
          continue;

        // Trk1: Kaon, Trk2: Proton
        // pT-dependent PID cut
        // TPC PID cut
        if ((trk1.pt() < 0.2) && (std::abs(trk1.tpcNSigmaKa()) > 4.0))
          continue;
        if ((trk1.pt() >= 0.2) && (trk1.pt() < 0.4) && (std::abs(trk1.tpcNSigmaKa()) > 3.0))
          continue;
        if ((trk1.pt() >= 0.4) && (trk1.pt() < 0.45) && (std::abs(trk1.tpcNSigmaKa()) > 2.0))
          continue;

        if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
          if ((trk1.pt() >= 0.45) && (trk1.pt() < 0.6) && (std::abs(trk1.tpcNSigmaKa()) > 2.0))
            continue;
          if (trk1.pt() > 0.6)
            continue;
        }

        if ((trk2.pt() < 0.25) && (std::abs(trk2.tpcNSigmaPr()) > 4.0))
          continue;
        if ((trk2.pt() >= 0.25) && (trk2.pt() < 0.7) && (std::abs(trk2.tpcNSigmaPr()) > 3.0))
          continue;
        if ((trk2.pt() >= 0.7) && (trk2.pt() < 0.8) && (std::abs(trk2.tpcNSigmaPr()) > 2.0))
          continue;

        if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) != aod::resodaughter::kHasTOF) {
          if ((trk2.pt() >= 0.8) && (trk2.pt() < 1.1) && (std::abs(trk2.tpcNSigmaPr()) > 2.0))
            continue;
          if (trk2.pt() > 1.1)
            continue;
        }

        // TOF PID cut
        if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
          if ((trk1.pt() >= 0.1) && (std::abs(trk1.tofNSigmaKa()) > static_cast<float_t>(pidnSigmaPreSelectionCut)))
            continue;
          if ((trk1.pt() >= 0.1) && (std::abs(trk1.tpcNSigmaKa()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF)))
            continue;
          // LOGF(info, "Inside the TOF PID: value for Track 1: %f, PIDflag: %d, check: %d", trk1.tofNSigmaKa(), trk1.tofPIDselectionFlag(), trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF);
        }

        if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
          if ((trk2.pt() >= 0.1) && (std::abs(trk2.tofNSigmaPr()) > static_cast<float_t>(pidnSigmaPreSelectionCut)))
            continue;
          if ((trk2.pt() >= 0.1) && (std::abs(trk2.tpcNSigmaPr()) > static_cast<float_t>(pidnSigmaTPCCuthasTOF)))
            continue;
          // LOGF(info, "Inside the TOF PID: value for Track 2: %f, PIDflag: %d, check: %d", trk2.tofNSigmaPr(), trk2.tofPIDselectionFlag(), trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF);
        }

        lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
        lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPr);
        lResonance = lDecayDaughter1 + lDecayDaughter2;

        if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
          continue;
        // LOGF(info, "track 1 (K) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk1.globalIndex(), trk1.tofNSigmaKa(), trk1.tpcNSigmaKa());  //trk1=Kaon
        // LOGF(info, "track 2 (p) ResoID: %d, NsigTOF: %f NSigTPC: %f", trk2.globalIndex(), trk2.tofNSigmaPr(), trk2.tpcNSigmaPr());  //trk2=Proton
        // LOGF(info, "ResoMass: %f", lResonance.M());

        histos.fill(HIST("Analysis/lambda1520invmassME"), lResonance.M());
        histos.fill(HIST("Analysis/h3lambda1520invmassME"), collision1.multTPCtemp(), lResonance.Pt(), lResonance.M());
      }
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processME, "Process EventMixing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambda1520analysis>(cfgc, TaskName{"lf-lambda1520analysis"})};
}
