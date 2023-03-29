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

/// \file k892analysis.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

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

struct k892analysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;

  framework::Service<o2::ccdb::BasicCCDBManager> ccdb; /// Accessing the CCDB
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::QAObject};
  // Configurables

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};

  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

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

    // Mass QA (quick check)
    histos.add("k892invmass", "Invariant mass of K(892)0", kTH1F, {{900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("k892invmassME", "Invariant mass of K(892)0 mixed event", kTH1F, {{900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"}});
    histos.add("trk1pT", "pT distribution of track1", kTH1F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("trk2pT", "pT distribution of track1", kTH1F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    histos.add("TOF_TPC_Map1", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_Nsigma1", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TPC_Nsigma1", "TPC NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_TPC_Map2", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TOF_Nsigma2", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});
    histos.add("TPC_Nsigma2", "TPC NSigma for Kaons;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {{1000, -10, 10}, {1000, -10, 10}}});

    // 3d histogram
    histos.add("h3k892invmass", "Invariant mass of K(892)0", kTH3F, {{300, 0, 3000}, {100, 0.0f, 10.0f}, {900, 0.6, 1.5}});
    histos.add("h3k892invmassME", "Invariant mass of K(892)0 mixed event", kTH3F, {{300, 0, 3000}, {100, 0.0f, 10.0f}, {900, 0.6, 1.5}});

    if (doprocessMC) {
      histos.add("h3recok892invmass", "Invariant mass of Reconstructed MC K(892)0", kTH3F, {{300, 0, 3000}, {100, 0.0f, 10.0f}, {900, 0.6, 1.5}});
      histos.add("truek892pt", "pT distribution of True MC K(892)0", kTH1F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      histos.add("reconk892pt", "pT distribution of Reconstructed MC K(892)0", kTH1F, {{100, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
      histos.add("reconk892invmass", "Inv mass distribution of Reconstructed MC Phi", kTH1F, {{900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"}});
    }
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks, dTracks))) {
      // Un-like sign pair only
      if (trk1.sign() * trk2.sign() > 0)
        continue;

      // Trk1: Pion, Trk2: Kaon
      // pT-dependent TPC PID cut
      if ((trk1.pt() < 0.3) && (std::abs(trk1.tpcNSigmaPi()) > 6.0))
        continue;
      if ((trk1.pt() >= 0.3) && (trk1.pt() < 0.4) && (std::abs(trk1.tpcNSigmaPi()) > 4.0))
        continue;
      if ((trk1.pt() >= 0.4) && (std::abs(trk1.tpcNSigmaPi()) > 2.0))
        continue;

      if ((trk2.pt() < 0.3) && (std::abs(trk2.tpcNSigmaKa()) > 6.0))
        continue;
      if ((trk2.pt() >= 0.3) && (trk2.pt() < 0.4) && (std::abs(trk2.tpcNSigmaKa()) > 4.0))
        continue;
      if ((trk2.pt() >= 0.4) && (std::abs(trk2.tpcNSigmaKa()) > 2.0))
        continue;

      //  --- PID QA Pion -
      histos.fill(HIST("TOF_Nsigma1"), trk1.pt(), trk1.tofNSigmaPi());
      histos.fill(HIST("TPC_Nsigma1"), trk1.pt(), trk1.tpcNSigmaPi());
      histos.fill(HIST("TOF_TPC_Map1"), trk1.tofNSigmaPi(), trk1.tpcNSigmaPi());
      //  --- PID QA Kaon +
      histos.fill(HIST("TOF_Nsigma2"), trk2.pt(), trk2.tofNSigmaKa());
      histos.fill(HIST("TPC_Nsigma2"), trk2.pt(), trk2.tpcNSigmaKa());
      histos.fill(HIST("TOF_TPC_Map2"), trk2.tofNSigmaKa(), trk2.tpcNSigmaKa());

      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);

      histos.fill(HIST("trk1pT"), trk1.pt());
      histos.fill(HIST("trk2pT"), trk2.pt());

      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;

      histos.fill(HIST("k892invmass"), lResonance.M());
      histos.fill(HIST("h3k892invmass"), collision.multV0M(), lResonance.Pt(), lResonance.M());

      if constexpr (IsMC) {
        if (abs(trk1.pdgCode()) != kPiPlus || abs(trk2.pdgCode()) != kKPlus)
          continue;
        auto mother1 = trk1.motherId();
        auto mother2 = trk2.motherId();
        if (mother1 == mother2) {             // Same mother
          if (abs(trk1.motherPDG()) == 313) { // k892(0)
            histos.fill(HIST("reconk892pt"), lResonance.Pt());
            histos.fill(HIST("reconk892invmass"), lResonance.M());
            histos.fill(HIST("h3recok892invmass"), collision.multV0M(), lResonance.Pt(), lResonance.M());
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
      Partition<aod::ResoTracks> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<false>(collision, colTracks);
    }
  }
  PROCESS_SWITCH(k892analysis, processData, "Process Event for data", true);

  void processMC(aod::ResoCollisions& collisions,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks, aod::McParticles const& mcParticles)
  {
    LOGF(debug, "[MC] MC events: %d", collisions.size());
    for (auto& collision : collisions) {
      Partition<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<true>(collision, colTracks);
    }

    // Not related to the real collisions
    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 313) {          // K892(0)
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          // LOGF(info, "[Rapidity cut] K892(0): %d, y: %f", part.pdgCode(), part.y());
          continue;
        }
        bool pass1 = false;
        bool pass2 = false;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (abs(dau.pdgCode()) == kKPlus) { // Decay to Kaons
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kPiPlus) { // Decay to Pions
            pass1 = true;
          }
        }
        if (!pass1 || !pass2)
          continue;
        histos.fill(HIST("truek892pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(k892analysis, processMC, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVetZTPCtemp = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultTPCtemp>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVetZTPCtemp colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVetZTPCtemp> pairs{colBinning, 10, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      Partition<aod::ResoTracks> selectedTracks1 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks1.bindTable(tracks1);
      Partition<aod::ResoTracks> selectedTracks2 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks2.bindTable(tracks2);

      for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(selectedTracks1, selectedTracks2))) {
        // Un-like sign pair only
        if (trk1.sign() * trk2.sign() > 0)
          continue;
        if ((trk1.pt() < 0.3) && (std::abs(trk1.tpcNSigmaPi()) > 6.0))
          continue;
        if ((trk1.pt() >= 0.3) && (trk1.pt() < 0.4) && (std::abs(trk1.tpcNSigmaPi()) > 4.0))
          continue;
        if ((trk1.pt() >= 0.4) && (std::abs(trk1.tpcNSigmaPi()) > 2.0))
          continue;

        if ((trk2.pt() < 0.3) && (std::abs(trk2.tpcNSigmaKa()) > 6.0))
          continue;
        if ((trk2.pt() >= 0.3) && (trk2.pt() < 0.4) && (std::abs(trk2.tpcNSigmaKa()) > 4.0))
          continue;
        if ((trk2.pt() >= 0.4) && (std::abs(trk2.tpcNSigmaKa()) > 2.0))
          continue;

        lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
        lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
        lResonance = lDecayDaughter1 + lDecayDaughter2;

        if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
          continue;

        histos.fill(HIST("k892invmassME"), lResonance.M());
        histos.fill(HIST("h3k892invmassME"), collision1.multV0M(), lResonance.Pt(), lResonance.M());
      }
    }
  };
  PROCESS_SWITCH(k892analysis, processME, "Process EventMixing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892analysis>(cfgc, TaskName{"lf-k892analysis"})};
}
