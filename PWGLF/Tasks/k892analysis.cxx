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
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ///// Configurables
  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};
  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 2.0, "TPC nSigma cut for Pion"}; // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 2.0, "TOF nSigma cut for Pion"}; // TOF
  // Kaon
  Configurable<std::vector<double>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100.};
    AxisSpec centAxis = {centBinning, "V0M (%)"};
    AxisSpec multAxis = {0, 0, 100, "V0M (%)"}; // for future
    AxisSpec ptAxis = {200, 0, 20, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {130, -6.5, 6.5};

    // Mass QA (quick check)
    histos.add("k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTH1F, {invMassAxis});
    histos.add("k892invmassLS", "Invariant mass of K(892)0 like sign", kTH1F, {invMassAxis});
    histos.add("k892invmassME", "Invariant mass of K(892)0 mixed event", kTH1F, {invMassAxis});
    // pT QA
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAbefore/trkpT_ka", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_ka", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_ka_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_ka_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});

    // 3d histogram
    histos.add("h3k892invmassDS", "Invariant mass of K(892)0 differnt sign", kTH3F, {{300, 0, 3000}, ptAxis, invMassAxis}); // TODO: multiplicity bin has to be updatde.
    histos.add("h3k892invmassLS", "Invariant mass of K(892)0 same sign", kTH3F, {{300, 0, 3000}, ptAxis, invMassAxis});     // TODO: multiplicity bin has to be updatde.
    histos.add("h3k892invmassME", "Invariant mass of K(892)0 mixed event", kTH3F, {{300, 0, 3000}, ptAxis, invMassAxis});   // TODO: multiplicity bin has to be updatde.

    if (doprocessMC) {
      // MC QA
      histos.add("h3recok892invmass", "Invariant mass of Reconstructed MC K(892)0", kTH3F, {{300, 0, 3000}, ptAxis, invMassAxis}); // TODO: multiplicity bin has to be updatde.
      histos.add("truek892pt", "pT distribution of True MC K(892)0", kTH1F, {ptAxis});
      histos.add("reconk892pt", "pT distribution of Reconstructed MC K(892)0", kTH1F, {ptAxis});
      histos.add("reconk892invmass", "Inv mass distribution of Reconstructed MC Phi", kTH1F, {invMassAxis});
    }
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    bool isTrk1Selected{true}, isTrk2Selected{true}, isTrk1hasTOF{false}, isTrk2hasTOF{false};
    auto vKaonTPCPIDpTintv = static_cast<std::vector<double>>(kaonTPCPIDpTintv);
    auto vKaonTPCPIDcuts = static_cast<std::vector<double>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<double>>(kaonTOFPIDpTintv);
    auto vKaonTOFPIDcuts = static_cast<std::vector<double>>(kaonTOFPIDcuts);
    auto lengthOfkaonTPCPIDpTintv = static_cast<int>(vKaonTPCPIDpTintv.size());
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      //// Initialize variables
      // Trk1: Pion, Trk2: Kaon
      isTrk1Selected = true;
      isTrk2Selected = true;
      isTrk1hasTOF = false;
      isTrk2hasTOF = false;
      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = -999.;

      // hasTOF?
      if ((trk1.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
        isTrk1hasTOF = true;
      }
      if ((trk2.tofPIDselectionFlag() & aod::resodaughter::kHasTOF) == aod::resodaughter::kHasTOF) {
        isTrk2hasTOF = true;
      }
      //// PID selections
      // For Pion candidate, we don't need to apply pT-dependent PID cuts
      if (std::abs(trk1NSigmaPiTPC) > cMaxTPCnSigmaPion)
        isTrk1Selected = false;
      if (isTrk1hasTOF) {
        trk1NSigmaPiTOF = trk1.tofNSigmaPi();
        if (std::abs(trk1NSigmaPiTOF) > cMaxTOFnSigmaPion)
          isTrk1Selected = false;
      }
      // For Kaon candidate, we need to apply pT-dependent PID cuts
      if (lengthOfkaonTPCPIDpTintv > 0) {
        for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
          if (trk2ptKa < vKaonTPCPIDpTintv[i]) {
            if (std::abs(trk2NSigmaKaTPC) > vKaonTPCPIDcuts[i])
              isTrk2Selected = false;
          }
        }
      }
      if (isTrk2hasTOF) {
        trk2NSigmaKaTOF = trk2.tofNSigmaKa();
        if (lengthOfkaonTPCPIDpTintv > 0) {
          for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
            if (trk2ptKa < vKaonTOFPIDpTintv[i]) {
              if (std::abs(trk2NSigmaKaTOF) > vKaonTOFPIDcuts[i])
                isTrk2Selected = false;
            }
          }
        }
      }

      //// QA plots before the selection
      //  --- PID QA Pion
      histos.fill(HIST("QAbefore/TPC_Nsigma_pi_all"), trk1ptPi, trk1NSigmaPiTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAbefore/TOF_Nsigma_pi_all"), trk1ptPi, trk1NSigmaPiTOF);
        histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
      }
      //  --- PID QA Kaon
      histos.fill(HIST("QAbefore/TPC_Nsigmaka_all"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAbefore/TOF_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }
      histos.fill(HIST("QAbefore/trkpT_pi"), trk1ptPi);
      histos.fill(HIST("QAbefore/trkpT_ka"), trk2ptKa);

      //// Apply the selection
      if (!isTrk1Selected || !isTrk2Selected)
        continue;

      //// QA plots after the selection
      //  --- PID QA Pion
      histos.fill(HIST("QAafter/TPC_Nsigma_pi_all"), trk1ptPi, trk1NSigmaPiTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAafter/TOF_Nsigma_pi_all"), trk1ptPi, trk1NSigmaPiTOF);
        histos.fill(HIST("QAafter/TOF_TPC_Map_pi_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
      }
      //  --- PID QA Kaon
      histos.fill(HIST("QAafter/TPC_Nsigmaka_all"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QAafter/TOF_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }
      histos.fill(HIST("QAafter/trkpT_pi"), trk1ptPi);
      histos.fill(HIST("QAafter/trkpT_ka"), trk2ptKa);

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;
      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          histos.fill(HIST("k892invmassDS"), lResonance.M());
          histos.fill(HIST("h3k892invmassDS"), collision.multV0M(), lResonance.Pt(), lResonance.M()); // TODO: multV0M has to be updatde.
        } else {
          histos.fill(HIST("k892invmassME"), lResonance.M());
          histos.fill(HIST("h3k892invmassME"), collision.multV0M(), lResonance.Pt(), lResonance.M()); // TODO: multV0M has to be updatde.
        }

        // MC
        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != kPiPlus || abs(trk2.pdgCode()) != kKPlus)
            continue;
          auto mother1 = trk1.motherId();
          auto mother2 = trk2.motherId();
          if (mother1 == mother2) {             // Same mother
            if (abs(trk1.motherPDG()) == 313) { // k892(0)
              histos.fill(HIST("reconk892pt"), lResonance.Pt());
              histos.fill(HIST("reconk892invmass"), lResonance.M());
              histos.fill(HIST("h3recok892invmass"), collision.multV0M(), lResonance.Pt(), lResonance.M()); // TODO: multV0M has to be updatde.
            }
          }
        }
      } else {
        if constexpr (!IsMix) {
          histos.fill(HIST("k892invmassLS"), lResonance.M());
          histos.fill(HIST("h3k892invmassLS"), collision.multV0M(), lResonance.Pt(), lResonance.M()); // TODO: multV0M has to be updatde.
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
      fillHistograms<false, false>(collision, colTracks, colTracks);
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
      fillHistograms<true, false>(collision, colTracks, colTracks);
    }

    // Not related to the real collisions
    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 313) {          // K892(0)
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          continue;
        }
        bool pass1 = false;
        bool pass2 = false;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (abs(dau.pdgCode()) == kKPlus) { // At least one decay to Kaon
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kPiPlus) { // At least one decay to Pion
            pass1 = true;
          }
        }
        if (!pass1 || !pass2) // If we have both decay products
          continue;
        histos.fill(HIST("truek892pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(k892analysis, processMC, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVetZTPCtemp = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultTPCtemp>; // TODO: MultTPCtemp has to be updatde.
  BinningTypeVetZTPCtemp colBinning{{CfgVtxBins, CfgMultBins}, true};
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(resotracks);
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVetZTPCtemp> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      Partition<aod::ResoTracks> selectedTracks1 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks1.bindTable(tracks1);
      Partition<aod::ResoTracks> selectedTracks2 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)); // Basic DCA cuts
      selectedTracks2.bindTable(tracks2);

      fillHistograms<false, true>(collision1, selectedTracks1, selectedTracks2);
    }
  };
  PROCESS_SWITCH(k892analysis, processME, "Process EventMixing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k892analysis>(cfgc, TaskName{"lf-k892analysis"})};
}
