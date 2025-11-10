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

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TLorentzVector.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct phianalysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  ///// Configurables
  /// Histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Binning of the centrality axis"};
  Configurable<float> cInvMassStart{"cInvMassStart", 0.6, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 1.5, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};
  Configurable<int> cPIDBins{"cPIDBins", 65, "PID binning"};
  Configurable<float> cPIDQALimit{"cPIDQALimit", 6.5, "PID QA limit"};
  Configurable<int> cDCABins{"cDCABins", 150, "DCA binning"};
  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0., 1., 5., 10., 30., 50., 70., 100., 110.}, "Mixing bins - multiplicity"};
  /// Pre-selection cuts
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  /// DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for PID selection
  /// TPC nCluster cut
  Configurable<int> cMinTPCNclsFound{"cMinTPCNclsFound", 70, "Minimum TPC cluster found"};
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", 3.0, "Combined nSigma cut for Kaon"}; // Combined
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "By pass TOF PID selection"};                          // By pass TOF PID selection
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "V0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};

    // Mass QA (quick check)
    histos.add("phiinvmassDS", "Invariant mass of phi differnt sign", kTH1F, {invMassAxis});
    histos.add("phiinvmassLS", "Invariant mass of phi like sign", kTH1F, {invMassAxis});
    histos.add("phiinvmassME", "Invariant mass of phi mixed event", kTH1F, {invMassAxis});
    // DCA QA
    histos.add("QAbefore/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    // pT QA
    histos.add("QAbefore/trkpT", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT", "pT distribution of kaon track candidates", kTH1F, {ptAxis});
    // PID QA before cuts
    histos.add("QAbefore/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QAbefore/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {ptAxisQA, pidQAAxis}});
    histos.add("QAbefore/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {ptAxisQA, pidQAAxis}});
    // PID QA after cuts
    histos.add("QAafter/TOF_TPC_Mapka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QAafter/TOF_Nsigma_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {ptAxisQA, pidQAAxis}});
    histos.add("QAafter/TPC_Nsigmaka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {ptAxisQA, pidQAAxis}});

    // 3d histogram
    histos.add("h3phiinvmassDS", "Invariant mass of phi differnt sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3phiinvmassLS", "Invariant mass of phi same sign", kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("h3phiinvmassME", "Invariant mass of phi mixed event", kTH3F, {centAxis, ptAxis, invMassAxis});

    if (doprocessMCLight) {
      // MC QA
      histos.add("QAMCTrue/trkDCAxy", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAz", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
      histos.add("h3Recphiinvmass", "Invariant mass of Reconstructed MC phi", kTH3F, {centAxis, ptAxis, invMassAxis});
      histos.add("phiGen", "pT distribution of True MC phi", kTH2F, {ptAxis, centAxis});
      histos.add("phiRec", "pT distribution of Reconstructed MC phi", kTH2F, {ptAxis, centAxis});
      histos.add("phiRecinvmass", "Inv mass distribution of Reconstructed MC Phi", kTH1F, {invMassAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in phi analysis:";
    histos.print();
  }

  double massKa = MassKaonCharged;

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsFound() < cMinTPCNclsFound)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }
  // PID selection tool
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
      tpcPIDPassed = true;
    }
    if (cByPassTOF && tpcPIDPassed) {
      return true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaKa()) < cMaxTOFnSigmaKaon) {
        tofPIDPassed = true;
      }
      if ((nsigmaCutCombinedKaon > 0) && (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa() < nsigmaCutCombinedKaon * nsigmaCutCombinedKaon)) {
        tofPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    auto multiplicity = collision.cent();
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsUpperIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      //// Initialize variables
      // Trk1: Kaon, Trk2: Kaon
      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();
      auto trk1ptKa = trk1.pt();
      auto trk1NSigmaKaTPC = trk1.tpcNSigmaKa();
      auto trk1NSigmaKaTOF = (isTrk1hasTOF) ? trk1.tofNSigmaKa() : -999.;
      // auto trk2ptKa = trk2.pt();
      // auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      // auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      if constexpr (!IsMix) {
        //// QA plots before the selection
        //  --- PID QA Kaon
        histos.fill(HIST("QAbefore/TPC_Nsigmaka_all"), trk1ptKa, trk1NSigmaKaTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QAbefore/TOF_Nsigma_all"), trk1ptKa, trk1NSigmaKaTOF);
          histos.fill(HIST("QAbefore/TOF_TPC_Mapka_all"), trk1NSigmaKaTOF, trk1NSigmaKaTPC);
        }
        histos.fill(HIST("QAbefore/trkpT"), trk1ptKa);
        histos.fill(HIST("QAbefore/trkDCAxy"), trk1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAz"), trk1.dcaZ());
      }

      //// Apply the selection
      if (cUseOnlyTOFTrackKa && (!isTrk1hasTOF || !isTrk2hasTOF))
        continue;
      if (!selectionPIDKaon(trk1) || !selectionPIDKaon(trk2))
        continue;

      if constexpr (!IsMix) {
        //// QA plots after the selection
        //  --- PID QA Kaon
        histos.fill(HIST("QAafter/TPC_Nsigmaka_all"), trk1ptKa, trk1NSigmaKaTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QAafter/TOF_Nsigma_all"), trk1ptKa, trk1NSigmaKaTOF);
          histos.fill(HIST("QAafter/TOF_TPC_Mapka_all"), trk1NSigmaKaTOF, trk1NSigmaKaTPC);
        }
        histos.fill(HIST("QAafter/trkpT"), trk1ptKa);
        histos.fill(HIST("QAafter/trkDCAxy"), trk1.dcaXY());
        histos.fill(HIST("QAafter/trkDCAz"), trk1.dcaZ());
      }

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massKa);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (abs(lResonance.Rapidity()) > 0.5)
        continue;
      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          if (trk1.sign() > 0) {
            histos.fill(HIST("phiinvmassDS"), lResonance.M());
            histos.fill(HIST("h3phiinvmassDS"), multiplicity, lResonance.Pt(), lResonance.M());
          } else {
          }
        } else {
          histos.fill(HIST("phiinvmassME"), lResonance.M());
          histos.fill(HIST("h3phiinvmassME"), multiplicity, lResonance.Pt(), lResonance.M());
        }

        // MC
        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != 321 || abs(trk2.pdgCode()) != 321)
            continue;
          if (trk1.motherId() != trk2.motherId()) // Same mother
            continue;
          if (abs(trk1.motherPDG()) != 333)
            continue;

          // Track selection check.
          histos.fill(HIST("QAMCTrue/trkDCAxy"), trk2.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAz"), trk2.dcaZ());

          // MC histograms
          histos.fill(HIST("phiRec"), lResonance.Pt(), multiplicity);
          histos.fill(HIST("phiRecinvmass"), lResonance.M());
          histos.fill(HIST("h3Recphiinvmass"), multiplicity, lResonance.Pt(), lResonance.M());
        }
      } else {
        if constexpr (!IsMix)
          continue;
        if (trk1.sign() > 0) {
          histos.fill(HIST("phiinvmassLS"), lResonance.M());
          histos.fill(HIST("h3phiinvmassLS"), multiplicity, lResonance.Pt(), lResonance.M());
        } else {
        }
      }
    }
  }

  void processDataLight(aod::ResoCollision& collision,
                        aod::ResoTracks const& resotracks)
  {
    // LOG(info) << "new collision, zvtx: " << collision.posZ();
    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(phianalysis, processDataLight, "Process Event for data", false);

  void processMCLight(aod::ResoCollision& collision,
                      soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(phianalysis, processMCLight, "Process Event for MC", false);

  void processMCTrue(aod::ResoCollision& collision, aod::ResoMCParents& resoParents)
  {
    auto multiplicity = collision.cent();
    for (auto& part : resoParents) {  // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 333) // phi(0)
        continue;
      if (abs(part.y()) > 0.5) { // rapidity cut
        continue;
      }
      if (abs(part.daughterPDG1()) != 321 || abs(part.daughterPDG2()) != 321) { // At least one decay to Kaon
        continue;
      }
      histos.fill(HIST("phiGen"), part.pt(), multiplicity);
    }
  }
  PROCESS_SWITCH(phianalysis, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMELight(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(phianalysis, processMELight, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phianalysis>(cfgc, TaskName{"lf-phianalysis"})};
}
