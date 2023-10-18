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

/// \file k1analysis.cxx
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

struct k1analysis {
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
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  /// PID Selections
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};              // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};              // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};  // Combined
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
  // Kaon
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon"};              // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon"};              // TOF
  Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"}; // Combined
  // Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // bachelor pion TOF PID?
  Configurable<int> cDoTOFPID{"cDoTOFPID", 1, "Do TOF PID"};

  // K(892)0 selection
  Configurable<double> cK892masswindow{"cK892masswindow", 0.1, "K(892)0 inv mass selection window"};
  Configurable<double> cPiPiMin{"cPiPiMin", 0, "Pion pair inv mass selection minimum"};
  Configurable<double> cPiPiMax{"cPiPiMax", 999, "Pion pair inv mass selection maximum"};
  Configurable<double> cPiKaMin{"cPiKaMin", 0, "bPion-Kaon pair inv mass selection minimum"};
  Configurable<double> cPiKaMax{"cPiKaMax", 999, "bPion-Kaon pair inv mass selection maximum"};
  Configurable<double> cMinAngle{"cMinAngle", 0, "Minimum angle between K(892)0 and bachelor pion"};
  Configurable<double> cMaxAngle{"cMaxAngle", 4, "Maximum angle between K(892)0 and bachelor pion"};
  Configurable<double> cMinPairAsym{"cMinPairAsym", -1, "Minimum pair asymmetry"};
  Configurable<double> cMaxPairAsym{"cMaxPairAsym", 1, "Maximum pair asymmetry"};

  // K1 selection
  Configurable<double> cK1MaxRap{"cK1MaxRap", 0.5, "K1 maximum rapidity"};
  Configurable<double> cK1MinRap{"cK1MinRap", -0.5, "K1 minimum rapidity"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {300, 0, 3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {500, 0, 5, "DCA_{#it{xy}} (cm)"};
    AxisSpec invMassAxis = {900, 0.6, 1.5, "Invariant Mass (GeV/#it{c}^2)"};        // K(892)0
    AxisSpec invMassAxisReso = {1600, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec invMassAxisScan = {250, 0, 2.5, "Invariant Mass (GeV/#it{c}^2)"};      // For selection
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};

    // Mass QA (quick check)
    histos.add("k892invmass", "Invariant mass of K(892)0", HistType::kTH1F, {invMassAxis});
    histos.add("k1invmass", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_LS", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_Mix", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    if (doprocessMC) {
      histos.add("k1invmass_MC", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    }
    // DCA QA
    histos.add("QA/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkDCAxy_ka", "DCAxy distribution of kaon track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkDCAxy_pi_bach", "DCAxy distribution of bachelor pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkDCAz_ka", "DCAz distribution of kaon track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkDCAz_pi_bach", "DCAz distribution of bachelor pion track candidates", HistType::kTH1F, {dcazAxis});

    // pT QA
    histos.add("QA/trkpT_pi", "pT distribution of pion track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkpT_ka", "pT distribution of kaon track candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkpT_pi_bach", "pT distribution of bachelor pion track candidates", HistType::kTH1F, {ptAxis});
    // PID QA after cuts
    histos.add("QA/TOF_TPC_Map_pi", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QA/TOF_Nsigma_pi", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/TPC_Nsigma_pi", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/TOF_TPC_Map_ka", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QA/TOF_Nsigma_ka", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/TPC_Nsigmaka", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/TOF_TPC_Map_pi_bach", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
    histos.add("QA/TOF_Nsigma_pi_bach", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/TPC_Nsigma_pi_bach", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH2F, {ptAxis, pidQAAxis}});
    histos.add("QA/InvMass_piK_pipi", "Invariant mass of pion + kaon and pion+pion;Invariant Mass (GeV/#it{c}^{2});Invariant Mass (GeV/#it{c}^{2});", {HistType::kTH2F, {invMassAxisScan, invMassAxisScan}});
    histos.add("QA/InvMass_piK_pika", "Invariant mass of pion + kaon and pion+kaon;Invariant Mass (GeV/#it{c}^{2});Invariant Mass (GeV/#it{c}^{2});", {HistType::kTH2F, {invMassAxisScan, invMassAxisScan}});
    histos.add("QA/K1OA", "Opening angle of K1(1270)pm", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle of K1(1270)pm"}});
    histos.add("QA/K1PairAsymm", "Pair asymmetry of K1(1270)pm", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry of K1(1270)pm"}});

    // Invariant mass histograms
    histos.add("hK892invmass_PP", "Invariant mass of K(892)0 (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_NP", "Invariant mass of K(892)0 (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_PN", "Invariant mass of K(892)0 (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK892invmass_NN", "Invariant mass of K(892)0 (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxis});
    histos.add("hK1invmass_NPP", "Invariant mass of K(892)0 + pion (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_NPN", "Invariant mass of K(892)0 + pion (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_PNP", "Invariant mass of K(892)0 + pion (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_PNN", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    // K892-LS bkg
    histos.add("hK1invmass_PPP", "Invariant mass of K(892)0 + pion (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_PPN", "Invariant mass of K(892)0 + pion (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_NNP", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_NNN", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    // Mixed event
    histos.add("hK1invmass_NPP_Mix", "Invariant mass of K(892)0 + pion (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_NPN_Mix", "Invariant mass of K(892)0 + pion (Matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_PNP_Mix", "Invariant mass of K(892)0 + pion (Anti-matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
    histos.add("hK1invmass_PNN_Mix", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});

    if (doprocessMC) {
      histos.add("QAMC/InvMass_piK_pipi", "Invariant mass of pion + kaon and pion+pion;Invariant Mass (GeV/#it{c}^{2});Invariant Mass (GeV/#it{c}^{2});", {HistType::kTH2F, {invMassAxisScan, invMassAxisScan}});
      histos.add("QAMC/InvMass_piK_pika", "Invariant mass of pion + kaon and pion+kaon;Invariant Mass (GeV/#it{c}^{2});Invariant Mass (GeV/#it{c}^{2});", {HistType::kTH2F, {invMassAxisScan, invMassAxisScan}});
      histos.add("QAMC/K1OA", "Opening angle of K1(1270)pm", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle of K1(1270)pm"}});
      histos.add("QAMC/K1PairAsymm", "Pair asymmetry of K1(1270)pm", HistType::kTH1F, {AxisSpec{100, 0, 1, "Pair asymmetry of K1(1270)pm"}});

      histos.add("hK1invmass_NPP_MC", "Invariant mass of K(892)0 + pion (Matter + Matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});
      histos.add("hK1invmass_PNN_MC", "Invariant mass of K(892)0 + pion (Anti-matter + Anti-matter)", HistType::kTH3F, {centAxis, ptAxis, invMassAxisReso});

      histos.add("hReconK892pt", "pT distribution of Reconstructed MC K(892)0", HistType::kTH1F, {ptAxis});
      histos.add("hTrueK1pt", "pT distribution of True MC K1", HistType::kTH1F, {ptAxis});
      histos.add("hReconK1pt", "pT distribution of Reconstructed MC K1", HistType::kTH1F, {ptAxis});

      histos.add("k1invmass_noK1", "Invariant mass of K1(1270)pm", HistType::kTH1F, {invMassAxisReso});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in spectraTOF";
    histos.print();
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();
  double massK892 = TDatabasePDG::Instance()->GetParticle(313)->Mass();

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cMinPtcut)
      return false;
    if (track.dcaXY() > cMaxDCArToPVcut)
      return false;
    if (track.dcaZ() < cMinDCAzToPVcut || track.dcaZ() > cMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    return true;
  }

  // PID selection tools
  template <typename T>
  bool selectionPIDPion(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
      tpcPIDPassed = true;
    }
    if (candidate.hasTOF()) {
      if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPion) {
        tofPIDPassed = true;
      }
      if ((nsigmaCutCombinedPion > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
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
  template <typename T>
  bool selectionPIDKaon(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    if (std::abs(candidate.tpcNSigmaKa()) < cMaxTPCnSigmaKaon) {
      tpcPIDPassed = true;
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

  template <typename T, typename T2>
  bool isTrueK1(const T& trk1, const T& trk2, const T2& bTrack)
  {
    if (abs(trk1.pdgCode()) != kPiPlus || abs(trk2.pdgCode()) != kKPlus)
      return false;
    auto mother1 = trk1.motherId();
    auto mother2 = trk2.motherId();
    if (mother1 != mother2)
      return false;
    if (abs(trk1.motherPDG()) != 313)
      return false;
    if (abs(bTrack.pdgCode()) != kPiPlus)
      return false;
    if (abs(bTrack.motherPDG()) != 10323)
      return false;
    auto siblings = bTrack.siblingIds();
    if (siblings[0] != mother1 && siblings[1] != mother1)
      return false;

    return true;
  }

  template <typename T>
  bool isTrueK892(const T& trk1, const T& trk2)
  {
    if (abs(trk1.pdgCode()) != kPiPlus || abs(trk2.pdgCode()) != kKPlus)
      return false;
    auto mother1 = trk1.motherId();
    auto mother2 = trk2.motherId();
    if (mother1 != mother2)
      return false;
    if (abs(trk1.motherPDG()) != 313)
      return false;
    return true;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonanceK892, lDecayDaughter_bach, lResonanceK1;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks2, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.
      // Trk1: Pion, Trk2: Kaon
      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      auto isTrk1hasTOF = trk1.hasTOF();
      auto isTrk2hasTOF = trk2.hasTOF();
      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;
      auto trk2ptKa = trk2.pt();
      auto trk2NSigmaKaTPC = trk2.tpcNSigmaKa();
      auto trk2NSigmaKaTOF = (isTrk2hasTOF) ? trk2.tofNSigmaKa() : -999.;

      //// PID selections
      if (cUseOnlyTOFTrackPi && !isTrk1hasTOF)
        continue;
      if (cUseOnlyTOFTrackKa && !isTrk2hasTOF)
        continue;
      if (!selectionPIDPion(trk1) || !selectionPIDKaon(trk2))
        continue;

      //// QA plots after the selection
      //  --- PID QA Pion
      if constexpr (!IsMix) {
        histos.fill(HIST("QA/TPC_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/TOF_Nsigma_pi"), trk1ptPi, trk1NSigmaPiTOF);
          histos.fill(HIST("QA/TOF_TPC_Map_pi"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
        }
        //  --- PID QA Kaon
        histos.fill(HIST("QA/TPC_Nsigmaka"), trk2ptKa, trk2NSigmaKaTPC);
        if (isTrk1hasTOF) {
          histos.fill(HIST("QA/TOF_Nsigma_ka"), trk2ptKa, trk2NSigmaKaTOF);
          histos.fill(HIST("QA/TOF_TPC_Map_ka"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
        }
        histos.fill(HIST("QA/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QA/trkpT_ka"), trk2ptKa);

        histos.fill(HIST("QA/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QA/trkDCAxy_ka"), trk2.dcaXY());
        histos.fill(HIST("QA/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QA/trkDCAz_ka"), trk2.dcaZ());
      }

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonanceK892 = lDecayDaughter1 + lDecayDaughter2;

      if constexpr (!IsMix) {
        histos.fill(HIST("k892invmass"), lResonanceK892.M()); // quick check
        if (trk1.sign() > 0) {                                // Positive pion
          if (trk2.sign() > 0)                                // Positive kaon
            histos.fill(HIST("hK892invmass_PP"), collision.multV0M(), lResonanceK892.Pt(), lResonanceK892.M());
          else                                                                                                  // Negative kaon
            histos.fill(HIST("hK892invmass_PN"), collision.multV0M(), lResonanceK892.Pt(), lResonanceK892.M()); // Anti-K(892)0
        } else {                                                                                                // Negative pion
          if (trk2.sign() > 0)                                                                                  // Positive kaon
            histos.fill(HIST("hK892invmass_NP"), collision.multV0M(), lResonanceK892.Pt(), lResonanceK892.M()); // K(892)0
          else                                                                                                  // Negative kaon
            histos.fill(HIST("hK892invmass_NN"), collision.multV0M(), lResonanceK892.Pt(), lResonanceK892.M());
        }
      }
      // Like-sign rejection for K(892)0 - disabled for further LS bkg study
      // if (trk1.sign() * trk2.sign() > 0)
      //   continue;

      if constexpr (IsMC) { // MC Check of K(892)0
        if (isTrueK892(trk1, trk2))
          histos.fill(HIST("hReconK892pt"), lResonanceK892.Pt());
      }

      // Mass window cut
      if (std::abs(lResonanceK892.M() - massK892) > cK892masswindow)
        continue;
      // Add one more track loop for K1 reconstruction
      for (auto bTrack : dTracks1) {
        // ID cut
        if (bTrack.index() == trk1.index() || bTrack.index() == trk2.index())
          continue;
        // Track cut
        if (!trackCut(bTrack))
          continue;

        auto bTrkPt = bTrack.pt();
        auto bTrkTPCnSigmaPi = bTrack.tpcNSigmaPi();
        auto isbTrkhasTOF = bTrack.hasTOF();
        auto bTrack_TOFnSigma = (isbTrkhasTOF) ? bTrack.tofNSigmaPi() : -999.;

        // PID selection
        if (!selectionPIDPion(bTrack))
          continue;

        if constexpr (!IsMix) {
          histos.fill(HIST("QA/trkpT_pi_bach"), bTrkPt);
          //  --- PID QA Pion
          histos.fill(HIST("QA/TPC_Nsigma_pi_bach"), bTrkPt, bTrkTPCnSigmaPi);
          if (isbTrkhasTOF) {
            histos.fill(HIST("QA/TOF_Nsigma_pi_bach"), bTrkPt, bTrack_TOFnSigma);
            histos.fill(HIST("QA/TOF_TPC_Map_pi_bach"), bTrack_TOFnSigma, bTrkTPCnSigmaPi);
          }
          histos.fill(HIST("QA/trkDCAxy_pi_bach"), bTrack.dcaXY());
          histos.fill(HIST("QA/trkDCAz_pi_bach"), bTrack.dcaZ());
        }

        // K1 reconstruction
        lDecayDaughter_bach.SetXYZM(bTrack.px(), bTrack.py(), bTrack.pz(), massPi);
        lResonanceK1 = lResonanceK892 + lDecayDaughter_bach;

        // Rapidity cut
        if (lResonanceK1.Rapidity() > cK1MaxRap || lResonanceK1.Rapidity() < cK1MinRap)
          continue;
        // Opening angle cut
        auto lK1Angle = lResonanceK892.Angle(lDecayDaughter_bach.Vect());
        // Pair asymmetry cut
        auto lPairAsym = (lResonanceK892.E() - lDecayDaughter_bach.E()) / (lResonanceK892.E() + lDecayDaughter_bach.E());
        // PiPi, PiKa mass range cut
        TLorentzVector tempPiPi = lDecayDaughter1 + lDecayDaughter_bach;
        TLorentzVector tempPiKa = lDecayDaughter2 + lDecayDaughter_bach;
        if constexpr (!IsMix) {
          histos.fill(HIST("QA/InvMass_piK_pipi"), lResonanceK892.M(), tempPiPi.M());
          histos.fill(HIST("QA/InvMass_piK_pika"), lResonanceK892.M(), tempPiKa.M());
          histos.fill(HIST("QA/K1OA"), lK1Angle);
          histos.fill(HIST("QA/K1PairAsymm"), lPairAsym);
        }
        if (tempPiPi.M() < cPiPiMin || tempPiPi.M() > cPiPiMax)
          continue;
        if (tempPiKa.M() < cPiKaMin || tempPiKa.M() > cPiKaMax)
          continue;
        if (lK1Angle < cMinAngle || lK1Angle > cMaxAngle)
          continue;
        if (lPairAsym < cMinPairAsym || lPairAsym > cMaxPairAsym)
          continue;

        if constexpr (!IsMix) {                                   // Same event pair
          if (trk1.sign() * trk2.sign() < 0) {                    // K892
            if (bTrack.sign() > 0) {                              // bachelor pi+
              if (trk2.sign() > 0) {                              // kaon + means K(892)0 is matter.
                histos.fill(HIST("k1invmass"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_NPP"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("k1invmass_LS"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_PNP"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            } else {                                                 // bachelor pi-
              if (trk2.sign() > 0) {                                 // kaon + means K(892)0 is matter.
                histos.fill(HIST("k1invmass_LS"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_NPN"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("k1invmass"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_PNN"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            }
          } else {                   // K892-LS (false)
            if (bTrack.sign() > 0) { // bachelor pi+
              if (trk2.sign() > 0) { // Kaon+
                histos.fill(HIST("hK1invmass_PPP"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("hK1invmass_PPN"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            } else {                 // bachelor pi-
              if (trk2.sign() > 0) { // Kaon_
                histos.fill(HIST("hK1invmass_NNN"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("hK1invmass_NNP"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            }
          }

          if constexpr (IsMC) {
            if (isTrueK1(trk1, trk2, bTrack)) {
              histos.fill(HIST("hReconK1pt"), lResonanceK1.Pt());
              histos.fill(HIST("QAMC/InvMass_piK_pipi"), lResonanceK892.M(), tempPiPi.M());
              histos.fill(HIST("QAMC/InvMass_piK_pika"), lResonanceK892.M(), tempPiKa.M());
              histos.fill(HIST("QAMC/K1OA"), lK1Angle);
              histos.fill(HIST("QAMC/K1PairAsymm"), lPairAsym);

              if ((bTrack.sign() > 0) && (trk2.sign() > 0)) { // Matter
                histos.fill(HIST("hK1invmass_NPP_MC"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
                histos.fill(HIST("k1invmass_MC"), lResonanceK1.M()); // quick check
              }
              if ((bTrack.sign() < 0) && (trk2.sign() < 0)) { // Anti-matter
                histos.fill(HIST("hK1invmass_PNN_MC"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
                histos.fill(HIST("k1invmass_MC"), lResonanceK1.M()); // quick check
              }
              histos.fill(HIST("hTrueK1pt"), lResonanceK1.Pt());
            } else {
              if (((bTrack.sign() > 0) && (trk2.sign() > 0)) || ((bTrack.sign() < 0) && (trk2.sign() < 0)))
                histos.fill(HIST("k1invmass_noK1"), lResonanceK1.M()); // quick check
            }
          }
        } else {                                                      // Mixed event pair
          if (trk1.sign() * trk2.sign() < 0) {                        // K892
            if (bTrack.sign() > 0) {                                  // bachelor pi+
              if (trk2.sign() > 0) {                                  // kaon + means K(892)0 is matter.
                histos.fill(HIST("k1invmass_Mix"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_NPP_Mix"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("hK1invmass_PNP_Mix"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            } else {                 // bachelor pi-
              if (trk2.sign() > 0) { // kaon + means K(892)0 is matter.
                histos.fill(HIST("hK1invmass_NPN_Mix"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              } else {
                histos.fill(HIST("k1invmass_Mix"), lResonanceK1.M()); // quick check
                histos.fill(HIST("hK1invmass_PNN_Mix"), collision.multV0M(), lResonanceK1.Pt(), lResonanceK1.M());
              }
            }
          }
        }
      }
    }
  }

  void processData(aod::ResoCollision& collision,
                   aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(k1analysis, processData, "Process Event for data without Partitioning", true);

  void processMC(aod::ResoCollision& collision,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(k1analysis, processMC, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents& resoParents)
  {
    for (auto& part : resoParents) {    // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 10323) // K892(0)
        continue;
      if (abs(part.y()) > 0.5) { // rapidity cut
        continue;
      }
      bool pass1 = false;
      bool pass2 = false;
      if (abs(part.daughterPDG1()) == 313 || abs(part.daughterPDG2()) == 313) { // At least one decay to Kaon
        pass2 = true;
      }
      if (abs(part.daughterPDG1()) == kPiPlus || abs(part.daughterPDG2()) == kPiPlus) { // At least one decay to Pion
        pass1 = true;
      }
      if (!pass1 || !pass2) // If we have both decay products
        continue;
      histos.fill(HIST("hTrueK1pt"), part.pt());
    }
  }
  PROCESS_SWITCH(k1analysis, processMCTrue, "Process Event for MC", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(k1analysis, processME, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<k1analysis>(cfgc, TaskName{"lf-k1analysis"})};
}
