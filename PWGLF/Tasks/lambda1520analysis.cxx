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

/// \file lambda1520analysis.cxx
/// \brief This task reconstructs track-track decay lambda(1520) resonance candidate
/// \author Hirak Kumar Koley <hirak.koley@cern.ch>

#include <TLorentzVector.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::soa;

struct lambda1520analysis {
  // Define slice per Resocollision
  SliceCache cache;
  Preslice<aod::ResoTracks> perResoCollision = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  // Eta-asymmetry switch
  Configurable<bool> isEtaAssym{"isEtaAssym", false, "isEtaAssym"};

  // Pre-selection Track cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> cMinPtcut{"cMinPtcut", 0.2f, "Minimal pT for tracks"};
  Configurable<float> cMaxPtcut{"cMaxPtcut", 10.0f, "Maximal pT for tracks"};
  Configurable<int> cMinTPCncr{"cMinTPCncr", 70, "Minimum number of TPC X rows"};
  // DCA Selections
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.12f, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0f, "Track DCAz cut to PV Minimum"};

  /// PID Selections
  Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 2.0, "TPC nSigma cut for Proton"};              // TPC
  Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 2.0, "TOF nSigma cut for Proton"};              // TOF
  Configurable<double> cMaxTPCnSigmaProtonVETO{"cMaxTPCnSigmaProtonVETO", 3.0, "TPC nSigma VETO cut for Proton"}; // TPC
  Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 2.0, "TPC nSigma cut for Kaon"};                    // TPC
  Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 2.0, "TOF nSigma cut for Kaon"};                    // TOF
  Configurable<double> cMaxTPCnSigmaKaonVETO{"cMaxTPCnSigmaKaonVETO", 3.0, "TPC nSigma VETO cut for Kaon"};       // TPC
  // Kaon
  Configurable<bool> IsptIndependentKaonCut{"IsptIndependentKaonCut", true, "Flag to proceed with pT dependent Kaon pid cuts"};
  Configurable<std::vector<double>> kaonTPCPIDpTintv{"kaonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTPCPIDcuts{"kaonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDpTintv{"kaonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> kaonTOFPIDcuts{"kaonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};
  // Proton
  Configurable<bool> IsptIndependentProtonCut{"IsptIndependentProtonCut", true, "Flag to proceed with pT dependent Proton pid cuts"};
  Configurable<std::vector<double>> protonTPCPIDpTintv{"protonTPCPIDpTintv", {999.}, "pT intervals for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTPCPIDcuts{"protonTPCPIDcuts", {2}, "nSigma list for Kaon TPC PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDpTintv{"protonTOFPIDpTintv", {999.}, "pT intervals for Kaon TOF PID cuts"};
  Configurable<std::vector<double>> protonTOFPIDcuts{"protonTOFPIDcuts", {2}, "nSigma list for Kaon TOF PID cuts"};

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f}, "Mixing bins - multiplicity"};

  /// Figures
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, ""};
  ConfigurableAxis binsMass{"binsMass", {500, 1.3, 3.0}, "Invariant Mass (GeV/#it{c}^2)"};
  ConfigurableAxis binsMult{"binsMult", {150, 0.0, 150.0}, "mult_{FT0M}"};
  ConfigurableAxis binsDCA{"binsDCA", {1000, -5, 5}, ""};
  ConfigurableAxis binsTPCXrows{"binsTPCXrows", {200, 0, 200}, ""};
  ConfigurableAxis binsnSigma{"binsnSigma", {130, -6.5, 6.5}, ""};

  void init(o2::framework::InitContext&)
  {
    // axes
    AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{binsEta, ""};
    AxisSpec axisMassLambda1520{binsMass, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec axisMult{binsMult, "mult_{V0M}"};
    AxisSpec axisDCA{binsDCA, ""};
    AxisSpec axisTPCXrow{binsTPCXrows, ""};
    AxisSpec pidQAAxis = {binsnSigma, ""};

    // Track QA before cuts
    //  --- Track
    histos.add("QA/QAbefore/Track/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QA/QAbefore/Track/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAbefore/Track/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAbefore/Track/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QA/QAbefore/Track/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAbefore/Track/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAbefore/Track/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAbefore/Track/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAbefore/Track/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
    histos.add("QA/QAbefore/Track/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
    histos.add("QA/QAbefore/Track/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});

    // PID QA after cuts
    //  --- Kaon
    histos.add("QA/QAafter/Kaon/TOF_TPC_Map_ka_all", "TOF + TPC Combined PID for Kaon;#sigma_{TOF}^{Kaon};#sigma_{TPC}^{Kaon}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QA/QAafter/Kaon/TOF_Nsigma_ka_all", "TOF NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Kaon};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_all", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly", "TPC NSigma for Kaon;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Kaon};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Kaon/dcaZ", "DCA_{Z} distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm); ", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAafter/Kaon/dcaXY", "DCA_{XY} momentum distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAafter/Kaon/TPC_CR", "# TPC Xrows distribution of selected Kaons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
    histos.add("QA/QAafter/Kaon/pT", "pT distribution of Kaons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
    histos.add("QA/QAafter/Kaon/eta", "#eta distribution of Kaons; #eta; Counts;", {HistType::kTH1F, {axisEta}});

    //  --- Proton
    histos.add("QA/QAafter/Proton/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2D, {pidQAAxis, pidQAAxis}});
    histos.add("QA/QAafter/Proton/TOF_Nsigma_pr_all", "TOF NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH2D, {axisPt, pidQAAxis}});
    histos.add("QA/QAafter/Proton/dcaZ", "DCA_{Z} distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{Z} (cm);", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAafter/Proton/dcaXY", "DCA_{XY} momentum distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); DCA_{XY} (cm);", HistType::kTH2F, {axisPt, axisDCA});
    histos.add("QA/QAafter/Proton/TPC_CR", "# TPC Xrows distribution of selected Protons; #it{p}_{T} (GeV/#it{c}); TPC X rows", HistType::kTH2F, {axisPt, axisTPCXrow});
    histos.add("QA/QAafter/Proton/pT", "pT distribution of Protons; #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
    histos.add("QA/QAafter/Proton/eta", "#eta distribution of Protons; #eta; Counts;", {HistType::kTH1F, {axisEta}});

    if (!doprocessMC || !doprocessMCLight || !doprocessMCTrue) {
      // Mass QA 1D for quick check
      histos.add("Result/Data/lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
      histos.add("Result/Data/lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});   // K+ + Pr
      histos.add("Result/Data/lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}}); // K- + anti-Pr

      // 3d histogram
      histos.add("Result/Data/h3lambda1520invmass", "Invariant mass of #Lambda(1520) K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/Data/h3lambda1520invmassLSPP", "Invariant mass of #Lambda(1520) Like Sign Method K^{#plus}p^{#plus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});   // K+ + Pr
      histos.add("Result/Data/h3lambda1520invmassLSMM", "Invariant mass of #Lambda(1520) Like Sign Method K^{#minus}p^{#minus}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520}); // K- + anti-Pr

      if (doprocessME || doprocessMELight) {
        histos.add("Result/Data/lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}; Invariant Mass (GeV/#it{c}^2); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/h3lambda1520invmassME", "Invariant mass of #Lambda(1520) mixed event K^{#pm}p^{#mp}", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      }

      if (isEtaAssym) {
        histos.add("Result/Data/hlambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/hlambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/hlambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", {HistType::kTH1F, {axisMassLambda1520}});
        histos.add("Result/Data/hlambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", {HistType::kTH1F, {axisMassLambda1520}});

        histos.add("Result/Data/h3lambda1520invmassUnlikeSignAside", "Invariant mass of #Lambda(1520) Unlike Sign A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassLikeSignAside", "Invariant mass of #Lambda(1520) Like Sign A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassUnlikeSignCside", "Invariant mass of #Lambda(1520) Unlike Sign C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        histos.add("Result/Data/h3lambda1520invmassLikeSignCside", "Invariant mass of #Lambda(1520) Like Sign C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        if (doprocessME || doprocessMELight) {
          histos.add("Result/Data/hlambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", {HistType::kTH1F, {axisMassLambda1520}});
          histos.add("Result/Data/hlambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", {HistType::kTH1F, {axisMassLambda1520}});

          histos.add("Result/Data/h3lambda1520invmassMixedAside", "Invariant mass of #Lambda(1520) Mixed A side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
          histos.add("Result/Data/h3lambda1520invmassMixedCside", "Invariant mass of #Lambda(1520) Mixed C side", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
        }
      }
    } else if (doprocessMC || doprocessMCLight || doprocessMCTrue) {
      histos.add("Result/MC/h3recolambda1520invmass", "Invariant mass of Reconstructed MC #Lambda(1520)", HistType::kTH3F, {axisMult, axisPt, axisMassLambda1520});
      histos.add("Result/MC/truelambda1520pt", "pT distribution of True MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("Result/MC/recolambda1520pt", "pT distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisPt}});
      histos.add("Result/MC/recolambda1520invmass", "Inv mass distribution of Reconstructed MC #Lambda(1520); #it{p}_{T} (GeV/#it{c}); Counts;", {HistType::kTH1F, {axisMassLambda1520}});
    }
  }

  double massKa = TDatabasePDG::Instance()->GetParticle(kKMinus)->Mass();
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();

  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    // basic track cuts
    if (track.pt() < cMinPtcut || track.pt() > cMaxPtcut)
      return false;
    if (fabs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (fabs(track.dcaZ()) < cMinDCAzToPVcut || fabs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsCrossedRows() < cMinTPCncr)
      return false;
    if (fabs(track.eta()) > cfgCutEta)
      return false;

    return true;
  }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    auto vKaonTPCPIDpTintv = static_cast<std::vector<double>>(kaonTPCPIDpTintv);
    auto vKaonTPCPIDcuts = static_cast<std::vector<double>>(kaonTPCPIDcuts);
    auto vKaonTOFPIDpTintv = static_cast<std::vector<double>>(kaonTOFPIDpTintv);
    auto vKaonTOFPIDcuts = static_cast<std::vector<double>>(kaonTOFPIDcuts);
    auto lengthOfkaonTPCPIDpTintv = static_cast<int>(vKaonTPCPIDpTintv.size());

    auto vProtonTPCPIDpTintv = static_cast<std::vector<double>>(protonTPCPIDpTintv);
    auto vProtonTPCPIDcuts = static_cast<std::vector<double>>(protonTPCPIDcuts);
    auto vProtonTOFPIDpTintv = static_cast<std::vector<double>>(protonTOFPIDpTintv);
    auto vProtonTOFPIDcuts = static_cast<std::vector<double>>(protonTOFPIDcuts);
    auto lengthOfprotonTPCPIDpTintv = static_cast<int>(vProtonTPCPIDpTintv.size());

    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {
      // Full index policy is needed to consider all possible combinations
      if (trk1.index() == trk2.index())
        continue; // We need to run (0,1), (1,0) pairs as well. but same id pairs are not needed.

      //// Initialize variables
      // Trk1: Proton, Trk2: Kaon
      bool isTrk1Selected{true}, isTrk2Selected{true}, isTrk1hasTOF{false}, isTrk2hasTOF{false};
      auto trk1ptPr = trk1.pt();
      auto trk1NSigmaPrTPC = trk1.tpcNSigmaPr();
      auto trk1NSigmaPrTOF = -999.;
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
      // For Proton candidate:
      if (IsptIndependentProtonCut) {
        // to apply pT-independent PID cuts
        if (isTrk1hasTOF) {
          trk1NSigmaPrTOF = trk1.tofNSigmaPr();
          if (std::abs(trk1NSigmaPrTOF) > cMaxTOFnSigmaProton)
            isTrk1Selected = false;
          if (std::abs(trk1NSigmaPrTPC) > cMaxTPCnSigmaProtonVETO)
            isTrk1Selected = false;
        } else {
          if (std::abs(trk1NSigmaPrTPC) > cMaxTPCnSigmaProton)
            isTrk1Selected = false;
        }
      } else {
        // to apply pT-dependent PID cuts
        if (isTrk1hasTOF) {
          trk1NSigmaPrTOF = trk1.tofNSigmaPr();
          if (lengthOfprotonTPCPIDpTintv > 0) {
            for (int i = 0; i < lengthOfprotonTPCPIDpTintv; i++) {
              if (trk1ptPr < vProtonTOFPIDpTintv[i]) {
                if (std::abs(trk1NSigmaPrTOF) > vProtonTOFPIDcuts[i])
                  isTrk1Selected = false;
                if (std::abs(trk1NSigmaPrTPC) > cMaxTPCnSigmaProtonVETO)
                  isTrk1Selected = false;
              }
            }
          }
        } else {
          if (lengthOfprotonTPCPIDpTintv > 0) {
            for (int i = 0; i < lengthOfprotonTPCPIDpTintv; i++) {
              if (trk1ptPr > vProtonTPCPIDpTintv[i]) {
                if (std::abs(trk1NSigmaPrTPC) > vProtonTPCPIDcuts[i])
                  isTrk1Selected = false;
              }
            }
          }
        }
      }

      // For Kaon candidate:
      // to apply pT-independent PID cuts
      if (IsptIndependentKaonCut) {
        if (isTrk2hasTOF) {
          trk2NSigmaKaTOF = trk2.tofNSigmaKa();
          if (std::abs(trk2NSigmaKaTOF) > cMaxTOFnSigmaKaon)
            isTrk2Selected = false;
          if (std::abs(trk2NSigmaKaTPC) > cMaxTPCnSigmaKaonVETO)
            isTrk2Selected = false;
        } else {
          if (std::abs(trk2NSigmaKaTPC) > cMaxTPCnSigmaKaon)
            isTrk2Selected = false;
        }
      } else {
        // to apply pT-dependent PID cuts
        if (isTrk2hasTOF) {
          trk2NSigmaKaTOF = trk2.tofNSigmaKa();
          if (lengthOfkaonTPCPIDpTintv > 0) {
            for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
              if (trk2ptKa < vKaonTOFPIDpTintv[i]) {
                if (std::abs(trk2NSigmaKaTOF) > vKaonTOFPIDcuts[i])
                  isTrk2Selected = false;
                if (std::abs(trk2NSigmaKaTPC) > cMaxTPCnSigmaKaonVETO)
                  isTrk2Selected = false;
              }
            }
          }
        } else {
          if (lengthOfkaonTPCPIDpTintv > 0) {
            for (int i = 0; i < lengthOfkaonTPCPIDpTintv; i++) {
              if (trk2ptKa > vKaonTPCPIDpTintv[i]) {
                if (std::abs(trk2NSigmaKaTPC) > vKaonTPCPIDcuts[i])
                  isTrk2Selected = false;
              }
            }
          }
        }
      }

      //// QA plots before the selection
      //  --- Track QA all
      histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_pr_all"), trk1ptPr, trk1NSigmaPrTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_pr_all"), trk1ptPr, trk1NSigmaPrTOF);
        histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
      }
      histos.fill(HIST("QA/QAbefore/Track/TPC_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk2hasTOF) {
        histos.fill(HIST("QA/QAbefore/Track/TOF_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QA/QAbefore/Track/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }

      histos.fill(HIST("QA/QAbefore/Track/dcaZ"), trk1ptPr, trk1.dcaZ());
      histos.fill(HIST("QA/QAbefore/Track/dcaXY"), trk1ptPr, trk1.dcaXY());
      histos.fill(HIST("QA/QAbefore/Track/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
      histos.fill(HIST("QA/QAbefore/Track/pT"), trk1ptPr);
      histos.fill(HIST("QA/QAbefore/Track/eta"), trk1.eta());

      // apply the track cut
      if (!trackCut(trk1) || !trackCut(trk2))
        continue;

      //// Apply the pid selection
      if (!isTrk1Selected || !isTrk2Selected)
        continue;

      //// QA plots after the selection
      //  --- PID QA Proton
      histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_all"), trk1ptPr, trk1NSigmaPrTPC);
      if (isTrk1hasTOF) {
        histos.fill(HIST("QA/QAafter/Proton/TOF_Nsigma_pr_all"), trk1ptPr, trk1NSigmaPrTOF);
        histos.fill(HIST("QA/QAafter/Proton/TOF_TPC_Map_pr_all"), trk1NSigmaPrTOF, trk1NSigmaPrTPC);
      }
      if (!isTrk1hasTOF) {
        histos.fill(HIST("QA/QAafter/Proton/TPC_Nsigma_pr_TPConly"), trk1ptPr, trk1NSigmaPrTPC);
      }
      histos.fill(HIST("QA/QAafter/Proton/dcaZ"), trk1ptPr, trk1.dcaZ());
      histos.fill(HIST("QA/QAafter/Proton/dcaXY"), trk1ptPr, trk1.dcaXY());
      histos.fill(HIST("QA/QAafter/Proton/TPC_CR"), trk1ptPr, trk1.tpcNClsCrossedRows());
      histos.fill(HIST("QA/QAafter/Proton/pT"), trk1ptPr);
      histos.fill(HIST("QA/QAafter/Proton/eta"), trk1.eta());

      //  --- PID QA Kaon
      histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTPC);
      if (isTrk2hasTOF) {
        histos.fill(HIST("QA/QAafter/Kaon/TOF_Nsigma_ka_all"), trk2ptKa, trk2NSigmaKaTOF);
        histos.fill(HIST("QA/QAafter/Kaon/TOF_TPC_Map_ka_all"), trk2NSigmaKaTOF, trk2NSigmaKaTPC);
      }
      if (!isTrk2hasTOF) {
        histos.fill(HIST("QA/QAafter/Kaon/TPC_Nsigma_ka_TPConly"), trk2ptKa, trk2NSigmaKaTPC);
      }
      histos.fill(HIST("QA/QAafter/Kaon/dcaZ"), trk2ptKa, trk2.dcaZ());
      histos.fill(HIST("QA/QAafter/Kaon/dcaXY"), trk2ptKa, trk2.dcaXY());
      histos.fill(HIST("QA/QAafter/Kaon/TPC_CR"), trk2ptKa, trk2.tpcNClsCrossedRows());
      histos.fill(HIST("QA/QAafter/Kaon/pT"), trk2ptKa);
      histos.fill(HIST("QA/QAafter/Kaon/eta"), trk2.eta());

      //// Resonance reconstruction
      lDecayDaughter1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPr);
      lDecayDaughter2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massKa);
      lResonance = lDecayDaughter1 + lDecayDaughter2;
      // Rapidity cut
      if (lResonance.Rapidity() > 0.5 || lResonance.Rapidity() < -0.5)
        continue;
      //// Un-like sign pair only
      if (trk1.sign() * trk2.sign() < 0) {
        if constexpr (!IsMix) {
          histos.fill(HIST("Result/Data/lambda1520invmass"), lResonance.M());
          histos.fill(HIST("Result/Data/h3lambda1520invmass"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          if (isEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignAside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          } else if (isEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassUnlikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassUnlikeSignCside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          }
        } else {
          histos.fill(HIST("Result/Data/lambda1520invmassME"), lResonance.M());
          histos.fill(HIST("Result/Data/h3lambda1520invmassME"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          if (isEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedAside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          } else if (isEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassMixedCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassMixedCside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          }
        }

        // MC
        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != kProton || abs(trk2.pdgCode()) != kKPlus)
            continue;
          auto mother1 = trk1.motherId();
          auto mother2 = trk2.motherId();
          if (mother1 == mother2) {              // Same mother
            if (abs(trk1.motherPDG()) == 3124) { // lambda1520(0)
              histos.fill(HIST("Result/MC/reconlambda1520pt"), lResonance.Pt());
              histos.fill(HIST("Result/MC/reconlambda1520invmass"), lResonance.M());
              histos.fill(HIST("Result/MC/h3recolambda1520invmass"), collision.multV0M(), lResonance.Pt(), lResonance.M());
            }
          }
        }
      } else {
        if constexpr (!IsMix) {
          if (isEtaAssym && trk1.eta() > 0.2 && trk1.eta() < 0.8 && trk2.eta() > 0.2 && trk2.eta() < 0.8) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignAside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignAside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          } else if (isEtaAssym && trk1.eta() > -0.6 && trk1.eta() < 0.0 && trk2.eta() > -0.6 && trk2.eta() < 0.0) { // Eta-range will be updated
            histos.fill(HIST("Result/Data/hlambda1520invmassLikeSignCside"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLikeSignCside"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          }
          // Like sign pair ++
          if (trk1.sign() > 0 && trk2.sign() > 0) {
            histos.fill(HIST("Result/Data/lambda1520invmassLSPP"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSPP"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          }

          // Like sign pair --
          if (trk1.sign() < 0 && trk2.sign() < 0) {
            histos.fill(HIST("Result/Data/lambda1520invmassLSMM"), lResonance.M());
            histos.fill(HIST("Result/Data/h3lambda1520invmassLSMM"), collision.multV0M(), lResonance.Pt(), lResonance.M());
          }
        }
      }
    }
  }

  void processData(aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "[DATA] Processing %d collisions", collisions.size());
    for (auto& collision : collisions) {
      Partition<aod::ResoTracks> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) <= static_cast<float_t>(cMaxDCArToPVcut)) /* ((0.0105f + 0.0350f / npow(o2::aod::track::pt, 1.1f)))) */ && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // Basic DCA cuts
      selectedTracks.bindTable(resotracks);
      auto colTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<false, false>(collision, colTracks, colTracks);
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processData, "Process Event for data with partition", true);

  void processDataLight(aod::ResoCollision& collision,
                        aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(lambda1520analysis, processDataLight, "Process Event for data without partition", false);

  void processMC(aod::ResoCollisions& collisions, soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resomctracks, aod::McParticles const& mcParticles)
  {
    LOGF(debug, "[MC] MC events: %d", collisions.size());
    for (auto& collision : collisions) {
      Partition<soa::Join<aod::ResoTracks, aod::ResoMCTracks>> selectedTracks = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) <= static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // Basic DCA cuts
      selectedTracks.bindTable(resomctracks);
      auto colMCTracks = selectedTracks->sliceByCached(aod::resodaughter::resoCollisionId, collision.globalIndex(), cache);
      fillHistograms<true, false>(collision, colMCTracks, colMCTracks);
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
          if (dau.pt() < cMinPtcut || fabs(dau.eta()) > cfgCutEta)
            continue;
          //  LOGF(info,"daughter pt: %f, eta: %f", dau.pt(), dau.eta());
          if (abs(dau.pdgCode()) == kKPlus) { // Decay to Kaons
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kProton) { // Decay to Protons
            pass1 = true;
          }
        }
        if (!pass1 || !pass2)
          continue;
        histos.fill(HIST("Result/MC/truelambda1520pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processMC, "Process Event for MC with partition", false);

  void processMCLight(aod::ResoCollision& collision,
                      soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(lambda1520analysis, processMCLight, "Process Event for MC Light without partition", false);

  void processMCTrue(aod::ResoCollisions& collisions, aod::McParticles const& mcParticles)
  {
    // Not related to the real collisions
    for (auto& part : mcParticles) {             // loop over all MC particles
      if (abs(part.pdgCode()) == 3124) {         // Lambda1520(0)
        if (part.y() > 0.5 || part.y() < -0.5) { // rapidity cut
          continue;
        }
        bool pass1 = false;
        bool pass2 = false;
        for (auto& dau : part.daughters_as<aod::McParticles>()) {
          if (abs(dau.pdgCode()) == kKPlus) { // At least one decay to Kaon
            pass2 = true;
          }
          if (abs(dau.pdgCode()) == kProton) { // At least one decay to Proton
            pass1 = true;
          }
        }
        if (!pass1 || !pass2) // If we have both decay products
          continue;
        histos.fill(HIST("Result/MC/truelambda1520pt"), part.pt());
      }
    }
  }
  PROCESS_SWITCH(lambda1520analysis, processMCTrue, "Process Event for MC only", false);

  // Processing Event Mixing
  using BinningTypeVetZTPCtemp = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::MultV0M>;
  BinningTypeVetZTPCtemp colBinning{{CfgVtxBins, CfgMultBins}, true};
  void processME(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    LOGF(debug, "Event Mixing Started");
    auto tracksTuple = std::make_tuple(resotracks);
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVetZTPCtemp> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      // Kaons
      Partition<aod::ResoTracks> selectedTracks1 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // && aod::tracktag::isKaon == true; // Basic DCA cuts
      selectedTracks1.bindTable(tracks1);
      // Protons
      Partition<aod::ResoTracks> selectedTracks2 = (o2::aod::track::pt > static_cast<float_t>(cMinPtcut)) && (nabs(o2::aod::track::eta) < static_cast<float_t>(cfgCutEta)) && (nabs(o2::aod::track::dcaZ) > static_cast<float_t>(cMinDCAzToPVcut)) && (nabs(o2::aod::track::dcaZ) < static_cast<float_t>(cMaxDCAzToPVcut)) && (nabs(o2::aod::track::dcaXY) < static_cast<float_t>(cMaxDCArToPVcut)) && (aod::resodaughter::tpcNClsCrossedRows > static_cast<uint8_t>(cMinTPCncr)); // && aod::tracktag::isProton == true; // Basic DCA cuts
      selectedTracks2.bindTable(tracks2);

      fillHistograms<false, true>(collision1, selectedTracks1, selectedTracks2);
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processME, "Process EventMixing with partition", true);

  void processMELight(o2::aod::ResoCollisions& collisions, aod::ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    SameKindPair<aod::ResoCollisions, aod::ResoTracks, BinningTypeVetZTPCtemp> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(lambda1520analysis, processMELight, "Process EventMixing light without partition", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambda1520analysis>(cfgc, TaskName{"lf-lambda1520analysis"})};
}
