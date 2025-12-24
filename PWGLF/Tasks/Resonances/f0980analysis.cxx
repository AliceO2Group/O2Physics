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
///
/// \file f0980analysis.cxx
/// \brief f0(980) analysis in pp 13.6 TeV
/// \author Yunseul Bae (ybae@cern.ch), Junlee Kim (jikim1290@gmail.com)
/// \since 01/07/2024

#include "PWGLF/DataModel/LFResonanceTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include "Math/LorentzVector.h"
#include "Math/Vector4D.h"
#include "TVector2.h"

#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;

struct f0980analysis {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  //  Event selections
  Configurable<float> cfgMinpT{"cfgMinpT", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Maximum longitudinal DCA"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};

  //  Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};                   // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType |
                                                                                                             // kTPCNCls | kTPCCrossedRows |
                                                                                                             // kTPCCrossedRowsOverNCls |
                                                                                                             // kTPCChi2NDF | kTPCRefit |
                                                                                                             // kITSNCls | kITSChi2NDF |
                                                                                                             // kITSRefit | kITSHits) |
                                                                                                             // kInAcceptanceTracks (kPtRange |
                                                                                                             // kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"};

  //  PID
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"};
  Configurable<double> cMaxTPCnSigmaPionWoTOF{"cMaxTPCnSigmaPionWoTOF", 2.0, "TPC nSigma cut without TOF for Pion"};
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};
  Configurable<int> selectType{"selectType", 0, "PID selection type"};

  //  Axis
  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis pTAxis{"pTAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 95.0, 100.0, 105.0, 110.0}, "Centrality Binning"};

  Configurable<bool> cfgFindRT{"cfgFindRT", false, "boolean for RT analysis"};
  Configurable<bool> cfgFindEP{"cfgFindEP", false, "boolean for Event plane analysis"};
  Configurable<bool> cfgQAEPLT{"cfgQAEPLT", false, "Fill QA histograms for Event Plane and Leading Track"};
  Configurable<bool> cfgQASelection{"cfgQASelection", true, "Fill QA histograms for Selection"};
  Configurable<bool> cfgQAMCTrue{"cfgQAMCTrue", false, "Fill QA histograms for MC True Selection"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> lptBinning = {0, 5.0, 13.0, 20.0, 50.0, 1000.0};

    AxisSpec rtAxis = {3, 0, 3};
    AxisSpec lptAxis = {lptBinning};                    //  Minimum leading hadron pT selection
    AxisSpec epAxis = {10, 0, o2::constants::math::PI}; //  Event Plane
    AxisSpec epqaAxis = {200, -o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec epResAxis = {200, -2, 2};

    AxisSpec pidqaAxis = {60, -6, 6, "#sigma"};
    AxisSpec pTqaAxis = {200, 0, 20, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec phiqaAxis = {72, 0, o2::constants::math::TwoPI, "#Phi"}; //  Azimuthal angle axis
    AxisSpec etaqaAxis = {150, -2, 2, "#eta"};                        //  Pseudorapidity axis
    AxisSpec rapqaAxis = {60, -1.5, 1.5, "#it{y}"};                   //  Rapidity axis
    AxisSpec dcaxyAxis = {200, -5.0, 5.0, "DCA_{xy} (cm)"};           //  DCAxy axis
    AxisSpec dcazAxis = {200, -5.0, 5.0, "DCA_{z} (cm)"};             //  DCAz axis

    AxisSpec collCutAxis = {4, -0.5, 3.5, "Collision cut index for MC"};

    if (cfgFindRT) {
      histos.add("hInvMass_f0980_US_RT", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSpp_RT", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSmm_RT", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
    } else if (cfgFindEP) {
      histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
    } else {
      histos.add("hInvMass_f0980_US", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
      histos.add("hInvMass_f0980_LSpp", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
      histos.add("hInvMass_f0980_LSmm", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis}});
    }

    if (cfgQAEPLT) {
      // Event Plane QA
      histos.add("QA/EPhist", "", {HistType::kTH2F, {centAxis, epqaAxis}});
      histos.add("QA/hEPResAB", "", {HistType::kTH2F, {centAxis, epResAxis}});
      histos.add("QA/hEPResBC", "", {HistType::kTH2F, {centAxis, epResAxis}});
      histos.add("QA/hEPResAC", "", {HistType::kTH2F, {centAxis, epResAxis}});
      //  Leading track pT QA
      histos.add("QA/LTpt", "", {HistType::kTH3F, {pTqaAxis, centAxis, phiqaAxis}});
    }
    if (cfgQASelection) {
      //  General QA
      histos.add("QA/TrackPt", "", {HistType::kTH1F, {pTqaAxis}});
      histos.add("QA/TrackEta", "", {HistType::kTH1F, {etaqaAxis}});
      histos.add("QA/TrackPhi", "", {HistType::kTH1F, {phiqaAxis}});
      //  Track selection QA
      histos.add("QA/trkDCAxy_BC", "DCA_{xy} for pion tracks (before cuts)", HistType::kTH2F, {pTqaAxis, dcaxyAxis});
      histos.add("QA/trkDCAz_BC", "DCA_{z} for pion tracks (before cuts)", HistType::kTH2F, {pTqaAxis, dcazAxis});
      histos.add("QA/trkDCAxy", "DCA_{xy} for pion tracks (after cuts)", HistType::kTH2F, {pTqaAxis, dcaxyAxis});
      histos.add("QA/trkDCAz", "DCA_{z} for pion tracks (after cuts)", HistType::kTH2F, {pTqaAxis, dcazAxis});
      //  PID QA
      histos.add("QA/Nsigma_TPC_BC", "TPC n#sigma^{#pi} (before PID cuts); p_{T} (GeV/c); n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QA/Nsigma_TOF_BC", "TOF n#sigma^{#pi} (before PID cuts); p_{T} (GeV/c); n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QA/Nsigma_TPC_TOF_BC", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});
      histos.add("QA/Nsigma_TPC", "TPC n#sigma^{#pi} (after PID cuts); p_{T} (GeV/c); n#sigma_{TPC}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QA/Nsigma_TOF", "TOF n#sigma^{#pi} (after PID cuts); p_{T} (GeV/c); n#sigma_{TOF}^{#pi}", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
      histos.add("QA/Nsigma_TPC_TOF", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});
    }

    if (doprocessMCRec) {
      histos.add("MCL/hpT_f0980_REC", "Reconstructed f0 signals", HistType::kTH3F, {massAxis, pTqaAxis, centAxis});
    }
    if (doprocessMCTrue) {
      // histos.add("MCL/hpT_f0980_GEN", "Generated f0 signals", HistType::kTH2F, {pTqaAxis, centAxis});
      histos.add("MCL/hpT_f0980_GEN", "Generated f0 signals; selIdx; p_{T} (GeV/c); Centrality (%)", HistType::kTH3F, {collCutAxis, pTqaAxis, centAxis});
      if (cfgQAMCTrue) {
        histos.add("QAMCTrue/f0_pt_y", "Generated f0 ; #it{p}_{T} (GeV/#it{c}) ; #it{y}", HistType::kTH2F, {pTqaAxis, rapqaAxis});
        histos.add("QAMCTrue/f0_pt_cent", "Generated f0 ; #it{p}_{T} (GeV/#it{c}); Centrality (%)", HistType::kTH2F, {pTqaAxis, centAxis});
      }
    }
    histos.print();
  }

  double massPi = MassPionCharged;

  static constexpr float OneThird = 1.0f / 3.0f;
  static constexpr float PIthird = o2::constants::math::PI * OneThird;
  static constexpr float TWOPIthird = o2::constants::math::TwoPI * OneThird;

  int rtIndex(double pairphi, double lhphi)
  {
    double dphi = std::fabs(TVector2::Phi_mpi_pi(lhphi - pairphi));

    if (dphi < PIthird)
      return 0;
    if (dphi < TWOPIthird && dphi > PIthird)
      return 1;
    if (dphi > TWOPIthird)
      return 2;

    return -1;
  }

  template <typename TrackType>
  bool selTrack(const TrackType track)
  {
    if (std::abs(track.pt()) < cfgMinpT)
      return false;
    if (std::fabs(track.eta()) > cfgMaxEta)
      return false;
    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;

    return true;
  }

  template <typename TrackType>
  bool selPion(const TrackType track)
  {
    switch (selectType) {
      case 0:
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
          return false;
        break;
      case 1:
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion)
          return false;
        break;
      case 2:
        if (track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi() >= nsigmaCutCombinedPion * nsigmaCutCombinedPion)
          return false;
        break;
      case 3:
        if (track.hasTOF()) {
          if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
            return false;
        } else {
          if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPionWoTOF)
            return false;
        }
        break;
    }
    return true;
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    double lhpT = 0.;
    double lhphi = 0.;
    double relphi = 0.;
    if (cfgFindRT) {
      for (const auto& trk : dTracks) {
        if (trk.pt() > lhpT) {
          lhpT = trk.pt();
          lhphi = trk.phi();
        }
      }
      if (cfgQAEPLT)
        histos.fill(HIST("QA/LTpt"), lhpT, collision.cent(), lhphi);
    } else if (cfgFindEP) {
      if (cfgQAEPLT) {
        histos.fill(HIST("QA/EPhist"), collision.cent(), collision.evtPl());
        histos.fill(HIST("QA/hEPResAB"), collision.cent(), collision.evtPlResAB());
        histos.fill(HIST("QA/hEPResBC"), collision.cent(), collision.evtPlResBC());
        histos.fill(HIST("QA/hEPResAC"), collision.cent(), collision.evtPlResAC());
      }
    }

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> pion1, pion2, reco;
    for (const auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(dTracks, dTracks))) {
      if (cfgQASelection) {
        histos.fill(HIST("QA/trkDCAxy_BC"), trk1.pt(), trk1.dcaXY());
        histos.fill(HIST("QA/trkDCAz_BC"), trk1.pt(), trk1.dcaZ());
      }
      if (!selTrack(trk1) || !selTrack(trk2))
        continue;
      if (cfgQASelection) {
        histos.fill(HIST("QA/trkDCAxy"), trk1.pt(), trk1.dcaXY());
        histos.fill(HIST("QA/trkDCAz"), trk1.pt(), trk1.dcaZ());

        histos.fill(HIST("QA/Nsigma_TPC_BC"), trk1.pt(), trk1.tpcNSigmaPi());
        if (trk1.hasTOF()) {
          histos.fill(HIST("QA/Nsigma_TOF_BC"), trk1.pt(), trk1.tofNSigmaPi());
          histos.fill(HIST("QA/Nsigma_TPC_TOF_BC"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        }
      }
      if (!selPion(trk1) || !selPion(trk2))
        continue;
      if (cfgQASelection) {
        histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
        if (trk1.hasTOF()) {
          histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
          histos.fill(HIST("QA/Nsigma_TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        }

        histos.fill(HIST("QA/TrackPt"), trk1.pt());
        histos.fill(HIST("QA/TrackEta"), trk1.eta());
        histos.fill(HIST("QA/TrackPhi"), trk1.phi());
      }

      pion1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      pion2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      reco = pion1 + pion2;
      if (reco.Rapidity() > cfgMaxRap || reco.Rapidity() < cfgMinRap)
        continue;

      if (cfgFindEP) {
        relphi = TVector2::Phi_0_2pi(reco.Phi() - collision.evtPl());
        if (relphi > o2::constants::math::PI) {
          relphi -= o2::constants::math::PI;
        }
      }

      if (trk1.sign() * trk2.sign() < 0) {
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_US_RT"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_US_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_US"), reco.M(), reco.Pt(), collision.cent());
        }

        if constexpr (IsMC) {
          if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kPiPlus)
            continue;
          if (trk1.motherId() != trk2.motherId())
            continue;
          if (std::abs(trk1.motherPDG()) != 9010221)
            continue;
          histos.fill(HIST("MCL/hpT_f0980_REC"), reco.M(), reco.Pt(), collision.cent());
        }
      } else if (trk1.sign() > 0 && trk2.sign() > 0) {
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_LSpp_RT"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_LSpp"), reco.M(), reco.Pt(), collision.cent());
        }
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_LSmm_RT"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        } else if (cfgFindEP) {
          histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
        } else {
          histos.fill(HIST("hInvMass_f0980_LSmm"), reco.M(), reco.Pt(), collision.cent());
        }
      }
    }
  }

  void processData(o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoEvtPlCollisions>::iterator const& collision,
                   o2::aod::ResoTracks const& resotracks)
  {
    fillHistograms<false>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processData, "Process Event for data", true);

  void processMCRec(o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoEvtPlCollisions>::iterator const& collision,
                    o2::soa::Join<o2::aod::ResoTracks, o2::aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processMCRec, "Process Event for MC", false);

  void processMCTrue(o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoMCCollisions>::iterator const& resoCollision,
                     o2::aod::ResoMCParents const& resoParents)
  {
    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != 9010221)
        continue;
      if (!part.producedByGenerator())
        continue;
      if (part.y() < cfgMinRap || part.y() > cfgMaxRap) {
        continue;
      }
      bool pass = false;
      if ((std::abs(part.daughterPDG1()) == kPiPlus && std::abs(part.daughterPDG2()) == kPiPlus)) {
        pass = true;
      }
      if (!pass) // If we have both decay products
        continue;

      // no event selection
      histos.fill(HIST("MCL/hpT_f0980_GEN"), 0, part.pt(), resoCollision.cent());
      // |zvtx|<10 cm
      if (resoCollision.isVtxIn10()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 1, part.pt(), resoCollision.cent());
      }
      // |zvtx|<10 cm & TVX trigger
      if (resoCollision.isVtxIn10() && resoCollision.isTriggerTVX()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 2, part.pt(), resoCollision.cent());
      }
      if (resoCollision.isInAfterAllCuts()) {
        histos.fill(HIST("MCL/hpT_f0980_GEN"), 3, part.pt(), resoCollision.cent());
      }
      if (cfgQAMCTrue) {
        histos.fill(HIST("QAMCTrue/f0_pt_y"), part.pt(), part.y());
        histos.fill(HIST("QAMCTrue/f0_pt_cent"), part.pt(), resoCollision.cent());
      }
    }
  };
  PROCESS_SWITCH(f0980analysis, processMCTrue, "Process Event for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<f0980analysis>(cfgc, TaskName{"lf-f0980analysis"})};
}
