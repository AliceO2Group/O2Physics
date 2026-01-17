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
#include <TLorentzVector.h>

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
  Configurable<bool> cfgFindRT{"cfgFindRT", false, "boolean for RT analysis"};

  //  Track selections
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
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
  Configurable<int> selectType{"SelectType", 0, "PID selection type"};

  //  Axis
  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis pTAxis{"pTAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 95.0, 100.0, 105.0, 110.0}, "Centrality Binning"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> lptBinning = {0, 5.0, 13.0, 20.0, 50.0, 1000.0};

    AxisSpec rtAxis = {3, 0, 3};
    AxisSpec lptAxis = {lptBinning}; //  Minimum leading hadron pT selection

    AxisSpec pidqaAxis = {60, -6, 6};
    AxisSpec pTqaAxis = {200, 0, 20};
    AxisSpec phiqaAxis = {72, 0, o2::constants::math::TwoPI}; //  Azimuthal angle axis

    AxisSpec epAxis = {10, 0, o2::constants::math::PI}; //  Event Plane
    AxisSpec epqaAxis = {200, -o2::constants::math::PI, o2::constants::math::PI};
    AxisSpec epResAxis = {200, -2, 2};

    if (cfgFindRT) {
      histos.add("hInvMass_f0980_US", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSpp", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
      histos.add("hInvMass_f0980_LSmm", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, rtAxis, lptAxis}});
    } else {
      histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
      histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, pTAxis, centAxis, epAxis}});
    }

    histos.add("QA/EPhist", "", {HistType::kTH2F, {centAxis, epqaAxis}});
    histos.add("QA/hEPResAB", "", {HistType::kTH2F, {centAxis, epResAxis}});
    histos.add("QA/hEPResBC", "", {HistType::kTH2F, {centAxis, epResAxis}});
    histos.add("QA/hEPResAC", "", {HistType::kTH2F, {centAxis, epResAxis}});

    histos.add("QA/LTpt", "", {HistType::kTH3F, {pTqaAxis, centAxis, phiqaAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/Nsigma_TPC_TOF", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});

    if (doprocessMCLight) {
      histos.add("MCL/hpT_f0980_GEN", "generated f0 signals", HistType::kTH1F, {pTqaAxis});
      histos.add("MCL/hpT_f0980_REC", "reconstructed f0 signals", HistType::kTH3F, {massAxis, pTqaAxis, centAxis});
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
    }
    histos.fill(HIST("QA/EPhist"), collision.cent(), collision.evtPl());
    histos.fill(HIST("QA/hEPResAB"), collision.cent(), collision.evtPlResAB());
    histos.fill(HIST("QA/hEPResBC"), collision.cent(), collision.evtPlResBC());
    histos.fill(HIST("QA/hEPResAC"), collision.cent(), collision.evtPlResAC());
    histos.fill(HIST("QA/LTpt"), lhpT, collision.cent(), lhphi);

    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>> pion1, pion2, reco;
    for (const auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(dTracks, dTracks))) {

      if (!selTrack(trk1) || !selTrack(trk2))
        continue;
      histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
      histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
      histos.fill(HIST("QA/Nsigma_TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());

      if (!selPion(trk1) || !selPion(trk2))
        continue;

      pion1 = ROOT::Math::PxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
      pion2 = ROOT::Math::PxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
      reco = pion1 + pion2;
      if (reco.Rapidity() > cfgMaxRap || reco.Rapidity() < cfgMinRap)
        continue;

      relphi = TVector2::Phi_0_2pi(reco.Phi() - collision.evtPl());
      if (relphi > o2::constants::math::PI) {
        relphi -= o2::constants::math::PI;
      }

      if (trk1.sign() * trk2.sign() < 0) {
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_US"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        }
        histos.fill(HIST("hInvMass_f0980_US_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
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
          histos.fill(HIST("hInvMass_f0980_LSpp"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        }
        histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        if (cfgFindRT) {
          histos.fill(HIST("hInvMass_f0980_LSmm"), reco.M(), reco.Pt(), collision.cent(), rtIndex(reco.Phi(), lhphi), lhpT);
        }
        histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), reco.M(), reco.Pt(), collision.cent(), relphi);
      }
    }
  }

  void processData(o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoEvtPlCollisions>::iterator const& collision,
                   o2::aod::ResoTracks const& resotracks)
  {
    fillHistograms<false>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processData, "Process Event for data", true);

  void processMCLight(
    o2::soa::Join<o2::aod::ResoCollisions, o2::aod::ResoEvtPlCollisions>::iterator const& collision,
    o2::soa::Join<o2::aod::ResoTracks, o2::aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processMCLight, "Process Event for MC", false);

  void processMCTrue(const o2::aod::ResoMCParents& resoParents)
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
      if ((std::abs(part.daughterPDG1()) == kPiPlus &&
           std::abs(part.daughterPDG2()) == kPiPlus)) {
        pass = true;
      }
      if (!pass) // If we have both decay products
        continue;

      histos.fill(HIST("MCL/hpT_f0980_GEN"), part.pt());
    }
  };
  PROCESS_SWITCH(f0980analysis, processMCTrue, "Process Event for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<f0980analysis>(cfgc, TaskName{"lf-f0980analysis"})};
}
