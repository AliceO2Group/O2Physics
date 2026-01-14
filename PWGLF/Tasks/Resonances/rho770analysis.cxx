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
/// \file rho770analysis.cxx
/// \brief rho(770)0 analysis in pp 13 & 13.6 TeV
/// \author Hyunji Lim (hyunji.lim@cern.ch)
/// \since 14/01/2026

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include "Math/Vector4D.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using LorentzVectorPxPyPzMVector = ROOT::Math::PxPyPzMVector;
using LorentzVectorPxPyPzEVector = ROOT::Math::PxPyPzEVector;

enum class TrackPIDMode : int {
  TPCOrTOF = 0,
  OnlyTPC = 1,
  Combined = 2,
  TPCWithTOFVeto = 3
};

struct rho770analysis {
  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.15, "Maximum transverse DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Maximum longitudinal DCA"};
  Configurable<float> cfgMaxTPC{"cfgMaxTPC", 5.0, "Maximum TPC PID with TOF"};
  Configurable<float> cfgMaxTOF{"cfgMaxTOF", 5.0, "Maximum TOF PID with TPC"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};

  // Track selection
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType |
                                                                                                             // kTPCNCls | kTPCCrossedRows |
                                                                                                             // kTPCCrossedRowsOverNCls |
                                                                                                             // kTPCChi2NDF | kTPCRefit |
                                                                                                             // kITSNCls | kITSChi2NDF |
                                                                                                             // kITSRefit | kITSHits) |
                                                                                                             // kInAcceptanceTracks (kPtRange |
                                                                                                             // kEtaRange)
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 1, "Number of TPC cluster"};
  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"}; // refit is included in global track selection
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
  Configurable<int> cfgTPCRows{"cfgTPCRows", 80, "Minimum Number of TPC Crossed Rows "};

  // PID
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};                          // TOF
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 5.0, "TPC nSigma cut for Pion"};                          // TPC
  Configurable<double> cMaxTPCnSigmaPionnoTOF{"cMaxTPCnSigmaPionnoTOF", 3.0, "TPC nSigma cut for Pion in no TOF case"}; // TPC
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", 3.0, "Combined nSigma cut for Pion"};
  Configurable<int> selectTypeInt{"selectType", static_cast<int>(TrackPIDMode::OnlyTPC), "Pion PID selection mode: 0=TPCOrTOF, 1=OnlyTPC, 2=Combined, 3=TPCWithTOFVeto"};

  // Axis
  ConfigurableAxis massAxis{"massAxis", {400, 0.2, 2.2}, "Invariant mass axis"};
  ConfigurableAxis massK0sAxis{"massK0sAxis", {200, 0.46, 0.54}, "K0s Invariant mass axis"};
  ConfigurableAxis massKstarAxis{"massKstarAxis", {200, 0.7, 1.1}, "Kstar Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0}, "Transverse momentum Binning"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 95.0, 100.0, 105.0, 110.0}, "Centrality  Binning"};

  void init(o2::framework::InitContext&)
  {
    AxisSpec pidqaAxis = {120, -6, 6};
    AxisSpec pTqaAxis = {200, 0, 20};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    histos.add("hInvMass_rho770_US", "unlike invariant mass", {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis}});
    histos.add("hInvMass_rho770_LSpp", "++ invariant mass", {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis}});
    histos.add("hInvMass_rho770_LSmm", "-- invariant mass", {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis}});

    histos.add("hInvMass_K0s_US", "K0s unlike invariant mass", {HistType::kTHnSparseF, {massK0sAxis, ptAxis, centAxis}});
    histos.add("hInvMass_K0s_LSpp", "K0s ++ invariant mass", {HistType::kTHnSparseF, {massK0sAxis, ptAxis, centAxis}});
    histos.add("hInvMass_K0s_LSmm", "K0s -- invariant mass", {HistType::kTHnSparseF, {massK0sAxis, ptAxis, centAxis}});

    histos.add("hInvMass_Kstar_US", "Kstar unlike invariant mass", {HistType::kTHnSparseF, {massKstarAxis, ptAxis, centAxis}});
    histos.add("hInvMass_Kstar_LSpp", "Kstar ++ invariant mass", {HistType::kTHnSparseF, {massKstarAxis, ptAxis, centAxis}});
    histos.add("hInvMass_Kstar_LSmm", "Kstar -- invariant mass", {HistType::kTHnSparseF, {massKstarAxis, ptAxis, centAxis}});

    histos.add("QA/hMultiplicityPercent", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});

    histos.add("QA/Nsigma_TPC_BF", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/Nsigma_TOF_BF", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/TPC_TOF_BF", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {pTqaAxis, pidqaAxis}});
    histos.add("QA/TPC_TOF", "", {HistType::kTH2F, {pidqaAxis, pidqaAxis}});

    if (doprocessMCLight) {
      histos.add("MCL/hpT_rho770_GEN", "generated rho770 signals", HistType::kTHnSparseF, {mcLabelAxis, massAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_rho770_REC", "reconstructed rho770 signals", HistType::kTHnSparseF, {massAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_omega_REC", "reconstructed omega signals", HistType::kTHnSparseF, {massAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_K0s_REC", "reconstructed K0s signals", HistType::kTHnSparseF, {massAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_Kstar_REC", "reconstructed Kstar signals", HistType::kTHnSparseF, {massAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_K0s_pipi_REC", "reconstructed K0s signals in pipi hist", HistType::kTHnSparseF, {massK0sAxis, ptAxis, centAxis});
      histos.add("MCL/hpT_Kstar_Kpi_REC", "reconstructed rho770 signals in KK hist", HistType::kTHnSparseF, {massKstarAxis, ptAxis, centAxis});
    }

    histos.print();
  }

  double massPi = MassPionCharged;
  double massKa = MassKaonCharged;

  template <typename TrackType>
  bool selTrack(const TrackType track)
  {
    if (std::abs(track.pt()) < cfgMinPt)
      return false;
    if (std::abs(track.eta()) > cfgMaxEta)
      return false;
    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;
    if (track.tpcNClsCrossedRows() < cfgTPCRows)
      return false;

    return true;
  }

  template <typename TrackType>
  bool selPion(const TrackType track)
  {
    const auto mode = static_cast<TrackPIDMode>(selectTypeInt.value);

    if (mode == TrackPIDMode::TPCOrTOF) { // TPC or TOF
      if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
        return false;
    }
    if (mode == TrackPIDMode::OnlyTPC) { // only TPC
      if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion)
        return false;
    }
    if (mode == TrackPIDMode::Combined) { // combining
      if (track.tpcNSigmaPi() * track.tpcNSigmaPi() + track.tofNSigmaPi() * track.tofNSigmaPi() >= nsigmaCutCombinedPion * nsigmaCutCombinedPion)
        return false;
    }
    if (mode == TrackPIDMode::TPCWithTOFVeto) { // TPC TOF veto
      if (track.hasTOF()) {
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaPi()) >= cMaxTOFnSigmaPion)
          return false;
      } else if (!track.hasTOF()) {
        if (std::fabs(track.tpcNSigmaPi()) >= cMaxTPCnSigmaPionnoTOF)
          return false;
      }
    }
    return true;
  }

  template <typename TrackType>
  bool selKaon(const TrackType track)
  {
    const auto mode = static_cast<TrackPIDMode>(selectTypeInt.value);

    if (mode == TrackPIDMode::TPCOrTOF) {
      if (std::fabs(track.tpcNSigmaKa()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaKa()) >= cMaxTOFnSigmaPion)
        return false;
    }
    if (mode == TrackPIDMode::OnlyTPC) {
      if (std::fabs(track.tpcNSigmaKa()) >= cMaxTPCnSigmaPion)
        return false;
    }
    if (mode == TrackPIDMode::Combined) {
      if (track.tpcNSigmaKa() * track.tpcNSigmaKa() + track.tofNSigmaKa() * track.tofNSigmaKa() >= nsigmaCutCombinedPion * nsigmaCutCombinedPion)
        return false;
    }
    if (mode == TrackPIDMode::TPCWithTOFVeto) {
      if (track.hasTOF()) {
        if (std::fabs(track.tpcNSigmaKa()) >= cMaxTPCnSigmaPion || std::fabs(track.tofNSigmaKa()) >= cMaxTOFnSigmaPion)
          return false;
      } else if (!track.hasTOF()) {
        if (std::fabs(track.tpcNSigmaKa()) >= cMaxTPCnSigmaPionnoTOF)
          return false;
      }
    }
    return true;
  }

  template <bool IsMC, typename CollisionType, typename CenMult, typename TracksType>
  void fillHistograms(const CollisionType& /*collision*/, const CenMult& multiplicity, const TracksType& dTracks)
  {
    histos.fill(HIST("QA/hMultiplicityPercent"), multiplicity);

    LorentzVectorPxPyPzMVector part1, part2, reco;
    for (const auto& [trk1, trk2] : combinations(CombinationsUpperIndexPolicy(dTracks, dTracks))) {

      if (trk1.index() == trk2.index()) {
        histos.fill(HIST("QA/Nsigma_TPC_BF"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF_BF"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF_BF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());

        if (!selTrack(trk1))
          continue;

        histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        continue;
      }

      if (!selTrack(trk1) || !selTrack(trk2))
        continue;

      if (selPion(trk1) && selPion(trk2)) {

        part1 = LorentzVectorPxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
        part2 = LorentzVectorPxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
        reco = part1 + part2;

        if (reco.Rapidity() > cfgMaxRap || reco.Rapidity() < cfgMinRap)
          continue;

        if (trk1.sign() * trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_rho770_US"), reco.M(), reco.Pt(), multiplicity);
          histos.fill(HIST("hInvMass_K0s_US"), reco.M(), reco.Pt(), multiplicity);

          if constexpr (IsMC) {
            if (trk1.motherId() != trk2.motherId())
              continue;
            if (std::abs(trk1.pdgCode()) == kPiPlus && std::abs(trk2.pdgCode()) == kPiPlus) {
              if (std::abs(trk1.motherPDG()) == kRho770_0) {
                histos.fill(HIST("MCL/hpT_rho770_REC"), reco.M(), reco.Pt(), multiplicity);
              } else if (std::abs(trk1.motherPDG()) == kOmega) {
                histos.fill(HIST("MCL/hpT_omega_REC"), reco.M(), reco.Pt(), multiplicity);
              } else if (std::abs(trk1.motherPDG()) == kK0Short) {
                histos.fill(HIST("MCL/hpT_K0s_REC"), reco.M(), reco.Pt(), multiplicity);
                histos.fill(HIST("MCL/hpT_K0s_pipi_REC"), reco.M(), reco.Pt(), multiplicity);
              }
            } else if ((std::abs(trk1.pdgCode()) == kPiPlus && std::abs(trk2.pdgCode()) == kKPlus) || (std::abs(trk1.pdgCode()) == kKPlus && std::abs(trk2.pdgCode()) == kPiPlus)) {
              if (std::abs(trk1.motherPDG()) == kK0Star892) {
                histos.fill(HIST("MCL/hpT_Kstar_REC"), reco.M(), reco.Pt(), multiplicity);
              }
            }
          }
        } else if (trk1.sign() > 0 && trk2.sign() > 0) {
          histos.fill(HIST("hInvMass_rho770_LSpp"), reco.M(), reco.Pt(), multiplicity);
          histos.fill(HIST("hInvMass_K0s_LSpp"), reco.M(), reco.Pt(), multiplicity);
        } else if (trk1.sign() < 0 && trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_rho770_LSmm"), reco.M(), reco.Pt(), multiplicity);
          histos.fill(HIST("hInvMass_K0s_LSmm"), reco.M(), reco.Pt(), multiplicity);
        }
      }

      if ((selPion(trk1) && selKaon(trk2)) || (selKaon(trk1) && selPion(trk2))) {
        if (selPion(trk1) && selKaon(trk2)) {
          part1 = LorentzVectorPxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massPi);
          part2 = LorentzVectorPxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massKa);
        } else if (selKaon(trk1) && selPion(trk2)) {
          part1 = LorentzVectorPxPyPzMVector(trk1.px(), trk1.py(), trk1.pz(), massKa);
          part2 = LorentzVectorPxPyPzMVector(trk2.px(), trk2.py(), trk2.pz(), massPi);
        }
        reco = part1 + part2;

        if (reco.Rapidity() > cfgMaxRap || reco.Rapidity() < cfgMinRap)
          continue;

        if (trk1.sign() * trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_Kstar_US"), reco.M(), reco.Pt(), multiplicity);

          if constexpr (IsMC) {
            if (trk1.motherId() != trk2.motherId())
              continue;
            if ((std::abs(trk1.pdgCode()) == kPiPlus && std::abs(trk2.pdgCode()) == kKPlus) || (std::abs(trk1.pdgCode()) == kKPlus && std::abs(trk2.pdgCode()) == kPiPlus)) {
              if (std::abs(trk1.motherPDG()) == kK0Star892) {
                histos.fill(HIST("MCL/hpT_Kstar_Kpi_REC"), reco.M(), reco.Pt(), multiplicity);
              }
            }
          }
        } else if (trk1.sign() > 0 && trk2.sign() > 0) {
          histos.fill(HIST("hInvMass_Kstar_LSpp"), reco.M(), reco.Pt(), multiplicity);
        } else if (trk1.sign() < 0 && trk2.sign() < 0) {
          histos.fill(HIST("hInvMass_Kstar_LSmm"), reco.M(), reco.Pt(), multiplicity);
        }
      }
    }
  }

  void processData(aod::ResoCollision const& collision, aod::ResoTracks const& resotracks)
  {
    auto multiplicity = collision.cent();
    fillHistograms<false>(collision, multiplicity, resotracks);
  }
  PROCESS_SWITCH(rho770analysis, processData, "Process Event for data", true);

  void processMCLight(ResoMCCols::iterator const& collision,
                      aod::ResoCollisionColls const& collisionIndex,
                      soa::Join<aod::ResoCollisionCandidatesMC, aod::PVMults> const& collisionsMC,
                      soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks,
                      soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&)
  {
    float multiplicity;
    auto linkRow = collisionIndex.iteratorAt(collision.globalIndex());
    auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

    auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

    if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;

    auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
    multiplicity = mcColl.centFT0M();

    if (!collision.isInAfterAllCuts())
      return;

    fillHistograms<true>(collision, multiplicity, resotracks);
  }
  PROCESS_SWITCH(rho770analysis, processMCLight, "Process Event for MC", false);

  void processMCTrue(ResoMCCols::iterator const& collision,
                     aod::ResoCollisionColls const& collisionIndex,
                     aod::ResoMCParents const& resoParents,
                     aod::ResoCollisionCandidatesMC const& collisionsMC,
                     soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&)
  {
    float multiplicity;
    auto linkRow = collisionIndex.iteratorAt(collision.globalIndex());
    auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

    auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

    if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;

    auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
    multiplicity = mcColl.centFT0M();

    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != kRho770_0)
        continue;
      if (!part.producedByGenerator())
        continue;
      if (part.y() < cfgMinRap || part.y() > cfgMaxRap)
        continue;
      if (!(std::abs(part.daughterPDG1()) == kPiPlus && std::abs(part.daughterPDG2()) == kPiPlus))
        continue;

      LorentzVectorPxPyPzEVector truthpar;
      truthpar = LorentzVectorPxPyPzEVector(part.px(), part.py(), part.pz(), part.e());
      auto mass = truthpar.M();

      histos.fill(HIST("MCL/hpT_rho770_GEN"), 0, mass, part.pt(), multiplicity);

      if (collision.isVtxIn10()) {
        histos.fill(HIST("MCL/hpT_rho770_GEN"), 1, mass, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isInSel8()) {
        histos.fill(HIST("MCL/hpT_rho770_GEN"), 2, mass, part.pt(), multiplicity);
      }
      if (collision.isVtxIn10() && collision.isTriggerTVX()) {
        histos.fill(HIST("MCL/hpT_rho770_GEN"), 3, mass, part.pt(), multiplicity);
      }
      if (collision.isInAfterAllCuts()) {
        histos.fill(HIST("MCL/hpT_rho770_GEN"), 4, mass, part.pt(), multiplicity);
      }
    }
  };
  PROCESS_SWITCH(rho770analysis, processMCTrue, "Process Event for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<rho770analysis>(cfgc)};
}
