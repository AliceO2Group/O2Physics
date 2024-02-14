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

/// \author Junlee Kim (jikim1290@gmail.com)

#include <TLorentzVector.h>
#include "TVector2.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct f0980analysis {
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15,
                               "Minimum transverse momentum for charged track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8,
                                "Maximum pseudorapidiy for charged track"};
  Configurable<float> cfgMinDCArToPVcut{"cfgMinDCArToPVcut", -0.5,
                                        "Minimum transverse DCA"};
  Configurable<float> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5,
                                        "Maximum transverse DCA"};
  Configurable<float> cfgMinDCAzToPVcut{"cfgMinDCAzToPVcut", -2.0,
                                        "Minimum longitudinal DCA"};
  Configurable<float> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0,
                                        "Maximum longitudinal DCA"};
  Configurable<float> cfgMaxTPCStandalone{"cfgMaxTPCStandalone", 2.0,
                                          "Maximum TPC PID as standalone"};
  Configurable<float> cfgMaxTPC{"cfgMaxTPC", 5.0, "Maximum TPC PID with TOF"};
  Configurable<float> cfgMaxTOF{"cfgMaxTOF", 3.0, "Maximum TOF PID with TPC"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5, "Minimum rapidity for pair"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5, "Maximum rapidity for pair"};
  Configurable<int> cfgMinTPCncr{"cfgMinTPCncr", 70, "minimum TPC cluster"};

  Configurable<bool> cfgPrimaryTrack{
    "cfgPrimaryTrack", true,
    "Primary track selection"}; // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{
    "cfgGlobalWoDCATrack", true,
    "Global track selection without DCA"}; // kQualityTracks (kTrackType |
                                           // kTPCNCls | kTPCCrossedRows |
                                           // kTPCCrossedRowsOverNCls |
                                           // kTPCChi2NDF | kTPCRefit |
                                           // kITSNCls | kITSChi2NDF |
                                           // kITSRefit | kITSHits) |
                                           // kInAcceptanceTracks (kPtRange |
                                           // kEtaRange)
  Configurable<bool> cfgPVContributor{
    "cfgPVContributor", true,
    "PV contributor track selection"}; // PV Contriuibutor
  Configurable<bool> cfgUseTOF{
    "cfgUseTOF", false,
    "Flag for the usage of TOF for PID"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> ptBinning = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8,
                                     1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5,
                                     5.0, 6.0, 7.0, 8.0, 10.0, 13.0, 20.0};
    std::vector<double> lptBinning = {0, 5.0, 13.0, 20.0, 50.0, 1000.0};

    AxisSpec centAxis = {22, 0, 110};
    AxisSpec ptAxis = {ptBinning};
    AxisSpec massAxis = {400, 0.2, 2.2};
    AxisSpec RTAxis = {3, 0, 3};
    AxisSpec LptAxis = {lptBinning}; // Minimum leading hadron pT selection

    AxisSpec PIDqaAxis = {120, -6, 6};
    AxisSpec pTqaAxis = {200, 0, 20};
    AxisSpec phiqaAxis = {72, 0., 2.0 * constants::math::PI};
    AxisSpec EPAxis = {10, 0, constants::math::PI};
    AxisSpec EPqaAxis = {200, -constants::math::PI, constants::math::PI};
    AxisSpec EPresAxis = {200, -2, 2};

    histos.add("hInvMass_f0980_US", "unlike invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, RTAxis, LptAxis}});
    histos.add("hInvMass_f0980_LSpp", "++ invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, RTAxis, LptAxis}});
    histos.add("hInvMass_f0980_LSmm", "-- invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, RTAxis, LptAxis}});

    histos.add("hInvMass_f0980_US_EPA", "unlike invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, EPAxis}});
    histos.add("hInvMass_f0980_LSpp_EPA", "++ invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, EPAxis}});
    histos.add("hInvMass_f0980_LSmm_EPA", "-- invariant mass",
               {HistType::kTHnSparseF, {massAxis, ptAxis, centAxis, EPAxis}});

    histos.add("QA/hEPResAB", "", {HistType::kTH2F, {centAxis, EPresAxis}});
    histos.add("QA/hEPResAC", "", {HistType::kTH2F, {centAxis, EPresAxis}});
    histos.add("QA/hEPResBC", "", {HistType::kTH2F, {centAxis, EPresAxis}});

    histos.add("QA/Nsigma_TPC", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/Nsigma_TOF", "", {HistType::kTH2F, {pTqaAxis, PIDqaAxis}});
    histos.add("QA/TPC_TOF", "", {HistType::kTH2F, {PIDqaAxis, PIDqaAxis}});

    histos.add("QA/LTpt", "", {HistType::kTH3F, {pTqaAxis, centAxis, phiqaAxis}});
    histos.add("QA/EPhist", "", {HistType::kTH2F, {centAxis, EPqaAxis}});

    if (doprocessMCLight) {
      histos.add("MCL/hpT_f0980_GEN", "generated f0 signals", HistType::kTH1F,
                 {pTqaAxis});
      histos.add("MCL/hpT_f0980_REC", "reconstructed f0 signals",
                 HistType::kTH3F, {massAxis, pTqaAxis, centAxis});
    }

    histos.print();
  }

  double massPi = MassPionCharged;

  int RTIndex(double pairphi, double lhphi)
  {
    double dphi = std::fabs(TVector2::Phi_mpi_pi(lhphi - pairphi));
    if (dphi < constants::math::PI / 3.0)
      return 0;
    if (dphi < 2.0 * constants::math::PI / 3.0 && dphi > constants::math::PI / 3.0)
      return 1;
    if (dphi > 2.0 * constants::math::PI / 3.0)
      return 2;

    return -1;
  }

  template <typename TrackType>
  bool SelTrack(const TrackType track)
  {
    if (track.pt() < cfgMinPt)
      return false;
    if (std::fabs(track.eta()) > cfgMaxEta)
      return false;
    if (track.dcaXY() < cfgMinDCArToPVcut || track.dcaXY() > cfgMaxDCArToPVcut)
      return false;
    if (track.dcaZ() < cfgMinDCAzToPVcut || track.dcaZ() > cfgMaxDCAzToPVcut)
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;

    return true;
  }

  template <typename TrackType>
  bool SelPion(const TrackType track)
  {
    if (track.hasTOF() || !cfgUseTOF) {
      if (std::fabs(track.tpcNSigmaPi()) > cfgMaxTPCStandalone) {
        return false;
      }
    } else {
      if (std::fabs(track.tpcNSigmaPi()) > cfgMaxTPC ||
          std::fabs(track.tofNSigmaPi()) > cfgMaxTOF) {
        return false;
      }
    }
    return true;
  }

  template <bool IsMC, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision,
                      const TracksType& dTracks)
  {
    double LHpt = 0.;
    double LHphi;
    double relPhi;
    for (auto& trk : dTracks) {
      if (trk.pt() > LHpt) {
        LHpt = trk.pt();
        LHphi = trk.phi();
      }
    }
    histos.fill(HIST("QA/EPhist"), collision.cent(), collision.evtPl());
    histos.fill(HIST("QA/hEPResAB"), collision.cent(), collision.evtPlResAB());
    histos.fill(HIST("QA/hEPResAC"), collision.cent(), collision.evtPlResBC());
    histos.fill(HIST("QA/hEPResBC"), collision.cent(), collision.evtPlResAC());
    histos.fill(HIST("QA/LTpt"), LHpt, collision.cent(), LHphi);

    TLorentzVector Pion1, Pion2, Reco;
    for (auto& [trk1, trk2] :
         combinations(CombinationsUpperIndexPolicy(dTracks, dTracks))) {
      if (trk1.index() == trk2.index()) {
        if (!SelTrack(trk1))
          continue;
        histos.fill(HIST("QA/Nsigma_TPC"), trk1.pt(), trk1.tpcNSigmaPi());
        histos.fill(HIST("QA/Nsigma_TOF"), trk1.pt(), trk1.tofNSigmaPi());
        histos.fill(HIST("QA/TPC_TOF"), trk1.tpcNSigmaPi(), trk1.tofNSigmaPi());
        continue;
      }

      if (!SelTrack(trk1) || !SelTrack(trk2))
        continue;
      if (!SelPion(trk1) || !SelPion(trk2))
        continue;

      Pion1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), massPi);
      Pion2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), massPi);
      Reco = Pion1 + Pion2;

      if (Reco.Rapidity() > cfgMaxRap || Reco.Rapidity() < cfgMinRap)
        continue;

      relPhi = TVector2::Phi_0_2pi(Reco.Phi() - collision.evtPl());
      if (relPhi > constants::math::PI) {
        relPhi -= constants::math::PI;
      }

      if (trk1.sign() * trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_US"), Reco.M(), Reco.Pt(),
                    collision.cent(), RTIndex(Reco.Phi(), LHphi), LHpt);
        histos.fill(HIST("hInvMass_f0980_US_EPA"), Reco.M(), Reco.Pt(),
                    collision.cent(), relPhi);
        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != 211 || abs(trk2.pdgCode()) != 211)
            continue;
          if (trk1.motherId() != trk2.motherId())
            continue;
          if (abs(trk1.motherPDG()) != 9010221)
            continue;
          histos.fill(HIST("MCL/hpT_f0980_REC"), Reco.M(), Reco.Pt(),
                      collision.cent());
        }
      } else if (trk1.sign() > 0 && trk2.sign() > 0) {
        histos.fill(HIST("hInvMass_f0980_LSpp"), Reco.M(), Reco.Pt(),
                    collision.cent(), RTIndex(Reco.Phi(), LHphi), LHpt);
        histos.fill(HIST("hInvMass_f0980_LSpp_EPA"), Reco.M(), Reco.Pt(),
                    collision.cent(), relPhi);
      } else if (trk1.sign() < 0 && trk2.sign() < 0) {
        histos.fill(HIST("hInvMass_f0980_LSmm"), Reco.M(), Reco.Pt(),
                    collision.cent(), RTIndex(Reco.Phi(), LHphi), LHpt);
        histos.fill(HIST("hInvMass_f0980_LSmm_EPA"), Reco.M(), Reco.Pt(),
                    collision.cent(), relPhi);
      }
    }
  }

  void processData(aod::ResoCollision& collision,
                   aod::ResoTracks const& resotracks)
  {
    fillHistograms<false>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processData, "Process Event for data", true);

  void processMCLight(
    aod::ResoCollision& collision,
    soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resotracks)
  {
    fillHistograms<true>(collision, resotracks);
  }
  PROCESS_SWITCH(f0980analysis, processMCLight, "Process Event for MC", false);

  void processMCTrue(aod::ResoMCParents& resoParents)
  {

    for (auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 9010221)
        continue;
      if (!part.producedByGenerator())
        continue;
      if (part.y() < cfgMinRap || part.y() > cfgMaxRap) {
        continue;
      }
      bool pass = false;
      if ((abs(part.daughterPDG1()) == 211 &&
           abs(part.daughterPDG2()) == 211)) {
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
