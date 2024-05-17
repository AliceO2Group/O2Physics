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

/// \file xi1530Analysis.cxx
/// \brief Invariant Mass Reconstruction of Xi(1530) Resonance
/// \author Yash Patley <yash.patley@cern.ch>

#include <TLorentzVector.h>
#include <TRandom.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct cascadeXiAnalysis {

  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;

  // Tracks
  Configurable<float> cPtMin{"cPtMin", 0.15, "Minimum Track pT"};
  Configurable<float> cEtaCut{"cEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<float> cDcaz{"cDcaz", 1., "Minimum DCAz"};
  Configurable<float> cDcaxy{"cDcaxy", 0.1, "Minimum DCAxy"};
  Configurable<float> cPIDprecut{"cPIDprecut", 5, "Preselection PID TPC TOF cut"};
  Configurable<bool> cPrimaryTrack{"cPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cGlobalWoDCATrack{"cGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cPVContributor{"cPVContributor", true, "PV contributor track selection"};           // PV Contriuibutor

  // PID Pions
  Configurable<float> cMaxTpcNsigmaPi{"cMaxTpcNsigmaPi", 3.0, "Max TPCNsigma Pions"};
  Configurable<float> cMaxTofNsigmaPi{"cMaxTofNsigmaPi", 3.0, "Max TOFNsigma Pions"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisP(400, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisTpcNsigma(401, -10.025, 10.025, "n#sigma(TPC)");
    const AxisSpec axisTofNsigma(401, -10.025, 10.025, "n#sigma(TOF)");
    const AxisSpec axisXiMass(3000, 0., 3., "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisRadius(1000, 0, 100, "r(cm)");
    const AxisSpec axisCosPA(1000, 0.90, 1.1, "cos(#theta_{PA})(rad)");
    const AxisSpec axisDca(1000, -10., 10., "dca (cm)");
    const AxisSpec axisDcaDau(1000, 0., 10., "Daug DCA (cm^{2})");
    const AxisSpec axisLambdaMass(1000, 0., 10., "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisPtXiStar(500, 0, 10, "p_{T} (GeV/#it{c})");
    const AxisSpec axisInvMassXiStar(400, 1.4, 1.8, "m_{#Xi#pi} (GeV/#it{c}^{2})");

    histos.add("QA/CascXi/h1d_mass_Xi", "#Xi^{+} mass", kTH1F, {axisXiMass});
    histos.add("QA/CascXi/h1d_v0_radius", "V0 Radius", kTH1F, {axisRadius});
    histos.add("QA/CascXi/h1d_casc_radius", "Cascade Radius", kTH1F, {axisRadius});
    histos.add("QA/CascXi/h1d_v0_cosPA", "V0 Cosine of PA", kTH1F, {axisCosPA});
    histos.add("QA/CascXi/h1d_casc_cosPA", "Casc Cosine of PA", kTH1F, {axisCosPA});
    histos.add("QA/CascXi/h1d_dca_postoPV", "DCA Positive to PV", kTH1F, {axisDca});
    histos.add("QA/CascXi/h1d_dca_negtoPV", "DCA Negative to PV", kTH1F, {axisDca});
    histos.add("QA/CascXi/h1d_dca_bachtoPV", "DCA Bachelor to PV", kTH1F, {axisDca});
    histos.add("QA/CascXi/h1d_dca_v0toPV", "DCA V0 to PV", kTH1F, {axisDca});
    histos.add("QA/CascXi/h1d_dca_v0_dau", "DCA V0 Daughter", kTH1F, {axisDcaDau});
    histos.add("QA/CascXi/h1d_dca_casc_dau", "DCA Cascade Daughter", kTH1F, {axisDcaDau});
    histos.add("QA/Pions/h2d_pi_tpc_nsigma_vs_p", "Pions", kTH2D, {axisP, axisTpcNsigma});
    histos.add("QA/Pions/h2d_pi_tof_nsigma_vs_p", "Pions", kTH2D, {axisP, axisTofNsigma});
    histos.add("QA/Pions/h2d_pi_nsigma_tof_vs_tpc", "Pions", kTH2D, {axisTpcNsigma, axisTofNsigma});
    histos.add("Analysis/h1d_mass_Xistar", "Inv Mass Xi(1530)", kTH1D, {axisInvMassXiStar});
    histos.add("Analysis/h2d_mass_vs_pt_Xistar", "Xi(1530)", kTH2D, {axisInvMassXiStar, axisPtXiStar});
    histos.add("Analysis/h1d_mass_Xistar_LS", "Inv Mass Xi(1530)", kTH1D, {axisInvMassXiStar});
    histos.add("Analysis/h2d_mass_vs_pt_Xistar_LS", "Xi(1530)", kTH2D, {axisInvMassXiStar, axisPtXiStar});
  }

  template <typename trackType>
  bool selBachTracks(trackType const& track)
  {

    if (track.pt() < cPtMin)
      return false;

    if (std::abs(track.eta()) > cEtaCut)
      return false;

    if (std::abs(track.dcaZ()) > cDcaz)
      return false;

    if (std::abs(track.dcaXY()) > cDcaxy)
      return false;

    if (cPrimaryTrack && !track.isPrimaryTrack())
      return false;

    if (cGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;

    if (cPVContributor && !track.isPVContributor())
      return false;

    return true;
  }

  template <typename trackType>
  bool selPIDPions(trackType const& track)
  {
    bool tpcPIDselFlag{false}, tofPIDselFlag{false};

    float tpcNsigmaPi = std::abs(track.tpcNSigmaPi());
    float tofNsigmaPi = std::abs(track.tofNSigmaPi());

    if (track.hasTOF()) {
      if (tofNsigmaPi < cMaxTofNsigmaPi) {
        tofPIDselFlag = true;
      }
      if (tpcNsigmaPi < cMaxTpcNsigmaPi) {
        tpcPIDselFlag = true;
      }
    } else {
      tofPIDselFlag = true;
      if (tpcNsigmaPi < cMaxTofNsigmaPi) {
        tpcPIDselFlag = true;
      }
    }
    if (tofPIDselFlag && tpcPIDselFlag) {
      return true;
    }
    return false;
  }

  template <bool mix, bool mc, typename cascType, typename trackType>
  void fillDataHisto(trackType const& bachTrks, cascType const& cascTrks, const float /*cent*/)
  {
    TLorentzVector p1, p2, p;
    float pi_p_tot = 0;
    // trk1 -> pi
    // trk2 -> Xi
    for (auto const& [trk1, trk2] : soa::combinations(soa::CombinationsFullIndexPolicy(cascTrks, bachTrks))) {
      if (trk1.index() == trk2.index()) {
        continue;
      }

      // select primary pions (one of the decay daughter of Xi*)
      if (!selBachTracks(trk2)) {
        continue;
      }

      if (!selPIDPions(trk2)) {
        continue;
      }

      pi_p_tot = TMath::Sqrt(trk2.px() * trk2.px() + trk2.py() * trk2.py() + trk2.pz() * trk2.pz());

      // fill QA for pion selection
      if constexpr (!mix) {
        histos.fill(HIST("QA/Pions/h2d_pi_tpc_nsigma_vs_p"), pi_p_tot, trk2.tpcNSigmaPi());
        if (trk2.hasTOF()) {
          histos.fill(HIST("QA/Pions/h2d_pi_tof_nsigma_vs_p"), pi_p_tot, trk2.tofNSigmaPi());
          histos.fill(HIST("QA/Pions/h2d_pi_nsigma_tof_vs_tpc"), trk2.tpcNSigmaPi(), trk2.tofNSigmaPi());
        }
      }

      // Make Lorentz 4-vector
      p1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), trk1.mXi());
      p2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), MassPionCharged);
      p = p1 + p2;

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      if constexpr (!mix && !mc) {
        if (trk1.sign() * trk2.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_mass_Xistar"), p.M());
          histos.fill(HIST("Analysis/h2d_mass_vs_pt_Xistar"), p.M(), p.Pt());
        } else {
          histos.fill(HIST("Analysis/h1d_mass_Xistar_LS"), p.M());
          histos.fill(HIST("Analysis/h2d_mass_vs_pt_Xistar_LS"), p.M(), p.Pt());
        }
      }
    }
  }

  void process(aod::ResoCollisions::iterator const& resoCollision, aod::ResoCascades const& cascTracks, aod::ResoTracks const& resoTracks)
  {
    fillDataHisto<false, false>(resoTracks, cascTracks, resoCollision.cent());
  }

  void processReso(aod::ResoCollisions::iterator const& /*col*/, aod::ResoCascades const& cascTracks)
  {

    for (auto const& casc : cascTracks) {
      histos.fill(HIST("QA/CascXi/h1d_mass_Xi"), casc.mXi());
      histos.fill(HIST("QA/CascXi/h1d_v0_radius"), casc.transRadius());
      histos.fill(HIST("QA/CascXi/h1d_casc_radius"), casc.casctransRadius());
      histos.fill(HIST("QA/CascXi/h1d_v0_cosPA"), casc.v0CosPA());
      histos.fill(HIST("QA/CascXi/h1d_casc_cosPA"), casc.cascCosPA());
      histos.fill(HIST("QA/CascXi/h1d_dca_postoPV"), casc.dcapostopv());
      histos.fill(HIST("QA/CascXi/h1d_dca_negtoPV"), casc.dcanegtopv());
      histos.fill(HIST("QA/CascXi/h1d_dca_bachtoPV"), casc.dcabachtopv());
      histos.fill(HIST("QA/CascXi/h1d_dca_v0toPV"), casc.dcav0topv());
      histos.fill(HIST("QA/CascXi/h1d_dca_v0_dau"), casc.daughDCA());
      histos.fill(HIST("QA/CascXi/h1d_dca_casc_dau"), casc.cascdaughDCA());
    }
  }

  PROCESS_SWITCH(cascadeXiAnalysis, processReso, "ResoProcess", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cascadeXiAnalysis>(cfgc)};
}
