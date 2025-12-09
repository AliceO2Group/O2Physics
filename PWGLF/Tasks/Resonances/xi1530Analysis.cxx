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

#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TLorentzVector.h>
#include <TRandom.h>

#include <vector>

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

  // Kinematics cuts
  Configurable<bool> cKinCuts{"cKinCuts", true, "Kinematic Cuts for p-K pair opening angle"};
  Configurable<std::vector<float>> cKinCutsPt{"cKinCutsPt", {0.0, 0.4, 0.6, 0.8, 1.0, 1.4, 1.8, 2.2, 2.6, 3.0, 4.0, 5.0, 6.0, 1e10}, "p_{T} of L* for kinematic cuts"};
  Configurable<std::vector<float>> cKinLowerCutsAlpha{"cKinLowerCutsAlpha", {1.5, 1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.08, 0.07, 0.06, 0.04, 0.02}, "Lower cut on Opening angle of p-K of L*"};
  Configurable<std::vector<float>> cKinUpperCutsAlpha{"cKinUpperCutsAlpha", {3.0, 2.0, 1.5, 1.4, 1.0, 0.8, 0.6, 0.5, 0.45, 0.35, 0.3, 0.25, 0.2}, "Upper cut on Opening angle of p-K of L*"};

  // PID Selections
  Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
  Configurable<bool> cUseTpcOnly{"cUseTpcOnly", false, "Use TPC Only selection"};                             // TPC only tracks
  Configurable<float> cRejNsigmaTpc{"cRejNsigmaTpc", 3.0, "Reject tracks to improve purity of TPC PID"};      // Reject missidentified particles when tpc bands merge
  Configurable<float> cRejNsigmaTof{"cRejNsigmaTof", 3.0, "Reject tracks to improve purity of TOF PID"};      // Reject missidentified particles when tpc bands merge
  // Pion
  Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 4.0, "TPC nSigma cut for Pion"};              // TPC
  Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};              // TOF
  Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", 3.0, "Combined nSigma cut for Pion"}; // Combined
  Configurable<std::vector<float>> pionTPCPIDp{"pionTPCPIDp", {0.15, 0.3, 0.35, 0.40, 0.54}, "p dependent TPC cuts Pions"};
  Configurable<std::vector<float>> pionTPCPIDcut{"pionTPCPIDcut", {5., 4., 3., 2.}, "TPC nsigma cuts protons"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisP_pid(600, 0., 6., "p (GeV/c)");
    const AxisSpec axisPt_pid(600, 0., 6., "p_{T} (GeV/c)");
    const AxisSpec axisDCAz(500, -0.5, 0.5, {"DCA_{z} (cm)"});
    const AxisSpec axisDCAxy(240, -0.12, 0.12, {"DCA_{xy} (cm)"});
    const AxisSpec axisTPCNsigma(401, -10.025, 10.025, {"n#sigma^{TPC}"});
    const AxisSpec axisTOFNsigma(401, -10.025, 10.025, {"n#sigma^{TOF}"});
    const AxisSpec axisdEdx(380, 10, 200, {"#frac{dE}{dx}"});
    const AxisSpec axisTpcNsigma(401, -10.025, 10.025, "n#sigma(TPC)");
    const AxisSpec axisTofNsigma(401, -10.025, 10.025, "n#sigma(TOF)");
    const AxisSpec axisXiMass(3000, 0., 3., "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisRadius(1000, 0, 100, "r(cm)");
    const AxisSpec axisCosPA(60, 0.95, 1.01, "cos(#theta_{PA})(rad)");
    const AxisSpec axisDca(200, -10., 10., "dca (cm)");
    const AxisSpec axisDcaDau(100, 0., 10., "Daug DCA (cm^{2})");
    const AxisSpec axisLambdaMass(400, 1.24, 1.64, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisAlpha(200, 0 - 0.15, 2 * TMath::Pi() + 0.15, "#alpha");
    const AxisSpec axisPtXiStar(500, 0, 10, "p_{T} (GeV/#it{c})");
    const AxisSpec axisInvMassXiStar(400, 1.4, 1.8, "m_{#Xi#pi} (GeV/#it{c}^{2})");

    // Create Histograms.
    // Event
    histos.add("Event/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});

    // QA Pions
    histos.add("QA_Bach_Pi/h2d_pi_dca_z", "Pions", kTH2F, {axisPt_pid, axisDCAz});
    histos.add("QA_Bach_Pi/h2d_pi_dca_xy", "Pions", kTH2F, {axisPt_pid, axisDCAxy});
    histos.add("QA_Bach_Pi/h2d_pi_dEdx_p", "TPC Signal Pion", kTH2F, {axisP_pid, axisdEdx});
    histos.add("QA_Bach_Pi/h2d_pi_nsigma_tpc_pt", "Pions", kTH2F, {axisPt_pid, axisTPCNsigma});
    histos.add("QA_Bach_Pi/h2d_pi_nsigma_tpc_p", "Pions", kTH2F, {axisP_pid, axisTPCNsigma});
    histos.add("QA_Bach_Pi/h2d_pi_nsigma_tof_pt", "Pions", kTH2F, {axisPt_pid, axisTOFNsigma});
    histos.add("QA_Bach_Pi/h2d_pi_nsigma_tof_p", "Pions", kTH2F, {axisP_pid, axisTOFNsigma});
    histos.add("QA_Bach_Pi/h2d_pi_nsigma_tof_vs_tpc", "n#sigma(TOF) vs n#sigma(TPC) Pions", kTH2F, {axisTPCNsigma, axisTOFNsigma});

    // QA checks for Pions
    histos.add("QA_Checks/h1d_pi_pt", "p_{T}-spectra Pions", kTH1F, {axisPt_pid});

    // QA Xi
    histos.add("QA_Casc_Xi/h1d_mass_Xi", "#Xi^{+} mass", kTH1F, {axisXiMass});
    histos.add("QA_Casc_Xi/h1d_v0_radius", "V0 Radius", kTH1F, {axisRadius});
    histos.add("QA_Casc_Xi/h1d_casc_radius", "Cascade Radius", kTH1F, {axisRadius});
    histos.add("QA_Casc_Xi/h1d_v0_cosPA", "V0 Cosine of PA", kTH1F, {axisCosPA});
    histos.add("QA_Casc_Xi/h1d_casc_cosPA", "Casc Cosine of PA", kTH1F, {axisCosPA});
    histos.add("QA_Casc_Xi/h1d_dca_postoPV", "DCA Positive to PV", kTH1F, {axisDca});
    histos.add("QA_Casc_Xi/h1d_dca_negtoPV", "DCA Negative to PV", kTH1F, {axisDca});
    histos.add("QA_Casc_Xi/h1d_dca_bachtoPV", "DCA Bachelor to PV", kTH1F, {axisDca});
    histos.add("QA_Casc_Xi/h1d_dca_v0toPV", "DCA V0 to PV", kTH1F, {axisDca});
    histos.add("QA_Casc_Xi/h1d_dca_v0_dau", "DCA V0 Daughter", kTH1F, {axisDcaDau});
    histos.add("QA_Casc_Xi/h1d_dca_casc_dau", "DCA Cascade Daughter", kTH1F, {axisDcaDau});

    // QA Kinematic Cuts
    histos.add("QA_Checks/h2d_xistar_alpha_vs_pt", "#alpha_{oa} vs p_{T} Before Cuts", kTH2F, {axisPtXiStar, axisAlpha});
    histos.add("QA_Checks/h2d_sel_kincuts_xistar_alpha_vs_pt", "#alpha_{oa} vs p_{T} After Cuts", kTH2F, {axisPtXiStar, axisAlpha});

    // Invariant Mass Analysis
    histos.add("Analysis/h1d_mass_Xistar", "Inv Mass Xi(1530)", kTH1D, {axisInvMassXiStar});
    histos.add("Analysis/h3d_mass_vs_pt_Xistar", "Xi(1530)", kTHnSparseD, {axisInvMassXiStar, axisPtXiStar, axisCent});
    histos.add("Analysis/h1d_mass_Xistar_LS", "Inv Mass Xi(1530)", kTH1D, {axisInvMassXiStar});
    histos.add("Analysis/h3d_mass_vs_pt_Xistar_LS", "Xi(1530)", kTHnSparseD, {axisInvMassXiStar, axisPtXiStar, axisCent});
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

  template <typename T>
  bool selectionPIDPions(const T& candidate)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};
    auto tpcPIDp = static_cast<std::vector<float>>(pionTPCPIDp);
    auto tpcPIDcut = static_cast<std::vector<float>>(pionTPCPIDcut);
    int nitr = static_cast<int>(tpcPIDp.size());
    float p = RecoDecay::p(candidate.px(), candidate.py(), candidate.pz());

    float tpcNsigmaPi = std::abs(candidate.tpcNSigmaPi());
    float tpcNsigmaKa = std::abs(candidate.tpcNSigmaKa());
    float tpcNsigmaPr = std::abs(candidate.tpcNSigmaPr());
    float tofNsigmaPi = std::abs(candidate.tofNSigmaPi());
    float tofNsigmaKa = std::abs(candidate.tofNSigmaKa());
    float tofNsigmaPr = std::abs(candidate.tofNSigmaPr());

    float tpcTofNsigmaPi = tpcNsigmaPi * tpcNsigmaPi + tofNsigmaPi * tofNsigmaPi;
    float tpcTofNsigmaKa = tpcNsigmaKa * tpcNsigmaKa + tofNsigmaKa * tofNsigmaKa;
    float tpcTofNsigmaPr = tpcNsigmaPr * tpcNsigmaPr + tofNsigmaPr * tofNsigmaPr;
    float combinedCut = nsigmaCutCombinedPion * nsigmaCutCombinedPion;
    float combinedRejCut = cRejNsigmaTof * cRejNsigmaTpc;

    if (!cUseTpcOnly && candidate.hasTOF()) {
      if (tofNsigmaPi < cMaxTOFnSigmaPion && tofNsigmaPr > cRejNsigmaTof && tofNsigmaKa > cRejNsigmaTof) {
        tofPIDPassed = true;
      }
      // square cut
      if ((nsigmaCutCombinedPion < 0) && (tpcNsigmaPi < cMaxTPCnSigmaPion)) {
        tpcPIDPassed = true;
      }
      // circular cut
      if ((nsigmaCutCombinedPion > 0) && (tpcTofNsigmaPi < combinedCut && tpcTofNsigmaPr > combinedRejCut && tpcTofNsigmaKa > combinedRejCut)) {
        tofPIDPassed = true;
        tpcPIDPassed = true;
      }
    } else {
      tofPIDPassed = true;
      if (cUseTpcOnly) {
        if (tpcNsigmaPi < cMaxTPCnSigmaPion && tpcNsigmaPr > cRejNsigmaTpc && tpcNsigmaKa > cRejNsigmaTpc) {
          tpcPIDPassed = true;
        }
      } else {
        for (int i = 0; i < nitr - 1; ++i) {
          if (p >= tpcPIDp[i] && p < tpcPIDp[i + 1] && (tpcNsigmaPi < tpcPIDcut[i] && tpcNsigmaPr > cRejNsigmaTpc && tpcNsigmaKa > cRejNsigmaTpc)) {
            tpcPIDPassed = true;
          }
        }
        if (tpcPIDPassed && ((tpcNsigmaPi > tpcNsigmaKa) || (tpcNsigmaPi > tpcNsigmaPr))) {
          tpcPIDPassed = false;
        }
      }
    }
    if (tpcPIDPassed && tofPIDPassed) {
      return true;
    }
    return false;
  }

  // kinematic cuts method
  template <typename U, typename K, typename T>
  bool kinCuts(U trk1, K trk2, T p, float& alpha)
  {
    // initialize
    std::vector<float> kinCutsPt = static_cast<std::vector<float>>(cKinCutsPt);
    std::vector<float> kinLowerCutsAlpha = static_cast<std::vector<float>>(cKinLowerCutsAlpha);
    std::vector<float> kinUpperCutsAlpha = static_cast<std::vector<float>>(cKinUpperCutsAlpha);
    int kinCutsSize = static_cast<int>(kinUpperCutsAlpha.size());

    TVector3 v1(trk1.px(), trk1.py(), trk1.pz());
    TVector3 v2(trk2.px(), trk2.py(), trk2.pz());
    alpha = v1.Angle(v2);

    for (int i = 0; i < kinCutsSize; ++i) {
      if ((p.Pt() > kinCutsPt[i] && p.Pt() <= kinCutsPt[i + 1]) && (alpha < kinLowerCutsAlpha[i] || alpha > kinUpperCutsAlpha[i])) {
        return false;
      }
    }

    return true;
  }

  template <typename T>
  void fillPionQAHistos(T const& track)
  {

    float p = RecoDecay::p(track.px(), track.py(), track.pz());

    // select particle (Pion)
    if (selectionPIDPions(track)) {
      histos.fill(HIST("QA_Checks/h1d_pi_pt"), track.pt());
      histos.fill(HIST("QA_Bach_Pi/h2d_pi_dca_z"), track.pt(), track.dcaZ());
      histos.fill(HIST("QA_Bach_Pi/h2d_pi_dca_xy"), track.pt(), track.dcaXY());
      histos.fill(HIST("QA_Bach_Pi/h2d_pi_dEdx_p"), p, track.tpcSignal());
      histos.fill(HIST("QA_Bach_Pi/h2d_pi_nsigma_tpc_p"), p, track.tpcNSigmaPi());
      histos.fill(HIST("QA_Bach_Pi/h2d_pi_nsigma_tpc_pt"), track.pt(), track.tpcNSigmaPi());
      if (!cUseTpcOnly && track.hasTOF()) {
        histos.fill(HIST("QA_Bach_Pi/h2d_pi_nsigma_tof_p"), p, track.tofNSigmaPi());
        histos.fill(HIST("QA_Bach_Pi/h2d_pi_nsigma_tof_pt"), track.pt(), track.tofNSigmaPi());
        histos.fill(HIST("QA_Bach_Pi/h2d_pi_nsigma_tof_vs_tpc"), track.tpcNSigmaPi(), track.tofNSigmaPi());
      }
    }
  }

  template <bool mix, bool mc, typename cascType, typename trackType>
  void fillDataHisto(cascType const& cascTrks, trackType const& bachTrks, float cent_class)
  {
    TLorentzVector p1, p2, p;

    // trk1 -> Xi | trk2 -> pi
    for (auto const& [trk1, trk2] : soa::combinations(soa::CombinationsFullIndexPolicy(cascTrks, bachTrks))) {

      // // do not analyze same index tracks
      // if (trk1.index() == trk2.index()) {
      //   continue;
      // }

      // select primary tracks for Pion selection
      if (!selBachTracks(trk2)) {
        continue;
      }

      // Apply Pion PID
      if (!selectionPIDPions(trk2)) {
        continue;
      }

      // Make Lorentz 4-vector
      p1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), trk1.mXi());
      p2.SetXYZM(trk2.px(), trk2.py(), trk2.pz(), MassPionCharged);
      p = p1 + p2;

      // Get the opening angle b/w Xi and Pion and kinematic cut flag.
      float alpha = 0;
      bool kinCutFlag{true};
      if (cKinCuts) {
        kinCutFlag = kinCuts(trk1, trk2, p, alpha);
        if constexpr (!mix && !mc) {
          histos.fill(HIST("QA_Checks/h2d_xistar_alpha_vs_pt"), p.Pt(), alpha);
        }
      }

      // apply kincuts
      if (cKinCuts && !kinCutFlag) {
        continue;
      }

      if (std::abs(p.Rapidity()) > 0.5)
        continue;

      if constexpr (!mix && !mc) {
        if (trk1.sign() * trk2.sign() < 0) {
          histos.fill(HIST("Analysis/h1d_mass_Xistar"), p.M());
          histos.fill(HIST("Analysis/h3d_mass_vs_pt_Xistar"), p.M(), p.Pt(), cent_class);
        } else {
          histos.fill(HIST("Analysis/h1d_mass_Xistar_LS"), p.M());
          histos.fill(HIST("Analysis/h3d_mass_vs_pt_Xistar_LS"), p.M(), p.Pt(), cent_class);
        }
      }
    }
  }

  void process(aod::ResoCollisions::iterator const& resoCollision, aod::ResoCascades const& cascTracks, aod::ResoTracks const& resoTracks)
  {

    histos.fill(HIST("Event/h1d_ft0m_mult_percentile"), resoCollision.cent());

    // QA for pion
    for (auto const& track : resoTracks) {

      // apply primary selection
      if (!selBachTracks(track)) {
        continue;
      }

      // QA histos
      fillPionQAHistos(track);
    }

    // QA Xi
    for (auto const& casc : cascTracks) {
      histos.fill(HIST("QA_Casc_Xi/h1d_mass_Xi"), casc.mXi());
      histos.fill(HIST("QA_Casc_Xi/h1d_v0_radius"), casc.transRadius());
      histos.fill(HIST("QA_Casc_Xi/h1d_casc_radius"), casc.cascTransRadius());
      histos.fill(HIST("QA_Casc_Xi/h1d_v0_cosPA"), casc.v0CosPA());
      histos.fill(HIST("QA_Casc_Xi/h1d_casc_cosPA"), casc.cascCosPA());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_postoPV"), casc.dcapostopv());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_negtoPV"), casc.dcanegtopv());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_bachtoPV"), casc.dcabachtopv());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_v0toPV"), casc.dcav0topv());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_v0_dau"), casc.daughDCA());
      histos.fill(HIST("QA_Casc_Xi/h1d_dca_casc_dau"), casc.cascDaughDCA());
    }

    fillDataHisto<false, false>(cascTracks, resoTracks, resoCollision.cent());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cascadeXiAnalysis>(cfgc)};
}
