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

/// \file lambdaCorrelationAnalysis.cxx
/// \brief R2 correlation of Lambda baryons.
/// \author Yash Patley <yash.patley@cern.ch>

#include <TLorentzVector.h>

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

enum PidType {
  kPion = 0,
  kProton
};

enum ParticleType {
  kLambda = 0,
  kAntiLambda,
  kDummy
};

enum ParticlePairType {
  kLambdaAntiLambda = 0,
  kLambdaLambda,
  kAntiLambdaAntiLambda
};

enum RecGenType {
  kRec = 0,
  kGen
};

struct lambdaCorrAnalysis {

  // Collisions
  Configurable<float> cfg_z_vtx{"cfg_z_vtx", 10.0, "z vertex cut"};
  Configurable<bool> cfg_trigger_tvx_sel{"cfg_trigger_tvx_sel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cfg_tf_border{"cfg_tf_border", false, "Timeframe Border Selection"};
  Configurable<bool> cfg_noitsro_border{"cfg_noitsro_border", false, "No ITSRO Border Cut"};
  Configurable<bool> cfg_sel8_sel{"cfg_sel8_sel", true, "Sel8 (T0A + T0C) Selection"};
  Configurable<bool> cfg_itstpc_vtx{"cfg_itstpc_vtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cfg_pileup_reject{"cfg_pileup_reject", false, "Pileup rejection"};
  Configurable<bool> cfg_zvtx_time_diff{"cfg_zvtx_time_diff", false, "z-vtx time diff selection"};

  // Tracks
  Configurable<float> cfg_pt_min{"cfg_pt_min", 0.2, "p_{T} minimum"};
  Configurable<float> cfg_pt_max{"cfg_pt_max", 4.0, "p_{T} minimum"};
  Configurable<float> cfg_eta_cut{"cfg_eta_cut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cfg_min_crossed_rows{"cfg_min_crossed_rows", 70, "min crossed rows"};
  Configurable<double> cfg_tpc_nsigma{"cfg_tpc_nsigma", 3.0, "TPC NSigma Selection Cut"};
  Configurable<double> cfg_tof_nsigma{"cfg_tof_nsigma", 3.0, "TOF NSigma Selection Cut"};

  // V0s
  Configurable<double> cfg_min_dca_V0_daughters{"cfg_min_dca_V0_daughters", 1.0, "min DCA between V0 daughters"};
  Configurable<double> cfg_min_dca_pos_to_PV{"cfg_min_dca_pos_to_PV", 0.1, "Minimum V0 Positive Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_neg_to_PV{"cfg_min_dca_neg_to_PV", 0.1, "Minimum V0 Negative Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_V0_to_PV{"cfg_min_dca_V0_to_PV", 0.6, "Minimum DCA V0 to PV"};
  Configurable<double> cfg_min_V0_radius{"cfg_min_V0_radius", 0.0, "Minimum V0 radius from PV"};
  Configurable<double> cfg_max_V0_radius{"cfg_max_V0_radius", 50.0, "Maximum V0 radius from PV"};
  Configurable<double> cfg_min_ctau{"cfg_min_ctau", 0.0, "Minimum ctau"};
  Configurable<double> cfg_max_ctau{"cfg_max_ctau", 50.0, "Maximum ctau"};
  Configurable<double> cfg_min_V0_cosPA{"cfg_min_V0_cosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cfg_lambda_mass_window{"cfg_lambda_mass_window", 0.01, "Mass Window to select Lambda"};
  Configurable<double> cfg_kshort_rej{"cfg_kshort_rej", 0.005, "Reject K0Short Candidates"};

  // V0s kinmatic acceptance
  Configurable<float> cfg_v0_pt_min{"cfg_v0_pt_min", 0.5, "Minimum V0 pT"};
  Configurable<float> cfg_v0_pt_max{"cfg_v0_pt_max", 2.5, "Minimum V0 pT"};
  Configurable<float> cfg_v0_rap_max{"cfg_v0_rap_max", 0.8, "|rap| cut"};

  // bool eta/rapidity
  Configurable<bool> cfg_do_eta_analysis{"cfg_do_eta_analysis", false, "Eta Analysis"};

  // V0s MC
  Configurable<bool> cfg_is_primary_lambda{"cfg_is_primary_lambda", true, "Primary Lambda"};
  Configurable<bool> cfg_secondary_lambda{"cfg_secondary_lambda", false, "Secondary Lambda"};
  Configurable<bool> cfg_sel_v0mcrec{"cfg_sel_v0mcrec", true, "MC Rec Selection"};

  // lambda mass window
  Configurable<float> cfg_v0mass_min{"cfg_v0mass_min", 1.11, "Minimum V0 Mass"};
  Configurable<float> cfg_v0mass_max{"cfg_v0mass_max", 1.12, "Maximum V0 Mass"};
  Configurable<bool> cfg_shared_dau{"cfg_shared_dau", true, "Remove shared proton and pion"};

  // Global Configurables
  Configurable<int> cfg_nRapBins{"cfg_nRapBins", 16, "N Rapidity Bins"};
  Configurable<float> cfg_Rap_Min{"cfg_Rap_Min", -0.8, "Minimum Rapidity"};
  Configurable<float> cfg_Rap_Max{"cfg_Rap_Max", 0.8, "Maximum Rapidity"};
  Configurable<int> cfg_nPhiBins{"cfg_nPhiBins", 64, "N Phi Bins"};
  Configurable<float> cfg_Phi_Min{"cfg_Phi_Min", 0, "Minimum Phi"};
  Configurable<float> cfg_Phi_Max{"cfg_Phi_Max", 2 * TMath::Pi(), "Maximum Phi"};

  // Initialize global variables
  float nrapbins = 0.;
  float kminrap = 0.;
  float kmaxrap = 0.;
  float nphibins = 0.;
  float kminphi = 0.;
  float kmaxphi = 0.;
  float rapbinwidth = 0.;
  float phibinwidth = 0.;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // global variable
    nrapbins = static_cast<float>(cfg_nRapBins);
    kminrap = static_cast<float>(cfg_Rap_Min);
    kmaxrap = static_cast<float>(cfg_Rap_Max);
    nphibins = static_cast<float>(cfg_nPhiBins);
    kminphi = static_cast<float>(cfg_Phi_Min);
    kmaxphi = static_cast<float>(cfg_Phi_Max);

    rapbinwidth = (kmaxrap - kminrap) / nrapbins;
    phibinwidth = (kmaxphi - kminphi) / nphibins;

    int knrapphibins = static_cast<int>(cfg_nRapBins) * static_cast<int>(cfg_nPhiBins);
    float kminrapphi = 0.;
    float kmaxrapphi = knrapphibins;

    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");

    const AxisSpec axisV0Mass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(200, 0., 5., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(16, -0.8, 0.8, "rap");
    const AxisSpec axisV0Phi(36, 0., 2. * TMath::Pi(), "#phi (rad)");

    const AxisSpec axisRadius(200, 0, 100, "r(cm)");
    const AxisSpec axisCosPA(120, 0.97, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(200, 0, 2., "dca (cm)");
    const AxisSpec axisDcaProngPV(200, 0, 20., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (cm^{2})");
    const AxisSpec axisCTau(200, 0, 100, "c#tau (cm/#it{c})");
    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    const AxisSpec axisMomPID(80, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    const AxisSpec axisPt(40, 0.5, 2.5, "p_{T} (GeV/#it{c})");
    const AxisSpec axisRap(cfg_nRapBins, cfg_Rap_Min, cfg_Rap_Max, "rap");
    const AxisSpec axisPhi(cfg_nPhiBins, cfg_Phi_Min, cfg_Phi_Max, "#phi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "rap #phi");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1d_collision_posz", "V_{Z} Distribution", kTH1F, {axisVz});
    histos.add("Event/Reco/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h1d_lambda_multiplicity", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_antilambda_multiplicity", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});

    // QA
    histos.add("QA_Checks/h2d_before_topo_cuts_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_after_topo_cuts_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_before_masswincut_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_after_masswincut_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});

    // QA Lambda
    histos.add("QA_Sel_Lambda/h1d_pos_prong_pt", "Pos-Prong p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_pt", "Neg-Prong p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_pos_prong_eta", "Pos-Prong #eta-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_eta", "Neg-Prong #eta-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_pos_prong_phi", "Pos-Prong #phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_phi", "Neg-Prong #phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA_Sel_Lambda/h1d_V0_inv_mass", "V_{0} mass", kTH1F, {axisV0Mass});
    histos.add("QA_Sel_Lambda/h1d_V0_pt", "V_{0} p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_V0_eta", "#eta-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_V0_rap", "y-distribution", kTH1F, {axisV0Rap});
    histos.add("QA_Sel_Lambda/h1d_V0_phi", "#phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA_Sel_Lambda/h1d_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA_Sel_Lambda/h1d_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA_Sel_Lambda/h1d_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA_Sel_Lambda/h1d_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA_Sel_Lambda/h1d_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});

    histos.add("QA_Sel_Lambda/h2d_pos_prong_dEdx_vs_p", "TPC Signal Pos-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_dEdx_vs_p", "TPC Signal Neg-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_nsigma_pr_tpc", "TPC n#sigma Pos-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_nsigma_pi_tpc", "TPC n#sigma Pos-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_nsigma_pr_tpc", "TPC n#sigma Neg-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_nsigma_pi_tpc", "TPC n#sigma Neg-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_nsigma_pr_tof", "TOF n#sigma Pos-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_nsigma_pi_tof", "TOF n#sigma Pos-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_nsigma_pr_tof", "TOF n#sigma Neg-Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_nsigma_pi_tof", "TOF n#sigma Neg-Prong", kTH2F, {axisMomPID, axisNsigma});

    histos.add("QA_Sel_Lambda/h2d_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});

    // QA Anti-Lambda
    histos.addClone("QA_Sel_Lambda/", "QA_Sel_AntiLambda/");

    // single and two particle densities
    histos.add("Reco/h2d_n2_pt1pt2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_pt1pt2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_pt1pt2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_eta1eta2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_eta1eta2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_eta1eta2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_phi1phi2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_phi1phi2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_phi1phi2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_pt1eta2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n2_pt1eta2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n2_pt1eta2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n2_pt1phi2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisPhi});
    histos.add("Reco/h2d_n2_pt1phi2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisPhi});
    histos.add("Reco/h2d_n2_pt1phi2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisPhi});

    histos.add("Reco/h2d_n1_pteta_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPt});
    histos.add("Reco/h2d_n1_pteta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPt});
    histos.add("Reco/h2d_n1_ptrap_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPt});
    histos.add("Reco/h2d_n1_ptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPt});
    histos.add("Reco/h2d_n1_ptphi_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisPhi, axisPt});
    histos.add("Reco/h2d_n1_ptphi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisPhi, axisPt});
    histos.add("Reco/h2d_n1_ptmass_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisV0Mass, axisPt});
    histos.add("Reco/h2d_n1_ptmass_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisV0Mass, axisPt});

    histos.add("Reco/h2d_n1_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n1_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n2_LaP_LaM", "#rho_{2}^{#Lambda - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_LaP_LaP", "#rho_{2}^{#Lambda - #Lambda}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_LaM_LaM", "#rho_{2}^{#bar{#Lambda} - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});

    // MC Generated Histograms
    if (doprocessMCGen) {
      // clone reco histograms
      histos.addClone("Event/Reco/", "Event/McGen/");
      histos.addClone("Reco/", "McGen/");

      // initialize MCGen histograms
      histos.add("McGen/h1d_pt_lambda", "#Lambda p_{T}", kTH1F, {axisV0Pt});
      histos.add("McGen/h1d_pt_antilambda", "#bar{#Lambda} p_{T}", kTH1F, {axisV0Pt});
      histos.add("McGen/h1d_eta_lambda", "#Lambda #eta-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_eta_antilambda", "#bar{#Lambda} #eta-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_y_lambda", "#Lambda y-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_y_antilambda", "#bar{#Lambda} y-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_phi_lambda", "#Lambda #phi-distribution", kTH1F, {axisV0Phi});
      histos.add("McGen/h1d_phi_antilambda", "#bar{#Lambda} #phi-distribution", kTH1F, {axisV0Phi});
    }
  }

  template <typename C>
  bool selCol(C const& col)
  {

    if (fabs(col.posZ()) > cfg_z_vtx) {
      return false;
    }

    if (!col.sel8() && cfg_sel8_sel) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsTriggerTVX) && cfg_trigger_tvx_sel) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) && cfg_tf_border) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoITSROFrameBorder) && cfg_noitsro_border) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsVertexITSTPC) && cfg_itstpc_vtx) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoSameBunchPileup) && cfg_pileup_reject) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && cfg_zvtx_time_diff) {
      return false;
    }

    return true;
  }

  template <PidType pid, typename T>
  bool selPIDTrack(T const& track)
  {

    bool selTPCv0type = false, selTOFv0type = false;
    float tpcNSigma = 0., tofNSigma = 0.;

    switch (pid) {

      case kPion:
        tpcNSigma = track.tpcNSigmaPi();
        tofNSigma = track.tofNSigmaPi();
        break;

      case kProton:
        tpcNSigma = track.tpcNSigmaPr();
        tofNSigma = track.tofNSigmaPr();
        break;
    }

    if (track.hasTOF()) {
      if (fabs(tofNSigma) < cfg_tof_nsigma) {
        selTOFv0type = true;
      }
      if (fabs(tpcNSigma) < cfg_tpc_nsigma) {
        selTPCv0type = true;
      }
    } else {
      selTOFv0type = true;
      if (fabs(tpcNSigma) < cfg_tpc_nsigma) {
        selTPCv0type = true;
      }
    }

    if (selTPCv0type && selTOFv0type) {
      return true;
    }

    return false;
  }

  template <typename V, typename T>
  ParticleType selLambda(V const& v0, T const&)
  {

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (selPIDTrack<kProton>(postrack) && selPIDTrack<kPion>(negtrack)) {
      return kLambda;
    } else if (selPIDTrack<kProton>(negtrack) && selPIDTrack<kPion>(postrack)) {
      return kAntiLambda;
    } else {
      return kDummy;
    }
  }

  template <typename T>
  bool selTrack(T const& track)
  {

    if (track.pt() < cfg_pt_min || track.pt() > cfg_pt_max) {
      return false;
    }

    if (fabs(track.eta()) > cfg_eta_cut) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < cfg_min_crossed_rows) {
      return false;
    }

    return true;
  }

  template <typename V, typename T>
  bool selTopoCuts(V const& v0, T const&)
  {
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // apply track selection on V0 daughters
    if (!selTrack(postrack) || !selTrack(negtrack)) {
      return false;
    }

    // apply topological cuts on V0s
    if (v0.dcaV0daughters() > cfg_min_dca_V0_daughters) {
      return false;
    }

    if (fabs(v0.dcapostopv()) < cfg_min_dca_pos_to_PV) {
      return false;
    }

    if (fabs(v0.dcanegtopv()) < cfg_min_dca_neg_to_PV) {
      return false;
    }

    if (v0.dcav0topv() > cfg_min_dca_V0_to_PV) {
      return false;
    }

    if ((v0.v0radius() > cfg_max_V0_radius) || (v0.v0radius() < cfg_min_V0_radius)) {
      return false;
    }

    if (v0.v0cosPA() < cfg_min_V0_cosPA) {
      return false;
    }

    return true;
  }

  template <typename V>
  bool lambdaKinCuts(V const& v0, float mass)
  {
    // apply Mass Window Selection || Armentros-Podolansky Selection
    if ((mass <= cfg_v0mass_min || mass >= cfg_v0mass_max) || (fabs(v0.mK0Short() - MassK0Short) <= cfg_kshort_rej)) {
      return false;
    }

    // apply kinematic cuts
    if (v0.pt() <= cfg_v0_pt_min || v0.pt() >= cfg_v0_pt_max) {
      return false;
    }

    float rap = 0.;
    if (cfg_do_eta_analysis) {
      rap = v0.eta();
    } else {
      rap = v0.yLambda();
    }

    if (fabs(rap) >= cfg_v0_rap_max) {
      return false;
    }

    return true;
  }

  template <ParticleType v0part, typename V>
  bool selV0MCParticle(V const& v0track)
  {
    auto v0mcpart = v0track.mcParticle();

    if (cfg_is_primary_lambda && !v0mcpart.isPhysicalPrimary()) {
      return false;
    } else if (cfg_secondary_lambda && v0mcpart.isPhysicalPrimary()) {
      return false;
    }

    if (v0part == kLambda && v0mcpart.pdgCode() != 3122) {
      return false;
    } else if (v0part == kAntiLambda && v0mcpart.pdgCode() != -3122) {
      return false;
    }

    return true;
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillQALambda(C const& col, V const& v0, T const&)
  {

    static constexpr std::string_view sub_dir[] = {"QA_Sel_Lambda/", "QA_Sel_AntiLambda/"};

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;

    // daugthers
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    float mass;

    if constexpr (part == kLambda) {
      mass = v0.mLambda();
    } else {
      mass = v0.mAntiLambda();
    }

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_inv_mass"), mass);
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_pt"), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_eta"), v0.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_rap"), v0.yLambda());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_phi"), v0.phi());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pt_vs_alpha"), v0.alpha(), v0.qtarm());

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_pos_to_PV"), fabs(v0.dcapostopv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_neg_to_PV"), fabs(v0.dcanegtopv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_to_PV"), fabs(v0.dcav0topv()));
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_radius"), v0.v0radius());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_ctau"), ctau);

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_pt"), postrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_phi"), postrack.phi());

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_phi"), negtrack.phi());

    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_nsigma_pr_tpc"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_nsigma_pi_tpc"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_nsigma_pr_tpc"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_nsigma_pi_tpc"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());

    if (postrack.hasTOF()) {
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_nsigma_pr_tof"), postrack.tofExpMom(), postrack.tofNSigmaPr());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_nsigma_pi_tof"), postrack.tofExpMom(), postrack.tofNSigmaPi());
    }

    if (negtrack.hasTOF()) {
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_nsigma_pr_tof"), negtrack.tofExpMom(), negtrack.tofNSigmaPr());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_nsigma_pi_tof"), negtrack.tofExpMom(), negtrack.tofNSigmaPi());
    }
  }

  float getPhi(float phi)
  {
    if (phi >= 0) {
      return phi;
    } else {
      return phi + 2 * TMath::Pi();
    }
  }

  template <ParticleType part, RecGenType rec_gen, typename U>
  void fillSingleHistos(U& p)
  {
    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir_hist[] = {"LaP", "LaM"};

    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_pteta_") + HIST(sub_dir_hist[part]), p.Eta(), p.Pt());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_ptrap_") + HIST(sub_dir_hist[part]), p.Rapidity(), p.Pt());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_ptphi_") + HIST(sub_dir_hist[part]), getPhi(p.Phi()), p.Pt());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_ptmass_") + HIST(sub_dir_hist[part]), p.M(), p.Pt());

    if (cfg_do_eta_analysis) {
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_") + HIST(sub_dir_hist[part]), p.Eta(), getPhi(p.Phi()));
    } else {
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_") + HIST(sub_dir_hist[part]), p.Rapidity(), getPhi(p.Phi()));
    }
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, typename U>
  void fillPairHistos(U& p1, U& p2)
  {

    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir_hist[] = {"LaP_LaM", "LaP_LaP", "LaM_LaM"};

    int rapbin1 = 0, rapbin2 = 0;

    if (cfg_do_eta_analysis) {
      rapbin1 = static_cast<int>((p1.Eta() - kminrap) / rapbinwidth);
      rapbin2 = static_cast<int>((p2.Eta() - kminrap) / rapbinwidth);
    } else {
      rapbin1 = static_cast<int>((p1.Rapidity() - kminrap) / rapbinwidth);
      rapbin2 = static_cast<int>((p1.Rapidity() - kminrap) / rapbinwidth);
    }

    int phibin1 = static_cast<int>(getPhi(p1.Phi()) / phibinwidth);
    int phibin2 = static_cast<int>(getPhi(p1.Phi()) / phibinwidth);

    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1pt2_") + HIST(sub_dir_hist[part_pair]), p1.Pt(), p2.Pt());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_eta1eta2_") + HIST(sub_dir_hist[part_pair]), p1.Eta(), p2.Eta());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_phi1phi2_") + HIST(sub_dir_hist[part_pair]), getPhi(p1.Phi()), getPhi(p2.Phi()));
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1eta2_") + HIST(sub_dir_hist[part_pair]), p1.Pt(), p2.Eta());
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1phi2_") + HIST(sub_dir_hist[part_pair]), p1.Pt(), getPhi(p2.Phi()));

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_") + HIST(sub_dir_hist[part_pair]), rapphix + 0.5, rapphiy + 0.5);
    }
  }

  template <RecGenType rg, typename T>
  void getR2CorrHists(T& p1, T& p2)
  {
    int n1 = p1.size(), n2 = p2.size();

    // Lambda-Lambda
    if (n1 != 0) {
      for (int i = 0; i < n1; ++i) {
        fillSingleHistos<kLambda, rg>(p1[i]);
        for (int j = 0; j < n1; ++j) {
          if (i != j) {
            fillPairHistos<kLambdaLambda, rg>(p1[i], p1[j]);
          }
        }
      }
    }

    // AntiLambda-AntiLambda
    if (n2 != 0) {
      for (int i = 0; i < n2; ++i) {
        fillSingleHistos<kAntiLambda, rg>(p2[i]);
        for (int j = 0; j < n2; ++j) {
          if (i != j) {
            fillPairHistos<kAntiLambdaAntiLambda, rg>(p2[i], p2[j]);
          }
        }
      }
    }

    if (n1 == 0 || n2 == 0) {
      return;
    }

    // Lambda-AntiLambda
    for (int i = 0; i < n1; ++i) {
      for (int j = 0; j < n2; ++j) {
        fillPairHistos<kLambdaAntiLambda, rg>(p1[i], p2[j]);
      }
    }
  }

  template <typename L, typename I>
  void updateLambdaList(L& v, I& pos, I& neg)
  {

    std::vector<TLorentzVector> t;
    int n = v.size();

    for (int i = 0; i < n; ++i) {

      for (int j = 0; j < n; ++j) {

        if (i == j) {
          continue;
        }

        if (pos[i] == pos[j] || neg[i] == neg[j]) {
          continue;
        }
      }
    }
  }

  template <bool mcrec, typename C, typename V, typename T>
  void analyzeTracks(C const& collision, V const& V0s, T const& tracks)
  {

    std::vector<TLorentzVector> vLaP, vLaM;
    std::vector<int> vLaPPosId, vLaPNegId, vLaMPosId, vLaMNegId;
    TLorentzVector p;

    for (auto const& v0 : V0s) {

      // sel Lambda or Anti-Lambda
      ParticleType part = selLambda(v0, tracks);

      // apply topological selections
      if (!selTopoCuts(v0, tracks)) {
        continue;
      }

      if (part == kLambda) {
        // apply kinematic cuts
        if (!lambdaKinCuts(v0, v0.mLambda())) {
          continue;
        }
        // McRec selection (PID, PrimaryLambda)
        if constexpr (mcrec == true) {
          if (cfg_sel_v0mcrec && !selV0MCParticle<kLambda>(v0)) {
            continue;
          }
        }
        fillQALambda<kLambda>(collision, v0, tracks);
        p.SetPtEtaPhiM(v0.pt(), v0.eta(), v0.phi(), v0.mLambda());
        vLaP.push_back(p);
        vLaPPosId.push_back(v0.posTrackId());
        vLaPNegId.push_back(v0.negTrackId());
      } else if (part == kAntiLambda) {
        // apply kinematic cuts
        if (!lambdaKinCuts(v0, v0.mAntiLambda())) {
          continue;
        }
        // McRec selection (PID, PrimaryLambda)
        if constexpr (mcrec == true) {
          if (cfg_sel_v0mcrec && !selV0MCParticle<kAntiLambda>(v0)) {
            continue;
          }
        }
        fillQALambda<kAntiLambda>(collision, v0, tracks);
        p.SetPtEtaPhiM(v0.pt(), v0.eta(), v0.phi(), v0.mAntiLambda());
        vLaM.push_back(p);
        vLaMPosId.push_back(v0.posTrackId());
        vLaMNegId.push_back(v0.negTrackId());
      }
    }

    // Update List of Lambdas That Shares Daughters (Remove Fake Lambdas)
    updateLambdaList(vLaP, vLaPPosId, vLaPNegId);
    updateLambdaList(vLaM, vLaMPosId, vLaMNegId);

    // fill density histograms
    getR2CorrHists<kRec>(vLaP, vLaM);

    // clear vectors
    vLaP.clear();
    vLaM.clear();
  }

  // ----------------- Data -------------------- //

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

  void processData(Collisions::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {

    // select collision
    if (!selCol(collision)) {
      return;
    }

    histos.fill(HIST("Event/Reco/h1d_collision_posz"), collision.posZ());

    analyzeTracks<false>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(lambdaCorrAnalysis, processData, "Process for DATA", true);

  // ----------------- McRec -------------------- //

  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  SliceCache cache1;

  // Service<o2::framework::O2DatabasePDG> pdgDB;

  using CollisionsWithMcLabels = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels, aod::PVMults>;
  // using McCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  using McCollisions = aod::McCollisions;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;

  void processMCReco(CollisionsWithMcLabels const& collisions, McCollisions const&, McV0Tracks const& V0s, aod::McParticles const& /*mcParticles*/, TracksMC const& tracks)
  {

    // Event Loop
    for (auto const& collision : collisions) {

      // check for corresponding MCGen Collision
      if (!collision.has_mcCollision()) {
        return;
      }

      // select collision
      if (!selCol(collision)) {
        return;
      }

      histos.fill(HIST("Event/Reco/h1d_collision_posz"), collision.posZ());

      // auto const& mcCollision = collision.mcCollision_as<aod::McCollisions::iterator>();

      // v0-track loop
      auto v0sThisCollision = V0s.sliceBy(perCol, collision.globalIndex());

      analyzeTracks<true>(collision, v0sThisCollision, tracks);
    }
  }

  PROCESS_SWITCH(lambdaCorrAnalysis, processMCReco, "Process for MC Reconstructed", false);

  // --------------------- McGen ---------------------- //

  void processMCGen(McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {

    // apply collision cuts
    if (fabs(mcCollision.posZ()) > cfg_z_vtx) {
      return;
    }

    histos.fill(HIST("Event/McGen/h1d_collision_posz"), mcCollision.posZ());

    TLorentzVector p;
    std::vector<TLorentzVector> vLaP, vLaM;

    for (auto const& mcpart : mcParticles) {

      // check for Primary Lambdas/AntiLambdas
      if (cfg_is_primary_lambda && !mcpart.isPhysicalPrimary()) {
        continue;
      } else if (cfg_secondary_lambda && mcpart.isPhysicalPrimary()) {
        continue;
      }

      // apply kinematic acceptance
      if (mcpart.pt() < cfg_v0_pt_min || mcpart.pt() > cfg_v0_pt_max) {
        continue;
      }

      float rap = 0.;
      if (cfg_do_eta_analysis) {
        rap = mcpart.eta();
      } else {
        rap = mcpart.y();
      }

      if (fabs(rap) > cfg_v0_rap_max) {
        continue;
      }

      p.SetPxPyPzE(mcpart.px(), mcpart.py(), mcpart.pz(), mcpart.e());

      // find daughter ids
      auto mcpart_daughters = mcpart.daughters_as<aod::McParticles>();

      // Fill histograms
      if (mcpart.pdgCode() == 3122) {
        histos.fill(HIST("McGen/h1d_pt_lambda"), mcpart.pt());
        histos.fill(HIST("McGen/h1d_eta_lambda"), mcpart.eta());
        histos.fill(HIST("McGen/h1d_y_lambda"), mcpart.y());
        histos.fill(HIST("McGen/h1d_phi_lambda"), mcpart.phi());
        vLaP.push_back(p);
      } else if (mcpart.pdgCode() == -3122) {
        histos.fill(HIST("McGen/h1d_pt_antilambda"), mcpart.pt());
        histos.fill(HIST("McGen/h1d_eta_antilambda"), mcpart.eta());
        histos.fill(HIST("McGen/h1d_y_antilambda"), mcpart.y());
        histos.fill(HIST("McGen/h1d_phi_antilambda"), mcpart.phi());
        vLaM.push_back(p);
      }
    }

    // fill density histograms
    getR2CorrHists<kGen>(vLaP, vLaM);

    // clear vectors
    vLaP.clear();
    vLaM.clear();
  }

  PROCESS_SWITCH(lambdaCorrAnalysis, processMCGen, "Process for MC Generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdaCorrAnalysis>(cfgc)};
}
