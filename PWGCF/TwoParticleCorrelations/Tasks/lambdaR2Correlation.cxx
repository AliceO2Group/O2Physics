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

#include <vector>

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
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace lambdacollision
{
DECLARE_SOA_COLUMN(Cent, cent, float);
}
DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ);
using LambdaCollision = LambdaCollisions::iterator;

namespace lambdamcgencollision
{
}
DECLARE_SOA_TABLE(LambdaMCGenCollisions, "AOD", "LMCGENCOLS", o2::soa::Index<>,
                  o2::aod::mccollision::PosX,
                  o2::aod::mccollision::PosY,
                  o2::aod::mccollision::PosZ);
using LambdaMCGenCollision = LambdaMCGenCollisions::iterator;

namespace lambdatrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaCollision, lambdaCollision);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(PosTrackId, postrackid, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negtrackid, int64_t);
DECLARE_SOA_COLUMN(V0Type, v0type, int8_t);
DECLARE_SOA_COLUMN(Cospa, cospa, float);
DECLARE_SOA_COLUMN(DcaDau, dcadau, float);
} // namespace lambdatrack
DECLARE_SOA_TABLE(LambdaTracks, "AOD", "LAMBDATRACKS", o2::soa::Index<>,
                  lambdatrack::LambdaCollisionId,
                  lambdatrack::Px,
                  lambdatrack::Py,
                  lambdatrack::Pz,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::V0Type,
                  lambdatrack::Cospa,
                  lambdatrack::DcaDau);
using LambdaTrack = LambdaTracks::iterator;

namespace lambdamcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMCGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(LambdaMCGenTracks, "AOD", "LMCGENTRACKS", o2::soa::Index<>,
                  lambdamcgentrack::LambdaMCGenCollisionId,
                  o2::aod::mcparticle::Px,
                  o2::aod::mcparticle::Py,
                  o2::aod::mcparticle::Pz,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::V0Type,
                  lambdatrack::Cospa,
                  lambdatrack::DcaDau);
using LambdaMCGenTrack = LambdaMCGenTracks::iterator;

} // namespace o2::aod

enum ParticleType {
  kLambda = 0,
  kAntiLambda
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

struct lambdaCorrTableProducer {

  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::LambdaMCGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMCGenTracks> lambdaMCGenTrackTable;

  // Collisions
  Configurable<float> cfg_min_z_vtx{"cfg_min_z_vtx", -10.0, "z vertex cut"};
  Configurable<float> cfg_max_z_vtx{"cfg_max_z_vtx", 10.0, "z vertex cut"};
  Configurable<bool> cfg_sel8_sel{"cfg_sel8_sel", true, "Sel8 (T0A + T0C) Selection"};
  Configurable<bool> cfg_trigger_tvx_sel{"cfg_trigger_tvx_sel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cfg_tf_border{"cfg_tf_border", false, "Timeframe Border Selection"};
  Configurable<bool> cfg_noitsro_border{"cfg_noitsro_border", false, "No ITSRO Border Cut"};
  Configurable<bool> cfg_itstpc_vtx{"cfg_itstpc_vtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cfg_pileup_reject{"cfg_pileup_reject", false, "Pileup rejection"};
  Configurable<bool> cfg_zvtx_time_diff{"cfg_zvtx_time_diff", false, "z-vtx time diff selection"};

  // Tracks
  Configurable<float> cfg_pt_min{"cfg_pt_min", 0.2, "p_{T} minimum"};
  Configurable<float> cfg_pt_max{"cfg_pt_max", 4.0, "p_{T} minimum"};
  Configurable<float> cfg_eta_cut{"cfg_eta_cut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cfg_min_crossed_rows{"cfg_min_crossed_rows", 70, "min crossed rows"};
  Configurable<double> cfg_tpc_nsigma{"cfg_tpc_nsigma", 2.0, "TPC NSigma Selection Cut"};
  Configurable<double> cfg_track_dcaXY_min{"cfg_track_dcaXY_min", 0.05, "Minimum DcaXY of Daughter Tracks"};

  // V0s
  Configurable<double> cfg_min_dca_V0_daughters{"cfg_min_dca_V0_daughters", 1.0, "min DCA between V0 daughters"};
  Configurable<double> cfg_min_dca_pos_to_PV{"cfg_min_dca_pos_to_PV", 0.1, "Minimum V0 Positive Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_neg_to_PV{"cfg_min_dca_neg_to_PV", 0.1, "Minimum V0 Negative Track DCAr cut to PV"};
  Configurable<double> cfg_min_dca_V0_to_PV{"cfg_min_dca_V0_to_PV", 0.6, "Minimum DCA V0 to PV"};
  Configurable<double> cfg_min_V0_radius{"cfg_min_V0_radius", 0.0, "Minimum V0 radius from PV"};
  Configurable<double> cfg_max_V0_radius{"cfg_max_V0_radius", 50.0, "Maximum V0 radius from PV"};
  Configurable<double> cfg_min_ctau{"cfg_min_ctau", 0.0, "Minimum ctau"};
  Configurable<double> cfg_max_ctau{"cfg_max_ctau", 50.0, "Maximum ctau"};
  Configurable<double> cfg_min_V0_cosPA{"cfg_min_V0_cosPA", 0.998, "Minimum V0 CosPA to PV"};
  Configurable<double> cfg_lambda_mass_window{"cfg_lambda_mass_window", 0.007, "Mass Window to select Lambda"};
  Configurable<double> cfg_kshort_rej{"cfg_kshort_rej", 0.005, "Reject K0Short Candidates"};
  Configurable<bool> cfg_kshort_rej_flag{"cfg_kshort_rej_flag", false, "K0short Mass Rej Flag"};
  Configurable<bool> cfg_armpod_flag{"cfg_armpod_flag", true, "Armentros-Podolanski Cut Flag"};
  Configurable<double> cfg_armpod_val{"cfg_armpod_val", 0.5, "Armentros-Podolanski Slope Parameter"};

  // V0s kinmatic acceptance
  Configurable<float> cfg_v0_pt_min{"cfg_v0_pt_min", 0.3, "Minimum V0 pT"};
  Configurable<float> cfg_v0_pt_max{"cfg_v0_pt_max", 3.0, "Minimum V0 pT"};
  Configurable<float> cfg_v0_rap_max{"cfg_v0_rap_max", 0.8, "|rap| cut"};

  // bool eta/rapidity
  Configurable<bool> cfg_do_eta_analysis{"cfg_do_eta_analysis", false, "Eta Analysis"};

  // V0s MC
  Configurable<bool> cfg_has_mc_flag{"cfg_has_mc_flag", true, "Has Mc Tag"};
  Configurable<bool> cfg_rec_primary_lambda{"cfg_rec_primary_lambda", false, "Primary Lambda"};
  Configurable<bool> cfg_rec_secondary_lambda{"cfg_rec_secondary_lambda", false, "Secondary Lambda"};
  Configurable<bool> cfg_rec_pid_flag{"cfg_rec_pid_flag", false, "PID Flag"};
  Configurable<bool> cfg_gen_primary_lambda{"cfg_gen_primary_lambda", true, "Primary Lambda"};
  Configurable<bool> cfg_gen_secondary_lambda{"cfg_gen_secondary_lambda", false, "Secondary Lambda"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisCol(6, 0, 6, "");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(200, 1.09, 1.14, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(64, 0.3, 3.5, "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(12, -0.6, 0.6, "y");
    const AxisSpec axisV0Eta(24, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., 2. * TMath::Pi(), "#phi (rad)");

    const AxisSpec axisRadius(400, 0, 200, "r(cm)");
    const AxisSpec axisCosPA(100, 0.99, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(1000, -5., 5., "dca (cm)");
    const AxisSpec axisDcaProngPV(1000, -50., 50., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisCTau(400, 0, 200, "c#tau (cm)");
    const AxisSpec axisGCTau(400, 0, 200, "#gammac#tau (cm)");
    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    const AxisSpec axisTrackPt(80, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackDCA(200, -1, 1, "dca_{XY} (cm)");
    const AxisSpec axisMomPID(80, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    // Create Histograms.
    // Event histograms
    histos.add("Event/h1d_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
    histos.add("Event/h1d_collisions_info", "# of Collisions", kTH1F, {axisCol});

    // QA
    histos.add("QA_Checks/h1d_tracks_info", "# of tracks", kTH1F, {axisCol});
    histos.add("QA_Checks/h1d_lambda_mass", "M_{#Lambda}", kTH1F, {axisV0Mass});
    histos.add("QA_Checks/h1d_antilambda_mass", "M_{#Lambda}", kTH1F, {axisV0Mass});
    histos.add("QA_Checks/h2d_lambda_pt_vs_mass", "p_{T} vs M_{#Lambda}", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("QA_Checks/h2d_antilambda_pt_vs_mass", "p_{T} vs M_{#Lambda}", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("QA_Checks/h2d_before_topo_cuts_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_after_topo_cuts_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_before_masswincut_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA_Checks/h2d_after_masswincut_pt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});

    // QA Lambda
    histos.add("QA_Sel_Lambda/h1d_V0_inv_mass", "V_{0} mass", kTH1F, {axisV0Mass});
    histos.add("QA_Sel_Lambda/h1d_V0_pt", "V_{0} p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_V0_eta", "#eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_V0_rap", "y-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_V0_phi", "#phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA_Sel_Lambda/h2d_V0_pt_vs_eta", "p_{T} vs #eta", kTH2F, {axisV0Eta, axisV0Pt});
    histos.add("QA_Sel_Lambda/h2d_V0_pt_vs_rap", "p_{T} vs y", kTH2F, {axisV0Eta, axisV0Pt});
    histos.add("QA_Sel_Lambda/h2d_V0_pt_vs_phi", "p_{T} vs #phi", kTH2F, {axisV0Phi, axisV0Pt});
    histos.add("QA_Sel_Lambda/h2d_V0_pt_vs_mass", "p_{T} vs M_{p#pi}", kTH2F, {axisV0Mass, axisV0Pt});

    histos.add("QA_Sel_Lambda/h1d_eta_check_V0_pt", "V_{0} p_{T}", kTH1F, {axisV0Pt});
    histos.add("QA_Sel_Lambda/h1d_eta_check_V0_eta", "#eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_eta_check_V0_rap", "y-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h2d_eta_check_V0_pt_vs_eta", "p_{T} vs #eta", kTH2F, {axisV0Eta, axisV0Pt});
    histos.add("QA_Sel_Lambda/h2d_eta_check_V0_pt_vs_rap", "p_{T} vs y", kTH2F, {axisV0Eta, axisV0Pt});

    histos.add("QA_Sel_Lambda/h1d_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA_Sel_Lambda/h1d_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA_Sel_Lambda/h1d_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA_Sel_Lambda/h1d_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA_Sel_Lambda/h1d_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA_Sel_Lambda/h1d_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});
    histos.add("QA_Sel_Lambda/h1d_V0_gctau", "V_{0} #gammac#tau", kTH1F, {axisGCTau});
    histos.add("QA_Sel_Lambda/h2d_qt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});

    histos.add("QA_Sel_Lambda/h1d_pos_prong_pt", "Pos-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_pt", "Neg-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA_Sel_Lambda/h1d_pos_prong_eta", "Pos-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_eta", "Neg-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_pos_prong_phi", "Pos-Prong #phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA_Sel_Lambda/h1d_neg_prong_phi", "Neg-Prong #phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA_Sel_Lambda/h1d_eta_check_pos_prong_pt", "Pos-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA_Sel_Lambda/h1d_eta_check_neg_prong_pt", "Neg-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA_Sel_Lambda/h1d_eta_check_pos_prong_eta", "Pos-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_eta_check_neg_prong_eta", "Neg-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA_Sel_Lambda/h1d_eta_check_pos_prong_phi", "Pos-Prong #phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA_Sel_Lambda/h1d_eta_check_neg_prong_phi", "Neg-Prong #phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA_Sel_Lambda/h2d_pos_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_dEdx_vs_p", "TPC Signal Pos-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_dEdx_vs_p", "TPC Signal Neg-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_pos_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA_Sel_Lambda/h2d_neg_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});

    // QA Anti-Lambda
    histos.addClone("QA_Sel_Lambda/", "QA_Sel_AntiLambda/");

    // MC Generated Histograms
    if (doprocessMCGen) {
      // McReco Histos
      histos.add("QA_Checks/h2d_lambda_pt_vs_eta_ps", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_lambda_pt_vs_eta_prim", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_lambda_pt_vs_eta_secd", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_antilambda_pt_vs_eta_ps", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_antilambda_pt_vs_eta_prim", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_antilambda_pt_vs_eta_secd", "", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("QA_Checks/h2d_tracks_pid_before_mccuts", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("QA_Checks/h2d_tracks_pid_phyprim", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("QA_Checks/h2d_tracks_pid_after_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("QA_Checks/h2d_lambda_from_sigma", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("QA_Checks/h2d_lambda_from_cascade", "PIDs", kTH2F, {axisPID, axisV0Pt});

      // McGen Histos
      histos.add("McGen/h1d_collisions_info", "# of collisions", kTH1F, {axisCol});
      histos.add("McGen/h1d_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
      histos.add("McGen/h1d_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1d_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1d_mass_lambda", "M_{#Lambda}", kTH1F, {axisV0Mass});
      histos.add("McGen/h1d_pt_lambda", "#Lambda p_{T}", kTH1F, {axisV0Pt});
      histos.add("McGen/h1d_eta_lambda", "#Lambda #eta-distribution", kTH1F, {axisV0Eta});
      histos.add("McGen/h1d_y_lambda", "#Lambda y-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_phi_lambda", "#Lambda #phi-distribution", kTH1F, {axisV0Phi});
      histos.add("McGen/h2d_pteta_lambda", "#Lambda p_{T} vs #eta", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("McGen/h2d_ptrap_lambda", "#Lambda p_{T} vs y", kTH2F, {axisV0Rap, axisV0Pt});
      histos.add("McGen/h2d_ptphi_lambda", "#Lambda p_{T} vs #phi", kTH2F, {axisV0Phi, axisV0Pt});
      histos.add("McGen/h1d_mass_antilambda", "M_{#bar{#Lambda}}", kTH1F, {axisV0Mass});
      histos.add("McGen/h1d_pt_antilambda", "#bar{#Lambda} p_{T}", kTH1F, {axisV0Pt});
      histos.add("McGen/h1d_eta_antilambda", "#bar{#Lambda} #eta-distribution", kTH1F, {axisV0Eta});
      histos.add("McGen/h1d_y_antilambda", "#bar{#Lambda} y-distribution", kTH1F, {axisV0Rap});
      histos.add("McGen/h1d_phi_antilambda", "#bar{#Lambda} #phi-distribution", kTH1F, {axisV0Phi});
      histos.add("McGen/h2d_pteta_antilambda", "#bar{#Lambda} p_{T} vs #eta", kTH2F, {axisV0Eta, axisV0Pt});
      histos.add("McGen/h2d_ptrap_antilambda", "#bar{#Lambda} p_{T} vs y", kTH2F, {axisV0Rap, axisV0Pt});
      histos.add("McGen/h2d_ptphi_antilambda", "#bar{#Lambda} p_{T} vs #phi", kTH2F, {axisV0Phi, axisV0Pt});
    }
  }

  template <typename C>
  bool selCol(C const& col)
  {

    if (col.posZ() < cfg_min_z_vtx || col.posZ() >= cfg_max_z_vtx) {
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

  template <typename T>
  bool dauTrackSelection(T const& track)
  {

    if (track.pt() < cfg_pt_min || track.pt() > cfg_pt_max) {
      return false;
    }

    if (fabs(track.eta()) >= cfg_eta_cut) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < cfg_min_crossed_rows) {
      return false;
    }

    if (fabs(track.dcaXY()) < cfg_track_dcaXY_min) {
      return false;
    }

    return true;
  }

  template <typename C, typename V, typename T>
  bool topologicalCutsV0(C const& col, V const& v0, T const&)
  {

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (!dauTrackSelection(postrack) || !dauTrackSelection(negtrack)) {
      return false;
    }

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

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau < cfg_min_ctau || ctau > cfg_max_ctau) {
      return false;
    }

    if (v0.v0cosPA() < cfg_min_V0_cosPA) {
      return false;
    }

    return true;
  }

  template <ParticleType part, typename T>
  bool selPIDTrack(T const& postrack, T const& negtrack)
  {
    bool return_flag = false;
    float tpcNSigmaPr = 0., tpcNSigmaPi = 0.;

    switch (part) {
      // postrack = Proton, negtrack = Pion
      case kLambda:
        tpcNSigmaPr = postrack.tpcNSigmaPr();
        tpcNSigmaPi = negtrack.tpcNSigmaPi();
        break;

      // negtrack = Proton, postrack = Pion
      case kAntiLambda:
        tpcNSigmaPr = negtrack.tpcNSigmaPr();
        tpcNSigmaPi = postrack.tpcNSigmaPi();
        break;
    }

    if (fabs(tpcNSigmaPr) < cfg_tpc_nsigma && fabs(tpcNSigmaPi) < cfg_tpc_nsigma) {
      return_flag = true;
    }

    return return_flag;
  }

  template <typename T>
  bool checkKinCuts(T const& v0track)
  {

    // pT cut
    if (v0track.pt() <= cfg_v0_pt_min || v0track.pt() >= cfg_v0_pt_max) {
      return false;
    }

    // apply rapidity/pseudorapidity acceptance
    float rap = 0.;
    if (cfg_do_eta_analysis) {
      rap = v0track.eta();
    } else {
      rap = v0track.yLambda();
    }

    if (fabs(rap) >= cfg_v0_rap_max) {
      return false;
    }

    return true;
  }

  template <typename T>
  void getPDGsIDs(T const& mcparts, std::vector<int>& PDGs, std::vector<int>& IDs)
  {

    for (auto mcpart : mcparts) {
      if (mcpart.globalIndex() != 0) {
        IDs.push_back(mcpart.globalIndex());
        PDGs.push_back(mcpart.pdgCode());
      }
    }

    while (IDs.size() > 2) {
      IDs.pop_back();
    }

    while (PDGs.size() > 2) {
      PDGs.pop_back();
    }
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillQALambda(C const& col, V const& v0, T const&)
  {

    static constexpr std::string_view sub_dir[] = {"QA_Sel_Lambda/", "QA_Sel_AntiLambda/"};

    // daugthers
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();
    float mass = 0.;

    if constexpr (part == kLambda) {
      mass = v0.mLambda();
    } else {
      mass = v0.mAntiLambda();
    }

    // ctau
    float e = RecoDecay::e(v0.px(), v0.py(), v0.pz(), mass);
    float gamma = e / mass;
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    float gctau = ctau * gamma;

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_inv_mass"), mass);
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_pt"), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_eta"), v0.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_rap"), v0.yLambda());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_phi"), v0.phi());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_V0_pt_vs_eta"), v0.eta(), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_V0_pt_vs_rap"), v0.yLambda(), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_V0_pt_vs_phi"), v0.phi(), v0.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_V0_pt_vs_mass"), mass, v0.pt());

    if (fabs(v0.eta()) > 0.8) {
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_V0_pt"), v0.pt());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_V0_eta"), v0.eta());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_V0_rap"), v0.yLambda());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_eta_check_V0_pt_vs_eta"), v0.eta(), v0.pt());
      histos.fill(HIST(sub_dir[part]) + HIST("h2d_eta_check_V0_pt_vs_rap"), v0.yLambda(), v0.pt());
    }

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_pos_to_PV"), v0.dcapostopv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_neg_to_PV"), v0.dcanegtopv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_dca_V0_to_PV"), v0.dcav0topv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_radius"), v0.v0radius());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_ctau"), ctau);
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_V0_gctau"), gctau);
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_qt_vs_alpha"), v0.alpha(), v0.qtarm());

    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_pt"), postrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_pos_prong_phi"), postrack.phi());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1d_neg_prong_phi"), negtrack.phi());

    if (fabs(v0.eta()) > 0.8) {
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_pos_prong_pt"), postrack.pt());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_pos_prong_eta"), postrack.eta());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_pos_prong_phi"), postrack.phi());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_neg_prong_pt"), negtrack.pt());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_neg_prong_eta"), negtrack.eta());
      histos.fill(HIST(sub_dir[part]) + HIST("h1d_eta_check_neg_prong_phi"), negtrack.phi());
    }

    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_dcaXY_vs_pt"), postrack.pt(), postrack.dcaXY());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_dcaXY_vs_pt"), negtrack.pt(), negtrack.dcaXY());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_tpc_nsigma_pr_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_tpc_nsigma_pr_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_pos_prong_tpc_nsigma_pi_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(sub_dir[part]) + HIST("h2d_neg_prong_tpc_nsigma_pi_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
  }

  template <bool data, bool reco, typename C, typename V, typename T>
  void selV0Particle(C const& collision, V const& v0track, T const& tracks)
  {

    // apply kinematic cuts
    if (!checkKinCuts(v0track)) {
      return;
    }

    // initialize daughter tracks
    auto postrack = v0track.template posTrack_as<T>();
    auto negtrack = v0track.template negTrack_as<T>();

    // initialize mass and v0lambda/v0antilambda
    float mass = 0.;
    ParticleType v0part;

    // apply daughter particle id
    // check for Lambda / Anti-Lambda
    if (selPIDTrack<kLambda>(postrack, negtrack)) {
      mass = v0track.mLambda();
      v0part = kLambda;
    } else if (selPIDTrack<kAntiLambda>(postrack, negtrack)) {
      mass = v0track.mAntiLambda();
      v0part = kAntiLambda;
    } else {
      return;
    }

    histos.fill(HIST("QA_Checks/h1d_tracks_info"), 3.5);
    histos.fill(HIST("QA_Checks/h2d_before_masswincut_pt_vs_alpha"), v0track.alpha(), v0track.qtarm());

    // apply Armentros-Podolanski Selection
    if (cfg_armpod_flag && (fabs(v0track.alpha()) < v0track.qtarm() / cfg_armpod_val)) {
      return;
    }

    // apply kshort rejection hypothesis
    if (cfg_kshort_rej_flag && (fabs(v0track.mK0Short() - MassK0Short) <= cfg_kshort_rej)) {
      return;
    }

    // fill mass histograms before applying mass window cut to get % purity
    if (v0part == kLambda) {
      histos.fill(HIST("QA_Checks/h1d_lambda_mass"), mass);
      histos.fill(HIST("QA_Checks/h2d_lambda_pt_vs_mass"), mass, v0track.pt());
    } else if (v0part == kAntiLambda) {
      histos.fill(HIST("QA_Checks/h1d_antilambda_mass"), mass);
      histos.fill(HIST("QA_Checks/h2d_antilambda_pt_vs_mass"), mass, v0track.pt());
    }

    // apply mass window cut (Selection of Lambda/AntiLambda)
    if (fabs(mass - MassLambda0) >= cfg_lambda_mass_window) {
      return;
    }

    histos.fill(HIST("QA_Checks/h2d_after_masswincut_pt_vs_alpha"), v0track.alpha(), v0track.qtarm());
    histos.fill(HIST("QA_Checks/h1d_tracks_info"), 4.5);

    // MC Reco Analysis
    if constexpr (reco) {
      auto v0mcpart = v0track.mcParticle();

      // Get information of all the reconstructed V0s
      histos.fill(HIST("QA_Checks/h2d_tracks_pid_before_mccuts"), v0mcpart.pdgCode(), v0mcpart.pt());

      // Get all lambdas/anti-lambdas before any selection
      if (v0mcpart.pdgCode() == kLambda0) {
        histos.fill(HIST("QA_Checks/h2d_lambda_pt_vs_eta_ps"), v0mcpart.eta(), v0mcpart.pt());
      } else if (v0mcpart.pdgCode() == kLambda0Bar) {
        histos.fill(HIST("QA_Checks/h2d_antilambda_pt_vs_eta_ps"), v0mcpart.eta(), v0mcpart.pt());
      }

      // get all primary reconstructed V0s
      if (v0mcpart.isPhysicalPrimary()) {
        histos.fill(HIST("QA_Checks/h2d_tracks_pid_phyprim"), v0mcpart.pdgCode(), v0mcpart.pt());
      }

      // Get Daughters and Mothers
      bool decay_channel_flag = false;
      std::vector<int> daughterPDGs{}, daughterIDs{}, motherPDGs{}, motherIDs{};
      auto mcpart_daughters = v0mcpart.template daughters_as<aod::McParticles>();
      auto mcpart_mothers = v0mcpart.template mothers_as<aod::McParticles>();

      if (cfg_rec_pid_flag) {

        if (v0part == kLambda && v0mcpart.pdgCode() != kLambda0) {
          return;
        } else if (v0part == kAntiLambda && v0mcpart.pdgCode() != kLambda0Bar) {
          return;
        }

        getPDGsIDs(mcpart_daughters, daughterPDGs, daughterIDs);
        getPDGsIDs(mcpart_mothers, motherPDGs, motherIDs);

        // Decay to Proton-Pion
        if (abs(daughterPDGs[0]) == kProton && abs(daughterPDGs[1]) == kPiPlus) {
          decay_channel_flag = true;
        }

        // check for correct decay channel
        if (!decay_channel_flag) {
          return;
        }

        // fill histograms to get secondary contaminations
        if (v0mcpart.pdgCode() == kLambda0) {
          if (v0mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("QA_Checks/h2d_lambda_pt_vs_eta_prim"), v0mcpart.eta(), v0mcpart.pt());
          } else {
            histos.fill(HIST("QA_Checks/h2d_lambda_pt_vs_eta_secd"), v0mcpart.eta(), v0mcpart.pt());
          }
        } else if (v0mcpart.pdgCode() == kLambda0Bar) {
          if (v0mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("QA_Checks/h2d_antilambda_pt_vs_eta_prim"), v0mcpart.eta(), v0mcpart.pt());
          } else {
            histos.fill(HIST("QA_Checks/h2d_antilambda_pt_vs_eta_secd"), v0mcpart.eta(), v0mcpart.pt());
          }
        }

        // check whether the selected lambda is a Physical Primary / Secondary
        if (cfg_rec_primary_lambda && !v0mcpart.isPhysicalPrimary()) {
          return;
        } else if (cfg_rec_secondary_lambda && v0mcpart.isPhysicalPrimary()) {
          return;
        }

        // check the secondary lambdas coming from Sigma, Cascades and Omegas
        if (abs(motherPDGs[0]) == kSigma0 || abs(motherPDGs[1]) == kSigma0Bar) {
          histos.fill(HIST("QA_Checks/h2d_lambda_from_sigma"), v0mcpart.pdgCode(), v0mcpart.pt());
        }

        if (abs(motherPDGs[0]) == kXiMinus || abs(motherPDGs[1]) == kXiMinus) {
          histos.fill(HIST("QA_Checks/h2d_lambda_from_cascade"), v0mcpart.pdgCode(), v0mcpart.pt());
        }
      }

      // Fill the counter for selected primary/secondary Lambda/AntiLambda
      histos.fill(HIST("QA_Checks/h1d_tracks_info"), 5.5);
      histos.fill(HIST("QA_Checks/h2d_tracks_pid_after_sel"), v0mcpart.pdgCode(), v0mcpart.pt());
    }

    if (v0part == kLambda) {
      fillQALambda<kLambda>(collision, v0track, tracks);
    } else {
      fillQALambda<kAntiLambda>(collision, v0track, tracks);
    }

    // Fill Lambda/AntiLambda Table
    lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0track.px(), v0track.py(), v0track.pz(), v0track.pt(), v0track.eta(), v0track.phi(), v0track.yLambda(), mass, postrack.index(), negtrack.index(), (int8_t)v0part, v0track.v0cosPA(), v0track.dcaV0daughters());
  }

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

  void processData(Collisions::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {

    histos.fill(HIST("Event/h1d_collisions_info"), 1.5);

    // select collision
    if (!selCol(collision)) {
      return;
    }

    histos.fill(HIST("Event/h1d_collisions_info"), 2.5);
    histos.fill(HIST("Event/h1d_collision_posZ"), collision.posZ());

    lambdaCollisionTable(collision.centFT0M(), collision.posX(), collision.posY(), collision.posZ());

    for (auto const& v0 : V0s) {

      histos.fill(HIST("QA_Checks/h1d_tracks_info"), 1.5);
      histos.fill(HIST("QA_Checks/h2d_before_topo_cuts_pt_vs_alpha"), v0.alpha(), v0.qtarm());

      // apply topological cuts on v0 candidates
      if (!topologicalCutsV0(collision, v0, tracks)) {
        continue;
      }

      histos.fill(HIST("QA_Checks/h1d_tracks_info"), 2.5);
      histos.fill(HIST("QA_Checks/h2d_after_topo_cuts_pt_vs_alpha"), v0.alpha(), v0.qtarm());

      selV0Particle<true, false>(collision, v0, tracks);
    }
  }

  PROCESS_SWITCH(lambdaCorrTableProducer, processData, "Process for DATA", true);

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

    for (auto const& collision : collisions) {

      histos.fill(HIST("Event/h1d_collisions_info"), 1.5);

      // select collision
      if (!selCol(collision)) {
        continue;
      }

      histos.fill(HIST("Event/h1d_collisions_info"), 2.5);

      // check for corresponding MCGen Collision
      if (cfg_has_mc_flag && !collision.has_mcCollision()) {
        continue;
      }

      histos.fill(HIST("Event/h1d_collisions_info"), 3.5);
      histos.fill(HIST("Event/h1d_collision_posZ"), collision.posZ());

      lambdaCollisionTable(collision.centFT0M(), collision.posX(), collision.posY(), collision.posZ());

      // auto const& mcCollision = collision.mcCollision_as<aod::McCollisions::iterator>();

      // v0-track loop
      auto v0sThisCollision = V0s.sliceBy(perCol, collision.globalIndex());
      for (auto const& v0 : v0sThisCollision) {

        histos.fill(HIST("QA_Checks/h1d_tracks_info"), 0.5);

        // check for corresponding MCGen Particle
        if (cfg_has_mc_flag && !v0.has_mcParticle()) {
          continue;
        }

        histos.fill(HIST("QA_Checks/h1d_tracks_info"), 1.5);
        histos.fill(HIST("QA_Checks/h2d_before_topo_cuts_pt_vs_alpha"), v0.alpha(), v0.qtarm());

        // apply topological cuts on v0 candidates
        if (!topologicalCutsV0(collision, v0, tracks)) {
          continue;
        }

        histos.fill(HIST("QA_Checks/h1d_tracks_info"), 2.5);
        histos.fill(HIST("QA_Checks/h2d_after_topo_cuts_pt_vs_alpha"), v0.alpha(), v0.qtarm());

        selV0Particle<false, true>(collision, v0, tracks);
      }
    }
  }

  PROCESS_SWITCH(lambdaCorrTableProducer, processMCReco, "Process for MC Reconstructed", false);

  void processMCGen(McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {

    histos.fill(HIST("McGen/h1d_collisions_info"), 1.5);

    // apply collision cuts
    if (mcCollision.posZ() < cfg_min_z_vtx || mcCollision.posZ() > cfg_max_z_vtx) {
      return;
    }

    histos.fill(HIST("McGen/h1d_collisions_info"), 2.5);
    histos.fill(HIST("McGen/h1d_collision_posZ"), mcCollision.posZ());
    lambdaMCGenCollisionTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());
    float mass = 0.;

    for (auto const& mcpart : mcParticles) {

      // check for Lambda first
      if (abs(mcpart.pdgCode()) != kLambda0) {
        continue;
      }

      // check for Primary Lambdas/AntiLambdas
      if (cfg_gen_primary_lambda && !mcpart.isPhysicalPrimary()) {
        continue;
      } else if (cfg_gen_secondary_lambda && mcpart.isPhysicalPrimary()) {
        continue;
      }

      // apply kinematic acceptance
      if (mcpart.pt() <= cfg_v0_pt_min || mcpart.pt() >= cfg_v0_pt_max) {
        continue;
      }

      float rap = 0.;
      if (cfg_do_eta_analysis) {
        rap = mcpart.eta();
      } else {
        rap = mcpart.y();
      }

      if (fabs(rap) >= cfg_v0_rap_max) {
        continue;
      }

      // Get Daughters and Mothers
      bool decay_channel_flag = false;
      std::vector<int> daughterPDGs{}, daughterIDs{}, motherPDGs{}, motherIDs{};
      auto mcpart_daughters = mcpart.template daughters_as<aod::McParticles>();
      auto mcpart_mothers = mcpart.template mothers_as<aod::McParticles>();
      getPDGsIDs(mcpart_daughters, daughterPDGs, daughterIDs);
      getPDGsIDs(mcpart_mothers, motherPDGs, motherIDs);

      // Decay to Proton-Pion
      if (abs(daughterPDGs[0]) == kProton && abs(daughterPDGs[1]) == kPiPlus) {
        decay_channel_flag = true;
      }

      if (!decay_channel_flag) {
        continue;
      }

      mass = RecoDecay::m(mcpart.p(), mcpart.e());

      // Fill histograms
      if (mcpart.pdgCode() == kLambda0) {
        histos.fill(HIST("McGen/h1d_lambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1d_lambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1d_mass_lambda"), mass);
        histos.fill(HIST("McGen/h1d_pt_lambda"), mcpart.pt());
        histos.fill(HIST("McGen/h1d_eta_lambda"), mcpart.eta());
        histos.fill(HIST("McGen/h1d_y_lambda"), mcpart.y());
        histos.fill(HIST("McGen/h1d_phi_lambda"), mcpart.phi());
        histos.fill(HIST("McGen/h2d_pteta_lambda"), mcpart.eta(), mcpart.pt());
        histos.fill(HIST("McGen/h2d_ptrap_lambda"), mcpart.y(), mcpart.pt());
        histos.fill(HIST("McGen/h2d_ptphi_lambda"), mcpart.phi(), mcpart.pt());
        lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(), mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), mass, daughterIDs[0], daughterIDs[1], (int8_t)kLambda, -999., -999.);
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        histos.fill(HIST("McGen/h1d_antilambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1d_antilambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1d_mass_antilambda"), mass);
        histos.fill(HIST("McGen/h1d_pt_antilambda"), mcpart.pt());
        histos.fill(HIST("McGen/h1d_eta_antilambda"), mcpart.eta());
        histos.fill(HIST("McGen/h1d_y_antilambda"), mcpart.y());
        histos.fill(HIST("McGen/h1d_phi_antilambda"), mcpart.phi());
        histos.fill(HIST("McGen/h2d_pteta_antilambda"), mcpart.eta(), mcpart.pt());
        histos.fill(HIST("McGen/h2d_ptrap_antilambda"), mcpart.y(), mcpart.pt());
        histos.fill(HIST("McGen/h2d_ptphi_antilambda"), mcpart.phi(), mcpart.pt());
        lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(), mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), mass, daughterIDs[1], daughterIDs[0], (int8_t)kAntiLambda, -999., -999.);
      }
    }
  }

  PROCESS_SWITCH(lambdaCorrTableProducer, processMCGen, "Process for MC Generated", false);
};

struct lambdaCorrelationAnalysis {

  // Global Configurables
  Configurable<int> cfg_nRapBins{"cfg_nRapBins", 12, "N Rapidity Bins"};
  Configurable<float> cfg_Rap_Min{"cfg_Rap_Min", -0.6, "Minimum Rapidity"};
  Configurable<float> cfg_Rap_Max{"cfg_Rap_Max", 0.6, "Maximum Rapidity"};
  Configurable<int> cfg_nPhiBins{"cfg_nPhiBins", 64, "N Phi Bins"};
  Configurable<float> cfg_Phi_Min{"cfg_Phi_Min", 0, "Minimum Phi"};
  Configurable<float> cfg_Phi_Max{"cfg_Phi_Max", 2 * TMath::Pi(), "Maximum Phi"};

  // remove lambda with shared daughters
  Configurable<bool> cfg_remove_lambda{"cfg_remove_lambda", true, "Flag to remove lambda"};

  // pt dependent efficiencies
  Configurable<bool> cfg_eff_corr_flag{"cfg_eff_corr_flag", true, "Efficiency Correction Flag"};
  Configurable<std::vector<float>> cfg_pt_bins{"cfg_pt_bins", {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5}, "pT bins for efficiencies"};
  Configurable<std::vector<float>> cfg_lambda_eff{"cfg_lambda_eff", {0.00480, 0.03529, 0.07232, 0.10678, 0.13147, 0.14560, 0.15536, 0.16519, 0.17139, 0.17817, 0.18111, 0.18234, 0.18192, 0.18211, 0.17968}, "Lambda Efficiencies"};
  Configurable<std::vector<float>> cfg_antilambda_eff{"cfg_antilambda_eff", {0.00414, 0.03219, 0.06641, 0.09892, 0.12297, 0.13752, 0.14767, 0.15738, 0.16584, 0.17065, 0.17625, 0.17741, 0.17681, 0.17844, 0.17354}, "AntiLambda Efficiencies"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Initialize global variables
  float nrapbins = 0.;
  float kminrap = 0.;
  float kmaxrap = 0.;
  float nphibins = 0.;
  float kminphi = 0.;
  float kmaxphi = 0.;
  float rapbinwidth = 0.;
  float phibinwidth = 0.;
  float q = 0., e = 0., qinv = 0.;

  std::vector<float> v_pt_bins = static_cast<std::vector<float>>(cfg_pt_bins);
  std::vector<float> v_lambda_eff = static_cast<std::vector<float>>(cfg_lambda_eff);
  std::vector<float> v_antilambda_eff = static_cast<std::vector<float>>(cfg_antilambda_eff);

  void init(InitContext const&)
  {
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

    const AxisSpec axisCheck(1, 0, 1, "");
    const AxisSpec axisPosZ(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisMass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisPt(64, 0.3, 3.5, "p_{T} (GeV/#it{c})");
    const AxisSpec axisEta(24, -1.2, 1.2, "#eta");
    const AxisSpec axisCPA(100, 0.99, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisRap(cfg_nRapBins, cfg_Rap_Min, cfg_Rap_Max, "y");
    const AxisSpec axisPhi(cfg_nPhiBins, cfg_Phi_Min, cfg_Phi_Max, "#phi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "y #phi");
    const AxisSpec axisQinv(100, 0, 10, "q_{inv} (GeV/#it{c})");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1d_collision_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/Reco/h1d_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h1d_lambda_multiplicity", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_antilambda_multiplicity", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_lambda_sdau", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_antilambda_sdau", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_lambda_totmult", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1d_antilambda_totmult", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});

    // InvMass, DcaDau and CosPA
    histos.add("Reco/QA_Lambda/h1d_V0_mass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/QA_Lambda/h1d_V0_cpa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/QA_Lambda/h1d_V0_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});

    histos.addClone("Reco/QA_Lambda/", "Reco/QA_AntiLambda/");

    // single and two particle densities
    // 1D Histograms
    histos.add("Reco/h1d_n1_pt_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisPt});
    histos.add("Reco/h1d_n1_pt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisPt});
    histos.add("Reco/h1d_n1_eta_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_eta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_rap_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_rap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_phi_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisPhi});
    histos.add("Reco/h1d_n1_phi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisPhi});
    histos.add("Reco/h2d_n1_pt_vs_eta_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisEta, axisPt});
    histos.add("Reco/h2d_n1_pt_vs_eta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisEta, axisPt});
    histos.add("Reco/h1d_n2_qinv_LaP_LaM", "#rho_{2}^{#Lambda-#bar{#Lambda}}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaP_LaP", "#rho_{2}^{#Lambda-#Lambda}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaM_LaM", "#rho_{2}^{#bar{#Lambda}-#bar{#Lambda}}", kTH1D, {axisQinv});

    // 2D Histograms
    histos.add("Reco/h2d_n1_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n1_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n2_LaP_LaM", "#rho_{2}^{#Lambda - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_LaP_LaP", "#rho_{2}^{#Lambda - #Lambda}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_LaM_LaM", "#rho_{2}^{#bar{#Lambda} - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_pt1pt2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_pt1pt2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_pt1pt2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisPt});
    histos.add("Reco/h2d_n2_eta1eta2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisEta, axisEta});
    histos.add("Reco/h2d_n2_eta1eta2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisEta, axisEta});
    histos.add("Reco/h2d_n2_eta1eta2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisEta, axisEta});
    histos.add("Reco/h2d_n2_phi1phi2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_phi1phi2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_phi1phi2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPhi, axisPhi});
    histos.add("Reco/h2d_n2_rap1rap2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_rap1rap2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_rap1rap2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisRap, axisRap});
    histos.add("Reco/h2d_n2_pt1eta2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisEta});
    histos.add("Reco/h2d_n2_pt1eta2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisEta});
    histos.add("Reco/h2d_n2_pt1eta2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisEta});
    histos.add("Reco/h2d_n2_pt1phi2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisPhi});
    histos.add("Reco/h2d_n2_pt1phi2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisPhi});
    histos.add("Reco/h2d_n2_pt1phi2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisPhi});
    histos.add("Reco/h2d_n2_pt1rap2_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n2_pt1rap2_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n2_pt1rap2_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisPt, axisRap});

    // MCGen
    if (doprocessMCGen) {
      histos.addClone("Event/Reco/", "Event/McGen/");
      histos.addClone("Reco/", "McGen/");
    }
  }

  template <bool fillHist, ParticleType part, RecGenType rec_gen, typename T, typename V>
  bool removeLambdaSharingDau(T const& v, V const& vs)
  {
    // check whether to remove lambda or not
    if (!cfg_remove_lambda) {
      return true;
    }

    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir[] = {"QA_Lambda/", "QA_AntiLambda/"};

    bool ret_flag = true;

    for (auto const& x : vs) {
      if ((v.index() != x.index()) && (v.postrackid() == x.postrackid() || v.negtrackid() == x.negtrackid())) {
        if (fillHist) {
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1d_V0_mass"), x.mass());
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1d_V0_cpa"), x.cospa());
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1d_V0_dcadau"), x.dcadau());
        }
        if (std::abs(v.mass() - MassLambda0) > std::abs(x.mass() - MassLambda0)) {
          ret_flag = false;
          break;
        }
      }
    }

    return ret_flag;
  }

  template <ParticleType part>
  void get_corr_factor(float& eff, float pt)
  {

    if (!cfg_eff_corr_flag) {
      eff = 1.;
      return;
    }

    int n = v_lambda_eff.size();

    for (int i = 0; i < n; ++i) {
      if (pt > v_pt_bins[i] && pt <= v_pt_bins[i + 1]) {
        if constexpr (part == kLambda) {
          eff = v_lambda_eff[i];
        } else if constexpr (part == kAntiLambda) {
          eff = v_antilambda_eff[i];
        }
      } else {
        eff = 1.;
      }
    }
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, typename U>
  void fillPairHistos(U& p1, U& p2)
  {

    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir_hist[] = {"LaP_LaM", "LaP_LaP", "LaM_LaM"};

    int rapbin1 = static_cast<int>((p1.rap() - kminrap) / rapbinwidth);
    int rapbin2 = static_cast<int>((p2.rap() - kminrap) / rapbinwidth);

    int phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    int phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    float eff_1 = 1., eff_2 = 1.;

    if constexpr (rec_gen == kRec) {
      if (part_pair == kLambdaAntiLambda) {
        get_corr_factor<kLambda>(eff_1, p1.pt());
        get_corr_factor<kAntiLambda>(eff_2, p2.pt());
      } else if (part_pair == kLambdaLambda) {
        get_corr_factor<kLambda>(eff_1, p1.pt());
        get_corr_factor<kLambda>(eff_2, p2.pt());
      } else if (part_pair == kAntiLambdaAntiLambda) {
        get_corr_factor<kAntiLambda>(eff_1, p1.pt());
        get_corr_factor<kAntiLambda>(eff_2, p2.pt());
      }
    }

    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1pt2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.pt(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_eta1eta2_") + HIST(sub_dir_hist[part_pair]), p1.eta(), p2.eta(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_phi1phi2_") + HIST(sub_dir_hist[part_pair]), p1.phi(), p2.phi(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_rap1rap2_") + HIST(sub_dir_hist[part_pair]), p1.rap(), p2.rap(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1eta2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.eta(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1phi2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.phi(), eff_1 * eff_2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1rap2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.rap(), eff_1 * eff_2);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_") + HIST(sub_dir_hist[part_pair]), rapphix + 0.5, rapphiy + 0.5, eff_1 * eff_2);
    }

    // qinv histos
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n2_qinv_") + HIST(sub_dir_hist[part_pair]), qinv, eff_1 * eff_2);
  }

  template <ParticleType part, RecGenType rec_gen, typename T>
  void analyzeSingles(T const& tracks)
  {

    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir_hist[] = {"LaP", "LaM"};

    int ntrk1 = 0, ntrk2 = 0, ntrk3 = 0;
    float eff = 1.;

    for (auto const& track : tracks) {
      ++ntrk1;
      if (!removeLambdaSharingDau<true, part, rec_gen>(track, tracks)) {
        ++ntrk2;
        continue;
      }
      if constexpr (rec_gen == kRec) {
        get_corr_factor<part>(eff, track.pt());
      }
      ++ntrk3;

      // QA and Efficiency Calculation Plots
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_pt_") + HIST(sub_dir_hist[part]), track.pt(), eff);
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_eta_") + HIST(sub_dir_hist[part]), track.eta(), eff);
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_phi_") + HIST(sub_dir_hist[part]), track.phi(), eff);
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_rap_") + HIST(sub_dir_hist[part]), track.rap(), eff);
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_pt_vs_eta_") + HIST(sub_dir_hist[part]), track.eta(), track.pt(), eff);

      // Rho1 for R2 Calculation
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_") + HIST(sub_dir_hist[part]), track.rap(), track.phi(), eff);
    }

    // fill multiplicity histograms
    if (ntrk1 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_lambda_totmult"), ntrk1);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_antilambda_totmult"), ntrk1);
      }
    }

    if (ntrk2 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_lambda_sdau"), ntrk2);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_antilambda_sdau"), ntrk2);
      }
    }

    if (ntrk3 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_lambda_multiplicity"), ntrk3);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_antilambda_multiplicity"), ntrk3);
      }
    }
  }

  template <ParticlePairType partpair, RecGenType rec_gen, bool samelambda, typename T>
  void analyzePairs(T const& trks_1, T const& trks_2)
  {
    for (auto const& trk_1 : trks_1) {
      if (!removeLambdaSharingDau<false, kLambda, rec_gen>(trk_1, trks_1)) {
        continue;
      }
      for (auto const& trk_2 : trks_2) {
        // check for same index for Lambda-Lambda / AntiLambda-AntiLambda
        if (samelambda && ((trk_1.index() == trk_2.index()))) {
          continue;
        }
        // check if Lambda shares a daughter and select the one closest to PDG Mass
        if (!removeLambdaSharingDau<false, kLambda, rec_gen>(trk_2, trks_2)) {
          continue;
        }
        fillPairHistos<partpair, rec_gen>(trk_1, trk_2);
      }
    }
  }

  using Lambda_Collisions = aod::LambdaCollisions;
  using Lambda_Tracks = aod::LambdaTracks;

  SliceCache cache;

  Partition<Lambda_Tracks> part_lambda_tracks = (aod::lambdatrack::v0type == (int8_t)kLambda);
  Partition<Lambda_Tracks> part_anti_lambda_tracks = (aod::lambdatrack::v0type == (int8_t)kAntiLambda);

  void processDataReco(Lambda_Collisions::iterator const& collision, Lambda_Tracks const&)
  {

    histos.fill(HIST("Event/Reco/h1d_collision_posz"), collision.posZ());
    histos.fill(HIST("Event/Reco/h1d_ft0m_mult_percentile"), collision.cent());

    auto lambda_tracks = part_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto anti_lambda_tracks = part_anti_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kRec>(lambda_tracks);
    analyzeSingles<kAntiLambda, kRec>(anti_lambda_tracks);
    analyzePairs<kLambdaAntiLambda, kRec, false>(lambda_tracks, anti_lambda_tracks);
    analyzePairs<kLambdaLambda, kRec, true>(lambda_tracks, lambda_tracks);
    analyzePairs<kAntiLambdaAntiLambda, kRec, true>(anti_lambda_tracks, anti_lambda_tracks);
  }

  PROCESS_SWITCH(lambdaCorrelationAnalysis, processDataReco, "Process for Data and MCReco", true);

  using LambdaMCGenCollisions = aod::LambdaMCGenCollisions;
  using LambdaMCGenTracks = aod::LambdaMCGenTracks;

  SliceCache cachemc;
  Partition<LambdaMCGenTracks> part_lambda_mcgen_tracks = aod::lambdatrack::v0type == (int8_t)kLambda;
  Partition<LambdaMCGenTracks> part_antilambda_mcgen_tracks = aod::lambdatrack::v0type == (int8_t)kAntiLambda;

  void processMCGen(LambdaMCGenCollisions::iterator const& mcgencol, LambdaMCGenTracks const&)
  {

    histos.fill(HIST("Event/McGen/h1d_collision_posz"), mcgencol.posZ());

    auto lambda_mcgen_tracks = part_lambda_mcgen_tracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto antilambda_mcgen_tracks = part_antilambda_mcgen_tracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);

    analyzeSingles<kLambda, kGen>(lambda_mcgen_tracks);
    analyzeSingles<kAntiLambda, kGen>(antilambda_mcgen_tracks);
    analyzePairs<kLambdaAntiLambda, kGen, false>(lambda_mcgen_tracks, antilambda_mcgen_tracks);
    analyzePairs<kLambdaLambda, kGen, true>(lambda_mcgen_tracks, lambda_mcgen_tracks);
    analyzePairs<kAntiLambdaAntiLambda, kGen, true>(antilambda_mcgen_tracks, antilambda_mcgen_tracks);
  }

  PROCESS_SWITCH(lambdaCorrelationAnalysis, processMCGen, "Process for MC Generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdaCorrTableProducer>(cfgc),
    adaptAnalysisTask<lambdaCorrelationAnalysis>(cfgc)};
}
