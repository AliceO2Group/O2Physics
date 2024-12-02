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
#include <string>

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
#include "CCDB/BasicCCDBManager.h"
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
DECLARE_SOA_COLUMN(CorrFact, corrfact, float);
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
                  lambdatrack::DcaDau,
                  lambdatrack::CorrFact);
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
                  lambdatrack::DcaDau,
                  lambdatrack::CorrFact);
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

enum DMCType {
  kData = 0,
  kMC
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

  // V0s MC
  Configurable<bool> cfg_has_mc_flag{"cfg_has_mc_flag", true, "Has Mc Tag"};
  Configurable<bool> cfg_sel_true_lambda{"cfg_sel_true_lambda", true, "Select True Lambda"};
  Configurable<bool> cfg_rec_primary_lambda{"cfg_rec_primary_lambda", true, "Primary Lambda"};
  Configurable<bool> cfg_gen_primary_lambda{"cfg_gen_primary_lambda", true, "Primary Lambda"};

  // Efficiency Correction
  Configurable<bool> cfg_corr_flag{"cfg_corr_flag", true, "Efficiency Correction Flag"};
  Configurable<int> cfg_cf_hist_type{"cfg_cf_hist_type", 0, "Histogram Type {TH1F, TH2F, TH3F}"};

  // CCDB
  Configurable<std::string> cfg_ccdb_url{"cfg_ccdb_url", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cfg_ccdb_path{"cfg_ccdb_path", "Users/y/ypatley/lambda_corr_fact", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> v_cf_str = {{"h1f_LaP_cf", "h1f_LaM_cf"},
                                                    {"h2f_LaP_cf", "h2f_LaM_cf"},
                                                    {"h3f_LaP_cf", "h3f_LaM_cf"}};

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cfg_ccdb_url.value);
    ccdb->setCaching(true);

    // initialize axis specifications
    const AxisSpec axisCol(20, 0, 20, "");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(200, 1.09, 1.14, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(32, 0.3, 3.5, "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(24, -1.2, 1.2, "y");
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
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCol});

    // QA
    histos.add("Tracks/h1f_tracks_info", "# of tracks", kTH1F, {axisCol});
    histos.add("Tracks/h2f_armpod_before_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h2f_armpod_after_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h1f_lambda_pt_vs_invm", "p_{T} vs M_{#Lambda}", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("Tracks/h1f_antilambda_pt_vs_invm", "p_{T} vs M_{#bar{#Lambda}}", kTH2F, {axisV0Mass, axisV0Pt});

    // QA Lambda
    histos.add("QA/Lambda/h2f_qt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/Lambda/h1f_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA/Lambda/h1f_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA/Lambda/h1f_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA/Lambda/h1f_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA/Lambda/h1f_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});
    histos.add("QA/Lambda/h1f_V0_gctau", "V_{0} #gammac#tau", kTH1F, {axisGCTau});

    histos.add("QA/Lambda/h1f_pos_prong_pt", "Pos-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_neg_prong_pt", "Neg-Prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_pos_prong_eta", "Pos-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_neg_prong_eta", "Neg-Prong #eta-distribution", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_pos_prong_phi", "Pos-Prong #phi-distribution", kTH1F, {axisV0Phi});
    histos.add("QA/Lambda/h1f_neg_prong_phi", "Neg-Prong #phi-distribution", kTH1F, {axisV0Phi});

    histos.add("QA/Lambda/h2f_pos_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_neg_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_pos_prong_dEdx_vs_p", "TPC Signal Pos-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_neg_prong_dEdx_vs_p", "TPC Signal Neg-Prong", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisNsigma});

    // QA Anti-Lambda
    histos.addClone("QA/Lambda/", "QA/AntiLambda/");

    // MC Generated Histograms
    if (doprocessMCGen) {
      // McReco Histos
      histos.add("Tracks/h2f_tracks_pid_before_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_tracks_pid_after_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_lambda_from_sigma", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_lambda_from_cascade", "PIDs", kTH2F, {axisPID, axisV0Pt});

      // McGen Histos
      histos.add("McGen/h1f_collisions_info", "# of collisions", kTH1F, {axisCol});
      histos.add("McGen/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
      histos.add("McGen/h1f_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
    }
  }

  template <typename C>
  bool sel_collision(C const& col)
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
  bool dau_track_selection(T const& track)
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
  bool lambda_topo_kin_sel(C const& col, V const& v0, T const&)
  {
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (!dau_track_selection(postrack) || !dau_track_selection(negtrack)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 2.5);

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

    // cosine of pointing angle
    if (v0.v0cosPA() < cfg_min_V0_cosPA) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 3.5);

    // pT cut
    if (v0.pt() <= cfg_v0_pt_min || v0.pt() >= cfg_v0_pt_max) {
      return false;
    }

    // rapidity cut
    if (fabs(v0.yLambda()) >= cfg_v0_rap_max) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 4.5);

    // Armentros-Podolanski Selection
    if (cfg_armpod_flag && (fabs(v0.alpha()) < v0.qtarm() / cfg_armpod_val)) {
      return false;
    }

    // Kshort mass rejection hypothesis
    if (cfg_kshort_rej_flag && (fabs(v0.mK0Short() - MassK0Short) <= cfg_kshort_rej)) {
      return false;
    }

    return true;
  }

  template <ParticleType part, typename T>
  bool lambda_dau_tpc_pid(T const& postrack, T const& negtrack)
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

  template <typename V, typename T>
  bool sel_lambda(V const& v0, T const&, ParticleType& v0type)
  {
    // initialize daughter tracks
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // initialize selection flags
    bool lambda_flag = false, antilambda_flag = false;

    // get v0 track as lambda
    if ((fabs(v0.mLambda() - MassLambda0) < cfg_lambda_mass_window) && (lambda_dau_tpc_pid<kLambda>(postrack, negtrack))) {
      lambda_flag = true;
      v0type = kLambda;
    }

    // get v0 track as anti-lambda
    if ((fabs(v0.mAntiLambda() - MassLambda0) < cfg_lambda_mass_window) && (lambda_dau_tpc_pid<kAntiLambda>(postrack, negtrack))) {
      antilambda_flag = true;
      v0type = kAntiLambda;
    }

    if (!lambda_flag && !antilambda_flag) { // neither Lambda nor Anti-Lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), 18.5);
      return false;
    } else if (lambda_flag && antilambda_flag) { // check if the track is identified as lambda and anti-lambda both (DISCARD THIS TRACK)
      histos.fill(HIST("Tracks/h1f_tracks_info"), 19.5);
      return false;
    }

    return true;
  }

  template <typename V, typename T>
  bool sel_true_lambda(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check for lambda/anti-lambda
    if (abs(mcpart.pdgCode()) != kLambda0) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 10.5);
      return false;
    }

    // check for primary/secondary lambda
    if (cfg_rec_primary_lambda && !mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 11.5);
      return false;
    } else if (!cfg_rec_primary_lambda && mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 12.5);
      return false;
    }

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // check if the daughters have corresponding mcparticle
    if (!postrack.has_mcParticle() || !negtrack.has_mcParticle()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 13.5);
      return false;
    }

    auto mcpostrack = postrack.template mcParticle_as<aod::McParticles>();
    auto mcnegtrack = negtrack.template mcParticle_as<aod::McParticles>();

    if (abs(mcpostrack.pdgCode()) != kProton || abs(mcnegtrack.pdgCode()) != kPiPlus) { // incorrect decay channel -> return false
      histos.fill(HIST("Tracks/h1f_tracks_info"), 14.5);
      return false;
    }

    return true;
  }

  template <ParticleType part, typename C, typename V>
  float get_corr_factor(C const& col, V const& v0)
  {
    // Check for efficiency correction flag and Rec/Gen Data
    if (!cfg_corr_flag) {
      return 1.;
    }

    // Get  from CCDB
    auto ccdb_obj = ccdb->getForTimeStamp<TList>(cfg_ccdb_path.value, -1);

    // Check CCDB Object
    if (!ccdb_obj) {
      LOGF(warning, "CCDB OBJECT NOT FOUND");
      return 1.;
    }

    // get ccdb object
    TObject* obj = reinterpret_cast<TObject*>(ccdb_obj->FindObject(Form("%s", v_cf_str[cfg_cf_hist_type][part].c_str())));

    if (obj->InheritsFrom("TH1F")) {
      TH1F* hist = reinterpret_cast<TH1F*> obj->Clone();
      int pt_bin = hist->GetXaxis()->FindBin(v0.pt());
      return hist->GetBinContent(pt_bin);
    } else if (obj->InheritsFrom("TH2F")) {
      TH2F* hist = reinterpret_cast<TH2F*> obj->Clone();
      int pt_bin = hist->GetXaxis()->FindBin(v0.pt());
      int rap_bin = hist->GetYaxis()->FindBin(v0.yLambda());
      return hist->GetBinContent(pt_bin, rap_bin);
    } else if (obj->InheritsFrom("TH3F")) {
      TH3F* hist = reinterpret_cast<TH3F*> obj->Clone();
      int pt_bin = hist->GetXaxis()->FindBin(v0.pt());
      int rap_bin = hist->GetYaxis()->FindBin(v0.yLambda());
      int vz_bin = hist->GetZaxis()->FindBin(col.posZ());
      return hist->GetBinContent(pt_bin, rap_bin, vz_bin);
    } else {
      return 1.;
    }
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fill_qa_lambda(C const& col, V const& v0, T const&)
  {
    static constexpr std::string_view sub_dir[] = {"QA/Lambda/", "QA/AntiLambda/"};

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

    histos.fill(HIST(sub_dir[part]) + HIST("h2f_qt_vs_alpha"), v0.alpha(), v0.qtarm());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_dca_pos_to_PV"), v0.dcapostopv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_dca_neg_to_PV"), v0.dcanegtopv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_dca_V0_to_PV"), v0.dcav0topv());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_V0_radius"), v0.v0radius());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_V0_ctau"), ctau);
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_V0_gctau"), gctau);

    histos.fill(HIST(sub_dir[part]) + HIST("h1f_pos_prong_pt"), postrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_pos_prong_phi"), postrack.phi());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(sub_dir[part]) + HIST("h1f_neg_prong_phi"), negtrack.phi());

    histos.fill(HIST(sub_dir[part]) + HIST("h2f_pos_prong_dcaXY_vs_pt"), postrack.pt(), postrack.dcaXY());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_neg_prong_dcaXY_vs_pt"), negtrack.pt(), negtrack.dcaXY());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pr_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pr_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pi_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(sub_dir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pi_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
  }

  template <DMCType dmc, typename C, typename V, typename T>
  void fill_lambda_tables(C const& collision, V const& v0tracks, T const& tracks)
  {
    if constexpr (dmc == kMC) {
      histos.fill(HIST("Events/h1f_collisions_info"), 0.5);
      if (!collision.has_mcCollision()) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), 1.5);

    // select collision
    if (!sel_collision(collision)) {
      return;
    }

    histos.fill(HIST("Events/h1f_collisions_info"), 2.5);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());

    lambdaCollisionTable(collision.centFT0M(), collision.posX(), collision.posY(), collision.posZ());

    // initialize v0track objects
    ParticleType v0type = kLambda;
    float mass = 0., corr_fact = 1.;

    for (auto const& v0 : v0tracks) {
      // check for corresponding MCGen Particle
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), 0.5);
        if (!v0.has_mcParticle()) {
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), 1.5);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // topological and kinematic selection
      // armeteros-podolanski selection | kshort mass rejection hypothesis
      if (!lambda_topo_kin_sel(collision, v0, tracks)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), 5.5);

      // select v0 as lambda / anti-lambda
      if (!sel_lambda(v0, tracks, v0type)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), 6.5);
      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

      // we have v0 as lambda
      // do MC analysis
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cfg_sel_true_lambda && !sel_true_lambda(v0, tracks)) {
          continue;
        }
        histos.fill(HIST("Tracks/h1f_tracks_info"), 7.5);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());
      }

      // get correction factors and mass
      corr_fact = (v0type == kLambda) ? get_corr_factor<kLambda>(collision, v0) : get_corr_factor<kAntiLambda>(collision, v0);
      mass = (v0type == kLambda) ? v0.mLambda() : v0.mAntiLambda();

      // fill lambda qa
      if (v0type == kLambda) {
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fill_qa_lambda<kLambda>(collision, v0, tracks);
      } else {
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fill_qa_lambda<kAntiLambda>(collision, v0, tracks);
      }

      // Fill Lambda/AntiLambda Table
      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0.px(), v0.py(), v0.pz(),
                       v0.pt(), v0.eta(), v0.phi(), v0.yLambda(), mass,
                       v0.template posTrack_as<T>().index(), v0.template negTrack_as<T>().index(),
                       (int8_t)v0type, v0.v0cosPA(), v0.dcaV0daughters(), corr_fact);
    }
  }

  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

  void processData(Collisions::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    fill_lambda_tables<kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(lambdaCorrTableProducer, processData, "Process for DATA", true);

  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;

  void processMCReco(soa::Join<Collisions, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                     McV0Tracks const& V0s, aod::McParticles const&, TracksMC const& tracks)
  {
    fill_lambda_tables<kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(lambdaCorrTableProducer, processMCReco, "Process for MC Reconstructed", false);

  void processMCGen(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("McGen/h1f_collisions_info"), 1.5);

    // apply collision cuts
    if (mcCollision.posZ() < cfg_min_z_vtx || mcCollision.posZ() > cfg_max_z_vtx) {
      return;
    }

    histos.fill(HIST("McGen/h1f_collisions_info"), 2.5);
    histos.fill(HIST("McGen/h1f_collision_posZ"), mcCollision.posZ());
    lambdaMCGenCollisionTable(mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

    // initialize track objects
    ParticleType v0type = kLambda;

    for (auto const& mcpart : mcParticles) {

      // check for Lambda first
      if (mcpart.pdgCode() == kLambda0) {
        v0type = kLambda;
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        v0type = kAntiLambda;
      } else {
        continue;
      }

      // check for Primary Lambdas/AntiLambdas
      if (cfg_gen_primary_lambda && !mcpart.isPhysicalPrimary()) {
        continue;
      }

      // apply kinematic acceptance
      if (mcpart.pt() <= cfg_v0_pt_min || mcpart.pt() >= cfg_v0_pt_max || fabs(mcpart.y()) >= cfg_v0_rap_max) {
        continue;
      }

      // check for correct decay channel
      auto dautracks = mcpart.template daughters_as<aod::McParticles>();
      std::vector<int> daughterPDGs, daughterIDs;
      for (auto dautrack : dautracks) {
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
      }
      if (abs(daughterPDGs[0]) != kProton || abs(daughterPDGs[1]) != kPiPlus) {
        continue;
      }

      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            daughterIDs[0], daughterIDs[1], (int8_t)v0type, -999., -999., 1.);
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

  void init(InitContext const&)
  {
    // Set Density Histogram Attributes
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
    const AxisSpec axisPt(64, 0.2, 3.5, "p_{T} (GeV/#it{c})");
    const AxisSpec axisEta(24, -1.2, 1.2, "#eta");
    const AxisSpec axisCPA(100, 0.99, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisRap(cfg_nRapBins, cfg_Rap_Min, cfg_Rap_Max, "y");
    const AxisSpec axisPhi(cfg_nPhiBins, cfg_Phi_Min, cfg_Phi_Max, "#phi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "y #phi");
    const AxisSpec axisQinv(100, 0, 10, "q_{inv} (GeV/#it{c})");

    const AxisSpec axisEfPt(19, 0.2, 4.0, "p_{T}");
    const AxisSpec axisEfEta(24, -1.2, 1.2, "#eta");
    const AxisSpec axisEfRap(24, -1.2, 1.2, "y");
    const AxisSpec axisEfPosZ(10, -10., 10., "V_{Z}");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1f_collision_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/Reco/h1f_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h1i_lambda_multiplicity", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_antilambda_multiplicity", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_lambda_sdau", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_antilambda_sdau", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_lambda_totmult", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_antilambda_totmult", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});

    // InvMass, DcaDau and CosPA
    histos.add("Reco/QA_Lambda/h1f_invmass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/QA_Lambda/h1f_cospa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/QA_Lambda/h1f_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});

    histos.addClone("Reco/QA_Lambda/", "Reco/QA_AntiLambda/");

    // Efficiency Histograms
    histos.add("Reco/Efficiency/h2f_n1_pteta_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisEfPt, axisEfEta});
    histos.add("Reco/Efficiency/h2f_n1_pteta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisEfPt, axisEfEta});
    histos.add("Reco/Efficiency/h2f_n1_ptrap_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisEfPt, axisEfRap});
    histos.add("Reco/Efficiency/h2f_n1_ptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisEfPt, axisEfRap});
    histos.add("Reco/Efficiency/h3f_n1_ptetaposz_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisEfPt, axisEfEta, axisEfPosZ});
    histos.add("Reco/Efficiency/h3f_n1_ptetaposz_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisEfPt, axisEfEta, axisEfPosZ});
    histos.add("Reco/Efficiency/h3f_n1_ptrapposz_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisEfPt, axisEfRap, axisEfPosZ});
    histos.add("Reco/Efficiency/h3f_n1_ptrapposz_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisEfPt, axisEfRap, axisEfPosZ});

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

    // rho1 for R2 deta dphi histograms
    histos.add("Reco/h2d_n1_rapphi_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n1_rapphi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPhi});

    // rho1 for R2 qinv histograms
    histos.add("Reco/h2d_n1_ptrap_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisPt, axisRap});
    histos.add("Reco/h2d_n1_ptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisPt, axisRap});

    // rho2 for R2 deta dphi histograms
    histos.add("Reco/h2d_n2_rapphi_LaP_LaM", "#rho_{2}^{#Lambda - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_rapphi_LaP_LaP", "#rho_{2}^{#Lambda - #Lambda}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_rapphi_LaM_LaM", "#rho_{2}^{#bar{#Lambda} - #bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});

    // rho2 for R2 qinv histograms
    histos.add("Reco/h1d_n2_qinv_LaP_LaM", "#rho_{2}^{#Lambda-#bar{#Lambda}}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaP_LaP", "#rho_{2}^{#Lambda-#Lambda}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaM_LaM", "#rho_{2}^{#bar{#Lambda}-#bar{#Lambda}}", kTH1D, {axisQinv});

    // rho2 for qa checks
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
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1f_invmass"), x.mass());
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1f_cospa"), x.cospa());
          histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST(sub_dir[part]) + HIST("h1f_dcadau"), x.dcadau());
        }
        if (std::abs(v.mass() - MassLambda0) > std::abs(x.mass() - MassLambda0)) {
          ret_flag = false;
          break;
        }
      }
    }

    return ret_flag;
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

    float corfac1 = p1.corrfact(), corfac2 = p2.corrfact();

    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1pt2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.pt(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_eta1eta2_") + HIST(sub_dir_hist[part_pair]), p1.eta(), p2.eta(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_phi1phi2_") + HIST(sub_dir_hist[part_pair]), p1.phi(), p2.phi(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_rap1rap2_") + HIST(sub_dir_hist[part_pair]), p1.rap(), p2.rap(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1eta2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.eta(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1phi2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.phi(), corfac1 * corfac2);
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_pt1rap2_") + HIST(sub_dir_hist[part_pair]), p1.pt(), p2.rap(), corfac1 * corfac2);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n2_rapphi_") + HIST(sub_dir_hist[part_pair]), rapphix + 0.5, rapphiy + 0.5, corfac1 * corfac2);
    }

    // qinv histos
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n2_qinv_") + HIST(sub_dir_hist[part_pair]), qinv, corfac1 * corfac2);
  }

  template <ParticleType part, RecGenType rec_gen, typename C, typename T>
  void analyzeSingles(C const& col, T const& tracks)
  {
    static constexpr std::string_view sub_dir_recgen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view sub_dir_hist[] = {"LaP", "LaM"};

    int ntrk1 = 0, ntrk2 = 0, ntrk3 = 0;

    for (auto const& track : tracks) {
      ++ntrk1;
      if (!removeLambdaSharingDau<true, part, rec_gen>(track, tracks)) {
        ++ntrk2;
        continue;
      }
      ++ntrk3;

      // QA Plots
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_pt_") + HIST(sub_dir_hist[part]), track.pt(), track.corrfact());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_eta_") + HIST(sub_dir_hist[part]), track.eta(), track.corrfact());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_phi_") + HIST(sub_dir_hist[part]), track.phi(), track.corrfact());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h1d_n1_rap_") + HIST(sub_dir_hist[part]), track.rap(), track.corrfact());

      // Efficiency Calculation Plots
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("Efficiency/h2f_n1_pteta_") + HIST(sub_dir_hist[part]), track.pt(), track.eta());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("Efficiency/h2f_n1_ptrap_") + HIST(sub_dir_hist[part]), track.pt(), track.rap());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("Efficiency/h3f_n1_ptetaposz_") + HIST(sub_dir_hist[part]), track.pt(), track.eta(), col.posZ());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("Efficiency/h3f_n1_ptrapposz_") + HIST(sub_dir_hist[part]), track.pt(), track.rap(), col.posZ());

      // Rho1 for R2 Calculation
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_ptrap_") + HIST(sub_dir_hist[part]), track.pt(), track.rap(), track.corrfact());
      histos.fill(HIST(sub_dir_recgen[rec_gen]) + HIST("h2d_n1_rapphi_") + HIST(sub_dir_hist[part]), track.rap(), track.phi(), track.corrfact());
    }

    // fill multiplicity histograms
    if (ntrk1 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_lambda_totmult"), ntrk1);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_antilambda_totmult"), ntrk1);
      }
    }

    if (ntrk2 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_lambda_sdau"), ntrk2);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_antilambda_sdau"), ntrk2);
      }
    }

    if (ntrk3 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_lambda_multiplicity"), ntrk3);
      } else {
        histos.fill(HIST("Event/") + HIST(sub_dir_recgen[rec_gen]) + HIST("h1i_antilambda_multiplicity"), ntrk3);
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
    histos.fill(HIST("Event/Reco/h1f_collision_posz"), collision.posZ());
    histos.fill(HIST("Event/Reco/h1f_ft0m_mult_percentile"), collision.cent());

    auto lambda_tracks = part_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto anti_lambda_tracks = part_anti_lambda_tracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kRec>(collision, lambda_tracks);
    analyzeSingles<kAntiLambda, kRec>(collision, anti_lambda_tracks);
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

    histos.fill(HIST("Event/McGen/h1f_collision_posz"), mcgencol.posZ());

    auto lambda_mcgen_tracks = part_lambda_mcgen_tracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto antilambda_mcgen_tracks = part_antilambda_mcgen_tracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);

    analyzeSingles<kLambda, kGen>(mcgencol, lambda_mcgen_tracks);
    analyzeSingles<kAntiLambda, kGen>(mcgencol, antilambda_mcgen_tracks);
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
