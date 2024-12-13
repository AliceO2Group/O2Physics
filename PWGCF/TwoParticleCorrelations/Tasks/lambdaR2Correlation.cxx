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

/// \file lambdaR2Correlation.cxx
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
using namespace o2::constants::math;

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
DECLARE_SOA_TABLE(LambdaMcGenCollisions, "AOD", "LMCGENCOLS", o2::soa::Index<>,
                  o2::aod::mccollision::PosX,
                  o2::aod::mccollision::PosY,
                  o2::aod::mccollision::PosZ);
using LambdaMcGenCollision = LambdaMcGenCollisions::iterator;

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
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t);
DECLARE_SOA_COLUMN(V0Type, v0Type, int8_t);
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(CorrFact, corrFact, float);
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
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::CorrFact);
using LambdaTrack = LambdaTracks::iterator;

namespace lambdamcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMcGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(LambdaMcGenTracks, "AOD", "LMCGENTRACKS", o2::soa::Index<>,
                  lambdamcgentrack::LambdaMcGenCollisionId,
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
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::CorrFact);
using LambdaMcGenTrack = LambdaMcGenTracks::iterator;

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

struct LambdaCorrTableProducer {

  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::LambdaMcGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMcGenTracks> lambdaMCGenTrackTable;

  // Collisions
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "z vertex cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "z vertex cut"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.2, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 4.0, "p_{T} minimum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "min crossed rows"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 2.0, "TPC NSigma Selection Cut"};
  Configurable<double> cTrackMinDcaXY{"cTrackMinDcaXY", 0.05, "Minimum DcaXY of Daughter Tracks"};

  // V0s
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 1.0, "min DCA between V0 daughters"};
  Configurable<double> cMinDcaPosToPV{"cMinDcaPosToPV", 0.1, "Minimum V0 Positive Track DCAr cut to PV"};
  Configurable<double> cMinDcaNegToPV{"cMinDcaNegToPV", 0.1, "Minimum V0 Negative Track DCAr cut to PV"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.6, "Minimum DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.0, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 50.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Minimum ctau"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 50.0, "Maximum ctau"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.998, "Minimum V0 CosPA to PV"};
  Configurable<double> cLambdaMassWindow{"cLambdaMassWindow", 0.007, "Mass Window to select Lambda"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.005, "Reject K0Short Candidates"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", false, "K0short Mass Rej Flag"};
  Configurable<double> cArmPodCutValue{"cArmPodCutValue", 0.5, "Armentros-Podolanski Slope Parameter"};
  Configurable<bool> cArmPodCutFlag{"cArmPodCutFlag", true, "Armentros-Podolanski Cut Flag"};

  // V0s kinmatic acceptance
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.8, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 3.2, "Minimum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.6, "|rap| cut"};

  // V0s MC
  Configurable<bool> cHasMcFlag{"cHasMcFlag", true, "Has Mc Tag"};
  Configurable<bool> cSelectTrueLambda{"cSelectTrueLambda", true, "Select True Lambda"};
  Configurable<bool> cRecPrimaryLambda{"cRecPrimaryLambda", true, "Primary Lambda"};
  Configurable<bool> cGenPrimaryLambda{"cGenPrimaryLambda", true, "Primary Lambda"};

  // Efficiency Correction
  Configurable<bool> cCorrectionFlag{"cCorrectionFlag", true, "Efficiency Correction Flag"};
  Configurable<int> cCorrFactHistType{"cCorrFactHistType", 0, "Histogram Type {TH1F, TH2F, TH3F}"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cPathCCDB{"cPathCCDB", "Users/y/ypatley/lambda_corr_fact", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vCorrFactStrings = {{"h1f_LaP_cf", "h1f_LaM_cf"},
                                                            {"h2f_LaP_cf", "h2f_LaM_cf"},
                                                            {"h3f_LaP_cf", "h3f_LaM_cf"}};

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cUrlCCDB.value);
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
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

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
  bool selCollision(C const& col)
  {
    if (col.posZ() < cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }

    if (!col.sel8() && cSel8Trig) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsTriggerTVX) && cTriggerTvxSel) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoTimeFrameBorder) && cTFBorder) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoITSROFrameBorder) && cNoItsROBorder) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsVertexITSTPC) && cItsTpcVtx) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kNoSameBunchPileup) && cPileupReject) {
      return false;
    }

    if (!col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && cZVtxTimeDiff) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selDaughterTracks(T const& track)
  {
    if (track.pt() < cTrackMinPt || track.pt() > cTrackMaxPt) {
      return false;
    }

    if (std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < cMinTpcCrossedRows) {
      return false;
    }

    if (std::abs(track.dcaXY()) < cTrackMinDcaXY) {
      return false;
    }

    return true;
  }

  template <typename C, typename V, typename T>
  bool selLambdaWithTopoKinCuts(C const& col, V const& v0, T const&)
  {
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    if (!selDaughterTracks(postrack) || !selDaughterTracks(negtrack)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 2.5);

    if (v0.dcaV0daughters() > cMinV0DcaDaughters) {
      return false;
    }

    if (std::abs(v0.dcapostopv()) < cMinDcaPosToPV) {
      return false;
    }

    if (std::abs(v0.dcanegtopv()) < cMinDcaNegToPV) {
      return false;
    }

    if (v0.dcav0topv() > cMinDcaV0ToPV) {
      return false;
    }

    if ((v0.v0radius() > cMaxV0TransRadius) || (v0.v0radius() < cMinV0TransRadius)) {
      return false;
    }

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau < cMinV0CTau || ctau > cMaxV0CTau) {
      return false;
    }

    // cosine of pointing angle
    if (v0.v0cosPA() < cMinV0CosPA) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 3.5);

    // pT cut
    if (v0.pt() <= cMinV0Pt || v0.pt() >= cMaxV0Pt) {
      return false;
    }

    // rapidity cut
    if (std::abs(v0.yLambda()) >= cMaxV0Rap) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 4.5);

    // Armentros-Podolanski Selection
    if (cArmPodCutFlag && (std::abs(v0.alpha()) < v0.qtarm() / cArmPodCutValue)) {
      return false;
    }

    // Kshort mass rejection hypothesis
    if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
      return false;
    }

    return true;
  }

  template <ParticleType part, typename T>
  bool selLambdaDauWithTpcPid(T const& postrack, T const& negtrack)
  {
    bool returnFlag = false;
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

    if (std::abs(tpcNSigmaPr) < cTpcNsigmaCut && std::abs(tpcNSigmaPi) < cTpcNsigmaCut) {
      returnFlag = true;
    }

    return returnFlag;
  }

  template <typename V, typename T>
  bool selLambdaMassWindow(V const& v0, T const&, ParticleType& v0type)
  {
    // initialize daughter tracks
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // initialize selection flags
    bool lambdaFlag = false, antiLambdaFlag = false;

    // get v0 track as lambda
    if ((std::abs(v0.mLambda() - MassLambda0) < cLambdaMassWindow) && (selLambdaDauWithTpcPid<kLambda>(postrack, negtrack))) {
      lambdaFlag = true;
      v0type = kLambda;
    }

    // get v0 track as anti-lambda
    if ((std::abs(v0.mAntiLambda() - MassLambda0) < cLambdaMassWindow) && (selLambdaDauWithTpcPid<kAntiLambda>(postrack, negtrack))) {
      antiLambdaFlag = true;
      v0type = kAntiLambda;
    }

    if (!lambdaFlag && !antiLambdaFlag) { // neither Lambda nor Anti-Lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), 18.5);
      return false;
    } else if (lambdaFlag && antiLambdaFlag) { // check if the track is identified as lambda and anti-lambda both (DISCARD THIS TRACK)
      histos.fill(HIST("Tracks/h1f_tracks_info"), 19.5);
      return false;
    }

    return true;
  }

  template <typename V, typename T>
  bool selTrueMcRecLambda(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check for lambda/anti-lambda
    if (std::abs(mcpart.pdgCode()) != kLambda0) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 10.5);
      return false;
    }

    // check for primary/secondary lambda
    if (cRecPrimaryLambda && !mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), 11.5);
      return false;
    } else if (!cRecPrimaryLambda && mcpart.isPhysicalPrimary()) {
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

    if (std::abs(mcpostrack.pdgCode()) != kProton || std::abs(mcnegtrack.pdgCode()) != kPiPlus) { // incorrect decay channel -> return false
      histos.fill(HIST("Tracks/h1f_tracks_info"), 14.5);
      return false;
    }

    return true;
  }

  template <ParticleType part, typename C, typename V>
  float getCorrectionFactors(C const& col, V const& v0)
  {
    // Check for efficiency correction flag and Rec/Gen Data
    if (!cCorrectionFlag) {
      return 1.;
    }

    // Get  from CCDB
    auto ccdbObj = ccdb->getForTimeStamp<TList>(cPathCCDB.value, -1);

    // Check CCDB Object
    if (!ccdbObj) {
      LOGF(warning, "CCDB OBJECT NOT FOUND");
      return 1.;
    }

    // get ccdb object
    TObject* obj = reinterpret_cast<TObject*>(ccdbObj->FindObject(Form("%s", vCorrFactStrings[cCorrFactHistType][part].c_str())));

    if (obj->InheritsFrom("TH1F")) {
      TH1F* hist = reinterpret_cast<TH1F*>(obj->Clone());
      int ptBin = hist->GetXaxis()->FindBin(v0.pt());
      return hist->GetBinContent(ptBin);
    } else if (obj->InheritsFrom("TH2F")) {
      TH2F* hist = reinterpret_cast<TH2F*>(obj->Clone());
      int ptBin = hist->GetXaxis()->FindBin(v0.pt());
      int rapBin = hist->GetYaxis()->FindBin(v0.yLambda());
      return hist->GetBinContent(ptBin, rapBin);
    } else if (obj->InheritsFrom("TH3F")) {
      TH3F* hist = reinterpret_cast<TH3F*>(obj->Clone());
      int ptBin = hist->GetXaxis()->FindBin(v0.pt());
      int rapBin = hist->GetYaxis()->FindBin(v0.yLambda());
      int vzBin = hist->GetZaxis()->FindBin(col.posZ());
      return hist->GetBinContent(ptBin, rapBin, vzBin);
    } else {
      return 1.;
    }
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillLambdaQAHistos(C const& col, V const& v0, T const&)
  {
    static constexpr std::string_view SubDir[] = {"QA/Lambda/", "QA/AntiLambda/"};

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

    histos.fill(HIST(SubDir[part]) + HIST("h2f_qt_vs_alpha"), v0.alpha(), v0.qtarm());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_pos_to_PV"), v0.dcapostopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_neg_to_PV"), v0.dcanegtopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_to_PV"), v0.dcav0topv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_radius"), v0.v0radius());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_ctau"), ctau);
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_gctau"), gctau);

    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_pt"), postrack.pt());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_phi"), postrack.phi());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_phi"), negtrack.phi());

    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dcaXY_vs_pt"), postrack.pt(), postrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dcaXY_vs_pt"), negtrack.pt(), negtrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pr_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pr_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pi_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pi_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
  }

  template <DMCType dmc, typename C, typename V, typename T>
  void fillLambdaTables(C const& collision, V const& v0tracks, T const& tracks)
  {
    if constexpr (dmc == kMC) {
      histos.fill(HIST("Events/h1f_collisions_info"), 0.5);
      if (!collision.has_mcCollision()) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), 1.5);

    // select collision
    if (!selCollision(collision)) {
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
      if (!selLambdaWithTopoKinCuts(collision, v0, tracks)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), 5.5);

      // select v0 as lambda / anti-lambda
      if (!selLambdaMassWindow(v0, tracks, v0type)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), 6.5);
      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

      // we have v0 as lambda
      // do MC analysis
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks)) {
          continue;
        }
        histos.fill(HIST("Tracks/h1f_tracks_info"), 7.5);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());
      }

      // get correction factors and mass
      corr_fact = (v0type == kLambda) ? getCorrectionFactors<kLambda>(collision, v0) : getCorrectionFactors<kAntiLambda>(collision, v0);
      mass = (v0type == kLambda) ? v0.mLambda() : v0.mAntiLambda();

      // fill lambda qa
      if (v0type == kLambda) {
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
      } else {
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
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
    fillLambdaTables<kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCorrTableProducer, processData, "Process for DATA", true);

  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;

  void processMCReco(soa::Join<Collisions, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                     McV0Tracks const& V0s, aod::McParticles const&, TracksMC const& tracks)
  {
    fillLambdaTables<kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaCorrTableProducer, processMCReco, "Process for MC Reconstructed", false);

  void processMCGen(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("McGen/h1f_collisions_info"), 1.5);

    // apply collision cuts
    if (mcCollision.posZ() < cMinZVtx || mcCollision.posZ() > cMaxZVtx) {
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
      if (cGenPrimaryLambda && !mcpart.isPhysicalPrimary()) {
        continue;
      }

      // apply kinematic acceptance
      if (mcpart.pt() <= cMinV0Pt || mcpart.pt() >= cMaxV0Pt || std::abs(mcpart.y()) >= cMaxV0Rap) {
        continue;
      }

      // check for correct decay channel
      auto dautracks = mcpart.template daughters_as<aod::McParticles>();
      std::vector<int> daughterPDGs, daughterIDs;
      for (auto const& dautrack : dautracks) {
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
      }
      if (std::abs(daughterPDGs[0]) != kProton || std::abs(daughterPDGs[1]) != kPiPlus) {
        continue;
      }

      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            daughterIDs[0], daughterIDs[1], (int8_t)v0type, -999., -999., 1.);
    }
  }

  PROCESS_SWITCH(LambdaCorrTableProducer, processMCGen, "Process for MC Generated", false);
};

struct LambdaR2Correlation {

  // Global Configurables
  Configurable<int> cNRapBins{"cNRapBins", 12, "N Rapidity Bins"};
  Configurable<float> cMinRap{"cMinRap", -0.6, "Minimum Rapidity"};
  Configurable<float> cMaxRap{"cMaxRap", 0.6, "Maximum Rapidity"};
  Configurable<int> cNPhiBins{"cNPhiBins", 64, "N Phi Bins"};

  // remove lambda with shared daughters
  Configurable<bool> cRemoveLambdaSharingDau{"cRemoveLambdaSharingDau", true, "Flag to remove lambda"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Initialize global variables
  float nrapbins = 0.;
  float kminrap = 0.;
  float kmaxrap = 0.;
  float nphibins = 0.;
  float kminphi = 0.;
  float kmaxphi = TwoPI;
  float rapbinwidth = 0.;
  float phibinwidth = 0.;
  float q = 0., e = 0., qinv = 0.;

  void init(InitContext const&)
  {
    // Set Density Histogram Attributes
    nrapbins = static_cast<float>(cNRapBins);
    kminrap = static_cast<float>(cMinRap);
    kmaxrap = static_cast<float>(cMaxRap);
    nphibins = static_cast<float>(cNPhiBins);

    rapbinwidth = (kmaxrap - kminrap) / nrapbins;
    phibinwidth = (kmaxphi - kminphi) / nphibins;

    int knrapphibins = static_cast<int>(cNRapBins) * static_cast<int>(cNPhiBins);
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
    const AxisSpec axisRap(cNRapBins, cMinRap, cMaxRap, "y");
    const AxisSpec axisPhi(cNPhiBins, 0., TwoPI, "#phi (rad)");
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
    if (!cRemoveLambdaSharingDau) {
      return true;
    }

    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirHist[] = {"QA_Lambda/", "QA_AntiLambda/"};

    bool retFlag = true;

    for (auto const& x : vs) {
      if ((v.index() != x.index()) && (v.posTrackId() == x.posTrackId() || v.negTrackId() == x.negTrackId())) {
        if (fillHist) {
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirHist[part]) + HIST("h1f_invmass"), x.mass());
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirHist[part]) + HIST("h1f_cospa"), x.cosPA());
          histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirHist[part]) + HIST("h1f_dcadau"), x.dcaDau());
        }
        if (std::abs(v.mass() - MassLambda0) > std::abs(x.mass() - MassLambda0)) {
          retFlag = false;
          break;
        }
      }
    }

    return retFlag;
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, typename U>
  void fillPairHistos(U& p1, U& p2)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirHist[] = {"LaP_LaM", "LaP_LaP", "LaM_LaM"};

    int rapbin1 = static_cast<int>((p1.rap() - kminrap) / rapbinwidth);
    int rapbin2 = static_cast<int>((p2.rap() - kminrap) / rapbinwidth);

    int phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    int phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    float corfac1 = p1.corrFact(), corfac2 = p2.corrFact();

    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_pt1pt2_") + HIST(SubDirHist[part_pair]), p1.pt(), p2.pt(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_eta1eta2_") + HIST(SubDirHist[part_pair]), p1.eta(), p2.eta(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_phi1phi2_") + HIST(SubDirHist[part_pair]), p1.phi(), p2.phi(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_rap1rap2_") + HIST(SubDirHist[part_pair]), p1.rap(), p2.rap(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_pt1eta2_") + HIST(SubDirHist[part_pair]), p1.pt(), p2.eta(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_pt1phi2_") + HIST(SubDirHist[part_pair]), p1.pt(), p2.phi(), corfac1 * corfac2);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_pt1rap2_") + HIST(SubDirHist[part_pair]), p1.pt(), p2.rap(), corfac1 * corfac2);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_rapphi_") + HIST(SubDirHist[part_pair]), rapphix + 0.5, rapphiy + 0.5, corfac1 * corfac2);
    }

    // qinv histos
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n2_qinv_") + HIST(SubDirHist[part_pair]), qinv, corfac1 * corfac2);
  }

  template <ParticleType part, RecGenType rec_gen, typename C, typename T>
  void analyzeSingles(C const& col, T const& tracks)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirHist[] = {"LaP", "LaM"};

    int ntrk1 = 0, ntrk2 = 0, ntrk3 = 0;

    for (auto const& track : tracks) {
      ++ntrk1;
      if (!removeLambdaSharingDau<true, part, rec_gen>(track, tracks)) {
        ++ntrk2;
        continue;
      }
      ++ntrk3;

      // QA Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_pt_") + HIST(SubDirHist[part]), track.pt(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_eta_") + HIST(SubDirHist[part]), track.eta(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_phi_") + HIST(SubDirHist[part]), track.phi(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_rap_") + HIST(SubDirHist[part]), track.rap(), track.corrFact());

      // Efficiency Calculation Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h2f_n1_pteta_") + HIST(SubDirHist[part]), track.pt(), track.eta());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h2f_n1_ptrap_") + HIST(SubDirHist[part]), track.pt(), track.rap());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h3f_n1_ptetaposz_") + HIST(SubDirHist[part]), track.pt(), track.eta(), col.posZ());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h3f_n1_ptrapposz_") + HIST(SubDirHist[part]), track.pt(), track.rap(), col.posZ());

      // Rho1 for R2 Calculation
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n1_ptrap_") + HIST(SubDirHist[part]), track.pt(), track.rap(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n1_rapphi_") + HIST(SubDirHist[part]), track.rap(), track.phi(), track.corrFact());
    }

    // fill multiplicity histograms
    if (ntrk1 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_lambda_totmult"), ntrk1);
      } else {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_antilambda_totmult"), ntrk1);
      }
    }

    if (ntrk2 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_lambda_sdau"), ntrk2);
      } else {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_antilambda_sdau"), ntrk2);
      }
    }

    if (ntrk3 != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_lambda_multiplicity"), ntrk3);
      } else {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_antilambda_multiplicity"), ntrk3);
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

  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaTracks = aod::LambdaTracks;

  SliceCache cache;

  Partition<LambdaTracks> partLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda);
  Partition<LambdaTracks> partAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda);

  void processDataReco(LambdaCollisions::iterator const& collision, LambdaTracks const&)
  {
    histos.fill(HIST("Event/Reco/h1f_collision_posz"), collision.posZ());
    histos.fill(HIST("Event/Reco/h1f_ft0m_mult_percentile"), collision.cent());

    auto lambdaTracks = partLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto antiLambdaTracks = partAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kRec>(collision, lambdaTracks);
    analyzeSingles<kAntiLambda, kRec>(collision, antiLambdaTracks);
    analyzePairs<kLambdaAntiLambda, kRec, false>(lambdaTracks, antiLambdaTracks);
    analyzePairs<kLambdaLambda, kRec, true>(lambdaTracks, lambdaTracks);
    analyzePairs<kAntiLambdaAntiLambda, kRec, true>(antiLambdaTracks, antiLambdaTracks);
  }

  PROCESS_SWITCH(LambdaR2Correlation, processDataReco, "Process for Data and MCReco", true);

  using LambdaMcGenCollisions = aod::LambdaMcGenCollisions;
  using LambdaMcGenTracks = aod::LambdaMcGenTracks;

  SliceCache cachemc;
  Partition<LambdaMcGenTracks> partLambdaMcGenTracks = aod::lambdatrack::v0Type == (int8_t)kLambda;
  Partition<LambdaMcGenTracks> partAntiLambdaMcGenTracks = aod::lambdatrack::v0Type == (int8_t)kAntiLambda;

  void processMCGen(LambdaMcGenCollisions::iterator const& mcgencol, LambdaMcGenTracks const&)
  {

    histos.fill(HIST("Event/McGen/h1f_collision_posz"), mcgencol.posZ());

    auto lambdaMcGenTracks = partLambdaMcGenTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto antiLambdaMcGenTracks = partAntiLambdaMcGenTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);

    analyzeSingles<kLambda, kGen>(mcgencol, lambdaMcGenTracks);
    analyzeSingles<kAntiLambda, kGen>(mcgencol, antiLambdaMcGenTracks);
    analyzePairs<kLambdaAntiLambda, kGen, false>(lambdaMcGenTracks, antiLambdaMcGenTracks);
    analyzePairs<kLambdaLambda, kGen, true>(lambdaMcGenTracks, lambdaMcGenTracks);
    analyzePairs<kAntiLambdaAntiLambda, kGen, true>(antiLambdaMcGenTracks, antiLambdaMcGenTracks);
  }

  PROCESS_SWITCH(LambdaR2Correlation, processMCGen, "Process for MC Generated", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LambdaCorrTableProducer>(cfgc),
    adaptAnalysisTask<LambdaR2Correlation>(cfgc)};
}
