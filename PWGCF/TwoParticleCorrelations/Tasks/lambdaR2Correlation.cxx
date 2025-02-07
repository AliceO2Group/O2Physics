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
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"
#include "CCDB/BasicCCDBManager.h"
#include "TPDGCode.h"
#include "TRandom.h"

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
                  lambdacollision::Cent,
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

namespace lambdatrackext
{
DECLARE_SOA_COLUMN(LambdaSharingDaughter, lambdaSharingDaughter, bool);
DECLARE_SOA_COLUMN(LambdaSharingDauIds, lambdaSharingDauIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(TrueLambdaFlag, trueLambdaFlag, bool);
} // namespace lambdatrackext
DECLARE_SOA_TABLE(LambdaTracksExt, "AOD", "LAMBDATRACKSEXT",
                  lambdatrackext::LambdaSharingDaughter,
                  lambdatrackext::LambdaSharingDauIds,
                  lambdatrackext::TrueLambdaFlag);

using LambdaTrackExt = LambdaTracksExt::iterator;

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

enum CollisionLabels {
  kTotColBeforeHasMcCollision = 1,
  kTotCol,
  kPassSelCol
};

enum TrackLabels {
  kTracksBeforeHasMcParticle = 1,
  kAllV0Tracks,
  kPassV0DauTrackSel,
  kPassV0KinCuts,
  kPassV0TopoSel,
  kNotLambdaNotAntiLambda,
  kV0AsLambdaAntiLambda,
  kPassV0MassWinCuts,
  kNotPrimaryLambda,
  kNotSecondaryLambda,
  kLambdaDauNotMcParticle,
  kLambdaNotPrPiMinus,
  kAntiLambdaNotAntiPrPiPlus,
  kPassTrueLambdaSel,
  kEffCorrPt,
  kEffCorrPtRap,
  kEffCorrPtRapVz,
  kNoEffCorr,
  kGenTotAccLambda,
  kGenLambdaNoDau,
  kGenLambdaToPrPi
};

enum RunType {
  kRun3 = 0,
  kRun2
};

enum ParticleType {
  kLambda = 0,
  kAntiLambda
};

enum ParticlePairType {
  kLambdaAntiLambda = 0,
  kLambdaLambda,
  kAntiLambdaAntiLambda
};

enum ShareDauLambda {
  kUniqueLambda = 0,
  kLambdaShareDau
};

enum RecGenType {
  kRec = 0,
  kGen
};

enum DMCType {
  kData = 0,
  kMC
};

struct LambdaTableProducer {

  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::LambdaMcGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMcGenTracks> lambdaMCGenTrackTable;

  // Collisions
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "z vertex cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "z vertex cut"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cInt7Trig{"cInt7Trig", false, "kINT7 MB Trigger"};
  Configurable<bool> cSel7Trig{"cSel7Trig", false, "Sel7 (V0A + V0C) Selection Run2"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", false, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", false, "Good ITS Layers All"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.2, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 999.0, "p_{T} minimum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "min crossed rows"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 2.0, "TPC NSigma Selection Cut"};
  Configurable<double> cTrackMinDcaXY{"cTrackMinDcaXY", 0.05, "Minimum DcaXY of Daughter Tracks"};
  Configurable<bool> cIsGlobalTrackWoDca{"cIsGlobalTrackWoDca", true, "Check for Global Track"};

  // V0s
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Minimum DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1., "Maximum DCA between V0 daughters"};
  Configurable<double> cMinDcaPosToPV{"cMinDcaPosToPV", 0.05, "Minimum V0 Positive Track DCAr cut to PV"};
  Configurable<double> cMinDcaNegToPV{"cMinDcaNegToPV", 0.05, "Minimum V0 Negative Track DCAr cut to PV"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Minimum DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 999.0, "Maximum DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Minimum ctau"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Maximum ctau"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cLambdaMassWindow{"cLambdaMassWindow", 0.005, "Mass Window to select Lambda"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.01, "Reject K0Short Candidates"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0short Mass Rej Flag"};
  Configurable<double> cArmPodCutValue{"cArmPodCutValue", 0.5, "Armentros-Podolanski Slope Parameter"};
  Configurable<bool> cArmPodCutFlag{"cArmPodCutFlag", false, "Armentros-Podolanski Cut Flag"};

  // V0s kinmatic acceptance
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.8, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 2.8, "Minimum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.6, "|rap| cut"};
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Do Eta Analysis"};

  // V0s MC
  Configurable<bool> cHasMcFlag{"cHasMcFlag", true, "Has Mc Tag"};
  Configurable<bool> cSelectTrueLambda{"cSelectTrueLambda", false, "Select True Lambda"};
  Configurable<bool> cRecPrimaryLambda{"cRecPrimaryLambda", false, "Primary Reconstructed Lambda"};
  Configurable<bool> cRecSecondaryLambda{"cRecSecondaryLambda", false, "Secondary Reconstructed Lambda"};
  Configurable<bool> cGenPrimaryLambda{"cGenPrimaryLambda", true, "Primary Generated Lambda"};
  Configurable<bool> cGenDecayChannel{"cGenDecayChannel", true, "Gen Level Decay Channel Flag"};

  // Efficiency Correction
  Configurable<bool> cCorrectionFlag{"cCorrectionFlag", false, "Efficiency Correction Flag"};
  Configurable<int> cCorrFactHist{"cCorrFactHist", 0, "Correction Factor Histogram"};
  Configurable<bool> cDoEtaCorr{"cDoEtaCorr", false, "Do Eta Corr"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cPathCCDB{"cPathCCDB", "Users/y/ypatley/lambda_corr_fact", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vCorrFactStrings = {{"hEffVsPtLambda", "hEffVsPtAntiLambda"},
                                                            {"hEffVsPtYLambda", "hEffVsPtYAntiLambda"},
                                                            {"hEffVsPtEtaLambda", "hEffVsPtEtaAntiLambda"},
                                                            {"hEffVsPtYVzLambda", "hEffVsPtYVzAntiLambda"},
                                                            {"hEffVsPtEtaVzLambda", "hEffVsPtEtaVzAntiLambda"}};

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cUrlCCDB.value);
    ccdb->setCaching(true);

    // initialize axis specifications
    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(25, 0.5, 25.5, "");
    const AxisSpec axisCent(105, 0, 105, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(200, 1.09, 1.14, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(38, 0.2, 4.0, "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(24, -1.2, 1.2, "y");
    const AxisSpec axisV0Eta(24, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

    const AxisSpec axisRadius(2000, 0, 200, "r(cm)");
    const AxisSpec axisCosPA(500, 0.995, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(1000, 0., 10., "dca (cm)");
    const AxisSpec axisDcaProngPV(5000, -50., 50., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisCTau(2000, 0, 200, "c#tau (cm)");
    const AxisSpec axisGCTau(2000, 0, 200, "#gammac#tau (cm)");
    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    const AxisSpec axisTrackPt(40, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackDCA(200, -1, 1, "dca_{XY} (cm)");
    const AxisSpec axisMomPID(80, 0, 4, "p (GeV/#it{c})");
    const AxisSpec axisNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    // Create Histograms.
    // Event histograms
    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCols});
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});

    // QA
    histos.add("Tracks/h1f_tracks_info", "# of tracks", kTH1F, {axisTrks});
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
      histos.add("Tracks/h2f_lambda_from_omega", "PIDs", kTH2F, {axisPID, axisV0Pt});

      // McGen Histos
      histos.add("McGen/h1f_collisions_info", "# of collisions", kTH1F, {axisCols});
      histos.add("McGen/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
      histos.add("McGen/h1f_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});

      // set bin lables specific to MC
      histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotColBeforeHasMcCollision, "kTotColBeforeHasMcCollision");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kTracksBeforeHasMcParticle, "kTracksBeforeHasMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotPrimaryLambda, "kNotPrimaryLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotSecondaryLambda, "kNotSecondaryLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kLambdaDauNotMcParticle, "kLambdaDauNotMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kLambdaNotPrPiMinus, "kLambdaNotPrPiMinus");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAntiLambdaNotAntiPrPiPlus, "kAntiLambdaNotAntiPrPiPlus");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassTrueLambdaSel, "kPassTrueLambdaSel");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenTotAccLambda, "kGenTotAccLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenLambdaNoDau, "kGenLambdaNoDau");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kGenLambdaToPrPi, "kGenLambdaToPrPi");
    }

    // set bin labels
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllV0Tracks, "kAllV0Tracks");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0DauTrackSel, "kPassV0DauTrackSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0KinCuts, "kPassV0KinCuts");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0TopoSel, "kPassV0TopoSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotLambdaNotAntiLambda, "kNotLambdaNotAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0AsLambdaAntiLambda, "kV0AsLambdaAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0MassWinCuts, "kPassV0MassWinCuts");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPt, "kEffCorrPt");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtRap, "kEffCorrPtRap");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtRapVz, "kEffCorrPtRapVz");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNoEffCorr, "kNoEffCorr");
  }

  template <RunType run, typename C>
  bool selCollision(C const& col)
  {
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }

    if constexpr (run == kRun3) {
      if (cSel8Trig && !col.sel8()) {
        return false;
      }
    } else {
      if (cInt7Trig && !col.alias_bit(kINT7)) {
        return false;
      }

      if (cSel7Trig && !col.sel7()) {
        return false;
      }
    }

    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }

    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }

    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }

    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }

    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }

    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool selDaughterTracks(T const& track)
  {
    if (track.pt() <= cTrackMinPt || track.pt() >= cTrackMaxPt) {
      return false;
    }

    if (std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    if (track.tpcNClsCrossedRows() <= cMinTpcCrossedRows) {
      return false;
    }

    if (std::abs(track.dcaXY()) <= cTrackMinDcaXY) {
      return false;
    }

    if (cIsGlobalTrackWoDca && !track.isGlobalTrackWoDCA()) {
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

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0DauTrackSel);

    // pT cut
    if (v0.pt() <= cMinV0Pt || v0.pt() >= cMaxV0Pt) {
      return false;
    }

    // rapidity cut
    if (std::abs(v0.yLambda()) >= cMaxV0Rap) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0KinCuts);

    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters) {
      return false;
    }

    if (std::abs(v0.dcapostopv()) < cMinDcaPosToPV) {
      return false;
    }

    if (std::abs(v0.dcanegtopv()) < cMinDcaNegToPV) {
      return false;
    }

    if (v0.dcav0topv() <= cMinDcaV0ToPV || v0.dcav0topv() >= cMaxDcaV0ToPV) {
      return false;
    }

    if (v0.v0radius() <= cMinV0TransRadius || v0.v0radius() >= cMaxV0TransRadius) {
      return false;
    }

    // ctau
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau <= cMinV0CTau || ctau >= cMaxV0CTau) {
      return false;
    }

    // cosine of pointing angle
    if (v0.v0cosPA() <= cMinV0CosPA) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0TopoSel);

    // all selection criterion passed (Return True)
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
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaNotAntiLambda);
      return false;
    } else if (lambdaFlag && antiLambdaFlag) { // check if the track is identified as lambda and anti-lambda both (DISCARD THIS TRACK)
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0AsLambdaAntiLambda);
      return false;
    }

    // Armentros-Podolanski Selection
    if (cArmPodCutFlag && (std::abs(v0.alpha()) < v0.qtarm() / cArmPodCutValue)) {
      return false;
    }

    // Kshort mass rejection hypothesis
    if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
      return false;
    }

    if (lambdaFlag || antiLambdaFlag) {
      return true;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), 23.5);

    return false;
  }

  template <typename V, typename T>
  bool selTrueMcRecLambda(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    if (std::abs(mcpart.pdgCode()) != kLambda0) {
      return false;
    }

    // check for primary/secondary lambda
    if (cRecPrimaryLambda && !mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotPrimaryLambda);
      return false;
    } else if (cRecSecondaryLambda && mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotSecondaryLambda);
      return false;
    }

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // check if the daughters have corresponding mcparticle
    if (!postrack.has_mcParticle() || !negtrack.has_mcParticle()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kLambdaDauNotMcParticle);
      return false;
    }

    auto mcpostrack = postrack.template mcParticle_as<aod::McParticles>();
    auto mcnegtrack = negtrack.template mcParticle_as<aod::McParticles>();

    if (mcpart.pdgCode() == kLambda0) {
      if (mcpostrack.pdgCode() != kProton || mcnegtrack.pdgCode() != kPiMinus) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kLambdaNotPrPiMinus);
        return false;
      }
    } else if (mcpart.pdgCode() == kLambda0Bar) {
      if (mcpostrack.pdgCode() != kPiPlus || mcnegtrack.pdgCode() != kProtonBar) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kAntiLambdaNotAntiPrPiPlus);
        return false;
      }
    }

    // get information about secondary lambdas
    if (cRecSecondaryLambda) {
      auto lambdaMothers = mcpart.template mothers_as<aod::McParticles>();
      if (std::abs(lambdaMothers[0].pdgCode()) == 3112 || std::abs(lambdaMothers[0].pdgCode()) == 3212 || std::abs(lambdaMothers[0].pdgCode()) == 3222) {
        histos.fill(HIST("Tracks/h2f_lambda_from_sigma"), mcpart.pdgCode(), mcpart.pt());
      } else if (std::abs(lambdaMothers[0].pdgCode()) == 3312 || std::abs(lambdaMothers[0].pdgCode()) == 3322) {
        histos.fill(HIST("Tracks/h2f_lambda_from_cascade"), mcpart.pdgCode(), mcpart.pt());
      } else if (std::abs(lambdaMothers[0].pdgCode()) == 3334) {
        histos.fill(HIST("Tracks/h2f_lambda_from_omega"), mcpart.pdgCode(), mcpart.pt());
      }
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
    TObject* obj = reinterpret_cast<TObject*>(ccdbObj->FindObject(Form("%s", vCorrFactStrings[cCorrFactHist][part].c_str())));
    TH1F* hist = reinterpret_cast<TH1F*>(obj->Clone());
    float retVal = 0.;
    float rap = (cDoEtaCorr) ? v0.eta() : v0.yLambda();

    if (std::string(obj->ClassName()) == "TH1F") {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPt);
      retVal = hist->GetBinContent(hist->FindBin(v0.pt()));
    } else if (std::string(obj->ClassName()) == "TH2F") {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtRap);
      retVal = hist->GetBinContent(hist->FindBin(v0.pt(), rap));
    } else if (std::string(obj->ClassName()) == "TH3F") {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtRapVz);
      retVal = hist->GetBinContent(hist->FindBin(v0.pt(), rap, col.posZ()));
    } else {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNoEffCorr);
      LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
      retVal = 1.;
    }

    delete hist;
    return retVal;
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

  template <RunType run, DMCType dmc, typename C, typename V, typename T>
  void fillLambdaTables(C const& collision, V const& v0tracks, T const& tracks)
  {
    if constexpr (dmc == kMC) {
      histos.fill(HIST("Events/h1f_collisions_info"), kTotColBeforeHasMcCollision);
      if (!collision.has_mcCollision()) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    // select collision
    if (!selCollision<run>(collision)) {
      return;
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());

    float cent = 0.;

    if constexpr (run == kRun3) {
      cent = collision.centFT0M();
    } else {
      cent = collision.centRun2V0M();
    }

    lambdaCollisionTable(cent, collision.posX(), collision.posY(), collision.posZ());

    // initialize v0track objects
    ParticleType v0type = kLambda;
    float mass = 0., corr_fact = 1.;

    for (auto const& v0 : v0tracks) {
      // check for corresponding MCGen Particle
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kTracksBeforeHasMcParticle);
        if (!v0.has_mcParticle() || !v0.template posTrack_as<T>().has_mcParticle() || !v0.template negTrack_as<T>().has_mcParticle()) {
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // topological and kinematic selection
      if (!selLambdaWithTopoKinCuts(collision, v0, tracks)) {
        continue;
      }

      // select v0 as lambda / anti-lambda
      // armeteros-podolanski selection | kshort mass rejection hypothesis
      if (!selLambdaMassWindow(v0, tracks, v0type)) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0MassWinCuts);

      // we have v0 as lambda
      // do MC analysis
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks)) {
          continue;
        }
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPassTrueLambdaSel);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());
      }

      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

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

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;

  void processDataRun3(CollisionsRun3::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    fillLambdaTables<kRun3, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processDataRun3, "Process for Run3 DATA", true);

  void processDataRun2(CollisionsRun2::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    fillLambdaTables<kRun2, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processDataRun2, "Process for Run2 DATA", false);

  void processMCRecoRun3(soa::Join<CollisionsRun3, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, aod::McParticles const&, TracksMC const& tracks)
  {
    fillLambdaTables<kRun3, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRecoRun3, "Process for Run3 MC Reconstructed", false);

  void processMCRecoRun2(soa::Join<CollisionsRun2, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, aod::McParticles const&, TracksMC const& tracks)
  {
    fillLambdaTables<kRun2, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRecoRun2, "Process for Run2 MC Reconstructed", false);

  void processMCGen(soa::Join<aod::McCollisions, aod::McCentFT0Ms>::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("McGen/h1f_collisions_info"), 1.5);

    // apply collision cuts
    if (mcCollision.posZ() <= cMinZVtx || mcCollision.posZ() >= cMaxZVtx) {
      return;
    }

    histos.fill(HIST("McGen/h1f_collisions_info"), 2.5);
    histos.fill(HIST("McGen/h1f_collision_posZ"), mcCollision.posZ());
    lambdaMCGenCollisionTable(mcCollision.centFT0M(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

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

      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenTotAccLambda);

      // get daughter track info and check for decay channel flag
      if (!mcpart.has_daughters()) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaNoDau);
        continue;
      }
      auto dautracks = mcpart.template daughters_as<aod::McParticles>();
      std::vector<int> daughterPDGs, daughterIDs;
      bool dauKinCutFlag = true;
      for (auto const& dautrack : dautracks) {
        // check kinematic selection on daughters as well
        if (dautrack.pt() <= cTrackMinPt || dautrack.pt() >= cTrackMaxPt || std::abs(dautrack.eta()) >= cTrackEtaCut) {
          dauKinCutFlag = false;
        }
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
      }
      if (!dauKinCutFlag) { // check daughter acceptance
        continue;
      }
      if (cGenDecayChannel && (std::abs(daughterPDGs[0]) != kProton || std::abs(daughterPDGs[1]) != kPiPlus)) { // check decay channel
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaToPrPi);

      if (v0type == kLambda) {
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), mcpart.pdgCode());
      } else {
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), mcpart.pdgCode());
      }

      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            daughterIDs[0], daughterIDs[1], (int8_t)v0type, -999., -999., 1.);
    }
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCGen, "Process for MC Generated", false);
};

struct LambdaTracksExtProducer {

  Produces<aod::LambdaTracksExt> lambdaTrackExtTable;

  // Configurables
  Configurable<bool> cAcceptAllLambda{"cAcceptAllLambda", false, "Accept all Lambda"};
  Configurable<bool> cRejAllLambdaShaDau{"cRejAllLambdaShaDau", true, "Reject all Lambda sharing daughters"};
  Configurable<bool> cSelLambdaMassPdg{"cSelLambdaMassPdg", false, "Select Lambda closest to Pdg Mass"};
  Configurable<bool> cSelLambdaTScore{"cSelLambdaTScore", false, "Select Lambda based on t-score"};
  Configurable<float> cA{"cA", 0.6, "a * |lambdaMass - lambdaPdgMass|"};
  Configurable<float> cB{"cB", 0.6, "b * DcaPrPi"};
  Configurable<float> cC{"cC", 0.6, "c * Cos(theta_{PA})"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // Axis Specifications
    const AxisSpec axisMult(10, 0, 10);
    const AxisSpec axisMass(100, 1.06, 1.16, "Inv Mass (GeV/#it{c}^{2})");
    const AxisSpec axisCPA(100, 0.995, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");

    // Histograms Booking
    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_lambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_antilambda_mult", "Multiplicity", kTH1I, {axisMult});

    // InvMass, DcaDau and CosPA
    histos.add("Reco/h1f_lambda_invmass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/h1f_lambda_cospa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/h1f_lambda_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});
    histos.add("Reco/h1f_antilambda_invmass", "M_{p#pi}", kTH1F, {axisMass});
    histos.add("Reco/h1f_antilambda_cospa", "cos(#theta_{PA})", kTH1F, {axisCPA});
    histos.add("Reco/h1f_antilambda_dcadau", "DCA_{p#pi} at V0 Decay Vertex", kTH1F, {axisDcaDau});

    histos.addClone("Reco/", "SharingDau/");
  }

  template <ShareDauLambda sd, typename T>
  void fillHistos(T const& track)
  {
    static constexpr std::string_view SubDir[] = {"Reco/", "SharingDau/"};

    if (track.v0Type() == kLambda) {
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_invmass"), track.mass());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_dcadau"), track.dcaDau());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_lambda_cospa"), track.cosPA());
    } else {
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_invmass"), track.mass());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_dcadau"), track.dcaDau());
      histos.fill(HIST(SubDir[sd]) + HIST("h1f_antilambda_cospa"), track.cosPA());
    }
  }

  void process(aod::LambdaCollisions::iterator const&, aod::LambdaTracks const& tracks)
  {

    int nTotLambda = 0, nTotAntiLambda = 0, nSelLambda = 0, nSelAntiLambda = 0;

    for (auto const& lambda : tracks) {
      bool lambdaMinDeltaMassFlag = true, lambdaMinTScoreFlag = true;
      bool lambdaSharingDauFlag = false, trueLambdaFlag = false;
      std::vector<int64_t> vSharedDauLambdaIndex;
      float tLambda = 0., tTrack = 0.;

      if (lambda.v0Type() == kLambda) {
        ++nTotLambda;
      } else if (lambda.v0Type() == kAntiLambda) {
        ++nTotAntiLambda;
      }

      tLambda = (cA * std::abs(lambda.mass() - MassLambda0)) + (cB * lambda.dcaDau()) + (cC * std::abs(lambda.cosPA() - 1.));

      for (auto const& track : tracks) {
        // check lambda index (don't analyze same lambda track !!!)
        if (lambda.index() == track.index()) {
          continue;
        }

        // check only lambda-lambda || antilambda-antilambda
        if (lambda.v0Type() != track.v0Type()) {
          continue;
        }

        // check if lambda shares daughters with any other track
        if (lambda.posTrackId() == track.posTrackId() || lambda.negTrackId() == track.negTrackId()) {
          vSharedDauLambdaIndex.push_back(track.index());
          lambdaSharingDauFlag = true;

          // decision based on mass closest to PdgMass of Lambda
          if (std::abs(lambda.mass() - MassLambda0) > std::abs(track.mass() - MassLambda0)) {
            lambdaMinDeltaMassFlag = false;
          }

          // decisions based on t-score
          tTrack = (cA * std::abs(track.mass() - MassLambda0)) + (cB * track.dcaDau()) + (cC * std::abs(track.cosPA() - 1.));
          if (tLambda > tTrack) {
            lambdaMinTScoreFlag = false;
          }
        }
      }

      // fill QA histograms
      if (lambdaSharingDauFlag) {
        fillHistos<kLambdaShareDau>(lambda);
      } else {
        fillHistos<kUniqueLambda>(lambda);
      }

      if (cAcceptAllLambda) { // Accept all lambda
        trueLambdaFlag = true;
      } else if (cRejAllLambdaShaDau && !lambdaSharingDauFlag) { // Reject all lambda sharing daughter
        trueLambdaFlag = true;
      } else if (cSelLambdaMassPdg && lambdaMinDeltaMassFlag) { // Select lambda closest to pdg mass
        trueLambdaFlag = true;
      } else if (cSelLambdaTScore && lambdaMinTScoreFlag) { // Select lambda based on t-score
        trueLambdaFlag = true;
      }

      // Multiplicity of selected lambda
      if (trueLambdaFlag) {
        if (lambda.v0Type() == kLambda) {
          ++nSelLambda;
        } else if (lambda.v0Type() == kAntiLambda) {
          ++nSelAntiLambda;
        }
      }

      // fill LambdaTrackExt table
      lambdaTrackExtTable(lambdaSharingDauFlag, vSharedDauLambdaIndex, trueLambdaFlag);
    }

    // fill multiplicity histograms
    if (nTotLambda != 0) {
      histos.fill(HIST("h1i_totlambda_mult"), nTotLambda);
    }

    if (nTotAntiLambda != 0) {
      histos.fill(HIST("h1i_totantilambda_mult"), nTotAntiLambda);
    }

    if (nSelLambda != 0) {
      histos.fill(HIST("h1i_lambda_mult"), nSelLambda);
    }

    if (nSelAntiLambda != 0) {
      histos.fill(HIST("h1i_antilambda_mult"), nSelAntiLambda);
    }
  }
};

struct LambdaR2Correlation {

  // Global Configurables
  Configurable<int> cNPtBins{"cNPtBins", 10, "N pT Bins"};
  Configurable<float> cMinPt{"cMinPt", 0.8, "pT Min"};
  Configurable<float> cMaxPt{"cMaxPt", 2.8, "pT Max"};
  Configurable<int> cNRapBins{"cNRapBins", 12, "N Rapidity Bins"};
  Configurable<float> cMinRap{"cMinRap", -0.6, "Minimum Rapidity"};
  Configurable<float> cMaxRap{"cMaxRap", 0.6, "Maximum Rapidity"};
  Configurable<int> cNPhiBins{"cNPhiBins", 36, "N Phi Bins"};

  // Eta/Rap Analysis
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Eta/Rap Analysis Flag"};

  // Rotation Angle Min/Max
  Configurable<float> cRotAngleMin{"cRotAngleMin", -0.12, "Rotation Angle Minimum"};
  Configurable<float> cRotAngleMax{"cRotAngleMax", 0.12, "Rotation Angle Minimum"};

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
  TRandom* ran = new TRandom();

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
    const AxisSpec axisPt(cNPtBins, cMinPt, cMaxPt, "p_{T} (GeV/#it{c})");
    const AxisSpec axisEta(cNRapBins, cMinRap, cMaxRap, "#eta");
    const AxisSpec axisRap(cNRapBins, cMinRap, cMaxRap, "y");
    const AxisSpec axisPhi(cNPhiBins, 0., TwoPI, "#phi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "y #phi");
    const AxisSpec axisQinv(100, 0, 10, "q_{inv} (GeV/#it{c})");

    const AxisSpec axisEfPt(cNPtBins, cMinPt, cMaxPt, "p_{T}");
    const AxisSpec axisEfEta(cNRapBins, cMinRap, cMaxRap, "#eta");
    const AxisSpec axisEfRap(cNRapBins, cMinRap, cMaxRap, "y");
    const AxisSpec axisEfPosZ(10, -10., 10., "V_{Z}");
    const AxisSpec axisEfCent(10, 0, 100, "FT0M(%)");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1f_collision_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/Reco/h1f_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h1i_lambda_mult", "#Lambda - Multiplicity", kTH1I, {axisMult});
    histos.add("Event/Reco/h1i_antilambda_mult", "#bar{#Lambda} - Multiplicity", kTH1I, {axisMult});

    // Efficiency Histograms
    // Single Particle Efficiencies
    histos.add("Reco/Efficiency/h1f_n1_pt_LaP", "#rho_{1}^{#Lambda}", kTH1F, {axisEfPt});
    histos.add("Reco/Efficiency/h1f_n1_pt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1F, {axisEfPt});
    histos.add("Reco/Efficiency/h2f_n1_pteta_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisEfPt, axisEfEta});
    histos.add("Reco/Efficiency/h2f_n1_pteta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisEfPt, axisEfEta});
    histos.add("Reco/Efficiency/h2f_n1_ptrap_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisEfPt, axisEfRap});
    histos.add("Reco/Efficiency/h2f_n1_ptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisEfPt, axisEfRap});
    histos.add("Reco/Efficiency/h4f_n1_ptetavzmult_LaP", "#rho_{1}^{#Lambda}", kTHnSparseF, {axisEfPt, axisEfEta, axisEfPosZ, axisEfCent});
    histos.add("Reco/Efficiency/h4f_n1_ptetavzmult_LaM", "#rho_{1}^{#bar{#Lambda}}", kTHnSparseF, {axisEfPt, axisEfEta, axisEfPosZ, axisEfCent});
    histos.add("Reco/Efficiency/h4f_n1_ptrapvzmult_LaP", "#rho_{1}^{#Lambda}", kTHnSparseF, {axisEfPt, axisEfRap, axisEfPosZ, axisEfCent});
    histos.add("Reco/Efficiency/h4f_n1_ptrapvzmult_LaM", "#rho_{1}^{#bar{#Lambda}}", kTHnSparseF, {axisEfPt, axisEfRap, axisEfPosZ, axisEfCent});

    // Single and Two Particle Densities
    // 1D Histograms
    histos.add("Reco/h1d_n1_mass_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisMass});
    histos.add("Reco/h1d_n1_mass_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisMass});
    histos.add("Reco/h1d_n1_pt_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisPt});
    histos.add("Reco/h1d_n1_pt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisPt});
    histos.add("Reco/h1d_n1_eta_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_eta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisEta});
    histos.add("Reco/h1d_n1_rap_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisRap});
    histos.add("Reco/h1d_n1_rap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisRap});
    histos.add("Reco/h1d_n1_phi_LaP", "#rho_{1}^{#Lambda}", kTH1D, {axisPhi});
    histos.add("Reco/h1d_n1_phi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH1D, {axisPhi});

    // rho1 for R2 RapPhi histograms
    histos.add("Reco/h2d_n1_rapphi_LaP", "#rho_{1}^{#Lambda}", kTH2D, {axisRap, axisPhi});
    histos.add("Reco/h2d_n1_rapphi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2D, {axisRap, axisPhi});

    // rho2 for R2 Rap1Phi1Rap2Phi2 histograms
    histos.add("Reco/h2d_n2_rapphi_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_rapphi_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2D, {axisRapPhi, axisRapPhi});
    histos.add("Reco/h2d_n2_rapphi_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2D, {axisRapPhi, axisRapPhi});

    // rho2 for R2 Qinv histograms
    histos.add("Reco/h1d_n2_qinv_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n2_qinv_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n1n1_qinv_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n1n1_qinv_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH1D, {axisQinv});
    histos.add("Reco/h1d_n1n1_qinv_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH1D, {axisQinv});

    // MCGen
    if (doprocessMCGen) {
      histos.addClone("Event/Reco/", "Event/McGen/");
      histos.addClone("Reco/", "McGen/");
    }
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, typename U>
  void fillPairHistos(U& p1, U& p2)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirHist[] = {"LaP_LaM", "LaP_LaP", "LaM_LaM"};

    float rap1 = (cDoEtaAnalysis) ? p1.eta() : p1.rap();
    float rap2 = (cDoEtaAnalysis) ? p2.eta() : p2.rap();

    int rapbin1 = static_cast<int>((rap1 - kminrap) / rapbinwidth);
    int rapbin2 = static_cast<int>((rap2 - kminrap) / rapbinwidth);

    int phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    int phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    float corfac = p1.corrFact() * p2.corrFact();

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n2_rapphi_") + HIST(SubDirHist[part_pair]), rapphix + 0.5, rapphiy + 0.5, corfac);
    }

    // qinv histograms
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n2_qinv_") + HIST(SubDirHist[part_pair]), qinv, corfac);

    // Rotate momentum vector about z-axis to get N1N1_Qinv Histograms
    float ranPhi = ran->Uniform(cRotAngleMin, cRotAngleMax);
    float p2x = p2.pt() * std::cos(RecoDecay::constrainAngle((p2.phi() + ranPhi), 0));
    float p2y = p2.pt() * std::sin(RecoDecay::constrainAngle((p2.phi() + ranPhi), 0));
    q = RecoDecay::p((p1.px() - p2x), (p1.py() - p2y), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2x, p2y, p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1n1_qinv_") + HIST(SubDirHist[part_pair]), qinv, corfac);
  }

  template <ParticleType part, RecGenType rec_gen, typename C, typename T>
  void analyzeSingles(C const& col, T const& tracks)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirHist[] = {"LaP", "LaM"};

    int ntrk = 0;

    for (auto const& track : tracks) {
      // count tracks
      ++ntrk;

      // QA Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_mass_") + HIST(SubDirHist[part]), track.mass());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_pt_") + HIST(SubDirHist[part]), track.pt(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_eta_") + HIST(SubDirHist[part]), track.eta(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_phi_") + HIST(SubDirHist[part]), track.phi(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h1d_n1_rap_") + HIST(SubDirHist[part]), track.rap(), track.corrFact());

      // Efficiency Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h1f_n1_pt_") + HIST(SubDirHist[part]), track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h2f_n1_pteta_") + HIST(SubDirHist[part]), track.pt(), track.eta());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h2f_n1_ptrap_") + HIST(SubDirHist[part]), track.pt(), track.rap());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h4f_n1_ptetavzmult_") + HIST(SubDirHist[part]), track.pt(), track.eta(), col.posZ(), col.cent());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h4f_n1_ptrapvzmult_") + HIST(SubDirHist[part]), track.pt(), track.rap(), col.posZ(), col.cent());

      // Rho1 for N1RapPhi
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h2d_n1_rapphi_") + HIST(SubDirHist[part]), track.rap(), track.phi(), track.corrFact());
    }

    // fill multiplicity histograms
    if (ntrk != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_lambda_mult"), ntrk);
      } else {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h1i_antilambda_mult"), ntrk);
      }
    }
  }

  template <ParticlePairType partpair, RecGenType rec_gen, bool samelambda, typename T>
  void analyzePairs(T const& trks_1, T const& trks_2)
  {
    for (auto const& trk_1 : trks_1) {
      for (auto const& trk_2 : trks_2) {
        // check for same index for Lambda-Lambda / AntiLambda-AntiLambda
        if (samelambda && ((trk_1.index() == trk_2.index()))) {
          continue;
        }
        fillPairHistos<partpair, rec_gen>(trk_1, trk_2);
      }
    }
  }

  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaTracks = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;

  SliceCache cache;
  Partition<LambdaTracks> partLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrackext::trueLambdaFlag == true);
  Partition<LambdaTracks> partAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrackext::trueLambdaFlag == true);

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
    histos.fill(HIST("Event/McGen/h1f_ft0m_mult_percentile"), mcgencol.cent());

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
    adaptAnalysisTask<LambdaTableProducer>(cfgc),
    adaptAnalysisTask<LambdaTracksExtProducer>(cfgc),
    adaptAnalysisTask<LambdaR2Correlation>(cfgc)};
}
