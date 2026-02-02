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

/// \file Lambdacascadecorrelation.cxx
/// \brief Correlation-balance functions of multistrange baryons
/// \author Oveis Sheibani <oveis.sheibani@cern.ch>

#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include <Math/Vector4D.h>
#include <TFile.h>
#include <TH2F.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TProfile.h>
#include <TTree.h> // Required for TTree output

#include <array>
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using namespace o2::soa;

// Zorro zorro;

// use parameters + cov mat non-propagated, aux info + (extension propagated)
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using FullTracksExtWithPID = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;
using FullTracksExtIUWithPID = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr>;

//

namespace o2::aod
{
namespace cascadeflags
{
DECLARE_SOA_COLUMN(IsSelected, isSelected, int); //~!
} // namespace cascadeflags
DECLARE_SOA_TABLE(CascadeFlags, "AOD", "CASCADEFLAGS", //!
                  cascadeflags::IsSelected);
using CascDataExtSelected = soa::Join<CascDataExt, CascadeFlags>;
} // namespace o2::aod

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
using MyCollisionsMult = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using MyCascades = soa::Filtered<aod::CascDataExtSelected>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

namespace o2::aod
{
namespace lambdacollision
{
DECLARE_SOA_COLUMN(Cent, cent, float);
DECLARE_SOA_COLUMN(Mult, mult, float);
DECLARE_SOA_COLUMN(RefCollId, refCollId, int64_t); // <--- 1. Add this line
} // namespace lambdacollision
DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
                  lambdacollision::RefCollId, // <--- 2. Add this line
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ);
using LambdaCollision = LambdaCollisions::iterator;

namespace lambdamcgencollision
{
}
DECLARE_SOA_TABLE(LambdaMcGenCollisions, "AOD", "LMCGENCOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
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
DECLARE_SOA_COLUMN(CosPA, cosPA, float);
DECLARE_SOA_COLUMN(DcaDau, dcaDau, float);
DECLARE_SOA_COLUMN(V0Type, v0Type, int8_t);
DECLARE_SOA_COLUMN(V0PrmScd, v0PrmScd, int8_t);
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
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::V0Type,
                  lambdatrack::V0PrmScd,
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
                  lambdatrack::V0PrmScd,
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
  kV0KShortMassRej,
  kNotLambdaNotAntiLambda,
  kV0IsBothLambdaAntiLambda,
  kNotLambdaAfterSel,
  kV0IsLambdaOrAntiLambda,
  kPassV0DauTrackSel,
  kPassV0KinCuts,
  kPassV0TopoSel,
  kAllSelPassed,
  kPrimaryLambda,
  kSecondaryLambda,
  kLambdaDauNotMcParticle,
  kLambdaNotPrPiMinus,
  kAntiLambdaNotAntiPrPiPlus,
  kPassTrueLambdaSel,
  kEffCorrPtCent,
  kEffCorrPtRapCent,
  kNoEffCorr,
  kPFCorrPtCent,
  kPFCorrPtRapCent,
  kNoPFCorr,
  kGenTotAccLambda,
  kGenLambdaNoDau,
  kGenLambdaToPrPi
};

enum CentEstType {
  kCentFT0M = 0,
  kCentFV0A
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

enum CorrHistDim {
  OneDimCorr = 1,
  TwoDimCorr,
  ThreeDimCorr
};

enum PrmScdType {
  kPrimary = 0,
  kSecondary
};

enum PrmScdPairType {
  kPP = 0,
  kPS,
  kSP,
  kSS
};

struct LambdaTableProducer {

  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::LambdaMcGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMcGenTracks> lambdaMCGenTrackTable;

  // Collisions
  Configurable<int> cCentEstimator{"cCentEstimator", 0, "Centrality Estimator : 0-FT0M, 1-FV0A"};
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinMult{"cMinMult", 0., "Minumum Multiplicity"};
  Configurable<float> cMaxMult{"cMaxMult", 100.0, "Maximum Multiplicity"};
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
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.15, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 999.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "TPC Min Crossed Rows"};
  Configurable<float> cMinTpcCROverCls{"cMinTpcCROverCls", 0.8, "Tpc Min Crossed Rows Over Findable Clusters"};
  Configurable<float> cMaxTpcSharedClusters{"cMaxTpcSharedClusters", 0.4, "Tpc Max Shared Clusters"};
  Configurable<float> cMaxChi2Tpc{"cMaxChi2Tpc", 4, "Max Chi2 Tpc"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 3.0, "TPC NSigma Selection Cut"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  // V0s
  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.02, "Minimum Proton DCAr to PV"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.06, "Minimum Pion DCAr to PV"};
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Minimum DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1., "Maximum DCA between V0 daughters"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Minimum DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 999.0, "Maximum DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Minimum ctau"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Maximum ctau"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.01, "Reject K0Short Candidates"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0short Mass Rej Flag"};

  // V0s kinmatic acceptance
  Configurable<float> cMinV0Mass{"cMinV0Mass", 1.10, "V0 Mass Min"};
  Configurable<float> cMaxV0Mass{"cMaxV0Mass", 1.12, "V0 Mass Min"};
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.8, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 4.2, "Minimum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.5, "|rap| cut"};
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Do Eta Analysis"};
  Configurable<bool> cV0TypeSelFlag{"cV0TypeSelFlag", false, "V0 Type Selection Flag"};
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};

  // V0s MC
  Configurable<bool> cHasMcFlag{"cHasMcFlag", true, "Has Mc Tag"};
  Configurable<bool> cSelectTrueLambda{"cSelectTrueLambda", true, "Select True Lambda"};
  Configurable<bool> cSelMCPSV0{"cSelMCPSV0", true, "Select Primary/Secondary V0"};
  Configurable<bool> cCheckRecoDauFlag{"cCheckRecoDauFlag", true, "Check for reco daughter PID"};
  Configurable<bool> cGenPrimaryLambda{"cGenPrimaryLambda", true, "Primary Generated Lambda"};
  Configurable<bool> cGenSecondaryLambda{"cGenSecondaryLambda", false, "Secondary Generated Lambda"};
  Configurable<bool> cGenDecayChannel{"cGenDecayChannel", true, "Gen Level Decay Channel Flag"};
  Configurable<bool> cRecoMomResoFlag{"cRecoMomResoFlag", false, "Check effect of momentum space smearing on balance function"};

  // Efficiency Correction
  Configurable<bool> cCorrectionFlag{"cCorrectionFlag", false, "Correction Flag"};
  Configurable<bool> cGetEffFact{"cGetEffFact", false, "Get Efficiency Factor Flag"};
  Configurable<bool> cGetPrimFrac{"cGetPrimFrac", false, "Get Primary Fraction Flag"};
  Configurable<int> cCorrFactHist{"cCorrFactHist", 0, "Efficiency Factor Histogram"};
  Configurable<int> cPrimFracHist{"cPrimFracHist", 0, "Primary Fraction Histogram"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cPathCCDB{"cPathCCDB", "Users/y/ypatley/lambda_corr_fact", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vCorrFactStrings = {{"hEffVsPtCentLambda", "hEffVsPtCentAntiLambda"},
                                                            {"hEffVsPtYCentLambda", "hEffVsPtYCentAntiLambda"},
                                                            {"hEffVsPtEtaCentLambda", "hEffVsPtEtaCentAntiLambda"}};

  // initialize corr_factor objects
  std::vector<std::vector<std::string>> vPrimFracStrings = {{"hPrimFracVsPtCentLambda", "hPrimFracVsPtCentAntiLambda"},
                                                            {"hPrimFracVsPtYCentLambda", "hPrimFracVsPtYCentAntiLambda"},
                                                            {"hPrimFracVsPtEtaCentLambda", "hPrimFracVsPtEtaCentAntiLambda"}};

  // Initialize Global Variables
  float cent = 0., mult = 0.;
  float pt = 0., eta = 0., rap = 0., phi = 0.;

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdb->setURL(cUrlCCDB.value);
    ccdb->setCaching(true);

    // initialize axis specifications
    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(30, 0.5, 30.5, "");
    const AxisSpec axisCent(100, 0, 100, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(140, 1.08, 1.15, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(100., 0., 10., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(48, -1.2, 1.2, "y");
    const AxisSpec axisV0Eta(48, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

    const AxisSpec axisRadius(2000, 0, 200, "r(cm)");
    const AxisSpec axisCosPA(300, 0.97, 1.0, "cos(#theta_{PA})");
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

    // Kinematic Histograms
    histos.add("McRec/Lambda/hPt", "Transverse Momentum", kTH1F, {axisV0Pt});
    histos.add("McRec/Lambda/hEta", "Pseudorapidity", kTH1F, {axisV0Eta});
    histos.add("McRec/Lambda/hRap", "Rapidity", kTH1F, {axisV0Rap});
    histos.add("McRec/Lambda/hPhi", "Azimuthal Angle", kTH1F, {axisV0Phi});

    // QA Anti-Lambda
    histos.addClone("QA/Lambda/", "QA/AntiLambda/");
    histos.addClone("McRec/Lambda/", "McRec/AntiLambda/");

    // MC Generated Histograms
    if (doprocessMCRun3 || doprocessMCRun2 || doprocessMCRecoRun3 || doprocessMCRecoRun2) {
      // McReco Histos
      histos.add("Tracks/h2f_tracks_pid_before_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_tracks_pid_after_sel", "PIDs", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_lambda_mothers_pdg", "PIDs", kTH2F, {axisPID, axisV0Pt});

      // McGen Histos
      histos.add("McGen/h1f_collision_recgen", "# of Reco Collision Associated to One Mc Generator Collision", kTH1F, {axisMult});
      histos.add("McGen/h1f_collisions_info", "# of collisions", kTH1F, {axisCols});
      histos.add("McGen/h2f_collision_posZ", "V_{z}-distribution", kTH2F, {axisVz, axisVz});
      histos.add("McGen/h2f_collision_cent", "FT0M Centrality", kTH2F, {axisCent, axisCent});
      histos.add("McGen/h1f_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});

      histos.addClone("McRec/", "McGen/");

      histos.add("McGen/Lambda/Proton/hPt", "Proton p_{T}", kTH1F, {axisTrackPt});
      histos.add("McGen/Lambda/Proton/hEta", "Proton #eta", kTH1F, {axisV0Eta});
      histos.add("McGen/Lambda/Proton/hRap", "Proton y", kTH1F, {axisV0Rap});
      histos.add("McGen/Lambda/Proton/hPhi", "Proton #phi", kTH1F, {axisV0Phi});

      histos.addClone("McGen/Lambda/Proton/", "McGen/Lambda/Pion/");
      histos.addClone("McGen/Lambda/Proton/", "McGen/AntiLambda/Proton/");
      histos.addClone("McGen/Lambda/Pion/", "McGen/AntiLambda/Pion/");

      // set bin lables specific to MC
      histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotColBeforeHasMcCollision, "kTotColBeforeHasMcCollision");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kTracksBeforeHasMcParticle, "kTracksBeforeHasMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPrimaryLambda, "kPrimaryLambda");
      histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kSecondaryLambda, "kSecondaryLambda");
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
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0KShortMassRej, "kV0KShortMassRej");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotLambdaNotAntiLambda, "kNotLambdaNotAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0IsBothLambdaAntiLambda, "kV0IsBothLambdaAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNotLambdaAfterSel, "kNotLambdaAfterSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kV0IsLambdaOrAntiLambda, "kV0IsLambdaOrAntiLambda");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0DauTrackSel, "kPassV0DauTrackSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0KinCuts, "kPassV0KinCuts");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPassV0TopoSel, "kPassV0TopoSel");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kAllSelPassed, "kAllSelPassed");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtCent, "kEffCorrPtCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kEffCorrPtRapCent, "kEffCorrPtRapCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNoEffCorr, "kNoEffCorr");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPFCorrPtCent, "kPFCorrPtCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kPFCorrPtRapCent, "kPFCorrPtRapCent");
    histos.get<TH1>(HIST("Tracks/h1f_tracks_info"))->GetXaxis()->SetBinLabel(TrackLabels::kNoPFCorr, "kNoPFCorr");
  }

  template <RunType run, typename C>
  bool selCollision(C const& col)
  {
    // VtxZ Selection
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }

    if constexpr (run == kRun3) { // Run3 Min-Bias Trigger
      // select centrality estimator
      if (cCentEstimator == kCentFT0M) {
        cent = col.centFT0M();
      } else if (cCentEstimator == kCentFV0A) {
        cent = col.centFV0A();
      }
      if (cSel8Trig && !col.sel8()) {
        return false;
      }
    } else { // Run2 Min-Bias Trigger
      cent = col.centRun2V0M();
      if (cInt7Trig && !col.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !col.sel7()) {
        return false;
      }
    }

    if (cent <= cMinMult || cent >= cMaxMult) { // select centrality percentile class
      return false;
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

    // Set Multiplicity
    mult = col.multNTracksPV();

    return true;
  }

  // Kinematic Selection
  bool kinCutSelection(float const& pt, float const& rap, float const& ptMin, float const& ptMax, float const& rapMax)
  {
    if (pt <= ptMin || pt >= ptMax || rap >= rapMax) {
      return false;
    }

    return true;
  }

  // Track Selection
  template <typename T>
  bool selTrack(T const& track)
  {
    if (!kinCutSelection(track.pt(), std::abs(track.eta()), cTrackMinPt, cTrackMaxPt, cTrackEtaCut)) {
      return false;
    }

    if (track.tpcNClsCrossedRows() <= cMinTpcCrossedRows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < cMinTpcCROverCls) {
      return false;
    }

    if (track.tpcNClsShared() > cMaxTpcSharedClusters) {
      return false;
    }

    if (track.tpcChi2NCl() > cMaxChi2Tpc) {
      return false;
    }

    return true;
  }

  // Daughter Track Selection
  template <typename V, typename T>
  bool selDaughterTracks(V const& v0, T const&, ParticleType const& v0Type)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    if (!selTrack(posTrack) || !selTrack(negTrack)) {
      return false;
    }

    // Apply DCA Selection on Daughter Tracks Based on Lambda/AntiLambda daughters
    float dcaProton = 0., dcaPion = 0.;
    if (v0Type == kLambda) {
      dcaProton = std::abs(v0.dcapostopv());
      dcaPion = std::abs(v0.dcanegtopv());
    } else if (v0Type == kAntiLambda) {
      dcaPion = std::abs(v0.dcapostopv());
      dcaProton = std::abs(v0.dcanegtopv());
    }

    if (dcaProton < cMinDcaProtonToPV || dcaPion < cMinDcaPionToPV) {
      return false;
    }

    return true;
  }

  template <typename C, typename V, typename T>
  bool topoCutSelection(C const& col, V const& v0, T const&)
  {
    // DCA
    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters) {
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
    // Kshort mass rejection hypothesis
    if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0KShortMassRej);
      return false;
    }

    // initialize daughter tracks
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // initialize selection flags
    bool lambdaFlag = false, antiLambdaFlag = false;

    // get v0 track as lambda
    if ((v0.mLambda() > cMinV0Mass && v0.mLambda() < cMaxV0Mass) && (selLambdaDauWithTpcPid<kLambda>(postrack, negtrack))) {
      lambdaFlag = true;
      v0type = kLambda;
    }

    // get v0 track as anti-lambda
    if ((v0.mAntiLambda() > cMinV0Mass && v0.mAntiLambda() < cMaxV0Mass) && (selLambdaDauWithTpcPid<kAntiLambda>(postrack, negtrack))) {
      antiLambdaFlag = true;
      v0type = kAntiLambda;
    }

    if (!lambdaFlag && !antiLambdaFlag) { // neither Lambda nor Anti-Lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaNotAntiLambda);
      return false;
    } else if (lambdaFlag && antiLambdaFlag) { // check if the track is identified as lambda and anti-lambda both (DISCARD THIS TRACK)
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsBothLambdaAntiLambda);
      return false;
    }

    if (lambdaFlag || antiLambdaFlag) {
      return true;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaAfterSel);

    return false;
  }

  template <typename C, typename V, typename T>
  bool selV0Particle(C const& col, V const& v0, T const& tracks, ParticleType& v0Type)
  {
    // Apply Lambda Mass Hypothesis
    if (!selLambdaMassWindow(v0, tracks, v0Type)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsLambdaOrAntiLambda);

    // Apply Daughter Track Selection
    if (!selDaughterTracks(v0, tracks, v0Type)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0DauTrackSel);

    // Apply Kinematic Selection
    float rap = 0.;
    if (!cDoEtaAnalysis) {
      rap = std::abs(v0.yLambda());
    } else {
      rap = std::abs(v0.eta());
    }

    if (!kinCutSelection(v0.pt(), rap, cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0KinCuts);

    // Apply Topological Selection
    if (!topoCutSelection(col, v0, tracks)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0TopoSel);

    // All Selection Criterion Passed
    return true;
  }

  template <typename V, typename T>
  bool hasAmbiguousDaughters(V const& v0, T const&)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();

    auto posTrackCompCols = posTrack.compatibleCollIds();
    auto negTrackCompCols = negTrack.compatibleCollIds();

    // Check if daughter tracks belongs to more than one collision (Ambiguous Tracks)
    if (posTrackCompCols.size() > 1 || negTrackCompCols.size() > 1) {
      return true;
    }

    // Check if compatible collision index matches the track collision index
    if (((posTrackCompCols.size() != 0) && (posTrackCompCols[0] != posTrack.collisionId())) ||
        ((negTrackCompCols.size() != 0) && (negTrackCompCols[0] != negTrack.collisionId()))) {
      return true;
    }

    // Pass as not ambiguous
    return false;
  }

  template <typename V>
  PrmScdType isPrimaryV0(V const& v0)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check for secondary lambda
    if (!mcpart.isPhysicalPrimary()) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kSecondaryLambda);
      return kSecondary;
    }

    histos.fill(HIST("Tracks/h1f_tracks_info"), kPrimaryLambda);
    return kPrimary;
  }

  template <typename V, typename T>
  bool selTrueMcRecLambda(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();

    // check if Lambda/AntiLambda
    if (std::abs(mcpart.pdgCode()) != kLambda0) {
      return false;
    }

    // Check for daughters
    if (cCheckRecoDauFlag) {
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
    }

    return true;
  }

  template <ParticleType part, typename V>
  float getCorrectionFactors(V const& v0)
  {
    // Check for efficiency correction flag
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

    // initialize efficiency factor and primary fraction values
    float effCorrFact = 1., primFrac = 1.;
    float rap = (cDoEtaAnalysis) ? v0.eta() : v0.yLambda();

    // Get Efficiency Factor
    if (cGetEffFact) {
      TObject* objEff = reinterpret_cast<TObject*>(ccdbObj->FindObject(Form("%s", vCorrFactStrings[cCorrFactHist][part].c_str())));
      TH1F* histEff = reinterpret_cast<TH1F*>(objEff->Clone());
      if (histEff->GetDimension() == TwoDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtCent);
        effCorrFact = histEff->GetBinContent(histEff->FindBin(cent, v0.pt()));
      } else if (histEff->GetDimension() == ThreeDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtRapCent);
        effCorrFact = histEff->GetBinContent(histEff->FindBin(cent, v0.pt(), rap));
      } else {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kNoEffCorr);
        LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
        effCorrFact = 1.;
      }
      delete histEff;
    }

    // Get Primary Fraction
    // (The dimension of this could be different than efficiency because of large errors !!!)
    if (cGetPrimFrac) {
      TObject* objPrm = reinterpret_cast<TObject*>(ccdbObj->FindObject(Form("%s", vPrimFracStrings[cPrimFracHist][part].c_str())));
      TH1F* histPrm = reinterpret_cast<TH1F*>(objPrm->Clone());
      if (histPrm->GetDimension() == TwoDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPFCorrPtCent);
        primFrac = histPrm->GetBinContent(histPrm->FindBin(cent, v0.pt()));
      } else if (histPrm->GetDimension() == ThreeDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPFCorrPtRapCent);
        primFrac = histPrm->GetBinContent(histPrm->FindBin(cent, v0.pt(), rap));
      } else {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kNoPFCorr);
        LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
        primFrac = 1.;
      }
      delete histPrm;
    }

    return primFrac * effCorrFact;
  }

  template <typename V, typename T>
  void fillLambdaMothers(V const& v0, T const&)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();
    auto lambdaMothers = mcpart.template mothers_as<aod::McParticles>();
    histos.fill(HIST("Tracks/h2f_lambda_mothers_pdg"), lambdaMothers[0].pdgCode(), v0.pt());
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

  // Fill Lambda Kinematic Histograms
  template <RecGenType rg, ParticleType part>
  void fillKinematicHists(float const& pt, float const& eta, float const& y, float const& phi)
  {
    static constexpr std::string_view SubDirRG[] = {"McRec/", "McGen/"};
    static constexpr std::string_view SubDirPart[] = {"Lambda/", "AntiLambda/"};

    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPt"), pt);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hEta"), eta);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hRap"), y);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPhi"), phi);
  }

  // Reconstructed Level Tables
  template <RunType run, DMCType dmc, typename C, typename V, typename T>
  void fillLambdaRecoTables(C const& collision, V const& v0tracks, T const& tracks)
  {
    // Total Collisions
    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    // Select Collision (Only for Data... McRec has been selected already !!!)
    if constexpr (dmc == kData) {
      if (!selCollision<run>(collision)) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());

    // Fill Collision Table
    // lambdaCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ());
    lambdaCollisionTable(cent, mult, collision.globalIndex(), collision.posX(), collision.posY(), collision.posZ());

    // initialize v0track objects
    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float mass = 0., corr_fact = 1.;

    for (auto const& v0 : v0tracks) {
      // check for corresponding MCGen Particle
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kTracksBeforeHasMcParticle);
        if (!v0.has_mcParticle()) {
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // Select V0 Particle as Lambda/AntiLambda
      if (!selV0Particle(collision, v0, tracks, v0Type)) {
        continue;
      }

      // Select V0 Type Selection
      if (cV0TypeSelFlag && v0.v0Type() != cV0TypeSelection) {
        continue;
      }

      // we have v0 as lambda
      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllSelPassed);

      // Remove lambda with ambiguous daughters (Only for run3)
      if constexpr (run == kRun3) {
        if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks)) {
          continue;
        }
      }

      // Get Lambda mass and kinematic variables
      mass = (v0Type == kLambda) ? v0.mLambda() : v0.mAntiLambda();
      pt = v0.pt();
      eta = v0.eta();
      rap = v0.yLambda();
      phi = v0.phi();

      // do MC analysis
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());

        // Get Primary/Secondary Lambda
        if (cSelMCPSV0) {
          v0PrmScdType = isPrimaryV0(v0);
        }

        // check for true Lambda/Anti-Lambda
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks)) {
          continue;
        }

        // get mothers information
        if (v0PrmScdType == kSecondary) {
          fillLambdaMothers(v0, tracks);
        }

        histos.fill(HIST("Tracks/h1f_tracks_info"), kPassTrueLambdaSel);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());

        if (cRecoMomResoFlag) {
          auto mc = v0.template mcParticle_as<aod::McParticles>();
          pt = mc.pt();
          eta = mc.eta();
          rap = mc.y();
          phi = mc.phi();
          float y = (cDoEtaAnalysis) ? eta : rap;
          // apply kinematic selection (On Truth)
          if (!kinCutSelection(pt, std::abs(y), cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
            continue;
          }
        }
      }

      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

      // get correction factors
      corr_fact = (v0Type == kLambda) ? getCorrectionFactors<kLambda>(v0) : getCorrectionFactors<kAntiLambda>(v0);

      // fill lambda qa
      if (v0Type == kLambda) {
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      } else {
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kAntiLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      }

      // Fill Lambda/AntiLambda Table
      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0.px(), v0.py(), v0.pz(),
                       pt, eta, phi, rap, mass, v0.template posTrack_as<T>().index(), v0.template negTrack_as<T>().index(),
                       v0.v0cosPA(), v0.dcaV0daughters(), (int8_t)v0Type, v0PrmScdType, corr_fact);
    }
  }

  // MC Generater Level Tables
  template <RunType run, typename C, typename M>
  void fillLambdaMcGenTables(C const& mcCollision, M const& mcParticles)
  {
    // Fill McGen Collision Table
    lambdaMCGenCollisionTable(cent, mult, mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

    // initialize track objects
    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float rap = 0.;

    for (auto const& mcpart : mcParticles) {
      // check for Lambda first
      if (mcpart.pdgCode() == kLambda0) {
        v0Type = kLambda;
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        v0Type = kAntiLambda;
      } else {
        continue;
      }

      // check for Primary Lambda/AntiLambda
      if (mcpart.isPhysicalPrimary()) {
        v0PrmScdType = kPrimary;
      } else {
        v0PrmScdType = kSecondary;
      }

      // Decide Eta/Rap
      if (!cDoEtaAnalysis) {
        rap = mcpart.y();
      } else {
        rap = mcpart.eta();
      }

      // Apply Kinematic Acceptance
      if (!kinCutSelection(mcpart.pt(), std::abs(rap), cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
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
      std::vector<float> vDauPt, vDauEta, vDauRap, vDauPhi;
      for (auto const& dautrack : dautracks) {
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
        vDauPt.push_back(dautrack.pt());
        vDauEta.push_back(dautrack.eta());
        vDauRap.push_back(dautrack.y());
        vDauPhi.push_back(dautrack.phi());
      }
      if (cGenDecayChannel) { // check decay channel
        if (v0Type == kLambda) {
          if (daughterPDGs[0] != kProton || daughterPDGs[1] != kPiMinus) {
            continue;
          }
        } else if (v0Type == kAntiLambda) {
          if (daughterPDGs[0] != kProtonBar || daughterPDGs[1] != kPiPlus) {
            continue;
          }
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaToPrPi);

      if (v0Type == kLambda) {
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/Lambda/Proton/hPt"), vDauPt[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hEta"), vDauEta[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hRap"), vDauRap[0]);
        histos.fill(HIST("McGen/Lambda/Proton/hPhi"), vDauPhi[0]);
        histos.fill(HIST("McGen/Lambda/Pion/hPt"), vDauPt[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hEta"), vDauEta[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hRap"), vDauRap[1]);
        histos.fill(HIST("McGen/Lambda/Pion/hPhi"), vDauPhi[1]);
        fillKinematicHists<kGen, kLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      } else {
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[0]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[1]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/AntiLambda/Pion/hPt"), vDauPt[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hEta"), vDauEta[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hRap"), vDauRap[0]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hPhi"), vDauPhi[0]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hPt"), vDauPt[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hEta"), vDauEta[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hRap"), vDauRap[1]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hPhi"), vDauPhi[1]);
        fillKinematicHists<kGen, kAntiLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      }

      // Fill Lambda McGen Table
      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            daughterIDs[0], daughterIDs[1], (int8_t)v0Type, -999., -999., v0PrmScdType, 1.);
    }
  }

  template <RunType run, DMCType dmc, typename M, typename C, typename V, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, V const& V0s, T const& tracks, P const& mcParticles)
  {
    // Number of Rec Collisions Associated to the McGen Collision
    int nRecCols = collisions.size();
    if (nRecCols != 0) {
      histos.fill(HIST("McGen/h1f_collision_recgen"), nRecCols);
    }
    // Do not analyze if more than one reco collision is accociated to one mc gen collision
    if (nRecCols != 1) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kTotCol);
    // Check the reco collision
    if (!collisions.begin().has_mcCollision() || !selCollision<run>(collisions.begin()) || collisions.begin().mcCollisionId() != mcCollision.globalIndex()) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(), collisions.begin().posZ());
    auto v0Tracks = V0s.sliceBy(perCollision, collisions.begin().globalIndex());
    fillLambdaRecoTables<run, dmc>(collisions.begin(), v0Tracks, tracks);
    fillLambdaMcGenTables<run>(mcCollision, mcParticles);
  }

  SliceCache cache;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCollision = aod::v0data::collisionId;

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFV0As, aod::PVMults>;
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::PVMults>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::TrackCompColls>;
  using TracksRun2 = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using TracksMCRun2 = soa::Join<TracksRun2, aod::McTrackLabels>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processDataRun3(CollisionsRun3::iterator const& collision, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    fillLambdaRecoTables<kRun3, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processDataRun3, "Process for Run3 DATA", true);

  void processDataRun2(CollisionsRun2::iterator const& collision, aod::V0Datas const& V0s, TracksRun2 const& tracks)
  {
    fillLambdaRecoTables<kRun2, kData>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processDataRun2, "Process for Run2 DATA", false);

  void processMCRecoRun3(soa::Join<CollisionsRun3, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, TracksMC const& tracks, aod::McParticles const&)
  {
    // check collision
    if (!selCollision<kRun3>(collision)) {
      return;
    }
    fillLambdaRecoTables<kRun3, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRecoRun3, "Process for Run3 McReco DATA", false);

  void processMCRecoRun2(soa::Join<CollisionsRun2, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&,
                         McV0Tracks const& V0s, TracksMCRun2 const& tracks, aod::McParticles const&)
  {
    // check collision
    if (!selCollision<kRun2>(collision)) {
      return;
    }
    fillLambdaRecoTables<kRun2, kMC>(collision, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRecoRun2, "Process for Run2 McReco DATA", false);

  void processMCRun3(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun3, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMC const& tracks,
                     aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun3, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRun3, "Process for Run3 MC RecoGen", false);

  void processMCRun2(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun2, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMCRun2 const& tracks,
                     aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun2, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRun2, "Process for Run2 MC RecoGen", false);
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
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    // Histograms Booking
    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_lambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_antilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h2d_n2_etaphi_LaP_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaP_LaP", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaM_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});

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

        // check if lambda shares daughters with any other track
        if (lambda.posTrackId() == track.posTrackId() || lambda.negTrackId() == track.negTrackId()) {
          vSharedDauLambdaIndex.push_back(track.index());
          lambdaSharingDauFlag = true;

          // Fill DEta-DPhi Histogram
          if ((lambda.v0Type() == kLambda && track.v0Type() == kAntiLambda) || (lambda.v0Type() == kAntiLambda && track.v0Type() == kLambda)) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.v0Type() == kLambda && track.v0Type() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaP"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          }

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
  Configurable<int> cNPtBins{"cNPtBins", 34, "N pT Bins"};
  Configurable<float> cMinPt{"cMinPt", 0.8, "pT Min"};
  Configurable<float> cMaxPt{"cMaxPt", 4.2, "pT Max"};
  Configurable<int> cNRapBins{"cNRapBins", 20, "N Rapidity Bins"};
  Configurable<float> cMinRap{"cMinRap", -0.5, "Minimum Rapidity"};
  Configurable<float> cMaxRap{"cMaxRap", 0.5, "Maximum Rapidity"};
  Configurable<int> cNPhiBins{"cNPhiBins", 36, "N Phi Bins"};
  Configurable<bool> cAnaSecondaries{"cAnaSecondaries", false, "Analysze Secondaries"};
  Configurable<bool> cAnaPairs{"cAnaPairs", false, "Analyze Pairs Flag"};
  Configurable<bool> cAnaSecondaryPairs{"cAnaSecondaryPairs", false, "Analyze Secondary Pairs Flag"};

  // Eta/Rap Analysis
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Eta/Rap Analysis Flag"};

  // Centrality Axis
  ConfigurableAxis cMultBins{"cMultBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 30.0f, 50.f, 80.0f, 100.f}, "Variable Mult-Bins"};

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
  float cent = 0.;

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
    const AxisSpec axisCent(cMultBins, "FT0M (%)");
    const AxisSpec axisChMult(200, 0, 200, "N_{ch}");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisMass(100, 1.06, 1.16, "M_{#Lambda} (GeV/#it{c}^{2})");
    const AxisSpec axisPt(cNPtBins, cMinPt, cMaxPt, "p_{T} (GeV/#it{c})");
    const AxisSpec axisEta(cNRapBins, cMinRap, cMaxRap, "#eta");
    const AxisSpec axisRap(cNRapBins, cMinRap, cMaxRap, "y");
    const AxisSpec axisPhi(cNPhiBins, 0., TwoPI, "#varphi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "y #varphi");
    const AxisSpec axisQinv(100, 0, 10, "q_{inv} (GeV/#it{c})");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1f_collision_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/Reco/h1f_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h2f_Mult_vs_Centrality", "N_{ch} vs FT0M(%)", kTH2F, {axisCent, axisChMult});
    histos.add("Event/Reco/h2f_lambda_mult", "#Lambda - Multiplicity", kTH2F, {axisCent, axisMult});
    histos.add("Event/Reco/h2f_antilambda_mult", "#bar{#Lambda} - Multiplicity", kTH2F, {axisCent, axisMult});

    // Efficiency Histograms
    // Single Particle Efficiencies
    histos.add("Reco/Primary/Efficiency/h2f_n1_centpt_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisPt});
    histos.add("Reco/Primary/Efficiency/h2f_n1_centpt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisPt});
    histos.add("Reco/Primary/Efficiency/h3f_n1_centpteta_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPt, axisEta});
    histos.add("Reco/Primary/Efficiency/h3f_n1_centpteta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPt, axisEta});
    histos.add("Reco/Primary/Efficiency/h3f_n1_centptrap_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPt, axisRap});
    histos.add("Reco/Primary/Efficiency/h3f_n1_centptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPt, axisRap});

    // Single and Two Particle Densities
    // 1D Histograms
    histos.add("Reco/Primary/h3f_n1_centmasspt_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisMass, axisPt});
    histos.add("Reco/Primary/h3f_n1_centmasspt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisMass, axisPt});
    histos.add("Reco/Primary/h2f_n1_pt_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisPt});
    histos.add("Reco/Primary/h2f_n1_pt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisPt});
    histos.add("Reco/Primary/h2f_n1_eta_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisEta});
    histos.add("Reco/Primary/h2f_n1_eta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisEta});
    histos.add("Reco/Primary/h2f_n1_rap_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisRap});
    histos.add("Reco/Primary/h2f_n1_rap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisRap});
    histos.add("Reco/Primary/h2f_n1_phi_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisPhi});
    histos.add("Reco/Primary/h2f_n1_phi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisPhi});

    // rho1 for R2 RapPhi
    histos.add("Reco/Primary/h3f_n1_rapphi_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisRap, axisPhi});
    histos.add("Reco/Primary/h3f_n1_rapphi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisRap, axisPhi});

    // rho1 for Q_{inv}
    histos.add("Reco/Primary/h3f_n1_pteta_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPt, axisEta});
    histos.add("Reco/Primary/h3f_n1_pteta_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPt, axisEta});

    // Clone Singles Primary/Secondary Histogram
    if (cAnaSecondaries) {
      histos.addClone("Reco/Primary/", "Reco/Secondary/");
    }

    if (cAnaPairs) {
      // rho2 for numerator of R2
      histos.add("Reco/PP/h3f_n2_raprap_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH3F, {axisCent, axisRap, axisRap});
      histos.add("Reco/PP/h3f_n2_raprap_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH3F, {axisCent, axisRap, axisRap});
      histos.add("Reco/PP/h3f_n2_raprap_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH3F, {axisCent, axisRap, axisRap});
      histos.add("Reco/PP/h3f_n2_phiphi_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH3F, {axisCent, axisPhi, axisPhi});
      histos.add("Reco/PP/h3f_n2_phiphi_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH3F, {axisCent, axisPhi, axisPhi});
      histos.add("Reco/PP/h3f_n2_phiphi_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH3F, {axisCent, axisPhi, axisPhi});

      // rho2 for R2 Rap1Phi1Rap2Phi2
      histos.add("Reco/PP/h3f_n2_rapphi_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/PP/h3f_n2_rapphi_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/PP/h3f_n2_rapphi_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});

      // rho2 for R2 Qinv
      histos.add("Reco/PP/h2f_n2_qinv_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH2F, {axisCent, axisQinv});
      histos.add("Reco/PP/h2f_n2_qinv_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH2F, {axisCent, axisQinv});
      histos.add("Reco/PP/h2f_n2_qinv_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH2F, {axisCent, axisQinv});

      // Clone Pairs Histograms
      if (cAnaSecondaryPairs) {
        histos.addClone("Reco/PP/", "Reco/PS/");
        histos.addClone("Reco/PP/", "Reco/SP/");
        histos.addClone("Reco/PP/", "Reco/SS/");
      }
    }

    // MCGen
    if (doprocessMCGen) {
      histos.addClone("Event/Reco/", "Event/McGen/");
      histos.addClone("Reco/", "McGen/");
    }
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, PrmScdPairType psp, typename U>
  void fillPairHistos(U& p1, U& p2)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirPrmScd[] = {"PP/", "PS/", "SP/", "SS/"};
    static constexpr std::string_view SubDirHist[] = {"LaP_LaM", "LaP_LaP", "LaM_LaM"};

    float rap1 = (cDoEtaAnalysis) ? p1.eta() : p1.rap();
    float rap2 = (cDoEtaAnalysis) ? p2.eta() : p2.rap();

    int rapbin1 = static_cast<int>((rap1 - kminrap) / rapbinwidth);
    int rapbin2 = static_cast<int>((rap2 - kminrap) / rapbinwidth);

    int phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    int phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    float corfac = p1.corrFact() * p2.corrFact();

    // fill rho2 histograms
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[psp]) + HIST("h3f_n2_raprap_") + HIST(SubDirHist[part_pair]), cent, rap1, rap2, corfac);
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[psp]) + HIST("h3f_n2_phiphi_") + HIST(SubDirHist[part_pair]), cent, p1.phi(), p2.phi(), corfac);

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      int rapphix = rapbin1 * nphibins + phibin1;
      int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[psp]) + HIST("h3f_n2_rapphi_") + HIST(SubDirHist[part_pair]), cent, rapphix + 0.5, rapphiy + 0.5, corfac);
    }

    // qinv histograms
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassLambda0);
    qinv = std::sqrt(-RecoDecay::m2(q, e));
    histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[psp]) + HIST("h2f_n2_qinv_") + HIST(SubDirHist[part_pair]), cent, qinv, corfac);
  }

  template <ParticleType part, RecGenType rec_gen, PrmScdType pst, typename T>
  void analyzeSingles(T const& tracks)
  {
    static constexpr std::string_view SubDirRecGen[] = {"Reco/", "McGen/"};
    static constexpr std::string_view SubDirPrmScd[] = {"Primary/", "Secondary/"};
    static constexpr std::string_view SubDirHist[] = {"LaP", "LaM"};

    int ntrk = 0;

    for (auto const& track : tracks) {
      // count tracks
      ++ntrk;

      // Efficiency Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("Efficiency/h2f_n1_centpt_") + HIST(SubDirHist[part]), cent, track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("Efficiency/h3f_n1_centpteta_") + HIST(SubDirHist[part]), cent, track.pt(), track.eta());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("Efficiency/h3f_n1_centptrap_") + HIST(SubDirHist[part]), cent, track.pt(), track.rap());

      // QA Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h3f_n1_centmasspt_") + HIST(SubDirHist[part]), cent, track.mass(), track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h2f_n1_pt_") + HIST(SubDirHist[part]), cent, track.pt(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h2f_n1_eta_") + HIST(SubDirHist[part]), cent, track.eta(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h2f_n1_phi_") + HIST(SubDirHist[part]), cent, track.phi(), track.corrFact());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h2f_n1_rap_") + HIST(SubDirHist[part]), cent, track.rap(), track.corrFact());

      // Rho1 for N1RapPhi
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h3f_n1_rapphi_") + HIST(SubDirHist[part]), cent, track.rap(), track.phi(), track.corrFact());

      // Rho1 for Q_{inv} Bkg Estimation
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST(SubDirPrmScd[pst]) + HIST("h3f_n1_pteta_") + HIST(SubDirHist[part]), cent, track.pt(), track.eta(), track.corrFact());
    }

    // fill multiplicity histograms
    if (ntrk != 0) {
      if (part == kLambda) {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h2f_lambda_mult"), cent, ntrk);
      } else {
        histos.fill(HIST("Event/") + HIST(SubDirRecGen[rec_gen]) + HIST("h2f_antilambda_mult"), cent, ntrk);
      }
    }
  }

  template <ParticlePairType partpair, RecGenType rec_gen, PrmScdPairType psp, bool samelambda, typename T>
  void analyzePairs(T const& trks_1, T const& trks_2)
  {
    for (auto const& trk_1 : trks_1) {
      for (auto const& trk_2 : trks_2) {
        // check for same index for Lambda-Lambda / AntiLambda-AntiLambda
        if (samelambda && ((trk_1.index() == trk_2.index()))) {
          continue;
        }
        fillPairHistos<partpair, rec_gen, psp>(trk_1, trk_2);
      }
    }
  }

  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaTracks = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;

  // using MyCascades = aod::CascDataExt;   // ← NOT CascDatas. NEVER CascDatas.
  using MyCascades = soa::Filtered<aod::CascDataExt>;

  SliceCache cache;
  Partition<LambdaTracks> partPrimLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrackext::trueLambdaFlag == true) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);
  Partition<LambdaTracks> partPrimAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrackext::trueLambdaFlag == true) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);
  Partition<LambdaTracks> partSecdLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrackext::trueLambdaFlag == true) && (aod::lambdatrack::v0PrmScd == (int8_t)kSecondary);
  Partition<LambdaTracks> partSecdAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrackext::trueLambdaFlag == true) && (aod::lambdatrack::v0PrmScd == (int8_t)kSecondary);

  void processDataReco(LambdaCollisions::iterator const& collision, LambdaTracks const&)
  {
    histos.fill(HIST("Event/Reco/h1f_collision_posz"), collision.posZ());
    histos.fill(HIST("Event/Reco/h1f_ft0m_mult_percentile"), collision.cent());
    histos.fill(HIST("Event/Reco/h2f_Mult_vs_Centrality"), collision.cent(), collision.mult());

    cent = collision.cent();

    auto lambdaPrimTracks = partPrimLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto antiLambdaPrimTracks = partPrimAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto lambdaSecdTracks = partSecdLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto antiLambdaSecdTracks = partSecdAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kRec, kPrimary>(lambdaPrimTracks);
    analyzeSingles<kAntiLambda, kRec, kPrimary>(antiLambdaPrimTracks);

    if (cAnaSecondaries) {
      analyzeSingles<kLambda, kRec, kSecondary>(lambdaSecdTracks);
      analyzeSingles<kAntiLambda, kRec, kSecondary>(antiLambdaSecdTracks);
    }

    if (cAnaPairs) {
      // Primary Pairs Only
      analyzePairs<kLambdaAntiLambda, kRec, kPP, false>(lambdaPrimTracks, antiLambdaPrimTracks);
      analyzePairs<kLambdaLambda, kRec, kPP, true>(lambdaPrimTracks, lambdaPrimTracks);
      analyzePairs<kAntiLambdaAntiLambda, kRec, kPP, true>(antiLambdaPrimTracks, antiLambdaPrimTracks);

      // Secondary Pairs
      if (cAnaSecondaryPairs) {
        analyzePairs<kLambdaAntiLambda, kRec, kPS, false>(lambdaPrimTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaLambda, kRec, kPS, true>(lambdaPrimTracks, lambdaSecdTracks);
        analyzePairs<kAntiLambdaAntiLambda, kRec, kPS, true>(antiLambdaPrimTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaAntiLambda, kRec, kSP, false>(lambdaSecdTracks, antiLambdaPrimTracks);
        analyzePairs<kLambdaLambda, kRec, kSP, true>(lambdaSecdTracks, lambdaPrimTracks);
        analyzePairs<kAntiLambdaAntiLambda, kRec, kSP, true>(antiLambdaSecdTracks, antiLambdaPrimTracks);
        analyzePairs<kLambdaAntiLambda, kRec, kSS, false>(lambdaSecdTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaLambda, kRec, kSS, true>(lambdaSecdTracks, lambdaSecdTracks);
        analyzePairs<kAntiLambdaAntiLambda, kRec, kSS, true>(antiLambdaSecdTracks, antiLambdaSecdTracks);
      }
    }
  }

  PROCESS_SWITCH(LambdaR2Correlation, processDataReco, "Process for Data and MCReco", true);

  using LambdaMcGenCollisions = aod::LambdaMcGenCollisions;
  using LambdaMcGenTracks = aod::LambdaMcGenTracks;

  SliceCache cachemc;
  Partition<LambdaMcGenTracks> partMcPrimLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);
  Partition<LambdaMcGenTracks> partMcPrimAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);
  Partition<LambdaMcGenTracks> partMcSecdLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kSecondary);
  Partition<LambdaMcGenTracks> partMcSecdAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kSecondary);

  void processMCGen(LambdaMcGenCollisions::iterator const& mcgencol, LambdaMcGenTracks const&)
  {
    histos.fill(HIST("Event/McGen/h1f_collision_posz"), mcgencol.posZ());
    histos.fill(HIST("Event/McGen/h1f_ft0m_mult_percentile"), mcgencol.cent());
    histos.fill(HIST("Event/McGen/h2f_Mult_vs_Centrality"), mcgencol.cent(), mcgencol.mult());

    cent = mcgencol.cent();

    auto lambdaPrimTracks = partMcPrimLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto antiLambdaPrimTracks = partMcPrimAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto lambdaSecdTracks = partMcSecdLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);
    auto antiLambdaSecdTracks = partMcSecdAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cachemc);

    analyzeSingles<kLambda, kGen, kPrimary>(lambdaPrimTracks);
    analyzeSingles<kAntiLambda, kGen, kPrimary>(antiLambdaPrimTracks);

    if (cAnaSecondaries) {
      analyzeSingles<kLambda, kGen, kSecondary>(lambdaSecdTracks);
      analyzeSingles<kAntiLambda, kGen, kSecondary>(antiLambdaSecdTracks);
    }

    if (cAnaPairs) {
      // Primary Pairs Only
      analyzePairs<kLambdaAntiLambda, kGen, kPP, false>(lambdaPrimTracks, antiLambdaPrimTracks);
      analyzePairs<kLambdaLambda, kGen, kPP, true>(lambdaPrimTracks, lambdaPrimTracks);
      analyzePairs<kAntiLambdaAntiLambda, kGen, kPP, true>(antiLambdaPrimTracks, antiLambdaPrimTracks);

      // Secondary Pairs
      if (cAnaSecondaryPairs) {
        analyzePairs<kLambdaAntiLambda, kGen, kPS, false>(lambdaPrimTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaLambda, kGen, kPS, true>(lambdaPrimTracks, lambdaSecdTracks);
        analyzePairs<kAntiLambdaAntiLambda, kGen, kPS, true>(antiLambdaPrimTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaAntiLambda, kGen, kSP, false>(lambdaSecdTracks, antiLambdaPrimTracks);
        analyzePairs<kLambdaLambda, kGen, kSP, true>(lambdaSecdTracks, lambdaPrimTracks);
        analyzePairs<kAntiLambdaAntiLambda, kGen, kSP, true>(antiLambdaSecdTracks, antiLambdaPrimTracks);
        analyzePairs<kLambdaAntiLambda, kGen, kSS, false>(lambdaSecdTracks, antiLambdaSecdTracks);
        analyzePairs<kLambdaLambda, kGen, kSS, true>(lambdaSecdTracks, lambdaSecdTracks);
        analyzePairs<kAntiLambdaAntiLambda, kGen, kSS, true>(antiLambdaSecdTracks, antiLambdaSecdTracks);
      }
    }
  }

  PROCESS_SWITCH(LambdaR2Correlation, processMCGen, "Process for MC Generated", false);
};

struct CascadeSelector {

  //   Zorro zorro;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB;

  Produces<aod::CascadeFlags> cascflags;

  // Configurables
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<bool> useTrigger{"useTrigger", false, "Use trigger selection on skimmed data"};
  Configurable<std::string> triggerList{"triggerList", "fDoubleXi, fDoubleOmega, fOmegaXi", "List of triggers used to select events"};
  Configurable<bool> doTFBorderCut{"doTFBorderCut", true, "Switch to apply TimeframeBorderCut event selection"};
  Configurable<bool> doSel8{"doSel8", true, "Switch to apply sel8 event selection"};
  Configurable<bool> doNoSameBunchPileUp{"doNoSameBunchPileUp", true, "Switch to apply NoSameBunchPileUp event selection"};
  Configurable<int> INEL{"INEL", 0, "Number of charged tracks within |eta| < 1 has to be greater than value"};
  Configurable<double> maxVertexZ{"maxVertexZ", 10., "Maximum value of z coordinate of PV"};
  Configurable<float> etaCascades{"etaCascades", 0.8, "min/max of eta for cascades"};
  Configurable<bool> doCompetingMassCut{"doCompetingMassCut", true, "Switch to apply a competing mass cut for the Omega's"};
  Configurable<float> competingMassWindow{"competingMassWindow", 0.01, "Mass window for the competing mass cut"};

  // Tracklevel
  Configurable<float> tpcNsigmaBachelor{"tpcNsigmaBachelor", 3, "TPC NSigma bachelor"};
  Configurable<float> tpcNsigmaProton{"tpcNsigmaProton", 3, "TPC NSigma proton <- lambda"};
  Configurable<float> tpcNsigmaPion{"tpcNsigmaPion", 3, "TPC NSigma pion <- lambda"};
  Configurable<int> minTPCCrossedRows{"minTPCCrossedRows", 80, "min N TPC crossed rows"}; // TODO: finetune! 80 > 159/2, so no split tracks?
  Configurable<int> minITSClusters{"minITSClusters", 4, "minimum number of ITS clusters"};
  Configurable<float> etaTracks{"etaTracks", 1.0, "min/max of eta for tracks"};
  Configurable<float> tpcChi2{"tpcChi2", 4, "TPC Chi2"};
  Configurable<float> itsChi2{"itsChi2", 36, "ITS Chi2"};

  // Selection criteria - compatible with core wagon autodetect - copied from cascadeanalysis.cxx
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  Configurable<double> v0setting_cospa{"v0setting_cospa", 0.995, "v0setting_cospa"};
  Configurable<float> v0setting_dcav0dau{"v0setting_dcav0dau", 1.0, "v0setting_dcav0dau"};
  Configurable<float> v0setting_dcapostopv{"v0setting_dcapostopv", 0.1, "v0setting_dcapostopv"};
  Configurable<float> v0setting_dcanegtopv{"v0setting_dcanegtopv", 0.1, "v0setting_dcanegtopv"};
  Configurable<float> v0setting_radius{"v0setting_radius", 0.9, "v0setting_radius"};
  Configurable<double> cascadesetting_cospa{"cascadesetting_cospa", 0.95, "cascadesetting_cospa"};
  Configurable<float> cascadesetting_dcacascdau{"cascadesetting_dcacascdau", 1.0, "cascadesetting_dcacascdau"};
  Configurable<float> cascadesetting_dcabachtopv{"cascadesetting_dcabachtopv", 0.05, "cascadesetting_dcabachtopv"};
  Configurable<float> cascadesetting_cascradius{"cascadesetting_cascradius", 0.9, "cascadesetting_cascradius"};
  Configurable<float> cascadesetting_v0masswindow{"cascadesetting_v0masswindow", 0.01, "cascadesetting_v0masswindow"};
  Configurable<float> cascadesetting_mindcav0topv{"cascadesetting_mindcav0topv", 0.01, "cascadesetting_mindcav0topv"};
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  // TODO: variables as function of Omega mass, only do Xi for now
  ConfigurableAxis radiusAxis = {"radiusAxis", {100, 0.0f, 50.0f}, "cm"};
  ConfigurableAxis cpaAxis = {"cpaAxis", {100, 0.95f, 1.0f}, "CPA"};
  ConfigurableAxis vertexAxis = {"vertexAxis", {100, -10.0f, 10.0f}, "cm"};
  ConfigurableAxis dcaAxis = {"dcaAxis", {100, 0.0f, 2.0f}, "cm"};
  ConfigurableAxis invXiMassAxis = {"invXiMassAxis", {100, 1.28f, 1.38f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis invOmegaMassAxis = {"invOmegaMassAxis", {100, 1.62f, 1.72f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis ptAxis = {"ptAxis", {150, 0, 15}, "#it{p}_{T}"};
  ConfigurableAxis rapidityAxis{"rapidityAxis", {100, -1.f, 1.f}, "y"};
  ConfigurableAxis invLambdaMassAxis{"invLambdaMassAxis", {100, 1.07f, 1.17f}, "Inv. Mass (GeV/c^{2})"};
  AxisSpec itsClustersAxis{8, -0.5, 7.5, "number of ITS clusters"};
  AxisSpec tpcRowsAxis{160, -0.5, 159.5, "TPC crossed rows"};
  HistogramRegistry registry{
    "registry",
    {
      // basic selection variables
      {"hV0Radius", "hV0Radius", {HistType::kTH3F, {radiusAxis, invXiMassAxis, ptAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH3F, {radiusAxis, invXiMassAxis, ptAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH3F, {cpaAxis, invXiMassAxis, ptAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH3F, {cpaAxis, invXiMassAxis, ptAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH3F, {vertexAxis, invXiMassAxis, ptAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH3F, {dcaAxis, invXiMassAxis, ptAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH3F, {invLambdaMassAxis, invXiMassAxis, ptAxis}}},

      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {invXiMassAxis, ptAxis, rapidityAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {invXiMassAxis, ptAxis, rapidityAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, rapidityAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {invOmegaMassAxis, ptAxis, rapidityAxis}}},

      // // invariant mass per cut, start with Xi
      // {"hMassXi0", "Xi inv mass before selections", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi1", "Xi inv mass after TPCnCrossedRows cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi2", "Xi inv mass after ITSnClusters cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi3", "Xi inv mass after topo cuts", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi4", "Xi inv mass after V0 daughters PID cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},
      // {"hMassXi5", "Xi inv mass after bachelor PID cut", {HistType::kTH2F, {invMassAxis, ptAxis}}},

      // ITS & TPC clusters, with Xi inv mass
      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH3F, {tpcRowsAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH3F, {itsClustersAxis, invXiMassAxis, ptAxis}}},
      {"hTPCChi2Pos", "hTPCChi2Pos", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Pos"}}}},
      {"hTPCChi2Neg", "hTPCChi2Neg", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Neg"}}}},
      {"hTPCChi2Bach", "hTPCChi2Bach", {HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Bach"}}}},
      {"hITSChi2Pos", "hITSChi2Pos", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Pos"}}}},
      {"hITSChi2Neg", "hITSChi2Neg", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Neg"}}}},
      {"hITSChi2Bach", "hITSChi2Bach", {HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Bach"}}}},

      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},
    },
  };

  // Keep track of which selections the candidates pass
  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);

    auto h = registry.add<TH1>("hSelectionStatus", "hSelectionStatus", HistType::kTH1I, {{10, 0, 10, "status"}});
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "nTPC OK");
    h->GetXaxis()->SetBinLabel(3, "nITS OK");
    h->GetXaxis()->SetBinLabel(4, "track Chi2 OK");
    h->GetXaxis()->SetBinLabel(5, "Topo OK");
    h->GetXaxis()->SetBinLabel(6, "Track eta OK");
    h->GetXaxis()->SetBinLabel(7, "Cascade eta OK");
    h->GetXaxis()->SetBinLabel(8, "V0 PID OK");
    h->GetXaxis()->SetBinLabel(9, "Bach PID OK");

    auto hEventSel = registry.add<TH1>("hEventSel", "hEventSel", HistType::kTH1I, {{10, 0, 10, "selection criteria"}});
    hEventSel->GetXaxis()->SetBinLabel(1, "All");
    hEventSel->GetXaxis()->SetBinLabel(2, "sel8");
    hEventSel->GetXaxis()->SetBinLabel(3, "INEL0");
    hEventSel->GetXaxis()->SetBinLabel(4, "V_z");
    hEventSel->GetXaxis()->SetBinLabel(5, "NoSameBunchPileUp");
    hEventSel->GetXaxis()->SetBinLabel(6, "Selected events");

    if (doprocessRecMC) {
      // only create the rec matched to gen histograms if relevant
      registry.add("truerec/hV0Radius", "hV0Radius", HistType::kTH1F, {radiusAxis});
      registry.add("truerec/hCascRadius", "hCascRadius", HistType::kTH1F, {radiusAxis});
      registry.add("truerec/hV0CosPA", "hV0CosPA", HistType::kTH1F, {cpaAxis});
      registry.add("truerec/hCascCosPA", "hCascCosPA", HistType::kTH1F, {cpaAxis});
      registry.add("truerec/hDCAPosToPV", "hDCAPosToPV", HistType::kTH1F, {vertexAxis});
      registry.add("truerec/hDCANegToPV", "hDCANegToPV", HistType::kTH1F, {vertexAxis});
      registry.add("truerec/hDCABachToPV", "hDCABachToPV", HistType::kTH1F, {vertexAxis});
      registry.add("truerec/hDCAV0ToPV", "hDCAV0ToPV", HistType::kTH1F, {vertexAxis});
      registry.add("truerec/hDCAV0Dau", "hDCAV0Dau", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hDCACascDau", "hDCACascDau", HistType::kTH1F, {dcaAxis});
      registry.add("truerec/hLambdaMass", "hLambdaMass", HistType::kTH1F, {invLambdaMassAxis});
      registry.add("truerec/hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", HistType::kTH1F, {tpcRowsAxis});
      registry.add("truerec/hITSnClustersPos", "hITSnClustersPos", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hITSnClustersNeg", "hITSnClustersNeg", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hITSnClustersBach", "hITSnClustersBach", HistType::kTH1F, {itsClustersAxis});
      registry.add("truerec/hTPCChi2Pos", "hTPCChi2Pos", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Pos"}});
      registry.add("truerec/hTPCChi2Neg", "hTPCChi2Neg", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Neg"}});
      registry.add("truerec/hTPCChi2Bach", "hTPCChi2Bach", HistType::kTH1F, {{100, 0, 10, "TPC Chi2 Bach"}});
      registry.add("truerec/hITSChi2Pos", "hITSChi2Pos", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Pos"}});
      registry.add("truerec/hITSChi2Neg", "hITSChi2Neg", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Neg"}});
      registry.add("truerec/hITSChi2Bach", "hITSChi2Bach", HistType::kTH1F, {{100, 0, 100, "ITS Chi2 Bach"}});
      registry.add("truerec/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("truerec/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("truerec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("truerec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, rapidityAxis});
    }

    if (doprocessGenMC) {
      // only create the MC gen histograms if relevant
      registry.add("gen/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("gen/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("gen/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("gen/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, rapidityAxis});

      registry.add("genwithrec/hXiMinus", "hXiMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("genwithrec/hXiPlus", "hXiPlus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("genwithrec/hOmegaMinus", "hOmegaMinus", HistType::kTH2F, {ptAxis, rapidityAxis});
      registry.add("genwithrec/hOmegaPlus", "hOmegaPlus", HistType::kTH2F, {ptAxis, rapidityAxis});

      registry.add("genwithrec/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events with reconstructed event"}});
      registry.add("gen/hNevents", "hNevents", HistType::kTH1F, {{1, 0, 1, "N generated events"}});
    }
  }

  template <typename TCollision>
  bool eventSelection(TCollision const& collision, bool fillHistos)
  {
    //    if (useTrigger) {
    //      auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    //      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
    //      bool eventTrigger = zorro.isSelected(bc.globalBC());
    //      if (eventTrigger) {
    //        if (fillHistos)
    //          registry.fill(HIST("hTriggerQA"), 1);
    //      } else {
    //        if (fillHistos)
    //          registry.fill(HIST("hTriggerQA"), 0);
    //        return false;
    //      }
    //    }
    // fill event selection based on which selection criteria are applied and passed
    if (fillHistos)
      registry.fill(HIST("hEventSel"), 0);
    if (doSel8 && !collision.sel8()) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 1);
      return false;
    } else if (collision.multNTracksPVeta1() <= INEL) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 2);
      return false;
    } else if (std::abs(collision.posZ()) > maxVertexZ) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 3);
      return false;
    } else if (doNoSameBunchPileUp && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      if (fillHistos)
        registry.fill(HIST("hEventSel"), 4);
      return false;
    }
    // passes all selections
    if (fillHistos)
      registry.fill(HIST("hEventSel"), 5);
    return true;
  }

  template <typename TCollision>
  void fillMatchedHistos(LabeledCascades::iterator rec, int flag, TCollision collision)
  {
    if (flag == 0)
      return;
    if (!rec.has_mcParticle())
      return;
    auto gen = rec.mcParticle();
    if (!gen.isPhysicalPrimary())
      return;
    int genpdg = gen.pdgCode();
    if ((flag < 3 && std::abs(genpdg) == 3312) || (flag > 1 && std::abs(genpdg) == 3334)) {
      // if casc is consistent with Xi and has matched gen Xi OR cand is consistent with Omega and has matched gen omega
      // have to do this in case we reco true Xi with only Omega hypothesis (or vice versa) (very unlikely)
      registry.fill(HIST("truerec/hV0Radius"), rec.v0radius());
      registry.fill(HIST("truerec/hCascRadius"), rec.cascradius());
      registry.fill(HIST("truerec/hV0CosPA"), rec.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hCascCosPA"), rec.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hDCAPosToPV"), rec.dcapostopv());
      registry.fill(HIST("truerec/hDCANegToPV"), rec.dcanegtopv());
      registry.fill(HIST("truerec/hDCABachToPV"), rec.dcabachtopv());
      registry.fill(HIST("truerec/hDCAV0ToPV"), rec.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("truerec/hDCAV0Dau"), rec.dcaV0daughters());
      registry.fill(HIST("truerec/hDCACascDau"), rec.dcacascdaughters());
      registry.fill(HIST("truerec/hLambdaMass"), rec.mLambda());
      registry.fill(HIST("truerec/hITSnClustersPos"), rec.posTrack_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hITSnClustersNeg"), rec.negTrack_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hITSnClustersBach"), rec.bachelor_as<FullTracksExtIUWithPID>().itsNCls());
      registry.fill(HIST("truerec/hTPCnCrossedRowsPos"), rec.posTrack_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hTPCnCrossedRowsNeg"), rec.negTrack_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hTPCnCrossedRowsBach"), rec.bachelor_as<FullTracksExtIUWithPID>().tpcNClsCrossedRows());
      registry.fill(HIST("truerec/hITSChi2Pos"), rec.posTrack_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hITSChi2Neg"), rec.negTrack_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hITSChi2Bach"), rec.bachelor_as<FullTracksExtIUWithPID>().itsChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Pos"), rec.posTrack_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Neg"), rec.negTrack_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      registry.fill(HIST("truerec/hTPCChi2Bach"), rec.bachelor_as<FullTracksExtIUWithPID>().tpcChi2NCl());
      switch (genpdg) { // is matched so we can use genpdg
        case 3312:
          registry.fill(HIST("truerec/hXiMinus"), rec.pt(), rec.yXi());
          break;
        case -3312:
          registry.fill(HIST("truerec/hXiPlus"), rec.pt(), rec.yXi());
          break;
        case 3334:
          registry.fill(HIST("truerec/hOmegaMinus"), rec.pt(), rec.yOmega());
          break;
        case -3334:
          registry.fill(HIST("truerec/hOmegaPlus"), rec.pt(), rec.yOmega());
          break;
      }
    }
  }

  template <typename TCascade, typename TCollision>
  int processCandidate(TCascade const& casc, TCollision const& collision)
  {
    // these are the tracks:
    auto bachTrack = casc.template bachelor_as<FullTracksExtIUWithPID>();
    auto posTrack = casc.template posTrack_as<FullTracksExtIUWithPID>();
    auto negTrack = casc.template negTrack_as<FullTracksExtIUWithPID>();

    // topo variables before cuts:
    registry.fill(HIST("hV0Radius"), casc.v0radius(), casc.mXi(), casc.pt());
    registry.fill(HIST("hCascRadius"), casc.cascradius(), casc.mXi(), casc.pt());
    registry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters(), casc.mXi(), casc.pt());
    registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters(), casc.mXi(), casc.pt());
    registry.fill(HIST("hLambdaMass"), casc.mLambda(), casc.mXi(), casc.pt());

    registry.fill(HIST("hITSnClustersPos"), posTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSnClustersNeg"), negTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSnClustersBach"), bachTrack.itsNCls(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsPos"), posTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsNeg"), negTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hTPCnCrossedRowsBach"), bachTrack.tpcNClsCrossedRows(), casc.mXi(), casc.pt());
    registry.fill(HIST("hITSChi2Pos"), posTrack.itsChi2NCl());
    registry.fill(HIST("hITSChi2Neg"), negTrack.itsChi2NCl());
    registry.fill(HIST("hITSChi2Bach"), bachTrack.itsChi2NCl());
    registry.fill(HIST("hTPCChi2Pos"), posTrack.tpcChi2NCl());
    registry.fill(HIST("hTPCChi2Neg"), negTrack.tpcChi2NCl());
    registry.fill(HIST("hTPCChi2Bach"), bachTrack.tpcChi2NCl());

    registry.fill(HIST("hSelectionStatus"), 0); // all the cascade before selections
    // registry.fill(HIST("hMassXi0"), casc.mXi(), casc.pt());

    // TPC N crossed rows todo: check if minTPCCrossedRows > 50
    if (posTrack.tpcNClsCrossedRows() < minTPCCrossedRows || negTrack.tpcNClsCrossedRows() < minTPCCrossedRows || bachTrack.tpcNClsCrossedRows() < minTPCCrossedRows)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 1); // passes nTPC crossed rows
    // registry.fill(HIST("hMassXi1"), casc.mXi(), casc.pt());

    // ITS N clusters todo: check if minITSClusters > 0
    if (posTrack.itsNCls() < minITSClusters || negTrack.itsNCls() < minITSClusters || bachTrack.itsNCls() < minITSClusters)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 2); // passes nITS clusters
    // registry.fill(HIST("hMassXi2"), casc.mXi(), casc.pt());

    // Chi2 cuts
    if (posTrack.itsChi2NCl() > itsChi2 || negTrack.itsChi2NCl() > itsChi2 || bachTrack.itsChi2NCl() > itsChi2)
      return 0;
    if (posTrack.tpcChi2NCl() > tpcChi2 || negTrack.tpcChi2NCl() > tpcChi2 || bachTrack.tpcChi2NCl() > tpcChi2)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 3); // passes Chi2 cuts

    //// TOPO CUTS //// TODO: improve!
    double pvx = collision.posX();
    double pvy = collision.posY();
    double pvz = collision.posZ();
    if (casc.v0radius() < v0setting_radius ||
        casc.cascradius() < cascadesetting_cascradius ||
        casc.v0cosPA(pvx, pvy, pvz) < v0setting_cospa ||
        casc.casccosPA(pvx, pvy, pvz) < cascadesetting_cospa ||
        casc.dcav0topv(pvx, pvy, pvz) < cascadesetting_mindcav0topv ||
        std::abs(casc.mLambda() - 1.115683) > cascadesetting_v0masswindow)
      return 0; // It failed at least one topo selection

    registry.fill(HIST("hSelectionStatus"), 4); // passes topo
    // registry.fill(HIST("hMassXi3"), casc.mXi(), casc.pt());

    if (std::abs(posTrack.eta()) > etaTracks || std::abs(negTrack.eta()) > etaTracks || std::abs(bachTrack.eta()) > etaTracks)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 5); // passes track eta

    if (std::abs(casc.eta()) > etaCascades)
      return 0;

    registry.fill(HIST("hSelectionStatus"), 6); // passes candidate eta

    // TODO: TOF (for pT > 2 GeV per track?)

    //// TPC PID ////
    // Lambda check
    if (casc.sign() < 0) {
      // Proton check:
      if (std::abs(posTrack.tpcNSigmaPr()) > tpcNsigmaProton)
        return 0;
      // Pion check:
      if (std::abs(negTrack.tpcNSigmaPi()) > tpcNsigmaPion)
        return 0;
    } else {
      // Proton check:
      if (std::abs(negTrack.tpcNSigmaPr()) > tpcNsigmaProton)
        return 0;
      // Pion check:
      if (std::abs(posTrack.tpcNSigmaPi()) > tpcNsigmaPion)
        return 0;
    }
    registry.fill(HIST("hSelectionStatus"), 7); // passes V0 daughters PID
    // registry.fill(HIST("hMassXi4"), casc.mXi(), casc.pt());

    // setting selection flag based on bachelor PID (and competing mass cut for omega's)
    int flag = 0;
    if (std::abs(bachTrack.tpcNSigmaPi()) < tpcNsigmaBachelor)
      flag = 1;
    if (std::abs(bachTrack.tpcNSigmaKa()) < tpcNsigmaBachelor && (!doCompetingMassCut || std::abs(pdgDB->Mass(3312) - casc.mXi()) > competingMassWindow))
      flag = 3 - flag; // 3 if only consistent with omega, 2 if consistent with both

    switch (flag) {
      case 1:                                       // only Xi
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
        }
        break;
      case 2:                                       // Xi or Omega
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
        }
        break;
      case 3:                                       // only Omega
        registry.fill(HIST("hSelectionStatus"), 8); // passes bach PID
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
        }
        break;
    }

    return flag;

  } // processCandidate

  void processGenMC(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, MyCollisions>> const& collisions, aod::McParticles const& mcParticles)
  {
    // evsel
    if (INEL >= 0 && !pwglf::isINELgtNmc(mcParticles, INEL, pdgDB))
      return;
    if (std::abs(mcCollision.posZ()) > maxVertexZ)
      return;

    registry.fill(HIST("gen/hNevents"), 0);

    for (auto const& mcPart : mcParticles) {
      if (!mcPart.isPhysicalPrimary())
        continue;
      if (std::abs(mcPart.eta()) > etaCascades)
        continue;

      switch (mcPart.pdgCode()) {
        case 3312:
          registry.fill(HIST("gen/hXiMinus"), mcPart.pt(), mcPart.y());
          break;
        case -3312:
          registry.fill(HIST("gen/hXiPlus"), mcPart.pt(), mcPart.y());
          break;
        case 3334:
          registry.fill(HIST("gen/hOmegaMinus"), mcPart.pt(), mcPart.y());
          break;
        case -3334:
          registry.fill(HIST("gen/hOmegaPlus"), mcPart.pt(), mcPart.y());
          break;
      }
    }

    // Do the same thing, but now making sure there is at least one matched reconstructed event:
    if (collisions.size() < 1) {
      return;
    } else {
      bool evSel = false; // will be true if at least one rec. collision passes evsel
      for (auto const& collision : collisions) {
        // can be more than 1 rec. collisions due to event splitting
        evSel = eventSelection(collision, false);
        if (evSel) // exit loop if we find 1 rec. event that passes evsel
          break;
      }
      if (evSel) {
        // N gen events with a reconstructed event
        registry.fill(HIST("genwithrec/hNevents"), 0);

        for (auto const& mcPart : mcParticles) {
          if (!mcPart.isPhysicalPrimary())
            continue;
          if (std::abs(mcPart.eta()) > etaCascades)
            continue;

          switch (mcPart.pdgCode()) {
            case 3312:
              registry.fill(HIST("genwithrec/hXiMinus"), mcPart.pt(), mcPart.y());
              break;
            case -3312:
              registry.fill(HIST("genwithrec/hXiPlus"), mcPart.pt(), mcPart.y());
              break;
            case 3334:
              registry.fill(HIST("genwithrec/hOmegaMinus"), mcPart.pt(), mcPart.y());
              break;
            case -3334:
              registry.fill(HIST("genwithrec/hOmegaPlus"), mcPart.pt(), mcPart.y());
              break;
          }
        }
      }
    }
  } // processGen

  // wrappers for data/MC processes on reco level
  void processRecData(MyCollisions::iterator const& collision, aod::CascDataExt const& Cascades, FullTracksExtIUWithPID const&, aod::BCsWithTimestamps const&)
  {
    bool evSel = eventSelection(collision, true);
    // do not skip the collision if event selection fails - this will lead to the cascadeFlag table having less entries than the Cascade table, and therefor not joinable.
    for (auto const& casc : Cascades) {
      if (!evSel) {
        cascflags(0);
        continue;
      }
      int flag = processCandidate(casc, collision);
      cascflags(flag);
    }
  }

  void processRecMC(MyCollisions::iterator const& collision, LabeledCascades const& Cascades, FullTracksExtIUWithPID const&, aod::BCsWithTimestamps const&, aod::McParticles const&)
  {
    bool evSel = eventSelection(collision, true);
    // do not skip the collision if event selection fails - this will lead to the cascadeFlag table having less entries than the Cascade table, and therefor not joinable.
    for (auto const& casc : Cascades) {
      if (!evSel) {
        cascflags(0);
        continue;
      }
      int flag = processCandidate(casc, collision);
      cascflags(flag);
      // do mc matching here
      fillMatchedHistos(casc, flag, collision); // if sign < 0 then pdg > 0
    }
  }

  PROCESS_SWITCH(CascadeSelector, processRecData, "Process rec data", true);
  PROCESS_SWITCH(CascadeSelector, processRecMC, "Process rec MC", false);
  PROCESS_SWITCH(CascadeSelector, processGenMC, "Process gen MC", false);
}; // struct

struct CascadeCorrelations {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Configurables
  Configurable<float> maxRapidity{"maxRapidity", 0.5, "|y| < maxRapidity"};
  Configurable<float> zVertexCut{"zVertexCut", 10, "Cut on PV position"};
  Configurable<int> nMixedEvents{"nMixedEvents", 10, "Number of events to be mixed"};
  Configurable<bool> doEfficiencyCorrection{"doEfficiencyCorrection", true, "flag to do efficiency corrections"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "CCDB url"};
  Configurable<bool> useTrigger{"useTrigger", false, "Use trigger selection on skimmed data"};
  Configurable<std::string> triggerList{"triggerList", "fDoubleXi, fDoubleOmega, fOmegaXi", "List of triggers used to select events"};
  Configurable<std::string> efficiencyCCDBPath{"efficiencyCCDBPath", "Users/r/rspijker/test/EffTest", "Path of the efficiency corrections"};
  Configurable<bool> doTFBorderCut{"doTFBorderCut", true, "Switch to apply TimeframeBorderCut event selection"};
  Configurable<bool> doSel8{"doSel8", true, "Switch to apply sel8 event selection"};

  ConfigurableAxis radiusAxis = {"radiusAxis", {100, 0.0f, 50.0f}, "cm"};
  ConfigurableAxis cpaAxis = {"cpaAxis", {100, 0.95f, 1.0f}, "CPA"};
  ConfigurableAxis invMassAxis = {"invMassAxis", {1000, 1.0f, 2.0f}, "Inv. Mass (GeV/c^{2})"};
  ConfigurableAxis deltaPhiAxis = {"deltaPhiAxis", {180, -PIHalf, 3 * PIHalf}, "#Delta#varphi"}; // 180 is divisible by 18 (tpc sectors) and 20 (run 2 binning)
  ConfigurableAxis ptAxis = {"ptAxis", {150, 0, 15}, "#it{p}_{T}"};
  ConfigurableAxis vertexAxis = {"vertexAxis", {200, -10.0f, 10.0f}, "cm"};
  ConfigurableAxis dcaAxis = {"dcaAxis", {100, 0.0f, 2.0f}, "cm"};
  ConfigurableAxis multiplicityAxis{"multiplicityAxis", {100, 0, 100}, "Multiplicity (centFT0M?)"};
  ConfigurableAxis invLambdaMassAxis{"invLambdaMassAxis", {100, 1.07f, 1.17f}, "Inv. Mass (GeV/c^{2})"};
  AxisSpec signAxis{3, -1.5, 1.5, "sign of cascade"};
  AxisSpec deltaYAxis{40, -2.f, 2.f, "#Delta y"};
  AxisSpec rapidityAxis{100, -1.f, 1.f, "y"};
  AxisSpec selectionFlagAxis{4, -0.5f, 3.5f, "Selection flag of casc candidate"};
  AxisSpec itsClustersAxis{8, -0.5, 7.5, "number of ITS clusters"};
  AxisSpec tpcRowsAxis{160, -0.5, 159.5, "TPC crossed rows"};

  // initialize efficiency maps
  TH1D* hEffXiMin;
  TH1D* hEffXiPlus;
  TH1D* hEffOmegaMin;
  TH1D* hEffOmegaPlus;

  // used in MC closure test
  Service<o2::framework::O2DatabasePDG> pdgDB;
  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    if (doEfficiencyCorrection) {
      TList* effList = ccdb->getForTimeStamp<TList>(efficiencyCCDBPath, 1);
      if (!effList) {
        LOGF(fatal, "null ptr in efficiency list!");
      }
      hEffXiMin = static_cast<TH1D*>(effList->FindObject("hXiMinEff"));
      hEffXiPlus = static_cast<TH1D*>(effList->FindObject("hXiPlusEff"));
      hEffOmegaMin = static_cast<TH1D*>(effList->FindObject("hOmegaMinEff"));
      hEffOmegaPlus = static_cast<TH1D*>(effList->FindObject("hOmegaPlusEff"));
    }

    // zorroSummary.setObject(zorro.getZorroSummary());

    mCounter.mPdgDatabase = pdgDB.service;
    mCounter.mSelectPrimaries = true;
  }

  double getEfficiency(TH1* h, double pT, double y = 0)
  {
    // This function returns 1 / eff
    double eff = h->GetBinContent(h->FindFixBin(pT, y));
    if (eff == 0)
      return 0;
    else
      return 1. / eff;
  }

  bool autoCorrelation(std::array<int, 3> triggerTracks, std::array<int, 3> assocTracks)
  {
    // function that loops over 2 arrays of track indices, checking for common elements
    for (int triggerTrack : triggerTracks) {
      for (int assocTrack : assocTracks) {
        if (triggerTrack == assocTrack)
          return true;
      }
    }
    return false;
  }

  HistogramRegistry registry{
    "registry",
    {
      // inv mass
      {"hMassXiMinus", "hMassXiMinus", {HistType::kTH3F, {{200, 1.24, 1.44, "Inv. Mass (GeV/c^{2})"}, ptAxis, rapidityAxis}}},
      {"hMassXiPlus", "hMassXiPlus", {HistType::kTH3F, {{200, 1.24, 1.44, "Inv. Mass (GeV/c^{2})"}, ptAxis, rapidityAxis}}},
      {"hMassOmegaMinus", "hMassOmegaMinus", {HistType::kTH3F, {{200, 1.6, 1.8, "Inv. Mass (GeV/c^{2})"}, ptAxis, rapidityAxis}}},
      {"hMassOmegaPlus", "hMassOmegaPlus", {HistType::kTH3F, {{200, 1.6, 1.8, "Inv. Mass (GeV/c^{2})"}, ptAxis, rapidityAxis}}},
      // efficiency corrected inv mass
      {"hMassXiEffCorrected", "hMassXiEffCorrected", {HistType::kTHnSparseF, {invMassAxis, signAxis, ptAxis, rapidityAxis, vertexAxis, multiplicityAxis}}, true},
      {"hMassOmegaEffCorrected", "hMassOmegaEffCorrected", {HistType::kTHnSparseF, {invMassAxis, signAxis, ptAxis, rapidityAxis, vertexAxis, multiplicityAxis}}, true},

      // trigger QA
      {"hTriggerQA", "hTriggerQA", {HistType::kTH1F, {{2, -0.5, 1.5, "Trigger y/n"}}}},

      // basic selection variables (after cuts)
      {"hV0Radius", "hV0Radius", {HistType::kTH1F, {radiusAxis}}},
      {"hCascRadius", "hCascRadius", {HistType::kTH1F, {radiusAxis}}},
      {"hV0CosPA", "hV0CosPA", {HistType::kTH1F, {cpaAxis}}},
      {"hCascCosPA", "hCascCosPA", {HistType::kTH1F, {cpaAxis}}},
      {"hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {vertexAxis}}},
      {"hDCAV0Dau", "hDCAV0Dau", {HistType::kTH1F, {dcaAxis}}},
      {"hDCACascDau", "hDCACascDau", {HistType::kTH1F, {dcaAxis}}},
      {"hLambdaMass", "hLambdaMass", {HistType::kTH1F, {invLambdaMassAxis}}},
      {"hTPCnCrossedRowsPos", "hTPCnCrossedRowsPos", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hTPCnCrossedRowsNeg", "hTPCnCrossedRowsNeg", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hTPCnCrossedRowsBach", "hTPCnCrossedRowsBach", {HistType::kTH1F, {tpcRowsAxis}}},
      {"hITSnClustersPos", "hITSnClustersPos", {HistType::kTH1F, {itsClustersAxis}}},
      {"hITSnClustersNeg", "hITSnClustersNeg", {HistType::kTH1F, {itsClustersAxis}}},
      {"hITSnClustersBach", "hITSnClustersBach", {HistType::kTH1F, {itsClustersAxis}}},

      {"hSelectionFlag", "hSelectionFlag", {HistType::kTH1I, {selectionFlagAxis}}},
      {"hAutoCorrelation", "hAutoCorrelation", {HistType::kTH1I, {{4, -0.5f, 3.5f, "Types of SS autocorrelation"}}}},
      {"hAutoCorrelationOS", "hAutoCorrelationOS", {HistType::kTH1I, {{2, -1.f, 1.f, "Charge of OS autocorrelated track"}}}},
      {"hPhi", "hPhi", {HistType::kTH1F, {{180, 0, TwoPI, "#varphi"}}}},
      {"hEta", "hEta", {HistType::kTH1F, {{100, -2, 2, "#eta"}}}},
      {"hRapidityXi", "hRapidityXi", {HistType::kTH1F, {rapidityAxis}}},
      {"hRapidityOmega", "hRapidityOmega", {HistType::kTH1F, {rapidityAxis}}},

      // correlation histos
      {"hDeltaPhiSS", "hDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"hDeltaPhiOS", "hDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},

      {"hXiXi", "hXiXi", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hXiOm", "hXiOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"hOmOm", "hOmOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},

      // Mixed events
      {"MixedEvents/hMEVz1", "hMEVz1", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEVz2", "hMEVz2", {HistType::kTH1F, {vertexAxis}}},
      {"MixedEvents/hMEDeltaPhiSS", "hMEDeltaPhiSS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEDeltaPhiOS", "hMEDeltaPhiOS", {HistType::kTH1F, {deltaPhiAxis}}},
      {"MixedEvents/hMEQA", "hMEQA", {HistType::kTH1I, {{2, 0, 2, "QA for exceptions in ME (this histogram should have 0 entries!)"}}}},
      {"MixedEvents/hMEAutoCorrelation", "hMEAutoCorrelation", {HistType::kTH1I, {{4, -0.5f, 3.5f, "Types of SS autocorrelation"}}}},
      {"MixedEvents/hMEAutoCorrelationOS", "hMEAutoCorrelationOS", {HistType::kTH1I, {{2, -1.f, 1.f, "Charge of OS autocorrelated track"}}}},

      {"MixedEvents/hMEXiXi", "hMEXiXi", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEXiOm", "hMEXiOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},
      {"MixedEvents/hMEOmOm", "hMEOmOm", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, signAxis, signAxis, ptAxis, ptAxis, invMassAxis, invMassAxis, vertexAxis, multiplicityAxis}}, true},

      // MC closure
      {"MC/hMCPlusMinus", "hMCPlusMinus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCPlusPlus", "hMCPlusPlus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCMinusPlus", "hMCMinusPlus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},
      {"MC/hMCMinusMinus", "hMCMinusMinus", {HistType::kTHnSparseF, {deltaPhiAxis, deltaYAxis, ptAxis, ptAxis, vertexAxis, multiplicityAxis}}, true},

      {"MC/hGenMultNoReco", "hGenMultNoReco", {HistType::kTH1I, {{100, 0, 100, "Number of generated charged primaries"}}}},
      {"MC/hGenMultOneReco", "hGenMultOneReco", {HistType::kTH1I, {{100, 0, 100, "Number of generated charged primaries"}}}},
      {"MC/hSplitEvents", "hSplitEvents", {HistType::kTH1I, {{10, 0, 10, "Number of rec. events per gen event"}}}},

      // debug
      {"MC/hPhi", "hPhi", {HistType::kTH1F, {{180, 0, TwoPI}}}},
      {"MC/hEta", "hEta", {HistType::kTH1F, {{100, -2, 2}}}},
      {"MC/hRapidity", "hRapidity", {HistType::kTH1F, {{100, -2, 2}}}},
    },
  };

  // cascade filter
  Filter cascadeSelector = aod::cascadeflags::isSelected > 0;

  // Warning: it is not possible to use this axis as configurable due to a bug - however, default values are sensible.
  SliceCache cache;
  ConfigurableAxis axisVtxZ{"axisVtxZ", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  // ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100, 1000}, "Mixing bins - multiplicity"};
  // using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::centFT0M>;
  // BinningType colBinning{{axisVtxZ, axisMult}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;
  BinningType colBinning{{axisVtxZ}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.

  void processSameEvent(MyCollisionsMult::iterator const& collision, MyCascades const& Cascades, FullTracksExtIU const&, aod::BCsWithTimestamps const&)
  {
    //    if (useTrigger) {
    //      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    //      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerList);
    //      bool eventTrigger = zorro.isSelected(bc.globalBC());
    //      if (eventTrigger) {
    //        registry.fill(HIST("hTriggerQA"), 1);
    //      } else {
    //        registry.fill(HIST("hTriggerQA"), 0);
    //        return;
    //      }
    //    }

    double weight;
    // Some QA on the cascades
    for (auto const& casc : Cascades) {
      if (casc.isSelected() <= 2) { // not exclusively an Omega --> consistent with Xi or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassXiMinus"), casc.mXi(), casc.pt(), casc.yXi());
          weight = getEfficiency(hEffXiMin, casc.pt(), casc.yXi());
        } else {
          registry.fill(HIST("hMassXiPlus"), casc.mXi(), casc.pt(), casc.yXi());
          weight = getEfficiency(hEffXiPlus, casc.pt(), casc.yXi());
        }
        // LOGF(info, "casc pt %f, weight %f", casc.pt(), weight);
        registry.fill(HIST("hMassXiEffCorrected"), casc.mXi(), casc.sign(), casc.pt(), casc.yXi(), collision.posZ(), collision.centFT0M(), weight);
        registry.fill(HIST("hRapidityXi"), casc.yXi());
      }
      if (casc.isSelected() >= 2) { // consistent with Omega or both
        if (casc.sign() < 0) {
          registry.fill(HIST("hMassOmegaMinus"), casc.mOmega(), casc.pt(), casc.yOmega());
          weight = getEfficiency(hEffOmegaMin, casc.pt(), casc.yOmega());
        } else {
          registry.fill(HIST("hMassOmegaPlus"), casc.mOmega(), casc.pt(), casc.yOmega());
          weight = getEfficiency(hEffOmegaPlus, casc.pt(), casc.yOmega());
        }
        registry.fill(HIST("hMassOmegaEffCorrected"), casc.mOmega(), casc.sign(), casc.pt(), casc.yOmega(), collision.posZ(), collision.centFT0M(), weight);
        registry.fill(HIST("hRapidityOmega"), casc.yOmega());
      }
      registry.fill(HIST("hV0Radius"), casc.v0radius());
      registry.fill(HIST("hCascRadius"), casc.cascradius());
      registry.fill(HIST("hV0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hCascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAPosToPV"), casc.dcapostopv());
      registry.fill(HIST("hDCANegToPV"), casc.dcanegtopv());
      registry.fill(HIST("hDCABachToPV"), casc.dcabachtopv());
      registry.fill(HIST("hDCAV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAV0Dau"), casc.dcaV0daughters());
      registry.fill(HIST("hDCACascDau"), casc.dcacascdaughters());
      registry.fill(HIST("hLambdaMass"), casc.mLambda());
      registry.fill(HIST("hITSnClustersPos"), casc.posTrack_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hITSnClustersNeg"), casc.negTrack_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hITSnClustersBach"), casc.bachelor_as<FullTracksExtIU>().itsNCls());
      registry.fill(HIST("hTPCnCrossedRowsPos"), casc.posTrack_as<FullTracksExtIU>().tpcNClsCrossedRows());
      registry.fill(HIST("hTPCnCrossedRowsNeg"), casc.negTrack_as<FullTracksExtIU>().tpcNClsCrossedRows());
      registry.fill(HIST("hTPCnCrossedRowsBach"), casc.bachelor_as<FullTracksExtIU>().tpcNClsCrossedRows());

      registry.fill(HIST("hSelectionFlag"), casc.isSelected());
      registry.fill(HIST("hPhi"), casc.phi());
      registry.fill(HIST("hEta"), casc.eta());
    } // casc loop

    for (auto& [c0, c1] : combinations(Cascades, Cascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT. As we can't swap the cascade tables themselves, we swap the addresses and later dereference them
      auto* triggerAddress = &c0;
      auto* assocAddress = &c1;
      if (assocAddress->pt() > triggerAddress->pt()) {
        std::swap(triggerAddress, assocAddress);
      }
      auto trigger = *triggerAddress;
      auto assoc = *assocAddress;

      // autocorrelation check
      std::array<int, 3> triggerTracks = {trigger.posTrackId(), trigger.negTrackId(), trigger.bachelorId()};
      std::array<int, 3> assocTracks = {assoc.posTrackId(), assoc.negTrackId(), assoc.bachelorId()};
      if (autoCorrelation(triggerTracks, assocTracks))
        continue;

      // calculate angular correlations
      double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

      double invMassXiTrigg = trigger.mXi();
      double invMassOmTrigg = trigger.mOmega();
      double invMassXiAssoc = assoc.mXi();
      double invMassOmAssoc = assoc.mOmega();

      double weightTrigg = 1.;
      double weightAssoc = 1.;

      if (trigger.isSelected() <= 2 && std::abs(trigger.yXi()) < maxRapidity) { // trigger Xi
        if (doEfficiencyCorrection)
          weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffXiMin, trigger.pt()) : getEfficiency(hEffXiPlus, trigger.pt());
        if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
          if (doEfficiencyCorrection)
            weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt()) : getEfficiency(hEffXiPlus, assoc.pt());
          registry.fill(HIST("hXiXi"), dphi, trigger.yXi() - assoc.yXi(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
        }
        if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
          if (doEfficiencyCorrection)
            weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt()) : getEfficiency(hEffOmegaPlus, assoc.pt());
          registry.fill(HIST("hXiOm"), dphi, trigger.yXi() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
        }
      }
      if (trigger.isSelected() >= 2 && std::abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
        if (doEfficiencyCorrection)
          weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffOmegaMin, trigger.pt()) : getEfficiency(hEffOmegaPlus, trigger.pt());
        if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
          if (doEfficiencyCorrection)
            weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt()) : getEfficiency(hEffXiPlus, assoc.pt());
          // if Omega-Xi, fill the Xi-Omega histogram (flip the trigger/assoc and dphy,dy signs)
          registry.fill(HIST("hXiOm"), RecoDecay::constrainAngle(assoc.phi() - trigger.phi(), -PIHalf), -(trigger.yOmega() - assoc.yXi()), assoc.sign(), trigger.sign(), assoc.pt(), trigger.pt(), invMassXiAssoc, invMassOmTrigg, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
        }
        if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
          if (doEfficiencyCorrection)
            weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt()) : getEfficiency(hEffOmegaPlus, assoc.pt());
          registry.fill(HIST("hOmOm"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, collision.posZ(), collision.centFT0M(), weightTrigg * weightAssoc);
        }
      }

      // QA plots
      if (trigger.sign() * assoc.sign() < 0) {
        registry.fill(HIST("hDeltaPhiOS"), dphi);
      } else {
        registry.fill(HIST("hDeltaPhiSS"), dphi);
      }
    } // correlations
  } // process same event

  void processMixedEvent(MyCollisionsMult const& collisions, MyCascades const& Cascades, FullTracksExtIU const&)
  {
    auto cascadesTuple = std::make_tuple(Cascades);
    SameKindPair<MyCollisionsMult, MyCascades, BinningType> pair{colBinning, nMixedEvents, -1, collisions, cascadesTuple, &cache};

    for (auto const& [col1, cascades1, col2, cascades2] : pair) {
      if (!col1.sel8() || !col2.sel8())
        continue;
      if (std::abs(col1.posZ()) > zVertexCut || std::abs(col2.posZ()) > zVertexCut)
        continue;
      if (col1.globalIndex() == col2.globalIndex()) {
        registry.fill(HIST("hMEQA"), 0.5);
        continue;
      }
      registry.fill(HIST("MixedEvents/hMEVz1"), col1.posZ());
      registry.fill(HIST("MixedEvents/hMEVz2"), col2.posZ());

      for (auto& [casc1, casc2] : combinations(CombinationsFullIndexPolicy(cascades1, cascades2))) {
        // specify FullIndexPolicy since the cascades are from different collisions
        auto* triggerAddress = &casc1;
        auto* assocAddress = &casc2;
        if (assocAddress->pt() > triggerAddress->pt()) {
          std::swap(triggerAddress, assocAddress);
        }
        auto trigger = *triggerAddress;
        auto assoc = *assocAddress;

        if (trigger.collisionId() == assoc.collisionId()) {
          registry.fill(HIST("hMEQA"), 1.5);
          continue;
        }

        std::array<int, 3> triggerTracks = {trigger.posTrackId(), trigger.negTrackId(), trigger.bachelorId()};
        std::array<int, 3> assocTracks = {assoc.posTrackId(), assoc.negTrackId(), assoc.bachelorId()};
        if (autoCorrelation(triggerTracks, assocTracks))
          continue;

        double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

        double invMassXiTrigg = trigger.mXi();
        double invMassOmTrigg = trigger.mOmega();
        double invMassXiAssoc = assoc.mXi();
        double invMassOmAssoc = assoc.mOmega();

        double weightTrigg = 1.;
        double weightAssoc = 1.;

        if (trigger.isSelected() <= 2 && std::abs(trigger.yXi()) < maxRapidity) { // trigger Xi
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffXiMin, trigger.pt()) : getEfficiency(hEffXiPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt()) : getEfficiency(hEffXiPlus, assoc.pt());
            registry.fill(HIST("MixedEvents/hMEXiXi"), dphi, trigger.yXi() - assoc.yXi(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassXiAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt()) : getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("MixedEvents/hMEXiOm"), dphi, trigger.yXi() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassXiTrigg, invMassOmAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
          }
        }
        if (trigger.isSelected() >= 2 && std::abs(trigger.yOmega()) < maxRapidity) { // trigger Omega
          if (doEfficiencyCorrection)
            weightTrigg = trigger.sign() < 0 ? getEfficiency(hEffOmegaMin, trigger.pt()) : getEfficiency(hEffOmegaPlus, trigger.pt());
          if (assoc.isSelected() <= 2 && std::abs(assoc.yXi()) < maxRapidity) { // assoc Xi
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffXiMin, assoc.pt()) : getEfficiency(hEffXiPlus, assoc.pt());
            // if Omega-Xi, fill the Xi-Omega histogram (flip the trigger/assoc and dphy,dy signs)
            registry.fill(HIST("MixedEvents/hMEXiOm"), RecoDecay::constrainAngle(assoc.phi() - trigger.phi(), -PIHalf), -(trigger.yOmega() - assoc.yXi()), assoc.sign(), trigger.sign(), assoc.pt(), trigger.pt(), invMassXiAssoc, invMassOmTrigg, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
          }
          if (assoc.isSelected() >= 2 && std::abs(assoc.yOmega()) < maxRapidity) { // assoc Omega
            if (doEfficiencyCorrection)
              weightAssoc = assoc.sign() < 0 ? getEfficiency(hEffOmegaMin, assoc.pt()) : getEfficiency(hEffOmegaPlus, assoc.pt());
            registry.fill(HIST("MixedEvents/hMEOmOm"), dphi, trigger.yOmega() - assoc.yOmega(), trigger.sign(), assoc.sign(), trigger.pt(), assoc.pt(), invMassOmTrigg, invMassOmAssoc, col1.posZ(), col1.centFT0M(), weightTrigg * weightAssoc);
          }
        }

        // QA plots
        if (trigger.sign() * assoc.sign() < 0) {
          registry.fill(HIST("MixedEvents/hMEDeltaPhiOS"), dphi);
        } else {
          registry.fill(HIST("MixedEvents/hMEDeltaPhiSS"), dphi);
        }
      } // correlations
    } // collisions
  } // process mixed events

  Configurable<float> etaGenCascades{"etaGenCascades", 0.8, "min/max of eta for generated cascades"};
  Filter genCascadesFilter = nabs(aod::mcparticle::pdgCode) == 3312;

  void processMC(aod::McCollision const&, soa::SmallGroups<soa::Join<aod::McCollisionLabels, MyCollisionsMult>> const& collisions, soa::Filtered<aod::McParticles> const& genCascades, aod::McParticles const& mcParticles)
  {
    // Let's do some logic on matched reconstructed collisions - if there less or more than one, fill some QA and skip the rest
    double FT0mult = -1; // non-sensible default value just in case
    double vtxz = -999.; // non-sensible default value just in case
    if (collisions.size() < 1) {
      registry.fill(HIST("MC/hSplitEvents"), 0);
      registry.fill(HIST("MC/hGenMultNoReco"), mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles));
      return;
    } else if (collisions.size() == 1) {
      registry.fill(HIST("MC/hSplitEvents"), 1);
      registry.fill(HIST("MC/hGenMultOneReco"), mCounter.countFT0A(mcParticles) + mCounter.countFT0C(mcParticles));
      for (auto const& collision : collisions) { // not really a loop, as there is only one collision
        FT0mult = collision.centFT0M();
        vtxz = collision.posZ();
      }
    } else if (collisions.size() > 1) {
      registry.fill(HIST("MC/hSplitEvents"), collisions.size());
      return;
    }

    // QA
    for (auto& casc : genCascades) {
      if (!casc.isPhysicalPrimary())
        continue;
      registry.fill(HIST("MC/hPhi"), casc.phi());
      registry.fill(HIST("MC/hEta"), casc.eta());
      registry.fill(HIST("MC/hRapidity"), casc.y());
    }

    for (auto& [c0, c1] : combinations(genCascades, genCascades)) { // combinations automatically applies strictly upper in case of 2 identical tables
      // Define the trigger as the particle with the highest pT. As we can't swap the cascade tables themselves, we swap the addresses and later dereference them
      auto* triggerAddress = &c0;
      auto* assocAddress = &c1;
      if (assocAddress->pt() > triggerAddress->pt()) {
        std::swap(triggerAddress, assocAddress);
      }
      auto trigger = *triggerAddress;
      auto assoc = *assocAddress;

      if (!trigger.isPhysicalPrimary() || !assoc.isPhysicalPrimary())
        continue; // require the cascades to be primaries
      if (std::abs(trigger.eta()) > etaGenCascades)
        continue; // only apply eta cut to trigger - trigger normalization still valid without introducing 2-particle-acceptance effects

      double dphi = RecoDecay::constrainAngle(trigger.phi() - assoc.phi(), -PIHalf);

      if (trigger.pdgCode() < 0) { // anti-trigg --> Plus
        if (assoc.pdgCode() < 0) { // anti-assoc --> Plus
          registry.fill(HIST("MC/hMCPlusPlus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
        } else { // assoc --> Minus
          registry.fill(HIST("MC/hMCPlusMinus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
        }
      } else {                     // trig --> Minus
        if (assoc.pdgCode() < 0) { // anti-assoc --> Plus
          registry.fill(HIST("MC/hMCMinusPlus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
        } else {
          registry.fill(HIST("MC/hMCMinusMinus"), dphi, trigger.y() - assoc.y(), trigger.pt(), assoc.pt(), vtxz, FT0mult);
        }
      }
    }
  }

  PROCESS_SWITCH(CascadeCorrelations, processSameEvent, "Process same events", true);
  PROCESS_SWITCH(CascadeCorrelations, processMixedEvent, "Process mixed events", true);
  PROCESS_SWITCH(CascadeCorrelations, processMC, "Process MC", false);

}; // struct


struct LambdaXiCorrelation {

  // --- Configurables ---
  Configurable<float> maxY{"maxY", 0.5, "Max |y| for Lambda and Xi"};
  Configurable<bool> useEff{"useEff", false, "Apply Lambda efficiency correction"};

  // --- Outputs ---
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // --- Data Slicing Definitions ---
  using GoodLambdas = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;
  Partition<GoodLambdas> goodLambda = aod::lambdatrackext::trueLambdaFlag == true;

  Preslice<GoodLambdas> lambdasPerCollision = aod::lambdatrack::lambdaCollisionId;
  Preslice<aod::CascDataExt> cascadesPerCollision = aod::cascdata::collisionId;

  using LambdaCollisionsExt = aod::LambdaCollisions;

  // --- R2 Calculation Helper ---
  // R2 = (N_events * Pair_Yield) / (Single_Yield_1 * Single_Yield_2) - 1
  static TH2* calculateR2(TH2* hPairs, TH1* hSinglesTrig, TH1* hSinglesAssoc, double nEvents)
  {
    if (!hPairs || !hSinglesTrig || !hSinglesAssoc || nEvents <= 0)
      return nullptr;

    TH2* hR2 = reinterpret_cast<TH2*>(hPairs->Clone(Form("%s_R2", hPairs->GetName())));
    hR2->Reset();

    double nS1 = hSinglesTrig->Integral();
    double nS2 = hSinglesAssoc->Integral();

    if (nS1 > 0 && nS2 > 0) {
      hR2->Add(hPairs);
      hR2->Scale(nEvents / (nS1 * nS2));

      for (int i = 1; i <= hR2->GetNbinsX(); i++) {
        for (int j = 1; j <= hR2->GetNbinsY(); j++) {
          double content = hR2->GetBinContent(i, j);
          hR2->SetBinContent(i, j, content - 1.0);
        }
      }
    }
    return hR2;
  }

  void init(InitContext const&)
  {
    // --- 1. Axis Definitions ---
    const AxisSpec dphi{72, -PIHalf, 3 * PIHalf, "#Delta#varphi"};
    const AxisSpec dy{40, -2.0f, 2.0f, "#Delta y"};
    const AxisSpec pt{100, 0, 10, "p_{T} (GeV/c)"};
    const AxisSpec rap{100, -1.0, 1.0, "y"};
    const AxisSpec massLam{100, 1.09, 1.14, "M_{p#pi} (GeV/c^{2})"};
    const AxisSpec massXi{100, 1.28, 1.36, "M_{#Lambda#pi} (GeV/c^{2})"};
    const AxisSpec radius{100, 0, 100, "Radius (cm)"};
    const AxisSpec cpa{100, 0.9, 1.0, "Cos(PA)"};
    const AxisSpec dca{100, 0.0, 5.0, "DCA (cm)"};
    const AxisSpec pvDca{100, -10.0, 10.0, "DCA to PV (cm)"};

    // --- 2. Histograms ---
    histos.add("Event/hEventCount", "Event Counter", kTH1F, {{1, 0, 1, "Count"}});

    // Singles: Lambda
    histos.add("Singles/Lambda/hPt", "Lambda p_{T}", kTH1F, {pt});
    histos.add("Singles/AntiLambda/hPt", "AntiLambda p_{T}", kTH1F, {pt});

    histos.add("Singles/Lambda/hPtVsMass", "Lambda p_{T} vs Mass", kTH2F, {massLam, pt});
    histos.add("Singles/AntiLambda/hPtVsMass", "AntiLambda p_{T} vs Mass", kTH2F, {massLam, pt});

    // Singles: Xi & QA
    histos.add("QA/Xi/hRadius", "Xi Radius", kTH1F, {radius});
    histos.add("QA/Xi/hCosPA", "Xi CosPA", kTH1F, {cpa});
    histos.add("QA/Xi/hDCAV0Dau", "DCA V0 Daughters", kTH1F, {dca});
    histos.add("QA/Xi/hDCACascDau", "DCA Casc Daughters", kTH1F, {dca});
    histos.add("QA/Xi/hDCAV0ToPV", "DCA V0 to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCAPosToPV", "DCA Pos to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCANegToPV", "DCA Neg to PV", kTH1F, {pvDca});
    histos.add("QA/Xi/hDCABachToPV", "DCA Bach to PV", kTH1F, {pvDca});

    histos.add("Singles/XiMinus/hPtVsMass", "Xi^{-} p_{T} vs Mass", kTH2F, {massXi, pt});
    histos.add("Singles/XiPlus/hPtVsMass", "Xi^{+} p_{T} vs Mass", kTH2F, {massXi, pt});
    histos.add("Singles/XiMinus/hRap", "Xi^{-} Rapidity", kTH1F, {rap});
    histos.add("Singles/XiPlus/hRap", "Xi^{+} Rapidity", kTH1F, {rap});

    // Pairs: Charge Separated (R2 Inputs)
    histos.add("Pairs/Lam_XiM/hDeltaPhiDeltaY", "L-Xi-", kTH2F, {dphi, dy});
    histos.add("Pairs/Lam_XiP/hDeltaPhiDeltaY", "L-Xi+", kTH2F, {dphi, dy});
    histos.add("Pairs/AntiLam_XiM/hDeltaPhiDeltaY", "AL-Xi-", kTH2F, {dphi, dy});
    histos.add("Pairs/AntiLam_XiP/hDeltaPhiDeltaY", "AL-Xi+", kTH2F, {dphi, dy});
  }

  // --- Analysis Functions ---

  template <typename T>
  void analyzeSinglesLambda(T const& tracks)
  {
    for (const auto& track : tracks) {
      if (std::abs(track.rap()) > maxY)
        continue;

      float w = useEff ? track.corrFact() : 1.0f;
      bool isAnti = (track.v0Type() == 1);

      if (!isAnti) {
        histos.fill(HIST("Singles/Lambda/hPt"), track.pt(), w);
        histos.fill(HIST("Singles/Lambda/hPtVsMass"), track.mass(), track.pt(), w);
      } else {
        histos.fill(HIST("Singles/AntiLambda/hPt"), track.pt(), w);
        histos.fill(HIST("Singles/AntiLambda/hPtVsMass"), track.mass(), track.pt(), w);
      }
    }
  }

  template <typename T, typename F>
  void analyzeSinglesXi(T const& cascades, F const& flagsStart, float pvX, float pvY, float pvZ)
  {
    for (const auto& casc : cascades) {
      if ((flagsStart + casc.globalIndex()).isSelected() == 0)
        continue;

      float xiY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, MassXi0);
      if (std::abs(xiY) > maxY)
        continue;

      // QA Filling
      histos.fill(HIST("QA/Xi/hRadius"), casc.cascradius());
      histos.fill(HIST("QA/Xi/hCosPA"), casc.casccosPA(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Xi/hDCAV0Dau"), casc.dcaV0daughters());
      histos.fill(HIST("QA/Xi/hDCACascDau"), casc.dcacascdaughters());
      histos.fill(HIST("QA/Xi/hDCAV0ToPV"), casc.dcav0topv(pvX, pvY, pvZ));
      histos.fill(HIST("QA/Xi/hDCAPosToPV"), casc.dcapostopv());
      histos.fill(HIST("QA/Xi/hDCANegToPV"), casc.dcanegtopv());
      histos.fill(HIST("QA/Xi/hDCABachToPV"), casc.dcabachtopv());

      if (casc.sign() < 0) {
        histos.fill(HIST("Singles/XiMinus/hPtVsMass"), casc.mXi(), casc.pt());
        histos.fill(HIST("Singles/XiMinus/hRap"), xiY);
      } else {
        histos.fill(HIST("Singles/XiPlus/hPtVsMass"), casc.mXi(), casc.pt());
        histos.fill(HIST("Singles/XiPlus/hRap"), xiY);
      }
    }
  }

  template <typename L, typename C, typename F>
  void analyzePairs(L const& lambdas, C const& cascades, F const& flagsStart)
  {
    for (const auto& lam : lambdas) {
      if (std::abs(lam.rap()) > maxY)
        continue;
      float wLam = useEff ? lam.corrFact() : 1.0f;
      bool isAntiLam = (lam.v0Type() == 1);

      for (const auto& casc : cascades) {
        if ((flagsStart + casc.globalIndex()).isSelected() == 0)
          continue;

        float xiY = RecoDecay::y(std::array{casc.px(), casc.py(), casc.pz()}, MassXi0);
        if (std::abs(xiY) > maxY)
          continue;

        float dphi = RecoDecay::constrainAngle(casc.phi() - lam.phi(), -PIHalf);
        float dy = xiY - lam.rap();

        bool isXiPlus = (casc.sign() > 0);

        if (!isAntiLam && !isXiPlus)
          histos.fill(HIST("Pairs/Lam_XiM/hDeltaPhiDeltaY"), dphi, dy, wLam);
        else if (!isAntiLam && isXiPlus)
          histos.fill(HIST("Pairs/Lam_XiP/hDeltaPhiDeltaY"), dphi, dy, wLam);
        else if (isAntiLam && !isXiPlus)
          histos.fill(HIST("Pairs/AntiLam_XiM/hDeltaPhiDeltaY"), dphi, dy, wLam);
        else if (isAntiLam && isXiPlus)
          histos.fill(HIST("Pairs/AntiLam_XiP/hDeltaPhiDeltaY"), dphi, dy, wLam);
      }
    }
  }

  void process(LambdaCollisionsExt::iterator const& lambdacoll,
               GoodLambdas const& /*lambdas*/,
               aod::CascDataExt const& cascades,
               aod::CascadeFlags const& cascflags)
  {
    histos.fill(HIST("Event/hEventCount"), 0.5);

    auto lambdasInThisEvent = goodLambda->sliceBy(lambdasPerCollision, lambdacoll.globalIndex());
    const int64_t refCollisionIndex = lambdacoll.refCollId();
    auto cascadesInThisEvent = cascades.sliceBy(cascadesPerCollision, refCollisionIndex);

    float pvX = lambdacoll.posX();
    float pvY = lambdacoll.posY();
    float pvZ = lambdacoll.posZ();

    auto flagsStart = cascflags.begin();

    analyzeSinglesLambda(lambdasInThisEvent);
    analyzeSinglesXi(cascadesInThisEvent, flagsStart, pvX, pvY, pvZ);
    analyzePairs(lambdasInThisEvent, cascadesInThisEvent, flagsStart);
  }

  PROCESS_SWITCH(LambdaXiCorrelation, process, "Λ–Ξ correlation (Final Complete)", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{

    adaptAnalysisTask<LambdaTableProducer>(cfgc),
    adaptAnalysisTask<LambdaTracksExtProducer>(cfgc),
    adaptAnalysisTask<LambdaR2Correlation>(cfgc),
    adaptAnalysisTask<CascadeSelector>(cfgc),
    adaptAnalysisTask<CascadeCorrelations>(cfgc),
    adaptAnalysisTask<LambdaXiCorrelation>(cfgc)

  };
}
