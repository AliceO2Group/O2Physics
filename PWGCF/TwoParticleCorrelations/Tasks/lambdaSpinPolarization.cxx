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

/// \file lambdaSpinPolarization.cxx
/// \brief Task to study the Lambda spin polarization
///\author Subhadeep Roy <subhadeep.roy@cern.ch>
/// \author Yash Patley <yash.patley@cern.ch>,

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TList.h>
#include <TObject.h>
#include <TPDGCode.h>
#include <TString.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <string_view>
#include <vector>

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
DECLARE_SOA_COLUMN(Mult, mult, float);
DECLARE_SOA_COLUMN(TimeStamp, timeStamp, uint64_t);
// DECALRE_SOA_
} // namespace lambdacollision

DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ,
                  lambdacollision::TimeStamp
                  //, lambdacollision::CollisionId
);
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
DECLARE_SOA_COLUMN(PrPx, prPx, float);
DECLARE_SOA_COLUMN(PrPy, prPy, float);
DECLARE_SOA_COLUMN(PrPz, prPz, float);
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
                  lambdatrack::PrPx,
                  lambdatrack::PrPy,
                  lambdatrack::PrPz,
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
                  lambdatrack::PrPx,
                  lambdatrack::PrPy,
                  lambdatrack::PrPz,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::V0Type,
                  lambdatrack::CosPA,
                  lambdatrack::DcaDau,
                  lambdatrack::V0PrmScd,
                  lambdatrack::CorrFact);
using LambdaMcGenTrack = LambdaMcGenTracks::iterator;

namespace lambdamixeventcollision
{
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int);
} // namespace lambdamixeventcollision

DECLARE_SOA_TABLE(LambdaMixEventCollisions, "AOD", "LAMBDAMIXCOLS", o2::soa::Index<>,
                  lambdamixeventcollision::CollisionIndex,
                  lambdacollision::Cent,
                  aod::collision::PosZ,
                  lambdacollision::TimeStamp);
using LambdaMixEventCollision = LambdaMixEventCollisions::iterator;

namespace lambdamixeventtracks
{
// DECLARE_SOA_INDEX_COLUMN(LambdaMixEventCollision, lambdaMixEventCollision);
DECLARE_SOA_COLUMN(LambdaMixEventCollisionIdx, lambdaMixEventCollisionIdx, int);
DECLARE_SOA_COLUMN(LambdaMixEventTrackIdx, lambdaMixEventTrackIdx, int);
DECLARE_SOA_COLUMN(LambdaMixEventTimeStamp, lambdaMixEventTimeStamp, uint64_t);

} // namespace lambdamixeventtracks

DECLARE_SOA_TABLE(LambdaMixEventTracks, "AOD", "LAMBDAMIXTRKS", o2::soa::Index<>,
                  // lambdamixeventtracks::LambdaMixEventCollisionId,
                  lambdamixeventtracks::LambdaMixEventCollisionIdx,
                  lambdamixeventtracks::LambdaMixEventTrackIdx,
                  lambdatrack::Px,
                  lambdatrack::Py,
                  lambdatrack::Pz,
                  lambdatrack::Mass,
                  lambdatrack::PrPx,
                  lambdatrack::PrPy,
                  lambdatrack::PrPz,
                  lambdatrack::V0Type,
                  lambdamixeventtracks::LambdaMixEventTimeStamp);

using LambdaMixEventTrack = LambdaMixEventTracks::iterator;
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
  kCentFT0C
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
  kAntiLambdaLambda = 1,
  kLambdaLambda = 2,
  kAntiLambdaAntiLambda = 3,

  kLambdaSBAntiLambda = 4,
  kAntiLambdaSBLambda = 5,
  kLambdaSBLambda = 6,
  kAntiLambdaSBAntiLambda = 7,

  kLambdaSBSBAntiLambda = 8,
  kAntiLambdaSBSBLambda = 9,
  kLambdaSBSBLambda = 10,
  kAntiLambdaSBSBAntiLambda = 11
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

  Configurable<int> cCentEstimator{"cCentEstimator", 0, "Centrality Estimator: 0=FT0M, 1=FT0C"};
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ (cm)"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ (cm)"};
  Configurable<float> cMinMult{"cMinMult", 0.0, "Min centrality percentile"};
  Configurable<float> cMaxMult{"cMaxMult", 100.0, "Max centrality percentile"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A+T0C) Run3"};
  Configurable<bool> cInt7Trig{"cInt7Trig", false, "kINT7 MB Trigger"};
  Configurable<bool> cSel7Trig{"cSel7Trig", false, "Sel7 (V0A+V0C) Run2"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "TVX Trigger Selection"};
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
  Configurable<float> cMinTpcCROverCls{"cMinTpcCROverCls", 0.8, "TPC Min CR/Findable Cls"};
  Configurable<float> cMaxTpcSharedClusters{"cMaxTpcSharedClusters", 0.4, "TPC Max Shared Clusters"};
  Configurable<float> cMaxChi2Tpc{"cMaxChi2Tpc", 4, "Max TPC Chi2/ndf"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 3.0, "TPC nSigma PID cut"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.02, "Min proton DCA to PV (cm)"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.06, "Min pion DCA to PV (cm)"};
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Min DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1., "Max DCA between V0 daughters"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Min DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 999.0, "Max DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Min V0 decay radius (cm)"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Max V0 decay radius (cm)"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Min cTau (cm)"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Max cTau (cm)"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Min V0 cos(PA)"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.01, "K0s mass rejection window"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0s mass rejection flag"};

  // V0s kinmatic acceptance
  Configurable<float> cMinV0Mass{"cMinV0Mass", 1.10, "V0 Mass Min"};
  Configurable<float> cMaxV0Mass{"cMaxV0Mass", 1.12, "V0 Mass Max"};
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

    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(30, 0.5, 30.5, "");
    const AxisSpec axisCent(100, 0, 100, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(120, 1.08, 1.20, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(100., 0., 10., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(48, -1.2, 1.2, "y");
    const AxisSpec axisV0Eta(48, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

    const AxisSpec axisRadius(2000, 0, 200, "r (cm)");
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
    const AxisSpec axisNsigma(401, -10.025, 10.025, "n#sigma");
    const AxisSpec axisdEdx(360, 20, 200, "#frac{dE}{dx}");

    // Create Histograms.
    // Event histograms
    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCols});
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});

    histos.add("Tracks/h1f_tracks_info", "# of tracks", kTH1F, {axisTrks});
    histos.add("Tracks/h2f_armpod_before_sel", "Armenteros-Podolanski (before)", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h2f_armpod_after_sel", "Armenteros-Podolanski (after)", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h1f_lambda_pt_vs_invm", "p_{T} vs M_{#Lambda}", kTH2F, {axisV0Mass, axisV0Pt});
    histos.add("Tracks/h1f_antilambda_pt_vs_invm", "p_{T} vs M_{#bar{#Lambda}}", kTH2F, {axisV0Mass, axisV0Pt});

    histos.add("QA/Lambda/h2f_qt_vs_alpha", "Armenteros-Podolanski", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/Lambda/h1f_dca_V0_daughters", "DCA V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA/Lambda/h1f_dca_pos_to_PV", "DCA pos-prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_neg_to_PV", "DCA neg-prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA/Lambda/h1f_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA/Lambda/h1f_V0_radius", "V0 decay radius", kTH1F, {axisRadius});
    histos.add("QA/Lambda/h1f_V0_ctau", "c#tau", kTH1F, {axisCTau});
    histos.add("QA/Lambda/h1f_V0_gctau", "#gammac#tau", kTH1F, {axisGCTau});
    histos.add("QA/Lambda/h1f_pos_prong_pt", "Pos-prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_neg_prong_pt", "Neg-prong p_{T}", kTH1F, {axisTrackPt});
    histos.add("QA/Lambda/h1f_pos_prong_eta", "Pos-prong #eta", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_neg_prong_eta", "Neg-prong #eta", kTH1F, {axisV0Eta});
    histos.add("QA/Lambda/h1f_pos_prong_phi", "Pos-prong #phi", kTH1F, {axisV0Phi});
    histos.add("QA/Lambda/h1f_neg_prong_phi", "Neg-prong #phi", kTH1F, {axisV0Phi});
    histos.add("QA/Lambda/h2f_pos_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_neg_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_pos_prong_dEdx_vs_p", "TPC dE/dx pos", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_neg_prong_dEdx_vs_p", "TPC dE/dx neg", kTH2F, {axisMomPID, axisdEdx});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma_{p} pos", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma_{p} neg", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma_{#pi} pos", kTH2F, {axisMomPID, axisNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma_{#pi} neg", kTH2F, {axisMomPID, axisNsigma});

    histos.add("McRec/Lambda/hPt", "p_{T}", kTH1F, {axisV0Pt});
    histos.add("McRec/Lambda/hEta", "#eta", kTH1F, {axisV0Eta});
    histos.add("McRec/Lambda/hRap", "y", kTH1F, {axisV0Rap});
    histos.add("McRec/Lambda/hPhi", "#phi", kTH1F, {axisV0Phi});

    // QA Anti-Lambda
    histos.addClone("QA/Lambda/", "QA/AntiLambda/");
    histos.addClone("McRec/Lambda/", "McRec/AntiLambda/");

    if (doprocessMCRun3 || doprocessMCRun2) {
      histos.add("Tracks/h2f_tracks_pid_before_sel", "PIDs before sel", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_tracks_pid_after_sel", "PIDs after sel", kTH2F, {axisPID, axisV0Pt});
      histos.add("Tracks/h2f_lambda_mothers_pdg", "Lambda mothers", kTH2F, {axisPID, axisV0Pt});

      histos.add("McGen/h1f_collision_recgen", "RecGen collisions", kTH1F, {axisMult});
      histos.add("McGen/h1f_collisions_info", "Collisions info", kTH1F, {axisCols});
      histos.add("McGen/h2f_collision_posZ", "V_{z} rec vs gen", kTH2F, {axisVz, axisVz});
      histos.add("McGen/h2f_collision_cent", "Centrality rec vs gen", kTH2F, {axisCent, axisCent});
      histos.add("McGen/h1f_lambda_daughter_PDG", "Lambda dau PDG", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "AntiLambda dau PDG", kTH1F, {axisPID});

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
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx)
      return false;

    if constexpr (run == kRun3) {
      if (cCentEstimator == kCentFT0M)
        cent = col.centFT0M();
      else if (cCentEstimator == kCentFT0C)
        cent = col.centFT0C();
      if (cSel8Trig && !col.sel8())
        return false;
    } else {
      cent = col.centRun2V0M();
      if (cInt7Trig && !col.alias_bit(kINT7))
        return false;
      if (cSel7Trig && !col.sel7())
        return false;
    }

    if (cent <= cMinMult || cent >= cMaxMult)
      return false;
    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll))
      return false;

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
    if (!kinCutSelection(track.pt(), std::abs(track.eta()), cTrackMinPt, cTrackMaxPt, cTrackEtaCut))
      return false;
    if (track.tpcNClsCrossedRows() <= cMinTpcCrossedRows)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cMinTpcCROverCls)
      return false;
    if (track.tpcNClsShared() > cMaxTpcSharedClusters)
      return false;
    if (track.tpcChi2NCl() > cMaxChi2Tpc)
      return false;
    return true;
  }

  // Daughter Track Selection
  template <typename V, typename T>
  bool selDaughterTracks(V const& v0, T const&, ParticleType const& v0Type)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    if (!selTrack(posTrack) || !selTrack(negTrack))
      return false;

    // Apply DCA Selection on Daughter Tracks Based on Lambda/AntiLambda daughters
    float dcaProton = 0., dcaPion = 0.;
    if (v0Type == kLambda) {
      dcaProton = std::abs(v0.dcapostopv());
      dcaPion = std::abs(v0.dcanegtopv());
    } else if (v0Type == kAntiLambda) {
      dcaPion = std::abs(v0.dcapostopv());
      dcaProton = std::abs(v0.dcanegtopv());
    }
    if (dcaProton < cMinDcaProtonToPV || dcaPion < cMinDcaPionToPV)
      return false;
    return true;
  }

  template <typename C, typename V, typename T>
  bool topoCutSelection(C const& col, V const& v0, T const&)
  {
    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters)
      return false;
    if (v0.dcav0topv() <= cMinDcaV0ToPV || v0.dcav0topv() >= cMaxDcaV0ToPV)
      return false;
    if (v0.v0radius() <= cMinV0TransRadius || v0.v0radius() >= cMaxV0TransRadius)
      return false;

    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau <= cMinV0CTau || ctau >= cMaxV0CTau)
      return false;
    if (v0.v0cosPA() <= cMinV0CosPA)
      return false;
    return true;
  }

  template <ParticleType part, typename T>
  bool selLambdaDauWithTpcPid(T const& postrack, T const& negtrack)
  {
    float tpcNSigmaPr = 0., tpcNSigmaPi = 0.;
    if constexpr (part == kLambda) {
      tpcNSigmaPr = postrack.tpcNSigmaPr();
      tpcNSigmaPi = negtrack.tpcNSigmaPi();
    } else {
      tpcNSigmaPr = negtrack.tpcNSigmaPr();
      tpcNSigmaPi = postrack.tpcNSigmaPi();
    }
    return (std::abs(tpcNSigmaPr) < cTpcNsigmaCut && std::abs(tpcNSigmaPi) < cTpcNsigmaCut);
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

    if ((v0.mLambda() > cMinV0Mass && v0.mLambda() < cMaxV0Mass) &&
        selLambdaDauWithTpcPid<kLambda>(postrack, negtrack)) {
      lambdaFlag = true;
      v0type = kLambda;
    }
    if ((v0.mAntiLambda() > cMinV0Mass && v0.mAntiLambda() < cMaxV0Mass) &&
        selLambdaDauWithTpcPid<kAntiLambda>(postrack, negtrack)) {
      antiLambdaFlag = true;
      v0type = kAntiLambda;
    }

    if (!lambdaFlag && !antiLambdaFlag) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kNotLambdaNotAntiLambda);
      return false;
    }
    if (lambdaFlag && antiLambdaFlag) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsBothLambdaAntiLambda);
      return false;
    }
    return true;
  }

  template <typename C, typename V, typename T>
  bool selV0Particle(C const& col, V const& v0, T const& tracks, ParticleType& v0Type)
  {
    if (!selLambdaMassWindow(v0, tracks, v0Type))
      return false;
    histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsLambdaOrAntiLambda);

    if (!selDaughterTracks(v0, tracks, v0Type))
      return false;
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0DauTrackSel);

    float rap = cDoEtaAnalysis ? std::abs(v0.eta()) : std::abs(v0.yLambda());
    if (!kinCutSelection(v0.pt(), rap, cMinV0Pt, cMaxV0Pt, cMaxV0Rap))
      return false;
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0KinCuts);

    if (!topoCutSelection(col, v0, tracks))
      return false;
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0TopoSel);

    // All Selection Criterion Passed
    return true;
  }

  template <typename V, typename T>
  bool hasAmbiguousDaughters(V const& v0, T const&)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    auto posCC = posTrack.compatibleCollIds();
    auto negCC = negTrack.compatibleCollIds();
    if (posCC.size() > 1 || negCC.size() > 1)
      return true;
    if ((posCC.size() != 0 && posCC[0] != posTrack.collisionId()) ||
        (negCC.size() != 0 && negCC[0] != negTrack.collisionId()))
      return true;
    return false;
  }

  template <typename V>
  PrmScdType isPrimaryV0(V const& v0)
  {
    auto mcpart = v0.template mcParticle_as<aod::McParticles>();
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
    if (std::abs(mcpart.pdgCode()) != kLambda0)
      return false;

    if (cCheckRecoDauFlag) {
      auto postrack = v0.template posTrack_as<T>();
      auto negtrack = v0.template negTrack_as<T>();
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
    if (!cCorrectionFlag)
      return 1.;

    // Get  from CCDB
    auto ccdbObj = ccdb->getForTimeStamp<TList>(cPathCCDB.value, -1);
    if (!ccdbObj) {
      LOGF(warning, "CCDB OBJECT NOT FOUND");
      return 1.;
    }

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
        LOGF(warning, "CCDB: not a histogram!");
        effCorrFact = 1.;
      }
      delete histEff;
    }
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
        LOGF(warning, "CCDB: not a histogram!");
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
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();
    float mass = (part == kLambda) ? v0.mLambda() : v0.mAntiLambda();

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
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_pt"), negtrack.pt());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_eta"), postrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_neg_prong_eta"), negtrack.eta());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_pos_prong_phi"), postrack.phi());
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

  template <RunType run, DMCType dmc, typename C, typename B, typename V, typename T>
  void fillLambdaRecoTables(C const& collision, B const& bc, V const& v0tracks, T const& tracks)
  {
    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    if constexpr (dmc == kData) {
      if (!selCollision<run>(collision))
        return;
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());

    // Fill Collision Table
    lambdaCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ(), bc.timestamp());

    // initialize v0track objects
    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float mass = 0., corr_fact = 1.;
    float prPx = 0., prPy = 0., prPz = 0.;

    for (auto const& v0 : v0tracks) {
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kTracksBeforeHasMcParticle);
        if (!v0.has_mcParticle())
          continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      if (!selV0Particle(collision, v0, tracks, v0Type))
        continue;
      if (cV0TypeSelFlag && v0.v0Type() != cV0TypeSelection)
        continue;

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllSelPassed);

      if constexpr (run == kRun3) {
        if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks))
          continue;
      }

      mass = (v0Type == kLambda) ? v0.mLambda() : v0.mAntiLambda();
      pt = v0.pt();
      eta = v0.eta();
      rap = v0.yLambda();
      phi = v0.phi();

      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cSelMCPSV0)
          v0PrmScdType = isPrimaryV0(v0);
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks))
          continue;
        if (v0PrmScdType == kSecondary)
          fillLambdaMothers(v0, tracks);
        histos.fill(HIST("Tracks/h1f_tracks_info"), kPassTrueLambdaSel);
        histos.fill(HIST("Tracks/h2f_tracks_pid_after_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cRecoMomResoFlag) {
          auto mc = v0.template mcParticle_as<aod::McParticles>();
          pt = mc.pt();
          eta = mc.eta();
          rap = mc.y();
          phi = mc.phi();
          float y = cDoEtaAnalysis ? eta : rap;
          if (!kinCutSelection(pt, std::abs(y), cMinV0Pt, cMaxV0Pt, cMaxV0Rap))
            continue;
        }
      }

      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());
      corr_fact = (v0Type == kLambda) ? getCorrectionFactors<kLambda>(v0) : getCorrectionFactors<kAntiLambda>(v0);

      if (v0Type == kLambda) {
        prPx = v0.template posTrack_as<T>().px();
        prPy = v0.template posTrack_as<T>().py();
        prPz = v0.template posTrack_as<T>().pz();
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      } else {
        prPx = v0.template negTrack_as<T>().px();
        prPy = v0.template negTrack_as<T>().py();
        prPz = v0.template negTrack_as<T>().pz();
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kAntiLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      }

      lambdaTrackTable(lambdaCollisionTable.lastIndex(),
                       v0.px(), v0.py(), v0.pz(), pt, eta, phi, rap, mass,
                       prPx, prPy, prPz,
                       v0.template posTrack_as<T>().index(), v0.template negTrack_as<T>().index(),
                       v0.v0cosPA(), v0.dcaV0daughters(), (int8_t)v0Type, v0PrmScdType, corr_fact);
    }
  }

  template <RunType run, typename C, typename M>
  void fillLambdaMcGenTables(C const& mcCollision, M const& mcParticles)
  {
    lambdaMCGenCollisionTable(cent, mult, mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float rap = 0.;
    float prPx = 0., prPy = 0., prPz = 0.;

    for (auto const& mcpart : mcParticles) {
      if (mcpart.pdgCode() == kLambda0)
        v0Type = kLambda;
      else if (mcpart.pdgCode() == kLambda0Bar)
        v0Type = kAntiLambda;
      else
        continue;

      v0PrmScdType = mcpart.isPhysicalPrimary() ? kPrimary : kSecondary;
      rap = cDoEtaAnalysis ? mcpart.eta() : mcpart.y();
      if (!kinCutSelection(mcpart.pt(), std::abs(rap), cMinV0Pt, cMaxV0Pt, cMaxV0Rap))
        continue;
      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenTotAccLambda);

      if (!mcpart.has_daughters()) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaNoDau);
        continue;
      }
      auto dautracks = mcpart.template daughters_as<aod::McParticles>();
      std::vector<int> daughterPDGs, daughterIDs;
      std::vector<float> vDauPt, vDauEta, vDauRap, vDauPhi;
      std::vector<float> vDauPx, vDauPy, vDauPz;
      for (auto const& dautrack : dautracks) {
        daughterPDGs.push_back(dautrack.pdgCode());
        daughterIDs.push_back(dautrack.globalIndex());
        vDauPt.push_back(dautrack.pt());
        vDauEta.push_back(dautrack.eta());
        vDauRap.push_back(dautrack.y());
        vDauPhi.push_back(dautrack.phi());
        vDauPx.push_back(dautrack.px());
        vDauPy.push_back(dautrack.py());
        vDauPz.push_back(dautrack.pz());
      }
      if (cGenDecayChannel) {
        if (v0Type == kLambda && (daughterPDGs[0] != kProton || daughterPDGs[1] != kPiMinus))
          continue;
        if (v0Type == kAntiLambda && (daughterPDGs[0] != kProtonBar || daughterPDGs[1] != kPiPlus))
          continue;
      }
      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaToPrPi);

      if (v0Type == kLambda) {
        prPx = vDauPx[0];
        prPy = vDauPy[0];
        prPz = vDauPz[0];
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
        prPx = vDauPx[1];
        prPy = vDauPy[1];
        prPz = vDauPz[1];
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

      lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(),
                            mcpart.px(), mcpart.py(), mcpart.pz(),
                            mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(),
                            RecoDecay::m(mcpart.p(), mcpart.e()),
                            prPx, prPy, prPz,
                            daughterIDs[0], daughterIDs[1],
                            (int8_t)v0Type, -999., -999., v0PrmScdType, 1.);
    }
  }

  template <RunType run, DMCType dmc, typename M, typename C, typename V, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, V const& /*V0s*/, T const& /*tracks*/, P const& mcParticles)
  {
    int nRecCols = collisions.size();
    if (nRecCols != 0)
      histos.fill(HIST("McGen/h1f_collision_recgen"), nRecCols);
    if (nRecCols != 1)
      return;

    histos.fill(HIST("McGen/h1f_collisions_info"), kTotCol);
    if (!collisions.begin().has_mcCollision() ||
        !selCollision<run>(collisions.begin()) ||
        collisions.begin().mcCollisionId() != mcCollision.globalIndex())
      return;
    histos.fill(HIST("McGen/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(), collisions.begin().posZ());
    fillLambdaMcGenTables<run>(mcCollision, mcParticles);
  }

  SliceCache cache;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCollision = aod::v0data::collisionId;

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>;
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::PVMults>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::TrackCompColls>;
  using TracksRun2 = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using TracksMCRun2 = soa::Join<TracksRun2, aod::McTrackLabels>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processDataRun3(CollisionsRun3::iterator const& collision, aod::BCsWithTimestamps const&,
                       aod::V0Datas const& V0s, Tracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    fillLambdaRecoTables<kRun3, kData>(collision, bc, V0s, tracks);
  }
  PROCESS_SWITCH(LambdaTableProducer, processDataRun3, "Process for Run3 DATA", true);

  void processMCRun3(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun3, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMC const& tracks, aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun3, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }
  PROCESS_SWITCH(LambdaTableProducer, processMCRun3, "Process for Run3 MC RecoGen", false);

  void processMCRun2(aod::McCollisions::iterator const& mcCollision,
                     soa::SmallGroups<soa::Join<CollisionsRun2, aod::McCollisionLabels>> const& collisions,
                     McV0Tracks const& V0s, TracksMCRun2 const& tracks, aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun2, kMC>(mcCollision, collisions, V0s, tracks, mcParticles);
  }
  PROCESS_SWITCH(LambdaTableProducer, processMCRun2, "Process for Run2 MC RecoGen", false);
};

struct LambdaTracksExtProducer {

  Produces<aod::LambdaTracksExt> lambdaTrackExtTable;

  Configurable<bool> cAcceptAllLambda{"cAcceptAllLambda", false, "Accept all lambda (ignore sharing)"};
  Configurable<bool> cRejAllLambdaShaDau{"cRejAllLambdaShaDau", true, "Reject lambda sharing daughters"};
  Configurable<bool> cSelLambdaMassPdg{"cSelLambdaMassPdg", false, "Select lambda closest to PDG mass"};
  Configurable<bool> cSelLambdaTScore{"cSelLambdaTScore", false, "Select lambda by t-score"};
  Configurable<float> cA{"cA", 0.6, "t-score weight: |mass - PDGmass|"};
  Configurable<float> cB{"cB", 0.6, "t-score weight: DCA daughters"};
  Configurable<float> cC{"cC", 0.6, "t-score weight: |cosPA - 1|"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisMult(10, 0, 10);
    const AxisSpec axisMass(120, 1.08, 1.20, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisCPA(100, 0.995, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daughter DCA (#sigma)");
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_lambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_antilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h2d_n2_etaphi_LaP_LaM", "#rho_{2}^{Share} #Lambda#bar{#Lambda}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaM_LaP", "#rho_{2}^{Share} #bar{#Lambda}#Lambda", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaP_LaP", "#rho_{2}^{Share} #Lambda#Lambda", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaM_LaM", "#rho_{2}^{Share} #bar{#Lambda}#bar{#Lambda}", kTH2D, {axisDEta, axisDPhi});

    histos.add("Reco/h1f_lambda_invmass", "M_{#Lambda}", kTH1F, {axisMass});
    histos.add("Reco/h1f_lambda_cospa", "cos(PA)", kTH1F, {axisCPA});
    histos.add("Reco/h1f_lambda_dcadau", "DCA daughters", kTH1F, {axisDcaDau});
    histos.add("Reco/h1f_antilambda_invmass", "M_{#bar{#Lambda}}", kTH1F, {axisMass});
    histos.add("Reco/h1f_antilambda_cospa", "cos(PA)", kTH1F, {axisCPA});
    histos.add("Reco/h1f_antilambda_dcadau", "DCA daughters", kTH1F, {axisDcaDau});
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

      if (lambda.v0Type() == kLambda)
        ++nTotLambda;
      else if (lambda.v0Type() == kAntiLambda)
        ++nTotAntiLambda;

      tLambda = (cA * std::abs(lambda.mass() - MassLambda0)) +
                (cB * lambda.dcaDau()) +
                (cC * std::abs(lambda.cosPA() - 1.));

      for (auto const& track : tracks) {
        if (lambda.index() == track.index())
          continue;

        if (lambda.posTrackId() == track.posTrackId() || lambda.negTrackId() == track.negTrackId()) {
          vSharedDauLambdaIndex.push_back(track.index());
          lambdaSharingDauFlag = true;

          if (lambda.v0Type() == kLambda && track.v0Type() == kAntiLambda)
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaM"), lambda.eta() - track.eta(),
                        RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kLambda)
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaP"), lambda.eta() - track.eta(),
                        RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          else if (lambda.v0Type() == kLambda && track.v0Type() == kLambda)
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaP"), lambda.eta() - track.eta(),
                        RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kAntiLambda)
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaM"), lambda.eta() - track.eta(),
                        RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));

          if (std::abs(lambda.mass() - MassLambda0) > std::abs(track.mass() - MassLambda0))
            lambdaMinDeltaMassFlag = false;

          tTrack = (cA * std::abs(track.mass() - MassLambda0)) +
                   (cB * track.dcaDau()) +
                   (cC * std::abs(track.cosPA() - 1.));
          if (tLambda > tTrack)
            lambdaMinTScoreFlag = false;
        }
      }

      if (lambdaSharingDauFlag)
        fillHistos<kLambdaShareDau>(lambda);
      else
        fillHistos<kUniqueLambda>(lambda);

      if (cAcceptAllLambda)
        trueLambdaFlag = true;
      else if (cRejAllLambdaShaDau && !lambdaSharingDauFlag)
        trueLambdaFlag = true;
      else if (cSelLambdaMassPdg && lambdaMinDeltaMassFlag)
        trueLambdaFlag = true;
      else if (cSelLambdaTScore && lambdaMinTScoreFlag)
        trueLambdaFlag = true;

      if (trueLambdaFlag) {
        if (lambda.v0Type() == kLambda)
          ++nSelLambda;
        else if (lambda.v0Type() == kAntiLambda)
          ++nSelAntiLambda;
      }

      lambdaTrackExtTable(lambdaSharingDauFlag, vSharedDauLambdaIndex, trueLambdaFlag);
    }

    if (nTotLambda != 0)
      histos.fill(HIST("h1i_totlambda_mult"), nTotLambda);
    if (nTotAntiLambda != 0)
      histos.fill(HIST("h1i_totantilambda_mult"), nTotAntiLambda);
    if (nSelLambda != 0)
      histos.fill(HIST("h1i_lambda_mult"), nSelLambda);
    if (nSelAntiLambda != 0)
      histos.fill(HIST("h1i_antilambda_mult"), nSelAntiLambda);
  }
};

struct LambdaSpinPolarization {

  Produces<aod::LambdaMixEventCollisions> lambdaMixEvtCol;
  Produces<aod::LambdaMixEventTracks> lambdaMixEvtTrk;

  Configurable<int> cNPtBins{"cNPtBins", 30, "N pT bins"};
  Configurable<float> cMinPt{"cMinPt", 0.5f, "pT min (GeV/c)"};
  Configurable<float> cMaxPt{"cMaxPt", 4.5f, "pT max (GeV/c)"};
  Configurable<int> cNRapBins{"cNRapBins", 10, "N rapidity bins"};
  Configurable<float> cMinRap{"cMinRap", -0.5f, "Rapidity min"};
  Configurable<float> cMaxRap{"cMaxRap", 0.5f, "Rapidity max"};
  Configurable<int> cNPhiBins{"cNPhiBins", 36, "N phi bins"};
  Configurable<int> cNBinsCosTS{"cNBinsCosTS", 20, "N costheta* bins"};
  Configurable<int> cNBinsDeltaR{"cNBinsDeltaR", 20, "N DeltaR bins"};

  Configurable<float> cMassHistMin{"cMassHistMin", 1.08f, "Mass histogram min (GeV/c2)"};
  Configurable<float> cMassHistMax{"cMassHistMax", 1.20f, "Mass histogram max (GeV/c2)"};
  Configurable<int> cNMassBins{"cNMassBins", 120, "Mass histogram N bins"};

  Configurable<float> cSigMinLambda{"cSigMinLambda", 1.108f, "Signal region min (GeV/c2)"};
  Configurable<float> cSigMaxLambda{"cSigMaxLambda", 1.123f, "Signal region max (GeV/c2)"};
  Configurable<float> cSbLeftMin{"cSbLeftMin", 1.080f, "Left sideband min (GeV/c2)"};
  Configurable<float> cSbLeftMax{"cSbLeftMax", 1.100f, "Left sideband max (GeV/c2)"};
  Configurable<float> cSbRightMin{"cSbRightMin", 1.135f, "Right sideband min (GeV/c2)"};
  Configurable<float> cSbRightMax{"cSbRightMax", 1.155f, "Right sideband max (GeV/c2)"};

  Configurable<bool> cInvBoostFlag{"cInvBoostFlag", true, "Inverse boost flag"};
  Configurable<bool> cDoAtlasMethod{"cDoAtlasMethod", false, "Fill pair-boost (ATLAS) histograms"};
  Configurable<bool> cDoStarMethod{"cDoStarMethod", true, "Fill lab-boost (STAR) histograms"};
  Configurable<int> mixingParameter{"mixingParameter", 5, "ME pool depth"};
  Configurable<int> cMEMode{"cMEMode", 1, "ME mode: 0=standard, 1=kinematicConstrained"};

  ConfigurableAxis cMultBins{"cMultBins", {VARIABLE_WIDTH, 0.f, 10.f, 30.f, 50.f, 80.f, 100.f}, "Multiplicity bins"};
  ConfigurableAxis axisCentME{"axisCentME", {VARIABLE_WIDTH, 0, 10, 30, 50, 100}, "ME centrality bins"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -7, -3, 0, 3, 7}, "ME vtxZ bins"};

  Configurable<float> cMaxDeltaPt{"cMaxDeltaPt", 0.1f, "Kinematic ME: max |deltaPt(SE)-deltaPt(ME)| (GeV/c)"};
  Configurable<float> cMaxDeltaPhi{"cMaxDeltaPhi", 0.1f, "Kinematic ME: max |deltaPhi(SE)-deltaPhi(ME)| (rad)"};
  Configurable<float> cMaxDeltaRap{"cMaxDeltaRap", 0.1f, "Kinematic ME: max |deltaRap(SE)-deltaRap(ME)|"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  float cent = 0.;
  struct PoolTrack {
    float _px, _py, _pz, _pt, _rap, _phi, _mass;
    float _prPx, _prPy, _prPz;
    float px() const { return _px; }
    float py() const { return _py; }
    float pz() const { return _pz; }
    float pt() const { return _pt; }
    float rap() const { return _rap; }
    float phi() const { return _phi; }
    float mass() const { return _mass; }
    float prPx() const { return _prPx; }
    float prPy() const { return _prPy; }
    float prPz() const { return _prPz; }
  };

  template <typename T>
  PoolTrack toPoolTrack(T const& trk)
  {
    return PoolTrack{trk.px(), trk.py(), trk.pz(), trk.pt(),
                     trk.rap(), trk.phi(), trk.mass(),
                     trk.prPx(), trk.prPy(), trk.prPz()};
  }

  void init(InitContext const&)
  {
    const AxisSpec axisCheck(1, 0, 1, "");
    const AxisSpec axisCent(cMultBins, "FT0M (%)");

    const AxisSpec axisMass(cNMassBins, cMassHistMin, cMassHistMax, "M_{#Lambda} (GeV/#it{c}^{2})");
    const AxisSpec axisPt(cNPtBins, cMinPt, cMaxPt, "p_{T} (GeV/#it{c})");
    const AxisSpec axisDRap(2 * cNRapBins, cMinRap - cMaxRap, cMaxRap - cMinRap, "#Deltay");
    const AxisSpec axisDPhi(cNPhiBins, -PI, PI, "#Delta#varphi");
    const AxisSpec axisCosTS(cNBinsCosTS, -1, 1, "cos(#theta*)");
    const AxisSpec axisDR(cNBinsDeltaR, 0, 4, "#DeltaR");
    const AxisSpec axisPosZ(220, -7, 7, "V_{z} (cm)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");

    histos.add("QA/ME/hPoolCentVz", "ME pool;cent (%);V_{z}", kTH2F, {axisCentME, axisVtxZME});
    histos.add("QA/ME/hLambdaMultVsCent", "ME #Lambda mult;cent;N", kTH2F, {axisCentME, {50, 0, 50}});
    histos.add("QA/ME/hAntiLambdaMultVsCent", "ME #bar{#Lambda} mult;cent;N", kTH2F, {axisCentME, {50, 0, 50}});

    histos.add("SE/Reco/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/Reco/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/Reco/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/Reco/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("SE/RecoBkgSigSB/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSigSB/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSigSB/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSigSB/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("SE/RecoBkgSBSB/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSBSB/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSBSB/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("SE/RecoBkgSBSB/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("ME/Reco/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/Reco/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/Reco/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/Reco/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} inclusive",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("ME/RecoBkgSigSB/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSigSB/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSigSB/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSigSB/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} Sig#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("ME/RecoBkgSBSB/h2f_n2_mass_LaPLaM", "M_{inv}: #Lambda#bar{#Lambda} SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSBSB/h2f_n2_mass_LaMLaP", "M_{inv}: #bar{#Lambda}#Lambda SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSBSB/h2f_n2_mass_LaPLaP", "M_{inv}: #Lambda#Lambda SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});
    histos.add("ME/RecoBkgSBSB/h2f_n2_mass_LaMLaM", "M_{inv}: #bar{#Lambda}#bar{#Lambda} SB#timesSB",
               kTHnSparseF, {axisMass, axisMass, axisPt, axisPt, axisDRap, axisDPhi, axisDR});

    histos.add("SE/RecoCorr/Star/h2f_n2_dltaR_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_dltaR_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_dltaR_LaPLaP", "#rho_{2} #Lambda#Lambda [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_dltaR_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});

    histos.add("SE/RecoCorr/Star/h2f_n2_ctheta_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_ctheta_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_ctheta_LaPLaP", "#rho_{2} #Lambda#Lambda [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Star/h2f_n2_ctheta_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});

    histos.add("SE/RecoCorr/Atlas/h2f_n2_dltaR_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_dltaR_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_dltaR_LaPLaP", "#rho_{2} #Lambda#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_dltaR_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});

    histos.add("SE/RecoCorr/Atlas/h2f_n2_ctheta_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_ctheta_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_ctheta_LaPLaP", "#rho_{2} #Lambda#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("SE/RecoCorr/Atlas/h2f_n2_ctheta_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});

    histos.add("ME/RecoCorr/Star/h2f_n2_dltaR_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_dltaR_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_dltaR_LaPLaP", "#rho_{2} #Lambda#Lambda [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_dltaR_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDR, axisCosTS});

    histos.add("ME/RecoCorr/Star/h2f_n2_ctheta_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_ctheta_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_ctheta_LaPLaP", "#rho_{2} #Lambda#Lambda [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Star/h2f_n2_ctheta_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Star]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});

    histos.add("ME/RecoCorr/Atlas/h2f_n2_dltaR_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_dltaR_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_dltaR_LaPLaP", "#rho_{2} #Lambda#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_dltaR_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDR, axisCosTS});

    histos.add("ME/RecoCorr/Atlas/h2f_n2_ctheta_LaPLaM", "#rho_{2} #Lambda#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_ctheta_LaMLaP", "#rho_{2} #bar{#Lambda}#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_ctheta_LaPLaP", "#rho_{2} #Lambda#Lambda [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});
    histos.add("ME/RecoCorr/Atlas/h2f_n2_ctheta_LaMLaM", "#rho_{2} #bar{#Lambda}#bar{#Lambda} [Atlas]", kTHnSparseF, {axisCent, axisDRap, axisDPhi, axisCosTS});

    histos.addClone("SE/RecoCorr/Star/", "SE/RecoCorrBkgSigSB/Star/");
    histos.addClone("SE/RecoCorr/Atlas/", "SE/RecoCorrBkgSigSB/Atlas/");
    histos.addClone("SE/RecoCorr/Star/", "SE/RecoCorrBkgSBSB/Star/");
    histos.addClone("SE/RecoCorr/Atlas/", "SE/RecoCorrBkgSBSB/Atlas/");

    histos.addClone("ME/RecoCorr/Star/", "ME/RecoCorrBkgSigSB/Star/");
    histos.addClone("ME/RecoCorr/Atlas/", "ME/RecoCorrBkgSigSB/Atlas/");
    histos.addClone("ME/RecoCorr/Star/", "ME/RecoCorrBkgSBSB/Star/");
    histos.addClone("ME/RecoCorr/Atlas/", "ME/RecoCorrBkgSBSB/Atlas/");
  }

  bool isSignal(float m) const
  {
    return (m >= cSigMinLambda.value && m <= cSigMaxLambda.value);
  }
  bool isSideband(float m) const
  {
    return ((m >= cSbLeftMin.value && m <= cSbLeftMax.value) ||
            (m >= cSbRightMin.value && m <= cSbRightMax.value));
  }
  bool isInclusive(float m) const
  {
    return (m >= cMassHistMin.value && m <= cMassHistMax.value);
  }

  void getBoostVector(std::array<float, 4> const& p, std::array<float, 3>& v, bool inverseBoostFlag = true)
  {
    int n = p.size();
    for (int i = 0; i < n - 1; ++i) {
      if (inverseBoostFlag) {
        v[i] = -p[i] / RecoDecay::e(p[0], p[1], p[2], p[3]);
      } else {
        v[i] = p[i] / RecoDecay::e(p[0], p[1], p[2], p[3]);
      }
    }
  }

  void boost(std::array<float, 4>& p, std::array<float, 3> const& b)
  {
    float e = RecoDecay::e(p[0], p[1], p[2], p[3]);
    float b2 = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
    float gam = 1.f / std::sqrt(1.f - b2);
    float bp = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
    float gam2 = (b2 > 0.f) ? (gam - 1.f) / b2 : 0.f;
    p[0] = p[0] + gam2 * bp * b[0] + gam * b[0] * e;
    p[1] = p[1] + gam2 * bp * b[1] + gam * b[1] * e;
    p[2] = p[2] + gam2 * bp * b[2] + gam * b[2] * e;
  }

  template <ParticlePairType part_pair, typename U>
  void fillPairHistos(U const& p1, U const& p2)
  {
    static constexpr std::string_view SubDir[] = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};

    constexpr bool IsSigSB = (part_pair == kLambdaSBAntiLambda ||
                              part_pair == kAntiLambdaSBLambda ||
                              part_pair == kLambdaSBLambda ||
                              part_pair == kAntiLambdaSBAntiLambda);
    constexpr bool IsSBSB = (part_pair == kLambdaSBSBAntiLambda ||
                             part_pair == kAntiLambdaSBSBLambda ||
                             part_pair == kLambdaSBSBLambda ||
                             part_pair == kAntiLambdaSBSBAntiLambda);

    constexpr int Idx = IsSigSB  ? static_cast<int>(part_pair) - 4
                        : IsSBSB ? static_cast<int>(part_pair) - 8
                                 : static_cast<int>(part_pair);

    float drap = p1.rap() - p2.rap();
    float dphi = RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PI);
    float dR = std::sqrt(drap * drap + dphi * dphi);

    if constexpr (!IsSigSB && !IsSBSB) {
      histos.fill(HIST("SE/Reco/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR);
    } else if constexpr (IsSigSB) {
      histos.fill(HIST("SE/RecoBkgSigSB/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR);
    } else {
      histos.fill(HIST("SE/RecoBkgSBSB/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR);
    }

    std::array<float, 4> l1 = {p1.px(), p1.py(), p1.pz(), MassLambda0};
    std::array<float, 4> l2 = {p2.px(), p2.py(), p2.pz(), MassLambda0};
    std::array<float, 4> pr1 = {p1.prPx(), p1.prPy(), p1.prPz(), MassProton};
    std::array<float, 4> pr2 = {p2.prPx(), p2.prPy(), p2.prPz(), MassProton};

    if (cDoAtlasMethod) {
      auto l1a = l1;
      auto l2a = l2;
      auto pr1a = pr1;
      auto pr2a = pr2;
      std::array<float, 4> llpair = {l1a[0] + l2a[0], l1a[1] + l2a[1], l1a[2] + l2a[2], l1a[3] + l2a[3]};
      std::array<float, 3> vPair;
      getBoostVector(llpair, vPair, cInvBoostFlag);
      boost(l1a, vPair);
      boost(l2a, vPair);
      boost(pr1a, vPair);
      boost(pr2a, vPair);
      std::array<float, 3> v1p, v2p;
      getBoostVector(l1a, v1p, cInvBoostFlag);
      getBoostVector(l2a, v2p, cInvBoostFlag);
      boost(pr1a, v1p);
      boost(pr2a, v2p);
      std::array<float, 3> n1 = {pr1a[0], pr1a[1], pr1a[2]}, n2 = {pr2a[0], pr2a[1], pr2a[2]};
      float ctheta = RecoDecay::dotProd(n1, n2) /
                     (RecoDecay::sqrtSumOfSquares(n1[0], n1[1], n1[2]) *
                      RecoDecay::sqrtSumOfSquares(n2[0], n2[1], n2[2]));
      if constexpr (!IsSigSB && !IsSBSB) {
        histos.fill(HIST("SE/RecoCorr/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorr/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      } else if constexpr (IsSigSB) {
        histos.fill(HIST("SE/RecoCorrBkgSigSB/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorrBkgSigSB/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      } else {
        histos.fill(HIST("SE/RecoCorrBkgSBSB/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorrBkgSBSB/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      }
    }

    if (cDoStarMethod) {
      auto pr1s = pr1;
      auto pr2s = pr2;
      std::array<float, 3> v1lab, v2lab;
      getBoostVector(l1, v1lab, cInvBoostFlag);
      getBoostVector(l2, v2lab, cInvBoostFlag);
      boost(pr1s, v1lab);
      boost(pr2s, v2lab);
      std::array<float, 3> n1 = {pr1s[0], pr1s[1], pr1s[2]}, n2 = {pr2s[0], pr2s[1], pr2s[2]};
      float ctheta = RecoDecay::dotProd(n1, n2) /
                     (RecoDecay::sqrtSumOfSquares(n1[0], n1[1], n1[2]) *
                      RecoDecay::sqrtSumOfSquares(n2[0], n2[1], n2[2]));
      if constexpr (!IsSigSB && !IsSBSB) {
        histos.fill(HIST("SE/RecoCorr/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorr/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      } else if constexpr (IsSigSB) {
        histos.fill(HIST("SE/RecoCorrBkgSigSB/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorrBkgSigSB/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      } else {
        histos.fill(HIST("SE/RecoCorrBkgSBSB/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]), cent, drap, dphi, ctheta);
        histos.fill(HIST("SE/RecoCorrBkgSBSB/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]), cent, dR, ctheta);
      }
    }
  }

  template <ParticlePairType part_pair>
  void fillPairHistosWeighted(PoolTrack const& p1, PoolTrack const& p2, float w)
  {
    static constexpr std::string_view SubDir[] = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};

    constexpr bool IsSigSB = (part_pair == kLambdaSBAntiLambda ||
                              part_pair == kAntiLambdaSBLambda ||
                              part_pair == kLambdaSBLambda ||
                              part_pair == kAntiLambdaSBAntiLambda);
    constexpr bool IsSBSB = (part_pair == kLambdaSBSBAntiLambda ||
                             part_pair == kAntiLambdaSBSBLambda ||
                             part_pair == kLambdaSBSBLambda ||
                             part_pair == kAntiLambdaSBSBAntiLambda);

    constexpr int Idx = IsSigSB  ? static_cast<int>(part_pair) - 4
                        : IsSBSB ? static_cast<int>(part_pair) - 8
                                 : static_cast<int>(part_pair);

    float drap = p1.rap() - p2.rap();
    float dphi = RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PI);
    float dR = std::sqrt(drap * drap + dphi * dphi);

    if constexpr (!IsSigSB && !IsSBSB) {
      histos.fill(HIST("ME/Reco/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR, w);
    } else if constexpr (IsSigSB) {
      histos.fill(HIST("ME/RecoBkgSigSB/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR, w);
    } else {
      histos.fill(HIST("ME/RecoBkgSBSB/h2f_n2_mass_") + HIST(SubDir[Idx]),
                  p1.mass(), p2.mass(), p1.pt(), p2.pt(), drap, dphi, dR, w);
    }

    std::array<float, 4> l1 = {p1.px(), p1.py(), p1.pz(), MassLambda0};
    std::array<float, 4> l2 = {p2.px(), p2.py(), p2.pz(), MassLambda0};
    std::array<float, 4> pr1 = {p1.prPx(), p1.prPy(), p1.prPz(), MassProton};
    std::array<float, 4> pr2 = {p2.prPx(), p2.prPy(), p2.prPz(), MassProton};

    if (cDoAtlasMethod) {
      auto l1a = l1;
      auto l2a = l2;
      auto pr1a = pr1;
      auto pr2a = pr2;
      std::array<float, 4> llpair = {l1a[0] + l2a[0], l1a[1] + l2a[1],
                                     l1a[2] + l2a[2], l1a[3] + l2a[3]};
      std::array<float, 3> vPair;
      getBoostVector(llpair, vPair, cInvBoostFlag);
      boost(l1a, vPair);
      boost(l2a, vPair);
      boost(pr1a, vPair);
      boost(pr2a, vPair);
      std::array<float, 3> v1p, v2p;
      getBoostVector(l1a, v1p, cInvBoostFlag);
      getBoostVector(l2a, v2p, cInvBoostFlag);
      boost(pr1a, v1p);
      boost(pr2a, v2p);
      std::array<float, 3> n1 = {pr1a[0], pr1a[1], pr1a[2]},
                           n2 = {pr2a[0], pr2a[1], pr2a[2]};
      float ctheta = RecoDecay::dotProd(n1, n2) /
                     (RecoDecay::sqrtSumOfSquares(n1[0], n1[1], n1[2]) *
                      RecoDecay::sqrtSumOfSquares(n2[0], n2[1], n2[2]));

      if constexpr (!IsSigSB && !IsSBSB) {
        histos.fill(HIST("ME/RecoCorr/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorr/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      } else if constexpr (IsSigSB) {
        histos.fill(HIST("ME/RecoCorrBkgSigSB/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorrBkgSigSB/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      } else {
        histos.fill(HIST("ME/RecoCorrBkgSBSB/Atlas/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorrBkgSBSB/Atlas/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      }
    }

    if (cDoStarMethod) {
      auto pr1s = pr1;
      auto pr2s = pr2;
      std::array<float, 3> v1lab, v2lab;
      getBoostVector(l1, v1lab, cInvBoostFlag);
      getBoostVector(l2, v2lab, cInvBoostFlag);
      boost(pr1s, v1lab);
      boost(pr2s, v2lab);
      std::array<float, 3> n1 = {pr1s[0], pr1s[1], pr1s[2]},
                           n2 = {pr2s[0], pr2s[1], pr2s[2]};
      float ctheta = RecoDecay::dotProd(n1, n2) /
                     (RecoDecay::sqrtSumOfSquares(n1[0], n1[1], n1[2]) *
                      RecoDecay::sqrtSumOfSquares(n2[0], n2[1], n2[2]));

      if constexpr (!IsSigSB && !IsSBSB) {
        histos.fill(HIST("ME/RecoCorr/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorr/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      } else if constexpr (IsSigSB) {
        histos.fill(HIST("ME/RecoCorrBkgSigSB/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorrBkgSigSB/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      } else {
        histos.fill(HIST("ME/RecoCorrBkgSBSB/Star/h2f_n2_ctheta_") + HIST(SubDir[Idx]),
                    cent, drap, dphi, ctheta, w);
        histos.fill(HIST("ME/RecoCorrBkgSBSB/Star/h2f_n2_dltaR_") + HIST(SubDir[Idx]),
                    cent, dR, ctheta, w);
      }
    }
  }

  template <ParticlePairType part_pair_sig,
            ParticlePairType part_pair_bkg_sigsb,
            ParticlePairType part_pair_bkg_sbsb,
            bool samelambda, typename T>
  void analyzePairsWithMassWindow(T const& trks_1, T const& trks_2)
  {
    for (auto const& trk_1 : trks_1) {
      if (!isInclusive(trk_1.mass()))
        continue;
      bool t1sig = isSignal(trk_1.mass());
      bool t1sb = isSideband(trk_1.mass());

      for (auto const& trk_2 : trks_2) {
        if constexpr (samelambda) {
          if (trk_1.index() == trk_2.index())
            continue;
        }
        if (!isInclusive(trk_2.mass()))
          continue;
        bool t2sig = isSignal(trk_2.mass());
        bool t2sb = isSideband(trk_2.mass());

        fillPairHistos<part_pair_sig>(trk_1, trk_2);
        if ((t1sig && t2sb) || (t1sb && t2sig))
          fillPairHistos<part_pair_bkg_sigsb>(trk_1, trk_2);
        if (t1sb && t2sb)
          fillPairHistos<part_pair_bkg_sbsb>(trk_1, trk_2);
      }
    }
  }

  static constexpr int MEModeStandard = 0;
  static constexpr int MEModeKinematic = 1;

  template <ParticlePairType part_pair_sig,
            ParticlePairType part_pair_bkg_sigsb,
            ParticlePairType part_pair_bkg_sbsb,
            bool samelambda, typename T>
  void analyzePairsME(T const& trks_1, T const& trks_2)
  {

    for (auto const& trk_1 : trks_1) {
      if (!isInclusive(trk_1.mass()))
        continue;
      const bool t1sig = isSignal(trk_1.mass());
      const bool t1sb = isSideband(trk_1.mass());
      PoolTrack p1 = toPoolTrack(trk_1);

      for (auto const& trk_2 : trks_2) {
        if constexpr (samelambda) {
          if (trk_1.index() == trk_2.index())
            continue;
        }
        if (!isInclusive(trk_2.mass()))
          continue;
        const bool t2sig = isSignal(trk_2.mass());
        const bool t2sb = isSideband(trk_2.mass());
        PoolTrack p2 = toPoolTrack(trk_2);

        fillPairHistosWeighted<part_pair_sig>(p1, p2, 1.0f);

        if ((t1sig && t2sb) || (t1sb && t2sig))
          fillPairHistosWeighted<part_pair_bkg_sigsb>(p1, p2, 1.0f);

        if (t1sb && t2sb)
          fillPairHistosWeighted<part_pair_bkg_sbsb>(p1, p2, 1.0f);
      }
    }
  }

  template <ParticlePairType part_pair_sig,
            ParticlePairType part_pair_bkg_sigsb,
            ParticlePairType part_pair_bkg_sbsb,
            bool samelambda, typename T>
  void analyzePairsMEKinematic(T const& se_trks_1, T const& se_trks_2,
                               T const& me_pool_1, T const& me_pool_2)
  {
    std::vector<PoolTrack> meVec1, meVec2;
    meVec1.reserve(me_pool_1.size());
    meVec2.reserve(me_pool_2.size());
    for (auto const& me : me_pool_1)
      if (isInclusive(me.mass()))
        meVec1.push_back(toPoolTrack(me));
    for (auto const& me : me_pool_2)
      if (isInclusive(me.mass()))
        meVec2.push_back(toPoolTrack(me));

    if (meVec1.empty() && meVec2.empty())
      return;

    for (auto const& trk1 : se_trks_1) {
      if (!isInclusive(trk1.mass()))
        continue;
      const bool se1sig = isSignal(trk1.mass());
      const bool se1sb = isSideband(trk1.mass());
      PoolTrack p1 = toPoolTrack(trk1);

      for (auto const& trk2 : se_trks_2) {
        if constexpr (samelambda) {
          if (trk1.index() == trk2.index())
            continue;
        }
        if (!isInclusive(trk2.mass()))
          continue;
        const bool se2sig = isSignal(trk2.mass());
        const bool se2sb = isSideband(trk2.mass());
        PoolTrack p2 = toPoolTrack(trk2);

        {
          std::vector<PoolTrack> matchA;
          for (auto const& meP : meVec2) {
            if (std::abs(meP.pt() - p2.pt()) < cMaxDeltaPt &&
                std::abs(meP.rap() - p2.rap()) < cMaxDeltaRap &&
                std::abs(RecoDecay::constrainAngle(meP.phi() - p2.phi(), -PI)) < cMaxDeltaPhi)
              matchA.push_back(meP);
          }
          if (!matchA.empty()) {
            const float wA = 1.0f / static_cast<float>(matchA.size());
            for (auto const& meP2 : matchA) {
              const bool me2sig = isSignal(meP2.mass());
              const bool me2sb = isSideband(meP2.mass());

              fillPairHistosWeighted<part_pair_sig>(p1, meP2, wA);

              if ((se1sig && me2sb) || (se1sb && me2sig))
                fillPairHistosWeighted<part_pair_bkg_sigsb>(p1, meP2, wA);

              if (se1sb && me2sb)
                fillPairHistosWeighted<part_pair_bkg_sbsb>(p1, meP2, wA);
            }
          }
        }

        {
          std::vector<PoolTrack> matchB;
          for (auto const& meP : meVec1) {
            if (std::abs(meP.pt() - p1.pt()) < cMaxDeltaPt &&
                std::abs(meP.rap() - p1.rap()) < cMaxDeltaRap &&
                std::abs(RecoDecay::constrainAngle(meP.phi() - p1.phi(), -PI)) < cMaxDeltaPhi)
              matchB.push_back(meP);
          }
          if (!matchB.empty()) {
            const float wB = 1.0f / static_cast<float>(matchB.size());
            for (auto const& meP1 : matchB) {
              const bool me1sig = isSignal(meP1.mass());
              const bool me1sb = isSideband(meP1.mass());

              fillPairHistosWeighted<part_pair_sig>(meP1, p2, wB);

              if ((me1sig && se2sb) || (me1sb && se2sig))
                fillPairHistosWeighted<part_pair_bkg_sigsb>(meP1, p2, wB);

              if (me1sb && se2sb)
                fillPairHistosWeighted<part_pair_bkg_sbsb>(meP1, p2, wB);
            }
          }
        }

      } // trk2
    } // trk1
  } // analyzePairsMEKinematic

  // =========================================================================
  // O2 framework declarations
  // =========================================================================
  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaTracks = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;

  Preslice<LambdaTracks> perCollisionLambda = aod::lambdatrack::lambdaCollisionId;
  SliceCache cache;

  Partition<LambdaTracks> partLambdaTracks =
    (aod::lambdatrack::v0Type == (int8_t)kLambda) &&
    (aod::lambdatrackext::trueLambdaFlag == true) &&
    (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  Partition<LambdaTracks> partAntiLambdaTracks =
    (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) &&
    (aod::lambdatrackext::trueLambdaFlag == true) &&
    (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  void processDummy(LambdaCollisions::iterator const&) {}
  PROCESS_SWITCH(LambdaSpinPolarization, processDummy, "Dummy", false);

  void processDataReco(LambdaCollisions::iterator const& collision, LambdaTracks const&)
  {
    cent = collision.cent();
    auto lTrks = partLambdaTracks->sliceByCached(
      aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto alTrks = partAntiLambdaTracks->sliceByCached(
      aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzePairsWithMassWindow<kLambdaAntiLambda, kLambdaSBAntiLambda, kLambdaSBSBAntiLambda, false>(lTrks, alTrks);
    analyzePairsWithMassWindow<kAntiLambdaLambda, kAntiLambdaSBLambda, kAntiLambdaSBSBLambda, false>(alTrks, lTrks);
    analyzePairsWithMassWindow<kLambdaLambda, kLambdaSBLambda, kLambdaSBSBLambda, true>(lTrks, lTrks);
    analyzePairsWithMassWindow<kAntiLambdaAntiLambda, kAntiLambdaSBAntiLambda, kAntiLambdaSBSBAntiLambda, true>(alTrks, alTrks);
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataReco, "SE only (data/MCReco)", false);

  struct GetMultiplicity {
    float operator()(auto const& col) const { return col.cent(); }
  };
  using MixedBinning = FlexibleBinningPolicy<std::tuple<GetMultiplicity>,
                                             o2::aod::collision::PosZ,
                                             GetMultiplicity>;
  MixedBinning binningOnVtxAndMult{{GetMultiplicity{}}, {axisVtxZME, axisCentME}, true};

  void processDataRecoMixed(LambdaCollisions const& col, LambdaTracks const&)
  {
    for (auto const& [col1, col2] :
         soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (col1.globalIndex() == col2.globalIndex())
        continue;
      cent = col1.cent();
      histos.fill(HIST("QA/ME/hPoolCentVz"), col1.cent(), col1.posZ());

      auto lTrks1 = partLambdaTracks->sliceByCached(
        aod::lambdatrack::lambdaCollisionId, col1.globalIndex(), cache);
      auto lTrks2 = partLambdaTracks->sliceByCached(
        aod::lambdatrack::lambdaCollisionId, col2.globalIndex(), cache);
      auto alTrks1 = partAntiLambdaTracks->sliceByCached(
        aod::lambdatrack::lambdaCollisionId, col1.globalIndex(), cache);
      auto alTrks2 = partAntiLambdaTracks->sliceByCached(
        aod::lambdatrack::lambdaCollisionId, col2.globalIndex(), cache);

      histos.fill(HIST("QA/ME/hLambdaMultVsCent"), col1.cent(), lTrks1.size());
      histos.fill(HIST("QA/ME/hAntiLambdaMultVsCent"), col1.cent(), alTrks1.size());

      if (cMEMode == MEModeStandard) {
        analyzePairsME<kLambdaAntiLambda, kLambdaSBAntiLambda,
                       kLambdaSBSBAntiLambda, false>(lTrks1, alTrks2);
        analyzePairsME<kAntiLambdaLambda, kAntiLambdaSBLambda,
                       kAntiLambdaSBSBLambda, false>(alTrks1, lTrks2);
        analyzePairsME<kLambdaLambda, kLambdaSBLambda,
                       kLambdaSBSBLambda, false>(lTrks1, lTrks2);
        analyzePairsME<kAntiLambdaAntiLambda, kAntiLambdaSBAntiLambda,
                       kAntiLambdaSBSBAntiLambda, false>(alTrks1, alTrks2);

      } else if (cMEMode == MEModeKinematic) {
        analyzePairsMEKinematic<kLambdaAntiLambda, kLambdaSBAntiLambda,
                                kLambdaSBSBAntiLambda, false>(
          lTrks1, alTrks1, lTrks2, alTrks2);
        analyzePairsMEKinematic<kAntiLambdaLambda, kAntiLambdaSBLambda,
                                kAntiLambdaSBSBLambda, false>(
          alTrks1, lTrks1, alTrks2, lTrks2);
        analyzePairsMEKinematic<kLambdaLambda, kLambdaSBLambda,
                                kLambdaSBSBLambda, true>(
          lTrks1, lTrks1, lTrks2, lTrks2);
        analyzePairsMEKinematic<kAntiLambdaAntiLambda, kAntiLambdaSBAntiLambda,
                                kAntiLambdaSBSBAntiLambda, true>(
          alTrks1, alTrks1, alTrks2, alTrks2);
      }
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataRecoMixed,
                 "ME (modes 0=standard, 1=kinematicConstrained)", false);

  void processDataRecoMixEvent(LambdaCollisions::iterator const& collision,
                               LambdaTracks const&)
  {
    auto lTrks = partLambdaTracks->sliceByCached(
      aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto alTrks = partAntiLambdaTracks->sliceByCached(
      aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    if (lTrks.size() == 0 && alTrks.size() == 0)
      return;

    lambdaMixEvtCol(collision.index(), collision.cent(),
                    collision.posZ(), collision.timeStamp());

    for (auto const& track : lTrks)
      lambdaMixEvtTrk(collision.index(), track.globalIndex(),
                      track.px(), track.py(), track.pz(), track.mass(),
                      track.prPx(), track.prPy(), track.prPz(),
                      track.v0Type(), collision.timeStamp());
    for (auto const& track : alTrks)
      lambdaMixEvtTrk(collision.index(), track.globalIndex(),
                      track.px(), track.py(), track.pz(), track.mass(),
                      track.prPx(), track.prPy(), track.prPz(),
                      track.v0Type(), collision.timeStamp());
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataRecoMixEvent,
                 "Mix-event table filling", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LambdaTableProducer>(cfgc),
    adaptAnalysisTask<LambdaTracksExtProducer>(cfgc),
    adaptAnalysisTask<LambdaSpinPolarization>(cfgc)};
}
