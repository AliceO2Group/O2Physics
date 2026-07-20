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

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TList.h>
#include <TObject.h>
#include <TPDGCode.h>

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstddef>
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
DECLARE_SOA_COLUMN(BField, bField, float);
} // namespace lambdacollision
DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent,
                  lambdacollision::Mult,
                  aod::collision::PosX,
                  aod::collision::PosY,
                  aod::collision::PosZ,
                  lambdacollision::BField);
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
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(PosTrackId, posTrackId, int64_t);
DECLARE_SOA_COLUMN(NegTrackId, negTrackId, int64_t);
DECLARE_SOA_COLUMN(PosTrackKin, posTrackKin, std::vector<float>);
DECLARE_SOA_COLUMN(NegTrackKin, negTrackKin, std::vector<float>);
DECLARE_SOA_COLUMN(PartType, partType, int8_t);
DECLARE_SOA_COLUMN(CorrFact, corrFact, float);
} // namespace lambdatrack
DECLARE_SOA_TABLE(LambdaTracks, "AOD", "LAMBDATRACKS", o2::soa::Index<>,
                  lambdatrack::LambdaCollisionId,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Px,
                  lambdatrack::Py,
                  lambdatrack::Pz,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::PosTrackKin,
                  lambdatrack::NegTrackKin,
                  lambdatrack::PartType,
                  lambdatrack::CorrFact);
using LambdaTrack = LambdaTracks::iterator;

namespace kaontrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaCollision, lambdaCollision);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Rap, rap, float);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(KaonTrackId, kaonTrackId, int64_t);
DECLARE_SOA_COLUMN(PartType, partType, int8_t);
DECLARE_SOA_COLUMN(CorrFact, corrFact, float);
} // namespace kaontrack
DECLARE_SOA_TABLE(KaonTracks, "AOD", "KAONTRACKS", o2::soa::Index<>,
                  kaontrack::LambdaCollisionId,
                  kaontrack::Pt,
                  kaontrack::Eta,
                  kaontrack::Phi,
                  kaontrack::Rap,
                  kaontrack::Px,
                  kaontrack::Py,
                  kaontrack::Pz,
                  kaontrack::Mass,
                  kaontrack::KaonTrackId,
                  kaontrack::PartType,
                  kaontrack::CorrFact);
using KaonTrack = KaonTracks::iterator;

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

namespace kaontrackext
{
DECLARE_SOA_COLUMN(KaonSharingLambdaDau, kaonSharingLambdaDau, bool);
DECLARE_SOA_COLUMN(KaonSharingLambdaDauIds, kaonSharingLambdaDauIds, std::vector<int64_t>);
DECLARE_SOA_COLUMN(TrueKaonFlag, trueKaonFlag, bool);
} // namespace kaontrackext
DECLARE_SOA_TABLE(KaonTracksExt, "AOD", "KAONTRACKSEXT",
                  kaontrackext::KaonSharingLambdaDau,
                  kaontrackext::KaonSharingLambdaDauIds,
                  kaontrackext::TrueKaonFlag);
using KaonTrackExt = KaonTracksExt::iterator;

namespace lambdamcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMcGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(LambdaMcGenTracks, "AOD", "LMCGENTRACKS", o2::soa::Index<>,
                  lambdamcgentrack::LambdaMcGenCollisionId,
                  lambdatrack::Pt,
                  lambdatrack::Eta,
                  lambdatrack::Phi,
                  lambdatrack::Rap,
                  lambdatrack::Px,
                  lambdatrack::Py,
                  lambdatrack::Pz,
                  lambdatrack::Mass,
                  lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId,
                  lambdatrack::PartType,
                  lambdatrack::CorrFact);
using LambdaMcGenTrack = LambdaMcGenTracks::iterator;

namespace kaonmcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMcGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(KaonMcGenTracks, "AOD", "KMCGENTRACKS", o2::soa::Index<>,
                  kaonmcgentrack::LambdaMcGenCollisionId,
                  kaontrack::Pt,
                  kaontrack::Eta,
                  kaontrack::Phi,
                  kaontrack::Rap,
                  kaontrack::Px,
                  kaontrack::Py,
                  kaontrack::Pz,
                  kaontrack::Mass,
                  kaontrack::KaonTrackId,
                  kaontrack::PartType,
                  kaontrack::CorrFact);
using KaonMcGenTrack = KaonMcGenTracks::iterator;
} // namespace o2::aod

enum CollisionLabels {
  kTotColBeforeHasMcCollision = 1,
  kTotCol,
  kPassSelCol
};

enum LambdaLabels {
  kLambdaBeforeHasMcParticle = 1,
  kAllV0Tracks,
  kPassV0KinCuts,
  kPassV0TopoSel,
  kPassK0ShortMassRej,
  kPassV0DauTrackSel,
  kAllSelPassed,
  kGenTotAccLambda,
  kGenLambdaNoDau,
};

enum KaonLabels {
  kKaonAllChargedTracks = 1,
  kKaonPassKinSel,
  kKaonPassGlobalSel,
  kKaonPassDcaSel,
  kKaonPassElRejSel,
  kKaonPassAllSel
};

enum EffCorrType {
  kEffCorrPtCent = 0,
  kEffCorrPtRapCent
};

enum CentEstType {
  kCentFT0M = 0,
  kCentFT0C
};

enum ParticleType {
  kLambda = 0,
  kAntiLambda,
  kKaonPlus,
  kKaonMinus
};

enum ParticlePairType {
  kLambdaAntiLambda = 0,
  kLambdaLambda,
  kAntiLambdaAntiLambda,
  kLambdaKaonPlus,
  kLambdaKaonMinus,
  kAntiLambdaKaonPlus,
  kAntiLambdaKaonMinus
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
  // Table Producers
  Produces<aod::LambdaCollisions> lambdaCollisionTable;
  Produces<aod::LambdaTracks> lambdaTrackTable;
  Produces<aod::KaonTracks> kaonTrackTable;
  Produces<aod::LambdaMcGenCollisions> lambdaMCGenCollisionTable;
  Produces<aod::LambdaMcGenTracks> lambdaMCGenTrackTable;
  Produces<aod::KaonMcGenTracks> kaonMCGenTrackTable;

  // Centrality Axis
  ConfigurableAxis cCentBins{"cCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Variable Centrality Bins"};

  // Collisions
  Configurable<int> cCentEstimator{"cCentEstimator", 1, "Centrality Estimator : 0-FT0M, 1-FT0C"};
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinCent{"cMinCent", 0., "Minumum Centrality"};
  Configurable<float> cMaxCent{"cMaxCent", 100.0, "Maximum Centrality"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", true, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.1, "p_{T} minimum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "TPC Min Crossed Rows"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 3.0, "TPC NSigma Selection Cut"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  // Kaon Tracks
  Configurable<float> cKaonMinPt{"cKaonMinPt", 0.3, "Kaon Min pT"};
  Configurable<float> cKaonMaxPt{"cKaonMaxPt", 2.4, "Kaon Max pT"};
  Configurable<float> cKaonRapCut{"cKaonRapCut", 0.5, "Kaon |y| cut"};
  Configurable<bool> cKaonGlobalSel{"cKaonGlobalSel", true, "Global Track"};
  Configurable<float> cKaonDcaXYCut{"cKaonDcaXYCut", 0.1, "DcaXY Cut"};
  Configurable<float> cKaonDcaZCut{"cKaonDcaZCut", 1., "DcaZ Cut"};
  Configurable<float> cTpcElRejCutMin{"cTpcElRejCutMin", -3., "Electron Rejection Cut Minimum"};
  Configurable<float> cTpcElRejCutMax{"cTpcElRejCutMax", 5., "Electron Rejection Cut Maximum"};
  Configurable<float> cKaonTpcNSigmaCut{"cKaonTpcNSigmaCut", 2, "TPC Kaon NSigma Cut"};
  Configurable<float> cTpcRejCut{"cTpcRejCut", 3, "TPC Rej Cut"};
  Configurable<float> cKaonTofNSigmaCut{"cKaonTofNSigmaCut", 2, "TOF Kaon NSigma Cut"};
  Configurable<float> cTofRejCut{"cTofRejCut", 3, "TOF Rej Cut"};
  Configurable<float> cKaonTpcPtSel{"cKaonTpcPtSel", 0.7, "Kaon TPC pT cutoff"};

  // V0s
  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.02, "Minimum Proton DCAr to PV"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.06, "Minimum Pion DCAr to PV"};
  Configurable<double> cMaxDcaV0Daughters{"cMaxDcaV0Daughters", 1., "Maximum DCA between V0 daughters"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 0.1, "Maximum DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Maximum ctau"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.01, "Reject K0Short Candidates"};
  Configurable<float> cArmPodSel{"cArmPodSel", 0.2, "Armentros-Podolanski Selection"};

  // DCA
  Configurable<int> cWeakDecayProcess{"cWeakDecayProcess", 4, "Weak decay process tag for DCA analysis"};

  // V0s acceptance
  Configurable<float> cLambdaMassWindow{"cLambdaMassWindow", 0.007, "Lambda Mass Window"};
  Configurable<float> cLambdaMinPt{"cLambdaMinPt", 0.8, "Minimum Lambda pT"};
  Configurable<float> cLambdaMaxPt{"cLambdaMaxPt", 3.8, "Maximum Lambda pT"};
  Configurable<float> cLambdaRapCut{"cLambdaRapCut", 0.5, "Lambda |rap| cut"};
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};

  // Efficiency Correction
  Configurable<bool> cGetCorrectionFlag{"cGetCorrectionFlag", false, "Apply correction flag"};
  Configurable<bool> cGetEffFact{"cGetEffFact", false, "Get Efficiency Factor Flag"};
  Configurable<int> cCorrFactHist{"cCorrFactHist", 0, "Efficiency Factor Histogram"};

  // MC
  Configurable<bool> cSelPrimaryParticle{"cSelPrimaryParticle", false, "Select primary particle"};
  Configurable<bool> cSelTrueParticle{"cSelTrueParticle", false, "Select MC PDG matched particle"};

  // CCDB
  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://alice-ccdb.cern.ch", "ALICE CCDB URL"};
  Configurable<std::string> cPathCCDBRecoEff{"cPathCCDBRecoEff", "Users/y/ypatley/LBF/PP/RecoEfficiency", "Path for ccdb-object for reco efficiency"};
  Configurable<int64_t> nolaterthan{"nolaterthan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::parameters::GRPMagField* grpo = nullptr;

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Initialize corr_factor objects
  std::vector<std::vector<std::string>> vCorrFactStrings = {{"hEffVsPtCentLambda", "hEffVsPtCentAntiLambda", "hEffVsPtCentKaonPlus", "hEffVsPtCentKaonMinus"}, {"hEffVsPtYCentLambda", "hEffVsPtYCentAntiLambda", "hEffVsPtYCentKaonPlus", "hEffVsPtYCentKaonMinus"}};

  // Store correction histograms
  struct CorrHist {
    std::array<TH2F*, 4> vPtCentCorrHists{};
    std::array<TH3F*, 4> vPtRapCentCorrHists{};
  } corrHist;

  // Initialize Global Variables
  float cent = 0., mult = 0.;
  TList* ccdbObjRecoEff = nullptr;
  TList* ccdbObjMatchEff = nullptr;
  static constexpr auto SubDir = std::array{"QA/Lambda/", "QA/AntiLambda/", "QA/KaonPlus/", "QA/KaonMinus/"};

  void init(InitContext const&)
  {
    // initialize axis specifications
    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(15, 0.5, 15.5, "");
    const AxisSpec axisEffChecks(3, 0, 3, "");
    const AxisSpec axisCent(100, 0, 100, "Centrality(%)");
    const AxisSpec axisVarCent(cCentBins, "FT0C%");
    const AxisSpec axisPVMults(1000, 0, 1000, "N_{PV}");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(140, 1.08, 1.15, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisV0Pt(100., 0., 10., "p_{T} (GeV/#it{c})");
    const AxisSpec axisV0Rap(48, -1.2, 1.2, "y");
    const AxisSpec axisV0Eta(48, -1.2, 1.2, "#eta");
    const AxisSpec axisV0Phi(36, 0., TwoPI, "#phi (rad)");

    const AxisSpec axisRadius(200, 0, 200, "r(cm)");
    const AxisSpec axisCosPA(100, 0.99, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaV0PV(100, 0., 0.1, "dca (cm)");
    const AxisSpec axisDcaProngPV(5000, -50., 50., "dca (cm)");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daug DCA (#sigma)");
    const AxisSpec axisCTau(80, 0, 40, "c#tau (cm)");
    const AxisSpec axisAlpha(40, -1, 1, "#alpha");
    const AxisSpec axisQtarm(40, 0, 0.4, "q_{T}");

    const AxisSpec axisDcaLambda(100, 0., cMaxDcaV0ToPV, "Dca_{V^{0}}");
    const AxisSpec axisDcaKaon(100, -cKaonDcaXYCut, cKaonDcaXYCut, "Dca_{K}");

    const AxisSpec axisITSTPCTrackPt(100, 0, 10, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackPt(40, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackDCA(200, -1, 1, "dca_{XY} (cm)");
    const AxisSpec axisMomPID(80, 0, 4, "p_{T} (GeV/#it{c})");
    const AxisSpec axisTrackNsigma(401, -10.025, 10.025, {"n#sigma"});
    const AxisSpec axisTrackdEdx(360, 20, 200, "#frac{dE}{dx}");
    const AxisSpec axisTrackTofSignal(240, 0, 1.2, "#beta");

    // Create Histograms.
    // Event histograms
    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCols});
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
    histos.add("Events/h1f_collision_cent", "Centrality distribution", kTH1F, {axisCent});

    // QA
    histos.add("Tracks/h1f_lambda_info", "Lambda selection info", kTH1F, {axisTrks});
    histos.add("Tracks/h1f_kaon_info", "Kaon selection info", kTH1F, {axisTrks});
    histos.add("Tracks/h1f_effcorr_info", "Efficiency correction info", kTH1F, {axisEffChecks});
    histos.add("Tracks/h2f_armpod_before_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("Tracks/h2f_armpod_after_sel", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});

    // QA Lambda
    histos.add("QA/Lambda/h3f_centmasspt", "Invariant Mass", kTH3F, {axisCent, axisV0Pt, axisV0Mass});
    histos.add("QA/Lambda/h2f_qt_vs_alpha", "Armentros-Podolanski Plot", kTH2F, {axisAlpha, axisQtarm});
    histos.add("QA/Lambda/h1f_dca_V0_daughters", "DCA between V0 daughters", kTH1F, {axisDcaDau});
    histos.add("QA/Lambda/h1f_dca_pos_to_PV", "DCA positive prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_neg_to_PV", "DCA negative prong to PV", kTH1F, {axisDcaProngPV});
    histos.add("QA/Lambda/h1f_dca_V0_to_PV", "DCA V0 to PV", kTH1F, {axisDcaV0PV});
    histos.add("QA/Lambda/h1f_V0_cospa", "cos(#theta_{PA})", kTH1F, {axisCosPA});
    histos.add("QA/Lambda/h1f_V0_radius", "V_{0} Decay Radius in XY plane", kTH1F, {axisRadius});
    histos.add("QA/Lambda/h1f_V0_ctau", "V_{0} c#tau", kTH1F, {axisCTau});
    histos.add("QA/Lambda/h2f_pos_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_neg_prong_dcaXY_vs_pt", "DCA vs p_{T}", kTH2F, {axisTrackPt, axisTrackDCA});
    histos.add("QA/Lambda/h2f_pos_prong_dEdx_vs_p", "TPC Signal Pos-Prong", kTH2F, {axisMomPID, axisTrackdEdx});
    histos.add("QA/Lambda/h2f_neg_prong_dEdx_vs_p", "TPC Signal Neg-Prong", kTH2F, {axisMomPID, axisTrackdEdx});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisTrackNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pr_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisTrackNsigma});
    histos.add("QA/Lambda/h2f_pos_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Pos Prong", kTH2F, {axisMomPID, axisTrackNsigma});
    histos.add("QA/Lambda/h2f_neg_prong_tpc_nsigma_pi_vs_p", "TPC n#sigma Neg Prong", kTH2F, {axisMomPID, axisTrackNsigma});

    // Dca analysis histogram
    histos.add("QA/Lambda/DCA/h3f_Dca_vs_pT_cent", "DCA", kTH3F, {axisCent, axisV0Pt, axisDcaLambda});
    histos.add("QA/KaonPlus/DCA/h3f_Dca_vs_pT_cent", "DCA", kTH3F, {axisCent, axisV0Pt, axisDcaKaon});

    // QA Kaons
    histos.add("QA/KaonPlus/hdEdX", "dE/dx vs pT", kTH2F, {axisMomPID, axisTrackdEdx});
    histos.add("QA/KaonPlus/hTOFSignal", "#beta_{TOF} vs p_{T}", kTH2F, {axisMomPID, axisTrackTofSignal});
    histos.add("QA/KaonPlus/hTPCNSigma", "n#sigma_{TPC} vs p_{T}", kTH2F, {axisMomPID, axisTrackNsigma});
    histos.add("QA/KaonPlus/hTOFNSigma", "n#sigma_{TOF} vs p_{T}", kTH2F, {axisMomPID, axisTrackNsigma});

    // MC Generated Histograms
    if (doprocessMCRecoGen || doprocessMCReco) {
      // DCA analysis
      histos.add("QA/Lambda/DCA/h3f_Prm_Dca_vs_pT_cent", "Primary DCA", kTH3F, {axisCent, axisV0Pt, axisDcaLambda});
      histos.add("QA/Lambda/DCA/h3f_Scd_Dca_vs_pT_cent", "Weak Decay DCA", kTH3F, {axisCent, axisV0Pt, axisDcaLambda});
      histos.add("QA/Lambda/DCA/h3f_Mat_Dca_vs_pT_cent", "Material DCA", kTH3F, {axisCent, axisV0Pt, axisDcaLambda});
      histos.add("QA/KaonPlus/DCA/h3f_Prm_Dca_vs_pT_cent", "Primary DCA", kTH3F, {axisCent, axisV0Pt, axisDcaKaon});
      histos.add("QA/KaonPlus/DCA/h3f_Scd_Dca_vs_pT_cent", "Weak Decay DCA", kTH3F, {axisCent, axisV0Pt, axisDcaKaon});
      histos.add("QA/KaonPlus/DCA/h3f_Mat_Dca_vs_pT_cent", "Material DCA", kTH3F, {axisCent, axisV0Pt, axisDcaKaon});

      // McGen Histos
      histos.add("McGen/h1f_collision_recgen", "# of Reco Collision Associated to One Mc Generator Collision", kTH1F, {axisMult});
      histos.add("McGen/h1f_collisions_info", "# of collisions", kTH1F, {axisCols});
      histos.add("McGen/h2f_collision_posZ", "V_{z}-distribution", kTH2F, {axisVz, axisVz});
      histos.add("McGen/h2f_collision_cent", "FT0M Centrality", kTH2F, {axisCent, axisCent});
      histos.add("McGen/h1f_lambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});
      histos.add("McGen/h1f_antilambda_daughter_PDG", "PDG Daughters", kTH1F, {axisPID});

      // set bin lables specific to MC
      histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotColBeforeHasMcCollision, "kTotColBeforeHasMcCollision");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
      histos.get<TH1>(HIST("McGen/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
      histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kLambdaBeforeHasMcParticle, "kLambdaBeforeHasMcParticle");
      histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kGenTotAccLambda, "kGenTotAccLambda");
      histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kGenLambdaNoDau, "kGenLambdaNoDau");
    }

    // Clone histograms
    // QA Anti-Lambda
    histos.addClone("QA/Lambda/", "QA/AntiLambda/");

    // QA KaonMinus
    histos.addClone("QA/KaonPlus/", "QA/KaonMinus/");

    // Set bin labels
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kTotCol, "kTotCol");
    histos.get<TH1>(HIST("Events/h1f_collisions_info"))->GetXaxis()->SetBinLabel(CollisionLabels::kPassSelCol, "kPassSelCol");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kAllV0Tracks, "kAllV0Tracks");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kPassK0ShortMassRej, "kPassK0ShortMassRej");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kPassV0DauTrackSel, "kPassV0DauTrackSel");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kPassV0KinCuts, "kPassV0KinCuts");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kPassV0TopoSel, "kPassV0TopoSel");
    histos.get<TH1>(HIST("Tracks/h1f_lambda_info"))->GetXaxis()->SetBinLabel(LambdaLabels::kAllSelPassed, "kAllSelPassed");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonAllChargedTracks, "kKaonAllChargedTracks");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonPassKinSel, "kKaonPassKinSel");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonPassGlobalSel, "kKaonPassGlobalSel");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonPassDcaSel, "kKaonPassDcaSel");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonPassElRejSel, "kKaonPassElRejSel");
    histos.get<TH1>(HIST("Tracks/h1f_kaon_info"))->GetXaxis()->SetBinLabel(KaonLabels::kKaonPassAllSel, "kKaonPassAllSel");

    // Load correction factor
    if (cGetCorrectionFlag) {
      // Set CCDB url
      ccdb->setURL(cUrlCCDB.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(nolaterthan.value);

      // Get CCDB object
      ccdbObjRecoEff = ccdb->getForTimeStamp<TList>(cPathCCDBRecoEff.value, nolaterthan.value);

      // Load reco eff corrections
      LoadRecoEfficiencyHistograms();
    }
  }

  template <typename T>
  void GetCorrFactHists(T& vhists, std::vector<std::string> const& strings)
  {
    using HistPtr = typename T::value_type;
    for (size_t i = 0; i < vhists.size(); ++i) {
      auto* obj = ccdbObjRecoEff->FindObject(strings[i].c_str());
      if (!obj) {
        LOGF(fatal, "CCDB object %s not found!", strings[i].c_str());
        continue;
      }
      vhists[i] = dynamic_cast<HistPtr>(obj->Clone());
      if (!vhists[i]) {
        LOGF(fatal, "CCDB object %s could not be cast!", strings[i].c_str());
      }
    }
  }

  // Load reco efficiency histograms
  void LoadRecoEfficiencyHistograms()
  {
    if (cCorrFactHist == kEffCorrPtCent) {
      GetCorrFactHists(corrHist.vPtCentCorrHists, vCorrFactStrings[kEffCorrPtCent]);
    } else if (cCorrFactHist == kEffCorrPtRapCent) {
      GetCorrFactHists(corrHist.vPtRapCentCorrHists, vCorrFactStrings[kEffCorrPtRapCent]);
    }
  }

  // Get magnetic field
  float getMagneticField(int64_t const& timestamp)
  {
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 1;
      }
    }
    auto field = std::lround(5.f * grpo->getL3Current() / 30000.f);
    return 0.1 * field;
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    // Vz Selection
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }

    // Run 3 Min-Bias Trigger
    if (cSel8Trig && !col.sel8()) {
      return false;
    }
    if (cCentEstimator == kCentFT0M) {
      cent = col.centFT0M();
    } else if (cCentEstimator == kCentFT0C) {
      cent = col.centFT0C();
    }

    if (cent <= cMinCent || cent >= cMaxCent) { // select centrality percentile class
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
    mult = col.multTPC();

    return true;
  }

  // Ambiguous track rejection
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

  template <typename V, typename T>
  bool selLambdaDauTracks(V const& v0, T const& postrack, T const& negtrack, ParticleType& partType)
  {
    // Kinematic selection
    if (postrack.pt() <= cTrackMinPt || negtrack.pt() <= cTrackMinPt || std::abs(postrack.eta()) >= cTrackEtaCut || std::abs(negtrack.eta()) >= cTrackEtaCut || postrack.tpcNClsCrossedRows() <= cMinTpcCrossedRows || negtrack.tpcNClsCrossedRows() <= cMinTpcCrossedRows) {
      return false;
    }

    // Apply DCA Selection on Daughter Tracks Based on Lambda/AntiLambda daughters
    if (partType == kLambda) {
      if (std::abs(v0.dcapostopv()) <= cMinDcaProtonToPV || std::abs(v0.dcanegtopv()) <= cMinDcaPionToPV) {
        return false;
      }
    } else if (partType == kAntiLambda) {
      if (std::abs(v0.dcapostopv()) <= cMinDcaPionToPV || std::abs(v0.dcanegtopv()) <= cMinDcaProtonToPV) {
        return false;
      }
    } else {
      return false;
    }

    // Daughter track PID [Dau1 = PosTrack, Dau2 = NegTrack]
    float tpcNSigmaDau1 = 0., tpcNSigmaDau2 = 0.;

    if (partType == kLambda) {
      tpcNSigmaDau1 = postrack.tpcNSigmaPr();
      tpcNSigmaDau2 = negtrack.tpcNSigmaPi();
    } else if (partType == kAntiLambda) {
      tpcNSigmaDau1 = postrack.tpcNSigmaPi();
      tpcNSigmaDau2 = negtrack.tpcNSigmaPr();
    } else {
      LOGF(fatal, "Particle call not Lambda !!!");
      return false;
    }

    if (std::abs(tpcNSigmaDau1) >= cTpcNsigmaCut || std::abs(tpcNSigmaDau2) >= cTpcNsigmaCut) {
      return false;
    }

    return true;
  }

  template <typename C, typename V, typename T>
  bool selLambda(C const& collision, V const& v0, T const& tracks, ParticleType& partType)
  {
    // Initialize daughter tracks
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // Kinematic selections
    if (v0.pt() <= cLambdaMinPt || v0.pt() >= cLambdaMaxPt || std::abs(v0.yLambda()) >= cLambdaRapCut) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_lambda_info"), kPassV0KinCuts);

    // Decay length
    float ctauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * MassLambda0;

    // Topological selections
    if (v0.dcaV0daughters() >= cMaxDcaV0Daughters || v0.dcav0topv() >= cMaxDcaV0ToPV || v0.v0radius() <= cMinV0TransRadius || v0.v0Type() != cV0TypeSelection || ctauLambda >= cMaxV0CTau || v0.v0cosPA() <= cMinV0CosPA) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_lambda_info"), kPassV0TopoSel);

    // K0s mass rejection
    if (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_lambda_info"), kPassK0ShortMassRej);

    // Select Lambda / AntiLambda based on Armenteros-Podolanski
    bool lambdaFlag = false, antiLambdaFlag = false;
    float mass = 0.;

    if (v0.qtarm() < cArmPodSel * std::abs(v0.alpha()) && v0.alpha() > 0) {
      lambdaFlag = true;
      partType = kLambda;
      mass = v0.mLambda();
    } else {
      antiLambdaFlag = true;
      partType = kAntiLambda;
      mass = v0.mAntiLambda();
    }

    // Select Lambda daughters
    if (!selLambdaDauTracks(v0, postrack, negtrack, partType)) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_lambda_info"), kPassV0DauTrackSel);

    // Remove lambda with ambiguous daughters
    if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks)) {
      return false;
    }

    // Fill invariant mass histograms after all selection
    if (lambdaFlag) {
      histos.fill(HIST("QA/Lambda/h3f_centmasspt"), cent, v0.pt(), mass);
    } else {
      histos.fill(HIST("QA/AntiLambda/h3f_centmasspt"), cent, v0.pt(), mass);
    }

    // Mass window selection
    if ((lambdaFlag && std::abs(v0.mLambda() - MassLambda0) >= cLambdaMassWindow) || (antiLambdaFlag && std::abs(v0.mAntiLambda() - MassLambda0) >= cLambdaMassWindow)) {
      return false;
    }

    // All Selection Criterion Passed
    return true;
  }

  template <typename T>
  bool selKaonTrack(T const& track, float const& rap)
  {
    // Kinematic selection
    if (track.pt() <= cKaonMinPt || track.pt() >= cKaonMaxPt || std::abs(rap) >= cKaonRapCut) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonPassKinSel);

    // Global track selection
    if (cKaonGlobalSel && !track.isGlobalTrackWoDCA()) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonPassGlobalSel);

    // Dca selection
    if (std::abs(track.dcaXY()) >= cKaonDcaXYCut || std::abs(track.dcaZ()) >= cKaonDcaZCut) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonPassDcaSel);

    // Electron rejection
    if (std::abs(track.tpcNSigmaPi()) > cTpcRejCut && std::abs(track.tpcNSigmaKa()) > cTpcRejCut && std::abs(track.tpcNSigmaPr()) > cTpcRejCut && track.tpcNSigmaEl() > cTpcElRejCutMin && track.tpcNSigmaEl() < cTpcElRejCutMax) {
      return false;
    }

    histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonPassElRejSel);

    // Kaon PID TPC + TOF
    if (track.hasTOF()) {
      if (std::abs(track.tofNSigmaKa()) >= cKaonTofNSigmaCut || std::abs(track.tofNSigmaPi()) < cTofRejCut || std::abs(track.tofNSigmaPr()) < cTofRejCut || std::abs(track.tpcNSigmaKa()) >= cKaonTpcNSigmaCut) {
        return false;
      }
    } else {
      if (track.pt() >= cKaonTpcPtSel || std::abs(track.tpcNSigmaKa()) >= cKaonTpcNSigmaCut || std::abs(track.tpcNSigmaPi()) < cTpcRejCut || std::abs(track.tpcNSigmaPr()) < cTpcRejCut) {
        return false;
      }
    }

    histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonPassAllSel);

    return true;
  }

  // Correction factors
  template <ParticleType part, typename V>
  float getCorrectionFactors(V const& v, float const& rap)
  {
    // Initialize efficiency factor
    float effCorrFact = 1.;

    // Get Efficiency Factor
    if (cGetEffFact) {
      if (cCorrFactHist == kEffCorrPtCent) {
        histos.fill(HIST("Tracks/h1f_effcorr_info"), kEffCorrPtCent);
        effCorrFact = corrHist.vPtCentCorrHists[part]->GetBinContent(corrHist.vPtCentCorrHists[part]->FindBin(cent, v.pt()));
      } else if (cCorrFactHist == kEffCorrPtRapCent) {
        histos.fill(HIST("Tracks/h1f_effcorr_info"), kEffCorrPtRapCent);
        effCorrFact = corrHist.vPtRapCentCorrHists[part]->GetBinContent(corrHist.vPtRapCentCorrHists[part]->FindBin(cent, v.pt(), rap));
      } else {
        LOGF(warning, "CCDB OBJECT IS NOT A HISTOGRAM !!!");
      }
    }

    // Return correction factor
    return effCorrFact;
  }

  // Lambda QA
  template <ParticleType part, typename C, typename V, typename T>
  void fillLambdaQAHistos(C const& col, V const& v0, T const&)
  {
    // Daugthers
    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    // Decay length
    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;

    histos.fill(HIST(SubDir[part]) + HIST("h2f_qt_vs_alpha"), v0.alpha(), v0.qtarm());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_daughters"), v0.dcaV0daughters());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_pos_to_PV"), v0.dcapostopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_neg_to_PV"), v0.dcanegtopv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_dca_V0_to_PV"), v0.dcav0topv());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_cospa"), v0.v0cosPA());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_radius"), v0.v0radius());
    histos.fill(HIST(SubDir[part]) + HIST("h1f_V0_ctau"), ctau);
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dcaXY_vs_pt"), postrack.pt(), postrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dcaXY_vs_pt"), negtrack.pt(), negtrack.dcaXY());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_dEdx_vs_p"), postrack.tpcInnerParam(), postrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_dEdx_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pr_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pr_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPr());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_pos_prong_tpc_nsigma_pi_vs_p"), postrack.tpcInnerParam(), postrack.tpcNSigmaPi());
    histos.fill(HIST(SubDir[part]) + HIST("h2f_neg_prong_tpc_nsigma_pi_vs_p"), negtrack.tpcInnerParam(), negtrack.tpcNSigmaPi());
  }

  // Kaon QA
  template <ParticleType part, typename T>
  void fillKaonQA(T const& track)
  {
    histos.fill(HIST(SubDir[part]) + HIST("hdEdX"), track.pt(), track.tpcSignal());
    histos.fill(HIST(SubDir[part]) + HIST("hTPCNSigma"), track.pt(), track.tpcNSigmaKa());
    if (track.hasTOF()) {
      histos.fill(HIST(SubDir[part]) + HIST("hTOFSignal"), track.pt(), track.beta());
      histos.fill(HIST(SubDir[part]) + HIST("hTOFNSigma"), track.pt(), track.tofNSigmaKa());
    }
  }

  // Dca analysis
  template <ParticleType part, DMCType dmc, typename T>
  void getDcaHist(T const& track, float const& dca)
  {
    // Fill histograms
    histos.fill(HIST(SubDir[part]) + HIST("DCA/h3f_Dca_vs_pT_cent"), cent, track.pt(), dca);

    if constexpr (dmc == kMC) {
      auto mcpart = track.mcParticle();
      if (mcpart.isPhysicalPrimary()) {
        histos.fill(HIST(SubDir[part]) + HIST("DCA/h3f_Prm_Dca_vs_pT_cent"), cent, track.pt(), dca);
      } else {
        if (mcpart.getProcess() == cWeakDecayProcess) { // Weak decay
          histos.fill(HIST(SubDir[part]) + HIST("DCA/h3f_Scd_Dca_vs_pT_cent"), cent, track.pt(), dca);
        } else { // Material
          histos.fill(HIST(SubDir[part]) + HIST("DCA/h3f_Mat_Dca_vs_pT_cent"), cent, track.pt(), dca);
        }
      }
    }
  }

  // Reconstructed Level Tables
  template <DMCType dmc, typename C, typename B, typename V, typename T>
  void fillLambdaRecoTables(C const& collision, B const&, V const& v0tracks, T const& tracks)
  {
    // Total Collisions
    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    // Select Collision (Only for Data... McRec has been selected already !!!)
    if constexpr (dmc == kData) {
      if (!selCollision(collision)) {
        return;
      }
    }

    // Fill Collision Histograms
    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());
    histos.fill(HIST("Events/h1f_collision_cent"), cent);

    // Fill Collision Table
    float magField = getMagneticField(collision.template foundBC_as<BCsRun3>().timestamp());
    lambdaCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ(), magField);

    // initialize v0track objects
    ParticleType partType = kLambda;
    float lambdaMass = 0., lambdaCorrFact = 1.;

    // Loop over V0s to select Lambda
    for (auto const& v0 : v0tracks) {
      // Daugthers
      auto postrack = v0.template posTrack_as<T>();
      auto negtrack = v0.template negTrack_as<T>();

      // Check for corresponding MCGen Particle
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_lambda_info"), kLambdaBeforeHasMcParticle);
        if (!v0.has_mcParticle() || !postrack.has_mcParticle() || !negtrack.has_mcParticle()) { // check corresponding MC particle
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_lambda_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      // Select V0 as Lambda/AntiLambda
      if (!selLambda(collision, v0, tracks, partType)) {
        continue;
      }

      // We have v0 as lambda
      histos.fill(HIST("Tracks/h1f_lambda_info"), kAllSelPassed);
      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());

      // MC matching
      if constexpr (dmc == kMC) {
        auto mcpart = v0.mcParticle();
        if (cSelPrimaryParticle && !mcpart.isPhysicalPrimary()) { // Primary lambda selection
          continue;
        }

        if (cSelTrueParticle && std::abs(mcpart.pdgCode()) != kLambda0) { // True lambda selection
          continue;
        }
      }

      // Get Lambda mass and correction factor
      lambdaMass = (partType == kLambda) ? v0.mLambda() : v0.mAntiLambda();
      if (cGetCorrectionFlag) {
        lambdaCorrFact = (partType == kLambda) ? getCorrectionFactors<kLambda>(v0, v0.yLambda()) : getCorrectionFactors<kAntiLambda>(v0, v0.yLambda());
      }

      // Fill lambda qa
      if (partType == kLambda) {
        getDcaHist<kLambda, dmc>(v0, v0.dcav0topv());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
      } else {
        getDcaHist<kAntiLambda, dmc>(v0, v0.dcav0topv());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
      }

      // Daughter kinematics
      std::vector<float> posTrackKin = {postrack.pt(), postrack.eta(), postrack.phi()};
      std::vector<float> negTrackKin = {negtrack.pt(), negtrack.eta(), negtrack.phi()};

      // Fill Lambda/AntiLambda Table
      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0.pt(), v0.eta(), v0.phi(), v0.yLambda(), v0.px(), v0.py(), v0.pz(), lambdaMass,
                       v0.template posTrack_as<T>().index(), v0.template negTrack_as<T>().index(), posTrackKin, negTrackKin,
                       (int8_t)partType, lambdaCorrFact);
    }

    // Loop over tracks to select Kaon
    float kaonCorrFact = 1.;
    for (auto const& track : tracks) {
      // Check corresponding MC particle
      if constexpr (dmc == kMC) {
        if (!track.has_mcParticle()) {
          continue;
        }
      }

      // All charged tracks
      histos.fill(HIST("Tracks/h1f_kaon_info"), kKaonAllChargedTracks);

      // Kaon rapidity
      std::array<float, 3> mom = {track.px(), track.py(), track.pz()};
      float rap = RecoDecay::y(mom, MassKPlus);
      if (!selKaonTrack(track, rap)) { // Kaon selection
        continue;
      }

      // MC matching
      if constexpr (dmc == kMC) {
        auto mcpart = track.mcParticle();
        if (cSelPrimaryParticle && !mcpart.isPhysicalPrimary()) { // Primary kaon selection
          continue;
        }

        if (cSelTrueParticle && std::abs(mcpart.pdgCode()) != kKPlus) { // True kaon selection
          continue;
        }
      }

      // K+ / K-
      if (track.sign() >= 0) {
        getDcaHist<kKaonPlus, dmc>(track, track.dcaXY());
        fillKaonQA<kKaonPlus>(track);
        partType = kKaonPlus;
      } else if (track.sign() <= 0) {
        getDcaHist<kKaonMinus, dmc>(track, track.dcaXY());
        fillKaonQA<kKaonMinus>(track);
        partType = kKaonMinus;
      } else {
        continue;
      }

      // Get Kaon correction factor
      if (cGetCorrectionFlag) {
        kaonCorrFact = (partType == kKaonPlus) ? getCorrectionFactors<kKaonPlus>(track, rap) : getCorrectionFactors<kKaonMinus>(track, rap);
      }

      // Fill table
      kaonTrackTable(lambdaCollisionTable.lastIndex(), track.pt(), track.eta(), track.phi(), rap, track.px(), track.py(), track.pz(), MassKaonCharged,
                     track.globalIndex(), (int8_t)partType, kaonCorrFact);
    }
  }

  // MC Generater Level Tables
  template <typename C, typename M>
  void fillLambdaMcGenTables(C const& mcCollision, M const& mcParticles)
  {
    // Fill McGen Collision Table
    lambdaMCGenCollisionTable(cent, mult, mcCollision.posX(), mcCollision.posY(), mcCollision.posZ());

    // initialize track objects
    ParticleType partType = kLambda;

    // Loop over MC particles
    for (auto const& mcpart : mcParticles) {
      // Check for Primary Lambda/Anti-Lambda/K+/K-
      if (mcpart.isPhysicalPrimary() && mcpart.pdgCode() == kLambda0) {
        partType = kLambda;
      } else if (mcpart.isPhysicalPrimary() && mcpart.pdgCode() == kLambda0Bar) {
        partType = kAntiLambda;
      } else if (mcpart.isPhysicalPrimary() && mcpart.pdgCode() == kKPlus) {
        partType = kKaonPlus;
      } else if (mcpart.isPhysicalPrimary() && mcpart.pdgCode() == kKMinus) {
        partType = kKaonMinus;
      } else {
        continue;
      }

      // Fill Lambda Table
      if (partType == kLambda || partType == kAntiLambda) {
        // Kinematic selection
        if (mcpart.pt() <= cLambdaMinPt || mcpart.pt() >= cLambdaMaxPt || std::abs(mcpart.y()) >= cLambdaRapCut) {
          continue;
        }

        histos.fill(HIST("Tracks/h1f_lambda_info"), kGenTotAccLambda);

        // get daughter track info and check for decay channel flag
        if (!mcpart.has_daughters()) {
          histos.fill(HIST("Tracks/h1f_lambda_info"), kGenLambdaNoDau);
          continue;
        }
        auto dautracks = mcpart.template daughters_as<aod::McParticles>();
        std::vector<int> daughterPDGs, daughterIDs;
        for (auto const& dautrack : dautracks) {
          daughterPDGs.push_back(dautrack.pdgCode());
          daughterIDs.push_back(dautrack.globalIndex());
        }

        if (partType == kLambda) {
          histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[0]);
          histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[1]);
          histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), mcpart.pdgCode());
        } else {
          histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[0]);
          histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[1]);
          histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), mcpart.pdgCode());
        }
        // Fill table
        lambdaMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), mcpart.px(), mcpart.py(), mcpart.pz(), RecoDecay::m(mcpart.p(), mcpart.e()),
                              daughterIDs[0], daughterIDs[1], (int8_t)partType, 1.);
      }

      // Fill Kaon Table
      if (partType == kKaonPlus || partType == kKaonMinus) {
        // Kinematic selection
        if (mcpart.pt() <= cKaonMinPt || mcpart.pt() >= cKaonMaxPt || std::abs(mcpart.y()) >= cKaonRapCut) {
          continue;
        }

        // Fill table
        kaonMCGenTrackTable(lambdaMCGenCollisionTable.lastIndex(), mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(), mcpart.px(), mcpart.py(), mcpart.pz(), RecoDecay::m(mcpart.p(), mcpart.e()),
                            mcpart.globalIndex(), (int8_t)partType, 1.);
      }
    }
  }

  template <DMCType dmc, typename M, typename C, typename B, typename V, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, B const& bc, V const& V0s, T const& tracks, P const& mcParticles)
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
    if (!collisions.begin().has_mcCollision() || !selCollision(collisions.begin()) || collisions.begin().mcCollisionId() != mcCollision.globalIndex()) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(), collisions.begin().posZ());
    auto v0Tracks = V0s.sliceBy(v0sPerCollision, collisions.begin().globalIndex());
    auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collisions.begin().globalIndex());
    fillLambdaRecoTables<dmc>(collisions.begin(), bc, v0Tracks, tracksThisCollision);
    fillLambdaMcGenTables(mcCollision, mcParticles);
  }

  // BC, Collision, tracks and V0s
  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using Collisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::TPCMults, aod::PVMults, aod::MultsGlobal>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksExtra, aod::TracksDCA, aod::pidTPCEl, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::TOFSignal, aod::pidTOFbeta, aod::TrackCompColls>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;

  SliceCache cache;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> v0sPerCollision = aod::v0data::collisionId;
  Preslice<TracksMC> tracksPerCollision = aod::track::collisionId;

  void processDummy(Collisions::iterator const&) {}

  PROCESS_SWITCH(LambdaTableProducer, processDummy, "Dummy Process", true);

  void processData(Collisions::iterator const& collision, BCsRun3 const& bc, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    fillLambdaRecoTables<kData>(collision, bc, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processData, "Process for DATA", false);

  void processMCReco(soa::Join<Collisions, aod::McCollisionLabels>::iterator const& collision, aod::McCollisions const&, BCsRun3 const& bc,
                     McV0Tracks const& V0s, TracksMC const& tracks, aod::McParticles const&)
  {
    // check collision
    if (!selCollision(collision)) {
      return;
    }
    fillLambdaRecoTables<kMC>(collision, bc, V0s, tracks);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCReco, "Process for McReco DATA", false);

  void processMCRecoGen(aod::McCollisions::iterator const& mcCollision,
                        soa::SmallGroups<soa::Join<Collisions, aod::McCollisionLabels>> const& collisions, BCsRun3 const& bc,
                        McV0Tracks const& V0s, TracksMC const& tracks,
                        aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kMC>(mcCollision, collisions, bc, V0s, tracks, mcParticles);
  }

  PROCESS_SWITCH(LambdaTableProducer, processMCRecoGen, "Process for MC RecoGen", false);
};

struct LambdaTracksExtProducer {
  // Tables
  Produces<aod::LambdaTracksExt> lambdaTrackExtTable;
  Produces<aod::KaonTracksExt> kaonTrackExtTable;

  // Configurables
  Configurable<bool> cAcceptAllLambda{"cAcceptAllLambda", false, "Accept all Lambda"};
  Configurable<bool> cRejAllLambdaShaDau{"cRejAllLambdaShaDau", true, "Reject all Lambda sharing daughters"};
  Configurable<bool> cAcceptAllKaon{"cAcceptAllKaon", false, "Accept all Kaons"};
  Configurable<bool> cRejAllKaonShaLaDau{"cRejAllKaonShaLaDau", true, "Reject all Kaons sharing Lambda daughters"};

  // Histogram Registry.
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    // Axis Specifications
    const AxisSpec axisMult(10, 0, 10);
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    // Histograms Booking
    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_lambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_antilambda_mult", "Multiplicity", kTH1I, {axisMult});

    histos.add("h1i_totkaplus_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totkaminus_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_kaplus_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_kaminus_mult", "Multiplicity", kTH1I, {axisMult});

    histos.add("h2d_n2_etaphi_LaP_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaP_LaP", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_LaM_LaM", "#rho_{2}^{SharePair}", kTH2D, {axisDEta, axisDPhi});

    histos.add("h2d_n2_etaphi_KaPLaP", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_KaPLaM", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_KaMLaP", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("h2d_n2_etaphi_KaMLaM", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
  }

  void processDummy(aod::LambdaCollisions::iterator const&) {}

  PROCESS_SWITCH(LambdaTracksExtProducer, processDummy, "Dummy Process", true);

  void processLambdaTrackExt(aod::LambdaCollisions::iterator const&, aod::LambdaTracks const& tracks)
  {
    int nTotLambda = 0, nTotAntiLambda = 0, nSelLambda = 0, nSelAntiLambda = 0;

    for (auto const& lambda : tracks) {
      bool lambdaSharingDauFlag = false, trueLambdaFlag = false;
      std::vector<int64_t> vSharedDauLambdaIndex;

      if (lambda.partType() == kLambda) {
        ++nTotLambda;
      } else if (lambda.partType() == kAntiLambda) {
        ++nTotAntiLambda;
      }

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
          if ((lambda.partType() == kLambda && track.partType() == kAntiLambda) || (lambda.partType() == kAntiLambda && track.partType() == kLambda)) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.partType() == kLambda && track.partType() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaP"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.partType() == kAntiLambda && track.partType() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          }
        }
      }

      // Accept / Reject
      trueLambdaFlag = cAcceptAllLambda || (cRejAllLambdaShaDau && !lambdaSharingDauFlag);

      // Multiplicity of selected lambda
      if (trueLambdaFlag) {
        if (lambda.partType() == kLambda) {
          ++nSelLambda;
        } else if (lambda.partType() == kAntiLambda) {
          ++nSelAntiLambda;
        }
      }

      // Fill LambdaTrackExt table
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

  PROCESS_SWITCH(LambdaTracksExtProducer, processLambdaTrackExt, "Process for lambda track extension", false);

  void processKaonTrackExt(aod::LambdaCollisions::iterator const&, aod::LambdaTracks const& lambdaTracks, aod::KaonTracks const& kaonTracks)
  {
    int nTotKaonPlus = 0, nTotKaonMinus = 0, nSelKaonPlus = 0, nSelKaonMinus = 0;

    for (auto const& kaonTrack : kaonTracks) {
      bool kaonSharingLambdaDauFlag = false, trueKaonFlag = false;
      std::vector<int64_t> vKaonShareDauLambdaIndex;

      if (kaonTrack.partType() == kKaonPlus) {
        ++nTotKaonPlus;
      } else if (kaonTrack.partType() == kKaonMinus) {
        ++nTotKaonMinus;
      }

      for (auto const& lambdaTrack : lambdaTracks) {
        // Removal based on shared track index
        if (kaonTrack.kaonTrackId() == lambdaTrack.posTrackId() || kaonTrack.kaonTrackId() == lambdaTrack.negTrackId()) {
          vKaonShareDauLambdaIndex.push_back(kaonTrack.kaonTrackId());
          kaonSharingLambdaDauFlag = true;

          // Fill Deta-Dphi Histogram
          if (kaonTrack.partType() == kKaonPlus && lambdaTrack.partType() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_KaPLaP"), kaonTrack.eta() - lambdaTrack.eta(), RecoDecay::constrainAngle(kaonTrack.phi() - lambdaTrack.phi(), -PIHalf));
          } else if (kaonTrack.partType() == kKaonPlus && lambdaTrack.partType() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_KaPLaM"), kaonTrack.eta() - lambdaTrack.eta(), RecoDecay::constrainAngle(kaonTrack.phi() - lambdaTrack.phi(), -PIHalf));
          } else if (kaonTrack.partType() == kKaonMinus && lambdaTrack.partType() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_KaMLaP"), kaonTrack.eta() - lambdaTrack.eta(), RecoDecay::constrainAngle(kaonTrack.phi() - lambdaTrack.phi(), -PIHalf));
          } else if (kaonTrack.partType() == kKaonMinus && lambdaTrack.partType() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_KaMLaM"), kaonTrack.eta() - lambdaTrack.eta(), RecoDecay::constrainAngle(kaonTrack.phi() - lambdaTrack.phi(), -PIHalf));
          }
        }
      }

      // Accept / Reject
      trueKaonFlag = cAcceptAllKaon || (cRejAllKaonShaLaDau && !kaonSharingLambdaDauFlag);

      // Multiplicity of selected kaons
      if (trueKaonFlag) {
        if (kaonTrack.partType() == kKaonPlus) {
          ++nSelKaonPlus;
        } else if (kaonTrack.partType() == kKaonMinus) {
          ++nSelKaonMinus;
        }
      }

      // Fill LambdaTrackExt table
      kaonTrackExtTable(kaonSharingLambdaDauFlag, vKaonShareDauLambdaIndex, trueKaonFlag);
    }

    // Fill multiplicity histograms
    if (nTotKaonPlus != 0) {
      histos.fill(HIST("h1i_totkaplus_mult"), nTotKaonPlus);
    }

    if (nTotKaonMinus != 0) {
      histos.fill(HIST("h1i_totkaminus_mult"), nTotKaonMinus);
    }

    if (nSelKaonPlus != 0) {
      histos.fill(HIST("h1i_kaplus_mult"), nSelKaonPlus);
    }

    if (nSelKaonMinus != 0) {
      histos.fill(HIST("h1i_kaminus_mult"), nSelKaonMinus);
    }
  }

  PROCESS_SWITCH(LambdaTracksExtProducer, processKaonTrackExt, "Process for kaon track extension", false);
};

struct LambdaR2Correlation {
  // Global Configurables
  Configurable<int> cLambdaNPtBins{"cLambdaNPtBins", 34, "N pT Bins"};
  Configurable<float> cLambdaPtMin{"cLambdaPtMin", 0.7, "Lambda pT Min"};
  Configurable<float> cLambdaPtMax{"cLambdaPtMax", 3.4, "Lambda pT Max"};
  Configurable<int> cKaonNPtBins{"cKaonNPtBins", 20, "N pT Bins"};
  Configurable<float> cKaonPtMin{"cKaonPtMin", 0.3, "Kaon pT Min"};
  Configurable<float> cKaonPtMax{"cKaonPtMax", 2.2, "Kaon pT Max"};

  Configurable<int> cNRapBins{"cNRapBins", 10, "N Rapidity Bins"};
  Configurable<float> cMinRap{"cMinRap", -0.5, "Minimum Rapidity"};
  Configurable<float> cMaxRap{"cMaxRap", 0.5, "Maximum Rapidity"};
  Configurable<int> cNPhiBins{"cNPhiBins", 36, "N Phi Bins"};
  Configurable<bool> cAnaPairs{"cAnaPairs", false, "Analyze Pairs Flag"};

  // Lambda Kaon femtoscopic correction
  Configurable<bool> cApplyFemtoSel{"cApplyFemtoSel", false, "Femto qinv selection"};
  Configurable<float> cFemtoCut{"cFemtoCut", 0.02, "Kaon--Lambda Femto qinv cut"};

  // Lambda Kaon two-track cuts
  Configurable<bool> cApplyTwoTrackCut{"cApplyTwoTrackCut", false, "Flag for two track cut"};
  Configurable<float> cDEtaCut{"cDEtaCut", 0.02, "DEta cut"};
  Configurable<float> cDPhiStarCut{"cDPhiStarCut", 0.02, "DPhiStar cut"};

  // Centrality Axis
  ConfigurableAxis cCentBins{"cCentBins", {VARIABLE_WIDTH, 0.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.f, 60.0f, 70.0f, 80.0f, 90.0f, 100.f}, "Variable Mult-Bins"};

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
  float magField = 0.;
  std::array<float, 9> tpcRadii = {0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4};

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
    const AxisSpec axisCent(cCentBins, "FT0C (%)");
    const AxisSpec axisChMult(200, 0, 200, "N_{ch}");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");
    const AxisSpec axisMass(100, 1.06, 1.16, "M_{#Lambda} (GeV/#it{c}^{2})");
    const AxisSpec axisPtLambda(cLambdaNPtBins, cLambdaPtMin, cLambdaPtMax, "p_{T} (GeV/#it{c})");
    const AxisSpec axisPtKaon(cKaonNPtBins, cKaonPtMin, cKaonPtMax, "p_{T} (GeV/#it{c})");
    const AxisSpec axisEta(cNRapBins, cMinRap, cMaxRap, "#eta");
    const AxisSpec axisRap(cNRapBins, cMinRap, cMaxRap, "y");
    const AxisSpec axisPhi(cNPhiBins, 0., TwoPI, "#varphi (rad)");
    const AxisSpec axisRapPhi(knrapphibins, kminrapphi, kmaxrapphi, "y #varphi");

    // Create Histograms.
    // Event
    histos.add("Event/Reco/h1f_collision_posz", "V_{Z} Distribution", kTH1F, {axisPosZ});
    histos.add("Event/Reco/h1f_ft0m_mult_percentile", "FT0M (%)", kTH1F, {axisCent});
    histos.add("Event/Reco/h2f_Mult_vs_Centrality", "N_{ch} vs FT0M(%)", kTProfile, {axisCent});

    // Two track cut
    histos.add("QA/TwoTrackCut/Before/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/TwoTrackCut/After/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/FemtoCut/Before/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});
    histos.add("QA/FemtoCut/After/h2d_n2_detadphi", "#rho_{2}", kTH2D, {axisDEta, axisDPhi});

    // Efficiency Histograms
    // Single Particle Efficiencies
    histos.add("Reco/Efficiency/h2f_n1_centpt_LaP", "#rho_{1}^{#Lambda}", kTH2F, {axisCent, axisPtLambda});
    histos.add("Reco/Efficiency/h2f_n1_centpt_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH2F, {axisCent, axisPtLambda});
    histos.add("Reco/Efficiency/h2f_n1_centpt_KaP", "#rho_{1}^{K^{#plus}}", kTH2F, {axisCent, axisPtKaon});
    histos.add("Reco/Efficiency/h2f_n1_centpt_KaM", "#rho_{1}^{K^{#minus}}", kTH2F, {axisCent, axisPtKaon});
    histos.add("Reco/Efficiency/h3f_n1_centptrap_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPtLambda, axisRap});
    histos.add("Reco/Efficiency/h3f_n1_centptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPtLambda, axisRap});
    histos.add("Reco/Efficiency/h3f_n1_centptrap_KaP", "#rho_{1}^{K^{#plus}}", kTH3F, {axisCent, axisPtKaon, axisRap});
    histos.add("Reco/Efficiency/h3f_n1_centptrap_KaM", "#rho_{1}^{K^{#minus}}", kTH3F, {axisCent, axisPtKaon, axisRap});

    // Single and Two Particle Densities
    // 1D Histograms
    histos.add("Reco/h3f_n1_centptmass_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPtLambda, axisMass});
    histos.add("Reco/h3f_n1_centptmass_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPtLambda, axisMass});
    histos.add("Reco/h3f_n1_centptrap_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisPtLambda, axisRap});
    histos.add("Reco/h3f_n1_centptrap_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisPtLambda, axisRap});
    histos.add("Reco/h3f_n1_centptrap_KaP", "#rho_{1}^{K^{#plus}}", kTH3F, {axisCent, axisPtKaon, axisRap});
    histos.add("Reco/h3f_n1_centptrap_KaM", "#rho_{1}^{K^{#minus}}", kTH3F, {axisCent, axisPtKaon, axisRap});

    // rho1 for R2 RapPhi
    histos.add("Reco/h3f_n1_rapphi_LaP", "#rho_{1}^{#Lambda}", kTH3F, {axisCent, axisRap, axisPhi});
    histos.add("Reco/h3f_n1_rapphi_LaM", "#rho_{1}^{#bar{#Lambda}}", kTH3F, {axisCent, axisRap, axisPhi});
    histos.add("Reco/h3f_n1_rapphi_KaP", "#rho_{1}^{K^{#plus}}", kTH3F, {axisCent, axisRap, axisPhi});
    histos.add("Reco/h3f_n1_rapphi_KaM", "#rho_{1}^{K^{#minus}}", kTH3F, {axisCent, axisRap, axisPhi});

    if (cAnaPairs) {
      // rho2 for R2 Rap1Phi1Rap2Phi2
      histos.add("Reco/h3f_n2_rapphi_LaP_LaM", "#rho_{2}^{#Lambda#bar{#Lambda}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaP_LaP", "#rho_{2}^{#Lambda#Lambda}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaM_LaM", "#rho_{2}^{#bar{#Lambda}#bar{#Lambda}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaP_KaM", "#rho_{2}^{#LambdaK^{#minus}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaP_KaP", "#rho_{2}^{#LambdaK^{#plus}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaM_KaM", "#rho_{2}^{#bar{#Lambda}K^{#plus}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
      histos.add("Reco/h3f_n2_rapphi_LaM_KaP", "#rho_{2}^{#bar{#Lambda}K^{#minus}}", kTH3F, {axisCent, axisRapPhi, axisRapPhi});
    }

    // MCGen
    if (doprocessMCGen) {
      histos.addClone("Event/Reco/", "Event/McGen/");
      histos.addClone("Reco/", "McGen/");
    }
  }

  template <typename A>
  bool checkClosePair(A const& v1, A const& v2, int charge1, int charge2)
  {
    // DPhiStar average over different TPC radius
    float dphistar = 0., n = 0.;
    for (auto const& radius : tpcRadii) {
      float arg1 = 0.15 * magField * radius / v1[0];
      float arg2 = 0.15 * magField * radius / v2[0];
      if (std::abs(arg1) < 1.0 && std::abs(arg2) < 1.0) {
        dphistar += v1[2] - v2[2] - (charge1 * std::abs(std::asin(arg1))) + (charge2 * std::abs(std::asin(arg2)));
        ++n;
      }
    }

    // Inf check
    if (n == 0) {
      return false;
    }

    // DPhiStar
    dphistar = RecoDecay::constrainAngle(dphistar / n, -PIHalf);

    // DEta
    float deta = v1[1] - v2[1];

    // Return
    return (std::abs(deta) < cDEtaCut && std::abs(dphistar) < cDPhiStarCut);
  }

  template <typename T, typename V>
  bool isClosePair(V const& lambda, T const& track)
  {
    // Before
    histos.fill(HIST("QA/TwoTrackCut/Before/h2d_n2_detadphi"), track.eta() - lambda.eta(), RecoDecay::constrainAngle(track.phi() - lambda.phi(), -PIHalf));

    // Close pair flag
    bool retFlag = false;

    // Assign kinematics
    std::array<float, 3> trackKin = {track.pt(), track.eta(), track.phi()};
    std::array<float, 3> lambdaPosTrackKin = {lambda.posTrackKin()[0], lambda.posTrackKin()[1], lambda.posTrackKin()[2]};
    std::array<float, 3> lambdaNegTrackKin = {lambda.negTrackKin()[0], lambda.negTrackKin()[1], lambda.negTrackKin()[2]};

    if (track.partType() == kKaonPlus) {
      retFlag = checkClosePair(trackKin, lambdaPosTrackKin, 1, 1) || checkClosePair(trackKin, lambdaNegTrackKin, 1, -1);
    } else if (track.partType() == kKaonMinus) {
      retFlag = checkClosePair(trackKin, lambdaPosTrackKin, -1, 1) || checkClosePair(trackKin, lambdaNegTrackKin, -1, -1);
    }

    // Fill QA
    if (!retFlag) { // Pair accept
      histos.fill(HIST("QA/TwoTrackCut/After/h2d_n2_detadphi"), track.eta() - lambda.eta(), RecoDecay::constrainAngle(track.phi() - lambda.phi(), -PIHalf));
    }

    // Pair reject
    return retFlag;
  }

  template <typename T, typename V>
  bool isCloseQinv(V const& p1, T const& p2)
  {
    // Before
    histos.fill(HIST("QA/FemtoCut/Before/h2d_n2_detadphi"), p1.eta() - p2.eta(), RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));

    // Calculate qinv
    q = RecoDecay::p((p1.px() - p2.px()), (p1.py() - p2.py()), (p1.pz() - p2.pz()));
    e = RecoDecay::e(p1.px(), p1.py(), p1.pz(), MassLambda0) - RecoDecay::e(p2.px(), p2.py(), p2.pz(), MassKaonCharged);
    qinv = std::sqrt(-RecoDecay::m2(q, e));

    // Pair accept
    if (qinv > cFemtoCut) {
      histos.fill(HIST("QA/FemtoCut/After/h2d_n2_detadphi"), p1.eta() - p2.eta(), RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PIHalf));
      return false;
    }

    // Pair reject
    return true;
  }

  template <ParticlePairType part_pair, RecGenType rec_gen, typename T1, typename T2>
  void fillPairHistos(T1& p1, T2& p2)
  {
    static constexpr auto SubDirRecGen = std::array{"Reco/", "McGen/"};
    static constexpr auto SubDirHist = std::array{"LaP_LaM", "LaP_LaP", "LaM_LaM", "LaP_KaP", "LaP_KaM", "LaM_KaP", "LaM_KaM", "KaP_KaM", "KaP_KaP", "KaM_KaM"};

    const auto rapbin1 = static_cast<int>((p1.rap() - kminrap) / rapbinwidth);
    const auto rapbin2 = static_cast<int>((p2.rap() - kminrap) / rapbinwidth);

    const auto phibin1 = static_cast<int>(p1.phi() / phibinwidth);
    const auto phibin2 = static_cast<int>(p2.phi() / phibinwidth);

    const float corfac = p1.corrFact() * p2.corrFact();

    if (rapbin1 >= 0 && rapbin2 >= 0 && phibin1 >= 0 && phibin2 >= 0 && rapbin1 < nrapbins && rapbin2 < nrapbins && phibin1 < nphibins && phibin2 < nphibins) {

      const int rapphix = rapbin1 * nphibins + phibin1;
      const int rapphiy = rapbin2 * nphibins + phibin2;

      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n2_rapphi_") + HIST(SubDirHist[part_pair]), cent, rapphix + 0.5, rapphiy + 0.5, corfac);
    }
  }

  template <ParticleType part, RecGenType rec_gen, typename T>
  void analyzeSingles(T const& tracks)
  {
    static constexpr auto SubDirRecGen = std::array{"Reco/", "McGen/"};
    static constexpr auto SubDirHist = std::array{"LaP", "LaM", "KaP", "KaM"};

    for (auto const& track : tracks) {
      // Efficiency Plots
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h2f_n1_centpt_") + HIST(SubDirHist[part]), cent, track.pt());
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("Efficiency/h3f_n1_centptrap_") + HIST(SubDirHist[part]), cent, track.pt(), track.rap());

      // QA Plots
      if (part == kLambda || part == kAntiLambda) {
        histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n1_centptmass_") + HIST(SubDirHist[part]), cent, track.pt(), track.mass());
      }
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n1_centptrap_") + HIST(SubDirHist[part]), cent, track.pt(), track.rap(), track.corrFact());

      // Rho1 for N1RapPhi
      histos.fill(HIST(SubDirRecGen[rec_gen]) + HIST("h3f_n1_rapphi_") + HIST(SubDirHist[part]), cent, track.rap(), track.phi(), track.corrFact());
    }
  }

  template <ParticlePairType partpair, RecGenType rec_gen, bool same, typename T1, typename T2>
  void analyzePairs(T1 const& trks_1, T2 const& trks_2)
  {
    for (auto const& trk_1 : trks_1) {
      for (auto const& trk_2 : trks_2) {
        // check for same index for Lambda-Lambda / AntiLambda-AntiLambda
        if (same && ((trk_1.index() == trk_2.index()))) {
          continue;
        }

        // Lambda-Kaon close pair rejection
        if constexpr (rec_gen == kRec) {
          if constexpr (partpair == kLambdaKaonPlus || partpair == kLambdaKaonMinus || partpair == kAntiLambdaKaonPlus || partpair == kAntiLambdaKaonMinus) {
            // Close pair
            if (cApplyTwoTrackCut && isClosePair(trk_1, trk_2)) {
              continue;
            }
            // Femto selection
            if (cApplyFemtoSel && isCloseQinv(trk_1, trk_2)) {
              continue;
            }
          }
        }

        // Fill pair histograms
        fillPairHistos<partpair, rec_gen>(trk_1, trk_2);
      }
    }
  }

  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaTracks = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;
  using KaonTracks = soa::Join<aod::KaonTracks, aod::KaonTracksExt>;

  SliceCache cache;
  Partition<LambdaTracks> partLambdaTracks = (aod::lambdatrack::partType == (int8_t)kLambda) && (aod::lambdatrackext::trueLambdaFlag == true);
  Partition<LambdaTracks> partAntiLambdaTracks = (aod::lambdatrack::partType == (int8_t)kAntiLambda) && (aod::lambdatrackext::trueLambdaFlag == true);
  Partition<KaonTracks> partKaonPlusTracks = (aod::kaontrack::partType == (int8_t)kKaonPlus) && (aod::kaontrackext::trueKaonFlag == true);
  Partition<KaonTracks> partKaonMinusTracks = (aod::kaontrack::partType == (int8_t)kKaonMinus) && (aod::kaontrackext::trueKaonFlag == true);

  void processDummy(aod::LambdaCollisions::iterator const&) {}

  PROCESS_SWITCH(LambdaR2Correlation, processDummy, "Dummy Process", true);

  void processDataReco(LambdaCollisions::iterator const& collision, LambdaTracks const&, KaonTracks const&)
  {
    histos.fill(HIST("Event/Reco/h1f_collision_posz"), collision.posZ());
    histos.fill(HIST("Event/Reco/h1f_ft0m_mult_percentile"), collision.cent());
    histos.fill(HIST("Event/Reco/h2f_Mult_vs_Centrality"), collision.cent(), collision.mult());

    cent = collision.cent();
    magField = collision.bField();

    auto lambdaTracks = partLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto antiLambdaTracks = partAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto kaonPlusTracks = partKaonPlusTracks->sliceByCached(aod::kaontrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto kaonMinusTracks = partKaonMinusTracks->sliceByCached(aod::kaontrack::lambdaCollisionId, collision.globalIndex(), cache);

    analyzeSingles<kLambda, kRec>(lambdaTracks);
    analyzeSingles<kAntiLambda, kRec>(antiLambdaTracks);
    analyzeSingles<kKaonPlus, kRec>(kaonPlusTracks);
    analyzeSingles<kKaonMinus, kRec>(kaonMinusTracks);

    if (cAnaPairs) {
      // Pairs Only
      analyzePairs<kLambdaAntiLambda, kRec, false>(lambdaTracks, antiLambdaTracks);
      analyzePairs<kLambdaLambda, kRec, true>(lambdaTracks, lambdaTracks);
      analyzePairs<kAntiLambdaAntiLambda, kRec, true>(antiLambdaTracks, antiLambdaTracks);
      analyzePairs<kLambdaKaonPlus, kRec, false>(lambdaTracks, kaonPlusTracks);
      analyzePairs<kLambdaKaonMinus, kRec, false>(lambdaTracks, kaonMinusTracks);
      analyzePairs<kAntiLambdaKaonPlus, kRec, false>(antiLambdaTracks, kaonPlusTracks);
      analyzePairs<kAntiLambdaKaonMinus, kRec, false>(antiLambdaTracks, kaonMinusTracks);
    }
  }

  PROCESS_SWITCH(LambdaR2Correlation, processDataReco, "Process for Data and MCReco", false);

  using LambdaMcGenCollisions = aod::LambdaMcGenCollisions;
  using LambdaMcGenTracks = aod::LambdaMcGenTracks;
  using KaonMcGenTracks = aod::KaonMcGenTracks;

  SliceCache cachemc;
  Partition<LambdaMcGenTracks> partMcLambdaTracks = (aod::lambdatrack::partType == (int8_t)kLambda);
  Partition<LambdaMcGenTracks> partMcAntiLambdaTracks = (aod::lambdatrack::partType == (int8_t)kAntiLambda);
  Partition<KaonMcGenTracks> partMcKaonPlusTracks = (aod::kaontrack::partType == (int8_t)kKaonPlus);
  Partition<KaonMcGenTracks> partMcKaonMinusTracks = (aod::kaontrack::partType == (int8_t)kKaonMinus);

  void processMCGen(LambdaMcGenCollisions::iterator const& mcgencol, LambdaMcGenTracks const&, KaonMcGenTracks const&)
  {
    histos.fill(HIST("Event/McGen/h1f_collision_posz"), mcgencol.posZ());
    histos.fill(HIST("Event/McGen/h1f_ft0m_mult_percentile"), mcgencol.cent());
    histos.fill(HIST("Event/McGen/h2f_Mult_vs_Centrality"), mcgencol.cent(), mcgencol.mult());

    cent = mcgencol.cent();

    auto lambdaTracks = partMcLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cache);
    auto antiLambdaTracks = partMcAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cache);
    auto kaonPlusTracks = partMcKaonPlusTracks->sliceByCached(aod::kaonmcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cache);
    auto kaonMinusTracks = partMcKaonMinusTracks->sliceByCached(aod::kaonmcgentrack::lambdaMcGenCollisionId, mcgencol.globalIndex(), cache);

    analyzeSingles<kLambda, kGen>(lambdaTracks);
    analyzeSingles<kAntiLambda, kGen>(antiLambdaTracks);
    analyzeSingles<kKaonPlus, kGen>(kaonPlusTracks);
    analyzeSingles<kKaonMinus, kGen>(kaonMinusTracks);

    if (cAnaPairs) {
      analyzePairs<kLambdaAntiLambda, kGen, false>(lambdaTracks, antiLambdaTracks);
      analyzePairs<kLambdaLambda, kGen, true>(lambdaTracks, lambdaTracks);
      analyzePairs<kAntiLambdaAntiLambda, kGen, true>(antiLambdaTracks, antiLambdaTracks);
      analyzePairs<kLambdaKaonPlus, kGen, false>(lambdaTracks, kaonPlusTracks);
      analyzePairs<kLambdaKaonMinus, kGen, false>(lambdaTracks, kaonMinusTracks);
      analyzePairs<kAntiLambdaKaonPlus, kGen, false>(antiLambdaTracks, kaonPlusTracks);
      analyzePairs<kAntiLambdaKaonMinus, kGen, false>(antiLambdaTracks, kaonMinusTracks);
    }
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
