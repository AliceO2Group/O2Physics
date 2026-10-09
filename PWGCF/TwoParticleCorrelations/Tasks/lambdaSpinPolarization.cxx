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

#include <TAxis.h>
#include <TH1.h>
#include <THnBase.h>
#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
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
} // namespace lambdacollision

DECLARE_SOA_TABLE(LambdaCollisions, "AOD", "LAMBDACOLS", o2::soa::Index<>,
                  lambdacollision::Cent, lambdacollision::Mult,
                  aod::collision::PosX, aod::collision::PosY,
                  aod::collision::PosZ, lambdacollision::TimeStamp);
using LambdaCollision = LambdaCollisions::iterator;

namespace lambdamcgencollision
{
}
DECLARE_SOA_TABLE(LambdaMcGenCollisions, "AOD", "LMCGENCOLS", o2::soa::Index<>,
                  lambdacollision::Cent, lambdacollision::Mult,
                  o2::aod::mccollision::PosX, o2::aod::mccollision::PosY,
                  o2::aod::mccollision::PosZ, lambdacollision::TimeStamp);
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
                  lambdatrack::LambdaCollisionId, lambdatrack::Px,
                  lambdatrack::Py, lambdatrack::Pz, lambdatrack::Pt,
                  lambdatrack::Eta, lambdatrack::Phi, lambdatrack::Rap,
                  lambdatrack::Mass, lambdatrack::PrPx, lambdatrack::PrPy,
                  lambdatrack::PrPz, lambdatrack::PosTrackId,
                  lambdatrack::NegTrackId, lambdatrack::CosPA,
                  lambdatrack::DcaDau, lambdatrack::V0Type,
                  lambdatrack::V0PrmScd, lambdatrack::CorrFact);
using LambdaTrack = LambdaTracks::iterator;

namespace lambdatrackext
{
DECLARE_SOA_COLUMN(LambdaSharingDaughter, lambdaSharingDaughter, bool);
} // namespace lambdatrackext

DECLARE_SOA_TABLE(LambdaTracksExt, "AOD", "LAMBDATRACKSEXT",
                  lambdatrackext::LambdaSharingDaughter);
using LambdaTrackExt = LambdaTracksExt::iterator;

namespace lambdamcgentrack
{
DECLARE_SOA_INDEX_COLUMN(LambdaMcGenCollision, lambdaMcGenCollision);
}
DECLARE_SOA_TABLE(LambdaMcGenTracks, "AOD", "LMCGENTRACKS", o2::soa::Index<>,
                  lambdamcgentrack::LambdaMcGenCollisionId,
                  o2::aod::mcparticle::Px, o2::aod::mcparticle::Py,
                  o2::aod::mcparticle::Pz, lambdatrack::Pt, lambdatrack::Eta,
                  lambdatrack::Phi, lambdatrack::Rap, lambdatrack::Mass,
                  lambdatrack::PrPx, lambdatrack::PrPy, lambdatrack::PrPz,
                  lambdatrack::PosTrackId, lambdatrack::NegTrackId,
                  lambdatrack::V0Type, lambdatrack::CosPA, lambdatrack::DcaDau,
                  lambdatrack::V0PrmScd, lambdatrack::CorrFact);
using LambdaMcGenTrack = LambdaMcGenTracks::iterator;

namespace lambdamixeventcollision
{
DECLARE_SOA_COLUMN(CollisionIndex, collisionIndex, int);
} // namespace lambdamixeventcollision

DECLARE_SOA_TABLE(LambdaMixEventCollisions, "AOD", "LAMBDAMIXCOLS",
                  o2::soa::Index<>, lambdamixeventcollision::CollisionIndex,
                  lambdacollision::Cent, aod::collision::PosZ,
                  lambdacollision::TimeStamp);
using LambdaMixEventCollision = LambdaMixEventCollisions::iterator;

namespace lambdamixeventmcgencollision
{
}
DECLARE_SOA_TABLE(LambdaMixEventMcGenCollisions, "AOD", "LAMBDAMIXMGCOLS",
                  o2::soa::Index<>, lambdamixeventcollision::CollisionIndex,
                  lambdacollision::Cent, aod::collision::PosZ,
                  lambdacollision::TimeStamp);
using LambdaMixEventMcGenCollision = LambdaMixEventMcGenCollisions::iterator;

namespace lambdamixeventtracks
{
DECLARE_SOA_COLUMN(LambdaMixEventCollisionIdx, lambdaMixEventCollisionIdx, int);
DECLARE_SOA_COLUMN(LambdaMixEventTrackIdx, lambdaMixEventTrackIdx, int);
DECLARE_SOA_COLUMN(LambdaMixEventTimeStamp, lambdaMixEventTimeStamp, uint64_t);
DECLARE_SOA_COLUMN(LambdaMixEventPosTrackIdx, lambdaMixEventPosTrackIdx, int);
DECLARE_SOA_COLUMN(LambdaMixEventNegTrackIdx, lambdaMixEventNegTrackIdx, int);
} // namespace lambdamixeventtracks

DECLARE_SOA_TABLE(LambdaMixEventTracks, "AOD", "LAMBDAMIXTRKS",
                  o2::soa::Index<>,
                  lambdamixeventtracks::LambdaMixEventCollisionIdx,
                  lambdamixeventtracks::LambdaMixEventTrackIdx, lambdatrack::Px,
                  lambdatrack::Py, lambdatrack::Pz, lambdatrack::Mass,
                  lambdatrack::PrPx, lambdatrack::PrPy, lambdatrack::PrPz,
                  lambdatrack::V0Type,
                  lambdamixeventtracks::LambdaMixEventTimeStamp,
                  lambdamixeventtracks::LambdaMixEventPosTrackIdx,
                  lambdamixeventtracks::LambdaMixEventNegTrackIdx);
using LambdaMixEventTrack = LambdaMixEventTracks::iterator;

namespace lambdamixeventmcgentracks
{
}
DECLARE_SOA_TABLE(LambdaMixEventMcGenTracks, "AOD", "LAMBDAMIXMGTRKS",
                  o2::soa::Index<>,
                  lambdamixeventtracks::LambdaMixEventCollisionIdx,
                  lambdamixeventtracks::LambdaMixEventTrackIdx, lambdatrack::Px,
                  lambdatrack::Py, lambdatrack::Pz, lambdatrack::Mass,
                  lambdatrack::PrPx, lambdatrack::PrPy, lambdatrack::PrPz,
                  lambdatrack::V0Type,
                  lambdamixeventtracks::LambdaMixEventTimeStamp);
using LambdaMixEventMcGenTrack = LambdaMixEventMcGenTracks::iterator;
} // namespace o2::aod

enum CollisionLabels { kTotColBeforeHasMcCollision = 1,
                       kTotCol,
                       kPassSelCol };

enum EventCutFlow { kEvAll = 1,
                    kEvTrigger,
                    kEvTvx,
                    kEvTFBorder,
                    kEvItsRofBorder,
                    kEvItsTpcVtx,
                    kEvNoSameBunchPileup,
                    kEvGoodZvtxFT0vsPV,
                    kEvGoodItsLayers,
                    kEvCentrality,
                    kEvVz,
                    kEvOneV0,
                    kEvTwoV0,
                    kNEventCutFlow };

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

enum CentEstType { kCentFT0M = 0,
                   kCentFT0C };

enum RunType { kRun3 = 0,
               kRun2 };

enum ParticleType { kLambda = 0,
                    kAntiLambda };

enum ParticlePairType { kLambdaAntiLambda = 0,
                        kAntiLambdaLambda,
                        kLambdaLambda,
                        kAntiLambdaAntiLambda };

enum PairTier { kSameEvent = 0,
                kSameEventInMixing,
                kMixedEvent,
                kMcGenSameEvent,
                kMcGenSameEventInMixing,
                kMcGenMixedEvent,
                kSameEventSharedDau,
                kSameEventInMixingSharedDau };

enum SharedDauRule { kSharedDauMaxSet = 0,
                     kSharedDauPairVeto,
                     kSharedDauDropAll };

enum PairEventCount { kPeAll = 1,
                      kPeOneCand,
                      kPeTwoCand,
                      kPeUnlikeSign,
                      kPeLambdaLambda,
                      kPeAntiLambdaAntiLambda,
                      kNPairEventCount };

enum MEEventStatus { kMEEventSeen = 1,
                     kMEEventOutsidePools,
                     kMEEventWarmUp,
                     kMEEventMixed,
                     kMEEventDonated,
                     kNMEEventStatus };

enum ShareDauLambda { kUniqueLambda = 0,
                      kLambdaShareDau };

enum RecGenType { kRec = 0,
                  kGen };

enum DMCType { kData = 0,
               kMC };

enum CorrHistDim { OneDimCorr = 1,
                   TwoDimCorr,
                   ThreeDimCorr };

enum PrmScdType { kPrimary = 0,
                  kSecondary };

enum PrmScdPairType { kPP = 0,
                      kPS,
                      kSP,
                      kSS };

static constexpr std::size_t NDaughtersTwoBody = 2;

static constexpr std::size_t NCandidatesForPair = 2;

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

  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.15, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 999.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<int> cMinTpcCrossedRows{"cMinTpcCrossedRows", 70, "TPC Min Crossed Rows"};
  Configurable<float> cMinTpcCROverCls{"cMinTpcCROverCls", -999, "TPC Min CR/Findable Cls"};
  Configurable<float> cMaxTpcSharedClusters{"cMaxTpcSharedClusters", 999, "TPC Max Shared Clusters"};
  Configurable<float> cMaxChi2Tpc{"cMaxChi2Tpc", 999, "Max TPC Chi2/ndf"};
  Configurable<double> cTpcNsigmaCut{"cTpcNsigmaCut", 5.0, "TPC nSigma PID cut"};
  Configurable<bool> cRemoveAmbiguousTracks{"cRemoveAmbiguousTracks", false, "Remove Ambiguous Tracks"};

  Configurable<double> cMinDcaProtonToPV{"cMinDcaProtonToPV", 0.05, "Min proton DCA to PV (cm)"};
  Configurable<double> cMinDcaPionToPV{"cMinDcaPionToPV", 0.1, "Min pion DCA to PV (cm)"};
  Configurable<double> cMinV0DcaDaughters{"cMinV0DcaDaughters", 0., "Min DCA between V0 daughters"};
  Configurable<double> cMaxV0DcaDaughters{"cMaxV0DcaDaughters", 1.4, "Max DCA between V0 daughters"};
  Configurable<double> cMinDcaV0ToPV{"cMinDcaV0ToPV", 0.0, "Min DCA V0 to PV"};
  Configurable<double> cMaxDcaV0ToPV{"cMaxDcaV0ToPV", 0.1, "Max DCA V0 to PV"};
  Configurable<double> cMinV0TransRadius{"cMinV0TransRadius", 0.5, "Min V0 decay radius (cm)"};
  Configurable<double> cMaxV0TransRadius{"cMaxV0TransRadius", 999.0, "Max V0 decay radius (cm)"};
  Configurable<double> cMinV0CTau{"cMinV0CTau", 0.0, "Min cTau (cm)"};
  Configurable<double> cMaxV0CTau{"cMaxV0CTau", 30.0, "Max cTau (cm)"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Min V0 cos(PA)"};
  Configurable<double> cKshortRejMassWindow{"cKshortRejMassWindow", 0.001, "K0s mass rejection window"};
  Configurable<bool> cKshortRejFlag{"cKshortRejFlag", true, "K0s mass rejection flag"};

  Configurable<float> cMinV0Mass{"cMinV0Mass", 1.09, "V0 Mass Min"};
  Configurable<float> cMaxV0Mass{"cMaxV0Mass", 1.14, "V0 Mass Max"};
  Configurable<float> cMinV0Pt{"cMinV0Pt", 0.6, "Minimum V0 pT"};
  Configurable<float> cMaxV0Pt{"cMaxV0Pt", 3.0, "Maximum V0 pT"};
  Configurable<float> cMaxV0Rap{"cMaxV0Rap", 0.5, "|rap| cut"};
  Configurable<bool> cDoEtaAnalysis{"cDoEtaAnalysis", false, "Do Eta Analysis"};
  Configurable<bool> cV0TypeSelFlag{"cV0TypeSelFlag", true, "V0 Type Selection Flag"};
  Configurable<int> cV0TypeSelection{"cV0TypeSelection", 1, "V0 Type Selection"};

  Configurable<bool> cHasMcFlag{"cHasMcFlag", true, "Has Mc Tag"};
  Configurable<bool> cSelectTrueLambda{"cSelectTrueLambda", false, "Select True Lambda"};
  Configurable<bool> cSelMCPSV0{"cSelMCPSV0", false, "Select Primary/Secondary V0"};
  Configurable<bool> cCheckRecoDauFlag{"cCheckRecoDauFlag", true, "Check for reco daughter PID"};
  Configurable<bool> cGenPrimaryLambda{"cGenPrimaryLambda", true, "Primary Generated Lambda"};
  Configurable<bool> cGenSecondaryLambda{"cGenSecondaryLambda", false, "Secondary Generated Lambda"};
  Configurable<bool> cGenDecayChannel{"cGenDecayChannel", false, "Gen Level Decay Channel Flag"};
  Configurable<bool> cRecoMomResoFlag{"cRecoMomResoFlag", false, "Check effect of momentum space smearing on balance function"};

  Configurable<bool> cCorrectionFlag{"cCorrectionFlag", false, "Correction Flag"};
  Configurable<bool> cGetEffFact{"cGetEffFact", true, "Get Efficiency Factor Flag"};
  Configurable<bool> cGetPrimFrac{"cGetPrimFrac", false, "Get Primary Fraction Flag"};
  Configurable<int> cCorrFactHist{"cCorrFactHist", 0, "Efficiency Factor Histogram"};
  Configurable<int> cPrimFracHist{"cPrimFracHist", 0, "Primary Fraction Histogram"};

  Configurable<std::string> cUrlCCDB{"cUrlCCDB", "http://alice-ccdb.cern.ch", "url of ccdb"};
  Configurable<std::string> cPathCCDB{"cPathCCDB", "Users/y/ypatley/SpinCorr/RecoEfficiency", "Path for ccdb-object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<std::vector<std::string>> vCorrFactStrings = {
    {"hEffVsPtCentLambda", "hEffVsPtCentAntiLambda"},
    {"hEffVsPtYCentLambda", "hEffVsPtYCentAntiLambda"},
    {"hEffVsPtEtaCentLambda", "hEffVsPtEtaCentAntiLambda"}};

  std::vector<std::vector<std::string>> vPrimFracStrings = {
    {"hPrimFracVsPtCentLambda", "hPrimFracVsPtCentAntiLambda"},
    {"hPrimFracVsPtYCentLambda", "hPrimFracVsPtYCentAntiLambda"},
    {"hPrimFracVsPtEtaCentLambda", "hPrimFracVsPtEtaCentAntiLambda"}};

  float cent = 0., mult = 0.;

  void init(InitContext const&)
  {
    ccdb->setURL(cUrlCCDB.value);
    ccdb->setCaching(true);

    const AxisSpec axisCols(5, 0.5, 5.5, "");
    const AxisSpec axisTrks(30, 0.5, 30.5, "");
    const AxisSpec axisCent(100, 0, 100, "FT0M (%)");
    const AxisSpec axisMult(10, 0, 10, "N_{#Lambda}");
    const AxisSpec axisVz(220, -11, 11, "V_{z} (cm)");
    const AxisSpec axisPID(8000, -4000, 4000, "PdgCode");

    const AxisSpec axisV0Mass(200, 1.08, 1.18, "M_{p#pi} (GeV/#it{c}^{2})");
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

    histos.add("Events/h1f_collisions_info", "# of Collisions", kTH1F, {axisCols});
    histos.add("Events/h1f_collision_posZ", "V_{z}-distribution", kTH1F, {axisVz});
    histos.add("Events/h1f_collision_cent", "Centrality of the selected collisions", kTH1F, {axisCent});
    auto hCutFlow = histos.add<TH1>("Events/hEventCutFlow", "event cut flow;;collisions", kTH1D, {{kNEventCutFlow - 1, 0.5, static_cast<double>(kNEventCutFlow) - 0.5}});
    hCutFlow->GetXaxis()->SetBinLabel(kEvAll, "all");
    hCutFlow->GetXaxis()->SetBinLabel(kEvTrigger, "trigger (sel8)");
    hCutFlow->GetXaxis()->SetBinLabel(kEvTvx, "TVX");
    hCutFlow->GetXaxis()->SetBinLabel(kEvTFBorder, "no TF border");
    hCutFlow->GetXaxis()->SetBinLabel(kEvItsRofBorder, "no ITS ROF border");
    hCutFlow->GetXaxis()->SetBinLabel(kEvItsTpcVtx, "ITS-TPC vertex");
    hCutFlow->GetXaxis()->SetBinLabel(kEvNoSameBunchPileup, "no same-bunch pileup");
    hCutFlow->GetXaxis()->SetBinLabel(kEvGoodZvtxFT0vsPV, "good z_{vtx} FT0 vs PV");
    hCutFlow->GetXaxis()->SetBinLabel(kEvGoodItsLayers, "good ITS layers");
    hCutFlow->GetXaxis()->SetBinLabel(kEvCentrality, "centrality range");
    hCutFlow->GetXaxis()->SetBinLabel(kEvVz, "V_{z} range");
    hCutFlow->GetXaxis()->SetBinLabel(kEvOneV0, "#geq 1 selected #Lambda/#bar{#Lambda}");
    hCutFlow->GetXaxis()->SetBinLabel(kEvTwoV0, "#geq 2 selected #Lambda/#bar{#Lambda}");

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

    histos.addClone("QA/Lambda/", "QA/AntiLambda/");
    histos.addClone("McRec/Lambda/", "McRec/AntiLambda/");

    if (doprocessMCRun3) {
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
    histos.fill(HIST("Events/hEventCutFlow"), kEvAll);
    if constexpr (run == kRun3) {
      if (cCentEstimator == kCentFT0M) {
        cent = col.centFT0M();
      } else if (cCentEstimator == kCentFT0C) {
        cent = col.centFT0C();
      }
      if (cSel8Trig && !col.sel8()) {
        return false;
      }
    } else {
      cent = col.centRun2V0M();
      if (cInt7Trig && !col.alias_bit(kINT7)) {
        return false;
      }
      if (cSel7Trig && !col.sel7()) {
        return false;
      }
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvTrigger);

    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvTvx);
    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvTFBorder);
    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvItsRofBorder);
    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvItsTpcVtx);
    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvNoSameBunchPileup);
    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvGoodZvtxFT0vsPV);
    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvGoodItsLayers);

    if (cent <= cMinMult || cent >= cMaxMult) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvCentrality);
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) {
      return false;
    }
    histos.fill(HIST("Events/hEventCutFlow"), kEvVz);

    mult = col.multNTracksPV();
    return true;
  }

  bool kinCutSelection(float const& pt, float const& rap, float const& ptMin, float const& ptMax, float const& rapMax)
  {
    return pt > ptMin && pt < ptMax && rap < rapMax;
  }

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

  template <typename V, typename T>
  bool selDaughterTracks(V const& v0, T const&, ParticleType const& v0Type)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    if (!selTrack(posTrack) || !selTrack(negTrack)) {
      return false;
    }

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
    if (v0.dcaV0daughters() <= cMinV0DcaDaughters || v0.dcaV0daughters() >= cMaxV0DcaDaughters) {
      return false;
    }
    if (v0.dcav0topv() <= cMinDcaV0ToPV || v0.dcav0topv() >= cMaxDcaV0ToPV) {
      return false;
    }
    if (v0.v0radius() <= cMinV0TransRadius || v0.v0radius() >= cMaxV0TransRadius) {
      return false;
    }

    float ctau = v0.distovertotmom(col.posX(), col.posY(), col.posZ()) * MassLambda0;
    if (ctau <= cMinV0CTau || ctau >= cMaxV0CTau) {
      return false;
    }
    if (v0.v0cosPA() <= cMinV0CosPA) {
      return false;
    }
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
    if (cKshortRejFlag && (std::abs(v0.mK0Short() - MassK0Short) <= cKshortRejMassWindow)) {
      histos.fill(HIST("Tracks/h1f_tracks_info"), kV0KShortMassRej);
      return false;
    }

    auto postrack = v0.template posTrack_as<T>();
    auto negtrack = v0.template negTrack_as<T>();

    bool lambdaFlag = false, antiLambdaFlag = false;

    if ((v0.mLambda() > cMinV0Mass && v0.mLambda() < cMaxV0Mass) && selLambdaDauWithTpcPid<kLambda>(postrack, negtrack)) {
      lambdaFlag = true;
      v0type = kLambda;
    }
    if ((v0.mAntiLambda() > cMinV0Mass && v0.mAntiLambda() < cMaxV0Mass) && selLambdaDauWithTpcPid<kAntiLambda>(postrack, negtrack)) {
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
    if (!selLambdaMassWindow(v0, tracks, v0Type)) {
      return false;
    }
    histos.fill(HIST("Tracks/h1f_tracks_info"), kV0IsLambdaOrAntiLambda);

    if (!selDaughterTracks(v0, tracks, v0Type)) {
      return false;
    }
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0DauTrackSel);

    float rap = cDoEtaAnalysis ? std::abs(v0.eta()) : std::abs(v0.yLambda());
    if (!kinCutSelection(v0.pt(), rap, cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
      return false;
    }
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0KinCuts);

    if (!topoCutSelection(col, v0, tracks)) {
      return false;
    }
    histos.fill(HIST("Tracks/h1f_tracks_info"), kPassV0TopoSel);

    return true;
  }

  template <typename V, typename T>
  bool hasAmbiguousDaughters(V const& v0, T const&)
  {
    auto posTrack = v0.template posTrack_as<T>();
    auto negTrack = v0.template negTrack_as<T>();
    auto posCC = posTrack.compatibleCollIds();
    auto negCC = negTrack.compatibleCollIds();
    if (posCC.size() > 1 || negCC.size() > 1) {
      return true;
    }
    if ((posCC.size() != 0 && posCC[0] != posTrack.collisionId()) ||
        (negCC.size() != 0 && negCC[0] != negTrack.collisionId())) {
      return true;
    }
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
    if (std::abs(mcpart.pdgCode()) != kLambda0) {
      return false;
    }

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
        if (mcpostrack.pdgCode() != kProton ||
            mcnegtrack.pdgCode() != kPiMinus) {
          histos.fill(HIST("Tracks/h1f_tracks_info"), kLambdaNotPrPiMinus);
          return false;
        }
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        if (mcpostrack.pdgCode() != kPiPlus ||
            mcnegtrack.pdgCode() != kProtonBar) {
          histos.fill(HIST("Tracks/h1f_tracks_info"),
                      kAntiLambdaNotAntiPrPiPlus);
          return false;
        }
      }
    }
    return true;
  }

  template <ParticleType part, typename V>
  float getCorrectionFactors(V const& v0)
  {
    if (!cCorrectionFlag) {
      return 1.;
    }

    auto ccdbObj = ccdb->getForTimeStamp<TList>(cPathCCDB.value, 1);
    if (!ccdbObj) {
      LOGF(warning, "CCDB OBJECT NOT FOUND");
      return 1.;
    }

    float effCorrFact = 1., primFrac = 1.;
    float rap = (cDoEtaAnalysis) ? v0.eta() : v0.yLambda();

    if (cGetEffFact) {
      auto* objEff = ccdbObj->FindObject(
        Form("%s", vCorrFactStrings[cCorrFactHist][part].c_str()));
      auto* histEff = dynamic_cast<TH1*>(objEff->Clone());
      if (histEff->GetDimension() == TwoDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtCent);
        effCorrFact = histEff->GetBinContent(histEff->FindBin(cent, v0.pt()));
      } else if (histEff->GetDimension() == ThreeDimCorr) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kEffCorrPtRapCent);
        effCorrFact =
          histEff->GetBinContent(histEff->FindBin(cent, v0.pt(), rap));
      } else {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kNoEffCorr);
        LOGF(warning, "CCDB: not a histogram!");
        effCorrFact = 1.;
      }
      delete histEff;
    }
    if (cGetPrimFrac) {
      auto* objPrm = ccdbObj->FindObject(
        Form("%s", vPrimFracStrings[cPrimFracHist][part].c_str()));
      auto* histPrm = dynamic_cast<TH1*>(objPrm->Clone());
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
    histos.fill(HIST("Tracks/h2f_lambda_mothers_pdg"),
                lambdaMothers[0].pdgCode(), v0.pt());
  }

  template <ParticleType part, typename C, typename V, typename T>
  void fillLambdaQAHistos(C const& col, V const& v0, T const&)
  {
    static constexpr std::array<std::string_view, 2> SubDir = {"QA/Lambda/", "QA/AntiLambda/"};
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
    static constexpr std::array<std::string_view, 2> SubDirRG = {"McRec/", "McGen/"};
    static constexpr std::array<std::string_view, 2> SubDirPart = {"Lambda/", "AntiLambda/"};

    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPt"), pt);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hEta"), eta);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hRap"), y);
    histos.fill(HIST(SubDirRG[rg]) + HIST(SubDirPart[part]) + HIST("hPhi"), phi);
  }

  template <RunType run, DMCType dmc, typename C, typename B, typename V,
            typename T>
  void fillLambdaRecoTables(C const& collision, B const& bc, V const& v0tracks, T const& tracks)
  {
    histos.fill(HIST("Events/h1f_collisions_info"), kTotCol);

    if constexpr (dmc == kData) {
      if (!selCollision<run>(collision)) {
        return;
      }
    }

    histos.fill(HIST("Events/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("Events/h1f_collision_posZ"), collision.posZ());
    histos.fill(HIST("Events/h1f_collision_cent"), cent);

    lambdaCollisionTable(cent, mult, collision.posX(), collision.posY(), collision.posZ(), bc.timestamp());

    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float mass = 0., corr_fact = 1.;
    float prPx = 0., prPy = 0., prPz = 0.;
    float pt = 0., eta = 0., rap = 0., phi = 0.;
    std::size_t nWritten = 0;

    for (auto const& v0 : v0tracks) {
      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h1f_tracks_info"), kTracksBeforeHasMcParticle);
        if (!v0.has_mcParticle()) {
          continue;
        }
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllV0Tracks);
      histos.fill(HIST("Tracks/h2f_armpod_before_sel"), v0.alpha(), v0.qtarm());

      if (!selV0Particle(collision, v0, tracks, v0Type)) {
        continue;
      }
      if (cV0TypeSelFlag && v0.v0Type() != cV0TypeSelection) {
        continue;
      }

      histos.fill(HIST("Tracks/h1f_tracks_info"), kAllSelPassed);

      if constexpr (run == kRun3) {
        if (cRemoveAmbiguousTracks && hasAmbiguousDaughters(v0, tracks)) {
          continue;
        }
      }

      mass = (v0Type == kLambda) ? v0.mLambda() : v0.mAntiLambda();
      pt = v0.pt();
      eta = v0.eta();
      rap = v0.yLambda();
      phi = v0.phi();

      if constexpr (dmc == kMC) {
        histos.fill(HIST("Tracks/h2f_tracks_pid_before_sel"), v0.mcParticle().pdgCode(), v0.pt());
        if (cSelMCPSV0) {
          v0PrmScdType = isPrimaryV0(v0);
        }
        if (cSelectTrueLambda && !selTrueMcRecLambda(v0, tracks)) {
          continue;
        }
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
          float y = cDoEtaAnalysis ? eta : rap;
          if (!kinCutSelection(pt, std::abs(y), cMinV0Pt, cMaxV0Pt, cMaxV0Rap)) {
            continue;
          }
        }
      }

      histos.fill(HIST("Tracks/h2f_armpod_after_sel"), v0.alpha(), v0.qtarm());
      corr_fact = (v0Type == kLambda) ? getCorrectionFactors<kLambda>(v0) : getCorrectionFactors<kAntiLambda>(v0);

      if (v0Type == kLambda) {
        prPx = v0.pxpos();
        prPy = v0.pypos();
        prPz = v0.pzpos();
        histos.fill(HIST("Tracks/h1f_lambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      } else {
        prPx = v0.pxneg();
        prPy = v0.pyneg();
        prPz = v0.pzneg();
        histos.fill(HIST("Tracks/h1f_antilambda_pt_vs_invm"), mass, v0.pt());
        fillLambdaQAHistos<kAntiLambda>(collision, v0, tracks);
        fillKinematicHists<kRec, kAntiLambda>(v0.pt(), v0.eta(), v0.yLambda(), v0.phi());
      }

      lambdaTrackTable(lambdaCollisionTable.lastIndex(), v0.px(), v0.py(),
                       v0.pz(), pt, eta, phi, rap, mass, prPx, prPy, prPz,
                       v0.template posTrack_as<T>().index(),
                       v0.template negTrack_as<T>().index(),
                       v0.v0cosPA(),
                       v0.dcaV0daughters(), (int8_t)v0Type,
                       v0PrmScdType,
                       corr_fact);
      ++nWritten;
    }

    if (nWritten >= 1) {
      histos.fill(HIST("Events/hEventCutFlow"), kEvOneV0);
    }
    if (nWritten >= NCandidatesForPair) {
      histos.fill(HIST("Events/hEventCutFlow"), kEvTwoV0);
    }
  }

  template <RunType run, typename B, typename C, typename M>
  void fillLambdaMcGenTables(B const& bc, C const& mcCollision, M const& mcParticles)
  {
    lambdaMCGenCollisionTable(cent, mult, mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
                              bc.timestamp());

    ParticleType v0Type = kLambda;
    PrmScdType v0PrmScdType = kPrimary;
    float rap = 0.;
    float prPx = 0., prPy = 0., prPz = 0.;

    for (auto const& mcpart : mcParticles) {
      if (mcpart.pdgCode() == kLambda0) {
        v0Type = kLambda;
      } else if (mcpart.pdgCode() == kLambda0Bar) {
        v0Type = kAntiLambda;
      } else {
        continue;
      }

      v0PrmScdType = mcpart.isPhysicalPrimary() ? kPrimary : kSecondary;
      rap = cDoEtaAnalysis ? mcpart.eta() : mcpart.y();
      if (!kinCutSelection(mcpart.pt(), std::abs(rap), cMinV0Pt, cMaxV0Pt,
                           cMaxV0Rap)) {
        continue;
      }
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

      int prIdx = -1, piIdx = -1;
      int const prPdg = (v0Type == kLambda) ? kProton : kProtonBar;
      int const piPdg = (v0Type == kLambda) ? kPiMinus : kPiPlus;
      for (std::size_t i = 0; i < daughterPDGs.size(); ++i) {
        if (daughterPDGs[i] == prPdg) {
          prIdx = static_cast<int>(i);
        } else if (daughterPDGs[i] == piPdg) {
          piIdx = static_cast<int>(i);
        }
      }
      if (prIdx < 0 || piIdx < 0) {
        continue;
      }
      if (cGenDecayChannel && daughterPDGs.size() != NDaughtersTwoBody) {
        continue;
      }
      histos.fill(HIST("Tracks/h1f_tracks_info"), kGenLambdaToPrPi);

      prPx = vDauPx[prIdx];
      prPy = vDauPy[prIdx];
      prPz = vDauPz[prIdx];

      int const posDauIdx = (v0Type == kLambda) ? prIdx : piIdx;
      int const negDauIdx = (v0Type == kLambda) ? piIdx : prIdx;

      if (v0Type == kLambda) {
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[prIdx]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), daughterPDGs[piIdx]);
        histos.fill(HIST("McGen/h1f_lambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/Lambda/Proton/hPt"), vDauPt[prIdx]);
        histos.fill(HIST("McGen/Lambda/Proton/hEta"), vDauEta[prIdx]);
        histos.fill(HIST("McGen/Lambda/Proton/hRap"), vDauRap[prIdx]);
        histos.fill(HIST("McGen/Lambda/Proton/hPhi"), vDauPhi[prIdx]);
        histos.fill(HIST("McGen/Lambda/Pion/hPt"), vDauPt[piIdx]);
        histos.fill(HIST("McGen/Lambda/Pion/hEta"), vDauEta[piIdx]);
        histos.fill(HIST("McGen/Lambda/Pion/hRap"), vDauRap[piIdx]);
        histos.fill(HIST("McGen/Lambda/Pion/hPhi"), vDauPhi[piIdx]);
        fillKinematicHists<kGen, kLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      } else {
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[prIdx]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), daughterPDGs[piIdx]);
        histos.fill(HIST("McGen/h1f_antilambda_daughter_PDG"), mcpart.pdgCode());
        histos.fill(HIST("McGen/AntiLambda/Proton/hPt"), vDauPt[prIdx]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hEta"), vDauEta[prIdx]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hRap"), vDauRap[prIdx]);
        histos.fill(HIST("McGen/AntiLambda/Proton/hPhi"), vDauPhi[prIdx]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hPt"), vDauPt[piIdx]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hEta"), vDauEta[piIdx]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hRap"), vDauRap[piIdx]);
        histos.fill(HIST("McGen/AntiLambda/Pion/hPhi"), vDauPhi[piIdx]);
        fillKinematicHists<kGen, kAntiLambda>(mcpart.pt(), mcpart.eta(), mcpart.y(), mcpart.phi());
      }

      lambdaMCGenTrackTable(
        lambdaMCGenCollisionTable.lastIndex(), mcpart.px(), mcpart.py(),
        mcpart.pz(), mcpart.pt(), mcpart.eta(), mcpart.phi(), mcpart.y(),
        RecoDecay::m(mcpart.p(), mcpart.e()), prPx, prPy, prPz,
        daughterIDs[posDauIdx], daughterIDs[negDauIdx], (int8_t)v0Type, -999.,
        -999., v0PrmScdType, 1.);
    }
  }

  template <RunType run, DMCType dmc, typename M, typename C, typename B, typename V, typename T, typename P>
  void analyzeMcRecoGen(M const& mcCollision, C const& collisions, B const&, V const& V0s, T const& tracks, P const& mcParticles)
  {
    int nRecCols = collisions.size();
    if (nRecCols != 0) {
      histos.fill(HIST("McGen/h1f_collision_recgen"), nRecCols);
    }
    if (nRecCols != 1) {
      return;
    }

    histos.fill(HIST("McGen/h1f_collisions_info"), kTotCol);
    if (!collisions.begin().has_mcCollision() ||
        !selCollision<run>(collisions.begin()) ||
        collisions.begin().mcCollisionId() != mcCollision.globalIndex()) {
      return;
    }
    histos.fill(HIST("McGen/h1f_collisions_info"), kPassSelCol);
    histos.fill(HIST("McGen/h2f_collision_posZ"), mcCollision.posZ(),
                collisions.begin().posZ());

    auto const& recCollision = collisions.begin();
    auto bc = recCollision.template bc_as<B>();
    auto v0sThisCollision =
      V0s.sliceBy(perCollision, recCollision.globalIndex());

    fillLambdaRecoTables<run, dmc>(recCollision, bc, v0sThisCollision, tracks);
    fillLambdaMcGenTables<run>(bc, mcCollision, mcParticles);
  }

  SliceCache cache;
  Preslice<soa::Join<aod::V0Datas, aod::McV0Labels>> perCollision = aod::v0data::collisionId;

  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::PVMults>;
  using CollisionsRun2 = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms, aod::PVMults>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                           aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr,
                           aod::TrackCompColls>;
  using TracksRun2 = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra,
                               aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr>;
  using TracksMC = soa::Join<Tracks, aod::McTrackLabels>;
  using TracksMCRun2 = soa::Join<TracksRun2, aod::McTrackLabels>;
  using McV0Tracks = soa::Join<aod::V0Datas, aod::McV0Labels>;

  void processDataRun3(CollisionsRun3::iterator const& collision, aod::BCsWithTimestamps const&, aod::V0Datas const& V0s, Tracks const& tracks)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    fillLambdaRecoTables<kRun3, kData>(collision, bc, V0s, tracks);
  }
  PROCESS_SWITCH(LambdaTableProducer, processDataRun3, "Process for Run3 DATA", true);

  void processMCRun3(
    aod::McCollisions::iterator const& mcCollision,
    soa::SmallGroups<soa::Join<CollisionsRun3, aod::McCollisionLabels>> const& collisions,
    aod::BCsWithTimestamps const& bcts, McV0Tracks const& V0s,
    TracksMC const& tracks, aod::McParticles const& mcParticles)
  {
    analyzeMcRecoGen<kRun3, kMC>(mcCollision, collisions, bcts, V0s, tracks, mcParticles);
  }
  PROCESS_SWITCH(LambdaTableProducer, processMCRun3, "Process for Run3 MC RecoGen", false);
};

struct LambdaTracksExtProducer {

  Produces<aod::LambdaTracksExt> lambdaTrackExtTable;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    const AxisSpec axisMult(10, 0, 10);
    const AxisSpec axisMass(200, 1.08, 1.18, "M_{p#pi} (GeV/#it{c}^{2})");
    const AxisSpec axisCPA(100, 0.995, 1.0, "cos(#theta_{PA})");
    const AxisSpec axisDcaDau(75, 0., 1.5, "Daughter DCA (#sigma)");
    const AxisSpec axisDEta(320, -1.6, 1.6, "#Delta#eta");
    const AxisSpec axisDPhi(640, -PIHalf, 3. * PIHalf, "#Delta#varphi");

    histos.add("h1i_totlambda_mult", "Multiplicity", kTH1I, {axisMult});
    histos.add("h1i_totantilambda_mult", "Multiplicity", kTH1I, {axisMult});
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
    static constexpr std::array<std::string_view, 2> SubDir = {"Reco/", "SharingDau/"};
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
    int nTotLambda = 0, nTotAntiLambda = 0;

    for (auto const& lambda : tracks) {
      bool lambdaSharingDauFlag = false;

      if (lambda.v0Type() == kLambda) {
        ++nTotLambda;
      } else if (lambda.v0Type() == kAntiLambda) {
        ++nTotAntiLambda;
      }

      for (auto const& track : tracks) {
        if (lambda.index() == track.index()) {
          continue;
        }

        if (lambda.posTrackId() == track.posTrackId() || lambda.negTrackId() == track.negTrackId()) {
          lambdaSharingDauFlag = true;

          if (lambda.v0Type() == kLambda && track.v0Type() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaP"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.v0Type() == kLambda && track.v0Type() == kLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaP_LaP"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          } else if (lambda.v0Type() == kAntiLambda && track.v0Type() == kAntiLambda) {
            histos.fill(HIST("h2d_n2_etaphi_LaM_LaM"), lambda.eta() - track.eta(), RecoDecay::constrainAngle((lambda.phi() - track.phi()), -PIHalf));
          }
        }
      }

      if (lambdaSharingDauFlag) {
        fillHistos<kLambdaShareDau>(lambda);
      } else {
        fillHistos<kUniqueLambda>(lambda);
      }

      lambdaTrackExtTable(lambdaSharingDauFlag);
    }

    if (nTotLambda != 0) {
      histos.fill(HIST("h1i_totlambda_mult"), nTotLambda);
    }
    if (nTotAntiLambda != 0) {
      histos.fill(HIST("h1i_totantilambda_mult"), nTotAntiLambda);
    }
  }
};

struct LambdaSpinPolarization {
  Produces<aod::LambdaMixEventCollisions> lambdaMixEvtCol;
  Produces<aod::LambdaMixEventTracks> lambdaMixEvtTrk;
  Produces<aod::LambdaMixEventMcGenCollisions> lambdaMixEvtMGCol;
  Produces<aod::LambdaMixEventMcGenTracks> lambdaMixEvtMGTrk;

  Configurable<float> cMassAccMin{"cMassAccMin", 1.0956f, "Accepted mass min (GeV/c2); candidates outside the window are not paired"};
  Configurable<float> cMassAccMax{"cMassAccMax", 1.1356f, "Accepted mass max (GeV/c2); candidates outside the window are not paired"};
  Configurable<float> cMassHistMin{"cMassHistMin", 1.0806f, "Mass axes min (GeV/c2), below cMassAccMin"};
  Configurable<float> cMassHistMax{"cMassHistMax", 1.1806f, "Mass axes max (GeV/c2), above cMassAccMax"};
  Configurable<int> cNMassBins{"cNMassBins", 200, "Mass axes N bins"};
  ConfigurableAxis axisDeltaR{"axisDeltaR", {VARIABLE_WIDTH, 0.0, 0.4, 0.8, 1.2, 1.8, 2.4, 3.1, 3.5}, "DeltaR bins (last bin: outside the analysis windows)"};
  Configurable<int> cNBinsCosTS{"cNBinsCosTS", 10, "N costheta* bins"};
  Configurable<int> cSharedDauRule{"cSharedDauRule", 0, "Candidates sharing a daughter track: 0 = keep the largest set without shared tracks per cluster (ties: closest to the PDG mass), 1 = veto only the pairs sharing a track, 2 = drop every candidate sharing a track"};
  Configurable<int> cMinCandPerEvent{"cMinCandPerEvent", 2, "Collisions with fewer reconstructed Lambda + AntiLambda candidates are skipped by the pair, mixing and derived-table processes (the event QA counts all collisions)"};

  Configurable<int> cMEPoolDepth{"cMEPoolDepth", 100000, "Rolling ME pool: depth in donor events per (centrality, Vz) bin"};
  Configurable<int> cMEPoolMinEvents{"cMEPoolMinEvents", 20, "Rolling ME pool: events a pool must hold before its bin is mixed (warm-up)"};
  Configurable<int> cMEPoolMinCand{"cMEPoolMinCand", 2, "Rolling ME pool: only events with at least this many accepted candidates donate"};
  ConfigurableAxis axisCentME{"axisCentME", {VARIABLE_WIDTH, 0, 10, 20, 40, 100}, "ME pool centrality bins"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10., -7., -4., -2., 0, 2., 4., 7., 10.}, "ME pool vtxZ bins"};

  Configurable<float> cMaxDeltaPt{"cMaxDeltaPt", 0.1f, "Kinematic ME: max |deltaPt(SE)-deltaPt(ME)| (GeV/c)"};
  Configurable<float> cMaxDeltaPhi{"cMaxDeltaPhi", 0.02f, "Kinematic ME: max |deltaPhi(SE)-deltaPhi(ME)| (rad)"};
  Configurable<float> cMaxDeltaRap{"cMaxDeltaRap", 0.03f, "Kinematic ME: max |deltaRap(SE)-deltaRap(ME)|"};
  Configurable<bool> cFillSEInMixing{"cFillSEInMixing", true, "Kinematic ME: fill every same-event pair the mixing processes, matched or not (SEproc)"};

  Configurable<bool> cFillKinWSpectra{"cFillKinWSpectra", true, "Kinematic ME: fill the leg-2 spectra of SEproc and ME (weight-map input, closure)"};
  Configurable<float> cMinPt{"cMinPt", 0.6f, "Leg spectra and ME pool grid: pT min (GeV/c)"};
  Configurable<float> cMaxPt{"cMaxPt", 3.0f, "Leg spectra and ME pool grid: pT max (GeV/c)"};
  Configurable<float> cMinRap{"cMinRap", -0.5f, "Leg spectra and ME pool grid: rapidity min"};
  Configurable<float> cMaxRap{"cMaxRap", 0.5f, "Leg spectra and ME pool grid: rapidity max"};
  Configurable<int> cNKinWPtBins{"cNKinWPtBins", 0, "Leg spectra: pT bins (0: no wider than cMaxDeltaPt)"};
  Configurable<int> cNKinWRapBins{"cNKinWRapBins", 0, "Leg spectra: rapidity bins (0: no wider than cMaxDeltaRap)"};
  Configurable<int> cNKinWPhiBins{"cNKinWPhiBins", 0, "Leg spectra: phi bins (0: no wider than cMaxDeltaPhi)"};
  Configurable<bool> cApplyKinWeights{"cApplyKinWeights", false, "Kinematic ME: weight the leg-2 replacements with the maps from CCDB"};
  Configurable<std::string> cKinWeightCcdbUrl{"cKinWeightCcdbUrl", "http://alice-ccdb.cern.ch", "Kinematic weights: CCDB url"};
  Configurable<std::string> cKinWeightCcdbPath{"cKinWeightCcdbPath", "", "Kinematic weights: CCDB path of a TList with hKinW_<pair> (axes of KinW/*/hLeg2_<pair>)"};
  Configurable<int64_t> cKinWeightCcdbTimestamp{"cKinWeightCcdbTimestamp", -1, "Kinematic weights: CCDB timestamp (-1: timestamp of the first mixed collision)"};

  Configurable<bool> cFillClosePairQA{"cFillClosePairQA", true, "Close-pair QA of the same-event pairs: SE (processDataReco) and SEproc (processDataRecoMixed)"};
  Configurable<float> cClosePairMassWindow{"cClosePairMassWindow", 0.004f, "Close-pair QA: both legs within |m - m_Lambda| < this (GeV/c2); <= 0: all pairs"};
  ConfigurableAxis axisClosePair{"axisClosePair", {100, -0.15, 0.15}, "Close-pair QA: bins of the daughter Delta eta, Delta y and Delta phi"};
  Configurable<bool> cFillSharedDauQA{"cFillSharedDauQA", true, "Shared-daughter QA: SE and SEproc pairs with a candidate sharing a daughter track with another candidate, and the cluster sizes"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  uint64_t pairOrderSeed = 0;
  std::array<std::unique_ptr<THnBase>, 4> kinWeightMaps{};
  std::vector<TAxis> legAxesRef;
  bool kinWeightsLoaded = false;

  struct SharedDauCand {
    int64_t globalIdx;
    int64_t posId;
    int64_t negId;
    float dMass;
  };
  static constexpr std::size_t MaxExactCluster = 16;
  static constexpr int NoSharedTrackIdx = -1;
  std::vector<SharedDauCand> sdCands;
  std::vector<int64_t> sdConflicted;
  std::vector<int64_t> sdDropped;

  struct PoolTrack {
    float mPx, mPy, mPz, mPt, mRap, mPhi, mMass;
    float mPrPx, mPrPy, mPrPz;
    int8_t mSpecies;
    [[nodiscard]] float px() const { return mPx; }
    [[nodiscard]] float py() const { return mPy; }
    [[nodiscard]] float pz() const { return mPz; }
    [[nodiscard]] float pt() const { return mPt; }
    [[nodiscard]] float rap() const { return mRap; }
    [[nodiscard]] float phi() const { return mPhi; }
    [[nodiscard]] float mass() const { return mMass; }
    [[nodiscard]] float prPx() const { return mPrPx; }
    [[nodiscard]] float prPy() const { return mPrPy; }
    [[nodiscard]] float prPz() const { return mPrPz; }
    [[nodiscard]] int8_t species() const { return mSpecies; }
  };

  template <typename T>
  PoolTrack toPoolTrack(T const& trk)
  {
    return PoolTrack{trk.px(), trk.py(), trk.pz(), trk.pt(), trk.rap(), trk.phi(), trk.mass(), trk.prPx(), trk.prPy(), trk.prPz(), static_cast<int8_t>(trk.v0Type())};
  }

  class RollingPool
  {
   public:
    struct Grid {
      float ptMin = 0.f;
      float rapMin = 0.f;
      float cellPt = 1.f;
      float cellRap = 1.f;
      float cellPhi = TwoPI;
      int nPt = 1;
      int nRap = 1;
      int nPhi = 1;
    };

    RollingPool() = default;
    RollingPool(Grid const& g, int depth) : grid(g), maxEvents(depth) {}

    [[nodiscard]] std::size_t nEvents() const { return nPerEvent.size(); }

    void push(std::vector<PoolTrack> const& event)
    {
      for (auto const& cand : event) {
        const int64_t seq = seqBase + static_cast<int64_t>(slots.size());
        slots.push_back(Slot{.cand = cand, .next = NoSeq});
        auto [it, isNew] = cells[cand.species()].try_emplace(cellKey(cand), Cell{.head = seq, .tail = seq});
        if (!isNew) {
          slot(it->second.tail).next = seq;
          it->second.tail = seq;
        }
      }
      nPerEvent.push_back(static_cast<int>(event.size()));
      while (static_cast<int>(nPerEvent.size()) > maxEvents) {
        const int n = nPerEvent.front();
        for (int i = 0; i < n; ++i) {
          retireOldest();
        }
        nPerEvent.pop_front();
      }
    }

    template <typename F>
    void findMatches(PoolTrack const& leg, F const& isMatch, std::vector<PoolTrack const*>& out) const
    {
      out.clear();
      int iPt = 0, iRap = 0, iPhi = 0;
      cellCoords(leg, iPt, iRap, iPhi);
      std::array<int, 3> phiCells{};
      int nPhiCells = 0;
      for (int d = -1; d <= 1; ++d) {
        const int k = ((iPhi + d) % grid.nPhi + grid.nPhi) % grid.nPhi;
        if (std::find(phiCells.begin(), phiCells.begin() + nPhiCells, k) == phiCells.begin() + nPhiCells) {
          phiCells[nPhiCells++] = k;
        }
      }
      auto const& speciesCells = cells[leg.species()];
      for (int jPt = std::max(iPt - 1, 0); jPt <= std::min(iPt + 1, grid.nPt - 1); ++jPt) {
        for (int jRap = std::max(iRap - 1, 0); jRap <= std::min(iRap + 1, grid.nRap - 1); ++jRap) {
          for (int ip = 0; ip < nPhiCells; ++ip) {
            const auto it = speciesCells.find(key(jPt, jRap, phiCells[ip]));
            if (it == speciesCells.end()) {
              continue;
            }
            for (int64_t seq = it->second.head; seq != NoSeq;) {
              Slot const& sl = slot(seq);
              if (isMatch(sl.cand)) {
                out.push_back(&sl.cand);
              }
              seq = sl.next;
            }
          }
        }
      }
    }

   private:
    static constexpr int64_t NoSeq = -1;
    struct Slot {
      PoolTrack cand;
      int64_t next;
    };
    struct Cell {
      int64_t head;
      int64_t tail;
    };

    Slot& slot(int64_t seq) { return slots[static_cast<std::size_t>(seq - seqBase)]; }
    [[nodiscard]] Slot const& slot(int64_t seq) const { return slots[static_cast<std::size_t>(seq - seqBase)]; }

    void retireOldest()
    {
      Slot const& old = slots.front();
      auto& speciesCells = cells[old.cand.species()];
      auto it = speciesCells.find(cellKey(old.cand));
      if (it == speciesCells.end() || it->second.head != seqBase) {
        LOGF(fatal, "RollingPool: grid index out of step with the candidate queue");
        return;
      }
      if (old.next == NoSeq) {
        speciesCells.erase(it);
      } else {
        it->second.head = old.next;
      }
      slots.pop_front();
      ++seqBase;
    }

    [[nodiscard]] int64_t key(int iPt, int iRap, int iPhi) const
    {
      return (static_cast<int64_t>(iPt) * grid.nRap + iRap) * grid.nPhi + iPhi;
    }

    void cellCoords(PoolTrack const& c, int& iPt, int& iRap, int& iPhi) const
    {
      iPt = std::clamp(static_cast<int>(std::floor((c.pt() - grid.ptMin) / grid.cellPt)), 0, grid.nPt - 1);
      iRap = std::clamp(static_cast<int>(std::floor((c.rap() - grid.rapMin) / grid.cellRap)), 0, grid.nRap - 1);
      iPhi = std::clamp(static_cast<int>(std::floor(RecoDecay::constrainAngle(c.phi(), 0.f) / grid.cellPhi)), 0, grid.nPhi - 1);
    }

    [[nodiscard]] int64_t cellKey(PoolTrack const& c) const
    {
      int iPt = 0, iRap = 0, iPhi = 0;
      cellCoords(c, iPt, iRap, iPhi);
      return key(iPt, iRap, iPhi);
    }

    Grid grid{};
    int maxEvents = 0;
    std::deque<Slot> slots;
    std::deque<int> nPerEvent;
    int64_t seqBase = 0;
    std::array<std::unordered_map<int64_t, Cell>, 2> cells{};
  };

  std::vector<RollingPool> mePools;
  std::vector<RollingPool> mePoolsGen;
  TAxis mePoolAxisCent;
  TAxis mePoolAxisVz;

  void init(InitContext const&)
  {
    const AxisSpec axisM1(cNMassBins, cMassHistMin, cMassHistMax, "M_{1} (GeV/#it{c}^{2})");
    const AxisSpec axisM2(cNMassBins, cMassHistMin, cMassHistMax, "M_{2} (GeV/#it{c}^{2})");
    const AxisSpec axisDR(axisDeltaR, "#DeltaR");
    const AxisSpec axisCosTS(cNBinsCosTS, -1, 1, "cos(#theta*)");
    const std::vector<AxisSpec> pairAxes = {axisM1, axisM2, axisDR, axisCosTS};

    auto nBins = [](int configured, float range, float window) {
      if (configured > 0) {
        return configured;
      }
      return window > 0.f ? static_cast<int>(std::ceil(range / window)) : 1;
    };
    const AxisSpec axisLegPt(nBins(cNKinWPtBins, cMaxPt - cMinPt, cMaxDeltaPt), cMinPt, cMaxPt, "p_{T} (GeV/#it{c})");
    const AxisSpec axisLegRap(nBins(cNKinWRapBins, cMaxRap - cMinRap, cMaxDeltaRap), cMinRap, cMaxRap, "y");
    const AxisSpec axisLegPhi(nBins(cNKinWPhiBins, TwoPI, cMaxDeltaPhi), 0., TwoPI, "#varphi (rad)");
    const std::vector<AxisSpec> legAxes = {axisLegPt, axisLegRap, axisLegPhi};
    legAxesRef.clear();
    for (auto const& spec : legAxes) {
      const std::vector<double> edges = axisEdges(spec);
      legAxesRef.emplace_back(static_cast<int>(edges.size()) - 1, edges.data());
    }

    const bool massAxesContainWindow = cMassHistMin.value < cMassAccMin.value && cMassAccMax.value < cMassHistMax.value;
    if (!massAxesContainWindow) {
      LOGF(fatal, "The mass axes [%.4f, %.4f] must be wider than the accepted mass window [%.4f, %.4f]", cMassHistMin.value, cMassHistMax.value, cMassAccMin.value, cMassAccMax.value);
    }
    for (auto const& edge : {cMassAccMin.value, cMassAccMax.value}) {
      if (!isOnMassBinEdge(edge)) {
        LOGF(warning, "Accepted mass %.5f is not a bin edge of the mass axes [%.4f, %.4f] / %d", edge, cMassHistMin.value, cMassHistMax.value, cNMassBins.value);
      }
    }
    if (cSharedDauRule.value < kSharedDauMaxSet || cSharedDauRule.value > kSharedDauDropAll) {
      LOGF(fatal, "cSharedDauRule must be %d (largest set without shared tracks), %d (veto pairs sharing a track) or %d (drop every candidate sharing a track)", static_cast<int>(kSharedDauMaxSet), static_cast<int>(kSharedDauPairVeto), static_cast<int>(kSharedDauDropAll));
    }
    if (cMinCandPerEvent.value < 1) {
      LOGF(fatal, "cMinCandPerEvent must be at least 1, it is %d", cMinCandPerEvent.value);
    }

    if (doprocessDataRecoMixed || doprocessMcGenMixed) {
      const std::vector<double> centEdges = axisEdges(AxisSpec(axisCentME));
      const std::vector<double> vzEdges = axisEdges(AxisSpec(axisVtxZME));
      mePoolAxisCent.Set(static_cast<int>(centEdges.size()) - 1, centEdges.data());
      mePoolAxisVz.Set(static_cast<int>(vzEdges.size()) - 1, vzEdges.data());

      static constexpr std::array<int64_t, 3> PoolSizeSteps = {1, 2, 5};
      static constexpr int64_t PoolSizeDecade = 10;
      const int64_t depth = std::max(1, cMEPoolDepth.value);
      std::vector<double> poolSizeEdges = {0.};
      for (int64_t decade = 1; decade < depth; decade *= PoolSizeDecade) {
        for (auto const& step : PoolSizeSteps) {
          if (step * decade < depth) {
            poolSizeEdges.push_back(static_cast<double>(step * decade));
          }
        }
      }
      poolSizeEdges.push_back(static_cast<double>(depth));
      const AxisSpec axisPoolSize(poolSizeEdges, "events in the pool");

      std::vector<std::string> mixQaDirs;
      if (doprocessDataRecoMixed) {
        mixQaDirs.emplace_back("QA/ME/");
      }
      if (doprocessMcGenMixed) {
        mixQaDirs.emplace_back("McGen/QA/ME/");
      }
      for (auto const& dir : mixQaDirs) {
        histos.add((dir + "hMixedCentVz").c_str(), "events mixed;cent (%);V_{z} (cm)", kTH2F, {axisCentME, axisVtxZME});
        auto hStatus = histos.add<TH1>((dir + "hEventStatus").c_str(), "rolling-pool mixing;;collisions", kTH1D, {{kNMEEventStatus - 1, 0.5, static_cast<double>(kNMEEventStatus) - 0.5}});
        hStatus->GetXaxis()->SetBinLabel(kMEEventSeen, "seen");
        hStatus->GetXaxis()->SetBinLabel(kMEEventOutsidePools, "outside the pool bins");
        hStatus->GetXaxis()->SetBinLabel(kMEEventWarmUp, "pairs, pool in warm-up");
        hStatus->GetXaxis()->SetBinLabel(kMEEventMixed, "mixed");
        hStatus->GetXaxis()->SetBinLabel(kMEEventDonated, "donated to the pool");
        histos.add((dir + "hPoolEventsAtMixing").c_str(), "pool size when a collision is mixed;events in the pool;collisions", kTH1F, {axisPoolSize});
        histos.add((dir + "hMatchesLeg2").c_str(), "matches per replaced leg 2;N matches;pair type", kTH2F, {{100, -0.5, 99.5}, {4, -0.5, 3.5}});
        histos.add((dir + "hDeltaRShift").c_str(), "ME pairs are filed at their own #DeltaR;#DeltaR_{seed};#DeltaR_{mixed pair} - #DeltaR_{seed}", kTH2F, {{35, 0., 3.5}, {100, -0.2, 0.2}});
        histos.add((dir + "hLambdaMultVsCent").c_str(), "ME #Lambda mult;cent;N", kTH2D, {axisCentME, {50, 0, 50}});
        histos.add((dir + "hAntiLambdaMultVsCent").c_str(), "ME #bar{#Lambda} mult;cent;N", kTH2D, {axisCentME, {50, 0, 50}});
      }

      if (cMaxDeltaPt.value <= 0.f || cMaxDeltaRap.value <= 0.f || cMaxDeltaPhi.value <= 0.f) {
        LOGF(fatal, "Rolling ME pool: the matching windows cMaxDeltaPt, cMaxDeltaRap and cMaxDeltaPhi must be positive");
      }
      RollingPool::Grid grid;
      grid.ptMin = cMinPt.value;
      grid.rapMin = cMinRap.value;
      grid.cellPt = cMaxDeltaPt.value;
      grid.cellRap = cMaxDeltaRap.value;
      grid.nPt = std::max(1, static_cast<int>(std::ceil((cMaxPt.value - cMinPt.value) / grid.cellPt)));
      grid.nRap = std::max(1, static_cast<int>(std::ceil((cMaxRap.value - cMinRap.value) / grid.cellRap)));
      grid.nPhi = std::max(1, static_cast<int>(std::floor(TwoPI / cMaxDeltaPhi.value)));
      grid.cellPhi = TwoPI / static_cast<float>(grid.nPhi);
      const auto nPools = static_cast<std::size_t>(mePoolAxisCent.GetNbins()) * static_cast<std::size_t>(mePoolAxisVz.GetNbins());
      if (doprocessDataRecoMixed) {
        mePools.assign(nPools, RollingPool(grid, cMEPoolDepth.value));
      }
      if (doprocessMcGenMixed) {
        mePoolsGen.assign(nPools, RollingPool(grid, cMEPoolDepth.value));
      }
      LOGF(info, "Rolling ME pools: %d (centrality) x %d (Vz) bins, depth %d donor events, grid %d x %d x %d (pT, y, phi)",
           mePoolAxisCent.GetNbins(), mePoolAxisVz.GetNbins(), cMEPoolDepth.value, grid.nPt, grid.nRap, grid.nPhi);
    }

    std::vector<std::string> evDirs;
    if (doprocessDataReco || doprocessDataRecoMixed) {
      evDirs.emplace_back("");
    }
    if (doprocessMcGen || doprocessMcGenMixed) {
      evDirs.emplace_back("McGen/");
    }
    const AxisSpec axisCandMass(cNMassBins, cMassHistMin, cMassHistMax, "M_{p#pi} (GeV/#it{c}^{2})");
    for (auto const& dir : evDirs) {
      auto hCount = histos.add<TH1>((dir + "Events/hEventCount").c_str(), "collisions of the pair analysis;;collisions", kTH1D, {{kNPairEventCount - 1, 0.5, static_cast<double>(kNPairEventCount) - 0.5}});
      hCount->GetXaxis()->SetBinLabel(kPeAll, "all");
      hCount->GetXaxis()->SetBinLabel(kPeOneCand, "#geq 1 accepted #Lambda/#bar{#Lambda}");
      hCount->GetXaxis()->SetBinLabel(kPeTwoCand, "#geq 2 accepted (enter pairs)");
      hCount->GetXaxis()->SetBinLabel(kPeUnlikeSign, "#geq 1 #Lambda#bar{#Lambda} pair");
      hCount->GetXaxis()->SetBinLabel(kPeLambdaLambda, "#geq 1 #Lambda#Lambda pair");
      hCount->GetXaxis()->SetBinLabel(kPeAntiLambdaAntiLambda, "#geq 1 #bar{#Lambda}#bar{#Lambda} pair");
      histos.add((dir + "Events/hNCandAccepted").c_str(), "accepted candidates per collision;N_{#Lambda};N_{#bar{#Lambda}}", kTH2D, {{21, -0.5, 20.5}, {21, -0.5, 20.5}});
      histos.add((dir + "Events/hCentPairEvents").c_str(), "collisions with a pair;cent (%);collisions", kTH1D, {{100, 0., 100.}});
      for (auto const& species : {"Lambda/", "AntiLambda/"}) {
        const std::string cand = dir + "QA/PairCand/" + species;
        histos.add((cand + "hMassVsPt").c_str(), "accepted candidates of collisions with a pair", kTH2F, {axisCandMass, axisLegPt});
        histos.add((cand + "hRapVsPhi").c_str(), "accepted candidates of collisions with a pair", kTH2F, {axisLegRap, axisLegPhi});
      }
    }

    static constexpr std::array<std::string_view, 4> PairTags = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};
    for (auto const& tag : PairTags) {
      const std::string pair{tag};
      if (doprocessDataReco) {
        histos.add(("SE/hPair_" + pair).c_str(), ("SE " + pair).c_str(), kTHnSparseF, pairAxes);
      }
      if (doprocessDataRecoMixed) {
        histos.add(("ME/hPair_" + pair).c_str(), ("ME " + pair).c_str(), kTHnSparseD, pairAxes, true);
      }
      if (doprocessDataRecoMixed && cFillSEInMixing) {
        histos.add(("SEproc/hPair_" + pair).c_str(), ("SE processed by the mixing " + pair).c_str(), kTHnSparseF, pairAxes);
      }
      if (doprocessDataRecoMixed && cFillKinWSpectra) {
        histos.add(("KinW/SEproc/hLeg2_" + pair).c_str(), ("leg 2 of SEproc " + pair).c_str(), kTHnSparseD, legAxes);
        histos.add(("KinW/ME/hLeg2_" + pair).c_str(), ("candidates replacing leg 2 " + pair).c_str(), kTHnSparseD, legAxes, true);
      }
      if (cFillSharedDauQA && doprocessDataReco) {
        histos.add(("QA/SharedDau/SE/hPair_" + pair).c_str(), ("SE pairs with a candidate sharing a daughter " + pair).c_str(), kTHnSparseF, pairAxes);
      }
      if (cFillSharedDauQA && doprocessDataRecoMixed) {
        histos.add(("QA/SharedDau/SEproc/hPair_" + pair).c_str(), ("SEproc pairs with a candidate sharing a daughter " + pair).c_str(), kTHnSparseF, pairAxes);
      }
      if (doprocessMcGen) {
        histos.add(("McGen/SE/hPair_" + pair).c_str(), ("MC generated SE " + pair).c_str(), kTHnSparseF, pairAxes);
      }
      if (doprocessMcGenMixed) {
        histos.add(("McGen/ME/hPair_" + pair).c_str(), ("MC generated ME " + pair).c_str(), kTHnSparseD, pairAxes, true);
      }
      if (doprocessMcGenMixed && cFillSEInMixing) {
        histos.add(("McGen/SEproc/hPair_" + pair).c_str(), ("MC generated SE processed by the mixing " + pair).c_str(), kTHnSparseF, pairAxes);
      }
    }
    if (cFillSharedDauQA && (doprocessDataReco || doprocessDataRecoMixed)) {
      histos.add("QA/SharedDau/hClusterSize", "clusters of candidates sharing daughter tracks;candidates in the cluster;clusters", kTH1F, {{20, 0.5, 20.5}});
    }

    if (cFillClosePairQA && (doprocessDataReco || doprocessDataRecoMixed)) {
      const AxisSpec axisCPdEta(axisClosePair, "#Delta#eta");
      const AxisSpec axisCPdY(axisClosePair, "#Deltay");
      const AxisSpec axisCPdPhi(axisClosePair, "#Delta#varphi (rad)");
      std::vector<std::pair<std::string, std::string>> cpDirs;
      if (doprocessDataReco) {
        cpDirs.emplace_back("QA/ClosePair/", "");
      }
      if (doprocessDataRecoMixed) {
        cpDirs.emplace_back("QA/ClosePairSEproc/", "SEproc, ");
      }
      for (auto const& [dir, tag] : cpDirs) {
        for (auto const& cls : {"LamLam", "ALamALam", "UnlikeSign"}) {
          for (auto const& combo : {"pp", "pipi", "p1pi2", "pi1p2"}) {
            const std::string name = dir + cls + "/h_" + combo;
            const std::string title = tag + cls + ", daughters " + combo;
            histos.add((name + "_dEta").c_str(), title.c_str(), kTH3F, {axisDR, axisCPdEta, axisCPdPhi});
            histos.add((name + "_dY").c_str(), title.c_str(), kTH3F, {axisDR, axisCPdY, axisCPdPhi});
          }
        }
      }
    }

    if (cApplyKinWeights) {
      if (cKinWeightCcdbPath.value.empty()) {
        LOGF(fatal, "cApplyKinWeights is set but cKinWeightCcdbPath is empty");
      }
      ccdb->setURL(cKinWeightCcdbUrl.value);
      ccdb->setCaching(true);
      ccdb->setFatalWhenNull(false);
    }
  }

  bool isAcceptedMass(float m) const
  {
    return m >= cMassAccMin.value && m < cMassAccMax.value;
  }

  static std::vector<double> axisEdges(AxisSpec const& axis)
  {
    if (!axis.nBins.has_value()) {
      return axis.binEdges;
    }
    const int n = axis.nBins.value();
    std::vector<double> edges(n + 1);
    for (int i = 0; i <= n; ++i) {
      edges[i] = axis.binEdges[0] + (axis.binEdges[1] - axis.binEdges[0]) * i / n;
    }
    return edges;
  }

  static constexpr double MassEdgeTolerance = 1e-3;
  [[nodiscard]] bool isOnMassBinEdge(double m) const
  {
    const double width = (static_cast<double>(cMassHistMax.value) - cMassHistMin.value) / cNMassBins.value;
    const double k = (m - cMassHistMin.value) / width;
    return std::abs(k - std::round(k)) < MassEdgeTolerance;
  }

  void getBoostVector(std::array<float, 4> const& p, std::array<float, 3>& v)
  {
    int n = p.size();
    for (int i = 0; i < n - 1; ++i) {
      v[i] = -p[i] / RecoDecay::e(p[0], p[1], p[2], p[3]);
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

  static float cosOpeningAngle(std::array<float, 4> const& a, std::array<float, 4> const& b)
  {
    std::array<float, 3> n1 = {a[0], a[1], a[2]}, n2 = {b[0], b[1], b[2]};
    return RecoDecay::dotProd(n1, n2) /
           (RecoDecay::sqrtSumOfSquares(n1[0], n1[1], n1[2]) *
            RecoDecay::sqrtSumOfSquares(n2[0], n2[1], n2[2]));
  }

  float cosThetaStarAtlas(std::array<float, 4> const& l1, std::array<float, 4> const& l2, std::array<float, 4> const& pr1, std::array<float, 4> const& pr2)
  {
    auto l1a = l1;
    auto l2a = l2;
    auto pr1a = pr1;
    auto pr2a = pr2;

    const float mPair = RecoDecay::m(std::array{std::array{l1a[0], l1a[1], l1a[2]}, std::array{l2a[0], l2a[1], l2a[2]}}, std::array{l1a[3], l2a[3]});
    std::array<float, 4> llpair = {l1a[0] + l2a[0], l1a[1] + l2a[1], l1a[2] + l2a[2], mPair};
    std::array<float, 3> vPair{};
    getBoostVector(llpair, vPair);
    boost(l1a, vPair);
    boost(l2a, vPair);
    boost(pr1a, vPair);
    boost(pr2a, vPair);
    std::array<float, 3> v1p{}, v2p{};
    getBoostVector(l1a, v1p);
    getBoostVector(l2a, v2p);
    boost(pr1a, v1p);
    boost(pr2a, v2p);
    return cosOpeningAngle(pr1a, pr2a);
  }

  static float deltaR(PoolTrack const& p1, PoolTrack const& p2)
  {
    const float drap = p1.rap() - p2.rap();
    const float dphi = RecoDecay::constrainAngle(p1.phi() - p2.phi(), -PI);
    return std::sqrt(drap * drap + dphi * dphi);
  }

  template <PairTier tier, ParticlePairType part_pair>
  void fillPair(PoolTrack const& p1, PoolTrack const& p2, float dR, float w)
  {
    static constexpr std::array<std::string_view, 8> TierDir = {"SE/hPair_", "SEproc/hPair_", "ME/hPair_",
                                                                "McGen/SE/hPair_", "McGen/SEproc/hPair_", "McGen/ME/hPair_",
                                                                "QA/SharedDau/SE/hPair_", "QA/SharedDau/SEproc/hPair_"};
    static constexpr std::array<std::string_view, 4> PairTag = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};

    std::array<float, 4> l1 = {p1.px(), p1.py(), p1.pz(), p1.mass()};
    std::array<float, 4> l2 = {p2.px(), p2.py(), p2.pz(), p2.mass()};
    std::array<float, 4> pr1 = {p1.prPx(), p1.prPy(), p1.prPz(), MassProton};
    std::array<float, 4> pr2 = {p2.prPx(), p2.prPy(), p2.prPz(), MassProton};

    const float ctheta = cosThetaStarAtlas(l1, l2, pr1, pr2);
    histos.fill(HIST(TierDir[tier]) + HIST(PairTag[part_pair]), p1.mass(), p2.mass(), dR, ctheta, w);
  }

  template <bool IsME, ParticlePairType part_pair>
  void fillLeg2Spectrum(PoolTrack const& leg, float w)
  {
    static constexpr std::array<std::string_view, 2> Dir = {"KinW/SEproc/hLeg2_", "KinW/ME/hLeg2_"};
    static constexpr std::array<std::string_view, 4> PairTag = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};
    histos.fill(HIST(Dir[IsME]) + HIST(PairTag[part_pair]), leg.pt(), leg.rap(), leg.phi(), w);
  }

  template <PairTier tier, int Cls, int Combo>
  void fillDaughterSeparation(std::array<float, 3> const& a, double ma, std::array<float, 3> const& b, double mb, float dR)
  {
    static_assert(tier == kSameEvent || tier == kSameEventInMixing, "close-pair QA is filled for SE and SEproc pairs only");
    static constexpr std::array<std::string_view, 2> TierDir = {"QA/ClosePair/", "QA/ClosePairSEproc/"};
    static constexpr std::array<std::string_view, 3> ClsDir = {"LamLam/", "ALamALam/", "UnlikeSign/"};
    static constexpr std::array<std::string_view, 4> ComboName = {"h_pp_", "h_pipi_", "h_p1pi2_", "h_pi1p2_"};
    const float dEta = std::asinh(a[2] / std::hypot(a[0], a[1])) - std::asinh(b[2] / std::hypot(b[0], b[1]));
    const float dY = RecoDecay::y(a, ma) - RecoDecay::y(b, mb);
    const float dPhi = RecoDecay::constrainAngle(std::atan2(a[1], a[0]) - std::atan2(b[1], b[0]), -PI);
    histos.fill(HIST(TierDir[tier]) + HIST(ClsDir[Cls]) + HIST(ComboName[Combo]) + HIST("dEta"), dR, dEta, dPhi);
    histos.fill(HIST(TierDir[tier]) + HIST(ClsDir[Cls]) + HIST(ComboName[Combo]) + HIST("dY"), dR, dY, dPhi);
  }

  template <PairTier tier, ParticlePairType part_pair>
  void fillClosePairQA(PoolTrack const& l1, PoolTrack const& l2, float dR)
  {
    constexpr int Cls = part_pair == kLambdaLambda ? 0 : (part_pair == kAntiLambdaAntiLambda ? 1 : 2);
    const std::array<float, 3> pr1 = {l1.prPx(), l1.prPy(), l1.prPz()};
    const std::array<float, 3> pr2 = {l2.prPx(), l2.prPy(), l2.prPz()};
    const std::array<float, 3> pi1 = {l1.px() - l1.prPx(), l1.py() - l1.prPy(), l1.pz() - l1.prPz()};
    const std::array<float, 3> pi2 = {l2.px() - l2.prPx(), l2.py() - l2.prPy(), l2.pz() - l2.prPz()};
    fillDaughterSeparation<tier, Cls, 0>(pr1, MassProton, pr2, MassProton, dR);
    fillDaughterSeparation<tier, Cls, 1>(pi1, MassPionCharged, pi2, MassPionCharged, dR);
    fillDaughterSeparation<tier, Cls, 2>(pr1, MassProton, pi2, MassPionCharged, dR);
    fillDaughterSeparation<tier, Cls, 3>(pi1, MassPionCharged, pr2, MassProton, dR);
  }

  bool isClosePairQAMass(float m) const
  {
    return cClosePairMassWindow.value <= 0.f || std::abs(m - MassLambda0) < cClosePairMassWindow.value;
  }

  bool isKinematicMatch(PoolTrack const& cand, PoolTrack const& leg) const
  {
    return std::abs(cand.pt() - leg.pt()) < cMaxDeltaPt.value &&
           std::abs(cand.rap() - leg.rap()) < cMaxDeltaRap.value &&
           std::abs(RecoDecay::constrainAngle(cand.phi() - leg.phi(), -PI)) < cMaxDeltaPhi.value;
  }

  static constexpr std::size_t NLegAxes = 3;
  template <ParticlePairType part_pair>
  float kinWeight(PoolTrack const& cand) const
  {
    THnBase const* map = kinWeightMaps[part_pair].get();
    if (!map) {
      return 1.f;
    }
    const std::array<double, NLegAxes> x = {cand.pt(), cand.rap(), cand.phi()};
    std::array<int, NLegAxes> idx{};
    for (std::size_t i = 0; i < NLegAxes; ++i) {
      TAxis const* axis = map->GetAxis(i);
      idx[i] = axis->FindFixBin(x[i]);
      if (idx[i] < 1 || idx[i] > axis->GetNbins()) {
        return 1.f;
      }
    }
    const int64_t bin = map->GetBin(idx.data());
    if (bin < 0) {
      return 1.f;
    }
    const double w = map->GetBinContent(bin);
    return w > 0. ? static_cast<float>(w) : 1.f;
  }

  static constexpr double AxisEdgeTolerance = 1e-6;
  [[nodiscard]] bool hasLegAxes(THnBase const& map) const
  {
    if (map.GetNdimensions() != static_cast<int>(legAxesRef.size())) {
      return false;
    }
    for (std::size_t i = 0; i < legAxesRef.size(); ++i) {
      TAxis const* axis = map.GetAxis(i);
      TAxis const& ref = legAxesRef[i];
      if (axis->GetNbins() != ref.GetNbins()) {
        return false;
      }
      for (int b = 1; b <= ref.GetNbins() + 1; ++b) {
        if (std::abs(axis->GetBinLowEdge(b) - ref.GetBinLowEdge(b)) > AxisEdgeTolerance) {
          return false;
        }
      }
    }
    return true;
  }

  void loadKinWeights(uint64_t collisionTimestamp)
  {
    static constexpr std::array<std::string_view, 4> PairTag = {"LaPLaM", "LaMLaP", "LaPLaP", "LaMLaM"};
    const int64_t ts = cKinWeightCcdbTimestamp.value >= 0 ? cKinWeightCcdbTimestamp.value : static_cast<int64_t>(collisionTimestamp);
    auto* list = ccdb->getForTimeStamp<TList>(cKinWeightCcdbPath.value, ts);
    if (!list) {
      LOGP(fatal, "Kinematic weights: no object at {} for timestamp {}", cKinWeightCcdbPath.value, ts);
      return;
    }
    auto const* guard = dynamic_cast<TNamed const*>(list->FindObject("kinWeightGuardStatus"));
    if (guard != nullptr && !TString(guard->GetTitle()).BeginsWith("OK")) {
      LOGF(fatal, "Kinematic weights at %s are flagged by their builder: %s", cKinWeightCcdbPath.value.c_str(), guard->GetTitle());
      return;
    }
    for (std::size_t i = 0; i < PairTag.size(); ++i) {
      const std::string name = "hKinW_" + std::string(PairTag[i]);
      auto* map = dynamic_cast<THnBase*>(list->FindObject(name.c_str()));
      if (!map || !hasLegAxes(*map)) {
        LOGF(fatal, "Kinematic weights: %s is missing, or its axes differ from KinW/*/hLeg2_%s (pT, y, phi) of this configuration", name.c_str(), std::string(PairTag[i]).c_str());
        return;
      }
      kinWeightMaps[i].reset(dynamic_cast<THnBase*>(map->Clone()));
    }
    auto const* form = dynamic_cast<TNamed const*>(list->FindObject("kinWeightForm"));
    kinWeightsLoaded = true;
    LOGF(info, "Kinematic weights loaded from %s (form: %s)", cKinWeightCcdbPath.value.c_str(), form ? form->GetTitle() : "not stated");
  }

  static uint64_t splitMix(uint64_t x)
  {
    x += 0x9E3779B97F4A7C15ULL;
    x = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
  }

  template <typename C>
  void setPairOrderSeed(C const& col)
  {
    const uint64_t vz = std::bit_cast<uint32_t>(col.posZ());
    pairOrderSeed = splitMix(((vz << 32) | static_cast<uint32_t>(col.globalIndex())) ^ (col.timeStamp() * 0xD6E8FEB86659FD93ULL));
  }

  template <typename T>
  bool isCountedOrder(T const& trk_1, T const& trk_2) const
  {
    const int64_t g1 = trk_1.globalIndex(), g2 = trk_2.globalIndex();
    if (g1 == g2) {
      return false;
    }
    const bool firstIsLo = g1 < g2;
    const uint64_t lo = static_cast<uint32_t>(firstIsLo ? g1 : g2);
    const uint64_t hi = static_cast<uint32_t>(firstIsLo ? g2 : g1);
    const bool swapped = splitMix(pairOrderSeed ^ ((lo << 32) | hi)) & 1ULL;
    return firstIsLo != swapped;
  }

  template <typename T1, typename T2>
  static bool sharesDaughter(T1 const& a, T2 const& b)
  {
    return a.posTrackId() == b.posTrackId() || a.negTrackId() == b.negTrackId();
  }

  [[nodiscard]] bool sharesDaughterCand(std::size_t a, std::size_t b) const
  {
    return a != b && (sdCands[a].posId == sdCands[b].posId || sdCands[a].negId == sdCands[b].negId);
  }

  [[nodiscard]] bool isDropped(int64_t globalIdx) const
  {
    return std::binary_search(sdDropped.begin(), sdDropped.end(), globalIdx);
  }

  [[nodiscard]] bool isConflicted(int64_t globalIdx) const
  {
    return std::binary_search(sdConflicted.begin(), sdConflicted.end(), globalIdx);
  }

  template <bool Reco, typename T>
  [[nodiscard]] bool isKept(T const& trk) const
  {
    return !Reco || !isDropped(trk.globalIndex());
  }

  template <bool Reco, typename T>
  [[nodiscard]] bool isUsable(T const& trk) const
  {
    return isAcceptedMass(trk.mass()) && isKept<Reco>(trk);
  }

  template <typename T>
  [[nodiscard]] bool hasMinCandidates(T const& lTrks, T const& alTrks) const
  {
    return std::cmp_greater_equal(lTrks.size() + alTrks.size(), cMinCandPerEvent.value);
  }

  void selectLargestFreeSet(std::vector<std::size_t> const& members, std::vector<char>& keep) const
  {
    const std::size_t s = members.size();
    keep.assign(s, 0);
    if (s <= MaxExactCluster) {
      std::array<uint32_t, MaxExactCluster> conf{};
      for (std::size_t a = 0; a < s; ++a) {
        for (std::size_t b = 0; b < s; ++b) {
          if (sharesDaughterCand(members[a], members[b])) {
            conf[a] |= static_cast<uint32_t>(1) << b;
          }
        }
      }
      uint32_t best = 0;
      int bestN = -1;
      float bestScore = 0.f;
      const uint32_t nMasks = static_cast<uint32_t>(1) << s;
      for (uint32_t mask = 1; mask < nMasks; ++mask) {
        bool isFree = true;
        int nIn = 0;
        float score = 0.f;
        for (std::size_t a = 0; a < s && isFree; ++a) {
          if (((mask >> a) & 1U) == 0U) {
            continue;
          }
          if ((conf[a] & mask) != 0U) {
            isFree = false;
            continue;
          }
          ++nIn;
          score += sdCands[members[a]].dMass;
        }
        if (isFree && (nIn > bestN || (nIn == bestN && score < bestScore))) {
          best = mask;
          bestN = nIn;
          bestScore = score;
        }
      }
      for (std::size_t a = 0; a < s; ++a) {
        keep[a] = static_cast<char>((best >> a) & 1U);
      }
      return;
    }
    std::vector<std::size_t> order(s);
    for (std::size_t a = 0; a < s; ++a) {
      order[a] = a;
    }
    std::sort(order.begin(), order.end(), [this, &members](std::size_t x, std::size_t y) {
      SharedDauCand const& cx = sdCands[members[x]];
      SharedDauCand const& cy = sdCands[members[y]];
      return cx.dMass < cy.dMass || (cx.dMass == cy.dMass && cx.globalIdx < cy.globalIdx);
    });
    for (auto const& a : order) {
      bool isFree = true;
      for (std::size_t b = 0; b < s && isFree; ++b) {
        if (keep[b] != 0 && sharesDaughterCand(members[a], members[b])) {
          isFree = false;
        }
      }
      if (isFree) {
        keep[a] = 1;
      }
    }
  }

  template <typename T>
  void resolveSharedDaughters(T const& lTrks, T const& alTrks, bool fillQA)
  {
    sdCands.clear();
    sdConflicted.clear();
    sdDropped.clear();
    auto collect = [this](auto const& trks) {
      for (auto const& trk : trks) {
        if (trk.lambdaSharingDaughter()) {
          sdCands.push_back(SharedDauCand{trk.globalIndex(), trk.posTrackId(), trk.negTrackId(), static_cast<float>(std::abs(trk.mass() - MassLambda0))});
        }
      }
    };
    collect(lTrks);
    collect(alTrks);
    const std::size_t n = sdCands.size();
    if (n == 0) {
      return;
    }
    std::sort(sdCands.begin(), sdCands.end(), [](SharedDauCand const& a, SharedDauCand const& b) { return a.globalIdx < b.globalIdx; });
    std::vector<int> cluster(n, -1);
    std::vector<std::size_t> members;
    std::vector<std::size_t> stack;
    std::vector<char> keep;
    int nClusters = 0;
    for (std::size_t seed = 0; seed < n; ++seed) {
      if (cluster[seed] >= 0) {
        continue;
      }
      members.clear();
      stack.assign(1, seed);
      cluster[seed] = nClusters;
      while (!stack.empty()) {
        const std::size_t a = stack.back();
        stack.pop_back();
        members.push_back(a);
        for (std::size_t b = 0; b < n; ++b) {
          if (cluster[b] < 0 && sharesDaughterCand(a, b)) {
            cluster[b] = nClusters;
            stack.push_back(b);
          }
        }
      }
      ++nClusters;
      std::sort(members.begin(), members.end());
      if (fillQA && cFillSharedDauQA) {
        histos.fill(HIST("QA/SharedDau/hClusterSize"), members.size());
      }
      if (members.size() < NCandidatesForPair) {
        continue;
      }
      for (auto const& a : members) {
        sdConflicted.push_back(sdCands[a].globalIdx);
      }
      if (cSharedDauRule.value == kSharedDauPairVeto) {
        continue;
      }
      if (cSharedDauRule.value == kSharedDauDropAll) {
        keep.assign(members.size(), 0);
      } else {
        selectLargestFreeSet(members, keep);
      }
      for (std::size_t a = 0; a < members.size(); ++a) {
        if (keep[a] == 0) {
          sdDropped.push_back(sdCands[members[a]].globalIdx);
        }
      }
    }
    std::sort(sdConflicted.begin(), sdConflicted.end());
    std::sort(sdDropped.begin(), sdDropped.end());
  }

  template <PairTier tier, ParticlePairType part_pair, typename T>
  void analyzePairsSE(T const& trks_1, T const& trks_2)
  {
    constexpr bool Reco = tier == kSameEvent;
    for (auto const& trk_1 : trks_1) {
      if (!isUsable<Reco>(trk_1)) {
        continue;
      }
      const PoolTrack p1 = toPoolTrack(trk_1);

      for (auto const& trk_2 : trks_2) {
        if (!isCountedOrder(trk_1, trk_2) || !isUsable<Reco>(trk_2)) {
          continue;
        }
        if constexpr (Reco) {
          if (sharesDaughter(trk_1, trk_2)) {
            continue;
          }
        }
        const PoolTrack p2 = toPoolTrack(trk_2);
        const float dR = deltaR(p1, p2);
        fillPair<tier, part_pair>(p1, p2, dR, 1.f);
        if constexpr (Reco) {
          if (cFillClosePairQA && isClosePairQAMass(p1.mass()) && isClosePairQAMass(p2.mass())) {
            fillClosePairQA<kSameEvent, part_pair>(p1, p2, dR);
          }
          if (cFillSharedDauQA && (isConflicted(trk_1.globalIndex()) || isConflicted(trk_2.globalIndex()))) {
            fillPair<kSameEventSharedDau, part_pair>(p1, p2, dR, 1.f);
          }
        }
      }
    }
  }

  template <bool Gen, typename T>
  void fillEventQA(float centrality, T const& lTrks, T const& alTrks)
  {
    static constexpr std::array<std::string_view, 2> EvDir = {"Events/", "McGen/Events/"};
    std::size_t nLambda = 0, nAntiLambda = 0;
    for (auto const& trk : lTrks) {
      if (isUsable<!Gen>(trk)) {
        ++nLambda;
      }
    }
    for (auto const& trk : alTrks) {
      if (isUsable<!Gen>(trk)) {
        ++nAntiLambda;
      }
    }
    histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeAll);
    histos.fill(HIST(EvDir[Gen]) + HIST("hNCandAccepted"), nLambda, nAntiLambda);
    if (nLambda + nAntiLambda >= 1) {
      histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeOneCand);
    }
    if (nLambda + nAntiLambda < NCandidatesForPair) {
      return;
    }
    histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeTwoCand);
    histos.fill(HIST(EvDir[Gen]) + HIST("hCentPairEvents"), centrality);
    if (nLambda >= 1 && nAntiLambda >= 1) {
      histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeUnlikeSign);
    }
    if (nLambda >= NCandidatesForPair) {
      histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeLambdaLambda);
    }
    if (nAntiLambda >= NCandidatesForPair) {
      histos.fill(HIST(EvDir[Gen]) + HIST("hEventCount"), kPeAntiLambdaAntiLambda);
    }
    for (auto const& trk : lTrks) {
      if (isUsable<!Gen>(trk)) {
        fillPairCandQA<Gen, kLambda>(trk);
      }
    }
    for (auto const& trk : alTrks) {
      if (isUsable<!Gen>(trk)) {
        fillPairCandQA<Gen, kAntiLambda>(trk);
      }
    }
  }

  template <bool Gen, ParticleType part, typename T>
  void fillPairCandQA(T const& trk)
  {
    static constexpr std::array<std::string_view, 2> CandDir = {"QA/PairCand/", "McGen/QA/PairCand/"};
    static constexpr std::array<std::string_view, 2> Species = {"Lambda/", "AntiLambda/"};
    histos.fill(HIST(CandDir[Gen]) + HIST(Species[part]) + HIST("hMassVsPt"), trk.mass(), trk.pt());
    histos.fill(HIST(CandDir[Gen]) + HIST(Species[part]) + HIST("hRapVsPhi"), trk.rap(), trk.phi());
  }

  [[nodiscard]] int poolBin(float centrality, float vz) const
  {
    const int iCent = mePoolAxisCent.FindFixBin(centrality);
    const int iVz = mePoolAxisVz.FindFixBin(vz);
    if (iCent < 1 || iCent > mePoolAxisCent.GetNbins() || iVz < 1 || iVz > mePoolAxisVz.GetNbins()) {
      return -1;
    }
    return (iCent - 1) * mePoolAxisVz.GetNbins() + (iVz - 1);
  }

  template <bool Gen, ParticlePairType part_pair, typename T>
  void analyzePairsMEKinematic(T const& se_trks_1, T const& se_trks_2, RollingPool const& pool)
  {
    static constexpr std::array<std::string_view, 2> MixQaDir = {"QA/ME/", "McGen/QA/ME/"};
    constexpr PairTier SeTier = Gen ? kMcGenSameEventInMixing : kSameEventInMixing;
    constexpr PairTier MeTier = Gen ? kMcGenMixedEvent : kMixedEvent;
    std::vector<PoolTrack const*> matches;
    for (auto const& trk1 : se_trks_1) {
      if (!isUsable<!Gen>(trk1)) {
        continue;
      }
      const PoolTrack p1 = toPoolTrack(trk1);

      for (auto const& trk2 : se_trks_2) {
        if (!isCountedOrder(trk1, trk2) || !isUsable<!Gen>(trk2)) {
          continue;
        }
        if constexpr (!Gen) {
          if (sharesDaughter(trk1, trk2)) {
            continue;
          }
        }
        const PoolTrack p2 = toPoolTrack(trk2);
        const float dRSeed = deltaR(p1, p2);

        if (cFillSEInMixing) {
          fillPair<SeTier, part_pair>(p1, p2, dRSeed, 1.f);
        }
        if constexpr (!Gen) {
          if (cFillClosePairQA && isClosePairQAMass(p1.mass()) && isClosePairQAMass(p2.mass())) {
            fillClosePairQA<kSameEventInMixing, part_pair>(p1, p2, dRSeed);
          }
          if (cFillSharedDauQA && (isConflicted(trk1.globalIndex()) || isConflicted(trk2.globalIndex()))) {
            fillPair<kSameEventInMixingSharedDau, part_pair>(p1, p2, dRSeed, 1.f);
          }
          if (cFillKinWSpectra) {
            fillLeg2Spectrum<false, part_pair>(p2, 1.f);
          }
        }

        pool.findMatches(p2, [this, &p2](PoolTrack const& cand) { return isKinematicMatch(cand, p2); }, matches);
        histos.fill(HIST(MixQaDir[Gen]) + HIST("hMatchesLeg2"), matches.size(), part_pair);
        for (auto const& meP2 : matches) {
          float kinW = 1.f;
          if constexpr (!Gen) {
            kinW = kinWeight<part_pair>(*meP2);
          }
          const float w = kinW / static_cast<float>(matches.size());
          const float dRMixed = deltaR(p1, *meP2);
          fillPair<MeTier, part_pair>(p1, *meP2, dRMixed, w);
          histos.fill(HIST(MixQaDir[Gen]) + HIST("hDeltaRShift"), dRSeed, dRMixed - dRSeed, w);
          if constexpr (!Gen) {
            if (cFillKinWSpectra) {
              fillLeg2Spectrum<true, part_pair>(*meP2, w);
            }
          }
        }
      }
    }
  }

  template <bool Gen, typename C, typename T>
  void mixCollision(C const& col, T const& lTrks, T const& alTrks, std::vector<RollingPool>& pools, std::vector<PoolTrack>& eventCands)
  {
    static constexpr std::array<std::string_view, 2> MixQaDir = {"QA/ME/", "McGen/QA/ME/"};
    histos.fill(HIST(MixQaDir[Gen]) + HIST("hEventStatus"), kMEEventSeen);
    const int bin = poolBin(col.cent(), col.posZ());
    if (bin < 0) {
      histos.fill(HIST(MixQaDir[Gen]) + HIST("hEventStatus"), kMEEventOutsidePools);
      return;
    }
    std::size_t nLambda = 0, nAntiLambda = 0;
    for (auto const& trk : lTrks) {
      if (isKept<!Gen>(trk)) {
        ++nLambda;
      }
    }
    for (auto const& trk : alTrks) {
      if (isKept<!Gen>(trk)) {
        ++nAntiLambda;
      }
    }
    histos.fill(HIST(MixQaDir[Gen]) + HIST("hLambdaMultVsCent"), col.cent(), nLambda);
    histos.fill(HIST(MixQaDir[Gen]) + HIST("hAntiLambdaMultVsCent"), col.cent(), nAntiLambda);

    eventCands.clear();
    for (auto const& trk : lTrks) {
      if (isUsable<!Gen>(trk)) {
        eventCands.push_back(toPoolTrack(trk));
      }
    }
    for (auto const& trk : alTrks) {
      if (isUsable<!Gen>(trk)) {
        eventCands.push_back(toPoolTrack(trk));
      }
    }

    RollingPool& pool = pools[bin];
    if (eventCands.size() >= NCandidatesForPair) {
      if (static_cast<int>(pool.nEvents()) >= cMEPoolMinEvents.value) {
        if constexpr (!Gen) {
          if (cApplyKinWeights && !kinWeightsLoaded) {
            loadKinWeights(col.timeStamp());
          }
        }
        setPairOrderSeed(col);
        histos.fill(HIST(MixQaDir[Gen]) + HIST("hEventStatus"), kMEEventMixed);
        histos.fill(HIST(MixQaDir[Gen]) + HIST("hMixedCentVz"), col.cent(), col.posZ());
        histos.fill(HIST(MixQaDir[Gen]) + HIST("hPoolEventsAtMixing"), pool.nEvents());
        analyzePairsMEKinematic<Gen, kLambdaAntiLambda>(lTrks, alTrks, pool);
        analyzePairsMEKinematic<Gen, kAntiLambdaLambda>(alTrks, lTrks, pool);
        analyzePairsMEKinematic<Gen, kLambdaLambda>(lTrks, lTrks, pool);
        analyzePairsMEKinematic<Gen, kAntiLambdaAntiLambda>(alTrks, alTrks, pool);
      } else {
        histos.fill(HIST(MixQaDir[Gen]) + HIST("hEventStatus"), kMEEventWarmUp);
      }
    }
    if (static_cast<int>(eventCands.size()) >= cMEPoolMinCand.value) {
      pool.push(eventCands);
      histos.fill(HIST(MixQaDir[Gen]) + HIST("hEventStatus"), kMEEventDonated);
    }
  }

  using LambdaCollisions = aod::LambdaCollisions;
  using LambdaMcGenCollisions = aod::LambdaMcGenCollisions;
  using LambdaTracks = soa::Join<aod::LambdaTracks, aod::LambdaTracksExt>;

  Preslice<LambdaTracks> perCollisionLambda = aod::lambdatrack::lambdaCollisionId;
  SliceCache cache;

  Partition<LambdaTracks> partLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  Partition<LambdaTracks> partAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  SliceCache cachemc;

  Partition<aod::LambdaMcGenTracks> partMcLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  Partition<aod::LambdaMcGenTracks> partMcAntiLambdaTracks = (aod::lambdatrack::v0Type == (int8_t)kAntiLambda) && (aod::lambdatrack::v0PrmScd == (int8_t)kPrimary);

  void processDummy(LambdaCollisions::iterator const&) {}
  PROCESS_SWITCH(LambdaSpinPolarization, processDummy, "Dummy", false);

  void processDataReco(LambdaCollisions::iterator const& collision, LambdaTracks const&)
  {
    setPairOrderSeed(collision);
    auto lTrks = partLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto alTrks = partAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    resolveSharedDaughters(lTrks, alTrks, !doprocessDataRecoMixed);
    if (!doprocessDataRecoMixed) {
      fillEventQA<false>(collision.cent(), lTrks, alTrks);
    }
    if (!hasMinCandidates(lTrks, alTrks)) {
      return;
    }

    analyzePairsSE<kSameEvent, kLambdaAntiLambda>(lTrks, alTrks);
    analyzePairsSE<kSameEvent, kAntiLambdaLambda>(alTrks, lTrks);
    analyzePairsSE<kSameEvent, kLambdaLambda>(lTrks, lTrks);
    analyzePairsSE<kSameEvent, kAntiLambdaAntiLambda>(alTrks, alTrks);
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataReco, "SE only (data/MCReco)", false);

  void processDataRecoMixed(LambdaCollisions const& cols, LambdaTracks const&)
  {
    std::vector<PoolTrack> eventCands;
    for (auto const& col : cols) {
      auto lTrks = partLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, col.globalIndex(), cache);
      auto alTrks = partAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, col.globalIndex(), cache);
      resolveSharedDaughters(lTrks, alTrks, true);
      fillEventQA<false>(col.cent(), lTrks, alTrks);
      if (!hasMinCandidates(lTrks, alTrks)) {
        continue;
      }
      mixCollision<false>(col, lTrks, alTrks, mePools, eventCands);
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataRecoMixed, "ME: kinematic mixing against rolling pools that cross DataFrames", false);

  void processMcGen(LambdaMcGenCollisions::iterator const& collision, aod::LambdaMcGenTracks const&)
  {
    setPairOrderSeed(collision);
    auto lTrks = partMcLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, collision.globalIndex(), cachemc);
    auto alTrks = partMcAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, collision.globalIndex(), cachemc);
    if (!doprocessMcGenMixed) {
      fillEventQA<true>(collision.cent(), lTrks, alTrks);
    }
    if (!hasMinCandidates(lTrks, alTrks)) {
      return;
    }

    analyzePairsSE<kMcGenSameEvent, kLambdaAntiLambda>(lTrks, alTrks);
    analyzePairsSE<kMcGenSameEvent, kAntiLambdaLambda>(alTrks, lTrks);
    analyzePairsSE<kMcGenSameEvent, kLambdaLambda>(lTrks, lTrks);
    analyzePairsSE<kMcGenSameEvent, kAntiLambdaAntiLambda>(alTrks, alTrks);
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processMcGen, "MC generated: SE pairs (McGen/SE), as processDataReco", false);

  void processMcGenMixed(LambdaMcGenCollisions const& cols, aod::LambdaMcGenTracks const&)
  {
    std::vector<PoolTrack> eventCands;
    for (auto const& col : cols) {
      auto lTrks = partMcLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, col.globalIndex(), cachemc);
      auto alTrks = partMcAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, col.globalIndex(), cachemc);
      fillEventQA<true>(col.cent(), lTrks, alTrks);
      if (!hasMinCandidates(lTrks, alTrks)) {
        continue;
      }
      mixCollision<true>(col, lTrks, alTrks, mePoolsGen, eventCands);
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processMcGenMixed, "MC generated: kinematic mixing against rolling pools, as processDataRecoMixed", false);

  void processDataRecoMixEvent(LambdaCollisions::iterator const& collision, LambdaTracks const&)
  {
    auto lTrks = partLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);
    auto alTrks = partAntiLambdaTracks->sliceByCached(aod::lambdatrack::lambdaCollisionId, collision.globalIndex(), cache);

    if (!hasMinCandidates(lTrks, alTrks)) {
      return;
    }

    lambdaMixEvtCol(collision.index(), collision.cent(), collision.posZ(), collision.timeStamp());

    auto fillMixEventTrack = [&](auto const& track) {
      const bool shared = track.lambdaSharingDaughter();
      lambdaMixEvtTrk(collision.index(), track.globalIndex(), track.px(),
                      track.py(), track.pz(), track.mass(), track.prPx(),
                      track.prPy(), track.prPz(), track.v0Type(),
                      collision.timeStamp(),
                      shared ? static_cast<int>(track.posTrackId()) : NoSharedTrackIdx,
                      shared ? static_cast<int>(track.negTrackId()) : NoSharedTrackIdx);
    };
    for (auto const& track : lTrks) {
      fillMixEventTrack(track);
    }
    for (auto const& track : alTrks) {
      fillMixEventTrack(track);
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processDataRecoMixEvent, "Mix-event table filling", false);

  void processMcGenMixEvent(LambdaMcGenCollisions::iterator const& collision, aod::LambdaMcGenTracks const&)
  {
    auto lTrks = partMcLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, collision.globalIndex(), cachemc);
    auto alTrks = partMcAntiLambdaTracks->sliceByCached(aod::lambdamcgentrack::lambdaMcGenCollisionId, collision.globalIndex(), cachemc);

    if (!hasMinCandidates(lTrks, alTrks)) {
      return;
    }

    lambdaMixEvtMGCol(collision.globalIndex(), collision.cent(),
                      collision.posZ(), collision.timeStamp());

    for (auto const& track : lTrks) {
      lambdaMixEvtMGTrk(collision.globalIndex(), track.globalIndex(),
                        track.px(), track.py(), track.pz(), track.mass(),
                        track.prPx(), track.prPy(), track.prPz(),
                        track.v0Type(), collision.timeStamp());
    }
    for (auto const& track : alTrks) {
      lambdaMixEvtMGTrk(collision.globalIndex(), track.globalIndex(),
                        track.px(), track.py(), track.pz(), track.mass(),
                        track.prPx(), track.prPy(), track.prPz(),
                        track.v0Type(), collision.timeStamp());
    }
  }
  PROCESS_SWITCH(LambdaSpinPolarization, processMcGenMixEvent, "Mix-event McGen table filling", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LambdaTableProducer>(cfgc),
                      adaptAnalysisTask<LambdaTracksExtProducer>(cfgc),
                      adaptAnalysisTask<LambdaSpinPolarization>(cfgc)};
}
