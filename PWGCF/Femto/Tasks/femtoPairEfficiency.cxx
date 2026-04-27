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

/// \file femtoPairEfficiency.cxx
/// \brief task for analysis of charged-particle pseudorapidity density at midrapidity in pp collisions triggered with a paritcle pair
/// \author Anton Riedel (anton.riedel@tum.ch)

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>

using o2::framework::Configurable;
using std::string;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection,
                                  aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe,
                                  aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe>;
using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;
using V0TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;

enum {
  kTrackTypebegin = 0,
  kGlobalplusITS = 1,
  kGlobalonly,
  kITSonly,
  kTrackTypeend
};

enum {
  kGenpTbegin = 0,
  kNoGenpTVar = 1,
  kGenpTup,
  kGenpTdown,
  kGenpTend
};

enum {
  kGenTrkTypebegin = 0,
  kGenAll = 1,
  kGenPion,
  kGenKaon,
  kGenProton,
  kGenOther,
  kGenTrkTypeend
};

enum {
  kRecTrkTypebegin = 0,
  kRecoAll = 1,
  kRecoPion,
  kRecoKaon,
  kRecoProton,
  kRecoOther,
  kRecoSecondary,
  kRecoWeakDecay,
  kRecoFake,
  kRecoBkg,
  kRecTrkTypeend
};

static constexpr TrackSelectionFlags::flagtype TrackSelectionIts =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;
static constexpr TrackSelectionFlags::flagtype TrackSelectionTpc =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDca =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;
static constexpr TrackSelectionFlags::flagtype TrackSelectionDcaxyOnly =
  TrackSelectionFlags::kDCAxy;

AxisSpec axisEvent{15, 0.5, 15.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisPhi{{0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "#phi", "PhiAxis"};
AxisSpec axisPhi2{629, 0, o2::constants::math::TwoPI, "#phi"};
AxisSpec axisCent{100, 0, 100, "#Cent"};
AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
AxisSpec axisGenTrkType = {kGenTrkTypeend - 1, +kGenTrkTypebegin + 0.5, +kGenTrkTypeend - 0.5, "", "GenTrackTypeAxis"};
AxisSpec axisRecTrkType = {kRecTrkTypeend - 1, +kRecTrkTypebegin + 0.5, +kRecTrkTypeend - 0.5, "", "RecTrackTypeAxis"};
AxisSpec axisMassK0s = {200, 0.4, 0.6, "K0sMass", "K0sMass"};
AxisSpec axisMassLambda = {200, 1.07, 1.17, "Lambda/AntiLamda Mass", "Lambda/AntiLamda Mass"};
auto static constexpr KminCharge = 3.f;
auto static constexpr KminPtCut = 0.1f;

struct FemtoPairEfficiency {

  o2::framework::SliceCache cache;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> v0radiusK0SCut{"v0radiusK0SCut", 1.2f, "K0S RadiusCut"};
  Configurable<float> dcapostopvK0SCut{"dcapostopvK0SCut", 0.05f, "K0S dcapostopvCut"};
  Configurable<float> dcanegtopvK0SCut{"dcanegtopvK0SCut", 0.05f, "K0S dcanegtopvCut"};
  Configurable<float> v0cospaK0SCut{"v0cospaK0SCut", 0.995f, "K0S v0cospaCut"};
  Configurable<float> dcav0daughterK0Scut{"dcav0daughterK0Scut", 1.0f, "K0S dcav0daughtercut"};
  Configurable<float> minTPCnClsK0SCut{"minTPCnClsK0SCut", 50.0f, "K0S minTPCnClsCut"};
  Configurable<float> nSigmaTpcK0SCut{"nSigmaTpcK0SCut", 5.0f, "K0S nSigmaTpcCut"};
  Configurable<float> v0etaK0SCut{"v0etaK0SCut", 0.9f, "K0S v0etaCut"};
  Configurable<float> v0radiusLambdaCut{"v0radiusLambdaCut", 1.2f, "Lambda RadiusCut"};
  Configurable<float> dcapostopvLambdaCut{"dcapostopvLambdaCut", 0.05f, "Lambda dcapostopvCut"};
  Configurable<float> dcanegtopvLambdaCut{"dcanegtopvLambdaCut", 0.05f, "Lambda dcanegtopvCut"};
  Configurable<float> v0cospaLambdaCut{"v0cospaLambdaCut", 0.995f, "Lambda v0cospaCut"};
  Configurable<float> dcav0daughterLambdacut{"dcav0daughterLambdacut", 1.0f, "Lambda dcav0daughtercut"};
  Configurable<float> minTPCnClsLambdaCut{"minTPCnClsLambdaCut", 50.0f, "Lambda minTPCnClsCut"};
  Configurable<float> nSigmaTpcLambdaCut{"nSigmaTpcLambdaCut", 5.0f, "Lambda nSigmaTpcCut"};
  Configurable<float> v0etaLambdaCut{"v0etaLambdaCut", 0.9f, "Lambda v0etaCut"};
  Configurable<float> extraphicut1{"extraphicut1", 3.07666f, "Extra Phi cut 1"};
  Configurable<float> extraphicut2{"extraphicut2", 3.12661f, "Extra Phi cut 2"};
  Configurable<float> extraphicut3{"extraphicut3", 0.03f, "Extra Phi cut 3"};
  Configurable<float> extraphicut4{"extraphicut4", 6.253f, "Extra Phi cut 4"};
  ConfigurableAxis multHistBin{"multHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis pvHistBin{"pvHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis fv0aMultHistBin{"fv0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0aMultHistBin{"ft0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0cMultHistBin{"ft0cMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ptHistBin{"ptHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis centralityBinning{"centralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};
  ConfigurableAxis binsMult{"binsMult", {500, 0.0f, +500.0f}, ""};
  ConfigurableAxis binsDCA{"binsDCA", {500, -10.0f, 10.0f}, ""};

  Configurable<bool> isApplyTFcut{"isApplyTFcut", true, "Enable TimeFrameBorder cut"};
  Configurable<bool> isApplyITSROcut{"isApplyITSROcut", true, "Enable ITS ReadOutFrameBorder cut"};
  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Enable INEL > 0 condition"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<bool> isApplyPairCut{"isApplyPairCut", false, "Apply pair cut"};
  Configurable<float> kStarMax{"kStarMax", 0.2, "k* threshold for pair cut"};

  struct : ConfigurableGroup {
    string prefix = string("TrackSel1");
    Configurable<float> tpcClustersMin{"tpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};
    Configurable<float> tpcCrossedRowsMin{"tpcCrossedRowsMin", {80.f}, "Minimum number of crossed rows in TPC"};
    Configurable<float> tpcClustersOverCrossedRows{"tpcClustersOverCrossedRows", {0.83f}, "Minimum fraction of clusters over crossed rows in TPC"};
    Configurable<float> tpcSharedClustersMax{"tpcSharedClustersMax", {160.f}, "Maximum number of shared clusters in TPC"};
    Configurable<float> tpcSharedClusterFractionMax{"tpcSharedClusterFractionMax", {1.f}, "Maximum fraction of shared clusters in TPC"};
    Configurable<float> itsClustersMin{"itsClustersMin", {5.f}, "Minimum number of clusters in ITS"};
    Configurable<float> itsIbClustersMin{"itsIbClustersMin", {3.f}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
    Configurable<std::vector<float>> dcaxyMax{"dcaxyMax", {0.004, 0.013, -1}, "Maximum |dca_xy| as a function of pT (v[0] + v[1]*pow(x, v[2])"};
    Configurable<std::vector<float>> dcazMax{"dcazMax", {0.004, 0.013, -1}, "Maximum |dca_z| as a function of pT (v[0] + v[1]*pow(x, v[2])"};
    Configurable<float> pidThres{"pidThres", 0.8, "Minimum momentum to required TOF PID"};
    Configurable<float> tpcNsigma{"tpcNsigma", 3, "Maximum NSigma"};
    Configurable<float> tpctofNsigma{"tpctofNsigma", 3, "Maximum NSigma TPCTOF"};
    Configurable<int> pdgCode{"pdgCode", 2212, "Pdg Code"};
  } TrackSel1;

  struct : ConfigurableGroup {
    string prefix = string("TrackSel2");
    Configurable<float> tpcClustersMin{"tpcClustersMin", {90.f}, "Minimum number of clusters in TPC"};
    Configurable<float> tpcCrossedRowsMin{"tpcCrossedRowsMin", {80.f}, "Minimum number of crossed rows in TPC"};
    Configurable<float> tpcClustersOverCrossedRows{"tpcClustersOverCrossedRows", {0.83f}, "Minimum fraction of clusters over crossed rows in TPC"};
    Configurable<float> tpcSharedClustersMax{"tpcSharedClustersMax", {160.f}, "Maximum number of shared clusters in TPC"};
    Configurable<float> tpcSharedClusterFractionMax{"tpcSharedClusterFractionMax", {1.f}, "Maximum fraction of shared clusters in TPC"};
    Configurable<float> itsClustersMin{"itsClustersMin", {5.f}, "Minimum number of clusters in ITS"};
    Configurable<float> itsIbClustersMin{"itsIbClustersMin", {3.f}, "Minimum number of clusters in inner barrel (max 3) of ITS"};
    Configurable<std::vector<float>> dcaxyMax{"dcaxyMax", {0.004, 0.013, -1}, "Maximum |dca_xy| as a function of pT (v[0] + v[1]*pow(x, v[2])"};
    Configurable<std::vector<float>> dcazMax{"dcazMax", {0.004, 0.013, -1}, "Maximum |dca_z| as a function of pT (v[0] + v[1]*pow(x, v[2])"};
    Configurable<float> pidThres{"pidThres", 0.8, "Minimum momentum to required TOF PID"};
    Configurable<float> tpcNsigma{"tpcNsigma", 3, "Maximum NSigma"};
    Configurable<float> tpctofNsigma{"tpctofNsigma", 3, "Maximum NSigma TPCTOF"};
    Configurable<int> pdgCode{"pdgCode", 2212, "Pdg Code"};
  } TrackSel2;

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin, "Mult", "MultAxis"};
    AxisSpec axisPV = {pvHistBin, "PV", "PVAxis"};
    AxisSpec axisFv0aMult = {fv0aMultHistBin, "fv0a", "FV0AMultAxis"};
    AxisSpec axisFt0aMult = {ft0aMultHistBin, "ft0a", "FT0AMultAxis"};
    AxisSpec axisFt0cMult = {ft0cMultHistBin, "ft0c", "FT0CMultAxis"};
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec dcaAxis = {binsDCA, "DCA vs PV"};
    AxisSpec multAxis = {binsMult, "Multiplicity #eta<0.5"};
    AxisSpec axisKstar = {200, 0., 1., "k* (GeV/c)", "KstarAxis"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "kIsTriggerTVX");
    x->SetBinLabel(3, "kNoTimeFrameBorder");
    x->SetBinLabel(4, "kNoITSROFrameBorder");
    x->SetBinLabel(5, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(6, "kIsGoodZvtxFT0vsPV"); // small difference between z-vertex from PV and from FT0
    x->SetBinLabel(7, "INEL > 0");
    x->SetBinLabel(8, "hasPair");

    histos.add("hTrack1Pt", "Track1 pT", kTH1D, {AxisSpec{ptHistBin, "pT", "pTAxis"}}, false);
    histos.add("hTrack1Eta", "Track1 eta", kTH1D, {axisEta}, false);
    histos.add("hTrack1Phi", "Track1 phi", kTH1D, {axisPhi2}, false);
    histos.add("hTrack2Pt", "Track2 pT", kTH1D, {AxisSpec{ptHistBin, "pT", "pTAxis"}}, false);
    histos.add("hTrack2Eta", "Track2 eta", kTH1D, {axisEta}, false);
    histos.add("hTrack2Phi", "Track2 phi", kTH1D, {axisPhi2}, false);
    histos.add("hPairKstar", "Pair k*", kTH1D, {axisKstar}, false);

    histos.add("hTrack1PtMC", "Track1 pT MC", kTH1D, {AxisSpec{ptHistBin, "pT", "pTAxis"}}, false);
    histos.add("hTrack1EtaMC", "Track1 eta MC", kTH1D, {axisEta}, false);
    histos.add("hTrack1PhiMC", "Track1 phi MC", kTH1D, {axisPhi2}, false);
    histos.add("hTrack2PtMC", "Track2 pT MC", kTH1D, {AxisSpec{ptHistBin, "pT", "pTAxis"}}, false);
    histos.add("hTrack2EtaMC", "Track2 eta MC", kTH1D, {axisEta}, false);
    histos.add("hTrack2PhiMC", "Track2 phi MC", kTH1D, {axisPhi2}, false);
    histos.add("hPairKstarMC", "Pair k* MC", kTH1D, {axisKstar}, false);

    if (doprocessData) {
      histos.add("hdcaxy", "dca to pv in the xy plane", kTH1D, {dcaAxis}, false);
      histos.add("hdcaz", "dca to pv in the z axis", kTH1D, {dcaAxis}, false);
      histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hdatadndeta", "hdatadndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisTrackType}, false);
      histos.add("hdatadndetaMB", "hdatadndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi}, false);
    }

    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisMult, axisFt0aMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFt0cMult}, true);
      histos.add("GlobalMult_vs_FV0A", "GlobalMult_vs_FV0A", kTH2F, {axisMult, axisFv0aMult}, true);
      histos.add("GlobalMult_vs_NPVtracks", "GlobalMult_vs_NPVtracks", kTH2F, {axisMult, axisPV}, true);
      histos.add("NPVtracks_vs_FT0C", "NPVtracks_vs_FT0C", kTH2F, {axisPV, axisFt0cMult}, true);
    }

    if (doprocessStrangeYield) {
      histos.add("hzvtxcent", "hzvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("K0sCentEtaMass", "K0sCentEtaMass", kTH3D, {centAxis, axisEta, axisMassK0s}, false);
      histos.add("LambdaCentEtaMass", "LambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
      histos.add("AntiLambdaCentEtaMass", "AntiLambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
    }

    if (doprocessMCeff) {
      histos.add("hmcdcaxy", "dca to pv in the xy plane", kTH1D, {dcaAxis}, false);
      histos.add("hmcdcaz", "dca to pv in the z axis", kTH1D, {dcaAxis}, false);
      histos.add("hGenMCvertexZ", "hGenMCvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hGenMCvtxzcent", "hGenMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hGenMCAssoRecvertexZ", "hGenMCAssoRecvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hGenMCAssoRecvtxzcent", "hGenMCAssoRecvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hGenMCdndeta", "hGenMCdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi}, false);
      histos.add("hGenMCAssoRecdndeta", "hGenMCAssoRecdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hRecMCvertexZ", "hRecMCvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hRecMCvtxzcent", "hRecMCvtxzcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("hRecMCcentrality", "hRecMCcentrality", kTH1D, {axisCent}, false);
      histos.add("hRecMCphivseta", "hRecMCphivseta", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hRecMCdndeta", "hRecMCdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisRecTrkType}, false);
      histos.add("hRecMCdndetaMB", "hRecMCdndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi, axisRecTrkType}, false);
      histos.add("hGenMCAssoRecdndetaMB", "hGenMCAssoRecdndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi, axisGenTrkType}, false);
    }

    if (doprocessEvtLossSigLossMC) {
      histos.add("MCEventHist", "MCEventHist", kTH1F, {axisEvent}, false);
      auto hstat = histos.get<TH1>(HIST("MCEventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All MC events");
      x->SetBinLabel(2, "MC events with reco event after event selection");
      histos.add("hMultEta05Gen", "multiplicity in eta<0.5 of generated MC events", kTH1F, {multAxis});
      histos.add("hMultEta05GenAssoRec", "multiplicity in eta<0.5 of selected MC events", kTH1F, {multAxis});
      histos.add("hgendndetaVsMultEta05BeforeEvtSel", "hgendndetaBeforeEvtSel vs multiplicity in eta<0.5", kTH2F, {axisEta, multAxis});
      histos.add("hgendndetaVsMultEta05AfterEvtSel", "hgendndetaAfterEvtSel vs multiplicity in eta<0.5", kTH2F, {axisEta, multAxis});
      histos.add("hGenCent", "Centrality of generated MC events", kTH1F, {axisCent});
      histos.add("hGenAssoRecCent", "Centrality of selected MC events", kTH1F, {axisCent});
      histos.add("hgendndetaBeforeEvtSel", "Eta of all generated particles", kTH1F, {axisEta});
      histos.add("hgendndetaAfterEvtSel", "Eta of generated particles after EvtSel", kTH1F, {axisEta});
      histos.add("hgendndetaVscentBeforeEvtSel", "hgendndetaBeforeEvtSel vs centrality", kTH2F, {axisEta, centAxis});
      histos.add("hgendndetaVscentAfterEvtSel", "hgendndetaAfterEvtSel vs centrality", kTH2F, {axisEta, centAxis});
    }
  }

  template <typename T, typename S>
  bool checkTrack(T const& track, S const& sel)
  {
    if (track.tpcNClsFound() < sel.tpcClustersMin.value)
      return false;
    if (track.tpcNClsCrossedRows() < sel.tpcCrossedRowsMin.value)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < sel.tpcClustersOverCrossedRows.value)
      return false;
    if (track.tpcNClsShared() > sel.tpcSharedClustersMax.value)
      return false;
    if (track.tpcFractionSharedCls() > sel.tpcSharedClusterFractionMax.value)
      return false;
    if (track.itsNCls() < sel.itsClustersMin.value)
      return false;
    if (track.itsNClsInnerBarrel() < sel.itsIbClustersMin.value)
      return false;

    const auto& dxy = sel.dcaxyMax.value;
    const auto& dz = sel.dcazMax.value;
    float dcaXyMax = dxy.at(0) + dxy.at(1) * std::pow(track.pt(), dxy.at(2));
    float dcaZMax = dz.at(0) + dz.at(1) * std::pow(track.pt(), dz.at(2));

    if (std::abs(track.dcaXY()) > dcaXyMax)
      return false;
    if (std::abs(track.dcaZ()) > dcaZMax)
      return false;

    return true;
  }

  template <typename T, typename S>
  bool checkTrackPid(T const& track, S const& sel)
  {
    float nsigmaTPC = -99;
    float nsigmaTPCTOF = -99;
    switch (std::abs(sel.pdgCode.value)) {
      case 211:
        nsigmaTPC = track.tpcNSigmaPi();
        nsigmaTPCTOF = std::hypot(track.tpcNSigmaPi(), track.tofNSigmaPi());
        break;
      case 321:
        nsigmaTPC = track.tpcNSigmaKa();
        nsigmaTPCTOF = std::hypot(track.tpcNSigmaKa(), track.tofNSigmaKa());
        break;
      case 2212:
        nsigmaTPC = track.tpcNSigmaPr();
        nsigmaTPCTOF = std::hypot(track.tpcNSigmaPr(), track.tofNSigmaPr());
        break;
      case 1000010020:
        nsigmaTPC = track.tpcNSigmaDe();
        nsigmaTPCTOF = std::hypot(track.tpcNSigmaDe(), track.tofNSigmaDe());
        break;
      default:
        LOG(fatal) << "Unsupported PDG code: " << std::abs(sel.pdgCode.value);
        return false;
    }
    if (track.p() < sel.pidThres.value) {
      if (std::abs(nsigmaTPC) > sel.tpcNsigma.value)
        return false;
    } else {
      if (nsigmaTPCTOF > sel.tpctofNsigma.value)
        return false;
    }
    return true;
  }

  float getKstar(ROOT::Math::PtEtaPhiMVector const& part1, ROOT::Math::PtEtaPhiMVector const& part2)
  {
    // compute pair momentum
    auto sum = part1 + part2;
    // Boost particle 1 to the pair rest frame (Prf) and calculate k* (would be equivalent using particle 2)
    // make a copy of particle 1
    auto particle1Prf = ROOT::Math::PtEtaPhiMVector(part1);
    // get lorentz boost into pair rest frame
    ROOT::Math::Boost boostPrf(sum.BoostToCM());
    // boost particle 1 into pair rest frame and calculate its momentum, which has the same value as k*
    return static_cast<float>(boostPrf(particle1Prf).P());
  }

  template <typename T>
  bool hasPair(const T& tracks)
  {
    bool hasPair = false;
    float mass1 = o2::analysis::femto::utils::getPdgMass(TrackSel1.pdgCode.value);
    float mass2 = o2::analysis::femto::utils::getPdgMass(TrackSel2.pdgCode.value);
    float kstar;

    ROOT::Math::PtEtaPhiMVector particle1;
    ROOT::Math::PtEtaPhiMVector particle2;

    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsUpperIndexPolicy(tracks, tracks))) {

      bool order1 = checkTrack(p1, TrackSel1) && checkTrackPid(p1, TrackSel1) &&
                    checkTrack(p2, TrackSel2) && checkTrackPid(p2, TrackSel2);
      bool order2 = checkTrack(p1, TrackSel2) && checkTrackPid(p1, TrackSel2) &&
                    checkTrack(p2, TrackSel1) && checkTrackPid(p2, TrackSel1);

      if (!order1 && !order2) {
        continue;
      }

      if (order1) {
        particle1 = ROOT::Math::PtEtaPhiMVector(p1.pt(), p1.eta(), p1.phi(), mass1);
        particle2 = ROOT::Math::PtEtaPhiMVector(p2.pt(), p2.eta(), p2.phi(), mass2);
      } else {
        particle1 = ROOT::Math::PtEtaPhiMVector(p2.pt(), p2.eta(), p2.phi(), mass1);
        particle2 = ROOT::Math::PtEtaPhiMVector(p1.pt(), p1.eta(), p1.phi(), mass2);
      }

      kstar = getKstar(particle1, particle2);

      if (kstar > kStarMax.value) {
        continue;
      }
      // fill histograms
      histos.fill(HIST("hTrack1Pt"), particle1.Pt());
      histos.fill(HIST("hTrack1Eta"), particle1.Eta());
      histos.fill(HIST("hTrack1Phi"), RecoDecay::constrainAngle(particle1.Phi()));
      histos.fill(HIST("hTrack2Pt"), particle2.Pt());
      histos.fill(HIST("hTrack2Eta"), particle2.Eta());
      histos.fill(HIST("hTrack2Phi"), RecoDecay::constrainAngle(particle2.Phi()));
      histos.fill(HIST("hPairKstar"), kstar);

      hasPair = true;
    }
    return hasPair;
  }

  template <typename T>
  bool hasPairMC(const T& tracks)
  {
    float mass1 = o2::analysis::femto::utils::getPdgMass(std::abs(TrackSel1.pdgCode.value));
    float mass2 = o2::analysis::femto::utils::getPdgMass(std::abs(TrackSel2.pdgCode.value));

    for (auto const& [p1, p2] : o2::soa::combinations(o2::soa::CombinationsUpperIndexPolicy(tracks, tracks))) {

      bool order1 = std::abs(p1.pdgCode()) == std::abs(TrackSel1.pdgCode.value) &&
                    std::abs(p2.pdgCode()) == std::abs(TrackSel2.pdgCode.value);
      bool order2 = std::abs(p1.pdgCode()) == std::abs(TrackSel2.pdgCode.value) &&
                    std::abs(p2.pdgCode()) == std::abs(TrackSel1.pdgCode.value);

      if (!order1 && !order2) {
        continue;
      }

      ROOT::Math::PtEtaPhiMVector particle1;
      ROOT::Math::PtEtaPhiMVector particle2;

      if (order1) {
        particle1 = ROOT::Math::PtEtaPhiMVector(p1.pt(), p1.eta(), p1.phi(), mass1);
        particle2 = ROOT::Math::PtEtaPhiMVector(p2.pt(), p2.eta(), p2.phi(), mass2);
      } else {
        particle1 = ROOT::Math::PtEtaPhiMVector(p2.pt(), p2.eta(), p2.phi(), mass1);
        particle2 = ROOT::Math::PtEtaPhiMVector(p1.pt(), p1.eta(), p1.phi(), mass2);
      }

      float kstar = getKstar(particle1, particle2);

      // fill for all pairs regardless of kstar cut
      histos.fill(HIST("hTrack1PtMC"), particle1.Pt());
      histos.fill(HIST("hTrack1EtaMC"), particle1.Eta());
      histos.fill(HIST("hTrack1PhiMC"), RecoDecay::constrainAngle(particle1.Phi()));
      histos.fill(HIST("hTrack2PtMC"), particle2.Pt());
      histos.fill(HIST("hTrack2EtaMC"), particle2.Eta());
      histos.fill(HIST("hTrack2PhiMC"), RecoDecay::constrainAngle(particle2.Phi()));
      histos.fill(HIST("hPairKstarMC"), kstar);

      if (kstar > kStarMax.value) {
        continue;
      }
      return true;
    }
    return false;
  }

  template <typename CheckCol, typename CheckTracks>
  bool isEventSelected(CheckCol const& col, CheckTracks const& tracks)
  {
    histos.fill(HIST("EventHist"), 1);
    if (!col.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);
    if (isApplyTFcut && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);
    if (isApplyITSROcut && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);
    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);
    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);
    if (isApplyInelgt0 && !col.isInelGt0()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);
    if (isApplyPairCut && !hasPair(tracks)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);
    return true;
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  template <typename CheckGenTrack>
  bool isGenTrackSelected(CheckGenTrack const& track)
  {
    if (!track.isPhysicalPrimary()) {
      return false;
    }
    if (!track.producedByGenerator()) {
      return false;
    }
    auto pdgTrack = pdg->GetParticle(track.pdgCode());
    if (pdgTrack == nullptr) {
      return false;
    }
    if (std::abs(pdgTrack->Charge()) < KminCharge) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > extraphicut1 && track.phi() < extraphicut2) || track.phi() <= extraphicut3 || track.phi() >= extraphicut4)) {
      return false;
    }
    return true;
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));

  void processData(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), collision.posZ());
    histos.fill(HIST("CentPercentileHist"), collision.centFT0M());
    histos.fill(HIST("hdatazvtxcent"), collision.posZ(), collision.centFT0M());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdcaxy"), track.dcaXY());
      histos.fill(HIST("hdcaz"), track.dcaZ());
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
      histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kGlobalplusITS);
      histos.fill(HIST("hdatadndetaMB"), collision.posZ(), track.eta(), track.phi());
      if (track.hasTPC()) {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kITSonly);
      }
    }
  }

  void processCorrelation(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    if (std::abs(collision.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("VtxZHist"), collision.posZ());

    auto nchTracks = 0;
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      nchTracks++;
    }

    histos.fill(HIST("GlobalMult_vs_FT0A"), nchTracks, collision.multFT0A());
    histos.fill(HIST("GlobalMult_vs_FT0C"), nchTracks, collision.multFT0C());
    histos.fill(HIST("GlobalMult_vs_FV0A"), nchTracks, collision.multFV0A());
    histos.fill(HIST("GlobalMult_vs_NPVtracks"), nchTracks, collision.multNTracksPV());
    histos.fill(HIST("NPVtracks_vs_FT0C"), collision.multNTracksPV(), collision.multFT0C());
  }

  void processStrangeYield(CollisionDataTable::iterator const& collision, V0TrackCandidates const&, aod::V0Datas const& v0data, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    if (std::abs(collision.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("hzvtxcent"), collision.posZ(), collision.centFT0M());
    for (const auto& v0track : v0data) {
      auto v0pTrack = v0track.template posTrack_as<V0TrackCandidates>();
      auto v0nTrack = v0track.template negTrack_as<V0TrackCandidates>();
      if (std::abs(v0pTrack.eta()) <= v0etaK0SCut && std::abs(v0nTrack.eta()) <= v0etaK0SCut && v0pTrack.tpcNClsFound() >= minTPCnClsK0SCut && v0nTrack.tpcNClsFound() >= minTPCnClsK0SCut && std::abs(v0track.dcapostopv()) >= dcapostopvK0SCut && std::abs(v0track.dcanegtopv()) >= dcanegtopvK0SCut && v0track.v0radius() >= v0radiusK0SCut && v0track.v0cosPA() >= v0cospaK0SCut && std::abs(v0track.dcaV0daughters()) <= dcav0daughterK0Scut && std::abs(v0pTrack.tpcNSigmaPi()) <= nSigmaTpcK0SCut && std::abs(v0nTrack.tpcNSigmaPi()) <= nSigmaTpcK0SCut) {

        histos.fill(HIST("K0sCentEtaMass"), collision.centFT0M(), v0track.eta(), v0track.mK0Short());
      }
      if (std::abs(v0pTrack.eta()) <= v0etaLambdaCut && std::abs(v0nTrack.eta()) <= v0etaLambdaCut && v0pTrack.tpcNClsFound() >= minTPCnClsLambdaCut && v0nTrack.tpcNClsFound() >= minTPCnClsLambdaCut && std::abs(v0track.dcapostopv()) >= dcapostopvLambdaCut && std::abs(v0track.dcanegtopv()) >= dcanegtopvLambdaCut && v0track.v0radius() >= v0radiusLambdaCut && v0track.v0cosPA() >= v0cospaLambdaCut && std::abs(v0track.dcaV0daughters()) <= dcav0daughterLambdacut) {

        if (std::abs(v0pTrack.tpcNSigmaPr()) <= nSigmaTpcLambdaCut && std::abs(v0nTrack.tpcNSigmaPi()) <= nSigmaTpcLambdaCut) {
          histos.fill(HIST("LambdaCentEtaMass"), collision.centFT0M(), v0track.eta(), v0track.mLambda());
        }
        if (std::abs(v0pTrack.tpcNSigmaPi()) <= nSigmaTpcLambdaCut && std::abs(v0nTrack.tpcNSigmaPr()) <= nSigmaTpcLambdaCut) {
          histos.fill(HIST("AntiLambdaCentEtaMass"), collision.centFT0M(), v0track.eta(), v0track.mAntiLambda());
        }
      }
    }
  }

  void processMCeff(soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    auto gencent = -999;
    bool atLeastOne = false;

    for (const auto& RecCol : RecCols) {
      auto recTracksPart = RecTracks.sliceByCached(aod::track::collisionId, RecCol.globalIndex(), cache);
      if (!isEventSelected(RecCol, recTracksPart)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
      gencent = RecCol.centFT0M();
    }

    histos.fill(HIST("hGenMCvertexZ"), mcCollision.posZ());
    histos.fill(HIST("hGenMCvtxzcent"), mcCollision.posZ(), gencent);

    if (atLeastOne) {
      histos.fill(HIST("hGenMCAssoRecvertexZ"), mcCollision.posZ());
      histos.fill(HIST("hGenMCAssoRecvtxzcent"), mcCollision.posZ(), gencent);
    }

    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("hGenMCdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        histos.fill(HIST("hGenMCAssoRecdndetaMB"), mcCollision.posZ(), particle.eta(), particle.phi(), static_cast<double>(kGenAll));
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
        }
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kGenPion;
            break;
          case PDG_t::kKPlus:
            pid = kGenKaon;
            break;
          case PDG_t::kProton:
            pid = kGenProton;
            break;
          default:
            pid = kGenOther;
            break;
        }
        histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, particle.eta(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
        histos.fill(HIST("hGenMCAssoRecdndetaMB"), mcCollision.posZ(), particle.eta(), particle.phi(), static_cast<double>(pid));
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {
      auto recTracksPart = RecTracks.sliceByCached(aod::track::collisionId, RecCol.globalIndex(), cache);
      if (!isEventSelected(RecCol, recTracksPart)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      histos.fill(HIST("hRecMCvertexZ"), RecCol.posZ());
      histos.fill(HIST("hRecMCcentrality"), RecCol.centFT0M());
      histos.fill(HIST("hRecMCvtxzcent"), RecCol.posZ(), RecCol.centFT0M());

      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        if (!Rectrack.has_mcParticle()) {
          histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoBkg));
          histos.fill(HIST("hRecMCdndetaMB"), RecCol.posZ(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoBkg));
          continue;
        }
        auto mcpart = Rectrack.mcParticle();
        if (RecCol.mcCollisionId() != mcpart.mcCollisionId()) {
          continue;
        }
        histos.fill(HIST("hmcdcaxy"), Rectrack.dcaXY());
        histos.fill(HIST("hmcdcaz"), Rectrack.dcaZ());
        histos.fill(HIST("hRecMCphivseta"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoAll));
        histos.fill(HIST("hRecMCdndetaMB"), RecCol.posZ(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoAll));

        int pid = 0;
        if (mcpart.isPhysicalPrimary()) {
          switch (std::abs(mcpart.pdgCode())) {
            case PDG_t::kPiPlus:
              pid = kRecoPion;
              break;
            case PDG_t::kKPlus:
              pid = kRecoKaon;
              break;
            case PDG_t::kProton:
              pid = kRecoProton;
              break;
            default:
              pid = kRecoOther;
              break;
          }
        } else {
          pid = kRecoSecondary;
        }
        if (mcpart.has_mothers()) {
          auto mcpartMother = mcpart.template mothers_as<aod::McParticles>().front();
          if (mcpartMother.pdgCode() == PDG_t::kK0Short || std::abs(mcpartMother.pdgCode()) == PDG_t::kLambda0) {
            pid = kRecoWeakDecay;
          }
        }
        if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
          pid = kRecoFake;
        }
        mclabels.push_back(Rectrack.mcParticleId());
        histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), RecCol.centFT0M(), mcpart.eta(), mcpart.phi(), static_cast<double>(pid));
        histos.fill(HIST("hRecMCdndetaMB"), RecCol.posZ(), mcpart.eta(), mcpart.phi(), static_cast<double>(pid));
      } // track (mcrec) loop
    } // collision loop
  }

  void processEvtLossSigLossMC(soa::Join<aod::McCollisions, aod::McCollsExtra, aod::MultMCExtras, aod::McCentFT0Ms>::iterator const& mcCollision, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {

    if ((isApplyInelgt0 && !mcCollision.isInelGt0()) || (isApplyPairCut && !hasPairMC(GenParticles))) {
      return;
    }
    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    histos.fill(HIST("hGenCent"), mcCollision.centFT0M());
    histos.fill(HIST("hMultEta05Gen"), mcCollision.multMCNParticlesEta05());
    bool atLeastOne = false;
    for (const auto& RecCol : RecCols) {
      auto recTracksPart = RecTracks.sliceByCached(aod::track::collisionId, RecCol.globalIndex(), cache);
      if (!isEventSelected(RecCol, recTracksPart)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
    }
    // Generated events with at least one reconstructed collision (event loss estimation)
    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      histos.fill(HIST("hGenAssoRecCent"), mcCollision.centFT0M());
      histos.fill(HIST("hMultEta05GenAssoRec"), mcCollision.multMCNParticlesEta05());
    }
    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      // All generated particles
      histos.fill(HIST("hgendndetaBeforeEvtSel"), particle.eta());
      histos.fill(HIST("hgendndetaVscentBeforeEvtSel"), particle.eta(), mcCollision.centFT0M());
      histos.fill(HIST("hgendndetaVsMultEta05BeforeEvtSel"), particle.eta(), mcCollision.multMCNParticlesEta05());
      if (atLeastOne) {
        // All generated particles with at least one reconstructed collision (signal loss estimation)
        histos.fill(HIST("hgendndetaAfterEvtSel"), particle.eta());
        histos.fill(HIST("hgendndetaVscentAfterEvtSel"), particle.eta(), mcCollision.centFT0M());
        histos.fill(HIST("hgendndetaVsMultEta05AfterEvtSel"), particle.eta(), mcCollision.multMCNParticlesEta05());
      }
    }
  }
  PROCESS_SWITCH(FemtoPairEfficiency, processData, "process data CentFT0C", false);
  PROCESS_SWITCH(FemtoPairEfficiency, processCorrelation, "do correlation study in data/MC", false);
  PROCESS_SWITCH(FemtoPairEfficiency, processStrangeYield, "Strange particle yield", false);
  PROCESS_SWITCH(FemtoPairEfficiency, processMCeff, "process MC efficiency function", false);
  PROCESS_SWITCH(FemtoPairEfficiency, processEvtLossSigLossMC, "process Signal Loss, Event Loss", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FemtoPairEfficiency>(cfgc)};
}
