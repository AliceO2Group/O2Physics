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
///
/// \file heavyionMultiplicity.cxx
///
/// \brief task for analysis of charged-particle multiplicity at midrapidity
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since September 15, 2023

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TPDGCode.h>

#include <cmath>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentFT0CVariant2s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>;
using ColDataTablepp = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Cs, aod::CentFV0As, aod::CentFT0CVariant1s, aod::CentFT0CVariant2s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>>;
using ColMCRecTablepp = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults, aod::CentFT0Ms>>;
using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection>;
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
  kSpeciesbegin = 0,
  kSpPion = 1,
  kSpKaon,
  kSpProton,
  kSpOther,
  kSpStrangeDecay,
  kBkg,
  kSpNotPrimary,
  kSpAll,
  kSpeciesend
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

AxisSpec axisEvent{10, 0.5, 10.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisEtaExtended{100, -5, 5, "#eta", "EtaAxisExtended"};
AxisSpec axisPhi{{0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "#phi", "PhiAxis"};
AxisSpec axisPhi2{629, 0, o2::constants::math::TwoPI, "#phi"};
AxisSpec axisCent{100, 0, 100, "#Cent"};
AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
AxisSpec axisSpecies = {kSpeciesend - 1, +kSpeciesbegin + 0.5, +kSpeciesend - 0.5, "", "SpeciesAxis"};
AxisSpec axisGenTrkType = {kGenTrkTypeend - 1, +kGenTrkTypebegin + 0.5, +kGenTrkTypeend - 0.5, "", "GenTrackTypeAxis"};
AxisSpec axisRecTrkType = {kRecTrkTypeend - 1, +kRecTrkTypebegin + 0.5, +kRecTrkTypeend - 0.5, "", "RecTrackTypeAxis"};
AxisSpec axisMassK0s = {200, 0.4, 0.6, "K0sMass", "K0sMass"};
AxisSpec axisMassLambda = {200, 1.07, 1.17, "Lambda/AntiLamda Mass", "Lambda/AntiLamda Mass"};
AxisSpec axisTracks{9, 0.5, 9.5, "#tracks", "TrackAxis"};
AxisSpec axisDeltaEta{50, -1.0, +1.0, "#Delta(#eta)"};
auto static constexpr KminCharge = 3.f;
auto static constexpr KminPtCut = 0.1f;
auto static constexpr KnItsLayers = 7;

struct HeavyionMultiplicity {

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
  ConfigurableAxis occupancyBin{"occupancyBin", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 5000, 10000}, ""};
  ConfigurableAxis centBinGen{"centBinGen", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 5000, 10000}, ""};
  ConfigurableAxis binsImpactPar{"binsImpactPar", {VARIABLE_WIDTH, 0.0, 3.00065, 4.28798, 6.14552, 7.6196, 8.90942, 10.0897, 11.2002, 12.2709, 13.3167, 14.4173, 23.2518}, "Binning of the impact parameter axis"};
  ConfigurableAxis binsMult{"binsMult", {500, 0.0f, +500.0f}, ""};
  ConfigurableAxis binsDCA{"binsDCA", {500, -10.0f, 10.0f}, ""};

  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyExtraCorrCut{"isApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<float> npvTracksCut{"npvTracksCut", 1.0f, "Apply extra NPVtracks cut"};
  Configurable<float> ft0cCut{"ft0cCut", 1.0f, "Apply extra FT0C cut"};
  Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", true, "Enable NoCollInTimeRangeStandard cut"};
  Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", false, "Enable NoCollInRofStandard cut"};
  Configurable<bool> isApplyNoHighMultCollInPrevRof{"isApplyNoHighMultCollInPrevRof", false, "Enable NoHighMultCollInPrevRof cut"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", false, "Enable FT0CbasedOccupancy cut"};
  Configurable<bool> isApplyCentFT0C{"isApplyCentFT0C", true, "Centrality based on FT0C"};
  Configurable<bool> isApplyCentFV0A{"isApplyCentFV0A", false, "Centrality based on FV0A"};
  Configurable<bool> isApplyCentFT0CVariant1{"isApplyCentFT0CVariant1", false, "Centrality based on FT0C variant1"};
  Configurable<bool> isApplyCentFT0CVariant2{"isApplyCentFT0CVariant2", false, "Centrality based on FT0C variant2 (Run2 like truncation)"};
  Configurable<bool> isApplyCentFT0M{"isApplyCentFT0M", false, "Centrality based on FT0A + FT0C"};
  Configurable<bool> isApplyCentNGlobal{"isApplyCentNGlobal", false, "Centrality based on global tracks"};
  Configurable<bool> isApplyCentMFT{"isApplyCentMFT", false, "Centrality based on MFT tracks"};
  Configurable<bool> isApplySplitRecCol{"isApplySplitRecCol", false, "Split MC reco collisions"};
  Configurable<bool> isApplyInelgt0{"isApplyInelgt0", false, "Enable INEL > 0 condition"};
  Configurable<bool> isApplyTVX{"isApplyTVX", false, "Enable TVX trigger sel"};

  void init(InitContext const&)
  {
    AxisSpec axisMult = {multHistBin, "Mult", "MultAxis"};
    AxisSpec axisPV = {pvHistBin, "PV", "PVAxis"};
    AxisSpec axisFv0aMult = {fv0aMultHistBin, "fv0a", "FV0AMultAxis"};
    AxisSpec axisFt0aMult = {ft0aMultHistBin, "ft0a", "FT0AMultAxis"};
    AxisSpec axisFt0cMult = {ft0cMultHistBin, "ft0c", "FT0CMultAxis"};
    AxisSpec centAxis = {centralityBinning, "Centrality", "CentralityAxis"};
    AxisSpec axisPt = {ptHistBin, "pT", "pTAxis"};
    AxisSpec axisOccupancy = {occupancyBin, "occupancy", "OccupancyAxis"};
    AxisSpec axisCentBinGen = {centBinGen, "GenCentrality", "CentGenAxis"};
    AxisSpec impactParAxis = {binsImpactPar, "Impact Parameter"};
    AxisSpec multAxis = {binsMult, "Multiplicity #eta<0.5"};
    AxisSpec dcaAxis = {binsDCA, "DCA vs PV"};

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoSameBunchPileup"); // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV"); // small difference between z-vertex from PV and from FT0
    x->SetBinLabel(5, "ApplyExtraCorrCut");
    x->SetBinLabel(6, "ApplyNoCollInTimeRangeStandard");
    x->SetBinLabel(7, "ApplyNoCollInRofStandard");
    x->SetBinLabel(8, "ApplyNoHighMultCollInPrevRof");
    x->SetBinLabel(9, "INEL > 0");

    if (doprocessData) {
      histos.add("hdcaxy", "dca to pv in the xy plane", kTH1D, {dcaAxis}, false);
      histos.add("hdcaz", "dca to pv in the z axis", kTH1D, {dcaAxis}, false);
      histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hdatadndeta", "hdatadndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi, axisTrackType}, false);
      histos.add("hdatadndetaMB", "hdatadndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi}, false);
    }

    if (doprocessMonteCarlo || doprocessMCpTefficiency || doprocessMCcheckFakeTracks) {
      histos.add("CentPercentileMCRecHist", "CentPercentileMCRecHist", kTH1D, {axisCent}, false);
      histos.add("hmczvtxcent", "hmczvtxcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
    }

    if (doprocessMonteCarlo) {
      histos.add("hmcdcaxy", "dca to pv in the xy plane", kTH1D, {dcaAxis}, false);
      histos.add("hmcdcaz", "dca to pv in the z axis", kTH1D, {dcaAxis}, false);
      histos.add("MCrecPhiVsEtaHist", "MCrecPhiVsEtaHist", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hmcrecdndeta", "hmcrecdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi, axisSpecies, axisTrackType}, false);
      histos.add("hmcrecdndetaMB", "hmcrecdndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi, axisSpecies}, false);
      histos.add("hmcgendndeta", "hmcgendndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisSpecies, axisGenPtVary}, false);
      histos.add("hmcgendndetaMB", "hmcgendndetaMB", kTHnSparseD, {axisVtxZ, axisEta, axisPhi, axisSpecies}, false);
    }

    if (doprocessMCpTefficiency) {
      histos.add("hmcrecdndpt", "hmcrecdndpt", kTHnSparseD, {centAxis, axisOccupancy, axisTrackType, axisPt}, false);
      histos.add("hmcgendndpt", "hmcgendndpt", kTHnSparseD, {centAxis, axisPt, axisGenPtVary}, false);
    }

    if (doprocessMCcheckFakeTracks) {
      histos.add("hTracksCount", "hTracksCount", kTHnSparseD, {centAxis, axisTracks}, false);
      auto htrack = histos.get<THnSparse>(HIST("hTracksCount"));
      auto* x2 = htrack->GetAxis(1);
      x2->SetBinLabel(1, "All tracks");
      x2->SetBinLabel(2, "Non-fake tracks");
      for (int i = 0; i < KnItsLayers; i++) {
        x2->SetBinLabel(i + 3, Form("layer %d", i));
      }
    }

    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisMult, axisFt0aMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFt0cMult}, true);
      histos.add("Centrality_vs_FT0C", "Centrality_vs_FT0C", kTH2F, {centAxis, axisFt0cMult}, true);
      histos.add("NPVtracks_vs_FT0C", "NPVtracks_vs_FT0C", kTH2F, {axisPV, axisFt0cMult}, true);
      histos.add("GlobalMult_vs_FV0A", "GlobalMult_vs_FV0A", kTH2F, {axisMult, axisFv0aMult}, true);
      histos.add("Centrality_vs_FV0A", "Centrality_vs_FV0A", kTH2F, {centAxis, axisFv0aMult}, true);
      histos.add("CentFT0Ccentrality_vs_GlobalMult", "CentFT0Ccentrality_vs_GlobalMult", kTH2F, {centAxis, axisMult}, true);
      histos.add("NPVtracks_vs_GlobalMult", "NPVtracks_vs_GlobalMult", kTH2F, {axisPV, axisMult}, true);
    }

    if (doprocessStrangeYield) {
      histos.add("hzvtxcent", "hzvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("K0sCentEtaMass", "K0sCentEtaMass", kTH3D, {centAxis, axisEta, axisMassK0s}, false);
      histos.add("LambdaCentEtaMass", "LambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
      histos.add("AntiLambdaCentEtaMass", "AntiLambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
    }

    if (doprocessppData) {
      histos.add("MultPercentileHist", "MultPercentileHist", kTH1D, {axisCent}, false);
      histos.add("hdatazvtxmultpp", "hdatazvtxmultpp", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("PhiVsEtaHistpp", "PhiVsEtaHistpp", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hdatadndetapp", "hdatadndetapp", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisTrackType}, false);
    }

    if (doprocessppMonteCarlo) {
      histos.add("MultPercentileMCRecHist", "MultPercentileMCRecHist", kTH1D, {axisCent}, false);
      histos.add("hmczvtxmultpp", "hmczvtxmultpp", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("MCrecPhiVsEtaHistpp", "MCrecPhiVsEtaHistpp", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hmcrecdndetapp", "hmcrecdndetapp", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisSpecies, axisTrackType}, false);
      histos.add("hmcgendndetapp", "hmcgendndetapp", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisSpecies, axisGenPtVary}, false);
    }

    if (doprocessGen) {
      histos.add("MultBarrelEta10_vs_FT0A", "MultBarrelEta10_vs_FT0A", kTH2F, {axisMult, axisFt0aMult}, true);
      histos.add("MultBarrelEta10_vs_FT0C", "MultBarrelEta10_vs_FT0C", kTH2F, {axisMult, axisFt0cMult}, true);
      histos.add("MultBarrelEta10", "MultBarrelEta10", kTH1F, {axisMult}, true);
      histos.add("MultFT0A", "MultFT0A", kTH1F, {axisFt0aMult}, true);
      histos.add("MultFT0C", "MultFT0C", kTH1F, {axisFt0cMult}, true);
      histos.add("dndeta10_vs_FT0C", "dndeta10_vs_FT0C", kTH2F, {axisEta, axisCentBinGen}, true);
      histos.add("dndeta10_vs_FT0A", "dndeta10_vs_FT0A", kTH2F, {axisEta, axisCentBinGen}, true);
      histos.add("mult10_vs_FT0C", "mult10_vs_FT0C", kTH2F, {axisMult, axisCentBinGen}, true);
      histos.add("mult10_vs_FT0A", "mult10_vs_FT0A", kTH2F, {axisMult, axisCentBinGen}, true);
    }

    if (doprocessEvtLossSigLossMC) {
      histos.add("MCEventHist", "MCEventHist", kTH1F, {axisEvent}, false);
      auto hstat = histos.get<TH1>(HIST("MCEventHist"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All MC events");
      x->SetBinLabel(2, "MC events with reco event after event selection");
      x->SetBinLabel(3, "MC events with no reco events");
      histos.add("hgendndetaVscentGenwithNOreco", "dndeta vs impact parameter, gen events with no reco", kTH2F, {axisEtaExtended, impactParAxis});
      histos.add("hgendndetaVscentGenwithReco", "dndeta vs impact parameter, gen events with at least one reco", kTH2F, {axisEtaExtended, impactParAxis});

      histos.add("hMultEta05GenwithNoreco", "multiplicity in eta<0.5 of generated MC events, with no recoevent", kTH1F, {multAxis});
      histos.add("hMultEta05Gen", "multiplicity in eta<0.5 of generated MC events", kTH1F, {multAxis});
      histos.add("hMultEta05Rec", "multiplicity in eta<0.5 of selected MC events", kTH1F, {multAxis});
      histos.add("hMultEta05vsCentrRec", "multiplicity in eta<0.5 of selected MC events vs centrality", kTH2F, {axisCent, multAxis});
      histos.add("hgendndetaVsMultEta05BeforeEvtSel", "hgendndetaBeforeEvtSel vs multiplicity in eta<0.5", kTH2F, {axisEta, multAxis});
      histos.add("hgendndetaVsMultEta05AfterEvtSel", "hgendndetaAfterEvtSel vs multiplicity in eta<0.5", kTH2F, {axisEta, multAxis});

      histos.add("hMultGen", "multiplicity of generated MC events", kTH1F, {axisFt0cMult});
      histos.add("hMultRec", "multiplicity of selected MC events", kTH1F, {axisFt0cMult});
      histos.add("hMultvsCentrRec", "multiplicity of selected MC events vs centrality", kTH2F, {axisCent, axisFt0cMult});
      histos.add("hgendndetaVsMultBeforeEvtSel", "hgendndetaBeforeEvtSel vs multiplicity", kTH2F, {axisEta, axisFt0cMult});
      histos.add("hgendndetaVsMultAfterEvtSel", "hgendndetaAfterEvtSel vs multiplicity", kTH2F, {axisEta, axisFt0cMult});

      histos.add("hImpactParameterGenwithNoreco", "Impact parameter of generated MC events, with no recoevent", kTH1F, {impactParAxis});
      histos.add("hImpactParameterGen", "Impact parameter of generated MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParameterRec", "Impact parameter of selected MC events", kTH1F, {impactParAxis});
      histos.add("hImpactParvsCentrRec", "Impact parameter of selected MC events vs centrality", kTH2F, {axisCent, impactParAxis});
      histos.add("hgendndetaBeforeEvtSel", "Eta of all generated particles", kTH1F, {axisEta});
      histos.add("hgendndetaAfterEvtSel", "Eta of generated particles after EvtSel", kTH1F, {axisEta});
      histos.add("hgendndetaVscentBeforeEvtSel", "hgendndetaBeforeEvtSel vs centrality", kTH2F, {axisEta, impactParAxis});
      histos.add("hgendndetaVscentAfterEvtSel", "hgendndetaAfterEvtSel vs centrality", kTH2F, {axisEta, impactParAxis});
    }

    if (doprocessMCeff) {
      histos.add("hGenMCvertexZ", "hGenMCvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hGenMCvtxzcent", "hGenMCvtxzcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
      histos.add("hGenMCAssoRecvertexZ", "hGenMCAssoRecvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hGenMCAssoRecvtxzcent", "hGenMCAssoRecvtxzcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
      histos.add("hGenMCdndeta", "hGenMCdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi}, false);
      histos.add("hGenMCAssoRecdndeta", "hGenMCAssoRecdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi, axisGenTrkType, axisGenPtVary}, false);

      histos.add("hRecMCvertexZ", "hRecMCvertexZ", kTH1D, {axisVtxZ}, false);
      histos.add("hRecMCvtxzcent", "hRecMCvtxzcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
      histos.add("hRecMCcentrality", "hRecMCcentrality", kTH1D, {axisCent}, false);
      histos.add("hRecMCphivseta", "hRecMCphivseta", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hRecMCdndeta", "hRecMCdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi, axisRecTrkType}, false);
      histos.add("etaResolution", "etaResolution", kTH2D, {axisEta, axisDeltaEta});
    }
  }

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);

    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);

    if (isApplyTVX && !col.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }

    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);

    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);

    if (isApplyExtraCorrCut && col.multNTracksPV() > npvTracksCut && col.multFT0C() < (10 * col.multNTracksPV() - ft0cCut)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);

    if (isApplyNoCollInTimeRangeStandard && !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);

    if (isApplyNoCollInRofStandard && !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);

    if (isApplyNoHighMultCollInPrevRof && !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);

    if (isApplyInelgt0 && !col.isInelGt0()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);
    return true;
  }

  template <typename CheckColCent>
  float selColCent(CheckColCent const& col)
  {
    auto cent = -1;
    if (isApplyCentFT0C) {
      cent = col.centFT0C();
    }
    if (isApplyCentFV0A) {
      cent = col.centFV0A();
    }
    if (isApplyCentFT0CVariant1) {
      cent = col.centFT0CVariant1();
    }
    if (isApplyCentFT0CVariant2) {
      cent = col.centFT0CVariant2();
    }
    if (isApplyCentFT0M) {
      cent = col.centFT0M();
    }
    if (isApplyCentNGlobal) {
      cent = col.centNGlobal();
    }
    if (isApplyCentMFT) {
      cent = col.centMFT();
    }
    return cent;
  }

  template <typename CheckColCent>
  float selColMultMC(CheckColCent const& col)
  {
    auto cent = -1;
    if (isApplyCentFT0C) {
      cent = col.multMCFT0C();
    }
    if (isApplyCentFV0A) {
      cent = col.multMCFV0A();
    }
    return cent;
  }

  template <typename CheckColOccu>
  float selColOccu(CheckColOccu const& col)
  {
    auto occu = isApplyFT0CbasedOccupancy ? col.ft0cOccupancyInTimeRange() : col.trackOccupancyInTimeRange();
    return occu;
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

  expressions::Filter trackSelectionProperMixed = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                                  ncheckbit(aod::track::trackCutFlag, TrackSelectionIts) &&
                                                  ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                         ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true) &&
                                                  ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                                         ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));

  void processData(CollisionDataTable::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());
    histos.fill(HIST("CentPercentileHist"), selColCent(cols));
    histos.fill(HIST("hdatazvtxcent"), cols.posZ(), selColCent(cols), selColOccu(cols));

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("hdcaxy"), track.dcaXY());
      histos.fill(HIST("hdcaz"), track.dcaZ());
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
      histos.fill(HIST("hdatadndeta"), cols.posZ(), selColCent(cols), selColOccu(cols), track.eta(), track.phi(), kGlobalplusITS);
      histos.fill(HIST("hdatadndetaMB"), cols.posZ(), track.eta(), track.phi());

      if (track.hasTPC()) {
        histos.fill(HIST("hdatadndeta"), cols.posZ(), selColCent(cols), selColOccu(cols), track.eta(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatadndeta"), cols.posZ(), selColCent(cols), selColOccu(cols), track.eta(), track.phi(), kITSonly);
      }
    }
  }

  void processCorrelation(CollisionDataTable::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    if (std::abs(cols.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("VtxZHist"), cols.posZ());

    auto nchTracks = 0;
    for (const auto& track : tracks) {
      if (std::abs(track.eta()) >= etaRange) {
        continue;
      }
      nchTracks++;
    }

    histos.fill(HIST("GlobalMult_vs_FT0A"), nchTracks, cols.multFT0A());
    histos.fill(HIST("GlobalMult_vs_FT0C"), nchTracks, cols.multFT0C());
    histos.fill(HIST("Centrality_vs_FT0C"), cols.centFT0C(), cols.multFT0C());
    histos.fill(HIST("NPVtracks_vs_FT0C"), cols.multNTracksPV(), cols.multFT0C());
    histos.fill(HIST("GlobalMult_vs_FV0A"), nchTracks, cols.multFV0A());
    histos.fill(HIST("Centrality_vs_FV0A"), cols.centFV0A(), cols.multFV0A());
    histos.fill(HIST("CentFT0Ccentrality_vs_GlobalMult"), cols.centFT0C(), nchTracks);
    histos.fill(HIST("NPVtracks_vs_GlobalMult"), cols.multNTracksPV(), nchTracks);
  }

  void processMonteCarlo(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {

    if (isApplySplitRecCol && (RecCols.size() == 0 || RecCols.size() > 1)) {
      return;
    }

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCol.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selColCent(RecCol));
      histos.fill(HIST("hmczvtxcent"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol));
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("hmcdcaxy"), Rectrack.dcaXY());
        histos.fill(HIST("hmcdcaz"), Rectrack.dcaZ());
        histos.fill(HIST("MCrecPhiVsEtaHist"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("hmcrecdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kGlobalplusITS);
        histos.fill(HIST("hmcrecdndetaMB"), RecCol.posZ(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll));
        if (Rectrack.hasTPC()) {
          histos.fill(HIST("hmcrecdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kGlobalonly);
        } else {
          histos.fill(HIST("hmcrecdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kITSonly);
        }

        if (Rectrack.has_mcParticle()) {
          int pid = kBkg;
          auto mcpart = Rectrack.template mcParticle_as<aod::McParticles>();
          if (mcpart.isPhysicalPrimary()) {
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kSpPion;
                break;
              case PDG_t::kKPlus:
                pid = kSpKaon;
                break;
              case PDG_t::kProton:
                pid = kSpProton;
                break;
              default:
                pid = kSpOther;
                break;
            }
          } else {
            pid = kSpNotPrimary;
          }
          if (mcpart.has_mothers()) {
            auto mcpartMother = mcpart.template mothers_as<aod::McParticles>().front();
            if (mcpartMother.pdgCode() == PDG_t::kK0Short || std::abs(mcpartMother.pdgCode()) == PDG_t::kLambda0) {
              pid = kSpStrangeDecay;
            }
          }
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kBkg;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("hmcrecdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(pid), kGlobalplusITS);
          histos.fill(HIST("hmcrecdndetaMB"), RecCol.posZ(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(pid));
        } else {
          histos.fill(HIST("hmcrecdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kBkg), kGlobalplusITS);
          histos.fill(HIST("hmcrecdndetaMB"), RecCol.posZ(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kBkg));
        }
      } // track (mcrec) loop

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kNoGenpTVar);
        histos.fill(HIST("hmcgendndetaMB"), RecCol.posZ(), particle.eta(), particle.phi(), static_cast<double>(kSpAll));
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTup);
          histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTdown);
        }

        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kSpPion;
            break;
          case PDG_t::kKPlus:
            pid = kSpKaon;
            break;
          case PDG_t::kProton:
            pid = kSpProton;
            break;
          default:
            pid = kSpOther;
            break;
        }
        histos.fill(HIST("hmcgendndeta"), RecCol.posZ(), selColCent(RecCol), particle.eta(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
        histos.fill(HIST("hmcgendndetaMB"), RecCol.posZ(), particle.eta(), particle.phi(), static_cast<double>(pid));
      } // track (mcgen) loop
    } // collision loop
  }

  void processMCpTefficiency(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (std::abs(RecCol.posZ()) >= vtxRange) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCol.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selColCent(RecCol));
      histos.fill(HIST("hmczvtxcent"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol));

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (Rectrack.has_mcParticle()) {
          auto mcpart = Rectrack.mcParticle();
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hmcrecdndpt"), selColCent(RecCol), selColOccu(RecCol), kGlobalplusITS, mcpart.pt());
            if (Rectrack.hasTPC()) {
              histos.fill(HIST("hmcrecdndpt"), selColCent(RecCol), selColOccu(RecCol), kGlobalonly, mcpart.pt());
            } else {
              histos.fill(HIST("hmcrecdndpt"), selColCent(RecCol), selColOccu(RecCol), kITSonly, mcpart.pt());
            }
          }
        }
      }

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("hmcgendndpt"), selColCent(RecCol), particle.pt(), kNoGenpTVar);
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hmcgendndpt"), selColCent(RecCol), particle.pt(), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndpt"), selColCent(RecCol), particle.pt(), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndpt"), selColCent(RecCol), particle.pt(), kGenpTup);
          histos.fill(HIST("hmcgendndpt"), selColCent(RecCol), particle.pt(), kGenpTdown);
        }
      }
    }
  }

  void processMCcheckFakeTracks(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCols, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (std::abs(RecCol.posZ()) >= vtxRange) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCol.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selColCent(RecCol));
      histos.fill(HIST("hmczvtxcent"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol));

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (!Rectrack.hasTPC()) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), selColCent(RecCol), 1);
        bool isFakeItsTracks = false;
        for (int i = 0; i < KnItsLayers; i++) {
          if (Rectrack.mcMask() & 1 << i) {
            isFakeItsTracks = true;
            histos.fill(HIST("hTracksCount"), selColCent(RecCol), i + 3);
            break;
          }
        }
        if (isFakeItsTracks) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), selColCent(RecCol), 2);
      }
    }
  }

  void processStrangeYield(CollisionDataTable::iterator const& cols, V0TrackCandidates const&, aod::V0Datas const& v0data)
  {
    if (!isEventSelected(cols)) {
      return;
    }
    if (std::abs(cols.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("hzvtxcent"), cols.posZ(), selColCent(cols));
    for (const auto& v0track : v0data) {
      auto v0pTrack = v0track.template posTrack_as<V0TrackCandidates>();
      auto v0nTrack = v0track.template negTrack_as<V0TrackCandidates>();
      if (std::abs(v0pTrack.eta()) <= v0etaK0SCut && std::abs(v0nTrack.eta()) <= v0etaK0SCut && v0pTrack.tpcNClsFound() >= minTPCnClsK0SCut && v0nTrack.tpcNClsFound() >= minTPCnClsK0SCut && std::abs(v0track.dcapostopv()) >= dcapostopvK0SCut && std::abs(v0track.dcanegtopv()) >= dcanegtopvK0SCut && v0track.v0radius() >= v0radiusK0SCut && v0track.v0cosPA() >= v0cospaK0SCut && std::abs(v0track.dcaV0daughters()) <= dcav0daughterK0Scut && std::abs(v0pTrack.tpcNSigmaPi()) <= nSigmaTpcK0SCut && std::abs(v0nTrack.tpcNSigmaPi()) <= nSigmaTpcK0SCut) {

        histos.fill(HIST("K0sCentEtaMass"), selColCent(cols), v0track.eta(), v0track.mK0Short());
      }
      if (std::abs(v0pTrack.eta()) <= v0etaLambdaCut && std::abs(v0nTrack.eta()) <= v0etaLambdaCut && v0pTrack.tpcNClsFound() >= minTPCnClsLambdaCut && v0nTrack.tpcNClsFound() >= minTPCnClsLambdaCut && std::abs(v0track.dcapostopv()) >= dcapostopvLambdaCut && std::abs(v0track.dcanegtopv()) >= dcanegtopvLambdaCut && v0track.v0radius() >= v0radiusLambdaCut && v0track.v0cosPA() >= v0cospaLambdaCut && std::abs(v0track.dcaV0daughters()) <= dcav0daughterLambdacut) {

        if (std::abs(v0pTrack.tpcNSigmaPr()) <= nSigmaTpcLambdaCut && std::abs(v0nTrack.tpcNSigmaPi()) <= nSigmaTpcLambdaCut) {
          histos.fill(HIST("LambdaCentEtaMass"), selColCent(cols), v0track.eta(), v0track.mLambda());
        }
        if (std::abs(v0pTrack.tpcNSigmaPi()) <= nSigmaTpcLambdaCut && std::abs(v0nTrack.tpcNSigmaPr()) <= nSigmaTpcLambdaCut) {
          histos.fill(HIST("AntiLambdaCentEtaMass"), selColCent(cols), v0track.eta(), v0track.mAntiLambda());
        }
      }
    }
  }

  void processppData(ColDataTablepp::iterator const& cols, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(cols)) {
      return;
    }

    histos.fill(HIST("VtxZHist"), cols.posZ());
    histos.fill(HIST("MultPercentileHist"), cols.centFT0M());
    histos.fill(HIST("hdatazvtxmultpp"), cols.posZ(), cols.centFT0M());

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("PhiVsEtaHistpp"), track.phi(), track.eta());
      histos.fill(HIST("hdatadndetapp"), cols.posZ(), cols.centFT0M(), track.eta(), track.phi(), kGlobalplusITS);
      if (track.hasTPC()) {
        histos.fill(HIST("hdatadndetapp"), cols.posZ(), cols.centFT0M(), track.eta(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatadndetapp"), cols.posZ(), cols.centFT0M(), track.eta(), track.phi(), kITSonly);
      }
    } // track loop
  }

  void processppMonteCarlo(CollisionMCTrueTable::iterator const&, ColMCRecTablepp const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    if (isApplySplitRecCol && (RecCols.size() == 0 || RecCols.size() > 1)) {
      return;
    }

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;

      histos.fill(HIST("VtxZHist"), RecCol.posZ());
      histos.fill(HIST("MultPercentileMCRecHist"), RecCol.centFT0M());
      histos.fill(HIST("hmczvtxmultpp"), RecCol.posZ(), RecCol.centFT0M());

      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("MCrecPhiVsEtaHistpp"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("hmcrecdndetapp"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kGlobalplusITS);
        if (Rectrack.hasTPC()) {
          histos.fill(HIST("hmcrecdndetapp"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kGlobalonly);
        } else {
          histos.fill(HIST("hmcrecdndetapp"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kSpAll), kITSonly);
        }

        if (Rectrack.has_mcParticle()) {
          int pid = kBkg;
          auto mcpart = Rectrack.template mcParticle_as<aod::McParticles>();
          if (mcpart.isPhysicalPrimary()) {
            switch (std::abs(mcpart.pdgCode())) {
              case PDG_t::kPiPlus:
                pid = kSpPion;
                break;
              case PDG_t::kKPlus:
                pid = kSpKaon;
                break;
              case PDG_t::kProton:
                pid = kSpProton;
                break;
              default:
                pid = kSpOther;
                break;
            }
          } else {
            pid = kSpNotPrimary;
          }
          if (mcpart.has_mothers()) {
            auto mcpartMother = mcpart.template mothers_as<aod::McParticles>().front();
            if (mcpartMother.pdgCode() == PDG_t::kK0Short || std::abs(mcpartMother.pdgCode()) == PDG_t::kLambda0) {
              pid = kSpStrangeDecay;
            }
          }
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kBkg;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("hmcrecdndetapp"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(pid), kGlobalplusITS);
        } else {
          histos.fill(HIST("hmcrecdndetapp"), RecCol.posZ(), RecCol.centFT0M(), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kBkg), kGlobalplusITS);
        }
      } // track (mcrec) loop

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kNoGenpTVar);
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTup);
          histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(kSpAll), kGenpTdown);
        }

        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case PDG_t::kPiPlus:
            pid = kSpPion;
            break;
          case PDG_t::kKPlus:
            pid = kSpKaon;
            break;
          case PDG_t::kProton:
            pid = kSpProton;
            break;
          default:
            pid = kSpOther;
            break;
        }
        histos.fill(HIST("hmcgendndetapp"), RecCol.posZ(), RecCol.centFT0M(), particle.eta(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // track (mcgen) loop
    } // collision loop
  }

  void processGen(aod::McCollisions::iterator const&, aod::McParticles const& GenParticles)
  {

    int multFT0A = 0;
    int multFT0C = 0;
    int multBarrelEta10 = 0;

    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      if (std::abs(particle.eta()) < 1.0) {
        multBarrelEta10++;
      }
      if (-3.3 < particle.eta() && particle.eta() < -2.1) {
        multFT0C++;
      }
      if (3.5 < particle.eta() && particle.eta() < 4.9) {
        multFT0A++;
      }
    }

    histos.fill(HIST("MultBarrelEta10_vs_FT0A"), multBarrelEta10, multFT0A);
    histos.fill(HIST("MultBarrelEta10_vs_FT0C"), multBarrelEta10, multFT0C);
    histos.fill(HIST("MultBarrelEta10"), multBarrelEta10);
    histos.fill(HIST("MultFT0A"), multFT0A);
    histos.fill(HIST("MultFT0C"), multFT0C);
    histos.fill(HIST("mult10_vs_FT0A"), multBarrelEta10, multFT0A);
    histos.fill(HIST("mult10_vs_FT0C"), multBarrelEta10, multFT0C);

    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("dndeta10_vs_FT0A"), particle.eta(), multFT0A);
      histos.fill(HIST("dndeta10_vs_FT0C"), particle.eta(), multFT0C);
    }
  }

  void processEvtLossSigLossMC(soa::Join<CollisionMCTrueTable, aod::MultMCExtras>::iterator const& mcCollision, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles)
  {
    if (isApplyInelgt0 && !mcCollision.isInelGt0()) {
      return;
    }

    if (isApplyTVX && !(mcCollision.multMCFT0C() > 0 && mcCollision.multMCFT0A() > 0)) {
      return;
    }

    if (std::abs(mcCollision.posZ()) >= vtxRange) {
      return;
    }
    // All generated events
    histos.fill(HIST("MCEventHist"), 1);
    histos.fill(HIST("hImpactParameterGen"), mcCollision.impactParameter());
    histos.fill(HIST("hMultEta05Gen"), mcCollision.multMCNParticlesEta05());
    histos.fill(HIST("hMultGen"), selColMultMC(mcCollision));

    if (RecCols.size() == 0) {
      histos.fill(HIST("MCEventHist"), 3);
      histos.fill(HIST("hImpactParameterGenwithNoreco"), mcCollision.impactParameter());
      histos.fill(HIST("hMultEta05GenwithNoreco"), mcCollision.multMCNParticlesEta05());
    }

    bool atLeastOne = false;
    auto centrality = -999.;
    auto numcontributors = -999;
    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (std::abs(RecCol.posZ()) >= vtxRange) {
        continue;
      }
      if (RecCol.numContrib() <= numcontributors) {
        continue;
      } else {
        numcontributors = RecCol.numContrib();
      }
      centrality = selColCent(RecCol);
      atLeastOne = true;
    }

    // Generated events with at least one reconstructed collision (event loss estimation)
    if (atLeastOne) {
      histos.fill(HIST("MCEventHist"), 2);
      histos.fill(HIST("hImpactParameterRec"), mcCollision.impactParameter());
      histos.fill(HIST("hMultEta05Rec"), mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("hMultRec"), selColMultMC(mcCollision));
      histos.fill(HIST("hImpactParvsCentrRec"), centrality, mcCollision.impactParameter());
      histos.fill(HIST("hMultEta05vsCentrRec"), centrality, mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("hMultvsCentrRec"), centrality, selColMultMC(mcCollision));
    }

    for (const auto& particle : GenParticles) {

      if (RecCols.size() == 0) {
        histos.fill(HIST("hgendndetaVscentGenwithNOreco"), particle.eta(), mcCollision.impactParameter());
      } else {
        histos.fill(HIST("hgendndetaVscentGenwithReco"), particle.eta(), mcCollision.impactParameter());
      }

      if (!isGenTrackSelected(particle)) {
        continue;
      }

      // All generated particles
      histos.fill(HIST("hgendndetaBeforeEvtSel"), particle.eta());
      histos.fill(HIST("hgendndetaVscentBeforeEvtSel"), particle.eta(), mcCollision.impactParameter());
      histos.fill(HIST("hgendndetaVsMultEta05BeforeEvtSel"), particle.eta(), mcCollision.multMCNParticlesEta05());
      histos.fill(HIST("hgendndetaVsMultBeforeEvtSel"), particle.eta(), selColMultMC(mcCollision));

      if (atLeastOne) {
        // All generated particles with at least one reconstructed collision (signal loss estimation)
        histos.fill(HIST("hgendndetaAfterEvtSel"), particle.eta());
        histos.fill(HIST("hgendndetaVscentAfterEvtSel"), particle.eta(), mcCollision.impactParameter());
        histos.fill(HIST("hgendndetaVsMultEta05AfterEvtSel"), particle.eta(), mcCollision.multMCNParticlesEta05());
        histos.fill(HIST("hgendndetaVsMultAfterEvtSel"), particle.eta(), selColMultMC(mcCollision));
      }
    }
  }

  void processMCeff(soa::Join<aod::McCollisions, aod::McCollsExtra>::iterator const& mcCollision, CollisionMCRecTable const& RecCols, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    auto gencent = -999;
    auto genoccu = -999;
    bool atLeastOne = false;

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      atLeastOne = true;
      gencent = selColCent(RecCol);
      genoccu = selColOccu(RecCol);
    }

    histos.fill(HIST("hGenMCvertexZ"), mcCollision.posZ());
    histos.fill(HIST("hGenMCvtxzcent"), mcCollision.posZ(), gencent, genoccu);

    if (atLeastOne) {
      histos.fill(HIST("hGenMCAssoRecvertexZ"), mcCollision.posZ());
      histos.fill(HIST("hGenMCAssoRecvtxzcent"), mcCollision.posZ(), gencent, genoccu);
    }

    for (const auto& particle : GenParticles) {
      if (!isGenTrackSelected(particle)) {
        continue;
      }
      histos.fill(HIST("hGenMCdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi());
      if (atLeastOne) {
        histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kNoGenpTVar);
        if (particle.pt() < KminPtCut) {
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTup);
          histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(kGenAll), kGenpTdown);
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
        histos.fill(HIST("hGenMCAssoRecdndeta"), mcCollision.posZ(), gencent, genoccu, particle.eta(), particle.phi(), static_cast<double>(pid), kNoGenpTVar);
      } // Associated with reco col
    } // track (mcgen) loop

    for (const auto& RecCol : RecCols) {
      if (!isEventSelected(RecCol)) {
        continue;
      }
      if (RecCol.globalIndex() != mcCollision.bestCollisionIndex()) {
        continue;
      }
      histos.fill(HIST("hRecMCvertexZ"), RecCol.posZ());
      histos.fill(HIST("hRecMCcentrality"), selColCent(RecCol));
      histos.fill(HIST("hRecMCvtxzcent"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol));

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCol.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("hRecMCphivseta"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoAll));
        if (Rectrack.has_mcParticle()) {
          int pid = 0;
          auto mcpart = Rectrack.mcParticle();
          histos.fill(HIST("etaResolution"), Rectrack.eta(), Rectrack.eta() - mcpart.eta());
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
          histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), mcpart.eta(), mcpart.phi(), static_cast<double>(pid));
        } else {
          histos.fill(HIST("hRecMCdndeta"), RecCol.posZ(), selColCent(RecCol), selColOccu(RecCol), Rectrack.eta(), Rectrack.phi(), static_cast<double>(kRecoBkg));
        }
      } // track (mcrec) loop
    } // collision loop
  }

  PROCESS_SWITCH(HeavyionMultiplicity, processData, "process data CentFT0C", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processCorrelation, "do correlation study in data", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processMonteCarlo, "process MC CentFT0C", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processMCpTefficiency, "process MC pTefficiency", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processMCcheckFakeTracks, "Check Fake tracks", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processStrangeYield, "Strange particle yield", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processppData, "process pp data", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processppMonteCarlo, "process pp MC", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processGen, "process pure MC gen", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processEvtLossSigLossMC, "process Signal Loss, Event Loss", false);
  PROCESS_SWITCH(HeavyionMultiplicity, processMCeff, "process extra efficiency function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyionMultiplicity>(cfgc)};
}
