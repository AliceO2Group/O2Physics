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

#include <cmath>
#include <cstdlib>
#include <TPDGCode.h>
#include <vector>

#include "PWGMM/Mult/DataModel/bestCollisionTable.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"
#include "PWGMM/Mult/DataModel/Index.h"
#include "Common/DataModel/PIDResponse.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

using CollisionDataTable = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>;
using TrackDataTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using FilTrackDataTable = soa::Filtered<TrackDataTable>;
using CollisionMCTrueTable = aod::McCollisions;
using TrackMCTrueTable = aod::McParticles;
using CollisionMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs, aod::CentFT0CVariant1s, aod::CentFT0Ms, aod::CentNGlobals, aod::CentMFTs>>;
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

AxisSpec axisEvent{12, 0.5, 12.5, "#Event", "EventAxis"};
AxisSpec axisVtxZ{40, -20, 20, "Vertex Z", "VzAxis"};
AxisSpec axisEta{40, -2, 2, "#eta", "EtaAxis"};
AxisSpec axisPhi{{0, o2::constants::math::PIQuarter, o2::constants::math::PIHalf, o2::constants::math::PIQuarter * 3., o2::constants::math::PI, o2::constants::math::PIQuarter * 5., o2::constants::math::PIHalf * 3., o2::constants::math::PIQuarter * 7., o2::constants::math::TwoPI}, "#phi", "PhiAxis"};
AxisSpec axisPhi2{629, 0, o2::constants::math::TwoPI, "#phi"};
AxisSpec axisCent{100, 0, 100, "#Cent"};
AxisSpec axisTrackType = {kTrackTypeend - 1, +kTrackTypebegin + 0.5, +kTrackTypeend - 0.5, "", "TrackTypeAxis"};
AxisSpec axisGenPtVary = {kGenpTend - 1, +kGenpTbegin + 0.5, +kGenpTend - 0.5, "", "GenpTVaryAxis"};
AxisSpec axisSpecies = {kSpeciesend - 1, +kSpeciesbegin + 0.5, +kSpeciesend - 0.5, "", "SpeciesAxis"};
AxisSpec axisMassK0s = {200, 0.4, 0.6, "K0sMass", "K0sMass"};
AxisSpec axisMassLambda = {200, 1.07, 1.17, "Lambda/AntiLamda Mass", "Lambda/AntiLamda Mass"};
AxisSpec axisTracks{9, 0.5, 9.5, "#tracks", "TrackAxis"};

struct HeavyionMultiplicity {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  Configurable<float> etaRange{"etaRange", 1.0f, "Eta range to consider"};
  Configurable<float> vtxRange{"vtxRange", 10.0f, "Vertex Z range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> v0radiusCut{"v0radiusCut", 1.2f, "RadiusCut"};
  Configurable<float> dcapostopvCut{"dcapostopvCut", 0.05f, "dcapostopvCut"};
  Configurable<float> dcanegtopvCut{"dcanegtopvCut", 0.05f, "dcanegtopvCut"};
  Configurable<float> v0cospaCut{"v0cospaCut", 0.995f, "v0cospaCut"};
  Configurable<float> dcav0daughtercut{"dcav0daughtercut", 1.0f, "dcav0daughtercut"};
  Configurable<float> minTPCnClsCut{"minTPCnClsCut", 50.0f, "minTPCnClsCut"};
  Configurable<float> nSigmaTpcCut{"nSigmaTpcCut", 5.0f, "nSigmaTpcCut"};
  ConfigurableAxis multHistBin{"multHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis pvHistBin{"pvHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis fv0aMultHistBin{"fv0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0aMultHistBin{"ft0aMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ft0cMultHistBin{"ft0cMultHistBin", {501, -0.5, 500.5}, ""};
  ConfigurableAxis ptHistBin{"ptHistBin", {200, 0., 20.}, ""};
  ConfigurableAxis centralityBinning{"centralityBinning", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, ""};
  ConfigurableAxis occupancyBin{"occupancyBin", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 5000, 10000}, ""};

  Configurable<bool> isApplySameBunchPileup{"isApplySameBunchPileup", true, "Enable SameBunchPileup cut"};
  Configurable<bool> isApplyGoodZvtxFT0vsPV{"isApplyGoodZvtxFT0vsPV", true, "Enable GoodZvtxFT0vsPV cut"};
  Configurable<bool> isApplyVertexITSTPC{"isApplyVertexITSTPC", true, "Enable VertexITSTPC cut"};
  Configurable<bool> isApplyVertexTOFmatched{"isApplyVertexTOFmatched", true, "Enable VertexTOFmatched cut"};
  Configurable<bool> isApplyVertexTRDmatched{"isApplyVertexTRDmatched", true, "Enable VertexTRDmatched cut"};
  Configurable<bool> isApplyExtraCorrCut{"isApplyExtraCorrCut", false, "Enable extra NPVtracks vs FTOC correlation cut"};
  Configurable<bool> isApplyExtraPhiCut{"isApplyExtraPhiCut", false, "Enable extra phi cut"};
  Configurable<float> npvTracksCut{"npvTracksCut", 1.0f, "Apply extra NPVtracks cut"};
  Configurable<float> ft0cCut{"ft0cCut", 1.0f, "Apply extra FT0C cut"};
  Configurable<bool> isApplyNoCollInTimeRangeStandard{"isApplyNoCollInTimeRangeStandard", true, "Enable NoCollInTimeRangeStandard cut"};
  Configurable<bool> isApplyNoCollInRofStandard{"isApplyNoCollInRofStandard", true, "Enable NoCollInRofStandard cut"};
  Configurable<bool> isApplyNoHighMultCollInPrevRof{"isApplyNoHighMultCollInPrevRof", true, "Enable NoHighMultCollInPrevRof cut"};
  Configurable<bool> isApplyFT0CbasedOccupancy{"isApplyFT0CbasedOccupancy", true, "Enable FT0CbasedOccupancy cut"};
  Configurable<bool> isApplyCentFT0C{"isApplyCentFT0C", false, "Centrality based on FT0C"};
  Configurable<bool> isApplyCentFT0CVariant1{"isApplyCentFT0CVariant1", false, "Centrality based on FT0C variant1"};
  Configurable<bool> isApplyCentFT0M{"isApplyCentFT0M", false, "Centrality based on FT0A + FT0C"};
  Configurable<bool> isApplyCentNGlobal{"isApplyCentNGlobal", false, "Centrality based on global tracks"};
  Configurable<bool> isApplyCentMFT{"isApplyCentMFT", false, "Centrality based on MFT tracks"};

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

    histos.add("EventHist", "EventHist", kTH1D, {axisEvent}, false);
    histos.add("VtxZHist", "VtxZHist", kTH1D, {axisVtxZ}, false);

    auto hstat = histos.get<TH1>(HIST("EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "kNoSameBunchPileup");  // reject collisions in case of pileup with another collision in the same foundBC
    x->SetBinLabel(4, "kIsGoodZvtxFT0vsPV");  // small difference between z-vertex from PV and from FT0
    x->SetBinLabel(5, "kIsVertexITSTPC");     // at least one ITS-TPC track (reject vertices built from ITS-only tracks)
    x->SetBinLabel(6, "kIsVertexTOFmatched"); // at least one of vertex contributors is matched to TOF
    x->SetBinLabel(7, "kIsVertexTRDmatched"); // at least one of vertex contributors is matched to TRD
    x->SetBinLabel(8, "Centrality");
    x->SetBinLabel(9, "ApplyExtraCorrCut");
    x->SetBinLabel(10, "ApplyNoCollInTimeRangeStandard");
    x->SetBinLabel(11, "ApplyNoCollInRofStandard");
    x->SetBinLabel(12, "ApplyNoHighMultCollInPrevRof");

    if (doprocessData) {
      histos.add("CentPercentileHist", "CentPercentileHist", kTH1D, {axisCent}, false);
      histos.add("hdatamult", "hdatamult", kTHnSparseD, {axisVtxZ, axisMult, centAxis}, false);
      histos.add("hdatadndeta", "hdatadndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi, axisTrackType}, false);
      histos.add("hdatazvtxcent", "hdatazvtxcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
      histos.add("PhiVsEtaHist", "PhiVsEtaHist", kTH2D, {axisPhi2, axisEta}, false);
    }

    if (doprocessMonteCarlo || doprocessMCpTefficiency || doprocessMCcheckFakeTracks || doprocessMCfillspecies) {
      histos.add("CentPercentileMCRecHist", "CentPercentileMCRecHist", kTH1D, {axisCent}, false);
      histos.add("hmczvtxcent", "hmczvtxcent", kTH3D, {axisVtxZ, centAxis, axisOccupancy}, false);
    }

    if (doprocessMonteCarlo) {
      histos.add("MCrecPhiVsEtaHist", "MCrecPhiVsEtaHist", kTH2D, {axisPhi2, axisEta}, false);
      histos.add("hmcrecdndeta", "hmcrecdndeta", kTHnSparseD, {axisVtxZ, centAxis, axisOccupancy, axisEta, axisPhi}, false);
      histos.add("hmcgendndeta", "hmcgendndeta", kTHnSparseD, {axisVtxZ, centAxis, axisEta, axisPhi, axisGenPtVary}, false);
    }

    if (doprocessMCfillspecies) {
      histos.add("FillMCrecSpecies", "FillMCrecSpecies", kTHnSparseD, {centAxis, axisOccupancy, axisEta, axisSpecies}, false);
      histos.add("FillMCgenSpecies", "FillMCgenSpecies", kTHnSparseD, {centAxis, axisEta, axisSpecies}, false);
    }

    if (doprocessMCpTefficiency) {
      histos.add("hmcrecdndpt", "hmcrecdndpt", kTHnSparseD, {centAxis, axisPt}, false);
      histos.add("hmcgendndpt", "hmcgendndpt", kTHnSparseD, {centAxis, axisPt, axisGenPtVary}, false);
    }

    if (doprocessMCcheckFakeTracks) {
      histos.add("hTracksCount", "hTracksCount", kTHnSparseD, {centAxis, axisTracks}, false);
      auto htrack = histos.get<THnSparse>(HIST("hTracksCount"));
      auto* x2 = htrack->GetAxis(1);
      x2->SetBinLabel(1, "All tracks");
      x2->SetBinLabel(2, "Non-fake tracks");
      for (int i = 0; i < 7; i++) {
        x2->SetBinLabel(i + 3, Form("layer %d", i));
      }
    }

    if (doprocessCorrelation) {
      histos.add("GlobalMult_vs_FT0A", "GlobalMult_vs_FT0A", kTH2F, {axisMult, axisFt0aMult}, true);
      histos.add("GlobalMult_vs_FT0C", "GlobalMult_vs_FT0C", kTH2F, {axisMult, axisFt0cMult}, true);
      histos.add("NPVtracks_vs_FT0C", "NPVtracks_vs_FT0C", kTH2F, {axisPV, axisFt0cMult}, true);
      histos.add("GlobalMult_vs_FV0A", "GlobalMult_vs_FV0A", kTH2F, {axisMult, axisFv0aMult}, true);
      histos.add("NPVtracks_vs_GlobalMult", "NPVtracks_vs_GlobalMult", kTH2F, {axisPV, axisMult}, true);
    }

    if (doprocessStrangeYield) {
      histos.add("hzvtxcent", "hzvtxcent", kTH2D, {axisVtxZ, centAxis}, false);
      histos.add("K0sCentEtaMass", "K0sCentEtaMass", kTH3D, {centAxis, axisEta, axisMassK0s}, false);
      histos.add("LambdaCentEtaMass", "LambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
      histos.add("AntiLambdaCentEtaMass", "AntiLambdaCentEtaMass", kTH3D, {centAxis, axisEta, axisMassLambda}, false);
    }
  }

  template <typename CheckTrack>
  bool isTrackSelected(CheckTrack const& track)
  {
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > 3.07666 && track.phi() < 3.12661) || track.phi() <= 0.03 || track.phi() >= 6.253)) {
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
    if (std::abs(pdgTrack->Charge()) < 3) {
      return false;
    }
    if (std::abs(track.eta()) >= etaRange) {
      return false;
    }
    if (isApplyExtraPhiCut && ((track.phi() > 3.07666 && track.phi() < 3.12661) || track.phi() <= 0.03 || track.phi() >= 6.253)) {
      return false;
    }
    return true;
  }
  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("EventHist"), 1);

    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("EventHist"), 2);

    if (isApplySameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 3);

    if (isApplyGoodZvtxFT0vsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 4);

    if (isApplyVertexITSTPC && !col.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 5);

    if (isApplyVertexTOFmatched && !col.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 6);

    if (isApplyVertexTRDmatched && !col.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 7);

    if (col.centFT0C() < 0. || col.centFT0C() > 100.) {
      return false;
    }
    histos.fill(HIST("EventHist"), 8);

    if (isApplyExtraCorrCut && col.multNTracksPV() > npvTracksCut && col.multFT0C() < (10 * col.multNTracksPV() - ft0cCut)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 9);

    if (isApplyNoCollInTimeRangeStandard && !col.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 10);

    if (isApplyNoCollInRofStandard && !col.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 11);

    if (isApplyNoHighMultCollInPrevRof && !col.selection_bit(o2::aod::evsel::kNoHighMultCollInPrevRof)) {
      return false;
    }
    histos.fill(HIST("EventHist"), 12);
    return true;
  }
  template <typename CheckColCent>
  float selectColCentrality(CheckColCent const& col)
  {
    auto cent = -1;
    if (isApplyCentFT0C) {
      cent = col.centFT0C();
    }
    if (isApplyCentFT0CVariant1) {
      cent = col.centFT0CVariant1();
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
  expressions::Filter trackSelectionProperMixed = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                                  ncheckbit(aod::track::trackCutFlag, TrackSelectionIts) &&
                                                  ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                         ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true) &&
                                                  ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                                         ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));

  void processData(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision)) {
      return;
    }
    histos.fill(HIST("VtxZHist"), collision.posZ());
    histos.fill(HIST("CentPercentileHist"), selectColCentrality(collision));
    auto occupancyValue = isApplyFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    histos.fill(HIST("hdatazvtxcent"), collision.posZ(), selectColCentrality(collision), occupancyValue);

    auto nchTracks = 0;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      histos.fill(HIST("PhiVsEtaHist"), track.phi(), track.eta());
      nchTracks++;
      histos.fill(HIST("hdatadndeta"), collision.posZ(), selectColCentrality(collision), occupancyValue, track.eta(), track.phi(), kGlobalplusITS);
      if (track.hasTPC()) {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), selectColCentrality(collision), occupancyValue, track.eta(), track.phi(), kGlobalonly);
      } else {
        histos.fill(HIST("hdatadndeta"), collision.posZ(), selectColCentrality(collision), occupancyValue, track.eta(), track.phi(), kITSonly);
      }
    }
    histos.fill(HIST("hdatamult"), collision.posZ(), nchTracks, selectColCentrality(collision));
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processData, "process data CentFT0C", false);

  void processCorrelation(CollisionDataTable::iterator const& collision, FilTrackDataTable const& tracks)
  {
    if (!isEventSelected(collision)) {
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
    histos.fill(HIST("NPVtracks_vs_FT0C"), collision.multNTracksPV(), collision.multFT0C());
    histos.fill(HIST("GlobalMult_vs_FV0A"), nchTracks, collision.multFV0A());
    histos.fill(HIST("NPVtracks_vs_GlobalMult"), collision.multNTracksPV(), nchTracks);
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processCorrelation, "do correlation study in data", false);

  void processMonteCarlo(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCollision : RecCollisions) {
      if (!isEventSelected(RecCollision)) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selectColCentrality(RecCollision));
      auto occupancyValue = isApplyFT0CbasedOccupancy ? RecCollision.ft0cOccupancyInTimeRange() : RecCollision.trackOccupancyInTimeRange();
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), selectColCentrality(RecCollision), occupancyValue);

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("MCrecPhiVsEtaHist"), Rectrack.phi(), Rectrack.eta());
        histos.fill(HIST("hmcrecdndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), occupancyValue, Rectrack.eta(), Rectrack.phi());
      } // track (mcrec) loop

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), particle.eta(), particle.phi(), kNoGenpTVar);
        if (particle.pt() < 0.1) {
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), particle.eta(), particle.phi(), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), particle.eta(), particle.phi(), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), particle.eta(), particle.phi(), kGenpTup);
          histos.fill(HIST("hmcgendndeta"), RecCollision.posZ(), selectColCentrality(RecCollision), particle.eta(), particle.phi(), kGenpTdown);
        }
      } // track (mcgen) loop
    } // collision loop
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processMonteCarlo, "process MC CentFT0C", false);

  void processMCpTefficiency(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCollision : RecCollisions) {
      if (!isEventSelected(RecCollision)) {
        continue;
      }
      if (std::abs(RecCollision.posZ()) >= vtxRange) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selectColCentrality(RecCollision));
      auto occupancyValue = isApplyFT0CbasedOccupancy ? RecCollision.ft0cOccupancyInTimeRange() : RecCollision.trackOccupancyInTimeRange();
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), selectColCentrality(RecCollision), occupancyValue);

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (Rectrack.has_mcParticle()) {
          auto mcpart = Rectrack.mcParticle();
          if (mcpart.isPhysicalPrimary()) {
            histos.fill(HIST("hmcrecdndpt"), selectColCentrality(RecCollision), mcpart.pt());
          }
        }
      }

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("hmcgendndpt"), selectColCentrality(RecCollision), particle.pt(), kNoGenpTVar);
        if (particle.pt() < 0.1) {
          histos.fill(HIST("hmcgendndpt"), selectColCentrality(RecCollision), particle.pt(), kGenpTup, -10.0 * particle.pt() + 2);
          histos.fill(HIST("hmcgendndpt"), selectColCentrality(RecCollision), particle.pt(), kGenpTdown, 5.0 * particle.pt() + 0.5);
        } else {
          histos.fill(HIST("hmcgendndpt"), selectColCentrality(RecCollision), particle.pt(), kGenpTup);
          histos.fill(HIST("hmcgendndpt"), selectColCentrality(RecCollision), particle.pt(), kGenpTdown);
        }
      }
    }
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processMCpTefficiency, "process MC pTefficiency", false);

  void processMCcheckFakeTracks(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCollision : RecCollisions) {
      if (!isEventSelected(RecCollision)) {
        continue;
      }
      if (std::abs(RecCollision.posZ()) >= vtxRange) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selectColCentrality(RecCollision));
      auto occupancyValue = isApplyFT0CbasedOccupancy ? RecCollision.ft0cOccupancyInTimeRange() : RecCollision.trackOccupancyInTimeRange();
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), selectColCentrality(RecCollision), occupancyValue);

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      for (const auto& Rectrack : recTracksPart) {
        if (std::abs(Rectrack.eta()) >= etaRange) {
          continue;
        }
        if (!Rectrack.hasTPC()) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), selectColCentrality(RecCollision), 1);
        bool isFakeItsTracks = false;
        for (int i = 0; i < 7; i++) {
          if (Rectrack.mcMask() & 1 << i) {
            isFakeItsTracks = true;
            histos.fill(HIST("hTracksCount"), selectColCentrality(RecCollision), i + 3);
            break;
          }
        }
        if (isFakeItsTracks) {
          continue;
        }
        histos.fill(HIST("hTracksCount"), selectColCentrality(RecCollision), 2);
      }
    }
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processMCcheckFakeTracks, "Check Fake tracks", false);

  void processMCfillspecies(CollisionMCTrueTable::iterator const&, CollisionMCRecTable const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    for (const auto& RecCollision : RecCollisions) {
      if (!isEventSelected(RecCollision)) {
        continue;
      }
      if (std::abs(RecCollision.posZ()) >= vtxRange) {
        continue;
      }
      histos.fill(HIST("VtxZHist"), RecCollision.posZ());
      histos.fill(HIST("CentPercentileMCRecHist"), selectColCentrality(RecCollision));
      auto occupancyValue = isApplyFT0CbasedOccupancy ? RecCollision.ft0cOccupancyInTimeRange() : RecCollision.trackOccupancyInTimeRange();
      histos.fill(HIST("hmczvtxcent"), RecCollision.posZ(), selectColCentrality(RecCollision), occupancyValue);

      auto recTracksPart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      std::vector<int> mclabels;
      for (const auto& Rectrack : recTracksPart) {
        if (!isTrackSelected(Rectrack)) {
          continue;
        }
        histos.fill(HIST("FillMCrecSpecies"), selectColCentrality(RecCollision), occupancyValue, Rectrack.eta(), static_cast<double>(kSpAll));
        if (Rectrack.has_mcParticle()) {
          int pid = kBkg;
          auto mcpart = Rectrack.template mcParticle_as<aod::McParticles>();
          if (mcpart.isPhysicalPrimary()) {
            switch (std::abs(mcpart.pdgCode())) {
              case 211:
                pid = kSpPion;
                break;
              case 321:
                pid = kSpKaon;
                break;
              case 2212:
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
            if (mcpartMother.pdgCode() == 310 || std::abs(mcpartMother.pdgCode()) == 3122) {
              pid = kSpStrangeDecay;
            }
          }
          if (find(mclabels.begin(), mclabels.end(), Rectrack.mcParticleId()) != mclabels.end()) {
            pid = kBkg;
          }
          mclabels.push_back(Rectrack.mcParticleId());
          histos.fill(HIST("FillMCrecSpecies"), selectColCentrality(RecCollision), occupancyValue, Rectrack.eta(), static_cast<double>(pid));
        } else {
          histos.fill(HIST("FillMCrecSpecies"), selectColCentrality(RecCollision), occupancyValue, Rectrack.eta(), static_cast<double>(kBkg));
        }
      } // rec track loop

      for (const auto& particle : GenParticles) {
        if (!isGenTrackSelected(particle)) {
          continue;
        }
        histos.fill(HIST("FillMCgenSpecies"), selectColCentrality(RecCollision), particle.eta(), static_cast<double>(kSpAll));
        int pid = 0;
        switch (std::abs(particle.pdgCode())) {
          case 211:
            pid = kSpPion;
            break;
          case 321:
            pid = kSpKaon;
            break;
          case 2212:
            pid = kSpProton;
            break;
          default:
            pid = kSpOther;
            break;
        }
        histos.fill(HIST("FillMCgenSpecies"), selectColCentrality(RecCollision), particle.eta(), static_cast<double>(pid));
      } // gen track loop
    } // collision loop
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processMCfillspecies, "Fill particle species in MC", false);

  void processStrangeYield(CollisionDataTable::iterator const& collision, V0TrackCandidates const&, aod::V0Datas const& v0data)
  {
    if (!isEventSelected(collision)) {
      return;
    }
    if (std::abs(collision.posZ()) >= vtxRange) {
      return;
    }
    histos.fill(HIST("hzvtxcent"), collision.posZ(), selectColCentrality(collision));
    for (const auto& v0track : v0data) {
      auto v0pTrack = v0track.template posTrack_as<V0TrackCandidates>();
      auto v0nTrack = v0track.template negTrack_as<V0TrackCandidates>();
      if (std::abs(v0pTrack.eta()) > 0.9 || std::abs(v0nTrack.eta()) > 0.9) {
        continue;
      }
      if (v0pTrack.tpcNClsFound() < minTPCnClsCut) {
        continue;
      }
      if (v0nTrack.tpcNClsFound() < minTPCnClsCut) {
        continue;
      }
      if (std::abs(v0pTrack.tpcNSigmaPi()) > nSigmaTpcCut) {
        continue;
      }
      if (std::abs(v0nTrack.tpcNSigmaPi()) > nSigmaTpcCut) {
        continue;
      }
      if (std::abs(v0pTrack.tpcNSigmaPr()) > nSigmaTpcCut) {
        continue;
      }
      if (std::abs(v0nTrack.tpcNSigmaPr()) > nSigmaTpcCut) {
        continue;
      }
      if (std::abs(v0track.dcapostopv()) < dcapostopvCut || std::abs(v0track.dcanegtopv()) < dcanegtopvCut || v0track.v0radius() < v0radiusCut || v0track.v0cosPA() < v0cospaCut || std::abs(v0track.dcaV0daughters()) > dcav0daughtercut) {
        continue;
      }
      histos.fill(HIST("K0sCentEtaMass"), selectColCentrality(collision), v0track.eta(), v0track.mK0Short());
      histos.fill(HIST("LambdaCentEtaMass"), selectColCentrality(collision), v0track.eta(), v0track.mLambda());
      histos.fill(HIST("AntiLambdaCentEtaMass"), selectColCentrality(collision), v0track.eta(), v0track.mAntiLambda());
    }
  }
  PROCESS_SWITCH(HeavyionMultiplicity, processStrangeYield, "Strange particle yield", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyionMultiplicity>(cfgc)};
}
