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
/// \file longrangeCorrelation.cxx
///
/// \brief task for long range correlation analysis
/// \author Abhi Modak (abhi.modak@cern.ch) and Debojit sarkar (debojit.sarkar@cern.ch)
/// \since April 22, 2025

#include <TH1F.h>
#include <chrono>
#include <string>
#include <vector>
#include <TComplex.h>
#include <TMath.h>
#include <cstdio>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/Core/RecoDecay.h"
#include "CommonConstants/MathConstants.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;
using namespace o2::constants::math;

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

AxisSpec axisEvent{10, 0.5, 9.5, "#Event", "EventAxis"};
AxisSpec amplitudeFT0{5000, 0, 10000, "FT0 amplitude"};
AxisSpec channelFT0Axis{96, 0.0, 96.0, "FT0 channel"};

struct LongrangeCorrelation {

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  SliceCache cache;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  std::vector<o2::detectors::AlignParam>* offsetFT0;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 1.0f, "Eta range to consider"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<float> cfgPtCutMin{"cfgPtCutMin", 0.2f, "minimum accepted track pT"};
  Configurable<float> cfgPtCutMax{"cfgPtCutMax", 10.0f, "maximum accepted track pT"};
  Configurable<int> mixingParameter{"mixingParameter", 5, "how many events are mixed"};
  Configurable<int> cfgMinMult{"cfgMinMult", 0, "Minimum multiplicity for collision"};
  Configurable<int> cfgMaxMult{"cfgMaxMult", 10, "Maximum multiplicity for collision"};
  Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for mixed event"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "Mixing bins - z-vertex"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
  ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
  ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
  ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 80, 100}, "multiplicity / centrality axis for histograms"};

  using CollTable = soa::Join<aod::Collisions, aod::EvSels>;
  using TrksTable = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  Preslice<TrksTable> perCollision = aod::track::collisionId;

  OutputObj<CorrelationContainer> same{Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};
  OutputObj<CorrelationContainer> mixed{Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult))};

  void init(InitContext const&)
  {
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    LOGF(info, "Getting alignment offsets from the CCDB...");
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", cfgCcdbParam.noLaterThan.value);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());

    // QA histos
    histos.add("QA/EventHist", "events", kTH1F, {axisEvent}, false);
    histos.add("QA/VtxZHist", "v_{z} (cm)", kTH1F, {axisVtxZ}, false);
    histos.add("QA/hMEpvz1", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("QA/hMEpvz2", ";pvz;Entries", kTH1F, {{30, -15, 15}});
    histos.add("QA/hMixingQA", "events", kTH1F, {axisEvent}, false);

    auto hstat = histos.get<TH1>(HIST("QA/EventHist"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "All events");
    x->SetBinLabel(2, "sel8");
    x->SetBinLabel(3, "|vz|<10");

    histos.add("SE/hMult", "event multiplicity", kTH1D, {axisMultiplicity});
    histos.add("SE/Trig_etavsphi", ";#eta;#phi", kTH2D, {axisPhi, axisEtaTrig});
    histos.add("SE/Trig_eta", "#eta", kTH1D, {axisEtaTrig});
    histos.add("SE/Trig_phi", "#eta", kTH1D, {axisPhi});
    histos.add("SE/Trig_pt", "p_{T}", kTH1D, {axisPtTrigger});
    histos.add("SE/hMult_used", "event multiplicity", kTH1F, {axisMultiplicity});
    histos.add("SE/Trig_hist", "trig hist", kTHnSparseF, {axisSample, axisVtxZ, axisPtTrigger});
    histos.add("SE/FT0Amp", "ftoamult", kTH2D, {channelFT0Axis, amplitudeFT0});
    histos.add("SE/FT0Aeta", "ft0a;#eta", kTH1D, {axisEtaAssoc});
    histos.add("SE/FT0Aphi", "ft0a;#phi", kTH1D, {axisPhi});
    histos.add("SE/FT0Aetavsphi", ";ft0a;#eta;#phi", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add("SE/deltaEta_deltaPhi", ";#delta#eta;#delta#phi", kTH2D, {axisDeltaPhi, axisDeltaEta});

    histos.add("ME/hMult", "event multiplicity", kTH1D, {axisMultiplicity});
    histos.add("ME/Trig_etavsphi", ";#eta;#phi", kTH2D, {axisPhi, axisEtaTrig});
    histos.add("ME/Trig_eta", "#eta", kTH1D, {axisEtaTrig});
    histos.add("ME/Trig_phi", "#eta", kTH1D, {axisPhi});
    histos.add("ME/Trig_pt", "p_{T}", kTH1D, {axisPtTrigger});
    histos.add("ME/FT0Amp", "ftoamult", kTH2D, {channelFT0Axis, amplitudeFT0});
    histos.add("ME/FT0Aeta", "ft0a;#eta", kTH1D, {axisEtaAssoc});
    histos.add("ME/FT0Aphi", "ft0a;#phi", kTH1D, {axisPhi});
    histos.add("ME/FT0Aetavsphi", ";ft0a;#eta;#phi", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add("ME/deltaEta_deltaPhi", ";#delta#eta;#delta#phi", kTH2D, {axisDeltaPhi, axisDeltaEta});

    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVtxZ, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{axisVertexEfficiency, "z-vtx (cm)"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisEtaEfficiency, "#eta"}};

    std::vector<AxisSpec> userAxis;

    same.setObject(new CorrelationContainer(Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("sameEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer(Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), Form("mixedEvent_%i_%i", static_cast<int>(cfgMinMult), static_cast<int>(cfgMaxMult)), corrAxis, effAxis, userAxis));
  }

  double getPhiFT0(int chno, double offsetX, double offsetY)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    return RecoDecay::phi(chPos.X() + offsetX, chPos.Y() + offsetY);
  }

  double getEtaFT0(int chno, double offsetX, double offsetY, double offsetZ)
  {
    o2::ft0::Geometry ft0Det;
    ft0Det.calculateChannelCenter();
    auto chPos = ft0Det.getChannelCenter(chno);
    auto x = chPos.X() + offsetX;
    auto y = chPos.Y() + offsetY;
    auto z = chPos.Z() + offsetZ;
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                              ncheckbit(aod::track::trackCutFlag, TrackSelectionIts);
  Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionTpc), true);
  Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, TrackSelectionDcaxyOnly),
                                     ncheckbit(aod::track::trackCutFlag, TrackSelectionDca));
  Filter fTracksEta = nabs(aod::track::eta) < cfgEtaCut;
  Filter fTracksPt = (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax);

  template <typename CheckCol>
  bool isEventSelected(CheckCol const& col)
  {
    histos.fill(HIST("QA/EventHist"), 1);
    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("QA/EventHist"), 2);
    if (std::abs(col.posZ()) >= cfgVtxCut) {
      return false;
    }
    histos.fill(HIST("QA/EventHist"), 3);
    histos.fill(HIST("QA/VtxZHist"), col.posZ());
    return true;
  }

  template <typename TTracks>
  void fillYield(TTracks tracks, bool mixing)
  {
    if (mixing) {
      histos.fill(HIST("ME/hMult"), tracks.size());
      for (auto const& triggerTrack : tracks) {
        histos.fill(HIST("ME/Trig_etavsphi"), triggerTrack.phi(), triggerTrack.eta());
        histos.fill(HIST("ME/Trig_eta"), triggerTrack.eta());
        histos.fill(HIST("ME/Trig_phi"), triggerTrack.phi());
        histos.fill(HIST("ME/Trig_pt"), triggerTrack.pt());
      }
    } else {
      histos.fill(HIST("SE/hMult"), tracks.size());
      for (auto const& triggerTrack : tracks) {
        histos.fill(HIST("SE/Trig_etavsphi"), triggerTrack.phi(), triggerTrack.eta());
        histos.fill(HIST("SE/Trig_eta"), triggerTrack.eta());
        histos.fill(HIST("SE/Trig_phi"), triggerTrack.phi());
        histos.fill(HIST("SE/Trig_pt"), triggerTrack.pt());
      }
    }
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TFT0s>
  void fillCorrelation(TTarget target, TTriggers const& triggers, TFT0s const& ft0, bool mixing, float vz)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    if (!mixing)
      histos.fill(HIST("SE/hMult_used"), triggers.size());
    for (auto const& triggerTrack : triggers) {
      if (!mixing)
        histos.fill(HIST("SE/Trig_hist"), fSampleIndex, vz, triggerTrack.pt());

      auto offsetFT0Ax = (*offsetFT0)[0].getX();
      auto offsetFT0Ay = (*offsetFT0)[0].getY();
      auto offsetFT0Az = (*offsetFT0)[0].getZ();
      for (std::size_t iChA = 0; iChA < ft0.channelA().size(); iChA++) {
        auto chanelid = ft0.channelA()[iChA];
        float ampl = ft0.amplitudeA()[iChA];
        if (ampl <= 0)
          continue;
        if (mixing)
          histos.fill(HIST("ME/FT0Amp"), chanelid, ampl);
        else
          histos.fill(HIST("SE/FT0Amp"), chanelid, ampl);

        auto phiA = getPhiFT0(chanelid, offsetFT0Ax, offsetFT0Ay);
        auto etaA = getEtaFT0(chanelid, offsetFT0Ax, offsetFT0Ay, offsetFT0Az);

        if (mixing) {
          histos.fill(HIST("ME/FT0Aeta"), etaA);
          histos.fill(HIST("ME/FT0Aphi"), phiA);
          histos.fill(HIST("ME/FT0Aetavsphi"), phiA, etaA);
        } else {
          histos.fill(HIST("SE/FT0Aeta"), etaA);
          histos.fill(HIST("SE/FT0Aphi"), phiA);
          histos.fill(HIST("SE/FT0Aetavsphi"), phiA, etaA);
        }
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - phiA, -PIHalf);
        float deltaEta = triggerTrack.eta() - etaA;
        if (mixing)
          histos.fill(HIST("ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        else
          histos.fill(HIST("SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta);
      } // associated ft0 tracks
    } // trigger tracks
  } // fillCorrelation

  void processSE(CollTable::iterator const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    if (!isEventSelected(col)) {
      return;
    }
    if (col.has_foundFT0()) {
      fillYield(tracks, false);
      const auto& ft0 = col.foundFT0();
      if (tracks.size() < cfgMinMult || tracks.size() >= cfgMaxMult) {
        return;
      }
      fillCorrelation<CorrelationContainer::kCFStepReconstructed>(same, tracks, ft0, false, col.posZ());
    }
  } // same event

  void processME(CollTable const& col, aod::FT0s const&, TrksTable const& tracks)
  {
    auto getTracksSize = [&tracks, this](CollTable::iterator const& collision) {
      auto associatedTracks = tracks.sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), this->cache);
      return associatedTracks.size();
    };
    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getTracksSize)>, aod::collision::PosZ, decltype(getTracksSize)>;
    MixedBinning binningOnVtxAndMult{{getTracksSize}, {axisVtxZME, axisMultME}, true};
    for (auto const& [col1, col2] : soa::selfCombinations(binningOnVtxAndMult, mixingParameter, -1, col, col)) {
      if (!isEventSelected(col1) || !isEventSelected(col2)) {
        continue;
      }
      if (col1.globalIndex() == col2.globalIndex()) {
        histos.fill(HIST("QA/hMixingQA"), 1.0); // same-collision pair counting
        continue;
      }
      if (col1.has_foundFT0() && col2.has_foundFT0()) {
        histos.fill(HIST("QA/hMEpvz1"), col1.posZ());
        histos.fill(HIST("QA/hMEpvz2"), col2.posZ());
        auto slicedTriggerTracks = tracks.sliceBy(perCollision, col1.globalIndex());
        fillYield(slicedTriggerTracks, true);
        const auto& ft0 = col2.foundFT0();
        if (slicedTriggerTracks.size() < cfgMinMult || slicedTriggerTracks.size() >= cfgMaxMult) {
          continue;
        }
        fillCorrelation<CorrelationContainer::kCFStepReconstructed>(mixed, slicedTriggerTracks, ft0, true, col1.posZ());
      }
    }
  } // mixed event

  PROCESS_SWITCH(LongrangeCorrelation, processSE, "process same event", false);
  PROCESS_SWITCH(LongrangeCorrelation, processME, "process mixed event", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LongrangeCorrelation>(cfgc)};
}
