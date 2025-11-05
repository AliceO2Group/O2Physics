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
/// \file longrangecorrDerived.cxx
///
/// \brief task for long range correlation analysis based on derived table
/// \author Abhi Modak (abhi.modak@cern.ch)
/// \since November 05, 2025

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/TwoParticleCorrelations/DataModel/LongRangeDerived.h"
#include "PWGMM/Mult/DataModel/bestCollisionTable.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "FV0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"

#include <TComplex.h>
#include <TH1F.h>
#include <TMath.h>
#include <TPDGCode.h>

#include <chrono>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::fwdtrack;
using namespace o2::aod::evsel;
using namespace o2::constants::math;

struct LongrangecorrDerived {

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for mixed event"};
  Configurable<int> cfgNmixedevent{"cfgNmixedevent", 5, "how many events are mixed"};
  Configurable<int> cfgPidMask{"cfgPidMask", 0, "Selection bitmask for the TPC particle"};
  Configurable<int> cfgV0Mask{"cfgV0Mask", 0, "Selection bitmask for the V0 particle"};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 15, 25, 50, 60, 1000}, "multiplicity axis"};
  ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
  ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt assoc axis for histograms"};
  ConfigurableAxis axisSample{"axisSample", {cfgSampleSize, 0, cfgSampleSize}, "sample axis for histograms"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
  ConfigurableAxis channelFt0aAxis{"channelFt0aAxis", {96, 0.0, 96.0}, "FT0A channel"};
  ConfigurableAxis amplitudeFt0a{"amplitudeFt0a", {5000, 0, 10000}, "FT0A amplitude"};
  ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
  ConfigurableAxis axisInvMass{"axisInvMass", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0}, "invariant mass axis"};
  ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 1000}, "Mixing bins - multiplicity"};
  ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "Mixing bins - z-vertex"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  void init(InitContext const&)
  {
    std::vector<AxisSpec> corrAxis = {{axisSample, "Sample"},
                                      {axisVtxZ, "z-vtx (cm)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{axisVertexEfficiency, "z-vtx (cm)"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisEtaEfficiency, "#eta"}};
    std::vector<AxisSpec> userAxis = {{axisMultiplicity, "multiplicity"},
                                      {axisInvMass, "m (GeV/c^2)"}};

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    histos.add("hMultiplicity", "hMultiplicity", kTH1D, {axisMultiplicity});
    histos.add("hCentrality", "hCentrality", kTH1D, {axisMultiplicity});
    histos.add("hVertexZ", "hVertexZ", kTH1D, {axisVtxZ});

    histos.add("Trig_eta", "Trig_eta", kTH1D, {axisEtaTrig});
    histos.add("Trig_phi", "Trig_phi", kTH1D, {axisPhi});
    histos.add("Trig_etavsphi", "Trig_etavsphi", kTH2D, {axisPhi, axisEtaTrig});
    histos.add("Trig_pt", "Trig_pt", kTH1D, {axisPtTrigger});
    histos.add("Trig_hist", "Trig_hist", kTHnSparseF, {axisSample, axisVtxZ, axisPtTrigger, axisMultiplicity, axisInvMass});

    histos.add("Assoc_eta", "Assoc_eta", kTH1D, {axisEtaAssoc});
    histos.add("Assoc_phi", "Assoc_phi", kTH1D, {axisPhi});
    histos.add("Assoc_etavsphi", "Assoc_etavsphi", kTH2D, {axisPhi, axisEtaAssoc});
    histos.add("Assoc_pt", "Assoc_pt", kTH1D, {axisPtAssoc});

    histos.add("deltaEta_deltaPhi_same", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
    histos.add("deltaEta_deltaPhi_mixed", "", kTH2D, {axisDeltaPhi, axisDeltaEta});
  }

  template <typename TCollision>
  void fillCollQA(TCollision const& col)
  {
    histos.fill(HIST("hMultiplicity"), col.multiplicity());
    histos.fill(HIST("hCentrality"), col.centrality());
    histos.fill(HIST("hVertexZ"), col.zvtx());
  }

  template <typename TTrack>
  void fillTrigTrackQA(TTrack const& track)
  {
    histos.fill(HIST("Trig_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Trig_eta"), track.eta());
    histos.fill(HIST("Trig_phi"), track.phi());
    histos.fill(HIST("Trig_pt"), track.pt());
  }

  template <typename TTrack>
  void fillAssocTrackQA(TTrack const& track)
  {
    histos.fill(HIST("Assoc_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Assoc_eta"), track.eta());
    histos.fill(HIST("Assoc_phi"), track.phi());
  }

  template <class T>
  using HasTpcTrack = decltype(std::declval<T&>().trackType());
  template <class T>
  using HasV0Track = decltype(std::declval<T&>().v0Type());
  template <class T>
  using HasInvMass = decltype(std::declval<T&>().invMass());

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TAssocs>
  void fillCorrHist(TTarget target, TTriggers const& triggers, TAssocs const& assocs, bool mixing, float vz, float multiplicity, float eventWeight)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      if constexpr (std::experimental::is_detected<HasTpcTrack, typename TTriggers::iterator>::value) {
        if (cfgPidMask != 0 && (cfgPidMask & (1u << static_cast<uint32_t>(triggerTrack.trackType()))) == 0u)
          continue;
      } else if constexpr (std::experimental::is_detected<HasV0Track, typename TTriggers::iterator>::value) {
        if (cfgV0Mask != 0 && (cfgV0Mask & (1u << static_cast<uint32_t>(triggerTrack.v0Type()))) == 0u)
          continue;
      }
      if (!mixing) {
        fillTrigTrackQA(triggerTrack);
        if constexpr (std::experimental::is_detected<HasInvMass, typename TTriggers::iterator>::value) {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity, triggerTrack.invMass(), eventWeight);
        } else {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, triggerTrack.pt(), multiplicity, 1.0, eventWeight);
        }
      }
      for (auto const& assoTrack : assocs) {
        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - assoTrack.phi(), -PIHalf);
        float deltaEta = triggerTrack.eta() - assoTrack.eta();
        if (!mixing) {
          fillAssocTrackQA(assoTrack);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta);
        } else {
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta);
        }
        if constexpr (std::experimental::is_detected<HasInvMass, typename TTriggers::iterator>::value) {
          target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta, multiplicity, triggerTrack.invMass(), eventWeight);
        } else {
          target->getPairHist()->Fill(step, fSampleIndex, vz, triggerTrack.pt(), triggerTrack.pt(), deltaPhi, deltaEta, multiplicity, 0., eventWeight);
        }
      } // associated tracks
    } // trigger tracks
  } // fill correlation

  template <typename TCollision, typename TTriggers, typename TAssocs>
  void processSame(TCollision const& col, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(col.zvtx()) >= cfgVtxCut) {
      return;
    }
    fillCollQA(col);
    fillCorrHist<CorrelationContainer::kCFStepReconstructed>(same, triggers, assocs, false, col.zvtx(), col.multiplicity(), 1.0);
  } // process same

  template <typename TCollision, typename... TrackTypes>
  void processMixed(TCollision const& col, TrackTypes&&... tracks)
  {
    auto getMultiplicity = [this](auto& collision) {
      (void)this;
      return collision.multiplicity();
    };
    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::lrcorrcolltable::Zvtx, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {axisVtxZME, axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<TCollision, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgNmixedevent, -1, col, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      float eventweight = 1.0f / it.currentWindowNeighbours();
      fillCorrHist<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, true, col1.zvtx(), col1.multiplicity(), eventweight);
    } // pair loop
  } // process mixed

  void processTpcft0aSE(aod::CollLRTables::iterator const& col, aod::TrkLRTables const& tracks, aod::Ft0aLRTables const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processTpcft0cSE(aod::CollLRTables::iterator const& col, aod::TrkLRTables const& tracks, aod::Ft0cLRTables const& ft0cs)
  {
    processSame(col, tracks, ft0cs);
  }

  void processTpcmftSE(aod::CollLRTables::iterator const& col, aod::TrkLRTables const& tracks, aod::MftTrkLRTables const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processMftft0aSE(aod::CollLRTables::iterator const& col, aod::MftTrkLRTables const& mfts, aod::Ft0aLRTables const& ft0as)
  {
    processSame(col, mfts, ft0as);
  }

  void processV0ft0aSE(aod::CollLRTables::iterator const& col, aod::V0TrkLRTables const& tracks, aod::Ft0aLRTables const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processV0mftSE(aod::CollLRTables::iterator const& col, aod::V0TrkLRTables const& tracks, aod::MftTrkLRTables const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processTpcmftbestSE(aod::CollLRTables::iterator const& col, aod::TrkLRTables const& tracks, aod::MftBestTrkLRTables const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processMftbestft0aSE(aod::CollLRTables::iterator const& col, aod::MftBestTrkLRTables const& mfts, aod::Ft0aLRTables const& ft0as)
  {
    processSame(col, mfts, ft0as);
  }

  void processV0mftbestSE(aod::CollLRTables::iterator const& col, aod::V0TrkLRTables const& tracks, aod::MftBestTrkLRTables const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processTpcft0aME(aod::CollLRTables const& col, aod::TrkLRTables const& tracks, aod::Ft0aLRTables const& ft0as)
  {
    processMixed(col, tracks, ft0as);
  }

  void processTpcft0cME(aod::CollLRTables const& col, aod::TrkLRTables const& tracks, aod::Ft0cLRTables const& ft0cs)
  {
    processMixed(col, tracks, ft0cs);
  }

  void processTpcmftME(aod::CollLRTables const& col, aod::TrkLRTables const& tracks, aod::MftTrkLRTables const& mfts)
  {
    processMixed(col, tracks, mfts);
  }

  void processMftft0aME(aod::CollLRTables const& col, aod::MftTrkLRTables const& mfts, aod::Ft0aLRTables const& ft0as)
  {
    processMixed(col, mfts, ft0as);
  }

  void processV0ft0aME(aod::CollLRTables const& col, aod::V0TrkLRTables const& tracks, aod::Ft0aLRTables const& ft0as)
  {
    processMixed(col, tracks, ft0as);
  }

  void processV0mftME(aod::CollLRTables const& col, aod::V0TrkLRTables const& tracks, aod::MftTrkLRTables const& mfts)
  {
    processMixed(col, tracks, mfts);
  }

  void processTpcmftbestME(aod::CollLRTables const& col, aod::TrkLRTables const& tracks, aod::MftBestTrkLRTables const& mfts)
  {
    processMixed(col, tracks, mfts);
  }

  void processMftbestft0aME(aod::CollLRTables const& col, aod::MftBestTrkLRTables const& mfts, aod::Ft0aLRTables const& ft0as)
  {
    processMixed(col, mfts, ft0as);
  }

  void processV0mftbestME(aod::CollLRTables const& col, aod::V0TrkLRTables const& tracks, aod::MftBestTrkLRTables const& mfts)
  {
    processMixed(col, tracks, mfts);
  }

  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0aSE, "same event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0aME, "mixed event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0cSE, "same event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcft0cME, "mixed event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcmftSE, "same event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcmftME, "mixed event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMftft0aSE, "same event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMftft0aME, "mixed event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0aSE, "same event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0ft0aME, "mixed event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0mftSE, "same event V0 vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0mftME, "mixed event V0 vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcmftbestSE, "same event TPC vs best MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processTpcmftbestME, "mixed event TPC vs best MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMftbestft0aSE, "same event best MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMftbestft0aME, "mixed event best MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0mftbestSE, "same event V0 vs best MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processV0mftbestME, "mixed event V0 vs best MFT", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongrangecorrDerived>(cfgc)};
}
