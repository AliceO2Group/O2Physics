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
#include "PWGCF/TwoParticleCorrelations/DataModel/LongRangeDerived.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TRandom.h>

#include <cstdint>
#include <cstdio>
#include <string>
#include <tuple>
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
  SGSelector sgSelector;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    Configurable<int> cfgNmixedevent{"cfgNmixedevent", 5, "how many events are mixed"};
    Configurable<double> cfgSampleSize{"cfgSampleSize", 10, "Sample size for bootstrapping"};
    Configurable<int> cfgPidMask{"cfgPidMask", 0, "Selection bitmask for the TPC particle"};
    Configurable<int> cfgV0Mask{"cfgV0Mask", 0, "Selection bitmask for the V0 particle"};
    Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0f, "Vertex Z range to consider"};
    Configurable<bool> isUseCentEst{"isUseCentEst", false, "Centrality based classification"};
    Configurable<int> isUseDataLikeMult{"isUseDataLikeMult", 0, "Data like mult/cent classification"};

    Configurable<float> cfgTpcMinNclsFound{"cfgTpcMinNclsFound", 50.0f, ""};
    Configurable<float> cfgTpcMinNCrossedRows{"cfgTpcMinNCrossedRows", 70.0f, ""};
    Configurable<float> cfgTpcMaxChi2PerCluster{"cfgTpcMaxChi2PerCluster", 4.0f, ""};
    Configurable<float> cfgTpcMaxDcaZ{"cfgTpcMaxDcaZ", 1.0f, ""};

    Configurable<int> cfgMftCluster{"cfgMftCluster", 5, "cut on MFT Cluster"};
    Configurable<float> cfgMftDcaxy{"cfgMftDcaxy", 2.0f, "cut on DCA xy for MFT tracks"};
    Configurable<float> cfgMftDcaz{"cfgMftDcaz", 2.0f, "cut on DCA z for MFT tracks"};
    Configurable<bool> cfgRejectAmbTrk{"cfgRejectAmbTrk", false, "Condition to reject Ambiguous tracks"};
    Configurable<bool> cfgRejectNonAmbTrk{"cfgRejectNonAmbTrk", false, "Condition to reject Non-Ambiguous tracks"};
  } cfgSel;

  struct : ConfigurableGroup {
    ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 10, 15, 25, 50, 60, 1000}, "multiplicity axis"};
    ConfigurableAxis axisPhi{"axisPhi", {96, 0, TwoPI}, "#phi axis"};
    ConfigurableAxis axisEtaTrig{"axisEtaTrig", {40, -1., 1.}, "#eta trig axis"};
    ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
    ConfigurableAxis axisVtxZ{"axisVtxZ", {40, -20, 20}, "vertex axis"};
    ConfigurableAxis axisEtaAssoc{"axisEtaAssoc", {96, 3.5, 4.9}, "#eta assoc axis"};
    ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
    ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -6, -2}, "delta eta axis for histograms"};
    ConfigurableAxis axisInvMass{"axisInvMass", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0}, "invariant mass axis"};
    ConfigurableAxis axisInvMassQA{"axisInvMassQA", {20, 0.45, 0.55}, "invariant mass axis for QA"};
    ConfigurableAxis axisAmplitude{"axisAmplitude", {5000, 0, 10000}, "FT0 amplitude"};
    ConfigurableAxis axisChannel{"axisChannel", {208, 0, 208}, "FT0 channel"};
    ConfigurableAxis axisMultME{"axisMultME", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 1000}, "Mixing bins - multiplicity"};
    ConfigurableAxis axisVtxZME{"axisVtxZME", {VARIABLE_WIDTH, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}, "Mixing bins - z-vertex"};
    ConfigurableAxis axisSample{"axisSample", {10, 0, 10}, "sample axis for histograms"};

    ConfigurableAxis axisTPCNClsFound{"axisTPCNClsFound", {200, -0.5, 199.5}, "TPC Cluster axis"};
    ConfigurableAxis axisTPCNClsCrossedRows{"axisTPCNClsCrossedRows", {200, -0.5, 199.5}, "TPC NCrossedRow axis"};
    ConfigurableAxis axisTPCChi2NCl{"axisTPCChi2NCl", {20, 0.0, 20.0}, "TPC Chi2/NCl axis"};
    ConfigurableAxis axisTPCdcaZ{"axisTPCdcaZ", {200, -10.0, 10.0}, "TPC dcaZ axis"};

    ConfigurableAxis axisMFTAmbDegree{"axisMFTAmbDegree", {50, -0.5, 49.5}, "Track Ambiguity axis"};
    ConfigurableAxis axisMFTNClusters{"axisMFTNClusters", {200, -0.5, 199.5}, "MFT Cluster axis"};
    ConfigurableAxis axisMFTbestDCAXY{"axisMFTbestDCAXY", {200, -10.0, 10.0}, "MFT dcaXY axis"};
    ConfigurableAxis axisMFTbestDCAZ{"axisMFTbestDCAZ", {200, -10.0, 10.0}, "MFT dcaZ axis"};

    ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
    ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
    ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {1, 0.5, 8.0}, "pt axis for efficiency histograms"};
  } cfgAxis;

  Configurable<float> cfgFv0Cut{"cfgFv0Cut", 50.0f, "FV0A threshold"};
  Configurable<float> cfgFt0aCut{"cfgFt0aCut", 100.0f, "FT0A threshold"};
  Configurable<float> cfgFt0cCut{"cfgFt0cCut", 50.0f, "FT0C threshold"};
  Configurable<float> cfgZdcCut{"cfgZdcCut", 0.1f, "ZDC threshold"};
  Configurable<int> cfgGapSideCut{"cfgGapSideCut", 0, "Gap-side A=0, C=1, AC = 2, No Gap = -1, All events = 3"};

  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  using CollsTable = aod::LRCollisions;
  using TrksTable = aod::LRMidTracks;
  using MftTrksTable = aod::LRMftTracks;
  using Ft0aTrksTable = aod::LRFt0aTracks;
  using Ft0cTrksTable = aod::LRFt0cTracks;
  using V0TrksTable = aod::LRV0Tracks;

  using McCollsTable = aod::LRMcCollisions;
  using McTrksTable = aod::LRMidMcTracks;
  using McMftTrksTable = aod::LRMftMcTracks;
  using McFt0aTrksTable = aod::LRFt0aMcTracks;
  using McFt0cTrksTable = aod::LRFt0cMcTracks;

  using UpcCollsTable = soa::Join<aod::UpcLRCollisions, aod::UpcSgLRCollisions, aod::LRZdcs>;
  using TrksUpcTable = aod::UpcLRMidTracks;
  using MftTrksUpcTable = aod::UpcLRMftTracks;
  using Ft0aTrksUpcTable = aod::UpcLRFt0aTracks;
  using Ft0cTrksUpcTable = aod::UpcLRFt0cTracks;
  using V0TrksUpcTable = aod::UpcLRV0Tracks;

  Preslice<TrksTable> perColTpc = aod::lrcorrtrktable::lrCollisionId;
  Preslice<MftTrksTable> perColMft = aod::lrcorrtrktable::lrCollisionId;
  Preslice<Ft0aTrksTable> perColFt0a = aod::lrcorrtrktable::lrCollisionId;
  Preslice<Ft0cTrksTable> perColFt0c = aod::lrcorrtrktable::lrCollisionId;
  Preslice<V0TrksTable> perColV0 = aod::lrcorrtrktable::lrCollisionId;

  Preslice<TrksUpcTable> perUpcColTpc = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<MftTrksUpcTable> perUpcColMft = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<Ft0aTrksUpcTable> perUpcColFt0a = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<Ft0cTrksUpcTable> perUpcColFt0c = aod::lrcorrtrktable::upcLRCollisionId;
  Preslice<V0TrksUpcTable> perUpcColV0 = aod::lrcorrtrktable::upcLRCollisionId;

  Preslice<McTrksTable> perMcColTpc = aod::lrcorrmctrktable::lrMcCollisionId;
  Preslice<McMftTrksTable> perMcColMft = aod::lrcorrmctrktable::lrMcCollisionId;
  Preslice<McFt0aTrksTable> perMcColFt0a = aod::lrcorrmctrktable::lrMcCollisionId;
  Preslice<McFt0cTrksTable> perMcColFt0c = aod::lrcorrmctrktable::lrMcCollisionId;

  void init(InitContext const&)
  {
    std::vector<AxisSpec> corrAxis = {{cfgAxis.axisSample, "Sample"},
                                      {cfgAxis.axisVtxZ, "z-vtx (cm)"},
                                      {cfgAxis.axisMultiplicity, "multiplicity"},
                                      {cfgAxis.axisPtTrigger, "p_{T} (GeV/c)"},
                                      {cfgAxis.axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {cfgAxis.axisDeltaEta, "#Delta#eta"}};
    std::vector<AxisSpec> effAxis = {{cfgAxis.axisVertexEfficiency, "z-vtx (cm)"},
                                     {cfgAxis.axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {cfgAxis.axisEtaEfficiency, "#eta"}};
    std::vector<AxisSpec> userAxis = {{cfgAxis.axisInvMass, "m (GeV/c^2)"}};

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    histos.add("hMultiplicity", "hMultiplicity", kTH1D, {cfgAxis.axisMultiplicity});
    histos.add("hCentrality", "hCentrality", kTH1D, {cfgAxis.axisMultiplicity});
    histos.add("hVertexZ", "hVertexZ", kTH1D, {cfgAxis.axisVtxZ});

    histos.add("hGapSide", "hGapSide", kTH1I, {{5, -0.5, 4.5}});
    histos.add("hTrueGapSide", "hTrueGapSide", kTH1I, {{6, -1.5, 4.5}});
    histos.add("hTrueGapSide_AfterSel", "hTrueGapSide_AfterSel", kTH1I, {{6, -1.5, 4.5}});

    histos.add("Trig_eta", "Trig_eta", kTH1D, {cfgAxis.axisEtaTrig});
    histos.add("Trig_phi", "Trig_phi", kTH1D, {cfgAxis.axisPhi});
    histos.add("Trig_etavsphi", "Trig_etavsphi", kTH2D, {cfgAxis.axisPhi, cfgAxis.axisEtaTrig});
    histos.add("Trig_pt", "Trig_pt", kTH1D, {cfgAxis.axisPtTrigger});
    histos.add("Trig_invMass", "Trig_invMass", kTH1D, {cfgAxis.axisInvMassQA});
    histos.add("Trig_hist", "Trig_hist", kTHnSparseF, {cfgAxis.axisSample, cfgAxis.axisVtxZ, cfgAxis.axisMultiplicity, cfgAxis.axisPtTrigger, cfgAxis.axisInvMass});
    histos.add("Trig_amp", "Trig_amp", kTH1D, {cfgAxis.axisAmplitude});
    histos.add("Channel_vs_Trig_amp", "Channel_vs_Trig_amp", kTH2D, {cfgAxis.axisChannel, cfgAxis.axisAmplitude});

    histos.add("Assoc_eta", "Assoc_eta", kTH1D, {cfgAxis.axisEtaAssoc});
    histos.add("Assoc_phi", "Assoc_phi", kTH1D, {cfgAxis.axisPhi});
    histos.add("Assoc_etavsphi", "Assoc_etavsphi", kTH2D, {cfgAxis.axisPhi, cfgAxis.axisEtaAssoc});
    histos.add("Assoc_amp", "Assoc_amp", kTH1D, {cfgAxis.axisAmplitude});
    histos.add("Channel_vs_Assoc_amp", "Channel_vs_Assoc_amp", kTH2D, {cfgAxis.axisChannel, cfgAxis.axisAmplitude});

    histos.add("deltaEta_deltaPhi_same", "deltaEta_deltaPhi_same", kTH2D, {cfgAxis.axisDeltaPhi, cfgAxis.axisDeltaEta});
    histos.add("deltaEta_deltaPhi_mixed", "deltaEta_deltaPhi_mixed", kTH2D, {cfgAxis.axisDeltaPhi, cfgAxis.axisDeltaEta});

    histos.add("TPCNClsFound", "TPCNClsFound", kTH1D, {cfgAxis.axisTPCNClsFound});
    histos.add("TPCNClsCrossedRows", "TPCNClsCrossedRows", kTH1D, {cfgAxis.axisTPCNClsCrossedRows});
    histos.add("TPCChi2NCl", "TPCChi2NCl", kTH1D, {cfgAxis.axisTPCChi2NCl});
    histos.add("TPCdcaZ", "TPCdcaZ", kTH1D, {cfgAxis.axisTPCdcaZ});

    histos.add("MFTAmbDegree", "MFTAmbDegree", kTH1D, {cfgAxis.axisMFTAmbDegree});
    histos.add("MFTNClusters", "MFTNClusters", kTH1D, {cfgAxis.axisMFTNClusters});
    histos.add("MFTbestDCAXY", "MFTbestDCAXY", kTH1D, {cfgAxis.axisMFTbestDCAXY});
    histos.add("MFTbestDCAZ", "MFTbestDCAZ", kTH1D, {cfgAxis.axisMFTbestDCAZ});
  }

  template <typename TTrack>
  bool isTrackSelected(TTrack const& track)
  {
    if constexpr (requires { track.tpcNClsFound(); }) {
      if (track.tpcNClsFound() < cfgSel.cfgTpcMinNclsFound)
        return false;
      if (track.tpcNClsCrossedRows() < cfgSel.cfgTpcMinNCrossedRows)
        return false;
      if (track.tpcChi2NCl() > cfgSel.cfgTpcMaxChi2PerCluster)
        return false;
      if (std::abs(track.dcaZ()) > cfgSel.cfgTpcMaxDcaZ)
        return false;
      return true;
    } else if constexpr (requires { track.nClusters(); }) {
      if (track.nClusters() < cfgSel.cfgMftCluster)
        return false;
      if (std::abs(track.bestDCAXY()) >= cfgSel.cfgMftDcaxy)
        return false;
      if (std::abs(track.bestDCAZ()) >= cfgSel.cfgMftDcaz)
        return false;
      if (cfgSel.cfgRejectAmbTrk && track.ambDegree() > 1)
        return false;
      if (cfgSel.cfgRejectNonAmbTrk && track.ambDegree() == 1)
        return false;
      return true;
    } else {
      return true;
    }
  }

  template <typename TCollision>
  void fillCollQA(TCollision const& col)
  {
    histos.fill(HIST("hMultiplicity"), col.multiplicity());
    if constexpr (requires { col.centrality(); }) {
      histos.fill(HIST("hCentrality"), col.centrality());
    }
    histos.fill(HIST("hVertexZ"), col.posZ());
  }

  template <typename TTrack>
  void fillTrigTrackQA(TTrack const& track)
  {
    histos.fill(HIST("Trig_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Trig_eta"), track.eta());
    histos.fill(HIST("Trig_phi"), track.phi());
    if constexpr (requires { track.channelID(); }) {
      histos.fill(HIST("Trig_amp"), track.amplitude());
      histos.fill(HIST("Channel_vs_Trig_amp"), track.channelID(), track.amplitude());
    } else {
      histos.fill(HIST("Trig_pt"), track.pt());
    }
    if constexpr (requires { track.invMass(); }) {
      histos.fill(HIST("Trig_invMass"), track.invMass());
    }
    if constexpr (requires { track.tpcNClsFound(); }) {
      histos.fill(HIST("TPCNClsFound"), track.tpcNClsFound());
      histos.fill(HIST("TPCNClsCrossedRows"), track.tpcNClsCrossedRows());
      histos.fill(HIST("TPCChi2NCl"), track.tpcChi2NCl());
      histos.fill(HIST("TPCdcaZ"), track.dcaZ());
    }
    if constexpr (requires { track.nClusters(); }) {
      histos.fill(HIST("MFTNClusters"), track.nClusters());
      histos.fill(HIST("MFTbestDCAXY"), track.bestDCAXY());
      histos.fill(HIST("MFTbestDCAZ"), track.bestDCAZ());
      histos.fill(HIST("MFTAmbDegree"), track.ambDegree());
    }
  }

  template <typename TTrack>
  void fillAssocTrackQA(TTrack const& track)
  {
    histos.fill(HIST("Assoc_etavsphi"), track.phi(), track.eta());
    histos.fill(HIST("Assoc_eta"), track.eta());
    histos.fill(HIST("Assoc_phi"), track.phi());
    if constexpr (requires { track.channelID(); }) {
      histos.fill(HIST("Assoc_amp"), track.amplitude());
      histos.fill(HIST("Channel_vs_Assoc_amp"), track.channelID(), track.amplitude());
    }
    if constexpr (requires { track.nClusters(); }) {
      histos.fill(HIST("MFTNClusters"), track.nClusters());
      histos.fill(HIST("MFTbestDCAXY"), track.bestDCAXY());
      histos.fill(HIST("MFTbestDCAZ"), track.bestDCAZ());
      histos.fill(HIST("MFTAmbDegree"), track.ambDegree());
    }
  }

  template <bool fillHist = true, typename CheckCol>
  bool isUpcEventSelected(CheckCol const& col)
  {
    if constexpr (fillHist) {
      histos.fill(HIST("hGapSide"), col.gapSide());
    }
    int truegapSide = sgSelector.trueGap(col, cfgFv0Cut, cfgFt0aCut, cfgFt0cCut, cfgZdcCut);
    if constexpr (fillHist) {
      histos.fill(HIST("hTrueGapSide"), truegapSide);
    }
    if (truegapSide != cfgGapSideCut)
      return false;
    if constexpr (fillHist) {
      histos.fill(HIST("hTrueGapSide_AfterSel"), truegapSide);
    }
    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TAssocs>
  void fillCorrHist(TTarget target, TTriggers const& triggers, TAssocs const& assocs, bool mixing, float vz, float multiplicity, float eventWeight)
  {
    int fSampleIndex = gRandom->Uniform(0, cfgSel.cfgSampleSize);
    for (auto const& triggerTrack : triggers) {
      auto trigAmpl = 1.0f;
      if constexpr (requires { triggerTrack.channelID(); }) {
        trigAmpl = triggerTrack.amplitude();
      } else {
        trigAmpl = 1.0;
      }

      if (!isTrackSelected(triggerTrack))
        continue;

      if constexpr (requires { triggerTrack.trackType(); }) {
        if (cfgSel.cfgPidMask != 0 && (cfgSel.cfgPidMask & (1u << static_cast<uint32_t>(triggerTrack.trackType()))) == 0u)
          continue;
      } else if constexpr (requires { triggerTrack.v0Type(); }) {
        if (cfgSel.cfgV0Mask != 0 && (cfgSel.cfgV0Mask & (1u << static_cast<uint32_t>(triggerTrack.v0Type()))) == 0u)
          continue;
      }
      if (!mixing) {
        fillTrigTrackQA(triggerTrack);
        if constexpr (requires { triggerTrack.channelID(); }) {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, 1.0, 1.0, eventWeight * trigAmpl);
        } else if constexpr (requires { triggerTrack.invMass(); }) {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, triggerTrack.pt(), triggerTrack.invMass(), eventWeight * trigAmpl);
        } else {
          histos.fill(HIST("Trig_hist"), fSampleIndex, vz, multiplicity, triggerTrack.pt(), 1.0, eventWeight * trigAmpl);
        }
      }
      for (auto const& assoTrack : assocs) {
        auto assoAmpl = 1.0f;
        if constexpr (requires { assoTrack.channelID(); }) {
          assoAmpl = assoTrack.amplitude();
        } else {
          assoAmpl = 1.0f;
        }

        if (!isTrackSelected(assoTrack))
          continue;

        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - assoTrack.phi(), -PIHalf);
        float deltaEta = triggerTrack.eta() - assoTrack.eta();
        if (!mixing) {
          fillAssocTrackQA(assoTrack);
          histos.fill(HIST("deltaEta_deltaPhi_same"), deltaPhi, deltaEta, eventWeight * trigAmpl * assoAmpl);
        } else {
          histos.fill(HIST("deltaEta_deltaPhi_mixed"), deltaPhi, deltaEta, eventWeight * trigAmpl * assoAmpl);
        }
        if constexpr (requires { triggerTrack.channelID(); }) {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, 1.0, deltaPhi, deltaEta, 1.0, eventWeight * trigAmpl * assoAmpl);
        } else if constexpr (requires { triggerTrack.invMass(); }) {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, triggerTrack.pt(), deltaPhi, deltaEta, triggerTrack.invMass(), eventWeight * trigAmpl * assoAmpl);
        } else {
          target->getPairHist()->Fill(step, fSampleIndex, vz, multiplicity, triggerTrack.pt(), deltaPhi, deltaEta, 1.0, eventWeight * trigAmpl * assoAmpl);
        }
      } // associated tracks
    } // trigger tracks
  } // fill correlation

  template <typename TCollision, typename TTriggers, typename TAssocs>
  void processSame(TCollision const& col, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(col.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    fillCollQA(col);
    auto multiplicity = 1.0f;
    if constexpr (requires { col.centrality(); }) {
      if (cfgSel.isUseCentEst)
        multiplicity = col.centrality();
      else
        multiplicity = col.multiplicity();
    } else {
      multiplicity = col.multiplicity();
    }
    fillCorrHist<CorrelationContainer::kCFStepReconstructed>(same, triggers, assocs, false, col.posZ(), multiplicity, 1.0);
  } // process same

  template <typename TCollision, typename... TrackTypes>
  void processMixed(TCollision const& cols, TrackTypes&&... tracks)
  {
    auto getMultiplicity = [this](auto& col) {
      if constexpr (requires { col.gapSide(); }) {
        if (!isUpcEventSelected<false>(col)) {
          return -1.0f;
        }
      } else {
        (void)this;
      }
      auto multiplicity = 1.0f;
      if constexpr (requires { col.centrality(); }) {
        if (cfgSel.isUseCentEst)
          multiplicity = col.centrality();
        else
          multiplicity = col.multiplicity();
      } else {
        multiplicity = col.multiplicity();
      }
      return multiplicity;
    };
    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {cfgAxis.axisVtxZME, cfgAxis.axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<TCollision, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgSel.cfgNmixedevent, -1, cols, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      if constexpr (requires { col1.gapSide(); } || requires { col2.gapSide(); }) {
        if (!isUpcEventSelected<false>(col1) || !isUpcEventSelected<false>(col2)) {
          continue;
        }
      }
      float eventweight = 1.0f / it.currentWindowNeighbours();
      auto multiplicity = getMultiplicity(col1);
      fillCorrHist<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
    } // pair loop
  } // process mixed

  template <typename TTriggers, typename TAssocs>
  void processMcSame(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, TTriggers const& triggers, TAssocs const& assocs)
  {
    if (std::abs(mccollision.posZ()) >= cfgSel.cfgVtxCut) {
      return;
    }
    fillCollQA(mccollision);
    auto multiplicity = mccollision.multiplicity();
    if (cfgSel.isUseDataLikeMult > 0) {
      for (const auto& collision : collisions) {
        if (cfgSel.isUseCentEst)
          multiplicity = collision.centrality();
        else
          multiplicity = collision.multiplicity();
      }
    }
    fillCorrHist<CorrelationContainer::kCFStepAll>(same, triggers, assocs, false, mccollision.posZ(), multiplicity, 1.0);
  } // process MC same

  template <typename... TrackTypes>
  void processMcMixed(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, TrackTypes&&... tracks)
  {
    bool useMCMultiplicity = (cfgSel.isUseDataLikeMult == 0);
    auto getMultiplicity =
      [&collisions, &useMCMultiplicity, this](auto& col) {
        if (useMCMultiplicity)
          return col.multiplicity();
        auto groupedCollisions = collisions.sliceByCached(aod::lrcorrcolltable::lrMcCollisionId, col.globalIndex(), this->cache);
        if (groupedCollisions.size() == 0)
          return -1.0f;
        if (cfgSel.isUseCentEst)
          return groupedCollisions.begin().centrality();
        else
          return groupedCollisions.begin().multiplicity();
      };

    using MixedBinning = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    MixedBinning binningOnVtxAndMult{{getMultiplicity}, {cfgAxis.axisVtxZME, cfgAxis.axisMultME}, true};
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TupleAtrack = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TupleBtrack = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<McCollsTable, TupleAtrack, TupleBtrack, MixedBinning> pairs{binningOnVtxAndMult, cfgSel.cfgNmixedevent, -1, mccollisions, tracksTuple, &cache};
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [col1, tracks1, col2, tracks2] = *it;
      float eventweight = 1.0f / it.currentWindowNeighbours();
      auto multiplicity = getMultiplicity(col1);
      fillCorrHist<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, true, col1.posZ(), multiplicity, eventweight);
    } // pair loop
  } // process MC mixed

  void processTpcft0aSE(CollsTable::iterator const& col, TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processTpcft0cSE(CollsTable::iterator const& col, TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processSame(col, tracks, ft0cs);
  }

  void processTpcmftSE(CollsTable::iterator const& col, TrksTable const& tracks, MftTrksTable const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processMftft0aSE(CollsTable::iterator const& col, MftTrksTable const& mfts, Ft0aTrksTable const& ft0as)
  {
    processSame(col, mfts, ft0as);
  }

  void processV0ft0aSE(CollsTable::iterator const& col, V0TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processSame(col, tracks, ft0as);
  }

  void processV0mftSE(CollsTable::iterator const& col, V0TrksTable const& tracks, MftTrksTable const& mfts)
  {
    processSame(col, tracks, mfts);
  }

  void processFt0aft0cSE(CollsTable::iterator const& col, Ft0aTrksTable const& ft0as, Ft0cTrksTable const& ft0cs)
  {
    processSame(col, ft0as, ft0cs);
  }

  void processTpcft0aME(CollsTable const& cols, TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processTpcft0cME(CollsTable const& cols, TrksTable const& tracks, Ft0cTrksTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processTpcmftME(CollsTable const& cols, TrksTable const& tracks, MftTrksTable const& mfts)
  {
    processMixed(cols, tracks, mfts);
  }

  void processMftft0aME(CollsTable const& cols, MftTrksTable const& mfts, Ft0aTrksTable const& ft0as)
  {
    processMixed(cols, mfts, ft0as);
  }

  void processV0ft0aME(CollsTable const& cols, V0TrksTable const& tracks, Ft0aTrksTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processV0mftME(CollsTable const& cols, V0TrksTable const& tracks, MftTrksTable const& mfts)
  {
    processMixed(cols, tracks, mfts);
  }

  void processFt0aft0cME(CollsTable const& cols, Ft0aTrksTable const& ft0as, Ft0cTrksTable const& ft0cs)
  {
    processMixed(cols, ft0as, ft0cs);
  }

  void processUpcTpcft0aSE(UpcCollsTable::iterator const& col, TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0as);
  }

  void processUpcTpcft0cSE(UpcCollsTable::iterator const& col, TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0cs);
  }

  void processUpcTpcmftSE(UpcCollsTable::iterator const& col, TrksUpcTable const& tracks, MftTrksUpcTable const& mfts)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, mfts);
  }

  void processUpcMftft0aSE(UpcCollsTable::iterator const& col, MftTrksUpcTable const& mfts, Ft0aTrksUpcTable const& ft0as)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, mfts, ft0as);
  }

  void processUpcV0ft0aSE(UpcCollsTable::iterator const& col, V0TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, ft0as);
  }

  void processUpcV0mftSE(UpcCollsTable::iterator const& col, V0TrksUpcTable const& tracks, MftTrksUpcTable const& mfts)
  {
    if (!isUpcEventSelected<true>(col)) {
      return;
    }
    processSame(col, tracks, mfts);
  }

  void processUpcTpcft0aME(UpcCollsTable const& cols, TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processUpcTpcft0cME(UpcCollsTable const& cols, TrksUpcTable const& tracks, Ft0cTrksUpcTable const& ft0cs)
  {
    processMixed(cols, tracks, ft0cs);
  }

  void processUpcTpcmftME(UpcCollsTable const& cols, TrksUpcTable const& tracks, MftTrksUpcTable const& mfts)
  {
    processMixed(cols, tracks, mfts);
  }

  void processUpcMftft0aME(UpcCollsTable const& cols, MftTrksUpcTable const& mfts, Ft0aTrksUpcTable const& ft0as)
  {
    processMixed(cols, mfts, ft0as);
  }

  void processUpcV0ft0aME(UpcCollsTable const& cols, V0TrksUpcTable const& tracks, Ft0aTrksUpcTable const& ft0as)
  {
    processMixed(cols, tracks, ft0as);
  }

  void processUpcV0mftME(UpcCollsTable const& cols, V0TrksUpcTable const& tracks, MftTrksUpcTable const& mfts)
  {
    processMixed(cols, tracks, mfts);
  }

  void processMcTpcft0aSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcSame(mccollision, collisions, tracks, ft0as);
  }

  void processMcTpcft0cSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcSame(mccollision, collisions, tracks, ft0cs);
  }

  void processMcTpcmftSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McTrksTable const& tracks, McMftTrksTable const& mfts)
  {
    processMcSame(mccollision, collisions, tracks, mfts);
  }

  void processMcMftft0aSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McMftTrksTable const& mfts, McFt0aTrksTable const& ft0as)
  {
    processMcSame(mccollision, collisions, mfts, ft0as);
  }

  void processMcFt0aft0cSE(McCollsTable::iterator const& mccollision, soa::SmallGroups<aod::LRCollisionsWithLabel> const& collisions, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcSame(mccollision, collisions, ft0as, ft0cs);
  }

  void processMcTpcft0aME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McTrksTable const& tracks, McFt0aTrksTable const& ft0as)
  {
    processMcMixed(mccollisions, collisions, tracks, ft0as);
  }

  void processMcTpcft0cME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McTrksTable const& tracks, McFt0cTrksTable const& ft0cs)
  {
    processMcMixed(mccollisions, collisions, tracks, ft0cs);
  }

  void processMcTpcmftME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McTrksTable const& tracks, McMftTrksTable const& mfts)
  {
    processMcMixed(mccollisions, collisions, tracks, mfts);
  }

  void processMcMftft0aME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McMftTrksTable const& mfts, McFt0aTrksTable const& ft0as)
  {
    processMcMixed(mccollisions, collisions, mfts, ft0as);
  }

  void processMcFt0aft0cME(McCollsTable const& mccollisions, aod::LRCollisionsWithLabel const& collisions, McFt0aTrksTable const& ft0as, McFt0cTrksTable const& ft0cs)
  {
    processMcMixed(mccollisions, collisions, ft0as, ft0cs);
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
  PROCESS_SWITCH(LongrangecorrDerived, processFt0aft0cSE, "same event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processFt0aft0cME, "mixed event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0aSE, "same UPC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0aME, "mixed UPC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0cSE, "same UPC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcft0cME, "mixed UPC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcmftSE, "same UPC event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcTpcmftME, "mixed UPC event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcMftft0aSE, "same UPC event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcMftft0aME, "mixed UPC event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0aSE, "same UPC event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0ft0aME, "mixed UPC event V0 vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0mftSE, "same UPC event V0 vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processUpcV0mftME, "mixed UPC event V0 vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0aSE, "same MC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0aME, "mixed MC event TPC vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0cSE, "same MC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcft0cME, "mixed MC event TPC vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcmftSE, "same MC event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcTpcmftME, "mixed MC event TPC vs MFT", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcMftft0aSE, "same MC event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcMftft0aME, "mixed MC event MFT vs FT0A", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcFt0aft0cSE, "same MC event FT0A vs FT0C", false);
  PROCESS_SWITCH(LongrangecorrDerived, processMcFt0aft0cME, "mixed MC event FT0A vs FT0C", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LongrangecorrDerived>(cfgc)};
}
