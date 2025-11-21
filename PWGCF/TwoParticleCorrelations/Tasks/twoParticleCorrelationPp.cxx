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

/// \file   twoParticleCorrelationPp.cxx
/// \brief  Task for two particle correlation in pp in order to calculate a baseline in a template fit for the flow coefficients
/// \author Josué Martínez García <josuem@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
constexpr float kThreeHalfPi = 1.5f * PI;

struct TwoParticleCorrelationPp {

  // Declare configurables on events/collisions
  Configurable<int> minMultiplicity{"minMultiplicity", 2, {"Range on multiplicity"}};
  Configurable<int> range1Max{"range1Max", 10, {"Range on multiplicity"}};
  Configurable<int> range2Min{"range2Min", 11, {"Range on multiplicity"}};
  Configurable<int> range2Max{"range2Max", 20, {"Range on multiplicity"}};
  Configurable<int> range3Min{"range3Min", 21, {"Range on multiplicity"}};
  Configurable<int> range3Max{"range3Max", 30, {"Range on multiplicity"}};
  Configurable<int> range4Min{"range4Min", 31, {"Range on multiplicity"}};
  Configurable<int> range4Max{"range4Max", 40, {"Range on multiplicity"}};
  Configurable<int> range5Min{"range5Min", 41, {"Range on multiplicity"}};
  Configurable<int> range5Max{"range5Max", 50, {"Range on multiplicity"}};
  Configurable<int> nEventsMixed{"nEventsMixed", 5, {"Events to be Mixed"}};
  Configurable<bool> evSel8{"evSel8", true, {"rejects collisions using sel8()"}};
  Configurable<float> cfgZVtxCut = {"cfgZVtxCut", 10.0, "Vertex z cut. Default 10 cm"};
  Configurable<bool> evSelkNoSameBunchPileup{"evSelkNoSameBunchPileup", true, {"rejects collisions which are associated with the same found-by-T0 bunch crossing"}};
  Configurable<bool> evSelkNoITSROFrameBorder{"evSelkNoITSROFrameBorder", true, {"reject events at ITS ROF border"}};
  Configurable<bool> evSelkNoTimeFrameBorder{"evSelkNoTimeFrameBorder", true, {"reject events at TF border"}};
  Configurable<bool> evSelkIsGoodZvtxFT0vsPV{"evSelkIsGoodZvtxFT0vsPV", true, {"removes collisions with large differences between z of PV by tracks and z of PV from FT0 A-C time difference, use this cut at low multiplicities with caution"}};
  Configurable<bool> evSelkNoCollInTimeRangeStandard{"evSelkNoCollInTimeRangeStandard", true, {"no collisions in specified time range"}};
  // Declare configurables on tracks
  Configurable<float> cutMyptMin{"cutMyptMin", 0.2, {"My Track cut"}};
  Configurable<float> cutMyptMax{"cutMyptMax", 3., {"My Track cut"}};
  Configurable<float> cutMyetaMin{"cutMyetaMin", -0.8, {"My Track cut"}};
  Configurable<float> cutMyetaMax{"cutMyetaMax", 0.8, {"My Track cut"}};
  Configurable<float> cutMydcaZmax{"cutMydcaZmax", 2.f, {"My Track cut"}};
  Configurable<float> cutMydcaXYmax{"cutMydcaXYmax", 1e0f, {"My Track cut"}};
  Configurable<bool> cutMydcaXYusePt{"cutMydcaXYusePt", false, {"My Track cut"}};
  Configurable<bool> cutMyHasITS{"cutMyHasITS", true, {"My Track cut"}};
  Configurable<int> cutMyITSNClsMin{"cutMyITSNClsMin", 1, {"My Track cut"}};
  Configurable<float> cutMyITSChi2NClMax{"cutMyITSChi2NClMax", 36.f, {"My Track cut"}};
  Configurable<bool> cutMyHasTPC{"cutMyHasTPC", true, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyTPCNClsCrossedRowsMin{"cutMyTPCNClsCrossedRowsMin", 70, {"My Track cut"}};
  Configurable<int> cutMyTPCNClsFindableMin{"cutMyTPCNClsFindableMin", 50, {"My Track cut"}};
  Configurable<int> cutMyTPCNClsMin{"cutMyTPCNClsMin", 1, {"My Track cut"}};
  Configurable<float> cutMyTPCNClsCrossedRowsOverNClsFindableMin{"cutMyTPCNClsCrossedRowsOverNClsFindableMin", 0.8f, {"My Track cut"}};
  Configurable<float> cutMyTPCNClsOverFindableNClsMin{"cutMyTPCNClsOverFindableNClsMin", 0.5f, {"My Track cut"}};
  Configurable<float> cutMyTPCChi2NclMax{"cutMyTPCChi2NclMax", 4.f, {"My Track cut"}};
  // Declare configurables for correlations

  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut",
                                               {cfgPairCutDefaults[0],
                                                5,
                                                {"Photon", "K0", "Lambda", "Phi", "Rho"}},
                                               "Pair cuts on various particles"};
  // Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, {"Two track cut"}};
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {32, -PIHalf, kThreeHalfPi}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {32, -1.6, 1.6}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 110.1}, "multiplicity / multiplicity axis for histograms"};
  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  PairCuts mPairCuts;
  bool doPairCuts = false;

  void init(InitContext&)
  {
    LOGF(info, "Starting init");

    const AxisSpec axisCountTracks{17, -0.5, 16.5};
    const AxisSpec axisCountEvents{8, -0.5, 7.5};

    histos.add("yields", "multiplicity vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "multiplicity"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    histos.add("etaphi", "multiplicity vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity"}, {100, -2, 2, "#eta"}, {64, 0., TwoPI, "#varphi"}}});
    histos.add("sameEvent2D", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_2_10", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_11_20", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_21_30", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_31_40", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("sameEvent_41_50", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent2D", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_2_10", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_11_20", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_21_30", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_31_40", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("mixedEvent_41_50", "#Delta #eta vs #Delta #phi", {HistType::kTH2F, {axisDeltaEta, axisDeltaPhi}});
    histos.add("Tracks/hTracksAfterCuts", " ; ; counts", kTH1F, {axisCountTracks});
    histos.add("Events/hEventsAfterCuts", " ; ; counts", kTH1F, {axisCountEvents});

    const int maxMixBin = axisMultiplicity->size() * axisVertex->size();
    histos.add("eventcount", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    mPairCuts.SetHistogramRegistry(&histos);
    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 ||
        cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      doPairCuts = true;
    }

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, {}));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, {}));

    LOGF(info, "Finishing init");
  }

  // Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgZVtxCut;
  // Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && (requireGlobalTrackInFilter() || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  using FullCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>;

  enum EventType {
    SameEvent = 1,
    MixedEvent = 2
  };

  SliceCache cache;
  std::vector<double> vtxBinsEdges{VARIABLE_WIDTH, -10.0f, -7.0f, -5.0f, -2.5f, 0.0f, 2.5f, 5.0f, 7.0f, 10.0f};
  using BinningType = ColumnBinningPolicy<aod::collision::PosZ>;

  template <typename TCollision>
  bool isEventSelected(TCollision collision)
  {
    if (evSel8 && !collision.sel8()) {
      return false;
    }
    if (collision.posZ() < -cfgZVtxCut || cfgZVtxCut < collision.posZ()) {
      return false;
    }
    if (evSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (evSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (evSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (evSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (evSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isTrackCut(T const& track)
  {
    if (track.sign() != 1 && track.sign() != -1) {
      return false;
    }
    if (track.pt() < cutMyptMin || track.pt() > cutMyptMax) {
      return false;
    }
    if (track.eta() < cutMyetaMin || track.eta() > cutMyetaMax) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cutMydcaZmax) {
      return false;
    }
    if (cutMydcaXYusePt) {
      float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
      if (std::abs(track.dcaXY()) > maxDCA) {
        return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > cutMydcaXYmax) {
        return false;
      }
    }
    if (track.isPVContributor() == false) {
      return false;
    }
    // Quality Track
    // ITS
    if (cutMyHasITS && !track.hasITS()) {
      return false; // ITS refit
    }
    if (track.itsNCls() < cutMyITSNClsMin) {
      return false;
    }
    if (track.itsChi2NCl() > cutMyITSChi2NClMax) {
      return false;
    }
    // TPC
    if (cutMyHasTPC && !track.hasTPC()) {
      return false; // TPC refit
    }
    if (track.tpcNClsCrossedRows() < cutMyTPCNClsCrossedRowsMin) {
      return false;
    }
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutMyTPCNClsMin) {
      return false; // tpcNClsFound()
    }
    if (track.tpcNClsFindable() < cutMyTPCNClsFindableMin) {
      return false;
    }
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
      return false; //
    }
    if ((static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
      return false; //
    }
    if (track.tpcChi2NCl() > cutMyTPCChi2NclMax) {
      return false; // TPC chi2
    }
    return true;
  }

  template <typename TTracks>
  void fillQA(TTracks tracks, float multiplicity)
  {
    for (const auto& track : tracks) {
      histos.fill(HIST("yields"), multiplicity, track.pt(), track.eta());
      histos.fill(HIST("etaphi"), multiplicity, track.eta(), track.phi());
    }
  }

  template <typename TTarget>
  bool fillCollision(TTarget target, float multiplicity)
  {
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    return true;
  }

  template <typename TTarget, typename TTracks>
  void fillCorrelations(TTarget target, TTracks tracks1, TTracks tracks2, float multiplicity, float posZ, int system)
  {
    for (const auto& track1 : tracks1) {
      if (isTrackCut(track1) == false) {
        return;
      }
      float phi1 = track1.phi();
      phi1 = RecoDecay::constrainAngle(phi1, 0.f);
      float eta1 = track1.eta();
      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, posZ, 1.0);
      for (const auto& track2 : tracks2) {
        if (track1 == track2) {
          continue;
        }
        if (isTrackCut(track2) == false) {
          return;
        }
        float phi2 = track2.phi();
        phi2 = RecoDecay::constrainAngle(phi2, 0.f);
        float eta2 = track2.eta();
        /*if (doPairCuts && mPairCuts.conversionCuts(track1, track2)) {
          continue;
        }*/
        float deltaPhi = phi1 - phi2;
        float deltaEta = eta1 - eta2;
        deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
        target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                    deltaEta,
                                    track2.pt(), track1.pt(),
                                    multiplicity,
                                    deltaPhi,
                                    posZ);

        if (system == SameEvent) {
          if (minMultiplicity <= multiplicity) {
            histos.fill(HIST("sameEvent2D"), deltaEta, deltaPhi);
          }
          if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
            histos.fill(HIST("sameEvent_2_10"), deltaEta, deltaPhi);
          } else if (range2Min <= multiplicity && multiplicity <= range2Max) {
            histos.fill(HIST("sameEvent_11_20"), deltaEta, deltaPhi);
          } else if (range3Min <= multiplicity && multiplicity <= range3Max) {
            histos.fill(HIST("sameEvent_21_30"), deltaEta, deltaPhi);
          } else if (range4Min <= multiplicity && multiplicity <= range4Max) {
            histos.fill(HIST("sameEvent_31_40"), deltaEta, deltaPhi);
          } else if (range5Min <= multiplicity && multiplicity <= range5Max) {
            histos.fill(HIST("sameEvent_41_50"), deltaEta, deltaPhi);
          } else {
            continue;
          }
        } else if (system == MixedEvent) {
          if (minMultiplicity <= multiplicity) {
            histos.fill(HIST("mixedEvent2D"), deltaEta, deltaPhi);
          }
          if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
            histos.fill(HIST("mixedEvent_2_10"), deltaEta, deltaPhi);
          } else if (range2Min <= multiplicity && multiplicity <= range2Max) {
            histos.fill(HIST("mixedEvent_11_20"), deltaEta, deltaPhi);
          } else if (range3Min <= multiplicity && multiplicity <= range3Max) {
            histos.fill(HIST("mixedEvent_21_30"), deltaEta, deltaPhi);
          } else if (range4Min <= multiplicity && multiplicity <= range4Max) {
            histos.fill(HIST("mixedEvent_31_40"), deltaEta, deltaPhi);
          } else if (range5Min <= multiplicity && multiplicity <= range5Max) {
            histos.fill(HIST("mixedEvent_41_50"), deltaEta, deltaPhi);
          } else {
            continue;
          }
        }
      }
    }
  }

  void processMixed(FullCollisions const& collisions, FullTracks const& tracks)
  {
    BinningType bindingOnVtx{{vtxBinsEdges}, true};
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<FullCollisions,
                 FullTracks,
                 BinningType>
      pairs{bindingOnVtx, nEventsMixed, -1, collisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (collision1.size() == 0 || collision2.size() == 0) {
        continue;
      }

      if (isEventSelected(collision1) == false || isEventSelected(collision2) == false) {
        continue;
      }

      float multiplicity = 0;
      for (const auto& track : tracks1) {
        if (isTrackCut(track) == false) {
          continue;
        }
        ++multiplicity;
      }
      if (fillCollision(mixed, multiplicity) == false) {
        return;
      }
      // histos.fill(HIST("eventcount"), bindingOnVtx.getBin({collision1.posZ()}));
      fillCorrelations(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), MixedEvent);
    }
  }

  PROCESS_SWITCH(TwoParticleCorrelationPp, processMixed, "Process mixed events", true);

  void processSame(FullCollisions::iterator const& collision, FullTracks const& tracks)
  {
    float multiplicity = 0;

    if (isEventSelected(collision) == false) {
      return;
    }

    for (const auto& track : tracks) {
      if (isTrackCut(track) == false) {
        continue;
      }
      ++multiplicity;
    }
    // LOGF(debug, "Filling same events");
    if (fillCollision(same, multiplicity) == false) {
      return;
    }
    histos.fill(HIST("eventcount"), -2);
    if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
      histos.fill(HIST("eventcount"), 1);
    } else if (range2Min <= multiplicity && multiplicity <= range2Max) {
      histos.fill(HIST("eventcount"), 2);
    } else if (range3Min <= multiplicity && multiplicity <= range3Max) {
      histos.fill(HIST("eventcount"), 3);
    } else if (range4Min <= multiplicity && multiplicity <= range4Max) {
      histos.fill(HIST("eventcount"), 4);
    } else if (range5Min <= multiplicity && multiplicity <= range5Max) {
      histos.fill(HIST("eventcount"), 5);
    }
    fillQA(tracks, multiplicity);
    fillCorrelations(same, tracks, tracks, multiplicity, collision.posZ(), SameEvent);
  }

  PROCESS_SWITCH(TwoParticleCorrelationPp, processSame, "Process same events", true);

  void process(FullCollisions::iterator const& collision, FullTracks const& tracks)
  {
    // Configure events flow histogram labels
    auto hFlowEvents = histos.get<TH1>(HIST("Events/hEventsAfterCuts"));
    hFlowEvents->GetXaxis()->SetBinLabel(1, "All tracks");
    hFlowEvents->GetXaxis()->SetBinLabel(2, "sel8");
    hFlowEvents->GetXaxis()->SetBinLabel(3, "Z vertex");
    hFlowEvents->GetXaxis()->SetBinLabel(4, "kNoSameBunchPileup");
    hFlowEvents->GetXaxis()->SetBinLabel(5, "kNoITSROFrameBorder");
    hFlowEvents->GetXaxis()->SetBinLabel(6, "kNoTimeFrameBorder");
    hFlowEvents->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
    hFlowEvents->GetXaxis()->SetBinLabel(8, "kNoCollInTimeRangeStandard");

    histos.fill(HIST("Events/hEventsAfterCuts"), 0);
    if (evSel8 && !collision.sel8()) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 1);
    if (collision.posZ() < -cfgZVtxCut || cfgZVtxCut < collision.posZ()) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 2);
    if (evSelkNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 3);
    if (evSelkNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 4);
    if (evSelkNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 5);
    if (evSelkIsGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 6);
    if (evSelkNoCollInTimeRangeStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    histos.fill(HIST("Events/hEventsAfterCuts"), 7);

    // Configure track flow histogram labels
    auto hFlowTracks = histos.get<TH1>(HIST("Tracks/hTracksAfterCuts"));
    hFlowTracks->GetXaxis()->SetBinLabel(1, "All tracks");
    hFlowTracks->GetXaxis()->SetBinLabel(2, "Track sign");
    hFlowTracks->GetXaxis()->SetBinLabel(3, "p_{T} range");
    hFlowTracks->GetXaxis()->SetBinLabel(4, "#eta range");
    hFlowTracks->GetXaxis()->SetBinLabel(5, "dcaZ");
    hFlowTracks->GetXaxis()->SetBinLabel(6, "dcaXY");
    hFlowTracks->GetXaxis()->SetBinLabel(7, "PV contrib cut");
    hFlowTracks->GetXaxis()->SetBinLabel(8, "has ITS cut");
    hFlowTracks->GetXaxis()->SetBinLabel(9, "N clusters ITS cut");
    hFlowTracks->GetXaxis()->SetBinLabel(10, "#chi^{2} N cluster ITS cut");
    hFlowTracks->GetXaxis()->SetBinLabel(11, "has TPC cut");
    hFlowTracks->GetXaxis()->SetBinLabel(12, "N clusters crossed row TPC cut");
    hFlowTracks->GetXaxis()->SetBinLabel(13, "(N cluster findable - N cluster minus findable) TPC cut");
    hFlowTracks->GetXaxis()->SetBinLabel(14, "N cluster findable TPC cut");
    hFlowTracks->GetXaxis()->SetBinLabel(15, "(N cluster crossed row / N cluster findable) TPC cut");
    hFlowTracks->GetXaxis()->SetBinLabel(16, "(N cluster findable - N cluster minus findable) / N cluster findable cut");
    hFlowTracks->GetXaxis()->SetBinLabel(17, "#chi^{2} N cluster TPC cut");

    for (const auto& track : tracks) {
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 0);
      if (track.sign() != 1 && track.sign() != -1) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 1);
      if (track.pt() < cutMyptMin || track.pt() > cutMyptMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 2);
      if (track.eta() < cutMyetaMin || track.eta() > cutMyetaMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 3);
      if (std::abs(track.dcaZ()) > cutMydcaZmax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 4);
      if (cutMydcaXYusePt) {
        float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
        if (std::abs(track.dcaXY()) > maxDCA) {
          continue;
        }
      } else {
        if (std::abs(track.dcaXY()) > cutMydcaXYmax) {
          continue;
        }
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 5);
      if (track.isPVContributor() == false) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 6);
      // Quality Track
      // ITS
      if (cutMyHasITS && !track.hasITS()) {
        continue; // ITS refit
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 7);
      if (track.itsNCls() < cutMyITSNClsMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 8);
      if (track.itsChi2NCl() > cutMyITSChi2NClMax) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 9);
      // TPC
      if (cutMyHasTPC && !track.hasTPC()) {
        continue; // TPC refit
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 10);
      if (track.tpcNClsCrossedRows() < cutMyTPCNClsCrossedRowsMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 11);
      if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutMyTPCNClsMin) {
        continue; // tpcNClsFound()
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 12);
      if (track.tpcNClsFindable() < cutMyTPCNClsFindableMin) {
        continue;
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 13);
      if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
        continue; //
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 14);
      if ((static_cast<float>(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) / static_cast<float>(track.tpcNClsFindable())) < cutMyTPCNClsCrossedRowsOverNClsFindableMin) {
        continue; //
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 15);
      if (track.tpcChi2NCl() > cutMyTPCChi2NclMax) {
        continue; // TPC chi2
      }
      histos.fill(HIST("Tracks/hTracksAfterCuts"), 16);
    }
  }

  PROCESS_SWITCH(TwoParticleCorrelationPp, process, "Process mixed events", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TwoParticleCorrelationPp>(cfgc),
  };
}
