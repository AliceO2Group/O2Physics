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

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr float cfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
constexpr float kThreeHalfPi = 1.5f * PI;

struct TwoParticleCorrelations {

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
	// Declare configurables for correlations
	Configurable<float> cfgZVtxCut = {"zvtxcut", 10.0, "Vertex z cut. Default 10 cm"};
	Configurable<float> cfgPtCutMin = {"minpt", 0.2, "Minimum accepted track pT. Default 0.2 GeV"};
	Configurable<float> cfgPtCutMax = {"maxpt", 3.0, "Maximum accepted track pT. Default 3.0 GeV"};
	Configurable<float> cfgEtaCut = {"etacut", 0.8, "Eta cut. Default 0.8"};
	Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut",
	                                            {cfgPairCutDefaults[0],
	                                            5,
	                                            {"Photon", "K0", "Lambda", "Phi", "Rho"}},
	                                            "Pair cuts on various particles"};
	//Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, {"Two track cut"}};
	ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
	ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {32, -constants::math::PIHalf, constants::math::PIHalf * 3}, "delta phi axis for histograms"};
	ConfigurableAxis axisDeltaEta{"axisDeltaEta", {32, -1.6, 1.6}, "delta eta axis for histograms"};
	ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
	ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
	ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1},"multiplicity / centrality axis for histograms"};
	ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
	ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
	ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25,2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

	// Output definitions
	OutputObj<CorrelationContainer> same{"sameEvent"};
	OutputObj<CorrelationContainer> mixed{"mixedEvent"};
	HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
	PairCuts mPairCuts;
	bool doPairCuts = false;

  void init(InitContext&)
  {
    LOGF(info, "Starting init");

    histos.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    histos.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, 2 * M_PI, "#varphi"}}});
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

  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgZVtxCut;
  Filter trackFilter = (nabs(aod::track::eta) < cfgEtaCut) && (aod::track::pt > cfgPtCutMin) && (aod::track::pt < cfgPtCutMax) && (requireGlobalTrackInFilter() || (aod::track::isGlobalTrackSDD == (uint8_t) true));

  struct SameEventTag {
  };
  struct MixedEventTag {
  };

  template <typename TTracks>
  void fillQA(TTracks tracks, float multiplicity)
  {
  	for (auto& track: tracks) {
  	  histos.fill(HIST("yields"), multiplicity, track.pt(), track.eta());
  	  histos.fill(HIST("etaphi"), multiplicity, track.eta(), track.phi());
  	}
  }

  template <typename TTarget, typename TCollision>
  bool fillCollision(TTarget target, TCollision collision, float multiplicity)
  {
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    if (!collision.alias_bit(kINT7) || !collision.sel7()){
      return false;
    }
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    return true;
  }

  template <typename TTarget, typename TTracks, typename TTag>
  void fillCorrelations(TTarget target, TTracks tracks1, TTracks tracks2, float multiplicity, float posZ)
  {
  	for (const auto& track1 : tracks1) {
  	  if (isTrackCut(track1) == false) {
  	    return;
  	  }
  	  float phi1 = phi(track1.px(), track1.py());
  	  phi1 = RecoDecay::constrainAngle(phi1, 0.f);
  	  float eta1 = eta(track1.px(), track1.py(), track1.pz());
  	  target->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, posZ, 1.0);
  	  for (const auto& track2 : tracks2) {
  	    if (track1 == track2) {
  	      continue;
  	    }
  	    if (isTrackCut(track2) == false) {
  	      return;
  	    }
  	    float phi2 = phi(track2.px(), track2.py());
  	    phi2 = RecoDecay::constrainAngle(phi2, 0.f);
  	    float eta2 = eta(track2.px(), track2.py(), track2.pz());
  	    if (doPairCuts && mPairCuts.conversionCuts(track1, track2)) {
  	      continue;
  	    }
  	    float deltaPhi = phi1 - phi2;
  	    float deltaEta = eta1 - eta2;
  	    deltaPhi = RecoDecay::constrainAngle(deltaPhi, -PIHalf);
  	    target->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
  	                                deltaEta,
  	                                track2.pt(), track1.pt(),
  	                                multiplicity,
  	                                deltaPhi,
  	                                posZ);
  	    if constexpr (std::is_same_v<TTag, SameEventTag>) {
  	      if (minMultiplicity <= multiplicity) {
  	        histos.fill(HIST("sameEvent2D"), deltaEta, deltaPhi);
  	      }
  	      if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
  	        histos.fill(HIST("sameEvent_2_10"), deltaEta, deltaPhi);
  	      }
  	      if (range2Min <= multiplicity && multiplicity <= range2Max) {
  	        histos.fill(HIST("sameEvent_11_20"), deltaEta, deltaPhi);
  	      }
  	      if (range3Min <= multiplicity && multiplicity <= range3Max) {
  	        histos.fill(HIST("sameEvent_21_30"), deltaEta, deltaPhi);
  	      }
  	      if (range4Min <= multiplicity && multiplicity <= range4Max) {
  	        histos.fill(HIST("sameEvent_31_40"), deltaEta, deltaPhi);
  	      }
  	      if (range5Min <= multiplicity && multiplicity <= range5Max) {
  	        histos.fill(HIST("sameEvent_41_50"), deltaEta, deltaPhi);
  	      }
  	    } else if constexpr (std::is_same_v<TTag, MixedEventTag>) {
  	      if (minMultiplicity <= multiplicity) {
  	        histos.fill(HIST("mixedEvent2D"), deltaEta, deltaPhi);
  	      }
  	      if (minMultiplicity <= multiplicity && multiplicity <= range1Max) {
  	        histos.fill(HIST("mixedEvent_2_10"), deltaEta, deltaPhi);
  	      }
  	      if (range2Min <= multiplicity && multiplicity <= range2Max) {
  	        histos.fill(HIST("mixedEvent_11_20"), deltaEta, deltaPhi);
  	      }
  	      if (range3Min <= multiplicity && multiplicity <= range3Max) {
  	        histos.fill(HIST("mixedEvent_21_30"), deltaEta, deltaPhi);
  	      }
  	      if (range4Min <= multiplicity && multiplicity <= range4Max) {
  	        histos.fill(HIST("mixedEvent_31_40"), deltaEta, deltaPhi);
  	      }
  	      if (range5Min <= multiplicity && multiplicity <= range5Max) {
  	        histos.fill(HIST("mixedEvent_41_50"), deltaEta, deltaPhi);
  	      }
  	    }
  	  }
  	}
  }

  void processSame()
  {
  }

  void processMixed()
  {
  }
  void process(aod::Collision const&, aod::Tracks const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TwoParticleCorrelations>(cfgc),
  };
}