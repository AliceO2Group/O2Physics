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
/// \file    flowJSPCAnalysis.cxx
/// \brief   Task for the calculation of SPC with filtered data.
/// \author  Maxim Virta (maxim.virta@cern.ch), Cindy Mordasini (cindy.mordasini@cern.ch), Neelkamal Mallick (neelkamal.mallick@cern.ch)

// Standard headers.
#include <TFormula.h>
#include <TRandom3.h>

#include <chrono>
#include <memory>
#include <string>
#include <vector>

// O2 headers. //
// The first two are mandatory.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

// O2 Physics headers. //
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"
#include "PWGCF/JCorran/Core/FlowJSPCAnalysis.h"
#include "PWGCF/JCorran/Core/FlowJSPCObservables.h"
#include "PWGCF/JCorran/Core/FlowJHistManager.h"

// Namespaces and definitions.
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                               aod::FT0sCorrected, aod::CentFT0Ms,
                               aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                               aod::CentFDDMs, aod::CentNTPVs>;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::JWeights>;

struct flowJSPCAnalysis {
  HistogramRegistry spcHistograms{"SPCResults", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  FlowJSPCAnalysis spcAnalysis;
  FlowJSPCAnalysis::JQVectorsT jqvecs;
  template <class T>
  using HasWeightNUA = decltype(std::declval<T&>().weightNUA());
  template <class T>
  using HasWeightEff = decltype(std::declval<T&>().weightEff());
  template <class T>
  using HasTrackType = decltype(std::declval<T&>().trackType());

  HistogramRegistry qaHistRegistry{"qaHistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  FlowJHistManager histManager;

  FlowJSPCObservables spcObservables;

  // Set Configurables here
  Configurable<bool> cfgFillQA{"cfgFillQA", true, "Fill QA plots"};

  Configurable<int> cfgWhichSPC{"cfgWhichSPC", 0, "Which SPC observables to compute."};

  struct : ConfigurableGroup {
    Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "Minimum pT used for track selection."};
    Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "Maximum pT used for track selection."};
    Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum eta used for track selection."};
  } cfgTrackCuts;

  struct : ConfigurableGroup {
    Configurable<int> cfgCentEst{"cfgCentEst", 2, "Centrality estimator."};
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.0f, "Maximum primary vertex cut applied for the events."};
    Configurable<int> cfgMultMin{"cfgMultMin", 10, "Minimum number of particles required for the event to have."};
  } cfgEventCuts;

  O2_DEFINE_CONFIGURABLE(cfgTrackBitMask, uint16_t, 0, "Track selection bitmask to use as defined in the filterCorrelations.cxx task");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrelationsMask, uint16_t, 0, "Selection bitmask for the multiplicity correlations. This should match the filter selection cfgEstimatorBitMask.")
  O2_DEFINE_CONFIGURABLE(cfgMultCutFormula, std::string, "", "Multiplicity correlations cut formula. A result greater than zero results in accepted event. Parameters: [cFT0C] FT0C centrality, [mFV0A] V0A multiplicity, [mGlob] global track multiplicity, [mPV] PV track multiplicity, [cFT0M] FT0M centrality")

  // // Filters to be applied to the received data.
  // // The analysis assumes the data has been subjected to a QA of its selection,
  // // and thus only the final distributions of the data for analysis are saved.
  Filter collFilter = (nabs(aod::collision::posZ) < cfgEventCuts.cfgZvtxMax);

  Filter trackFilter = (aod::track::pt > cfgTrackCuts.cfgPtMin) && (aod::track::pt < cfgTrackCuts.cfgPtMax) && (nabs(aod::track::eta) < cfgTrackCuts.cfgEtaMax);
  Filter cftrackFilter = (nabs(aod::cftrack::eta) < cfgTrackCuts.cfgEtaMax) && (aod::cftrack::pt > cfgTrackCuts.cfgPtMin) && (aod::cftrack::pt < cfgTrackCuts.cfgPtMax) && ncheckbit(aod::track::trackType, as<uint8_t>(cfgTrackBitMask));

  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint, aod::cfmultset::NMultiplicityEstimators> multCutFormulaParamIndex;

  void init(InitContext const&)
  {
    // Add histomanager here
    spcAnalysis.setHistRegistry(&spcHistograms);
    spcAnalysis.createHistos();

    spcObservables.setSPCObservables(cfgWhichSPC);
    spcAnalysis.setFullCorrSet(spcObservables.harmonicArray);

    histManager.setHistRegistryQA(&qaHistRegistry);
    histManager.setDebugLog(false);
    histManager.createHistQA();

    if (!cfgMultCutFormula.value.empty()) {
      multCutFormula = std::make_unique<TFormula>("multCutFormula", cfgMultCutFormula.value.c_str());
      std::fill_n(multCutFormulaParamIndex.begin(), std::size(multCutFormulaParamIndex), ~0u);
      std::array<std::string, aod::cfmultset::NMultiplicityEstimators> pars = {"cFT0C", "mFV0A", "mPV", "mGlob", "cFT0M"}; // must correspond the order of MultiplicityEstimators
      for (uint i = 0, n = multCutFormula->GetNpar(); i < n; ++i) {
        auto m = std::find(pars.begin(), pars.end(), multCutFormula->GetParName(i));
        if (m == pars.end()) {
          LOGF(warning, "Unknown parameter in cfgMultCutFormula: %s", multCutFormula->GetParName(i));
          continue;
        }
        const uint estIdx = std::distance(pars.begin(), m);
        if ((cfgMultCorrelationsMask.value & (1u << estIdx)) == 0) {
          LOGF(warning, "The centrality/multiplicity estimator %s is not available to be used in cfgMultCutFormula. Ensure cfgMultCorrelationsMask is correct and matches the CFMultSets in derived data.", m->c_str());
        } else {
          multCutFormulaParamIndex[estIdx] = i;
          LOGF(info, "Multiplicity cut parameter %s in use.", m->c_str());
        }
      }
    }
  }

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  // void process(soa::Filtered<MyCollisions>::iterator const& coll, soa::Filtered<soa::Join<aod::MyTracks, aod::JWeights>> const& tracks)
  {
    if (tracks.size() < cfgEventCuts.cfgMultMin)
      return;

    float cent = collision.multiplicity();
    if (cent < 0. || cent > 100.) {
      return;
    }
    int cBin = histManager.getCentBin(cent);
    spcHistograms.fill(HIST("FullCentrality"), cent);
    int nTracks = tracks.size();

    double wNUA = 1.0;
    double wEff = 1.0;
    for (const auto& track : tracks) {
      if (cfgFillQA) {
        // histManager.FillTrackQA<0>(track, cBin, collision.posZ());

        using JInputClassIter = typename TrackT::iterator;
        if constexpr (std::experimental::is_detected<HasWeightNUA, const JInputClassIter>::value)
          wNUA = track.weightNUA();
        if constexpr (std::experimental::is_detected<HasWeightEff, const JInputClassIter>::value)
          wEff = track.weightEff();
        histManager.fillTrackQA<1>(track, cBin, wEff, wNUA, collision.posZ());

        if constexpr (std::experimental::is_detected<HasWeightNUA, const JInputClassIter>::value) {
          spcAnalysis.fillQAHistograms(cBin, track.phi(), 1. / track.weightNUA());
        }
        if constexpr (std::experimental::is_detected<HasTrackType, const JInputClassIter>::value) {
          if (track.trackType() != cfgTrackBitMask.value) {
            LOGF(warning, "trackType %d (expected %d) is passed to the analysis", track.trackType(), cfgTrackBitMask.value);
          }
        }
      }
    }

    if (cfgFillQA)
      histManager.fillEventQA<1>(collision, cBin, cent, nTracks);

    jqvecs.Calculate(tracks, 0.0, cfgTrackCuts.cfgEtaMax);
    spcAnalysis.setQvectors(&jqvecs);
    spcAnalysis.calculateCorrelators(cBin);
  }

  template <class CollType>
  bool passOutlier(CollType const& collision)
  {
    if (cfgMultCutFormula.value.empty())
      return true;
    for (uint i = 0; i < aod::cfmultset::NMultiplicityEstimators; ++i) {
      if ((cfgMultCorrelationsMask.value & (1u << i)) == 0 || multCutFormulaParamIndex[i] == ~0u)
        continue;
      auto estIndex = std::popcount(static_cast<uint32_t>(cfgMultCorrelationsMask.value & ((1u << i) - 1)));
      multCutFormula->SetParameter(multCutFormulaParamIndex[i], collision.multiplicities()[estIndex]);
    }
    return multCutFormula->Eval() > 0.0f;
  }

  void processJDerived(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processJDerived, "Process derived data", false);

  void processJDerivedCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processJDerivedCorrected, "Process derived data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processCFDerived, "Process CF derived data", false);

  void processCFDerivedCorrected(soa::Filtered<aod::CFCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processCFDerivedCorrected, "Process CF derived data with corrections", true);

  void processCFDerivedMultSetCorrected(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>>::iterator const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    if (std::popcount(static_cast<uint32_t>(cfgMultCorrelationsMask.value)) != static_cast<int>(collision.multiplicities().size()))
      LOGF(fatal, "Multiplicity selections (cfgMultCorrelationsMask = 0x%x) do not match the size of the table column (%ld).", cfgMultCorrelationsMask.value, collision.multiplicities().size());
    if (!passOutlier(collision))
      return;
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(flowJSPCAnalysis, processCFDerivedMultSetCorrected, "Process CF derived data with corrections and multiplicity sets", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<flowJSPCAnalysis>(cfgc)};
}
