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
/// \author Jasper Parkkila (jparkkil@cern.ch)
/// \author Dong Jo Kim (djkim@jyu.fi)
/// \since Sep 2022

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "ReconstructionDataFormats/V0.h"

#include <TFormula.h>

#include <array>
#include <deque>
#include <memory>
#include <string>

// #include "CCDB/BasicCCDBManager.h"

#include "JFFlucAnalysis.h"
#include "JFFlucAnalysisO2Hist.h"

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/JCatalyst.h"

#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct jflucAnalysisTask {
  ~jflucAnalysisTask()
  {
    if (pcf)
      delete pcf;
    if (pcf2Prong)
      delete pcf2Prong;
  }

  O2_DEFINE_CONFIGURABLE(etamin, float, 0.4, "Minimum eta for tracks");
  O2_DEFINE_CONFIGURABLE(etamax, float, 0.8, "Maximum eta for tracks");
  O2_DEFINE_CONFIGURABLE(ptmin, float, 0.2, "Minimum pt for tracks");
  O2_DEFINE_CONFIGURABLE(ptmax, float, 5.0, "Maximum pt for tracks");
  O2_DEFINE_CONFIGURABLE(cfgCentBinsForMC, int, 0, "0 = OFF and 1 = ON for data like multiplicity/centrality bins for MC process");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrelationsMask, uint16_t, 0, "Selection bitmask for the multiplicity correlations. This should match the filter selection cfgEstimatorBitMask.")
  O2_DEFINE_CONFIGURABLE(cfgMultCutFormula, std::string, "", "Multiplicity correlations cut formula. A result greater than zero results in accepted event. Parameters: [cFT0C] FT0C centrality, [mFV0A] V0A multiplicity, [mGlob] global track multiplicity, [mPV] PV track multiplicity")

  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis phiAxis{"axisPhi", {50, 0.0, o2::constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis etaAxis{"axisEta", {40, -2.0, 2.0}, "eta axis for histograms"};
  ConfigurableAxis zvtAxis{"axisZvt", {20, -10.0, 10.0}, "zvertex axis for histograms"};
  ConfigurableAxis ptAxis{"axisPt", {60, 0.0, 300.0}, "pt axis for histograms"};
  ConfigurableAxis massAxis{"axisMass", {1, 0.0, 10.0}, "mass axis for histograms"};

  ConfigurableAxis vnCorrAxis{"vnCorrAxis", {2048, -0.1, 0.1}, "vn correlation axis"};
  ConfigurableAxis fourCorrSCAxis{"4pCorrAxisSC", {2048, -0.001, 0.001}, "4-particle correlation axis for SC"};
  ConfigurableAxis twoCorrSCAxis{"2pCorrAxisSC", {2048, -0.1, 0.1}, "2-particle correlation axis for SC"};
  ConfigurableAxis mixedCorrAxis{"mixedCorrAxis", {2048, -3.0, 3.0}, "N-particle correlation axis"};

  Filter jtrackFilter = (aod::jtrack::pt > ptmin) && (aod::jtrack::pt < ptmax);                                                     // eta cuts done by jfluc
  Filter cftrackFilter = (aod::cftrack::pt > ptmin) && (aod::cftrack::pt < ptmax);                                                  // eta cuts done by jfluc
  Filter cfmcparticleFilter = (aod::cfmcparticle::pt > ptmin) && (aod::cfmcparticle::pt < ptmax) && (aod::cfmcparticle::sign != 0); // eta cuts done by jfluc
  Filter cf2pFilter = (aod::cf2prongtrack::pt > ptmin) && (aod::cf2prongtrack::pt < ptmax);

  HistogramRegistry registry{"registry"};

  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint, 4> multCutFormulaParamIndex;

  void init(InitContext const&)
  {
    auto axisSpecMult = AxisSpec(axisMultiplicity);
    auto axisSpecPhi = AxisSpec(phiAxis);
    auto axisSpecEta = AxisSpec(etaAxis);
    auto axisSpecZvt = AxisSpec(zvtAxis);
    auto axisSpecPt = AxisSpec(ptAxis);
    auto axisSpecMass = AxisSpec(massAxis);
    auto axisSpecVn = AxisSpec(vnCorrAxis);
    auto axisSpec4pSC = AxisSpec(fourCorrSCAxis);
    auto axisSpec2pSC = AxisSpec(twoCorrSCAxis);
    auto axisSpecMixed = AxisSpec(mixedCorrAxis);
    if (doprocessJDerived || doprocessJDerivedCorrected || doprocessCFDerived || doprocessCFDerivedCorrected || doprocessCFDerivedMultSet || doprocessCFDerivedMultSetCorrected || doprocessMCCFDerived) {
      pcf = new JFFlucAnalysisO2Hist(registry, axisSpecMult, axisSpecPhi, axisSpecEta, axisSpecZvt, axisSpecPt, axisSpecMass, axisSpecVn, axisSpec4pSC, axisSpec2pSC, axisSpecMixed, cfgMultCorrelationsMask, "jfluc");
      pcf->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);
      pcf->UserCreateOutputObjects();
    } else {
      pcf = 0;
    }
    if (doprocessCF2ProngDerived || doprocessCF2ProngDerivedCorrected) {
      pcf2Prong = new JFFlucAnalysisO2Hist(registry, axisSpecMult, axisSpecPhi, axisSpecEta, axisSpecZvt, axisSpecPt, axisSpecMass, axisSpecVn, axisSpec4pSC, axisSpec2pSC, axisSpecMixed, cfgMultCorrelationsMask, "jfluc2prong");
      pcf2Prong->AddFlags(JFFlucAnalysis::kFlucEbEWeighting);
      pcf2Prong->UserCreateOutputObjects();

      ConfigurableAxis axisInvMassHistogram{"axisInvMassHistogram", {1000, 1.0, 3.0}, "invariant mass histogram binning"};
      registry.add("invMass", "2-prong invariant mass (GeV/c^2)", {HistType::kTH3F, {axisInvMassHistogram, {8, 0.0, 8.0, "p_{T}"}, axisMultiplicity}});
    } else {
      pcf2Prong = 0;
    }
    if ((doprocessCFDerivedMultSet || doprocessCFDerivedMultSetCorrected) && cfgMultCorrelationsMask == 0)
      LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when MultSet process functions are in use.");

    if (!cfgMultCutFormula.value.empty()) {
      if (cfgMultCorrelationsMask == 0)
        LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when outlier cuts are enabled.");
      multCutFormula = std::make_unique<TFormula>("multCutFormula", cfgMultCutFormula.value.c_str());
      std::fill_n(multCutFormulaParamIndex.begin(), std::size(multCutFormulaParamIndex), ~0u);
      std::array<std::string, 4> pars = {"cFT0C", "mFV0A", "mPV", "mGlob"}; // must correspond the order of MultiplicityEstimators
      for (uint i = 0, n = multCutFormula->GetNpar(); i < n; ++i) {
        auto m = std::find(pars.begin(), pars.end(), multCutFormula->GetParName(i));
        if (m == pars.end()) {

          LOGF(warning, "Unknown parameter in cfgMultCutFormula: %s", multCutFormula->GetParName(i));
          continue;
        }
        if ((cfgMultCorrelationsMask.value & (1u << i)) == 0) {
          LOGF(warning, "The centrality/multiplicity estimator %s is not available to be used in cfgMultCutFormula. Ensure cfgMultCorrelationsMask is correct and matches the CFMultSets in derived data.");
        } else {
          multCutFormulaParamIndex[std::distance(pars.begin(), m)] = i;
          LOGF(info, "Multiplicity cut parameter %s in use.", m->c_str());
        }
      }
    }
  }

  template <class T>
  using hasInvMass = decltype(std::declval<T&>().invMass());

  template <class CollisionT, class TrackT>
  void analyze(CollisionT const& collision, TrackT const& tracks)
  {
    pcf->Init();
    pcf->SetEventCentrality(collision.multiplicity());
    pcf->SetEventVertex(collision.posZ());
    pcf->FillQA(tracks);
    pcf->FillMultSet(collision);
    qvecs.Calculate(tracks, etamin, etamax);
    pcf->SetJQVectors(&qvecs);
    pcf->UserExec("");
  }

  template <class CollisionT, class POITrackT, class REFTrackT>
  void analyze(CollisionT const& collision, POITrackT const& poiTracks, REFTrackT const& refTracks)
  {
    if constexpr (std::experimental::is_detected<hasInvMass, typename POITrackT::iterator>::value) {
      for (auto& track : poiTracks)
        registry.fill(HIST("invMass"), track.invMass(), track.pt(), collision.multiplicity());
    }
    pcf2Prong->Init();
    pcf2Prong->SetEventCentrality(collision.multiplicity());
    pcf2Prong->SetEventVertex(collision.posZ());
    pcf2Prong->FillQA(refTracks, 0u);
    pcf2Prong->FillQA(poiTracks, 1u); // type = 1, all POI tracks in this list are of the same type
    qvecsRef.Calculate(refTracks, etamin, etamax);
    pcf2Prong->SetJQVectors(&qvecs, &qvecsRef);
    const AxisSpec& a = AxisSpec(massAxis);
    for (uint i = 0; i < a.getNbins(); ++i) {
      qvecs.Calculate(poiTracks, etamin, etamax, a.binEdges[i], a.binEdges[i + 1]);
      pcf2Prong->SetAverageInvariantMass(0.5f * (a.binEdges[i] + a.binEdges[i + 1]));
      pcf2Prong->UserExec(""); // The analysis needs to be called many times, once for each mass bin. For each of the bins, SetInvariantMass is used
    }
  }

  template <class CollType>
  bool passOutlier(CollType const& collision)
  {
    if (cfgMultCutFormula.value.empty())
      return true;
    for (uint i = 0; i < 4; ++i) {
      if ((cfgMultCorrelationsMask.value & (1u << i)) == 0 || multCutFormulaParamIndex[i] == ~0u)
        continue;
      auto estIndex = std::popcount(cfgMultCorrelationsMask.value & ((1u << i) - 1));
      multCutFormula->SetParameter(multCutFormulaParamIndex[i], collision.multiplicities()[estIndex]);
    }
    return multCutFormula->Eval() > 0.0f;
  }

  void processJDerived(aod::JCollision const& collision, soa::Filtered<aod::JTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerived, "Process derived data", false);

  void processJDerivedCorrected(aod::JCollision const& collision, soa::Filtered<soa::Join<aod::JTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processJDerivedCorrected, "Process derived data with corrections", false);

  void processCFDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerived, "Process CF derived data", false);

  void processCFDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedCorrected, "Process CF derived data with corrections", true);

  void processCFDerivedMultSet(soa::Join<aod::CFCollisions, aod::CFMultSets>::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    if (std::popcount(cfgMultCorrelationsMask.value) != static_cast<int>(collision.multiplicities().size()))
      LOGF(fatal, "Multiplicity selections (cfgMultCorrelationsMask = 0x%x) do not match the size of the table column (%ld). The histogram filling relies on the preservation of order.", cfgMultCorrelationsMask.value, collision.multiplicities().size());
    if (!passOutlier(collision))
      return;
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedMultSet, "Process CF derived data with multiplicity sets", false);

  void processCFDerivedMultSetCorrected(soa::Join<aod::CFCollisions, aod::CFMultSets>::iterator const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks)
  {
    if (std::popcount(cfgMultCorrelationsMask.value) != static_cast<int>(collision.multiplicities().size()))
      LOGF(fatal, "Multiplicity selections (cfgMultCorrelationsMask = 0x%x) do not match the size of the table column (%ld). The histogram filling relies on the preservation of order.", cfgMultCorrelationsMask.value, collision.multiplicities().size());
    if (!passOutlier(collision))
      return;
    analyze(collision, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCFDerivedMultSetCorrected, "Process CF derived data with corrections and multiplicity sets", false);

  void processCF2ProngDerived(aod::CFCollision const& collision, soa::Filtered<aod::CFTracks> const& tracks, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    analyze(collision, p2tracks, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCF2ProngDerived, "Process CF derived data with 2-prongs as POI and charged particles as REF.", false);

  void processCF2ProngDerivedCorrected(aod::CFCollision const& collision, soa::Filtered<soa::Join<aod::CFTracks, aod::JWeights>> const& tracks, soa::Filtered<soa::Join<aod::CF2ProngTracks, aod::J2ProngWeights>> const& p2tracks)
  {
    analyze(collision, p2tracks, tracks);
  }
  PROCESS_SWITCH(jflucAnalysisTask, processCF2ProngDerivedCorrected, "Process CF derived data with 2-prongs as POI and charged particles as REF with corrections.", false);

  void processMCCFDerived(aod::CFMcCollision const& mcCollision, soa::Filtered<aod::CFMcParticles> const& particles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    auto multiplicity = mcCollision.multiplicity();
    if (cfgCentBinsForMC > 0) {
      if (collisions.size() == 0) {
        return;
      }
      for (const auto& collision : collisions) {
        multiplicity = collision.multiplicity();
      }
    }
    pcf->Init();
    pcf->SetEventCentrality(multiplicity);
    pcf->SetEventVertex(mcCollision.posZ());
    pcf->FillQA(particles);
    qvecs.Calculate(particles, etamin, etamax);
    pcf->SetJQVectors(&qvecs);
    pcf->UserExec("");
  }
  PROCESS_SWITCH(jflucAnalysisTask, processMCCFDerived, "Process CF derived MC data", false);

  JFFlucAnalysis::JQVectorsT qvecs;
  JFFlucAnalysis::JQVectorsT qvecsRef;
  JFFlucAnalysisO2Hist* pcf;
  JFFlucAnalysisO2Hist* pcf2Prong;
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    // adaptAnalysisTask<jflucTypeLoader>(cfgc),
    adaptAnalysisTask<jflucAnalysisTask>(cfgc)};
}
