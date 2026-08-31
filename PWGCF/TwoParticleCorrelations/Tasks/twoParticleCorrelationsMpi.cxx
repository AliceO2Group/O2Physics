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

/// \file twoParticleCorrelationsMpi.cxx
/// \brief task for the MPI-proxy classification of correlation calculations with CF-filtered tracks for O2 analysis
/// \author Emil Gorm Dahlbæk Nielsen <emil.gorm.nielsen@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TFormula.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>
#include <TTree.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <experimental/type_traits>
#include <functional>
#include <iterator>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

static constexpr std::array<std::array<float, 5>, 1> CfgPairCutDefaults = {{{-1, -1, -1, -1, -1}}};

struct TwoParticleCorrelationsMpi {
  SliceCache cache;

  // Configuration
  Configurable<float> cfgCutVertex{"cfgCutVertex", 7.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.5f, "Minimal pT for tracks"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};

  Configurable<int> cfgPtOrder{"cfgPtOrder", 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)"};
  Configurable<int> cfgTriggerCharge{"cfgTriggerCharge", 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative"};
  Configurable<int> cfgAssociatedCharge{"cfgAssociatedCharge", 0, "Select on charge of associated particle: 0 = all charged; 1 = positive; -1 = negative"};
  Configurable<int> cfgPairCharge{"cfgPairCharge", 0, "Select on charge of particle pair: 0 = all; 1 = like sign; -1 = unlike sign"};
  Configurable<int> cfgCorrelationMethod{"cfgCorrelationMethod", 0, "Correlation method, 0 = all, 1 = dd, 2 = ddbar"};

  Configurable<float> cfgTwoTrackCut{"cfgTwoTrackCut", -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)"};
  Configurable<float> cfgTwoTrackCutMinRadius{"cfgTwoTrackCutMinRadius", 0.8f, "Two track cut: radius in m from which two track cuts are applied"};
  ;
  Configurable<int> cfgLocalEfficiency{"cfgLocalEfficiency", 0, "0 = OFF and 1 = ON for local efficiency"};
  Configurable<bool> cfgDropStepRECO{"cfgDropStepRECO", false, "choice to drop step RECO if efficiency correction is used"};
  Configurable<int> cfgCentBinsForMC{"cfgCentBinsForMC", 0, "0 = OFF and 1 = ON for data like multiplicity/centrality bins for MC steps"};
  Configurable<uint16_t> cfgTrackBitMask{"cfgTrackBitMask", 0, "BitMask for track selection systematics; refer to the enum TrackSelectionCuts in filtering task"};
  Configurable<uint16_t> cfgMultCorrelationsMask{"cfgMultCorrelationsMask", 0, "Selection bitmask for the multiplicity correlations. This should match the filter selection cfgEstimatorBitMask."};
  Configurable<std::string> cfgMultCutFormula{"cfgMultCutFormula", "", "Multiplicity correlations cut formula. A result greater than zero results in accepted event. Parameters: [cFT0C] FT0C centrality, [mFV0A] V0A multiplicity, [mGlob] global track multiplicity, [mPV] PV track multiplicity, [cFT0M] FT0M centrality"};

  // Suggested values: Photon: 0.004; K0 and Lambda: 0.005
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {CfgPairCutDefaults.front().data(), 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  Configurable<std::string> cfgEfficiencyTrigger{"cfgEfficiencyTrigger", "", "CCDB path to efficiency object for trigger particles"};
  Configurable<std::string> cfgEfficiencyAssociated{"cfgEfficiencyAssociated", "", "CCDB path to efficiency object for associated particles"};
  Configurable<std::string> cfgNuncSeedsTemplateFile{"cfgNuncSeedsTemplateFile", "", "Local ROOT file containing ensembleYieldTemplates"};
  Configurable<std::string> cfgNuncSeedsTemplate{"cfgNuncSeedsTemplate", "", "CCDB path to the ensemble-yield template ccdb_object"};
  Configurable<float> cfgMinPairAcceptance{"cfgMinPairAcceptance", 0.05f, "Minimum pair acceptance used by the event seed estimator"};
  Configurable<float> cfgMaxPairAcceptanceWeight{"cfgMaxPairAcceptanceWeight", 20.f, "Maximum allowed inverse pair-acceptance weight"};
  Configurable<int> cfgEventSeedEstimatorMethod{"cfgEventSeedEstimatorMethod", 0, "Event seed estimator: 0 = fixed ensemble pair probabilities, 1 = Gamma-regularized weighted MAP-EM"};
  Configurable<float> cfgEventSeedPriorExposure{"cfgEventSeedPriorExposure", 1.f, "Gamma-prior strength in equivalent event exposures for MAP-EM"};
  Configurable<int> cfgEventSeedEMMaxIterations{"cfgEventSeedEMMaxIterations", 10, "Maximum number of weighted MAP-EM iterations per event"};
  Configurable<float> cfgEventSeedEMTolerance{"cfgEventSeedEMTolerance", 1.e-5f, "Relative component-yield convergence tolerance for MAP-EM"};
  Configurable<std::vector<float>> cfgEventSeedPercentileLower{"cfgEventSeedPercentileLower", {}, "Lower event-seed percentile boundary in each axisMultiplicity bin; empty keeps the raw Seed user axis"};
  Configurable<std::vector<float>> cfgEventSeedPercentileUpper{"cfgEventSeedPercentileUpper", {}, "Upper event-seed percentile boundary in each axisMultiplicity bin; empty keeps the raw Seed user axis"};

  Configurable<int> cfgNumMixedEvents{"cfgNumMixedEvents", 5, "Number of mixed events per event"};

  Configurable<int> cfgVerbosity{"cfgVerbosity", 1, "Verbosity level (0 = major, 1 = per collision)"};

  Configurable<int> cfgDecayParticleMask{"cfgDecayParticleMask", 0, "Selection bitmask for the decay particles: 0 = no selection"};
  Configurable<float> cfgV0RapidityMax{"cfgV0RapidityMax", 0.8, "Maximum rapidity for the decay particles (0 = no selection)"};
  Configurable<int> cfgUserAxis{"cfgUserAxis", 0, "Additional user axis: 0 = OFF, 1 = invariant mass, 2 = event seed"};
  Configurable<std::vector<int>> cfgMcTriggerPDGs{"cfgMcTriggerPDGs", {}, "MC PDG codes to use exclusively as trigger particles and exclude from associated particles. Empty = no selection."};

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  ConfigurableAxis axisUser{"axisUser", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 5.0}, "additional user axis for histograms"};

  ConfigurableAxis axisMultCorrCent{"axisMultCorrCent", {100, 0, 100}, "multiplicity correlation axis for centralities"};
  ConfigurableAxis axisMultCorrV0{"axisMultCorrV0", {1000, 0, 100000}, "multiplicity correlation axis for V0 multiplicities"};
  ConfigurableAxis axisMultCorrMult{"axisMultCorrMult", {1000, 0, 1000}, "multiplicity correlation axis for track multiplicities"};

  // This filter is applied to AOD and derived data (column names are identical)
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // This filter is only applied to AOD
  Filter collisionVertexTypeFilter = (aod::collision::flags & static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks)) == static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks);

  // Track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t)true));
  Filter cfTrackFilter = (nabs(aod::cftrack::eta) < cfgCutEta) && (aod::cftrack::pt > cfgCutPt) && ncheckbit(aod::track::trackType, as<uint8_t>(cfgTrackBitMask));

  // MC filters
  Filter cfMCCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  Filter cfMCParticleFilter = (nabs(aod::cfmcparticle::eta) < cfgCutEta) && (aod::cfmcparticle::pt > cfgCutPt); // && (aod::cfmcparticle::sign != 0); //check the sign manually, some specials may be neutral
  Filter mcParticleFilter = (nabs(aod::mcparticle::eta) < cfgCutEta) && (aod::mcparticle::pt > cfgCutPt);

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  // persistent caches
  std::vector<float> efficiencyAssociatedCache;

  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint, aod::cfmultset::NMultiplicityEstimators> multCutFormulaParamIndex{};

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  struct YieldTemplate {
    int trigBin = -1;
    int assocBin = -1;
    int multBin = -1;
    double nchLow = 0.0;
    double nchHigh = 0.0;
    double trigPtLow = 0.0;
    double trigPtHigh = 0.0;
    double assocPtLow = 0.0;
    double assocPtHigh = 0.0;
    std::array<double, 10> parameters{};
  };

  struct EventSeedPairObservation {
    std::size_t templateIndex = 0;
    std::array<double, 3> normalizedShapes{}; // near, away, baseline densities normalized by their fitted yields
    double acceptanceWeight = 0.0;
  };

  struct EventSeedEstimate {
    int nTriggers = 0;
    int nCandidatePairs = 0;
    int nPairs = 0;
    int nPairsWithoutTemplate = 0;
    double nearPairs = 0.0;
    double awayPairs = 0.0;
    double baselinePairs = 0.0;
    double probabilityNearPairs = 0.0;
    double probabilityAwayPairs = 0.0;
    double probabilityBaselinePairs = 0.0;
    double rawNearPairs = 0.0;
    double rawAwayPairs = 0.0;
    double rawBaselinePairs = 0.0;
    int nAcceptanceCorrectedPairs = 0;
    int nPairsRejectedByAcceptance = 0;
    double sumAcceptanceWeights = 0.0;
    std::vector<EventSeedPairObservation> pairObservations;
    std::vector<std::array<double, 3>> templatePriorCounts;
    std::array<double, 3> priorComponentCounts{};
    int emIterations = 0;
    bool emConverged = false;
    bool usedMapEM = false;

    [[nodiscard]] bool hasTriggers() const { return nTriggers > 0; }
    [[nodiscard]] bool hasAcceptanceCorrectedPairs() const { return nAcceptanceCorrectedPairs > 0; }
    [[nodiscard]] bool isValid() const { return hasTriggers() && hasAcceptanceCorrectedPairs(); }
    [[nodiscard]] double nearYield() const { return hasTriggers() ? nearPairs / nTriggers : 0.0; }
    [[nodiscard]] double awayYield() const { return hasTriggers() ? awayPairs / nTriggers : 0.0; }
    [[nodiscard]] double probabilityNearYield() const { return hasTriggers() ? probabilityNearPairs / nTriggers : 0.0; }
    [[nodiscard]] double probabilityAwayYield() const { return hasTriggers() ? probabilityAwayPairs / nTriggers : 0.0; }
    [[nodiscard]] double probabilityNuncSeeds() const
    {
      const double denominator = 1.0 + probabilityNearYield() + probabilityAwayYield();
      return isValid() && denominator > 0.0 ? nTriggers / denominator : -1.0;
    }
    [[nodiscard]] double nuncSeeds() const
    {
      const double denominator = 1.0 + nearYield() + awayYield();
      return isValid() && denominator > 0.0 ? nTriggers / denominator : -1.0;
    }
  };

  struct PendingSeedTriggerFill {
    float pt;
    float multiplicity;
    float posZ;
    float weight;
  };

  struct PendingSeedPairFill {
    float deltaEta;
    float assocPt;
    float triggerPt;
    float multiplicity;
    float deltaPhi;
    float posZ;
    float weight;
  };

  std::vector<YieldTemplate> yieldTemplates;
  std::vector<std::unique_ptr<TH3D>> pairAcceptanceMaps;
  std::vector<std::unique_ptr<TH2D>> pairAcceptanceEtaVertexMaps;
  int pairAcceptanceSchemaVersion = 0;
  const TList* loadedCcdbYieldTemplateObject = nullptr;
  bool eventSeedEstimatorEnabled = false;
  bool eventSeedPercentileAxisEnabled = false;

  enum CorrelationMethod {
    All = 0,
    Dd,
    Ddbar
  };
  enum UserAxisMode {
    NoUserAxis = 0,
    InvariantMassAxis,
    EventSeedAxis
  };
  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  Service<o2::framework::O2DatabasePDG> pdg{};

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  using DerivedCollisions = soa::Filtered<aod::CFCollisions>;
  using DerivedTracks = soa::Filtered<aod::CFTracks>;

  void init(o2::framework::InitContext&)
  {
    if (cfgUserAxis < NoUserAxis || cfgUserAxis > EventSeedAxis) {
      LOGF(fatal, "Unsupported cfgUserAxis=%d; use 0 (off), 1 (invariant mass), or 2 (event seed)", cfgUserAxis.value);
    }
    if (doprocessMCSameDerived && (doprocessSameDerived || doprocessSameDerivedMultSet)) {
      LOGF(fatal, "processMCSameDerived is mutually exclusive with the reconstructed derived same-event processes because it also fills those outputs");
    }
    if (doprocessSameGenMC && doprocessMCSameDerived) {
      LOGF(fatal, "processSameGenMC and processMCSameDerived are mutually exclusive because both fill the generated same-event outputs");
    }
    if (!cfgNuncSeedsTemplateFile.value.empty() && !cfgNuncSeedsTemplate.value.empty()) {
      LOGF(fatal, "Configure only one template source: cfgNuncSeedsTemplateFile or cfgNuncSeedsTemplate");
    }
    if (cfgEventSeedEstimatorMethod < 0 || cfgEventSeedEstimatorMethod > 1) {
      LOGF(fatal, "Unsupported cfgEventSeedEstimatorMethod=%d; use 0 (fixed probabilities) or 1 (weighted MAP-EM)", cfgEventSeedEstimatorMethod.value);
    }
    if (cfgEventSeedPriorExposure < 0.f || cfgEventSeedEMMaxIterations < 1 || cfgEventSeedEMTolerance <= 0.f) {
      LOGF(fatal, "MAP-EM configuration requires non-negative prior exposure, at least one iteration, and positive tolerance");
    }
    eventSeedEstimatorEnabled = !cfgNuncSeedsTemplateFile.value.empty() || !cfgNuncSeedsTemplate.value.empty();
    if (cfgUserAxis == EventSeedAxis && !eventSeedEstimatorEnabled) {
      LOGF(fatal, "cfgUserAxis=2 requires an event-seed template configured through cfgNuncSeedsTemplateFile or cfgNuncSeedsTemplate");
    }
    const bool hasLowerPercentileBoundaries = !cfgEventSeedPercentileLower->empty();
    const bool hasUpperPercentileBoundaries = !cfgEventSeedPercentileUpper->empty();
    if (hasLowerPercentileBoundaries != hasUpperPercentileBoundaries) {
      LOGF(fatal, "Configure both cfgEventSeedPercentileLower and cfgEventSeedPercentileUpper, or leave both empty");
    }
    eventSeedPercentileAxisEnabled = hasLowerPercentileBoundaries;
    if (eventSeedPercentileAxisEnabled) {
      const int nMultiplicityBins = AxisSpec(axisMultiplicity).getNbins();
      if (static_cast<int>(cfgEventSeedPercentileLower->size()) != nMultiplicityBins ||
          static_cast<int>(cfgEventSeedPercentileUpper->size()) != nMultiplicityBins) {
        LOGF(fatal, "Event-seed percentile boundary vectors must each contain exactly %d values, one per axisMultiplicity bin", nMultiplicityBins);
      }
      for (int multBin = 0; multBin < nMultiplicityBins; ++multBin) {
        if (!std::isfinite(cfgEventSeedPercentileLower->at(multBin)) ||
            !std::isfinite(cfgEventSeedPercentileUpper->at(multBin)) ||
            cfgEventSeedPercentileLower->at(multBin) >= cfgEventSeedPercentileUpper->at(multBin)) {
          LOGF(fatal, "Invalid event-seed percentile boundaries in multiplicity bin %d: lower=%g upper=%g",
               multBin, cfgEventSeedPercentileLower->at(multBin), cfgEventSeedPercentileUpper->at(multBin));
        }
      }
      if (cfgUserAxis == EventSeedAxis) {
        const std::vector<double> expectedEdges{-1.5, -0.5, 0.5, 1.5, 2.5};
        const auto& configuredEdges = AxisSpec(axisUser).binEdges;
        const float threshold = 1.e-6;
        if (configuredEdges.size() != expectedEdges.size() ||
            !std::equal(configuredEdges.begin(), configuredEdges.end(), expectedEdges.begin(),
                        [&threshold](double lhs, double rhs) { return std::abs(lhs - rhs) < threshold; })) {
          LOGF(fatal, "Percentile-class Seed axis requires axisUser edges {-1.5,-0.5,0.5,1.5,2.5}");
        }
      }
    }
    LOGF(info, "Event seed estimator histogram booking: %s (local template='%s', CCDB template='%s')",
         eventSeedEstimatorEnabled ? "enabled" : "disabled",
         cfgNuncSeedsTemplateFile.value.c_str(), cfgNuncSeedsTemplate.value.c_str());

    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, o2::constants::math::TwoPI, "#varphi"}}});
    if (doprocessSameDerivedMultSet) {
      if (cfgMultCorrelationsMask == 0) {
        LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when MultSet process functions are in use.");
      }
      std::vector<AxisSpec> multAxes;
      if ((cfgMultCorrelationsMask & aod::cfmultset::CentFT0C) != 0) {
        multAxes.emplace_back(axisMultCorrCent, "FT0C centrality");
      }
      if ((cfgMultCorrelationsMask & aod::cfmultset::MultFV0A) != 0) {
        multAxes.emplace_back(axisMultCorrV0, "V0A multiplicity");
      }
      if ((cfgMultCorrelationsMask & aod::cfmultset::MultNTracksPV) != 0) {
        multAxes.emplace_back(axisMultCorrMult, "Nch PV");
      }
      if ((cfgMultCorrelationsMask & aod::cfmultset::MultNTracksGlobal) != 0) {
        multAxes.emplace_back(axisMultCorrMult, "Nch Global");
      }
      if ((cfgMultCorrelationsMask & aod::cfmultset::CentFT0M) != 0) {
        multAxes.emplace_back(axisMultCorrCent, "FT0M centrality");
      }
      registry.add("multCorrelations", "Multiplicity correlations", {HistType::kTHnSparseF, multAxes});
    }
    registry.add("multiplicity", "event multiplicity", {HistType::kTH1F, {{1000, 0, 100, "/multiplicity/centrality"}}});
    if (eventSeedEstimatorEnabled) {
      registry.add("eventSeedEstimator", "event-level template estimator", {HistType::kTHnSparseF, {{100, 0, 100, "multiplicity"}, {100, -0.5, 99.5, "N_{trig}"}, {200, 0, 20, "Y_{near}"}, {200, 0, 20, "Y_{away}"}, {200, 0, 100, "N_{uncorrelated seeds}"}}});
      registry.add("eventSeedPairProbabilities", "summed pair probabilities", {HistType::kTH3F, {{200, 0, 200, "#Sigma P_{baseline}"}, {200, 0, 200, "#Sigma P_{near}"}, {200, 0, 200, "#Sigma P_{away}"}}});
      registry.add("eventSeedEstimatorStatus", "event-estimator status", {HistType::kTH1F, {{5, -0.5, 4.5, "status"}}});
      registry.add("eventSeedTemplateCoverage", "template coverage vs multiplicity", {HistType::kTH2F, {{100, 0, 100, "multiplicity"}, {102, -0.01, 1.01, "matched/candidate pairs"}}});
      registry.add("eventSeedEstimateVsMultiplicity", "event-level seed estimate vs multiplicity", {HistType::kTH2F, {{100, 0, 100, "multiplicity"}, {200, 0, 100, "N_{uncorrelated seeds}"}}});
      registry.add("eventSeedFixedProbabilityVsMultiplicity", "fixed-probability event Seed versus multiplicity;multiplicity;N_{seed}^{probability sum}", {HistType::kTH2F, {axisMultiplicity, {1000, 0, 100}}});
      registry.add("eventNearYieldVsMultiplicity", "event-level near yield vs multiplicity", {HistType::kTH2F, {{100, 0, 100, "multiplicity"}, {200, 0, 20, "Y_{near}"}}});
      registry.add("eventAwayYieldVsMultiplicity", "event-level away yield vs multiplicity", {HistType::kTH2F, {{100, 0, 100, "multiplicity"}, {200, 0, 20, "Y_{away}"}}});
      registry.add("profileEventNTriggers", "mean event trigger count", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventNearPairs", "mean summed near-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventAwayPairs", "mean summed away-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventBaselinePairs", "mean summed baseline-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventRawNearPairs", "mean raw summed near-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventRawAwayPairs", "mean raw summed away-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventRawBaselinePairs", "mean raw summed baseline-pair probability", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventAcceptanceCoverage", "fraction of template-matched pairs with valid acceptance", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventMeanAcceptanceWeight", "mean inverse acceptance weight for corrected pairs", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("eventSeedAcceptanceWeight", "inverse pair-acceptance weight", {HistType::kTH1F, {{200, 0, 20, "1/A"}}});
      registry.add("profileEventNearYield", "mean event near yield", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventAwayYield", "mean event away yield", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventNuncSeeds", "mean event uncorrelated-seed estimate", {HistType::kTProfile, {axisMultiplicity}});
      registry.add("profileEventEstimatorValidity", "fraction of events with a valid seed estimate", {HistType::kTProfile, {axisMultiplicity}});
      if (cfgEventSeedEstimatorMethod == 1) {
        registry.add("eventSeedEMIterations", "weighted MAP-EM iterations;iterations;events", {HistType::kTH1F, {{cfgEventSeedEMMaxIterations + 1, -0.5, cfgEventSeedEMMaxIterations + 0.5}}});
        registry.add("profileEventEMConvergence", "weighted MAP-EM convergence fraction", {HistType::kTProfile, {axisMultiplicity}});
        registry.add("profileEventEMNearPrior", "mean MAP-EM near-component prior mode", {HistType::kTProfile, {axisMultiplicity}});
        registry.add("profileEventEMAwayPrior", "mean MAP-EM away-component prior mode", {HistType::kTProfile, {axisMultiplicity}});
        registry.add("profileEventEMBaselinePrior", "mean MAP-EM baseline-component prior mode", {HistType::kTProfile, {axisMultiplicity}});
        registry.add("eventSeedMAPVsProbability", "MAP-EM versus fixed-probability seed estimate;N_{seed}^{probability sum};N_{seed}^{MAP-EM}", {HistType::kTH2F, {{200, 0, 100}, {200, 0, 100}}});
        registry.add("eventSeedMAPEMVsMultiplicity", "MAP-EM event Seed versus multiplicity;multiplicity;N_{seed}^{MAP-EM}", {HistType::kTH2F, {axisMultiplicity, {1000, 0, 100}}});
        registry.add("profileEventProbabilityNuncSeeds", "mean fixed-probability seed estimate", {HistType::kTProfile, {axisMultiplicity}});
      }
      auto* estimatorStatus = registry.get<TH1>(HIST("eventSeedEstimatorStatus")).get();
      estimatorStatus->GetXaxis()->SetBinLabel(1, "unused");
      estimatorStatus->GetXaxis()->SetBinLabel(2, "no template-covered triggers");
      estimatorStatus->GetXaxis()->SetBinLabel(3, "no template-covered pairs");
      estimatorStatus->GetXaxis()->SetBinLabel(4, "no acceptance-corrected pairs");
      estimatorStatus->GetXaxis()->SetBinLabel(5, "valid estimate");
      registry.add("mcValidation/estimatedSeedsVsTrueNMPI", "template estimator response;N_{MPI}^{true};N_{seed}^{estimated}", {HistType::kTH2F, {{101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      registry.add("mcValidation/estimatedSeedsVsTrueNMPIVsMultiplicity", "event-level template estimator response;multiplicity;N_{MPI}^{true};N_{seed}^{estimated}", {HistType::kTH3F, {axisMultiplicity, {101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      if (cfgEventSeedEstimatorMethod == 1) {
        registry.add("mcValidation/probabilityEstimatedSeedsVsTrueNMPI", "fixed-probability estimator response;N_{MPI}^{true};N_{seed}^{probability sum}", {HistType::kTH2F, {{101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      }
      registry.add("mcValidation/profileEstimatedSeedsVsTrueNMPI", "mean template estimate;N_{MPI}^{true};#LT N_{seed}^{estimated} #GT", {HistType::kTProfile, {{101, -0.5, 100.5}}});
      registry.add("mcValidation/profileBiasVsTrueNMPI", "mean estimator bias;N_{MPI}^{true};#LT N_{seed}^{estimated} - N_{MPI}^{true} #GT", {HistType::kTProfile, {{101, -0.5, 100.5}}});
      registry.add("mcValidation/relativeResidualVsTrueNMPI", "relative estimator residual;N_{MPI}^{true};(N_{seed}^{estimated} - N_{MPI}^{true}) / N_{MPI}^{true}", {HistType::kTH2F, {{101, -0.5, 100.5}, {240, -3., 3.}}});
      registry.add("mcValidation/trueNMPIVsMultiplicity", "true MPI count versus reconstructed multiplicity;multiplicity;N_{MPI}^{true}", {HistType::kTH2F, {{100, 0., 100.}, {101, -0.5, 100.5}}});
      registry.add("mcValidation/templateCoverageVsTrueNMPI", "template pair coverage versus true MPI count;N_{MPI}^{true};matched / candidate pairs", {HistType::kTH2F, {{101, -0.5, 100.5}, {102, -0.01, 1.01}}});
      registry.add("mcValidation/status", "MC template-estimator validation status;status;events", {HistType::kTH1F, {{7, -0.5, 6.5}}});
      auto* mcValidationStatus = registry.get<TH1>(HIST("mcValidation/status")).get();
      mcValidationStatus->GetXaxis()->SetBinLabel(1, "no selected reconstructed collision");
      mcValidationStatus->GetXaxis()->SetBinLabel(2, "invalid N MPI");
      mcValidationStatus->GetXaxis()->SetBinLabel(3, "outside template multiplicity");
      mcValidationStatus->GetXaxis()->SetBinLabel(4, "no template-covered triggers");
      mcValidationStatus->GetXaxis()->SetBinLabel(5, "no template-covered pairs");
      mcValidationStatus->GetXaxis()->SetBinLabel(6, "no acceptance-corrected pairs");
      mcValidationStatus->GetXaxis()->SetBinLabel(7, "valid response");
      registry.add("mcValidation/generated/estimatedSeedsVsTrueNMPI", "generated-level template estimator response;N_{MPI}^{true};N_{seed,gen}^{estimated}", {HistType::kTH2F, {{101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      registry.add("mcValidation/generated/estimatedSeedsVsTrueNMPIVsMultiplicity", "generated event-level template estimator response;N_{ch}^{gen};N_{MPI}^{true};N_{seed,gen}^{estimated}", {HistType::kTH3F, {axisMultiplicity, {101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      if (cfgEventSeedEstimatorMethod == 1) {
        registry.add("mcValidation/generated/probabilityEstimatedSeedsVsTrueNMPI", "generated-level fixed-probability estimator response;N_{MPI}^{true};N_{seed,gen}^{probability sum}", {HistType::kTH2F, {{101, -0.5, 100.5}, {202, -0.5, 100.5}}});
      }
      registry.add("mcValidation/generated/profileEstimatedSeedsVsTrueNMPI", "mean generated-level template estimate;N_{MPI}^{true};#LT N_{seed,gen}^{estimated} #GT", {HistType::kTProfile, {{101, -0.5, 100.5}}});
      registry.add("mcValidation/generated/profileBiasVsTrueNMPI", "mean generated-level estimator bias;N_{MPI}^{true};#LT N_{seed,gen}^{estimated} - N_{MPI}^{true} #GT", {HistType::kTProfile, {{101, -0.5, 100.5}}});
      registry.add("mcValidation/generated/relativeResidualVsTrueNMPI", "generated-level relative estimator residual;N_{MPI}^{true};(N_{seed,gen}^{estimated} - N_{MPI}^{true}) / N_{MPI}^{true}", {HistType::kTH2F, {{101, -0.5, 100.5}, {240, -3., 3.}}});
      registry.add("mcValidation/generated/trueNMPIVsMultiplicity", "true MPI count versus generated charged multiplicity;N_{ch}^{gen};N_{MPI}^{true}", {HistType::kTH2F, {{101, -0.5, 100.5}, {101, -0.5, 100.5}}});
      registry.add("mcValidation/generated/templateCoverageVsTrueNMPI", "generated-level template pair coverage versus true MPI count;N_{MPI}^{true};matched / candidate pairs", {HistType::kTH2F, {{101, -0.5, 100.5}, {102, -0.01, 1.01}}});
      registry.add("mcValidation/generated/status", "generated-level MC template-estimator validation status;status;events", {HistType::kTH1F, {{6, -0.5, 5.5}}});
      auto* generatedValidationStatus = registry.get<TH1>(HIST("mcValidation/generated/status")).get();
      generatedValidationStatus->GetXaxis()->SetBinLabel(1, "invalid N MPI");
      generatedValidationStatus->GetXaxis()->SetBinLabel(2, "outside template multiplicity");
      generatedValidationStatus->GetXaxis()->SetBinLabel(3, "no template-covered triggers");
      generatedValidationStatus->GetXaxis()->SetBinLabel(4, "no template-covered pairs");
      generatedValidationStatus->GetXaxis()->SetBinLabel(5, "no acceptance-corrected pairs");
      generatedValidationStatus->GetXaxis()->SetBinLabel(6, "valid response");
    }
    registry.add("yvspt", "y vs pT", {HistType::kTH2F, {{100, -1, 1, "y"}, {100, 0, 20, "p_{T}"}}}); // y vs pT for all tracks (control histogram)

    const int maxMixBin = AxisSpec(axisMultiplicity).getNbins() * AxisSpec(axisVertex).getNbins();
    // The bin numbers for the control histograms (eventcount_*) come from getBin(...) and are the following: #mult_bin * #number_of_z_bins + #zbin
    registry.add("eventcount_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_same", "bin", {HistType::kTH2F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}}});
    registry.add("trackcount_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}, {10, -0.5, 9.5}}});

    mPairCuts.SetHistogramRegistry(&registry);

    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 || cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      cfg.mPairCuts = true;
    }

    if (cfgTwoTrackCut > 0) {
      mPairCuts.SetTwoTrackCuts(cfgTwoTrackCut, cfgTwoTrackCutMinRadius);
    }

    // --- OBJECT INIT ---

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
        if ((cfgMultCorrelationsMask.value & (1u << i)) == 0) {
          LOGF(warning, "The centrality/multiplicity estimator %s is not available to be used in cfgMultCutFormula. Ensure cfgMultCorrelationsMask is correct and matches the CFMultSets in derived data.");
        } else {
          multCutFormulaParamIndex[std::distance(pars.begin(), m)] = i;
          LOGF(info, "Multiplicity cut parameter %s in use.", m->c_str());
        }
      }
    }

    std::vector<AxisSpec> corrAxis;
    corrAxis.reserve(6);
    corrAxis.emplace_back(axisDeltaEta, "#Delta#eta");
    corrAxis.emplace_back(axisPtAssoc, "p_{T} (GeV/c)");
    corrAxis.emplace_back(axisPtTrigger, "p_{T} (GeV/c)");
    corrAxis.emplace_back(axisMultiplicity, "multiplicity / centrality");
    corrAxis.emplace_back(axisDeltaPhi, "#Delta#varphi (rad)");
    corrAxis.emplace_back(axisVertex, "z-vtx (cm)");
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis;
    std::vector<AxisSpec> userMixingAxis;

    if (cfgUserAxis == InvariantMassAxis) {
      userAxis.emplace_back(axisUser, "m (GeV/c^2)");
      userMixingAxis.emplace_back(axisUser, "m (GeV/c^2)");
    } else if (cfgUserAxis == EventSeedAxis) {
      const char* title = eventSeedPercentileAxisEnabled ? "N_{seed} percentile class" : "N_{seed}";
      userAxis.emplace_back(axisUser, title);
      userMixingAxis.emplace_back(axisUser, title);
    }
    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userMixingAxis));

    same->setTrackEtaCut(cfgCutEta);
    mixed->setTrackEtaCut(cfgCutEta);

    if (!cfgEfficiencyAssociated.value.empty()) {
      efficiencyAssociatedCache.reserve(512);
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time

    loadLocalYieldTemplates();
  }

  int getMagneticField(uint64_t timestamp)
  {
    // TODO done only once (and not per run). Will be replaced by CCDBConfigurable
    static o2::parameters::GRPObject* grpo = nullptr;
    // static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      // grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <class T>
  using HasMultSet = decltype(std::declval<T&>().multiplicities());

  template <typename TCollision, typename TTracks>
  void fillQA(const TCollision& collision, float multiplicity, const TTracks& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    if constexpr (std::experimental::is_detected<HasMultSet, TCollision>::value) {
      if (std::popcount(cfgMultCorrelationsMask.value) != static_cast<int>(collision.multiplicities().size())) {
        LOGF(fatal, "Multiplicity selections (cfgMultCorrelationsMask = 0x%x) do not match the size of the table column (%ld). The histogram filling relies on the preservation of order.", cfgMultCorrelationsMask.value, collision.multiplicities().size());
      }
      // need to convert to vec of doubles since THnSparse has no way to fill vec of floats directly
      std::vector<double> v(collision.multiplicities().begin(), collision.multiplicities().end());
      registry.get<THnSparse>(HIST("multCorrelations")).get()->Fill(v.data());
    }
    for (const auto& track1 : tracks) {
      registry.fill(HIST("yields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), multiplicity, track1.eta(), track1.phi());
    }
  }

  template <class T>
  using HasInvMass = decltype(std::declval<T&>().invMass());
  template <class T>
  using HasPDGCode = decltype(std::declval<T&>().pdgCode());

  template <typename TCollision, typename TTracks1, typename TTracks2>
  void fillQA(const TCollision& collision, float multiplicity, float posZ, const TTracks1& tracks1, const TTracks2& tracks2)
  {
    for (const auto& track1 : tracks1) {
      if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value && std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
        if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track1.decay()))) == 0u) {
          continue;
        }
        registry.fill(HIST("invMass"), track1.invMass(), track1.pt(), multiplicity, posZ);
      }
      if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
        if (!cfgMcTriggerPDGs->empty() && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track1.pdgCode()) == cfgMcTriggerPDGs->end()) {
          continue;
        }
      }
      registry.fill(HIST("yieldsTrigger"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphiTrigger"), multiplicity, track1.eta(), track1.phi());
    }
    fillQA(collision, multiplicity, tracks2);
  }

  template <typename TTarget>
  bool fillContainerEvent(TTarget target, float multiplicity, CorrelationContainer::CFStep step)
  {
    const float containerMultiplicity = getCorrelationContainerMultiplicity(multiplicity);
    if (containerMultiplicity < 0.f) {
      return false;
    }
    target->fillEvent(containerMultiplicity, step);
    return true;
  }

  template <typename TTarget, typename TCollision>
  bool fillCollisionAOD(const TTarget& target, const TCollision& collision, const float& multiplicity)
  {
    if (!fillContainerEvent(target, multiplicity, CorrelationContainer::kCFStepAll)) {
      return false;
    }

    if (!collision.alias_bit(kINT7) || !collision.sel7()) {
      return false;
    }

    fillContainerEvent(target, multiplicity, CorrelationContainer::kCFStepReconstructed);

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool checkObject(TTrack& track)
  {
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
      // If using MC trigger PDGs, allow ONLY those PDGs to bypass isPhysicalPrimary
      if (!cfgMcTriggerPDGs->empty()) {
        // track has pdgCode in this compilation branch (you only call checkObject where that is true)
        const bool isWantedTrigger =
          std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track.pdgCode()) != cfgMcTriggerPDGs->end();
        if (isWantedTrigger) {
          return true; // allow phi, K*, etc. even if not physical primary
        }
        // For everything else keep original definition
        return track.isPhysicalPrimary();
      }
      return track.isPhysicalPrimary();
    } else if constexpr (step == CorrelationContainer::kCFStepTrackedOnlyPrim) {
      return track.isPhysicalPrimary() && (track.flags() & aod::cfmcparticle::kReconstructed);
    } else if constexpr (step == CorrelationContainer::kCFStepTracked) {
      return (track.flags() & aod::cfmcparticle::kReconstructed);
    }

    return true;
  }

  template <class T>
  using HasSign = decltype(std::declval<T&>().sign());
  template <class T>
  using HasDecay = decltype(std::declval<T&>().decay());
  template <class T>
  using HasPartDaugh0Id = decltype(std::declval<T&>().cfParticleDaugh0Id());
  template <class T>
  using HasPartDaugh1Id = decltype(std::declval<T&>().cfParticleDaugh1Id());

  template <typename TParticle>
  int getParticleSign(const TParticle& particle)
  {
    if constexpr (std::experimental::is_detected<HasSign, TParticle>::value) {
      return particle.sign();
    } else if constexpr (std::experimental::is_detected<HasPDGCode, TParticle>::value) {
      const auto* pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (pdgParticle) {
        return (pdgParticle->Charge() > 0.0) - (pdgParticle->Charge() < 0.0);
      }
    }
    return 0;
  }

  template <class CollType>
  bool passOutlier(CollType const& collision)
  {
    if (cfgMultCutFormula.value.empty()) {
      return true;
    }
    for (uint i = 0; i < aod::cfmultset::NMultiplicityEstimators; ++i) {
      if ((cfgMultCorrelationsMask.value & (1u << i)) == 0 || multCutFormulaParamIndex[i] == ~0u) {
        continue;
      }
      auto estIndex = std::popcount(cfgMultCorrelationsMask.value & ((1u << i) - 1));
      multCutFormula->SetParameter(multCutFormulaParamIndex[i], collision.multiplicities()[estIndex]);
    }
    return multCutFormula->Eval() > 0.0f;
  }

  template <typename T>
  std::tuple<bool, float> getV0Rapidity(const T& track)
  {
    if constexpr (!std::experimental::is_detected<HasDecay, T>::value) {
      return {false, 0.0f}; // no decay type, return dummy rapidity
    }
    const auto decayType = track.decay();
    float mass = 0.f;

    if (decayType == aod::cf2prongtrack::K0stoPiPi) {
      mass = o2::constants::physics::MassK0Short;
    } else if (decayType == aod::cf2prongtrack::LambdatoPPi || decayType == aod::cf2prongtrack::AntiLambdatoPiP) {
      mass = o2::constants::physics::MassLambda;
    } else if (decayType == aod::cf2prongtrack::PhiToKKPID1 || decayType == aod::cf2prongtrack::PhiToKKPID2 || decayType == aod::cf2prongtrack::PhiToKKPID3) {
      mass = o2::constants::physics::MassPhi;
    } else {
      return {false, 0.0f}; // unsupported decay type, return dummy rapidity
    }

    const float pt = track.pt();
    const float eta = track.eta();
    const float phi = track.phi();

    const float px = pt * std::cos(phi);
    const float py = pt * std::sin(phi);
    const float pz = pt * std::sinh(eta);

    const float p2 = px * px + py * py + pz * pz;

    const float E = std::sqrt(p2 + mass * mass); // o2-linter: disable=name/function-variable (E is standard variable name for energy)
    return {true, 0.5f * std::log((E + pz) / (E - pz))};
  }

  void loadYieldTemplatesFromTree(TTree* tree, const TNamed* schemaVersion, const TNamed* correlationStep, const std::string& source)
  {
    if (!tree) {
      LOGF(fatal, "Missing ensembleYieldTemplates in %s", source.c_str());
      return;
    }
    if (schemaVersion == nullptr || TString(schemaVersion->GetTitle()) != "1") {
      LOGF(fatal, "Unsupported or missing ensemble-yield template schema version in %s", source.c_str());
      return;
    }
    if (correlationStep == nullptr || TString(correlationStep->GetTitle()) != "kCFStepReconstructed") {
      LOGF(fatal, "Ensemble-yield templates must be derived at kCFStepReconstructed");
      return;
    }

    YieldTemplate value;
    int fitStatus = -1;
    std::array<double, 10> parameters{};
    tree->SetBranchAddress("trigBin", &value.trigBin);
    tree->SetBranchAddress("assocBin", &value.assocBin);
    tree->SetBranchAddress("multBin", &value.multBin);
    tree->SetBranchAddress("nchLow", &value.nchLow);
    tree->SetBranchAddress("nchHigh", &value.nchHigh);
    tree->SetBranchAddress("trigPtLow", &value.trigPtLow);
    tree->SetBranchAddress("trigPtHigh", &value.trigPtHigh);
    tree->SetBranchAddress("assocPtLow", &value.assocPtLow);
    tree->SetBranchAddress("assocPtHigh", &value.assocPtHigh);
    tree->SetBranchAddress("parameters", parameters.data());
    tree->SetBranchAddress("fitStatus", &fitStatus);

    yieldTemplates.clear();
    yieldTemplates.reserve(tree->GetEntries());
    for (int64_t iEntry = 0; iEntry < tree->GetEntries(); ++iEntry) {
      tree->GetEntry(iEntry);
      if (fitStatus != 0) {
        LOGF(warning, "Skipping failed yield template (%d, %d, %d), fit status %d", value.trigBin, value.assocBin, value.multBin, fitStatus);
        continue;
      }
      value.parameters = parameters;
      if (value.parameters[2] <= 0.0 || value.parameters[5] <= 0.0 || value.parameters[8] <= 0.0) {
        LOGF(warning, "Skipping yield template (%d, %d, %d) with non-positive Gaussian width", value.trigBin, value.assocBin, value.multBin);
        continue;
      }
      const auto intervalMatchesAxis = [](const AxisSpec& axis, double low, double high) {
        constexpr double Tolerance = 1e-6;
        const auto& edges = axis.binEdges;
        for (std::size_t index = 0; index + 1 < edges.size(); ++index) {
          if (std::abs(edges[index] - low) < Tolerance && std::abs(edges[index + 1] - high) < Tolerance) {
            return true;
          }
        }
        return false;
      };
      if (!intervalMatchesAxis(AxisSpec(axisMultiplicity), value.nchLow, value.nchHigh) ||
          !intervalMatchesAxis(AxisSpec(axisPtTrigger), value.trigPtLow, value.trigPtHigh) ||
          !intervalMatchesAxis(AxisSpec(axisPtAssoc), value.assocPtLow, value.assocPtHigh)) {
        LOGF(fatal, "Template binning does not match the configured axes for template (%d, %d, %d)", value.trigBin, value.assocBin, value.multBin);
        tree->ResetBranchAddresses();
        return;
      }
      yieldTemplates.push_back(value);
    }
    tree->ResetBranchAddresses();

    if (yieldTemplates.empty()) {
      LOGF(fatal, "No valid ensemble-yield templates were loaded from %s", source.c_str());
      return;
    }
    LOGF(info, "Loaded %zu ensemble-yield templates from %s", yieldTemplates.size(), source.c_str());
  }

  void loadPairAcceptanceMaps(const std::function<TObject*(const char*)>& findObject, const std::string& source)
  {
    const auto* schemaVersion = dynamic_cast<const TNamed*>(findObject("pairAcceptanceSchemaVersion"));
    const auto* normalization = dynamic_cast<const TNamed*>(findObject("pairAcceptanceNormalization"));
    const auto* axes = dynamic_cast<const TNamed*>(findObject("pairAcceptanceAxes"));
    const TString schema = schemaVersion != nullptr ? schemaVersion->GetTitle() : "";
    if (schemaVersion == nullptr || (schema != "2" && schema != "3") || normalization == nullptr || axes == nullptr) {
      LOGF(fatal, "Unsupported or missing pair-acceptance metadata in %s", source.c_str());
      return;
    }

    const int nMultiplicityBins = AxisSpec(axisMultiplicity).getNbins();
    pairAcceptanceMaps.clear();
    pairAcceptanceEtaVertexMaps.clear();
    pairAcceptanceSchemaVersion = schema.Atoi();
    int multiplicityDependentOnly = 2;
    if (pairAcceptanceSchemaVersion == multiplicityDependentOnly) {
      pairAcceptanceMaps.resize(nMultiplicityBins);
    } else {
      pairAcceptanceEtaVertexMaps.resize(nMultiplicityBins);
    }
    for (int multBin = 0; multBin < nMultiplicityBins; ++multBin) {
      if (pairAcceptanceSchemaVersion == multiplicityDependentOnly) {
        auto* inputMap = dynamic_cast<TH3D*>(findObject(Form("pairAcceptance_mult_%d", multBin)));
        auto* clone = inputMap != nullptr ? dynamic_cast<TH3D*>(inputMap->Clone(Form("loadedPairAcceptance_mult_%d", multBin))) : nullptr;
        if (clone == nullptr) {
          LOGF(fatal, "Missing or invalid pairAcceptance_mult_%d in %s", multBin, source.c_str());
          pairAcceptanceMaps.clear();
          return;
        }
        clone->SetDirectory(nullptr);
        pairAcceptanceMaps[multBin].reset(clone);
      } else {
        auto* inputMap = dynamic_cast<TH2D*>(findObject(Form("pairAcceptanceEtaVertex_mult_%d", multBin)));
        auto* clone = inputMap != nullptr ? dynamic_cast<TH2D*>(inputMap->Clone(Form("loadedPairAcceptanceEtaVertex_mult_%d", multBin))) : nullptr;
        if (clone == nullptr) {
          LOGF(fatal, "Missing or invalid pairAcceptanceEtaVertex_mult_%d in %s", multBin, source.c_str());
          pairAcceptanceEtaVertexMaps.clear();
          return;
        }
        clone->SetDirectory(nullptr);
        pairAcceptanceEtaVertexMaps[multBin].reset(clone);
      }
    }
    LOGF(info, "Loaded %d schema-%d pair-acceptance maps from %s", nMultiplicityBins, pairAcceptanceSchemaVersion, source.c_str());
  }

  void loadLocalYieldTemplates()
  {
    if (cfgNuncSeedsTemplateFile.value.empty()) {
      return;
    }

    std::unique_ptr<TFile> input(TFile::Open(cfgNuncSeedsTemplateFile.value.c_str(), "READ"));
    if (!input || input->IsZombie()) {
      LOGF(fatal, "Could not open ensemble-yield template file: %s", cfgNuncSeedsTemplateFile.value.c_str());
      return;
    }

    auto* calibration = dynamic_cast<TList*>(input->Get("ccdb_object"));
    auto findObject = [&](const char* name) -> TObject* {
      if (auto* object = input->Get(name)) {
        return object;
      }
      return calibration ? calibration->FindObject(name) : nullptr;
    };
    loadYieldTemplatesFromTree(
      dynamic_cast<TTree*>(findObject("ensembleYieldTemplates")),
      dynamic_cast<TNamed*>(findObject("ensembleYieldTemplateSchemaVersion")),
      dynamic_cast<TNamed*>(findObject("ensembleYieldTemplateCorrelationStep")),
      cfgNuncSeedsTemplateFile.value);
    loadPairAcceptanceMaps(findObject, cfgNuncSeedsTemplateFile.value);
  }

  void loadCcdbYieldTemplates(uint64_t timestamp)
  {
    if (cfgNuncSeedsTemplate.value.empty()) {
      return;
    }

    auto* calibration = ccdb->getForTimeStamp<TList>(cfgNuncSeedsTemplate.value, timestamp);
    if (!calibration) {
      LOGF(fatal, "Could not load ensemble-yield templates from CCDB path %s at timestamp %llu", cfgNuncSeedsTemplate.value.c_str(), timestamp);
      return;
    }
    if (calibration == loadedCcdbYieldTemplateObject) {
      return;
    }

    loadYieldTemplatesFromTree(
      dynamic_cast<TTree*>(calibration->FindObject("ensembleYieldTemplates")),
      dynamic_cast<TNamed*>(calibration->FindObject("ensembleYieldTemplateSchemaVersion")),
      dynamic_cast<TNamed*>(calibration->FindObject("ensembleYieldTemplateCorrelationStep")),
      cfgNuncSeedsTemplate.value);
    loadPairAcceptanceMaps(
      [calibration](const char* name) -> TObject* { return calibration->FindObject(name); },
      cfgNuncSeedsTemplate.value);
    loadedCcdbYieldTemplateObject = calibration;
  }

  const YieldTemplate* findYieldTemplate(double multiplicity, double trigPt, double assocPt) const
  {
    for (const auto& yieldTemplate : yieldTemplates) {
      if (multiplicity >= yieldTemplate.nchLow && multiplicity < yieldTemplate.nchHigh &&
          trigPt >= yieldTemplate.trigPtLow && trigPt < yieldTemplate.trigPtHigh &&
          assocPt >= yieldTemplate.assocPtLow && assocPt < yieldTemplate.assocPtHigh) {
        return &yieldTemplate;
      }
    }
    return nullptr;
  }

  bool hasTriggerTemplate(double multiplicity, double trigPt) const
  {
    return std::any_of(yieldTemplates.begin(), yieldTemplates.end(), [multiplicity, trigPt](const auto& yieldTemplate) {
      return multiplicity >= yieldTemplate.nchLow && multiplicity < yieldTemplate.nchHigh &&
             trigPt >= yieldTemplate.trigPtLow && trigPt < yieldTemplate.trigPtHigh;
    });
  }

  bool hasMultiplicityTemplate(double multiplicity) const
  {
    return std::any_of(yieldTemplates.begin(), yieldTemplates.end(), [multiplicity](const auto& yieldTemplate) {
      return multiplicity >= yieldTemplate.nchLow && multiplicity < yieldTemplate.nchHigh;
    });
  }

  static double evaluateGaussian(double x, const double* parameters)
  {
    const double pull = (x - parameters[1]) / parameters[2];
    return parameters[0] * std::exp(-0.5 * pull * pull);
  }

  static std::array<double, 3> getTemplateComponentIntegrals(const YieldTemplate& yieldTemplate)
  {
    const auto& p = yieldTemplate.parameters;
    const double gaussianIntegralFactor = std::sqrt(o2::constants::math::TwoPI);
    return {
      std::max(0.0, gaussianIntegralFactor * (p[0] * p[2] + p[3] * p[5])),
      std::max(0.0, gaussianIntegralFactor * p[6] * p[8]),
      std::max(0.0, o2::constants::math::TwoPI * p[9])};
  }

  void addTriggerPriorExpectations(EventSeedEstimate& estimate, double multiplicity, double trigPt)
  {
    if (estimate.templatePriorCounts.empty()) {
      estimate.templatePriorCounts.resize(yieldTemplates.size());
    }
    for (std::size_t templateIndex = 0; templateIndex < yieldTemplates.size(); ++templateIndex) {
      const auto& yieldTemplate = yieldTemplates[templateIndex];
      if (multiplicity < yieldTemplate.nchLow || multiplicity >= yieldTemplate.nchHigh ||
          trigPt < yieldTemplate.trigPtLow || trigPt >= yieldTemplate.trigPtHigh) {
        continue;
      }
      const auto componentIntegrals = getTemplateComponentIntegrals(yieldTemplate);
      for (std::size_t component = 0; component < componentIntegrals.size(); ++component) {
        // The fit integrals are per-trigger yields. Summing one copy for every
        // templated trigger makes the Gamma-prior mode follow the event's actual
        // trigger-pT composition without using its observed pair assignments.
        estimate.templatePriorCounts[templateIndex][component] += componentIntegrals[component];
        estimate.priorComponentCounts[component] += componentIntegrals[component];
      }
    }
  }

  int findPairAcceptanceMultiplicityBin(double multiplicity) const
  {
    const auto& edges = AxisSpec(axisMultiplicity).binEdges;
    const auto upper = std::upper_bound(edges.begin(), edges.end(), multiplicity);
    return static_cast<int>(std::distance(edges.begin(), upper)) - 1;
  }

  double getEventSeedUserAxisValue(double multiplicity, double eventSeed) const
  {
    if (!eventSeedPercentileAxisEnabled || eventSeed < 0.0) {
      return eventSeed;
    }
    const int multBin = findPairAcceptanceMultiplicityBin(multiplicity);
    if (multBin < 0 || multBin >= static_cast<int>(cfgEventSeedPercentileLower->size())) {
      return -1.0;
    }
    if (eventSeed < cfgEventSeedPercentileLower->at(multBin)) {
      return 0.0;
    }
    if (eventSeed >= cfgEventSeedPercentileUpper->at(multBin)) {
      return 2.0;
    }
    return 1.0;
  }

  double getPairAcceptance(double multiplicity, double deltaPhi, double deltaEta, double posZ) const
  {
    const int multBin = findPairAcceptanceMultiplicityBin(multiplicity);
    int etaVertexMultiplicityDependentOnly = 3;
    if (pairAcceptanceSchemaVersion == etaVertexMultiplicityDependentOnly) {
      if (multBin < 0 || multBin >= static_cast<int>(pairAcceptanceEtaVertexMaps.size()) || pairAcceptanceEtaVertexMaps[multBin] == nullptr) {
        return 0.0;
      }
      const auto* map = pairAcceptanceEtaVertexMaps[multBin].get();
      const int etaBin = map->GetXaxis()->FindFixBin(deltaEta);
      const int vertexBin = map->GetYaxis()->FindFixBin(posZ);
      if (etaBin < 1 || etaBin > map->GetNbinsX() || vertexBin < 1 || vertexBin > map->GetNbinsY()) {
        return 0.0;
      }
      return map->GetBinContent(etaBin, vertexBin);
    }
    if (multBin < 0 || multBin >= static_cast<int>(pairAcceptanceMaps.size()) || pairAcceptanceMaps[multBin] == nullptr) {
      return 0.0;
    }
    const auto* map = pairAcceptanceMaps[multBin].get();
    const int phiBin = map->GetXaxis()->FindFixBin(deltaPhi);
    const int etaBin = map->GetYaxis()->FindFixBin(deltaEta);
    const int vertexBin = map->GetZaxis()->FindFixBin(posZ);
    if (phiBin < 1 || phiBin > map->GetNbinsX() || etaBin < 1 || etaBin > map->GetNbinsY() ||
        vertexBin < 1 || vertexBin > map->GetNbinsZ()) {
      return 0.0;
    }
    return map->GetBinContent(phiBin, etaBin, vertexBin);
  }

  void addPairProbabilities(EventSeedEstimate& estimate, const YieldTemplate& yieldTemplate,
                            double multiplicity, double deltaPhi, double deltaEta, double posZ, bool fillAcceptanceQA = true)
  {
    const auto& parameters = yieldTemplate.parameters;
    const double near = std::max(0.0, evaluateGaussian(deltaPhi, parameters.data()) + evaluateGaussian(deltaPhi, parameters.data() + 3));
    const double away = std::max(0.0, evaluateGaussian(deltaPhi, parameters.data() + 6));
    const double baseline = std::max(0.0, parameters[9]);
    const double total = baseline + near + away;
    if (total <= 0.0) {
      return;
    }
    const double nearProbability = near / total;
    const double awayProbability = away / total;
    const double baselineProbability = baseline / total;
    estimate.rawNearPairs += nearProbability;
    estimate.rawAwayPairs += awayProbability;
    estimate.rawBaselinePairs += baselineProbability;
    ++estimate.nPairs;

    const double acceptance = getPairAcceptance(multiplicity, deltaPhi, deltaEta, posZ);
    if (!std::isfinite(acceptance) || acceptance < cfgMinPairAcceptance || acceptance <= 0.0) {
      ++estimate.nPairsRejectedByAcceptance;
      return;
    }
    const double acceptanceWeight = 1.0 / acceptance;
    if (!std::isfinite(acceptanceWeight) || acceptanceWeight > cfgMaxPairAcceptanceWeight) {
      ++estimate.nPairsRejectedByAcceptance;
      return;
    }
    estimate.probabilityNearPairs += acceptanceWeight * nearProbability;
    estimate.probabilityAwayPairs += acceptanceWeight * awayProbability;
    estimate.probabilityBaselinePairs += acceptanceWeight * baselineProbability;
    // Method 0 uses these values directly. Method 1 overwrites them with the
    // event-wide MAP-EM component yields after all pairs have been collected.
    estimate.nearPairs = estimate.probabilityNearPairs;
    estimate.awayPairs = estimate.probabilityAwayPairs;
    estimate.baselinePairs = estimate.probabilityBaselinePairs;

    if (cfgEventSeedEstimatorMethod == 1) {
      const auto componentIntegrals = getTemplateComponentIntegrals(yieldTemplate);
      EventSeedPairObservation observation;
      observation.templateIndex = static_cast<std::size_t>(&yieldTemplate - yieldTemplates.data());
      observation.normalizedShapes = {
        componentIntegrals[0] > 0.0 ? near / componentIntegrals[0] : 0.0,
        componentIntegrals[1] > 0.0 ? away / componentIntegrals[1] : 0.0,
        componentIntegrals[2] > 0.0 ? baseline / componentIntegrals[2] : 0.0};
      observation.acceptanceWeight = acceptanceWeight;
      estimate.pairObservations.push_back(observation);
    }
    estimate.sumAcceptanceWeights += acceptanceWeight;
    ++estimate.nAcceptanceCorrectedPairs;
    if (fillAcceptanceQA) {
      registry.fill(HIST("eventSeedAcceptanceWeight"), acceptanceWeight);
    }
  }

  void finalizeEventSeedEstimate(EventSeedEstimate& estimate)
  {
    if (cfgEventSeedEstimatorMethod == 0 || !estimate.hasAcceptanceCorrectedPairs()) {
      return;
    }

    estimate.usedMapEM = true;
    auto componentYields = estimate.priorComponentCounts;
    const double priorExposure = cfgEventSeedPriorExposure;
    constexpr double Tiny = 1.e-15;
    bool hasPositivePrior = false;
    for (const auto& value : componentYields) {
      hasPositivePrior = hasPositivePrior || value > 0.0;
    }
    if (!hasPositivePrior) {
      return;
    }

    for (int iteration = 0; iteration < cfgEventSeedEMMaxIterations; ++iteration) {
      std::array<double, 3> weightedComponentCounts{};

      for (const auto& observation : estimate.pairObservations) {
        if (observation.templateIndex >= estimate.templatePriorCounts.size()) {
          continue;
        }
        const auto& templatePrior = estimate.templatePriorCounts[observation.templateIndex];
        std::array<double, 3> numerators{};
        double denominator = 0.0;
        for (std::size_t component = 0; component < numerators.size(); ++component) {
          if (estimate.priorComponentCounts[component] <= 0.0) {
            continue;
          }
          // templatePrior/priorComponentCounts is the ensemble-predicted pT-bin
          // fraction for this component. At the prior mode, the responsibility
          // therefore reduces to the existing fixed-template probability.
          const double templateFraction = templatePrior[component] / estimate.priorComponentCounts[component];
          numerators[component] = componentYields[component] * templateFraction * observation.normalizedShapes[component];
          denominator += numerators[component];
        }
        if (!std::isfinite(denominator) || denominator <= 0.0) {
          continue;
        }
        for (std::size_t component = 0; component < numerators.size(); ++component) {
          weightedComponentCounts[component] += observation.acceptanceWeight * numerators[component] / denominator;
        }
      }

      std::array<double, 3> updatedYields{};
      double maximumRelativeChange = 0.0;
      for (std::size_t component = 0; component < updatedYields.size(); ++component) {
        // Gamma(shape=1+s*m, rate=s) has mode m. Combining it with one
        // event exposure gives the closed MAP update below. With s=0 this is
        // the unregularized weighted EM update.
        updatedYields[component] =
          (priorExposure * estimate.priorComponentCounts[component] + weightedComponentCounts[component]) /
          (priorExposure + 1.0);
        const double scale = std::max({std::abs(componentYields[component]), std::abs(updatedYields[component]), Tiny});
        maximumRelativeChange = std::max(maximumRelativeChange, std::abs(updatedYields[component] - componentYields[component]) / scale);
      }
      componentYields = updatedYields;
      estimate.emIterations = iteration + 1;
      if (maximumRelativeChange < cfgEventSeedEMTolerance) {
        estimate.emConverged = true;
        break;
      }
    }

    estimate.nearPairs = componentYields[0];
    estimate.awayPairs = componentYields[1];
    estimate.baselinePairs = componentYields[2];
  }

  void fillEventSeedEstimatorQA(double multiplicity, const EventSeedEstimate& estimate)
  {
    if (yieldTemplates.empty()) {
      return;
    }
    if (!hasMultiplicityTemplate(multiplicity)) {
      return;
    }
    registry.fill(HIST("profileEventNTriggers"), multiplicity, estimate.nTriggers);
    registry.fill(HIST("profileEventNearPairs"), multiplicity, estimate.nearPairs);
    registry.fill(HIST("profileEventAwayPairs"), multiplicity, estimate.awayPairs);
    registry.fill(HIST("profileEventBaselinePairs"), multiplicity, estimate.baselinePairs);
    registry.fill(HIST("profileEventRawNearPairs"), multiplicity, estimate.rawNearPairs);
    registry.fill(HIST("profileEventRawAwayPairs"), multiplicity, estimate.rawAwayPairs);
    registry.fill(HIST("profileEventRawBaselinePairs"), multiplicity, estimate.rawBaselinePairs);
    const double acceptanceCoverage = estimate.nPairs > 0 ? static_cast<double>(estimate.nAcceptanceCorrectedPairs) / estimate.nPairs : 0.0;
    const double meanAcceptanceWeight = estimate.nAcceptanceCorrectedPairs > 0 ? estimate.sumAcceptanceWeights / estimate.nAcceptanceCorrectedPairs : 0.0;
    registry.fill(HIST("profileEventAcceptanceCoverage"), multiplicity, acceptanceCoverage);
    registry.fill(HIST("profileEventMeanAcceptanceWeight"), multiplicity, meanAcceptanceWeight);
    registry.fill(HIST("profileEventEstimatorValidity"), multiplicity, estimate.isValid() ? 1.0 : 0.0);
    if (!estimate.hasTriggers()) {
      registry.fill(HIST("eventSeedEstimatorStatus"), 1.0);
      return;
    }
    const double coverage = estimate.nCandidatePairs > 0 ? static_cast<double>(estimate.nPairs) / estimate.nCandidatePairs : 0.0;
    registry.fill(HIST("eventSeedTemplateCoverage"), multiplicity, coverage);
    if (estimate.nPairs == 0) {
      registry.fill(HIST("eventSeedEstimatorStatus"), 2.0);
      return;
    }
    if (!estimate.hasAcceptanceCorrectedPairs()) {
      registry.fill(HIST("eventSeedEstimatorStatus"), 3.0);
      return;
    }
    registry.fill(HIST("eventSeedEstimatorStatus"), 4.0);
    registry.fill(HIST("eventSeedEstimator"), multiplicity, estimate.nTriggers, estimate.nearYield(), estimate.awayYield(), estimate.nuncSeeds());
    registry.fill(HIST("eventSeedPairProbabilities"), estimate.probabilityBaselinePairs, estimate.probabilityNearPairs, estimate.probabilityAwayPairs);
    registry.fill(HIST("eventSeedEstimateVsMultiplicity"), multiplicity, estimate.nuncSeeds());
    registry.fill(HIST("eventSeedFixedProbabilityVsMultiplicity"), multiplicity, estimate.probabilityNuncSeeds());
    registry.fill(HIST("eventNearYieldVsMultiplicity"), multiplicity, estimate.nearYield());
    registry.fill(HIST("eventAwayYieldVsMultiplicity"), multiplicity, estimate.awayYield());
    registry.fill(HIST("profileEventNearYield"), multiplicity, estimate.nearYield());
    registry.fill(HIST("profileEventAwayYield"), multiplicity, estimate.awayYield());
    registry.fill(HIST("profileEventNuncSeeds"), multiplicity, estimate.nuncSeeds());
    if (estimate.usedMapEM) {
      registry.fill(HIST("profileEventProbabilityNuncSeeds"), multiplicity, estimate.probabilityNuncSeeds());
      registry.fill(HIST("eventSeedEMIterations"), estimate.emIterations);
      registry.fill(HIST("profileEventEMConvergence"), multiplicity, estimate.emConverged ? 1.0 : 0.0);
      registry.fill(HIST("profileEventEMNearPrior"), multiplicity, estimate.priorComponentCounts[0]);
      registry.fill(HIST("profileEventEMAwayPrior"), multiplicity, estimate.priorComponentCounts[1]);
      registry.fill(HIST("profileEventEMBaselinePrior"), multiplicity, estimate.priorComponentCounts[2]);
      registry.fill(HIST("eventSeedMAPVsProbability"), estimate.probabilityNuncSeeds(), estimate.nuncSeeds());
      registry.fill(HIST("eventSeedMAPEMVsMultiplicity"), multiplicity, estimate.nuncSeeds());
    }
  }

  void fillMCValidation(double multiplicity, const EventSeedEstimate& estimate, int trueNMPI)
  {
    if (!eventSeedEstimatorEnabled || yieldTemplates.empty()) {
      return;
    }
    if (trueNMPI < 0) {
      registry.fill(HIST("mcValidation/status"), 1.0);
      return;
    }
    registry.fill(HIST("mcValidation/trueNMPIVsMultiplicity"), multiplicity, trueNMPI);
    if (!hasMultiplicityTemplate(multiplicity)) {
      registry.fill(HIST("mcValidation/status"), 2.0);
      return;
    }
    if (!estimate.hasTriggers()) {
      registry.fill(HIST("mcValidation/status"), 3.0);
      return;
    }
    const double coverage = estimate.nCandidatePairs > 0 ? static_cast<double>(estimate.nPairs) / estimate.nCandidatePairs : 0.0;
    registry.fill(HIST("mcValidation/templateCoverageVsTrueNMPI"), trueNMPI, coverage);
    if (estimate.nPairs == 0) {
      registry.fill(HIST("mcValidation/status"), 4.0);
      return;
    }
    if (!estimate.hasAcceptanceCorrectedPairs()) {
      registry.fill(HIST("mcValidation/status"), 5.0);
      return;
    }

    const double estimatedSeeds = estimate.nuncSeeds();
    const double bias = estimatedSeeds - trueNMPI;
    registry.fill(HIST("mcValidation/status"), 6.0);
    registry.fill(HIST("mcValidation/estimatedSeedsVsTrueNMPI"), trueNMPI, estimatedSeeds);
    registry.fill(HIST("mcValidation/estimatedSeedsVsTrueNMPIVsMultiplicity"), multiplicity, trueNMPI, estimatedSeeds);
    if (estimate.usedMapEM) {
      registry.fill(HIST("mcValidation/probabilityEstimatedSeedsVsTrueNMPI"), trueNMPI, estimate.probabilityNuncSeeds());
    }
    registry.fill(HIST("mcValidation/profileEstimatedSeedsVsTrueNMPI"), trueNMPI, estimatedSeeds);
    registry.fill(HIST("mcValidation/profileBiasVsTrueNMPI"), trueNMPI, bias);
    if (trueNMPI > 0) {
      registry.fill(HIST("mcValidation/relativeResidualVsTrueNMPI"), trueNMPI, bias / trueNMPI);
    }
  }

  void fillGeneratedMCValidation(double multiplicity, const EventSeedEstimate& estimate, int trueNMPI)
  {
    if (!eventSeedEstimatorEnabled || yieldTemplates.empty()) {
      return;
    }
    if (trueNMPI < 0) {
      registry.fill(HIST("mcValidation/generated/status"), 0.0);
      return;
    }
    registry.fill(HIST("mcValidation/generated/trueNMPIVsMultiplicity"), multiplicity, trueNMPI);
    if (!hasMultiplicityTemplate(multiplicity)) {
      registry.fill(HIST("mcValidation/generated/status"), 1.0);
      return;
    }
    if (!estimate.hasTriggers()) {
      registry.fill(HIST("mcValidation/generated/status"), 2.0);
      return;
    }
    const double coverage = estimate.nCandidatePairs > 0 ? static_cast<double>(estimate.nPairs) / estimate.nCandidatePairs : 0.0;
    registry.fill(HIST("mcValidation/generated/templateCoverageVsTrueNMPI"), trueNMPI, coverage);
    if (estimate.nPairs == 0) {
      registry.fill(HIST("mcValidation/generated/status"), 3.0);
      return;
    }
    if (!estimate.hasAcceptanceCorrectedPairs()) {
      registry.fill(HIST("mcValidation/generated/status"), 4.0);
      return;
    }

    const double estimatedSeeds = estimate.nuncSeeds();
    const double bias = estimatedSeeds - trueNMPI;
    registry.fill(HIST("mcValidation/generated/status"), 5.0);
    registry.fill(HIST("mcValidation/generated/estimatedSeedsVsTrueNMPI"), trueNMPI, estimatedSeeds);
    registry.fill(HIST("mcValidation/generated/estimatedSeedsVsTrueNMPIVsMultiplicity"), multiplicity, trueNMPI, estimatedSeeds);
    if (estimate.usedMapEM) {
      registry.fill(HIST("mcValidation/generated/probabilityEstimatedSeedsVsTrueNMPI"), trueNMPI, estimate.probabilityNuncSeeds());
    }
    registry.fill(HIST("mcValidation/generated/profileEstimatedSeedsVsTrueNMPI"), trueNMPI, estimatedSeeds);
    registry.fill(HIST("mcValidation/generated/profileBiasVsTrueNMPI"), trueNMPI, bias);
    if (trueNMPI > 0) {
      registry.fill(HIST("mcValidation/generated/relativeResidualVsTrueNMPI"), trueNMPI, bias / trueNMPI);
    }
  }

  template <CorrelationContainer::CFStep step, typename TTarget>
  void flushSeedAxisFills(TTarget target, const std::vector<PendingSeedTriggerFill>& triggerFills, const std::vector<PendingSeedPairFill>& pairFills, double eventSeed)
  {
    for (const auto& fill : triggerFills) {
      target->getTriggerHist()->Fill(step, fill.pt, fill.multiplicity, fill.posZ, getEventSeedUserAxisValue(fill.multiplicity, eventSeed), fill.weight);
    }
    for (const auto& fill : pairFills) {
      target->getPairHist()->Fill(step, fill.deltaEta, fill.assocPt, fill.triggerPt, fill.multiplicity, fill.deltaPhi, fill.posZ, getEventSeedUserAxisValue(fill.multiplicity, eventSeed), fill.weight);
    }
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks>
  double estimateEventSeedWithoutFilling(TTarget target, TTracks& tracks, float multiplicity, float posZ, int magField)
  {
    EventSeedEstimate estimate;
    std::vector<PendingSeedTriggerFill> discardedTriggerFills;
    std::vector<PendingSeedPairFill> discardedPairFills;
    discardedTriggerFills.reserve(tracks.size());
    discardedPairFills.reserve(tracks.size() * tracks.size());
    fillCorrelations<step>(target, tracks, tracks, multiplicity, posZ, magField, 1.0f, &estimate, &discardedTriggerFills, &discardedPairFills, -1.f, false, false);
    finalizeEventSeedEstimate(estimate);
    return estimate.nuncSeeds();
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelations(TTarget target, TTracks1& tracks1, TTracks2& tracks2, float multiplicity, float posZ, int magField, float eventWeight, EventSeedEstimate* seedEstimate = nullptr, std::vector<PendingSeedTriggerFill>* pendingTriggerFills = nullptr, std::vector<PendingSeedPairFill>* pendingPairFills = nullptr, double eventSeed = -1.0, bool fillEstimatorAcceptanceQA = true, bool fillLoopQA = true)
  {
    const float containerMultiplicity = getCorrelationContainerMultiplicity(multiplicity);
    if (containerMultiplicity < 0.f) {
      return;
    }

    // Cache efficiency for particles (too many FindBin lookups)
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      if (cfg.mEfficiencyAssociated) {
        efficiencyAssociatedCache.clear();
        efficiencyAssociatedCache.reserve(tracks2.size());
        for (const auto& track : tracks2) {
          efficiencyAssociatedCache.push_back(getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ));
        }
      }
    }

    for (const auto& track1 : tracks1) {
      // LOGF(info, "Track %f | %f | %f  %d %d", track1.eta(), track1.phi(), track1.pt(), track1.isGlobalTrack(), track1.isGlobalTrackSDD());

      if constexpr (step <= CorrelationContainer::kCFStepTracked && !std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
        if (!checkObject<step>(track1)) {
          continue;
        }
      }

      // sign check and PDG code special cases
      if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
        // If the MC trigger particle is on the trigger PDG code list, we will accept them regardless of their charge.
        if (!cfgMcTriggerPDGs->empty()) {
          if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track1.pdgCode()) == cfgMcTriggerPDGs->end()) {
            continue;
          }
        } else { // otherwise check the sign against the configuration
          const int sign = getParticleSign(track1);
          if (cfgTriggerCharge != 0) {
            if (cfgTriggerCharge * sign < 0) {
              continue;
            }
          } else if (sign == 0) {
            continue; // reject neutral MC particles
          }
        }
      } else if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value) {
        // Check reco objects that have the sign attribute. There are no neutrals to deal with.
        if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
          continue;
        }
      }

      if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
        if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track1.decay()))) == 0u) {
          continue; // skip particles that do not match the decay mask
        }
        if (cfgV0RapidityMax > 0) {
          auto [t, y] = getV0Rapidity(track1);
          if (t && std::abs(y) > cfgV0RapidityMax) {
            continue; // V0s are not allowed to be outside the rapidity range
          }
          if (fillLoopQA) {
            registry.fill(HIST("yvspt"), y, track1.pt());
          }
        }
      }

      if constexpr (std::experimental::is_detected<HasPartDaugh0Id, typename TTracks1::iterator>::value) {
        if (track1.cfParticleDaugh0Id() < 0 && track1.cfParticleDaugh1Id() < 0) {
          continue; // these we could not match
        }
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ);
        }
      }

      if (cfgUserAxis == EventSeedAxis) {
        if (pendingTriggerFills) {
          pendingTriggerFills->push_back({track1.pt(), containerMultiplicity, posZ, triggerWeight});
        } else {
          target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, eventSeed, triggerWeight);
        }
      } else if (cfgUserAxis == InvariantMassAxis) {
        if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value) {
          target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, track1.invMass(), triggerWeight);
        } else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
          // TParticlePDG *p = pdg->GetParticle(track1.pdgCode());
          // target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, p->Mass(), triggerWeight);
          target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, 1.8, triggerWeight);
        } else {
          LOGF(fatal, "Can not fill invariant-mass user axis without invMass column. Disable cfgUserAxis or select another mode.");
        }
      } else {
        target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, triggerWeight);
      }

      const bool triggerHasTemplate = seedEstimate && hasTriggerTemplate(multiplicity, track1.pt());
      if (triggerHasTemplate) {
        ++seedEstimate->nTriggers;
        if (cfgEventSeedEstimatorMethod == 1) {
          addTriggerPriorExpectations(*seedEstimate, multiplicity, track1.pt());
        }
      }

      for (const auto& track2 : tracks2) {
        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if (track1.globalIndex() == track2.globalIndex()) {
            // LOGF(info, "Track identical: %f | %f | %f || %f | %f | %f", track1.eta(), track1.phi(), track1.pt(),  track2.eta(), track2.phi(), track2.pt());
            continue;
          }
        }

        if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks2::iterator>::value) { // skip those that are specifically chosen to be triggers
          if (!cfgMcTriggerPDGs->empty() && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track2.pdgCode()) != cfgMcTriggerPDGs->end()) {
            continue;
          }
        }

        // Daughter particle checks
        if constexpr (std::experimental::is_detected<HasPartDaugh0Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfParticleDaugh0Id()) { // do not correlate daughter particles of the same event
            continue;
          }
        }
        if constexpr (std::experimental::is_detected<HasPartDaugh1Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfParticleDaugh1Id()) { // do not correlate daughter particles of the same event
            continue;
          }
        }

        if constexpr (step <= CorrelationContainer::kCFStepTracked && !std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (!checkObject<step>(track2)) {
            continue;
          }
        }

        if constexpr (std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track2.decay()))) == 0u) {
            continue; // skip particles that do not match the decay mask
          }
        }

        if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value && std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (cfgCorrelationMethod == CorrelationMethod::Dd && track1.decay() != track2.decay()) {
            continue;
          }
          if (cfgCorrelationMethod == CorrelationMethod::Ddbar && track1.decay() == track2.decay()) {
            continue;
          }
        }

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        if constexpr (std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value || std::experimental::is_detected<HasPDGCode, typename TTracks2::iterator>::value) {
          const int associatedSign = getParticleSign(track2);
          if (cfgAssociatedCharge != 0) {
            if (cfgAssociatedCharge * associatedSign < 0) {
              continue;
            }
          } else if (associatedSign == 0) { // mc particles come in neutrals, need to check explicitly
            continue;
          }
        }

        if constexpr ((std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value || std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) &&
                      (std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value || std::experimental::is_detected<HasPDGCode, typename TTracks2::iterator>::value)) {
          if (cfgPairCharge != 0 && cfgPairCharge * getParticleSign(track1) * getParticleSign(track2) < 0) {
            continue;
          }
        }

        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if constexpr (step >= CorrelationContainer::kCFStepReconstructed) {
            if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value && std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
              if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
                continue;
              }
              if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, magField)) {
                continue;
              }
            }
          }
        }

        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= efficiencyAssociatedCache[track2.filteredIndex()];
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PIHalf);
        const float deltaEta = track1.eta() - track2.eta();

        if (triggerHasTemplate) {
          ++seedEstimate->nCandidatePairs;
          if (const auto* yieldTemplate = findYieldTemplate(multiplicity, track1.pt(), track2.pt())) {
            addPairProbabilities(*seedEstimate, *yieldTemplate, multiplicity, deltaPhi, deltaEta, posZ, fillEstimatorAcceptanceQA);
          } else {
            ++seedEstimate->nPairsWithoutTemplate;
          }
        }

        // last param is the weight
        if (cfgUserAxis == EventSeedAxis) {
          if (pendingPairFills) {
            pendingPairFills->push_back({deltaEta, track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, associatedWeight});
          } else {
            target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, eventSeed, associatedWeight);
          }
        } else if (cfgUserAxis == InvariantMassAxis) {
          if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value) {
            target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, track1.invMass(), associatedWeight);
          } else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
            target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, 1.8, associatedWeight); // p->Mass()
          } else {
            LOGF(fatal, "Can not fill invariant-mass user axis without invMass column. Disable cfgUserAxis or select another mode.");
          }
        } else {
          target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, associatedWeight);
        }
      }
    }
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (!cfgEfficiencyTrigger.value.empty()) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyTrigger.value.c_str(), "READ");
        cfg.mEfficiencyTrigger = dynamic_cast<THn*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      }
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
    }
    if (!cfgEfficiencyAssociated.value.empty()) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyAssociated = TFile::Open(cfgEfficiencyAssociated.value.c_str(), "READ");
        cfg.mEfficiencyAssociated = dynamic_cast<THn*>(fEfficiencyAssociated->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      }
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
    }
    cfg.efficiencyLoaded = true;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    std::array<int, 4> effVars{};
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars.data());
  }

  float getCorrelationContainerMultiplicity(float multiplicity) const
  {
    return multiplicity;
  }

  template <typename TCollision, typename TTracks>
  void processSameAODT(TCollision const& collision, TTracks const& tracks, const int* trueNMPI = nullptr)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    if (cfgVerbosity > 0) {
      LOGF(info, "processSameAOD: Tracks for collision: %d | Vertex: %.1f | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.sel7(), collision.centRun2V0M());
    }

    // TODO will go to CCDBConfigurable
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    loadEfficiency(bc.timestamp());
    loadCcdbYieldTemplates(bc.timestamp());

    const auto multiplicity = collision.centRun2V0M();

    if (!fillCollisionAOD(same, collision, multiplicity)) {
      return;
    }
    registry.fill(HIST("eventcount_same"), -2);
    fillQA(collision, multiplicity, tracks);
    EventSeedEstimate seedEstimate;
    std::vector<PendingSeedTriggerFill> pendingTriggerFills;
    std::vector<PendingSeedPairFill> pendingPairFills;
    if (cfgUserAxis == EventSeedAxis) {
      pendingTriggerFills.reserve(tracks.size());
      pendingPairFills.reserve(tracks.size() * tracks.size());
    }
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, collision.posZ(), getMagneticField(bc.timestamp()), 1.0f, &seedEstimate,
                                                                 cfgUserAxis == EventSeedAxis ? &pendingTriggerFills : nullptr,
                                                                 cfgUserAxis == EventSeedAxis ? &pendingPairFills : nullptr);
    finalizeEventSeedEstimate(seedEstimate);
    if (cfgUserAxis == EventSeedAxis) {
      flushSeedAxisFills<CorrelationContainer::kCFStepReconstructed>(same, pendingTriggerFills, pendingPairFills, seedEstimate.nuncSeeds());
    }
    fillEventSeedEstimatorQA(multiplicity, seedEstimate);
    if (trueNMPI) {
      fillMCValidation(multiplicity, seedEstimate, *trueNMPI);
    }
  }

  // Version with explicit nested loop
  void processSameAOD(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    processSameAODT(collision, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processSameAOD, "Process same event on AOD", true);

  void processSameGenMC(soa::Filtered<aod::CFMcCollisionsWithExtra>::iterator const& mcCollision,
                        soa::Filtered<aod::CFMcParticles> const& mcParticles,
                        soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    if (!cfgNuncSeedsTemplate.value.empty()) {
      for (const auto& collision : collisions) {
        loadCcdbYieldTemplates(collision.timestamp());
        break;
      }
    }

    const auto generatedMultiplicity = mcCollision.multiplicity();
    fillContainerEvent(same, generatedMultiplicity, CorrelationContainer::kCFStepAll);
    EventSeedEstimate seedEstimate;
    std::vector<PendingSeedTriggerFill> pendingTriggerFills;
    std::vector<PendingSeedPairFill> pendingPairFills;
    if (cfgUserAxis == EventSeedAxis) {
      pendingTriggerFills.reserve(mcParticles.size());
      pendingPairFills.reserve(mcParticles.size() * mcParticles.size());
    }
    fillCorrelations<CorrelationContainer::kCFStepAll>(same, mcParticles, mcParticles, generatedMultiplicity, mcCollision.posZ(), 0, 1.0f, &seedEstimate,
                                                       cfgUserAxis == EventSeedAxis ? &pendingTriggerFills : nullptr,
                                                       cfgUserAxis == EventSeedAxis ? &pendingPairFills : nullptr);
    finalizeEventSeedEstimate(seedEstimate);
    if (cfgUserAxis == EventSeedAxis) {
      flushSeedAxisFills<CorrelationContainer::kCFStepAll>(same, pendingTriggerFills, pendingPairFills, seedEstimate.nuncSeeds());
    }
    fillGeneratedMCValidation(generatedMultiplicity, seedEstimate, mcCollision.nMPI());
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processSameGenMC, "Process generated MC events from derived data and validate against the stored HepMC N MPI", false);

  template <class CollType, class TTracks1, class TTracks2>
  void processSameDerivedT(CollType const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2, const int* trueNMPI = nullptr)
  {
    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerivedT: Tracks for collision: %d/%d | Vertex: %.1f | Multiplicity/Centrality: %.1f", tracks1.size(), tracks2.size(), collision.posZ(), collision.multiplicity());
    }
    loadEfficiency(collision.timestamp());
    loadCcdbYieldTemplates(collision.timestamp());

    const auto multiplicity = collision.multiplicity();

    int field = 0;
    if (cfgTwoTrackCut > 0) {
      field = getMagneticField(collision.timestamp());
    }

    int bin = configurableBinningDerived.getBin({collision.posZ(), collision.multiplicity()});
    registry.fill(HIST("eventcount_same"), bin);
    registry.fill(HIST("trackcount_same"), bin, tracks1.size());
    if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
      fillQA(collision, multiplicity, collision.posZ(), tracks1, tracks2);
    } else {
      fillQA(collision, multiplicity, tracks1);
    }

    const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
    const bool fillReco = !(cfgDropStepRECO && hasEfficiency);
    EventSeedEstimate seedEstimate;

    if (cfgUserAxis == EventSeedAxis) {
      std::vector<PendingSeedTriggerFill> pendingTriggerFills;
      std::vector<PendingSeedPairFill> pendingPairFills;
      pendingTriggerFills.reserve(tracks1.size());
      pendingPairFills.reserve(tracks1.size() * tracks2.size());

      if (fillReco) {
        fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepReconstructed);
        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, &seedEstimate, &pendingTriggerFills, &pendingPairFills);
        finalizeEventSeedEstimate(seedEstimate);
        flushSeedAxisFills<CorrelationContainer::kCFStepReconstructed>(same, pendingTriggerFills, pendingPairFills, seedEstimate.nuncSeeds());
      } else if (hasEfficiency) {
        fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepCorrected);
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, &seedEstimate, &pendingTriggerFills, &pendingPairFills);
        finalizeEventSeedEstimate(seedEstimate);
        flushSeedAxisFills<CorrelationContainer::kCFStepCorrected>(same, pendingTriggerFills, pendingPairFills, seedEstimate.nuncSeeds());
      }

      if (fillReco && hasEfficiency) {
        fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepCorrected);
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, nullptr, nullptr, nullptr, getEventSeedUserAxisValue(multiplicity, seedEstimate.nuncSeeds()));
      }
      fillEventSeedEstimatorQA(multiplicity, seedEstimate);
      if (trueNMPI) {
        fillMCValidation(multiplicity, seedEstimate, *trueNMPI);
      }
      return;
    }

    if (fillReco) {
      fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, &seedEstimate);
    }
    if (hasEfficiency) {
      fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, fillReco ? nullptr : &seedEstimate);
    }
    finalizeEventSeedEstimate(seedEstimate);
    fillEventSeedEstimatorQA(multiplicity, seedEstimate);
    if (trueNMPI) {
      fillMCValidation(multiplicity, seedEstimate, *trueNMPI);
    }
  }

  void processSameDerived(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processSameDerived, "Process same event on derived data", false);

  void processSameDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>>::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    if (!passOutlier(collision)) {
      return;
    }
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processSameDerivedMultSet, "Process same event on derived data with multiplicity sets", false);

  using BinningTypeAOD = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentRun2V0M>;
  void processMixedAOD(AodCollisions const& collisions, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    // Strictly upper categorised collisions, for cfgNumMixedEvents combinations per bin, skipping those in entry -1
    BinningTypeAOD configurableBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<AodCollisions, AodTracks, BinningTypeAOD> pairs{configurableBinning, cfgNumMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    int skipID = -1;
    double triggerEventSeed = -1.0;
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      int bin = configurableBinning.getBin({collision1.posZ(), collision1.centRun2V0M()});
      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedAOD: Mixed collisions bin: %d pair: %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, collision1.globalIndex(), collision1.posZ(), collision1.centRun2V0M(), collision2.globalIndex(), collision2.posZ(), collision2.centRun2V0M());
      }
      if (collision1.globalIndex() == skipID) {
        continue;
      }

      if (it.isNewWindow()) {
        skipID = -1;
        if (!fillCollisionAOD(mixed, collision1, collision1.centRun2V0M())) {
          skipID = collision1.globalIndex();
          continue;
        }
        if (cfgUserAxis == EventSeedAxis) {
          auto bc = collision1.bc_as<aod::BCsWithTimestamps>();
          loadCcdbYieldTemplates(bc.timestamp());
          triggerEventSeed = getEventSeedUserAxisValue(collision1.centRun2V0M(), estimateEventSeedWithoutFilling<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp())));
        }
      }
      if (!collision2.alias_bit(kINT7) || !collision2.sel7()) {
        continue;
      }

      registry.fill(HIST("eventcount_mixed"), bin);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()), 1.0f / it.currentWindowNeighbours(), nullptr, nullptr, nullptr, triggerEventSeed);
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMixedAOD, "Process mixed events on AOD", false);

  template <class CollType, typename... TrackTypes>
  void processMixedDerivedT(CollType const& collisions, TrackTypes&&... tracks)
  {
    auto getMultiplicity =
      [this](auto& col) {
        if constexpr (std::experimental::is_detected<HasMultSet, CollType>::value) {
          if (!passOutlier(col)) {
            return -1.0f;
          }
        } else {
          (void)this; // fix compile error on unused 'this' capture
        }
        return col.multiplicity();
      };

    using BinningTypeDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    BinningTypeDerived configurableBinningDerived{{getMultiplicity}, {axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    //  Strictly upper categorised collisions, for cfgNumMixedEvents combinations per bin, skipping those in entry -1
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TA = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<CollType, TA, TB, BinningTypeDerived> pairs{configurableBinningDerived, cfgNumMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    double triggerEventSeed = -1.0;
    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      float multiplicity = getMultiplicity(collision1);
      int bin = configurableBinningDerived.getBin(std::tuple(collision1.posZ(), multiplicity));
      float eventWeight = 1.0f / it.currentWindowNeighbours();
      int field = 0;
      if (cfgTwoTrackCut > 0) {
        field = getMagneticField(collision1.timestamp());
      }

      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(), collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      bool hasEfficiencyMixed = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
      bool fillRecoMixed = !(cfgDropStepRECO && hasEfficiencyMixed);

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());
        hasEfficiencyMixed = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
        fillRecoMixed = !(cfgDropStepRECO && hasEfficiencyMixed);

        if (cfgUserAxis == EventSeedAxis) {
          loadCcdbYieldTemplates(collision1.timestamp());
          if constexpr (std::is_same_v<std::remove_cvref_t<TA>, std::remove_cvref_t<TB>>) {
            triggerEventSeed = getEventSeedUserAxisValue(collision1.multiplicity(), estimateEventSeedWithoutFilling<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, collision1.multiplicity(), collision1.posZ(), field));
          } else {
            LOGF(fatal, "Event-seed user axis for mixed events requires the same trigger and associated track table so the trigger event can be estimated independently");
          }
        }

        if (fillRecoMixed) {
          fillContainerEvent(mixed, collision1.multiplicity(), CorrelationContainer::kCFStepReconstructed);
        }
      }

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());
      if (fillRecoMixed) {
        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight, nullptr, nullptr, nullptr, triggerEventSeed);
      }

      if (hasEfficiencyMixed) {
        if (it.isNewWindow()) {
          fillContainerEvent(mixed, collision1.multiplicity(), CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight, nullptr, nullptr, nullptr, triggerEventSeed);
      }
    }
  }

  void processMixedDerived(DerivedCollisions const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMixedDerived, "Process mixed events on derived data", false);

  void processMixedDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>> const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMixedDerivedMultSet, "Process mixed events on derived data with multiplicity sets", false);

  int getSpecies(int pdgCode)
  {
    switch (pdgCode) {
      case 211: // pion
      case -211:
        return 0;
      case 321: // Kaon
      case -321:
        return 1;
      case 2212: // proton
      case -2212:
        return 2;
      default:
        break;
    }
    if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), pdgCode) != cfgMcTriggerPDGs->end()) {
      return 4;
    }
    // The efficiency histogram is hardcoded to contain 5 species. Anything special will have the 4th slot.
    return 3;
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracksWithLabel> perCollision = aod::cftrack::cfCollisionId;
  void processMCEfficiency(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, aod::CFMcParticles const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const& tracks)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    }

    auto multiplicity = mcCollision.multiplicity();
    if (cfgCentBinsForMC > 0) {
      if (collisions.size() == 0) {
        return;
      }
      for (const auto& collision : collisions) {
        multiplicity = collision.multiplicity();
      }
    }
    // Primaries
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() && mcParticle.sign() != 0 && !(std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), mcParticle.pdgCode()) != cfgMcTriggerPDGs->end())) {
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
      }
    }
    for (const auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());
        LOGF(info, "  which has %d tracks", groupedTracks.size());
      }

      for (const auto& track : groupedTracks) {
        if (cfgTrackBitMask > 0 && (track.trackType() & (uint8_t)cfgTrackBitMask) != (uint8_t)cfgTrackBitMask) {
          continue;
        }
        if (track.has_cfMCParticle()) {
          const auto& mcParticle = track.cfMCParticle();
          if (mcParticle.isPhysicalPrimary()) {
            same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          }
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, mcParticle.eta(), mcParticle.pt(), getSpecies(mcParticle.pdgCode()), multiplicity, mcCollision.posZ());
          // LOGF(info, "Filled track %d", track.globalIndex());
        } else {
          // fake track
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ());
        }
      }
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMCEfficiency, "MC: Extract efficiencies", false);

  template <class McCollision, class Particles1, class Particles2>
  void processMCSameDerivedT(McCollision const& mcCollision, Particles1 const& mcParticles1, Particles2 const& mcParticles2, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    if (cfgVerbosity > 0) {
      LOGF(info, "processMCSameDerivedT. MC collision: %d, particles1: %d, particles2: %d, collisions: %d", mcCollision.globalIndex(), mcParticles1.size(), mcParticles2.size(), collisions.size());
    }

    auto multiplicity = mcCollision.multiplicity();
    if (cfgCentBinsForMC > 0) {
      if (collisions.size() == 0) {
        return;
      }
      for (const auto& collision : collisions) {
        multiplicity = collision.multiplicity();
      }
    }

    if (!(doprocessMCSameDerived || doprocessSameDerived || doprocessSameDerivedMultSet)) {
      if constexpr (std::experimental::is_detected<HasDecay, typename Particles1::iterator>::value) {
        fillQA(mcCollision, multiplicity, mcCollision.posZ(), mcParticles1, mcParticles2);
      } else {
        fillQA(mcCollision, multiplicity, mcParticles1);
      }
    }

    fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    if (collisions.size() == 0) {
      return;
    }

    fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
    fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelations<CorrelationContainer::kCFStepTracked>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    // kCFStepReconstructed and kCFStepCorrected are filled below for every
    // reconstructed collision associated with this MC collision.
  }

  // NOTE SmallGroups includes soa::Filtered always
  Preslice<aod::CFTracks> derivedTracksPerCollision = aod::cftrack::cfCollisionId;
  void processMCSameDerived(soa::Filtered<aod::CFMcCollisionsWithExtra>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, soa::Filtered<aod::CFTracks> const& tracks) // TODO. For mixed no need to check the daughters since the events are different
  {
    processMCSameDerivedT(mcCollision, mcParticles, mcParticles, collisions);
    if (eventSeedEstimatorEnabled && collisions.size() == 0) {
      registry.fill(HIST("mcValidation/status"), 0.0);
    }
    const int trueNMPI = mcCollision.nMPI();
    for (const auto& collision : collisions) {
      auto collisionTracks = tracks.sliceBy(derivedTracksPerCollision, collision.globalIndex());
      processSameDerivedT(collision, collisionTracks, collisionTracks, &trueNMPI);
    }
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMCSameDerived, "Process generated and reconstructed MC same events on derived data and validate against HepMC N MPI", false);

  PresliceUnsorted<aod::CFCollisionsWithLabel> collisionPerMCCollision = aod::cfcollision::cfMcCollisionId;
  template <typename... ParticleTypes>
  void processMCMixedDerivedT(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions, ParticleTypes&&... particles)
  {
    bool useMCMultiplicity = (cfgCentBinsForMC == 0);
    auto getMultiplicity =
      [&collisions, &useMCMultiplicity, this](auto& col) {
        if (useMCMultiplicity) {
          return col.multiplicity();
        }
        auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, col.globalIndex());
        if (groupedCollisions.size() == 0) {
          return -1.0f;
        }
        return groupedCollisions.begin().multiplicity();
      };

    using BinningTypeMCDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    BinningTypeMCDerived configurableBinning{{getMultiplicity}, {axisVertex, axisMultiplicity}, true};

    // Strictly upper categorised collisions, for cfgNumMixedEvents combinations per bin, skipping those in entry -1
    auto tuple = std::make_tuple(std::forward<ParticleTypes>(particles)...);
    using TA = std::tuple_element<0, decltype(tuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tuple)> - 1, decltype(tuple)>::type;
    Pair<soa::Filtered<aod::CFMcCollisions>, TA, TB, BinningTypeMCDerived> pairs{configurableBinning, cfgNumMixedEvents, -1, mcCollisions, tuple, &cache}; // -1 is the number of the bin to skip

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      float multiplicity = getMultiplicity(collision1);
      if (cfgVerbosity > 0) {
        int bin = configurableBinning.getBin(std::tuple(collision1.posZ(), multiplicity));
        LOGF(info, "processMCMixedDerived: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), getMultiplicity(collision1), collision2.globalIndex(), collision2.posZ(), getMultiplicity(collision2));
      }

      // STEP 0
      if (it.isNewWindow()) {
        fillContainerEvent(mixed, multiplicity, CorrelationContainer::kCFStepAll);
      }
      fillCorrelations<CorrelationContainer::kCFStepAll>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      // check if collision1 has at least one reconstructed collision
      auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, collision1.globalIndex());
      if (cfgVerbosity > 0) {
        LOGF(info, "Found %d related collisions", groupedCollisions.size());
      }
      if (groupedCollisions.size() == 0) {
        continue;
      }

      // STEP 2, 4, 5
      if (it.isNewWindow()) {
        fillContainerEvent(mixed, multiplicity, CorrelationContainer::kCFStepVertex);
        fillContainerEvent(mixed, multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
        fillContainerEvent(mixed, multiplicity, CorrelationContainer::kCFStepTracked);
      }
      fillCorrelations<CorrelationContainer::kCFStepVertex>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);
      fillCorrelations<CorrelationContainer::kCFStepTracked>(mixed, tracks1, tracks2, multiplicity, collision1.posZ(), 0, eventWeight);

      // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processMixedDerived
      //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
    }
  }

  void processMCMixedDerived(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions)
  {
    processMCMixedDerivedT(mcCollisions, collisions, mcParticles);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMCMixedDerived,
                 "Process MC mixed events on derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TwoParticleCorrelationsMpi>(cfgc)};
}
