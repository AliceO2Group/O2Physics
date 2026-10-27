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

/// \file mftcorrelations.cxx
/// \brief task for two-particle (hadron-hadron) angular correlations built from MFT tracks

/// \author Mintu Haldar <mintu.haldar@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DataFormatsParameters/GRPObject.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h> // needed for Pair<> (two different track tables mixed together)
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TFormula.h>
#include <THn.h>
#include <TMath.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <experimental/type_traits>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

// Bitmask describing which multiplicity/centrality estimators feed the "multCorrelations" QA
// histogram and, when set, the outlier-rejection formula.
namespace o2::aod
{
namespace mftmultset
{
DECLARE_SOA_COLUMN(Multiplicities, multiplicities, std::vector<float>); //! List of auxiliary multiplicities
enum MultiplicityEstimators : uint16_t {
  CentFT0C = 1 << 0,
  MultFV0A = 1 << 1,
  MultNTracksPV = 1 << 2,
  MultNTracksGlobal = 1 << 3,
  CentFT0M = 1 << 4
};
inline constexpr uint32_t NMultiplicityEstimators = __builtin_ctz(CentFT0M) + 1;
} // namespace mftmultset

DECLARE_SOA_TABLE(MFTMultSets, "AOD", "MFTMULTSET", mftmultset::Multiplicities); //! Auxilary multiplicity set table
using MFTMultSet = MFTMultSets::iterator;
} // namespace o2::aod

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct Mftcorrelations {
  SliceCache cache;

  // ------------------------------------------------------------------ Configurables
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")

  // Global acceptance
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrackEtaCut, float, 0.8f, "Eta cut for global (central-barrel) tracks used as MFT-Global trigger")
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrackPtMin, float, 0.15f, "Minimal pT for global tracks")
  O2_DEFINE_CONFIGURABLE(cfgGlobalTrackPtMax, float, 10.0f, "Maximal pT for global tracks")

  // MFT acceptance is roughly -3.6 < eta < -2.5; keep as configurables since it is tag/setup dependent.
  O2_DEFINE_CONFIGURABLE(cfgMFTCutEtaMin, float, -3.6f, "Minimal eta for MFT tracks")
  O2_DEFINE_CONFIGURABLE(cfgMFTCutEtaMax, float, -2.5f, "Maximal eta for MFT tracks")
  O2_DEFINE_CONFIGURABLE(cfgMFTCutPtMin, float, 0.5f, "Minimal pT for MFT tracks")
  O2_DEFINE_CONFIGURABLE(cfgMFTCutPtMax, float, 20.0f, "Maximal pT for MFT tracks")
  O2_DEFINE_CONFIGURABLE(cfgMFTMinNClusters, int, 5, "Minimal number of MFT clusters on track")

  O2_DEFINE_CONFIGURABLE(cfgPtOrder, int, 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)")
  O2_DEFINE_CONFIGURABLE(cfgTriggerCharge, int, 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative")
  O2_DEFINE_CONFIGURABLE(cfgAssociatedCharge, int, 0, "Select on charge of associated particle: 0 = all; 1 = positive; -1 = negative")
  O2_DEFINE_CONFIGURABLE(cfgPairCharge, int, 0, "Select on charge of the pair: 0 = all; 1 = like sign; -1 = unlike sign")

  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, int, 0, "0 = OFF and 1 = ON: read efficiency THn from a local ROOT file instead of CCDB")
  O2_DEFINE_CONFIGURABLE(cfgDropStepRECO, bool, false, "drop step RECO if efficiency correction is used")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyTrigger, std::string, "", "CCDB path (or local file if cfgLocalEfficiency=1) to efficiency object for trigger particles")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyAssociated, std::string, "", "CCDB path (or local file if cfgLocalEfficiency=1) to efficiency object for associated particles")
  // FIX #6d: efficiency for central-barrel (global) tracks, used in MFT-Global (trigger) and Global-Global channels.
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyGlobal, std::string, "", "CCDB path (or local file if cfgLocalEfficiency=1) to efficiency object for central-barrel (global) tracks")

  O2_DEFINE_CONFIGURABLE(cfgMultCorrelationsMask, uint16_t, 0, "Bitmask of multiplicity/centrality estimators to fill in multCorrelations QA / use for outlier rejection (see mftmultset enum)")
  O2_DEFINE_CONFIGURABLE(cfgMultCutFormula, std::string, "", "Multiplicity outlier-rejection formula. A result > 0 accepts the event. Parameters: [cFT0C], [mFV0A], [mPV], [mGlob], [cFT0M]")

  O2_DEFINE_CONFIGURABLE(cfgNoMixedEvents, int, 5, "Number of mixed events per event")
  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")

  // ------------------------------------------------------------------ Phi weight
  // NOTE: keep OFF (empty file) when the efficiency object already has a phi axis, otherwise phi is corrected twice.
  O2_DEFINE_CONFIGURABLE(cfgPhiWeightFile, std::string, "", "Local ROOT file with the phi-weight TH1 (empty = disabled)")
  O2_DEFINE_CONFIGURABLE(cfgPhiWeightName, std::string, "", "Name of the phi-weight TH1 object inside cfgPhiWeightFile")

  // ------------------------------------------------------------------ Axes
  ConfigurableAxis axisVertex{"axisVertex", {10, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2.2, 2.2}, "delta eta axis for histograms (MFT-MFT so somewhat narrower than TPC-TPC)"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -3.6, -2.5}, "eta axis for efficiency histograms (MFT acceptance)"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0}, "pt axis for efficiency histograms"};
  ConfigurableAxis axisPhiEfficiency{"axisPhiEfficiency", {36, 0, TwoPI}, "phi axis for efficiency histograms"};

  ConfigurableAxis axisEtaEfficiencyGlobal{"axisEtaEfficiencyGlobal", {16, -0.8, 0.8}, "eta axis for efficiency histograms (global tracks)"};
  ConfigurableAxis axisPtEfficiencyGlobal{"axisPtEfficiencyGlobal", {VARIABLE_WIDTH, 0.0, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 20.0}, "pt axis for efficiency histograms (global tracks)"};

  ConfigurableAxis axisMultCorrCent{"axisMultCorrCent", {100, 0, 100}, "multiplicity correlation axis for centralities"};
  ConfigurableAxis axisMultCorrV0{"axisMultCorrV0", {1000, 0, 100000}, "multiplicity correlation axis for V0 multiplicities"};
  ConfigurableAxis axisMultCorrMult{"axisMultCorrMult", {1000, 0, 1000}, "multiplicity correlation axis for track multiplicities"};

  // ------------------------------------------------------------------ QA Axes
  ConfigurableAxis axisMultQA{"axisMultQA", {100, 0, 100}, "multiplicity/centrality axis for QA histograms"};
  ConfigurableAxis axisPtQA{"axisPtQA", {VARIABLE_WIDTH, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 15.0, 20.0}, "pT axis for QA histograms (yields)"};
  ConfigurableAxis axisEtaQA{"axisEtaQA", {50, -4, -2}, "eta axis for QA histograms"};
  ConfigurableAxis axisPhiQA{"axisPhiQA", {200, 0, TwoPI}, "phi axis for QA histograms"};

  // Global MFT
  ConfigurableAxis axisEtaGlobalQA{"axisEtaGlobalQA", {40, -1, 1}, "eta axis for global-track QA histograms"};
  ConfigurableAxis axisDeltaEtaMftGlobal{"axisDeltaEtaMftGlobal", {60, 1.0, 5.0}, "delta eta axis for MFT-Global correlations (deltaEta = etaGlobal - etaMFT is POSITIVE, roughly 1.7-4.4); override via JSON if needed"};
  // FIX #4: delta-eta axis for the global-global (central-barrel hadron-hadron) channel.
  ConfigurableAxis axisDeltaEtaGlobal{"axisDeltaEtaGlobal", {32, -1.6, 1.6}, "delta eta axis for global-global (central-barrel hadron-hadron) correlations"};

  // ------------------------------------------------------------------ Event selection
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;

  // ------------------------------------------------------------------ MFT track selection
  Filter mftTrackFilter = (aod::fwdtrack::eta > cfgMFTCutEtaMin) && (aod::fwdtrack::eta < cfgMFTCutEtaMax) &&
                          (aod::fwdtrack::pt > cfgMFTCutPtMin) && (aod::fwdtrack::pt < cfgMFTCutPtMax);

  // ------------------------------------------------------------------ Global track selection
  Filter globalTrackFilter = (nabs(aod::track::eta) < cfgGlobalTrackEtaCut) &&
                             (aod::track::pt > cfgGlobalTrackPtMin) && (aod::track::pt < cfgGlobalTrackPtMax) &&
                             requireGlobalTrackInFilter();

  // ------------------------------------------------------------------ Output
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  OutputObj<CorrelationContainer> sameMftGlobal{"sameEventMftGlobal"};
  OutputObj<CorrelationContainer> mixedMftGlobal{"mixedEventMftGlobal"};

  // FIX #4: global-global (central-barrel hadron-hadron) output objects.
  OutputObj<CorrelationContainer> sameGlobal{"sameEventGlobalGlobal"};
  OutputObj<CorrelationContainer> mixedGlobal{"mixedEventGlobalGlobal"};

  HistogramRegistry registry{"registry"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  struct Config {
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    THn* mEfficiencyGlobal = nullptr; // FIX #6d
    TH1* mPhiWeight = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  std::vector<float> efficiencyAssociatedCache;
  std::vector<float> phiWeightCache;
  std::vector<float> phiWeightCache1;
  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint32_t, aod::mftmultset::NMultiplicityEstimators> multCutFormulaParamIndex;

  // Collisions are expected to already carry the joined centrality/multiplicity estimator tables.
  using ExtendedCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::Mults>;
  using FilteredCollisions = soa::Filtered<ExtendedCollisions>;
  using FilteredMFTTracks = soa::Filtered<aod::MFTTracks>;
  using FilteredGlobalTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection>>;

  // For MFT MC Efficiency
  // NOTE SmallGroups includes soa::Filtered always
  using McEffCollisions = soa::SmallGroups<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::Mults>>;
  using McEffMFTTracks = soa::Filtered<soa::Join<aod::MFTTracks, aod::McMFTTrackLabels>>;

  // For Global MC Efficiency
  using McEffGlobalTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::McTrackLabels>>;

  void init(o2::framework::InitContext&)
  {
    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{axisMultQA, "multiplicity/centrality"}, {axisPtQA, "p_{T}"}, {axisEtaQA, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{axisMultQA, "multiplicity/centrality"}, {axisEtaQA, "#eta"}, {axisPhiQA, "#varphi"}}});
    // FIX #6b: efficiency-corrected QA histograms (weight = efficiency correction factor)
    registry.add("yields_corrected", "efficiency-corrected multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{axisMultQA, "multiplicity/centrality"}, {axisPtQA, "p_{T}"}, {axisEtaQA, "#eta"}}});
    registry.add("etaphi_corrected", "efficiency-corrected multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{axisMultQA, "multiplicity/centrality"}, {axisEtaQA, "#eta"}, {axisPhiQA, "#varphi"}}});
    registry.add("multiplicity", "event multiplicity", {HistType::kTH1F, {{axisMultQA, "multiplicity/centrality"}}});
    registry.add("nClusters", "MFT track number of clusters", {HistType::kTH1F, {{12, -0.5, 11.5, "n clusters"}}});
    registry.add("multiplicity_raw", "event multiplicity (no cuts)", {HistType::kTH1F, {{axisMultQA, "multiplicity/centrality"}}});

    // MC Efficiency
    registry.add("MCEff/hNRecoColl", "N reco collisions per MC collision", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    registry.add("MCEff/hCutFlow", "1=all MC coll, 2=size==1, 3=sel8", {HistType::kTH1F, {{3, 0.5, 3.5}}});

    registry.add("MCEffGlobal/hNRecoColl", "N reco collisions per MC collision (global)", {HistType::kTH1F, {{6, -0.5, 5.5}}});
    registry.add("MCEffGlobal/hCutFlow", "1=all MC coll, 2=size==1, 3=sel8", {HistType::kTH1F, {{3, 0.5, 3.5}}});

    if (cfgMultCorrelationsMask != 0) {
      std::vector<AxisSpec> multAxes;
      if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0C)
        multAxes.emplace_back(axisMultCorrCent, "FT0C centrality");
      if (cfgMultCorrelationsMask & aod::mftmultset::MultFV0A)
        multAxes.emplace_back(axisMultCorrV0, "V0A multiplicity");
      if (cfgMultCorrelationsMask & aod::mftmultset::MultNTracksPV)
        multAxes.emplace_back(axisMultCorrMult, "Nch PV");
      if (cfgMultCorrelationsMask & aod::mftmultset::MultNTracksGlobal)
        multAxes.emplace_back(axisMultCorrMult, "Nch Global");
      if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0M)
        multAxes.emplace_back(axisMultCorrCent, "FT0M centrality");
      registry.add("multCorrelations", "Multiplicity correlations", {HistType::kTHnSparseF, multAxes});

      // FIX #3
      if (!(cfgMultCorrelationsMask & aod::mftmultset::CentFT0C)) {
        LOGF(warning,
             "cfgMultCorrelationsMask=%d does not include CentFT0C. Mixed-event binning follows "
             "getMultiplicity(), so make sure the efficiency objects were produced with the same estimator.",
             (int)cfgMultCorrelationsMask);
      }
    }

    const int maxMixBin = AxisSpec(axisMultiplicity).getNbins() * AxisSpec(axisVertex).getNbins();
    registry.add("eventcount_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_same", "bin", {HistType::kTH2F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}}});
    registry.add("trackcount_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}, {50, -0.5, 49.5}}});

    // FIX #1
    registry.add("eventcount_mftglobal_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mftglobal_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_mftglobal_same", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}, {50, -0.5, 49.5}}});
    registry.add("trackcount_mftglobal_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}, {50, -0.5, 49.5}}});

    if (!cfgMultCutFormula.value.empty()) {
      multCutFormula = std::make_unique<TFormula>("multCutFormula", cfgMultCutFormula.value.c_str());
      std::fill_n(multCutFormulaParamIndex.begin(), std::size(multCutFormulaParamIndex), ~0u);
      std::array<std::string, aod::mftmultset::NMultiplicityEstimators> pars = {"cFT0C", "mFV0A", "mPV", "mGlob", "cFT0M"};
      for (uint32_t i = 0, n = multCutFormula->GetNpar(); i < n; ++i) {
        auto m = std::find(pars.begin(), pars.end(), multCutFormula->GetParName(i));
        if (m == pars.end()) {
          LOGF(warning, "Unknown parameter in cfgMultCutFormula: %s", multCutFormula->GetParName(i));
          continue;
        }
        multCutFormulaParamIndex[std::distance(pars.begin(), m)] = i;
      }
    }

    registry.add("mftmft/SE/deltaEta_deltaPhi", "same event #Delta#eta-#Delta#varphi (MFT-MFT)", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEta, "#Delta#eta"}}});
    registry.add("mftmft/ME/deltaEta_deltaPhi", "mixed event #Delta#eta-#Delta#varphi (MFT-MFT)", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEta, "#Delta#eta"}}});

    registry.add("mftGlobal/SE/Trig_etavsphi", "trigger (global) eta vs phi", {HistType::kTH2F, {{axisPhiQA, "#varphi"}, {axisEtaGlobalQA, "#eta"}}});
    registry.add("mftGlobal/SE/Trig_pt", "trigger (global) pT", {HistType::kTH1F, {{axisPtTrigger, "p_{T}"}}});
    registry.add("mftGlobal/SE/Assoc_etavsphi", "associated (MFT) eta vs phi", {HistType::kTH2F, {{axisPhiQA, "#varphi"}, {axisEtaQA, "#eta"}}});
    registry.add("mftGlobal/SE/Assoc_pt", "associated (MFT) pT", {HistType::kTH1F, {{axisPtAssoc, "p_{T}"}}});
    registry.add("mftGlobal/SE/deltaEta_deltaPhi", "same event #Delta#eta-#Delta#varphi", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEtaMftGlobal, "#Delta#eta"}}});
    registry.add("mftGlobal/ME/deltaEta_deltaPhi", "mixed event #Delta#eta-#Delta#varphi", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEtaMftGlobal, "#Delta#eta"}}});

    // FIX #4: QA histograms for the global-global channel.
    registry.add("eventcount_global_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_global_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_global_same", "bin", {HistType::kTH2F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}}});
    registry.add("trackcount_global_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {50, -0.5, 49.5}, {50, -0.5, 49.5}}});

    registry.add("global/SE/etavsphi", "global track eta vs phi", {HistType::kTH2F, {{axisPhiQA, "#varphi"}, {axisEtaGlobalQA, "#eta"}}});
    registry.add("global/SE/pt", "global track pT", {HistType::kTH1F, {{axisPtQA, "p_{T}"}}});
    registry.add("global/SE/deltaEta_deltaPhi", "same event #Delta#eta-#Delta#varphi (global-global)", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEtaGlobal, "#Delta#eta"}}});
    registry.add("global/ME/deltaEta_deltaPhi", "mixed event #Delta#eta-#Delta#varphi (global-global)", {HistType::kTH2F, {{axisDeltaPhi, "#Delta#varphi"}, {axisDeltaEtaGlobal, "#Delta#eta"}}});

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"},
                                     {axisPhiEfficiency, "#varphi (rad)"}};

    std::vector<AxisSpec> effAxisGlobal = {{axisEtaEfficiencyGlobal, "#eta"},
                                           {axisPtEfficiencyGlobal, "p_{T} (GeV/c)"},
                                           {axisVertexEfficiency, "z-vtx (cm)"},
                                           {axisPhiEfficiency, "#varphi (rad)"}};

    std::vector<AxisSpec> userAxis; // no extra axis needed for hadron-hadron

    std::vector<AxisSpec> corrAxisMftGlobal = {{axisDeltaEtaMftGlobal, "#Delta#eta"},
                                               {axisPtAssoc, "p_{T} (GeV/c)"},
                                               {axisPtTrigger, "p_{T} (GeV/c)"},
                                               {axisMultiplicity, "multiplicity / centrality"},
                                               {axisDeltaPhi, "#Delta#varphi (rad)"},
                                               {axisVertex, "z-vtx (cm)"}};

    // FIX #4
    std::vector<AxisSpec> corrAxisGlobal = {{axisDeltaEtaGlobal, "#Delta#eta"},
                                            {axisPtAssoc, "p_{T} (GeV/c)"},
                                            {axisPtTrigger, "p_{T} (GeV/c)"},
                                            {axisMultiplicity, "multiplicity / centrality"},
                                            {axisDeltaPhi, "#Delta#varphi (rad)"},
                                            {axisVertex, "z-vtx (cm)"}};

    sameMftGlobal.setObject(new CorrelationContainer("sameEventMftGlobal", "sameEventMftGlobal", corrAxisMftGlobal, effAxis, userAxis));
    mixedMftGlobal.setObject(new CorrelationContainer("mixedEventMftGlobal", "mixedEventMftGlobal", corrAxisMftGlobal, effAxis, userAxis));

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    same->setTrackEtaCut(cfgMFTCutEtaMax); // upper bound of |eta| cut used internally by CorrelationContainer

    sameGlobal.setObject(new CorrelationContainer("sameEventGlobalGlobal", "sameEventGlobalGlobal", corrAxisGlobal, effAxisGlobal, userAxis));
    mixedGlobal.setObject(new CorrelationContainer("mixedEventGlobalGlobal", "mixedEventGlobalGlobal", corrAxisGlobal, effAxisGlobal, userAxis));
    sameGlobal->setTrackEtaCut(cfgGlobalTrackEtaCut);

    if (!cfgEfficiencyAssociated.value.empty())
      efficiencyAssociatedCache.reserve(512);

    // FIX #5: load the phi weight once at init().
    loadPhiWeight();

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  // ------------------------------------------------------------------ helpers

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
    }
    return grpo->getNominalL3Field();
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (!cfgEfficiencyTrigger.value.empty()) {
      if (cfgLocalEfficiency > 0) {
        TFile* f = TFile::Open(cfgEfficiencyTrigger.value.c_str(), "READ");
        cfg.mEfficiencyTrigger = reinterpret_cast<THn*>(f->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      }
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load trigger efficiency from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded MFT trigger efficiency from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (!cfgEfficiencyAssociated.value.empty()) {
      if (cfgLocalEfficiency > 0) {
        TFile* f = TFile::Open(cfgEfficiencyAssociated.value.c_str(), "READ");
        cfg.mEfficiencyAssociated = reinterpret_cast<THn*>(f->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      }
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load associated efficiency from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded MFT associated efficiency from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }
    // FIX #6d
    if (!cfgEfficiencyGlobal.value.empty()) {
      if (cfgLocalEfficiency > 0) {
        TFile* f = TFile::Open(cfgEfficiencyGlobal.value.c_str(), "READ");
        cfg.mEfficiencyGlobal = reinterpret_cast<THn*>(f->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyGlobal = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyGlobal, timestamp);
      }
      if (cfg.mEfficiencyGlobal == nullptr) {
        LOGF(fatal, "Could not load global-track efficiency from %s", cfgEfficiencyGlobal.value.c_str());
      }
      LOGF(info, "Loaded global-track efficiency from %s (%p)", cfgEfficiencyGlobal.value.c_str(), (void*)cfg.mEfficiencyGlobal);
    }
    cfg.efficiencyLoaded = true;
  }

  // MFT phi is in [-pi,pi], McParticle/global-track phi in [0,2pi). Use [0,2pi) everywhere.
  static float wrapPhi(float phi)
  {
    return RecoDecay::constrainAngle(phi, 0.f);
  }

  // FIX #6a: expected axis order of the efficiency object: (eta, pt, multiplicity, posZ [, phi]).
  //  - phi is wrapped to [0, 2pi) (MFT track.phi() is in [-pi, pi])
  //  - 4-dim and 5-dim objects are supported
  //  - bins outside the axis range (under/overflow) return 0 = track outside efficiency acceptance
  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ, float phi)
  {
    if (eff == nullptr) {
      return 1.0;
    }
    constexpr int NDimWithoutPhi = 4; // eta, pt, mult, z
    constexpr int NDimWithPhi = 5;    // eta, pt, mult, z, phi
    const int nDim = eff->GetNdimensions();
    if (nDim != NDimWithoutPhi && nDim != NDimWithPhi) {
      LOGF(fatal, "Efficiency object has %d dimensions, expected %d (eta,pt,mult,z) or %d (+phi)", nDim, NDimWithoutPhi, NDimWithPhi);
    }
    const double vals[NDimWithPhi] = {eta, pt, multiplicity, posZ, static_cast<double>(RecoDecay::constrainAngle(phi, 0.f))};
    int bins[NDimWithPhi];
    for (int i = 0; i < nDim; ++i) {
      TAxis* ax = eff->GetAxis(i);
      bins[i] = ax->FindBin(vals[i]);
      if (bins[i] < 1 || bins[i] > ax->GetNbins()) {
        return 0.0;
      }
    }
    return eff->GetBinContent(bins);
  }

  template <class T>
  using HasSign = decltype(std::declval<T&>().sign());
  template <class T>
  using HasPdgCode = decltype(std::declval<T&>().pdgCode());
  template <class T>
  using HasIsPhysicalPrimary = decltype(std::declval<T&>().isPhysicalPrimary());

  // Retrieve the multiplicity/centrality estimator selected by cfgMultCorrelationsMask.
  // Falls back to FT0C centrality if no mask bit is set.
  template <typename TCollision>
  float getMultiplicity(TCollision const& collision)
  {
    if (cfgMultCorrelationsMask == 0) {
      return collision.centFT0C();
    }
    if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0C)
      return collision.centFT0C();
    if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0M)
      return collision.centFT0M();
    return collision.centFT0C();
  }

  template <typename TCollision>
  bool passOutlier(TCollision const& collision)
  {
    if (cfgMultCutFormula.value.empty() || multCutFormula == nullptr)
      return true;
    if (multCutFormulaParamIndex[0] != ~0u)
      multCutFormula->SetParameter(multCutFormulaParamIndex[0], collision.centFT0C());
    if (multCutFormulaParamIndex[1] != ~0u)
      multCutFormula->SetParameter(multCutFormulaParamIndex[1], collision.multFV0A());
    if (multCutFormulaParamIndex[2] != ~0u)
      multCutFormula->SetParameter(multCutFormulaParamIndex[2], collision.multNTracksPV());
    if (multCutFormulaParamIndex[3] != ~0u)
      multCutFormula->SetParameter(multCutFormulaParamIndex[3], collision.multNTracksPV()); // substitute a "global" estimator if/when available
    if (multCutFormulaParamIndex[4] != ~0u)
      multCutFormula->SetParameter(multCutFormulaParamIndex[4], collision.centFT0M());
    return multCutFormula->Eval() > 0.0f;
  }

  template <typename TCollision>
  void fillMultCorrelations(TCollision const& collision)
  {
    if (cfgMultCorrelationsMask == 0)
      return;
    std::vector<double> v;
    if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0C)
      v.push_back(collision.centFT0C());
    if (cfgMultCorrelationsMask & aod::mftmultset::MultFV0A)
      v.push_back(collision.multFV0A());
    if (cfgMultCorrelationsMask & aod::mftmultset::MultNTracksPV)
      v.push_back(collision.multNTracksPV());
    if (cfgMultCorrelationsMask & aod::mftmultset::MultNTracksGlobal)
      v.push_back(collision.multNTracksPV()); // substitute a "global" estimator if/when available
    if (cfgMultCorrelationsMask & aod::mftmultset::CentFT0M)
      v.push_back(collision.centFT0M());
    registry.get<THnSparse>(HIST("multCorrelations")).get()->Fill(v.data());
  }

  void loadPhiWeight()
  {
    if (cfg.mPhiWeight != nullptr || cfgPhiWeightFile.value.empty()) {
      return; // already loaded, or feature disabled
    }
    TFile* f = TFile::Open(cfgPhiWeightFile.value.c_str(), "READ");
    if (!f || f->IsZombie()) {
      LOGF(fatal, "Could not open phi weight file %s", cfgPhiWeightFile.value.c_str());
      return;
    }
    cfg.mPhiWeight = reinterpret_cast<TH1*>(f->Get(cfgPhiWeightName.value.c_str()));
    if (cfg.mPhiWeight) {
      cfg.mPhiWeight->SetDirectory(nullptr); // detach so it survives f->Close()
    } else {
      LOGF(fatal, "Could not find object %s in %s", cfgPhiWeightName.value.c_str(), cfgPhiWeightFile.value.c_str());
    }
    f->Close();
    LOGF(info, "Loaded phi weight '%s' from %s (%p)", cfgPhiWeightName.value.c_str(), cfgPhiWeightFile.value.c_str(), (void*)cfg.mPhiWeight);
  }

  // Looks up a bin-content-based weight from a 1-D phi map.
  // Returns 0 if the histogram is missing so the >0 check upstream skips it cleanly.
  float getPhiWeight(TH1* hPhiWeight, float phi)
  {
    if (!hPhiWeight) {
      return 0.f;
    }
    int binX = hPhiWeight->GetXaxis()->FindBin(wrapPhi(phi));
    return hPhiWeight->GetBinContent(binX);
  }

  // ------------------------------------------------------------------ QA
  // FIX #6b: "yields"/"etaphi" are RAW. "yields_corrected"/"etaphi_corrected" are weighted with the
  // efficiency correction factor (associated-particle efficiency object).
  template <typename TCollision, typename TTracks>
  void fillQA(TCollision const& collision, float multiplicity, TTracks const& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    fillMultCorrelations(collision);
    const float posZ = collision.posZ();

    for (const auto& track : tracks) {
      if (track.nClusters() < cfgMFTMinNClusters) {
        continue;
      }
      registry.fill(HIST("yields"), multiplicity, track.pt(), track.eta());

      float phiWeight = 1.f;
      if (cfg.mPhiWeight) {
        float phiWgt = getPhiWeight(cfg.mPhiWeight, track.phi());
        phiWeight = (phiWgt > 0) ? 1.f / phiWgt : 1.f;
      }
      registry.fill(HIST("etaphi"), multiplicity, track.eta(), wrapPhi(track.phi()), phiWeight);
      registry.fill(HIST("nClusters"), track.nClusters());

      float effWeight = 1.f;
      if (cfg.mEfficiencyAssociated) {
        effWeight = getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ, track.phi());
      }
      registry.fill(HIST("yields_corrected"), multiplicity, track.pt(), track.eta(), effWeight);
      registry.fill(HIST("etaphi_corrected"), multiplicity, track.eta(), wrapPhi(track.phi()), effWeight);
    }
  }

  // ------------------------------------------------------------------ event-level acceptance
  template <typename TCollision>
  bool acceptCollision(TCollision const& collision)
  {
    if (!collision.sel8()) {
      return false;
    }
    if (!passOutlier(collision)) {
      return false;
    }
    return true;
  }

  // ------------------------------------------------------------------ correlation filling
  // Hadron-hadron (MFT-MFT). tracks1 = trigger collection, tracks2 = associated collection.
  // Same-event: called as fillCorrelations(..., tracks, tracks, ...) -> self-correlation with index check.
  // Mixed-event: called as fillCorrelations(..., tracks1, tracks2, ...) from two different collisions.

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelations(TTarget target, TTracks1 const& tracks1, TTracks2 const& tracks2, float multiplicity, bool mixing, float posZ, float eventWeight)
  {
    // ---------------------------------------------------------------- per-call caches
    // Indexed by plain position in the loop (NOT filteredIndex()), and built for BOTH
    // collections. This is what makes the mixed event correct: tracks1 and tracks2 have the same
    // C++ type in SameKindPair, so a type-based "sameCollection" is true in ME as well.
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      efficiencyAssociatedCache.clear();
      if (cfg.mEfficiencyAssociated) {
        efficiencyAssociatedCache.reserve(tracks2.size());
        for (const auto& track : tracks2) {
          efficiencyAssociatedCache.push_back(getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ, track.phi()));
        }
      }
    }

    phiWeightCache.clear();  // associated side (tracks2)
    phiWeightCache1.clear(); // trigger side (tracks1)
    if (cfg.mPhiWeight) {
      phiWeightCache.reserve(tracks2.size());
      for (const auto& track : tracks2) {
        phiWeightCache.push_back(getPhiWeight(cfg.mPhiWeight, track.phi()));
      }
      phiWeightCache1.reserve(tracks1.size());
      for (const auto& track : tracks1) {
        phiWeightCache1.push_back(getPhiWeight(cfg.mPhiWeight, track.phi()));
      }
    }

    // ---------------------------------------------------------------- trigger loop
    int i1 = -1;
    for (const auto& track1 : tracks1) {
      ++i1; // must stay before any 'continue'
      if (track1.nClusters() < cfgMFTMinNClusters) {
        continue;
      }

      if constexpr (std::experimental::is_detected<HasPdgCode, typename TTracks1::iterator>::value) {
        if constexpr (std::experimental::is_detected<HasIsPhysicalPrimary, typename TTracks1::iterator>::value) {
          if (!track1.isPhysicalPrimary())
            continue;
        }
      } else if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value) {
        if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0)
          continue;
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ, track1.phi());
        }
      }
      // phi weight of the trigger: applied BEFORE the trigger histogram is filled, so the
      // per-trigger normalisation and the pair weights stay consistent.
      if (cfg.mPhiWeight) {
        const float w1 = phiWeightCache1[i1];
        if (w1 > 0.f) {
          triggerWeight *= 1.f / w1;
        }
      }
      target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, triggerWeight);

      // ---------------------------------------------------------------- associated loop
      int i2 = -1;
      for (const auto& track2 : tracks2) {
        ++i2; // must stay before any 'continue'

        if (!mixing && track1.globalIndex() == track2.globalIndex()) {
          continue; // no self-correlation (same event only)
        }
        if (track2.nClusters() < cfgMFTMinNClusters) {
          continue;
        }

        if constexpr (std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
          if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.sign() < 0)
            continue;
        }
        if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value &&
                      std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
          if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0)
            continue;
        }

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        float associatedWeight = triggerWeight; // already contains event, trigger eff, trigger phi weight
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= efficiencyAssociatedCache[i2];
          }
        }
        if (cfg.mPhiWeight) {
          const float w2 = phiWeightCache[i2];
          if (w2 > 0.f) {
            associatedWeight *= 1.f / w2;
          }
        }

        const float deltaEta = track1.eta() - track2.eta();
        const float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PIHalf);

        if (mixing) {
          registry.fill(HIST("mftmft/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        } else {
          registry.fill(HIST("mftmft/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        }

        target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
      }
    }
  }
  // FIX #2 + FIX #6d: MFT-Global. Trigger = central-barrel track (cfg.mEfficiencyGlobal),
  // associated = MFT track (cfg.mEfficiencyAssociated).
  template <CorrelationContainer::CFStep step, typename TTarget, typename TTriggers, typename TMFTs>
  void fillCorrMftGlobal(TTarget target, TTriggers const& triggers, TMFTs const& mftTracks, bool mixing, float posZ, float multiplicity, float eventWeight)
  {
    constexpr bool TriggerHasSign = std::experimental::is_detected<HasSign, typename TTriggers::iterator>::value;
    constexpr bool AssocHasSign = std::experimental::is_detected<HasSign, typename TMFTs::iterator>::value;

    if (!mixing) {
      for (const auto& assocTrack : mftTracks) {
        if (assocTrack.nClusters() < cfgMFTMinNClusters) {
          continue;
        }
        registry.fill(HIST("mftGlobal/SE/Assoc_etavsphi"), assocTrack.phi(), assocTrack.eta());
        registry.fill(HIST("mftGlobal/SE/Assoc_pt"), assocTrack.pt());
      }
    }

    for (const auto& triggerTrack : triggers) {
      if constexpr (TriggerHasSign) {
        if (cfgTriggerCharge != 0 && cfgTriggerCharge * triggerTrack.sign() < 0)
          continue;
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyGlobal) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyGlobal, triggerTrack.eta(), triggerTrack.pt(), multiplicity, posZ, triggerTrack.phi());
        }
      }
      if (!mixing) {
        registry.fill(HIST("mftGlobal/SE/Trig_etavsphi"), triggerTrack.phi(), triggerTrack.eta());
        registry.fill(HIST("mftGlobal/SE/Trig_pt"), triggerTrack.pt());
      }
      target->getTriggerHist()->Fill(step, triggerTrack.pt(), multiplicity, posZ, triggerWeight);

      for (const auto& assocTrack : mftTracks) {
        if (assocTrack.nClusters() < cfgMFTMinNClusters) {
          continue;
        }
        if constexpr (AssocHasSign) {
          if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * assocTrack.sign() < 0)
            continue;
        }
        if constexpr (TriggerHasSign && AssocHasSign) {
          if (cfgPairCharge != 0 && cfgPairCharge * triggerTrack.sign() * assocTrack.sign() < 0)
            continue;
        }
        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= getEfficiencyCorrection(cfg.mEfficiencyAssociated, assocTrack.eta(), assocTrack.pt(), multiplicity, posZ, assocTrack.phi());
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(triggerTrack.phi() - assocTrack.phi(), -o2::constants::math::PIHalf);
        float deltaEta = triggerTrack.eta() - assocTrack.eta();
        if (mixing) {
          registry.fill(HIST("mftGlobal/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        } else {
          registry.fill(HIST("mftGlobal/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        }
        target->getPairHist()->Fill(step, deltaEta, assocTrack.pt(), triggerTrack.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
      }
    }
  }

  // FIX #4 + FIX #6d: Global-Global (central-barrel hadron-hadron). Both tracks use cfg.mEfficiencyGlobal.
  // No nClusters() cut: global tracks are quality-selected via FilteredGlobalTracks.
  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrGlobalGlobal(TTarget target, TTracks1 const& tracks1, TTracks2 const& tracks2, bool mixing, float posZ, float multiplicity, float eventWeight)
  {
    constexpr bool Track1HasSign = std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value;
    constexpr bool Track2HasSign = std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value;

    for (const auto& track1 : tracks1) {
      if constexpr (Track1HasSign) {
        if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0)
          continue;
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyGlobal) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyGlobal, track1.eta(), track1.pt(), multiplicity, posZ, track1.phi());
        }
      }
      if (!mixing) {
        registry.fill(HIST("global/SE/etavsphi"), track1.phi(), track1.eta());
        registry.fill(HIST("global/SE/pt"), track1.pt());
      }
      target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, triggerWeight);

      for (const auto& track2 : tracks2) {
        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if (track1.globalIndex() == track2.globalIndex()) {
            continue; // no self-correlation
          }
        }

        if constexpr (Track2HasSign) {
          if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.sign() < 0)
            continue;
        }
        if constexpr (Track1HasSign && Track2HasSign) {
          if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0)
            continue;
        }
        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyGlobal) {
            associatedWeight *= getEfficiencyCorrection(cfg.mEfficiencyGlobal, track2.eta(), track2.pt(), multiplicity, posZ, track2.phi());
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PIHalf);
        float deltaEta = track1.eta() - track2.eta();
        if (mixing) {
          registry.fill(HIST("global/ME/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        } else {
          registry.fill(HIST("global/SE/deltaEta_deltaPhi"), deltaPhi, deltaEta);
        }
        target->getPairHist()->Fill(step, deltaEta, track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
      }
    }
  }

  // ------------------------------------------------------------------ SAME EVENT (data)
  void processSameMFT(FilteredCollisions::iterator const& collision, aod::BCsWithTimestamps const&, FilteredMFTTracks const& tracks)
  {
    registry.fill(HIST("multiplicity_raw"), collision.centFT0C()); // unconditional, no cuts
    if (!acceptCollision(collision)) {
      return;
    }

    loadEfficiency(collision.template bc_as<aod::BCsWithTimestamps>().timestamp());
    const float multiplicity = getMultiplicity(collision);

    using BinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    BinningType binning{{axisVertex, axisMultiplicity}, true};
    int bin = binning.getBin({collision.posZ(), multiplicity});

    if (cfgVerbosity > 0) {
      LOGF(info, "processSameMFT: %d MFT tracks | Vertex: %.2f | Mult: %.2f", tracks.size(), collision.posZ(), multiplicity);
    }

    registry.fill(HIST("eventcount_same"), bin);
    registry.fill(HIST("trackcount_same"), bin, tracks.size());
    fillQA(collision, multiplicity, tracks);

    const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
    const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    if (fillReco) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, false, collision.posZ(), 1.0f);
    }
    if (hasEfficiency) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks, tracks, multiplicity, false, collision.posZ(), 1.0f);
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processSameMFT, "Process same event MFT-MFT correlations", true);

  void processMftGlobalSE(FilteredCollisions::iterator const& collision, FilteredGlobalTracks const& globalTracks, FilteredMFTTracks const& mftTracks, aod::BCsWithTimestamps const&)
  {
    if (!acceptCollision(collision)) {
      return;
    }
    loadEfficiency(collision.template bc_as<aod::BCsWithTimestamps>().timestamp());
    const float multiplicity = getMultiplicity(collision);

    using MixingBinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    MixingBinningType binning{{axisVertex, axisMultiplicity}, true};
    int bin = binning.getBin({collision.posZ(), multiplicity});
    registry.fill(HIST("eventcount_mftglobal_same"), bin);
    registry.fill(HIST("trackcount_mftglobal_same"), bin, globalTracks.size(), mftTracks.size());

    // FIX #6d: trigger = global track, associated = MFT track
    const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyGlobal != nullptr);
    const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

    sameMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    if (fillReco) {
      sameMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(sameMftGlobal, globalTracks, mftTracks, false, collision.posZ(), multiplicity, 1.0f);
    }
    if (hasEfficiency) {
      sameMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrMftGlobal<CorrelationContainer::kCFStepCorrected>(sameMftGlobal, globalTracks, mftTracks, false, collision.posZ(), multiplicity, 1.0f);
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processMftGlobalSE, "Process same event MFT-vs-global-track correlations", false);

  // FIX #4: Same-event global-global (central-barrel hadron-hadron) correlations.
  void processGlobalGlobalSE(FilteredCollisions::iterator const& collision, FilteredGlobalTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    if (!acceptCollision(collision)) {
      return;
    }
    loadEfficiency(collision.template bc_as<aod::BCsWithTimestamps>().timestamp());
    const float multiplicity = getMultiplicity(collision);

    using MixingBinningType = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
    MixingBinningType binning{{axisVertex, axisMultiplicity}, true};
    int bin = binning.getBin({collision.posZ(), multiplicity});

    if (cfgVerbosity > 0) {
      LOGF(info, "processGlobalGlobalSE: %d global tracks | Vertex: %.2f | Mult: %.2f", tracks.size(), collision.posZ(), multiplicity);
    }

    registry.fill(HIST("eventcount_global_same"), bin);
    registry.fill(HIST("trackcount_global_same"), bin, tracks.size());

    // FIX #6d
    const bool hasEfficiency = (cfg.mEfficiencyGlobal != nullptr);
    const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

    sameGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    if (fillReco) {
      sameGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrGlobalGlobal<CorrelationContainer::kCFStepReconstructed>(sameGlobal, tracks, tracks, false, collision.posZ(), multiplicity, 1.0f);
    }
    if (hasEfficiency) {
      sameGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrGlobalGlobal<CorrelationContainer::kCFStepCorrected>(sameGlobal, tracks, tracks, false, collision.posZ(), multiplicity, 1.0f);
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processGlobalGlobalSE, "Process same event global-global (central-barrel) correlations", false);

  void processSameMFTMultSet(FilteredCollisions::iterator const& collision, aod::BCsWithTimestamps const& bcs, FilteredMFTTracks const& tracks)
  {
    if (cfgMultCorrelationsMask == 0) {
      LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when the MultSet process function is in use.");
    }
    processSameMFT(collision, bcs, tracks);
  }
  PROCESS_SWITCH(Mftcorrelations, processSameMFTMultSet, "Process same event MFT-MFT correlations with multi-estimator QA", false);

  // ------------------------------------------------------------------ MIXED EVENT (data)
  void processMixedMFT(FilteredCollisions const& collisions, FilteredMFTTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // Bin mixed events on the SAME multiplicity/centrality estimator that getMultiplicity() returns
    // (driven by cfgMultCorrelationsMask), instead of a hardcoded CentFT0C column. This keeps the
    // mixing pool synchronized with whatever estimator is actually filled into the CorrelationContainer,
    // and lets you switch estimators (e.g. FT0C -> FT0M) purely from the JSON config.
    auto getMult = [this](auto const& c) { return getMultiplicity(c); };
    using BinningType = FlexibleBinningPolicy<std::tuple<decltype(getMult)>, aod::collision::PosZ, decltype(getMult)>;
    BinningType binning{{getMult}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<FilteredCollisions, FilteredMFTTracks, BinningType> pairs{binning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!acceptCollision(collision1) || !collision2.sel8()) {
        continue;
      }

      const float multiplicity = getMultiplicity(collision1);
      int bin = binning.getBin(std::tuple(collision1.posZ(), multiplicity));
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      if (cfgVerbosity > 0) {
        LOGF(info, "processMixedMFT: bin %d | pair (%d,%d) | Vertex %.2f/%.2f | Mult %.2f", bin, collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), multiplicity);
      }

      if (it.isNewWindow()) {
        loadEfficiency(collision1.template bc_as<aod::BCsWithTimestamps>().timestamp());
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      }

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());

      const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
      const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

      if (it.isNewWindow() && fillReco) {
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      }
      if (it.isNewWindow() && hasEfficiency) {
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      }

      if (fillReco) {
        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, multiplicity, true, collision1.posZ(), eventWeight);
      }
      if (hasEfficiency) {
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, multiplicity, true, collision1.posZ(), eventWeight);
      }
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processMixedMFT, "Process mixed event MFT-MFT correlations", false);

  void processMftGlobalME(FilteredCollisions const& collisions, FilteredGlobalTracks const& globalTracks, FilteredMFTTracks const& mftTracks, aod::BCsWithTimestamps const&)
  {
    // Same fix as processMixedMFT: bin on getMultiplicity() (driven by cfgMultCorrelationsMask)
    // instead of a hardcoded CentFT0C column, so the mixing pool matches whatever estimator is
    // actually filled into the CorrelationContainer.
    auto getMult = [this](auto const& c) { return getMultiplicity(c); };
    using MixingBinningType = FlexibleBinningPolicy<std::tuple<decltype(getMult)>, aod::collision::PosZ, decltype(getMult)>;
    MixingBinningType binning{{getMult}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(globalTracks, mftTracks);
    Pair<FilteredCollisions, FilteredGlobalTracks, FilteredMFTTracks, MixingBinningType> pairs{binning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    int pairCounter = 0;
    int acceptedPairCounter = 0;
    int newWindowCounter = 0;

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, globalTracks1, collision2, mftTracks2] = *it;
      pairCounter++;

      if (!acceptCollision(collision1) || !collision2.sel8()) {
        continue;
      }
      acceptedPairCounter++;

      const float multiplicity = getMultiplicity(collision1);
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyGlobal != nullptr);
      const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

      int bin = binning.getBin(std::tuple(collision1.posZ(), multiplicity));

      if (it.isNewWindow()) {
        newWindowCounter++;
        loadEfficiency(collision1.template bc_as<aod::BCsWithTimestamps>().timestamp());
        mixedMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
        if (fillReco) {
          mixedMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        }
        if (hasEfficiency) {
          mixedMftGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
        }
        registry.fill(HIST("eventcount_mftglobal_mixed"), bin);
      }
      registry.fill(HIST("trackcount_mftglobal_mixed"), bin, globalTracks1.size(), mftTracks2.size());

      if (fillReco) {
        fillCorrMftGlobal<CorrelationContainer::kCFStepReconstructed>(mixedMftGlobal, globalTracks1, mftTracks2, true, collision1.posZ(), multiplicity, eventWeight);
      }
      if (hasEfficiency) {
        fillCorrMftGlobal<CorrelationContainer::kCFStepCorrected>(mixedMftGlobal, globalTracks1, mftTracks2, true, collision1.posZ(), multiplicity, eventWeight);
      }
    }

    if (cfgVerbosity > 0) {
      LOGF(info, "processMftGlobalME: collisions %d | pairs iterated %d | accepted %d | new windows %d",
           collisions.size(), pairCounter, acceptedPairCounter, newWindowCounter);
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processMftGlobalME, "Process mixed event MFT-vs-global-track correlations", false);

  // FIX #4: Mixed-event global-global (central-barrel hadron-hadron) correlations.
  void processGlobalGlobalME(FilteredCollisions const& collisions, FilteredGlobalTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    auto getMult = [this](auto const& c) { return getMultiplicity(c); };
    using MixingBinningType = FlexibleBinningPolicy<std::tuple<decltype(getMult)>, aod::collision::PosZ, decltype(getMult)>;
    MixingBinningType binning{{getMult}, {axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<FilteredCollisions, FilteredGlobalTracks, MixingBinningType> pairs{binning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;

      if (!acceptCollision(collision1) || !collision2.sel8()) {
        continue;
      }

      const float multiplicity = getMultiplicity(collision1);
      int bin = binning.getBin(std::tuple(collision1.posZ(), multiplicity));
      float eventWeight = 1.0f / it.currentWindowNeighbours();

      if (cfgVerbosity > 0) {
        LOGF(info, "processGlobalGlobalME: bin %d | pair (%d,%d) | Vertex %.2f/%.2f | Mult %.2f", bin, collision1.globalIndex(), collision2.globalIndex(), collision1.posZ(), collision2.posZ(), multiplicity);
      }

      const bool hasEfficiency = (cfg.mEfficiencyGlobal != nullptr);
      const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

      if (it.isNewWindow()) {
        loadEfficiency(collision1.template bc_as<aod::BCsWithTimestamps>().timestamp());
        mixedGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
        if (fillReco) {
          mixedGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
        }
        if (hasEfficiency) {
          mixedGlobal->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
        }
        registry.fill(HIST("eventcount_global_mixed"), bin);
      }
      registry.fill(HIST("trackcount_global_mixed"), bin, tracks1.size(), tracks2.size());

      if (fillReco) {
        fillCorrGlobalGlobal<CorrelationContainer::kCFStepReconstructed>(mixedGlobal, tracks1, tracks2, true, collision1.posZ(), multiplicity, eventWeight);
      }
      if (hasEfficiency) {
        fillCorrGlobalGlobal<CorrelationContainer::kCFStepCorrected>(mixedGlobal, tracks1, tracks2, true, collision1.posZ(), multiplicity, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processGlobalGlobalME, "Process mixed event global-global (central-barrel) correlations", false);

  void processMixedMFTMultSet(FilteredCollisions const& collisions, FilteredMFTTracks const& tracks, aod::BCsWithTimestamps const& bcs)
  {
    if (cfgMultCorrelationsMask == 0) {
      LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when the MultSet process function is in use.");
    }
    processMixedMFT(collisions, tracks, bcs);
  }
  PROCESS_SWITCH(Mftcorrelations, processMixedMFTMultSet, "Process mixed event MFT-MFT correlations with multi-estimator QA", false);

  template <typename TParticle>
  bool isChargedPrimary(TParticle const& p)
  {
    if (!p.isPhysicalPrimary()) {
      return false;
    }
    auto* entry = pdg->GetParticle(p.pdgCode());
    return entry != nullptr && entry->Charge() != 0;
  }

  // ------------------------------------------------------------------ MC efficiency extraction
  Preslice<McEffMFTTracks> perCollisionMFT = aod::fwdtrack::collisionId;

  void processMCEfficiencyMFT(aod::McCollisions::iterator const& mcCollision,
                              aod::McParticles const& mcParticles,
                              McEffCollisions const& collisions,
                              McEffMFTTracks const& tracks)
  {
    registry.fill(HIST("MCEff/hNRecoColl"), collisions.size());
    registry.fill(HIST("MCEff/hCutFlow"), 1);
    if (collisions.size() != 1) {
      return;
    }
    registry.fill(HIST("MCEff/hCutFlow"), 2);
    if (cfgVerbosity > 0) {
      LOGF(info, "processMCEfficiencyMFT: ENTER | mcCollision=%lld | nMcParticles=%d | nCollisions=%d | nTracks=%d",
           mcCollision.globalIndex(), mcParticles.size(), collisions.size(), tracks.size());
    }
    float multiplicity = -1;
    bool passSel8 = false;
    int64_t collGlobalIndex = -1;
    for (const auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("MCEff/hCutFlow"), 3);
      passSel8 = true;
      multiplicity = getMultiplicity(collision); // driven by cfgMultCorrelationsMask
      collGlobalIndex = collision.globalIndex();
    }
    if (!passSel8) {
      return;
    }

    for (const auto& mcParticle : mcParticles) {
      if (!isChargedPrimary(mcParticle) || mcParticle.eta() < cfgMFTCutEtaMin || mcParticle.eta() > cfgMFTCutEtaMax || mcParticle.pt() < cfgMFTCutPtMin || mcParticle.pt() > cfgMFTCutPtMax) {
        continue;
      }
      same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), wrapPhi(mcParticle.phi()));
    }

    auto groupedTracks = tracks.sliceBy(perCollisionMFT, collGlobalIndex);

    std::unordered_set<int64_t> countedPrimaries;
    std::unordered_set<int64_t> countedAll;

    for (const auto& track : groupedTracks) {
      if (track.nClusters() < cfgMFTMinNClusters) {
        continue;
      }
      if (track.has_mcParticle()) {
        const auto& mcParticle = track.mcParticle();
        const int64_t idx = mcParticle.globalIndex();

        if (mcParticle.isPhysicalPrimary() && countedPrimaries.insert(idx).second) {
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries,
                                               mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), wrapPhi(mcParticle.phi()));
        }
        if (countedAll.insert(idx).second) {
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll,
                                               mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), wrapPhi(mcParticle.phi()));
        }
      } else {
        // FIX #6e: wrap phi to [0,2pi) like the MC-particle phi used for the other steps
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake,
                                             track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ(),
                                             RecoDecay::constrainAngle(track.phi(), 0.f));
      }
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processMCEfficiencyMFT, "MC: Extract MFT tracking efficiency", false);

  //-------------------------------------------------------------------MC Efficiency Extraction for Global

  Preslice<McEffGlobalTracks> perCollisionGlobal = aod::track::collisionId;

  void processMCEfficiencyGlobal(aod::McCollisions::iterator const& mcCollision,
                                 aod::McParticles const& mcParticles,
                                 McEffCollisions const& collisions,
                                 McEffGlobalTracks const& tracks)
  {
    registry.fill(HIST("MCEffGlobal/hNRecoColl"), collisions.size());
    registry.fill(HIST("MCEffGlobal/hCutFlow"), 1);
    if (collisions.size() != 1) { // exactly one reconstructed collision per MC collision (as in the MFT efficiency)
      return;
    }
    registry.fill(HIST("MCEffGlobal/hCutFlow"), 2);
    if (cfgVerbosity > 0) {
      LOGF(info, "processMCEfficiencyGlobal: mcCollision=%lld | nMcParticles=%d | nCollisions=%d | nTracks=%d",
           mcCollision.globalIndex(), mcParticles.size(), collisions.size(), tracks.size());
    }

    float multiplicity = -1;
    bool passSel8 = false;
    int64_t collGlobalIndex = -1;
    for (const auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      registry.fill(HIST("MCEffGlobal/hCutFlow"), 3);
      passSel8 = true;
      multiplicity = getMultiplicity(collision); // same estimator as in data (cfgMultCorrelationsMask)
      collGlobalIndex = collision.globalIndex();
    }
    if (!passSel8) {
      return;
    }

    // Generated: charged physical primaries in the global-track acceptance
    for (const auto& mcParticle : mcParticles) {
      if (!isChargedPrimary(mcParticle) ||
          std::abs(mcParticle.eta()) >= cfgGlobalTrackEtaCut ||
          mcParticle.pt() < cfgGlobalTrackPtMin || mcParticle.pt() > cfgGlobalTrackPtMax) {
        continue;
      }
      sameGlobal->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), mcParticle.phi());
    }

    // Reconstructed: global tracks of the selected reco collision (globalTrackFilter already applied)
    auto groupedTracks = tracks.sliceBy(perCollisionGlobal, collGlobalIndex);

    std::unordered_set<int64_t> countedPrimaries; // count each MC particle once (split tracks)
    std::unordered_set<int64_t> countedAll;

    for (const auto& track : groupedTracks) {
      if (track.has_mcParticle()) {
        const auto& mcParticle = track.mcParticle();
        const int64_t idx = mcParticle.globalIndex();

        if (mcParticle.isPhysicalPrimary() && countedPrimaries.insert(idx).second) {
          sameGlobal->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries,
                                                     mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), mcParticle.phi());
        }
        if (countedAll.insert(idx).second) {
          sameGlobal->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll,
                                                     mcParticle.eta(), mcParticle.pt(), 0, multiplicity, mcCollision.posZ(), mcParticle.phi());
        }
      } else {
        // fake track
        sameGlobal->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake,
                                                   track.eta(), track.pt(), 0, multiplicity, mcCollision.posZ(),
                                                   RecoDecay::constrainAngle(track.phi(), 0.f));
      }
    }
  }
  PROCESS_SWITCH(Mftcorrelations, processMCEfficiencyGlobal, "MC: Extract central-barrel (global track) tracking efficiency", false);

  // ------------------------------------------------------------------ MC generator-level correlations (truth reference / closure test)
  using McRecoCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>;

  // -1, 0, +1 from the PDG database charge (0 = neutral or unknown code)
  template <typename TParticle>
  int getMcCharge(TParticle const& p)
  {
    auto* entry = pdg->GetParticle(p.pdgCode());
    if (entry == nullptr) {
      return 0;
    }
    const double q = entry->Charge();
    return (q > 0) - (q < 0);
  }

  // cheap cuts first, PDG lookup last
  template <typename TParticle>
  bool acceptMftTruth(TParticle const& p)
  {
    return p.isPhysicalPrimary() && p.eta() > cfgMFTCutEtaMin && p.eta() < cfgMFTCutEtaMax &&
           p.pt() > cfgMFTCutPtMin && p.pt() < cfgMFTCutPtMax && getMcCharge(p) != 0;
  }

  template <typename TParticle>
  bool acceptGlobalTruth(TParticle const& p)
  {
    return p.isPhysicalPrimary() && std::abs(p.eta()) < cfgGlobalTrackEtaCut &&
           p.pt() > cfgGlobalTrackPtMin && p.pt() < cfgGlobalTrackPtMax && getMcCharge(p) != 0;
  }

  template <typename TParticles, typename TSelector>
  std::vector<int64_t> collectAccepted(TParticles const& particles, TSelector accept)
  {
    std::vector<int64_t> idx;
    idx.reserve(64);
    for (const auto& p : particles) {
      if (accept(p)) {
        idx.push_back(p.globalIndex());
      }
    }
    return idx;
  }

  template <typename TTarget, typename TParticles1, typename TParticles2>
  void fillTruthPairs(TTarget& target, TParticles1 const& particles1, TParticles2 const& particles2,
                      std::vector<int64_t> const& acc1, std::vector<int64_t> const& acc2,
                      float multiplicity, float posZ, float weight, bool sameCollection)
  {
    for (const auto& i1 : acc1) {
      const auto& p1 = particles1.iteratorAt(i1 - particles1.begin().globalIndex());
      target->getTriggerHist()->Fill(CorrelationContainer::kCFStepAll, p1.pt(), multiplicity, posZ, weight);
      for (const auto& i2 : acc2) {
        if (sameCollection && i1 == i2) {
          continue;
        }
        const auto& p2 = particles2.iteratorAt(i2 - particles2.begin().globalIndex());
        if (cfgPtOrder != 0 && p2.pt() >= p1.pt()) {
          continue;
        }
        float deltaPhi = RecoDecay::constrainAngle(p1.phi() - p2.phi(), -o2::constants::math::PIHalf);
        target->getPairHist()->Fill(CorrelationContainer::kCFStepAll, p1.eta() - p2.eta(), p2.pt(), p1.pt(), multiplicity, deltaPhi, posZ, weight);
      }
    }
  }

  // ---------------- same event
  template <typename TTarget, typename TSelector>
  void fillTruthSame(TTarget& target, aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles,
                     soa::SmallGroups<McRecoCollisions> const& collisions, TSelector accept)
  {
    if (collisions.size() != 1) { // same reco-collision requirement as processMCEfficiencyMFT
      return;
    }
    auto collision = collisions.begin();
    if (!collision.sel8()) {
      return;
    }
    const float multiplicity = getMultiplicity(collision);
    const auto acc = collectAccepted(mcParticles, accept);

    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillTruthPairs(target, mcParticles, mcParticles, acc, acc, multiplicity, mcCollision.posZ(), 1.0f, true);
  }

  void processMCSameMFT(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles,
                        soa::SmallGroups<McRecoCollisions> const& collisions)
  {
    fillTruthSame(same, mcCollision, mcParticles, collisions, [this](auto const& p) { return acceptMftTruth(p); });
  }
  PROCESS_SWITCH(Mftcorrelations, processMCSameMFT, "Process MC-truth same event MFT-acceptance correlations (charged primaries)", false);

  void processMCSameGlobal(aod::McCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles,
                           soa::SmallGroups<McRecoCollisions> const& collisions)
  {
    fillTruthSame(sameGlobal, mcCollision, mcParticles, collisions, [this](auto const& p) { return acceptGlobalTruth(p); });
  }
  PROCESS_SWITCH(Mftcorrelations, processMCSameGlobal, "Process MC-truth same event global-acceptance correlations (charged primaries)", false);

  // ---------------- mixed event
  PresliceUnsorted<McRecoCollisions> perMcCollision = aod::mccollisionlabel::mcCollisionId;

  template <typename TTarget, typename TSelector>
  void fillTruthMixed(TTarget& target, aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles,
                      McRecoCollisions const& collisions, TSelector accept)
  {
    // multiplicity of an MC collision = that of its single, sel8-passing reco collision; otherwise -1 (skipped by the binning)
    auto getMult = [&collisions, this](auto const& mcCol) {
      auto grouped = collisions.sliceBy(perMcCollision, mcCol.globalIndex());
      if (grouped.size() != 1) {
        return -1.0f;
      }
      auto c = grouped.begin();
      if (!c.sel8()) {
        return -1.0f;
      }
      return static_cast<float>(getMultiplicity(c));
    };

    using BinningTypeMC = FlexibleBinningPolicy<std::tuple<decltype(getMult)>, aod::mccollision::PosZ, decltype(getMult)>;
    BinningTypeMC binning{{getMult}, {axisVertex, axisMultiplicity}, true};
    auto tuple = std::make_tuple(mcParticles);
    SameKindPair<aod::McCollisions, aod::McParticles, BinningTypeMC> pairs{binning, cfgNoMixedEvents, -1, mcCollisions, tuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [mcCollision1, particles1, mcCollision2, particles2] = *it;
      const float eventWeight = 1.0f / it.currentWindowNeighbours();
      const float multiplicity = getMult(mcCollision1);

      if (it.isNewWindow()) {
        target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
      }
      const auto acc1 = collectAccepted(particles1, accept);
      const auto acc2 = collectAccepted(particles2, accept); // once per pair, not once per trigger
      fillTruthPairs(target, particles1, particles2, acc1, acc2, multiplicity, mcCollision1.posZ(), eventWeight, false);
    }
  }

  void processMCMixedMFT(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, McRecoCollisions const& collisions)
  {
    fillTruthMixed(mixed, mcCollisions, mcParticles, collisions, [this](auto const& p) { return acceptMftTruth(p); });
  }
  PROCESS_SWITCH(Mftcorrelations, processMCMixedMFT, "Process MC-truth mixed event MFT-acceptance correlations (charged primaries)", false);

  void processMCMixedGlobal(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, McRecoCollisions const& collisions)
  {
    fillTruthMixed(mixedGlobal, mcCollisions, mcParticles, collisions, [this](auto const& p) { return acceptGlobalTruth(p); });
  }
  PROCESS_SWITCH(Mftcorrelations, processMCMixedGlobal, "Process MC-truth mixed event global-acceptance correlations (charged primaries)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Mftcorrelations>(cfgc)};
}
