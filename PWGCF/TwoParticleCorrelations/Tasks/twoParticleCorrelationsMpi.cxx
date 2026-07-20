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
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TDirectory.h>
#include <TFile.h>
#include <TFormula.h>
#include <TGraph.h>
#include <THn.h>
#include <TList.h>
#include <TString.h>

#include <sys/types.h>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <experimental/type_traits>
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
  Configurable<int> cfgNuncSeedEstimator{"cfgNuncSeedEstimator", 0, "Estimator for number of uncorrelated seeds, 0 = bin look up, 1 = interpolation"};

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
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {CfgPairCutDefaults.front().data(), 5, 1, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  Configurable<std::string> cfgEfficiencyTrigger{"cfgEfficiencyTrigger", "", "CCDB path to efficiency object for trigger particles"};
  Configurable<std::string> cfgEfficiencyAssociated{"cfgEfficiencyAssociated", "", "CCDB path to efficiency object for associated particles"};
  Configurable<std::string> cfgNuncSeedsCalibration{"cfgNuncSeedsCalibration", "", "CCDB path to calibration object for number of uncorrelated seeds classification"};

  Configurable<int> cfgNumMixedEvents{"cfgNumMixedEvents", 5, "Number of mixed events per event"};

  Configurable<int> cfgVerbosity{"cfgVerbosity", 1, "Verbosity level (0 = major, 1 = per collision)"};

  Configurable<int> cfgDecayParticleMask{"cfgDecayParticleMask", 0, "Selection bitmask for the decay particles: 0 = no selection"};
  Configurable<float> cfgV0RapidityMax{"cfgV0RapidityMax", 0.8, "Maximum rapidity for the decay particles (0 = no selection)"};
  Configurable<int> cfgMassAxis{"cfgMassAxis", 0, "Use invariant mass axis (0 = OFF, 1 = ON)"};
  Configurable<std::vector<int>> cfgMcTriggerPDGs{"cfgMcTriggerPDGs", {}, "MC PDG codes to use exclusively as trigger particles and exclude from associated particles. Empty = no selection."};

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis axisUncSeeds{"axisUncSeeds", {VARIABLE_WIDTH, 0, 0.2, 0.4, 0.6, 0.8, 1}, "uncorrelated seeds axes in quantiles"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  ConfigurableAxis axisInvMass{"axisInvMass", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 5.0}, "invariant mass axis for histograms"};

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

  struct EventClassifier {
    std::vector<double> nuncSeedsByMultBin;
    std::vector<double> quantileBoundaries;
    std::unique_ptr<TGraph> interpolation;
    bool isLoaded = false;
    int nClasses = 0;
  };

  EventClassifier ec;

  enum CorrelationMethod {
    All = 0,
    Dd,
    Ddbar
  };
  enum EventClassEstimatorMethod {
    BinLookup = 0,
    Interpolation
  };

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  using DerivedCollisions = soa::Filtered<aod::CFCollisions>;
  using DerivedTracks = soa::Filtered<aod::CFTracks>;

  void init(o2::framework::InitContext&)
  {
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
    registry.add("yvspt", "y vs pT", {HistType::kTH2F, {{100, -1, 1, "y"}, {100, 0, 20, "p_{T}"}}}); // y vs pT for all tracks (control histogram)

    const bool useEventClassifier = !cfgNuncSeedsCalibration.value.empty();
    AxisSpec eventClassAxis = useEventClassifier ? AxisSpec(axisUncSeeds) : AxisSpec(axisMultiplicity);
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

    std::string eventClassString = useEventClassifier ? "uncorrelated seed quantile" : "multiplicity / centrality";
    std::vector<AxisSpec> corrAxis;
    corrAxis.reserve(6);
    corrAxis.emplace_back(axisDeltaEta, "#Delta#eta");
    corrAxis.emplace_back(axisPtAssoc, "p_{T} (GeV/c)");
    corrAxis.emplace_back(axisPtTrigger, "p_{T} (GeV/c)");
    corrAxis.emplace_back(eventClassAxis.binEdges, eventClassString);
    corrAxis.emplace_back(axisDeltaPhi, "#Delta#varphi (rad)");
    corrAxis.emplace_back(axisVertex, "z-vtx (cm)");
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis;
    std::vector<AxisSpec> userMixingAxis;

    if (cfgMassAxis != 0) {
      userAxis.emplace_back(axisInvMass, "m (GeV/c^2)");
      userMixingAxis.emplace_back(axisInvMass, "m (GeV/c^2)");
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

    loadEventClassifier(now);
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

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelations(TTarget target, TTracks1& tracks1, TTracks2& tracks2, float multiplicity, float posZ, int magField, float eventWeight)
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
          if (cfgTriggerCharge != 0) {
            if (cfgTriggerCharge * track1.sign() < 0) {
              continue;
            }
          } else if (track1.sign() == 0) {
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
          registry.fill(HIST("yvspt"), y, track1.pt());
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

      if (cfgMassAxis) {
        if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value) {
          target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, track1.invMass(), triggerWeight);
        } else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
          // TParticlePDG *p = pdg->GetParticle(track1.pdgCode());
          // target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, p->Mass(), triggerWeight);
          target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, 1.8, triggerWeight);
        } else {
          LOGF(fatal, "Can not fill mass axis without invMass column. Disable cfgMassAxis.");
        }
      } else {
        target->getTriggerHist()->Fill(step, track1.pt(), containerMultiplicity, posZ, triggerWeight);
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

        if constexpr (std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
          if (cfgAssociatedCharge != 0) {
            if (cfgAssociatedCharge * track2.sign() < 0) {
              continue;
            }
          } else if (track2.sign() == 0) { // mc particles come in neutrals, need to check explicitly
            continue;
          }
        }

        if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value && std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
          if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0) {
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

        // last param is the weight
        if (cfgMassAxis) {
          if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value) {
            target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, track1.invMass(), associatedWeight);
          } else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
            target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, 1.8, associatedWeight); // p->Mass()
          } else {
            LOGF(fatal, "Can not fill mass axis without invMass column. Disable cfgMassAxis.");
          }
        } else {
          target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), containerMultiplicity, deltaPhi, posZ, associatedWeight);
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

  void loadEventClassifier(int64_t timestamp)
  {
    if (cfgNuncSeedsCalibration.value.empty()) {
      return;
    }

    constexpr int TrigBin = 0;
    constexpr int AssocBin = 0;
    auto* list = ccdb->getForTimeStamp<TList>(cfgNuncSeedsCalibration.value, timestamp);

    if (!list) {
      LOGF(fatal, "Could not load event-classifier calibration from CCDB path: %s", cfgNuncSeedsCalibration.value.c_str());
      return;
    }

    auto* hNunc = dynamic_cast<TH1D*>(list->FindObject(Form("calibNuncSeedsVsNch_%d_%d", TrigBin, AssocBin)));
    auto* hQuantiles = dynamic_cast<TH1D*>(list->FindObject(Form("calibNuncSeedsQuantileBoundaries_%d_%d", TrigBin, AssocBin)));

    if (!hNunc) {
      LOGF(fatal, "Missing CCDB object calibNuncSeedsVsNch_%d_%d", TrigBin, AssocBin);
      return;
    }

    if (!hQuantiles) {
      LOGF(fatal, "Missing CCDB object calibNuncSeedsQuantileBoundaries_%d_%d", TrigBin, AssocBin);
      return;
    }

    const AxisSpec multAxis(axisMultiplicity);
    if (hNunc->GetNbinsX() != multAxis.getNbins()) {
      LOGF(fatal,
           "Multiplicity axis mismatch for event classifier: calibration has %d bins, axisMultiplicity has %d bins",
           hNunc->GetNbinsX(),
           multAxis.getNbins());
      return;
    }

    ec.nuncSeedsByMultBin.clear();
    ec.quantileBoundaries.clear();
    ec.nuncSeedsByMultBin.reserve(hNunc->GetNbinsX());
    ec.quantileBoundaries.reserve(hQuantiles->GetNbinsX());

    ec.interpolation = std::make_unique<TGraph>(hNunc->GetNbinsX());
    for (int i = 1; i <= hNunc->GetNbinsX(); ++i) {
      ec.nuncSeedsByMultBin.push_back(hNunc->GetBinContent(i));
      ec.interpolation->SetPoint(i - 1, hNunc->GetXaxis()->GetBinCenter(i), hNunc->GetBinContent(i));
    }

    for (int i = 1; i <= hQuantiles->GetNbinsX(); ++i) {
      ec.quantileBoundaries.push_back(hQuantiles->GetBinContent(i));
    }

    ec.nClasses = static_cast<int>(ec.quantileBoundaries.size()) - 1;
    ec.isLoaded = ec.nClasses > 0;

    if (!ec.isLoaded) {
      LOGF(fatal, "Event-classifier calibration has fewer than two quantile boundaries");
      return;
    }

    LOGF(info,
         "Loaded multiplicity-only event classifier from %s using calibNuncSeedsVsNch_%d_%d and calibNuncSeedsQuantileBoundaries_%d_%d",
         cfgNuncSeedsCalibration.value.c_str(),
         TrigBin,
         AssocBin,
         TrigBin,
         AssocBin);
  }

  int findMultiplicityBin(double multiplicity) const
  {
    const AxisSpec multAxis(axisMultiplicity);
    const auto& edges = multAxis.binEdges;
    const int minEdges = 2;
    if (edges.size() < minEdges || multiplicity < edges.front() || multiplicity >= edges.back()) {
      return -1;
    }

    const auto upper = std::upper_bound(edges.begin(), edges.end(), multiplicity);
    return static_cast<int>(std::distance(edges.begin(), upper)) - 1;
  }

  int classify(double s) const
  {
    const int minEdges = 2;
    if (ec.quantileBoundaries.size() < minEdges) {
      return -1;
    }

    for (int i = 0; i < static_cast<int>(ec.quantileBoundaries.size()) - 1; ++i) {
      if (s >= ec.quantileBoundaries[i] && s < ec.quantileBoundaries[i + 1]) {
        return i;
      }
    }

    return s < ec.quantileBoundaries.front() ? 0 : static_cast<int>(ec.quantileBoundaries.size()) - 2;
  }

  int getEventClass(double multiplicity) const
  {
    if (!ec.isLoaded) {
      return -1;
    }

    double s = -1.0;

    if (cfgNuncSeedEstimator.value == EventClassEstimatorMethod::BinLookup) {
      const int multBin = findMultiplicityBin(multiplicity);

      if (multBin < 0 || multBin >= static_cast<int>(ec.nuncSeedsByMultBin.size())) {
        return -1;
      }

      s = ec.nuncSeedsByMultBin[multBin];
    }

    if (cfgNuncSeedEstimator.value == EventClassEstimatorMethod::Interpolation) {
      if (!ec.interpolation) {
        return -1;
      }

      s = ec.interpolation->Eval(multiplicity);
    }

    if (cfgNuncSeedEstimator.value != EventClassEstimatorMethod::BinLookup && cfgNuncSeedEstimator.value != EventClassEstimatorMethod::Interpolation) {
      return -1;
    }

    return classify(s);
  }

  float getCorrelationContainerMultiplicity(float multiplicity) const
  {
    if (!ec.isLoaded) {
      return multiplicity;
    }

    const int eventClass = getEventClass(multiplicity);
    if (eventClass < 0 || ec.nClasses <= 0) {
      return -1.f;
    }

    return (static_cast<float>(eventClass) + 0.5f) / static_cast<float>(ec.nClasses);
  }

  // Version with explicit nested loop
  void processSameAOD(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracks const& tracks)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    if (cfgVerbosity > 0) {
      LOGF(info, "processSameAOD: Tracks for collision: %d | Vertex: %.1f | INT7: %d | V0M: %.1f", tracks.size(), collision.posZ(), collision.sel7(), collision.centRun2V0M());
    }

    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    loadEfficiency(bc.timestamp());

    const auto multiplicity = collision.centRun2V0M();

    if (!fillCollisionAOD(same, collision, multiplicity)) {
      return;
    }
    registry.fill(HIST("eventcount_same"), -2);
    fillQA(collision, multiplicity, tracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, collision.posZ(), getMagneticField(bc.timestamp()), 1.0f);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processSameAOD, "Process same event on AOD", true);

  template <class CollType, class TTracks1, class TTracks2>
  void processSameDerivedT(CollType const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerivedT: Tracks for collision: %d/%d | Vertex: %.1f | Multiplicity/Centrality: %.1f", tracks1.size(), tracks2.size(), collision.posZ(), collision.multiplicity());
    }
    loadEfficiency(collision.timestamp());

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

    if (fillReco) {
      fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);
    }
    if (hasEfficiency) {
      fillContainerEvent(same, multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);
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
      }
      if (!collision2.alias_bit(kINT7) || !collision2.sel7()) {
        continue;
      }

      registry.fill(HIST("eventcount_mixed"), bin);

      auto bc = collision1.bc_as<aod::BCsWithTimestamps>();

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.centRun2V0M(), collision1.posZ(), getMagneticField(bc.timestamp()), 1.0f / it.currentWindowNeighbours());
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

        if (fillRecoMixed) {
          fillContainerEvent(mixed, collision1.multiplicity(), CorrelationContainer::kCFStepReconstructed);
        }
      }

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());
      if (fillRecoMixed) {
        fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);
      }

      if (hasEfficiencyMixed) {
        if (it.isNewWindow()) {
          fillContainerEvent(mixed, collision1.multiplicity(), CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);
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

  template <class Particles1, class Particles2>
  void processMCSameDerivedT(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, Particles1 const& mcParticles1, Particles2 const& mcParticles2, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
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

    if (!(doprocessSameDerived || doprocessSameDerivedMultSet)) {
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

    // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processSameDerived
    //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
  }

  // NOTE SmallGroups includes soa::Filtered always
  void processMCSameDerived(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions) // TODO. For mixed no need to check the daughters since the events are different
  {
    processMCSameDerivedT(mcCollision, mcParticles, mcParticles, collisions);
  }
  PROCESS_SWITCH(TwoParticleCorrelationsMpi, processMCSameDerived, "Process MC same event on derived data", false);

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
