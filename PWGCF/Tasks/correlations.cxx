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

/// \file correlations.cxx
/// \brief task for the correlation calculations with CF-filtered tracks for O2 analysis
/// \author Jan Fiete Grosse-Oetringhaus <jan.fiete.grosse-oetringhaus@cern.ch>, Jasper Parkkila <jasper.parkkila@cern.ch>

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"

#include <TDirectory.h>
#include <TFile.h>
#include <TFormula.h>
#include <TH1F.h>
#include <THn.h>
#include <TVector2.h>

#include <cmath>
#include <experimental/type_traits>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

// NOTE This is a nice idea but will again make it impossible to use subwagon configurations...
// namespace o2::aod
// {
// namespace cfcorr
// {
// DECLARE_SOA_COLUMN(Correction, correction, float);          //! Correction factor for this track
// } // namespace cfcorr
// DECLARE_SOA_TABLE(CFCorrections, "AOD", "CFCORRECTIONS", //! Table which attaches efficiency correction factor to tracks
//                   cfcorreff::Correction);
// } // namespace o2::aod

static constexpr float kCfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};

struct CorrelationTask {
  SliceCache cache;

  // Configuration
  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")

  O2_DEFINE_CONFIGURABLE(cfgPtOrder, int, 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)");
  O2_DEFINE_CONFIGURABLE(cfgTriggerCharge, int, 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgAssociatedCharge, int, 0, "Select on charge of associated particle: 0 = all charged; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgPairCharge, int, 0, "Select on charge of particle pair: 0 = all; 1 = like sign; -1 = unlike sign");
  O2_DEFINE_CONFIGURABLE(cfgCorrelationMethod, int, 0, "Correlation method, 0 = all, 1 = dd, 2 = ddbar");

  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCut, float, -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)");
  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCutMinRadius, float, 0.8f, "Two track cut: radius in m from which two track cuts are applied");
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, int, 0, "0 = OFF and 1 = ON for local efficiency");
  O2_DEFINE_CONFIGURABLE(cfgCentBinsForMC, int, 0, "0 = OFF and 1 = ON for data like multiplicity/centrality bins for MC steps");
  O2_DEFINE_CONFIGURABLE(cfgTrackBitMask, uint16_t, 0, "BitMask for track selection systematics; refer to the enum TrackSelectionCuts in filtering task");
  O2_DEFINE_CONFIGURABLE(cfgMultCorrelationsMask, uint16_t, 0, "Selection bitmask for the multiplicity correlations. This should match the filter selection cfgEstimatorBitMask.")
  O2_DEFINE_CONFIGURABLE(cfgMultCutFormula, std::string, "", "Multiplicity correlations cut formula. A result greater than zero results in accepted event. Parameters: [cFT0C] FT0C centrality, [mFV0A] V0A multiplicity, [mGlob] global track multiplicity, [mPV] PV track multiplicity, [cFT0M] FT0M centrality")

  // Suggested values: Photon: 0.004; K0 and Lambda: 0.005
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {kCfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  O2_DEFINE_CONFIGURABLE(cfgEfficiencyTrigger, std::string, "", "CCDB path to efficiency object for trigger particles")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyAssociated, std::string, "", "CCDB path to efficiency object for associated particles")

  O2_DEFINE_CONFIGURABLE(cfgNoMixedEvents, int, 5, "Number of mixed events per event")

  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")

  O2_DEFINE_CONFIGURABLE(cfgDecayParticleMask, int, 0, "Selection bitmask for the decay particles: 0 = no selection")
  O2_DEFINE_CONFIGURABLE(cfgV0RapidityMax, float, 0.8, "Maximum rapidity for the decay particles (0 = no selection)")
  O2_DEFINE_CONFIGURABLE(cfgMassAxis, int, 0, "Use invariant mass axis (0 = OFF, 1 = ON)")
  O2_DEFINE_CONFIGURABLE(cfgMcTriggerPDGs, std::vector<int>, {}, "MC PDG codes to use exclusively as trigger particles and exclude from associated particles. Empty = no selection.")

  O2_DEFINE_CONFIGURABLE(cfgPtDepMLbkg, std::vector<float>, {}, "pT interval for ML training")
  O2_DEFINE_CONFIGURABLE(cfgPtCentDepMLbkgSel, std::vector<float>, {}, "Bkg ML selection")

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  ConfigurableAxis axisInvMass{"axisInvMass", {VARIABLE_WIDTH, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 5.0}, "invariant mass axis for histograms"};

  // This filter is applied to AOD and derived data (column names are identical)
  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  // This filter is only applied to AOD
  Filter collisionVertexTypeFilter = (aod::collision::flags & static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks)) == static_cast<uint16_t>(aod::collision::CollisionFlagsRun2::Run2VertexerTracks);

  // Track filters
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPt) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true));
  Filter cfTrackFilter = (nabs(aod::cftrack::eta) < cfgCutEta) && (aod::cftrack::pt > cfgCutPt) && ((aod::track::trackType & (uint8_t)cfgTrackBitMask) == (uint8_t)cfgTrackBitMask);

  // MC filters
  Filter cfMCCollisionFilter = nabs(aod::mccollision::posZ) < cfgCutVertex;
  Filter cfMCParticleFilter = (nabs(aod::cfmcparticle::eta) < cfgCutEta) && (aod::cfmcparticle::pt > cfgCutPt); // && (aod::cfmcparticle::sign != 0); //check the sign manually, some specials may be neutral

  // HF filters
  Filter track2pFilter = (nabs(aod::cf2prongtrack::eta) < cfgCutEta) && (aod::cf2prongtrack::pt > cfgCutPt);

  // Output definitions
  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  // persistent caches
  std::vector<float> efficiencyAssociatedCache;
  std::vector<int> p2indexCache;

  std::unique_ptr<TFormula> multCutFormula;
  std::array<uint, aod::cfmultset::NMultiplicityEstimators> multCutFormulaParamIndex;

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>>;
  using AodTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>>;

  using DerivedCollisions = soa::Filtered<aod::CFCollisions>;
  using DerivedTracks = soa::Filtered<aod::CFTracks>;

  void init(o2::framework::InitContext&)
  {
    if ((doprocessSame2ProngDerivedML || doprocessSame2Prong2ProngML || doprocessMixed2ProngDerivedML || doprocessMixed2Prong2ProngML) && (cfgPtDepMLbkg->empty() || cfgPtCentDepMLbkgSel->empty()))
      LOGF(fatal, "cfgPtDepMLbkg or cfgPtCentDepMLbkgSel can not be empty when ML 2-prong selections are used.");
    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, o2::constants::math::TwoPI, "#varphi"}}});
    if (doprocessSame2ProngDerived || doprocessSame2ProngDerivedML || doprocessSame2Prong2Prong || doprocessSame2Prong2ProngML || doprocessMCSameDerived2Prong) {
      registry.add("yieldsTrigger", "multiplicity/centrality vs pT vs eta (triggers)", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
      registry.add("etaphiTrigger", "multiplicity/centrality vs eta vs phi (triggers)", {HistType::kTH3F, {{100, 0, 100, "multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, o2::constants::math::TwoPI, "#varphi"}}});
      const AxisSpec& a = AxisSpec(axisInvMass);
      AxisSpec axisSpecMass = {1000, a.binEdges[0], a.binEdges[a.getNbins()]};
      registry.add("invMass", "2-prong invariant mass (GeV/c^2)", {HistType::kTH3F, {axisSpecMass, axisPtTrigger, axisMultiplicity}});
      if (doprocessSame2Prong2Prong || doprocessSame2Prong2ProngML) {
        registry.add("invMassTwoPart", "2D 2-prong invariant mass (GeV/c^2)", {HistType::kTHnSparseF, {axisSpecMass, axisSpecMass, axisPtTrigger, axisPtAssoc, axisMultiplicity}});
        registry.add("invMassTwoPartDPhi", "2D 2-prong invariant mass (GeV/c^2)", {HistType::kTHnSparseF, {axisSpecMass, axisSpecMass, axisPtTrigger, axisPtAssoc, axisDeltaPhi}});
        registry.add("invMassTwoPartDEta", "2D 2-prong invariant mass (GeV/c^2)", {HistType::kTHnSparseF, {axisSpecMass, axisSpecMass, axisPtTrigger, axisPtAssoc, axisDeltaEta}});
      }
    }
    if (doprocessSameDerivedMultSet) {
      if (cfgMultCorrelationsMask == 0)
        LOGF(fatal, "cfgMultCorrelationsMask can not be 0 when MultSet process functions are in use.");
      std::vector<AxisSpec> multAxes;
      if (cfgMultCorrelationsMask & aod::cfmultset::CentFT0C)
        multAxes.emplace_back(100, 0, 100, "FT0C centrality");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultFV0A)
        multAxes.emplace_back(1000, 0, 100000, "V0A multiplicity");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultNTracksPV)
        multAxes.emplace_back(100, 0, 1000, "Nch PV");
      if (cfgMultCorrelationsMask & aod::cfmultset::MultNTracksGlobal)
        multAxes.emplace_back(100, 0, 1000, "Nch Global");
      registry.add("multCorrelations", "Multiplicity correlations", {HistType::kTHnSparseF, multAxes});
    }
    registry.add("multiplicity", "event multiplicity", {HistType::kTH1F, {{1000, 0, 100, "/multiplicity/centrality"}}});
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

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis;
    std::vector<AxisSpec> userMixingAxis;

    if (cfgMassAxis != 0) {
      userAxis.emplace_back(axisInvMass, "m (GeV/c^2)");
      userMixingAxis.emplace_back(axisInvMass, "m (GeV/c^2)");
    }
    if (doprocessSame2Prong2Prong || doprocessSame2Prong2ProngML)
      userAxis.emplace_back(axisInvMass, "m (GeV/c^2)");
    if (doprocessMixed2Prong2Prong || doprocessMixed2Prong2ProngML)
      userMixingAxis.emplace_back(axisInvMass, "m (GeV/c^2)");

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userMixingAxis));

    same->setTrackEtaCut(cfgCutEta);
    mixed->setTrackEtaCut(cfgCutEta);

    if (!cfgEfficiencyAssociated.value.empty())
      efficiencyAssociatedCache.reserve(512);
    if (doprocessMCEfficiency2Prong || doprocessMCEfficiency2ProngML) {
      p2indexCache.reserve(16);
      if (cfgMcTriggerPDGs->empty())
        LOGF(fatal, "At least one PDG code in {} is to be selected to process 2-prong efficiency.", cfgMcTriggerPDGs.name);
    }

    // o2-ccdb-upload -p Users/jgrosseo/correlations/LHC15o -f /tmp/correction_2011_global.root -k correction

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
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

  template <class p2typeIterator>
  bool passMLScore(const p2typeIterator& track)
  {
    auto it = std::lower_bound(cfgPtDepMLbkg->begin(), cfgPtDepMLbkg->end(), track.pt());
    int idx = std::distance(cfgPtDepMLbkg->begin(), it) - 1;
    return !((track.decay() == 0 && track.mlProbD0()[0] > cfgPtCentDepMLbkgSel->at(idx)) || (track.decay() == 1 && track.mlProbD0bar()[0] > cfgPtCentDepMLbkgSel->at(idx)));
  }

  template <class T>
  using HasMultSet = decltype(std::declval<T&>().multiplicities());

  template <typename TCollision, typename TTracks>
  void fillQA(const TCollision& collision, float multiplicity, const TTracks& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    if constexpr (std::experimental::is_detected<HasMultSet, TCollision>::value) {
      if (std::popcount(cfgMultCorrelationsMask.value) != static_cast<int>(collision.multiplicities().size()))
        LOGF(fatal, "Multiplicity selections (cfgMultCorrelationsMask = 0x%x) do not match the size of the table column (%ld). The histogram filling relies on the preservation of order.", cfgMultCorrelationsMask.value, collision.multiplicities().size());
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
  void fillQA(const TCollision& collision, float multiplicity, const TTracks1& tracks1, const TTracks2& tracks2)
  {
    for (const auto& track1 : tracks1) {
      if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value && std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
        if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track1.decay()))) == 0u)
          continue;
        if constexpr (std::experimental::is_detected<HasMlProbD0, typename TTracks1::iterator>::value) {
          if (!passMLScore(track1))
            continue;
        }
        registry.fill(HIST("invMass"), track1.invMass(), track1.pt(), multiplicity);
        for (const auto& track2 : tracks2) {
          if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks2::iterator>::value && std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
            if (doprocessSame2Prong2Prong || doprocessMixed2Prong2Prong || doprocessSame2Prong2ProngML || doprocessMixed2Prong2ProngML) {
              if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track1.decay()))) == 0u)
                continue;
              if constexpr (std::experimental::is_detected<HasMlProbD0, typename TTracks2::iterator>::value) {
                if (!passMLScore(track2))
                  continue;
              }

              if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks1::iterator>::value) {
                if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks2::iterator>::value) {
                  if (track1.cfTrackProng0Id() == track2.cfTrackProng0Id()) {
                    continue;
                  }
                }
                if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks2::iterator>::value) {
                  if (track1.cfTrackProng0Id() == track2.cfTrackProng1Id()) {
                    continue;
                  }
                }
              }

              if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks1::iterator>::value) {
                if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks2::iterator>::value) {
                  if (track1.cfTrackProng1Id() == track2.cfTrackProng0Id()) {
                    continue;
                  }
                }
                if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks2::iterator>::value) {
                  if (track1.cfTrackProng1Id() == track2.cfTrackProng1Id()) {
                    continue;
                  }
                }
              } // no shared prong for two mothers

              if (cfgCorrelationMethod == 1 && track1.decay() != track2.decay())
                continue;
              if (cfgCorrelationMethod == 2 && track1.decay() == track2.decay())
                continue;
              registry.fill(HIST("invMassTwoPart"), track1.invMass(), track2.invMass(), track1.pt(), track2.pt(), multiplicity);
              registry.fill(HIST("invMassTwoPartDPhi"), track1.invMass(), track2.invMass(), track1.pt(), track2.pt(), TVector2::Phi_0_2pi(track1.phi() - track2.phi() + TMath::Pi() / 2.0) - TMath::Pi() / 2.0);
              if (std::abs(track1.phi() - track2.phi()) < constants::math::PI * 0.5) {
                registry.fill(HIST("invMassTwoPartDEta"), track1.invMass(), track2.invMass(), track1.pt(), track2.pt(), track1.eta() - track2.eta());
              }
            }
          }
        }
      }
      if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
        if (!cfgMcTriggerPDGs->empty() && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track1.pdgCode()) == cfgMcTriggerPDGs->end())
          continue;
      }
      registry.fill(HIST("yieldsTrigger"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphiTrigger"), multiplicity, track1.eta(), track1.phi());
    }
    fillQA(collision, multiplicity, tracks2);
  }

  template <typename TTarget, typename TCollision>
  bool fillCollisionAOD(TTarget target, TCollision collision, float multiplicity)
  {
    target->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);

    if (!collision.alias_bit(kINT7) || !collision.sel7()) {
      return false;
    }

    target->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);

    return true;
  }

  template <CorrelationContainer::CFStep step, typename TTrack>
  bool checkObject(TTrack& track)
  {
    if constexpr (step <= CorrelationContainer::kCFStepAnaTopology) {
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
  using HasMcDecay = decltype(std::declval<T&>().mcDecay());
  template <class T>
  using HasProng0Id = decltype(std::declval<T&>().cfTrackProng0Id());
  template <class T>
  using HasProng1Id = decltype(std::declval<T&>().cfTrackProng1Id());
  template <class T>
  using HasMlProbD0 = decltype(std::declval<T&>().mlProbD0());
  template <class T>
  using HasPartDaugh0Id = decltype(std::declval<T&>().cfParticleDaugh0Id());
  template <class T>
  using HasPartDaugh1Id = decltype(std::declval<T&>().cfParticleDaugh1Id());

  template <class CollType>
  bool passOutlier(CollType const& collision)
  {
    if (cfgMultCutFormula.value.empty())
      return true;
    for (uint i = 0; i < aod::cfmultset::NMultiplicityEstimators; ++i) {
      if ((cfgMultCorrelationsMask.value & (1u << i)) == 0 || multCutFormulaParamIndex[i] == ~0u)
        continue;
      auto estIndex = std::popcount(cfgMultCorrelationsMask.value & ((1u << i) - 1));
      multCutFormula->SetParameter(multCutFormulaParamIndex[i], collision.multiplicities()[estIndex]);
    }
    return multCutFormula->Eval() > 0.0f;
  }

  template <typename T>
  std::tuple<bool, float> getV0Rapidity(const T& track)
  {
    if constexpr (!std::experimental::is_detected<HasDecay, T>::value)
      return {false, 0.0f}; // no decay type, return dummy rapidity
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

    const float E = std::sqrt(p2 + mass * mass);
    return {true, 0.5f * std::log((E + pz) / (E - pz))};
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelations(TTarget target, TTracks1& tracks1, TTracks2& tracks2, float multiplicity, float posZ, int magField, float eventWeight)
  {
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
          if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track1.pdgCode()) == cfgMcTriggerPDGs->end())
            continue;
        } else { // otherwise check the sign against the configuration
          if (cfgTriggerCharge != 0) {
            if (cfgTriggerCharge * track1.sign() < 0)
              continue;
          } else if (track1.sign() == 0) {
            continue; // reject neutral MC particles
          }
        }
      } else if constexpr (std::experimental::is_detected<HasSign, typename TTracks1::iterator>::value) {
        // Check reco objects that have the sign attribute. There are no neutrals to deal with.
        if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0)
          continue;
      }

      if constexpr (std::experimental::is_detected<HasMcDecay, typename TTracks1::iterator>::value) {
        if (((track1.mcDecay() != aod::cf2prongtrack::D0ToPiK) && (track1.mcDecay() != aod::cf2prongtrack::D0barToKPi)) || (track1.decay() & aod::cf2prongmcpart::Prompt) == 0)
          continue;
      } else if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value) {
        if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track1.decay()))) == 0u) {
          continue; // skip particles that do not match the decay mask
        }
        if (cfgV0RapidityMax > 0) {
          auto [t, y] = getV0Rapidity(track1);
          if (t && std::abs(y) > cfgV0RapidityMax)
            continue; // V0s are not allowed to be outside the rapidity range
          registry.fill(HIST("yvspt"), y, track1.pt());
        }
      }

      if constexpr (std::experimental::is_detected<HasPartDaugh0Id, typename TTracks1::iterator>::value) {
        if (track1.cfParticleDaugh0Id() < 0 && track1.cfParticleDaugh1Id() < 0)
          continue; // these we could not match
      }

      if constexpr (std::experimental::is_detected<HasMlProbD0, typename TTracks1::iterator>::value) {
        if (!passMLScore(track1))
          continue;
      } // ML selection

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ);
        }
      }

      if (cfgMassAxis) {
        if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value)
          target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, track1.invMass(), triggerWeight);
        else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
          // TParticlePDG *p = pdg->GetParticle(track1.pdgCode());
          // target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, p->Mass(), triggerWeight);
          target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, 1.8, triggerWeight);
        } else {
          LOGF(fatal, "Can not fill mass axis without invMass column. Disable cfgMassAxis.");
        }
      } else {
        target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, triggerWeight);
      }

      for (const auto& track2 : tracks2) {
        if constexpr (std::is_same<TTracks1, TTracks2>::value) {
          if (track1.globalIndex() == track2.globalIndex()) {
            // LOGF(info, "Track identical: %f | %f | %f || %f | %f | %f", track1.eta(), track1.phi(), track1.pt(),  track2.eta(), track2.phi(), track2.pt());
            continue;
          }
        }
        if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks2::iterator>::value) { // skip those that are specifically chosen to be triggers
          if (!cfgMcTriggerPDGs->empty() && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), track2.pdgCode()) != cfgMcTriggerPDGs->end())
            continue; // TODO: fix cases like MC D0-D0
        }

        // Daughter track and particle checks
        if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfTrackProng0Id()) // do not correlate daughter tracks of the same event
            continue;
        }
        if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfTrackProng1Id()) // do not correlate daughter tracks of the same event
            continue;
        }
        if constexpr (std::experimental::is_detected<HasPartDaugh0Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfParticleDaugh0Id()) // do not correlate daughter particles of the same event
            continue;
        }
        if constexpr (std::experimental::is_detected<HasPartDaugh1Id, typename TTracks1::iterator>::value) {
          if (track2.globalIndex() == track1.cfParticleDaugh1Id()) // do not correlate daughter particles of the same event
            continue;
        }

        if constexpr (step <= CorrelationContainer::kCFStepTracked && !std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (!checkObject<step>(track2)) {
            continue;
          }
        }

        // If decay attributes are found for the second track/particle, we assume 2p-2p correlation
        if constexpr (std::experimental::is_detected<HasMcDecay, typename TTracks2::iterator>::value) {
          if ((((track2.mcDecay()) != aod::cf2prongtrack::D0ToPiK) && ((track2.mcDecay()) != aod::cf2prongtrack::D0barToKPi)) || (track2.decay() & aod::cf2prongmcpart::Prompt) == 0)
            continue;
        } else if constexpr (std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(track2.decay()))) == 0u) {
            continue; // skip particles that do not match the decay mask
          }

          // track2 here is charged hadron so we don't need rapidity cut for this track...this rapidity is only needed for V0
          /*
            if (cfgV0RapidityMax > 0) {
            auto [t, y] = getV0Rapidity(track2);
                  if (t && std::abs(y) > cfgV0RapidityMax)
            continue;
            }*/
        }

        if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value && std::experimental::is_detected<HasDecay, typename TTracks2::iterator>::value) {
          if (cfgCorrelationMethod == 1 && track1.decay() != track2.decay())
            continue;
          if (cfgCorrelationMethod == 2 && track1.decay() == track2.decay())
            continue;
        }

        if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks1::iterator>::value) {
          if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks2::iterator>::value) {
            if (track1.cfTrackProng0Id() == track2.cfTrackProng0Id()) {
              continue;
            }
          }
          if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks2::iterator>::value) {
            if (track1.cfTrackProng0Id() == track2.cfTrackProng1Id()) {
              continue;
            }
          }
        }

        if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks1::iterator>::value) {
          if constexpr (std::experimental::is_detected<HasProng0Id, typename TTracks2::iterator>::value) {
            if (track1.cfTrackProng1Id() == track2.cfTrackProng0Id()) {
              continue;
            }
          }
          if constexpr (std::experimental::is_detected<HasProng1Id, typename TTracks2::iterator>::value) {
            if (track1.cfTrackProng1Id() == track2.cfTrackProng1Id()) {
              continue;
            }
          }
        } // no shared prong for two mothers
        // TODO MC daughters check ^^

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        if constexpr (std::experimental::is_detected<HasSign, typename TTracks2::iterator>::value) {
          // TODO: support for MC D0-D0 case
          if (cfgAssociatedCharge != 0) {
            if (cfgAssociatedCharge * track2.sign() < 0)
              continue;
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

        if constexpr (std::experimental::is_detected<HasMlProbD0, typename TTracks2::iterator>::value) {
          if (!passMLScore(track2))
            continue;
        } // ML selection

        // last param is the weight
        if (cfgMassAxis && (doprocessSame2Prong2Prong || doprocessMixed2Prong2Prong || doprocessSame2Prong2ProngML || doprocessMixed2Prong2ProngML) && !(doprocessSame2ProngDerived || doprocessSame2ProngDerivedML || doprocessMixed2ProngDerived || doprocessMixed2ProngDerivedML)) {
          if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value && std::experimental::is_detected<HasInvMass, typename TTracks2::iterator>::value)
            target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, track2.invMass(), track1.invMass(), associatedWeight);
          else
            LOGF(fatal, "Can not fill mass axis without invMass column. \n no mass for two particles");
        } else if (cfgMassAxis) {
          if constexpr (std::experimental::is_detected<HasInvMass, typename TTracks1::iterator>::value)
            target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, track1.invMass(), associatedWeight);
          else if constexpr (std::experimental::is_detected<HasPDGCode, typename TTracks1::iterator>::value) {
            // TParticlePDG *p = pdg->GetParticle(track1.pdgCode()); //TODO: get the mass for the PDG properly
            target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, 1.8, associatedWeight); // p->Mass()
          } else {
            LOGF(fatal, "Can not fill mass axis without invMass column. Disable cfgMassAxis.");
          }
        } else {
          target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, associatedWeight);
        }
      }
    }
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (cfgEfficiencyTrigger.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyTrigger.value.c_str(), "READ");
        cfg.mEfficiencyTrigger = reinterpret_cast<THn*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      }
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (cfgEfficiencyAssociated.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyAssociated = TFile::Open(cfgEfficiencyAssociated.value.c_str(), "READ");
        cfg.mEfficiencyAssociated = reinterpret_cast<THn*>(fEfficiencyAssociated->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      }
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }
    cfg.efficiencyLoaded = true;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars);
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

    if (fillCollisionAOD(same, collision, multiplicity) == false) {
      return;
    }
    registry.fill(HIST("eventcount_same"), -2);
    fillQA(collision, multiplicity, tracks);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks, tracks, multiplicity, collision.posZ(), getMagneticField(bc.timestamp()), 1.0f);
  }
  PROCESS_SWITCH(CorrelationTask, processSameAOD, "Process same event on AOD", true);

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
    if constexpr (std::experimental::is_detected<HasDecay, typename TTracks1::iterator>::value)
      fillQA(collision, multiplicity, tracks1, tracks2);
    else
      fillQA(collision, multiplicity, tracks1);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
    fillCorrelations<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);

    if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelations<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f);
    }
  }

  void processSameDerived(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSameDerived, "Process same event on derived data", false);

  void processSameDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>>::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks)
  {
    if (!passOutlier(collision))
      return;
    processSameDerivedT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSameDerivedMultSet, "Process same event on derived data with multiplicity sets", false);

  void processSame2ProngDerived(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    processSameDerivedT(collision, p2tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSame2ProngDerived, "Process same event on derived data", false);

  void processSame2ProngDerivedML(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CFTracks> const& tracks, soa::Filtered<soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls>> const& p2tracks)
  {
    processSameDerivedT(collision, p2tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSame2ProngDerivedML, "Process same event on derived data with ML scores", false);

  void processSame2Prong2Prong(DerivedCollisions::iterator const& collision, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    processSameDerivedT(collision, p2tracks, p2tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSame2Prong2Prong, "Process same event on derived data", false);

  void processSame2Prong2ProngML(DerivedCollisions::iterator const& collision, soa::Filtered<soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls>> const& p2tracks)
  {
    processSameDerivedT(collision, p2tracks, p2tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processSame2Prong2ProngML, "Process same event on derived data with ML scores", false);

  using BinningTypeAOD = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentRun2V0M>;
  void processMixedAOD(AodCollisions const& collisions, AodTracks const& tracks, aod::BCsWithTimestamps const&)
  {
    // NOTE legacy function for O2 integration tests. Full version needs derived data

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    BinningTypeAOD configurableBinning{{axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    auto tracksTuple = std::make_tuple(tracks);
    SameKindPair<AodCollisions, AodTracks, BinningTypeAOD> pairs{configurableBinning, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

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
        if (fillCollisionAOD(mixed, collision1, collision1.centRun2V0M()) == false) {
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
  PROCESS_SWITCH(CorrelationTask, processMixedAOD, "Process mixed events on AOD", false);

  template <class CollType, typename... TrackTypes>
  void processMixedDerivedT(CollType const& collisions, TrackTypes&&... tracks)
  {
    auto getMultiplicity =
      [this](auto& col) {
        if constexpr (std::experimental::is_detected<HasMultSet, CollType>::value) {
          if (!passOutlier(col))
            return -1.0f;
        } else {
          (void)this; // fix compile error on unused 'this' capture
        }
        return col.multiplicity();
      };

    using BinningTypeDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::collision::PosZ, decltype(getMultiplicity)>;
    BinningTypeDerived configurableBinningDerived{{getMultiplicity}, {axisVertex, axisMultiplicity}, true}; // true is for 'ignore overflows' (true by default). Underflows and overflows will have bin -1.
    //  Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TA = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<CollType, TA, TB, BinningTypeDerived> pairs{configurableBinningDerived, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

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

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());

        mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepReconstructed);
      }

      // LOGF(info, "Tracks: %d and %d entries", tracks1.size(), tracks2.size());

      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());
      fillCorrelations<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);

      if (cfg.mEfficiencyAssociated || cfg.mEfficiencyTrigger) {
        if (it.isNewWindow()) {
          mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelations<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight);
      }
    }
  }

  void processMixedDerived(DerivedCollisions const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixedDerived, "Process mixed events on derived data", false);

  void processMixedDerivedMultSet(soa::Filtered<soa::Join<aod::CFCollisions, aod::CFMultSets>> const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedT(collisions, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixedDerivedMultSet, "Process mixed events on derived data with multiplicity sets", false);

  void processMixed2ProngDerived(DerivedCollisions const& collisions, DerivedTracks const& tracks, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    processMixedDerivedT(collisions, p2tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixed2ProngDerived, "Process mixed events on derived data", false);

  void processMixed2ProngDerivedML(DerivedCollisions const& collisions, DerivedTracks const& tracks, soa::Filtered<soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls>> const& p2tracks)
  {
    processMixedDerivedT(collisions, p2tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixed2ProngDerivedML, "Process mixed events on derived data with ML scores", false);

  void processMixed2Prong2Prong(DerivedCollisions const& collisions, soa::Filtered<aod::CF2ProngTracks> const& p2tracks)
  {
    processMixedDerivedT(collisions, p2tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixed2Prong2Prong, "Process mixed events on derived data", false);

  void processMixed2Prong2ProngML(DerivedCollisions const& collisions, soa::Filtered<soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls>> const& p2tracks)
  {
    processMixedDerivedT(collisions, p2tracks);
  }
  PROCESS_SWITCH(CorrelationTask, processMixed2Prong2ProngML, "Process mixed events on derived data with ML scores", false);

  // Version with combinations
  /*void processWithCombinations(soa::Join<aod::Collisions, aod::CentRun2V0Ms>::iterator const& collision, aod::BCsWithTimestamps const&, soa::Filtered<aod::Tracks> const& tracks)
  {
    LOGF(info, "Tracks for collision (Combination run): %d", tracks.size());

    // TODO will go to CCDBConfigurable
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    const auto multiplicity = collision.centRun2V0M();

    for (auto track1 = tracks.begin(); track1 != tracks.end(); ++track1) {

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      //       LOGF(info, "TRACK %f %f | %f %f | %f %f", track1.eta(), track1.eta(), track1.phi(), track1.phi2(), track1.pt(), track1.pt());

      same->getTriggerHist()->Fill(CorrelationContainer::kCFStepReconstructed, track1.pt(), multiplicity, collision.posZ());
      // mixed->getTriggerHist()->Fill(eventValues, CorrelationContainer::kCFStepReconstructed);
    }

    for (const auto& [track1, track2] : combinations(tracks, tracks)) {
      // LOGF(info, "Combination %d %d", track1.index(), track2.index());

      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }
      if (cfgAssociatedCharge != 0 && cfgAssociatedCharge * track2.sign() < 0) {
        continue;
      }
      if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0) {
        continue;
      }

      if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
        continue;
      }

      if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, getMagneticField(bc.timestamp()))) {
        continue;
      }

      float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math:    :PIHalf);
      same->getPairHist()->Fill(CorrelationContainer::kCFStepReconstructed,
                                track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, collision.posZ());
      // mixed->getPairHist()->Fill(values, CorrelationContainer::kCFStepReconstructed);
    }
  }
  PROCESS_SWITCH(CorrelationTask, processWithCombinations, "Process same event on AOD with combinations", false);*/

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
    }
    if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), pdgCode) != cfgMcTriggerPDGs->end())
      return 4; // NOTE - if changed, the number in processMCEfficiency2Prong needs to be changed too since we skip the getSpecies call
    else        // The efficiency histogram is hardcoded to contain 5 species. Anything special will have the 4th slot.
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
      if (mcParticle.isPhysicalPrimary() && mcParticle.sign() != 0 && !((doprocessMCEfficiency2Prong || doprocessMCEfficiency2ProngML) && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), mcParticle.pdgCode()) != cfgMcTriggerPDGs->end())) {
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
        if (track.has_cfMCParticle()) {
          const auto& mcParticle = track.cfMCParticle();
          if ((doprocessMCEfficiency2Prong || doprocessMCEfficiency2ProngML) && std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), mcParticle.pdgCode()) != cfgMcTriggerPDGs->end())
            continue; // properly booked by the 2Prong efficiency function, ignore here
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
  PROCESS_SWITCH(CorrelationTask, processMCEfficiency, "MC: Extract efficiencies", false);

  template <class p2type>
  void processMCEfficiency2ProngT(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Join<aod::CFMcParticles, aod::CF2ProngMcParts> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const&, p2type const& p2tracks, Preslice<p2type>& perCollision2Prong)
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
    // Primaries
    p2indexCache.clear();
    for (const auto& mcParticle : mcParticles) {
      if (std::find(cfgMcTriggerPDGs->begin(), cfgMcTriggerPDGs->end(), mcParticle.pdgCode()) != cfgMcTriggerPDGs->end()) {
        if (((mcParticle.mcDecay() != aod::cf2prongtrack::D0ToPiK) && (mcParticle.mcDecay() != aod::cf2prongtrack::D0barToKPi)) || (mcParticle.decay() & aod::cf2prongmcpart::Prompt) == 0)
          continue; // wrong decay channel
        if (mcParticle.cfParticleDaugh0Id() < 0 && mcParticle.cfParticleDaugh1Id() < 0)
          continue; // daughters not found
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::MC, mcParticle.eta(), mcParticle.pt(), 4, multiplicity, mcCollision.posZ());
        p2indexCache.push_back(mcParticle.globalIndex());
      }
    }
    for (const auto& collision : collisions) {
      auto grouped2ProngTracks = p2tracks.sliceBy(perCollision2Prong, collision.globalIndex());

      for (const auto& p2track : grouped2ProngTracks) {
        if (cfgDecayParticleMask != 0 && (cfgDecayParticleMask & (1u << static_cast<uint32_t>(p2track.decay()))) == 0u)
          continue;
        // Check if the mc particles of the prongs are found.
        if constexpr (std::experimental::is_detected<HasMlProbD0, typename p2type::iterator>::value) {
          if (!passMLScore(p2track))
            continue;
        }
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoAll, p2track.eta(), p2track.pt(), 4, multiplicity, mcCollision.posZ());
        auto fillMC2p = [&](const aod::CFTracksWithLabel::iterator& p) -> bool {
          if (!p.has_cfMCParticle())
            return false;
          auto m = std::find_if(p2indexCache.begin(), p2indexCache.end(), [&](const auto& t) -> bool {
            const auto& mcParticle = mcParticles.iteratorAt(t - mcParticles.begin().globalIndex());
            return (p.cfMCParticleId() == mcParticle.cfParticleDaugh0Id() || p.cfMCParticleId() == mcParticle.cfParticleDaugh1Id());
          });
          if (m == p2indexCache.end())
            return false;
          const auto& mcParticle = mcParticles.iteratorAt(*m - mcParticles.begin().globalIndex());
          same->getTrackHistEfficiency()->Fill(CorrelationContainer::RecoPrimaries, mcParticle.eta(), mcParticle.pt(), 4, multiplicity, mcCollision.posZ());
          return true;
        };
        if (p2track.has_cfTrackProng0()) {
          //
          if (const auto& p0 = p2track.template cfTrackProng0_as<aod::CFTracksWithLabel>(); fillMC2p(p0))
            continue;
        }
        if (p2track.has_cfTrackProng1()) {
          if (const auto& p1 = p2track.template cfTrackProng1_as<aod::CFTracksWithLabel>(); fillMC2p(p1))
            continue;
        }

        // fake track
        same->getTrackHistEfficiency()->Fill(CorrelationContainer::Fake, p2track.eta(), p2track.pt(), 4, multiplicity, mcCollision.posZ());
      }
    }
  }

  Preslice<aod::CF2ProngTracks> perCollision2Prong = aod::cftrack::cfCollisionId;
  void processMCEfficiency2Prong(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Join<aod::CFMcParticles, aod::CF2ProngMcParts> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const& tracks, aod::CF2ProngTracks const& p2tracks)
  {
    processMCEfficiency2ProngT(mcCollision, mcParticles, collisions, tracks, p2tracks, perCollision2Prong);
  }
  PROCESS_SWITCH(CorrelationTask, processMCEfficiency2Prong, "MC: Extract efficiencies for 2-prong particles", false);

  Preslice<soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls>> perCollision2ProngML = aod::cftrack::cfCollisionId;
  void processMCEfficiency2ProngML(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Join<aod::CFMcParticles, aod::CF2ProngMcParts> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions, aod::CFTracksWithLabel const& tracks, soa::Join<aod::CF2ProngTracks, aod::CF2ProngTrackmls> const& p2tracks)
  {
    processMCEfficiency2ProngT(mcCollision, mcParticles, collisions, tracks, p2tracks, perCollision2ProngML);
  }
  PROCESS_SWITCH(CorrelationTask, processMCEfficiency2ProngML, "MC: Extract efficiencies for 2-prong particles with ML scores", false);

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

    if (!(doprocessSameDerived || doprocessSameDerivedMultSet || doprocessSame2ProngDerived || doprocessSame2ProngDerivedML || doprocessSame2Prong2Prong || doprocessSame2Prong2ProngML)) {
      if constexpr (std::experimental::is_detected<HasDecay, typename Particles1::iterator>::value)
        fillQA(mcCollision, multiplicity, mcParticles1, mcParticles2);
      else
        fillQA(mcCollision, multiplicity, mcParticles1);
    }

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
    fillCorrelations<CorrelationContainer::kCFStepAll>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    if (collisions.size() == 0) {
      return;
    }

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
    fillCorrelations<CorrelationContainer::kCFStepVertex>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
    fillCorrelations<CorrelationContainer::kCFStepTrackedOnlyPrim>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    same->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
    fillCorrelations<CorrelationContainer::kCFStepTracked>(same, mcParticles1, mcParticles2, multiplicity, mcCollision.posZ(), 0, 1.0f);

    // NOTE kCFStepReconstructed and kCFStepCorrected are filled in processSameDerived
    //      This also means that if a MC collision had several reconstructed vertices (collisions), all of them are filled
  }

  // NOTE SmallGroups includes soa::Filtered always
  void processMCSameDerived(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions) // TODO. For mixed no need to check the daughters since the events are different
  {
    processMCSameDerivedT(mcCollision, mcParticles, mcParticles, collisions);
  }
  PROCESS_SWITCH(CorrelationTask, processMCSameDerived, "Process MC same event on derived data", false);

  void processMCSameDerived2Prong(soa::Filtered<aod::CFMcCollisions>::iterator const& mcCollision, soa::Filtered<soa::Join<aod::CFMcParticles, aod::CF2ProngMcParts>> const& mcParticles2Prong, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions)
  {
    // Subscribe to the MC particles table twice, once joined with the 2prongs and another time without.
    // This is to avoid triggering any 2p-2p specific cases in the templated fillCorrelations().
    processMCSameDerivedT(mcCollision, mcParticles2Prong, mcParticles, collisions);
  }
  PROCESS_SWITCH(CorrelationTask, processMCSameDerived2Prong, "Process MC same event on derived data", false);

  // TODO: add MC 2Prong2Prong functions when needed

  PresliceUnsorted<aod::CFCollisionsWithLabel> collisionPerMCCollision = aod::cfcollision::cfMcCollisionId;
  template <typename... ParticleTypes>
  void processMCMixedDerivedT(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions, ParticleTypes&&... particles)
  {
    bool useMCMultiplicity = (cfgCentBinsForMC == 0);
    auto getMultiplicity =
      [&collisions, &useMCMultiplicity, this](auto& col) {
        if (useMCMultiplicity)
          return col.multiplicity();
        auto groupedCollisions = collisions.sliceBy(collisionPerMCCollision, col.globalIndex());
        if (groupedCollisions.size() == 0)
          return -1.0f;
        return groupedCollisions.begin().multiplicity();
      };

    using BinningTypeMCDerived = FlexibleBinningPolicy<std::tuple<decltype(getMultiplicity)>, aod::mccollision::PosZ, decltype(getMultiplicity)>;
    BinningTypeMCDerived configurableBinning{{getMultiplicity}, {axisVertex, axisMultiplicity}, true};

    // Strictly upper categorised collisions, for cfgNoMixedEvents combinations per bin, skipping those in entry -1
    auto tuple = std::make_tuple(std::forward<ParticleTypes>(particles)...);
    using TA = std::tuple_element<0, decltype(tuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tuple)> - 1, decltype(tuple)>::type;
    Pair<soa::Filtered<aod::CFMcCollisions>, TA, TB, BinningTypeMCDerived> pairs{configurableBinning, cfgNoMixedEvents, -1, mcCollisions, tuple, &cache}; // -1 is the number of the bin to skip

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
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepAll);
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
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepVertex);
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepTrackedOnlyPrim);
        mixed->fillEvent(multiplicity, CorrelationContainer::kCFStepTracked);
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
  PROCESS_SWITCH(CorrelationTask, processMCMixedDerived, "Process MC mixed events on derived data", false);

  void processMCMixedDerived2Prong(soa::Filtered<aod::CFMcCollisions> const& mcCollisions, soa::Filtered<soa::Join<aod::CFMcParticles, aod::CF2ProngMcParts>> const& mcParticles2Prong, soa::Filtered<aod::CFMcParticles> const& mcParticles, soa::Filtered<aod::CFCollisionsWithLabel> const& collisions)
  {
    // Subscribe to the MC particles table twice, once joined with the 2prongs and another time without.
    // This is to avoid triggering any 2p-2p specific cases in the templated fillCorrelations().
    processMCMixedDerivedT(mcCollisions, collisions, mcParticles2Prong, mcParticles);
  }
  PROCESS_SWITCH(CorrelationTask, processMCMixedDerived2Prong, "Process MC mixed events on derived data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrelationTask>(cfgc)};
}
