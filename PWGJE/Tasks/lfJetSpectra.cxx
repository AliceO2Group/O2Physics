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

/// \file lfJetSpectra.cxx
/// \brief pT spectra of K0s and Lambda/AntiLambda in jets and in the underlying event
///
/// \author Francesca Ercolessi <francesca.ercolessi@cern.ch>

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetCandidateUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetReducedDataV0.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>
#include <TProfile2D.h>

#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Jets reconstructed by the PWGJE V0 jet finders (jet-finder-v0-{data,mcd,mcp}-charged): charged tracks (particles) + V0 candidates,
// where the V0 daughters are removed from the track (particle) list and the V0 enters the clustering as a single constituent
using V0ChargedJetsWithConstituents = soa::Join<aod::V0ChargedJets, aod::V0ChargedJetConstituents>;
using MCDV0JetsWithConstituents = soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetConstituents>;
using MCPV0JetsWithConstituents = soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents>;
// Jets geometrically matched between detector and particle level (PWGJE jet-matching-mc-v0-ch), for the jet pT response
using MatchedMCDV0Jets = soa::Join<aod::V0ChargedMCDetectorLevelJets, aod::V0ChargedMCDetectorLevelJetsMatchedToV0ChargedMCParticleLevelJets>;
using MatchedMCPV0JetsWithConstituents = soa::Join<aod::V0ChargedMCParticleLevelJets, aod::V0ChargedMCParticleLevelJetConstituents, aod::V0ChargedMCParticleLevelJetsMatchedToV0ChargedMCDetectorLevelJets>;

struct LfJetSpectra {
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  enum V0Species { kK0s = 0,
                   kLambda,
                   kAntiLambda };

  // Steps of the V0 selection, in the order in which they are applied (bins of the cut flow histogram)
  enum V0SelectionStep { kStepAll = 0,
                         kStepJetFinderAcceptance,
                         kStepNotRejected,
                         kStepSignalFlag,
                         kStepRadius,
                         kStepDcaDauToPV,
                         kStepCosPA,
                         kStepDcaV0Dau,
                         kStepDauEta,
                         kStepCompetingMass,
                         kNSelectionSteps };
  static constexpr std::array<std::string_view, kNSelectionSteps> StepLabels = {"all", "jet finder acc.", "not rejected", "signal flag", "radius", "DCA dau-PV", "cos PA", "DCA V0 dau", "dau #eta", "competing mass"};

  // Folders of the V0 QA histograms: all the stored V0s (before selections) and the selected V0s of each species.
  // The first three entries are also the species names used in the MC histograms
  enum V0QAFolder { kQAK0s = 0,
                    kQALambda,
                    kQAAntiLambda,
                    kQAAllV0s };
  static constexpr std::array<std::string_view, 4> QAFolders = {"K0s", "Lambda", "AntiLambda", "AllV0s"};

  // MC histograms (one per species, axes: multiplicity, pT).
  // Reconstructed V0s: true V0s (PDG of the V0 and of the daughters) passing the selections, at reconstructed pT.
  // Generated V0s: physical primary V0s, at generated pT. In particle-level jets (kGenJet, kGenUE) they are the V0s
  // clustered by the particle-level jet finder (aod::CandidatesV0MCP, all decay channels with v0ChargedDecaysOnly = false),
  // with the same in-jet / UE definition as in data. Otherwise they are taken from aod::JetParticles (all decay channels,
  // so the efficiency includes the branching ratio) and, for the reco jets, "in jet" means within R from the jet axis
  enum MCHistogram { kRecJet = 0,      // reco primary V0s, constituents of a reco jet (efficiency numerator, jet)
                     kRecUE,           // reco primary V0s in the UE of a reco jet (efficiency numerator, UE)
                     kRecJetIncl,      // as kRecJet, primaries + secondaries (feed-down / primary fraction)
                     kRecUEIncl,       // as kRecUE, primaries + secondaries (feed-down / primary fraction)
                     kRecFullEvent,    // reco primary V0s in all the reconstructed events (inclusive efficiency numerator)
                     kRecEventWithJet, // reco primary V0s in the reconstructed events with at least one selected jet
                     kGenRecoJetJet,   // generated V0s within R from the axis of a reco jet (efficiency denominator, jet)
                     kGenRecoJetUE,    // generated V0s in the UE of a reco jet (efficiency denominator, UE)
                     kGenRecoEvent,    // generated V0s in the reconstructed events (inclusive efficiency denominator)
                     kGenEventWithJet, // generated V0s in the reconstructed events with at least one selected jet
                     kGenJet,          // generated V0s, constituents of a particle-level jet
                     kGenUE,           // generated V0s in the UE of a particle-level jet
                     kGenFullEvent,    // generated V0s in all the generated events (signal loss)
                     kNMCHistograms };
  static constexpr std::array<std::string_view, kNMCHistograms> MCHistSuffixes = {"_reconstructed_jet", "_reconstructed_ue", "_reconstructed_jet_incl", "_reconstructed_ue_incl", "_reconstructed_fullevent", "_reconstructed_eventwjet",
                                                                                  "_generated_recojet_jet", "_generated_recojet_ue", "_generated_recoevent", "_generated_eventwjet",
                                                                                  "_generated_jet", "_generated_ue", "_generated_fullevent"};

  // MC origin of a reconstructed V0, from the PDG code of the common mother of the two daughters (0 if there is none).
  // "Fake" V0s (not a K0s, Lambda or AntiLambda) are clustered in the jets as the real ones and add their momentum to the jet
  enum V0Origin { kOriginK0s = 0,
                  kOriginLambda,
                  kOriginAntiLambda,
                  kOriginGammaConversion,
                  kOriginOtherMother,
                  kOriginNoCommonMother,
                  kNV0Origins };
  static constexpr std::array<std::string_view, kNV0Origins> OriginLabels = {"K^{0}_{S}", "#Lambda", "#bar{#Lambda}", "#gamma conversion", "other common mother", "no common mother"};

  // Source of a real V0 (K0s, Lambda, AntiLambda) clustered in a reconstructed jet: physical primary or secondary, classified
  // by the PDG code of the mother. At particle level only the physical primary V0s are clustered: the momentum
  // of the V0s from Xi+- and Omega+- is carried by their mother, a charged primary, while the V0s from Xi0 and from the
  // interactions in the material have no counterpart in the particle-level jets
  enum V0Source { kSourcePrimary = 0,
                  kSourceXiCharged,
                  kSourceXi0,
                  kSourceOmega,
                  kSourceOtherSecondary,
                  kNV0Sources };
  static constexpr std::array<std::string_view, kNV0Sources> SourceLabels = {"primary", "from #Xi^{#pm}", "from #Xi^{0}", "from #Omega^{#pm}", "other secondary"};

  Configurable<bool> fillV0QAHistograms{"fillV0QAHistograms", true, "fill the V0 QA histograms (cut flow and topological variables before/after selections)"};

  // Event selection
  Configurable<std::string> evSel{"evSel", "sel8WithoutTimeFrameBorderCut", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.f, "vertex z cut"};

  // Jet selection. The anti-kT algorithm and the resolution parameter are set in the V0 jet finders
  // (jetAlgorithm = 2, jetRadius = {0.4}); here only the jets with R = rJet are analysed
  Configurable<std::string> jetFinderTaskName{"jetFinderTaskName", "jet-finder-v0-data-charged", "name of the V0 jet finder task in the workflow for data, from which the event, track and V0 selections are read"};
  Configurable<std::string> jetFinderMCDTaskName{"jetFinderMCDTaskName", "jet-finder-v0-mcd-charged", "name of the V0 jet finder task in the workflow for MC detector level"};
  Configurable<std::string> jetFinderMCPTaskName{"jetFinderMCPTaskName", "jet-finder-v0-mcp-charged", "name of the V0 jet finder task in the workflow for MC particle level"};
  Configurable<std::string> derivedDataProducerTaskName{"derivedDataProducerTaskName", "jet-deriveddata-producer", "name of the PWGJE derived data producer in the MC workflow (producing the particle-level V0 table)"};
  Configurable<float> rJet{"rJet", 0.4, "Jet resolution parameter (R)"};
  Configurable<float> minJetPt{"minJetPt", 10.0, "Minimum background-subtracted pt of the jet (GeV/c)"};
  Configurable<float> maxJetEta{"maxJetEta", 0.5, "Maximum |eta| of the jet axis (should be <= 0.9 - R so that jet and UE cones are fully in acceptance)"};

  // V0 selections, applied in the same way to the V0s in jets and in the UE
  struct : ConfigurableGroup {
    std::string prefix = "v0Cuts";
    Configurable<bool> requireV0SignalFlags{"requireV0SignalFlags", false, "require also the species flag of the v0-selector (K0S/Lambda/AntiLambda candidate)"};
    // topological selections
    Configurable<float> minimumV0Radius{"minimumV0Radius", 0.9, "Minimum V0 transverse decay radius (cm)"};
    Configurable<float> maximumV0Radius{"maximumV0Radius", 40.0, "Maximum V0 transverse decay radius (cm)"};
    Configurable<float> dcaDaughterToPVmin{"dcaDaughterToPVmin", 0.06, "Minimum DCA of the V0 daughters to the primary vertex (cm)"};
    Configurable<float> v0cospaMin{"v0cospaMin", 0.99, "Minimum V0 cosine of pointing angle"};
    Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5, "Maximum DCA between the V0 daughters (cm)"};
    Configurable<float> daughterEtaMax{"daughterEtaMax", 0.8, "Maximum |eta| of the V0 daughters"};
    // competing V0 mass rejection
    Configurable<bool> doCompetingMassRejection{"doCompetingMassRejection", true, "Apply the competing V0 mass rejection"};
    Configurable<float> lambdaRejWindow{"lambdaRejWindow", 0.005, "Lambda mass window rejected for K0s candidates (GeV/c^2)"};
    Configurable<float> k0sRejWindow{"k0sRejWindow", 0.010, "K0s mass window rejected for (Anti)Lambda candidates (GeV/c^2)"};
  } v0Cuts;

  // MC options
  Configurable<float> inclusiveEtaMax{"inclusiveEtaMax", 0.8, "Maximum |eta| (generated) of the V0s in the inclusive (full event) MC histograms"};

  // Multiplicity binning
  ConfigurableAxis multBinning{"multBinning", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100}, "Binning of multiplicity axis"};

  Service<o2::framework::O2DatabasePDG> pdgDatabase;

  // phiExclusionMin at or below this value: no phi exclusion region (same convention as the V0 jet finder)
  static constexpr float PhiExclusionUnset = -998.f;

  // Selections of a V0 jet finder, read from its configuration in init() so that the jet finder input
  // (used for the rho estimate), the selected events and the V0s tested for the UE are exactly the ones of the clustering
  struct JetFinderConfig {
    std::string eventSelections = "sel8";
    bool skipMBGapEvents = true;
    bool applyRCTSelections = true;
    std::string trackSelections = "globalTracks";
    std::string particleSelections = "PhysicalPrimary";
    float trackPtMin = 0.15;
    float trackPtMax = 1000.0;
    float trackEtaMin = -0.9;
    float trackEtaMax = 0.9;
    float trackPhiMin = -999;
    float trackPhiMax = 999;
    float phiExclusionMin = -999;
    float phiExclusionMax = 999;
    bool isEmbedding = false;
    float candPtMin = 0.0;
    float candPtMax = 100.0;
    float candYMin = -0.8;
    float candYMax = 0.8;
    float candPDG = static_cast<float>(PDG_t::kK0Short);
    bool useV0SignalFlags = true;
    std::vector<int> eventSelectionBits;
  };
  JetFinderConfig jetFinderConfigDet; // data or MC detector level
  JetFinderConfig jetFinderConfigGen; // MC particle level
  int trackSelection = -1;
  int candIndex = 0;

  std::vector<int> eventSelectionBits;
  std::vector<fastjet::PseudoJet> inputParticles;

  Filter jetCollisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter jetMcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  // Particle-level jet finder input and jets of the MC collision associated to a reconstructed collision (jet pT response)
  Preslice<aod::CandidatesV0MCP> v0MCPPerMcCollision = aod::jv0indices::mcCollisionId;
  Preslice<aod::JetParticles> particlesPerMcCollision = aod::jmcparticle::mcCollisionId;
  Preslice<MatchedMCPV0JetsWithConstituents> mcpJetsPerMcCollision = aod::jet::mcCollisionId;

  // Read the configuration of a V0 jet finder running in the same workflow
  void readJetFinderConfig(InitContext& initContext, const std::string& taskName, JetFinderConfig& config, bool isParticleLevel)
  {
    auto readJetFinderOption = [&](const std::string& optName, auto& value) {
      if (!o2::common::core::getTaskOptionValue(initContext, taskName, optName, value, false)) {
        LOG(fatal) << "Option '" << optName << "' not found in task '" << taskName << "': the V0 jet finder must run in the same workflow";
      }
    };
    readJetFinderOption("eventSelections", config.eventSelections);
    readJetFinderOption("skipMBGapEvents", config.skipMBGapEvents);
    readJetFinderOption("applyRCTSelections", config.applyRCTSelections);
    readJetFinderOption("trackPtMin", config.trackPtMin);
    readJetFinderOption("trackPtMax", config.trackPtMax);
    readJetFinderOption("trackEtaMin", config.trackEtaMin);
    readJetFinderOption("trackEtaMax", config.trackEtaMax);
    readJetFinderOption("trackPhiMin", config.trackPhiMin);
    readJetFinderOption("trackPhiMax", config.trackPhiMax);
    readJetFinderOption("phiExclusionMin", config.phiExclusionMin);
    readJetFinderOption("phiExclusionMax", config.phiExclusionMax);
    readJetFinderOption("candPtMin", config.candPtMin);
    readJetFinderOption("candPtMax", config.candPtMax);
    readJetFinderOption("candYMin", config.candYMin);
    readJetFinderOption("candYMax", config.candYMax);
    if (isParticleLevel) {
      readJetFinderOption("particleSelections", config.particleSelections);
    } else {
      readJetFinderOption("trackSelections", config.trackSelections);
      readJetFinderOption("isEmbedding", config.isEmbedding);
      readJetFinderOption("candPDG", config.candPDG);
      readJetFinderOption("useV0SignalFlags", config.useV0SignalFlags);
    }
    config.eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(config.eventSelections);

    // Sanity checks on the jet finder settings needed by this analysis
    int jetAlgorithm = -1;
    readJetFinderOption("jetAlgorithm", jetAlgorithm);
    if (jetAlgorithm != static_cast<int>(fastjet::antikt_algorithm)) {
      LOG(fatal) << taskName << " must use the anti-kT algorithm (jetAlgorithm = 2), found " << jetAlgorithm;
    }
    // jetRadius is read by hand: an array option is stored in the Variant as a pointer + size, so
    // getTaskOptionValue with a std::vector<double> reinterprets the pointer as a vector (std::bad_alloc)
    std::vector<double> jetRadii;
    for (const auto& device : initContext.services().get<RunningWorkflowInfo const>().devices) {
      if (device.name != taskName) {
        continue;
      }
      for (const auto& option : device.options) {
        if (option.name == "jetRadius") {
          const double* radii = option.defaultValue.get<double*>();
          jetRadii.assign(radii, radii + option.defaultValue.size());
        }
      }
    }
    if (jetRadii.empty()) {
      LOG(fatal) << "Option 'jetRadius' not found in task '" << taskName << "': the V0 jet finder must run in the same workflow";
    }
    if (std::none_of(jetRadii.begin(), jetRadii.end(), [&](double r) { return std::round(r * 100.0) == std::round(rJet * 100.0f); })) {
      LOG(fatal) << "R = " << rJet.value << " is not among the jetRadius values of " << taskName;
    }
    bool saveJetsWithCandidatesOnly = true;
    readJetFinderOption("saveJetsWithCandidatesOnly", saveJetsWithCandidatesOnly);
    if (saveJetsWithCandidatesOnly) {
      LOG(fatal) << taskName << " must run with saveJetsWithCandidatesOnly = false, otherwise jets without V0s are not stored and the jet and UE yields are biased";
    }
    int ghostRepeat = 0;
    readJetFinderOption("ghostRepeat", ghostRepeat);
    if (ghostRepeat == 0) {
      LOG(fatal) << taskName << " must run with ghostRepeat > 0, the jet area is needed for the background subtraction";
    }
  }

  void init(InitContext& initContext)
  {
    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(evSel));

    if (doprocessData && doprocessMCReco) {
      LOG(fatal) << "processData and processMCReco cannot be enabled at the same time";
    }
    if (doprocessData) {
      readJetFinderConfig(initContext, jetFinderTaskName.value, jetFinderConfigDet, false);
    }
    if (doprocessData && (doprocessMCResponse || doprocessMCResponseGen)) {
      LOG(fatal) << "processData and processMCResponse(Gen) cannot be enabled at the same time";
    }
    if (doprocessMCReco || doprocessMCResponse || doprocessMCResponseGen) {
      readJetFinderConfig(initContext, jetFinderMCDTaskName.value, jetFinderConfigDet, false);
    }
    if (doprocessMCGen || doprocessMCResponse || doprocessMCResponseGen) {
      readJetFinderConfig(initContext, jetFinderMCPTaskName.value, jetFinderConfigGen, true);
      // The particle-level V0 table must contain all the decay channels, so that all the generated V0s are clustered by
      // the particle-level jet finder and the generated V0 constituents of the jets include the branching ratio
      bool v0ChargedDecaysOnly = true;
      if (!o2::common::core::getTaskOptionValue(initContext, derivedDataProducerTaskName.value, "v0ChargedDecaysOnly", v0ChargedDecaysOnly, false)) {
        LOG(fatal) << "Option 'v0ChargedDecaysOnly' not found in task '" << derivedDataProducerTaskName.value << "': the PWGJE derived data producer must run in the same workflow";
      }
      if (v0ChargedDecaysOnly) {
        LOG(fatal) << derivedDataProducerTaskName.value << " must run with v0ChargedDecaysOnly = false, otherwise the V0s decaying to neutral particles are not clustered in the particle-level jets";
      }
    }
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(jetFinderConfigDet.trackSelections);
    candIndex = (jetFinderConfigDet.candPDG == static_cast<float>(PDG_t::kLambda0)) ? 1 : 0; // same convention as the V0 jet finder

    const AxisSpec multAxis{multBinning, "FT0M percentile"};
    const AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassK0sAxis{200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassLambdaAxis{200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"};
    const AxisSpec rhoAxis{200, 0., 20., "#rho (GeV/#it{c})"};

    if (doprocessData) {
      // Normalisation and jet QA
      registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {multAxis});
      registryData.add("number_of_jets_vsmultiplicity", "number of jets in data vs multiplicity", HistType::kTH1D, {multAxis});
      registryData.add("rhoPerpCone", "Perpendicular cone #rho", HistType::kTH2F, {multAxis, rhoAxis});
      registryData.add("jetPtRawVsSub", "Jets; #it{p}_{T}^{raw} (GeV/#it{c}); #it{p}_{T}^{sub} (GeV/#it{c})", HistType::kTH2F, {{200, 0., 200.}, {220, -20., 200.}});
      registryData.add("jetPtEtaPhi", "Jets; #it{p}_{T}^{sub} (GeV/#it{c}); #eta; #phi", HistType::kTH3D, {{200, 0., 200.}, {20, -1.f, 1.f}, {18 * 8, 0.f, o2::constants::math::TwoPI}});

      // V0s in jets and in the underlying event
      registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multAxis, ptAxis, invMassLambdaAxis});
      registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multAxis, ptAxis, invMassLambdaAxis});
      registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multAxis, ptAxis, invMassLambdaAxis});
      registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multAxis, ptAxis, invMassLambdaAxis});
      registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multAxis, ptAxis, invMassK0sAxis});
      registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multAxis, ptAxis, invMassK0sAxis});
    }

    // V0 content of the selected jets (number of V0 constituents per jet vs background-subtracted jet pT)
    const AxisSpec jetPtAxis{200, 0., 200., "#it{p}_{T,jet}^{sub} (GeV/#it{c})"};
    const AxisSpec nV0Axis{30, -0.5, 29.5, "number of V0s per jet"};
    if (doprocessData || doprocessMCReco) {
      auto& registry = doprocessData ? registryData : registryMC;
      registry.add("jetV0s/nV0Constituents", "V0 constituents (all V0s clustered by the jet finder)", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registry.add("jetV0s/nK0sSelected", "V0 constituents passing the K^{0}_{S} selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registry.add("jetV0s/nLambdaSelected", "V0 constituents passing the #Lambda selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registry.add("jetV0s/nAntiLambdaSelected", "V0 constituents passing the #bar{#Lambda} selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
    }
    if (doprocessMCReco) {
      const AxisSpec originAxis{kNV0Origins, -0.5, static_cast<double>(kNV0Origins) - 0.5, "MC origin"};
      const AxisSpec hypothesisAxis{3, -0.5, 2.5, "selection hypothesis"};
      const AxisSpec fractionAxis{100, 0., 1., "#Sigma#it{p}_{T}^{fake V0} / #it{p}_{T,jet}^{raw}"};
      registryMC.add("jetV0s/nRealV0Constituents", "real V0 constituents (K^{0}_{S}, #Lambda, #bar{#Lambda})", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registryMC.add("jetV0s/nFakeV0Constituents", "fake V0 constituents", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registryMC.add("jetV0s/nTrueK0sSelected", "true K^{0}_{S} constituents passing the K^{0}_{S} selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registryMC.add("jetV0s/nTrueLambdaSelected", "true #Lambda constituents passing the #Lambda selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registryMC.add("jetV0s/nTrueAntiLambdaSelected", "true #bar{#Lambda} constituents passing the #bar{#Lambda} selection", HistType::kTH2F, {jetPtAxis, nV0Axis});
      registryMC.add("jetV0s/fakeV0PtFraction", "fraction of the jet pT carried by fake V0 constituents", HistType::kTH2F, {jetPtAxis, fractionAxis});
      registryMC.add("jetV0s/fakeSelectedV0PtFraction", "fraction of the jet pT carried by fake V0 constituents passing at least one selection", HistType::kTH2F, {jetPtAxis, fractionAxis});
      // Mean fraction vs jet pT, filled with the exact value (not limited by the binning of the fraction axis)
      registryMC.add("jetV0s/fakeV0PtFractionProfile", "mean fraction of the jet pT carried by fake V0 constituents;#it{p}_{T,jet}^{sub} (GeV/#it{c});#LT #Sigma#it{p}_{T}^{fake V0} / #it{p}_{T,jet}^{raw} #GT", HistType::kTProfile, {jetPtAxis});
      registryMC.add("jetV0s/fakeSelectedV0PtFractionProfile", "mean fraction of the jet pT carried by fake V0 constituents passing at least one selection;#it{p}_{T,jet}^{sub} (GeV/#it{c});#LT #Sigma#it{p}_{T}^{fake V0} / #it{p}_{T,jet}^{raw} #GT", HistType::kTProfile, {jetPtAxis});
      registryMC.add("jetV0s/V0ConstituentOrigin", "MC origin of the V0 constituents", HistType::kTH3F, {jetPtAxis, ptAxis, originAxis});
      registryMC.add("jetV0s/SelectedV0ConstituentOrigin", "MC origin of the V0 constituents passing the selection", HistType::kTH3F, {hypothesisAxis, ptAxis, originAxis});
      auto hOrigin = registryMC.get<TH3>(HIST("jetV0s/V0ConstituentOrigin"));
      auto hSelectedOrigin = registryMC.get<TH3>(HIST("jetV0s/SelectedV0ConstituentOrigin"));
      hSelectedOrigin->GetXaxis()->SetBinLabel(1, "K^{0}_{S}");
      hSelectedOrigin->GetXaxis()->SetBinLabel(2, "#Lambda");
      hSelectedOrigin->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}");
      for (int i = 0; i < kNV0Origins; i++) {
        hOrigin->GetZaxis()->SetBinLabel(i + 1, std::string(OriginLabels[i]).c_str());
        hSelectedOrigin->GetZaxis()->SetBinLabel(i + 1, std::string(OriginLabels[i]).c_str());
      }

      // Real V0 constituents by source (primary / secondary from Xi+-, Xi0, Omega+-, other) and fraction of the (raw) jet pT
      // they carry: the secondary ones are the detector-level counterpart of the momentum of their mother
      const AxisSpec sourceAxis{kNV0Sources, -0.5, static_cast<double>(kNV0Sources) - 0.5, "V0 source"};
      const AxisSpec sourceFractionAxis{100, 0., 1., "#Sigma#it{p}_{T}^{V0} / #it{p}_{T,jet}^{raw}"};
      registryMC.add("jetV0s/V0ConstituentSource", "source of the real V0 constituents", HistType::kTH3F, {jetPtAxis, ptAxis, sourceAxis});
      registryMC.add("jetV0s/V0SourcePtFraction", "fraction of the jet pT carried by the real V0 constituents of each source", HistType::kTH3F, {jetPtAxis, sourceAxis, sourceFractionAxis});
      registryMC.add("jetV0s/V0SourcePtFractionProfile", "mean fraction of the jet pT carried by the real V0 constituents of each source;#it{p}_{T,jet}^{sub} (GeV/#it{c});V0 source;#LT #Sigma#it{p}_{T}^{V0} / #it{p}_{T,jet}^{raw} #GT", HistType::kTProfile2D, {jetPtAxis, sourceAxis});
      for (const auto& axis : {registryMC.get<TH3>(HIST("jetV0s/V0ConstituentSource"))->GetZaxis(), registryMC.get<TH3>(HIST("jetV0s/V0SourcePtFraction"))->GetYaxis(),
                               registryMC.get<TProfile2D>(HIST("jetV0s/V0SourcePtFractionProfile"))->GetYaxis()}) {
        for (int i = 0; i < kNV0Sources; i++) {
          axis->SetBinLabel(i + 1, std::string(SourceLabels[i]).c_str());
        }
      }
    }

    // Jet pT response (detector vs particle level, geometrically matched jets with R = rJet and |eta| <= maxJetEta), with the
    // background subtracted at both levels as in the analysis. The residual is also studied vs the fraction of the
    // particle-level jet pT carried by Xi+- and Omega+-, whose momentum is carried by secondary V0s at detector level
    const AxisSpec jetPtSubDetAxis{220, -20., 200., "#it{p}_{T,jet}^{sub, det} (GeV/#it{c})"};
    const AxisSpec jetPtSubGenAxis{220, -20., 200., "#it{p}_{T,jet}^{sub, gen} (GeV/#it{c})"};
    const AxisSpec residualAxis{200, -2., 2., "(#it{p}_{T,jet}^{det} - #it{p}_{T,jet}^{gen}) / #it{p}_{T,jet}^{gen}"};
    if (doprocessMCResponse) {
      registryMC.add("jetResponse/ptSubDetVsGen", "jet pT response, background subtracted", HistType::kTH2F, {jetPtSubDetAxis, jetPtSubGenAxis});
      registryMC.add("jetResponse/ptRawDetVsGen", "jet pT response, no background subtraction;#it{p}_{T,jet}^{raw, det} (GeV/#it{c});#it{p}_{T,jet}^{raw, gen} (GeV/#it{c})", HistType::kTH2F, {{200, 0., 200.}, {200, 0., 200.}});
      registryMC.add("jetResponse/residualSub", "jet pT residual, background subtracted", HistType::kTH2F, {jetPtSubGenAxis, residualAxis});
      registryMC.add("jetResponse/residualRaw", "jet pT residual, no background subtraction;#it{p}_{T,jet}^{raw, gen} (GeV/#it{c})", HistType::kTH2F, {{200, 0., 200.}, residualAxis});
      registryMC.add("jetResponse/residualRawVsMultiStrangeFraction", "jet pT residual (no background subtraction) vs fraction of the particle-level jet pT carried by #Xi^{#pm} and #Omega^{#pm};#it{p}_{T,jet}^{raw, gen} (GeV/#it{c});#Sigma#it{p}_{T}^{#Xi^{#pm}, #Omega^{#pm}} / #it{p}_{T,jet}^{raw, gen}",
                     HistType::kTH3F, {{20, 0., 200.}, {50, 0., 0.5}, residualAxis});
      registryMC.add("jetResponse/detJets", "detector-level jets in acceptance", HistType::kTH1F, {jetPtSubDetAxis});
      registryMC.add("jetResponse/detJetsMatched", "detector-level jets in acceptance matched to a particle-level jet in acceptance", HistType::kTH1F, {jetPtSubDetAxis});
    }
    if (doprocessMCResponseGen) {
      registryMC.add("jetResponse/genJets", "particle-level jets in acceptance", HistType::kTH1F, {jetPtSubGenAxis});
      registryMC.add("jetResponse/genJetsMatched", "particle-level jets in acceptance matched to a detector-level jet in acceptance", HistType::kTH1F, {jetPtSubGenAxis});
    }

    // MC histograms for efficiency and corrections. The multiplicity is the generated FT0M percentile of the MC collision
    if (doprocessMCReco) {
      registryMC.add("number_of_events_vsmultiplicity_rec", "number of reconstructed events vs multiplicity", HistType::kTH1D, {multAxis});
      registryMC.add("number_of_events_wjet_vsmultiplicity_rec", "number of reconstructed events with at least one selected jet vs multiplicity", HistType::kTH1D, {multAxis});
      registryMC.add("number_of_jets_vsmultiplicity_rec", "number of reconstructed jets vs multiplicity", HistType::kTH1D, {multAxis});
      registryMC.add("rhoPerpCone_rec", "Perpendicular cone #rho, detector level", HistType::kTH2F, {multAxis, rhoAxis});
    }
    if (doprocessMCGen) {
      registryMC.add("number_of_events_vsmultiplicity_gen", "number of generated events vs multiplicity", HistType::kTH1D, {multAxis});
      registryMC.add("number_of_jets_vsmultiplicity_gen", "number of generated jets vs multiplicity", HistType::kTH1D, {multAxis});
      registryMC.add("rhoPerpCone_gen", "Perpendicular cone #rho, particle level", HistType::kTH2F, {multAxis, rhoAxis});
    }
    for (int iHist = 0; iHist < kNMCHistograms; iHist++) {
      const bool isGenLevelHistogram = (iHist == kGenJet || iHist == kGenUE || iHist == kGenFullEvent);
      if ((isGenLevelHistogram && !doprocessMCGen) || (!isGenLevelHistogram && !doprocessMCReco)) {
        continue;
      }
      for (int iSpecies = kK0s; iSpecies <= kAntiLambda; iSpecies++) {
        const std::string name = std::string(QAFolders[iSpecies]) + std::string(MCHistSuffixes[iHist]);
        registryMC.add(name.c_str(), name.c_str(), HistType::kTH2F, {multAxis, ptAxis});
      }
    }

    // V0 QA
    if (fillV0QAHistograms && (doprocessData || doprocessMCReco)) {
      const std::array<std::string, 5> flagLabels = {"all", "K0S", "#Lambda", "#bar{#Lambda}", "rejected"};
      const int nFlagBins = static_cast<int>(flagLabels.size());
      registryQC.add("QA/V0SignalFlags", "v0-selector flags of the stored V0s", HistType::kTH1D, {{nFlagBins, -0.5, nFlagBins - 0.5}});
      auto hFlags = registryQC.get<TH1>(HIST("QA/V0SignalFlags"));
      for (int i = 0; i < nFlagBins; i++) {
        hFlags->GetXaxis()->SetBinLabel(i + 1, flagLabels[i].c_str());
      }

      registryQC.add("QA/V0CutFlow", "V0 cut flow", HistType::kTH2D, {{3, -0.5, 2.5}, {kNSelectionSteps, -0.5, static_cast<double>(kNSelectionSteps) - 0.5}});
      auto hCutFlow = registryQC.get<TH2>(HIST("QA/V0CutFlow"));
      hCutFlow->GetXaxis()->SetBinLabel(1, "K^{0}_{S}");
      hCutFlow->GetXaxis()->SetBinLabel(2, "#Lambda");
      hCutFlow->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}");
      for (int step = 0; step < kNSelectionSteps; step++) {
        hCutFlow->GetYaxis()->SetBinLabel(step + 1, std::string(StepLabels[step]).c_str());
      }

      const AxisSpec radiusAxis{240, 0., 60., "V0 radius (cm)"};
      const AxisSpec dcaToPVAxis{200, 0., 1., "|DCA| to PV (cm)"};
      const AxisSpec cosPAAxis{250, 0.95, 1., "cos PA"};
      const AxisSpec dcaDauAxis{200, 0., 2., "DCA V0 daughters (cm)"};
      const AxisSpec etaAxis{100, -1., 1., "#eta_{daughter}"};
      const AxisSpec alphaAxis{200, -1., 1., "#alpha^{Arm}"};
      const AxisSpec qtAxis{150, 0., 0.3, "#it{q}_{T}^{Arm} (GeV/#it{c})"};
      const AxisSpec deltaMassAxis{200, 0., 0.1, "|m_{competing} - m_{PDG}| (GeV/#it{c}^{2})"};
      for (const auto& folder : QAFolders) {
        const std::string dir = "QA/" + std::string(folder) + "/";
        registryQC.add((dir + "V0Radius").c_str(), "V0 radius", HistType::kTH2F, {ptAxis, radiusAxis});
        registryQC.add((dir + "DcaPosToPV").c_str(), "DCA positive daughter to PV", HistType::kTH2F, {ptAxis, dcaToPVAxis});
        registryQC.add((dir + "DcaNegToPV").c_str(), "DCA negative daughter to PV", HistType::kTH2F, {ptAxis, dcaToPVAxis});
        registryQC.add((dir + "V0CosPA").c_str(), "V0 cosine of pointing angle", HistType::kTH2F, {ptAxis, cosPAAxis});
        registryQC.add((dir + "DcaV0Daughters").c_str(), "DCA between V0 daughters", HistType::kTH2F, {ptAxis, dcaDauAxis});
        registryQC.add((dir + "PosDaughterEta").c_str(), "#eta of the positive daughter", HistType::kTH2F, {ptAxis, etaAxis});
        registryQC.add((dir + "NegDaughterEta").c_str(), "#eta of the negative daughter", HistType::kTH2F, {ptAxis, etaAxis});
        registryQC.add((dir + "Armenteros").c_str(), "Armenteros-Podolanski", HistType::kTH2F, {alphaAxis, qtAxis});
      }
      registryQC.add("QA/K0s/Mass", "K^{0}_{S} invariant mass", HistType::kTH2F, {ptAxis, invMassK0sAxis});
      registryQC.add("QA/Lambda/Mass", "#Lambda invariant mass", HistType::kTH2F, {ptAxis, invMassLambdaAxis});
      registryQC.add("QA/AntiLambda/Mass", "#bar{#Lambda} invariant mass", HistType::kTH2F, {ptAxis, invMassLambdaAxis});
      registryQC.add("QA/K0s/CompetingMassDelta", "K^{0}_{S}: distance of m_{#Lambda} from the #Lambda mass", HistType::kTH2F, {ptAxis, deltaMassAxis});
      registryQC.add("QA/Lambda/CompetingMassDelta", "#Lambda: distance of m_{K^{0}_{S}} from the K^{0}_{S} mass", HistType::kTH2F, {ptAxis, deltaMassAxis});
      registryQC.add("QA/AntiLambda/CompetingMassDelta", "#bar{#Lambda}: distance of m_{K^{0}_{S}} from the K^{0}_{S} mass", HistType::kTH2F, {ptAxis, deltaMassAxis});
    }
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Jet finder input, background and jet selection (common to data, MC detector level and MC particle level)

  // Event selection: the one of this task plus the one of the jet finder, so that events without jets are only
  // the ones in which the jet finder did not find any jet
  template <typename T>
  bool isSelectedCollision(T const& collision, JetFinderConfig const& config)
  {
    return jetderiveddatautilities::selectCollision(collision, eventSelectionBits) &&
           jetderiveddatautilities::selectCollision(collision, config.eventSelectionBits, config.skipMBGapEvents, config.applyRCTSelections);
  }

  template <typename T>
  bool passedJetFinderKinematics(T const& particle, JetFinderConfig const& config)
  {
    if (particle.pt() < config.trackPtMin || particle.pt() >= config.trackPtMax ||
        particle.eta() < config.trackEtaMin || particle.eta() > config.trackEtaMax ||
        particle.phi() < config.trackPhiMin || particle.phi() > config.trackPhiMax) {
      return false;
    }
    if (config.phiExclusionMin > PhiExclusionUnset && particle.phi() > config.phiExclusionMin && particle.phi() < config.phiExclusionMax) {
      return false;
    }
    return true;
  }

  // Track selection of the V0 jet finder (kinematic filter + track selection + removal of the V0 daughters)
  template <typename T, typename U>
  bool isJetFinderTrack(T const& track, U const& v0s)
  {
    if (!passedJetFinderKinematics(track, jetFinderConfigDet)) {
      return false;
    }
    if (!jetderiveddatautilities::selectTrack(track, trackSelection, jetFinderConfigDet.isEmbedding)) {
      return false;
    }
    for (const auto& v0 : v0s) {
      if (jetcandidateutilities::isDaughterTrack(track, v0)) {
        return false;
      }
    }
    return true;
  }

  // Collect the indices of the V0 pseudojets (status candidate) in the jet finder input
  std::unordered_set<int> getJetFinderV0s()
  {
    std::unordered_set<int> jetFinderV0s;
    for (const auto& particle : inputParticles) {
      if (particle.user_info<fastjetutilities::fastjet_user_info>().getStatus() == JetConstituentStatus::candidate) {
        jetFinderV0s.insert(particle.user_info<fastjetutilities::fastjet_user_info>().getIndex());
      }
    }
    return jetFinderV0s;
  }

  // Rebuild the detector-level jet finder input (V0s + tracks without V0 daughters), with the same selections used in the clustering.
  // Returns the V0s accepted by the jet finder: only these can be jet constituents, so only these are tested for the UE
  template <typename TV0s, typename TTracks>
  std::unordered_set<int> buildDetectorLevelInput(TV0s const& v0s, TTracks const& tracks)
  {
    inputParticles.clear();
    jetfindingutilities::analyseV0s(inputParticles, v0s, jetFinderConfigDet.candPtMin, jetFinderConfigDet.candPtMax, jetFinderConfigDet.candYMin, jetFinderConfigDet.candYMax, candIndex, jetFinderConfigDet.useV0SignalFlags);
    const auto jetFinderV0s = getJetFinderV0s();
    for (const auto& track : tracks) {
      if (isJetFinderTrack(track, v0s)) {
        fastjetutilities::fillTracks(track, inputParticles, track.globalIndex());
      }
    }
    return jetFinderV0s;
  }

  // Rebuild the particle-level jet finder input (V0s + charged particles without V0 daughters), with the same selections used in the clustering.
  // The secondary V0s (e.g. from Xi/Omega decays) are not clustered (analyseV0s keeps only the physical primary MC V0s): their
  // momentum is already carried by the charged primary mother
  template <typename TV0s, typename TParticles>
  std::unordered_set<int> buildParticleLevelInput(TV0s const& v0s, TParticles const& particles)
  {
    inputParticles.clear();
    jetfindingutilities::analyseV0s(inputParticles, v0s, jetFinderConfigGen.candPtMin, jetFinderConfigGen.candPtMax, jetFinderConfigGen.candYMin, jetFinderConfigGen.candYMax, 0, false);
    const auto jetFinderV0s = getJetFinderV0s();
    std::vector<fastjet::PseudoJet> chargedParticles;
    jetfindingutilities::analyseParticles<true>(chargedParticles, jetFinderConfigGen.particleSelections, static_cast<int>(JetType::charged), particles, pdgDatabase, &v0s);
    for (const auto& particle : chargedParticles) {
      if (passedJetFinderKinematics(particle, jetFinderConfigGen)) {
        inputParticles.push_back(particle);
      }
    }
    return jetFinderV0s;
  }

  // Background: rho from the two cones perpendicular to the leading jet (PWGJE jetutilities::estimateRhoPerpCone),
  // computed on the jet finder input. Returns a negative value if there are no jets with R = rJet
  template <typename TJets>
  float estimateRho(TJets const& jets)
  {
    float leadingJetPt = -1.;
    float leadingJetEta = 0.;
    float leadingJetPhi = 0.;
    for (const auto& jet : jets) {
      if (jet.r() == std::round(rJet * 100.0f) && jet.pt() > leadingJetPt) {
        leadingJetPt = jet.pt();
        leadingJetEta = jet.eta();
        leadingJetPhi = jet.phi();
      }
    }
    if (leadingJetPt < 0.) {
      return -1.;
    }
    const auto [rho, rhoM] = jetutilities::estimateRhoPerpCone(inputParticles, fastjet::PtYPhiM(leadingJetPt, leadingJetEta, leadingJetPhi), rJet);
    return rho;
  }

  // Jet pT corrected for the background as pT,sub = pT,raw - rho * A
  template <typename TJet>
  float subtractedPt(TJet const& jet, float rho)
  {
    return jet.pt() - rho * jet.area();
  }

  template <typename TJet>
  bool isJetAboveThreshold(TJet const& jet, float rho)
  {
    return jet.r() == std::round(rJet * 100.0f) && subtractedPt(jet, rho) >= minJetPt;
  }

  template <typename TJet>
  bool isSelectedJet(TJet const& jet, float rho)
  {
    return isJetAboveThreshold(jet, rho) && std::abs(jet.eta()) <= maxJetEta;
  }

  // Jet acceptance of the jet pT response (no pT threshold)
  template <typename TJet>
  bool isInJetAcceptance(TJet const& jet)
  {
    return jet.r() == std::round(rJet * 100.0f) && std::abs(jet.eta()) <= maxJetEta;
  }

  // Jets above the pT threshold after subtraction: their V0 constituents and their axes are used to exclude from the UE
  // of all the other jets the V0s inside another jet (independently of the jet eta cut, to also catch jets partially
  // overlapping with the UE cones)
  template <typename TJets>
  void collectJetsAboveThreshold(TJets const& jets, float rho, std::unordered_set<int>& v0sInJets, std::vector<std::array<float, 2>>& jetAxes)
  {
    for (const auto& jet : jets) {
      if (!isJetAboveThreshold(jet, rho)) {
        continue;
      }
      jetAxes.push_back({jet.eta(), jet.phi()});
      for (const auto& v0Id : jet.candidatesIds()) {
        v0sInJets.insert(v0Id);
      }
    }
  }

  // UE region: two cones of radius R centred on the axes perpendicular to the jet in azimuth (same eta, phi_jet +- pi/2)
  template <typename TJet>
  bool isInUE(float eta, float phi, TJet const& jet)
  {
    return jetutilities::deltaR(eta, phi, jet.eta(), jet.phi() + o2::constants::math::PIHalf) < rJet ||
           jetutilities::deltaR(eta, phi, jet.eta(), jet.phi() - o2::constants::math::PIHalf) < rJet;
  }

  bool isInAnyJetCone(float eta, float phi, std::vector<std::array<float, 2>> const& jetAxes)
  {
    return std::any_of(jetAxes.begin(), jetAxes.end(), [&](const auto& axis) { return jetutilities::deltaR(eta, phi, axis[0], axis[1]) < rJet; });
  }

  // ---------------------------------------------------------------------------------------------------------------
  // V0 selection

  // Distance of the competing mass hypothesis from the PDG mass, used in the competing V0 mass rejection
  template <typename T>
  float competingMassDelta(T const& v0, V0Species species)
  {
    if (species == kK0s) {
      return std::min(std::abs(v0.mLambda() - o2::constants::physics::MassLambda), std::abs(v0.mAntiLambda() - o2::constants::physics::MassLambda));
    }
    return std::abs(v0.mK0Short() - o2::constants::physics::MassK0Short);
  }

  // V0 selection for a given species hypothesis, applied step by step. Returns the last step passed
  // (kStepCompetingMass = all the selections passed). It is used both for the selection and for the cut flow QA
  template <typename T>
  int lastPassedSelectionStep(T const& v0, V0Species species, bool isJetFinderV0)
  {
    if (!isJetFinderV0) {
      return kStepAll;
    }
    if (v0.isRejectedCandidate()) {
      return kStepJetFinderAcceptance;
    }
    if (v0Cuts.requireV0SignalFlags) {
      const bool hasFlag = (species == kK0s) ? v0.isK0SCandidate() : (species == kLambda) ? v0.isLambdaCandidate()
                                                                                          : v0.isAntiLambdaCandidate();
      if (!hasFlag) {
        return kStepNotRejected;
      }
    }
    if (v0.v0radius() < v0Cuts.minimumV0Radius || v0.v0radius() > v0Cuts.maximumV0Radius) {
      return kStepSignalFlag;
    }
    if (std::abs(v0.dcapostopv()) < v0Cuts.dcaDaughterToPVmin || std::abs(v0.dcanegtopv()) < v0Cuts.dcaDaughterToPVmin) {
      return kStepRadius;
    }
    if (v0.v0cosPA() < v0Cuts.v0cospaMin) {
      return kStepDcaDauToPV;
    }
    if (std::abs(v0.dcaV0daughters()) > v0Cuts.dcaV0DaughtersMax) {
      return kStepCosPA;
    }
    if (std::abs(v0.positiveeta()) > v0Cuts.daughterEtaMax || std::abs(v0.negativeeta()) > v0Cuts.daughterEtaMax) {
      return kStepDcaV0Dau;
    }
    if (v0Cuts.doCompetingMassRejection && competingMassDelta(v0, species) < ((species == kK0s) ? v0Cuts.lambdaRejWindow : v0Cuts.k0sRejWindow)) {
      return kStepDauEta;
    }
    return kStepCompetingMass;
  }

  // Full selection of a V0 (accepted by the jet finder) for a given species hypothesis
  template <typename T>
  bool passedV0Selection(T const& v0, V0Species species)
  {
    return lastPassedSelectionStep(v0, species, true) == kStepCompetingMass;
  }

  // ---------------------------------------------------------------------------------------------------------------
  // V0 QA

  // QA of the topological variables, for all the V0s (folder "AllV0s") or for the selected V0s of one species
  template <int iFolder, typename T>
  void fillV0QAVariables(T const& v0)
  {
    constexpr auto Folder = HIST("QA/") + HIST(QAFolders[iFolder]);
    registryQC.fill(Folder + HIST("/V0Radius"), v0.pt(), v0.v0radius());
    registryQC.fill(Folder + HIST("/DcaPosToPV"), v0.pt(), std::abs(v0.dcapostopv()));
    registryQC.fill(Folder + HIST("/DcaNegToPV"), v0.pt(), std::abs(v0.dcanegtopv()));
    registryQC.fill(Folder + HIST("/V0CosPA"), v0.pt(), v0.v0cosPA());
    registryQC.fill(Folder + HIST("/DcaV0Daughters"), v0.pt(), std::abs(v0.dcaV0daughters()));
    registryQC.fill(Folder + HIST("/PosDaughterEta"), v0.pt(), v0.positiveeta());
    registryQC.fill(Folder + HIST("/NegDaughterEta"), v0.pt(), v0.negativeeta());
    registryQC.fill(Folder + HIST("/Armenteros"), v0.alpha(), v0.qtarm());
    if constexpr (iFolder == kQAK0s) {
      registryQC.fill(Folder + HIST("/Mass"), v0.pt(), v0.mK0Short());
      registryQC.fill(Folder + HIST("/CompetingMassDelta"), v0.pt(), competingMassDelta(v0, kK0s));
    } else if constexpr (iFolder == kQALambda) {
      registryQC.fill(Folder + HIST("/Mass"), v0.pt(), v0.mLambda());
      registryQC.fill(Folder + HIST("/CompetingMassDelta"), v0.pt(), competingMassDelta(v0, kLambda));
    } else if constexpr (iFolder == kQAAntiLambda) {
      registryQC.fill(Folder + HIST("/Mass"), v0.pt(), v0.mAntiLambda());
      registryQC.fill(Folder + HIST("/CompetingMassDelta"), v0.pt(), competingMassDelta(v0, kAntiLambda));
    }
  }

  // Inclusive V0 QA (once per V0 and per event, independently of the jets): v0-selector flags, topological variables
  // before the selections, cut flow and topological variables after the selections for each species
  template <typename T>
  void fillV0QA(T const& v0, bool isJetFinderV0)
  {
    registryQC.fill(HIST("QA/V0SignalFlags"), 0);
    if (v0.isK0SCandidate()) {
      registryQC.fill(HIST("QA/V0SignalFlags"), 1);
    }
    if (v0.isLambdaCandidate()) {
      registryQC.fill(HIST("QA/V0SignalFlags"), 2);
    }
    if (v0.isAntiLambdaCandidate()) {
      registryQC.fill(HIST("QA/V0SignalFlags"), 3);
    }
    if (v0.isRejectedCandidate()) {
      registryQC.fill(HIST("QA/V0SignalFlags"), 4);
    }
    fillV0QAVariables<kQAAllV0s>(v0);

    for (const auto& species : {kK0s, kLambda, kAntiLambda}) {
      const int lastStep = lastPassedSelectionStep(v0, species, isJetFinderV0);
      for (int step = kStepAll; step <= lastStep; step++) {
        registryQC.fill(HIST("QA/V0CutFlow"), species, step);
      }
      if (lastStep != kStepCompetingMass) {
        continue;
      }
      if (species == kK0s) {
        fillV0QAVariables<kQAK0s>(v0);
      } else if (species == kLambda) {
        fillV0QAVariables<kQALambda>(v0);
      } else {
        fillV0QAVariables<kQAAntiLambda>(v0);
      }
    }
  }

  // ---------------------------------------------------------------------------------------------------------------
  // Data

  // Fill the V0 histograms for one region (in-jet or UE). The same selection is applied in the two regions.
  // A V0 can satisfy more than one species hypothesis
  template <bool isInJet, typename T>
  void fillV0(T const& v0, float multiplicity)
  {
    if (passedV0Selection(v0, kK0s)) {
      if constexpr (isInJet) {
        registryData.fill(HIST("K0s_in_jet"), multiplicity, v0.pt(), v0.mK0Short());
      } else {
        registryData.fill(HIST("K0s_in_ue"), multiplicity, v0.pt(), v0.mK0Short());
      }
    }
    if (passedV0Selection(v0, kLambda)) {
      if constexpr (isInJet) {
        registryData.fill(HIST("Lambda_in_jet"), multiplicity, v0.pt(), v0.mLambda());
      } else {
        registryData.fill(HIST("Lambda_in_ue"), multiplicity, v0.pt(), v0.mLambda());
      }
    }
    if (passedV0Selection(v0, kAntiLambda)) {
      if constexpr (isInJet) {
        registryData.fill(HIST("AntiLambda_in_jet"), multiplicity, v0.pt(), v0.mAntiLambda());
      } else {
        registryData.fill(HIST("AntiLambda_in_ue"), multiplicity, v0.pt(), v0.mAntiLambda());
      }
    }
  }

  void processData(soa::Filtered<aod::JetCollisions>::iterator const& jcoll, V0ChargedJetsWithConstituents const& v0jets, aod::CandidatesV0Data const& v0s, aod::JetTracks const& tracks)
  {
    if (!isSelectedCollision(jcoll, jetFinderConfigDet)) {
      return;
    }
    const float multiplicity = jcoll.centFT0M();
    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);

    const auto jetFinderV0s = buildDetectorLevelInput(v0s, tracks);

    // Inclusive V0 QA, filled once per V0 for all the selected events (also without jets)
    if (fillV0QAHistograms) {
      for (const auto& v0 : v0s) {
        fillV0QA(v0, jetFinderV0s.contains(v0.globalIndex()));
      }
    }

    const float rho = estimateRho(v0jets);
    if (rho < 0.) {
      return;
    }
    registryData.fill(HIST("rhoPerpCone"), multiplicity, rho);

    std::unordered_set<int> v0sInJets;
    std::vector<std::array<float, 2>> jetAxes;
    collectJetsAboveThreshold(v0jets, rho, v0sInJets, jetAxes);

    for (const auto& jet : v0jets) {
      if (jet.r() == std::round(rJet * 100.0f)) {
        registryData.fill(HIST("jetPtRawVsSub"), jet.pt(), subtractedPt(jet, rho));
      }
      if (!isSelectedJet(jet, rho)) {
        continue;
      }
      registryData.fill(HIST("number_of_jets_vsmultiplicity"), multiplicity);
      registryData.fill(HIST("jetPtEtaPhi"), subtractedPt(jet, rho), jet.eta(), jet.phi());
      fillJetV0Content<false, aod::CandidatesV0Data>(jet, subtractedPt(jet, rho));

      // V0-jet matching, in-jet region: a V0 is "in jet" if it is a constituent of the jet.
      // The V0 jet finder clusters each V0 as a single pseudojet and removes its daughters from the track list,
      // so the V0 is assigned to at most one jet and its daughters are never counted as separate jet constituents
      for (const auto& v0 : jet.candidates_as<aod::CandidatesV0Data>()) {
        fillV0<true>(v0, multiplicity);
      }

      // V0-jet matching, UE region: V0s in the perpendicular cones. The V0s are those accepted by the jet finder
      // (same selection as the in-jet V0s), excluding the V0s that are constituents of another jet above threshold
      for (const auto& v0 : v0s) {
        if (!jetFinderV0s.contains(v0.globalIndex()) || v0sInJets.contains(v0.globalIndex())) {
          continue;
        }
        if (isInUE(v0.eta(), v0.phi(), jet)) {
          fillV0<false>(v0, multiplicity);
        }
      }
    }
  }
  PROCESS_SWITCH(LfJetSpectra, processData, "V0s in jets and in the UE in data", true);

  // ---------------------------------------------------------------------------------------------------------------
  // MC

  template <int iSpecies, int iHist>
  void fillMCHistogram(float multiplicity, float pt, float weight)
  {
    registryMC.fill(HIST(QAFolders[iSpecies]) + HIST(MCHistSuffixes[iHist]), multiplicity, pt, weight);
  }

  template <int iHist>
  void fillMCHistogram(V0Species species, float multiplicity, float pt, float weight)
  {
    if (species == kK0s) {
      fillMCHistogram<kK0s, iHist>(multiplicity, pt, weight);
    } else if (species == kLambda) {
      fillMCHistogram<kLambda, iHist>(multiplicity, pt, weight);
    } else {
      fillMCHistogram<kAntiLambda, iHist>(multiplicity, pt, weight);
    }
  }

  // Reconstructed V0 matched to a true V0 of the given species: PDG code of the V0 and of the two daughters
  template <typename T>
  bool isTrueV0(T const& v0, V0Species species)
  {
    switch (species) {
      case kK0s:
        return v0.pdgCode() == PDG_t::kK0Short && v0.pdgCodePositive() == PDG_t::kPiPlus && v0.pdgCodeNegative() == PDG_t::kPiMinus;
      case kLambda:
        return v0.pdgCode() == PDG_t::kLambda0 && v0.pdgCodePositive() == PDG_t::kProton && v0.pdgCodeNegative() == PDG_t::kPiMinus;
      case kAntiLambda:
        return v0.pdgCode() == PDG_t::kLambda0Bar && v0.pdgCodePositive() == PDG_t::kPiPlus && v0.pdgCodeNegative() == PDG_t::kProtonBar;
    }
    return false;
  }

  // Species of a generated V0 (-1 if not a K0s, Lambda or AntiLambda)
  template <typename T>
  int genV0Species(T const& v0)
  {
    switch (v0.pdgCode()) {
      case PDG_t::kK0Short:
        return kK0s;
      case PDG_t::kLambda0:
        return kLambda;
      case PDG_t::kLambda0Bar:
        return kAntiLambda;
    }
    return -1;
  }

  // MC origin of a reconstructed V0, from the PDG code of the common mother of the two daughters (0 if there is none)
  template <typename T>
  int v0Origin(T const& v0)
  {
    switch (v0.pdgCode()) {
      case PDG_t::kK0Short:
        return kOriginK0s;
      case PDG_t::kLambda0:
        return kOriginLambda;
      case PDG_t::kLambda0Bar:
        return kOriginAntiLambda;
      case PDG_t::kGamma:
        return kOriginGammaConversion;
      case 0:
        return kOriginNoCommonMother; // the two daughters have no common mother, e.g. one of them is a primary particle and the other is a secondary from a weak decay
    }
    return kOriginOtherMother; // the two daughters have a common mother, but it is not a K0s, Lambda, AntiLambda or gamma, e.g. material or multi-particle decay or
  }

  // Source of a real V0 (see V0Source), from the PDG code of the mother of the V0
  template <typename T>
  int v0Source(T const& v0)
  {
    if (v0.isPhysicalPrimary()) {
      return kSourcePrimary;
    }
    switch (std::abs(v0.pdgCodeMother())) {
      case PDG_t::kXiMinus:
        return kSourceXiCharged;
      case o2::constants::physics::Pdg::kXi0:
        return kSourceXi0;
      case PDG_t::kOmegaMinus:
        return kSourceOmega;
    }
    return kSourceOtherSecondary;
  }

  // V0 content of a selected jet: number of V0 constituents (all, and passing the selection of each species). In MC also
  // the number of real / fake constituents, the true selected constituents, the MC origin of the constituents and the
  // fraction of the (raw) jet pT carried by fake V0s, which are clustered by the jet finder as the real ones
  template <bool isMC, typename TV0Table, typename TJet>
  void fillJetV0Content(TJet const& jet, float ptSub, float weight = 1.)
  {
    int nV0s = 0;
    int nReal = 0;
    int nFake = 0;
    std::array<int, 3> nSelected{0, 0, 0};
    std::array<int, 3> nTrueSelected{0, 0, 0};
    float fakePt = 0.;
    float fakeSelectedPt = 0.;
    std::array<float, kNV0Sources> sourcePt{};
    for (const auto& v0 : jet.template candidates_as<TV0Table>()) {
      nV0s++;
      bool passedAnySelection = false;
      for (const auto& species : {kK0s, kLambda, kAntiLambda}) {
        if (!passedV0Selection(v0, species)) {
          continue;
        }
        passedAnySelection = true;
        nSelected[species]++;
        if constexpr (isMC) {
          if (isTrueV0(v0, species)) {
            nTrueSelected[species]++;
          }
          registryMC.fill(HIST("jetV0s/SelectedV0ConstituentOrigin"), species, v0.pt(), v0Origin(v0), weight);
        }
      }
      if constexpr (isMC) {
        const int origin = v0Origin(v0);
        registryMC.fill(HIST("jetV0s/V0ConstituentOrigin"), ptSub, v0.pt(), origin, weight);
        if (origin <= kOriginAntiLambda) {
          nReal++;
          const int source = v0Source(v0);
          sourcePt[source] += v0.pt();
          registryMC.fill(HIST("jetV0s/V0ConstituentSource"), ptSub, v0.pt(), source, weight);
        } else {
          nFake++;
          fakePt += v0.pt();
          if (passedAnySelection) {
            fakeSelectedPt += v0.pt();
          }
        }
      }
    }

    if constexpr (isMC) {
      registryMC.fill(HIST("jetV0s/nV0Constituents"), ptSub, nV0s, weight);
      registryMC.fill(HIST("jetV0s/nK0sSelected"), ptSub, nSelected[kK0s], weight);
      registryMC.fill(HIST("jetV0s/nLambdaSelected"), ptSub, nSelected[kLambda], weight);
      registryMC.fill(HIST("jetV0s/nAntiLambdaSelected"), ptSub, nSelected[kAntiLambda], weight);
      registryMC.fill(HIST("jetV0s/nRealV0Constituents"), ptSub, nReal, weight);
      registryMC.fill(HIST("jetV0s/nFakeV0Constituents"), ptSub, nFake, weight);
      registryMC.fill(HIST("jetV0s/nTrueK0sSelected"), ptSub, nTrueSelected[kK0s], weight);
      registryMC.fill(HIST("jetV0s/nTrueLambdaSelected"), ptSub, nTrueSelected[kLambda], weight);
      registryMC.fill(HIST("jetV0s/nTrueAntiLambdaSelected"), ptSub, nTrueSelected[kAntiLambda], weight);
      registryMC.fill(HIST("jetV0s/fakeV0PtFraction"), ptSub, fakePt / jet.pt(), weight);
      registryMC.fill(HIST("jetV0s/fakeSelectedV0PtFraction"), ptSub, fakeSelectedPt / jet.pt(), weight);
      registryMC.fill(HIST("jetV0s/fakeV0PtFractionProfile"), ptSub, fakePt / jet.pt(), weight);
      registryMC.fill(HIST("jetV0s/fakeSelectedV0PtFractionProfile"), ptSub, fakeSelectedPt / jet.pt(), weight);
      for (int source = 0; source < kNV0Sources; source++) {
        registryMC.fill(HIST("jetV0s/V0SourcePtFraction"), ptSub, source, sourcePt[source] / jet.pt(), weight);
        registryMC.fill(HIST("jetV0s/V0SourcePtFractionProfile"), ptSub, source, sourcePt[source] / jet.pt(), weight);
      }
    } else {
      registryData.fill(HIST("jetV0s/nV0Constituents"), ptSub, nV0s);
      registryData.fill(HIST("jetV0s/nK0sSelected"), ptSub, nSelected[kK0s]);
      registryData.fill(HIST("jetV0s/nLambdaSelected"), ptSub, nSelected[kLambda]);
      registryData.fill(HIST("jetV0s/nAntiLambdaSelected"), ptSub, nSelected[kAntiLambda]);
    }
  }

  // Reconstructed true V0s passing the selections: primaries in iHist, primaries + secondaries in iHistIncl (if >= 0)
  template <int iHist, int iHistIncl, typename T>
  void fillMCRecoV0(T const& v0, float multiplicity, float weight)
  {
    for (const auto& species : {kK0s, kLambda, kAntiLambda}) {
      if (!isTrueV0(v0, species) || !passedV0Selection(v0, species)) {
        continue;
      }
      if constexpr (iHistIncl >= 0) {
        fillMCHistogram<iHistIncl>(species, multiplicity, v0.pt(), weight);
      }
      if (v0.isPhysicalPrimary()) {
        fillMCHistogram<iHist>(species, multiplicity, v0.pt(), weight);
      }
    }
  }

  // Generated physical primary V0s from the PWGJE reduced MC particle table
  template <int iHist, typename T>
  void fillMCGenV0(T const& v0, float multiplicity, float weight)
  {
    const int species = genV0Species(v0);
    if (species < 0 || !v0.isPhysicalPrimary()) {
      return;
    }
    fillMCHistogram<iHist>(static_cast<V0Species>(species), multiplicity, v0.pt(), weight);
  }

  // MC detector level: reconstructed V0s (efficiency numerators, feed-down) and generated V0s in the reconstructed
  // events (efficiency denominators), in jets and in the UE of the reconstructed jets
  void processMCReco(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, aod::JetMcCollisions const&, MCDV0JetsWithConstituents const& v0jets,
                     aod::CandidatesV0MCD const& v0s, aod::JetTracks const& tracks, aod::JetParticles const& jetParticles)
  {
    if (!jcoll.has_mcCollision() || !isSelectedCollision(jcoll, jetFinderConfigDet)) {
      return;
    }
    const auto mcCollision = jcoll.mcCollision_as<aod::JetMcCollisions>();
    const float multiplicity = mcCollision.centFT0M();
    const float weight = mcCollision.weight();
    registryMC.fill(HIST("number_of_events_vsmultiplicity_rec"), multiplicity, weight);

    // Generated particles of the MC collision associated to the reconstructed collision
    const auto genV0s = jetParticles.sliceBy(particlesPerMcCollision, mcCollision.globalIndex());

    const auto jetFinderV0s = buildDetectorLevelInput(v0s, tracks);

    if (fillV0QAHistograms) {
      for (const auto& v0 : v0s) {
        fillV0QA(v0, jetFinderV0s.contains(v0.globalIndex()));
      }
    }

    // Inclusive efficiency: reconstructed and generated V0s with |eta_gen| < inclusiveEtaMax, in all the reconstructed events
    for (const auto& v0 : v0s) {
      if (std::abs(RecoDecay::eta(std::array{v0.pxMC(), v0.pyMC(), v0.pzMC()})) < inclusiveEtaMax) {
        fillMCRecoV0<kRecFullEvent, -1>(v0, multiplicity, weight);
      }
    }
    for (const auto& genV0 : genV0s) {
      if (std::abs(genV0.eta()) < inclusiveEtaMax) {
        fillMCGenV0<kGenRecoEvent>(genV0, multiplicity, weight);
      }
    }

    const float rho = estimateRho(v0jets);
    if (rho < 0.) {
      return;
    }
    registryMC.fill(HIST("rhoPerpCone_rec"), multiplicity, rho, weight);

    std::unordered_set<int> v0sInJets;
    std::vector<std::array<float, 2>> jetAxes;
    collectJetsAboveThreshold(v0jets, rho, v0sInJets, jetAxes);

    bool hasSelectedJet = false;
    for (const auto& jet : v0jets) {
      if (!isSelectedJet(jet, rho)) {
        continue;
      }
      hasSelectedJet = true;
      registryMC.fill(HIST("number_of_jets_vsmultiplicity_rec"), multiplicity, weight);
      fillJetV0Content<true, aod::CandidatesV0MCD>(jet, subtractedPt(jet, rho), weight);

      // Reconstructed V0s, in-jet region: constituents of the reconstructed jet (same definition as in data)
      for (const auto& v0 : jet.candidates_as<aod::CandidatesV0MCD>()) {
        fillMCRecoV0<kRecJet, kRecJetIncl>(v0, multiplicity, weight);
      }
      // Reconstructed V0s, UE region: same definition as in data
      for (const auto& v0 : v0s) {
        if (!jetFinderV0s.contains(v0.globalIndex()) || v0sInJets.contains(v0.globalIndex())) {
          continue;
        }
        if (isInUE(v0.eta(), v0.phi(), jet)) {
          fillMCRecoV0<kRecUE, kRecUEIncl>(v0, multiplicity, weight);
        }
      }

      // Generated V0s placed with respect to the reconstructed jet axis, using the generated momentum. The generated V0s
      // are not jet constituents, so "in jet" means within R from the jet axis. In the UE, the generated V0s within R
      // from the axis of another jet above threshold are excluded (as the reconstructed V0s constituents of another jet)
      for (const auto& genV0 : genV0s) {
        if (jetutilities::deltaR(genV0.eta(), genV0.phi(), jet.eta(), jet.phi()) < rJet) {
          fillMCGenV0<kGenRecoJetJet>(genV0, multiplicity, weight);
        } else if (isInUE(genV0.eta(), genV0.phi(), jet) && !isInAnyJetCone(genV0.eta(), genV0.phi(), jetAxes)) {
          fillMCGenV0<kGenRecoJetUE>(genV0, multiplicity, weight);
        }
      }
    }

    // Inclusive efficiency in the events with at least one selected jet
    if (hasSelectedJet) {
      registryMC.fill(HIST("number_of_events_wjet_vsmultiplicity_rec"), multiplicity, weight);
      for (const auto& v0 : v0s) {
        if (std::abs(RecoDecay::eta(std::array{v0.pxMC(), v0.pyMC(), v0.pzMC()})) < inclusiveEtaMax) {
          fillMCRecoV0<kRecEventWithJet, -1>(v0, multiplicity, weight);
        }
      }
      for (const auto& genV0 : genV0s) {
        if (std::abs(genV0.eta()) < inclusiveEtaMax) {
          fillMCGenV0<kGenEventWithJet>(genV0, multiplicity, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(LfJetSpectra, processMCReco, "V0s in jets and in the UE, MC detector level (efficiency and feed-down)", false);

  // MC particle level: generated V0s in particle-level jets and in their UE (defined as in data), and in the full event
  void processMCGen(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision, MCPV0JetsWithConstituents const& v0jets,
                    aod::CandidatesV0MCP const& v0s, aod::JetParticles const& particles)
  {
    if (!jetderiveddatautilities::selectCollision(mcCollision, jetFinderConfigGen.eventSelectionBits, jetFinderConfigGen.skipMBGapEvents, jetFinderConfigGen.applyRCTSelections)) {
      return;
    }
    const float multiplicity = mcCollision.centFT0M();
    const float weight = mcCollision.weight();
    registryMC.fill(HIST("number_of_events_vsmultiplicity_gen"), multiplicity, weight);

    const auto genV0s = particles.sliceBy(particlesPerMcCollision, mcCollision.globalIndex());
    for (const auto& genV0 : genV0s) {
      if (std::abs(genV0.eta()) < inclusiveEtaMax) {
        fillMCGenV0<kGenFullEvent>(genV0, multiplicity, weight);
      }
    }

    const auto jetFinderV0s = buildParticleLevelInput(v0s, particles);
    const float rho = estimateRho(v0jets);
    if (rho < 0.) {
      return;
    }
    registryMC.fill(HIST("rhoPerpCone_gen"), multiplicity, rho, weight);

    std::unordered_set<int> v0sInJets;
    std::vector<std::array<float, 2>> jetAxes;
    collectJetsAboveThreshold(v0jets, rho, v0sInJets, jetAxes);

    for (const auto& jet : v0jets) {
      if (!isSelectedJet(jet, rho)) {
        continue;
      }
      registryMC.fill(HIST("number_of_jets_vsmultiplicity_gen"), multiplicity, weight);

      // Generated V0s, in-jet region: constituents of the particle-level jet (anti-kT clustering of the V0s as single
      // pseudojets, daughters removed from the particle list, as at detector level)
      for (const auto& v0 : jet.candidates_as<aod::CandidatesV0MCP>()) {
        fillMCGenV0<kGenJet>(v0, multiplicity, weight);
      }
      // Generated V0s, UE region: accepted by the particle-level jet finder and not constituents of another jet above threshold
      for (const auto& v0 : v0s) {
        if (!jetFinderV0s.contains(v0.globalIndex()) || v0sInJets.contains(v0.globalIndex())) {
          continue;
        }
        if (isInUE(v0.eta(), v0.phi(), jet)) {
          fillMCGenV0<kGenUE>(v0, multiplicity, weight);
        }
      }
    }
  }
  PROCESS_SWITCH(LfJetSpectra, processMCGen, "V0s in jets and in the UE, MC particle level", false);

  // ---------------------------------------------------------------------------------------------------------------
  // Jet pT response (requires jet-matching-mc-v0-ch)

  // Fraction of the (raw) pT of a particle-level jet carried by Xi+- and Omega+- (charged primaries, clustered as single
  // particles at particle level, while at detector level their momentum is carried by the secondary V0 and the bachelor)
  template <typename TJet>
  float multiStrangePtFraction(TJet const& mcpJet)
  {
    float pt = 0.;
    for (const auto& particle : mcpJet.template tracks_as<aod::JetParticles>()) {
      const int pdg = std::abs(particle.pdgCode());
      if (pdg == PDG_t::kXiMinus || pdg == PDG_t::kOmegaMinus) {
        pt += particle.pt();
      }
    }
    return pt / mcpJet.pt();
  }

  // Detector-level jets of the reconstructed collision and particle-level jets matched to them. The background is
  // estimated at each level from the jet finder input of that level (particle level: MC collision associated to the
  // reconstructed collision), as in processMCReco and processMCGen. Detector-level side: response and matching purity
  void processMCResponse(soa::Filtered<aod::JetCollisionsMCD>::iterator const& jcoll, aod::JetMcCollisions const&, MatchedMCDV0Jets const& mcdJets,
                         aod::CandidatesV0MCD const& v0s, aod::JetTracks const& tracks,
                         MatchedMCPV0JetsWithConstituents const& mcpJets, aod::CandidatesV0MCP const& v0sMCP, aod::JetParticles const& particles)
  {
    if (!jcoll.has_mcCollision() || !isSelectedCollision(jcoll, jetFinderConfigDet)) {
      return;
    }
    const auto mcCollision = jcoll.mcCollision_as<aod::JetMcCollisions>();
    if (std::abs(mcCollision.posZ()) >= vertexZCut ||
        !jetderiveddatautilities::selectCollision(mcCollision, jetFinderConfigGen.eventSelectionBits, jetFinderConfigGen.skipMBGapEvents, jetFinderConfigGen.applyRCTSelections)) {
      return;
    }
    const float weight = mcCollision.weight();

    buildDetectorLevelInput(v0s, tracks);
    const float rhoDet = estimateRho(mcdJets);
    if (rhoDet < 0.) {
      return;
    }
    const auto mcpJetsOfMcCollision = mcpJets.sliceBy(mcpJetsPerMcCollision, jcoll.mcCollisionId());
    buildParticleLevelInput(v0sMCP.sliceBy(v0MCPPerMcCollision, jcoll.mcCollisionId()), particles.sliceBy(particlesPerMcCollision, jcoll.mcCollisionId()));
    const float rhoGen = estimateRho(mcpJetsOfMcCollision); // < 0 only if there are no particle-level jets, i.e. no matches

    for (const auto& mcdJet : mcdJets) {
      if (!isInJetAcceptance(mcdJet)) {
        continue;
      }
      const float ptSubDet = subtractedPt(mcdJet, rhoDet);
      registryMC.fill(HIST("jetResponse/detJets"), ptSubDet, weight);
      for (const auto& mcpJet : mcdJet.template matchedJetGeo_as<MatchedMCPV0JetsWithConstituents>()) {
        if (!isInJetAcceptance(mcpJet)) {
          continue;
        }
        const float ptSubGen = subtractedPt(mcpJet, rhoGen);
        const float residualRaw = (mcdJet.pt() - mcpJet.pt()) / mcpJet.pt();
        registryMC.fill(HIST("jetResponse/detJetsMatched"), ptSubDet, weight);
        registryMC.fill(HIST("jetResponse/ptSubDetVsGen"), ptSubDet, ptSubGen, weight);
        registryMC.fill(HIST("jetResponse/ptRawDetVsGen"), mcdJet.pt(), mcpJet.pt(), weight);
        if (ptSubGen > 0.) {
          registryMC.fill(HIST("jetResponse/residualSub"), ptSubGen, (ptSubDet - ptSubGen) / ptSubGen, weight);
        }
        registryMC.fill(HIST("jetResponse/residualRaw"), mcpJet.pt(), residualRaw, weight);
        registryMC.fill(HIST("jetResponse/residualRawVsMultiStrangeFraction"), mcpJet.pt(), multiStrangePtFraction(mcpJet), residualRaw, weight);
      }
    }
  }
  PROCESS_SWITCH(LfJetSpectra, processMCResponse, "jet pT response (detector-level side), requires jet-matching-mc-v0-ch", false);

  // Particle-level side of the jet pT response: matching efficiency of the particle-level jets in acceptance (matched to a
  // detector-level jet in acceptance, in a reconstructed collision passing the event selection of processMCResponse)
  void processMCResponseGen(soa::Filtered<aod::JetMcCollisions>::iterator const& mcCollision, MatchedMCPV0JetsWithConstituents const& mcpJets,
                            aod::CandidatesV0MCP const& v0s, aod::JetParticles const& particles,
                            soa::SmallGroups<aod::JetCollisionsMCD> const& collisions, MatchedMCDV0Jets const&)
  {
    if (!jetderiveddatautilities::selectCollision(mcCollision, jetFinderConfigGen.eventSelectionBits, jetFinderConfigGen.skipMBGapEvents, jetFinderConfigGen.applyRCTSelections)) {
      return;
    }
    const float weight = mcCollision.weight();

    // SmallGroups contains the reconstructed collisions linked to this JetMcCollision.
    std::unordered_set<int> selectedCollisions;
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) < vertexZCut && isSelectedCollision(collision, jetFinderConfigDet)) {
        selectedCollisions.insert(collision.globalIndex());
      }
    }

    buildParticleLevelInput(v0s, particles);
    const float rho = estimateRho(mcpJets);
    if (rho < 0.) {
      return;
    }
    for (const auto& mcpJet : mcpJets) {
      if (!isInJetAcceptance(mcpJet)) {
        continue;
      }
      const float ptSub = subtractedPt(mcpJet, rho);
      registryMC.fill(HIST("jetResponse/genJets"), ptSub, weight);
      for (const auto& mcdJet : mcpJet.template matchedJetGeo_as<MatchedMCDV0Jets>()) {
        if (isInJetAcceptance(mcdJet) && selectedCollisions.contains(mcdJet.collisionId())) {
          registryMC.fill(HIST("jetResponse/genJetsMatched"), ptSub, weight);
          break;
        }
      }
    }
  }
  PROCESS_SWITCH(LfJetSpectra, processMCResponseGen, "jet pT response (particle-level side: matching efficiency), requires jet-matching-mc-v0-ch", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<LfJetSpectra>(cfgc)};
}
