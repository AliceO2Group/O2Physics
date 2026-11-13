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

/// \file promptphotonQA.cxx
/// \brief Generator and reconstructed photon QA with shared MC/data selections.
/// \author Adrian Fereydon Nassirpour

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/JetReducedData.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include <TH1.h>
#include <TH2.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>

struct PromptphotonQA {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry extendedHistos{"extendedHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  // Separate the large energy-stack collection from the other QA registries.
  HistogramRegistry clusterEnergyHistos{"clusterEnergyHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry dataClusterEnergyHistos{"dataClusterEnergyHistos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Service<o2::framework::O2DatabasePDG> pdg{};

  // Shared counts, generator status codes and numerical validation bounds.
  static constexpr int NSamples = 2;
  static constexpr int NNeutralBranches = 3;
  static constexpr int NAnalysisStages = 7;
  static constexpr int NClusterSelections = 4;
  static constexpr int NEventCuts = 7;
  static constexpr int NClusterCuts = 9;
  static constexpr int NTruthGroups = 3;
  static constexpr int NTruthLabelStatuses = 4;
  static constexpr int NPhotonGroups = 2;
  static constexpr int NAncestryClasses = 5;
  static constexpr int HadronCodeLowerBound = 100;
  static constexpr int NucleusCodeLowerBound = 1000000000;
  static constexpr int NTrackingTruthSets = 5;
  static constexpr int NSignalSamples = 3;
  static constexpr int NTrackingMatches = 3;
  static constexpr int FinalGeneratorStatusEnd = 81;
  static constexpr int HardProcessStatus = -23;
  static constexpr float MaxContributorFraction = 1.0001f;
  static constexpr double MaxContributorFractionSum = 1.001;
  static constexpr double MaxExtrapolationAbsEta = 10.;
  static constexpr double MinMatchingExtrapolationPhi = -10.;

  // basic enum for cateogrization of prompt photons
  enum class Origin : int { Unclassified = 0,
                            Direct = 1,
                            Fragmentation = 2,
                            Decay = 3,
                            Excluded = 4,
                            Mixed = 5 };
  struct OriginalResult {
    Origin origin = Origin::Unclassified;
    bool leptonBranch = false;
    std::vector<int64_t> directHardIndices;
    std::vector<int64_t> fragmentationHardIndices;
  };
  // Configuration and framework filters
  Configurable<double> cfgMaxDCArToPVcut{"cfgMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};
  Configurable<double> cfgMaxDCAzToPVcut{"cfgMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgConnectedToPV{"cfgConnectedToPV", true, "PV contributor track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<double> cfgnFindableTPCClusters{"cfgnFindableTPCClusters", 50, "nFindable TPC Clusters"};
  Configurable<double> cfgnTPCCrossedRows{"cfgnTPCCrossedRows", 70, "nCrossed TPC Rows"};
  Configurable<double> cfgnRowsOverFindable{"cfgnRowsOverFindable", 1.2, "Maximum TPC crossed rows / findable clusters (old selection)"};
  Configurable<double> cfgTpcChi2{"cfgTpcChi2", 4.0, "nTPC Chi2 per Cluster"};
  Configurable<double> cfgItsChi2{"cfgItsChi2", 36.0, "nITS Chi2 per Cluster"};
  Configurable<int> cfgClusterDefinition{"cfgClusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default, 0 = kV1Default"};
  Configurable<float> cfgMinTime{"cfgMinTime", -30., "Minimum cluster time for time cut"};
  Configurable<float> cfgMaxTime{"cfgMaxTime", +35., "Maximum cluster time for time cut"};
  Configurable<float> cfgMinClusterEnergy{"cfgMinClusterEnergy", 0.7f, "Minimal cluster energy"};
  Configurable<int> cfgMinNCells{"cfgMinNCells", 2, "Minimal amount of cells per cluster"};
  Configurable<int> cfgMaxNLM{"cfgMaxNLM", 2, "Maximal amount of local Maxima per cluster"};
  Configurable<bool> cfgExoticContribution{"cfgExoticContribution", false, "Exotic cluster in the data"};
  Configurable<double> cfgtrkMinPt{"cfgtrkMinPt", 0.15, "set track min pT"};
  Configurable<double> cfgtrkMaxEta{"cfgtrkMaxEta", 0.6, "set track max Eta"};
  Configurable<float> cfgVtxCut{"cfgVtxCut", 10.0, "V_z cut selection"};
  Configurable<float> cfgLowM02{"cfgLowM02", 0.1, "Lower-bound M02 cut"};
  Configurable<float> cfgHighM02{"cfgHighM02", 0.3, "Higher-bound M02 cut"};
  Configurable<float> cfgLowClusterE{"cfgLowClusterE", 0.5, "Higher-bound Cluster E cut"};
  Configurable<float> cfgHighClusterE{"cfgHighClusterE", 500, "Lower-bound Cluster E cut"};
  Configurable<bool> cfgEmcTrigger{"cfgEmcTrigger", true, "Require EMC readout for event"};
  Configurable<bool> cfgPtClusterCut{"cfgPtClusterCut", true, "Performs Pt-dependent cluster-track matching"};
  Configurable<std::string> cfgTrackFilter{"cfgTrackFilter", "globalTracks", "set track selections"};
  Configurable<bool> cfgJETracks{"cfgJETracks", false, "Use JE track-selection bits instead of original track quality cuts"};
  Configurable<bool> cfgReqRecPsRec{"cfgReqRecPsRec", true, "Apply EMCal/DCal fiducial acceptance to reconstructed clusters in both modes"};
  Configurable<float> cfgEMClowPSphi{"cfgEMClowPSphi", static_cast<float>(80.f * o2::constants::math::Deg2Rad), "lower limit of the EMC acceptance if Rec PS is required"};
  Configurable<float> cfgEMChighPSphi{"cfgEMChighPSphi", static_cast<float>(187.f * o2::constants::math::Deg2Rad), "higher limit of the EMC acceptance if Rec PS is required"};
  Configurable<float> cfgEMChighPSeta{"cfgEMChighPSeta", 0.7, "symmetric eta cut if Rec PS is required"};
  Configurable<float> cfgDClowPSphi{"cfgDClowPSphi", static_cast<float>(260.f * o2::constants::math::Deg2Rad), "lower limit of the DCal acceptance if Rec PS is required"};
  Configurable<float> cfgDChighPSphi{"cfgDChighPSphi", static_cast<float>(327.f * o2::constants::math::Deg2Rad), "higher limit of the DCal acceptance if Rec PS is required"};
  Configurable<float> cfgDCGapEndPhi{"cfgDCGapEndPhi", static_cast<float>(320. * o2::constants::math::Deg2Rad), "End of full-size DCal modules (rad)"};
  Configurable<float> cfgDCGapEta{"cfgDCGapEta", 0.22f, "Minimum |eta| in full-size DCal modules"};
  Configurable<float> cfgMatchMinEOverP{"cfgMatchMinEOverP", 1.7f, "Minimum cluster E / scalar sum of matched track momenta for recovery study"};
  Configurable<std::string> cfgTriggerMasks{"cfgTriggerMasks", "", "possible JE Trigger masks: fJetChLowPt,fJetChHighPt,fTrackLowPt,fTrackHighPt,fJetD0ChLowPt,fJetD0ChHighPt,fJetLcChLowPt,fJetLcChHighPt,fEMCALReadout,fJetFullHighPt,fJetFullLowPt,fJetNeutralHighPt,fJetNeutralLowPt,fGammaVeryHighPtEMCAL,fGammaVeryHighPtDCAL,fGammaHighPtEMCAL,fGammaHighPtDCAL,fGammaLowPtEMCAL,fGammaLowPtDCAL,fGammaVeryLowPtEMCAL,fGammaVeryLowPtDCAL"};
  Configurable<float> cfgMinPromptFraction{"cfgMinPromptFraction", 0.5f, "Exclusive minimum single-contributor energy fraction for truth dominance; also used for prompt tagging"};
  Configurable<bool> cfgRecoverHighEOverSumPt{"cfgRecoverHighEOverSumPt", false, "Allow matched cluster with E/sumPt >= threshold to pass charged veto"};
  Configurable<float> cfgEOverSumPtMin{"cfgEOverSumPtMin", 1.7f, "E/sumPt threshold for optional recovery"};
  // Counterpart means any reconstructed collision carrying this MC label, before REC event cuts.
  Configurable<bool> cfgGenReqRec{"cfgGenReqRec", false, "Require a reconstructed counterpart for GEN and selected-event efficiency denominators"};
  Configurable<bool> cfgRequireEventTrack{"cfgRequireEventTrack", false, "Require a selected track in both MC and data (old data-only cut)"};

  Configurable<float> cfgIsoRadius{"cfgIsoRadius", 0.4f, "Charged isolation cone radius"};
  Configurable<float> cfgIsoMaxSumPt{"cfgIsoMaxSumPt", 1.5f, "Upper exclusive raw cone sum pT (GeV/c)"};
  Configurable<float> cfgIsoMinEOverP{"cfgIsoMinEOverP", 1.7f, "Minimum candidate E / scalar sum of charged momenta in cone"};
  Configurable<bool> cfgIsoExcludeMatchedTracks{"cfgIsoExcludeMatchedTracks", false, "Exclude tracks matched to the candidate cluster from reconstructed cone sums"};
  Configurable<bool> cfgIsoFullCone{"cfgIsoFullCone", true, "Require the cone inside constituent eta acceptance"};

  // Common baseline for both reconstructed samples. Keep acceptance and M02
  // downstream so the before/after shower-shape distributions stay distinct.
  Filter clusterDefinitionSelection =
    (aod::jcluster::definition == cfgClusterDefinition) &&
    (aod::jcluster::time >= cfgMinTime) &&
    (aod::jcluster::time <= cfgMaxTime) &&
    (nabs(aod::jcluster::time) <= std::numeric_limits<float>::max()) &&
    (aod::jcluster::energy > cfgMinClusterEnergy) &&
    (aod::jcluster::energy <= std::numeric_limits<float>::max()) &&
    (aod::jcluster::nCells >= cfgMinNCells) &&
    (aod::jcluster::nlm <= cfgMaxNLM) &&
    (aod::jcluster::isExotic == cfgExoticContribution);

  // AOD merged tables. We of course require the correct task for all of these tables.
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using RecoTracks = soa::Join<aod::JTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using DataCollisions = soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::EvSels, aod::JEMCCollisionLbs>;
  using MCCollisions = soa::Join<aod::JCollisions, aod::JCollisionPIs, aod::EvSels, aod::JEMCCollisionLbs, aod::JMcCollisionLbs>;
  using DataClusters = soa::Join<aod::JClusters, aod::JClusterTracks>;
  using MCClusters = soa::Join<aod::JClusters, aod::JClusterTracks, aod::JMcClusterLbs>;
  using FilteredDataClusters = soa::Filtered<DataClusters>;
  using FilteredMCClusters = soa::Filtered<MCClusters>;
  using BcCandidates = soa::Join<aod::JBCs, aod::JBCPIs>;
  // init presciles and filters
  PresliceUnsorted<aod::JEMCTracks> emcTrackPerTrack = aod::jemctrack::trackId;
  Preslice<aod::JClusters> clustersPerCollision = aod::jcluster::collisionId;
  Preslice<aod::JTracks> tracksPerCollision = aod::jtrack::collisionId;
  int trackFilter = -1;
  std::vector<int> triggerMaskBits;

  // WRAPPERS to fill histograms WITH and WITHOUT the pThat weight from the MC. QA1 is for 1D histos, QA2 is for 2d histos.
  // by having these wrappers we ensure that the weighted and unweighted histos are always filed simultaneously and identically.
  struct QA1 {
    std::shared_ptr<TH1> weighted, unweighted;
    void fill(double x, double weight)
    {
      weighted->Fill(x, weight);
      if (unweighted) {
        unweighted->Fill(x, 1.0);
      }
    }
    void setBinLabel(int bin, const char* label)
    {
      weighted->GetXaxis()->SetBinLabel(bin, label);
      if (unweighted) {
        unweighted->GetXaxis()->SetBinLabel(bin, label);
      }
    }
  };
  struct QA2 {
    std::shared_ptr<TH2> weighted, unweighted;
    void fill(double x, double y, double weight)
    {
      weighted->Fill(x, y, weight);
      if (unweighted) {
        unweighted->Fill(x, y, 1.0);
      }
    }
  };

  // Histogram storage and selection stages, as well as several helper structs
  enum PerformanceStage { BeforeM02,
                          AfterM02,
                          Neutral,
                          Isolated,
                          IsolatedEOverP,
                          StrictNeutral,
                          MatchRecovery,
                          MatchRecoveryIsolated,
                          nPerformanceStages };
  inline static const std::array<std::string, nPerformanceStages> performanceStages{
    "beforeM02", "afterM02", "neutral", "isolated", "isolatedEOverP", "strictNeutral", "matchedEOverPRecovery", "matchedEOverPRecoveryIsolated"};

  // Exclusive truth categories for energy stacks. Inclusive prompt is an overlay,
  // since direct, fragmentation and mixed prompt are already separate stack entries.
  enum class ClusterOrigin : int { Direct,
                                   Fragmentation,
                                   MixedPrompt,
                                   Decay,
                                   Charged,
                                   Other,
                                   NoDominant,
                                   UnknownTruth,
                                   Count };
  static constexpr int NClusterOrigins = static_cast<int>(ClusterOrigin::Count);
  inline static const std::array<std::string, NClusterOrigins> clusterOrigins{
    "direct", "fragmentation", "mixedPrompt", "decay", "charged",
    "other", "noDominant", "unknownTruth"};
  inline static const std::array<std::string, 7> analysisEnergyStages{
    "beforeM02", "afterM02", "validMatching", "neutral", "coneEligible", "isolated", "isolatedEOverP"};

  struct ClusterEnergyQA {
    std::shared_ptr<QA1> total, prompt;
    std::array<std::shared_ptr<QA1>, NClusterOrigins> exclusive;
  };

  struct AnalysisCutQA {
    std::shared_ptr<QA1> candidates, goodSignal;
    std::shared_ptr<QA2> candidatesVsPt, goodSignalVsPt;
    std::array<ClusterEnergyQA, 7> energy;
  };

  struct ConeParticle {
    double pt = 0., eta = 0., phi = 0.;
    int64_t trackId = -1; // Reconstructed track identity; unused for generator constituents.
  };

  struct IsolationQA {
    std::shared_ptr<QA2> sumPt, eOverP;
    std::shared_ptr<QA1> eligible, isolated, eOverPPass, combined;
  };

  struct TrackKinematicsQA {
    std::shared_ptr<QA1> pt, eta, phi;
    std::shared_ptr<QA2> etaPhi;
  };
  struct TrackingQA {
    TrackKinematicsQA global, extrapolated, inCalo;
    std::shared_ptr<QA1> surfaceStatus;
    std::shared_ptr<QA2> etaShift, phiShift;
  };

  struct SignalValidationQA {
    std::shared_ptr<QA1> allCandidates, validTruth, anyPrompt, goodSignal, badSignal;
    std::shared_ptr<QA1> noPrompt, promptBelowThreshold, unknownTruth;
    std::shared_ptr<QA2> maxPromptFraction, totalPromptFraction;
  };

  struct ClusterQA {
    std::shared_ptr<QA1> energy, pt, eta, phi, m02, m20, cells, time, nlm;
    std::shared_ptr<QA2> m02Energy, m20Energy, etaPhi;
  };

  struct QAState {
    std::array<std::array<AnalysisCutQA, 3>, 2> analysisCuts;
    std::array<std::array<ClusterEnergyQA, 4>, 2> clusterEnergy;
    // Neutral, cone-eligible neutral, and isolated neutral candidate samples.
    std::array<SignalValidationQA, 3> signalValidation;
    std::array<TrackingQA, 2> tracks;
    // All tracking performance histograms use truth kinematics and unique truth IDs.
    std::array<TrackKinematicsQA, 5> trackingTruth;
    std::shared_ptr<QA1> trackingLabelStatus;
    IsolationQA genIsolation;
    std::array<IsolationQA, 2> recoIsolation;
    std::shared_ptr<QA2> genParticleMap, allPrimaryMap, allPhotonMap;
    std::array<std::shared_ptr<QA2>, 2> allClusterMap;
    std::array<std::shared_ptr<QA2>, 2> genPhotonMap, trackMap, matchedEOverP;
    std::array<std::shared_ptr<QA1>, 2> nLinked, nMatched, nLinkedValid, invalidMatchingPt, strictPt, recoveryPt, recoveredIsoPt;

    // Efficiency sets are reset once per dataframe, not per reconstructed collision.
    // This also deduplicates photons across split reconstructed collisions.
    std::array<std::array<std::set<int64_t>, nPerformanceStages>, 2> efficiencyMatched;
    std::array<std::shared_ptr<QA1>, 2> efficiencyAll, efficiencyAccepted, coneAll, coneAccepted;
    std::array<std::array<std::shared_ptr<QA1>, nPerformanceStages>, 2> efficiencyNumerator, purityNumerator;
    std::array<std::shared_ptr<QA1>, nPerformanceStages> purityDenominator, validLabels;
    bool recoEventAccepted = false;

    // Stages: quality+acceptance before M02; after M02; charged-veto pass; veto fail.
    // Truth groups are orthogonal to these blind reconstructed selections.
    std::array<std::array<ClusterQA, 4>, 2> blindQA;
    std::array<std::array<ClusterQA, 4>, 3> truthQA; // prompt, direct, fragmentation (overlapping)
    std::array<std::shared_ptr<QA1>, 2> eventCuts, clusterCuts, sumPt, eOverSumPt;
    std::array<std::shared_ptr<QA2>, 2> matchResidual;
    std::shared_ptr<QA1> truthStatus, maxPromptFraction, totalPromptFraction;
    std::shared_ptr<QA2> promptFractionEnergy;
    std::shared_ptr<QA1> genPromptEta, genPromptPhi;
    std::array<std::shared_ptr<QA1>, 2> genPhi, genEta, genPt;
    std::shared_ptr<QA1> genEvents, ancestryStatus, directBeforeCalo, directInCalo;
    std::array<std::shared_ptr<QA1>, 2> selectedEventTruth, selectedConeTruth;
    std::array<std::shared_ptr<QA1>, nPerformanceStages> directFiducialTagged, fiducialCandidates;
  };
  QAState qaState;

  //**************************************************
  // HELPER STRUCTS DONE//////////////////////////////
  //**************************************************

  // NOW FOR HELPER FUNCTIONS FOR HISTOGRAM BOOKING
  ////////////////////////////////////////////////////
  // checks whether isolation is required
  static bool isolatedStage(int stage)
  {
    return stage == Isolated || stage == IsolatedEOverP || stage == MatchRecoveryIsolated;
  }

  // simple phi02pi wrapper
  static double phi02pi(double phi)
  {
    return RecoDecay::constrainAngle(phi);
  }

  // Checks whether origin is prompt for a cluster
  static bool isPrompt(Origin origin)
  {
    return origin == Origin::Direct || origin == Origin::Fragmentation || origin == Origin::Mixed;
  }

  // Each Fill makes the weighted and unit-weight fill together, with identical selection.
  static std::string unweightedTitle(std::string title)
  {
    for (auto const& word : {std::string("Weighted "), std::string("weighted ")}) {
      std::size_t pos = 0;
      while ((pos = title.find(word)) != std::string::npos) {
        title.erase(pos, word.size());
      }
    }
    auto pos = title.find("Sum weights");
    if (pos != std::string::npos) {
      title.replace(pos, 11, "Events");
    }
    return title;
  }

  // folder bookkeeping in the root file
  static bool isMCPath(std::string const& name)
  {
    return name.starts_with("REC/") || name.starts_with("GEN/");
  }
  static std::string weightPath(std::string const& name, bool weighted)
  {
    const auto slash = name.find('/');
    return name.substr(0, slash + 1) + (weighted ? "weighted/" : "unweighted/") + name.substr(slash + 1);
  }

  // NOW WE BOOKKEEP ALL OF THE HISTOGRAMS
  std::shared_ptr<QA1> h1(std::string const& name, std::string const& title, int n, double lo, double hi, bool extended = false)
  {
    auto& registry = extended ? extendedHistos : histos;
    auto result = std::make_shared<QA1>();
    const bool paired = isMCPath(name);
    const auto mainTitle = paired ? title : unweightedTitle(title);
    result->weighted = registry.add<TH1>(paired ? weightPath(name, true) : name, mainTitle.c_str(), kTH1D, {{n, lo, hi}}, true);
    if (paired) {
      result->unweighted = registry.add<TH1>(weightPath(name, false), unweightedTitle(title).c_str(), kTH1D, {{n, lo, hi}}, true);
    }
    return result;
  }

  std::shared_ptr<QA2> h2(std::string const& name, std::string const& title, AxisSpec const& x, AxisSpec const& y, bool extended = false)
  {
    auto& registry = extended ? extendedHistos : histos;
    auto result = std::make_shared<QA2>();
    const bool paired = isMCPath(name);
    result->weighted = registry.add<TH2>(paired ? weightPath(name, true) : name, title.c_str(), kTH2D, {x, y}, true);
    if (paired) {
      result->unweighted = registry.add<TH2>(weightPath(name, false), unweightedTitle(title).c_str(), kTH2D, {x, y}, true);
    }
    return result;
  }

  std::shared_ptr<QA2> bookMap(std::string const& path)
  {
    return h2(path, ";#eta;#varphi (rad)", {100, -1, 1}, {144, 0, o2::constants::math::TwoPI}, true);
  }

  IsolationQA bookIsolation(std::string const& path)
  {
    IsolationQA qa;
    qa.sumPt = h2(path + "/coneSumPt", ";Candidate pT (GeV/c);Charged cone sum pT (GeV/c)", {250, 0, 500}, {250, 0, 50}, true);
    qa.eOverP = h2(path + "/coneEOverP", ";Candidate pT (GeV/c);E / cone sum |p|", {250, 0, 500}, {200, 0, 10}, true);
    const std::array<std::pair<const char*, std::shared_ptr<QA1>*>, 4> entries{{{"/eligible", &qa.eligible}, {"/isolated", &qa.isolated}, {"/eOverPPass", &qa.eOverPPass}, {"/combined", &qa.combined}}};
    for (auto const& entry : entries) {
      *entry.second = h1(path + entry.first, ";Candidate pT (GeV/c);Candidates", 250, 0, 500, true);
    }
    return qa;
  }

  TrackKinematicsQA bookTrackKinematics(std::string const& path)
  {
    TrackKinematicsQA q;
    q.pt = h1(path + "/pt", ";pT (GeV/c);Tracks or particles", 250, 0, 50, true);
    q.eta = h1(path + "/eta", ";eta;Tracks or particles", 100, -1, 1, true);
    q.phi = h1(path + "/phi", ";phi (rad);Tracks or particles", 144, 0, o2::constants::math::TwoPI, true);
    q.etaPhi = bookMap(path + "/etaPhi");
    return q;
  }

  void fillTrackKinematics(TrackKinematicsQA& q, double pt, double eta, double phi, double w)
  {
    q.pt->fill(pt, w);
    q.eta->fill(eta, w);
    q.phi->fill(phi02pi(phi), w);
    q.etaPhi->fill(eta, phi02pi(phi), w);
  }

  void bookCluster(ClusterQA& q, std::string const& p)
  {
    q.etaPhi = bookMap(p + "/etaPhi");
    q.energy = h1(p + "/energy", ";E_{cluster} (GeV);Weighted clusters", 250, 0, 500);
    q.pt = h1(p + "/pt", ";E/cosh(eta) (GeV/c);Weighted clusters", 250, 0, 500);
    q.eta = h1(p + "/eta", ";eta;Weighted clusters", 80, -0.8, 0.8);
    q.phi = h1(p + "/phi", ";phi;Weighted clusters", 90, 0, o2::constants::math::TwoPI);
    q.m02 = h1(p + "/m02", ";M02;Weighted clusters", 100, 0, 2);
    q.m20 = h1(p + "/m20", ";M20;Weighted clusters", 100, 0, 2);
    q.cells = h1(p + "/nCells", ";N cells;Weighted clusters", 100, -0.5, 99.5);
    q.time = h1(p + "/time", ";Time (ns);Weighted clusters", 200, -100, 100);
    q.nlm = h1(p + "/nlm", ";N local maxima;Weighted clusters", 20, -0.5, 19.5);
    q.m02Energy = h2(p + "/m02VsEnergy", ";E (GeV);M02", {250, 0, 500}, {100, 0, 2});
    q.m20Energy = h2(p + "/m20VsEnergy", ";E (GeV);M20", {250, 0, 500}, {100, 0, 2});
  }

  template <typename C>
  void fillCluster(ClusterQA& q, C const& c, double w)
  {
    q.etaPhi->fill(c.eta(), phi02pi(c.phi()), w);
    q.energy->fill(c.energy(), w);
    q.pt->fill(c.energy() / std::cosh(c.eta()), w);
    q.eta->fill(c.eta(), w);
    q.phi->fill(phi02pi(c.phi()), w);
    q.m02->fill(c.m02(), w);
    q.m20->fill(c.m20(), w);
    q.cells->fill(c.nCells(), w);
    q.time->fill(c.time(), w);
    q.nlm->fill(c.nlm(), w);
    q.m02Energy->fill(c.energy(), c.m02(), w);
    q.m20Energy->fill(c.energy(), c.m20(), w);
  }

  //*********************************************
  // DONE WITH HISTO BOOKING. NOW FOR THE MAIN INIT
  //*********************************************
  std::shared_ptr<QA1> energyH1(std::string const& name, std::string const& title)
  {
    const bool paired = isMCPath(name);
    auto& registry = paired ? clusterEnergyHistos : dataClusterEnergyHistos;
    auto result = std::make_shared<QA1>();
    result->weighted = registry.add<TH1>(paired ? weightPath(name, true) : name,
                                         title.c_str(), kTH1D, {{250, 0, 500}}, true);
    if (paired) {
      result->unweighted = registry.add<TH1>(weightPath(name, false),
                                             unweightedTitle(title).c_str(), kTH1D, {{250, 0, 500}}, true);
    }
    return result;
  }

  ClusterEnergyQA bookClusterEnergy(std::string const& path, bool truth)
  {
    ClusterEnergyQA q;
    q.total = energyH1(path + "/total", ";Cluster E (GeV);Clusters");
    if (truth) {
      q.prompt = energyH1(path + "/prompt", ";Cluster E (GeV);Prompt-dominated clusters (inclusive overlay)");
      for (int i = 0; i < NClusterOrigins; ++i) {
        q.exclusive[i] = energyH1(path + "/exclusive/" + clusterOrigins[i], ";Cluster E (GeV);Clusters");
      }
    }
    return q;
  }

  void bookAnalysisCuts()
  {
    const std::array<std::string, 3> branches{"neutral", "strictNeutral", "matchedEOverPRecovery"};
    const std::array<const char*, 7> labels{
      "quality+energy+configured acceptance", "M02", "valid track matching",
      "neutral branch selection", "cone eligible", "charged isolation", "cone E/sum|p|"};
    for (int sample = 0; sample < NSamples; ++sample) {
      for (int branch = 0; branch < NNeutralBranches; ++branch) {
        auto& q = qaState.analysisCuts[sample][branch];
        const std::string base = (sample == 0 ? "REC/analysisCuts/" : "DATA/analysisCuts/") + branches[branch] + "/";
        for (int stage = 0; stage < NAnalysisStages; ++stage) {
          q.energy[stage] = bookClusterEnergy(base + "energy/" + analysisEnergyStages[stage], sample == 0);
        }
        q.candidates = h1(base + "candidates", ";Cumulative selection;Candidates", 7, -0.5, 6.5, true);
        q.candidatesVsPt = h2(base + "candidatesVsPt", ";Cumulative selection;Cluster pT (GeV/c)",
                              {7, -0.5, 6.5}, {250, 0, 500}, true);
        if (sample == 0) {
          q.goodSignal = h1(base + "goodSignal", ";Cumulative selection;Prompt-dominated candidates", 7, -0.5, 6.5, true);
          q.goodSignalVsPt = h2(base + "goodSignalVsPt", ";Cumulative selection;Cluster pT (GeV/c)",
                                {7, -0.5, 6.5}, {250, 0, 500}, true);
        }
        for (int bin = 1; bin <= NAnalysisStages; ++bin) {
          q.candidates->setBinLabel(bin, labels[bin - 1]);
          if (q.goodSignal) {
            q.goodSignal->setBinLabel(bin, labels[bin - 1]);
          }
          for (auto const& hist : {q.candidatesVsPt, q.goodSignalVsPt}) {
            if (!hist) {
              continue;
            }
            hist->weighted->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
            if (hist->unweighted) {
              hist->unweighted->GetXaxis()->SetBinLabel(bin, labels[bin - 1]);
            }
          }
        }
      }
    }
  }

  void init(InitContext const&)
  {
    bookAnalysisCuts();

    if (doprocessMCGen && doprocessMCGenJE) {
      throw std::runtime_error("Select only one GEN input: processMCGen (truth only) or processMCGenJE (JE truth and collision labels)");
    }
    if (doprocessMCGen && cfgGenReqRec) {
      throw std::runtime_error("cfgGenReqRec requires JE labels: disable processMCGen and enable processMCGenJE, or run REC alone");
    }
    if (cfgMinPromptFraction <= 0.f || cfgMinPromptFraction > 1.f) {
      throw std::runtime_error("cfgMinPromptFraction must be in (0,1]");
    }
    if (cfgLowM02 > cfgHighM02 || cfgMinTime > cfgMaxTime || cfgLowClusterE > cfgHighClusterE) {
      throw std::runtime_error("Invalid cluster selection bounds");
    }
    const bool validIsolation = cfgIsoRadius > 0.f && cfgIsoMaxSumPt >= 0.f && cfgIsoMinEOverP >= 0.f;
    if (!validIsolation) {
      throw std::runtime_error("Invalid isolation radius or threshold");
    }
    const bool validMatchingGeometry = cfgMatchMinEOverP > 0.f && cfgDCGapEta >= 0.f && cfgDCGapEta < cfgEMChighPSeta &&
                                       cfgDClowPSphi < cfgDCGapEndPhi && cfgDCGapEndPhi <= cfgDChighPSphi;
    if (!validMatchingGeometry) {
      throw std::runtime_error("Invalid matched E/P or DCal geometry");
    }

    for (int sample = 0; sample < NSamples; ++sample) {
      auto& q = qaState.tracks[sample];
      const std::string base = sample == 0 ? "REC/tracking/" : "DATA/tracking/";
      q.global = bookTrackKinematics(base + "global");
      q.extrapolated = bookTrackKinematics(base + "surfaceExtrapolated");
      q.inCalo = bookTrackKinematics(base + "surfaceInCalo");
      q.surfaceStatus = h1(base + "surfaceStatus", ";Status;Selected tracks", 3, -0.5, 2.5, true);
      q.surfaceStatus->setBinLabel(1, "missing/invalid extrapolation");
      q.surfaceStatus->setBinLabel(2, "valid outside calo");
      q.surfaceStatus->setBinLabel(3, "valid inside calo");
      q.etaShift = h2(base + "etaShift", ";Global pT;eta(surface)-eta(global)", {250, 0, 50}, {200, -1, 1}, true);
      q.phiShift = h2(base + "phiShift", ";Global pT;phi(surface)-phi(global)", {250, 0, 50}, {200, -o2::constants::math::PI, o2::constants::math::PI}, true);
    }

    const std::array<std::string, 5> trackTruthNames{
      "allGenerated", "acceptedEvents", "globalMatched", "surfaceExtrapolatedMatched", "surfaceInCaloMatched"};

    for (int i = 0; i < NTrackingTruthSets; ++i) {
      qaState.trackingTruth[i] = bookTrackKinematics("REC/tracking/performance/" + trackTruthNames[i]);
    }

    qaState.trackingLabelStatus = h1("REC/tracking/labelStatus", ";Status;Selected tracks in accepted events", 5, -0.5, 4.5, true);
    const std::array<const char*, NTrackingTruthSets> trackLabelNames{"missing label", "invalid label", "different MC collision", "not selected charged primary", "matched charged primary"};
    for (int i = 0; i < NTrackingTruthSets; ++i) {
      qaState.trackingLabelStatus->setBinLabel(i + 1, trackLabelNames[i]);
    }
    qaState.ancestryStatus = h1("GEN/ancestry", ";Classification;Photons in fiducial acceptance", 5, -0.5, 4.5, true);
    const std::array<const char*, NAncestryClasses> ancestryLabels{"other", "direct", "fragmentation", "mixed", "lepton ancestry only or unresolved"};
    for (int i = 0; i < NAncestryClasses; ++i) {
      qaState.ancestryStatus->setBinLabel(i + 1, ancestryLabels[i]);
    }
    qaState.genIsolation = bookIsolation("GEN/isolation/prompt");
    qaState.allPrimaryMap = bookMap("GEN/allPrimaryParticleEtaPhi");
    qaState.allPhotonMap = bookMap("GEN/allPrimaryPhotonEtaPhi");
    qaState.genParticleMap = bookMap("GEN/chargedPrimaryEtaPhi");
    triggerMaskBits = jetderiveddatautilities::initialiseTriggerMaskBits(cfgTriggerMasks);
    trackFilter = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(cfgTrackFilter));
    const std::array<std::string, 4> stages{"beforeM02", "afterM02", "neutral", "charged"};
    for (int sample = 0; sample < NSamples; ++sample) {
      for (int stage = 0; stage < NClusterSelections; ++stage) {
        qaState.clusterEnergy[sample][stage] = bookClusterEnergy(
          (sample == 0 ? "REC/clusterEnergy/" : "DATA/clusterEnergy/") + stages[stage], sample == 0);
      }
    }

    for (int sample = 0; sample < NSamples; ++sample) {
      const std::string p = sample == 0 ? "REC" : "DATA";
      qaState.recoIsolation[sample] = bookIsolation(p + "/isolation/neutral");
      qaState.allClusterMap[sample] = bookMap(p + "/allClusterEtaPhi");
      qaState.trackMap[sample] = bookMap(p + "/selectedTrackEtaPhi");
      qaState.matchedEOverP[sample] = h2(p + "/association/eOverP", ";Track p (GeV/c);Cluster E / track p", {250, 0, 250}, {200, 0, 10}, true);
      qaState.nLinked[sample] = h1(p + "/association/nLinked", ";Linked tracks;Clusters after M02", 32, -0.5, 31.5, true);
      qaState.nLinkedValid[sample] = h1(p + "/association/nLinkedValid", ";Linked tracks;Clusters with valid matching", 32, -0.5, 31.5, true);
      qaState.invalidMatchingPt[sample] = h1(p + "/association/invalidMatchingPt", ";Cluster pT;Clusters with unusable matching", 250, 0, 500, true);
      qaState.recoveredIsoPt[sample] = h1(p + "/association/recoveryIsolatedPt", ";Cluster pT;Isolated unmatched or recovered clusters", 250, 0, 500, true);
      qaState.strictPt[sample] = h1(p + "/association/strictNeutralPt", ";Cluster pT;Unmatched clusters", 250, 0, 500, true);
      qaState.recoveryPt[sample] = h1(p + "/association/recoveryPt", ";Cluster pT;Unmatched or recovered clusters", 250, 0, 500, true);
      qaState.nMatched[sample] = h1(p + "/association/nMatched", ";Tracks passing matching window;Clusters with valid matching", 32, -0.5, 31.5, true);
      qaState.eventCuts[sample] = h1(p + "/eventCutFlow", ";Selection;Weighted events", 7, -0.5, 6.5);
      const std::array<const char*, NEventCuts> labels{"valid weight", "vertex", "sel8", "EMC readout", "trigger", "event track", "missing EMC extrapolation"};
      for (int b = 0; b < NEventCuts; ++b) {
        qaState.eventCuts[sample]->setBinLabel(b + 1, labels[b]);
      }
      qaState.clusterCuts[sample] = h1(p + "/clusterCutFlow", ";Selection;Weighted clusters", 9, -0.5, 8.5);
      const std::array<const char*, NClusterCuts> clabels{"all", "definition", "time", "minimum energy", "cells", "NLM", "exotic", "energy+acceptance", "M02"};
      for (int b = 0; b < NClusterCuts; ++b) {
        qaState.clusterCuts[sample]->setBinLabel(b + 1, clabels[b]);
      }
      qaState.sumPt[sample] = h1(p + "/matchedSumPt", ";Sum matched track pT (GeV/c);Weighted clusters", 250, 0, 250);
      qaState.eOverSumPt[sample] = h1(p + "/energyOverSumPt", ";E / sum pT;Weighted matched clusters", 200, 0, 10);
      qaState.matchResidual[sample] = h2(p + "/trackMatchResidual", ";Delta phi at EMCal;Delta eta at EMCal", {200, -0.1, 0.1}, {200, -0.1, 0.1});
      for (int st = 0; st < NClusterSelections; ++st) {
        bookCluster(qaState.blindQA[sample][st], p + "/" + stages[st]);
      }
    }

    const std::array<std::string, 3> truthNames{"prompt", "direct", "fragmentation"};
    for (int g = 0; g < NTruthGroups; ++g) {
      for (int st = 0; st < NClusterSelections; ++st) {
        bookCluster(qaState.truthQA[g][st], "REC/truth_" + truthNames[g] + "/" + stages[st]);
      }
    }

    qaState.truthStatus = h1("REC/truthStatus", ";Status;Weighted clusters before M02", 4, -0.5, 3.5);
    const std::array<const char*, NTruthLabelStatuses> tlabels{"valid labels", "no labels", "size mismatch", "invalid label/fraction"};
    for (int b = 0; b < NTruthLabelStatuses; ++b) {
      qaState.truthStatus->setBinLabel(b + 1, tlabels[b]);
    }
    qaState.maxPromptFraction = h1("REC/maxPromptFraction", ";Largest single prompt fraction;Weighted clusters", 101, 0, 1.01);
    qaState.totalPromptFraction = h1("REC/totalPromptFraction", ";Sum prompt energy fractions;Weighted clusters", 101, 0, 1.01);
    qaState.promptFractionEnergy = h2("REC/maxPromptFractionVsEnergy", ";Cluster E (GeV);Largest prompt fraction", {250, 0, 500}, {101, 0, 1.01});
    const std::array<std::string, 2> performanceGroups{"prompt", "direct"};

    for (int st = 0; st < nPerformanceStages; ++st) {
      qaState.validLabels[st] = h1("REC/performance/purity/" + performanceStages[st] + "/validLabels",
                                   ";Cluster pT;Candidates with valid contributor labels", 250, 0, 500, true);
    }

    for (int st = 0; st < nPerformanceStages; ++st) {
      qaState.purityDenominator[st] = h1("REC/performance/purity/" + performanceStages[st] + "/allCandidates",
                                         ";Reconstructed cluster pT (GeV/c);Weighted clusters", 250, 0, 500, st >= Isolated);
    }

    for (int g = 0; g < NPhotonGroups; ++g) {
      const auto base = "REC/performance/" + performanceGroups[g];
      qaState.coneAll[g] = h1(base + "/efficiency/coneEligibleGenerated", ";True photon pT;Cone-eligible photons", 250, 0, 500, true);
      qaState.coneAccepted[g] = h1(base + "/efficiency/coneEligibleAcceptedEvents", ";True photon pT;Cone-eligible photons in accepted events", 250, 0, 500, true);
      qaState.selectedEventTruth[g] = h1(base + "/efficiency/selectedEventPopulation", ";True photon pT;Fiducial photons in configured GEN event population", 250, 0, 500, true);
      qaState.selectedConeTruth[g] = h1(base + "/efficiency/selectedConePopulation", ";True photon pT;Cone-eligible photons in configured GEN event population", 250, 0, 500, true);
      qaState.efficiencyAll[g] = h1(base + "/efficiency/allGenerated", ";True photon pT (GeV/c);Weighted photons", 250, 0, 500);
      qaState.efficiencyAccepted[g] = h1(base + "/efficiency/acceptedEvents", ";True photon pT (GeV/c);Weighted photons", 250, 0, 500);
      for (int st = 0; st < nPerformanceStages; ++st) {
        qaState.efficiencyNumerator[g][st] = h1(base + "/efficiency/" + performanceStages[st], ";True photon pT (GeV/c);Weighted matched photons", 250, 0, 500, st >= Isolated);
        qaState.purityNumerator[g][st] = h1(base + "/purity/" + performanceStages[st], ";Reconstructed cluster pT (GeV/c);Weighted tagged clusters", 250, 0, 500, st >= Isolated);
      }
    }
    // Use one event ensemble on both sides of the acceptance ratio.
    const std::string acceptance = "REC/performance/direct/acceptance/";
    qaState.directBeforeCalo = h1(acceptance + "beforeCalo", ";True photon pT;Direct photons before geometry", 250, 0, 500, true);
    qaState.directInCalo = h1(acceptance + "inCalo", ";True photon pT;Direct photons within geometry", 250, 0, 500, true);
    for (int st = 0; st < nPerformanceStages; ++st) {
      const auto base = "REC/performance/direct/fiducialPurity/" + performanceStages[st];
      qaState.directFiducialTagged[st] = h1(base + "/tagged", ";Cluster pT;Same-collision fiducial direct tags", 250, 0, 500, true);
      qaState.fiducialCandidates[st] = h1(base + "/candidates", ";Cluster pT;Fiducial candidates", 250, 0, 500, true);
    }
    // Experimental selection defines the denominator; truth only grades its signal content.
    const std::array<std::string, 3> signalSamples{"neutral", "neutralConeEligible", "isolated"};
    for (int i = 0; i < NSignalSamples; ++i) {
      auto& q = qaState.signalValidation[i];
      const std::string base = "REC/signalValidation/" + signalSamples[i] + "/";
      const std::array<std::pair<const char*, std::shared_ptr<QA1>*>, 8> entries{{{"allCandidates", &q.allCandidates}, {"validTruth", &q.validTruth}, {"anyPrompt", &q.anyPrompt}, {"goodSignal", &q.goodSignal}, {"badSignal", &q.badSignal}, {"noPrompt", &q.noPrompt}, {"promptBelowThreshold", &q.promptBelowThreshold}, {"unknownTruth", &q.unknownTruth}}};
      for (auto const& entry : entries) {
        *entry.second = h1(base + entry.first, ";Cluster pT (GeV/c);Candidates", 250, 0, 500, true);
      }
      q.maxPromptFraction = h2(base + "maxPromptFractionVsPt",
                               ";Cluster pT (GeV/c);Largest single prompt deposited-energy fraction",
                               {250, 0, 500}, {101, 0, 1.01}, true);
      q.totalPromptFraction = h2(base + "totalPromptFractionVsPt",
                                 ";Cluster pT (GeV/c);Sum of prompt deposited-energy fractions",
                                 {250, 0, 500}, {101, 0, 1.01}, true);
    }
    // Match the reconstructed truth-tagged cluster angular binning.
    qaState.genPromptEta = h1("GEN/prompt/eta", ";eta;Prompt photons", 80, -0.8, 0.8);
    qaState.genPromptPhi = h1("GEN/prompt/phi", ";phi;Prompt photons", 90, 0, o2::constants::math::TwoPI);
    qaState.genEvents = h1("GEN/nEvents", ";MC events;Sum weights", 1, 0, 1);
    for (int g = 0; g < NPhotonGroups; ++g) {
      std::string p = g == 0 ? "GEN/direct" : "GEN/frag";
      qaState.genPhotonMap[g] = bookMap(p + "_etaPhiQA");
      qaState.genPhi[g] = h1(p + "_phiQA", ";phi;Weighted photons", 320, 0, o2::constants::math::TwoPI);
      qaState.genEta[g] = h1(p + "_etaQA", ";eta;Weighted photons", 50, -1, 1);
      qaState.genPt[g] = h1(p + "_ptQA", ";pT (GeV/c);Weighted photons", 250, 0, 500);
    }
  }
  //*********************************************
  // INIT COMPLETE
  // NOW TIME FOR ANALYSIS-LEVEL HELPER FUNCTIONS
  //*********************************************

  // Truth classification: follow photon copies, then only quark/gluon ancestors
  // of a parton emission. Stop each branch at its nearest hard endpoint.
  template <typename ParticleTable, typename Particle, typename Table>
  OriginalResult classifyOriginal(Particle const& photon, Table const& particles) const
  {
    OriginalResult result;
    const int64_t first = particles.offset();
    const int64_t count = particles.size();
    if (!photon.isPhysicalPrimary() || photon.getGenStatusCode() < 0 || photon.getGenStatusCode() >= FinalGeneratorStatusEnd) {
      return result;
    }
    std::vector<std::pair<int64_t, bool>> pending;
    std::set<std::pair<int64_t, bool>> visited;
    auto addMothers = [&](auto const& row, bool partonBranch) {
      for (auto const& mother : row.template mothers_as<ParticleTable>()) {
        pending.emplace_back(mother.globalIndex(), partonBranch);
      }
    };
    addMothers(photon, false);
    while (!pending.empty()) {
      auto [id, partonBranch] = pending.back();
      pending.pop_back();
      const auto local = id - first;
      if (local < 0 || local >= count ||
          !visited.emplace(id, partonBranch).second) {
        continue;
      }
      auto mother = particles.rawIteratorAt(local);
      const int code = std::abs(mother.pdgCode());
      const bool parton = (code >= kDown && code <= kTop) || code == kGluon;
      if ((!parton && code != kGamma) || (partonBranch && code == kGamma)) {
        result.leptonBranch |= code == kElectron || code == kMuonMinus || code == kTauMinus;
        continue;
      }
      if (mother.getGenStatusCode() == HardProcessStatus) {
        auto& endpoints = parton ? result.fragmentationHardIndices : result.directHardIndices;
        endpoints.push_back(id);
        continue;
      }
      addMothers(mother, parton);
    }
    const bool direct = !result.directHardIndices.empty();
    const bool frag = !result.fragmentationHardIndices.empty();
    result.origin = direct ? (frag ? Origin::Mixed : Origin::Direct)
                           : (frag ? Origin::Fragmentation : Origin::Unclassified);
    return result;
  }

  // Track, calorimeter and cone selections
  template <typename T>
  static bool finiteTrack(T const& t)
  {
    return std::isfinite(t.pt()) && std::isfinite(t.eta()) && std::isfinite(t.phi());
  }

  template <typename TrackType>
  bool trackSelection(const TrackType& track)
  {
    if (!finiteTrack(track) || !std::isfinite(track.dcaXY()) || !std::isfinite(track.dcaZ()) ||
        !std::isfinite(track.tpcCrossedRowsOverFindableCls()) || !std::isfinite(track.tpcChi2NCl()) ||
        !std::isfinite(track.itsChi2NCl())) {
      return false;
    }
    // Basic track cuts
    if (track.pt() < cfgtrkMinPt) {
      return false;
    }

    if (std::abs(track.eta()) > cfgtrkMaxEta) {
      return false;
    }

    if (std::abs(track.dcaXY()) > cfgMaxDCArToPVcut) {
      return false;
    }

    if (std::abs(track.dcaZ()) > cfgMaxDCAzToPVcut) {
      return false;
    }

    if (cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }

    if (track.tpcNClsFindable() < cfgnFindableTPCClusters) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < cfgnTPCCrossedRows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() > cfgnRowsOverFindable) {
      return false;
    }

    if (track.tpcChi2NCl() > cfgTpcChi2) {
      return false;
    }

    if (track.itsChi2NCl() > cfgItsChi2) {
      return false;
    }

    if (cfgConnectedToPV && !track.isPVContributor()) {
      return false;
    }

    return true;
  }; // end of track selection

  template <typename T>
  bool selectedTrack(T const& track)
  {
    if (!finiteTrack(track) || track.pt() < cfgtrkMinPt || std::abs(track.eta()) > cfgtrkMaxEta) {
      return false;
    }
    if (cfgJETracks) {
      return jetderiveddatautilities::selectTrack(track, trackFilter);
    }
    auto original = track.template track_as<TrackCandidates>();
    return trackSelection(original) && original.isGlobalTrack() &&
           (!cfgGlobalWoDCATrack || original.isGlobalTrackWoDCA());
  }

  bool coneEligible(double eta)
  {
    return std::abs(eta) <= cfgtrkMaxEta &&
           (!cfgIsoFullCone || std::abs(eta) + cfgIsoRadius <= cfgtrkMaxEta);
  }

  bool inCalo(double eta, double phi)
  {
    if (!(std::abs(eta) < cfgEMChighPSeta)) {
      return false;
    }
    phi = phi02pi(phi);
    if (phi > cfgEMClowPSphi && phi < cfgEMChighPSphi) {
      return true;
    }
    return phi > cfgDClowPSphi && phi < cfgDChighPSphi &&
           (phi >= cfgDCGapEndPhi || std::abs(eta) > cfgDCGapEta);
  }

  // Charged isolation optionally excludes tracks matched to this candidate.
  // Cone E/sum|p| uses the same constituents as sumPt: after a tight
  // sumPt cut its extra rejection vanishes above threshold*cosh(etaMax)*maxSumPt.
  std::pair<bool, bool> measureIsolation(IsolationQA* qa, double pt, double eta, double phi,
                                         double energy, std::vector<ConeParticle> const& particles, double w,
                                         std::set<int64_t> const* excludedTracks = nullptr)
  {
    if (!coneEligible(eta)) {
      return {false, false};
    }
    double sumPt = 0., sumP = 0.;
    for (auto const& particle : particles) {
      if (excludedTracks && excludedTracks->contains(particle.trackId)) {
        continue;
      }
      const double dEta = particle.eta - eta;
      const double dPhi = TVector2::Phi_mpi_pi(particle.phi - phi);
      if (dEta * dEta + dPhi * dPhi >= cfgIsoRadius * cfgIsoRadius) {
        continue;
      }
      sumPt += particle.pt;
      sumP += particle.pt * std::cosh(particle.eta);
    }
    const bool isolated = sumPt < cfgIsoMaxSumPt;
    const bool passesRatio = sumP == 0. || energy >= cfgIsoMinEOverP * sumP;
    if (qa) {
      qa->eligible->fill(pt, w);
      qa->sumPt->fill(pt, sumPt, w);
      if (sumP > 0.) {
        qa->eOverP->fill(pt, energy / sumP, w);
      }
      if (isolated) {
        qa->isolated->fill(pt, w);
      }
      if (passesRatio) {
        qa->eOverPPass->fill(pt, w);
      }
      if (isolated && passesRatio) {
        qa->combined->fill(pt, w);
      }
    }
    return {isolated, passesRatio};
  }

  // Cluster truth tagging and performance counting
  template <typename Particle>
  bool isDecayPhoton(Particle const& photon, aod::JMcParticles const& particles)
  {
    const int64_t first = particles.offset();
    const int64_t count = particles.size();
    const int64_t end = first + count;
    std::vector<int64_t> pending;
    std::set<int64_t> visited;
    for (auto const& mother : photon.template mothers_as<aod::JMcParticles>()) {
      pending.push_back(mother.globalIndex());
    }
    while (!pending.empty()) {
      const int64_t id = pending.back();
      pending.pop_back();
      if (id < first || id >= end || !visited.insert(id).second) {
        continue;
      }
      const auto mother = particles.rawIteratorAt(id - first);
      const int code = std::abs(mother.pdgCode());
      if (code > HadronCodeLowerBound && code < NucleusCodeLowerBound) {
        return true;
      }
      if (code != kGamma) {
        continue;
      }
      for (auto const& ancestor : mother.template mothers_as<aod::JMcParticles>()) {
        pending.push_back(ancestor.globalIndex());
      }
    }
    return false;
  }

  template <typename Particle>
  ClusterOrigin contributorOrigin(Particle const& particle, aod::JMcParticles const& particles)
  {
    if (particle.pdgCode() == kGamma) {
      const auto result = classifyOriginal<aod::JMcParticles>(particle, particles);
      if (result.origin == Origin::Direct) {
        return ClusterOrigin::Direct;
      }
      if (result.origin == Origin::Fragmentation) {
        return ClusterOrigin::Fragmentation;
      }
      if (result.origin == Origin::Mixed) {
        return ClusterOrigin::MixedPrompt;
      }
      return isDecayPhoton(particle, particles) ? ClusterOrigin::Decay : ClusterOrigin::Other;
    }
    // Classify the labelled contributor itself. Do not double-count its ancestors.
    auto* species = pdg->GetParticle(particle.pdgCode());
    return species && species->Charge() != 0. ? ClusterOrigin::Charged : ClusterOrigin::Other;
  }

  struct TruthTag {
    ClusterOrigin stackOrigin = ClusterOrigin::UnknownTruth;
    bool fiducialDirect = false, coneFiducialDirect = false;
    bool labelsValid = false;
    bool anyPrompt = false;
    double maxPromptFraction = 0., totalPromptFraction = 0.;
    bool prompt = false, direct = false, fragmentation = false;
    std::array<std::vector<int64_t>, 2> efficiencyIds;
  };

  template <typename C>
  TruthTag tagCluster(C const& c, aod::JMcParticles const& particles, double w)
  {
    TruthTag tag;
    const auto ids = c.mcParticlesIds();
    const auto fractions = c.amplitudeA(); // EMCAL producer stores getEnergyFractions(), NOT energies.
    if (ids.size() != fractions.size()) {
      qaState.truthStatus->fill(2, w);
      return tag;
    }
    if (ids.empty()) {
      qaState.truthStatus->fill(1, w);
      return tag;
    }
    // Validate the complete association before assigning any tag. Never index invalid labels.
    const int64_t first = particles.offset();
    const int64_t end = first + particles.size();
    double fractionSum = 0.;
    for (std::size_t i = 0; i < ids.size(); ++i) {
      if (ids[i] < first || ids[i] >= end || !std::isfinite(fractions[i]) || fractions[i] < 0.f || fractions[i] > MaxContributorFraction) {
        qaState.truthStatus->fill(3, w);
        return tag;
      }
      fractionSum += fractions[i];
    }
    if (fractionSum <= 0. || fractionSum > MaxContributorFractionSum) {
      qaState.truthStatus->fill(3, w);
      return tag;
    }
    tag.labelsValid = true;
    tag.stackOrigin = ClusterOrigin::NoDominant;
    double leadingFraction = -1.;
    int64_t leadingId = std::numeric_limits<int64_t>::max();
    qaState.truthStatus->fill(0.0, w);
    double maximum = 0., total = 0.;
    std::size_t i = 0;
    for (auto const& particle : c.template mcParticles_as<aod::JMcParticles>()) {
      const double fraction = fractions[i++];
      // The largest qualifying contributor gets exactly one stack category.
      // Below a 50% threshold several contributors can qualify; ties use the lowest ID.
      if (fraction > cfgMinPromptFraction &&
          (fraction > leadingFraction || (fraction == leadingFraction && particle.globalIndex() < leadingId))) {
        leadingFraction = fraction;
        leadingId = particle.globalIndex();
        tag.stackOrigin = contributorOrigin(particle, particles);
      }
      if (particle.pdgCode() != kGamma) {
        continue;
      }
      const auto result = classifyOriginal<aod::JMcParticles>(particle, particles);
      if (!isPrompt(result.origin)) {
        continue;
      }
      tag.anyPrompt |= fraction > 0.;
      maximum = std::max(maximum, fraction);
      total += fraction;
      if (fraction <= cfgMinPromptFraction) {
        continue;
      }
      tag.efficiencyIds[0].push_back(particle.globalIndex());
      if (!result.directHardIndices.empty()) {
        tag.efficiencyIds[1].push_back(particle.globalIndex());
      }
      tag.prompt = true;
      tag.direct |= !result.directHardIndices.empty();
      tag.fragmentation |= !result.fragmentationHardIndices.empty();
    }
    tag.maxPromptFraction = maximum;
    tag.totalPromptFraction = total;
    qaState.maxPromptFraction->fill(maximum, w);
    qaState.totalPromptFraction->fill(total, w);
    qaState.promptFractionEnergy->fill(c.energy(), maximum, w);
    return tag;
  }

  void fillSignalValidation(int sample, double pt, double w, TruthTag const& tag)
  {
    auto& q = qaState.signalValidation[sample];
    q.allCandidates->fill(pt, w);
    // Untaggable candidates are not proven signal; retain them explicitly as unknown.
    if (!tag.labelsValid) {
      q.unknownTruth->fill(pt, w);
      q.badSignal->fill(pt, w);
      return;
    }
    q.validTruth->fill(pt, w);
    q.maxPromptFraction->fill(pt, tag.maxPromptFraction, w);
    q.totalPromptFraction->fill(pt, tag.totalPromptFraction, w);
    if (tag.anyPrompt) {
      q.anyPrompt->fill(pt, w);
    }
    if (tag.prompt) {
      q.goodSignal->fill(pt, w);
    } else {
      q.badSignal->fill(pt, w);
      if (tag.anyPrompt) {
        q.promptBelowThreshold->fill(pt, w);
      } else {
        q.noPrompt->fill(pt, w);
      }
    }
  }

  void fillClusterEnergy(ClusterEnergyQA& q, double energy, double w, TruthTag const& tag)
  {
    q.total->fill(energy, w);
    if (!q.prompt) {
      return;
    } // Data has no truth decomposition.
    const int category = static_cast<int>(tag.stackOrigin);
    q.exclusive[category]->fill(energy, w);
    if (tag.stackOrigin == ClusterOrigin::Direct || tag.stackOrigin == ClusterOrigin::Fragmentation ||
        tag.stackOrigin == ClusterOrigin::MixedPrompt) {
      q.prompt->fill(energy, w);
    }
  }

  void fillAnalysisCut(int sample, int branch, int stage, double pt, double energy, double w, TruthTag const& tag)
  {
    auto& q = qaState.analysisCuts[sample][branch];
    q.candidates->fill(stage, w);
    q.candidatesVsPt->fill(stage, pt, w);
    fillClusterEnergy(q.energy[stage], energy, w, tag);
    if (tag.prompt && q.goodSignal) {
      q.goodSignal->fill(stage, w);
      q.goodSignalVsPt->fill(stage, pt, w);
    }
  }

  void fillAnalysisCommonCut(int sample, int stage, double pt, double energy, double w, TruthTag const& tag)
  {
    for (int branch = 0; branch < NNeutralBranches; ++branch) {
      fillAnalysisCut(sample, branch, stage, pt, energy, w, tag);
    }
  }

  void fillAnalysisFinalCuts(int sample, int branch, double pt, double energy, double w, TruthTag const& tag,
                             bool neutral, bool eligible, bool isolated, bool passesRatio)
  {
    if (!neutral) {
      return;
    }
    fillAnalysisCut(sample, branch, 3, pt, energy, w, tag);
    if (!eligible) {
      return;
    }
    fillAnalysisCut(sample, branch, 4, pt, energy, w, tag);
    if (!isolated) {
      return;
    }
    fillAnalysisCut(sample, branch, 5, pt, energy, w, tag);
    if (passesRatio) {
      fillAnalysisCut(sample, branch, 6, pt, energy, w, tag);
    }
  }

  template <typename C>
  void fillPerformance(int stage, C const& c, double w, TruthTag const& tag)
  {
    const double pt = c.energy() / std::cosh(c.eta());
    qaState.purityDenominator[stage]->fill(pt, w);
    if (stage == Neutral) {
      fillSignalValidation(0, pt, w, tag);
      if (coneEligible(c.eta())) {
        fillSignalValidation(1, pt, w, tag);
      }
    } else if (stage == Isolated) {
      fillSignalValidation(2, pt, w, tag);
    }
    // Only matches inside the reconstructed fiducial region enter efficiency.
    // Inclusive purity and blind QA can still cover the larger REC sample.
    if (inCalo(c.eta(), c.phi()) && (!isolatedStage(stage) || coneEligible(c.eta()))) {
      for (int g = 0; g < NPhotonGroups; ++g) {
        qaState.efficiencyMatched[g][stage].insert(tag.efficiencyIds[g].begin(), tag.efficiencyIds[g].end());
      }
      qaState.fiducialCandidates[stage]->fill(pt, w);
      if (isolatedStage(stage) ? tag.coneFiducialDirect : tag.fiducialDirect) {
        qaState.directFiducialTagged[stage]->fill(pt, w);
      }
    }
    if (tag.labelsValid) {
      qaState.validLabels[stage]->fill(pt, w);
    }
    if (tag.prompt) {
      qaState.purityNumerator[0][stage]->fill(pt, w);
    }
    if (tag.direct) {
      qaState.purityNumerator[1][stage]->fill(pt, w);
    }
  }

  template <typename C>
  void fillStage(int sample, int stage, C const& c, double w, TruthTag const& tag)
  {
    fillCluster(qaState.blindQA[sample][stage], c, w);
    fillClusterEnergy(qaState.clusterEnergy[sample][stage], c.energy(), w, tag);
    if (sample != 0) {
      return;
    }
    if (stage <= Neutral) {
      fillPerformance(stage, c, w, tag);
    }
    if (tag.prompt) {
      fillCluster(qaState.truthQA[0][stage], c, w);
    }
    if (tag.direct) {
      fillCluster(qaState.truthQA[1][stage], c, w);
    }
    if (tag.fragmentation) {
      fillCluster(qaState.truthQA[2][stage], c, w);
    }
  }

  // The first six cut-flow bins include clusters rejected by the framework
  // filter. Read only the bare cluster table for these diagnostic counts.
  void fillBaselineCutFlow(aod::JClusters const& rawClusters, int sample, double w)
  {
    for (auto const& c : rawClusters) {
      qaState.clusterCuts[sample]->fill(0.0, w);
      if (c.definition() != cfgClusterDefinition) {
        continue;
      }
      if (std::isfinite(c.eta()) && std::isfinite(c.phi())) {
        qaState.allClusterMap[sample]->fill(c.eta(), phi02pi(c.phi()), w);
      }
      qaState.clusterCuts[sample]->fill(1, w);
      if (!std::isfinite(c.time()) || c.time() < cfgMinTime || c.time() > cfgMaxTime) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(2, w);
      if (!std::isfinite(c.energy()) || c.energy() <= cfgMinClusterEnergy) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(3, w);
      if (c.nCells() < cfgMinNCells) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(4, w);
      if (c.nlm() > cfgMaxNLM) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(5, w);
    }
  }

  // Use the same extrapolation source as the cluster-track matching branch.
  template <typename Track>
  bool trackSurface(Track const& track, aod::JEMCTracks const& emctracks, double& eta, double& phi)
  {
    if (cfgJETracks) {
      auto rows = emctracks.sliceBy(emcTrackPerTrack, track.globalIndex());
      if (rows.size() == 0) {
        return false;
      }
      eta = rows.iteratorAt(0).etaEmcal();
      phi = rows.iteratorAt(0).phiEmcal();
    } else {
      auto original = track.template track_as<TrackCandidates>();
      eta = original.trackEtaEmcal();
      phi = original.trackPhiEmcal();
    }
    // Reject extrapolation sentinels before wrapping phi.
    if (!std::isfinite(eta) || !std::isfinite(phi) || std::abs(eta) > MaxExtrapolationAbsEta ||
        phi < -o2::constants::math::PI || phi > o2::constants::math::TwoPI) {
      return false;
    }
    phi = phi02pi(phi);
    return true;
  }

  template <typename Particle>
  bool trackingTruthSelected(Particle const& p)
  {
    if (!p.isPhysicalPrimary() || !finiteTrack(p) || p.pt() < cfgtrkMinPt ||
        std::abs(p.eta()) > cfgtrkMaxEta) {
      return false;
    }
    auto* species = pdg->GetParticle(p.pdgCode());
    return species && species->Charge() != 0.;
  }

  template <typename Tracks>
  void collectTrackingMatches(Tracks const& tracks, aod::JMcTrackLbs const& labels,
                              aod::JMcParticles const& particles, int64_t collisionId,
                              aod::JEMCTracks const& emctracks, double w,
                              std::array<std::set<int64_t>, 3>& matched)
  {
    const int64_t first = particles.offset();
    const int64_t count = particles.size();
    const int64_t end = first + count;
    for (auto const& track : tracks) {
      if (!selectedTrack(track)) {
        continue;
      }
      auto label = labels.rawIteratorAt(track.globalIndex() - labels.offset());
      if (!label.has_mcParticle()) {
        qaState.trackingLabelStatus->fill(0., w);
        continue;
      }
      const int64_t id = label.mcParticleId();
      if (id < first || id >= end) {
        qaState.trackingLabelStatus->fill(1., w);
        continue;
      }
      auto truth = particles.rawIteratorAt(id - first);
      if (!truth.has_mcCollision() || truth.mcCollisionId() != collisionId) {
        qaState.trackingLabelStatus->fill(2., w);
        continue;
      }
      if (!trackingTruthSelected(truth)) {
        qaState.trackingLabelStatus->fill(3., w);
        continue;
      }
      qaState.trackingLabelStatus->fill(4., w);
      matched[0].insert(id);
      double eta = 0., phi = 0.;
      if (!trackSurface(track, emctracks, eta, phi)) {
        continue;
      }
      matched[1].insert(id);
      if (inCalo(eta, phi)) {
        matched[2].insert(id);
      }
    }
  }

  struct ClusterMatch {
    bool valid, matched;
    double sumPt, sumP;
  };

  template <typename Cluster>
  ClusterMatch matchCluster(Cluster const& c, aod::JEMCTracks const& emctracks, int sample,
                            double w, double clusterPt, bool& missingExtrapolation,
                            std::set<int64_t>& matchedTrackIds)
  {
    double matchedPt = 0., matchedP = 0.;
    bool matched = false, validMatching = true;
    auto linkedTracks = c.template matchedTracks_as<RecoTracks>();
    qaState.nLinked[sample]->fill(linkedTracks.size(), w);
    int matchedCount = 0;
    for (auto const& track : linkedTracks) {
      if (!finiteTrack(track)) {
        validMatching = false;
        continue;
      }
      if (!selectedTrack(track)) {
        continue;
      }
      double etaT = 0., phiT = 0.;
      if (cfgJETracks) {
        auto extrapolated = emctracks.sliceBy(emcTrackPerTrack, track.globalIndex());
        if (extrapolated.size() == 0) {
          missingExtrapolation = true;
          validMatching = false;
          continue;
        }
        etaT = extrapolated.iteratorAt(0).etaEmcal();
        phiT = extrapolated.iteratorAt(0).phiEmcal();
      } else {
        auto original = track.template track_as<TrackCandidates>();
        etaT = original.trackEtaEmcal();
        phiT = original.trackPhiEmcal();
      }
      if (!std::isfinite(etaT) || !std::isfinite(phiT) || std::abs(etaT) > MaxExtrapolationAbsEta || phiT < MinMatchingExtrapolationPhi) {
        missingExtrapolation = true;
        validMatching = false;
        continue;
      }
      const double dEta = etaT - c.eta(), dPhi = TVector2::Phi_mpi_pi(phiT - c.phi());
      qaState.matchResidual[sample]->fill(dPhi, dEta, w);
      const double etaWindow = cfgPtClusterCut ? 0.010 + std::pow(track.pt() + 4.07, -2.5) : 0.05;
      const double phiWindow = cfgPtClusterCut ? 0.015 + std::pow(track.pt() + 3.65, -2.0) : 0.05;
      if (std::abs(dEta) < etaWindow && std::abs(dPhi) < phiWindow) {
        ++matchedCount;
        matchedTrackIds.insert(track.globalIndex());
        const double momentum = track.pt() * std::cosh(track.eta());
        if (momentum > 0.) {
          qaState.matchedEOverP[sample]->fill(momentum, c.energy() / momentum, w);
        }
        matchedP += momentum;
        matched = true;
        matchedPt += track.pt();
      }
    }

    // Incomplete matching cannot establish that a cluster is neutral.
    if (!validMatching) {
      qaState.invalidMatchingPt[sample]->fill(clusterPt, w);
      return {.valid = false, .matched = matched, .sumPt = matchedPt, .sumP = matchedP};
    }
    qaState.nLinkedValid[sample]->fill(linkedTracks.size(), w);
    qaState.nMatched[sample]->fill(matchedCount, w);
    return {.valid = true, .matched = matched, .sumPt = matchedPt, .sumP = matchedP};
  }

  // Shared reconstructed-event processing for MC and data.
  template <typename Collision, typename Clusters, typename Tagger>
  void processReco(Collision const& collision, Clusters const& clusters, aod::JClusters const& rawClusters,
                   RecoTracks const& tracks,
                   aod::JEMCTracks const& emctracks, int sample, double w, Tagger tagger)
  {
    qaState.recoEventAccepted = false;
    if (!std::isfinite(w)) {
      throw std::runtime_error("Non-finite MC event weight");
    }
    qaState.eventCuts[sample]->fill(0.0, w);
    if (!std::isfinite(collision.posZ()) || std::abs(collision.posZ()) >= cfgVtxCut) {
      return;
    }
    qaState.eventCuts[sample]->fill(1, w);
    if (!collision.sel8()) {
      return;
    }
    qaState.eventCuts[sample]->fill(2, w);
    if (cfgEmcTrigger && !collision.isEmcalReadout()) {
      return;
    }
    qaState.eventCuts[sample]->fill(3, w);
    if (!jetderiveddatautilities::selectTrigger(collision, triggerMaskBits)) {
      return;
    }
    qaState.eventCuts[sample]->fill(4, w);
    if (cfgRequireEventTrack) {
      bool found = false;
      for (auto const& track : tracks) {
        if (selectedTrack(track)) {
          found = true;
          break;
        }
      }
      if (!found) {
        return;
      }
    }
    qaState.eventCuts[sample]->fill(5, w);
    qaState.recoEventAccepted = true;
    bool missingExtrapolation = false;
    std::vector<ConeParticle> constituents;
    for (auto const& track : tracks) {
      if (!selectedTrack(track)) {
        continue;
      }
      constituents.push_back({track.pt(), track.eta(), track.phi(), track.globalIndex()});
      qaState.trackMap[sample]->fill(track.eta(), phi02pi(track.phi()), w);
      auto& tq = qaState.tracks[sample];
      fillTrackKinematics(tq.global, track.pt(), track.eta(), track.phi(), w);
      double surfaceEta = 0., surfacePhi = 0.;
      if (!trackSurface(track, emctracks, surfaceEta, surfacePhi)) {
        tq.surfaceStatus->fill(0., w);
      } else {
        // pT remains the global-track pT; only eta/phi are evaluated at the surface.
        fillTrackKinematics(tq.extrapolated, track.pt(), surfaceEta, surfacePhi, w);
        tq.etaShift->fill(track.pt(), surfaceEta - track.eta(), w);
        tq.phiShift->fill(track.pt(), TVector2::Phi_mpi_pi(surfacePhi - track.phi()), w);
        if (inCalo(surfaceEta, surfacePhi)) {
          fillTrackKinematics(tq.inCalo, track.pt(), surfaceEta, surfacePhi, w);
          tq.surfaceStatus->fill(2., w);
        } else {
          tq.surfaceStatus->fill(1., w);
        }
      }
    }
    fillBaselineCutFlow(rawClusters, sample, w);
    for (auto const& c : clusters) {
      qaState.clusterCuts[sample]->fill(6, w);
      if (!std::isfinite(c.eta()) || !std::isfinite(c.phi()) || !std::isfinite(c.m02()) || !std::isfinite(c.m20())) {
        continue;
      }
      if (c.energy() < cfgLowClusterE || c.energy() > cfgHighClusterE || (cfgReqRecPsRec && !inCalo(c.eta(), c.phi()))) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(7, w);
      auto tag = tagger(c);
      fillStage(sample, 0, c, w, tag);
      const double analysisPt = c.energy() / std::cosh(c.eta());
      fillAnalysisCommonCut(sample, 0, analysisPt, c.energy(), w, tag);
      if (c.m02() < cfgLowM02 || c.m02() > cfgHighM02) {
        continue;
      }
      qaState.clusterCuts[sample]->fill(8, w);
      fillStage(sample, 1, c, w, tag);
      fillAnalysisCommonCut(sample, 1, analysisPt, c.energy(), w, tag);
      const double clusterPt = c.energy() / std::cosh(c.eta());
      std::set<int64_t> matchedTrackIds;
      const auto [validMatching, matched, matchedPt, matchedP] =
        matchCluster(c, emctracks, sample, w, clusterPt, missingExtrapolation, matchedTrackIds);
      auto const* excludedTracks = cfgIsoExcludeMatchedTracks ? &matchedTrackIds : nullptr;
      if (!validMatching) {
        continue;
      }
      fillAnalysisCommonCut(sample, 2, clusterPt, c.energy(), w, tag);
      qaState.sumPt[sample]->fill(matchedPt, w);
      bool charged = matched;
      if (matched && matchedPt > 0.) {
        const double ratio = c.energy() / matchedPt;
        qaState.eOverSumPt[sample]->fill(ratio, w);
        if (cfgRecoverHighEOverSumPt && ratio >= cfgEOverSumPtMin) {
          charged = false;
        }
      }
      // Compare strict matching and E/sum|p| recovery on the same after-M02 population.
      // The legacy neutral branch above uses its independently configurable E/sumPt rule.
      const bool recoveryPass = !matched || (matchedP > 0. && c.energy() >= cfgMatchMinEOverP * matchedP);
      if (!matched) {
        qaState.strictPt[sample]->fill(clusterPt, w);
        if (sample == 0) {
          fillPerformance(StrictNeutral, c, w, tag);
        }
      }
      if (recoveryPass) {
        qaState.recoveryPt[sample]->fill(clusterPt, w);
        if (sample == 0) {
          fillPerformance(MatchRecovery, c, w, tag);
        }
      }
      // Cumulative diagnostics compare the branches without changing existing selections.
      const auto [analysisIsolated, analysisRatio] = measureIsolation(nullptr,
                                                                      clusterPt, c.eta(), c.phi(), c.energy(), constituents, w, excludedTracks);
      const bool eligible = coneEligible(c.eta());
      fillAnalysisFinalCuts(sample, 0, clusterPt, c.energy(), w, tag, !charged, eligible, analysisIsolated, analysisRatio);
      fillAnalysisFinalCuts(sample, 1, clusterPt, c.energy(), w, tag, !matched, eligible, analysisIsolated, analysisRatio);
      fillAnalysisFinalCuts(sample, 2, clusterPt, c.energy(), w, tag, recoveryPass, eligible, analysisIsolated, analysisRatio);
      fillStage(sample, charged ? 3 : 2, c, w, tag);
      if (!charged || recoveryPass) {
        const auto [isolated, passesRatio] = measureIsolation(charged ? nullptr : &qaState.recoIsolation[sample],
                                                              c.energy() / std::cosh(c.eta()), c.eta(), c.phi(), c.energy(), constituents, w, excludedTracks);
        if (isolated && recoveryPass) {
          qaState.recoveredIsoPt[sample]->fill(clusterPt, w);
          if (sample == 0) {
            fillPerformance(MatchRecoveryIsolated, c, w, tag);
          }
        }
        if (sample == 0 && isolated && !charged) {
          fillPerformance(Isolated, c, w, tag);
          if (passesRatio) {
            fillPerformance(IsolatedEOverP, c, w, tag);
          }
        }
      }
    }
    if (missingExtrapolation) {
      qaState.eventCuts[sample]->fill(6, w);
    }
  }
  // MC bookkeeping helpers: keep truth classification and efficiency counting together.
  template <typename Cluster>
  TruthTag tagClusterForCollision(Cluster const& c, aod::JMcParticles const& particles,
                                  int64_t mcCollisionId, double weight)
  {
    auto tag = tagCluster(c, particles, weight);
    // Purity retains all contributors. Efficiency matches only photons from
    // this reconstructed collision's associated MC collision.
    for (auto idsIt = tag.efficiencyIds.begin(); idsIt != tag.efficiencyIds.end(); ++idsIt) {
      auto& ids = *idsIt;
      std::vector<int64_t> sameCollision;
      for (auto const& id : ids) {
        const auto photon = particles.rawIteratorAt(id - particles.offset());
        if (photon.mcCollisionId() == mcCollisionId) {
          sameCollision.push_back(id);
        }
      }
      ids = std::move(sameCollision);
    }
    // efficiencyIds[1] now contains only direct contributors from this MC event.
    for (auto const& id : tag.efficiencyIds[1]) {
      const auto direct = particles.rawIteratorAt(id - particles.offset());
      if (!std::isfinite(direct.pt()) || !inCalo(direct.eta(), direct.phi())) {
        continue;
      }
      tag.fiducialDirect = true;
      tag.coneFiducialDirect |= coneEligible(direct.eta());
    }
    return tag;
  }

  void fillTrackingPerformance(aod::JMcParticles const& particles,
                               std::set<int64_t> const& acceptedMC,
                               std::array<std::set<int64_t>, 3> const& trackingMatched)
  {
    for (auto const& truth : particles) {
      if (!truth.has_mcCollision() || !trackingTruthSelected(truth)) {
        continue;
      }
      const auto mcCollision = truth.mcCollision_as<aod::JMcCollisions>();
      const double w = mcCollision.weight();
      if (!std::isfinite(w)) {
        throw std::runtime_error("Non-finite tracking truth weight");
      }
      fillTrackKinematics(qaState.trackingTruth[0], truth.pt(), truth.eta(), truth.phi(), w);
      if (acceptedMC.contains(mcCollision.globalIndex())) {
        fillTrackKinematics(qaState.trackingTruth[1], truth.pt(), truth.eta(), truth.phi(), w);
      }
      for (int i = 0; i < NTrackingMatches; ++i) {
        if (trackingMatched[i].contains(truth.globalIndex())) {
          fillTrackKinematics(qaState.trackingTruth[i + 2], truth.pt(), truth.eta(), truth.phi(), w);
        }
      }
    }
  }

  void fillPhotonPerformance(aod::JMcParticles const& particles,
                             std::set<int64_t> const& reconstructedMC,
                             std::set<int64_t> const& acceptedMC)
  {
    // Truth fiducial region is explicit and shared by numerator/denominator.
    // No reconstructed energy or shower-shape requirement is applied to truth.
    for (auto const& photon : particles) {
      if (photon.pdgCode() != kGamma || !std::isfinite(photon.pt()) || !std::isfinite(photon.eta()) || !std::isfinite(photon.phi())) {
        continue;
      }
      if (!photon.has_mcCollision()) {
        continue;
      }
      const auto mcCollision = photon.mcCollision_as<aod::JMcCollisions>();
      const double weight = mcCollision.weight();
      if (!std::isfinite(weight)) {
        throw std::runtime_error("Non-finite efficiency MC weight");
      }
      const auto result = classifyOriginal<aod::JMcParticles>(photon, particles);
      const bool selectedEvent = !cfgGenReqRec || reconstructedMC.contains(mcCollision.globalIndex());
      const bool direct = !result.directHardIndices.empty();
      if (direct && selectedEvent) {
        qaState.directBeforeCalo->fill(photon.pt(), weight);
      }
      if (!inCalo(photon.eta(), photon.phi())) {
        continue;
      }
      if (direct && selectedEvent) {
        qaState.directInCalo->fill(photon.pt(), weight);
      }
      const std::array<bool, 2> belongs{isPrompt(result.origin), !result.directHardIndices.empty()};
      for (int g = 0; g < NPhotonGroups; ++g) {
        if (!belongs[g]) {
          continue;
        }
        const bool eligible = coneEligible(photon.eta());
        if (selectedEvent) {
          qaState.selectedEventTruth[g]->fill(photon.pt(), weight);
          if (eligible) {
            qaState.selectedConeTruth[g]->fill(photon.pt(), weight);
          }
        }
        if (eligible) {
          qaState.coneAll[g]->fill(photon.pt(), weight);
          if (acceptedMC.contains(mcCollision.globalIndex())) {
            qaState.coneAccepted[g]->fill(photon.pt(), weight);
          }
        }
        qaState.efficiencyAll[g]->fill(photon.pt(), weight);
        if (acceptedMC.contains(mcCollision.globalIndex())) {
          qaState.efficiencyAccepted[g]->fill(photon.pt(), weight);
        }
        for (int st = 0; st < nPerformanceStages; ++st)
        // Pair isolated numerators with the truth cone-eligible denominators.
        {
          if ((eligible || !isolatedStage(st)) &&
              qaState.efficiencyMatched[g][st].contains(photon.globalIndex())) {
            qaState.efficiencyNumerator[g][st]->fill(photon.pt(), weight);
          }
        }
      }
    }
  }

  void fillMCReco(MCCollisions const& collisions, FilteredMCClusters const& clusters,
                  aod::JClusters const& rawClusters, RecoTracks const& tracks,
                  aod::JMcParticles const& particles, aod::JEMCTracks const& emctracks,
                  aod::JMcTrackLbs const& trackLabels)
  {
    if (trackLabels.size() != tracks.size() || trackLabels.offset() != tracks.offset()) {
      throw std::runtime_error("JMcTrackLbs must be row-aligned with JTracks; enable the JE MC track-label producer");
    }
    std::array<std::set<int64_t>, 3> trackingMatched;
    for (auto group = qaState.efficiencyMatched.begin(); group != qaState.efficiencyMatched.end(); ++group) {
      for (auto stage = group->begin(); stage != group->end(); ++stage) {
        stage->clear();
      }
    }
    std::set<int64_t> reconstructedMC, acceptedMC;
    for (auto const& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto mcCollision = collision.mcCollision_as<aod::JMcCollisions>();
      reconstructedMC.insert(mcCollision.globalIndex());
      const double weight = mcCollision.weight();
      auto eventClusters = clusters.sliceBy(clustersPerCollision, collision.globalIndex());
      auto rawEventClusters = rawClusters.sliceBy(clustersPerCollision, collision.globalIndex());
      auto eventTracks = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      processReco(collision, eventClusters, rawEventClusters, eventTracks, emctracks, 0, weight,
                  [&](auto const& c) {
                    return tagClusterForCollision(c, particles, mcCollision.globalIndex(), weight);
                  });
      if (qaState.recoEventAccepted) {
        acceptedMC.insert(mcCollision.globalIndex());
        collectTrackingMatches(eventTracks, trackLabels, particles, mcCollision.globalIndex(), emctracks, weight, trackingMatched);
      }
    }
    fillTrackingPerformance(particles, acceptedMC, trackingMatched);
    fillPhotonPerformance(particles, reconstructedMC, acceptedMC);
  }

  template <typename ParticleTable, typename Collision>
  void fillGenCollision(Collision const& collision, ParticleTable const& particles)
  {
    const double weight = collision.weight();
    if (!std::isfinite(weight)) {
      throw std::runtime_error("Non-finite generator weight");
    }
    qaState.genEvents->fill(0.5, weight);
    std::vector<ConeParticle> constituents;
    for (auto const& particle : particles) {
      if (!particle.isPhysicalPrimary() || !std::isfinite(particle.eta()) || !std::isfinite(particle.phi())) {
        continue;
      }
      const auto phi = phi02pi(particle.phi());
      qaState.allPrimaryMap->fill(particle.eta(), phi, weight);
      if (particle.pdgCode() == kGamma) {
        qaState.allPhotonMap->fill(particle.eta(), phi, weight);
      }
      if (!std::isfinite(particle.pt()) || particle.pt() < cfgtrkMinPt || std::abs(particle.eta()) > cfgtrkMaxEta) {
        continue;
      }
      auto* species = pdg->GetParticle(particle.pdgCode());
      if (!species || species->Charge() == 0.) {
        continue;
      }
      constituents.push_back({particle.pt(), particle.eta(), particle.phi()});
      qaState.genParticleMap->fill(particle.eta(), phi02pi(particle.phi()), weight);
    }
    // Apply the same truth fiducial and ancestry rules as the efficiency loop.
    for (auto const& photon : particles) {
      if (photon.pdgCode() != kGamma || !photon.isPhysicalPrimary() || photon.getGenStatusCode() < 0 ||
          !std::isfinite(photon.pt()) || !std::isfinite(photon.e()) || !std::isfinite(photon.phi()) ||
          !inCalo(photon.eta(), photon.phi())) {
        continue;
      }
      const auto phi = phi02pi(photon.phi());
      const auto result = classifyOriginal<ParticleTable>(photon, particles);
      if (isPrompt(result.origin)) {
        qaState.genPromptEta->fill(photon.eta(), weight);
        qaState.genPromptPhi->fill(phi, weight);
        measureIsolation(&qaState.genIsolation, photon.pt(), photon.eta(), phi, photon.e(), constituents, weight);
      }
      const bool direct = !result.directHardIndices.empty();
      const bool frag = !result.fragmentationHardIndices.empty();
      qaState.ancestryStatus->fill(direct ? (frag ? 3 : 1) : (frag ? 2 : (result.leptonBranch ? 4 : 0)), weight);
      // Groups overlap only for mixed ancestry. Each group counts the photon once.
      for (int g = 0; g < NPhotonGroups; ++g) {
        if (!(g == 0 ? direct : frag)) {
          continue;
        }
        qaState.genPhotonMap[g]->fill(photon.eta(), phi, weight);
        qaState.genPhi[g]->fill(phi, weight);
        qaState.genEta[g]->fill(photon.eta(), weight);
        qaState.genPt[g]->fill(photon.pt(), weight);
      }
    }
  }

  //******************************************
  // ANALYSIS HELPERS COMPLETE!!!
  // NOW FOR THE MAIN PROCESS FUNCTIONS
  //******************************************

  void processMCRec(MCCollisions const& collisions, FilteredMCClusters const& clusters,
                    aod::JClusters const& rawClusters,
                    RecoTracks const& tracks, TrackCandidates const&, aod::JMcParticles const& particles,
                    BcCandidates const&, aod::JEMCTracks const& emctracks, aod::JMcCollisions const&,
                    aod::JMcTrackLbs const& trackLabels)
  {
    fillMCReco(collisions, clusters, rawClusters, tracks, particles, emctracks, trackLabels);
  }

  PROCESS_SWITCH(PromptphotonQA, processMCRec, "Weighted reconstructed MC QA", false);

  void processData(DataCollisions::iterator const& collision, FilteredDataClusters const& clusters,
                   aod::JClusters const& rawClusters,
                   RecoTracks const& tracks, TrackCandidates const&, BcCandidates const&,
                   aod::JEMCTracks const& emctracks)
  {
    processReco(collision, clusters, rawClusters, tracks, emctracks, 1, 1., [](auto const&) { return TruthTag{}; });
  }
  PROCESS_SWITCH(PromptphotonQA, processData, "Real-data reconstructed QA", false);

  // Standard generator-only input has no reconstructed-table dependencies.
  void processMCGen(aod::McCollisions::iterator const& collision, aod::McParticles const& particles)
  {
    fillGenCollision(collision, particles);
  }
  PROCESS_SWITCH(PromptphotonQA, processMCGen, "Generator QA on standard MC truth only", true);

  // O2 groups the same JE collision rows used by processMCRec by their JE MC
  // label. Test for existence before event selection; split vertices count once.
  void processMCGenJE(aod::JMcCollisions::iterator const& collision,
                      soa::SmallGroups<MCCollisions> const& recoCollisions,
                      aod::JMcParticles const& particles)
  {
    if (cfgGenReqRec && recoCollisions.size() == 0) {
      return;
    }
    fillGenCollision(collision, particles);
  }
  PROCESS_SWITCH(PromptphotonQA, processMCGenJE, "Generator QA on JE truth with optional reconstructed counterpart", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<PromptphotonQA>(context)};
}
