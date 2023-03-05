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

// jet finder task
//
// Author: Jochen Klein, Nima Zardoshti, Raymond Ehlers

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetFinder.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec jetData = {
    "jetInputData",
    VariantType::String,
    "",
    {"Jet input data type. Options include Data, MCParticleLevel, MCDetectorLevel, and HybridIntermediate."},
  };
  workflowOptions.push_back(jetData);
  ConfigParamSpec jetType = {
    "jetType",
    VariantType::Int,
    1,
    {"Type of jet. Options include full = 0, charged = 1, neutral = 2."},
  };
  workflowOptions.push_back(jetType);
}

#include "Framework/runDataProcessing.h"

enum class JetType_t {
  full = 0,
  charged = 1,
  neutral = 2,
};

template <typename JetTable, typename TrackConstituentTable, typename ClusterConstituentTable, typename ConstituentSubTable>
struct JetFinderTask {
  Produces<JetTable> jetsTable;
  Produces<TrackConstituentTable> trackConstituentsTable;
  Produces<ClusterConstituentTable> clusterConstituentsTable;
  Produces<ConstituentSubTable> constituentsSubTable;
  OutputObj<TH2F> hJetPt{"h_jet_pt"};
  OutputObj<TH2F> hJetPhi{"h_jet_phi"};
  OutputObj<TH2F> hJetEta{"h_jet_eta"};
  OutputObj<TH2F> hJetN{"h_jet_n"};

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> trackPtCut{"trackPtCut", 0.1, "minimum constituent pT"};
  Configurable<float> trackEtaCut{"trackEtaCut", 0.9, "constituent eta cut"};
  Configurable<float> trackPhiMinCut{"trackPhiMinCut", -999, "track min phi cut"};
  Configurable<float> trackPhiMaxCut{"trackPhiMaxCut", 999, "track max phi cut"};
  Configurable<float> clusterEtaCut{"clusterEtaCut", 0.7, "cluster eta cut"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMinCut{"clusterPhiMinCut", -999, "cluster min phi cut"};
  Configurable<float> clusterPhiMaxCut{"clusterPhiMaxCut", 999, "cluster max phi cut"};
  Configurable<bool> DoRhoAreaSub{"DoRhoAreaSub", false, "do rho area subtraction"};
  Configurable<bool> DoConstSub{"DoConstSub", false, "do constituent subtraction"};
  Configurable<float> jetPtMin{"jetPtMin", 10.0, "minimum jet pT"};
  Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  Configurable<std::vector<double>> jetR{"jetR", {0.4}, "jet resolution parameters"};
  // FIXME: This should be named jetType. However, as of Aug 2021, it doesn't appear possible
  //        to set both global and task level options. This should be resolved when workflow
  //        level customization is available
  Configurable<int> jetType2{"jetType2", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<std::string> mClusterDefinition{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};

  Filter collisionFilter = nabs(aod::collision::posZ) < vertexZCut;
  // Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);
  Filter trackFilter = (nabs(aod::track::eta) < trackEtaCut) && (aod::track::phi > trackPhiMinCut) && (aod::track::phi < trackPhiMaxCut) && (requireGlobalTrackInFilter()) && (aod::track::pt > trackPtCut);

  o2::aod::EMCALClusterDefinition clusDef = o2::aod::emcalcluster::getClusterDefinitionFromString(mClusterDefinition.value);
  // Filter clusterDefinitionSelection = o2::aod::emcalcluster::definition == static_cast<int>(clusDef);
  Filter clusterFilter = (o2::aod::emcalcluster::definition == static_cast<int>(clusDef)) && (nabs(aod::emcalcluster::eta) < clusterEtaCut) && (aod::emcalcluster::phi > clusterPhiMinCut) && (aod::emcalcluster::phi < clusterPhiMaxCut);

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder; // should be a configurable but for now this cant be changed on hyperloop
  // FIXME: Once configurables support enum, ideally we can
  JetType_t _jetType;

  void init(InitContext const&)
  {
    hJetPt.setObject(new TH2F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100., 10, 0.05, 1.05));
    hJetPhi.setObject(new TH2F("h_jet_phi", "jet #phi;#phi",
                               80, -1., 7., 10, 0.05, 1.05));
    hJetEta.setObject(new TH2F("h_jet_eta", "jet #eta;#eta",
                               70, -0.7, 0.7, 10, 0.05, 1.05));
    hJetN.setObject(new TH2F("h_jet_n", "jet n;n constituents",
                             30, 0., 30., 10, 0.05, 1.05));
    if (DoRhoAreaSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::rhoAreaSub);
    }
    if (DoConstSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::constSub);
    }
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.ghostRepeatN = ghostRepeat;
  }

  template <typename T>
  bool processInit(T const& collision)
  {
    // FIXME: Reenable event selection when working (disabled Aug 2021)
    /*
    if (!collision.alias()[kINT7]) {
      return false; //remove hard code
    }
    if (!collision.sel7()) {
      return false; //remove hard code
    }
    */

    jets.clear();
    inputParticles.clear();

    return true;
  }

  template <typename T>
  void processImplementation(T const& collision)
  {
    LOG(debug) << "Process Implementation";
    // NOTE: Can't just iterate directly - we have to cast first
    auto jetRValues = static_cast<std::vector<double>>(jetR);
    for (auto R : jetRValues) {
      // Update jet finder R and find jets
      jetFinder.jetR = R;
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));

      for (const auto& jet : jets) {
        jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                  jet.E(), jet.m(), jet.area(), std::round(R * 100));
        hJetPt->Fill(jet.pt(), R);
        hJetPhi->Fill(jet.phi(), R);
        hJetEta->Fill(jet.eta(), R);
        hJetN->Fill(jet.constituents().size(), R);
        for (const auto& constituent : jet.constituents()) { // event or jetwise
          if (DoConstSub) {
            // Since we're copying the consituents, we can combine the tracks and clusters together
            // We only have to keep the uncopied versions separated due to technical constraints.
            constituentsSubTable(jetsTable.lastIndex(), constituent.pt(), constituent.eta(), constituent.phi(),
                                 constituent.E(), constituent.m(), constituent.user_index());
          }
          if (constituent.user_index() < 0) {
            // Cluster
            // -1 to account for the convention of negative indices for clusters.
            clusterConstituentsTable(jetsTable.lastIndex(), -1 * constituent.user_index());
          } else {
            // Tracks
            trackConstituentsTable(jetsTable.lastIndex(), constituent.user_index());
          }
        }
      }
    }
  }

  void processParticleLevel(aod::McCollision const& collision, aod::McParticles const& particles)
  {
    // Setup
    // As of June 2021, I don't think enums are supported as configurables, so we have to handle the conversion here.
    // TODO: Double cast is to work around conversion failure.
    _jetType = static_cast<JetType_t>(static_cast<int>(jetType2));

    // Initialziation and event selection
    // TODO: MC event selection?
    jets.clear();
    inputParticles.clear();

    // As of June 2021, how best to check for charged particles? It doesn't seem to be in
    // the McParticles table, so for now we select by PID.
    // charged hadron (e-, mu-, pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-)
    // We define full jets as taking everything, and neutral jets as not charged.
    std::vector<unsigned int> selectedPIDs{11, 13, 211, 321, 2212, 3222, 3112, 3312, 3334};
    auto selectedBegin = selectedPIDs.begin();
    auto selectedEnd = selectedPIDs.end();
    for (auto& particle : particles) {
      if (_jetType != JetType_t::full) {
        bool foundChargedParticle = (std::find(selectedBegin, selectedEnd, particle.pdgCode()) == selectedEnd);
        if (_jetType == JetType_t::charged && foundChargedParticle == false) {
          continue;
        }
        if (_jetType == JetType_t::neutral && foundChargedParticle == true) {
          continue;
        }
      }

      inputParticles.emplace_back(
        fastjet::PseudoJet(
          particle.px(), particle.py(), particle.pz(), particle.e()));
      inputParticles.back().set_user_index(particle.globalIndex());
    }

    processImplementation(collision);
  }

  PROCESS_SWITCH(JetFinderTask, processParticleLevel, "Particle level jet finding", false);

  template <typename T, typename U>
  void processData(T const& collision, U const& tracks, aod::EMCALClusters const* clusters = nullptr)
  {
    LOG(debug) << "Process data!";
    // Setup
    // As of June 2021, I don't think enums are supported as configurables, so we have to handle the conversion here.
    // FIXME: Double cast is to work around conversion failure.
    _jetType = static_cast<JetType_t>(static_cast<int>(jetType2));

    // Initialziation and event selection
    bool accepted = processInit(collision);
    if (!accepted) {
      return;
    }
    LOG(debug) << "Accepted event!";

    if (_jetType == JetType_t::full || _jetType == JetType_t::charged) {
      for (auto& track : tracks) {
        fillConstituents(track, inputParticles);
        inputParticles.back().set_user_index(track.globalIndex());
      }
    }
    if (_jetType == JetType_t::full || _jetType == JetType_t::neutral) {
      if (clusters) {
        for (auto& cluster : *clusters) {
          // The right thing to do here would be to fully calculate the momentum correcting for the vertex position.
          // However, it's not clear that this functionality exists yet (21 June 2021)
          double pt = cluster.energy() / std::cosh(cluster.eta());
          inputParticles.emplace_back(
            fastjet::PseudoJet(
              pt * std::cos(cluster.phi()),
              pt * std::sin(cluster.phi()),
              pt * std::sinh(cluster.eta()),
              cluster.energy()));
          // Clusters are denoted with negative indices.
          inputParticles.back().set_user_index(-1 * cluster.globalIndex());
        }
      } else {
        throw std::runtime_error("Requested clusters, but they're not available!");
      }
    }

    processImplementation(collision);
  }

  void processDataCharged(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                          soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    LOG(debug) << "Process data charged!";
    processData(collision, tracks);
  }

  PROCESS_SWITCH(JetFinderTask, processDataCharged, "Data jet finding for charged jets", true);

  void processDataFull(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                       soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks,
                       soa::Filtered<aod::EMCALClusters> const& clusters)
  {
    LOG(debug) << "Process data full!";
    processData(collision, tracks, &clusters);
  }

  PROCESS_SWITCH(JetFinderTask, processDataFull, "Data jet finding for full and neutral jets", false);
};

using JetFinderData = JetFinderTask<o2::aod::Jets, o2::aod::JetTrackConstituents, o2::aod::JetClusterConstituents, o2::aod::JetConstituentsSub>;
using JetFinderMCParticleLevel = JetFinderTask<o2::aod::MCParticleLevelJets, o2::aod::MCParticleLevelJetTrackConstituents, o2::aod::MCParticleLevelJetClusterConstituents, o2::aod::MCParticleLevelJetConstituentsSub>;
using JetFinderMCDetectorLevel = JetFinderTask<o2::aod::MCDetectorLevelJets, o2::aod::MCDetectorLevelJetTrackConstituents, o2::aod::MCDetectorLevelJetClusterConstituents, o2::aod::MCDetectorLevelJetConstituentsSub>;
using JetFinderHybridIntermediate = JetFinderTask<o2::aod::HybridIntermediateJets, o2::aod::HybridIntermediateJetTrackConstituents, o2::aod::HybridIntermediateJetClusterConstituents, o2::aod::HybridIntermediateJetConstituentsSub>;

enum class JetInputData_t {
  Data,
  MCParticleLevel,
  MCDetectorLevel,
  HybridIntermediate
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Determine the jet input data types.
  // NOTE: It's possible to create multiple jet finders by passing multiple jet input data types in the string.
  //       Each type is separated by a comma.
  // FIXME: String validation and normalization. There's too much room for error with the strings. (Is there a better way to do this?)
  auto jetInputData = cfgc.options().get<std::string>("jetInputData");
  const std::map<std::string, JetInputData_t> jetInputDataTypes = {
    {"Data", JetInputData_t::Data},
    {"MCParticleLevel", JetInputData_t::MCParticleLevel},
    {"MCDetectorLevel", JetInputData_t::MCDetectorLevel},
    {"HybridIntermediate", JetInputData_t::HybridIntermediate},
    {"", JetInputData_t::Data}, // Default to data
  };
  // Tokenize using stringstream
  std::vector<std::string> jetInputDataOptions;
  std::stringstream ss;
  ss << jetInputData;
  while (ss.good()) {
    std::string substring;
    getline(ss, substring, ',');
    jetInputDataOptions.push_back(substring);
  }
  // Jet type
  auto jetType = static_cast<JetType_t>(cfgc.options().get<int>("jetType"));
  std::vector<o2::framework::DataProcessorSpec> tasks;
  for (auto opt : jetInputDataOptions) {
    auto jetData = jetInputDataTypes.at(opt);
    switch (jetData) {
      case JetInputData_t::MCParticleLevel:
        // We don't need a switch on jet type here because the arguments for particle level jets are the same
        // for charged and full jets
        tasks.emplace_back(
          adaptAnalysisTask<JetFinderMCParticleLevel>(cfgc,
                                                      SetDefaultProcesses{{{"processParticleLevel", true}, {"processDataCharged", false}, {"processDataFull", false}}}, TaskName{"jet-finder-MC"}));
        break;
      case JetInputData_t::MCDetectorLevel:
        if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderMCDetectorLevel>(cfgc,
                                                        SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-MC-detector-level"}));
        } else {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderMCDetectorLevel>(cfgc, TaskName{"jet-finder-MC-detector-level"}));
        }
        break;
      case JetInputData_t::HybridIntermediate:
        if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc,
                                                           SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-hybrid-intermedaite"}));
        } else {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderHybridIntermediate>(cfgc, TaskName{"jet-finder-hybrid-intermedaite"}));
        }
        break;
      case JetInputData_t::Data: // Intentionally fall through to the default for data.
      default:
        if (jetType == JetType_t::full || jetType == JetType_t::neutral) {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderData>(cfgc,
                                             SetDefaultProcesses{{{"processParticleLevel", false}, {"processDataCharged", false}, {"processDataFull", true}}}, TaskName{"jet-finder-data"}));
        } else {
          tasks.emplace_back(
            adaptAnalysisTask<JetFinderData>(cfgc, TaskName{"jet-finder-data"}));
        }
        break;
    }
  }
  return WorkflowSpec{tasks};
}
