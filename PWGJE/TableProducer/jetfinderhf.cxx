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
// Authors: Nima Zardoshti, Jochen Klein

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TDatabasePDG.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "PWGJE/DataModel/JetHF.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// NB: runDataProcessing.h must be included after customize!
#include "Framework/runDataProcessing.h"

template <typename JetTable, typename ConstituentTable, typename ConstituentSubTable>
struct JetFinderHFTask {
  Produces<JetTable> jetsTable;
  Produces<ConstituentTable> constituentsTable;
  Produces<ConstituentSubTable> constituentsSubTable;
  OutputObj<TH2F> h2JetPt{"h2_jet_pt"};
  OutputObj<TH2F> h2JetPhi{"h2_jet_phi"};
  OutputObj<TH2F> h2JetEta{"h2_jet_eta"};
  OutputObj<TH2F> h2JetNTracks{"h2_jet_ntracks"};
  OutputObj<TH1F> hJetPt{"h_jet_pt"};
  OutputObj<TH1F> hJetPhi{"h_jet_phi"};
  OutputObj<TH1F> hJetEta{"h_jet_eta"};
  OutputObj<TH1F> hJetNTracks{"h_jet_ntracks"};
  OutputObj<TH1F> hCandPt{"h_cand_pt"};

  Service<O2DatabasePDG> pdg;
  TrackSelection globalTracks;

  std::vector<fastjet::PseudoJet> jets;
  std::vector<fastjet::PseudoJet> inputParticles;
  JetFinder jetFinder;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.8, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.8, "maximum track eta"};
  Configurable<float> trackPhiMin{"trackPhiMin", -999, "minimum track phi"};
  Configurable<float> trackPhiMax{"trackPhiMax", 999, "maximum track phi"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  // cluster level configurables
  Configurable<std::string> clusterDefinitionS{"clusterDefinition", "kV3Default", "cluster definition to be selected, e.g. V3Default"};
  Configurable<float> clusterEtaMin{"clusterEtaMin", -0.7, "minimum cluster eta"}; // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterEtaMax{"clusterEtaMax", 0.7, "maximum cluster eta"};  // For ECMAL: |eta| < 0.7, phi = 1.40 - 3.26
  Configurable<float> clusterPhiMin{"clusterPhiMin", -999, "minimum cluster phi"};
  Configurable<float> clusterPhiMax{"clusterPhiMax", 999, "maximum cluster phi"};

  // HF candidate level configurables
  Configurable<std::string> candSpecie_s{"candSpecie_s", "D0", "options are D0, Lc, BPlus"};
  Configurable<std::string> candDecayChannel_s{"candDecayChannel_s", "default", "look up in task"};
  Configurable<float> candPtMin{"candPtMin", 0.0, "minimum candidate pT"};
  Configurable<float> candPtMax{"candPtMax", 100.0, "maximum candidate pT"};
  Configurable<float> candYMin{"candYMin", -0.8, "minimum candidate eta"};
  Configurable<float> candYMax{"candYMax", 0.8, "maximum candidate eta"};
  // HF candidiate selection configurables
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};
  Configurable<int> selectionFlagLcToPKPi{"selectionFlagLcToPKPi", 1, "Selection Flag for Lc->PKPi"};
  Configurable<int> selectionFlagLcToPiPK{"selectionFlagLcToPiPK", 1, "Selection Flag for Lc->PiPK"};
  Configurable<int> selectionFlagBPlus{"selectionFlagBPlus", 1, "Selection Flag for B+"};

  // jet level configurables
  Configurable<std::vector<double>> jetRadius{"jetRadius", {0.4}, "jet resolution parameters"};
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<int> jetTypeParticleLevel{"jetTypeParticleLevel", 1, "Type of stored jets. 0 = full, 1 = charged, 2 = neutral"};
  Configurable<int> jetAlgorithm{"jetAlgorithm", 2, "jet clustering algorithm. 0 = kT, 1 = C/A, 2 = Anti-kT"};
  Configurable<int> jetRecombScheme{"jetRecombScheme", 0, "jet recombination scheme. 0 = E-scheme, 1 = pT-scheme, 2 = pT2-scheme"};
  Configurable<float> jetGhostArea{"jetGhostArea", 0.005, "jet ghost area"};
  Configurable<int> ghostRepeat{"ghostRepeat", 1, "set to 0 to gain speed if you dont need area calculation"};
  Configurable<bool> DoRhoAreaSub{"DoRhoAreaSub", false, "do rho area subtraction"};
  Configurable<bool> DoConstSub{"DoConstSub", false, "do constituent subtraction"};

  int candPDG;
  int candDecay;

  void init(InitContext const&)
  {
    if (static_cast<std::string>(trackSelections) == "globalTracks") {
      globalTracks = getGlobalTrackSelection();
      globalTracks.SetEtaRange(trackEtaMin, trackEtaMax);
    }

    h2JetPt.setObject(new TH2F("h2_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                               100, 0., 100., 10, 0.05, 1.05));
    h2JetPhi.setObject(new TH2F("h2_jet_phi", "jet #phi;#phi",
                                80, -1., 7., 10, 0.05, 1.05));
    h2JetEta.setObject(new TH2F("h2_jet_eta", "jet #eta;#eta",
                                70, -0.7, 0.7, 10, 0.05, 1.05));
    h2JetNTracks.setObject(new TH2F("h2_jet_ntracks", "jet n;n constituents",
                                    30, 0., 30., 10, 0.05, 1.05));

    hJetPt.setObject(new TH1F("h_jet_pt", "jet p_{T};p_{T} (GeV/#it{c})",
                              100, 0., 100.));
    hJetPhi.setObject(new TH1F("h_jet_phi", "jet #phi; #phi",
                               140, -7.0, 7.0));
    hJetEta.setObject(new TH1F("h_jet_eta", "jet #eta; #eta",
                               30, -1.5, 1.5));
    hJetNTracks.setObject(new TH1F("h_jet_ntracks", "jet N tracks ; N tracks",
                                   150, -0.5, 99.5));
    hCandPt.setObject(new TH1F("h_cand_pt", "jet p_{T,cand};p_{T,cand} (GeV/#it{c})",
                               100, 0., 100.));

    if (DoRhoAreaSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::rhoAreaSub);
    }
    if (DoConstSub) {
      jetFinder.setBkgSubMode(JetFinder::BkgSubMode::constSub);
    }

    jetFinder.etaMin = trackEtaMin;
    jetFinder.etaMax = trackEtaMax;
    jetFinder.jetPtMin = jetPtMin;
    jetFinder.jetPtMax = jetPtMax;
    jetFinder.algorithm = static_cast<fastjet::JetAlgorithm>(static_cast<int>(jetAlgorithm));
    jetFinder.recombScheme = static_cast<fastjet::RecombinationScheme>(static_cast<int>(jetRecombScheme));
    jetFinder.ghostArea = jetGhostArea;
    jetFinder.ghostRepeatN = ghostRepeat;

    auto candSpecie = static_cast<std::string>(candSpecie_s);
    auto candDecayChannel = static_cast<std::string>(candDecayChannel_s);
    if (candSpecie == "D0") {
      candPDG = static_cast<int>(pdg::Code::kD0);
      candDecay = static_cast<int>(aod::hf_cand_2prong::DecayType::D0ToPiK);
    }
    if (candSpecie == "BPlus") {
      candPDG = static_cast<int>(pdg::Code::kBPlus);
      candDecay = static_cast<int>(aod::hf_cand_bplus::DecayType::BplusToD0Pi);
    }
    if (candSpecie == "Lc") {
      candPDG = static_cast<int>(pdg::Code::kLambdaCPlus);
      candDecay = static_cast<int>(aod::hf_cand_3prong::DecayType::LcToPKPi);
    }
    if (candSpecie == "JPsi") {
      candPDG = static_cast<int>(pdg::Code::kJPsi);
      if (candDecayChannel == "default" || candDecayChannel == "ee") {
        candDecay = static_cast<int>(aod::hf_cand_2prong::DecayType::JpsiToEE);
      }
      if (candDecayChannel == "mumu") {
        candDecay = static_cast<int>(aod::hf_cand_2prong::DecayType::JpsiToMuMu);
      }
    }
  }

  using JetTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;
  using JetParticles2Prong = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand2ProngMcGen>>;
  using JetParticles3Prong = soa::Filtered<soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>>;
  using JetParticlesBPlus = soa::Filtered<soa::Join<aod::McParticles, aod::HfCandBplusMcGen>>;

  o2::aod::EMCALClusterDefinition clusterDefinition = o2::aod::emcalcluster::getClusterDefinitionFromString(clusterDefinitionS.value);
  Filter collisionFilter = (nabs(aod::collision::posZ) < vertexZCut);
  Filter trackCuts = (aod::track::pt >= trackPtMin && aod::track::pt < trackPtMax && aod::track::eta > trackEtaMin && aod::track::eta < trackEtaMax && aod::track::phi >= trackPhiMin && aod::track::phi <= trackPhiMax);
  Filter partCuts = (aod::mcparticle::pt >= trackPtMin && aod::mcparticle::pt < trackPtMax);
  Filter clusterFilter = (o2::aod::emcalcluster::definition == static_cast<int>(clusterDefinition) && aod::emcalcluster::eta > clusterEtaMin && aod::emcalcluster::eta < clusterEtaMax && aod::emcalcluster::phi >= clusterPhiMin && aod::emcalcluster::phi <= clusterPhiMax);
  // Filter candidateCuts = (aod::hf_cand::Pt >= candPtMin && aod::hf_cand::Pt < candPtMax); FIXME: why wont this work?
  Filter candidateCutsD0 = (aod::hf_sel_candidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_sel_candidate_d0::isSelD0bar >= selectionFlagD0bar);
  Filter candidateCutsLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLcToPKPi || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLcToPiPK);
  Filter candidateCutsBPlus = (aod::hf_sel_candidate_bplus::isSelBplusToD0Pi >= selectionFlagBPlus);

  // function that takes any generic candidate, performs selections and adds the candidate to the fastjet list
  template <typename T>
  bool processCandidate(T const& candidate)
  {
    if (candidate.y(RecoDecay::getMassPDG(candPDG)) < candYMin || candidate.y(RecoDecay::getMassPDG(candPDG)) > candYMax) {
      return false;
    }
    if (candidate.pt() < candPtMin || candidate.pt() >= candPtMax) {
      return false;
    }
    FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candPDG));
    return true;
  }

  // function that checks the MC status of a candidate and then calls the function to processCandidates
  template <typename T>
  bool processCandidateMC(T const& candidate)
  {
    if (!(std::abs(candidate.flagMcMatchRec()) == 1 << candDecay)) {
      return false;
    }
    return processCandidate(candidate);
  }

  // function that performs track selections on each track
  template <typename T>
  bool processTrackSelection(T const& track)
  {
    if (static_cast<std::string>(trackSelections) == "globalTracks" && !globalTracks.IsSelected(track)) {
      return false;
    } else if (static_cast<std::string>(trackSelections) == "QualityTracks" && !track.isQualityTrack()) {
      return false;
    } else {
      return true;
    }
  }

  // function that adds tracks to the fastjet list, removing daughters of 2Prong candidates
  template <typename T, typename U>
  void processTracks2Prong(T const& tracks, U const& candidate)
  {
    for (auto& track : tracks) {
      if (!processTrackSelection(track)) {
        continue;
      }
      if (candidate.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
        continue;
      }
      FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
    }
  }
  // function that adds tracks to the fastjet list, removing daughters of 3Prong candidates
  template <typename T, typename U>
  void processTracks3Prong(T const& tracks, U const& candidate)
  {
    for (auto& track : tracks) {
      if (!processTrackSelection(track)) {
        continue;
      }
      if (candidate.template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.template prong2_as<JetTracks>().globalIndex() == track.globalIndex()) {
        continue;
      }
      FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
    }
  }
  // function that adds tracks to the fastjet list, removing daughters of B+ candidates
  template <typename T, typename U>
  void processTracksBPlus(T const& tracks, U const& candidate)
  {
    for (auto& track : tracks) {
      if (!processTrackSelection(track)) {
        continue;
      }
      if (candidate.template prong0_as<aod::HfCand2Prong>().template prong0_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.template prong0_as<aod::HfCand2Prong>().template prong1_as<JetTracks>().globalIndex() == track.globalIndex() || candidate.template prong1_as<JetTracks>().globalIndex() == track.globalIndex()) {
        continue;
      }
      FastJetUtilities::fillTracks(track, inputParticles, track.globalIndex());
    }
  }

  // function that adds clusters to the fastjet list
  template <typename T>
  void processClusters(T const& clusters)
  {
    for (auto& cluster : *clusters) {
      // add cluster selections
      FastJetUtilities::fillClusters(cluster, inputParticles, cluster.globalIndex());
    }
  }

  // function that calls the jet finding and fills the relevant tables
  template <typename T>
  void jetFinding(T const& collision)
  {
    auto candidatepT = 0.0;
    auto jetRValues = static_cast<std::vector<double>>(jetRadius);
    for (auto R : jetRValues) {
      jetFinder.jetR = R;
      jets.clear();
      fastjet::ClusterSequenceArea clusterSeq(jetFinder.findJets(inputParticles, jets));
      for (const auto& jet : jets) {
        bool isHFJet = false;
        for (const auto& constituent : jet.constituents()) {
          if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
            isHFJet = true;
            candidatepT = constituent.pt();
            break;
          }
        }
        if (!isHFJet) {
          continue;
        }
        std::vector<int> trackconst;
        std::vector<int> candconst;
        std::vector<int> clusterconst;

        jetsTable(collision, jet.pt(), jet.eta(), jet.phi(),
                  jet.E(), jet.m(), jet.area(), std::round(R * 100));
        for (const auto& constituent : sorted_by_pt(jet.constituents())) {
          // need to add seperate thing for constituent subtraction
          if (DoConstSub) { // FIXME: needs to be addressed in Haadi's PR
            constituentsSubTable(jetsTable.lastIndex(), constituent.pt(), constituent.eta(), constituent.phi(),
                                 constituent.E(), constituent.m(), constituent.user_index());
          }

          if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::track)) {
            trackconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
          }
          if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::cluster)) {
            clusterconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
          }
          if (constituent.template user_info<FastJetUtilities::fastjet_user_info>().getStatus() == static_cast<int>(JetConstituentStatus::candidateHF)) {
            candconst.push_back(constituent.template user_info<FastJetUtilities::fastjet_user_info>().getIndex());
          }
        }
        constituentsTable(jetsTable.lastIndex(), trackconst, clusterconst, candconst);
        h2JetPt->Fill(jet.pt(), R);
        h2JetPhi->Fill(jet.phi(), R);
        h2JetEta->Fill(jet.rap(), R);
        h2JetNTracks->Fill(jet.constituents().size(), R);
        hJetPt->Fill(jet.pt());
        hJetPhi->Fill(jet.phi());
        hJetEta->Fill(jet.rap());
        hJetNTracks->Fill(jet.constituents().size());
        hCandPt->Fill(candidatepT);
        break;
      }
    }
  }

  // function that processes data for all 2Prong candidates
  template <typename T, typename U, typename M>
  void processData2Prong(T const& collision, U const& tracks, M const& candidates)
  {
    if (!collision.sel8())
      return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidate(candidate)) {
        continue;
      }
      processTracks2Prong(tracks, candidate);
      jetFinding(collision);
    }
  }

  // function that processes MC det for all 2Prong candidates
  template <typename T, typename U, typename M>
  void processMCD2Prong(T const& collision, U const& tracks, M const& candidates)
  {
    // if (!collision.sel8())
    // return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidateMC(candidate)) {
        continue;
      }
      processTracks2Prong(tracks, candidate);
      jetFinding(collision);
    }
  }

  // function that processes data for all 3Prong candidates
  template <typename T, typename U, typename M>
  void processData3Prong(T const& collision, U const& tracks, M const& candidates)
  {
    if (!collision.sel8())
      return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidate(candidate)) {
        continue;
      }
      processTracks3Prong(tracks, candidate);
      jetFinding(collision);
    }
  }
  // function that processes MC det for all 3Prong candidates
  template <typename T, typename U, typename M>
  void processMCD3Prong(T const& collision, U const& tracks, M const& candidates)
  {
    // if (!collision.sel8())
    // return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidateMC(candidate)) {
        continue;
      }
      processTracks3Prong(tracks, candidate);
      jetFinding(collision);
    }
  }

  // function that processes data for all B+ candidates
  template <typename T, typename U, typename M>
  void processDataBPlus(T const& collision, U const& tracks, M const& candidates)
  {
    if (!collision.sel8())
      return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidate(candidate)) {
        continue;
      }
      processTracksBPlus(tracks, candidate);
      jetFinding(collision);
    }
  }

  // function that processes MC det for all B+ candidates
  template <typename T, typename U, typename M>
  void processMCDBPlus(T const& collision, U const& tracks, M const& candidates)
  {
    // if (!collision.sel8())
    // return;

    for (auto& candidate : candidates) {
      inputParticles.clear();
      if (!processCandidateMC(candidate)) {
        continue;
      }
      processTracksBPlus(tracks, candidate);
      jetFinding(collision);
    }
  }

  // function that checks if a candidate has any daughters that need to be removed from the event at gen level
  template <typename T, typename U> // CHECK if this works for everything
  bool checkDaughters(T const& particle, U const globalIndex)
  {
    for (auto daughter : particle.template daughters_as<aod::McParticles>()) {
      if (daughter.globalIndex() == globalIndex) {
        return true;
      }
      if (checkDaughters(daughter, globalIndex)) {
        return true;
      }
    }
    return false;
  }

  // function that generalically processes gen level events
  template <typename T, typename U, typename M>
  void processParticlesMCGen(T const& collision, U const& particles, M& candidates)
  {

    for (auto const& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) & (1 << candDecay)) {
        auto particleY = RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode()));
        if (particleY < candYMin || particleY > candYMax) {
          continue;
        }
        if (particle.pt() < candPtMin || particle.pt() >= candPtMax) {
          continue;
        }
        candidates.push_back(particle);
      }
    }
    for (auto& candidate : candidates) {
      inputParticles.clear();
      for (auto& particle : particles) {
        // TODO: can we do this through the filter?
        if (particle.eta() < trackEtaMin || particle.eta() > trackEtaMax) {
          continue;
        }
        if (particle.getGenStatusCode() != 1) { // CHECK : Does this include HF hadrons that decay?
          continue;
        }
        auto pdgParticle = pdg->GetParticle(particle.pdgCode());
        auto pdgCharge = pdgParticle ? std::abs(pdgParticle->Charge()) : -1.0;
        if (jetTypeParticleLevel == static_cast<int>(JetType::charged) && pdgCharge < 3.0) {
          continue;
        }
        if (jetTypeParticleLevel == static_cast<int>(JetType::neutral) && pdgCharge != 0.0) {
          continue;
        }
        if (checkDaughters(candidate, particle.globalIndex())) {
          continue;
        }
        FastJetUtilities::fillTracks(particle, inputParticles, particle.globalIndex(), static_cast<int>(JetConstituentStatus::track), RecoDecay::getMassPDG(particle.pdgCode()));
      }
      FastJetUtilities::fillTracks(candidate, inputParticles, candidate.globalIndex(), static_cast<int>(JetConstituentStatus::candidateHF), RecoDecay::getMassPDG(candidate.pdgCode()));
      jetFinding(collision);
    }
  }

  // check if type JetParticles2Prong can be templated. then you can just use one function everywhere
  // function that is called for gen level events with 2 prong candidates
  template <typename T, typename U>
  void processHFGen2Prong(T const& collision, U const& particles)
  {
    jets.clear();
    inputParticles.clear();
    LOG(debug) << "Per Event MCP";
    std::vector<JetParticles2Prong::iterator> candidates;
    candidates.clear();
    processParticlesMCGen(collision, particles, candidates);
  }
  // function that is called for gen level events with 3 prong candidates
  template <typename T, typename U>
  void processHFGen3Prong(T const& collision, U const& particles)
  {
    LOG(debug) << "Per Event MCP";
    std::vector<JetParticles3Prong::iterator> candidates;
    processParticlesMCGen(collision, particles, candidates);
  }
  // function that is called for gen level events with B+ candidates
  template <typename T, typename U>
  void processHFGenBPlus(T const& collision, U const& particles)
  {
    LOG(debug) << "Per Event MCP";
    std::vector<JetParticlesBPlus::iterator> candidates;
    processParticlesMCGen(collision, particles, candidates);
  }

  void processDummy(aod::Collisions const& collision)
  {
  }
  PROCESS_SWITCH(JetFinderHFTask, processDummy, "Dummy process function turned on by default", true);

  void processD0Data(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                     JetTracks const& tracks,
                     soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0>> const& candidates)
  {
    processData2Prong(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processD0Data, "D0 jet finding on data", false);

  void processD0MCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                    JetTracks const& tracks,
                    soa::Filtered<soa::Join<aod::HfCand2Prong, aod::HfSelD0, aod::HfCand2ProngMcRec>> const& candidates)
  {
    processMCD2Prong(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processD0MCD, "D0 finding on MC detector level", false);

  void processBPlusData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                        JetTracks const& tracks,
                        soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi>> const& candidates,
                        aod::HfCand2Prong const& HFdaughters)
  {
    processDataBPlus(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processBPlusData, "B+ jet finding on data", false);

  void processBPlusMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                       JetTracks const& tracks,
                       soa::Filtered<soa::Join<aod::HfCandBplus, aod::HfSelBplusToD0Pi, aod::HfCandBplusMcRec>> const& candidates,
                       aod::HfCand2Prong const& HFdaughters)
  {
    processMCDBPlus(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processBPlusMCD, "B+ finding on MC detector level", false);

  void processLcData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
                     JetTracks const& tracks,
                     soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>> const& candidates)
  {
    processData3Prong(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processLcData, "Lc jet finding on data", false);

  void processLcMCD(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                    JetTracks const& tracks,
                    soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>> const& candidates)
  {
    processMCD3Prong(collision, tracks, candidates);
  }
  PROCESS_SWITCH(JetFinderHFTask, processLcMCD, "Lc finding on MC detector level", false);

  void process2ProngMCP(aod::McCollision const& collision,
                        JetParticles2Prong const& particles)
  {
    processHFGen2Prong(collision, particles);
  }
  PROCESS_SWITCH(JetFinderHFTask, process2ProngMCP, "2-prong HF jet finding on MC particle level", false);

  void process3ProngMCP(aod::McCollision const& collision,
                        JetParticles3Prong const& particles)
  {
    processHFGen3Prong(collision, particles);
  }
  PROCESS_SWITCH(JetFinderHFTask, process3ProngMCP, "3-prong HF jet finding on MC particle level", false);

  void processBPlusMCP(aod::McCollision const& collision,
                       JetParticlesBPlus const& particles)
  {
    processHFGenBPlus(collision, particles);
  }
  PROCESS_SWITCH(JetFinderHFTask, processBPlusMCP, "B+ HF jet finding on MC particle level", false);
};

using JetFinderD0 = JetFinderHFTask<o2::aod::D0Jets, o2::aod::D0JetConstituents, o2::aod::D0JetConstituentsSub>;
using MCDetectorLevelJetFinderD0 = JetFinderHFTask<o2::aod::MCDetectorLevelD0Jets, o2::aod::MCDetectorLevelD0JetConstituents, o2::aod::MCDetectorLevelD0JetConstituentsSub>;
using MCParticleLevelJetFinderD0 = JetFinderHFTask<o2::aod::MCParticleLevelD0Jets, o2::aod::MCParticleLevelD0JetConstituents, o2::aod::MCParticleLevelD0JetConstituentsSub>;

using JetFinderBPlus = JetFinderHFTask<o2::aod::BPlusJets, o2::aod::BPlusJetConstituents, o2::aod::BPlusJetConstituentsSub>;
using MCDetectorLevelJetFinderBPlus = JetFinderHFTask<o2::aod::MCDetectorLevelBPlusJets, o2::aod::MCDetectorLevelBPlusJetConstituents, o2::aod::MCDetectorLevelBPlusJetConstituentsSub>;
using MCParticleLevelJetFinderBPlus = JetFinderHFTask<o2::aod::MCParticleLevelBPlusJets, o2::aod::MCParticleLevelBPlusJetConstituents, o2::aod::MCParticleLevelBPlusJetConstituentsSub>;

using JetFinderLc = JetFinderHFTask<o2::aod::LcJets, o2::aod::LcJetConstituents, o2::aod::LcJetConstituentsSub>;
using MCDetectorLevelJetFinderLc = JetFinderHFTask<o2::aod::MCDetectorLevelLcJets, o2::aod::MCDetectorLevelLcJetConstituents, o2::aod::MCDetectorLevelLcJetConstituentsSub>;
using MCParticleLevelJetFinderLc = JetFinderHFTask<o2::aod::MCParticleLevelLcJets, o2::aod::MCParticleLevelLcJetConstituents, o2::aod::MCParticleLevelLcJetConstituentsSub>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<JetFinderD0>(cfgc,
                                                    SetDefaultProcesses{},
                                                    TaskName{"jet-finder-D0-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderD0>(cfgc,
                                                                   SetDefaultProcesses{},
                                                                   TaskName{"jet-finder-D0-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderD0>(cfgc,
                                                                   SetDefaultProcesses{},
                                                                   TaskName{"jet-finder-D0-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderBPlus>(cfgc,
                                                       SetDefaultProcesses{},
                                                       TaskName{"jet-finder-BPlus-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderBPlus>(cfgc,
                                                                      SetDefaultProcesses{},
                                                                      TaskName{"jet-finder-BPlus-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderBPlus>(cfgc,
                                                                      SetDefaultProcesses{},
                                                                      TaskName{"jet-finder-BPlus-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<JetFinderLc>(cfgc,
                                                    SetDefaultProcesses{},
                                                    TaskName{"jet-finder-Lc-data"}));

  tasks.emplace_back(adaptAnalysisTask<MCDetectorLevelJetFinderLc>(cfgc,
                                                                   SetDefaultProcesses{},
                                                                   TaskName{"jet-finder-Lc-mcd"}));

  tasks.emplace_back(adaptAnalysisTask<MCParticleLevelJetFinderLc>(cfgc,
                                                                   SetDefaultProcesses{},
                                                                   TaskName{"jet-finder-Lc-mcp"}));

  return WorkflowSpec{tasks};
}
