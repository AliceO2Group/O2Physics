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

/// \file jetmatching.cxx
/// \brief jet nSubjettiness calculation task

/// \author Aimeric Landou <aimeric.landou@cern.ch>

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/Logger.h"
#include "Framework/HistogramRegistry.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/FastJetUtilities.h"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/AxesDefinition.hh"
#include "fastjet/contrib/MeasureDefinition.hh"
#include "fastjet/contrib/SoftDrop.hh"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
using Tracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

template <typename CollisionTable, typename JetTable, typename JetEventWeightTable, typename TrackTable>
struct nSubJettiness {

  HistogramRegistry registry{
    "registry",
    {
      {"hErrorControl", "hErrorControl", {HistType::kTH1F, {{2, -999.0f, -997.0f, ""}}}},
    },
  };

  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};

  Configurable<float> SD_beta{"SD_beta", 0.0, "SoftDrop beta"};
  Configurable<float> SD_z_cut{"SD_z_cut", 0.10, "SoftDrop z cut"};

  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 100.0f}, ""};
  ConfigurableAxis DeltaRBinning{"Delta-R-binning", {50, 0.0, 0.5}, ""};
  ConfigurableAxis NSubRatioBinning{"NSub-Ratio-binning", {50, 0.0f, 1.2f}, ""};

  std::vector<fastjet::PseudoJet> jetConstituents;
  fastjet::ClusterSequence clusterSeq_pseudoJet;

  void init(InitContext const&)
  {
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec DeltaRAxis = {DeltaRBinning, "#Delta R"};
    AxisSpec NSubRatioAxis = {NSubRatioBinning, "#tau_{2}/#tau_{1}"};

    registry.add("hNSubRatio21_Kt", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_Kt"))->Sumw2();
    registry.add("hNSubRatio21_CA", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_CA"))->Sumw2();
    registry.add("hNSubRatio21_SD", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_SD"))->Sumw2();

    registry.add("hNSubRatio21VsPt_Kt", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_Kt"))->Sumw2();
    registry.add("hNSubRatio21VsPt_CA", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_CA"))->Sumw2();
    registry.add("hNSubRatio21VsPt_SD", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_SD"))->Sumw2();

    registry.add("hDeltaR_Kt", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_Kt"))->Sumw2();
    registry.add("hDeltaR_CA", "#Delta R distribution with C/A subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_CA"))->Sumw2();
    registry.add("hDeltaR_SD", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_SD"))->Sumw2();

    registry.add("hDeltaRVsPt_Kt", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_Kt"))->Sumw2();
    registry.add("hDeltaRVsPt_CA", "#Delta R distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_CA"))->Sumw2();
    registry.add("hDeltaRVsPt_SD", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_SD"))->Sumw2();

    registry.get<TH1>(HIST("hErrorControl"))->GetXaxis()->SetBinLabel(1, "Single track jet");
    registry.get<TH1>(HIST("hErrorControl"))->GetXaxis()->SetBinLabel(2, "Single track groomed jet");
  }

  // function that converts a jet from the O2Physics jet table into a pseudojet; the fastjet cluster sequence of the jet needs to be given as input and will be modified by the function to save the clustering information
  template <typename JetTableElement>
  void jetToPseudoJet(JetTableElement const& jet, fastjet::PseudoJet& pseudoJet, std::vector<fastjet::PseudoJet> jetConstituents, fastjet::ClusterSequence& clusterSeqInput)
  {
    for (auto& jetConstituent : jet.template tracks_as<TrackTable>()) {
      FastJetUtilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }

    fastjet::JetDefinition jet_def(fastjet::JetAlgorithm::kt_algorithm, 2.5 * jet.r() / 100);
    fastjet::ClusterSequence clusterSeq(jetConstituents, jet_def);

    clusterSeqInput.transfer_from_sequence(clusterSeq);
    pseudoJet = sorted_by_pt(clusterSeqInput.inclusive_jets())[0]; // returns the jet found by the reclustering with the largest pT; should be the one corresponding to the jet table element
  }

  // function that returns the N-subjettiness ratio and the distance betewwen the two axes considered for tau2, in the form of a vector
  template <typename AxesTypeArg>
  std::vector<float> getNsubRatio21(fastjet::PseudoJet pseudoJet, float const& jetR, AxesTypeArg const& axesType, std::string const& axesTypeString)
  {
    std::vector<float> nSub_results;

    fastjet::contrib::Nsubjettiness nSub_1(1, axesType, fastjet::contrib::NormalizedMeasure(1.0, jetR));
    fastjet::contrib::Nsubjettiness nSub_2(2, axesType, fastjet::contrib::NormalizedMeasure(1.0, jetR));

    if (axesTypeString == "SD") {
      fastjet::contrib::SoftDrop softdropAlgo(SD_beta, SD_z_cut);
      pseudoJet = softdropAlgo(pseudoJet);
      if (pseudoJet.constituents().size() < 2) { // this analsysis requires at least 2 subjets in the jet which in turn requires at least two constituents
        nSub_results = {-998., -998.};           // error values
        return nSub_results;
      }
    }
    float tau1 = nSub_1.result(pseudoJet);
    float tau2 = nSub_2.result(pseudoJet);
    std::vector<fastjet::PseudoJet> axes_nSub_2 = nSub_2.currentAxes(); // gets the two axes used in the 2-subjettiness calculation

    float nSubRatio_21 = tau2 / tau1;
    float DeltaR = axes_nSub_2[0].delta_R(axes_nSub_2[1]); // distance between axes

    nSub_results = {nSubRatio_21, DeltaR};
    return nSub_results;
  }

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);

  template <typename T>
  void processOneJet(T const& jet, float weight = 1)
  {
    if (jet.tracksIds().size() < 2) { // this analsysis requires at least 2 subjets in the jet which in turn requires at least two constituents
      registry.fill(HIST("hErrorControl"), -999);
    } else {
      fastjet::PseudoJet pseudoJet;
      jetConstituents.clear();
      jetToPseudoJet(jet, pseudoJet, jetConstituents, clusterSeq_pseudoJet);
      std::vector<float> nSub_Kt_results = getNsubRatio21(pseudoJet, jet.r() / 100., fastjet::contrib::KT_Axes(), "Kt");
      std::vector<float> nSub_CA_results = getNsubRatio21(pseudoJet, jet.r() / 100., fastjet::contrib::CA_Axes(), "CA");
      std::vector<float> nSub_SD_results = getNsubRatio21(pseudoJet, jet.r() / 100., fastjet::contrib::CA_Axes(), "SD");

      registry.fill(HIST("hNSubRatio21_Kt"), nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21_CA"), nSub_CA_results[0], weight);
      registry.fill(HIST("hDeltaR_Kt"), nSub_Kt_results[1], weight);
      registry.fill(HIST("hDeltaR_CA"), nSub_CA_results[1], weight);

      registry.fill(HIST("hNSubRatio21VsPt_Kt"), jet.pt(), nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21VsPt_CA"), jet.pt(), nSub_CA_results[0], weight);
      registry.fill(HIST("hDeltaRVsPt_Kt"), jet.pt(), nSub_Kt_results[1], weight);
      registry.fill(HIST("hDeltaRVsPt_CA"), jet.pt(), nSub_CA_results[1], weight);

      if (nSub_SD_results[0] > -900) { // makes sure the grooming found at least one splitting that satisfied the Soft Drop condition before filling the result histograms for SD, like it is done for kT and CA
        registry.fill(HIST("hNSubRatio21_SD"), nSub_SD_results[0], weight);
        registry.fill(HIST("hDeltaR_SD"), nSub_SD_results[1], weight);

        registry.fill(HIST("hNSubRatio21VsPt_SD"), jet.pt(), nSub_SD_results[0], weight);
        registry.fill(HIST("hDeltaRVsPt_SD"), jet.pt(), nSub_SD_results[1], weight);
      } else {
        registry.fill(HIST("hErrorControl"), -998);
      }
    }
  }

  void processJets(typename CollisionTable::iterator const& collision, soa::Filtered<JetTable> const& jets, TrackTable const& tracks)
  {
    for (auto& jet : jets) {
      processOneJet(jet);
    }
  }
  PROCESS_SWITCH(nSubJettiness, processJets, "Process function, turned off by default", false);

  void processJetsWeighted(typename CollisionTable::iterator const& collision, soa::Filtered<soa::Join<JetTable, JetEventWeightTable>> const& jets, TrackTable const& tracks)
  {
    float weight;
    for (auto& jet : jets) {
      weight = jet.eventWeight();
      processOneJet(jet, weight);
    }
  }
  PROCESS_SWITCH(nSubJettiness, processJetsWeighted, "Process function with weighted events, turned off by default", false);

  void processDummy(aod::Tracks const& track)
  {
  }
  PROCESS_SWITCH(nSubJettiness, processDummy, "Dummy process function, turned on by default", true);
};
using NSubjettinessChargedJetMCParticleLevel = nSubJettiness<aod::McCollisions, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>, aod::ChargedMCParticleLevelJetEventWeights, aod::McParticles>;
using NSubjettinessChargedJetMCDetectorLevel = nSubJettiness<aod::Collisions, soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::ChargedMCDetectorLevelJetEventWeights, aod::Tracks>;
// using NSubjettinessChargedJetMCParticleLevel = nSubJettiness<aod::McCollisions, soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>, aod::McParticles>;
// using NSubjettinessChargedJetMCDetectorLevel = nSubJettiness<aod::Collisions,   soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>, aod::Tracks>;
// using NSubjettinessChargedJetDataLevel       = nSubJettiness<aod::Collisions,   soa::Join<aod::ChargedJets,                aod::ChargedJetConstituents>,                aod::Tracks>;

// using NSubjettinessD0ChargedJetMCParticleLevel = nSubJettiness<aod::McCollisions, soa::Join<aod::D0ChargedMCParticleLevelJets, aod::D0ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights>, aod::McParticles>;
// using NSubjettinessD0ChargedJetMCDetectorLevel = nSubJettiness<aod::Collisions,   soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>, aod::Tracks>;
// // using NSubjettinessD0ChargedJetDataLevel       = nSubJettiness<aod::Collisions,   soa::Join<aod::D0ChargedJets,                aod::D0ChargedJetConstituents>,                                                            aod::Tracks>;

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  std::vector<o2::framework::DataProcessorSpec> tasks;

  tasks.emplace_back(adaptAnalysisTask<NSubjettinessChargedJetMCParticleLevel>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-nsubjettiness-charged-mcp"}));

  tasks.emplace_back(adaptAnalysisTask<NSubjettinessChargedJetMCDetectorLevel>(cfgc,
                                                                               SetDefaultProcesses{},
                                                                               TaskName{"jet-nsubjettiness-charged-mcd"}));

  // tasks.emplace_back(adaptAnalysisTask<NSubjettinessChargedJetDataLevel>(cfgc,
  //                                                                SetDefaultProcesses{},
  //                                                                TaskName{"jet-nsubjettiness-charged-data"}));

  // tasks.emplace_back(adaptAnalysisTask<NSubjettinessD0ChargedJetMCParticleLevel>(cfgc,
  //                                                                      SetDefaultProcesses{},
  //                                                                      TaskName{"jet-nsubjettiness-d0-mcp"}));

  // tasks.emplace_back(adaptAnalysisTask<NSubjettinessD0ChargedJetMCDetectorLevel>(cfgc,
  //                                                                      SetDefaultProcesses{},
  //                                                                      TaskName{"jet-nsubjettiness-d0-mcd"}));

  // tasks.emplace_back(adaptAnalysisTask<NSubjettinessD0ChargedJetDataLevel>(cfgc,
  //                                                                SetDefaultProcesses{},
  //                                                                TaskName{"jet-nsubjettiness-d0-data"}));

  return WorkflowSpec{tasks};
}
