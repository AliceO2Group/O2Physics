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
/// \author Nima Zardoshti <nima.zardoshti@cern.ch>

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

namespace o2::aod
{
namespace nsubjettiness
{
DECLARE_SOA_COLUMN(JetPt, jetPt, float);
DECLARE_SOA_COLUMN(JetEta, jetEta, float);
DECLARE_SOA_COLUMN(JetPhi, jetPhi, float);
DECLARE_SOA_COLUMN(Nsub1Kt, nsub1kt, float);
DECLARE_SOA_COLUMN(Nsub2Kt, nsub2kt, float);
DECLARE_SOA_COLUMN(DeltaRKt, deltaRKt, float);
DECLARE_SOA_COLUMN(Nsub1CA, nsub1CA, float);
DECLARE_SOA_COLUMN(Nsub2CA, nsub2CA, float);
DECLARE_SOA_COLUMN(DeltaRCA, deltaRCA, float);
DECLARE_SOA_COLUMN(Nsub1CASD, nsub1CASD, float);
DECLARE_SOA_COLUMN(Nsub2CASD, nsub2CASD, float);
DECLARE_SOA_COLUMN(DeltaRCASD, deltaRCASD, float);
DECLARE_SOA_COLUMN(HfPt, hfPt, float);
DECLARE_SOA_COLUMN(HfEta, hfEta, float);
DECLARE_SOA_COLUMN(HfPhi, hfPhi, float);
DECLARE_SOA_COLUMN(HfY, hfY, float);
DECLARE_SOA_COLUMN(HfMass, hfMass, float);
DECLARE_SOA_DYNAMIC_COLUMN(MC, mc, []() -> int { return 0; });
} // namespace nsubjettiness
DECLARE_SOA_TABLE(JetTable, "AOD", "JETTABLE",
                  nsubjettiness::JetPt,
                  nsubjettiness::JetEta,
                  nsubjettiness::JetPhi,
                  nsubjettiness::Nsub1Kt,
                  nsubjettiness::Nsub2Kt,
                  nsubjettiness::DeltaRKt,
                  nsubjettiness::Nsub1CA,
                  nsubjettiness::Nsub2CA,
                  nsubjettiness::DeltaRCA,
                  nsubjettiness::Nsub1CASD,
                  nsubjettiness::Nsub2CASD,
                  nsubjettiness::DeltaRCASD);

DECLARE_SOA_TABLE(JetTableMC, "AOD", "JETTABLEMC",
                  nsubjettiness::JetPt,
                  nsubjettiness::JetEta,
                  nsubjettiness::JetPhi,
                  nsubjettiness::Nsub1Kt,
                  nsubjettiness::Nsub2Kt,
                  nsubjettiness::DeltaRKt,
                  nsubjettiness::Nsub1CA,
                  nsubjettiness::Nsub2CA,
                  nsubjettiness::DeltaRCA,
                  nsubjettiness::Nsub1CASD,
                  nsubjettiness::Nsub2CASD,
                  nsubjettiness::DeltaRCASD,
                  nsubjettiness::MC<>);

DECLARE_SOA_TABLE(CandTable, "AOD", "CANDTABLE",
                  nsubjettiness::HfPt,
                  nsubjettiness::HfEta,
                  nsubjettiness::HfPhi,
                  nsubjettiness::HfY,
                  nsubjettiness::HfMass);

DECLARE_SOA_TABLE(CandTableMC, "AOD", "CANDTABLEMC",
                  nsubjettiness::HfPt,
                  nsubjettiness::HfEta,
                  nsubjettiness::HfPhi,
                  nsubjettiness::HfY,
                  nsubjettiness::MC<>);
} // namespace o2::aod

struct NSubjettinessTask {

  Produces<aod::JetTable> jetTable;
  Produces<aod::JetTableMC> jetTableMC;
  Produces<aod::CandTable> candidateTable;
  Produces<aod::CandTableMC> candidateTableMC;
  HistogramRegistry registry;

  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};

  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> SD_beta{"SD_beta", 0.0, "SoftDrop beta"};
  Configurable<float> SD_z_cut{"SD_z_cut", 0.10, "SoftDrop z cut"};

  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 100.0f}, ""};
  ConfigurableAxis DeltaRBinning{"Delta-R-binning", {50, 0.0, 0.5}, ""};
  ConfigurableAxis NSubRatioBinning{"NSub-Ratio-binning", {50, 0.0f, 1.2f}, ""};

  JetFinder jetReclusterer;
  std::vector<float> nSub_Kt_results;
  std::vector<float> nSub_CA_results;
  std::vector<float> nSub_CASD_results;

  std::vector<fastjet::PseudoJet> dummyVector;
  fastjet::JetDefinition dummyJetDef{fastjet::antikt_algorithm, 0.1, fastjet::E_scheme, fastjet::Best};
  fastjet::GhostedAreaSpec ghostareaspec{0.0, 1, 0.05};
  fastjet::AreaDefinition dummyAreaDef{fastjet::active_area, ghostareaspec};
  fastjet::ClusterSequenceArea clusterSeq{dummyVector, dummyJetDef, dummyAreaDef};

  void init(InitContext const&)
  {
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/c)"};
    AxisSpec DeltaRAxis = {DeltaRBinning, "#Delta R"};
    AxisSpec NSubRatioAxis = {NSubRatioBinning, "#tau_{2}/#tau_{1}"};

    registry.add("hNSubRatio21_Kt", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_Kt"))->Sumw2();
    registry.add("hNSubRatio21_CA", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_CA"))->Sumw2();
    registry.add("hNSubRatio21_CASD", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_CASD"))->Sumw2();

    registry.add("hNSubRatio21_Kt_Part", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_Kt_Part"))->Sumw2();
    registry.add("hNSubRatio21_CA_Part", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_CA_Part"))->Sumw2();
    registry.add("hNSubRatio21_CASD_Part", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH1F, {NSubRatioAxis}});
    registry.get<TH1>(HIST("hNSubRatio21_CASD_Part"))->Sumw2();

    registry.add("hNSubRatio21VsPt_Kt", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_Kt"))->Sumw2();
    registry.add("hNSubRatio21VsPt_CA", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_CA"))->Sumw2();
    registry.add("hNSubRatio21VsPt_CASD", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_CASD"))->Sumw2();

    registry.add("hNSubRatio21VsPt_Kt_Part", "#tau_{2}/#tau_{1} distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_Kt_Part"))->Sumw2();
    registry.add("hNSubRatio21VsPt_CA_Part", "#tau_{2}/#tau_{1} distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_CA_Part"))->Sumw2();
    registry.add("hNSubRatio21VsPt_CASD_Part", "#tau_{2}/#tau_{1} distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, NSubRatioAxis}});
    registry.get<TH2>(HIST("hNSubRatio21VsPt_CASD_Part"))->Sumw2();

    registry.add("hDeltaR_Kt", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_Kt"))->Sumw2();
    registry.add("hDeltaR_CA", "#Delta R distribution with C/A subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_CA"))->Sumw2();
    registry.add("hDeltaR_CASD", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_CASD"))->Sumw2();

    registry.add("hDeltaR_Kt_Part", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_Kt_Part"))->Sumw2();
    registry.add("hDeltaR_CA_Part", "#Delta R distribution with C/A subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_CA_Part"))->Sumw2();
    registry.add("hDeltaR_CASD_Part", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH1F, {DeltaRAxis}});
    registry.get<TH1>(HIST("hDeltaR_CASD_Part"))->Sumw2();

    registry.add("hDeltaRVsPt_Kt", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_Kt"))->Sumw2();
    registry.add("hDeltaRVsPt_CA", "#Delta R distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_CA"))->Sumw2();
    registry.add("hDeltaRVsPt_CASD", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_CASD"))->Sumw2();

    registry.add("hDeltaRVsPt_Kt_Part", "#Delta R distribution with k_{T} subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_Kt_Part"))->Sumw2();
    registry.add("hDeltaRVsPt_CA_Part", "#Delta R distribution with C/A subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_CA_Part"))->Sumw2();
    registry.add("hDeltaRVsPt_CASD_Part", "#Delta R distribution with Softdrop subjet axes", {HistType::kTH2F, {ptAxis, DeltaRAxis}});
    registry.get<TH2>(HIST("hDeltaRVsPt_CASD_Part"))->Sumw2();

    jetReclusterer.isReclustering = true;
    jetReclusterer.algorithm = fastjet::JetAlgorithm::antikt_algorithm;
    jetReclusterer.jetR = jetR;
  }

  // function that converts a jet from the O2Physics jet table into a pseudojet; the fastjet cluster sequence of the jet needs to be given as input and will be modified by the function to save the clustering information
  template <bool isHF, typename T, typename U, typename V>
  fastjet::PseudoJet jetToPseudoJet(T const& jet, U const& tracks, V const candidates)
  {
    std::vector<fastjet::PseudoJet> jetConstituents;
    for (auto& jetConstituent : jet.template tracks_as<U>()) {
      fastjetutilities::fillTracks(jetConstituent, jetConstituents, jetConstituent.globalIndex());
    }
    if constexpr (isHF) {
      for (auto& jetHFConstituent : jet.template hfcandidates_as<V>()) {
        fastjetutilities::fillTracks(jetHFConstituent, jetConstituents, jetHFConstituent.globalIndex());
      }
    }
    std::vector<fastjet::PseudoJet> jetReclustered;
    clusterSeq = jetReclusterer.findJets(jetConstituents, jetReclustered);
    jetReclustered = sorted_by_pt(jetReclustered);
    return jetReclustered[0];
  }

  // function that returns the N-subjettiness ratio and the distance betewwen the two axes considered for tau2, in the form of a vector
  template <typename AxesTypeArg>
  std::vector<float> getNSubjettiness(fastjet::PseudoJet jet, float const& jetR, AxesTypeArg const& axesType, bool softDrop)
  {
    if (jet.constituents().size() < 2) { // this analsysis requires at least 2 subjets in the jet which in turn requires at least two constituents
      return {-1.0, -1.0, -1.0};         // error values
    }
    fastjet::contrib::Nsubjettiness nSub_1(1, axesType, fastjet::contrib::NormalizedMeasure(1.0, jetR));
    fastjet::contrib::Nsubjettiness nSub_2(2, axesType, fastjet::contrib::NormalizedMeasure(1.0, jetR));

    if (softDrop) {
      fastjet::contrib::SoftDrop softdropAlgo(SD_beta, SD_z_cut);
      jet = softdropAlgo(jet);
      if (jet.constituents().size() < 2) { // this analsysis requires at least 2 subjets in the jet which in turn requires at least two constituents
        return {-2.0, -2.0, -2.0};         // error values
      }
    }
    float tau1 = nSub_1.result(jet);
    float tau2 = nSub_2.result(jet);
    std::vector<fastjet::PseudoJet> axes_nSub_2 = nSub_2.currentAxes(); // gets the two axes used in the 2-subjettiness calculation

    float DeltaR = axes_nSub_2[0].delta_R(axes_nSub_2[1]); // distance between axes

    return {tau1, tau2, DeltaR};
  }

  Filter jetCuts = aod::jet::r == nround(jetR.node() * 100.0f);
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);

  template <bool isMCP, bool isHF, typename T, typename U, typename V>
  void processJet(T const& jet, U const& tracks, V const& candidates, float weight = 1.0)
  {
    fastjet::PseudoJet pseudoJet(jetToPseudoJet<isHF>(jet, tracks, candidates));
    nSub_Kt_results = getNSubjettiness(pseudoJet, jet.r() / 100., fastjet::contrib::KT_Axes(), false);
    nSub_CA_results = getNSubjettiness(pseudoJet, jet.r() / 100., fastjet::contrib::CA_Axes(), false);
    nSub_CASD_results = getNSubjettiness(pseudoJet, jet.r() / 100., fastjet::contrib::CA_Axes(), true);

    if constexpr (isMCP) {

      registry.fill(HIST("hNSubRatio21_Kt_Part"), nSub_Kt_results[1] / nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21_CA_Part"), nSub_CA_results[1] / nSub_CA_results[0], weight);
      registry.fill(HIST("hNSubRatio21_CASD_Part"), nSub_CASD_results[1] / nSub_CASD_results[0], weight);
      registry.fill(HIST("hDeltaR_Kt_Part"), nSub_Kt_results[2], weight);
      registry.fill(HIST("hDeltaR_CA_Part"), nSub_CA_results[2], weight);
      registry.fill(HIST("hDeltaR_CASD_Part"), nSub_CASD_results[2], weight);

      registry.fill(HIST("hNSubRatio21VsPt_Kt_Part"), jet.pt(), nSub_Kt_results[1] / nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21VsPt_CA_Part"), jet.pt(), nSub_CA_results[1] / nSub_CA_results[0], weight);
      registry.fill(HIST("hNSubRatio21VsPt_CASD_Part"), jet.pt(), nSub_CASD_results[1] / nSub_CASD_results[0], weight);
      registry.fill(HIST("hDeltaRVsPt_Kt_Part"), jet.pt(), nSub_Kt_results[2], weight);
      registry.fill(HIST("hDeltaRVsPt_CA_Part"), jet.pt(), nSub_CA_results[2], weight);
      registry.fill(HIST("hDeltaRVsPt_CASD_Part"), jet.pt(), nSub_CASD_results[2], weight);

    } else {
      registry.fill(HIST("hNSubRatio21_Kt"), nSub_Kt_results[1] / nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21_CA"), nSub_CA_results[1] / nSub_CA_results[0], weight);
      registry.fill(HIST("hNSubRatio21_CASD"), nSub_CASD_results[1] / nSub_CASD_results[0], weight);
      registry.fill(HIST("hDeltaR_Kt"), nSub_Kt_results[2], weight);
      registry.fill(HIST("hDeltaR_CA"), nSub_CA_results[2], weight);
      registry.fill(HIST("hDeltaR_CASD"), nSub_CASD_results[2], weight);

      registry.fill(HIST("hNSubRatio21VsPt_Kt"), jet.pt(), nSub_Kt_results[1] / nSub_Kt_results[0], weight);
      registry.fill(HIST("hNSubRatio21VsPt_CA"), jet.pt(), nSub_CA_results[1] / nSub_CA_results[0], weight);
      registry.fill(HIST("hNSubRatio21VsPt_CASD"), jet.pt(), nSub_CASD_results[1] / nSub_CASD_results[0], weight);
      registry.fill(HIST("hDeltaRVsPt_Kt"), jet.pt(), nSub_Kt_results[2], weight);
      registry.fill(HIST("hDeltaRVsPt_CA"), jet.pt(), nSub_CA_results[2], weight);
      registry.fill(HIST("hDeltaRVsPt_CASD"), jet.pt(), nSub_CASD_results[2], weight);
    }
  }

  template <typename T, typename U>
  void fillJetTable(T const& jet, U& table)
  {
    table(jet.pt(), jet.eta(), jet.phi(), nSub_Kt_results[0], nSub_Kt_results[1], nSub_Kt_results[2], nSub_CA_results[0], nSub_CA_results[1], nSub_CA_results[2], nSub_CASD_results[0], nSub_CASD_results[1], nSub_CASD_results[2]);
  }

  template <bool isMC, typename T, typename U>
  void fillCandidateTable(T const& candidate, U& table)
  {
    if constexpr (isMC) {
      table(candidate.pt(), candidate.eta(), candidate.phi(), candidate.y());
    } else {
      table(candidate.pt(), candidate.eta(), candidate.phi(), candidate.y(), candidate.m());
    }
  }

  void processJetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& jets, JetTracks const& tracks)
  {
    for (auto& jet : jets) {
      processJet<false, false>(jet, tracks, tracks);
      fillJetTable(jet, jetTable);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsData, "Process function for inclusive jets in data", true);

  void processJetsDataEWS(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedEventWiseSubtractedJets, aod::ChargedEventWiseSubtractedJetConstituents>> const& jets, JetTracksSub const& tracks)
  {
    for (auto& jet : jets) {
      processJet<false, false>(jet, tracks, tracks);
      fillJetTable(jet, jetTable);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsDataEWS, "Process function for inclusive jets with eventwise subtraction in data", true);

  void processJetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents>> const& jets, JetTracks const& tracks)
  {
    for (auto& jet : jets) {
      processJet<false, false>(jet, tracks, tracks);
      fillJetTable(jet, jetTable);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsMCD, "Process function for inclusive jets in mcd", false);

  void processJetsMCDWeighted(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetEventWeights>> const& jets, JetTracks const& tracks)
  {
    for (auto& jet : jets) {
      processJet<false, false>(jet, tracks, tracks, jet.eventWeight());
      fillJetTable(jet, jetTable);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsMCDWeighted, "Process function for inclusive jets in weighted mcd", false);

  void processJetsMCP(JetMcCollision const& mcCollision, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents>> const& jets, JetParticles const& particles)
  {
    for (auto& jet : jets) {
      processJet<true, false>(jet, particles, particles);
      fillJetTable(jet, jetTableMC);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsMCP, "Process function for inclusive jets in mcp", false);

  void processJetsMCPWeighted(JetMcCollision const& mcCollision, soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetEventWeights>> const& jets, JetParticles const& particles)
  {
    for (auto& jet : jets) {
      processJet<true, false>(jet, particles, particles, jet.eventWeight());
      fillJetTable(jet, jetTableMC);
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processJetsMCPWeighted, "Process function for inclusive jets in weighted mcp", false);

  void processD0JetsData(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::D0ChargedJets, aod::D0ChargedJetConstituents>> const& jets, JetTracks const& tracks, CandidatesD0Data const& candidates)
  {
    for (auto& jet : jets) {
      processJet<false, true>(jet, tracks, candidates);
      fillJetTable(jet, jetTable);
      for (auto& candidate : jet.hfcandidates_as<CandidatesD0Data>()) {
        fillCandidateTable<false>(candidate, candidateTable);
        break;
      }
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processD0JetsData, "Process function for D0 jets in data", false);

  void processD0JetsMCD(soa::Filtered<JetCollisions>::iterator const& collision, soa::Filtered<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>> const& jets, JetTracks const& tracks, CandidatesD0MCD const& candidates)
  {
    for (auto& jet : jets) {
      processJet<false, true>(jet, tracks, candidates);
      fillJetTable(jet, jetTable);
      for (auto& candidate : jet.hfcandidates_as<CandidatesD0MCD>()) {
        fillCandidateTable<false>(candidate, candidateTable);
        break;
      }
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processD0JetsMCD, "Process function for D0 jets in mcd", false);

  void processD0JetsMCP(JetMcCollision const& mcCollision, soa::Filtered<soa::Join<aod::D0ChargedMCDetectorLevelJets, aod::D0ChargedMCDetectorLevelJetConstituents>> const& jets, JetParticles const& particles, CandidatesD0MCP const& hfParticles)
  {
    for (auto& jet : jets) {
      processJet<true, true>(jet, particles, hfParticles);
      fillJetTable(jet, jetTable);
      for (auto& hfParticle : jet.hfcandidates_as<CandidatesD0MCP>()) {
        fillCandidateTable<true>(hfParticle, candidateTableMC);
        break;
      }
    }
  }
  PROCESS_SWITCH(NSubjettinessTask, processD0JetsMCP, "Process function for D0 jets in mcp", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<NSubjettinessTask>(
    cfgc, TaskName{"jet-nsubjettiness"})};
}
