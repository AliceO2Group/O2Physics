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

/// \file   bjetTaggingGnn.cxx
/// \brief  b-jet tagging using GNN
///
/// \author Changhwan Choi <changhwan.choi@cern.ch>, Pusan National University

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Framework/ASoA.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include <Framework/Configurable.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BjetTaggingGnn {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> eventSelectionsSel{"eventSelectionsSel", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};
  Configurable<float> trackPtMinGnn{"trackPtMinGnn", 0.5, "minimum track pT for GNN inputs"};

  Configurable<float> maxIPxy{"maxIPxy", 10, "maximum track DCA in xy plane"};
  Configurable<float> maxIPz{"maxIPz", 10, "maximum track DCA in z direction"};

  Configurable<float> trackNppCrit{"trackNppCrit", 0.95, "track not physical primary ratio"};

  // sv level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<double> dbMin{"dbMin", -10., "minimum GNN Db"};
  Configurable<double> dbMax{"dbMax", 20., "maximum GNN Db"};
  Configurable<int> dbNbins{"dbNbins", 3000, "number of bins in axisDbFine"};

  Configurable<bool> doDataDriven{"doDataDriven", false, "Flag whether to use fill THnSpase for data driven methods"};
  Configurable<bool> callSumw2{"callSumw2", false, "Flag whether to call THnSparse::Sumw2() for error calculation"};

  Configurable<int> trainingDatasetRatioParam{"trainingDatasetRatioParam", 0, "Parameter for splitting training/evaluation datasets by collisionId"};

  std::vector<int> eventSelectionBits;
  std::vector<int> eventSelectionBitsSel;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    eventSelectionBitsSel = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelectionsSel));

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
    // registry.add("h_vertexZ_truth", "Vertex Z truth;#it{Z} (cm)", {HistType::kTH1F, {{100, -20.0, 20.0}}}, callSumw2);
    registry.add("h_event_counter", ";Event counter", {HistType::kTH1F, {{6, 0.0, 6.0}}}, callSumw2);
    auto hEventCounter = registry.get<TH1>(HIST("h_event_counter"));
    hEventCounter->GetXaxis()->SetBinLabel(1, "Coll+TVX");
    hEventCounter->GetXaxis()->SetBinLabel(2, "Coll+TVX+Sel8");
    hEventCounter->GetXaxis()->SetBinLabel(3, "Coll+TVX+Sel8+");
    hEventCounter->GetXaxis()->SetBinLabel(4, "McColl(INEL)");
    hEventCounter->GetXaxis()->SetBinLabel(5, "McColl(-> Coll+TVX+Sel8)");
    hEventCounter->GetXaxis()->SetBinLabel(6, "McColl(-> Coll+TVX+Sel8+)");
    registry.add("hBCCounter", "", {HistType::kTH1F, {{10, 0.0, 10.}}}, callSumw2);
    auto hBCCounter = registry.get<TH1>(HIST("hBCCounter"));
    hBCCounter->GetXaxis()->SetBinLabel(1, "AllBC");
    hBCCounter->GetXaxis()->SetBinLabel(2, "BC+TVX");
    hBCCounter->GetXaxis()->SetBinLabel(3, "BC+TVX+NoTFB");
    hBCCounter->GetXaxis()->SetBinLabel(4, "BC+TVX+NoTFB+NoITSROFB");
    hBCCounter->GetXaxis()->SetBinLabel(5, "CollinBC");
    hBCCounter->GetXaxis()->SetBinLabel(6, "CollinBC+Sel8");
    hBCCounter->GetXaxis()->SetBinLabel(7, "CollinBC+Sel8+VtxZ");
    hBCCounter->GetXaxis()->SetBinLabel(8, "CollinBC+Sel8Full");
    hBCCounter->GetXaxis()->SetBinLabel(9, "CollinBC+Sel8Full+GoodZvtx");
    hBCCounter->GetXaxis()->SetBinLabel(10, "CollinBC+Sel8Full+VtxZ+GoodZvtx");

    const AxisSpec axisTrackpT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTrackpTFine{1000, 0., 10., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisJetpT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDb{200, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisDbFine{dbNbins, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisJetMass{200, 0., 50., "#it{m}_{jet} (GeV/#it{c}^{2})"};
    const AxisSpec axisJetProb{200, 0., 40., "-ln(JP)"};
    const AxisSpec axisNTracks{42, 0, 42, "#it{n}_{tracks}"};

    registry.add("h_jetpT", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    registry.add("h_Db", "", {HistType::kTH1F, {axisDbFine}});
    registry.add("h2_jetpT_Db", "", {HistType::kTH2F, {axisJetpT, axisDb}});

    if (doprocessDataTracks || doprocessMCDTracks) {
      registry.add("h_trackpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_tracketa", "", {HistType::kTH1F, {{100, trackEtaMin, trackEtaMax, "#it{#eta}"}}}, callSumw2);
      registry.add("h_trackphi", "", {HistType::kTH1F, {{100, 0.0, 2.0 * M_PI, "#it{#phi}"}}}, callSumw2);
    }

    if (doprocessMCDTracks) {
      registry.add("h2_trackpT_partpT", "", {HistType::kTH2F, {axisTrackpT, axisTrackpT}}, callSumw2);
      registry.add("h_partpT_matched_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
      registry.add("h_partpT", "", {HistType::kTH1F, {axisTrackpT}}, callSumw2);
      registry.add("h_partpT_fine", "", {HistType::kTH1F, {axisTrackpTFine}}, callSumw2);
    }

    if (doprocessDataJetsSel || doprocessMCDJetsSel) {
      registry.add("h_jetpT_sel", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_tvx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    }

    if (doprocessMCDJets) {
      registry.add("h_jetpT_b", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_c", "c-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_lf", "lf-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_Db_b", "b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_c", "c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_lf", "lf-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h2_jetpT_Db_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_Response_DetjetpT_PartjetpT", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}, callSumw2});
      registry.add("h2_jetpT_Db_lf_none", "lf-jet (none)", {HistType::kTH2F, {axisJetpT, axisDb}}, callSumw2);
      registry.add("h2_jetpT_Db_lf_matched", "lf-jet (matched)", {HistType::kTH2F, {axisJetpT, axisDb}}, callSumw2);
      registry.add("h2_jetpT_Db_npp", "NotPhysPrim", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h_Db_npp", "NotPhysPrim", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_jetpT_matched", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_matched", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_matched", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_matched", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("hSparse_overflow_b", "", {HistType::kTHnSparseF, {axisJetpT, axisJetpT, axisNTracks}});
    }

    if (doprocessMCDJetsSel) {
      registry.add("h_jetpT_b_sel", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_sel", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_sel", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h_jetpT_b_tvx", "b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_tvx", "", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
      registry.add("h2_Response_DetjetpT_PartjetpT_b_tvx", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}}, callSumw2);
    }

    if (doprocessMCPJets) {
      registry.add("h_jetpT_particle", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_c", "particle c-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_lf", "particle lf-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_sel", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_sel", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_tvx", "", {HistType::kTH1F, {axisJetpT}}, callSumw2);
      registry.add("h_jetpT_particle_b_tvx", "particle b-jet", {HistType::kTH1F, {axisJetpT}}, callSumw2);
    }

    if (doDataDriven) {
      registry.add("hSparse_Incljets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
      if (doprocessMCDJets) {
        registry.add("hSparse_bjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        registry.add("hSparse_cjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets_none", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets_matched", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisNTracks}}, callSumw2);
      }
    }
  }

  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackFilter = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt < jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using AnalysisCollisions = soa::Join<aod::JetCollisions, aod::JCollisionPIs>;
  using FilteredCollisions = soa::Filtered<AnalysisCollisions>;
  using DataJets = soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetTags>;
  using FilteredDataJets = soa::Filtered<DataJets>;
  using AnalysisTracks = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracks = soa::Filtered<AnalysisTracks>;

  using MCDJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags>;
  using FilteredMCDJets = soa::Filtered<MCDJets>;
  using AnalysisTracksMCD = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  using FilteredTracksMCD = soa::Filtered<AnalysisTracksMCD>;

  using AnalysisCollisionsMCD = soa::Join<aod::JetCollisionsMCD, aod::JCollisionPIs>;
  using FilteredCollisionsMCD = soa::Filtered<AnalysisCollisionsMCD>;

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionsMCP = soa::Filtered<aod::JetMcCollisions>;
  using MCPJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef>;
  using FilteredMCPJets = soa::Filtered<MCPJets>;

  Service<o2::framework::O2DatabasePDG> pdg;

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks>
  int analyzeJetTrackInfo(AnyCollision const& /*collision*/, AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/ /*, int8_t jetFlavor = 0, double weight = 1.0*/)
  {
    int nTracks = 0;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin || !jettaggingutilities::trackAcceptanceWithDca(constituent, maxIPxy, maxIPz)) {
        continue;
      }

      // ...

      ++nTracks;
    }
    return nTracks;
  }

  void processDummy(FilteredCollisions::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollisions::iterator const& collision, FilteredDataJets const& alljets, FilteredTracks const& allTracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 2.5); // Coll+TVX+Sel8+...

    registry.fill(HIST("h_vertexZ"), collision.posZ());

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      int nTracks = analyzeJetTrackInfo(collision, analysisJet, allTracks);

      registry.fill(HIST("h_jetpT"), analysisJet.pt());
      registry.fill(HIST("h_Db"), analysisJet.scoreML());
      registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML());

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJets, "jet information in Data", false);

  void processDataJetsSel(AnalysisCollisions::iterator const& collision, FilteredDataJets const& alljets)
  {
    registry.fill(HIST("h_event_counter"), 0.5); // Coll+TVX

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      registry.fill(HIST("h_jetpT_tvx"), analysisJet.pt());
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSel)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 1.5); // Coll+TVX+Sel8

    for (const auto& analysisJet : alljets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      registry.fill(HIST("h_jetpT_sel"), analysisJet.pt());
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJetsSel, "jet information in Data (sel8)", false);

  void processDataTracks(FilteredCollisions::iterator const& collision, AnalysisTracks const& tracks)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    for (const auto& track : tracks) {
      if (track.eta() <= trackEtaMin || track.eta() >= trackEtaMax) {
        continue;
      }
      registry.fill(HIST("h_trackpT"), track.pt());
      registry.fill(HIST("h_tracketa"), track.eta());
      registry.fill(HIST("h_trackphi"), track.phi());
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataTracks, "track information in Data", false);

  void processMCDJets(FilteredCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredTracksMCD const& /*allTracks*/, FilteredMCPJets const& /*MCPjets*/, aod::JetParticles const& /*MCParticles*/, FilteredCollisionsMCP const& /*mcCollisions*/)
  {
    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    // Uses only collisionId % trainingDatasetRaioParam != 0 for evaluation dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam == 0) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 2.5, weightEvt);

    registry.fill(HIST("h_vertexZ"), collision.posZ(), weightEvt);

    bool matchedMcColl = collision.has_mcCollision() && std::fabs(collision.template mcCollision_as<FilteredCollisionsMCP>().posZ()) < vertexZCut;

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      // int nTracks = analyzeJetTrackInfo(collision, analysisJet, allTracks /*, jetFlavor, weight*/);
      int nTracks = 0;

      int nNppTracks = 0;
      for (const auto& constituent : analysisJet.template tracks_as<FilteredTracksMCD>()) {
        if (constituent.pt() < trackPtMinGnn) {
          continue;
        }
        if (!constituent.has_mcParticle() || !constituent.template mcParticle_as<aod::JetParticles>().isPhysicalPrimary()) {
          ++nNppTracks;
        }
        ++nTracks;
      }

      registry.fill(HIST("h_jetpT"), analysisJet.pt(), weightEvt);
      registry.fill(HIST("h_Db"), analysisJet.scoreML(), weightEvt);
      registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_b"), analysisJet.pt(), weightEvt);
        registry.fill(HIST("h_Db_b"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_b"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_c"), analysisJet.pt(), weightEvt);
        registry.fill(HIST("h_Db_c"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_c"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
      } else {
        registry.fill(HIST("h_jetpT_lf"), analysisJet.pt(), weightEvt);
        registry.fill(HIST("h_Db_lf"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_lf"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::none) {
          registry.fill(HIST("h2_jetpT_Db_lf_none"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h2_jetpT_Db_lf_matched"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        }
      }

      if (static_cast<float>(nNppTracks) / nTracks > trackNppCrit) {
        registry.fill(HIST("h_Db_npp"), analysisJet.scoreML(), weightEvt);
        registry.fill(HIST("h2_jetpT_Db_npp"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h_Db_npp_b"), analysisJet.scoreML(), weightEvt);
          registry.fill(HIST("h2_jetpT_Db_npp_b"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h_Db_npp_c"), analysisJet.scoreML(), weightEvt);
          registry.fill(HIST("h2_jetpT_Db_npp_c"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        } else {
          registry.fill(HIST("h_Db_npp_lf"), analysisJet.scoreML(), weightEvt);
          registry.fill(HIST("h2_jetpT_Db_npp_lf"), analysisJet.pt(), analysisJet.scoreML(), weightEvt);
        }
      }

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("hSparse_bjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("hSparse_cjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
        } else {
          registry.fill(HIST("hSparse_lfjets"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          if (jetFlavor == JetTaggingSpecies::none) {
            registry.fill(HIST("hSparse_lfjets_none"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          } else {
            registry.fill(HIST("hSparse_lfjets_matched"), analysisJet.pt(), analysisJet.scoreML(), nTracks, weightEvt);
          }
        }
      }

      if (!matchedMcColl) {
        continue;
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        registry.fill(HIST("h_jetpT_matched"), analysisJet.pt(), weightEvt);
        registry.fill(HIST("h_jetpT_particle_matched"), mcpjet.pt(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b"), analysisJet.pt(), mcpjet.pt(), weightEvt);
          registry.fill(HIST("h_jetpT_b_matched"), analysisJet.pt(), weightEvt);
          registry.fill(HIST("h_jetpT_particle_b_matched"), mcpjet.pt(), weightEvt);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_c"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lf"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJets, "jet information in MC", false);

  void processMCDJetsSel(AnalysisCollisionsMCD::iterator const& collision, FilteredMCDJets const& MCDjets, FilteredMCPJets const& /*MCPjets*/)
  {
    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    registry.fill(HIST("h_event_counter"), 0.5, weightEvt);

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      registry.fill(HIST("h_jetpT_tvx"), analysisJet.pt(), weightEvt);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_b_tvx"), analysisJet.pt(), weightEvt);
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_tvx"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_tvx"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        }
      }
    }

    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBitsSel)) {
      return;
    }
    registry.fill(HIST("h_event_counter"), 1.5, weightEvt);

    for (const auto& analysisJet : MCDjets) {

      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (analysisJet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      registry.fill(HIST("h_jetpT_sel"), analysisJet.pt(), weightEvt);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_b_sel"), analysisJet.pt(), weightEvt);
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<FilteredMCPJets>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_sel"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b_sel"), analysisJet.pt(), mcpjet.pt(), weightEvt);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDJetsSel, "jet information in MC (sel8)", false);

  PresliceUnsorted<AnalysisCollisionsMCD> collisionsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;
  Preslice<FilteredMCPJets> mcpjetsPerMCPCollision = aod::jmccollisionlb::mcCollisionId;

  void processMCPJets(aod::McCollisions::iterator const& mcCollision, FilteredMCPJets const& mcpjets, AnalysisCollisionsMCD const& collisions, aod::JetParticles const& /*MCParticles*/)
  {
    float weightEvt = useEventWeight ? mcCollision.weight() : 1.f;
    registry.fill(HIST("h_event_counter"), 3.5, weightEvt); // McColl(INEL)
    auto collisionspermccollision = collisions.sliceBy(collisionsPerMCPCollision, mcCollision.globalIndex());
    if (collisionspermccollision.size() >= 1) {
      if (jetderiveddatautilities::selectCollision(collisionspermccollision.begin(), eventSelectionBitsSel)) {
        registry.fill(HIST("h_event_counter"), 4.5, weightEvt); // McColl(-> Coll+TVX+Sel8)
      }
      if (jetderiveddatautilities::selectCollision(collisionspermccollision.begin(), eventSelectionBits) && std::fabs(collisionspermccollision.begin().posZ()) < vertexZCut) {
        registry.fill(HIST("h_event_counter"), 5.5, weightEvt); // McColl(-> Coll+TVX+Sel8+...)
      }
    }

    auto mcpjetspermcpcollision = mcpjets.sliceBy(mcpjetsPerMCPCollision, mcCollision.globalIndex());
    for (const auto& mcpjet : mcpjetspermcpcollision) {
      bool jetIncluded = false;
      for (const auto& jetR : jetRadiiValues) {
        if (mcpjet.r() == static_cast<int>(jetR * 100)) {
          jetIncluded = true;
          break;
        }
      }

      if (!jetIncluded) {
        continue;
      }

      float pTHat = 10. / (std::pow(weightEvt, 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int8_t jetFlavor = mcpjet.origin();

      registry.fill(HIST("h_jetpT_particle_tvx"), mcpjet.pt(), weightEvt);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_b_tvx"), mcpjet.pt(), weightEvt);
      }

      if (collisionspermccollision.size() >= 1) {
        if (jetderiveddatautilities::selectCollision(collisionspermccollision.begin(), eventSelectionBitsSel)) {
          registry.fill(HIST("h_jetpT_particle_sel"), mcpjet.pt(), weightEvt);

          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("h_jetpT_particle_b_sel"), mcpjet.pt(), weightEvt);
          }
        }

        if (jetderiveddatautilities::selectCollision(collisionspermccollision.begin(), eventSelectionBits) && std::fabs(collisionspermccollision.begin().posZ()) < vertexZCut && std::fabs(mcCollision.posZ()) < vertexZCut) {
          registry.fill(HIST("h_jetpT_particle"), mcpjet.pt(), weightEvt);

          if (jetFlavor == JetTaggingSpecies::beauty) {
            registry.fill(HIST("h_jetpT_particle_b"), mcpjet.pt(), weightEvt);
          } else if (jetFlavor == JetTaggingSpecies::charm) {
            registry.fill(HIST("h_jetpT_particle_c"), mcpjet.pt(), weightEvt);
          } else {
            registry.fill(HIST("h_jetpT_particle_lf"), mcpjet.pt(), weightEvt);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCPJets, "mc collision information", false);

  Preslice<aod::JetParticles> mcparticlesPerMCPCollision = aod::jmcparticle::mcCollisionId;

  void processMCDTracks(FilteredCollisionsMCD::iterator const& collision, AnalysisTracksMCD const& tracks, FilteredCollisionsMCP const& /*mcCollisions*/, aod::JetParticles const& allParticles)
  {
    float weightEvt = useEventWeight ? collision.weight() : 1.f;
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    // Uses only collisionId % trainingDatasetRaioParam != 0 for evaluation dataset
    if (trainingDatasetRatioParam && collision.collisionId() % trainingDatasetRatioParam == 0) {
      return;
    }

    bool matchedMcColl = collision.has_mcCollision() && std::fabs(collision.template mcCollision_as<FilteredCollisionsMCP>().posZ()) < vertexZCut;

    for (const auto& track : tracks) {
      if (track.eta() <= trackEtaMin || track.eta() >= trackEtaMax) {
        continue;
      }
      registry.fill(HIST("h_trackpT"), track.pt(), weightEvt);
      registry.fill(HIST("h_tracketa"), track.eta(), weightEvt);
      registry.fill(HIST("h_trackphi"), track.phi(), weightEvt);

      if (!matchedMcColl || !track.has_mcParticle()) {
        continue;
      }
      auto particle = track.template mcParticle_as<aod::JetParticles>();
      if (particle.isPhysicalPrimary() && particle.eta() > trackEtaMin && particle.eta() < trackEtaMax) {
        registry.fill(HIST("h2_trackpT_partpT"), track.pt(), particle.pt(), weightEvt);
        registry.fill(HIST("h_partpT_matched_fine"), particle.pt(), weightEvt);
      }
    }

    if (!matchedMcColl) {
      return;
    }

    auto const particles = allParticles.sliceBy(mcparticlesPerMCPCollision, collision.mcCollisionId());

    for (const auto& particle : particles) {
      auto pdgParticle = pdg->GetParticle(particle.pdgCode());
      if (!pdgParticle || pdgParticle->Charge() == 0.0) {
        continue;
      }
      if (particle.isPhysicalPrimary() && particle.eta() > trackEtaMin && particle.eta() < trackEtaMax) {
        registry.fill(HIST("h_partpT"), particle.pt(), weightEvt);
        registry.fill(HIST("h_partpT_fine"), particle.pt(), weightEvt);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCDTracks, "track information in MCD", false);

  PresliceUnsorted<o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>> perFoundBC = aod::evsel::foundBCId;

  void processBCs(soa::Join<aod::BCs, aod::BcSels> const& bcs, soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    if (bcs.size() == 0) {
      return;
    }
    for (const auto& bc : bcs) {
      registry.fill(HIST("hBCCounter"), 0.5); // All BC
      if (bc.selection_bit(aod::evsel::kIsTriggerTVX)) {
        registry.fill(HIST("hBCCounter"), 1.5); // BC+TVX
        if (bc.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
          registry.fill(HIST("hBCCounter"), 2.5); // BC+TVX+NoTFB
          if (bc.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
            registry.fill(HIST("hBCCounter"), 3.5); // BC+TVX+NoTFB+NoITSROFB ----> this goes to Lumi i.e. hLumiAfterBCcuts in eventSelection task
          }
        }
      }
      auto collisionsInBC = collisions.sliceBy(perFoundBC, bc.globalIndex());
      for (const auto& collision : collisionsInBC) {
        registry.fill(HIST("hBCCounter"), 4.5); // CollinBC
        if (collision.sel8()) {
          registry.fill(HIST("hBCCounter"), 5.5); // CollinBC+sel8
          if (std::fabs(collision.posZ()) < vertexZCut) {
            registry.fill(HIST("hBCCounter"), 6.5); // CollinBC+sel8+VtxZ
          }
          if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
            registry.fill(HIST("hBCCounter"), 7.5); // CollinBC+sel8Full
            if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
              registry.fill(HIST("hBCCounter"), 8.5); // CollinBC+sel8Full+GoodZvtx
              if (std::fabs(collision.posZ()) < vertexZCut) {
                registry.fill(HIST("hBCCounter"), 9.5); // CollinBC+sel8Full+VtxZ+GoodZvtx ----> this goes to my analysis task for jet events selection
              }
            }
          }
        }
      } // collision loop
    } // bc loop
  }
  PROCESS_SWITCH(BjetTaggingGnn, processBCs, "BCs for 0 vertex QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetTaggingGnn>(cfgc)};
}
