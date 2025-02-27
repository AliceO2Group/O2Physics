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

#include <algorithm>
#include <cmath>
#include <vector>
#include <string>

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/DataModel/JetTagging.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"

#include "Framework/Logger.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BjetTaggingGnn {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  Configurable<float> pTHatMaxMCD{"pTHatMaxMCD", 999.0, "maximum fraction of hard scattering for jet acceptance in detector MC"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum fraction of hard scattering for jet acceptance in particle MC"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the calculation of pTHat"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.5, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};

  Configurable<float> trackNppCrit{"trackNppCrit", 0.95, "track not physical primary ratio"};

  // track level configurables
  Configurable<float> svPtMin{"svPtMin", 0.5, "minimum SV pT"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -99.0, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 99.0, "maximum jet pseudorapidity"};

  Configurable<std::vector<double>> jetRadii{"jetRadii", std::vector<double>{0.4}, "jet resolution parameters"};

  Configurable<bool> doDataDriven{"doDataDriven", false, "Flag whether to use fill THnSpase for data driven methods"};
  Configurable<bool> callSumw2{"callSumw2", false, "Flag whether to call THnSparse::Sumw2() for error calculation"};

  std::vector<int> eventSelectionBits;

  std::vector<double> jetRadiiValues;

  void init(InitContext const&)
  {
    jetRadiiValues = (std::vector<double>)jetRadii;

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    registry.add("h_vertexZ", "Vertex Z;#it{Z} (cm)", {HistType::kTH1F, {{40, -20.0, 20.0}}});

    const AxisSpec axisJetpT{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisDb{200, -10., 20., "#it{D}_{b}"};
    const AxisSpec axisDbFine{3000, -10., 20., "#it{D}_{b}"};
    const AxisSpec axisSVMass{200, 0., 10., "#it{m}_{SV} (GeV/#it{c}^{2})"};
    const AxisSpec axisSVEnergy{200, 0., 100., "#it{E}_{SV} (GeV)"};
    const AxisSpec axisSLxy{200, 0., 100., "#it{SL}_{xy}"};
    const AxisSpec axisJetMass{200, 0., 50., "#it{m}_{jet} (GeV/#it{c}^{2})"};
    const AxisSpec axisJetProb{200, 0., 40., "-ln(JP)"};
    const AxisSpec axisNTracks{42, 0, 42, "#it{n}_{tracks}"};

    registry.add("h_jetpT", "", {HistType::kTH1F, {axisJetpT}});
    registry.add("h_Db", "", {HistType::kTH1F, {axisDbFine}});
    registry.add("h2_jetpT_Db", "", {HistType::kTH2F, {axisJetpT, axisDb}});
    registry.add("h2_jetpT_SVMass", "", {HistType::kTH2F, {axisJetpT, axisSVMass}});
    registry.add("h2_jetpT_jetMass", "", {HistType::kTH2F, {axisJetpT, axisJetMass}});
    registry.add("h2_jetpT_jetProb", "", {HistType::kTH2F, {axisJetpT, axisJetProb}});
    registry.add("h2_jetpT_nTracks", "", {HistType::kTH2F, {axisJetpT, axisNTracks}});

    if (doprocessMCJets) {
      registry.add("h_jetpT_b", "b-jet", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_jetpT_c", "c-jet", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_jetpT_lf", "lf-jet", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_Db_b", "b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_c", "c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_lf", "lf-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h2_jetpT_Db_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_SVMass_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisSVMass}});
      registry.add("h2_jetpT_SVMass_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisSVMass}});
      registry.add("h2_jetpT_SVMass_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisSVMass}});
      registry.add("h2_jetpT_jetMass_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetMass}});
      registry.add("h2_jetpT_jetMass_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisJetMass}});
      registry.add("h2_jetpT_jetMass_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisJetMass}});
      registry.add("h2_jetpT_jetProb_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetProb}});
      registry.add("h2_jetpT_jetProb_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisJetProb}});
      registry.add("h2_jetpT_jetProb_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisJetProb}});
      registry.add("h2_jetpT_nTracks_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisNTracks}});
      registry.add("h2_jetpT_nTracks_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisNTracks}});
      registry.add("h2_jetpT_nTracks_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisNTracks}});
      registry.add("h2_Response_DetjetpT_PartjetpT", "", {HistType::kTH2F, {axisJetpT, axisJetpT}});
      registry.add("h2_Response_DetjetpT_PartjetpT_b", "b-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}});
      registry.add("h2_Response_DetjetpT_PartjetpT_c", "c-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}});
      registry.add("h2_Response_DetjetpT_PartjetpT_lf", "lf-jet", {HistType::kTH2F, {axisJetpT, axisJetpT}});
      registry.add("h2_jetpT_Db_lf_none", "lf-jet (none)", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_lf_matched", "lf-jet (matched)", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp", "NotPhysPrim", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h2_jetpT_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH2F, {axisJetpT, axisDb}});
      registry.add("h_Db_npp", "NotPhysPrim", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_b", "NotPhysPrim b-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_c", "NotPhysPrim c-jet", {HistType::kTH1F, {axisDbFine}});
      registry.add("h_Db_npp_lf", "NotPhysPrim lf-jet", {HistType::kTH1F, {axisDbFine}});
      // registry.add("h2_pT_dcaXY_pp", "tracks", {HistType::kTH2F, {axisJetpT, {200, 0., 1.}}});
      // registry.add("h2_pT_dcaXY_npp", "NotPhysPrim tracks", {HistType::kTH2F, {axisJetpT, {200, 0., 1.}}});
      // registry.add("h2_pT_dcaZ_pp", "tracks", {HistType::kTH2F, {axisJetpT, {200, 0., 2.}}});
      // registry.add("h2_pT_dcaZ_npp", "NotPhysPrim tracks", {HistType::kTH2F, {axisJetpT, {200, 0., 2.}}});
    }

    if (doprocessMCTruthJets) {
      registry.add("h_jetpT_particle", "", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_jetpT_particle_b", "particle b-jet", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_jetpT_particle_c", "particle c-jet", {HistType::kTH1F, {axisJetpT}});
      registry.add("h_jetpT_particle_lf", "particle lf-jet", {HistType::kTH1F, {axisJetpT}});
    }

    if (doDataDriven) {
      registry.add("hSparse_Incljets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
      if (doprocessMCJets) {
        registry.add("hSparse_bjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
        registry.add("hSparse_cjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets_none", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
        registry.add("hSparse_lfjets_matched", "", {HistType::kTHnSparseF, {axisJetpT, axisDbFine, axisSVMass, axisJetMass, axisNTracks}}, callSumw2);
      }
    }
  }

  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackCuts = (aod::jtrack::pt > trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt <= jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  using FilteredCollision = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs>>;
  using DataJets = soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents, aod::ChargedJetTags, aod::DataSecondaryVertex3ProngIndices>>;
  using JetTrackswID = soa::Filtered<soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>>;
  using JetTracksMCDwID = soa::Filtered<soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>>;
  using SVTable = aod::DataSecondaryVertex3Prongs;
  using MCDSVTable = aod::MCDSecondaryVertex3Prongs;

  template <typename AnyCollision, typename AnalysisJet, typename AnyTracks>
  int analyzeJetTrackInfo(AnyCollision const& /*collision*/, AnalysisJet const& analysisJet, AnyTracks const& /*allTracks*/ /*, int8_t jetFlavor = 0, double weight = 1.0*/)
  {
    int nTracks = 0;
    for (const auto& constituent : analysisJet.template tracks_as<AnyTracks>()) {

      if (constituent.pt() < trackPtMin) {
        continue;
      }

      // ...

      ++nTracks;
    }
    return nTracks;
  }

  template <typename AnalysisJet, typename SecondaryVertices>
  SecondaryVertices::iterator analyzeJetSVInfo(AnalysisJet const& analysisJet, SecondaryVertices const& allSVs, bool& checkSV /*, int8_t jetFlavor = 0, double weight = 1.0*/)
  {
    using SVType = typename SecondaryVertices::iterator;

    auto compare = [](SVType& sv1, SVType& sv2) {
      return (sv1.decayLengthXY() / sv1.errorDecayLengthXY()) > (sv2.decayLengthXY() / sv2.errorDecayLengthXY());
    };

    auto svs = analysisJet.template secondaryVertices_as<SecondaryVertices>();

    std::sort(svs.begin(), svs.end(), compare);

    checkSV = false;
    for (const auto& candSV : svs) {

      if (candSV.pt() < svPtMin) {
        continue;
      }

      checkSV = true;
      return candSV;
    }

    // No SV found
    return *allSVs.begin();
  }

  void processDummy(FilteredCollision::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDummy, "Dummy process function turned on by default", true);

  void processDataJets(FilteredCollision::iterator const& collision, DataJets const& alljets, JetTrackswID const& allTracks, SVTable const& allSVs)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

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

      float mSV = -1.f;
      // float eSV = -1.f;
      // float slXY = -1.f;

      bool checkSV;
      // auto sv = jettaggingutilities::jetFromProngMaxDecayLength<MCDSVTable>(analysisJet, prongChi2PCAMin, prongChi2PCAMax, prongsigmaLxyMax, prongIPxyMin, prongIPxyMax, false, &checkSV);
      auto sv = analyzeJetSVInfo(analysisJet, allSVs, checkSV);

      if (checkSV) {
        mSV = sv.m();
        // eSV = sv.e();
        // slXY = sv.decayLengthXY() / sv.errorDecayLengthXY();
      }

      registry.fill(HIST("h_jetpT"), analysisJet.pt());
      registry.fill(HIST("h_Db"), analysisJet.scoreML());
      registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML());
      registry.fill(HIST("h2_jetpT_SVMass"), analysisJet.pt(), mSV);
      registry.fill(HIST("h2_jetpT_jetMass"), analysisJet.pt(), analysisJet.mass());
      registry.fill(HIST("h2_jetpT_jetProb"), analysisJet.pt(), analysisJet.jetProb());
      registry.fill(HIST("h2_jetpT_nTracks"), analysisJet.pt(), nTracks);

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processDataJets, "jet information in Data", false);

  using MCDJetTable = soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags, aod::ChargedMCDetectorLevelJetEventWeights, aod::MCDSecondaryVertex3ProngIndices>>;
  using MCPJetTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef, aod::ChargedMCParticleLevelJetEventWeights>>;
  using FilteredCollisionMCD = soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::JMcCollisionLbs>>;

  void processMCJets(FilteredCollisionMCD::iterator const& collision, MCDJetTable const& MCDjets, MCPJetTable const& /*MCPjets*/, JetTracksMCDwID const& /*allTracks*/, MCDSVTable const& allSVs, aod::JetParticles const& /*MCParticles*/)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    registry.fill(HIST("h_vertexZ"), collision.posZ());

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

      float weight = useEventWeight ? analysisJet.eventWeight() : 1.f;
      float pTHat = 10. / (std::pow(analysisJet.eventWeight(), 1.0 / pTHatExponent));
      if (analysisJet.pt() > pTHatMaxMCD * pTHat) {
        continue;
      }

      int8_t jetFlavor = analysisJet.origin();

      // int nTracks = analyzeJetTrackInfo(collision, analysisJet, allTracks /*, jetFlavor, weight*/);
      int nTracks = 0;

      int nNppTracks = 0;
      for (const auto& constituent : analysisJet.template tracks_as<JetTracksMCDwID>()) {
        if (constituent.pt() < trackPtMin) {
          continue;
        }
        if (!constituent.has_mcParticle() || !constituent.template mcParticle_as<aod::JetParticles>().isPhysicalPrimary()) {
          // registry.fill(HIST("h2_pT_dcaXY_npp"), constituent.pt(), constituent.dcaXY());
          // registry.fill(HIST("h2_pT_dcaZ_npp"), constituent.pt(), constituent.dcaZ());
          ++nNppTracks;
        } else {
          // registry.fill(HIST("h2_pT_dcaXY_pp"), constituent.pt(), constituent.dcaXY());
          // registry.fill(HIST("h2_pT_dcaZ_pp"), constituent.pt(), constituent.dcaZ());
        }
        ++nTracks;
      }

      float mSV = -1.f;
      // float eSV = -1.f;
      // float slXY = -1.f;

      bool checkSV;
      auto sv = analyzeJetSVInfo(analysisJet, allSVs, checkSV /*, jetFlavor, weight*/);

      if (checkSV) {
        mSV = sv.m();
        // eSV = sv.e();
        // slXY = sv.decayLengthXY() / sv.errorDecayLengthXY();
      }

      registry.fill(HIST("h_jetpT"), analysisJet.pt(), weight);
      registry.fill(HIST("h_Db"), analysisJet.scoreML(), weight);
      registry.fill(HIST("h2_jetpT_Db"), analysisJet.pt(), analysisJet.scoreML(), weight);
      registry.fill(HIST("h2_jetpT_SVMass"), analysisJet.pt(), mSV, weight);
      registry.fill(HIST("h2_jetpT_jetMass"), analysisJet.pt(), analysisJet.mass(), weight);
      registry.fill(HIST("h2_jetpT_jetProb"), analysisJet.pt(), analysisJet.jetProb(), weight);
      registry.fill(HIST("h2_jetpT_nTracks"), analysisJet.pt(), nTracks, weight);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_b"), analysisJet.pt(), weight);
        registry.fill(HIST("h_Db_b"), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_Db_b"), analysisJet.pt(), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_SVMass_b"), analysisJet.pt(), mSV, weight);
        registry.fill(HIST("h2_jetpT_jetMass_b"), analysisJet.pt(), analysisJet.mass(), weight);
        registry.fill(HIST("h2_jetpT_jetProb_b"), analysisJet.pt(), analysisJet.jetProb(), weight);
        registry.fill(HIST("h2_jetpT_nTracks_b"), analysisJet.pt(), nTracks, weight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_c"), analysisJet.pt(), weight);
        registry.fill(HIST("h_Db_c"), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_Db_c"), analysisJet.pt(), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_SVMass_c"), analysisJet.pt(), mSV, weight);
        registry.fill(HIST("h2_jetpT_jetMass_c"), analysisJet.pt(), analysisJet.mass(), weight);
        registry.fill(HIST("h2_jetpT_jetProb_c"), analysisJet.pt(), analysisJet.jetProb(), weight);
        registry.fill(HIST("h2_jetpT_nTracks_c"), analysisJet.pt(), nTracks, weight);
      } else {
        registry.fill(HIST("h_jetpT_lf"), analysisJet.pt(), weight);
        registry.fill(HIST("h_Db_lf"), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_Db_lf"), analysisJet.pt(), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_SVMass_lf"), analysisJet.pt(), mSV, weight);
        registry.fill(HIST("h2_jetpT_jetMass_lf"), analysisJet.pt(), analysisJet.mass(), weight);
        registry.fill(HIST("h2_jetpT_jetProb_lf"), analysisJet.pt(), analysisJet.jetProb(), weight);
        registry.fill(HIST("h2_jetpT_nTracks_lf"), analysisJet.pt(), nTracks, weight);
        if (jetFlavor == JetTaggingSpecies::none) {
          registry.fill(HIST("h2_jetpT_Db_lf_none"), analysisJet.pt(), analysisJet.scoreML(), weight);
        } else {
          registry.fill(HIST("h2_jetpT_Db_lf_matched"), analysisJet.pt(), analysisJet.scoreML(), weight);
        }
      }

      if (static_cast<float>(nNppTracks) / nTracks > trackNppCrit) {
        registry.fill(HIST("h_Db_npp"), analysisJet.scoreML(), weight);
        registry.fill(HIST("h2_jetpT_Db_npp"), analysisJet.pt(), analysisJet.scoreML(), weight);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h_Db_npp_b"), analysisJet.scoreML(), weight);
          registry.fill(HIST("h2_jetpT_Db_npp_b"), analysisJet.pt(), analysisJet.scoreML(), weight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h_Db_npp_c"), analysisJet.scoreML(), weight);
          registry.fill(HIST("h2_jetpT_Db_npp_c"), analysisJet.pt(), analysisJet.scoreML(), weight);
        } else {
          registry.fill(HIST("h_Db_npp_lf"), analysisJet.scoreML(), weight);
          registry.fill(HIST("h2_jetpT_Db_npp_lf"), analysisJet.pt(), analysisJet.scoreML(), weight);
        }
      }

      if (doDataDriven) {
        registry.fill(HIST("hSparse_Incljets"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("hSparse_bjets"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("hSparse_cjets"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
        } else {
          registry.fill(HIST("hSparse_lfjets"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
          if (jetFlavor == JetTaggingSpecies::none) {
            registry.fill(HIST("hSparse_lfjets_none"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
          } else {
            registry.fill(HIST("hSparse_lfjets_matched"), analysisJet.pt(), analysisJet.scoreML(), mSV, analysisJet.mass(), nTracks, weight);
          }
        }
      }

      for (const auto& mcpjet : analysisJet.template matchedJetGeo_as<MCPJetTable>()) {
        if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
          continue;
        }

        registry.fill(HIST("h2_Response_DetjetpT_PartjetpT"), analysisJet.pt(), mcpjet.pt(), weight);
        if (jetFlavor == JetTaggingSpecies::beauty) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_b"), analysisJet.pt(), mcpjet.pt(), weight);
        } else if (jetFlavor == JetTaggingSpecies::charm) {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_c"), analysisJet.pt(), mcpjet.pt(), weight);
        } else {
          registry.fill(HIST("h2_Response_DetjetpT_PartjetpT_lf"), analysisJet.pt(), mcpjet.pt(), weight);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCJets, "jet information in MC", false);

  Filter mccollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  using FilteredCollisionMCP = soa::Filtered<aod::JMcCollisions>;

  void processMCTruthJets(FilteredCollisionMCP::iterator const& /*collision*/, MCPJetTable const& MCPjets, aod::JetParticles const& /*MCParticles*/)
  {

    for (const auto& mcpjet : MCPjets) {

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

      float weight = useEventWeight ? mcpjet.eventWeight() : 1.0;
      float pTHat = 10. / (std::pow(mcpjet.eventWeight(), 1.0 / pTHatExponent));
      if (mcpjet.pt() > pTHatMaxMCP * pTHat) {
        continue;
      }

      int8_t jetFlavor = mcpjet.origin();

      registry.fill(HIST("h_jetpT_particle"), mcpjet.pt(), weight);

      if (jetFlavor == JetTaggingSpecies::beauty) {
        registry.fill(HIST("h_jetpT_particle_b"), mcpjet.pt(), weight);
      } else if (jetFlavor == JetTaggingSpecies::charm) {
        registry.fill(HIST("h_jetpT_particle_c"), mcpjet.pt(), weight);
      } else {
        registry.fill(HIST("h_jetpT_particle_lf"), mcpjet.pt(), weight);
      }
    }
  }
  PROCESS_SWITCH(BjetTaggingGnn, processMCTruthJets, "truth jet information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetTaggingGnn>(cfgc)};
}
