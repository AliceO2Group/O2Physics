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

/// \file   bjetSubstructureGnn.cxx
/// \brief  b-jet fragmentation and groomed substructure analysis using GNN tagging
///
/// \author Yungwang Hwang <yungwang.hwang@cern.ch>, Yonsei University

#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/Core/JetTaggingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/JetSubstructure.h"
#include "PWGJE/DataModel/JetTagging.h"

#include "Common/DataModel/EventSelection.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisDataModelHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct BjetSubstructureGnn {

  HistogramRegistry registry;

  // event level configurables
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<bool> useEventWeight{"useEventWeight", true, "Flag whether to scale histograms with the event weight"};

  // track level configurables
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum track pT"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum track pT"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum track eta"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum track eta"};

  // jet level configurables
  Configurable<float> jetPtMin{"jetPtMin", 0.0, "minimum jet pT"};
  Configurable<float> jetPtMax{"jetPtMax", 1000.0, "maximum jet pT"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.9, "minimum jet eta"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.9, "maximum jet eta"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};

  Configurable<double> dbMin{"dbMin", -15.0, "minimum GNN db"};
  Configurable<double> dbMax{"dbMax", 15.0, "maximum GNN db"};
  Configurable<int> dbNbins{"dbNbins", 300, "number of bins in axisDbFine"};

  Configurable<bool> callSumw2{"callSumw2", true, "Flag whether to call THnSparse::Sumw2() for error calculation"};

  // fragmentation observables
  ConfigurableAxis binFragmentationZ{"binFragmentationZ", {VARIABLE_WIDTH, 0.0, 0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060, 0.080, 0.100, 0.125, 0.150, 0.200, 0.300, 0.400, 0.500, 0.750, 1.000, 1.001}, "binning for longitudinal momentum fraction z"};
  ConfigurableAxis binJt{"binJt", {VARIABLE_WIDTH, 0., 0.005, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.625, 0.750, 1.0, 1.25, 1.5, 2.0, 3.5}, "binning for transverse momentum to jet axis jT"};
  ConfigurableAxis binTrackJetR{"binTrackJetR", {10, 0.0, 0.5}, "binning for track-jet radial distance r"};

  // groomed substructure observables
  Configurable<float> softDropZCut{"softDropZCut", 0.1f, "Soft Drop Z cut"};
  Configurable<float> softDropBeta{"softDropBeta", 0.f, "Soft Drop beta"};
  ConfigurableAxis binZg{"binZg", {20, 0.0, 0.5}, "binning for groomed momentum sharing zg"};
  ConfigurableAxis binRg{"binRg", {20, 0.0, 0.5}, "binning for groomed radius Rg"};
  ConfigurableAxis binNSD{"binNSD", {15, -0.5, 14.5}, "binning for number of Soft Drop splittings"};

  // event selection bits
  std::vector<int> eventSelectionBits;

  // Jet flavor categories: b, c, lf, none
  // MCP flavour comes from processMCPByDistance, which does not return "none".
  static constexpr std::size_t kNFlavors = 4;

  // MCP histograms indexed
  std::array<std::shared_ptr<TH1>, kNFlavors> hMCPFragmentationJetPt;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPZ;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPJt;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPR;

  std::array<std::shared_ptr<TH1>, kNFlavors> hMCPGroomedJetPt;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPZg;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPRg;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCPNSD;

  // MCD histograms indexed
  std::array<std::shared_ptr<TH1>, kNFlavors> hMCDFragmentationJetPt;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCDFragmentationJetPtDb;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDZ;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDJt;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDR;

  std::array<std::shared_ptr<TH1>, kNFlavors> hMCDGroomedJetPt;
  std::array<std::shared_ptr<TH2>, kNFlavors> hMCDGroomedJetPtDb;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDZg;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDRg;
  std::array<std::shared_ptr<TH3>, kNFlavors> hMCDNSD;

  void init(InitContext const&)
  {

    eventSelectionBits = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));

    const AxisSpec axisJetpT{200, 0., 200., "#it{p}_{T, ch jet} (GeV/#it{c})"};
    const AxisSpec axisJetEta{200, -0.8, 0.8, "#it{#eta}_{jet}"};
    const AxisSpec axisJetPhi{200, 0., o2::constants::math::TwoPI, "#it{phi}_{jet}"};
    const AxisSpec axisDb{dbNbins, dbMin, dbMax, "#it{D}_{b}"};
    const AxisSpec axisJetFlavorCat{4, -0.5, 3.5, "flavour"};
    const AxisSpec axisNConstituents{101, -0.5, 100.5, "#it{N}_{constituents}"};
    const AxisSpec axisSoftDropStatus{2, -0.5, 1.5, "Soft Drop status"};

    const AxisSpec axisFragmentationZ{binFragmentationZ, "#it{z}"};
    const AxisSpec axisJt{binJt, "#it{j}_{T} (GeV/#it{c})"};
    const AxisSpec axisTrackJetR{binTrackJetR, "#it{r}"};

    const AxisSpec axisZg{binZg, "#it{z}_g"};
    const AxisSpec axisRg{binRg, "#it{R}_{g}"};
    const AxisSpec axisNSD{binNSD, "#it{n}_{SD}"};

    const std::array<std::string, kNFlavors> flavorNames{"b", "c", "lf", "none"};

    auto addJetQAHistograms = [&](const std::string& basePath) {
      registry.add(basePath + "/h_jetPt", "", HistType::kTH1F, {axisJetpT}, callSumw2);
      registry.add(basePath + "/h_jetEta", "", HistType::kTH1F, {axisJetEta}, callSumw2);
      registry.add(basePath + "/h_jetPhi", "", HistType::kTH1F, {axisJetPhi}, callSumw2);
      registry.add(basePath + "/h_event_counter", "", HistType::kTH1F, {{1, 0.0, 1.0}}, callSumw2);
      registry.add(basePath + "/h_vertexZ", "", HistType::kTH1F, {{100, -20.0, 20.0}}, callSumw2);
      auto hJetFlavor = registry.add<TH1>(basePath + "/h_jetFlavor", "", HistType::kTH1F, {axisJetFlavorCat}, callSumw2);
      setJetFlavorCatAxisLabels(hJetFlavor->GetXaxis());
    };

    if (doprocessMCDFragmentation) {
      addJetQAHistograms("fragmentation/mcd");
      registry.add("fragmentation/mcd/h2_jetPt_nConstituents", "", HistType::kTH2F, {axisJetpT, axisNConstituents}, callSumw2);

      for (std::size_t flav = 0; flav < kNFlavors; flav++) {
        const auto& flavor = flavorNames[flav];

        hMCDFragmentationJetPt[flav] = registry.add<TH1>("fragmentation/mcd/" + flavor + "/h_jetPt", "", HistType::kTH1F, {axisJetpT}, callSumw2);
        hMCDFragmentationJetPtDb[flav] = registry.add<TH2>("fragmentation/mcd/" + flavor + "/h2_jetPt_Db", "", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
        hMCDZ[flav] = registry.add<TH3>("fragmentation/mcd/" + flavor + "/h3_jetPt_z_Db", "", HistType::kTH3F, {axisJetpT, axisFragmentationZ, axisDb}, callSumw2);
        hMCDJt[flav] = registry.add<TH3>("fragmentation/mcd/" + flavor + "/h3_jetPt_jT_Db", "", HistType::kTH3F, {axisJetpT, axisJt, axisDb}, callSumw2);
        hMCDR[flav] = registry.add<TH3>("fragmentation/mcd/" + flavor + "/h3_jetPt_r_Db", "", HistType::kTH3F, {axisJetpT, axisTrackJetR, axisDb}, callSumw2);
      }

      registry.add("fragmentation/mcd/h2_jetPt_Db", "", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
      registry.add("fragmentation/mcd/h3_jetPt_z_Db", "", {HistType::kTH3F, {axisJetpT, axisFragmentationZ, axisDb}}, callSumw2);
      registry.add("fragmentation/mcd/h3_jetPt_r_Db", "", {HistType::kTH3F, {axisJetpT, axisTrackJetR, axisDb}}, callSumw2);
      registry.add("fragmentation/mcd/h3_jetPt_jT_Db", "", {HistType::kTH3F, {axisJetpT, axisJt, axisDb}}, callSumw2);
    }

    if (doprocessMCPFragmentation) {
      addJetQAHistograms("fragmentation/mcp");
      registry.add("fragmentation/mcp/h2_jetPt_nConstituents", "", HistType::kTH2F, {axisJetpT, axisNConstituents}, callSumw2);

      for (std::size_t flav = 0; flav < kNFlavors; flav++) {
        const auto& flavor = flavorNames[flav];

        hMCPFragmentationJetPt[flav] = registry.add<TH1>("fragmentation/mcp/" + flavor + "/h_jetPt", "", HistType::kTH1F, {axisJetpT}, callSumw2);
        hMCPJt[flav] = registry.add<TH2>("fragmentation/mcp/" + flavor + "/h2_jetPt_jT", "", HistType::kTH2F, {axisJetpT, axisJt}, callSumw2);
        hMCPZ[flav] = registry.add<TH2>("fragmentation/mcp/" + flavor + "/h2_jetPt_z", "", HistType::kTH2F, {axisJetpT, axisFragmentationZ}, callSumw2);
        hMCPR[flav] = registry.add<TH2>("fragmentation/mcp/" + flavor + "/h2_jetPt_r", "", HistType::kTH2F, {axisJetpT, axisTrackJetR}, callSumw2);
      }

      registry.add("fragmentation/mcp/h2_jetPt_z", "", HistType::kTH2F, {axisJetpT, axisFragmentationZ}, callSumw2);
      registry.add("fragmentation/mcp/h2_jetPt_r", "", HistType::kTH2F, {axisJetpT, axisTrackJetR}, callSumw2);
      registry.add("fragmentation/mcp/h2_jetPt_jT", "", HistType::kTH2F, {axisJetpT, axisJt}, callSumw2);
    }

    if (doprocessMCDGroomedSubstructure) {
      addJetQAHistograms("groomed/mcd");
      auto hSoftDropStatus = registry.add<TH2>("groomed/mcd/h2_jetPt_softDropStatus", "", HistType::kTH2F, {axisJetpT, axisSoftDropStatus}, callSumw2);
      hSoftDropStatus->GetYaxis()->SetBinLabel(1, "no accepted splitting");
      hSoftDropStatus->GetYaxis()->SetBinLabel(2, "accepted splitting");

      for (std::size_t flav = 0; flav < kNFlavors; flav++) {
        const auto& flavor = flavorNames[flav];

        hMCDGroomedJetPt[flav] = registry.add<TH1>("groomed/mcd/" + flavor + "/h_jetPt", "", HistType::kTH1F, {axisJetpT}, callSumw2);
        hMCDGroomedJetPtDb[flav] = registry.add<TH2>("groomed/mcd/" + flavor + "/h2_jetPt_Db", "", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
        hMCDZg[flav] = registry.add<TH3>("groomed/mcd/" + flavor + "/h3_jetPt_zg_Db", "", HistType::kTH3F, {axisJetpT, axisZg, axisDb}, callSumw2);
        hMCDRg[flav] = registry.add<TH3>("groomed/mcd/" + flavor + "/h3_jetPt_rg_Db", "", HistType::kTH3F, {axisJetpT, axisRg, axisDb}, callSumw2);
        hMCDNSD[flav] = registry.add<TH3>("groomed/mcd/" + flavor + "/h3_jetPt_nSD_Db", "", HistType::kTH3F, {axisJetpT, axisNSD, axisDb}, callSumw2);
      }

      registry.add("groomed/mcd/h2_jetPt_Db", "", HistType::kTH2F, {axisJetpT, axisDb}, callSumw2);
      registry.add("groomed/mcd/h3_jetPt_zg_Db", "", HistType::kTH3F, {axisJetpT, axisZg, axisDb}, callSumw2);
      registry.add("groomed/mcd/h3_jetPt_rg_Db", "", HistType::kTH3F, {axisJetpT, axisRg, axisDb}, callSumw2);
      registry.add("groomed/mcd/h3_jetPt_nSD_Db", "", HistType::kTH3F, {axisJetpT, axisNSD, axisDb}, callSumw2);
    }

    if (doprocessMCPGroomedSubstructure) {
      addJetQAHistograms("groomed/mcp");
      auto hSoftDropStatus = registry.add<TH2>("groomed/mcp/h2_jetPt_softDropStatus", "", HistType::kTH2F, {axisJetpT, axisSoftDropStatus}, callSumw2);
      hSoftDropStatus->GetYaxis()->SetBinLabel(1, "no accepted splitting");
      hSoftDropStatus->GetYaxis()->SetBinLabel(2, "accepted splitting");

      for (std::size_t flav = 0; flav < kNFlavors; flav++) {
        const auto& flavor = flavorNames[flav];

        hMCPGroomedJetPt[flav] = registry.add<TH1>("groomed/mcp/" + flavor + "/h_jetPt", "", HistType::kTH1F, {axisJetpT}, callSumw2);
        hMCPZg[flav] = registry.add<TH2>("groomed/mcp/" + flavor + "/h2_jetPt_zg", "", HistType::kTH2F, {axisJetpT, axisZg}, callSumw2);
        hMCPRg[flav] = registry.add<TH2>("groomed/mcp/" + flavor + "/h2_jetPt_rg", "", HistType::kTH2F, {axisJetpT, axisRg}, callSumw2);
        hMCPNSD[flav] = registry.add<TH2>("groomed/mcp/" + flavor + "/h2_jetPt_nSD", "", HistType::kTH2F, {axisJetpT, axisNSD}, callSumw2);
      }

      registry.add("groomed/mcp/h2_jetPt_zg", "", HistType::kTH2F, {axisJetpT, axisZg}, callSumw2);
      registry.add("groomed/mcp/h2_jetPt_rg", "", HistType::kTH2F, {axisJetpT, axisRg}, callSumw2);
      registry.add("groomed/mcp/h2_jetPt_nSD", "", HistType::kTH2F, {axisJetpT, axisNSD}, callSumw2);
    }
  }

  Filter collisionFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter trackFilter = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter jetFilter = (aod::jet::pt >= jetPtMin && aod::jet::pt < jetPtMax && aod::jet::eta < jetEtaMax - aod::jet::r / 100.f && aod::jet::eta > jetEtaMin + aod::jet::r / 100.f);

  // MCD
  using AnalysisCollisionsMCD = soa::Join<aod::JetCollisionsMCD, aod::JCollisionOutliers>;
  using FilteredCollisionsMCD = soa::Filtered<AnalysisCollisionsMCD>;

  using MCDFragmentationJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags>;
  using MCDGroomedJets = soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets, aod::ChargedMCDetectorLevelJetFlavourDef, aod::ChargedMCDetectorLevelJetTags, aod::CMCDJetSSs>;

  using FilteredMCDFragmentationJets = soa::Filtered<MCDFragmentationJets>;
  using FilteredMCDGroomedJets = soa::Filtered<MCDGroomedJets>;

  using AnalysisTracksMCD = aod::JetTracksMCD;

  // Apply the task-level pT and eta selection to detector-level fragmentation constituents.
  using FilteredTracksMCD = soa::Filtered<AnalysisTracksMCD>;

  // MCP
  using AnalysisCollisionsMCP = soa::Join<aod::JetMcCollisions, aod::JMcCollisionPIs, aod::JMcCollisionOutliers>;

  using MCPFragmentationJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef>;
  using MCPGroomedJets = soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets, aod::ChargedMCParticleLevelJetFlavourDef, aod::CMCPJetSSs>;

  using FilteredMCPFragmentationJets = soa::Filtered<MCPFragmentationJets>;
  using FilteredMCPGroomedJets = soa::Filtered<MCPGroomedJets>;

  // Retrieve particle-level jets associated with each MC collision.
  Preslice<FilteredMCPFragmentationJets> mcpFragmentationJetsPerMCPCollision = aod::jet::mcCollisionId;
  Preslice<FilteredMCPGroomedJets> mcpGroomedJetsPerMCPCollision = aod::jet::mcCollisionId;

  struct FragmentationObservables {
    float z;
    float jT;
    float r;
  };

  struct GroomedObservables {
    float zg{-1.f};
    float rg{-1.f};
    int nSD{0};
    bool hasAcceptedSplitting{false};
  };

  template <typename AnalysisJet, typename Constituent>
  FragmentationObservables calculateFragmentation(AnalysisJet const& jet, Constituent const& constituent)
  {
    const float jetPt = jet.pt();

    const float z = (constituent.px() * jet.px() + constituent.py() * jet.py()) / (jetPt * jetPt);
    const float jT = std::abs(constituent.px() * jet.py() - constituent.py() * jet.px()) / jetPt;
    const float r = jetutilities::deltaR(jet, constituent);

    return {z, jT, r};
  }

  template <typename AnalysisJet>
  GroomedObservables calculateGroomedSubstructure(AnalysisJet const& jet)
  {
    GroomedObservables result;
    const auto ptLeading = jet.ptLeading();
    const auto ptSubLeading = jet.ptSubLeading();
    const auto theta = jet.theta();

    const float jetR = jet.r() / 100.f;
    if (jetR <= 0.f) {
      return result;
    }

    for (std::size_t i = 0; i < theta.size(); i++) {
      const float sumPt = ptLeading[i] + ptSubLeading[i];
      if (sumPt <= 0.f) {
        continue;
      }

      const float z = ptSubLeading[i] / sumPt;
      const float threshold = softDropZCut * std::pow(theta[i] / jetR, softDropBeta);

      if (z < threshold) {
        continue;
      }

      if (!result.hasAcceptedSplitting) {
        result.zg = z;
        result.rg = theta[i];
        result.hasAcceptedSplitting = true;
      }
      result.nSD++;
    }
    return result;
  }

  static std::size_t getJetFlavorCat(int8_t jetFlavor)
  {
    if (jetFlavor == JetTaggingSpecies::beauty) {
      return 0;
    }
    if (jetFlavor == JetTaggingSpecies::charm) {
      return 1;
    }
    if (jetFlavor == JetTaggingSpecies::none) {
      return 3;
    }
    return 2;
  }

  static void setJetFlavorCatAxisLabels(TAxis* axis)
  {
    axis->SetBinLabel(1, "b");
    axis->SetBinLabel(2, "c");
    axis->SetBinLabel(3, "lf");
    axis->SetBinLabel(4, "none");
  }

  template <typename Jet>
  void fillMCDFragmentationJetHistograms(Jet const& mcdjet, std::size_t flavorIndex, float Db, float eventWeight)
  {
    registry.fill(HIST("fragmentation/mcd/h_jetFlavor"), static_cast<double>(flavorIndex), eventWeight);
    registry.fill(HIST("fragmentation/mcd/h2_jetPt_Db"), mcdjet.pt(), Db, eventWeight);
    registry.fill(HIST("fragmentation/mcd/h_jetPt"), mcdjet.pt(), eventWeight);
    registry.fill(HIST("fragmentation/mcd/h_jetEta"), mcdjet.eta(), eventWeight);
    registry.fill(HIST("fragmentation/mcd/h_jetPhi"), mcdjet.phi(), eventWeight);

    hMCDFragmentationJetPt[flavorIndex]->Fill(mcdjet.pt(), eventWeight);
    hMCDFragmentationJetPtDb[flavorIndex]->Fill(mcdjet.pt(), Db, eventWeight);
  }

  template <typename Jet>
  void fillMCDFragmentationHistograms(Jet const& mcdjet, FragmentationObservables const& observable, std::size_t flavorIndex, float Db, float eventWeight)
  {
    registry.fill(HIST("fragmentation/mcd/h3_jetPt_z_Db"), mcdjet.pt(), observable.z, Db, eventWeight);
    registry.fill(HIST("fragmentation/mcd/h3_jetPt_jT_Db"), mcdjet.pt(), observable.jT, Db, eventWeight);
    registry.fill(HIST("fragmentation/mcd/h3_jetPt_r_Db"), mcdjet.pt(), observable.r, Db, eventWeight);

    hMCDZ[flavorIndex]->Fill(mcdjet.pt(), observable.z, Db, eventWeight);
    hMCDJt[flavorIndex]->Fill(mcdjet.pt(), observable.jT, Db, eventWeight);
    hMCDR[flavorIndex]->Fill(mcdjet.pt(), observable.r, Db, eventWeight);
  }

  template <typename Jet>
  void fillMCPFragmentationJetHistograms(Jet const& mcpjet, std::size_t flavorIndex, float eventWeight)
  {
    registry.fill(HIST("fragmentation/mcp/h_jetFlavor"), static_cast<double>(flavorIndex), eventWeight);
    registry.fill(HIST("fragmentation/mcp/h_jetPt"), mcpjet.pt(), eventWeight);
    registry.fill(HIST("fragmentation/mcp/h_jetEta"), mcpjet.eta(), eventWeight);
    registry.fill(HIST("fragmentation/mcp/h_jetPhi"), mcpjet.phi(), eventWeight);

    hMCPFragmentationJetPt[flavorIndex]->Fill(mcpjet.pt(), eventWeight);
  }

  template <typename Jet>
  void fillMCPFragmentationHistograms(Jet const& mcpjet, FragmentationObservables const& observable, std::size_t flavorIndex, float eventWeight)
  {
    registry.fill(HIST("fragmentation/mcp/h2_jetPt_z"), mcpjet.pt(), observable.z, eventWeight);
    registry.fill(HIST("fragmentation/mcp/h2_jetPt_jT"), mcpjet.pt(), observable.jT, eventWeight);
    registry.fill(HIST("fragmentation/mcp/h2_jetPt_r"), mcpjet.pt(), observable.r, eventWeight);

    hMCPZ[flavorIndex]->Fill(mcpjet.pt(), observable.z, eventWeight);
    hMCPJt[flavorIndex]->Fill(mcpjet.pt(), observable.jT, eventWeight);
    hMCPR[flavorIndex]->Fill(mcpjet.pt(), observable.r, eventWeight);
  }

  template <typename Jet>
  void fillMCDGroomedJetHistograms(Jet const& mcdjet, std::size_t flavorIndex, float Db, float eventWeight)
  {
    registry.fill(HIST("groomed/mcd/h_jetFlavor"), static_cast<double>(flavorIndex), eventWeight);
    registry.fill(HIST("groomed/mcd/h2_jetPt_Db"), mcdjet.pt(), Db, eventWeight);
    registry.fill(HIST("groomed/mcd/h_jetPt"), mcdjet.pt(), eventWeight);
    registry.fill(HIST("groomed/mcd/h_jetEta"), mcdjet.eta(), eventWeight);
    registry.fill(HIST("groomed/mcd/h_jetPhi"), mcdjet.phi(), eventWeight);

    hMCDGroomedJetPt[flavorIndex]->Fill(mcdjet.pt(), eventWeight);
    hMCDGroomedJetPtDb[flavorIndex]->Fill(mcdjet.pt(), Db, eventWeight);
  }

  template <typename Jet>
  void fillMCDGroomedSubstructureHistograms(Jet const& mcdjet, GroomedObservables const& observable, std::size_t flavorIndex, float Db, float eventWeight)
  {
    const double softDropStatus = observable.hasAcceptedSplitting ? 1. : 0.;
    registry.fill(HIST("groomed/mcd/h2_jetPt_softDropStatus"), mcdjet.pt(), softDropStatus, eventWeight);
    registry.fill(HIST("groomed/mcd/h3_jetPt_nSD_Db"), mcdjet.pt(), observable.nSD, Db, eventWeight);
    hMCDNSD[flavorIndex]->Fill(mcdjet.pt(), observable.nSD, Db, eventWeight);

    if (!observable.hasAcceptedSplitting) {
      return;
    }

    registry.fill(HIST("groomed/mcd/h3_jetPt_zg_Db"), mcdjet.pt(), observable.zg, Db, eventWeight);
    registry.fill(HIST("groomed/mcd/h3_jetPt_rg_Db"), mcdjet.pt(), observable.rg, Db, eventWeight);
    hMCDZg[flavorIndex]->Fill(mcdjet.pt(), observable.zg, Db, eventWeight);
    hMCDRg[flavorIndex]->Fill(mcdjet.pt(), observable.rg, Db, eventWeight);
  }

  template <typename Jet>
  void fillMCPGroomedJetHistograms(Jet const& mcpjet, std::size_t flavorIndex, float eventWeight)
  {
    registry.fill(HIST("groomed/mcp/h_jetFlavor"), static_cast<double>(flavorIndex), eventWeight);
    registry.fill(HIST("groomed/mcp/h_jetPt"), mcpjet.pt(), eventWeight);
    registry.fill(HIST("groomed/mcp/h_jetEta"), mcpjet.eta(), eventWeight);
    registry.fill(HIST("groomed/mcp/h_jetPhi"), mcpjet.phi(), eventWeight);

    hMCPGroomedJetPt[flavorIndex]->Fill(mcpjet.pt(), eventWeight);
  }

  template <typename Jet>
  void fillMCPGroomedSubstructureHistograms(Jet const& mcpjet, GroomedObservables const& observable, std::size_t flavorIndex, float eventWeight)
  {
    const double softDropStatus = observable.hasAcceptedSplitting ? 1. : 0.;
    registry.fill(HIST("groomed/mcp/h2_jetPt_softDropStatus"), mcpjet.pt(), softDropStatus, eventWeight);
    registry.fill(HIST("groomed/mcp/h2_jetPt_nSD"), mcpjet.pt(), observable.nSD, eventWeight);
    hMCPNSD[flavorIndex]->Fill(mcpjet.pt(), observable.nSD, eventWeight);

    if (!observable.hasAcceptedSplitting) {
      return;
    }

    registry.fill(HIST("groomed/mcp/h2_jetPt_zg"), mcpjet.pt(), observable.zg, eventWeight);
    registry.fill(HIST("groomed/mcp/h2_jetPt_rg"), mcpjet.pt(), observable.rg, eventWeight);
    hMCPZg[flavorIndex]->Fill(mcpjet.pt(), observable.zg, eventWeight);
    hMCPRg[flavorIndex]->Fill(mcpjet.pt(), observable.rg, eventWeight);
  }

  template <typename Tracks, typename Jet>
  bool isAcceptedJet(Jet const& jet)
  {
    static constexpr float largeNegativeNumber = -98.0f;
    static constexpr float largePositiveNumber = 9999.0f;

    const bool applyMin = leadingConstituentPtMin > largeNegativeNumber;
    const bool applyMax = leadingConstituentPtMax < largePositiveNumber;

    if (!applyMin && !applyMax) {
      return true;
    }

    bool passesMin = !applyMin;

    for (const auto& constituent : jet.template tracks_as<Tracks>()) {
      const float pt = constituent.pt();

      if (applyMin && pt >= leadingConstituentPtMin) {
        passesMin = true;
      }

      if (applyMax && pt > leadingConstituentPtMax) {
        return false;
      }
    }

    return passesMin;
  }

  void processDummy(FilteredCollisionsMCD::iterator const& /*collision*/)
  {
  }
  PROCESS_SWITCH(BjetSubstructureGnn, processDummy, "Dummy process function turned on by default", true);

  void processMCPFragmentation(AnalysisCollisionsMCP const& mcCollisions, FilteredMCPFragmentationJets const& MCPjets, aod::JetParticles const& /*particles*/)
  {
    for (const auto& mcCollision : mcCollisions) {
      if (mcCollision.isOutlier() || std::abs(mcCollision.posZ()) > vertexZCut) {
        continue;
      }

      const float eventWeight = useEventWeight ? mcCollision.weight() : 1.f;

      registry.fill(HIST("fragmentation/mcp/h_event_counter"), 0., eventWeight);
      registry.fill(HIST("fragmentation/mcp/h_vertexZ"), mcCollision.posZ(), eventWeight);

      const auto& jetsInCollision = MCPjets.sliceBy(mcpFragmentationJetsPerMCPCollision, mcCollision.mcCollisionId());

      for (const auto& mcpjet : jetsInCollision) {
        int8_t jetFlavor = mcpjet.origin();
        const auto flavorCat = getJetFlavorCat(jetFlavor);

        if (mcpjet.pt() <= 0.0f) {
          continue;
        }

        fillMCPFragmentationJetHistograms(mcpjet, flavorCat, eventWeight);

        const auto constituents = mcpjet.template tracks_as<aod::JetParticles>();
        registry.fill(HIST("fragmentation/mcp/h2_jetPt_nConstituents"), mcpjet.pt(), static_cast<double>(constituents.size()), eventWeight);

        for (const auto& constituent : constituents) {
          const auto observable = calculateFragmentation(mcpjet, constituent);

          fillMCPFragmentationHistograms(mcpjet, observable, flavorCat, eventWeight);
        }
      }
    }
  }
  PROCESS_SWITCH(BjetSubstructureGnn, processMCPFragmentation, "Fill particle-level MC jet fragmentation histograms", false);

  void processMCDFragmentation(FilteredCollisionsMCD::iterator const& collision, FilteredMCDFragmentationJets const& MCDjets, FilteredTracksMCD const& /*tracks*/)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }
    if (collision.isOutlier()) {
      return;
    }

    float eventWeight = useEventWeight ? collision.weight() : 1.f;

    registry.fill(HIST("fragmentation/mcd/h_event_counter"), 0.0, eventWeight);
    registry.fill(HIST("fragmentation/mcd/h_vertexZ"), collision.posZ(), eventWeight);

    for (const auto& mcdjet : MCDjets) {
      int8_t jetFlavor = mcdjet.origin();
      const auto flavorCat = getJetFlavorCat(jetFlavor);

      if (mcdjet.pt() <= 0.0f) {
        continue;
      }

      if (!isAcceptedJet<FilteredTracksMCD>(mcdjet)) {
        continue;
      }

      const float Db = mcdjet.scoreML();

      fillMCDFragmentationJetHistograms(mcdjet, flavorCat, Db, eventWeight);

      const auto constituents = mcdjet.template tracks_as<FilteredTracksMCD>();
      registry.fill(HIST("fragmentation/mcd/h2_jetPt_nConstituents"), mcdjet.pt(), static_cast<double>(constituents.size()), eventWeight);

      for (const auto& constituent : constituents) {
        const auto observable = calculateFragmentation(mcdjet, constituent);

        fillMCDFragmentationHistograms(mcdjet, observable, flavorCat, Db, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BjetSubstructureGnn, processMCDFragmentation, "Fill detector-level MC jet fragmentation histograms", false);

  void processMCPGroomedSubstructure(AnalysisCollisionsMCP const& mcCollisions, FilteredMCPGroomedJets const& MCPJets)
  {
    for (const auto& mcCollision : mcCollisions) {
      if (mcCollision.isOutlier() || std::abs(mcCollision.posZ()) > vertexZCut) {
        continue;
      }

      const float eventWeight = useEventWeight ? mcCollision.weight() : 1.f;

      registry.fill(HIST("groomed/mcp/h_event_counter"), 0., eventWeight);
      registry.fill(HIST("groomed/mcp/h_vertexZ"), mcCollision.posZ(), eventWeight);

      const auto jetsInCollision = MCPJets.sliceBy(mcpGroomedJetsPerMCPCollision, mcCollision.mcCollisionId());

      for (const auto& mcpjet : jetsInCollision) {
        if (mcpjet.pt() <= 0.f) {
          continue;
        }

        const int8_t jetFlavor = mcpjet.origin();
        const auto flavorCat = getJetFlavorCat(jetFlavor);

        fillMCPGroomedJetHistograms(mcpjet, flavorCat, eventWeight);

        const auto observable = calculateGroomedSubstructure(mcpjet);

        fillMCPGroomedSubstructureHistograms(mcpjet, observable, flavorCat, eventWeight);
      }
    }
  }
  PROCESS_SWITCH(BjetSubstructureGnn, processMCPGroomedSubstructure, "Fill particle-level MC groomed jet substructure histograms", false);

  void processMCDGroomedSubstructure(FilteredCollisionsMCD::iterator const& collision, FilteredMCDGroomedJets const& MCDJets, FilteredTracksMCD const& /*tracks*/)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelectionBits)) {
      return;
    }

    if (collision.isOutlier()) {
      return;
    }

    const float eventWeight = useEventWeight ? collision.weight() : 1.f;

    registry.fill(HIST("groomed/mcd/h_event_counter"), 0., eventWeight);
    registry.fill(HIST("groomed/mcd/h_vertexZ"), collision.posZ(), eventWeight);

    for (const auto& mcdjet : MCDJets) {
      if (mcdjet.pt() <= 0.f) {
        continue;
      }

      if (!isAcceptedJet<FilteredTracksMCD>(mcdjet)) {
        continue;
      }

      const int8_t jetFlavor = mcdjet.origin();
      const auto flavorCat = getJetFlavorCat(jetFlavor);
      const float Db = mcdjet.scoreML();

      fillMCDGroomedJetHistograms(mcdjet, flavorCat, Db, eventWeight);

      const auto observable = calculateGroomedSubstructure(mcdjet);

      fillMCDGroomedSubstructureHistograms(mcdjet, observable, flavorCat, Db, eventWeight);
    }
  }
  PROCESS_SWITCH(BjetSubstructureGnn, processMCDGroomedSubstructure, "Fill detector-level MC groomed jet substructure histograms", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<BjetSubstructureGnn>(cfgc)};
}
