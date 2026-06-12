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

/// \file jetSpectraChargedGen.cxx
/// \brief Charged-particle jet spectra at MC-particle level for generator-only
///        AODs (no reconstruction tables required). Companion to
///        jetSpectraCharged.cxx for productions such as EPOS 4 on-the-fly
///        where only mc-particle and mc-collision tables exist.
/// \author Joonsuk Bae <joonsuk.bae@cern.ch>

#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/Core/JetUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetSubtraction.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct JetSpectraChargedGen {

  using JetBkgRhoMcCollisions = soa::Join<aod::JetMcCollisions, aod::BkgChargedMcRhos>;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<float> selectedJetsRadius{"selectedJetsRadius", 0.2, "resolution parameter for histograms without radius"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.7, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.7, "maximum jet pseudorapidity"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum particle pseudorapidity"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum particle pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "minimum area fraction relative to pi R^2"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT of the leading constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT of the leading constituent"};
  Configurable<float> pTHatMaxMCP{"pTHatMaxMCP", 999.0, "maximum jet pT in units of pTHat (outlier rejection)"};
  Configurable<float> pTHatAbsoluteMin{"pTHatAbsoluteMin", -99.0, "minimum pTHat"};
  Configurable<float> pTHatExponent{"pTHatExponent", 6.0, "exponent of the event weight for the back-calculation of pTHat"};
  Configurable<float> simPtRef{"simPtRef", 10.0, "reference pT for the back-calculation of pTHat from the event weight"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "truth vertex-z cut; <=0 disables"};
  Configurable<float> kappa{"kappa", 1.0, "angularity kappa"};
  Configurable<float> alpha{"alpha", 1.0, "angularity alpha"};
  Configurable<double> jetPtMax{"jetPtMax", 200.0, "max of the jet pT axis"};

  static constexpr float kConfigSwitchLow = -98.0f;
  static constexpr float kConfigSwitchHigh = 9998.0f;

  void init(InitContext const&)
  {
    AxisSpec jetPtAxis = {static_cast<int>(jetPtMax), 0.0, jetPtMax, "#it{p}_{T,jet} (GeV/#it{c})"};
    AxisSpec etaAxis = {100, -1.0, 1.0, "#eta_{jet}"};
    AxisSpec phiAxis = {160, -1.0, 7.0, "#varphi_{jet}"};
    AxisSpec rhoAxis = {200, 0.0, 200.0, "#rho (GeV/#it{c})"};

    registry.add("h_mccollisions", "MC collisions counter;;entries", {HistType::kTH1F, {{1, 0.0, 1.0}}});
    registry.add("h_mccollisions_rho", "#rho distribution;#rho (GeV/#it{c});entries", {HistType::kTH1F, {rhoAxis}});
    registry.add("p_xSection", ";;<#sigma_{gen}> from HepMC (pb)", {HistType::kTProfile, {{1, 0.0, 1.0}}});

    registry.add("h_jet_pt_part", ";#it{p}_{T,jet}^{part} (GeV/#it{c});entries", {HistType::kTH1F, {jetPtAxis}});
    registry.add("h_jet_eta_part", ";#eta_{jet}^{part};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_jet_phi_part", ";#varphi_{jet}^{part};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h3_jet_pt_jet_eta_jet_phi_part", ";#it{p}_{T,jet}^{part};#eta_{jet}^{part};#varphi_{jet}^{part}",
                 {HistType::kTH3F, {jetPtAxis, etaAxis, phiAxis}});
    registry.add("h2_jet_pt_part_jet_area_part", ";#it{p}_{T,jet}^{part};jet area",
                 {HistType::kTH2F, {jetPtAxis, {200, 0.0, 4.0}}});
    registry.add("h2_jet_pt_part_jet_ntracks_part", ";#it{p}_{T,jet}^{part};#it{N}_{tracks}",
                 {HistType::kTH2F, {jetPtAxis, {100, 0.0, 100.0}}});
    registry.add("h2_jet_pt_part_track_pt_part", ";#it{p}_{T,jet}^{part};#it{p}_{T,track}^{part}",
                 {HistType::kTH2F, {jetPtAxis, jetPtAxis}});
    registry.add("h2_jet_pt_jet_angularity_part", ";#it{p}_{T,jet}^{part};angularity",
                 {HistType::kTH2F, {jetPtAxis, {200, 0.0, 4.0}}});

    registry.add("h_jet_pt_part_rhoareasubtracted", ";#it{p}_{T,jet}^{part,sub} (GeV/#it{c});entries",
                 {HistType::kTH1F, {{static_cast<int>(2 * jetPtMax), -jetPtMax, jetPtMax}}});
    registry.add("h_jet_eta_part_rhoareasubtracted", ";#eta_{jet}^{part};entries", {HistType::kTH1F, {etaAxis}});
    registry.add("h_jet_phi_part_rhoareasubtracted", ";#varphi_{jet}^{part};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("h3_jet_pt_jet_eta_jet_phi_part_rhoareasubtracted",
                 ";#it{p}_{T,jet}^{part,sub};#eta_{jet}^{part};#varphi_{jet}^{part}",
                 {HistType::kTH3F, {{static_cast<int>(2 * jetPtMax), -jetPtMax, jetPtMax}, etaAxis, phiAxis}});
    registry.add("h2_jet_pt_part_jet_area_part_rhoareasubtracted", ";#it{p}_{T,jet}^{part,sub};jet area",
                 {HistType::kTH2F, {{static_cast<int>(2 * jetPtMax), -jetPtMax, jetPtMax}, {200, 0.0, 4.0}}});
    registry.add("h2_jet_pt_part_jet_ntracks_part_rhoareasubtracted", ";#it{p}_{T,jet}^{part,sub};#it{N}_{tracks}",
                 {HistType::kTH2F, {{static_cast<int>(2 * jetPtMax), -jetPtMax, jetPtMax}, {100, 0.0, 100.0}}});
  }

  // Truth-only zvtx cut: returns true if accepted.
  template <typename TMcColl>
  bool passVertexZ(TMcColl const& mccollision)
  {
    if (vertexZCut <= 0.0f) {
      return true;
    }
    return std::abs(mccollision.posZ()) < vertexZCut;
  }

  template <typename TJet>
  bool isAcceptedJet(TJet const& jet)
  {
    if (jetAreaFractionMin > kConfigSwitchLow) {
      if (jet.area() < jetAreaFractionMin * o2::constants::math::PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkMin = (leadingConstituentPtMin > kConfigSwitchLow);
    bool checkMax = (leadingConstituentPtMax < kConfigSwitchHigh);
    if (!checkMin && !checkMax) {
      return true;
    }
    bool passMin = !checkMin;
    bool passMax = checkMax;
    for (const auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
      double pt = constituent.pt();
      if (checkMin && pt >= leadingConstituentPtMin) {
        passMin = true;
      }
      if (checkMax && pt > leadingConstituentPtMax) {
        passMax = false;
      }
    }
    return passMin && passMax;
  }

  template <typename TJet>
  void fillRawHistograms(TJet const& jet, float weight, float pTHat)
  {
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    if (jet.r() == std::round(selectedJetsRadius * 100.0f)) {
      registry.fill(HIST("h_jet_pt_part"), jet.pt(), weight);
      registry.fill(HIST("h_jet_eta_part"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part"), jet.phi(), weight);
      registry.fill(HIST("h3_jet_pt_jet_eta_jet_phi_part"), jet.pt(), jet.eta(), jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_area_part"), jet.pt(), jet.area(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_ntracks_part"), jet.pt(), jet.tracksIds().size(), weight);
    }
    float angularity = 0.0f;
    for (const auto& constituent : jet.template tracks_as<aod::JetParticles>()) {
      registry.fill(HIST("h2_jet_pt_part_track_pt_part"), jet.pt(), constituent.pt(), weight);
      angularity += std::pow(constituent.pt(), kappa) * std::pow(jetutilities::deltaR(jet, constituent), alpha);
    }
    angularity /= (jet.pt() * (jet.r() / 100.0f));
    registry.fill(HIST("h2_jet_pt_jet_angularity_part"), jet.pt(), angularity, weight);
  }

  template <typename TJet>
  void fillAreaSubHistograms(TJet const& jet, float rho, float weight, float pTHat)
  {
    if (jet.pt() > pTHatMaxMCP * pTHat || pTHat < pTHatAbsoluteMin) {
      return;
    }
    if (jet.r() != std::round(selectedJetsRadius * 100.0f)) {
      return;
    }
    double ptSub = jet.pt() - (rho * jet.area());
    registry.fill(HIST("h_jet_pt_part_rhoareasubtracted"), ptSub, weight);
    registry.fill(HIST("h3_jet_pt_jet_eta_jet_phi_part_rhoareasubtracted"), ptSub, jet.eta(), jet.phi(), weight);
    if (ptSub > 0.0) {
      registry.fill(HIST("h_jet_eta_part_rhoareasubtracted"), jet.eta(), weight);
      registry.fill(HIST("h_jet_phi_part_rhoareasubtracted"), jet.phi(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_area_part_rhoareasubtracted"), ptSub, jet.area(), weight);
      registry.fill(HIST("h2_jet_pt_part_jet_ntracks_part_rhoareasubtracted"), ptSub, jet.tracksIds().size(), weight);
    }
  }

  void processDummy(aod::JDummys const&) {}
  PROCESS_SWITCH(JetSpectraChargedGen, processDummy, "dummy process", true);

  void processChargedGen(aod::JetMcCollisions::iterator const& mccollision,
                         soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                         aod::JetParticles const&)
  {
    if (!passVertexZ(mccollision)) {
      return;
    }
    registry.fill(HIST("h_mccollisions"), 0.5);
    registry.fill(HIST("p_xSection"), 0.5, mccollision.xsectGen());
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet(jet)) {
        continue;
      }
      fillRawHistograms(jet, /*weight*/ 1.0f, /*pTHat*/ 999.0f);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedGen, processChargedGen,
                 "charged-jet spectra at MC particle level (generator-only AOD)", false);

  void processChargedGenWeighted(aod::JetMcCollisions::iterator const& mccollision,
                                 soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                 aod::JetParticles const&)
  {
    if (!passVertexZ(mccollision)) {
      return;
    }
    float weight = mccollision.weight();
    float pTHat = simPtRef / std::pow(weight, 1.0f / pTHatExponent);
    registry.fill(HIST("h_mccollisions"), 0.5, weight);
    registry.fill(HIST("p_xSection"), 0.5, mccollision.xsectGen(), weight);
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet(jet)) {
        continue;
      }
      fillRawHistograms(jet, weight, pTHat);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedGen, processChargedGenWeighted,
                 "charged-jet spectra at MC particle level, weighted (generator-only AOD)", false);

  void processChargedAreaSubGen(JetBkgRhoMcCollisions::iterator const& mccollision,
                                soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                aod::JetParticles const&)
  {
    if (!passVertexZ(mccollision)) {
      return;
    }
    registry.fill(HIST("h_mccollisions"), 0.5);
    registry.fill(HIST("p_xSection"), 0.5, mccollision.xsectGen());
    registry.fill(HIST("h_mccollisions_rho"), mccollision.rho());
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet(jet)) {
        continue;
      }
      fillAreaSubHistograms(jet, mccollision.rho(), /*weight*/ 1.0f, /*pTHat*/ 999.0f);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedGen, processChargedAreaSubGen,
                 "charged-jet area-subtracted spectra at MC particle level (generator-only AOD)", false);

  void processChargedAreaSubGenWeighted(JetBkgRhoMcCollisions::iterator const& mccollision,
                                        soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents> const& jets,
                                        aod::JetParticles const&)
  {
    if (!passVertexZ(mccollision)) {
      return;
    }
    float weight = mccollision.weight();
    float pTHat = simPtRef / std::pow(weight, 1.0f / pTHatExponent);
    registry.fill(HIST("h_mccollisions"), 0.5, weight);
    registry.fill(HIST("p_xSection"), 0.5, mccollision.xsectGen(), weight);
    registry.fill(HIST("h_mccollisions_rho"), mccollision.rho(), weight);
    for (auto const& jet : jets) {
      if (!jetfindingutilities::isInEtaAcceptance(jet, jetEtaMin, jetEtaMax, trackEtaMin, trackEtaMax)) {
        continue;
      }
      if (!isAcceptedJet(jet)) {
        continue;
      }
      fillAreaSubHistograms(jet, mccollision.rho(), weight, pTHat);
    }
  }
  PROCESS_SWITCH(JetSpectraChargedGen, processChargedAreaSubGenWeighted,
                 "charged-jet area-subtracted spectra at MC particle level, weighted (generator-only AOD)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JetSpectraChargedGen>(cfgc)};
}
