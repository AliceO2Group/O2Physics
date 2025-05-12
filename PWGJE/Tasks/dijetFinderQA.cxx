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

// dijet finder QA task
//
/// \author Dongguk Kim <dongguk.kim@cern.ch>

#include <cmath>
#include <TRandom3.h>
#include <string>
#include <vector>

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetFinder.h"
#include "PWGJE/Core/JetFindingUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "PWGJE/Core/JetDerivedDataUtilities.h"

#include "EventFiltering/filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DijetFinderQATask {

  HistogramRegistry registry;

  Configurable<float> centralityMin{"centralityMin", -999.0, "minimum centrality"};
  Configurable<float> centralityMax{"centralityMax", 999.0, "maximum centrality"};
  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0f, "Accepted z-vertex range"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};
  Configurable<float> trackPtMin{"trackPtMin", 0.15, "minimum pT acceptance for tracks"};
  Configurable<float> trackPtMax{"trackPtMax", 1000.0, "maximum pT acceptance for tracks"};
  Configurable<float> trackEtaMin{"trackEtaMin", -0.9, "minimum eta acceptance for tracks"};
  Configurable<float> trackEtaMax{"trackEtaMax", 0.9, "maximum eta acceptance for tracks"};
  Configurable<float> leadingConstituentPtMin{"leadingConstituentPtMin", -99.0, "minimum pT selection on jet constituent"};
  Configurable<float> leadingConstituentPtMax{"leadingConstituentPtMax", 9999.0, "maximum pT selection on jet constituent"};
  Configurable<float> setJetPtCut{"setJetPtCut", 20., "set jet pt minimum cut"};
  Configurable<float> setPhiCut{"setPhiCut", 0.5, "set phicut"};
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> jetPtMin{"jetPtMin", 20.0, "minimum jet pT cut"};
  Configurable<double> jetPtMax{"jetPtMax", 200., "set jet pT bin max"};
  Configurable<float> jetEtaMin{"jetEtaMin", -0.5, "minimum jet pseudorapidity"};
  Configurable<float> jetEtaMax{"jetEtaMax", 0.5, "maximum jet pseudorapidity"};
  Configurable<float> jetAreaFractionMin{"jetAreaFractionMin", -99.0, "used to make a cut on the jet areas"};

  std::vector<int> eventSelection;
  int trackSelection = -1;

  std::vector<double> dijetMassBins;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
    // Add histogram for event counts

    auto dijetMassTemp = 0.0;
    while (dijetMassTemp <= 2 * jetPtMax) {
      dijetMassBins.push_back(dijetMassTemp);
      dijetMassTemp += 5.0;
    }

    AxisSpec dijetMassAxis = {dijetMassBins, "M_{jj} (GeV/#it{c}^2)"};
    AxisSpec eventCountAxis = {{0.5, 1.5}, "events"};

    if (doprocessDijetMCP) {
      registry.add("h_part_dijet_mass", "Dijet invariant mass;;entries", {HistType::kTH1F, {dijetMassAxis}});
      registry.add("hColCounterFinal_MCP", "Event count;;entries", {HistType::kTH1F, {eventCountAxis}});
    }

    if (doprocessDijetMCD) {
      registry.add("h_detec_dijet_mass", "Dijet invariant mass;;entries", {HistType::kTH1F, {dijetMassAxis}});
      registry.add("hColCounterFinal_MCD", "Event count;;entries", {HistType::kTH1F, {eventCountAxis}});
    }

    if (doprocessDijetData) {
      registry.add("h_data_dijet_mass", "Dijet invariant mass;;entries", {HistType::kTH1F, {dijetMassAxis}});
      registry.add("hColCounterFinal_Data", "Event count;;entries", {HistType::kTH1F, {eventCountAxis}});
    }

    if (doprocessDijetMCMatched) {
      registry.add("h_matched_dijet_mass", "M_{jj matched};M_{jj part}; M_{jj det}", {HistType::kTH2F, {dijetMassAxis, dijetMassAxis}});
    }
  }

  /****************************************************************************************************************************************************************/
  Filter trackCuts = (aod::jtrack::pt >= trackPtMin && aod::jtrack::pt < trackPtMax && aod::jtrack::eta > trackEtaMin && aod::jtrack::eta < trackEtaMax);
  Filter eventCuts = (nabs(aod::jcollision::posZ) < vertexZCut && aod::jcollision::centrality >= centralityMin && aod::jcollision::centrality < centralityMax);
  Filter mcCollisionsFilter = nabs(aod::jmccollision::posZ) < vertexZCut;
  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f);
  /****************************************************************************************************************************************************************/

  template <typename T, typename U>
  bool isAcceptedJet(U const& jet)
  {

    if (jetAreaFractionMin > -98.0) {
      if (jet.area() < jetAreaFractionMin * M_PI * (jet.r() / 100.0) * (jet.r() / 100.0)) {
        return false;
      }
    }
    bool checkConstituentPt = true;
    bool checkConstituentMinPt = (leadingConstituentPtMin > -98.0);
    bool checkConstituentMaxPt = (leadingConstituentPtMax < 9998.0);
    if (!checkConstituentMinPt && !checkConstituentMaxPt) {
      checkConstituentPt = false;
    }

    if (checkConstituentPt) {
      bool isMinLeadingConstituent = !checkConstituentMinPt;
      bool isMaxLeadingConstituent = true;

      for (const auto& constituent : jet.template tracks_as<T>()) {
        double pt = constituent.pt();

        if (checkConstituentMinPt && pt >= leadingConstituentPtMin) {
          isMinLeadingConstituent = true;
        }
        if (checkConstituentMaxPt && pt > leadingConstituentPtMax) {
          isMaxLeadingConstituent = false;
        }
      }
      return isMinLeadingConstituent && isMaxLeadingConstituent;
    }

    return true;
  }

  template <typename T>
  void fillMassHistogramsMCP(T const& mass)
  {
    registry.fill(HIST("h_part_dijet_mass"), mass);
  }

  template <typename T>
  void fillMassHistogramsMCD(T const& mass)
  {
    registry.fill(HIST("h_detec_dijet_mass"), mass);
  }

  template <typename T>
  void fillMassHistogramsData(T const& mass)
  {
    registry.fill(HIST("h_data_dijet_mass"), mass);
  }

  template <typename T>
  void fillMassHistogramsMCMatched(T const& mass_P, T const& mass_D)
  {
    registry.fill(HIST("h_matched_dijet_mass"), mass_P, mass_D);
  }

  void processDummy(aod::JDummys const&)
  {
  }
  PROCESS_SWITCH(DijetFinderQATask, processDummy, "dummy", false);

  void processDijetMCP(soa::Filtered<aod::JetMcCollisions>::iterator const&, soa::Filtered<aod::ChargedMCParticleLevelJets> const& jets, soa::SmallGroups<aod::JetCollisionsMCD> const& collisions)
  {
    if (collisions.size() == 0) {
      return;
    }
    for (auto& collision : collisions) {
      if (fabs(collision.posZ()) > vertexZCut || !jetderiveddatautilities::selectCollision(collision, eventSelection))
        return;
    }

    registry.fill(HIST("hColCounterFinal_MCP"), 1);

    std::vector<std::array<double, 3>> jetPtcuts;
    for (auto& jet : jets) {
      jetPtcuts.push_back({jet.pt(), jet.eta(), jet.phi()});
    }

    if (jetPtcuts.size() >= 2) {
      auto& leading_jet = jetPtcuts[0];
      bool found_pair = false;

      for (size_t i = 1; i < jetPtcuts.size() && !found_pair; i++) {
        auto& candidate_jet = jetPtcuts[i];
        Double_t dphi = fabs(leading_jet[2] - candidate_jet[2]);
        Double_t deta = fabs(leading_jet[1] - candidate_jet[1]);
        Double_t condition = fabs(dphi - M_PI);

        if (condition < setPhiCut * M_PI) {
          Double_t pt1 = leading_jet[0];
          Double_t pt2 = candidate_jet[0];
          Double_t dijet_mass = sqrt(2 * pt1 * pt2 * (cosh(deta) - cos(dphi)));
          fillMassHistogramsMCP(dijet_mass);
          found_pair = true;
        }
      }
    }
  }
  PROCESS_SWITCH(DijetFinderQATask, processDijetMCP, "QA for invariant mass of dijet in particle level MC", false);

  void processDijetMCD(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::ChargedMCDetectorLevelJets> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    registry.fill(HIST("hColCounterFinal_MCD"), 1);

    std::vector<std::array<double, 3>> jetPtcuts;
    for (auto& jet : jets) {
      jetPtcuts.push_back({jet.pt(), jet.eta(), jet.phi()});
    }

    if (jetPtcuts.size() >= 2) {
      auto& leading_jet = jetPtcuts[0];
      bool found_pair = false;

      for (size_t i = 1; i < jetPtcuts.size() && !found_pair; i++) {
        auto& candidate_jet = jetPtcuts[i];
        Double_t dphi = fabs(leading_jet[2] - candidate_jet[2]);
        Double_t deta = fabs(leading_jet[1] - candidate_jet[1]);
        Double_t condition = fabs(dphi - M_PI);

        if (condition < setPhiCut * M_PI) {
          Double_t pt1 = leading_jet[0];
          Double_t pt2 = candidate_jet[0];
          Double_t dijet_mass = sqrt(2 * pt1 * pt2 * (cosh(deta) - cos(dphi)));
          fillMassHistogramsMCD(dijet_mass);
          found_pair = true;
        }
      }
    }
  }
  PROCESS_SWITCH(DijetFinderQATask, processDijetMCD, "QA for invariant mass of dijet in detector level MC", false);

  void processDijetData(soa::Filtered<aod::JetCollisions>::iterator const& collision, soa::Filtered<aod::ChargedJets> const& jets)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    // Fill event count histogram
    registry.fill(HIST("hColCounterFinal_Data"), 1);

    std::vector<std::array<double, 3>> jetPtcuts;
    for (auto& jet : jets) {
      jetPtcuts.push_back({jet.pt(), jet.eta(), jet.phi()});
    }

    if (jetPtcuts.size() >= 2) {
      auto& leading_jet = jetPtcuts[0];
      bool found_pair = false;

      for (size_t i = 1; i < jetPtcuts.size() && !found_pair; i++) {
        auto& candidate_jet = jetPtcuts[i];
        Double_t dphi = fabs(leading_jet[2] - candidate_jet[2]);
        Double_t deta = fabs(leading_jet[1] - candidate_jet[1]);
        Double_t condition = fabs(dphi - M_PI);

        if (condition < setPhiCut * M_PI) {
          Double_t pt1 = leading_jet[0];
          Double_t pt2 = candidate_jet[0];
          Double_t dijet_mass = sqrt(2 * pt1 * pt2 * (cosh(deta) - cos(dphi)));
          fillMassHistogramsData(dijet_mass);
          found_pair = true;
        }
      }
    }
  }
  PROCESS_SWITCH(DijetFinderQATask, processDijetData, "QA for invariant mass of dijet in data", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processDijetMCMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                             soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                             JetMCPTable const&, aod::JetTracks const&, aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }

    std::vector<std::array<double, 3>> jetPtcuts_D;
    std::vector<std::array<double, 3>> jetPtcuts_P;

    for (auto& jet : mcdjets) {
      if (jet.has_matchedJetGeo()) {
        for (auto& matchedJet : jet.template matchedJetPt_as<JetMCPTable>()) {
          if (matchedJet.pt() > setJetPtCut) {
            jetPtcuts_D.push_back({jet.pt(), jet.eta(), jet.phi()});
            jetPtcuts_P.push_back({matchedJet.pt(), matchedJet.eta(), matchedJet.phi()});
            break;
          }
        }
      }
    }

    if (jetPtcuts_D.size() >= 2 && jetPtcuts_P.size() >= 2) {
      auto& leading_jet_D = jetPtcuts_D[0];
      auto& leading_jet_P = jetPtcuts_P[0];
      bool found_pair = false;

      for (size_t i = 1; i < jetPtcuts_D.size() && !found_pair; i++) {
        auto& candidate_jet_D = jetPtcuts_D[i];
        auto& candidate_jet_P = jetPtcuts_P[i];

        Double_t dphi_D = fabs(leading_jet_D[2] - candidate_jet_D[2]);
        Double_t deta_D = fabs(leading_jet_D[1] - candidate_jet_D[1]);
        Double_t dphi_P = fabs(leading_jet_P[2] - candidate_jet_P[2]);
        Double_t deta_P = fabs(leading_jet_P[1] - candidate_jet_P[1]);
        Double_t condition = fabs(dphi_D - M_PI);

        if (condition < setPhiCut * M_PI) {
          double pt1_D = leading_jet_D[0];
          double pt2_D = candidate_jet_D[0];
          double dijet_mass_D = sqrt(2 * pt1_D * pt2_D * (cosh(deta_D) - cos(dphi_D)));

          double pt1_P = leading_jet_P[0];
          double pt2_P = candidate_jet_P[0];
          double dijet_mass_P = sqrt(2 * pt1_P * pt2_P * (cosh(deta_P) - cos(dphi_P)));

          fillMassHistogramsMCMatched(dijet_mass_P, dijet_mass_D);
          found_pair = true;
        }
      }
    }
  }
  PROCESS_SWITCH(DijetFinderQATask, processDijetMCMatched, "QA for invariant mass of dijet in mcmactched", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<DijetFinderQATask>(cfgc, TaskName{"dijet-finder-charged-qa"})};
}
