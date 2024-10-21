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

/// \file jetSpectraEseTask.cxx
/// \brief jet spectra analysis framework with ESE (19/08/2024)
///
/// \author Joachim C. K. B. Hansen, Lund University

#include "Framework/ASoA.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGJE/Core/FastJetUtilities.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"

#include "Common/DataModel/EseTable.h"
#include "Common/Core/FFitWeights.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

using ColWqVecFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QVecFT0Cs, aod::QPercentileFT0Cs>;
using JColwESE = soa::Join<aod::JCollisions, aod::CentFT0Cs, aod::QVecFT0Cs, aod::QPercentileFT0Cs>;

struct JetSpectraEseTask {
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis bindPhi{"bindPhi", {100, -TMath::Pi() - 1, TMath::Pi() + 1}, ""};
  ConfigurableAxis binESE{"binESE", {100, 0, 100}, ""};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.2, "jet resolution parameter"};
  Configurable<std::vector<float>> CentRange{"CentRange", {30, 50}, "centrality region of interest"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T,jet}"};
  AxisSpec dPhiAxis = {bindPhi, "#Delta#phi"};
  AxisSpec eseAxis = {binESE, "#it{q}_{2}"};

  HistogramRegistry registry{"registry",
                             {{"h_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}}},
                              {"h_jet_pt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_pt_bkgsub", "jet pT background sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{200, 0., 200.}}}},
                              {"h_jet_eta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}}},
                              {"h_jet_phi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}}},
                              {"h_rho", ";#rho;entries", {HistType::kTH1F, {{100, 0, 200.}}}},
                              {"h_jet_area", ";area_{jet};entries", {HistType::kTH1F, {{100, 0, 10.}}}},
                              {"h_Psi2", "#Psi_{2};entries;", {HistType::kTH1F, {{150, -2.5, 2.5}}}},
                              {"h_dPhi", "#Delta#phi;entries;", {HistType::kTH1F, {{dPhiAxis}}}},
                              {"jet_pt_dPhi_q2", "", {HistType::kTH3F, {{jetPtAxis}, {dPhiAxis}, {eseAxis}}}}}};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && nabs(aod::jet::eta) < 0.9f - jetR;

  using JColPI = soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::BkgChargedRhos>::iterator;
  using ColESE = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QVecFT0Cs, aod::QPercentileFT0Cs, aod::FEseCols>;
  using JJets = soa::Filtered<aod::ChargedJets>;
  void processESEDataCharged(JColPI const& collision,
                             ColESE const&,
                             JJets const& jets,
                             aod::JetTracks const& tracks)
  {
    auto originalCollision = collision.collision_as<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QVecFT0Cs, aod::QPercentileFT0Cs, aod::FEseCols>>();
    registry.fill(HIST("h_collisions"), 0.5);
    if (originalCollision.fESECOL()[0] != 1)
      return;
    registry.fill(HIST("h_collisions"), 1.5);
    if (originalCollision.centFT0C() < CentRange->at(0) || originalCollision.centFT0C() > CentRange->at(1))
      return;
    registry.fill(HIST("h_collisions"), 2.5);
    float vPsi2 = FFitWeights::EventPlane(originalCollision, 2);
    auto qPerc = originalCollision.qPERCFT0C();
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection)) {
      return;
    }
    registry.fill(HIST("h_collisions"), 3.5);
    double leadingTrackpT = 0.0;
    for (auto& track : tracks) {
      if (track.pt() > 5.0) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return;

    registry.fill(HIST("h_collisions"), 4.5);
    registry.fill(HIST("h_rho"), collision.rho());
    for (auto const& jet : jets) {
      float jetpT_bkgsub = jet.pt() - (collision.rho() * jet.area());
      registry.fill(HIST("h_jet_pt"), jet.pt());
      registry.fill(HIST("h_jet_pt_bkgsub"), jetpT_bkgsub);
      registry.fill(HIST("h_jet_eta"), jet.eta());
      registry.fill(HIST("h_jet_phi"), jet.phi());
      registry.fill(HIST("h_Psi2"), vPsi2);
      registry.fill(HIST("h_jet_area"), jet.area());

      float dPhi = RecoDecay::constrainAngle(jet.phi() - vPsi2, -o2::constants::math::PI);
      registry.fill(HIST("h_dPhi"), dPhi);
      if (qPerc[0] < 0)
        continue;
      registry.fill(HIST("jet_pt_dPhi_q2"), jetpT_bkgsub, dPhi, qPerc[0]);
    }
    registry.fill(HIST("h_collisions"), 5.5);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process ese collisions", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc)}; }
