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

#include <string>

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
#include "Common/DataModel/Qvectors.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#include "Framework/runDataProcessing.h"

struct JetSpectraEseTask {
  ConfigurableAxis binJetPt{"binJetPt", {200, 0., 200.}, ""};
  ConfigurableAxis bindPhi{"bindPhi", {100, -TMath::Pi() - 1, TMath::Pi() + 1}, ""};
  ConfigurableAxis binESE{"binESE", {100, 0, 100}, ""};
  ConfigurableAxis binCos{"binCos", {100, -1.05, 1.05}, ""};

  Configurable<float> jetPtMin{"jetPtMin", 5.0, "minimum jet pT cut"};
  Configurable<float> jetR{"jetR", 0.2, "jet resolution parameter"};
  Configurable<float> vertexZCut{"vertexZCut", 10.0, "vertex z cut"};
  Configurable<std::vector<float>> CentRange{"CentRange", {30, 50}, "centrality region of interest"};
  Configurable<double> leadingJetPtCut{"fLeadingJetPtCut", 5.0, "leading jet pT cut"};

  Configurable<std::string> eventSelections{"eventSelections", "sel8", "choose event selection"};
  Configurable<std::string> trackSelections{"trackSelections", "globalTracks", "set track selections"};

  Configurable<int> fColSwitch{"fColSwitch", 0, "collision switch"};

  AxisSpec jetPtAxis = {binJetPt, "#it{p}_{T,jet}"};
  AxisSpec dPhiAxis = {bindPhi, "#Delta#phi"};
  AxisSpec eseAxis = {binESE, "#it{q}_{2}"};

  AxisSpec cosAxis = {binCos, ""};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  int eventSelection = -1;
  int trackSelection = -1;

  void init(o2::framework::InitContext&)
  {
    eventSelection = jetderiveddatautilities::initialiseEventSelection(static_cast<std::string>(eventSelections));
    trackSelection = jetderiveddatautilities::initialiseTrackSelection(static_cast<std::string>(trackSelections));

    LOGF(info, "jetSpectraEse::init()");

    switch (fColSwitch) {
      case 0:
        LOGF(info, "JetSpectraEseTask::init() - using data");
        registry.add("hEventCounter", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("hJetPt", "jet pT;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetPt_bkgsub", "jet pT background sub;#it{p}_{T,jet} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("hJetEta", "jet #eta;#eta_{jet};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("hJetPhi", "jet #phi;#phi_{jet};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("hRho", ";#rho;entries", {HistType::kTH1F, {{100, 0, 200.}}});
        registry.add("hJetArea", ";area_{jet};entries", {HistType::kTH1F, {{100, 0, 10.}}});
        registry.add("hdPhi", "#Delta#phi;entries;", {HistType::kTH1F, {{dPhiAxis}}});
        registry.add("hJetPtdPhiq2", "", {HistType::kTH3F, {{jetPtAxis}, {dPhiAxis}, {eseAxis}}});
        registry.add("hPsi2FT0C", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2FT0A", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2FV0A", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2TPCpos", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hPsi2TPCneg", ";Centrality; #Psi_{2}", {HistType::kTH2F, {{100, 0, 100}, {150, -2.5, 2.5}}});
        registry.add("hCosPsi2FT0CmFT0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FT0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0CmFV0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FV0AmFT0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmFT0C", ";Centrality;cos(2(#Psi_{2}^{FT0A}-#Psi_{2}^{FT0C}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FT0AmFV0A", ";Centrality;cos(2(#Psi_{2}^{FT0C}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2FV0AmFT0C", ";Centrality;cos(2(#Psi_{2}^{FV0A}-#Psi_{2}^{FT0C}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCposmTPCneg", ";Centrality;cos(2(#Psi_{2}^{TPCpos}-#Psi_{2}^{TPCneg}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCposmFV0A", ";Centrality;cos(2(#Psi_{2}^{TPCpos}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});
        registry.add("hCosPsi2TPCnegmFV0A", ";Centrality;cos(2(#Psi_{2}^{TPCneg}-#Psi_{2}^{FV0A}));#it{q}_{2}", {HistType::kTH3F, {{100, 0, 100}, {cosAxis}, {eseAxis}}});

        break;
      case 1:
        LOGF(info, "JetSpectraEseTask::init() - using MC");
        registry.add("h_mc_collisions", "event status;event status;entries", {HistType::kTH1F, {{10, 0.0, 10.0}}});
        registry.add("h_part_jet_pt", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_part_jet_eta", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_part_jet_phi", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_part_jet_pt_match", "particle level jet pT;#it{p}_{T,jet part} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_part_jet_eta_match", "particle level jet #eta;#eta_{jet part};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_part_jet_phi_match", "particle level jet #phi;#phi_{jet part};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_detector_jet_pt", "detector level jet pT;#it{p}_{T,jet det} (GeV/#it{c});entries", {HistType::kTH1F, {{jetPtAxis}}});
        registry.add("h_detector_jet_eta", "detector level jet #eta;#eta_{jet det};entries", {HistType::kTH1F, {{100, -1.0, 1.0}}});
        registry.add("h_detector_jet_phi", "detector level jet #phi;#phi_{jet det};entries", {HistType::kTH1F, {{80, -1.0, 7.}}});
        registry.add("h_matched_jets_pt_delta", "#it{p}_{T,jet part}; det - part", {HistType::kTH2F, {{jetPtAxis}, {200, -20., 20.0}}});
        registry.add("h_matched_jets_eta_delta", "#eta_{jet part}; det - part", {HistType::kTH2F, {{100, -1.0, 1.0}, {200, -20.0, 20.0}}});
        registry.add("h_matched_jets_phi_delta", "#phi_{jet part}; det - part", {HistType::kTH2F, {{80, -1.0, 7.}, {200, -20.0, 20.}}});
        registry.add("h_response_mat_match", "#it{p}_{T, jet det}; #it{p}_{T, jet part}", HistType::kTH2F, {jetPtAxis, jetPtAxis});
        break;
    }
  }

  Filter jetCuts = aod::jet::pt > jetPtMin&& aod::jet::r == nround(jetR.node() * 100.0f) && nabs(aod::jet::eta) < 0.9f - jetR;
  Filter colFilter = nabs(aod::jcollision::posZ) < vertexZCut;
  Filter mcCollisionFilter = nabs(aod::jmccollision::posZ) < vertexZCut;

  void processESEDataCharged(soa::Join<aod::JetCollisions, aod::JCollisionPIs, aod::BkgChargedRhos>::iterator const& collision,
                             soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFV0AVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QPercentileFT0Cs> const&,
                             soa::Filtered<aod::ChargedJets> const& jets,
                             aod::JetTracks const& tracks)
  {
    float counter{0.5f};
    registry.fill(HIST("hEventCounter"), counter++);
    const auto originalCollision = collision.collision_as<soa::Join<aod::Collisions, aod::CentFT0Cs, aod::QvectorFT0CVecs, aod::QvectorFT0AVecs, aod::QvectorFV0AVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QPercentileFT0Cs>>();
    registry.fill(HIST("hEventCounter"), counter++);
    if (originalCollision.centFT0C() < CentRange->at(0) || originalCollision.centFT0C() > CentRange->at(1))
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    const auto vPsi2 = procEP(originalCollision);
    const auto qPerc = originalCollision.qPERCFT0C();
    if (qPerc[0] < 0)
      return;
    registry.fill(HIST("hEventCounter"), counter++);

    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    registry.fill(HIST("hEventCounter"), counter++);

    if (!isAcceptedLeadTrack(tracks))
      return;

    registry.fill(HIST("hEventCounter"), counter++);
    registry.fill(HIST("hRho"), collision.rho());
    for (auto const& jet : jets) {
      float jetpT_bkgsub = jet.pt() - (collision.rho() * jet.area());
      registry.fill(HIST("hJetPt"), jet.pt());
      registry.fill(HIST("hJetPt_bkgsub"), jetpT_bkgsub);
      registry.fill(HIST("hJetEta"), jet.eta());
      registry.fill(HIST("hJetPhi"), jet.phi());
      registry.fill(HIST("hJetArea"), jet.area());

      float dPhi = RecoDecay::constrainAngle(jet.phi() - vPsi2, -o2::constants::math::PI);
      registry.fill(HIST("hdPhi"), dPhi);
      registry.fill(HIST("hJetPtdPhiq2"), jetpT_bkgsub, dPhi, qPerc[0]); /* check the dphi */
    }
    registry.fill(HIST("hEventCounter"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processESEDataCharged, "process ese collisions", true);

  void processMCParticleLevel(soa::Filtered<aod::ChargedMCParticleLevelJets>::iterator const& jet)
  {
    registry.fill(HIST("h_part_jet_pt"), jet.pt());
    registry.fill(HIST("h_part_jet_eta"), jet.eta());
    registry.fill(HIST("h_part_jet_phi"), jet.phi());
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCParticleLevel, "jets on particle level MC", false);

  using JetMCPTable = soa::Filtered<soa::Join<aod::ChargedMCParticleLevelJets, aod::ChargedMCParticleLevelJetConstituents, aod::ChargedMCParticleLevelJetsMatchedToChargedMCDetectorLevelJets>>;
  void processMCChargedMatched(soa::Filtered<aod::JetCollisionsMCD>::iterator const& collision,
                               soa::Filtered<soa::Join<aod::ChargedMCDetectorLevelJets, aod::ChargedMCDetectorLevelJetConstituents, aod::ChargedMCDetectorLevelJetsMatchedToChargedMCParticleLevelJets>> const& mcdjets,
                               JetMCPTable const&,
                               aod::JetTracks const&,
                               aod::JetParticles const&)
  {
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return;

    float counter{0.5f};
    registry.fill(HIST("h_mc_collisions"), counter++);
    for (const auto& mcdjet : mcdjets) {

      registry.fill(HIST("h_detector_jet_pt"), mcdjet.pt());
      registry.fill(HIST("h_detector_jet_eta"), mcdjet.eta());
      registry.fill(HIST("h_detector_jet_phi"), mcdjet.phi());
      for (auto& mcpjet : mcdjet.template matchedJetGeo_as<JetMCPTable>()) {

        registry.fill(HIST("h_part_jet_pt_match"), mcpjet.pt());
        registry.fill(HIST("h_part_jet_eta_match"), mcpjet.eta());
        registry.fill(HIST("h_part_jet_phi_match"), mcpjet.phi());

        registry.fill(HIST("h_matched_jets_pt_delta"), mcpjet.pt(), mcdjet.pt() - mcpjet.pt());
        registry.fill(HIST("h_matched_jets_phi_delta"), mcpjet.phi(), mcdjet.phi() - mcpjet.phi());
        registry.fill(HIST("h_matched_jets_eta_delta"), mcpjet.eta(), mcdjet.eta() - mcpjet.eta());

        registry.fill(HIST("h_response_mat_match"), mcdjet.pt(), mcpjet.pt());
      }
    }
    registry.fill(HIST("h_mc_collisions"), counter++);
  }
  PROCESS_SWITCH(JetSpectraEseTask, processMCChargedMatched, "jet matched mcp and mcd", false);

  template <typename T>
  bool isAcceptedLeadTrack(T const& tracks)
  {
    double leadingTrackpT{0.0};
    for (const auto& track : tracks) {
      if (track.pt() > leadingJetPtCut) {
        if (track.pt() > leadingTrackpT) {
          leadingTrackpT = track.pt();
        }
      }
    }
    if (leadingTrackpT == 0.0)
      return false;
    else
      return true;
  }

  template <typename qVectors>
  float procEP(qVectors const& vec)
  {
    const auto epFT0A = 1 / 2.0 * TMath::ATan2(vec.qvecFT0AImVec()[0], vec.qvecFT0AReVec()[0]);
    const auto epFV0A = 1 / 2.0 * std::atan2(vec.qvecFV0AImVec()[0], vec.qvecFV0AReVec()[0]);
    const auto epFT0C = 1 / 2.0 * std::atan2(vec.qvecFT0CImVec()[0], vec.qvecFT0CReVec()[0]);
    const auto epTPCpos = 1 / 2.0 * std::atan2(vec.qvecTPCposImVec()[0], vec.qvecTPCposReVec()[0]);
    const auto epTPCneg = 1 / 2.0 * std::atan2(vec.qvecTPCnegImVec()[0], vec.qvecTPCnegReVec()[0]);

    registry.fill(HIST("hPsi2FT0C"), vec.centFT0C(), epFT0C);
    registry.fill(HIST("hPsi2FT0A"), vec.centFT0C(), epFT0A);
    registry.fill(HIST("hPsi2FV0A"), vec.centFT0C(), epFV0A);
    registry.fill(HIST("hPsi2TPCpos"), vec.centFT0C(), epTPCpos);
    registry.fill(HIST("hPsi2TPCneg"), vec.centFT0C(), epTPCneg);

    const auto cosPsi2FT0CmFT0A = cosPsiXY(epFT0C, epFT0A);
    const auto cosPsi2FT0CmFV0A = cosPsiXY(epFT0C, epFV0A);
    const auto cosPsi2FV0AmFT0A = cosPsiXY(epFV0A, epFT0A);
    const auto cosPsi2FT0AmFT0C = cosPsiXY(epFT0A, epFT0C);
    const auto cosPsi2FT0AmFV0A = cosPsiXY(epFT0A, epFV0A);
    const auto cosPsi2FV0AmFT0C = cosPsiXY(epFV0A, epFT0C);
    const auto cosPsi2TPCposmTPCneg = cosPsiXY(epTPCpos, epTPCneg);
    const auto cosPsi2TPCposmFV0A = cosPsiXY(epTPCpos, epFV0A);
    const auto cosPsi2TPCnegmFV0A = cosPsiXY(epTPCneg, epFV0A);

    registry.fill(HIST("hCosPsi2FT0CmFT0A"), vec.centFT0C(), cosPsi2FT0CmFT0A, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2FT0CmFV0A"), vec.centFT0C(), cosPsi2FT0CmFV0A, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2FV0AmFT0A"), vec.centFT0C(), cosPsi2FV0AmFT0A, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2FT0AmFT0C"), vec.centFT0C(), cosPsi2FT0AmFT0C, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2FT0AmFV0A"), vec.centFT0C(), cosPsi2FT0AmFV0A, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2FV0AmFT0C"), vec.centFT0C(), cosPsi2FV0AmFT0C, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2TPCposmTPCneg"), vec.centFT0C(), cosPsi2TPCposmTPCneg, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2TPCposmFV0A"), vec.centFT0C(), cosPsi2TPCposmFV0A, vec.qPERCFT0C()[0]);
    registry.fill(HIST("hCosPsi2TPCnegmFV0A"), vec.centFT0C(), cosPsi2TPCnegmFV0A, vec.qPERCFT0C()[0]);

    return epFT0A;
  }

  template <typename Psi>
  float cosPsiXY(Psi const& psiX, Psi const& psiY)
  {
    return std::cos(2.0 * (psiX - psiY));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<JetSpectraEseTask>(cfgc)}; }
