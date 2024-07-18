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
///
/// \brief  task for an analysis of rho photoproduction in UPCs, intended usage is with UD tables
///         includes event tagging based on ZN information, track selection, reconstruction,
///         and also some basic stuff for phi anisotropy studies
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

// #include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h" // has some useful funtions for stuff not available from the tables

// ROOT headers
#include <Math/Vector4D.h> // these should apparently be used instead of TLorentzVector
#include <Math/Vector2D.h>
//#include "TEfficiency.h"

#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <chrono>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct upcMcRhoAnalysis {
  // configurables
  double PcEtaCut = 0.9; // cut on track eta as per physics coordination recommendation
  Configurable<bool> tracksRequireTOF{"tracksRequireTOF", false, "requireTOF"};
  Configurable<double> tracksTpcNSigmaPiCut{"tracksTpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
//   Configurable<double> tracksTofNSigmaPiCut{"tracksTofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
  Configurable<double> tracksPtMaxCut{"tracksPtMaxCut", 1.0, "ptMaxCut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "dcaMaxCut"};
  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "posZMaxCut"};
  Configurable<double> ZNcommonEnergyCut{"ZNcommonEnergyCut", 0.0, "ZNcommonEnergyCut"};
  Configurable<double> ZNtimeCut{"ZNtimeCut", 2.0, "ZNtimeCut"};
  Configurable<double> systemMassMinCut{"2systemMassMinCut", 0.5, "2systemMassMinCut"};
  Configurable<double> systemMassMaxCut{"2systemMassMaxCut", 1.1, "2systemMassMaxCut"};
  Configurable<double> systemPtCut{"2systemPtMaxCut", 0.12, "2systemPtMaxCut"};
  Configurable<double> systemYCut{"2systemYCut", 0.9, "2systemYCut"};

  ConfigurableAxis mAxis{"mAxis", {800, 0.0, 8.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {500, 0.0, 5.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis pt2Axis{"pt2Axis", {5000, 0.0, 0.05}, "#it{p}_{T}^{2} (GeV^{2}/#it{c}^{2})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, 2.0*o2::constants::math::PI}, "#phi"};
  ConfigurableAxis phiAssymAxis{"phiAssymAxis", {180, 0, o2::constants::math::PI}, "#phi"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&) {
    registry.add("QC/hMcParticlePdgCode", ";PDG code;counts", kTH1D, {{1001, -500.5, 500.5}});
    registry.add("QC/hMcParticleEta", ";#eta;counts", kTH1D, {etaAxis});
    registry.add("QC/hMcParticlePt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("sim/hRhoMass", ";m_{#rho} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("sim/hRhoPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("sim/hRhoY", ";y;counts", kTH1D, {yAxis});
    registry.add("sim/hRhoPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/hRhoMass", ";m_{#rho} (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/hRhoPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/hRhoY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/hRhoPhi", ";#phi;counts", kTH1D, {phiAxis});
  }

  template <typename T>
  bool trackPassesCuts(T const& track) { // check if track passes preliminary cuts (PID done separately)
    if (!track.isPVContributor()) return false;
    if (!track.hasITS()) return false;
    if (std::abs(track.dcaZ()) > tracksDcaMaxCut || std::abs(track.dcaXY()) > tracksDcaMaxCut) return false;
    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut) return false;
    if (std::abs(track.pt()) > tracksPtMaxCut) return false;
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) { // check if pre-cut tracks pass PID cut (fall within n-dimensional hypersphere)
    double radius = 0.0;
    for (const auto& track : cutTracks) radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  double tracksTotalCharge(const T& cutTracks) { // calculate total charge of selected tracks
    double charge = 0.0;
    for (const auto& track : cutTracks) charge += track.sign();
    return charge;
  }

  template <typename T>
  bool systemPassCuts(const T& system) { // check if system passes system cuts
    if (system.M() < systemMassMinCut || system.M() > systemMassMaxCut) return false;
    if (system.Pt() > systemPtCut) return false;
    if (std::abs(system.Rapidity()) > systemYCut) return false;
    return true;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) { // reconstruct system from 4-vectors
    ROOT::Math::PxPyPzMVector system;
    for (const auto& track4Vec : cutTracks4Vecs) system += track4Vec;
    return system;
  }

  template <typename T>
  double getPhiAssym(const T& cutTracks) { // phi anisotropy based on STAR + ALICE papers
    // assign the tracks randomly
    std::vector<int> indices = {0, 1};
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); // get time-based seed
    std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed)); // shuffle indices
    // calculate phi
    ROOT::Math::XYVector pOne(cutTracks[indices[0]].px(), cutTracks[indices[0]].py());
    ROOT::Math::XYVector pTwo(cutTracks[indices[1]].px(), cutTracks[indices[1]].py());
    auto pPlus = pOne + pTwo;
    auto pMinus = pOne - pTwo;
    // no method for direct calculation of angle -> use dot product formula
    double cosPhi = (pPlus.Dot(pMinus))/(std::sqrt(pPlus.Mag2()) * std::sqrt(pMinus.Mag2()));
    return std::acos(cosPhi);
  }

  void processSim(aod::McCollision const& mccollision, aod::McParticles const& mcparticles) {
    std::vector<ROOT::Math::PxPyPzMVector> mcPions;
    for (const auto& mcparticle : mcparticles) {
        registry.fill(HIST("QC/hMcParticlePdgCode"), mcparticle.pdgCode());
        if (!mcparticle.isPhysicalPrimary()) continue;
        if (!mcparticle.producedByGenerator()) continue; // I think this is the same as isPhysicalPrimary
        if (std::abs(mcparticle.pdgCode()) != 211) continue;
        registry.fill(HIST("QC/hMcParticleEta"), mcparticle.eta());
        registry.fill(HIST("QC/hMcParticlePt"), mcparticle.pt());
        mcPions.push_back(ROOT::Math::PxPyPzMVector(mcparticle.px(), mcparticle.py(), mcparticle.pz(), o2::constants::physics::MassPionCharged));
    }
    if (mcPions.size() != 2) return;
    auto rho = reconstructSystem(mcPions);
    if (std::abs(rho.Rapidity()) > 0.9) return;
    registry.fill(HIST("sim/hRhoMass"), rho.M());
    registry.fill(HIST("sim/hRhoPt"), rho.Pt());
    registry.fill(HIST("sim/hRhoY"), rho.Rapidity());
    registry.fill(HIST("sim/hRhoPhi"), rho.Phi() + o2::constants::math::PI);
  } PROCESS_SWITCH(upcMcRhoAnalysis, processSim, "analyse simulated tracks", true);

  void processReco(aod::Collision const& collision, soa::Join<aod::Tracks, aod::pidTPCPi, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection> const& tracks) {
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut) return;
    std::vector<decltype(tracks.begin())> cutTracks;
    for (const auto& track : tracks) {
      if (!trackPassesCuts(track)) continue;
      cutTracks.push_back(track);
    }
  } PROCESS_SWITCH(upcMcRhoAnalysis, processReco, "analyse reconstructed tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcMcRhoAnalysis>(cfgc)
  };
}