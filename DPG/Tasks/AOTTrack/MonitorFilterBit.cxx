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

/// \file   MonitorFilterBit.cxx
/// \author Andrea Rossi <andrea.rossi@cern.ch>
///
/// \brief Task performing basic checks on filter-bit selections.
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::track;

struct CheckFilterBit {

  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis binsEta{"binsEta", {30, -1.5, 1.5}, ""};
  Configurable<float> zVtxCut{"zVtxCut", 10., "Primary Vtx z cut"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0., 2 * M_PI}, "Phi binning"};

  HistogramRegistry histos;
  Int_t ncollisionCounter = 0;

  void init(InitContext const&)
  {

    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axisEta{binsEta, "#it{#eta}"};
    const AxisSpec axisPhi{binsPhi, "#it{#eta}"};
    histos.add("EventProp/histMCcollZ", "MC coll Z (cm); #it{z_{MCcoll}} (cm)", kTH1D, {{100, -20., 20.}});
    histos.add("EventProp/histDatacollZ", "MC coll Z (cm); #it{z_{MCcoll}} (cm)", kTH1D, {{100, -20., 20.}});
    histos.add("EventProp/histPtTrackNegCollID", "pt", kTH1D, {axisPt});

    histos.add("Tracks/Reco/histptAll", "pt", kTH1D, {axisPt});
    histos.add("Tracks/Reco/histpt3DAll", "All tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB0", "FB0 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB1", "FB1 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB2", "FB2 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB3", "FB3 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/Reco/histpt3DFB4", "FB4 tracks;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/MCgen/histMCgenpt", "pt", kTH1D, {axisPt});
    histos.add("Tracks/MCgen/histMCgen3dPhysPrimary", "MC Phys. Prim.;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedProdRad1to15cm", "MC Prod Rad_xy 1 to 15 cm;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedProdRad1mumto5mm", "MC Prod Rad_xy 1#mum to 5 mm ;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/MCgen/histMCgen3dChargedfromHFdecay", "MC Phys. Prim from HF decay ;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCPhysPrimCollMatch/histpt", "pt", kTH1D, {axisPt});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCPhysPrimCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1to15cmCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});

    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB0", "FB0;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB1", "FB1;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB2", "FB2;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB3", "FB3;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
    histos.add("Tracks/RecoMCfromHFdecayCollMatch/histptFB4", "FB4;#it{p}_{T} (GeV/#it{c});#it{#eta};#it{#varphi}", kTH3D, {axisPt, axisEta, axisPhi});
  }

  void processData(soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension> const& tracks)
  { // test whether the name can be changed using the PROCESS_SWITCH method (see qaEventTrack class)

    for (auto& track : tracks) {

      if (std::abs(track.eta()) < 0.9)
        histos.fill(HIST("Tracks/Reco/histptAll"), track.pt());
      histos.fill(HIST("Tracks/Reco/histpt3DAll"), track.pt(), track.eta(), track.phi());
      if (track.isGlobalTrack())
        histos.fill(HIST("Tracks/Reco/histpt3DFB0"), track.pt(), track.eta(), track.phi());
      if (track.trackCutFlagFb1())
        histos.fill(HIST("Tracks/Reco/histpt3DFB1"), track.pt(), track.eta(), track.phi());
      if (track.trackCutFlagFb2())
        histos.fill(HIST("Tracks/Reco/histpt3DFB2"), track.pt(), track.eta(), track.phi());
      if (track.trackCutFlagFb3())
        histos.fill(HIST("Tracks/Reco/histpt3DFB3"), track.pt(), track.eta(), track.phi());
      if (track.trackCutFlagFb4())
        histos.fill(HIST("Tracks/Reco/histpt3DFB4"), track.pt(), track.eta(), track.phi());
      if (track.trackCutFlagFb5())
        histos.fill(HIST("Tracks/Reco/histpt3DFB5"), track.pt(), track.eta(), track.phi());
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processData, "process data", true);

  void processRecoMC(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions)
  { // this will loop over data (PV) collisions

    histos.fill(HIST("EventProp/histDatacollZ"), collision.posZ());
    if (std::abs(collision.posZ()) > 10.)
      return;
    ncollisionCounter++;
    for (auto& track : tracks) {

      if (track.collisionId() < 0)
        histos.fill(HIST("EventProp/histPtTrackNegCollID"), track.pt());
      bool has_MCparticle = track.has_mcParticle();
      if (has_MCparticle) {
        /// the track is not fake
        auto mcparticle = track.mcParticle();
        // auto collReco = track.collision_as<CollisionTableMC>();
        // auto collMC = mcparticle.mcCollision();
        auto mcCollID_recoColl = collision.mcCollisionId();
        auto mcCollID_particle = mcparticle.mcCollisionId();
        bool indexMatchOK = (mcCollID_recoColl == mcCollID_particle);
        // double pvZdiff = collision.posZ() - collMC.posZ();
        if (indexMatchOK) {
          double prodRadius2 = mcparticle.vx() * mcparticle.vx() + mcparticle.vy() * mcparticle.vy();
          int partpdg = std::abs(mcparticle.pdgCode());
          if (partpdg == 211 || partpdg == 321 || partpdg == 2212 || partpdg == 11 || partpdg == 13) {
            if (mcparticle.isPhysicalPrimary()) {
              int isHF = RecoDecay::getCharmHadronOrigin(mcParticles, mcparticle, false);
              if (std::abs(track.eta()) < 0.9) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histpt"), mcparticle.pt()); // note: one needs to avoid double counting of tracks reco both in TPC and ITS but not matched
              }
              if (track.isGlobalTrack()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb1()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb2()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb3()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb4()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
              }
              if (track.trackCutFlagFb5()) {
                histos.fill(HIST("Tracks/RecoMCPhysPrimCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
              }
              if (isHF == RecoDecay::OriginType::Prompt || isHF == RecoDecay::OriginType::NonPrompt) {
                if (track.isGlobalTrack()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb1()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb2()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb3()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb4()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
                }
                if (track.trackCutFlagFb5()) {
                  histos.fill(HIST("Tracks/RecoMCfromHFdecayCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
                }
              }
            } else if (prodRadius2 > 1. && prodRadius2 < 225.) {
              if (track.isGlobalTrack())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb1())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb2())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb3())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb4())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb5())
                histos.fill(HIST("Tracks/RecoMCRad1to15cmCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
            }
            if (prodRadius2 > 1.e-8 && prodRadius2 < 0.25) {
              if (track.isGlobalTrack())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB0"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb1())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB1"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb2())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB2"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb3())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB3"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb4())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB4"), track.pt(), track.eta(), track.phi());
              if (track.trackCutFlagFb5())
                histos.fill(HIST("Tracks/RecoMCRad1mumto5mmCollMatch/histptFB5"), track.pt(), track.eta(), track.phi());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processRecoMC, "process mc", true);

  void processMC(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {

    histos.fill(HIST("EventProp/histMCcollZ"), mcCollision.posZ());
    if (std::abs(mcCollision.posZ()) > 10.)
      return;
    ncollisionCounter++;
    for (auto& mcpart : mcParticles) {
      double prodRadius2 = mcpart.vx() * mcpart.vx() + mcpart.vy() * mcpart.vy();
      if (std::abs(mcpart.pdgCode()) == 211 || std::abs(mcpart.pdgCode()) == 321 || std::abs(mcpart.pdgCode()) == 2212 || std::abs(mcpart.pdgCode()) == 11 || std::abs(mcpart.pdgCode()) == 13) {
        int isHF = RecoDecay::getCharmHadronOrigin(mcParticles, mcpart, false);
        if (mcpart.isPhysicalPrimary()) {
          if (std::abs(mcpart.eta()) < 0.9) {
            histos.fill(HIST("Tracks/MCgen/histMCgenpt"), mcpart.pt());
          }
          histos.fill(HIST("Tracks/MCgen/histMCgen3dPhysPrimary"), mcpart.pt(), mcpart.eta(), mcpart.phi());
          if (isHF == RecoDecay::OriginType::Prompt || isHF == RecoDecay::OriginType::NonPrompt)
            histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedfromHFdecay"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        } else if (prodRadius2 > 1. && prodRadius2 < 225.) {
          histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedProdRad1to15cm"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        }
        if (prodRadius2 > 1.e-8 && prodRadius2 < 0.25) {
          histos.fill(HIST("Tracks/MCgen/histMCgen3dChargedProdRad1mumto5mm"), mcpart.pt(), mcpart.eta(), mcpart.phi());
        }
      }
    }
  }
  PROCESS_SWITCH(CheckFilterBit, processMC, "process mc", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CheckFilterBit>(cfgc)

  };
}

// void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
//{
//   ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Use MC info"}};
//   workflowOptions.push_back(optionDoMC);
// }
