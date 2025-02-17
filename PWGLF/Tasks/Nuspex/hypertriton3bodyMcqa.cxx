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
//
/// \file hypertriton3bodyMcqa.cxx
/// \brief QA for MC productions which contain hypertriton 3body decay process, including special checks for TOF PID
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include <cmath>
#include <array>
#include <cstdlib>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>

#include "CommonDataFormat/InteractionRecord.h"
#include "CommonDataFormat/IRFrame.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/pidTOFGeneric.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using ColwithEvTimes = o2::soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels, aod::EvTimeTOFFT0>;
using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::pidTOFFullDe, aod::TOFEvTime, aod::TOFSignal, aod::EvTimeTOFFT0ForTrack>;
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

template <class TMCTrackTo, typename TMCParticle>
bool is3bodyDecayedH3L(TMCParticle const& particle)
{
  if (std::abs(particle.pdgCode()) != 1010010030) {
    return false;
  }
  bool haveProton = false, havePion = false, haveDeuteron = false;
  bool haveAntiProton = false, haveAntiPion = false, haveAntiDeuteron = false;
  for (const auto& mcDaughter : particle.template daughters_as<TMCTrackTo>()) {
    if (mcDaughter.pdgCode() == 2212)
      haveProton = true;
    if (mcDaughter.pdgCode() == -2212)
      haveAntiProton = true;
    if (mcDaughter.pdgCode() == 211)
      havePion = true;
    if (mcDaughter.pdgCode() == -211)
      haveAntiPion = true;
    if (mcDaughter.pdgCode() == 1000010020)
      haveDeuteron = true;
    if (mcDaughter.pdgCode() == -1000010020)
      haveAntiDeuteron = true;
  }
  if (haveProton && haveAntiPion && haveDeuteron && particle.pdgCode() > 0) {
    return true;
  } else if (haveAntiProton && havePion && haveAntiDeuteron && particle.pdgCode() < 0) {
    return true;
  }
  return false;
}

template <class TMCTrackTo, typename TMCParticle>
bool isPairedH3LDaughters(TMCParticle const& mctrack0, TMCParticle const& mctrack1, TMCParticle const& mctrack2)
{
  for (const auto& particleMother : mctrack0.template mothers_as<TMCTrackTo>()) {
    if (!(particleMother.pdgCode() == 1010010030 && mctrack0.pdgCode() == 2212 && mctrack1.pdgCode() == -211 && mctrack2.pdgCode() == 1000010020) &&
        !(particleMother.pdgCode() == -1010010030 && mctrack0.pdgCode() == -2212 && mctrack1.pdgCode() == 211 && mctrack2.pdgCode() == -1000010020)) {
      continue;
    }
    bool flag1 = false, flag2 = false;
    for (const auto& mcDaughter : particleMother.template daughters_as<TMCTrackTo>()) {
      if (mcDaughter.globalIndex() == mctrack1.globalIndex())
        flag1 = true;
      if (mcDaughter.globalIndex() == mctrack2.globalIndex())
        flag2 = true;
    }
    if (!flag1 || !flag2)
      continue;
    // move the requirement in mass region into the loop to draw a histogram
    // double hypertritonMCMass = RecoDecay::m(array{array{mctrack0.px(), mctrack0.py(), mctrack0.pz()}, array{mctrack1.px(), mctrack1.py(), mctrack1.pz()}, array{mctrack2.px(), mctrack2.py(), mctrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
    // if (hypertritonMCMass > 2.990 && hypertritonMCMass < 2.993)
    return true;
  }
  return false;
}

// check the properties of daughters candidates and true daughters
struct Hypertriton3bodyMcqa {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Preslice<MCLabeledTracksIU> perCollisionTracks = aod::track::collisionId;

  int mRunNumber;

  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hParticleCounter", "hParticleCounter", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},

      {"hTPCNCls", "hTPCNCls", {HistType::kTH1F, {{160, 0.0f, 160.0f}}}},
      {"hTPCNClsCrossedRows", "hTPCNClsCrossedRows", {HistType::kTH1F, {{160, 0.0f, 160.0f}}}},
      {"hTrackEta", "hTrackEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hTrackITSNcls", "hTrackITSNcls", {HistType::kTH1F, {{10, 0.0f, 10.0f}}}},
      {"hTrackMcRapidity", "hTrackMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hTrackNsigmaProton", "hTrackNsigmaProton", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hTrackNsigmaPion", "hTrackNsigmaPion", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hTrackNsigmaDeuteron", "hTrackNsigmaDeuteron", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},

      {"hHypertritonEta", "hHypertritomEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hHypertritonMcRapidity", "hHypertritonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hHypertritonMcPt", "hHypertritonMcPt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},

      {"hProtonCounter", "hProtonCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hProtonPt", "hProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonP", "hProtonP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonMcPt", "hProtonMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonMcP", "hProtonMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hProtonEta", "hProtonEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hProtonMcRapidity", "hProtonMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hProtonNsigmaProton", "hProtonNsigmaProton", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hProtonTPCNCls", "hProtonTPCNCls", {HistType::kTH1F, {{120, 0.0f, 120.0f}}}},
      {"hProtonTPCBB", "hProtonTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hProtonTPCBBAfterTPCNclsCut", "hProtonTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDauProtonPt", "hDauProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauProtonMcPt", "hDauProtonMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauProtonNsigmaProton", "hDauProtonNsigmaProton", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hDauProtonTPCVsPt", "hDauProtonTPCVsPt", {HistType::kTH2F, {{50, 0.0f, 5.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},

      {"hPionCounter", "hPionCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hPionPt", "hPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionP", "hPionP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionMcPt", "hPionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionMcP", "hPionMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPionEta", "hPionEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hPionMcRapidity", "hPionMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hPionNsigmaPion", "hPionNsigmaPion", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hPionTPCNCls", "hPionTPCNCls", {HistType::kTH1F, {{160, 0.0f, 160.0f}}}},
      {"hPionTPCBB", "hPionTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hPionTPCBBAfterTPCNclsCut", "hPionTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDauPionPt", "hDauPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauPionMcPt", "hDauPionMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauPionNsigmaPion", "hDauPionNsigmaPion", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hDauPionTPCVsPt", "hDauPionTPCVsPt", {HistType::kTH2F, {{20, 0.0f, 2.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDauPionDcaXY", "hDauPionDcaXY", {HistType::kTH1F, {{100, -10.0f, 10.0f}}}},

      {"hDeuteronCounter", "hDeuteronCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},
      {"hDeuteronPt", "hDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronP", "hDeuteronP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronMcPt", "hDeuteronMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronMcP", "hDeuteronMcP", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDeuteronEta", "hDeuteronEta", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hDeuteronMcRapidity", "hDeuteronMcRapidity", {HistType::kTH1F, {{200, -10.0f, 10.0f}}}},
      {"hDeuteronNsigmaDeuteron", "hDeuteronNsigmaDeuteron", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hDeuteronTPCNCls", "hDeuteronTPCNCls", {HistType::kTH1F, {{120, 0.0f, 120.0f}}}},
      {"hDeuteronTPCBB", "hDeuteronTPCBB", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDeuteronTPCBBAfterTPCNclsCut", "hDeuteronTPCBBAfterTPCNclsCut", {HistType::kTH2F, {{320, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDauDeuteronPt", "hDauDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauDeuteronMcPt", "hDauDeuteronMcPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDauDeuteronTPCVsPt", "hDauDeuteronTPCVsPt", {HistType::kTH2F, {{80, 0.0f, 8.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsP", "hDauDeuteronTOFNSigmaVsP", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsPHasTOF", "hDauDeuteronTOFNSigmaVsPHasTOF", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronMatchCounter", "hDauDeuteronMatchCounter", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},

      {"hTPCBB", "hTPCBB", {HistType::kTH2F, {{120, -8.0f, 8.0f, "p/z(GeV/c)"}, {100, 0.0f, 1000.0f, "TPCSignal"}}}},

      {"hPairedH3LDaughers", "hPairedH3LDaughers", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hPairedH3LDaughersInvMass", "hPairedH3LDaughersInvMass", {HistType::kTH1F, {{300, 2.9f, 3.2f}}}},
      {"hDuplicatedH3LDaughers", "hDuplicatedH3LDaughers", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},

      // Diff checks always requir hasTOF
      {"hDiffTrackTOFSignal", "hDiffTrackTOFSignal", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"hDiffEvTimeForTrack", "hDiffEvTimeForTrack", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"hDiffTrackTOFNSigmaDe", "hDiffTrackTOFNSigmaDe", {HistType::kTH1F, {{200, -1.0f, 1.0f}}}},
      {"hDauDeuteronNewTOFNSigmaVsP", "hDauDeuteronNewTOFNSigmaVsP", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hWrongDeuteronTOFNSigmaVsP", "hWrongDeuteronTOFNSigmaVsP", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hWrongDeuteronNewTOFNSigmaVsP", "hWrongDeuteronNewTOFNSigmaVsP", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDiffColTime", "hDiffColTime", {HistType::kTH1F, {{200, -100.0f, 100.0f}}}},
      {"hDauDeuteronDiffTOFNsigmaDeHasTOF", "hDauDeuteronDiffTOFNsigmaDeHasTOF", {HistType::kTH1F, {{200, -100.0f, 100.0f}}}},

      // _v2 for using relinked collision
      {"hDauDeuteronTOFNSigmaVsP_CorrectCol", "hDauDeuteronTOFNSigmaVsP_CorrectCol", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronNewTOFNSigmaVsP_CorrectCol", "hDauDeuteronNewTOFNSigmaVsP_CorrectCol", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsP_v2", "hDauDeuteronTOFNSigmaVsP_v2", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronNewTOFNSigmaVsP_v2_AO2D", "hDauDeuteronNewTOFNSigmaVsP_v2 AO2D", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronNewTOFNSigmaVsP_v2_EvSel", "hDauDeuteronNewTOFNSigmaVsP_v2 EvSel", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsColTimeRes_v2", "hDauDeuteronTOFNSigmaVsColTimeRes_v2", {HistType::kTH2F, {{100, 0.0f, 400.0f, "CollisionTimeRes(ns)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsColTimeRes_v2_AO2D", "hDauDeuteronTOFNSigmaVsColTimeRes_v2 AO2D", {HistType::kTH2F, {{100, 0.0f, 400.0f, "CollisionTimeRes(ns)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFNSigmaVsColTimeRes_v2_EvSel", "hDauDeuteronTOFNSigmaVsColTimeRes_v2 EvSel", {HistType::kTH2F, {{100, 0.0f, 400.0f, "CollisionTimeRes(ns)"}, {600, -300.0f, 300.0f, "TOF n#sigma"}}}},
      {"hDauDeuteronTOFPIDCounter", "hDauDeuteronTOFPIDCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
      {"hDauDeuteronTOFPIDCounter_CloseBC", "hDauDeuteronTOFPIDCounter CloseBC", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
    },
  };

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(1, "Readin");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(2, "Has_mcparticle");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(3, "Rapidity Cut");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(4, "McisHypertriton");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(5, "McisProton");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(6, "McisPion");
    registry.get<TH1>(HIST("hParticleCounter"))->GetXaxis()->SetBinLabel(7, "McisDeuteron");

    TString trackCounterbinLabel[2] = {"hasMom", "FromHypertriton"};
    for (int i{0}; i < 2; i++) {
      registry.get<TH1>(HIST("hProtonCounter"))->GetXaxis()->SetBinLabel(i + 1, trackCounterbinLabel[i]);
      registry.get<TH1>(HIST("hPionCounter"))->GetXaxis()->SetBinLabel(i + 1, trackCounterbinLabel[i]);
      registry.get<TH1>(HIST("hDeuteronCounter"))->GetXaxis()->SetBinLabel(i + 1, trackCounterbinLabel[i]);
    }
    registry.get<TH1>(HIST("hDuplicatedH3LDaughers"))->GetXaxis()->SetBinLabel(1, "proton");
    registry.get<TH1>(HIST("hDuplicatedH3LDaughers"))->GetXaxis()->SetBinLabel(2, "pion");
    registry.get<TH1>(HIST("hDuplicatedH3LDaughers"))->GetXaxis()->SetBinLabel(3, "deuteron");

    registry.get<TH1>(HIST("hDauDeuteronMatchCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hDauDeuteronMatchCounter"))->GetXaxis()->SetBinLabel(2, "correct collision");
    registry.get<TH1>(HIST("hDauDeuteronMatchCounter"))->GetXaxis()->SetBinLabel(3, "hasTOF");
    registry.get<TH1>(HIST("hDauDeuteronMatchCounter"))->GetXaxis()->SetBinLabel(4, "hasTOF & correct collsion");

    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter"))->GetXaxis()->SetBinLabel(1, "Origin |n#sigma| >= 5");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter"))->GetXaxis()->SetBinLabel(2, "BothBC work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter"))->GetXaxis()->SetBinLabel(3, "Only BCAO2D work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter"))->GetXaxis()->SetBinLabel(4, "Only BCEvSel work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter"))->GetXaxis()->SetBinLabel(5, "BothBC not work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter_CloseBC"))->GetXaxis()->SetBinLabel(1, "Origin |n#sigma| < 6");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter_CloseBC"))->GetXaxis()->SetBinLabel(2, "BothBC work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter_CloseBC"))->GetXaxis()->SetBinLabel(3, "Only BCAO2D work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter_CloseBC"))->GetXaxis()->SetBinLabel(4, "Only BCEvSel work");
    registry.get<TH1>(HIST("hDauDeuteronTOFPIDCounter_CloseBC"))->GetXaxis()->SetBinLabel(5, "BothBC not work");
  }

  Configurable<float> dcapiontopv{"dcapiontopv", .05, "DCA Pion To PV"};
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
  Configurable<bool> mc_event_selection{"mc_event_selection", true, "mc event selection count post kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};

  // CCDB TOF PID paras
  Configurable<int64_t> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};
  Configurable<std::string> paramFileName{"paramFileName", "", "Path to the parametrization object. If empty the parametrization is not taken from file"};
  Configurable<std::string> parametrizationPath{"parametrizationPath", "TOF/Calib/Params", "Path of the TOF parametrization on the CCDB or in the file, if the paramFileName is not empty"};
  Configurable<std::string> passName{"passName", "", "Name of the pass inside of the CCDB parameter collection. If empty, the automatically deceted from metadata (to be implemented!!!)"};
  Configurable<std::string> timeShiftCCDBPath{"timeShiftCCDBPath", "", "Path of the TOF time shift vs eta. If empty none is taken"};
  Configurable<bool> loadResponseFromCCDB{"loadResponseFromCCDB", false, "Flag to load the response from the CCDB"};
  Configurable<bool> fatalOnPassNotAvailable{"fatalOnPassNotAvailable", true, "Flag to throw a fatal if the pass is not available in the retrieved CCDB object"};

  o2::aod::pidtofgeneric::TofPidNewCollision<MCLabeledTracksIU::iterator> bachelorTOFPID;
  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    // Initial TOF PID Paras, copied from PIDTOF.h
    timestamp.value = bc.timestamp();
    ccdb->setTimestamp(timestamp.value);
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    // TODO: implement the automatic pass name detection from metadata
    if (passName.value == "") {
      passName.value = "unanchored"; // temporary default
      LOG(warning) << "Passed autodetect mode for pass, not implemented yet, waiting for metadata. Taking '" << passName.value << "'";
    }
    LOG(info) << "Using parameter collection, starting from pass '" << passName.value << "'";

    const std::string fname = paramFileName.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file " << fname << ", using param: " << parametrizationPath.value;
      if (1) {
        o2::tof::ParameterCollection paramCollection;
        paramCollection.loadParamFromFile(fname, parametrizationPath.value);
        LOG(info) << "+++ Loaded parameter collection from file +++";
        if (!paramCollection.retrieveParameters(mRespParamsV2, passName.value)) {
          if (fatalOnPassNotAvailable) {
            LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          } else {
            LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
          }
        } else {
          mRespParamsV2.setShiftParameters(paramCollection.getPars(passName.value));
          mRespParamsV2.printShiftParameters();
        }
      } else {
        mRespParamsV2.loadParamFromFile(fname.data(), parametrizationPath.value);
      }
    } else if (loadResponseFromCCDB) { // Loading it from CCDB
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << parametrizationPath.value << " for timestamp " << timestamp.value;
      o2::tof::ParameterCollection* paramCollection = ccdb->getForTimeStamp<o2::tof::ParameterCollection>(parametrizationPath.value, timestamp.value);
      paramCollection->print();
      if (!paramCollection->retrieveParameters(mRespParamsV2, passName.value)) { // Attempt at loading the parameters with the pass defined
        if (fatalOnPassNotAvailable) {
          LOGF(fatal, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        } else {
          LOGF(warning, "Pass '%s' not available in the retrieved CCDB object", passName.value.data());
        }
      } else { // Pass is available, load non standard parameters
        mRespParamsV2.setShiftParameters(paramCollection->getPars(passName.value));
        mRespParamsV2.printShiftParameters();
      }
    }
    mRespParamsV2.print();
    if (timeShiftCCDBPath.value != "") {
      if (timeShiftCCDBPath.value.find(".root") != std::string::npos) {
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Pos", true);
        mRespParamsV2.setTimeShiftParameters(timeShiftCCDBPath.value, "gmean_Neg", false);
      } else {
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/pos", timeShiftCCDBPath.value.c_str()), timestamp.value), true);
        mRespParamsV2.setTimeShiftParameters(ccdb->getForTimeStamp<TGraph>(Form("%s/neg", timeShiftCCDBPath.value.c_str()), timestamp.value), false);
      }
    }

    bachelorTOFPID.SetParams(mRespParamsV2);
  }

  struct Indexdaughters { // check duplicated paired daughters
    int64_t index0;
    int64_t index1;
    int64_t index2;
    bool operator==(const Indexdaughters& t) const
    {
      return (this->index0 == t.index0 && this->index1 == t.index1 && this->index2 == t.index2);
    }
  };

  void process(ColwithEvTimes const& collisions, MCLabeledTracksIU const& tracks, aod::McParticles const& /*particlesMC*/, aod::McCollisions const& /*mcCollisions*/, aod::BCsWithTimestamps const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      registry.fill(HIST("hEventCounter"), 0.5);
      if (mc_event_selection && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);
      if (event_posZ_selection && std::abs(collision.posZ()) > 10.f) { // 10cm
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2.5);

      std::vector<int> protons, pions, deuterons;                     // index for daughter tracks
      std::unordered_set<int64_t> set_proton, set_pion, set_deuteron; // check duplicated daughters
      int itrack = -1;

      auto coltracks = tracks.sliceBy(perCollisionTracks, collision.globalIndex());

      for (const auto& track : coltracks) {

        ++itrack;
        registry.fill(HIST("hParticleCounter"), 0.5);
        registry.fill(HIST("hTrackITSNcls"), track.itsNCls());
        registry.fill(HIST("hTPCNCls"), track.tpcNClsFound());
        registry.fill(HIST("hTPCNClsCrossedRows"), track.tpcNClsCrossedRows());
        registry.fill(HIST("hTrackNsigmaDeuteron"), track.tpcNSigmaDe());
        registry.fill(HIST("hTrackNsigmaProton"), track.tpcNSigmaPr());
        registry.fill(HIST("hTrackNsigmaPion"), track.tpcNSigmaPi());

        if (!track.has_mcParticle()) {
          continue;
        }
        auto mcparticle = track.mcParticle_as<aod::McParticles>();
        registry.fill(HIST("hTPCBB"), track.p() * track.sign(), track.tpcSignal());

        registry.fill(HIST("hParticleCounter"), 1.5);

        // if (TMath::Abs(mcparticle.y()) > 0.9) {continue;}
        registry.fill(HIST("hParticleCounter"), 2.5);
        registry.fill(HIST("hTrackEta"), track.eta());
        registry.fill(HIST("hTrackMcRapidity"), mcparticle.y());

        // Hypertriton detected directly
        if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hParticleCounter"), 3.5);
          registry.fill(HIST("hHypertritonMcPt"), mcparticle.pt());
          registry.fill(HIST("hHypertritonEta"), track.eta());
          registry.fill(HIST("hHypertritonMcRapidity"), mcparticle.y());
        }

        // Proton
        if (mcparticle.pdgCode() == 2212 || mcparticle.pdgCode() == -2212) {
          registry.fill(HIST("hParticleCounter"), 4.5);
          if (track.tpcNClsFound() > 70) {
            registry.fill(HIST("hProtonTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
          }

          if (mcparticle.has_mothers()) {
            registry.fill(HIST("hProtonCounter"), 0.5);
            for (const auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
              bool flag_H3L = is3bodyDecayedH3L<aod::McParticles>(particleMother);
              if (!flag_H3L) {
                continue;
              }
              protons.push_back(itrack);
              auto p = set_proton.insert(mcparticle.globalIndex());
              if (p.second == false)
                registry.fill(HIST("hDuplicatedH3LDaughers"), 0);
              registry.fill(HIST("hProtonCounter"), 1.5);
              registry.fill(HIST("hDauProtonPt"), track.pt());
              registry.fill(HIST("hDauProtonMcPt"), mcparticle.pt());
              registry.fill(HIST("hDauProtonNsigmaProton"), track.tpcNSigmaPr());
              registry.fill(HIST("hDauProtonTPCVsPt"), track.pt(), track.tpcNSigmaPr());
            }
          }

          registry.fill(HIST("hProtonMcPt"), mcparticle.pt());
          registry.fill(HIST("hProtonMcP"), mcparticle.p());
          registry.fill(HIST("hProtonPt"), track.pt());
          registry.fill(HIST("hProtonP"), track.p());

          registry.fill(HIST("hProtonNsigmaProton"), track.tpcNSigmaPr());
          registry.fill(HIST("hProtonTPCNCls"), track.tpcNClsFound());
          registry.fill(HIST("hProtonEta"), track.eta());
          registry.fill(HIST("hProtonMcRapidity"), mcparticle.y());
          registry.fill(HIST("hProtonTPCBB"), track.p() * track.sign(), track.tpcSignal());
        }

        // Pion
        if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
          registry.fill(HIST("hParticleCounter"), 5.5);
          if (track.tpcNClsFound() > 70) {
            registry.fill(HIST("hPionTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
          }

          if (mcparticle.has_mothers()) {
            registry.fill(HIST("hPionCounter"), 0.5);
            for (const auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
              bool flag_H3L = is3bodyDecayedH3L<aod::McParticles>(particleMother);
              if (!flag_H3L) {
                continue;
              }
              pions.push_back(itrack);
              auto p = set_pion.insert(mcparticle.globalIndex());
              if (p.second == false) {
                registry.fill(HIST("hDuplicatedH3LDaughers"), 1);
              }
              registry.fill(HIST("hPionCounter"), 1.5);
              registry.fill(HIST("hDauPionPt"), track.pt());
              registry.fill(HIST("hDauPionMcPt"), mcparticle.pt());
              registry.fill(HIST("hDauPionTPCVsPt"), track.pt(), track.tpcNSigmaPi());
              registry.fill(HIST("hDauPionDcaXY"), track.dcaXY());
            }
          }

          registry.fill(HIST("hPionMcPt"), mcparticle.pt());
          registry.fill(HIST("hPionMcP"), mcparticle.p());
          registry.fill(HIST("hPionPt"), track.pt());
          registry.fill(HIST("hPionP"), track.p());

          registry.fill(HIST("hPionNsigmaPion"), track.tpcNSigmaPi());
          registry.fill(HIST("hPionTPCNCls"), track.tpcNClsFound());
          registry.fill(HIST("hPionEta"), track.eta());
          registry.fill(HIST("hPionMcRapidity"), mcparticle.y());
          registry.fill(HIST("hPionTPCBB"), track.p() * track.sign(), track.tpcSignal());
        }

        float tofNsigmaDe = -999;
        static constexpr float kCSPEED = TMath::C() * 1.0e2f * 1.0e-12f; // c in cm/ps

        if (track.hasTOF() && track.has_collision()) {
          auto responseDe = o2::pid::tof::ExpTimes<MCLabeledTracksIU::iterator, o2::track::PID::Deuteron>();
          // float bachExpTime = track.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v

          float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
          float bachExpTime = track.length() * std::sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v
          float tofsignal = track.trackTime() * 1000 + bachExpTime;                                                                                        // in ps

          float expSigma = responseDe.GetExpectedSigma(mRespParamsV2, track, tofsignal, track.tofEvTimeErr());
          // tofNsigmaDe = (track.tofSignal() - track.tofEvTime() - responseDe.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
          tofNsigmaDe = (tofsignal - track.tofEvTime() - responseDe.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;
          // tofNsigmaDe = (tofsignal - track.evTimeForTrack() - responseDe.GetCorrectedExpectedSignal(mRespParamsV2, track)) / expSigma;

          if (collision.bcId() == collision.foundBCId()) {
            registry.fill(HIST("hDiffColTime"), track.tofEvTime() - collision.collisionTime());
          }

          // Assume deuterons linked to the correct collision, result of new TOF PID should be same as the default one
          registry.fill(HIST("hDiffTrackTOFSignal"), track.tofSignal() - tofsignal);
          registry.fill(HIST("hDiffEvTimeForTrack"), track.tofEvTime() - track.evTimeForTrack());
          registry.fill(HIST("hDiffTrackTOFNSigmaDe"), track.tofNSigmaDe() - bachelorTOFPID.GetTOFNSigma(o2::track::PID::Deuteron, track, collision, collision));
          // registry.fill(HIST("hDiffTrackTOFNSigmaDe"), track.tofExpSigmaDe() - bachelorTOFPID.GetTOFNSigma(o2::track::PID::Deuteron, track, collision, collision));
        }

        // Deuteron
        if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
          registry.fill(HIST("hParticleCounter"), 6.5);
          if (track.tpcNClsFound() > 70) {
            registry.fill(HIST("hDeuteronTPCBBAfterTPCNclsCut"), track.p() * track.sign(), track.tpcSignal());
          }

          if (mcparticle.has_mothers()) {
            registry.fill(HIST("hDeuteronCounter"), 0.5);
            for (const auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
              bool flag_H3L = is3bodyDecayedH3L<aod::McParticles>(particleMother);
              if (!flag_H3L) {
                continue;
              }
              deuterons.push_back(itrack);
              auto p = set_deuteron.insert(mcparticle.globalIndex());
              if (p.second == false)
                registry.fill(HIST("hDuplicatedH3LDaughers"), 2);
              registry.fill(HIST("hDeuteronCounter"), 1.5);
              registry.fill(HIST("hDauDeuteronPt"), track.pt());
              registry.fill(HIST("hDauDeuteronMcPt"), mcparticle.pt());
              registry.fill(HIST("hDauDeuteronTPCVsPt"), track.pt(), track.tpcNSigmaDe());
              registry.fill(HIST("hDauDeuteronTOFNSigmaVsP"), track.sign() * track.p(), track.tofNSigmaDe());

              registry.fill(HIST("hDauDeuteronNewTOFNSigmaVsP"), track.sign() * track.p(), tofNsigmaDe);
              if (track.hasTOF()) {
                registry.fill(HIST("hDauDeuteronTOFNSigmaVsPHasTOF"), track.sign() * track.p(), track.tofNSigmaDe());
                registry.fill(HIST("hDauDeuteronDiffTOFNsigmaDeHasTOF"), track.tofNSigmaDe() - tofNsigmaDe);
              }
              registry.fill(HIST("hDauDeuteronMatchCounter"), 0.5);
              if (mcparticle.mcCollisionId() == collision.mcCollisionId()) {
                registry.fill(HIST("hDauDeuteronMatchCounter"), 1.5);
              }
              if (track.hasTOF()) {
                registry.fill(HIST("hDauDeuteronMatchCounter"), 2.5);
                if (mcparticle.mcCollisionId() == collision.mcCollisionId()) {
                  registry.fill(HIST("hDauDeuteronMatchCounter"), 3.5);
                }
              }
            }
          }

          registry.fill(HIST("hDeuteronMcPt"), mcparticle.pt());
          registry.fill(HIST("hDeuteronMcP"), mcparticle.p());
          registry.fill(HIST("hDeuteronPt"), track.pt());
          registry.fill(HIST("hDeuteronP"), track.p());

          registry.fill(HIST("hDeuteronNsigmaDeuteron"), track.tpcNSigmaDe());
          registry.fill(HIST("hDeuteronTPCNCls"), track.tpcNClsFound());
          registry.fill(HIST("hDeuteronEta"), track.eta());
          registry.fill(HIST("hDeuteronMcRapidity"), mcparticle.y());
          registry.fill(HIST("hDeuteronTPCBB"), track.p() * track.sign(), track.tpcSignal());
        } else {
          if (track.hasTOF()) {
            registry.fill(HIST("hWrongDeuteronTOFNSigmaVsP"), track.sign() * track.p(), track.tofNSigmaDe());
            registry.fill(HIST("hWrongDeuteronNewTOFNSigmaVsP"), track.sign() * track.p(), tofNsigmaDe);
          }
        }
      }

      std::vector<Indexdaughters> set_pair;
      for (size_t iproton = 0; iproton < protons.size(); iproton++) {
        auto track0 = tracks.iteratorAt(protons[iproton]);
        auto mctrack0 = track0.mcParticle_as<aod::McParticles>();
        for (size_t ipion = 0; ipion < pions.size(); ipion++) {
          auto track1 = tracks.iteratorAt(pions[ipion]);
          auto mctrack1 = track1.mcParticle_as<aod::McParticles>();
          for (size_t ideuteron = 0; ideuteron < deuterons.size(); ideuteron++) {
            auto track2 = tracks.iteratorAt(deuterons[ideuteron]);
            auto mctrack2 = track2.mcParticle_as<aod::McParticles>();
            if (isPairedH3LDaughters<aod::McParticles>(mctrack0, mctrack1, mctrack2)) {
              registry.fill(HIST("hPairedH3LDaughers"), 0);
              // MC mass cut, to check if the daughters are from materials
              double hypertritonMCMass = RecoDecay::m(std::array{std::array{mctrack0.px(), mctrack0.py(), mctrack0.pz()}, std::array{mctrack1.px(), mctrack1.py(), mctrack1.pz()}, std::array{mctrack2.px(), mctrack2.py(), mctrack2.pz()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
              registry.fill(HIST("hPairedH3LDaughersInvMass"), hypertritonMCMass);
              if (hypertritonMCMass < 2.990 || hypertritonMCMass > 2.993)
                continue;
              registry.fill(HIST("hPairedH3LDaughers"), 1);
              // duplicated daughters check
              Indexdaughters temp = {mctrack0.globalIndex(), mctrack1.globalIndex(), mctrack2.globalIndex()};
              auto p = std::find(set_pair.begin(), set_pair.end(), temp);
              if (p == set_pair.end()) {
                set_pair.push_back(temp);
                registry.fill(HIST("hPairedH3LDaughers"), 2);
              }
            }
          }
        }
      }
    }

    // Check for recalculated TOF PID for secondary deuterons

    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }

    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcparticle = track.mcParticle_as<aod::McParticles>();
      if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
        if (!mcparticle.has_mothers()) {
          continue;
        }
        const auto evtReconstructed = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcparticle.mcCollision_as<aod::McCollisions>().globalIndex());
        if (evtReconstructed == SelectedEvents.end() || !track.has_collision()) {
          continue;
        }
        if (!track.has_collision()) {
          continue;
        }
        auto collision = collisions.iteratorAt(evtReconstructed - SelectedEvents.begin());
        auto originalcollision = track.collision_as<ColwithEvTimes>();

        for (const auto& particleMother : mcparticle.mothers_as<aod::McParticles>()) {
          bool flag_H3L = is3bodyDecayedH3L<aod::McParticles>(particleMother);
          if (!flag_H3L) {
            continue;
          }

          auto bc = collision.bc_as<aod::BCsWithTimestamps>();
          initCCDB(bc);
          float tofNsigmaDeAO2D = -999;
          float tofNsigmaDeEvSel = -999;

          if (track.hasTOF()) {
            /*auto responseDe = o2::pid::tof::ExpTimes<MCLabeledTracksIU::iterator, o2::track::PID::Deuteron>();
            //float bachExpTime = track.length() * sqrt((o2::constants::physics::MassDeuteron * o2::constants::physics::MassDeuteron) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v
            float mMassHyp = o2::track::pid_constants::sMasses2Z[track.pidForTracking()];
            float bachExpTime = track.length() * std::sqrt((mMassHyp * mMassHyp) + (track.tofExpMom() * track.tofExpMom())) / (kCSPEED * track.tofExpMom()); // L*E/(p*c) = L/v
             */

            tofNsigmaDeAO2D = bachelorTOFPID.GetTOFNSigma(o2::track::PID::Deuteron, track, originalcollision, collision);
            tofNsigmaDeEvSel = bachelorTOFPID.GetTOFNSigma(o2::track::PID::Deuteron, track, originalcollision, collision, false);

            if (collision.globalIndex() == originalcollision.globalIndex()) {
              registry.fill(HIST("hDauDeuteronTOFNSigmaVsP_CorrectCol"), track.sign() * track.p(), track.tofNSigmaDe());
              registry.fill(HIST("hDauDeuteronNewTOFNSigmaVsP_CorrectCol"), track.sign() * track.p(), tofNsigmaDeAO2D);
              continue;
            }

            /*if (originalcollision.collisionTimeRes() > 40){
              continue;
            }*/
            registry.fill(HIST("hDauDeuteronTOFNSigmaVsP_v2"), track.sign() * track.p(), track.tofNSigmaDe());
            registry.fill(HIST("hDauDeuteronNewTOFNSigmaVsP_v2_AO2D"), track.sign() * track.p(), tofNsigmaDeAO2D);
            registry.fill(HIST("hDauDeuteronNewTOFNSigmaVsP_v2_EvSel"), track.sign() * track.p(), tofNsigmaDeEvSel);
            registry.fill(HIST("hDauDeuteronTOFNSigmaVsColTimeRes_v2"), collision.collisionTimeRes(), track.tofNSigmaDe());
            registry.fill(HIST("hDauDeuteronTOFNSigmaVsColTimeRes_v2_AO2D"), originalcollision.collisionTimeRes(), tofNsigmaDeAO2D);
            registry.fill(HIST("hDauDeuteronTOFNSigmaVsColTimeRes_v2_EvSel"), originalcollision.collisionTimeRes(), tofNsigmaDeEvSel);

            if (std::abs(track.tofNSigmaDe()) >= 5) {
              registry.fill(HIST("hDauDeuteronTOFPIDCounter"), 0.5);
              if (std::abs(tofNsigmaDeAO2D) < 5 && std::abs(tofNsigmaDeEvSel) < 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter"), 1.5);
              } else if (std::abs(tofNsigmaDeAO2D) < 5 && std::abs(tofNsigmaDeEvSel) >= 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter"), 2.5);
              } else if (std::abs(tofNsigmaDeAO2D) >= 5 && std::abs(tofNsigmaDeEvSel) < 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter"), 3.5);
              } else if (std::abs(tofNsigmaDeAO2D) >= 5 && std::abs(tofNsigmaDeEvSel) >= 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter"), 4.5);
              }
            } else if (std::abs(track.tofNSigmaDe()) < 5) {
              registry.fill(HIST("hDauDeuteronTOFPIDCounter_CloseBC"), 0.5);
              if (std::abs(tofNsigmaDeAO2D) < 5 && std::abs(tofNsigmaDeEvSel) < 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter_CloseBC"), 1.5);
              } else if (std::abs(tofNsigmaDeAO2D) < 5 && std::abs(tofNsigmaDeEvSel) >= 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter_CloseBC"), 2.5);
              } else if (std::abs(tofNsigmaDeAO2D) >= 5 && std::abs(tofNsigmaDeEvSel) < 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter_CloseBC"), 3.5);
              } else if (std::abs(tofNsigmaDeAO2D) >= 5 && std::abs(tofNsigmaDeEvSel) >= 5) {
                registry.fill(HIST("hDauDeuteronTOFPIDCounter_CloseBC"), 4.5);
              }
            }
          }
        }
      }
    }
  }
};

// check the performance of mcparticle
struct Hypertriton3bodyMcParticleCheck {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hMcCollCounter", "hMcCollCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}}},

      {"h3dMCDecayedHypertriton", "h3dMCDecayedHypertriton", {HistType::kTH3F, {{20, -1.0f, 1.0f, "Rapidity"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {50, 0.0f, 50.0f, "ct(cm)"}}}},

      {"hMcHypertritonCounter", "hMcHypertritonCounter", {HistType::kTH1F, {{9, 0.0f, 9.0f}}}},
      {"hMcHypertritonPt", "hMcHypertritonPt", {HistType::kTH1F, {{300, 0.0f, 15.0f}}}},
      {"hMcProtonPt", "hMcProtonPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hMcPionPt", "hMcPionPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hMcDeuteronPt", "hMcDeuteronPt", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},

      {"hMcRecoInvMass", "hMcRecoInvMass", {HistType::kTH1F, {{100, 2.95, 3.05f}}}},

      {"hDiffDaughterR", "hDiffDaughterR", {HistType::kTH1F, {{10000, -100, 100}}}}, // difference between minR of pion&proton and R of deuteron(bachelor)
      {"hTrackX", "hTrackX", {HistType::kTH1F, {{10000, -100, 100}}}},
      {"hTrackY", "hTrackY", {HistType::kTH1F, {{10000, -100, 100}}}},
      {"hTrackZ", "hTrackZ", {HistType::kTH1F, {{10000, -100, 100}}}},
    },
  };

  o2::pid::tof::TOFResoParamsV2 mRespParamsV2;

  void init(InitContext&)
  {
    registry.get<TH1>(HIST("hMcCollCounter"))->GetXaxis()->SetBinLabel(1, "Total Counter");
    registry.get<TH1>(HIST("hMcCollCounter"))->GetXaxis()->SetBinLabel(2, "Reconstructed");

    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "Hypertriton All");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "Matter All");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(3, "AntiMatter All");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(4, "confirm to 3-body decay");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(5, "Matter");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(6, "AntiMatter");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(7, "Rapidity");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(8, "Lifetime");
    registry.get<TH1>(HIST("hMcHypertritonCounter"))->GetXaxis()->SetBinLabel(9, "PtCut");
  }

  Configurable<bool> mc_event_selection{"mc_event_selection", true, "mc event selection count post kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};

  Preslice<aod::McParticles> permcCollision = o2::aod::mcparticle::mcCollisionId;

  std::vector<int64_t> mcPartIndices;
  template <typename TTrackTable>
  void SetTrackIDForMC(aod::McParticles const& particlesMC, TTrackTable const& tracks)
  {
    mcPartIndices.clear();
    mcPartIndices.resize(particlesMC.size());
    std::fill(mcPartIndices.begin(), mcPartIndices.end(), -1);
    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcparticle = track.template mcParticle_as<aod::McParticles>();
        if (mcPartIndices[mcparticle.globalIndex()] == -1) {
          mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
        } else {
          auto candTrack = tracks.rawIteratorAt(mcPartIndices[mcparticle.globalIndex()]);
          // Use the track which has innest information (also best quality?
          if (track.x() < candTrack.x()) {
            mcPartIndices[mcparticle.globalIndex()] = track.globalIndex();
          }
        }

        // Checks for TrackR
        registry.fill(HIST("hTrackX"), track.x());
        registry.fill(HIST("hTrackY"), track.y());
        registry.fill(HIST("hTrackZ"), track.z());
      }
    }
  }

  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& particlesMC, const o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>& collisions, MCLabeledTracksIU const& tracks)
  {
    SetTrackIDForMC(particlesMC, tracks);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (mc_event_selection && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      if (event_posZ_selection && std::abs(collision.posZ()) > 10.f) { // 10cm
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);

    for (const auto& mcCollision : mcCollisions) {
      registry.fill(HIST("hMcCollCounter"), 0.5);
      const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();
      if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
        continue;
      }
      registry.fill(HIST("hMcCollCounter"), 1.5);

      const auto& dparticlesMC = particlesMC.sliceBy(permcCollision, mcCollision.globalIndex());

      for (const auto& mcparticle : dparticlesMC) {

        if (mcparticle.pdgCode() == 2212 || mcparticle.pdgCode() == -2212) {
          registry.fill(HIST("hMcProtonPt"), mcparticle.pt());
        }
        if (mcparticle.pdgCode() == 211 || mcparticle.pdgCode() == -211) {
          registry.fill(HIST("hMcPionPt"), mcparticle.pt());
        }
        if (mcparticle.pdgCode() == 1000010020 || mcparticle.pdgCode() == -1000010020) {
          registry.fill(HIST("hMcDeuteronPt"), mcparticle.pt());
        }
        if (mcparticle.pdgCode() == 1010010030) {
          registry.fill(HIST("hMcHypertritonCounter"), 1.5);
        } else if (mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hMcHypertritonCounter"), 2.5);
        }
        if (mcparticle.pdgCode() == 1010010030 || mcparticle.pdgCode() == -1010010030) {
          registry.fill(HIST("hMcHypertritonCounter"), 0.5);
          registry.fill(HIST("hMcHypertritonPt"), mcparticle.pt());

          double dauDeuteronPos[3] = {-999, -999, -999};
          double dauProtonMom[3] = {-999, -999, -999};
          double dauPionMom[3] = {-999, -999, -999};
          double dauDeuteronMom[3] = {-999, -999, -999};
          double MClifetime = 999;
          double dauProtonTrackR = 9999, dauPionTrackR = 99999, dauDeuteronTrackR = 999999;
          bool flag_H3L = is3bodyDecayedH3L<aod::McParticles>(mcparticle);
          if (!flag_H3L) {
            continue;
          }
          for (const auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
            if (std::abs(mcparticleDaughter.pdgCode()) == 2212) {
              dauProtonMom[0] = mcparticleDaughter.px();
              dauProtonMom[1] = mcparticleDaughter.py();
              dauProtonMom[2] = mcparticleDaughter.pz();
              if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
                auto trackProton = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
                dauProtonTrackR = trackProton.x();
              }
            }
            if (std::abs(mcparticleDaughter.pdgCode()) == 211) {
              dauPionMom[0] = mcparticleDaughter.px();
              dauPionMom[1] = mcparticleDaughter.py();
              dauPionMom[2] = mcparticleDaughter.pz();
              if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
                auto trackPion = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
                dauPionTrackR = trackPion.x();
              }
            }
            if (std::abs(mcparticleDaughter.pdgCode()) == 1000010020) {
              dauDeuteronPos[0] = mcparticleDaughter.vx();
              dauDeuteronPos[1] = mcparticleDaughter.vy();
              dauDeuteronPos[2] = mcparticleDaughter.vz();
              dauDeuteronMom[0] = mcparticleDaughter.px();
              dauDeuteronMom[1] = mcparticleDaughter.py();
              dauDeuteronMom[2] = mcparticleDaughter.pz();
              if (mcPartIndices[mcparticleDaughter.globalIndex()] != -1) {
                auto trackDeuteron = tracks.rawIteratorAt(mcPartIndices[mcparticleDaughter.globalIndex()]);
                dauDeuteronTrackR = trackDeuteron.x();
              }
            }
          }
          if (mcparticle.pdgCode() == 1010010030) {
            registry.fill(HIST("hMcHypertritonCounter"), 3.5);
            registry.fill(HIST("hMcHypertritonCounter"), 4.5);
          }
          if (mcparticle.pdgCode() == -1010010030) {
            registry.fill(HIST("hMcHypertritonCounter"), 3.5);
            registry.fill(HIST("hMcHypertritonCounter"), 5.5);
          }
          double hypertritonMCMass = RecoDecay::m(std::array{std::array{dauProtonMom[0], dauProtonMom[1], dauProtonMom[2]}, std::array{dauPionMom[0], dauPionMom[1], dauPionMom[2]}, std::array{dauDeuteronMom[0], dauDeuteronMom[1], dauDeuteronMom[2]}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
          registry.fill(HIST("hMcRecoInvMass"), hypertritonMCMass);

          if (hypertritonMCMass > 2.990 && hypertritonMCMass < 2.993) {
            MClifetime = RecoDecay::sqrtSumOfSquares(dauDeuteronPos[0] - mcparticle.vx(), dauDeuteronPos[1] - mcparticle.vy(), dauDeuteronPos[2] - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
            registry.fill(HIST("h3dMCDecayedHypertriton"), mcparticle.y(), mcparticle.pt(), MClifetime);

            double diffTrackR = dauDeuteronTrackR - std::min(dauPionTrackR, dauProtonTrackR);
            registry.fill(HIST("hDiffDaughterR"), diffTrackR);

            // int daughterPionCount = 0;
            // for (auto& mcparticleDaughter : mcparticle.daughters_as<aod::McParticles>()) {
            //   if (std::abs(mcparticleDaughter.pdgCode()) == 211) {
            //     daughterPionCount++;
            //   }
            // }

            // Counter for hypertriton N_gen
            if (std::abs(mcparticle.y()) < 1) {
              registry.fill(HIST("hMcHypertritonCounter"), 6.5);
              if (MClifetime < 40) {
                registry.fill(HIST("hMcHypertritonCounter"), 7.5);
                if (mcparticle.pt() > 1 && mcparticle.pt() < 10) {
                  registry.fill(HIST("hMcHypertritonCounter"), 8.5);
                }
              }
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Hypertriton3bodyMcqa>(cfgc),
    adaptAnalysisTask<Hypertriton3bodyMcParticleCheck>(cfgc),
  };
}
