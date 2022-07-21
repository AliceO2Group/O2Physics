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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGUD/DataModel/UDTables.h"

#include <TH1D.h>
#include "TLorentzVector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandAnalyzer {
  Preslice<o2::aod::UDMcParticles> perMcCollision = o2::aod::udmcparticle::udMcCollisionId;

  float ft0DummyTime = 32.767f;
  float mmuon = 0.1056583745;
  float muonPDG = 13;

  HistogramRegistry registry{
    "registry",
    {{"TracksWithFT0/Eta", ";#eta;", {HistType::kTH1D, {{80, -4., 4.}}}},
     {"TracksWithFT0/TimeFT0A", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     {"TracksWithFT0/TimeFT0C", ";time C, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     //
     {"PairsMC/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairsMC/EtaMu", ";#eta_{#mu};", {HistType::kTH1D, {{100, -6., 6.}}}},
     //
     {"PairSelection/MassRaw", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtRaw", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassMC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/Pt", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/EtaMu", ";#eta_{#mu};", {HistType::kTH1D, {{100, -6., 6.}}}},
     //
     {"PairSelection/MassNoFT0", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0Contam", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0MC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtNoFT0", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}}}};

  void collectMC(o2::aod::UDMcCollisions const& mcCollisions, o2::aod::UDMcParticles const& mcParticles)
  {
    int32_t nMCEvents = mcCollisions.size();
    // collect MC distributions
    for (int32_t i = 0; i < nMCEvents; i++) {
      auto mcPartsGroup = mcParticles.sliceBy(perMcCollision, i);
      std::vector<aod::UDMcParticle> mcPartsFiltered;
      for (const auto& mcPart : mcPartsGroup) {
        if (!mcPart.isPhysicalPrimary() || std::abs(mcPart.pdgCode()) != muonPDG) {
          continue;
        }
        mcPartsFiltered.emplace_back(mcPart);
      }
      // sanity check
      if (mcPartsFiltered.size() != 2) {
        continue;
      }
      const auto& part1 = mcPartsFiltered[0];
      const auto& part2 = mcPartsFiltered[1];
      TLorentzVector p1, p2, p;
      p1.SetXYZM(part1.px(), part1.py(), part1.pz(), mmuon);
      p2.SetXYZM(part2.px(), part2.py(), part2.pz(), mmuon);
      p = p1 + p2;
      registry.fill(HIST("PairsMC/Mass"), p.M());
      registry.fill(HIST("PairsMC/EtaMu"), p1.Eta());
      registry.fill(HIST("PairsMC/EtaMu"), p2.Eta());
      mcPartsFiltered.clear();
    }
  }

  template <typename TTrack1, typename TTrack2>
  void processCandidate(o2::aod::UDCollision const& cand, o2::aod::UDMcParticles const& mcParticles, TTrack1& tr1, TTrack2& tr2)
  {
    TLorentzVector p1, p2, pPair;
    p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mmuon);
    p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mmuon);
    pPair = p1 + p2;
    // distributions without selection
    registry.fill(HIST("PairSelection/MassRaw"), pPair.M());
    registry.fill(HIST("PairSelection/PtRaw"), pPair.Pt());
    // unlike-charged tracks requirement
    if (tr1.sign() * tr2.sign() >= 0) {
      return;
    }
    int mcPartId1 = tr1.udMcParticleId();
    int mcPartId2 = tr2.udMcParticleId();
    const auto& mcPart1 = mcParticles.iteratorAt(mcPartId1);
    const auto& mcPart2 = mcParticles.iteratorAt(mcPartId2);
    TLorentzVector mcP1, mcP2, pPairMC;
    mcP1.SetXYZM(mcPart1.px(), mcPart1.py(), mcPart1.pz(), mmuon);
    mcP2.SetXYZM(mcPart2.px(), mcPart2.py(), mcPart2.pz(), mmuon);
    pPairMC = mcP1 + mcP2;
    // reconstructed distributions and "MC" distributions for reconstructed tracks
    registry.fill(HIST("PairSelection/Mass"), pPair.M());
    registry.fill(HIST("PairSelection/Pt"), pPair.Pt());
    registry.fill(HIST("PairSelection/EtaMu"), p1.Eta());
    registry.fill(HIST("PairSelection/EtaMu"), p2.Eta());
    registry.fill(HIST("PairSelection/MassMC"), pPairMC.M());
    // check FT0 signal
    if (cand.hasFT0()) {
      registry.fill(HIST("TracksWithFT0/Eta"), p1.Eta());
      registry.fill(HIST("TracksWithFT0/Eta"), p2.Eta());
      registry.fill(HIST("TracksWithFT0/TimeFT0A"), cand.timeFT0A());
      registry.fill(HIST("TracksWithFT0/TimeFT0C"), cand.timeFT0C());
      // if there is a signal, candidate passes if timeA is dummy
      // and timeC is between +/- 1 ns
      bool checkA = std::abs(cand.timeFT0A() - ft0DummyTime) < 1e-3;
      bool checkC = cand.timeFT0C() > -1. && cand.timeFT0C() < 1.;
      if (!(checkA && checkC)) {
        return;
      }
    }
    registry.fill(HIST("PairSelection/MassNoFT0"), pPair.M());
    registry.fill(HIST("PairSelection/MassNoFT0MC"), pPairMC.M());
    registry.fill(HIST("PairSelection/PtNoFT0"), pPair.Pt());
    if (mcPartId1 == -1 || mcPartId2 == -1) {
      registry.fill(HIST("PairSelection/MassNoFT0Contam"), pPair.M());
    }
  }

  // process candidates with 2 muon tracks
  void processFwdMC(o2::aod::UDCollisions const& eventCandidates,
                    soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTrackCollisionIDs, o2::aod::UDFwdTracksExtra, o2::aod::UDMcFwdTrackLabels> const& muonTracks,
                    o2::aod::UDMcCollisions const& mcCollisions,
                    o2::aod::UDMcParticles const& mcParticles)
  {
    collectMC(mcCollisions, mcParticles);

    // assuming that candidates have exatly 2 muon tracks and 0 barrel tracks
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto muonTracksPerCand = muonTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      const auto& tr1 = muonTracksPerCand.iteratorAt(0);
      const auto& tr2 = muonTracksPerCand.iteratorAt(1);
      processCandidate(cand, mcParticles, tr1, tr2);
    }
  }

  // process candidates with 1 muon and 1 barrel tracks
  void processSemiFwdMC(o2::aod::UDCollisions const& eventCandidates,
                        soa::Join<o2::aod::UDFwdTracks, o2::aod::UDFwdTrackCollisionIDs, o2::aod::UDFwdTracksExtra, o2::aod::UDMcFwdTrackLabels> const& muonTracks,
                        soa::Join<o2::aod::UDTracks, o2::aod::UDTrackCollisionIDs, o2::aod::UDTracksExtra, o2::aod::UDMcTrackLabels> const& barTracks,
                        o2::aod::UDMcCollisions const& mcCollisions,
                        o2::aod::UDMcParticles const& mcParticles)
  {
    collectMC(mcCollisions, mcParticles);

    // assuming that candidates have exatly 1 muon track and 1 barrel track
    for (const auto& cand : eventCandidates) {
      auto candID = cand.globalIndex();
      auto muonTracksPerCand = muonTracks.select(o2::aod::udfwdtrack::udCollisionId == candID);
      auto barTracksPerCand = barTracks.select(o2::aod::udtrack::udCollisionId == candID);
      const auto& tr1 = muonTracksPerCand.iteratorAt(0);
      const auto& tr2 = barTracksPerCand.iteratorAt(0);
      processCandidate(cand, mcParticles, tr1, tr2);
    }
  }

  // todo: add processes without MC information
  PROCESS_SWITCH(UpcCandAnalyzer, processFwdMC, "Analyse forward candidates", false);
  PROCESS_SWITCH(UpcCandAnalyzer, processSemiFwdMC, "Analyse semiforward candidates", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcCandAnalyzer>(cfgc)};
}
