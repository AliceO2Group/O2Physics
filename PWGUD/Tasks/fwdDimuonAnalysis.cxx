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

#include <TH1D.h>
#include <TH2D.h>
#include "TLorentzVector.h"

#include "Common/DataModel/EventSelection.h"
#include "PWGUD/DataModel/UDTables.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define mmuon 0.1056583745

uint64_t relativeTime_to_LocalBC(double relativeTimeStampInNS)
{
  return relativeTimeStampInNS > 0. ? std::round(relativeTimeStampInNS / o2::constants::lhc::LHCBunchSpacingNS) : 0;
}

struct AnalysisTracksMc {
  Configurable<int> signalSourceID{"signalSourceID", 1, "signal source ID"};
  Configurable<int> selPDG{"selPDG", 13, "abs(pdg) of primary particles of interest"};

  // mc pt distribution from mcparticle
  OutputObj<TH1D> hPtMC{TH1D("hPt", ";pt;", 100, 0., 10.)};
  // mc eta distribution from mcparticle
  OutputObj<TH1D> hEtaMC{TH1D("hEta", ";#eta;", 3000, -4., 4.)};
  // mc pair distributions
  OutputObj<TH1D> hPairMassMC{TH1D("hPairMass", ";m;", 100, 0., 10.)};
  OutputObj<TH1D> hPairPtMC{TH1D("hPairPt", ";m;", 100, 0., 10.)};

  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

  void process(aod::McCollisions const& mcCollisions,
               aod::McParticles_001 const& mcParticles,
               Tracks const& fwdTracks,
               aod::AmbiguousFwdTracks const& ambFwdTracks,
               aod::BCs const& bcs)
  {
    // mc collisions
    // -----------------------------------------------------------
    for (const auto& mcCol : mcCollisions) {
      if (mcCol.generatorsID() != signalSourceID)
        continue;
      auto mcParts = mcParticles.sliceBy(aod::mcparticle::mcCollisionId, mcCol.globalIndex());
      std::vector<aod::McParticles_001::iterator> mcPartsFiltered;
      for (const auto& mcPart : mcParts) {
        if (!mcPart.isPhysicalPrimary() || std::abs(mcPart.pdgCode()) != selPDG)
          continue;
        mcPartsFiltered.emplace_back(mcPart);
      }
      if (mcPartsFiltered.size() != 2)
        continue;
      const auto& part1 = mcPartsFiltered[0];
      const auto& part2 = mcPartsFiltered[1];
      TLorentzVector p1, p2, p;
      p1.SetXYZM(part1.px(), part1.py(), part1.pz(), mmuon);
      p2.SetXYZM(part2.px(), part2.py(), part2.pz(), mmuon);
      p = p1 + p2;
      hPairMassMC->Fill(p.M());
      hPairPtMC->Fill(p.Pt());
      mcPartsFiltered.clear();
    }

    // mc particles
    // -----------------------------------------------------------
    for (const auto& mcPart : mcParticles) {
      float pt = mcPart.pt();
      float eta = mcPart.eta();
      // signals
      if (mcPart.mcCollision().generatorsID() == signalSourceID &&
          (mcPart.isPhysicalPrimary() && std::abs(mcPart.pdgCode()) == selPDG)) {
        hPtMC->Fill(pt);
        hEtaMC->Fill(eta);
      }
    }
  }
};

const char* fwdTrackLabels[5] = {"GlobalMuonTrack",
                                 "GlobalMuonTrackOtherMatch",
                                 "GlobalForwardTrack",
                                 "MuonStandaloneTrack",
                                 "MCHStandaloneTrack"};

struct AnalysisTracks {
  Configurable<int> signalSourceID{"signalSourceID", 1, "signal source ID"};
  Configurable<int> filterFwdType{"filterFwdType", -1, "fwd track type for skimming"};
  Configurable<int> selPDG{"selPDG", 13, "abs(pdg) of primary particles of interest"};

  HistogramRegistry registry{
    "registry",
    {{"TrackTypes/TrackTypesAll", ";;", {HistType::kTH1I, {{5, 0, 5}}}},
     {"TrackTypes/AmbTrackTypesAll", ";;", {HistType::kTH1I, {{5, 0, 5}}}},
     {"TrackTypes/TrackTypesSig", ";;", {HistType::kTH1I, {{5, 0, 5}}}},
     {"TrackTypes/AmbTrackTypesSig", ";;", {HistType::kTH1I, {{5, 0, 5}}}},
     //
     {"TrackTypes/NTracksPerMcPart", ";;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"TrackTypes/NGlobalMuonTracksPerMcPart", ";;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"TrackTypes/NMuonStandaloneTracksPerMcPart", ";;", {HistType::kTH1I, {{100, 0, 100}}}},
     //
     {"WronglyAttached/FwdTracksAll", "; N tracks attached to min-bias;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"WronglyAttached/FwdTracksUnambig", "; N tracks attached to min-bias;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"WronglyAttached/MIDTracksAll", "; N tracks attached to min-bias;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"WronglyAttached/MIDTracksUnambig", "; N tracks attached to min-bias;", {HistType::kTH1I, {{100, 0, 100}}}},
     {"WronglyAttached/MIDTracksUnambigDeltaT", "; Time difference (track time - col. time), ns;", {HistType::kTH1D, {{1000, -500., 500.}}}},
     {"WronglyAttached/ColTimeResolution", ";Time resolution, ns;", {HistType::kTH1D, {{500, 0., 500.}}}},
     //
     {"Timings/TrTimeRes", "; Time Resolustion, ns;", {HistType::kTH1D, {{10000, 0., 1000.}}}},
     {"Timings/GlobalMuonTimeRes", "; Time Resolustion, ns;", {HistType::kTH1D, {{10000, 0., 1000.}}}},
     {"Timings/GlobalForwardTimeRes", "; Time Resolustion, ns;", {HistType::kTH1D, {{10000, 0., 1000.}}}},
     {"Timings/AmbMuonStandaloneTrackTime", "; Time, ns;", {HistType::kTH1D, {{5000, -50., 50.}}}},
     //
     {"BCs/AmbFwdTracks", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     {"BCs/FT0Signals", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     {"BCs/SliceSize", ";BC slice size for amb. tracks;", {HistType::kTH1I, {{100, 0, 100}}}},
     //
     {"TracksWithFT0/Eta", ";#eta;", {HistType::kTH1D, {{80, -4., 4.}}}},
     {"TracksWithFT0/TimeFT0A", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     {"TracksWithFT0/TimeFT0C", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     //
     {"PairSelection/MassRaw", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtRaw", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassMC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/Pt", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0MC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtNoFT0", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/NTracksPerBC", ";N tracks per BC;", {HistType::kTH1I, {{50, 0, 50}}}},
     //
     {"PairSelection/TracksChi2", ";chi2;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/TracksChi2Raw", ";chi2;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/NMatchesExact", ";N pairs;", {HistType::kTH1I, {{2, 0, 2}}}},
     {"PairSelection/NMatchesTrueExact", ";N true pairs;", {HistType::kTH1I, {{2, 0, 2}}}},
     {"PairSelection/NMatchesTrueSigExact", ";N true signal pairs;", {HistType::kTH1I, {{2, 0, 2}}}},
     {"PairSelection/NMatchesFakeExact", ";N fakes;", {HistType::kTH1I, {{2, 0, 2}}}}}};

  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;
  // using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  void init(InitContext&)
  {
    for (int i = 0; i < 5; i++) {
      registry.get<TH1>(HIST("TrackTypes/TrackTypesAll"))->GetXaxis()->SetBinLabel(i + 1, fwdTrackLabels[i]);
      registry.get<TH1>(HIST("TrackTypes/AmbTrackTypesAll"))->GetXaxis()->SetBinLabel(i + 1, fwdTrackLabels[i]);
      registry.get<TH1>(HIST("TrackTypes/TrackTypesSig"))->GetXaxis()->SetBinLabel(i + 1, fwdTrackLabels[i]);
      registry.get<TH1>(HIST("TrackTypes/AmbTrackTypesSig"))->GetXaxis()->SetBinLabel(i + 1, fwdTrackLabels[i]);
    }
  }

  void process(Collisions const& collisions,
               aod::McCollisions const& mcCollisions,
               aod::McParticles_001 const& mcParticles,
               Tracks const& fwdTracks,
               aod::AmbiguousFwdTracks const& ambFwdTracks,
               aod::BCs const& bcs,
               aod::FT0s const& ft0s)
  {
    using namespace o2::aod::fwdtrack;

    // bcs with ft0 signals
    std::map<uint64_t, std::pair<int64_t, int64_t>> BCsWithFT0;
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      registry.fill(HIST("BCs/FT0Signals"), bc);
      BCsWithFT0[bc] = std::make_pair(ft0.bc().globalIndex(), ft0.globalIndex());
    }

    // vector for duplicate checks
    std::vector<std::pair<int, int>> mcPartIdsToNFwdTracks;            // pairs <mcParticleId, counter of trackIDs>
    std::vector<std::pair<int, int>> mcPartIdsToNGlobalMuonTracks;     // pairs <mcParticleId, counter of trackIDs>
    std::vector<std::pair<int, int>> mcPartIdsToNMuonStandaloneTracks; // pairs <mcParticleId, counter of trackIDs>

    for (const auto& tr : fwdTracks) {
      registry.fill(HIST("TrackTypes/TrackTypesAll"), tr.trackType());
      if (tr.has_mcParticle()) {
        if (tr.mcParticle().mcCollision().generatorsID() == signalSourceID)
          registry.fill(HIST("TrackTypes/TrackTypesSig"), tr.trackType());
      }
      //
      registry.fill(HIST("Timings/TrTimeRes"), tr.trackTimeRes());
      if (tr.trackType() == ForwardTrackTypeEnum::GlobalMuonTrack)
        registry.fill(HIST("Timings/GlobalMuonTimeRes"), tr.trackTimeRes());
      if (tr.trackType() == ForwardTrackTypeEnum::GlobalForwardTrack)
        registry.fill(HIST("Timings/GlobalForwardTimeRes"), tr.trackTimeRes());
      //
      // n fwd tracks per mc particle
      auto mcPartId = tr.mcParticleId();
      auto it = std::find_if(mcPartIdsToNFwdTracks.begin(), mcPartIdsToNFwdTracks.end(),
                             [=](const std::pair<int, int>& element) { return element.first == mcPartId; });
      if (it != mcPartIdsToNFwdTracks.end())
        it->second++; // add 1 to counter corresponding to `mcPartId`
      else
        mcPartIdsToNFwdTracks.emplace_back(mcPartId, 1);
      //
      // n GlobalMuonTracks per mc particle
      if (tr.trackType() == ForwardTrackTypeEnum::GlobalMuonTrack) {
        auto it = std::find_if(mcPartIdsToNGlobalMuonTracks.begin(), mcPartIdsToNGlobalMuonTracks.end(),
                               [=](const std::pair<int, int>& element) { return element.first == mcPartId; });
        if (it != mcPartIdsToNGlobalMuonTracks.end())
          it->second++; // add 1 to counter corresponding to `mcPartId`
        else
          mcPartIdsToNGlobalMuonTracks.emplace_back(mcPartId, 1);
      }
      //
      // n MuonStandaloneTracks per mc particle
      if (tr.trackType() == ForwardTrackTypeEnum::MuonStandaloneTrack) {
        auto it = std::find_if(mcPartIdsToNMuonStandaloneTracks.begin(), mcPartIdsToNMuonStandaloneTracks.end(),
                               [=](const std::pair<int, int>& element) { return element.first == mcPartId; });
        if (it != mcPartIdsToNMuonStandaloneTracks.end())
          it->second++; // add 1 to counter corresponding to `mcPartId`
        else
          mcPartIdsToNMuonStandaloneTracks.emplace_back(mcPartId, 1);
      }
    }

    // counters of tracks per mc particle
    for (const auto& item : mcPartIdsToNFwdTracks) {
      registry.fill(HIST("TrackTypes/NTracksPerMcPart"), item.second);
    }
    for (const auto& item : mcPartIdsToNGlobalMuonTracks) {
      registry.fill(HIST("TrackTypes/NGlobalMuonTracksPerMcPart"), item.second);
    }
    for (const auto& item : mcPartIdsToNMuonStandaloneTracks) {
      registry.fill(HIST("TrackTypes/NMuonStandaloneTracksPerMcPart"), item.second);
    }

    // vectors for pair matching
    std::vector<std::pair<uint64_t, std::vector<int>>> ambTracksToGlobalBCs; // pairs <globalBC, vector of trackIDs>
    std::vector<std::tuple<int, int, int64_t>> ambTrackMatchesExact;         // tuples <trackID, trackID, globalBC> for exactly 2 tracks per BC
    // ids of ambiguous forward tracks
    std::vector<int> ambigFwdTracks;

    // collect (ambiguous) trackIds belonging to same BCs
    for (const auto& ambFwdTr : ambFwdTracks) {
      auto trId = ambFwdTr.fwdtrackId();
      ambigFwdTracks.emplace_back(trId);
      const auto& tr = fwdTracks.iteratorAt(trId);
      auto trType = tr.trackType();
      registry.fill(HIST("TrackTypes/AmbTrackTypesAll"), trType);
      if (tr.has_mcParticle()) {
        if (tr.mcParticle().mcCollision().generatorsID() == signalSourceID)
          registry.fill(HIST("TrackTypes/AmbTrackTypesSig"), tr.trackType());
      }
      //
      if (filterFwdType != -1 && trType != filterFwdType) // if needed, filter tracks used for matching by type
        continue;
      //
      const auto& bcSlice = ambFwdTr.bc();
      auto sliceSize = bcSlice.size();
      auto first = bcSlice.begin();
      auto firstBC = first.globalBC();
      auto bc = relativeTime_to_LocalBC(firstBC * o2::constants::lhc::LHCBunchSpacingNS + (double)tr.trackTime());
      //
      registry.fill(HIST("Timings/AmbMuonStandaloneTrackTime"), tr.trackTime());
      registry.fill(HIST("BCs/SliceSize"), sliceSize);
      registry.fill(HIST("BCs/AmbFwdTracks"), bc);
      //
      auto it = std::find_if(ambTracksToGlobalBCs.begin(), ambTracksToGlobalBCs.end(),
                             [=](const std::pair<uint64_t, std::vector<int>>& element) { return element.first == bc; });
      if (it != ambTracksToGlobalBCs.end())
        it->second.emplace_back(ambFwdTr.fwdtrackId()); // add Id to vector corresponding to `bc`
      else
        ambTracksToGlobalBCs.emplace_back(bc, std::vector(1, ambFwdTr.fwdtrackId())); // add vector corresponding to `bc` with 1 element (trId)
    }

    // collecting matches
    // considering exactly same BCs
    for (const auto& item : ambTracksToGlobalBCs) {
      auto bc = item.first;
      auto trackIDs = item.second;
      auto nTracks = trackIDs.size();
      registry.fill(HIST("PairSelection/NTracksPerBC"), nTracks);
      if (nTracks == 2)
        ambTrackMatchesExact.emplace_back(trackIDs[0], trackIDs[1], bc);
    }

    // check "exact" matches
    for (const auto& item : ambTrackMatchesExact) {
      const auto& tr1 = fwdTracks.iteratorAt(std::get<0>(item));
      const auto& tr2 = fwdTracks.iteratorAt(std::get<1>(item));
      // has ft0 signal?
      auto bc = std::get<2>(item);
      auto bcFT0 = BCsWithFT0.find(bc);
      bool hasFT0 = bcFT0 != BCsWithFT0.end();
      if (hasFT0) {
        registry.fill(HIST("TracksWithFT0/Eta"), tr1.eta());
        registry.fill(HIST("TracksWithFT0/Eta"), tr2.eta());
        const auto& ft0 = ft0s.iteratorAt(bcFT0->second.second);
        registry.fill(HIST("TracksWithFT0/TimeFT0A"), ft0.timeA());
        registry.fill(HIST("TracksWithFT0/TimeFT0C"), ft0.timeC());
      }
      //
      int mcPartId1 = tr1.mcParticleId();
      int mcPartId2 = tr2.mcParticleId();
      if (mcPartId1 < 0 || mcPartId2 < 0)
        continue;
      const auto& mcPart1 = mcParticles.iteratorAt(mcPartId1);
      const auto& mcPart2 = mcParticles.iteratorAt(mcPartId2);
      //
      int mcId1 = tr1.mcParticle().mcCollisionId();
      int mcId2 = tr2.mcParticle().mcCollisionId();
      int genId1 = tr1.mcParticle().mcCollision().generatorsID();
      int genId2 = tr2.mcParticle().mcCollision().generatorsID();
      //
      registry.fill(HIST("PairSelection/NMatchesExact"), 1);
      if (mcId1 == mcId2)
        registry.fill(HIST("PairSelection/NMatchesTrueExact"), 1);
      if (mcId1 != mcId2)
        registry.fill(HIST("PairSelection/NMatchesFakeExact"), 1);
      if (mcId1 == mcId2 && genId1 == genId2 && genId1 == signalSourceID)
        registry.fill(HIST("PairSelection/NMatchesTrueSigExact"), 1);
      //
      TLorentzVector mcP1, mcP2, pPairMC;
      mcP1.SetXYZM(mcPart1.px(), mcPart1.py(), mcPart1.pz(), mmuon);
      mcP2.SetXYZM(mcPart2.px(), mcPart2.py(), mcPart2.pz(), mmuon);
      pPairMC = mcP1 + mcP2;
      //
      TLorentzVector p1, p2, pPair;
      p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mmuon);
      p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mmuon);
      pPair = p1 + p2;
      registry.fill(HIST("PairSelection/MassRaw"), pPair.M());
      registry.fill(HIST("PairSelection/PtRaw"), pPair.Pt());
      registry.fill(HIST("PairSelection/TracksChi2Raw"), tr1.chi2());
      registry.fill(HIST("PairSelection/TracksChi2Raw"), tr2.chi2());
      // some simple selection
      if (tr1.sign() * tr2.sign() >= 0)
        continue;
      if (pPair.Pt() > 0.25) {
        // comment out for debugging
        // LOGF(info, "mcPartId1 = %d", mcPartId1);
        // p1.Print();
        // mcP1.Print();
        // LOGF(info, "mcPartId2 = %d", mcPartId2);
        // p2.Print();
        // mcP2.Print();
        continue;
      }
      //
      registry.fill(HIST("PairSelection/Mass"), pPair.M());
      registry.fill(HIST("PairSelection/MassMC"), pPairMC.M());
      registry.fill(HIST("PairSelection/Pt"), pPair.Pt());
      registry.fill(HIST("PairSelection/TracksChi2"), tr1.chi2());
      registry.fill(HIST("PairSelection/TracksChi2"), tr2.chi2());
      //
      bool passFT0 = false;
      if (hasFT0) {
        const auto& ft0 = ft0s.iteratorAt(bcFT0->second.second);
        passFT0 = ft0.timeC() > -1. && ft0.timeC() < 1.;
      } else {
        passFT0 = true;
      }
      if (passFT0) {
        registry.fill(HIST("PairSelection/MassNoFT0"), pPair.M());
        registry.fill(HIST("PairSelection/MassNoFT0MC"), pPairMC.M());
        registry.fill(HIST("PairSelection/PtNoFT0"), pPair.Pt());
      }
    }

    // signal tracks wrongly attached to hadronic collisions
    for (const auto& col : collisions) {
      auto colT = (double)col.collisionTime(); // relative time
      auto colTRes = (double)col.collisionTimeRes();
      int genId = col.mcCollision().generatorsID();
      if (genId != signalSourceID) {
        int nSignalsAll = 0;
        int nSignalsUnambig = 0;
        int nSignalsMIDAll = 0;
        int nSignalsMIDUnambig = 0;
        auto groupedTracks = fwdTracks.sliceBy(aod::track::collisionId, col.globalIndex());
        for (const auto& tr : groupedTracks) {
          if (!tr.has_mcParticle())
            continue;
          int src = tr.mcParticle().mcCollision().generatorsID();
          if (src != signalSourceID)
            continue;
          auto it = std::find(ambigFwdTracks.begin(), ambigFwdTracks.end(), tr.globalIndex());
          bool isAmbiguous = it != ambigFwdTracks.end();
          nSignalsAll++;
          if (!isAmbiguous)
            nSignalsUnambig++;
          if (tr.trackType() == ForwardTrackTypeEnum::MuonStandaloneTrack)
            nSignalsMIDAll++;
          if (!isAmbiguous && tr.trackType() == ForwardTrackTypeEnum::MuonStandaloneTrack) {
            nSignalsMIDUnambig++;
            auto trackT = (double)tr.trackTime();
            double deltaT = trackT - colT;
            registry.fill(HIST("WronglyAttached/MIDTracksUnambigDeltaT"), deltaT);
          }
        }
        if (nSignalsAll != 0)
          registry.fill(HIST("WronglyAttached/FwdTracksAll"), nSignalsAll);
        if (nSignalsMIDAll != 0)
          registry.fill(HIST("WronglyAttached/MIDTracksAll"), nSignalsMIDAll);
        if (nSignalsUnambig != 0)
          registry.fill(HIST("WronglyAttached/FwdTracksUnambig"), nSignalsUnambig);
        if (nSignalsMIDUnambig != 0) {
          registry.fill(HIST("WronglyAttached/MIDTracksUnambig"), nSignalsMIDUnambig);
          registry.fill(HIST("WronglyAttached/ColTimeResolution"), colTRes);
        }
      }
    }

    BCsWithFT0.clear();
    ambTracksToGlobalBCs.clear();
    ambTrackMatchesExact.clear();
    mcPartIdsToNFwdTracks.clear();
    mcPartIdsToNGlobalMuonTracks.clear();
    mcPartIdsToNMuonStandaloneTracks.clear();
    ambigFwdTracks.clear();
  }
};

struct AnalysisCollisions {
  HistogramRegistry registry{
    "registry",
    {{"BCs/McCollisionsSig", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     {"BCs/McCollisionsBkg", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     //
     {"Timings/Resolution", ";Time resolution, ns;", {HistType::kTH1D, {{500, 0., 500.}}}}}};

  using Collisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
  using Tracks = soa::Join<aod::FwdTracks, aod::FwdTracksCov, aod::McFwdTrackLabels>;

  void process(Collisions const& collisions,
               aod::McCollisions const& mcCollsions,
               aod::BCs const& bcs,
               aod::FT0s const& ft0s)
  {
    // time resolution
    for (const auto& col : collisions) {
      auto timeRes = col.collisionTimeRes();
      registry.fill(HIST("Timings/Resolution"), timeRes);
    }

    // mc collisions bcs
    for (const auto& mcCol : mcCollsions) {
      if (mcCol.generatorsID() == 1)
        registry.fill(HIST("BCs/McCollisionsSig"), mcCol.bc().globalBC());
      if (mcCol.generatorsID() == 0)
        registry.fill(HIST("BCs/McCollisionsBkg"), mcCol.bc().globalBC());
    }
  }
};

struct FwdDimuonsAnalysis {
  HistogramRegistry registry{
    "registry",
    {{"BCs/MuonTrackPairs", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     {"BCs/FT0Signals", ";BC;", {HistType::kTH1D, {{1000000, 0., 1000000.}}}},
     //
     {"TracksWithFT0/Eta", ";#eta;", {HistType::kTH1D, {{80, -4., 4.}}}},
     {"TracksWithFT0/TimeFT0A", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     {"TracksWithFT0/TimeFT0C", ";time A, ns;", {HistType::kTH1D, {{1000, -50., 50.}}}},
     //
     {"PairSelection/MassRaw", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtRaw", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/Mass", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassMC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/Pt", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/MassNoFT0", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/MassNoFT0MC", ";#it{m}_{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/PtNoFT0", ";#it{p}_{T}^{#mu#mu}, GeV;", {HistType::kTH1D, {{100, 0., 10.}}}},
     //
     {"PairSelection/TracksChi2Raw", ";Chi2;", {HistType::kTH1D, {{100, 0., 10.}}}},
     {"PairSelection/TracksChi2", ";Chi2;", {HistType::kTH1D, {{100, 0., 10.}}}}}};

  void process(o2::aod::EventCandidates const& eventCandidates,
               soa::Join<o2::aod::SkimmedMuons, o2::aod::SkimmedMuonsExtra, o2::aod::SkimmedMuonTrackLabels> const& muonTracks,
               o2::aod::SkimmedMCParticles const& mcParticles,
               o2::aod::BCs const& bcs,
               o2::aod::FT0s const& ft0s)
  {
//    // bcs with ft0 signals
//    std::map<uint64_t, std::pair<int64_t, int64_t>> BCsWithFT0;
//    for (const auto& ft0 : ft0s) {
//      uint64_t bc = ft0.bc().globalBC();
//      registry.fill(HIST("BCs/FT0Signals"), bc);
//      BCsWithFT0[bc] = std::make_pair(ft0.bc().globalIndex(), ft0.globalIndex());
//    }
//
//    LOGF(info, "N candidates = %d", eventCandidates.size());
//
//    for (const auto& cand : eventCandidates) {
//      const auto& trackIDs = cand.matchedFwdTracksIds();
//      const auto& tr1 = muonTracks.iteratorAt(trackIDs[0]);
//      const auto& tr2 = muonTracks.iteratorAt(trackIDs[1]);
//      // has ft0 signal?
//      int64_t bc = cand.globalBC();
//      registry.fill(HIST("BCs/MuonTrackPairs"), bc);
//      auto bcFT0 = BCsWithFT0.find(bc);
//      bool hasFT0 = bcFT0 != BCsWithFT0.end();
//      if (hasFT0) {
//        registry.fill(HIST("TracksWithFT0/Eta"), tr1.eta());
//        registry.fill(HIST("TracksWithFT0/Eta"), tr2.eta());
//        const auto& ft0 = ft0s.iteratorAt(bcFT0->second.second);
//        registry.fill(HIST("TracksWithFT0/TimeFT0A"), ft0.timeA());
//        registry.fill(HIST("TracksWithFT0/TimeFT0C"), ft0.timeC());
//      }
//      //
//      int mcPartId1 = tr1.mcParticleId();
//      int mcPartId2 = tr2.mcParticleId();
//      const auto& mcPart1 = mcParticles.iteratorAt(mcPartId1);
//      const auto& mcPart2 = mcParticles.iteratorAt(mcPartId2);
//      TLorentzVector mcP1, mcP2, pPairMC;
//      mcP1.SetXYZM(mcPart1.px(), mcPart1.py(), mcPart1.pz(), mmuon);
//      mcP2.SetXYZM(mcPart2.px(), mcPart2.py(), mcPart2.pz(), mmuon);
//      pPairMC = mcP1 + mcP2;
//      //
//      TLorentzVector p1, p2, pPair;
//      p1.SetXYZM(tr1.px(), tr1.py(), tr1.pz(), mmuon);
//      p2.SetXYZM(tr2.px(), tr2.py(), tr2.pz(), mmuon);
//      pPair = p1 + p2;
//      registry.fill(HIST("PairSelection/MassRaw"), pPair.M());
//      registry.fill(HIST("PairSelection/PtRaw"), pPair.Pt());
//      registry.fill(HIST("PairSelection/TracksChi2Raw"), tr1.chi2());
//      registry.fill(HIST("PairSelection/TracksChi2Raw"), tr2.chi2());
//      // some simple selection
//      if (tr1.sign() * tr2.sign() >= 0)
//        continue;
//      if (pPair.Pt() > 0.25) {
//        continue;
//      }
//      //
//      registry.fill(HIST("PairSelection/Mass"), pPair.M());
//      registry.fill(HIST("PairSelection/MassMC"), pPairMC.M());
//      registry.fill(HIST("PairSelection/Pt"), pPair.Pt());
//      registry.fill(HIST("PairSelection/TracksChi2"), tr1.chi2());
//      registry.fill(HIST("PairSelection/TracksChi2"), tr2.chi2());
//      //
//      bool passFT0 = false;
//      if (hasFT0) {
//        const auto& ft0 = ft0s.iteratorAt(bcFT0->second.second);
//        passFT0 = ft0.timeC() > -1. && ft0.timeC() < 1.;
//      } else {
//        passFT0 = true;
//      }
//      if (passFT0) {
//        registry.fill(HIST("PairSelection/MassNoFT0"), pPair.M());
//        registry.fill(HIST("PairSelection/MassNoFT0MC"), pPairMC.M());
//        registry.fill(HIST("PairSelection/PtNoFT0"), pPair.Pt());
//      }
//    }
//    BCsWithFT0.clear();
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec
  {
    // debug
    // adaptAnalysisTask<AnalysisTracksMc>(cfgc),
    // adaptAnalysisTask<AnalysisTracks>(cfgc),
    // adaptAnalysisTask<AnalysisCollisions>(cfgc),
    // release?
    adaptAnalysisTask<FwdDimuonsAnalysis>(cfgc)
  };
}
