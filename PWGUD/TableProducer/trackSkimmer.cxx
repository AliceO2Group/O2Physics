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
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackSkimmer {
  int32_t signalGenID{1};

  Produces<o2::aod::SkimmedMuons> muonTracks;
  Produces<o2::aod::SkimmedMuonsExtra> muonsExtra;
  Produces<o2::aod::SkimmedMuonsCov> muonsCov;
  Produces<o2::aod::SkimmedMuonTrackLabels> muonLabels;
  Produces<o2::aod::SkimmedMCEvents> skMCEvents;
  Produces<o2::aod::SkimmedMCParticles> skMCParticles;

  // cuts for forward tracks
  Configurable<int> fUseFwdCuts{"useFwdCuts", 1, "Use cuts for forward tracks"};
  Configurable<float> fFwdPtLow{"fwdPtLow", 0.5, "Minimal Pt for forward tracks"};
  Configurable<float> fFwdPtHigh{"fwdPtHigh", 4., "Maximal Pt for forward tracks"};
  Configurable<float> fFwdEtaLow{"fwdEtaLow", -4.0, "Maximal Eta for forward tracks"};
  Configurable<float> fFwdEtaHigh{"fwdEtaHigh", -2.5, "Maximal Eta for forward tracks"};
  Configurable<float> fMuonRAtAbsorberEndLow{"muonRAtAbsorberEndLow", 17.6, "Minimal muon R at absorber end"};
  Configurable<float> fMuonRAtAbsorberEndHigh{"muonRAtAbsorberEndHigh", 89.5, "Maximal muon R at absorber end"};
  Configurable<float> fMuonPDcaHighFirst{"fMuonPDcaHighFirst", 594.0, "Primary PDCA cut: Maximal value for R < 26.5"};
  Configurable<float> fMuonPDcaHighSecond{"fMuonPDcaHighSecond", 324.0, "Additional PDCA cut: Maximal value for R >= 26.5"};
  Configurable<float> fChi2Low{"fwdChi2Low", 0.0, "Minimal Chi2 for forward tracks"};
  Configurable<float> fChi2High{"fwdChi2High", 10000.0, "Maximal Chi2 for forward tracks"};

  template <typename TFwdTrack>
  bool applyFwdCuts(TFwdTrack const& track)
  {
    // using any cuts at all?
    if (!fUseFwdCuts) {
      return true;
    }
    // check Pt cuts
    float pt = track.pt();
    bool checkPt = pt > fFwdPtLow && pt < fFwdPtHigh;
    if (!checkPt) {
      return false;
    }
    // check pseudorapidity cuts
    float eta = track.eta();
    bool checkEta = eta > fFwdEtaLow && eta < fFwdEtaHigh;
    if (!checkEta) {
      return false;
    }
    // check muon R
    float r = track.rAtAbsorberEnd();
    bool checkR = r > fMuonRAtAbsorberEndLow && r < fMuonRAtAbsorberEndHigh;
    if (!checkR) {
      return false;
    }
    // check pDCA
    float pDCA = track.pDca();
    bool checkPDCA = r < 26.5 ? pDCA < fMuonPDcaHighFirst : pDCA < fMuonPDcaHighSecond;
    if (!checkPDCA) {
      return false;
    }
    // check Chi2
    float chi2 = track.chi2();
    bool checkChi2 = chi2 > fChi2Low && chi2 < fChi2High;
    if (!checkChi2) {
      return false;
    }
    return true;
  }

  template <uint32_t TTrackType, uint32_t TUseMC, typename TFwdTracks, typename TAmbFwdTracks, typename TBCs,
            typename TMcFwdTrackLabels, typename TMcCollisions, typename TMcParticles>
  void skimFwdTracks(TFwdTracks const& tracks,
                     TAmbFwdTracks const& ambTracks,
                     TBCs const& bcs,
                     TMcFwdTrackLabels* mcFwdTrackLabels,
                     TMcCollisions* mcCollisions,
                     TMcParticles* mcParticles)
  {
    std::vector<int32_t> newPartIDs;
    std::vector<int32_t> newEventIDs;

    if constexpr (static_cast<bool>(TUseMC)) {
      newPartIDs.resize(mcParticles->size(), -1);
      newEventIDs.resize(mcCollisions->size(), -1);
      int32_t newPartID = 0;
      int32_t newEventID = 0;
      for (const auto& label : *mcFwdTrackLabels) {
        int32_t mcPartID = label.mcParticleId();
        const auto& mcPart = mcParticles->iteratorAt(mcPartID);
        int32_t mcEventID = mcPart.mcCollisionId();
        const auto& mcEvent = mcCollisions->iteratorAt(mcEventID);
        bool isSignal = mcEvent.generatorsID() == signalGenID;
        if (isSignal) {
          newPartIDs[mcPartID] = newPartID;
          newPartID++;
          if (newEventIDs[mcEventID] == -1) {
            newEventIDs[mcEventID] = newEventID;
            newEventID++;
          }
        }
      }

      // storing MC particles
      for (int32_t i = 0; i < mcParticles->size(); i++) {
        if (newPartIDs[i] == -1) {
          continue;
        }
        const auto& mcPart = mcParticles->iteratorAt(i);
        int32_t mcEventID = mcPart.mcCollisionId();
        int32_t newEventID = newEventIDs[mcEventID];
        // collecting new mother IDs
        const auto& motherIDs = mcPart.mothersIds();
        std::vector<int32_t> newMotherIDs;
        newMotherIDs.reserve(motherIDs.size());
        for (auto motherID : motherIDs) {
          newMotherIDs.push_back(newPartIDs[motherID]);
        }
        // collecting new daughter IDs
        const auto& daughterIDs = mcPart.daughtersIds();
        int32_t newDaughterIDs[2] = {-1, -1};
        if (daughterIDs.size() > 0) {
          newDaughterIDs[0] = daughterIDs.front();
          newDaughterIDs[1] = daughterIDs.back();
        }
        skMCParticles(newEventID, mcPart.pdgCode(), mcPart.statusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                      mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
      }

      // storing MC events
      for (int32_t i = 0; i < mcCollisions->size(); i++) {
        if (newEventIDs[i] == -1) {
          continue;
        }
        const auto& mcEvent = mcCollisions->iteratorAt(i);
        skMCEvents(mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                   mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
      }

      newEventIDs.clear();
    }

    for (const auto& ambTr : ambTracks) {
      auto trId = ambTr.fwdtrackId();
      const auto& tr = tracks.iteratorAt(trId);
      // filter only interesting type of tracks if needed
      if constexpr (TTrackType != -1) {
        auto trType = tr.trackType();
        if (trType != TTrackType) {
          continue;
        }
      }
      // skip if doesn't pass cuts
      if (!applyFwdCuts(tr)) {
        continue;
      }
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      uint64_t firstBC = first.globalBC();
      double trTime = firstBC * o2::constants::lhc::LHCBunchSpacingNS + (double)tr.trackTime(); // using absolute time
      muonTracks(tr.px(), tr.py(), tr.pz(), tr.sign(), trTime, tr.trackTimeRes());
      muonsCov(tr.x(), tr.y(), tr.z(), tr.tgl(), tr.signed1Pt(), tr.cXX(), tr.cXY(),
               tr.cYY(), tr.cPhiX(), tr.cPhiY(), tr.cPhiPhi(), tr.cTglX(), tr.cTglY(),
               tr.cTglPhi(), tr.cTglTgl(), tr.c1PtX(), tr.c1PtY(), tr.c1PtPhi(), tr.c1PtTgl(),
               tr.c1Pt21Pt2());
      muonsExtra(tr.nClusters(), tr.pDca(), tr.rAtAbsorberEnd(), tr.chi2(), tr.chi2MatchMCHMID(),
                 tr.mchBitMap(), tr.midBitMap(), tr.midBoards());
      if constexpr (static_cast<bool>(TUseMC)) {
        const auto& label = mcFwdTrackLabels->iteratorAt(trId);
        int32_t mcPartID = label.mcParticleId();
        uint16_t mcMask = label.mcMask();
        int32_t newPartID = newPartIDs[mcPartID];
        muonLabels(newPartID, mcMask);
      }
    }

    newPartIDs.clear();
  }

  void processFwdMC(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& tracks,
                    o2::aod::AmbiguousFwdTracks const& ambTracks,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
                    o2::aod::McCollisions const& mcCollisions,
                    o2::aod::McParticles const& mcParticles,
                    o2::aod::BCs const& bcs)
  {
    using namespace o2::aod::fwdtrack;
    const uint32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    const uint32_t useMC = 1u;
    skimFwdTracks<trType, useMC>(tracks, ambTracks, bcs, &mcFwdTrackLabels, &mcCollisions, &mcParticles);
  }

  void processFwd(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& tracks,
                  o2::aod::AmbiguousFwdTracks const& ambTracks,
                  o2::aod::BCs const& bcs)
  {
    using namespace o2::aod::fwdtrack;
    const uint32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    const uint32_t useMC = 0u;
    skimFwdTracks<trType, useMC>(tracks, ambTracks, bcs,
                                 (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McCollisions*)nullptr, (o2::aod::McParticles_001*)nullptr);
  }

  PROCESS_SWITCH(TrackSkimmer, processFwdMC, "Produce only muon tracks with MC information", false);
  PROCESS_SWITCH(TrackSkimmer, processFwd, "Produce only muon tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<TrackSkimmer>(cfgc)};
}