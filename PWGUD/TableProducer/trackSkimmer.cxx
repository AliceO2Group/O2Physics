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
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGUD/DataModel/UDTables.h"

#include "TH1I.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct TrackSkimmer {
  int32_t fSignalGenID{1};
  bool fDoMC;

  OutputObj<TH1I> hNSig{TH1I("hNSig", ";n tracks;", 2, 0, 2)};
  OutputObj<TH1I> hNSigITSTPC{TH1I("hNSigITSTPC", ";n tracks;", 2, 0, 2)};

  Produces<o2::aod::SkimmedMCEvents> skMCEvents;
  Produces<o2::aod::SkimmedMCParticles> skMCParticles;

  Produces<o2::aod::SkimmedMuons> muonTracks;
  Produces<o2::aod::SkimmedMuonsExtra> muonsExtra;
  Produces<o2::aod::SkimmedMuonsCov> muonsCov;
  Produces<o2::aod::SkimmedMuonTrackLabels> muonLabels;

  Produces<o2::aod::SkimmedBarTracks> barTracks;
  Produces<o2::aod::SkimmedBarTracksCov> barsCov;
  Produces<o2::aod::SkimmedBarTracksExtra> barsExtra;
  Produces<o2::aod::SkimmedBarTrackLabels> barLabels;

  // cuts for forward tracks
  Configurable<int> fUseFwdCuts{"useFwdCuts", 1, "Use cuts for forward tracks"};
  // basic
  Configurable<float> fFwdPtLow{"fwdPtLow", 0.5, "Minimal Pt for forward tracks"};
  Configurable<float> fFwdPtHigh{"fwdPtHigh", 4., "Maximal Pt for forward tracks"};
  Configurable<float> fFwdEtaLow{"fwdEtaLow", -4.0, "Maximal Eta for forward tracks"};
  Configurable<float> fFwdEtaHigh{"fwdEtaHigh", -2.5, "Maximal Eta for forward tracks"};
  // quality
  Configurable<float> fMuonRAtAbsorberEndLow{"muonRAtAbsorberEndLow", 17.6, "Minimal muon R at absorber end"};
  Configurable<float> fMuonRAtAbsorberEndHigh{"muonRAtAbsorberEndHigh", 89.5, "Maximal muon R at absorber end"};
  Configurable<float> fMuonPDcaHighFirst{"fMuonPDcaHighFirst", 594.0, "Primary PDCA cut: Maximal value for R < 26.5"};
  Configurable<float> fMuonPDcaHighSecond{"fMuonPDcaHighSecond", 324.0, "Additional PDCA cut: Maximal value for R >= 26.5"};
  Configurable<float> fFwdChi2Low{"fwdChi2Low", 0.0, "Minimal Chi2 for forward tracks"};
  Configurable<float> fFwdChi2High{"fwdChi2High", 10000.0, "Maximal Chi2 for forward tracks"};

  // cuts for central-barrel tracks
  Configurable<int> fUseBarCuts{"useBarCuts", 1, "Use cuts for barrel tracks"};
  // basic
  Configurable<float> fBarPtLow{"barPtLow", 0., "Minimal Pt for barrel tracks"};
  Configurable<float> fBarPtHigh{"barPtHigh", 1000., "Maximal Pt for barrel tracks"};
  Configurable<float> fBarEtaLow{"barEtaLow", -0.9, "Maximal Eta for barrel tracks"};
  Configurable<float> fBarEtaHigh{"barEtaHigh", 0.9, "Maximal Eta for barrel tracks"};
  // quality: ITS
  Configurable<int> fITSNClusLow{"ITSNClusLow", 4, "Minimal number of ITS clusters"};
  Configurable<int> fITSNClusHigh{"ITSNClusHigh", 9, "Maximal number of ITS clusters"};
  Configurable<float> fITSChi2Low{"ITSChi2Low", 0., "Minimal Chi2 in ITS per cluster"};
  Configurable<float> fITSChi2High{"ITSChi2High", 5., "Maximal Chi2 in ITS per cluster"};
  // quality: TPC
  Configurable<int> fTPCNClusCRLow{"TPCNClusCRLow", 70, "Minimal number of TPC clusters (crossed rows)"};
  Configurable<int> fTPCNClusCRHigh{"TPCNClusCRHigh", 161, "Maximal number of TPC clusters (crossed rows)"};
  Configurable<float> fTPCChi2Low{"TPCChi2Low", 0., "Minimal Chi2 in TPC per cluster"};
  Configurable<float> fTPCChi2High{"TPCChi2High", 4., "Maximal Chi2 in TPC per cluster"};
  // todo: quality: DCA
  // Configurable<int> fCheckMaxDcaXY{"checkMaxDCACut", 1, "Apply cut on maximal DCA_xy"};
  // Configurable<float> fDcaZLow{"dcaZLow", -3., "Minimal DCA_z for barrel tracks"};
  // Configurable<float> fDcaZHigh{"dcaZHigh", 3., "Maximal DCA_z for barrel tracks"};

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
    bool checkChi2 = chi2 > fFwdChi2Low && chi2 < fFwdChi2High;
    if (!checkChi2) {
      return false;
    }
    return true;
  }

  template <typename TBarTrack>
  bool applyBarCuts(TBarTrack const& track)
  {
    // using any cuts at all?
    if (!fUseBarCuts) {
      return true;
    }
    // check Pt cuts
    float pt = track.pt();
    bool checkPt = pt > fBarPtLow && pt < fBarPtHigh;
    if (!checkPt) {
      return false;
    }
    // check pseudorapidity cuts
    float eta = track.eta();
    bool checkEta = eta > fBarEtaLow && eta < fBarEtaHigh;
    if (!checkEta) {
      return false;
    }
    // check ITS cuts
    bool checkITSNClus = track.itsNCls() >= static_cast<uint8_t>(fITSNClusLow) && track.itsNCls() <= static_cast<uint8_t>(fITSNClusHigh);
    if (!checkITSNClus) {
      return false;
    }
    bool checkChi2ITS = track.itsChi2NCl() > fITSChi2Low && track.itsChi2NCl() < fITSChi2High;
    if (!checkChi2ITS) {
      return false;
    }
    // check TPC cuts
    int checkTPCNClus = track.tpcNClsCrossedRows() > static_cast<int16_t>(fTPCNClusCRLow) && track.tpcNClsCrossedRows() < static_cast<int16_t>(fTPCNClusCRHigh);
    if (!checkTPCNClus) {
      return false;
    }
    bool checkChi2TPC = track.tpcChi2NCl() > fTPCChi2Low && track.tpcChi2NCl() < fTPCChi2High;
    if (!checkChi2TPC) {
      return false;
    }
    // todo: check DCA
    // if (fCheckMaxDcaXY) {
    //   float dca = track.dcaXY();
    //   float maxDCA = 0.0105f + 0.0350f / pow(pt, 1.1f);
    //   if (dca < maxDCA) {
    //     return false;
    //   }
    // }
    // bool checkDCAZ = track.dcaZ() > fDcaZLow && track.dcaZ() < fDcaZHigh;
    // if (!checkDCAZ) {
    //   return false;
    // }
    return true;
  }

  template <typename TMcCollisions, typename TMcParticles, typename TBCs>
  void skimMCInfo(TMcCollisions const& mcCollisions,
                  TMcParticles const& mcParticles,
                  TBCs const& bcs,
                  std::map<int32_t, int32_t>& newPartIDs)
  {
    std::vector<int32_t> newEventIDs(mcCollisions.size(), -1);

    int32_t newPartID = 0;
    int32_t newEventID = 0;
    int32_t nMCParticles = mcParticles.size();
    for (int32_t mcPartID = 0; mcPartID < nMCParticles; mcPartID++) {
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      const auto& mcEvent = mcCollisions.iteratorAt(mcEventID);
      bool isSignal = mcEvent.generatorsID() == fSignalGenID;
      if (!isSignal || !mcPart.producedByGenerator()) {
        continue;
      }
      newPartIDs[mcPartID] = newPartID;
      newPartID++;
      if (newEventIDs[mcEventID] == -1) {
        newEventIDs[mcEventID] = newEventID;
        newEventID++;
      }
    }

    // storing MC particles
    for (const auto& item : newPartIDs) {
      int32_t mcPartID = item.first;
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      int32_t newEventID = newEventIDs[mcEventID];
      // collecting new mother IDs
      const auto& motherIDs = mcPart.mothersIds();
      std::vector<int32_t> newMotherIDs;
      newMotherIDs.reserve(motherIDs.size());
      for (auto motherID : motherIDs) {
        if (motherID >= nMCParticles) {
          continue;
        }
        auto it = newPartIDs.find(motherID);
        if (it != newPartIDs.end()) {
          newMotherIDs.push_back(it->second);
        }
      }
      // collecting new daughter IDs
      const auto& daughterIDs = mcPart.daughtersIds();
      int32_t newDaughterIDs[2] = {-1, -1};
      if (daughterIDs.size() > 0) {
        int32_t firstDaughter = daughterIDs.front();
        int32_t lastDaughter = daughterIDs.back();
        if (firstDaughter >= nMCParticles || lastDaughter >= nMCParticles) {
          continue;
        }
        auto itFirst = newPartIDs.find(firstDaughter);
        auto itLast = newPartIDs.find(lastDaughter);
        if (itFirst != newPartIDs.end() && itLast != newPartIDs.end()) {
          newDaughterIDs[0] = newPartIDs.at(daughterIDs.front());
          newDaughterIDs[1] = newPartIDs.at(daughterIDs.back());
        }
      }
      skMCParticles(newEventID, mcPart.pdgCode(), mcPart.statusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
    }

    // storing MC events
    for (int32_t i = 0; i < mcCollisions.size(); i++) {
      if (newEventIDs[i] == -1) {
        continue;
      }
      const auto& mcEvent = mcCollisions.iteratorAt(i);
      skMCEvents(mcEvent.bc().globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                 mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  template <int32_t TTrackType, typename TFwdTracks, typename TAmbFwdTracks, typename TBCs, typename TMcFwdTrackLabels>
  void skimFwdTracks(TFwdTracks const& tracks,
                     TAmbFwdTracks const& ambTracks,
                     TBCs const& bcs,
                     TMcFwdTrackLabels* mcTrackLabels,
                     std::map<int32_t, int32_t> const& newPartIDs)
  {
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
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trId);
        uint16_t mcMask = label.mcMask();
        auto it = newPartIDs.find(label.mcParticleId());
        int32_t newPartID = it != newPartIDs.end() ? it->second : -1;
        muonLabels(newPartID, mcMask);
      }
    }
  }

  template <typename TBarTracks, typename TAmbBarTracks, typename TBCs, typename TMcBarTrackLabels>
  void skimBarTracks(TBarTracks const& tracks,
                     TAmbBarTracks const& ambTracks,
                     TBCs const& bcs,
                     TMcBarTrackLabels* mcTrackLabels,
                     std::map<int32_t, int32_t> const& newPartIDs,
                     o2::aod::McParticles const& mcParticles,
                     o2::aod::McCollisions const& mcCollisions)
  {
    // debug
    LOGF(info, "n barrel tracks = %d", tracks.size());

    for (const auto& ambTr : ambTracks) {
      auto trId = ambTr.trackId();
      const auto& tr = tracks.iteratorAt(trId);
      // using only tracks with TOF match
      if (!tr.hasTOF()) {
        continue;
      }
      // skip if doesn't pass cuts
      if (!applyBarCuts(tr)) {
        continue;
      }
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      uint64_t firstBC = first.globalBC();
      double trTime = firstBC * o2::constants::lhc::LHCBunchSpacingNS + (double)tr.trackTime(); // using absolute time
      barTracks(tr.px(), tr.py(), tr.pz(), tr.sign(), trTime, tr.trackTimeRes());
      barsCov(tr.x(), tr.alpha(), tr.y(), tr.z(), tr.snp(), tr.tgl(), tr.signed1Pt(),
              tr.cYY(), tr.cZY(), tr.cZZ(), tr.cSnpY(), tr.cSnpZ(), tr.cSnpSnp(), tr.cTglY(),
              tr.cTglZ(), tr.cTglSnp(), tr.cTglTgl(), tr.c1PtY(), tr.c1PtZ(), tr.c1PtSnp(),
              tr.c1PtTgl(), tr.c1Pt21Pt2());
      barsExtra(tr.flags(), tr.itsClusterMap(), tr.tpcNClsFindable(),
                tr.tpcNClsFindableMinusFound(), tr.tpcNClsFindableMinusCrossedRows(),
                tr.tpcNClsShared(), tr.itsChi2NCl(), tr.tpcChi2NCl(), tr.tofChi2(),
                tr.tpcSignal(), tr.length(), tr.tofExpMom());
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trId);
        uint16_t mcMask = label.mcMask();
        int32_t mcPartID = label.mcParticleId();
        // debug
        if (mcPartID >= 0) {
          const auto& mcPart = mcParticles.iteratorAt(mcPartID);
          bool isSignal = mcPart.mcCollision().generatorsID() == fSignalGenID;
          if (isSignal && mcPart.producedByGenerator()) {
            hNSig->Fill(1);
            if (tr.hasITS() && tr.hasTPC())
              hNSigITSTPC->Fill(1);
          }
        }
        //
        auto it = newPartIDs.find(mcPartID);
        int32_t newPartID = it != newPartIDs.end() ? it->second : -1;
        barLabels(newPartID, mcMask);
      }
    }
  }

  void processFwdMC(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& tracks,
                    o2::aod::AmbiguousFwdTracks const& ambTracks,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
                    o2::aod::McCollisions const& mcCollisions,
                    o2::aod::McParticles const& mcParticles,
                    o2::aod::BCs const& bcs)
  {
    fDoMC = true;
    using namespace o2::aod::fwdtrack;
    const int32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    std::map<int32_t, int32_t> newPartIDs;
    skimMCInfo(mcCollisions, mcParticles, bcs, newPartIDs);
    skimFwdTracks<trType>(tracks, ambTracks, bcs, &mcFwdTrackLabels, newPartIDs);
    newPartIDs.clear();
  }

  void processFwd(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& tracks,
                  o2::aod::AmbiguousFwdTracks const& ambTracks,
                  o2::aod::BCs const& bcs)
  {
    fDoMC = false;
    using namespace o2::aod::fwdtrack;
    const int32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    std::map<int32_t, int32_t> dummyMap;
    skimFwdTracks<trType>(tracks, ambTracks, bcs, (o2::aod::McFwdTrackLabels*)nullptr, dummyMap);
  }

  // todo: switch to propagated tracks?
  void processAllMC(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& fwdTracks,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
                    o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                    o2::soa::Join<o2::aod::TracksIU, o2::aod::TracksCovIU, o2::aod::TracksExtra> const& barTracks,
                    o2::aod::McTrackLabels const& mcBarTrackLabels,
                    o2::aod::AmbiguousTracks const& ambBarTracks,
                    o2::aod::McCollisions const& mcCollisions,
                    o2::aod::McParticles const& mcParticles,
                    o2::aod::BCs const& bcs)
  {
    fDoMC = true;
    using namespace o2::aod::fwdtrack;
    const int32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    std::map<int32_t, int32_t> newPartIDs;
    skimMCInfo(mcCollisions, mcParticles, bcs, newPartIDs);
    skimFwdTracks<trType>(fwdTracks, ambFwdTracks, bcs, &mcFwdTrackLabels, newPartIDs);
    skimBarTracks(barTracks, ambBarTracks, bcs, &mcBarTrackLabels, newPartIDs, mcParticles, mcCollisions);
    newPartIDs.clear();
  }

  PROCESS_SWITCH(TrackSkimmer, processFwdMC, "Produce only muon tracks with MC information", false);
  PROCESS_SWITCH(TrackSkimmer, processFwd, "Produce only muon tracks", false);
  PROCESS_SWITCH(TrackSkimmer, processAllMC, "Produce barrel and muon tracks with MC information", false);
  // todo: process without MC info
  // PROCESS_SWITCH(TrackSkimmer, processAll, "Produce barrel and muon tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<TrackSkimmer>(cfgc)};
}