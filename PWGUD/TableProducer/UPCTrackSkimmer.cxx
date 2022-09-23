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

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcTrackSkimmer {
  int32_t fSignalGenID{1};
  bool fDoMC;

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;

  Produces<o2::aod::UDFwdTracks> fwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> fwdTracksExtra;
  Produces<o2::aod::UDMcFwdTrackLabels> fwdTrackLabels;

  Produces<o2::aod::UDTracks> udTracks;
  Produces<o2::aod::UDTracksExtra> udTracksExtra;
  Produces<o2::aod::UDTracksPID> udTracksPID;
  Produces<o2::aod::UDMcTrackLabels> udTrackLabels;

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
  // quality: DCA
  Configurable<int> fCheckMaxDcaXY{"checkMaxDCACut", 1, "Apply cut on maximal DCA_xy"};
  Configurable<float> fDcaZLow{"dcaZLow", -3., "Minimal DCA_z for barrel tracks"};
  Configurable<float> fDcaZHigh{"dcaZHigh", 3., "Maximal DCA_z for barrel tracks"};

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
    // check DCA
    if (fCheckMaxDcaXY) {
      float dca = track.dcaXY();
      float maxDCA = 0.0105f + 0.0350f / pow(pt, 1.1f);
      if (dca < maxDCA) {
        return false;
      }
    }
    bool checkDCAZ = track.dcaZ() > fDcaZLow && track.dcaZ() < fDcaZHigh;
    if (!checkDCAZ) {
      return false;
    }
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
    // loop over MC particles to select only the ones from signal events
    // and calculate new MC table IDs
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
      udMCParticles(newEventID, mcPart.pdgCode(), mcPart.statusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
    }

    // storing MC events
    for (int32_t i = 0; i < mcCollisions.size(); i++) {
      if (newEventIDs[i] == -1) {
        continue;
      }
      const auto& mcEvent = mcCollisions.iteratorAt(i);
      udMCCollisions(mcEvent.bc().globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
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
      double realTime = firstBC * o2::constants::lhc::LHCBunchSpacingNS + (double)tr.trackTime(); // "real" = relative to start of TF
      uint64_t bc = realTime > 0. ? std::round(realTime / o2::constants::lhc::LHCBunchSpacingNS) : 0;
      double trTime = bc * o2::constants::lhc::LHCBunchSpacingNS - realTime;
      fwdTracks(tr.px(), tr.py(), tr.pz(), tr.sign(), bc, trTime, tr.trackTimeRes());
      fwdTracksExtra(tr.nClusters(), tr.pDca(), tr.rAtAbsorberEnd(), tr.chi2(), tr.chi2MatchMCHMID(),
                     tr.mchBitMap(), tr.midBitMap(), tr.midBoards());
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trId);
        uint16_t mcMask = label.mcMask();
        auto it = newPartIDs.find(label.mcParticleId());
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != newPartIDs.end() ? it->second : -1;
        fwdTrackLabels(newPartID, mcMask);
      }
    }
  }

  template <typename TBarTracks, typename TAmbBarTracks, typename TBCs, typename TMcBarTrackLabels>
  void skimBarTracks(TBarTracks const& tracks,
                     TAmbBarTracks const& ambTracks,
                     TBCs const& bcs,
                     TMcBarTrackLabels* mcTrackLabels,
                     std::map<int32_t, int32_t> const& newPartIDs)
  {
    for (const auto& ambTr : ambTracks) {
      auto trId = ambTr.trackId();
      const auto& tr = tracks.iteratorAt(trId);
      // using only tracks with TOF match
      // todo: make this selection optional?
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
      double realTime = firstBC * o2::constants::lhc::LHCBunchSpacingNS + (double)tr.trackTime(); // "real" = relative to start of TF
      uint64_t bc = realTime > 0. ? std::round(realTime / o2::constants::lhc::LHCBunchSpacingNS) : 0;
      double trTime = bc * o2::constants::lhc::LHCBunchSpacingNS - realTime;
      udTracks(tr.px(), tr.py(), tr.pz(), tr.sign(), bc, trTime, tr.trackTimeRes());
      udTracksExtra(tr.itsClusterMap(), tr.tpcNClsFindable(), tr.tpcNClsFindableMinusFound(), tr.tpcNClsFindableMinusCrossedRows(),
                    tr.tpcNClsShared(), tr.trdPattern(), tr.itsChi2NCl(), tr.tpcChi2NCl(), tr.trdChi2(), tr.tofChi2(),
                    tr.tpcSignal(), tr.tofSignal(), tr.trdSignal(), tr.length(), tr.tofExpMom(), tr.detectorMap());
      udTracksPID(tr.tpcNSigmaEl(), tr.tpcNSigmaMu(), tr.tpcNSigmaPi(), tr.tpcNSigmaKa(), tr.tpcNSigmaPr(),
                  tr.tofNSigmaEl(), tr.tofNSigmaMu(), tr.tofNSigmaPi(), tr.tofNSigmaKa(), tr.tofNSigmaPr());
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trId);
        uint16_t mcMask = label.mcMask();
        int32_t mcPartID = label.mcParticleId();
        auto it = newPartIDs.find(mcPartID);
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != newPartIDs.end() ? it->second : -1;
        udTrackLabels(newPartID, mcMask);
      }
    }
  }

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

  // process only MCH-MID tracks with MC information
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

  // process only MCH-MID tracks without MC information
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

  // process both barrel and muon tracks with MC information
  void processAllMC(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& fwdTracks,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
                    o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                    BarrelTracks const& barTracks,
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
    skimBarTracks(barTracks, ambBarTracks, bcs, &mcBarTrackLabels, newPartIDs);
    newPartIDs.clear();
  }

  // process both barrel and muon tracks without MC information
  void processAll(o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov> const& fwdTracks,
                  o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                  BarrelTracks const& barTracks,
                  o2::aod::AmbiguousTracks const& ambBarTracks,
                  o2::aod::BCs const& bcs)
  {
    fDoMC = true;
    using namespace o2::aod::fwdtrack;
    const int32_t trType = ForwardTrackTypeEnum::MuonStandaloneTrack;
    std::map<int32_t, int32_t> dummyMap;
    skimFwdTracks<trType>(fwdTracks, ambFwdTracks, bcs, (o2::aod::McFwdTrackLabels*)nullptr, dummyMap);
    skimBarTracks(barTracks, ambBarTracks, bcs, (o2::aod::McTrackLabels*)nullptr, dummyMap);
  }

  PROCESS_SWITCH(UpcTrackSkimmer, processFwdMC, "Produce only muon tracks with MC information", false);
  PROCESS_SWITCH(UpcTrackSkimmer, processFwd, "Produce only muon tracks", false);
  PROCESS_SWITCH(UpcTrackSkimmer, processAllMC, "Produce barrel and muon tracks with MC information", false);
  PROCESS_SWITCH(UpcTrackSkimmer, processAll, "Produce barrel and muon tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<UpcTrackSkimmer>(cfgc)};
}
