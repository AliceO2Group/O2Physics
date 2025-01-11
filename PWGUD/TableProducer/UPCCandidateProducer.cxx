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
/// \author Nazar Burmasov, nazar.burmasov@cern.ch
/// \author Diana Krupova, diana.krupova@cern.ch
/// \since 04.06.2024

#include <limits>
#include <unordered_set>
#include <utility>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <vector>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "DataFormatsFIT/Triggers.h"
#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"
#include "DataFormatsITSMFT/ROFRecord.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducer {
  bool fDoMC{false};

  std::map<int32_t, int32_t> fNewPartIDs;
  uint64_t fMaxBC{0}; // max BC for ITS-TPC search

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;

  Produces<o2::aod::UDFwdTracks> udFwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> udFwdTracksExtra;
  Produces<o2::aod::UDFwdIndices> udFwdIndices;
  Produces<o2::aod::UDFwdTracksCls> udFwdTrkClusters;
  Produces<o2::aod::UDMcFwdTrackLabels> udFwdTrackLabels;

  Produces<o2::aod::UDTracks> udTracks;
  Produces<o2::aod::UDTracksExtra> udTracksExtra;
  Produces<o2::aod::UDTracksDCA> udTracksDCA;
  Produces<o2::aod::UDTracksPID> udTracksPID;
  Produces<o2::aod::UDMcTrackLabels> udTrackLabels;
  Produces<o2::aod::UDTracksFlags> udTracksFlags;

  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSels> eventCandidatesSels;
  Produces<o2::aod::UDCollisionsSelsCent> eventCandidatesSelsCent;
  Produces<o2::aod::UDCollisionsSelsFwd> eventCandidatesSelsFwd;
  Produces<o2::aod::UDCollisionSelExtras> eventCandidatesSelExtras;

  Produces<o2::aod::UDZdcsReduced> udZdcsReduced;

  std::vector<bool> fwdSelectors;
  std::vector<bool> barrelSelectors;

  // skimmer flags
  // choose a source of signal MC events
  Configurable<int> fSignalGenID{"signalGenID", 1, "Signal generator ID"};

  // load cuts
  UPCCutparHolder upcCuts = UPCCutparHolder();
  MutableConfigurable<UPCCutparHolder> inputCuts{"UPCCuts", {}, "UPC event cuts"};

  // candidate producer flags
  Configurable<int> fFilterFT0{"filterFT0", -1, "Filter candidates by FT0 TOR(central) or T0A(fwd)"};
  Configurable<int> fFilterTSC{"filterTSC", -1, "Filter candidates by FT0 TSC"};
  Configurable<int> fFilterTVX{"filterTVX", -1, "Filter candidates by FT0 TVX"};
  Configurable<int> fFilterFV0{"filterFV0", -1, "Filter candidates by FV0A"};

  Configurable<uint64_t> fBCWindowFITAmps{"bcWindowFITAmps", 20, "BC range for T0A/V0A amplitudes array [-range, +(range-1)]"};
  Configurable<int> fBcWindowMCH{"bcWindowMCH", 20, "Time window for MCH-MID to MCH-only matching for Muon candidates"};
  Configurable<int> fBcWindowITSTPC{"bcWindowITSTPC", 20, "Time window for TOF/ITS-TPC to ITS-TPC matching for Central candidates"};

  Configurable<int> fMuonTrackTShift{"muonTrackTShift", 0, "Time shift for Muon tracks"};
  Configurable<int> fBarrelTrackTShift{"barrelTrackTShift", 0, "Time shift for Central Barrel tracks"};

  Configurable<uint32_t> fNFwdProngs{"nFwdProngs", 0, "Matched forward tracks per candidate"};
  Configurable<uint32_t> fNBarProngs{"nBarProngs", 2, "Matched barrel tracks per candidate"};

  Configurable<int> fFilterRangeFT0{"filterRangeFT0", 0, "BC range (+/-) for filtration by FT0 signals"};
  Configurable<int> fSearchITSTPC{"searchITSTPC", 0, "Search for ITS-TPC tracks near candidates"};
  Configurable<int> fSearchRangeITSTPC{"searchRangeITSTPC", 50, "BC range for ITS-TPC tracks search wrt TOF tracks"};

  Configurable<float> fMinEtaMFT{"minEtaMFT", -3.6, "Minimum eta for MFT tracks"};
  Configurable<float> fMaxEtaMFT{"maxEtaMFT", -2.5, "Maximum eta for MFT tracks"};

  // QA histograms
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels, o2::aod::BCFlags>;

  using ForwardTracks = o2::soa::Join<o2::aod::UDFwdTracksProp, o2::aod::UDFwdTracksCovProp>;

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFbeta,
                                     o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

  typedef std::pair<uint64_t, std::vector<int64_t>> BCTracksPair;

  void init(InitContext&)
  {
    fwdSelectors.resize(upchelpers::kNFwdSels - 1, false);
    barrelSelectors.resize(upchelpers::kNBarrelSels - 1, false);

    upcCuts = (UPCCutparHolder)inputCuts;

    const AxisSpec axisTrgCounters{10, 0.5, 10.5, ""};
    histRegistry.add("hCountersTrg", "", kTH1F, {axisTrgCounters});
    histRegistry.get<TH1>(HIST("hCountersTrg"))->GetXaxis()->SetBinLabel(1, "TCE");
    histRegistry.get<TH1>(HIST("hCountersTrg"))->GetXaxis()->SetBinLabel(2, "ZNA");
    histRegistry.get<TH1>(HIST("hCountersTrg"))->GetXaxis()->SetBinLabel(3, "ZNC");

    const AxisSpec axisBcDist{201, 0.5, 200.5, ""};
    histRegistry.add("hDistToITSTPC", "", kTH1F, {axisBcDist});

    const AxisSpec axisSelFwd{upchelpers::kNFwdSels, 0., static_cast<double>(upchelpers::kNFwdSels), ""};
    histRegistry.add("MuonsSelCounter", "", kTH1F, {axisSelFwd});
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelRabs + 1, "Rabs");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelpDCA + 1, "pDCA");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelChi2 + 1, "Chi2");

    const AxisSpec axisSelBar{upchelpers::kNBarrelSels, 0., static_cast<double>(upchelpers::kNBarrelSels), ""};
    histRegistry.add("BarrelsSelCounter", "", kTH1F, {axisSelBar});
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelHasTOF + 1, "HasTOF");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelITSNCls + 1, "ITSNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelITSChi2 + 1, "ITSChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelTPCNCls + 1, "TPCNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelTPCChi2 + 1, "TPCChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelDCAXY + 1, "DCAXY");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelDCAZ + 1, "DCAZ");
  }

  template <typename T>
  bool applyFwdCuts(const T& track)
  {
    histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelAll, 1);

    // using any cuts at all?
    if (!upcCuts.getUseFwdCuts()) {
      return true;
    }

    upchelpers::applyFwdCuts(upcCuts, track, fwdSelectors);

    if (fwdSelectors[upchelpers::kFwdSelPt])
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelPt, 1);
    if (fwdSelectors[upchelpers::kFwdSelEta])
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelEta, 1);
    if (fwdSelectors[upchelpers::kFwdSelRabs])
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelRabs, 1);
    if (fwdSelectors[upchelpers::kFwdSelpDCA])
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelpDCA, 1);
    if (fwdSelectors[upchelpers::kFwdSelChi2])
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelChi2, 1);
    bool pass = fwdSelectors[upchelpers::kFwdSelPt] &&
                fwdSelectors[upchelpers::kFwdSelEta] &&
                fwdSelectors[upchelpers::kFwdSelRabs] &&
                fwdSelectors[upchelpers::kFwdSelpDCA] &&
                fwdSelectors[upchelpers::kFwdSelChi2];
    return pass;
  }

  template <typename T>
  bool applyBarCuts(const T& track)
  {
    // using any cuts at all?
    if (!upcCuts.getUseBarCuts())
      return true;
    if (upcCuts.getAmbigSwitch() == 1 && !track.isPVContributor())
      return false;
    if (upcCuts.getRequireTOF() && !track.hasTOF())
      return false;
    if (track.pt() < upcCuts.getBarPtLow())
      return false;
    if (track.pt() > upcCuts.getBarPtHigh())
      return false;
    if (track.eta() < upcCuts.getBarEtaLow())
      return false;
    if (track.eta() > upcCuts.getBarEtaHigh())
      return false;
    if (track.itsNCls() < static_cast<uint8_t>(upcCuts.getITSNClusLow()))
      return false;
    if (track.itsNCls() > static_cast<uint8_t>(upcCuts.getITSNClusHigh()))
      return false;
    if (track.itsChi2NCl() < upcCuts.getITSChi2Low())
      return false;
    if (track.itsChi2NCl() > upcCuts.getITSChi2High())
      return false;
    if (track.tpcNClsFound() < static_cast<int16_t>(upcCuts.getTPCNClsLow()))
      return false;
    if (track.tpcNClsFound() > static_cast<int16_t>(upcCuts.getTPCNClsHigh()))
      return false;
    if (track.tpcChi2NCl() < upcCuts.getTPCChi2Low())
      return false;
    if (track.tpcChi2NCl() > upcCuts.getTPCChi2High())
      return false;
    if (track.dcaZ() < upcCuts.getDcaZLow())
      return false;
    if (track.dcaZ() > upcCuts.getDcaZHigh())
      return false;
    if (upcCuts.getCheckMaxDcaXY()) {
      float dca = track.dcaXY();
      float maxDCA = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
      if (dca > maxDCA)
        return false;
    }
    return true;
  }

  auto findClosestBC(uint64_t globalBC, std::map<uint64_t, int32_t>& bcs)
  {
    auto it = bcs.lower_bound(globalBC);
    auto bc1 = it->first;
    if (it != bcs.begin())
      --it;
    auto bc2 = it->first;
    auto dbc1 = bc1 >= globalBC ? bc1 - globalBC : globalBC - bc1;
    auto dbc2 = bc2 >= globalBC ? bc2 - globalBC : globalBC - bc2;
    auto bc = (dbc1 <= dbc2) ? bc1 : bc2;
    return bc;
  }

  auto findClosestTrackBCiter(uint64_t globalBC, std::vector<BCTracksPair>& bcs)
  {
    auto it = std::lower_bound(bcs.begin(), bcs.end(), globalBC,
                               [](const BCTracksPair& p, uint64_t bc) {
                                 return p.first < bc;
                               });
    auto bc1 = it->first;
    auto it1 = it;
    if (it != bcs.begin())
      --it;
    auto it2 = it;
    auto bc2 = it->first;
    auto dbc1 = bc1 >= globalBC ? bc1 - globalBC : globalBC - bc1;
    auto dbc2 = bc2 >= globalBC ? bc2 - globalBC : globalBC - bc2;
    return (dbc1 <= dbc2) ? it1 : it2;
  }

  auto findClosestTrackBCiterNotEq(uint64_t globalBC, std::vector<BCTracksPair>& bcs)
  {
    auto it = std::find_if(
      bcs.begin(),
      bcs.end(),
      [globalBC](const auto& v) {
        return v.first > globalBC;
      });
    auto bc1 = it->first;
    auto it1 = it;
    if (it != bcs.begin())
      --it;
    if (it->first == globalBC)
      --it;
    auto it2 = it;
    auto bc2 = it->first;
    auto dbc1 = bc1 >= globalBC ? bc1 - globalBC : globalBC - bc1;
    auto dbc2 = bc2 >= globalBC ? bc2 - globalBC : globalBC - bc2;
    return (dbc1 <= dbc2) ? it1 : it2;
  }

  template <typename TBCs>
  void skimMCInfo(o2::aod::McCollisions const& mcCollisions,
                  o2::aod::McParticles const& mcParticles,
                  TBCs const& /*bcs*/)
  {
    std::vector<int32_t> newEventIDs(mcCollisions.size(), -1);

    int32_t newPartID = 0;
    int32_t newEventID = 0;
    int32_t nMCParticles = mcParticles.size();
    // loop over MC particles to select only the ones from signal events
    // and calculate new MC table IDs
    for (int32_t mcPartID = 0; mcPartID < nMCParticles; mcPartID++) {
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      if (!mcPart.has_mcCollision())
        continue;
      int32_t mcEventID = mcPart.mcCollisionId();
      const auto& mcEvent = mcCollisions.iteratorAt(mcEventID);
      bool isSignal = mcEvent.getSourceId() == fSignalGenID;
      if (!isSignal) {
        continue;
      }
      fNewPartIDs[mcPartID] = newPartID;
      newPartID++;
      if (newEventIDs[mcEventID] == -1) {
        newEventIDs[mcEventID] = newEventID;
        newEventID++;
      }
    }

    std::vector<int32_t> newMotherIDs{};

    // storing MC particles
    for (const auto& item : fNewPartIDs) {
      int32_t mcPartID = item.first;
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      int32_t newEventID = newEventIDs[mcEventID];
      // collecting new mother IDs
      if (mcPart.has_mothers()) {
        const auto& motherIDs = mcPart.mothersIds();
        for (auto motherID : motherIDs) {
          if (motherID >= nMCParticles) {
            continue;
          }
          auto it = fNewPartIDs.find(motherID);
          if (it != fNewPartIDs.end()) {
            newMotherIDs.push_back(it->second);
          }
        }
      }
      // collecting new daughter IDs
      int32_t newDaughterIDs[2] = {-1, -1};
      if (mcPart.has_daughters()) {
        const auto& daughterIDs = mcPart.daughtersIds();
        int32_t firstDaughter = daughterIDs.front();
        int32_t lastDaughter = daughterIDs.back();
        if (firstDaughter >= nMCParticles || lastDaughter >= nMCParticles) {
          continue;
        }
        auto itFirst = fNewPartIDs.find(firstDaughter);
        auto itLast = fNewPartIDs.find(lastDaughter);
        if (itFirst != fNewPartIDs.end() && itLast != fNewPartIDs.end()) {
          newDaughterIDs[0] = fNewPartIDs.at(daughterIDs.front());
          newDaughterIDs[1] = fNewPartIDs.at(daughterIDs.back());
        }
      }
      udMCParticles(newEventID, mcPart.pdgCode(), mcPart.getHepMCStatusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
      newMotherIDs.clear();
    }

    // storing MC events
    for (int32_t i = 0; i < mcCollisions.size(); i++) {
      if (newEventIDs[i] == -1) {
        continue;
      }
      const auto& mcEvent = mcCollisions.iteratorAt(i);
      const auto& bc = mcEvent.bc_as<TBCs>();
      udMCCollisions(bc.globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                     mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  void fillFwdTracks(ForwardTracks const& tracks,
                     std::vector<int64_t> const& trackIDs,
                     int32_t candID,
                     uint64_t globalBC, uint64_t closestBcMCH,
                     const o2::aod::McFwdTrackLabels* mcTrackLabels)
  {
    for (auto trackID : trackIDs) {
      const auto& track = tracks.iteratorAt(trackID);
      double trTime = track.trackTime();
      double mchmidChi2 = track.chi2MatchMCHMID();
      if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack) {
        trTime = (static_cast<int64_t>(globalBC) - static_cast<int64_t>(closestBcMCH)) * o2::constants::lhc::LHCBunchSpacingNS; // track time relative to MCH-MID track
        mchmidChi2 = -999.;                                                                                                     // no MID match
      }
      double mchmftChi2 = track.chi2MatchMCHMFT();
      int64_t globalIndex = track.globalIndex();
      int64_t mchIndex = track.matchMCHTrackId();
      int64_t mftIndex = track.matchMFTTrackId();
      udFwdTracks(candID, track.px(), track.py(), track.pz(), track.sign(), globalBC, trTime, track.trackTimeRes());
      udFwdTracksExtra(track.trackType(), track.nClusters(), track.pDca(), track.rAtAbsorberEnd(), track.chi2(), mchmidChi2, mchmftChi2,
                       track.mchBitMap(), track.midBitMap(), track.midBoards());
      udFwdIndices(candID, globalIndex, mchIndex, mftIndex);
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trackID);
        uint16_t mcMask = label.mcMask();
        auto it = fNewPartIDs.find(label.mcParticleId());
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
        udFwdTrackLabels(newPartID, mcMask);
      }
    }
  }

  void fillFwdClusters(const std::vector<int>& trackIds,
                       o2::aod::FwdTrkCls const& fwdTrkCls)
  {
    std::map<int, std::vector<int>> clustersPerTrack;
    for (const auto& cls : fwdTrkCls) {
      clustersPerTrack[cls.fwdtrackId()].push_back(cls.globalIndex());
    }
    int newId = 0;
    for (auto trackId : trackIds) {
      const auto& clusters = clustersPerTrack.at(trackId);
      for (auto clsId : clusters) {
        const auto& clsInfo = fwdTrkCls.iteratorAt(clsId);
        udFwdTrkClusters(newId, clsInfo.x(), clsInfo.y(), clsInfo.z(), clsInfo.clInfo());
      }
      newId++;
    }
  }

  void fillBarrelTracks(BarrelTracks const& tracks,
                        std::vector<int64_t> const& trackIDs,
                        int32_t candID,
                        uint64_t globalBC,
                        uint64_t closestBcITSTPC,
                        const o2::aod::McTrackLabels* mcTrackLabels,
                        std::unordered_map<int64_t, uint64_t>& /*ambBarrelTrBCs*/)
  {
    for (auto trackID : trackIDs) {
      const auto& track = tracks.iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      int64_t colId = track.collisionId() >= 0 ? track.collisionId() : -1;
      if (!track.hasTOF() && closestBcITSTPC != std::numeric_limits<uint64_t>::max())
        trTime = (static_cast<int64_t>(globalBC) - static_cast<int64_t>(closestBcITSTPC)) * o2::constants::lhc::LHCBunchSpacingNS; // track time relative to TOF track
      udTracks(candID, track.px(), track.py(), track.pz(), track.sign(), globalBC, trTime, track.trackTimeRes());
      udTracksExtra(track.tpcInnerParam(), track.itsClusterSizes(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.trdPattern(), track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                    track.tpcSignal(), track.tofSignal(), track.trdSignal(), track.length(), track.tofExpMom(), track.detectorMap());
      udTracksPID(track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                  track.beta(), track.betaerror(),
                  track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
      udTracksDCA(track.dcaZ(), track.dcaXY());
      udTracksFlags(colId, track.isPVContributor());
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trackID);
        uint16_t mcMask = label.mcMask();
        int32_t mcPartID = label.mcParticleId();
        auto it = fNewPartIDs.find(mcPartID);
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
        udTrackLabels(newPartID, mcMask);
      }
    }
  }

  bool checkFT0(upchelpers::FITInfo& info, bool isCentral)
  {
    const uint64_t presBitNum = 16;
    bool hasNoFT0 = true;
    for (uint64_t ibit = presBitNum - fFilterRangeFT0; ibit <= presBitNum + fFilterRangeFT0; ibit++) {
      bool check = false;
      if (isCentral) {
        bool isBB = TESTBIT(info.BBFT0Apf, ibit) || TESTBIT(info.BBFT0Cpf, ibit);
        bool isBG = TESTBIT(info.BGFT0Apf, ibit) || TESTBIT(info.BGFT0Cpf, ibit);
        check = !isBB && !isBG;
      } else {
        bool checkA = TESTBIT(info.BGFT0Apf, ibit);
        bool checkC = TESTBIT(info.BBFT0Cpf, ibit);
        check = checkA && checkC;
      }
      if (!check) {
        hasNoFT0 = false;
        break;
      }
    }
    return hasNoFT0;
  }

  template <typename TBCs>
  void processFITInfo(upchelpers::FITInfo& fitInfo,
                      uint64_t midbc,
                      std::vector<std::pair<uint64_t, int64_t>>& v,
                      TBCs const& bcs,
                      o2::aod::FT0s const& /*ft0s*/,
                      o2::aod::FDDs const& /*fdds*/,
                      o2::aod::FV0As const& /*fv0as*/)
  {
    auto it = std::find_if(v.begin(),
                           v.end(),
                           [midbc](const std::pair<uint64_t, int64_t>& p) { return p.first == midbc; });

    if (it != v.end()) {
      auto bcId = it->second;
      auto bcEntry = bcs.iteratorAt(bcId);
      if (bcEntry.has_foundFT0()) {
        auto ft0 = bcEntry.foundFT0();
        fitInfo.timeFT0A = ft0.timeA();
        fitInfo.timeFT0C = ft0.timeC();
        const auto& ampsA = ft0.amplitudeA();
        const auto& ampsC = ft0.amplitudeC();
        fitInfo.ampFT0A = 0.;
        for (auto amp : ampsA)
          fitInfo.ampFT0A += amp;
        fitInfo.ampFT0C = 0.;
        for (auto amp : ampsC)
          fitInfo.ampFT0C += amp;
        fitInfo.triggerMaskFT0 = ft0.triggerMask();
      }
      if (bcEntry.has_foundFV0()) {
        auto fv0a = bcEntry.foundFV0();
        fitInfo.timeFV0A = fv0a.time();
        const auto& amps = fv0a.amplitude();
        fitInfo.ampFV0A = 0.;
        for (auto amp : amps)
          fitInfo.ampFV0A += amp;
        fitInfo.triggerMaskFV0A = fv0a.triggerMask();
      }
      if (bcEntry.has_foundFDD()) {
        auto fdd = bcEntry.foundFDD();
        fitInfo.timeFDDA = fdd.timeA();
        fitInfo.timeFDDC = fdd.timeC();
        const auto& ampsA = fdd.chargeA();
        const auto& ampsC = fdd.chargeC();
        fitInfo.ampFDDA = 0.;
        for (auto amp : ampsA) {
          fitInfo.ampFDDA += amp;
        }
        fitInfo.ampFDDC = 0.;
        for (auto amp : ampsC) {
          fitInfo.ampFDDC += amp;
        }
        fitInfo.triggerMaskFDD = fdd.triggerMask();
      }
    }

    const uint64_t range = 16;
    uint64_t left = midbc >= range ? midbc - range : 0;
    uint64_t right = fMaxBC >= midbc + range ? midbc + range : fMaxBC;

    std::pair<uint64_t, int64_t> dummyPair(left, 0);
    auto curit = std::lower_bound(v.begin(), v.end(), dummyPair,
                                  [](const std::pair<uint64_t, int64_t>& left, const std::pair<uint64_t, int64_t>& right) { return left.first < right.first; });

    if (curit == v.end()) // no BCs with FT0 info at all
      return;

    uint64_t curbc = curit->first;
    while (curbc <= right) {
      uint64_t bit = curbc - (midbc - range);
      int64_t bcGlId = curit->second;
      const auto& bc = bcs.iteratorAt(bcGlId);
      if (!bc.selection_bit(o2::aod::evsel::kNoBGT0A))
        SETBIT(fitInfo.BGFT0Apf, bit);
      if (!bc.selection_bit(o2::aod::evsel::kNoBGT0C))
        SETBIT(fitInfo.BGFT0Cpf, bit);
      if (bc.selection_bit(o2::aod::evsel::kIsBBT0A))
        SETBIT(fitInfo.BBFT0Apf, bit);
      if (bc.selection_bit(o2::aod::evsel::kIsBBT0C))
        SETBIT(fitInfo.BBFT0Cpf, bit);
      if (!bc.selection_bit(o2::aod::evsel::kNoBGV0A))
        SETBIT(fitInfo.BGFV0Apf, bit);
      if (bc.selection_bit(o2::aod::evsel::kIsBBV0A))
        SETBIT(fitInfo.BBFV0Apf, bit);
      if (!bc.selection_bit(o2::aod::evsel::kNoBGFDA))
        SETBIT(fitInfo.BGFDDApf, bit);
      if (!bc.selection_bit(o2::aod::evsel::kNoBGFDC))
        SETBIT(fitInfo.BGFDDCpf, bit);
      if (bc.selection_bit(o2::aod::evsel::kIsBBFDA))
        SETBIT(fitInfo.BBFDDApf, bit);
      if (bc.selection_bit(o2::aod::evsel::kIsBBFDC))
        SETBIT(fitInfo.BBFDDCpf, bit);
      ++curit;
      if (curit == v.end())
        break;
      curbc = curit->first;
    }
  }

  template <int32_t tracksSwitch, typename TAmbTrack>
  int64_t getAmbTrackId(TAmbTrack ambTrack)
  {
    int64_t trkId = -1;
    if constexpr (tracksSwitch == 0) { // central barrel
      trkId = ambTrack.trackId();
    }
    if constexpr (tracksSwitch == 1) { // forward tracks
      trkId = ambTrack.fwdtrackId();
    }
    return trkId;
  }

  // "uncorrected" bcs
  template <int32_t tracksSwitch, typename TBCs, typename TAmbTracks>
  void collectAmbTrackBCs(std::unordered_map<int64_t, uint64_t>& ambTrIds,
                          TAmbTracks ambTracks)
  {
    for (const auto& ambTrk : ambTracks) {
      auto trkId = getAmbTrackId<tracksSwitch>(ambTrk);
      const auto& bcSlice = ambTrk.template bc_as<TBCs>();
      uint64_t trackBC = -1;
      if (bcSlice.size() != 0) {
        auto first = bcSlice.begin();
        trackBC = first.globalBC();
      }
      ambTrIds[trkId] = trackBC;
    }
  }

  void addTrack(std::vector<BCTracksPair>& v, uint64_t bc, int64_t trkId)
  {
    auto it = std::find_if(v.begin(), v.end(),
                           [bc](const BCTracksPair& element) { return element.first == bc; });
    if (it != v.end())
      it->second.push_back(trkId);
    else
      v.emplace_back(std::make_pair(bc, std::vector<int64_t>({trkId})));
  }

  // trackType == 0 -> hasTOF
  // trackType == 1 -> hasITS and not hasTOF
  template <typename TBCs>
  void collectBarrelTracks(std::vector<BCTracksPair>& bcsMatchedTrIds,
                           int trackType,
                           TBCs const& /*bcs*/,
                           o2::aod::Collisions const& /*collisions*/,
                           BarrelTracks const& barrelTracks,
                           o2::aod::AmbiguousTracks const& /*ambBarrelTracks*/,
                           std::unordered_map<int64_t, uint64_t>& ambBarrelTrBCs)
  {
    for (const auto& trk : barrelTracks) {
      if (!trk.hasTPC())
        continue;
      if (trackType == 0 && !trk.hasTOF())
        continue;
      if (trackType == 1 && !(trk.hasITS() && !trk.hasTOF()))
        continue;
      if (!applyBarCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int32_t nContrib = -1;
      uint64_t trackBC = 0;
      if (trk.has_collision()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
        trackBC = col.bc_as<TBCs>().globalBC();
      } else {
        auto ambIter = ambBarrelTrBCs.find(trkId);
        if (ambIter != ambBarrelTrBCs.end())
          trackBC = ambIter->second;
      }
      int64_t tint = TMath::FloorNint(trk.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + static_cast<float>(fBarrelTrackTShift));
      uint64_t bc = trackBC + tint;
      if (nContrib > upcCuts.getMaxNContrib())
        continue;
      addTrack(bcsMatchedTrIds, bc, trkId);
    }
  }

  template <typename TBCs>
  void collectForwardTracks(std::vector<BCTracksPair>& bcsMatchedTrIds,
                            int typeFilter,
                            TBCs const& /*bcs*/,
                            o2::aod::Collisions const& /*collisions*/,
                            ForwardTracks const& fwdTracks,
                            o2::aod::AmbiguousFwdTracks const& /*ambFwdTracks*/,
                            std::unordered_map<int64_t, uint64_t>& ambFwdTrBCs)
  {
    for (const auto& trk : fwdTracks) {
      if (trk.trackType() != typeFilter)
        continue;
      if (!applyFwdCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int32_t nContrib = -1;
      uint64_t trackBC = 0;
      auto ambIter = ambFwdTrBCs.find(trkId);
      if (ambIter == ambFwdTrBCs.end()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
        trackBC = col.bc_as<TBCs>().globalBC();
      } else {
        trackBC = ambIter->second;
      }
      int64_t tint = TMath::FloorNint(trk.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + static_cast<float>(fMuonTrackTShift));
      uint64_t bc = trackBC + tint;
      if (nContrib > upcCuts.getMaxNContrib())
        continue;
      addTrack(bcsMatchedTrIds, bc, trkId);
    }
  }

  int32_t searchTracks(uint64_t midbc, uint64_t range, uint32_t tracksToFind,
                       std::vector<int64_t>& tracks,
                       std::vector<BCTracksPair>& v,
                       std::unordered_set<int64_t>& matchedTracks,
                       bool skipMidBC = false)
  {
    uint32_t count = 0;
    uint64_t left = midbc >= range ? midbc - range : 0;
    uint64_t right = fMaxBC >= midbc + range ? midbc + range : fMaxBC;
    BCTracksPair dummyPair(left, {});
    auto curit = std::lower_bound(v.begin(), v.end(), dummyPair,
                                  [](const BCTracksPair& left, const BCTracksPair& right) { return left.first < right.first; });
    if (curit == v.end()) // no ITS-TPC tracks nearby at all -> near last BCs
      return -1;
    uint64_t curbc = curit->first;
    while (curbc <= right) { // moving forward to midbc+range
      if (skipMidBC && curbc == midbc) {
        ++curit;
        if (curit == v.end())
          break;
        curbc = curit->first;
      }
      uint32_t size = curit->second.size();
      if (size > 1) // too many tracks per BC -> possibly another event
        return -2;
      count += size;
      if (count > tracksToFind) // too many tracks nearby
        return -3;
      if (matchedTracks.find(curit->second[0]) == matchedTracks.end()) {
        tracks.push_back(curit->second[0]);
        matchedTracks.insert(curit->second[0]);
      }
      ++curit;
      if (curit == v.end())
        break;
      curbc = curit->first;
    }
    if (count != tracksToFind)
      return -4;
    return 0; // found exactly tracksToFind tracks in [midbc - range, midbc + range]
  }

  template <typename TBCs>
  void createCandidatesCentral(BarrelTracks const& barrelTracks,
                               o2::aod::AmbiguousTracks const& ambBarrelTracks,
                               TBCs const& bcs,
                               o2::aod::Collisions const& collisions,
                               o2::aod::FT0s const& ft0s,
                               o2::aod::FDDs const& /*fdds*/,
                               o2::aod::FV0As const& fv0as,
                               o2::aod::Zdcs const& zdcs,
                               const o2::aod::McTrackLabels* mcBarrelTrackLabels)
  {
    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsTOF;
    std::vector<BCTracksPair> bcsMatchedTrIdsITSTPC;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, uint64_t> ambBarrelTrBCs;
    if (upcCuts.getAmbigSwitch() != 1)
      collectAmbTrackBCs<0, BCsWithBcSels>(ambBarrelTrBCs, ambBarrelTracks);

    collectBarrelTracks(bcsMatchedTrIdsTOF,
                        0,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrBCs);

    collectBarrelTracks(bcsMatchedTrIdsITSTPC,
                        1,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrBCs);

    std::sort(bcsMatchedTrIdsTOF.begin(), bcsMatchedTrIdsTOF.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::sort(bcsMatchedTrIdsITSTPC.begin(), bcsMatchedTrIdsITSTPC.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::map<uint64_t, int32_t> mapGlobalBcWithTOR{};
    std::map<uint64_t, int32_t> mapGlobalBcWithTVX{};
    std::map<uint64_t, int32_t> mapGlobalBcWithTSC{};
    for (const auto& ft0 : ft0s) {
      uint64_t globalBC = ft0.bc_as<TBCs>().globalBC();
      int32_t globalIndex = ft0.globalIndex();
      if (!(std::abs(ft0.timeA()) > 2.f && std::abs(ft0.timeC()) > 2.f))
        mapGlobalBcWithTOR[globalBC] = globalIndex;
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex)) { // TVX
        mapGlobalBcWithTVX[globalBC] = globalIndex;
      }
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen)) { // TVX & TCE
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("TCE", 1);
      }
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex) &&
          (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen) ||
           TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitSCen))) { // TVX & (TSC | TCE)
        mapGlobalBcWithTSC[globalBC] = globalIndex;
      }
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithV0A{};
    for (const auto& fv0a : fv0as) {
      if (std::abs(fv0a.time()) > 15.f)
        continue;
      uint64_t globalBC = fv0a.bc_as<TBCs>().globalBC();
      mapGlobalBcWithV0A[globalBC] = fv0a.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithZdc{};
    for (const auto& zdc : zdcs) {
      if (std::abs(zdc.timeZNA()) > 2.f && std::abs(zdc.timeZNC()) > 2.f)
        continue;
      if (!(std::abs(zdc.timeZNA()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNA", 1);
      if (!(std::abs(zdc.timeZNC()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNC", 1);
      auto globalBC = zdc.bc_as<TBCs>().globalBC();
      mapGlobalBcWithZdc[globalBC] = zdc.globalIndex();
    }

    auto nTORs = mapGlobalBcWithTOR.size();
    auto nTSCs = mapGlobalBcWithTSC.size();
    auto nTVXs = mapGlobalBcWithTVX.size();
    auto nFV0As = mapGlobalBcWithV0A.size();
    auto nZdcs = mapGlobalBcWithZdc.size();
    auto nBcsWithITSTPC = bcsMatchedTrIdsITSTPC.size();

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    auto updateFitInfo = [&](uint64_t globalBC, upchelpers::FITInfo& fitInfo) {
      fitInfo.timeFT0A = -999.f;
      fitInfo.timeFT0C = -999.f;
      fitInfo.timeFV0A = -999.f;
      fitInfo.ampFT0A = 0.f;
      fitInfo.ampFT0C = 0.f;
      fitInfo.ampFV0A = 0.f;
      fitInfo.BBFT0Apf = -999;
      fitInfo.BBFV0Apf = -999;
      fitInfo.distClosestBcTOR = 999;
      fitInfo.distClosestBcTSC = 999;
      fitInfo.distClosestBcTVX = 999;
      fitInfo.distClosestBcV0A = 999;
      if (nTORs > 0) {
        uint64_t closestBcTOR = findClosestBC(globalBC, mapGlobalBcWithTOR);
        fitInfo.distClosestBcTOR = globalBC - static_cast<int64_t>(closestBcTOR);
        if (std::abs(fitInfo.distClosestBcTOR) <= fFilterFT0)
          return false;
        auto ft0Id = mapGlobalBcWithTOR.at(closestBcTOR);
        auto ft0 = ft0s.iteratorAt(ft0Id);
        fitInfo.timeFT0A = ft0.timeA();
        fitInfo.timeFT0C = ft0.timeC();
        const auto& t0AmpsA = ft0.amplitudeA();
        const auto& t0AmpsC = ft0.amplitudeC();
        for (auto amp : t0AmpsA)
          fitInfo.ampFT0A += amp;
        for (auto amp : t0AmpsC)
          fitInfo.ampFT0C += amp;
      }
      if (nTSCs > 0) {
        uint64_t closestBcTSC = findClosestBC(globalBC, mapGlobalBcWithTSC);
        fitInfo.distClosestBcTSC = globalBC - static_cast<int64_t>(closestBcTSC);
        if (std::abs(fitInfo.distClosestBcTSC) <= fFilterTSC)
          return false;
      }
      if (nTVXs > 0) {
        uint64_t closestBcTVX = findClosestBC(globalBC, mapGlobalBcWithTVX);
        fitInfo.distClosestBcTVX = globalBC - static_cast<int64_t>(closestBcTVX);
        if (std::abs(fitInfo.distClosestBcTVX) <= fFilterTVX)
          return false;
      }
      if (nFV0As > 0) {
        uint64_t closestBcV0A = findClosestBC(globalBC, mapGlobalBcWithV0A);
        fitInfo.distClosestBcV0A = globalBC - static_cast<int64_t>(closestBcV0A);
        if (std::abs(fitInfo.distClosestBcV0A) <= fFilterFV0)
          return false;
        auto fv0aId = mapGlobalBcWithV0A.at(closestBcV0A);
        auto fv0a = fv0as.iteratorAt(fv0aId);
        fitInfo.timeFV0A = fv0a.time();
        const auto& v0Amps = fv0a.amplitude();
        for (auto amp : v0Amps)
          fitInfo.ampFV0A += amp;
      }
      return true;
    };

    // candidates with TOF
    int32_t candID = 0;
    for (auto& pair : bcsMatchedTrIdsTOF) {
      auto globalBC = pair.first;
      auto& barrelTrackIDs = pair.second;
      uint32_t nTOFs = barrelTrackIDs.size();
      if (nTOFs > fNBarProngs) // too many tracks
        continue;
      auto closestBcITSTPC = std::numeric_limits<uint64_t>::max();
      if (nTOFs < fNBarProngs && nBcsWithITSTPC > 0) { // adding ITS-TPC tracks
        auto itClosestBcITSTPC = findClosestTrackBCiter(globalBC, bcsMatchedTrIdsITSTPC);
        if (itClosestBcITSTPC == bcsMatchedTrIdsITSTPC.end())
          continue;
        closestBcITSTPC = itClosestBcITSTPC->first;
        int64_t distClosestBcITSTPC = globalBC - static_cast<int64_t>(closestBcITSTPC);
        histRegistry.fill(HIST("hDistToITSTPC"), std::abs(distClosestBcITSTPC));
        if (std::abs(distClosestBcITSTPC) > fBcWindowITSTPC)
          continue;
        auto& itstpcTracks = itClosestBcITSTPC->second;
        uint32_t nITSTPCs = itstpcTracks.size();
        if ((nTOFs + nITSTPCs) != fNBarProngs)
          continue;
        barrelTrackIDs.insert(barrelTrackIDs.end(), itstpcTracks.begin(), itstpcTracks.end());
        itClosestBcITSTPC->second.clear(); // BC is matched to BC with TOF, removing tracks, but leaving BC
      }
      upchelpers::FITInfo fitInfo{};
      if (!updateFitInfo(globalBC, fitInfo))
        continue;
      if (nZdcs > 0) {
        auto itZDC = mapGlobalBcWithZdc.find(globalBC);
        if (itZDC != mapGlobalBcWithZdc.end()) {
          const auto& zdc = zdcs.iteratorAt(itZDC->second);
          float timeZNA = zdc.timeZNA();
          float timeZNC = zdc.timeZNC();
          float eComZNA = zdc.energyCommonZNA();
          float eComZNC = zdc.energyCommonZNC();
          udZdcsReduced(candID, timeZNA, timeZNC, eComZNA, eComZNC);
        }
      }
      uint16_t numContrib = fNBarProngs;
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      for (auto id : barrelTrackIDs) {
        const auto& tr = barrelTracks.iteratorAt(id);
        netCharge += tr.sign();
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = RgtrwTOF / static_cast<float>(numContrib);
      // store used tracks
      int upc_flag = 0;
      // TODO: introduce better check on association of collision and reconstruction mode
      if (bcs.iteratorAt(0).flags() == o2::itsmft::ROFRecord::VtxUPCMode)
        upc_flag = 1;
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, globalBC, closestBcITSTPC, mcBarrelTrackLabels, ambBarrelTrBCs);
      eventCandidates(globalBC, runNumber, dummyX, dummyY, dummyZ, upc_flag, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      eventCandidatesSelsCent(fitInfo.distClosestBcTOR,
                              fitInfo.distClosestBcTSC,
                              fitInfo.distClosestBcTVX,
                              fitInfo.distClosestBcV0A);
      candID++;
    }

    // candidates without TOF
    for (auto& pair : bcsMatchedTrIdsITSTPC) {
      auto globalBC = pair.first;
      auto& barrelTrackIDs = pair.second;
      uint32_t nThisITSTPCs = barrelTrackIDs.size();
      if (nThisITSTPCs > fNBarProngs || nThisITSTPCs == 0) // too many tracks / already matched to TOF
        continue;
      auto closestBcITSTPC = std::numeric_limits<uint64_t>::max();
      if (nThisITSTPCs < fNBarProngs) { // adding ITS-TPC tracks
        auto itClosestBcITSTPC = findClosestTrackBCiterNotEq(globalBC, bcsMatchedTrIdsITSTPC);
        if (itClosestBcITSTPC == bcsMatchedTrIdsITSTPC.end())
          continue;
        closestBcITSTPC = itClosestBcITSTPC->first;
        int64_t distClosestBcITSTPC = globalBC - static_cast<int64_t>(closestBcITSTPC);
        histRegistry.fill(HIST("hDistToITSTPC"), std::abs(distClosestBcITSTPC));
        if (std::abs(distClosestBcITSTPC) > fBcWindowITSTPC)
          continue;
        auto& itstpcTracks = itClosestBcITSTPC->second;
        uint32_t nITSTPCs = itstpcTracks.size();
        if ((nThisITSTPCs + nITSTPCs) != fNBarProngs)
          continue;
        barrelTrackIDs.insert(barrelTrackIDs.end(), itstpcTracks.begin(), itstpcTracks.end());
        itClosestBcITSTPC->second.clear();
      }
      upchelpers::FITInfo fitInfo{};
      if (!updateFitInfo(globalBC, fitInfo))
        continue;
      if (nZdcs > 0) {
        auto itZDC = mapGlobalBcWithZdc.find(globalBC);
        if (itZDC != mapGlobalBcWithZdc.end()) {
          const auto& zdc = zdcs.iteratorAt(itZDC->second);
          float timeZNA = zdc.timeZNA();
          float timeZNC = zdc.timeZNC();
          float eComZNA = zdc.energyCommonZNA();
          float eComZNC = zdc.energyCommonZNC();
          udZdcsReduced(candID, timeZNA, timeZNC, eComZNA, eComZNC);
        }
      }
      uint16_t numContrib = fNBarProngs;
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      for (auto id : barrelTrackIDs) {
        const auto& tr = barrelTracks.iteratorAt(id);
        netCharge += tr.sign();
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = RgtrwTOF / static_cast<float>(numContrib);
      // store used tracks
      int upc_flag = 0;
      // TODO: introduce better check on association of collision and reconstruction mode
      if (bcs.iteratorAt(0).flags() == o2::itsmft::ROFRecord::VtxUPCMode)
        upc_flag = 1;
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, globalBC, closestBcITSTPC, mcBarrelTrackLabels, ambBarrelTrBCs);
      eventCandidates(globalBC, runNumber, dummyX, dummyY, dummyZ, upc_flag, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      eventCandidatesSelsCent(fitInfo.distClosestBcTOR,
                              fitInfo.distClosestBcTSC,
                              fitInfo.distClosestBcTVX,
                              fitInfo.distClosestBcV0A);
      barrelTrackIDs.clear();
      candID++;
    }

    ambBarrelTrBCs.clear();
    bcsMatchedTrIdsITSTPC.clear();
    bcsMatchedTrIdsTOF.clear();
  }

  template <typename TBCs>
  void createCandidatesSemiFwd(BarrelTracks const& barrelTracks,
                               o2::aod::AmbiguousTracks const& ambBarrelTracks,
                               ForwardTracks const& fwdTracks,
                               o2::aod::FwdTrkCls const& /*fwdTrkClusters*/,
                               o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                               TBCs const& bcs,
                               o2::aod::Collisions const& collisions,
                               o2::aod::FT0s const& ft0s,
                               o2::aod::FDDs const& fdds,
                               o2::aod::FV0As const& fv0as,
                               const o2::aod::McTrackLabels* mcBarrelTrackLabels,
                               const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    LOGP(debug, "barrelTracks.size()={}", barrelTracks.size());
    LOGP(debug, "ambBarrelTracks.size()={}", ambBarrelTracks.size());
    LOGP(debug, "fwdTracks.size()={}", fwdTracks.size());
    LOGP(debug, "ambFwdTracks.size()={}", ambFwdTracks.size());

    fMaxBC = bcs.iteratorAt(bcs.size() - 1).globalBC(); // restrict ITS-TPC track search to [0, fMaxBC]

    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsTOF;
    std::vector<BCTracksPair> bcsMatchedTrIdsITSTPC;
    std::vector<BCTracksPair> bcsMatchedTrIdsMID;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, uint64_t> ambBarrelTrBCs;
    collectAmbTrackBCs<0, BCsWithBcSels>(ambBarrelTrBCs, ambBarrelTracks);

    std::unordered_map<int64_t, uint64_t> ambFwdTrBCs;
    collectAmbTrackBCs<1, BCsWithBcSels>(ambFwdTrBCs, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    collectBarrelTracks(bcsMatchedTrIdsTOF,
                        0,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrBCs);

    collectBarrelTracks(bcsMatchedTrIdsITSTPC,
                        1,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrBCs);

    LOGP(debug, "bcsMatchedTrIdsMID.size()={}", bcsMatchedTrIdsMID.size());
    LOGP(debug, "bcsMatchedTrIdsTOF.size()={}", bcsMatchedTrIdsTOF.size());

    uint32_t nBCsWithITSTPC = bcsMatchedTrIdsITSTPC.size();
    uint32_t nBCsWithMID = bcsMatchedTrIdsMID.size();

    std::sort(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::vector<BCTracksPair> bcsMatchedTrIdsTOFTagged(nBCsWithMID);
    for (const auto& pair : bcsMatchedTrIdsTOF) {
      uint64_t bc = pair.first;
      auto it = std::find_if(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
                             [bc](const auto& item) { return item.first == bc; });
      if (it != bcsMatchedTrIdsMID.end()) {
        uint32_t ibc = it - bcsMatchedTrIdsMID.begin();
        bcsMatchedTrIdsTOFTagged[ibc].second = pair.second;
      }
    }

    bcsMatchedTrIdsTOF.clear();

    std::sort(bcsMatchedTrIdsITSTPC.begin(), bcsMatchedTrIdsITSTPC.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    if (nBCsWithITSTPC > 0 && fSearchITSTPC == 1) {
      std::unordered_set<int64_t> matchedTracks;
      for (uint32_t ibc = 0; ibc < nBCsWithMID; ++ibc) {
        auto& pairMID = bcsMatchedTrIdsMID[ibc];
        auto& pairTOF = bcsMatchedTrIdsTOFTagged[ibc];
        uint64_t bc = pairMID.first;
        auto& trackIdsMID = pairMID.second;
        auto& trackIdsTOF = pairTOF.second;
        uint32_t nMIDtracks = trackIdsMID.size();
        uint32_t nTOFtracks = trackIdsTOF.size();
        if (nMIDtracks > fNFwdProngs || nTOFtracks > fNBarProngs) // too many MID and/or TOF tracks?!
          continue;
        if (nMIDtracks == fNFwdProngs && nTOFtracks == fNBarProngs) { // check for ITS-TPC tracks
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, 0, tracks, bcsMatchedTrIdsITSTPC, matchedTracks, false);
          if (res < 0) { // too many tracks nearby -> rejecting
            trackIdsTOF.push_back(0);
            continue;
          }
        }
        if (nMIDtracks == fNFwdProngs && nTOFtracks < fNBarProngs && !upcCuts.getRequireTOF()) { // add ITS-TPC track if needed
          uint32_t tracksToFind = fNBarProngs - nTOFtracks;
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, tracksToFind, tracks, bcsMatchedTrIdsITSTPC, matchedTracks, true);
          if (res < 0) // too many or not enough tracks nearby -> rejecting
            continue;
          trackIdsTOF.insert(trackIdsTOF.end(), tracks.begin(), tracks.end());
        }
      }
    }

    bcsMatchedTrIdsITSTPC.clear();

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    std::vector<std::pair<uint64_t, int64_t>> indexBCglId;
    indexBCglId.reserve(bcs.size());
    for (const auto& bc : bcs) {
      if (bc.has_foundFT0() || bc.has_foundFV0() || bc.has_foundFDD())
        indexBCglId.emplace_back(std::make_pair(bc.globalBC(), bc.globalIndex()));
    }

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    // storing n-prong matches
    int32_t candID = 0;
    for (uint32_t ibc = 0; ibc < nBCsWithMID; ++ibc) {
      auto& pairMID = bcsMatchedTrIdsMID[ibc];
      auto& pairTOF = bcsMatchedTrIdsTOFTagged[ibc];
      auto& fwdTrackIDs = pairMID.second;
      auto& barrelTrackIDs = pairTOF.second;
      uint32_t nMIDtracks = fwdTrackIDs.size();
      uint32_t nBarrelTracks = barrelTrackIDs.size(); // TOF + ITS-TPC tracks
      uint16_t numContrib = nBarrelTracks + nMIDtracks;
      uint64_t bc = pairMID.first;
      // sanity check
      if (nBarrelTracks != fNBarProngs || nMIDtracks != fNFwdProngs) {
        continue;
      }
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      upchelpers::FITInfo fitInfo{};
      processFITInfo(fitInfo, bc, indexBCglId, bcs, ft0s, fdds, fv0as);
      if (fFilterFT0) {
        if (!checkFT0(fitInfo, false))
          continue;
      }
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      for (auto id : barrelTrackIDs) {
        const auto& tr = barrelTracks.iteratorAt(id);
        netCharge += tr.sign();
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = RgtrwTOF / static_cast<float>(numContrib);
      // store used tracks
      int upc_flag = 0;
      // TODO: introduce better check on association of collision and reconstruction mode
      if (bcs.iteratorAt(0).flags() == o2::itsmft::ROFRecord::VtxUPCMode)
        upc_flag = 1;
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, bc, mcFwdTrackLabels);
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, bc, mcBarrelTrackLabels, ambBarrelTrBCs);
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, upc_flag, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      candID++;
    }

    indexBCglId.clear();
    ambFwdTrBCs.clear();
    bcsMatchedTrIdsMID.clear();
    ambBarrelTrBCs.clear();
    bcsMatchedTrIdsTOFTagged.clear();
  }

  template <typename T>
  void fillAmplitudes(const T& t,
                      const std::map<uint64_t, int32_t>& mapBCs,
                      std::vector<float>& amps,
                      std::vector<int8_t>& relBCs,
                      uint64_t gbc)
  {
    auto s = gbc - fBCWindowFITAmps;
    auto e = gbc + (fBCWindowFITAmps - 1);
    auto it = mapBCs.lower_bound(s);
    while (it->first <= e && it != mapBCs.end()) {
      int i = it->first - s;
      auto id = it->second;
      const auto& row = t.iteratorAt(id);
      float totalAmp = 0.f;
      if constexpr (std::is_same_v<T, o2::aod::FT0s>) {
        const auto& itAmps = row.amplitudeA();
        totalAmp = std::accumulate(itAmps.begin(), itAmps.end(), 0.f);
      }
      if constexpr (std::is_same_v<T, o2::aod::FV0As>) {
        const auto& itAmps = row.amplitude();
        totalAmp = std::accumulate(itAmps.begin(), itAmps.end(), 0.f);
      }
      if (totalAmp > 0.f) {
        amps.push_back(totalAmp);
        relBCs.push_back(gbc - (i + s));
      }
      ++it;
    }
  }

  template <typename TBCs>
  void createCandidatesFwd(ForwardTracks const& fwdTracks,
                           o2::aod::FwdTrkCls const& fwdTrkClusters,
                           o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                           TBCs const& bcs,
                           o2::aod::Collisions const& collisions,
                           o2::aod::FT0s const& ft0s,
                           o2::aod::FDDs const& fdds,
                           o2::aod::FV0As const& fv0as,
                           o2::aod::Zdcs const& zdcs,
                           const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsMID;
    std::vector<BCTracksPair> bcsMatchedTrIdsMCH;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, uint64_t> ambFwdTrBCs;
    collectAmbTrackBCs<1, BCsWithBcSels>(ambFwdTrBCs, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    collectForwardTracks(bcsMatchedTrIdsMCH,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    std::sort(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::sort(bcsMatchedTrIdsMCH.begin(), bcsMatchedTrIdsMCH.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::map<uint64_t, int32_t> mapGlobalBcWithT0A{};
    for (const auto& ft0 : ft0s) {
      if (!TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex))
        continue;
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen)) { // TVX & TCE
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("TCE", 1);
      }
      if (std::abs(ft0.timeA()) > 2.f)
        continue;
      uint64_t globalBC = ft0.bc_as<TBCs>().globalBC();
      mapGlobalBcWithT0A[globalBC] = ft0.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithV0A{};
    for (const auto& fv0a : fv0as) {
      if (!TESTBIT(fv0a.triggerMask(), o2::fit::Triggers::bitA))
        continue;
      if (std::abs(fv0a.time()) > 15.f)
        continue;
      uint64_t globalBC = fv0a.bc_as<TBCs>().globalBC();
      mapGlobalBcWithV0A[globalBC] = fv0a.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithZdc{};
    for (const auto& zdc : zdcs) {
      if (std::abs(zdc.timeZNA()) > 2.f && std::abs(zdc.timeZNC()) > 2.f)
        continue;
      if (!(std::abs(zdc.timeZNA()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNA", 1);
      if (!(std::abs(zdc.timeZNC()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNC", 1);
      auto globalBC = zdc.bc_as<TBCs>().globalBC();
      mapGlobalBcWithZdc[globalBC] = zdc.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithFDD{};
    uint8_t twoLayersA = 0;
    uint8_t twoLayersC = 0;
    for (const auto& fdd : fdds) {
      // get signal coincidence
      for (int i = 0; i < 4; i++) {
        if (fdd.chargeA()[i + 4] > 0 && fdd.chargeA()[i] > 0)
          twoLayersA++;
        if (fdd.chargeC()[i + 4] > 0 && fdd.chargeC()[i] > 0)
          twoLayersC++;
      }
      // if no signal, continue
      if ((twoLayersA == 0) && (twoLayersC == 0))
        continue;
      uint64_t globalBC = fdd.bc_as<TBCs>().globalBC();
      mapGlobalBcWithFDD[globalBC] = fdd.globalIndex();
    }

    auto nFT0s = mapGlobalBcWithT0A.size();
    auto nFV0As = mapGlobalBcWithV0A.size();
    auto nZdcs = mapGlobalBcWithZdc.size();
    auto nBcsWithMCH = bcsMatchedTrIdsMCH.size();
    auto nFDDs = mapGlobalBcWithFDD.size();

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    std::vector<int> selTrackIds{};

    // storing n-prong matches
    int32_t candID = 0;

    for (auto& pair : bcsMatchedTrIdsMID) { // candidates without MFT
      auto globalBC = static_cast<int64_t>(pair.first);
      const auto& fwdTrackIDs = pair.second; // only MID-matched tracks at the moment
      uint32_t nMIDs = fwdTrackIDs.size();
      if (nMIDs > fNFwdProngs) // too many tracks
        continue;
      std::vector<int64_t> trkCandIDs{};
      if (nMIDs == fNFwdProngs) {
        trkCandIDs.insert(trkCandIDs.end(), fwdTrackIDs.begin(), fwdTrackIDs.end());
      }
      uint64_t closestBcMCH = 0;
      if (nMIDs < fNFwdProngs && nBcsWithMCH > 0) { // adding MCH tracks
        auto itClosestBcMCH = findClosestTrackBCiter(globalBC, bcsMatchedTrIdsMCH);
        closestBcMCH = itClosestBcMCH->first;
        int64_t distClosestBcMCH = globalBC - static_cast<int64_t>(closestBcMCH);
        if (std::abs(distClosestBcMCH) > fBcWindowMCH)
          continue;
        auto& mchTracks = itClosestBcMCH->second;
        uint32_t nMCHs = mchTracks.size();
        if ((nMCHs + nMIDs) != fNFwdProngs)
          continue;
        trkCandIDs.insert(trkCandIDs.end(), fwdTrackIDs.begin(), fwdTrackIDs.end());
        trkCandIDs.insert(trkCandIDs.end(), mchTracks.begin(), mchTracks.end());
      }
      upchelpers::FITInfo fitInfo{};
      fitInfo.timeFT0A = -999.f;
      fitInfo.timeFT0C = -999.f;
      fitInfo.timeFV0A = -999.f;
      fitInfo.ampFT0A = 0.f;
      fitInfo.ampFT0C = 0.f;
      fitInfo.ampFV0A = 0.f;
      std::vector<float> amplitudesT0A{};
      std::vector<float> amplitudesV0A{};
      std::vector<int8_t> relBCsT0A{};
      std::vector<int8_t> relBCsV0A{};
      uint8_t chFT0A = 0;
      uint8_t chFT0C = 0;
      if (nFT0s > 0) {
        uint64_t closestBcT0A = findClosestBC(globalBC, mapGlobalBcWithT0A);
        int64_t distClosestBcT0A = globalBC - static_cast<int64_t>(closestBcT0A);
        if (std::abs(distClosestBcT0A) <= fFilterFT0)
          continue;
        fitInfo.distClosestBcT0A = distClosestBcT0A;
        auto ft0Id = mapGlobalBcWithT0A.at(closestBcT0A);
        auto ft0 = ft0s.iteratorAt(ft0Id);
        fitInfo.timeFT0A = ft0.timeA();
        fitInfo.timeFT0C = ft0.timeC();
        const auto& t0AmpsA = ft0.amplitudeA();
        const auto& t0AmpsC = ft0.amplitudeC();
        fitInfo.ampFT0A = std::accumulate(t0AmpsA.begin(), t0AmpsA.end(), 0.f);
        fitInfo.ampFT0C = std::accumulate(t0AmpsC.begin(), t0AmpsC.end(), 0.f);
        chFT0A = ft0.amplitudeA().size();
        chFT0C = ft0.amplitudeC().size();
        fillAmplitudes(ft0s, mapGlobalBcWithT0A, amplitudesT0A, relBCsT0A, globalBC);
      }
      uint8_t chFV0A = 0;
      if (nFV0As > 0) {
        uint64_t closestBcV0A = findClosestBC(globalBC, mapGlobalBcWithV0A);
        int64_t distClosestBcV0A = globalBC - static_cast<int64_t>(closestBcV0A);
        if (std::abs(distClosestBcV0A) <= fFilterFV0)
          continue;
        fitInfo.distClosestBcV0A = distClosestBcV0A;
        auto fv0aId = mapGlobalBcWithV0A.at(closestBcV0A);
        auto fv0a = fv0as.iteratorAt(fv0aId);
        fitInfo.timeFV0A = fv0a.time();
        const auto& v0Amps = fv0a.amplitude();
        fitInfo.ampFV0A = std::accumulate(v0Amps.begin(), v0Amps.end(), 0.f);
        chFV0A = fv0a.amplitude().size();
        fillAmplitudes(fv0as, mapGlobalBcWithV0A, amplitudesV0A, relBCsV0A, globalBC);
      }
      uint8_t chFDDA = 0;
      uint8_t chFDDC = 0;
      if (nFDDs > 0) {
        uint64_t closestBcFDD = findClosestBC(globalBC, mapGlobalBcWithFDD);
        auto fddId = mapGlobalBcWithFDD.at(closestBcFDD);
        auto fdd = fdds.iteratorAt(fddId);
        fitInfo.timeFDDA = fdd.timeA();
        fitInfo.timeFDDC = fdd.timeC();
        fitInfo.ampFDDA = 0;
        for (int i = 0; i < 8; i++)
          fitInfo.ampFDDA += fdd.chargeA()[i];
        fitInfo.ampFDDC = 0;
        for (int i = 0; i < 8; i++)
          fitInfo.ampFDDC += fdd.chargeC()[i];
        fitInfo.triggerMaskFDD = fdd.triggerMask();
        // get signal coincidence
        for (int i = 0; i < 4; i++) {
          if (fdd.chargeA()[i + 4] > 0 && fdd.chargeA()[i] > 0)
            chFDDA++;
          if (fdd.chargeC()[i + 4] > 0 && fdd.chargeC()[i] > 0)
            chFDDC++;
        }
      }
      if (nZdcs > 0) {
        auto itZDC = mapGlobalBcWithZdc.find(globalBC);
        if (itZDC != mapGlobalBcWithZdc.end()) {
          const auto& zdc = zdcs.iteratorAt(itZDC->second);
          float timeZNA = zdc.timeZNA();
          float timeZNC = zdc.timeZNC();
          float eComZNA = zdc.energyCommonZNA();
          float eComZNC = zdc.energyCommonZNC();
          udZdcsReduced(candID, timeZNA, timeZNC, eComZNA, eComZNC);
        }
      }
      uint16_t numContrib = fNFwdProngs;
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      for (auto id : trkCandIDs) {
        auto tr = fwdTracks.iteratorAt(id);
        netCharge += tr.sign();
        selTrackIds.push_back(id);
      }
      // store used tracks
      int upc_flag = 0;
      // TODO: introduce better check on association of collision and reconstruction mode
      if (bcs.iteratorAt(0).flags() == o2::itsmft::ROFRecord::VtxUPCMode)
        upc_flag = 1;
      fillFwdTracks(fwdTracks, trkCandIDs, candID, globalBC, closestBcMCH, mcFwdTrackLabels);
      eventCandidates(globalBC, runNumber, dummyX, dummyY, dummyZ, upc_flag, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      eventCandidatesSelExtras(chFT0A, chFT0C, chFDDA, chFDDC, chFV0A, 0, 0, 0, 0, 0);
      eventCandidatesSelsFwd(fitInfo.distClosestBcV0A,
                             fitInfo.distClosestBcT0A,
                             amplitudesT0A,
                             relBCsT0A,
                             amplitudesV0A,
                             relBCsV0A);
      candID++;
      trkCandIDs.clear();
    }

    fillFwdClusters(selTrackIds, fwdTrkClusters);

    selTrackIds.clear();
    ambFwdTrBCs.clear();
    bcsMatchedTrIdsMID.clear();
    bcsMatchedTrIdsMCH.clear();
    mapGlobalBcWithT0A.clear();
    mapGlobalBcWithV0A.clear();
    mapGlobalBcWithFDD.clear();
  }

  template <typename TBCs>
  void createCandidatesFwdGlobal(ForwardTracks const& fwdTracks,
                                 o2::aod::FwdTrkCls const& /*fwdTrkClusters*/,
                                 o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                                 TBCs const& bcs,
                                 o2::aod::Collisions const& collisions,
                                 o2::aod::FT0s const& ft0s,
                                 o2::aod::FDDs const& fdds,
                                 o2::aod::FV0As const& fv0as,
                                 o2::aod::Zdcs const& zdcs,
                                 const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsMID;
    std::vector<BCTracksPair> bcsMatchedTrIdsMCH;
    std::vector<BCTracksPair> bcsMatchedTrIdsGlobal;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, uint64_t> ambFwdTrBCs;
    collectAmbTrackBCs<1, BCsWithBcSels>(ambFwdTrBCs, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    collectForwardTracks(bcsMatchedTrIdsMCH,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    collectForwardTracks(bcsMatchedTrIdsGlobal,
                         o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    std::sort(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::sort(bcsMatchedTrIdsMCH.begin(), bcsMatchedTrIdsMCH.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::sort(bcsMatchedTrIdsGlobal.begin(), bcsMatchedTrIdsGlobal.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    std::map<uint64_t, int32_t> mapGlobalBcWithT0A{};
    for (const auto& ft0 : ft0s) {
      if (!TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitVertex))
        continue;
      if (TESTBIT(ft0.triggerMask(), o2::fit::Triggers::bitCen)) { // TVX & TCE
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("TCE", 1);
      }
      if (std::abs(ft0.timeA()) > 2.f)
        continue;
      uint64_t globalBC = ft0.bc_as<TBCs>().globalBC();
      mapGlobalBcWithT0A[globalBC] = ft0.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithV0A{};
    for (const auto& fv0a : fv0as) {
      if (!TESTBIT(fv0a.triggerMask(), o2::fit::Triggers::bitA))
        continue;
      if (std::abs(fv0a.time()) > 15.f)
        continue;
      uint64_t globalBC = fv0a.bc_as<TBCs>().globalBC();
      mapGlobalBcWithV0A[globalBC] = fv0a.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithZdc{};
    for (const auto& zdc : zdcs) {
      if (std::abs(zdc.timeZNA()) > 2.f && std::abs(zdc.timeZNC()) > 2.f)
        continue;
      if (!(std::abs(zdc.timeZNA()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNA", 1);
      if (!(std::abs(zdc.timeZNC()) > 2.f))
        histRegistry.get<TH1>(HIST("hCountersTrg"))->Fill("ZNC", 1);
      auto globalBC = zdc.bc_as<TBCs>().globalBC();
      mapGlobalBcWithZdc[globalBC] = zdc.globalIndex();
    }

    std::map<uint64_t, int32_t> mapGlobalBcWithFDD{};
    uint8_t twoLayersA = 0;
    uint8_t twoLayersC = 0;
    for (const auto& fdd : fdds) {
      // get signal coincidence
      for (int i = 0; i < 4; i++) {
        if (fdd.chargeA()[i + 4] > 0 && fdd.chargeA()[i] > 0)
          twoLayersA++;
        if (fdd.chargeC()[i + 4] > 0 && fdd.chargeC()[i] > 0)
          twoLayersC++;
      }
      // if no signal, continue
      if ((twoLayersA == 0) && (twoLayersC == 0))
        continue;
      uint64_t globalBC = fdd.bc_as<TBCs>().globalBC();
      mapGlobalBcWithFDD[globalBC] = fdd.globalIndex();
    }

    auto nFT0s = mapGlobalBcWithT0A.size();
    auto nFV0As = mapGlobalBcWithV0A.size();
    auto nZdcs = mapGlobalBcWithZdc.size();
    auto nFDDs = mapGlobalBcWithFDD.size();

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    std::vector<int> selTrackIdsGlobal{};

    int32_t candID = 0;

    for (auto& pair : bcsMatchedTrIdsGlobal) {
      auto globalBC = static_cast<int64_t>(pair.first);
      const auto& fwdTrackIDs = pair.second; // Forward tracks (Global with MFT)
      if (fwdTrackIDs.size() != 2) {         // ensure we have two MFT tracks
        continue;
      }

      // find the corresponding MCH-MID tracks
      auto midIt = std::find_if(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
                                [globalBC](const auto& midPair) {
                                  return midPair.first == static_cast<uint64_t>(globalBC);
                                });
      const auto* midTrackIDs = (midIt != bcsMatchedTrIdsMID.end()) ? &midIt->second : nullptr;

      // ensure MCH-MID tracks are available
      if (!midTrackIDs || midTrackIDs->size() != 2) {
        continue;
      }

      std::vector<int64_t> trkCandIDs;

      // retrieve global track eta and apply the logic
      bool firstInAcceptance = false, secondInAcceptance = false;

      const auto& trk1 = fwdTracks.iteratorAt(fwdTrackIDs[0]);
      const auto& trk2 = fwdTracks.iteratorAt(fwdTrackIDs[1]);

      if (trk1.eta() > fMinEtaMFT && trk1.eta() < fMaxEtaMFT) {
        firstInAcceptance = true;
      }
      if (trk2.eta() > fMinEtaMFT && trk2.eta() < fMaxEtaMFT) {
        secondInAcceptance = true;
      }

      // handle the four cases
      if (!firstInAcceptance && !secondInAcceptance) {
        // Case 1: Both outside MFT acceptance
        trkCandIDs.insert(trkCandIDs.end(), midTrackIDs->begin(), midTrackIDs->end());
      } else if (firstInAcceptance && !secondInAcceptance) {
        // Case 2: First inside, second outside
        trkCandIDs.push_back(fwdTrackIDs[0]);    // Keep first global
        trkCandIDs.push_back((*midTrackIDs)[1]); // Replace second with MCH-MID
      } else if (!firstInAcceptance && secondInAcceptance) {
        // Case 3: First outside, second inside
        trkCandIDs.push_back((*midTrackIDs)[0]); // Replace first with MCH-MID
        trkCandIDs.push_back(fwdTrackIDs[1]);    // Keep second global
      } else {
        // Case 4: Both inside MFT acceptance
        trkCandIDs.insert(trkCandIDs.end(), fwdTrackIDs.begin(), fwdTrackIDs.end());
      }

      uint64_t closestBcMCH = 0;
      upchelpers::FITInfo fitInfo{};
      fitInfo.timeFT0A = -999.f;
      fitInfo.timeFT0C = -999.f;
      fitInfo.timeFV0A = -999.f;
      fitInfo.ampFT0A = 0.f;
      fitInfo.ampFT0C = 0.f;
      fitInfo.ampFV0A = 0.f;
      std::vector<float> amplitudesT0A{};
      std::vector<float> amplitudesV0A{};
      std::vector<int8_t> relBCsT0A{};
      std::vector<int8_t> relBCsV0A{};
      uint8_t chFT0A = 0;
      uint8_t chFT0C = 0;
      if (nFT0s > 0) {
        uint64_t closestBcT0A = findClosestBC(globalBC, mapGlobalBcWithT0A);
        int64_t distClosestBcT0A = globalBC - static_cast<int64_t>(closestBcT0A);
        if (std::abs(distClosestBcT0A) <= fFilterFT0)
          continue;
        fitInfo.distClosestBcT0A = distClosestBcT0A;
        auto ft0Id = mapGlobalBcWithT0A.at(closestBcT0A);
        auto ft0 = ft0s.iteratorAt(ft0Id);
        fitInfo.timeFT0A = ft0.timeA();
        fitInfo.timeFT0C = ft0.timeC();
        const auto& t0AmpsA = ft0.amplitudeA();
        const auto& t0AmpsC = ft0.amplitudeC();
        fitInfo.ampFT0A = std::accumulate(t0AmpsA.begin(), t0AmpsA.end(), 0.f);
        fitInfo.ampFT0C = std::accumulate(t0AmpsC.begin(), t0AmpsC.end(), 0.f);
        chFT0A = ft0.amplitudeA().size();
        chFT0C = ft0.amplitudeC().size();
        fillAmplitudes(ft0s, mapGlobalBcWithT0A, amplitudesT0A, relBCsT0A, globalBC);
      }
      uint8_t chFV0A = 0;
      if (nFV0As > 0) {
        uint64_t closestBcV0A = findClosestBC(globalBC, mapGlobalBcWithV0A);
        int64_t distClosestBcV0A = globalBC - static_cast<int64_t>(closestBcV0A);
        if (std::abs(distClosestBcV0A) <= fFilterFV0)
          continue;
        fitInfo.distClosestBcV0A = distClosestBcV0A;
        auto fv0aId = mapGlobalBcWithV0A.at(closestBcV0A);
        auto fv0a = fv0as.iteratorAt(fv0aId);
        fitInfo.timeFV0A = fv0a.time();
        const auto& v0Amps = fv0a.amplitude();
        fitInfo.ampFV0A = std::accumulate(v0Amps.begin(), v0Amps.end(), 0.f);
        chFV0A = fv0a.amplitude().size();
        fillAmplitudes(fv0as, mapGlobalBcWithV0A, amplitudesV0A, relBCsV0A, globalBC);
      }
      uint8_t chFDDA = 0;
      uint8_t chFDDC = 0;
      if (nFDDs > 0) {
        uint64_t closestBcFDD = findClosestBC(globalBC, mapGlobalBcWithFDD);
        auto fddId = mapGlobalBcWithFDD.at(closestBcFDD);
        auto fdd = fdds.iteratorAt(fddId);
        fitInfo.timeFDDA = fdd.timeA();
        fitInfo.timeFDDC = fdd.timeC();
        fitInfo.ampFDDA = 0;
        for (int i = 0; i < 8; i++)
          fitInfo.ampFDDA += fdd.chargeA()[i];
        fitInfo.ampFDDC = 0;
        for (int i = 0; i < 8; i++)
          fitInfo.ampFDDC += fdd.chargeC()[i];
        fitInfo.triggerMaskFDD = fdd.triggerMask();
        // get signal coincidence
        for (int i = 0; i < 4; i++) {
          if (fdd.chargeA()[i + 4] > 0 && fdd.chargeA()[i] > 0)
            chFDDA++;
          if (fdd.chargeC()[i + 4] > 0 && fdd.chargeC()[i] > 0)
            chFDDC++;
        }
      }
      if (nZdcs > 0) {
        auto itZDC = mapGlobalBcWithZdc.find(globalBC);
        if (itZDC != mapGlobalBcWithZdc.end()) {
          const auto& zdc = zdcs.iteratorAt(itZDC->second);
          float timeZNA = zdc.timeZNA();
          float timeZNC = zdc.timeZNC();
          float eComZNA = zdc.energyCommonZNA();
          float eComZNC = zdc.energyCommonZNC();
          udZdcsReduced(candID, timeZNA, timeZNC, eComZNA, eComZNC);
        }
      }
      uint16_t numContrib = fNFwdProngs;
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      for (auto id : trkCandIDs) {
        auto tr = fwdTracks.iteratorAt(id);
        netCharge += tr.sign();
        selTrackIdsGlobal.push_back(id);
      }
      // store used tracks
      int upc_flag = 0;
      // TODO: introduce better check on association of collision and reconstruction mode
      if (bcs.iteratorAt(0).flags() == o2::itsmft::ROFRecord::VtxUPCMode)
        upc_flag = 1;
      fillFwdTracks(fwdTracks, trkCandIDs, candID, globalBC, closestBcMCH, mcFwdTrackLabels);
      eventCandidates(globalBC, runNumber, dummyX, dummyY, dummyZ, upc_flag, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      eventCandidatesSelExtras(chFT0A, chFT0C, chFDDA, chFDDC, chFV0A, 0, 0, 0, 0, 0);
      eventCandidatesSelsFwd(fitInfo.distClosestBcV0A,
                             fitInfo.distClosestBcT0A,
                             amplitudesT0A,
                             relBCsT0A,
                             amplitudesV0A,
                             relBCsV0A);
      candID++;
      trkCandIDs.clear();
    }

    selTrackIdsGlobal.clear();
    ambFwdTrBCs.clear();
    bcsMatchedTrIdsMID.clear();
    bcsMatchedTrIdsMCH.clear();
    bcsMatchedTrIdsGlobal.clear();
    mapGlobalBcWithT0A.clear();
    mapGlobalBcWithV0A.clear();
    mapGlobalBcWithFDD.clear();
  }

  // data processors
  // _________________________________________________

  // create candidates for forward/semiforward region
  // forward:     n fwd tracks + 0 barrel tracks
  // semiforward: n fwd tracks + m barrel tracks
  void processSemiFwd(ForwardTracks const& fwdTracks,
                      o2::aod::FwdTrkCls const& fwdTrkClusters,
                      o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                      BarrelTracks const& barrelTracks,
                      o2::aod::AmbiguousTracks const& ambTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    createCandidatesSemiFwd(barrelTracks, ambTracks,
                            fwdTracks, fwdTrkClusters, ambFwdTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as,
                            (o2::aod::McTrackLabels*)nullptr, (o2::aod::McFwdTrackLabels*)nullptr);
  }

  // create candidates for central region
  void processCentral(BarrelTracks const& barrelTracks,
                      o2::aod::AmbiguousTracks const& ambBarrelTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      o2::aod::Zdcs const& zdcs)
  {
    fDoMC = false;
    createCandidatesCentral(barrelTracks, ambBarrelTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as, zdcs,
                            (o2::aod::McTrackLabels*)nullptr);
  }

  // MC processors
  // _________________________________________________

  // create candidates for forward/semiforward region
  // forward:     n fwd tracks + 0 barrel tracks
  // semiforward: n fwd tracks + m barrel tracks
  void processSemiFwdMC(ForwardTracks const& fwdTracks,
                        o2::aod::FwdTrkCls const& fwdTrkClusters,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambTracks,
                        BCsWithBcSels const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McFwdTrackLabels const& mcFwdTrackLabels, o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesSemiFwd(barrelTracks, ambTracks,
                            fwdTracks, fwdTrkClusters, ambFwdTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as,
                            &mcBarrelTrackLabels, &mcFwdTrackLabels);
    fNewPartIDs.clear();
  }

  // create candidates for central region
  void processCentralMC(BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambBarrelTracks,
                        BCsWithBcSels const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        o2::aod::Zdcs const& zdcs,
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesCentral(barrelTracks, ambBarrelTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as, zdcs,
                            &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  // create candidates for forward region
  // forward: n fwd tracks
  void processForward(ForwardTracks const& fwdTracks,
                      o2::aod::FwdTrkCls const& fwdTrkClusters,
                      o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      o2::aod::Zdcs const& zdcs)
  {
    fDoMC = false;
    createCandidatesFwd(fwdTracks, fwdTrkClusters, ambFwdTracks,
                        bcs, collisions,
                        ft0s, fdds, fv0as, zdcs,
                        (o2::aod::McFwdTrackLabels*)nullptr);
  }

  void processForwardMC(ForwardTracks const& fwdTracks,
                        o2::aod::FwdTrkCls const& fwdTrkClusters,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        BCsWithBcSels const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        o2::aod::Zdcs const& zdcs,
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McFwdTrackLabels const& mcFwdTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesFwd(fwdTracks, fwdTrkClusters, ambFwdTracks,
                        bcs, collisions,
                        ft0s, fdds, fv0as, zdcs,
                        &mcFwdTrackLabels);
    fNewPartIDs.clear();
  }

  // create candidates for forward region including MFT
  // forward: n fwd tracks
  void processForwardGlobal(ForwardTracks const& fwdTracks,
                            o2::aod::FwdTrkCls const& fwdTrkClusters,
                            o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                            BCsWithBcSels const& bcs,
                            o2::aod::Collisions const& collisions,
                            o2::aod::FT0s const& ft0s,
                            o2::aod::FDDs const& fdds,
                            o2::aod::FV0As const& fv0as,
                            o2::aod::Zdcs const& zdcs)
  {
    fDoMC = false;
    createCandidatesFwdGlobal(fwdTracks, fwdTrkClusters, ambFwdTracks,
                              bcs, collisions,
                              ft0s, fdds, fv0as, zdcs,
                              (o2::aod::McFwdTrackLabels*)nullptr);
  }

  void processForwardGlobalMC(ForwardTracks const& fwdTracks,
                              o2::aod::FwdTrkCls const& fwdTrkClusters,
                              o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                              BCsWithBcSels const& bcs,
                              o2::aod::Collisions const& collisions,
                              o2::aod::FT0s const& ft0s,
                              o2::aod::FDDs const& fdds,
                              o2::aod::FV0As const& fv0as,
                              o2::aod::Zdcs const& zdcs,
                              o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                              o2::aod::McFwdTrackLabels const& mcFwdTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesFwdGlobal(fwdTracks, fwdTrkClusters, ambFwdTracks,
                              bcs, collisions,
                              ft0s, fdds, fv0as, zdcs,
                              &mcFwdTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward/forward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwdMC, "Produce candidates in semiforward/forward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processCentralMC, "Produce candidates in central region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processForward, "Produce candidates in forward region", false);
  PROCESS_SWITCH(UpcCandProducer, processForwardGlobal, "Produce candidates in forward region with MFT", true);
  PROCESS_SWITCH(UpcCandProducer, processForwardMC, "Produce candidates in forward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processForwardGlobalMC, "Produce candidates in forward region with MFT and MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
