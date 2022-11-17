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
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducerQa {
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisPt{500, 0., 5., ""};
    const AxisSpec axisEta{400, -2., 2., ""};
    const AxisSpec axisPhi{628, 0., 6.28, ""};
    const AxisSpec axisColNContrib{1001, -1., 1000., ""};

    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TRD", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_TRD_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TOF", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD", "", kTH2F, {axisPt, axisColNContrib});
    histRegistry.add("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD_TOF", "", kTH2F, {axisPt, axisColNContrib});

    histRegistry.add("TracksQA/Barrel/Eta/TPC", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TRD", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_TRD_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TOF", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TRD", "", kTH1F, {axisEta});
    histRegistry.add("TracksQA/Barrel/Eta/TPC_ITS_TRD_TOF", "", kTH1F, {axisEta});

    histRegistry.add("TracksQA/Barrel/Phi/TPC", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TRD", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_TRD_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TOF", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TRD", "", kTH1F, {axisPhi});
    histRegistry.add("TracksQA/Barrel/Phi/TPC_ITS_TRD_TOF", "", kTH1F, {axisPhi});
  }

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra>;

  void updateBarrelTrackQA(const BarrelTracks::iterator& track, int32_t colNContrib)
  {
    // check basic cuts
    // only tracks in TOF acceptance
    if (track.pt() < 0.3 || std::abs(track.eta()) > 0.8 || track.tpcNClsCrossedRows() < 70)
      return;

    int8_t mask = 0;
    if (track.hasTPC())
      SETBIT(mask, 0);
    if (track.hasITS())
      SETBIT(mask, 1);
    if (track.hasTOF())
      SETBIT(mask, 2);
    if (track.hasTRD())
      SETBIT(mask, 3);

    float pt = track.pt();
    float eta = track.eta();
    float phi = track.phi();

    if (mask == 1) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC"), phi);
    }
    if (mask == 3) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS"), phi);
    }
    if (mask == 5) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TOF"), phi);
    }
    if (mask == 7) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TOF"), phi);
    }
    if (mask == 9) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TRD"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TRD"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TRD"), phi);
    }
    if (mask == 11) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TRD"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TRD"), phi);
    }
    if (mask == 13) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_TRD_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_TRD_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_TRD_TOF"), phi);
    }
    if (mask == 15) {
      histRegistry.fill(HIST("TracksQA/Barrel/PtVsColNContrib/TPC_ITS_TRD_TOF"), pt, colNContrib);
      histRegistry.fill(HIST("TracksQA/Barrel/Eta/TPC_ITS_TRD_TOF"), eta);
      histRegistry.fill(HIST("TracksQA/Barrel/Phi/TPC_ITS_TRD_TOF"), phi);
    }
  }

  void process(o2::aod::Collisions const& collisions, BarrelTracks const& tracks, o2::aod::AmbiguousTracks const& ambTracks)
  {
    std::unordered_set<int64_t> ambTrIds;
    for (const auto& ambTrk : ambTracks) {
      auto trkId = ambTrk.trackId();
      ambTrIds.insert(trkId);
    }

    for (const auto& track : tracks) {
      int32_t nContrib = -1;
      if (ambTrIds.find(track.globalIndex()) == ambTrIds.end()) {
        const auto& col = track.collision();
        nContrib = col.numContrib();
      }
      updateBarrelTrackQA(track, nContrib);
    }

    ambTrIds.clear();
  }
};

struct UpcCandProducer {
  bool fDoMC{false};

  std::map<int32_t, int32_t> fNewPartIDs;
  uint64_t fMaxBC{0}; // max BC for ITS-TPC search

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;

  Produces<o2::aod::UDFwdTracks> udFwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> udFwdTracksExtra;
  Produces<o2::aod::UDMcFwdTrackLabels> udFwdTrackLabels;

  Produces<o2::aod::UDTracks> udTracks;
  Produces<o2::aod::UDTracksExtra> udTracksExtra;
  Produces<o2::aod::UDTracksDCA> udTracksDCA;
  Produces<o2::aod::UDTracksPID> udTracksPID;
  Produces<o2::aod::UDMcTrackLabels> udTrackLabels;
  Produces<o2::aod::UDTracksFlags> udTracksFlags;

  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSels> eventCandidatesSels;

  std::vector<bool> fwdSelectors;
  std::vector<bool> barrelSelectors;

  // skimmer flags
  // choose a source of signal MC events
  Configurable<int> fSignalGenID{"signalGenID", 1, "Signal generator ID"};

  // load cuts
  UPCCutparHolder upcCuts = UPCCutparHolder();
  MutableConfigurable<UPCCutparHolder> inputCuts{"UPCCuts", {}, "UPC event cuts"};

  // candidate producer flags
  Configurable<int> fFilterFT0{"filterFT0", 0, "Filter candidates by FT0 signals"};
  Configurable<int> fFilterRangeFT0{"filterRangeFT0", 0, "BC range (+/-) for filtration by FT0 signals"};
  Configurable<int> fSearchITSTPC{"searchITSTPC", 0, "Search for ITS-TPC tracks near candidates"};
  Configurable<int> fSearchRangeITSTPC{"searchRangeITSTPC", 50, "BC range for ITS-TPC tracks search wrt TOF tracks"};
  Configurable<uint32_t> fNFwdProngs{"nFwdProngs", 0, "Matched forward tracks per candidate"};
  Configurable<uint32_t> fNBarProngs{"nBarProngs", 2, "Matched barrel tracks per candidate"};

  // QA histograms
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

  typedef std::pair<uint64_t, std::vector<int64_t>> BCTracksPair;

  void init(InitContext&)
  {
    fwdSelectors.resize(upchelpers::kNFwdSels - 1, false);
    barrelSelectors.resize(upchelpers::kNBarrelSels - 1, false);

    upcCuts = (UPCCutparHolder)inputCuts;

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

  bool applyFwdCuts(const ForwardTracks::iterator& track)
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

  bool applyBarCuts(const BarrelTracks::iterator& track)
  {
    histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelAll, 1);

    // using any cuts at all?
    if (!upcCuts.getUseBarCuts()) {
      return true;
    }

    upchelpers::applyBarrelCuts(upcCuts, track, barrelSelectors);

    if (barrelSelectors[upchelpers::kBarrelSelHasTOF])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelHasTOF, 1);
    if (barrelSelectors[upchelpers::kBarrelSelPt])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelPt, 1);
    if (barrelSelectors[upchelpers::kBarrelSelEta])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelEta, 1);
    if (barrelSelectors[upchelpers::kBarrelSelITSNCls])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelITSNCls, 1);
    if (barrelSelectors[upchelpers::kBarrelSelITSChi2])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelITSChi2, 1);
    if (barrelSelectors[upchelpers::kBarrelSelTPCNCls])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelTPCNCls, 1);
    if (barrelSelectors[upchelpers::kBarrelSelTPCChi2])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelTPCChi2, 1);
    if (barrelSelectors[upchelpers::kBarrelSelDCAXY])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelDCAXY, 1);
    if (barrelSelectors[upchelpers::kBarrelSelDCAZ])
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelDCAZ, 1);

    bool pass = barrelSelectors[upchelpers::kBarrelSelPt] &&
                barrelSelectors[upchelpers::kBarrelSelEta] &&
                barrelSelectors[upchelpers::kBarrelSelITSNCls] &&
                barrelSelectors[upchelpers::kBarrelSelITSChi2] &&
                barrelSelectors[upchelpers::kBarrelSelTPCNCls] &&
                barrelSelectors[upchelpers::kBarrelSelTPCChi2] &&
                barrelSelectors[upchelpers::kBarrelSelDCAXY] &&
                barrelSelectors[upchelpers::kBarrelSelDCAZ];
    return pass;
  }

  void skimMCInfo(o2::aod::McCollisions const& mcCollisions,
                  o2::aod::McParticles const& mcParticles,
                  BCsWithBcSels const& bcs)
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

    // storing MC particles
    for (const auto& item : fNewPartIDs) {
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
        auto it = fNewPartIDs.find(motherID);
        if (it != fNewPartIDs.end()) {
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
        auto itFirst = fNewPartIDs.find(firstDaughter);
        auto itLast = fNewPartIDs.find(lastDaughter);
        if (itFirst != fNewPartIDs.end() && itLast != fNewPartIDs.end()) {
          newDaughterIDs[0] = fNewPartIDs.at(daughterIDs.front());
          newDaughterIDs[1] = fNewPartIDs.at(daughterIDs.back());
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
      auto bc = mcEvent.bc_as<BCsWithBcSels>();
      udMCCollisions(bc.globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                     mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  template <typename TTrack, typename TAmbTracks>
  uint64_t getTrackBC(TTrack track,
                      TAmbTracks const& ambTracks,
                      int64_t ambTrId,
                      o2::aod::Collisions const& collisions,
                      BCsWithBcSels const& bcs)
  {
    uint64_t trackBC = 0;
    if (ambTrId >= 0) {
      const auto& ambTr = ambTracks.iteratorAt(ambTrId);
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      trackBC = first.globalBC();
    } else {
      const auto& col = track.collision();
      trackBC = col.template bc_as<BCsWithBcSels>().globalBC();
    }
    int64_t tint = std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
    uint64_t bc = trackBC + tint;
    return bc;
  }

  void fillFwdTracks(ForwardTracks const& tracks,
                     std::vector<int64_t> const& trackIDs,
                     int32_t candID,
                     uint64_t bc,
                     const o2::aod::McFwdTrackLabels* mcTrackLabels)
  {
    for (auto trackID : trackIDs) {
      const auto& track = tracks.iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      udFwdTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udFwdTracksExtra(track.nClusters(), track.pDca(), track.rAtAbsorberEnd(), track.chi2(), track.chi2MatchMCHMID(),
                       track.mchBitMap(), track.midBitMap(), track.midBoards());
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

  void fillBarrelTracks(BarrelTracks const& tracks,
                        std::vector<int64_t> const& trackIDs,
                        int32_t candID,
                        uint64_t bc,
                        const o2::aod::McTrackLabels* mcTrackLabels,
                        std::unordered_map<int64_t, int64_t>& ambBarrelTrIds)
  {
    for (auto trackID : trackIDs) {
      const auto& track = tracks.iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      int64_t colId = -1;
      if (ambBarrelTrIds.find(trackID) == ambBarrelTrIds.end()) {
        colId = track.collisionId();
      }
      udTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udTracksExtra(track.itsClusterMap(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.trdPattern(), track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                    track.tpcSignal(), track.tofSignal(), track.trdSignal(), track.length(), track.tofExpMom(), track.detectorMap());
      udTracksPID(track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
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

  void processFITInfo(upchelpers::FITInfo& fitInfo,
                      uint64_t midbc,
                      std::vector<std::pair<uint64_t, int64_t>>& v,
                      BCsWithBcSels const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
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
      if (bc.has_foundFT0()) {
        const auto& ft0 = bc.foundFT0();
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
        if (!bc.selection()[evsel::kNoBGT0A])
          SETBIT(fitInfo.BGFT0Apf, bit);
        if (!bc.selection()[evsel::kNoBGT0C])
          SETBIT(fitInfo.BGFT0Cpf, bit);
        if (bc.selection()[evsel::kIsBBT0A])
          SETBIT(fitInfo.BBFT0Apf, bit);
        if (bc.selection()[evsel::kIsBBT0C])
          SETBIT(fitInfo.BBFT0Cpf, bit);
      }
      if (bc.has_foundFV0()) {
        const auto& fv0a = bc.foundFV0();
        fitInfo.timeFV0A = fv0a.time();
        const auto& amps = fv0a.amplitude();
        fitInfo.ampFV0A = 0.;
        for (auto amp : amps)
          fitInfo.ampFV0A += amp;
        fitInfo.triggerMaskFV0A = fv0a.triggerMask();
        if (!bc.selection()[evsel::kNoBGV0A])
          SETBIT(fitInfo.BGFV0Apf, bit);
        if (bc.selection()[evsel::kIsBBV0A])
          SETBIT(fitInfo.BBFV0Apf, bit);
      }
      if (bc.has_foundFDD()) {
        const auto& fdd = bc.foundFDD();
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
        if (!bc.selection()[evsel::kNoBGFDA])
          SETBIT(fitInfo.BGFDDApf, bit);
        if (!bc.selection()[evsel::kNoBGFDC])
          SETBIT(fitInfo.BGFDDCpf, bit);
        if (bc.selection()[evsel::kIsBBFDA])
          SETBIT(fitInfo.BBFDDApf, bit);
        if (bc.selection()[evsel::kIsBBFDC])
          SETBIT(fitInfo.BBFDDCpf, bit);
      }
      ++curit;
      if (curit == v.end())
        break;
      curbc = curit->first;
    }
  }

  template <int32_t tracksSwitch, typename TAmbTrack>
  int64_t getAmbTrackId(TAmbTrack& ambTrack)
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

  template <int32_t tracksSwitch, typename TAmbTracks>
  void collectAmbTracks(std::unordered_map<int64_t, int64_t>& ambTrIds,
                        TAmbTracks const& ambTracks)
  {
    for (const auto& ambTrk : ambTracks) {
      auto trkId = getAmbTrackId<tracksSwitch>(ambTrk);
      ambTrIds[trkId] = ambTrk.globalIndex();
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

  void collectBarrelTracks(std::vector<BCTracksPair>& bcsMatchedTrIdsTOF,
                           std::vector<BCTracksPair>& bcsMatchedTrIdsITSTPC,
                           BCsWithBcSels const& bcs,
                           o2::aod::Collisions const& collisions,
                           BarrelTracks const& barrelTracks,
                           o2::aod::AmbiguousTracks const& ambBarrelTracks,
                           std::unordered_map<int64_t, int64_t>& ambBarrelTrIds)
  {
    for (const auto& trk : barrelTracks) {
      if (!applyBarCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int64_t ambTrId = -1;
      int32_t nContrib = -1;
      auto ambIter = ambBarrelTrIds.find(trkId);
      if (ambIter == ambBarrelTrIds.end()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
      } else {
        ambTrId = ambIter->second;
      }
      uint64_t bc = getTrackBC(trk, ambBarrelTracks, ambTrId, collisions, bcs);
      if (bc > fMaxBC)
        continue;
      if (!upcCuts.getRequireITSTPC() && trk.hasTOF() && nContrib <= upcCuts.getMaxNContrib())
        addTrack(bcsMatchedTrIdsTOF, bc, trkId);
      if (upcCuts.getRequireITSTPC() && trk.hasTOF() && trk.hasITS() && trk.hasTPC() && nContrib <= upcCuts.getMaxNContrib())
        addTrack(bcsMatchedTrIdsTOF, bc, trkId);
      if (fSearchITSTPC == 1 && !trk.hasTOF() && trk.hasITS() && trk.hasTPC())
        addTrack(bcsMatchedTrIdsITSTPC, bc, trkId);
    }
  }

  void collectForwardTracks(std::vector<BCTracksPair>& bcsMatchedTrIdsMID,
                            BCsWithBcSels const& bcs,
                            o2::aod::Collisions const& collisions,
                            ForwardTracks const& fwdTracks,
                            o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                            std::unordered_map<int64_t, int64_t>& ambFwdTrIds)
  {
    for (const auto& trk : fwdTracks) {
      if (!applyFwdCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int64_t ambTrId = -1;
      int32_t nContrib = -1;
      auto ambIter = ambFwdTrIds.find(trkId);
      if (ambIter == ambFwdTrIds.end()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
      } else {
        ambTrId = ambIter->second;
      }
      uint64_t bc = getTrackBC(trk, ambFwdTracks, ambTrId, collisions, bcs);
      if (bc > fMaxBC)
        continue;
      auto trkType = trk.trackType();
      if (trkType == o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack && nContrib <= upcCuts.getMaxNContrib())
        addTrack(bcsMatchedTrIdsMID, bc, trkId);
    }
  }

  int32_t searchTracks(uint64_t midbc, uint64_t range, uint32_t tracksToFind,
                       std::vector<int64_t>& tracks,
                       std::vector<BCTracksPair>& v)
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
      uint32_t size = curit->second.size();
      if (size > 1) // too many tracks per BC -> possibly another event
        return -2;
      count += size;
      if (count > tracksToFind) // too many tracks nearby
        return -3;
      tracks.push_back(curit->second[0]);
      ++curit;
      if (curit == v.end())
        break;
      curbc = curit->first;
    }
    if (count != tracksToFind)
      return -4;
    return 0; // found exactly tracksToFind tracks in [midbc - range, midbc + range]
  }

  void createCandidatesCentral(BarrelTracks const& barrelTracks,
                               o2::aod::AmbiguousTracks const& ambBarrelTracks,
                               BCsWithBcSels const& bcs,
                               o2::aod::Collisions const& collisions,
                               o2::aod::FT0s const& ft0s,
                               o2::aod::FDDs const& fdds,
                               o2::aod::FV0As const& fv0as,
                               const o2::aod::McTrackLabels* mcBarrelTrackLabels)
  {
    fMaxBC = bcs.iteratorAt(bcs.size() - 1).globalBC(); // restrict ITS-TPC track search to [0, fMaxBC]

    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsTOF;
    std::vector<BCTracksPair> bcsMatchedTrIdsITSTPC;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, int64_t> ambBarrelTrIds;
    collectAmbTracks<0>(ambBarrelTrIds, ambBarrelTracks);

    collectBarrelTracks(bcsMatchedTrIdsTOF, bcsMatchedTrIdsITSTPC,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrIds);

    uint32_t nBCsWithITSTPC = bcsMatchedTrIdsITSTPC.size();

    std::sort(bcsMatchedTrIdsTOF.begin(), bcsMatchedTrIdsTOF.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::sort(bcsMatchedTrIdsITSTPC.begin(), bcsMatchedTrIdsITSTPC.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    if (nBCsWithITSTPC > 0 && fSearchITSTPC == 1) {
      for (auto& pair : bcsMatchedTrIdsTOF) {
        uint64_t bc = pair.first;
        auto& trackIds = pair.second;
        uint32_t nTOFtracks = trackIds.size();
        if (nTOFtracks > fNBarProngs) // too many TOF tracks?!
          continue;
        if (nTOFtracks == fNBarProngs) { // check for ITS-TPC tracks
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, 0, tracks, bcsMatchedTrIdsITSTPC);
          if (res < 0) { // too many tracks nearby -> rejecting
            trackIds.push_back(0);
            continue;
          }
        }
        if (nTOFtracks < fNBarProngs && !upcCuts.getRequireTOF()) { // add ITS-TPC track if needed
          uint32_t tracksToFind = fNBarProngs - nTOFtracks;
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, tracksToFind, tracks, bcsMatchedTrIdsITSTPC);
          if (res < 0) // too many or not enough tracks nearby -> rejecting
            continue;
          trackIds.insert(trackIds.end(), tracks.begin(), tracks.end());
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
    for (const auto& item : bcsMatchedTrIdsTOF) {
      auto& barrelTrackIDs = item.second;
      uint16_t numContrib = barrelTrackIDs.size();
      // sanity check
      if (numContrib != fNBarProngs)
        continue;
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      uint64_t bc = item.first;
      upchelpers::FITInfo fitInfo{};
      processFITInfo(fitInfo, bc, indexBCglId, bcs, ft0s, fdds, fv0as);
      if (fFilterFT0) {
        if (!checkFT0(fitInfo, true))
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
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels, ambBarrelTrIds);
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      candID++;
    }

    indexBCglId.clear();
    ambBarrelTrIds.clear();
    bcsMatchedTrIdsTOF.clear();
  }

  void createCandidatesSemiFwd(BarrelTracks const& barrelTracks,
                               o2::aod::AmbiguousTracks const& ambBarrelTracks,
                               ForwardTracks const& fwdTracks,
                               o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                               BCsWithBcSels const& bcs,
                               o2::aod::Collisions const& collisions,
                               o2::aod::FT0s const& ft0s,
                               o2::aod::FDDs const& fdds,
                               o2::aod::FV0As const& fv0as,
                               const o2::aod::McTrackLabels* mcBarrelTrackLabels,
                               const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    fMaxBC = bcs.iteratorAt(bcs.size() - 1).globalBC(); // restrict ITS-TPC track search to [0, fMaxBC]

    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsTOF;
    std::vector<BCTracksPair> bcsMatchedTrIdsITSTPC;
    std::vector<BCTracksPair> bcsMatchedTrIdsMID;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, int64_t> ambBarrelTrIds;
    collectAmbTracks<0>(ambBarrelTrIds, ambBarrelTracks);

    std::unordered_map<int64_t, int64_t> ambFwdTrIds;
    collectAmbTracks<1>(ambFwdTrIds, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrIds);

    collectBarrelTracks(bcsMatchedTrIdsTOF, bcsMatchedTrIdsITSTPC,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrIds);

    uint32_t nBCsWithITSTPC = bcsMatchedTrIdsITSTPC.size();
    uint32_t nBCsWithMID = bcsMatchedTrIdsMID.size();

    // tag TOF tracks to MID tracks
    std::unordered_set<uint64_t> bcsWithMID;
    for (const auto& pair : bcsMatchedTrIdsMID) {
      bcsWithMID.insert(pair.first);
    }

    std::vector<BCTracksPair> bcsMatchedTrIdsTOFTagged;
    bcsMatchedTrIdsTOFTagged.reserve(nBCsWithMID);
    for (const auto& pair : bcsMatchedTrIdsTOF) {
      if (bcsWithMID.find(pair.first) != bcsWithMID.end())
        bcsMatchedTrIdsTOFTagged.emplace_back(pair);
    }

    bcsWithMID.clear();
    bcsMatchedTrIdsTOF.clear();

    std::sort(bcsMatchedTrIdsTOFTagged.begin(), bcsMatchedTrIdsTOFTagged.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::sort(bcsMatchedTrIdsITSTPC.begin(), bcsMatchedTrIdsITSTPC.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::sort(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    if (nBCsWithITSTPC > 0 && fSearchITSTPC == 1) {
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
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, 0, tracks, bcsMatchedTrIdsITSTPC);
          if (res < 0) { // too many tracks nearby -> rejecting
            trackIdsTOF.push_back(0);
            continue;
          }
        }
        if (nMIDtracks == fNFwdProngs && nTOFtracks < fNBarProngs && !upcCuts.getRequireTOF()) { // add ITS-TPC track if needed
          uint32_t tracksToFind = fNBarProngs - nTOFtracks;
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, tracksToFind, tracks, bcsMatchedTrIdsITSTPC);
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
      // sanity check
      if (nBarrelTracks != fNBarProngs || nMIDtracks != fNFwdProngs)
        continue;
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      uint64_t bc = pairMID.first;
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
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, mcFwdTrackLabels);
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels, ambBarrelTrIds);
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      candID++;
    }

    indexBCglId.clear();
    ambFwdTrIds.clear();
    bcsMatchedTrIdsMID.clear();
    ambBarrelTrIds.clear();
    bcsMatchedTrIdsTOFTagged.clear();
  }

  // data processors
  // _________________________________________________

  // create candidates for forward/semiforward region
  // forward:     n fwd tracks + 0 barrel tracks
  // semiforward: n fwd tracks + m barrel tracks
  void processSemiFwd(ForwardTracks const& fwdTracks,
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
                            fwdTracks, ambFwdTracks,
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
                      o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    createCandidatesCentral(barrelTracks, ambBarrelTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as,
                            (o2::aod::McTrackLabels*)nullptr);
  }

  // MC processors
  // _________________________________________________

  // create candidates for forward/semiforward region
  // forward:     n fwd tracks + 0 barrel tracks
  // semiforward: n fwd tracks + m barrel tracks
  void processSemiFwdMC(ForwardTracks const& fwdTracks,
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
                            fwdTracks, ambFwdTracks,
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
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesCentral(barrelTracks, ambBarrelTracks,
                            bcs, collisions,
                            ft0s, fdds, fv0as,
                            &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward/forward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwdMC, "Produce candidates in semiforward/forward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processCentralMC, "Produce candidates in central region with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducerQa>(cfgc),
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
