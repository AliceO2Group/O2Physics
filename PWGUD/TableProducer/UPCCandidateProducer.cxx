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

struct UpcCandProducer {
  bool fDoMC{false};
  bool fDoSemiFwd{false};

  std::map<int32_t, int32_t> fNewPartIDs;

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
  Configurable<int> fCheckTPCPID{"checkTPCPID", 0, "Check TPC PID. Useful for central selection -- see `tpcPIDSwitch` option"};
  Configurable<int> fTPCPIDSwitch{"tpcPIDSwitch", 0, "PID switch: 0 -- two muons/pions, 1 -- two electrons, 2 -- electron + muon/pion"};
  Configurable<int> fFilterFT0{"filterFT0", 0, "Filter candidates by FT0 signals"};
  Configurable<int> fFilterRangeFT0{"filterRangeFT0", 0, "BC range (+/-) for filtration by FT0 signals"};
  Configurable<int> fNFwdProngs{"nFwdProngs", 2, "Matched forward tracks per candidate"};
  Configurable<int> fNBarProngs{"nBarProngs", 0, "Matched barrel tracks per candidate"};

  // QA histograms
  Configurable<int> fCollectBTracksQA{"collectBTracksQA", 0, "Collect kinematic distributions for all filtered central barrel tracks"};
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

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

    if (fCollectBTracksQA) {
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

    bool pass = barrelSelectors[upchelpers::kBarrelSelHasTOF] &&
                barrelSelectors[upchelpers::kBarrelSelPt] &&
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

  // trackSwitch: 0 -> forward, 1 -> barrel
  template <int32_t trackSwitch, typename TTracks>
  void filterTracks(TTracks* tracks,
                    o2::aod::Collisions const& collisions,
                    std::unordered_map<int32_t, int32_t>& ambTrIDs,
                    std::vector<int32_t>& filteredTrackIDs)
  {
    for (const auto& tr : *tracks) {
      // skip if track doesn't pass selection
      bool pass = false;
      if constexpr (trackSwitch == 0) {
        if (upcCuts.getTrackType() != -1) {
          auto trType = tr.trackType();
          if (trType != upcCuts.getTrackType()) {
            continue;
          }
        }
        pass = applyFwdCuts(tr);
      }
      if constexpr (trackSwitch == 1) {
        pass = applyBarCuts(tr);
      }
      int32_t colNContrib = -1;
      auto ambIter = ambTrIDs.find(tr.globalIndex());
      if (ambIter == ambTrIDs.end()) {
        int32_t colId = tr.collisionId();
        const auto& col = collisions.iteratorAt(colId);
        colNContrib = col.numContrib();
      }
      if (pass) {
        if ((colNContrib >= 0 && colNContrib <= upcCuts.getMaxNContrib()) || (colNContrib == -1 && upcCuts.getAmbigSwitch() != 1))
          filteredTrackIDs.push_back(tr.globalIndex());
      }
      if constexpr (trackSwitch == 1) {
        if (fCollectBTracksQA)
          updateBTrackQA(tr, colNContrib);
      }
    }
  }

  void updateBTrackQA(const BarrelTracks::iterator& track, int32_t colNContrib)
  {
    // check basic cuts
    // only tracks in TOF acceptance
    if (!barrelSelectors[upchelpers::kBarrelSelPt] ||
        std::abs(track.eta()) > 0.8 ||
        !barrelSelectors[upchelpers::kBarrelSelTPCNCls])
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

  template <typename TTrack, typename TAmbTracks>
  uint64_t getTrackBC(TTrack track,
                      TAmbTracks* ambTracks,
                      std::unordered_map<int32_t, int32_t>& ambTrIDs,
                      o2::aod::Collisions const& collisions,
                      BCsWithBcSels const& bcs)
  {
    uint64_t trackBC = 0;
    auto ambIter = ambTrIDs.find(track.globalIndex());
    if (ambIter == ambTrIDs.end()) {
      const auto& col = track.collision();
      trackBC = col.template bc_as<BCsWithBcSels>().globalBC();
    } else if (upcCuts.getAmbigSwitch() != 1) {
      const auto& ambTr = ambTracks->iteratorAt(ambIter->second);
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      trackBC = first.globalBC();
    }
    int64_t tint = std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
    uint64_t bc = trackBC + tint;
    return bc;
  }

  template <typename TFwdTracks, typename TMcFwdTrackLabels>
  void fillFwdTracks(TFwdTracks* tracks,
                     std::vector<int32_t> const& trackIDs,
                     int32_t candID,
                     uint64_t bc,
                     TMcFwdTrackLabels* mcTrackLabels)
  {
    for (int32_t trackID : trackIDs) {
      const auto& track = tracks->iteratorAt(trackID);
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

  template <typename TBarrelTracks, typename TMcTrackLabels>
  void fillBarrelTracks(TBarrelTracks* tracks,
                        std::vector<int32_t> const& trackIDs,
                        int32_t candID,
                        uint64_t bc,
                        TMcTrackLabels* mcTrackLabels,
                        std::unordered_map<int32_t, int32_t>& ambBarrelTrIds)
  {
    for (int32_t trackID : trackIDs) {
      const auto& track = tracks->iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      bool isAmbiguous = false;
      if (upcCuts.getAmbigSwitch() != 1) {
        isAmbiguous = ambBarrelTrIds.find(trackID) != ambBarrelTrIds.end();
      }
      udTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udTracksExtra(track.itsClusterMap(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.trdPattern(), track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                    track.tpcSignal(), track.tofSignal(), track.trdSignal(), track.length(), track.tofExpMom(), track.detectorMap());
      udTracksPID(track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                  track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
      udTracksDCA(track.dcaZ(), track.dcaXY());
      udTracksFlags(isAmbiguous, track.isPVContributor());
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

  // filter tracks in the central barrel using TPC PID
  // return `true` if a candidate passes "number of tracks" requirement from `nBarProngs`, `false` otherwise
  // if candidate passes, `filteredTrackIDs` contain track IDs of passed tracks
  template <typename TBarrelTracks>
  bool checkTPCPID(TBarrelTracks* tracks, std::vector<int32_t> const& trackIDs, std::vector<int32_t>& filteredTrackIDs)
  {
    int32_t nIDs = trackIDs.size();
    std::vector<bool> pidFlagsEl(nIDs, false);
    std::vector<bool> pidFlagsPi(nIDs, false);
    int32_t nEl = 0;
    int32_t nPi = 0;

    for (int32_t itr = 0; itr < nIDs; itr++) {
      const auto& tr = tracks->iteratorAt(trackIDs[itr]);
      float tpcNSigmaEl = tr.tpcNSigmaEl();
      float tpcNSigmaPi = tr.tpcNSigmaPi();

      bool isEl = std::abs(tpcNSigmaEl) < 3.0f && std::abs(tpcNSigmaEl) < std::abs(tpcNSigmaPi);
      bool isPi = std::abs(tpcNSigmaPi) < 3.0f && std::abs(tpcNSigmaEl) > std::abs(tpcNSigmaPi);

      pidFlagsEl[itr] = isEl;
      if (isEl) {
        nEl++;
      }

      pidFlagsPi[itr] = isPi;
      if (isPi) {
        nPi++;
      }
    }

    bool pass = false;

    // two muons/pions
    if (fTPCPIDSwitch == 0 && nPi == 2 && nEl == 0) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // two electrons
    if (fTPCPIDSwitch == 1 && nEl == 2) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // electron + muon/pion
    if (fTPCPIDSwitch == 2 && nEl == 1 && nPi == 1) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    return pass;
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
  void processFITInfo(TBCs const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>>& bcsMatchedTrIds,
                      std::unordered_map<uint64_t, upchelpers::FITInfo>& bcsWithFIT)
  {
    const uint64_t presBitNum = 16; // bit for present BC

    std::unordered_set<uint64_t> addedBCs;
    std::vector<std::pair<uint64_t, upchelpers::FITInfo>> vbcsWithFIT;

    // gather FIT information for each BC
    // even if there is no FIT info at all -> store dummy "info{}"
    for (const auto& bc : bcs) {
      int32_t ft0Id = bc.foundFT0Id();
      upchelpers::FITInfo info{};
      if (ft0Id != -1) {
        const auto& ft0 = ft0s.iteratorAt(ft0Id);
        info.timeFT0A = ft0.timeA();
        info.timeFT0C = ft0.timeC();
        const auto& ampsA = ft0.amplitudeA();
        const auto& ampsC = ft0.amplitudeC();
        info.ampFT0A = 0.;
        for (auto amp : ampsA) {
          info.ampFT0A += amp;
        }
        info.ampFT0C = 0.;
        for (auto amp : ampsC) {
          info.ampFT0C += amp;
        }
        info.triggerMaskFT0 = ft0.triggerMask();
        if (!bc.selection()[evsel::kNoBGT0A])
          SETBIT(info.BGFT0Apf, presBitNum);
        if (!bc.selection()[evsel::kNoBGT0C])
          SETBIT(info.BGFT0Cpf, presBitNum);
        if (bc.selection()[evsel::kIsBBT0A])
          SETBIT(info.BBFT0Apf, presBitNum);
        if (bc.selection()[evsel::kIsBBT0C])
          SETBIT(info.BBFT0Cpf, presBitNum);
      }
      int32_t fddId = bc.foundFDDId();
      if (fddId != -1) {
        const auto& fdd = fdds.iteratorAt(fddId);
        info.timeFDDA = fdd.timeA();
        info.timeFDDC = fdd.timeC();
        const auto& ampsA = fdd.chargeA();
        const auto& ampsC = fdd.chargeC();
        info.ampFDDA = 0.;
        for (auto amp : ampsA) {
          info.ampFDDA += amp;
        }
        info.ampFDDC = 0.;
        for (auto amp : ampsC) {
          info.ampFDDC += amp;
        }
        info.triggerMaskFDD = fdd.triggerMask();
        if (!bc.selection()[evsel::kNoBGFDA])
          SETBIT(info.BGFDDApf, presBitNum);
        if (!bc.selection()[evsel::kNoBGFDC])
          SETBIT(info.BGFDDCpf, presBitNum);
        if (bc.selection()[evsel::kIsBBFDA])
          SETBIT(info.BBFDDApf, presBitNum);
        if (bc.selection()[evsel::kIsBBFDC])
          SETBIT(info.BBFDDCpf, presBitNum);
      }
      int32_t fv0aId = bc.foundFV0Id();
      if (fv0aId != -1) {
        const auto& fv0a = fv0as.iteratorAt(fv0aId);
        info.timeFV0A = fv0a.time();
        const auto& amps = fv0a.amplitude();
        info.ampFV0A = 0.;
        for (auto amp : amps) {
          info.ampFV0A += amp;
        }
        info.triggerMaskFV0A = fv0a.triggerMask();
        if (!bc.selection()[evsel::kNoBGV0A])
          SETBIT(info.BGFV0Apf, presBitNum);
        if (bc.selection()[evsel::kIsBBV0A])
          SETBIT(info.BBFV0Apf, presBitNum);
      }
      vbcsWithFIT.push_back({bc.globalBC(), info});
      addedBCs.insert(bc.globalBC());
    }

    // add BCs from event candidates if not yet there
    for (const auto& item : bcsMatchedTrIds) {
      uint64_t bc = item.first;
      if (addedBCs.find(bc) == addedBCs.end()) {
        addedBCs.insert(bc);
        upchelpers::FITInfo info{};
        vbcsWithFIT.push_back({bc, info});
      }
    }

    addedBCs.clear();

    std::sort(vbcsWithFIT.begin(), vbcsWithFIT.end(), [](const auto& left, const auto& right) { return left.first < right.first; });

    int32_t nBCsWithFIT = vbcsWithFIT.size();
    for (int32_t ibc = 0; ibc < nBCsWithFIT; ++ibc) {
      uint64_t bc = vbcsWithFIT[ibc].first;
      auto& info = vbcsWithFIT[ibc].second;
      // scrolling back on the original vector
      for (int32_t jbc = ibc - 1; jbc >= 0; --jbc) {
        uint64_t pastbc = vbcsWithFIT[jbc].first;
        uint64_t deltaBC = bc - pastbc;
        if (deltaBC > 16) // range [1, 16]
          break;
        const auto& pastInfo = vbcsWithFIT[jbc].second;
        uint64_t pastBitNum = deltaBC - 1; // bits range in [0, 15]; bit == 16 -> present
        if (TESTBIT(pastInfo.BGFT0Apf, presBitNum))
          SETBIT(info.BGFT0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BGFT0Cpf, presBitNum))
          SETBIT(info.BGFT0Cpf, pastBitNum);
        if (TESTBIT(pastInfo.BBFT0Apf, presBitNum))
          SETBIT(info.BBFT0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BBFT0Cpf, presBitNum))
          SETBIT(info.BBFT0Cpf, pastBitNum);
        if (TESTBIT(pastInfo.BGFDDApf, presBitNum))
          SETBIT(info.BGFDDApf, pastBitNum);
        if (TESTBIT(pastInfo.BGFDDCpf, presBitNum))
          SETBIT(info.BGFDDCpf, pastBitNum);
        if (TESTBIT(pastInfo.BBFDDApf, presBitNum))
          SETBIT(info.BBFDDApf, pastBitNum);
        if (TESTBIT(pastInfo.BBFDDCpf, presBitNum))
          SETBIT(info.BBFDDCpf, pastBitNum);
        if (TESTBIT(pastInfo.BGFV0Apf, presBitNum))
          SETBIT(info.BGFV0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BBFV0Apf, presBitNum))
          SETBIT(info.BBFV0Apf, pastBitNum);
      }
      // scrolling forward on the original vector
      for (int32_t jbc = ibc + 1; jbc <= nBCsWithFIT; ++jbc) {
        uint64_t futbc = vbcsWithFIT[jbc].first;
        uint64_t deltaBC = futbc - bc;
        if (deltaBC > 15) // [1, 15]
          break;
        const auto& futInfo = vbcsWithFIT[jbc].second;
        uint64_t futBitNum = deltaBC + 16; // bits range in [17, 31]
        if (TESTBIT(futInfo.BGFT0Apf, presBitNum))
          SETBIT(info.BGFT0Apf, futBitNum);
        if (TESTBIT(futInfo.BGFT0Cpf, presBitNum))
          SETBIT(info.BGFT0Cpf, futBitNum);
        if (TESTBIT(futInfo.BBFT0Apf, presBitNum))
          SETBIT(info.BBFT0Apf, futBitNum);
        if (TESTBIT(futInfo.BBFT0Cpf, presBitNum))
          SETBIT(info.BBFT0Cpf, futBitNum);
        if (TESTBIT(futInfo.BGFDDApf, presBitNum))
          SETBIT(info.BGFDDApf, futBitNum);
        if (TESTBIT(futInfo.BGFDDCpf, presBitNum))
          SETBIT(info.BGFDDCpf, futBitNum);
        if (TESTBIT(futInfo.BBFDDApf, presBitNum))
          SETBIT(info.BBFDDApf, futBitNum);
        if (TESTBIT(futInfo.BBFDDCpf, presBitNum))
          SETBIT(info.BBFDDCpf, futBitNum);
        if (TESTBIT(futInfo.BGFV0Apf, presBitNum))
          SETBIT(info.BGFV0Apf, futBitNum);
        if (TESTBIT(futInfo.BBFV0Apf, presBitNum))
          SETBIT(info.BBFV0Apf, futBitNum);
      }
      // add current pair to map
      bcsWithFIT.insert({bc, info});
    }

    vbcsWithFIT.clear();
  }

  template <typename TFwdTracks, typename TAmbiguousFwdTracks,
            typename TBarrelTracks, typename TAmbiguousTracks,
            typename TBCs,
            typename TMcFwdTrackLabels, typename TMcTrackLabels>
  void createCandidates(TFwdTracks* fwdTracks,
                        TBarrelTracks* barrelTracks,
                        TAmbiguousFwdTracks* ambFwdTracks,
                        TAmbiguousTracks* ambBarrelTracks,
                        TBCs const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        TMcFwdTrackLabels* mcFwdTrackLabels,
                        TMcTrackLabels* mcBarrelTrackLabels)
  {
    // pairs of global BCs and vectors of matched track IDs:
    // global BC <-> <vector of fwd. trackIDs, vector of barrel trackIDs>
    std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> bcsMatchedTrIds;

    // track IDs after cuts
    std::vector<int32_t> filteredFwdTracks;
    std::vector<int32_t> filteredBarrelTracks;

    // trackID -> index in amb. track table
    std::unordered_map<int32_t, int32_t> ambFwdTrIds;
    std::unordered_map<int32_t, int32_t> ambBarrelTrIds;

    // forward matching
    if (fwdTracks != nullptr) {
      if (upcCuts.getAmbigSwitch() != 1) {
        for (const auto& ambTr : *ambFwdTracks) {
          auto trId = ambTr.fwdtrackId();
          ambFwdTrIds[trId] = ambTr.globalIndex();
        }
      }
      filterTracks<0>(fwdTracks, collisions, ambFwdTrIds, filteredFwdTracks);
      for (int32_t fwdTrID : filteredFwdTracks) {
        const auto& fwdTr = fwdTracks->iteratorAt(fwdTrID);
        uint64_t bc = getTrackBC(fwdTr, ambFwdTracks, ambFwdTrIds, collisions, bcs);
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.first.emplace_back(fwdTr.globalIndex());
        } else {
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(1, fwdTrID), std::vector<int32_t>());
        }
      }
    }

    // central barrel tracks
    if (barrelTracks != nullptr) {
      if (upcCuts.getAmbigSwitch() != 1 || fCollectBTracksQA) {
        for (const auto& ambTr : *ambBarrelTracks) {
          auto trId = ambTr.trackId();
          ambBarrelTrIds[trId] = ambTr.globalIndex();
        }
      }
      filterTracks<1>(barrelTracks, collisions, ambBarrelTrIds, filteredBarrelTracks);
      for (int32_t barTrID : filteredBarrelTracks) {
        const auto& barTr = barrelTracks->iteratorAt(barTrID);
        uint64_t bc = getTrackBC(barTr, ambBarrelTracks, ambBarrelTrIds, collisions, bcs);
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.second.emplace_back(barTr.globalIndex());
        } else if (!fDoSemiFwd) { // tag central-barrel tracks to forward tracks in semiforward case
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(), std::vector<int32_t>(1, barTrID));
        }
      }
    }

    std::unordered_map<uint64_t, upchelpers::FITInfo> bcsWithFIT;
    processFITInfo(bcs, ft0s, fdds, fv0as, bcsMatchedTrIds, bcsWithFIT);

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    // storing n-prong matches
    int32_t candID = 0;
    for (const auto& item : bcsMatchedTrIds) {
      uint64_t bc = item.first;
      std::vector<int32_t> fwdTrackIDs = item.second.first;
      std::vector<int32_t> barrelTrackIDs = item.second.second;
      int32_t nFwdTracks = fwdTrackIDs.size();
      int32_t nBarTracks = barrelTrackIDs.size();
      // check number of tracks in a candidate
      bool checkForward = nFwdTracks == fNFwdProngs;
      bool checkCentral = nBarTracks == fNBarProngs;
      // check TPC PID if needed
      std::vector<int32_t> filteredTrackIDs;
      if (fCheckTPCPID) {
        checkCentral = checkTPCPID(barrelTracks, barrelTrackIDs, filteredTrackIDs);
        if (checkCentral) {
          barrelTrackIDs.swap(filteredTrackIDs);
          nBarTracks = barrelTrackIDs.size();
        }
      }
      if (!checkForward || !checkCentral) {
        continue;
      }
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      uint16_t numContrib = nFwdTracks + nBarTracks;
      for (auto id : barrelTrackIDs) {
        const auto& tr = barrelTracks->iteratorAt(id);
        netCharge += tr.sign();
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = nBarTracks != 0 ? RgtrwTOF / static_cast<float>(nBarTracks) : 0.;
      if (RgtrwTOF == 0 && fNBarProngs != 0) { // require at least 1 TOF track in central and semiforward cases
        continue;
      }
      for (auto id : fwdTrackIDs) {
        const auto& tr = fwdTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      upchelpers::FITInfo fitInfo{};
      auto it = bcsWithFIT.find(bc);
      if (it != bcsWithFIT.end()) {
        fitInfo = it->second;
        if (fFilterFT0) { // check FT0 signal in the same BC
          bool hasNoFT0 = checkFT0(fitInfo, barrelTracks != nullptr && fwdTracks == nullptr);
          if (!hasNoFT0)
            continue;
        }
      }
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

    ambFwdTrIds.clear();
    ambBarrelTrIds.clear();
    filteredFwdTracks.clear();
    filteredBarrelTracks.clear();
    bcsMatchedTrIds.clear();
  }

  // data processors
  // _________________________________________________

  // create candidates for forward region
  void processFwd(ForwardTracks const& fwdTracks,
                  o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                  BCsWithBcSels const& bcs,
                  o2::aod::Collisions const& collisions,
                  o2::aod::FT0s const& ft0s,
                  o2::aod::FDDs const& fdds,
                  o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    fDoSemiFwd = false;
    createCandidates(&fwdTracks, static_cast<BarrelTracks*>(nullptr),
                     &ambFwdTracks, (o2::aod::AmbiguousTracks*)nullptr,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  // create candidates for semiforward region
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
    fDoSemiFwd = true;
    createCandidates(&fwdTracks, &barrelTracks,
                     &ambFwdTracks, &ambTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
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
    fDoSemiFwd = false;
    createCandidates(static_cast<ForwardTracks*>(nullptr), &barrelTracks,
                     (o2::aod::AmbiguousFwdTracks*)nullptr, &ambBarrelTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  // MC processors
  // _________________________________________________

  // create candidates for forward region
  void processFwdMC(ForwardTracks const& fwdTracks,
                    o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                    BCsWithBcSels const& bcs,
                    o2::aod::Collisions const& collisions,
                    o2::aod::FT0s const& ft0s,
                    o2::aod::FDDs const& fdds,
                    o2::aod::FV0As const& fv0as,
                    o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels)
  {
    fDoMC = true;
    fDoSemiFwd = false;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates(&fwdTracks, static_cast<BarrelTracks*>(nullptr),
                     &ambFwdTracks, (o2::aod::AmbiguousTracks*)nullptr,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     &mcFwdTrackLabels, (o2::aod::McTrackLabels*)nullptr);
    fNewPartIDs.clear();
  }

  // create candidates for semiforward region
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
    fDoSemiFwd = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates(&fwdTracks, &barrelTracks,
                     &ambFwdTracks, &ambTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     &mcFwdTrackLabels, &mcBarrelTrackLabels);
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
    fDoSemiFwd = false;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates(static_cast<ForwardTracks*>(nullptr), &barrelTracks,
                     (o2::aod::AmbiguousFwdTracks*)nullptr, &ambBarrelTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducer, processFwd, "Produce candidates for forward rapidities", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
  PROCESS_SWITCH(UpcCandProducer, processFwdMC, "Produce candidates for forward rapidities with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwdMC, "Produce candidates in semiforward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processCentralMC, "Produce candidates in central region with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
