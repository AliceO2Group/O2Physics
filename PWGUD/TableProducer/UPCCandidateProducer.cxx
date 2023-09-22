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
                                     o2::aod::TOFSignal, o2::aod::pidTOFbeta,
                                     o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

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
      if (!mcPart.has_mcCollision())
        continue;
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
      udMCParticles(newEventID, mcPart.pdgCode(), mcPart.getHepMCStatusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
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
                        std::unordered_map<int64_t, uint64_t>& ambBarrelTrBCs)
  {
    for (auto trackID : trackIDs) {
      const auto& track = tracks.iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      int64_t colId = -1;
      if (ambBarrelTrBCs.find(trackID) == ambBarrelTrBCs.end()) {
        colId = track.collisionId();
      }
      udTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udTracksExtra(track.tpcInnerParam(), track.itsClusterMap(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
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

  void processFITInfo(upchelpers::FITInfo& fitInfo,
                      uint64_t midbc,
                      std::vector<std::pair<uint64_t, int64_t>>& v,
                      BCsWithBcSels const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
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
  template <int32_t tracksSwitch, typename TAmbTracks>
  void collectAmbTrackBCs(std::unordered_map<int64_t, uint64_t>& ambTrIds,
                          TAmbTracks ambTracks)
  {
    for (const auto& ambTrk : ambTracks) {
      auto trkId = getAmbTrackId<tracksSwitch>(ambTrk);
      const auto& bcSlice = ambTrk.bc();
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

  void collectBarrelTracks(std::vector<BCTracksPair>& bcsMatchedTrIdsA,
                           std::vector<BCTracksPair>& bcsMatchedTrIdsB,
                           BCsWithBcSels const& bcs,
                           o2::aod::Collisions const& collisions,
                           BarrelTracks const& barrelTracks,
                           o2::aod::AmbiguousTracks const& ambBarrelTracks,
                           std::unordered_map<int64_t, uint64_t>& ambBarrelTrBCs)
  {
    for (const auto& trk : barrelTracks) {
      if (!applyBarCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int32_t nContrib = -1;
      uint64_t trackBC = 0;
      auto ambIter = ambBarrelTrBCs.find(trkId);
      if (ambIter == ambBarrelTrBCs.end()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
        trackBC = col.bc_as<BCsWithBcSels>().globalBC();
      } else {
        trackBC = ambIter->second;
      }
      int64_t tint = TMath::FloorNint(trk.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
      uint64_t bc = trackBC + tint;
      if (bc > fMaxBC)
        continue;
      bool checkNContrib = nContrib <= upcCuts.getMaxNContrib();
      if (!checkNContrib)
        continue;
      bool needITSITS = upcCuts.getProduceITSITS() && (trk.hasTOF() || (trk.hasITS() && trk.hasTPC()));
      bool needAllTOF = !upcCuts.getProduceITSITS() && !upcCuts.getRequireITSTPC() && trk.hasTOF();
      bool needTOFWithITS = !upcCuts.getProduceITSITS() && upcCuts.getRequireITSTPC() && trk.hasTOF() && trk.hasITS() && trk.hasTPC();
      bool addToA = needITSITS || needAllTOF || needTOFWithITS;
      if (addToA)
        addTrack(bcsMatchedTrIdsA, bc, trkId);
      if (fSearchITSTPC == 1 && !trk.hasTOF() && trk.hasITS() && trk.hasTPC())
        addTrack(bcsMatchedTrIdsB, bc, trkId);
    }
  }

  void collectForwardTracks(std::vector<BCTracksPair>& bcsMatchedTrIdsMID,
                            BCsWithBcSels const& bcs,
                            o2::aod::Collisions const& collisions,
                            ForwardTracks const& fwdTracks,
                            o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                            std::unordered_map<int64_t, uint64_t>& ambFwdTrBCs)
  {
    for (const auto& trk : fwdTracks) {
      if (upcCuts.getTrackType() != -1) {  // need to filter by type ?
        if (upcCuts.getTrackType() == 5) { // is MCH-MID or MCH ?
          if (trk.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack &&
              trk.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack)
            continue;
        }
        if (trk.trackType() != upcCuts.getTrackType()) // getTrackType < 5 -> check exact track type
          continue;
      }
      if (!applyFwdCuts(trk))
        continue;
      int64_t trkId = trk.globalIndex();
      int32_t nContrib = -1;
      uint64_t trackBC = 0;
      auto ambIter = ambFwdTrBCs.find(trkId);
      if (ambIter == ambFwdTrBCs.end()) {
        const auto& col = trk.collision();
        nContrib = col.numContrib();
        trackBC = col.bc_as<BCsWithBcSels>().globalBC();
      } else {
        trackBC = ambIter->second;
      }
      int64_t tint = TMath::FloorNint(trk.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
      uint64_t bc = trackBC + tint;
      if (nContrib <= upcCuts.getMaxNContrib())
        addTrack(bcsMatchedTrIdsMID, bc, trkId);
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
    std::unordered_map<int64_t, uint64_t> ambBarrelTrBCs;
    collectAmbTrackBCs<0>(ambBarrelTrBCs, ambBarrelTracks);

    collectBarrelTracks(bcsMatchedTrIdsTOF, bcsMatchedTrIdsITSTPC,
                        bcs, collisions,
                        barrelTracks, ambBarrelTracks, ambBarrelTrBCs);

    uint32_t nBCsWithITSTPC = bcsMatchedTrIdsITSTPC.size();

    std::sort(bcsMatchedTrIdsTOF.begin(), bcsMatchedTrIdsTOF.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });
    std::sort(bcsMatchedTrIdsITSTPC.begin(), bcsMatchedTrIdsITSTPC.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

    if (nBCsWithITSTPC > 0 && fSearchITSTPC == 1) {
      std::unordered_set<int64_t> matchedTracks;
      for (auto& pair : bcsMatchedTrIdsTOF) {
        uint64_t bc = pair.first;
        auto& trackIds = pair.second;
        uint32_t nTOFtracks = trackIds.size();
        if (nTOFtracks > fNBarProngs) // too many TOF tracks?!
          continue;
        if (nTOFtracks == fNBarProngs) { // check for ITS-TPC tracks
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, 0, tracks, bcsMatchedTrIdsITSTPC, matchedTracks, true);
          if (res < 0) { // too many tracks nearby -> rejecting
            trackIds.push_back(0);
            continue;
          }
        }
        if (nTOFtracks < fNBarProngs && !upcCuts.getRequireTOF()) { // add ITS-TPC track if needed
          uint32_t tracksToFind = fNBarProngs - nTOFtracks;
          std::vector<int64_t> tracks;
          tracks.reserve(fNBarProngs * 2); // precautions
          int32_t res = searchTracks(bc, fSearchRangeITSTPC, tracksToFind, tracks, bcsMatchedTrIdsITSTPC, matchedTracks, true);
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
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels, ambBarrelTrBCs);
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
    ambBarrelTrBCs.clear();
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
    collectAmbTrackBCs<0>(ambBarrelTrBCs, ambBarrelTracks);

    std::unordered_map<int64_t, uint64_t> ambFwdTrBCs;
    collectAmbTrackBCs<1>(ambFwdTrBCs, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    collectBarrelTracks(bcsMatchedTrIdsTOF, bcsMatchedTrIdsITSTPC,
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
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, mcFwdTrackLabels);
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels, ambBarrelTrBCs);
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
    ambFwdTrBCs.clear();
    bcsMatchedTrIdsMID.clear();
    ambBarrelTrBCs.clear();
    bcsMatchedTrIdsTOFTagged.clear();
  }

  void createCandidatesFwd(ForwardTracks const& fwdTracks,
                           o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                           BCsWithBcSels const& bcs,
                           o2::aod::Collisions const& collisions,
                           o2::aod::FT0s const& ft0s,
                           o2::aod::FDDs const& fdds,
                           o2::aod::FV0As const& fv0as,
                           const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    fMaxBC = bcs.iteratorAt(bcs.size() - 1).globalBC(); // restrict ITS-TPC track search to [0, fMaxBC]

    // pairs of global BCs and vectors of matched track IDs:
    std::vector<BCTracksPair> bcsMatchedTrIdsMID;

    // trackID -> index in amb. track table
    std::unordered_map<int64_t, uint64_t> ambFwdTrBCs;
    collectAmbTrackBCs<1>(ambFwdTrBCs, ambFwdTracks);

    collectForwardTracks(bcsMatchedTrIdsMID,
                         bcs, collisions,
                         fwdTracks, ambFwdTracks, ambFwdTrBCs);

    uint32_t nBCsWithMID = bcsMatchedTrIdsMID.size();

    std::sort(bcsMatchedTrIdsMID.begin(), bcsMatchedTrIdsMID.end(),
              [](const auto& left, const auto& right) { return left.first < right.first; });

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
      auto& fwdTrackIDs = pairMID.second;
      uint32_t nMIDtracks = fwdTrackIDs.size();
      uint16_t numContrib = nMIDtracks;
      uint64_t bc = pairMID.first;
      // sanity check
      if (nMIDtracks != fNFwdProngs) {
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
      for (auto id : fwdTrackIDs) {
        const auto& tr = fwdTracks.iteratorAt(id);
        netCharge += tr.sign();
      }
      // store used tracks
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, mcFwdTrackLabels);
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
    ambFwdTrBCs.clear();
    bcsMatchedTrIdsMID.clear();
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

  // create candidates for forward region
  // forward: n fwd tracks
  void processForward(ForwardTracks const& fwdTracks,
                      o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    createCandidatesFwd(fwdTracks, ambFwdTracks,
                        bcs, collisions,
                        ft0s, fdds, fv0as,
                        (o2::aod::McFwdTrackLabels*)nullptr);
  }

  void processForwardMC(ForwardTracks const& fwdTracks,
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
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidatesFwd(fwdTracks, ambFwdTracks,
                        bcs, collisions,
                        ft0s, fdds, fv0as,
                        &mcFwdTrackLabels);
  }

  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward/forward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwdMC, "Produce candidates in semiforward/forward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processCentralMC, "Produce candidates in central region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processForward, "Produce caniddates in forward region", false);
  PROCESS_SWITCH(UpcCandProducer, processForwardMC, "Produce caniddates in forward region with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
