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

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducer {
  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDTrackCollisionIDs> barrelCandIds;
  Produces<o2::aod::UDFwdTrackCollisionIDs> muonCandIds;

  bool fDoSemiFwd{false};

  Configurable<int> fNFwdProngs{"nFwdProngs", 2, "Matched forward tracks per candidate"};
  Configurable<int> fNBarProngs{"nBarProngs", 0, "Matched barrel tracks per candidate"};

  // helper struct
  struct FT0Info {
    float amplitudeA = -1;
    float amplitudeC = -1;
    float timeA = -999.;
    float timeC = -999.;
    uint8_t triggerMask = 0;
  };

  template <typename TFwdTracks, typename TBarrelTracks,
            typename TBCs, typename TFT0s>
  void createCandidates(TFwdTracks* fwdTracks,
                        TBarrelTracks* barTracks,
                        TBCs const& bcs,
                        TFT0s const& ft0s)
  {
    // map track IDs to the respective event candidate IDs
    std::vector<int32_t> barTrackCandIds;
    std::vector<int32_t> fwdTrackCandIds;

    std::map<uint64_t, int32_t> BCsWithFT0;
    // collect BCs with FT0 signals
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      BCsWithFT0[bc] = ft0.globalIndex();
    }

    // pairs of global BCs and vectors of matched track IDs:
    // global BC <-> <vector of fwd. trackIDs, vector of barrel trackIDs>
    std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> bcsMatchedFwdTrIds;

    // forward matching
    if (fwdTracks != nullptr) {
      fwdTrackCandIds.resize(fwdTracks->size(), -1);
      for (const auto& fwdTr : *fwdTracks) {
        uint64_t bc = fwdTr.globalBC();
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedFwdTrIds.find(bc);
        if (it != bcsMatchedFwdTrIds.end()) {
          it->second.first.emplace_back(fwdTr.globalIndex());
        } else {
          bcsMatchedFwdTrIds[bc] = std::make_pair(std::vector<int32_t>(1, fwdTr.globalIndex()), std::vector<int32_t>());
        }
      }
    }

    // central barrel tracks
    if (barTracks != nullptr) {
      barTrackCandIds.resize(barTracks->size(), -1);
      for (const auto& barTr : *barTracks) {
        uint64_t bc = barTr.globalBC();
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedFwdTrIds.find(bc);
        if (it != bcsMatchedFwdTrIds.end()) {
          it->second.second.emplace_back(barTr.globalIndex());
        } else if (!fDoSemiFwd) { // tag central-barrel tracks to forward tracks in semiforward case
          bcsMatchedFwdTrIds[bc] = std::make_pair(std::vector<int32_t>(), std::vector<int32_t>(1, barTr.globalIndex()));
        }
      }
    }

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    // storing n-prong matches
    int32_t candID = 0;
    for (const auto& item : bcsMatchedFwdTrIds) {
      uint64_t bc = item.first;
      std::vector<int32_t> fwdTrackIDs = item.second.first;
      std::vector<int32_t> barTrackIDs = item.second.second;
      int32_t nFwdTracks = fwdTrackIDs.size();
      int32_t nBarTracks = barTrackIDs.size();
      // check number of tracks in a candidate
      bool checkForward = nFwdTracks == fNFwdProngs;
      bool checkCentral = nBarTracks == fNBarProngs;
      // check TPC PID if needed
      if (!checkForward || !checkCentral) {
        continue;
      }
      float RgtrwTOF = 0.;
      for (auto id : barTrackIDs) {
        const auto& tr = barTracks->iteratorAt(id);
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = nBarTracks != 0 ? RgtrwTOF / (float)nBarTracks : 0.;
      if (RgtrwTOF == 0 && fNBarProngs != 0) { // require at least 1 TOF track in central and semiforward cases
        continue;
      }
      int8_t netCharge = 0;
      uint16_t numContrib = nFwdTracks + nBarTracks;
      for (auto id : fwdTrackIDs) {
        fwdTrackCandIds[id] = candID;
        const auto& tr = fwdTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      for (auto id : barTrackIDs) {
        barTrackCandIds[id] = candID;
        const auto& tr = barTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      // fetching FT0 information
      // if there is no FT0 signal, dummy info will be used
      FT0Info ft0Info;
      auto ft0Iter = BCsWithFT0.find(bc);
      if (ft0Iter != BCsWithFT0.end()) {
        const auto& ft0 = ft0s.iteratorAt(ft0Iter->second);
        const auto& ampsA = ft0.amplitudeA();
        const auto& ampsC = ft0.amplitudeC();
        ft0Info.amplitudeA = 0.;
        for (auto amp : ampsA) {
          ft0Info.amplitudeA += amp;
        }
        ft0Info.amplitudeC = 0.;
        for (auto amp : ampsC) {
          ft0Info.amplitudeC += amp;
        }
        ft0Info.timeA = ft0.timeA();
        ft0Info.timeC = ft0.timeC();
        ft0Info.triggerMask = ft0.triggerMask();
      }
      int32_t runNumber = bcs.iteratorAt(0).runNumber();
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF,
                      ft0Info.amplitudeA, ft0Info.amplitudeC, ft0Info.timeA, ft0Info.timeC, ft0Info.triggerMask);
      candID++;
    }
    bcsMatchedFwdTrIds.clear();
    BCsWithFT0.clear();

    if (fwdTracks != nullptr) {
      for (const auto& fwdTr : *fwdTracks) {
        muonCandIds(fwdTrackCandIds[fwdTr.globalIndex()]);
      }
    }

    fwdTrackCandIds.clear();

    if (barTracks != nullptr) {
      for (const auto& barTr : *barTracks) {
        barrelCandIds(barTrackCandIds[barTr.globalIndex()]);
      }
    }

    barTrackCandIds.clear();
  }

  using BarrelTracks = o2::soa::Join<o2::aod::UDTracks, o2::aod::UDTracksExtra, o2::aod::UDTracksPID>;

  // create candidates for forward region
  void processFwd(o2::aod::UDFwdTracks const& muonTracks,
                  o2::aod::BCs const& bcs,
                  o2::aod::FT0s const& ft0s)
  {
    fDoSemiFwd = false;
    createCandidates(&muonTracks, (BarrelTracks*)nullptr, bcs, ft0s);
  }

  // create candidates for semiforward region
  void processSemiFwd(o2::aod::UDFwdTracks const& muonTracks,
                      BarrelTracks const& barTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::FT0s const& ft0s)
  {
    fDoSemiFwd = true;
    createCandidates(&muonTracks, &barTracks, bcs, ft0s);
  }

  // create candidates for central region
  void processCentral(o2::aod::UDFwdTracks const& muonTracks,
                      BarrelTracks const& barTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::FT0s const& ft0s)
  {
    fDoSemiFwd = false;
    createCandidates((o2::aod::UDFwdTracks*)nullptr, &barTracks, bcs, ft0s);
  }

  PROCESS_SWITCH(UpcCandProducer, processFwd, "Produce candidates for forward rapidities", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
