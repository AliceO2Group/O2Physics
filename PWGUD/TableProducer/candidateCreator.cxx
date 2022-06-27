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

#include "CommonConstants/LHCConstants.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct CandidateCreator {
  Produces<o2::aod::EventCandidates> eventCandidates;

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

  template <typename TFwdTracks, typename TBarTracks,
            typename TBCs, typename TFT0s>
  void createCandidates(TFwdTracks* fwdTracks,
                        TBarTracks* barTracks,
                        TBCs const& bcs,
                        TFT0s const& ft0s)
  {
    std::map<uint64_t, int32_t> BCsWithFT0;
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      BCsWithFT0[bc] = ft0.globalIndex();
    }

    std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> bcsMatchedFwdTrIds; // pairs of global BCs and vectors of matched track IDs

    // forward matching
    if (fwdTracks != nullptr) {
      for (const auto& fwdTr : *fwdTracks) {
        uint64_t bc = std::trunc(fwdTr.trackTime() / o2::constants::lhc::LHCBunchSpacingNS); // absolute time to BC
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
      for (const auto& barTr : *barTracks) {
        uint64_t bc = std::trunc(barTr.trackTime() / o2::constants::lhc::LHCBunchSpacingNS); // absolute time to BC
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

    // storing n-prong matches
    for (const auto& item : bcsMatchedFwdTrIds) {
      uint64_t bc = item.first;
      const std::vector<int32_t>& fwdTrackIDs = item.second.first;
      const std::vector<int32_t>& barTrackIDs = item.second.second;
      int32_t nFwdTracks = fwdTrackIDs.size();
      int32_t nBarTracks = barTrackIDs.size();
      if (!(nFwdTracks == fNFwdProngs && nBarTracks == fNBarProngs)) {
        continue;
      }
      // fetching FT0 information
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
      eventCandidates(bc, runNumber, fwdTrackIDs, barTrackIDs,
                      ft0Info.amplitudeA, ft0Info.amplitudeC, ft0Info.timeA, ft0Info.timeC, ft0Info.triggerMask);
    }
    bcsMatchedFwdTrIds.clear();
    BCsWithFT0.clear();
  }

  void processFwd(o2::aod::SkimmedMuons const& muonTracks,
                  o2::aod::BCs const& bcs,
                  o2::aod::FT0s const& ft0s)
  {
    createCandidates(&muonTracks, (o2::aod::SkimmedBarTracks*)nullptr, bcs, ft0s);
  }

  void processSemiFwd(o2::aod::SkimmedMuons const& muonTracks,
                      o2::aod::SkimmedBarTracks const& barTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::FT0s const& ft0s)
  {
    fDoSemiFwd = true;
    createCandidates(&muonTracks, &barTracks, bcs, ft0s);
  }

  PROCESS_SWITCH(CandidateCreator, processFwd, "Produce candidates for forward rapidities", false);
  PROCESS_SWITCH(CandidateCreator, processSemiFwd, "Produce candidates in semiforward region", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CandidateCreator>(cfgc)};
}
