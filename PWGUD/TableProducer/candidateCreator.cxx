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

  template <uint32_t TNFwdProngs, uint32_t TNBarProngs, typename TFwdTracks, typename TBarTracks,
            typename TBCs, typename TFT0s>
  void createCandidates(TFwdTracks* fwdTracks, TBarTracks* barTracks,
                        TBCs const& bcs, TFT0s const& ft0s)
  {
    std::map<uint64_t, int32_t> BCsWithFT0;
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      BCsWithFT0[bc] = ft0.globalIndex();
    }

    std::vector<std::pair<int64_t, std::vector<int32_t>>> bcsMatchedFwdTrIds; // pairs of global BCs and vectors of matched track IDs

    // forward matching
    if constexpr (TNFwdProngs > 0 && TNBarProngs == 0) {
      for (const auto& fwdTr : *fwdTracks) {
        uint64_t bc = std::trunc(fwdTr.trackTime() / o2::constants::lhc::LHCBunchSpacingNS); // absolute time to BC
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = std::find_if(bcsMatchedFwdTrIds.begin(), bcsMatchedFwdTrIds.end(),
                               [=](const std::pair<uint64_t, std::vector<int32_t>>& element) { return element.first == bc; });
        if (it != bcsMatchedFwdTrIds.end()) {
          it->second.emplace_back(fwdTr.globalIndex());
        } else {
          bcsMatchedFwdTrIds.emplace_back(bc, std::vector<int32_t>(1, fwdTr.globalIndex()));
        }
      }
      // storing n-prong matches
      std::vector<int32_t> dummyBarTrackIDs;
      for (const auto& item : bcsMatchedFwdTrIds) {
        uint64_t bc = item.first;
        const std::vector<int32_t>& fwdTrackIDs = item.second;
        size_t nTracks = fwdTrackIDs.size();
        if (nTracks != TNFwdProngs) {
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
        eventCandidates(bc, runNumber, fwdTrackIDs, dummyBarTrackIDs,
                        ft0Info.amplitudeA, ft0Info.amplitudeC, ft0Info.timeA, ft0Info.timeC, ft0Info.triggerMask);
      }
      bcsMatchedFwdTrIds.clear();
      BCsWithFT0.clear();
    }

    if constexpr (TNFwdProngs > 0 && TNBarProngs > 0) {
      // todo: semiforward matching
    }

    if constexpr (TNFwdProngs == 0 && TNBarProngs > 0) {
      // todo: central barrel matching
    }
  }

  void process(o2::aod::SkimmedMuons const& muonTracks,
               o2::aod::BCs const& bcs,
               o2::aod::FT0s const& ft0s)
  {
    if (fNBarProngs <= 0 && fNFwdProngs > 0) {
      switch (fNFwdProngs) {
        case 1:
          createCandidates<1, 0>(&muonTracks, (o2::aod::Tracks*)nullptr, bcs, ft0s);
          break;
        case 2:
          createCandidates<2, 0>(&muonTracks, (o2::aod::Tracks*)nullptr, bcs, ft0s);
          break;
        case 3:
          createCandidates<3, 0>(&muonTracks, (o2::aod::Tracks*)nullptr, bcs, ft0s);
          break;
        case 4:
          createCandidates<4, 0>(&muonTracks, (o2::aod::Tracks*)nullptr, bcs, ft0s);
          break;
        default:
          LOGF(warn, "Unsupported number of forward prongs! Available numbers: 1, 2, 3, 4");
          return;
      }
    }
    if (fNBarProngs > 0 && fNFwdProngs <= 0) {
      LOGF(warn, "Central barrel is not supported yet");
      return;
    }
    if (fNBarProngs > 0 && fNFwdProngs > 0) {
      LOGF(warn, "Semiforward analysis is not supported yet");
      return;
    }
    if (fNBarProngs <= 0 && fNFwdProngs <= 0) {
      LOGF(warn, "Wrong number of prongs! Check inputs");
      return;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CandidateCreator>(cfgc)};
}
