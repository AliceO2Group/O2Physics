// Copyright 2019-2023 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtoDreamCollisionMasker.cxx
/// \brief Tasks creates bitmasks for collisions
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <Framework/Configurable.h>
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include <cstdint>
#include <vector>
#include <algorithm>
#include <bitset>
#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2::aod;
using namespace o2::framework;

namespace CollisionMasks
{
enum Parts {
  kPartOne,
  kPartTwo,
  kPartThree,
  kNParts
};
}

struct femtoDreamCollisionMasker {
  Produces<FDColMasks> Masks;
  /// track selection bits
  std::array<std::vector<uint32_t>, CollisionMasks::kNParts> TrackCutBits = {};
  /// track tpc pid bits
  std::array<std::vector<uint32_t>, CollisionMasks::kNParts> TrackPIDTPCBits = {};
  /// track tpctof pid bits
  std::array<std::vector<uint32_t>, CollisionMasks::kNParts> TrackPIDTPCTOFBits = {};
  /// track momemtum threshold for PID
  std::array<std::vector<float>, CollisionMasks::kNParts> TrackPIDThreshold = {};

  // Configurable<bool> ConfSelfConfiguration{"ConfSelfConfiguration", true, "Flag to active self-configuration"};

  // particle 1
  // Configurable<std::vector<uint32_t>> ConfTrackCutBitsPart1{"ConfTrackCutsBitsPart1", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackCutBitsPart2{"ConfTrackCutsBitsPart2", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackCutBitsPart3{"ConfTrackCutsBitsPart3", {1}, "Cut bits for collision masking"};
  //
  // // particle 2
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCBitsPart1{"ConfTrackPIDTPCBitsPart1", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCBitsPart2{"ConfTrackPIDTPCBitsPart2", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCBitsPart3{"ConfTrackPIDTPCBitsPart3", {1}, "Cut bits for collision masking"};
  //
  // // particle 3
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCTOFBitsPart1{"ConfTrackPIDTPCTOFBitsPart1", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCTOFBitsPart2{"ConfTrackPIDTPCTOFBitsPart2", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<uint32_t>> ConfTrackPIDTPCTOFBitsPart3{"ConfTrackPIDTPCTOFBitsPart3", {1}, "Cut bits for collision masking"};
  // Configurable<std::vector<float>> ConfTrackPIDThresholdsPart1{"ConTrackPIDThresholdPart1", {0.75}, "Cut bits for collision masking"};
  // Configurable<std::vector<float>> ConfTrackPIDThresholdsPart2{"ConTrackPIDThresholdPart2", {0.75}, "Cut bits for collision masking"};
  // Configurable<std::vector<float>> ConfTrackPIDThresholdsPart3{"ConTrackPIDThresholdPart3", {0.75}, "Cut bits for collision masking"};

  void init(InitContext& context)
  {
		// asdf
    // if (ConfSelfConfiguration.value) {
      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      LOGF(info, " Collision masker self-configuration ");
      LOGF(info, "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*");
      auto& workflows = context.services().get<RunningWorkflowInfo const>();
      for (DeviceSpec const& device : workflows.devices) {
        if (device.name.compare("femto-dream-pair-task-track-track") == 0) {
          for (auto const& option : device.options) {
            LOG(info) << option.name;
            if (option.name.compare(std::string("ConfCutPartOne")) == 0) {
              TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDTPCPartOne")) == 0) {
              TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDTPCTOFPartOne")) == 0) {
              TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDThresPartOne")) == 0) {
              TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
            } else if (option.name.compare(std::string("ConfCutPartTwo")) == 0) {
              TrackCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDTPCPartTwo")) == 0) {
              TrackPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDTPCTOFPartTwo")) == 0) {
              TrackPIDTPCTOFBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<uint32_t>());
            } else if (option.name.compare(std::string("ConfPIDThresPartTwo")) == 0) {
              TrackPIDThreshold.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
            }
          }
        }
      }
    // } else {
    //   TrackCutBits = std::array<std::vector<uint32_t>, CollisionMasks::kNParts>{ConfTrackCutBitsPart1.value, ConfTrackCutBitsPart2.value, ConfTrackCutBitsPart3.value};
    //   TrackPIDTPCBits = std::array<std::vector<uint32_t>, CollisionMasks::kNParts>{ConfTrackPIDTPCBitsPart1.value, ConfTrackPIDTPCBitsPart2.value, ConfTrackPIDTPCBitsPart3.value};
    //   TrackPIDTPCTOFBits = std::array<std::vector<uint32_t>, CollisionMasks::kNParts>{ConfTrackPIDTPCTOFBitsPart1.value, ConfTrackPIDTPCTOFBitsPart2.value, ConfTrackPIDTPCTOFBitsPart3.value};
    //   TrackPIDThreshold = std::array<std::vector<float>, CollisionMasks::kNParts>{ConfTrackPIDThresholdsPart1.value, ConfTrackPIDThresholdsPart2.value, ConfTrackPIDThresholdsPart3.value};
    // }
  };

  void processTrack(FDCollision const& col, FDParticles const& parts)
  {
    // create a bit mask for particle one, particle two and particle three
    std::array<std::bitset<64>, CollisionMasks::kNParts> Mask = {{0}};
    // iterate over all tracks in this collision
    for (auto& part : parts) {
      // iterate over particleOne/Two/Three
      for (size_t nPart = 0; nPart < CollisionMasks::kNParts; nPart++) {
        // iterate over all selections for particleOne/Two/Three
        for (size_t index = 0; index < TrackCutBits.at(nPart).size(); index++) {
        // set the bit at the index of the selection equal to one if the track passes
          // check track cuts
          if ((part.cut() & TrackCutBits.at(nPart).at(index)) == TrackCutBits.at(nPart).at(index)) {
            // check pid cuts
            if (part.p() < TrackPIDThreshold.at(nPart).at(index)) {
              if ((part.cut() & TrackPIDTPCBits.at(nPart).at(index)) == TrackPIDTPCBits.at(nPart).at(index)) {
                Mask.at(nPart).set(index);
              }
            } else {
              if ((part.cut() & TrackPIDTPCTOFBits.at(nPart).at(index)) == TrackPIDTPCTOFBits.at(nPart).at(index)) {
                Mask.at(nPart).set(index);
              }
            }
          }
        }
      }
    }
    // LOG(info) << "Part1:" << Mask.at(CollisionMasks::kPartOne).to_string();
    // LOG(info) << "Part2:" << Mask.at(CollisionMasks::kPartTwo).to_string();
    // LOG(info) << "Part3:" << Mask.at(CollisionMasks::kPartThree).to_string();
    // fill bitmask for each collision
    Masks(static_cast<uint32_t>(Mask.at(CollisionMasks::kPartOne).to_ulong()),
          static_cast<uint32_t>(Mask.at(CollisionMasks::kPartTwo).to_ulong()),
          static_cast<uint32_t>(Mask.at(CollisionMasks::kPartThree).to_ulong()));
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femtoDreamCollisionMasker>(cfgc)};
  return workflow;
};
