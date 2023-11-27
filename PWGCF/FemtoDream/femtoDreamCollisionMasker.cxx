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

  void init(InitContext& context)
  {
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
    // for(auto i : TrackCutBits.at(0)){
    // 	LOG(info) << i;
    // }
  };

  void process(FDCollision const& col, FDParticles const& parts)
  {
    std::array<std::bitset<64>, CollisionMasks::kNParts> Mask = {{0}};
    for (auto& part : parts) {
      for (size_t nPart = 0; nPart < CollisionMasks::kNParts; nPart++) {
        for (size_t index = 0; index < TrackCutBits.at(nPart).size(); index++) {
          if ((part.cut() & TrackCutBits.at(nPart).at(index)) == TrackCutBits.at(nPart).at(index)) {
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
