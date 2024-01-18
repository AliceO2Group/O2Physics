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
/// \brief Tasks creates bitmasks for femtodream collisions
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include <vector>
#include <bitset>
#include <algorithm>

#include "fairlogger/Logger.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGCF/DataModel/FemtoDerived.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

namespace CollisionMasks
{
enum Parts {
  kPartOne,
  kPartTwo,
  kPartThree,
  kNParts,
};
enum Tasks {
  kTrackTrack,
  kTrackV0,
  kNTasks,
};
} // namespace CollisionMasks

struct femoDreamCollisionMasker {
  Produces<FDColMasks> Masks;

  // particle selection bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackCutBits;
  // particle tpc pid bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackPIDTPCBits;
  // particle tpctof pid bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackPIDTPCTOFBits;
  // particle momemtum threshold for PID
  std::array<std::vector<float>, CollisionMasks::kNParts> TrackPIDThreshold;

  // particle selection for v0
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> V0CutBits;

  // particle selection bits for v0 children
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> PosChildCutBits;
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> NegChildCutBits;

  // particle tpc pid bits for v0 children
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> PosChildPIDTPCBits;
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> NegChildPIDTPCBits;

  // cuts from filters on the pair task
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterPtMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterPtMax;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterEtaMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterEtaMax;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassMax;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassAntiMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassAntiMax;

  int TaskFinder = -1;

  void init(InitContext& context)
  {
    std::vector<std::string> MatchedWorkflows;
    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";
    LOG(info) << " Collision masker self-configuration ";
    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-track") != std::string::npos) {
        LOG(info) << "Matched workflow: " << device.name;
        TaskFinder = CollisionMasks::kTrackTrack;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("ConfTrk1_CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_minPt")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_maxPt")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_minEta")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_maxEta")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk2_CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk2_TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk2_TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk2_PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk2_minPt")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk2_maxPt")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk2_minEta")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk2_maxEta")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          }
        }
      } else if (device.name.find(std::string("femto-dream-pair-task-track-v0")) != std::string::npos) {
        LOG(info) << "Matched workflow: " << device.name;
        TaskFinder = CollisionMasks::kTrackV0;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("ConfTrk1_CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTrk1_PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_minPt")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_maxPt")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_minEta")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfTrk1_maxEta")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_CutBit")) == 0) {
            V0CutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfV02_minPt")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_maxPt")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_minEta")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_maxEta")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_minInvMass")) == 0) {
            FilterInvMassMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_maxInvMass")) == 0) {
            FilterInvMassMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_minInvMassAnti")) == 0) {
            FilterInvMassAntiMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_maxInvMassAnti")) == 0) {
            FilterInvMassAntiMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfV02_ChildPos_CutBit")) == 0) {
            PosChildCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfV02_ChildPos_TPCBit")) == 0) {
            PosChildPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfV02_ChildNeg_CutBit")) == 0) {
            NegChildCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfV02_ChildNeg_TPCBit")) == 0) {
            NegChildPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          }
        }
      }
    }

    if (8 * sizeof(femtodreamcollision::BitMaskType) < TrackCutBits.at(CollisionMasks::kPartOne).size() ||
        8 * sizeof(femtodreamcollision::BitMaskType) < TrackCutBits.at(CollisionMasks::kPartTwo).size() ||
        8 * sizeof(femtodreamcollision::BitMaskType) < TrackCutBits.at(CollisionMasks::kPartThree).size()) {
      LOG(fatal) << "Too many variations!";
      LOG(fatal) << "Collision masker only supports up to " << 8 * sizeof(femtodreamcollision::BitMaskType) << " variations!";
    }
  }

  // make bitmask for a track
  template <typename T, typename R>
  void MaskForTrack(T& BitSet, CollisionMasks::Parts P, R& track)
  {
    if (track.partType() != static_cast<uint8_t>(femtodreamparticle::kTrack)) {
      return;
    }
    for (size_t index = 0; index < TrackCutBits.at(P).size(); index++) {
      // check filter cuts
      if (track.pt() < FilterPtMin.at(P).at(index) || track.pt() > FilterPtMax.at(P).at(index) ||
          track.eta() < FilterEtaMin.at(P).at(index) || track.eta() > FilterEtaMax.at(P).at(index)) {
        // if they are not passed, skip the particle
        continue;
      }
      // set the bit at the index of the selection equal to one if the track passes all selections
      // check track cuts
      if ((track.cut() & TrackCutBits.at(P).at(index)) == TrackCutBits.at(P).at(index)) {
        // check pid cuts
        if (track.p() <= TrackPIDThreshold.at(P).at(index)) {
          if ((track.pidcut() & TrackPIDTPCBits.at(P).at(index)) == TrackPIDTPCBits.at(P).at(index)) {
            BitSet.at(P).set(index);
          }
        } else {
          if ((track.pidcut() & TrackPIDTPCTOFBits.at(P).at(index)) == TrackPIDTPCTOFBits.at(P).at(index)) {
            BitSet.at(P).set(index);
          }
        }
      }
    }
  }

  // make bit mask for v0
  template <typename T, typename R, typename S>
  void MaskForV0(T& BitSet, CollisionMasks::Parts P, R& v0, S& parts)
  {
    if (v0.partType() != static_cast<uint8_t>(femtodreamparticle::kV0)) {
      return;
    }
    for (size_t index = 0; index < V0CutBits.at(P).size(); index++) {
      // check filter cuts
      if (v0.pt() < FilterPtMin.at(P).at(index) || v0.pt() > FilterPtMax.at(P).at(index) ||
          v0.eta() < FilterEtaMin.at(P).at(index) || v0.eta() > FilterEtaMax.at(P).at(index) ||
          v0.mLambda() < FilterInvMassMin.at(P).at(index) || v0.mLambda() > FilterInvMassMax.at(P).at(index) ||
          v0.mAntiLambda() < FilterInvMassAntiMin.at(P).at(index) || v0.mAntiLambda() > FilterInvMassAntiMax.at(P).at(index)) {
        // if they are not passed, skip the particle
        continue;
      }
      // check cut bit of v0
      if ((v0.cut() & V0CutBits.at(P).at(index)) == V0CutBits.at(P).at(index)) {
        const auto& posChild = parts.iteratorAt(v0.index() - 2);
        const auto& negChild = parts.iteratorAt(v0.index() - 1);
        // This is how it is supposed to work but there seems to be an issue
        // for now, keep in sync with femtodreampairtasktrackv0
        // auto posChild = v0.template children_as<FDParticles>().front();
        // auto negChild = v0.template children_as<FDParticles>().back();
        // check cut on v0 children
        if ((posChild.cut() & PosChildCutBits.at(P).at(index)) == PosChildCutBits.at(P).at(index) &&
            (posChild.pidcut() & PosChildPIDTPCBits.at(P).at(index)) == PosChildPIDTPCBits.at(P).at(index) &&
            (negChild.cut() & NegChildCutBits.at(P).at(index)) == NegChildCutBits.at(P).at(index) &&
            (negChild.pidcut() & NegChildPIDTPCBits.at(P).at(index)) == NegChildPIDTPCBits.at(P).at(index)) {
          BitSet.at(P).set(index);
        }
      }
    }
  }

  void process(FDCollision const& col, FDParticles const& parts)
  {
    // create a bit mask for particle one, particle two and particle three
    std::array<std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)>, CollisionMasks::kNParts> Mask = {{0}};

    switch (TaskFinder) {
      case CollisionMasks::kTrackTrack:
        // pair-track-track task
        // create mask for track 1 and track 2
        for (auto const& part : parts) {
          MaskForTrack(Mask, CollisionMasks::kPartOne, part);
          MaskForTrack(Mask, CollisionMasks::kPartTwo, part);
        }
        break;
      case CollisionMasks::kTrackV0:
        // pair-track-v0 task
        // create mask for track 1 and v0 2
        for (auto const& part : parts) {
          MaskForTrack(Mask, CollisionMasks::kPartOne, part);
          MaskForV0(Mask, CollisionMasks::kPartTwo, part, parts);
        }
        // TODO: add all supported pair/triplet tasks
        break;
      default:
        LOG(fatal) << "No femtodream pair task found!";
    }
    // fill bitmask for each collision
    Masks(static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartOne).to_ulong()),
          static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartTwo).to_ulong()),
          static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartThree).to_ulong()));
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femoDreamCollisionMasker>(cfgc)};
  return workflow;
};
