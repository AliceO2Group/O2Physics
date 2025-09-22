// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
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
/// \author Laura Serksnyte, TU München, laura.serksnyte@tum.de

#include <cstdint>
#include <vector>
#include <string>
#include <bitset>
#include <algorithm>
#include <random>
#include <chrono>

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
  kTrackTrackTrack,
  kTrackTrackV0,
  kNTasks,
};
} // namespace CollisionMasks

struct femoDreamCollisionMasker {
  Produces<FDColMasks> Masks;
  Produces<FDDownSample> DownSample;

  // configurable for downsampling
  Configurable<float> ConfDownsampling{"ConfDownsampling", -1., "Fraction of events to be used in mixed event sample. Factor should be between 0 and 1. Deactivate with negative value"};
  Configurable<uint64_t> ConfSeed{"ConfSeed", 0, "Seed for downsampling. Set to 0 for using a seed unique in time."};

  std::mt19937* rng = nullptr;

  // particle selection bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackCutBits;
  // particle tpc pid bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackPIDTPCBits;
  // particle tpc pid bits for rejection
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackPIDTPCBitsReject;
  // particle tpctof pid bits
  std::array<std::vector<femtodreamparticle::cutContainerType>, CollisionMasks::kNParts> TrackPIDTPCTOFBits;
  // particle momemtum threshold for PID
  std::array<std::vector<float>, CollisionMasks::kNParts> TrackPIDThreshold;
  // bool if same particle species are mixed
  // bits indicate number of particles in the event
  // kPartOne -> 1 track, kPartTwo -> 2 tracks, kPartThree -> 3 tracks
  std::vector<bool> TrackSameSpecies;
  std::vector<bool> TrackDCACutPtDep;

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
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterTempFitVarMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterTempFitVarMax;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassMax;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassAntiMin;
  std::array<std::vector<float>, CollisionMasks::kNParts> FilterInvMassAntiMax;

  int TaskFinder = -1;

  void init(InitContext& context)
  {

    // seed rng for downsampling
    if (ConfDownsampling.value > 0) {
      uint64_t randomSeed = 0;
      if (ConfSeed.value == 0) {
        randomSeed = static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
      } else {
        randomSeed = ConfSeed.value;
      }
      rng = new std::mt19937(randomSeed);
    }

    int counter = 0;
    std::bitset<8 * sizeof(aod::femtodreamcollision::BitMaskType)> bitmask;
    std::vector<std::string> MatchedWorkflows;
    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";
    LOG(info) << " Collision masker self-configuration ";
    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";
    auto& workflows = context.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      if (device.name.find("femto-dream-pair-task-track-track") != std::string::npos) {
        LOG(info) << "Matched workflow: " << device.name;
        bitmask.set(counter);
        LOG(info) << "Assign bitmask: " << bitmask.to_string();
        bitmask.reset();
        counter++;
        TaskFinder = CollisionMasks::kTrackTrack;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("Track1.CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCBit_Reject")) == 0) {
            TrackPIDTPCBitsReject.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.PtMin")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.PtMax")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.EtaMin")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.EtaMax")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.TempFitVarMin")) == 0) {
            FilterTempFitVarMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.TempFitVarMax")) == 0) {
            FilterTempFitVarMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track2.TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track2.TPCBit_Reject")) == 0) {
            TrackPIDTPCBitsReject.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track2.TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track2.PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.PtMin")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.PtMax")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.EtaMin")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.EtaMax")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.TempFitVarMin")) == 0) {
            FilterTempFitVarMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track2.TempFitVarMax")) == 0) {
            FilterTempFitVarMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Option.SameSpecies")) == 0) {
            TrackSameSpecies.push_back(option.defaultValue.get<bool>());
          } else if (option.name.compare(std::string("Option.DCACutPtDep")) == 0) {
            TrackDCACutPtDep.push_back(option.defaultValue.get<bool>());
          }
        }
      } else if (device.name.find(std::string("femto-dream-pair-task-track-v0")) != std::string::npos) {
        LOG(info) << "Matched workflow: " << device.name;
        bitmask.set(counter);
        LOG(info) << "Assign bitmask: " << bitmask.to_string();
        bitmask.reset();
        counter++;
        TaskFinder = CollisionMasks::kTrackV0;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("Option.DCACutPtDep")) == 0) {
            TrackDCACutPtDep.push_back(option.defaultValue.get<bool>());
          } else if (option.name.compare(std::string("Track1.CutBit")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCBit_Reject")) == 0) {
            TrackPIDTPCBitsReject.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.TPCTOFBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Track1.PIDThres")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.PtMin")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.PtMax")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.EtaMin")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.EtaMax")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.TempFitVarMin")) == 0) {
            FilterTempFitVarMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Track1.TempFitVarMax")) == 0) {
            FilterTempFitVarMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.CutBit")) == 0) {
            V0CutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("V02.PtMin")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.PtMax")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.EtaMin")) == 0) {
            FilterEtaMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.EtaMax")) == 0) {
            FilterEtaMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.InvMassMin")) == 0) {
            FilterInvMassMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.InvMassMax")) == 0) {
            FilterInvMassMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.InvMassAntiMin")) == 0) {
            FilterInvMassAntiMin.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.InvMassAntiMax")) == 0) {
            FilterInvMassAntiMax.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("V02.ChildPos_CutBit")) == 0) {
            PosChildCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("V02.ChildPos_TPCBit")) == 0) {
            PosChildPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("V02.ChildNeg_CutBit")) == 0) {
            NegChildCutBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("V02.ChildNeg_TPCBit")) == 0) {
            NegChildPIDTPCBits.at(CollisionMasks::kPartTwo).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          }
        }
      } else if ((device.name.find("femto-dream-triplet-task-track-track-track") != std::string::npos) || (device.name.find("femto-dream-triplet-task-track-track-track-pb-pb") != std::string::npos)) {
        LOG(info) << "Matched workflow: " << device.name;
        TaskFinder = CollisionMasks::kTrackTrackTrack;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("ConfCutPart")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTPCPIDBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTPCTOFPIDBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfPIDthrMom")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMaxpT")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMinpT")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMaxDCAxy")) == 0) {
            FilterTempFitVarMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMinDCAxy")) == 0) {
            FilterTempFitVarMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfDCACutPtDep")) == 0) {
            TrackDCACutPtDep.push_back(option.defaultValue.get<bool>());
          }
        }
      } else if ((device.name.find("femto-dream-triplet-task-track-track-v0") != std::string::npos) || (device.name.find("femto-dream-triplet-task-track-track-v0-pb-pb") != std::string::npos)) {
        LOG(info) << "Matched workflow: " << device.name;
        TaskFinder = CollisionMasks::kTrackTrackV0;
        for (auto const& option : device.options) {
          if (option.name.compare(std::string("ConfCutPart")) == 0) {
            TrackCutBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTPCPIDBit")) == 0) {
            TrackPIDTPCBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfTPCTOFPIDBit")) == 0) {
            TrackPIDTPCTOFBits.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("ConfPIDthrMom")) == 0) {
            TrackPIDThreshold.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMaxpT")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMinpT")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMaxDCAxy")) == 0) {
            FilterTempFitVarMax.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfMinDCAxy")) == 0) {
            FilterTempFitVarMin.at(CollisionMasks::kPartOne).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfCutV0")) == 0) {
            V0CutBits.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Conf_ChildPos_CutV0")) == 0) {
            PosChildCutBits.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Conf_ChildPos_TPCBitV0")) == 0) {
            PosChildPIDTPCBits.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Conf_ChildNeg_CutV0")) == 0) {
            NegChildCutBits.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Conf_ChildNeg_TPCBitV0")) == 0) {
            NegChildPIDTPCBits.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<femtodreamparticle::cutContainerType>());
          } else if (option.name.compare(std::string("Conf_minInvMass_V0")) == 0) {
            FilterInvMassMin.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Conf_maxInvMass_V0")) == 0) {
            FilterInvMassMax.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Conf_minInvMassAnti_V0")) == 0) {
            FilterInvMassAntiMin.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Conf_maxInvMassAnti_V0")) == 0) {
            FilterInvMassAntiMax.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Conf_minPt_V0")) == 0) {
            FilterPtMin.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("Conf_maxPt_V0")) == 0) {
            FilterPtMax.at(CollisionMasks::kPartThree).push_back(option.defaultValue.get<float>());
          } else if (option.name.compare(std::string("ConfDCACutPtDep")) == 0) {
            TrackDCACutPtDep.push_back(option.defaultValue.get<bool>());
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

  // make bitmask for a track for two body task
  template <typename T, typename R, typename S>
  void MaskForTrack(T& BitSet, CollisionMasks::Parts P, R& track, S& counter)
  {
    if (track.partType() != static_cast<uint8_t>(femtodreamparticle::kTrack)) {
      return;
    }
    for (size_t index = 0; index < TrackCutBits.at(P).size(); index++) {
      // check filter cuts
      // if they are not passed, skip the particle
      if (track.pt() < FilterPtMin.at(P).at(index) || track.pt() > FilterPtMax.at(P).at(index) ||
          track.eta() < FilterEtaMin.at(P).at(index) || track.eta() > FilterEtaMax.at(P).at(index)) {
        continue;
      }
      // check if we apply pt dependend dca cut
      // if they do not pass this cut, skip particle as well
      if (TrackDCACutPtDep.at(index)) {
        if (std::fabs(track.tempFitVar()) > 0.0105f + (0.035f / std::pow(track.pt(), 1.1f))) {
          continue;
        }
      } else {
        // or cut on the DCA directly
        if (track.tempFitVar() < FilterTempFitVarMin.at(P).at(index) || track.tempFitVar() > FilterTempFitVarMax.at(P).at(index)) {
          continue;
        }
      }
      // set the bit at the index of the selection equal to one if the track passes all selections
      // check track cuts
      if ((track.cut() & TrackCutBits.at(P).at(index)) == TrackCutBits.at(P).at(index)) {
        // check pid cuts
        if (track.p() <= TrackPIDThreshold.at(P).at(index)) {
          if ((track.pidcut() & TrackPIDTPCBits.at(P).at(index)) == TrackPIDTPCBits.at(P).at(index) && ((track.pidcut() & TrackPIDTPCBitsReject.at(P).at(index)) == 0u)) {
            BitSet.at(P).set(index);
            counter.at(P).at(index)++;
          }
        } else {
          if ((track.pidcut() & TrackPIDTPCTOFBits.at(P).at(index)) == TrackPIDTPCTOFBits.at(P).at(index)) {
            BitSet.at(P).set(index);
            counter.at(P).at(index)++;
          }
        }
      }
    }
  }

  // Make bitmask for a track for three body task
  // This function ALWAYS checks Track P; if howManyTracksToSetInMask is set to more than 1, it also checks how many
  // tracks there are which pass Track P selection and accordingly sets the bits in the mask for Tracks P+1 or both P+1 and P+2
  template <typename T>
  void MaskForTrack_ThreeBody(T& BitSet, CollisionMasks::Parts P, FDCollision const& /*col*/, FDParticles const& parts, int howManyTracksToSetInMask)
  {
    if (howManyTracksToSetInMask == 2) {
      if (P == 2) {
        LOG(fatal) << "You are checking third particle out of three but asking to set two new bits, not possible!";
      }
    }
    if (howManyTracksToSetInMask == 3) {
      if (P >= 1) {
        LOG(fatal) << "You are checking second or third particle out of three but asking to set three new bits, not possible!";
      }
    }
    for (size_t index = 0; index < TrackCutBits.at(P).size(); index++) {
      int countTracksWhichPassSelection = 0;
      for (auto const& track : parts) {
        if (track.partType() != static_cast<uint8_t>(femtodreamparticle::kTrack)) {
          continue;
        }

        // check filter cuts
        if (track.pt() < FilterPtMin.at(P).at(index) || track.pt() > FilterPtMax.at(P).at(index)) {
          // if they are not passed, skip the particle
          continue;
        }

        if (TrackDCACutPtDep.at(index)) {
          if (std::fabs(track.tempFitVar()) > 0.0105f + (0.035f / std::pow(track.pt(), 1.1f))) {
            continue;
          }
        } else {
          // or cut on the DCA directly
          if (track.tempFitVar() < FilterTempFitVarMin.at(P).at(index) || track.tempFitVar() > FilterTempFitVarMax.at(P).at(index)) {
            continue;
          }
        }

        // set the bit at the index of the selection equal to one if the track passes all selections
        // check track cuts
        if ((track.cut() & TrackCutBits.at(P).at(index)) == TrackCutBits.at(P).at(index)) {
          // check pid cuts
          if (track.p() <= TrackPIDThreshold.at(P).at(index)) {
            if ((track.pidcut() & TrackPIDTPCBits.at(P).at(index)) == TrackPIDTPCBits.at(P).at(index)) {
              countTracksWhichPassSelection = countTracksWhichPassSelection + 1;
            }
          } else {
            if ((track.pidcut() & TrackPIDTPCTOFBits.at(P).at(index)) == TrackPIDTPCTOFBits.at(P).at(index)) {
              countTracksWhichPassSelection = countTracksWhichPassSelection + 1;
            }
          }
        }
      }
      if (countTracksWhichPassSelection >= 1)
        BitSet.at(P).set(index);
      if (howManyTracksToSetInMask == 2) {
        if (countTracksWhichPassSelection >= 2) {
          if (P == CollisionMasks::kPartOne)
            BitSet.at(CollisionMasks::kPartTwo).set(index);
          if (P == CollisionMasks::kPartTwo)
            BitSet.at(CollisionMasks::kPartThree).set(index);
        }
      }
      if (howManyTracksToSetInMask == 3) {
        if (countTracksWhichPassSelection >= 2)
          BitSet.at(CollisionMasks::kPartTwo).set(index);
        if (countTracksWhichPassSelection >= 3)
          BitSet.at(CollisionMasks::kPartThree).set(index);
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

  // Make bitmask for a V0 for three body task
  // This function ALWAYS checks V0 P; if howManyVoaToSetInMask is set to more than 1, it also checks how many
  // V0s there are which pass V0 P selection and accordingly sets the bits in the mask for V0s P+1 or both P+1 and P+2
  template <typename T>
  void MaskForV0_ThreeBody(T& BitSet, CollisionMasks::Parts P, FDCollision const& /*col*/, FDParticles const& parts, int howManyV0sToSetInMask)
  {
    if (howManyV0sToSetInMask == 2) {
      if (P == 2) {
        LOG(fatal) << "You are checking third particle out of three but asking to set two new bits, not possible!";
      }
    }
    if (howManyV0sToSetInMask == 3) {
      if (P >= 1) {
        LOG(fatal) << "You are checking second or third particle out of three but asking to set three new bits, not possible!";
      }
    }
    for (size_t index = 0; index < V0CutBits.at(P).size(); index++) {
      int countV0sWhichPassSelection = 0;
      for (auto const& V0 : parts) {
        if (V0.partType() != static_cast<uint8_t>(femtodreamparticle::kV0)) {
          continue;
        }
        // check filter cuts for pt
        if (V0.pt() < FilterPtMin.at(P).at(index) || V0.pt() > FilterPtMax.at(P).at(index)) {
          continue;
        }
        // check filter cuts for Minv
        if (V0.mLambda() < FilterInvMassMin.at(P).at(index) || V0.mLambda() > FilterInvMassMax.at(P).at(index) ||
            V0.mAntiLambda() < FilterInvMassAntiMin.at(P).at(index) || V0.mAntiLambda() > FilterInvMassAntiMax.at(P).at(index)) {
          continue;
        }
        // set the bit at the index of the selection equal to one if the track passes all selections
        // check track cuts
        if ((V0.cut() & V0CutBits.at(P).at(index)) == V0CutBits.at(P).at(index)) {
          // check daughter cuts and pid cuts
          const auto& posChild = parts.iteratorAt(V0.index() - 2);
          const auto& negChild = parts.iteratorAt(V0.index() - 1);
          if ((posChild.cut() & PosChildCutBits.at(P).at(index)) == PosChildCutBits.at(P).at(index) &&
              (posChild.pidcut() & PosChildPIDTPCBits.at(P).at(index)) == PosChildPIDTPCBits.at(P).at(index) &&
              (negChild.cut() & NegChildCutBits.at(P).at(index)) == NegChildCutBits.at(P).at(index) &&
              (negChild.pidcut() & NegChildPIDTPCBits.at(P).at(index)) == NegChildPIDTPCBits.at(P).at(index)) {
            countV0sWhichPassSelection = countV0sWhichPassSelection + 1;
          }
        }
      }
      if (countV0sWhichPassSelection >= 1)
        BitSet.at(P).set(index);

      if (howManyV0sToSetInMask == 2) {
        if (countV0sWhichPassSelection >= 2) {
          if (P == CollisionMasks::kPartOne)
            BitSet.at(CollisionMasks::kPartTwo).set(index);
          if (P == CollisionMasks::kPartTwo)
            BitSet.at(CollisionMasks::kPartThree).set(index);
        }
      }

      if (howManyV0sToSetInMask == 3) {
        if (countV0sWhichPassSelection >= 2)
          BitSet.at(CollisionMasks::kPartTwo).set(index);
        if (countV0sWhichPassSelection >= 3)
          BitSet.at(CollisionMasks::kPartThree).set(index);
      }
    }
  }

  void process(FDCollision const& col, FDParticles const& parts)
  {
    // create a bit mask for particle one, particle two and particle three
    std::array<std::bitset<8 * sizeof(femtodreamcollision::BitMaskType)>, CollisionMasks::kNParts> Mask = {{0}};

    // array to keep track of found particles
    std::array<std::array<int, 8 * sizeof(femtodreamcollision::BitMaskType)>, CollisionMasks::kNParts> FoundParts = {{{0}}};

    switch (TaskFinder) {
      case CollisionMasks::kTrackTrack:
        // pair-track-track task
        // create mask for track 1 and track 2
        for (auto const& part : parts) {
          MaskForTrack(Mask, CollisionMasks::kPartOne, part, FoundParts);
          MaskForTrack(Mask, CollisionMasks::kPartTwo, part, FoundParts);
        }
        // check if SameSpecies option is set
        // if so, only set the bit for kPartTwo to true if there are at least to tracks in the collision that pass the selections
        // this allows to (optionally) select only events with at least two particles for the mixing
        for (size_t index = 0; index < TrackSameSpecies.size(); index++) {
          if (TrackSameSpecies.at(index) == true) {
            Mask.at(CollisionMasks::kPartTwo).reset(index);
            if (FoundParts.at(CollisionMasks::kPartOne).at(index) >= 2) {
              Mask.at(CollisionMasks::kPartTwo).set(index);
            }
          }
        }
        break;
      case CollisionMasks::kTrackV0:
        // pair-track-v0 task
        // create mask for track 1 and v0 2
        for (auto const& part : parts) {
          MaskForTrack(Mask, CollisionMasks::kPartOne, part, FoundParts);
          MaskForV0(Mask, CollisionMasks::kPartTwo, part, parts);
        }
        break;
      case CollisionMasks::kTrackTrackTrack:
        // triplet-track-track-track task
        // this assumes three identical tracks are used in the task
        MaskForTrack_ThreeBody(Mask, CollisionMasks::kPartOne, col, parts, 3);
        break;
      case CollisionMasks::kTrackTrackV0:
        // triplet-track-track-v0 task
        // this assumes two identical tracks and one V0 is used in the task; V0 is last particle
        MaskForTrack_ThreeBody(Mask, CollisionMasks::kPartOne, col, parts, 2);
        MaskForV0_ThreeBody(Mask, CollisionMasks::kPartThree, col, parts, 1);
        break;
      // TODO: add all supported pair/triplet tasks
      default:
        LOG(fatal) << "No femtodream task found!";
    }
    // fill bitmask for each collision
    Masks(static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartOne).to_ulong()),
          static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartTwo).to_ulong()),
          static_cast<femtodreamcollision::BitMaskType>(Mask.at(CollisionMasks::kPartThree).to_ulong()));

    bool UseInMixedEvent = true;
    std::uniform_real_distribution<> dist(0, 1);

    if (ConfDownsampling.value > 0 && (1 - dist(*rng)) > ConfDownsampling.value) {
      UseInMixedEvent = false;
    }

    DownSample(UseInMixedEvent);
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<femoDreamCollisionMasker>(cfgc)};
  return workflow;
};
