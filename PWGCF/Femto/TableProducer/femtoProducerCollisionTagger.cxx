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

/// \file femtoProducerCollisionTagger.cxx
/// \brief Self-configuring task that tags femto collisions with a bit per (pair/triplet task
///        instance, pair type) found in the same workflow, set if the collision has what it
///        takes to form at least one pair of that type: both the particle-1-role and
///        particle-2-role candidate present (or, for a same-species pair, at least 2 distinct
///        candidates passing the one shared selection). Used to pre-filter the mixed-event pool
///        instead of discovering empty collisions only after slicing. Supports any number of
///        same-type "subwagons" (the same task run several times with different cuts in one
///        train) -- bits are assigned in discovery order.
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/Femto/Core/collisionTagger.h"
#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <bitset>
#include <string>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoProducerCollisionTagger {
  o2::framework::Produces<o2::aod::FColTags> producedTags;

  o2::framework::SliceCache cache;

  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;

  o2::framework::Preslice<FemtoTracks> perColTracks = o2::aod::femtobase::stored::fColId;
  o2::framework::Preslice<FemtoLambdas> perColLambdas = o2::aod::femtobase::stored::fColId;
  o2::framework::Preslice<FemtoK0shorts> perColK0shorts = o2::aod::femtobase::stored::fColId;

  // --- femto-pair-track-track: one bit per instance, "a track-track pair is possible here" ---
  std::vector<collisiontagger::TrackTrackVariant> mTrackTrackVariants;

  // --- femto-pair-track-v0: one bit per instance per enabled v0 species ---
  std::vector<collisiontagger::TrackV0Variant> mTrackLambdaVariants;
  std::vector<collisiontagger::TrackV0Variant> mTrackK0shortVariants;

  void init(o2::framework::InitContext& context)
  {
    // FColTags gets exactly one row written per collision; if more than one process function were
    // enabled, each would try to write its own row for the same collision. Until FColTags is
    // extended to carry one column per wagon type, only one process function may run at a time --
    // run separate tagger instances (with separate output files/merges) for separate wagon types.
    int activeProcesses = static_cast<int>(doprocessTrackTrack) + static_cast<int>(doprocessTrackV0);
    if (activeProcesses != 1) {
      LOG(fatal) << "Exactly one process function must be enabled on the collision tagger, found " << activeProcesses << ". Breaking...";
    }

    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";
    LOG(info) << " Collision tagger self-configuration  ";
    LOG(info) << "*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*";

    // discoverTaggedRoles() is the single canonical bit-assignment walk, shared with every
    // consuming pair/triplet task (see Core/collisionTagger.h) -- it assigns bits to every
    // recognized task instance found in the workflow, regardless of which of THIS tagger
    // instance's own process functions is enabled, so that bit numbers stay globally consistent
    // even when multiple tagger instances (one per wagon type) run side by side.
    auto& workflows = context.services().get<o2::framework::RunningWorkflowInfo const>();
    auto roles = collisiontagger::discoverTaggedRoles(workflows);

    auto findDevice = [&workflows](std::string const& name) -> o2::framework::DeviceSpec const& {
      for (auto const& device : workflows.devices) {
        if (device.name == name) {
          return device;
        }
      }
      LOG(fatal) << "Device " << name << " vanished between discovery and lookup. Breaking...";
      return workflows.devices.front(); // unreachable, silences -Wreturn-type
    };

    for (auto const& r : roles) {
      LOG(info) << "  bit " << r.bit << " = " << r.deviceName << " " << r.role;
      if (doprocessTrackTrack && r.role == collisiontagger::role::TrackTrackPair && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackTrack) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackTrackVariant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readTrackSelection(device, "TrackSelection2");
        v.sameSpecies = collisiontagger::hasEnabledOption(device, "Mixing.sameSpecies");
        v.bit = r.bit;
        mTrackTrackVariants.push_back(v);
      } else if (doprocessTrackV0 && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackV0) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == collisiontagger::role::TrackLambdaPair) {
          collisiontagger::TrackV0Variant v;
          v.trackSel = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.v0Sel = collisiontagger::readV0Selection(device, "LambdaSelection1");
          v.bit = r.bit;
          mTrackLambdaVariants.push_back(v);
        } else if (r.role == collisiontagger::role::TrackK0shortPair) {
          collisiontagger::TrackV0Variant v;
          v.trackSel = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.v0Sel = collisiontagger::readV0Selection(device, "K0shortSelection1");
          v.bit = r.bit;
          mTrackK0shortVariants.push_back(v);
        }
      }
      // TODO: add the remaining pair/triplet task types (track-cascade, track-kink, v0-v0, d0-d0,
      // resonances, triplets, ...), following the same pattern.
    }

    if (doprocessTrackTrack && mTrackTrackVariants.empty()) {
      LOG(fatal) << "processTrackTrack is enabled, but no femto-pair-track-track device was found in this workflow.";
    }
    if (doprocessTrackV0 && mTrackLambdaVariants.empty() && mTrackK0shortVariants.empty()) {
      LOG(fatal) << "processTrackV0 is enabled, but no femto-pair-track-v0 device (with processLambdaSameEvent/MixedEvent "
                 << "or processK0shortSameEvent/MixedEvent enabled) was found in this workflow.";
    }
    if (!roles.empty() && roles.back().bit >= collisiontagger::MaxTagBits) {
      LOG(fatal) << "Too many tag bits requested (" << (roles.back().bit + 1) << "); collision tagger supports up to " << collisiontagger::MaxTagBits << ".";
    }
  }

  void processTrackTrack(FemtoCollisions const& cols, FemtoTracks const& tracks)
  {
    for (auto const& col : cols) {
      std::bitset<collisiontagger::MaxTagBits> tag;
      auto tracksThisCol = tracks.sliceBy(perColTracks, col.globalIndex());
      for (auto const& v : mTrackTrackVariants) {
        bool found1 = false;
        int needed2 = v.sameSpecies ? 2 : 1;
        int count2 = 0;
        for (auto const& track : tracksThisCol) {
          found1 = found1 || collisiontagger::trackPasses(track, v.sel1);
          if (count2 < needed2 && collisiontagger::trackPasses(track, v.sel2)) {
            count2++;
          }
          if (found1 && count2 >= needed2) {
            break;
          }
        }
        if (found1 && count2 >= needed2) {
          tag.set(v.bit);
        }
      }
      producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
    }
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackTrack, "Tag collisions for femto-pair-track-track wagon(s)", false);

  void processTrackV0(FemtoCollisions const& cols, FemtoTracks const& tracks, FemtoLambdas const& lambdas, FemtoK0shorts const& k0shorts)
  {
    for (auto const& col : cols) {
      std::bitset<collisiontagger::MaxTagBits> tag;
      auto tracksThisCol = tracks.sliceBy(perColTracks, col.globalIndex());

      if (!mTrackLambdaVariants.empty()) {
        auto lambdasThisCol = lambdas.sliceBy(perColLambdas, col.globalIndex());
        for (auto const& v : mTrackLambdaVariants) {
          bool hasTrack = false;
          for (auto const& track : tracksThisCol) {
            if (collisiontagger::trackPasses(track, v.trackSel)) {
              hasTrack = true;
              break;
            }
          }
          if (!hasTrack) {
            continue;
          }
          for (auto const& lambda : lambdasThisCol) {
            if (collisiontagger::lambdaPasses(lambda, v.v0Sel)) {
              tag.set(v.bit);
              break;
            }
          }
        }
      }

      if (!mTrackK0shortVariants.empty()) {
        auto k0shortsThisCol = k0shorts.sliceBy(perColK0shorts, col.globalIndex());
        for (auto const& v : mTrackK0shortVariants) {
          bool hasTrack = false;
          for (auto const& track : tracksThisCol) {
            if (collisiontagger::trackPasses(track, v.trackSel)) {
              hasTrack = true;
              break;
            }
          }
          if (!hasTrack) {
            continue;
          }
          for (auto const& k0short : k0shortsThisCol) {
            if (collisiontagger::k0shortPasses(k0short, v.v0Sel)) {
              tag.set(v.bit);
              break;
            }
          }
        }
      }

      producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
    }
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackV0, "Tag collisions for femto-pair-track-v0 wagon(s)", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    o2::framework::adaptAnalysisTask<FemtoProducerCollisionTagger>(context),
  };
  return workflow;
}
