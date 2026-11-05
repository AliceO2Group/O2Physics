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

/// \file collisionTagger.h
/// \brief Shared machinery for the self-configuring collision tagger (TableProducer/femtoProducerCollisionTagger.cxx)
///        and for pair/triplet tasks that want to read back their own assigned tag bit(s).
///
/// The tagger walks the running workflow's device list once, and for every matching pair/triplet
/// task instance it finds, assigns the next free bit(s) of FColTags::collisionTag to that
/// instance's particle role(s), in discovery order. A consuming task cannot read the tagger's
/// bit assignment back out of it directly: RunningWorkflowInfo only exposes each device's
/// CLI/JSON-resolved option values, a static snapshot taken once at workflow construction -- it is
/// never updated by what another device's init() computes at runtime, and there is no ordering
/// guarantee between different devices' init() calls regardless. So instead, discoverTaggedRoles()
/// is the single canonical implementation of the bit-assignment walk, included and called
/// identically by the tagger (to assign bits and fill FColTags) and by every consuming task (to
/// recompute, from the same static device list, which bit(s) ended up belonging to its own device
/// name). Keeping this walk in one place is what keeps the two sides from silently drifting apart.
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTO_CORE_COLLISIONTAGGER_H_
#define PWGCF_FEMTO_CORE_COLLISIONTAGGER_H_

#include "PWGCF/Femto/Core/dataTypes.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/RunningWorkflowInfo.h>

#include <cmath>
#include <string>
#include <vector>

namespace o2::analysis::femto::collisiontagger
{

constexpr int MaxTagBits = 8 * sizeof(datatypes::CollisionTagType);

// Supporting an arbitrary number of same-type "subwagons" (same task, different cuts, all running
// in one train) rules out reusing Partition<T> for the cut check itself: it is a compile-time
// construct, permanently bound to one specific Configurable object's fields. So instead, each
// matching sibling device's configured cut VALUES are copied into a plain value struct (no
// Configurable, no Partition), one push_back per discovered subwagon, and checked by hand against
// each particle -- the same approach femtoDreamCollisionMasker.cxx uses for its own multi-variant
// support. The check functions below are a direct, faithful port of the corresponding
// MAKE_*_PARTITION macro in Core/partitions.h, just expressed as a runtime check instead of a
// compiled Gandiva expression.
struct TrackSelectionValues {
  int chargeAbs = 1;
  int chargeSign = 1;
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  datatypes::TrackMaskType maskLowMomentum = 0;
  datatypes::TrackMaskType maskHighMomentum = 0;
  datatypes::TrackMaskType rejectionMaskLowMomentum = 0;
  datatypes::TrackMaskType rejectionMaskHighMomentum = 0;
  float pidThres = 1e6f;
};

struct V0SelectionValues {
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  float massMin = 0.f;
  float massMax = 1e6f;
  datatypes::V0MaskType mask = 0;
  int sign = 0; // only used for lambdas; k0shorts carry a dummy value for compatibility
};

// one bit = one pair type for one task instance: both sides' selections, read together, since the
// bit means "both are present", not "particle 1 is present" or "particle 2 is present" separately
struct TrackTrackVariant {
  TrackSelectionValues sel1;
  TrackSelectionValues sel2;
  bool sameSpecies;
  int bit;
};

struct TrackV0Variant {
  TrackSelectionValues trackSel;
  V0SelectionValues v0Sel;
  int bit;
};

inline TrackSelectionValues readTrackSelection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  TrackSelectionValues v;
  for (auto const& option : device.options) {
    if (option.name == prefix + ".chargeAbs") {
      v.chargeAbs = option.defaultValue.get<int>();
    } else if (option.name == prefix + ".chargeSign") {
      v.chargeSign = option.defaultValue.get<int>();
    } else if (option.name == prefix + ".ptMin") {
      v.ptMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".ptMax") {
      v.ptMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".etaMin") {
      v.etaMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".etaMax") {
      v.etaMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".phiMin") {
      v.phiMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".phiMax") {
      v.phiMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".maskLowMomentum") {
      v.maskLowMomentum = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".maskHighMomentum") {
      v.maskHighMomentum = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".rejectionMaskLowMomentum") {
      v.rejectionMaskLowMomentum = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".rejectionMaskHighMomentum") {
      v.rejectionMaskHighMomentum = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".pidThres") {
      v.pidThres = option.defaultValue.get<float>();
    }
  }
  return v;
}

inline V0SelectionValues readV0Selection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  V0SelectionValues v;
  for (auto const& option : device.options) {
    if (option.name == prefix + ".ptMin") {
      v.ptMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".ptMax") {
      v.ptMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".etaMin") {
      v.etaMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".etaMax") {
      v.etaMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".phiMin") {
      v.phiMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".phiMax") {
      v.phiMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".massMin") {
      v.massMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".massMax") {
      v.massMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".mask") {
      v.mask = option.defaultValue.get<datatypes::V0MaskType>();
    } else if (option.name == prefix + ".sign") {
      v.sign = option.defaultValue.get<int>();
    }
  }
  return v;
}

// reads the current value of a plain (non-group-prefixed) boolean option on a sibling device,
// e.g. one of its doprocessXXX-backed PROCESS_SWITCH flags
inline bool hasEnabledOption(o2::framework::DeviceSpec const& device, std::string const& name)
{
  for (auto const& option : device.options) {
    if (option.name == name) {
      return option.defaultValue.get<bool>();
    }
  }
  return false;
}

// mirrors MAKE_TRACK_PARTITION(selection) in Core/partitions.h
template <typename T>
bool trackPasses(T const& track, TrackSelectionValues const& sel)
{
  if (sel.chargeSign != 0 && ((sel.chargeSign > 0) != (track.signedPt() > 0.f))) {
    return false;
  }
  float absPt = std::fabs(sel.chargeAbs * track.signedPt());
  if (!(absPt > sel.ptMin && absPt < sel.ptMax)) {
    return false;
  }
  if (!(track.eta() > sel.etaMin && track.eta() < sel.etaMax)) {
    return false;
  }
  if (!(track.phi() > sel.phiMin && track.phi() < sel.phiMax)) {
    return false;
  }
  float momentum = absPt * std::cosh(track.eta());
  if (momentum <= sel.pidThres) {
    return (track.mask() & sel.maskLowMomentum) == sel.maskLowMomentum &&
           (track.mask() & sel.rejectionMaskLowMomentum) == static_cast<datatypes::TrackMaskType>(0);
  }
  return (track.mask() & sel.maskHighMomentum) == sel.maskHighMomentum &&
         (track.mask() & sel.rejectionMaskHighMomentum) == static_cast<datatypes::TrackMaskType>(0);
}

// mirrors MAKE_LAMBDA_PARTITION(selection) in Core/partitions.h
template <typename T>
bool lambdaPasses(T const& v0, V0SelectionValues const& sel)
{
  if (sel.sign != 0 && ((sel.sign > 0) != (v0.signedPt() > 0.f))) {
    return false;
  }
  float absPt = std::fabs(v0.signedPt());
  if (!(absPt > sel.ptMin && absPt < sel.ptMax)) {
    return false;
  }
  if (!(v0.eta() > sel.etaMin && v0.eta() < sel.etaMax)) {
    return false;
  }
  if (!(v0.phi() > sel.phiMin && v0.phi() < sel.phiMax)) {
    return false;
  }
  if (!(v0.mass() > sel.massMin && v0.mass() < sel.massMax)) {
    return false;
  }
  return (v0.mask() & sel.mask) == sel.mask;
}

// mirrors MAKE_K0SHORT_PARTITION(selection) in Core/partitions.h
template <typename T>
bool k0shortPasses(T const& v0, V0SelectionValues const& sel)
{
  if (!(v0.pt() > sel.ptMin && v0.pt() < sel.ptMax)) {
    return false;
  }
  if (!(v0.eta() > sel.etaMin && v0.eta() < sel.etaMax)) {
    return false;
  }
  if (!(v0.phi() > sel.phiMin && v0.phi() < sel.phiMax)) {
    return false;
  }
  if (!(v0.mass() > sel.massMin && v0.mass() < sel.massMax)) {
    return false;
  }
  return (v0.mask() & sel.mask) == sel.mask;
}

// Role names and device-name substrings, declared once: a consuming task needs the exact same
// role spelling discoverTaggedRoles() below uses, to look its own bit back up via findOwnBit() --
// fed in as ConfTagging's template argument below, so a task's role is fixed at compile time
// rather than a loose string repeated at every call site. Declared as named `char` arrays, not
// `const char*` pointing at string literals: a pointer non-type template argument must be the
// address of a named object with static storage duration, and the standard specifically
// disallows forming that address from a string literal directly (two identical literals aren't
// guaranteed the same address). A named array has none of that ambiguity, and decays to
// `const char*` on its own wherever a plain pointer is expected (ConfTagging<role::TrackTrackPair>
// below, the std::string construction in discoverTaggedRoles(), etc.). `inline` (not just
// `constexpr`) matters too: without it, a namespace-scope variable has internal linkage, so each
// translation unit would get its own copy at its own address -- fine for a value, but the same
// template instantiation must denote the same entity everywhere it's used. `inline` guarantees
// one shared address. The device-name substrings are, unavoidably, the literal
// o2physics_add_dpl_workflow() names from the respective CMakeLists.txt -- there is no way to
// derive those at compile time, so they are hardcoded here, exactly once, rather than repeated at
// every call site.
namespace role
{
inline constexpr char TrackTrackPair[] = "TrackTrackPair";
inline constexpr char TrackLambdaPair[] = "TrackLambdaPair";
inline constexpr char TrackK0shortPair[] = "TrackK0shortPair";
} // namespace role

namespace devicename
{
constexpr const char* FemtoPairTrackTrack = "femto-pair-track-track";
constexpr const char* FemtoPairTrackV0 = "femto-pair-track-v0";
} // namespace devicename

// One assigned bit of FColTags::collisionTag: `deviceName` is the matched sibling device's own
// (unique) DPL device name, `role` names which particle-selection group of that device the bit
// covers (e.g. "TrackSelection1", "LambdaSelection1"), and `bit` is its position in collisionTag.
// One assigned bit means "this collision has what it takes to form a pair of this TYPE": both the
// particle-1-role candidate and the particle-2-role candidate (or, for a same-species pair, at
// least 2 distinct candidates passing the one shared selection) are present. One bit per pair type
// per task instance -- not one bit per particle role -- since during mixing the same bit is
// checked identically on both the anchor and the partner collision (either one could end up
// supplying either role across the full mixing sweep).
struct TaggedRole {
  std::string deviceName;
  std::string role;
  int bit;
};

// The single canonical bit-assignment walk: scans the workflow's device list in order and, for
// every recognized pair/triplet task instance found, consumes the next free bit(s) for its
// role(s) -- exactly what femtoProducerCollisionTagger.cxx uses to decide what to tag, and what a
// consuming task calls to recompute (from the same static device list) which bit(s) ended up
// being its own. Keep this in sync with femtoProducerCollisionTagger.cxx's supported task types: whichever
// task types are added there must be added here identically, in the same order, or the tagger and
// its consumers will silently disagree about which bit means what.
inline std::vector<TaggedRole> discoverTaggedRoles(o2::framework::RunningWorkflowInfo const& workflows)
{
  std::vector<TaggedRole> roles;
  int counter = 0;
  for (o2::framework::DeviceSpec const& device : workflows.devices) {
    if (device.name.find(devicename::FemtoPairTrackTrack) != std::string::npos) {
      // Gated on the one process function that actually consumes this bit: an instance that
      // never enables processMixedEventWithTagging will never look this bit up, so assigning it
      // one anyway would just waste a slot out of the hard-capped budget (MaxTagBits below) for
      // nothing.
      if (hasEnabledOption(device, "processMixedEventWithTagging")) {
        roles.push_back({device.name, role::TrackTrackPair, counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackV0) != std::string::npos) {
      // Each species gates independently on its own *MixedEventWithTagging switch -- Lambda and
      // K0short are activated completely independently of each other (e.g. same-event-only for
      // one species, tagged-mixed-event for the other is a legitimate combination), and, same as
      // TrackTrackPair above, same-event/Mc/WithFitVars variants are irrelevant here since none
      // of them ever consume a tag bit.
      if (hasEnabledOption(device, "processLambdaMixedEventWithTagging")) {
        roles.push_back({device.name, role::TrackLambdaPair, counter++});
      }
      if (hasEnabledOption(device, "processK0shortMixedEventWithTagging")) {
        roles.push_back({device.name, role::TrackK0shortPair, counter++});
      }
    }
    // TODO: add the remaining pair/triplet task types (track-cascade, track-kink, v0-v0, d0-d0,
    // resonances, triplets, ...), following the same pattern. Whatever is added here must also be
    // added to femtoProducerCollisionTagger.cxx's init(), in the same order, since that is where the
    // selection VALUES behind each of these bits actually get read and checked.
  }
  return roles;
}

// Convenience for a consuming task: find the bit this device's own `role` was assigned, or -1 if
// it was never reached (e.g. this device's own process switch wasn't enabled on the tagger, so it
// never requested this role at all).
inline int findOwnBit(std::vector<TaggedRole> const& roles, std::string const& ownDeviceName, std::string const& role)
{
  for (auto const& r : roles) {
    if (r.deviceName == ownDeviceName && r.role == role) {
      return r.bit;
    }
  }
  return -1;
}

// Carries a consuming task's own, fixed pair-type role (Role, one of the role:: constants above,
// supplied at the instantiation site -- e.g. ConfTagging<role::TrackTrackPair> confTagging;) and
// the bit that role was assigned, discovered once in init() into `mask` (not read from the JSON
// config). Bundling both together means a call like findOwnBit(roles, ownDevice.name,
// confTagging.role) always uses the exact role this same object's mask gets set for, instead of
// the role being a second, independently-spelled literal at the call site. Every task that opts
// into tagger-filtered mixing instantiates one of these; wrapped in a ConfigurableGroup (rather
// than a bare Configurable) so more fields can be added here later without having to touch every
// consumer.
template <char const* role>
struct ConfTagging : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Tagging");
  static constexpr char const* Role = role;
  o2::framework::Configurable<datatypes::CollisionTagType> internalTagMask{"internalTagMask", 0, "internal: discovered collision tagger bitmask, set automatically in init() -- do not set manually"};
};

} // namespace o2::analysis::femto::collisiontagger

#endif // PWGCF_FEMTO_CORE_COLLISIONTAGGER_H_
