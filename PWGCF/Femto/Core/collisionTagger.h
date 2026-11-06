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

#include <bitset>
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

// Despite the name, this shape (sign + pt + eta + phi + mass window + flat mask) is NOT
// V0-specific: MAKE_CASCADE_PARTITION, MAKE_SIGMA_PARTITION, MAKE_SIGMAPLUS_PARTITION,
// MAKE_D0_PARTITION, and MAKE_CHARM3PRONG_PARTITION in Core/partitions.h are all byte-for-byte
// identical to MAKE_LAMBDA_PARTITION, and CascadeMaskType/KinkMaskType/CharmHadronMaskType are
// all, like V0MaskType, uint32_t -- so this one struct (and lambdaPasses()/readV0Selection()
// below) is reused verbatim for Xi/Omega (cascade), Sigma/SigmaPlus (kink), and D0/Lc
// (charm hadron) selections too, not just Lambda. Only SigmaPlusPcm needs its own struct
// (KinkPcmSelectionValues below), for the one extra mlScoreMin field it has beyond this shape.
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
  bool sameSpecies = false;
  int bit = 0;
};

// sel1/sel2 (not trackSel/v0Sel) for the same reason TrackTrackVariant above uses sel1/sel2: so
// tagCrossVariants() below can be written once, generically, for any cross-type pair -- naming
// the fields after what they ARE (track selection, v0 selection) would tie the generic helper's
// code to this one particle-type combination.
struct TrackV0Variant {
  TrackSelectionValues sel1;
  V0SelectionValues sel2;
  int bit = 0;
};

// SigmaPlusPcm's selection is V0SelectionValues plus one extra field (mlScoreMin), so it needs
// its own value struct, check function (kinkPcmPasses, below readV0Selection), and variant type
// -- everything else about it (sign/pt/eta/phi/mass/mask) is identical to V0SelectionValues.
struct KinkPcmSelectionValues {
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  float massMin = 0.f;
  float massMax = 1e6f;
  datatypes::KinkMaskType mask = 0;
  int sign = 0;
  float mlScoreMin = -1.f;
};

struct TrackKinkPcmVariant {
  TrackSelectionValues sel1;
  KinkPcmSelectionValues sel2;
  int bit = 0;
};

// Same-type pair where BOTH legs use the V0SelectionValues shape: covers v0-v0 (lambda-lambda,
// k0short-k0short, lambda-k0short) and d0-d0 (same shape, see V0SelectionValues above), checked
// with tagPairVariants() when sel1/sel2 are the same species (sameSpecies-aware) or tagCrossVariants()
// when they are not (e.g. lambda-k0short: sameSpecies is meaningless there and left false/unused).
struct V0V0Variant {
  V0SelectionValues sel1;
  V0SelectionValues sel2;
  bool sameSpecies = false;
  int bit = 0;
};

// Phi -> K+K-, Rho0 -> pi+pi-: both daughters are the same species, so MAKE_RESONANCE_0_PARTITON
// (Core/partitions.h) only needs one pair of momentum-threshold-gated daughter masks, applied to
// whichever of the two stored columns (maskPosDau/maskNegDau) belongs to the daughter that was
// above/below threshold at build time (posDauHasHighMomentum/negDauHasHighMomentum, read straight
// off the candidate -- these switches are baked in at production, not re-derived here).
struct TwoTrackResonanceSelectionValues {
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  float massMin = 0.f;
  float massMax = 1e6f;
  datatypes::TrackMaskType dauMaskBelowThres = 0;
  datatypes::TrackMaskType dauMaskAboveThres = 0;
};

// Kstar0 -> K+pi-, Kstar0bar -> K-pi+: the physical pos/neg daughter's species flips with the
// candidate's charge sign (see MAKE_RESONANCE_1_PARTITON in Core/partitions.h), so the daughter
// masks are named by species (kaon/pion) rather than by charge, and sign picks which hypothesis a
// given candidate's signedPt must match.
struct Kstar0SelectionValues {
  int sign = 0;
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  float massMin = 0.f;
  float massMax = 1e6f;
  datatypes::TrackMaskType kaonDauMaskBelowThres = 0;
  datatypes::TrackMaskType kaonDauMaskAboveThres = 0;
  datatypes::TrackMaskType pionDauMaskBelowThres = 0;
  datatypes::TrackMaskType pionDauMaskAboveThres = 0;
};

// Pi0/Eta -> gamma gamma: photon daughters have no momentum-dependent PID switch and no pos/neg
// distinction, so (unlike TwoTrackResonanceSelectionValues above) there is only one flat mask,
// already the AND of both legs' own masks at production time (see femtotwophotonresonances::MaskDau
// in FemtoTables.h) -- same shape as V0SelectionValues/k0shortPasses, just plain pt (no sign) and a
// PhotonMaskType mask instead of V0MaskType.
struct TwoPhotonResonanceSelectionValues {
  float ptMin = 0.f;
  float ptMax = 1e6f;
  float etaMin = -1e6f;
  float etaMax = 1e6f;
  float phiMin = 0.f;
  float phiMax = 1.f * o2::constants::math::TwoPI;
  float massMin = 0.f;
  float massMax = 1e6f;
  datatypes::PhotonMaskType mask = 0;
};

// Cross-type pair variants for the three resonance-pairing tasks: track-twoTrackResonance (Phi/Rho0
// share one shape, Kstar0 needs its own because of the kaon/pion mask swap), track-twoPhotonResonance
// (Pi0/Eta), and v0-twoTrackResonance (lambda/k0short x phi/rho0/kstar0 -- same split as track-side).
struct TrackTwoTrackResonanceVariant {
  TrackSelectionValues sel1;
  TwoTrackResonanceSelectionValues sel2;
  int bit = 0;
};

struct TrackKstar0Variant {
  TrackSelectionValues sel1;
  Kstar0SelectionValues sel2;
  int bit = 0;
};

struct TrackTwoPhotonResonanceVariant {
  TrackSelectionValues sel1;
  TwoPhotonResonanceSelectionValues sel2;
  int bit = 0;
};

struct V0TwoTrackResonanceVariant {
  V0SelectionValues sel1;
  TwoTrackResonanceSelectionValues sel2;
  int bit = 0;
};

struct V0Kstar0Variant {
  V0SelectionValues sel1;
  Kstar0SelectionValues sel2;
  int bit = 0;
};

// Triplet variants, mirroring the "123 same, 12 same, or all different" pattern
// TripletTrackTrackTrackBuilder/TripletTrackTrackV0Builder/TripletTrackTrackCascadeBuilder use
// (Core/tripletBuilder.h, confMixing.particle123AreSameSpecies / particle12AreSameSpecies -- only
// these two combinations are supported, never e.g. "1&3 same"). Same-type (track-track-track):
// all three roles share one table, so same123 needs 3 distinct entries passing sel1 (sel2/sel3
// ignored, mirrored onto sel1 at discovery time -- same trick as TrackTrackVariant's sameSpecies
// above), same12 needs 2 distinct passing sel1 plus 1 more passing sel3 (sel2 ignored).
struct TrackTrackTrackVariant {
  TrackSelectionValues sel1;
  TrackSelectionValues sel2;
  TrackSelectionValues sel3;
  bool same123 = false;
  bool same12 = false;
  int bit = 0;
};

// Cross-type (track-track-v0, track-track-cascade): roles 1&2 share the track table, role 3 is a
// different table (v0 or cascade, hence V0SelectionValues shape) -- only same12 is possible here
// (no same123 since role 3 is never the same type as roles 1&2).
struct TrackTrackV0Variant {
  TrackSelectionValues sel1;
  TrackSelectionValues sel2;
  V0SelectionValues sel3;
  bool same12 = false;
  int bit = 0;
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

inline KinkPcmSelectionValues readKinkPcmSelection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  KinkPcmSelectionValues v;
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
      v.mask = option.defaultValue.get<datatypes::KinkMaskType>();
    } else if (option.name == prefix + ".sign") {
      v.sign = option.defaultValue.get<int>();
    } else if (option.name == prefix + ".mlScoreMin") {
      v.mlScoreMin = option.defaultValue.get<float>();
    }
  }
  return v;
}

inline TwoTrackResonanceSelectionValues readTwoTrackResonanceSelection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  TwoTrackResonanceSelectionValues v;
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
    } else if (option.name == prefix + ".dauMaskBelowThres") {
      v.dauMaskBelowThres = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".dauMaskAboveThres") {
      v.dauMaskAboveThres = option.defaultValue.get<datatypes::TrackMaskType>();
    }
  }
  return v;
}

inline Kstar0SelectionValues readKstar0Selection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  Kstar0SelectionValues v;
  for (auto const& option : device.options) {
    if (option.name == prefix + ".sign") {
      v.sign = option.defaultValue.get<int>();
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
    } else if (option.name == prefix + ".massMin") {
      v.massMin = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".massMax") {
      v.massMax = option.defaultValue.get<float>();
    } else if (option.name == prefix + ".kaonDauMaskBelowThres") {
      v.kaonDauMaskBelowThres = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".kaonDauMaskAboveThres") {
      v.kaonDauMaskAboveThres = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".pionDauMaskBelowThres") {
      v.pionDauMaskBelowThres = option.defaultValue.get<datatypes::TrackMaskType>();
    } else if (option.name == prefix + ".pionDauMaskAboveThres") {
      v.pionDauMaskAboveThres = option.defaultValue.get<datatypes::TrackMaskType>();
    }
  }
  return v;
}

inline TwoPhotonResonanceSelectionValues readTwoPhotonResonanceSelection(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  TwoPhotonResonanceSelectionValues v;
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
    } else if (option.name == prefix + ".dauMask") {
      v.mask = option.defaultValue.get<datatypes::PhotonMaskType>();
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

// true if ANY enabled option on this device has a name starting with `prefix` -- used instead of
// hasEnabledOption() where a pair/triplet task has several *WithTagging* process functions that
// all consume the same tag bit (processMixedEventWithTagging, processMixedEventWithTaggingAndMass,
// ...AndEventShape, ...AndFitVars, ...): matching the shared prefix means enabling any one of
// those variants correctly triggers bit assignment, instead of only the base one happening to
// have the exact name this check was written against.
inline bool hasAnyEnabledOptionStartingWith(o2::framework::DeviceSpec const& device, std::string const& prefix)
{
  for (auto const& option : device.options) {
    if (option.name.starts_with(prefix) && option.defaultValue.get<bool>()) {
      return true;
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
  if (absPt <= sel.ptMin || absPt >= sel.ptMax) {
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

// mirrors MAKE_LAMBDA_PARTITION(selection) in Core/partitions.h -- also mirrors
// MAKE_CASCADE_PARTITION/MAKE_SIGMA_PARTITION/MAKE_SIGMAPLUS_PARTITION/MAKE_D0_PARTITION/
// MAKE_CHARM3PRONG_PARTITION, which are byte-for-byte identical, so this same function is reused
// as-is for Xi/Omega, Sigma/SigmaPlus, and D0/Lc candidates too (see the comment on
// V0SelectionValues above).
template <typename T>
bool lambdaPasses(T const& v0, V0SelectionValues const& sel)
{
  if (sel.sign != 0 && ((sel.sign > 0) != (v0.signedPt() > 0.f))) {
    return false;
  }
  float absPt = std::fabs(v0.signedPt());
  if (absPt <= sel.ptMin || absPt >= sel.ptMax) {
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

// mirrors MAKE_SIGMAPLUSPCM_PARTITION(selection) in Core/partitions.h -- lambdaPasses() above
// plus the one extra mlScoreMin check.
template <typename T>
bool kinkPcmPasses(T const& kink, KinkPcmSelectionValues const& sel)
{
  if (sel.sign != 0 && ((sel.sign > 0) != (kink.signedPt() > 0.f))) {
    return false;
  }
  float absPt = std::fabs(kink.signedPt());
  if (absPt <= sel.ptMin || absPt >= sel.ptMax) {
    return false;
  }
  if (!(kink.eta() > sel.etaMin && kink.eta() < sel.etaMax)) {
    return false;
  }
  if (!(kink.phi() > sel.phiMin && kink.phi() < sel.phiMax)) {
    return false;
  }
  if (!(kink.mass() > sel.massMin && kink.mass() < sel.massMax)) {
    return false;
  }
  if (!(kink.mlScore() >= sel.mlScoreMin)) {
    return false;
  }
  return (kink.mask() & sel.mask) == sel.mask;
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

// mirrors MAKE_RESONANCE_0_PARTITON(selection) in Core/partitions.h (Phi, Rho0)
template <typename T>
bool twoTrackResonancePasses(T const& res, TwoTrackResonanceSelectionValues const& sel)
{
  if (!(res.pt() > sel.ptMin && res.pt() < sel.ptMax)) {
    return false;
  }
  if (!(res.eta() > sel.etaMin && res.eta() < sel.etaMax)) {
    return false;
  }
  if (!(res.phi() > sel.phiMin && res.phi() < sel.phiMax)) {
    return false;
  }
  if (!(res.mass() > sel.massMin && res.mass() < sel.massMax)) {
    return false;
  }
  datatypes::TrackMaskType posMask = res.posDauHasHighMomentum() ? sel.dauMaskAboveThres : sel.dauMaskBelowThres;
  if ((res.maskPosDau() & posMask) != posMask) {
    return false;
  }
  datatypes::TrackMaskType negMask = res.negDauHasHighMomentum() ? sel.dauMaskAboveThres : sel.dauMaskBelowThres;
  return (res.maskNegDau() & negMask) == negMask;
}

// mirrors MAKE_RESONANCE_1_PARTITON(selection) in Core/partitions.h (Kstar0/Kstar0bar) -- pos/neg
// daughter <-> kaon/pion resolved from the candidate's own signedPt, same as the real builder.
template <typename T>
bool kstar0Passes(T const& kstar0, Kstar0SelectionValues const& sel)
{
  if (sel.sign != 0 && ((sel.sign > 0) != (kstar0.signedPt() > 0.f))) {
    return false;
  }
  float absPt = std::fabs(kstar0.signedPt());
  if (absPt <= sel.ptMin || absPt >= sel.ptMax) {
    return false;
  }
  if (!(kstar0.eta() > sel.etaMin && kstar0.eta() < sel.etaMax)) {
    return false;
  }
  if (!(kstar0.phi() > sel.phiMin && kstar0.phi() < sel.phiMax)) {
    return false;
  }
  if (!(kstar0.mass() > sel.massMin && kstar0.mass() < sel.massMax)) {
    return false;
  }
  bool posIsKaon = kstar0.signedPt() > 0.f;
  datatypes::TrackMaskType posMaskThres = posIsKaon
                                            ? (kstar0.posDauHasHighMomentum() ? sel.kaonDauMaskAboveThres : sel.kaonDauMaskBelowThres)
                                            : (kstar0.posDauHasHighMomentum() ? sel.pionDauMaskAboveThres : sel.pionDauMaskBelowThres);
  if ((kstar0.maskPosDau() & posMaskThres) != posMaskThres) {
    return false;
  }
  datatypes::TrackMaskType negMaskThres = posIsKaon
                                            ? (kstar0.negDauHasHighMomentum() ? sel.pionDauMaskAboveThres : sel.pionDauMaskBelowThres)
                                            : (kstar0.negDauHasHighMomentum() ? sel.kaonDauMaskAboveThres : sel.kaonDauMaskBelowThres);
  return (kstar0.maskNegDau() & negMaskThres) == negMaskThres;
}

// mirrors MAKE_TWOPHOTONRESONANCE_PARTITION(selection) in Core/partitions.h (Pi0, Eta)
template <typename T>
bool twoPhotonResonancePasses(T const& res, TwoPhotonResonanceSelectionValues const& sel)
{
  if (!(res.pt() > sel.ptMin && res.pt() < sel.ptMax)) {
    return false;
  }
  if (!(res.eta() > sel.etaMin && res.eta() < sel.etaMax)) {
    return false;
  }
  if (!(res.phi() > sel.phiMin && res.phi() < sel.phiMax)) {
    return false;
  }
  if (!(res.mass() > sel.massMin && res.mass() < sel.massMax)) {
    return false;
  }
  return (res.maskDau() & sel.mask) == sel.mask;
}

using CollisionTagBitset = std::bitset<MaxTagBits>;

// true if at least `needed` entries in `range` satisfy `passes(entry, sel)`. The one primitive
// every per-collision tagging check below is built from; a new particle type only ever needs a
// new `passes` function (trackPasses/lambdaPasses/k0shortPasses above are examples), never a new
// copy of this loop.
template <typename Range, typename Selection, typename PassesFn>
inline bool hasAtLeast(Range const& range, Selection const& sel, PassesFn const& passes, int needed = 1)
{
  int count = 0;
  for (auto const& entry : range) {
    if (passes(entry, sel) && ++count >= needed) {
      return true;
    }
  }
  return false;
}

// Sets tag.set(v.bit) for every variant where `range` has what it takes to form a SAME-TYPE pair
// (both v.sel1 and v.sel2 satisfied -- or, for a same-species pair, 2 distinct entries passing
// the one shared selection v.sel2). Reusable for any future same-type pair (V0V0, D0D0, ...) --
// just a different `variants`/`range`/`passes` combination, the counting logic doesn't change.
template <typename Variant, typename Range, typename PassesFn>
inline void tagPairVariants(CollisionTagBitset& tag, std::vector<Variant> const& variants, Range const& range, PassesFn const& passes)
{
  for (auto const& v : variants) {
    bool ok1 = hasAtLeast(range, v.sel1, passes, 1);
    bool ok2 = hasAtLeast(range, v.sel2, passes, v.sameSpecies ? 2 : 1);
    if (ok1 && ok2) {
      tag.set(v.bit);
    }
  }
}

// Sets tag.set(v.bit) for every variant where range1 has what it takes for role 1 (v.sel1,
// checked with passes1) AND range2 has what it takes for role 2 (v.sel2, checked with passes2).
// Reusable for any future cross-type pair (track-cascade, track-kink, track-D0, ...) -- just a
// different `variants`/ranges/passes combination.
template <typename Variant, typename Range1, typename Range2, typename PassesFn1, typename PassesFn2>
inline void tagCrossVariants(CollisionTagBitset& tag, std::vector<Variant> const& variants, Range1 const& range1, PassesFn1 const& passes1, Range2 const& range2, PassesFn2 const& passes2)
{
  for (auto const& v : variants) {
    if (hasAtLeast(range1, v.sel1, passes1, 1) && hasAtLeast(range2, v.sel2, passes2, 1)) {
      tag.set(v.bit);
    }
  }
}

// One-call convenience for the common case of a single tagPairVariants() check immediately
// followed by producing the tag row -- covers any task type whose process function only ever
// needs one same-type-pair check per collision (e.g. processTrackTrack).
template <typename Produces, typename Variant, typename Range, typename PassesFn>
inline void tagAndProduce(Produces& producedTags, std::vector<Variant> const& variants, Range const& range, PassesFn const& passes)
{
  CollisionTagBitset tag;
  tagPairVariants(tag, variants, range, passes);
  producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
}

// Sets tag.set(v.bit) for every variant where `range` has what it takes to form a SAME-TYPE
// triplet (track-track-track): same123 needs 3 distinct entries passing sel1; same12 needs 2
// distinct passing sel1 plus 1 more passing sel3; otherwise sel1/sel2/sel3 each need their own
// independent match. Mirrors TripletTrackTrackTrackBuilder exactly (see TrackTrackTrackVariant).
template <typename Variant, typename Range, typename PassesFn>
inline void tagTripletVariants(CollisionTagBitset& tag, std::vector<Variant> const& variants, Range const& range, PassesFn const& passes)
{
  for (auto const& v : variants) {
    bool ok = false;
    if (v.same123) {
      ok = hasAtLeast(range, v.sel1, passes, 3);
    } else if (v.same12) {
      ok = hasAtLeast(range, v.sel1, passes, 2) && hasAtLeast(range, v.sel3, passes, 1);
    } else {
      ok = hasAtLeast(range, v.sel1, passes, 1) && hasAtLeast(range, v.sel2, passes, 1) && hasAtLeast(range, v.sel3, passes, 1);
    }
    if (ok) {
      tag.set(v.bit);
    }
  }
}

// Sets tag.set(v.bit) for every variant where range1 (roles 1&2, one table) and range2 (role 3, a
// different table) have what it takes to form a CROSS-TYPE triplet (track-track-v0,
// track-track-cascade): same12 needs 2 distinct entries in range1 passing sel1; otherwise role 1
// and role 2 each need their own independent match in range1. Mirrors
// TripletTrackTrackV0Builder/TripletTrackTrackCascadeBuilder exactly (see TrackTrackV0Variant).
template <typename Variant, typename Range1, typename Range2, typename PassesFn1, typename PassesFn2>
inline void tagCrossTripletVariants(CollisionTagBitset& tag, std::vector<Variant> const& variants, Range1 const& range1, PassesFn1 const& passes1, Range2 const& range2, PassesFn2 const& passes2)
{
  for (auto const& v : variants) {
    bool ok = v.same12 ? hasAtLeast(range1, v.sel1, passes1, 2)
                       : (hasAtLeast(range1, v.sel1, passes1, 1) && hasAtLeast(range1, v.sel2, passes1, 1));
    if (ok && hasAtLeast(range2, v.sel3, passes2, 1)) {
      tag.set(v.bit);
    }
  }
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
// below, the explicit static_cast<char const*>(...) in discoverTaggedRoles(), etc.). `inline` (not just
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
inline constexpr char TrackXiPair[] = "TrackXiPair";
inline constexpr char TrackOmegaPair[] = "TrackOmegaPair";
inline constexpr char TrackSigmaPair[] = "TrackSigmaPair";
inline constexpr char TrackSigmaPlusPair[] = "TrackSigmaPlusPair";
inline constexpr char TrackSigmaPlusPcmPair[] = "TrackSigmaPlusPcmPair";
inline constexpr char TrackD0Pair[] = "TrackD0Pair";
inline constexpr char TrackLcPair[] = "TrackLcPair";
inline constexpr char LambdaLambdaPair[] = "LambdaLambdaPair";
inline constexpr char K0shortK0shortPair[] = "K0shortK0shortPair";
inline constexpr char LambdaK0shortPair[] = "LambdaK0shortPair";
inline constexpr char D0D0Pair[] = "D0D0Pair";
inline constexpr char TrackPhiPair[] = "TrackPhiPair";
inline constexpr char TrackKstar0Pair[] = "TrackKstar0Pair";
inline constexpr char TrackRho0Pair[] = "TrackRho0Pair";
inline constexpr char TrackPi0Pair[] = "TrackPi0Pair";
inline constexpr char TrackEtaPair[] = "TrackEtaPair";
inline constexpr char LambdaPhiPair[] = "LambdaPhiPair";
inline constexpr char LambdaKstar0Pair[] = "LambdaKstar0Pair";
inline constexpr char LambdaRho0Pair[] = "LambdaRho0Pair";
inline constexpr char K0shortPhiPair[] = "K0shortPhiPair";
inline constexpr char K0shortKstar0Pair[] = "K0shortKstar0Pair";
inline constexpr char K0shortRho0Pair[] = "K0shortRho0Pair";
inline constexpr char TrackTrackTrackTriplet[] = "TrackTrackTrackTriplet";
inline constexpr char TrackTrackLambdaTriplet[] = "TrackTrackLambdaTriplet";
inline constexpr char TrackTrackXiTriplet[] = "TrackTrackXiTriplet";
inline constexpr char TrackTrackOmegaTriplet[] = "TrackTrackOmegaTriplet";
} // namespace role

namespace devicename
{
constexpr const char* FemtoPairTrackTrack = "femto-pair-track-track";
constexpr const char* FemtoPairTrackV0 = "femto-pair-track-v0";
constexpr const char* FemtoPairTrackCascade = "femto-pair-track-cascade";
constexpr const char* FemtoPairTrackKink = "femto-pair-track-kink";
constexpr const char* FemtoPairTrackCharmHadron = "femto-pair-track-charm-hadron";
constexpr const char* FemtoPairV0V0 = "femto-pair-v0-v0";
constexpr const char* FemtoPairD0D0 = "femto-pair-d0-d0";
constexpr const char* FemtoPairTrackTwoTrackResonance = "femto-pair-track-two-track-resonance";
constexpr const char* FemtoPairTrackTwoPhotonResonance = "femto-pair-track-two-photon-resonance";
constexpr const char* FemtoPairV0TwoTrackResonance = "femto-pair-v0-two-track-resonance";
constexpr const char* FemtoTripletTrackTrackTrack = "femto-triplet-track-track-track";
constexpr const char* FemtoTripletTrackTrackV0 = "femto-triplet-track-track-v0";
constexpr const char* FemtoTripletTrackTrackCascade = "femto-triplet-track-track-cascade";
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
  int bit = 0;
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
      // Gated on whether ANY of the *WithTagging* process functions that actually consume this
      // bit is enabled (processMixedEventWithTagging itself, or one of the
      // ...AndMass/...AndEventShape/...AndFitVars combined variants -- all share the same role
      // and tag bit, since the tag-bit check and the extra joined columns are independent
      // concerns). An instance that enables none of them will never look this bit up, so
      // assigning it one anyway would just waste a slot out of the hard-capped budget
      // (MaxTagBits below) for nothing.
      if (hasAnyEnabledOptionStartingWith(device, "processMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackTrackPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackV0) != std::string::npos) {
      // Each species gates independently on its own *MixedEventWithTagging* switch family --
      // Lambda and K0short are activated completely independently of each other (e.g.
      // same-event-only for one species, tagged-mixed-event for the other is a legitimate
      // combination), and, same as TrackTrackPair above, same-event/Mc/WithFitVars-without-tagging
      // variants are irrelevant here since none of them ever consume a tag bit.
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackLambdaPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processK0shortMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackK0shortPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackCascade) != std::string::npos) {
      // Xi and Omega are mutually exclusive on a single instance (the task itself fatals if
      // both are enabled), but gating independently costs nothing and needs no special-casing
      // here.
      if (hasAnyEnabledOptionStartingWith(device, "processXiMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackXiPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processOmegaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackOmegaPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackKink) != std::string::npos) {
      if (hasAnyEnabledOptionStartingWith(device, "processSigmaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackSigmaPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processSigmaPlusMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackSigmaPlusPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processSigmaPlusPcmMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackSigmaPlusPcmPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackCharmHadron) != std::string::npos) {
      if (hasAnyEnabledOptionStartingWith(device, "processD0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackD0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processLcMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackLcPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairV0V0) != std::string::npos) {
      // Lambda-lambda, k0short-k0short and lambda-k0short are mutually exclusive on a single
      // instance (the task itself fatals if more than one is enabled), but gating independently
      // costs nothing and needs no special-casing here.
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaLambdaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::LambdaLambdaPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processK0shortK0shortMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::K0shortK0shortPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaK0shortMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::LambdaK0shortPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairD0D0) != std::string::npos) {
      // Only one pair type possible per instance (sign on D0Selection1/2 picks the hypothesis), so
      // this task's WithTagging process function is just "processMixedEventWithTagging", same local
      // name as femtoPairTrackTrack.cxx's -- fine since option names are scoped per-device.
      if (hasAnyEnabledOptionStartingWith(device, "processMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::D0D0Pair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackTwoTrackResonance) != std::string::npos) {
      // Phi, Kstar0 and Rho0 are mutually exclusive on a single instance (the task itself fatals
      // if more than one is enabled), but gating independently costs nothing here either.
      if (hasAnyEnabledOptionStartingWith(device, "processPhiMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackPhiPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processKstar0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackKstar0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processRho0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackRho0Pair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairTrackTwoPhotonResonance) != std::string::npos) {
      // Pi0 and Eta are mutually exclusive on a single instance, same as above.
      if (hasAnyEnabledOptionStartingWith(device, "processPi0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackPi0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processEtaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackEtaPair), counter++});
      }
    } else if (device.name.find(devicename::FemtoPairV0TwoTrackResonance) != std::string::npos) {
      // All six v0 x resonance combinations are mutually exclusive on a single instance (the task
      // itself fatals if more than one is enabled), same as above.
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaPhiMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::LambdaPhiPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaKstar0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::LambdaKstar0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processLambdaRho0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::LambdaRho0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processK0shortPhiMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::K0shortPhiPair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processK0shortKstar0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::K0shortKstar0Pair), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processK0shortRho0MixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::K0shortRho0Pair), counter++});
      }
    } else if (device.name.find(devicename::FemtoTripletTrackTrackTrack) != std::string::npos) {
      if (hasAnyEnabledOptionStartingWith(device, "processMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackTrackTrackTriplet), counter++});
      }
    } else if (device.name.find(devicename::FemtoTripletTrackTrackV0) != std::string::npos) {
      if (hasAnyEnabledOptionStartingWith(device, "processMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackTrackLambdaTriplet), counter++});
      }
    } else if (device.name.find(devicename::FemtoTripletTrackTrackCascade) != std::string::npos) {
      // Xi and Omega are mutually exclusive on a single instance (the task itself fatals if both
      // are enabled), but gating independently costs nothing here either.
      if (hasAnyEnabledOptionStartingWith(device, "processXiMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackTrackXiTriplet), counter++});
      }
      if (hasAnyEnabledOptionStartingWith(device, "processOmegaMixedEventWithTagging")) {
        roles.push_back({device.name, static_cast<char const*>(role::TrackTrackOmegaTriplet), counter++});
      }
    }
    // TODO: add the remaining triplet task types (mc-particle triplets, ...), following the same
    // pattern. Whatever is added here must also be added to femtoProducerCollisionTagger.cxx's
    // init(), in the same order, since that is where the selection VALUES behind each of these
    // bits actually get read and checked.
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
// the bit that role was assigned, discovered once in init() into `internalTagMask` (not read from
// the JSON config). Bundling both together means a call like findOwnBit(roles, ownDevice.name,
// confTagging.Role) always uses the exact role this same object's mask gets set for, instead of
// the role being a second, independently-spelled literal at the call site. Every task that opts
// into tagger-filtered mixing instantiates one of these; wrapped in a ConfigurableGroup (rather
// than a bare Configurable) so more fields can be added here later without having to touch every
// consumer. `prefix` folds in Role itself: a task with more than one pair type (e.g.
// femtoPairTrackV0.cxx's Lambda and K0short, each needing their own tag bit) instantiates
// ConfTagging twice, and without this, both would register the exact same CLI/JSON option name
// ("Tagging.internalTagMask" twice over), which DPL rejects at startup as an ambiguous option.
template <char const* role>
struct ConfTagging : o2::framework::ConfigurableGroup {
  std::string prefix = std::string("Tagging") + role;
  static constexpr char const* Role = role;
  o2::framework::Configurable<datatypes::CollisionTagType> internalTagMask{"internalTagMask", 0, "internal: discovered collision tagger bitmask, set automatically in init() -- do not set manually"};
};

} // namespace o2::analysis::femto::collisiontagger

#endif // PWGCF_FEMTO_CORE_COLLISIONTAGGER_H_
