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

#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>

#include <string>
#include <vector>

using namespace o2::analysis::femto;

struct FemtoProducerCollisionTagger {
  o2::framework::Produces<o2::aod::FColTags> producedTags;

  // Process functions below take a single collision ITERATOR (not the whole table) as their
  // first argument, with the track/v0 tables as later arguments -- the same auto-grouping-by-
  // fColId mechanism every processSameEvent() in this codebase already relies on, so tracks/
  // lambdas/k0shorts arrive pre-sliced to that one collision with no manual Preslice/SliceCache/
  // sliceBy() needed here at all.
  using FemtoCollisions = o2::soa::Join<o2::aod::FCols, o2::aod::FColMasks>;
  using FemtoCollision = FemtoCollisions::iterator;
  using FemtoTracks = o2::soa::Join<o2::aod::FTracks, o2::aod::FTrackMasks>;
  using FemtoLambdas = o2::soa::Join<o2::aod::FLambdas, o2::aod::FLambdaMasks>;
  using FemtoK0shorts = o2::soa::Join<o2::aod::FK0shorts, o2::aod::FK0shortMasks>;
  using FemtoXis = o2::soa::Join<o2::aod::FXis, o2::aod::FXiMasks>;
  using FemtoOmegas = o2::soa::Join<o2::aod::FOmegas, o2::aod::FOmegaMasks>;
  using FemtoSigmas = o2::soa::Join<o2::aod::FSigmas, o2::aod::FSigmaMasks>;
  using FemtoSigmaPlus = o2::soa::Join<o2::aod::FSigmaPlus, o2::aod::FSigmaPlusMasks>;
  using FemtoSigmaPlusPcm = o2::soa::Join<o2::aod::FSigmaPlusPcm, o2::aod::FSigmaPlusPcmMasks, o2::aod::FSigmaPlusPcmMlScores>;
  using FemtoD0s = o2::soa::Join<o2::aod::FD0s, o2::aod::FD0Masks>;
  using FemtoLcs = o2::soa::Join<o2::aod::FLcs, o2::aod::FLcMasks>;
  using FemtoPhis = o2::soa::Join<o2::aod::FPhis, o2::aod::FPhiMasks>;
  using FemtoKstar0s = o2::soa::Join<o2::aod::FKstar0s, o2::aod::FKstar0Masks>;
  using FemtoRho0s = o2::soa::Join<o2::aod::FRho0s, o2::aod::FRho0Masks>;
  using FemtoPi0s = o2::soa::Join<o2::aod::FPi0s, o2::aod::FPi0Masks>;
  using FemtoEtas = o2::soa::Join<o2::aod::FEtas, o2::aod::FEtaMasks>;

  // --- femto-pair-track-track: one bit per instance, "a track-track pair is possible here" ---
  std::vector<collisiontagger::TrackTrackVariant> mTrackTrackVariants;

  // --- femto-pair-track-v0: one bit per instance per enabled v0 species ---
  std::vector<collisiontagger::TrackV0Variant> mTrackLambdaVariants;
  std::vector<collisiontagger::TrackV0Variant> mTrackK0shortVariants;

  // --- femto-pair-track-cascade: one bit per instance per enabled cascade species ---
  std::vector<collisiontagger::TrackV0Variant> mTrackXiVariants;
  std::vector<collisiontagger::TrackV0Variant> mTrackOmegaVariants;

  // --- femto-pair-track-kink: one bit per instance per enabled kink species ---
  std::vector<collisiontagger::TrackV0Variant> mTrackSigmaVariants;
  std::vector<collisiontagger::TrackV0Variant> mTrackSigmaPlusVariants;
  std::vector<collisiontagger::TrackKinkPcmVariant> mTrackSigmaPlusPcmVariants;

  // --- femto-pair-track-charm-hadron: one bit per instance per enabled charm hadron species ---
  std::vector<collisiontagger::TrackV0Variant> mTrackD0Variants;
  std::vector<collisiontagger::TrackV0Variant> mTrackLcVariants;

  // --- femto-pair-v0-v0: one bit per instance per enabled v0-v0 combination ---
  std::vector<collisiontagger::V0V0Variant> mLambdaLambdaVariants;
  std::vector<collisiontagger::V0V0Variant> mK0shortK0shortVariants;
  std::vector<collisiontagger::V0V0Variant> mLambdaK0shortVariants;

  // --- femto-pair-d0-d0: one bit per instance ---
  std::vector<collisiontagger::V0V0Variant> mD0D0Variants;

  // --- femto-pair-track-two-track-resonance: one bit per instance per enabled species ---
  std::vector<collisiontagger::TrackTwoTrackResonanceVariant> mTrackPhiVariants;
  std::vector<collisiontagger::TrackKstar0Variant> mTrackKstar0Variants;
  std::vector<collisiontagger::TrackTwoTrackResonanceVariant> mTrackRho0Variants;

  // --- femto-pair-track-two-photon-resonance: one bit per instance per enabled species ---
  std::vector<collisiontagger::TrackTwoPhotonResonanceVariant> mTrackPi0Variants;
  std::vector<collisiontagger::TrackTwoPhotonResonanceVariant> mTrackEtaVariants;

  // --- femto-pair-v0-two-track-resonance: one bit per instance per enabled combination ---
  std::vector<collisiontagger::V0TwoTrackResonanceVariant> mLambdaPhiVariants;
  std::vector<collisiontagger::V0Kstar0Variant> mLambdaKstar0Variants;
  std::vector<collisiontagger::V0TwoTrackResonanceVariant> mLambdaRho0Variants;
  std::vector<collisiontagger::V0TwoTrackResonanceVariant> mK0shortPhiVariants;
  std::vector<collisiontagger::V0Kstar0Variant> mK0shortKstar0Variants;
  std::vector<collisiontagger::V0TwoTrackResonanceVariant> mK0shortRho0Variants;

  // --- femto-triplet-track-track-track: one bit per instance ---
  std::vector<collisiontagger::TrackTrackTrackVariant> mTrackTrackTrackVariants;

  // --- femto-triplet-track-track-v0: one bit per instance ---
  std::vector<collisiontagger::TrackTrackV0Variant> mTrackTrackLambdaVariants;

  // --- femto-triplet-track-track-cascade: one bit per instance per enabled cascade species ---
  std::vector<collisiontagger::TrackTrackV0Variant> mTrackTrackXiVariants;
  std::vector<collisiontagger::TrackTrackV0Variant> mTrackTrackOmegaVariants;

  void init(o2::framework::InitContext& context)
  {
    // FColTags gets exactly one row written per collision; if more than one process function were
    // enabled, each would try to write its own row for the same collision. Until FColTags is
    // extended to carry one column per wagon type, only one process function may run at a time --
    // run separate tagger instances (with separate output files/merges) for separate wagon types.
    int activeProcesses = static_cast<int>(doprocessTrackTrack) + static_cast<int>(doprocessTrackLambda) + static_cast<int>(doprocessTrackK0short) + static_cast<int>(doprocessTrackXi) + static_cast<int>(doprocessTrackOmega) + static_cast<int>(doprocessTrackSigma) + static_cast<int>(doprocessTrackSigmaPlus) + static_cast<int>(doprocessTrackSigmaPlusPcm) + static_cast<int>(doprocessTrackD0) + static_cast<int>(doprocessTrackLc) + static_cast<int>(doprocessLambdaLambda) + static_cast<int>(doprocessK0shortK0short) + static_cast<int>(doprocessLambdaK0short) + static_cast<int>(doprocessD0D0) + static_cast<int>(doprocessTrackTwoTrackResonance) + static_cast<int>(doprocessTrackTwoPhotonResonance) + static_cast<int>(doprocessLambdaTwoTrackResonance) + static_cast<int>(doprocessK0shortTwoTrackResonance) + static_cast<int>(doprocessTripletTrackTrackTrack) + static_cast<int>(doprocessTripletTrackTrackV0) + static_cast<int>(doprocessTripletTrackTrackXi) + static_cast<int>(doprocessTripletTrackTrackOmega);
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
      if (doprocessTrackTrack && r.role == static_cast<char const*>(collisiontagger::role::TrackTrackPair) && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackTrack) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackTrackVariant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sameSpecies = collisiontagger::hasEnabledOption(device, "Mixing.sameSpecies");
        // When sameSpecies is on, the real pair task's own dispatcher (PairTrackTrackBuilder::
        // processMixedEvent in pairBuilder.h) uses partition1 for BOTH roles and never consults
        // TrackSelection2 at all -- mirror that here instead of reading TrackSelection2 off the
        // device, since the user is only expected to configure TrackSelection1 in that mode and
        // TrackSelection2 may be left at an unrelated default.
        v.sel2 = v.sameSpecies ? v.sel1 : collisiontagger::readTrackSelection(device, "TrackSelection2");
        v.bit = r.bit;
        mTrackTrackVariants.push_back(v);
      } else if (doprocessTrackLambda && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackV0) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackLambdaPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "LambdaSelection1");
        v.bit = r.bit;
        mTrackLambdaVariants.push_back(v);
      } else if (doprocessTrackK0short && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackV0) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackK0shortPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "K0shortSelection1");
        v.bit = r.bit;
        mTrackK0shortVariants.push_back(v);
      } else if (doprocessTrackXi && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackCascade) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackXiPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "XiSelection");
        v.bit = r.bit;
        mTrackXiVariants.push_back(v);
      } else if (doprocessTrackOmega && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackCascade) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackOmegaPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "OmegaSelection");
        v.bit = r.bit;
        mTrackOmegaVariants.push_back(v);
      } else if (doprocessTrackSigma && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackKink) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackSigmaPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "SigmaSelection1");
        v.bit = r.bit;
        mTrackSigmaVariants.push_back(v);
      } else if (doprocessTrackSigmaPlus && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackKink) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackSigmaPlusPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "SigmaPlusSelection1");
        v.bit = r.bit;
        mTrackSigmaPlusVariants.push_back(v);
      } else if (doprocessTrackSigmaPlusPcm && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackKink) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::TrackSigmaPlusPcmPair)) {
          collisiontagger::TrackKinkPcmVariant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readKinkPcmSelection(device, "SigmaPlusPcmSelection1");
          v.bit = r.bit;
          mTrackSigmaPlusPcmVariants.push_back(v);
        }
      } else if ((doprocessTrackD0 || doprocessTrackLc) && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackCharmHadron) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::TrackD0Pair)) {
          collisiontagger::TrackV0Variant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readV0Selection(device, "D0Selection1");
          v.bit = r.bit;
          mTrackD0Variants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::TrackLcPair)) {
          collisiontagger::TrackV0Variant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readV0Selection(device, "LcSelection1");
          v.bit = r.bit;
          mTrackLcVariants.push_back(v);
        }
      } else if (doprocessLambdaLambda && r.deviceName.find(collisiontagger::devicename::FemtoPairV0V0) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::LambdaLambdaPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::V0V0Variant v;
        v.sel1 = collisiontagger::readV0Selection(device, "LambdaSelection1");
        v.sameSpecies = collisiontagger::hasEnabledOption(device, "Mixing.sameSpecies");
        // Mirror PairV0V0Builder::processMixedEvent(), which uses partition1 for BOTH roles and
        // never consults LambdaSelection2 at all when sameSpecies -- same override as
        // TrackTrackPair above.
        v.sel2 = v.sameSpecies ? v.sel1 : collisiontagger::readV0Selection(device, "LambdaSelection2");
        v.bit = r.bit;
        mLambdaLambdaVariants.push_back(v);
      } else if (doprocessK0shortK0short && r.deviceName.find(collisiontagger::devicename::FemtoPairV0V0) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::K0shortK0shortPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::V0V0Variant v;
        v.sel1 = collisiontagger::readV0Selection(device, "K0shortSelection1");
        v.sameSpecies = collisiontagger::hasEnabledOption(device, "Mixing.sameSpecies");
        v.sel2 = v.sameSpecies ? v.sel1 : collisiontagger::readV0Selection(device, "K0shortSelection2");
        v.bit = r.bit;
        mK0shortK0shortVariants.push_back(v);
      } else if (doprocessLambdaK0short && r.deviceName.find(collisiontagger::devicename::FemtoPairV0V0) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::LambdaK0shortPair)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::V0V0Variant v;
        v.sel1 = collisiontagger::readV0Selection(device, "LambdaSelection1");
        v.sel2 = collisiontagger::readV0Selection(device, "K0shortSelection2");
        v.bit = r.bit;
        mLambdaK0shortVariants.push_back(v);
      } else if (doprocessD0D0 && r.deviceName.find(collisiontagger::devicename::FemtoPairD0D0) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::V0V0Variant v;
        v.sel1 = collisiontagger::readV0Selection(device, "D0Selection1");
        collisiontagger::V0SelectionValues sel2Raw = collisiontagger::readV0Selection(device, "D0Selection2");
        bool confSameSpecies = collisiontagger::hasEnabledOption(device, "Mixing.sameSpecies");
        // PairD0D0Builder::init() forces sameSpecies off when the two signs disagree (see
        // pairBuilder.h) -- mirror that override here too, otherwise a sign-flip-only d0-d0bar
        // config would be mis-tagged as a same-species pair.
        v.sameSpecies = confSameSpecies && (v.sel1.sign == sel2Raw.sign);
        v.sel2 = v.sameSpecies ? v.sel1 : sel2Raw;
        v.bit = r.bit;
        mD0D0Variants.push_back(v);
      } else if (doprocessTrackTwoTrackResonance && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackTwoTrackResonance) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::TrackPhiPair)) {
          collisiontagger::TrackTwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "PhiSelection");
          v.bit = r.bit;
          mTrackPhiVariants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::TrackKstar0Pair)) {
          collisiontagger::TrackKstar0Variant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readKstar0Selection(device, "Kstar0Selection");
          v.bit = r.bit;
          mTrackKstar0Variants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::TrackRho0Pair)) {
          collisiontagger::TrackTwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "Rho0Selection");
          v.bit = r.bit;
          mTrackRho0Variants.push_back(v);
        }
      } else if (doprocessTrackTwoPhotonResonance && r.deviceName.find(collisiontagger::devicename::FemtoPairTrackTwoPhotonResonance) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::TrackPi0Pair)) {
          collisiontagger::TrackTwoPhotonResonanceVariant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readTwoPhotonResonanceSelection(device, "Pi0Selection1");
          v.bit = r.bit;
          mTrackPi0Variants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::TrackEtaPair)) {
          collisiontagger::TrackTwoPhotonResonanceVariant v;
          v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
          v.sel2 = collisiontagger::readTwoPhotonResonanceSelection(device, "EtaSelection1");
          v.bit = r.bit;
          mTrackEtaVariants.push_back(v);
        }
      } else if (doprocessLambdaTwoTrackResonance && r.deviceName.find(collisiontagger::devicename::FemtoPairV0TwoTrackResonance) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::LambdaPhiPair)) {
          collisiontagger::V0TwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readV0Selection(device, "LambdaSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "PhiSelection");
          v.bit = r.bit;
          mLambdaPhiVariants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::LambdaKstar0Pair)) {
          collisiontagger::V0Kstar0Variant v;
          v.sel1 = collisiontagger::readV0Selection(device, "LambdaSelection1");
          v.sel2 = collisiontagger::readKstar0Selection(device, "Kstar0Selection");
          v.bit = r.bit;
          mLambdaKstar0Variants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::LambdaRho0Pair)) {
          collisiontagger::V0TwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readV0Selection(device, "LambdaSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "Rho0Selection");
          v.bit = r.bit;
          mLambdaRho0Variants.push_back(v);
        }
      } else if (doprocessK0shortTwoTrackResonance && r.deviceName.find(collisiontagger::devicename::FemtoPairV0TwoTrackResonance) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        if (r.role == static_cast<char const*>(collisiontagger::role::K0shortPhiPair)) {
          collisiontagger::V0TwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readV0Selection(device, "K0shortSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "PhiSelection");
          v.bit = r.bit;
          mK0shortPhiVariants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::K0shortKstar0Pair)) {
          collisiontagger::V0Kstar0Variant v;
          v.sel1 = collisiontagger::readV0Selection(device, "K0shortSelection1");
          v.sel2 = collisiontagger::readKstar0Selection(device, "Kstar0Selection");
          v.bit = r.bit;
          mK0shortKstar0Variants.push_back(v);
        } else if (r.role == static_cast<char const*>(collisiontagger::role::K0shortRho0Pair)) {
          collisiontagger::V0TwoTrackResonanceVariant v;
          v.sel1 = collisiontagger::readV0Selection(device, "K0shortSelection1");
          v.sel2 = collisiontagger::readTwoTrackResonanceSelection(device, "Rho0Selection");
          v.bit = r.bit;
          mK0shortRho0Variants.push_back(v);
        }
      } else if (doprocessTripletTrackTrackTrack && r.deviceName.find(collisiontagger::devicename::FemtoTripletTrackTrackTrack) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackTrackTrackVariant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        collisiontagger::TrackSelectionValues sel2Raw = collisiontagger::readTrackSelection(device, "TrackSelection2");
        collisiontagger::TrackSelectionValues sel3Raw = collisiontagger::readTrackSelection(device, "TrackSelection3");
        v.same123 = collisiontagger::hasEnabledOption(device, "Mixing.particle123AreSameSpecies");
        v.same12 = collisiontagger::hasEnabledOption(device, "Mixing.particle12AreSameSpecies");
        // Mirror TripletTrackTrackTrackBuilder::processMixedEvent(), which ignores
        // TrackSelection2/3 entirely when same123 (uses partition1 for all three roles), and
        // ignores TrackSelection2 (but not 3) when same12 -- same override trick as
        // TrackTrackVariant's sameSpecies above.
        v.sel2 = (v.same123 || v.same12) ? v.sel1 : sel2Raw;
        v.sel3 = v.same123 ? v.sel1 : sel3Raw;
        v.bit = r.bit;
        mTrackTrackTrackVariants.push_back(v);
      } else if (doprocessTripletTrackTrackV0 && r.deviceName.find(collisiontagger::devicename::FemtoTripletTrackTrackV0) != std::string::npos) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackTrackV0Variant v;
        v.sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        collisiontagger::TrackSelectionValues sel2Raw = collisiontagger::readTrackSelection(device, "TrackSelection2");
        v.same12 = collisiontagger::hasEnabledOption(device, "Mixing.particle12AreSameSpecies");
        v.sel2 = v.same12 ? v.sel1 : sel2Raw;
        v.sel3 = collisiontagger::readV0Selection(device, "LambdaSelection1");
        v.bit = r.bit;
        mTrackTrackLambdaVariants.push_back(v);
      } else if (doprocessTripletTrackTrackXi && r.deviceName.find(collisiontagger::devicename::FemtoTripletTrackTrackCascade) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackTrackXiTriplet)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackSelectionValues sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        collisiontagger::TrackSelectionValues sel2Raw = collisiontagger::readTrackSelection(device, "TrackSelection2");
        bool same12 = collisiontagger::hasEnabledOption(device, "Mixing.particle12AreSameSpecies");
        collisiontagger::TrackTrackV0Variant v;
        v.sel1 = sel1;
        v.sel2 = same12 ? sel1 : sel2Raw;
        v.same12 = same12;
        v.sel3 = collisiontagger::readV0Selection(device, "XiSelection");
        v.bit = r.bit;
        mTrackTrackXiVariants.push_back(v);
      } else if (doprocessTripletTrackTrackOmega && r.deviceName.find(collisiontagger::devicename::FemtoTripletTrackTrackCascade) != std::string::npos && r.role == static_cast<char const*>(collisiontagger::role::TrackTrackOmegaTriplet)) {
        o2::framework::DeviceSpec const& device = findDevice(r.deviceName);
        collisiontagger::TrackSelectionValues sel1 = collisiontagger::readTrackSelection(device, "TrackSelection1");
        collisiontagger::TrackSelectionValues sel2Raw = collisiontagger::readTrackSelection(device, "TrackSelection2");
        bool same12 = collisiontagger::hasEnabledOption(device, "Mixing.particle12AreSameSpecies");
        collisiontagger::TrackTrackV0Variant v;
        v.sel1 = sel1;
        v.sel2 = same12 ? sel1 : sel2Raw;
        v.same12 = same12;
        v.sel3 = collisiontagger::readV0Selection(device, "OmegaSelection");
        v.bit = r.bit;
        mTrackTrackOmegaVariants.push_back(v);
      }
      // TODO: add the remaining triplet task types (mc-particle triplets, ...), following the
      // same pattern.
    }

    if (doprocessTrackTrack && mTrackTrackVariants.empty()) {
      LOG(fatal) << "processTrackTrack is enabled, but no femto-pair-track-track device was found in this workflow.";
    }
    if (doprocessTrackLambda && mTrackLambdaVariants.empty()) {
      LOG(fatal) << "processTrackLambda is enabled, but no femto-pair-track-v0 device (with processLambdaMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackK0short && mTrackK0shortVariants.empty()) {
      LOG(fatal) << "processTrackK0short is enabled, but no femto-pair-track-v0 device (with processK0shortMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackXi && mTrackXiVariants.empty()) {
      LOG(fatal) << "processTrackXi is enabled, but no femto-pair-track-cascade device (with processXiMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackOmega && mTrackOmegaVariants.empty()) {
      LOG(fatal) << "processTrackOmega is enabled, but no femto-pair-track-cascade device (with processOmegaMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackSigma && mTrackSigmaVariants.empty()) {
      LOG(fatal) << "processTrackSigma is enabled, but no femto-pair-track-kink device (with processSigmaMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackSigmaPlus && mTrackSigmaPlusVariants.empty()) {
      LOG(fatal) << "processTrackSigmaPlus is enabled, but no femto-pair-track-kink device (with processSigmaPlusMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackSigmaPlusPcm && mTrackSigmaPlusPcmVariants.empty()) {
      LOG(fatal) << "processTrackSigmaPlusPcm is enabled, but no femto-pair-track-kink device (with "
                 << "processSigmaPlusPcmMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTrackD0 && mTrackD0Variants.empty()) {
      LOG(fatal) << "processTrackD0 is enabled, but no femto-pair-track-charm-hadron device (with processD0MixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessTrackLc && mTrackLcVariants.empty()) {
      LOG(fatal) << "processTrackLc is enabled, but no femto-pair-track-charm-hadron device (with processLcMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessLambdaLambda && mLambdaLambdaVariants.empty()) {
      LOG(fatal) << "processLambdaLambda is enabled, but no femto-pair-v0-v0 device (with processLambdaLambdaMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessK0shortK0short && mK0shortK0shortVariants.empty()) {
      LOG(fatal) << "processK0shortK0short is enabled, but no femto-pair-v0-v0 device (with processK0shortK0shortMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessLambdaK0short && mLambdaK0shortVariants.empty()) {
      LOG(fatal) << "processLambdaK0short is enabled, but no femto-pair-v0-v0 device (with processLambdaK0shortMixedEventWithTagging "
                 << "enabled) was found in this workflow.";
    }
    if (doprocessD0D0 && mD0D0Variants.empty()) {
      LOG(fatal) << "processD0D0 is enabled, but no femto-pair-d0-d0 device (with processMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTrackTwoTrackResonance && mTrackPhiVariants.empty() && mTrackKstar0Variants.empty() && mTrackRho0Variants.empty()) {
      LOG(fatal) << "processTrackTwoTrackResonance is enabled, but no femto-pair-track-two-track-resonance device (with "
                 << "processPhiMixedEventWithTagging/processKstar0MixedEventWithTagging/processRho0MixedEventWithTagging enabled) "
                 << "was found in this workflow.";
    }
    if (doprocessTrackTwoPhotonResonance && mTrackPi0Variants.empty() && mTrackEtaVariants.empty()) {
      LOG(fatal) << "processTrackTwoPhotonResonance is enabled, but no femto-pair-track-two-photon-resonance device (with "
                 << "processPi0MixedEventWithTagging/processEtaMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessLambdaTwoTrackResonance && mLambdaPhiVariants.empty() && mLambdaKstar0Variants.empty() && mLambdaRho0Variants.empty()) {
      LOG(fatal) << "processLambdaTwoTrackResonance is enabled, but no femto-pair-v0-two-track-resonance device (with one of "
                 << "processLambdaPhiMixedEventWithTagging/processLambdaKstar0MixedEventWithTagging/"
                 << "processLambdaRho0MixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessK0shortTwoTrackResonance && mK0shortPhiVariants.empty() && mK0shortKstar0Variants.empty() && mK0shortRho0Variants.empty()) {
      LOG(fatal) << "processK0shortTwoTrackResonance is enabled, but no femto-pair-v0-two-track-resonance device (with one of "
                 << "processK0shortPhiMixedEventWithTagging/processK0shortKstar0MixedEventWithTagging/"
                 << "processK0shortRho0MixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTripletTrackTrackTrack && mTrackTrackTrackVariants.empty()) {
      LOG(fatal) << "processTripletTrackTrackTrack is enabled, but no femto-triplet-track-track-track device (with "
                 << "processMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTripletTrackTrackV0 && mTrackTrackLambdaVariants.empty()) {
      LOG(fatal) << "processTripletTrackTrackV0 is enabled, but no femto-triplet-track-track-v0 device (with "
                 << "processMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTripletTrackTrackXi && mTrackTrackXiVariants.empty()) {
      LOG(fatal) << "processTripletTrackTrackXi is enabled, but no femto-triplet-track-track-cascade device (with "
                 << "processXiMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (doprocessTripletTrackTrackOmega && mTrackTrackOmegaVariants.empty()) {
      LOG(fatal) << "processTripletTrackTrackOmega is enabled, but no femto-triplet-track-track-cascade device (with "
                 << "processOmegaMixedEventWithTagging enabled) was found in this workflow.";
    }
    if (!roles.empty() && roles.back().bit >= collisiontagger::MaxTagBits) {
      LOG(fatal) << "Too many tag bits requested (" << (roles.back().bit + 1) << "); collision tagger supports up to " << collisiontagger::MaxTagBits << ".";
    }
  }

  // `tracksThisCol` arrives already sliced to `col` (see FemtoCollision comment above) -- the
  // entire per-collision tagging step is one call to tagAndProduce(), which builds the tag
  // bitset, runs the same-type-pair check (both track roles' selections, or 2 distinct entries
  // for a same-species pair) against every variant, and produces the row.
  void processTrackTrack(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol)
  {
    collisiontagger::tagAndProduce(producedTags, mTrackTrackVariants, tracksThisCol,
                                   [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); });
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackTrack, "Tag collisions for femto-pair-track-track wagon(s)", false);

  // Lambda and K0short are split into two separate process functions (rather than one combined
  // processTrackV0 reading both FemtoLambdas and FemtoK0shorts unconditionally) because the
  // consuming femto-pair-track-v0 task itself enforces that a single instance only ever tags
  // Lambda OR K0short, never both (see its init(): "Both lambda-track and k0short-track
  // processing is enabled" is fatal there) -- so requiring both tables as DPL inputs here would
  // make the tagger unusable against any input file that only has one of the two V0 species
  // derived (e.g. Lambda-only productions with no FK0shorts table at all: DPL subscribes to a
  // process function's full argument list regardless of whether the branch using it is actually
  // reached at runtime, so the AOD reader fails outright looking for a tree that was never
  // produced, instead of this cleanly tagging zero K0short pairs).
  void processTrackLambda(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoLambdas const& lambdasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto trackPassesFn = [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); };
    collisiontagger::tagCrossVariants(tag, mTrackLambdaVariants, tracksThisCol, trackPassesFn, lambdasThisCol,
                                      [](auto const& lambda, auto const& sel) { return collisiontagger::signedV0Passes(lambda, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackLambda, "Tag collisions for femto-pair-track-v0 wagon(s) (Lambda)", false);

  void processTrackK0short(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoK0shorts const& k0shortsThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto trackPassesFn = [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); };
    collisiontagger::tagCrossVariants(tag, mTrackK0shortVariants, tracksThisCol, trackPassesFn, k0shortsThisCol,
                                      [](auto const& k0short, auto const& sel) { return collisiontagger::k0shortPasses(k0short, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackK0short, "Tag collisions for femto-pair-track-v0 wagon(s) (K0short)", false);

  // Xi and Omega are split into separate process functions (rather than one combined
  // processTrackCascade reading both FemtoXis and FemtoOmegas unconditionally): same reasoning as
  // processTrackLambda/processTrackK0short above -- femtoProducerDerivedToDerived.cxx's
  // processTracksOmegas has no corresponding "always comes with Xi" guarantee (Xi and Omega are
  // each their own independently-switchable, default-off producer process function), so a
  // Xi-only (or Omega-only) derived file would otherwise make the AOD reader fail looking for a
  // tree that was never produced.
  void processTrackXi(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoXis const& xisThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackXiVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, xisThisCol, [](auto const& cascade, auto const& sel) { return collisiontagger::signedV0Passes(cascade, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackXi, "Tag collisions for femto-pair-track-cascade wagon(s) (Xi)", false);

  void processTrackOmega(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoOmegas const& omegasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackOmegaVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, omegasThisCol, [](auto const& cascade, auto const& sel) { return collisiontagger::signedV0Passes(cascade, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackOmega, "Tag collisions for femto-pair-track-cascade wagon(s) (Omega)", false);

  // Sigma and SigmaPlus are likewise split: femtoProducerDerivedToDerived.cxx's processTracksSigma
  // and processTracksSigmaPlus are each their own independently-switchable, default-off process
  // function, not bundled -- same risk as Xi/Omega and Lambda/K0short above. SigmaPlusPcm already
  // had its own process function below (comes from a genuinely separate, mutually exclusive
  // producer process function), unaffected by this split.
  void processTrackSigma(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoSigmas const& sigmasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackSigmaVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, sigmasThisCol, [](auto const& sigma, auto const& sel) { return collisiontagger::signedV0Passes(sigma, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackSigma, "Tag collisions for femto-pair-track-kink wagon(s) (Sigma)", false);

  void processTrackSigmaPlus(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoSigmaPlus const& sigmaPlusThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackSigmaPlusVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, sigmaPlusThisCol, [](auto const& sigmaPlus, auto const& sel) { return collisiontagger::signedV0Passes(sigmaPlus, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackSigmaPlus, "Tag collisions for femto-pair-track-kink wagon(s) (SigmaPlus)", false);

  void processTrackSigmaPlusPcm(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoSigmaPlusPcm const& sigmaPlusPcmThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackSigmaPlusPcmVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, sigmaPlusPcmThisCol, [](auto const& sigmaPlusPcm, auto const& sel) { return collisiontagger::kinkPcmPasses(sigmaPlusPcm, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackSigmaPlusPcm, "Tag collisions for femto-pair-track-kink wagon(s) (SigmaPlusPcm)", false);

  // D0 and Lc are also mutually exclusive on the producer side (separate process functions, counted
  // together in its own "exactly one process enabled" check), so each gets its own process function
  // here too, for the same reason as SigmaPlusPcm above.
  void processTrackD0(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoD0s const& d0sThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackD0Variants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, d0sThisCol, [](auto const& d0, auto const& sel) { return collisiontagger::signedV0Passes(d0, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackD0, "Tag collisions for femto-pair-track-charm-hadron wagon(s) (D0)", false);

  void processTrackLc(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoLcs const& lcsThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mTrackLcVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, lcsThisCol, [](auto const& lc, auto const& sel) { return collisiontagger::signedV0Passes(lc, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackLc, "Tag collisions for femto-pair-track-charm-hadron wagon(s) (Lc)", false);

  // No track table involved here at all -- v0-v0 is checked purely between the two v0 tables
  // themselves (lambda/k0short), same as D0D0 below is checked purely within FemtoD0s.
  // Split by which table(s) are actually needed: LambdaLambda only needs FemtoLambdas,
  // K0shortK0short only needs FemtoK0shorts -- same reasoning as the other splits above
  // (femtoProducerDerivedToDerived.cxx's processTracksLambdas/processTracksK0shorts are each
  // independently switchable, default off). LambdaK0short is a genuine cross-species pair, so it
  // legitimately needs both tables and stays as-is.
  void processLambdaLambda(FemtoCollision const& /*col*/, FemtoLambdas const& lambdasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagPairVariants(tag, mLambdaLambdaVariants, lambdasThisCol,
                                     [](auto const& v0, auto const& sel) { return collisiontagger::signedV0Passes(v0, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processLambdaLambda, "Tag collisions for femto-pair-v0-v0 wagon(s) (Lambda-Lambda)", false);

  void processK0shortK0short(FemtoCollision const& /*col*/, FemtoK0shorts const& k0shortsThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagPairVariants(tag, mK0shortK0shortVariants, k0shortsThisCol,
                                     [](auto const& v0, auto const& sel) { return collisiontagger::k0shortPasses(v0, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processK0shortK0short, "Tag collisions for femto-pair-v0-v0 wagon(s) (K0short-K0short)", false);

  void processLambdaK0short(FemtoCollision const& /*col*/, FemtoLambdas const& lambdasThisCol, FemtoK0shorts const& k0shortsThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossVariants(tag, mLambdaK0shortVariants, lambdasThisCol, [](auto const& v0, auto const& sel) { return collisiontagger::signedV0Passes(v0, sel); }, k0shortsThisCol, [](auto const& v0, auto const& sel) { return collisiontagger::k0shortPasses(v0, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processLambdaK0short, "Tag collisions for femto-pair-v0-v0 wagon(s) (Lambda-K0short)", false);

  void processD0D0(FemtoCollision const& /*col*/, FemtoD0s const& d0sThisCol)
  {
    collisiontagger::tagAndProduce(producedTags, mD0D0Variants, d0sThisCol,
                                   [](auto const& d0, auto const& sel) { return collisiontagger::signedV0Passes(d0, sel); });
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processD0D0, "Tag collisions for femto-pair-d0-d0 wagon(s)", false);

  // femto-producer-resonances has no mutual-exclusivity gate between Phi/Kstar0/Rho0 (all three
  // default to enabled, built from the same track pool with different PID hypotheses), so they can
  // safely be combined into one process function here too.
  void processTrackTwoTrackResonance(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoPhis const& phisThisCol, FemtoKstar0s const& kstar0sThisCol, FemtoRho0s const& rho0sThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto trackPassesFn = [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); };
    if (!mTrackPhiVariants.empty()) {
      collisiontagger::tagCrossVariants(tag, mTrackPhiVariants, tracksThisCol, trackPassesFn, phisThisCol,
                                        [](auto const& phi, auto const& sel) { return collisiontagger::twoTrackResonancePasses(phi, sel); });
    }
    if (!mTrackKstar0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mTrackKstar0Variants, tracksThisCol, trackPassesFn, kstar0sThisCol,
                                        [](auto const& kstar0, auto const& sel) { return collisiontagger::kstar0Passes(kstar0, sel); });
    }
    if (!mTrackRho0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mTrackRho0Variants, tracksThisCol, trackPassesFn, rho0sThisCol,
                                        [](auto const& rho0, auto const& sel) { return collisiontagger::twoTrackResonancePasses(rho0, sel); });
    }
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackTwoTrackResonance, "Tag collisions for femto-pair-track-two-track-resonance wagon(s)", false);

  // femto-producer-two-photon-resonances also has no mutual-exclusivity gate between Pi0/Eta.
  void processTrackTwoPhotonResonance(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoPi0s const& pi0sThisCol, FemtoEtas const& etasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto trackPassesFn = [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); };
    if (!mTrackPi0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mTrackPi0Variants, tracksThisCol, trackPassesFn, pi0sThisCol,
                                        [](auto const& pi0, auto const& sel) { return collisiontagger::twoPhotonResonancePasses(pi0, sel); });
    }
    if (!mTrackEtaVariants.empty()) {
      collisiontagger::tagCrossVariants(tag, mTrackEtaVariants, tracksThisCol, trackPassesFn, etasThisCol,
                                        [](auto const& eta, auto const& sel) { return collisiontagger::twoPhotonResonancePasses(eta, sel); });
    }
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTrackTwoPhotonResonance, "Tag collisions for femto-pair-track-two-photon-resonance wagon(s)", false);

  // No track table involved here at all -- checked purely between the v0 tables (lambda/k0short)
  // and the resonance tables (phi/kstar0/rho0) themselves. Split by V0 species (Lambda vs
  // K0short) for the same reason as processLambdaLambda/processK0shortK0short above -- each
  // needs only its own V0 table. Phi/Kstar0/Rho0 stay combined within each (that axis is safe:
  // femtoProducerResonances.cxx builds all three by default with no mutual-exclusivity gate, see
  // processTrackTwoTrackResonance below).
  void processLambdaTwoTrackResonance(FemtoCollision const& /*col*/, FemtoLambdas const& lambdasThisCol, FemtoPhis const& phisThisCol, FemtoKstar0s const& kstar0sThisCol, FemtoRho0s const& rho0sThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto signedV0PassesFn = [](auto const& lambda, auto const& sel) { return collisiontagger::signedV0Passes(lambda, sel); };
    if (!mLambdaPhiVariants.empty()) {
      collisiontagger::tagCrossVariants(tag, mLambdaPhiVariants, lambdasThisCol, signedV0PassesFn, phisThisCol,
                                        [](auto const& phi, auto const& sel) { return collisiontagger::twoTrackResonancePasses(phi, sel); });
    }
    if (!mLambdaKstar0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mLambdaKstar0Variants, lambdasThisCol, signedV0PassesFn, kstar0sThisCol,
                                        [](auto const& kstar0, auto const& sel) { return collisiontagger::kstar0Passes(kstar0, sel); });
    }
    if (!mLambdaRho0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mLambdaRho0Variants, lambdasThisCol, signedV0PassesFn, rho0sThisCol,
                                        [](auto const& rho0, auto const& sel) { return collisiontagger::twoTrackResonancePasses(rho0, sel); });
    }
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processLambdaTwoTrackResonance, "Tag collisions for femto-pair-v0-two-track-resonance wagon(s) (Lambda)", false);

  void processK0shortTwoTrackResonance(FemtoCollision const& /*col*/, FemtoK0shorts const& k0shortsThisCol, FemtoPhis const& phisThisCol, FemtoKstar0s const& kstar0sThisCol, FemtoRho0s const& rho0sThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    auto k0shortPassesFn = [](auto const& k0short, auto const& sel) { return collisiontagger::k0shortPasses(k0short, sel); };
    if (!mK0shortPhiVariants.empty()) {
      collisiontagger::tagCrossVariants(tag, mK0shortPhiVariants, k0shortsThisCol, k0shortPassesFn, phisThisCol,
                                        [](auto const& phi, auto const& sel) { return collisiontagger::twoTrackResonancePasses(phi, sel); });
    }
    if (!mK0shortKstar0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mK0shortKstar0Variants, k0shortsThisCol, k0shortPassesFn, kstar0sThisCol,
                                        [](auto const& kstar0, auto const& sel) { return collisiontagger::kstar0Passes(kstar0, sel); });
    }
    if (!mK0shortRho0Variants.empty()) {
      collisiontagger::tagCrossVariants(tag, mK0shortRho0Variants, k0shortsThisCol, k0shortPassesFn, rho0sThisCol,
                                        [](auto const& rho0, auto const& sel) { return collisiontagger::twoTrackResonancePasses(rho0, sel); });
    }
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processK0shortTwoTrackResonance, "Tag collisions for femto-pair-v0-two-track-resonance wagon(s) (K0short)", false);

  void processTripletTrackTrackTrack(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagTripletVariants(tag, mTrackTrackTrackVariants, tracksThisCol,
                                        [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTripletTrackTrackTrack, "Tag collisions for femto-triplet-track-track-track wagon(s)", false);

  void processTripletTrackTrackV0(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoLambdas const& lambdasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossTripletVariants(tag, mTrackTrackLambdaVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, lambdasThisCol, [](auto const& lambda, auto const& sel) { return collisiontagger::signedV0Passes(lambda, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTripletTrackTrackV0, "Tag collisions for femto-triplet-track-track-v0 wagon(s)", false);

  // Split by Xi/Omega for the same reason as processTrackXi/processTrackOmega above.
  void processTripletTrackTrackXi(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoXis const& xisThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossTripletVariants(tag, mTrackTrackXiVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, xisThisCol, [](auto const& cascade, auto const& sel) { return collisiontagger::signedV0Passes(cascade, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTripletTrackTrackXi, "Tag collisions for femto-triplet-track-track-cascade wagon(s) (Xi)", false);

  void processTripletTrackTrackOmega(FemtoCollision const& /*col*/, FemtoTracks const& tracksThisCol, FemtoOmegas const& omegasThisCol)
  {
    collisiontagger::CollisionTagBitset tag;
    collisiontagger::tagCrossTripletVariants(tag, mTrackTrackOmegaVariants, tracksThisCol, [](auto const& track, auto const& sel) { return collisiontagger::trackPasses(track, sel); }, omegasThisCol, [](auto const& cascade, auto const& sel) { return collisiontagger::signedV0Passes(cascade, sel); });
    producedTags(static_cast<datatypes::CollisionTagType>(tag.to_ullong()));
  }
  PROCESS_SWITCH(FemtoProducerCollisionTagger, processTripletTrackTrackOmega, "Tag collisions for femto-triplet-track-track-cascade wagon(s) (Omega)", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{
    o2::framework::adaptAnalysisTask<FemtoProducerCollisionTagger>(context),
  };
  return workflow;
}
