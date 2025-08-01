// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file femtodreamConverter.cxx
/// \brief converter task for femtoDream to femtoUnited
/// \author Anton Riedel, TU München, anton.riedel@tum.de

#include "PWGCF/DataModel/FemtoDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoCollisionsDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct FemtodreamConverter {

  // femto collisions
  Produces<FUCols> outputCollisions;

  // femto tracks
  Produces<FUTracks> outputTracks;
  Produces<FUTrackMasks> outputTrackMasks;
  // Produces<FTrackDCAs> outputTrackDCAs;
  // Produces<FTrackExtras> outputTrackExtras;
  //
  // // femto vzeros
  Produces<FULambdas> outputLambdas;
  // Produces<FLambdaExtras> outputLambdaExtras;
  //
  // Produces<FULambdaDaus> outputLambdaDaus;

  // void init(o2::framework::InitContext&){};

  void process(FDCollision const& col, FDParticles const& parts)
  {
    outputCollisions(col.posZ(), col.multNtr(), col.multV0M(), col.sphericity(), col.magField());
    for (auto const& part : parts) {
      // fill tracks
      if (part.partType() == femtodreamparticle::ParticleType::kTrack) {
        // use last index for the collision index
        // the index stored in old format might be different, but it is only about having a
        // unique index for every collision
        int sign = 1;
        if (part.cut() & 1) {
          sign = -1;
        }
        outputTracks(outputCollisions.lastIndex(), sign * part.pt(), part.eta(), part.phi());
        // in old format tpc and tpctof information is stored in the same bitmask
        // so we only fill the TPC mask in the FemtoUnited framework
        // for analysis the bitmask has to be computed with the old tools
        uint64_t combined = (static_cast<uint64_t>(part.cut()) << 32) | part.pidcut();
        outputTrackMasks(combined);

        // fill Lambdas
        // when Lambdas are filled the daughters will also be filled, so they are skipped in the outer loop
      } else if (part.partType() == femtodreamparticle::ParticleType::kV0) {
        // vzero
        const auto& posChild = parts.iteratorAt(part.index() - 2);
        const auto& negChild = parts.iteratorAt(part.index() - 1);
        int sign = 1;
        if (part.cut() & 1) {
          sign = -1;
        }
        if (sign > 0) {
          outputLambdas(outputCollisions.lastIndex(), part.pt(), part.eta(), part.phi(), part.mLambda(), posChild.childrenIds()[0], negChild.childrenIds()[1]);
        } else {
          outputLambdas(outputCollisions.lastIndex(), -part.pt(), part.eta(), part.phi(), part.mAntiLambda(), posChild.childrenIds()[0], negChild.childrenIds()[1]);
        }
        // vzero daughters
      }
    }
  };
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FemtodreamConverter>(cfgc)};
  return workflow;
}
