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

/// \file HFPvRefitSkimConverted.cxx
/// \brief Dummy workflow to add PVrefit tables to Run 2 converted data with original values
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN Padova, Italy

#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "Common/Core/trackUtilities.h"

using namespace o2::framework;

#include "Framework/runDataProcessing.h"

struct HFPvRefitSkimConverted {

  // Tables to be added
  Produces<aod::HfPvRefitProng2> rowProng2PVrefit;
  Produces<aod::HfPvRefitProng3> rowProng3PVrefit;

  // process function
  void process(aod::Collisions const& collisions,
               aod::Hf2Prongs const& rowsTrackIndexProng2,
               aod::Hf3Prongs const& rowsTrackIndexProng3,
               aod::Tracks const&)
  {

    // dummy tables for 2 prong candidates
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      // original PV information
      auto track0 = rowTrackIndexProng2.index0_as<aod::Tracks>();
      auto collision = track0.collision();
      auto primaryVertex = getPrimaryVertex(collision);

      // fill table row with coordinates of origial PV (dummy)
      rowProng2PVrefit(primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(), primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2());
    }

    // dummy tables for 3 prong candidates
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

      // original PV information
      auto track0 = rowTrackIndexProng3.index0_as<aod::Tracks>();
      auto collision = track0.collision();
      auto primaryVertex = getPrimaryVertex(collision);

      // fill table row with coordinates of origial PV (dummy)
      rowProng3PVrefit(primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(), primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2());
    }
  } // end of process function
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<HFPvRefitSkimConverted>(cfgc)};
  return workflow;
}