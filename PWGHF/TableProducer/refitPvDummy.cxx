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

/// \file refitPvDummy.cxx
/// \brief Dummy workflow to add PVrefit tables to Run 2 converted data with original values
///
/// \author Mattia Faggin <mattia.faggin@cern.ch>, University and INFN Padova, Italy

#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/trackUtilities.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct HfRefitPvDummy {

  // Tables to be added
  Produces<aod::HfPvRefit2Prong> rowProng2PVrefit;
  Produces<aod::HfPvRefit3Prong> rowProng3PVrefit;
  Produces<aod::Hf2Prongs_001> rowProng2WithCollId;
  Produces<aod::Hf3Prongs_001> rowProng3WithCollId;
  Produces<aod::HfCascades_001> rowCascadesWithCollId;

  // process function for 2-prong candidates
  void process2Prongs(aod::Collisions const&,
                      aod::Hf2Prongs_000 const& rowsTrackIndexProng2,
                      aod::Tracks const&)
  {
    // dummy tables for 2 prong candidates
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {

      // original PV information
      auto track0 = rowTrackIndexProng2.prong0();
      auto collision = track0.collision();
      auto primaryVertex = getPrimaryVertex(collision);

      // fill table row with coordinates of origial PV (dummy)
      rowProng2PVrefit(primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(), primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2());

      // fill version 001 of candidate table, containing also the collision ID
      rowProng2WithCollId(collision.globalIndex(), rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id(), rowTrackIndexProng2.hfflag());
    }
  } // end of process function for 2 prong candidates

  PROCESS_SWITCH(HfRefitPvDummy, process2Prongs, "Produce tables for 2-prong candidates", true);

  // process function for 3-prong candidates
  void process3Prongs(aod::Collisions const&,
                      aod::Hf3Prongs_000 const& rowsTrackIndexProng3,
                      aod::Tracks const&)
  {
    // dummy tables for 3 prong candidates
    for (const auto& rowTrackIndexProng3 : rowsTrackIndexProng3) {

      // original PV information
      auto track0 = rowTrackIndexProng3.prong0();
      auto collision = track0.collision();
      auto primaryVertex = getPrimaryVertex(collision);

      // fill table row with coordinates of origial PV (dummy)
      rowProng3PVrefit(primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ(),
                       primaryVertex.getSigmaX2(), primaryVertex.getSigmaXY(), primaryVertex.getSigmaY2(), primaryVertex.getSigmaXZ(), primaryVertex.getSigmaYZ(), primaryVertex.getSigmaZ2());

      // fill version 001 of candidate table, containing also the collision ID
      rowProng3WithCollId(collision.globalIndex(), rowTrackIndexProng3.prong0Id(), rowTrackIndexProng3.prong1Id(), rowTrackIndexProng3.prong2Id(), rowTrackIndexProng3.hfflag());
    }
  } // end of process function for 3 prong candidates

  PROCESS_SWITCH(HfRefitPvDummy, process3Prongs, "Produce tables for 3-prong candidates", true);

  // process function for cascade candidates
  void processCascades(aod::Collisions const&,
                       aod::HfCascades_000 const& rowsTrackIndexCascade,
                       aod::Tracks const&)
  {
    // dummy tables for cascade candidates
    for (const auto& rowTrackIndexCascade : rowsTrackIndexCascade) {

      // original PV information
      auto track0 = rowTrackIndexCascade.prong0();
      auto collision = track0.collision();
      // fill version 001 of candidate table, containing also the collision ID
      rowCascadesWithCollId(collision.globalIndex(), rowTrackIndexCascade.prong0Id(), rowTrackIndexCascade.v0Id());
    }
  } // end of process function for cascade candidates

  PROCESS_SWITCH(HfRefitPvDummy, processCascades, "Produce tables for hf-cascade candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<HfRefitPvDummy>(cfgc)};
  return workflow;
}
