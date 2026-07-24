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

/// \file femtoProducerLiteConverter.cxx
/// \brief Task that converts FLiteTracks (binned kinematics) back to FTracks (float kinematics)
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerLiteConverter {
  o2::framework::Produces<o2::aod::FCols> producedCols;
  o2::framework::Produces<o2::aod::FTracks> producedTracks;

  void init(o2::framework::InitContext&)
  {
  }

  void processLiteCols(o2::aod::FLiteCols::iterator const& liteCols)
  {
    producedCols(liteCols.posZ(),
                 liteCols.mult(),
                 liteCols.cent(),
                 liteCols.magField());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteCols, "Convert FLiteCols to FCols", true);

  void processLiteTracks(o2::aod::FLiteTracks::iterator const& liteTrack)
  {
    producedTracks(liteTrack.fColId(),
                   liteTrack.signedPt(),
                   liteTrack.eta(),
                   liteTrack.phi());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteTracks, "Convert FLiteTracks to FTracks", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{o2::framework::adaptAnalysisTask<FemtoProducerLiteConverter>(context)};
  return workflow;
}
