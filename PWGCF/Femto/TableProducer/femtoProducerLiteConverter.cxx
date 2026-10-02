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
/// \brief Task that converts FLite* tables (binned kinematics) back to their full-precision counterparts
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerLiteConverter {
  o2::framework::Produces<o2::aod::FCols> producedCols;
  o2::framework::Produces<o2::aod::FColShapes> producedColShapes;
  o2::framework::Produces<o2::aod::FTracks> producedTracks;
  o2::framework::Produces<o2::aod::FLambdas> producedLambdas;
  o2::framework::Produces<o2::aod::FK0shorts> producedK0shorts;
  o2::framework::Produces<o2::aod::FXis> producedXis;
  o2::framework::Produces<o2::aod::FOmegas> producedOmegas;
  o2::framework::Produces<o2::aod::FSigmas> producedSigmas;
  o2::framework::Produces<o2::aod::FSigmaPlus> producedSigmaPlus;

  void init(o2::framework::InitContext&)
  {
    if (!doprocessLiteCols) {
      LOG(fatal) << "At least processLiteCols needs to be enabled, otherwise the collision index points no where!";
    }
  }

  void processLiteCols(o2::aod::FLiteCols::iterator const& liteCol)
  {
    producedCols(liteCol.posZ(),
                 liteCol.mult(),
                 liteCol.cent(),
                 liteCol.magField());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteCols, "Convert FLiteCols to FCols", true);

  void processLiteColShapes(o2::aod::FLiteColShapes::iterator const& liteColShape)
  {
    producedColShapes(liteColShape.qvec(),
                      liteColShape.eventPlaneAngle());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteColShapes, "Convert FLiteColShapes to FColShapes", false);

  void processLiteTracks(o2::aod::FLiteTracks::iterator const& liteTrack)
  {
    producedTracks(liteTrack.fLiteColId(),
                   liteTrack.signedPt(),
                   liteTrack.eta(),
                   liteTrack.phi());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteTracks, "Convert FLiteTracks to FTracks", true);

  void processLiteLambdas(o2::aod::FLiteLambdas::iterator const& liteLambda)
  {
    producedLambdas(liteLambda.fLiteColId(),
                    liteLambda.signedPt(),
                    liteLambda.eta(),
                    liteLambda.phi(),
                    liteLambda.lambdaMass(),
                    liteLambda.posDauId(),
                    liteLambda.negDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteLambdas, "Convert FLiteLambdas to FLambdas", false);

  void processLiteK0shorts(o2::aod::FLiteK0shorts::iterator const& liteK0short)
  {
    producedK0shorts(liteK0short.fLiteColId(),
                     liteK0short.pt(),
                     liteK0short.eta(),
                     liteK0short.phi(),
                     liteK0short.k0shortMass(),
                     liteK0short.posDauId(),
                     liteK0short.negDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteK0shorts, "Convert FLiteK0shorts to FK0shorts", false);

  void processLiteXis(o2::aod::FLiteXis::iterator const& liteXi)
  {
    producedXis(liteXi.fLiteColId(),
                liteXi.signedPt(),
                liteXi.eta(),
                liteXi.phi(),
                liteXi.xiMass(),
                liteXi.bachelorId(),
                liteXi.posDauId(),
                liteXi.negDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteXis, "Convert FLiteXis to FXis", false);

  void processLiteOmegas(o2::aod::FLiteOmegas::iterator const& liteOmega)
  {
    producedOmegas(liteOmega.fLiteColId(),
                   liteOmega.signedPt(),
                   liteOmega.eta(),
                   liteOmega.phi(),
                   liteOmega.omegaMass(),
                   liteOmega.bachelorId(),
                   liteOmega.posDauId(),
                   liteOmega.negDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteOmegas, "Convert FLiteOmegas to FOmegas", false);

  void processLiteSigmas(o2::aod::FLiteSigmas::iterator const& liteSigma)
  {
    producedSigmas(liteSigma.fLiteColId(),
                   liteSigma.signedPt(),
                   liteSigma.eta(),
                   liteSigma.phi(),
                   liteSigma.sigmaMass(),
                   liteSigma.chaDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteSigmas, "Convert FLiteSigmas to FSigmas", false);

  void processLiteSigmaPlus(o2::aod::FLiteSigmaPlus::iterator const& liteSigmaPlus)
  {
    producedSigmaPlus(liteSigmaPlus.fLiteColId(),
                      liteSigmaPlus.signedPt(),
                      liteSigmaPlus.eta(),
                      liteSigmaPlus.phi(),
                      liteSigmaPlus.sigmaMass(),
                      liteSigmaPlus.chaDauId());
  }
  PROCESS_SWITCH(FemtoProducerLiteConverter, processLiteSigmaPlus, "Convert FLiteSigmaPlus to FSigmaPlus", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{o2::framework::adaptAnalysisTask<FemtoProducerLiteConverter>(context)};
  return workflow;
}
