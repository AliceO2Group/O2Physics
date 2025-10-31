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

/// \file converterDstarIndices.cxx
/// \brief Task for conversion of HfDstars to version 001, using the collision index from the D0 daughter
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

struct HfConverterDstarIndices {
  Produces<aod::HfDstars_001> dstarIndices;

  void process(aod::HfDstars_000::iterator const& candDstar,
               aod::Hf2Prongs const&)
  {
    auto candDzero = candDstar.prongD0_as<aod::Hf2Prongs>();
    dstarIndices(candDzero.collisionId(), candDstar.prong0Id(), candDstar.prongD0Id());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfConverterDstarIndices>(cfgc),
  };
}
