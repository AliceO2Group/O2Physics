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

/// \file converterReduced3ProngsMl.cxx
/// \brief Task for conversion of HfRed3ProngsMl to version 001
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/runDataProcessing.h>

using namespace o2;
using namespace o2::framework;

// Swaps covariance matrix elements if the data is known to be bogus (collision_000 is bogus)
struct HfConverterReduced3ProngsMl {
  Produces<aod::HfRed3ProngsMl_001> ml3Prongs;

  void process(aod::HfRed3ProngsMl_000::iterator const& mlScoreTable)
  {
    ml3Prongs(mlScoreTable.mlScoreBkgMassHypo0(), mlScoreTable.mlScorePromptMassHypo0(), mlScoreTable.mlScoreNonpromptMassHypo0(), -1.f, -1.f, -1.f);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfConverterReduced3ProngsMl>(cfgc),
  };
}
