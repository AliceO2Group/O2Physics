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

/// \file femtoProducerMaskConverter001To002.cxx
/// \brief Task that converts *Masks_001 tables to *Masks_002
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#include "PWGCF/Femto/Core/dataTypes.h"
#include "PWGCF/Femto/DataModel/FemtoTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

using namespace o2::analysis::femto;

struct FemtoProducerMaskConverter001To002 {
  o2::framework::Produces<o2::aod::FLambdaMasks_002> producedLambdaMasks;
  o2::framework::Produces<o2::aod::FK0shortMasks_002> producedK0shortMasks;
  o2::framework::Produces<o2::aod::FXiMasks_002> producedXiMasks;
  o2::framework::Produces<o2::aod::FOmegaMasks_002> producedOmegaMasks;

  void init(o2::framework::InitContext&)
  {
  }

  void processLambdaMasks001(o2::aod::FLambdaMasks_001::iterator const& lambdaMaskOld)
  {
    producedLambdaMasks(static_cast<datatypes::V0MaskType>(lambdaMaskOld.mask()));
  }
  PROCESS_SWITCH(FemtoProducerMaskConverter001To002, processLambdaMasks001, "Convert LambdaMasks_001 to LambdaMasks_002", false);

  void processK0shortMasks001(o2::aod::FK0shortMasks_001::iterator const& k0shortOld)
  {
    producedK0shortMasks(static_cast<datatypes::V0MaskType>(k0shortOld.mask()));
  }
  PROCESS_SWITCH(FemtoProducerMaskConverter001To002, processK0shortMasks001, "Convert K0shortMasks_001 to K0shortMasks_002", false);

  void processXiMasks001(o2::aod::FXiMasks_001::iterator const& xiMaskOld)
  {
    producedXiMasks(static_cast<datatypes::CascadeMaskType>(xiMaskOld.mask()));
  }
  PROCESS_SWITCH(FemtoProducerMaskConverter001To002, processXiMasks001, "Convert XiMasks_001 to XiMasks_002", false);

  void processOmegaMasks001(o2::aod::FOmegaMasks_001::iterator const& omegaMaskOld)
  {
    producedOmegaMasks(static_cast<datatypes::CascadeMaskType>(omegaMaskOld.mask()));
  }
  PROCESS_SWITCH(FemtoProducerMaskConverter001To002, processOmegaMasks001, "Convert OmegaMasks_001 to OmegaMasks_002", false);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& context)
{
  o2::framework::WorkflowSpec workflow{o2::framework::adaptAnalysisTask<FemtoProducerMaskConverter001To002>(context)};
  return workflow;
}
