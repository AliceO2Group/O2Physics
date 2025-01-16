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

/// \file candidateCreatorMcGenOnly.cxx
/// \brief McGen only selection of heavy-flavour particles
///
/// \author Nima Zardoshti, nima.zardoshti@cern.ch, CERN

#include <memory>
#include <string>
#include <vector>

#include <TPDGCode.h>

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"

#include "PWGHF/Core/HfMcGenHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreatorMcGenOnly {

  Produces<aod::HfCand2ProngMcGen> rowMcMatchGen2Prong;
  Produces<aod::HfCand3ProngMcGen> rowMcMatchGen3Prong;
  Produces<aod::HfCandBplusMcGen> rowMcMatchGenBplus;
  Produces<aod::HfCandB0McGen> rowMcMatchGenB0;
  Configurable<bool> fill2Prong{"fill2Prong", false, "fill table for 2 prong candidates"};
  Configurable<bool> fill3Prong{"fill3Prong", false, "fill table for 3 prong candidates"};
  Configurable<bool> fillBplus{"fillBplus", false, "fill table for for B+ candidates"};
  Configurable<bool> fillB0{"fillB0", false, "fill table for B0 candidates"};
  Configurable<bool> rejectBackground2Prong{"rejectBackground2Prong", false, "Reject particles from PbPb background for 2 prong candidates"};
  Configurable<bool> rejectBackground3Prong{"rejectBackground3Prong", false, "Reject particles from PbPb background for 3 prong candidates"};
  Configurable<bool> rejectBackgroundBplus{"rejectBackgroundBplus", false, "Reject particles from PbPb background for B+ candidates"};
  Configurable<bool> rejectBackgroundB0{"rejectBackgroundB0", false, "Reject particles from PbPb background for B0 candidates"};
  Configurable<bool> createDplus{"createDplus", false, "Create D+ in 3 prong"};
  Configurable<bool> createDs{"createDs", false, "Create Ds in 3 prong"};
  Configurable<bool> createLc{"createLc", false, "Create Lc in 3 prong"};
  Configurable<bool> createXic{"createXic", false, "Create Xic in 3 prong"};
  Configurable<float> minPosZ{"minPosZ", 10.0, "min z position of promary vertex"};

  // Filter filterSelectMcCollisions = nabs(aod::mcCollision::posZ) < minPosZ;

  void process(aod::McCollision const&,
               aod::McParticles const& mcParticles)
  {
    if (fill2Prong) {
      hf_mcgen_helper::fill2ProngMcMatchGen(mcParticles, rowMcMatchGen2Prong, rejectBackground2Prong);
    }
    if (fill3Prong) {
      hf_mcgen_helper::fill3ProngMcMatchGen(mcParticles, rowMcMatchGen3Prong, rejectBackground3Prong, createDplus, createDs, createLc, createXic);
    }
    if (fillBplus) {
      hf_mcgen_helper::fillBplusMcMatchGen(mcParticles, rowMcMatchGenBplus);
    }
    if (fillB0) {
      hf_mcgen_helper::fillB0McMatchGen(mcParticles, rowMcMatchGenB0);
    }
  }
  PROCESS_SWITCH(HfCandidateCreatorMcGenOnly, process, "Process candidates", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorMcGenOnly>(cfgc, TaskName{"hf-candidate-creator-mcgen-only"})};
}
