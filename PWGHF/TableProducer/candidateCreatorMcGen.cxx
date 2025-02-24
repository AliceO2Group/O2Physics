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

/// \file candidateCreatorMcGen.cxx
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

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/Utils/utilsMcGen.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::framework;

/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreatorMcGen {

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

  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  void process(aod::McCollisions const& mcCollisions,
               aod::McParticles const& mcParticles)
  {

    for (const auto& mcCollision : mcCollisions) {
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      if (fill2Prong) {
        hf_mc_gen::fillMcMatchGen2Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen2Prong, rejectBackground2Prong);
      }
      if (fill3Prong) {
        hf_mc_gen::fillMcMatchGen3Prong(mcParticles, mcParticlesPerMcColl, rowMcMatchGen3Prong, rejectBackground3Prong);
      }
    }
    if (fillBplus) {
      hf_mc_gen::fillMcMatchGenBplus(mcParticles, rowMcMatchGenBplus);
    }
    if (fillB0) {
      hf_mc_gen::fillMcMatchGenB0(mcParticles, rowMcMatchGenB0);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfCandidateCreatorMcGen>(cfgc)};
}
