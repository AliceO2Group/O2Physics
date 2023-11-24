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
//

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct lambdakzerocollisionbuilder {
  // TPC pid (copied over from central services for reference)
  Produces<aod::V0CollRefs> v0collref; // raw table for checks
  Produces<aod::V0Collision> v0coll;   // table with Nsigmas

  Configurable<bool> fillEmptyCollisions{"fillEmptyCollisions", false, "fill collision entries without candidates"};

  // For manual sliceBy
  Preslice<aod::V0Datas> perCollision = o2::aod::v0data::collisionId;

  void init(InitContext& context)
  {
  }

  void process(soa::Join<aod::Collisions, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As> const& collisions, aod::V0Datas const& V0s)
  {
    int currentCollIdx = -1;
    for (const auto& collision : collisions) {
      const uint64_t collIdx = collision.globalIndex();
      auto V0Table_thisCollision = V0s.sliceBy(perCollision, collIdx);
      // V0 table sliced
      if (V0Table_thisCollision.size() > 0 || fillEmptyCollisions) {
        if (currentCollIdx != collIdx) {
          v0coll(collision.posX(), collision.posY(), collision.posZ(),
                 collision.centFT0M(), collision.centFT0A(),
                 collision.centFT0C(), collision.centFV0A());
          currentCollIdx = collIdx;
        }
      }
      for (int i = 0; i < V0Table_thisCollision.size(); i++) {
        v0collref(v0coll.lastIndex());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdakzerocollisionbuilder>(cfgc)};
}
