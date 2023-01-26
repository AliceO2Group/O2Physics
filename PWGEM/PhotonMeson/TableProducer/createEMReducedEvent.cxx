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
// ========================
//
// This code produces reduced events for photon analyses.
//    Please write to: daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
//#include "Common/Core/trackUtilities.h"
//#include "Common/Core/TrackSelection.h"
//#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
//#include "Common/DataModel/PIDResponse.h"
//#include "Common/Core/RecoDecay.h"
//#include "DetectorsVertexing/DCAFitterN.h"
//#include "DetectorsBase/Propagator.h"
//#include "DetectorsBase/GeometryManager.h"
//#include "DataFormatsParameters/GRPObject.h"
//#include "DataFormatsParameters/GRPMagField.h"
//#include "CCDB/BasicCCDBManager.h"
//#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

struct createEMReducedEvent {
  Produces<o2::aod::EMReducedEvents> event;

  enum SubSystem {
    kPCM = 0x1,
    kPHOS = 0x2,
    kEMC = 0x4,
    kUndef = -1,
  };

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{5, 0.5f, 5.5f}}}},
    },
  };

  Preslice<aod::V0Photons> perCollision_pcm = aod::v0photon::collisionId;
  Preslice<aod::PHOSClusters> perCollision_phos = aod::phoscluster::collisionId;

  template <uint8_t system, typename TPCMs, typename TPHOSs, typename TEMCs>
  void process(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const&, TPCMs const& v0photons, TPHOSs const& phosclusters, TEMCs const& emcclusters)
  // void process(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const&, aod::V0Photons const& v0photons)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      int ng_pcm = 0;
      int ng_phos = 0;
      int ng_emc = 0;

      if constexpr (static_cast<bool>(system & kPCM)) {
        auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
        ng_pcm = v0photons_coll.size();
      }
      if constexpr (static_cast<bool>(system & kPHOS)) {
        auto phos_coll = phosclusters.sliceBy(perCollision_phos, collision.globalIndex());
        ng_phos = phos_coll.size();
      }

      uint64_t tag = 0;
      // store event selection decisions
      for (int i = 0; i < kNsel; i++) {
        if (collision.selection()[i] > 0) {
          tag |= (uint64_t(1) << i);
        }
      }
      event(collision.globalIndex(), tag, collision.bc().runNumber(), collision.sel8(),
            collision.posX(), collision.posY(), collision.posZ(),
            collision.numContrib(), collision.collisionTime(), collision.collisionTimeRes(),
            ng_pcm, ng_phos, ng_emc); // ng is needed for event mixing to filter events that contain at least 1 photon.

    } // end of collision loop

  } // end of process

  void process_PCM(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const& bcs, aod::V0Photons const& v0photons)
  {
    const uint8_t sysflag = kPCM;
    ;
    process<sysflag>(collisions, bcs, v0photons, nullptr, nullptr);
  }
  void process_PHOS(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const& bcs, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPHOS;
    process<sysflag>(collisions, bcs, nullptr, phosclusters, nullptr);
  }
  void process_PCM_PHOS(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::BCs const& bcs, aod::V0Photons const& v0photons, aod::PHOSClusters const& phosclusters)
  {
    const uint8_t sysflag = kPCM | kPHOS;
    process<sysflag>(collisions, bcs, v0photons, phosclusters, nullptr);
  }

  void processDummy(soa::Join<aod::Collisions, aod::EvSels> const& collisions) {}

  PROCESS_SWITCH(createEMReducedEvent, process_PCM, "create em event table for PCM", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PHOS, "create em event table for PHOS", false);
  PROCESS_SWITCH(createEMReducedEvent, process_PCM_PHOS, "create em event table for PCM, PHOS", false);
  PROCESS_SWITCH(createEMReducedEvent, processDummy, "processDummy", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<createEMReducedEvent>(cfgc, TaskName{"create-emreduced-event"}),
  };
}
