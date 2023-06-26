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
// O2 includes

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/CaloClusters.h"
#include "DataFormatsPHOS/TriggerRecord.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "../filterTables.h"
#include "Framework/HistogramRegistry.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0Recalculation>;
using MyV0Photon = MyV0Photons::iterator;

using MyCollisions = soa::Join<aod::Collisions, aod::PhotonFilters>;
using MyCollision = MyCollisions::iterator;

struct EMPhotonFilterQC {
  HistogramRegistry registry{"registry"};

  void init(o2::framework::InitContext&)
  {
    addhistograms();
  }

  void addhistograms()
  {
    auto hEventCounter = registry.add<TH1>("hEventCounter", "hEventCounter", kTH1F, {{6, 0.5f, 6.5f}});
    hEventCounter->GetXaxis()->SetBinLabel(1, "all");
    hEventCounter->GetXaxis()->SetBinLabel(2, "PCM W wire");
    registry.add<TH2>("hGammaConvXY", "Conversion point XY;X (cm);Y (cm)", kTH2F, {{200, -100.f, 100.f}, {200, -100.f, 100.f}});
    registry.add<TH2>("hGammaConvRZ", "Conversion point RZ;Z (cm);R_{xy} (cm)", kTH2F, {{200, -100.f, 100.f}, {100, 0.f, 100.f}});
  }

  Preslice<MyV0Photons> perCollision_pcm = aod::v0photon::collisionId;
  void processPCM(MyCollisions const& collisions, MyV0Photons const& v0photons, aod::V0Legs const& v0legs)
  {
    for (auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 1);

      if (!collision.hasPCMWwire()) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2);

      auto v0photons_coll = v0photons.sliceBy(perCollision_pcm, collision.globalIndex());
      for (auto& v0 : v0photons_coll) {
        registry.fill(HIST("hGammaConvXY"), v0.recalculatedVtxX(), v0.recalculatedVtxY());
        registry.fill(HIST("hGammaConvRZ"), v0.recalculatedVtxZ(), v0.recalculatedVtxR());
      }
    } // end of collision loop
  }

  void processPHOS(MyCollisions const& collisions) {}
  void processEMC(MyCollisions const& collisions) {}

  PROCESS_SWITCH(EMPhotonFilterQC, processPCM, "Process PCM software trigger QC", true);
  PROCESS_SWITCH(EMPhotonFilterQC, processPHOS, "Process PHOS software trigger QC", false);
  PROCESS_SWITCH(EMPhotonFilterQC, processEMC, "Process EMC software trigger QC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{
    adaptAnalysisTask<EMPhotonFilterQC>(cfg, TaskName{"em-photon-filter-qc"})};
}
