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

/// \brief skim PHOS cluster info
/// dependencies: o2-analysis-calo-cluster
/// \author daiki.sekihata@cern.ch

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/CaloClusters.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct skimmerPHOS {
  Produces<aod::PHOSClusters> tablePHOS;

  // Configurable for filter/cuts
  Configurable<float> minM20{"minM20", 0.1, "Minimum M20 cut"};
  Configurable<float> minE{"minE", 0.2, "Minimum energy cut in GeV"};
  Configurable<int> minNcell{"minNcell", 0, "Minimum Ncell cut"};

  HistogramRegistry registry;
  void init(o2::framework::InitContext&)
  {
    auto hPHOSClusterFilter = registry.add<TH1>("hPHOSClusterFilter", "hPHOSClusterFilter", kTH1I, {{5, 0.5, 5.5}});
    hPHOSClusterFilter->GetXaxis()->SetBinLabel(1, "in");
    hPHOSClusterFilter->GetXaxis()->SetBinLabel(2, "energy cut");
    hPHOSClusterFilter->GetXaxis()->SetBinLabel(3, "M20 cut");
    hPHOSClusterFilter->GetXaxis()->SetBinLabel(4, "Ncell cut");
    hPHOSClusterFilter->GetXaxis()->SetBinLabel(5, "out");
  }

  void process(aod::Collisions const&, aod::BCs const&, aod::CaloClusters const& clusters)
  {
    for (auto& cluster : clusters) {
      if (cluster.caloType() != 0) { // 0 is PHOS, 1 is EMCal
        continue;
      }

      registry.fill(HIST("hPHOSClusterFilter"), 1);

      if (cluster.e() < minE) {
        continue;
        registry.fill(HIST("hPHOSClusterFilter"), 2);
      }
      if (cluster.m20() < minM20) {
        continue;
        registry.fill(HIST("hPHOSClusterFilter"), 3);
      }
      if (cluster.ncell() < minNcell) {
        continue;
        registry.fill(HIST("hPHOSClusterFilter"), 4);
      }

      registry.fill(HIST("hPHOSClusterFilter"), 5);
      tablePHOS(
        cluster.colId(), cluster.bcId(), cluster.trackIndex(),
        cluster.e(), cluster.globalx(), cluster.globaly(), cluster.globalz(),
        cluster.m02(), cluster.m20(), cluster.ncell(), cluster.time(), cluster.distBad(), cluster.nlm());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<skimmerPHOS>(cfgc, TaskName{"skimmer-phos"})};
  return workflow;
}
