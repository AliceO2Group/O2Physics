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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"

/// \struct PHOS pi0 analysis
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosPi0 {

  Configurable<float> mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -50.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"mMinCluTime", 100.e-9, "Max. cluster time"};
  Configurable<int> mMinCluNcell{"mMinCluNcell", 2, "min cells in cluster"};
  Configurable<int> mMixedEvents{"mixedEvents", 10, "number of events to mix"};

  Configurable<float> mOccE{"minOccE", 0.6, "Minimum cluster energy to fill occupancy"};

  using FilteredClusters = soa::Filtered<aod::CaloClusters>;

  HistogramRegistry mHistManager{"phosPi0Histograms"};

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS Cluster QA monitor task ...";

    const AxisSpec
      cluPhiAxis{256, 4.3633231, 5.5850536, "phi", ""},
      cluZAxis{56, -60., 60., "z", ""},
      amplitudeAxisLarge{100, 0., 20., "amplitude", "Amplutude (GeV)"},
      timeAxisLarge{500, -1500.e-9, 3500.e-9, "celltime", "cell time (ns)"},
      multAxis{100, 0., 100.},
      mggAxis{250, 0., 1., "mgg", "m_{#gamma#gamma} (GeV/c^{2})"},
      vertexAxis{100, -20., 20., "z", "z (cm)"},
      centAxis{10, 0., 10.},
      centralityAxis{100, 0., 100., "centrality", "centrality"};

    mHistManager.add("eventsAll", "Number of events", HistType::kTH1F, {{1, 0., 1.}});
    mHistManager.add("contributors", "Centrality", HistType::kTH1F, {centralityAxis});
    mHistManager.add("vertex", "vertex", HistType::kTH1F, {vertexAxis});
    mHistManager.add("vertex2", "PHOS with vertex", HistType::kTH1F, {vertexAxis});

    mHistManager.add("cluSp", "Cluster spectrum per module",
                     HistType::kTH2F, {amplitudeAxisLarge, {4, 1., 5., "module", "module"}});
    mHistManager.add("cluETime", "Cluster time vs E",
                     HistType::kTH3F, {amplitudeAxisLarge, timeAxisLarge, {4, 1., 5., "module", "module"}});
    mHistManager.add("mggRe", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggMi", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});

    mHistManager.add("cluOcc", "Cluster occupancy ", HistType::kTH2F, {cluPhiAxis, cluZAxis});
    mHistManager.add("cluE", "Cluster energy", HistType::kTH2F, {cluPhiAxis, cluZAxis});
    mHistManager.add("cluTime", "Cluster time", HistType::kTH2F, {cluPhiAxis, cluZAxis});
  }

  /// \brief Process PHOS data
  void process(aod::BCs const& bcs,
               aod::Collisions const& collisions,
               aod::CaloClusters const& clusters)
  {
    for (const auto& col : collisions) {
      mHistManager.fill(HIST("contributors"), col.numContrib());
      mHistManager.fill(HIST("vertex"), col.posZ());
      //   mHistManager.fill(HIST("centralityFT0M"), collision.centFT0M);
      //   mHistManager.fill(HIST("centralityFDDM"), collision.centFDDM);
      //   mHistManager.fill(HIST("centralityNTPV"), collision.centNTPV);
    }
    uint64_t bcevent = 0;
    for (const auto& clu : clusters) {
      if (clu.collision().bc().globalBC() != bcevent) {
        mHistManager.fill(HIST("vertex2"), clu.collision().posZ());
        mHistManager.fill(HIST("eventsAll"), 0.);
        bcevent = clu.collision().bc().globalBC();
      }

      mHistManager.fill(HIST("cluETime"), clu.e(), clu.time(), clu.mod());

      if (clu.caloType() != 0 || clu.e() < mMinCluE || clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime) {
        continue;
      }

      float cen1 = log(1. + clu.collision_as<o2::aod::Collisions>().numContrib()) / log(2.);
      int iCenBin1 = static_cast<int>(cen1 / 10);
      if (iCenBin1 < 0 || iCenBin1 > 9) {
        continue;
      }

      mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
      if (clu.e() > mOccE) {
        double phi = 6.2831853 + atan2(clu.y(), clu.x()); // Only negative phi from tan2
        mHistManager.fill(HIST("cluOcc"), phi, clu.z());
        mHistManager.fill(HIST("cluE"), phi, clu.z(), clu.e());
        mHistManager.fill(HIST("cluTime"), phi, clu.z(), clu.time());
      }

      // inv mass
      auto clu2 = clu;
      ++clu2;
      int nMix = mMixedEvents; // Number of events to mix
      uint64_t bcurrent = 0;
      for (; clu2 != clusters.end() && nMix > 0; clu2++) {
        if (clu2.caloType() != 0 || clu2.e() < mMinCluE || clu2.ncell() < mMinCluNcell ||
            clu2.time() > mMaxCluTime || clu2.time() < mMinCluTime) {
          continue;
        }
        float cen2 = log(1. + clu2.collision_as<o2::aod::Collisions>().numContrib()) / log(2.);
        int iCenBin2 = static_cast<int>(cen2 / 10);
        if (iCenBin1 != iCenBin2) {
          continue;
        }

        double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                   pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
        if (m > 0) {
          m = sqrt(m);
        }
        double pt = sqrt(pow(clu.px() + clu2.px(), 2) +
                         pow(clu.py() + clu2.py(), 2));
        if (clu.collision().bc().globalBC() == clu2.collision().bc().globalBC()) { // Real
          mHistManager.fill(HIST("mggRe"), m, pt, cen1);
        } else { // Mixed
          if (clu2.collision().bc().globalBC() != bcurrent) {
            --nMix;
            bcurrent = clu2.collision().bc().globalBC();
          }
          mHistManager.fill(HIST("mggMi"), m, pt, cen1);
        }
      }
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosPi0>(cfgc)};
}
