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
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"},
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
                     HistType::kTH2F, {amplitudeAxisLarge, modAxis});
    mHistManager.add("cluETime", "Cluster time vs E",
                     HistType::kTH3F, {amplitudeAxisLarge, timeAxisLarge, modAxis});
    mHistManager.add("mggRe", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggMi", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggReCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggMiCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggReDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});
    mHistManager.add("mggMiDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, centAxis});

    mHistManager.add("cluOcc", "Cluster occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluCPVOcc", "Cluster with CPV occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluDispOcc", "Cluster with Disp occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluBothOcc", "Cluster with Both occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluE", "Cluster energy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
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

      if (clu.e() < mMinCluE || clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime) {
        continue;
      }

      mHistManager.fill(HIST("cluSp"), clu.e(), clu.mod());
      if (clu.e() > mOccE) {
        mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
        if (clu.trackdist() > 2.) {
          mHistManager.fill(HIST("cluCPVOcc"), clu.x(), clu.z(), clu.mod());
          if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
            mHistManager.fill(HIST("cluBothOcc"), clu.x(), clu.z(), clu.mod());
          }
        }
        if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
          mHistManager.fill(HIST("cluDispOcc"), clu.x(), clu.z(), clu.mod());
        }
        mHistManager.fill(HIST("cluE"), clu.x(), clu.z(), clu.mod(), clu.e());
      }

      float cen1 = 1.; // TODO: To be extended to use centrality
      // inv mass
      auto clu2 = clu;
      ++clu2;
      int nMix = mMixedEvents; // Number of events to mix
      uint64_t bcurrent = 0;
      for (; clu2 != clusters.end() && nMix > 0; clu2++) {
        if (clu2.e() < mMinCluE || clu2.ncell() < mMinCluNcell ||
            clu2.time() > mMaxCluTime || clu2.time() < mMinCluTime) {
          continue;
        }
        double m = pow(clu.e() + clu2.e(), 2) - pow(clu.px() + clu2.px(), 2) -
                   pow(clu.py() + clu2.py(), 2) - pow(clu.pz() + clu2.pz(), 2);
        if (m > 0) {
          m = sqrt(m);
        }
        double pt = sqrt(pow(clu.px() + clu2.px(), 2) +
                         pow(clu.py() + clu2.py(), 2));
        if (clu.collision() == clu2.collision()) { // Real
          mHistManager.fill(HIST("mggRe"), m, pt, cen1);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggReCPV"), m, pt, cen1);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggReDisp"), m, pt, cen1);
          }
        } else { // Mixed
          if (clu2.collision().bc().globalBC() != bcurrent) {
            --nMix;
            bcurrent = clu2.collision().bc().globalBC();
          }
          mHistManager.fill(HIST("mggMi"), m, pt, cen1);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggMiCPV"), m, pt, cen1);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggMiDisp"), m, pt, cen1);
          }
        }
      }
    }
  }

  //_____________________________________________________________________________
  bool TestLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    // Parameterizatino for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosPi0>(cfgc)};
}
