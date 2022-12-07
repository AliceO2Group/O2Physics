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

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"

/// \struct PHOS QA
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since March, 2022
///
/// This task monitors simplу cluster quantities like
/// - Energy distribution
/// - Time distribution
/// - Count rate in 2D representation
struct phosCluQA {

  o2::framework::Configurable<double> mMinCluE{"minCluE", 0.3, "Minimum cluster energy for histograms."};
  //   o2::framework::Configurable<double> mMinCellTimeMain{"minCellTimeMain", -50, "Min. cell time of main bunch selection"};
  //   o2::framework::Configurable<double> mMaxCellTimeMain{"maxCellTimeMain", 100, "Max. cell time of main bunch selection"};
  o2::framework::Configurable<int> mVetoBCID{"vetoBCID", -1, "BC ID to be excluded"};

  o2::framework::HistogramRegistry mHistManager{"phosCluQAHistograms"};

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    LOG(info) << "Initializing PHOS Cluster QA monitor task ...";

    const o2Axis
      cluPhiAxis{256, 4.3633231, 5.5850536, "phi", ""},
      cluZAxis{56, -60., 60., "z", ""},
      amplitudeAxisLarge{1000, 0., 10., "amplitude", "Amplutude (GeV)"},
      timeAxisLarge{1000, -1500.e-9, 3500.e-9, "celltime", "cell time (ns)"},
      multAxis{100, 0., 100.},
      mggAxis{250, 0., 1.},
      bcAxis{3501, -0.5, 3500.5};
    mHistManager.add("eventsAll", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventsSelected", "Number of events", o2HistType::kTH1F, {{1, 0.5, 1.5}});
    mHistManager.add("eventBCAll", "Bunch crossing ID of event (all events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("eventBCSelected", "Bunch crossing ID of event (selected events)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCAll", "Bunch crossing ID of cell (all cells)", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCSelected", "Bunch crossing ID of cell (selected cells)", o2HistType::kTH1F, {{bcAxis}});

    for (int i = 1; i < 5; ++i) {
      mHistManager.add(Form("cluSpM%d", i), Form("Cluster spectrum for module %d", i),
                       o2HistType::kTH1F, {amplitudeAxisLarge});
      mHistManager.add(Form("cluMultM%d", i), Form("Cluster multiplicity for module %d", i),
                       o2HistType::kTH2F, {amplitudeAxisLarge, multAxis});
      mHistManager.add(Form("cluETimeM%d", i), Form("Correlation between cell amplitude and time in module %d", i),
                       o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge});
    }
    mHistManager.add("cluOcc", "Cluster occupancy ", o2HistType::kTH2F, {cluPhiAxis, cluZAxis});
    mHistManager.add("cluE", "Cluster energy", o2HistType::kTH2F, {cluPhiAxis, cluZAxis});
    mHistManager.add("cluTime", "Cluster time", o2HistType::kTH2F, {cluPhiAxis, cluZAxis});
    for (int i = 1; i < 5; ++i) {
      for (int j = i; j < 5; ++j) {
        mHistManager.add(Form("mggReM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                         o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge});
        mHistManager.add(Form("mggMiM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                         o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge});
      }
    }
  }

  /// \brief Process PHOS data
  void process(o2::aod::BCs const& bcs,
               o2::aod::Collisions const&,
               o2::aod::CaloClusters const& clusters)
  {
    for (const auto& bc : bcs) {
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(bc.globalBC());
      mHistManager.fill(HIST("eventsAll"), 1);
      mHistManager.fill(HIST("eventBCAll"), eventIR.bc);
      if (mVetoBCID >= 0 && eventIR.bc == mVetoBCID)
        continue;
      mHistManager.fill(HIST("eventsSelected"), 1);
      mHistManager.fill(HIST("eventBCSelected"), eventIR.bc);
    }
    for (const auto& clu : clusters) {
      if (clu.caloType() != 0)
        continue;
      clu.collision_as<o2::aod::Collisions>();
      clu.collision_as<o2::aod::Collisions>().bc_as<o2::aod::BCs>();
      o2::InteractionRecord ir;
      ir.setFromLong(clu.collision_as<o2::aod::Collisions>().bc_as<o2::aod::BCs>().globalBC());
      mHistManager.fill(HIST("cluBCAll"), ir.bc);

      if (mVetoBCID >= 0 && ir.bc == mVetoBCID)
        continue;
      mHistManager.fill(HIST("cluBCSelected"), ir.bc);
      if (clu.mod() == 1) {
        mHistManager.fill(HIST("cluSpM1"), clu.e());
        mHistManager.fill(HIST("cluMultM1"), clu.e(), clu.ncell());
        mHistManager.fill(HIST("cluETimeM1"), clu.time(), clu.e());
      }
      if (clu.mod() == 2) {
        mHistManager.fill(HIST("cluSpM2"), clu.e());
        mHistManager.fill(HIST("cluMultM2"), clu.e(), clu.ncell());
        mHistManager.fill(HIST("cluETimeM2"), clu.time(), clu.e());
      }
      if (clu.mod() == 3) {
        mHistManager.fill(HIST("cluSpM3"), clu.e());
        mHistManager.fill(HIST("cluMultM3"), clu.e(), clu.ncell());
        mHistManager.fill(HIST("cluETimeM3"), clu.time(), clu.e());
      }
      if (clu.mod() == 4) {
        mHistManager.fill(HIST("cluSpM4"), clu.e());
        mHistManager.fill(HIST("cluMultM4"), clu.e(), clu.ncell());
        mHistManager.fill(HIST("cluETimeM4"), clu.time(), clu.e());
      }
      if (clu.e() > mMinCluE) {
        double phi = 6.2831853 + atan2(clu.y(), clu.x()); // Only negative phi from tan2
        mHistManager.fill(HIST("cluOcc"), phi, clu.z());
        mHistManager.fill(HIST("cluE"), phi, clu.z(), clu.e());
        mHistManager.fill(HIST("cluTime"), phi, clu.z(), clu.time());
      }
    }
    // inv mass
    for (const auto& clu1 : clusters) {
      auto clu2 = clu1;
      ++clu2;
      int nMix = 100; // Number of photons to mix
      for (; clu2 != clusters.end() && nMix > 0; clu2++) {
        double m = pow(clu1.e() + clu2.e(), 2) - pow(clu1.px() + clu2.px(), 2) -
                   pow(clu1.py() + clu2.py(), 2) - pow(clu1.pz() + clu2.pz(), 2);
        if (m > 0)
          m = sqrt(m);
        double pt = sqrt(pow(clu1.px() + clu2.px(), 2) +
                         pow(clu1.py() + clu2.py(), 2));
        if (clu1.collision().bc().globalBC() == clu2.collision().bc().globalBC()) { // Real
          if (clu1.mod() == 1 && clu2.mod() == 1) {
            mHistManager.fill(HIST("mggReM11"), m, pt);
          }
          if (clu1.mod() == 2 && clu2.mod() == 2) {
            mHistManager.fill(HIST("mggReM22"), m, pt);
          }
          if (clu1.mod() == 3 && clu2.mod() == 3) {
            mHistManager.fill(HIST("mggReM33"), m, pt);
          }
          if (clu1.mod() == 4 && clu2.mod() == 4) {
            mHistManager.fill(HIST("mggReM44"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 2) || (clu1.mod() == 2 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggReM12"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 3) || (clu1.mod() == 3 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggReM13"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggReM14"), m, pt);
          }
          if ((clu1.mod() == 2 && clu2.mod() == 3) || (clu1.mod() == 3 && clu2.mod() == 2)) {
            mHistManager.fill(HIST("mggReM23"), m, pt);
          }
          if ((clu1.mod() == 2 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 2)) {
            mHistManager.fill(HIST("mggReM24"), m, pt);
          }
          if ((clu1.mod() == 3 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 3)) {
            mHistManager.fill(HIST("mggReM34"), m, pt);
          }
        } else { // Mixed
          --nMix;
          if (clu1.mod() == 1 && clu2.mod() == 1) {
            mHistManager.fill(HIST("mggMiM11"), m, pt);
          }
          if (clu1.mod() == 2 && clu2.mod() == 2) {
            mHistManager.fill(HIST("mggMiM22"), m, pt);
          }
          if (clu1.mod() == 3 && clu2.mod() == 3) {
            mHistManager.fill(HIST("mggMiM33"), m, pt);
          }
          if (clu1.mod() == 4 && clu2.mod() == 4) {
            mHistManager.fill(HIST("mggMiM44"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 2) || (clu1.mod() == 2 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggMiM12"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 3) || (clu1.mod() == 3 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggMiM13"), m, pt);
          }
          if ((clu1.mod() == 1 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 1)) {
            mHistManager.fill(HIST("mggMiM14"), m, pt);
          }
          if ((clu1.mod() == 2 && clu2.mod() == 3) || (clu1.mod() == 3 && clu2.mod() == 2)) {
            mHistManager.fill(HIST("mggMiM23"), m, pt);
          }
          if ((clu1.mod() == 2 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 2)) {
            mHistManager.fill(HIST("mggMiM24"), m, pt);
          }
          if ((clu1.mod() == 3 && clu2.mod() == 4) || (clu1.mod() == 4 && clu2.mod() == 3)) {
            mHistManager.fill(HIST("mggMiM34"), m, pt);
          }
        }
      }
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosCluQA>(cfgc)};
}
