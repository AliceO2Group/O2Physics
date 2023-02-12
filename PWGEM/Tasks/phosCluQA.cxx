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
  o2::framework::Configurable<double> mCpvMinE{"cpvMinE", 200, "Min CPV amplitude"};

  o2::framework::HistogramRegistry mHistManager{"phosCluQAHistograms"};

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    LOG(info) << "Initializing PHOS Cluster QA monitor task ...";

    const o2Axis
      amplitudeAxisLarge{1000, 0., 10., "amplitude", "Amplutude (GeV)"},
      timeAxisLarge{1000, -1500.e-9, 3500.e-9, "celltime", "cell time (ns)"},
      multAxis{100, 0., 100.},
      mggAxis{250, 0., 1.},
      bcAxis{3501, -0.5, 3500.5},
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"};
    mHistManager.add("bcAll", "Number of BC", o2HistType::kTH1F, {{4, 0.5, 4.5}});
    mHistManager.add("cpvBCAll", "Bunch crossing ID of event with CPV", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCAll", "Bunch crossing ID of event with PHOS", o2HistType::kTH1F, {bcAxis});

    mHistManager.add("cpvPhosEvents", "Number of common PHOS and CPV events", o2HistType::kTH2F, {{3, 0., 3.}, {100, 0., 100}});
    mHistManager.add("cpvAmbPhosEvents", "Number of common PHOS and CPV events", o2HistType::kTH2F, {{3, 0., 3.}, {100, 0., 100}});

    mHistManager.add("CPVSp", "CPV spectrum", o2HistType::kTH2F, {{200, 0., 10000}, {3, 2., 5.}});
    mHistManager.add("CPVOcc", "CPV occupancy", o2HistType::kTH3F, {{128, -72., 72.}, {56, -63., 63.}, {3, 2., 5.}});

    for (int i = 2; i < 5; ++i) {
      mHistManager.add(Form("CPVPHOSDistReMod%d", i), Form("PHOS-CPV distance for module %d", i),
                       o2HistType::kTH3F, {{100, -100., 100.}, {100, -100., 100.}, {100, 0., 10.}});
      mHistManager.add(Form("CPVPHOSDistMiMod%d", i), Form("PHOS-CPV distance for module %d", i),
                       o2HistType::kTH3F, {{100, -100., 100.}, {100, -100., 100.}, {100, 0., 10.}});
    }

    for (int i = 1; i < 5; ++i) {
      mHistManager.add(Form("cluSpM%d", i), Form("Cluster spectrum for module %d", i),
                       o2HistType::kTH1F, {amplitudeAxisLarge});
      mHistManager.add(Form("cluMultM%d", i), Form("Cluster multiplicity for module %d", i),
                       o2HistType::kTH2F, {amplitudeAxisLarge, multAxis});
      mHistManager.add(Form("cluETimeM%d", i), Form("Correlation between cell amplitude and time in module %d", i),
                       o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge});
    }
    mHistManager.add("cluOcc", "Cluster occupancy ", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluE", "Cluster energy", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
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
               o2::aod::Collisions const& colls,
               o2::aod::CaloClusters const& clusters,
               o2::aod::CaloAmbiguousClusters const& ambclusters,
               o2::aod::CPVClusters const& cpvs)
  {

    // If several collisions appear in BC, choose one with largers number of contributors
    std::map<int64_t, int> colMap;
    for (auto cl : colls) {
      auto colbc = colMap.find(cl.bc().globalBC());
      if (colbc == colMap.end()) { // single collision per BC
        colMap[cl.bc().globalBC()] = 1;
      } else { // not unique collision per BC
        colbc->second++;
      }
    }
    for (const auto& bc : bcs) {
      o2::InteractionRecord eventIR;
      eventIR.setFromLong(bc.globalBC());
      mHistManager.fill(HIST("bcAll"), 1);
      auto colbc = colMap.find(bc.globalBC());

      if (colbc == colMap.end()) { // no Collision
        mHistManager.fill(HIST("bcAll"), 2);
      } else {
        if (colbc->second == 1) {
          mHistManager.fill(HIST("bcAll"), 3);
        } else {
          mHistManager.fill(HIST("bcAll"), 4);
        }
      }
    }

    std::map<uint64_t, int> ncpvClu;
    for (const auto& cpvclu : cpvs) {
      ncpvClu[cpvclu.bcId()]++;
      o2::InteractionRecord ir;
      ir.setFromLong(cpvclu.bc().globalBC());
      mHistManager.fill(HIST("cpvBCAll"), ir.bc);
    }
    std::map<uint64_t, int> nphosClu;
    for (const auto& clu : clusters) {
      nphosClu[clu.collision().bcId()]++;
    }
    std::map<uint64_t, int> nphosAmbClu;
    for (const auto& clu : ambclusters) {
      nphosAmbClu[clu.bcId()]++;
    }
    for (const auto& [bcid, n] : ncpvClu) {
      if (nphosClu.find(bcid) != nphosClu.end()) {
        mHistManager.fill(HIST("cpvPhosEvents"), 0., n);
      } else {
        mHistManager.fill(HIST("cpvPhosEvents"), 1., n);
      }
    }
    for (const auto& [bcid, n] : ncpvClu) {
      if (nphosAmbClu.find(bcid) != nphosAmbClu.end()) {
        mHistManager.fill(HIST("cpvAmbPhosEvents"), 0., n);
      } else {
        mHistManager.fill(HIST("cpvAmbPhosEvents"), 1., n);
      }
    }

    for (const auto& [bcid, n] : nphosClu) {
      if (ncpvClu.find(bcid) == ncpvClu.end()) {
        mHistManager.fill(HIST("cpvPhosEvents"), 2., n);
      }
    }

    for (const auto& [bcid, n] : nphosAmbClu) {
      if (ncpvClu.find(bcid) == ncpvClu.end()) {
        mHistManager.fill(HIST("cpvAmbPhosEvents"), 2., n);
      }
    }

    for (const auto& cpvclu : cpvs) {
      int bcCPV = cpvclu.bcId();
      mHistManager.fill(HIST("CPVSp"), cpvclu.amplitude(), cpvclu.moduleNumber());
      if (cpvclu.amplitude() < mCpvMinE) {
        continue;
      }
      mHistManager.fill(HIST("CPVOcc"), cpvclu.posX(), cpvclu.posZ(), cpvclu.moduleNumber());

      // CPV-PHOS matching
      int currentPHOSBc = 0;
      double dist = 9999.;
      double dphi = 0., dz = 0., eClose = 0;
      for (const auto& clu : clusters) {
        if (clu.mod() != cpvclu.moduleNumber()) {
          continue;
        }
        int bcPHOS = clu.collision().bcId();
        if (currentPHOSBc == 0) {
          currentPHOSBc = bcPHOS;
        }
        if (bcPHOS != currentPHOSBc) {  // switched to new event
          if (currentPHOSBc == bcCPV) { // Real
            if (clu.mod() == 2) {
              mHistManager.fill(HIST("CPVPHOSDistReMod2"), dphi, dz, eClose);
            } else {
              if (clu.mod() == 3) {
                mHistManager.fill(HIST("CPVPHOSDistReMod3"), dphi, dz, eClose);
              } else {
                if (clu.mod() == 4) {
                  mHistManager.fill(HIST("CPVPHOSDistReMod4"), dphi, dz, eClose);
                }
              }
            }
          } else {
            if (clu.mod() == 2) {
              mHistManager.fill(HIST("CPVPHOSDistMiMod2"), dphi, dz, eClose);
            } else {
              if (clu.mod() == 3) {
                mHistManager.fill(HIST("CPVPHOSDistMiMod3"), dphi, dz, eClose);
              } else {
                if (clu.mod() == 4) {
                  mHistManager.fill(HIST("CPVPHOSDistMiMod4"), dphi, dz, eClose);
                }
              }
            }
          }
          currentPHOSBc = bcPHOS;
          dist = 9999.;
          dphi = 999.;
          dz = 999.;
        }
        double dr2 = pow(clu.x() - cpvclu.posX(), 2) + pow(clu.z() - cpvclu.posZ(), 2);
        if (dr2 < dist) {
          dist = dr2;
          dz = clu.z() - cpvclu.posZ();
          dphi = clu.x() - cpvclu.posX();
          eClose = clu.e();
        }
      } // phos clusters
      // last event
      if (currentPHOSBc == bcCPV) { // Real
        if (cpvclu.moduleNumber() == 2) {
          mHistManager.fill(HIST("CPVPHOSDistReMod2"), dphi, dz, eClose);
        } else {
          if (cpvclu.moduleNumber() == 3) {
            mHistManager.fill(HIST("CPVPHOSDistReMod3"), dphi, dz, eClose);
          } else {
            if (cpvclu.moduleNumber() == 4) {
              mHistManager.fill(HIST("CPVPHOSDistReMod4"), dphi, dz, eClose);
            }
          }
        }
      } else {
        if (cpvclu.moduleNumber() == 2) {
          mHistManager.fill(HIST("CPVPHOSDistMiMod2"), dphi, dz, eClose);
        } else {
          if (cpvclu.moduleNumber() == 3) {
            mHistManager.fill(HIST("CPVPHOSDistMiMod3"), dphi, dz, eClose);
          } else {
            if (cpvclu.moduleNumber() == 4) {
              mHistManager.fill(HIST("CPVPHOSDistMiMod4"), dphi, dz, eClose);
            }
          }
        }
      }
    }

    for (const auto& clu : clusters) {
      mHistManager.fill(HIST("cluBCAll"), clu.collision().bc().globalBC());

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
        mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
        mHistManager.fill(HIST("cluE"), clu.x(), clu.z(), clu.mod(), clu.e());
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
        if (clu1.collision() == clu2.collision()) { // Real
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
