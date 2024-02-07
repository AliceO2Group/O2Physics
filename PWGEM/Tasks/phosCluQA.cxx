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
#include "Common/DataModel/EventSelection.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/HistogramRegistry.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

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

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

using mcClusters = soa::Join<aod::CaloClusters, aod::PHOSCluLabels>;
using mcAmpClusters = soa::Join<aod::CaloAmbiguousClusters, aod::PHOSAmbCluLabels>;

struct phosCluQA {

  ConfigurableAxis amplitudeAxisLarge{"amplitude", {1000, 0., 100.}, "Amplutude (GeV)"};
  ConfigurableAxis timeAxisLarge{"celltime", {1000, -1500.e-9, 3500.e-9}, "cell time (ns)"};
  ConfigurableAxis mggAxis{"mgg", {250, 0., 1.}, "m_{#gamma#gamma} (GeV/c^{2})"};
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<int> mEvSelTrigAmb{"mEvSelTrigAmb", kTVXinPHOS, "Select events (for ambigious clusters) with this trigger"};
  Configurable<double> mMinCluE{"minCluE", 0.3, "Minimum cluster energy for histograms."};
  Configurable<double> mCpvMinE{"cpvMinE", 200, "Min CPV amplitude"};

  o2::framework::HistogramRegistry mHistManager{"phosCluQAHistograms"};
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  std::unique_ptr<o2::phos::Geometry> geomPHOS;

  bool fillBCmap = true;
  std::array<TH2*, 25> hRe, hMi, hReAmb, hMiAmb;
  std::array<TH2*, 4> hCluMul, hCluMulAmb, hCluETime, hCluETimeAmb;
  std::array<TH1*, 4> hCluSp, hCluSpAmb;

  /// \brief Create output histograms
  void init(o2::framework::InitContext const&)
  {
    using o2HistType = o2::framework::HistType;
    using o2Axis = o2::framework::AxisSpec;
    LOG(info) << "Initializing PHOS Cluster QA monitor task ...";

    const o2Axis
      multAxis{100, 0., 100.},
      bcAxis{3501, -0.5, 3500.5},
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"};
    auto h{std::get<std::shared_ptr<TH1>>(mHistManager.add("bcAll", "Number of BC", o2HistType::kTH1F, {{4, 0.5, 4.5}}))};
    h->GetXaxis()->SetBinLabel(1, "All BC");
    h->GetXaxis()->SetBinLabel(2, "BC w/o coll.");
    h->GetXaxis()->SetBinLabel(3, "BC w one coll.");
    h->GetXaxis()->SetBinLabel(4, "BC w few coll.");

    mHistManager.add("BCA", "Bunch crossing schedule A only", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("BCC", "Bunch crossing schedule C only", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("BCB", "Bunch crossing schedule Both", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cpvBCAll", "Bunch crossing ID of event with CPV", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCAll", "Bunch crossing ID of event with PHOS", o2HistType::kTH1F, {bcAxis});
    mHistManager.add("ambcluBCAll", "Bunch crossing ID of event with PHOS", o2HistType::kTH1F, {bcAxis});

    auto h2{std::get<std::shared_ptr<TH2>>(mHistManager.add("cpvPhosEvents", "Number of common PHOS and CPV events", o2HistType::kTH2F, {{3, 0., 3.}, {100, 0., 100}}))};
    h2->GetXaxis()->SetBinLabel(1, "CPV with PHOS");
    h2->GetXaxis()->SetBinLabel(2, "CPV w/o PHOS");
    h2->GetXaxis()->SetBinLabel(3, "PHOS w/o CPV");

    auto h3{std::get<std::shared_ptr<TH2>>(mHistManager.add("cpvAmbPhosEvents", "Number of common PHOS and CPV events", o2HistType::kTH2F, {{3, 0., 3.}, {100, 0., 100}}))};
    h3->GetXaxis()->SetBinLabel(1, "CPV with PHOS");
    h3->GetXaxis()->SetBinLabel(2, "CPV w/o PHOS");
    h3->GetXaxis()->SetBinLabel(3, "PHOS w/o CPV");

    mHistManager.add("eventsTrig", "Number of trigger events", HistType::kTH1F, {{2, 0., 2.}});
    mHistManager.add("CPVSp", "CPV spectrum", o2HistType::kTH2F, {{200, 0., 10000}, {3, 2., 5.}});
    mHistManager.add("CPVOcc", "CPV occupancy", o2HistType::kTH3F, {{128, -72., 72.}, {56, -63., 63.}, {3, 2., 5.}});

    for (int i = 2; i < 5; ++i) {
      mHistManager.add(Form("CPVPHOSDistReMod%d", i), Form("PHOS-CPV distance for module %d", i),
                       o2HistType::kTH3F, {{100, -100., 100.}, {100, -100., 100.}, {100, 0., 10.}});
      mHistManager.add(Form("CPVPHOSDistMiMod%d", i), Form("PHOS-CPV distance for module %d", i),
                       o2HistType::kTH3F, {{100, -100., 100.}, {100, -100., 100.}, {100, 0., 10.}});
    }

    for (int i = 1; i < 5; ++i) {
      hCluSp[i - 1] = (std::get<std::shared_ptr<TH1>>(mHistManager.add(Form("cluSpM%d", i), Form("Cluster spectrum for module %d", i),
                                                                       o2HistType::kTH1F, {amplitudeAxisLarge})))
                        .get();
      hCluMul[i - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("cluMultM%d", i), Form("Cluster multiplicity for module %d", i),
                                                                        o2HistType::kTH2F, {amplitudeAxisLarge, multAxis})))
                         .get();
      hCluETime[i - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("cluETimeM%d", i), Form("Correlation between cell amplitude and time in module %d", i),
                                                                          o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge})))
                           .get();

      hCluSpAmb[i - 1] = (std::get<std::shared_ptr<TH1>>(mHistManager.add(Form("ambcluSpM%d", i), Form("Cluster spectrum for module %d", i),
                                                                          o2HistType::kTH1F, {amplitudeAxisLarge})))
                           .get();
      hCluMulAmb[i - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("ambcluMultM%d", i), Form("Cluster multiplicity for module %d", i),
                                                                           o2HistType::kTH2F, {amplitudeAxisLarge, multAxis})))
                            .get();
      hCluETimeAmb[i - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("ambcluETimeM%d", i), Form("Correlation between cell amplitude and time in module %d", i),
                                                                             o2HistType::kTH2F, {timeAxisLarge, amplitudeAxisLarge})))
                              .get();
    }
    mHistManager.add("cluOcc", "Cluster occupancy ", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluE", "Cluster energy", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluOcc", "Cluster occupancy ", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluE", "Cluster energy", o2HistType::kTH3F, {phiAxis, zAxis, modAxis});
    for (int i = 1; i < 5; ++i) {
      for (int j = i; j < 5; ++j) {
        hRe[5 * (i - 1) + j - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("mggReM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                                                                                    o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge})))
                                     .get();
        hMi[5 * (i - 1) + j - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("mggMiM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                                                                                    o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge})))
                                     .get();
        hReAmb[5 * (i - 1) + j - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("ambmggReM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                                                                                       o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge})))
                                        .get();
        hMiAmb[5 * (i - 1) + j - 1] = (std::get<std::shared_ptr<TH2>>(mHistManager.add(Form("ambmggMiM%d%d", i, j), Form("inv mass for module %d%d", i, j),
                                                                                       o2HistType::kTH2F, {mggAxis, amplitudeAxisLarge})))
                                        .get();
      }
    }
    if (isMC) {
      mHistManager.add("cluNprim", "Number of primaries", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {10, 0., 10., "N_{prim}", "N_{prim}"}});
      mHistManager.add("cluEdep", "Deposited energy", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {20, 0., 2., "E_{dep}/E_{rec}", "E_{dep}/E_{rec}"}});
      mHistManager.add("cluEdepN", "Deposited energy vs iparent", HistType::kTH2F, {{10, 0., 10., "i_{parent}", "i_{parent}"}, {20, 0., 2., "E_{dep}/E_{rec}", "E_{dep}/E_{rec}"}});
      mHistManager.add("mcEdepAll", "Erec vs true E", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {100, 0., 10., "E_{true}", "E_{true} (GeV)"}});
      mHistManager.add("mcEdepGamma", "Erec vs true E", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {100, 0., 10., "E_{true}", "E_{true} (GeV)"}});
      mHistManager.add("cluDxe", "dx vs E", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {100, -50., 50., "x_{rec}-x_{true}", "x_{rec}-x_{true} (cm)"}});
      mHistManager.add("cluDze", "dz vs E", HistType::kTH2F, {{100, 0., 10., "E_{rec}", "E_{rec} (GeV)"}, {100, -50., 50., "z_{rec}-z_{true}", "z_{rec}-z_{true} (cm)"}});
    }

    ccdb->setURL(o2::base::NameConf::getCCDBServer());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
  }

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  /// \brief Process PHOS data
  void processData(BCsWithBcSels const& bcs,
                   SelCollisions const& colls,
                   o2::aod::CaloClusters const& clusters,
                   o2::aod::CaloAmbiguousClusters const& ambclusters,
                   o2::aod::CPVClusters const& cpvs)
  {
    FillQCHistos(bcs, colls);

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

    o2::InteractionRecord ir;
    for (const auto& clu : clusters) {
      ir.setFromLong(clu.collision().bc().globalBC());
      mHistManager.fill(HIST("cluBCAll"), ir.bc);

      if (!clu.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().alias_bit(mEvSelTrig))
        continue;
      mHistManager.fill(HIST("eventsTrig"), 1.);

      hCluSp[clu.mod() - 1]->Fill(clu.e());
      hCluMul[clu.mod() - 1]->Fill(clu.e(), clu.ncell());
      hCluETime[clu.mod() - 1]->Fill(clu.time(), clu.e());
      if (clu.e() > mMinCluE) {
        mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
        mHistManager.fill(HIST("cluE"), clu.x(), clu.z(), clu.mod(), clu.e());
      }
    }
    for (const auto& clu : ambclusters) {
      ir.setFromLong(clu.bc().globalBC());
      mHistManager.fill(HIST("ambcluBCAll"), ir.bc);

      if (!clu.bc_as<BCsWithBcSels>().alias_bit(mEvSelTrigAmb))
        continue;
      mHistManager.fill(HIST("eventsTrigAmb"), 1.);

      hCluSpAmb[clu.mod() - 1]->Fill(clu.e());
      hCluMulAmb[clu.mod() - 1]->Fill(clu.e(), clu.ncell());
      hCluETimeAmb[clu.mod() - 1]->Fill(clu.time(), clu.e());

      if (clu.e() > mMinCluE) {
        mHistManager.fill(HIST("ambcluOcc"), clu.x(), clu.z(), clu.mod());
        mHistManager.fill(HIST("ambcluE"), clu.x(), clu.z(), clu.mod(), clu.e());
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
          if (clu1.mod() < clu2.mod()) {
            hRe[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hRe[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
        } else { // Mixed
          if (clu1.mod() < clu2.mod()) {
            hMi[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hMi[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
          --nMix;
        }
      }
    }

    // inv mass
    for (const auto& clu1 : ambclusters) {
      auto clu2 = clu1;
      ++clu2;
      int nMix = 100; // Number of photons to mix
      for (; clu2 != ambclusters.end() && nMix > 0; clu2++) {
        double m = pow(clu1.e() + clu2.e(), 2) - pow(clu1.px() + clu2.px(), 2) -
                   pow(clu1.py() + clu2.py(), 2) - pow(clu1.pz() + clu2.pz(), 2);
        if (m > 0)
          m = sqrt(m);
        double pt = sqrt(pow(clu1.px() + clu2.px(), 2) +
                         pow(clu1.py() + clu2.py(), 2));
        if (clu1.bc() == clu2.bc()) { // Real
          if (clu1.mod() < clu2.mod()) {
            hReAmb[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hReAmb[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
        } else { // Mixed
          if (clu1.mod() < clu2.mod()) {
            hMiAmb[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hMiAmb[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
          --nMix;
        }
      }
    }
  }
  PROCESS_SWITCH(phosCluQA, processData, "Process real data", true);

  void processMC(BCsWithBcSels const& bcs,
                 SelCollisions const& colls,
                 aod::McParticles const& mcParticles,
                 mcClusters const& clusters,
                 mcAmpClusters const& ambclusters,
                 o2::aod::CPVClusters const& cpvs)
  {
    FillQCHistos(bcs, colls);

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

    o2::InteractionRecord ir;
    for (const auto& clu : clusters) {
      ir.setFromLong(clu.collision().bc().globalBC());
      mHistManager.fill(HIST("cluBCAll"), ir.bc);

      if (!clu.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().alias_bit(mEvSelTrig))
        continue;
      mHistManager.fill(HIST("eventsTrig"), 1.);

      hCluSp[clu.mod() - 1]->Fill(clu.e());
      hCluMul[clu.mod() - 1]->Fill(clu.e(), clu.ncell());
      hCluETime[clu.mod() - 1]->Fill(clu.time(), clu.e());
      if (clu.e() > mMinCluE) {
        mHistManager.fill(HIST("cluOcc"), clu.x(), clu.z(), clu.mod());
        mHistManager.fill(HIST("cluE"), clu.x(), clu.z(), clu.mod(), clu.e());
      }

      // MC energy, position resolution etc
      auto mcList = clu.labels();     // const std::vector<int>
      auto mcEdep = clu.amplitides(); // const std::vector<float>&
      mHistManager.fill(HIST("cluNprim"), clu.e(), mcList.size());
      for (size_t iii = 0; iii < mcList.size(); iii++) {
        mHistManager.fill(HIST("cluEdep"), clu.e(), mcEdep[iii] / clu.e());
        mHistManager.fill(HIST("cluEdepN"), static_cast<float>(iii), mcEdep[iii] / clu.e());
      }

      if (mcList.size() > 0) {
        float mce = mcParticles.iteratorAt(mcList[0]).e();
        mHistManager.fill(HIST("mcEdepAll"), clu.e(), mce);
        if (mcParticles.iteratorAt(mcList[0]).pdgCode() == 22) {
          mHistManager.fill(HIST("mcEdepGamma"), clu.e(), mce);
          // if photon from veretx, compare position resolution
          if (pow(mcParticles.iteratorAt(mcList[0]).vx(), 2) + pow(mcParticles.iteratorAt(mcList[0]).vy(), 2) < 1.) {
            // calculate impact position on PHOS
            TVector3 vtx(mcParticles.iteratorAt(mcList[0]).vx(), mcParticles.iteratorAt(mcList[0]).vy(), mcParticles.iteratorAt(mcList[0]).vz());
            TVector3 p(mcParticles.iteratorAt(mcList[0]).px(), mcParticles.iteratorAt(mcList[0]).py(), mcParticles.iteratorAt(mcList[0]).pz());
            int16_t mod;
            float x, z;
            if (geomPHOS->impactOnPHOS(vtx, p, mod, z, x)) { // photon should hit PHOS
              if (mod == clu.mod()) {
                mHistManager.fill(HIST("cluDxe"), clu.e(), clu.x() - x);
                mHistManager.fill(HIST("cluDze"), clu.e(), clu.z() - z);
              }
            }
          }
        }
      }
    }
    for (const auto& clu : ambclusters) {
      ir.setFromLong(clu.bc().globalBC());
      mHistManager.fill(HIST("ambcluBCAll"), ir.bc);

      if (!clu.bc_as<BCsWithBcSels>().alias_bit(mEvSelTrigAmb))
        continue;
      mHistManager.fill(HIST("eventsTrigAmb"), 1.);

      hCluSpAmb[clu.mod() - 1]->Fill(clu.e());
      hCluMulAmb[clu.mod() - 1]->Fill(clu.e(), clu.ncell());
      hCluETimeAmb[clu.mod() - 1]->Fill(clu.time(), clu.e());

      if (clu.e() > mMinCluE) {
        mHistManager.fill(HIST("ambcluOcc"), clu.x(), clu.z(), clu.mod());
        mHistManager.fill(HIST("ambcluE"), clu.x(), clu.z(), clu.mod(), clu.e());
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
          if (clu1.mod() < clu2.mod()) {
            hRe[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hRe[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
        } else { // Mixed
          if (clu1.mod() < clu2.mod()) {
            hMi[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hMi[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
          --nMix;
        }
      }
    }

    // inv mass
    for (const auto& clu1 : ambclusters) {
      auto clu2 = clu1;
      ++clu2;
      int nMix = 100; // Number of photons to mix
      for (; clu2 != ambclusters.end() && nMix > 0; clu2++) {
        double m = pow(clu1.e() + clu2.e(), 2) - pow(clu1.px() + clu2.px(), 2) -
                   pow(clu1.py() + clu2.py(), 2) - pow(clu1.pz() + clu2.pz(), 2);
        if (m > 0)
          m = sqrt(m);
        double pt = sqrt(pow(clu1.px() + clu2.px(), 2) +
                         pow(clu1.py() + clu2.py(), 2));
        if (clu1.bc() == clu2.bc()) { // Real
          if (clu1.mod() < clu2.mod()) {
            hReAmb[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hReAmb[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
        } else { // Mixed
          if (clu1.mod() < clu2.mod()) {
            hMiAmb[5 * (clu1.mod() - 1) + clu2.mod() - 1]->Fill(m, pt);
          } else {
            hMiAmb[5 * (clu2.mod() - 1) + clu1.mod() - 1]->Fill(m, pt);
          }
          --nMix;
        }
      }
    }
  }
  PROCESS_SWITCH(phosCluQA, processMC, "Process MC data", false);

  void FillQCHistos(BCsWithBcSels const& bcs,
                    SelCollisions const& colls)
  {
    // Filll BC map
    if (fillBCmap && bcs.begin() != bcs.end()) {
      auto rl = ccdb->getRunDuration(bcs.begin().runNumber());
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", rl.first);
      constexpr int nBCsPerOrbit = 3564;
      std::bitset<nBCsPerOrbit> beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
      std::bitset<nBCsPerOrbit> beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
      std::bitset<nBCsPerOrbit> bcPatternA = beamPatternA & ~beamPatternC;
      std::bitset<nBCsPerOrbit> bcPatternC = ~beamPatternA & beamPatternC;
      std::bitset<nBCsPerOrbit> bcPatternB = beamPatternA & beamPatternC;
      for (int i = 0; i < nBCsPerOrbit; i++) {
        if (bcPatternB[i])
          mHistManager.fill(HIST("BCB"), i);
        if (bcPatternA[i])
          mHistManager.fill(HIST("BCA"), i);
        if (bcPatternC[i])
          mHistManager.fill(HIST("BCC"), i);
      }
      fillBCmap = false;
    }

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

    // Fill histograms to see energy resolution, position resolution
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return o2::framework::WorkflowSpec{
    o2::framework::adaptAnalysisTask<phosCluQA>(cfgc)};
}
