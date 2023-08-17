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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "PHOSBase/Geometry.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

/// \struct PHOS pi0 analysis
/// \brief Monitoring task for PHOS related quantities
/// \author Dmitri Peresunko, NRC "Kurchatov institute"
/// \since Nov, 2022
///

using namespace o2;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct phosPi0 {

  Configurable<float> mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -50.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"mMinCluTime", 100.e-9, "Max. cluster time"};
  Configurable<int> mMinCluNcell{"mMinCluNcell", 2, "min cells in cluster"};
  Configurable<int> mMixedEvents{"mixedEvents", 10, "number of events to mix"};
  Configurable<bool> mEventSelection{"mEventSelection", true, "to apply event selection"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};

  Configurable<float> mOccE{"minOccE", 0.6, "Minimum cluster energy to fill occupancy"};

  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  using FilteredClusters = soa::Filtered<aod::CaloClusters>;

  HistogramRegistry mHistManager{"phosPi0Histograms"};

  bool fillBCmap = true; // fill BC map once
  static constexpr int nBCsPerOrbit = 3564;
  std::bitset<nBCsPerOrbit> mbcPatternB;
  std::bitset<nBCsPerOrbit> mbcPatternA;
  std::bitset<nBCsPerOrbit> mbcPatternC;

  /// \brief Create output histograms
  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS pi0 analysis task ...";

    const AxisSpec
      bcAxis{3501, -0.5, 3500.5},
      zAxis{56, -63., 63., "z", "z (cm)"},
      phiAxis{64, -72., 72., "x", "x (cm)"},
      modAxis{4, 1., 5., "module", "Module"},
      amplitudeAxisLarge{100, 0., 20., "amplitude", "Amplutude (GeV)"},
      timeAxisLarge{500, -1500.e-9, 3500.e-9, "celltime", "cell time (ns)"},
      multAxis{100, 0., 100.},
      mggAxis{250, 0., 1., "mgg", "m_{#gamma#gamma} (GeV/c^{2})"},
      vertexAxis{100, -20., 20., "z", "z (cm)"},
      modCombAxis{10, 0., 10.},
      centAxis{10, 0., 10.},
      centralityAxis{100, 0., 100., "centrality", "centrality"};

    auto h{std::get<std::shared_ptr<TH1>>(mHistManager.add("eventsCol", "Number of events", HistType::kTH1F, {{9, 0., 9.}}))};
    h->GetXaxis()->SetBinLabel(1, "All");
    h->GetXaxis()->SetBinLabel(2, "T0a||T0c");
    h->GetXaxis()->SetBinLabel(3, "T0a&&T0c");
    h->GetXaxis()->SetBinLabel(4, "V0A");
    h->GetXaxis()->SetBinLabel(5, "kIsTriggerTVX");
    h->GetXaxis()->SetBinLabel(6, "kTVXinPHOS");
    h->GetXaxis()->SetBinLabel(7, "PHOSClu");
    h->GetXaxis()->SetBinLabel(8, "PHOSClu&&Trig");
    h->GetXaxis()->SetBinLabel(9, "PHOSClu&&Trig&&BB");

    auto h2{std::get<std::shared_ptr<TH1>>(mHistManager.add("eventsBC", "Number of events per trigger", HistType::kTH1F, {{9, 0., 9.}}))};
    h2->GetXaxis()->SetBinLabel(1, "All");
    h2->GetXaxis()->SetBinLabel(2, "T0a||T0c");
    h2->GetXaxis()->SetBinLabel(3, "T0a&&T0c");
    h2->GetXaxis()->SetBinLabel(4, "V0A");
    h2->GetXaxis()->SetBinLabel(5, "kIsTriggerTVX");
    h2->GetXaxis()->SetBinLabel(6, "kTVXinPHOS");
    h2->GetXaxis()->SetBinLabel(7, "PHOSClu");
    h2->GetXaxis()->SetBinLabel(8, "PHOSClu&&Trig");
    h2->GetXaxis()->SetBinLabel(9, "PHOSClu&&Trig&&BB");

    mHistManager.add("contributors", "Centrality", HistType::kTH1F, {centralityAxis});
    mHistManager.add("vertex", "vertex", HistType::kTH1F, {vertexAxis});
    mHistManager.add("vertex2", "PHOS with vertex", HistType::kTH1F, {vertexAxis});
    mHistManager.add("BCA", "Bunch crossing schedule A only", HistType::kTH1F, {bcAxis});
    mHistManager.add("BCC", "Bunch crossing schedule C only", HistType::kTH1F, {bcAxis});
    mHistManager.add("BCB", "Bunch crossing schedule Both", HistType::kTH1F, {bcAxis});
    mHistManager.add("cpvBCAll", "Bunch crossing ID of event with CPV", HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCAll", "Bunch crossing ID of event with PHOS", HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCAll2", "Bunch crossing ID of event with PHOS clu", HistType::kTH1F, {bcAxis});
    mHistManager.add("ambcluBCAll", "Bunch crossing ID of event with PHOS", HistType::kTH1F, {bcAxis});
    mHistManager.add("evBCkTVXinPHOS", "Bunch crossing ID of event with PHOS", HistType::kTH1F, {bcAxis});
    mHistManager.add("cluBCkTVXinPHOS", "Bunch crossing ID of event with PHOS", HistType::kTH1F, {bcAxis});
    mHistManager.add("ambCluBCkTVXinPHOS", "Bunch crossing ID of event with PHOS", HistType::kTH1F, {bcAxis});

    mHistManager.add("cluSp", "Cluster spectrum per module",
                     HistType::kTH2F, {amplitudeAxisLarge, modAxis});
    mHistManager.add("ambcluSp", "Amb. Cluster spectrum per module",
                     HistType::kTH2F, {amplitudeAxisLarge, modAxis});
    mHistManager.add("cluETime", "Cluster time vs E",
                     HistType::kTH3F, {amplitudeAxisLarge, timeAxisLarge, modAxis});
    mHistManager.add("ambcluETime", "Amb. cluster time vs E",
                     HistType::kTH3F, {amplitudeAxisLarge, timeAxisLarge, modAxis});
    mHistManager.add("mggRe", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMi", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReBoth", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiBoth", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});

    mHistManager.add("mggReAmb", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiAmb", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReAmbCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiAmbCPV", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReAmbDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiAmbDisp", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggReAmbBoth", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});
    mHistManager.add("mggMiAmbBoth", "inv mass for centrality",
                     HistType::kTH3F, {mggAxis, amplitudeAxisLarge, modCombAxis});

    mHistManager.add("cluOcc", "Cluster occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluCPVOcc", "Cluster with CPV occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluDispOcc", "Cluster with Disp occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluBothOcc", "Cluster with Both occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("cluE", "Cluster energy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluOcc", "Cluster occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluCPVOcc", "Cluster with CPV occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluDispOcc", "Cluster with Disp occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluBothOcc", "Cluster with Both occupancy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
    mHistManager.add("ambcluE", "Cluster energy", HistType::kTH3F, {phiAxis, zAxis, modAxis});
  }

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using BCsWithBcSels = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;

  /// \brief Process PHOS data
  void process(BCsWithBcSels const& bcs,
               SelCollisions const& collisions,
               aod::CaloClusters const& clusters,
               aod::CaloAmbiguousClusters const& ambclusters)
  {

    // Filll BC map
    if (fillBCmap && bcs.begin() != bcs.end()) {
      auto rl = ccdb->getRunDuration(bcs.begin().runNumber());
      auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", rl.first);
      std::bitset<nBCsPerOrbit> beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
      std::bitset<nBCsPerOrbit> beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
      mbcPatternB = beamPatternA & beamPatternC;
      mbcPatternA = beamPatternA & ~beamPatternC;
      mbcPatternC = ~beamPatternA & beamPatternC;
      for (int i = 0; i < nBCsPerOrbit; i++) {
        if (mbcPatternB[i])
          mHistManager.fill(HIST("BCB"), i);
        if (mbcPatternA[i])
          mHistManager.fill(HIST("BCA"), i);
        if (mbcPatternC[i])
          mHistManager.fill(HIST("BCC"), i);
      }
      fillBCmap = false;
    }
    o2::InteractionRecord ir;
    for (const auto& col : collisions) {
      mHistManager.fill(HIST("contributors"), col.numContrib());
      mHistManager.fill(HIST("vertex"), col.posZ());
      ir.setFromLong(col.bc_as<BCsWithBcSels>().globalBC());
      mHistManager.fill(HIST("cluBCAll"), ir.bc);
      if (col.alias_bit(kTVXinPHOS)) {
        mHistManager.fill(HIST("evBCkTVXinPHOS"), ir.bc);
      }
      mHistManager.fill(HIST("eventsCol"), 0.);
      if (col.selection_bit(kIsBBT0A) || col.selection_bit(kIsBBT0C)) {
        mHistManager.fill(HIST("eventsCol"), 1);
      }
      if (col.selection_bit(kIsBBT0A) && col.selection_bit(kIsBBT0C)) {
        mHistManager.fill(HIST("eventsCol"), 2);
      }
      if (col.selection_bit(kIsBBV0A)) {
        mHistManager.fill(HIST("eventsCol"), 3);
      }
      if (col.selection_bit(kIsTriggerTVX)) {
        mHistManager.fill(HIST("eventsCol"), 4);
      }
      if (col.alias_bit(kTVXinPHOS)) {
        mHistManager.fill(HIST("eventsCol"), 5);
      }
    }

    for (const auto& bc : bcs) {
      mHistManager.fill(HIST("eventsBC"), 0);
      if (bc.selection_bit(kIsBBT0A) || bc.selection_bit(kIsBBT0C)) {
        mHistManager.fill(HIST("eventsBC"), 1);
      }
      if (bc.selection_bit(kIsBBT0A) && bc.selection_bit(kIsBBT0C)) {
        mHistManager.fill(HIST("eventsBC"), 2);
      }
      if (bc.selection_bit(kIsBBV0A)) {
        mHistManager.fill(HIST("eventsBC"), 3);
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mHistManager.fill(HIST("eventsBC"), 4);
      }
      if (bc.alias_bit(kTVXinPHOS)) {
        mHistManager.fill(HIST("eventsBC"), 5);
      }
    }

    uint64_t bcevent = 0;
    for (const auto& clu : clusters) {
      if (clu.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().globalBC() != bcevent) { // New BC
        bcevent = clu.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().globalBC();
        ir.setFromLong(bcevent);
        mHistManager.fill(HIST("eventsCol"), 6.);
        if (clu.collision_as<SelCollisions>().alias_bit(mEvSelTrig)) {
          mHistManager.fill(HIST("eventsCol"), 7.);
          if (mbcPatternB[ir.bc]) {
            mHistManager.fill(HIST("eventsCol"), 8.);
            mHistManager.fill(HIST("vertex2"), clu.collision_as<SelCollisions>().posZ());
          }
        }
        mHistManager.fill(HIST("cluBCAll2"), ir.bc);
        if (clu.collision_as<SelCollisions>().alias_bit(kTVXinPHOS)) {
          mHistManager.fill(HIST("cluBCkTVXinPHOS"), ir.bc);
        }
      }

      if (mEventSelection) {
        if (!clu.collision_as<SelCollisions>().alias_bit(mEvSelTrig) || !mbcPatternB[ir.bc]) {
          continue;
        }
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

      // inv mass
      auto clu2 = clu;
      ++clu2;
      int nMix = mMixedEvents; // Number of events to mix
      bool skipMix = false;
      uint64_t bcurrentMix = 0;
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
        int modComb = ModuleCombination(clu.mod(), clu2.mod());
        if (clu.collision() == clu2.collision()) { // Real
          mHistManager.fill(HIST("mggRe"), m, pt, modComb);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggReCPV"), m, pt, modComb);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggReDisp"), m, pt, modComb);
            if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
              mHistManager.fill(HIST("mggReBoth"), m, pt, modComb);
            }
          }
        } else { // Mixed
          if (clu2.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().globalBC() != bcurrentMix) {
            bcurrentMix = clu2.collision_as<SelCollisions>().bc_as<BCsWithBcSels>().globalBC();
            o2::InteractionRecord irMix;
            irMix.setFromLong(bcurrentMix);
            if (mEventSelection) {
              skipMix = !clu2.collision_as<SelCollisions>().alias_bit(mEvSelTrig) || !mbcPatternB[irMix.bc];
            }
            if (!skipMix) {
              --nMix;
            }
          }
          if (skipMix) {
            continue;
          }
          mHistManager.fill(HIST("mggMi"), m, pt, modComb);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggMiCPV"), m, pt, modComb);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggMiDisp"), m, pt, modComb);
            if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
              mHistManager.fill(HIST("mggMiBoth"), m, pt, modComb);
            }
          }
        }
      }
    }
    // }

    // same for amb clusters
    bcevent = 0;
    for (const auto& clu : ambclusters) {
      if (clu.bc_as<BCsWithBcSels>().globalBC() != bcevent) {
        bcevent = clu.bc_as<BCsWithBcSels>().globalBC();
        ir.setFromLong(bcevent);
        mHistManager.fill(HIST("eventsBC"), 6.);
        if (clu.bc_as<BCsWithBcSels>().selection_bit(mEvSelTrig)) {
          mHistManager.fill(HIST("eventsBC"), 7.);
          if (mbcPatternB[ir.bc]) {
            mHistManager.fill(HIST("eventsBC"), 8.);
          }
        }
        mHistManager.fill(HIST("ambcluBCAll"), ir.bc);
        if (clu.bc_as<BCsWithBcSels>().selection_bit(kIsTriggerTVX)) {
          mHistManager.fill(HIST("ambCluBCkTVXinPHOS"), ir.bc);
        }
      }

      if (mEventSelection && (!clu.bc_as<BCsWithBcSels>().selection_bit(mEvSelTrig) || !mbcPatternB[ir.bc])) {
        continue;
      }

      mHistManager.fill(HIST("ambcluETime"), clu.e(), clu.time(), clu.mod());
      if (clu.e() < mMinCluE || clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime) {
        continue;
      }

      mHistManager.fill(HIST("ambcluSp"), clu.e(), clu.mod());
      if (clu.e() > mOccE) {
        mHistManager.fill(HIST("ambcluOcc"), clu.x(), clu.z(), clu.mod());
        if (clu.trackdist() > 2.) {
          mHistManager.fill(HIST("ambcluCPVOcc"), clu.x(), clu.z(), clu.mod());
          if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
            mHistManager.fill(HIST("ambcluBothOcc"), clu.x(), clu.z(), clu.mod());
          }
        }
        if (TestLambda(clu.e(), clu.m02(), clu.m20())) {
          mHistManager.fill(HIST("ambcluDispOcc"), clu.x(), clu.z(), clu.mod());
        }
        mHistManager.fill(HIST("ambcluE"), clu.x(), clu.z(), clu.mod(), clu.e());
      }

      // inv mass
      auto clu2 = clu;
      ++clu2;
      int nMix = mMixedEvents; // Number of events to mix
      uint64_t bcurrentMix = 0;
      bool skipMix = false;
      for (; clu2 != ambclusters.end() && nMix > 0; clu2++) {
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
        int modComb = ModuleCombination(clu.mod(), clu2.mod());
        if (clu.bc_as<BCsWithBcSels>() == clu2.bc_as<BCsWithBcSels>()) { // Real
          mHistManager.fill(HIST("mggReAmb"), m, pt, modComb);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggReAmbCPV"), m, pt, modComb);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggReAmbDisp"), m, pt, modComb);
            if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
              mHistManager.fill(HIST("mggReAmbBoth"), m, pt, modComb);
            }
          }
        } else { // Mixed
          if (clu2.bc_as<BCsWithBcSels>().globalBC() != bcurrentMix) {
            bcurrentMix = clu2.bc_as<BCsWithBcSels>().globalBC();
            o2::InteractionRecord irMix;
            irMix.setFromLong(bcevent);
            if (mEventSelection) {
              skipMix = !clu2.bc_as<BCsWithBcSels>().selection_bit(mEvSelTrig) || !mbcPatternB[irMix.bc];
            }
            if (!skipMix) {
              --nMix;
            }
          }
          if (skipMix) {
            continue;
          }
          mHistManager.fill(HIST("mggMiAmb"), m, pt, modComb);
          if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
            mHistManager.fill(HIST("mggMiAmbCPV"), m, pt, modComb);
          }
          if (TestLambda(clu.e(), clu.m02(), clu.m20()) && TestLambda(clu2.e(), clu2.m02(), clu2.m20())) {
            mHistManager.fill(HIST("mggMiAmbDisp"), m, pt, modComb);
            if (clu.trackdist() > 2. && clu2.trackdist() > 2.) {
              mHistManager.fill(HIST("mggMiAmbBoth"), m, pt, modComb);
            }
          }
        }
      }
    }
  }

  //_____________________________________________________________________________
  int ModuleCombination(int m1, int m2)
  {
    // enumerates possible module combinations
    // (1,1)=0, (2,2)=1, (3,3)=2, (4,4)=3, (1,2)=(2,1)=4, (2,3)=(3,2)=5, (3,4)=(4,3)=6, (1,3)=(3,1)=7,
    // (2,4)=(4,2)=8, (1,4)=(4,1)=9
    int d = TMath::Abs(m1 - m2);
    if (d == 0) {
      return m1 - 1;
    }
    if (d == 1) {
      return 3 + TMath::Min(m1, m2);
    }
    if (d == 2) {
      return 6 + TMath::Min(m1, m2);
    }
    return 9;
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
