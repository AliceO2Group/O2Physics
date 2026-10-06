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
/// \file sgDeuteronSpectra.cxx
/// \brief Analysis for the (anti)deuteron production in UPC
/// \author Marika Rasa
/// \since September 2026

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/SGTrackSelector.h"
#include "PWGUD/DataModel/UDTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SgDeuteronSpectra {
  // UPC cuts
  SGSelector sgSelector;
  Configurable<float> fv0Cut{"fv0Cut", 50., "FV0A threshold"};
  Configurable<float> zdcCut{"zdcCut", 10., "ZDC threshold"};
  Configurable<float> ft0aCut{"ft0aCut", 150., "FT0A threshold"};
  Configurable<float> ft0cCut{"ft0cCut", 50., "FT0C threshold"};
  Configurable<float> fddaCut{"fddaCut", 10000., "FDDA threshold"};
  Configurable<float> fddcCut{"fddcCut", 10000., "FDDC threshold"};

  // Track cuts
  Configurable<float> pvCut{"pvCut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZCut{"dcaZCut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXYCut{"dcaXYCut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> ptCut{"ptCut", 0.1, "Track Pt"};
  Configurable<float> nsigmatpcPreselCut{"nsigmatpcPreselCut", 3.0, "nsigma TPC preselection cut for TOF analysis"};

  // Configurable axis for histograms
  ConfigurableAxis ptAxis{"ptAxis", {200, 0.0, 10.0}, "p_{T}"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {800, -20.0, 20.0}, "nSigma axis for TPC and TOF"};

  // Initialize histogram registry
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec axispt{ptAxis, "p_{T}"};
    const AxisSpec axistpc{nSigmaAxis, "n#sigma_{TPC}"};
    const AxisSpec axistof{nSigmaAxis, "n#sigma_{TOF}"};

    // Collision histograms
    registry.add("collisions/GapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{3, -0.5, 2.5}}});
    registry.add("collisions/TrueGapSide", "Gap Side: A, C, A+C", {HistType::kTH1F, {{4, -1.5, 2.5}}});
    registry.add("tracks/Deut_Pt_TPC_GapA", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Deut_Pt_TOF_GapA", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Deut_Pt_TOF_GapA_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TPC_GapA", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Antideut_Pt_TOF_GapA", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TOF_GapA_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Deut_Pt_TPC_GapC", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Deut_Pt_TOF_GapC", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Deut_Pt_TOF_GapC_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TPC_GapC", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Antideut_Pt_TOF_GapC", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TOF_GapC_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Deut_Pt_TPC_DoubleGap", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Deut_Pt_TOF_DoubleGap", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Deut_Pt_TOF_DoubleGap_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TPC_DoubleGap", "", {HistType::kTH2F, {axispt, axistpc}});
    registry.add("tracks/Antideut_Pt_TOF_DoubleGap", "", {HistType::kTH2F, {axispt, axistof}});
    registry.add("tracks/Antideut_Pt_TOF_DoubleGap_TPCpresel", "", {HistType::kTH2F, {axispt, axistof}});
  }

  // Define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::SGCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksPIDExtra, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    registry.fill(HIST("collisions/GapSide"), coll.gapSide(), 1.);
    std::vector<float> fitCut = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    int truegapSide = sgSelector.trueGap(coll, fitCut[0], fitCut[1], fitCut[2], zdcCut);
    int gapA = 0;
    int gapC = 1;
    int doubleGap = 2;
    registry.fill(HIST("collisions/TrueGapSide"), truegapSide, 1.);

    std::vector<float> parameters = {pvCut, dcaZCut, dcaXYCut, tpcChi2Cut, tpcNClsFindableCut, itsChi2Cut, etaCut, ptCut};

    for (const auto& t : tracks) {
      if (trackselector(t, parameters) != 0) {
        if (truegapSide == gapA) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_GapA"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_GapA"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_GapA_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_GapA"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_GapA"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Antideut_Pt_TOF_GapA_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          }
        }

        if (truegapSide == gapC) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_GapC"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_GapC"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_GapC_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_GapC"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_GapC"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Antideut_Pt_TOF_GapC_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          }
        }

        if (truegapSide == doubleGap) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_DoubleGap"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_DoubleGap"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_DoubleGap_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_DoubleGap"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_DoubleGap"), t.pt(), t.tofNSigmaDe());
            if (std::abs(t.tpcNSigmaDe()) < nsigmatpcPreselCut) {
              registry.fill(HIST("tracks/Antideut_Pt_TOF_DoubleGap_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SgDeuteronSpectra>(cfgc),
  };
}
