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
// \Single Gap Event Analyzer for (anti)deuteron production
// \author Marika Rasa, marika.rasa@cern.ch
// \since  September 2026

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

#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct SGDeuteronSpectra {
  // UPC cuts
  SGSelector sgSelector;
  Configurable<float> FV0_cut{"FV0", 50., "FV0A threshold"};
  Configurable<float> ZDC_cut{"ZDC", 10., "ZDC threshold"};
  Configurable<float> FT0A_cut{"FT0A", 150., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> FDDA_cut{"FDDA", 10000., "FDDA threshold"};
  Configurable<float> FDDC_cut{"FDDC", 10000., "FDDC threshold"};

  // Track cuts
  Configurable<float> PV_cut{"PV_cut", 1.0, "Use Only PV tracks"};
  Configurable<float> dcaZ_cut{"dcaZ_cut", 2.0, "dcaZ cut"};
  Configurable<float> dcaXY_cut{"dcaXY_cut", 0.0, "dcaXY cut (0 for Pt-function)"};
  Configurable<float> tpcChi2_cut{"tpcChi2_cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindable_cut{"tpcNClsFindable_cut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2_cut{"itsChi2_cut", 36, "Max itsChi2NCl"};
  Configurable<float> eta_cut{"eta_cut", 0.9, "Track Pseudorapidity"};
  Configurable<float> pt_cut{"pt_cut", 0.1, "Track Pt"};

  // Configurable axis for histograms
  ConfigurableAxis ptAxis{"ptAxis", {200, 0.0, 10.0}, "p_{T}"};
  ConfigurableAxis nsigmaAxis{"nSigmaAxis", {800, -20.0, 20.0}, "nSigma axis for TPC and TOF"};

  // Initialize histogram registry
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec axispt{ptAxis, "p_{T}"};
    const AxisSpec axistpc{nsigmaAxis, "n#sigma_{TPC}"};
    const AxisSpec axistof{nsigmaAxis, "n#sigma_{TOF}"};

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
    float FIT_cut[5] = {FV0_cut, FT0A_cut, FT0C_cut, FDDA_cut, FDDC_cut};
    int truegapSide = sgSelector.trueGap(coll, FIT_cut[0], FIT_cut[1], FIT_cut[2], ZDC_cut);
    registry.fill(HIST("collisions/TrueGapSide"), truegapSide, 1.);

    std::vector<float> parameters = {PV_cut, dcaZ_cut, dcaXY_cut, tpcChi2_cut, tpcNClsFindable_cut, itsChi2_cut, eta_cut, pt_cut};

    for (const auto& t : tracks) {
      if (trackselector(t, parameters)) {
        if (truegapSide == 0) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_GapA"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_GapA"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_GapA_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_GapA"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_GapA"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
              registry.fill(HIST("tracks/Antideut_Pt_TOF_GapA_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          }
        }

        if (truegapSide == 1) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_GapC"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_GapC"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_GapC_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_GapC"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_GapC"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
              registry.fill(HIST("tracks/Antideut_Pt_TOF_GapC_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          }
        }

        if (truegapSide == 2) {
          if (t.sign() > 0) {
            registry.fill(HIST("tracks/Deut_Pt_TPC_DoubleGap"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Deut_Pt_TOF_DoubleGap"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
              registry.fill(HIST("tracks/Deut_Pt_TOF_DoubleGap_TPCpresel"), t.pt(), t.tofNSigmaDe());
            }
          } else {
            registry.fill(HIST("tracks/Antideut_Pt_TPC_DoubleGap"), t.pt(), t.tpcNSigmaDe());
            registry.fill(HIST("tracks/Antideut_Pt_TOF_DoubleGap"), t.pt(), t.tofNSigmaDe());
            if (TMath::Abs(t.tpcNSigmaDe()) < 3.0) {
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
    adaptAnalysisTask<SGDeuteronSpectra>(cfgc, TaskName{"sgdeuteronspectra"}),
  };
}
