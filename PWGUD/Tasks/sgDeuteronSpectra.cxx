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
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/RecoDecay.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

constexpr std::array<std::array<float, 4>, 1> DCAxyDefault = {{{0.0105, 0.035, 1.1, 0.001}}};
constexpr std::array<std::array<float, 4>, 1> DCAzDefault = {{{0.0105, 0.035, 1.1, 0.001}}};
const std::vector<std::string> dcaParNames{"par0", "par1", "par2", "mean"};

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
  Configurable<LabeledArray<float>> dcaxyParameters{"dcaxyParameters", {DCAxyDefault[0].data(), 4, dcaParNames}, "DCAxy parametrization"};
  Configurable<LabeledArray<float>> dcazParameters{"dcazParameters", {DCAzDefault[0].data(), 4, dcaParNames}, "DCAz parametrization"};
  Configurable<float> nsigmaDCAxyCut{"nsigmaDCAxyCut", 5.0, "n#sigma cut for DCA_{XY}"};
  Configurable<float> nsigmaDCAzCut{"nsigmaDCAzCut", 5.0, "n#sigma cut for DCA_{Z}"};
  Configurable<float> tpcChi2Cut{"tpcChi2Cut", 4, "Max tpcChi2NCl"};
  Configurable<float> tpcNClsFindableCut{"tpcNClsFindableCut", 70, "Min tpcNClsFindable"};
  Configurable<float> itsChi2Cut{"itsChi2Cut", 36, "Max itsChi2NCl"};
  Configurable<float> etaCut{"etaCut", 0.9, "Track Pseudorapidity"};
  Configurable<float> ptCut{"ptCut", 0.1, "Track Pt"};
  Configurable<float> nsigmatpcPreselCut{"nsigmatpcPreselCut", 3.0, "nsigma TPC preselection cut for TOF analysis"};
  Configurable<float> occupancyCut{"occupancyCut", 1000, "Occupancy cut"};
  Configurable<float> vertexZCut{"vertexZCut", 10, "VertexZ cut"};
  Configurable<float> maxTracks{"maxTracks", 100, "Maximum number of tracks in the event"};
  Configurable<bool> useVtxFT0Cut{"useVtxFT0Cut", false, "Use zVtxFT0vPV cut in event selection"};

  // Configurable axis for histograms
  ConfigurableAxis ptAxis{"ptAxis", {200, 0.0, 10.0}, "p_{T}"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {800, -20.0, 20.0}, "nSigma axis for TPC and TOF"};
  ConfigurableAxis dcaAxis{"dcaAxis", {600, -3.0, 3.0}, "DCA axis"};

  // Initialize histogram registry
  HistogramRegistry registry{"registry", {}};

  void init(InitContext&)
  {
    const AxisSpec axispt{ptAxis, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axistpc{nSigmaAxis, "n#sigma_{TPC}"};
    const AxisSpec axistof{nSigmaAxis, "n#sigma_{TOF}"};
    const AxisSpec axisdcaxy{dcaAxis, "DCA_{XY} (cm)"};
    const AxisSpec axisdcaz{dcaAxis, "DCA_{Z} (cm)"};

    // Collision histograms
    if (doprocessData || doprocessDCA) {
      registry.add("collisions/selectionCuts", "", {HistType::kTH1D, {{9, 0, 9}}});
      registry.add("collisions/VertexZ_before", "", {HistType::kTH1D, {{800, -20, 20}}});
      registry.add("collisions/VertexZ_after", "", {HistType::kTH1D, {{800, -20, 20}}});
      registry.add("collisions/GapSide", "Gap Side: Gap A, Gap C, Double Gap", {HistType::kTH1D, {{3, -0.5, 2.5}}});
      registry.add("collisions/TrueGapSide", "Gap Side: No Gap, Gap A, Gap C, Double Gap, No UPC, Trk Out Range, Bad Double Gap", {HistType::kTH1D, {{7, -1.5, 5.5}}});
    }
    if (doprocessData) {
      registry.add("tracks/DCAxy_Pt_postracks_postcuts", "DCAxy vs Pt for positive tracks after cuts", {HistType::kTH2D, {axispt, axisdcaxy}});
      registry.add("tracks/DCAxy_Pt_negtracks_postcuts", "DCAxy vs Pt for negative tracks after cuts", {HistType::kTH2D, {axispt, axisdcaxy}});
      registry.add("tracks/DCAz_Pt_postracks_postcuts", "DCAz vs Pt for positive tracks after cuts", {HistType::kTH2D, {axispt, axisdcaz}});
      registry.add("tracks/DCAz_Pt_negtracks_postcuts", "DCAz vs Pt for negative tracks after cuts", {HistType::kTH2D, {axispt, axisdcaz}});
      registry.add("tracks/Deut_Pt_TPC_GapA", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Deut_Pt_TOF_GapA", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Deut_Pt_TOF_GapA_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TPC_GapA", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Antideut_Pt_TOF_GapA", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TOF_GapA_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Deut_Pt_TPC_GapC", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Deut_Pt_TOF_GapC", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Deut_Pt_TOF_GapC_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TPC_GapC", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Antideut_Pt_TOF_GapC", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TOF_GapC_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Deut_Pt_TPC_DoubleGap", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Deut_Pt_TOF_DoubleGap", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Deut_Pt_TOF_DoubleGap_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TPC_DoubleGap", "", {HistType::kTH2D, {axispt, axistpc}});
      registry.add("tracks/Antideut_Pt_TOF_DoubleGap", "", {HistType::kTH2D, {axispt, axistof}});
      registry.add("tracks/Antideut_Pt_TOF_DoubleGap_TPCpresel", "", {HistType::kTH2D, {axispt, axistof}});
    }
    if (doprocessDCA) {
      registry.add("tracks/DCAxy_Pt_postracks_precuts", "DCAxy vs Pt for positive tracks before cuts", {HistType::kTH2D, {axispt, axisdcaxy}});
      registry.add("tracks/DCAxy_Pt_negtracks_precuts", "DCAxy vs Pt for negative tracks before cuts", {HistType::kTH2D, {axispt, axisdcaxy}});
      registry.add("tracks/DCAz_Pt_postracks_precuts", "DCAz vs Pt for positive tracks before cuts", {HistType::kTH2D, {axispt, axisdcaz}});
      registry.add("tracks/DCAz_Pt_negtracks_precuts", "DCAz vs Pt for negative tracks before cuts", {HistType::kTH2D, {axispt, axisdcaz}});
    }
  }

  // Define data types
  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions, aod::UDZdcsReduced>; // UDCollisions
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksPIDExtra, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  // define function for event selection
  bool passEventSelection(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    registry.fill(HIST("collisions/VertexZ_before"), coll.posZ());
    registry.fill(HIST("collisions/selectionCuts"), 0.5);
    // applying event cuts
    if (coll.tfb() == 0) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 1.5);
    if (coll.itsROFb() == 0) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 2.5);
    if (coll.sbp() == 0) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 3.5);
    if (useVtxFT0Cut) {
      if (coll.zVtxFT0vPV() == 0) {
        return false;
      }
      registry.fill(HIST("collisions/selectionCuts"), 4.5);
    } else {
      registry.fill(HIST("collisions/selectionCuts"), 4.5);
    }
    if (coll.vtxITSTPC() == 0) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 5.5);
    if (std::abs(coll.posZ()) > vertexZCut) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 6.5);
    if (coll.occupancyInTime() > occupancyCut) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 7.5);
    if (tracks.size() > maxTracks) {
      return false;
    }
    registry.fill(HIST("collisions/selectionCuts"), 8.5);
    registry.fill(HIST("collisions/VertexZ_after"), coll.posZ());
    return true;
  }

  // fill gapside after quality cuts
  int fillGapSide(UDCollisionFull const& coll)
  {
    registry.fill(HIST("collisions/GapSide"), coll.gapSide(), 1.);
    std::vector<float> fitCut = {fv0Cut, ft0aCut, ft0cCut, fddaCut, fddcCut};
    int truegapSide = sgSelector.trueGap(coll, fitCut[0], fitCut[1], fitCut[2], zdcCut);
    registry.fill(HIST("collisions/TrueGapSide"), truegapSide, 1.);

    return truegapSide;
  }

  float computeNsigmaDCA(const float pt, const float dca, const char* dcaType = "xy")
  {
    const int nParams = 4;
    std::array<float, 4> params{{0., 0., 0., 0.}};
    if (std::strcmp(dcaType, "xy") == 0) {
      for (int i = 0; i < nParams; ++i) {
        params[i] = dcaxyParameters->get(i);
      }
    } else if (std::strcmp(dcaType, "z") == 0) {
      for (int i = 0; i < nParams; ++i) {
        params[i] = dcazParameters->get(i);
      }
    } else {
      LOG(error) << "Invalid dcaType. Use 'xy' or 'z'.";
    }
    const float sigma = params[0] + params[1] / std::pow(pt, params[2]);
    const float mean = params[3];
    return (dca - mean) / sigma;
  }

  // define function for track selection
  bool passTrackSelection(UDTracksFull::iterator const& tracks)
  {
    if (!tracks.isPVContributor()) {
      return false;
    }
    if (tracks.tpcChi2NCl() > tpcChi2Cut) {
      return false;
    }
    if (tracks.tpcNClsFindable() < tpcNClsFindableCut) {
      return false;
    }
    if (tracks.itsChi2NCl() > itsChi2Cut) {
      return false;
    }
    auto momentum = std::array<double, 3>{tracks.px(), tracks.py(), tracks.pz()};
    if (std::abs(RecoDecay::eta(momentum)) > etaCut) {
      return false;
    }
    if (tracks.pt() < ptCut) {
      return false;
    }
    const float nsigmaDCAxy = computeNsigmaDCA(tracks.pt(), tracks.dcaXY(), "xy");
    const float nsigmaDCAz = computeNsigmaDCA(tracks.pt(), tracks.dcaZ(), "z");
    if (std::abs(nsigmaDCAxy) > nsigmaDCAxyCut) {
      return false;
    }
    if (std::abs(nsigmaDCAz) > nsigmaDCAzCut) {
      return false;
    }
    return true;
  }

  void processData(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    if (!passEventSelection(coll, tracks)) {
      return;
    }
    auto truegapSide = fillGapSide(coll);

    for (const auto& t : tracks) {
      if (passTrackSelection(t)) {
        if (t.sign() > 0) {
          registry.fill(HIST("tracks/DCAxy_Pt_postracks_postcuts"), t.pt(), t.dcaXY());
          registry.fill(HIST("tracks/DCAz_Pt_postracks_postcuts"), t.pt(), t.dcaZ());
        } else {
          registry.fill(HIST("tracks/DCAxy_Pt_negtracks_postcuts"), t.pt(), t.dcaXY());
          registry.fill(HIST("tracks/DCAz_Pt_negtracks_postcuts"), t.pt(), t.dcaZ());
        }
        if (truegapSide == o2::aod::sgselector::SingleGapA) {
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

        if (truegapSide == o2::aod::sgselector::SingleGapC) {
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

        if (truegapSide == o2::aod::sgselector::DoubleGap) {
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
  PROCESS_SWITCH(SgDeuteronSpectra, processData, "Process data for spectra analysis", true)

  void processDCA(UDCollisionFull const& coll, UDTracksFull const& tracks)
  {
    passEventSelection(coll, tracks);
    for (const auto& t : tracks) {
      if (t.sign() > 0) {
        registry.fill(HIST("tracks/DCAxy_Pt_postracks_precuts"), t.pt(), t.dcaXY());
        registry.fill(HIST("tracks/DCAz_Pt_postracks_precuts"), t.pt(), t.dcaZ());
      } else {
        registry.fill(HIST("tracks/DCAxy_Pt_negtracks_precuts"), t.pt(), t.dcaXY());
        registry.fill(HIST("tracks/DCAz_Pt_negtracks_precuts"), t.pt(), t.dcaZ());
      }
    }
  }
  PROCESS_SWITCH(SgDeuteronSpectra, processDCA, "Process data for DCA fits", false)
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<SgDeuteronSpectra>(cfgc),
  };
}
