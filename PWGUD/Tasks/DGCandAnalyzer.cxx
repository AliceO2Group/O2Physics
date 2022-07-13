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
// \brief Analyses reduced tables (DGCandidates, DGTracks) of DG candidates produced with DGCandProducer
//
//     options:
//           anaPars.mNCombine(2)
//           anaPars.mTPCnSigmas(120, 0.)
//
//           mTPCnSigmas contains 10 blocks (particles) of 12 elements:
//              0: PID
//              1: sign
//           2, 3: min/max nsigma for e
//           4, 5: min/max nsigma for pi
//           6, 7: min/max nsigma for mu
//           8, 9: min/max nsigma for Ka
//          10,11: min/max nsigma for Pr
//          In test for particle with PID it is required: min < nsigma < max
//          In test for all other particles it is required: nsigam < min || nsigam > max
//
//     usage: copts="--configuration json://DGCandAnalyzerConfig.json -b"
//
//           o2-analysis-ud-dgcand-analyzer $copts > DGCandAnalyzer.log
//
// \author Paul Buehler, paul.buehler@oeaw.ac.at
// \since  06.06.2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Common/DataModel/PIDResponse.h"
#include "EventFiltering/PWGUD/DGCutparHolder.h"
#include "PWGUD/DataModel/DGCandidates.h"
#include "PWGUD/Core/DGPIDSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // get a DGCutparHolder and DGAnaparHolder
  DGCutparHolder diffCuts = DGCutparHolder();
  MutableConfigurable<DGCutparHolder> DGCuts{"DGCuts", {}, "DG event cuts"};
  DGAnaparHolder anaPars = DGAnaparHolder();
  MutableConfigurable<DGAnaparHolder> DGPars{"AnaPars", {}, "Analysis parameters"};

  // PID selector
  DGPIDSelector pidsel = DGPIDSelector();

  // define histograms
  HistogramRegistry registry{
    "registry",
    {
      {"nIVMs", "#nIVMs", {HistType::kTH1F, {{36, -0.5, 35.5}}}},
      {"IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      {"IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
    }};

  void init(InitContext&)
  {
    diffCuts = (DGCutparHolder)DGCuts;
    anaPars = (DGAnaparHolder)DGPars;
    pidsel.init(anaPars);
  }

  void process(aod::DGCandidate const& dgcand, aod::DGTracks const& dgtracks)
  {

    // skip events with too few/many tracks
    if (dgcand.numContrib() < diffCuts.minNTracks() || dgcand.numContrib() > diffCuts.maxNTracks()) {
      return;
    }

    // skip events with out-of-range net charge
    if (dgcand.netCharge() < diffCuts.minNetCharge() || dgcand.netCharge() > diffCuts.maxNetCharge()) {
      return;
    }

    // skip events with out-of-range rgtrwTOF
    if (dgcand.rgtrwTOF() < diffCuts.minRgtrwTOF()) {
      return;
    }

    // find track combinations which are compatible with anaPars.TPCnSigmas()
    auto nIVMs = pidsel.computeIVMs(anaPars.nCombine(), dgtracks);

    // update histograms
    registry.get<TH1>(HIST("nIVMs"))->Fill(nIVMs, 1.);
    for (auto ivm : pidsel.IVMs()) {
      registry.get<TH2>(HIST("IVMptSysDG"))->Fill(ivm.M(), ivm.Perp());
      for (auto ind : ivm.trkinds()) {
        auto track = dgtracks.rawIteratorAt(ind);
        registry.get<TH2>(HIST("IVMptTrkDG"))->Fill(ivm.M(), track.pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DGCandAnalyzer>(cfgc, TaskName{"dgcandanalyzer"}),
  };
}
