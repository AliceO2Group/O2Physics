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
// \brief Saves relevant information of DG candidates
//
//     options:
//           anaPars.mMinNTracks(0)
//           anaPars.mMaxNTracks(10000)
//           anaPars.mNCombine(2)
//           anaPars.mPIDinfo(120, 0.)
//
//           mPIDinfo contains 10 blocks (particles) of 12 elements:
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

#include "Common/Core/PID/PIDResponse.h"
#include "PWGUD/DataModel/DGCandidates.h"
#include "PWGUD/Core/pidSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct DGCandAnalyzer {

  // get a cutHolder
  anaparHolder anaPars = anaparHolder();
  MutableConfigurable<anaparHolder> DGPars{"AnaPars", {}, "Analysis parameters"};

  // PID selector
  pidSelector pidsel = pidSelector();
  TList* IVMs = new TList();

  HistogramRegistry registry{
    "registry",
    {
      {"IVMptSysDG", "#IVMptSysDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
      {"IVMptTrkDG", "#IVMptTrkDG", {HistType::kTH2F, {{100, 0., 5.}, {350, 0., 3.5}}}},
    }};

  void init(InitContext&)
  {
    anaPars = (anaparHolder)DGPars;
    pidsel.init(anaPars);
  }

  void process(aod::DGCandidate const& dgcand, aod::DGTracks const& dgtracks)
  {

    // skip events with too few/many tracks
    if (dgcand.numContrib() < anaPars.minNTracks() || dgcand.numContrib() > anaPars.maxNTracks()) {
      return;
    }

    // find track combinations which are compatible with anaPars.PIDinfo()
    auto nCombs = pidsel.computeIVMs(anaPars.nCombine(), dgtracks);

    // update histograms
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
