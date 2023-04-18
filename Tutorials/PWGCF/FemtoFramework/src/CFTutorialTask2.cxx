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
///
/// \brief Femtodream Tutorial 2
/// \author Luca Barioglio, Anton Riedel

// O2 includes
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// STEP 2
// Example task illustrating how to use Configurables and Filters

namespace o2::aod
{
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyTracks =
  soa::Join<aod::FullTracks, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
            aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe>;
using MyCollision = MyCollisions::iterator;
using MyTrack = MyTracks::iterator;
} // namespace o2::aod

struct CFTutorialTask2 {
  HistogramRegistry histos{
    "Histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  // TODO
  // define configurables for zvtx, eta, pt and nsigma_TPC cut

  // TODO
  // define filters for collisions and tracks

  // TODO
  // apply filters to collisions and tracks

  // Equivalent of the AliRoot task UserCreateOutputObjects
  void init(o2::framework::InitContext&)
  {
    // Define your axes
    // Constant bin width axis
    AxisSpec vtxZAxis = {100, -20, 20};
    // Variable bin width axis
    std::vector<double> ptBinning = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1,
                                     1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0,
                                     2.2, 2.4, 2.8, 3.2, 3.6, 4.};
    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};

    // Add histograms to histogram manager (as in the output object of in AliPhysics)
    histos.add("hZvtx", ";Z (cm)", kTH1F, {vtxZAxis});
    histos.add("hP", ";#it{p} (GeV/#it{c})", kTH1F, {{35, 0.5, 4.}});
    histos.add("hEta", ";#it{p} (GeV/#it{c})", kTH1F, {{100, -1.5, 1.5}});
    histos.add("hPt", ";#it{p}_{T} (GeV/#it{c})", kTH1F, {ptAxis});
    histos.add("hNsigmaTPC", ";#it{p} (GeV/#it{c}); n#sigma_{TPC}^{proton}", kTH2F,
               {{35, 0.5, 4.}, {100, -5., 5.}});
  }

  // TODO
  // use filtered collisions and track tables
  void process(aod::Collision const& coll, aod::Tracks const& inputTracks)
  {
    histos.fill(HIST("hZvtx"), coll.posZ());

    for (auto track : inputTracks) {
      // TODO
      // add cut on TPCNSigmaPr
      // Note: TPCNSigmaPr is a dynamic column and as such it is not compatible with filters

      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hP"), track.p());
      histos.fill(HIST("hPt"), track.pt());
      histos.fill(HIST("hEta"), track.eta());
      // histos.fill(HIST("hNsigmaTPCP"), track.tpcInnerParam(), track.tpcNSigmaPr());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Equivalent to the AddTask in AliPhysics
  WorkflowSpec workflow{adaptAnalysisTask<CFTutorialTask2>(cfgc)};
  return workflow;
}
