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

/// \file jV02DmesonToy.cxx
/// \brief One timeframe of synthetic collisions for JV02DmesonAnalysis.
/// \author Dong Jo Kim (djkim@jyu.fi)
///
/// 200 events at centrality 25. N_ch = 20.
/// Even events: XY = 0.08, 6 tracks in pT [1.0, 1.2), 4 D0 at mass 1.864 and pT 3.5.
/// Odd events:  XY = 0.02, 2 tracks in pT [1.0, 1.2), 1 D0 at the same mass and pT.
/// The remaining charged tracks sit at pT 0.3.
/// Expected charged-hadron v02 in the 1.1 GeV/c bin is 0.3.
/// Expected D0 v02, with an empty sideband, is 0.36.

#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AlgorithmSpec.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/ControlService.h>
#include <Framework/DataProcessorSpec.h>
#include <Framework/TableBuilder.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

struct JV02DmesonToy {
  Produces<aod::BCs> bcs;
  Produces<aod::Collisions> collisions;
  Produces<aod::EvSels> evsels;
  Produces<aod::CentFT0Cs> centsC;
  Produces<aod::CentFT0Ms> centsM;
  Produces<aod::Qvectors> qvectors;
  Produces<aod::CFCollisions> cfCollisions;
  Produces<aod::CFCollRefs> cfRefs;
  Produces<aod::CFTracks> cfTracks;
  Produces<aod::CF2ProngTracks> cfCandidates;

  void fill()
  {
    constexpr int NEvents = 200;
    constexpr int NDetectors = 7;
    constexpr int NCorrections = 4;
    constexpr int FullyCorrected = 3;
    constexpr int Ft0C = 0;
    constexpr int Ft0A = 1;
    constexpr int EventPair = 2;
    constexpr float Centrality = 25.f;
    constexpr float Nch = 20.f;
    const uint64_t selection =
      (uint64_t{1} << aod::evsel::kNoSameBunchPileup) |
      (uint64_t{1} << aod::evsel::kIsGoodZvtxFT0vsPV) |
      (uint64_t{1} << aod::evsel::kNoTimeFrameBorder) |
      (uint64_t{1} << aod::evsel::kNoITSROFrameBorder);

    std::vector<float> re(NDetectors * NCorrections, 0.f);
    std::vector<float> im(NDetectors * NCorrections, 0.f);
    std::vector<float> amp(NDetectors, 0.f);
    amp[Ft0A] = 1.f;
    amp[Ft0C] = 1.f;

    for (int ievt = 0; ievt < NEvents; ++ievt) {
      const bool high = ievt % EventPair == 0;
      const float xy = high ? 0.08f : 0.02f;
      const float q = std::sqrt(xy);
      re[Ft0C * NCorrections + FullyCorrected] = q;
      re[Ft0A * NCorrections + FullyCorrected] = q;

      bcs(544116, static_cast<uint64_t>(ievt), uint64_t{0}, uint64_t{0});
      collisions(ievt, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, uint16_t{0}, 1.f, uint16_t{20}, 0.f, 0.f);
      evsels(uint32_t{0}, selection, uint32_t{0}, false, true, ievt, -1, -1, -1, -1, 100, 0.f, 0.f);
      centsC(Centrality);
      centsM(Centrality);
      qvectors(Centrality, true, re, im, amp);
      cfCollisions(544116, 0.f, Nch, uint64_t{0});
      cfRefs(ievt);

      const int nSignalTracks = high ? 6 : 2;
      for (int itrk = 0; itrk < static_cast<int>(Nch); ++itrk) {
        const float pt = itrk < nSignalTracks ? 1.1f : 0.3f;
        const int8_t sign = itrk % EventPair == 0 ? int8_t{1} : int8_t{-1};
        cfTracks(ievt, pt, 0.f, 0.f, sign, uint8_t{1});
      }
      const int nD0 = high ? 4 : 1;
      for (int icand = 0; icand < nD0; ++icand) {
        cfCandidates(ievt, 0, 1, 3.5f, 0.f, 0.f, static_cast<float>(o2::constants::physics::MassD0), static_cast<uint8_t>(aod::cf2prongtrack::D0ToPiK));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    DataProcessorSpec{
      "jv02-dmeson-toy",
      Inputs{},
      Outputs{
        OutputForTable<aod::BCs>::spec(),
        OutputForTable<aod::Collisions>::spec(),
        OutputForTable<aod::EvSels>::spec(),
        OutputForTable<aod::CentFT0Cs>::spec(),
        OutputForTable<aod::CentFT0Ms>::spec(),
        OutputForTable<aod::Qvectors>::spec(),
        OutputForTable<aod::CFCollisions>::spec(),
        OutputForTable<aod::CFCollRefs>::spec(),
        OutputForTable<aod::CFTracks>::spec(),
        OutputForTable<aod::CF2ProngTracks>::spec(),
        OutputSpec{"TFN", "TFNumber"},
        OutputSpec{"TFF", "TFFilename"},
      },
      AlgorithmSpec{adaptStateless([](ProcessingContext& pc) {
        JV02DmesonToy toy;
        toy.bcs.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::BCs>::ref())));
        toy.collisions.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::Collisions>::ref())));
        toy.evsels.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::EvSels>::ref())));
        toy.centsC.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CentFT0Cs>::ref())));
        toy.centsM.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CentFT0Ms>::ref())));
        toy.qvectors.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::Qvectors>::ref())));
        toy.cfCollisions.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CFCollisions>::ref())));
        toy.cfRefs.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CFCollRefs>::ref())));
        toy.cfTracks.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CFTracks>::ref())));
        toy.cfCandidates.resetCursor(std::move(pc.outputs().make<TableBuilder>(OutputForTable<aod::CF2ProngTracks>::ref())));
        toy.fill();
        pc.outputs().make<uint64_t>(Output{"TFN", "TFNumber"}) = 1;
        pc.outputs().make<std::string>(Output{"TFF", "TFFilename"}) = std::string{"toy"};
        pc.services().get<ControlService>().endOfStream();
        pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
      })},
    },
  };
}
