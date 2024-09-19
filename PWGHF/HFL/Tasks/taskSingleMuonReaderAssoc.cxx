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

/// \file taskSingleMuonReaderAssoc.cxx
/// \brief Task used to read the derived table produced by the tableMaker-association of DQ framework and extract observables on single muons needed for the HF-muon analysis.
/// \author Maolin Zhang <maolin.zhang@cern.ch>, CCNU

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

using MyCollisions = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra>;
using MyMcMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsLabels>;

namespace o2::aod
{
namespace single_muon
{
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(DcaXY, dcaXY, float);
DECLARE_SOA_COLUMN(DeltaPt, deltaPt, float);
DECLARE_SOA_COLUMN(Chi2, chi2, float);
} // namespace single_muon
DECLARE_SOA_TABLE(HfSingleMuon, "AOD", "SINGLEMUON", single_muon::Pt, single_muon::DcaXY, single_muon::DeltaPt, single_muon::Chi2);
} // namespace o2::aod

struct HfTaskSingleMuonReaderAssoc {
  Produces<aod::HfSingleMuon> singleMuon;

  Configurable<int> trkType{"trkType", 0, "Muon track type, valid values are 0, 1, 2, 3 and 4"};
  Configurable<float> etaMin{"etaMin", -3.6, "eta minimum value"};
  Configurable<float> etaMax{"etaMax", -2.5, "eta maximum value"};
  Configurable<float> pDcaMax{"pDcaMax", 594., "p*DCA maximum value"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5, "R at absorber end maximum value"};
  Configurable<float> rAbsMin{"rAbsMin", 26.5, "R at absorber end minimum value"};
  Configurable<bool> fillMcHist{"fillMcHist", false, "fill MC-related histograms"};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    AxisSpec axisPt{200, 0., 100., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{100, -4., -2., "#it{#eta}"};
    AxisSpec axisDCA{2000, 0., 2., "#it{DCA}_{xy} (cm)"};
    AxisSpec axisChi2MatchMCHMFT{100, 0., 100., "MCH-MFT matching #chi^{2}"};
    AxisSpec axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec axisRabs{1000, 0, 100, "R at Absorber End (cm)"};
    AxisSpec axisDeltaPt{10000, -50, 50, "#Delta #it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisVtxZ{80, -20., 20., "#it{z}_{vtx} (cm)"};

    HistogramConfigSpec hTHnMu{HistType::kTHnSparseF, {axisPt, axisEta, axisDCA, axisRabs, axisSign, axisChi2MatchMCHMFT, axisDeltaPt}, 7};
    HistogramConfigSpec hVtxZ{HistType::kTH1F, {axisVtxZ}};

    registry.add("hMuAfterCuts", "", hTHnMu);
    if (fillMcHist) {
      registry.add("hMuAfterCutsTrue", "", hTHnMu);
      registry.add("hMuAfterCutsFake", "", hTHnMu);
    }
    registry.add("hVtxZ", "", hVtxZ);
  }

  template <typename TCollisions, typename TMuons>
  void runMuonSelAssoc(ReducedMuonsAssoc const& assocs, TCollisions const& collisions, TMuons const&)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hVtxZ"), collision.posZ());
    }
    for (const auto& assoc : assocs) {
      auto muon = assoc.template reducedmuon_as<TMuons>();
      if (muon.trackType() != trkType) {
        continue;
      }

      const auto eta(muon.eta()), pDca(muon.pDca());
      const auto rAbs(muon.rAtAbsorberEnd());
      const auto dcaXY(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));
      const auto pt(muon.pt());
      const auto charge(muon.sign());
      const auto chi2(muon.chi2MatchMCHMFT());

      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      if ((rAbs >= rAbsMax) || (rAbs < rAbsMin)) {
        continue;
      }
      if (pDca >= pDcaMax) {
        continue;
      }

      // histograms after acceptance cuts
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<TMuons>();
        auto Dpt = muonType3.pt() - pt;

        singleMuon(pt, dcaXY, Dpt, chi2);
        registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, rAbs, charge, chi2, Dpt);
      }
    }
  }

  template <typename TCollisions, typename TMuons>
  void runMuonSelMcAssoc(ReducedMuonsAssoc const& assocs, TCollisions const& collisions, TMuons const&)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("hVtxZ"), collision.posZ());
    }
    for (const auto& assoc : assocs) {
      auto muon = assoc.template reducedmuon_as<TMuons>();
      if (muon.trackType() != trkType) {
        continue;
      }

      const auto eta(muon.eta()), pDca(muon.pDca());
      const auto rAbs(muon.rAtAbsorberEnd());
      const auto dcaXY(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));
      const auto pt(muon.pt());
      const auto charge(muon.sign());
      const auto chi2(muon.chi2MatchMCHMFT());

      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      if ((rAbs >= rAbsMax) || (rAbs < rAbsMin)) {
        continue;
      }
      if (pDca >= pDcaMax) {
        continue;
      }

      // histograms after acceptance cuts
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<TMuons>();
        auto Dpt = muonType3.pt() - pt;

        singleMuon(pt, dcaXY, Dpt, chi2);
        registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, rAbs, charge, chi2, Dpt);
        if (muon.mcMask() == 0) {
          registry.fill(HIST("hMuAfterCutsTrue"), pt, eta, dcaXY, rAbs, charge, chi2, Dpt);
        }
        if (muon.mcMask() == 128) {
          registry.fill(HIST("hMuAfterCutsFake"), pt, eta, dcaXY, rAbs, charge, chi2, Dpt);
        }
      }
    }
  }

  void processMuon(ReducedMuonsAssoc const& assocs, MyCollisions const& collisions, MyMuons const& muons)
  {
    runMuonSelAssoc(assocs, collisions, muons);
  }
  PROCESS_SWITCH(HfTaskSingleMuonReaderAssoc, processMuon, "run muon selection with real data", true);

  void processMuonMc(ReducedMuonsAssoc const& assocs, MyCollisions const& collisions, MyMcMuons const& muons)
  {
    runMuonSelMcAssoc(assocs, collisions, muons);
  }
  PROCESS_SWITCH(HfTaskSingleMuonReaderAssoc, processMuonMc, "run muon selection with MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSingleMuonReaderAssoc>(cfgc)};
}
