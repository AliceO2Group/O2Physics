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

/// \file taskSingleMuon.cxx
/// \brief Task used to extract the observables on single muons needed for the HF-muon analysis.
/// \author Maolin Zhang <maolin.zhang@cern.ch>, CCNU

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <sys/types.h>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace ambiguous_mft
{ // for MA2T
DECLARE_SOA_INDEX_COLUMN(MFTTrack, track);
DECLARE_SOA_INDEX_COLUMN(AmbiguousMFTTrack, ambMftTrack);
} // namespace ambiguous_mft

DECLARE_SOA_INDEX_TABLE_USER(MA2T, MFTTracks, "MA2T", ambiguous_mft::MFTTrackId, ambiguous_mft::AmbiguousMFTTrackId);
} // namespace o2::aod

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMcMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MFTTracksExtra = soa::Join<aod::MFTTracks, aod::MA2T>; // MFT track + index of ambiguous track

struct HfTaskSingleMuonSelectionAmbiguousMftIndexBuilder {
  // build the index table MA2T
  Builds<aod::MA2T> idx;

  void init(InitContext const&) {}
  // Additionnal methods provided:
  // mfttrack.has_ambMftTrack()
  // mfttrack.ambMftTrack()
};

struct HfTaskSingleMuon {
  Configurable<uint> trkType{"trkType", 0u, "Muon track type, valid values are 0, 1, 2, 3 and 4"};
  Configurable<uint> mcMaskSelection{"mcMaskSelection", 0u, "McMask for correct match, valid values are 0 and 128"};
  Configurable<float> etaMin{"etaMin", -3.6, "eta minimum value"};
  Configurable<float> etaMax{"etaMax", -2.5, "eta maximum value"};
  Configurable<float> pDcaMin{"pDcaMin", 324., "p*DCA maximum value for small Rabs"};
  Configurable<float> pDcaMax{"pDcaMax", 594., "p*DCA maximum value for large Rabs"};
  Configurable<float> rAbsMin{"rAbsMin", 17.6, "R at absorber end minimum value"};
  Configurable<float> rAbsMax{"rAbsMax", 89.5, "R at absorber end maximum value"};
  Configurable<float> rAbsMid{"rAbsMid", 26.5, "R at absorber end split point for different p*DCA selections"};
  Configurable<float> zVtx{"zVtx", 10., "Z edge of primary vertex [cm]"};
  Configurable<bool> fillMcHist{"fillMcHist", false, "fill MC-related muon histograms"};
  Configurable<bool> reduceAmbMft{"reduceAmbMft", false, "reduce ambiguous MFT tracks"};
  Configurable<bool> reduceOrphMft{"reduceOrphMft", true, "reduce orphan MFT tracks"};

  Filter posZfilter = nabs(aod::collision::posZ) < zVtx;
  Filter mcMaskFilter = aod::mcfwdtracklabel::mcMask == mcMaskSelection;

  o2::framework::HistogramRegistry registry{
    "registry",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  void init(InitContext&)
  {
    AxisSpec const axisPt{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisEta{250, -5., 0., "#it{#eta}"};
    AxisSpec const axisDCA{500, 0., 5., "#it{DCA}_{xy} (cm)"};
    AxisSpec const axisChi2MatchMCHMFT{100, 0., 100., "MCH-MFT matching #chi^{2}"};
    AxisSpec const axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec const axisPDca{100000, 0, 100000, "#it{p} #times DCA (GeV/#it{c} * cm)"};
    AxisSpec const axisVtxZ{80, -20., 20., "#it{z}_{vtx} (cm)"};
    AxisSpec const axisDCAx{1000, -5., 5., "#it{DCA}_{x or y} (cm)"};
    AxisSpec const axisPtDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};
    AxisSpec const axisEtaDif{200, -2., 2., "#it{#eta} diff"};
    AxisSpec const axisDeltaPt{60, -30, 30, "#Delta #it{p}_{T} (GeV/#it{c})"};

    HistogramConfigSpec const hTHnMu{HistType::kTHnSparseF, {axisPt, axisEta, axisDCA, axisPDca, axisSign, axisChi2MatchMCHMFT}, 6};
    HistogramConfigSpec const h2PtMc{HistType::kTH2F, {axisPt, axisPtDif}};
    HistogramConfigSpec const h2EtaMc{HistType::kTH2F, {axisEta, axisEtaDif}};
    HistogramConfigSpec const h2DCA{HistType::kTH2F, {axisDCAx, axisDCAx}};
    HistogramConfigSpec const h3DeltaPt{HistType::kTH3F, {axisPt, axisEta, axisDeltaPt}};
    HistogramConfigSpec const hVtxZ{HistType::kTH1F, {axisVtxZ}};

    registry.add("hMuBeforeCuts", "", hTHnMu);
    registry.add("hMuAfterCuts", "", hTHnMu);
    registry.add("h3DeltaPtBeforeCuts", "", h3DeltaPt);
    registry.add("h3DeltaPtAfterCuts", "", h3DeltaPt);
    registry.add("h2DCA", "", h2DCA);
    registry.add("hVtxZ", "", hVtxZ);
    if (fillMcHist) {
      registry.add("h2PtMcBeforeCuts", "", h2PtMc);
      registry.add("h2PtMcAfterCuts", "", h2PtMc);
      registry.add("h2EtaMcBeforeCuts", "", h2EtaMc);
      registry.add("h2EtaMcAfterCuts", "", h2EtaMc);
    }
  }

  template <typename TCollision, typename TMFT, typename TMuons>
  void runMuonSel(TCollision const& collision, TMFT const& /*tracksMFT*/, TMuons const& muons)
  {
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hVtxZ"), collision.posZ());

    for (const auto& muon : muons) {
      if (muon.trackType() != trkType) {
        continue;
      }
      const auto eta(muon.eta()), pDca(muon.pDca()), rAbs(muon.rAtAbsorberEnd());
      const auto dcaXY(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));
      const auto pt(muon.pt());
      const auto charge(muon.sign());
      const auto chi2(muon.chi2MatchMCHMFT());

      if (muon.has_matchMFTTrack()) {
        auto trkMFT = muon.template matchMFTTrack_as<MFTTracksExtra>();
        if (reduceAmbMft && trkMFT.has_ambMftTrack()) {
          continue;
        }
        if (reduceOrphMft && (!reduceAmbMft) && trkMFT.has_ambMftTrack()) {
          continue;
        }
      }
      // histograms before the acceptance cuts
      registry.fill(HIST("hMuBeforeCuts"), pt, eta, dcaXY, pDca, charge, chi2);
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        registry.fill(HIST("h3DeltaPtBeforeCuts"), pt, eta, muonType3.pt() - pt);
      }
      // standard muon acceptance cuts
      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      if ((rAbs >= rAbsMax) || (rAbs < rAbsMin)) {
        continue;
      }
      if ((rAbs < rAbsMid) && (pDca >= pDcaMin)) {
        continue;
      }
      if ((rAbs >= rAbsMid) && (pDca >= pDcaMax)) {
        continue;
      }
      if ((muon.chi2() >= 1e6) || (muon.chi2() < 0)) {
        continue;
      }
      // histograms after acceptance cuts
      registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, pDca, charge, chi2);
      registry.fill(HIST("h2DCA"), muon.fwdDcaX(), muon.fwdDcaY());
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        registry.fill(HIST("h3DeltaPtAfterCuts"), pt, eta, muonType3.pt() - pt);
      }
    }
  }

  template <typename TCollision, typename TMFT, typename TMuons, typename TMC>
  void runMuonSelMc(TCollision const& collision, TMFT const& /*tracksMFT*/, TMuons const& muons, TMC const& /*mc*/)
  {
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hVtxZ"), collision.posZ());

    for (const auto& muon : muons) {
      if (muon.trackType() != trkType) {
        continue;
      }
      const auto eta(muon.eta()), pDca(muon.pDca());
      const auto dcaXY(RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY()));
      const auto pt(muon.pt());
      const auto charge(muon.sign());
      const auto chi2(muon.chi2MatchMCHMFT());

      if (muon.has_matchMFTTrack()) {
        auto trkMFT = muon.template matchMFTTrack_as<MFTTracksExtra>();
        if (reduceAmbMft && trkMFT.has_ambMftTrack()) {
          continue;
        }
        if (reduceOrphMft && (!reduceAmbMft) && trkMFT.has_ambMftTrack()) {
          continue;
        }
      }
      // histograms before the acceptance cuts
      registry.fill(HIST("hMuBeforeCuts"), pt, eta, dcaXY, pDca, charge, chi2);
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMcMuons>();
        registry.fill(HIST("h3DeltaPtBeforeCuts"), pt, eta, muonType3.pt() - pt);
      }
      if (fillMcHist) {
        if (muon.has_mcParticle()) {
          auto mcMuon = muon.mcParticle();
          registry.fill(HIST("h2PtMcBeforeCuts"), mcMuon.pt(), mcMuon.pt() - pt);
          registry.fill(HIST("h2EtaMcBeforeCuts"), mcMuon.eta(), mcMuon.eta() - eta);
        }
      }

      // standard muon acceptance cuts
      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      if ((muon.chi2() >= 1e6) || (muon.chi2() < 0)) {
        continue;
      }
      // histograms after acceptance cuts
      registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, pDca, charge, chi2);
      registry.fill(HIST("h2DCA"), muon.fwdDcaX(), muon.fwdDcaY());
      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMcMuons>();
        registry.fill(HIST("h3DeltaPtAfterCuts"), pt, eta, muonType3.pt() - pt);
      }
      if (fillMcHist) {
        if (muon.has_mcParticle()) {
          auto mcMuon = muon.mcParticle();
          registry.fill(HIST("h2PtMcAfterCuts"), mcMuon.pt(), mcMuon.pt() - pt);
          registry.fill(HIST("h2EtaMcAfterCuts"), mcMuon.eta(), mcMuon.eta() - eta);
        }
      }
    }
  }

  void processMuon(soa::Filtered<MyCollisions>::iterator const& collision,
                   MFTTracksExtra const& tracksMFT,
                   MyMuons const& muons)
  {
    runMuonSel(collision, tracksMFT, muons);
  }
  void processMuonMc(soa::Filtered<MyCollisions>::iterator const& collision,
                     MFTTracksExtra const& tracksMFT,
                     soa::Filtered<MyMcMuons> const& muons,
                     aod::McParticles const& mc)
  {
    runMuonSelMc(collision, tracksMFT, muons, mc);
  }

  PROCESS_SWITCH(HfTaskSingleMuon, processMuon, "run muon selection with real data", true);
  PROCESS_SWITCH(HfTaskSingleMuon, processMuonMc, "run muon selection with MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskSingleMuonSelectionAmbiguousMftIndexBuilder>(cfgc),
    adaptAnalysisTask<HfTaskSingleMuon>(cfgc),
  };
}
