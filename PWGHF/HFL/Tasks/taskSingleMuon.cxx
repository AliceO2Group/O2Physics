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

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace ambii
{ // for MA2T
DECLARE_SOA_INDEX_COLUMN(MFTTrack, track);
DECLARE_SOA_INDEX_COLUMN(AmbiguousMFTTrack, ambMFTtrack);
} // namespace ambii

DECLARE_SOA_INDEX_TABLE_USER(MA2T, MFTTracks, "MA2T", ambii::MFTTrackId, ambii::AmbiguousMFTTrackId);
} // namespace o2::aod

using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
using MyMcMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
using MFTTracksExtra = soa::Join<aod::MFTTracks, aod::MA2T>; // MFT track + index of ambiguous track

struct ProduceMFTAmbii {
  // build the index table MA2T
  Builds<aod::MA2T> idx;
  void init(InitContext const&) {}
  // Additionnal methods provided:
  // mfttrack.has_ambMFTtrack()
  // mfttrack.ambMFTtrack()
};

struct HfTaskSingleMuonSelection {
  Configurable<uint8_t> trkType{"trkType", 0, "Muon track type, validate values are 0, 1, 2, 3 and 4"};
  Configurable<uint8_t> mcMaskSelection{"mcMaskSelection", 0, "McMask for correct match, validate values are 0 and 128"};
  Configurable<float> etaLow{"etaLow", -3.6, "Eta acceptance low edge"};
  Configurable<float> etaUp{"etaUp", -2.5, "Eta acceptance up edge"};
  Configurable<float> pDcaLow{"pDcaLow", 324., "pDCA acceptance low edge"};
  Configurable<float> pDcaUp{"pDcaUp", 594., "pDCA acceptance up edge"};
  Configurable<float> rAbsLow{"rAbsLow", 17.6, "Rabs acceptance low edge"};
  Configurable<float> rAbsUp{"rAbsUp", 89.5, "Rabs acceptance up edge"};
  Configurable<float> rAbsMid{"rAbsMid", 26.5, "Rabs acceptance middle edge"};
  Configurable<float> zVtx{"zVtx", 10., "Z edge of primary vertex [cm]"};
  Configurable<bool> fillMcHist{"fillMcHist", false, "fill MC-related muon histograms"};
  Configurable<bool> reduceAmbMFT{"reduceAmbMFT", false, "reduce ambiguous MFT tracks"};

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
    AxisSpec axisPt{200, 0., 200., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec axisEta{250, -5., 0., "#eta"};
    AxisSpec axisDCA{500, 0., 5., "DCA_{xy} (cm)"};
    AxisSpec axisChi2MatchMCHMFT{100, 0., 100., "MCH-MFT matching #chi^{2}"};
    AxisSpec axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec axisPDca{100000, 0, 100000, "p #times DCA (GeV/#it{c} * cm)"};
    AxisSpec axisVtxZ{80, -20., 20., "Z_{vtx} (cm)"};
    AxisSpec axisDCAx{1000, -5., 5., "DCA_{x or y} (cm)"};
    AxisSpec axisPtDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};
    AxisSpec axisEtaDif{200, -2., 2., "#eta diff"};
    AxisSpec axisDeltaPt{60, -30, 30, "#Delta #it{p}_{T} (GeV/#it{c})"};

    HistogramConfigSpec hTHnMu{HistType::kTHnSparseF, {axisPt, axisEta, axisDCA, axisPDca, axisSign, axisChi2MatchMCHMFT}, 6};
    HistogramConfigSpec h2PtMc{HistType::kTH2F, {axisPt, axisPtDif}};
    HistogramConfigSpec h2EtaMc{HistType::kTH2F, {axisEta, axisEtaDif}};
    HistogramConfigSpec h2DCA{HistType::kTH2F, {axisDCAx, axisDCAx}};
    HistogramConfigSpec h3DeltaPt{HistType::kTH3F, {axisPt, axisEta, axisDeltaPt}};
    HistogramConfigSpec hVtxZ{HistType::kTH1F, {axisVtxZ}};

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
  void runMuonSel(TCollision const& collision, TMFT const& tracksMFT, TMuons const& muons)
  {
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hVtxZ"), collision.posZ());

    for (const auto& muon : muons) {
      if (muon.trackType() == trkType) {
        const auto eta(muon.eta()), pDca(muon.pDca()), rAbs(muon.rAtAbsorberEnd());
        const auto dcaXY(std::sqrt(std::pow(muon.fwdDcaX(), 2.) + std::pow(muon.fwdDcaY(), 2.)));
        const auto pt(muon.pt());
        const auto charge(muon.sign());
        const auto chi2(muon.chi2MatchMCHMFT());

        auto trkMFT = muon.template matchMFTTrack_as<MFTTracksExtra>();
        if (reduceAmbMFT && trkMFT.has_ambMFTtrack()) {
          continue;
        }
        // histograms before the acceptance cuts
        registry.fill(HIST("hMuBeforeCuts"), pt, eta, dcaXY, pDca, charge, chi2);
        if (muon.matchMCHTrackId() > 0.) {
          auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
          registry.fill(HIST("h3DeltaPtBeforeCuts"), pt, eta, muonType3.pt() - pt);
        }
        // standard muon acceptance cuts
        if ((eta >= etaUp) || (eta < etaLow)) {
          continue;
        }
        if ((rAbs >= rAbsUp) || (rAbs < rAbsLow)) {
          continue;
        }
        if ((rAbs < rAbsMid) && (pDca >= pDcaLow)) {
          continue;
        }
        if ((rAbs >= rAbsMid) && (pDca >= pDcaUp)) {
          continue;
        }
        if ((muon.chi2() >= 1e6) || (muon.chi2() < 0.)) {
          continue;
        }
        // histograms after acceptance cuts
        registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, pDca, charge, chi2);
        registry.fill(HIST("h2DCA"), muon.fwdDcaX(), muon.fwdDcaY());
        if (muon.matchMCHTrackId() > 0.) {
          auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
          registry.fill(HIST("h3DeltaPtAfterCuts"), pt, eta, muonType3.pt() - pt);
        }
      }
    }
  }

  template <typename TCollision, typename TMFT, typename TMuons, typename TMC>
  void runMuonSelMc(TCollision const& collision, TMFT const& tracksMFT, TMuons const& muons, TMC const& mc)
  {
    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hVtxZ"), collision.posZ());

    for (const auto& muon : muons) {
      if (muon.trackType() == trkType) {
        const auto eta(muon.eta()), pDca(muon.pDca());
        const auto dcaXY(std::sqrt(std::pow(muon.fwdDcaX(), 2.) + std::pow(muon.fwdDcaY(), 2.)));
        const auto pt(muon.pt());
        const auto charge(muon.sign());
        const auto chi2(muon.chi2MatchMCHMFT());

        auto trkMFT = muon.template matchMFTTrack_as<MFTTracksExtra>();
        if (reduceAmbMFT && trkMFT.has_ambMFTtrack()) {
          continue;
        }
        // histograms before the acceptance cuts
        registry.fill(HIST("hMuBeforeCuts"), pt, eta, dcaXY, pDca, charge, chi2);
        if (muon.matchMCHTrackId() > 0.) {
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
        if ((eta >= etaUp) || (eta < etaLow)) {
          continue;
        }
        if ((muon.chi2() >= 1e6) || (muon.chi2() < 0.)) {
          continue;
        }
        // histograms after acceptance cuts
        registry.fill(HIST("hMuAfterCuts"), pt, eta, dcaXY, pDca, charge, chi2);
        registry.fill(HIST("h2DCA"), muon.fwdDcaX(), muon.fwdDcaY());
        if (muon.matchMCHTrackId() > 0.) {
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
  }

  void processMuon(soa::Filtered<MyCollisions>::iterator const& collision, MFTTracksExtra const& tracksMFT, MyMuons const& muons)
  {
    runMuonSel(collision, tracksMFT, muons);
  }
  void processMuonMc(soa::Filtered<MyCollisions>::iterator const& collision, MFTTracksExtra const& tracksMFT, soa::Filtered<MyMcMuons> const& muons, aod::McParticles const& mc)
  {
    runMuonSelMc(collision, tracksMFT, muons, mc);
  }

  PROCESS_SWITCH(HfTaskSingleMuonSelection, processMuon, "run muon selection with real data", true);
  PROCESS_SWITCH(HfTaskSingleMuonSelection, processMuonMc, "run muon selection with MC data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ProduceMFTAmbii>(cfgc),
    adaptAnalysisTask<HfTaskSingleMuonSelection>(cfgc),
  };
}
