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

/// \file taskSingleMuonMult.cxx
/// \brief Task used to study the Open heavy flavour decay muon production as a function of multiplicity.
/// \author Md Samsul Islam <md.samsul.islam@cern.ch>, IITB

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::fwdtrack;

struct HfTaskSingleMuonMult {

  // enum for event selection bins
  enum EventSelection {
    AllEvents = 0,
    Sel8,
    VtxZAfterSel,
    NEventSelection
  };

  // enum for muon track selection bins
  enum MuonSelection {
    NoCut = 0,
    EtaCut,
    RAbsorbCut,
    PDcaCut,
    Chi2Cut,
    NMuonSelection
  };

  Configurable<float> zVtxMax{"zVtxMax", 10., "maxium z of primary vertex [cm]"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.15, "minimum pt of tracks"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "maximum pseudorapidity of tracks"};
  Configurable<float> etaMin{"etaMin", -3.6, "minimum pseudorapidity"};
  Configurable<float> etaMax{"etaMax", -2.5, "maximum pseudorapidity"};
  Configurable<float> pDcaMin{"pDcaMin", 324., "p*DCA value for small RAbsorb"};
  Configurable<float> pDcaMax{"pDcaMax", 594., "p*DCA value for large RAbsorb"};
  Configurable<float> rAbsorbMin{"rAbsorbMin", 17.6, "R at absorber end minimum value"};
  Configurable<float> rAbsorbMax{"rAbsorbMax", 89.5, "R at absorber end maximum value"};
  Configurable<float> rAbsorbMid{"rAbsorbMid", 26.5, "R at absorber end split point for different p*DCA selections"};
  Configurable<float> chi2Max{"chi2Max", 1e6f, "MCH-MFT matching chi2 maximum value"};
  Configurable<bool> reduceOrphMft{"reduceOrphMft", true, "reduce orphan MFT tracks"};

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;
  using MyMuons = soa::Join<aod::FwdTracks, aod::FwdTracksDCA>;
  using MyMcMuons = soa::Join<aod::FwdTracks, aod::McFwdTrackLabels, aod::FwdTracksDCA>;
  using MyTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksIU, aod::TracksDCA, aod::TrackSelection>>;

  // Filter Global Track for Multiplicty
  Filter trackFilter = ((nabs(aod::track::eta) < etaTrackMax) && (aod::track::pt > ptTrackMin));

  // Number the types of muon tracks
  static constexpr uint8_t NTrackTypes{ForwardTrackTypeEnum::MCHStandaloneTrack + 1};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    AxisSpec const axisCent = {101, -0.5, 100.5, "centrality"};
    AxisSpec const axisEvent{NEventSelection, 0, NEventSelection, "Event Selection"};
    AxisSpec const axisVtxZ{80, -20., 20., "#it{z}_{vtx} (cm)"};
    AxisSpec const axisMuon{NMuonSelection, 0, NMuonSelection, "Muon Selection"};
    AxisSpec const axisNCh{500, 0.5, 500.5, "#it{N}_{ch}"};
    AxisSpec const axisNMu{20, -0.5, 19.5, "#it{N}_{#mu}"};
    AxisSpec const axisPt{1000, 0., 500., "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisEta{1000, -5., 5., "#it{#eta}"};
    AxisSpec const axisTheta{500, 170., 180., "#it{#theta}"};
    AxisSpec const axisRAbsorb{1000, 0., 100., "#it{R}_{Absorb} (cm)"};
    AxisSpec const axisDCA{500, 0., 5., "#it{DCA}_{xy} (cm)"};
    AxisSpec const axisChi2MatchMCHMFT{1000, 0., 1000., "MCH-MFT matching #chi^{2}"};
    AxisSpec const axisSign{5, -2.5, 2.5, "Charge"};
    AxisSpec const axisPDca{100000, 0, 100000, "#it{p} #times DCA (GeV/#it{c} * cm)"};
    AxisSpec const axisDCAx{1000, -5., 5., "#it{DCA}_{x or y} (cm)"};
    AxisSpec const axisEtaDif{200, -2., 2., "#it{#eta} diff"};
    AxisSpec const axisDeltaPt{10000, -50, 50, "#Delta #it{p}_{T} (GeV/#it{c})"};
    AxisSpec const axisTrackType{8, -1.5, 6.5, "TrackType"};
    AxisSpec const axisPtDif{200, -2., 2., "#it{p}_{T} diff (GeV/#it{c})"};

    registry.add("hCentrality", "Centrality Percentile", {HistType::kTH1F, {axisCent}});
    registry.add("hEventSel", " Number of Events", {HistType::kTH1F, {axisEvent}});
    registry.add("hNch", "Charged Particle Multiplicity", {HistType::kTH1F, {axisNCh}});
    registry.add("hVtxZBeforeSel", "Z-vertex distribution before zVtx Cut", {HistType::kTH1F, {axisVtxZ}});
    registry.add("hVtxZAfterSel", "Z-vertex distribution after zVtx Cut", {HistType::kTH1F, {axisVtxZ}});

    registry.add("hMuonSel", "Selection of muon tracks at various kinematic cuts", {HistType::kTH1F, {axisMuon}});
    registry.add("hMuBeforeMatchMFT", "Muon information before any Kinemeatic cuts applied", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("hMuBeforeAccCuts", "Muon information before applying Acceptance cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("h3DCABeforeAccCuts", "DCAx,DCAy,DCAz information before Acceptance cuts", {HistType::kTH3F, {axisDCAx, axisDCAx, axisTrackType}});
    registry.add("hMuDeltaPtBeforeAccCuts", "Muon information with DeltaPt before applying Acceptance cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDeltaPt, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("hMuAfterEtaCuts", "Muon information after applying Eta cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("hMuAfterRAbsorbCuts", "Muon information after applying RAbsorb cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("hMuAfterPdcaCuts", "Muon information after applying Pdca cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("hMuAfterAccCuts", "Muon information after applying all Kinematic cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDCA, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});
    registry.add("h3DCAAfterAccCuts", "DCAx,DCAy,DCAz information after Acceptance cuts", {HistType::kTH3F, {axisDCAx, axisDCAx, axisTrackType}});
    registry.add("hMuDeltaPtAfterAccCuts", "Muon information with DeltaPt after applying Acceptance cuts", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisTheta, axisRAbsorb, axisDeltaPt, axisPDca, axisChi2MatchMCHMFT, axisTrackType}, 10});

    registry.add("hTHnTrk", "Muon information with multiplicity", {HistType::kTHnSparseF, {axisCent, axisNCh, axisPt, axisEta, axisSign}, 5});
    registry.add("h3MultNchNmu", "Number of muons and multiplicity", {HistType::kTH3F, {axisCent, axisNCh, axisNMu}});
    registry.add("hMultNchNmuTrackType", "Number of muons with different types and multiplicity", {HistType::kTHnSparseF, {axisCent, axisNCh, axisNMu, axisTrackType}, 4});

    auto hEvstat = registry.get<TH1>(HIST("hEventSel"));
    auto* xEv = hEvstat->GetXaxis();
    xEv->SetBinLabel(AllEvents + 1, "All events");
    xEv->SetBinLabel(Sel8 + 1, "sel8");
    xEv->SetBinLabel(VtxZAfterSel + 1, "VtxZAfterSel");

    auto hMustat = registry.get<TH1>(HIST("hMuonSel"));
    auto* xMu = hMustat->GetXaxis();
    xMu->SetBinLabel(NoCut + 1, "noCut");
    xMu->SetBinLabel(EtaCut + 1, "etaCut");
    xMu->SetBinLabel(RAbsorbCut + 1, "RAbsorbCut");
    xMu->SetBinLabel(PDcaCut + 1, "pDcaCut");
    xMu->SetBinLabel(Chi2Cut + 1, "chi2Cut");
  }

  void process(MyCollisions::iterator const& collision,
               MyTracks const& tracks,
               MyMuons const& muons)
  {
    registry.fill(HIST("hEventSel"), AllEvents);

    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventSel"), Sel8);
    registry.fill(HIST("hVtxZBeforeSel"), collision.posZ());

    if (std::abs(collision.posZ()) > zVtxMax) {
      return;
    }
    registry.fill(HIST("hEventSel"), VtxZAfterSel);
    registry.fill(HIST("hVtxZAfterSel"), collision.posZ());

    // T0M centrality
    const auto cent = collision.centFT0M();
    registry.fill(HIST("hCentrality"), cent);

    // Charged particles
    for (const auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
    }

    auto nCh{tracks.size()};
    if (nCh < 1) {
      return;
    }
    registry.fill(HIST("hNch"), nCh);

    for (const auto& track : tracks) {
      registry.fill(HIST("hTHnTrk"), cent, nCh, track.pt(), track.eta(), track.sign());
    }

    // muons per event
    int nMu{0};
    int nMuType[NTrackTypes] = {0};

    for (const auto& muon : muons) {
      const auto pt{muon.pt()}, eta{muon.eta()}, theta{90.0f - ((std::atan(muon.tgl())) * constants::math::Rad2Deg)}, pDca{muon.pDca()}, rAbsorb{muon.rAtAbsorberEnd()}, chi2{muon.chi2MatchMCHMFT()};
      const auto dcaXY{RecoDecay::sqrtSumOfSquares(muon.fwdDcaX(), muon.fwdDcaY())};
      const auto muTrackType{muon.trackType()};
      float dptBefore{0.}, dptAfter{0.};

      registry.fill(HIST("hMuBeforeMatchMFT"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);

      // histograms before the acceptance cuts
      registry.fill(HIST("hMuonSel"), NoCut);
      registry.fill(HIST("hMuBeforeAccCuts"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);
      registry.fill(HIST("h3DCABeforeAccCuts"), muon.fwdDcaX(), muon.fwdDcaY(), muTrackType);

      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        dptBefore = muonType3.pt() - pt;
      }
      registry.fill(HIST("hMuDeltaPtBeforeAccCuts"), cent, nCh, pt, eta, theta, rAbsorb, dptBefore, pDca, chi2, muTrackType);

      // Apply various standard muon acceptance cuts
      // eta cuts
      if ((eta >= etaMax) || (eta < etaMin)) {
        continue;
      }
      registry.fill(HIST("hMuonSel"), EtaCut);
      registry.fill(HIST("hMuAfterEtaCuts"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);

      // Rabsorb cuts
      if ((rAbsorb < rAbsorbMin) || (rAbsorb >= rAbsorbMax)) {
        continue;
      }
      registry.fill(HIST("hMuonSel"), RAbsorbCut);
      registry.fill(HIST("hMuAfterRAbsorbCuts"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);

      if ((rAbsorb < rAbsorbMid) && (pDca >= pDcaMin)) {
        continue;
      }
      if ((rAbsorb >= rAbsorbMid) && (pDca >= pDcaMax)) {
        continue;
      }
      registry.fill(HIST("hMuonSel"), PDcaCut);
      registry.fill(HIST("hMuAfterPdcaCuts"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);

      //  MCH-MFT matching chi2
      if (muon.chi2() >= chi2Max) {
        continue;
      }
      registry.fill(HIST("hMuonSel"), Chi2Cut);

      // histograms after acceptance cuts
      registry.fill(HIST("hMuAfterAccCuts"), cent, nCh, pt, eta, theta, rAbsorb, dcaXY, pDca, chi2, muTrackType);
      registry.fill(HIST("h3DCAAfterAccCuts"), muon.fwdDcaX(), muon.fwdDcaY(), muTrackType);

      if (muon.has_matchMCHTrack()) {
        auto muonType3 = muon.template matchMCHTrack_as<MyMuons>();
        dptAfter = muonType3.pt() - pt;
      }
      registry.fill(HIST("hMuDeltaPtAfterAccCuts"), cent, nCh, pt, eta, theta, rAbsorb, dptAfter, pDca, chi2, muTrackType);

      nMu++;
      nMuType[muTrackType]++;
    }

    registry.fill(HIST("h3MultNchNmu"), cent, nCh, nMu);

    // Fill number of muons of various types with multiplicity
    for (auto indexType{0u}; indexType < NTrackTypes; ++indexType) {
      if (nMuType[indexType] > 0) {
        registry.fill(HIST("hMultNchNmuTrackType"), cent, nCh, nMuType[indexType], indexType);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskSingleMuonMult>(cfgc)};
}
