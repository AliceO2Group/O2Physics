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
/// \file hypertriton3bodyanalysis.cxx
/// \brief Standard analysis workflow for hypertriton 3-body decay
/// \author Yuanzhe Wang <yuanzhe.wang@cern.ch>

#include <cmath>
#include <array>
#include <cstdlib>
#include <algorithm>
#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
// #include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
// using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe, aod::pidTOFFullDe>; // For TOF PID check
using MCLabeledTracksIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct hypertriton3bodyQa {
  // Basic checks
  HistogramRegistry registry{
    "registry",
    {
      {"hVtxRadius", "hVtxRadius", {HistType::kTH1F, {{1000, 0.0f, 100.0f, "cm"}}}},
      {"hVtxCosPA", "hVtxCosPA", {HistType::kTH1F, {{1000, 0.9f, 1.0f}}}},
      {"hPtProton", "hPtProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtPionMinus", "hPtPionMinus", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtDeuteron", "hPtDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiProton", "hPtAntiProton", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtPionPlus", "hPtPionPlus", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hPtAntiDeuteron", "hPtAntiDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f}}}},
      {"hDCAXYProtonToPV", "hDCAXYProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAXYPionToPV", "hDCAXYPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAXYDeuteronToPV", "hDCAXYDeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAProtonToPV", "hDCAProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAPionToPV", "hDCAPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCADeuteronToPV", "hDCADeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hProtonTPCNcls", "hProtonTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hPionTPCNcls", "hPionTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hDeuteronTPCNcls", "hDeuteronTPCNcls", {HistType::kTH1F, {{300, 0, 300, "TPC cluster"}}}},
      {"hDCAVtxDau", "hDCAVtxDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hVtxPt", "hVtxPt", {HistType::kTH1F, {{200, 0.0f, 10.0f, "p_{T}"}}}},
      {"hTOFPIDDeuteron", "hTOFPIDDeuteron", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hDeuTOFNsigma", "Deuteron TOF Nsigma distribution", {HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {2000, -100, 100, "TOF n#sigma"}}}},
      {"hDeuTOFNsigmaWithTPC", "Deuteron TOF Nsigma distribution", {HistType::kTH2F, {{1200, -6, 6, "#it{p} (GeV/#it{c})"}, {1000, -100, 100, "TOF n#sigma"}}}},
    },
  };

  void init(InitContext const&)
  {
    AxisSpec massAxis = {120, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"};
    registry.add("hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {massAxis}});
    registry.add("hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {massAxis}});
    // Check for selection criteria
    registry.add("hDiffRVtxProton", "hDiffRVtxProton", HistType::kTH1F, {{100, -10, 10}});     // difference between the radius of decay vertex and minR of proton
    registry.add("hDiffRVtxPion", "hDiffRVtxPion", HistType::kTH1F, {{100, -10, 10}});         // difference between the radius of decay vertex and minR of pion
    registry.add("hDiffRVtxDeuteron", "hDiffRVtxDeuteron", HistType::kTH1F, {{100, -10, 10}}); // difference between the radius of decay vertex and minR of deuteron
    registry.add("hDiffDaughterR", "hDiffDaughterR", HistType::kTH1F, {{10000, -100, 100}});   // difference between minR of pion&proton and R of deuteron(bachelor)
  }

  void process(aod::Collision const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, FullTracksExtIU const& /*tracks*/)
  {
    for (const auto& vtx : vtx3bodydatas) {
      auto track0 = vtx.track0_as<FullTracksExtIU>();
      auto track1 = vtx.track1_as<FullTracksExtIU>();
      auto track2 = vtx.track2_as<FullTracksExtIU>();

      registry.fill(HIST("hVtxRadius"), vtx.vtxradius());
      registry.fill(HIST("hVtxCosPA"), vtx.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
      registry.fill(HIST("hDCAVtxDau"), vtx.dcaVtxdaughters());
      registry.fill(HIST("hVtxPt"), vtx.pt());
      registry.fill(HIST("hMassHypertriton"), vtx.mHypertriton());
      registry.fill(HIST("hMassAntiHypertriton"), vtx.mAntiHypertriton());
      if (std::abs(track2.tpcNSigmaDe()) < 5) {
        registry.fill(HIST("hDeuTOFNsigmaWithTPC"), track2.tpcInnerParam() * track2.sign(), vtx.tofNSigmaBachDe());
      }
      if (track2.sign() > 0) {
        registry.fill(HIST("hPtProton"), track0.pt());
        registry.fill(HIST("hPtPionMinus"), track1.pt());
        registry.fill(HIST("hPtDeuteron"), track2.pt());
        registry.fill(HIST("hDCAXYProtonToPV"), vtx.dcaXYtrack0topv());
        registry.fill(HIST("hDCAXYPionToPV"), vtx.dcaXYtrack1topv());
        registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack0topv());
        registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack1topv());
        registry.fill(HIST("hProtonTPCNcls"), track0.tpcNClsCrossedRows());
        registry.fill(HIST("hPionTPCNcls"), track1.tpcNClsCrossedRows());
        registry.fill(HIST("hDiffRVtxProton"), track0.x() - vtx.vtxradius());
        registry.fill(HIST("hDiffRVtxPion"), track1.x() - vtx.vtxradius());
      } else {
        registry.fill(HIST("hPtPionPlus"), track0.pt());
        registry.fill(HIST("hPtAntiProton"), track1.pt());
        registry.fill(HIST("hPtAntiDeuteron"), track2.pt());
        registry.fill(HIST("hDCAXYProtonToPV"), vtx.dcaXYtrack1topv());
        registry.fill(HIST("hDCAXYPionToPV"), vtx.dcaXYtrack0topv());
        registry.fill(HIST("hDCAProtonToPV"), vtx.dcatrack1topv());
        registry.fill(HIST("hDCAPionToPV"), vtx.dcatrack0topv());
        registry.fill(HIST("hProtonTPCNcls"), track1.tpcNClsCrossedRows());
        registry.fill(HIST("hPionTPCNcls"), track0.tpcNClsCrossedRows());
        registry.fill(HIST("hDiffRVtxProton"), track1.x() - vtx.vtxradius());
        registry.fill(HIST("hDiffRVtxPion"), track0.x() - vtx.vtxradius());
      }
      registry.fill(HIST("hDCAXYDeuteronToPV"), vtx.dcaXYtrack2topv());
      registry.fill(HIST("hDCADeuteronToPV"), vtx.dcatrack2topv());
      registry.fill(HIST("hDeuteronTPCNcls"), track2.tpcNClsCrossedRows());
      registry.fill(HIST("hTOFPIDDeuteron"), vtx.tofNSigmaBachDe());
      registry.fill(HIST("hDeuTOFNsigma"), track2.tpcInnerParam() * track2.sign(), vtx.tofNSigmaBachDe());
      registry.fill(HIST("hDiffRVtxDeuteron"), track2.x() - vtx.vtxradius());
      float diffTrackR = track2.x() - std::min(track0.x(), track1.x());
      registry.fill(HIST("hDiffDaughterR"), diffTrackR);
    }
  }
};

struct hypertriton3bodyAnalysis {

  Preslice<aod::Vtx3BodyDatas> perCollisionVtx3BodyDatas = o2::aod::vtx3body::collisionId;

  // Selection criteria
  Configurable<double> vtxcospa{"vtxcospa", 0.99, "Vtx CosPA"};         // double -> N.B. dcos(x)/dx = 0 at x=0)
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"}; // loose cut
  Configurable<float> dcapiontopv{"dcapiontopv", .05, "DCA Pion To PV"};
  Configurable<float> etacut{"etacut", 0.9, "etacut"};
  Configurable<float> rapiditycut{"rapiditycut", 1, "rapiditycut"};
  Configurable<float> tofPIDNSigmaMin{"tofPIDNSigmaMin", -5, "tofPIDNSigmaMin"};
  Configurable<float> tofPIDNSigmaMax{"tofPIDNSigmaMax", 5, "tofPIDNSigmaMax"};
  Configurable<float> tpcPIDNSigmaCut{"tpcPIDNSigmaCut", 5, "tpcPIDNSigmaCut"};
  Configurable<bool> event_sel8_selection{"event_sel8_selection", true, "event selection count post sel8 cut"};
  Configurable<bool> mc_event_selection{"mc_event_selection", true, "mc event selection count post kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", true, "event selection count post poZ cut"};
  Configurable<float> lifetimecut{"lifetimecut", 40., "lifetimecut"}; // ct
  Configurable<float> minProtonPt{"minProtonPt", 0.3, "minProtonPt"};
  Configurable<float> maxProtonPt{"maxProtonPt", 5, "maxProtonPt"};
  Configurable<float> minPionPt{"minPionPt", 0.1, "minPionPt"};
  Configurable<float> maxPionPt{"maxPionPt", 1.2, "maxPionPt"};
  Configurable<float> minDeuteronPt{"minDeuteronPt", 0.6, "minDeuteronPt"};
  Configurable<float> maxDeuteronPt{"maxDeuteronPt", 10, "maxDeuteronPt"};
  Configurable<float> minDeuteronPUseTOF{"minDeuteronPUseTOF", 1, "minDeuteronPt Enable TOF PID"};
  Configurable<float> h3LMassLowerlimit{"h3LMassLowerlimit", 2.96, "Hypertriton mass lower limit"};
  Configurable<float> h3LMassUpperlimit{"h3LMassUpperlimit", 3.04, "Hypertriton mass upper limit"};
  Configurable<int> mintpcNClsproton{"mintpcNClsproton", 90, "min tpc Nclusters for proton"};
  Configurable<int> mintpcNClspion{"mintpcNClspion", 70, "min tpc Nclusters for pion"};
  Configurable<int> mintpcNClsdeuteron{"mintpcNClsdeuteron", 100, "min tpc Nclusters for deuteron"};

  Configurable<float> mcsigma{"mcsigma", 0.0015, "sigma of mc invariant mass fit"}; // obtained from MC
  Configurable<int> bachelorPdgCode{"bachelorPdgCode", 1000010020, "pdgCode of bachelor daughter"};
  Configurable<int> motherPdgCode{"motherPdgCode", 1010010030, "pdgCode of mother track"};

  // 3sigma region for Dalitz plot
  float lowersignallimit = o2::constants::physics::MassHyperTriton - 3 * mcsigma;
  float uppersignallimit = o2::constants::physics::MassHyperTriton + 3 * mcsigma;

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
      {"hCandidatesCounter", "hCandidatesCounter", {HistType::kTH1F, {{12, 0.0f, 12.0f}}}},
      {"hMassHypertriton", "hMassHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassAntiHypertriton", "hMassAntiHypertriton", {HistType::kTH1F, {{80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hMassHypertritonTotal", "hMassHypertritonTotal", {HistType::kTH1F, {{300, 2.9f, 3.2f, "Inv. Mass (GeV/c^{2})"}}}},
      {"hPtProton", "hPtProton", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hPtPionMinus", "hPtPionMinus", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hPtDeuteron", "hPtDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hPtAntiProton", "hPtAntiProton", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hPtPionPlus", "hPtPionPlus", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hPtAntiDeuteron", "hPtAntiDeuteron", {HistType::kTH1F, {{200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}}}},
      {"hDCAXYProtonToPV", "hDCAXYProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAXYPionToPV", "hDCAXYPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAXYDeuteronToPV", "hDCAXYDeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAProtonToPV", "hDCAProtonToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCAPionToPV", "hDCAPionToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hDCADeuteronToPV", "hDCADeuteronToPV", {HistType::kTH1F, {{1000, -10.0f, 10.0f, "cm"}}}},
      {"hProtonTPCNcls", "hProtonTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hPionTPCNcls", "hPionTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hDeuteronTPCNcls", "hDeuteronTPCNcls", {HistType::kTH1F, {{180, 0, 180, "TPC cluster"}}}},
      {"hVtxCosPA", "hVtxCosPA", {HistType::kTH1F, {{1000, 0.9f, 1.0f}}}},
      {"hDCAVtxDau", "hDCAVtxDau", {HistType::kTH1F, {{1000, 0.0f, 10.0f, "cm^{2}"}}}},
      {"hTOFPIDDeuteron", "hTOFPIDDeuteron", {HistType::kTH1F, {{2000, -100.0f, 100.0f}}}},
      {"hTPCPIDProton", "hTPCPIDProton", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hTPCPIDPion", "hTPCPIDPion", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hTPCPIDDeuteron", "hTPCPIDDeuteron", {HistType::kTH1F, {{120, -6.0f, 6.0f}}}},
      {"hProtonTPCBB", "hProtonTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hPionTPCBB", "hPionTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hDeuteronTPCBB", "hDeuteronTPCBB", {HistType::kTH2F, {{160, -8.0f, 8.0f, "p/z(GeV/c)"}, {200, 0.0f, 1000.0f, "TPCSignal"}}}},
      {"hProtonTPCVsPt", "hProtonTPCVsPt", {HistType::kTH2F, {{50, 0.0f, 5.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hPionTPCVsPt", "hPionTPCVsPt", {HistType::kTH2F, {{20, 0.0f, 2.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDeuteronTPCVsPt", "hDeuteronTPCVsPt", {HistType::kTH2F, {{80, 0.0f, 8.0f, "#it{p}_{T} (GeV/c)"}, {120, -6.0f, 6.0f, "TPC n#sigma"}}}},
      {"hDeuteronTOFVsPBeforeTOFCut", "hDeuteronTOFVsPBeforeTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronTOFVsPAfterTOFCut", "hDeuteronTOFVsPAfterTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},

      {"hDalitz", "hDalitz", {HistType::kTH2F, {{120, 7.85, 8.45, "M^{2}(dp) (GeV^{2}/c^{4})"}, {60, 1.1, 1.4, "M^{2}(p#pi) (GeV^{2}/c^{4})"}}}},
      {"h3dMassHypertriton", "h3dMassHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dMassAntiHypertriton", "h3dMassAntiHypertriton", {HistType::kTH3F, {{20, 0.0f, 100.0f, "Cent (%)"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      {"h3dTotalHypertriton", "h3dTotalHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},

      {"hDeuteronTOFVsPBeforeTOFCutSig", "hDeuteronTOFVsPBeforeTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronTOFVsPAfterTOFCutSig", "hDeuteronTOFVsPAfterTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"h3dTotalTrueHypertriton", "h3dTotalTrueHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {80, 2.96f, 3.04f, "Inv. Mass (GeV/c^{2})"}}}},
      // For TOF PID check
      /*{"hDeuteronDefaultTOFVsPBeforeTOFCut", "hDeuteronDefaultTOFVsPBeforeTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronDefaultTOFVsPAtferTOFCut", "hDeuteronDefaultTOFVsPAtferTOFCut", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronDefaultTOFVsPBeforeTOFCutSig", "hDeuteronDefaultTOFVsPBeforeTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},
      {"hDeuteronDefaultTOFVsPAfterTOFCutSig", "hDeuteronDefaultTOFVsPAfterTOFCutSig", {HistType::kTH2F, {{40, -10.0f, 10.0f, "p/z (GeV/c)"}, {40, -10.0f, 10.0f, "TOF n#sigma"}}}},*/
    },
  };

  //------------------------------------------------------------------
  // Fill stats histograms
  enum vtxstep { kCandAll = 0,
                 kCandDauEta,
                 kCandDauPt,
                 kCandTPCNcls,
                 kCandTPCPID,
                 kCandTOFPID,
                 kCandDcaToPV,
                 kCandRapidity,
                 kCandct,
                 kCandCosPA,
                 kCandDcaDau,
                 kCandInvMass,
                 kNCandSteps };

  struct {
    std::array<int32_t, kNCandSteps> candstats;
    std::array<int32_t, kNCandSteps> truecandstats;
  } statisticsRegistry;

  void resetHistos()
  {
    for (int ii = 0; ii < kNCandSteps; ii++) {
      statisticsRegistry.candstats[ii] = 0;
      statisticsRegistry.truecandstats[ii] = 0;
    }
  }
  void FillCandCounter(int kn, bool istrue = false)
  {
    statisticsRegistry.candstats[kn]++;
    if (istrue) {
      statisticsRegistry.truecandstats[kn]++;
    }
  }
  void fillHistos()
  {
    for (int ii = 0; ii < kNCandSteps; ii++) {
      registry.fill(HIST("hCandidatesCounter"), ii, statisticsRegistry.candstats[ii]);
      if (doprocessMC == true) {
        registry.fill(HIST("hTrueHypertritonCounter"), ii, statisticsRegistry.truecandstats[ii]);
      }
    }
  }

  ConfigurableAxis dcaBinning{"dca-binning", {200, 0.0f, 1.0f}, ""};
  ConfigurableAxis ptBinning{"pt-binning", {200, 0.0f, 10.0f}, ""};

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(1, "total");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(2, "sel8");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(3, "vertexZ");
    registry.get<TH1>(HIST("hEventCounter"))->GetXaxis()->SetBinLabel(4, "has Candidate");

    if (doprocessMC == true) {
      registry.add("hTrueHypertritonCounter", "hTrueHypertritonCounter", HistType::kTH1F, {{12, 0.0f, 12.0f}});
      auto hGeneratedHypertritonCounter = registry.add<TH1>("hGeneratedHypertritonCounter", "hGeneratedHypertritonCounter", HistType::kTH1F, {{2, 0.0f, 2.0f}});
      hGeneratedHypertritonCounter->GetXaxis()->SetBinLabel(1, "Total");
      hGeneratedHypertritonCounter->GetXaxis()->SetBinLabel(2, "3-body decay");
      registry.add("hPtGeneratedHypertriton", "hPtGeneratedHypertriton", HistType::kTH1F, {{200, 0.0f, 10.0f}});
      registry.add("hctGeneratedHypertriton", "hctGeneratedHypertriton", HistType::kTH1F, {{50, 0, 50, "ct(cm)"}});
      registry.add("hEtaGeneratedHypertriton", "hEtaGeneratedHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hRapidityGeneratedHypertriton", "hRapidityGeneratedHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hPtGeneratedAntiHypertriton", "hPtGeneratedAntiHypertriton", HistType::kTH1F, {{200, 0.0f, 10.0f}});
      registry.add("hctGeneratedAntiHypertriton", "hctGeneratedAntiHypertriton", HistType::kTH1F, {{50, 0, 50, "ct(cm)"}});
      registry.add("hEtaGeneratedAntiHypertriton", "hEtaGeneratedAntiHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
      registry.add("hRapidityGeneratedAntiHypertriton", "hRapidityGeneratedAntiHypertriton", HistType::kTH1F, {{40, -2.0f, 2.0f}});
    }

    TString CandCounterbinLabel[kNCandSteps] = {"Total", "TrackEta", "DauPt", "TPCNcls", "TPCPID", "d TOFPID", "PionDcatoPV", "MomRapidity", "Lifetime", "VtxCosPA", "VtxDcaDau", "InvMass"};
    for (int i{0}; i < kNCandSteps; i++) {
      registry.get<TH1>(HIST("hCandidatesCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
      if (doprocessMC == true) {
        registry.get<TH1>(HIST("hTrueHypertritonCounter"))->GetXaxis()->SetBinLabel(i + 1, CandCounterbinLabel[i]);
      }
    }
  }

  //------------------------------------------------------------------
  // Selections for candidates
  template <typename TCollisionTable, typename TTrackTable, typename TCandTable>
  bool SelectCand(TCollisionTable const& collision, TCandTable const& candData, TTrackTable const& trackProton, TTrackTable const& trackPion, TTrackTable const& trackDeuteron, bool isMatter, bool isTrueCand = false, double MClifetime = -1, double lPt = -1)
  {
    FillCandCounter(kCandAll, isTrueCand);

    // Selection on daughters
    if (std::abs(trackProton.eta()) > etacut || std::abs(trackPion.eta()) > etacut || std::abs(trackDeuteron.eta()) > etacut) {
      return false;
    }
    FillCandCounter(kCandDauEta, isTrueCand);

    if (trackProton.pt() < minProtonPt || trackProton.pt() > maxProtonPt || trackPion.pt() < minPionPt || trackPion.pt() > maxPionPt || trackDeuteron.pt() < minDeuteronPt || trackDeuteron.pt() > maxDeuteronPt) {
      return false;
    }
    FillCandCounter(kCandDauPt, isTrueCand);

    if (trackProton.tpcNClsFound() < mintpcNClsproton || trackPion.tpcNClsFound() < mintpcNClspion || trackDeuteron.tpcNClsFound() < mintpcNClsdeuteron) {
      return false;
    }
    FillCandCounter(kCandTPCNcls, isTrueCand);

    if (std::abs(trackProton.tpcNSigmaPr()) > tpcPIDNSigmaCut || std::abs(trackPion.tpcNSigmaPi()) > tpcPIDNSigmaCut || std::abs(trackDeuteron.tpcNSigmaDe()) > tpcPIDNSigmaCut) {
      return false;
    }
    FillCandCounter(kCandTPCPID, isTrueCand);

    // registry.fill(HIST("hDeuteronDefaultTOFVsPBeforeTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tofNSigmaDe());
    registry.fill(HIST("hDeuteronTOFVsPBeforeTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());
    if (isTrueCand) {
      // registry.fill(HIST("hDeuteronDefaultTOFVsPBeforeTOFCutSig"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tofNSigmaDe());
      registry.fill(HIST("hDeuteronTOFVsPBeforeTOFCutSig"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());
    }
    if ((candData.tofNSigmaBachDe() < tofPIDNSigmaMin || candData.tofNSigmaBachDe() > tofPIDNSigmaMax) && trackDeuteron.p() > minDeuteronPUseTOF) {
      return false;
    }
    FillCandCounter(kCandTOFPID, isTrueCand);
    // registry.fill(HIST("hDeuteronDefaultTOFVsPAtferTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tofNSigmaDe());
    registry.fill(HIST("hDeuteronTOFVsPAfterTOFCut"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());
    if (isTrueCand) {
      // registry.fill(HIST("hDeuteronDefaultTOFVsPAfterTOFCutSig"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tofNSigmaDe());
      registry.fill(HIST("hDeuteronTOFVsPAfterTOFCutSig"), trackDeuteron.sign() * trackDeuteron.p(), candData.tofNSigmaBachDe());
    }

    double dcapion = isMatter ? candData.dcatrack1topv() : candData.dcatrack0topv();
    if (std::abs(dcapion) < dcapiontopv) {
      return false;
    }
    FillCandCounter(kCandDcaToPV, isTrueCand);

    // Selection on candidate hypertriton
    if (std::abs(candData.yHypertriton()) > rapiditycut) {
      return false;
    }
    FillCandCounter(kCandRapidity, isTrueCand);

    double ct = candData.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassHyperTriton;
    if (ct > lifetimecut) {
      return false;
    }
    FillCandCounter(kCandct, isTrueCand);

    double cospa = candData.vtxcosPA(collision.posX(), collision.posY(), collision.posZ());
    if (cospa < vtxcospa) {
      return false;
    }
    FillCandCounter(kCandCosPA, isTrueCand);

    if (candData.dcaVtxdaughters() > dcavtxdau) {
      return false;
    }
    FillCandCounter(kCandDcaDau, isTrueCand);

    if ((isMatter && candData.mHypertriton() > h3LMassLowerlimit && candData.mHypertriton() < h3LMassUpperlimit)) {
      // Hypertriton
      registry.fill(HIST("hPtProton"), trackProton.pt());
      registry.fill(HIST("hPtPionMinus"), trackPion.pt());
      registry.fill(HIST("hPtDeuteron"), trackDeuteron.pt());
      registry.fill(HIST("hDCAXYProtonToPV"), candData.dcaXYtrack0topv());
      registry.fill(HIST("hDCAXYPionToPV"), candData.dcaXYtrack1topv());
      registry.fill(HIST("hDCAProtonToPV"), candData.dcatrack0topv());
      registry.fill(HIST("hDCAPionToPV"), candData.dcatrack1topv());

      registry.fill(HIST("hMassHypertriton"), candData.mHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mHypertriton());
      registry.fill(HIST("h3dMassHypertriton"), 0., candData.pt(), candData.mHypertriton()); // collision.centV0M() instead of 0. once available
      registry.fill(HIST("h3dTotalHypertriton"), ct, candData.pt(), candData.mHypertriton());
      if (candData.mHypertriton() > lowersignallimit && candData.mHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(std::array{std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, std::array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(std::array{std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}, std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
      if (isTrueCand) {
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, candData.mHypertriton());
      }
    } else if ((!isMatter && candData.mAntiHypertriton() > h3LMassLowerlimit && candData.mAntiHypertriton() < h3LMassUpperlimit)) {
      // AntiHypertriton
      registry.fill(HIST("hPtAntiProton"), trackProton.pt());
      registry.fill(HIST("hPtPionPlus"), trackPion.pt());
      registry.fill(HIST("hPtAntiDeuteron"), trackDeuteron.pt());
      registry.fill(HIST("hDCAXYProtonToPV"), candData.dcaXYtrack1topv());
      registry.fill(HIST("hDCAXYPionToPV"), candData.dcaXYtrack0topv());
      registry.fill(HIST("hDCAProtonToPV"), candData.dcatrack1topv());
      registry.fill(HIST("hDCAPionToPV"), candData.dcatrack0topv());

      registry.fill(HIST("hMassAntiHypertriton"), candData.mAntiHypertriton());
      registry.fill(HIST("hMassHypertritonTotal"), candData.mAntiHypertriton());
      registry.fill(HIST("h3dMassAntiHypertriton"), 0., candData.pt(), candData.mAntiHypertriton()); // collision.centV0M() instead of 0. once available
      registry.fill(HIST("h3dTotalHypertriton"), ct, candData.pt(), candData.mAntiHypertriton());
      if (candData.mAntiHypertriton() > lowersignallimit && candData.mAntiHypertriton() < uppersignallimit) {
        registry.fill(HIST("hDalitz"), RecoDecay::m2(std::array{std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, std::array{candData.pxtrack2(), candData.pytrack2(), candData.pztrack2()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron}), RecoDecay::m2(std::array{std::array{candData.pxtrack1(), candData.pytrack1(), candData.pztrack1()}, std::array{candData.pxtrack0(), candData.pytrack0(), candData.pztrack0()}}, std::array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged}));
      }
      if (isTrueCand) {
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, candData.mHypertriton());
      }
    } else {
      return false;
    }

    FillCandCounter(kCandInvMass, isTrueCand);

    registry.fill(HIST("hDCAXYDeuteronToPV"), candData.dcaXYtrack2topv());
    registry.fill(HIST("hDCADeuteronToPV"), candData.dcatrack2topv());
    registry.fill(HIST("hVtxCosPA"), candData.vtxcosPA(collision.posX(), collision.posY(), collision.posZ()));
    registry.fill(HIST("hDCAVtxDau"), candData.dcaVtxdaughters());
    registry.fill(HIST("hProtonTPCNcls"), trackProton.tpcNClsCrossedRows());
    registry.fill(HIST("hPionTPCNcls"), trackPion.tpcNClsCrossedRows());
    registry.fill(HIST("hDeuteronTPCNcls"), trackDeuteron.tpcNClsCrossedRows());
    registry.fill(HIST("hTPCPIDProton"), trackProton.tpcNSigmaPr());
    registry.fill(HIST("hTPCPIDPion"), trackPion.tpcNSigmaPi());
    registry.fill(HIST("hTPCPIDDeuteron"), trackDeuteron.tpcNSigmaDe());
    registry.fill(HIST("hProtonTPCBB"), trackProton.sign() * trackProton.p(), trackProton.tpcSignal());
    registry.fill(HIST("hPionTPCBB"), trackPion.sign() * trackPion.p(), trackPion.tpcSignal());
    registry.fill(HIST("hDeuteronTPCBB"), trackDeuteron.sign() * trackDeuteron.p(), trackDeuteron.tpcSignal());
    registry.fill(HIST("hProtonTPCVsPt"), trackProton.pt(), trackProton.tpcNSigmaPr());
    registry.fill(HIST("hPionTPCVsPt"), trackProton.pt(), trackPion.tpcNSigmaPi());
    registry.fill(HIST("hDeuteronTPCVsPt"), trackDeuteron.pt(), trackDeuteron.tpcNSigmaDe());
    registry.fill(HIST("hTOFPIDDeuteron"), candData.tofNSigmaBachDe());

    return true;
  }

  //------------------------------------------------------------------
  // Analysis process for a single candidate
  template <class TTrackClass, typename TCollisionTable, typename TCandTable>
  void CandidateAnalysis(TCollisionTable const& collision, TCandTable const& candData, bool& if_hasvtx, bool isTrueCand = false, double MClifetime = -1, double lPt = -1)
  {

    auto track0 = candData.template track0_as<TTrackClass>();
    auto track1 = candData.template track1_as<TTrackClass>();
    auto track2 = candData.template track2_as<TTrackClass>();

    bool isMatter = track2.sign() > 0;

    auto& trackProton = isMatter ? track0 : track1;
    auto& trackPion = isMatter ? track1 : track0;
    auto& trackDeuteron = track2;

    if (SelectCand(collision, candData, trackProton, trackPion, trackDeuteron, isMatter, isTrueCand, MClifetime, lPt)) {
      if_hasvtx = true;
    }
  }

  //------------------------------------------------------------------
  // collect information for generated hypertriton (should be called after event selection)
  void GetGeneratedH3LInfo(aod::McParticles const& particlesMC)
  {
    for (const auto& mcparticle : particlesMC) {
      if (std::abs(mcparticle.pdgCode()) != motherPdgCode) {
        continue;
      }
      registry.fill(HIST("hGeneratedHypertritonCounter"), 0.5);

      bool haveProton = false, havePionPlus = false, haveDeuteron = false;
      bool haveAntiProton = false, havePionMinus = false, haveAntiDeuteron = false;
      double MClifetime = -1;
      for (const auto& mcparticleDaughter : mcparticle.template daughters_as<aod::McParticles>()) {
        if (mcparticleDaughter.pdgCode() == 2212)
          haveProton = true;
        if (mcparticleDaughter.pdgCode() == -2212)
          haveAntiProton = true;
        if (mcparticleDaughter.pdgCode() == 211)
          havePionPlus = true;
        if (mcparticleDaughter.pdgCode() == -211)
          havePionMinus = true;
        if (mcparticleDaughter.pdgCode() == bachelorPdgCode) {
          haveDeuteron = true;
          MClifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
        if (mcparticleDaughter.pdgCode() == -bachelorPdgCode) {
          haveAntiDeuteron = true;
          MClifetime = RecoDecay::sqrtSumOfSquares(mcparticleDaughter.vx() - mcparticle.vx(), mcparticleDaughter.vy() - mcparticle.vy(), mcparticleDaughter.vz() - mcparticle.vz()) * o2::constants::physics::MassHyperTriton / mcparticle.p();
        }
      }
      if (haveProton && havePionMinus && haveDeuteron && mcparticle.pdgCode() == motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hPtGeneratedHypertriton"), mcparticle.pt());
        registry.fill(HIST("hctGeneratedHypertriton"), MClifetime);
        registry.fill(HIST("hEtaGeneratedHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedHypertriton"), mcparticle.y());
      } else if (haveAntiProton && havePionPlus && haveAntiDeuteron && mcparticle.pdgCode() == -motherPdgCode) {
        registry.fill(HIST("hGeneratedHypertritonCounter"), 1.5);
        registry.fill(HIST("hPtGeneratedAntiHypertriton"), mcparticle.pt());
        registry.fill(HIST("hctGeneratedAntiHypertriton"), MClifetime);
        registry.fill(HIST("hEtaGeneratedAntiHypertriton"), mcparticle.eta());
        registry.fill(HIST("hRapidityGeneratedAntiHypertriton"), mcparticle.y());
      }
    }
  }

  //------------------------------------------------------------------
  // process real data analysis
  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, aod::Vtx3BodyDatas const& vtx3bodydatas, FullTracksExtIU const& /*tracks*/)
  {
    registry.fill(HIST("hEventCounter"), 0.5);
    if (event_sel8_selection && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 1.5);
    if (event_posZ_selection && std::abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }
    registry.fill(HIST("hEventCounter"), 2.5);

    bool if_hasvtx = false;

    for (const auto& vtx : vtx3bodydatas) {
      CandidateAnalysis<FullTracksExtIU>(collision, vtx, if_hasvtx);
    }

    if (if_hasvtx)
      registry.fill(HIST("hEventCounter"), 3.5);
    fillHistos();
    resetHistos();
  }
  PROCESS_SWITCH(hypertriton3bodyAnalysis, processData, "Real data analysis", true);

  //------------------------------------------------------------------
  // process mc analysis
  void processMC(soa::Join<aod::Collisions, aod::EvSels> const& collisions, aod::Vtx3BodyDatas const& vtx3bodydatas, aod::McParticles const& particlesMC, MCLabeledTracksIU const& /*tracks*/)
  {
    GetGeneratedH3LInfo(particlesMC);

    for (const auto& collision : collisions) {
      registry.fill(HIST("hEventCounter"), 0.5);
      if (mc_event_selection && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
        continue;
      }
      registry.fill(HIST("hEventCounter"), 1.5);
      if (event_posZ_selection && std::abs(collision.posZ()) > 10.f) { // 10cm
        continue;
      }
      registry.fill(HIST("hEventCounter"), 2.5);

      bool if_hasvtx = false;
      auto vtxsthiscol = vtx3bodydatas.sliceBy(perCollisionVtx3BodyDatas, collision.globalIndex());

      for (const auto& vtx : vtxsthiscol) {
        // int lLabel = -1;
        int lPDG = -1;
        float lPt = -1;
        double MClifetime = -1;
        bool isTrueCand = false;
        auto track0 = vtx.track0_as<MCLabeledTracksIU>();
        auto track1 = vtx.track1_as<MCLabeledTracksIU>();
        auto track2 = vtx.track2_as<MCLabeledTracksIU>();
        if (track0.has_mcParticle() && track1.has_mcParticle() && track2.has_mcParticle()) {
          auto lMCTrack0 = track0.mcParticle_as<aod::McParticles>();
          auto lMCTrack1 = track1.mcParticle_as<aod::McParticles>();
          auto lMCTrack2 = track2.mcParticle_as<aod::McParticles>();
          if (lMCTrack0.has_mothers() && lMCTrack1.has_mothers() && lMCTrack2.has_mothers()) {
            for (const auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
              for (const auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
                for (const auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
                  if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
                    // lLabel = lMother1.globalIndex();
                    lPt = lMother1.pt();
                    lPDG = lMother1.pdgCode();
                    if ((lPDG == motherPdgCode && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == bachelorPdgCode) ||
                        (lPDG == -motherPdgCode && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -bachelorPdgCode)) {
                      isTrueCand = true;
                      MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
                    }
                  }
                }
              }
            }
          }
        }

        CandidateAnalysis<MCLabeledTracksIU>(collision, vtx, if_hasvtx, isTrueCand, MClifetime, lPt);
      }

      if (if_hasvtx)
        registry.fill(HIST("hEventCounter"), 3.5);
      fillHistos();
      resetHistos();
    }
  }
  PROCESS_SWITCH(hypertriton3bodyAnalysis, processMC, "MC analysis", false);
};

// check vtx3body with mclabels
struct hypertriton3bodyLabelCheck {

  Configurable<bool> mc_event_selection{"mc_event_selection", true, "mc event selection count post kIsTriggerTVX and kNoTimeFrameBorder"};
  Configurable<bool> event_posZ_selection{"event_posZ_selection", false, "event selection count post poZ cut"};
  Configurable<float> tpcPIDNSigmaCut{"tpcPIDNSigmaCut", 5, "tpcPIDNSigmaCut"};
  Configurable<int> motherPdgCode{"motherPdgCode", 1010010030, "pdgCode of mother track"};

  HistogramRegistry registry{"registry", {}};

  void init(InitContext const&)
  {
    if (doprocessData == false) {
      auto hLabeledVtxCounter = registry.add<TH1>("hLabeledVtxCounter", "hLabeledVtxCounter", HistType::kTH1F, {{3, 0.0f, 3.0f}});
      hLabeledVtxCounter->GetXaxis()->SetBinLabel(1, "Readin");
      hLabeledVtxCounter->GetXaxis()->SetBinLabel(2, "TrueMCH3L");
      hLabeledVtxCounter->GetXaxis()->SetBinLabel(3, "Nonrepetitive");
      registry.add("hMassTrueH3L", "hMassTrueH3L", HistType::kTH1F, {{80, 2.96f, 3.04f}});
      registry.add("hMassTrueH3LMatter", "hMassTrueH3LMatter", HistType::kTH1F, {{80, 2.96f, 3.04f}});
      registry.add("hMassTrueH3LAntiMatter", "hMassTrueH3LAntiMatter", HistType::kTH1F, {{80, 2.96f, 3.04f}});
      auto hPIDCounter = registry.add<TH1>("hPIDCounter", "hPIDCounter", HistType::kTH1F, {{6, 0.0f, 6.0f}});
      hPIDCounter->GetXaxis()->SetBinLabel(1, "H3L Proton PID > 5");
      hPIDCounter->GetXaxis()->SetBinLabel(2, "H3L Pion PID > 5");
      hPIDCounter->GetXaxis()->SetBinLabel(3, "H3L Deuteron PID > 5");
      hPIDCounter->GetXaxis()->SetBinLabel(4, "#bar{H3L} Proton PID > 5");
      hPIDCounter->GetXaxis()->SetBinLabel(5, "#bar{H3L} Pion PID > 5");
      hPIDCounter->GetXaxis()->SetBinLabel(6, "#bar{H3L} Deuteron PID > 5");
      auto hHypertritonCounter = registry.add<TH1>("hHypertritonCounter", "hHypertritonCounter", HistType::kTH1F, {{4, 0.0f, 4.0f}});
      hHypertritonCounter->GetXaxis()->SetBinLabel(1, "H3L");
      hHypertritonCounter->GetXaxis()->SetBinLabel(2, "H3L daughters pass PID");
      hHypertritonCounter->GetXaxis()->SetBinLabel(3, "#bar{H3L}");
      hHypertritonCounter->GetXaxis()->SetBinLabel(4, "#bar{H3L} daughters pass PID");
      auto hDecay3BodyCounter = registry.add<TH1>("hDecay3BodyCounter", "hDecay3BodyCounter", HistType::kTH1F, {{5, 0.0f, 5.0f}});
      hDecay3BodyCounter->GetXaxis()->SetBinLabel(1, "Total");
      hDecay3BodyCounter->GetXaxis()->SetBinLabel(2, "True H3L");
      hDecay3BodyCounter->GetXaxis()->SetBinLabel(3, "Unduplicated H3L");
      hDecay3BodyCounter->GetXaxis()->SetBinLabel(4, "Correct collision");
      hDecay3BodyCounter->GetXaxis()->SetBinLabel(5, "Same ColID for daughters");
      registry.add("hDiffRVtxProton", "hDiffRVtxProton", HistType::kTH1F, {{100, -10, 10}});     // difference between the radius of decay vertex and minR of proton
      registry.add("hDiffRVtxPion", "hDiffRVtxPion", HistType::kTH1F, {{100, -10, 10}});         // difference between the radius of decay vertex and minR of pion
      registry.add("hDiffRVtxDeuteron", "hDiffRVtxDeuteron", HistType::kTH1F, {{100, -10, 10}}); // difference between the radius of decay vertex and minR of deuteron
    }
  }

  struct Indexdaughters { // check duplicated paired daughters
    int64_t index0;
    int64_t index1;
    int64_t index2;
    bool operator==(const Indexdaughters& t) const
    {
      return (this->index0 == t.index0 && this->index1 == t.index1 && this->index2 == t.index2);
    }
  };

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const&)
  {
    // dummy function
  }
  PROCESS_SWITCH(hypertriton3bodyLabelCheck, processData, "Donot check MC label tables", true);

  void processCheckLabel(soa::Join<aod::Collisions, o2::aod::McCollisionLabels, aod::EvSels>::iterator const& collision, aod::Decay3Bodys const& decay3bodys, soa::Join<aod::Vtx3BodyDatas, aod::McVtx3BodyLabels> const& vtx3bodydatas, MCLabeledTracksIU const& /*tracks*/, aod::McParticles const& /*particlesMC*/, aod::McCollisions const& /*mcCollisions*/)
  {
    // check the decay3body table
    std::vector<Indexdaughters> set_pair;
    for (const auto& d3body : decay3bodys) {
      registry.fill(HIST("hDecay3BodyCounter"), 0.5);
      auto lTrack0 = d3body.track0_as<MCLabeledTracksIU>();
      auto lTrack1 = d3body.track1_as<MCLabeledTracksIU>();
      auto lTrack2 = d3body.track2_as<MCLabeledTracksIU>();
      if (!lTrack0.has_mcParticle() || !lTrack1.has_mcParticle() || !lTrack2.has_mcParticle()) {
        continue;
      }
      auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
      auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
      auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
      if (!lMCTrack0.has_mothers() || !lMCTrack1.has_mothers() || !lMCTrack2.has_mothers()) {
        continue;
      }

      for (const auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (const auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (const auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
            if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
              registry.fill(HIST("hDecay3BodyCounter"), 1.5);
              // duplicated daughters check
              Indexdaughters temp = {lMCTrack0.globalIndex(), lMCTrack1.globalIndex(), lMCTrack2.globalIndex()};
              auto p = std::find(set_pair.begin(), set_pair.end(), temp);
              if (p == set_pair.end()) {
                set_pair.push_back(temp);
                registry.fill(HIST("hDecay3BodyCounter"), 2.5);
                if (lMother0.mcCollisionId() == collision.mcCollisionId()) {
                  registry.fill(HIST("hDecay3BodyCounter"), 3.5);
                  if (lTrack0.collisionId() == lTrack1.collisionId() && lTrack0.collisionId() == lTrack2.collisionId()) {
                    registry.fill(HIST("hDecay3BodyCounter"), 4.5);
                  }
                }
              }
            }
          }
        }
      }
    }

    if (mc_event_selection && (!collision.selection_bit(aod::evsel::kIsTriggerTVX) || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))) {
      return;
    }

    if (event_posZ_selection && std::abs(collision.posZ()) > 10.f) { // 10cm
      return;
    }

    std::vector<int64_t> set_mothertrack;
    for (const auto& vtx : vtx3bodydatas) {
      registry.fill(HIST("hLabeledVtxCounter"), 0.5);
      if (vtx.mcParticleId() != -1) {
        auto mcparticle = vtx.mcParticle_as<aod::McParticles>();
        auto lTrack0 = vtx.track0_as<MCLabeledTracksIU>();
        auto lTrack1 = vtx.track1_as<MCLabeledTracksIU>();
        auto lTrack2 = vtx.track2_as<MCLabeledTracksIU>();
        if (std::abs(mcparticle.pdgCode()) != motherPdgCode) {
          continue;
        }
        registry.fill(HIST("hLabeledVtxCounter"), 1.5);
        registry.fill(HIST("hDiffRVtxDeuteron"), lTrack2.x() - vtx.vtxradius());
        if (mcparticle.pdgCode() > 0) {
          registry.fill(HIST("hHypertritonCounter"), 0.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mHypertriton());
          registry.fill(HIST("hMassTrueH3LMatter"), vtx.mHypertriton());
          registry.fill(HIST("hDiffRVtxProton"), lTrack0.x() - vtx.vtxradius());
          registry.fill(HIST("hDiffRVtxPion"), lTrack1.x() - vtx.vtxradius());
          auto p = std::find(set_mothertrack.begin(), set_mothertrack.end(), mcparticle.globalIndex());
          if (p == set_mothertrack.end()) {
            set_mothertrack.push_back(mcparticle.globalIndex());
            registry.fill(HIST("hLabeledVtxCounter"), 2.5);
          }
          if (std::abs(lTrack0.tpcNSigmaPr()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 0.5);
          }
          if (std::abs(lTrack1.tpcNSigmaPi()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 1.5);
          }
          if (std::abs(lTrack2.tpcNSigmaDe()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 2.5);
          }
          if (std::abs(lTrack0.tpcNSigmaPr()) < tpcPIDNSigmaCut && std::abs(lTrack1.tpcNSigmaPi()) < tpcPIDNSigmaCut && std::abs(lTrack2.tpcNSigmaDe()) < tpcPIDNSigmaCut) {
            registry.fill(HIST("hHypertritonCounter"), 1.5);
          }
        } else {
          registry.fill(HIST("hHypertritonCounter"), 2.5);
          registry.fill(HIST("hMassTrueH3L"), vtx.mAntiHypertriton());
          registry.fill(HIST("hMassTrueH3LAntiMatter"), vtx.mAntiHypertriton());
          registry.fill(HIST("hDiffRVtxProton"), lTrack1.x() - vtx.vtxradius());
          registry.fill(HIST("hDiffRVtxPion"), lTrack0.x() - vtx.vtxradius());
          auto p = std::find(set_mothertrack.begin(), set_mothertrack.end(), mcparticle.globalIndex());
          if (p == set_mothertrack.end()) {
            set_mothertrack.push_back(mcparticle.globalIndex());
            registry.fill(HIST("hLabeledVtxCounter"), 2.5);
          }
          if (std::abs(lTrack0.tpcNSigmaPi()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 4.5);
          }
          if (std::abs(lTrack1.tpcNSigmaPr()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 3.5);
          }
          if (std::abs(lTrack2.tpcNSigmaDe()) > tpcPIDNSigmaCut) {
            registry.fill(HIST("hPIDCounter"), 5.5);
          }
          if (std::abs(lTrack0.tpcNSigmaPi()) < tpcPIDNSigmaCut && std::abs(lTrack1.tpcNSigmaPr()) < tpcPIDNSigmaCut && std::abs(lTrack2.tpcNSigmaDe()) < tpcPIDNSigmaCut) {
            registry.fill(HIST("hHypertritonCounter"), 3.5);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(hypertriton3bodyLabelCheck, processCheckLabel, "Check MC label tables", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyAnalysis>(cfgc),
    adaptAnalysisTask<hypertriton3bodyQa>(cfgc),
    adaptAnalysisTask<hypertriton3bodyLabelCheck>(cfgc),
  };
}
