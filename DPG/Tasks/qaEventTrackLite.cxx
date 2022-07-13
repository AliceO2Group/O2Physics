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
/// \file   qaEventTrackLite.cxx
/// \author Mario Krüger <mario.kruger@cern.ch>
/// \author Mattia Faggin <mattia.faggin@cern.ch>
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \brief  Light version of the task to produce QA objects for the track and the event properties.
///         This task runs on prefiltered data
///

#include "qaEventTrack.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

// Lite version of the QA task to run on skimmed dataset
struct qaEventTrackLite {
  // Binning
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 5.0, 10.0, 20.0, 50.0}, ""};
  ConfigurableAxis bins1overPt{"bins1overPt", {100, 0, 10}, "1/pt binning (c/GeV)"};
  ConfigurableAxis binsImpPar{"binsImpPar", {200, -0.15, 0.15}, "Impact parameter binning (cm)"};
  ConfigurableAxis binsEta{"binsEta", {200, -2., 2.}, "Eta binning"};
  ConfigurableAxis binsPhi{"binsPhi", {180, 0., 2 * M_PI}, "Phi binning"};

  HistogramRegistry histos;

  // Event selections
  Configurable<float> vtxZMax{"vtxZMax", 10.f, "Max VTX Z position"};
  // Track selections
  Configurable<bool> bHasITS{"bHasITS", false, "Select only DPG tracks with ITS (no forcing on TPC information)"};
  Configurable<bool> bItsStandalone{"bItsStandalone", false, "Select only ITS standalone DPG tracks"};
  Configurable<bool> bTpcOnly{"bTpcOnly", false, "Select only TPC only DPG tracks"};
  Configurable<bool> bItsTpcMatched{"bItsTpcMatched", false, "Select ITS-TPC matched DPG tracks"};
  // Kinematic selections
  Configurable<float> ptMin{"ptMin", 0., "Minimum track pt"};
  Configurable<float> etaMin{"etaMin", -10., "Minimum eta for DPG tracks"};
  Configurable<float> etaMax{"etaMax", 10., "Maximum eta for DPG tracks"};
  // ITS selections
  Configurable<float> chi2ItsMin{"chi2ItsMin", -1001.f, "Max ITS chi2"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 1000.f, "Max ITS chi2"};
  // TPC selections
  Configurable<int> nClusterTpcMin{"nClusterTpcMin", -1001, "Minimum number of TPC clusters"};
  Configurable<int> nCrossedRowsTpcMin{"nCrossedRowsTpcMin", -1001, "Minimum number of TPC crossed rows"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMin{"nCrossedRowsTpcOverFindableClustersTpcMin", -1, "Minimum ratio between TPC crossed rows and findable clusters"};
  Configurable<float> chi2TpcMin{"chi2TpcMin", -1001.f, "Max TPC chi2"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 1000.f, "Max TPC chi2"};
  // TOF selections
  Configurable<float> chi2TofMin{"chi2TofMin", -1001.f, "Max TOF chi2"};
  Configurable<float> lengthMin{"lengthMin", -1001.f, "Min length"};
  Configurable<float> dcaXYmax{"dcaXYMax", 999., "Max dca XY"};

  // Selection 1
  // ITS selections
  Configurable<float> chi2ItsMaxSel1{"chi2ItsMaxSel1", 1000.f, "Max ITS chi2 Sel1"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel1{"nClusterTpcMinSel1", -1001, "Minimum number of TPC clusters Sel1"};
  Configurable<int> nCrossedRowsTpcMinSel1{"nCrossedRowsTpcMinSel1", -1001, "Minimum number of TPC crossed rows Sel1"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel1{"nCrossedRowsTpcOverFindableClustersTpcMinSel1", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel1"};
  Configurable<float> chi2TpcMaxSel1{"chi2TpcMaxSel1", 1000.f, "Max TPC chi2 Sel1"};
  Configurable<bool> bItsTpcMatchedSel1{"bItsTpcMatchedSel1", false, "Select ITS-TPC matched DPG tracks, sel1"};
  Configurable<float> dcaXYmaxSel1{"dcaXYMaxSel1", 999., "Max dca XY sel1"};

  // Selection 2
  // ITS selections
  Configurable<float> chi2ItsMaxSel2{"chi2ItsMaxSel2", 1000.f, "Max ITS chi2 Sel2"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel2{"nClusterTpcMinSel2", -1001, "Minimum number of TPC clusters Sel2"};
  Configurable<int> nCrossedRowsTpcMinSel2{"nCrossedRowsTpcMinSel2", -1001, "Minimum number of TPC crossed rows Sel2"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel2{"nCrossedRowsTpcOverFindableClustersTpcMinSel2", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel2"};
  Configurable<float> chi2TpcMaxSel2{"chi2TpcMaxSel2", 1000.f, "Max TPC chi2 Sel2"};
  Configurable<bool> bItsTpcMatchedSel2{"bItsTpcMatchedSel2", false, "Select ITS-TPC matched DPG tracks, sel2"};
  Configurable<float> dcaXYmaxSel2{"dcaXYMaxSel2", 999., "Max dca XY sel2"};

  // Selection 3
  // ITS selections
  Configurable<float> chi2ItsMaxSel3{"chi2ItsMaxSel3", 1000.f, "Max ITS chi2 Sel3"};
  // TPC selections
  Configurable<int> nClusterTpcMinSel3{"nClusterTpcMinSel3", -1001, "Minimum number of TPC clusters Sel3"};
  Configurable<int> nCrossedRowsTpcMinSel3{"nCrossedRowsTpcMinSel3", -1001, "Minimum number of TPC crossed rows Sel3"};
  Configurable<float> nCrossedRowsTpcOverFindableClustersTpcMinSel3{"nCrossedRowsTpcOverFindableClustersTpcMinSel3", -1, "Minimum ratio between TPC crossed rows and findable clusters Sel3"};
  Configurable<float> chi2TpcMaxSel3{"chi2TpcMaxSel3", 1000.f, "Max TPC chi2 Sel3"};
  Configurable<bool> bItsTpcMatchedSel3{"bItsTpcMatchedSel3", false, "Select ITS-TPC matched DPG tracks, sel3"};
  Configurable<float> dcaXYmaxSel3{"dcaXYMaxSel3", 999., "Max dca XY sel3"};

  // MC selections
  Configurable<int> pdgCodeSel{"pdgCodeSel", 2, "pdgCode based particle selection, 1 defines pi,K,p,mu,e, 2 all final-state charged particles including light (hyper)nuclei"};
  Configurable<bool> checkPdgAtReco{"checkPdgAtReco", false, "check pdg code also at reco levo for data-like reference"};

  void init(InitContext const&)
  {
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/c)"};
    const AxisSpec axis1overPt{bins1overPt, "1/#it{p}_{T} (GeV/c)^{-1}"};
    const AxisSpec axisEta{binsEta, "#it{#eta}"};
    const AxisSpec axisPhi{binsPhi, "#it{#phi} (rad)"};

    // kine histograms
    histos.add("Tracks/VertexPositionZ", "", kTH1D, {{100, -20.f, 20.f, "Vertex Z (cm)"}});
    histos.add("Tracks/VertexPositionZvsEta", "", kTH2D, {axisEta, {100, -20.f, 20.f, "Vertex Z (cm)"}});
    histos.add("Tracks/VertexPositionZvsEtaHasITS", "", kTH2D, {axisEta, {100, -20.f, 20.f, "Vertex Z (cm)"}});
    histos.add("Tracks/Kine/pt", "#it{p}_{T}", kTH1D, {axisPt});
    histos.add("Tracks/Kine/eta", "#eta", kTH1D, {axisEta});
    histos.add("Tracks/Kine/phi", "#phi", kTH1D, {axisPhi});
    histos.add("Tracks/Kine/phiVsEta", "#phi", kTH2D, {axisPhi, axisEta});
    histos.add("Tracks/length", "track length in cm;#it{Length} (cm);", kTH1D, {{400, 0, 1000}});
    const AxisSpec axisImpParRPhi{binsImpPar, "#it{d}_{r#it{#varphi}} (#cm)"};
    const AxisSpec axisImpParZAxis{binsImpPar, "#it{d}_{z} (#cm)"};
    histos.add("Tracks/dcaXY", "distance of closest approach in #it{xy} plane", kTH1D, {axisImpParRPhi});
    histos.add("Tracks/dcaZ", "distance of closest approach in #it{z}", kTH1D, {axisImpParZAxis});
    histos.add("Tracks/dcaXYvsPt", "d_#it{xy} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});
    histos.add("Tracks/dcaZvsPt", "d_#it{z} vs. #it{p}_{T}", kTH2D, {axisPt, axisImpParRPhi});
    histos.add("Tracks/dcaXYvsPtvsEta", "d_#it{xy} vs. #it{p}_{T} vs. #eta", kTH3D, {axisPt, axisEta, axisImpParRPhi});
    histos.add("Tracks/dcaZvsPtvsEta", "d_#it{z} vs. #it{p}_{T} vs. #eta", kTH3D, {axisPt, axisEta, axisImpParRPhi});

    // its histograms
    histos.add("Tracks/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {{100, 0, 40}});
    histos.add("Tracks/ITS/itsNCl", "ITS number of clusters;# clusters ITS", kTH1D, {{8, -0.5, 7.5}});
    histos.add("Tracks/ITS/itsNClvsItsHitmap", "ITS number of clusters vs. ITS hitmap;# clusters ITS; ITS hitmap", kTH2D, {{8, -0.5, 7.5, "# clusters ITS"}, {128, 0, 128, "ITS hitmap"}});
    histos.add("Tracks/ITS/itsNClstvsEtavsPt", "profile2D;", kTProfile2D, {axisEta, axisPt});
    // tpc histograms
    histos.add("Tracks/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {{100, 0, 10}});
    histos.add("Tracks/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {{165, -0.5, 164.5}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {{60, 0.7, 1.3}});
    histos.add("Tracks/TPC/tpcNClsFoundvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# clusters TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsvsPt", "", kTH2D, {axisPt, {165, -0.5, 164.5, "# crossed rows TPC"}});
    histos.add("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt", "", kTH2D, {axisPt, {60, 0.7, 1.3, "crossed rows / findable clusters TPC"}});
    histos.add("Tracks/TPC/TPCnClstvsEtavsPt", "profile2D;", kTProfile2D, {axisEta, axisPt});
    histos.add("Tracks/TPC/dEdxvsP", "", kTH2D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {500, 0, 1000, "d#it{E}/d#it{x} (a.u.)"}});
    histos.add("Tracks/TPC/dEdxvsPvsEta", "", kTH3D, {{5000, 0, 10, "#it{p} (GeV/#it{c})"}, {20, -2, 2, "#it{#eta}"}, {500, 0, 1000, "d#it{E}/d#it{x} (a.u.)"}});
    // trd histograms
    histos.add("Tracks/TRD/trdChi2", "chi2 in TRD", kTH1D, {{100, 0, 10, "chi2 / cluster TRD"}});
    // tof histograms
    histos.add("Tracks/TOF/tofChi2", "chi2 in TOF", kTH1D, {{100, 0, 10, "chi2 / cluster TOF"}});
    // matching histogram
    histos.add("Tracks/matchedDet", "matched detectors", kTH1D, {{4, 0.5, 4.5, ""}});
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(1, "hasTPC");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(2, "hasITS");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(3, "hasTRD");
    histos.get<TH1>(HIST("Tracks/matchedDet"))->GetXaxis()->SetBinLabel(4, "hasTOF");
    // kinematics
    histos.add("Tracks/relativeResoPt", "relative #it{p}_{T} resolution;#sigma(#it{p}_{T})/#it{p}_{T};#it{p}_{T};", kTH2D, {{axisPt, {500, 0., 1, "#sigma{#it{p}}/#it{p}_{T}"}}});
    histos.add("Tracks/relativeResoPtMean", "mean relative #it{p}_{T} resolution;;average #sigma_{#it{p}_{T}}/#it{p}_{T};", kTProfile, {axisPt});
    histos.add("Tracks/relativeResoPtMeanvsEtavsPt", "mean relative #it{p}_{T} resolution;", kTProfile2D, {axisPt, axisEta});
    histos.add("Tracks/reso1overPtMeanvsEtavs1overPt", "mean 1/#it{p}_{T} resolution;", kTProfile2D, {axis1overPt, axisEta});
    // track cuts map
    histos.add("TrackCuts/selPtEtaPhiNoSel", "pt eta phi map no sel; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
    histos.add("TrackCuts/selPtEtaPhiSel1", "pt eta phi map sel 1; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
    histos.add("TrackCuts/selPtEtaPhiSel2", "pt eta phi map sel 2; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
    histos.add("TrackCuts/selPtEtaPhiSel3", "pt eta phi map sel 3; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});

    // MC histograms
    if (doprocessMCLite) {
      histos.add("Particles/PDGs", "Particle PDGs;PDG Code", kTH1D, {{100, 0.f, 100.f}});
      histos.add("Particles/Kine/pt", "Particle #it{p}_{T}", kTH1D, {axisPt});
      histos.add("Particles/Kine/eta", "Particle #eta", kTH1D, {axisEta});
      histos.add("Particles/Kine/phi", "Particle #phi", kTH1D, {axisPhi});
      histos.add("Particle/selPtEtaPhiMCGenPrimary", "pt eta phi map MC gen Primary; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
      histos.add("Particle/selPtEtaPhiMCRecoNoSelPrimary", "pt eta phi map MC gen Primary sel0; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
      histos.add("Particle/selPtEtaPhiMCRecoSel1Primary", "pt eta phi map MC gen Primary sel1; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
      histos.add("Particle/selPtEtaPhiMCRecoSel2Primary", "pt eta phi map MC gen Primary sel2; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
      histos.add("Particle/selPtEtaPhiMCRecoSel3Primary", "pt eta phi map MC gen Primary sel3; pt,eta,phi", kTH3D, {axisPt, {50, -1.2, 1.2, "#eta"}, {30, 0., 2 * M_PI, "#varphi"}});
      histos.add("Tracks/resoPhivsPtvsEta", "#varphi(reco)-#varphi(gen);", kTH3D, {axisPt, axisEta, {180, -M_PI, M_PI, "#varphi(reco)-#varphi(gen)"}});
      histos.add("Tracks/phiRecovsphiGen", "#varphi(reco) vs. #varphi(gen);", kTH2D, {axisPhi, axisPhi});
      histos.get<TH2>(HIST("Tracks/phiRecovsphiGen"))->GetXaxis()->SetTitle("#varphi(reco)");
      histos.get<TH2>(HIST("Tracks/phiRecovsphiGen"))->GetYaxis()->SetTitle("#varphi(gen)");
    }
  }

  ///////////////
  /// Filters ///
  ///////////////
  // Kinematics
  Filter ptCut = o2::aod::dpgtrack::pt > ptMin;
  Filter etaCut = etaMin < o2::aod::dpgtrack::eta && o2::aod::dpgtrack::eta < etaMax;
  // Detector matching
  Filter filterHasITS = (bHasITS.node() == false) || (bItsStandalone.node() == false && bTpcOnly.node() == false && bItsTpcMatched.node() == false && o2::aod::dpgtrack::hasITS == true);
  Filter itsStandaloneTracks = (bItsStandalone.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == false);
  Filter tpcOnlyTracks = (bTpcOnly.node() == false) || (o2::aod::dpgtrack::hasITS == false && o2::aod::dpgtrack::hasTPC == true);
  Filter itsTpcMatchedTracks = (bItsTpcMatched.node() == false) || (o2::aod::dpgtrack::hasITS == true && o2::aod::dpgtrack::hasTPC == true);
  // ITS
  Filter itsChi2 = (bTpcOnly.node() == true) || (chi2ItsMin < o2::aod::track::itsChi2NCl && o2::aod::track::itsChi2NCl < chi2ItsMax);
  // TPC
  Filter tpcChi2s = (bItsStandalone.node() == true) || (chi2TpcMin < o2::aod::track::tpcChi2NCl && o2::aod::track::tpcChi2NCl < chi2TpcMax);
  Filter tpcNclusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsFound > (int16_t)nClusterTpcMin);
  Filter tpcNcrossedRows = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcNClsCrossedRows > (int16_t)nCrossedRowsTpcMin);
  Filter tpcNcrossedRowsOverFindableClusters = (bItsStandalone.node() == true) || (o2::aod::dpgtrack::tpcCrossedRowsOverFindableCls > nCrossedRowsTpcOverFindableClustersTpcMin);
  // TOF
  Filter tofChi = o2::aod::track::tofChi2 > chi2TofMin;
  Filter length = o2::aod::track::length > lengthMin;

  template <bool isMC, typename trackType>
  void fillHistograms(trackType const& track)
  {

    if (TMath::Abs(track.dpgCollision().posZ()) > vtxZMax)
      return;
    if constexpr (isMC) {
      if (track.productionMode() == 0 && isPdgSelected(track.pdgCode())) {
        histos.fill(HIST("Particle/selPtEtaPhiMCGenPrimary"), track.ptMC(), track.etaMC(), track.phiMC());
      }
    }
    // temporary additional selections
    if (!applyTrackSelectionsNotFiltered(track)) {
      return;
    }
    Bool_t sel1 = applyTrackSelectionsSel1(track);
    Bool_t sel2 = applyTrackSelectionsSel2(track);
    Bool_t sel3 = applyTrackSelectionsSel3(track);
    Bool_t isPdgOk = true;
    if constexpr (isMC) {
      isPdgOk = isPdgSelected(track.pdgCode());
      if (track.productionMode() == 0 && isPdgOk) {
        histos.get<TH1>(HIST("Particles/PDGs"))->Fill(Form("%i", track.pdgCode()), 1);
        histos.fill(HIST("Particle/selPtEtaPhiMCRecoNoSelPrimary"), track.ptMC(), track.etaMC(), track.phiMC());

        if (sel1) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel1Primary"), track.ptMC(), track.etaMC(), track.phiMC());
        }
        if (sel2) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel2Primary"), track.ptMC(), track.etaMC(), track.phiMC());
        }
        if (sel3) {
          histos.fill(HIST("Particle/selPtEtaPhiMCRecoSel3Primary"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      }

      histos.fill(HIST("Particles/Kine/pt"), track.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), track.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), track.phiMC());
    }
    histos.fill(HIST("Tracks/VertexPositionZ"), track.dpgCollision().posZ());
    histos.fill(HIST("Tracks/VertexPositionZvsEta"), track.eta(), track.dpgCollision().posZ());
    if (track.hasITS()) {
      histos.fill(HIST("Tracks/VertexPositionZvsEtaHasITS"), track.eta(), track.dpgCollision().posZ());
    }
    histos.fill(HIST("Tracks/Kine/pt"), track.pt());
    histos.fill(HIST("Tracks/Kine/eta"), track.eta());
    histos.fill(HIST("Tracks/Kine/phi"), track.phi());
    histos.fill(HIST("Tracks/Kine/phiVsEta"), track.phi(), track.eta());
    histos.fill(HIST("Tracks/dcaXY"), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZ"), track.dcaZ());
    histos.fill(HIST("Tracks/dcaXYvsPt"), track.pt(), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZvsPt"), track.pt(), track.dcaZ());
    histos.fill(HIST("Tracks/dcaXYvsPtvsEta"), track.pt(), track.eta(), track.dcaXY());
    histos.fill(HIST("Tracks/dcaZvsPtvsEta"), track.pt(), track.eta(), track.dcaZ());
    histos.fill(HIST("Tracks/length"), track.length());
    histos.fill(HIST("Tracks/ITS/itsChi2NCl"), track.itsChi2NCl());
    histos.fill(HIST("Tracks/ITS/itsNCl"), track.itsNCls());
    histos.fill(HIST("Tracks/ITS/itsNClvsItsHitmap"), track.itsNCls(), track.itsClusterMap());
    histos.fill(HIST("Tracks/TPC/tpcChi2NCl"), track.tpcChi2NCl());
    histos.fill(HIST("Tracks/TPC/tpcNClsFound"), track.tpcNClsFound());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRows"), track.tpcNClsCrossedRows());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableCls"), track.tpcCrossedRowsOverFindableCls());
    histos.fill(HIST("Tracks/TPC/tpcNClsFoundvsPt"), track.pt(), track.tpcNClsFound());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsvsPt"), track.pt(), track.tpcNClsCrossedRows());
    histos.fill(HIST("Tracks/TPC/tpcCrossedRowsOverFindableClsvsPt"), track.pt(), track.tpcCrossedRowsOverFindableCls());
    histos.fill(HIST("Tracks/TPC/dEdxvsP"), track.pt() / (sin(2 * atan2(1, exp(track.eta())))), track.tpcSignal());
    histos.fill(HIST("Tracks/TPC/dEdxvsPvsEta"), track.pt() / (sin(2 * atan2(1, exp(track.eta())))), track.eta(), track.tpcSignal());
    histos.fill(HIST("Tracks/TRD/trdChi2"), track.trdChi2());
    histos.fill(HIST("Tracks/TOF/tofChi2"), track.tofChi2());
    if (track.hasTPC()) {
      histos.fill(HIST("Tracks/matchedDet"), 1);
    }
    if (track.hasITS()) {
      histos.fill(HIST("Tracks/matchedDet"), 2);
    }
    if (track.hasTRD()) {
      histos.fill(HIST("Tracks/matchedDet"), 3);
    }
    if (track.hasTOF()) {
      histos.fill(HIST("Tracks/matchedDet"), 4);
    }
    histos.fill(HIST("Tracks/relativeResoPt"), track.pt(), track.ptReso());
    histos.fill(HIST("Tracks/relativeResoPtMean"), track.pt(), track.ptReso());
    histos.fill(HIST("Tracks/relativeResoPtMeanvsEtavsPt"), track.pt(), track.eta(), track.ptReso());
    histos.fill(HIST("Tracks/reso1overPtMeanvsEtavs1overPt"), 1. / track.pt(), track.eta(), track.ptReso() / (track.pt())); // reso(1/pt) = reso(pt)/pt^2 and here ptReso==reso(pt)/pt
    histos.fill(HIST("Tracks/TPC/TPCnClstvsEtavsPt"), track.eta(), track.pt(), track.tpcNClsFound());
    histos.fill(HIST("Tracks/ITS/itsNClstvsEtavsPt"), track.eta(), track.pt(), track.itsNCls());
    if constexpr (isMC) {
      histos.fill(HIST("Tracks/resoPhivsPtvsEta"), track.pt(), track.eta(), track.phi() - track.phiMC());
      histos.fill(HIST("Tracks/phiRecovsphiGen"), track.phi(), track.phiMC());
    }

    if constexpr (isMC) {
      if (isPdgOk || !checkPdgAtReco) {
        histos.fill(HIST("TrackCuts/selPtEtaPhiNoSel"), track.ptMC(), track.etaMC(), track.phiMC());
      }
    } else {
      histos.fill(HIST("TrackCuts/selPtEtaPhiNoSel"), track.pt(), track.eta(), track.phi());
    }

    if (sel1) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel1"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel1"), track.pt(), track.eta(), track.phi());
      }
    }
    if (sel2) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel2"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel2"), track.pt(), track.eta(), track.phi());
      }
    }
    if (sel3) {
      if constexpr (isMC) {
        if (isPdgOk || !checkPdgAtReco) {
          histos.fill(HIST("TrackCuts/selPtEtaPhiSel3"), track.ptMC(), track.etaMC(), track.phiMC());
        }
      } else {
        histos.fill(HIST("TrackCuts/selPtEtaPhiSel3"), track.pt(), track.eta(), track.phi());
      }
    }
  }

  // Process data
  void processDataLite(o2::soa::Filtered<aod::DPGTracks> const& tracks, aod::DPGCollisions const&)
  {
    for (const auto& track : tracks) {
      fillHistograms<false>(track);
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processDataLite, "process data lite", true);

  // Process MC
  void processMCLite(o2::soa::Filtered<soa::Join<aod::DPGTracks, aod::DPGRecoParticles>> const& tracks, aod::DPGCollisions const&, aod::DPGNonRecoParticles const& particles)
  {
    for (const auto& track : tracks) {
      fillHistograms<true>(track);
    }

    for (const auto& particle : particles) {
      if (TMath::Abs(particle.dpgCollision().posZ()) > vtxZMax)
        continue;

      histos.fill(HIST("Particles/Kine/pt"), particle.ptMC());
      histos.fill(HIST("Particles/Kine/eta"), particle.etaMC());
      histos.fill(HIST("Particles/Kine/phi"), particle.phiMC());

      if (particle.productionMode() == 0 && isPdgSelected(particle.pdgCode()))
        histos.fill(HIST("Particle/selPtEtaPhiMCGenPrimary"), particle.ptMC(), particle.etaMC(), particle.phiMC());
    }
  }
  PROCESS_SWITCH(qaEventTrackLite, processMCLite, "process MC lite", false);

  template <typename T>
  bool applyTrackSelectionsNotFiltered(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMin)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMin)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMin)
      return false;
    if (bItsStandalone == true && (track.hasITS() == false || track.hasTPC() == true))
      return false;
    if (bTpcOnly == true && (track.hasITS() == true || track.hasTPC() == false))
      return false;
    if (bItsTpcMatched == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.itsChi2NCl() > chi2ItsMax)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMax)
      return false;
    if (track.tpcNClsFound() < nClusterTpcMin)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmax)
      return false;
    return true;
  }

  template <typename T>
  bool applyTrackSelectionsSel1(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel1)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel1)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel1)
      return false;
    if (bItsTpcMatchedSel1 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel1)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel1)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel1)
      return false;

    return true;
  }
  template <typename T>
  bool applyTrackSelectionsSel2(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel2)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel2)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel2)
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel2)
      return false;
    if (bItsTpcMatchedSel2 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel2)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel2)
      return false;

    return true;
  }
  template <typename T>
  bool applyTrackSelectionsSel3(const T& track)
  {
    if (track.tpcNClsFound() < nClusterTpcMinSel3)
      return false;
    if (track.tpcNClsCrossedRows() < nCrossedRowsTpcMinSel3)
      return false;
    if (bItsTpcMatchedSel3 == true && (track.hasITS() == false || track.hasTPC() == false))
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < nCrossedRowsTpcOverFindableClustersTpcMinSel3)
      return false;
    if (track.itsChi2NCl() > chi2ItsMaxSel3)
      return false;
    if (track.tpcChi2NCl() > chi2TpcMaxSel3)
      return false;
    if (TMath::Abs(track.dcaXY()) > dcaXYmaxSel3)
      return false;

    return true;
  }

  bool isPdgSelected(const Int_t pdgcode)
  { // mimics selection of charged particles or id particles
    Int_t abspdgcode = TMath::Abs(pdgcode);
    if (abspdgcode == pdgCodeSel)
      return true;
    if (pdgCodeSel == 1 || pdgCodeSel == 2) {
      if (abspdgcode == 211 || abspdgcode == 321 || abspdgcode == 2212 || abspdgcode == 11 || abspdgcode == 13)
        return true;
      if (pdgCodeSel == 2) {
        if (abspdgcode == 3222 || abspdgcode == 3112 || abspdgcode == 3312 || abspdgcode == 3334 || abspdgcode == 1000010020 || abspdgcode == 1000010030 || abspdgcode == 1000020030 || abspdgcode == 1000020040 || abspdgcode == 1010010030 || abspdgcode == 1010020040)
          return true;
      }
    }
    return false;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaEventTrackLite>(cfgc)};
}
