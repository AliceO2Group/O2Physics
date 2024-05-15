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
/// \brief  task for the analysis of rho photoproduction in UPCs
/// \author Jakub Juracka, jakub.juracka@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"

#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h" // has some useful funtions for stuff not available from the tables

// ROOT headers
#include <Math/Vector4D.h> // this should apparently be used instead of TLorentzVector
//#include "TEfficiency.h" // for eventual MC studies

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDZdcsReduced>::iterator;
using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;

struct upcRhoAnalysis {
  // configurables
  double PcEtaCut = 0.9; // cut on track eta as per physics coordination recommendation
  Configurable<bool> tracksRequireTOF{"tracksRequireTOF", false, "requireTOF"};
  Configurable<double> tracksTpcNSigmaPiCut{"treacksTpcNSigmaPiCut", 3.0, "tpcNSigmaPiCut"};
  Configurable<double> tracksTpcNSigmaElCut{"treacksTpcNSigmaElCut", 3.0, "tpcNSigmaElCut"};
  Configurable<double> tracksTofNSigmaPiCut{"treacksTofNSigmaPiCut", 3.0, "tofNSigmaPiCut"};
  Configurable<double> tracksPtMaxCut{"tracksPtMaxCut", 2.0, "ptMaxCut"};
  Configurable<double> tracksDcaMaxCut{"tracksDcaMaxCut", 1.0, "dcaMaxCut"};
  Configurable<double> collisionsPosZMaxCut{"collisionsPosZMaxCut", 10.0, "posZMaxCut"};

  ConfigurableAxis mAxis{"mAxis", {800, 0.0, 8.0}, "m (GeV/#it{c}^{2})"};
  ConfigurableAxis ptAxis{"ptAxis", {500, 0.0, 5.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis etaAxis{"etaAxis", {180, -0.9, 0.9}, "#eta"};
  ConfigurableAxis yAxis{"yAxis", {180, -0.9, 0.9}, "y"};
  ConfigurableAxis phiAxis{"phiAxis", {180, 0.0, 2.0*o2::constants::math::PI}, "#phi"};
  ConfigurableAxis nTracksAxis{"nTracksAxis", {101, -0.5, 100.5}, "N_{tracks}"};
  ConfigurableAxis tpcNSigmaPiAxis{"tpcNSigmaPiAxis", {400, -10.0, 30.0}, "TPC n#sigma_{#pi}"};
  ConfigurableAxis tofNSigmaPiAxis{"tofNSigmaPiAxis", {400, -20.0, 20.0}, "TOF n#sigma_{#pi}"};
  ConfigurableAxis dcaAxis{"dcaXYAxis", {1000, -5.0, 5.0}, "DCA (cm)"};

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&){
    // selection counter
    std::vector<std::string> selectionNames = {"all tracks", "PV tracks", "DCA cut", "ITS required", "#eta cut", "p_{T} cut", "2 tracks in event", "opposite charge"};
    int nSelections = selectionNames.size();
    registry.add("hSelectionCounter", ";;counts", kTH1D, {{nSelections, 0.5, nSelections + 0.5}});
    for (int i = 0; i < nSelections; i++) registry.get<TH1>(HIST("hSelectionCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionNames[i].c_str());
    registry.get<TH1>(HIST("hSelectionCounter"))->SetStats(0);
    // all collisions
    registry.add("QC/collisions/uncut/hNetCharge", ";net charge;counts", kTH1D, {{11, -5.5, 5.5}});
    registry.add("QC/collisions/uncut/hNumContributors", ";N_{contributors};counts", kTH1D, {{11, -0.5, 10.5}});
    registry.add("QC/collisions/uncut/hRgtrwTOF", ";fraction of tracks with TOF;counts", kTH1D, {{101, -0.005, 1.005}});
    registry.add("QC/collisions/uncut/hFt0Amplitude", ";A;C;counts", kTH2D, {{201, -1.5, 200.5}, {201, -1.5, 200.5}});
    registry.add("QC/collisions/uncut/hFt0Time", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
    registry.add("QC/collisions/uncut/hFddAmplitude", ";A;C;counts", kTH2D, {{300, 0.0, 3000.0}, {300, 0.0, 3000.0}});
    registry.add("QC/collisions/uncut/hFddTime", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
    registry.add("QC/collisions/uncut/hFv0Amplitude", ";A;counts", kTH1D, {{72, -1.5, 70.5}});
    registry.add("QC/collisions/uncut/hFv0Time", ";t;counts", kTH1D, {{1050, -999.5, 50.5}});
    registry.add("QC/collisions/uncut/hPosX", ";x (cm);counts", kTH1D, {{200, -1.0, 1.0}});
    registry.add("QC/collisions/uncut/hPosY", ";y (cm);counts", kTH1D, {{200, -1.0, 1.0}});
    registry.add("QC/collisions/uncut/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/collisions/uncut/hZDCcommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0},{250, -5.0, 20.0}});
    registry.add("QC/collisions/uncut/hZDCtime", ";ZNA time;ZNC time;counts", kTH2D, {{200, -10.0, 10.0},{200, -10.0, 10.0}});
    // collisions passing cuts
    registry.add("QC/collisions/cut/hNetCharge", ";net charge;counts", kTH1D, {{11, -5.5, 5.5}});
    registry.add("QC/collisions/cut/hNumContributors", ";N_{contributors};counts", kTH1D, {{11, -0.5, 10.5}});
    registry.add("QC/collisions/cut/hRgtrwTOF", ";fraction of tracks with TOF;counts", kTH1D, {{101, -0.005, 1.005}});
    registry.add("QC/collisions/cut/hFt0Amplitude", ";A;C;counts", kTH2D, {{201, -1.5, 200.5}, {201, -1.5, 200.5}});
    registry.add("QC/collisions/cut/hFt0Time", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
    registry.add("QC/collisions/cut/hFddAmplitude", ";A;C;counts", kTH2D, {{300, 0.0, 3000.0}, {300, 0.0, 3000.0}});
    registry.add("QC/collisions/cut/hFddTime", ";t_{A};t_{C};counts", kTH2D, {{1050, -999.5, 50.5}, {1050, -999.5, 50.5}});
    registry.add("QC/collisions/cut/hFv0Amplitude", ";A;counts", kTH1D, {{72, -1.5, 70.5}});
    registry.add("QC/collisions/cut/hFv0Time", ";t;counts", kTH1D, {{1050, -999.5, 50.5}});
    registry.add("QC/collisions/cut/hPosX", ";x (cm);counts", kTH1D, {{200, -1.0, 1.0}});
    registry.add("QC/collisions/cut/hPosY", ";y (cm);counts", kTH1D, {{200, -1.0, 1.0}});
    registry.add("QC/collisions/cut/hPosZ", ";z (cm);counts", kTH1D, {{400, -20.0, 20.0}});
    registry.add("QC/collisions/cut/hZDCcommonEnergy", ";ZNA common energy;ZNC common energy;counts", kTH2D, {{250, -5.0, 20.0},{250, -5.0, 20.0}});
    registry.add("QC/collisions/cut/hZDCtime", ";ZNA time;ZNC time;counts", kTH2D, {{200, -10.0, 10.0},{200, -10.0, 10.0}});
    // all tracks
    registry.add("QC/allTracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("QC/allTracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
    registry.add("QC/allTracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("QC/allTracks/hCharge", ";charge;counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("QC/allTracks/hNTracks", ";N_{tracks};counts", kTH1D, {nTracksAxis});
    registry.add("QC/allTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
    registry.add("QC/allTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
    registry.add("QC/allTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
    registry.add("QC/allTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {100, 0.0, 250.0}});
    registry.add("QC/allTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi^{+}};TPC n#sigma_{#pi^{-}};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
    registry.add("QC/allTracks/hTpcNClsFindable", ";N_{findable};counts", kTH1D, {{171, -0.5, 170.5}});
    registry.add("QC/allTracks/hTpcNClsFindableMinusFound", ";N_{findable} - N_{found};counts", kTH1D, {{11, -0.5, 10.5}});
    // tracks passing selections
    registry.add("QC/cutTracks/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("QC/cutTracks/hEta", ";#eta;counts", kTH1D, {etaAxis});
    registry.add("QC/cutTracks/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("QC/cutTracks/hCharge", ";charge;counts", kTH1D, {{5, -2.5, 2.5}});
    registry.add("QC/cutTracks/hNTracks", ";N_{tracks};counts", kTH1D, {nTracksAxis});
    registry.add("QC/cutTracks/hTpcNSigmaPi", ";TPC n#sigma_{#pi};counts", kTH1D, {tpcNSigmaPiAxis});
    registry.add("QC/cutTracks/hTofNSigmaPi", ";TOF n#sigma_{#pi};counts", kTH1D, {tofNSigmaPiAxis});
    registry.add("QC/cutTracks/hDcaXYZ", ";DCA_{z} (cm);DCA_{xy} (cm);counts", kTH2D, {dcaAxis, dcaAxis});
    registry.add("QC/cutTracks/hTpcSignalVsPt", ";p_{T} (GeV/#it{c});TPC signal;counts", kTH2D, {ptAxis, {500, 0.0, 500.0}});
    registry.add("QC/cutTracks/hTpcNsigmaPi2D", ";TPC n#sigma_{#pi^{+}};TPC n#sigma_{#pi^{-}};counts", kTH2D, {tpcNSigmaPiAxis, tpcNSigmaPiAxis});
    registry.add("QC/cutTracks/hTpcNClsFindable", ";N_{findable};counts", kTH1D, {{171, -0.5, 170.5}});
    registry.add("QC/cutTracks/hTpcNClsFindableMinusFound", ";N_{findable} - N_{found};counts", kTH1D, {{11, -0.5, 10.5}});
    // reco pions
    registry.add("reco/pions/unlike-sign/hPt", ";p_{T}(#pi^{+}) (GeV/#it{c});p_{T}(#pi^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/unlike-sign/hEta", ";#eta(#pi^{+});#eta(#pi^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/unlike-sign/hPhi", ";#phi(#pi^{+});#phi(#pi^{-});counts", kTH2D, {phiAxis, phiAxis});
    registry.add("reco/pions/like-sign/hPt", ";p_{T}(#pi_{1}) (GeV/#it{c});p_{T}(#pi_{2}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/pions/like-sign/hEta", ";#eta(#pi_{1});#eta(#pi_{2});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/pions/like-sign/hPhi", ";#phi(#pi_{1});#phi(#pi_{2});counts", kTH2D, {phiAxis, phiAxis});
    // reco electron-positrons
    registry.add("reco/electrons/hPt", ";p_{T}(e^{+}) (GeV/#it{c});p_{T}(e^{-}) (GeV/#it{c});counts", kTH2D, {ptAxis, ptAxis});
    registry.add("reco/electrons/hEta", ";#eta(e^{+});#eta(e^{-});counts", kTH2D, {etaAxis, etaAxis});
    registry.add("reco/electrons/hPhi", ";#phi(e^{+});#phi(e^{-});counts", kTH2D, {phiAxis, phiAxis});
    // reco rhos
    registry.add("reco/system/2pi/unlike-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/unlike-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/unlike-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/unlike-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/unlike-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/2pi/like-sign/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2pi/like-sign/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2pi/like-sign/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2pi/like-sign/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2pi/like-sign/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    // el-pos system
    registry.add("reco/system/2el/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/2el/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/2el/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/2el/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/2el/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    // 4pi system
    registry.add("reco/system/4pi/net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/4pi/net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/4pi/net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/4pi/net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/4pi/net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/4pi/non-net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/4pi/non-net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/4pi/non-net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/4pi/non-net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/4pi/non-net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    // 6pi system
    registry.add("reco/system/6pi/net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/6pi/net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/6pi/net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/6pi/net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/6pi/net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
    registry.add("reco/system/6pi/non-net-zero/hM", ";m (GeV/#it{c}^{2});counts", kTH1D, {mAxis});
    registry.add("reco/system/6pi/non-net-zero/hPt", ";p_{T} (GeV/#it{c});counts", kTH1D, {ptAxis});
    registry.add("reco/system/6pi/non-net-zero/hPtVsM", ";m (GeV/#it{c}^{2});p_{T} (GeV/#it{c});counts", kTH2D, {mAxis, ptAxis});
    registry.add("reco/system/6pi/non-net-zero/hY", ";y;counts", kTH1D, {yAxis});
    registry.add("reco/system/6pi/non-net-zero/hPhi", ";#phi;counts", kTH1D, {phiAxis});
  }

  template <typename T>
  bool trackPassesCuts(T const& track) { // check if track passes preliminary cuts (PID done separately)
    registry.fill(HIST("hSelectionCounter"), 1);
    if (!track.isPVContributor()) return false;
    registry.fill(HIST("hSelectionCounter"), 2);
    if (track.dcaZ() > tracksDcaMaxCut || track.dcaXY() > tracksDcaMaxCut) return false;
    registry.fill(HIST("hSelectionCounter"), 3);
    if (!track.hasITS()) return false;
    registry.fill(HIST("hSelectionCounter"), 4);
    if (std::abs(eta(track.px(), track.py(), track.pz())) > PcEtaCut) return false;
    registry.fill(HIST("hSelectionCounter"), 5);
    if (std::abs(track.pt()) > tracksPtMaxCut) return false;
    registry.fill(HIST("hSelectionCounter"), 6);
    return true;
  }

  template <typename T>
  bool tracksPassPiPID(const T& cutTracks) { // check if pre-cut tracks pass PID cut (fall within n-dimensional hypersphere)
    double radius = 0.0;
    for (const auto& track : cutTracks) radius += std::pow(track.tpcNSigmaPi(), 2);
    return radius < std::pow(tracksTpcNSigmaPiCut, 2);
  }

  template <typename T>
  double tracksTotalCharge(const T& cutTracks) { // calculate total charge of selected tracks
    double charge = 0.0;
    for (const auto& track : cutTracks) charge += track.sign();
    return charge;
  }

  ROOT::Math::PxPyPzMVector reconstructSystem(const std::vector<ROOT::Math::PxPyPzMVector>& cutTracks4Vecs) { // reconstruct system from 4-vectors
    ROOT::Math::PxPyPzMVector system;
    for (const auto& track4Vec : cutTracks4Vecs) system += track4Vec;
    return system;
  }

  void processReco(FullUDCollision const& collision, FullUDTracks const& tracks, aod::UDZdcsReduced const&) {
    registry.fill(HIST("QC/collisions/uncut/hNetCharge"), collision.netCharge());
    registry.fill(HIST("QC/collisions/uncut/hNumContributors"), collision.numContrib());
    registry.fill(HIST("QC/collisions/uncut/hRgtrwTOF"), collision.rgtrwTOF());
    registry.fill(HIST("QC/collisions/uncut/hFt0Amplitude"), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC());
    registry.fill(HIST("QC/collisions/uncut/hFt0Time"), collision.timeFT0A(), collision.timeFT0C());
    registry.fill(HIST("QC/collisions/uncut/hFddAmplitude"), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC());
    registry.fill(HIST("QC/collisions/uncut/hFddTime"), collision.timeFDDA(), collision.timeFDDC());
    registry.fill(HIST("QC/collisions/uncut/hFv0Amplitude"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("QC/collisions/uncut/hFv0Time"), collision.timeFV0A());
    registry.fill(HIST("QC/collisions/uncut/hPosX"), collision.posX());
    registry.fill(HIST("QC/collisions/uncut/hPosY"), collision.posY());
    registry.fill(HIST("QC/collisions/uncut/hPosZ"), collision.posZ());
    registry.fill(HIST("QC/collisions/uncut/hZDCcommonEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());
    registry.fill(HIST("QC/collisions/uncut/hZDCtime"), collision.timeZNA(), collision.timeZNC());

    // apply some cuts on collisions
    if (std::abs(collision.posZ()) > collisionsPosZMaxCut) return;
    // if (collision.has_<aod::UDZdcsReduced>()) return; // this doesn't work for some reason

    registry.fill(HIST("QC/collisions/cut/hNetCharge"), collision.netCharge());
    registry.fill(HIST("QC/collisions/cut/hNumContributors"), collision.numContrib());
    registry.fill(HIST("QC/collisions/cut/hRgtrwTOF"), collision.rgtrwTOF());
    registry.fill(HIST("QC/collisions/cut/hFt0Amplitude"), collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC());
    registry.fill(HIST("QC/collisions/cut/hFt0Time"), collision.timeFT0A(), collision.timeFT0C());
    registry.fill(HIST("QC/collisions/cut/hFddAmplitude"), collision.totalFDDAmplitudeA(), collision.totalFDDAmplitudeC());
    registry.fill(HIST("QC/collisions/cut/hFddTime"), collision.timeFDDA(), collision.timeFDDC());
    registry.fill(HIST("QC/collisions/cut/hFv0Amplitude"), collision.totalFV0AmplitudeA());
    registry.fill(HIST("QC/collisions/cut/hFv0Time"), collision.timeFV0A());
    registry.fill(HIST("QC/collisions/cut/hPosX"), collision.posX());
    registry.fill(HIST("QC/collisions/cut/hPosY"), collision.posY());
    registry.fill(HIST("QC/collisions/cut/hPosZ"), collision.posZ());

    // vectors for storing selected tracks and their 4-vectors
    std::vector<decltype(tracks.begin())> cutTracks;
    std::vector<ROOT::Math::PxPyPzMVector> cutTracks4Vecs;

    for (const auto& track : tracks) {
      registry.fill(HIST("QC/allTracks/hPt"), track.pt());
      registry.fill(HIST("QC/allTracks/hEta"), eta(track.px(), track.py(), track.pz()));
      registry.fill(HIST("QC/allTracks/hPhi"), phi(track.px(), track.py()));
      registry.fill(HIST("QC/allTracks/hCharge"), track.sign());
      registry.fill(HIST("QC/allTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/allTracks/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/allTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      registry.fill(HIST("QC/allTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("QC/allTracks/hTpcNClsFindable"), track.tpcNClsFindable());
      registry.fill(HIST("QC/allTracks/hTpcNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());

      if (!trackPassesCuts(track)) continue; // apply cuts
      cutTracks.push_back(track);
      cutTracks4Vecs.push_back(ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), o2::constants::physics::MassPionCharged)); // apriori assume pion mass
      
      registry.fill(HIST("QC/cutTracks/hPt"), track.pt());
      registry.fill(HIST("QC/cutTracks/hEta"), eta(track.px(), track.py(), track.pz()));
      registry.fill(HIST("QC/cutTracks/hPhi"), phi(track.px(), track.py()));
      registry.fill(HIST("QC/cutTracks/hCharge"), track.sign());
      registry.fill(HIST("QC/cutTracks/hTpcNSigmaPi"), track.tpcNSigmaPi());
      registry.fill(HIST("QC/cutTracks/hTofNSigmaPi"), track.tofNSigmaPi());
      registry.fill(HIST("QC/cutTracks/hDcaXYZ"), track.dcaZ(), track.dcaXY());
      registry.fill(HIST("QC/cutTracks/hTpcSignalVsPt"), track.pt(), track.tpcSignal());
      registry.fill(HIST("QC/cutTracks/hTpcNClsFindable"), track.tpcNClsFindable());
      registry.fill(HIST("QC/cutTracks/hTpcNClsFindableMinusFound"), track.tpcNClsFindableMinusFound());
    }
    registry.fill(HIST("QC/allTracks/hNTracks"), tracks.size());
    registry.fill(HIST("QC/cutTracks/hNTracks"), cutTracks.size());

    auto system = reconstructSystem(cutTracks4Vecs); // reconstruct system from 4-vectors

    if (cutTracks.size() == 2) {
      ROOT::Math::PxPyPzMVector piPos, piNeg;
      decltype(tracks.begin()) posTrack, negTrack;
      
      registry.fill(HIST("hSelectionCounter"), 7, 2.); // the weight is 2. because we are looking at pairs of tracks

      // unlike-sign tracks
      if (tracksTotalCharge(cutTracks) == 0) {
        registry.fill(HIST("hSelectionCounter"), 8, 2.);
        if (cutTracks[0].sign() > 0) {
          piPos = cutTracks4Vecs[0];
          posTrack = cutTracks[0];
          piNeg = cutTracks4Vecs[1];
          negTrack = cutTracks[1];
        } else {
          piNeg = cutTracks4Vecs[0];
          negTrack = cutTracks[0];
          piPos = cutTracks4Vecs[1];
          posTrack = cutTracks[1];
        }

        registry.fill(HIST("QC/cutTracks/hTpcNsigmaPi2D"), posTrack.tpcNSigmaPi(), negTrack.tpcNSigmaPi());
        // pi pairs
        if (tracksPassPiPID(cutTracks)) { // PID performed here because I also need to check for electron-positron pairs
          registry.fill(HIST("reco/pions/unlike-sign/hPt"), piPos.Pt(), piNeg.Pt());
          registry.fill(HIST("reco/pions/unlike-sign/hEta"), piPos.Eta(), piNeg.Eta());
          registry.fill(HIST("reco/pions/unlike-sign/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI); // shift by pi to get to the range [0, 2pi]
          registry.fill(HIST("reco/system/2pi/unlike-sign/hM"), system.M());
          registry.fill(HIST("reco/system/2pi/unlike-sign/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2pi/unlike-sign/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2pi/unlike-sign/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2pi/unlike-sign/hPhi"), system.Phi() + o2::constants::math::PI);
        }
        // electron-positron pairs
        if (std::pow(posTrack.tpcNSigmaEl(), 2) + std::pow(negTrack.tpcNSigmaEl(), 2) < std::pow(tracksTpcNSigmaElCut, 2)) {
          registry.fill(HIST("reco/electrons/hPt"), piPos.Pt(), piNeg.Pt());
          registry.fill(HIST("reco/electrons/hEta"), piPos.Eta(), piNeg.Eta());
          registry.fill(HIST("reco/electrons/hPhi"), piPos.Phi() + o2::constants::math::PI, piNeg.Phi() + o2::constants::math::PI);
          registry.fill(HIST("reco/system/2el/hM"), system.M());
          registry.fill(HIST("reco/system/2el/hPt"), system.Pt());
          registry.fill(HIST("reco/system/2el/hPtVsM"), system.M(), system.Pt());
          registry.fill(HIST("reco/system/2el/hY"), system.Rapidity());
          registry.fill(HIST("reco/system/2el/hPhi"), system.Phi() + o2::constants::math::PI);
        }
      }

      // like-sign tracks PIDed as pions
      if (tracksTotalCharge(cutTracks) != 0 && tracksPassPiPID(cutTracks)) {
        registry.fill(HIST("reco/pions/like-sign/hPt"), cutTracks4Vecs[0].Pt(), cutTracks4Vecs[1].Pt());
        registry.fill(HIST("reco/pions/like-sign/hEta"), cutTracks4Vecs[0].Eta(), cutTracks4Vecs[1].Eta());
        registry.fill(HIST("reco/pions/like-sign/hPhi"), cutTracks4Vecs[0].Phi() + o2::constants::math::PI, cutTracks4Vecs[1].Phi() + o2::constants::math::PI);
        registry.fill(HIST("reco/system/2pi/like-sign/hM"), system.M());
        registry.fill(HIST("reco/system/2pi/like-sign/hPt"), system.Pt());
        registry.fill(HIST("reco/system/2pi/like-sign/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/2pi/like-sign/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/2pi/like-sign/hPhi"), system.Phi() + o2::constants::math::PI);
      }
    }

    if (cutTracks.size() == 4 && tracksPassPiPID(cutTracks)) {
      if (tracksTotalCharge(cutTracks) == 0) {
        registry.fill(HIST("reco/system/4pi/net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/4pi/net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/4pi/net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/4pi/net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/4pi/net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      } else {
        registry.fill(HIST("reco/system/4pi/non-net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/4pi/non-net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      }
    }

    if (cutTracks.size() == 6 && tracksPassPiPID(cutTracks)) {
      if (tracksTotalCharge(cutTracks) == 0) {
        registry.fill(HIST("reco/system/6pi/net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/6pi/net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/6pi/net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/6pi/net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/6pi/net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      } else {
        registry.fill(HIST("reco/system/6pi/non-net-zero/hM"), system.M());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPt"), system.Pt());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPtVsM"), system.M(), system.Pt());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hY"), system.Rapidity());
        registry.fill(HIST("reco/system/6pi/non-net-zero/hPhi"), system.Phi() + o2::constants::math::PI);
      }
    }
  } PROCESS_SWITCH(upcRhoAnalysis, processReco, "analyse reco tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<upcRhoAnalysis>(cfgc)
  };
}