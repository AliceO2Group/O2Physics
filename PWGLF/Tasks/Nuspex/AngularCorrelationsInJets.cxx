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
// Authors: Lars Jörgensen
// Date: 29.05.2024

#include <cmath>
#include <vector>
#include <TMath.h>
#include <TObjArray.h>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "TVector2.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCLfFullPr, aod::pidTPCLfFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::TOFSignal, aod::pidTOFmass, aod::pidTOFbeta>;

struct AngularCorrelationsInJets {

  HistogramRegistry registryData{"jetOutput", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(InitContext&)
  {
    // AxisSpec specName = {binningInfo, "axisLabel"};
    AxisSpec ptAxis = {1000, 0, 100, "#it{p}_{T} [GeV/#it{c}]"};
    AxisSpec particleTypeAxis = {4, 1, 5, "[p, ap, d, ad]"};
    AxisSpec nsigmapTAxis = {1000, -50, 50, "#it{p}_{T} [GeV/#it{c}]"};
    AxisSpec nsigmaAxis = {1000, -5, 5, "n#sigma"};
    AxisSpec dcazAxis = {200, -3, 3, "DCA_{z} [cm]"};
    AxisSpec dcaxyAxis = {200, -2, 2, "DCA_{xy} [cm]"};
    AxisSpec angDistPhiAxis = {1000, -2, 5, "#Delta#varphi"};
    AxisSpec angDistEtaAxis = {1000, -2, 2, "#Delta#eta"};

    // Counters
    registryData.add("hNumberOfEvents", "Number of events", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hNumberOfJets", "Total number of jets", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hEventProtocol", "Event protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hTrackProtocol", "Track protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hNumPartInJet", "Number of particles in a jet", HistType::kTH1I, {{200, 0, 200}});

    // (Pseudo)Rapidity
    registryData.add("hEtaFullEvent", "Particle pseudorapidity;#eta", HistType::kTH1F, {{200, -1, 1}});
    registryData.add("hJetRapidity", "Jet rapidity;#it{y}", HistType::kTH1F, {{200, -1, 1}});

    // pT
    registryData.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {ptAxis});
    registryData.add("hPtJetParticle", "p_{T} of particles in jets", HistType::kTH1D, {ptAxis});
    registryData.add("hPtSubtractedJet", "Subtracted jet p_{T}", HistType::kTH1D, {ptAxis});
    registryData.add("hPtJetProtonDeuteron", "p_{T} of (anti)p, (anti)d", HistType::kTH2D, {particleTypeAxis, ptAxis});
    registryData.add("hPtTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{2000, 0, 500}});
    registryData.add("hPtDiff", "pT difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{100, -0.0000005, 0.0000005}});

    // nSigma
    registryData.add("hTPCnsigma", "TPC n#sigma for full event", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
    registryData.add("hTOFnsigma", "TOF n#sigma for full event", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
    registryData.add("hTPCnsigmaProton", "TPC n#sigma for (anti)proton", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
    registryData.add("hTOFnsigmaProton", "TOF n#sigma for (anti)proton", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
    registryData.add("hTPCnsigmaDeuteron", "TPC n#sigma for (anti)deuteron", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});
    registryData.add("hTOFnsigmaDeuteron", "TOF n#sigma for (anti)deuteron", HistType::kTH2F, {nsigmapTAxis, nsigmaAxis});

    // DCA
    registryData.add("hDCAxyFullJet", "DCA_{xy} of full jet", HistType::kTH2F, {ptAxis, dcaxyAxis});
    registryData.add("hDCAzFullJet", "DCA_{z} of full jet", HistType::kTH2F, {ptAxis, dcazAxis});
    registryData.add("hDCAzJetProton", "DCA_{z} of protons after TPC cut", HistType::kTH2F, {ptAxis, dcazAxis});
    registryData.add("hDCAzJetAntiproton", "DCA_{z} of antiprotons after TPC cut", HistType::kTH2F, {ptAxis, dcazAxis});
    registryData.add("hDCAzJetDeuteron", "DCA_{z} of deuterons after TPC cut", HistType::kTH2F, {ptAxis, dcazAxis});
    registryData.add("hDCAzJetAntideuteron", "DCA_{z} of antideuterons after TPC cut", HistType::kTH2F, {ptAxis, dcazAxis});

    // Angular Distributions
    registryData.add("hDeltaPhiSE", "#Delta#varphi of (anti)p, (anti)d in single event", HistType::kTH2D, {particleTypeAxis, angDistPhiAxis});
    registryData.add("hDeltaPhiME", "#Delta#varphi of (anti)p, (anti)d in mixed events", HistType::kTH2D, {particleTypeAxis, angDistPhiAxis});
    registryData.add("hDeltaPhiEtaSE", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in single event", HistType::kTH3D, {particleTypeAxis, angDistPhiAxis, angDistEtaAxis});
    registryData.add("hDeltaPhiEtaME", "#Delta#varphi vs #Delta#eta of (anti)p, (anti)d in mixed events", HistType::kTH3D, {particleTypeAxis, angDistPhiAxis, angDistEtaAxis});

    registryData.add("hJetConeRadius", "Jet Radius;#it{R}", HistType::kTH1F, {{100, 0, 1}});
  }

  // Preliminary Cuts
  // Configurable<int> fTPCRefit{"TPCRefit", 0, "Require TPC refit"};
  Configurable<float> fMinNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> fMinReqClusterITS{"minReqClusterITS", 2.0, "min number of clusters required in ITS"};
  Configurable<float> fMinRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.7f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> fMaxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> fMaxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> fMaxDCAxy{"maxDCA_xy", 0.5f, "max DCA to vertex xy"};
  Configurable<float> fMaxDCAz{"maxDCA_z", 2.4f, "max DCA to vertex z"};
  Configurable<float> fMaxEta{"maxEta", 0.8, "max pseudorapidity"};

  // Jet Cuts
  Configurable<float> fJetRadius{"jetRadius", 0.4, "jet radius R"};
  Configurable<float> fMinJetPt{"minJetPt", 10.0, "minimum total pT to accept jet"};
  Configurable<float> fMinJetParticlePt{"minJetParticlePt", 0.0, "minimum pT to accept jet particle"};
  Configurable<float> fMinLeadingPt{"minLeadingPt", 5.0, "minimum pT for leading track"};

  // Proton Cuts
  Configurable<float> fProtonDCAxy{"protonDCAxy", 0.5, "[proton] DCAxy cut"};
  Configurable<float> fProtonDCAz{"protonDCAz", 1.0, "[proton] DCAz cut"};
  Configurable<float> fProtonTPCTOFpT{"protonTPCTOFswitchpT", 0.7, "[proton] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fProtonTPCnsigLow{"protonTPCnsigmapTLow", 5.0, "[proton] max TPC nsigma with low pT"};
  Configurable<float> fProtonTPCnsigHigh{"protonTPCnsigmapTHigh", 4.0, "[proton] max TPC nsigma with high pT"};
  Configurable<float> fProtonTOFnsigLow{"protonTOFnsigmapTLow", 10.0, "[proton] max TOF nsigma with low pT"};
  Configurable<float> fProtonTOFnsigHigh{"protonTOFnsigmapTHigh", 10.0, "[proton] max TOF nsigma with high pT"};

  // Antiproton Cuts
  Configurable<float> fAntiprotonDCAxy{"antiprotonDCAxy", 0.5, "[antiproton] DCAxy cut"};
  Configurable<float> fAntiprotonDCAz{"antiprotonDCAz", 1.0, "[antiproton] DCAz cut"};
  Configurable<float> fAntiprotonTPCTOFpT{"antiprotonTPCTOFswitchpT", 0.7, "[antiproton] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fAntiprotonTPCnsigLow{"antiprotonTPCnsigmapTLow", 5.0, "[antiproton] max TPC nsigma with low pT"};
  Configurable<float> fAntiprotonTPCnsigHigh{"antiprotonTPCnsigmapTHigh", 4.0, "[antiproton] max TPC nsigma with high pT"};
  Configurable<float> fAntiprotonTOFnsigLow{"antiprotonTOFnsigmapTLow", 10.0, "[antiproton] max TOF nsigma with low pT"};
  Configurable<float> fAntiprotonTOFnsigHigh{"antiprotonTOFnsigmapTHigh", 10.0, "[antiproton] max TOF nsigma with high pT"};

  // Deuteron Cuts
  Configurable<float> fDeuteronDCAxy{"deuteronDCAxy", 0.5, "[deuteron] DCAxy cut"};
  Configurable<float> fDeuteronDCAz{"deuteronDCAz", 1.0, "[deuteron] DCAz cut"};
  Configurable<float> fDeuteronTPCTOFpT{"deuteronTPCTOFswitchpT", 0.7, "[deuteron] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fDeuteronTPCnsigLow{"deuteronTPCnsigmapTLow", 5.0, "[deuteron] max TPC nsigma with low pT"};
  Configurable<float> fDeuteronTPCnsigHigh{"deuteronTPCnsigmapTHigh", 4.0, "[deuteron] max TPC nsigma with high pT"};
  Configurable<float> fDeuteronTOFnsigLow{"deuteronTOFnsigmapTLow", 10.0, "[deuteron] max TOF nsigma with low pT"};
  Configurable<float> fDeuteronTOFnsigHigh{"deuteronTOFnsigmapTHigh", 10.0, "[deuteron] max TOF nsigma with high pT"};

  // Antideuteron Cuts
  Configurable<float> fAntideuteronDCAxy{"antideuteronDCAxy", 0.5, "[antideuteron] DCAxy cut"};
  Configurable<float> fAntideuteronDCAz{"antideuteronDCAz", 1.0, "[antideuteron] DCAz cut"};
  Configurable<float> fAntideuteronTPCTOFpT{"antideuteronTPCTOFswitchpT", 0.7, "[antideuteron] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fAntideuteronTPCnsigLow{"antideuteronTPCnsigmapTLow", 5.0, "[antideuteron] max TPC nsigma with low pT"};
  Configurable<float> fAntideuteronTPCnsigHigh{"antideuteronTPCnsigmapTHigh", 4.0, "[antideuteron] max TPC nsigma with high pT"};
  Configurable<float> fAntideuteronTOFnsigLow{"antideuteronTOFnsigmapTLow", 10.0, "[antideuteron] max TOF nsigma with low pT"};
  Configurable<float> fAntideuteronTOFnsigHigh{"antideuteronTOFnsigmapTHigh", 10.0, "[antideuteron] max TOF nsigma with high pT"};

  Configurable<int> fTrackBufferSize{"trackBufferSize", 2000, "Number of mixed-event tracks being stored"};

  std::vector<typename FullTracks::iterator> fTrackBufferProton;
  std::vector<typename FullTracks::iterator> fTrackBufferAntiproton;
  std::vector<typename FullTracks::iterator> fTrackBufferDeuteron;
  std::vector<typename FullTracks::iterator> fTrackBufferAntideuteron;

  //****************************************************************************************************

  template <typename CollisionType, typename TracksType>
  void angCorrData(const CollisionType& /*collision*/, const TracksType& tracks)
  {
    registryData.fill(HIST("hNumberOfEvents"), 0);
    registryData.fill(HIST("hEventProtocol"), 0);

    std::vector<fastjet::PseudoJet> jetInput;
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> constituents;
    fastjet::PseudoJet hardestJet(0., 0., 0., 0.);
    fastjet::PseudoJet subtractedJet(0., 0., 0., 0.);
    jetInput.clear();
    jets.clear();
    constituents.clear();

    for (const auto& track : tracks) {
      registryData.fill(HIST("hTrackProtocol"), 0);
      registryData.fill(HIST("hEtaFullEvent"), track.eta());
      registryData.fill(HIST("hPtFullEvent"), track.pt());

      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(track.mass())); // is this fine? just TOF for invariant mass?
      inputPseudoJet.set_user_index(track.globalIndex());
      jetInput.emplace_back(inputPseudoJet);
    }

    if (jetInput.size() < 2)
      return;
    registryData.fill(HIST("hEventProtocol"), 1);

    double ghost_maxrap = 1.0;
    double ghost_area = 0.005;
    int ghost_repeat = 1;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, fJetRadius);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap, ghost_repeat, ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef); // or CSActiveArea?
    jets = sorted_by_pt(clusterSeq.inclusive_jets());
    if (jets.size() == 0)
      return;
    registryData.fill(HIST("hEventProtocol"), 2);

    hardestJet = jets[0];

    if (hardestJet.pt() < fMinJetPt)
      return;
    registryData.fill(HIST("hEventProtocol"), 3);
    if (hardestJet.constituents().size() < 2)
      return;
    registryData.fill(HIST("hEventProtocol"), 4);
    registryData.fill(HIST("hNumberOfJets"), 0);
    registryData.fill(HIST("hPtTotalJet"), hardestJet.pt());
    registryData.fill(HIST("hJetRapidity"), hardestJet.rap());

    constituents = hardestJet.constituents();

    for (int i = 0; i < constituents.size(); i++) {
      registryData.fill(HIST("hPtJetParticle"), constituents[i].pt());
      double DeltaPhi = TVector2::Phi_0_2pi(constituents[i].phi() - hardestJet.phi());
      double DeltaEta = constituents[i].eta() - hardestJet.eta();
      double Delta = TMath::Sqrt(DeltaPhi * DeltaPhi + DeltaEta * DeltaEta) / (constituents[i].pt() * constituents[i].pt()); // need 1/pT^2?
      registryData.fill(HIST("hJetConeRadius"), Delta);
    }

    fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2)); // fix!
    fastjet::JetMedianBackgroundEstimator bkgEst(selector, jetDefBkg, areaDefBkg);
    fastjet::Subtractor subtractor(&bkgEst);
    subtractor.set_use_rho_m(true);
    bkgEst.set_particles(jetInput);

    subtractedJet = subtractor(hardestJet);
    if (subtractedJet.has_constituents()) {
      for (int i = 0; i < subtractedJet.constituents().size(); i++) {
        registryData.fill(HIST("hPtSubtractedJet"), subtractedJet.constituents()[i].pt());
      }
    }

    std::vector<typename TracksType::iterator> jetProtons;
    std::vector<typename TracksType::iterator> jetAntiprotons;
    std::vector<typename TracksType::iterator> jetDeuterons;
    std::vector<typename TracksType::iterator> jetAntideuterons;

    for (int i = 0; i < constituents.size(); i++) {
      fastjet::PseudoJet pseudoParticle = constituents[i];
      int id = pseudoParticle.user_index();
      typename TracksType::iterator jetParticle = tracks.iteratorAt(id);
      registryData.fill(HIST("hDCAxyFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("hDCAzFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("hTPCnsigma"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
      registryData.fill(HIST("hTOFnsigma"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());

      double ptDiff = pseudoParticle.pt() - jetParticle.pt();
      registryData.fill(HIST("hPtDiff"), ptDiff);
      int particleType = 0;

      if (jetParticle.pt() < fMinJetParticlePt)
        continue;
      if (IsProton(jetParticle) || IsAntiproton(jetParticle)) { // collect (anti)protons in jet
        registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
        registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());
        if (IsProton(jetParticle)) {
          particleType = 1;
          registryData.fill(HIST("hTrackProtocol"), 6); // # protons
          jetProtons.emplace_back(jetParticle);
        } else {
          particleType = 2;
          registryData.fill(HIST("hTrackProtocol"), 7); // # antiprotons
          jetAntiprotons.emplace_back(jetParticle);
        }
      } else if (IsDeuteron(jetParticle) || IsAntideuteron(jetParticle)) { // collect (anti)deuterons in jet
        registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaDe());
        registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaDe());
        if (IsDeuteron(jetParticle)) {
          particleType = 3;
          registryData.fill(HIST("hTrackProtocol"), 8); // # deuterons
          jetDeuterons.emplace_back(jetParticle);
        } else {
          particleType = 4;
          registryData.fill(HIST("hTrackProtocol"), 9); // # antideuterons
          jetAntideuterons.emplace_back(jetParticle);
        }
      }
      registryData.fill(HIST("hPtJetProtonDeuteron"), particleType, jetParticle.pt());
    } // for (int i=0; i<(int)constituents.size(); i++)

    if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2) && (jetDeuterons.size() < 2) && (jetAntideuterons.size() < 2))
      return;
    registryData.fill(HIST("hEventProtocol"), 5);

    if (jetProtons.size() > 1) {
      for (int i = 0; i < jetProtons.size(); i++) {
        for (int j = i + 1; j < jetProtons.size(); j++) {
          double DeltaPhi = TVector2::Phi_0_2pi(jetProtons[i].phi() - jetProtons[j].phi());
          double DeltaEta = TMath::Abs(jetProtons[i].eta() - jetProtons[j].eta());
          if (DeltaPhi > (1.5 * TMath::Pi())) {
            DeltaPhi = DeltaPhi - 2 * TMath::Pi();
          }
          registryData.fill(HIST("hDeltaPhiSE"), 1, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaSE"), 1, DeltaPhi, DeltaEta);
        }
        FillMixedEventDeltas(jetProtons[i], 1);
        SetTrackBuffer(jetProtons[i], 1);
      }
    }
    if (jetAntiprotons.size() > 1) {
      for (int i = 0; i < jetAntiprotons.size(); i++) {
        for (int j = i + 1; j < jetAntiprotons.size(); j++) {
          double DeltaPhi = TVector2::Phi_0_2pi(jetAntiprotons[i].phi() - jetAntiprotons[j].phi());
          double DeltaEta = TMath::Abs(jetAntiprotons[i].eta() - jetAntiprotons[j].eta());
          if (DeltaPhi > (1.5 * TMath::Pi())) {
            DeltaPhi = DeltaPhi - 2 * TMath::Pi();
          }
          registryData.fill(HIST("hDeltaPhiSE"), 2, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaSE"), 2, DeltaPhi, DeltaEta);
        }
        FillMixedEventDeltas(jetAntiprotons[i], 2);
        SetTrackBuffer(jetAntiprotons[i], 2);
      }
    }
    if (jetDeuterons.size() > 1) {
      for (int i = 0; i < jetDeuterons.size(); i++) {
        for (int j = i + 1; j < jetDeuterons.size(); j++) {
          double DeltaPhi = TVector2::Phi_0_2pi(jetDeuterons[i].phi() - jetDeuterons[j].phi());
          double DeltaEta = TMath::Abs(jetDeuterons[i].eta() - jetDeuterons[j].eta());
          if (DeltaPhi > (1.5 * TMath::Pi())) {
            DeltaPhi = DeltaPhi - 2 * TMath::Pi();
          }
          registryData.fill(HIST("hDeltaPhiSE"), 3, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaSE"), 3, DeltaPhi, DeltaEta);
        }
        FillMixedEventDeltas(jetDeuterons[i], 3);
        SetTrackBuffer(jetDeuterons[i], 3);
      }
    }
    if (jetAntideuterons.size() > 1) {
      for (int i = 0; i < jetAntideuterons.size(); i++) {
        for (int j = i + 1; j < jetAntideuterons.size(); j++) {
          double DeltaPhi = TVector2::Phi_0_2pi(jetAntideuterons[i].phi() - jetAntideuterons[j].phi());
          double DeltaEta = TMath::Abs(jetAntideuterons[i].eta() - jetAntideuterons[j].eta());
          if (DeltaPhi > (1.5 * TMath::Pi())) {
            DeltaPhi = DeltaPhi - 2 * TMath::Pi();
          }
          registryData.fill(HIST("hDeltaPhiSE"), 4, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaSE"), 4, DeltaPhi, DeltaEta);
        }
        FillMixedEventDeltas(jetAntideuterons[i], 4);
        SetTrackBuffer(jetAntideuterons[i], 4);
      }
    }
  }

  template <typename T1>
  bool IsProton(const T1& track)
  {
    bool isProton = false;
    if (track.sign() < 0)
      return isProton;

    if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fProtonTPCnsigLow)
      isProton = true;
    if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fProtonTPCnsigHigh)
      isProton = true;

    registryData.fill(HIST("hDCAzJetProton"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fProtonDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fProtonDCAz)
      return false;

    if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fProtonTOFnsigLow)
      isProton = true;
    if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fProtonTOFnsigHigh)
      isProton = true;

    return isProton;
  }

  template <typename T2>
  bool IsAntiproton(const T2& track)
  {
    bool isAntiproton = false;
    if (track.sign() < 0)
      return isAntiproton;

    if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fAntiprotonTPCnsigLow)
      isAntiproton = true;
    if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fAntiprotonTPCnsigHigh)
      isAntiproton = true;

    registryData.fill(HIST("hDCAzJetAntiproton"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fAntiprotonDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fAntiprotonDCAz)
      return false;

    if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fAntiprotonTOFnsigLow)
      isAntiproton = true;
    if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fAntiprotonTOFnsigHigh)
      isAntiproton = true;

    return isAntiproton;
  }

  template <typename T3>
  bool IsDeuteron(const T3& track)
  {
    bool isDeuteron = false;
    if (track.sign() < 0)
      return isDeuteron;

    if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fDeuteronTPCnsigLow)
      isDeuteron = true;
    if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fDeuteronTPCnsigHigh)
      isDeuteron = true;

    registryData.fill(HIST("hDCAzJetDeuteron"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fDeuteronDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fDeuteronDCAz)
      return false;

    if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fDeuteronTOFnsigLow)
      isDeuteron = true;
    if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fDeuteronTOFnsigHigh)
      isDeuteron = true;

    return isDeuteron;
  }

  template <typename T4>
  bool IsAntideuteron(const T4& track)
  {
    bool isAntideuteron = false;
    if (track.sign() < 0)
      return isAntideuteron;

    if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fAntideuteronTPCnsigLow)
      isAntideuteron = true;
    if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fAntideuteronTPCnsigHigh)
      isAntideuteron = true;

    registryData.fill(HIST("hDCAzJetAntideuteron"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fAntideuteronDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fAntideuteronDCAz)
      return false;

    if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fAntideuteronTOFnsigLow)
      isAntideuteron = true;
    if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fAntideuteronTOFnsigHigh)
      isAntideuteron = true;

    return isAntideuteron;
  }

  template <typename T5>
  void SetTrackBuffer(const T5& track, int particleType)
  {
    switch (particleType) {
      case 1:
        if (fTrackBufferProton.size() == fTrackBufferSize) {
          fTrackBufferProton.insert(fTrackBufferProton.begin(), track);
          fTrackBufferProton.resize(fTrackBufferSize);
        } else if (fTrackBufferProton.size() < fTrackBufferSize) {
          fTrackBufferProton.emplace_back(track);
        }
        break;
      case 2:
        if (fTrackBufferAntiproton.size() == fTrackBufferSize) {
          fTrackBufferAntiproton.insert(fTrackBufferAntiproton.begin(), track);
          fTrackBufferAntiproton.resize(fTrackBufferSize);
        } else if (fTrackBufferAntiproton.size() < fTrackBufferSize) {
          fTrackBufferAntiproton.emplace_back(track);
        }
        break;
      case 3:
        if (fTrackBufferDeuteron.size() == fTrackBufferSize) {
          fTrackBufferDeuteron.insert(fTrackBufferDeuteron.begin(), track);
          fTrackBufferDeuteron.resize(fTrackBufferSize);
        } else if (fTrackBufferDeuteron.size() < fTrackBufferSize) {
          fTrackBufferDeuteron.emplace_back(track);
        }
        break;
      case 4:
        if (fTrackBufferAntideuteron.size() == fTrackBufferSize) {
          fTrackBufferAntideuteron.insert(fTrackBufferAntideuteron.begin(), track);
          fTrackBufferAntideuteron.resize(fTrackBufferSize);
        } else if (fTrackBufferAntideuteron.size() < fTrackBufferSize) {
          fTrackBufferAntideuteron.emplace_back(track);
        }
        break;
      default:
        LOG(warn) << "SetTrackBuffer: invalid particle ID!";
    }
  }

  template <typename T6>
  void FillMixedEventDeltas(const T6& track, int particleType)
  {
    switch (particleType) {
      case 1:
        if (fTrackBufferProton.size() == 0)
          return;
        for (int i = 0; i < fTrackBufferProton.size(); i++) { // can I do this even if the track buffer isn't even full yet?
          double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferProton[i].phi());
          double DeltaEta = TMath::Abs(track.eta() - fTrackBufferProton[i].eta());
          registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
        }
        break;
      case 2:
        if (fTrackBufferAntiproton.size() == 0)
          return;
        for (int i = 0; i < fTrackBufferAntiproton.size(); i++) {
          double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferAntiproton[i].phi());
          double DeltaEta = TMath::Abs(track.eta() - fTrackBufferAntiproton[i].eta());
          registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
        }
        break;
      case 3:
        if (fTrackBufferDeuteron.size() == 0)
          return;
        for (int i = 0; i < fTrackBufferDeuteron.size(); i++) {
          double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferDeuteron[i].phi());
          double DeltaEta = TMath::Abs(track.eta() - fTrackBufferDeuteron[i].eta());
          registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
        }
        break;
      case 4:
        if (fTrackBufferAntideuteron.size() == 0)
          return;
        for (int i = 0; i < fTrackBufferAntideuteron.size(); i++) {
          double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - fTrackBufferAntideuteron[i].phi());
          double DeltaEta = TMath::Abs(track.eta() - fTrackBufferAntideuteron[i].eta());
          registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
        }
        break;
      default:
        LOG(warn) << "FillMixedEventDeltas: invalid particle ID!";
    }
  }

  //****************************************************************************************************

  Filter prelimTrackCuts = (/* aod::track::TPCrefit == fTPCRefit && */ aod::track::itsChi2NCl < fMaxChi2ITS && aod::track::tpcChi2NCl < fMaxChi2TPC && nabs(aod::track::dcaXY) < fMaxDCAxy && nabs(aod::track::dcaZ) < fMaxDCAz && nabs(aod::track::eta) < fMaxEta);

  void process_ang_corr_data(aod::Collision const& collision, soa::Filtered<FullTracks> const& tracks)
  {
    angCorrData(collision, tracks);
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, process_ang_corr_data, "ang correlations in reco'ed jets", true);
};

//****************************************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc, TaskName{"angular-correlations-in-jets"})};
}
