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
// author: Lars Jörgensen

#include <vector>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/BasicCCDBManager.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/TableProducer/PID/pidTOFBase.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "TVector2.h"
#include "TVector3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AxisSpecs {
  AxisSpec ptAxisPos = {1000, 0, 100, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec ptAxisFull = {2000, -100, 100, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec nsigmapTAxis = {1000, -50, 50, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec nsigmaAxis = {1000, -15, 15, "n#sigma"};
  AxisSpec dcazAxis = {1000, -1, 1, "DCA_{z} [cm]"};
  AxisSpec dcaxyAxis = {1000, -0.5, 0.5, "DCA_{xy} [cm]"};
  AxisSpec angDistPhiAxis = {1000, -2, 5, "#Delta#varphi"};
  AxisSpec angDistEtaAxis = {1000, -2, 2, "#Delta#eta"};
};

struct AngularCorrelationsInJets {
  // Preliminary Cuts
  Configurable<int> fMinNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of crossed rows TPC"};
  Configurable<int> fMinReqClusterITS{"minReqClusterITS", 2, "min number of clusters required in ITS"};
  Configurable<int> fMinReqClusterTPC{"minReqClusterTPC", 70, "min number of clusters required in TPC"};
  Configurable<float> fMinRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.7, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> fMaxChi2ITS{"maxChi2ITS", 36.0, "max chi2 per cluster ITS"};
  Configurable<float> fMaxChi2TPC{"maxChi2TPC", 4.0, "max chi2 per cluster TPC"};
  Configurable<float> fMaxDCAxy{"maxDCA_xy", 0.5, "max DCA to vertex xy"};
  Configurable<float> fMaxDCAz{"maxDCA_z", 1.0, "max DCA to vertex z"};
  Configurable<float> fMaxEta{"maxEta", 0.8, "max pseudorapidity"}; // consider jet cone?

  // Jet Cuts
  Configurable<float> fJetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> fMinJetPt{"minJetPt", 10.0, "minimum total pT to accept jet"};
  Configurable<float> fMinJetParticlePt{"minJetParticlePt", 0.0, "minimum pT to accept jet particle"};
  Configurable<float> fMinLeadingPt{"minLeadingPt", 5.0, "minimum pT for leading track"};

  // Proton Cuts
  Configurable<float> fProtonDCAxyYield{"protonDCAxyYield", 0.05, "[proton] DCAxy cut for yield"};
  Configurable<float> fProtonDCAzYield{"protonDCAzYield", 1.0, "[proton] DCAz cut for yield"};
  Configurable<float> fProtonDCAxyCF{"protonDCAxyCF", 0.05, "[proton] DCAxy cut for CF"};
  Configurable<float> fProtonDCAzCF{"protonDCAzCF", 1.0, "[proton] DCAz cut for CF"};
  Configurable<float> fProtonTPCTOFpT{"protonTPCTOFswitchpT", 0.7, "[proton] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fProtonTPCnsigLowYield{"protonTPCnsigmaLowPtYield", 4.0, "[proton] max TPC nsigma with low pT for yield"};
  Configurable<float> fProtonTPCnsigHighYield{"protonTPCnsigmaHighPtYield", 4.0, "[proton] max TPC nsigma with high pT for yield"};
  Configurable<float> fProtonTPCnsigLowCF{"protonTPCnsigmaLowPtCF", 4.0, "[proton] max TPC nsigma with low pT for CF"};
  Configurable<float> fProtonTPCnsigHighCF{"protonTPCnsigmaHighPtCF", 4.0, "[proton] max TPC nsigma with high pT for CF"};
  Configurable<float> fProtonTOFnsigLowYield{"protonTOFnsigmaLowPtYield", 15.0, "[proton] max TOF nsigma with low pT for yield"};
  Configurable<float> fProtonTOFnsigHighYield{"protonTOFnsigmaHighPtYield", 15.0, "[proton] max TOF nsigma with high pT yield"};
  Configurable<float> fProtonTOFnsigLowCF{"protonTOFnsigmaLowPtCF", 15.0, "[proton] max TOF nsigma with low pT for CF"};
  Configurable<float> fProtonTOFnsigHighCF{"protonTOFnsigmaHighPtCF", 15.0, "[proton] max TOF nsigma with high pT for CF"};

  // Antiproton Cuts
  Configurable<float> fAntiprotonDCAxyYield{"antiprotonDCAxyYield", 0.05, "[antiproton] DCAxy cut for yield"};
  Configurable<float> fAntiprotonDCAzYield{"antiprotonDCAzYield", 1.0, "[antiproton] DCAz cut for yield"};
  Configurable<float> fAntiprotonDCAxyCF{"antiprotonDCAxyCF", 0.05, "[antiproton] DCAxy cut for CF"};
  Configurable<float> fAntiprotonDCAzCF{"antiprotonDCAzCF", 1.0, "[antiproton] DCAz cut for CF"};
  Configurable<float> fAntiprotonTPCTOFpT{"antiprotonTPCTOFswitchpT", 0.7, "[antiproton] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fAntiprotonTPCnsigLowYield{"antiprotonTPCnsigmaLowPtYield", 4.0, "[antiproton] max TPC nsigma with low pT for yield"};
  Configurable<float> fAntiprotonTPCnsigHighYield{"antiprotonTPCnsigmaHighPtYield", 4.0, "[antiproton] max TPC nsigma with high pT for yield"};
  Configurable<float> fAntiprotonTPCnsigLowCF{"antiprotonTPCnsigmaLowPtCF", 4.0, "[antiproton] max TPC nsigma with low pT for CF"};
  Configurable<float> fAntiprotonTPCnsigHighCF{"antiprotonTPCnsigmaHighPtCF", 4.0, "[antiproton] max TPC nsigma with high pT for CF"};
  Configurable<float> fAntiprotonTOFnsigLowYield{"antiprotonTOFnsigmaLowPtYield", 15.0, "[antiproton] min TOF nsigma with low pT for yield"};
  Configurable<float> fAntiprotonTOFnsigHighYield{"antiprotonTOFnsigmaHighPtYield", 15.0, "[antiproton] min TOF nsigma with high pT for yield"};
  Configurable<float> fAntiprotonTOFnsigLowCF{"antiprotonTOFnsigmaLowPtCF", 15.0, "[antiproton] max TOF nsigma with low pT for CF"};
  Configurable<float> fAntiprotonTOFnsigHighCF{"antiprotonTOFnsigmaHighPtCF", 15.0, "[antiproton] max TOF nsigma with high pT for CF"};

  // Deuteron Cuts
  Configurable<float> fDeuteronDCAxyYield{"deuteronDCAxyYield", 0.05, "[deuteron] DCAxy cut for yield"};
  Configurable<float> fDeuteronDCAzYield{"deuteronDCAzYield", 1.0, "[deuteron] DCAz cut for yield"};
  Configurable<float> fDeuteronDCAxyCF{"deuteronDCAxyCF", 0.05, "[deuteron] DCAxy cut for CF"};
  Configurable<float> fDeuteronDCAzCF{"deuteronDCAzCF", 1.0, "[deuteron] DCAz cut for CF"};
  Configurable<float> fDeuteronTPCTOFpT{"deuteronTPCTOFswitchpT", 0.7, "[deuteron] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fDeuteronTPCnsigLowYield{"deuteronTPCnsigmaLowPtYield", 4.0, "[deuteron] max TPC nsigma with low pT for yield"};
  Configurable<float> fDeuteronTPCnsigHighYield{"deuteronTPCnsigmaHighPtYield", 4.0, "[deuteron] max TPC nsigma with high pT for yield"};
  Configurable<float> fDeuteronTPCnsigLowCF{"deuteronTPCnsigmaLowPtCF", 4.0, "[deuteron] max TPC nsigma with low pT for CF"};
  Configurable<float> fDeuteronTPCnsigHighCF{"deuteronTPCnsigmaHighPtCF", 4.0, "[deuteron] max TPC nsigma with high pT for CF"};
  Configurable<float> fDeuteronTOFnsigLowYield{"deuteronTOFnsigmaLowPtYield", 15.0, "[deuteron] min TOF nsigma with low pT for yield"};
  Configurable<float> fDeuteronTOFnsigHighYield{"deuteronTOFnsigmaHighPtYield", 15.0, "[deuteron] min TOF nsigma with high pT for yield"};
  Configurable<float> fDeuteronTOFnsigLowCF{"deuteronTOFnsigmaLowPtCF", 15.0, "[deuteron] max TOF nsigma with low pT for CF"};
  Configurable<float> fDeuteronTOFnsigHighCF{"deuteronTOFnsigmaHighPtCF", 15.0, "[deuteron] max TOF nsigma with high pT for CF"};

  // Antideuteron Cuts
  Configurable<float> fAntideuteronDCAxyYield{"antideuteronDCAxyYield", 0.05, "[antideuteron] DCAxy cut for yield"};
  Configurable<float> fAntideuteronDCAzYield{"antideuteronDCAzYield", 1.0, "[antideuteron] DCAz cut for yield"};
  Configurable<float> fAntideuteronDCAxyCF{"antideuteronDCAxyCF", 0.05, "[antideuteron] DCAxy cut for CF"};
  Configurable<float> fAntideuteronDCAzCF{"antideuteronDCAzCF", 1.0, "[antideuteron] DCAz cut for CF"};
  Configurable<float> fAntideuteronTPCTOFpT{"antideuteronTPCTOFswitchpT", 0.7, "[antideuteron] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fAntideuteronTPCnsigLowYield{"antideuteronTPCnsigmaLowPtYield", 4.0, "[antideuteron] max TPC nsigma with low pT for yield"};
  Configurable<float> fAntideuteronTPCnsigHighYield{"antideuteronTPCnsigmaHighPtYield", 4.0, "[antideuteron] max TPC nsigma with high pT for yield"};
  Configurable<float> fAntideuteronTPCnsigLowCF{"antideuteronTPCnsigmaLowPtCF", 4.0, "[antideuteron] max TPC nsigma with low pT for CF"};
  Configurable<float> fAntideuteronTPCnsigHighCF{"antideuteronTPCnsigmaHighPtCF", 4.0, "[antideuteron] max TPC nsigma with high pT for CF"};
  Configurable<float> fAntideuteronTOFnsigLowYield{"antideuteronTOFnsigmaLowPtYield", 15.0, "[antideuteron] min TOF nsigma with low pT for yield"};
  Configurable<float> fAntideuteronTOFnsigHighYield{"antideuteronTOFnsigmaHighPtYield", 15.0, "[antideuteron] min TOF nsigma with high pT for yield"};
  Configurable<float> fAntideuteronTOFnsigLowCF{"antideuteronTOFnsigmaLowPtCF", 15.0, "[antideuteron] max TOF nsigma with low pT for CF"};
  Configurable<float> fAntideuteronTOFnsigHighCF{"antideuteronTOFnsigmaHighPtCF", 15.0, "[antideuteron] max TOF nsigma with high pT for CF"};

  // Helium-3 Cuts
  Configurable<float> fHeliumDCAxy{"heliumDCAxy", 0.5, "[helium] DCAxy cut"};
  Configurable<float> fHeliumDCAz{"heliumDCAz", 1.0, "[helium] DCAz cut"};
  Configurable<float> fHeliumTPCTOFpT{"heliumTPCTOFswitchpT", 0.7, "[helium] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fHeliumTPCnsigLowYield{"heliumTPCnsigmaLowPtYield", 4.0, "[helium] max TPC nsigma with low pT for yield"};
  Configurable<float> fHeliumTPCnsigHighYield{"heliumTPCnsigmaHighPtYield", 4.0, "[helium] max TPC nsigma with high pT for yield"};
  Configurable<float> fHeliumTOFnsigLowYield{"heliumTOFnsigmaLowPtYield", 15.0, "[helium] min TOF nsigma with low pT for yield"};
  Configurable<float> fHeliumTOFnsigHighYield{"heliumTOFnsigmaHighPtYield", 15.0, "[helium] min TOF nsigma with high pT for yield"};

  // Antihelium-3 Cuts
  Configurable<float> fAntiheliumDCAxy{"antiheliumDCAxy", 0.5, "[antihelium] DCAxy cut"};
  Configurable<float> fAntiheliumDCAz{"antiheliumDCAz", 1.0, "[antihelium] DCAz cut"};
  Configurable<float> fAntiheliumTPCTOFpT{"antiheliumTPCTOFswitchpT", 0.7, "[antihelium] pT for switch in TPC/TOF nsigma"};
  Configurable<float> fAntiheliumTPCnsigLowYield{"antiheliumTPCnsigmaLowPtYield", 4.0, "[antihelium] max TPC nsigma with low pT for yield"};
  Configurable<float> fAntiheliumTPCnsigHighYield{"antiheliumTPCnsigmaHighPtYield", 4.0, "[antihelium] max TPC nsigma with high pT for yield"};
  Configurable<float> fAntiheliumTOFnsigLowYield{"antiheliumTOFnsigmaLowPtYield", 15.0, "[antihelium] min TOF nsigma with low pT for yield"};
  Configurable<float> fAntiheliumTOFnsigHighYield{"antiheliumTOFnsigmaHighPtYield", 15.0, "[antihelium] min TOF nsigma with high pT for yield"};

  Configurable<int> fBufferSize{"trackBufferSize", 2000, "Number of mixed-event tracks being stored"};

  // QC Configurables
  Configurable<float> fZVtx{"zVtx", 9999, "max zVertex"};
  Configurable<float> fRmax{"Rmax", 0.3, "Maximum radius for jet and UE regions"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  using FullTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                   aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta>;
  using FullTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                   aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  Filter prelimTrackCuts = (aod::track::itsChi2NCl < fMaxChi2ITS &&
                            aod::track::tpcChi2NCl < fMaxChi2TPC &&
                            nabs(aod::track::dcaXY) < fMaxDCAxy &&
                            nabs(aod::track::dcaZ) < fMaxDCAz &&
                            nabs(aod::track::eta) < fMaxEta);

  Preslice<FullTracksRun2> perCollisionFullTracksRun2 = o2::aod::track::collisionId;
  Preslice<FullTracksRun3> perCollisionFullTracksRun3 = o2::aod::track::collisionId;

  AxisSpecs axisSpecs;

  HistogramRegistry registryData{"dataOutput", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryQA{"dataQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Counters
    registryData.add("hNumberOfEvents", "Number of events", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hNumberOfJets", "Total number of jets", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hEventProtocol", "Event protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hTrackProtocol", "Track protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hNumPartInJet", "Number of particles in a jet", HistType::kTH1I, {{200, 0, 200}});

    // (Pseudo)Rapidity
    registryData.add("hJetRapidity", "Jet rapidity;#it{y}", HistType::kTH1F, {{200, -1, 1}});

    // pT
    registryData.add("hPtJetParticle", "p_{T} of particles in jets", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtSubtractedJet", "Subtracted jet p_{T}", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetProton", "p_{T} of p", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetAntiproton", "p_{T} of antip", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetDeuteron", "p_{T} of d", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetAntideuteron", "p_{T} of antid", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetHelium", "p_{T} of He", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetAntihelium", "p_{T} of anti-He", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{2000, 0, 500}});
    registryData.add("hPtDiff", "pT difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{100, -5, 5}});

    // nSigma
    registryData.add("hTPCsignal", "TPC signal", HistType::kTH2F, {{1000, 0, 100, "#it{p} [GeV/#it{c}]"}, {5000, 0, 5000, "d#it{E}/d#it{X} (a.u.)"}});
    registryData.add("hTOFsignal", "TOF signal", HistType::kTH2F, {{1000, 0, 100, "#it{p} [GeV/#it{c}]"}, {550, 0, 1.1, "#beta (TOF)"}});
    registryData.add("hTPCnsigmaProton", "TPC n#sigma for (anti)proton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaProton", "TOF n#sigma for (anti)proton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaDeuteron", "TPC n#sigma for (anti)deuteron", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaDeuteron", "TOF n#sigma for (anti)deuteron", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaHelium", "TPC n#sigma for (anti)helium", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaHelium", "TOF n#sigma for (anti)helium", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});

    // DCA
    registryData.add("hDCAxyFullJet", "DCA_{xy} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcaxyAxis});
    registryData.add("hDCAzFullJet", "DCA_{z} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetProton", "DCA_{z} of protons after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetAntiproton", "DCA_{z} of antiprotons after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetDeuteron", "DCA_{z} of deuterons after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetAntideuteron", "DCA_{z} of antideuterons after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetHelium", "DCA_{z} of helium after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetAntihelium", "DCA_{z} of antihelium after TPC cut", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});

    // Angular Distributions
    registryData.add("hDeltaPhiSEFull", "#Delta#varphi of particles in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEJet", "#Delta#varphi of jet particles in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEProton", "#Delta#varphi of protons in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEAntiproton", "#Delta#varphi of antiprotons in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEDeuteron", "#Delta#varphi of deuterons in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEAntideuteron", "#Delta#varphi of antideuterons in single event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEFull", "#Delta#varphi of particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEJet", "#Delta#varphi of jet particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEProton", "#Delta#varphi of protons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEAntiproton", "#Delta#varphi of antiprotons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEDeuteron", "#Delta#varphi of deuterons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEAntideuteron", "#Delta#varphi of antideuterons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiEtaSEFull", "#Delta#varphi vs #Delta#eta of full particles in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEJet", "#Delta#varphi vs #Delta#eta of jet particles in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEProton", "#Delta#varphi vs #Delta#eta of protons in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEDeuteron", "#Delta#varphi vs #Delta#eta of deuterons in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEAntideuteron", "#Delta#varphi vs #Delta#eta of antideuterons in single event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEFull", "#Delta#varphi vs #Delta#eta of particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEJet", "#Delta#varphi vs #Delta#eta of jet particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEProton", "#Delta#varphi vs #Delta#eta of protons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEDeuteron", "#Delta#varphi vs #Delta#eta of deuterons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEAntideuteron", "#Delta#varphi vs #Delta#eta of antideuterons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});

    registryData.add("hJetConeRadius", "Jet Radius;#it{R}", HistType::kTH1F, {{100, 0, 1}});

    // QA
    registryQA.add("hTOFmass", "TOF mass", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 5, "#it{m} [GeV/#it{c}^{2}]"}});
    registryQA.get<TH2>(HIST("hTOFmass"))->Sumw2();
    registryQA.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {axisSpecs.ptAxisPos});
    registryQA.add("hEtaFullEvent", "Particle pseudorapidity;#eta", HistType::kTH1F, {{200, -1, 1}});
    registryQA.add("hCrossedRowsTPC", "Crossed rows TPC", HistType::kTH2I, {axisSpecs.ptAxisPos, {135, 65, 200}});
    registryQA.add("hClusterITS", "ITS clusters", HistType::kTH2I, {axisSpecs.ptAxisPos, {10, 0, 10}});
    registryQA.add("hClusterTPC", "TPC clusters", HistType::kTH2I, {axisSpecs.ptAxisPos, {135, 65, 200}});
    registryQA.add("hRatioCrossedRowsTPC", "Ratio crossed rows/findable TPC", HistType::kTH2F, {axisSpecs.ptAxisPos, {100, 0.5, 1.5}});
    registryQA.add("hChi2ITS", "ITS #chi^{2}", HistType::kTH2F, {axisSpecs.ptAxisPos, {400, 0, 40}});
    registryQA.add("hChi2TPC", "TPC #chi^{2}", HistType::kTH2F, {axisSpecs.ptAxisPos, {50, 0, 5}});
    registryQA.add("hDCAxyFullEvent", "DCA_{xy} of full event", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcaxyAxis});
    registryQA.add("hDCAzFullEvent", "DCA_{z} of full event", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryQA.add("hJetPtVsNumPart", "Total jet p_{T} vs number of constituents", HistType::kTH2F, {axisSpecs.ptAxisPos, {100, 0, 100}});

    // QA Histograms for Comparison with nuclei_in_jets.cxx
    registryQA.add("hMultiplicityJetPlusUE", "hMultiplicityJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQA.add("hMultiplicityJet", "hMultiplicityJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQA.add("hMultiplicityUE", "hMultiplicityUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
    registryQA.add("hPtLeading", "hPtLeading", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hEtaLeading", "hEtaLeading", HistType::kTH1F, {{100, -0.8, 0.8, "#eta"}});
    registryQA.add("hPhiLeading", "hPhiLeading", HistType::kTH1F, {{100, 0, TMath::TwoPi(), "#phi"}});
    registryQA.add("hRJet", "hRJet", HistType::kTH1F, {{100, 0.0, 0.5, "#it{R}"}});
    registryQA.add("hRUE", "hRUE", HistType::kTH1F, {{100, 0.0, 0.5, "#it{R}"}});
    registryQA.add("hAngleJetLeadingTrack", "hAngleJetLeadingTrack", HistType::kTH1F, {{200, 0.0, 50.0, "#theta"}});
    registryQA.add("hPtJetPlusUE", "hPtJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hPtJet", "hPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hPtUE", "hPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hDeltaEtadeltaPhiJet", "hDeltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaEtadeltaPhiUE", "hDeltaEtadeltaPhiUE", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaEtadeltaPhiLeadingJet", "hDeltaEtadeltaPhiLeadingJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaJetPt", "hDeltaJetPt", HistType::kTH1F, {{200, -2, 2, "#Delta#it{p}_{T} (GeV/#it{c})"}});

    // QC Histograms for ptJet < ptLeading
    registryQA.add("hNParticlesClusteredInJet", "hNParticlesClusteredInJet", HistType::kTH1F, {{50, 0, 50, "#it{N}_{ch}"}});
    registryQA.add("hPtParticlesClusteredInJet", "hPtParticlesClusteredInJet", HistType::kTH1F, {{200, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hDeltaEtaDeltaPhiJetAxis", "hDeltaEtaDeltaPhiJetAxis", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaEtaDeltaPhiJetAxisLeading", "hDeltaEtaDeltaPhiJetAxisLeading", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
  }

  std::vector<std::pair<double, double>> fBufferProton;
  std::vector<std::pair<double, double>> fBufferAntiproton;
  std::vector<std::pair<double, double>> fBufferDeuteron;
  std::vector<std::pair<double, double>> fBufferAntideuteron;
  std::vector<std::pair<double, double>> fBufferJet;
  std::vector<std::pair<double, double>> fBufferFull;

  template <class Bc>
  void initCCDB(Bc const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
  }

  template <class T>
  bool selectTrack(T const& track)
  {
    if (track.tpcNClsCrossedRows() < fMinRatioCrossedRowsTPC * track.tpcNClsFindable() ||
        track.tpcNClsCrossedRows() < fMinNCrossedRowsTPC ||
        track.tpcNClsFound() < fMinReqClusterTPC ||
        track.itsNCls() < fMinReqClusterITS) {
      return false;
    }
    if (doprocessRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit) ||
          !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool isProton(const T& track, bool tightCuts)
  {
    if (track.sign() < 0)
      return false;

    if (tightCuts) { // for correlation function
      // TPC
      if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fProtonTPCnsigLowCF)
        return false;
      if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fProtonTPCnsigHighCF)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fProtonDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > fProtonDCAzCF)
        return false;

      // TOF
      if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fProtonTOFnsigLowCF)
        return false; // drop for low pt if no TOF
      if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fProtonTOFnsigHighCF)
        return false;
    } else { // for yields
      // TPC
      if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fProtonTPCnsigLowYield)
        return false;
      if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fProtonTPCnsigHighYield)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fProtonDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > fProtonDCAzYield)
        return false;

      // TOF
      if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fProtonTOFnsigLowYield)
        return false;
      if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fProtonTOFnsigHighYield)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isAntiproton(const T& track, bool tightCuts)
  {
    if (track.sign() > 0)
      return false;

    if (tightCuts) { // for correlation function
      // TPC
      if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fAntiprotonTPCnsigLowCF)
        return false;
      if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fAntiprotonTPCnsigHighCF)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fAntiprotonDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > fAntiprotonDCAzCF)
        return false;

      // TOF
      if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fAntiprotonTOFnsigLowCF)
        return false;
      if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fAntiprotonTOFnsigHighCF)
        return false;
    } else { // for yields
      // TPC
      if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fAntiprotonTPCnsigLowYield)
        return false;
      if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > fAntiprotonTPCnsigHighYield)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fAntiprotonDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > fAntiprotonDCAzYield)
        return false;

      // TOF
      if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fAntiprotonTOFnsigLowYield)
        return false;
      if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > fAntiprotonTOFnsigHighYield)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isDeuteron(const T& track, bool tightCuts)
  {
    if (track.sign() < 0)
      return false;

    if (tightCuts) { // for correlation function
      // TPC
      if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fDeuteronTPCnsigLowCF)
        return false;
      if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fDeuteronTPCnsigHighCF)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fDeuteronDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > fDeuteronDCAzCF)
        return false;

      // TOF
      if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fDeuteronTOFnsigLowCF)
        return false;
      if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fDeuteronTOFnsigHighCF)
        return false;
    } else { // for yields
      // TPC
      if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fDeuteronTPCnsigLowYield)
        return false;
      if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fDeuteronTPCnsigHighYield)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fDeuteronDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > fDeuteronDCAzYield)
        return false;

      // TOF
      if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fDeuteronTOFnsigLowYield)
        return false;
      if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fDeuteronTOFnsigHighYield)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isAntideuteron(const T& track, bool tightCuts)
  {
    if (track.sign() > 0)
      return false;

    if (tightCuts) { // for correlation function
      // TPC
      if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fAntideuteronTPCnsigLowCF)
        return false;
      if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fAntideuteronTPCnsigHighCF)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fAntideuteronDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > fAntideuteronDCAzCF)
        return false;

      // TOF
      if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fAntideuteronTOFnsigLowCF)
        return false;
      if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fAntideuteronTOFnsigHighCF)
        return false;
    } else { // for yields
      // TPC
      if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fAntideuteronTPCnsigLowYield)
        return false;
      if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > fAntideuteronTPCnsigHighYield)
        return false;

      // DCA
      if (TMath::Abs(track.dcaXY()) > fAntideuteronDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > fAntideuteronDCAzYield)
        return false;

      // TOF
      if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fAntideuteronTOFnsigLowYield)
        return false;
      if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > fAntideuteronTOFnsigHighYield)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isHelium(const T& track)
  {
    if (track.sign() < 0)
      return false;

    // TPC
    if (track.pt() < fHeliumTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > fHeliumTPCnsigLowYield)
      return false;
    if (track.pt() > fHeliumTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > fHeliumTPCnsigHighYield)
      return false;

    // DCA
    if (TMath::Abs(track.dcaXY()) > fHeliumDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fHeliumDCAz)
      return false;

    // TOF
    if (track.pt() < fHeliumTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > fHeliumTOFnsigLowYield)
      return false;
    if (track.pt() > fHeliumTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > fHeliumTOFnsigHighYield)
      return false;

    return true;
  }

  template <typename T>
  bool isAntihelium(const T& track)
  {
    if (track.sign() > 0)
      return false;

    // TPC
    if (track.pt() < fAntiheliumTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > fAntiheliumTPCnsigLowYield)
      return false;
    if (track.pt() > fAntiheliumTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > fAntiheliumTPCnsigHighYield)
      return false;

    // DCA
    if (TMath::Abs(track.dcaXY()) > fAntiheliumDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fAntiheliumDCAz)
      return false;

    // TOF
    if (track.pt() < fAntiheliumTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > fAntiheliumTOFnsigLowYield)
      return false;
    if (track.pt() > fAntiheliumTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > fAntiheliumTOFnsigHighYield)
      return false;

    return true;
  }

  void setTrackBuffer(const auto& tempBuffer, auto& buffer)
  {
    for (const auto& pair : tempBuffer) {
      if (static_cast<int>(buffer.size()) == fBufferSize) {
        buffer.insert(buffer.begin(), pair);
        buffer.resize(fBufferSize);
      } else if (static_cast<int>(buffer.size()) < fBufferSize) {
        buffer.emplace_back(pair);
      }
    }
  }

  void fillMixedEventDeltas(const auto& track, const auto& buffer, int particleType, const TVector3 jetAxis)
  {
    if (buffer.size() == 0)
      return;
    if (std::isnan(track.phi()) || std::isnan(jetAxis.Phi()))
      return;
    for (int i = 0; i < static_cast<int>(buffer.size()); i++) { // loop over tracks in buffer
      if (std::isnan(buffer.at(i).first))
        continue;
      if (buffer.at(i).first > 2 * TMath::Pi() || buffer.at(i).first < -2 * TMath::Pi()) {
        registryData.fill(HIST("hTrackProtocol"), 16);
        continue;
      }

      double phiToAxis = TVector2::Phi_0_2pi(track.phi() - jetAxis.Phi());
      double etaToAxis = track.eta() - jetAxis.Eta();
      double DeltaPhi = TVector2::Phi_0_2pi(phiToAxis - buffer.at(i).first);
      if (DeltaPhi > (1.5 * TMath::Pi())) { // ensure range of [-pi/2, 3/2 pi]
        DeltaPhi = DeltaPhi - 2 * TMath::Pi();
      }
      double DeltaEta = etaToAxis - buffer.at(i).second;

      switch (particleType) {
        case -1:
          registryData.fill(HIST("hDeltaPhiMEFull"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEFull"), DeltaPhi, DeltaEta);
          break;
        case 0:
          registryData.fill(HIST("hDeltaPhiMEJet"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEJet"), DeltaPhi, DeltaEta);
          break;
        case 1:
          registryData.fill(HIST("hDeltaPhiMEProton"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEProton"), DeltaPhi, DeltaEta);
          break;
        case 2:
          registryData.fill(HIST("hDeltaPhiMEAntiproton"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEAntiproton"), DeltaPhi, DeltaEta);
          break;
        case 3:
          registryData.fill(HIST("hDeltaPhiMEDeuteron"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEDeuteron"), DeltaPhi, DeltaEta);
          break;
        case 4:
          registryData.fill(HIST("hDeltaPhiMEAntideuteron"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEAntideuteron"), DeltaPhi, DeltaEta);
          break;
      }
    } // for (int i = 0; i < static_cast<int>(buffer.size()); i++)
  }

  void doCorrelations(const auto& particleVector, const auto& buffer, auto& tempBuffer, int particleType, const TVector3 jetAxis)
  {
    if (std::isnan(jetAxis.Phi()))
      return;
    for (int i = 0; i < static_cast<int>(particleVector.size()); i++) {
      if (std::isnan(particleVector.at(i).phi()))
        continue;
      double phiToAxis = TVector2::Phi_0_2pi(particleVector.at(i).phi() - jetAxis.Phi());
      double etaToAxis = particleVector.at(i).eta() - jetAxis.Eta();
      if (TMath::Abs(particleVector.at(i).phi()) > 2 * TMath::Pi()) {
        registryData.fill(HIST("hTrackProtocol"), 14);
        continue;
      }
      for (int j = i + 1; j < static_cast<int>(particleVector.size()); j++) {
        if ((j == static_cast<int>(particleVector.size())) || std::isnan(particleVector.at(j).phi()))
          continue;
        if (TMath::Abs(particleVector.at(j).phi()) > 2 * TMath::Pi()) {
          registryData.fill(HIST("hTrackProtocol"), 15);
          continue;
        }

        double DeltaPhi = TVector2::Phi_0_2pi(particleVector.at(i).phi() - particleVector.at(j).phi());
        double DeltaEta = particleVector.at(i).eta() - particleVector[j].eta();
        if (DeltaPhi > (1.5 * TMath::Pi())) {
          DeltaPhi = DeltaPhi - 2 * TMath::Pi();
        }
        switch (particleType) {
          case -1:
            registryData.fill(HIST("hDeltaPhiSEFull"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEFull"), DeltaPhi, DeltaEta);
            break;
          case 0:
            registryData.fill(HIST("hDeltaPhiSEJet"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEJet"), DeltaPhi, DeltaEta);
            break;
          case 1:
            registryData.fill(HIST("hDeltaPhiSEProton"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEProton"), DeltaPhi, DeltaEta);
            break;
          case 2:
            registryData.fill(HIST("hDeltaPhiSEAntiproton"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEAntiproton"), DeltaPhi, DeltaEta);
            break;
          case 3:
            registryData.fill(HIST("hDeltaPhiSEDeuteron"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEDeuteron"), DeltaPhi, DeltaEta);
            break;
          case 4:
            registryData.fill(HIST("hDeltaPhiSEAntideuteron"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEAntideuteron"), DeltaPhi, DeltaEta);
            break;
        }
      }
      fillMixedEventDeltas(particleVector.at(i), buffer, particleType, jetAxis);
      tempBuffer.emplace_back(std::make_pair(phiToAxis, etaToAxis));
    }
  }

  double getDeltaPhi(double a1, double a2)
  {
    if (std::isnan(a1) || std::isnan(a2) || a1 == -999 || a2 == -999)
      return -999;
    double deltaPhi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = TMath::Abs(phi1 - phi2);

    if (diff <= TMath::Pi())
      deltaPhi = diff;
    if (diff > TMath::Pi())
      deltaPhi = TMath::TwoPi() - diff;

    return deltaPhi;
  }

  void getPerpendicularAxis(TVector3 p, TVector3& u, double sign)
  {
    // Initialization
    double ux(0), uy(0), uz(0);

    // Components of Vector p
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    // Protection 1
    if (px == 0 && py != 0) {

      uy = -(pz * pz) / py;
      ux = sign * sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Protection 2
    if (py == 0 && px != 0) {

      ux = -(pz * pz) / px;
      uy = sign * sqrt(px * px - (pz * pz * pz * pz) / (px * px));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Equation Parameters
    double a = px * px + py * py;
    double b = 2.0 * px * pz * pz;
    double c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    double delta = b * b - 4.0 * a * c;

    // Protection agains delta<0
    if (delta < 0) {
      return;
    }

    // Solutions
    ux = (-b + sign * sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }

  template <typename U>
  void fillHistograms(U const& tracks)
  {
    std::vector<std::pair<double, double>> fTempBufferProton;
    std::vector<std::pair<double, double>> fTempBufferAntiproton;
    std::vector<std::pair<double, double>> fTempBufferDeuteron;
    std::vector<std::pair<double, double>> fTempBufferAntideuteron;
    std::vector<std::pair<double, double>> fTempBufferJet;
    std::vector<std::pair<double, double>> fTempBufferFull;
    fTempBufferProton.clear();
    fTempBufferAntiproton.clear();
    fTempBufferDeuteron.clear();
    fTempBufferAntideuteron.clear();
    fTempBufferJet.clear();
    fTempBufferFull.clear();
    std::vector<fastjet::PseudoJet> jetInput;
    std::map<int, typename U::iterator> particles;
    std::vector<typename U::iterator> particlesForCF;
    jetInput.clear();
    particles.clear();
    int index = 0;
    int leadingID = 0;
    fastjet::PseudoJet hardestJet(0., 0., 0., 0.);
    fastjet::PseudoJet subtractedJet(0., 0., 0., 0.);
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> constituents;
    jets.clear();
    constituents.clear();

    for (const auto& track : tracks) {
      if (!selectTrack(track))
        continue;

      double mass;
      if (track.hasTOF()) {
        mass = track.mass(); // check reliability, maybe use only pion mass
        registryQA.fill(HIST("hTOFmass"), track.pt(), track.mass());
        registryData.fill(HIST("hTrackProtocol"), 1);
      } else {
        mass = 0.139; // pion mass as default, ~80% are pions
        registryData.fill(HIST("hTrackProtocol"), 2);
      }

      if (track.pt() > fMinLeadingPt) {
        leadingID = track.globalIndex();
      }

      if (track.tpcNClsFindable() != 0) {
        registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
      }
      registryQA.fill(HIST("hPtFullEvent"), track.pt());
      registryQA.fill(HIST("hEtaFullEvent"), track.eta());
      registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
      registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
      registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
      registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
      registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(mass));
      inputPseudoJet.set_user_index(index);
      particles[index] = track;
      particlesForCF.emplace_back(track);
      jetInput.emplace_back(inputPseudoJet);

      index++;
    } // for (const auto& track : tracks)

    if (jetInput.size() < 2)
      return;
    registryData.fill(HIST("hEventProtocol"), 2);

    // Reconstruct Jets
    double ghost_maxrap = 1.0;
    double ghost_area = 0.005;
    int ghost_repeat = 1;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, fJetR);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap, ghost_repeat, ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());

    if (jets.size() == 0)
      return;

    registryData.fill(HIST("hEventProtocol"), 3);

    hardestJet = jets[0];

    if (hardestJet.pt() < fMinJetPt)
      return;

    registryData.fill(HIST("hEventProtocol"), 4);
    if (hardestJet.constituents().size() < 2) // unlikely but maybe yields before this point?
      return;

    registryData.fill(HIST("hEventProtocol"), 5);
    registryData.fill(HIST("hNumberOfJets"), 0);
    registryData.fill(HIST("hPtTotalJet"), hardestJet.pt());
    registryData.fill(HIST("hJetRapidity"), hardestJet.rap());
    registryData.fill(HIST("hNumPartInJet"), hardestJet.constituents().size());
    registryQA.fill(HIST("hJetPtVsNumPart"), hardestJet.pt(), hardestJet.constituents().size());

    for (const auto& constituent : hardestJet.constituents()) {
      registryData.fill(HIST("hPtJetParticle"), constituent.pt());
      if (std::isnan(constituent.phi()) || std::isnan(hardestJet.phi()))
        continue;
      double DeltaPhi = TVector2::Phi_0_2pi(constituent.phi() - hardestJet.phi());
      double DeltaEta = constituent.eta() - hardestJet.eta();
      double Delta = TMath::Sqrt(DeltaPhi * DeltaPhi + DeltaEta * DeltaEta);
      registryData.fill(HIST("hJetConeRadius"), Delta);
    }

    fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2)); // TODO: fix subtraction -> check PWGJE
    fastjet::JetMedianBackgroundEstimator bkgEst(selector, jetDefBkg, areaDefBkg);
    fastjet::Subtractor subtractor(&bkgEst);
    subtractor.set_use_rho_m(true);
    bkgEst.set_particles(jetInput);

    subtractedJet = subtractor(hardestJet);
    if (subtractedJet.has_constituents()) {
      for (const auto& subConstituent : subtractedJet.constituents()) {
        registryData.fill(HIST("hPtSubtractedJet"), subConstituent.pt());
      }
    }

    if (!hardestJet.has_constituents())
      return;
    constituents = hardestJet.constituents();

    // QA for comparison with nuclei_in_jets
    const auto& leadingTrack = tracks.iteratorAt(leadingID);
    TVector3 pLeading(leadingTrack.px(), leadingTrack.py(), leadingTrack.pz());
    TVector3 pJet(hardestJet.px(), hardestJet.py(), hardestJet.pz());
    TVector3 UEAxis1(0.0, 0.0, 0.0);
    TVector3 UEAxis2(0.0, 0.0, 0.0);
    getPerpendicularAxis(pJet, UEAxis1, +1.0);
    getPerpendicularAxis(pJet, UEAxis2, -1.0);

    doCorrelations(particlesForCF, fBufferFull, fTempBufferFull, -1, pJet);
    setTrackBuffer(fTempBufferFull, fBufferFull);

    if (UEAxis1.Mag() == 0 || UEAxis2.Mag() == 0)
      return;
    double deltaEta = pLeading.Eta() - pJet.Eta();
    double deltaPhi = getDeltaPhi(pLeading.Phi(), pJet.Phi());
    registryQA.fill(HIST("hPtLeading"), leadingTrack.pt());
    registryQA.fill(HIST("hPhiLeading"), leadingTrack.phi());
    registryQA.fill(HIST("hEtaLeading"), leadingTrack.eta());
    registryQA.fill(HIST("hAngleJetLeadingTrack"), (180.0 / TMath::Pi()) * pLeading.Angle(pJet));
    if (deltaPhi != -999)
      registryQA.fill(HIST("hDeltaEtadeltaPhiLeadingJet"), deltaEta, deltaPhi);

    double NchJetPlusUE(0);
    double NchJet(0);
    double NchUE(0);
    double ptJetPlusUE(0);
    double ptJet(0);
    double ptUE(0);

    for (const auto& [index, track] : particles) {
      TVector3 particleDir(track.px(), track.py(), track.pz());
      double deltaEtaJet = particleDir.Eta() - pJet.Eta();
      double deltaPhiJet = getDeltaPhi(particleDir.Phi(), pJet.Phi());
      double deltaRJet = sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
      double deltaEtaUE1 = particleDir.Eta() - UEAxis1.Eta();
      double deltaPhiUE1 = getDeltaPhi(particleDir.Phi(), UEAxis1.Phi());
      double deltaRUE1 = sqrt(deltaEtaUE1 * deltaEtaUE1 + deltaPhiUE1 * deltaPhiUE1);
      double deltaEtaUE2 = particleDir.Eta() - UEAxis2.Eta();
      double deltaPhiUE2 = getDeltaPhi(particleDir.Phi(), UEAxis2.Phi());
      double deltaRUE2 = sqrt(deltaEtaUE2 * deltaEtaUE2 + deltaPhiUE2 * deltaPhiUE2);

      if (deltaRJet < fRmax) {
        if (deltaPhiJet != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiJet"), deltaEtaJet, deltaPhiJet);
        registryQA.fill(HIST("hRJet"), deltaRJet);
        NchJetPlusUE++;
        ptJetPlusUE = ptJetPlusUE + track.pt();
      }
      if (deltaRUE1 < fRmax) {
        if (deltaPhiUE1 != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiUE"), deltaEtaUE1, deltaPhiUE1);
        registryQA.fill(HIST("hRUE"), deltaRUE1);
        NchUE++;
        ptUE = ptUE + track.pt();
      }
      if (deltaRUE2 < fRmax) {
        if (deltaPhiUE2 != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiUE"), deltaEtaUE2, deltaPhiUE2);
        registryQA.fill(HIST("hRUE"), deltaRUE2);
        NchUE++;
        ptUE = ptUE + track.pt();
      }
    } // for (const auto& [index, track] : particles)

    NchJet = NchJetPlusUE - 0.5 * NchUE;
    ptJet = ptJetPlusUE - 0.5 * ptUE;
    registryQA.fill(HIST("hMultiplicityJetPlusUE"), NchJetPlusUE);
    registryQA.fill(HIST("hMultiplicityJet"), NchJet);
    registryQA.fill(HIST("hMultiplicityUE"), 0.5 * NchUE);
    registryQA.fill(HIST("hPtJetPlusUE"), ptJetPlusUE);
    registryQA.fill(HIST("hPtJet"), ptJet);
    registryQA.fill(HIST("hPtUE"), 0.5 * ptUE);
    registryQA.fill(HIST("hDeltaJetPt"), hardestJet.pt() - ptJetPlusUE);

    int nPartClusteredJet = static_cast<int>(constituents.size());

    // Fill QA Histograms
    if (ptJetPlusUE < fMinJetPt) {

      registryQA.fill(HIST("hNParticlesClusteredInJet"), nPartClusteredJet);
      double dEta = pLeading.Eta() - pJet.Eta();
      double dPhi = getDeltaPhi(pLeading.Phi(), pJet.Phi());
      if (dPhi != -999)
        registryQA.fill(HIST("hDeltaEtaDeltaPhiJetAxisLeading"), dEta, dPhi);

      for (const auto& track : constituents) {
        TVector3 particleDir(track.px(), track.py(), track.pz());
        double dEta = particleDir.Eta() - pJet.Eta();
        double dPhi = getDeltaPhi(particleDir.Phi(), pJet.Phi());
        registryQA.fill(HIST("hPtParticlesClusteredInJet"), track.pt());
        if (dPhi != -999)
          registryQA.fill(HIST("hDeltaEtaDeltaPhiJetAxis"), dEta, dPhi);
      }
    }

    std::vector<typename U::iterator> jetProtons;
    std::vector<typename U::iterator> jetAntiprotons;
    std::vector<typename U::iterator> jetDeuterons;
    std::vector<typename U::iterator> jetAntideuterons;
    std::vector<typename U::iterator> jetHelium;
    std::vector<typename U::iterator> jetAntihelium;
    std::vector<typename U::iterator> jetAll;

    for (int i = 0; i < static_cast<int>(constituents.size()); i++) {
      registryData.fill(HIST("hTrackProtocol"), 3);
      fastjet::PseudoJet pseudoParticle = constituents.at(i);
      int id = pseudoParticle.user_index();
      const auto& jetParticle = particles[id];
      jetAll.emplace_back(jetParticle);

      registryData.fill(HIST("hDCAxyFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("hDCAzFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("hTPCsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
      if (jetParticle.hasTOF())
        registryData.fill(HIST("hTOFsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());

      double ptDiff = pseudoParticle.pt() - jetParticle.pt();
      registryData.fill(HIST("hPtDiff"), ptDiff);

      if (jetParticle.pt() < fMinJetParticlePt)
        continue;
      if (isProton(jetParticle, false)) { // collect protons in jet
        registryData.fill(HIST("hPtJetProton"), jetParticle.pt());
        registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());
        registryData.fill(HIST("hTrackProtocol"), 4); // # protons
        if (isProton(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 5); // # high purity protons
          jetProtons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetProton"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isAntiproton(jetParticle, false)) { // collect antiprotons in jet
        registryData.fill(HIST("hPtJetAntiproton"), jetParticle.pt());
        registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());
        registryData.fill(HIST("hTrackProtocol"), 6); // # antiprotons
        if (isAntiproton(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 7); // # high purity antiprotons
          jetAntiprotons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetAntiproton"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isDeuteron(jetParticle, false)) { // collect deuterons in jet
        registryData.fill(HIST("hPtJetDeuteron"), jetParticle.pt());
        registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaDe());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaDe());
        registryData.fill(HIST("hTrackProtocol"), 8); // # deuterons
        if (isDeuteron(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 9); // # high purity deuterons
          jetDeuterons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetDeuteron"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isAntideuteron(jetParticle, false)) {
        registryData.fill(HIST("hPtJetAntideuteron"), jetParticle.pt());
        registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaDe());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaDe());
        registryData.fill(HIST("hTrackProtocol"), 10); // # antideuterons
        if (isAntideuteron(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 11); // # high purity antideuterons
          jetAntideuterons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetAntideuteron"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isHelium(jetParticle) || isAntihelium(jetParticle)) { // collect (anti)helium in jet
        registryData.fill(HIST("hPtJetHelium"), jetParticle.pt());
        registryData.fill(HIST("hTPCnsigmaHelium"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaHe());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaHelium"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaHe());
        if (isHelium(jetParticle)) {
          registryData.fill(HIST("hTrackProtocol"), 12); // # helium
          jetDeuterons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetHelium"), jetParticle.pt(), jetParticle.dcaZ());
        } else {
          registryData.fill(HIST("hTrackProtocol"), 13); // # antihelium
          jetAntideuterons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetAntihelium"), jetParticle.pt(), jetParticle.dcaZ());
        }
      }
    } // for (int i=0; i<static_cast<int>(constituents.size()); i++)

    if (jetAll.size() > 1) { // general correlation function
      doCorrelations(jetAll, fBufferJet, fTempBufferJet, 0, pJet);
      setTrackBuffer(fTempBufferJet, fBufferJet);
    }

    if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2) && (jetDeuterons.size() < 2) && (jetAntideuterons.size() < 2))
      return;
    registryData.fill(HIST("hEventProtocol"), 6);

    if (jetProtons.size() > 1) {
      doCorrelations(jetProtons, fBufferProton, fTempBufferProton, 1, pJet);
      setTrackBuffer(fTempBufferProton, fBufferProton);
    }
    if (jetAntiprotons.size() > 1) {
      doCorrelations(jetAntiprotons, fBufferAntiproton, fTempBufferAntiproton, 2, pJet);
      setTrackBuffer(fTempBufferAntiproton, fBufferAntiproton);
    }
    if (jetDeuterons.size() > 1) {
      doCorrelations(jetDeuterons, fBufferDeuteron, fTempBufferDeuteron, 3, pJet);
      setTrackBuffer(fTempBufferDeuteron, fBufferDeuteron);
    }
    if (jetAntideuterons.size() > 1) {
      doCorrelations(jetAntideuterons, fBufferAntideuteron, fTempBufferAntideuteron, 4, pJet);
      setTrackBuffer(fTempBufferAntideuteron, fBufferAntideuteron);
    }
  }

  void processRun2(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun2> const& tracks,
                   BCsWithRun2Info const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      registryData.fill(HIST("hEventProtocol"), 0);
      registryData.fill(HIST("hNumberOfEvents"), 0);
      if (!collision.alias_bit(kINT7))
        continue;
      registryData.fill(HIST("hEventProtocol"), 1);

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun2, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun2, "process Run 2 data", true);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun3> const& tracks)
  {
    for (const auto& collision : collisions) {
      registryData.fill(HIST("hEventProtocol"), 0);
      registryData.fill(HIST("hNumberOfEvents"), 0);
      if (!collision.sel8())
        continue;
      registryData.fill(HIST("hEventProtocol"), 1);
      if (TMath::Abs(collision.posZ()) > fZVtx)
        continue;

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun3, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun3, "process Run 3 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc)};
}
