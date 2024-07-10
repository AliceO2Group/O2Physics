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
// #include "PWGJE/Core/JetFinder.h"
// #include "PWGJE/Core/JetFindingUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AngularCorrelationsInJets {
  // Preliminary Cuts
  Configurable<int> fMinNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of crossed rows TPC"};
  Configurable<int> fMinReqClusterITS{"minReqClusterITS", 2, "min number of clusters required in ITS"};
  Configurable<int> fMinReqClusterTPC{"minReqClusterTPC", 70, "min number of clusters required in TPC"};
  Configurable<float> fMinRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.7f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> fMaxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  Configurable<float> fMaxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> fMaxDCAxy{"maxDCA_xy", 0.5f, "max DCA to vertex xy"};
  Configurable<float> fMaxDCAz{"maxDCA_z", 2.4f, "max DCA to vertex z"};
  Configurable<float> fMaxEta{"maxEta", 0.8, "max pseudorapidity"}; // consider jet cone?

  // Jet Cuts
  Configurable<float> fJetR{"jetR", 0.4, "jet resolution parameter"};
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

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  using FullTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                   aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFmass, aod::pidTOFbeta>;
  using FullTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                   aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFmass, aod::pidTOFbeta>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  Filter prelimTrackCuts = (aod::track::itsChi2NCl < fMaxChi2ITS &&
                            aod::track::tpcChi2NCl < fMaxChi2TPC &&
                            nabs(aod::track::dcaXY) < fMaxDCAxy &&
                            nabs(aod::track::dcaZ) < fMaxDCAz &&
                            nabs(aod::track::eta) < fMaxEta);

  Preslice<FullTracksRun2> perCollisionFullTracksRun2 = o2::aod::track::collisionId;
  Preslice<FullTracksRun3> perCollisionFullTracksRun3 = o2::aod::track::collisionId;

  AxisSpec ptAxis = {1000, 0, 100, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec particleTypeAxis = {4, 1, 5, "[p, ap, d, ad]"};
  AxisSpec nsigmapTAxis = {1000, -50, 50, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec nsigmaAxis = {1000, -5, 5, "n#sigma"};
  AxisSpec dcazAxis = {200, -3, 3, "DCA_{z} [cm]"};
  AxisSpec dcaxyAxis = {200, -2, 2, "DCA_{xy} [cm]"};
  AxisSpec angDistPhiAxis = {1000, -2, 5, "#Delta#varphi"};
  AxisSpec angDistEtaAxis = {1000, -2, 2, "#Delta#eta"};

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
    registryData.add("hPtJetParticle", "p_{T} of particles in jets", HistType::kTH1D, {ptAxis});
    registryData.add("hPtSubtractedJet", "Subtracted jet p_{T}", HistType::kTH1D, {ptAxis});
    registryData.add("hPtJetProtonDeuteron", "p_{T} of (anti)p, (anti)d", HistType::kTH2D, {particleTypeAxis, ptAxis});
    registryData.add("hPtTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{2000, 0, 500}});
    registryData.add("hPtDiff", "pT difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{100, -5, 5}});

    // nSigma
    registryData.add("hTPCsignal", "TPC signal", HistType::kTH2F, {{1000, 0, 100, "#it{p} [GeV/#it{c}]"}, {5000, 0, 5000, "d#it{E}/d#it{X} (a.u.)"}});
    registryData.add("hTOFsignal", "TOF signal", HistType::kTH2F, {{1000, 0, 100, "#it{p} [GeV/#it{c}]"}, {550, 0, 1.1, "#beta (TOF)"}});
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

    // QA
    registryQA.add("hTOFmass", "TOF mass", HistType::kTH2F, {ptAxis, {1000, 0, 5, "#it{m} [GeV/#it{c}^{2}]"}});
    registryQA.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {ptAxis});
    registryQA.add("hEtaFullEvent", "Particle pseudorapidity;#eta", HistType::kTH1F, {{200, -1, 1}});
    registryQA.get<TH2>(HIST("hTOFmass"))->Sumw2();
    registryQA.add("hCrossedRowsTPC", "Crossed rows TPC", HistType::kTH2I, {ptAxis, {135, 65, 200}});
    registryQA.add("hClusterITS", "ITS clusters", HistType::kTH2I, {ptAxis, {10, 0, 10}});
    registryQA.add("hClusterTPC", "TPC clusters", HistType::kTH2I, {ptAxis, {135, 65, 200}});
    registryQA.add("hRatioCrossedRowsTPC", "Ratio crossed rows/findable TPC", HistType::kTH2F, {ptAxis, {100, 0.5, 1.5}});
    registryQA.add("hChi2ITS", "ITS #chi^{2}", HistType::kTH2F, {ptAxis, {400, 0, 40}});
    registryQA.add("hChi2TPC", "TPC #chi^{2}", HistType::kTH2F, {ptAxis, {50, 0, 5}});
    registryQA.add("hDCAxyFullEvent", "DCA_{xy} of full event", HistType::kTH2F, {ptAxis, dcaxyAxis});
    registryQA.add("hDCAzFullEvent", "DCA_{z} of full event", HistType::kTH2F, {ptAxis, dcazAxis});
  }

  std::vector<typename FullTracksRun2::iterator> fTrackBufferProtonRun2;
  std::vector<typename FullTracksRun2::iterator> fTrackBufferAntiprotonRun2;
  std::vector<typename FullTracksRun2::iterator> fTrackBufferDeuteronRun2;
  std::vector<typename FullTracksRun2::iterator> fTrackBufferAntideuteronRun2;
  std::vector<typename FullTracksRun3::iterator> fTrackBufferProtonRun3;
  std::vector<typename FullTracksRun3::iterator> fTrackBufferAntiprotonRun3;
  std::vector<typename FullTracksRun3::iterator> fTrackBufferDeuteronRun3;
  std::vector<typename FullTracksRun3::iterator> fTrackBufferAntideuteronRun3;
  // TODO: check if FullTracksRun2 works for Run3 too or add Run3 track buffers

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
      if (!(track.trackType() & o2::aod::track::Run2Track) //||
                                                           //!(track.flags() & o2::aod::track::TPCrefit) ||
                                                           //!(track.flags() & o2::aod::track::ITSrefit)
      ) {
        return false;
      }
    }
    return true;
  }

  std::vector<fastjet::PseudoJet> findJets(std::vector<fastjet::PseudoJet> jetInput)
  {
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> constituents;
    jets.clear();
    constituents.clear();
    fastjet::PseudoJet hardestJet(0., 0., 0., 0.);
    fastjet::PseudoJet subtractedJet(0., 0., 0., 0.);

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
      return constituents;

    registryData.fill(HIST("hEventProtocol"), 3);

    hardestJet = jets[0];

    if (hardestJet.pt() < fMinJetPt)
      return constituents;
    registryData.fill(HIST("hEventProtocol"), 4);
    if (hardestJet.constituents().size() < 2)
      return constituents;
    registryData.fill(HIST("hEventProtocol"), 5);
    registryData.fill(HIST("hNumberOfJets"), 0);
    registryData.fill(HIST("hPtTotalJet"), hardestJet.pt());
    registryData.fill(HIST("hJetRapidity"), hardestJet.rap());
    registryData.fill(HIST("hNumPartInJet"), hardestJet.constituents().size());

    for (const auto& constituent : hardestJet.constituents()) {
      registryData.fill(HIST("hPtJetParticle"), constituent.pt());
      double DeltaPhi = TVector2::Phi_0_2pi(constituent.phi() - hardestJet.phi());
      double DeltaEta = constituent.eta() - hardestJet.eta();
      double Delta = TMath::Sqrt(DeltaPhi * DeltaPhi + DeltaEta * DeltaEta);
      registryData.fill(HIST("hJetConeRadius"), Delta);
    }

    fastjet::Selector selector = fastjet::SelectorAbsEtaMax(1.0) * (!fastjet::SelectorNHardest(2)); // TODO: fix subtraction
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
    return hardestJet.constituents();
  }

  template <typename T>
  bool isProton(const T& track)
  {
    bool isProton = false;
    if (track.sign() < 0)
      return isProton;

    // TPC
    if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fProtonTPCnsigLow)
      isProton = true;
    if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fProtonTPCnsigHigh)
      isProton = true;

    registryData.fill(HIST("hDCAzJetProton"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fProtonDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fProtonDCAz)
      return false;

    // TOF
    if (track.pt() < fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fProtonTOFnsigLow)
      isProton = true;
    if (track.pt() > fProtonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fProtonTOFnsigHigh)
      isProton = true;

    return isProton;
  }

  template <typename T>
  bool isAntiproton(const T& track)
  {
    bool isAntiproton = false;
    if (track.sign() < 0)
      return isAntiproton;

    // TPC
    if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fAntiprotonTPCnsigLow)
      isAntiproton = true;
    if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) < fAntiprotonTPCnsigHigh)
      isAntiproton = true;

    registryData.fill(HIST("hDCAzJetAntiproton"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fAntiprotonDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fAntiprotonDCAz)
      return false;

    // TOF
    if (track.pt() < fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fAntiprotonTOFnsigLow)
      isAntiproton = true;
    if (track.pt() > fAntiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) < fAntiprotonTOFnsigHigh)
      isAntiproton = true;

    return isAntiproton;
  }

  template <typename T>
  bool isDeuteron(const T& track)
  {
    bool isDeuteron = false;
    if (track.sign() < 0)
      return isDeuteron;

    // TPC
    if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fDeuteronTPCnsigLow)
      isDeuteron = true;
    if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fDeuteronTPCnsigHigh)
      isDeuteron = true;

    registryData.fill(HIST("hDCAzJetDeuteron"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fDeuteronDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fDeuteronDCAz)
      return false;

    // TOF
    if (track.pt() < fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fDeuteronTOFnsigLow)
      isDeuteron = true;
    if (track.pt() > fDeuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fDeuteronTOFnsigHigh)
      isDeuteron = true;

    return isDeuteron;
  }

  template <typename T>
  bool isAntideuteron(const T& track)
  {
    bool isAntideuteron = false;
    if (track.sign() < 0)
      return isAntideuteron;

    // TPC
    if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fAntideuteronTPCnsigLow)
      isAntideuteron = true;
    if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) < fAntideuteronTPCnsigHigh)
      isAntideuteron = true;

    registryData.fill(HIST("hDCAzJetAntideuteron"), track.pt() * track.sign(), track.dcaZ());

    if (TMath::Abs(track.dcaXY()) > fAntideuteronDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > fAntideuteronDCAz)
      return false;

    // TOF
    if (track.pt() < fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fAntideuteronTOFnsigLow)
      isAntideuteron = true;
    if (track.pt() > fAntideuteronTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) < fAntideuteronTOFnsigHigh)
      isAntideuteron = true;

    return isAntideuteron;
  }

  void setTrackBuffer(const auto& tempBuffer, auto& buffer)
  {
    for (const auto& track : tempBuffer) {
      if (static_cast<int>(buffer.size()) == fTrackBufferSize) {
        buffer.insert(buffer.begin(), track);
        buffer.resize(fTrackBufferSize);
      } else if (static_cast<int>(buffer.size()) < fTrackBufferSize) {
        buffer.emplace_back(track);
      }
    }
  }

  void fillMixedEventDeltas(const auto& track, const auto& buffer, int particleType)
  {
    if (buffer.size() == 0)
      return;
    for (int i = 0; i < static_cast<int>(buffer.size()); i++) {
      if (std::isnan(buffer.at(i).phi()))
        continue;
      if (buffer.at(i).phi() > 2 * TMath::Pi() || buffer.at(i).phi() < -2 * TMath::Pi()) {
        registryData.fill(HIST("hTrackProtocol"), 13);
        continue;
      }

      double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - buffer.at(i).phi());
      if (DeltaPhi > (1.5 * TMath::Pi())) {
        DeltaPhi = DeltaPhi - 2 * TMath::Pi();
      }
      double DeltaEta = TMath::Abs(track.eta() - buffer.at(i).eta());
      registryData.fill(HIST("hDeltaPhiME"), particleType, DeltaPhi);
      registryData.fill(HIST("hDeltaPhiEtaME"), particleType, DeltaPhi, DeltaEta);
    }
  }

  void doCorrelations(const auto& particleVector, const auto& buffer, auto& tempBuffer, int particleType)
  {
    for (int i = 0; i < static_cast<int>(particleVector.size()); i++) { // maybe simply introduce phi cut?
      if (particleVector.at(i).phi() > 2 * TMath::Pi() || particleVector.at(i).phi() < -2 * TMath::Pi()) {
        registryData.fill(HIST("hTrackProtocol"), 11);
        continue;
      }
      for (int j = i + 1; j < static_cast<int>(particleVector.size()); j++) {
        if ((j == static_cast<int>(particleVector.size())) || std::isnan(particleVector.at(j).phi()))
          continue;
        if (particleVector.at(j).phi() > 2 * TMath::Pi() || particleVector.at(j).phi() < -2 * TMath::Pi()) {
          registryData.fill(HIST("hTrackProtocol"), 12);
          continue;
        }

        double DeltaPhi = TVector2::Phi_0_2pi(particleVector.at(i).phi() - particleVector.at(j).phi());
        double DeltaEta = TMath::Abs(particleVector.at(i).eta() - particleVector[j].eta());
        if (DeltaPhi > (1.5 * TMath::Pi())) {
          DeltaPhi = DeltaPhi - 2 * TMath::Pi();
        }
        registryData.fill(HIST("hDeltaPhiSE"), 1, DeltaPhi);
        registryData.fill(HIST("hDeltaPhiEtaSE"), 1, DeltaPhi, DeltaEta);
      }
      fillMixedEventDeltas(particleVector.at(i), buffer, particleType);
      tempBuffer.emplace_back(particleVector.at(i));
    }
  }

  template <typename T, typename U>
  void fillHistogramsRun2(T const& collision, U const& allTracks)
  {
    std::vector<typename FullTracksRun2::iterator> fTempBufferProton;
    std::vector<typename FullTracksRun2::iterator> fTempBufferAntiproton;
    std::vector<typename FullTracksRun2::iterator> fTempBufferDeuteron;
    std::vector<typename FullTracksRun2::iterator> fTempBufferAntideuteron;
    fTempBufferProton.clear();
    fTempBufferAntiproton.clear();
    fTempBufferDeuteron.clear();
    fTempBufferAntideuteron.clear();
    std::vector<fastjet::PseudoJet> jetInput;
    std::map<int, typename U::iterator> particles;
    jetInput.clear();
    particles.clear();
    int index = 0;

    auto tracks = allTracks.sliceBy(perCollisionFullTracksRun2, collision.globalIndex());

    for (const auto& track : tracks) {
      if (!selectTrack(track))
        continue;

      double mass;
      if (track.hasTOF()) {
        mass = track.mass(); // check reliability, maybe use only pion mass
        registryQA.fill(HIST("hTOFmass"), track.pt(), track.mass());
        registryData.fill(HIST("hTrackProtocol"), 4);
      } else {
        mass = 0.139; // pion mass as default, ~80% are pions
        registryData.fill(HIST("hTrackProtocol"), 5);
      }

      // double ratioCrossedRowsTPC = track.tpcNClsCrossedRows()/track.tpcNClsFindable();
      registryQA.fill(HIST("hPtFullEvent"), track.pt());
      registryQA.fill(HIST("hEtaFullEvent"), track.eta());
      registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
      registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
      // registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), ratioCrossedRowsTPC);
      registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
      registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
      registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(mass));
      inputPseudoJet.set_user_index(index);
      particles[index] = track;
      jetInput.emplace_back(inputPseudoJet);

      index++;
    } // for (const auto& track : tracks)

    if (jetInput.size() < 2)
      return;
    registryData.fill(HIST("hEventProtocol"), 2);

    std::vector<fastjet::PseudoJet> constituents = findJets(jetInput);
    if (constituents.empty())
      return;

    std::vector<typename U::iterator> jetProtons; // replace with IDs?
    std::vector<typename U::iterator> jetAntiprotons;
    std::vector<typename U::iterator> jetDeuterons;
    std::vector<typename U::iterator> jetAntideuterons;

    for (int i = 0; i < static_cast<int>(constituents.size()); i++) {
      registryData.fill(HIST("hTrackProtocol"), 10);
      fastjet::PseudoJet pseudoParticle = constituents.at(i);
      int id = pseudoParticle.user_index();
      const auto& jetParticle = particles[id];

      registryData.fill(HIST("hDCAxyFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("hDCAzFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("hTPCsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
      if (jetParticle.hasTOF())
        registryData.fill(HIST("hTOFsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());

      double ptDiff = pseudoParticle.pt() - jetParticle.pt();
      registryData.fill(HIST("hPtDiff"), ptDiff);
      int particleType = 0;

      if (jetParticle.pt() < fMinJetParticlePt)
        continue;
      if (isProton(jetParticle) || isAntiproton(jetParticle)) { // collect (anti)protons in jet
        registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());
        if (isProton(jetParticle)) {
          particleType = 1;
          registryData.fill(HIST("hTrackProtocol"), 6); // # protons
          jetProtons.emplace_back(jetParticle);
        } else {
          particleType = 2;
          registryData.fill(HIST("hTrackProtocol"), 7); // # antiprotons
          jetAntiprotons.emplace_back(jetParticle);
        }
      } else if (isDeuteron(jetParticle) || isAntideuteron(jetParticle)) { // collect (anti)deuterons in jet
        registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaDe());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaDe());
        if (isDeuteron(jetParticle)) {
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
    } // for (int i=0; i<static_cast<int>(constituents.size()); i++)

    if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2) && (jetDeuterons.size() < 2) && (jetAntideuterons.size() < 2))
      return;
    registryData.fill(HIST("hEventProtocol"), 6);

    if (jetProtons.size() > 1) {
      doCorrelations(jetProtons, fTrackBufferProtonRun2, fTempBufferProton, 1);
      setTrackBuffer(fTempBufferProton, fTrackBufferProtonRun2);
    }
    if (jetAntiprotons.size() > 1) {
      doCorrelations(jetAntiprotons, fTrackBufferAntiprotonRun2, fTempBufferAntiproton, 2);
      setTrackBuffer(fTempBufferAntiproton, fTrackBufferAntiprotonRun2);
    }
    if (jetDeuterons.size() > 1) {
      doCorrelations(jetDeuterons, fTrackBufferDeuteronRun2, fTempBufferDeuteron, 3);
      setTrackBuffer(fTempBufferDeuteron, fTrackBufferDeuteronRun2);
    }
    if (jetAntideuterons.size() > 1) {
      doCorrelations(jetAntideuterons, fTrackBufferAntideuteronRun2, fTempBufferAntideuteron, 4);
      setTrackBuffer(fTempBufferAntideuteron, fTrackBufferAntideuteronRun2);
    }
  }

  template <typename T, typename U>
  void fillHistogramsRun3(T const& collision, U const& allTracks)
  {
    std::vector<typename FullTracksRun3::iterator> fTempBufferProton;
    std::vector<typename FullTracksRun3::iterator> fTempBufferAntiproton;
    std::vector<typename FullTracksRun3::iterator> fTempBufferDeuteron;
    std::vector<typename FullTracksRun3::iterator> fTempBufferAntideuteron;
    fTempBufferProton.clear();
    fTempBufferAntiproton.clear();
    fTempBufferDeuteron.clear();
    fTempBufferAntideuteron.clear();
    std::vector<fastjet::PseudoJet> jetInput;
    std::map<int, typename U::iterator> particles;
    jetInput.clear();
    particles.clear();
    int index = 0;

    auto tracks = allTracks.sliceBy(perCollisionFullTracksRun2, collision.globalIndex());

    for (const auto& track : tracks) {
      if (!selectTrack(track))
        continue;

      double mass;
      if (track.hasTOF()) {
        mass = track.mass(); // check reliability, maybe use only pion mass
        registryQA.fill(HIST("hTOFmass"), track.pt(), track.mass());
        registryData.fill(HIST("hTrackProtocol"), 4);
      } else {
        mass = 0.139; // pion mass as default, ~80% are pions
        registryData.fill(HIST("hTrackProtocol"), 5);
      }

      // double ratioCrossedRowsTPC = track.tpcNClsCrossedRows()/track.tpcNClsFindable();
      registryQA.fill(HIST("hPtFullEvent"), track.pt());
      registryQA.fill(HIST("hEtaFullEvent"), track.eta());
      registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
      registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
      // registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), ratioCrossedRowsTPC);
      registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
      registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
      registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(mass));
      inputPseudoJet.set_user_index(index);
      particles[index] = track;
      jetInput.emplace_back(inputPseudoJet);

      index++;
    } // for (const auto& track : tracks)

    if (jetInput.size() < 2)
      return;
    registryData.fill(HIST("hEventProtocol"), 2);

    std::vector<fastjet::PseudoJet> constituents = findJets(jetInput);
    if (constituents.empty())
      return;

    std::vector<typename U::iterator> jetProtons; // replace with IDs?
    std::vector<typename U::iterator> jetAntiprotons;
    std::vector<typename U::iterator> jetDeuterons;
    std::vector<typename U::iterator> jetAntideuterons;

    for (int i = 0; i < static_cast<int>(constituents.size()); i++) {
      registryData.fill(HIST("hTrackProtocol"), 10);
      fastjet::PseudoJet pseudoParticle = constituents.at(i);
      int id = pseudoParticle.user_index();
      const auto& jetParticle = particles[id];

      registryData.fill(HIST("hDCAxyFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("hDCAzFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("hTPCsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
      if (jetParticle.hasTOF())
        registryData.fill(HIST("hTOFsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());

      double ptDiff = pseudoParticle.pt() - jetParticle.pt();
      registryData.fill(HIST("hPtDiff"), ptDiff);
      int particleType = 0;

      if (jetParticle.pt() < fMinJetParticlePt)
        continue;
      if (isProton(jetParticle) || isAntiproton(jetParticle)) { // collect (anti)protons in jet
        registryData.fill(HIST("hTPCnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaPr());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaProton"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaPr());
        if (isProton(jetParticle)) {
          particleType = 1;
          registryData.fill(HIST("hTrackProtocol"), 6); // # protons
          jetProtons.emplace_back(jetParticle);
        } else {
          particleType = 2;
          registryData.fill(HIST("hTrackProtocol"), 7); // # antiprotons
          jetAntiprotons.emplace_back(jetParticle);
        }
      } else if (isDeuteron(jetParticle) || isAntideuteron(jetParticle)) { // collect (anti)deuterons in jet
        registryData.fill(HIST("hTPCnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcNSigmaDe());
        if (jetParticle.hasTOF())
          registryData.fill(HIST("hTOFnsigmaDeuteron"), jetParticle.pt() * jetParticle.sign(), jetParticle.tofNSigmaDe());
        if (isDeuteron(jetParticle)) {
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
    } // for (int i=0; i<static_cast<int>(constituents.size()); i++)

    if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2) && (jetDeuterons.size() < 2) && (jetAntideuterons.size() < 2))
      return;
    registryData.fill(HIST("hEventProtocol"), 6);

    if (jetProtons.size() > 1) {
      doCorrelations(jetProtons, fTrackBufferProtonRun3, fTempBufferProton, 1);
      setTrackBuffer(fTempBufferProton, fTrackBufferProtonRun3);
    }
    if (jetAntiprotons.size() > 1) {
      doCorrelations(jetAntiprotons, fTrackBufferAntiprotonRun3, fTempBufferAntiproton, 2);
      setTrackBuffer(fTempBufferAntiproton, fTrackBufferAntiprotonRun3);
    }
    if (jetDeuterons.size() > 1) {
      doCorrelations(jetDeuterons, fTrackBufferDeuteronRun3, fTempBufferDeuteron, 3);
      setTrackBuffer(fTempBufferDeuteron, fTrackBufferDeuteronRun3);
    }
    if (jetAntideuterons.size() > 1) {
      doCorrelations(jetAntideuterons, fTrackBufferAntideuteronRun3, fTempBufferAntideuteron, 4);
      setTrackBuffer(fTempBufferAntideuteron, fTrackBufferAntideuteronRun3);
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

      fillHistogramsRun2(collision, tracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun2, "process Run 2 data", true);

  void processRun3(aod::Collision const& collision, soa::Filtered<FullTracksRun3> const& tracks)
  {
    fillHistogramsRun3(collision, tracks);
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun3, "process Run 3 data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc)};
}
