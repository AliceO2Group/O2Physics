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
/// \file angularCorrelationsInJets.cxx
///
/// \author Lars Jörgensen (lars.christian.joergensen@cern.ch)
/// \brief task for analysis of angular correlations in jets using Fastjet

#include <vector>
#include <utility>
#include <map>
#include <string>

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
#include "Common/Core/RecoDecay.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "TVector3.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct AxisSpecs {
  AxisSpec ptAxisPos = {1000, 0, 100, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec ptAxisFull = {2000, -100, 100, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec nsigmapTAxis = {500, 0, 50, "#it{p}_{T} [GeV/#it{c}]"};
  AxisSpec nsigmaAxis = {300, -15, 15, "n#sigma"};
  AxisSpec dcazAxis = {1000, -1, 1, "DCA_{z} [cm]"};
  AxisSpec dcaxyAxis = {1000, -0.5, 0.5, "DCA_{xy} [cm]"};
  AxisSpec angDistPhiAxis = {1000, -2, 5, "#Delta#varphi"};
  AxisSpec angDistEtaAxis = {1000, -2, 2, "#Delta#eta"};
};

struct AngularCorrelationsInJets {
  // Switches
  Configurable<bool> useRejectionCut{"useRejectionCut", true, "use nsigmaRejection for correlations"};
  Configurable<bool> outputQC{"outputQC", true, "add QC output"};
  Configurable<bool> doppCorrelations{"doppCorrelations", true, "measure correlations for p-p"};
  Configurable<bool> doapapCorrelations{"doapapCorrelations", false, "measure correlations for pbar-pbar"};
  Configurable<bool> dopapCorrelations{"dopapCorrelations", false, "measure correlations for p-pbar"};
  Configurable<bool> dopipiCorrelations{"dopipiCorrelations", false, "measure correlations for pi+-p+, pi--pi-"};
  Configurable<bool> doJetCorrelations{"doJetCorrelations", false, "measure correlations for all particles inside jets"};
  Configurable<bool> doFullCorrelations{"doFullCorrelations", false, "measure correlations for all particles in an event"};
  Configurable<bool> measureKaons{"measureKaons", false, "measure correlations for K-K"};

  // Track Cuts
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80, "min number of crossed rows TPC"};
  Configurable<int> minReqClusterITS{"minReqClusterITS", 5, "min number of clusters required in ITS"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.05, "max DCA to vertex xy"};
  Configurable<float> maxDCAz{"maxDCAz", 0.05, "max DCA to vertex z"};
  Configurable<float> maxEta{"maxEta", 0.8, "max pseudorapidity"};
  Configurable<float> deltaEtaEdge{"deltaEtaEdge", 0.05, "min eta distance of jet from acceptance edge"};
  Configurable<float> minTrackPt{"minTrackPt", 0.3, "minimum track pT"};
  Configurable<float> requirePVContributor{"requirePVContributor", false, "require track to be PV contributor"};

  // Jet Cuts
  Configurable<float> jetR{"jetR", 0.3, "jet resolution parameter"};
  Configurable<float> minJetPt{"minJetPt", 10.0, "minimum total pT to accept jet"};

  // Proton Cuts
  Configurable<float> protonDCAxyCF{"protonDCAxyCF", 0.05, "[proton] DCAxy cut for CF"};
  Configurable<float> protonDCAzCF{"protonDCAzCF", 0.02, "[proton] DCAz cut for CF"};
  Configurable<float> protonTPCTOFpT{"protonTPCTOFpT", 0.7, "[proton] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> protonTPCnsigma{"protonTPCnsigma", 4.0, "[proton] max TPC nsigma for pt > 0/1.5/3.0 GeV"};
  Configurable<float> protonTOFnsigma{"protonTOFnsigma", 3.0, "[proton] max TOF nsigma for pt > 0/1.5/3.0 GeV"};

  // Antiproton Cuts
  Configurable<float> antiprotonDCAxyCF{"antiprotonDCAxyCF", 0.05, "[antiproton] DCAxy cut for CF"};
  Configurable<float> antiprotonDCAzCF{"antiprotonDCAzCF", 0.02, "[antiproton] DCAz cut for CF"};
  Configurable<float> antiprotonTPCTOFpT{"antiprotonTPCTOFpT", 0.7, "[antiproton] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> antiprotonTPCnsigma{"antiprotonTPCnsigma", 4.0, "[antiproton] max TPC nsigma for pt > 0/1.5/3.0 GeV"};
  Configurable<float> antiprotonTOFnsigma{"antiprotonTOFnsigma", 3.0, "[antiproton] max TOF nsigma for pt > 0/1.5/3.0 GeV"};

  // Pion & Kaon PID
  Configurable<float> pionDCAxy{"pionDCAxy", 0.05, "[pion] DCAxy cut"};
  Configurable<float> pionDCAz{"pionDCAz", 0.05, "[pion] DCAz cut"};
  Configurable<float> pionTPCTOFpT{"pionTPCTOFpT", 0.7, "[pion] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> pionTPCnsigmaLowPt{"pionTPCnsigmaLowPt", 3.0, "[pion] max TPC nsigma with low pT"};
  Configurable<float> pionTPCnsigmaHighPt{"pionTPCnsigmaHighPt", 3.0, "[pion] max TPC nsigma with high pT"};
  Configurable<float> pionTOFnsigma{"pionTOFnsigma", 3.0, "[pion] max TOF nsigma"};
  Configurable<float> kaonDCAxy{"kaonDCAxy", 0.05, "[kaon] DCAxy cut"};
  Configurable<float> kaonDCAz{"kaonDCAz", 0.05, "[kaon] DCAz cut"};
  Configurable<float> kaonTPCTOFpT{"kaonTPCTOFpT", 0.7, "[kaon] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> kaonTPCnsigmaLowPt{"kaonTPCnsigmaLowPt", 3.0, "[kaon] max TPC nsigma with low pT"};
  Configurable<float> kaonTPCnsigmaHighPt{"kaonTPCnsigmaHighPt", 3.0, "[kaon] max TPC nsigma with high pT"};
  Configurable<float> kaonTOFnsigma{"kaonTOFnsigma", 3.0, "[kaon] max TOF nsigma"};

  Configurable<float> nsigmaRejection{"nsigmaRejection", 1.0, "reject tracks with nsigma < nsigmaRejection for >1 species"};
  Configurable<int> trackBufferSize{"trackBufferSize", 200, "Number of mixed-event tracks being stored"};

  // QC Configurables
  Configurable<float> zVtx{"zVtx", 10.0, "max zVertex"};
  Configurable<float> rMax{"rMax", 0.3, "Maximum radius for jet and UE regions"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  using FullTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                   aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa>;
  using FullTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                   aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa>;
  using McTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                 aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa, aod::McTrackLabels>;
  using McTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                 aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa, aod::McTrackLabels>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  using McCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Filter prelimTrackCuts = (aod::track::itsChi2NCl < maxChi2ITS &&
                            aod::track::tpcChi2NCl < maxChi2TPC &&
                            nabs(aod::track::dcaXY) < maxDCAxy);

  Preslice<FullTracksRun2> perCollisionFullTracksRun2 = o2::aod::track::collisionId;
  Preslice<FullTracksRun3> perCollisionFullTracksRun3 = o2::aod::track::collisionId;
  Preslice<McTracksRun2> perCollisionMcTracksRun2 = o2::aod::track::collisionId;
  Preslice<McTracksRun3> perCollisionMcTracksRun3 = o2::aod::track::collisionId;

  AxisSpecs axisSpecs;

  HistogramRegistry registryData{"data", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMC{"MC", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryQC{"QC", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  JetBkgSubUtils bkgSub;
  std::vector<int> eventSelection;

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    // Counters
    registryData.add("numberOfEvents", "Number of events", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("numberOfJets", "Total number of jets", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("eventProtocol", "Event protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("trackProtocol", "Track protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("numPartInJet", "Number of particles in a jet", HistType::kTH1I, {{200, 0, 200}});
    registryData.add("numJetsInEvent", "Number of jets selected", HistType::kTH1I, {{10, 0, 10}});
    registryData.add("jetRapidity", "Jet rapidity;#it{y}", HistType::kTH1F, {{200, -1, 1}});

    registryData.add("ptJetProton", "p_{T} of protons", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("ptJetAntiproton", "p_{T} of antiprotons", HistType::kTH1D, {axisSpecs.ptAxisPos});

    registryQC.add("ptJetProtonVsTotalJet", "Proton p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQC.add("ptJetAntiprotonVsTotalJet", "Antiproton p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});

    registryData.add("tpcNSigmaProtonCF", "TPC n#sigma for proton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tofNSigmaProtonCF", "TOF n#sigma for proton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tpcNSigmaAntiprotonCF", "TPC n#sigma for antiproton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tofNSigmaAntiprotonCF", "TOF n#sigma for antiproton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tpcNSigmaPion", "TPC n#sigma for pion", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tofNSigmaPion", "TOF n#sigma for pion", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});

    registryData.add("deltaPhiSEFull", "#Delta#varphi of particles in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiSEJet", "#Delta#varphi of jet particles in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiSEProton", "#Delta#varphi of protons in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiSEAntiproton", "#Delta#varphi of antiprotons in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiSEPion", "#Delta#varphi of pions in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiSEProtonAntiproton", "#Delta#varphi of proton-antiproton in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});

    registryData.add("deltaPhiMEFull", "#Delta#varphi of particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiMEJet", "#Delta#varphi of jet particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiMEProton", "#Delta#varphi of protons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiMEAntiproton", "#Delta#varphi of antiprotons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiMEPion", "#Delta#varphi of pions in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("deltaPhiMEProtonAntiproton", "#Delta#varphi of proton-antiproton in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});

    registryData.add("deltaPhiEtaSEFull", "#Delta#varphi vs #Delta#eta of full particles in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaSEJet", "#Delta#varphi vs #Delta#eta of jet particles in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaSEProton", "#Delta#varphi vs #Delta#eta of protons in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaSEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaSEPion", "#Delta#varphi vs #Delta#eta of pions in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaSEProtonAntiproton", "#Delta#varphi vs #Delta#eta of proton-antiproton in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});

    registryData.add("deltaPhiEtaMEFull", "#Delta#varphi vs #Delta#eta of particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaMEJet", "#Delta#varphi vs #Delta#eta of jet particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaMEProton", "#Delta#varphi vs #Delta#eta of protons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaMEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaMEPion", "#Delta#varphi vs #Delta#eta of pions in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("deltaPhiEtaMEProtonAntiproton", "#Delta#varphi vs #Delta#eta of proton-antiproton in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});

    registryData.add("dcaZJetProton", "DCA_{z} of high purity protons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("dcaZJetAntiproton", "DCA_{z} of high purity antiprotons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("dcaZJetPion", "DCA_{z} of high purity pions", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    // pT
    registryData.add("ptJetParticle", "p_{T} of particles in jets", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("ptTotalSubJetArea", "Subtracted full jet p_{T} (area)", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("ptTotalSubJetPerp", "Subtracted full jet p_{T} (perpendicular)", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("ptTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{1000, 0, 500}});
    registryQC.add("ptJetPionVsTotalJet", "Pion p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQC.add("ptJetKaonVsTotalJet", "Kaon p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});

    // nSigma
    registryData.add("tpcSignal", "TPC signal", HistType::kTH2F, {{400, -10, 10, "#it{p} [GeV/#it{c}]"}, {1000, 0, 500, "d#it{E}/d#it{X} (a.u.)"}});
    registryData.add("tofSignal", "TOF signal", HistType::kTH2F, {{400, -20, 20, "#it{p} [GeV/#it{c}]"}, {400, 0.75, 1.05, "#beta (TOF)"}});
    registryData.add("tpcNSigmaKaon", "TPC n#sigma for kaon", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("tofNSigmaKaon", "TOF n#sigma for kaon", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});

    // DCA
    registryData.add("dcaXYFullJet", "DCA_{xy} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcaxyAxis});
    registryData.add("dcaZFullJet", "DCA_{z} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcazAxis});
    registryData.add("dcaZJetKaon", "DCA_{z} of high purity kaons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});

    // Angular Distributions
    registryQC.add("phiFullEvent", "#varphi in full event", HistType::kTH1F, {{1000, 0, 6.3}});
    registryQC.add("phiPtFullEvent", "#varphi vs. p_{T} in full event", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 6.3}});
    registryQC.add("phiJet", "#varphi in jet", HistType::kTH1F, {{1000, 0, 6.3}});
    registryQC.add("phiPtJet", "#varphi vs. p_{T} in jet", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 6.3}});
    registryQC.add("etaFullEvent", "#eta in full event", HistType::kTH1F, {{1000, -1, 1}});
    registryQC.add("etaPtFullEvent", "#eta vs. p_{T} in full event", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, -1, 1}});
    registryQC.add("etaJet", "#eta in jet", HistType::kTH1F, {{1000, -1, 1}});
    registryQC.add("etaPtJet", "#eta vs. p_{T} in jet", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, -1, 1}});

    // QA
    registryQC.add("rhoEstimatePerp", "Background #rho (perp)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQC.add("rhoMEstimatePerp", "Background #rho_{m} (perp)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQC.add("rhoEstimateArea", "Background #rho (area)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQC.add("rhoMEstimateArea", "Background #rho_{m} (area)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQC.add("jetPtVsNumPart", "Total jet p_{T} vs number of constituents", HistType::kTH2F, {axisSpecs.ptAxisPos, {100, 0, 100}});

    if (doprocessMCRun2 || doprocessMCRun3) {
      registryMC.add("ptJetProtonMC", "Truth jet proton p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
      registryMC.add("ptJetAntiprotonMC", "Truth jet antiproton p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
      registryMC.add("numberOfTruthParticles", "Truth yields (anti)p, (anti)d, (anti)He-3", HistType::kTH1I, {{6, 0, 6}});
    }

    if (outputQC) {
      registryQC.add("ptDiff", "p_{T} difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{2000, -0.0001, 0.0001}});
      registryQC.add("jetConeRadius", "Jet Radius;#it{R}", HistType::kTH1F, {{100, 0, 1}});
      registryQC.add("maxRadiusVsPt", "Max Cone Radius vs p_{T}", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {100, 0, 1}});
      registryQC.add("jetBkgDeltaPt", "#Delta p_{T} Clustered Cone - Pure Jet", HistType::kTH1F, {{200, 0, 10}});

      registryQC.add("ptFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {axisSpecs.ptAxisPos});
      registryQC.add("crossedRowsTPC", "Crossed rows TPC", HistType::kTH2I, {axisSpecs.ptAxisPos, {135, 65, 200}});
      registryQC.add("clusterITS", "ITS clusters", HistType::kTH2I, {axisSpecs.ptAxisPos, {10, 0, 10}});
      registryQC.add("clusterTPC", "TPC clusters", HistType::kTH2I, {axisSpecs.ptAxisPos, {135, 65, 200}});
      registryQC.add("ratioCrossedRowsTPC", "Ratio crossed rows/findable TPC", HistType::kTH2F, {axisSpecs.ptAxisPos, {100, 0.5, 1.5}});
      registryQC.add("chi2ITS", "ITS #chi^{2}", HistType::kTH2F, {axisSpecs.ptAxisPos, {400, 0, 40}});
      registryQC.add("chi2TPC", "TPC #chi^{2}", HistType::kTH2F, {axisSpecs.ptAxisPos, {50, 0, 5}});
      registryQC.add("dcaXYFullEvent", "DCA_{xy} of full event", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcaxyAxis});
      registryQC.add("dcaZFullEvent", "DCA_{z} of full event", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});

      registryQC.add("multiplicityJetPlusUE", "multiplicityJetPlusUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("multiplicityJet", "multiplicityJet", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("multiplicityUE", "multiplicityUE", HistType::kTH1F, {{100, 0, 100, "#it{N}_{ch}"}});
      registryQC.add("ptJetPlusUE", "ptJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("ptJet", "ptJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("ptUE", "ptUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("deltaEtadeltaPhiJet", "deltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, constants::math::PIHalf, "#Delta#phi"}});
      registryQC.add("deltaEtadeltaPhiUE", "deltaEtadeltaPhiUE", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, constants::math::PIHalf, "#Delta#phi"}});
      registryQC.add("deltaJetPt", "deltaJetPt", HistType::kTH1F, {{200, -2, 2, "#Delta#it{p}_{T} (GeV/#it{c})"}});
      registryQC.add("nParticlesClusteredInJet", "nParticlesClusteredInJet", HistType::kTH1F, {{50, 0, 50, "#it{N}_{ch}"}});
      registryQC.add("ptParticlesClusteredInJet", "ptParticlesClusteredInJet", HistType::kTH1F, {{200, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});
    }
  }

  std::vector<std::pair<double, double>> fBufferProton;
  std::vector<std::pair<double, double>> fBufferAntiproton;
  std::vector<std::pair<double, double>> fBufferPiPlus;
  std::vector<std::pair<double, double>> fBufferPiMinus;
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

  template <typename T>
  bool hasITSHit(const T& track, int layer)
  {
    int bit = layer - 1;
    return (track.itsClusterMap() & (1 << bit));
  }

  template <typename T>
  bool selectTrackForJetReco(const T& track)
  {
    if (track.dcaZ() > maxDCAz)
      return false;
    double maxEtaForJetReco = 0.8;
    if (track.eta() > maxEtaForJetReco)
      return false;
    double minTrackPtForJetReco = 0.1;
    if (track.pt() < minTrackPtForJetReco)
      return false;
    if (!track.hasITS())
      return false;
    if (!track.hasTPC())
      return false;
    int minCrossedRowsForJetReco = 70;
    if (track.tpcNClsCrossedRows() < minCrossedRowsForJetReco)
      return false;
    if ((!hasITSHit(track, 1)) && (!hasITSHit(track, 2)) && (!hasITSHit(track, 3)))
      return false;
    double minRatioCrRowsFindableJetReco = 0.8;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < minRatioCrRowsFindableJetReco)
      return false;
    if (std::fabs(track.dcaXY()) > (0.0105 + 0.035 / std::pow(track.pt(), 1.1)))
      return false;
    if (doprocessRun2 || doprocessMCRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit) ||
          !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selectTrack(const T& track)
  {
    if (requirePVContributor && !(track.isPVContributor()))
      return false;
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < minReqClusterITS)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if ((static_cast<double>(track.tpcNClsCrossedRows()) / static_cast<double>(track.tpcNClsFindable())) < minRatioCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.itsChi2NCl() > maxChi2ITS)
      return false;
    if (std::fabs(track.eta()) > maxEta)
      return false;
    if (track.pt() < minTrackPt)
      return false;

    return true;
  }

  template <class T>
  bool singleSpeciesTPCNSigma(T const& track) // reject any track that has TPC nsigma < 3 for more than 1 species
  {
    if (useRejectionCut && (track.tpcNSigmaStoreEl() < nsigmaRejection || track.tpcNSigmaStoreMu() < nsigmaRejection || track.tpcNSigmaStorePi() < nsigmaRejection || track.tpcNSigmaStoreKa() < nsigmaRejection || track.tpcNSigmaStoreTr() < nsigmaRejection || track.tpcNSigmaStoreAl() < nsigmaRejection || track.tpcNSigmaDe() < nsigmaRejection || track.tpcNSigmaHe() < nsigmaRejection))
      return false;
    return true;
  }

  template <typename T>
  bool isProton(const T& track)
  {
    if (track.sign() < 0)
      return false;

    double pt = track.pt();

    // DCA
    double maxDCApt = 1.2;
    if (pt < maxDCApt) {
      if (std::abs(track.dcaXY()) > protonDCAxyCF)
        return false;
      if (std::abs(track.dcaZ()) > protonDCAzCF)
        return false;
    }

    // nsigma
    double midPt = 1.5;
    double highPt = 3.0;

    double maxTPCnsigma = protonTPCnsigma;
    double maxTOFnsigma = protonTOFnsigma;
    if (pt > midPt) {
      maxTPCnsigma = protonTPCnsigma - 1;
      maxTOFnsigma = protonTOFnsigma - 1;
    }
    if (pt > highPt) {
      maxTPCnsigma = protonTPCnsigma - 2;
      maxTOFnsigma = protonTOFnsigma - 2;
    }

    registryData.fill(HIST("tpcNSigmaProtonCF"), track.pt(), track.tpcNSigmaPr());
    if (pt < protonTPCTOFpT && (std::abs(track.tpcNSigmaPr()) > maxTPCnsigma))
      return false;

    double tofNSigma = 999;
    if (track.hasTOF()) {
      registryData.fill(HIST("tofNSigmaProtonCF"), track.pt(), track.tofNSigmaPr());
      tofNSigma = track.tofNSigmaPr();
    }

    if (pt > protonTPCTOFpT && ((std::abs(tofNSigma) > maxTOFnsigma) || std::abs(track.tpcNSigmaPr()) > maxTPCnsigma))
      return false;

    if (useRejectionCut && !singleSpeciesTPCNSigma(track))
      return false;

    return true;
  }

  template <typename T>
  bool isAntiproton(const T& track)
  {
    if (track.sign() > 0)
      return false;

    double pt = track.pt();

    // DCA
    double maxDCApt = 1.2;
    if (pt < maxDCApt) {
      if (std::abs(track.dcaXY()) > antiprotonDCAxyCF)
        return false;
      if (std::abs(track.dcaZ()) > antiprotonDCAzCF)
        return false;
    }

    // nsigma
    double midPt = 1.5;
    double highPt = 3.0;

    double maxTPCnsigma = antiprotonTPCnsigma;
    double maxTOFnsigma = antiprotonTOFnsigma;
    if (pt > midPt) {
      maxTPCnsigma = antiprotonTPCnsigma - 1;
      maxTOFnsigma = antiprotonTOFnsigma - 1;
    }
    if (pt > highPt) {
      maxTPCnsigma = antiprotonTPCnsigma - 2;
      maxTOFnsigma = antiprotonTOFnsigma - 2;
    }

    registryData.fill(HIST("tpcNSigmaAntiprotonCF"), track.pt(), track.tpcNSigmaPr());
    if (pt < antiprotonTPCTOFpT && (std::abs(track.tpcNSigmaPr()) > maxTPCnsigma))
      return false;

    double tofNSigma = 999;
    if (track.hasTOF()) {
      registryData.fill(HIST("tofNSigmaAntiprotonCF"), track.pt(), track.tofNSigmaPr());
      tofNSigma = track.tofNSigmaPr();
    }

    if (pt > antiprotonTPCTOFpT && ((std::abs(tofNSigma) > maxTOFnsigma) || std::abs(track.tpcNSigmaPr()) > maxTPCnsigma))
      return false;

    if (useRejectionCut && !singleSpeciesTPCNSigma(track))
      return false;

    return true;
  }

  template <typename T>
  bool isPion(const T& track)
  {
    // DCA
    if (std::abs(track.dcaXY()) > pionDCAxy)
      return false;
    if (std::abs(track.dcaZ()) > pionDCAz)
      return false;

    registryData.fill(HIST("tpcNSigmaPion"), track.pt(), track.tpcNSigmaStorePi());

    // TPC
    if (track.pt() < pionTPCTOFpT && std::abs(track.tpcNSigmaStorePi()) > pionTPCnsigmaLowPt)
      return false;
    if (track.pt() > pionTPCTOFpT && std::abs(track.tpcNSigmaStorePi()) > pionTPCnsigmaHighPt)
      return false;

    // TOF
    if (track.hasTOF()) {
      registryData.fill(HIST("tofNSigmaPion"), track.pt(), track.tofNSigmaStorePi());
      if (track.pt() > pionTPCTOFpT && std::abs(track.tofNSigmaStorePi()) > pionTOFnsigma)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isKaon(const T& track)
  {
    // DCA
    if (std::abs(track.dcaXY()) > kaonDCAxy)
      return false;
    if (std::abs(track.dcaZ()) > kaonDCAz)
      return false;

    registryData.fill(HIST("tpcNSigmaKaon"), track.pt(), track.tpcNSigmaStoreKa());

    // TPC
    if (track.pt() < kaonTPCTOFpT && std::abs(track.tpcNSigmaStoreKa()) > kaonTPCnsigmaLowPt)
      return false;
    if (track.pt() > kaonTPCTOFpT && std::abs(track.tpcNSigmaStoreKa()) > kaonTPCnsigmaHighPt)
      return false;

    // TOF
    if (track.hasTOF()) {
      registryData.fill(HIST("tofNSigmaKaon"), track.pt(), track.tofNSigmaStoreKa());
      if (track.pt() > kaonTPCTOFpT && std::abs(track.tofNSigmaStoreKa()) > kaonTOFnsigma)
        return false;
    }

    return true;
  }

  void setTrackBuffer(const auto& tempBuffer, auto& buffer) // refresh track buffer
  {
    for (const auto& pair : tempBuffer) {
      if (static_cast<int>(buffer.size()) == trackBufferSize) {
        buffer.insert(buffer.begin(), pair);
        buffer.resize(trackBufferSize);
      } else if (static_cast<int>(buffer.size()) < trackBufferSize) {
        buffer.emplace_back(pair);
      }
    }
  }

  void fillMixedEventDeltas(const auto& track, const auto& buffer, int particleType, const TVector3 jetAxis) // correlate tracks from current event with tracks from buffer, i.e. other events
  {
    if (buffer.size() == 0)
      return;
    if (std::isnan(track.phi()) || std::isnan(jetAxis.Phi()))
      return;
    for (int i = 0; i < static_cast<int>(buffer.size()); i++) { // loop over tracks in buffer
      if (std::isnan(buffer.at(i).first))
        continue;
      if (buffer.at(i).first > constants::math::TwoPI || buffer.at(i).first < constants::math::TwoPI) {
        registryData.fill(HIST("trackProtocol"), 13); // # buffer tracks failed with phi > 2 pi
        continue;
      }

      double phiToAxis = RecoDecay::constrainAngle(track.phi() - jetAxis.Phi(), 0);
      double etaToAxis = track.eta() - jetAxis.Eta();
      double deltaPhi = RecoDecay::constrainAngle(phiToAxis - buffer.at(i).first, -constants::math::PIHalf);
      double deltaEta = etaToAxis - buffer.at(i).second;

      switch (particleType) {
        case -1:
          registryData.fill(HIST("deltaPhiMEFull"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEFull"), deltaPhi, deltaEta);
          break;
        case 0:
          registryData.fill(HIST("deltaPhiMEJet"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEJet"), deltaPhi, deltaEta);
          break;
        case 1:
          registryData.fill(HIST("deltaPhiMEProton"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEProton"), deltaPhi, deltaEta);
          break;
        case 2:
          registryData.fill(HIST("deltaPhiMEAntiproton"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEAntiproton"), deltaPhi, deltaEta);
          break;
        case 3:
          registryData.fill(HIST("deltaPhiMEPion"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEPion"), deltaPhi, deltaEta);
          break;
        case 4:
          registryData.fill(HIST("deltaPhiMEProtonAntiproton"), deltaPhi);
          registryData.fill(HIST("deltaPhiEtaMEProtonAntiproton"), deltaPhi, deltaEta);
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
      double phiToAxis = RecoDecay::constrainAngle(particleVector.at(i).phi() - jetAxis.Phi(), 0);
      double etaToAxis = particleVector.at(i).eta() - jetAxis.Eta();
      if (std::abs(particleVector.at(i).phi()) > constants::math::TwoPI) {
        registryData.fill(HIST("trackProtocol"), 11); // # tracks failed with phi > 2 pi
        continue;
      }
      for (int j = i + 1; j < static_cast<int>(particleVector.size()); j++) {
        if ((j == static_cast<int>(particleVector.size())) || std::isnan(particleVector.at(j).phi()))
          continue;
        if (std::abs(particleVector.at(j).phi()) > constants::math::TwoPI) {
          registryData.fill(HIST("trackProtocol"), 12); // # tracks failed with phi > 2 pi
          continue;
        }

        double deltaPhi = RecoDecay::constrainAngle(particleVector.at(i).phi() - particleVector.at(j).phi(), -constants::math::PIHalf);
        double deltaEta = particleVector.at(i).eta() - particleVector.at(j).eta();
        switch (particleType) {
          case -1:
            registryData.fill(HIST("deltaPhiSEFull"), deltaPhi);
            registryData.fill(HIST("deltaPhiEtaSEFull"), deltaPhi, deltaEta);
            break;
          case 0:
            registryData.fill(HIST("deltaPhiSEJet"), deltaPhi);
            registryData.fill(HIST("deltaPhiEtaSEJet"), deltaPhi, deltaEta);
            break;
          case 1:
            registryData.fill(HIST("deltaPhiSEProton"), deltaPhi);
            registryData.fill(HIST("deltaPhiEtaSEProton"), deltaPhi, deltaEta);
            break;
          case 2:
            registryData.fill(HIST("deltaPhiSEAntiproton"), deltaPhi);
            registryData.fill(HIST("deltaPhiEtaSEAntiproton"), deltaPhi, deltaEta);
            break;
          case 3:
            registryData.fill(HIST("deltaPhiSEPion"), deltaPhi);
            registryData.fill(HIST("deltaPhiEtaSEPion"), deltaPhi, deltaEta);
            break;
        }
      }
      fillMixedEventDeltas(particleVector.at(i), buffer, particleType, jetAxis);
      tempBuffer.emplace_back(std::make_pair(phiToAxis, etaToAxis));
    }
  }

  void doCorrelationsAnti(const auto& particleVector, const auto& particleVectorAnti, const auto& bufferAnti, auto& tempBuffer, const TVector3 jetAxis) // correlations between particle/antiparticle
  {
    if (std::isnan(jetAxis.Phi()))
      return;
    for (int i = 0; i < static_cast<int>(particleVector.size()); i++) {
      if (std::isnan(particleVector.at(i).phi()))
        continue;
      double phiToAxis = RecoDecay::constrainAngle(particleVector.at(i).phi() - jetAxis.Phi(), 0);
      double etaToAxis = particleVector.at(i).eta() - jetAxis.Eta();
      if (std::abs(particleVector.at(i).phi()) > constants::math::TwoPI) {
        registryData.fill(HIST("trackProtocol"), 14); // # tracks failed with phi > 2 pi
        continue;
      }
      for (int j = 0; j < static_cast<int>(particleVectorAnti.size()); j++) {
        if (std::isnan(particleVectorAnti.at(j).phi()))
          continue;
        if (std::abs(particleVectorAnti.at(j).phi()) > constants::math::TwoPI) {
          registryData.fill(HIST("trackProtocol"), 15); // # tracks failed with phi > 2 pi
          continue;
        }

        double deltaPhi = RecoDecay::constrainAngle(particleVector.at(i).phi() - particleVectorAnti.at(j).phi(), -constants::math::PIHalf);
        double deltaEta = particleVector.at(i).eta() - particleVectorAnti.at(j).eta();
        registryData.fill(HIST("deltaPhiSEProtonAntiproton"), deltaPhi);
        registryData.fill(HIST("deltaPhiEtaSEProtonAntiproton"), deltaPhi, deltaEta);
        break;
      }
      fillMixedEventDeltas(particleVector.at(i), bufferAnti, 4, jetAxis);
      tempBuffer.emplace_back(std::make_pair(phiToAxis, etaToAxis));
    }
  }

  double getDeltaPhi(double a1, double a2)
  {
    double failedPhi = -999;
    if (std::isnan(a1) || std::isnan(a2) || a1 == failedPhi || a2 == failedPhi)
      return -999;
    double deltaPhi(0);
    double phi1 = RecoDecay::constrainAngle(a1, 0);
    double phi2 = RecoDecay::constrainAngle(a2, 0);
    double diff = std::abs(phi1 - phi2);

    if (diff <= constants::math::PI)
      deltaPhi = diff;
    if (diff > constants::math::PI)
      deltaPhi = constants::math::TwoPI - diff;

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
      ux = sign * std::abs(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // Protection 2
    if (py == 0 && px != 0) {
      ux = -(pz * pz) / px;
      uy = sign * std::abs(px * px - (pz * pz * pz * pz) / (px * px));
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
    ux = (-b + sign * std::abs(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }

  int analyseJet(int jetCounter, fastjet::PseudoJet jet, const auto& particles, auto& jetProtons, auto& jetAntiprotons, auto& jetPiPlus, auto& jetPiMinus, auto& jetAll, double rho, double rhoM, double rhoPerp, double rhoMPerp)
  {
    if (!jet.has_constituents())
      return jetCounter;
    fastjet::PseudoJet subtractedJetPerp(0., 0., 0., 0.);
    subtractedJetPerp = bkgSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
    fastjet::PseudoJet subtractedJetArea(0., 0., 0., 0.);
    subtractedJetArea = bkgSub.doRhoAreaSub(jet, rho, rhoM);

    if (subtractedJetPerp.pt() < minJetPt) // cut on jet w/o bkg
      return jetCounter;
    if ((std::fabs(jet.eta()) + jetR) > (maxEta - deltaEtaEdge))
      return jetCounter;
    jetCounter++;
    registryData.fill(HIST("ptTotalSubJetPerp"), subtractedJetPerp.pt());
    registryData.fill(HIST("ptTotalSubJetArea"), subtractedJetArea.pt());
    registryQC.fill(HIST("rhoEstimateArea"), jet.pt(), rho); // switch to subtracted jet pt
    registryQC.fill(HIST("rhoMEstimateArea"), jet.pt(), rhoM);
    registryQC.fill(HIST("rhoEstimatePerp"), jet.pt(), rhoPerp);
    registryQC.fill(HIST("rhoMEstimatePerp"), jet.pt(), rhoMPerp);
    double jetBkgDeltaPt = jet.pt() - subtractedJetPerp.pt();

    registryData.fill(HIST("eventProtocol"), 4);
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    registryData.fill(HIST("eventProtocol"), 5);
    registryData.fill(HIST("numberOfJets"), 0);
    registryData.fill(HIST("ptTotalJet"), jet.pt());
    registryData.fill(HIST("jetRapidity"), jet.rap());
    registryData.fill(HIST("numPartInJet"), jet.constituents().size());

    TVector3 pJet(0., 0., 0.);
    pJet.SetXYZ(jet.px(), jet.py(), jet.pz());

    if (outputQC) {
      registryQC.fill(HIST("jetBkgDeltaPt"), jetBkgDeltaPt);
      registryQC.fill(HIST("jetPtVsNumPart"), jet.pt(), jet.constituents().size());

      double maxRadius = 0;
      for (const auto& constituent : constituents) {
        registryData.fill(HIST("ptJetParticle"), constituent.pt());
        registryQC.fill(HIST("phiJet"), constituent.phi());
        registryQC.fill(HIST("phiPtJet"), constituent.pt(), constituent.phi());
        registryQC.fill(HIST("etaJet"), constituent.eta());
        registryQC.fill(HIST("etaPtJet"), constituent.pt(), constituent.eta());

        if (std::isnan(constituent.phi()) || std::isnan(jet.phi())) // geometric jet cone
          continue;
        double deltaPhi = RecoDecay::constrainAngle(constituent.phi() - jet.phi(), -constants::math::PIHalf);
        double deltaEta = constituent.eta() - jet.eta();
        double delta = std::abs(deltaPhi * deltaPhi + deltaEta * deltaEta);
        registryQC.fill(HIST("jetConeRadius"), delta);
        if (delta > maxRadius)
          maxRadius = delta;
      }
      registryQC.fill(HIST("maxRadiusVsPt"), jet.pt(), maxRadius);

      TVector3 ueAxis1(0.0, 0.0, 0.0);
      TVector3 ueAxis2(0.0, 0.0, 0.0);
      getPerpendicularAxis(pJet, ueAxis1, +1.0);
      getPerpendicularAxis(pJet, ueAxis2, -1.0);

      double nchJetPlusUE(0);
      double nchJet(0);
      double nchUE(0);
      double ptJetPlusUE(0);
      double ptJet(0);
      double ptUE(0);

      for (const auto& [index, track] : particles) {
        TVector3 particleDir(track.px(), track.py(), track.pz());
        double deltaEtaJet = particleDir.Eta() - pJet.Eta();
        double deltaPhiJet = getDeltaPhi(particleDir.Phi(), pJet.Phi());
        double deltaRJet = std::abs(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
        double deltaEtaUE1 = particleDir.Eta() - ueAxis1.Eta();
        double deltaPhiUE1 = getDeltaPhi(particleDir.Phi(), ueAxis1.Phi());
        double deltaRUE1 = std::abs(deltaEtaUE1 * deltaEtaUE1 + deltaPhiUE1 * deltaPhiUE1);
        double deltaEtaUE2 = particleDir.Eta() - ueAxis2.Eta();
        double deltaPhiUE2 = getDeltaPhi(particleDir.Phi(), ueAxis2.Phi());
        double deltaRUE2 = std::abs(deltaEtaUE2 * deltaEtaUE2 + deltaPhiUE2 * deltaPhiUE2);
        double failedPhi = -999;
        if (deltaRJet < rMax) {
          if (deltaPhiJet != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiJet"), deltaEtaJet, deltaPhiJet);
          nchJetPlusUE++;
          ptJetPlusUE = ptJetPlusUE + track.pt();
        }
        if (deltaRUE1 < rMax) {
          if (deltaPhiUE1 != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEtaUE1, deltaPhiUE1);
          nchUE++;
          ptUE = ptUE + track.pt();
        }
        if (deltaRUE2 < rMax) {
          if (deltaPhiUE2 != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEtaUE2, deltaPhiUE2);
          nchUE++;
          ptUE = ptUE + track.pt();
        }
      } // for (const auto& [index, track] : particles)

      nchJet = nchJetPlusUE - 0.5 * nchUE;
      ptJet = ptJetPlusUE - 0.5 * ptUE;
      registryQC.fill(HIST("multiplicityJetPlusUE"), nchJetPlusUE);
      registryQC.fill(HIST("multiplicityJet"), nchJet);
      registryQC.fill(HIST("multiplicityUE"), 0.5 * nchUE);
      registryQC.fill(HIST("ptJetPlusUE"), ptJetPlusUE);
      registryQC.fill(HIST("ptJet"), ptJet);
      registryQC.fill(HIST("ptUE"), 0.5 * ptUE);
      registryQC.fill(HIST("deltaJetPt"), jet.pt() - ptJetPlusUE);

      int nPartClusteredJet = static_cast<int>(constituents.size());

      // Fill QA Histograms
      if (ptJetPlusUE < minJetPt) { // swap for sub pt?

        registryQC.fill(HIST("nParticlesClusteredInJet"), nPartClusteredJet);

        for (const auto& track : constituents) {
          registryQC.fill(HIST("ptParticlesClusteredInJet"), track.pt());
        }
      }
    }

    std::vector<std::pair<double, double>> fTempBufferProton;
    std::vector<std::pair<double, double>> fTempBufferAntiproton;
    std::vector<std::pair<double, double>> fTempBufferPiPlus;
    std::vector<std::pair<double, double>> fTempBufferPiMinus;
    std::vector<std::pair<double, double>> fTempBufferJet;
    fTempBufferProton.clear();
    fTempBufferPiPlus.clear();
    fTempBufferPiMinus.clear();
    fTempBufferAntiproton.clear();
    fTempBufferJet.clear();

    for (const auto& pseudoParticle : constituents) { // analyse jet constituents - this is where the magic happens
      registryData.fill(HIST("trackProtocol"), 2);
      int id = pseudoParticle.user_index();
      const auto& jetParticle = particles.at(id);
      if (!selectTrack(jetParticle))
        continue;
      jetAll.emplace_back(jetParticle);

      registryData.fill(HIST("dcaXYFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("dcaZFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("tpcSignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
      if (jetParticle.hasTOF()) {
        registryData.fill(HIST("tofSignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());
      }
      if (outputQC) {
        double ptDiff = pseudoParticle.pt() - jetParticle.pt();
        registryQC.fill(HIST("ptDiff"), ptDiff);
      }

      if (isProton(jetParticle)) {
        registryData.fill(HIST("trackProtocol"), 3); // # high purity protons
        registryData.fill(HIST("ptJetProton"), jetParticle.pt());
        registryData.fill(HIST("dcaZJetProton"), jetParticle.pt(), jetParticle.dcaZ());
        registryQC.fill(HIST("ptJetProtonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        jetProtons.emplace_back(jetParticle);
      } else if (isAntiproton(jetParticle)) {
        registryData.fill(HIST("trackProtocol"), 4); // # high purity antiprotons
        registryData.fill(HIST("ptJetAntiproton"), jetParticle.pt());
        registryData.fill(HIST("dcaZJetAntiproton"), jetParticle.pt(), jetParticle.dcaZ());
        registryQC.fill(HIST("ptJetAntiprotonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        jetAntiprotons.emplace_back(jetParticle);
      } else if (isPion(jetParticle)) {
        registryQC.fill(HIST("ptJetPionVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("trackProtocol"), 5);
        registryData.fill(HIST("dcaZJetPion"), jetParticle.pt(), jetParticle.dcaZ());
        if (jetParticle.sign() > 0) {
          jetPiPlus.emplace_back(jetParticle);
        } else if (jetParticle.sign() < 0) {
          jetPiMinus.emplace_back(jetParticle);
        }
      }
      if (measureKaons && isKaon(jetParticle)) {
        registryQC.fill(HIST("ptJetKaonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("trackProtocol"), 6);
        registryData.fill(HIST("dcaZJetKaon"), jetParticle.pt(), jetParticle.dcaZ());
      }
    } // for (const auto& pseudoParticle : constituents)

    if (doJetCorrelations && jetAll.size() > 1) { // general correlation function
      doCorrelations(jetAll, fBufferJet, fTempBufferJet, 0, pJet);
      setTrackBuffer(fTempBufferJet, fBufferJet);
    }

    if (dopapCorrelations && (jetProtons.size() > 0) && (jetAntiprotons.size() > 0)) {
      doCorrelationsAnti(jetProtons, jetAntiprotons, fBufferAntiproton, fTempBufferProton, pJet);
      doCorrelationsAnti(jetAntiprotons, jetProtons, fBufferProton, fTempBufferAntiproton, pJet); // divide SE distributions by 2 in post
    }
    int minNumPartForCorrelations = 2;
    if ((static_cast<int>(jetProtons.size()) < minNumPartForCorrelations) && (static_cast<int>(jetAntiprotons.size()) < minNumPartForCorrelations) && (static_cast<int>(jetPiPlus.size()) < minNumPartForCorrelations) && (static_cast<int>(jetPiMinus.size()) < minNumPartForCorrelations))
      return jetCounter;
    registryData.fill(HIST("eventProtocol"), 6);

    if (doppCorrelations && jetProtons.size() > 1) {
      doCorrelations(jetProtons, fBufferProton, fTempBufferProton, 1, pJet);
      setTrackBuffer(fTempBufferProton, fBufferProton);
    }
    if (doapapCorrelations && jetAntiprotons.size() > 1) {
      doCorrelations(jetAntiprotons, fBufferAntiproton, fTempBufferAntiproton, 2, pJet);
      setTrackBuffer(fTempBufferAntiproton, fBufferAntiproton);
    }
    if (dopipiCorrelations && jetPiPlus.size() > 1) {
      doCorrelations(jetPiPlus, fBufferPiPlus, fTempBufferPiPlus, 3, pJet);
      setTrackBuffer(fTempBufferPiPlus, fBufferPiPlus);
    }
    if (dopipiCorrelations && jetPiMinus.size() > 1) {
      doCorrelations(jetPiMinus, fBufferPiMinus, fTempBufferPiMinus, 3, pJet);
      setTrackBuffer(fTempBufferPiMinus, fBufferPiMinus);
    }
    return jetCounter;
  }

  template <typename U>
  void fillHistograms(U const& tracks)
  {
    std::vector<typename U::iterator> jetProtons;
    std::vector<typename U::iterator> jetAntiprotons;
    std::vector<typename U::iterator> jetPiPlus;
    std::vector<typename U::iterator> jetPiMinus;
    std::vector<typename U::iterator> jetAll;
    jetProtons.clear();
    jetAntiprotons.clear();
    jetPiPlus.clear();
    jetPiMinus.clear();
    jetAll.clear();
    std::vector<std::pair<double, double>> fTempBufferFull;
    fTempBufferFull.clear();
    std::vector<fastjet::PseudoJet> jetInput;         // input for jet finder
    std::map<int, typename U::iterator> particles;    // all selected particles in event
    std::vector<typename U::iterator> particlesForCF; // particles for full event angular correlations
    jetInput.clear();
    particles.clear();
    int index = 0;
    int jetCounter = 0;
    std::vector<fastjet::PseudoJet> jets;
    jets.clear();

    for (const auto& track : tracks) {
      registryData.fill(HIST("trackProtocol"), 0); // # all tracks
      if (!selectTrackForJetReco(track))
        continue;

      registryData.fill(HIST("trackProtocol"), 1); // # tracks selected for jet reconstruction
      double mass = 0.139;

      if (outputQC && (track.tpcNClsFindable() != 0)) {
        registryQC.fill(HIST("ratioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
      }

      if (outputQC) {
        registryQC.fill(HIST("ptFullEvent"), track.pt());
        registryQC.fill(HIST("crossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
        registryQC.fill(HIST("clusterITS"), track.pt(), track.itsNCls());
        registryQC.fill(HIST("clusterTPC"), track.pt(), track.tpcNClsFound());
        registryQC.fill(HIST("chi2ITS"), track.pt(), track.itsChi2NCl());
        registryQC.fill(HIST("chi2TPC"), track.pt(), track.tpcChi2NCl());
        registryQC.fill(HIST("dcaXYFullEvent"), track.pt(), track.dcaXY());
        registryQC.fill(HIST("dcaZFullEvent"), track.pt(), track.dcaZ());
        registryQC.fill(HIST("phiFullEvent"), track.phi());
        registryQC.fill(HIST("phiPtFullEvent"), track.pt(), track.phi());
        registryQC.fill(HIST("etaFullEvent"), track.eta());
        registryQC.fill(HIST("etaPtFullEvent"), track.pt(), track.eta());
      }
      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(mass));
      inputPseudoJet.set_user_index(index);
      particles[index] = track;
      particlesForCF.emplace_back(track);
      jetInput.emplace_back(inputPseudoJet);

      index++;
    } // for (const auto& track : tracks)

    int minNumPartForJetReco = 2;
    if (static_cast<int>(jetInput.size()) < minNumPartForJetReco)
      return;
    registryData.fill(HIST("eventProtocol"), 2);

    // Reconstruct Jets
    double ghostMaxrap = 1.0;
    double ghostArea = 0.005;
    int ghostRepeat = 1;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, jetR);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghostMaxrap, ghostRepeat, ghostArea));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghostMaxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());

    if (jets.size() == 0)
      return;

    registryData.fill(HIST("eventProtocol"), 3);

    bool doSparse = true;
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(jetInput, doSparse);
    auto [rhoPerp, rhoMPerp] = bkgSub.estimateRhoPerpCone(jetInput, jets);

    for (const auto& jet : jets) {
      jetCounter = analyseJet(jetCounter, jet, particles, jetProtons, jetAntiprotons, jetPiPlus, jetPiMinus, jetAll, rho, rhoM, rhoPerp, rhoMPerp);
    }
    registryData.fill(HIST("numJetsInEvent"), jetCounter);

    TVector3 hardestJetAxis(jets.at(0).px(), jets.at(0).py(), jets.at(0).pz()); // for full event, use hardest jet as orientation
    doCorrelations(particlesForCF, fBufferFull, fTempBufferFull, -1, hardestJetAxis);
    setTrackBuffer(fTempBufferFull, fBufferFull);
  }

  template <typename U>
  void fillHistogramsMC(U const& tracks)
  {
    std::vector<fastjet::PseudoJet> jetInput;      // input for jet finder
    std::map<int, typename U::iterator> particles; // all selected particles in event
    jetInput.clear();
    particles.clear();
    int index = 0;
    double mass = 0.139; // pion mass for input because 80% is pions anyway
    std::vector<fastjet::PseudoJet> jets;
    jets.clear();

    for (const auto& track : tracks) {
      if (outputQC && (track.tpcNClsFindable() != 0)) {
        registryQC.fill(HIST("ratioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
      }
      registryQC.fill(HIST("ptFullEvent"), track.pt());
      registryQC.fill(HIST("crossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQC.fill(HIST("clusterITS"), track.pt(), track.itsNCls());
      registryQC.fill(HIST("clusterTPC"), track.pt(), track.tpcNClsFound());
      registryQC.fill(HIST("chi2ITS"), track.pt(), track.itsChi2NCl());
      registryQC.fill(HIST("chi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQC.fill(HIST("dcaXYFullEvent"), track.pt(), track.dcaXY());
      registryQC.fill(HIST("dcaZFullEvent"), track.pt(), track.dcaZ());
      registryQC.fill(HIST("phiFullEvent"), track.phi());
      registryQC.fill(HIST("phiPtFullEvent"), track.pt(), track.phi());
      registryQC.fill(HIST("etaFullEvent"), track.eta());
      registryQC.fill(HIST("etaPtFullEvent"), track.pt(), track.eta());

      fastjet::PseudoJet inputPseudoJet(track.px(), track.py(), track.pz(), track.energy(mass));
      inputPseudoJet.set_user_index(index);
      particles[index] = track;
      jetInput.emplace_back(inputPseudoJet);

      index++;
    } // for (const auto& track : tracks)

    int minNumPartForJetReco = 2;
    if (static_cast<int>(jetInput.size()) < minNumPartForJetReco)
      return;
    registryData.fill(HIST("eventProtocol"), 2);

    // Reconstruct Jets
    double ghostMaxrap = 1.0;
    double ghostArea = 0.005;
    int ghostRepeat = 1;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, jetR);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghostMaxrap, ghostRepeat, ghostArea));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghostMaxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());

    if (jets.size() == 0)
      return;

    registryData.fill(HIST("eventProtocol"), 3);

    bool doSparse = true;
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(jetInput, doSparse);
    auto [rhoPerp, rhoMPerp] = bkgSub.estimateRhoPerpCone(jetInput, jets);

    for (auto& jet : jets) { // o2-linter: disable=const-ref-in-for-loop (jets are modified)
      if (!jet.has_constituents())
        continue;
      fastjet::PseudoJet subtractedJetPerp(0., 0., 0., 0.);
      subtractedJetPerp = bkgSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
      fastjet::PseudoJet subtractedJetArea(0., 0., 0., 0.);
      subtractedJetArea = bkgSub.doRhoAreaSub(jet, rho, rhoM);

      if (subtractedJetPerp.pt() < minJetPt) // cut on jet w/o bkg
        continue;
      registryData.fill(HIST("ptTotalSubJetPerp"), subtractedJetPerp.pt());
      registryQC.fill(HIST("rhoEstimateArea"), jet.pt(), rho); // switch to subtracted jet pt
      registryQC.fill(HIST("rhoMEstimateArea"), jet.pt(), rhoM);
      registryQC.fill(HIST("rhoEstimatePerp"), jet.pt(), rhoPerp);
      registryQC.fill(HIST("rhoMEstimatePerp"), jet.pt(), rhoMPerp);
      double jetBkgDeltaPt = jet.pt() - subtractedJetPerp.pt();
      registryQC.fill(HIST("jetBkgDeltaPt"), jetBkgDeltaPt);

      registryData.fill(HIST("eventProtocol"), 4);
      std::vector<fastjet::PseudoJet> constituents = jet.constituents();

      registryData.fill(HIST("eventProtocol"), 5);
      registryData.fill(HIST("numberOfJets"), 0);
      registryData.fill(HIST("ptTotalJet"), jet.pt());
      registryData.fill(HIST("jetRapidity"), jet.rap());
      registryData.fill(HIST("numPartInJet"), jet.constituents().size());
      registryQC.fill(HIST("jetPtVsNumPart"), jet.pt(), jet.constituents().size());

      double maxRadius = 0;
      for (const auto& constituent : constituents) {
        registryData.fill(HIST("ptJetParticle"), constituent.pt());
        registryQC.fill(HIST("phiJet"), constituent.phi());
        registryQC.fill(HIST("phiPtJet"), constituent.pt(), constituent.phi());
        registryQC.fill(HIST("etaJet"), constituent.eta());
        registryQC.fill(HIST("etaPtJet"), constituent.pt(), constituent.eta());

        if (std::isnan(constituent.phi()) || std::isnan(jet.phi())) // geometric jet cone
          continue;
        double deltaPhi = RecoDecay::constrainAngle(constituent.phi() - jet.phi(), -constants::math::PIHalf);
        double deltaEta = constituent.eta() - jet.eta();
        double delta = std::abs(deltaPhi * deltaPhi + deltaEta * deltaEta);
        registryQC.fill(HIST("jetConeRadius"), delta);
        if (delta > maxRadius)
          maxRadius = delta;
      }
      registryQC.fill(HIST("maxRadiusVsPt"), jet.pt(), maxRadius);

      TVector3 pJet(0., 0., 0.);
      pJet.SetXYZ(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0.0, 0.0, 0.0);
      TVector3 ueAxis2(0.0, 0.0, 0.0);
      getPerpendicularAxis(pJet, ueAxis1, +1.0);
      getPerpendicularAxis(pJet, ueAxis2, -1.0);

      double nchJetPlusUE(0);
      double nchJet(0);
      double nchUE(0);

      for (const auto& [index, track] : particles) {
        TVector3 particleDir(track.px(), track.py(), track.pz());
        double deltaEtaJet = particleDir.Eta() - pJet.Eta();
        double deltaPhiJet = getDeltaPhi(particleDir.Phi(), pJet.Phi());
        double deltaRJet = std::abs(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
        double deltaEtaUE1 = particleDir.Eta() - ueAxis1.Eta();
        double deltaPhiUE1 = getDeltaPhi(particleDir.Phi(), ueAxis1.Phi());
        double deltaRUE1 = std::abs(deltaEtaUE1 * deltaEtaUE1 + deltaPhiUE1 * deltaPhiUE1);
        double deltaEtaUE2 = particleDir.Eta() - ueAxis2.Eta();
        double deltaPhiUE2 = getDeltaPhi(particleDir.Phi(), ueAxis2.Phi());
        double deltaRUE2 = std::abs(deltaEtaUE2 * deltaEtaUE2 + deltaPhiUE2 * deltaPhiUE2);

        double failedPhi = -999;
        if (deltaRJet < rMax) {
          if (deltaPhiJet != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiJet"), deltaEtaJet, deltaPhiJet);
          nchJetPlusUE++;
        }
        if (deltaRUE1 < rMax) {
          if (deltaPhiUE1 != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEtaUE1, deltaPhiUE1);
          nchUE++;
        }
        if (deltaRUE2 < rMax) {
          if (deltaPhiUE2 != failedPhi)
            registryQC.fill(HIST("deltaEtadeltaPhiUE"), deltaEtaUE2, deltaPhiUE2);
          nchUE++;
        }
      } // for (const auto& [index, track] : particles)

      nchJet = nchJetPlusUE - 0.5 * nchUE;
      registryQC.fill(HIST("multiplicityJetPlusUE"), nchJetPlusUE);
      registryQC.fill(HIST("multiplicityJet"), nchJet);
      registryQC.fill(HIST("multiplicityUE"), 0.5 * nchUE);

      for (const auto& pseudoParticle : constituents) { // analyse jet constituents - this is where the magic happens
        registryData.fill(HIST("trackProtocol"), 3);
        int id = pseudoParticle.user_index();
        const auto& jetParticle = particles.at(id);
        if (!selectTrack(jetParticle))
          continue;

        registryData.fill(HIST("dcaXYFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
        registryData.fill(HIST("dcaZFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
        registryData.fill(HIST("tpcSignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
        if (jetParticle.hasTOF()) {
          registryData.fill(HIST("tofSignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());
        }
        if (outputQC) {
          double ptDiff = pseudoParticle.pt() - jetParticle.pt();
          registryQC.fill(HIST("ptDiff"), ptDiff);
        }

        // if (jetParticle.pt() < minJetParticlePt)
        //   continue;
        if (!jetParticle.has_mcParticle())
          continue;
        switch (jetParticle.mcParticle().pdgCode()) {
          case kProton:
            registryMC.fill(HIST("numberOfTruthParticles"), 0);
            registryMC.fill(HIST("ptJetProtonMC"), jetParticle.pt());
            break;
          case kProtonBar:
            registryMC.fill(HIST("numberOfTruthParticles"), 1);
            registryMC.fill(HIST("ptJetAntiprotonMC"), jetParticle.pt());
            break;
          default:
            continue;
        }

        if (!selectTrackForJetReco(jetParticle))
          continue;
        switch (jetParticle.mcParticle().pdgCode()) {
          case kProton:
            registryData.fill(HIST("ptJetProton"), jetParticle.pt());
            registryQC.fill(HIST("ptJetProtonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
            registryData.fill(HIST("trackProtocol"), 4); // # protons
            break;
          case kProtonBar:
            registryData.fill(HIST("ptJetAntiproton"), jetParticle.pt());
            registryQC.fill(HIST("ptJetAntiprotonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
            registryData.fill(HIST("trackProtocol"), 6); // # antiprotons
            break;
          default:
            continue;
        }
      } // for (const auto& pseudoParticle : constituents)
    } // for (auto& jet : jets)
  }

  void processRun2(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun2> const& tracks,
                   BCsWithRun2Info const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      registryData.fill(HIST("eventProtocol"), 0);
      if (!collision.alias_bit(kINT7))
        continue;
      registryData.fill(HIST("numberOfEvents"), 0);
      registryData.fill(HIST("eventProtocol"), 1);

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun2, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun2, "process Run 2 data w/o jet tables", false);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun3> const& tracks)
  {
    for (const auto& collision : collisions) {
      registryData.fill(HIST("eventProtocol"), 0);
      if (!collision.sel8())
        continue;
      registryData.fill(HIST("numberOfEvents"), 0);
      registryData.fill(HIST("eventProtocol"), 1);
      if (std::abs(collision.posZ()) > zVtx)
        continue;

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun3, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun3, "process Run 3 data w/o jet tables", false);

  void processMCRun2(McCollisions const& collisions, soa::Filtered<McTracksRun2> const& tracks, BCsWithRun2Info const&, aod::McParticles const&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      registryData.fill(HIST("eventProtocol"), 0);
      if (!collision.alias_bit(kINT7))
        continue;
      registryData.fill(HIST("numberOfEvents"), 0);
      registryData.fill(HIST("eventProtocol"), 1);

      auto slicedTracks = tracks.sliceBy(perCollisionMcTracksRun2, collision.globalIndex());

      fillHistogramsMC(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processMCRun2, "process Run 2 MC w/o jet tables, not currently usable", false);

  void processMCRun3(McCollisions const& collisions, soa::Filtered<McTracksRun3> const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      registryData.fill(HIST("eventProtocol"), 0);
      if (!collision.sel8())
        continue;
      registryData.fill(HIST("numberOfEvents"), 0);
      registryData.fill(HIST("eventProtocol"), 1);
      if (std::abs(collision.posZ()) > zVtx)
        continue;

      auto slicedTracks = tracks.sliceBy(perCollisionMcTracksRun3, collision.globalIndex());

      fillHistogramsMC(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processMCRun3, "process Run 3 MC w/o jet tables, not currently usable", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc)};
}
