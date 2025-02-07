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
#include <utility>
#include <map>

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
#include "Common/DataModel/McCollisionExtra.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGJE/DataModel/Jet.h"

#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "TVector2.h"
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
  // Preliminary Cuts
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of crossed rows TPC"};
  Configurable<int> minReqClusterITS{"minReqClusterITS", 2, "min number of clusters required in ITS"};
  Configurable<int> minReqClusterTPC{"minReqClusterTPC", 70, "min number of clusters required in TPC"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.7, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0, "max chi2 per cluster ITS"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0, "max chi2 per cluster TPC"};
  Configurable<float> maxDCAxy{"maxDCAxy", 0.5, "max DCA to vertex xy"};
  Configurable<float> maxDCAz{"maxDCAz", 1.0, "max DCA to vertex z"};
  Configurable<float> maxEta{"maxEta", 0.8, "max pseudorapidity"}; // consider jet cone?
  Configurable<float> minTrackPt{"minTrackPt", 0.1, "minimum track pT"};

  // Jet Cuts
  Configurable<float> jetR{"jetR", 0.4, "jet resolution parameter"};
  Configurable<float> minJetPt{"minJetPt", 5.0, "minimum total pT to accept jet"};
  Configurable<float> minJetParticlePt{"minJetParticlePt", 0.0, "minimum pT to accept jet particle"};

  // Proton Cuts
  Configurable<float> protonDCAxyYield{"protonDCAxyYield", 0.05, "[proton] DCAxy cut for yield"};
  Configurable<float> protonDCAzYield{"protonDCAzYield", 1.0, "[proton] DCAz cut for yield"};
  Configurable<float> protonDCAxyCF{"protonDCAxyCF", 0.05, "[proton] DCAxy cut for CF"};
  Configurable<float> protonDCAzCF{"protonDCAzCF", 1.0, "[proton] DCAz cut for CF"};
  Configurable<float> protonTPCTOFpT{"protonTPCTOFpT", 0.7, "[proton] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> protonTPCnsigmaLowPtYield{"protonTPCnsigmaLowPtYield", 4.0, "[proton] max TPC nsigma with low pT for yield"};
  Configurable<float> protonTPCnsigmaHighPtYield{"protonTPCnsigmaHighPtYield", 4.0, "[proton] max TPC nsigma with high pT for yield"};
  Configurable<float> protonTOFnsigmaHighPtYield{"protonTOFnsigmaHighPtYield", 4.0, "[proton] max TOF nsigma with high pT yield"};
  Configurable<float> protonNsigma{"protonNsigma", 2.0, "[proton] max combined nsigma for CF (sqrt(nsigTPC^2 + nsigTOF^2))"};

  // Antiproton Cuts
  Configurable<float> antiprotonDCAxyYield{"antiprotonDCAxyYield", 0.05, "[antiproton] DCAxy cut for yield"};
  Configurable<float> antiprotonDCAzYield{"antiprotonDCAzYield", 1.0, "[antiproton] DCAz cut for yield"};
  Configurable<float> antiprotonDCAxyCF{"antiprotonDCAxyCF", 0.05, "[antiproton] DCAxy cut for CF"};
  Configurable<float> antiprotonDCAzCF{"antiprotonDCAzCF", 1.0, "[antiproton] DCAz cut for CF"};
  Configurable<float> antiprotonTPCTOFpT{"antiprotonTPCTOFpT", 0.7, "[antiproton] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> antiprotonTPCnsigmaLowPtYield{"antiprotonTPCnsigmaLowPtYield", 4.0, "[antiproton] max TPC nsigma with low pT for yield"};
  Configurable<float> antiprotonTPCnsigmaHighPtYield{"antiprotonTPCnsigmaHighPtYield", 4.0, "[antiproton] max TPC nsigma with high pT for yield"};
  Configurable<float> antiprotonTOFnsigmaHighPtYield{"antiprotonTOFnsigmaHighPtYield", 4.0, "[antiproton] min TOF nsigma with high pT for yield"};
  Configurable<float> antiprotonNsigma{"antiprotonNsigma", 2.0, "[antiproton] max combined nsigma for CF (sqrt(nsigTPC^2 + nsigTOF^2))"};

  // Nuclei Cuts
  Configurable<float> nucleiDCAxyYield{"nucleiDCAxyYield", 0.05, "[nuclei] DCAxy cut for yield"};
  Configurable<float> nucleiDCAzYield{"nucleiDCAzYield", 0.02, "[nuclei] DCAz cut for yield"};
  Configurable<float> nucleiTPCTOFpT{"nucleiTPCTOFpT", 0.7, "[nuclei] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> nucleiTPCnsigmaLowPtYield{"nucleiTPCnsigmaLowPtYield", 4.0, "[nuclei] max TPC nsigma with low pT for yield"};
  Configurable<float> nucleiTPCnsigmaHighPtYield{"nucleiTPCnsigmaHighPtYield", 4.0, "[nuclei] max TPC nsigma with high pT for yield"};
  Configurable<float> nucleiTOFnsigmaHighPtYield{"nucleiTOFnsigmaHighPtYield", 4.0, "[nuclei] min TOF nsigma with high pT for yield"};

  // Antinuclei Cuts
  Configurable<float> antinucleiDCAxyYield{"antinucleiDCAxyYield", 0.05, "[antinuclei] DCAxy cut for yield"};
  Configurable<float> antinucleiDCAzYield{"antinucleiDCAzYield", 0.02, "[antinuclei] DCAz cut for yield"};
  Configurable<float> antinucleiTPCTOFpT{"antinucleiTPCTOFpT", 0.7, "[antinuclei] pT for switch in TPC/TPC+TOF nsigma"};
  Configurable<float> antinucleiTPCnsigmaLowPtYield{"antinucleiTPCnsigmaLowPtYield", 4.0, "[antinuclei] max TPC nsigma with low pT for yield"};
  Configurable<float> antinucleiTPCnsigmaHighPtYield{"antinucleiTPCnsigmaHighPtYield", 4.0, "[antinuclei] max TPC nsigma with high pT for yield"};
  Configurable<float> antinucleiTOFnsigmaHighPtYield{"antinucleiTOFnsigmaHighPtYield", 4.0, "[antinuclei] min TOF nsigma with high pT for yield"};

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

  Configurable<bool> useRejectionCut{"useRejectionCut", true, "use nsigmaRejection"};
  Configurable<float> nsigmaRejection{"nsigmaRejection", 1.0, "reject tracks with nsigma < nsigmaRejection for >1 species"};
  Configurable<bool> deuteronAnalysis{"deuteronAnalysis", true, "true [false]: analyse (anti)deuterons [(anti)helium-3]"};
  Configurable<bool> useTOFmass{"useTOFmass", false, "use TOF mass instead of pion mass if available"};

  Configurable<int> trackBufferSize{"trackBufferSize", 200, "Number of mixed-event tracks being stored"};

  // QC Configurables
  Configurable<float> zVtx{"zVtx", 10.0, "max zVertex"};
  Configurable<float> Rmax{"Rmax", 0.4, "Maximum radius for jet and UE regions"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber;

  using FullTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                   aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa>;
  using FullTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                   aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa>;
  using McTracksRun2 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCov, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection,
                                 aod::TrackSelectionExtension, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa, aod::McTrackLabels>;
  using McTracksRun3 = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection, aod::TrackSelectionExtension,
                                 aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullHe, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullHe, aod::pidTOFmass, aod::pidTOFbeta, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCTr, aod::pidTPCAl, aod::pidTOFPi, aod::pidTOFKa, aod::McTrackLabels>;
  using JTracksRun3 = soa::Join<aod::JetTracks, aod::JTrackExtras, aod::JTrackPIs>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;
  using McCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

  Filter prelimTrackCuts = (aod::track::itsChi2NCl < maxChi2ITS &&
                            aod::track::tpcChi2NCl < maxChi2TPC &&
                            nabs(aod::track::dcaXY) < maxDCAxy &&
                            nabs(aod::track::dcaZ) < maxDCAz &&
                            nabs(aod::track::eta) < maxEta &&
                            aod::track::pt > minTrackPt); // add more preliminary cuts to filter if possible
  Filter collisionFilter = (nabs(aod::jcollision::posZ) < zVtx);
  Filter trackCuts = (nabs(aod::jtrack::eta) > maxEta && aod::jtrack::pt > minJetParticlePt);
  // Filter partCuts = (aod::jmcparticle::pt >= trackPtMin && aod::jmcparticle::pt < trackPtMax);
  Filter jetFilter = (aod::jet::pt >= minJetPt && nabs(aod::jet::eta) < nabs(maxEta - aod::jet::r / 100.f));

  Preslice<FullTracksRun2> perCollisionFullTracksRun2 = o2::aod::track::collisionId;
  Preslice<FullTracksRun3> perCollisionFullTracksRun3 = o2::aod::track::collisionId;
  Preslice<McTracksRun2> perCollisionMcTracksRun2 = o2::aod::track::collisionId;
  Preslice<McTracksRun3> perCollisionMcTracksRun3 = o2::aod::track::collisionId;

  AxisSpecs axisSpecs;

  HistogramRegistry registryData{"dataOutput", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryMC{"MCOutput", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registryQA{"dataQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  JetBkgSubUtils bkgSub;
  std::vector<int> eventSelection;

  void init(o2::framework::InitContext&)
  {
    mRunNumber = 0;
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    eventSelection = jetderiveddatautilities::initialiseEventSelectionBits(static_cast<std::string>("sel8"));

    // Counters
    registryData.add("hNumberOfEvents", "Number of events", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hNumberOfJets", "Total number of jets", HistType::kTH1I, {{1, 0, 1}});
    registryData.add("hEventProtocol", "Event protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hTrackProtocol", "Track protocol", HistType::kTH1I, {{20, 0, 20}});
    registryData.add("hNumPartInJet", "Number of particles in a jet", HistType::kTH1I, {{200, 0, 200}});
    registryData.add("hNumJetsInEvent", "Number of jets selected", HistType::kTH1I, {{10, 0, 10}});

    // (Pseudo)Rapidity
    registryData.add("hJetRapidity", "Jet rapidity;#it{y}", HistType::kTH1F, {{200, -1, 1}});

    // pT
    registryData.add("hPtJetParticle", "p_{T} of particles in jets", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtTotalSubJetArea", "Subtracted full jet p_{T} (area)", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtTotalSubJetPerp", "Subtracted full jet p_{T} (perpendicular)", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetProton", "p_{T} of protons", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetAntiproton", "p_{T} of antiprotons", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetNuclei", "p_{T} of nuclei", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtJetAntinuclei", "p_{T} of antinuclei", HistType::kTH1D, {axisSpecs.ptAxisPos});
    registryData.add("hPtTotalJet", "p_{T} of entire jet;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1F, {{1000, 0, 500}});
    registryQA.add("hPtJetProtonVsTotalJet", "Proton p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQA.add("hPtJetAntiprotonVsTotalJet", "Antiproton p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQA.add("hPtJetNucleiVsTotalJet", "Nuclei p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQA.add("hPtJetAntinucleiVsTotalJet", "Antinuclei p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQA.add("hPtJetPionVsTotalJet", "Pion p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});
    registryQA.add("hPtJetKaonVsTotalJet", "Kaon p_{T} vs. jet p_{T}", HistType::kTH2D, {axisSpecs.ptAxisPos, {1000, 0, 500, "jet p_{T} [GeV/#it{c}]"}});

    // nSigma
    registryData.add("hTPCsignal", "TPC signal", HistType::kTH2F, {{1000, -100, 100, "#it{p} [GeV/#it{c}]"}, {1000, 0, 5000, "d#it{E}/d#it{X} (a.u.)"}});
    registryData.add("hTOFsignal", "TOF signal", HistType::kTH2F, {{1000, -100, 100, "#it{p} [GeV/#it{c}]"}, {550, 0, 1.1, "#beta (TOF)"}});
    registryData.add("hTPCnsigmaProton", "TPC n#sigma for proton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaProton", "TOF n#sigma for proton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaAntiproton", "TPC n#sigma for antiproton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaAntiproton", "TOF n#sigma for antiproton", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaNuclei", "TPC n#sigma for nuclei", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaNuclei", "TOF n#sigma for nuclei", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaAntinuclei", "TPC n#sigma for antinuclei", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaAntinuclei", "TOF n#sigma for antinuclei", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaProtonCF", "TPC n#sigma for proton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaProtonCF", "TOF n#sigma for proton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaAntiprotonCF", "TPC n#sigma for antiproton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaAntiprotonCF", "TOF n#sigma for antiproton CF", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaPion", "TPC n#sigma for pion", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaPion", "TOF n#sigma for pion", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTPCnsigmaKaon", "TPC n#sigma for kaon", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});
    registryData.add("hTOFnsigmaKaon", "TOF n#sigma for kaon", HistType::kTH2F, {axisSpecs.nsigmapTAxis, axisSpecs.nsigmaAxis});

    // DCA
    registryData.add("hDCAxyFullJet", "DCA_{xy} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcaxyAxis});
    registryData.add("hDCAzFullJet", "DCA_{z} of full jet", HistType::kTH2F, {axisSpecs.ptAxisFull, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetProton", "DCA_{z} of high purity protons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetAntiproton", "DCA_{z} of high purity antiprotons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetNuclei", "DCA_{z} of high purity nuclei", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetAntinuclei", "DCA_{z} of high purity antinuclei", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetPion", "DCA_{z} of high purity pions", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});
    registryData.add("hDCAzJetKaon", "DCA_{z} of high purity kaons", HistType::kTH2F, {axisSpecs.ptAxisPos, axisSpecs.dcazAxis});

    // Angular Distributions
    registryQA.add("hPhiFullEvent", "#varphi in full event", HistType::kTH1F, {{1000, 0, 6.3}});
    registryQA.add("hPhiPtFullEvent", "#varphi vs. p_{T} in full event", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 6.3}});
    registryQA.add("hPhiJet", "#varphi in jet", HistType::kTH1F, {{1000, 0, 6.3}});
    registryQA.add("hPhiPtJet", "#varphi vs. p_{T} in jet", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 6.3}});
    registryQA.add("hEtaFullEvent", "#eta in full event", HistType::kTH1F, {{1000, -1, 1}});
    registryQA.add("hEtaPtFullEvent", "#eta vs. p_{T} in full event", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, -1, 1}});
    registryQA.add("hEtaJet", "#eta in jet", HistType::kTH1F, {{1000, -1, 1}});
    registryQA.add("hEtaPtJet", "#eta vs. p_{T} in jet", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, -1, 1}});

    registryData.add("hDeltaPhiSEFull", "#Delta#varphi of particles in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEJet", "#Delta#varphi of jet particles in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEProton", "#Delta#varphi of protons in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiSEAntiproton", "#Delta#varphi of antiprotons in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEFull", "#Delta#varphi of particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEJet", "#Delta#varphi of jet particles in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEProton", "#Delta#varphi of protons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEAntiproton", "#Delta#varphi of antiprotons in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiEtaSEFull", "#Delta#varphi vs #Delta#eta of full particles in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEJet", "#Delta#varphi vs #Delta#eta of jet particles in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEProton", "#Delta#varphi vs #Delta#eta of protons in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaSEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEFull", "#Delta#varphi vs #Delta#eta of particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEJet", "#Delta#varphi vs #Delta#eta of jet particles in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEProton", "#Delta#varphi vs #Delta#eta of protons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEAntiproton", "#Delta#varphi vs #Delta#eta of antiprotons in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});

    registryData.add("hDeltaPhiSEProtonAntiproton", "#Delta#varphi of proton-antiproton in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEProtonAntiproton", "#Delta#varphi of proton-antiproton in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiEtaSEProtonAntiproton", "#Delta#varphi vs #Delta#eta of proton-antiproton in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEProtonAntiproton", "#Delta#varphi vs #Delta#eta of proton-antiproton in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiSEPion", "#Delta#varphi of pions in same event", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiMEPion", "#Delta#varphi of pions in mixed events", HistType::kTH1D, {axisSpecs.angDistPhiAxis});
    registryData.add("hDeltaPhiEtaSEPion", "#Delta#varphi vs #Delta#eta of pions in same event", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});
    registryData.add("hDeltaPhiEtaMEPion", "#Delta#varphi vs #Delta#eta of pions in mixed events", HistType::kTH2D, {axisSpecs.angDistPhiAxis, axisSpecs.angDistEtaAxis});

    // QA
    registryQA.add("hPtDiff", "p_{T} difference PseudoJet/original track;#it{p}_{T} [GeV/#it{c}]", HistType::kTH1D, {{100, -5, 5}});
    registryQA.add("hJetConeRadius", "Jet Radius;#it{R}", HistType::kTH1F, {{100, 0, 1}});
    registryQA.add("hMaxRadiusVsPt", "Max Cone Radius vs p_{T}", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {100, 0, 1}});
    registryQA.add("hRhoEstimatePerp", "Background #rho (perp)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQA.add("hRhoMEstimatePerp", "Background #rho_{m} (perp)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQA.add("hRhoEstimateArea", "Background #rho (area)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQA.add("hRhoMEstimateArea", "Background #rho_{m} (area)", HistType::kTH2F, {{axisSpecs.ptAxisPos}, {200, 0, 20}});
    registryQA.add("hJetBkgDeltaPt", "#Delta p_{T} Clustered Cone - Pure Jet", HistType::kTH1F, {{200, 0, 10}});

    registryQA.add("hTOFmass", "TOF mass vs p_{T}", HistType::kTH2F, {axisSpecs.ptAxisPos, {1000, 0, 5, "#it{m} [GeV/#it{c}^{2}]"}});
    registryQA.get<TH2>(HIST("hTOFmass"))->Sumw2();
    registryQA.add("hPtFullEvent", "p_{T} after basic cuts", HistType::kTH1F, {axisSpecs.ptAxisPos});
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
    registryQA.add("hPtJetPlusUE", "hPtJetPlusUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hPtJet", "hPtJet", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hPtUE", "hPtUE", HistType::kTH1F, {{500, 0, 50, "#it{p}_{T} (GeV/#it{c})"}});
    registryQA.add("hDeltaEtadeltaPhiJet", "hDeltaEtadeltaPhiJet", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaEtadeltaPhiUE", "hDeltaEtadeltaPhiUE", HistType::kTH2F, {{200, -0.5, 0.5, "#Delta#eta"}, {200, 0, 0.5 * TMath::Pi(), "#Delta#phi"}});
    registryQA.add("hDeltaJetPt", "hDeltaJetPt", HistType::kTH1F, {{200, -2, 2, "#Delta#it{p}_{T} (GeV/#it{c})"}});

    registryQA.add("hNParticlesClusteredInJet", "hNParticlesClusteredInJet", HistType::kTH1F, {{50, 0, 50, "#it{N}_{ch}"}});
    registryQA.add("hPtParticlesClusteredInJet", "hPtParticlesClusteredInJet", HistType::kTH1F, {{200, 0, 10, "#it{p}_{T} (GeV/#it{c})"}});

    // MC
    registryMC.add("hPtJetProtonMC", "Truth jet proton p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
    registryMC.add("hPtJetAntiprotonMC", "Truth jet antiproton p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
    registryMC.add("hPtJetNucleiMC", "Truth jet nuclei p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
    registryMC.add("hPtJetAntinucleiMC", "Truth jet antinuclei p_{T}", HistType::kTH1F, {axisSpecs.ptAxisPos});
    registryMC.add("hNumberOfTruthParticles", "Truth yields (anti)p, (anti)d, (anti)He-3", HistType::kTH1I, {{6, 0, 6}});
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

  template <class T>
  bool selectTrack(T const& track) // preliminary track selections
  {
    if (track.tpcNClsCrossedRows() < minRatioCrossedRowsTPC * track.tpcNClsFindable() ||
        track.tpcNClsCrossedRows() < minNCrossedRowsTPC ||
        track.tpcNClsFound() < minReqClusterTPC ||
        track.itsNCls() < minReqClusterITS) {
      return false;
    }
    if (doprocessRun2 || doprocessMCRun2) {
      if (!(track.trackType() & o2::aod::track::Run2Track) ||
          !(track.flags() & o2::aod::track::TPCrefit) ||
          !(track.flags() & o2::aod::track::ITSrefit)) {
        return false;
      }
    }
    return true;
  }

  template <class T>
  bool singleSpeciesTPCNSigma(T const& track)              // make cut configurable
  {                                                        // reject any track that has nsigma < 3 for more than 1 species
    if (useRejectionCut && (track.tpcNSigmaStoreEl() < nsigmaRejection || track.tpcNSigmaStoreMu() < nsigmaRejection || track.tpcNSigmaStorePi() < nsigmaRejection || track.tpcNSigmaStoreKa() < nsigmaRejection || track.tpcNSigmaStoreTr() < nsigmaRejection || track.tpcNSigmaStoreAl() < nsigmaRejection || track.tpcNSigmaDe() < nsigmaRejection || track.tpcNSigmaHe() < nsigmaRejection))
      return false;
    return true;
  }

  template <typename T>
  bool isProton(const T& track, bool tightCuts)
  {
    if (track.sign() < 0)
      return false;

    if (tightCuts) { // for correlation function
      // DCA
      if (TMath::Abs(track.dcaXY()) > protonDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > protonDCAzCF)
        return false;

      registryData.fill(HIST("hTPCnsigmaProtonCF"), track.pt(), track.tpcNSigmaPr());
      if (track.hasTOF())
        registryData.fill(HIST("hTOFnsigmaProtonCF"), track.pt(), track.tofNSigmaPr());

      // nsigma
      if (!track.hasTOF())
        return false;
      if ((track.pt() < protonTPCTOFpT && (TMath::Abs(track.tpcNSigmaPr()) > protonNsigma)) ||
          (track.pt() > protonTPCTOFpT && (TMath::Sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr()) > protonNsigma)))
        return false;
      if (useRejectionCut && !singleSpeciesTPCNSigma(track))
        return false;
    } else { // for yields
      // DCA
      if (TMath::Abs(track.dcaXY()) > protonDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > protonDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaProton"), track.pt(), track.tpcNSigmaPr());

      // TPC
      if (track.pt() < protonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > protonTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > protonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > protonTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaProton"), track.pt(), track.tofNSigmaPr());
        if (track.pt() > protonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > protonTOFnsigmaHighPtYield)
          return false;
      }
    }

    return true;
  }

  template <typename T>
  bool isAntiproton(const T& track, bool tightCuts)
  {
    if (track.sign() > 0)
      return false;

    if (tightCuts) { // for correlation function
      // DCA
      if (TMath::Abs(track.dcaXY()) > antiprotonDCAxyCF)
        return false;
      if (TMath::Abs(track.dcaZ()) > antiprotonDCAzCF)
        return false;

      registryData.fill(HIST("hTPCnsigmaAntiprotonCF"), track.pt(), track.tpcNSigmaPr());
      if (track.hasTOF())
        registryData.fill(HIST("hTOFnsigmaAntiprotonCF"), track.pt(), track.tofNSigmaPr());

      // nsigma
      if (!track.hasTOF())
        return false;
      if ((track.pt() < antiprotonTPCTOFpT && (TMath::Abs(track.tpcNSigmaPr()) > antiprotonNsigma)) ||
          (track.pt() > antiprotonTPCTOFpT && (TMath::Sqrt(track.tpcNSigmaPr() * track.tpcNSigmaPr() + track.tofNSigmaPr() * track.tofNSigmaPr()) > antiprotonNsigma)))
        return false;
      if (useRejectionCut && !singleSpeciesTPCNSigma(track))
        return false;
    } else { // for yields
      // DCA
      if (TMath::Abs(track.dcaXY()) > antiprotonDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > antiprotonDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaAntiproton"), track.pt(), track.tpcNSigmaPr());

      // TPC
      if (track.pt() < antiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > antiprotonTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > antiprotonTPCTOFpT && TMath::Abs(track.tpcNSigmaPr()) > antiprotonTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaAntiproton"), track.pt(), track.tofNSigmaPr());
        if (track.pt() > antiprotonTPCTOFpT && TMath::Abs(track.tofNSigmaPr()) > antiprotonTOFnsigmaHighPtYield)
          return false;
      }
    }

    return true;
  }

  template <typename T>
  bool isNucleus(const T& track)
  {
    if (track.sign() < 0)
      return false;
    if (deuteronAnalysis) {
      // DCA
      if (TMath::Abs(track.dcaXY()) > nucleiDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > nucleiDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaNuclei"), track.pt(), track.tpcNSigmaDe());

      // TPC
      if (track.pt() < nucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > nucleiTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > nucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > nucleiTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaNuclei"), track.pt(), track.tofNSigmaDe());
        if (track.pt() > nucleiTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > nucleiTOFnsigmaHighPtYield)
          return false;
      }
    } else {
      // DCA
      if (TMath::Abs(track.dcaXY()) > nucleiDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > nucleiDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaNuclei"), track.pt(), track.tpcNSigmaHe());

      // TPC
      if (track.pt() < nucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > nucleiTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > nucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > nucleiTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaNuclei"), track.pt(), track.tofNSigmaHe());
        if (track.pt() > nucleiTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > nucleiTOFnsigmaHighPtYield)
          return false;
      }
    }

    return true;
  }

  template <typename T>
  bool isAntinucleus(const T& track)
  {
    if (track.sign() > 0)
      return false;

    if (deuteronAnalysis) {
      // DCA
      if (TMath::Abs(track.dcaXY()) > antinucleiDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > antinucleiDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaAntinuclei"), track.pt(), track.tpcNSigmaDe());

      // TPC
      if (track.pt() < antinucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > antinucleiTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > antinucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaDe()) > antinucleiTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaAntinuclei"), track.pt(), track.tofNSigmaDe());
        if (track.pt() > antinucleiTPCTOFpT && TMath::Abs(track.tofNSigmaDe()) > antinucleiTOFnsigmaHighPtYield)
          return false;
      }
    } else {
      // DCA
      if (TMath::Abs(track.dcaXY()) > antinucleiDCAxyYield)
        return false;
      if (TMath::Abs(track.dcaZ()) > antinucleiDCAzYield)
        return false;

      registryData.fill(HIST("hTPCnsigmaAntinuclei"), track.pt(), track.tpcNSigmaHe());

      // TPC
      if (track.pt() < antinucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > antinucleiTPCnsigmaLowPtYield)
        return false;
      if (track.pt() > antinucleiTPCTOFpT && TMath::Abs(track.tpcNSigmaHe()) > antinucleiTPCnsigmaHighPtYield)
        return false;

      // TOF
      if (track.hasTOF()) {
        registryData.fill(HIST("hTOFnsigmaAntinuclei"), track.pt(), track.tofNSigmaHe());
        if (track.pt() > antinucleiTPCTOFpT && TMath::Abs(track.tofNSigmaHe()) > antinucleiTOFnsigmaHighPtYield)
          return false;
      }
    }

    return true;
  }

  template <typename T>
  bool isPion(const T& track)
  {
    // DCA
    if (TMath::Abs(track.dcaXY()) > pionDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > pionDCAz)
      return false;

    registryData.fill(HIST("hTPCnsigmaPion"), track.pt(), track.tpcNSigmaStorePi());

    // TPC
    if (track.pt() < pionTPCTOFpT && TMath::Abs(track.tpcNSigmaStorePi()) > pionTPCnsigmaLowPt)
      return false;
    if (track.pt() > pionTPCTOFpT && TMath::Abs(track.tpcNSigmaStorePi()) > pionTPCnsigmaHighPt)
      return false;

    // TOF
    if (track.hasTOF()) {
      registryData.fill(HIST("hTOFnsigmaPion"), track.pt(), track.tofNSigmaStorePi());
      if (track.pt() > pionTPCTOFpT && TMath::Abs(track.tofNSigmaStorePi()) > pionTOFnsigma)
        return false;
    }

    return true;
  }

  template <typename T>
  bool isKaon(const T& track)
  {
    // DCA
    if (TMath::Abs(track.dcaXY()) > kaonDCAxy)
      return false;
    if (TMath::Abs(track.dcaZ()) > kaonDCAz)
      return false;

    registryData.fill(HIST("hTPCnsigmaKaon"), track.pt(), track.tpcNSigmaStoreKa());

    // TPC
    if (track.pt() < kaonTPCTOFpT && TMath::Abs(track.tpcNSigmaStoreKa()) > kaonTPCnsigmaLowPt)
      return false;
    if (track.pt() > kaonTPCTOFpT && TMath::Abs(track.tpcNSigmaStoreKa()) > kaonTPCnsigmaHighPt)
      return false;

    // TOF
    if (track.hasTOF()) {
      registryData.fill(HIST("hTOFnsigmaKaon"), track.pt(), track.tofNSigmaStoreKa());
      if (track.pt() > kaonTPCTOFpT && TMath::Abs(track.tofNSigmaStoreKa()) > kaonTOFnsigma)
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
          registryData.fill(HIST("hDeltaPhiMEPion"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEPion"), DeltaPhi, DeltaEta);
          break;
        case 4:
          registryData.fill(HIST("hDeltaPhiMEProtonAntiproton"), DeltaPhi);
          registryData.fill(HIST("hDeltaPhiEtaMEProtonAntiproton"), DeltaPhi, DeltaEta);
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
        double DeltaEta = particleVector.at(i).eta() - particleVector.at(j).eta();
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
            registryData.fill(HIST("hDeltaPhiSEPion"), DeltaPhi);
            registryData.fill(HIST("hDeltaPhiEtaSEPion"), DeltaPhi, DeltaEta);
            break;
        }
      }
      fillMixedEventDeltas(particleVector.at(i), buffer, particleType, jetAxis);
      tempBuffer.emplace_back(std::make_pair(phiToAxis, etaToAxis));
    }
  }

  void doCorrelationsAnti(const auto& particleVector, const auto& particleVectorAnti, const auto& bufferAnti, auto& tempBuffer, const TVector3 jetAxis)
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
      for (int j = 0; j < static_cast<int>(particleVectorAnti.size()); j++) {
        if (std::isnan(particleVectorAnti.at(j).phi()))
          continue;
        if (TMath::Abs(particleVectorAnti.at(j).phi()) > 2 * TMath::Pi()) {
          registryData.fill(HIST("hTrackProtocol"), 15);
          continue;
        }

        double DeltaPhi = TVector2::Phi_0_2pi(particleVector.at(i).phi() - particleVectorAnti.at(j).phi());
        double DeltaEta = particleVector.at(i).eta() - particleVectorAnti.at(j).eta();
        if (DeltaPhi > (1.5 * TMath::Pi())) {
          DeltaPhi = DeltaPhi - 2 * TMath::Pi();
        }
        registryData.fill(HIST("hDeltaPhiSEProtonAntiproton"), DeltaPhi);
        registryData.fill(HIST("hDeltaPhiEtaSEProtonAntiproton"), DeltaPhi, DeltaEta);
        break;
      }
      fillMixedEventDeltas(particleVector.at(i), bufferAnti, 4, jetAxis);
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
    registryData.fill(HIST("hPtTotalSubJetPerp"), subtractedJetPerp.pt());
    registryData.fill(HIST("hPtTotalSubJetArea"), subtractedJetArea.pt());
    registryQA.fill(HIST("hRhoEstimateArea"), jet.pt(), rho); // switch to subtracted jet pt
    registryQA.fill(HIST("hRhoMEstimateArea"), jet.pt(), rhoM);
    registryQA.fill(HIST("hRhoEstimatePerp"), jet.pt(), rhoPerp);
    registryQA.fill(HIST("hRhoMEstimatePerp"), jet.pt(), rhoMPerp);
    double jetBkgDeltaPt = jet.pt() - subtractedJetPerp.pt();
    registryQA.fill(HIST("hJetBkgDeltaPt"), jetBkgDeltaPt);

    registryData.fill(HIST("hEventProtocol"), 4);
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    registryData.fill(HIST("hEventProtocol"), 5);
    registryData.fill(HIST("hNumberOfJets"), 0);
    registryData.fill(HIST("hPtTotalJet"), jet.pt());
    registryData.fill(HIST("hJetRapidity"), jet.rap());
    registryData.fill(HIST("hNumPartInJet"), jet.constituents().size());
    registryQA.fill(HIST("hJetPtVsNumPart"), jet.pt(), jet.constituents().size());

    double maxRadius = 0;
    for (const auto& constituent : constituents) {
      registryData.fill(HIST("hPtJetParticle"), constituent.pt());
      registryQA.fill(HIST("hPhiJet"), constituent.phi());
      registryQA.fill(HIST("hPhiPtJet"), constituent.pt(), constituent.phi());
      registryQA.fill(HIST("hEtaJet"), constituent.eta());
      registryQA.fill(HIST("hEtaPtJet"), constituent.pt(), constituent.eta());

      if (std::isnan(constituent.phi()) || std::isnan(jet.phi())) // geometric jet cone
        continue;
      double DeltaPhi = TVector2::Phi_0_2pi(constituent.phi() - jet.phi());
      if (DeltaPhi > TMath::Pi())
        DeltaPhi = DeltaPhi - 2 * TMath::Pi();
      double DeltaEta = constituent.eta() - jet.eta();
      double Delta = TMath::Sqrt(DeltaPhi * DeltaPhi + DeltaEta * DeltaEta);
      registryQA.fill(HIST("hJetConeRadius"), Delta);
      if (Delta > maxRadius)
        maxRadius = Delta;
    }
    registryQA.fill(HIST("hMaxRadiusVsPt"), jet.pt(), maxRadius);

    // QA for comparison with nuclei_in_jets
    TVector3 pJet(0., 0., 0.);
    pJet.SetXYZ(jet.px(), jet.py(), jet.pz());
    TVector3 UEAxis1(0.0, 0.0, 0.0);
    TVector3 UEAxis2(0.0, 0.0, 0.0);
    getPerpendicularAxis(pJet, UEAxis1, +1.0);
    getPerpendicularAxis(pJet, UEAxis2, -1.0);

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

      if (deltaRJet < Rmax) {
        if (deltaPhiJet != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiJet"), deltaEtaJet, deltaPhiJet);
        NchJetPlusUE++;
        ptJetPlusUE = ptJetPlusUE + track.pt();
      }
      if (deltaRUE1 < Rmax) {
        if (deltaPhiUE1 != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiUE"), deltaEtaUE1, deltaPhiUE1);
        NchUE++;
        ptUE = ptUE + track.pt();
      }
      if (deltaRUE2 < Rmax) {
        if (deltaPhiUE2 != -999)
          registryQA.fill(HIST("hDeltaEtadeltaPhiUE"), deltaEtaUE2, deltaPhiUE2);
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
    registryQA.fill(HIST("hDeltaJetPt"), jet.pt() - ptJetPlusUE);

    int nPartClusteredJet = static_cast<int>(constituents.size());

    // Fill QA Histograms
    if (ptJetPlusUE < minJetPt) { // swap for sub pt?

      registryQA.fill(HIST("hNParticlesClusteredInJet"), nPartClusteredJet);

      for (const auto& track : constituents) {
        registryQA.fill(HIST("hPtParticlesClusteredInJet"), track.pt());
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

    for (int i = 0; i < static_cast<int>(constituents.size()); i++) { // analyse jet constituents - this is where the magic happens
      registryData.fill(HIST("hTrackProtocol"), 3);
      fastjet::PseudoJet pseudoParticle = constituents.at(i);
      int id = pseudoParticle.user_index();
      const auto& jetParticle = particles.at(id);
      jetAll.emplace_back(jetParticle);

      registryData.fill(HIST("hDCAxyFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaXY());
      registryData.fill(HIST("hDCAzFullJet"), jetParticle.pt() * jetParticle.sign(), jetParticle.dcaZ());
      registryData.fill(HIST("hTPCsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.tpcSignal());
      if (jetParticle.hasTOF()) {
        registryData.fill(HIST("hTOFsignal"), jetParticle.pt() * jetParticle.sign(), jetParticle.beta());
      }
      double ptDiff = pseudoParticle.pt() - jetParticle.pt();
      registryQA.fill(HIST("hPtDiff"), ptDiff);

      if (jetParticle.pt() < minJetParticlePt)
        continue;
      if (isProton(jetParticle, false)) { // collect protons in jet
        registryData.fill(HIST("hPtJetProton"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetProtonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 4); // # protons
        if (isProton(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 5); // # high purity protons
          jetProtons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetProton"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isAntiproton(jetParticle, false)) { // collect antiprotons in jet
        registryData.fill(HIST("hPtJetAntiproton"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetAntiprotonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 6); // # antiprotons
        if (isAntiproton(jetParticle, true)) {
          registryData.fill(HIST("hTrackProtocol"), 7); // # high purity antiprotons
          jetAntiprotons.emplace_back(jetParticle);
          registryData.fill(HIST("hDCAzJetAntiproton"), jetParticle.pt(), jetParticle.dcaZ());
        }
      } else if (isNucleus(jetParticle)) { // collect nuclei in jet
        registryData.fill(HIST("hPtJetNuclei"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetNucleiVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 8); // # nuclei
        registryData.fill(HIST("hDCAzJetNuclei"), jetParticle.pt(), jetParticle.dcaZ());
      } else if (isAntinucleus(jetParticle)) {
        registryData.fill(HIST("hPtJetAntinuclei"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetAntinucleiVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 10); // # antinuclei
        registryData.fill(HIST("hDCAzJetAntinuclei"), jetParticle.pt(), jetParticle.dcaZ());
      } else if (isPion(jetParticle)) {
        registryData.fill(HIST("hPtJetPion"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetPionVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 11); // # antinuclei
        registryData.fill(HIST("hDCAzJetPion"), jetParticle.pt(), jetParticle.dcaZ());
        if (jetParticle.sign() > 0) {
          jetPiPlus.emplace_back(jetParticle);
        } else if (jetParticle.sign() < 0) {
          jetPiMinus.emplace_back(jetParticle);
        }
      } else if (isKaon(jetParticle)) {
        registryData.fill(HIST("hPtJetKaon"), jetParticle.pt());
        registryQA.fill(HIST("hPtJetKaonVsTotalJet"), jetParticle.pt(), subtractedJetPerp.pt());
        registryData.fill(HIST("hTrackProtocol"), 12); // # antinuclei
        registryData.fill(HIST("hDCAzJetKaon"), jetParticle.pt(), jetParticle.dcaZ());
      }
    } // for (int i=0; i<static_cast<int>(constituents.size()); i++)

    if (jetAll.size() > 1) { // general correlation function
      doCorrelations(jetAll, fBufferJet, fTempBufferJet, 0, pJet);
      setTrackBuffer(fTempBufferJet, fBufferJet);
    }
    jetCounter++;

    if ((jetProtons.size() > 0) && (jetAntiprotons.size() > 0)) {
      doCorrelationsAnti(jetProtons, jetAntiprotons, fBufferAntiproton, fTempBufferProton, pJet);
      doCorrelationsAnti(jetAntiprotons, jetProtons, fBufferProton, fTempBufferAntiproton, pJet); // divide SE distributions by 2 in post
    }
    if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2))
      return jetCounter;
    registryData.fill(HIST("hEventProtocol"), 6);

    if (jetProtons.size() > 1) {
      doCorrelations(jetProtons, fBufferProton, fTempBufferProton, 1, pJet);
      setTrackBuffer(fTempBufferProton, fBufferProton);
    }
    if (jetAntiprotons.size() > 1) {
      doCorrelations(jetAntiprotons, fBufferAntiproton, fTempBufferAntiproton, 2, pJet);
      setTrackBuffer(fTempBufferAntiproton, fBufferAntiproton);
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
      if (!selectTrack(track))
        continue;

      double mass;
      if (useTOFmass) {
        if (track.hasTOF()) {
          mass = track.mass(); // check reliability, maybe use only pion mass
          registryQA.fill(HIST("hTOFmass"), track.pt(), track.mass());
          registryData.fill(HIST("hTrackProtocol"), 1);
        } else {
          mass = 0.139; // pion mass as default, ~80% are pions
          registryData.fill(HIST("hTrackProtocol"), 2);
        }
      } else {
        mass = 0.139;
      }

      if (track.tpcNClsFindable() != 0) {
        registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
      }
      registryQA.fill(HIST("hPtFullEvent"), track.pt());
      registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
      registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
      registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
      registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
      registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
      registryQA.fill(HIST("hPhiFullEvent"), track.phi());
      registryQA.fill(HIST("hPhiPtFullEvent"), track.pt(), track.phi());
      registryQA.fill(HIST("hEtaFullEvent"), track.eta());
      registryQA.fill(HIST("hEtaPtFullEvent"), track.pt(), track.eta());

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
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap, ghost_repeat, ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());

    if (jets.size() == 0)
      return;

    registryData.fill(HIST("hEventProtocol"), 3);

    bool doSparse = true;
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(jetInput, doSparse);
    auto [rhoPerp, rhoMPerp] = bkgSub.estimateRhoPerpCone(jetInput, jets);

    for (const auto& jet : jets) {
      jetCounter = analyseJet(jetCounter, jet, particles, jetProtons, jetAntiprotons, jetPiPlus, jetPiMinus, jetAll, rho, rhoM, rhoPerp, rhoMPerp);
    }
    registryData.fill(HIST("hNumJetsInEvent"), jetCounter);

    TVector3 hardestJetAxis(jets.at(0).px(), jets.at(0).py(), jets.at(0).pz()); // for full event, use hardest jet as orientation
    doCorrelations(particlesForCF, fBufferFull, fTempBufferFull, -1, hardestJetAxis);
    setTrackBuffer(fTempBufferFull, fBufferFull);
  }

  template <typename U>
  void fillHistogramsMC(U const& tracks)
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
      if (!selectTrack(track))
        continue;

      double mass;
      if (useTOFmass) {
        if (track.hasTOF()) {
          mass = track.mass(); // check reliability, maybe use only pion mass
          registryQA.fill(HIST("hTOFmass"), track.pt(), track.mass());
          registryData.fill(HIST("hTrackProtocol"), 1);
        } else {
          mass = 0.139; // pion mass as default, ~80% are pions
          registryData.fill(HIST("hTrackProtocol"), 2);
        }
      } else {
        mass = 0.139;
      }

      if (track.tpcNClsFindable() != 0) {
        registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
      }
      registryQA.fill(HIST("hPtFullEvent"), track.pt());
      registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
      registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
      registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
      registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
      registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
      registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
      registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
      registryQA.fill(HIST("hPhiFullEvent"), track.phi());
      registryQA.fill(HIST("hPhiPtFullEvent"), track.pt(), track.phi());
      registryQA.fill(HIST("hEtaFullEvent"), track.eta());
      registryQA.fill(HIST("hEtaPtFullEvent"), track.pt(), track.eta());

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
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, jetR);
    fastjet::JetDefinition jetDefBkg(fastjet::kt_algorithm, 0.5);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(ghost_maxrap, ghost_repeat, ghost_area));
    fastjet::AreaDefinition areaDefBkg(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(ghost_maxrap));
    fastjet::ClusterSequenceArea clusterSeq(jetInput, jetDef, areaDef);
    jets = sorted_by_pt(clusterSeq.inclusive_jets());

    if (jets.size() == 0)
      return;

    registryData.fill(HIST("hEventProtocol"), 3);

    bool doSparse = true;
    auto [rho, rhoM] = bkgSub.estimateRhoAreaMedian(jetInput, doSparse);
    auto [rhoPerp, rhoMPerp] = bkgSub.estimateRhoPerpCone(jetInput, jets);

    for (auto& jet : jets) {
      jetCounter = analyseJet(jetCounter, jet, particles, jetProtons, jetAntiprotons, jetPiPlus, jetPiMinus, jetAll, rho, rhoM, rhoPerp, rhoMPerp);

      // MC Truth Particles
      fastjet::PseudoJet subtractedJetPerp(0., 0., 0., 0.);
      subtractedJetPerp = bkgSub.doRhoAreaSub(jet, rhoPerp, rhoMPerp);
      if (subtractedJetPerp.pt() < minJetPt) // cut on jet w/o bkg
        continue;
      for (const auto& constituent : jet.constituents()) {
        const auto& jetParticle = particles.at(constituent.user_index());
        if (!jetParticle.has_mcParticle())
          continue;
        switch (jetParticle.mcParticle().pdgCode()) {
          case 2212:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 0);
            registryMC.fill(HIST("hPtJetProtonMC"), jetParticle.pt());
            break;
          case -2212:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 1);
            registryMC.fill(HIST("hPtJetAntiprotonMC"), jetParticle.pt());
            break;
          case 1000010020:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 2);
            registryMC.fill(HIST("hPtJetNucleiMC"), jetParticle.pt());
            break;
          case -1000010020:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 3);
            registryMC.fill(HIST("hPtJetAntinucleiMC"), jetParticle.pt());
            break;
          case 1000020030:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 4);
            registryMC.fill(HIST("hPtJetNucleiMC"), jetParticle.pt());
            break;
          case -1000020030:
            registryMC.fill(HIST("hNumberOfTruthParticles"), 5);
            registryMC.fill(HIST("hPtJetAntinucleiMC"), jetParticle.pt());
            break;
          default:
            continue;
        }
      } // for (const auto& constituent : jet.constituents())
    } // for (const auto& jet : jets)
    registryData.fill(HIST("hNumJetsInEvent"), jetCounter);

    TVector3 hardestJetAxis(jets.at(0).px(), jets.at(0).py(), jets.at(0).pz()); // for full event, use hardest jet as orientation
    doCorrelations(particlesForCF, fBufferFull, fTempBufferFull, -1, hardestJetAxis);
    setTrackBuffer(fTempBufferFull, fBufferFull);
  }

  void processRun2(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun2> const& tracks,
                   BCsWithRun2Info const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      registryData.fill(HIST("hEventProtocol"), 0);
      if (!collision.alias_bit(kINT7))
        continue;
      registryData.fill(HIST("hNumberOfEvents"), 0);
      registryData.fill(HIST("hEventProtocol"), 1);

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun2, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun2, "process Run 2 data", false);

  void processRun3(soa::Join<aod::Collisions, aod::EvSels> const& collisions,
                   soa::Filtered<FullTracksRun3> const& tracks)
  {
    for (const auto& collision : collisions) {
      registryData.fill(HIST("hEventProtocol"), 0);
      if (!collision.sel8())
        continue;
      registryData.fill(HIST("hNumberOfEvents"), 0);
      registryData.fill(HIST("hEventProtocol"), 1);
      if (TMath::Abs(collision.posZ()) > zVtx)
        continue;

      auto slicedTracks = tracks.sliceBy(perCollisionFullTracksRun3, collision.globalIndex());

      fillHistograms(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun3, "process Run 3 data", false);

  // using JetTracksMCDwID = soa::Join<aod::JetTracksMCD, aod::JTrackExtras, aod::JTrackPIs>;
  void processRun3revised(soa::Filtered<soa::Join<aod::JetCollisions, aod::JCollisionPIs>>::iterator const& collision, soa::Filtered<soa::Join<aod::ChargedJets, aod::ChargedJetConstituents>> const& allJets, /* JTracksRun3 const& jtracks, soa::Join<aod::Collisions, aod::EvSels> const&, */ soa::Filtered<FullTracksRun3> const&) // check how to use bkg sub jets --- table or recluster + bkg sub?
  {
    registryData.fill(HIST("hEventProtocol"), 0);
    if (!jetderiveddatautilities::selectCollision(collision, eventSelection))
      return registryData.fill(HIST("hNumberOfEvents"), 0);
    registryData.fill(HIST("hEventProtocol"), 1);

    int jetCounter = 0;

    for (const auto& jet : allJets) { // loop over jets in event
      jetCounter++;
      std::vector<FullTracksRun3::iterator> jetProtons;
      std::vector<FullTracksRun3::iterator> jetAntiprotons;
      std::vector<FullTracksRun3::iterator> jetPiPlus;
      std::vector<FullTracksRun3::iterator> jetPiMinus;
      std::vector<FullTracksRun3::iterator> jetAll;
      std::vector<std::pair<double, double>> fTempBufferProton;
      std::vector<std::pair<double, double>> fTempBufferAntiproton;
      std::vector<std::pair<double, double>> fTempBufferPiPlus;
      std::vector<std::pair<double, double>> fTempBufferPiMinus;
      std::vector<std::pair<double, double>> fTempBufferJet;
      jetProtons.clear();
      jetAntiprotons.clear();
      jetPiPlus.clear();
      jetPiMinus.clear();
      jetAll.clear();
      fTempBufferProton.clear();
      fTempBufferAntiproton.clear();
      fTempBufferPiPlus.clear();
      fTempBufferPiMinus.clear();
      fTempBufferJet.clear();
      TVector3 pJet(0., 0., 0.);
      pJet.SetXYZ(jet.px(), jet.py(), jet.pz());

      registryData.fill(HIST("hNumberOfJets"), 0);
      registryData.fill(HIST("hPtTotalJet"), jet.pt());
      registryData.fill(HIST("hJetRapidity"), jet.eta());
      registryData.fill(HIST("hNumPartInJet"), jet.tracksIds().size());
      registryQA.fill(HIST("hJetPtVsNumPart"), jet.pt(), jet.tracksIds().size());
      registryQA.fill(HIST("hMaxRadiusVsPt"), jet.pt(), jet.r());

      for (const auto& track : jet.template tracks_as<FullTracksRun3>()) { // slice on jets?
        if (!selectTrack(track))
          continue;

        if (track.tpcNClsFindable() != 0) {
          registryQA.fill(HIST("hRatioCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows() / track.tpcNClsFindable());
        }
        registryQA.fill(HIST("hPtJetParticle"), track.pt());
        registryQA.fill(HIST("hCrossedRowsTPC"), track.pt(), track.tpcNClsCrossedRows());
        registryQA.fill(HIST("hClusterITS"), track.pt(), track.itsNCls());
        registryQA.fill(HIST("hClusterTPC"), track.pt(), track.tpcNClsFound());
        registryQA.fill(HIST("hChi2ITS"), track.pt(), track.itsChi2NCl());
        registryQA.fill(HIST("hChi2TPC"), track.pt(), track.tpcChi2NCl());
        registryQA.fill(HIST("hDCAxyFullEvent"), track.pt(), track.dcaXY());
        registryQA.fill(HIST("hDCAzFullEvent"), track.pt(), track.dcaZ());
        registryQA.fill(HIST("hPhiJet"), track.phi());
        registryQA.fill(HIST("hPhiPtJet"), track.pt(), track.phi());
        registryQA.fill(HIST("hEtaJet"), track.eta());
        registryQA.fill(HIST("hEtaPtJet"), track.pt(), track.eta());

        if (!std::isnan(track.phi()) && !std::isnan(jet.phi())) { // geometric jet cone
          double DeltaPhi = TVector2::Phi_0_2pi(track.phi() - jet.phi());
          if (DeltaPhi > TMath::Pi())
            DeltaPhi = DeltaPhi - 2 * TMath::Pi();
          double DeltaEta = track.eta() - jet.eta();
          double Delta = TMath::Sqrt(DeltaPhi * DeltaPhi + DeltaEta * DeltaEta);
          registryQA.fill(HIST("hJetConeRadius"), Delta);
        }

        // analyse jet constituents - this is where the magic happens
        registryData.fill(HIST("hTrackProtocol"), 3);
        jetAll.emplace_back(track);

        registryData.fill(HIST("hDCAxyFullJet"), track.pt() * track.sign(), track.dcaXY());
        registryData.fill(HIST("hDCAzFullJet"), track.pt() * track.sign(), track.dcaZ());
        registryData.fill(HIST("hTPCsignal"), track.pt() * track.sign(), track.tpcSignal());
        if (track.hasTOF()) {
          registryData.fill(HIST("hTOFsignal"), track.pt() * track.sign(), track.beta());
        }
        // double ptDiff = pseudoParticle.pt() - track.pt();
        // registryQA.fill(HIST("hPtDiff"), ptDiff);

        if (isProton(track, false)) { // collect protons in jet
          registryData.fill(HIST("hPtJetProton"), track.pt());
          registryQA.fill(HIST("hPtJetProtonVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 4); // # protons
          if (isProton(track, true)) {
            registryData.fill(HIST("hTrackProtocol"), 5); // # high purity protons
            jetProtons.emplace_back(track);
            registryData.fill(HIST("hDCAzJetProton"), track.pt(), track.dcaZ());
          }
        } else if (isAntiproton(track, false)) { // collect antiprotons in jet
          registryData.fill(HIST("hPtJetAntiproton"), track.pt());
          registryQA.fill(HIST("hPtJetAntiprotonVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 6); // # antiprotons
          if (isAntiproton(track, true)) {
            registryData.fill(HIST("hTrackProtocol"), 7); // # high purity antiprotons
            jetAntiprotons.emplace_back(track);
            registryData.fill(HIST("hDCAzJetAntiproton"), track.pt(), track.dcaZ());
          }
        } else if (isNucleus(track)) { // collect nuclei in jet
          registryData.fill(HIST("hPtJetNuclei"), track.pt());
          registryQA.fill(HIST("hPtJetNucleiVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 8); // # nuclei
          registryData.fill(HIST("hDCAzJetNuclei"), track.pt(), track.dcaZ());
        } else if (isAntinucleus(track)) {
          registryData.fill(HIST("hPtJetAntinuclei"), track.pt());
          registryQA.fill(HIST("hPtJetAntinucleiVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 10); // # antinuclei
          registryData.fill(HIST("hDCAzJetAntinuclei"), track.pt(), track.dcaZ());
        } else if (isPion(track)) {
          registryData.fill(HIST("hPtJetPion"), track.pt());
          registryQA.fill(HIST("hPtJetPionVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 11); // # antinuclei
          registryData.fill(HIST("hDCAzJetPion"), track.pt(), track.dcaZ());
          if (track.sign() > 0) {
            jetPiPlus.emplace_back(track);
          } else if (track.sign() < 0) {
            jetPiMinus.emplace_back(track);
          }
        } else if (isKaon(track)) {
          registryData.fill(HIST("hPtJetKaon"), track.pt());
          registryQA.fill(HIST("hPtJetKaonVsTotalJet"), track.pt(), jet.pt());
          registryData.fill(HIST("hTrackProtocol"), 12); // # antinuclei
          registryData.fill(HIST("hDCAzJetKaon"), track.pt(), track.dcaZ());
        }
      } // for (const auto& jtrack : jtracks)

      if (jetAll.size() > 1) { // general correlation function
        doCorrelations(jetAll, fBufferJet, fTempBufferJet, 0, pJet);
        setTrackBuffer(fTempBufferJet, fBufferJet);
      }

      if ((jetProtons.size() > 0) && (jetAntiprotons.size() > 0)) {
        doCorrelationsAnti(jetProtons, jetAntiprotons, fBufferAntiproton, fTempBufferProton, pJet);
        doCorrelationsAnti(jetAntiprotons, jetProtons, fBufferProton, fTempBufferAntiproton, pJet); // divide SE distributions by 2 in post
      }
      if ((jetProtons.size() < 2) && (jetAntiprotons.size() < 2) && jetPiPlus.size() < 2 && jetPiMinus.size() < 2)
        continue;
      registryData.fill(HIST("hEventProtocol"), 6);

      if (jetProtons.size() > 1) {
        doCorrelations(jetProtons, fBufferProton, fTempBufferProton, 1, pJet);
        setTrackBuffer(fTempBufferProton, fBufferProton);
      }
      if (jetAntiprotons.size() > 1) {
        doCorrelations(jetAntiprotons, fBufferAntiproton, fTempBufferAntiproton, 2, pJet);
        setTrackBuffer(fTempBufferAntiproton, fBufferAntiproton);
      }
      if (jetPiPlus.size() > 1) {
        doCorrelations(jetPiPlus, fBufferPiPlus, fTempBufferPiPlus, 1, pJet);
        setTrackBuffer(fTempBufferPiPlus, fBufferPiPlus);
      }
      if (jetPiMinus.size() > 1) {
        doCorrelations(jetPiMinus, fBufferPiMinus, fTempBufferPiMinus, 1, pJet);
        setTrackBuffer(fTempBufferPiMinus, fBufferPiMinus);
      }
    } // for (const auto& jet : allJets)

    registryData.fill(HIST("hNumJetsInEvent"), jetCounter);
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processRun3revised, "process Run 3 data w jet tables", true);

  void processMCRun2(McCollisions const& collisions, soa::Filtered<McTracksRun2> const& tracks, BCsWithRun2Info const&, aod::McParticles&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<BCsWithRun2Info>();
      initCCDB(bc);

      registryData.fill(HIST("hEventProtocol"), 0);
      if (!collision.alias_bit(kINT7))
        continue;
      registryData.fill(HIST("hNumberOfEvents"), 0);
      registryData.fill(HIST("hEventProtocol"), 1);

      auto slicedTracks = tracks.sliceBy(perCollisionMcTracksRun2, collision.globalIndex());

      fillHistogramsMC(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processMCRun2, "process Run 2 MC", false);

  void processMCRun3(McCollisions const& collisions, soa::Filtered<McTracksRun3> const& tracks, aod::McParticles&, aod::McCollisions const&)
  {
    for (const auto& collision : collisions) {
      registryData.fill(HIST("hEventProtocol"), 0);
      if (!collision.sel8())
        continue;
      registryData.fill(HIST("hNumberOfEvents"), 0);
      registryData.fill(HIST("hEventProtocol"), 1);
      if (TMath::Abs(collision.posZ()) > zVtx)
        continue;

      auto slicedTracks = tracks.sliceBy(perCollisionMcTracksRun3, collision.globalIndex());

      fillHistogramsMC(slicedTracks);
    }
  }
  PROCESS_SWITCH(AngularCorrelationsInJets, processMCRun3, "process Run 3 MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<AngularCorrelationsInJets>(cfgc)};
}
