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

/// \file taskElectronWeakBoson.cxx
/// \brief task for WeakBoson (W/Z) based on electron in mid-rapidity
/// \author S. Sakai & S. Ito (Univ. of Tsukuba)

#ifndef HomogeneousField
#define HomogeneousField // o2-linter: disable=name/macro (required by KFParticle)
#endif

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGJE/DataModel/EMCALClusters.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/KFparticle/KFUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <KFPTrack.h>
#include <KFParticle.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_centrality;

struct HfTaskElectronWeakBoson {

  // configurable parameters
  Configurable<int> nBinsPt{"nBinsPt", 100, "N bins in pt registry"};
  Configurable<float> binPtmax{"binPtmax", 100.0, "maximum pt registry"};
  Configurable<int> nBinsE{"nBinsE", 100, "N bins in E registry"};
  Configurable<float> binEmax{"binEmax", 100.0, "maximum E registry"};

  Configurable<float> vtxZ{"vtxZ", 10.f, ""};

  Configurable<float> etaTrMin{"etaTrMin", -1.0f, "minimun track eta"};
  Configurable<float> etaTrMax{"etaTrMax", 1.0f, "maximum track eta"};
  Configurable<float> etaEmcMax{"etaEmcMax", 0.6f, "maximum track eta"};
  Configurable<float> dcaxyMax{"dcaxyMax", 2.0f, "mximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> ptMin{"ptMin", 3.0f, "minimum pT cut"};
  Configurable<float> ptAssMin{"ptAssMin", 0.15, "minimum pT cut for associated hadrons"};
  Configurable<float> ptMatch{"ptMatch", 0.001, "pT match in Z->ee and associated tracks"};
  Configurable<float> ptZeeMin{"ptZeeMin", 20.0f, "minimum pT cut for Zee"};
  Configurable<float> chi2TpcMax{"chi2TpcMax", 4.0f, "tpc chi2 cut"};
  Configurable<float> nclItsMin{"nclItsMin", 2.0f, "its # of cluster cut"};
  Configurable<float> nclTpcMin{"nclTpcMin", 100.0f, "tpc # if cluster cut"};
  Configurable<float> nclcrossTpcMin{"nclcrossTpcMin", 100.0f, "tpc # of crossedRows cut"};
  Configurable<float> nsigTpcMinLose{"nsigTpcMinLose", -3.0, "tpc Nsig lose lower cut"};
  Configurable<float> nsigTpcMin{"nsigTpcMin", -1.0, "tpc Nsig lower cut"};
  Configurable<float> nsigTpcMax{"nsigTpcMax", 3.0, "tpc Nsig upper cut"};

  Configurable<float> phiEmcMin{"phiEmcMin", 1.39, "EMC phi acc min"};
  Configurable<float> phiEmcMax{"phiEmcMax", 3.36, "EMC phi acc max"};
  Configurable<int> clusterDefinition{"clusterDefinition", 10, "cluster definition to be selected, e.g. 10=kV3Default"};
  Configurable<float> timeEmcMin{"timeEmcMin", -25., "Minimum EMCcluster timing"};
  Configurable<float> timeEmcMax{"timeEmcMax", +20., "Maximum EMCcluster timing"};
  Configurable<float> m02Min{"m02Min", 0.1, "Minimum M02"};
  Configurable<float> m02Max{"m02Max", 0.9, "Maximum M02"};
  Configurable<float> rMatchMax{"rMatchMax", 0.05, "cluster - track matching cut"};
  Configurable<float> eopMin{"eopMin", 0.9, "Minimum eop"};
  Configurable<float> eopMax{"eopMax", 1.3, "Maximum eop"};

  Configurable<float> rIsolation{"rIsolation", 0.3, "cone radius for isolation cut"};
  Configurable<float> energyIsolationMax{"energyIsolationMax", 0.1, "isolation cut on energy"};
  Configurable<float> momentumIsolationMax{"momentumIsolationMax", 0.1, "isolation cut on momentum"};
  Configurable<int> trackIsolationMax{"trackIsolationMax", 3, "Maximum number of tracks in isolation cone"};

  Configurable<float> massZMin{"massZMin", 60.0, "Minimum Z mass (GeV/c^2)"};
  Configurable<float> massZMax{"massZMax", 120.0, "Maximum Z mass (GeV/c^2)"};
  Configurable<float> correctionPtElectron{"correctionPtElectron", 1.0, "momentum correction factor for decay electrons from Z boson"};

  // flag for THn
  Configurable<bool> isTHnElectron{"isTHnElectron", true, "Enables THn for electrons"};
  Configurable<float> ptTHnThresh{"ptTHnThresh", 5.0, "Threshold for THn make"};

  // Skimmed (trigger) dataset processing configurations
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "Enables processing of skimmed datasets"};
  Configurable<std::string> cfgTriggerName{"cfgTriggerName", "fGammaHighPtEMCAL", "Trigger of interest (comma separated for multiple)"};
  Configurable<bool> applySel8{"applySel8", true, "Apply sel8 filter or not"};

  // CCDB service configurations
  Configurable<std::string> cfgCCDBPath{"cfgCCDBPath", "Users/m/mpuccio/EventFiltering/OTS/", "Path to CCDB for trigger data"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // KFParticle
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  Configurable<int> chiSqNdfMax{"chiSqNdfMax", 10, "Chi2 Max for mass reco by KF particle"};

  // Centrality estimator configuration
  Configurable<int> centralityEstimator{"centralityEstimator", CentralityEstimator::FT0M, "Centrality estimator. See CentralityEstimator for valid values."};
  Configurable<bool> enableCentralityAnalysis{"enableCentralityAnalysis", true, "Enable centrality-dependent analysis"};
  Configurable<float> centralityMin{"centralityMin", -1, "minimum cut on centrality selection"};
  Configurable<float> centralityMax{"centralityMax", 101, "maximum cut on centrality selection"};
  Configurable<std::vector<double>> centralityBins{"centralityBins", {0, 20, 60, 100}, "centrality bins"};

  // QA for Z->ee
  Configurable<bool> enableZeeRecoQA{"enableZeeRecoQA", false, "Enable QA for Z->ee reconstruction"};
  Configurable<float> massZMinQA{"massZMinQA", 0.1, "minimum mass cut for Zee Reco QA"};
  // CCDB service object
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct HfElectronCandidate {
    float pt, eta, phi, dcaxyTrk, dcazTrk, eop, energyIso, momIso;
    int ntrackIso, nclusterTPC, nclusterITS;
    HfElectronCandidate(float ptr, float e, float ph, float dcaxy, float dcaz, float ep, float eiso, float piso, int ntrkiso, int nclstpc, int nclsits)
      : pt(ptr), eta(e), phi(ph), dcaxyTrk(dcaxy), dcazTrk(dcaz), eop(ep), energyIso(eiso), momIso(piso), ntrackIso(ntrkiso), nclusterTPC(nclstpc), nclusterITS(nclsits) {}
  };
  std::vector<HfElectronCandidate> selectedElectronsIso;
  std::vector<HfElectronCandidate> selectedPositronsIso;
  std::vector<HfElectronCandidate> selectedElectronsAss;

  struct HfZeeCandidate {
    float pt, eta, phi, mass, ptchild0, ptchild1;
    int charge;
    HfZeeCandidate(float ptr, float e, float ph, float m, int ch, float ptzee0, float ptzee1)
      : pt(ptr), eta(e), phi(ph), mass(m), ptchild0(ptzee0), ptchild1(ptzee1), charge(ch) {}
  };
  std::vector<HfZeeCandidate> reconstructedZ;
  using CollisionsWithCent = soa::Join<aod::Collisions, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
  using SelectedClusters = o2::aod::EMCALClusters;
  // PbPb
  // using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>;
  using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::FullTracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>;

  // pp
  // using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCEl, o2::aod::pidTOFEl>>;

  Filter eventFilter = (applySel8 ? (o2::aod::evsel::sel8 == true) : (o2::aod::evsel::sel8 == o2::aod::evsel::sel8));
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < vtxZ);

  Filter etafilter = (aod::track::eta < etaTrMax) && (aod::track::eta > etaTrMin);
  Filter dcaxyfilter = (nabs(aod::track::dcaXY) < dcaxyMax);
  Filter filterGlobalTr = requireGlobalTrackInFilter();

  Filter clusterDefinitionSelection = (o2::aod::emcalcluster::definition == clusterDefinition) && (o2::aod::emcalcluster::time >= timeEmcMin) && (o2::aod::emcalcluster::time <= timeEmcMax) && (o2::aod::emcalcluster::m02 > m02Min) && (o2::aod::emcalcluster::m02 < m02Max);

  // Data Handling Objects
  Preslice<o2::aod::EMCALClusterCells> perCluster = o2::aod::emcalclustercell::emcalclusterId;
  Preslice<o2::aod::EMCALAmbiguousClusterCells> perClusterAmb = o2::aod::emcalclustercell::emcalambiguousclusterId;
  PresliceUnsorted<o2::aod::EMCALMatchedTracks> perClusterMatchedTracks = o2::aod::emcalmatchedtrack::trackId;

  // config axis
  ConfigurableAxis confaxisPt{"confaxisPt", {100, 0, 100}, "p_{T}"};
  ConfigurableAxis confaxisPtZee{"confaxisPtZee", {60, 20, 80}, "p_{T}"};
  ConfigurableAxis confaxisPtZneg{"confaxisPtZneg", {60, 20, 80}, "p_{T,neg} (GeV/c)"};
  ConfigurableAxis confaxisPtZpos{"confaxisPtZpos", {60, 20, 80}, "p_{T,pos} (GeV/c)"};
  ConfigurableAxis confaxisEop{"confaxisEop", {300, -0.01, 1.49}, "E/p"};
  ConfigurableAxis confaxisEopZneg{"confaxisEopZneg", {300, -0.01, 1.49}, "E/p_{neg}"};
  ConfigurableAxis confaxisEopZpos{"confaxisEopZpos", {300, -0.01, 1.49}, "E/p_{pos}"};
  ConfigurableAxis confaxisIsoEnergy{"confaxisIsoEnergy", {255, 0, 2.0}, "E_{iso}"};
  ConfigurableAxis confaxisIsoEnergyZneg{"confaxisIsoEnergyZneg", {255, 0, 2.0}, "E_{iso,neg}"};
  ConfigurableAxis confaxisIsoEnergyZpos{"confaxisIsoEnergyZpos", {255, 0, 2.0}, "E_{iso,pos}"};
  ConfigurableAxis confaxisIsoMomentum{"confaxisIsoMomentum", {255, 0, 2.0}, "E_{iso}"};
  ConfigurableAxis confaxisIsoMomentumZneg{"confaxisIsoMomentumZneg", {255, 0, 1.5}, "E_{iso,neg}"};
  ConfigurableAxis confaxisIsoMomentumZpos{"confaxisIsoMomentumZpos", {255, 0, 1.5}, "E_{iso,pos}"};
  ConfigurableAxis confaxisIsoTrack{"confaxisIsoTrack", {25, -0.5, 24.5}, "Isolation Track"};
  ConfigurableAxis confaxisIsoTrackZneg{"confaxisIsoTrackZneg", {25, -0.5, 24.5}, "N_{isotrk,neg}"};
  ConfigurableAxis confaxisIsoTrackZpos{"confaxisIsoTrackZpos", {25, -0.5, 24.5}, "N_{isotrk,pos}"};
  ConfigurableAxis confaxisInvMassZgamma{"confaxisInvMassZgamma", {150, 0, 150}, "M_{ee} (GeV/c^{2})"};
  ConfigurableAxis confaxisInvMassZ{"confaxisInvMassZ", {130, 20, 150}, "M_{ee} (GeV/c^{2})"};
  ConfigurableAxis confaxisZfrag{"confaxisZfrag", {200, 0, 2.0}, "p_{T,h}/p_{T,Z}"};

  // Histogram registry: an object to hold your registrygrams
  HistogramRegistry registry{"registry"};

  // Zorro objects for skimmed data processing
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  void init(InitContext const&)
  {
    // Configure CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // CCDB path for debug
    LOGF(info, "CCDB path for Zorro: %s", cfgCCDBPath.value.c_str());

    // Setup Zorro Summary
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }
    // check centrality
    if (centralityEstimator < CentralityEstimator::FT0A || centralityEstimator > CentralityEstimator::FV0A) {
      LOGF(fatal, "Invalid centrality estimator: %d", static_cast<int>(centralityEstimator.value));
    }

    // add configurable for CCDB path
    zorro.setBaseCCDBPath(cfgCCDBPath.value);

    // define axes you want to use
    const AxisSpec axisZvtx{40, -20, 20, "Zvtx"};
    const AxisSpec axisCounter{1, 0, 1, "events"};
    const AxisSpec axisEta{20, -1.0, 1.0, "#eta"};
    const AxisSpec axisDCAxyneg{150, 0, 0.3, "DCAxy_{neg}"};
    const AxisSpec axisDCAxypos{150, 0, 0.3, "DCAxy_{pos}"};
    const AxisSpec axisDCAzneg{150, 0, 0.3, "DCAz_{neg}"};
    const AxisSpec axisDCAzpos{150, 0, 0.3, "DCAz_{neg}"};
    const AxisSpec axisNclsTPCneg{20, 79.5, 159.5, "nClsTpc_{neg}"};
    const AxisSpec axisNclsTPCpos{20, 79.5, 159.5, "nClsTpc_{neg}"};
    const AxisSpec axisNclsITSneg{9, -0.5, 8.5, "nClsIts_{neg}"};
    const AxisSpec axisNclsITSpos{9, -0.5, 8.5, "nClsIts_{neg}"};
    const AxisSpec axisSectorTPCneg{360, 0, 18, "TPCsector_{neg}"};
    const AxisSpec axisSectorTPCpos{360, 0, 18, "TPCsector_{pos}"};
    const AxisSpec axisNsigma{100, -5, 5, "N#sigma"};
    const AxisSpec axisDedx{150, 0, 150, "dEdx"};
    const AxisSpec axisE{nBinsE, 0, binEmax, "Energy"};
    const AxisSpec axisM02{100, 0, 1, "M02"};
    const AxisSpec axisdPhi{100, -0.5, 0.5, "dPhi"};
    const AxisSpec axisdEta{100, -0.5, 0.5, "dEta"};
    const AxisSpec axisdR{20, 0.0, 0.2, "dR"};
    const AxisSpec axisNcell{50, 0.0, 50.0, "Ncell"};
    const AxisSpec axisPhi{350, 0, 7, "Phi"};
    const AxisSpec axisChi2{250, 0.0, 25.0, "#chi^{2}"};
    const AxisSpec axisCluster{100, 0.0, 200.0, "counts"};
    const AxisSpec axisITSNCls{10, 0.0, 10, "counts"};
    const AxisSpec axisEMCtime{100, -50.0, 50, "EMC time"};
    const AxisSpec axisTrigger{3, -0.5, 2.5, "Trigger status of zorro"};
    const AxisSpec axisDPhiZh{64, -o2::constants::math::PIHalf, 3 * o2::constants::math::PIHalf, "#Delta#phi(Z-h)"};
    const AxisSpec axisPtHadron{nBinsPt, 0, binPtmax, "p_{T,hadron} (GeV/c)"};
    const AxisSpec axisPtZ{nBinsPt, 0, binPtmax, "p_{T,Z} (GeV/c)"};
    const AxisSpec axisSign{2, -2, 2, "charge sign"};
    const AxisSpec axisCentrality{centralityBins, "centrality"};
    const AxisSpec axisEop{confaxisEop, "E/p"};
    const AxisSpec axisEopZneg{confaxisEopZneg, "E/p_{neg}"};
    const AxisSpec axisEopZpos{confaxisEopZpos, "E/p_{[pos}"};
    const AxisSpec axisPt{confaxisPt, "p_{T}"};
    const AxisSpec axisPtZee{confaxisPtZee, "p_{T}"};
    const AxisSpec axisPtZneg{confaxisPtZneg, "p_{T,neg} (GeV/c)"};
    const AxisSpec axisPtZpos{confaxisPtZpos, "p_{T,pos} (GeV/c)"};
    const AxisSpec axisIsoEnergy{confaxisIsoEnergy, "E_{iso}"};
    const AxisSpec axisIsoEnergyZneg{confaxisIsoEnergyZneg, "E_{iso}"};
    const AxisSpec axisIsoEnergyZpos{confaxisIsoEnergyZpos, "E_{iso}"};
    const AxisSpec axisIsoMomentum{confaxisIsoMomentum, "Isolation momentum(GeV/C)"};
    const AxisSpec axisIsoMomentumZneg{confaxisIsoMomentumZneg, "Isolation momentum(GeV/C)"};
    const AxisSpec axisIsoMomentumZpos{confaxisIsoMomentumZpos, "Isolation momentum(GeV/C)"};
    const AxisSpec axisIsoTrack{confaxisIsoTrack, "Isolation Track"};
    const AxisSpec axisIsoTrackZneg{confaxisIsoTrackZneg, "N_{isotrk,neg}"};
    const AxisSpec axisIsoTrackZpos{confaxisIsoTrackZpos, "N_{isotrk,pos}"};
    const AxisSpec axisInvMassZgamma{confaxisInvMassZgamma, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisInvMassZ{confaxisInvMassZ, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisZfrag{confaxisZfrag, "p_{T,h}/p_{T,Z}"};

    // create registrygrams
    registry.add("hZvtx", "Z vertex", kTH1D, {axisZvtx});
    registry.add("hEventCounterInit", "hEventCounterInit", kTH1D, {axisCounter});
    registry.add("hEventCounter", "hEventCounter", kTH1D, {axisCounter});
    registry.add("hCentrality", "Centrality distribution", kTH1D, {axisCentrality});
    registry.add("hITSchi2", "ITS #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCchi2", "TPC #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCnCls", "TPC NCls", kTH1F, {axisCluster});
    registry.add("hITSnCls", "ITS NCls", kTH1F, {axisITSNCls});
    registry.add("hTPCnClsCrossedRows", "TPC NClsCrossedRows", kTH1F, {axisCluster});
    registry.add("hEta", "track eta", kTH1F, {axisEta});
    registry.add("hPt", "track pt", kTH1F, {axisPt});
    registry.add("hTPCNsigma", "TPC electron Nsigma", kTH2F, {{axisPt}, {axisNsigma}});
    registry.add("hEnergy", "EMC cluster energy", kTH1F, {axisE});
    registry.add("hEnergyNcell", "EMC cluster energy and cell", kTH2F, {{axisE}, {axisNcell}});
    registry.add("hTrMatchR", "Track EMC Match in radius", kTH2F, {{axisPt}, {axisdR}});
    registry.add("hTrMatch_mim", "Track EMC Match minimu minimumm", kTH2F, {{axisdPhi}, {axisdEta}});
    registry.add("hMatchPhi", "Match in Phi", kTH2F, {{axisPhi}, {axisPhi}});
    registry.add("hMatchEta", "Match in Eta", kTH2F, {{axisEta}, {axisEta}});
    registry.add("hEop", "energy momentum match", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopIsolation", "energy momentum match after isolation", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopIsolationTr", "energy momentum match after isolationTr", kTH2F, {{axisPt}, {axisEop}});
    registry.add("hEopNsigTPC", "Eop vs. Nsigma", kTH2F, {{axisNsigma}, {axisEop}});
    registry.add("hEMCtime", "EMC timing", kTH1F, {axisEMCtime});
    registry.add("hIsolationEnergy", "Isolation Energy", kTH2F, {{axisE}, {axisIsoEnergy}});
    registry.add("hInvMassZee", "invariant mass for Z ULS pair", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPt, axisInvMassZgamma});
    registry.add("hKfInvMassZee", "invariant mass for Z ULS pair KFp", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPt, axisInvMassZgamma});
    registry.add("hInvMassZeeQA", "QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisDCAxyneg, axisDCAxypos, axisDCAzpos, axisNclsTPCneg, axisNclsTPCpos, axisNclsITSneg, axisNclsITSpos, axisSectorTPCneg, axisSectorTPCneg, axisEopZneg, axisEopZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hInvMassZeeQAbg", "QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisDCAxyneg, axisDCAxypos, axisDCAzpos, axisNclsTPCneg, axisNclsTPCpos, axisNclsITSneg, axisNclsITSpos, axisSectorTPCneg, axisSectorTPCneg, axisEopZneg, axisEopZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hTHnElectrons", "electron info", HistType::kTHnSparseF, {axisPt, axisNsigma, axisM02, axisEop, axisIsoEnergy, axisIsoMomentum, axisIsoTrack, axisEta, axisDedx});
    registry.add("hTHnTrMatch", "Track EMC Match", HistType::kTHnSparseF, {axisPt, axisdPhi, axisdEta});

    // Z-hadron correlation histograms
    registry.add("hZHadronDphi", "Z-hadron #Delta#phi correlation", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPtZ, axisDPhiZh, axisZfrag, axisPtHadron});
    registry.add("hZptSpectrum", "Z boson p_{T} spectrum", kTH2F, {{axisSign}, {axisPtZ}});

    // hisotgram for EMCal trigger
    registry.add("hEMCalTrigger", "EMCal trigger", kTH1D, {axisTrigger});
  }

  double getIsolatedCluster(const o2::aod::EMCALCluster& cluster,
                            const SelectedClusters& clusters)
  {
    double energySum = 0.0;
    double isoEnergy = 10.0;
    double const etaAssCluster = cluster.eta();
    double const phiAssCluster = cluster.phi();

    for (const auto& associateCluster : clusters) {
      // Calculate angular distances
      double const dEta = associateCluster.eta() - etaAssCluster;
      double dPhi = associateCluster.phi() - phiAssCluster;

      // Normalize φ difference
      dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

      // Calculate ΔR
      double const deltaR = std::sqrt(dEta * dEta + dPhi * dPhi);

      // Sum energy within isolation cone
      if (deltaR < rIsolation) {
        energySum += associateCluster.energy();
      }
    }

    if (energySum > 0) {
      isoEnergy = energySum / cluster.energy() - 1.0;
    }

    registry.fill(HIST("hIsolationEnergy"), cluster.energy(), isoEnergy);

    return (isoEnergy);
  }
  std::pair<int, double> getIsolatedTrack(double etaEle,
                                          double phiEle,
                                          float pEle,
                                          TrackEle const& tracks)
  {
    int trackCount = 0;
    double isoMomentum = 10;
    double pSum = 0.0;
    // LOG(info) << "track p = " << pEle;

    for (const auto& track : tracks) {

      double const dEta = track.eta() - etaEle;
      double dPhi = track.phi() - phiEle;
      dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

      double const deltaR = std::sqrt(dEta * dEta + dPhi * dPhi);

      if (deltaR < rIsolation) {
        trackCount++;
        pSum += track.p();
      }
    }

    // LOG(info) << "momSun = " << pSum;
    if (pSum > 0) {
      isoMomentum = pSum / pEle - 1.0;
    }

    // LOG(info) << "isop = " << isoMomentum;
    return std::make_pair(trackCount - 1, isoMomentum);
  }

  void recoMassZee(const KFParticle& kfpIsoEle,
                   int charge,
                   float centrality,
                   TrackEle const& tracks)
  {
    // LOG(info) << "Invarimass cal by KF particle ";
    for (const auto& track : tracks) {

      if (std::abs(track.pt() - kfpIsoEle.GetPt()) < ptMatch) {
        continue;
      }
      if (track.pt() < ptZeeMin) {
        continue;
      }
      if (std::abs(track.tpcNSigmaEl()) > nsigTpcMax) {
        continue;
      }
      if (std::abs(track.eta()) > etaTrMax) {
        continue;
      }
      int pdgAss = kElectron;
      if (track.sign() > 0) {
        pdgAss = kPositron;
      }
      auto [trackCount, isoMomentum] = getIsolatedTrack(track.eta(), track.phi(), track.p(), tracks);
      if (isoMomentum > momentumIsolationMax) {
        continue;
      }

      KFPTrack const kfpTrackAssEle = createKFPTrackFromTrack(track);
      KFParticle const kfpAssEle(kfpTrackAssEle, pdgAss);
      // reco by RecoDecay
      auto child1 = RecoDecayPtEtaPhi::pVector(kfpIsoEle.GetPt() * correctionPtElectron, kfpIsoEle.GetEta(), kfpIsoEle.GetPhi());
      auto child2 = RecoDecayPtEtaPhi::pVector(kfpAssEle.GetPt() * correctionPtElectron, kfpAssEle.GetEta(), kfpAssEle.GetPhi());
      double const invMassEE = RecoDecay::m(std::array{child1, child2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});

      registry.fill(HIST("hInvMassZee"), centrality, track.sign() * charge, kfpIsoEle.GetPt(), invMassEE);

      // reco by KFparticle
      const KFParticle* electronPairs[2] = {&kfpIsoEle, &kfpAssEle};
      KFParticle zeeKF;
      zeeKF.SetConstructMethod(kfConstructMethod);
      zeeKF.Construct(electronPairs, 2);
      // LOG(info) << "Invarimass cal by KF particle Chi2/NDF = " << zeeKF.GetChi2()/zeeKF.GetNDF();
      float const chiSqNdf = zeeKF.GetChi2() / zeeKF.GetNDF();
      if (zeeKF.GetNDF() < 1) {
        continue;
      }
      if (zeeKF.GetChi2() < 0) {
        continue;
      }
      if (chiSqNdf > chiSqNdfMax) {
        continue;
      }
      float massZee, massZeeErr;
      zeeKF.GetMass(massZee, massZeeErr);
      registry.fill(HIST("hKfInvMassZee"), centrality, track.sign() * charge, kfpIsoEle.GetPt(), massZee);
      // LOG(info) << "Invarimass cal by KF particle mass = " << massZee;
      // LOG(info) << "Invarimass cal by RecoDecay = " << invMassEE;
      reconstructedZ.emplace_back(
        zeeKF.GetPt(),
        zeeKF.GetEta(),
        zeeKF.GetPhi(),
        massZee,
        track.sign() * charge,
        kfpIsoEle.GetPt(),
        kfpAssEle.GetPt());
    }
  }

  // void process(soa::Filtered<aod::Collisions>::iterator const& collision,
  void process(soa::Filtered<CollisionsWithCent>::iterator const& collision,
               aod::BCsWithTimestamps const&,
               SelectedClusters const& emcClusters,
               TrackEle const& tracks,
               o2::aod::EMCALMatchedTracks const& matchedtracks)
  {
    registry.fill(HIST("hEventCounterInit"), 0.5);

    // Get BC for this collision
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    uint64_t const globalBC = bc.globalBC();
    int const runNumber = bc.runNumber();

    // Initialize Zorro for the first event (once per run)
    static bool isFirstEvent = true;
    static int lastRunNumber = -1;

    if ((isFirstEvent || runNumber != lastRunNumber) && cfgSkimmedProcessing) {
      LOGF(info, "Initializing Zorro for run %d", runNumber);
      uint64_t const currentTimestamp = bc.timestamp();

      // debug for timestamp
      LOGF(info, "Using CCDB path: %s, timestamp: %llu", cfgCCDBPath.value.c_str(), currentTimestamp);

      // initialize Zorro
      zorro.initCCDB(ccdb.service, runNumber, currentTimestamp, cfgTriggerName);
      isFirstEvent = false;
      lastRunNumber = runNumber;

      // initialize magnetic field
      auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrpMag, currentTimestamp);
      o2::base::Propagator::initFieldFromGRP(grpo);
      double const magneticField = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << "magneticField = " << magneticField;
      if (magneticField != 0.0) {
        KFParticle::SetField(magneticField);
      }
    }

    // Check if this is a triggered event using Zorro
    bool isTriggered = true;
    if (cfgSkimmedProcessing) {
      isTriggered = zorro.isSelected(globalBC);
      registry.fill(HIST("hEMCalTrigger"), isTriggered ? 1 : 0);

      // Skip event if not triggered and we're processing skimmed data
      if (!isTriggered) {
        return;
      }
    }
    // initialze for inclusive-electron
    selectedElectronsIso.clear();
    selectedPositronsIso.clear();
    selectedElectronsAss.clear();
    reconstructedZ.clear();

    registry.fill(HIST("hEventCounter"), 0.5);

    // LOGF(info, "Collision index : %d", collision.index());
    // LOGF(info, "Number of tracks: %d", tracks.size());
    // LOGF(info, "Number of clusters: %d", clusters.size());

    registry.fill(HIST("hZvtx"), collision.posZ());

    // Calculate centrality
    float centrality = 1.0;
    if (enableCentralityAnalysis) {
      centrality = o2::hf_centrality::getCentralityColl(collision, centralityEstimator);
      // LOG(info) << centrality;
      if (centrality < centralityMin || centrality > centralityMax) {
        return;
      }
      registry.fill(HIST("hCentrality"), centrality);
    }

    for (const auto& track : tracks) {

      if (std::abs(track.eta()) > etaTrMax) {
        continue;
      }
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin) {
        continue;
      }
      if (std::abs(track.dcaXY()) > dcaxyMax) {
        continue;
      }
      if (track.itsChi2NCl() > chi2ItsMax) {
        continue;
      }
      if (track.tpcChi2NCl() > chi2TpcMax) {
        continue;
      }
      if (track.tpcNClsFound() < nclTpcMin) {
        continue;
      }
      if (track.itsNCls() < nclItsMin) {
        continue;
      }

      registry.fill(HIST("hEta"), track.eta());
      registry.fill(HIST("hITSchi2"), track.itsChi2NCl());
      registry.fill(HIST("hTPCchi2"), track.tpcChi2NCl());
      registry.fill(HIST("hTPCnCls"), track.tpcNClsFound());
      registry.fill(HIST("hITSnCls"), track.itsNCls());
      registry.fill(HIST("hTPCnClsCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("hPt"), track.pt());
      registry.fill(HIST("hTPCNsigma"), track.p(), track.tpcNSigmaEl());

      float eop = -0.01;
      float isoEnergy = 1.0;
      // track isolation
      auto [trackCount, isoMomentum] = getIsolatedTrack(track.eta(), track.phi(), track.p(), tracks);
      // LOG(info) << "isoMomentum = " << isoMomentum;

      if (track.pt() > ptAssMin) {
        selectedElectronsAss.emplace_back(
          track.pt(),
          track.eta(),
          track.phi(),
          track.dcaXY(),
          track.dcaZ(),
          eop,
          isoEnergy,
          isoMomentum,
          trackCount,
          track.tpcNClsCrossedRows(),
          track.itsNCls());
      }

      if (track.pt() < ptMin) {
        continue;
      }

      // LOG(info) << "tr phi, eta = " << track.phi() << " ; " << track.eta();
      // EMC acc
      bool isEMCacceptance = true;
      if (track.phi() < phiEmcMin || track.phi() > phiEmcMax) {
        isEMCacceptance = false;
      }
      if (std::abs(track.eta()) > etaEmcMax) {
        isEMCacceptance = false;
      }
      // LOG(info) << "EMC acc  = " << isEMCacceptance;
      auto tracksofcluster = matchedtracks.sliceBy(perClusterMatchedTracks, track.globalIndex());

      double rMin = 999.9;
      double dPhiMin = 999.9;
      double dEtaMin = 999.9;
      bool isIsolated = false;
      bool isIsolatedTr = false;

      if ((tracksofcluster.size() != 0) && isEMCacceptance) {
        int nMatch = 0;
        for (const auto& match : tracksofcluster) {
          if (match.emcalcluster_as<SelectedClusters>().time() < timeEmcMin || match.emcalcluster_as<SelectedClusters>().time() > timeEmcMax) {
            continue;
          }

          float const m02Emc = match.emcalcluster_as<SelectedClusters>().m02();
          float const energyEmc = match.emcalcluster_as<SelectedClusters>().energy();
          double const phiEmc = match.emcalcluster_as<SelectedClusters>().phi();
          double const etaEmc = match.emcalcluster_as<SelectedClusters>().eta();
          double const timeEmc = match.emcalcluster_as<SelectedClusters>().time();
          // LOG(info) << "tr phi0 = " << match.track_as<TrackEle>().phi();
          // LOG(info) << "tr phi1 = " << track.phi();
          // LOG(info) << "emc phi = " << phiEmc;

          if (nMatch == 0) {
            double const dEta = match.track_as<TrackEle>().trackEtaEmcal() - etaEmc;
            double dPhi = match.track_as<TrackEle>().trackPhiEmcal() - phiEmc;
            dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);

            registry.fill(HIST("hMatchPhi"), phiEmc, match.track_as<TrackEle>().trackPhiEmcal());
            registry.fill(HIST("hMatchEta"), etaEmc, match.track_as<TrackEle>().trackEtaEmcal());

            double const r = RecoDecay::sqrtSumOfSquares(dPhi, dEta);
            // LOG(info) << "r match = " << r;
            if (r < rMin) {
              rMin = r;
              dPhiMin = dPhi;
              dEtaMin = dEta;
            }
            registry.fill(HIST("hTHnTrMatch"), match.track_as<TrackEle>().pt(), dPhi, dEta);
            registry.fill(HIST("hEMCtime"), timeEmc);
            registry.fill(HIST("hEnergy"), energyEmc);

            if (std::abs(dPhi) > rMatchMax || std::abs(dEta) > rMatchMax) {
              continue;
            }

            registry.fill(HIST("hTrMatchR"), match.track_as<TrackEle>().pt(), r);
            registry.fill(HIST("hEnergyNcell"), energyEmc, match.emcalcluster_as<SelectedClusters>().nCells());

            const auto& cluster = match.emcalcluster_as<SelectedClusters>();

            eop = energyEmc / match.track_as<TrackEle>().p();
            // LOG(info) << "eop = " << eop;

            isoEnergy = getIsolatedCluster(cluster, emcClusters);

            if (match.track_as<TrackEle>().pt() > ptTHnThresh && isTHnElectron) {
              registry.fill(HIST("hTHnElectrons"), match.track_as<TrackEle>().pt(), match.track_as<TrackEle>().tpcNSigmaEl(), m02Emc, eop, isoEnergy, isoMomentum, trackCount, track.eta(), track.tpcSignal());
            }
            // LOG(info) << "E/p" << eop;
            registry.fill(HIST("hEopNsigTPC"), match.track_as<TrackEle>().tpcNSigmaEl(), eop);
            if (match.emcalcluster_as<SelectedClusters>().m02() < m02Min || match.emcalcluster_as<SelectedClusters>().m02() > m02Max) {
              continue;
            }

            if (match.track_as<TrackEle>().tpcNSigmaEl() > nsigTpcMin && match.track_as<TrackEle>().tpcNSigmaEl() < nsigTpcMax) {
              registry.fill(HIST("hEop"), match.track_as<TrackEle>().pt(), eop);
              if (eop > eopMin && eop < eopMax && isoEnergy < energyIsolationMax) {
                isIsolated = true;
              }
              if (eop > eopMin && eop < eopMax && trackCount < trackIsolationMax && isoMomentum < momentumIsolationMax) {
                isIsolatedTr = true;
              }

              if (isIsolated && isIsolatedTr) {
                registry.fill(HIST("hEopIsolation"), match.track_as<TrackEle>().pt(), eop);

                if (match.track_as<TrackEle>().pt() > ptZeeMin) {
                  int pdgIso = kElectron;
                  if (match.track_as<TrackEle>().sign() > 0) {
                    pdgIso = kPositron;
                  }
                  KFPTrack const kfpTrackIsoEle = createKFPTrackFromTrack(match.track_as<TrackEle>());
                  KFParticle const kfpIsoEle(kfpTrackIsoEle, pdgIso);
                  recoMassZee(kfpIsoEle, match.track_as<TrackEle>().sign(), centrality, tracks);

                } // end of pt cut for e from Z
              } // end if isolation cut
              if (isIsolatedTr) {
                registry.fill(HIST("hEopIsolationTr"), match.track_as<TrackEle>().pt(), eop);
              }
            } // end of PID cut
          } // end of nmatch == 0
          nMatch++;
        } // end of cluster match
      } // end of cluster

      if (rMin < rMatchMax) {
        // LOG(info) << "R mim = " << rMin;
        registry.fill(HIST("hTrMatch_mim"), dPhiMin, dEtaMin);
      }
      if (enableZeeRecoQA && track.pt() > ptZeeMin) {
        if (track.sign() < 0) {
          selectedElectronsIso.emplace_back(
            track.pt(),
            track.eta(),
            track.phi(),
            track.dcaXY(),
            track.dcaZ(),
            eop,
            isoEnergy,
            isoMomentum,
            trackCount,
            track.tpcNClsFound(),
            track.itsNCls());
        } else {
          selectedPositronsIso.emplace_back(
            track.pt(),
            track.eta(),
            track.phi(),
            track.dcaXY(),
            track.dcaZ(),
            eop,
            isoEnergy,
            isoMomentum,
            trackCount,
            track.tpcNClsFound(),
            track.itsNCls());
        }
      }

    } // end of track loop
    // Z-hadron
    if (!reconstructedZ.empty()) {
      for (const auto& zBoson : reconstructedZ) {
        // Z boson selection
        if (zBoson.mass < massZMin || zBoson.mass > massZMax) {
          continue;
        }
        registry.fill(HIST("hZptSpectrum"), zBoson.charge, zBoson.pt);
        for (const auto& trackAss : selectedElectronsAss) {
          if (std::abs(trackAss.pt - zBoson.ptchild0) < ptMatch) {
            continue;
          }
          if (std::abs(trackAss.pt - zBoson.ptchild1) < ptMatch) {
            continue;
          }
          // calculate Z-h correlation
          double const deltaPhi = RecoDecay::constrainAngle(trackAss.phi - zBoson.phi, -o2::constants::math::PIHalf);
          double const ptRatio = trackAss.pt / zBoson.pt;
          registry.fill(HIST("hZHadronDphi"), centrality, zBoson.charge, zBoson.pt, deltaPhi, ptRatio, trackAss.pt);
        }
      }
    } // end of Z-hadron correlation
    // Z->ee QA
    if (enableZeeRecoQA) {
      if (!selectedElectronsIso.empty() && !selectedPositronsIso.empty()) {
        // signal
        for (const auto& trackEle : selectedElectronsIso) {
          auto child1 = RecoDecayPtEtaPhi::pVector(trackEle.pt, trackEle.eta, trackEle.phi);
          float const sectorneg = trackEle.phi / o2::constants::math::SectorSpanRad;
          for (const auto& trackPos : selectedPositronsIso) {
            auto child2 = RecoDecayPtEtaPhi::pVector(trackPos.pt, trackPos.eta, trackPos.phi);
            double const invMass = RecoDecay::m(std::array{child1, child2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
            float const sectorpos = trackPos.phi / o2::constants::math::SectorSpanRad;
            if (invMass > massZMinQA) {
              registry.fill(HIST("hInvMassZeeQA"), invMass, trackEle.pt, trackPos.pt, trackEle.dcaxyTrk, trackPos.dcaxyTrk, trackPos.dcazTrk, trackEle.nclusterTPC, trackPos.nclusterTPC, trackEle.nclusterITS, trackPos.nclusterITS, sectorneg, sectorpos, trackEle.eop, trackPos.eop, trackEle.energyIso, trackPos.energyIso, trackEle.momIso, trackPos.momIso, trackEle.ntrackIso, trackPos.ntrackIso);
            }
          }
        }
        // bg e-e-
        for (size_t i = 0; i < selectedElectronsIso.size(); ++i) {
          const auto& trackEle = selectedElectronsIso[i];
          auto child1 = RecoDecayPtEtaPhi::pVector(trackEle.pt, trackEle.eta, trackEle.phi);
          float const sectorneg = trackEle.phi / o2::constants::math::SectorSpanRad;
          for (size_t j = i + 1; j < selectedElectronsIso.size(); ++j) {
            const auto& trackEle2 = selectedElectronsIso[j];
            auto child2 = RecoDecayPtEtaPhi::pVector(trackEle2.pt, trackEle2.eta, trackEle2.phi);
            float const sectorpos = trackEle2.phi / o2::constants::math::SectorSpanRad;
            double const invMass = RecoDecay::m(std::array{child1, child2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
            if (invMass > massZMinQA) {
              registry.fill(HIST("hInvMassZeeQAbg"), invMass, trackEle.pt, trackEle2.pt, trackEle.dcaxyTrk, trackEle2.dcaxyTrk, trackEle2.dcazTrk, trackEle.nclusterTPC, trackEle2.nclusterTPC, trackEle.nclusterITS, trackEle2.nclusterITS, sectorneg, sectorpos, trackEle.eop, trackEle2.eop, trackEle.energyIso, trackEle2.energyIso, trackEle.momIso, trackEle2.momIso, trackEle.ntrackIso, trackEle2.ntrackIso);
            }
          }
        }
        // bg e+e+
        for (size_t i = 0; i < selectedPositronsIso.size(); ++i) {
          const auto& trackPos = selectedPositronsIso[i];
          auto child1 = RecoDecayPtEtaPhi::pVector(trackPos.pt, trackPos.eta, trackPos.phi);
          float const sectorneg = trackPos.phi / o2::constants::math::SectorSpanRad;
          for (size_t j = i + 1; j < selectedPositronsIso.size(); ++j) {
            const auto& trackPos2 = selectedPositronsIso[j];
            auto child2 = RecoDecayPtEtaPhi::pVector(trackPos2.pt, trackPos2.eta, trackPos2.phi);
            float const sectorpos = trackPos2.phi / o2::constants::math::SectorSpanRad;
            double const invMass = RecoDecay::m(std::array{child1, child2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
            if (invMass > massZMinQA) {
              registry.fill(HIST("hInvMassZeeQAbg"), invMass, trackPos.pt, trackPos2.pt, trackPos.dcaxyTrk, trackPos2.dcaxyTrk, trackPos2.dcazTrk, trackPos.nclusterTPC, trackPos2.nclusterTPC, trackPos.nclusterITS, trackPos2.nclusterITS, sectorneg, sectorpos, trackPos.eop, trackPos2.eop, trackPos.energyIso, trackPos2.energyIso, trackPos.momIso, trackPos2.momIso, trackPos.ntrackIso, trackPos2.ntrackIso);
            }
          }
        }
      } // check arrays aren't empty
    } // end of Z->ee QA
  } // process
}; // struct HfTaskElectronWeakBoson

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskElectronWeakBoson>(cfgc)};
}
