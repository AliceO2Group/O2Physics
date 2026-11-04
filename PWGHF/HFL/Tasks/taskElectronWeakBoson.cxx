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
#include "Common/DataModel/Multiplicity.h"
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
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>
#include <TParticlePDG.h>
#include <TRandom3.h>

#include <KFPTrack.h>
#include <KFParticle.h>

#include <algorithm>
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

  Configurable<bool> applyTrackSys{"applyTrackSys", false, "apply additional track cuts for systematic study"};
  Configurable<float> etaTrMin{"etaTrMin", -1.0f, "minimun track eta"};
  Configurable<float> etaTrMax{"etaTrMax", 1.0f, "maximum track eta"};
  Configurable<float> etaEmcMax{"etaEmcMax", 0.6f, "maximum track eta"};
  Configurable<float> dcaxyMax{"dcaxyMax", 2.0f, "mximum DCA xy"};
  Configurable<float> chi2ItsMax{"chi2ItsMax", 15.0f, "its chi2 cut"};
  Configurable<float> ptMin{"ptMin", 3.0f, "minimum pT cut"};
  Configurable<float> ptAssMin{"ptAssMin", 0.15, "minimum pT cut for associated hadrons"};
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
  // Zorro's own default (Common/Core/Zorro.h: mBaseCCDBPath = "EventFiltering/Zorro/").
  // The old value "Users/m/mpuccio/EventFiltering/OTS/" is a personal area from early Run 3
  // and no longer holds FilterCounters, so Zorro::initCCDB aborts the job on it.
  Configurable<std::string> cfgCCDBPath{"cfgCCDBPath", "EventFiltering/Zorro/", "Path to CCDB for trigger data"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};

  // KFParticle
  Configurable<int> kfConstructMethod{"kfConstructMethod", 2, "KF Construct Method"};
  Configurable<int> chiSqNdfMax{"chiSqNdfMax", 10, "Chi2 Max for mass reco by KF particle"};

  // Centrality estimator configuration
  Configurable<int> centralityEstimator{"centralityEstimator", CentralityEstimator::FT0M, "Centrality estimator. See CentralityEstimator for valid values."};
  Configurable<bool> enableCentralityAnalysis{"enableCentralityAnalysis", true, "Enable centrality-dependent analysis"};
  Configurable<bool> enableMultiplicityPVAnalysis{"enableMultiplicityPVAnalysis", true, "Enable centrality-dependent analysis"};
  Configurable<bool> enableMultiplicityFT0MAnalysis{"enableMultiplicityFT0MAnalysis", true, "Enable centrality-dependent analysis"};
  Configurable<float> centralityMin{"centralityMin", -1, "minimum cut on centrality selection"};
  Configurable<float> centralityMax{"centralityMax", 101, "maximum cut on centrality selection"};
  Configurable<std::vector<double>> centralityBins{"centralityBins", {0, 20, 60, 100}, "centrality bins"};

  // QA for Z->ee.
  // Grouped in a ConfigurableGroup to stay below the O2 limit of 100 task members.
  struct : ConfigurableGroup {
    Configurable<bool> enableZeeTrkRecoQA{"enableZeeTrkRecoQA", false, "Enable QA for track info Z->ee reconstruction"};
    Configurable<bool> enableZeeEmcRecoQA{"enableZeeEmcRecoQA", false, "Enable QA for Emc info Z->ee reconstruction"};
    Configurable<float> massZMinQA{"massZMinQA", 0.1, "minimum mass cut for Zee Reco QA"};
  } zeeQaCfg;

  // UE (same remark: grouped, no prefix, option names unchanged)
  struct : ConfigurableGroup {
    Configurable<int> nRandomCones{"nRandomCones", 100, "number of random cones"};
    Configurable<float> rcHardE{"rcHardE", 5.0, "hard cluster veto energy"};
    Configurable<float> rcVetoR{"rcVetoR", 0.4, "veto radius"};
    Configurable<bool> useUEsub{"useUEsub", true, "apply UE subtraction in isolation"};
  } ueCfg;

  // W->e - h correlation. Grouped so that the three options cost one task member instead of
  // three: the struct is close to the 99 the framework allows.
  struct : ConfigurableGroup {
    Configurable<bool> enableWeHadron{"enableWeHadron", false, "fill the isolated-electron - hadron correlation"};
    Configurable<float> ptWeMin{"ptWeMin", 20.0f, "minimum pT of the trigger electron"};
    Configurable<bool> vetoZeeElectron{"vetoZeeElectron", true, "reject trigger electrons that build a Z->ee candidate"};
  } weCfg;

  // Systematics by cut variation (method C): one extra THnSparse axis whose bin 0 is the
  // nominal selection and bins 1..N the variations.
  struct : ConfigurableGroup {
    Configurable<bool> enableSystematics{"enableSystematics", false, "add the cut-variation axis"};
    // associated hadron, checked against dcaxyTrk / nclusterTPC (crossed rows) / nclusterITS
    Configurable<std::vector<float>> varDcaxyMax{"varDcaxyMax", {}, "variations of dcaxyMax"};
    Configurable<std::vector<float>> varNclCrossTpcMin{"varNclCrossTpcMin", {}, "variations of nclcrossTpcMin"};
    Configurable<std::vector<float>> varNclItsMin{"varNclItsMin", {}, "variations of nclItsMin"};
    // W trigger electron identification
    Configurable<std::vector<float>> varEopMin{"varEopMin", {}, "variations of eopMin"};
    Configurable<std::vector<float>> varEopMax{"varEopMax", {}, "variations of eopMax"};
    Configurable<std::vector<float>> varNsigTpcMin{"varNsigTpcMin", {}, "variations of nsigTpcMin"};
    Configurable<std::vector<float>> varNsigTpcMax{"varNsigTpcMax", {}, "variations of nsigTpcMax"};
    Configurable<std::vector<float>> varM02Min{"varM02Min", {}, "variations of m02Min"};
    Configurable<std::vector<float>> varM02Max{"varM02Max", {}, "variations of m02Max"};
    // Z candidate. Two-sided: reconstructedZ holds every KF chi2-passing pair at any mass,
    // so a wider window really adds candidates instead of being a no-op.
    Configurable<std::vector<float>> varMassZMin{"varMassZMin", {}, "variations of massZMin"};
    Configurable<std::vector<float>> varMassZMax{"varMassZMax", {}, "variations of massZMax"};
  } sysCfg;

  // MC : efficiency of the associated hadrons in the Z-h correlation.
  struct : ConfigurableGroup {
    std::string prefix = "mcEfficiency"; // JSON group name
    Configurable<bool> requirePhysicalPrimary{"requirePhysicalPrimary", true, "require isPhysicalPrimary for the generated associated hadrons"};
    Configurable<bool> requireSameMcCollision{"requireSameMcCollision", true, "require the matched MC particle to belong to the MC collision of the reconstructed one"};
    ConfigurableAxis confaxisPtHadron{"confaxisPtHadron", {VARIABLE_WIDTH, 0.0, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 15.0, 20.0, 30.0, 50.0, 100.0}, "p_{T,hadron} (GeV/c)"};
    ConfigurableAxis confaxisEtaHadron{"confaxisEtaHadron", {20, -1.0, 1.0}, "#eta_{hadron}"};
  } mcCfg;

  // CCDB service object
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  // PDG database, used to get the charge of the generated particles
  Service<o2::framework::O2DatabasePDG> pdg{};

  // TParticlePDG::Charge() is given in units of |e|/3, a charged particle has |Charge()| >= 1
  static constexpr double ChargeTolerance = 0.001;

  struct HfElectronCandidate {
    float pt, eta, phi, dcaxyTrk, dcazTrk, eop, energyIso, momIso, nSigmaTpc, m02Emc;
    int ntrackIso, nclusterTPC, nclusterITS, sign;
    int64_t indexTrk; // globalIndex of the track, used to remove the Z daughters exactly
    HfElectronCandidate(float ptr, float e, float ph, float dcaxy, float dcaz, float ep, float eiso, float piso, float nsig, float m02, int ntrkiso, int nclstpc, int nclsits, int sgn, int64_t idxtrk)
      : pt(ptr), eta(e), phi(ph), dcaxyTrk(dcaxy), dcazTrk(dcaz), eop(ep), energyIso(eiso), momIso(piso), nSigmaTpc(nsig), m02Emc(m02), ntrackIso(ntrkiso), nclusterTPC(nclstpc), nclusterITS(nclsits), sign(sgn), indexTrk(idxtrk) {}
  };
  std::vector<HfElectronCandidate> selectedElectronsIso;
  std::vector<HfElectronCandidate> selectedPositronsIso;
  std::vector<HfElectronCandidate> selectedElectronsAss;

  struct HfZeeCandidate {
    float pt, eta, phi, mass, ptchild0, ptchild1;
    int charge;
    int64_t indexChild0, indexChild1; // globalIndex of the two tracks used to build the Z
    HfZeeCandidate(float ptr, float e, float ph, float m, int ch, float ptzee0, float ptzee1, int64_t idxzee0, int64_t idxzee1)
      : pt(ptr), eta(e), phi(ph), mass(m), ptchild0(ptzee0), ptchild1(ptzee1), charge(ch), indexChild0(idxzee0), indexChild1(idxzee1) {}
  };
  std::vector<HfZeeCandidate> reconstructedZ;
  using CollisionsWithCent = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::Mults>;
  using SelectedClusters = o2::aod::EMCALClusters;
  // PbPb
  // using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl>;
  using TrackEle = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::FullTracks, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi>;

  // MC : same collisions/tracks as above, plus the MC labels. aod::EvSels already comes with
  // CollisionsWithCent (joining it again would duplicate the table), so the data and the MC
  // paths see exactly the same event selection.
  using CollisionsWithCentMc = soa::Join<CollisionsWithCent, o2::aod::McCollisionLabels>;
  using TrackEleMc = o2::soa::Join<TrackEle, o2::aod::McTrackLabels>;

  // pp
  // using TrackEle = o2::soa::Filtered<o2::soa::Join<o2::aod::Tracks, o2::aod::FullTracks, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCEl, o2::aod::pidTOFEl>>;

  Filter eventFilter = ifnode(as<bool>(applySel8), o2::aod::evsel::sel8 == true, true);
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
  ConfigurableAxis confaxisIsoEnergy{"confaxisIsoEnergy", {255, 0, 2.0}, "E_{iso}"};
  ConfigurableAxis confaxisIsoMomentum{"confaxisIsoMomentum", {255, 0, 2.0}, "p_{iso}"};
  ConfigurableAxis confaxisIsoTrack{"confaxisIsoTrack", {25, -0.5, 24.5}, "Isolation Track"};
  ConfigurableAxis confaxisEtCluster{"confaxisEtCluster", {100, 0, 100}, "E_{T} (GeV)"};
  ConfigurableAxis confaxisNsigmaHadron{"confaxisNsigmaHadron", {200, -10, 10}, "n#sigma^{#pi}_{TPC}"};
  ConfigurableAxis confaxisInvMassZgamma{"confaxisInvMassZgamma", {150, 0, 150}, "M_{ee} (GeV/c^{2})"};
  ConfigurableAxis confaxisInvMassZ{"confaxisInvMassZ", {130, 20, 150}, "M_{ee} (GeV/c^{2})"};
  ConfigurableAxis confaxisZfrag{"confaxisZfrag", {200, 0, 2.0}, "p_{T,h}/p_{T,Z}"};
  ConfigurableAxis confaxisDEtaZh{"confaxisDEtaZh", {200, -4.0, 4.0}, "#Delta#eta"};
  ConfigurableAxis confaxisMultPV{"confaxisMultPV", {200, 0, 200.0}, "multiplicity"};
  ConfigurableAxis confaxisMultFT0{"confaxisMultFT0", {1000, 0, 1000.0}, "multiplicity"};

  // Histogram registry: an object to hold your registrygrams
  HistogramRegistry registry{"registry"};

  // Zorro objects for skimmed data processing
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // defined rnd
  TRandom3* rnd = nullptr;

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

    // init random
    rnd = new TRandom3(0);

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
    const AxisSpec axisNsigmaZneg{100, -5, 5, "N#sigma_{pos}"};
    const AxisSpec axisNsigmaZpos{100, -5, 5, "N#sigma_{neg}"};
    const AxisSpec axisEtCluster{confaxisEtCluster, "E_{T} (GeV)"};
    const AxisSpec axisNsigmaHadron{confaxisNsigmaHadron, "n#sigma^{#pi}_{TPC}"};
    const AxisSpec axisSignHadron{2, -2, 2, "charge sign (hadron)"};
    const AxisSpec axisNMatch{10, -0.5, 9.5, "N_{cluster} within rMatchMax"};
    const AxisSpec axisE{nBinsE, 0, binEmax, "Energy"};
    const AxisSpec axisM02{100, 0, 1, "M02"};
    const AxisSpec axisM02neg{100, 0, 1, "M02(neg)"};
    const AxisSpec axisM02pos{100, 0, 1, "M02(pos)"};
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
    const AxisSpec axisDEtaZh{confaxisDEtaZh, "#Delta#Eta(Z-h)"};
    const AxisSpec axisPtHadron{nBinsPt, 0, binPtmax, "p_{T,hadron} (GeV/c)"};
    const AxisSpec axisPtZ{nBinsPt, 0, binPtmax, "p_{T,Z} (GeV/c)"};
    const AxisSpec axisSign{2, -2, 2, "charge sign"};
    const AxisSpec axisCentrality{centralityBins, "centrality"};
    const AxisSpec axisEop{confaxisEop, "E/p"};
    const AxisSpec axisEopZneg{confaxisEop, "E/p_{neg}"};
    const AxisSpec axisEopZpos{confaxisEop, "E/p_{[pos}"};
    const AxisSpec axisPt{confaxisPt, "p_{T}"};
    const AxisSpec axisPtZee{confaxisPtZee, "p_{T}"};
    const AxisSpec axisPtZneg{confaxisPtZneg, "p_{T,neg} (GeV/c)"};
    const AxisSpec axisPtZpos{confaxisPtZpos, "p_{T,pos} (GeV/c)"};
    const AxisSpec axisIsoEnergy{confaxisIsoEnergy, "E_{iso}"};
    const AxisSpec axisIsoEnergyZneg{confaxisIsoEnergy, "E_{iso, neg} (GeV)"};
    const AxisSpec axisIsoEnergyZpos{confaxisIsoEnergy, "E_{iso, pos} (GeV)"};
    const AxisSpec axisIsoMomentum{confaxisIsoMomentum, "p_{iso} (GeV/C)"};
    const AxisSpec axisIsoMomentumZneg{confaxisIsoMomentum, "p_{iso, neg}(GeV/C)"};
    const AxisSpec axisIsoMomentumZpos{confaxisIsoMomentum, "p_{iso, pos}(GeV/C)"};
    const AxisSpec axisIsoTrack{confaxisIsoTrack, "Isolation Track"};
    const AxisSpec axisIsoTrackZneg{confaxisIsoTrack, "N_{isotrk,neg}"};
    const AxisSpec axisIsoTrackZpos{confaxisIsoTrack, "N_{isotrk,pos}"};
    const AxisSpec axisInvMassZgamma{confaxisInvMassZgamma, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisInvMassZ{confaxisInvMassZ, "M_{ee} (GeV/c^{2})"};
    const AxisSpec axisZfrag{confaxisZfrag, "p_{T,h}/p_{T,Z}"};
    const AxisSpec axisMultPV{confaxisMultPV, "multiplicity"};
    const AxisSpec axisMultFT0{confaxisMultFT0, "multiplicity"};

    // create registrygrams
    registry.add("hZvtx", "Z vertex", kTH1D, {axisZvtx});
    registry.add("hEventCounterInit", "hEventCounterInit", kTH1D, {axisCounter});
    registry.add("hEventCounter", "hEventCounter", kTH1D, {axisCounter});
    registry.add("hCentrality", "Centrality distribution", kTH1D, {axisCentrality});
    registry.add("hCentMultFT0Corr", "Centrality distribution vs. FT0 Mult", kTH2D, {{axisCentrality}, {axisMultFT0}});
    registry.add("hCentMultPVCorr", "Centrality distribution vs. PV Mult", kTH2D, {{axisCentrality}, {axisMultPV}});
    registry.add("hMultPV", "multiplicity  distribution for PV", kTH2D, {{axisZvtx}, {axisMultPV}});
    registry.add("hMultFT0", "multiplicity distribution for FT0", kTH2D, {{axisZvtx}, {axisMultFT0}});
    registry.add("hMultFT0PV", "multiplicity distribution", kTH2D, {{axisMultFT0}, {axisMultPV}});
    registry.add("hITSchi2", "ITS #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCchi2", "TPC #chi^{2}", kTH1F, {axisChi2});
    registry.add("hTPCnCls", "TPC NCls", kTH1F, {axisCluster});
    registry.add("hITSnCls", "ITS NCls", kTH1F, {axisITSNCls});
    registry.add("hTPCnClsCrossedRows", "TPC NClsCrossedRows", kTH1F, {axisCluster});
    registry.add("hEta", "track eta", kTH1F, {axisEta});
    registry.add("hPt", "track pt", kTH1F, {axisPt});
    registry.add("hTPCNsigma", "TPC electron Nsigma", kTH2F, {{axisPt}, {axisNsigma}});
    registry.add("hEnergy", "EMC cluster energy", kTH1F, {axisE});
    registry.add("hEnergyT", "EMC cluster transverse energy E/cosh(#eta)", kTH1F, {axisE});
    registry.add("hEnergyEt", "EMC cluster energy vs. transverse energy", kTH2F, {{axisE}, {axisE}});
    registry.add("hEnergyMult", "EMC cluster energy vs Multiplicity", kTH2F, {{axisCentrality}, {axisE}});
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
    registry.add("hInvMassZeeTrkQA", "Track QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisDCAxyneg, axisDCAxypos, axisDCAzpos, axisNclsTPCneg, axisNclsTPCpos, axisNclsITSneg, axisNclsITSpos, axisSectorTPCneg, axisSectorTPCneg, axisEopZneg, axisEopZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hInvMassZeeTrkQAbg", "Track QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisDCAxyneg, axisDCAxypos, axisDCAzpos, axisNclsTPCneg, axisNclsTPCpos, axisNclsITSneg, axisNclsITSpos, axisSectorTPCneg, axisSectorTPCneg, axisEopZneg, axisEopZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hInvMassZeeEmcQA", "EMCal QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisEopZneg, axisEopZpos, axisM02neg, axisM02pos, axisNsigmaZneg, axisNsigmaZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hInvMassZeeEmcQAbg", "EMCal QA for invariant mass for Z", HistType::kTHnSparseF, {axisInvMassZ, axisPtZneg, axisPtZpos, axisEopZneg, axisEopZpos, axisM02neg, axisM02pos, axisNsigmaZneg, axisNsigmaZpos, axisIsoEnergyZneg, axisIsoEnergyZpos, axisIsoMomentumZneg, axisIsoMomentumZpos, axisIsoTrackZneg, axisIsoTrackZpos});
    registry.add("hTHnElectrons", "electron info", HistType::kTHnSparseF, {axisPt, axisNsigma, axisM02, axisEop, axisIsoEnergy, axisIsoMomentum, axisIsoTrack, axisEta, axisEtCluster});
    registry.add("hTHnTrMatch", "Track EMC Match", HistType::kTHnSparseF, {axisPt, axisdPhi, axisdEta});
    registry.add("hNMatchPerTrack", "number of matched clusters within rMatchMax per track", kTH2F, {{axisPt}, {axisNMatch}});

    // Z-hadron correlation histograms
    const int nSysVar = nCutVariations();
    if (sysCfg.enableSystematics) {
      const size_t nVar = static_cast<size_t>(nSysVar - 1);
      if (nVar == 0) {
        LOGF(warn, "systematics: enableSystematics is on but every variation vector is empty, so the axis has one bin and nothing is varied.");
      }
      for (const auto& [values, name] : variationVectorsNamed()) {
        if (!values->empty() && values->size() != nVar) {
          LOGF(fatal, "systematics: %s has %zu entries but there are %zu variations. Give it either %zu entries or none at all.", name, values->size(), nVar, nVar);
        }
      }
      for (size_t i = 0; i < nVar; ++i) {
        const int iv = static_cast<int>(i) + 1;
        if (variedValue(sysCfg.varDcaxyMax.value, iv, dcaxyMax) > dcaxyMax ||
            variedValue(sysCfg.varNclCrossTpcMin.value, iv, nclcrossTpcMin) < nclcrossTpcMin ||
            variedValue(sysCfg.varNclItsMin.value, iv, nclItsMin) < nclItsMin) {
          LOGF(warn, "systematics: variation %zu loosens an associated-hadron cut. The container only holds tracks passing the nominal selection, so that variation is identical to the nominal there. Loosen the nominal configurable instead.", i + 1);
        }
      }
      LOGF(info, "systematics: %d bins on the variation axis (bin 0 = nominal)", nSysVar);
    }
    const AxisSpec axisSysVar{nSysVar, -0.5, nSysVar - 0.5, "cut variation"};

    registry.add("hZHadronDphi", "Z-hadron #Delta#phi correlation", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPtZ, axisDPhiZh, axisDEtaZh, axisZfrag, axisPtHadron, axisNsigmaHadron, axisSignHadron, axisSysVar});
    registry.add("hZptSpectrum", "Z boson p_{T} spectrum", kTH3F, {{axisSign}, {axisPtZ}, {axisSysVar}});

    if (weCfg.enableWeHadron) {
      registry.add("hWeHadronDphi", "isolated electron - hadron #Delta#phi correlation", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPtZ, axisDPhiZh, axisDEtaZh, axisIsoEnergy, axisIsoMomentum, axisPtHadron, axisNsigmaHadron, axisSignHadron, axisSysVar});
      registry.add("hWeSpectrum", "trigger electrons for the e-h correlation", HistType::kTHnSparseF, {axisCentrality, axisSign, axisPtZ, axisIsoEnergy, axisIsoMomentum, axisSysVar});
    }

    // hisotgram for EMCal trigger
    registry.add("hEMCalTrigger", "EMCal trigger", kTH1D, {axisTrigger});

    // histogram for UE
    registry.add("hRho", "rho UE density", kTH1F, {axisE});
    registry.add("hSumERC", "RC sumE", kTH1F, {axisE});
    registry.add("hEnergyUE", "UE vs. centrality", kTH2F, {{axisCentrality}, {axisE}});

    // MC : associated-hadron efficiency
    const AxisSpec axisPtHadronMc{mcCfg.confaxisPtHadron, "p_{T,hadron} (GeV/c)"};
    const AxisSpec axisEtaHadronMc{mcCfg.confaxisEtaHadron, "#eta_{hadron}"};
    const AxisSpec axisMcCollCounter{4, -0.5, 3.5, "MC collision selection step"};

    if (doprocessMcGen) {
      registry.add("hMcGenCollCounter", "generated collisions (0: all, 1: |zvtx| cut, 2: reconstructed & selected)", kTH1D, {axisMcCollCounter});
      registry.add("hMcGenZvtx", "generated Z vertex", kTH1D, {axisZvtx});
      // denominator of the efficiency
      registry.add("hMcGenHadron", "generated primary charged hadrons", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc});
    }

    if (doprocessMcRec) {
      registry.add("hMcRecCollCounter", "reconstructed collisions in MC", kTH1D, {axisCounter});
      // numerator of the efficiency, in generated kinematics (use with hMcGenHadron)
      registry.add("hMcRecHadronMatchedGenKine", "reconstructed primary charged hadrons, generated kinematics", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc, axisSysVar});
      // same tracks, but in reconstructed kinematics (use to correct in measured bins)
      registry.add("hMcRecHadronMatchedRecKine", "reconstructed primary charged hadrons, reconstructed kinematics", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc, axisSysVar});
      // purity : all selected tracks and the non-primary / unlabelled sub-samples
      registry.add("hMcRecHadronAll", "all selected associated tracks", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc, axisSysVar});
      registry.add("hMcRecHadronSecondary", "selected tracks matched to a non-primary particle", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc, axisSysVar});
      registry.add("hMcRecHadronFake", "selected tracks without MC particle", HistType::kTHnSparseF, {axisCentrality, axisPtHadronMc, axisEtaHadronMc, axisSysVar});
      // momentum resolution
      registry.add("hMcPtGenRec", "generated vs. reconstructed p_{T}", kTH2F, {{axisPtHadronMc}, {axisPtHadronMc}});
    }
  }

  // Kinematic + quality selection shared by the data and the MC track loops,
  // so that the efficiency is evaluated with exactly the same cuts.
  template <typename TTrack>
  bool isSelectedTrackKine(TTrack const& track)
  {
    if (!track.isGlobalTrackWoPtEta()) {
      return false;
    }
    if (std::abs(track.eta()) > etaTrMax) {
      return false;
    }
    if (applyTrackSys) {
      if (track.tpcNClsCrossedRows() < nclcrossTpcMin) {
        return false;
      }
      if (std::abs(track.dcaXY()) > dcaxyMax) {
        return false;
      }
      if (track.itsChi2NCl() > chi2ItsMax) {
        return false;
      }
      if (track.tpcChi2NCl() > chi2TpcMax) {
        return false;
      }
      if (track.tpcNClsFound() < nclTpcMin) {
        return false;
      }
      if (track.itsNCls() < nclItsMin) {
        return false;
      }
    }
    return true;
  }

  // Centrality selection, applied on the centrality estimator itself (as in processData,
  // where the cut is done before the variable is possibly overwritten by the multiplicity)
  template <typename TCollision>
  bool isSelectedCentrality(TCollision const& collision)
  {
    if (!enableCentralityAnalysis) {
      return true;
    }
    float const centrality = o2::hf_centrality::getCentralityColl(collision, centralityEstimator);
    return (centrality >= centralityMin && centrality <= centralityMax);
  }

  // Value filled on the centrality axis, following the same priority as in processData
  template <typename TCollision>
  float getCentralityValue(TCollision const& collision)
  {
    float centrality = 1.0;
    if (enableCentralityAnalysis) {
      centrality = o2::hf_centrality::getCentralityColl(collision, centralityEstimator);
    }
    if (enableMultiplicityFT0MAnalysis) {
      centrality = collision.multFT0M();
    }
    if (enableMultiplicityPVAnalysis) {
      centrality = collision.multNTracksPV();
    }
    return centrality;
  }

  // Generated particle counted as an associated hadron
  template <typename TParticle>
  bool isChargedPrimary(TParticle const& particle)
  {
    if (mcCfg.requirePhysicalPrimary && !particle.isPhysicalPrimary()) {
      return false;
    }
    const TParticlePDG* pdgParticle = pdg->GetParticle(particle.pdgCode());
    if (!pdgParticle) {
      return false;
    }
    if (std::abs(pdgParticle->Charge()) < ChargeTolerance) {
      return false;
    }
    return true;
  }

  double getIsolatedCluster(const o2::aod::EMCALCluster& cluster,
                            const SelectedClusters& clusters,
                            float energyUE)
  {
    double energySum = 0.0;
    double energySumExcl = 0.0;
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
    energySumExcl = energySum - cluster.energy();
    if (energySum > 0) {
      isoEnergy = (energySumExcl - energyUE) / cluster.energy();
    }

    // LOG(info) <<"clustE = " << cluster.energy() << " ; energySum = " << energySum << " ; nclust in Cone = " <<  nclustSum - 1  << " ; UE = " << energyUE << " ; isoEnergy = " << isoEnergy;
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
  float estimateRhoRC(const SelectedClusters& clusters)
  {
    const float randomConeR = rIsolation;
    const float randomConeArea = o2::constants::math::PI * randomConeR * randomConeR;

    std::vector<float> sumErc;
    sumErc.reserve(ueCfg.nRandomCones);

    for (int i = 0; i < ueCfg.nRandomCones; i++) {

      float etarc = rnd->Uniform(-etaEmcMax, etaEmcMax); // in EMCal acceptance
      float phirc = rnd->Uniform(phiEmcMin, phiEmcMax);  // in EMCal acceptance

      float energySumRC = 0;

      for (const auto& c : clusters) {
        if (c.energy() > ueCfg.rcHardE) {
          continue;
        }
        double dEtarc = etarc - c.eta();
        double dPhirc = phirc - c.phi();
        dPhirc = RecoDecay::constrainAngle(dPhirc, -o2::constants::math::PI);
        double const deltaRrc = std::sqrt(dEtarc * dEtarc + dPhirc * dPhirc);
        if (deltaRrc < randomConeR) {
          energySumRC += c.energy();
        }
      }

      registry.fill(HIST("hSumERC"), energySumRC);
      sumErc.push_back(energySumRC);
    }

    if (sumErc.empty()) {
      return 0;
    }
    std::nth_element(sumErc.begin(),
                     sumErc.begin() + sumErc.size() / 2,
                     sumErc.end());

    float median = sumErc[sumErc.size() / 2];
    // LOG(info) << "median = " << median;
    registry.fill(HIST("hRho"), median / randomConeArea);

    return median / randomConeArea;
  }

  // Outcome of the track - EMCal cluster matching for one track
  struct HfEmcMatchResult {
    int64_t posBestMatch = -1;
    double rMin = 999.9;
    double dPhiMin = 999.9;
    double dEtaMin = 999.9;
    int nMatchInRadius = 0;
    bool isMatched = false;
  };

  template <typename TMatchedTracks>
  HfEmcMatchResult findBestEmcMatch(TMatchedTracks const& matchesOfTrack, float centrality)
  {
    HfEmcMatchResult result;
    int64_t posMatch = -1;

    for (const auto& match : matchesOfTrack) {
      ++posMatch;
      const auto& cluster = match.template emcalcluster_as<SelectedClusters>();
      if (cluster.time() < timeEmcMin || cluster.time() > timeEmcMax) {
        continue;
      }

      double const dEta = match.template track_as<TrackEle>().trackEtaEmcal() - cluster.eta();
      double dPhi = match.template track_as<TrackEle>().trackPhiEmcal() - cluster.phi();
      dPhi = RecoDecay::constrainAngle(dPhi, -o2::constants::math::PI);
      double const r = RecoDecay::sqrtSumOfSquares(dPhi, dEta);

      registry.fill(HIST("hMatchPhi"), cluster.phi(), match.template track_as<TrackEle>().trackPhiEmcal());
      registry.fill(HIST("hMatchEta"), cluster.eta(), match.template track_as<TrackEle>().trackEtaEmcal());
      registry.fill(HIST("hTHnTrMatch"), match.template track_as<TrackEle>().pt(), dPhi, dEta);
      registry.fill(HIST("hEMCtime"), cluster.time());
      registry.fill(HIST("hEnergy"), cluster.energy());
      double const energyT = cluster.energy() / std::cosh(cluster.eta());
      registry.fill(HIST("hEnergyT"), energyT);
      registry.fill(HIST("hEnergyEt"), cluster.energy(), energyT);
      registry.fill(HIST("hEnergyMult"), centrality, cluster.energy());

      if (std::abs(dPhi) <= rMatchMax && std::abs(dEta) <= rMatchMax) {
        result.nMatchInRadius++;
      }

      // LOG(info) << "r match = " << r;
      if (r < result.rMin) {
        result.rMin = r;
        result.dPhiMin = dPhi;
        result.dEtaMin = dEta;
        result.posBestMatch = posMatch;
        // LOG(info) << "match cl info : " << r << " ; " << cluster.energy();
      }
    }

    result.isMatched = (result.posBestMatch >= 0) && (std::abs(result.dPhiMin) <= rMatchMax) && (std::abs(result.dEtaMin) <= rMatchMax);
    return result;
  }

  void recoMassZee(const KFParticle& kfpIsoEle,
                   int charge,
                   int64_t indexIsoEle,
                   float centrality,
                   TrackEle const& tracks)
  {
    // LOG(info) << "Invarimass cal by KF particle ";
    for (const auto& track : tracks) {

      if (track.globalIndex() == indexIsoEle) {
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
      std::array<const KFParticle*, 2> electronPairs = {&kfpIsoEle, &kfpAssEle};
      KFParticle zeeKF;
      zeeKF.SetConstructMethod(kfConstructMethod);
      zeeKF.Construct(electronPairs.data(), 2);
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
      float massZee{}, massZeeErr{};
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
        kfpAssEle.GetPt(),
        indexIsoEle,
        track.globalIndex());
    }
  }

  // -------------------------------------------------------------------------
  // Cut variations (systematics). Bin 0 is always the nominal configurables; bin iv>0 reads
  // entry iv-1 of each variation vector. Member FUNCTIONS, so none of this costs a task member.
  // -------------------------------------------------------------------------
  int nCutVariations()
  {
    if (!sysCfg.enableSystematics) {
      return 1;
    }
    size_t n = 0;
    for (const auto& [values, name] : variationVectorsNamed()) {
      n = std::max(n, values->size());
    }
    return 1 + static_cast<int>(n);
  }
  // The eleven variation vectors, paired with their names for the messages in init().
  std::array<std::pair<const std::vector<float>*, const char*>, 11> variationVectorsNamed()
  {
    return {{{&sysCfg.varDcaxyMax.value, "varDcaxyMax"},
             {&sysCfg.varNclCrossTpcMin.value, "varNclCrossTpcMin"},
             {&sysCfg.varNclItsMin.value, "varNclItsMin"},
             {&sysCfg.varEopMin.value, "varEopMin"},
             {&sysCfg.varEopMax.value, "varEopMax"},
             {&sysCfg.varNsigTpcMin.value, "varNsigTpcMin"},
             {&sysCfg.varNsigTpcMax.value, "varNsigTpcMax"},
             {&sysCfg.varM02Min.value, "varM02Min"},
             {&sysCfg.varM02Max.value, "varM02Max"},
             {&sysCfg.varMassZMin.value, "varMassZMin"},
             {&sysCfg.varMassZMax.value, "varMassZMax"}}};
  }
  // One variation's value of a cut, falling back to the nominal for bin 0 or a short vector.
  static float variedValue(const std::vector<float>& values, int iv, float nominal)
  {
    if (iv <= 0 || static_cast<size_t>(iv - 1) >= values.size()) {
      return nominal;
    }
    return values[iv - 1];
  }

  // Associated hadron. dcaxyTrk, nclusterTPC (= crossed rows, see where the container is
  // filled) and nclusterITS are the three quantities the candidate carries.
  bool passHadronVariation(const HfElectronCandidate& hadron, int iv)
  {
    if (iv == 0) {
      return true; // the container already holds the nominal selection
    }
    if (std::abs(hadron.dcaxyTrk) > variedValue(sysCfg.varDcaxyMax.value, iv, dcaxyMax)) {
      return false;
    }
    if (hadron.nclusterTPC < variedValue(sysCfg.varNclCrossTpcMin.value, iv, nclcrossTpcMin)) {
      return false;
    }
    if (hadron.nclusterITS < variedValue(sysCfg.varNclItsMin.value, iv, nclItsMin)) {
      return false;
    }
    return true;
  }

  // Same three cuts on a track, for the MC efficiency: the numerator has to be recomputed with
  // the very selection the correlation used, or the correction does not match the data.
  template <typename TTrack>
  bool passHadronVariationTrack(TTrack const& track, int iv)
  {
    if (iv == 0) {
      return true;
    }
    if (std::abs(track.dcaXY()) > variedValue(sysCfg.varDcaxyMax.value, iv, dcaxyMax)) {
      return false;
    }
    if (track.tpcNClsCrossedRows() < variedValue(sysCfg.varNclCrossTpcMin.value, iv, nclcrossTpcMin)) {
      return false;
    }
    if (track.itsNCls() < variedValue(sysCfg.varNclItsMin.value, iv, nclItsMin)) {
      return false;
    }
    return true;
  }

  // The Z mass window for one variation. Unlike the associated-hadron cuts this one is
  // two-sided: reconstructedZ holds every pair that passed the KF chi2 cut, at any mass, so a
  // wider window genuinely adds candidates. iv == 0 uses the nominal configurables.
  bool passZVariation(const HfZeeCandidate& zBoson, int iv)
  {
    if (zBoson.mass < variedValue(sysCfg.varMassZMin.value, iv, massZMin) ||
        zBoson.mass > variedValue(sysCfg.varMassZMax.value, iv, massZMax)) {
      return false;
    }
    return true;
  }

  // W trigger electron. selectedElectronsIso carries no identification cut, so these three can
  // be varied in both directions.
  bool passElectronVariation(const HfElectronCandidate& ele, int iv)
  {
    if (ele.eop < variedValue(sysCfg.varEopMin.value, iv, eopMin) ||
        ele.eop > variedValue(sysCfg.varEopMax.value, iv, eopMax)) {
      return false;
    }
    if (ele.nSigmaTpc < variedValue(sysCfg.varNsigTpcMin.value, iv, nsigTpcMin) ||
        ele.nSigmaTpc > variedValue(sysCfg.varNsigTpcMax.value, iv, nsigTpcMax)) {
      return false;
    }
    if (ele.m02Emc < variedValue(sysCfg.varM02Min.value, iv, m02Min) ||
        ele.m02Emc > variedValue(sysCfg.varM02Max.value, iv, m02Max)) {
      return false;
    }
    return true;
  }

  // -------------------------------------------------------------------------
  // Z->ee QA.
  // -------------------------------------------------------------------------
  template <bool IsSignal, typename TCandidate>
  void fillZeeQaPair(const TCandidate& trackNeg, const TCandidate& trackPos)
  {
    auto child1 = RecoDecayPtEtaPhi::pVector(trackNeg.pt, trackNeg.eta, trackNeg.phi);
    auto child2 = RecoDecayPtEtaPhi::pVector(trackPos.pt, trackPos.eta, trackPos.phi);
    double const invMass = RecoDecay::m(std::array{child1, child2}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
    if (invMass <= zeeQaCfg.massZMinQA) {
      return;
    }
    float const sectorneg = trackNeg.phi / o2::constants::math::SectorSpanRad;
    float const sectorpos = trackPos.phi / o2::constants::math::SectorSpanRad;

    if (zeeQaCfg.enableZeeTrkRecoQA) {
      if constexpr (IsSignal) {
        registry.fill(HIST("hInvMassZeeTrkQA"), invMass, trackNeg.pt, trackPos.pt, trackNeg.dcaxyTrk, trackPos.dcaxyTrk, trackPos.dcazTrk, trackNeg.nclusterTPC, trackPos.nclusterTPC, trackNeg.nclusterITS, trackPos.nclusterITS, sectorneg, sectorpos, trackNeg.eop, trackPos.eop, trackNeg.energyIso, trackPos.energyIso, trackNeg.momIso, trackPos.momIso, trackNeg.ntrackIso, trackPos.ntrackIso);
      } else {
        registry.fill(HIST("hInvMassZeeTrkQAbg"), invMass, trackNeg.pt, trackPos.pt, trackNeg.dcaxyTrk, trackPos.dcaxyTrk, trackPos.dcazTrk, trackNeg.nclusterTPC, trackPos.nclusterTPC, trackNeg.nclusterITS, trackPos.nclusterITS, sectorneg, sectorpos, trackNeg.eop, trackPos.eop, trackNeg.energyIso, trackPos.energyIso, trackNeg.momIso, trackPos.momIso, trackNeg.ntrackIso, trackPos.ntrackIso);
      }
    }
    if (zeeQaCfg.enableZeeEmcRecoQA) {
      if constexpr (IsSignal) {
        registry.fill(HIST("hInvMassZeeEmcQA"), invMass, trackNeg.pt, trackPos.pt, trackNeg.eop, trackPos.eop, trackNeg.m02Emc, trackPos.m02Emc, trackNeg.nSigmaTpc, trackPos.nSigmaTpc, trackNeg.energyIso, trackPos.energyIso, trackNeg.momIso, trackPos.momIso, trackNeg.ntrackIso, trackPos.ntrackIso);
      } else {
        registry.fill(HIST("hInvMassZeeEmcQAbg"), invMass, trackNeg.pt, trackPos.pt, trackNeg.eop, trackPos.eop, trackNeg.m02Emc, trackPos.m02Emc, trackNeg.nSigmaTpc, trackPos.nSigmaTpc, trackNeg.energyIso, trackPos.energyIso, trackNeg.momIso, trackPos.momIso, trackNeg.ntrackIso, trackPos.ntrackIso);
      }
    }
  }

  void fillZeeQa()
  {
    if (!zeeQaCfg.enableZeeTrkRecoQA && !zeeQaCfg.enableZeeEmcRecoQA) {
      return;
    }

    for (const auto& trackEle : selectedElectronsIso) {
      for (const auto& trackPos : selectedPositronsIso) {
        fillZeeQaPair<true>(trackEle, trackPos);
      }
    }
    // each unordered like-sign pair once
    for (size_t i = 0; i < selectedElectronsIso.size(); ++i) {
      for (size_t j = i + 1; j < selectedElectronsIso.size(); ++j) {
        fillZeeQaPair<false>(selectedElectronsIso[i], selectedElectronsIso[j]);
      }
    }
    for (size_t i = 0; i < selectedPositronsIso.size(); ++i) {
      for (size_t j = i + 1; j < selectedPositronsIso.size(); ++j) {
        fillZeeQaPair<false>(selectedPositronsIso[i], selectedPositronsIso[j]);
      }
    }
  }

  // -------------------------------------------------------------------------
  // Z - h correlation.
  // -------------------------------------------------------------------------
  void fillZHadron(float centrality)
  {
    const int nVar = nCutVariations();
    for (const auto& zBoson : reconstructedZ) {
      // The Z counter is filled per variation, so that a varied mass window moves the
      // denominator of the per-trigger yield together with its numerator.
      bool anyVariation = false;
      for (int iv = 0; iv < nVar; ++iv) {
        if (!passZVariation(zBoson, iv)) {
          continue;
        }
        anyVariation = true;
        registry.fill(HIST("hZptSpectrum"), zBoson.charge, zBoson.pt, iv);
      }
      if (!anyVariation) { // outside the widest window of any variation
        continue;
      }

      for (const auto& trackAss : selectedElectronsAss) {
        // remove the two tracks used to build this Z, matched exactly by globalIndex
        if (trackAss.indexTrk == zBoson.indexChild0 || trackAss.indexTrk == zBoson.indexChild1) {
          continue;
        }
        double const deltaPhi = RecoDecay::constrainAngle(trackAss.phi - zBoson.phi, -o2::constants::math::PIHalf);
        double const deltaEta = zBoson.eta - trackAss.eta;
        double const ptRatio = trackAss.pt / zBoson.pt;
        // one entry per cut variation the pair passes; bin 0 is the nominal selection
        for (int iv = 0; iv < nVar; ++iv) {
          if (!passZVariation(zBoson, iv) || !passHadronVariation(trackAss, iv)) {
            continue;
          }
          registry.fill(HIST("hZHadronDphi"), centrality, zBoson.charge, zBoson.pt, deltaPhi, deltaEta,
                        ptRatio, trackAss.pt, trackAss.nSigmaTpc, trackAss.sign, iv);
        }
      }
    }
  }

  // -------------------------------------------------------------------------
  // W->e - h correlation.
  // -------------------------------------------------------------------------
  bool isZeeDaughter(int64_t indexTrk)
  {
    for (const auto& zBoson : reconstructedZ) {
      if (zBoson.charge >= 0) { // like-sign pairs are the combinatorial background, not a Z
        continue;
      }
      if (zBoson.mass < massZMin || zBoson.mass > massZMax) {
        continue;
      }
      if (indexTrk == zBoson.indexChild0 || indexTrk == zBoson.indexChild1) {
        return true;
      }
    }
    return false;
  }

  void fillWeHadron(const std::vector<HfElectronCandidate>& electrons, float centrality)
  {
    for (const auto& ele : electrons) {
      if (ele.pt < weCfg.ptWeMin) {
        continue;
      }
      if (ele.eop <= 0.f) {
        continue;
      }
      if (weCfg.vetoZeeElectron && isZeeDaughter(ele.indexTrk)) {
        continue;
      }

      const int nVar = nCutVariations();
      for (int iv = 0; iv < nVar; ++iv) {
        if (!passElectronVariation(ele, iv)) {
          continue;
        }
        registry.fill(HIST("hWeSpectrum"), centrality, ele.sign, ele.pt, ele.energyIso, ele.momIso, iv);

        for (const auto& trackAss : selectedElectronsAss) {
          // the trigger electron itself, matched exactly by globalIndex
          if (trackAss.indexTrk == ele.indexTrk) {
            continue;
          }
          if (!passHadronVariation(trackAss, iv)) {
            continue;
          }
          double const deltaPhi = RecoDecay::constrainAngle(trackAss.phi - ele.phi, -o2::constants::math::PIHalf);
          double const deltaEta = ele.eta - trackAss.eta;
          registry.fill(HIST("hWeHadronDphi"), centrality, ele.sign, ele.pt, deltaPhi, deltaEta,
                        ele.energyIso, ele.momIso, trackAss.pt, trackAss.nSigmaTpc, trackAss.sign, iv);
        }
      }
    }
  }

  void processData(soa::Filtered<CollisionsWithCent>::iterator const& collision,
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
      // LOG(info) << "centrality = " << o2::hf_centrality::getCentralityColl(collision, centralityEstimator) << " ; FTC = " << collision.multFT0C();
      if (centrality < centralityMin || centrality > centralityMax) {
        return;
      }
      registry.fill(HIST("hCentMultFT0Corr"), centrality, collision.multFT0M());
      registry.fill(HIST("hCentMultPVCorr"), centrality, collision.multNTracksPV());
    }

    if (enableMultiplicityFT0MAnalysis || enableMultiplicityPVAnalysis) {
      if (enableMultiplicityFT0MAnalysis) {
        centrality = collision.multFT0M();
      }
      if (enableMultiplicityPVAnalysis) {
        centrality = collision.multNTracksPV();
      }
      // LOG(info) << "raw mult PV = " << collision.multNTracksPV();
      // LOG(info) << "raw mult FT0M = " << collision.multFT0M();
      registry.fill(HIST("hMultPV"), collision.posZ(), collision.multNTracksPV());
      registry.fill(HIST("hMultFT0"), collision.posZ(), collision.multFT0M());
      registry.fill(HIST("hMultFT0PV"), collision.multFT0M(), collision.multNTracksPV());
    }

    registry.fill(HIST("hCentrality"), centrality);

    // UE estimate
    float rho = 0.f;
    float energyUE = 0.f;

    if (ueCfg.useUEsub) {
      rho = estimateRhoRC(emcClusters);
      energyUE = rho * (o2::constants::math::PI * rIsolation * rIsolation);
      registry.fill(HIST("hEnergyUE"), centrality, energyUE);
      // LOG(info) << "UE = " << energyUE;
    }

    const float ptIsoStoreMin = weCfg.enableWeHadron
                                  ? std::min(static_cast<float>(ptZeeMin), static_cast<float>(weCfg.ptWeMin))
                                  : static_cast<float>(ptZeeMin);

    // track loop
    for (const auto& track : tracks) {

      if (!isSelectedTrackKine(track)) {
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
      float m02 = -0.01;
      float isoEnergy = 99.0;
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
          track.tpcNSigmaPi(),
          m02,
          trackCount,
          track.tpcNClsCrossedRows(),
          track.itsNCls(),
          track.sign(),
          track.globalIndex());
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

        // Pass 1 : QA over every matched cluster, and find the closest one
        auto const emcMatch = findBestEmcMatch(tracksofcluster, centrality);
        rMin = emcMatch.rMin;
        dPhiMin = emcMatch.dPhiMin;
        dEtaMin = emcMatch.dEtaMin;
        registry.fill(HIST("hNMatchPerTrack"), track.pt(), emcMatch.nMatchInRadius);

        // Pass 2 : only the closest cluster is used for E/p, isolation and the Z reconstruction
        if (emcMatch.isMatched) {
          auto const match = tracksofcluster.iteratorAt(emcMatch.posBestMatch);
          const auto& cluster = match.emcalcluster_as<SelectedClusters>();
          float const m02Emc = cluster.m02();
          float const energyEmc = cluster.energy();
          // transverse energy of the cluster; cosh(eta) >= 1, so energyTEmc <= energyEmc
          float const energyTEmc = energyEmc / std::cosh(cluster.eta());

          registry.fill(HIST("hTrMatchR"), match.track_as<TrackEle>().pt(), rMin);
          registry.fill(HIST("hEnergyNcell"), energyEmc, cluster.nCells());

          eop = energyEmc / match.track_as<TrackEle>().p();
          m02 = m02Emc;
          // LOG(info) << "check match cl info : " << cluster.energy();
          isoEnergy = getIsolatedCluster(cluster, emcClusters, energyUE);

          if (match.track_as<TrackEle>().pt() > ptTHnThresh && isTHnElectron) {
            registry.fill(HIST("hTHnElectrons"), match.track_as<TrackEle>().pt(), match.track_as<TrackEle>().tpcNSigmaEl(), m02Emc, eop, isoEnergy, isoMomentum, trackCount, track.eta(), energyTEmc);
          }
          // LOG(info) << "E/p" << eop;
          registry.fill(HIST("hEopNsigTPC"), match.track_as<TrackEle>().tpcNSigmaEl(), eop);

          if (m02Emc >= m02Min && m02Emc <= m02Max &&
              match.track_as<TrackEle>().tpcNSigmaEl() > nsigTpcMin && match.track_as<TrackEle>().tpcNSigmaEl() < nsigTpcMax) {
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
                recoMassZee(kfpIsoEle, match.track_as<TrackEle>().sign(), match.track_as<TrackEle>().globalIndex(), centrality, tracks);

              } // end of pt cut for e from Z
            } // end if isolation cut
            if (isIsolatedTr) {
              registry.fill(HIST("hEopIsolationTr"), match.track_as<TrackEle>().pt(), eop);
            }
          } // end of M02 + PID cut
        } // end of isMatched
      } // end of cluster

      if (rMin < rMatchMax) {
        // LOG(info) << "R mim = " << rMin;
        registry.fill(HIST("hTrMatch_mim"), dPhiMin, dEtaMin);
      }
      // The e-h correlation reads the same two containers, so they must also be filled when
      // only that analysis is enabled, and down to the looser of the two pT thresholds.
      if ((zeeQaCfg.enableZeeTrkRecoQA || zeeQaCfg.enableZeeEmcRecoQA || weCfg.enableWeHadron) && track.pt() > ptIsoStoreMin) {
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
            track.tpcNSigmaEl(),
            m02,
            trackCount,
            track.tpcNClsFound(),
            track.itsNCls(),
            track.sign(),
            track.globalIndex());
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
            track.tpcNSigmaEl(),
            m02,
            trackCount,
            track.tpcNClsFound(),
            track.itsNCls(),
            track.sign(),
            track.globalIndex());
        }
      }

    } // end of track loop
    // Z - hadron correlation
    fillZHadron(centrality);

    // W->e - h correlation. Electrons and positrons are kept in two containers, so the same
    // routine runs twice; the electron charge goes on the sign axis (W- / W+).
    if (weCfg.enableWeHadron) {
      fillWeHadron(selectedElectronsIso, centrality);
      fillWeHadron(selectedPositronsIso, centrality);
    } // end of W->e - hadron correlation
    // Z->ee QA
    fillZeeQa();
  } // processData
  PROCESS_SWITCH(HfTaskElectronWeakBoson, processData, "process data (Z->ee and Z-h correlation)", true);

  // -------------------------------------------------------------------------
  // MC : denominator of the associated-hadron efficiency.
  // Looping over the generated collisions (with the reconstructed ones collected by
  // soa::SmallGroups) also gives access to the event-loss / split-collision effects
  // through hMcGenCollCounter.
  // -------------------------------------------------------------------------
  void processMcGen(aod::McCollision const& mcCollision,
                    soa::SmallGroups<CollisionsWithCentMc> const& collisions,
                    aod::McParticles const& mcParticles)
  {
    registry.fill(HIST("hMcGenCollCounter"), 0.5);

    if (std::abs(mcCollision.posZ()) > vtxZ) {
      return;
    }
    registry.fill(HIST("hMcGenCollCounter"), 1.5);
    registry.fill(HIST("hMcGenZvtx"), mcCollision.posZ());

    // look for a reconstructed collision passing the same event selection as in processData
    bool isReconstructed = false;
    float centrality = 1.0;
    for (const auto& collision : collisions) {
      if (applySel8 && !collision.sel8()) {
        continue;
      }
      if (std::abs(collision.posZ()) > vtxZ) {
        continue;
      }
      if (!isSelectedCentrality(collision)) {
        continue;
      }
      isReconstructed = true;
      centrality = getCentralityValue(collision);
      break;
    }
    if (!isReconstructed) {
      return;
    }
    registry.fill(HIST("hMcGenCollCounter"), 2.5);

    for (const auto& particle : mcParticles) {
      if (!isChargedPrimary(particle)) {
        continue;
      }
      // same acceptance as isSelectedTrackKine() applied to the reconstructed tracks
      if (std::abs(particle.eta()) > etaTrMax) {
        continue;
      }
      if (particle.pt() < ptAssMin) {
        continue;
      }
      registry.fill(HIST("hMcGenHadron"), centrality, particle.pt(), particle.eta());
    }
  } // processMcGen
  PROCESS_SWITCH(HfTaskElectronWeakBoson, processMcGen, "process MC generated level for the associated-hadron efficiency", false);

  // -------------------------------------------------------------------------
  // MC : numerator of the associated-hadron efficiency (and purity).
  // The track selection is identical to the one applied to selectedElectronsAss
  // in processData.
  // -------------------------------------------------------------------------
  void processMcRec(soa::Filtered<CollisionsWithCentMc>::iterator const& collision,
                    TrackEleMc const& tracks,
                    aod::McParticles const&,
                    aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }

    if (!isSelectedCentrality(collision)) {
      return;
    }
    float const centrality = getCentralityValue(collision);
    registry.fill(HIST("hMcRecCollCounter"), 0.5);

    for (const auto& track : tracks) {
      if (!isSelectedTrackKine(track)) {
        continue;
      }
      if (track.pt() < ptAssMin) {
        continue;
      }

      // The efficiency has to be recomputed with the same cut variation the correlation used,
      // so the numerator histograms carry the variation axis too. hMcGenHadron does not: no
      // track cut enters the denominator, so it is common to every variation.
      const int nVar = nCutVariations();

      for (int iv = 0; iv < nVar; ++iv) {
        if (passHadronVariationTrack(track, iv)) {
          registry.fill(HIST("hMcRecHadronAll"), centrality, track.pt(), track.eta(), iv);
        }
      }

      if (!track.has_mcParticle()) {
        for (int iv = 0; iv < nVar; ++iv) {
          if (passHadronVariationTrack(track, iv)) {
            registry.fill(HIST("hMcRecHadronFake"), centrality, track.pt(), track.eta(), iv);
          }
        }
        continue;
      }

      auto const mcParticle = track.mcParticle();

      // reject pileup : the matched particle must come from this MC collision
      const bool isSecondary = (mcCfg.requireSameMcCollision && mcParticle.mcCollisionId() != collision.mcCollisionId()) || !isChargedPrimary(mcParticle);
      if (isSecondary) {
        for (int iv = 0; iv < nVar; ++iv) {
          if (passHadronVariationTrack(track, iv)) {
            registry.fill(HIST("hMcRecHadronSecondary"), centrality, track.pt(), track.eta(), iv);
          }
        }
        continue;
      }

      for (int iv = 0; iv < nVar; ++iv) {
        if (!passHadronVariationTrack(track, iv)) {
          continue;
        }
        registry.fill(HIST("hMcRecHadronMatchedGenKine"), centrality, mcParticle.pt(), mcParticle.eta(), iv);
        registry.fill(HIST("hMcRecHadronMatchedRecKine"), centrality, track.pt(), track.eta(), iv);
      }
      registry.fill(HIST("hMcPtGenRec"), mcParticle.pt(), track.pt());
    }
  } // processMcRec
  PROCESS_SWITCH(HfTaskElectronWeakBoson, processMcRec, "process MC reconstructed level for the associated-hadron efficiency", false);

}; // struct HfTaskElectronWeakBoson

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskElectronWeakBoson>(cfgc)};
}
