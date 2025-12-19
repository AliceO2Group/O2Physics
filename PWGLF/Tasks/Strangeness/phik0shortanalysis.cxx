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
/// \file phik0shortanalysis.cxx
/// \brief Analysis task for the Phi and K0S rapidity correlations analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/TableHelper.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

enum {
  kGlobalplusITSonly = 0,
  kGlobalonly,
  kITSonly
};

enum {
  kSpAll = 0,
  kSpPion,
  kSpKaon,
  kSpProton,
  kSpOther,
  kSpStrangeDecay,
  kSpNotPrimary
};

enum {
  kNoGenpTVar = 0,
  kGenpTup,
  kGenpTdown
};

struct Phik0shortanalysis {
  // Histograms are defined with HistogramRegistry
  HistogramRegistry dataEventHist{"dataEventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcEventHist{"mcEventHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry dataPhiHist{"dataPhiHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcPhiHist{"mcPhiHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhiHist{"closureMCPhiHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry dataK0SHist{"dataK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcK0SHist{"mcK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry dataPhiK0SHist{"dataPhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcPhiK0SHist{"mcPhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhiK0SHist{"closureMCPhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry dataPionHist{"dataPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcPionHist{"mcPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry dataPhiPionHist{"dataPhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mcPhiPionHist{"mcPhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry closureMCPhiPionHist{"closureMCPhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mePhiK0SHist{"mePhiK0SHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry mePhiPionHist{"mePhiPionHist", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  struct : ConfigurableGroup {
    Configurable<bool> isData{"isData", true, "Data"};
    Configurable<bool> isMC{"isMC", true, "MC"};
    Configurable<bool> isClosure{"isClosure", true, "MC Closure"};

    Configurable<bool> isDataNewProc{"isDataNewProc", true, "New procedure  for Data"};
    Configurable<bool> isMCNewProc{"isMCNewProc", true, "New procedure for MC"};
    Configurable<bool> isClosureNewProc{"isClosureNewProc", true, "New procedure for MC Closure"};
    Configurable<bool> isMENewProc{"isMENewProc", true, "New procedure for ME"};
  } analysisModeConfigs;

  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for track selection (not necessarily common for trigger and the two associated particles)
  struct : ConfigurableGroup {
    Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on charge"};
    Configurable<float> cfgMinAbsCharge{"cfgMinAbsCharge", 3.0f, "Cut on absolute charge"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinChargedParticlePtcut{"cMinChargedParticlePtcut", 0.1f, "Track minimum pt cut"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
    Configurable<float> pTToUseTOF{"pTToUseTOF", 0.5f, "pT above which use TOF"};
    Configurable<float> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0f, "Track DCAz cut to PV Maximum"};
    Configurable<float> cMaxDCArToPV1Phi{"cMaxDCArToPV1Phi", 0.004f, "Track DCAr cut to PV config 1 for Phi"};
    Configurable<float> cMaxDCArToPV2Phi{"cMaxDCArToPV2Phi", 0.013f, "Track DCAr cut to PV config 2 for Phi"};
    Configurable<float> cMaxDCArToPV3Phi{"cMaxDCArToPV3Phi", 1.0f, "Track DCAr cut to PV config 3 for Phi"};
    Configurable<float> cMaxDCArToPV1Pion{"cMaxDCArToPV1Pion", 0.004f, "Track DCAr cut to PV config 1 for Pions"};
    Configurable<float> cMaxDCArToPV2Pion{"cMaxDCArToPV2Pion", 0.013f, "Track DCAr cut to PV config 2 for Pions"};
    Configurable<float> cMaxDCArToPV3Pion{"cMaxDCArToPV3Pion", 1.0f, "Track DCAr cut to PV config 3 for Pions"};

    Configurable<bool> cfgIsDCAzParameterized{"cfgIsDCAzParameterized", false, "IsDCAzParameterized"};
    Configurable<float> cMaxDCAzToPV1Pion{"cMaxDCAzToPV1Pion", 0.004f, "Track DCAz cut to PV config 1 for Pion"};
    Configurable<float> cMaxDCAzToPV2Pion{"cMaxDCAzToPV2Pion", 0.013f, "Track DCAz cut to PV config 2 for Pion"};
    Configurable<float> cMaxDCAzToPV3Pion{"cMaxDCAzToPV3Pion", 1.0f, "Track DCAz cut to PV config 3 for Pion"};

    Configurable<bool> isNoTOF{"isNoTOF", false, "isNoTOF"};
    Configurable<float> nSigmaCutTPCKa{"nSigmaCutTPCKa", 3.0f, "Value of the TPC Nsigma cut for Kaons"};
    Configurable<float> nSigmaCutCombinedKa{"nSigmaCutCombinedKa", 3.0f, "Value of the TOF Nsigma cut for Kaons"};

    Configurable<float> nSigmaCutTPCPion{"nSigmaCutTPCPion", 4.0f, "Value of the TPC Nsigma cut for Pions"};
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.2f, "Track minimum pt cut"};
    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
    Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};
    Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};

    Configurable<bool> applyExtraPhiCuts{"applyExtraPhiCuts", false, "Enable extra phi cut"};
    Configurable<std::vector<float>> extraPhiCuts{"extraPhiCuts", {3.07666f, 3.12661f, 0.03f, 6.253f}, "Extra phi cuts"};
  } trackConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};

  // Configurables on phi selection
  struct : ConfigurableGroup {
    Configurable<int> nBinsMPhi{"nBinsMPhi", 13, "N bins in cfgmassPhiaxis"};
    Configurable<float> lowMPhi{"lowMPhi", 1.0095f, "Upper limits on Phi mass for signal extraction"};
    Configurable<float> upMPhi{"upMPhi", 1.029f, "Upper limits on Phi mass for signal extraction"};

    Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi"};
    Configurable<float> maxPhiPt{"maxPhiPt", 10.0f, "Maximum pT for Phi"};
  } phiConfigs;

  // Configurables for V0 selection
  struct : ConfigurableGroup {
    Configurable<float> v0SettingCosPA{"v0SettingCosPA", 0.98f, "V0 CosPA"};
    Configurable<float> v0SettingRadius{"v0SettingRadius", 0.5f, "v0radius"};
    Configurable<float> v0SettingDCAV0Dau{"v0SettingDCAV0Dau", 1.0f, "DCA V0 Daughters"};
    Configurable<float> v0SettingDCAPosToPV{"v0SettingDCAPosToPV", 0.1f, "DCA Pos To PV"};
    Configurable<float> v0SettingDCANegToPV{"v0SettingDCANegToPV", 0.1f, "DCA Neg To PV"};
    Configurable<float> v0SettingMinPt{"v0SettingMinPt", 0.1f, "V0 min pt"};

    Configurable<bool> cfgisV0ForData{"cfgisV0ForData", true, "isV0ForData"};

    Configurable<bool> cfgFurtherV0Selection{"cfgFurtherV0Selection", false, "Further V0 selection"};
    Configurable<float> ctauK0s{"ctauK0s", 20.0f, "C tau K0s(cm)"};
    Configurable<float> paramArmenterosCut{"paramArmenterosCut", 0.2f, "parameter Armenteros Cut"};
    Configurable<float> v0rejK0s{"v0rejK0s", 0.005f, "V0 rej K0s"};

    Configurable<float> lowMK0S{"lowMK0S", 0.48f, "Lower limit on K0Short mass"};
    Configurable<float> upMK0S{"upMK0S", 0.52f, "Upper limit on K0Short mass"};
  } v0Configs;

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  // Configurables for delta y selection
  struct : ConfigurableGroup {
    Configurable<int> nBinsY{"nBinsY", 20, "Number of bins in y axis"};
    Configurable<int> nBinsDeltaY{"nBinsDeltaY", 20, "Number of bins in deltay axis"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
    Configurable<float> cfgFCutOnDeltaY{"cfgFCutOnDeltaY", 0.5f, "First upper bound on Deltay selection"};
    Configurable<float> cfgSCutOnDeltaY{"cfgSCutOnDeltaY", 0.1f, "Second upper bound on Deltay selection"};
    Configurable<std::vector<float>> cfgDeltaYAcceptanceBins{"cfgDeltaYAcceptanceBins", {0.5f}, "Rapidity acceptance bins"};
  } deltaYConfigs;

  // Configurable for RecMC
  Configurable<bool> cfgiskNoITSROFrameBorder{"cfgiskNoITSROFrameBorder", false, "kNoITSROFrameBorder request on RecMC collisions"};

  // Configurables for MC closure
  Configurable<bool> cfgisRecMCWPDGForClosure1{"cfgisRecMCWPDGForClosure1", false, "RecoMC with PDG Codes for Closure only for Associated particles"};
  Configurable<bool> cfgisRecMCWPDGForClosure2{"cfgisRecMCWPDGForClosure2", false, "RecoMC with PDG Codes for Closure"};
  Configurable<bool> cfgisGenMCForClosure{"cfgisGenMCForClosure", false, "GenMC for Closure"};

  // Configurables to choose the filling method
  Configurable<bool> fillMethodMultipleWeights{"fillMethodMultipleWeights", true, "Fill method Multiple Weights"};
  Configurable<bool> fillMethodSingleWeight{"fillMethodSingleWeight", false, "Fill method Single Weight"};
  Configurable<bool> applyEfficiency{"applyEfficiency", false, "Use efficiency for filling histograms"};

  // Configurables for dN/deta with phi computation
  Configurable<int> filterOnGenPhi{"filterOnGenPhi", 1, "Filter on Gen Phi (0: K+K- pair like Phi, 1: proper Phi)"};
  Configurable<int> filterOnRecoPhi{"filterOnRecoPhi", 1, "Filter on Reco Phi (0: without PDG, 1: with PDG)"};
  Configurable<bool> fillMcPartsForAllReco{"fillMcPartsForAllReco", false, "Fill MC particles for all associated reco collisions"};

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurable for CCDB
  Configurable<bool> useCCDB{"useCCDB", false, "Use CCDB for corrections"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbPurityPath{"ccdbPurityPath", "Users/s/scannito/PhiPuritiesData", "Correction path to file"};
  Configurable<std::string> ccdbEfficiencyPath{"ccdbEfficiencyPath", "Users/s/scannito/Efficiencies", "Correction path to file"};

  // Constants
  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;
  double massK0S = o2::constants::physics::MassK0Short;
  double massLambda = o2::constants::physics::MassLambda0;

  // Defining track flags
  static constexpr TrackSelectionFlags::flagtype TrackSelectionITS = TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF | TrackSelectionFlags::kITSHits;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionTPC = TrackSelectionFlags::kTPCNCls | TrackSelectionFlags::kTPCCrossedRowsOverNCls | TrackSelectionFlags::kTPCChi2NDF;
  static constexpr TrackSelectionFlags::flagtype TrackSelectionDCA = TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZVertex);

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0Configs.v0SettingDCAPosToPV && nabs(aod::v0data::dcanegtopv) > v0Configs.v0SettingDCANegToPV && aod::v0data::dcaV0daughters < v0Configs.v0SettingDCAV0Dau);

  // Defining filters on tracks
  Filter trackFilter = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionITS) &&
                       ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC), ncheckbit(aod::track::trackCutFlag, TrackSelectionTPC), true) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionDCA) &&
                       ncheckbit(aod::track::trackCutFlag, TrackSelectionFlags::kInAcceptanceTracks);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the V0s
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Join<FullV0s, aod::McV0Labels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;
  using FilteredTracks = soa::Filtered<FullTracks>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;
  using FilteredMCTracks = soa::Filtered<FullMCTracks>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  // Defining binning policy and axis for mixed event
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {20, -10, 10}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {20, 0, 100}, "Multiplicity percentil binning for mixing"};

  using BinningTypeVertexCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeVertexCent binningOnVertexAndCent{{axisVertexMixing, axisCentralityMixing}, true};

  // Cache for manual slicing
  SliceCache cache;

  // Preslice for manual slicing
  struct : PresliceGroup {
    Preslice<aod::Tracks> trackPerCollision = aod::track::collisionId;
    Preslice<aod::McParticles> mcPartPerMCCollision = aod::mcparticle::mcCollisionId;
    Preslice<aod::V0Datas> v0PerCollision = aod::v0::collisionId;
  } preslices;

  // Positive and negative tracks partitions
  Partition<FullTracks> posTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FilteredTracks> posFiltTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FilteredTracks> negFiltTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FilteredMCTracks> posFiltMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FilteredMCTracks> negFiltMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  // Necessary to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Necessary to get the CCDB for phi purities
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Set of functions for phi purity
  std::vector<std::vector<TF1*>> phiPurityFunctions = std::vector<std::vector<TF1*>>(binsMult->size(), std::vector<TF1*>(binspTPhi->size(), nullptr));

  // Efficiecy maps
  TH3F* effMapPhi;
  TH3F* effMapK0S;
  TH3F* effMapPionTPC;
  TH3F* effMapPionTPCTOF;

  void init(InitContext&)
  {
    // Axes
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec sigmassPhiAxis = {phiConfigs.nBinsMPhi, phiConfigs.lowMPhi, phiConfigs.upMPhi, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec massK0SAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec nSigmaPiAxis = {100, -10.0f, 10.0f, "N#sigma #pi"};
    AxisSpec vertexZAxis = {100, -cutZVertex, cutZVertex, "vrtx_{Z} [cm]"};
    AxisSpec etaAxis = {16, -trackConfigs.etaMax, trackConfigs.etaMax, "#eta"};
    AxisSpec yAxis = {deltaYConfigs.nBinsY, -deltaYConfigs.cfgYAcceptance, deltaYConfigs.cfgYAcceptance, "#it{y}"};
    AxisSpec deltayAxis = {deltaYConfigs.nBinsDeltaY, -1.0f, 1.0f, "#Delta#it{y}"};
    AxisSpec deltaphiAxis = {72, -o2::constants::math::PIHalf, o2::constants::math::PIHalf * 3, "#Delta#varphi"};
    AxisSpec phiAxis = {629, 0, o2::constants::math::TwoPI, "#phi"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec pTPhiAxis = {120, 0.0f, 12.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPhiAxis{(std::vector<double>)binspTPhi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTK0SAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPiAxis = {50, 0.0f, 5.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};

    // Histograms
    // Number of events per selection
    dataEventHist.add("hEventSelection", "hEventSelection", kTH1F, {{6, -0.5f, 5.5f}});
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "INEL>0 cut");
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a #phi cand");
    dataEventHist.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "With at least a #phi");

    // Event information
    dataEventHist.add("hVertexZ", "hVertexZ", kTH1F, {vertexZAxis});
    dataEventHist.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {multAxis});
    dataEventHist.add("hMultiplicityPercentWithPhi", "Multiplicity Percentile in Events with a Phi Candidate", kTH1F, {multAxis});
    dataEventHist.add("h2VertexZvsMult", "Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});

    // Eta distribution for dN/deta values estimation in Data
    dataEventHist.add("h5EtaDistribution", "Eta vs multiplicity in Data", kTHnSparseF, {vertexZAxis, binnedmultAxis, etaAxis, phiAxis, {3, -0.5f, 2.5f}});

    // Number of MC events per selection for Rec and Gen
    mcEventHist.add("hRecMCEventSelection", "hRecMCEventSelection", kTH1F, {{9, -0.5f, 8.5f}});
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(2, "kIsTriggerTVX");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(3, "kNoTimeFrameBorder");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(5, "posZ cut");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(6, "INEL>0 cut");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(7, "With at least a gen coll");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(8, "With at least a #phi cand");
    mcEventHist.get<TH1>(HIST("hRecMCEventSelection"))->GetXaxis()->SetBinLabel(9, "With at least a #phi");

    mcEventHist.add("hGenMCEventSelection", "hGenMCEventSelection", kTH1F, {{5, -0.5f, 4.5f}});
    mcEventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    mcEventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(2, "posZ cut");
    mcEventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(3, "INEL>0 cut");
    mcEventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(4, "With at least a #phi");
    mcEventHist.get<TH1>(HIST("hGenMCEventSelection"))->GetXaxis()->SetBinLabel(5, "With at least a reco coll");

    // MC Event information for Rec and Gen
    mcEventHist.add("hRecoMCVertexZ", "hRecoMCVertexZ", kTH1F, {vertexZAxis});
    mcEventHist.add("hUnbinnedRecoMCMultiplicityPercent", "RecoMC Multiplicity Percentile", kTH1F, {multAxis});
    mcEventHist.add("hRecoMCMultiplicityPercent", "RecoMC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    mcEventHist.add("hRecoMCMultiplicityPercentWithPhi", "RecoMC Multiplicity Percentile in Events with a Phi Candidate", kTH1F, {binnedmultAxis});
    mcEventHist.add("h2RecoMCVertexZvsMult", "RecoMC Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    mcEventHist.add("hSplitVertexZ", "Split in z-vtx", kTH1F, {{100, -5.0f, 5.0f}});

    mcEventHist.add("hGenMCVertexZ", "hGenMCVertexZ", kTH1F, {vertexZAxis});
    mcEventHist.add("hGenMCMultiplicityPercent", "GenMC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    mcEventHist.add("hGenMCAssocRecoMultiplicityPercent", "GenMC AssocReco Multiplicity Percentile", kTH1F, {binnedmultAxis});
    mcEventHist.add("h2GenMCAssocRecoVertexZvsMult", "GenMC AssocReco Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    mcEventHist.add("hGenMCRecoMultiplicityPercent", "GenMCReco Multiplicity Percentile", kTH1F, {binnedmultAxis});

    // Eta distribution for dN/deta values estimation in MC
    mcEventHist.add("h6RecoMCEtaDistribution", "Eta vs multiplicity in MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, etaAxis, phiAxis, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}});

    mcEventHist.add("h5GenMCEtaDistribution", "Eta vs multiplicity in MCGen", kTHnSparseF, {binnedmultAxis, etaAxis, phiAxis, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}});
    mcEventHist.add("h6GenMCAssocRecoEtaDistribution", "Eta vs multiplicity in MCGen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, etaAxis, phiAxis, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}});
    mcEventHist.add("h6GenMCAllAssocRecoEtaDistribution", "Eta vs multiplicity in MCGen Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, etaAxis, phiAxis, {6, -0.5f, 5.5f}, {3, -0.5f, 2.5f}});

    // Phi topological/PID cuts
    dataPhiHist.add("h2DauTracksPhiDCAxyPreCutData", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAzPreCutData", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAxyPostCutData", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAzPostCutData", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    dataPhiHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    dataPhiHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    dataPhiHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});

    if (analysisModeConfigs.isData) {
      // Phi invariant mass for computing purities and normalisation
      dataPhiHist.add("h3PhipurData", "Invariant mass of Phi for Purity (no K0S/Pi) in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

      dataPhiHist.add("h4PhipurK0SData", "Invariant mass of Phi for Purity (K0S) in Data", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      dataPhiHist.get<THnSparse>(HIST("h4PhipurK0SData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        dataPhiHist.get<THnSparse>(HIST("h4PhipurK0SData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      dataPhiHist.add("h4PhipurPiData", "Invariant mass of Phi for Purity (Pi) in Data", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      dataPhiHist.get<THnSparse>(HIST("h4PhipurPiData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        dataPhiHist.get<THnSparse>(HIST("h4PhipurPiData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }
    }

    // DCA plots for phi daughters in MCReco
    mcPhiHist.add("h2DauTracksPhiDCAxyPreCutMCReco", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAzPreCutMCReco", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAxyPostCutMCReco", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAzPostCutMCReco", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    if (analysisModeConfigs.isClosure) {
      // MCPhi invariant mass for computing purities
      closureMCPhiHist.add("h3PhipurMCClosure", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

      closureMCPhiHist.add("h4PhipurK0SMCClosure", "Invariant mass of Phi for Purity (K0S) in MCClosure", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      closureMCPhiHist.get<THnSparse>(HIST("h4PhipurK0SMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        closureMCPhiHist.get<THnSparse>(HIST("h4PhipurK0SMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      closureMCPhiHist.add("h4PhipurPiMCClosure", "Invariant mass of Phi for Purity (Pi) in MCClosure", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      closureMCPhiHist.get<THnSparse>(HIST("h4PhipurPiMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        closureMCPhiHist.get<THnSparse>(HIST("h4PhipurPiMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }
    }

    // K0S topological/PID cuts
    dataK0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    dataK0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    dataK0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    dataK0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});

    if (analysisModeConfigs.isData) {
      // 2D mass of Phi and K0S for Data
      dataPhiK0SHist.add("h5PhiK0SData", "2D Invariant mass of Phi and K0Short for Data", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, sigmassPhiAxis});
      dataPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        dataPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // 1D mass of K0S for Data
      dataPhiK0SHist.add("h3PhiK0SSEIncNew", "Invariant mass of K0Short for Same Event Inclusive", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
      dataPhiK0SHist.add("h3PhiK0SSEFCutNew", "Invariant mass of K0Short for Same Event Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
      dataPhiK0SHist.add("h3PhiK0SSESCutNew", "Invariant mass of K0Short for Same Event Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
    }

    // K0S rapidity in Data
    dataK0SHist.add("h3K0SRapidityData", "K0Short rapidity for Data", kTH3F, {binnedmultAxis, binnedpTK0SAxis, yAxis});

    if (analysisModeConfigs.isMC) {
      // RecMC K0S coupled to Phi
      mcPhiK0SHist.add("h4PhiK0SMCReco", "K0S coupled to Phi in MCReco", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
      mcPhiK0SHist.get<THnSparse>(HIST("h4PhiK0SMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiK0SHist.get<THnSparse>(HIST("h4PhiK0SMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // GenMC K0S coupled to Phi
      mcPhiK0SHist.add("h3PhiK0SMCGen", "K0S coupled toPhi in MCGen", kTH3F, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTK0SAxis});
      mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGen"))->GetXaxis()->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGen"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      mcPhiK0SHist.add("h3PhiK0SMCGenAssocReco", "K0S coupled toPhi in MCGen Associated MCReco Collision", kTH3F, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTK0SAxis});
      mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGenAssocReco"))->GetXaxis()->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGenAssocReco"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }
    }

    if (analysisModeConfigs.isClosure) {
      // 2D mass of Phi and K0S for Closure Test
      closureMCPhiK0SHist.add("h5PhiK0SMCClosure", "2D Invariant mass of Phi and K0Short for MC Closure Test", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, sigmassPhiAxis});
      closureMCPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        closureMCPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // 1D mass of K0S for Closure Test
      closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSEIncNew", "Invariant mass of K0Short for Inclusive for Closure Test", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
      closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSEFCutNew", "Invariant mass of K0Short for Deltay < FirstCut for Closure Test", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
      closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSESCutNew", "Invariant mass of K0Short for Deltay < SecondCut for Closure Test", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});
    }

    if (analysisModeConfigs.isData) {
      // Phi mass vs Pion NSigma dE/dx for Data
      dataPhiPionHist.add("h6PhiPiData", "Phi Invariant mass vs Pion nSigma TPC/TOF for Data", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigmassPhiAxis});
      dataPhiPionHist.get<THnSparse>(HIST("h6PhiPiData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        dataPhiPionHist.get<THnSparse>(HIST("h6PhiPiData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // Pion NSigma dE/dx for Data
      dataPhiPionHist.add("h4PhiPiSEIncNew", "Pion nSigma TPC/TOF for Same Event Inclusive", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
      dataPhiPionHist.add("h4PhiPiSEFCutNew", "Pion nSigma TPC/TOF for Same Event Deltay < FirstCut", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
      dataPhiPionHist.add("h4PhiPiSESCutNew", "Pion nSigma TPC/TOF for Same Event Deltay < SecondCut", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    }

    // Pion rapidity in Data
    dataPionHist.add("h3PiRapidityData", "Pion rapidity for Data", kTH3F, {binnedmultAxis, binnedpTPiAxis, yAxis});

    // DCA plots for pions in Data
    dataPionHist.add("h2TracksPiDCAxyPreCutData", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPionHist.add("h2TracksPiDCAzPreCutData", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    dataPionHist.add("h2TracksPiDCAxyPostCutData", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPionHist.add("h2TracksPiDCAzPostCutData", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    // DCA plots for pions in MCReco
    mcPionHist.add("h2TracksPiDCAxyPreCutMCReco", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h2TracksPiDCAzPreCutMCReco", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    mcPionHist.add("h2TracksPiDCAxyPostCutMCReco", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h2TracksPiDCAzPostCutMCReco", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    // DCA plots for pions in MCReco distinguishing Primaries, Secondaries from Weak Decay and Secondaries from Material
    mcPionHist.add("h3RecMCDCAxyPrimPi", "Dcaxy distribution vs pt for Primary Pions", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h3RecMCDCAxySecWeakDecayPi", "Dcaz distribution vs pt for Secondary Pions from Weak Decay", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h3RecMCDCAxySecMaterialPi", "Dcaxy distribution vs pt for Secondary Pions from Material", kTH2F, {binnedpTPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});

    if (analysisModeConfigs.isMC) {
      // RecMC Pion coupled to Phi with TPC
      mcPhiPionHist.add("h4PhiPiTPCMCReco", "Pion coupled to Phi in MCReco (TPC)", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}});
      mcPhiPionHist.get<THnSparse>(HIST("h4PhiPiTPCMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiPionHist.get<THnSparse>(HIST("h4PhiPiTPCMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // RecMC Pion coupled to Phi with TPC and TOF
      mcPhiPionHist.add("h5PhiPiTPCTOFMCReco", "Pion coupled to Phi in MCReco (TPC and TOF)", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
      mcPhiPionHist.get<THnSparse>(HIST("h5PhiPiTPCTOFMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiPionHist.get<THnSparse>(HIST("h5PhiPiTPCTOFMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      mcPhiPionHist.add("h3PhiPiMCGen", "Pion coupled to Phi in MCGen", kTH3F, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis});
      mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGen"))->GetXaxis()->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGen"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      mcPhiPionHist.add("h3PhiPiMCGenAssocReco", "Pion coupled to Phi in MCGen Associated Reco Collision", kTH3F, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis});
      mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGenAssocReco"))->GetXaxis()->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGenAssocReco"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }
    }

    if (analysisModeConfigs.isClosure) {
      // Phi mass vs Pion NSigma dE/dx for Closure Test
      closureMCPhiPionHist.add("h6PhiPiMCClosure", "Phi Invariant mass vs Pion nSigma TPC/TOF for MC Closure Test", kTHnSparseF, {{static_cast<int>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1), -0.5f, static_cast<float>(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f)}, binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigmassPhiAxis});
      closureMCPhiPionHist.get<THnSparse>(HIST("h6PhiPiMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
        closureMCPhiPionHist.get<THnSparse>(HIST("h6PhiPiMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", deltaYConfigs.cfgDeltaYAcceptanceBins->at(i)));
      }

      // Phi mass vs Pion NSigma dE/dx for Closure Test
      closureMCPhiPionHist.add("h4ClosureMCPhiPiSEIncNew", "Pion nSigma TPC/TOF for Inclusive for Closure Test", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
      closureMCPhiPionHist.add("h4ClosureMCPhiPiSEFCutNew", "Pion nSigma TPC/TOF for Deltay < FirstCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
      closureMCPhiPionHist.add("h4ClosureMCPhiPiSESCutNew", "Pion nSigma TPC/TOF for Deltay < SecondCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    }

    if (analysisModeConfigs.isMC) {
      // MCPhi invariant mass for computing efficiencies and MCnormalisation
      mcPhiHist.add("h2PhieffInvMass", "Invariant mass of Phi for Efficiency (no K0S/Pi)", kTH2F, {binnedmultAxis, massPhiAxis});

      mcPhiHist.add("h3PhieffK0SInvMassInc", "Invariant mass of Phi for Efficiency (K0S) Inclusive", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massPhiAxis});
      mcPhiHist.add("h3PhieffK0SInvMassFCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massPhiAxis});
      mcPhiHist.add("h3PhieffK0SInvMassSCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massPhiAxis});

      mcPhiHist.add("h3PhieffPiInvMassInc", "Invariant mass of Phi for Efficiency (Pi) Inclusive", kTH3F, {binnedmultAxis, binnedpTPiAxis, massPhiAxis});
      mcPhiHist.add("h3PhieffPiInvMassFCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedpTPiAxis, massPhiAxis});
      mcPhiHist.add("h3PhieffPiInvMassSCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedpTPiAxis, massPhiAxis});

      // GenMC Phi and Phi coupled to K0S and Pion
      mcPhiHist.add("h1PhiGenMC", "Phi for GenMC", kTH1F, {binnedmultAxis});
      mcPhiHist.add("h1PhiGenMCAssocReco", "Phi for GenMC Associated Reco Collision", kTH1F, {binnedmultAxis});

      mcPhiHist.add("h2PhieffK0SGenMCInc", "Phi coupled to K0Short for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedpTK0SAxis});
      mcPhiHist.add("h2PhieffK0SGenMCFCut", "Phi coupled to K0Short for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedpTK0SAxis});
      mcPhiHist.add("h2PhieffK0SGenMCSCut", "Phi coupled to K0Short for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedpTK0SAxis});

      mcPhiHist.add("h2PhieffK0SGenMCIncAssocReco", "Phi coupled to K0Short for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTK0SAxis});
      mcPhiHist.add("h2PhieffK0SGenMCFCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTK0SAxis});
      mcPhiHist.add("h2PhieffK0SGenMCSCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTK0SAxis});

      mcPhiHist.add("h2PhieffPiGenMCInc", "Phi coupled to Pion for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedpTPiAxis});
      mcPhiHist.add("h2PhieffPiGenMCFCut", "Phi coupled to Pion for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedpTPiAxis});
      mcPhiHist.add("h2PhieffPiGenMCSCut", "Phi coupled to Pion for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedpTPiAxis});

      mcPhiHist.add("h2PhieffPiGenMCIncAssocReco", "Phi coupled to Pion for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTPiAxis});
      mcPhiHist.add("h2PhieffPiGenMCFCutAssocReco", "Phi coupled to Pion for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTPiAxis});
      mcPhiHist.add("h2PhieffPiGenMCSCutAssocReco", "Phi coupled to Pion for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTPiAxis});

      // Rapidity smearing matrix for Phi
      mcPhiHist.add("h3PhiRapiditySmearing", "Rapidity Smearing Matrix for Phi", kTH3F, {binnedmultAxis, yAxis, yAxis});

      // MCK0S invariant mass and GenMC K0S for computing efficiencies
      mcK0SHist.add("h3K0SMCReco", "K0S for MCReco", kTH3F, {binnedmultAxis, binnedpTK0SAxis, massK0SAxis});

      mcK0SHist.add("h2K0SMCGen", "K0S for MCGen", kTH2F, {binnedmultAxis, binnedpTK0SAxis});
      mcK0SHist.add("h2K0SMCGenAssocReco", "K0S for MCGen Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTK0SAxis});

      // Rapidity smearing matrix for K0S and rapidity in GenMC
      mcK0SHist.add("h4K0SRapiditySmearing", "Rapidity Smearing Matrix for K0Short", kTHnSparseF, {binnedmultAxis, binnedpTK0SAxis, yAxis, yAxis});

      mcK0SHist.add("h3K0SRapidityGenMC", "Rapidity for K0Short for GenMC", kTH3F, {binnedmultAxis, binnedpTK0SAxis, yAxis});

      // MCPion invariant mass and GenMC Pion for computing efficiencies
      mcPionHist.add("h3PiTPCMCReco", "Pion for MCReco (TPC)", kTH3F, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}});
      mcPionHist.add("h4PiTPCTOFMCReco", "Pion for MCReco (TPC and TOF)", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});

      mcPionHist.add("h2PiMCGen", "Pion for GenMC", kTH2F, {binnedmultAxis, binnedpTPiAxis});
      mcPionHist.add("h2PiMCGenAssocReco", "Pion for GenMC Associated Reco Collision", kTH2F, {binnedmultAxis, binnedpTPiAxis});

      // Rapidity smearing matrix for Pion and rapidity in GenMC
      mcPionHist.add("h4PiRapiditySmearing", "Rapidity Smearing Matrix for Pion", kTHnSparseF, {binnedmultAxis, binnedpTPiAxis, yAxis, yAxis});

      mcPionHist.add("h3PiRapidityGenMC", "Rapidity for Pion for GenMC", kTH3F, {binnedmultAxis, binnedpTPiAxis, yAxis});
    }

    if (analysisModeConfigs.isDataNewProc) {
      // Histograms for new analysis procedure (to be finalized and renamed deleting other histograms)
      dataPhiHist.add("h3PhiDataNewProc", "Invariant mass of Phi in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      dataPhiK0SHist.add("h5PhiK0SDataNewProc", "2D Invariant mass of Phi and K0Short in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, massPhiAxis});
      dataPhiPionHist.add("h5PhiPiTPCDataNewProc", "Phi Invariant mass vs Pion nSigma TPC in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
      dataPhiPionHist.add("h5PhiPiTOFDataNewProc", "Phi Invariant mass vs Pion nSigma TOF in Data", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});

      dataPhiK0SHist.add("h5PhiK0SData2PartCorr", "Deltay vs deltaphi for Phi and K0Short in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTK0SAxis, deltayAxis, deltaphiAxis});
      dataPhiPionHist.add("h5PhiPiData2PartCorr", "Deltay vs deltaphi for Phi and Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTPiAxis, deltayAxis, deltaphiAxis});
    }

    if (analysisModeConfigs.isClosureNewProc) {
      closureMCPhiHist.add("h3PhiMCClosureNewProc", "Invariant mass of Phi in MC Closure test", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
      closureMCPhiK0SHist.add("h5PhiK0SMCClosureNewProc", "2D Invariant mass of Phi and K0Short in MC Closure Test", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, massPhiAxis});
      closureMCPhiPionHist.add("h5PhiPiTPCMCClosureNewProc", "Phi Invariant mass vs Pion nSigma TPC in MC Closure Test", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
      closureMCPhiPionHist.add("h5PhiPiTOFMCClosureNewProc", "Phi Invariant mass vs Pion nSigma TOF in MC Closure Test", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
    }

    if (analysisModeConfigs.isMENewProc) {
      mePhiK0SHist.add("h5PhiK0SMENewProc", "2D Invariant mass of Phi and K0Short in Mixed Event", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTK0SAxis, massK0SAxis, massPhiAxis});
      mePhiPionHist.add("h5PhiPiTPCMENewProc", "Phi Invariant mass vs Pion nSigma TPC in Mixed Event", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
      mePhiPionHist.add("h5PhiPiTOFMENewProc", "Phi Invariant mass vs Pion nSigma TOF in Mixed Event", kTHnSparseF, {deltayAxis, binnedmultAxis, binnedpTPiAxis, nSigmaPiAxis, massPhiAxis});
    }

    if (analysisModeConfigs.isMCNewProc) {
      mcPhiHist.add("h4PhiMCRecoNewProc", "Phi in MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTPhiAxis, yAxis});
      mcK0SHist.add("h4K0SMCRecoNewProc", "K0S in MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTK0SAxis, yAxis});
      mcPionHist.add("h4PiMCRecoNewProc", "Pion in MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTPiAxis, yAxis});
      mcPionHist.add("h4PiMCReco2NewProc", "Pion in MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTPiAxis, yAxis});

      mcPhiHist.add("h3PhiMCGenNewProc", "Phi in MCGen", kTH3F, {binnedmultAxis, pTPhiAxis, yAxis});
      mcK0SHist.add("h3K0SMCGenNewProc", "K0S in MCGen", kTH3F, {binnedmultAxis, pTK0SAxis, yAxis});
      mcPionHist.add("h3PiMCGenNewProc", "Pion in MCGen", kTH3F, {binnedmultAxis, pTPiAxis, yAxis});

      mcPhiHist.add("h4PhiMCGenAssocRecoNewProc", "Phi in MCGen Associated MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTPhiAxis, yAxis});
      mcK0SHist.add("h4K0SMCGenAssocRecoNewProc", "K0S in MCGen Associated MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTK0SAxis, yAxis});
      mcPionHist.add("h4PiMCGenAssocRecoNewProc", "Pion in MCGen Associated MCReco", kTHnSparseF, {vertexZAxis, binnedmultAxis, pTPiAxis, yAxis});
    }

    // Initialize CCDB only if purity or efficiencies are requested in the task
    if (useCCDB) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      if (fillMethodSingleWeight)
        getPhiPurityFunctionsFromCCDB();

      if (applyEfficiency) {
        getEfficiencyMapsFromCCDB();
      } else {
        effMapPhi = nullptr;
        effMapK0S = nullptr;
        effMapPionTPC = nullptr;
        effMapPionTPCTOF = nullptr;
      }
    }
  }

  // Event selection and QA filling
  template <bool isMC, typename T>
  bool acceptEventQA(const T& collision, bool QA)
  {
    if constexpr (!isMC) { // data event
      if (QA)
        dataEventHist.fill(HIST("hEventSelection"), 0); // all collisions
      if (!collision.sel8())
        return false;
      if (QA)
        dataEventHist.fill(HIST("hEventSelection"), 1); // sel8 collisions
      if (std::abs(collision.posZ()) >= cutZVertex)
        return false;
      if (QA) {
        dataEventHist.fill(HIST("hEventSelection"), 2); // vertex-Z selected
        dataEventHist.fill(HIST("hVertexZ"), collision.posZ());
      }
      if (!collision.isInelGt0())
        return false;
      if (QA)
        dataEventHist.fill(HIST("hEventSelection"), 3); // INEL>0 collisions
      return true;
    } else { // RecMC event
      if (QA)
        mcEventHist.fill(HIST("hRecMCEventSelection"), 0); // all collisions
      if (!collision.selection_bit(aod::evsel::kIsTriggerTVX))
        return false;
      if (QA)
        mcEventHist.fill(HIST("hRecMCEventSelection"), 1); // kIsTriggerTVX collisions
      if (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        return false;
      if (QA)
        mcEventHist.fill(HIST("hRecMCEventSelection"), 2); // kNoTimeFrameBorder collisions
      if (cfgiskNoITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        return false;
      if (QA)
        mcEventHist.fill(HIST("hRecMCEventSelection"), 3); // kNoITSROFrameBorder collisions (by default not requested by the selection)
      if (std::abs(collision.posZ()) > cutZVertex)
        return false;
      if (QA) {
        mcEventHist.fill(HIST("hRecMCEventSelection"), 4); // vertex-Z selected
        mcEventHist.fill(HIST("hRecoMCVertexZ"), collision.posZ());
      }
      if (!collision.isInelGt0())
        return false;
      if (QA)
        mcEventHist.fill(HIST("hRecMCEventSelection"), 5); // INEL>0 collisions
      return true;
    }
  }

  // Single track selection for strangeness sector
  template <typename T>
  bool selectionTrackStrangeness(const T& track)
  {
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC)
      return false;

    if (std::abs(track.eta()) > trackConfigs.etaMax)
      return false;
    return true;
  }

  // V0 selection
  template <typename T1, typename T2>
  bool selectionV0(const T1& v0, const T2& daughter1, const T2& daughter2)
  {
    if (!selectionTrackStrangeness(daughter1) || !selectionTrackStrangeness(daughter2))
      return false;

    if (v0.v0cosPA() < v0Configs.v0SettingCosPA)
      return false;
    if (v0.v0radius() < v0Configs.v0SettingRadius)
      return false;
    if (v0.pt() < v0Configs.v0SettingMinPt)
      return false;

    if (v0Configs.cfgisV0ForData) {
      if (std::abs(daughter1.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCPion)
        return false;
      if (std::abs(daughter2.tpcNSigmaPi()) > trackConfigs.nSigmaCutTPCPion)
        return false;
    }
    return true;
  }

  // Further V0 selection
  template <typename T1, typename T2>
  bool furtherSelectionV0(const T1& v0, const T2& collision)
  {
    if (v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massK0S > v0Configs.ctauK0s)
      return false;
    if (v0.qtarm() < (v0Configs.paramArmenterosCut * std::abs(v0.alpha())))
      return false;
    if (std::abs(v0.mLambda() - massLambda) < v0Configs.v0rejK0s)
      return false;
    return true;
  }

  // Topological track selection
  template <bool isMC, typename T>
  bool selectionTrackResonance(const T& track, bool doQA)
  {
    if (trackConfigs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (trackConfigs.cfgPVContributor && !track.isPVContributor())
      return false;

    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinKaonPtcut)
      return false;

    if (doQA) {
      if constexpr (!isMC) {
        dataPhiHist.fill(HIST("h2DauTracksPhiDCAxyPreCutData"), track.pt(), track.dcaXY());
        dataPhiHist.fill(HIST("h2DauTracksPhiDCAzPreCutData"), track.pt(), track.dcaZ());
      } else {
        mcPhiHist.fill(HIST("h2DauTracksPhiDCAxyPreCutMCReco"), track.pt(), track.dcaXY());
        mcPhiHist.fill(HIST("h2DauTracksPhiDCAzPreCutMCReco"), track.pt(), track.dcaZ());
      }
    }
    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPV1Phi + (trackConfigs.cMaxDCArToPV2Phi / std::pow(track.pt(), trackConfigs.cMaxDCArToPV3Phi)))
      return false;
    if (doQA) {
      if constexpr (!isMC) {
        dataPhiHist.fill(HIST("h2DauTracksPhiDCAxyPostCutData"), track.pt(), track.dcaXY());
        dataPhiHist.fill(HIST("h2DauTracksPhiDCAzPostCutData"), track.pt(), track.dcaZ());
      } else {
        mcPhiHist.fill(HIST("h2DauTracksPhiDCAxyPostCutMCReco"), track.pt(), track.dcaXY());
        mcPhiHist.fill(HIST("h2DauTracksPhiDCAzPostCutMCReco"), track.pt(), track.dcaZ());
      }
    }
    if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
      return false;
    return true;
  }

  // PIDKaon track selection
  template <typename T>
  bool selectionPIDKaon(const T& track)
  {
    if (!trackConfigs.isNoTOF && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2))
      return true;
    if (!trackConfigs.isNoTOF && !track.hasTOF() && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    if (trackConfigs.isNoTOF && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    return false;
  }

  template <typename T>
  bool selectionPIDKaonpTdependent(const T& track)
  {
    if (track.pt() < trackConfigs.pTToUseTOF && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    if (track.pt() >= trackConfigs.pTToUseTOF && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2))
      return true;
    return false;
  }

  // Reconstruct the Phi
  template <typename T1, typename T2>
  ROOT::Math::PxPyPzMVector recMother(const T1& track1, const T2& track2, float masscand1, float masscand2)
  {
    ROOT::Math::PxPyPzMVector daughter1(track1.px(), track1.py(), track1.pz(), masscand1); // set the daughter1 4-momentum
    ROOT::Math::PxPyPzMVector daughter2(track2.px(), track2.py(), track2.pz(), masscand2); // set the daughter2 4-momentum
    ROOT::Math::PxPyPzMVector mother = daughter1 + daughter2;                              // calculate the mother 4-momentum

    return mother;
  }

  // Topological selection for pions
  template <bool isTOFChecked, bool isMC, typename T>
  bool selectionPion(const T& track, bool doQA)
  {
    if (!track.isGlobalTrackWoDCA())
      return false;

    if (track.itsNCls() < trackConfigs.minITSnCls)
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinPionPtcut)
      return false;

    if constexpr (isTOFChecked) {
      if (track.pt() >= trackConfigs.pTToUseTOF && !track.hasTOF())
        return false;
    }

    if (doQA) {
      if constexpr (!isMC) {
        dataPionHist.fill(HIST("h2TracksPiDCAxyPreCutData"), track.pt(), track.dcaXY());
        dataPionHist.fill(HIST("h2TracksPiDCAzPreCutData"), track.pt(), track.dcaZ());
      } else {
        mcPionHist.fill(HIST("h2TracksPiDCAxyPreCutMCReco"), track.pt(), track.dcaXY());
        mcPionHist.fill(HIST("h2TracksPiDCAzPreCutMCReco"), track.pt(), track.dcaZ());
      }
    }
    if (std::abs(track.dcaXY()) > trackConfigs.cMaxDCArToPV1Pion + (trackConfigs.cMaxDCArToPV2Pion / std::pow(track.pt(), trackConfigs.cMaxDCArToPV3Pion)))
      return false;
    if (doQA) {
      if constexpr (!isMC) {
        dataPionHist.fill(HIST("h2TracksPiDCAxyPostCutData"), track.pt(), track.dcaXY());
        dataPionHist.fill(HIST("h2TracksPiDCAzPostCutData"), track.pt(), track.dcaZ());
      } else {
        mcPionHist.fill(HIST("h2TracksPiDCAxyPostCutMCReco"), track.pt(), track.dcaXY());
        mcPionHist.fill(HIST("h2TracksPiDCAzPostCutMCReco"), track.pt(), track.dcaZ());
      }
    }
    if (trackConfigs.cfgIsDCAzParameterized) {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPV1Pion + (trackConfigs.cMaxDCAzToPV2Pion / std::pow(track.pt(), trackConfigs.cMaxDCAzToPV3Pion)))
        return false;
    } else {
      if (std::abs(track.dcaZ()) > trackConfigs.cMaxDCAzToPVcut)
        return false;
    }
    return true;
  }

  template <typename T1, typename T2>
  bool eventHasRecoPhi(const T1& posTracks, const T2& negTracks)
  {
    int nPhi = 0;

    for (const auto& track1 : posTracks) {
      if (!selectionTrackResonance<false>(track1, false) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      for (const auto& track2 : negTracks) {
        if (!selectionTrackResonance<false>(track2, false) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        nPhi++;
      }
    }

    if (nPhi > 0)
      return true;
    return false;
  }

  template <typename T1, typename T2, typename T3>
  bool eventHasRecoPhiWPDG(const T1& posTracks, const T2& negTracks, const T3& mcParticles)
  {
    int nPhi = 0;

    for (const auto& track1 : posTracks) {
      if (!selectionTrackResonance<true>(track1, false) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      if (!track1.has_mcParticle())
        continue;
      auto mcTrack1 = mcParticles.rawIteratorAt(track1.mcParticleId());
      if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
        continue;

      for (const auto& track2 : negTracks) {
        if (!selectionTrackResonance<true>(track2, false) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        if (!track2.has_mcParticle())
          continue;
        auto mcTrack2 = mcParticles.rawIteratorAt(track2.mcParticleId());
        if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
          continue;

        const auto mcTrack1MotherIndexes = mcTrack1.mothersIds();
        const auto mcTrack2MotherIndexes = mcTrack2.mothersIds();

        float pTMother = -1.0f;
        float yMother = -1.0f;
        bool isMCMotherPhi = false;

        for (const auto& mcTrack1MotherIndex : mcTrack1MotherIndexes) {
          for (const auto& mcTrack2MotherIndex : mcTrack2MotherIndexes) {
            if (mcTrack1MotherIndex != mcTrack2MotherIndex)
              continue;

            const auto mother = mcParticles.rawIteratorAt(mcTrack1MotherIndex);
            if (mother.pdgCode() != o2::constants::physics::Pdg::kPhi)
              continue;

            pTMother = mother.pt();
            yMother = mother.y();
            isMCMotherPhi = true;
          }
        }

        if (!isMCMotherPhi)
          continue;
        if (pTMother < phiConfigs.minPhiPt || std::abs(yMother) > deltaYConfigs.cfgYAcceptance)
          continue;

        nPhi++;
      }
    }

    if (nPhi > 0)
      return true;
    return false;
  }

  template <typename T>
  bool eventHasGenKPair(const T& mcParticles)
  {
    int nKPair = 0;

    for (const auto& mcParticle1 : mcParticles) {
      if (!mcParticle1.isPhysicalPrimary() || std::abs(mcParticle1.eta()) > trackConfigs.etaMax)
        continue;

      for (const auto& mcParticle2 : mcParticles) {
        if (!mcParticle2.isPhysicalPrimary() || std::abs(mcParticle2.eta()) > trackConfigs.etaMax)
          continue;

        if ((mcParticle1.pdgCode() != PDG_t::kKPlus || mcParticle2.pdgCode() != PDG_t::kKMinus) &&
            (mcParticle1.pdgCode() != PDG_t::kKMinus || mcParticle2.pdgCode() != PDG_t::kKPlus))
          continue;

        ROOT::Math::PxPyPzMVector genKPair = recMother(mcParticle1, mcParticle2, massKa, massKa);
        if (genKPair.Pt() < phiConfigs.minPhiPt)
          continue;
        if (genKPair.M() < phiConfigs.lowMPhi || genKPair.M() > phiConfigs.upMPhi)
          continue;
        if (std::abs(genKPair.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        nKPair++;
      }
    }

    if (nKPair > 0)
      return true;
    return false;
  }

  template <typename T>
  bool eventHasGenPhi(const T& mcParticles)
  {
    int nPhi = 0;

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.pdgCode() != o2::constants::physics::Pdg::kPhi)
        continue;
      if (mcParticle.pt() < phiConfigs.minPhiPt)
        continue;
      if (std::abs(mcParticle.y()) > deltaYConfigs.cfgYAcceptance)
        continue;

      auto kDaughters = mcParticle.template daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2)
        continue;
      bool isPosKaon = false, isNegKaon = false;
      for (const auto& kDaughter : kDaughters) {
        if (kDaughter.pdgCode() == PDG_t::kKPlus)
          isPosKaon = true;
        if (kDaughter.pdgCode() == PDG_t::kKMinus)
          isNegKaon = true;
      }
      if (!isPosKaon || !isNegKaon)
        continue;

      nPhi++;
    }

    if (nPhi > 0)
      return true;
    return false;
  }

  template <typename T>
  bool selectionChargedGenParticle(const T& mcParticle)
  {
    if (!mcParticle.isPhysicalPrimary() || std::abs(mcParticle.eta()) > trackConfigs.etaMax)
      return false;

    auto pdgTrack = pdgDB->GetParticle(mcParticle.pdgCode());
    if (pdgTrack == nullptr)
      return false;
    if (std::abs(pdgTrack->Charge()) < trackConfigs.cfgMinAbsCharge)
      return false;

    return true;
  }

  int fromPDGToEnum(int pdgCode)
  {
    int pid = kSpAll;
    switch (std::abs(pdgCode)) {
      case PDG_t::kPiPlus:
        pid = kSpPion;
        break;
      case PDG_t::kKPlus:
        pid = kSpKaon;
        break;
      case PDG_t::kProton:
        pid = kSpProton;
        break;
      default:
        pid = kSpOther;
        break;
    }

    return pid;
  }

  // Get phi-meson purity functions from CCDB
  void getPhiPurityFunctionsFromCCDB()
  {
    TList* listPhiPurityFunctions = ccdb->get<TList>(ccdbPurityPath);
    if (!listPhiPurityFunctions)
      LOG(error) << "Problem getting TList object with phi purity functions!";

    for (size_t multIdx = 0; multIdx < binsMult->size() - 1; multIdx++) {
      for (size_t ptIdx = 0; ptIdx < binspTPhi->size() - 1; ptIdx++) {
        phiPurityFunctions[multIdx][ptIdx] = static_cast<TF1*>(listPhiPurityFunctions->FindObject(Form("funcFitPhiPur_%zu_%zu", multIdx, ptIdx)));
      }
    }
  }

  // Get the phi purity choosing the correct purity function according to the multiplicity and pt of the phi
  double getPhiPurity(float multiplicity, const ROOT::Math::PxPyPzMVector& Phi)
  {
    // Check if multiplicity is out of range
    if (multiplicity < binsMult->front() || multiplicity >= binsMult->back()) {
      LOG(info) << "Multiplicity out of range: " << multiplicity;
      return 0;
    }

    // Find the multiplicity bin using upper_bound which finds the first element strictly greater than 'multiplicity'
    // Subtract 1 to get the correct bin index
    auto multIt = std::upper_bound(binsMult->begin(), binsMult->end(), multiplicity);
    int multIdx = std::distance(binsMult->begin(), multIt) - 1;

    // Check if pT is out of range
    if (Phi.Pt() < binspTPhi->front() || Phi.Pt() >= binspTPhi->back()) {
      LOG(info) << "pT out of range: " << Phi.Pt();
      return 0;
    }

    // Find the pT bin using upper_bound
    // The logic is the same as for multiplicity
    auto pTIt = std::upper_bound(binspTPhi->begin(), binspTPhi->end(), Phi.Pt());
    int pTIdx = std::distance(binspTPhi->begin(), pTIt) - 1;

    return phiPurityFunctions[multIdx][pTIdx]->Eval(Phi.M());
  }

  void getEfficiencyMapsFromCCDB()
  {
    TList* listEfficiencyMaps = ccdb->get<TList>(ccdbEfficiencyPath);
    if (!listEfficiencyMaps)
      LOG(error) << "Problem getting TList object with efficiency maps!";

    effMapPhi = static_cast<TH3F*>(listEfficiencyMaps->FindObject("h3EfficiencyPhi"));
    if (!effMapPhi) {
      LOG(error) << "Problem getting efficiency map for Phi!";
      return;
    }
    effMapK0S = static_cast<TH3F*>(listEfficiencyMaps->FindObject("h3EfficiencyK0S"));
    if (!effMapK0S) {
      LOG(error) << "Problem getting efficiency map for K0S!";
      return;
    }
    effMapPionTPC = static_cast<TH3F*>(listEfficiencyMaps->FindObject("h3EfficiencyPionTPC"));
    if (!effMapPionTPC) {
      LOG(error) << "Problem getting efficiency map for Pion with TPC!";
      return;
    }
    effMapPionTPCTOF = static_cast<TH3F*>(listEfficiencyMaps->FindObject("h3EfficiencyPionTPCTOF"));
    if (!effMapPionTPCTOF) {
      LOG(error) << "Problem getting efficiency map for Pion with TPC and TOF!";
      return;
    }
  }

  // Fill 2D invariant mass histogram for V0 and Phi
  template <bool isMC, typename T>
  void fillInvMass2D(const T& V0, const std::vector<ROOT::Math::PxPyPzMVector>& listPhi, float multiplicity, const std::vector<float>& weights)
  {
    for (const auto& Phi : listPhi) {
      if constexpr (!isMC) { // same event
        dataPhiK0SHist.fill(HIST("h5PhiK0SData"), 0, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(V0.yK0Short() - Phi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          dataPhiK0SHist.fill(HIST("h5PhiK0SData"), i + 1, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(i + 1));
        }
      } else { // MC event
        closureMCPhiK0SHist.fill(HIST("h5PhiK0SMCClosure"), 0, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(V0.yK0Short() - Phi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          closureMCPhiK0SHist.fill(HIST("h5PhiK0SMCClosure"), i + 1, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(i + 1));
        }
      }
    }
  }

  // Fill Phi invariant mass vs Pion nSigmaTPC/TOF histogram
  template <bool isMC, typename T>
  void fillInvMassNSigma(const T& Pi, const std::vector<ROOT::Math::PxPyPzMVector>& listPhi, float multiplicity, const std::vector<float>& weights)
  {
    float nSigmaTOFPi = (Pi.hasTOF() ? Pi.tofNSigmaPi() : -999);

    for (const auto& Phi : listPhi) {
      if constexpr (!isMC) { // same event
        dataPhiPionHist.fill(HIST("h6PhiPiData"), 0, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(0));
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          dataPhiPionHist.fill(HIST("h6PhiPiData"), i + 1, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(i + 1));
        }
      } else { // MC event
        closureMCPhiPionHist.fill(HIST("h6PhiPiMCClosure"), 0, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(0));
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          closureMCPhiPionHist.fill(HIST("h6PhiPiMCClosure"), i + 1, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(i + 1));
        }
      }
    }
  }

  // Fill invariant mass histogram for V0
  template <bool isMC, typename T>
  void fillInvMass(const T& V0, float multiplicity, const std::vector<float>& weights)
  {
    if constexpr (!isMC) { // same event
      dataPhiK0SHist.fill(HIST("h3PhiK0SSEIncNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(0));
      dataPhiK0SHist.fill(HIST("h3PhiK0SSEFCutNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(1));
      dataPhiK0SHist.fill(HIST("h3PhiK0SSESCutNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(2));
    } else { // MC event
      closureMCPhiK0SHist.fill(HIST("h3ClosureMCPhiK0SSEIncNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(0));
      closureMCPhiK0SHist.fill(HIST("h3ClosureMCPhiK0SSEFCutNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(1));
      closureMCPhiK0SHist.fill(HIST("h3ClosureMCPhiK0SSESCutNew"), multiplicity, V0.pt(), V0.mK0Short(), weights.at(2));
    }
  }

  // Fill nSigmaTPC/TOF histogram for Pion
  template <bool isMC, typename T>
  void fillNSigma(const T& Pi, float multiplicity, const std::vector<float>& weights)
  {
    float nSigmaTOFPi = (Pi.hasTOF() ? Pi.tofNSigmaPi() : -999);

    if constexpr (!isMC) { // same event
      dataPhiPionHist.fill(HIST("h4PhiPiSEIncNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(0));
      dataPhiPionHist.fill(HIST("h4PhiPiSEFCutNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(1));
      dataPhiPionHist.fill(HIST("h4PhiPiSESCutNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(2));
    } else { // MC event
      closureMCPhiPionHist.fill(HIST("h4ClosureMCPhiPiSEIncNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(0));
      closureMCPhiPionHist.fill(HIST("h4ClosureMCPhiPiSEFCutNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(1));
      closureMCPhiPionHist.fill(HIST("h4ClosureMCPhiPiSESCutNew"), multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, weights.at(2));
    }
  }

  void processPhiPurityData(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    // Check if the event selection is passed
    if (!acceptEventQA<false>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    dataEventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;
    bool isFilledhV0 = false;

    double weight{1.0};

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<false>(track1, true) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      dataPhiHist.fill(HIST("hEta"), track1.eta());
      dataPhiHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      dataPhiHist.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());

      auto track1ID = track1.globalIndex();

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<false>(track2, true) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          dataEventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event
          dataEventHist.fill(HIST("hMultiplicityPercentWithPhi"), multiplicity);
          isCountedPhi = true;
        }

        if (fillMethodSingleWeight)
          weight *= (1 - getPhiPurity(multiplicity, recPhi));

        dataPhiHist.fill(HIST("h3PhipurData"), multiplicity, recPhi.Pt(), recPhi.M());

        std::vector<int> countsK0S(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (!isFilledhV0) {
            dataK0SHist.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
            dataK0SHist.fill(HIST("hV0CosPA"), v0.v0cosPA());

            // Filling the PID of the V0 daughters in the region of the K0 peak
            if (v0Configs.lowMK0S < v0.mK0Short() && v0.mK0Short() < v0Configs.upMK0S) {
              dataK0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi());
              dataK0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
            }
          }

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsK0S.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsK0S.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsK0S.at(i) > 0)
            dataPhiHist.fill(HIST("h4PhipurK0SData"), i, multiplicity, recPhi.Pt(), recPhi.M());
        }

        isFilledhV0 = true;

        std::vector<int> countsPi(deltaYConfigs.cfgDeltaYAcceptanceBins->size(), 0);

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion<true, false>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsPi.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsPi.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsPi.at(i) > 0)
            dataPhiHist.fill(HIST("h4PhipurPiData"), i, multiplicity, recPhi.Pt(), recPhi.M());
        }
      }
    }

    weight = 1 - weight;
    dataEventHist.fill(HIST("hEventSelection"), 5, weight); // at least a Phi in the event
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPurityData, "Process function for Phi Purities in Data", true);

  void processPhiK0SData(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const&, FullV0s const& V0s, V0DauTracks const&)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    dataEventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

      // Cut on V0 dynamic columns
      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;
      if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
        continue;

      dataK0SHist.fill(HIST("h3K0SRapidityData"), multiplicity, v0.pt(), v0.yK0Short());

      if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 1);

      // Phi reconstruction
      // Loop over positive tracks
      for (const auto& track1 : posThisColl) {
        if (!selectionTrackResonance<false>(track1, false) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative tracks
        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance<false>(track2, false) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
            continue;
          if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(multiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMass2D<false>(v0, listrecPhi, multiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillInvMass<false>(v0, multiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SData, "Process function for Phi-K0S Correlations in Data", false);

  void processPhiPionData(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const& fullTracks)
  {
    if (!collision.isInelGt0())
      return;

    float multiplicity = collision.centFT0M();
    dataEventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullTracks) {

      // Pion selection
      if (!selectionPion<true, false>(track, true))
        continue;

      dataPionHist.fill(HIST("h3PiRapidityData"), multiplicity, track.pt(), track.rapidity(massPi));

      if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 1);

      // Phi reconstruction
      // Loop over positive tracks
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance<false>(track1, false) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        // Loop over all negative tracks
        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance<false>(track2, false) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
            continue;
          if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(multiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMassNSigma<false>(track, listrecPhi, multiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillNSigma<false>(track, multiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPionData, "Process function for Phi-Pion Correlations in Data", false);

  void processPhiMCReco(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    mcEventHist.fill(HIST("hUnbinnedRecoMCMultiplicityPercent"), multiplicity);

    if (!collision.has_mcCollision())
      return;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hRecoMCMultiplicityPercent"), genmultiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<true>(track1, false) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      if (!track1.has_mcParticle())
        continue;
      auto mcTrack1 = track1.mcParticle_as<aod::McParticles>();
      if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
        continue;

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<true>(track2, false) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        if (!track2.has_mcParticle())
          continue;
        auto mcTrack2 = track2.mcParticle_as<aod::McParticles>();
        if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
          continue;

        bool isMCMotherPhi = false;
        auto mcMotherPhi = mcTrack1.mothers_as<aod::McParticles>()[0];
        for (const auto& MotherOfmcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
          for (const auto& MotherOfmcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
            if (MotherOfmcTrack1 == MotherOfmcTrack2 && MotherOfmcTrack1.pdgCode() == o2::constants::physics::Pdg::kPhi) {
              mcMotherPhi = MotherOfmcTrack1;
              isMCMotherPhi = true;
            }
          }
        }

        if (!isMCMotherPhi)
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
          continue;

        mcPhiHist.fill(HIST("h3PhiRapiditySmearing"), genmultiplicity, recPhi.Rapidity(), mcMotherPhi.y());

        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          mcEventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi in the event
          isCountedPhi = true;
        }

        mcPhiHist.fill(HIST("h2PhieffInvMass"), genmultiplicity, recPhi.M());

        std::array<bool, 3> isCountedK0S{false, false, false};

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          if (!v0.has_mcParticle()) {
            continue;
          }

          auto v0mcparticle = v0.mcParticle();
          if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
            continue;

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;
          if (!isCountedK0S.at(0)) {
            mcPhiHist.fill(HIST("h3PhieffK0SInvMassInc"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(0) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgFCutOnDeltaY)
            continue;
          if (!isCountedK0S.at(1)) {
            mcPhiHist.fill(HIST("h3PhieffK0SInvMassFCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(1) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgSCutOnDeltaY)
            continue;
          if (!isCountedK0S.at(2)) {
            mcPhiHist.fill(HIST("h3PhieffK0SInvMassSCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(2) = true;
          }
        }

        std::array<bool, 3> isCountedPi{false, false, false};

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (!track.has_mcParticle())
            continue;

          auto mcTrack = track.mcParticle_as<aod::McParticles>();
          if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus || !mcTrack.isPhysicalPrimary())
            continue;

          if (!selectionPion<true, true>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;
          if (!isCountedPi.at(0)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassInc"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(0) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgFCutOnDeltaY)
            continue;
          if (!isCountedPi.at(1)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassFCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(1) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgSCutOnDeltaY)
            continue;
          if (!isCountedPi.at(2)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassSCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(2) = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiMCReco, "Process function for Phi in MCReco (to be removed)", false);

  void processPhiK0SMCReco(SimCollisions const& collisions, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!acceptEventQA<true>(collision, false))
        continue;

      if (!collision.has_mcCollision())
        continue;

      const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
      float genmultiplicity = mcCollision.centFT0M();

      // Defining V0s in the collision
      auto v0sThisColl = V0s.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      // Defining McParticles in the collision
      auto mcParticlesThisColl = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

      // V0 already reconstructed by the builder
      for (const auto& v0 : v0sThisColl) {
        if (!v0.has_mcParticle())
          continue;

        auto v0mcparticle = v0.mcParticle();
        if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
          continue;

        const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

        if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
          continue;
        if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
          continue;

        mcK0SHist.fill(HIST("h4K0SRapiditySmearing"), genmultiplicity, v0.pt(), v0.yK0Short(), v0mcparticle.y());

        if (std::abs(v0mcparticle.y()) > deltaYConfigs.cfgYAcceptance)
          continue;

        mcK0SHist.fill(HIST("h3K0SMCReco"), genmultiplicity, v0mcparticle.pt(), v0.mK0Short());

        std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != o2::constants::physics::Pdg::kPhi)
            continue;
          if (mcParticle.pt() < phiConfigs.minPhiPt || mcParticle.pt() > phiConfigs.maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > deltaYConfigs.cfgYAcceptance)
            continue;

          counts.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0mcparticle.y() - mcParticle.y()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (counts.at(i) > 0)
            mcPhiK0SHist.fill(HIST("h4PhiK0SMCReco"), i, genmultiplicity, v0mcparticle.pt(), v0.mK0Short());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SMCReco, "Process function for Phi-K0S Correlations Efficiency correction in MCReco", false);

  void processPhiPionMCReco(SimCollisions const& collisions, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      if (!acceptEventQA<true>(collision, false))
        continue;

      if (!collision.has_mcCollision())
        continue;

      const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
      float genmultiplicity = mcCollision.centFT0M();

      // Defining tracks in the collision
      auto mcTracksThisColl = fullMCTracks.sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      // Defining McParticles in the collision
      auto mcParticlesThisColl = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

      // Loop over all primary pion candidates
      for (const auto& track : mcTracksThisColl) {
        // Pion selection
        if (!selectionPion<false, true>(track, false))
          continue;

        if (!track.has_mcParticle())
          continue;

        auto mcTrack = track.mcParticle_as<aod::McParticles>();
        if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus)
          continue;

        if (std::abs(mcTrack.y()) > deltaYConfigs.cfgYAcceptance)
          continue;

        // Primary pion selection
        if (mcTrack.isPhysicalPrimary()) {
          mcPionHist.fill(HIST("h3RecMCDCAxyPrimPi"), track.pt(), track.dcaXY());
        } else {
          if (mcTrack.getProcess() == 4) { // Selection of secondary pions from weak decay
            mcPionHist.fill(HIST("h3RecMCDCAxySecWeakDecayPi"), track.pt(), track.dcaXY());
          } else { // Selection of secondary pions from material interactions
            mcPionHist.fill(HIST("h3RecMCDCAxySecMaterialPi"), track.pt(), track.dcaXY());
          }
          continue;
        }

        mcPionHist.fill(HIST("h4PiRapiditySmearing"), genmultiplicity, track.pt(), track.rapidity(massPi), mcTrack.y());

        mcPionHist.fill(HIST("h3PiTPCMCReco"), genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi());

        std::vector<int> countsTPC(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != o2::constants::physics::Pdg::kPhi)
            continue;
          if (mcParticle.pt() < phiConfigs.minPhiPt || mcParticle.pt() > phiConfigs.maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsTPC.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(mcTrack.y() - mcParticle.y()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsTPC.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsTPC.at(i) > 0)
            mcPhiPionHist.fill(HIST("h4PhiPiTPCMCReco"), i, genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi());
        }

        if (track.pt() >= trackConfigs.pTToUseTOF && !track.hasTOF())
          continue;

        mcPionHist.fill(HIST("h4PiTPCTOFMCReco"), genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi(), track.tofNSigmaPi());

        std::vector<int> countsTPCTOF(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != o2::constants::physics::Pdg::kPhi)
            continue;
          if (mcParticle.pt() < phiConfigs.minPhiPt || mcParticle.pt() > phiConfigs.maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsTPCTOF.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(mcTrack.y() - mcParticle.y()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsTPCTOF.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsTPCTOF.at(i) > 0)
            mcPhiPionHist.fill(HIST("h5PhiPiTPCTOFMCReco"), i, genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPionMCReco, "Process function for Phi-Pion Correlations Efficiency correction in MCReco", false);

  void processPhiPurityMCClosure(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    if (!collision.has_mcCollision())
      return;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hRecoMCMultiplicityPercent"), genmultiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    double weight{1.0};

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<true>(track1, true) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<true>(track2, true) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          mcEventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi candidate in the event
          mcEventHist.fill(HIST("hRecoMCMultiplicityPercentWithPhi"), genmultiplicity);
          isCountedPhi = true;
        }

        if (fillMethodSingleWeight)
          weight *= (1 - getPhiPurity(genmultiplicity, recPhi));

        closureMCPhiHist.fill(HIST("h3PhipurMCClosure"), genmultiplicity, recPhi.Pt(), recPhi.M());

        std::vector<int> countsK0S(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!v0.has_mcParticle())
              continue;
            auto v0mcparticle = v0.mcParticle();
            if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
              continue;
          }

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsK0S.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsK0S.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsK0S.at(i) > 0)
            closureMCPhiHist.fill(HIST("h4PhipurK0SInvMass"), i, genmultiplicity, recPhi.Pt(), recPhi.M());
        }

        std::vector<int> countsPi(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!track.has_mcParticle())
              continue;
            auto mcTrack = track.mcParticle_as<aod::McParticles>();
            if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus || !mcTrack.isPhysicalPrimary())
              continue;
          }

          if (!selectionPion<true, true>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;

          countsPi.at(0)++;
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsPi.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsPi.at(i) > 0)
            closureMCPhiHist.fill(HIST("h4PhipurPiInvMass"), i, genmultiplicity, recPhi.Pt(), recPhi.M());
        }
      }
    }

    weight = 1 - weight;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 8, weight); // at least a Phi in the event
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPurityMCClosure, "Process function for Phi Purities in MCClosure", false);

  void processPhiK0SMCClosure(SimCollisions::iterator const& collision, FullMCTracks const&, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // V0 already reconstructed by the builder
    for (const auto& v0 : V0s) {
      if (cfgisRecMCWPDGForClosure1) {
        if (!v0.has_mcParticle())
          continue;
        auto v0mcparticle = v0.mcParticle();
        if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
          continue;
      }

      const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;
      if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
        continue;

      if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 1);

      // Phi reconstruction
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance<true>(track1, false) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance<true>(track2, false) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          if (cfgisRecMCWPDGForClosure2) {
            if (!track1.has_mcParticle())
              continue;
            auto mcTrack1 = track1.mcParticle_as<aod::McParticles>();
            if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
              continue;

            if (!track2.has_mcParticle())
              continue;
            auto mcTrack2 = track2.mcParticle_as<aod::McParticles>();
            if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
              continue;

            bool isMCMotherPhi = false;
            for (const auto& motherOfMcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
              for (const auto& motherOfMcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
                if (motherOfMcTrack1.pdgCode() != motherOfMcTrack2.pdgCode())
                  continue;
                if (motherOfMcTrack1.globalIndex() != motherOfMcTrack2.globalIndex())
                  continue;
                if (motherOfMcTrack1.pdgCode() != o2::constants::physics::Pdg::kPhi)
                  continue;
                isMCMotherPhi = true;
              }
            }
            if (!isMCMotherPhi)
              continue;
          }

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
            continue;
          if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(genmultiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMass2D<true>(v0, listrecPhi, genmultiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillInvMass<true>(v0, genmultiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SMCClosure, "Process function for Phi-K0S Correlations in MCClosure", false);

  void processPhiPionMCClosure(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, false))
      return;

    if (!collision.has_mcCollision())
      return;

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Loop over all primary pion candidates
    for (const auto& track : fullMCTracks) {
      if (cfgisRecMCWPDGForClosure1) {
        if (!track.has_mcParticle())
          continue;
        auto mcTrack = track.mcParticle_as<aod::McParticles>();
        if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus || !mcTrack.isPhysicalPrimary())
          continue;
      }

      // Pion selection
      if (!selectionPion<true, true>(track, true))
        continue;

      if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 1);

      // Phi reconstruction
      for (const auto& track1 : posThisColl) { // loop over all selected tracks
        if (!selectionTrackResonance<true>(track1, false) || !selectionPIDKaonpTdependent(track1))
          continue; // topological and PID selection

        auto track1ID = track1.globalIndex();

        for (const auto& track2 : negThisColl) {
          if (!selectionTrackResonance<true>(track2, false) || !selectionPIDKaonpTdependent(track2))
            continue; // topological and PID selection

          auto track2ID = track2.globalIndex();
          if (track2ID == track1ID)
            continue; // condition to avoid double counting of pair

          if (cfgisRecMCWPDGForClosure2) {
            if (!track1.has_mcParticle())
              continue;
            auto mcTrack1 = track1.mcParticle_as<aod::McParticles>();
            if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
              continue;

            if (!track2.has_mcParticle())
              continue;
            auto mcTrack2 = track2.mcParticle_as<aod::McParticles>();
            if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
              continue;

            bool isMCMotherPhi = false;
            for (const auto& motherOfMcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
              for (const auto& motherOfMcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
                if (motherOfMcTrack1.pdgCode() != motherOfMcTrack2.pdgCode())
                  continue;
                if (motherOfMcTrack1.globalIndex() != motherOfMcTrack2.globalIndex())
                  continue;
                if (motherOfMcTrack1.pdgCode() != o2::constants::physics::Pdg::kPhi)
                  continue;
                isMCMotherPhi = true;
              }
            }
            if (!isMCMotherPhi)
              continue;
          }

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
            continue;
          if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(genmultiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMassNSigma<true>(track, listrecPhi, genmultiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillNSigma<true>(track, genmultiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPionMCClosure, "Process function for Phi-Pion Correlations in MCClosure", false);

  void processPhiMCGen(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    mcEventHist.fill(HIST("hGenMCEventSelection"), 0); // all collisions
    if (std::abs(mcCollision.posZ()) > cutZVertex)
      return;
    mcEventHist.fill(HIST("hGenMCEventSelection"), 1); // vertex-Z selected
    mcEventHist.fill(HIST("hGenMCVertexZ"), mcCollision.posZ());
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;
    mcEventHist.fill(HIST("hGenMCEventSelection"), 2); // INEL>0 collisions

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    bool isCountedPhi = false;

    for (const auto& mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != o2::constants::physics::Pdg::kPhi)
        continue;
      auto kDaughters = mcParticle1.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2)
        continue;
      bool isPosKaon = false, isNegKaon = false;
      for (const auto& kDaughter : kDaughters) {
        if (kDaughter.pdgCode() == PDG_t::kKPlus)
          isPosKaon = true;
        if (kDaughter.pdgCode() == PDG_t::kKMinus)
          isNegKaon = true;
      }
      if (!isPosKaon || !isNegKaon)
        continue;
      if (std::abs(mcParticle1.y()) > deltaYConfigs.cfgYAcceptance)
        continue;

      if (!isCountedPhi) {
        mcEventHist.fill(HIST("hGenMCEventSelection"), 3); // at least a Phi in the event
        if (isAssocColl)
          mcEventHist.fill(HIST("hGenMCEventSelection"), 4); // with at least a reco collision
        isCountedPhi = true;
      }

      mcPhiHist.fill(HIST("h1PhiGenMC"), genmultiplicity);
      if (isAssocColl)
        mcPhiHist.fill(HIST("h1PhiGenMCAssocReco"), genmultiplicity);

      std::array<bool, 3> isCountedK0S = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != PDG_t::kK0Short)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;

        if (std::abs(mcParticle2.y()) > deltaYConfigs.cfgYAcceptance)
          continue;
        if (!isCountedK0S.at(0)) {
          mcPhiHist.fill(HIST("h2PhieffK0SGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgFCutOnDeltaY)
          continue;
        if (!isCountedK0S.at(1)) {
          mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgSCutOnDeltaY)
          continue;
        if (!isCountedK0S.at(2)) {
          mcPhiHist.fill(HIST("h2PhieffK0SGenMCSCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffK0SGenMCSCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(2) = true;
        }
      }

      std::array<bool, 3> isCountedPi = {false, false, false};

      for (const auto& mcParticle2 : mcParticles) {
        if (std::abs(mcParticle2.pdgCode()) != PDG_t::kPiPlus)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;

        if (std::abs(mcParticle2.y()) > deltaYConfigs.cfgYAcceptance)
          continue;
        if (!isCountedPi.at(0)) {
          mcPhiHist.fill(HIST("h2PhieffPiGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffPiGenMCIncAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgFCutOnDeltaY)
          continue;
        if (!isCountedPi.at(1)) {
          mcPhiHist.fill(HIST("h2PhieffPiGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffPiGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgSCutOnDeltaY)
          continue;
        if (!isCountedPi.at(2)) {
          mcPhiHist.fill(HIST("h2PhieffPiGenMCSCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffPiGenMCSCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(2) = true;
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiMCGen, "Process function for Phi in MCGen (to be removed)", false);

  void processPhiK0SMCGen(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > cutZVertex)
      return;
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    for (const auto& mcParticle1 : mcParticles) {
      if (mcParticle1.pdgCode() != PDG_t::kK0Short)
        continue;
      if (!mcParticle1.isPhysicalPrimary() || mcParticle1.pt() < v0Configs.v0SettingMinPt)
        continue;

      mcK0SHist.fill(HIST("h3K0SRapidityGenMC"), genmultiplicity, mcParticle1.pt(), mcParticle1.y());

      if (std::abs(mcParticle1.y()) > deltaYConfigs.cfgYAcceptance)
        continue;

      mcK0SHist.fill(HIST("h2K0SMCGen"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        mcK0SHist.fill(HIST("h2K0SMCGenAssocReco"), genmultiplicity, mcParticle1.pt());

      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != o2::constants::physics::Pdg::kPhi)
          continue;
        if (cfgisGenMCForClosure) {
          auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
          if (kDaughters.size() != 2)
            continue;
          bool isPosKaon = false, isNegKaon = false;
          for (const auto& kDaughter : kDaughters) {
            if (kDaughter.pdgCode() == PDG_t::kKPlus)
              isPosKaon = true;
            if (kDaughter.pdgCode() == PDG_t::kKMinus)
              isNegKaon = true;
          }
          if (!isPosKaon || !isNegKaon)
            continue;
        }
        if (mcParticle2.pt() < phiConfigs.minPhiPt || mcParticle2.pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(mcParticle2.y()) > deltaYConfigs.cfgYAcceptance)
          continue;

        counts.at(0)++;
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          counts.at(i + 1)++;
        }
      }

      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
        if (counts.at(i) > 0) {
          mcPhiK0SHist.fill(HIST("h3PhiK0SMCGen"), i, genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            mcPhiK0SHist.fill(HIST("h3PhiK0SMCGenAssocReco"), i, genmultiplicity, mcParticle1.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SMCGen, "Process function for Phi-K0S Correlations Efficiency correction in MCGen", false);

  void processPhiPionMCGen(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
  {
    if (std::abs(mcCollision.posZ()) > cutZVertex)
      return;
    if (!pwglf::isINELgtNmc(mcParticles, 0, pdgDB))
      return;

    bool isAssocColl = false;
    for (const auto& collision : collisions) {
      if (acceptEventQA<true>(collision, false)) {
        isAssocColl = true;
        break;
      }
    }

    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hGenMCMultiplicityPercent"), genmultiplicity);

    for (const auto& mcParticle1 : mcParticles) {
      if (std::abs(mcParticle1.pdgCode()) != PDG_t::kPiPlus)
        continue;
      if (!mcParticle1.isPhysicalPrimary() || mcParticle1.pt() < trackConfigs.cMinPionPtcut)
        continue;

      mcPionHist.fill(HIST("h3PiRapidityGenMC"), genmultiplicity, mcParticle1.pt(), mcParticle1.y());

      if (std::abs(mcParticle1.y()) > deltaYConfigs.cfgYAcceptance)
        continue;

      mcPionHist.fill(HIST("h2PiMCGen"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        mcPionHist.fill(HIST("h2PiMCGenAssocReco"), genmultiplicity, mcParticle1.pt());

      std::vector<int> counts(deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1, 0);

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != o2::constants::physics::Pdg::kPhi)
          continue;
        if (cfgisGenMCForClosure) {
          auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
          if (kDaughters.size() != 2)
            continue;
          bool isPosKaon = false, isNegKaon = false;
          for (const auto& kDaughter : kDaughters) {
            if (kDaughter.pdgCode() == PDG_t::kKPlus)
              isPosKaon = true;
            if (kDaughter.pdgCode() == PDG_t::kKMinus)
              isNegKaon = true;
          }
          if (!isPosKaon || !isNegKaon)
            continue;
        }
        if (mcParticle2.pt() < phiConfigs.minPhiPt || mcParticle2.pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(mcParticle2.y()) > deltaYConfigs.cfgYAcceptance)
          continue;

        counts.at(0)++;
        for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > deltaYConfigs.cfgDeltaYAcceptanceBins->at(i))
            continue;
          counts.at(i + 1)++;
        }
      }

      for (size_t i = 0; i < deltaYConfigs.cfgDeltaYAcceptanceBins->size() + 1; i++) {
        if (counts.at(i) > 0) {
          mcPhiPionHist.fill(HIST("h3PhiPiMCGen"), i, genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            mcPhiPionHist.fill(HIST("h3PhiPiMCGenAssocReco"), i, genmultiplicity, mcParticle1.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPionMCGen, "Process function for Phi-Pion Correlations Efficiency correction in MCGen", false);

  // dN/deta procedure
  void processdNdetaWPhiData(SelCollisions::iterator const& collision, FilteredTracks const& filteredTracks)
  {
    // Check if the event selection is passed
    if (!acceptEventQA<false>(collision, true))
      return;

    auto posThisColl = posFiltTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negFiltTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    // Check if the event contains a phi candidate
    if (!eventHasRecoPhi(posThisColl, negThisColl))
      return;

    dataEventHist.fill(HIST("hMultiplicityPercent"), collision.centFT0M());
    dataEventHist.fill(HIST("h2VertexZvsMult"), collision.posZ(), collision.centFT0M());

    for (const auto& track : filteredTracks) {
      if (trackConfigs.applyExtraPhiCuts && ((track.phi() > trackConfigs.extraPhiCuts->at(0) && track.phi() < trackConfigs.extraPhiCuts->at(1)) ||
                                             track.phi() <= trackConfigs.extraPhiCuts->at(2) || track.phi() >= trackConfigs.extraPhiCuts->at(3)))
        continue;

      dataEventHist.fill(HIST("h5EtaDistribution"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kGlobalplusITSonly);
      if (track.hasTPC()) {
        dataEventHist.fill(HIST("h5EtaDistribution"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kGlobalonly);
      } else {
        dataEventHist.fill(HIST("h5EtaDistribution"), collision.posZ(), collision.centFT0M(), track.eta(), track.phi(), kITSonly);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processdNdetaWPhiData, "Process function for dN/deta values in Data", false);

  void processdNdetaWPhiMC(MCCollisions const& mcCollisions, SimCollisions const& collisions, FilteredMCTracks const& filteredMCTracks, aod::McParticles const& mcParticles)
  {
    std::vector<std::vector<int>> collsGrouped(mcCollisions.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
      collsGrouped[mcCollision.globalIndex()].push_back(collision.globalIndex());
    }

    for (const auto& mcCollision : mcCollisions) {
      auto mcParticlesThisMcColl = mcParticles.sliceBy(preslices.mcPartPerMCCollision, mcCollision.globalIndex());

      if (!pwglf::isINELgtNmc(mcParticlesThisMcColl, 0, pdgDB))
        continue;
      switch (filterOnGenPhi) {
        case 0:
          if (!eventHasGenKPair(mcParticlesThisMcColl))
            continue;
          break;
        case 1:
          if (!eventHasGenPhi(mcParticlesThisMcColl))
            continue;
          break;
        default:
          break;
      }

      uint64_t numberAssocColl = 0;
      std::vector<float> zVtxs;

      auto& collIndexesThisMcColl = collsGrouped[mcCollision.globalIndex()];

      for (const auto& collisionIndex : collIndexesThisMcColl) {
        auto collision = collisions.rawIteratorAt(collisionIndex);

        if (acceptEventQA<true>(collision, false)) {
          auto filteredMCTracksThisColl = filteredMCTracks.sliceBy(preslices.trackPerCollision, collision.globalIndex());

          posFiltMCTracks.bindTable(filteredMCTracksThisColl);
          negFiltMCTracks.bindTable(filteredMCTracksThisColl);

          switch (filterOnRecoPhi) {
            case 0:
              if (!eventHasRecoPhi(posFiltMCTracks, negFiltMCTracks))
                continue;
              break;
            case 1:
              if (!eventHasRecoPhiWPDG(posFiltMCTracks, negFiltMCTracks, mcParticles))
                continue;
              break;
            default:
              break;
          }

          mcEventHist.fill(HIST("hRecoMCMultiplicityPercent"), mcCollision.centFT0M());
          mcEventHist.fill(HIST("h2RecoMCVertexZvsMult"), collision.posZ(), mcCollision.centFT0M());

          zVtxs.push_back(collision.posZ());

          for (const auto& track : filteredMCTracksThisColl) {
            if (trackConfigs.applyExtraPhiCuts && ((track.phi() > trackConfigs.extraPhiCuts->at(0) && track.phi() < trackConfigs.extraPhiCuts->at(1)) ||
                                                   track.phi() <= trackConfigs.extraPhiCuts->at(2) || track.phi() >= trackConfigs.extraPhiCuts->at(3)))
              continue;
            if (!track.has_mcParticle())
              continue;

            auto mcTrack = track.mcParticle();
            if (!mcTrack.isPhysicalPrimary() || std::abs(mcTrack.eta()) > trackConfigs.etaMax)
              continue;

            mcEventHist.fill(HIST("h6RecoMCEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcTrack.eta(), mcTrack.phi(), kSpAll, kGlobalplusITSonly);
            if (track.hasTPC()) {
              mcEventHist.fill(HIST("h6RecoMCEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcTrack.eta(), mcTrack.phi(), kSpAll, kGlobalonly);
            } else {
              mcEventHist.fill(HIST("h6RecoMCEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcTrack.eta(), mcTrack.phi(), kSpAll, kITSonly);
            }

            int pid = fromPDGToEnum(mcTrack.pdgCode());
            mcEventHist.fill(HIST("h6RecoMCEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcTrack.eta(), mcTrack.phi(), pid, kGlobalplusITSonly);
          }

          if (fillMcPartsForAllReco) {
            for (const auto& mcParticle : mcParticlesThisMcColl) {
              if (!selectionChargedGenParticle(mcParticle))
                continue;

              mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kNoGenpTVar);
              if (mcParticle.pt() < trackConfigs.cMinChargedParticlePtcut) {
                mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup, -10.0f * mcParticle.pt() + 2.0f);
                mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown, 5.0f * mcParticle.pt() + 0.5f);
              } else {
                mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup);
                mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown);
              }

              int pid = fromPDGToEnum(mcParticle.pdgCode());
              mcEventHist.fill(HIST("h6GenMCAllAssocRecoEtaDistribution"), collision.posZ(), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), pid, kNoGenpTVar);
            }
          }

          numberAssocColl++;
        }
      }

      mcEventHist.fill(HIST("hGenMCMultiplicityPercent"), mcCollision.centFT0M());

      if (numberAssocColl > 0) {
        float zVtxRef = zVtxs[0];
        if (zVtxs.size() > 1) {
          for (size_t i = 1; i < zVtxs.size(); ++i) {
            mcEventHist.fill(HIST("hSplitVertexZ"), zVtxs[i] - zVtxRef);
          }
        }

        mcEventHist.fill(HIST("hGenMCAssocRecoMultiplicityPercent"), mcCollision.centFT0M());
        mcEventHist.fill(HIST("h2GenMCAssocRecoVertexZvsMult"), zVtxRef, mcCollision.centFT0M());
      }

      for (const auto& mcParticle : mcParticlesThisMcColl) {
        if (!selectionChargedGenParticle(mcParticle))
          continue;

        int pid = fromPDGToEnum(mcParticle.pdgCode());

        mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kNoGenpTVar);
        if (mcParticle.pt() < trackConfigs.cMinChargedParticlePtcut) {
          mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup, -10.0f * mcParticle.pt() + 2.0f);
          mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown, 5.0f * mcParticle.pt() + 0.5f);
        } else {
          mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup);
          mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown);
        }
        mcEventHist.fill(HIST("h5GenMCEtaDistribution"), mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), pid, kNoGenpTVar);

        if (numberAssocColl > 0) {
          float zVtxRef = zVtxs[0];

          mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kNoGenpTVar);
          if (mcParticle.pt() < trackConfigs.cMinChargedParticlePtcut) {
            mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup, -10.0f * mcParticle.pt() + 2.0f);
            mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown, 5.0f * mcParticle.pt() + 0.5f);
          } else {
            mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTup);
            mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), kSpAll, kGenpTdown);
          }
          mcEventHist.fill(HIST("h6GenMCAssocRecoEtaDistribution"), zVtxRef, mcCollision.centFT0M(), mcParticle.eta(), mcParticle.phi(), pid, kNoGenpTVar);
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processdNdetaWPhiMC, "Process function for dN/deta values in MC", false);

  // New 2D analysis procedure
  void processPhiK0SPionData2D(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    // Check if the event selection is passed
    if (!acceptEventQA<false>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    dataEventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;
    bool isFilledhV0 = false;

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<false>(track1, true) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      dataPhiHist.fill(HIST("hEta"), track1.eta());
      dataPhiHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      dataPhiHist.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());

      auto track1ID = track1.globalIndex();

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<false>(track2, true) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          dataEventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event
          dataEventHist.fill(HIST("hMultiplicityPercentWithPhi"), multiplicity);
          isCountedPhi = true;
        }

        float efficiencyPhi = 1.0f;
        if (applyEfficiency) {
          efficiencyPhi = effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity());
          if (efficiencyPhi == 0)
            efficiencyPhi = 1.0f;
        }
        float weightPhi = applyEfficiency ? 1.0f / efficiencyPhi : 1.0f;
        dataPhiHist.fill(HIST("h3PhiDataNewProc"), multiplicity, recPhi.Pt(), recPhi.M(), weightPhi);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (!isFilledhV0) {
            dataK0SHist.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
            dataK0SHist.fill(HIST("hV0CosPA"), v0.v0cosPA());

            // Filling the PID of the V0 daughters in the region of the K0 peak
            if (v0Configs.lowMK0S < v0.mK0Short() && v0.mK0Short() < v0Configs.upMK0S) {
              dataK0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi());
              dataK0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
            }
          }

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiK0S = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiK0S = effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapK0S->Interpolate(multiplicity, v0.pt(), v0.yK0Short());
            if (efficiencyPhiK0S == 0)
              efficiencyPhiK0S = 1.0f;
          }
          float weightPhiK0S = applyEfficiency ? 1.0f / efficiencyPhiK0S : 1.0f;
          dataPhiK0SHist.fill(HIST("h5PhiK0SDataNewProc"), v0.yK0Short() - recPhi.Rapidity(), multiplicity, v0.pt(), v0.mK0Short(), recPhi.M(), weightPhiK0S);
        }

        isFilledhV0 = true;

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion<true, false>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiPion = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiPion = track.pt() < trackConfigs.pTToUseTOF ? effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPC->Interpolate(multiplicity, track.pt(), track.rapidity(massPi)) : effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPCTOF->Interpolate(multiplicity, track.pt(), track.rapidity(massPi));
            if (efficiencyPhiPion == 0)
              efficiencyPhiPion = 1.0f;
          }
          float weightPhiPion = applyEfficiency ? 1.0f / efficiencyPhiPion : 1.0f;
          dataPhiPionHist.fill(HIST("h5PhiPiTPCDataNewProc"), track.rapidity(massPi) - recPhi.Rapidity(), multiplicity, track.pt(), track.tpcNSigmaPi(), recPhi.M(), weightPhiPion);
          if (track.hasTOF())
            dataPhiPionHist.fill(HIST("h5PhiPiTOFDataNewProc"), track.rapidity(massPi) - recPhi.Rapidity(), multiplicity, track.pt(), track.tofNSigmaPi(), recPhi.M(), weightPhiPion);
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SPionData2D, "Process function for Phi-K0S and Phi-Pion Correlations in Data2D", false);

  void processPhiK0SPionMCClosure2D(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    if (!collision.has_mcCollision())
      return;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hRecoMCMultiplicityPercent"), genmultiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negMCTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    bool isCountedPhi = false;

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<true>(track1, true) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      auto track1ID = track1.globalIndex();

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<true>(track2, true) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt || recPhi.Pt() > phiConfigs.maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          mcEventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi candidate in the event
          mcEventHist.fill(HIST("hRecoMCMultiplicityPercentWithPhi"), genmultiplicity);
          isCountedPhi = true;
        }

        float efficiencyPhi = 1.0f;
        if (applyEfficiency) {
          efficiencyPhi = effMapPhi->Interpolate(genmultiplicity, recPhi.Pt(), recPhi.Rapidity());
          if (efficiencyPhi == 0)
            efficiencyPhi = 1.0f;
        }
        float weightPhi = applyEfficiency ? 1.0f / efficiencyPhi : 1.0f;
        closureMCPhiHist.fill(HIST("h3PhiMCClosureNewProc"), genmultiplicity, recPhi.Pt(), recPhi.M(), weightPhi);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!v0.has_mcParticle())
              continue;
            auto v0mcparticle = v0.mcParticle();
            if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
              continue;
          }

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiK0S = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiK0S = effMapPhi->Interpolate(genmultiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapK0S->Interpolate(genmultiplicity, v0.pt(), v0.yK0Short());
            if (efficiencyPhiK0S == 0)
              efficiencyPhiK0S = 1.0f;
          }
          float weightPhiK0S = applyEfficiency ? 1.0f / efficiencyPhiK0S : 1.0f;
          closureMCPhiK0SHist.fill(HIST("h5PhiK0SMCClosureNewProc"), v0.yK0Short() - recPhi.Rapidity(), genmultiplicity, v0.pt(), v0.mK0Short(), recPhi.M(), weightPhiK0S);
        }

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!track.has_mcParticle())
              continue;
            auto mcTrack = track.mcParticle_as<aod::McParticles>();
            if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus || !mcTrack.isPhysicalPrimary())
              continue;
          }

          if (!selectionPion<true, true>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiPion = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiPion = track.pt() < trackConfigs.pTToUseTOF ? effMapPhi->Interpolate(genmultiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPC->Interpolate(genmultiplicity, track.pt(), track.rapidity(massPi)) : effMapPhi->Interpolate(genmultiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPCTOF->Interpolate(genmultiplicity, track.pt(), track.rapidity(massPi));
            if (efficiencyPhiPion == 0)
              efficiencyPhiPion = 1.0f;
          }
          float weightPhiPion = applyEfficiency ? 1.0f / efficiencyPhiPion : 1.0f;
          closureMCPhiPionHist.fill(HIST("h5PhiPiTPCMCClosureNewProc"), track.rapidity(massPi) - recPhi.Rapidity(), genmultiplicity, track.pt(), track.tpcNSigmaPi(), recPhi.M(), weightPhiPion);
          if (track.hasTOF())
            closureMCPhiPionHist.fill(HIST("h5PhiPiTOFMCClosureNewProc"), track.rapidity(massPi) - recPhi.Rapidity(), genmultiplicity, track.pt(), track.tofNSigmaPi(), recPhi.M(), weightPhiPion);
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SPionMCClosure2D, "Process function for Phi-K0S and Phi-Pion Correlations in MCClosure2D", false);

  void processPhiK0SMixingEvent2D(SelCollisions const& collisions, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    auto tracksV0sTuple = std::make_tuple(fullTracks, V0s);
    Pair<SelCollisions, FullTracks, FullV0s, BinningTypeVertexCent> pairPhiK0S{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tracksV0sTuple, &cache};

    for (auto const& [collision1, tracks1, collision2, v0s2] : pairPhiK0S) {
      float multiplicity = collision1.centFT0M();

      Partition<FullTracks> posMixTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
      posMixTracks.bindTable(tracks1);
      Partition<FullTracks> negMixTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;
      negMixTracks.bindTable(tracks1);

      for (const auto& [posTrack1, negTrack1, v0] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posMixTracks, negMixTracks, v0s2))) {
        if (!selectionTrackResonance<false>(posTrack1, true) || !selectionPIDKaonpTdependent(posTrack1))
          continue;
        if (!selectionTrackResonance<false>(negTrack1, true) || !selectionPIDKaonpTdependent(negTrack1))
          continue;
        if (posTrack1.globalIndex() == negTrack1.globalIndex())
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(posTrack1, negTrack1, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

        if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
          continue;
        if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision2))
          continue;
        if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
          continue;

        float efficiencyPhiK0S = 1.0f;
        if (applyEfficiency) {
          efficiencyPhiK0S = effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapK0S->Interpolate(multiplicity, v0.pt(), v0.yK0Short());
          if (efficiencyPhiK0S == 0)
            efficiencyPhiK0S = 1.0f;
        }
        float weightPhiK0S = applyEfficiency ? 1.0f / efficiencyPhiK0S : 1.0f;
        mePhiK0SHist.fill(HIST("h5PhiK0SMENewProc"), v0.yK0Short() - recPhi.Rapidity(), multiplicity, v0.pt(), v0.mK0Short(), recPhi.M(), weightPhiK0S);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SMixingEvent2D, "Process Mixed Event for Phi-K0S Analysis 2D", false);

  void processPhiPionMixingEvent2D(SelCollisions const& collisions, FullTracks const& fullTracks)
  {
    auto tracksTuple = std::make_tuple(fullTracks);
    SameKindPair<SelCollisions, FullTracks, BinningTypeVertexCent> pairPhiPion{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto const& [collision1, tracks1, collision2, tracks2] : pairPhiPion) {
      float multiplicity = collision1.centFT0M();

      Partition<FullTracks> posMixTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
      posMixTracks.bindTable(tracks1);
      Partition<FullTracks> negMixTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;
      negMixTracks.bindTable(tracks1);

      for (const auto& [posTrack1, negTrack1, track] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posMixTracks, negMixTracks, tracks2))) {
        if (!selectionTrackResonance<false>(posTrack1, true) || !selectionPIDKaonpTdependent(posTrack1))
          continue;
        if (!selectionTrackResonance<false>(negTrack1, true) || !selectionPIDKaonpTdependent(negTrack1))
          continue;
        if (posTrack1.globalIndex() == negTrack1.globalIndex())
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(posTrack1, negTrack1, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!selectionPion<true, false>(track, false))
          continue;
        if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
          continue;

        float efficiencyPhiPion = 1.0f;
        if (applyEfficiency) {
          efficiencyPhiPion = track.pt() < trackConfigs.pTToUseTOF ? effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPC->Interpolate(multiplicity, track.pt(), track.rapidity(massPi)) : effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPCTOF->Interpolate(multiplicity, track.pt(), track.rapidity(massPi));
          if (efficiencyPhiPion == 0)
            efficiencyPhiPion = 1.0f;
        }
        float weightPhiPion = applyEfficiency ? 1.0f / efficiencyPhiPion : 1.0f;
        mePhiPionHist.fill(HIST("h5PhiPiTPCMENewProc"), track.rapidity(massPi) - recPhi.Rapidity(), multiplicity, track.pt(), track.tpcNSigmaPi(), recPhi.M(), weightPhiPion);
        if (track.hasTOF())
          mePhiPionHist.fill(HIST("h5PhiPiTOFMENewProc"), track.rapidity(massPi) - recPhi.Rapidity(), multiplicity, track.pt(), track.tofNSigmaPi(), recPhi.M(), weightPhiPion);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiPionMixingEvent2D, "Process Mixed Event for Phi-Pion Analysis 2D", false);

  void processAllPartMC(MCCollisions const& mcCollisions, SimCollisions const& collisions, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, aod::McParticles const& mcParticles)
  {

    std::vector<std::vector<int>> collsGrouped(mcCollisions.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision())
        continue;
      const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
      collsGrouped[mcCollision.globalIndex()].push_back(collision.globalIndex());
    }

    for (const auto& mcCollision : mcCollisions) {
      auto mcParticlesThisMcColl = mcParticles.sliceBy(preslices.mcPartPerMCCollision, mcCollision.globalIndex());

      if (!pwglf::isINELgtNmc(mcParticlesThisMcColl, 0, pdgDB))
        continue;
      switch (filterOnGenPhi) {
        case 0:
          if (!eventHasGenKPair(mcParticlesThisMcColl))
            continue;
          break;
        case 1:
          if (!eventHasGenPhi(mcParticlesThisMcColl))
            continue;
          break;
        default:
          break;
      }

      uint64_t numberAssocColl = 0;
      std::vector<float> zVtxs;

      auto& collIndexesThisMcColl = collsGrouped[mcCollision.globalIndex()];

      for (const auto& collisionIndex : collIndexesThisMcColl) {
        auto collision = collisions.rawIteratorAt(collisionIndex);

        if (acceptEventQA<true>(collision, false)) {
          auto fullMCTracksThisColl = fullMCTracks.sliceBy(preslices.trackPerCollision, collision.globalIndex());
          auto v0sThisColl = V0s.sliceBy(preslices.v0PerCollision, collision.globalIndex());

          posMCTracks.bindTable(fullMCTracksThisColl);
          negMCTracks.bindTable(fullMCTracksThisColl);

          switch (filterOnRecoPhi) {
            case 0:
              if (!eventHasRecoPhi(posMCTracks, negMCTracks))
                continue;
              break;
            case 1:
              if (!eventHasRecoPhiWPDG(posMCTracks, negMCTracks, mcParticles))
                continue;
              break;
            default:
              break;
          }

          mcEventHist.fill(HIST("hRecoMCMultiplicityPercent"), mcCollision.centFT0M());
          mcEventHist.fill(HIST("h2RecoMCVertexZvsMult"), collision.posZ(), mcCollision.centFT0M());

          zVtxs.push_back(collision.posZ());

          if ((filterOnGenPhi != 0 && filterOnGenPhi != 1) && (filterOnRecoPhi != 0 && filterOnRecoPhi != 1)) {
            for (const auto& track1 : posMCTracks) { // loop over all selected tracks
              if (!selectionTrackResonance<true>(track1, false) || !selectionPIDKaonpTdependent(track1))
                continue; // topological and PID selection

              auto track1ID = track1.globalIndex();

              if (!track1.has_mcParticle())
                continue;
              auto mcTrack1 = mcParticles.rawIteratorAt(track1.mcParticleId());
              if (mcTrack1.pdgCode() != PDG_t::kKPlus || !mcTrack1.isPhysicalPrimary())
                continue;

              for (const auto& track2 : negMCTracks) {
                if (!selectionTrackResonance<true>(track2, false) || !selectionPIDKaonpTdependent(track2))
                  continue; // topological and PID selection

                auto track2ID = track2.globalIndex();
                if (track2ID == track1ID)
                  continue; // condition to avoid double counting of pair

                if (!track2.has_mcParticle())
                  continue;
                auto mcTrack2 = mcParticles.rawIteratorAt(track2.mcParticleId());
                if (mcTrack2.pdgCode() != PDG_t::kKMinus || !mcTrack2.isPhysicalPrimary())
                  continue;

                const auto mcTrack1MotherIndexes = mcTrack1.mothersIds();
                const auto mcTrack2MotherIndexes = mcTrack2.mothersIds();

                float pTMother = -1.0f;
                float yMother = -1.0f;
                bool isMCMotherPhi = false;

                for (const auto& mcTrack1MotherIndex : mcTrack1MotherIndexes) {
                  for (const auto& mcTrack2MotherIndex : mcTrack2MotherIndexes) {
                    if (mcTrack1MotherIndex != mcTrack2MotherIndex)
                      continue;

                    const auto mother = mcParticles.rawIteratorAt(mcTrack1MotherIndex);
                    if (mother.pdgCode() != o2::constants::physics::Pdg::kPhi)
                      continue;

                    pTMother = mother.pt();
                    yMother = mother.y();
                    isMCMotherPhi = true;
                  }
                }

                if (!isMCMotherPhi)
                  continue;
                if (pTMother < phiConfigs.minPhiPt || std::abs(yMother) > deltaYConfigs.cfgYAcceptance)
                  continue;

                mcPhiHist.fill(HIST("h4PhiMCRecoNewProc"), collision.posZ(), mcCollision.centFT0M(), pTMother, yMother);
              }
            }
          }

          for (const auto& v0 : v0sThisColl) {
            if (!v0.has_mcParticle())
              continue;

            auto v0mcparticle = mcParticles.rawIteratorAt(v0.mcParticleId());
            if (v0mcparticle.pdgCode() != PDG_t::kK0Short || !v0mcparticle.isPhysicalPrimary())
              continue;

            const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
            const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

            if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
              continue;
            if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
              continue;
            if (std::abs(v0mcparticle.y()) > deltaYConfigs.cfgYAcceptance)
              continue;

            mcK0SHist.fill(HIST("h4K0SMCRecoNewProc"), collision.posZ(), mcCollision.centFT0M(), v0mcparticle.pt(), v0mcparticle.y());
          }

          for (const auto& track : fullMCTracksThisColl) {
            // Pion selection
            if (!selectionPion<false, true>(track, false))
              continue;

            if (!track.has_mcParticle())
              continue;

            auto mcTrack = mcParticles.rawIteratorAt(track.mcParticleId());
            if (std::abs(mcTrack.pdgCode()) != PDG_t::kPiPlus)
              continue;

            if (std::abs(mcTrack.y()) > deltaYConfigs.cfgYAcceptance)
              continue;

            // Primary pion selection
            if (mcTrack.isPhysicalPrimary()) {
              mcPionHist.fill(HIST("h3RecMCDCAxyPrimPi"), track.pt(), track.dcaXY());
            } else {
              if (mcTrack.getProcess() == 4) { // Selection of secondary pions from weak decay
                mcPionHist.fill(HIST("h3RecMCDCAxySecWeakDecayPi"), track.pt(), track.dcaXY());
              } else { // Selection of secondary pions from material interactions
                mcPionHist.fill(HIST("h3RecMCDCAxySecMaterialPi"), track.pt(), track.dcaXY());
              }
              continue;
            }

            mcPionHist.fill(HIST("h4PiMCRecoNewProc"), collision.posZ(), mcCollision.centFT0M(), mcTrack.pt(), mcTrack.y());

            if (track.pt() >= trackConfigs.pTToUseTOF && !track.hasTOF())
              continue;

            mcPionHist.fill(HIST("h4PiMCReco2NewProc"), collision.posZ(), mcCollision.centFT0M(), mcTrack.pt(), mcTrack.y());
          }

          numberAssocColl++;
        }
      }

      mcEventHist.fill(HIST("hGenMCMultiplicityPercent"), mcCollision.centFT0M());

      if (numberAssocColl > 0) {
        float zVtxRef = zVtxs[0];
        if (zVtxs.size() > 1) {
          for (size_t i = 1; i < zVtxs.size(); ++i) {
            mcEventHist.fill(HIST("hSplitVertexZ"), zVtxs[i] - zVtxRef);
          }
        }

        mcEventHist.fill(HIST("hGenMCAssocRecoMultiplicityPercent"), mcCollision.centFT0M());
        mcEventHist.fill(HIST("h2GenMCAssocRecoVertexZvsMult"), zVtxRef, mcCollision.centFT0M());
      }

      for (const auto& mcParticle : mcParticlesThisMcColl) {
        if (std::abs(mcParticle.y()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (filterOnGenPhi != 0 && filterOnGenPhi != 1) {
          // Phi selection
          if (mcParticle.pdgCode() == o2::constants::physics::Pdg::kPhi && mcParticle.pt() >= phiConfigs.minPhiPt) {
            mcPhiHist.fill(HIST("h3PhiMCGenNewProc"), mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
            if (numberAssocColl > 0) {
              float zVtxRef = zVtxs[0];
              mcPhiHist.fill(HIST("h4PhiMCGenAssocRecoNewProc"), zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
            }
          }
        }

        // K0S selection
        if (mcParticle.pdgCode() == PDG_t::kK0Short && mcParticle.isPhysicalPrimary() && mcParticle.pt() >= v0Configs.v0SettingMinPt) {
          mcK0SHist.fill(HIST("h3K0SMCGenNewProc"), mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          if (numberAssocColl > 0) {
            float zVtxRef = zVtxs[0];
            mcK0SHist.fill(HIST("h4K0SMCGenAssocRecoNewProc"), zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          }
        }

        // Pion selection
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus && mcParticle.isPhysicalPrimary() && mcParticle.pt() >= trackConfigs.cMinPionPtcut) {
          mcPionHist.fill(HIST("h3PiMCGenNewProc"), mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          if (numberAssocColl > 0) {
            float zVtxRef = zVtxs[0];
            mcPionHist.fill(HIST("h4PiMCGenAssocRecoNewProc"), zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          }
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processAllPartMC, "Process function for all particles (not for phi if triggered on it) in MC", false);

  // New 2D analysis procedure
  void processPhiK0SPionDeltayDeltaphiData2D(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
  {
    // Check if the event selection is passed
    if (!acceptEventQA<false>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    dataEventHist.fill(HIST("hMultiplicityPercent"), multiplicity);

    // Defining positive and negative tracks for phi reconstruction
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (!eventHasRecoPhi(posThisColl, negThisColl))
      return;

    dataEventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event

    bool isCountedPhi = false;
    bool isFilledhV0 = false;

    // Loop over all positive tracks
    for (const auto& track1 : posThisColl) {
      if (!selectionTrackResonance<false>(track1, true) || !selectionPIDKaonpTdependent(track1))
        continue; // topological and PID selection

      dataPhiHist.fill(HIST("hEta"), track1.eta());
      dataPhiHist.fill(HIST("hNsigmaKaonTPC"), track1.tpcInnerParam(), track1.tpcNSigmaKa());
      dataPhiHist.fill(HIST("hNsigmaKaonTOF"), track1.tpcInnerParam(), track1.tofNSigmaKa());

      auto track1ID = track1.globalIndex();

      // Loop over all negative tracks
      for (const auto& track2 : negThisColl) {
        if (!selectionTrackResonance<false>(track2, true) || !selectionPIDKaonpTdependent(track2))
          continue; // topological and PID selection

        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID)
          continue; // condition to avoid double counting of pair

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < phiConfigs.minPhiPt)
          continue;
        if (recPhi.M() < phiConfigs.lowMPhi || recPhi.M() > phiConfigs.upMPhi)
          continue;
        if (std::abs(recPhi.Rapidity()) > deltaYConfigs.cfgYAcceptance)
          continue;

        if (!isCountedPhi)
          isCountedPhi = true;

        float efficiencyPhi = 1.0f;
        if (applyEfficiency) {
          efficiencyPhi = effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity());
          if (efficiencyPhi == 0)
            efficiencyPhi = 1.0f;
        }
        float weightPhi = applyEfficiency ? 1.0f / efficiencyPhi : 1.0f;
        dataPhiHist.fill(HIST("h3PhiDataNewProc"), multiplicity, recPhi.Pt(), recPhi.M(), weightPhi);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          const auto& posDaughterTrack = v0.posTrack_as<V0DauTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (!isFilledhV0) {
            dataK0SHist.fill(HIST("hDCAV0Daughters"), v0.dcaV0daughters());
            dataK0SHist.fill(HIST("hV0CosPA"), v0.v0cosPA());

            // Filling the PID of the V0 daughters in the region of the K0 peak
            if (v0Configs.lowMK0S < v0.mK0Short() && v0.mK0Short() < v0Configs.upMK0S) {
              dataK0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi());
              dataK0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
            }
          }

          if (std::abs(v0.yK0Short()) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiK0S = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiK0S = effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapK0S->Interpolate(multiplicity, v0.pt(), v0.yK0Short());
            if (efficiencyPhiK0S == 0)
              efficiencyPhiK0S = 1.0f;
          }
          float weightPhiK0S = applyEfficiency ? 1.0f / efficiencyPhiK0S : 1.0f;
          dataPhiK0SHist.fill(HIST("h5PhiK0SData2PartCorr"), multiplicity, recPhi.Pt(), v0.pt(), recPhi.Rapidity() - v0.yK0Short(), recPhi.Phi() - v0.phi(), weightPhiK0S);
        }

        isFilledhV0 = true;

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion<true, false>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > deltaYConfigs.cfgYAcceptance)
            continue;

          float efficiencyPhiPion = 1.0f;
          if (applyEfficiency) {
            efficiencyPhiPion = track.pt() < trackConfigs.pTToUseTOF ? effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPC->Interpolate(multiplicity, track.pt(), track.rapidity(massPi)) : effMapPhi->Interpolate(multiplicity, recPhi.Pt(), recPhi.Rapidity()) * effMapPionTPCTOF->Interpolate(multiplicity, track.pt(), track.rapidity(massPi));
            if (efficiencyPhiPion == 0)
              efficiencyPhiPion = 1.0f;
          }
          float weightPhiPion = applyEfficiency ? 1.0f / efficiencyPhiPion : 1.0f;
          dataPhiPionHist.fill(HIST("h5PhiPiData2PartCorr"), multiplicity, recPhi.Pt(), track.pt(), recPhi.Rapidity() - track.rapidity(massPi), recPhi.Phi() - track.phi(), weightPhiPion);
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processPhiK0SPionDeltayDeltaphiData2D, "Process function for Phi-K0S and Phi-Pion Deltay and Deltaphi 2D Correlations in Data", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Phik0shortanalysis>(cfgc)};
}
