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

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TRandom.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TList.h>
#include <TF1.h>
#include <TPDGCode.h>
#include <Math/Vector4D.h>

#include <cstdlib>
#include <cmath>
#include <array>
#include <vector>
#include <algorithm>
#include <string>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Framework/ASoAHelpers.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

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

  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables for track selection (not necessarily common for trigger and the two associated particles)
  struct : ConfigurableGroup {
    Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0f, "Cut on charge"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
    Configurable<float> cMinKaonPtcut{"cMinKaonPtcut", 0.15f, "Track minimum pt cut"};
    Configurable<float> etaMax{"etaMax", 0.8f, "eta max"};
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
    Configurable<float> cMinPionPtcut{"cMinPionPtcut", 0.3f, "Track minimum pt cut"};
    Configurable<int> minTPCnClsFound{"minTPCnClsFound", 70, "min number of found TPC clusters"};
    Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80, "min number of TPC crossed rows"};
    Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
    Configurable<int> minITSnCls{"minITSnCls", 4, "min number of ITS clusters"};
    Configurable<float> maxChi2ITS{"maxChi2ITS", 36.0f, "max chi2 per cluster ITS"};
  } trackConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};
  Configurable<float> minPhiPt{"minPhiPt", 0.4f, "Minimum pT for Phi"};
  Configurable<float> maxPhiPt{"maxPhiPt", 10.0f, "Maximum pT for Phi"};

  // Configurables on phi mass
  Configurable<int> nBinsMPhi{"nBinsMPhi", 13, "N bins in cfgmassPhiaxis"};
  Configurable<float> lowMPhi{"lowMPhi", 1.0095f, "Upper limits on Phi mass for signal extraction"};
  Configurable<float> upMPhi{"upMPhi", 1.029f, "Upper limits on Phi mass for signal extraction"};

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
  } v0Configs;

  // Configurables on K0S mass
  Configurable<float> lowMK0S{"lowMK0S", 0.48f, "Lower limit on K0Short mass"};
  Configurable<float> upMK0S{"upMK0S", 0.52f, "Upper limit on K0Short mass"};

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for pions"};

  // Configurables for delta y selection
  Configurable<int> nBinsY{"nBinsY", 80, "Number of bins in y axis"};
  Configurable<int> nBinsDeltaY{"nBinsDeltaY", 24, "Number of bins in deltay axis"};
  Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
  Configurable<float> cfgYAcceptanceSmear{"cfgYAcceptanceSmear", 0.8f, "Rapidity acceptance for smearing matrix study"};
  Configurable<float> cfgFCutOnDeltaY{"cfgFCutOnDeltaY", 0.5f, "First upper bound on Deltay selection"};
  Configurable<float> cfgSCutOnDeltaY{"cfgSCutOnDeltaY", 0.1f, "Second upper bound on Deltay selection"};
  Configurable<std::vector<float>> cfgDeltaYAcceptanceBins{"cfgDeltaYAcceptanceBins", {1.0f, 0.8f, 0.6f, 0.5f, 0.3f, 0.2f}, "Rapidity acceptance bins"};

  // Configurable for RecMC
  Configurable<bool> cfgiskNoITSROFrameBorder{"cfgiskNoITSROFrameBorder", false, "kNoITSROFrameBorder request on RecMC collisions"};

  // Configurables for MC closure
  Configurable<bool> cfgisRecMCWPDGForClosure1{"cfgisRecMCWPDGForClosure1", false, "RecoMC with PDG Codes for Closure only for Associated particles"};
  Configurable<bool> cfgisRecMCWPDGForClosure2{"cfgisRecMCWPDGForClosure2", false, "RecoMC with PDG Codes for Closure"};
  Configurable<bool> cfgisGenMCForClosure{"cfgisGenMCForClosure", false, "GenMC for Closure"};

  // Configurables to choose the filling method
  Configurable<bool> fillMethodMultipleWeights{"fillMethodMultipleWeights", true, "Fill method Multiple Weights"};
  Configurable<bool> fillMethodSingleWeight{"fillMethodSingleWeight", false, "Fill method Single Weight"};

  // Configurable for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbPurityPath{"ccdbPurityPath", "Users/s/scannito/PhiPuritiesData", "Correction path to file"};

  // Constants
  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiPlus;
  double massK0S = o2::constants::physics::MassK0Short;
  double massLambda = o2::constants::physics::MassLambda0;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZVertex);

  // Defining filters on V0s (cannot filter on dynamic columns)
  Filter preFilterV0 = (nabs(aod::v0data::dcapostopv) > v0Configs.v0SettingDCAPosToPV && nabs(aod::v0data::dcanegtopv) > v0Configs.v0SettingDCANegToPV && aod::v0data::dcaV0daughters < v0Configs.v0SettingDCAV0Dau);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults>;
  using SimCollisions = soa::Join<SelCollisions, aod::McCollisionLabels>;
  using MCCollisions = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;

  // Defining the type of the V0s
  using FullV0s = soa::Filtered<aod::V0Datas>;
  using FullMCV0s = soa::Join<FullV0s, aod::McV0Labels>;

  // Defining the type of the tracks for data and MC
  using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTOFFullPi, aod::pidTOFFullKa>;
  using FullMCTracks = soa::Join<FullTracks, aod::McTrackLabels>;

  using V0DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCFullPi>;
  using V0DauMCTracks = soa::Join<V0DauTracks, aod::McTrackLabels>;

  // Defining the binning policy for mixed event
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  SliceCache cache;

  Partition<FullTracks> posTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullTracks> negTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  Partition<FullMCTracks> posMCTracks = aod::track::signed1Pt > trackConfigs.cfgCutCharge;
  Partition<FullMCTracks> negMCTracks = aod::track::signed1Pt < trackConfigs.cfgCutCharge;

  // Necessary to flag INEL>0 events in GenMC
  Service<o2::framework::O2DatabasePDG> pdgDB;

  // Necessary to get the CCDB for phi purities
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Set of functions for phi purity
  std::vector<std::vector<TF1*>> phiPurityFunctions = std::vector<std::vector<TF1*>>(binsMult->size(), std::vector<TF1*>(binspTPhi->size(), nullptr));

  void init(InitContext&)
  {
    // Axes
    AxisSpec massK0SAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec sigmassPhiAxis = {nBinsMPhi, lowMPhi, upMPhi, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec vertexZAxis = {100, -15.f, 15.f, "vrtx_{Z} [cm]"};
    AxisSpec yAxis = {nBinsY, -cfgYAcceptanceSmear, cfgYAcceptanceSmear, "#it{y}"};
    AxisSpec deltayAxis = {nBinsDeltaY, -1.2f, 1.2f, "#Delta#it{y}"};
    AxisSpec multAxis = {120, 0.0f, 120.0f, "centFT0M"};
    AxisSpec binnedmultAxis{(std::vector<double>)binsMult, "centFT0M"};
    AxisSpec binnedpTPhiAxis{(std::vector<double>)binspTPhi, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptK0SAxis = {60, 0.0f, 6.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedptK0SAxis{(std::vector<double>)binspTK0S, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptPiAxis = {30, 0.0f, 3.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedptPiAxis{(std::vector<double>)binspTPi, "#it{p}_{T} (GeV/#it{c})"};

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
    mcEventHist.add("hRecMCVertexZ", "hRecMCVertexZ", kTH1F, {vertexZAxis});
    mcEventHist.add("hRecMCMultiplicityPercent", "RecMC Multiplicity Percentile", kTH1F, {multAxis});
    mcEventHist.add("hRecMCGenMultiplicityPercent", "RecMC Gen Multiplicity Percentile", kTH1F, {binnedmultAxis});
    mcEventHist.add("hRecMCGenMultiplicityPercentWithPhi", "RecMC Gen Multiplicity Percentile in Events with a Phi Candidate", kTH1F, {binnedmultAxis});

    mcEventHist.add("hGenMCVertexZ", "hGenMCVertexZ", kTH1F, {vertexZAxis});
    mcEventHist.add("hGenMCMultiplicityPercent", "GenMC Multiplicity Percentile", kTH1F, {binnedmultAxis});

    // Phi topological/PID cuts
    dataPhiHist.add("h2DauTracksPhiDCAxyPreCutData", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAzPreCutData", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAxyPostCutData", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    dataPhiHist.add("h2DauTracksPhiDCAzPostCutData", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    dataPhiHist.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    dataPhiHist.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    dataPhiHist.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});

    // Phi invariant mass for computing purities and normalisation
    dataPhiHist.add("h3PhipurData", "Invariant mass of Phi for Purity (no K0S/Pi) in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

    dataPhiHist.add("h4PhipurK0SData", "Invariant mass of Phi for Purity (K0S) in Data", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    dataPhiHist.get<THnSparse>(HIST("h4PhipurK0SData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (int i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      dataPhiHist.get<THnSparse>(HIST("h4PhipurK0SData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    dataPhiHist.add("h4PhipurPiData", "Invariant mass of Phi for Purity (Pi) in Data", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    dataPhiHist.get<THnSparse>(HIST("h4PhipurPiData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (int i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      dataPhiHist.get<THnSparse>(HIST("h4PhipurPiData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // DCA plots for phi daughters in MCReco
    mcPhiHist.add("h2DauTracksPhiDCAxyPreCutMCReco", "Dcaxy distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAzPreCutMCReco", "Dcaz distribution vs pt before DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAxyPostCutMCReco", "Dcaxy distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPhiHist.add("h2DauTracksPhiDCAzPostCutMCReco", "Dcaz distribution vs pt after DCAxy cut", kTH2F, {{100, 0.0, 5.0, "#it{p}_{T} (GeV/#it{c})"}, {2000, -0.05, 0.05, "DCA_{z} (cm)"}});

    // MCPhi invariant mass for computing purities
    closureMCPhiHist.add("h3PhipurMCClosure", "Invariant mass of Phi for Purity (no K0S/Pi)", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis});

    closureMCPhiHist.add("h4PhipurK0SMCClosure", "Invariant mass of Phi for Purity (K0S) in MCClosure", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    closureMCPhiHist.get<THnSparse>(HIST("h4PhipurK0SMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (int i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      closureMCPhiHist.get<THnSparse>(HIST("h4PhipurK0SMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    closureMCPhiHist.add("h4PhipurPiMCClosure", "Invariant mass of Phi for Purity (Pi) in MCClosure", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedpTPhiAxis, massPhiAxis});
    closureMCPhiHist.get<THnSparse>(HIST("h4PhipurPiMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (int i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      closureMCPhiHist.get<THnSparse>(HIST("h4PhipurPiMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // K0S topological/PID cuts
    dataK0SHist.add("hDCAV0Daughters", "hDCAV0Daughters", kTH1F, {{55, 0.0f, 2.2f}});
    dataK0SHist.add("hV0CosPA", "hV0CosPA", kTH1F, {{100, 0.95f, 1.f}});
    dataK0SHist.add("hNSigmaPosPionFromK0S", "hNSigmaPosPionFromK0Short", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});
    dataK0SHist.add("hNSigmaNegPionFromK0S", "hNSigmaNegPionFromK0Short", kTH2F, {{100, 0.0, 5.0, "#it{p} (GeV/#it{c})"}, {100, -10.0f, 10.0f}});

    // 2D mass of Phi and K0S for Data
    dataPhiK0SHist.add("h5PhiK0SData", "2D Invariant mass of Phi and K0Short for Data", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptK0SAxis, massK0SAxis, sigmassPhiAxis});
    dataPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      dataPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // 1D mass of K0S for Data
    dataPhiK0SHist.add("h3PhiK0SSEIncNew", "Invariant mass of K0Short for Same Event Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});
    dataPhiK0SHist.add("h3PhiK0SSEFCutNew", "Invariant mass of K0Short for Same Event Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});
    dataPhiK0SHist.add("h3PhiK0SSESCutNew", "Invariant mass of K0Short for Same Event Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});

    // K0S rapidity in Data
    dataK0SHist.add("h3K0SRapidityData", "K0Short rapidity for Data", kTH3F, {binnedmultAxis, binnedptK0SAxis, yAxis});

    // RecMC K0S coupled to Phi
    mcPhiK0SHist.add("h4PhiK0SMCReco", "K0S coupled to Phi in MCReco", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptK0SAxis, massK0SAxis});
    mcPhiK0SHist.get<THnSparse>(HIST("h4PhiK0SMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiK0SHist.get<THnSparse>(HIST("h4PhiK0SMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // GenMC K0S coupled to Phi
    mcPhiK0SHist.add("h3PhiK0SMCGen", "K0S coupled toPhi in MCGen", kTH3F, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptK0SAxis});
    mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGen"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGen"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    mcPhiK0SHist.add("h3PhiK0SMCGenAssocReco", "K0S coupled toPhi in MCGen Associated MCReco Collision", kTH3F, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptK0SAxis});
    mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGenAssocReco"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiK0SHist.get<TH3>(HIST("h3PhiK0SMCGenAssocReco"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // 2D mass of Phi and K0S for Closure Test
    closureMCPhiK0SHist.add("h5PhiK0SMCClosure", "2D Invariant mass of Phi and K0Short for MC Closure Test", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptK0SAxis, massK0SAxis, sigmassPhiAxis});
    closureMCPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      closureMCPhiK0SHist.get<THnSparse>(HIST("h5PhiK0SMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // 1D mass of K0S for Closure Test
    closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSEIncNew", "Invariant mass of K0Short for Inclusive for Closure Test", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});
    closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSEFCutNew", "Invariant mass of K0Short for Deltay < FirstCut for Closure Test", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});
    closureMCPhiK0SHist.add("h3ClosureMCPhiK0SSESCutNew", "Invariant mass of K0Short for Deltay < SecondCut for Closure Test", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});

    // Phi mass vs Pion NSigma dE/dx for Data
    dataPhiPionHist.add("h6PhiPiData", "Phi Invariant mass vs Pion nSigma TPC/TOF for Data", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigmassPhiAxis});
    dataPhiPionHist.get<THnSparse>(HIST("h6PhiPiData"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      dataPhiPionHist.get<THnSparse>(HIST("h6PhiPiData"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // Pion NSigma dE/dx for Data
    dataPhiPionHist.add("h4PhiPiSEIncNew", "Pion nSigma TPC/TOF for Same Event Inclusive", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    dataPhiPionHist.add("h4PhiPiSEFCutNew", "Pion nSigma TPC/TOF for Same Event Deltay < FirstCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    dataPhiPionHist.add("h4PhiPiSESCutNew", "Pion nSigma TPC/TOF for Same Event Deltay < SecondCut", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});

    // Pion rapidity in Data
    dataPionHist.add("h3PiRapidityData", "Pion rapidity for Data", kTH3F, {binnedmultAxis, binnedptPiAxis, yAxis});

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
    mcPionHist.add("h3RecMCDCAxyPrimPi", "Dcaxy distribution vs pt for Primary Pions", kTH2F, {binnedptPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h3RecMCDCAxySecWeakDecayPi", "Dcaz distribution vs pt for Secondary Pions from Weak Decay", kTH2F, {binnedptPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});
    mcPionHist.add("h3RecMCDCAxySecMaterialPi", "Dcaxy distribution vs pt for Secondary Pions from Material", kTH2F, {binnedptPiAxis, {2000, -0.05, 0.05, "DCA_{xy} (cm)"}});

    // RecMC Pion coupled to Phi with TPC
    mcPhiPionHist.add("h4PhiPiTPCMCReco", "Pion coupled to Phi in MCReco (TPC)", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}});
    mcPhiPionHist.get<THnSparse>(HIST("h4PhiPiTPCMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiPionHist.get<THnSparse>(HIST("h4PhiPiTPCMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // RecMC Pion coupled to Phi with TPC and TOF
    mcPhiPionHist.add("h5PhiPiTPCTOFMCReco", "Pion coupled to Phi in MCReco (TPC and TOF)", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    mcPhiPionHist.get<THnSparse>(HIST("h5PhiPiTPCTOFMCReco"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiPionHist.get<THnSparse>(HIST("h5PhiPiTPCTOFMCReco"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    mcPhiPionHist.add("h3PhiPiMCGen", "Pion coupled to Phi in MCGen", kTH3F, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis});
    mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGen"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGen"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    mcPhiPionHist.add("h3PhiPiMCGenAssocReco", "Pion coupled to Phi in MCGen Associated Reco Collision", kTH3F, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis});
    mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGenAssocReco"))->GetXaxis()->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      mcPhiPionHist.get<TH3>(HIST("h3PhiPiMCGenAssocReco"))->GetXaxis()->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // Phi mass vs Pion NSigma dE/dx for Closure Test
    closureMCPhiPionHist.add("h6PhiPiMCClosure", "Phi Invariant mass vs Pion nSigma TPC/TOF for MC Closure Test", kTHnSparseF, {{cfgDeltaYAcceptanceBins->size(), -0.5f, cfgDeltaYAcceptanceBins->size() + 1.0f - 0.5f}, binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}, sigmassPhiAxis});
    closureMCPhiPionHist.get<THnSparse>(HIST("h6PhiPiMCClosure"))->GetAxis(0)->SetBinLabel(1, "Inclusive");
    for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
      closureMCPhiPionHist.get<THnSparse>(HIST("h6PhiPiMCClosure"))->GetAxis(0)->SetBinLabel(i + 2, Form("|Delta#it{y}| < %.1f", cfgDeltaYAcceptanceBins->at(i)));
    }

    // Phi mass vs Pion NSigma dE/dx for Closure Test
    closureMCPhiPionHist.add("h4ClosureMCPhiPiSEIncNew", "Pion nSigma TPC/TOF for Inclusive for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    closureMCPhiPionHist.add("h4ClosureMCPhiPiSEFCutNew", "Pion nSigma TPC/TOF for Deltay < FirstCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});
    closureMCPhiPionHist.add("h4ClosureMCPhiPiSESCutNew", "Pion nSigma TPC/TOF for Deltay < SecondCut for Closure Test", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});

    // MCPhi invariant mass for computing efficiencies and MCnormalisation
    mcPhiHist.add("h2PhieffInvMass", "Invariant mass of Phi for Efficiency (no K0S/Pi)", kTH2F, {binnedmultAxis, massPhiAxis});

    mcPhiHist.add("h3PhieffK0SInvMassInc", "Invariant mass of Phi for Efficiency (K0S) Inclusive", kTH3F, {binnedmultAxis, binnedptK0SAxis, massPhiAxis});
    mcPhiHist.add("h3PhieffK0SInvMassFCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, massPhiAxis});
    mcPhiHist.add("h3PhieffK0SInvMassSCut", "Invariant mass of Phi for Efficiency (K0S) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptK0SAxis, massPhiAxis});

    mcPhiHist.add("h3PhieffPiInvMassInc", "Invariant mass of Phi for Efficiency (Pi) Inclusive", kTH3F, {binnedmultAxis, binnedptPiAxis, massPhiAxis});
    mcPhiHist.add("h3PhieffPiInvMassFCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < FirstCut", kTH3F, {binnedmultAxis, binnedptPiAxis, massPhiAxis});
    mcPhiHist.add("h3PhieffPiInvMassSCut", "Invariant mass of Phi for Efficiency (Pi) Deltay < SecondCut", kTH3F, {binnedmultAxis, binnedptPiAxis, massPhiAxis});

    // GenMC Phi and Phi coupled to K0S and Pion
    mcPhiHist.add("h1PhiGenMC", "Phi for GenMC", kTH1F, {binnedmultAxis});
    mcPhiHist.add("h1PhiGenMCAssocReco", "Phi for GenMC Associated Reco Collision", kTH1F, {binnedmultAxis});

    mcPhiHist.add("h2PhieffK0SGenMCInc", "Phi coupled to K0Short for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    mcPhiHist.add("h2PhieffK0SGenMCFCut", "Phi coupled to K0Short for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    mcPhiHist.add("h2PhieffK0SGenMCSCut", "Phi coupled to K0Short for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    mcPhiHist.add("h2PhieffK0SGenMCIncAssocReco", "Phi coupled to K0Short for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    mcPhiHist.add("h2PhieffK0SGenMCFCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    mcPhiHist.add("h2PhieffK0SGenMCSCutAssocReco", "Phi coupled to K0Short for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    mcPhiHist.add("h2PhieffPiGenMCInc", "Phi coupled to Pion for GenMC Inclusive", kTH2F, {binnedmultAxis, binnedptPiAxis});
    mcPhiHist.add("h2PhieffPiGenMCFCut", "Phi coupled to Pion for GenMC Deltay < FirstCut", kTH2F, {binnedmultAxis, binnedptPiAxis});
    mcPhiHist.add("h2PhieffPiGenMCSCut", "Phi coupled to Pion for GenMC Deltay < SecondCut", kTH2F, {binnedmultAxis, binnedptPiAxis});

    mcPhiHist.add("h2PhieffPiGenMCIncAssocReco", "Phi coupled to Pion for GenMC Inclusive Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    mcPhiHist.add("h2PhieffPiGenMCFCutAssocReco", "Phi coupled to Pion for GenMC Deltay < FirstCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});
    mcPhiHist.add("h2PhieffPiGenMCSCutAssocReco", "Phi coupled to Pion for GenMC Deltay < SecondCut Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});

    // Rapidity smearing matrix for Phi
    mcPhiHist.add("h3PhiRapiditySmearing", "Rapidity Smearing Matrix for Phi", kTH3F, {binnedmultAxis, yAxis, yAxis});

    // MCK0S invariant mass and GenMC K0S for computing efficiencies
    mcK0SHist.add("h3K0SMCReco", "K0S for MCReco", kTH3F, {binnedmultAxis, binnedptK0SAxis, massK0SAxis});

    mcK0SHist.add("h2K0SMCGen", "K0S for MCGen", kTH2F, {binnedmultAxis, binnedptK0SAxis});
    mcK0SHist.add("h2K0SMCGenAssocReco", "K0S for MCGen Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptK0SAxis});

    // Rapidity smearing matrix for K0S and rapidity in GenMC
    mcK0SHist.add("h4K0SRapiditySmearing", "Rapidity Smearing Matrix for K0Short", kTHnSparseF, {binnedmultAxis, binnedptK0SAxis, yAxis, yAxis});

    mcK0SHist.add("h3K0SRapidityGenMC", "Rapidity for K0Short for GenMC", kTH3F, {binnedmultAxis, binnedptK0SAxis, yAxis});

    // MCPion invariant mass and GenMC Pion for computing efficiencies
    mcPionHist.add("h3PiTPCMCReco", "Pion for MCReco (TPC)", kTH3F, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}});
    mcPionHist.add("h4PiTPCTOFMCReco", "Pion for MCReco (TPC and TOF)", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, {100, -10.0f, 10.0f}, {100, -10.0f, 10.0f}});

    mcPionHist.add("h2PiMCGen", "Pion for GenMC", kTH2F, {binnedmultAxis, binnedptPiAxis});
    mcPionHist.add("h2PiMCGenAssocReco", "Pion for GenMC Associated Reco Collision", kTH2F, {binnedmultAxis, binnedptPiAxis});

    // Rapidity smearing matrix for Pion and rapidity in GenMC
    mcPionHist.add("h4PiRapiditySmearing", "Rapidity Smearing Matrix for Pion", kTHnSparseF, {binnedmultAxis, binnedptPiAxis, yAxis, yAxis});

    mcPionHist.add("h3PiRapidityGenMC", "Rapidity for Pion for GenMC", kTH3F, {binnedmultAxis, binnedptPiAxis, yAxis});

    // Initialize CCDB only if purity is requested in the task
    if (fillMethodSingleWeight) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      getPhiPurityFunctionsFromCCDB();
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
        mcEventHist.fill(HIST("hRecMCVertexZ"), collision.posZ());
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
  bool selectionTrackResonance(const T& track, bool isQA)
  {
    if (trackConfigs.cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (trackConfigs.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (trackConfigs.cfgPVContributor && !track.isPVContributor())
      return false;

    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;

    if (track.pt() < trackConfigs.cMinKaonPtcut)
      return false;
    if (std::abs(track.eta()) > trackConfigs.etaMax)
      return false;

    if (isQA) {
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
    if (isQA) {
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
    if (track.pt() < 0.5 && std::abs(track.tpcNSigmaKa()) < trackConfigs.nSigmaCutTPCKa)
      return true;
    if (track.pt() >= 0.5 && track.hasTOF() && (std::pow(track.tofNSigmaKa(), 2) + std::pow(track.tpcNSigmaKa(), 2)) < std::pow(trackConfigs.nSigmaCutCombinedKa, 2))
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
  bool selectionPion(const T& track, bool isQA)
  {
    if (!track.hasITS())
      return false;
    if (track.itsNCls() < trackConfigs.minITSnCls)
      return false;
    if (track.itsChi2NCl() > trackConfigs.maxChi2ITS)
      return false;

    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < trackConfigs.minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < trackConfigs.minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > trackConfigs.maxChi2TPC)
      return false;

    if (track.pt() < trackConfigs.cMinPionPtcut)
      return false;
    if (std::abs(track.eta()) > trackConfigs.etaMax)
      return false;

    if constexpr (isTOFChecked) {
      if (track.pt() >= 0.5 && !track.hasTOF())
        return false;
    }

    if (isQA) {
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
    if (isQA) {
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

  // Get phi-meson purity functions from CCDB
  void getPhiPurityFunctionsFromCCDB()
  {
    TList* listPhiPurityFunctions = ccdb->get<TList>(ccdbPurityPath);
    if (!listPhiPurityFunctions)
      LOG(fatal) << "Problem getting TList object with phi purity functions!";

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

  // Fill 2D invariant mass histogram for V0 and Phi
  template <bool isMC, typename T>
  void fillInvMass2D(const T& V0, const std::vector<ROOT::Math::PxPyPzMVector>& listPhi, float multiplicity, const std::vector<float>& weights)
  {
    for (const auto& Phi : listPhi) {
      if constexpr (!isMC) { // same event
        dataPhiK0SHist.fill(HIST("h5PhiK0SData"), 0, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
            continue;
          dataPhiK0SHist.fill(HIST("h5PhiK0SData"), i + 1, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(i + 1));
        }
      } else { // MC event
        closureMCPhiK0SHist.fill(HIST("h5PhiK0SMCClosure"), 0, multiplicity, V0.pt(), V0.mK0Short(), Phi.M(), weights.at(0));
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(V0.yK0Short() - Phi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
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
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
            continue;
          dataPhiPionHist.fill(HIST("h6PhiPiData"), i + 1, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(i + 1));
        }
      } else { // MC event
        closureMCPhiPionHist.fill(HIST("h6PhiPiMCClosure"), 0, multiplicity, Pi.pt(), Pi.tpcNSigmaPi(), nSigmaTOFPi, Phi.M(), weights.at(0));
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(Pi.rapidity(massPi) - Phi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
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

  void processQAPurity(SelCollisions::iterator const& collision, FullTracks const& fullTracks, FullV0s const& V0s, V0DauTracks const&)
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
        if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          dataEventHist.fill(HIST("hEventSelection"), 4); // at least a Phi candidate in the event
          dataEventHist.fill(HIST("hMultiplicityPercentWithPhi"), multiplicity);
          isCountedPhi = true;
        }

        if (fillMethodSingleWeight)
          weight *= (1 - getPhiPurity(multiplicity, recPhi));

        dataPhiHist.fill(HIST("h3PhipurData"), multiplicity, recPhi.Pt(), recPhi.M());

        std::vector<int> countsK0S(cfgDeltaYAcceptanceBins->size() + 1, 0);

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
            if (lowMK0S < v0.mK0Short() && v0.mK0Short() < upMK0S) {
              dataK0SHist.fill(HIST("hNSigmaPosPionFromK0S"), posDaughterTrack.tpcInnerParam(), posDaughterTrack.tpcNSigmaPi());
              dataK0SHist.fill(HIST("hNSigmaNegPionFromK0S"), negDaughterTrack.tpcInnerParam(), negDaughterTrack.tpcNSigmaPi());
            }
          }

          if (std::abs(v0.yK0Short()) > cfgYAcceptance)
            continue;

          countsK0S.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsK0S.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsK0S.at(i) > 0)
            dataPhiHist.fill(HIST("h4PhipurK0SData"), i, multiplicity, recPhi.Pt(), recPhi.M());
        }

        isFilledhV0 = true;

        std::vector<int> countsPi(cfgDeltaYAcceptanceBins->size(), 0);

        // Loop over all primary pion candidates
        for (const auto& track : fullTracks) {
          if (!selectionPion<true, false>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgYAcceptance)
            continue;

          countsPi.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsPi.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsPi.at(i) > 0)
            dataPhiHist.fill(HIST("h4PhipurPiData"), i, multiplicity, recPhi.Pt(), recPhi.M());
        }
      }
    }

    weight = 1 - weight;
    dataEventHist.fill(HIST("hEventSelection"), 5, weight); // at least a Phi in the event
  }

  PROCESS_SWITCH(Phik0shortanalysis, processQAPurity, "Process for QA and Phi Purities", true);

  void processSEPhiK0S(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const&, FullV0s const& V0s, V0DauTracks const&)
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

      if (std::abs(v0.yK0Short()) > cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(cfgDeltaYAcceptanceBins->size() + 1, 1);

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
          if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
            continue;
          if (recPhi.M() < lowMPhi || recPhi.M() > upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(multiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMass2D<false>(v0, listrecPhi, multiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillInvMass<false>(v0, multiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processSEPhiK0S, "Process Same Event for Phi-K0S Analysis", false);

  void processSEPhiPion(soa::Filtered<SelCollisions>::iterator const& collision, FullTracks const& fullTracks)
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

      if (std::abs(track.rapidity(massPi)) > cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(cfgDeltaYAcceptanceBins->size() + 1, 1);

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
          if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
            continue;
          if (recPhi.M() < lowMPhi || recPhi.M() > upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(multiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMassNSigma<false>(track, listrecPhi, multiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillNSigma<false>(track, multiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processSEPhiPion, "Process Same Event for Phi-Pion Analysis", false);

  void processRecMCPhiQA(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    float multiplicity = collision.centFT0M();
    mcEventHist.fill(HIST("hRecMCMultiplicityPercent"), multiplicity);

    if (!collision.has_mcCollision())
      return;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hRecMCGenMultiplicityPercent"), genmultiplicity);

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
      if (mcTrack1.pdgCode() != 321 || !mcTrack1.isPhysicalPrimary())
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
        if (mcTrack2.pdgCode() != -321 || !mcTrack2.isPhysicalPrimary())
          continue;

        bool isMCMotherPhi = false;
        auto mcMotherPhi = mcTrack1.mothers_as<aod::McParticles>()[0];
        for (const auto& MotherOfmcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
          for (const auto& MotherOfmcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
            if (MotherOfmcTrack1 == MotherOfmcTrack2 && MotherOfmcTrack1.pdgCode() == 333) {
              mcMotherPhi = MotherOfmcTrack1;
              isMCMotherPhi = true;
            }
          }
        }

        if (!isMCMotherPhi)
          continue;

        ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
        if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
          continue;

        mcPhiHist.fill(HIST("h3PhiRapiditySmearing"), genmultiplicity, recPhi.Rapidity(), mcMotherPhi.y());

        if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
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
          if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
            continue;

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          // Cut on V0 dynamic columns
          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (std::abs(v0.yK0Short()) > cfgYAcceptance)
            continue;
          if (!isCountedK0S.at(0)) {
            mcPhiHist.fill(HIST("h3PhieffK0SInvMassInc"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(0) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgFCutOnDeltaY)
            continue;
          if (!isCountedK0S.at(1)) {
            mcPhiHist.fill(HIST("h3PhieffK0SInvMassFCut"), genmultiplicity, v0.pt(), recPhi.M());
            isCountedK0S.at(1) = true;
          }
          if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgSCutOnDeltaY)
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
          if (std::abs(mcTrack.pdgCode()) != 211 || !mcTrack.isPhysicalPrimary())
            continue;

          if (!selectionPion<true, true>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgYAcceptance)
            continue;
          if (!isCountedPi.at(0)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassInc"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(0) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgFCutOnDeltaY)
            continue;
          if (!isCountedPi.at(1)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassFCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(1) = true;
          }
          if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgSCutOnDeltaY)
            continue;
          if (!isCountedPi.at(2)) {
            mcPhiHist.fill(HIST("h3PhieffPiInvMassSCut"), genmultiplicity, track.pt(), recPhi.M());
            isCountedPi.at(2) = true;
          }
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCPhiQA, "Process for ReCMCQA and Phi in RecMC", false);

  void processRecMCPhiK0S(SimCollisions const& collisions, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const& mcParticles)
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
        if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
          continue;

        const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
        const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

        if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
          continue;
        if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
          continue;

        mcK0SHist.fill(HIST("h4K0SRapiditySmearing"), genmultiplicity, v0.pt(), v0.yK0Short(), v0mcparticle.y());

        if (std::abs(v0mcparticle.y()) > cfgYAcceptance)
          continue;

        mcK0SHist.fill(HIST("h3K0SMCReco"), genmultiplicity, v0mcparticle.pt(), v0.mK0Short());

        std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != 333)
            continue;
          if (mcParticle.pt() < minPhiPt || mcParticle.pt() > maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > cfgYAcceptance)
            continue;

          counts.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0mcparticle.y() - mcParticle.y()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (counts.at(i) > 0)
            mcPhiK0SHist.fill(HIST("h4PhiK0SMCReco"), i, genmultiplicity, v0mcparticle.pt(), v0.mK0Short());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCPhiK0S, "Process RecMC for Phi-K0S Analysis", false);

  void processRecMCPhiPion(SimCollisions const& collisions, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const& mcParticles)
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
        if (std::abs(mcTrack.pdgCode()) != 211)
          continue;

        if (std::abs(mcTrack.y()) > cfgYAcceptance)
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

        std::vector<int> countsTPC(cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != 333)
            continue;
          if (mcParticle.pt() < minPhiPt || mcParticle.pt() > maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > cfgYAcceptance)
            continue;

          countsTPC.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(mcTrack.y() - mcParticle.y()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsTPC.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsTPC.at(i) > 0)
            mcPhiPionHist.fill(HIST("h4PhiPiTPCMCReco"), i, genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi());
        }

        if (track.pt() >= 0.5 && !track.hasTOF())
          continue;

        mcPionHist.fill(HIST("h4PiTPCTOFMCReco"), genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi(), track.tofNSigmaPi());

        std::vector<int> countsTPCTOF(cfgDeltaYAcceptanceBins->size() + 1, 0);

        for (const auto& mcParticle : mcParticlesThisColl) {
          if (mcParticle.pdgCode() != 333)
            continue;
          if (mcParticle.pt() < minPhiPt || mcParticle.pt() > maxPhiPt)
            continue;
          if (std::abs(mcParticle.y()) > cfgYAcceptance)
            continue;

          countsTPCTOF.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(mcTrack.y() - mcParticle.y()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsTPCTOF.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsTPCTOF.at(i) > 0)
            mcPhiPionHist.fill(HIST("h5PhiPiTPCTOFMCReco"), i, genmultiplicity, mcTrack.pt(), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCPhiPion, "Process RecMC for Phi-Pion Analysis", false);

  void processRecMCClosurePhiQA(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
  {
    if (!acceptEventQA<true>(collision, true))
      return;

    if (!collision.has_mcCollision())
      return;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 6); // with at least a gen collision

    const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
    float genmultiplicity = mcCollision.centFT0M();
    mcEventHist.fill(HIST("hRecMCGenMultiplicityPercent"), genmultiplicity);

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
        if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
          continue;
        if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
          continue;

        if (!isCountedPhi) {
          mcEventHist.fill(HIST("hRecMCEventSelection"), 7); // at least a Phi candidate in the event
          mcEventHist.fill(HIST("hRecMCGenMultiplicityPercentWithPhi"), genmultiplicity);
          isCountedPhi = true;
        }

        if (fillMethodSingleWeight)
          weight *= (1 - getPhiPurity(genmultiplicity, recPhi));

        closureMCPhiHist.fill(HIST("h3PhipurMCClosure"), genmultiplicity, recPhi.Pt(), recPhi.M());

        std::vector<int> countsK0S(cfgDeltaYAcceptanceBins->size() + 1, 0);

        // V0 already reconstructed by the builder
        for (const auto& v0 : V0s) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!v0.has_mcParticle())
              continue;
            auto v0mcparticle = v0.mcParticle();
            if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
              continue;
          }

          const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
          const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

          if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
            continue;
          if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
            continue;

          if (std::abs(v0.yK0Short()) > cfgYAcceptance)
            continue;

          countsK0S.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsK0S.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsK0S.at(i) > 0)
            closureMCPhiHist.fill(HIST("h4PhipurK0SInvMass"), i, genmultiplicity, recPhi.Pt(), recPhi.M());
        }

        std::vector<int> countsPi(cfgDeltaYAcceptanceBins->size() + 1, 0);

        // Loop over all primary pion candidates
        for (const auto& track : fullMCTracks) {
          if (cfgisRecMCWPDGForClosure1) {
            if (!track.has_mcParticle())
              continue;
            auto mcTrack = track.mcParticle_as<aod::McParticles>();
            if (std::abs(mcTrack.pdgCode()) != 211 || !mcTrack.isPhysicalPrimary())
              continue;
          }

          if (!selectionPion<true, true>(track, false))
            continue;

          if (std::abs(track.rapidity(massPi)) > cfgYAcceptance)
            continue;

          countsPi.at(0)++;
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            countsPi.at(i + 1)++;
          }
        }

        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          if (countsPi.at(i) > 0)
            closureMCPhiHist.fill(HIST("h4PhipurPiInvMass"), i, genmultiplicity, recPhi.Pt(), recPhi.M());
        }
      }
    }

    weight = 1 - weight;
    mcEventHist.fill(HIST("hRecMCEventSelection"), 8, weight); // at least a Phi in the event
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCClosurePhiQA, "Process for ReCMCQA and Phi in RecMCClosure", false);

  void processRecMCClosurePhiK0S(SimCollisions::iterator const& collision, FullMCTracks const&, FullMCV0s const& V0s, V0DauMCTracks const&, MCCollisions const&, aod::McParticles const&)
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
        if (v0mcparticle.pdgCode() != 310 || !v0mcparticle.isPhysicalPrimary())
          continue;
      }

      const auto& posDaughterTrack = v0.posTrack_as<V0DauMCTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<V0DauMCTracks>();

      if (!selectionV0(v0, posDaughterTrack, negDaughterTrack))
        continue;
      if (v0Configs.cfgFurtherV0Selection && !furtherSelectionV0(v0, collision))
        continue;

      if (std::abs(v0.yK0Short()) > cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(cfgDeltaYAcceptanceBins->size() + 1, 1);

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
            if (mcTrack1.pdgCode() != 321 || !mcTrack1.isPhysicalPrimary())
              continue;

            if (!track2.has_mcParticle())
              continue;
            auto mcTrack2 = track2.mcParticle_as<aod::McParticles>();
            if (mcTrack2.pdgCode() != -321 || !mcTrack2.isPhysicalPrimary())
              continue;

            bool isMCMotherPhi = false;
            for (const auto& motherOfMcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
              for (const auto& motherOfMcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
                if (motherOfMcTrack1.pdgCode() != motherOfMcTrack2.pdgCode())
                  continue;
                if (motherOfMcTrack1.globalIndex() != motherOfMcTrack2.globalIndex())
                  continue;
                if (motherOfMcTrack1.pdgCode() != 333)
                  continue;
                isMCMotherPhi = true;
              }
            }
            if (!isMCMotherPhi)
              continue;
          }

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
            continue;
          if (recPhi.M() < lowMPhi || recPhi.M() > upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(genmultiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(v0.yK0Short() - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMass2D<true>(v0, listrecPhi, genmultiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillInvMass<true>(v0, genmultiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCClosurePhiK0S, "Process RecMC for MCClosure Phi-K0S Analysis", false);

  void processRecMCClosurePhiPion(SimCollisions::iterator const& collision, FullMCTracks const& fullMCTracks, MCCollisions const&, aod::McParticles const&)
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
        if (std::abs(mcTrack.pdgCode()) != 211 || !mcTrack.isPhysicalPrimary())
          continue;
      }

      // Pion selection
      if (!selectionPion<true, true>(track, true))
        continue;

      if (std::abs(track.rapidity(massPi)) > cfgYAcceptance)
        continue;

      std::vector<ROOT::Math::PxPyPzMVector> listrecPhi;
      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);
      std::vector<float> weights(cfgDeltaYAcceptanceBins->size() + 1, 1);

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
            if (mcTrack1.pdgCode() != 321 || !mcTrack1.isPhysicalPrimary())
              continue;

            if (!track2.has_mcParticle())
              continue;
            auto mcTrack2 = track2.mcParticle_as<aod::McParticles>();
            if (mcTrack2.pdgCode() != -321 || !mcTrack2.isPhysicalPrimary())
              continue;

            bool isMCMotherPhi = false;
            for (const auto& motherOfMcTrack1 : mcTrack1.mothers_as<aod::McParticles>()) {
              for (const auto& motherOfMcTrack2 : mcTrack2.mothers_as<aod::McParticles>()) {
                if (motherOfMcTrack1.pdgCode() != motherOfMcTrack2.pdgCode())
                  continue;
                if (motherOfMcTrack1.globalIndex() != motherOfMcTrack2.globalIndex())
                  continue;
                if (motherOfMcTrack1.pdgCode() != 333)
                  continue;
                isMCMotherPhi = true;
              }
            }
            if (!isMCMotherPhi)
              continue;
          }

          ROOT::Math::PxPyPzMVector recPhi = recMother(track1, track2, massKa, massKa);
          if (recPhi.Pt() < minPhiPt || recPhi.Pt() > maxPhiPt)
            continue;
          if (recPhi.M() < lowMPhi || recPhi.M() > upMPhi)
            continue;
          if (std::abs(recPhi.Rapidity()) > cfgYAcceptance)
            continue;

          double phiPurity{};
          if (fillMethodSingleWeight)
            phiPurity = getPhiPurity(genmultiplicity, recPhi);

          if (fillMethodMultipleWeights)
            listrecPhi.push_back(recPhi);

          counts.at(0)++;
          weights.at(0) *= (1 - phiPurity);
          for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
            if (std::abs(track.rapidity(massPi) - recPhi.Rapidity()) > cfgDeltaYAcceptanceBins->at(i))
              continue;
            counts.at(i + 1)++;
            weights.at(i + 1) *= (1 - phiPurity);
          }
        }
      }

      if (fillMethodMultipleWeights) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1. / static_cast<float>(counts.at(i)) : 0);
        }
        fillInvMassNSigma<true>(track, listrecPhi, genmultiplicity, weights);
      } else if (fillMethodSingleWeight) {
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
          weights.at(i) = (counts.at(i) > 0 ? 1 - weights.at(i) : 0);
        }
        fillNSigma<true>(track, genmultiplicity, weights);
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processRecMCClosurePhiPion, "Process RecMC for MCClosure Phi-Pion Analysis", false);

  void processGenMCPhiQA(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
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
      if (mcParticle1.pdgCode() != 333)
        continue;
      auto kDaughters = mcParticle1.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2)
        continue;
      bool isPosKaon = false, isNegKaon = false;
      for (const auto& kDaughter : kDaughters) {
        if (kDaughter.pdgCode() == 321)
          isPosKaon = true;
        if (kDaughter.pdgCode() == -321)
          isNegKaon = true;
      }
      if (!isPosKaon || !isNegKaon)
        continue;
      if (std::abs(mcParticle1.y()) > cfgYAcceptance)
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
        if (mcParticle2.pdgCode() != 310)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;

        if (std::abs(mcParticle2.y()) > cfgYAcceptance)
          continue;
        if (!isCountedK0S.at(0)) {
          mcPhiHist.fill(HIST("h2PhieffK0SGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFCutOnDeltaY)
          continue;
        if (!isCountedK0S.at(1)) {
          mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffK0SGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedK0S.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSCutOnDeltaY)
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
        if (std::abs(mcParticle2.pdgCode()) != 211)
          continue;
        if (!mcParticle2.isPhysicalPrimary())
          continue;

        if (std::abs(mcParticle2.y()) > cfgYAcceptance)
          continue;
        if (!isCountedPi.at(0)) {
          mcPhiHist.fill(HIST("h2PhieffPiGenMCInc"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffPiGenMCIncAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(0) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgFCutOnDeltaY)
          continue;
        if (!isCountedPi.at(1)) {
          mcPhiHist.fill(HIST("h2PhieffPiGenMCFCut"), genmultiplicity, mcParticle2.pt());
          if (isAssocColl)
            mcPhiHist.fill(HIST("h2PhieffPiGenMCFCutAssocReco"), genmultiplicity, mcParticle2.pt());
          isCountedPi.at(1) = true;
        }
        if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgSCutOnDeltaY)
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

  PROCESS_SWITCH(Phik0shortanalysis, processGenMCPhiQA, "Process for ReCMCQA and Phi in RecMC", false);

  void processGenMCPhiK0S(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
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
      if (mcParticle1.pdgCode() != 310)
        continue;
      if (!mcParticle1.isPhysicalPrimary() || mcParticle1.pt() < v0Configs.v0SettingMinPt)
        continue;

      mcK0SHist.fill(HIST("h3K0SRapidityGenMC"), genmultiplicity, mcParticle1.pt(), mcParticle1.y());

      if (std::abs(mcParticle1.y()) > cfgYAcceptance)
        continue;

      mcK0SHist.fill(HIST("h2K0SMCGen"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        mcK0SHist.fill(HIST("h2K0SMCGenAssocReco"), genmultiplicity, mcParticle1.pt());

      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        if (cfgisGenMCForClosure) {
          auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
          if (kDaughters.size() != 2)
            continue;
          bool isPosKaon = false, isNegKaon = false;
          for (const auto& kDaughter : kDaughters) {
            if (kDaughter.pdgCode() == 321)
              isPosKaon = true;
            if (kDaughter.pdgCode() == -321)
              isNegKaon = true;
          }
          if (!isPosKaon || !isNegKaon)
            continue;
        }
        if (mcParticle2.pt() < minPhiPt || mcParticle2.pt() > maxPhiPt)
          continue;
        if (std::abs(mcParticle2.y()) > cfgYAcceptance)
          continue;

        counts.at(0)++;
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgDeltaYAcceptanceBins->at(i))
            continue;
          counts.at(i + 1)++;
        }
      }

      for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
        if (counts.at(i) > 0) {
          mcPhiK0SHist.fill(HIST("h3PhiK0SMCGen"), i, genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            mcPhiK0SHist.fill(HIST("h3PhiK0SMCGenAssocReco"), i, genmultiplicity, mcParticle1.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processGenMCPhiK0S, "Process GenMC for Phi-K0S Analysis", false);

  void processGenMCPhiPion(MCCollisions::iterator const& mcCollision, soa::SmallGroups<SimCollisions> const& collisions, aod::McParticles const& mcParticles)
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
      if (std::abs(mcParticle1.pdgCode()) != 211)
        continue;
      if (!mcParticle1.isPhysicalPrimary() || mcParticle1.pt() < trackConfigs.cMinPionPtcut)
        continue;

      mcPionHist.fill(HIST("h3PiRapidityGenMC"), genmultiplicity, mcParticle1.pt(), mcParticle1.y());

      if (std::abs(mcParticle1.y()) > cfgYAcceptance)
        continue;

      mcPionHist.fill(HIST("h2PiMCGen"), genmultiplicity, mcParticle1.pt());
      if (isAssocColl)
        mcPionHist.fill(HIST("h2PiMCGenAssocReco"), genmultiplicity, mcParticle1.pt());

      std::vector<int> counts(cfgDeltaYAcceptanceBins->size() + 1, 0);

      for (const auto& mcParticle2 : mcParticles) {
        if (mcParticle2.pdgCode() != 333)
          continue;
        if (cfgisGenMCForClosure) {
          auto kDaughters = mcParticle2.daughters_as<aod::McParticles>();
          if (kDaughters.size() != 2)
            continue;
          bool isPosKaon = false, isNegKaon = false;
          for (const auto& kDaughter : kDaughters) {
            if (kDaughter.pdgCode() == 321)
              isPosKaon = true;
            if (kDaughter.pdgCode() == -321)
              isNegKaon = true;
          }
          if (!isPosKaon || !isNegKaon)
            continue;
        }
        if (mcParticle2.pt() < minPhiPt || mcParticle2.pt() > maxPhiPt)
          continue;
        if (std::abs(mcParticle2.y()) > cfgYAcceptance)
          continue;

        counts.at(0)++;
        for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size(); i++) {
          if (std::abs(mcParticle1.y() - mcParticle2.y()) > cfgDeltaYAcceptanceBins->at(i))
            continue;
          counts.at(i + 1)++;
        }
      }

      for (size_t i = 0; i < cfgDeltaYAcceptanceBins->size() + 1; i++) {
        if (counts.at(i) > 0) {
          mcPhiPionHist.fill(HIST("h3PhiPiMCGen"), i, genmultiplicity, mcParticle1.pt());
          if (isAssocColl)
            mcPhiPionHist.fill(HIST("h3PhiPiMCGenAssocReco"), i, genmultiplicity, mcParticle1.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(Phik0shortanalysis, processGenMCPhiPion, "Process GenMC for Phi-Pion Analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Phik0shortanalysis>(cfgc)};
}
