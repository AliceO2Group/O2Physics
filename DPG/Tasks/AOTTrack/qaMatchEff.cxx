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

/// \file qaMatchEff.cxx
/// \brief ITS-TPC track matching checks, per charge and with PID
///
/// \author Rosario Turrisi  <rosario.turrisi@pd.infn.it>, INFN-PD
/// \author Mattia Faggin <mattia.faggin@ts.infn.it>, UniTs & INFN-TS

//
//  Internal version number: 6.3
//
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "CommonConstants/MathConstants.h"
#include "CCDB/BasicCCDBManager.h"
//
#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

//
namespace extConfPar
{
static constexpr int nParDCA = 1;
static constexpr int nParVaDCA = 2;
static const std::vector<std::string> parClassDCA{"TrVtx"};
static const std::vector<std::string> parNameDCA{"dcaXY", "dcaZ"};
static const float parTableDCA[nParDCA][nParVaDCA]{{9999.f, 99999.f}};
static constexpr int nParPID = 2;
static constexpr int nParVaPID = 6;
static const std::vector<std::string> parClassPID{"TPC", "TOF"};
static const std::vector<std::string> parNamePID{"nSigPionMin", "nSigPionMax", "nSigKaonMin", "nSigKaonMax", "nSigProtonMin", "nSigProtonMax"};
static const float parTablePID[nParPID][nParVaPID]{
  {-99999.f, 999999.f, -99999.f, 999999.f, -99999.f, 999999.f},
  {-99999.f, 999999.f, -99999.f, 999999.f, -99999.f, 999999.f}};
} // namespace extConfPar
//
// base namespaces
using namespace o2;
using namespace o2::constants::math;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace extConfPar;
//
//
struct qaMatchEff {
  int lastRunNumber = -1;
  bool timeMonitorSetUp = false;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  using BCsWithTimeStamp = soa::Join<aod::BCs, aod::Timestamps>;

  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> enableMonitorVsTime{"enableMonitorVsTime", false, "Enable the storage of ITS-TPC matching efficiency vs. time"};
  Configurable<bool> enableTHnSparseMonitorVsTime{"enableTHnSparseMonitorVsTime", false, "Enable the storage of ITS-TPC matching efficiency vs. time"};
  //
  // histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  //
  // Track selections
  Configurable<bool> b_useTrackSelections{"b_useTrackSelections", false, "Boolean to switch the track selections on/off."};
  Configurable<int> filterbitTrackSelections{"filterbitTrackSelections", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  // kinematics
  Configurable<float> ptMinCutInnerWallTPC{"ptMinCutInnerWallTPC", 0.1f, "Minimum transverse momentum calculated at the inner wall of TPC (GeV/c)"};
  Configurable<float> ptMinCut{"ptMinCut", 0.1f, "Minimum transverse momentum (GeV/c)"};
  Configurable<float> ptMaxCut{"ptMaxCut", 100.f, "Maximum transverse momentum (GeV/c)"};
  Configurable<float> etaMinCut{"etaMinCut", -2.0f, "Minimum pseudorapidity"};
  Configurable<float> etaMaxCut{"etaMaxCut", 2.0f, "Maximum pseudorapidity"};
  Configurable<bool> b_useTPCinnerWallPt{"b_useTPCinnerWallPt", false, "Boolean to switch the usage of pt calculated at the inner wall of TPC on/off."};
  // DCA and PID cuts
  Configurable<LabeledArray<float>> dcaMaxCut{"dcaMaxCut", {parTableDCA[0], nParDCA, nParVaDCA, parClassDCA, parNameDCA}, "Track DCA cuts"};
  Configurable<LabeledArray<float>> nSigmaPID{"nSigmaPID", {parTablePID[0], nParPID, nParVaPID, parClassPID, parNamePID}, "PID nSigma cuts TPC and TOF"};
  // TPC
  Configurable<int> tpcNClusterMin{"tpcNClusterMin", 0, "Minimum number of clusters in TPC"};
  Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 70, "Minimum number of crossed rows in TPC"};
  Configurable<float> tpcNCrossedRowsOverFindableClstMin{"tpcNCrossedRowsOverFindableClstMin", 0.8f, "Minimum fracion of crossed rows over findable custers in TPC"};
  Configurable<float> tpcChi2Max{"tpcChi2Max", 4.0f, "Maximum chi2 in TPC"};
  // ITS
  Configurable<float> itsChi2Max{"itsChi2Max", 36.0f, "Maximum chi2 in ITS"};
  Configurable<int> customITShitmap{"customITShitmap", 3, "ITS hitmap (think to the binary representation)"};
  Configurable<int> customMinITShits{"customMinITShits", 1, "Minimum number of layers crossed by a track among those in \"customITShitmap\""};
  // Other track settings
  //  TRD presence
  Configurable<int> isTRDThere{"isTRDThere", 2, "Integer to turn the presence of TRD off, on, don't care (0,1,anything else)"};
  //
  Configurable<bool> isitMC{"isitMC", false, "Reading MC files, data if false"};
  Configurable<bool> doDebug{"doDebug", false, "Flag of debug information"};
  // Histogram configuration
  //
  // histos axes limits
  Configurable<float> etaMin{"eta-min", -2.0f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 2.0f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.0f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 1.0f * TwoPI, "Upper limit in phi"};
  // histos bins
  Configurable<int> etaBins{"eta-bins", 40, "Number of eta bins"};
  Configurable<int> phiBins{"phi-bins", 18, "Number of phi bins"};
  Configurable<int> qoptBins{"qopt-bins", 500, "Number of Q/pt bins"};
  //
  // special histo, few particles explicitly stored, then pdg>3000
  Configurable<int> pdgBins{"pdg-bins", 14, "Number of pdg values counted"};
  //
  // histo axes
  //
  ConfigurableAxis ptBins{"ptBins", {100, 0.f, 20.f}, "pT binning"};
  ConfigurableAxis ptBinsVsTime{"ptBinsVsTime", {VARIABLE_WIDTH, 0.1, 0.5, 1.0, 2.0, 5.0}, "pT binning for monitorning vs time"};
  ConfigurableAxis ptInvserseBinsVsTime{"ptInverseBinsVsTime", {VARIABLE_WIDTH, 0, 0.2, 0.5, 1, 2, 10}, "1/pT binning for monitorning vs time"};
  ConfigurableAxis etaBinsVsTime{"etaBinsVsTime", {14, -1.4, 1.4}, "eta binning for monitoring vs time"};
  ConfigurableAxis posZBinsVsTime{"posZBinsVsTime", {2, -100, 100}, "posZ primary vertex binning for monitoring vs time"};
  ConfigurableAxis tpcClstBinsVsTime{"tpcClstBinsVsTime", {40, 0, 160}, "TPC cluster binning for monitoring vs time"};
  ConfigurableAxis itsClstBinsVsTime{"itsClstBinsVsTime", {9, 0, 9}, "ITS cluster binning for monitoring vs time"};
  //
  AxisSpec axisPDG{pdgBins, 0, pdgBins + 1.000, "pdgclass"};
  //
  AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
  //
  AxisSpec axisQoPt{qoptBins, -20, 20, "#Q/it{p}_{T} (GeV/#it{c})^{-1}"};
  //
  AxisSpec axisEta{etaBins, etaMin, etaMax, "#eta"};
  AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};
  AxisSpec axisDEta{etaBins, etaMin, etaMax, "D#eta"};
  AxisSpec axisDPh{phiBins, -PI, PI, "D#it{#varphi} (rad)"};
  //
  // pdg codes vector
  std::vector<int> pdgChoice = {211, 213, 215, 217, 219, 221, 223, 321, 411, 521, 2212, 1114, 2214};
  //
  // configuration for THnSparse's
  //
  Configurable<bool> makethn{"makethn", false, "choose if produce thnsparse"};
  ConfigurableAxis thnd0{"thnd0", {600, -3.0f, 3.0f}, "impact parameter in xy [cm]"};
  ConfigurableAxis thnPt{"thnPt", {30, 0.0f, 15.0f}, "pt [GeV/c]"};
  ConfigurableAxis thnPhi{"thnPhi", {18, 0.0f, TMath::TwoPi()}, "phi"};
  ConfigurableAxis thnEta{"thnEta", {20, -2.0f, 2.0f}, "eta"};
  ConfigurableAxis thnType{"thnType", {3, -0.5f, 2.5f}, "0: primary, 1: physical secondary, 2: sec. from material"};
  ConfigurableAxis thnLabelSign{"thnLabelSign", {3, -1.5f, 1.5f}, "-1/+1 antip./particle"};
  ConfigurableAxis thnSpec{"thnSpec", {5, 0.5f, 5.5f}, "particle from MC (1,2,3,4,5 -> e,pi,K,P,other)"};
  AxisSpec thnd0Axis{thnd0, "#it{d}_{r#it{#varphi}} [cm]"};
  AxisSpec thnPtAxis{thnPt, "#it{p}_{T}^{reco} [GeV/#it{c}]"};
  AxisSpec thnPhiAxis{thnPhi, "#varphi"};
  AxisSpec thnEtaAxis{thnEta, "#it{#eta}"};
  AxisSpec thnTypeAxis{thnType, "0:prim-1:sec-2:matsec"};
  AxisSpec thnLabelSignAxis{thnLabelSign, "+/- 1 for part./antipart."};
  AxisSpec thnSpecAxis{thnSpec, "particle from MC (1,2,3,4,5 -> e,pi,K,P,other)"};
  //
  // Tracks selection object
  TrackSelection cutObject;
  //
  // do you want pt comparison 2d's ?
  Configurable<bool> makept2d{"makept2d", false, "choose if produce pt reco/TPC derived pt 2dims "};
  //
  // common flags for PID
  Configurable<bool> isPIDPionRequired{"isPIDPionRequired", false, "choose if apply pion PID"};
  Configurable<bool> isPIDKaonRequired{"isPIDKaonRequired", false, "choose if apply kaon PID"};
  Configurable<bool> isPIDProtonRequired{"isPIDProtonRequired", false, "choose if apply proton PID"};
  //
  // limit for z position of primary vertex
  Configurable<float> zPrimVtxMax{"zPrimVtxax", 999.f, "Maximum asbolute value of z of primary vertex"};
  //
  //
  //      ******     BE VERY CAREFUL!   --  FILTERS !!!  *****
  //
  Filter zPrimVtxLim = nabs(aod::collision::posZ) < zPrimVtxMax;
  Filter trackFilter = (filterbitTrackSelections.node() == 0) ||
                       ((filterbitTrackSelections.node() == 1) && requireGlobalTrackInFilter()) || /// filterbit 4 track selections + tight DCA cuts
                       ((filterbitTrackSelections.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((filterbitTrackSelections.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((filterbitTrackSelections.node() == 4) && requireQualityTracksInFilter()) ||
                       ((filterbitTrackSelections.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  //
  //
  //
  //  Init function
  //
  void init(o2::framework::InitContext&)
  {
    if (doDebug)
      LOG(info) << "===========================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  is it MC? = " << isitMC;
    //
    // let's know if it's MC or data
    if (isitMC)
      initMC();
    else
      initData();

    if ((!isitMC && (doprocessMCFilteredTracks || doprocessMC || doprocessMCNoColl || doprocessTrkIUMC)) || (isitMC && (doprocessDataFilteredTracks || doprocessData || doprocessDataNoColl || doprocessTrkIUMC)))
      LOGF(fatal, "Initialization set for MC and processData function flagged  (or viceversa)! Fix the configuration.");
    if ((doprocessMCFilteredTracks && doprocessMC && doprocessMCNoColl && doprocessTrkIUMC) || (doprocessDataFilteredTracks && doprocessData && doprocessDataNoColl && doprocessTrkIUData))
      LOGF(fatal, "Cannot process for both without collision tag and with collision tag at the same time! Fix the configuration.");
    if (doprocessTrkIUMC && makethn) {
      LOGF(fatal, "No DCA for IU tracks. Put makethn = false.");
    }
    if (doprocessTrkIUData && makethn) {
      LOGF(fatal, "No DCA for IU tracks. Put makethn = false.");
    }
    //
    /// initialize the track selections
    if (b_useTrackSelections) {
      // kinematics
      cutObject.SetEtaRange(etaMinCut, etaMaxCut);
      cutObject.SetPtRange(ptMinCut, ptMaxCut);
      cutObject.SetMaxDcaXY(dcaMaxCut->get("TrVtx", "dcaXY")); /// max for dca implementend by hand in isTrackSelectedKineCuts
      cutObject.SetMaxDcaZ(dcaMaxCut->get("TrVtx", "dcaZ"));   /// max for dca implementend by hand in isTrackSelectedKineCuts
      // TPC
      cutObject.SetMinNClustersTPC(tpcNClusterMin);
      cutObject.SetMinNCrossedRowsTPC(tpcNCrossedRowsMin);
      cutObject.SetMinNCrossedRowsOverFindableClustersTPC(tpcNCrossedRowsOverFindableClstMin);
      cutObject.SetMaxChi2PerClusterTPC(tpcChi2Max);
      // ITS
      cutObject.SetMaxChi2PerClusterITS(itsChi2Max);
      // ITS hitmap
      std::set<uint8_t> set_customITShitmap; // = {};
      for (int index_ITSlayer = 0; index_ITSlayer < 7; index_ITSlayer++) {
        if ((customITShitmap & (1 << index_ITSlayer)) > 0) {
          set_customITShitmap.insert(static_cast<uint8_t>(index_ITSlayer));
        }
      }
      LOG(info) << "### customITShitmap: " << customITShitmap;
      LOG(info) << "### customMinITShits: " << customMinITShits;
      LOG(info) << "### set_customITShitmap.size(): " << set_customITShitmap.size();
      LOG(info) << "### Custom ITS hitmap checked: ";
      for (std::set<uint8_t>::iterator it = set_customITShitmap.begin(); it != set_customITShitmap.end(); it++) {
        LOG(info) << "Layer " << static_cast<int>(*it) << " ";
      }
      LOG(info) << "############";
      cutObject.SetRequireHitsInITSLayers(customMinITShits, set_customITShitmap);
    }
  }
  // end Init function
  //
  //
  // Init Data function - define data histograms
  void initData()
  {
    if (doDebug)
      LOGF(info, "*********************************************************** DATA  ***************************************************");
    //
    // data histos
    //
    // thnsparse for fractions - only if selected
    if (makethn)
      histos.add("data/sparse/thnsforfrac", "Sparse histo for imp. par. fraction analysis - data", kTHnSparseF, {thnd0Axis, thnPtAxis, thnPhiAxis, thnEtaAxis, thnTypeAxis, thnLabelSignAxis, thnSpecAxis});

    /// control plots
    histos.add("data/control/itsHitsMatched", "No. of hits vs ITS layer for ITS-TPC matched tracks;layer ITS", kTH2D, {{8, -1.5, 6.5}, {8, -0.5, 7.5, "No. of hits"}});
    histos.add("data/control/zPrimary", "Position of primary vertex along beam axis;z position [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/control/yPrimary", "Position of primary vertex along y axis;y position [cm]", kTH1D, {{200, -0.1, 0.1}}, true);
    histos.add("data/control/xPrimary", "Position of primary vertex along x axis;x position [cm]", kTH1D, {{200, -0.1, 0.1}}, true);
    histos.add("data/control/chi2Prim", "#chi^2 of primary vertex fit;#chi^2", kTH1D, {{200, 0., 100.0}}, true);
    histos.add("data/control/zDCA_tpc", "DCA along z TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/control/xyDCA_tpc", "DCA in x-y plane TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/control/zDCA_tpcits", "DCA along z TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/control/xyDCA_tpcits", "DCA in x-y plane TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    //
    histos.add("data/PID/zDCA_tpc_pi", "DCA along z - pions TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpc_pi", "DCA in x-y plane - pions TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/zDCA_tpcits_pi", "DCA along z - pions TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpcits_pi", "DCA in x-y plane - pions TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    //
    histos.add("data/PID/zDCA_tpc_ka", "DCA along z - kaons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpc_ka", "DCA in x-y plane - kaons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/zDCA_tpcits_ka", "DCA along z - kaons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpcits_ka", "DCA in x-y plane - kaons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    //
    histos.add("data/PID/zDCA_tpc_pr", "DCA along z - protons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpc_pr", "DCA in x-y plane - protons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/zDCA_tpcits_pr", "DCA along z - protons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpcits_pr", "DCA in x-y plane - protons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    //
    histos.add("data/PID/zDCA_tpc_noid", "DCA along z - no pi/K/P TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpc_noid", "DCA in x-y plane - no pi/K/P TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/zDCA_tpcits_noid", "DCA along z - no pi/K/P TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("data/PID/xyDCA_tpcits_noid", "DCA in x-y plane - no pi/K/P TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    /// compare pt's (tracking and innerParamTPC)
    if (makept2d) {
      histos.add("data/control/ptptconfTPCall", "Tracking pt vs TPC inner wall pt - TPC tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("data/control/ptptconfITSall", "Tracking pt vs TPC inner wall pt - ITS tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("data/control/ptptconfTPCITS", "Tracking pt vs TPC inner wall pt - TPC & ITS tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("data/control/ptptconfITSo", "Tracking pt vs TPC inner wall pt - ITS-only tracks", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("data/control/ptptconfTPCo", "Tracking pt vs TPC inner wall pt - TPC-only tracks", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
    }
    //
    // tpc, its and tpc+its request for all, positive and negative charges vs

    // Q/pt
    histos.add("data/qopthist_tpc", "Q/#it{p}_{T} distribution - data TPC tag", kTH1D, {axisQoPt}, true);
    histos.add("data/qopthist_tpcits", "Q/#it{p}_{T} distribution - data TPC+ITS tag", kTH1D, {axisQoPt}, true);

    // pt, phi, eta
    histos.add("data/pthist_tpc", "#it{p}_{T} distribution - data TPC tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpc", "#eta distribution - data TPC tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpc", "#phi distribution - data TPC tag", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_tpcits", "#it{p}_{T} distribution - data TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpcits", "#eta distribution - data TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpcits", "#phi distribution - data TPC+ITS tag", kTH1D, {axisPhi}, true);

    // pt, phi, eta TOF tagged
    histos.add("data/pthist_toftpc", "#it{p}_{T} distribution - data TOF+TPC tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_toftpc", "#eta distribution - data TOF+TPC tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_toftpc", "#phi distribution - data TOF+TPC tag", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_toftpcits", "#it{p}_{T} distribution - data TOF+TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_toftpcits", "#eta distribution - data TOF+TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_toftpcits", "#phi distribution - data TOF+TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    // if you want just pions
    if (isPIDPionRequired) {
      histos.add("data/PID/pthist_tpc_pi", "#it{p}_{T} distribution - data TPC tag - pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_pi", "#eta distribution - data TPC tag - pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_pi", "#phi distribution - data TPC tag - pions", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_pi", "#it{p}_{T} distribution - data TPC+ITS tag - pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_pi", "#eta distribution - data TPC+ITS tag - pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_pi", "#phi distribution - data TPC+ITS tag - pions", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_pi_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - pions PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pi_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pi_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_pi_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pi_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_piplus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - pions+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piplus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_piplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_piminus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - pions- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piminus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - pions- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_piminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_piminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - pions- PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_pi_PIDTPC", "#it{p}_{T} distribution - data TPC tag - pions PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pi_PIDTOF", "#it{p}_{T} distribution - data TPC tag - pions PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pi_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - pions PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_pi_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - pions PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pi_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - pions PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_piplus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - pions+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piplus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - pions+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - pions+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_piplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - pions+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - pions+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_piminus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - pions- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piminus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - pions- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - pions- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_piminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - pions- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_piminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - pions- PID TOF only", kTH1D, {axisPt}, true);

      // plus
      histos.add("data/PID/pthist_tpc_piplus", "#it{p}_{T} distribution - data TPC tag - pos. pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_piplus", "#eta distribution - data TPC tag - pos. pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_piplus", "#phi distribution - data TPC tag - pos. pions", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_piplus", "#it{p}_{T} distribution - data TPC+ITS tag - pos. pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_piplus", "#eta distribution - data TPC+ITS tag - pos. pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_piplus", "#phi distribution - data TPC+ITS tag - pos. pions", kTH1D, {axisPhi}, true);
      // minus
      histos.add("data/PID/pthist_tpc_piminus", "#it{p}_{T} distribution - data TPC tag - neg. pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_piminus", "#eta distribution - data TPC tag - neg. pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_piminus", "#phi distribution - data TPC tag - neg. pions", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_piminus", "#it{p}_{T} distribution - data TPC+ITS tag - neg. pions", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_piminus", "#eta distribution - data TPC+ITS tag - neg. pions", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_piminus", "#phi distribution - data TPC+ITS tag - neg. pions", kTH1D, {axisPhi}, true);
    }
    //
    // if you want just kaons
    if (isPIDKaonRequired) {
      histos.add("data/PID/pthist_tpc_ka", "#it{p}_{T} distribution - data TPC tag - kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_ka", "#eta distribution - data TPC tag - kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_ka", "#phi distribution - data TPC tag - kaons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_ka", "#it{p}_{T} distribution - data TPC+ITS tag - kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_ka", "#eta distribution - data TPC+ITS tag - kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_ka", "#phi distribution - data TPC+ITS tag - kaons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_ka_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - kaons PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_ka_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_ka_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_ka_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_ka_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_ka_PIDTPC", "#it{p}_{T} distribution - data TPC tag - kaons PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_ka_PIDTOF", "#it{p}_{T} distribution - data TPC tag - kaons PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_ka_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - kaons PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_ka_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - kaons PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_ka_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - kaons PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_kaplus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - kaons+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaplus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_kaplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_kaplus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - kaons+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaplus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - kaons+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - kaons+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_kaplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - kaons+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - kaons+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_kaminus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - kaons- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaminus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - kaons- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_kaminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_kaminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - kaons- PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_kaminus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - kaons- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaminus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - kaons- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - kaons- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_kaminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - kaons- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_kaminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - kaons- PID TOF only", kTH1D, {axisPt}, true);

      // plus
      histos.add("data/PID/pthist_tpc_kaplus", "#it{p}_{T} distribution - data TPC tag - pos. kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_kaplus", "#eta distribution - data TPC tag - pos. kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_kaplus", "#phi distribution - data TPC tag - pos. kaons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_kaplus", "#it{p}_{T} distribution - data TPC+ITS tag - pos. kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_kaplus", "#eta distribution - data TPC+ITS tag - pos. kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_kaplus", "#phi distribution - data TPC+ITS tag - pos. kaons", kTH1D, {axisPhi}, true);
      // minus
      histos.add("data/PID/pthist_tpc_kaminus", "#it{p}_{T} distribution - data TPC tag - neg. kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_kaminus", "#eta distribution - data TPC tag - neg. kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_kaminus", "#phi distribution - data TPC tag - neg. kaons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_kaminus", "#it{p}_{T} distribution - data TPC+ITS tag - neg. kaons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_kaminus", "#eta distribution - data TPC+ITS tag - neg. kaons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_kaminus", "#phi distribution - data TPC+ITS tag - neg. kaons", kTH1D, {axisPhi}, true);
    }
    //
    // if you want just protons
    if (isPIDProtonRequired) {
      histos.add("data/PID/pthist_tpc_pr", "#it{p}_{T} distribution - data TPC tag - protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_pr", "#eta distribution - data TPC tag - protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_pr", "#phi distribution - data TPC tag - protons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_pr", "#it{p}_{T} distribution - data TPC+ITS tag - protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_pr", "#eta distribution - data TPC+ITS tag - protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_pr", "#phi distribution - data TPC+ITS tag - protons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_pr_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - protons PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pr_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pr_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_pr_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_pr_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_pr_PIDTPC", "#it{p}_{T} distribution - data TPC tag - protons PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pr_PIDTOF", "#it{p}_{T} distribution - data TPC tag - protons PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pr_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - protons PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_pr_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - protons PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_pr_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - protons PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_prplus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - protons+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prplus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_prplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_prplus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - protons+ PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prplus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - protons+ PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prplus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - protons+ PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_prplus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - protons+ PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prplus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - protons+ PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_prminus_PIDTPC", "#it{p}_{T} distribution - data TPC+ITS tag - protons- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prminus_PIDTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC+ITS tag - protons- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpcits_prminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpcits_prminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC+ITS tag - protons- PID TOF only", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_prminus_PIDTPC", "#it{p}_{T} distribution - data TPC tag - protons- PID w/TPC at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prminus_PIDTOF", "#it{p}_{T} distribution - data TPC tag - protons- PID w/TOF at least", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prminus_PIDTPCTOF", "#it{p}_{T} distribution - data TPC tag - protons- PID TPC+TOF", kTH1D, {axisPt}, true);

      histos.add("data/PID/pthist_tpc_prminus_PIDTPC_O", "#it{p}_{T} distribution - data TPC tag - protons- PID TPC only", kTH1D, {axisPt}, true);
      histos.add("data/PID/pthist_tpc_prminus_PIDTOF_O", "#it{p}_{T} distribution - data TPC tag - protons- PID TOF only", kTH1D, {axisPt}, true);

      // plus
      histos.add("data/PID/pthist_tpc_prplus", "#it{p}_{T} distribution - data TPC tag - pos. protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_prplus", "#eta distribution - data TPC tag - pos. protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_prplus", "#phi distribution - data TPC tag - pos. protons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_prplus", "#it{p}_{T} distribution - data TPC+ITS tag - pos. protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_prplus", "#eta distribution - data TPC+ITS tag - pos. protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_prplus", "#phi distribution - data TPC+ITS tag - pos. protons", kTH1D, {axisPhi}, true);
      // minus
      histos.add("data/PID/pthist_tpc_prminus", "#it{p}_{T} distribution - data TPC tag - neg. protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_prminus", "#eta distribution - data TPC tag - neg. protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_prminus", "#phi distribution - data TPC tag - neg. protons", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_prminus", "#it{p}_{T} distribution - data TPC+ITS tag - neg. protons", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_prminus", "#eta distribution - data TPC+ITS tag - neg. protons", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_prminus", "#phi distribution - data TPC+ITS tag - neg. protons ", kTH1D, {axisPhi}, true);
    }
    //
    // if PID is required, build also non-identified spectra
    if (isPIDPionRequired || isPIDKaonRequired || isPIDProtonRequired) {
      histos.add("data/PID/pthist_tpc_noid", "#it{p}_{T} distribution - data TPC tag - no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_noid", "#eta distribution - data TPC tag - no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_noid", "#phi distribution - data TPC tag - no ident.", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_noid", "#it{p}_{T} distribution - data TPC+ITS tag - no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_noid", "#eta distribution - data TPC+ITS tag - no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_noid", "#phi distribution - data TPC+ITS tag - no ident.", kTH1D, {axisPhi}, true);
      // plus
      histos.add("data/PID/pthist_tpc_noidplus", "#it{p}_{T} distribution - data TPC tag - pos. no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_noidplus", "#eta distribution - data TPC tag - pos. no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_noidplus", "#phi distribution - data TPC tag - pos. no ident.", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_noidplus", "#it{p}_{T} distribution - data TPC+ITS tag - pos. no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_noidplus", "#eta distribution - data TPC+ITS tag - pos. no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_noidplus", "#phi distribution - data TPC+ITS tag - pos. no ident.", kTH1D, {axisPhi}, true);
      // minus
      histos.add("data/PID/pthist_tpc_noidminus", "#it{p}_{T} distribution - data TPC tag - neg. no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpc_noidminus", "#eta distribution - data TPC tag - neg. no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpc_noidminus", "#phi distribution - data TPC tag - neg. no ident.", kTH1D, {axisPhi}, true);

      histos.add("data/PID/pthist_tpcits_noidminus", "#it{p}_{T} distribution - data TPC+ITS tag - neg. no ident.", kTH1D, {axisPt}, true);
      histos.add("data/PID/etahist_tpcits_noidminus", "#eta distribution - data TPC+ITS tag - neg. no ident.", kTH1D, {axisEta}, true);
      histos.add("data/PID/phihist_tpcits_noidminus", "#phi distribution - data TPC+ITS tag - neg. no ident. ", kTH1D, {axisPhi}, true);
    }
    //

    histos.add("data/pthist_tpc_pos", "#it{p}_{T} distribution - data q>0 TPC tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpc_pos", "#eta distribution - data q>0 TPC tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpc_pos", "#phi distribution - data q>0 TPC tag", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_tpcits_pos", "#it{p}_{T} distribution - data q>0 TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpcits_pos", "#eta distribution - data q>0 TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpcits_pos", "#phi distribution - data q>0 TPC+ITS tag", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_tpc_neg", "#it{p}_{T} distribution - data q<0 TPC tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpc_neg", "#eta distribution - data q<0 TPC tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpc_neg", "#phi distribution - data q<0 TPC tag", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_tpcits_neg", "#it{p}_{T} distribution - data q<0 TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpcits_neg", "#eta distribution - data q<0 TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpcits_neg", "#phi distribution - data q<0 TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    // pt>0.5 GeV/c threshold
    histos.add("data/pthist_tpc_05", "#it{p}_{T} distribution - data TPC tag, #it{p}_{T}>0.5", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpc_05", "#eta distribution - data TPC tag, #it{p}_{T}>0.5", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpc_05", "#phi distribution - data TPC tag, #it{p}_{T}>0.5", kTH1D, {axisPhi}, true);

    histos.add("data/pthist_tpcits_05", "#it{p}_{T} distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1D, {axisPt}, true);
    histos.add("data/etahist_tpcits_05", "#eta distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1D, {axisEta}, true);
    histos.add("data/phihist_tpcits_05", "#phi distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1D, {axisPhi}, true);
  }
  //
  // Init MC function
  void initMC()
  {
    if (doDebug)
      LOGF(info, " +++++++++++++++++++++++  MC  ++++++++++++++++++++++++");

    //
    // adding histos to the registry
    // data histos
    // tpc request and tpc+its request for all, positive and negative charges
    // and for phys. primaries, decay secondaries and mat. secondaries (both
    // charges) vs pt, phi, eta (36 histos tot) pions only, also split in prim
    // secd secm
    //
    // thnsparse for fractions
    if (makethn)
      histos.add("MC/sparse/thnsforfrac", "Sparse histo for imp. par. fraction analysis - MC", kTHnSparseF, {thnd0Axis, thnPtAxis, thnPhiAxis, thnEtaAxis, thnTypeAxis, thnLabelSignAxis, thnSpecAxis});

    /// control plots
    histos.add("MC/control/itsHitsMatched", "No. of hits vs ITS layer for ITS-TPC matched tracks;layer ITS", kTH2D, {{8, -1.5, 6.5}, {8, -0.5, 7.5, "No. of hits"}});
    histos.add("MC/control/zPrimary", "Position of primary vertex along beam axis;z position [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/control/yPrimary", "Position of primary vertex along y axis;y position [cm]", kTH1D, {{200, -0.1, 0.1}}, true);
    histos.add("MC/control/xPrimary", "Position of primary vertex along x axis;x position [cm]", kTH1D, {{200, -0.1, 0.1}}, true);
    histos.add("MC/control/chi2Prim", "#chi^2 of primary vertex fit;#chi^2", kTH1D, {{200, 0., 100.0}}, true);
    histos.add("MC/control/zDCA_tpc", "DCA along z TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/control/xyDCA_tpc", "DCA in x-y plane TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/control/zDCA_tpcits", "DCA along z TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/control/xyDCA_tpcits", "DCA in x-y plane TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    /// compare pt's (tracking and innerParamTPC)
    if (makept2d) {
      histos.add("MC/control/ptptconfTPCall", "Tracking pt vs TPC inner wall pt - TPC tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("MC/control/ptptconfITSall", "Tracking pt vs TPC inner wall pt - ITS tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("MC/control/ptptconfTPCITS", "Tracking pt vs TPC inner wall pt - TPC & ITS tag", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("MC/control/ptptconfITSo", "Tracking pt vs TPC inner wall pt - ITS-only tracks", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
      histos.add("MC/control/ptptconfTPCo", "Tracking pt vs TPC inner wall pt - TPC-only tracks", kTH2D, {{100, 0.0, 10.0, "tracking #it{p}_{T}"}, {100, 0.0, 10.0, "TPC #it{p}_{T}"}});
    }
    //
    // all, positive, negative

    // Q/pt
    histos.add("MC/qopthist_tpc", "Q/#it{p}_{T} distribution - MC TPC tag", kTH1D, {axisQoPt}, true);
    histos.add("MC/qopthist_tpcits", "Q/#it{p}_{T} distribution - MC TPC+ITS tag", kTH1D, {axisQoPt}, true);
    //
    //  TOF tag
    histos.add("MC/pthist_toftpc", "#it{p}_{T} distribution - MC TOF+TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_toftpc", "#eta distribution - MC TOF+TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_toftpc", "#phi distribution - MC TOF+TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_toftpcits", "#it{p}_{T} distribution - MC TOF+TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_toftpcits", "#eta distribution - MC TOF+TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_toftpcits", "#phi distribution - MC TOF+TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    //
    histos.add("MC/pthist_tpc", "#it{p}_{T} distribution - MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpc", "#eta distribution - MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpc", "#phi distribution - MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpcits", "#it{p}_{T} distribution - MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpcits", "#eta distribution - MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpcits", "#phi distribution - MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpc_pos", "#it{p}_{T} distribution - MC q>0 TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpc_pos", "#eta distribution - MC q>0 TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpc_pos", "#phi distribution - MC q>0 TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpcits_pos", "#it{p}_{T} distribution - MC q>0 TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pos", "#eta distribution - MC q>0 TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pos", "#phi distribution - MC q>0 TPC+ITS tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpc_neg", "#it{p}_{T} distribution - MC q<0 TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpc_neg", "#eta distribution - MC q<0 TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpc_neg", "#phi distribution - MC q<0 TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpcits_neg", "#it{p}_{T} distribution - MC q<0 TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpcits_neg", "#eta distribution - MC q<0 TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpcits_neg", "#phi distribution - MC q<0 TPC+ITS tag", kTH1D, {axisPhi}, true);

    //
    // primaries, secondaries

    // Q/pt
    histos.add("MC/primsec/qopthist_tpc_prim", "Q/#it{p}_{T} distribution - MC prim TPC tag", kTH1D, {axisQoPt}, true);
    histos.add("MC/primsec/qopthist_tpcits_prim", "Q/#it{p}_{T} distribution - MC prim TPC+ITS tag", kTH1D, {axisQoPt}, true);

    histos.add("MC/primsec/zDCA_tpc_prim", "DCA along z TPC tag - primaries;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/zDCA_tpcits_prim", "DCA along z TPC+ITS tag - primaries;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/xyDCA_tpc_prim", "DCA in x-y plane TPC tag - primaries;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/xyDCA_tpcits_prim", "DCA in x-y plane TPC+ITS tag - primaries;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/pthist_tpc_prim", "#it{p}_{T} distribution - MC prim TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpc_prim", "#eta distribution - MC prim TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpc_prim", "#phi distribution - MC prim TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/primsec/pthist_tpcits_prim", "#it{p}_{T} distribution - MC prim TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpcits_prim", "#eta distribution - MC prim TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpcits_prim", "#phi distribution - MC prim TPC+ITS tag", kTH1D, {axisPhi}, true);

    // Q/pt
    histos.add("MC/primsec/qopthist_tpc_secd", "Q/#it{p}_{T} distribution - MC dec. sec. TPC tag", kTH1D, {axisQoPt}, true);
    histos.add("MC/primsec/qopthist_tpcits_secd", "Q/#it{p}_{T} distribution - MC dec. sec. TPC+ITS tag", kTH1D, {axisQoPt}, true);

    histos.add("MC/primsec/zDCA_tpc_secd", "DCA along z TPC tag - secd;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/zDCA_tpcits_secd", "DCA along z TPC+ITS tag - secd;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/xyDCA_tpc_secd", "DCA in x-y plane TPC tag - secd;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/xyDCA_tpcits_secd", "DCA in x-y plane TPC+ITS tag - secd;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/pthist_tpc_secd", "#it{p}_{T} distribution - MC dec. sec. TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpc_secd", "#eta distribution - MC dec. sec. TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpc_secd", "#phi distribution - MC dec. sec. TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/primsec/pthist_tpcits_secd", "#it{p}_{T} distribution - MC dec.sec. TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpcits_secd", "#eta distribution - MC dec. sec. TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpcits_secd", "#phi distribution - MC dec. sec. TPC+ITS tag", kTH1D, {axisPhi}, true);

    // Q/pt
    histos.add("MC/primsec/qopthist_tpc_secm", "Q/#it{p}_{T} distribution - MC mat. sec. TPC tag", kTH1D, {axisQoPt}, true);
    histos.add("MC/primsec/qopthist_tpcits_secm", "Q/#it{p}_{T} distribution - MC mat. sec. TPC+ITS tag", kTH1D, {axisQoPt}, true);

    histos.add("MC/primsec/zDCA_tpc_secm", "DCA along z TPC tag - secm;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/zDCA_tpcits_secm", "DCA along z TPC+ITS tag - secm;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/xyDCA_tpc_secm", "DCA in x-y plane TPC tag - secm;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/primsec/xyDCA_tpcits_secm", "DCA in x-y plane TPC+ITS tag - secm;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/primsec/pthist_tpc_secm", "#it{p}_{T} distribution - MC mat. sec. TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpc_secm", "#eta distribution - MC mat. sec. TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpc_secm", "#phi distribution - MC mat. sec. TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/primsec/pthist_tpcits_secm", "#it{p}_{T} distribution - MC mat.sec. TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/primsec/etahist_tpcits_secm", "#eta distribution - MC mat. sec. TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/primsec/phihist_tpcits_secm", "#phi distribution - MC mat. sec. TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    // pions only
    // all
    histos.add("MC/PID/zDCA_tpc_pi", "DCA along z - pions TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpc_pi", "DCA in x-y plane - pions TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/zDCA_tpcits_pi", "DCA along z - pions TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpcits_pi", "DCA in x-y plane - pions TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/PID/zDCA_tpc_ka", "DCA along z - kaons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpc_ka", "DCA in x-y plane - kaons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/zDCA_tpcits_ka", "DCA along z - kaons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpcits_ka", "DCA in x-y plane - kaons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/PID/zDCA_tpc_pr", "DCA along z - protons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpc_pr", "DCA in x-y plane - protons TPC tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/zDCA_tpcits_pr", "DCA along z - protons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);
    histos.add("MC/PID/xyDCA_tpcits_pr", "DCA in x-y plane - protons TPC+ITS tag;dca [cm]", kTH1D, {{200, -20.0, 20.0}}, true);

    histos.add("MC/PID/pthist_tpc_pi", "#it{p}_{T} distribution - #pi MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_pi", "#eta distribution - #pi MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_pi", "#phi distribution - #pi MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_pi", "#it{p}_{T} distribution - #pi MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_pi", "#eta distribution - #pi MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_pi", "#phi distribution - #pi MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // plus
    histos.add("MC/PID/pthist_tpc_piplus", "#it{p}_{T} distribution -pos. #pi MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_piplus", "#eta distribution -pos. #pi MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_piplus", "#phi distribution -pos. #pi MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_piplus", "#it{p}_{T} distribution -pos. #pi MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_piplus", "#eta distribution -pos. #pi MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_piplus", "#phi distribution -pos. #pi MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // minus
    histos.add("MC/PID/pthist_tpc_piminus", "#it{p}_{T} distribution -neg. #pi MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_piminus", "#eta distribution -neg. #pi MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_piminus", "#phi distribution -neg. #pi MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_piminus", "#it{p}_{T} distribution -neg. #pi MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_piminus", "#eta distribution -neg. #pi MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_piminus", "#phi distribution -neg. #pi MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // pions only
    // split in prim secd secm
    histos.add("MC/PID/pthist_tpc_pi_prim", "#it{p}_{T} distribution - #pi MC prim TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_pi_prim", "#eta distribution - #pi MC prim TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_pi_prim", "#phi distribution - #pi MC prim TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_pi_prim", "#it{p}_{T} distribution - #pi MC prim TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_pi_prim", "#eta distribution - #pi MC prim TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_pi_prim", "#phi distribution - #pi MC prim TPC+ITS tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpc_pi_secd", "#it{p}_{T} distribution - #pi MC dec. sec. TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_pi_secd", "#eta distribution - #pi MC dec. sec. TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_pi_secd", "#phi distribution - #pi MC dec. sec. TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_pi_secd", "#it{p}_{T} distribution - #pi MC dec.sec. TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_pi_secd", "#eta distribution - #pi MC dec. sec. TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_pi_secd", "#phi distribution - #pi MC dec. sec. TPC+ITS tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpc_pi_secm", "#it{p}_{T} distribution - #pi MC mat. sec. TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_pi_secm", "#eta distribution - #pi MC mat. sec. TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_pi_secm", "#phi distribution - #pi MC mat. sec. TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_pi_secm", "#it{p}_{T} distribution - #pi MC mat.sec. TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_pi_secm", "#eta distribution - #pi MC mat. sec. TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_pi_secm", "#phi distribution - #pi MC mat. sec. TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    // protons only
    // all
    histos.add("MC/PID/pthist_tpc_pr", "#it{p}_{T} distribution - prot MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_pr", "#eta distribution - prot MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_pr", "#phi distribution - prot MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_pr", "#it{p}_{T} distribution - prot MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_pr", "#eta distribution - prot MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_pr", "#phi distribution - prot MC TPC+ITS tag", kTH1D, {axisPhi}, true);
    // plus
    histos.add("MC/PID/pthist_tpc_prplus", "#it{p}_{T} distribution - pos. prot MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_prplus", "#eta distribution - pos. prot MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_prplus", "#phi distribution - pos. prot MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_prplus", "#it{p}_{T} distribution - pos. prot MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_prplus", "#eta distribution - pos. prot MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_prplus", "#phi distribution - pos. prot MC TPC+ITS tag", kTH1D, {axisPhi}, true);
    // minus
    histos.add("MC/PID/pthist_tpc_prminus", "#it{p}_{T} distribution - neg. prot MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_prminus", "#eta distribution - neg. prot MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_prminus", "#phi distribution - neg. prot MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_prminus", "#it{p}_{T} distribution - neg. prot MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_prminus", "#eta distribution - neg. prot MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_prminus", "#phi distribution - neg. prot MC TPC+ITS tag", kTH1D, {axisPhi}, true);
    //
    // kaons only
    // all
    histos.add("MC/PID/pthist_tpc_ka", "#it{p}_{T} distribution - kaons MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_ka", "#eta distribution - kaons MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_ka", "#phi distribution - kaons MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_ka", "#it{p}_{T} distribution - kaons MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_ka", "#eta distribution - kaons MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_ka", "#phi distribution - kaons MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // plus
    histos.add("MC/PID/pthist_tpc_kaplus", "#it{p}_{T} distribution - pos. kaons MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_kaplus", "#eta distribution - pos. kaons MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_kaplus", "#phi distribution - pos. kaons MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_kaplus", "#it{p}_{T} distribution - pos. kaons MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_kaplus", "#eta distribution - pos. kaons MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_kaplus", "#phi distribution - pos. kaons MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // minus
    histos.add("MC/PID/pthist_tpc_kaminus", "#it{p}_{T} distribution - neg. kaons MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_kaminus", "#eta distribution - neg. kaons MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_kaminus", "#phi distribution - neg. kaons MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_kaminus", "#it{p}_{T} distribution - neg. kaons MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_kaminus", "#eta distribution - neg. kaons MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_kaminus", "#phi distribution - neg. kaons MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // pions+kaons
    // all
    histos.add("MC/PID/pthist_tpc_piK", "#it{p}_{T} distribution - #pi+kaons MC TPC tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_piK", "#eta distribution - #pi+kaons MC TPC tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_piK", "#phi distribution - #pi+kaons MC TPC tag", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_piK", "#it{p}_{T} distribution - #pi+kaons MC TPC+ITS tag", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_piK", "#eta distribution - #pi+kaons MC TPC+ITS tag", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_piK", "#phi distribution - #pi+kaons MC TPC+ITS tag", kTH1D, {axisPhi}, true);

    // pt>0.5 GeV/c threshold
    histos.add("MC/pthist_tpc_05", "#it{p}_{T} distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpc_05", "#eta distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpc_05", "#phi distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1D, {axisPhi}, true);

    histos.add("MC/pthist_tpcits_05", "#it{p}_{T} distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1D, {axisPt}, true);
    histos.add("MC/etahist_tpcits_05", "#eta distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1D, {axisEta}, true);
    histos.add("MC/phihist_tpcits_05", "#phi distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1D, {axisPhi}, true);

    //
    // all but primary/secondary pions
    histos.add("MC/PID/pthist_tpc_nopi", "#it{p}_{T} distribution - MC TPC tag ! prim/secd #pi", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpc_nopi", "#eta distribution - MC TPC tag ! prim/secd #pi", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpc_nopi", "#phi distribution - MC TPC tag ! prim/secd #pi", kTH1D, {axisPhi}, true);

    histos.add("MC/PID/pthist_tpcits_nopi", "#it{p}_{T} distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1D, {axisPt}, true);
    histos.add("MC/PID/etahist_tpcits_nopi", "#eta distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1D, {axisEta}, true);
    histos.add("MC/PID/phihist_tpcits_nopi", "#phi distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1D, {axisPhi}, true);
    //
    // extras: difference between reconstructed and MC truth for eta, phi
    histos.add("MC/control/etahist_diff", "#eta difference track-MC ", kTH1D, {axisDEta}, true);
    histos.add("MC/control/phihist_diff", "#phi difference track-MC", kTH1D, {axisDPh}, true);
    //
    // hist sorting out PDG codes in wide bins
    histos.add("MC/PID/pdghist_num", "PDG code - when non primary #pi TPC+ITS tag", kTH1D, {axisPDG}, true);
    histos.add("MC/PID/pdghist_den", "PDG code - when non primary #pi TPC tag", kTH1D, {axisPDG}, true);
    histos.add("MC/PID/pdghist_denits", "PDG code - when non primary #pi ITS tag", kTH1D, {axisPDG}, true);

  } // end initMC

  /// Function calculating the pt at inner wall of TPC
  template <typename T>
  float computePtInParamTPC(T& track)
  {
    /// Using pt calculated at the inner wall of TPC
    /// Caveat: tgl still from tracking: this is not the value of tgl at the
    /// inner wall of TPC
    return track.tpcInnerParam() / sqrt(1.f + track.tgl() * track.tgl());
  }

  /// Function applying the kinematic selections
  template <typename T>
  bool isTrackSelectedKineCuts(T& track)
  {
    if (!b_useTrackSelections)
      return true; // no track selections applied
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kPtRange))
      return false;
    if (b_useTPCinnerWallPt && computePtInParamTPC(track) < ptMinCutInnerWallTPC) {
      return false; // pt selection active only if the required pt is that calculated at the inner wall of TPC
    }
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kEtaRange))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kDCAxy))
      return false;
    // dcaZ selection to simulate the dca cut in QC ()
    if (abs(track.dcaZ()) > dcaMaxCut->get("TrVtx", "dcaZ"))
      return false;
    return true;
  }
  /// Function applying the TPC selections
  template <typename T>
  bool isTrackSelectedTPCCuts(T& track)
  {
    if (!b_useTrackSelections)
      return true; // no track selections applied
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kTPCNCls))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kTPCCrossedRows))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kTPCCrossedRowsOverNCls))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kTPCChi2NDF))
      return false;
    return true;
  }
  /// Function applying the ITS selections
  template <typename T>
  bool isTrackSelectedITSCuts(T& track)
  {
    if (!b_useTrackSelections)
      return true; // no track selections applied
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kITSChi2NDF))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kITSHits))
      return false;
    return true;
  }
  //
  //
  // fill collision control plots
  //
  template <bool IS_MC, typename T>
  void fillGeneralHistos(T& coll)
  {
    if constexpr (IS_MC) {
      histos.fill(HIST("MC/control/zPrimary"), coll.posZ());
      histos.fill(HIST("MC/control/yPrimary"), coll.posY());
      histos.fill(HIST("MC/control/xPrimary"), coll.posX());
      histos.fill(HIST("MC/control/chi2Prim"), coll.chi2());
    } else {
      histos.fill(HIST("data/control/zPrimary"), coll.posZ());
      histos.fill(HIST("data/control/yPrimary"), coll.posY());
      histos.fill(HIST("data/control/xPrimary"), coll.posX());
      histos.fill(HIST("data/control/chi2Prim"), coll.chi2());
    }
    return;
  }
  //
  // define global variables
  int count = 0;
  int countData = 0;
  int countNoMC = 0;
  int siPDGCode = 0;
  int tpPDGCode = 0;
  std::vector<int>::iterator itr_pdg;
  float pdg_fill = 0.0;
  //
  //
  //
  /******************************************************************/
  //
  //
  /////////////////////////////////////////////////////
  ///   Template function to perform the analysis   ///
  /////////////////////////////////////////////////////
  template <bool IS_MC, typename T, typename P, typename B>
  void fillHistograms(T& tracks, P& mcParticles, B const& bcs)
  {
    //
    float trackPt = 0; //, ITStrackPt = 0;
    //
    //
    float tpcNSigmaPion = -999.f;
    float tpcNSigmaKaon = -999.f;
    float tpcNSigmaProton = -999.f;
    float tofNSigmaPion = -999.f;
    float tofNSigmaKaon = -999.f;
    float tofNSigmaProton = -999.f;
    //
    //
    for (auto& track : tracks) {
      // choose if we keep the track according to the TRD presence requirement
      if ((isTRDThere == 1) && !track.hasTRD())
        continue;
      if ((isTRDThere == 0) && track.hasTRD())
        continue;

      // kinematic track seletions for all tracks
      if (!isTrackSelectedKineCuts(track))
        continue;

      if constexpr (IS_MC) {
        if (!track.has_mcParticle()) {
          countNoMC++;
          if (doDebug)
            LOGF(warning, " N. %d track without MC particle, skipping...", countNoMC);
          continue;
        }
      }

      //
      // pt from full tracking or from TPCinnerWallPt
      float reco_pt = track.pt();
      float tpcinner_pt = computePtInParamTPC(track);

      /// Using pt calculated at the inner wall of TPC
      /// Caveat: tgl still from tracking: this is not the value of tgl at the
      /// inner wall of TPC
      if (b_useTPCinnerWallPt)
        trackPt = tpcinner_pt;
      else
        trackPt = reco_pt;

      /// special case for ITS tracks
      /// Using pt calculated at the inner wall of TPC
      /// Caveat: tgl still from tracking: this is not the value of tgl at the
      /// inner wall of TPC
      // if (b_useTPCinnerWallPtForITS)
      //   ITStrackPt = tpcinner_pt;
      // else
      //   ITStrackPt = reco_pt;

      countData++;
      //
      //  keep sign of track as boolean
      bool positiveTrack = false;
      if (track.signed1Pt() > 0)
        positiveTrack = true;
      //
      // PID sigmas
      if constexpr (!IS_MC) {
        tpcNSigmaPion = track.tpcNSigmaPi();
        tpcNSigmaKaon = track.tpcNSigmaKa();
        tpcNSigmaProton = track.tpcNSigmaPr();
        tofNSigmaPion = track.tofNSigmaPi();
        tofNSigmaKaon = track.tofNSigmaKa();
        tofNSigmaProton = track.tofNSigmaPr();
      }
      const bool trkWTOF = track.hasTOF();
      const bool trkWTPC = track.hasTPC();
      const bool trkWITS = track.hasITS();
      bool pionPIDwithTPC = (nSigmaPID->get("TPC", "nSigPionMin") < tpcNSigmaPion && tpcNSigmaPion < nSigmaPID->get("TPC", "nSigPionMax"));
      bool pionPIDwithTOF = (nSigmaPID->get("TOF", "nSigPionMin") < tofNSigmaPion && tofNSigmaPion < nSigmaPID->get("TOF", "nSigPionMax"));
      bool kaonPIDwithTPC = (nSigmaPID->get("TPC", "nSigKaonMin") < tpcNSigmaKaon && tpcNSigmaKaon < nSigmaPID->get("TPC", "nSigKaonMax"));
      bool kaonPIDwithTOF = (nSigmaPID->get("TOF", "nSigKaonMin") < tofNSigmaKaon && tofNSigmaKaon < nSigmaPID->get("TOF", "nSigKaonMax"));
      bool protonPIDwithTPC = (nSigmaPID->get("TPC", "nSigProtonMin") < tpcNSigmaProton && tpcNSigmaProton < nSigmaPID->get("TPC", "nSigProtonMax"));
      bool protonPIDwithTOF = (nSigmaPID->get("TOF", "nSigProtonMin") < tofNSigmaProton && tofNSigmaProton < nSigmaPID->get("TOF", "nSigProtonMax"));
      // isPion
      bool isPion = false;
      if (isPIDPionRequired && pionPIDwithTPC && ((!trkWTOF) || pionPIDwithTOF))
        isPion = true;
      // isKaon
      bool isKaon = false;
      if (isPIDKaonRequired && kaonPIDwithTPC && ((!trkWTOF) || kaonPIDwithTOF))
        isKaon = true;
      // isProton
      bool isProton = false;
      if (isPIDProtonRequired && protonPIDwithTPC && ((!trkWTOF) || protonPIDwithTOF))
        isProton = true;
      //
      int sayPrim = -1, signPDGCode = -2, specind = 0;
      if constexpr (IS_MC) {
        auto mcpart = track.mcParticle();
        siPDGCode = mcpart.pdgCode();
        tpPDGCode = TMath::Abs(siPDGCode);
        if (mcpart.isPhysicalPrimary()) {
          histos.get<TH1>(HIST("MC/control/etahist_diff"))->Fill(mcpart.eta() - track.eta());
          auto delta = mcpart.phi() - track.phi();
          if (delta > PI) {
            delta -= TwoPI;
          }
          if (delta < -PI) {
            delta += TwoPI;
          }
          histos.get<TH1>(HIST("MC/control/phihist_diff"))->Fill(delta);
        }

        /// MC info for THnSparse filling
        sayPrim = 0;
        specind = 0;
        if (mcpart.isPhysicalPrimary())
          sayPrim = 0;
        else if (mcpart.getProcess() == 4)
          sayPrim = 1;
        else
          sayPrim = 2;
        signPDGCode = siPDGCode / tpPDGCode;
        switch (tpPDGCode) {
          case 11:
            specind = 1;
            break;
          case 211:
            specind = 2;
            break;
          case 321:
            specind = 3;
            break;
          case 2212:
            specind = 4;
            break;
          default:
            specind = 5;
        }
      }
      //
      //
      // fill thnsparse for fraction analysis
      if (makethn) {
        if constexpr (IS_MC) {
          histos.fill(HIST("MC/sparse/thnsforfrac"), track.dcaXY(), trackPt, track.phi(), track.eta(), sayPrim, signPDGCode, specind);
        } else {
          histos.fill(HIST("data/sparse/thnsforfrac"), track.dcaXY(), trackPt, track.phi(), track.eta(), sayPrim, signPDGCode, specind);
        }
      }
      //
      // all tracks, no conditions
      //
      if (trkWTPC && isTrackSelectedTPCCuts(track)) {
        if constexpr (IS_MC) { ////////////////////////   MC
          //
          histos.fill(HIST("MC/control/zDCA_tpc"), track.dcaZ());
          histos.fill(HIST("MC/control/xyDCA_tpc"), track.dcaXY());
          //
          if (makept2d) {
            histos.fill(HIST("MC/control/ptptconfTPCall"), reco_pt, tpcinner_pt);
            if (!trkWITS)
              histos.fill(HIST("MC/control/ptptconfTPCo"), reco_pt, tpcinner_pt);
          }
          histos.get<TH1>(HIST("MC/qopthist_tpc"))->Fill(track.signed1Pt());
          histos.get<TH1>(HIST("MC/pthist_tpc"))->Fill(trackPt);
          histos.get<TH1>(HIST("MC/phihist_tpc"))->Fill(track.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc"))->Fill(track.eta());
          if (trkWTOF) {
            histos.get<TH1>(HIST("MC/pthist_toftpc"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/phihist_toftpc"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/etahist_toftpc"))->Fill(track.eta());
          }
        } else { ////////////////////////   DATA
          //
          histos.fill(HIST("data/control/zDCA_tpc"), track.dcaZ());
          histos.fill(HIST("data/control/xyDCA_tpc"), track.dcaXY());
          //
          if (makept2d) {
            histos.fill(HIST("data/control/ptptconfTPCall"), reco_pt, tpcinner_pt);
            if (!trkWITS)
              histos.fill(HIST("data/control/ptptconfTPCo"), reco_pt, tpcinner_pt);
          }
          histos.get<TH1>(HIST("data/qopthist_tpc"))->Fill(track.signed1Pt());
          histos.get<TH1>(HIST("data/pthist_tpc"))->Fill(trackPt);
          histos.get<TH1>(HIST("data/phihist_tpc"))->Fill(track.phi());
          histos.get<TH1>(HIST("data/etahist_tpc"))->Fill(track.eta());
          //
          // monitoring vs. time (debug reasons)
          if (enableMonitorVsTime && timeMonitorSetUp) {
            if (track.has_collision()) {
              const auto timestamp = track.collision().template bc_as<BCsWithTimeStamp>().timestamp(); /// NB: in ms
              histos.get<TH1>(HIST("data/hTrkTPCvsTime"))->Fill(timestamp);
              if (enableTHnSparseMonitorVsTime) {
                histos.get<THnSparse>(HIST("data/hTrkTPCvsTimePtEtaPosZ"))->Fill(timestamp, trackPt, track.eta(), track.collision().posZ(), 1. / trackPt, positiveTrack ? 0.5 : -0.5, track.tpcNClsFound(), track.itsNCls());
              }
            }
          }
          //
          // with TOF tag
          if (trkWTOF) {
            histos.get<TH1>(HIST("data/pthist_toftpc"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/phihist_toftpc"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/etahist_toftpc"))->Fill(track.eta());
          }
          //
          // PID is applied
          if (isPion) {
            histos.get<TH1>(HIST("data/PID/zDCA_tpc_pi"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("data/PID/xyDCA_tpc_pi"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("data/PID/pthist_tpc_pi"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/PID/phihist_tpc_pi"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/PID/etahist_tpc_pi"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_piplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_piplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_piminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_piminus"))->Fill(track.eta());
            }
          }
          if (isPIDPionRequired) {
            if (pionPIDwithTPC) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pi_PIDTPC"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus_PIDTPC"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus_PIDTPC"))->Fill(trackPt);
              if (!trkWTOF || !pionPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_pi_PIDTPC_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus_PIDTPC_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus_PIDTPC_O"))->Fill(trackPt);
              }
            }
            if (trkWTOF && pionPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pi_PIDTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus_PIDTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus_PIDTOF"))->Fill(trackPt);
              if (!pionPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_pi_PIDTOF_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus_PIDTOF_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus_PIDTOF_O"))->Fill(trackPt);
              }
            }
            if (pionPIDwithTPC && trkWTOF && pionPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pi_PIDTPCTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piplus_PIDTPCTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_piminus_PIDTPCTOF"))->Fill(trackPt);
            }
          }
          // end pions
          if (isKaon) {
            histos.get<TH1>(HIST("data/PID/zDCA_tpc_ka"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("data/PID/xyDCA_tpc_ka"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("data/PID/pthist_tpc_ka"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/PID/phihist_tpc_ka"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/PID/etahist_tpc_ka"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_kaplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_kaplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_kaminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_kaminus"))->Fill(track.eta());
            }
          }
          if (isPIDKaonRequired) {
            if (kaonPIDwithTPC) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_ka_PIDTPC"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus_PIDTPC"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus_PIDTPC"))->Fill(trackPt);
              if (!trkWTOF || !kaonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_ka_PIDTPC_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus_PIDTPC_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus_PIDTPC_O"))->Fill(trackPt);
              }
            }
            if (trkWTOF && kaonPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_ka_PIDTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus_PIDTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus_PIDTOF"))->Fill(trackPt);
              if (!kaonPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_ka_PIDTOF_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus_PIDTOF_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus_PIDTOF_O"))->Fill(trackPt);
              }
            }
            if (kaonPIDwithTPC && trkWTOF && kaonPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_ka_PIDTPCTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaplus_PIDTPCTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_kaminus_PIDTPCTOF"))->Fill(trackPt);
            }
          }
          // end kaons
          if (isProton) {
            histos.get<TH1>(HIST("data/PID/zDCA_tpc_pr"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("data/PID/xyDCA_tpc_pr"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("data/PID/pthist_tpc_pr"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/PID/phihist_tpc_pr"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/PID/etahist_tpc_pr"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_prplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_prplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_prminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_prminus"))->Fill(track.eta());
            }
          }
          if (isPIDProtonRequired) {
            if (protonPIDwithTPC) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pr_PIDTPC"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus_PIDTPC"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus_PIDTPC"))->Fill(trackPt);
              if (!trkWTOF || !protonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_pr_PIDTPC_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus_PIDTPC_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus_PIDTPC_O"))->Fill(trackPt);
              }
            }
            if (trkWTOF && protonPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pr_PIDTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus_PIDTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus_PIDTOF"))->Fill(trackPt);
              if (!protonPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpc_pr_PIDTOF_O"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus_PIDTOF_O"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus_PIDTOF_O"))->Fill(trackPt);
              }
            }
            if (protonPIDwithTPC && trkWTOF && protonPIDwithTOF) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_pr_PIDTPCTOF"))->Fill(trackPt);
              if (positiveTrack)
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prplus_PIDTPCTOF"))->Fill(trackPt);
              else
                histos.get<TH1>(HIST("data/PID/pthist_tpc_prminus_PIDTPCTOF"))->Fill(trackPt);
            }
          }
          // end protons
          if (!isPion && !isKaon && !isProton && (isPIDPionRequired || isPIDKaonRequired || isPIDProtonRequired)) {
            histos.get<TH1>(HIST("data/PID/zDCA_tpc_noid"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("data/PID/xyDCA_tpc_noid"))->Fill(track.dcaXY());
            histos.get<TH1>(HIST("data/PID/pthist_tpc_noid"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/PID/phihist_tpc_noid"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/PID/etahist_tpc_noid"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_noidplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_noidplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_noidplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/PID/pthist_tpc_noidminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpc_noidminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpc_noidminus"))->Fill(track.eta());
            }
          } // not pions, nor kaons, nor protons
        }   // end if DATA
        //
        if (trkWITS && isTrackSelectedITSCuts(track)) { ////////////////////////////////////////////   ITS tag inside TPC tagged
          if constexpr (IS_MC) {                        ////////////////////////   MC
            histos.get<TH1>(HIST("MC/control/zDCA_tpcits"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/control/xyDCA_tpcits"))->Fill(track.dcaXY());
            if (makept2d)
              histos.fill(HIST("MC/control/ptptconfTPCITS"), reco_pt, tpcinner_pt);
            histos.get<TH1>(HIST("MC/qopthist_tpcits"))->Fill(track.signed1Pt());
            histos.get<TH1>(HIST("MC/pthist_tpcits"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/phihist_tpcits"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits"))->Fill(track.eta());
            if (trkWTOF) {
              histos.get<TH1>(HIST("MC/pthist_toftpcits"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/phihist_toftpcits"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/etahist_toftpcits"))->Fill(track.eta());
            }
          } else { ////////////////////////   DATA
            histos.get<TH1>(HIST("data/control/zDCA_tpcits"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("data/control/xyDCA_tpcits"))->Fill(track.dcaXY());
            if (makept2d)
              histos.fill(HIST("data/control/ptptconfTPCITS"), reco_pt, tpcinner_pt);
            histos.get<TH1>(HIST("data/qopthist_tpcits"))->Fill(track.signed1Pt());
            //
            //  PID is applied
            if (isPion) {
              histos.get<TH1>(HIST("data/PID/zDCA_tpcits_pi"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("data/PID/xyDCA_tpcits_pi"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpcits_pi"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpcits_pi"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_piplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_piplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_piminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_piminus"))->Fill(track.eta());
              }
            }
            if (isPIDPionRequired) {
              if (pionPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi_PIDTPC"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus_PIDTPC"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus_PIDTPC"))->Fill(trackPt);
                if (!trkWTOF || !pionPIDwithTOF) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi_PIDTPC_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus_PIDTPC_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus_PIDTPC_O"))->Fill(trackPt);
                }
              }
              if (trkWTOF && pionPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi_PIDTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus_PIDTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus_PIDTOF"))->Fill(trackPt);
                if (!pionPIDwithTPC) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi_PIDTOF_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus_PIDTOF_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus_PIDTOF_O"))->Fill(trackPt);
                }
              }
              if (pionPIDwithTPC && trkWTOF && pionPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pi_PIDTPCTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piplus_PIDTPCTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_piminus_PIDTPCTOF"))->Fill(trackPt);
              }
            }
            // end pions
            if (isKaon) {
              histos.get<TH1>(HIST("data/PID/zDCA_tpcits_ka"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("data/PID/xyDCA_tpcits_ka"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpcits_ka"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpcits_ka"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_kaplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_kaplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_kaminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_kaminus"))->Fill(track.eta());
              }
            }
            if (isPIDKaonRequired) {
              if (kaonPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka_PIDTPC"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus_PIDTPC"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus_PIDTPC"))->Fill(trackPt);
                if (!trkWTOF || !kaonPIDwithTOF) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka_PIDTPC_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus_PIDTPC_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus_PIDTPC_O"))->Fill(trackPt);
                }
              }
              if (trkWTOF && kaonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka_PIDTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus_PIDTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus_PIDTOF"))->Fill(trackPt);
                if (!kaonPIDwithTPC) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka_PIDTOF_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus_PIDTOF_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus_PIDTOF_O"))->Fill(trackPt);
                }
              }
              if (kaonPIDwithTPC && trkWTOF && kaonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_ka_PIDTPCTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaplus_PIDTPCTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_kaminus_PIDTPCTOF"))->Fill(trackPt);
              }
            }
            // end kaons
            if (isProton) {
              histos.get<TH1>(HIST("data/PID/zDCA_tpcits_pr"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("data/PID/xyDCA_tpcits_pr"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpcits_pr"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpcits_pr"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_prplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_prplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_prminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_prminus"))->Fill(track.eta());
              }
            }
            if (isPIDProtonRequired) {
              if (protonPIDwithTPC) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr_PIDTPC"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus_PIDTPC"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus_PIDTPC"))->Fill(trackPt);
                if (!trkWTOF || !protonPIDwithTOF) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr_PIDTPC_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus_PIDTPC_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus_PIDTPC_O"))->Fill(trackPt);
                }
              }
              if (trkWTOF && protonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr_PIDTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus_PIDTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus_PIDTOF"))->Fill(trackPt);
                if (!protonPIDwithTPC) {
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr_PIDTOF_O"))->Fill(trackPt);
                  if (positiveTrack)
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus_PIDTOF_O"))->Fill(trackPt);
                  else
                    histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus_PIDTOF_O"))->Fill(trackPt);
                }
              }
              if (protonPIDwithTPC && trkWTOF && protonPIDwithTOF) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_pr_PIDTPCTOF"))->Fill(trackPt);
                if (positiveTrack)
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prplus_PIDTPCTOF"))->Fill(trackPt);
                else
                  histos.get<TH1>(HIST("data/PID/pthist_tpcits_prminus_PIDTPCTOF"))->Fill(trackPt);
              }
            }
            // end protons
            histos.get<TH1>(HIST("data/pthist_tpcits"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/phihist_tpcits"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits"))->Fill(track.eta());
            //
            // monitoring vs. time (debug reasons)
            if (enableMonitorVsTime && timeMonitorSetUp) {
              if (track.has_collision()) {
                const auto timestamp = track.collision().template bc_as<BCsWithTimeStamp>().timestamp(); /// NB: in ms
                histos.get<TH1>(HIST("data/hTrkITSTPCvsTime"))->Fill(timestamp);
                if (enableTHnSparseMonitorVsTime) {
                  histos.get<THnSparse>(HIST("data/hTrkITSTPCvsTimePtEtaPosZ"))->Fill(timestamp, trackPt, track.eta(), track.collision().posZ(), 1. / trackPt, positiveTrack ? 0.5 : -0.5, track.tpcNClsFound(), track.itsNCls());
                }
              }
            }
            //
            // not identified
            if (!isPion && !isKaon && !isProton) {
              histos.get<TH1>(HIST("data/PID/zDCA_tpcits_noid"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("data/PID/xyDCA_tpcits_noid"))->Fill(track.dcaXY());
              histos.get<TH1>(HIST("data/PID/pthist_tpcits_noid"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/PID/phihist_tpcits_noid"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/PID/etahist_tpcits_noid"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_noidplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_noidplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_noidplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("data/PID/pthist_tpcits_noidminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("data/PID/phihist_tpcits_noidminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("data/PID/etahist_tpcits_noidminus"))->Fill(track.eta());
              }
            }
            //
            // with TOF tag
            if (trkWTOF) {
              histos.get<TH1>(HIST("data/pthist_toftpcits"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/phihist_toftpcits"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/etahist_toftpcits"))->Fill(track.eta());
            }
          }
          /// control plot: correlation # ITS its vs ITS layer
          int itsNhits = 0;
          for (unsigned int i = 0; i < 7; i++) {
            if (track.itsClusterMap() & (1 << i)) {
              itsNhits += 1;
            }
          }
          bool trkHasITS = false;
          for (unsigned int i = 0; i < 7; i++) {
            if (track.itsClusterMap() & (1 << i)) {
              trkHasITS = true;
              if (IS_MC) { ////////////////////////   MC
                histos.fill(HIST("MC/control/itsHitsMatched"), i, itsNhits);
              } else { ////////////////////////   DATA
                histos.fill(HIST("data/control/itsHitsMatched"), i, itsNhits);
              }
            }
          }
          if (!trkHasITS) {
            if (IS_MC) { ////////////////////////   MC
              histos.fill(HIST("MC/control/itsHitsMatched"), -1, itsNhits);
            } else { ////////////////////////   DATA
              histos.fill(HIST("data/control/itsHitsMatched"), -1, itsNhits);
            }
          }
        } //  end if ITS
      }   //  end if TPC
      //
      // all tracks with pt>0.5
      if (trackPt > 0.5) {
        if (trkWTPC && isTrackSelectedTPCCuts(track)) {
          if constexpr (IS_MC) { ////////////////////////   MC
            histos.get<TH1>(HIST("MC/pthist_tpc_05"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/phihist_tpc_05"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_05"))->Fill(track.eta());
          } else { ////////////////////////   DATA
            histos.get<TH1>(HIST("data/pthist_tpc_05"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/phihist_tpc_05"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/etahist_tpc_05"))->Fill(track.eta());
          }
          if (trkWITS && isTrackSelectedITSCuts(track)) {
            if constexpr (IS_MC) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_05"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/phihist_tpcits_05"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_05"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/pthist_tpcits_05"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/phihist_tpcits_05"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/etahist_tpcits_05"))->Fill(track.eta());
            }
          } //  end if ITS
        }   //  end if TPC
      }     //  end if pt > 0.5
      //
      // positive only
      if (track.signed1Pt() > 0) {
        if (trkWTPC && isTrackSelectedTPCCuts(track)) {
          if constexpr (IS_MC) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pos"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/phihist_tpc_pos"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pos"))->Fill(track.eta());
          } else {
            histos.get<TH1>(HIST("data/pthist_tpc_pos"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/phihist_tpc_pos"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/etahist_tpc_pos"))->Fill(track.eta());
          }
          if (trkWITS && isTrackSelectedITSCuts(track)) {
            if constexpr (IS_MC) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pos"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/phihist_tpcits_pos"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pos"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/pthist_tpcits_pos"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/phihist_tpcits_pos"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/etahist_tpcits_pos"))->Fill(track.eta());
            }
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end positive
      //
      // negative only
      if (track.signed1Pt() < 0) {
        if (trkWTPC && isTrackSelectedTPCCuts(track)) {
          if constexpr (IS_MC) {
            histos.get<TH1>(HIST("MC/pthist_tpc_neg"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/phihist_tpc_neg"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_neg"))->Fill(track.eta());
          } else {
            histos.get<TH1>(HIST("data/pthist_tpc_neg"))->Fill(trackPt);
            histos.get<TH1>(HIST("data/phihist_tpc_neg"))->Fill(track.phi());
            histos.get<TH1>(HIST("data/etahist_tpc_neg"))->Fill(track.eta());
          }
          if (trkWITS && isTrackSelectedITSCuts(track)) {
            if constexpr (IS_MC) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_neg"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/phihist_tpcits_neg"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_neg"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("data/pthist_tpcits_neg"))->Fill(trackPt);
              histos.get<TH1>(HIST("data/phihist_tpcits_neg"))->Fill(track.phi());
              histos.get<TH1>(HIST("data/etahist_tpcits_neg"))->Fill(track.eta());
            }
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end negative

      if constexpr (IS_MC) { // MC
        auto mcpart = track.mcParticle();
        //
        // only primaries
        if (mcpart.isPhysicalPrimary()) {
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/primsec/zDCA_tpc_prim"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/primsec/xyDCA_tpc_prim"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/primsec/qopthist_tpc_prim"))->Fill(track.signed1Pt());
            histos.get<TH1>(HIST("MC/primsec/pthist_tpc_prim"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/primsec/phihist_tpc_prim"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/primsec/etahist_tpc_prim"))->Fill(track.eta());
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/primsec/zDCA_tpcits_prim"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/primsec/xyDCA_tpcits_prim"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/primsec/qopthist_tpcits_prim"))->Fill(track.signed1Pt());
              histos.get<TH1>(HIST("MC/primsec/pthist_tpcits_prim"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/primsec/phihist_tpcits_prim"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/primsec/etahist_tpcits_prim"))->Fill(track.eta());
            } //  end if ITS
          }   //  end if TPC
          //  end if primaries
        } else if (mcpart.getProcess() == 4) {
          //
          // only secondaries from decay
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/primsec/zDCA_tpc_secd"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/primsec/xyDCA_tpc_secd"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/primsec/qopthist_tpc_secd"))->Fill(track.signed1Pt());
            histos.get<TH1>(HIST("MC/primsec/pthist_tpc_secd"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/primsec/phihist_tpc_secd"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/primsec/etahist_tpc_secd"))->Fill(track.eta());
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/primsec/zDCA_tpcits_secd"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/primsec/xyDCA_tpcits_secd"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/primsec/qopthist_tpcits_secd"))->Fill(track.signed1Pt());
              histos.get<TH1>(HIST("MC/primsec/pthist_tpcits_secd"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/primsec/phihist_tpcits_secd"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/primsec/etahist_tpcits_secd"))->Fill(track.eta());
            } //  end if ITS
          }   //  end if TPC
          // end if secondaries from decay
        } else {
          //
          // only secondaries from material
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/primsec/zDCA_tpc_secm"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/primsec/xyDCA_tpc_secm"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/primsec/qopthist_tpc_secm"))->Fill(track.signed1Pt());
            histos.get<TH1>(HIST("MC/primsec/pthist_tpc_secm"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/primsec/phihist_tpc_secm"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/primsec/etahist_tpc_secm"))->Fill(track.eta());
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/primsec/zDCA_tpcits_secm"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/primsec/xyDCA_tpcits_secm"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/primsec/qopthist_tpcits_secm"))->Fill(track.signed1Pt());
              histos.get<TH1>(HIST("MC/primsec/pthist_tpcits_secm"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/primsec/phihist_tpcits_secm"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/primsec/etahist_tpcits_secm"))->Fill(track.eta());
            } //  end if ITS
          }   //  end if TPC
        }     // end if secondaries from material
        //
        // protons only
        if (tpPDGCode == 2212) {
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/PID/zDCA_tpc_pr"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/PID/xyDCA_tpc_pr"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/PID/pthist_tpc_pr"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/PID/phihist_tpc_pr"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/PID/etahist_tpc_pr"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_prplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_prplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_prplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_prminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_prminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_prminus"))->Fill(track.eta());
            }
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/zDCA_tpcits_pr"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/PID/xyDCA_tpcits_pr"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/PID/pthist_tpcits_pr"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpcits_pr"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpcits_pr"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_prplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_prplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_prplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_prminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_prminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_prminus"))->Fill(track.eta());
              }
            } //  end if ITS
          }   //  end if TPC
        }
        //
        // pions only
        if (tpPDGCode == 211) {
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/PID/zDCA_tpc_pi"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/PID/xyDCA_tpc_pi"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/PID/pthist_tpc_pi"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/PID/phihist_tpc_pi"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/PID/etahist_tpc_pi"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_piplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_piplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_piplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_piminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_piminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_piminus"))->Fill(track.eta());
            }
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/zDCA_tpcits_pi"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/PID/xyDCA_tpcits_pi"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/PID/pthist_tpcits_pi"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpcits_pi"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpcits_pi"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_piplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_piplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_piplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_piminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_piminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_piminus"))->Fill(track.eta());
              }
            } //  end if ITS
          }   //  end if TPC
          //
          // only primary pions
          if (mcpart.isPhysicalPrimary()) {
            if (trkWTPC && isTrackSelectedTPCCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_pi_prim"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_pi_prim"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_pi_prim"))->Fill(track.eta());
              if (trkWITS && isTrackSelectedITSCuts(track)) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_pi_prim"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_pi_prim"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_pi_prim"))->Fill(track.eta());
              } //  end if ITS
            }   //  end if TPC
            //  end if primaries
          } else if (mcpart.getProcess() == 4) {
            //
            // only secondary pions from decay
            if (trkWTPC && isTrackSelectedTPCCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_pi_secd"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_pi_secd"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_pi_secd"))->Fill(track.eta());
              if (trkWITS && isTrackSelectedITSCuts(track)) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_pi_secd"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_pi_secd"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_pi_secd"))->Fill(track.eta());
              } //  end if ITS
            }   //  end if TPC
            // end if secondaries from decay
          } else {
            //
            // only secondary pions from material
            if (trkWTPC && isTrackSelectedTPCCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_pi_secm"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_pi_secm"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_pi_secm"))->Fill(track.eta());
              if (trkWITS && isTrackSelectedITSCuts(track)) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_pi_secm"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_pi_secm"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_pi_secm"))->Fill(track.eta());
              } //  end if ITS
            }   //  end if TPC
          }     // end if secondaries from material  //
        }       // end pions only
        //
        // no primary/sec-d pions
        if (!((tpPDGCode == 211) && (mcpart.isPhysicalPrimary()))) {
          // gets the pdg code and finds its index in our vector
          itr_pdg = std::find(pdgChoice.begin(), pdgChoice.end(), tpPDGCode);
          if (itr_pdg != pdgChoice.cend())
            // index from zero, so increase by 1 to put in the right bin (and
            // 0.5 not needed but just not to sit in the edge)
            pdg_fill = static_cast<float>(std::distance(pdgChoice.begin(), itr_pdg)) + 1.5;
          else
            pdg_fill = -10.0;
          //
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/PID/pthist_tpc_nopi"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/PID/phihist_tpc_nopi"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/PID/etahist_tpc_nopi"))->Fill(track.eta());
            histos.get<TH1>(HIST("MC/PID/pdghist_den"))->Fill(pdg_fill);
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpcits_nopi"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpcits_nopi"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpcits_nopi"))->Fill(track.eta());
              histos.get<TH1>(HIST("MC/PID/pdghist_num"))->Fill(pdg_fill);
            } //  end if ITS
          }   //  end if TPC
        }     // end if not prim/sec-d pi
        //
        // kaons only
        if (tpPDGCode == 321) {
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/PID/zDCA_tpc_ka"))->Fill(track.dcaZ());
            histos.get<TH1>(HIST("MC/PID/xyDCA_tpc_ka"))->Fill(track.dcaXY());
            //
            histos.get<TH1>(HIST("MC/PID/pthist_tpc_ka"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/PID/phihist_tpc_ka"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/PID/etahist_tpc_ka"))->Fill(track.eta());
            if (positiveTrack) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_kaplus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_kaplus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_kaplus"))->Fill(track.eta());
            } else {
              histos.get<TH1>(HIST("MC/PID/pthist_tpc_kaminus"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpc_kaminus"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpc_kaminus"))->Fill(track.eta());
            }
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/zDCA_tpcits_ka"))->Fill(track.dcaZ());
              histos.get<TH1>(HIST("MC/PID/xyDCA_tpcits_ka"))->Fill(track.dcaXY());
              //
              histos.get<TH1>(HIST("MC/PID/pthist_tpcits_ka"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpcits_ka"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpcits_ka"))->Fill(track.eta());
              if (positiveTrack) {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_kaplus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_kaplus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_kaplus"))->Fill(track.eta());
              } else {
                histos.get<TH1>(HIST("MC/PID/pthist_tpcits_kaminus"))->Fill(trackPt);
                histos.get<TH1>(HIST("MC/PID/phihist_tpcits_kaminus"))->Fill(track.phi());
                histos.get<TH1>(HIST("MC/PID/etahist_tpcits_kaminus"))->Fill(track.eta());
              }
            } //  end if ITS
          }   //  end if TPC
        }
        //
        // pions and kaons together
        if (tpPDGCode == 211 || tpPDGCode == 321) {
          if (trkWTPC && isTrackSelectedTPCCuts(track)) {
            histos.get<TH1>(HIST("MC/PID/pthist_tpc_piK"))->Fill(trackPt);
            histos.get<TH1>(HIST("MC/PID/phihist_tpc_piK"))->Fill(track.phi());
            histos.get<TH1>(HIST("MC/PID/etahist_tpc_piK"))->Fill(track.eta());
            if (trkWITS && isTrackSelectedITSCuts(track)) {
              histos.get<TH1>(HIST("MC/PID/pthist_tpcits_piK"))->Fill(trackPt);
              histos.get<TH1>(HIST("MC/PID/phihist_tpcits_piK"))->Fill(track.phi());
              histos.get<TH1>(HIST("MC/PID/etahist_tpcits_piK"))->Fill(track.eta());
            } //  end if ITS
          }   //  end if TPC
        }
      }
      //
      //
    } //  end loop on tracks
    //
    //
    if (doDebug) {
      LOGF(info, "Selected tracks: %d ", countData);
      LOGF(info, "Selected tracks with MC: %d, tracks w/o MC: %d ", countData, countNoMC);
    }
  }

  /// Function to add histograms for time-based monitoring of matching efficiency (debug purposes)
  /// BEWARE: these hitograms are ok only if looked run-by-run
  void setUpTimeMonitoring(BCsWithTimeStamp const& bcs)
  {
    int runNumber = bcs.rawIteratorAt(0).runNumber();
    if (lastRunNumber != runNumber) {
      /// execute the code in this scope only once, i.e. when the current run is considered for the first time in this DF
      lastRunNumber = runNumber;
      int64_t tsSOR = 0;
      int64_t tsEOR = 0;

      /// reject AO2Ds for which no CCDB access is possible
      if (runNumber < 500000) {
        LOG(warning) << ">>> run number " << runNumber << " < 500000. access to CCDB not possible. Exiting";
        return;
      }

      /// If we are here, the current run was never considered before.
      /// Let's add the TH2 that we need for the monitoring
      /// Let's define the x-axis according to the start-of-run (SOR) and end-of-run (EOR) times
      o2::ccdb::CcdbApi ccdb_api;
      ccdb_api.init(ccdburl);
      std::map<string, string> metadataRCT, headers;
      headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
      tsSOR = atol(headers["SOR"].c_str());
      tsEOR = atol(headers["EOR"].c_str());
      double minMilliSec = floor(tsSOR); /// round tsSOR to the highest integer lower than tsSOR
      double maxMilliSec = ceil(tsEOR);  /// round tsEOR to the lowest integer higher than tsEOR
      const AxisSpec axisSeconds{static_cast<int>((maxMilliSec - minMilliSec) * 1. / 10.), minMilliSec, maxMilliSec, "time from January 1st, 1970 at UTC (unit: 10 ms)"};

      /// add histograms now
      const AxisSpec axisPtVsTime{ptBinsVsTime, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec axis1overPtVsTime{ptInvserseBinsVsTime, "1/#it{p}_{T} (#it{c}/GeV)"};
      const AxisSpec axisChargeVsTime{2, -1, 1, "charge"};
      const AxisSpec axisEtaVsTime{etaBinsVsTime, "#eta"};
      const AxisSpec axisPosZVsTime{posZBinsVsTime, "posZ (cm)"};
      const AxisSpec axisTpcClstVsTime{tpcClstBinsVsTime, "TPC clusters"};
      const AxisSpec axisItsClstVsTime{itsClstBinsVsTime, "TPC clusters"};
      // TPC tracks
      histos.add("data/hTrkTPCvsTime", "", kTH1D, {axisSeconds});
      if (enableTHnSparseMonitorVsTime) {
        histos.add("data/hTrkTPCvsTimePtEtaPosZ", "", kTHnSparseF, {axisSeconds, axisPtVsTime, axisEtaVsTime, axisPosZVsTime, axis1overPtVsTime, axisChargeVsTime, axisTpcClstVsTime, axisItsClstVsTime});
      }
      // ITS-TPC tracks
      histos.add("data/hTrkITSTPCvsTime", "", kTH1D, {axisSeconds});
      if (enableTHnSparseMonitorVsTime) {
        histos.add("data/hTrkITSTPCvsTimePtEtaPosZ", "", kTHnSparseF, {axisSeconds, axisPtVsTime, axisEtaVsTime, axisPosZVsTime, axis1overPtVsTime, axisChargeVsTime, axisTpcClstVsTime, axisItsClstVsTime});
      }

      /// time monitoring correctly set up
      timeMonitorSetUp = true;
    }
  }

  //////////////////////////////////////////////
  ///   Process MC with collision grouping   ///
  //////////////////////////////////////////////
  void processMC(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true>(tracks, mcParticles, mcParticles); /// 3rd argument non-sense in this case
    fillGeneralHistos<true>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processMC, "process MC", false);

  void processMCFilteredTracks(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::McTrackLabels>> const& tracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true>(tracks, mcParticles, mcParticles); /// 3rd argument non-sense in this case
    fillGeneralHistos<true>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processMCFilteredTracks, "process MC with filtered tracks with filterbit selections", false);

  ////////////////////////////////////////////////////////////
  ///   Process MC with collision grouping and IU tracks   ///
  ////////////////////////////////////////////////////////////
  void processTrkIUMC(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true>(tracks, mcParticles, mcParticles); /// 3rd argument non-sense in this case
    fillGeneralHistos<true>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processTrkIUMC, "process MC for IU tracks", false);

  /////////////////////////////////////////////
  ///   Process MC w/o collision grouping   ///
  /////////////////////////////////////////////
  void processMCNoColl(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& tracks, aod::McParticles const& mcParticles)
  {
    fillHistograms<true>(tracks, mcParticles, mcParticles); /// 3rd argument non-sense in this case
  }
  PROCESS_SWITCH(qaMatchEff, processMCNoColl, "process MC - no collision grouping", false);

  ////////////////////////////////////////////////
  ///   Process data with collision grouping   ///
  ////////////////////////////////////////////////
  void processData(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks, BCsWithTimeStamp const& bcs)
  {
    if (enableMonitorVsTime) {
      // tracks.rawIteratorAt(0).collision().bc_as<BCsWithTimeStamp>().timestamp(); /// NB: in ms
      setUpTimeMonitoring(bcs);
    }
    fillHistograms<false>(tracks, tracks, bcs); // 2nd argument not used in this case
    fillGeneralHistos<false>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processData, "process data", true);

  void processDataFilteredTracks(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>> const& tracks, BCsWithTimeStamp const& bcs)
  {
    if (enableMonitorVsTime) {
      // tracks.rawIteratorAt(0).collision().bc_as<BCsWithTimeStamp>().timestamp(); /// NB: in ms
      setUpTimeMonitoring(bcs);
    }
    fillHistograms<false>(tracks, tracks, bcs); // 2nd argument not used in this case
    fillGeneralHistos<false>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processDataFilteredTracks, "process data with filtered tracks with filterbit selections", false);

  /////////////////////////////////////////////////////////////
  ///   Process data with collision grouping and IU tracks  ///
  /////////////////////////////////////////////////////////////
  void processTrkIUData(soa::Filtered<aod::Collisions>::iterator const& collision, soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks)
  {
    fillHistograms<false>(tracks, tracks, tracks); // 2nd and 3rd arguments not used in this case
    fillGeneralHistos<false>(collision);
  }
  PROCESS_SWITCH(qaMatchEff, processTrkIUData, "process data", false);

  ///////////////////////////////////////////////
  ///   Process data w/o collision grouping   ///
  ///////////////////////////////////////////////
  void processDataNoColl(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr> const& tracks, BCsWithTimeStamp const& bcs)
  {
    if (enableMonitorVsTime) {
      setUpTimeMonitoring(bcs);
    }
    fillHistograms<false>(tracks, tracks, bcs); // 2nd argument not used in this case
  }
  PROCESS_SWITCH(qaMatchEff, processDataNoColl, "process data - no collision grouping", true);

}; // end of structure

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<qaMatchEff>(cfgc, TaskName{"qa-match-eff"})};
}
