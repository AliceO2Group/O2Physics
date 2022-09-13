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
/// \brief A task for basic checks on ITS-TPC track matching efficiency
/// \author Rosario Turrisi (rosario.turrisi@pd.infn.it)
/// \author Mattia Faggin (mfaggin@cern.ch)
/// \since 2022

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

//
// base namespaces
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
//
struct qaMatchEff {
  //
  // histogram registry
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  //
  // Track selections
  Configurable<bool> b_useTrackSelections{"b_useTrackSelections", false, "Boolean to switch the track selections on/off."};
  // kinematics
  Configurable<float> ptMinCut{"ptMinCut", 0.1f, "Minimum transverse momentum (GeV/c)"};
  Configurable<float> ptMaxCut{"ptMaxCut", 100.f, "Maximum transverse momentum (GeV/c)"};
  Configurable<float> etaMinCut{"etaMinCut", -2.0f, "Minimum pseudorapidity"};
  Configurable<float> etaMaxCut{"etaMaxCut", 2.0f, "Maximum pseudorapidity"};
  // TPC
  Configurable<int> tpcNClusterMin{"tpcNClusterMin", 0, "Minimum number of clusters in TPC"};
  Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 70, "Minimum number of crossed rows in TPC"};
  Configurable<float> tpcNCrossedRowsOverFindableClstMin{"tpcNCrossedRowsOverFindableClstMin", 0.8f, "Minimum fracion of crossed rows over findable custers in TPC"};
  Configurable<float> tpcChi2Max{"tpcChi2Max", 4.0f, "Maximum chi2 in TPC"};
  // ITS
  Configurable<float> itsChi2Max{"itsChi2Max", 36.0f, "Maximum chi2 in ITS"};
  Configurable<int> customITShitmap{"customITShitmap", 3, "ITS hitmap (think to the binary representation)"};
  Configurable<int> customMinITShits{"customMinITShits", 1, "Minimum number of layers crossed by a track among those in \"customITShitmap\""}; //
  //
  Configurable<bool> isitMC{"isitMC", false, "Reading MC files, data if false"};
  Configurable<bool> doDebug{"doDebug", false, "Flag of debug information"};
  // Histogram configuration
  //
  // histo x axes limits
  Configurable<float> etaMin{"eta-min", -2.0f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 2.0f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.0f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 1.0f * TwoPI, "Upper limit in phi"};
  // histos bins
  Configurable<int> etaBins{"eta-bins", 40, "Number of eta bins"};
  Configurable<int> phiBins{"phi-bins", 18, "Number of phi bins"};
  //
  // special histo, few particles explicitly stored, then pdg>3000
  Configurable<int> pdgBins{"pdg-bins", 14, "Number of pdg values counted"};
  //
  // histo axes
  //
  // non uniform pt binning
  std::vector<double> ptBinning = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
  //
  AxisSpec axisPDG{pdgBins, 0, pdgBins + 1.000, "pdgclass"};
  //
  AxisSpec axisPt{ptBinning, "#it{p}_{T} (GeV/#it{c})"};
  //
  AxisSpec axisEta{etaBins, etaMin, etaMax, "#eta"};
  AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};
  AxisSpec axisDEta{etaBins, etaMin, etaMax, "D#eta"};
  AxisSpec axisDPh{phiBins, -PI, PI, "D#it{#varphi} (rad)"};
  //
  // pdg codes vector
  std::vector<int> pdgChoice = {211, 213, 215, 217, 219, 221, 223, 321, 411, 521, 2212, 1114, 2214};
  //
  //
  // Track selection object
  TrackSelection cutObject;
  // Init function
  //
  void init(InitContext&)
  {
    if (doDebug)
      LOG(info) << "===========================================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  is it MC? = " << isitMC;
    //
    // let's know if it's MC or data
    if (isitMC)
      initMC();
    else
      initData();

    if ((!isitMC && (doprocessMC || doprocessMCNoColl)) || (isitMC && (doprocessData && doprocessDataNoColl)))
      LOGF(fatal, "Initialization set for MC and processData function flagged (or viceversa)! Fix the configuration.");
    if ((doprocessMC && doprocessMCNoColl) || (doprocessData && doprocessDataNoColl))
      LOGF(fatal, "Cannot process for both without collision tag and with collision tag at the same time! Fix the configuration.");

    /// initialize the track selections
    if (b_useTrackSelections) {
      // kinematics
      cutObject.SetEtaRange(etaMinCut, etaMaxCut);
      cutObject.SetPtRange(ptMinCut, ptMaxCut);
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
        LOG(info) << "Layer " << (int)(*it) << " ";
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
    // tpc request and tpc+its request for all, positive and negative charges vs pt, phi, eta (18 histos tot)
    histos.add("data/pthist_tpc", "#it{p}_{T} distribution - data TPC tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc", "#eta distribution - data TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc", "#phi distribution - data TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits", "#it{p}_{T} distribution - data TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits", "#eta distribution - data TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits", "#phi distribution - data TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpc_pos", "#it{p}_{T} distribution - data q>0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc_pos", "#eta distribution - data q>0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc_pos", "#phi distribution - data q>0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits_pos", "#it{p}_{T} distribution - data q>0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits_pos", "#eta distribution - data q>0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits_pos", "#phi distribution - data q>0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpc_neg", "#it{p}_{T} distribution - data q<0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc_neg", "#eta distribution - data q<0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc_neg", "#phi distribution - data q<0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits_neg", "#it{p}_{T} distribution - data q<0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits_neg", "#eta distribution - data q<0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits_neg", "#phi distribution - data q<0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // pt>0.5 GeV/c threshold
    histos.add("data/pthist_tpc_05", "#it{p}_{T} distribution - data TPC tag, #it{p}_{T}>0.5", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc_05", "#eta distribution - data TPC tag, #it{p}_{T}>0.5", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc_05", "#phi distribution - data TPC tag, #it{p}_{T}>0.5", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits_05", "#it{p}_{T} distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits_05", "#eta distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits_05", "#phi distribution - data TPC+ITS tag #it{p}_{T}>0.5", kTH1F, {axisPhi}, true);
  }
  //
  // Init MC function
  void initMC()
  {
    if (doDebug)
      LOGF(info, " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    MC  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

    //
    // adding histos to the registry
    // data histos
    // tpc request and tpc+its request for all, positive and negative charges
    // and for phys. primaries, decay secondaries and mat. secondaries (both charges) vs pt, phi, eta (36 histos tot)
    // pions only, also split in prim secd secm
    //
    // all, positive, negative
    histos.add("MC/pthist_tpc", "#it{p}_{T} distribution - MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc", "#eta distribution - MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc", "#phi distribution - MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits", "#it{p}_{T} distribution - MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits", "#eta distribution - MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits", "#phi distribution - MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pos", "#it{p}_{T} distribution - MC q>0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pos", "#eta distribution - MC q>0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pos", "#phi distribution - MC q>0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pos", "#it{p}_{T} distribution - MC q>0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pos", "#eta distribution - MC q>0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pos", "#phi distribution - MC q>0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_neg", "#it{p}_{T} distribution - MC q<0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_neg", "#eta distribution - MC q<0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_neg", "#phi distribution - MC q<0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_neg", "#it{p}_{T} distribution - MC q<0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_neg", "#eta distribution - MC q<0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_neg", "#phi distribution - MC q<0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // primaries, secondaries
    histos.add("MC/pthist_tpc_prim", "#it{p}_{T} distribution - MC prim TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_prim", "#eta distribution - MC prim TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_prim", "#phi distribution - MC prim TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_prim", "#it{p}_{T} distribution - MC prim TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_prim", "#eta distribution - MC prim TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_prim", "#phi distribution - MC prim TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_secd", "#it{p}_{T} distribution - MC dec. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_secd", "#eta distribution - MC dec. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_secd", "#phi distribution - MC dec. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_secd", "#it{p}_{T} distribution - MC dec.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_secd", "#eta distribution - MC dec. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_secd", "#phi distribution - MC dec. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_secm", "#it{p}_{T} distribution - MC mat. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_secm", "#eta distribution - MC mat. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_secm", "#phi distribution - MC mat. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_secm", "#it{p}_{T} distribution - MC mat.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_secm", "#eta distribution - MC mat. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_secm", "#phi distribution - MC mat. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // pions only
    // all
    histos.add("MC/pthist_tpc_pi", "#it{p}_{T} distribution - #pi MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi", "#eta distribution - #pi MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi", "#phi distribution - #pi MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi", "#it{p}_{T} distribution - #pi MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi", "#eta distribution - #pi MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi", "#phi distribution - #pi MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // split in prim secd secm
    histos.add("MC/pthist_tpc_pi_prim", "#it{p}_{T} distribution - #pi MC prim TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_prim", "#eta distribution - #pi MC prim TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_prim", "#phi distribution - #pi MC prim TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_prim", "#it{p}_{T} distribution - #pi MC prim TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_prim", "#eta distribution - #pi MC prim TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_prim", "#phi distribution - #pi MC prim TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pi_secd", "#it{p}_{T} distribution - #pi MC dec. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_secd", "#eta distribution - #pi MC dec. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_secd", "#phi distribution - #pi MC dec. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_secd", "#it{p}_{T} distribution - #pi MC dec.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_secd", "#eta distribution - #pi MC dec. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_secd", "#phi distribution - #pi MC dec. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pi_secm", "#it{p}_{T} distribution - #pi MC mat. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_secm", "#eta distribution - #pi MC mat. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_secm", "#phi distribution - #pi MC mat. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_secm", "#it{p}_{T} distribution - #pi MC mat.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_secm", "#eta distribution - #pi MC mat. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_secm", "#phi distribution - #pi MC mat. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    // protons only
    // all
    histos.add("MC/pthist_tpc_P", "#it{p}_{T} distribution - prot MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_P", "#eta distribution - prot MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_P", "#phi distribution - prot MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_P", "#it{p}_{T} distribution - prot MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_P", "#eta distribution - prot MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_P", "#phi distribution - prot MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // kaons only
    // all
    histos.add("MC/pthist_tpc_K", "#it{p}_{T} distribution - kaons MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_K", "#eta distribution - kaons MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_K", "#phi distribution - kaons MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_K", "#it{p}_{T} distribution - kaons MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_K", "#eta distribution - kaons MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_K", "#phi distribution - kaons MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // pions+kaons
    // all
    histos.add("MC/pthist_tpc_piK", "#it{p}_{T} distribution - #pi+kaons MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_piK", "#eta distribution - #pi+kaons MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_piK", "#phi distribution - #pi+kaons MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_piK", "#it{p}_{T} distribution - #pi+kaons MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_piK", "#eta distribution - #pi+kaons MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_piK", "#phi distribution - #pi+kaons MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // pt>0.5 GeV/c threshold
    histos.add("MC/pthist_tpc_05", "#it{p}_{T} distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_05", "#eta distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_05", "#phi distribution - MC TPC tag, #it{p}_{T}>0.5", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_05", "#it{p}_{T} distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_05", "#eta distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_05", "#phi distribution - MC TPC+ITS tag, #it{p}_{T}>0.5", kTH1F, {axisPhi}, true);
    //
    // all but primary/secondary pions
    histos.add("MC/pthist_tpc_nopi", "#it{p}_{T} distribution - MC TPC tag ! prim/secd #pi", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_nopi", "#eta distribution - MC TPC tag ! prim/secd #pi", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_nopi", "#phi distribution - MC TPC tag ! prim/secd #pi", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_nopi", "#it{p}_{T} distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_nopi", "#eta distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_nopi", "#phi distribution - MC TPC+ITS tag ! prim/secd #pi", kTH1F, {axisPhi}, true);
    //
    // extras: difference between reconstructed and MC truth for eta, phi
    histos.add("MC/etahist_diff", "#eta difference track-MC ", kTH1F, {axisDEta}, true);
    histos.add("MC/phihist_diff", "#phi difference track-MC", kTH1F, {axisDPh}, true);
    //
    // hist sorting out PDG codes in wide bins
    histos.add("MC/pdghist_num", "PDG code - when non primary #pi TPC+ITS tag", kTH1F, {axisPDG}, true);
    histos.add("MC/pdghist_den", "PDG code - when non primary #pi TPC tag", kTH1F, {axisPDG}, true);

  } // end initMC

  /// Function applying the kinematic selections
  template <typename T>
  bool isTrackSelectedKineCuts(T& track)
  {
    if (!b_useTrackSelections)
      return true; // no track selections applied
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kPtRange))
      return false;
    if (!cutObject.IsSelected(track, TrackSelection::TrackCuts::kEtaRange))
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
  // //
  // // fill histos for TPC (all tracks)
  // void fillAllTPC(){
  //   histos.get<TH1>(HIST("MC/pthist_tpc"))->Fill(jT.pt());
  //   histos.get<TH1>(HIST("MC/phihist_tpc"))->Fill(jT.phi());
  //   histos.get<TH1>(HIST("MC/etahist_tpc"))->Fill(jT.eta());
  // }
  // //
  // // fill histos for TPC+ITS (all tracks)
  // void fillAllTPCITS(){
  //   histos.get<TH1>(HIST("MC/pthist_tpcits"))->Fill(jT.pt());
  //   histos.get<TH1>(HIST("MC/phihist_tpcits"))->Fill(jT.phi());
  //   histos.get<TH1>(HIST("MC/etahist_tpcits"))->Fill(jT.eta());
  // }
  //
  // define global variables
  int count = 0;
  int countData = 0;
  int countNoMC = 0;
  int tpPDGCode = 0;
  std::vector<int>::iterator itr_pdg;
  float pdg_fill = 0.0;
  //
  //////////////////////////////////////////////// PROCESS FUNCTIONS //////////////////////////////////////////////////
  //
  //
  void processMC(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& jTracks, aod::McParticles const& mcParticles)
  {
    //
    //
    for (auto& jT : jTracks) {
      // kinematic track seletions for all tracks
      if (!isTrackSelectedKineCuts(jT))
        continue;

      if (!jT.has_mcParticle()) {
        countNoMC++;
        if (doDebug)
          LOGF(warning, " N.%d track without MC particle, skipping...", countNoMC);
        continue;
      }
      auto mcpart = jT.mcParticle();
      tpPDGCode = TMath::Abs(mcpart.pdgCode());
      if (mcpart.isPhysicalPrimary()) {
        histos.get<TH1>(HIST("MC/etahist_diff"))->Fill(mcpart.eta() - jT.eta());
        auto delta = mcpart.phi() - jT.phi();
        if (delta > PI) {
          delta -= TwoPI;
        }
        if (delta < -PI) {
          delta += TwoPI;
        }
        histos.get<TH1>(HIST("MC/phihist_diff"))->Fill(delta);
      }
      // count the tracks contained in the input file if they have MC counterpart
      count++;
      //
      // all tracks, no conditions
      if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
        histos.get<TH1>(HIST("MC/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("MC/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("MC/etahist_tpc"))->Fill(jT.eta());
        if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpcits"))->Fill(jT.eta());
        } //  end if ITS
      }   //  end if TPC
      //
      // all tracks with pt>0.5
      if (jT.pt() > 0.5) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_05"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_05"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_05"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_05"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_05"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_05"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }     //  end if pt > 0.5
      //
      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_pos"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_pos"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_pos"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end positive
      //
      // negative only
      if (jT.signed1Pt() < 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_neg"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_neg"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_neg"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end negative
      //
      // only primaries
      if (mcpart.isPhysicalPrimary()) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_prim"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_prim"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_prim"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_prim"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }     //  end if primaries
      //
      // only secondaries from decay
      else if (mcpart.getProcess() == 4) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_secd"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_secd"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_secd"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_secd"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only secondaries from material
        else {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_secm"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_secm"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_secm"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_secm"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_secm"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_secm"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
        }     // end if secondaries from material
              //
      }       // end if secondaries from decay
      //
      // protons only
      if (tpPDGCode == 2212) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_P"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_P"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_P"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_P"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_P"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_P"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions only
      if (tpPDGCode == 211) {
        //
        // all tracks
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pi"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pi"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pi"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_pi"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_pi"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_pi"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only primary pions
        if (mcpart.isPhysicalPrimary()) {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_prim"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_prim"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_prim"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_prim"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
        }     //  end if primaries
        //
        // only secondary pions from decay
        else if (mcpart.getProcess() == 4) {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_secd"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_secd"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_secd"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_secd"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
          //
          // only secondary pions from material
          else {
            if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpc_pi_secm"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpc_pi_secm"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpc_pi_secm"))->Fill(jT.eta());
              if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
                histos.get<TH1>(HIST("MC/pthist_tpcits_pi_secm"))->Fill(jT.pt());
                histos.get<TH1>(HIST("MC/phihist_tpcits_pi_secm"))->Fill(jT.phi());
                histos.get<TH1>(HIST("MC/etahist_tpcits_pi_secm"))->Fill(jT.eta());
              } //  end if ITS
            }   //  end if TPC
          }     // end if secondaries from material
          //
        } // end if secondaries from decay
          //
      }   // end pions only
      //
      // no primary/sec-d pions
      if (!((tpPDGCode == 211) && (mcpart.isPhysicalPrimary()))) {
        // gets the pdg code and finds its index in our vector
        itr_pdg = std::find(pdgChoice.begin(), pdgChoice.end(), tpPDGCode);
        if (itr_pdg != pdgChoice.cend())
          // index from zero, so increase by 1 to put in the right bin (and 0.5 not needed but just not to sit in the edge)
          pdg_fill = (float)std::distance(pdgChoice.begin(), itr_pdg) + 1.5;
        else
          pdg_fill = -10.0;
        //
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_nopi"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_nopi"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_nopi"))->Fill(jT.eta());
          histos.get<TH1>(HIST("MC/pdghist_den"))->Fill(pdg_fill);
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_nopi"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_nopi"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_nopi"))->Fill(jT.eta());
            histos.get<TH1>(HIST("MC/pdghist_num"))->Fill(pdg_fill);
          } //  end if ITS
        }   //  end if TPC
      }     // end if not prim/sec-d pi
      //
      // kaons only
      if (tpPDGCode == 321) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_K"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_K"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_K"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_K"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_K"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_K"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions and kaons together
      if (tpPDGCode == 211 || tpPDGCode == 321) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_piK"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_piK"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_piK"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_piK"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_piK"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_piK"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      //
    } //  end loop on tracks
    //
    //
    if (doDebug)
      LOGF(info, "Tracks: %d, w/out MC: %d ", count, countNoMC);
  } // end processMC
  //
  PROCESS_SWITCH(qaMatchEff, processMC, "process MC", false);
  //
  //
  void processMCNoColl(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels> const& jTracks, aod::McParticles const& mcParticles)
  {
    //
    //
    for (auto& jT : jTracks) {
      // kinematic track seletions for all tracks
      if (!isTrackSelectedKineCuts(jT))
        continue;

      if (!jT.has_mcParticle()) {
        countNoMC++;
        if (doDebug)
          LOGF(warning, " N.%d track without MC particle, skipping...", countNoMC);
        continue;
      }
      auto mcpart = jT.mcParticle();
      tpPDGCode = TMath::Abs(mcpart.pdgCode());
      if (mcpart.isPhysicalPrimary()) {
        histos.get<TH1>(HIST("MC/etahist_diff"))->Fill(mcpart.eta() - jT.eta());
        auto delta = mcpart.phi() - jT.phi();
        if (delta > PI) {
          delta -= TwoPI;
        }
        if (delta < -PI) {
          delta += TwoPI;
        }
        histos.get<TH1>(HIST("MC/phihist_diff"))->Fill(delta);
      }
      // count the tracks contained in the input file if they have MC counterpart
      count++;
      //
      // all tracks, no conditions
      if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
        histos.get<TH1>(HIST("MC/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("MC/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("MC/etahist_tpc"))->Fill(jT.eta());
        if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpcits"))->Fill(jT.eta());
        } //  end if ITS
      }   //  end if TPC
      //
      // all tracks with pt>0.5
      if (jT.pt() > 0.5) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_05"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_05"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_05"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_05"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_05"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_05"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }     //  end if pt > 0.5
      //
      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_pos"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_pos"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_pos"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end positive
      //
      // negative only
      if (jT.signed1Pt() < 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_neg"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_neg"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_neg"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end negative
      //
      // only primaries
      if (mcpart.isPhysicalPrimary()) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_prim"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_prim"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_prim"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_prim"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }     //  end if primaries
      //
      // only secondaries from decay
      else if (mcpart.getProcess() == 4) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_secd"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_secd"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_secd"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_secd"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only secondaries from material
        else {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_secm"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_secm"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_secm"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_secm"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_secm"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_secm"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
        }     // end if secondaries from material
              //
      }       // end if secondaries from decay
      //
      // protons only
      if (tpPDGCode == 2212) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_P"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_P"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_P"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_P"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_P"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_P"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions only
      if (tpPDGCode == 211) {
        //
        // all tracks
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pi"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pi"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pi"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_pi"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_pi"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_pi"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only primary pions
        if (mcpart.isPhysicalPrimary()) {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_prim"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_prim"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_prim"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_prim"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
        }     //  end if primaries
        //
        // only secondary pions from decay
        else if (mcpart.getProcess() == 4) {
          if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_secd"))->Fill(jT.eta());
            if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_secd"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_secd"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_secd"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
          //
          // only secondary pions from material
          else {
            if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
              histos.get<TH1>(HIST("MC/pthist_tpc_pi_secm"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpc_pi_secm"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpc_pi_secm"))->Fill(jT.eta());
              if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
                histos.get<TH1>(HIST("MC/pthist_tpcits_pi_secm"))->Fill(jT.pt());
                histos.get<TH1>(HIST("MC/phihist_tpcits_pi_secm"))->Fill(jT.phi());
                histos.get<TH1>(HIST("MC/etahist_tpcits_pi_secm"))->Fill(jT.eta());
              } //  end if ITS
            }   //  end if TPC
          }     // end if secondaries from material
          //
        } // end if secondaries from decay
          //
      }   // end pions only
      //
      // no primary/sec-d pions
      if (!((tpPDGCode == 211) && (mcpart.isPhysicalPrimary()))) {
        // gets the pdg code and finds its index in our vector
        itr_pdg = std::find(pdgChoice.begin(), pdgChoice.end(), tpPDGCode);
        if (itr_pdg != pdgChoice.cend())
          // index from zero, so increase by 1 to put in the right bin (and 0.5 not needed but just not to sit in the edge)
          pdg_fill = (float)std::distance(pdgChoice.begin(), itr_pdg) + 1.5;
        else
          pdg_fill = -10.0;
        //
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_nopi"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_nopi"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_nopi"))->Fill(jT.eta());
          histos.get<TH1>(HIST("MC/pdghist_den"))->Fill(pdg_fill);
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_nopi"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_nopi"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_nopi"))->Fill(jT.eta());
            histos.get<TH1>(HIST("MC/pdghist_num"))->Fill(pdg_fill);
          } //  end if ITS
        }   //  end if TPC
      }     // end if not prim/sec-d pi
      //
      // kaons only
      if (tpPDGCode == 321) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_K"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_K"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_K"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_K"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_K"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_K"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions and kaons together
      if (tpPDGCode == 211 || tpPDGCode == 321) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("MC/pthist_tpc_piK"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_piK"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_piK"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_piK"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_piK"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_piK"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      //
    } //  end loop on tracks
    //
    //
    if (doDebug)
      LOGF(info, "Tracks: %d, w/out MC: %d ", count, countNoMC);
  } // end processMCNoColl
  //
  PROCESS_SWITCH(qaMatchEff, processMCNoColl, "process MC - no loop on collisions", false);
  //
  //
  void processData(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& jTracks)
  {
    //
    //
    for (auto& jT : jTracks) {
      // kinematic track seletions for all tracks
      if (!isTrackSelectedKineCuts(jT))
        continue;
      //
      countData++;
      //
      // all tracks, no conditions
      if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
        histos.get<TH1>(HIST("data/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("data/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("data/etahist_tpc"))->Fill(jT.eta());
        if (jT.pt() > 0.5) {
          histos.get<TH1>(HIST("data/pthist_tpc_05"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_05"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_05"))->Fill(jT.eta());
        }
        if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpcits"))->Fill(jT.eta());
          if (jT.pt() > 0.5) {
            histos.get<TH1>(HIST("data/pthist_tpcits_05"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_05"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_05"))->Fill(jT.eta());
          }
        } //  end if ITS
      }   //  end if TPC
      //

      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("data/pthist_tpcits_pos"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_pos"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_pos"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end positive
      //
      // negative only
      if (jT.signed1Pt() < 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("data/pthist_tpcits_neg"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_neg"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_neg"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end negative
      //
      //
    } //  end loop on tracks
    //
    //
    if (doDebug)
      LOGF(info, "Tracks: %d ", countData);
    //
  } // end processData
  //
  PROCESS_SWITCH(qaMatchEff, processData, "process data", true);
  //
  //
  void processDataNoColl(soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& jTracks)
  {
    //
    //
    for (auto& jT : jTracks) {
      // kinematic track seletions for all tracks
      if (!isTrackSelectedKineCuts(jT))
        continue;
      //
      countData++;
      //
      // all tracks, no conditions
      if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
        histos.get<TH1>(HIST("data/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("data/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("data/etahist_tpc"))->Fill(jT.eta());
        if (jT.pt() > 0.5) {
          histos.get<TH1>(HIST("data/pthist_tpc_05"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_05"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_05"))->Fill(jT.eta());
        }
        if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpcits"))->Fill(jT.eta());
          if (jT.pt() > 0.5) {
            histos.get<TH1>(HIST("data/pthist_tpcits_05"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_05"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_05"))->Fill(jT.eta());
          }
        } //  end if ITS
      }   //  end if TPC
      //

      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("data/pthist_tpcits_pos"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_pos"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_pos"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end positive
      //
      // negative only
      if (jT.signed1Pt() < 0) {
        if (jT.hasTPC() && isTrackSelectedTPCCuts(jT)) {
          histos.get<TH1>(HIST("data/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS() && isTrackSelectedITSCuts(jT)) {
            histos.get<TH1>(HIST("data/pthist_tpcits_neg"))->Fill(jT.pt());
            histos.get<TH1>(HIST("data/phihist_tpcits_neg"))->Fill(jT.phi());
            histos.get<TH1>(HIST("data/etahist_tpcits_neg"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
            //
      }     // end negative
      //
      //
    } //  end loop on tracks
    //
    //
    if (doDebug)
      LOGF(info, "Tracks: %d ", countData);
    //
  } // end processDataNoColl
  //
  PROCESS_SWITCH(qaMatchEff, processDataNoColl, "process data - no collision dependence", true);
  //
  //
}; // end of structure

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaMatchEff>(cfgc, TaskName{"qa-match-eff"})};
}
