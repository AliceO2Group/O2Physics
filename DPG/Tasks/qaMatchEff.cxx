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
/// \brief A task for basic checks on track matching efficiency
/// \author Rosario Turrisi (rosario.turrisi@pd.infn.it)
/// \since 1887

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CommonConstants/MathConstants.h"
#include "Common/DataModel/EventSelection.h"
//
// centrality
//  ????

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
  //
  //
  Configurable<bool> isitMC{"isitMC", false, "Reading MC files, data if false"};
  Configurable<bool> doDebug{"doDebug", false, "Flag of debug information"};
  // Histogram configuration
  //
  // histo x axes limits
  Configurable<float> etaMin{"eta-min", -0.8f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 0.8f, "Upper limit in eta"};
  Configurable<float> phiMin{"phi-min", 0.0f, "Lower limit in phi"};
  Configurable<float> phiMax{"phi-max", 1.0f * TwoPI, "Upper limit in phi"};
  Configurable<float> ptMin{"pt-min", 0.0f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 15.f, "Upper limit in pT"};
  // histos bins
  Configurable<int> ptBins{"pt-bins", 16, "Number of pT bins"};
  Configurable<int> etaBins{"eta-bins", 16, "Number of eta bins"};
  Configurable<int> phiBins{"phi-bins", 18, "Number of phi bins"};
  // histo axes
  //
  // non uniform pt binning
  std::vector<double> ptBinning = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0};
  AxisSpec axisPt{ptBinning, "#it{p}_{T} (GeV/#it{c})"};
  // const AxisSpec axisPt{ptBins, ptMin, ptMax, "#it{p}_{T} (GeV/#it{c})"};
  //
  AxisSpec axisEta{etaBins, etaMin, etaMax, "#it{#eta}"};
  AxisSpec axisPhi{phiBins, phiMin, phiMax, "#it{#varphi} (rad)"};
  AxisSpec axisDEta{etaBins, etaMin, etaMax, "D#it{#eta}"};
  AxisSpec axisDPh{phiBins, -PI, PI, "D#it{#varphi} (rad)"};
  //
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
    histos.add("data/etahist_tpc", "#it{#eta} distribution - data TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc", "#it{#phi} distribution - data TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits", "#it{p}_{T} distribution - data TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits", "#it{#eta} distribution - data TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits", "#it{#phi} distribution - data TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpc_pos", "#it{p}_{T} distribution - data q>0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc_pos", "#it{#eta} distribution - data q>0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc_pos", "#it{#phi} distribution - data q>0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits_pos", "#it{p}_{T} distribution - data q>0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits_pos", "#it{#eta} distribution - data q>0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits_pos", "#it{#phi} distribution - data q>0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpc_neg", "#it{p}_{T} distribution - data q<0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpc_neg", "#it{#eta} distribution - data q<0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpc_neg", "#it{#phi} distribution - data q<0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("data/pthist_tpcits_neg", "#it{p}_{T} distribution - data q<0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("data/etahist_tpcits_neg", "#it{#eta} distribution - data q<0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("data/phihist_tpcits_neg", "#it{#phi} distribution - data q<0 TPC+ITS tag", kTH1F, {axisPhi}, true);
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
    histos.add("MC/etahist_tpc", "#it{#eta} distribution - MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc", "#it{#phi} distribution - MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits", "#it{p}_{T} distribution - MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits", "#it{#eta} distribution - MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits", "#it{#phi} distribution - MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pos", "#it{p}_{T} distribution - MC q>0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pos", "#it{#eta} distribution - MC q>0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pos", "#it{#phi} distribution - MC q>0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pos", "#it{p}_{T} distribution - MC q>0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pos", "#it{#eta} distribution - MC q>0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pos", "#it{#phi} distribution - MC q>0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_neg", "#it{p}_{T} distribution - MC q<0 TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_neg", "#it{#eta} distribution - MC q<0 TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_neg", "#it{#phi} distribution - MC q<0 TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_neg", "#it{p}_{T} distribution - MC q<0 TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_neg", "#it{#eta} distribution - MC q<0 TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_neg", "#it{#phi} distribution - MC q<0 TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // primaries, secondaries
    histos.add("MC/pthist_tpc_prim", "#it{p}_{T} distribution - MC prim TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_prim", "#it{#eta} distribution - MC prim TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_prim", "#it{#phi} distribution - MC prim TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_prim", "#it{p}_{T} distribution - MC prim TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_prim", "#it{#eta} distribution - MC prim TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_prim", "#it{#phi} distribution - MC prim TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_secd", "#it{p}_{T} distribution - MC dec. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_secd", "#it{#eta} distribution - MC dec. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_secd", "#it{#phi} distribution - MC dec. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_secd", "#it{p}_{T} distribution - MC dec.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_secd", "#it{#eta} distribution - MC dec. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_secd", "#it{#phi} distribution - MC dec. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_secm", "#it{p}_{T} distribution - MC mat. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_secm", "#it{#eta} distribution - MC mat. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_secm", "#it{#phi} distribution - MC mat. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_secm", "#it{p}_{T} distribution - MC mat.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_secm", "#it{#eta} distribution - MC mat. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_secm", "#it{#phi} distribution - MC mat. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // pions only
    // all
    histos.add("MC/pthist_tpc_pi", "#it{p}_{T} distribution - {#pi} MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi", "#it{#eta} distribution - {#pi} MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi", "#it{#phi} distribution - {#pi} MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi", "#it{p}_{T} distribution - {#pi} MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi", "#it{#eta} distribution - {#pi} MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi", "#it{#phi} distribution - {#pi} MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // split in prim secd secm
    histos.add("MC/pthist_tpc_pi_prim", "#it{p}_{T} distribution - {#pi} MC prim TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_prim", "#it{#eta} distribution - {#pi} MC prim TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_prim", "#it{#phi} distribution - {#pi} MC prim TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_prim", "#it{p}_{T} distribution - {#pi} MC prim TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_prim", "#it{#eta} distribution - {#pi} MC prim TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_prim", "#it{#phi} distribution - {#pi} MC prim TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pi_secd", "#it{p}_{T} distribution - {#pi} MC dec. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_secd", "#it{#eta} distribution - {#pi} MC dec. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_secd", "#it{#phi} distribution - {#pi} MC dec. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_secd", "#it{p}_{T} distribution - {#pi} MC dec.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_secd", "#it{#eta} distribution - {#pi} MC dec. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_secd", "#it{#phi} distribution - {#pi} MC dec. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpc_pi_secm", "#it{p}_{T} distribution - {#pi} MC mat. sec. TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_pi_secm", "#it{#eta} distribution - {#pi} MC mat. sec. TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_pi_secm", "#it{#phi} distribution - {#pi} MC mat. sec. TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_pi_secm", "#it{p}_{T} distribution - {#pi} MC mat.sec. TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_pi_secm", "#it{#eta} distribution - {#pi} MC mat. sec. TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_pi_secm", "#it{#phi} distribution - {#pi} MC mat. sec. TPC+ITS tag", kTH1F, {axisPhi}, true);
    // protons only
    // all
    histos.add("MC/pthist_tpc_P", "#it{p}_{T} distribution - prot MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_P", "#it{#eta} distribution - prot MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_P", "#it{#phi} distribution - prot MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_P", "#it{p}_{T} distribution - prot MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_P", "#it{#eta} distribution - prot MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_P", "#it{#phi} distribution - prot MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // kaons only
    // all
    histos.add("MC/pthist_tpc_K", "#it{p}_{T} distribution - kaons MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_K", "#it{#eta} distribution - kaons MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_K", "#it{#phi} distribution - kaons MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_K", "#it{p}_{T} distribution - kaons MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_K", "#it{#eta} distribution - kaons MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_K", "#it{#phi} distribution - kaons MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    // pions+kaons
    // all
    histos.add("MC/pthist_tpc_piK", "#it{p}_{T} distribution - {#pi}+kaons MC TPC tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpc_piK", "#it{#eta} distribution - {#pi}+kaons MC TPC tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpc_piK", "#it{#phi} distribution - {#pi}+kaons MC TPC tag", kTH1F, {axisPhi}, true);
    histos.add("MC/pthist_tpcits_piK", "#it{p}_{T} distribution - {#pi}+kaons MC TPC+ITS tag", kTH1F, {axisPt}, true);
    histos.add("MC/etahist_tpcits_piK", "#it{#eta} distribution - {#pi}+kaons MC TPC+ITS tag", kTH1F, {axisEta}, true);
    histos.add("MC/phihist_tpcits_piK", "#it{#phi} distribution - {#pi}+kaons MC TPC+ITS tag", kTH1F, {axisPhi}, true);
    //
    // extras: difference between reconstructed and MC truth for eta, phi
    histos.add("MC/etahist_diff", "#it{#eta} difference track-MC ", kTH1F, {axisDEta}, true);
    histos.add("MC/phihist_diff", "#it{#phi} difference track-MC", kTH1F, {axisDPh}, true);
  }
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
  //
  //////////////////////////////////////////////// PROCESS FUNCTIONS //////////////////////////////////////////////////
  //
  //
  void processMC(soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels> const& jTracks, aod::McParticles const& mcParticles)
  {
    //
    //
    for (auto& jT : jTracks) {
      if (!jT.has_mcParticle()) {
        countNoMC++;
        if (doDebug)
          LOGF(warning, " N.%d track without MC particle, skipping...", countNoMC);
        continue;
      }
      auto mcpart = jT.mcParticle();
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
      if (jT.hasTPC()) {
        histos.get<TH1>(HIST("MC/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("MC/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("MC/etahist_tpc"))->Fill(jT.eta());
        if (jT.hasITS()) {
          histos.get<TH1>(HIST("MC/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpcits"))->Fill(jT.eta());
        } //  end if ITS
      }   //  end if TPC
      //
      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS()) {
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
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS()) {
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
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_prim"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_prim"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_prim"))->Fill(jT.eta());
          if (jT.hasITS()) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_prim"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }     //  end if primaries
      //
      // only secondaries from decay
      else if (mcpart.getProcess() == 4) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_secd"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_secd"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_secd"))->Fill(jT.eta());
          if (jT.hasITS()) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_secd"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only secondaries from material
        else {
          if (jT.hasTPC()) {
            histos.get<TH1>(HIST("MC/pthist_tpc_secm"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_secm"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_secm"))->Fill(jT.eta());
            if (jT.hasITS()) {
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
      if (TMath::Abs(mcpart.pdgCode()) == 2212) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_P"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_P"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_P"))->Fill(jT.eta());
          if (jT.hasITS()) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_P"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_P"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_P"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions only
      if (TMath::Abs(mcpart.pdgCode()) == 211) {
        //
        // all tracks
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_pi"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_pi"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_pi"))->Fill(jT.eta());
          if (jT.hasITS()) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_pi"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_pi"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_pi"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
        //
        // only primary pions
        if (mcpart.isPhysicalPrimary()) {
          if (jT.hasTPC()) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_prim"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_prim"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_prim"))->Fill(jT.eta());
            if (jT.hasITS()) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_prim"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_prim"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_prim"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
        }     //  end if primaries
        //
        // only secondary pions from decay
        else if (mcpart.getProcess() == 4) {
          if (jT.hasTPC()) {
            histos.get<TH1>(HIST("MC/pthist_tpc_pi_secd"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpc_pi_secd"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpc_pi_secd"))->Fill(jT.eta());
            if (jT.hasITS()) {
              histos.get<TH1>(HIST("MC/pthist_tpcits_pi_secd"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpcits_pi_secd"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpcits_pi_secd"))->Fill(jT.eta());
            } //  end if ITS
          }   //  end if TPC
          //
          // only secondary pions from material
          else {
            if (jT.hasTPC()) {
              histos.get<TH1>(HIST("MC/pthist_tpc_pi_secm"))->Fill(jT.pt());
              histos.get<TH1>(HIST("MC/phihist_tpc_pi_secm"))->Fill(jT.phi());
              histos.get<TH1>(HIST("MC/etahist_tpc_pi_secm"))->Fill(jT.eta());
              if (jT.hasITS()) {
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
      // kaons only
      if (TMath::Abs(mcpart.pdgCode()) == 321) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_K"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_K"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_K"))->Fill(jT.eta());
          if (jT.hasITS()) {
            histos.get<TH1>(HIST("MC/pthist_tpcits_K"))->Fill(jT.pt());
            histos.get<TH1>(HIST("MC/phihist_tpcits_K"))->Fill(jT.phi());
            histos.get<TH1>(HIST("MC/etahist_tpcits_K"))->Fill(jT.eta());
          } //  end if ITS
        }   //  end if TPC
      }
      //
      // pions and kaons together
      if (TMath::Abs(mcpart.pdgCode()) == 211 || TMath::Abs(mcpart.pdgCode()) == 321) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("MC/pthist_tpc_piK"))->Fill(jT.pt());
          histos.get<TH1>(HIST("MC/phihist_tpc_piK"))->Fill(jT.phi());
          histos.get<TH1>(HIST("MC/etahist_tpc_piK"))->Fill(jT.eta());
          if (jT.hasITS()) {
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
  void processData(soa::Join<aod::Tracks, aod::TracksExtra> const& jTracks)
  {
    //
    //
    for (auto& jT : jTracks) {
      //
      countData++;
      //
      // all tracks, no conditions
      if (jT.hasTPC()) {
        histos.get<TH1>(HIST("data/pthist_tpc"))->Fill(jT.pt());
        histos.get<TH1>(HIST("data/phihist_tpc"))->Fill(jT.phi());
        histos.get<TH1>(HIST("data/etahist_tpc"))->Fill(jT.eta());
        if (jT.hasITS()) {
          histos.get<TH1>(HIST("data/pthist_tpcits"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpcits"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpcits"))->Fill(jT.eta());
        } //  end if ITS
      }   //  end if TPC
      //

      // positive only
      if (jT.signed1Pt() > 0) {
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("data/pthist_tpc_pos"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_pos"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_pos"))->Fill(jT.eta());
          if (jT.hasITS()) {
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
        if (jT.hasTPC()) {
          histos.get<TH1>(HIST("data/pthist_tpc_neg"))->Fill(jT.pt());
          histos.get<TH1>(HIST("data/phihist_tpc_neg"))->Fill(jT.phi());
          histos.get<TH1>(HIST("data/etahist_tpc_neg"))->Fill(jT.eta());
          if (jT.hasITS()) {
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
}; // end of structure

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaMatchEff>(cfgc, TaskName{"qa-match-eff"})};
}
