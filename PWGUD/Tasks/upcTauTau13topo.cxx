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
// \brief tau tau analysis 1e+3pi topology
// \author Adam Matyja, adam.tomasz.matyja@cern.ch, adam.matryja@ifj.edu.pl
// \since  January 2024
// to run it execute:
// copts="--configuration json://tautauConfig.json -b"
// o2-analysis-ud-tautau13topo $copts > output.log


#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"
//#include "Common/DataModel/EventSelection.h"
//#include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/DGPIDSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TauTau13topo {

  // configurables
  ConfigurableAxis ptAxis{"pAxis", {100, 0., 5.}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis etaAxis{"etaAxis", {100, -2., 2.}, "#eta"};
  ConfigurableAxis dedxAxis{"dedxAxis", {100, 20., 160.}, "dE/dx"};
  ConfigurableAxis minvAxis{"MinvAxis", {100, 0., 2.5}, "M_{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis phiAxis{"phiAxis", {100, 0., 3.2}, "#phi"};
  ConfigurableAxis vectorAxis{"vectorAxis", {100, 0., 2.}, "A_{V}"};
  ConfigurableAxis scalarAxis{"scalarAxis", {100, -1., 1.}, "A_{S}"};

  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};
  //ConfigurableAxis IVMAxis{"IVMAxis", {350, 0.0, 3.5}, "Invariant mass axis"};
  //ConfigurableAxis ptAxis{"ptAxis", {250, 0.0, 2.5}, "p_T axis"};
  //ConfigurableAxis nsTPCAxis{"nsTPCAxis", {100, -20.0, 20.0}, "nSigma TPC axis"};
  //ConfigurableAxis nsTOFAxis{"nsTOFAxis", {100, -100.0, 100.0}, "nSigma TOF axis"};

  //cut selection configurables
  Configurable<float> zvertexcut{"Zvertexcut", 15, "Z vertex cut"};
  Configurable<float> trkEtacut{"TrkEtacut", 15, "max track eta cut"};
  Configurable<float> minAnglecut{"minAnglecut", 0.05, "min angle between tracks cut"};
  Configurable<float> minNsigmaElcut{"minNsigmaElcut", -2, "min Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaElcut{"maxNsigmaElcut",  3, "max Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaPiVetocut{"maxNsigmaPiVetocut",  4, "max Nsigma for Pion veto cut"};
  Configurable<float> minPtEtrkcut{"minPtEtrkcut", 0.25, "min Pt for El track cut"};
  Configurable<bool> FITvetoFlag{"FITvetoFlag", {}, "To apply FIT veto"};
  Configurable<int> FITvetoWindow{"FITvetoWindow", 1, "FIT veto window"};
  // a pdg object
  TDatabasePDG* pdg = nullptr;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    pdg = TDatabasePDG::Instance();

    // dgcandidates histograms
    //    const AxisSpec axisIVM{IVMAxis, "IVM axis"};
    const AxisSpec axispt{ptAxis, "p_{T} axis"};
    const AxisSpec axiseta{etaAxis, "#eta - pseudo rapidity axis"};
    const AxisSpec axisdedx{dedxAxis, "dEdx axis"};
    const AxisSpec axisminv{minvAxis, "invariant mass axis"};
    const AxisSpec axisphi{phiAxis, "phi axis"};
    const AxisSpec axisav{vectorAxis, "AV axis"};
    const AxisSpec axisas{scalarAxis, "AS axis"};

    registry.add("global/hVertexXY", "Vertex position in x and y direction; V_x; V_y; Collisions", {HistType::kTH2F, {{50, -0.05, 0.05}, {50, -0.05, 0.05}}});
    registry.add("global/hVertexZ",  "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{ 100, -25., 25.}}});
    registry.add("global/hVertexZ15", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{ 100, -25., 25.}}});
    registry.add("global/hVertexZ10", "Vertex position in z direction; V_z; Collisions", {HistType::kTH1F, {{ 100, -25., 25.}}});
    registry.add("global/hNTracks", ";N_{tracks};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global/hNTracksGlobal", ";N_{tracks,global};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global/hNTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{ 5, 0., 5.}}});
    registry.add("global/hTrackPtPV", ";p_T^{trk}; Entries", {HistType::kTH1F, {axispt}});
    registry.add("global/hTrackPVTotCharge", "Q_{Tot};Q_{Tot}; Entries", {HistType::kTH1F, {{10,-5,5}}});
    registry.add("global/hTrackEtaPhiPV", ";Eta;Phi;",{HistType::kTH2D,{axiseta,{140, -3.5, 3.5}}});
    registry.add("global/hSignalTPCvsPtPV", ";Pt;TPC Signal",{HistType::kTH2F,{axispt,{200, 0., 200}}});
    registry.add("global/hITSbitPVtrk", "ITS bit for PV tracks; Layer hit;Entries", {HistType::kTH1F, {{ 10, 0., 10.}}});
    registry.add("global/hITSnbitsVsEtaPVtrk", "n ITS bits vs #eta for PV tracks; #eta;Layer hit;Entries", {HistType::kTH2F, {axiseta,{ 8, -1.,7.}}});
    registry.add("global/hITSbitVsEtaPVtrk", "ITS bit vs #eta for PV tracks; #eta;Layer hit;Entries", {HistType::kTH2F, {axiseta,{ 8, 0.,8.}}});
    registry.add("global/hEventEff", "Event cut efficiency: 0-All,1-PV=4,2-Qtot=0,3-El;Cut;entries", {HistType::kTH1F, {{ 25, 0., 25.} }});
    registry.add("global/hNCombAfterCut", "Combinations after cut: 0-All,5-M3pi,10-Dphi,15-N_{e},20-N_{v#pi},25-Pt,30-Vcal,35-Tot;N_{comb};entries", {HistType::kTH1F, {{ 40, 0., 40.} }});
    //registry.add("global/hInvMassElTrack", ";M_{inv}^{2};entries", {HistType::kTH1F, {{100, -0.01, 0.49 } }});
    registry.add("global/hDeltaAngleTrackPV", ";#Delta#alpha;entries", {HistType::kTH1F, {{100, -0.01, 0.49 } }});

    //cut0
    registry.add("control/h0Cut3piMassComb", "3#pi mass, 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis }});
    registry.add("control/h0Cut3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis }});
    registry.add("control/h0CutDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, { phiAxis}});
    registry.add("control/h0Cut13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.} }});
    registry.add("control/h0Cut13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis }});
    registry.add("control/h0Cut13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis }});
    registry.add("control/h0Cut13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.} }});
    registry.add("control/h0Cut4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt }});
    registry.add("control/h0Cut4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.} }});
    registry.add("control/h0Cut3piMassVsPt", "3#pi mass vs Pt, 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis,axispt }});
    registry.add("control/h0Cut4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100,1,5.},axispt }});
    //cut1
    registry.add("control/h1Cut3piMassComb", "3#pi mass, 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis }});
    registry.add("control/h1Cut3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis }});
    registry.add("control/h1CutDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, { phiAxis}});
    registry.add("control/h1Cut13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.} }});
    registry.add("control/h1Cut13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis }});
    registry.add("control/h1Cut13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis }});
    registry.add("control/h1Cut13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.} }});
    registry.add("control/h1Cut4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt }});
    registry.add("control/h1Cut4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.} }});
    registry.add("control/h1Cut3piMassVsPt", "3#pi mass vs Pt, 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis,axispt }});
    registry.add("control/h1Cut4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100,1,5.},axispt }});
    registry.add("control/h1CutDcaZ", "All 4 tracks dca ;dca_{Z};entries", {HistType::kTH1F, {{100,-0.05,0.05}}});
    registry.add("control/h1CutDcaXY", "All 4 tracks dca ;dca_{XY};entries", {HistType::kTH1F, {{100,-0.05,0.05} }});
    registry.add("control/h1CutChi2TPC", "All 4 tracks Chi2 ;Chi2_{TPC};entries", {HistType::kTH1F, {{48,-2,10.}}});
    registry.add("control/h1CutChi2ITS", "All 4 tracks Chi2 ;Chi2_{ITS};entries", {HistType::kTH1F, {{44,-2,20.}}});
    registry.add("control/h1CutTPCnclsFindable", "All 4 tracks NclFind ;N_{TPC,cl,findable};entries", {HistType::kTH1F, {{160,0,160.}}});

    //pid
    registry.add("pidTPC/hpvsdedxElHipCut0","In hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut1","All hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    //pid separately for each cut
    registry.add("pidTPC/hpvsdedxElHipCut2","IM hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut3","DP hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut4","El hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut5","Pi hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut6","Pt hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut7","Vc hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    //pid sequentialy 
    registry.add("pidTPC/hpvsdedxElHipCut10","El hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut11","Pi+10 hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut12","Vc+11 hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});
    registry.add("pidTPC/hpvsdedxElHipCut13","Pt+12 hip;p_{trk};dEdx_{trk}", {HistType::kTH2F, {ptAxis,dedxAxis }});


    registry.add("global/hFinalPtSpectrumEl",  ";p_{T}^{e} (GeV/c);entries",  {HistType::kTH1F, {{40, 0.,5.} }});

    registry.add("fit/bbFT0Abit",  "FT0A bits;bit;entries",  {HistType::kTH1F, {{32, 0.,32.} }});
    registry.add("fit/bbFT0Cbit",  "FT0C bits;bit;entries",  {HistType::kTH1F, {{32, 0.,32.} }});
    registry.add("fit/bbFV0Abit",  "FV0A bits;bit;entries",  {HistType::kTH1F, {{32, 0.,32.} }});
    registry.add("fit/bbFDDAbit",  "FDDA bits;bit;entries",  {HistType::kTH1F, {{32, 0.,32.} }});
    registry.add("fit/bbFDDCbit",  "FDDC bits;bit;entries",  {HistType::kTH1F, {{32, 0.,32.} }});
    registry.add("fit/bbFT0Aamplitude",  "FT0A amplitude;Amplitude;entries",  {HistType::kTH1F, {{100, -5.,95.} }});
    registry.add("fit/bbFT0Camplitude",  "FT0C amplitude;Amplitude;entries",  {HistType::kTH1F, {{100, -5.,95.} }});
    registry.add("fit/bbFV0Aamplitude",  "FV0A amplitude;Amplitude;entries",  {HistType::kTH1F, {{100, -5.,95.} }});
    registry.add("fit/bbFDDAamplitude",  "FDDA amplitude;Amplitude;entries",  {HistType::kTH1F, {{100, -5.,95.} }});
    registry.add("fit/bbFDDCamplitude",  "FDDC amplitude;Amplitude;entries",  {HistType::kTH1F, {{100, -5.,95.} }});

    registry.add("fit/timeFT0",  "FT0 time;time FT0A; time FT0C;entries",  {HistType::kTH2F, {{100, -5.,35.} ,{100, -5.,35.}}});
    registry.add("fit/timeFDD",  "FDD time;time FDDA; time FDDC;entries",  {HistType::kTH2F, {{100, -5.,35.} ,{100, -5.,35.}}});
  }

  float CalculateDeltaPhi(TLorentzVector p, TLorentzVector p1){
    float delta=p.Phi();
    if(delta < 0) delta += TMath::TwoPi();
    if(p1.Phi()<0) delta -= (p1.Phi()+TMath::TwoPi());
    else delta -= p1.Phi();
    if(delta<0) delta+=TMath::TwoPi();
    if(delta>TMath::Pi()) delta = TMath::TwoPi() - delta;
    return delta;
  }

  //fill control histograms per track
  template <int mode, typename T>  
  void FillControlHistos(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float pi3etasum){
    static constexpr std::string_view histoname[] = {"0","1","2","3","4","5","6","7","8","9",
						     "10","11","12","13","14","15","16","17","18","19",
						     "20","21","22","23","24","25","26","27","28","29",
						     "30","31","32","33","34","35","36","37","38","39"
    };
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut3piMassComb"))        ->Fill( pi3invMass);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut3trkPtTot"))          ->Fill(      pi3pt);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("CutDeltaPhi13topo"))     ->Fill(pi3deltaPhi);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut13AssymPt1ProngAver"))->Fill( pi3assymav);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut13Vector"))           ->Fill(  pi3vector);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut13Scalar"))           ->Fill(  pi3scalar);
    registry.get<TH1>(HIST("control/h") + HIST(histoname[mode]) + HIST("Cut13EtaSum"))           ->Fill(  pi3etasum);
  }



  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;

  void process(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  {
    //global checks
    registry.get<TH2>(HIST("global/hVertexXY"))->Fill(dgcand.posX(), dgcand.posY());
    registry.get<TH1>(HIST("global/hVertexZ"))->Fill(dgcand.posZ());
    if(TMath::Abs(dgcand.posZ())<15) registry.get<TH1>(HIST("global/hVertexZ15"))->Fill(dgcand.posZ());
    if(TMath::Abs(dgcand.posZ())<10) registry.get<TH1>(HIST("global/hVertexZ10"))->Fill(dgcand.posZ());

    registry.get<TH1>(HIST("global/hNTracks"))->Fill(dgtracks.size());

    //setup PV tracks partition
    Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    PVContributors.bindTable(dgtracks);

    registry.get<TH1>(HIST("global/hNTracksPV"))->Fill(PVContributors.size());

    UChar_t clustermap1;
    int nTofTrk=0;
    int nEtaIn15=0;
    int nITSbits=0;
    TLorentzVector p;
    TParticlePDG* pion = pdg->GetParticle(211);
    //loop over PV contributors
    for (auto trk : PVContributors) {
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), pion->Mass());
      registry.get<TH1>(HIST("global/hTrackPtPV"))->Fill(p.Pt());
      if(TMath::Abs(p.Eta()) < trkEtacut) nEtaIn15++;//1.5 is a default
      registry.get<TH2>(HIST("global/hTrackEtaPhiPV"))->Fill(p.Eta(),p.Phi());
      nITSbits=-1;
      if(trk.hasITS()){//ITS track
	clustermap1 = trk.itsClusterMap();
	for(int bitNo=0;bitNo<7;bitNo++){
	  if(TESTBIT(clustermap1, bitNo)) {//check ITS bits/layers for each PV track
	    registry.get<TH1>(HIST("global/hITSbitPVtrk"))->Fill(bitNo,1.);
	    registry.get<TH2>(HIST("global/hITSbitVsEtaPVtrk"))->Fill(p.Eta(),bitNo,1.);
	    nITSbits++;
	  }
	}//end of loop over ITS bits
      }//has ITS
      registry.get<TH2>(HIST("global/hITSnbitsVsEtaPVtrk"))->Fill(p.Eta(), nITSbits);
      if(trk.hasTPC())
	registry.get<TH2>(HIST("global/hSignalTPCvsPtPV"))->Fill(p.Pt(), trk.tpcSignal());
      if (trk.hasTOF()) nTofTrk++;
    }
    registry.get<TH1>(HIST("global/hNtofTrk"))->Fill(nTofTrk);

    //
    //selection
    //
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(0.,1.);

    // skip events with too few/many tracks
    if (dgcand.numContrib() != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Number of PV contributors is %d", dgcand.numContrib());
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(1.,1.);
    registry.get<TH1>(HIST("global/hTrackPVTotCharge"))->Fill(dgcand.netCharge());

    //if vz<15
    if(TMath::Abs(dgcand.posZ()) >= zvertexcut ){//default =15
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: VertexZ is %d", dgcand.posZ());
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(2.,1.);

    //if eta tracks <1.5
    if(nEtaIn15 != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Ntrk inside |eta|<1.5 is %d", nEtaIn15);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(3.,1.);

    // skip events with net charge != 0
    if (dgcand.netCharge() != 0) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Net charge is %d", dgcand.netCharge());
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(4.,1.);


    //    // skip events with out-of-range rgtrwTOF (fraction-of-good-tracks-with-TOF-hit)
    //    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    //    if (rtrwTOF < 0.25) {
    //      if (verbose) {
    //        LOGF(debug, "<tautau13topo> Candidate rejected: rtrwTOF is %f", rtrwTOF);
    //      }
    //      return;
    //    }
    
    //FIT informaton
    for (auto bit = 0; bit <= 32; bit++) {
      registry.get<TH1>(HIST("fit/bbFT0Abit"))->Fill(bit,TESTBIT(dgcand.bbFT0Apf(),bit) ) ;
      registry.get<TH1>(HIST("fit/bbFT0Cbit"))->Fill(bit,TESTBIT(dgcand.bbFT0Cpf(),bit) ) ;
      registry.get<TH1>(HIST("fit/bbFV0Abit"))->Fill(bit,TESTBIT(dgcand.bbFV0Apf(),bit) ) ;
      registry.get<TH1>(HIST("fit/bbFDDAbit"))->Fill(bit,TESTBIT(dgcand.bbFDDApf(),bit) ) ;
      registry.get<TH1>(HIST("fit/bbFDDCbit"))->Fill(bit,TESTBIT(dgcand.bbFDDCpf(),bit) ) ;
    }
    registry.get<TH1>(HIST("fit/bbFT0Aamplitude"))->Fill(dgcand.totalFT0AmplitudeA()) ;
    registry.get<TH1>(HIST("fit/bbFT0Camplitude"))->Fill(dgcand.totalFT0AmplitudeC()) ;
    registry.get<TH1>(HIST("fit/bbFV0Aamplitude"))->Fill(dgcand.totalFV0AmplitudeA()) ;
    registry.get<TH1>(HIST("fit/bbFDDAamplitude"))->Fill(dgcand.totalFDDAmplitudeA()) ;
    registry.get<TH1>(HIST("fit/bbFDDCamplitude"))->Fill(dgcand.totalFDDAmplitudeC()) ;
    
    registry.get<TH2>(HIST("fit/timeFT0"))->Fill(dgcand.timeFT0A(),dgcand.timeFT0C()) ;
    registry.get<TH2>(HIST("fit/timeFDD"))->Fill(dgcand.timeFDDA(),dgcand.timeFDDC()) ;

    // check FIT information
//    auto bitMin = -1 + 16;
//    auto bitMax = 1 + 16;
//    for (auto bit = bitMin; bit <= bitMax; bit++) {
//      if (TESTBIT(dgcand.bbFT0Apf(), bit) ||
//          TESTBIT(dgcand.bbFT0Cpf(), bit) ||
//          TESTBIT(dgcand.bbFV0Apf(), bit) ||
//          TESTBIT(dgcand.bbFDDApf(), bit) ||
//          TESTBIT(dgcand.bbFDDCpf(), bit)) {
//        return;
//      }
//    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(5.,1.);

    //
    //here PID from TPC starts to be
    //
    //temporary control variables per event with combinatorics
    float tmpMomentum[4];//ok
    float tmpPt[4];//ok
    float tmpDedx[4];//ok
    float pi3invMass[4];//ok
    float pi3pt[4];//ok
    float pi3deltaPhi[4];
    float pi3assymav[4];
    float pi3vector[4];//ok
    float pi3scalar[4];//ok
    float pi3etasum[4];
    float deltaPhiTmp=0;
    float scalarPtsum=0;
    float nSigmaEl[4];
    float nSigmaPi[4];
    float dcaZ[4];
    float dcaXY[4];
    float chi2TPC[4];
    float chi2ITS[4];
    float nclTPCfind[4];

    //first loop to add all the tracks together
    TLorentzVector p1;
    p = TLorentzVector(0., 0., 0., 0.);
    for (auto trk : PVContributors) {
      p1.SetXYZM(trk.px(),trk.py(),trk.pz(),pion->Mass());
      p+=p1;
      scalarPtsum+=trk.pt();
    }
    float pttot  =  p.Pt();
    float mass4pi = p.Mag();

    TVector3 v1(0,0,0);
    TVector3 vtmp(0,0,0);
    float deltaphi=0;
    //remove combinatoric
    bool flagVcalPV[4]={false,false,false,false};

    //second loop to calculate 1 by 1 each combinatorial variable
    int counterTmp=0;
    for (auto trk : PVContributors) {
      v1.SetXYZ(trk.px(),trk.py(),trk.pz());
      for (auto trk1 : PVContributors) {
	if(trk.index() == trk1.index()) continue;
	vtmp.SetXYZ(trk1.px(),trk1.py(),trk1.pz());
	deltaphi = v1.Angle(vtmp);
	registry.get<TH1>(HIST("global/hDeltaAngleTrackPV"))->Fill(deltaphi);
	if(deltaphi < minAnglecut) {//default 0.05
	  flagVcalPV[counterTmp] = true;
	}
      }
      nSigmaEl[counterTmp] = trk.tpcNSigmaEl();
      nSigmaPi[counterTmp] = trk.tpcNSigmaPi();
      dcaZ[counterTmp] = trk.dcaZ();
      dcaXY[counterTmp] = trk.dcaXY();
      chi2TPC[counterTmp] = trk.tpcChi2NCl();
      chi2ITS[counterTmp] = trk.itsChi2NCl();
      nclTPCfind[counterTmp] = trk.tpcNClsFindable();

      p1.SetXYZM(trk.px(),trk.py(),trk.pz(),pion->Mass());
      tmpMomentum[counterTmp]=p1.P();
      tmpPt[counterTmp]=p1.Pt();
      tmpDedx[counterTmp]=trk.tpcSignal();

      deltaPhiTmp = CalculateDeltaPhi(p-p1,p1);
      pi3invMass[counterTmp]=(p-p1).Mag();
      pi3pt[counterTmp]=(p-p1).Pt();
      pi3deltaPhi[counterTmp]=deltaPhiTmp;
      pi3assymav[counterTmp]=(p1.Pt()-(scalarPtsum-p1.Pt())/3.)/(p1.Pt()+(scalarPtsum-p1.Pt())/3.);
      pi3vector[counterTmp]=(p+p1).Pt()/(p-p1).Pt();
      pi3scalar[counterTmp]=(p.Pt()-p1.Pt())/(p.Pt()+p1.Pt());
      pi3etasum[counterTmp]=(p-p1).Eta()+p1.Eta();

      counterTmp++;
    }

    //control histos, max 4 per event
    for (int i=0;i<4;i++){
      FillControlHistos<0>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], pi3etasum[i]);
      registry.get<TH2>(HIST("control/h0Cut3piMassVsPt"))->Fill(pi3invMass[i],pi3pt[i]);
      registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut0"))->Fill(tmpMomentum[i],tmpDedx[i]);
    }
    //control, 1 per event
    registry.get<TH1>(HIST("control/h0Cut4trkPtTot"))->Fill(pttot);
    registry.get<TH1>(HIST("control/h0Cut4piMass"))->Fill(mass4pi);
    registry.get<TH2>(HIST("control/h0Cut4trkMassVsPt"))->Fill(mass4pi,pttot);

    //remove combinatoric
    bool flagTotal[4]= {false,false,false,false};
    bool flagIM[4]   = {false,false,false,false};
    bool flagDP[4]   = {false,false,false,false};
    bool flagEl[4]   = {false,false,false,false};
    bool flagPi[4]   = {false,false,false,false};
    bool flagPt[4]   = {false,false,false,false};

    //bool flagVcalPV[4]={false,false,false,false};
    //float deltaphi=0;

    for (int i=0;i<4;i++){
      if(pi3invMass[i]<1.8)  { 
	flagIM[i]=true; 
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut2"))->Fill(tmpMomentum[i],tmpDedx[i]);
      }
      if(pi3deltaPhi[i]>1.6) { 
	flagDP[i]=true; 
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut3"))->Fill(tmpMomentum[i],tmpDedx[i]);
      }
      if(minNsigmaElcut < nSigmaEl[i] && nSigmaEl[i] < maxNsigmaElcut) { //default (-2,3)
	flagEl[i]=true; 
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut4"))->Fill(tmpMomentum[i],tmpDedx[i]);
      } 
      if(TMath::Abs(nSigmaPi[i]) > maxNsigmaPiVetocut) { //default is 4
	flagPi[i]=true;
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut5"))->Fill(tmpMomentum[i],tmpDedx[i]);
      } 
      if(tmpPt[i] > minPtEtrkcut) { //0.25
	flagPt[i]=true;
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut6"))->Fill(tmpMomentum[i],tmpDedx[i]);
      } 
      if(!flagVcalPV[i]) {
	registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut7"))->Fill(tmpMomentum[i],tmpDedx[i]);
      }
      flagTotal[i] = flagEl[i] && flagPi[i] && flagPt[i] && !flagVcalPV[i];
    }

    int counterM3pi = flagIM[0] + flagIM[1] + flagIM[2] + flagIM[3];
    int counterDphi = flagDP[0] + flagDP[1] + flagDP[2] + flagDP[3];
    int counterEl   = flagEl[0] + flagEl[1] + flagEl[2] + flagEl[3];
    int counterPi   = flagPi[0] + flagPi[1] + flagPi[2] + flagPi[3];
    int counterPt   = flagPt[0] + flagPt[1] + flagPt[2] + flagPt[3];
    int counterVcal = !flagVcalPV[0] + !flagVcalPV[1] + !flagVcalPV[2] + !flagVcalPV[3];
    int counterTotal = flagTotal[0] + flagTotal[1] + flagTotal[2] + flagTotal[3];

    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(5.+counterM3pi,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(10.+counterDphi,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(15.+counterEl,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(20.+counterPi,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(25.+counterPt,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(30.+counterVcal,1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(35.+counterTotal,1.);


    //draw control histograms
    if(counterEl > 0){//Nelectrons>0
      for (int i=0;i<4;i++){
	if(flagEl[i]) registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut10"))->Fill(tmpMomentum[i],tmpDedx[i]);
      }
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(6.,1.);
      if(flagEl[0] * flagPi[0] + flagEl[1] * flagPi[1] + flagEl[2] * flagPi[2] + flagEl[3] * flagPi[3] > 0){//pi veto
	for (int i=0;i<4;i++){
	  if(flagEl[i] && flagPi[i]) registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut11"))->Fill(tmpMomentum[i],tmpDedx[i]);
	}
	registry.get<TH1>(HIST("global/hEventEff"))->Fill(7.,1.);
	if(flagEl[0] * flagPi[0] * !flagVcalPV[0] + 
	   flagEl[1] * flagPi[1] * !flagVcalPV[1] + 
	   flagEl[2] * flagPi[2] * !flagVcalPV[2] + 
	   flagEl[3] * flagPi[3] * !flagVcalPV[3] > 0){//vcal veto
	  for (int i=0;i<4;i++){
	    if(flagEl[i] && flagPi[i] && !flagVcalPV[i]) registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut12"))->Fill(tmpMomentum[i],tmpDedx[i]);
	  }
	  registry.get<TH1>(HIST("global/hEventEff"))->Fill(8.,1.);
	  if(flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] + 
	     flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] + 
	     flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] + 
	     flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] > 0){//pT veto
	    for (int i=0;i<4;i++){
	      if(flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i]) registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut13"))->Fill(tmpMomentum[i],tmpDedx[i]);
	    }
	    registry.get<TH1>(HIST("global/hEventEff"))->Fill(9.,1.);
	    
	  } else {
	    if (verbose) {
	      LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT");
	    }
	  }
	} else {
	  if (verbose) {
	    LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal");
	  }
	}
      } else {
	if (verbose) {
	  LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by pi PID");
	}
      }
    } else {//no electron
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: no electron PID among 4 tracks");
      }
    }//end of Nelectrons check


    // check FIT information
    if(FITvetoFlag){
      auto bitMin = 16 - FITvetoWindow;//default is +- 1 bc (1 bit)
      auto bitMax = 16 + FITvetoWindow;
      for (auto bit = bitMin; bit <= bitMax; bit++) {
	if (TESTBIT(dgcand.bbFT0Apf(), bit) ||
	    TESTBIT(dgcand.bbFT0Cpf(), bit) ||
	    TESTBIT(dgcand.bbFV0Apf(), bit) ||
	    TESTBIT(dgcand.bbFDDApf(), bit) ||
	    TESTBIT(dgcand.bbFDDCpf(), bit)) {
	  return;
	}
      }
    }

    if(counterTotal==1){
      for (int i=0;i<4;i++){
        registry.get<TH1>(HIST("control/h1CutDcaZ"))->Fill(dcaZ[i]);
        registry.get<TH1>(HIST("control/h1CutDcaXY"))->Fill(dcaXY[i]);
        registry.get<TH1>(HIST("control/h1CutChi2TPC"))->Fill(chi2TPC[i]);
        registry.get<TH1>(HIST("control/h1CutChi2ITS"))->Fill(chi2ITS[i]);

	if(flagTotal[i]) {
	  FillControlHistos<1>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], pi3etasum[i]);
	  registry.get<TH2>(HIST("control/h1Cut3piMassVsPt"))->Fill(pi3invMass[i],pi3pt[i]);
	  registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut1"))->Fill(tmpMomentum[i],tmpDedx[i]);
	  registry.get<TH1>(HIST("global/hFinalPtSpectrumEl"))->Fill(tmpPt[i]);
	  registry.get<TH1>(HIST("control/h1CutTPCnclsFindable"))->Fill(nclTPCfind[i]);
	}
      }
      registry.get<TH1>(HIST("control/h1Cut4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/h1Cut4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/h1Cut4trkMassVsPt"))->Fill(mass4pi,pttot);

      registry.get<TH1>(HIST("global/hEventEff"))->Fill(10.,1.);
    } else {//more than 1 electron candidate
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: more than one electron candidate");
      }
    }//end of Nelectrons check

    //trk ITS : eta vs n its layers - done earlier



  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauTau13topo>(cfgc, TaskName{"TauTau13topo"})
  };
}
