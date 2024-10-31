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
#include "Framework/AnalysisDataModel.h"

// #include "TDatabasePDG.h" // not recommended in o2
#include "Framework/O2DatabasePDGPlugin.h"

#include "TLorentzVector.h"
// #include "Common/DataModel/EventSelection.h"
// #include "Common/DataModel/TrackSelectionTables.h"

#include "Common/DataModel/PIDResponse.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/Core/DGPIDSelector.h"
#include "PWGUD/Core/SGSelector.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct TauTau13topo {
  SGSelector sgSelector;
  // configurables
  Configurable<float> FV0_cut{"FV0", 1000., "FV0A threshold"};
  Configurable<float> FT0A_cut{"FT0A", 150., "FT0A threshold"};
  Configurable<float> FT0C_cut{"FT0C", 50., "FT0C threshold"};
  Configurable<float> ZDC_cut{"ZDC", 100., "ZDC threshold"};
  Configurable<float> gap_Side{"gap", 2, "gap selection"};
  //  ConfigurableAxis ptAxis{"pAxis", {100, 0., 5.}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis ptAxis{"ptAxis", {120, 0., 4.}, "#it{p} (GeV/#it{c})"};
  //  ConfigurableAxis etaAxis{"etaAxis", {100, -2., 2.}, "#eta"};
  ConfigurableAxis dedxAxis{"dedxAxis", {100, 20., 160.}, "dE/dx"};
  ConfigurableAxis minvAxis{"MinvAxis", {100, 0.4, 3.5}, "M_{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis phiAxis{"phiAxis", {120, 0., 3.2}, "#phi"};
  //  ConfigurableAxis vectorAxis{"vectorAxis", {100, 0., 2.}, "A_{V}"};
  //  ConfigurableAxis scalarAxis{"scalarAxis", {100, -1., 1.}, "A_{S}"};
  Configurable<bool> verbose{"Verbose", {}, "Additional print outs"};

  // cut selection configurables
  Configurable<float> zvertexcut{"Zvertexcut", 15., "Z vertex cut"};
  Configurable<float> trkEtacut{"TrkEtacut", 1.5, "max track eta cut"};
  Configurable<bool> sameSign{"sameSign", {}, "Switch: same(true) or opposite(false) sign"};
  Configurable<float> ptTotcut{"PtTotcut", 0.15, "min pt of all 4 tracks cut"};
  Configurable<float> minAnglecut{"minAnglecut", 0.05, "min angle between tracks cut"};
  Configurable<float> minNsigmaElcut{"minNsigmaElcut", -2., "min Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaElcut{"maxNsigmaElcut", 3., "max Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaPiVetocut{"maxNsigmaPiVetocut", 4., "max Nsigma for Pion veto cut"};
  Configurable<float> maxNsigmaPrVetocut{"maxNsigmaPrVetocut", 3., "max Nsigma for Proton veto cut"};
  Configurable<float> maxNsigmaKaVetocut{"maxNsigmaKaVetocut", 3., "max Nsigma for Kaon veto cut"};
  Configurable<float> minPtEtrkcut{"minPtEtrkcut", 0.25, "min Pt for El track cut"};
  Configurable<bool> FITvetoFlag{"FITvetoFlag", {}, "To apply FIT veto"};
  Configurable<int> FITvetoWindow{"FITvetoWindow", 1, "FIT veto window"};
  Configurable<bool> useFV0ForVeto{"useFV0ForVeto", 0, "use FV0 for veto"};
  Configurable<bool> useFDDAForVeto{"useFDDAForVeto", 0, "use FDDA for veto"};
  Configurable<bool> useFDDCForVeto{"useFDDCForVeto", 0, "use FDDC for veto"};
  Configurable<int> nTofTrkMinCut{"nTofTrkMinCut", 0, "min TOF tracks"};

  Configurable<bool> invMass3piSignalRegion{"invMass3piSignalRegion", 1, "1-use inv mass 3pi in signal region, 0-in background region"};
  Configurable<float> invMass3piMaxcut{"invMass3piMaxcut", 20., "Z invariant mass of 3 pi cut"};
  Configurable<float> deltaPhiMincut{"deltaPhiMincut", 0., "delta phi electron - 3 pi direction cut"};
  Configurable<int> nTPCcrossedRowsMinCut{"nTPCcrossedRowsMinCut", 50, "min N_crossed TPC rows for electron candidate"};
  Configurable<float> nSigma3piMaxCut{"nSigma3piMaxCut", 20., "n sigma 3 pi max cut"};

  // Configurable<bool> DGactive{"DGactive", false, "Switch on DGproducer"};
  // Configurable<bool> SGactive{"SGactive", true, "Switch on SGproducer"};

  // a pdg object
  // TDatabasePDG* pdg = nullptr; //not recommended
  Service<o2::framework::O2DatabasePDG> pdg;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    // pdg = TDatabasePDG::Instance();

    // dgcandidates histograms
    const AxisSpec axisp{100, 0., 5., "#it{p} (GeV/#it{c})"};
    const AxisSpec axispt{ptAxis, "p_{T} axis"};
    // const AxisSpec axiseta{etaAxis, "#eta - pseudo rapidity axis"};
    const AxisSpec axiseta{100, -2., 2., "#eta"};
    const AxisSpec axisdedx{dedxAxis, "dEdx axis"};
    const AxisSpec axisminv{minvAxis, "invariant mass axis"};
    const AxisSpec axisphi{phiAxis, "phi axis"};
    //    const AxisSpec axisav{vectorAxis, "AV axis"};
    //    const AxisSpec axisas{scalarAxis, "AS axis"};
    const AxisSpec vectorAxis{100, 0., 2., "A_{V}"};
    const AxisSpec scalarAxis{100, -1., 1., "A_{S}"};
    const AxisSpec axisZDC{50, -1., 14., "#it{E} (TeV)"};

    registry.add("global/GapSide", "Associated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});
    registry.add("global/GapSideTrue", "Recalculated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});

    registry.add("global/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry.add("global/hZNACtime", "ZNA vs ZNC time; #it{time}_{ZNA} (ns); #it{time}_{ZNC} (ns); Collisions", {HistType::kTH2F, {{100, -10., 10.}, {100, -10., 10.}}});
    // registry.add("global/hZNACenergyTest", "ZNA or ZNC energy; #it{E}_{ZNA,ZNC} (GeV); Collisions", {HistType::kTH1F, {{100,-1000,0}}});

    registry.add("global/hVertexXY", "Vertex position in x and y direction; #it{V}_{x} (cm); #it{V}_{y} (cm); Collisions", {HistType::kTH2F, {{50, -0.05, 0.05}, {50, -0.02, 0.02}}});
    registry.add("global/hVertexZ", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    // registry.add("global/hVertexZ15", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    // registry.add("global/hVertexZ10", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    registry.add("global/hNTracks", ";N_{tracks};events", {HistType::kTH1D, {{20, 0., 20.}}});
    // registry.add("global/hNTracksGlobal", ";N_{tracks,global};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global/hNTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("global/hTrackPtPV", ";p_T^{trk}; Entries", {HistType::kTH1F, {axispt}});
    registry.add("global/hTrackEtaPhiPV", ";Eta;Phi;", {HistType::kTH2D, {axiseta, {140, -3.5, 3.5}}});
    registry.add("global/hTrackEfficiencyPVGlobal", "0-All,1-ntpc,2-rat,3-chitpc,4chiits,5-dcaz,6-dcaxy,7pt,8eta;Track efficiency; Entries", {HistType::kTH1F, {{15, 0, 15}}});
    registry.add("global/hTrackEtaPhiPVGlobal", ";Eta;Phi;", {HistType::kTH2D, {axiseta, {140, -3.5, 3.5}}});

    registry.add("global/hSignalTPCvsPtPV", ";Pt;TPC Signal", {HistType::kTH2F, {axispt, {200, 0., 200}}});
    registry.add("global/hITSbitPVtrk", "ITS bit for PV tracks; Layer hit;Entries", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.add("global/hITSnbitsVsEtaPVtrk", "n ITS bits vs #eta for PV tracks; #eta;Layer hit;Entries", {HistType::kTH2F, {axiseta, {8, -1., 7.}}});
    registry.add("global/hITSbitVsEtaPVtrk", "ITS bit vs #eta for PV tracks; #eta;Layer hit;Entries", {HistType::kTH2F, {axiseta, {8, 0., 8.}}});
    registry.add("global/hEventEff", "Event cut efficiency: 0-All,1-PV=4,2-Qtot=0,3-El;Cut;entries", {HistType::kTH1F, {{25, 0., 25.}}});
    registry.add("global/hNCombAfterCut", "Combinations after cut: 0-All,5-M3pi,10-Dphi,15-N_{e},20-N_{v#pi},25-Pt,30-Vcal,35-N_{vp},40-N_{vK},45-Tot;N_{comb};entries", {HistType::kTH1F, {{60, 0., 60.}}});
    // registry.add("global/hInvMassElTrack", ";M_{inv}^{2};entries", {HistType::kTH1F, {{100, -0.01, 0.49}}});
    registry.add("global/hDeltaAngleTrackPV", ";#Delta#alpha;entries", {HistType::kTH1F, {{136, 0., 3.2}}}); // 0.49
    registry.add("global/hTrkCheck", ";track type;entries", {HistType::kTH1F, {{16, -1, 15}}});

    registry.add("global1/hVertexZ", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    registry.add("global1/hNTracks", ";N_{tracks};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global1/hNTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registry.add("global1/hTrackPtPV", ";p_T^{trk}; Entries", {HistType::kTH1F, {axispt}});
    registry.add("global1/hTrackEtaPhiPV", ";Eta;Phi;", {HistType::kTH2D, {axiseta, {140, -3.5, 3.5}}});
    registry.add("global1/hTrackPVTotCharge", "Q_{Tot};Q_{Tot}; Entries", {HistType::kTH1F, {{11, -5, 6}}});

    // cut0
    registry.add("control/cut0/h3piMassComb", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut0/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut0/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut0/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut0/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut0/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut0/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut0/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut0/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut0/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut0/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut0/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut0/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut0/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut0/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut0/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    registry.add("control/cut0/hInvMass2ElAll", "Inv Mass of 2 Electrons from coherent peak;M_{inv}^{2e};entries", {HistType::kTH1F, {{150, -0.1, 9.}}});
    registry.add("control/cut0/hInvMass2ElCoh", "Inv Mass of 2 Electrons from coherent peak;M_{inv}^{2e};entries", {HistType::kTH1F, {{150, -0.1, 4.}}});
    registry.add("control/cut0/hGamPtCoh", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut0/hGamPtCohIM0", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut0/hN2gamma", "Number of gamma pairs among 3 comb;N_{#gamma#gamma};entries", {HistType::kTH1F, {{20, 0., 20.}}});
    registry.add("control/cut0/hGamAS", ";A_{S};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut0/hGamAV", ";A_{V};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut0/hInvMass2GamCoh", "Inv Mass of 2 Gamma from coherent peak;M_{inv}^{2#gamma};entries", {HistType::kTH1F, {{160, 0.5, 4.5}}});
    registry.add("control/cut0/hDeltaPhi2GamCoh", "Delta Phi of 2 Gamma from coherent peak;#Delta#phi^{2#gamma};entries", {HistType::kTH1F, {phiAxis}});

    // cut1
    registry.add("control/cut1/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut1/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut1/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut1/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut1/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut1/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut1/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut1/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut1/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut1/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut1/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut1/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut1/hDcaZ", "All 4 tracks dca ;dca_{Z};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    registry.add("control/cut1/hDcaXY", "All 4 tracks dca ;dca_{XY};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    registry.add("control/cut1/hChi2TPC", "All 4 tracks Chi2 ;Chi2_{TPC};entries", {HistType::kTH1F, {{48, -2, 10.}}});
    registry.add("control/cut1/hChi2ITS", "All 4 tracks Chi2 ;Chi2_{ITS};entries", {HistType::kTH1F, {{44, -2, 20.}}});
    registry.add("control/cut1/hTPCnclsFindable", "All 4 tracks NclFind ;N_{TPC,cl,findable};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut1/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut1/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut1/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut1/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    //   // cut1a for 2<m(4pi)<2.3 GeV
    //    registry.add("control/cut1/cut1a/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registry.add("control/cut1/cut1a/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registry.add("control/cut1/cut1a/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registry.add("control/cut1/cut1a/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registry.add("control/cut1/cut1a/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registry.add("control/cut1/cut1a/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut1/cut1a/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    //    registry.add("control/cut1/cut1a/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //    registry.add("control/cut1/cut1a/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut1/cut1a/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //    registry.add("control/cut1/cut1a/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //    registry.add("control/cut1/cut1a/hDcaZ", "All 4 tracks dca ;dca_{Z};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    //    registry.add("control/cut1/cut1a/hDcaXY", "All 4 tracks dca ;dca_{XY};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    //    registry.add("control/cut1/cut1a/hChi2TPC", "All 4 tracks Chi2 ;Chi2_{TPC};entries", {HistType::kTH1F, {{48, -2, 10.}}});
    //    registry.add("control/cut1/cut1a/hChi2ITS", "All 4 tracks Chi2 ;Chi2_{ITS};entries", {HistType::kTH1F, {{44, -2, 20.}}});
    //    registry.add("control/cut1/cut1a/hTPCnclsFindable", "All 4 tracks NclFind ;N_{TPC,cl,findable};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registry.add("control/cut1/cut1a/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut1/cut1a/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //    registry.add("control/cut1/cut1a/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});

    // cut20
    registry.add("control/cut20/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut20/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut20/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut20/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut20/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut20/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut20/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut20/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut20/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut20/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut20/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut20/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut20/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut20/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut20/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut20/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    registry.add("control/cut20/hInvMass2ElCoh", "Inv Mass of 2 Electrons from coherent peak;M_{inv}^{2e};entries", {HistType::kTH1F, {{150, -0.1, 4.}}});
    registry.add("control/cut20/hGamPtCohIM0", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut20/hGamAS", ";A_{S};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut20/hGamAV", ";A_{V};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut20/hInvMass2GamCoh", "Inv Mass of 2 Gamma from coherent peak;M_{inv}^{2#gamma};entries", {HistType::kTH1F, {{160, 0., 4.}}});
    registry.add("control/cut20/hDeltaPhi2GamCoh", "Delta Phi of 2 Gamma from coherent peak;#Delta#phi^{2#gamma};entries", {HistType::kTH1F, {phiAxis}});

    // cut21
    registry.add("control/cut21/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut21/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut21/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut21/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut21/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut21/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut21/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut21/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut21/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut21/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut21/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut21/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut21/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut21/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut21/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut21/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut22
    registry.add("control/cut22/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut22/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut22/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut22/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut22/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut22/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut22/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut22/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut22/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut22/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut22/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut22/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut22/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut22/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut22/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut22/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut23
    registry.add("control/cut23/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut23/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut23/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut23/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut23/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut23/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut23/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut23/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut23/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut23/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut23/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut23/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut23/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut23/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut23/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut23/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut24
    registry.add("control/cut24/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut24/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut24/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut24/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut24/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut24/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut24/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut24/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut24/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut24/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut24/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut24/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut24/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut24/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut24/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut24/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut25
    registry.add("control/cut25/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut25/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut25/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut25/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut25/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut25/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut25/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut25/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut25/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut25/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut25/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut25/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut25/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut25/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut25/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut25/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut26
    registry.add("control/cut26/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut26/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut26/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut26/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut26/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut26/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut26/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut26/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut26/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut26/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut26/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut26/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut26/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut26/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut26/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut26/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut27
    registry.add("control/cut27/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut27/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut27/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut27/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut27/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut27/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut27/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut27/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut27/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut27/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut27/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut27/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut27/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut27/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut27/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut27/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut28
    registry.add("control/cut28/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut28/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut28/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut28/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut28/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut28/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut28/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut28/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut28/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut28/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut28/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut28/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut28/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut28/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut28/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut28/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut29
    registry.add("control/cut29/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut29/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut29/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut29/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut29/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut29/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut29/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut29/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut29/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut29/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut29/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut29/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut29/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut29/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut29/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut29/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut30
    registry.add("control/cut30/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut30/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut30/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut30/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut30/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut30/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut30/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut30/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut30/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut30/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut30/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut30/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut30/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut30/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut30/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut30/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut31
    registry.add("control/cut31/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut31/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut31/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut31/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut31/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut31/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut31/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut31/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut31/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut31/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut31/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut31/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut31/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut31/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut31/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut31/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});

    // cut32
    registry.add("control/cut32/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut32/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut32/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut32/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut32/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut32/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut32/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut32/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut32/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut32/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut32/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut32/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut32/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut32/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut32/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut32/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry.add("control/cut32/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});

    // cut33
    registry.add("control/cut33/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut33/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut33/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut33/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut33/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut33/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut33/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut33/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut33/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut33/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut33/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut33/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut33/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut33/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut33/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut33/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry.add("control/cut33/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});

    // cut34
    registry.add("control/cut34/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registry.add("control/cut34/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut34/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registry.add("control/cut34/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registry.add("control/cut34/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registry.add("control/cut34/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registry.add("control/cut34/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    registry.add("control/cut34/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registry.add("control/cut34/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut34/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut34/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registry.add("control/cut34/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registry.add("control/cut34/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("control/cut34/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registry.add("control/cut34/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry.add("control/cut34/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry.add("control/cut34/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});

    // pid El
    registry.add("pidTPC/hpvsdedxElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut0CohPsi2s", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    registry.add("pidTPC/hpvsdedxElHipCut1", "All hip;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // pid separately for each cut (what we reject)
    registry.add("pidTPC/hpvsdedxElHipCut2", "rejected, IM hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut3", "rejected, DP hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut4", "rejected, El hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut5", "rejected, vPi hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut6", "rejected, Pt hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut7", "rejected, vVc hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut8", "rejected, pTot hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut9", "rejected, vPr hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut10", "rejected, vKa hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut11", "rejected, nCR hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut12", "rejected, s3pi hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut13", "rejected, s3pi hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // pid sequentialy
    registry.add("pidTPC/hpvsdedxElHipCut20", "El hip;    #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut21", "vPi+20 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut22", "vVc+21 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut23", "Pt+22 hip; #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut24", "vPr+23 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut26", "IM+25 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut27", "DP+26 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut28", "CR+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut29", "s3pi+28 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut30", "ptTot+29 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut31", "FIT+30 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut32", "TOF+31 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut33", "eTOF+1 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut34", "piTOF+33 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    registry.add("pidTPC/hpvsdedxElHipCut40", "All from gamma hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // pid 3pi
    // registry.add("pidTPC3pi/hpvsdedxElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC3pi/hpvsdedxElHipCut1", "All hip;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // pid sequentialy 3pi
    // registry.add("pidTPC3pi/hpvsdedxElHipCut20", "El hip;    #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut21", "vPi+20 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut22", "vVc+21 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut23", "Pt+22 hip; #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut24", "vPr+23 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut26", "ptTot+25 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut27", "FIT+26 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registry.add("pidTPC3pi/hpvsdedxElHipCut28", "TOF+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // final electron spectrum
    registry.add("global/hFinalPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});

    // fit histos
    registry.add("fit/bbFT0Abit", "FT0A bits;bit;entries", {HistType::kTH1F, {{32, 0., 32.}}});
    registry.add("fit/bbFT0Cbit", "FT0C bits;bit;entries", {HistType::kTH1F, {{32, 0., 32.}}});
    registry.add("fit/bbFV0Abit", "FV0A bits;bit;entries", {HistType::kTH1F, {{32, 0., 32.}}});
    registry.add("fit/bbFDDAbit", "FDDA bits;bit;entries", {HistType::kTH1F, {{32, 0., 32.}}});
    registry.add("fit/bbFDDCbit", "FDDC bits;bit;entries", {HistType::kTH1F, {{32, 0., 32.}}});
    registry.add("fit/bbFT0Aamplitude", "FT0A amplitude;Amplitude;entries", {HistType::kTH1F, {{100, -5., 95.}}});
    registry.add("fit/bbFT0Camplitude", "FT0C amplitude;Amplitude;entries", {HistType::kTH1F, {{100, -5., 95.}}});
    registry.add("fit/bbFT0ACamplitude", "FT0A vs FT0C amplitude;Amplitude FT0A;Amplitude FT0C;entries", {HistType::kTH2F, {{100, -5., 95.}, {100, -5., 95.}}});
    registry.add("fit/bbFV0Aamplitude", "FV0A amplitude;Amplitude;entries", {HistType::kTH1F, {{100, -5., 95.}}});
    registry.add("fit/bbFDDAamplitude", "FDDA amplitude;Amplitude;entries", {HistType::kTH1F, {{100, -5., 95.}}});
    registry.add("fit/bbFDDCamplitude", "FDDC amplitude;Amplitude;entries", {HistType::kTH1F, {{100, -5., 95.}}});
    registry.add("fit/bbFDDACamplitude", "FDDA vs FDDC amplitude;Amplitude FDDA;Amplitude FDDC;entries", {HistType::kTH2F, {{100, -5., 95.}, {100, -5., 95.}}});
    registry.add("fit/timeFT0", "FT0 time;time FT0A; time FT0C;entries", {HistType::kTH2F, {{100, -5., 35.}, {100, -5., 35.}}});
    registry.add("fit/timeFDD", "FDD time;time FDDA; time FDDC;entries", {HistType::kTH2F, {{100, -5., 35.}, {100, -5., 35.}}});
  }

  float CalculateDeltaPhi(TLorentzVector p, TLorentzVector p1)
  {
    float delta = p.Phi();
    if (delta < 0)
      delta += o2::constants::math::TwoPI;
    if (p1.Phi() < 0)
      delta -= (p1.Phi() + o2::constants::math::TwoPI);
    else
      delta -= p1.Phi();
    if (delta < 0)
      delta += o2::constants::math::TwoPI;
    if (delta > o2::constants::math::PI)
      delta = o2::constants::math::TwoPI - delta;
    return delta;
  }

  // fill control histograms per track
  template <int mode, typename T>
  //   void FillControlHistos(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float pi3etasum, float nCRtpc)
  void FillControlHistos(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float nCRtpc)
  {
    static constexpr std::string_view histoname[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                     "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                     "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                                                     "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"};
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h3piMassComb"))->Fill(pi3invMass);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h3trkPtTot"))->Fill(pi3pt);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/hDeltaPhi13topo"))->Fill(pi3deltaPhi);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h13AssymPt1ProngAver"))->Fill(pi3assymav);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h13Vector"))->Fill(pi3vector);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h13Scalar"))->Fill(pi3scalar);
    // registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/h13EtaSum"))->Fill(pi3etasum);
    registry.get<TH1>(HIST("control/cut") + HIST(histoname[mode]) + HIST("/hTPCnCrossedRows"))->Fill(nCRtpc);
  }

  template <typename T>
  int TrackCheck(T track)
  {
    // 1
    if (track.hasITS() && !track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 0;
    else if (!track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 1;
    else if (!track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 2;
    else if (!track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 3;
    // 2
    else if (track.hasITS() && track.hasTPC() && !track.hasTRD() && !track.hasTOF())
      return 4;
    else if (track.hasITS() && !track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 5;
    else if (track.hasITS() && !track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 6;
    else if (!track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 7;
    else if (!track.hasITS() && track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 8;
    else if (!track.hasITS() && !track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 9;
    // 3
    else if (track.hasITS() && track.hasTPC() && track.hasTRD() && !track.hasTOF())
      return 10;
    else if (track.hasITS() && track.hasTPC() && !track.hasTRD() && track.hasTOF())
      return 11;
    else if (track.hasITS() && !track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 12;
    else if (!track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 13;
    // 4
    else if (track.hasITS() && track.hasTPC() && track.hasTRD() && track.hasTOF())
      return 14;
    return -1;
  }

  //  using UDCollisionsFull = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>;
  //  using UDCollisionFull = UDCollisionsFull::iterator;
  using UDTracksFull = soa::Join<aod::UDTracks, aod::UDTracksPID, aod::UDTracksExtra, aod::UDTracksFlags, aod::UDTracksDCA>;
  using UDCollisionsFull2 = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>;
  using UDCollisionFull2 = UDCollisionsFull2::iterator;

  // PVContributors
  Filter PVContributorFilter = aod::udtrack::isPVContributor == true;
  using PVTracks = soa::Filtered<UDTracksFull>;

  //  void processDG(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  //  {
  //    int gapSide = 2;
  //    mainTask(gapSide,dgtracks);
  //  }
  //  PROCESS_SWITCH(TauTau13topo, processDG, "Process DG data", DGactive);

  // void processSG(UDCollisionFull2 const& dgcand, UDTracksFull const& dgtracks)
  void process(UDCollisionFull2 const& dgcand, UDTracksFull const& dgtracks, PVTracks const& PVContributors)
  {
    int gapSide = dgcand.gapSide();
    int truegapSide = sgSelector.trueGap(dgcand, FV0_cut, FT0A_cut, FT0C_cut, ZDC_cut);
    registry.fill(HIST("global/GapSide"), gapSide);
    registry.fill(HIST("global/GapSideTrue"), truegapSide);
    if (gapSide < 0 || gapSide > 2)
      return;
    gapSide = truegapSide;
    //    mainTask(gapSide,dgtracks);
    //  }
    //  PROCESS_SWITCH(TauTau13topo, processSG, "Process SG data", SGactive);
    //
    //  void mainTask(int gapSide, UDTracksFull const& dgtracks)
    //  {
    if (gapSide != gap_Side)
      return;
    // global checks
    registry.get<TH2>(HIST("global/hVertexXY"))->Fill(dgcand.posX(), dgcand.posY());
    registry.get<TH1>(HIST("global/hVertexZ"))->Fill(dgcand.posZ());
    //    if (TMath::Abs(dgcand.posZ()) < 10)
    //      registry.get<TH1>(HIST("global/hVertexZ10"))->Fill(dgcand.posZ());

    registry.get<TH1>(HIST("global/hNTracks"))->Fill(dgtracks.size());

    // setup PV tracks partition
    // Partition<UDTracksFull> PVContributors = aod::udtrack::isPVContributor == true;
    // PVContributors.bindTable(dgtracks);

    registry.get<TH1>(HIST("global/hNTracksPV"))->Fill(PVContributors.size());

    // zdc information
    float ZNAenergy = dgcand.energyCommonZNA();
    float ZNCenergy = dgcand.energyCommonZNC();
    // if (ZNAenergy < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(ZNAenergy);
    // if (ZNCenergy < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(ZNCenergy);
    if (ZNAenergy < 0)
      ZNAenergy = -1.;
    if (ZNCenergy < 0)
      ZNCenergy = -1.;
    registry.get<TH2>(HIST("global/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
    registry.get<TH2>(HIST("global/hZNACtime"))->Fill(dgcand.timeZNA(), dgcand.timeZNC());

    uint32_t clusterSizes;
    //    UChar_t clustermap1;
    //    bool isInnerITS = false;
    int nTofTrk = 0;
    int nEtaIn15 = 0;
    int nITSbits = -1;
    int npT100 = 0;
    TLorentzVector p;
    auto pionMass = pdg->Mass(211);
    auto electronMass = pdg->Mass(11);
    // TParticlePDG* pion = pdg->GetParticle(211);
    // TParticlePDG* electron = pdg->GetParticle(11);
    bool flagGlobalCheck = true;
    bool isGlobalTrack = true;
    int qtot = 0;
    // loop over PV contributors
    for (const auto& trk : PVContributors) {
      qtot += trk.sign();
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), pionMass);
      registry.get<TH1>(HIST("global/hTrackPtPV"))->Fill(p.Pt());
      if (std::abs(p.Eta()) < trkEtacut)
        nEtaIn15++; // 1.5 is a default
      registry.get<TH2>(HIST("global/hTrackEtaPhiPV"))->Fill(p.Eta(), p.Phi());

      if (trk.pt() > 0.1)
        npT100++;

      if (flagGlobalCheck) {
        registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(0., 1.);
        if (trk.tpcNClsCrossedRows() > 70) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(1., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.tpcNClsFindable() == 0) {
          isGlobalTrack = false;
        } else {
          if (trk.tpcNClsCrossedRows() / trk.tpcNClsFindable() > 0.8) {
            registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(2., 1.);
          } else {
            isGlobalTrack = false;
          }
        }

        if (trk.tpcChi2NCl() < 4.) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(3., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.itsChi2NCl() < 36.) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(4., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.dcaZ() < 2.) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(5., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.dcaXY() < 0.0105 * 0.035 / std::pow(trk.pt(), 1.1)) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(6., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.pt() > 0.1) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(7., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (std::abs(p.Eta()) < 0.8) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(8., 1.);
        } else {
          isGlobalTrack = false;
        }

        if (trk.hasITS()) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(9., 1.);

          // old version
          // clustermap1 = trk.itsClusterMap();
          // for (int bitNo = 0; bitNo < 7; bitNo++) {
          //   if (TESTBIT(clustermap1, bitNo)) { // check ITS bits/layers for each PV track
          //     registry.get<TH1>(HIST("global/hITSbitPVtrk"))->Fill(bitNo, 1.);
          //     registry.get<TH2>(HIST("global/hITSbitVsEtaPVtrk"))->Fill(p.Eta(), bitNo, 1.);
          //     nITSbits++;
          //   }
          // } // end of loop over ITS bits
          //
          // isInnerITS = TESTBIT(clustermap1, 0) || TESTBIT(clustermap1, 1) || TESTBIT(clustermap1, 2);
          // if (isInnerITS) {
          //   registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(10., 1.);
          // } else {
          //   isGlobalTrack = false;
          // }
          //
        } else {
          isGlobalTrack = false;
        }

        if (trk.hasTPC()) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(11., 1.);
        } else {
          isGlobalTrack = false;
        }

        // final global track
        if (isGlobalTrack) {
          registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(13., 1.);
          registry.get<TH2>(HIST("global/hTrackEtaPhiPVGlobal"))->Fill(p.Eta(), p.Phi());
        }
      } // end of flag check global

      // new version
      if (trk.hasITS()) { // ITS track
        nITSbits = -1;
        clusterSizes = trk.itsClusterSizes();
        // LOGF(info, "<tautau13topo> clistersize: %d", clusterSizes);
        for (int layer = 0; layer < 7; layer++) {
          if ((clusterSizes >> (layer * 4)) & 0xf) {
            registry.get<TH1>(HIST("global/hITSbitPVtrk"))->Fill(layer, 1.);
            registry.get<TH2>(HIST("global/hITSbitVsEtaPVtrk"))->Fill(p.Eta(), layer, 1.);
            nITSbits++;
          }
        } // end of loop over ITS bits
        // isInnerITS = TESTBIT(clustermap1, 0) || TESTBIT(clustermap1, 1) || TESTBIT(clustermap1, 2);
        // if (isInnerITS) {
        //   registry.get<TH1>(HIST("global/hTrackEfficiencyPVGlobal"))->Fill(10., 1.);
        // } else {
        //   isGlobalTrack = false;
        // }
      } // has ITS
      registry.get<TH2>(HIST("global/hITSnbitsVsEtaPVtrk"))->Fill(p.Eta(), nITSbits);
      if (trk.hasTPC())
        registry.get<TH2>(HIST("global/hSignalTPCvsPtPV"))->Fill(p.Pt(), trk.tpcSignal());
      if (trk.hasTOF())
        nTofTrk++;
    } // end of loop over PV tracks
    registry.get<TH1>(HIST("global/hNtofTrk"))->Fill(nTofTrk);

    //
    // selection
    //
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(0., 1.);

    // skip events with too few/many tracks
    // if (PVContributors.size() != 4 || dgcand.numContrib() != 4) {
    if (PVContributors.size() != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Number of PV contributors is %d, candContriv %d", PVContributors.size(), dgcand.numContrib());
      }
      return;
    }
    // old version in DG producer
    // if (dgcand.numContrib() != 4) {
    //   if (verbose) {
    //     LOGF(info, "<tautau13topo> Candidate rejected: Number of PV contributors is %d", dgcand.numContrib());
    //   }
    //   return;
    // }

    registry.get<TH1>(HIST("global/hEventEff"))->Fill(1., 1.);
    // registry.get<TH1>(HIST("global1/hTrackPVTotCharge"))->Fill(dgcand.netCharge());
    registry.get<TH1>(HIST("global1/hTrackPVTotCharge"))->Fill(qtot);
    registry.get<TH1>(HIST("global1/hVertexZ"))->Fill(dgcand.posZ());
    registry.get<TH1>(HIST("global1/hNTracks"))->Fill(dgtracks.size());
    registry.get<TH1>(HIST("global1/hNTracksPV"))->Fill(PVContributors.size());
    for (const auto& trk : PVContributors) {
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), pionMass);
      registry.get<TH1>(HIST("global1/hTrackPtPV"))->Fill(p.Pt());
      registry.get<TH2>(HIST("global1/hTrackEtaPhiPV"))->Fill(p.Eta(), p.Phi());
    }

    // if vz < 15
    if (std::abs(dgcand.posZ()) >= zvertexcut) { // default = 15
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: VertexZ is %f", dgcand.posZ());
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(2., 1.);

    // if eta tracks <1.5
    if (nEtaIn15 != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Ntrk inside |eta|<1.5 is %d", nEtaIn15);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(3., 1.);

    // skip events with net charge != 0
    if (!sameSign) { // opposite sign is signal
      // if (dgcand.netCharge() != 0) {
      if (qtot != 0) {
        if (verbose) {
          LOGF(info, "<tautau13topo> Candidate rejected: Net charge is %d (dgcand %d), while should be 0", qtot, dgcand.netCharge());
        }
        return;
      }
    } else { // same sign is background
      // if (dgcand.netCharge() == 0) {
      if (qtot == 0) {
        if (verbose) {
          LOGF(info, "<tautau13topo> Candidate rejected: Net charge is %d (dgcand %d), while should be not 0", qtot, dgcand.netCharge());
        }
        return;
      }
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(4., 1.);

    //    // skip events with out-of-range rgtrwTOF (fraction-of-good-tracks-with-TOF-hit)
    //    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    //    if (rtrwTOF < 0.25) {
    //      if (verbose) {
    //        LOGF(debug, "<tautau13topo> Candidate rejected: rtrwTOF is %f", rtrwTOF);
    //      }
    //      return;
    //    }
    //
    // FIT informaton
    for (auto bit = 0; bit <= 32; bit++) {
      registry.get<TH1>(HIST("fit/bbFT0Abit"))->Fill(bit, TESTBIT(dgcand.bbFT0Apf(), bit));
      registry.get<TH1>(HIST("fit/bbFT0Cbit"))->Fill(bit, TESTBIT(dgcand.bbFT0Cpf(), bit));
      registry.get<TH1>(HIST("fit/bbFV0Abit"))->Fill(bit, TESTBIT(dgcand.bbFV0Apf(), bit));
      registry.get<TH1>(HIST("fit/bbFDDAbit"))->Fill(bit, TESTBIT(dgcand.bbFDDApf(), bit));
      registry.get<TH1>(HIST("fit/bbFDDCbit"))->Fill(bit, TESTBIT(dgcand.bbFDDCpf(), bit));
    }
    registry.get<TH1>(HIST("fit/bbFT0Aamplitude"))->Fill(dgcand.totalFT0AmplitudeA());
    registry.get<TH1>(HIST("fit/bbFT0Camplitude"))->Fill(dgcand.totalFT0AmplitudeC());
    registry.get<TH2>(HIST("fit/bbFT0ACamplitude"))->Fill(dgcand.totalFT0AmplitudeA(), dgcand.totalFT0AmplitudeC());
    registry.get<TH1>(HIST("fit/bbFV0Aamplitude"))->Fill(dgcand.totalFV0AmplitudeA());
    registry.get<TH1>(HIST("fit/bbFDDAamplitude"))->Fill(dgcand.totalFDDAmplitudeA());
    registry.get<TH1>(HIST("fit/bbFDDCamplitude"))->Fill(dgcand.totalFDDAmplitudeC());
    registry.get<TH2>(HIST("fit/bbFDDACamplitude"))->Fill(dgcand.totalFDDAmplitudeA(), dgcand.totalFDDAmplitudeC());

    registry.get<TH2>(HIST("fit/timeFT0"))->Fill(dgcand.timeFT0A(), dgcand.timeFT0C());
    registry.get<TH2>(HIST("fit/timeFDD"))->Fill(dgcand.timeFDDA(), dgcand.timeFDDC());

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
    // registry.get<TH1>(HIST("global/hEventEff"))->Fill(5., 1.);

    // if pt tracks >0.100 GeV/c
    if (npT100 != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: number of tracks with pT>0.1GeV/c is %d", npT100);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(5., 1.);

    //
    // here PID from TPC starts to be
    //
    // temporary control variables per event with combinatorics
    float tmpMomentum[4];
    float tmpPt[4];
    float tmpDedx[4];
    float pi3invMass[4];
    float pi3pt[4];
    float pi3deltaPhi[4];
    float pi3assymav[4];
    float pi3vector[4];
    float pi3scalar[4];
    //    float pi3etasum[4];
    float deltaPhiTmp = 0;
    float scalarPtsum = 0;
    float nSigmaEl[4];
    float nSigmaPi[4];
    float nSigma3Pi[4] = {0., 0., 0., 0.};
    float nSigmaPr[4];
    float nSigmaKa[4];
    float dcaZ[4];
    float dcaXY[4];
    float chi2TPC[4];
    float chi2ITS[4];
    float nclTPCfind[4];
    float nclTPCcrossedRows[4];
    bool tmpHasTOF[4];
    float mass3pi1e[4];
    double trkTime[4];
    float trkTimeRes[4];
    double trkTimeTot = 0.;
    float trkTimeResTot = 10000.;

    // 2 gamma from 4 electrons
    // 12 34 | 01 23 |//1 //6 | 0 5 |counter<3?counter:5-counter counter<3?0:1
    // 13 24 | 02 13 |//2 //5 | 1 4 |
    // 14 23 | 03 12 |//3 //4 | 2 3 |
    TLorentzVector p1, p2;
    TLorentzVector gammaPair[3][2];
    float invMass2El[3][2];
    int counterTmp = 0;
    bool flagIMGam2ePV[4] = {true, true, true, true};

    for (const auto& trk : PVContributors) {
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), electronMass);
      for (const auto& trk1 : PVContributors) {
        if (trk.index() >= trk1.index())
          continue;
        p1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), electronMass);
        invMass2El[(counterTmp < 3 ? counterTmp : 5 - counterTmp)][(counterTmp < 3 ? 0 : 1)] = (p + p1).Mag2();
        gammaPair[(counterTmp < 3 ? counterTmp : 5 - counterTmp)][(counterTmp < 3 ? 0 : 1)] = (p + p1);
        registry.get<TH1>(HIST("control/cut0/hInvMass2ElAll"))->Fill((p + p1).Mag2());
        counterTmp++;
        if ((p + p1).M() < 0.015) {
          flagIMGam2ePV[trk.index()] = false;
          flagIMGam2ePV[trk1.index()] = false;
        }
      } // end of loop over PVContributors
    } // end of loop over PVContributors

    // first loop to add all the tracks together
    p = TLorentzVector(0., 0., 0., 0.);
    for (const auto& trk : PVContributors) {
      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), pionMass);
      p += p1;
      scalarPtsum += trk.pt();
    } // end of loop over PVContributors

    float pttot = p.Pt();
    float mass4pi = p.Mag();

    TVector3 v1(0, 0, 0);
    TVector3 vtmp(0, 0, 0);
    float deltaphi = 0;
    // remove combinatoric
    bool flagVcalPV[4] = {false, false, false, false};

    // second loop to calculate 1 by 1 each combinatorial variable
    counterTmp = 0;
    int tmpTrkCheck = -1;
    for (const auto& trk : PVContributors) {
      tmpTrkCheck = TrackCheck(trk); // check detectors associated to track
      registry.get<TH1>(HIST("global/hTrkCheck"))->Fill(tmpTrkCheck);

      // inv mass of 3pi + 1e
      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), pionMass);
      p2.SetXYZM(trk.px(), trk.py(), trk.pz(), electronMass);
      mass3pi1e[counterTmp] = (p - p1 + p2).Mag();

      v1.SetXYZ(trk.px(), trk.py(), trk.pz());
      for (const auto& trk1 : PVContributors) {
        if (trk.index() == trk1.index())
          continue;
        vtmp.SetXYZ(trk1.px(), trk1.py(), trk1.pz());
        deltaphi = v1.Angle(vtmp);
        registry.get<TH1>(HIST("global/hDeltaAngleTrackPV"))->Fill(deltaphi);
        if (deltaphi < minAnglecut) { // default 0.05
          flagVcalPV[counterTmp] = true;
        }
      } // end of loop over PVContributors
      nSigmaEl[counterTmp] = trk.tpcNSigmaEl();
      nSigmaPi[counterTmp] = trk.tpcNSigmaPi();
      nSigma3Pi[3] += (nSigmaPi[counterTmp] * nSigmaPi[counterTmp]);
      nSigmaPr[counterTmp] = trk.tpcNSigmaPr();
      nSigmaKa[counterTmp] = trk.tpcNSigmaKa();
      dcaZ[counterTmp] = trk.dcaZ();
      dcaXY[counterTmp] = trk.dcaXY();
      chi2TPC[counterTmp] = trk.tpcChi2NCl();
      chi2ITS[counterTmp] = trk.itsChi2NCl();
      nclTPCfind[counterTmp] = trk.tpcNClsFindable();
      nclTPCcrossedRows[counterTmp] = trk.tpcNClsCrossedRows();
      tmpHasTOF[counterTmp] = trk.hasTOF();
      trkTime[counterTmp] = trk.trackTime();
      trkTimeRes[counterTmp] = trk.trackTimeRes();

      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), pionMass);
      tmpMomentum[counterTmp] = p1.P();
      tmpPt[counterTmp] = p1.Pt();
      tmpDedx[counterTmp] = trk.tpcSignal();

      deltaPhiTmp = CalculateDeltaPhi(p - p1, p1);
      pi3invMass[counterTmp] = (p - p1).Mag();
      pi3pt[counterTmp] = (p - p1).Pt();
      pi3deltaPhi[counterTmp] = deltaPhiTmp;
      pi3assymav[counterTmp] = (p1.Pt() - (scalarPtsum - p1.Pt()) / 3.) / (p1.Pt() + (scalarPtsum - p1.Pt()) / 3.);
      pi3vector[counterTmp] = (p + p1).Pt() / (p - p1).Pt();
      pi3scalar[counterTmp] = (p.Pt() - p1.Pt()) / (p.Pt() + p1.Pt());
      //      pi3etasum[counterTmp] = (p - p1).Eta() + p1.Eta();

      counterTmp++;
    } // end of loop over PVContributors

    // calculate trk time and resolution total
    // 1. find best resolution
    int iTmpBest = -1;
    for (int i = 0; i < 4; i++) {
      if (trkTimeRes[i] < trkTimeResTot) {
        trkTimeResTot = trkTimeRes[i];
        iTmpBest = i;
      }
    }
    // 2. use best resol to calculate total time
    for (int i = 0; i < 4; i++) {
      if (i == iTmpBest)
        continue;
      trkTimeTot += fabs(trkTime[iTmpBest] - trkTime[i]);
    }
    trkTimeResTot = std::sqrt(trkTimeRes[0] * trkTimeRes[0] +
                              trkTimeRes[1] * trkTimeRes[1] +
                              trkTimeRes[2] * trkTimeRes[2] +
                              trkTimeRes[3] * trkTimeRes[3]);

    // control histos, max 4 per event, cut0
    for (int i = 0; i < 4; i++) {
      FillControlHistos<0>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
      registry.get<TH2>(HIST("control/cut0/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
      registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut0"))->Fill(tmpMomentum[i], tmpDedx[i]);
      // nsigma3Pi calculation
      nSigma3Pi[i] = nSigma3Pi[3] - (nSigmaPi[i] * nSigmaPi[i]);
      nSigma3Pi[i] = std::sqrt(nSigma3Pi[i]);
      registry.get<TH1>(HIST("control/cut0/hsigma3Pi"))->Fill(nSigma3Pi[i]);
      registry.get<TH1>(HIST("control/cut0/h3pi1eMass"))->Fill(mass3pi1e[i]);
    } // end of loop over 4 tracks

    // control, 1 per event
    registry.get<TH1>(HIST("control/cut0/h4trkPtTot"))->Fill(pttot);
    registry.get<TH1>(HIST("control/cut0/h4piMass"))->Fill(mass4pi);
    registry.get<TH2>(HIST("control/cut0/h4trkMassVsPt"))->Fill(mass4pi, pttot);
    registry.get<TH1>(HIST("control/cut0/hNtofTrk"))->Fill(nTofTrk);
    registry.get<TH2>(HIST("control/cut0/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);

    if (pttot < 0.150) {
      // give all the gg combinations
      // 12 34
      // 13 24
      // 14 23
      int nGamPair = 0;
      int nGamma = 0;
      int whichPair = -1;
      float scalarAsym = -1;
      float vectorAsym = -1;
      bool electronCheck = true;
      for (int i = 0; i < 4; i++) {
        if (std::abs(nSigmaEl[i]) > 5)
          electronCheck = false;
      } // end of loop over 4 tracks

      for (int i = 0; i < 3; i++) {
        registry.get<TH1>(HIST("control/cut0/hInvMass2ElCoh"))->Fill(invMass2El[i][0]);
        registry.get<TH1>(HIST("control/cut0/hInvMass2ElCoh"))->Fill(invMass2El[i][1]);
        registry.get<TH1>(HIST("control/cut0/hGamPtCoh"))->Fill(gammaPair[i][0].Pt());
        registry.get<TH1>(HIST("control/cut0/hGamPtCoh"))->Fill(gammaPair[i][1].Pt());
        if (invMass2El[i][0] < 0.15) {
          nGamma++;
          registry.get<TH1>(HIST("control/cut0/hGamPtCohIM0"))->Fill(gammaPair[i][0].Pt());
        }
        if (invMass2El[i][1] < 0.15) {
          nGamma++;
          registry.get<TH1>(HIST("control/cut0/hGamPtCohIM0"))->Fill(gammaPair[i][1].Pt());
        }
        if (invMass2El[i][0] < 0.15 && invMass2El[i][1] < 0.15) {
          nGamPair++;
          whichPair = i;
        }
      } // end of loop over 3 tracks
      registry.get<TH1>(HIST("control/cut0/hN2gamma"))->Fill(nGamma);
      registry.get<TH1>(HIST("control/cut0/hN2gamma"))->Fill(10 + nGamPair);
      if (nGamPair == 1) {
        scalarAsym = std::abs(gammaPair[whichPair][1].Pt() - gammaPair[whichPair][0].Pt()) / (gammaPair[whichPair][1].Pt() + gammaPair[whichPair][0].Pt());
        if ((gammaPair[whichPair][1] - gammaPair[whichPair][0]).Pt() != 0)
          vectorAsym = (gammaPair[whichPair][1] + gammaPair[whichPair][0]).Pt() / (gammaPair[whichPair][1] - gammaPair[whichPair][0]).Pt();
        registry.get<TH1>(HIST("control/cut0/hGamAS"))->Fill(scalarAsym);
        registry.get<TH1>(HIST("control/cut0/hGamAV"))->Fill(vectorAsym);
        registry.get<TH1>(HIST("control/cut0/hInvMass2GamCoh"))->Fill((gammaPair[whichPair][1] + gammaPair[whichPair][0]).M());
        registry.get<TH1>(HIST("control/cut0/hDeltaPhi2GamCoh"))->Fill(CalculateDeltaPhi(gammaPair[whichPair][1], gammaPair[whichPair][0]));
        for (int j = 0; j < 4; j++)
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut40"))->Fill(tmpMomentum[j], tmpDedx[j]);
        if ((gammaPair[whichPair][1] + gammaPair[whichPair][0]).M() > 3. &&
            (gammaPair[whichPair][1] + gammaPair[whichPair][0]).M() < 4.) {
          for (int j = 0; j < 4; j++)
            registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut0CohPsi2s"))->Fill(tmpMomentum[j], tmpDedx[j]);
        }

        if (electronCheck) {
          registry.get<TH1>(HIST("control/cut20/hInvMass2ElCoh"))->Fill(gammaPair[whichPair][0].M());
          registry.get<TH1>(HIST("control/cut20/hInvMass2ElCoh"))->Fill(gammaPair[whichPair][1].M());
          registry.get<TH1>(HIST("control/cut20/hGamPtCohIM0"))->Fill(gammaPair[whichPair][0].Pt());
          registry.get<TH1>(HIST("control/cut20/hGamPtCohIM0"))->Fill(gammaPair[whichPair][1].Pt());
          registry.get<TH1>(HIST("control/cut20/hGamAS"))->Fill(scalarAsym);
          registry.get<TH1>(HIST("control/cut20/hGamAV"))->Fill(vectorAsym);
          registry.get<TH1>(HIST("control/cut20/hInvMass2GamCoh"))->Fill((gammaPair[whichPair][1] + gammaPair[whichPair][0]).M());
          registry.get<TH1>(HIST("control/cut20/hDeltaPhi2GamCoh"))->Fill(CalculateDeltaPhi(gammaPair[whichPair][1], gammaPair[whichPair][0]));
        }

      } // ngam = 1

    } // end of check ptot<0.15

    // remove combinatorics
    bool flagTotal[4] = {false, false, false, false};
    bool flagIM[4] = {false, false, false, false};
    bool flagDP[4] = {false, false, false, false};
    bool flagEl[4] = {false, false, false, false};
    bool flagPi[4] = {false, false, false, false};
    bool flagPt[4] = {false, false, false, false};
    bool flagPr[4] = {false, false, false, false};
    bool flagKa[4] = {false, false, false, false};
    bool flagCR[4] = {false, false, false, false};
    bool flagS3pi[4] = {false, false, false, false};

    // bool flagVcalPV[4]={false,false,false,false};
    // float deltaphi=0;

    for (int i = 0; i < 4; i++) {
      if (pi3invMass[i] < invMass3piMaxcut) { // default should be 1.8
        if (invMass3piSignalRegion) {
          flagIM[i] = true;
        } else {
          flagIM[i] = false;
        }
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut2"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (pi3deltaPhi[i] > deltaPhiMincut) { // default should be 1.6
        flagDP[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut3"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (minNsigmaElcut < nSigmaEl[i] && nSigmaEl[i] < maxNsigmaElcut) { // default (-2,3)
        flagEl[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut4"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (std::abs(nSigmaPi[i]) > maxNsigmaPiVetocut) { // default is 4
        flagPi[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut5"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (tmpPt[i] > minPtEtrkcut) { // 0.25
        flagPt[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut6"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (flagVcalPV[i]) {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut7"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (pttot < 0.15) {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut8"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (std::abs(nSigmaPr[i]) > maxNsigmaPrVetocut) { // default is 3
        flagPr[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut9"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (std::abs(nSigmaKa[i]) > maxNsigmaKaVetocut) { // default is 3
        flagKa[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut10"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (nclTPCcrossedRows[i] > nTPCcrossedRowsMinCut) { // default is 50
        flagCR[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut11"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (nSigma3Pi[i] < nSigma3piMaxCut) { // default is 5
        flagS3pi[i] = true;
      } else {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut12"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      if (flagIMGam2ePV[i]) {
        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut13"))->Fill(tmpMomentum[i], tmpDedx[i]);
      }

      flagTotal[i] = flagEl[i] && flagPi[i] && flagPt[i] && !flagVcalPV[i] && flagPr[i] && flagKa[i] && flagIM[i] && flagDP[i] && flagCR[i] && flagS3pi[i];
    } // end of loop over 4 tracks

    int counterM3pi = flagIM[0] + flagIM[1] + flagIM[2] + flagIM[3];
    int counterDphi = flagDP[0] + flagDP[1] + flagDP[2] + flagDP[3];
    int counterEl = flagEl[0] + flagEl[1] + flagEl[2] + flagEl[3];
    int counterPi = flagPi[0] + flagPi[1] + flagPi[2] + flagPi[3];
    int counterPr = flagPr[0] + flagPr[1] + flagPr[2] + flagPr[3];
    int counterKa = flagKa[0] + flagKa[1] + flagKa[2] + flagKa[3];
    int counterPt = flagPt[0] + flagPt[1] + flagPt[2] + flagPt[3];
    int counterCR = flagCR[0] + flagCR[1] + flagCR[2] + flagCR[3];
    int counterVcal = !flagVcalPV[0] + !flagVcalPV[1] + !flagVcalPV[2] + !flagVcalPV[3];
    int counterS3pi = flagS3pi[0] + flagS3pi[1] + flagS3pi[2] + flagS3pi[3];
    int counterTotal = flagTotal[0] + flagTotal[1] + flagTotal[2] + flagTotal[3];

    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(5. + counterM3pi, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(10. + counterDphi, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(15. + counterEl, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(20. + counterPi, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(25. + counterPt, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(30. + counterVcal, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(35. + counterPr, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(40. + counterKa, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(45. + counterCR, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(50. + counterS3pi, 1.);
    registry.get<TH1>(HIST("global/hNCombAfterCut"))->Fill(55. + counterTotal, 1.);

    // draw PID histograms
    if (counterEl > 0) { // Nelectrons>0, cut20
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(6., 1.);
      registry.get<TH1>(HIST("control/cut20/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut20/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut20/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut20/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut20/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut20"))->Fill(tmpMomentum[i], tmpDedx[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut20"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          FillControlHistos<20>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
          registry.get<TH2>(HIST("control/cut20/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut20/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut20/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }

      if (flagEl[0] * flagPi[0] + flagEl[1] * flagPi[1] + flagEl[2] * flagPi[2] + flagEl[3] * flagPi[3] > 0) { // pi veto, cut21
        registry.get<TH1>(HIST("global/hEventEff"))->Fill(7., 1.);
        registry.get<TH1>(HIST("control/cut21/h4trkPtTot"))->Fill(pttot);
        registry.get<TH1>(HIST("control/cut21/h4piMass"))->Fill(mass4pi);
        registry.get<TH2>(HIST("control/cut21/h4trkMassVsPt"))->Fill(mass4pi, pttot);
        registry.get<TH1>(HIST("control/cut21/hNtofTrk"))->Fill(nTofTrk);
        registry.get<TH2>(HIST("control/cut21/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
        for (int i = 0; i < 4; i++) {
          if (flagEl[i] && flagPi[i]) {
            registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut21"))->Fill(tmpMomentum[i], tmpDedx[i]);
            // for (int j = 0; j < 4; j++) {
            //   if (i == j) continue;
            //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut21"))->Fill(tmpMomentum[j], tmpDedx[j]);
            // }
            FillControlHistos<21>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
            registry.get<TH2>(HIST("control/cut21/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
            registry.get<TH1>(HIST("control/cut21/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            registry.get<TH1>(HIST("control/cut21/h3pi1eMass"))->Fill(mass3pi1e[i]);
          }
        }

        if (flagEl[0] * flagPi[0] * !flagVcalPV[0] +
              flagEl[1] * flagPi[1] * !flagVcalPV[1] +
              flagEl[2] * flagPi[2] * !flagVcalPV[2] +
              flagEl[3] * flagPi[3] * !flagVcalPV[3] >
            0) { // vcal veto, cut22
          registry.get<TH1>(HIST("global/hEventEff"))->Fill(8., 1.);
          registry.get<TH1>(HIST("control/cut22/h4trkPtTot"))->Fill(pttot);
          registry.get<TH1>(HIST("control/cut22/h4piMass"))->Fill(mass4pi);
          registry.get<TH2>(HIST("control/cut22/h4trkMassVsPt"))->Fill(mass4pi, pttot);
          registry.get<TH1>(HIST("control/cut22/hNtofTrk"))->Fill(nTofTrk);
          registry.get<TH2>(HIST("control/cut22/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
          for (int i = 0; i < 4; i++) {
            if (flagEl[i] && flagPi[i] && !flagVcalPV[i]) {
              registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut22"))->Fill(tmpMomentum[i], tmpDedx[i]);
              // for (int j = 0; j < 4; j++) {
              //  if (i == j) continue;
              //  registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut22"))->Fill(tmpMomentum[j], tmpDedx[j]);
              // }
              FillControlHistos<22>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
              registry.get<TH2>(HIST("control/cut22/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registry.get<TH1>(HIST("control/cut22/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry.get<TH1>(HIST("control/cut22/h3pi1eMass"))->Fill(mass3pi1e[i]);
            }
          }

          if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] +
                flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] +
                flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] +
                flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] >
              0) { // pT veto, cut23
            registry.get<TH1>(HIST("global/hEventEff"))->Fill(9., 1.);
            registry.get<TH1>(HIST("control/cut23/h4trkPtTot"))->Fill(pttot);
            registry.get<TH1>(HIST("control/cut23/h4piMass"))->Fill(mass4pi);
            registry.get<TH2>(HIST("control/cut23/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registry.get<TH1>(HIST("control/cut23/hNtofTrk"))->Fill(nTofTrk);
            registry.get<TH2>(HIST("control/cut23/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
            for (int i = 0; i < 4; i++) {
              if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i]) {
                registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut23"))->Fill(tmpMomentum[i], tmpDedx[i]);
                // for (int j = 0; j < 4; j++) {
                //   if (i == j) continue;
                //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut23"))->Fill(tmpMomentum[j], tmpDedx[j]);
                // }
                FillControlHistos<23>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                registry.get<TH2>(HIST("control/cut23/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                registry.get<TH1>(HIST("control/cut23/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                registry.get<TH1>(HIST("control/cut23/h3pi1eMass"))->Fill(mass3pi1e[i]);
              }
            }

            if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] +
                  flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] +
                  flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] +
                  flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] >
                0) { // proton veto, cut24
              registry.get<TH1>(HIST("global/hEventEff"))->Fill(10., 1.);
              registry.get<TH1>(HIST("control/cut24/h4trkPtTot"))->Fill(pttot);
              registry.get<TH1>(HIST("control/cut24/h4piMass"))->Fill(mass4pi);
              registry.get<TH2>(HIST("control/cut24/h4trkMassVsPt"))->Fill(mass4pi, pttot);
              registry.get<TH1>(HIST("control/cut24/hNtofTrk"))->Fill(nTofTrk);
              registry.get<TH2>(HIST("control/cut24/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
              for (int i = 0; i < 4; i++) {
                if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i]) {
                  registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut24"))->Fill(tmpMomentum[i], tmpDedx[i]);
                  // for (int j = 0; j < 4; j++) {
                  //   if (i == j) continue;
                  //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut24"))->Fill(tmpMomentum[j], tmpDedx[j]);
                  // }
                  FillControlHistos<24>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                  registry.get<TH2>(HIST("control/cut24/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                  registry.get<TH1>(HIST("control/cut24/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                  registry.get<TH1>(HIST("control/cut24/h3pi1eMass"))->Fill(mass3pi1e[i]);
                }
              }

              if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] * flagKa[0] +
                    flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] * flagKa[1] +
                    flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] * flagKa[2] +
                    flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] * flagKa[3] >
                  0) { // kaon veto, cut25
                registry.get<TH1>(HIST("global/hEventEff"))->Fill(11., 1.);
                registry.get<TH1>(HIST("control/cut25/h4trkPtTot"))->Fill(pttot);
                registry.get<TH1>(HIST("control/cut25/h4piMass"))->Fill(mass4pi);
                registry.get<TH2>(HIST("control/cut25/h4trkMassVsPt"))->Fill(mass4pi, pttot);
                registry.get<TH1>(HIST("control/cut25/hNtofTrk"))->Fill(nTofTrk);
                registry.get<TH2>(HIST("control/cut25/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
                for (int i = 0; i < 4; i++) {
                  if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i] && flagKa[i]) {
                    registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut25"))->Fill(tmpMomentum[i], tmpDedx[i]);
                    // for (int j = 0; j < 4; j++) {
                    //   if (i == j) continue;
                    //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut25"))->Fill(tmpMomentum[j], tmpDedx[j]);
                    // }
                    FillControlHistos<25>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                    registry.get<TH2>(HIST("control/cut25/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                    registry.get<TH1>(HIST("control/cut25/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                    registry.get<TH1>(HIST("control/cut25/h3pi1eMass"))->Fill(mass3pi1e[i]);
                  }
                }

                if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] * flagKa[0] * flagIM[0] +
                      flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] * flagKa[1] * flagIM[1] +
                      flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] * flagKa[2] * flagIM[2] +
                      flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] * flagKa[3] * flagIM[3] >
                    0) { // 3pi cut, cut26
                  registry.get<TH1>(HIST("global/hEventEff"))->Fill(12., 1.);
                  registry.get<TH1>(HIST("control/cut26/h4trkPtTot"))->Fill(pttot);
                  registry.get<TH1>(HIST("control/cut26/h4piMass"))->Fill(mass4pi);
                  registry.get<TH2>(HIST("control/cut26/h4trkMassVsPt"))->Fill(mass4pi, pttot);
                  registry.get<TH1>(HIST("control/cut26/hNtofTrk"))->Fill(nTofTrk);
                  registry.get<TH2>(HIST("control/cut26/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
                  for (int i = 0; i < 4; i++) {
                    if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i] && flagKa[i] && flagIM[i]) {
                      registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut26"))->Fill(tmpMomentum[i], tmpDedx[i]);
                      // for (int j = 0; j < 4; j++) {
                      //   if (i == j) continue;
                      //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut26"))->Fill(tmpMomentum[j], tmpDedx[j]);
                      // }
                      FillControlHistos<26>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                      registry.get<TH2>(HIST("control/cut26/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                      registry.get<TH1>(HIST("control/cut26/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                      registry.get<TH1>(HIST("control/cut26/h3pi1eMass"))->Fill(mass3pi1e[i]);
                    }
                  }

                  if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] * flagKa[0] * flagIM[0] * flagDP[0] +
                        flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] * flagKa[1] * flagIM[1] * flagDP[1] +
                        flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] * flagKa[2] * flagIM[2] * flagDP[2] +
                        flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] * flagKa[3] * flagIM[3] * flagDP[3] >
                      0) { // delta phi cut, cut27
                    registry.get<TH1>(HIST("global/hEventEff"))->Fill(13., 1.);
                    registry.get<TH1>(HIST("control/cut27/h4trkPtTot"))->Fill(pttot);
                    registry.get<TH1>(HIST("control/cut27/h4piMass"))->Fill(mass4pi);
                    registry.get<TH2>(HIST("control/cut27/h4trkMassVsPt"))->Fill(mass4pi, pttot);
                    registry.get<TH1>(HIST("control/cut27/hNtofTrk"))->Fill(nTofTrk);
                    registry.get<TH2>(HIST("control/cut27/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
                    for (int i = 0; i < 4; i++) {
                      if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i] && flagKa[i] && flagIM[i] && flagDP[i]) {
                        registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut27"))->Fill(tmpMomentum[i], tmpDedx[i]);
                        // for (int j = 0; j < 4; j++) {
                        //   if (i == j) continue;
                        //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut27"))->Fill(tmpMomentum[j], tmpDedx[j]);
                        // }
                        FillControlHistos<27>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                        registry.get<TH2>(HIST("control/cut27/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                        registry.get<TH1>(HIST("control/cut27/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                        registry.get<TH1>(HIST("control/cut27/h3pi1eMass"))->Fill(mass3pi1e[i]);
                      }
                    }

                    if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] * flagKa[0] * flagIM[0] * flagDP[0] * flagCR[0] +
                          flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] * flagKa[1] * flagIM[1] * flagDP[1] * flagCR[1] +
                          flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] * flagKa[2] * flagIM[2] * flagDP[2] * flagCR[2] +
                          flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] * flagKa[3] * flagIM[3] * flagDP[3] * flagCR[3] >
                        0) { // Nc-rTPC cut, cut28
                      registry.get<TH1>(HIST("global/hEventEff"))->Fill(14., 1.);
                      registry.get<TH1>(HIST("control/cut28/h4trkPtTot"))->Fill(pttot);
                      registry.get<TH1>(HIST("control/cut28/h4piMass"))->Fill(mass4pi);
                      registry.get<TH2>(HIST("control/cut28/h4trkMassVsPt"))->Fill(mass4pi, pttot);
                      registry.get<TH1>(HIST("control/cut28/hNtofTrk"))->Fill(nTofTrk);
                      registry.get<TH2>(HIST("control/cut28/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
                      for (int i = 0; i < 4; i++) {
                        if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i] && flagKa[i] && flagIM[i] && flagDP[i] && flagCR[i]) {
                          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut28"))->Fill(tmpMomentum[i], tmpDedx[i]);
                          FillControlHistos<28>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                          registry.get<TH2>(HIST("control/cut28/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                          registry.get<TH1>(HIST("control/cut28/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                          registry.get<TH1>(HIST("control/cut28/h3pi1eMass"))->Fill(mass3pi1e[i]);
                        }
                      }

                      if (flagEl[0] * flagPi[0] * !flagVcalPV[0] * flagPt[0] * flagPr[0] * flagKa[0] * flagIM[0] * flagDP[0] * flagCR[0] * flagS3pi[0] +
                            flagEl[1] * flagPi[1] * !flagVcalPV[1] * flagPt[1] * flagPr[1] * flagKa[1] * flagIM[1] * flagDP[1] * flagCR[1] * flagS3pi[1] +
                            flagEl[2] * flagPi[2] * !flagVcalPV[2] * flagPt[2] * flagPr[2] * flagKa[2] * flagIM[2] * flagDP[2] * flagCR[2] * flagS3pi[2] +
                            flagEl[3] * flagPi[3] * !flagVcalPV[3] * flagPt[3] * flagPr[3] * flagKa[3] * flagIM[3] * flagDP[3] * flagCR[3] * flagS3pi[3] >
                          0) { // nsigma 3pi cut, cut29
                        registry.get<TH1>(HIST("global/hEventEff"))->Fill(15., 1.);
                        registry.get<TH1>(HIST("control/cut29/h4trkPtTot"))->Fill(pttot);
                        registry.get<TH1>(HIST("control/cut29/h4piMass"))->Fill(mass4pi);
                        registry.get<TH2>(HIST("control/cut29/h4trkMassVsPt"))->Fill(mass4pi, pttot);
                        registry.get<TH1>(HIST("control/cut29/hNtofTrk"))->Fill(nTofTrk);
                        registry.get<TH2>(HIST("control/cut29/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
                        for (int i = 0; i < 4; i++) {
                          if (flagEl[i] && flagPi[i] && !flagVcalPV[i] && flagPt[i] && flagPr[i] && flagKa[i] && flagIM[i] && flagDP[i] && flagCR[i] && flagS3pi[i]) {
                            registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut29"))->Fill(tmpMomentum[i], tmpDedx[i]);
                            FillControlHistos<29>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
                            registry.get<TH2>(HIST("control/cut29/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
                            registry.get<TH1>(HIST("control/cut29/hsigma3Pi"))->Fill(nSigma3Pi[i]);
                            registry.get<TH1>(HIST("control/cut29/h3pi1eMass"))->Fill(mass3pi1e[i]);
                            if (verbose) {
                              LOGF(info, "cut29 timeTot %f, resTot %f, trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTimeTot, trkTimeResTot, trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
                            }
                          }
                        }

                      } else {
                        if (verbose) {
                          LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID+Dphi+IM+CR");
                        }
                      } // end of nsigma 3pi cut
                    } else {
                      if (verbose) {
                        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID+Dphi+IM+CR");
                      }
                    } // end of TPC crossed rows for electron cut
                  } else {
                    if (verbose) {
                      LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID+Dphi+IM");
                    }
                  } // end of delta phi cut
                } else {
                  if (verbose) {
                    LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID+Dphi");
                  }
                } // end of inv mass 3 pi cut
              } else {
                if (verbose) {
                  LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID");
                }
              } // end of kaon veto
            } else {
              if (verbose) {
                LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID");
              }
            } // end of proton veto
          } else {
            if (verbose) {
              LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT");
            }
          } // end of pT veto
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal");
          }
        } // end of vcal veto
      } else {
        if (verbose) {
          LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by pi PID");
        }
      } // end of pi veto
    } else { // no electron
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: no electron PID among 4 tracks");
      }
    } // end of Nelectrons check

    // skip events with pttot<0.15
    if (pttot < ptTotcut) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: pt tot is %f", pttot);
      }
      return;
    }
    if (counterTotal > 0) {
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(16., 1.);
      registry.get<TH1>(HIST("control/cut30/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut30/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut30/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut30/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut30/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
      for (int i = 0; i < 4; i++) {
        if (flagTotal[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut30"))->Fill(tmpMomentum[i], tmpDedx[i]);
          FillControlHistos<30>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
          registry.get<TH2>(HIST("control/cut30/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut30/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut30/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    }

    // check FIT information
    if (FITvetoFlag) {
      auto bitMin = 16 - FITvetoWindow; // default is +- 1 bc (1 bit)
      auto bitMax = 16 + FITvetoWindow;
      for (auto bit = bitMin; bit <= bitMax; bit++) {
        if (TESTBIT(dgcand.bbFT0Apf(), bit))
          return;
        if (TESTBIT(dgcand.bbFT0Cpf(), bit))
          return;
        if (useFV0ForVeto && TESTBIT(dgcand.bbFV0Apf(), bit))
          return;
        if (useFDDAForVeto && TESTBIT(dgcand.bbFDDApf(), bit))
          return;
        if (useFDDCForVeto && TESTBIT(dgcand.bbFDDCpf(), bit))
          return;
      } // end of loop over bits
    } // end of check emptyness around given BC in FIT detectors
    if (counterTotal > 0) {
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(17., 1.);
      registry.get<TH1>(HIST("control/cut31/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut31/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut31/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut31/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut31/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
      for (int i = 0; i < 4; i++) {
        if (flagTotal[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut31"))->Fill(tmpMomentum[i], tmpDedx[i]);
          FillControlHistos<31>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
          registry.get<TH2>(HIST("control/cut31/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut31/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut31/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    }

    // n TOF tracks cut
    if (nTofTrk < nTofTrkMinCut) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: nTOFtracks is %d", nTofTrk);
      }
      return;
    }
    if (counterTotal > 0) {
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(18., 1.);
      registry.get<TH1>(HIST("control/cut32/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut32/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut32/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut32/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut32/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
      for (int i = 0; i < 4; i++) {
        if (flagTotal[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut32"))->Fill(tmpMomentum[i], tmpDedx[i]);
          FillControlHistos<32>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
          registry.get<TH2>(HIST("control/cut32/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut32/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut32/h3pi1eMass"))->Fill(mass3pi1e[i]);
          registry.get<TH1>(HIST("control/cut32/hPtSpectrumEl"))->Fill(tmpPt[i]);
          if (verbose) {
            LOGF(info, "cut32 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
          }
        }
      }
    }

    // electron has TOF hit

    // only 1 electron
    if (counterTotal == 1) {
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(19., 1.);
      for (int i = 0; i < 4; i++) {
        registry.get<TH1>(HIST("control/cut1/hDcaZ"))->Fill(dcaZ[i]);
        registry.get<TH1>(HIST("control/cut1/hDcaXY"))->Fill(dcaXY[i]);
        registry.get<TH1>(HIST("control/cut1/hChi2TPC"))->Fill(chi2TPC[i]);
        registry.get<TH1>(HIST("control/cut1/hChi2ITS"))->Fill(chi2ITS[i]);

        if (flagTotal[i]) {
          FillControlHistos<1>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
          registry.get<TH2>(HIST("control/cut1/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut1"))->Fill(tmpMomentum[i], tmpDedx[i]);
          for (int j = 0; j < 4; j++) {
            if (i == j)
              continue;
            registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut1"))->Fill(tmpMomentum[j], tmpDedx[j]);
          }
          registry.get<TH1>(HIST("global/hFinalPtSpectrumEl"))->Fill(tmpPt[i]);
          registry.get<TH1>(HIST("control/cut1/hTPCnclsFindable"))->Fill(nclTPCfind[i]);
          registry.get<TH1>(HIST("control/cut1/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut1/h3pi1eMass"))->Fill(mass3pi1e[i]);
          if (verbose) {
            LOGF(info, "cut1 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
          }

          // one electron with tof hit (cut33)
          if (tmpHasTOF[i]) {
            registry.get<TH1>(HIST("global/hEventEff"))->Fill(20., 1.);
            // registry.get<TH1>(HIST("control/cut33/hDcaZ"))->Fill(dcaZ[i]);
            // registry.get<TH1>(HIST("control/cut33/hDcaXY"))->Fill(dcaXY[i]);
            // registry.get<TH1>(HIST("control/cut33/hChi2TPC"))->Fill(chi2TPC[i]);
            // registry.get<TH1>(HIST("control/cut33/hChi2ITS"))->Fill(chi2ITS[i]);
            FillControlHistos<33>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
            registry.get<TH2>(HIST("control/cut33/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
            registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut33"))->Fill(tmpMomentum[i], tmpDedx[i]);
            registry.get<TH1>(HIST("control/cut33/h4trkPtTot"))->Fill(pttot);
            registry.get<TH1>(HIST("control/cut33/h4piMass"))->Fill(mass4pi);
            registry.get<TH2>(HIST("control/cut33/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registry.get<TH1>(HIST("control/cut33/hNtofTrk"))->Fill(nTofTrk);
            registry.get<TH2>(HIST("control/cut33/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
            registry.get<TH1>(HIST("control/cut33/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            registry.get<TH1>(HIST("control/cut33/h3pi1eMass"))->Fill(mass3pi1e[i]);
            registry.get<TH1>(HIST("control/cut33/hPtSpectrumEl"))->Fill(tmpPt[i]);
            if (verbose) {
              LOGF(info, "cut33 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
            }

            int otherTOFtracks = 0;
            for (int j = 0; j < 4; j++) {
              if (i == j)
                continue;
              if (tmpHasTOF[j])
                otherTOFtracks++;
            }
            // at least one pion with tof hit (cut34)
            if (otherTOFtracks >= 1) {
              registry.get<TH1>(HIST("global/hEventEff"))->Fill(21., 1.);
              // registry.get<TH1>(HIST("control/cut34/hDcaZ"))->Fill(dcaZ[i]);
              // registry.get<TH1>(HIST("control/cut34/hDcaXY"))->Fill(dcaXY[i]);
              // registry.get<TH1>(HIST("control/cut34/hChi2TPC"))->Fill(chi2TPC[i]);
              // registry.get<TH1>(HIST("control/cut34/hChi2ITS"))->Fill(chi2ITS[i]);
              FillControlHistos<34>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i]);
              registry.get<TH2>(HIST("control/cut34/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut34"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registry.get<TH1>(HIST("control/cut34/h4trkPtTot"))->Fill(pttot);
              registry.get<TH1>(HIST("control/cut34/h4piMass"))->Fill(mass4pi);
              registry.get<TH2>(HIST("control/cut34/h4trkMassVsPt"))->Fill(mass4pi, pttot);
              registry.get<TH1>(HIST("control/cut34/hNtofTrk"))->Fill(nTofTrk);
              registry.get<TH2>(HIST("control/cut34/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
              registry.get<TH1>(HIST("control/cut34/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry.get<TH1>(HIST("control/cut34/h3pi1eMass"))->Fill(mass3pi1e[i]);
              registry.get<TH1>(HIST("control/cut34/hPtSpectrumEl"))->Fill(tmpPt[i]);
              if (verbose) {
                LOGF(info, "cut34 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
              }
            } // end of at least one pion with tof hit (cut34)

          } // end of one electron with tof hit (cut33)
        }
      } // end of loop over 4 tracks
      registry.get<TH1>(HIST("control/cut1/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut1/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut1/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut1/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut1/hZNACenergy"))->Fill(ZNAenergy, ZNCenergy);
      // special case invmass 4pi (2,2.3)
      // if (mass4pi<2.3 && mass4pi>2) {
      // for (int i = 0; i < 4; i++) {
      //   registry.get<TH1>(HIST("control/cut1/cut1a/hDcaZ"))->Fill(dcaZ[i]);
      //   registry.get<TH1>(HIST("control/cut1/cut1a/hDcaXY"))->Fill(dcaXY[i]);
      //   registry.get<TH1>(HIST("control/cut1/cut1a/hChi2TPC"))->Fill(chi2TPC[i]);
      //   registry.get<TH1>(HIST("control/cut1/cut1a/hChi2ITS"))->Fill(chi2ITS[i]);
      //
      //   if (flagTotal[i]) {
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h3piMassComb"))->Fill(pi3invMass[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h3trkPtTot"))->Fill(pi3pt[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/hDeltaPhi13topo"))->Fill(pi3deltaPhi[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h13AssymPt1ProngAver"))->Fill(pi3assymav[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h13Vector"))->Fill(pi3vector[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h13Scalar"))->Fill(pi3scalar[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/h13EtaSum"))->Fill(pi3etasum[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/hTPCnCrossedRows"))->Fill(nclTPCcrossedRows[i]);
      //
      //     registry.get<TH2>(HIST("control/cut1/cut1a/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/hTPCnclsFindable"))->Fill(nclTPCfind[i]);
      //     registry.get<TH1>(HIST("control/cut1/cut1a/hsigma3Pi"))->Fill(nSigma3Pi[i]);
      //   }
      // }
      // registry.get<TH1>(HIST("control/cut1/cut1a/h4trkPtTot"))->Fill(pttot);
      // registry.get<TH1>(HIST("control/cut1/cut1a/h4piMass"))->Fill(mass4pi);
      // registry.get<TH2>(HIST("control/cut1/cut1a/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      // registry.get<TH1>(HIST("control/cut1/cut1a/hNtofTrk"))->Fill(nTofTrk);
      // } // end of mass window for 4pi case

    } else { // more than 1 electron candidate
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: more than one electron candidate");
      }
    } // end of 1electrons check
  } // end of process
  // check ntracks-4PVtracks
  // check pt of remaining (ntracks-4PVtracks) tracks
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauTau13topo>(cfgc, TaskName{"TauTau13topo"})};
}
