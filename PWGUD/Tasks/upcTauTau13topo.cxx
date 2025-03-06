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
/// \file upcTauTau13topo.cxx
/// \brief tau tau analysis 1e+3pi topology
/// \author Adam Matyja, adam.tomasz.matyja@cern.ch, adam.matryja@ifj.edu.pl
/// \since  January 2024
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

#include "Common/Core/RecoDecay.h"
// #include <CommonUtils/EnumFlags.h>
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct TauTau13topo {
  SGSelector sgSelector;
  // configurables
  Configurable<float> cutFV0{"cutFV0", 10000., "FV0A threshold"};
  Configurable<float> cutFT0A{"cutFT0A", 150., "FT0A threshold"};
  Configurable<float> cutFT0C{"cutFT0C", 50., "FT0C threshold"};
  Configurable<float> cutZDC{"cutZDC", 10000., "ZDC threshold"};
  Configurable<float> mGapSide{"mGapSide", 2, "gap selection"};
  //  ConfigurableAxis ptAxis{"pAxis", {100, 0., 5.}, "#it{p} (GeV/#it{c})"};
  ConfigurableAxis ptAxis{"ptAxis", {120, 0., 4.}, "#it{p} (GeV/#it{c})"};
  //  ConfigurableAxis etaAxis{"etaAxis", {100, -2., 2.}, "#eta"};
  ConfigurableAxis dedxAxis{"dedxAxis", {100, 20., 160.}, "dE/dx"};
  ConfigurableAxis minvAxis{"minvAxis", {100, 0.4, 3.5}, "M_{inv} (GeV/#it{c}^{2})"};
  ConfigurableAxis phiAxis{"phiAxis", {120, 0., 3.2}, "#phi"};
  //  ConfigurableAxis vectorAxis{"vectorAxis", {100, 0., 2.}, "A_{V}"};
  //  ConfigurableAxis scalarAxis{"scalarAxis", {100, -1., 1.}, "A_{S}"};
  Configurable<bool> verbose{"verbose", {}, "Additional print outs"};

  // cut selection configurables
  Configurable<float> zvertexcut{"zvertexcut", 10., "Z vertex cut"};
  Configurable<float> trkEtacut{"trkEtacut", 0.9, "max track eta cut"};
  Configurable<bool> sameSign{"sameSign", {}, "Switch: same(true) or opposite(false) sign"};
  Configurable<float> ptTotcut{"ptTotcut", 0.15, "min pt of all 4 tracks cut"};
  Configurable<float> minAnglecut{"minAnglecut", 0.05, "min angle between tracks cut"};
  Configurable<float> minNsigmaElcut{"minNsigmaElcut", -2., "min Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaElcut{"maxNsigmaElcut", 3., "max Nsigma for Electrons cut"};
  Configurable<float> maxNsigmaPiVetocut{"maxNsigmaPiVetocut", 4., "max Nsigma for Pion veto cut"};
  Configurable<float> maxNsigmaPrVetocut{"maxNsigmaPrVetocut", 3., "max Nsigma for Proton veto cut"};
  Configurable<float> maxNsigmaKaVetocut{"maxNsigmaKaVetocut", 3., "max Nsigma for Kaon veto cut"};
  Configurable<float> minPtEtrkcut{"minPtEtrkcut", 0.25, "min Pt for El track cut"};
  Configurable<bool> mFITvetoFlag{"mFITvetoFlag", {}, "To apply FIT veto"};
  Configurable<int> mFITvetoWindow{"mFITvetoWindow", 6, "FIT veto window"};
  Configurable<bool> useFV0ForVeto{"useFV0ForVeto", 0, "use FV0 for veto"};
  Configurable<bool> useFDDAForVeto{"useFDDAForVeto", 0, "use FDDA for veto"};
  Configurable<bool> useFDDCForVeto{"useFDDCForVeto", 0, "use FDDC for veto"};
  Configurable<int> nTofTrkMinCut{"nTofTrkMinCut", 1, "min TOF tracks"};

  Configurable<bool> invMass3piSignalRegion{"invMass3piSignalRegion", 1, "1-use inv mass 3pi in signal region, 0-in background region"};
  Configurable<float> invMass3piMaxcut{"invMass3piMaxcut", 1.8, "Z invariant mass of 3 pi cut"};
  Configurable<float> deltaPhiMincut{"deltaPhiMincut", 0., "delta phi electron - 3 pi direction cut"};
  Configurable<int> nTPCcrossedRowsMinCut{"nTPCcrossedRowsMinCut", 50, "min N_crossed TPC rows for electron candidate"};
  Configurable<float> nSigma3piMaxCut{"nSigma3piMaxCut", 5., "n sigma 3 pi max cut"};

  Configurable<int> generatorIDMC{"generatorIDMC", -1, "MC generator ID"};

  // Configurable<bool> DGactive{"DGactive", false, "Switch on DGproducer"};
  // Configurable<bool> SGactive{"SGactive", true, "Switch on SGproducer"};

  // a pdg object
  // TDatabasePDG* pdg = nullptr; //not recommended
  // Service<o2::framework::O2DatabasePDG> pdg;

  // initialize histogram registry
  HistogramRegistry registry{
    "registry",
    {}};

  HistogramRegistry registryMC{
    "registryMC",
    {}};
  HistogramRegistry registry1MC{
    "registry1MC",
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
    const AxisSpec axisInvMass4trk{160, 0.5, 8.5, "#it{M}^{4trk}_{inv} (GeV/#it{c}^{2})"};

    registry.add("global/RunNumber", "Run number; Run; Collisions", {HistType::kTH1F, {{150, 544013, 545367}}});
    registry.add("global/GapSide", "Associated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});
    registry.add("global/GapSideTrue", "Recalculated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});

    registry.add("global/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry.add("global/hZNACtime", "ZNA vs ZNC time; #it{time}_{ZNA} (ns); #it{time}_{ZNC} (ns); Collisions", {HistType::kTH2F, {{100, -10., 10.}, {100, -10., 10.}}});
    // registry.add("global/hZNACenergyTest", "ZNA or ZNC energy; #it{E}_{ZNA,ZNC} (GeV); Collisions", {HistType::kTH1F, {{100,-1000,0}}});

    registry.add("global/hVertexXY", "Vertex position in x and y direction; #it{V}_{x} (cm); #it{V}_{y} (cm); Collisions", {HistType::kTH2F, {{50, -0.05, 0.05}, {50, -0.02, 0.02}}});
    registry.add("global/hVertexZ", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
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
    registry.add("global/hEventEff", "Event cut efficiency: 0-All,1-PV=4,2-Qtot=0,3-El;Cut;entries", {HistType::kTH1F, {{27, -2., 25.}}});
    registry.add("global/hNCombAfterCut", "Combinations after cut: 0-All,5-M3pi,10-Dphi,15-N_{e},20-N_{v#pi},25-Pt,30-Vcal,35-N_{vp},40-N_{vK},45-Tot;N_{comb};entries", {HistType::kTH1F, {{60, 0., 60.}}});
    // registry.add("global/hInvMassElTrack", ";M_{inv}^{2};entries", {HistType::kTH1F, {{100, -0.01, 0.49}}});
    registry.add("global/hDeltaAngleTrackPV", ";#Delta#alpha;entries", {HistType::kTH1F, {{136, 0., 3.2}}}); // 0.49
    registry.add("global/hTrkCheck", ";track type;entries", {HistType::kTH1F, {{16, -1, 15}}});

    registry.add("global/hRecFlag", ";Reconstruction Flag;events", {HistType::kTH1F, {{10, 0., 10.}}});
    registry.add("global/hOccupancyInTime", ";Occupancy;events", {HistType::kTH1F, {{100, 0., 10000.}}});

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
    registry.add("control/cut0/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut0/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    registry.add("control/cut0/hInvMass2ElAll", "Inv Mass of 2 Electrons from coherent peak;M_{inv}^{2e};entries", {HistType::kTH1F, {{150, -0.1, 9.}}});
    registry.add("control/cut0/hInvMass2ElCoh", "Inv Mass of 2 Electrons from coherent peak;M_{inv}^{2e};entries", {HistType::kTH1F, {{150, -0.1, 4.}}});
    registry.add("control/cut0/hGamPtCoh", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut0/hGamPtCohIM0", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registry.add("control/cut0/hN2gamma", "Number of gamma pairs among 3 comb;N_{#gamma#gamma};entries", {HistType::kTH1F, {{20, 0., 20.}}});
    registry.add("control/cut0/hGamAS", ";A_{S};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut0/hGamAV", ";A_{V};entries", {HistType::kTH1F, {{100, 0, 0.2}}});
    registry.add("control/cut0/hInvMass2GamCoh", "Inv Mass of 2 Gamma from coherent peak;M_{inv}^{2#gamma};entries", {HistType::kTH1F, {{160, 0.5, 4.5}}});
    registry.add("control/cut0/hDeltaPhi2GamCoh", "Delta Phi of 2 Gamma from coherent peak;#Delta#phi^{2#gamma};entries", {HistType::kTH1F, {phiAxis}});

    // // cut1
    // registry.add("control/cut1/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    // registry.add("control/cut1/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    // registry.add("control/cut1/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    // registry.add("control/cut1/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    // registry.add("control/cut1/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    // registry.add("control/cut1/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    // //    registry.add("control/cut1/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    // registry.add("control/cut1/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    // registry.add("control/cut1/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registry.add("control/cut1/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registry.add("control/cut1/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    // registry.add("control/cut1/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    // registry.add("control/cut1/hDcaZ", "All 4 tracks dca ;dca_{Z};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    // registry.add("control/cut1/hDcaXY", "All 4 tracks dca ;dca_{XY};entries", {HistType::kTH1F, {{100, -0.05, 0.05}}});
    // registry.add("control/cut1/hChi2TPC", "All 4 tracks Chi2 ;Chi2_{TPC};entries", {HistType::kTH1F, {{48, -2, 10.}}});
    // registry.add("control/cut1/hChi2ITS", "All 4 tracks Chi2 ;Chi2_{ITS};entries", {HistType::kTH1F, {{44, -2, 20.}}});
    // registry.add("control/cut1/hChi2TOF", "All 4 tracks Chi2 ;Chi2_{TOF};entries", {HistType::kTH1F, {{48, -2, 10.}}});
    // registry.add("control/cut1/hTPCnclsFindable", "All 4 tracks NclFind ;N_{TPC,cl,findable};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    // registry.add("control/cut1/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registry.add("control/cut1/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registry.add("control/cut1/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    // registry.add("control/cut1/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    // registry.add("control/cut1/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    // registry.add("control/cut1/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut20/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut20/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut21/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut21/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut22/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut22/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //     // cut23
    //     registry.add("control/cut23/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //     registry.add("control/cut23/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //     registry.add("control/cut23/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //     registry.add("control/cut23/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //     registry.add("control/cut23/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //     registry.add("control/cut23/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //     //    registry.add("control/cut23/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    //     registry.add("control/cut23/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //     registry.add("control/cut23/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //     registry.add("control/cut23/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //     registry.add("control/cut23/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //     registry.add("control/cut23/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //     registry.add("control/cut23/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //     registry.add("control/cut23/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //     registry.add("control/cut23/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //     registry.add("control/cut23/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //     registry.add("control/cut23/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    //     registry.add("control/cut23/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut24/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut24/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut25/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut25/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut26/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut26/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut27/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut27/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut28/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut28/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut29/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut29/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut30/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    registry.add("control/cut30/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //    // cut31
    //    registry.add("control/cut31/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registry.add("control/cut31/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registry.add("control/cut31/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registry.add("control/cut31/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registry.add("control/cut31/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registry.add("control/cut31/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    //    registry.add("control/cut31/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    //    registry.add("control/cut31/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //    registry.add("control/cut31/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut31/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut31/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //    registry.add("control/cut31/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //    registry.add("control/cut31/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut31/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //    registry.add("control/cut31/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registry.add("control/cut31/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //    registry.add("control/cut31/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    //    registry.add("control/cut31/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //    // cut32
    //    registry.add("control/cut32/h3piMassComb", "3#pi mass, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registry.add("control/cut32/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registry.add("control/cut32/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registry.add("control/cut32/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registry.add("control/cut32/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registry.add("control/cut32/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    //    registry.add("control/cut32/h13EtaSum", ";#eta^{1-prong}+#eta^{3-prong};entries", {HistType::kTH1F, {{100, -4., 4.}}});
    //    registry.add("control/cut32/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //    registry.add("control/cut32/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut32/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut32/h3piMassVsPt", "3#pi mass vs Pt, 1 entry per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //    registry.add("control/cut32/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //    registry.add("control/cut32/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry.add("control/cut32/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //    registry.add("control/cut32/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registry.add("control/cut32/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //    registry.add("control/cut32/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {{40, 0., 5.}}});
    //    registry.add("control/cut32/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
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
    registry.add("control/cut33/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

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
    registry.add("control/cut34/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // pid El
    registry.add("pidTPC/hpvsdedxElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registry.add("pidTPC/hpvsdedxElHipCut0CohPsi2s", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // registry.add("pidTPC/hpvsdedxElHipCut1", "All hip;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
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
    //     registry.add("pidTPC3pi/hpvsdedxElHipCut1", "All hip;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
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

    // tof histos
    registry.add("pidTOF/hpvsNsigmaElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut20", "El hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut33", "eTOF+20 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut21", "vPi+33 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut24", "vPr+21 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut28", "CR+25 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut22", "vVc+28 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut29", "s3pi+22 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut26", "IM+29 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut34", "piTOF+26 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut30", "ptTot+34 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry.add("pidTOF/hpvsNsigmaElHipCut27", "DP+30 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    // registry.add("pidTOF/hpvsNsigmaElHipCut35", "ZDC+27 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100,-5.,5.}}});

    registry.add("pidTOF/h3piTOFchi2", "tof chi2;chi2 TOF;events", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("pidTOF/h3piTOFchi2Cut34", "tof chi2;chi2 TOF;events", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("pidTOF/h3piTOFchi2Cut30", "tof chi2;chi2 TOF;events", {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("pidTOF/h3piTOFchi2Cut27", "tof chi2;chi2 TOF;events", {HistType::kTH1F, {{100, 0., 10.}}});

    // MC part
    // histograms filled by processSimpleMCSG
    // CollisionMC histograms
    registry1MC.add("globalMC/hGeneratorID", ";Generator ID;events", {HistType::kTH1F, {{100, 0., 1000.}}});
    registryMC.add("globalMC/hMCZvertex", ";V_{Z}^{MC} (cm);events", {HistType::kTH1F, {{100, -25., 25.}}});
    registryMC.add("globalMC/hMCefficiency", ";Cut Number;events", {HistType::kTH1F, {{12, 0., 12.}}});
    registryMC.add("globalMC/hMCnPart", ";N_{part};Type;events", {HistType::kTH2F, {{25, 0., 25.}, {10, 0, 10}}});
    registryMC.add("globalMC/hMCetaGen", ";#eta^{gen};N^{MC particles}", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("globalMC/hMCphiGen", ";#phi^{gen};N^{MC particles}", {HistType::kTH1F, {{100, 0., 6.4}}});
    registryMC.add("globalMC/hMCyGen", ";y^{gen};N^{MC particles}", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("globalMC/hMCptGen", ";p_{T}^{gen};N^{MC particles}", {HistType::kTH1F, {{100, 0., 4.}}});

    // MC reconstructed with information from MC true
    registryMC.add("globalMCrec/hMCetaGenCol", ";#eta^{genCol};events", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("globalMCrec/hMCphiGenCol", ";#phi^{genCol};events", {HistType::kTH1F, {{100, 0., 6.4}}});
    registryMC.add("globalMCrec/hMCyGenCol", ";y^{genCol};events", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("globalMCrec/hMCptGenCol", ";p_{T}^{genCol};events", {HistType::kTH1F, {{100, 0., 4.}}});

    registryMC.add("globalMCrec/GapSide", "Associated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});
    registryMC.add("globalMCrec/GapSideTrue", "Recalculated gap side; gap index; Collisions", {HistType::kTH1F, {{5, -1, 4}}});
    registryMC.add("globalMCrec/hVertexXY", "Vertex position in x and y direction; #it{V}_{x} (cm); #it{V}_{y} (cm); Collisions", {HistType::kTH2F, {{50, -0.05, 0.05}, {50, -0.02, 0.02}}});
    registryMC.add("globalMCrec/hVertexZ", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    registryMC.add("globalMCrec/hNTracks", ";N_{tracks};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registryMC.add("globalMCrec/hNTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registryMC.add("globalMCrec/hNGhostTracks", ";N_{tracks};events", {HistType::kTH1D, {{10, 0., 10.}}});
    registryMC.add("globalMCrec/hNGhostTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{10, 0., 10.}}});
    registryMC.add("globalMCrec/hQtot", ";Q_{tot};events", {HistType::kTH1F, {{10, -5., 5.}}});
    registryMC.add("globalMCrec/hTrackToMCMatch", ";Match};events", {HistType::kTH1F, {{2, 0., 2.}}});
    registryMC.add("globalMCrec/hZNACenergy", "ZNA vs ZNC energy; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("globalMCrec/hZNACtime", "ZNA vs ZNC time; #it{time}_{ZNA} (ns); #it{time}_{ZNC} (ns); Collisions", {HistType::kTH2F, {{100, -10., 10.}, {100, -10., 10.}}});

    registry1MC.add("globalMCrec/hRecFlag", ";Reconstruction Flag;events", {HistType::kTH1F, {{10, 0., 10.}}});
    registry1MC.add("globalMCrec/hOccupancyInTime", ";Occupancy;events", {HistType::kTH1F, {{100, 0., 10000.}}});

    // registryMC.add("globalMCrec/hPtSpectrumElGen0", "Gen.;p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});//effiEl = 3 // hpTelec
    // registryMC.add("globalMCrec/hPtSpectrumElGen1", "Gen.;p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});//effiEl = 4, but still on MC sample
    registryMC.add("globalMCrec/hPtSpectrumElRec0", "Rec0;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 4, reconstruction
    registryMC.add("globalMCrec/hPtSpectrumElRec1", "Rec1;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 7, SGProducer
    registryMC.add("globalMCrec/hPtSpectrumElRec2", "Rec2;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 9, DoubleGap
    registryMC.add("globalMCrec/hPtSpectrumElRec3", "Rec3;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 11, PVtracks=4
    registryMC.add("globalMCrec/hPtSpectrumElRec4", "Rec4;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 13, Zvertex
    registryMC.add("globalMCrec/hPtSpectrumElRec5", "Rec5;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 15, eta acceptance
    registryMC.add("globalMCrec/hPtSpectrumElRec6", "Rec6;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 17, pT threshold
    registryMC.add("globalMCrec/hPtSpectrumElRec7", "Rec7;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 19, Qtot
    registryMC.add("globalMCrec/hPtSpectrumElRec8", "Rec8;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 21, NTOF>0
    registryMC.add("globalMCrec/hPtSpectrumElRec9", "Rec9;#it{p}_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}}); // effiEl = 23, FIT empty

    // global1 when we require SGProducer + double gap + nPVtracks=4
    registryMC.add("global1MCrec/hVertexZ", "Vertex position in z direction; #it{V}_{z} (cm); Collisions", {HistType::kTH1F, {{100, -25., 25.}}});
    registryMC.add("global1MCrec/hNTracks", ";N_{tracks};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registryMC.add("global1MCrec/hNTracksPV", ";N_{tracks,PV};events", {HistType::kTH1D, {{20, 0., 20.}}});
    registryMC.add("global1MCrec/hTrackPtPV", ";p_T^{trk}; Entries", {HistType::kTH1F, {axispt}});
    registryMC.add("global1MCrec/hTrackEtaPhiPV", ";Eta;Phi;", {HistType::kTH2D, {axiseta, {128, -0.05, 6.35}}});
    registryMC.add("global1MCrec/hTrackPVTotCharge", "Q_{Tot};Q_{Tot}; Entries", {HistType::kTH1F, {{11, -5, 6}}});

    registryMC.add("global1MCrec/hpTGenRecTracksPV", ";p_{T}^{Rec. tracks,PV} (GeV/c);p_{T}^{Gen} (GeV/c);events", {HistType::kTH2D, {{100, 0., 4.}, {100, 0., 4.}}});
    registryMC.add("global1MCrec/hDeltapTGenRecVsRecpTTracksPV", ";#Delta p_{T}^{Rec.-Gen. tracks,PV} (GeV/c);p_{T}^{Rec. tracks,PV} (GeV/c);events", {HistType::kTH2D, {{100, -4., 4.}, {100, 0., 4.}}});
    registryMC.add("global1MCrec/hEtaGenRecTracksPV", ";#eta^{Rec. tracks,PV} (GeV/c);#eta^{Gen} (GeV/c);events", {HistType::kTH2D, {{100, -2., 2.}, {100, -2., 2.}}});
    registryMC.add("global1MCrec/hDeltaEtaGenRecVsRecpTTracksPV", ";#Delta #eta^{Rec.-Gen. tracks,PV} (GeV/c);p_{T}^{Rec. tracks,PV} (GeV/c);events", {HistType::kTH2D, {{100, -0.25, 0.25}, {100, 0., 4.}}});
    registryMC.add("global1MCrec/hPhiGenRecTracksPV", ";#phi^{Rec. tracks,PV} (GeV/c);#phi^{Gen} (GeV/c);events", {HistType::kTH2D, {{100, 0., 6.4}, {100, 0., 6.4}}});
    registryMC.add("global1MCrec/hDeltaPhiGenRecVsRecpTTracksPV", ";#Delta #phi^{Rec.-Gen. tracks,PV} (GeV/c);p_{T}^{Rec. tracks,PV} (GeV/c);events", {HistType::kTH2D, {{100, -0.5, 0.5}, {100, 0., 4.}}});

    // pid El in MC true
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut1", "All hip;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // pid separately for each cut (what we reject)
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut2", "rejected, IM hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut3", "rejected, DP hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut4", "rejected, El hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut5", "rejected, vPi hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut6", "rejected, Pt hip; #it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut7", "rejected, vVc hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut8", "rejected, pTot hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut9", "rejected, vPr hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut10", "rejected, vKa hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut11", "rejected, nCR hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut12", "rejected, s3pi hip;#it{p}_{trk} (GeV/#it{c}); d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // pid sequentialy
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut20", "El hip;    #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut33", "eTOF+20 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut21", "vPi+33 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut24", "vPr+21 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut28", "CR+25 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut22", "vVc+28 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    //   registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut23", "Pt+22 hip; #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut29", "s3pi+22 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut26", "IM+29 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut34", "piTOF+26 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut30", "ptTot+34 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut27", "DP+30 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    //    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut31", "FIT+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    //    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut32", "TOF+31 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCEltrue/hpvsdedxElHipCut35", "ZDC+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // pid Pi in MC true
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});dE/dx_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut20", "El hip;    #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut33", "eTOF+20 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut21", "vPi+33 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut24", "vPr+21 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut28", "CR+25 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut22", "vVc+28 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut23", "Pt+22 hip; #it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut29", "s3pi+22 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut26", "IM+29 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut34", "piTOF+26 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut30", "ptTot+34 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut27", "DP+30 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut31", "FIT+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    // registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut32", "TOF+31 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});
    registryMC.add("pidTPCMCPitrue/hpvsdedxElHipCut35", "ZDC+27 hip;#it{p}_{trk} (GeV/#it{c});d#it{E}/d#it{x}_{trk}", {HistType::kTH2F, {axisp, dedxAxis}});

    // El PID in TOF MC true
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut0", "In hip ;#it{p}_{trk}(GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut20", "El hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut33", "eTOF+20 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut21", "vPi+33 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut24", "vPr+21 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut25", "vKa+24 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut28", "CR+25 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut22", "vVc+28 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut29", "s3pi+22 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut26", "IM+29 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut34", "piTOF+26 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut30", "ptTot+34 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut27", "DP+30 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});
    registry1MC.add("pidTOFMCEltrue/hpvsNsigmaElHipCut35", "ZDC+27 hip;#it{p}_{trk} (GeV/#it{c});N#sigma El^{TOF}_{trk}", {HistType::kTH2F, {axisp, {100, -5., 5.}}});

    // cut0
    registryMC.add("controlMCtrue/cut0/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut0/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut0/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut0/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut0/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut0/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut0/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut0/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut0/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut0/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut0/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut0/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut0/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    registryMC.add("controlMCtrue/cut0/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut0/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registry1MC.add("controlMCtrue/cut0/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut0/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut0/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut0/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut0/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut0/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut0/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut0/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCcomb/cut0/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registry1MC.add("controlMCcomb/cut0/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut20 MC
    registryMC.add("controlMCtrue/cut20/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut20/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut20/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut20/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut20/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut20/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut20/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut20/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut20/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut20/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut20/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut20/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut20/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut20/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut20/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut20/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut20/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut20/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut20/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut20/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut20/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut20/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut20/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut20/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut20/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut21 MC
    registryMC.add("controlMCtrue/cut21/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut21/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut21/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut21/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut21/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut21/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut21/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut21/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut21/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut21/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut21/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut21/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut21/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut21/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut21/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut21/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut21/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut21/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut21/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut21/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut21/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut21/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut21/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut21/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut21/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut24 MC
    registryMC.add("controlMCtrue/cut24/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut24/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut24/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut24/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut24/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut24/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut24/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut24/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut24/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut24/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut24/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut24/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut24/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut24/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut24/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut24/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut24/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut24/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut24/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut24/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut24/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut24/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut24/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut24/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut24/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut25 MC
    registryMC.add("controlMCtrue/cut25/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut25/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut25/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut25/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut25/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut25/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut25/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut25/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut25/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut25/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut25/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut25/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut25/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut25/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut25/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut25/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut25/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut25/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut25/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut25/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut25/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut25/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut25/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut25/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut25/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut28 MC
    registryMC.add("controlMCtrue/cut28/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut28/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut28/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut28/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut28/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut28/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut28/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut28/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut28/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut28/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut28/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut28/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut28/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut28/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut28/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut28/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut28/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut28/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut28/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut28/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut28/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut28/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut28/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut28/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut28/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut22 MC
    registryMC.add("controlMCtrue/cut22/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut22/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut22/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut22/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut22/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut22/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut22/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut22/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut22/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut22/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut22/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut22/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut22/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut22/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut22/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut22/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut22/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut22/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut22/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut22/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut22/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut22/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut22/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut22/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut22/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //    //cut23 MC
    //    registryMC.add("controlMCtrue/cut23/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registryMC.add("controlMCtrue/cut23/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //    registryMC.add("controlMCtrue/cut23/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //    registryMC.add("controlMCtrue/cut23/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //    registryMC.add("controlMCtrue/cut23/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registryMC.add("controlMCtrue/cut23/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registryMC.add("controlMCtrue/cut23/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registryMC.add("controlMCtrue/cut23/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registryMC.add("controlMCtrue/cut23/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registryMC.add("controlMCtrue/cut23/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registryMC.add("controlMCtrue/cut23/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registryMC.add("controlMCtrue/cut23/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //    registryMC.add("controlMCtrue/cut23/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry1MC.add("controlMCtrue/cut23/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    // registryMC.add("controlMCtrue/cut23/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //    // registryMC.add("controlMCtrue/cut23/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    //
    //    registryMC.add("controlMCcomb/cut23/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registryMC.add("controlMCcomb/cut23/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registryMC.add("controlMCcomb/cut23/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registryMC.add("controlMCcomb/cut23/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registryMC.add("controlMCcomb/cut23/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registryMC.add("controlMCcomb/cut23/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registryMC.add("controlMCcomb/cut23/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registryMC.add("controlMCcomb/cut23/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry1MC.add("controlMCcomb/cut23/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut29 MC
    registryMC.add("controlMCtrue/cut29/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut29/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut29/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut29/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut29/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut29/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut29/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut29/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut29/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut29/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut29/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut29/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut29/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut29/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut29/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut29/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut29/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut29/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut29/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut29/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut29/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut29/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut29/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut29/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut29/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut26 MC
    registryMC.add("controlMCtrue/cut26/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut26/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut26/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut26/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut26/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut26/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut26/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut26/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut26/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut26/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut26/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut26/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut26/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut26/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut26/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut26/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut26/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut26/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut26/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut26/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut26/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut26/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut26/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut26/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut26/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut30 MC
    registryMC.add("controlMCtrue/cut30/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut30/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut30/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut30/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut30/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut30/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut30/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut30/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut30/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut30/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut30/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut30/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut30/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut30/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut30/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut30/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut30/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut30/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut30/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut30/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut30/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut30/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut30/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut30/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut30/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut27 MC
    registryMC.add("controlMCtrue/cut27/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut27/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut27/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut27/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut27/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut27/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut27/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut27/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut27/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut27/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut27/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut27/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut27/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut27/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut27/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut27/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut27/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut27/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut27/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut27/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut27/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut27/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut27/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut27/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut27/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //   //cut31 MC
    //   registryMC.add("controlMCtrue/cut31/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //   registryMC.add("controlMCtrue/cut31/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //   registryMC.add("controlMCtrue/cut31/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //   registryMC.add("controlMCtrue/cut31/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //   registryMC.add("controlMCtrue/cut31/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //   registryMC.add("controlMCtrue/cut31/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //   registryMC.add("controlMCtrue/cut31/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //   registryMC.add("controlMCtrue/cut31/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //   registryMC.add("controlMCtrue/cut31/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //   registryMC.add("controlMCtrue/cut31/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //   registryMC.add("controlMCtrue/cut31/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //   registryMC.add("controlMCtrue/cut31/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //   registryMC.add("controlMCtrue/cut31/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //   registry1MC.add("controlMCtrue/cut31/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //   // registryMC.add("controlMCtrue/cut31/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //   // registryMC.add("controlMCtrue/cut31/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //   //
    //   registryMC.add("controlMCcomb/cut31/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //   registryMC.add("controlMCcomb/cut31/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //   registryMC.add("controlMCcomb/cut31/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //   registryMC.add("controlMCcomb/cut31/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //   registryMC.add("controlMCcomb/cut31/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //   registryMC.add("controlMCcomb/cut31/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //   registryMC.add("controlMCcomb/cut31/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //   registryMC.add("controlMCcomb/cut31/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //   registry1MC.add("controlMCcomb/cut31/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    //    //cut32 MC
    //    registryMC.add("controlMCtrue/cut32/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registryMC.add("controlMCtrue/cut32/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    //    registryMC.add("controlMCtrue/cut32/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    //    registryMC.add("controlMCtrue/cut32/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    //    registryMC.add("controlMCtrue/cut32/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registryMC.add("controlMCtrue/cut32/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registryMC.add("controlMCtrue/cut32/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registryMC.add("controlMCtrue/cut32/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registryMC.add("controlMCtrue/cut32/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registryMC.add("controlMCtrue/cut32/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registryMC.add("controlMCtrue/cut32/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registryMC.add("controlMCtrue/cut32/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    //    registryMC.add("controlMCtrue/cut32/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry1MC.add("controlMCtrue/cut32/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    // registryMC.add("controlMCtrue/cut32/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    //    // registryMC.add("controlMCtrue/cut32/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    //
    //    registryMC.add("controlMCcomb/cut32/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    //    registryMC.add("controlMCcomb/cut32/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    //    registryMC.add("controlMCcomb/cut32/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    //    registryMC.add("controlMCcomb/cut32/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    //    registryMC.add("controlMCcomb/cut32/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    //    registryMC.add("controlMCcomb/cut32/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    //    registryMC.add("controlMCcomb/cut32/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    //    registryMC.add("controlMCcomb/cut32/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //    registry1MC.add("controlMCcomb/cut32/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut33 MC
    registryMC.add("controlMCtrue/cut33/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut33/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut33/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut33/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut33/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut33/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut33/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut33/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut33/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut33/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut33/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut33/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut33/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut33/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut33/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut33/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut33/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut33/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut33/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut33/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut33/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut33/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut33/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut33/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut33/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut34 MC
    //     registryMC.add("controlMCtrue/cut34/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut34/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {axisInvMass4trk}});
    registryMC.add("controlMCtrue/cut34/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut34/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut34/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut34/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut34/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut34/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut34/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut34/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut34/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut34/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut34/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut34/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut34/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut34/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut34/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut34/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut34/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut34/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut34/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut34/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut34/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut34/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut34/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut34/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // cut35 MC
    //     registryMC.add("controlMCtrue/cut34/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registryMC.add("controlMCtrue/cut35/h4piMass", "4#pi mass;M_{inv}^{4#pi} (GeV/c^{2});entries", {HistType::kTH1F, {axisInvMass4trk}});
    registryMC.add("controlMCtrue/cut35/h4trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut35/h4trkMassVsPt", "4-track mass vs Pt;M_{inv}^{4track} (GeV/c^{2});p_{T}^{4track} (GeV/c);entries", {HistType::kTH2F, {{100, 1, 5.}, axispt}});
    registryMC.add("controlMCtrue/cut35/hZNACenergy", "ZNA vs ZNC energy, cut0; #it{E}_{ZNA} (GeV); #it{E}_{ZNC} (GeV); Collisions", {HistType::kTH2F, {axisZDC, axisZDC}});
    registryMC.add("controlMCtrue/cut35/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCtrue/cut35/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCtrue/cut35/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCtrue/cut35/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCtrue/cut35/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCtrue/cut35/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCtrue/cut35/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCtrue/cut35/h3piMassVsPt", "3#pi mass vs Pt, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});p_{T}^{3#pi} (GeV/c);entries", {HistType::kTH2F, {minvAxis, axispt}});
    registryMC.add("controlMCtrue/cut35/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCtrue/cut35/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    // registryMC.add("controlMCtrue/cut35/hNtofTrk", ";N_{TOF trk}; Entries", {HistType::kTH1F, {{7, 0., 7.}}});
    // registryMC.add("controlMCtrue/cut35/h3pi1eMass", "3#pi+e mass;M_{inv}^{3#pi+e} (GeV/c^{2});entries", {HistType::kTH1F, {{100, 0., 10.}}});
    //
    registryMC.add("controlMCcomb/cut35/h3piMass", "3#pi mass, up to 4 entries per event ;M_{inv}^{3#pi} (GeV/c^{2});entries", {HistType::kTH1F, {minvAxis}});
    registryMC.add("controlMCcomb/cut35/h3trkPtTot", ";p_{T} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut35/hDeltaPhi13topo", "#Delta #varphi 1+3 topo, 4 entries/event;#Delta#varphi^{1+3};entries", {HistType::kTH1F, {phiAxis}});
    registryMC.add("controlMCcomb/cut35/h13AssymPt1ProngAver", ";Delta Pt/Sum Pt (1Prong,Aver Pt)", {HistType::kTH1F, {{100, -1., 1.}}});
    registryMC.add("controlMCcomb/cut35/h13Vector", ";A_{V};entries", {HistType::kTH1F, {vectorAxis}});
    registryMC.add("controlMCcomb/cut35/h13Scalar", ";A_{S};entries", {HistType::kTH1F, {scalarAxis}});
    registryMC.add("controlMCcomb/cut35/hTPCnCrossedRows", "N crossed rows ;N_{TPC,crossed rows};entries", {HistType::kTH1F, {{160, 0, 160.}}});
    registryMC.add("controlMCcomb/cut35/hsigma3Pi", "#sqrt{#sigma_{1}^{2 }+#sigma_{2}^{2}+#sigma_{3}^{2}};#sigma^{3#pi};entries", {HistType::kTH1F, {{100, 0., 10.}}});
    registry1MC.add("controlMCcomb/cut35/hTofChi2El", ";TOF #chi^{2};entries", {HistType::kTH1F, {{100, 0., 10.}}});

    // ptSpectrum of electron for MC true and combinatorics
    registryMC.add("controlMCtrue/cut0/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut20/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut21/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut22/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut23/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut24/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut25/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut26/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut27/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut28/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut29/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut30/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut31/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    // registryMC.add("controlMCtrue/cut32/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut33/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut34/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCtrue/cut35/hPtSpectrumEl", ";p_{T}^{e} (GeV/c);entries", {HistType::kTH1F, {axispt}});

    registryMC.add("controlMCcomb/cut0/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {ptAxis}});
    registryMC.add("controlMCcomb/cut20/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut21/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut22/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut23/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut24/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut25/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut26/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut27/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut28/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut29/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut30/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut31/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut32/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut33/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut34/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    registryMC.add("controlMCcomb/cut35/hPtSpectrumEl", ";p_{T}^{comb} (GeV/c);entries", {HistType::kTH1F, {axispt}});
    // zrobic hpTspectrumEl dla cut 0,20-35: registry.get<TH1>(HIST("global/hFinalPtSpectrumEl"))->Fill(tmpPt[i]);

    // tau
    registryMC.add("tauMC/hMCeta", ";#eta^{#tau};N^{#tau} ", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("tauMC/hMCy", ";y^{#tau};N^{#tau}", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("tauMC/hMCphi", ";#phi^{#tau};N^{#tau}", {HistType::kTH1F, {{100, 0., 6.4}}});
    registryMC.add("tauMC/hMCpt", ";#it{p}_{T}^{#tau};N^{#tau}", {HistType::kTH1F, {{100, 0., 10.}}});

    registryMC.add("tauMC/hMCdeltaeta", ";#Delta#eta^{#tau};events ", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("tauMC/hMCdeltaphi", ";#Delta#phi^{#tau}(deg.);events", {HistType::kTH1F, {{100, 131., 181}}});

    // electron
    registryMC.add("electronMC/hMCeta", ";#eta^{e};N^{e}", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("electronMC/hMCy", ";y^{e};N^{e}", {HistType::kTH1F, {{100, -5., 5.}}});
    registryMC.add("electronMC/hMCphi", ";#phi^{e};N^{e}", {HistType::kTH1F, {{100, 0., 6.4}}});
    registryMC.add("electronMC/hMCpt", ";#it{p}_{T}^{e};N^{e}", {HistType::kTH1F, {{400, 0., 10.}}});

    // efficiency mu
    registryMC.add("efficiencyMCMu/effiMu", ";Efficiency #mu3#pi;events", {HistType::kTH1F, {{20, 0., 20.}}});
    registryMC.add("efficiencyMCMu/hpTmuon", ";p_{T}^{#mu, gen} (GeV/c);events", {HistType::kTH1F, {{200, 0., 5.}}});

    // efficiency pi
    registryMC.add("efficiencyMCPi/effiPi", ";Efficiency #pi3#pi;events", {HistType::kTH1F, {{20, 0., 20.}}});
    registryMC.add("efficiencyMCPi/hpTpi", ";p_{T}^{#pi, gen} (GeV/c);events", {HistType::kTH1F, {{200, 0., 5.}}});

    // efficiency el
    registryMC.add("efficiencyMCEl/effiEl", ";Efficiency e3#pi;events", {HistType::kTH1F, {{70, 0., 70.}}});
    registryMC.add("efficiencyMCEl/hpTelec", ";p_{T}^{e, gen} (GeV/c);events", {HistType::kTH1F, {axispt}});

    // efficiency el
    registryMC.add("efficiencyMCEl/hMCdeltaAlphaEpi", ";#Delta#alpha^{e-#pi(3x), gen} (rad);events", {HistType::kTH1F, {{100, -0.1, 3.2}}});
    registryMC.add("efficiencyMCEl/hMCdeltaAlphaPiPi", ";#Delta#alpha^{#pi-#pi(3x), gen} (rad);events", {HistType::kTH1F, {{100, -0.1, 3.2}}});
    registryMC.add("efficiencyMCEl/hMCdeltaPhiEpi", ";#Delta#phi^{e-#pi(3x), gen} (rad);events", {HistType::kTH1F, {{100, -0.1, 3.2}}});
    registryMC.add("efficiencyMCEl/hMCdeltaPhiPipi", ";#Delta#phi^{#pi-#pi(3x), gen} (rad);events", {HistType::kTH1F, {{100, -0.1, 3.2}}});
    registryMC.add("efficiencyMCEl/hMCvirtCal", ";virt Cal. #Delta #alpha^{#pi-#pi(3x), gen} (rad);events", {HistType::kTH1F, {{100, 0., 2.}}});
    registryMC.add("efficiencyMCEl/hMCScalar", ";A_{S}^{#pi-#pi(3x), gen};events", {HistType::kTH1F, {{100, 0., 1.}}});
    registryMC.add("efficiencyMCEl/hMCVector", ";A_{V}^{#pi-#pi(3x), gen};events", {HistType::kTH1F, {{100, 0., 2.}}});

    // missing eta vs phi
    registryMC.add("efficiencyMCEl/hMCptEl", ";p_{T}^{e, true} (GeV/c) ;events", {HistType::kTH1F, {{200, 0., 5.}}});
    // registryMC.add("efficiencyMCEl/hMCpt4trk",";p_{T}^{4 MC part.} (GeV/c) ;events",{HistType::kTH1F,{{100, 0., 5.}}});
    registryMC.add("efficiencyMCEl/hMCpt4trk", ";p_{T}^{4 MC part.} (GeV/c) ;events", {HistType::kTH1F, {axispt}});
    // registryMC.add("efficiencyMCEl/hMCinvmass4pi",";M_{inv}^{4#pi true} (GeV/c^{2}) ;events",{HistType::kTH1F,{{100, 0.4, 5.4}}});
    registryMC.add("efficiencyMCEl/hMCinvmass4pi", ";M_{inv}^{4#pi true} (GeV/c^{2}) ;events", {HistType::kTH1F, {axisInvMass4trk}});
    registryMC.add("efficiencyMCEl/hMCinvmass3pi", ";M_{inv}^{3#pi true} (GeV/c^{2}) ;events", {HistType::kTH1F, {{100, 0.4, 2.4}}});
    registryMC.add("efficiencyMCEl/hMCdeltaphi13", ";#Delta#phi^{1-3 true} ;events", {HistType::kTH1F, {{100, 0., 3.2}}});
  }

  float eta(float px, float py, float pz)
  // Just a simple function to return pseudorapidity
  {
    float arg = -2.; // outside valid range for std::atanh
    float mom = std::sqrt(px * px + py * py + pz * pz);
    if (mom != 0)
      arg = pz / mom;
    if (-1. < arg && arg < 1.)
      return std::atanh(arg); // definition of eta
    return -999.;
  }

  float phi(float px, float py)
  // Just a simple function to return azimuthal angle from 0 to 2pi
  {
    if (px != 0)
      return (std::atan2(py, px) + o2::constants::math::PI);
    return -999.;
  }

  float pt(float px, float py)
  // Just a simple function to return pt
  {
    return std::sqrt(px * px + py * py);
  }

  float rapidity(float energy, float pz)
  // Just a simple function to return track rapidity
  {
    return 0.5 * std::log((energy + pz) / (energy - pz));
  }

  // helper function to calculate delta alpha
  float deltaAlpha(auto particle1, auto particle2)
  {

    TVector3 vtmp(particle1.px(), particle1.py(), particle1.pz());
    TVector3 v1(particle2.px(), particle2.py(), particle2.pz());
    auto angle = v1.Angle(vtmp);

    return angle;
  }

  float invariantMass(float E, float px, float py, float pz)
  // Just a simple function to return invariant mass
  {
    return std::sqrt(E * E - px * px - py * py - pz * pz);
  }

  float calculateDeltaPhi(TLorentzVector p, TLorentzVector p1)
  {
    //    float delta = p.Phi();
    float delta = RecoDecay::constrainAngle(p.Phi());
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    // if (p1.Phi() < 0)
    //   delta -= (p1.Phi() + o2::constants::math::TwoPI);
    // else
    delta -= RecoDecay::constrainAngle(p1.Phi());
    delta = RecoDecay::constrainAngle(delta);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    if (delta > o2::constants::math::PI)
      delta = o2::constants::math::TwoPI - delta;
    return delta;
  }

  float calculateDeltaPhi(float p, float p1)
  {
    float delta = RecoDecay::constrainAngle(p);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    // if (p1 < 0)
    //   delta -= (p1 + o2::constants::math::TwoPI);
    // else
    delta -= RecoDecay::constrainAngle(p1);
    delta = RecoDecay::constrainAngle(delta);
    // if (delta < 0)
    //   delta += o2::constants::math::TwoPI;
    if (delta > o2::constants::math::PI)
      delta = o2::constants::math::TwoPI - delta;
    return delta;
  }

  //  // helper function to calculate scalar asymmetry
  //  float scalarAsym(auto particle1, auto particle2){
  //  auto delta = particle1.pt() - particle2.pt();
  //  return TMath::Abs(delta)/(particle1.pt() + particle2.pt());
  //  }
  //
  //  // helper function to calculate vector asymmetry
  //  float vectorAsym(auto particle1, auto particle2){
  //   auto delta = TMath::Sqrt((particle1.px() - particle2.px()) * (particle1.px() - particle2.px()) +
  //   (particle1.py() - particle2.py()) * (particle1.py() - particle2.py()));
  //   auto sum =  TMath::Sqrt((particle1.px() + particle2.px()) * (particle1.px() + particle2.px()) +
  //   (particle1.py() + particle2.py()) * (particle1.py() + particle2.py()));
  //   return sum/delta;
  //  }

  // helper function to calculate scalar asymmetry
  float scalarAsymMC(auto particle1, auto particle2)
  {
    auto pt1 = pt(particle1.px(), particle1.py());
    auto pt2 = pt(particle2.px(), particle2.py());
    auto delta = pt1 - pt2;
    return std::abs(delta) / (pt1 + pt2);
  }

  // helper function to calculate vector asymmetry
  float vectorAsym(auto particle1, auto particle2)
  {
    auto delta = std::sqrt((particle1.px() - particle2.px()) * (particle1.px() - particle2.px()) +
                           (particle1.py() - particle2.py()) * (particle1.py() - particle2.py()));
    auto sum = std::sqrt((particle1.px() + particle2.px()) * (particle1.px() + particle2.px()) +
                         (particle1.py() + particle2.py()) * (particle1.py() + particle2.py()));
    return sum / delta;
  }

  // fill control histograms per track
  template <int mode, typename T>
  //   void fillControlHistos(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float pi3etasum, float nCRtpc)
  void fillControlHistos(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float nCRtpc, float ptelec, float tofchi2)
  {
    static constexpr std::string_view kHistoname[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                      "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                      "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                                                      "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"};
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h3piMassComb"))->Fill(pi3invMass);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h3trkPtTot"))->Fill(pi3pt);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/hDeltaPhi13topo"))->Fill(pi3deltaPhi);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h13AssymPt1ProngAver"))->Fill(pi3assymav);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h13Vector"))->Fill(pi3vector);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h13Scalar"))->Fill(pi3scalar);
    // registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/h13EtaSum"))->Fill(pi3etasum);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/hTPCnCrossedRows"))->Fill(nCRtpc);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/hPtSpectrumEl"))->Fill(ptelec);
    registry.get<TH1>(HIST("control/cut") + HIST(kHistoname[mode]) + HIST("/hTofChi2El"))->Fill(tofchi2);
  }

  // fill control histograms per track in MC true
  template <int mode, typename T>
  void fillControlHistosMCtrue(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float nCRtpc, float ptelec, float tofchi2)
  {
    static constexpr std::string_view kHistoname[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                      "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                      "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                                                      "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"};
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h3piMass"))->Fill(pi3invMass);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h3trkPtTot"))->Fill(pi3pt);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/hDeltaPhi13topo"))->Fill(pi3deltaPhi);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h13AssymPt1ProngAver"))->Fill(pi3assymav);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h13Vector"))->Fill(pi3vector);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h13Scalar"))->Fill(pi3scalar);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/hTPCnCrossedRows"))->Fill(nCRtpc);
    registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/hPtSpectrumEl"))->Fill(ptelec);
    // registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h13EtaSum"))->Fill(pi3etasum);
    registry1MC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/hTofChi2El"))->Fill(tofchi2);
  }

  // fill control histograms per track in MC true
  template <int mode, typename T>
  void fillControlHistosMCcomb(T pi3invMass, float pi3pt, float pi3deltaPhi, float pi3assymav, float pi3vector, float pi3scalar, float nCRtpc, float ptelec, float tofchi2)
  {
    static constexpr std::string_view kHistoname[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
                                                      "10", "11", "12", "13", "14", "15", "16", "17", "18", "19",
                                                      "20", "21", "22", "23", "24", "25", "26", "27", "28", "29",
                                                      "30", "31", "32", "33", "34", "35", "36", "37", "38", "39"};
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/h3piMass"))->Fill(pi3invMass);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/h3trkPtTot"))->Fill(pi3pt);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/hDeltaPhi13topo"))->Fill(pi3deltaPhi);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/h13AssymPt1ProngAver"))->Fill(pi3assymav);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/h13Vector"))->Fill(pi3vector);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/h13Scalar"))->Fill(pi3scalar);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/hTPCnCrossedRows"))->Fill(nCRtpc);
    registryMC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/hPtSpectrumEl"))->Fill(ptelec);
    // registryMC.get<TH1>(HIST("controlMCtrue/cut") + HIST(kHistoname[mode]) + HIST("/h13EtaSum"))->Fill(pi3etasum);
    registry1MC.get<TH1>(HIST("controlMCcomb/cut") + HIST(kHistoname[mode]) + HIST("/hTofChi2El"))->Fill(tofchi2);
  }

  template <typename T>
  int trackCheck(T track)
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
  Filter pVContributorFilter = aod::udtrack::isPVContributor == true;
  using PVTracks = soa::Filtered<UDTracksFull>;

  //  using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  using LabeledTracks = soa::Join<aod::UDTracks, aod::UDMcTrackLabels, aod::UDTracksExtra, aod::UDTracksPID, aod::UDTracksFlags, aod::UDTracksDCA>;
  Preslice<aod::UDTracks> perCollision = aod::udtrack::udCollisionId;
  // PVContributors in MC handling
  Filter pVContributorFilterMC = aod::udtrack::isPVContributor == true;
  using PVTracksMC = soa::Filtered<LabeledTracks>;

  //  void processDG(UDCollisionFull const& dgcand, UDTracksFull const& dgtracks)
  //  {
  //    int gapSide = 2;
  //    mainTask(gapSide,dgtracks);
  //  }
  //  PROCESS_SWITCH(TauTau13topo, processDG, "Process DG data", DGactive);

  // void processSG(UDCollisionFull2 const& dgcand, UDTracksFull const& dgtracks)
  void processDataSG(UDCollisionFull2 const& dgcand, UDTracksFull const& dgtracks, PVTracks const& PVContributors)
  {
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(0., 1.);
    registry.get<TH1>(HIST("global/RunNumber"))->Fill(dgcand.runNumber());
    registry.get<TH1>(HIST("global/hRecFlag"))->Fill(dgcand.flags()); // reconstruction with upc settings flag
    // registry.get<TH1>(HIST("global/hOccupancyInTime"))->Fill(dgcand.occupancyInTime());

    int gapSide = dgcand.gapSide();
    int truegapSide = sgSelector.trueGap(dgcand, cutFV0, cutFT0A, cutFT0C, cutZDC);
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
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(1., 1.);
    if (gapSide != mGapSide)
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
    float energyZNA = dgcand.energyCommonZNA();
    float energyZNC = dgcand.energyCommonZNC();
    // if (energyZNA < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(energyZNA);
    // if (energyZNC < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(energyZNC);
    if (energyZNA < 0)
      energyZNA = -1.;
    if (energyZNC < 0)
      energyZNC = -1.;
    registry.get<TH2>(HIST("global/hZNACenergy"))->Fill(energyZNA, energyZNC);
    registry.get<TH2>(HIST("global/hZNACtime"))->Fill(dgcand.timeZNA(), dgcand.timeZNC());

    uint32_t clusterSizes;
    //    UChar_t clustermap1;
    //    bool isInnerITS = false;
    int nTofTrk = 0;
    int nEtaIn15 = 0;
    int nITSbits = -1;
    int npT100 = 0;
    TLorentzVector p;
    // auto const pionMass = MassPiPlus;
    // auto const electronMass = MassElectron;
    bool flagGlobalCheck = true;
    bool isGlobalTrack = true;
    int qtot = 0;
    // loop over PV contributors
    for (const auto& trk : PVContributors) {
      qtot += trk.sign();
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
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
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(2., 1.);

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

    registry.get<TH1>(HIST("global/hEventEff"))->Fill(3., 1.);
    // registry.get<TH1>(HIST("global1/hTrackPVTotCharge"))->Fill(dgcand.netCharge());
    registry.get<TH1>(HIST("global1/hTrackPVTotCharge"))->Fill(qtot);
    registry.get<TH1>(HIST("global1/hVertexZ"))->Fill(dgcand.posZ());
    registry.get<TH1>(HIST("global1/hNTracks"))->Fill(dgtracks.size());
    registry.get<TH1>(HIST("global1/hNTracksPV"))->Fill(PVContributors.size());
    for (const auto& trk : PVContributors) {
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
      registry.get<TH1>(HIST("global1/hTrackPtPV"))->Fill(p.Pt());
      registry.get<TH2>(HIST("global1/hTrackEtaPhiPV"))->Fill(p.Eta(), p.Phi());
    }

    // if vz < 10
    if (std::abs(dgcand.posZ()) >= zvertexcut) { // default = 10
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: VertexZ is %f", dgcand.posZ());
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(4., 1.);

    // if eta tracks <1.5
    if (nEtaIn15 != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: Ntrk inside |eta|<0.9 is %d", nEtaIn15);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(5., 1.);

    // if pt tracks >0.100 GeV/c
    if (npT100 != 4) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: number of tracks with pT>0.1GeV/c is %d", npT100);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(6., 1.);

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
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(7., 1.);

    //    // skip events with out-of-range rgtrwTOF (fraction-of-good-tracks-with-TOF-hit)
    //    auto rtrwTOF = udhelpers::rPVtrwTOF<false>(dgtracks, PVContributors.size());
    //    if (rtrwTOF < 0.25) {
    //      if (verbose) {
    //        LOGF(debug, "<tautau13topo> Candidate rejected: rtrwTOF is %f", rtrwTOF);
    //      }
    //      return;
    //    }

    // n TOF tracks cut
    if (nTofTrk < nTofTrkMinCut) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: nTOFtracks is %d", nTofTrk);
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(8., 1.);

    //
    // FIT informaton
    //
    auto bitMin = 16 - mFITvetoWindow; // default is +- 1 bc (1 bit)
    auto bitMax = 16 + mFITvetoWindow;
    bool flagFITveto = false;
    // check FIT information
    for (auto bit = bitMin; bit <= bitMax; bit++) {
      if (TESTBIT(dgcand.bbFT0Apf(), bit))
        flagFITveto = true;
      if (TESTBIT(dgcand.bbFT0Cpf(), bit))
        flagFITveto = true;
      if (useFV0ForVeto && TESTBIT(dgcand.bbFV0Apf(), bit))
        flagFITveto = true;
      if (useFDDAForVeto && TESTBIT(dgcand.bbFDDApf(), bit))
        flagFITveto = true;
      if (useFDDCForVeto && TESTBIT(dgcand.bbFDDCpf(), bit))
        flagFITveto = true;
    } // end of loop over FIT bits
    // FIT histos
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

    // FIT empty
    if (mFITvetoFlag && flagFITveto) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: FIT is not empty");
      }
      return;
    }
    registry.get<TH1>(HIST("global/hEventEff"))->Fill(9., 1.);

    //
    // here PID from TPC starts to be
    //
    // temporary control variables per event with combinatorics
    float tmpMomentum[4];
    float tmpPt[4];
    float tmpDedx[4];
    float tmpTofNsigmaEl[4];
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
    // float dcaZ[4];
    // float dcaXY[4];
    // float chi2TPC[4];
    // float chi2ITS[4];
    float chi2TOF[4];
    // float nclTPCfind[4];
    float nclTPCcrossedRows[4];
    bool trkHasTof[4];
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
      p.SetXYZM(trk.px(), trk.py(), trk.pz(), MassElectron);
      for (const auto& trk1 : PVContributors) {
        if (trk.index() >= trk1.index())
          continue;
        p1.SetXYZM(trk1.px(), trk1.py(), trk1.pz(), MassElectron);
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
      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
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
      tmpTrkCheck = trackCheck(trk); // check detectors associated to track
      registry.get<TH1>(HIST("global/hTrkCheck"))->Fill(tmpTrkCheck);

      // inv mass of 3pi + 1e
      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
      p2.SetXYZM(trk.px(), trk.py(), trk.pz(), MassElectron);
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
      // nSigmaPr[counterTmp] = trk.tpcNSigmaPr();
      nSigmaPr[counterTmp] = std::sqrt(trk.tofNSigmaPr() * trk.tofNSigmaPr() + trk.tpcNSigmaPr() * trk.tpcNSigmaPr());
      //      nSigmaKa[counterTmp] = trk.tpcNSigmaKa();
      nSigmaKa[counterTmp] = std::sqrt(trk.tofNSigmaKa() * trk.tofNSigmaKa() + trk.tpcNSigmaKa() * trk.tpcNSigmaKa());
      // dcaZ[counterTmp] = trk.dcaZ();
      // dcaXY[counterTmp] = trk.dcaXY();
      // chi2TPC[counterTmp] = trk.tpcChi2NCl();
      // chi2ITS[counterTmp] = trk.itsChi2NCl();
      chi2TOF[counterTmp] = trk.tofChi2();
      // nclTPCfind[counterTmp] = trk.tpcNClsFindable();
      nclTPCcrossedRows[counterTmp] = trk.tpcNClsCrossedRows();
      trkHasTof[counterTmp] = trk.hasTOF();
      trkTime[counterTmp] = trk.trackTime();
      trkTimeRes[counterTmp] = trk.trackTimeRes();

      p1.SetXYZM(trk.px(), trk.py(), trk.pz(), MassPiPlus);
      tmpMomentum[counterTmp] = p1.P();
      tmpPt[counterTmp] = p1.Pt();
      tmpDedx[counterTmp] = trk.tpcSignal();
      tmpTofNsigmaEl[counterTmp] = trk.tofNSigmaEl();

      deltaPhiTmp = calculateDeltaPhi(p - p1, p1);
      pi3invMass[counterTmp] = (p - p1).Mag();
      pi3pt[counterTmp] = (p - p1).Pt();
      pi3deltaPhi[counterTmp] = deltaPhiTmp;
      pi3assymav[counterTmp] = (p1.Pt() - (scalarPtsum - p1.Pt()) / 3.) / (p1.Pt() + (scalarPtsum - p1.Pt()) / 3.);
      //      pi3vector[counterTmp] = (p + p1).Pt() / (p - p1).Pt();
      pi3vector[counterTmp] = p.Pt() / (p - p1 - p1).Pt();
      // pi3scalar[counterTmp] = (p.Pt() - p1.Pt()) / (p.Pt() + p1.Pt());
      pi3scalar[counterTmp] = ((p - p1).Pt() - p1.Pt()) / ((p - p1).Pt() + p1.Pt());

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
      trkTimeTot += std::abs(trkTime[iTmpBest] - trkTime[i]);
    }
    trkTimeResTot = std::sqrt(trkTimeRes[0] * trkTimeRes[0] +
                              trkTimeRes[1] * trkTimeRes[1] +
                              trkTimeRes[2] * trkTimeRes[2] +
                              trkTimeRes[3] * trkTimeRes[3]);

    // control histos, max 4 per event, cut0
    for (int i = 0; i < 4; i++) {
      fillControlHistos<0>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
      registry.get<TH2>(HIST("control/cut0/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
      registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut0"))->Fill(tmpMomentum[i], tmpDedx[i]);
      registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut0"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
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
    registry.get<TH2>(HIST("control/cut0/hZNACenergy"))->Fill(energyZNA, energyZNC);

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
        registry.get<TH1>(HIST("control/cut0/hDeltaPhi2GamCoh"))->Fill(calculateDeltaPhi(gammaPair[whichPair][1], gammaPair[whichPair][0]));
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
          registry.get<TH1>(HIST("control/cut20/hDeltaPhi2GamCoh"))->Fill(calculateDeltaPhi(gammaPair[whichPair][1], gammaPair[whichPair][0]));
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
    //
    // electron
    //
    if (counterEl > 0) { // Nelectrons>0, cut20
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(10., 1.);
      registry.get<TH1>(HIST("control/cut20/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut20/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut20/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut20/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut20/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut20"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut20"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut20"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<20>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut20/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut20/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut20/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      // no electron
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: no electron PID among 4 tracks");
      }
      return;
    } // end of Nelectrons check

    //
    // electron with tof hit (cut33)
    //
    if (flagEl[0] * trkHasTof[0] +
          flagEl[1] * trkHasTof[1] +
          flagEl[2] * trkHasTof[2] +
          flagEl[3] * trkHasTof[3] >
        0) { // electron has tof hit cut 33
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(11., 1.);
      // registry.get<TH1>(HIST("control/cut33/hDcaZ"))->Fill(dcaZ[i]);
      // registry.get<TH1>(HIST("control/cut33/hDcaXY"))->Fill(dcaXY[i]);
      // registry.get<TH1>(HIST("control/cut33/hChi2TPC"))->Fill(chi2TPC[i]);
      // registry.get<TH1>(HIST("control/cut33/hChi2ITS"))->Fill(chi2ITS[i]);
      registry.get<TH1>(HIST("control/cut33/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut33/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut33/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut33/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut33/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i]) {
          fillControlHistos<33>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut33/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut33"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut33"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          registry.get<TH1>(HIST("control/cut33/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut33/h3pi1eMass"))->Fill(mass3pi1e[i]);
          // registry.get<TH1>(HIST("control/cut33/hPtSpectrumEl"))->Fill(tmpPt[i]);
        } // only for electron
      }
    } else {
      if (verbose) {
        LOGF(info, "cut33 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
        LOGF(debug, "<tautau13topo> Candidate rejected: no TOF hit for electron");
      }
      return;
    } // end of tof hit for electron

    //
    // pi veto cut21
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] >
        0) { // pi veto, cut21
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(12., 1.);
      registry.get<TH1>(HIST("control/cut21/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut21/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut21/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut21/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut21/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut21"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut21"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          fillControlHistos<21>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut21/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut21/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut21/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by pi PID");
      }
      return;
    } // end of pi veto

    //
    // proton veto
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] >
        0) { // proton veto, cut24
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(13., 1.);
      registry.get<TH1>(HIST("control/cut24/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut24/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut24/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut24/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut24/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut24"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut24"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut24"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<24>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut24/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut24/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut24/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID");
      }
      return;
    } // end of proton veto

    //
    // kaon veto
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] >
        0) { // kaon veto, cut25
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(14., 1.);
      registry.get<TH1>(HIST("control/cut25/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut25/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut25/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut25/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut25/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut25"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut25"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut25"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<25>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut25/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut25/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut25/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by KaPID");
      }
      return;
    } // end of kaon veto

    //
    // number of crossed rows in TPC for electron
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] >
        0) { // Nc-rTPC cut, cut28
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(15., 1.);
      registry.get<TH1>(HIST("control/cut28/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut28/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut28/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut28/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut28/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut28"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut28"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          fillControlHistos<28>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut28/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut28/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut28/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by CR");
      }
      return;
    } // end of TPC crossed rows for electron cut

    //
    // virtal calorimeter
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] >
        0) { // vcal veto, cut22
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(16., 1.);
      registry.get<TH1>(HIST("control/cut22/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut22/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut22/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut22/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut22/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut22"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut22"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //  if (i == j) continue;
          //  registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut22"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<22>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut22/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut22/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut22/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by Vcal");
      }
      return;
    } // end of vcal veto

    //
    // 3pi nsigma cut29
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] >
        0) { // nsigma 3pi cut, cut29
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(17., 1.);
      registry.get<TH1>(HIST("control/cut29/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut29/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut29/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut29/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut29/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut29"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut29"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          fillControlHistos<29>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
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
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by 3piPID");
      }
      return;
    } // end of nsigma 3pi cut

    //
    // IM cut
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] * flagIM[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] * flagIM[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] * flagIM[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] * flagIM[3] >
        0) { // 3pi cut, cut26
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(18., 1.);
      registry.get<TH1>(HIST("control/cut26/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut26/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut26/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut26/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut26/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut26"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut26"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut26"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<26>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut26/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut26/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut26/h3pi1eMass"))->Fill(mass3pi1e[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by IM");
      }
      return;
    } // end of inv mass 3 pi cut

    //
    // at least one pion with tof hit (cut34)
    //
    int otherTOFtracks[4];
    for (int i = 0; i < 4; i++) {
      otherTOFtracks[i] = 0;
      if (flagEl[i] && trkHasTof[i]) {
        for (int j = 0; j < 4; j++) {
          if (i == j)
            continue;
          if (trkHasTof[j]) {
            otherTOFtracks[i]++;
            registry.get<TH1>(HIST("pidTOF/h3piTOFchi2"))->Fill(chi2TOF[j]);
          }
        } // second loop over tracks
      }
    } // first loop over tracks
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] * flagIM[0] * (otherTOFtracks[0] >= 1) +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] * flagIM[1] * (otherTOFtracks[1] >= 1) +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] * flagIM[2] * (otherTOFtracks[2] >= 1) +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] * flagIM[3] * (otherTOFtracks[3] >= 1) >
        0) {                                                                // at lest 1 pi with tof hit, cut34
      registry.get<TH1>(HIST("global/hRecFlag"))->Fill(5 + dgcand.flags()); // reconstruction with upc settings flag
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(19., 1.);
      // registry.get<TH1>(HIST("control/cut34/hDcaZ"))->Fill(dcaZ[i]);
      // registry.get<TH1>(HIST("control/cut34/hDcaXY"))->Fill(dcaXY[i]);
      // registry.get<TH1>(HIST("control/cut34/hChi2TPC"))->Fill(chi2TPC[i]);
      // registry.get<TH1>(HIST("control/cut34/hChi2ITS"))->Fill(chi2ITS[i]);
      registry.get<TH1>(HIST("control/cut34/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut34/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut34/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut34/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut34/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (otherTOFtracks[i] >= 1)) {
          fillControlHistos<34>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut34/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut34"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut34"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          registry.get<TH1>(HIST("control/cut34/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut34/h3pi1eMass"))->Fill(mass3pi1e[i]);
          // registry.get<TH1>(HIST("control/cut34/hPtSpectrumEl"))->Fill(tmpPt[i]);
        } else if (!flagEl[i] && trkHasTof[i]) {
          registry.get<TH1>(HIST("pidTOF/h3piTOFchi2Cut34"))->Fill(chi2TOF[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by lack of TOF hit in 3pi");
      }
      return;
    } // end of at least one pion with tof hit (cut34)

    //
    // skip events with pttot<0.15
    //
    if (pttot < ptTotcut) {
      if (verbose) {
        LOGF(info, "<tautau13topo> Candidate rejected: pt tot is %f", pttot);
      }
      return;
    } else {
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(20., 1.);
      registry.get<TH1>(HIST("control/cut30/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut30/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut30/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut30/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut30/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (otherTOFtracks[i] >= 1)) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut30"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut30"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          fillControlHistos<30>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut30/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut30/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut30/h3pi1eMass"))->Fill(mass3pi1e[i]);
        } else if (!flagEl[i] && trkHasTof[i]) {
          registry.get<TH1>(HIST("pidTOF/h3piTOFchi2Cut30"))->Fill(chi2TOF[i]);
        }
      }
    } // end of pttot<0.15 cut30

    //
    // delta phi
    //
    if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] * flagIM[0] * (otherTOFtracks[0] >= 1) * flagDP[0] +
          flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] * flagIM[1] * (otherTOFtracks[1] >= 1) * flagDP[1] +
          flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] * flagIM[2] * (otherTOFtracks[2] >= 1) * flagDP[2] +
          flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] * flagIM[3] * (otherTOFtracks[3] >= 1) * flagDP[3] >
        0) { // delta phi cut, cut27
      registry.get<TH1>(HIST("global/hEventEff"))->Fill(21., 1.);
      registry.get<TH1>(HIST("control/cut27/h4trkPtTot"))->Fill(pttot);
      registry.get<TH1>(HIST("control/cut27/h4piMass"))->Fill(mass4pi);
      registry.get<TH2>(HIST("control/cut27/h4trkMassVsPt"))->Fill(mass4pi, pttot);
      registry.get<TH1>(HIST("control/cut27/hNtofTrk"))->Fill(nTofTrk);
      registry.get<TH2>(HIST("control/cut27/hZNACenergy"))->Fill(energyZNA, energyZNC);
      for (int i = 0; i < 4; i++) {
        if (flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (otherTOFtracks[i] >= 1) && flagDP[i]) {
          registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut27"))->Fill(tmpMomentum[i], tmpDedx[i]);
          registry.get<TH2>(HIST("pidTOF/hpvsNsigmaElHipCut27"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
          // for (int j = 0; j < 4; j++) {
          //   if (i == j) continue;
          //   registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut27"))->Fill(tmpMomentum[j], tmpDedx[j]);
          // }
          fillControlHistos<27>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
          registry.get<TH2>(HIST("control/cut27/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
          registry.get<TH1>(HIST("control/cut27/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          registry.get<TH1>(HIST("control/cut27/h3pi1eMass"))->Fill(mass3pi1e[i]);
        } else if (!flagEl[i] && trkHasTof[i]) {
          registry.get<TH1>(HIST("pidTOF/h3piTOFchi2Cut27"))->Fill(chi2TOF[i]);
          // LOGF(info, "<tautau13topo> chi2TOF %f", chi2TOF[i]);
        }
      }
    } else {
      if (verbose) {
        LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by Dphi");
      }
      return;
    } // end of Dphi

    // // only 1 electron
    // if (counterTotal == 1) {
    //   registry.get<TH1>(HIST("global/hEventEff"))->Fill(19., 1.);
    //   for (int i = 0; i < 4; i++) {
    //     registry.get<TH1>(HIST("control/cut1/hDcaZ"))->Fill(dcaZ[i]);
    //     registry.get<TH1>(HIST("control/cut1/hDcaXY"))->Fill(dcaXY[i]);
    //     registry.get<TH1>(HIST("control/cut1/hChi2TPC"))->Fill(chi2TPC[i]);
    //     registry.get<TH1>(HIST("control/cut1/hChi2ITS"))->Fill(chi2ITS[i]);
    //     registry.get<TH1>(HIST("control/cut1/hChi2TOF"))->Fill(chi2TOF[i]);
    //     if (flagTotal[i]) {
    //       fillControlHistos<1>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
    //       registry.get<TH2>(HIST("control/cut1/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
    //       registry.get<TH2>(HIST("pidTPC/hpvsdedxElHipCut1"))->Fill(tmpMomentum[i], tmpDedx[i]);
    //       for (int j = 0; j < 4; j++) {
    //         if (i == j)
    //           continue;
    //         registry.get<TH2>(HIST("pidTPC3pi/hpvsdedxElHipCut1"))->Fill(tmpMomentum[j], tmpDedx[j]);
    //       }
    //       registry.get<TH1>(HIST("global/hFinalPtSpectrumEl"))->Fill(tmpPt[i]);
    //       registry.get<TH1>(HIST("control/cut1/hTPCnclsFindable"))->Fill(nclTPCfind[i]);
    //       registry.get<TH1>(HIST("control/cut1/hsigma3Pi"))->Fill(nSigma3Pi[i]);
    //       registry.get<TH1>(HIST("control/cut1/h3pi1eMass"))->Fill(mass3pi1e[i]);
    //       if (verbose) {
    //         LOGF(info, "cut1 trackTime %f, %f, %f, %f Res %f, %f, %f, %f", trkTime[0], trkTime[1], trkTime[2], trkTime[3], trkTimeRes[0], trkTimeRes[1], trkTimeRes[2], trkTimeRes[3]);
    //       }
    //     }
    //   } // end of loop over 4 tracks
    //   registry.get<TH1>(HIST("control/cut1/h4trkPtTot"))->Fill(pttot);
    //   registry.get<TH1>(HIST("control/cut1/h4piMass"))->Fill(mass4pi);
    //   registry.get<TH2>(HIST("control/cut1/h4trkMassVsPt"))->Fill(mass4pi, pttot);
    //   registry.get<TH1>(HIST("control/cut1/hNtofTrk"))->Fill(nTofTrk);
    //   registry.get<TH2>(HIST("control/cut1/hZNACenergy"))->Fill(energyZNA, energyZNC);
    //   // special case invmass 4pi (2,2.3)
    //   // if (mass4pi<2.3 && mass4pi>2) {
    //   // for (int i = 0; i < 4; i++) {
    //   //   registry.get<TH1>(HIST("control/cut1/cut1a/hDcaZ"))->Fill(dcaZ[i]);
    //   //   registry.get<TH1>(HIST("control/cut1/cut1a/hDcaXY"))->Fill(dcaXY[i]);
    //   //   registry.get<TH1>(HIST("control/cut1/cut1a/hChi2TPC"))->Fill(chi2TPC[i]);
    //   //   registry.get<TH1>(HIST("control/cut1/cut1a/hChi2ITS"))->Fill(chi2ITS[i]);
    //   //
    //   //   if (flagTotal[i]) {
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h3piMassComb"))->Fill(pi3invMass[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h3trkPtTot"))->Fill(pi3pt[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/hDeltaPhi13topo"))->Fill(pi3deltaPhi[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h13AssymPt1ProngAver"))->Fill(pi3assymav[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h13Vector"))->Fill(pi3vector[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h13Scalar"))->Fill(pi3scalar[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/h13EtaSum"))->Fill(pi3etasum[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/hTPCnCrossedRows"))->Fill(nclTPCcrossedRows[i]);
    //   //
    //   //     registry.get<TH2>(HIST("control/cut1/cut1a/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/hTPCnclsFindable"))->Fill(nclTPCfind[i]);
    //   //     registry.get<TH1>(HIST("control/cut1/cut1a/hsigma3Pi"))->Fill(nSigma3Pi[i]);
    //   //   }
    //   // }
    //   // registry.get<TH1>(HIST("control/cut1/cut1a/h4trkPtTot"))->Fill(pttot);
    //   // registry.get<TH1>(HIST("control/cut1/cut1a/h4piMass"))->Fill(mass4pi);
    //   // registry.get<TH2>(HIST("control/cut1/cut1a/h4trkMassVsPt"))->Fill(mass4pi, pttot);
    //   // registry.get<TH1>(HIST("control/cut1/cut1a/hNtofTrk"))->Fill(nTofTrk);
    //   // } // end of mass window for 4pi case
    //
    // } else { // more than 1 electron candidate
    //   if (verbose) {
    //     LOGF(debug, "<tautau13topo> Candidate rejected: more than one electron candidate");
    //   }
    // } // end of 1electrons check
  } // end of processDataSG
  // check ntracks-4PVtracks
  // check pt of remaining (ntracks-4PVtracks) tracks

  //
  // basic distributions from MC related to tau, electron and MC particles
  //
  void processSimpleMCSG(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    registryMC.get<TH1>(HIST("globalMC/hMCZvertex"))->Fill(mcCollision.posZ());
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(0., 1.);
    registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(0., 1.);
    registryMC.get<TH1>(HIST("efficiencyMCMu/effiMu"))->Fill(0., 1.);
    registryMC.get<TH1>(HIST("efficiencyMCPi/effiPi"))->Fill(0., 1.);

    // check how many physical primaries
    // int countPrim = 0;
    int countGen = 0;
    int countBoth = 0;
    int countCharged = 0;
    int countChargedFromTau = 0;
    int countTau = 0;

    float etaTau[2];
    float phiTau[2];

    int pionCounter = 0;
    int singlePionIndex = -1;
    int tmpPionIndex = -1;
    //    int tmpPionGlobalIndex=-1;
    //    int singlePionGlobalIndex=-1;
    //    int singleElectronGlobalIndex=-1;
    //    int threePionGlobalIndex[3]={-1,-1,-1};
    //    int threePionIndex[3]={-1,-1,-1};
    //    int motherOfSinglePionIndex=-1;
    //    int motherOfThreePionIndex=-1;
    bool electronFound = false;
    bool muonFound = false;
    bool threePionsFound = false;
    bool singlePionFound = false;
    bool tauInRapidity = true;
    bool partFromTauInEta = true;
    float partPt = 0.;
    //    bool flagE3pi = false;
    //    bool flagMu3pi = false;
    //    bool flagPi3pi = false;
    bool flagElPlusElMinus = false; // e+ = 0, e- =1
    bool flagMuPlusMuMinus = false; // mu+ = 0, mu- =1
    bool flagPiPlusPiMinus = false; // pi+ = 0, pi- =1

    // loop over MC particles
    for (const auto& mcParticle : mcParticles) {
      // primaries
      // if (mcParticle.isPhysicalPrimary()) {
      // countPrim++;
      // }
      //
      // MC particles produced by generator only
      //
      if (mcParticle.producedByGenerator()) {
        countGen++;
        if (mcParticle.isPhysicalPrimary()) {
          countBoth++;
          if (mcParticle.pdgCode() != 22 && std::abs(mcParticle.pdgCode()) != 12 && std::abs(mcParticle.pdgCode()) != 14 && std::abs(mcParticle.pdgCode()) != 16 && mcParticle.pdgCode() != 130 && mcParticle.pdgCode() != 111) {
            countCharged++;

            registryMC.get<TH1>(HIST("globalMC/hMCetaGen"))->Fill(mcParticle.eta());
            registryMC.get<TH1>(HIST("globalMC/hMCphiGen"))->Fill(mcParticle.phi());
            registryMC.get<TH1>(HIST("globalMC/hMCyGen"))->Fill(mcParticle.y());
            registryMC.get<TH1>(HIST("globalMC/hMCptGen"))->Fill(mcParticle.pt());

            if (mcParticle.has_mothers()) {
              auto const& mother = mcParticle.mothers_first_as<aod::McParticles>();
              if (std::abs(mother.pdgCode()) == 15) {
                countChargedFromTau++;
              } // mother is tau
            } // mc particle has mother
          } // veto neutral particles
        } // physicsl primary
      } // generator produced by

      //
      // tau+/-
      //
      if (std::abs(mcParticle.pdgCode()) == 15) { // tau+/-
        countTau++;
        if (countTau <= 2) {
          etaTau[countTau - 1] = mcParticle.eta();
          phiTau[countTau - 1] = mcParticle.phi();
        }

        registryMC.get<TH1>(HIST("tauMC/hMCeta"))->Fill(mcParticle.eta());
        registryMC.get<TH1>(HIST("tauMC/hMCphi"))->Fill(mcParticle.phi());
        registryMC.get<TH1>(HIST("tauMC/hMCy"))->Fill(mcParticle.y());
        registryMC.get<TH1>(HIST("tauMC/hMCpt"))->Fill(mcParticle.pt());
        if (std::abs(mcParticle.y()) > 0.9)
          tauInRapidity = false;
        pionCounter = 0;
        if (mcParticle.has_daughters()) {
          for (const auto& daughter : mcParticle.daughters_as<aod::McParticles>()) {
            // pions from tau
            if (std::abs(daughter.pdgCode()) == 211) { // 211 = pi+
              pionCounter++;
              tmpPionIndex = daughter.index(); // returns index of daughter of tau, not in the event, not in the MC particles
              if (std::abs(daughter.eta()) > 0.9)
                partFromTauInEta = false;
            } // end of pion check
            // electron from tau
            if (std::abs(daughter.pdgCode()) == 11) { // 11 = electron
              if (daughter.pdgCode() == 11)
                flagElPlusElMinus = true;
              registryMC.get<TH1>(HIST("electronMC/hMCeta"))->Fill(daughter.eta());
              registryMC.get<TH1>(HIST("electronMC/hMCphi"))->Fill(daughter.phi());
              registryMC.get<TH1>(HIST("electronMC/hMCy"))->Fill(daughter.y());
              registryMC.get<TH1>(HIST("electronMC/hMCpt"))->Fill(daughter.pt());

              electronFound = !electronFound;
              partPt = static_cast<float>(daughter.pt());
              // singleElectronGlobalIndex = daughter.globalIndex();
              //  LOGF(info,"e pt %f",daughter.pt());
              if (std::abs(daughter.eta()) > 0.9)
                partFromTauInEta = false;
            } // end of electron check
            // muon from tau
            if (std::abs(daughter.pdgCode()) == 13) {
              if (daughter.pdgCode() == 13)
                flagMuPlusMuMinus = true;
              muonFound = !muonFound;
              partPt = static_cast<float>(daughter.pt());
              // LOGF(info,"mu pt %f",daughter.pt());
              if (std::abs(daughter.eta()) > 0.9)
                partFromTauInEta = false;
            } // end of muon check
          } // end of loop over daughters
          if (pionCounter == 3) {
            threePionsFound = true;
          } // end of 3pi check
          if (pionCounter == 1) {
            singlePionFound = true;
            singlePionIndex = tmpPionIndex;
            auto mcPartTmp = mcParticle.daughters_as<aod::McParticles>().begin() + singlePionIndex;
            if (mcPartTmp.pdgCode() == -211)
              flagPiPlusPiMinus = true;
            partPt = static_cast<float>(mcPartTmp.pt());
            // motherOfSinglePionIndex = mcParticle.index();
            if (std::abs(mcPartTmp.eta()) > 0.9)
              partFromTauInEta = false;
            // LOGF(info,"size %d; tau ID %d GID %d (pdg %d); pi ID %d LID %d GID %d, pt %f", mcParticle.size(), motherOfSinglePionIndex, mcParticle.globalIndex(), mcParticle.pdgCode(),singlePionIndex, mcPartTmp.index(), mcPartTmp.globalIndex(), mcPartTmp.pt());
          } // end of 1 pi check
        } // end check tau has daughter
      } // end of tau
    } // end of loop over MC particles
    // LOGF(info,"pt after %f",partPt);

    // tau related things
    if (countTau == 2) {
      registryMC.get<TH1>(HIST("tauMC/hMCdeltaeta"))->Fill(etaTau[0] - etaTau[1]);
      registryMC.get<TH1>(HIST("tauMC/hMCdeltaphi"))->Fill(calculateDeltaPhi(phiTau[0], phiTau[1]) * 180. / o2::constants::math::PI);
    }

    if (threePionsFound && electronFound) {
      // LOGF(info,"3pi + e found");
      registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(1., 1.);
      if (tauInRapidity) {
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(2., 1.);
        if (partFromTauInEta) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(3., 1.);
          registryMC.get<TH1>(HIST("efficiencyMCEl/hpTelec"))->Fill(partPt, 1.);
          // flagE3pi = true;
        } // particles from tau in eta
      } // tau in y
    } // el + 3pi

    if (threePionsFound && muonFound) {
      // LOGF(info,"3pi + mu found");
      registryMC.get<TH1>(HIST("efficiencyMCMu/effiMu"))->Fill(1., 1.);
      if (tauInRapidity) {
        registryMC.get<TH1>(HIST("efficiencyMCMu/effiMu"))->Fill(2., 1.);
        if (partFromTauInEta) {
          registryMC.get<TH1>(HIST("efficiencyMCMu/effiMu"))->Fill(3., 1.);
          registryMC.get<TH1>(HIST("efficiencyMCMu/hpTmuon"))->Fill(partPt, 1.);
          // flagMu3pi = true;
        } // particles from tau in eta
      } // tau in y
    } // el + 3pi

    if (singlePionFound && threePionsFound) {
      // LOGF(info,"3pi + pi found in MC");
      // flagPi3pi = true;
      registryMC.get<TH1>(HIST("efficiencyMCPi/effiPi"))->Fill(1., 1.);
      if (tauInRapidity) {
        registryMC.get<TH1>(HIST("efficiencyMCPi/effiPi"))->Fill(2., 1.);
        if (partFromTauInEta) {
          registryMC.get<TH1>(HIST("efficiencyMCPi/effiPi"))->Fill(3., 1.);
          registryMC.get<TH1>(HIST("efficiencyMCPi/hpTpi"))->Fill(partPt, 1.);
        } // particles from tau in eta
      } // tau in y
    } // el + 3pi

    registryMC.get<TH2>(HIST("globalMC/hMCnPart"))->Fill(mcParticles.size(), 0);
    registryMC.get<TH2>(HIST("globalMC/hMCnPart"))->Fill(countGen, 1);
    registryMC.get<TH2>(HIST("globalMC/hMCnPart"))->Fill(countBoth, 2);
    registryMC.get<TH2>(HIST("globalMC/hMCnPart"))->Fill(countCharged, 3);
    registryMC.get<TH2>(HIST("globalMC/hMCnPart"))->Fill(countChargedFromTau, 4);
    if (countChargedFromTau != 4)
      return;
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(1., 1.);
    if (electronFound && flagElPlusElMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(2., 1.); // e-
    else if (electronFound && !flagElPlusElMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(3., 1.); // e+
    if (muonFound && flagMuPlusMuMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(4., 1.); // mu-
    else if (muonFound && !flagMuPlusMuMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(5., 1.); // mu+
    if (singlePionFound && flagPiPlusPiMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(6., 1.); // pi-
    else if (singlePionFound && !flagPiPlusPiMinus)
      registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(7., 1.); // pi+

    if (!tauInRapidity)
      return;
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(8., 1.);
    if (!partFromTauInEta)
      return;
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(9., 1.);

  } // end of processSimpleMCSG

  void processEfficiencyMCSG(aod::UDMcCollision const& mcCollision,
                             soa::SmallGroups<soa::Join<aod::UDMcCollsLabels, aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>> const& collisions,
                             LabeledTracks const& tracks,
                             aod::UDMcParticles const& mcParticles)
  {
    if (verbose) {
      LOGF(info, "<tautau13topo_MC> GeneratorIDtot %d", mcCollision.generatorsID());
      // below is not implemented in UDMcCollisions
      // LOGF(info,"<tautau13topo_MC> GeneratorIDtot %d, GenID %d, subGenID %d, source %d", mcCollision.generatorsID(), mcCollision.getGeneratorId(), mcCollision.getSubGeneratorId(), mcCollision.getSourceId());
    }
    registry1MC.get<TH1>(HIST("globalMC/hGeneratorID"))->Fill(mcCollision.generatorsID());
    if (!(generatorIDMC < 0)) { // do not check generatorsID process if generatorIDMC < 0
      if (mcCollision.generatorsID() != generatorIDMC)
        return;
    }

    int indexProngMC[4];
    int index1ProngMC = -1;
    bool is1ProngElectronMC = false;
    // bool is1ProngMuonMC = false;
    // bool is1ProngPionMC = false;
    bool is3prong3PiMC = false;
    int motherIndex[4];

    int count = 0;

    bool tauInRapidity = true;
    bool partFromTauInEta = true;

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary()) {
        if (mcParticle.pdgCode() != 22 && std::abs(mcParticle.pdgCode()) != 12 && std::abs(mcParticle.pdgCode()) != 14 && std::abs(mcParticle.pdgCode()) != 16 && mcParticle.pdgCode() != 130 && mcParticle.pdgCode() != 111) {
          if (mcParticle.has_mothers()) {
            auto const& mother = mcParticle.mothers_first_as<aod::UDMcParticles>();
            if (std::abs(mother.pdgCode()) == 15) {
              if (std::abs(rapidity(mother.e(), mother.pz())) > 0.9)
                tauInRapidity = false;
              if (std::abs(eta(mcParticle.px(), mcParticle.py(), mcParticle.pz())) > 0.9)
                partFromTauInEta = false;

              if (std::abs(mcParticle.pdgCode()) == 11) {
                index1ProngMC = mcParticle.index();
                is1ProngElectronMC = true;
              } else if (std::abs(mcParticle.pdgCode()) == 13) {
                index1ProngMC = mcParticle.index();
                // is1ProngMuonMC = true;
              }

              if (count < 4) {
                indexProngMC[count] = mcParticle.index();
                motherIndex[count] = mother.globalIndex();
              } else {
                indexProngMC[3] = mcParticle.index();
                motherIndex[3] = mother.globalIndex();
              }
              count++;
              if (collisions.size() > 0) {
                registryMC.get<TH1>(HIST("globalMCrec/hMCetaGenCol"))->Fill(eta(mcParticle.px(), mcParticle.py(), mcParticle.pz()));
                registryMC.get<TH1>(HIST("globalMCrec/hMCphiGenCol"))->Fill(phi(mcParticle.px(), mcParticle.py()));
                registryMC.get<TH1>(HIST("globalMCrec/hMCyGenCol"))->Fill(rapidity(mcParticle.e(), mcParticle.pz()));
                registryMC.get<TH1>(HIST("globalMCrec/hMCptGenCol"))->Fill(pt(mcParticle.px(), mcParticle.py()));
              }
            } // mother is tau
          } // has mothers
        } // charged particles
      } // end if isPhysicalPrimary

    } // end loop over mcParticle

    if (count != 4)
      return;
    if (!tauInRapidity)
      return;
    if (!partFromTauInEta)
      return;
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(10., 1.); // just to confirm there is exactly the same selection

    if (index1ProngMC < 0) { // pion case + 3pi
      // bool onlyPi = true;
      // // int motherIndex[4];
      // for (int i = 0; i < 4; i++) {
      // auto const& tmpMC = mcParticles.begin() + indexProngMC[i];
      // if (std::abs(tmpMC.pdgCode()) != 211) onlyPi = false;
      // // // mother's check already done in a loop before
      // // // auto const& mother = tmpMC.mothers_first_as<aod::UDMcParticles>();
      // // // motherIndex[i] = mother.globalIndex();
      // // motherIndex[i] = (tmpMC.mothers_first_as<aod::UDMcParticles>()).globalIndex();
      // }
      int motherIndex1Pi = motherIndex[0];
      int motherIndexNew = -1;
      int nDifferences = 0;
      for (int i = 1; i < 4; i++) {
        if (motherIndex1Pi != motherIndex[i]) { // the same mother index
          nDifferences++;
          motherIndexNew = i;
        }
      }
      if (nDifferences == 3)
        index1ProngMC = indexProngMC[0];
      else
        index1ProngMC = indexProngMC[motherIndexNew];
      // is1ProngPionMC = true;
      // if (!onlyPi) LOGF(info, "ERROR: should be 4 pions, but they are not!");
    } // end of special check for pi + 3pi

    int index3ProngMC[3];
    if (index1ProngMC > 0) { // electron or muon case + 3pi
      int index3pi = 0;
      for (int i = 0; i < 4; i++) {
        if (index1ProngMC == indexProngMC[i])
          continue;
        index3ProngMC[index3pi] = indexProngMC[i];
        index3pi++;
      }
    }

    // create 1 prong and 3 prong MC references
    auto const& tmp1ProngMC = mcParticles.begin() + index1ProngMC;
    // LOGF(info,"tmp1ProngMC ID %d, GID %d", tmp1ProngMC.index(), tmp1ProngMC.globalIndex());

    auto const& tmpPion1MC = mcParticles.begin() + index3ProngMC[0];
    auto const& tmpPion2MC = mcParticles.begin() + index3ProngMC[1];
    auto const& tmpPion3MC = mcParticles.begin() + index3ProngMC[2];

    if (std::abs(tmpPion1MC.pdgCode()) == 211 && std::abs(tmpPion2MC.pdgCode()) == 211 && std::abs(tmpPion3MC.pdgCode()) == 211)
      is3prong3PiMC = true;

    //
    // here it comes e+3pi topology in MC
    //
    if (!(is1ProngElectronMC && is3prong3PiMC))
      return;

    // LOGF(info,"ID3pi1 %d, GID3pi1 %d",tmpPion1MC.index(),tmpPion1MC.globalIndex());
    // LOGF(info,"ID3pi2 %d, GID3pi2 %d",tmpPion2MC.index(),tmpPion2MC.globalIndex());
    // LOGF(info,"ID3pi3 %d, GID3pi3 %d",tmpPion3MC.index(),tmpPion3MC.globalIndex());

    auto deltaAlpha1 = deltaAlpha(tmp1ProngMC, tmpPion1MC);
    auto deltaAlpha2 = deltaAlpha(tmp1ProngMC, tmpPion2MC);
    auto deltaAlpha3 = deltaAlpha(tmp1ProngMC, tmpPion3MC);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaEpi"))->Fill(deltaAlpha1);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaEpi"))->Fill(deltaAlpha2);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaEpi"))->Fill(deltaAlpha3);
    //
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiEpi"))->Fill(calculateDeltaPhi(phi(tmp1ProngMC.px(), tmp1ProngMC.py()), phi(tmpPion1MC.px(), tmpPion1MC.py())));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiEpi"))->Fill(calculateDeltaPhi(phi(tmp1ProngMC.px(), tmp1ProngMC.py()), phi(tmpPion2MC.px(), tmpPion2MC.py())));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiEpi"))->Fill(calculateDeltaPhi(phi(tmp1ProngMC.px(), tmp1ProngMC.py()), phi(tmpPion3MC.px(), tmpPion3MC.py())));
    //
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiPipi"))->Fill(calculateDeltaPhi(phi(tmpPion1MC.px(), tmpPion1MC.py()), phi(tmpPion2MC.px(), tmpPion2MC.py())));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiPipi"))->Fill(calculateDeltaPhi(phi(tmpPion1MC.px(), tmpPion1MC.py()), phi(tmpPion3MC.px(), tmpPion3MC.py())));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaPhiPipi"))->Fill(calculateDeltaPhi(phi(tmpPion2MC.px(), tmpPion2MC.py()), phi(tmpPion3MC.px(), tmpPion3MC.py())));

    //
    auto deltaAlphaPi1 = deltaAlpha(tmpPion1MC, tmpPion2MC);
    auto deltaAlphaPi2 = deltaAlpha(tmpPion1MC, tmpPion3MC);
    auto deltaAlphaPi3 = deltaAlpha(tmpPion2MC, tmpPion3MC);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaPiPi"))->Fill(deltaAlphaPi1);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaPiPi"))->Fill(deltaAlphaPi2);
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaAlphaPiPi"))->Fill(deltaAlphaPi3);
    //
    float energyInCone = 0;
    float angleLimit = 0.5;
    if (deltaAlpha1 < angleLimit) {
      energyInCone += pt(tmpPion1MC.px(), tmpPion1MC.py());
    }
    if (deltaAlpha2 < angleLimit) {
      energyInCone += pt(tmpPion2MC.px(), tmpPion2MC.py());
    }
    if (deltaAlpha3 < angleLimit) {
      energyInCone += pt(tmpPion3MC.px(), tmpPion3MC.py());
    }
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCvirtCal"))->Fill(energyInCone);
    //
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCScalar"))->Fill(scalarAsymMC(tmp1ProngMC, tmpPion1MC));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCScalar"))->Fill(scalarAsymMC(tmp1ProngMC, tmpPion2MC));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCScalar"))->Fill(scalarAsymMC(tmp1ProngMC, tmpPion3MC));
    //
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCVector"))->Fill(vectorAsym(tmp1ProngMC, tmpPion1MC));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCVector"))->Fill(vectorAsym(tmp1ProngMC, tmpPion2MC));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCVector"))->Fill(vectorAsym(tmp1ProngMC, tmpPion3MC));

    // add eta phi
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCptEl"))->Fill(pt(tmp1ProngMC.px(), tmp1ProngMC.py()));

    float px3pi = tmpPion1MC.px() + tmpPion2MC.px() + tmpPion3MC.px();
    float py3pi = tmpPion1MC.py() + tmpPion2MC.py() + tmpPion3MC.py();
    float pz3pi = tmpPion1MC.pz() + tmpPion2MC.pz() + tmpPion3MC.pz();
    float en3pi = tmpPion1MC.e() + tmpPion2MC.e() + tmpPion3MC.e();

    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCpt4trk"))->Fill(pt(tmp1ProngMC.px() + px3pi, tmp1ProngMC.py() + py3pi));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCinvmass4pi"))->Fill(invariantMass(tmp1ProngMC.e() + en3pi, tmp1ProngMC.px() + px3pi, tmp1ProngMC.py() + py3pi, tmp1ProngMC.pz() + pz3pi));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCinvmass3pi"))->Fill(invariantMass(en3pi, px3pi, py3pi, pz3pi));
    registryMC.get<TH1>(HIST("efficiencyMCEl/hMCdeltaphi13"))->Fill(calculateDeltaPhi(phi(tmp1ProngMC.px(), tmp1ProngMC.py()), phi(px3pi, py3pi)));

    // reconstructed event
    if (collisions.size() < 1)
      return;
    registryMC.get<TH1>(HIST("globalMC/hMCefficiency"))->Fill(11., 1.); // there is at least 1 collision associated to MC collision
    if (is1ProngElectronMC && is3prong3PiMC) {
      registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(4., 1.);
      if (collisions.size() == 1)
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(5., 1.);

      // event selection flags
      bool zVertexFlag = false;
      bool allInEtaAcceptance = false;
      int nTrkInEtaRange = 0;
      // bool allAbovePtThreshold = false;
      int nTrkAbovePtThreshold = 0;

      int nPVTracks = 0;
      int nGhostTracks = 0;
      int nGhostPVTracks = 0;

      // int nGenTracks=0;
      // int nGenPrimTracks=0;
      // int nGenPVTracks=0;
      // int nGenPrimPVTracks=0;

      int trackId[4]; // local index in collision
      int trackCharge = 0;
      bool trackToMCmatch[4]; // match between data track and corresponding MC particle; true = match, false = track not found in this collision
      int trackMCId[4];       // when MC found, the global index from MC Particle is stored here

      int matchedElIndexToData = -1;
      int gapSide = -2;
      int truegapSide = -2;
      float reconstructedPtElMatchedToMC = -1;

      bool flagGapSideSGP = false;
      bool flagDoubleGap = false;
      bool tracksMatchedToMC = false;

      // FIT checks
      auto bitMin = 16 - mFITvetoWindow; // default is +- 1 bc (1 bit)
      auto bitMax = 16 + mFITvetoWindow;
      bool flagFITveto = false;

      for (const auto& collision : collisions) {
        // FIT flag set
        flagFITveto = false;
        for (auto bit = bitMin; bit <= bitMax; bit++) {
          if (TESTBIT(collision.bbFT0Apf(), bit))
            flagFITveto = true;
          if (TESTBIT(collision.bbFT0Cpf(), bit))
            flagFITveto = true;
          if (useFV0ForVeto && TESTBIT(collision.bbFV0Apf(), bit))
            flagFITveto = true;
          if (useFDDAForVeto && TESTBIT(collision.bbFDDApf(), bit))
            flagFITveto = true;
          if (useFDDCForVeto && TESTBIT(collision.bbFDDCpf(), bit))
            flagFITveto = true;
        } // end of loop over FIT bits

        registry1MC.get<TH1>(HIST("globalMCrec/hRecFlag"))->Fill(collision.flags()); // reconstruction with upc settings flag
        // registry1MC.get<TH1>(HIST("globalMCrec/hOccupancyInTime"))->Fill(collision.occupancyInTime());

        matchedElIndexToData = -1;
        reconstructedPtElMatchedToMC = -1;
        gapSide = collision.gapSide();
        truegapSide = sgSelector.trueGap(collision, cutFV0, cutFT0A, cutFT0C, cutZDC);
        registryMC.fill(HIST("globalMCrec/GapSide"), gapSide);
        registryMC.fill(HIST("globalMCrec/GapSideTrue"), truegapSide);
        // if (gapSide < 0 || gapSide > 2) continue; //old way
        if (gapSide >= 0 && gapSide <= 2)
          flagGapSideSGP = true;

        gapSide = truegapSide;
        // if (flagGapSideSGP) registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(6., 1./collisions.size());
        // if (flagGapSideSGP && tracksMatchedToMC) registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(7., 1./collisions.size()); // with true information

        // if (gapSide != mGapSide) continue; //old way
        if (gapSide == mGapSide)
          flagDoubleGap = true;
        // if (flagDoubleGap) registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(8., 1./collisions.size());
        // if (flagDoubleGap && tracksMatchedToMC) registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(9., 1./collisions.size()); // with true information
        registryMC.get<TH2>(HIST("globalMCrec/hVertexXY"))->Fill(collision.posX(), collision.posY());
        registryMC.get<TH1>(HIST("globalMCrec/hVertexZ"))->Fill(collision.posZ());

        zVertexFlag = true;
        if (std::abs(collision.posZ()) >= zvertexcut)
          zVertexFlag = false;

        auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
        registryMC.get<TH1>(HIST("globalMCrec/hNTracks"))->Fill(groupedTracks.size());
        nPVTracks = 0;
        nGhostTracks = 0;
        nGhostPVTracks = 0;

        for (int i = 0; i < 4; i++) {
          trackToMCmatch[i] = false;
          trackMCId[i] = -1;
        }

        // nGenTracks=0;
        // nGenPrimTracks=0;
        // nGenPVTracks=0;
        // nGenPrimPVTracks=0;

        trackCharge = 0;
        allInEtaAcceptance = false;
        nTrkInEtaRange = 0;
        // allAbovePtThreshold = false;
        nTrkAbovePtThreshold = 0;

        //
        // first loop over grouped Tracks
        //
        for (auto const& track : groupedTracks) {
          // ghost track
          if (!track.has_udMcParticle()) {
            nGhostTracks++;
          }

          if (track.isPVContributor()) {
            if (nPVTracks < 4) {
              trackId[nPVTracks] = track.index();
            }
            nPVTracks++;
            trackCharge += track.sign();
            // if (std::abs(eta(track.px(),track.py(),track.pz())) >= trkEtacut) allInEtaAcceptance = false;
            if (std::abs(eta(track.px(), track.py(), track.pz())) < trkEtacut)
              nTrkInEtaRange++;
            if (track.pt() > 0.1)
              nTrkAbovePtThreshold++;
            // // if (track.tpcNSigmaEl() > -2 && track.tpcNSigmaEl() < 3) atLeast1ElectronPID = true;
            // ptmp.SetXYZM(track.px(), track.py(), track.pz(), mpion);
            // // hPt->Fill(p.Pt());
            // ptot += ptmp;

            if (track.has_udMcParticle()) {
              // LOGF(info, "track ID %d match to MC (1p,3p0,3p1,3p2) (%d, %d, %d, %d)", track.udMcParticle().globalIndex(),  tmp1ProngMC.globalIndex(), tmpPion1MC.globalIndex(), tmpPion2MC.globalIndex(), tmpPion3MC.globalIndex());
              if (nPVTracks < 5) {
                trackMCId[nPVTracks - 1] = track.udMcParticle().globalIndex();
                if (trackMCId[nPVTracks - 1] == tmp1ProngMC.globalIndex()) {
                  matchedElIndexToData = nPVTracks - 1;
                  reconstructedPtElMatchedToMC = track.pt();
                }

                if (trackMCId[nPVTracks - 1] == tmp1ProngMC.globalIndex() ||
                    trackMCId[nPVTracks - 1] == tmpPion1MC.globalIndex() ||
                    trackMCId[nPVTracks - 1] == tmpPion2MC.globalIndex() ||
                    trackMCId[nPVTracks - 1] == tmpPion3MC.globalIndex())
                  trackToMCmatch[nPVTracks - 1] = true; // flag, we have a match data <=> MC
              }

            } else { // end of case where track has MC Particle associated
              // ghost PV track
              nGhostPVTracks++;
            }
          } else { // PV contributor
            if (track.has_udMcParticle()) {
              // LOGF(info,"non-PV trk: pid %4.0d (%f, %f, %f) gid %d pt %f",track.udMcParticle().pdgCode(), track.udMcParticle().vx(),track.udMcParticle().vy(),track.udMcParticle().vz(),track.udMcParticle().globalIndex(),track.pt());
              if (verbose) {
                LOGF(info, "non-PV trk: pid %4.0d gid %d", track.udMcParticle().pdgCode(), track.udMcParticle().globalIndex());
              }
            }
          }
          // if (track.isGlobalTrack()) nGlobalTracks++;
        } // end of loop over tracks
        registryMC.get<TH1>(HIST("globalMCrec/hNTracksPV"))->Fill(nPVTracks);
        registryMC.get<TH1>(HIST("globalMCrec/hNGhostTracks"))->Fill(nGhostTracks);
        registryMC.get<TH1>(HIST("globalMCrec/hNGhostTracksPV"))->Fill(nGhostPVTracks);
        registryMC.get<TH1>(HIST("globalMCrec/hQtot"))->Fill(trackCharge);

        // check whether tracks match to MC particles
        if (trackToMCmatch[0] && trackToMCmatch[1] && trackToMCmatch[2] && trackToMCmatch[3])
          tracksMatchedToMC = true;
        registryMC.get<TH1>(HIST("globalMCrec/hTrackToMCMatch"))->Fill(tracksMatchedToMC);

        if (nTrkInEtaRange >= 4)
          allInEtaAcceptance = true;
        if (nTrkAbovePtThreshold >= 4)
          nTrkAbovePtThreshold = true;

        // zdc information
        float energyZNA = collision.energyCommonZNA();
        float energyZNC = collision.energyCommonZNC();
        // if (energyZNA < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(energyZNA);
        // if (energyZNC < 0) registry.get<TH1>(HIST("global/hZNACenergyTest"))->Fill(energyZNC);
        if (energyZNA < 0)
          energyZNA = -1.;
        if (energyZNC < 0)
          energyZNC = -1.;
        registryMC.get<TH2>(HIST("globalMCrec/hZNACenergy"))->Fill(energyZNA, energyZNC);
        registryMC.get<TH2>(HIST("globalMCrec/hZNACtime"))->Fill(collision.timeZNA(), collision.timeZNC());

        //
        // here analysis event selection comes and track eta phi, pt comparison after that
        //
        // SG producer: flagGapSideSGP ok
        // Double gap: flagDoubleGap ok
        // npvtracks: nPVTracks ok
        // Zvertex: zVertexFlag
        // tracks in eta: allInEtaAcceptance
        // tracks charge : trackCharge
        // MC to data matching: tracksMatchedToMC

        // it is after reconstruction
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec0"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // skip events wrongly reconstructed by SG producernot
        if (!flagGapSideSGP) {
          if (verbose) {
            LOGF(info, "<tautau13topo MC> Candidate rejected: DGproducer flag is %d", gapSide);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(6., 1. / collisions.size());
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(7., 1. / collisions.size());           // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec1"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // skip not Double Gap events
        if (!flagDoubleGap) {
          if (verbose) {
            LOGF(info, "<tautau13topo MC> Candidate rejected: not double gapevent, gap value is %d", gapSide);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(8., 1. / collisions.size());
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(9., 1. / collisions.size());           // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec2"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // // skip events with too few/many tracks
        if (nPVTracks != 4) {
          if (verbose) {
            LOGF(info, "<tautau13topo MC> Candidate rejected: Number of PV contributors is %d", nPVTracks);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(10., 1.);
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(11., 1.);                              // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec3"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // //four reconstructed track with MC track link
        // auto const track1 = groupedTracks.begin()+trackId[0];
        // auto const track2 = groupedTracks.begin()+trackId[1];
        // auto const track3 = groupedTracks.begin()+trackId[2];
        // auto const track4 = groupedTracks.begin()+trackId[3];

        // here comes histos of global1 but in MC
        registryMC.get<TH1>(HIST("global1MCrec/hTrackPVTotCharge"))->Fill(trackCharge);
        registryMC.get<TH1>(HIST("global1MCrec/hVertexZ"))->Fill(collision.posZ());
        registryMC.get<TH1>(HIST("global1MCrec/hNTracks"))->Fill(groupedTracks.size());
        registryMC.get<TH1>(HIST("global1MCrec/hNTracksPV"))->Fill(nPVTracks);

        TLorentzVector p, p1;
        p.SetXYZM(0., 0., 0., 0.);
        TVector3 v1(0, 0, 0);
        TVector3 v2(0, 0, 0);
        float scalarPtsum = 0;
        bool flagVcalPV[4] = {false, false, false, false};
        float deltaphi = 0;
        bool trkHasTof[4] = {false, false, false, false};

        //
        // second loop, only over PV tracks
        //
        for (int i = 0; i < 4; i++) {
          auto const tmptrack = groupedTracks.begin() + trackId[i];
          if (tmptrack.hasTOF())
            trkHasTof[i] = true;
          v1.SetXYZ(tmptrack.px(), tmptrack.py(), tmptrack.pz());
          // second loop to calculate virtual calorimeter
          for (int j = 0; j < 4; j++) {
            if (i == j)
              continue;
            auto const tmptrack2 = groupedTracks.begin() + trackId[j];
            v2.SetXYZ(tmptrack2.px(), tmptrack2.py(), tmptrack2.pz());
            deltaphi = v1.Angle(v2);
            if (deltaphi < minAnglecut) { // default 0.05
              flagVcalPV[i] = true;
            }
          } // end of second loop
          float tmpEtaData = eta(tmptrack.px(), tmptrack.py(), tmptrack.pz());
          float tmpPhiData = phi(tmptrack.px(), tmptrack.py());
          registryMC.get<TH2>(HIST("global1MCrec/hTrackEtaPhiPV"))->Fill(tmpEtaData, tmpPhiData);
          registryMC.get<TH1>(HIST("global1MCrec/hTrackPtPV"))->Fill(tmptrack.pt());
          p1.SetXYZM(v1.X(), v1.Y(), v1.Z(), MassPiPlus); // in case of ghost

          if (trackMCId[i] >= 0) {
            p1.SetXYZM(v1.X(), v1.Y(), v1.Z(), (std::abs(tmptrack.udMcParticle().pdgCode()) == 211 ? MassPiPlus : MassElectron));
            float tmpPt = pt(tmptrack.udMcParticle().px(), tmptrack.udMcParticle().py());
            float tmpEta = eta(tmptrack.udMcParticle().px(), tmptrack.udMcParticle().py(), tmptrack.udMcParticle().pz());
            float tmpPhi = phi(tmptrack.udMcParticle().px(), tmptrack.udMcParticle().py());
            registryMC.get<TH2>(HIST("global1MCrec/hpTGenRecTracksPV"))->Fill(tmptrack.pt(), tmpPt);
            registryMC.get<TH2>(HIST("global1MCrec/hDeltapTGenRecVsRecpTTracksPV"))->Fill(tmptrack.pt() - tmpPt, tmptrack.pt());
            registryMC.get<TH2>(HIST("global1MCrec/hEtaGenRecTracksPV"))->Fill(tmpEtaData, tmpEta);
            registryMC.get<TH2>(HIST("global1MCrec/hDeltaEtaGenRecVsRecpTTracksPV"))->Fill(tmpEtaData - tmpEta, tmptrack.pt());
            registryMC.get<TH2>(HIST("global1MCrec/hPhiGenRecTracksPV"))->Fill(tmpPhiData, tmpPhi);
            registryMC.get<TH2>(HIST("global1MCrec/hDeltaPhiGenRecVsRecpTTracksPV"))->Fill(calculateDeltaPhi(tmpPhiData, tmpPhi), tmptrack.pt());
          } // MC infor exists
          p += p1;
          scalarPtsum += p1.Pt();
        } // end of short loop over tracks

        int nTofTracks = trkHasTof[0] + trkHasTof[1] + trkHasTof[2] + trkHasTof[3];

        // if vz < 10
        if (!zVertexFlag) { // default = 10
          if (verbose) {
            LOGF(info, "<tautau13topo> Candidate rejected: VertexZ is %f", collision.posZ());
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(12., 1.);
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(13., 1.);                              // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec4"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // if eta tracks <0.9 default
        if (!allInEtaAcceptance) {
          if (verbose) {
            LOGF(info, "<tautau13topo> Candidate rejected: Ntrk inside |eta|<0.9 is %d", nTrkInEtaRange);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(14., 1.);
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(15., 1.);                              // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec5"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // if pt of tracks >100 MeV/c
        if (!nTrkAbovePtThreshold) {
          if (verbose) {
            LOGF(info, "<tautau13topo> Candidate rejected: Ntrk with pT >100 MeV/c is %d", nTrkAbovePtThreshold);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(16., 1.);
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(17., 1.);                              // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec6"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        // skip events with net charge != 0
        if (!sameSign) { // opposite sign is signal
          if (trackCharge != 0) {
            if (verbose) {
              LOGF(info, "<tautau13topo> Candidate rejected: Net charge is %d (dgcand %d), while should be 0", trackCharge, collision.netCharge());
            }
            return;
          }
        } else { // same sign is background
          if (trackCharge == 0) {
            if (verbose) {
              LOGF(info, "<tautau13topo> Candidate rejected: Net charge is %d (dgcand %d), while should be not 0", trackCharge, collision.netCharge());
            }
            return;
          }
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(18., 1.);
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(19., 1.);                              // with true information
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec7"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        //
        // n TOF tracks cut 32
        //
        if (nTofTracks < nTofTrkMinCut) {
          if (verbose) {
            LOGF(info, "<tautau13topo> Candidate rejected: nTOFtracks is %d", nTofTracks);
          }
          return;
        }
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(56., 1.); // TOF tracks > Ntoftracks
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(57., 1.);                              // electron identified, tracks match to MC Particles
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec8"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        //
        // check FIT information
        //
        if (mFITvetoFlag) {
          if (flagFITveto) {
            if (verbose) {
              LOGF(info, "<tautau13topo> Candidate rejected: FIT not empty");
            }
            return;
          }
        } // end of check emptiness around given BC in FIT detectors
        registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(53., 1.); // electron identified
        if (tracksMatchedToMC) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(54., 1.);                              // electron identified, tracks match to MC Particles
          registryMC.get<TH1>(HIST("globalMCrec/hPtSpectrumElRec9"))->Fill(reconstructedPtElMatchedToMC); // pt El confirmed with true information
        }

        //
        // here PID from TPC starts to be
        //
        // temporary control variables per event with combinatorics
        float pttot = p.Pt();
        float mass4pi = p.Mag();
        int counterTmp = 0;

        float nSigmaEl[4];
        float nSigmaPi[4];
        float nSigma3Pi[4] = {0., 0., 0., 0.};
        float nSigmaPr[4];
        float nSigmaKa[4];
        // float dcaZ[4];
        // float dcaXY[4];
        // float chi2TPC[4];
        // float chi2ITS[4];
        float chi2TOF[4];
        // float nclTPCfind[4];
        float nclTPCcrossedRows[4];
        // bool tmpHasTOF[4];
        // double trkTime[4];
        // float trkTimeRes[4];

        float tmpMomentum[4];
        float tmpPt[4];
        float tmpDedx[4];
        float tmpTofNsigmaEl[4];

        float deltaPhiTmp = 0;
        float pi3invMass[4];
        float pi3pt[4];
        float pi3deltaPhi[4];
        float pi3assymav[4];
        float pi3vector[4];
        float pi3scalar[4];

        // bool trkHasTof[4] = {false, false, false, false};

        //    float mass3pi1e[4];
        //    double trkTimeTot = 0.;
        //    float trkTimeResTot = 10000.;
        for (int i = 0; i < 4; i++) {
          auto const tmptrack = groupedTracks.begin() + trackId[i];
          // if (tmptrack.hasTOF()) trkHasTof[i] = true;
          v1.SetXYZ(tmptrack.px(), tmptrack.py(), tmptrack.pz());
          p1.SetXYZM(v1.X(), v1.Y(), v1.Z(), MassPiPlus); // in case of ghost
          if (trackMCId[i] >= 0) {
            p1.SetXYZM(v1.X(), v1.Y(), v1.Z(), (i == matchedElIndexToData ? MassElectron : MassPiPlus));
          }

          nSigmaEl[counterTmp] = tmptrack.tpcNSigmaEl();
          nSigmaPi[counterTmp] = tmptrack.tpcNSigmaPi();
          nSigma3Pi[3] += (nSigmaPi[counterTmp] * nSigmaPi[counterTmp]);
          // nSigmaPr[counterTmp] = tmptrack.tpcNSigmaPr();
          // nSigmaPr[counterTmp] = (tmptrack.pt() < 1.5 ? tmptrack.tofNSigmaPr() : tmptrack.tpcNSigmaPr() );
          nSigmaPr[counterTmp] = std::sqrt(tmptrack.tofNSigmaPr() * tmptrack.tofNSigmaPr() + tmptrack.tpcNSigmaPr() * tmptrack.tpcNSigmaPr());
          // nSigmaKa[counterTmp] = tmptrack.tpcNSigmaKa();
          // nSigmaKa[counterTmp] = (tmptrack.pt() < 1.3 ? tmptrack.tofNSigmaKa() : tmptrack.tpcNSigmaKa() );
          nSigmaKa[counterTmp] = std::sqrt(tmptrack.tofNSigmaKa() * tmptrack.tofNSigmaKa() + tmptrack.tpcNSigmaKa() * tmptrack.tpcNSigmaKa());

          // dcaZ[counterTmp] =    tmptrack.dcaZ();
          // dcaXY[counterTmp] =   tmptrack.dcaXY();
          // chi2TPC[counterTmp] = tmptrack.tpcChi2NCl();
          // chi2ITS[counterTmp] = tmptrack.itsChi2NCl();
          chi2TOF[counterTmp] = tmptrack.tofChi2();
          // nclTPCfind[counterTmp] =        tmptrack.tpcNClsFindable();
          nclTPCcrossedRows[counterTmp] = tmptrack.tpcNClsCrossedRows();

          // tmpHasTOF[counterTmp] =  tmptrack.hasTOF();
          // trkTime[counterTmp] =    tmptrack.trackTime();
          // trkTimeRes[counterTmp] = tmptrack.trackTimeRes();

          tmpMomentum[counterTmp] = p1.P();
          tmpPt[counterTmp] = p1.Pt();
          tmpDedx[counterTmp] = tmptrack.tpcSignal();
          tmpTofNsigmaEl[counterTmp] = tmptrack.tofNSigmaEl();

          deltaPhiTmp = calculateDeltaPhi(p - p1, p1);
          pi3invMass[counterTmp] = (p - p1).Mag();
          pi3pt[counterTmp] = (p - p1).Pt();
          pi3deltaPhi[counterTmp] = deltaPhiTmp;
          pi3assymav[counterTmp] = (p1.Pt() - (scalarPtsum - p1.Pt()) / 3.) / (p1.Pt() + (scalarPtsum - p1.Pt()) / 3.);
          pi3vector[counterTmp] = p.Pt() / (p - p1 - p1).Pt();
          pi3scalar[counterTmp] = ((p - p1).Pt() - p1.Pt()) / ((p - p1).Pt() + p1.Pt());

          counterTmp++;
        } // end of second loop over PVtracks

        // fill the histograms with true information
        for (int i = 0; i < 4; i++) {
          nSigma3Pi[i] = nSigma3Pi[3] - (nSigmaPi[i] * nSigmaPi[i]);
          nSigma3Pi[i] = std::sqrt(nSigma3Pi[i]);
          if (i == matchedElIndexToData) {
            fillControlHistosMCtrue<0>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
            registryMC.get<TH2>(HIST("controlMCtrue/cut0/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut0"))->Fill(tmpMomentum[i], tmpDedx[i]);
            registryMC.get<TH1>(HIST("controlMCtrue/cut0/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut0"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            // registryMC.get<TH1>(HIST("control/cut0/h3pi1eMass"))->Fill(mass3pi1e[i]);
          } else { // only for 1prong = electron true
            fillControlHistosMCcomb<0>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
            registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut0"))->Fill(tmpMomentum[i], tmpDedx[i]);
            registryMC.get<TH1>(HIST("controlMCcomb/cut0/hsigma3Pi"))->Fill(nSigma3Pi[i]);
          }
        } // end of loop over PV tracks' informations
        // global variables
        registryMC.get<TH1>(HIST("controlMCtrue/cut0/h3pi1eMass"))->Fill(mass4pi); // 3pi + 1e mass
        registryMC.get<TH1>(HIST("controlMCtrue/cut0/h4trkPtTot"))->Fill(pttot);
        registryMC.get<TH1>(HIST("controlMCtrue/cut0/h4piMass"))->Fill(mass4pi);
        registryMC.get<TH2>(HIST("controlMCtrue/cut0/h4trkMassVsPt"))->Fill(mass4pi, pttot);
        // registryMC.get<TH1>(HIST("controlMCtrue/cut0/hNtofTrk"))->Fill(nTofTrk);
        registryMC.get<TH2>(HIST("controlMCtrue/cut0/hZNACenergy"))->Fill(energyZNA, energyZNC);

        // remove combinatorics
        // bool flagTotal[4] = {false, false, false, false};
        bool flagIM[4] = {false, false, false, false};
        bool flagDP[4] = {false, false, false, false};
        bool flagEl[4] = {false, false, false, false};
        bool flagPi[4] = {false, false, false, false};
        bool flagPr[4] = {false, false, false, false};
        bool flagKa[4] = {false, false, false, false};
        bool flagCR[4] = {false, false, false, false};
        bool flagS3pi[4] = {false, false, false, false};

        for (int i = 0; i < 4; i++) {
          if (pi3invMass[i] < invMass3piMaxcut) { // default should be 1.8
            if (invMass3piSignalRegion) {
              flagIM[i] = true;
            } else {
              flagIM[i] = false;
            }
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut2"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (pi3deltaPhi[i] > deltaPhiMincut) { // default should be 1.5
            flagDP[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut3"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (minNsigmaElcut < nSigmaEl[i] && nSigmaEl[i] < maxNsigmaElcut) { // default (-2,3)
            flagEl[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut4"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (std::abs(nSigmaPi[i]) > maxNsigmaPiVetocut) { // default is 4
            flagPi[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut5"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          // if (tmpPt[i] > minPtEtrkcut) { // 0.25
          //   flagPt[i] = true;
          // } else {
          //   registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut6"))->Fill(tmpMomentum[i], tmpDedx[i]);
          // }

          if (flagVcalPV[i]) {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut7"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (pttot < 0.15) {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut8"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (std::abs(nSigmaPr[i]) > maxNsigmaPrVetocut) { // default is 3
            flagPr[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut9"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (std::abs(nSigmaKa[i]) > maxNsigmaKaVetocut) { // default is 3
            flagKa[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut10"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (nclTPCcrossedRows[i] > nTPCcrossedRowsMinCut) { // default is 50
            flagCR[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut11"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          if (nSigma3Pi[i] < nSigma3piMaxCut) { // default is 5
            flagS3pi[i] = true;
          } else {
            registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut12"))->Fill(tmpMomentum[i], tmpDedx[i]);
          }

          // flagTotal[i] = flagEl[i] && flagPi[i] && flagPt[i] && !flagVcalPV[i] && flagPr[i] && flagKa[i] && flagIM[i] && flagDP[i] && flagCR[i] && flagS3pi[i];
        } // end of loop over 4 tracks

        int counterEl = flagEl[0] + flagEl[1] + flagEl[2] + flagEl[3];

        //
        // draw PID and control histograms
        //

        //
        // Nelectrons in TPC PID nsigma > 0, cut20
        //
        if (counterEl > 0) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(20., 1.); // at least 1 electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(21., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(22., 1.); // electron identified, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut20/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut20/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut20/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut20/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut20/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (flagEl[i] && tracksMatchedToMC && (i == matchedElIndexToData)) { // only for 1prong = electron true
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut20"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<20>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut20/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut20/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut20"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
              // registry1MC.get<TH1>(HIST("controlMCtrue/cut20/hTofChi2El"))->Fill(chi2TOF[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<20>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut20"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut20/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              // registry1MC.get<TH1>(HIST("controlMCcomb/cut20/hTofChi2El"))->Fill(chi2TOF[i]);
            }
          } // end of loop over 4 PV tracks
          // end of electron found cut20
        } else { // no electron
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: no electron PID among 4 tracks");
          }
          return;
        } // end of Nelectrons check

        //
        // electron has TOF hit
        //
        if (flagEl[0] * trkHasTof[0] +
              flagEl[1] * trkHasTof[1] +
              flagEl[2] * trkHasTof[2] +
              flagEl[3] * trkHasTof[3] >
            0) {                                                             // electron has tof hit cut 33
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(59., 1.); // electron identified + TOF hit
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(60., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] &&
              trkHasTof[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(61., 1.); // El PID +TOF, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut33/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut33/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut33/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut33/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut33/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] &&
                trkHasTof[i] && (i == matchedElIndexToData)) { // only for 1prong = electron true
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut33"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<33>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut33/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut33/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut33"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<33>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut33"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut33/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: electron has no tof hit ");
          }
          return;
        } // end of electron has tof hit cut 33

        //
        // electron survived pi veto, cut21
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] >
            0) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(23., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(24., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] &&
              trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(25., 1.); // pion veto, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut21/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut21/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut21/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut21/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut21/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && (i == matchedElIndexToData)) { // only for 1prong = electron true
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut21"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<21>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut21/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut21/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut21"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<21>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut21"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut21/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by pi PID");
          }
          return;
        } // end of pi veto, cut21

        //
        // additional proton veto on electron, cut24
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] >
            0) {                                                             // proton veto, cut24
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(26., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(27., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] && flagPr[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(28., 1.); // proton veto, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut24/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut24/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut24/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut24/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut24/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && flagEl[i] && flagPi[i] && flagPr[i] && (i == matchedElIndexToData)) { // only for 1prong = electron true
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut24"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<24>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut24/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut24/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut24"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<24>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut24"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut24/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by proton PID");
          }
          return;
        } // end of proton veto, cut24

        //
        // additional kaon veto on electron, cut25
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] >
            0) {                                                             // kaon veto, cut25
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(29., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(30., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] &&
              trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] && flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(31., 1.); // kaon veto, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut25/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut25/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut25/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut25/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut25/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut25"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<25>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut25/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut25/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut25"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<25>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut25"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut25/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by K PID");
          }
          return;
        } // end of proton veto, cut25

        //
        // crossd rows in TPC
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] >
            0) {                                                             // Nc-rTPC cut, cut28
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(32., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(33., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && trkHasTof[matchedElIndexToData] &&
              flagEl[matchedElIndexToData] && flagPi[matchedElIndexToData] && flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(34., 1.); // CR, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut28/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut28/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut28/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut28/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut28/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut28"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<28>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut28/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut28/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut28"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<28>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut28"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut28/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by CR in TPC");
          }
          return;
        } // end of TPC crossed rows for electron cut, cut28

        //
        // virtual calorimeter for electron
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] >
            0) {                                                             // vcal veto on electron, cut22
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(35., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(36., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(37., 1.); // kaon veto, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut22/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut22/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut22/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut22/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut22/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut22"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<22>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut22/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut22/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut22"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<22>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut22"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut22/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by Vcal");
          }
          return;
        } // end of vcal cut on electron, cut22

        //
        // nsigma 3pi cut, cut29
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] >
            0) {                                                             // nsigma 3pi cut, cut29
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(41., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(42., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(43., 1.); // kaon veto, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut29/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut29/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut29/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut29/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut29/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut29"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<29>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut29/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut29/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut29"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<29>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut29"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut29/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by 3piSigma");
          }
          return;
        } // end of sigma 3pi cut, cut29

        //
        // invariant mass of 3 pi <1.8
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] * flagIM[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] * flagIM[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] * flagIM[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] * flagIM[3] >
            0) {                                                             // IM 3pi cut, cut26
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(44., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(45., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData] &&
              flagIM[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(46., 1.); // IM 3pi, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut26/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut26/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut26/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut26/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut26/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut26"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<26>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut26/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut26/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut26"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<26>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut26"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut26/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: all electrons vetoed by IM");
          }
          return;
        } // end of IM 3pi cut, cut26

        //
        // at least one pion with tof hit (cut34)
        //
        int otherTOFtracks = 0;
        for (int j = 0; j < 4; j++) {
          if (j == matchedElIndexToData)
            continue;
          otherTOFtracks += trkHasTof[j];
        }

        if (otherTOFtracks >= 1) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(62., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(63., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData] &&
              flagIM[matchedElIndexToData] && trkHasTof[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(64., 1.); // 1tof hit in 3pi, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut34/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut34/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut34/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut34/hZNACenergy"))->Fill(energyZNA, energyZNC);
            registry1MC.get<TH1>(HIST("globalMCrec/hRecFlag"))->Fill(5 + collision.flags()); // reconstruction with upc settings flag
            // registryMC.get<TH1>(HIST("controlMCtrue/cut34/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut34"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<34>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut34/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut34/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut34"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<34>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut34"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut34/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo_MC> Candidate rejected: none of 3pi has no tof hit ");
          }
          return;
        } // end of at least one tof hit in 3pi, cut34

        //
        // skip events with pttot<0.15
        //
        if (pttot >= ptTotcut) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(47., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(48., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData] &&
              flagIM[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(49., 1.); // IM 3pi, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut30/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut30/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut30/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut30/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut30/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut30"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<30>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut30/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut30/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut30"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<30>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut30"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut30/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(info, "<tautau13topo_MC> Candidate rejected: pt tot is %f", pttot);
          }
          return;
        } // end of pttot cut, cut30

        //
        //  delta phi cut, cut27
        //
        if (flagEl[0] * trkHasTof[0] * flagPi[0] * flagPr[0] * flagKa[0] * flagCR[0] * !flagVcalPV[0] * flagS3pi[0] * flagIM[0] * flagDP[0] +
              flagEl[1] * trkHasTof[1] * flagPi[1] * flagPr[1] * flagKa[1] * flagCR[1] * !flagVcalPV[1] * flagS3pi[1] * flagIM[1] * flagDP[1] +
              flagEl[2] * trkHasTof[2] * flagPi[2] * flagPr[2] * flagKa[2] * flagCR[2] * !flagVcalPV[2] * flagS3pi[2] * flagIM[2] * flagDP[2] +
              flagEl[3] * trkHasTof[3] * flagPi[3] * flagPr[3] * flagKa[3] * flagCR[3] * !flagVcalPV[3] * flagS3pi[3] * flagIM[3] * flagDP[3] >
            0) {                                                             // delta phi cut, cut27
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(50., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(51., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && trkHasTof[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData] &&
              flagIM[matchedElIndexToData] && flagDP[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(52., 1.); // IM 3pi, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut27/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut27/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut27/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut27/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut27/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && trkHasTof[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && flagDP[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut27"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<27>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut27/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut27/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut27"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<27>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut27"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut27/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo> Candidate rejected: all electrons vetoed by piPID+Vcal+pT+prPID+KaPID+IM");
          }
          return;
        } // end of dphi 3pi-e cut, cut27

        //
        // ZDC cut Energy < 1 TeV on both sides
        //
        if (energyZNA < 1. && energyZNC < 1.) {
          registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(65., 1.); // electron identified
          if (tracksMatchedToMC)
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(66., 1.); // electron identified, tracks match to MC Particles
          if (tracksMatchedToMC && flagEl[matchedElIndexToData] && flagPi[matchedElIndexToData] &&
              flagPr[matchedElIndexToData] && flagKa[matchedElIndexToData] && flagCR[matchedElIndexToData] &&
              !flagVcalPV[matchedElIndexToData] && flagS3pi[matchedElIndexToData] &&
              flagIM[matchedElIndexToData] && flagDP[matchedElIndexToData] && trkHasTof[matchedElIndexToData]) {
            registryMC.get<TH1>(HIST("efficiencyMCEl/effiEl"))->Fill(67., 1.); // zdc cut, tracks match to MC Particles, El is MC El true

            registryMC.get<TH1>(HIST("controlMCtrue/cut35/h4piMass"))->Fill(mass4pi);
            registryMC.get<TH1>(HIST("controlMCtrue/cut35/h4trkPtTot"))->Fill(pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut35/h4trkMassVsPt"))->Fill(mass4pi, pttot);
            registryMC.get<TH2>(HIST("controlMCtrue/cut35/hZNACenergy"))->Fill(energyZNA, energyZNC);
            // registryMC.get<TH1>(HIST("controlMCtrue/cut35/hNtofTrk"))->Fill(nTofTrk);
          }
          for (int i = 0; i < 4; i++) {
            if (tracksMatchedToMC && flagEl[i] && flagPi[i] && flagPr[i] && flagKa[i] && flagCR[i] && !flagVcalPV[i] && flagS3pi[i] && flagIM[i] && flagDP[i] && trkHasTof[i] && (i == matchedElIndexToData)) {
              registryMC.get<TH2>(HIST("pidTPCMCEltrue/hpvsdedxElHipCut35"))->Fill(tmpMomentum[i], tmpDedx[i]);
              fillControlHistosMCtrue<35>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("controlMCtrue/cut35/h3piMassVsPt"))->Fill(pi3invMass[i], pi3pt[i]);
              registryMC.get<TH1>(HIST("controlMCtrue/cut35/hsigma3Pi"))->Fill(nSigma3Pi[i]);
              registry1MC.get<TH2>(HIST("pidTOFMCEltrue/hpvsNsigmaElHipCut35"))->Fill(tmpMomentum[i], tmpTofNsigmaEl[i]);
            } else if (tracksMatchedToMC && (i != matchedElIndexToData)) {
              fillControlHistosMCcomb<35>(pi3invMass[i], pi3pt[i], pi3deltaPhi[i], pi3assymav[i], pi3vector[i], pi3scalar[i], nclTPCcrossedRows[i], tmpPt[i], chi2TOF[i]);
              registryMC.get<TH2>(HIST("pidTPCMCPitrue/hpvsdedxElHipCut35"))->Fill(tmpMomentum[i], tmpDedx[i]);
              registryMC.get<TH1>(HIST("controlMCcomb/cut35/hsigma3Pi"))->Fill(nSigma3Pi[i]);
            }
          } // end of loop over 4 PV tracks
        } else {
          if (verbose) {
            LOGF(debug, "<tautau13topo> Candidate rejected: zdc cut ");
          }
          return;
        } // end of zdc, cut35

      } // end of loop over collisions
    } // end of electron + 3pi

  } // end of processEfficiencyMCSG

  PROCESS_SWITCH(TauTau13topo, processDataSG, "Run over SG Producer tables in reco level (reconstructed data or MC)", true);
  PROCESS_SWITCH(TauTau13topo, processSimpleMCSG, "Run over SG Producer tables in true level (MC true only)", false);
  PROCESS_SWITCH(TauTau13topo, processEfficiencyMCSG, "Run over SG Producer tables in true and reco level (MC true and reconstructed)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauTau13topo>(cfgc, TaskName{"TauTau13topo"})};
}
