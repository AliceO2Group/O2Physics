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
/// \brief
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  12.07.2022

// #include <algorithm>
// #include <iterator>

// O2 headers
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TF1.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcTauCentralBarrelRL {

  // Global varialbes
  bool isFirstReconstructedCollisions;
  int countCollisions;
  Service<o2::framework::O2DatabasePDG> pdg;
  SGSelector sgSelector;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", true, {"Print general info to terminal; default it true."}};
  Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};
  Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
  Configurable<float> cutMyGapSideFV0{"FV0", 100, "FV0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0A{"FT0A", 200., "FT0A threshold for SG selector"};
  Configurable<float> cutMyGapSideFT0C{"FT0C", 100., "FT0C threshold for SG selector"};
  Configurable<float> cutMyGapSideZDC{"ZDC", 10., "ZDC threshold for SG selector"};
  Configurable<bool> usePIDwithTOF{"usePIDwithTOF", false, {"Determine whether also TOF should be used in testPIDhypothesis"}};
  Configurable<bool> usePIDwithTOFsigmaAfterTPC{"usePIDwithTOFsigmaAfterTPC", true, {"Determine whether cut on TOF n sigma should be used after TPC-based decision in testPIDhypothesis"}};
  Configurable<float> cutMyTPCnSigmaEl{"cutMyTPCnSigmaEl", 3.f, {"n sigma cut on el in absolut values"}};
  Configurable<float> cutMyTPCnSigmaMu{"cutMyTPCnSigmaMu", 3.f, {"n sigma cut on mu in absolut values"}};
  Configurable<float> cutMyTPCnSigmaPi{"cutMyTPCnSigmaPi", 3.f, {"n sigma cut on pi in absolut values"}};
  Configurable<float> cutMyNsigmaTPCPIDselector{"cutMyNsigmaTPCPIDselector", 35.f, {"n sigma TPC cut on all particles in absolut values for testPIDhypothesis"}};
  Configurable<float> cutMyNsigmaTOFPIDselector{"cutMyNsigmaTOFPIDselector", 35.f, {"n sigma TOF cut on all particles in absolut values for testPIDhypothesis"}};
  Configurable<float> cutAvgITSclusterSize{"cutAvgITSclusterSize", 2.05f, {"specific study"}};
  Configurable<float> cutPtAvgITSclusterSize{"cutPtAvgITSclusterSize", 0.7f, {"specific study"}};
  Configurable<bool> cutMyGlobalTracksOnly{"cutMyGlobalTracksOnly", false, {"Applies cut on here defined global tracks"}};
  Configurable<float> cutMyGTptMin{"cutMyGTptMin", 0.1f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTptMax{"cutMyGTptMax", 1e10f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTetaMin{"cutMyGTetaMin", -0.8f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTetaMax{"cutMyGTetaMax", 0.8f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTdcaZmax{"cutMyGTdcaZmax", 2.f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTdcaXYmax{"cutMyGTdcaXYmax", 1e10f, {"MyGlobalTrack cut"}};
  Configurable<bool> cutMyGTdcaXYusePt{"cutMyGTdcaXYusePt", false, {"MyGlobalTrack cut"}};
  Configurable<bool> cutMyHasITS{"cutMyHasITS", true, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyGTitsNClsMin{"cutMyGTitsNClsMin", 1, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTitsChi2NclMax{"cutMyGTitsChi2NclMax", 36.f, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyGTitsHitsRule{"cutMyGTitsHitsRule", 0, {"MyGlobalTrack cut"}};
  Configurable<bool> cutMyHasTPC{"cutMyHasTPC", true, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyGTtpcNClsMin{"cutMyGTtpcNClsMin", 1, {"MyGlobalTrack cut"}};
  Configurable<int> cutMyGTtpcNClsCrossedRowsMin{"cutMyGTtpcNClsCrossedRowsMin", 70, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTtpcNClsCrossedRowsOverNClsMin{"cutMyGTtpcNClsCrossedRowsOverNClsMin", 0.8f, {"MyGlobalTrack cut"}};
  Configurable<float> cutMyGTtpcChi2NclMax{"cutMyGTtpcChi2NclMax", 4.f, {"MyGlobalTrack cut"}};
  Configurable<bool> applyTauEventSelection{"applyTauEventSelection", true, {"Select event"}};
  Configurable<bool> doMainHistos{"doMainHistos", true, {"Fill main histos"}};
  Configurable<bool> doPIDhistos{"doPIDhistos", true, {"Fill PID histos"}};
  Configurable<bool> doTwoTracks{"doTwoTracks", true, {"Define histos for two tracks and allow to fill them"}};
  Configurable<bool> doPionStudy{"doPionStudy", false, {"Define histos for two pions and allow to fill them"}};
  Configurable<bool> doMuonStudy{"doMuonStudy", false, {"Define histos for two muons and allow to fill them"}};
  Configurable<bool> doJpsiMuMuTests{"doJpsiMuMuTests", false, {"Define specific-tests histos for two muons and allow to fill them"}};
  Configurable<bool> doFourTracks{"doFourTracks", false, {"Define histos for four tracks and allow to fill them"}};
  Configurable<bool> doFourTrackPsi2S{"doFourTrackPsi2S", true, {"Define histos for Psi2S into four charged tracks (pi/mu) and allow to fill them"}};
  Configurable<bool> doSixTracks{"doSixTracks", false, {"Define histos for six tracks and allow to fill them"}};

  ConfigurableAxis axisNtracks{"axisNtracks", {30, -0.5, 29.5}, "Number of tracks in collision"};
  ConfigurableAxis axisZvtx{"axisZvtx", {40, -20., 20.}, "Z-vertex position (cm)"};
  ConfigurableAxis axisInvMass{"axisInvMass", {400, 1., 5.}, "Invariant mass (GeV/c^{2})"};
  ConfigurableAxis axisInvMassWide{"axisInvMassWide", {1000, 0., 10.}, "Invariant mass (GeV/c^{2}), wider range"};
  ConfigurableAxis axisMom{"axisMom", {400, 0., 2.}, "Momentum (GeV/c)"};
  ConfigurableAxis axisMomWide{"axisMomWide", {1000, 0., 10.}, "Momentum (GeV/c), wider range"};
  ConfigurableAxis axisMomSigned{"axisMomSigned", {800, -2., 2.}, "Signed momentum (GeV/c)"};
  ConfigurableAxis axisPt{"axisPt", {400, 0., 2.}, "Transversal momentum (GeV/c)"};
  ConfigurableAxis axisPhi{"axisPhi", {64, -2 * o2::constants::math::PI, 2 * o2::constants::math::PI}, "Azimuthal angle (a.y.)"};
  ConfigurableAxis axisModPhi{"axisModPhi", {400, 0., .4}, "Track fmod(#phi,#pi/9)"};
  ConfigurableAxis axisEta{"axisEta", {50, -1.2, 1.2}, "Pseudorapidity (a.u.)"};
  ConfigurableAxis axisRap{"axisRap", {50, -1.2, 1.2}, "Rapidity (a.u.)"};
  ConfigurableAxis axisAcoplanarity{"axisAcoplanarity", {32, 0.0, o2::constants::math::PI}, "Acoplanarity (rad)"};
  ConfigurableAxis axisTPCdEdx{"axisTPCdEdx", {2000, 0., 200.}, "TPC dE/dx (a.u.)"};
  ConfigurableAxis axisTOFsignal{"axisTOFsignal", {2500, -10000., 40000.}, "TOF signal (a.u.)"};
  ConfigurableAxis axisNsigma{"axisNsigma", {200, -10., 10.}, "n sigma"};
  ConfigurableAxis axisDCA{"axisDCA", {100, -0.5, 0.5}, "DCA (cm)"};
  ConfigurableAxis axisAvgITSclsSizes{"axisAvgITSclsSizes", {500, 0., 10.}, "ITS average cluster size"};
  ConfigurableAxis axisITSnCls{"axisITSnCls", {8, -0.5, 7.5}, "ITS n clusters"};
  ConfigurableAxis axisITSchi2{"axisITSchi2", {100, 0, 50}, "UTS chi2"};
  ConfigurableAxis axisTPCnCls{"axisTPCnCls", {165, -0.5, 164.5}, "TPC n clusters"};
  ConfigurableAxis axisTPCxRwsFrac{"axisTPCxRwsFrac", {200, 0.0, 2.0}, "TPC fraction of crossed raws"};
  ConfigurableAxis axisTPCchi2{"axisTPCchi2", {100, 0, 10}, "TPC chi2"};

  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
  using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;
  using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;

  TF1* funcPhiCutL = nullptr;
  TF1* funcPhiCutH = nullptr;

  // init
  void init(InitContext&)
  {
    mySetITShitsRule(cutMyGTitsHitsRule);

    if (verboseInfo)
      printLargeMessage("INIT METHOD");
    countCollisions = 0;
    isFirstReconstructedCollisions = true;

    histos.add("Events/hCountCollisions", ";;Number of  analysed collision (-)", HistType::kTH1D, {{1, 0.5, 1.5}});
    histos.add("Events/UDtableGapSide", ";GapSide value from UD table (-);Number of events (-)", HistType::kTH1D, {{4, -1.5, 2.5}});
    histos.add("Events/TrueGapSideDiffToTableValue", ";Difference trueGapSide from SGselector and gapSide from UD table (-);Number of events (-)", HistType::kTH1D, {{7, -3.5, 3.5}});
    histos.add("Events/hNreconstructedTracks", ";Number of tracks in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedPVGT", ";Number of good track particles from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedNotPVGT", ";Number of good track particles from NOT primary vertex in a collision (-);Number of events (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedPVGTelectrons", ";Number of good track identified electrons from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedPVGTmuons", ";Number of good track identified muons from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedPVGTpions", ";Number of good track identified pions from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hNreconstructedPVGTothers", ";Number of good track NOT identified electron/muon/pion particles from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {axisNtracks});
    histos.add("Events/hChannelsRatio", ";Channels (-);Number of events (-)", HistType::kTH1D, {{10, -0.5, 9.5}});

    histos.add("Tracks/raw/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {axisZvtx});
    histos.add("Tracks/raw/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("Tracks/raw/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("Tracks/raw/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("Tracks/raw/hTrackPtvsModPhi", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
    histos.add("Tracks/raw/hTrackPtvsModPhiTOF", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
    histos.add("Tracks/raw/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {axisEta});
    histos.add("Tracks/raw/hTrackDcaXY", ";Track DCA_{XY} (cm);Number of events (-)", HistType::kTH1D, {axisDCA});
    histos.add("Tracks/raw/hTrackPtvsDcaXY", ";Track #it{p_{T}} (GeV/c);Track DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
    histos.add("Tracks/raw/hTrackDcaZ", ";Track DCA_{Z} (cm);Number of events (-)", HistType::kTH1D, {axisDCA});
    histos.add("Tracks/raw/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH1D, {axisITSnCls});
    histos.add("Tracks/raw/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {axisITSchi2});
    histos.add("Tracks/raw/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/raw/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/raw/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/raw/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {axisTPCxRwsFrac});
    histos.add("Tracks/raw/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {axisTPCchi2});
    histos.add("Tracks/raw/PID/hTPCsignalVsZ", "All tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/raw/PID/hTPCsignalVsP", "All tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/raw/PID/hTPCsignalVsPt", "All tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/raw/PID/hTPCsignalVsEta", "All tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/raw/PID/hTPCsignalVsPhi", "All tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/raw/PID/hTOFsignalVsP", "All tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsZ", "Positively charged tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsP", "Positively charged tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsPt", "Positively charged tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsEta", "Positively charged tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsPhi", "Positively charged tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/raw/PID/PosCharge/hTOFsignalVsP", "Positively charged tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsZ", "Negatively charged tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsP", "Negatively charged tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsPt", "Negatively charged tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsEta", "Negatively charged tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/raw/PID/NegCharge/hTOFsignalVsP", "Negatively charged tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});

    histos.add("Tracks/GoodTrack/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {axisZvtx});
    histos.add("Tracks/GoodTrack/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
    histos.add("Tracks/GoodTrack/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
    histos.add("Tracks/GoodTrack/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
    histos.add("Tracks/GoodTrack/hTrackPtvsModPhi", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
    histos.add("Tracks/GoodTrack/hTrackPtvsModPhiTOF", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
    histos.add("Tracks/GoodTrack/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {axisEta});
    histos.add("Tracks/GoodTrack/hTrackDcaXY", ";Track DCA_{XY} (cm);Number of events (-)", HistType::kTH1D, {axisDCA});
    histos.add("Tracks/GoodTrack/hTrackPtvsDcaXY", ";Track #it{p_{T}} (GeV/c);Track DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
    histos.add("Tracks/GoodTrack/hTrackDcaZ", ";Track DCA_{Z} (cm);Number of events (-)", HistType::kTH1D, {axisDCA});
    histos.add("Tracks/GoodTrack/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH1D, {axisITSnCls});
    histos.add("Tracks/GoodTrack/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {axisITSchi2});
    histos.add("Tracks/GoodTrack/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/GoodTrack/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/GoodTrack/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {axisTPCnCls});
    histos.add("Tracks/GoodTrack/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {axisTPCxRwsFrac});
    histos.add("Tracks/GoodTrack/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {axisTPCchi2});
    histos.add("Tracks/GoodTrack/PID/hTPCsignalVsZ", "All good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/hTPCsignalVsP", "All good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/hTPCsignalVsPt", "All good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/hTPCsignalVsEta", "All good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/hTPCsignalVsPhi", "All good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/hTOFsignalVsP", "All good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsZ", "Positively charged good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP", "Positively charged good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt", "Positively charged good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta", "Positively charged good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi", "Positively charged good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/PosCharge/hTOFsignalVsP", "Positively charged good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsZ", "Negatively charged good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP", "Negatively charged good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt", "Negatively charged good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta", "Negatively charged good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/NegCharge/hTOFsignalVsP", "Negatively charged good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsZ", "Identified electrons;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP", "Identified electrons;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt", "Identified electrons;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta", "Identified electrons;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi", "Identified electrons;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFsignalVsP", "Identified electrons;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaVsP", "Identified electrons;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaVsP", "Identified electrons;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsMu", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsMu", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#mu}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPi", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPi", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsKa", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsKa", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPr", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPr", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsZ", "Identified muons;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsP", "Identified muons;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt", "Identified muons;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta", "Identified muons;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi", "Identified muons;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFsignalVsP", "Identified muons;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaVsP", "Identified muons;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaVsP", "Identified muons;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsEl", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#it{e}}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsEl", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#it{e}}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPi", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPi", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsKa", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsKa", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPr", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPr", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsZ", "Identified pions;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP", "Identified pions;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt", "Identified pions;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta", "Identified pions;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi", "Identified pions;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFsignalVsP", "Identified pions;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaVsP", "Identified pions;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaVsP", "Identified pions;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsEl", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#it{e}}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsEl", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#it{e}}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsMu", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsMu", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#mu}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsKa", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsKa", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsPr", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsPr", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {axisNsigma, axisNsigma});
    histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsZ", "Identified NOT electron/Muon/Pion;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisZvtx, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsP", "Identified NOT electron/Muon/Pion;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsPt", "Identified NOT electron/Muon/Pion;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPt, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsEta", "Identified NOT electron/Muon/Pion;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisEta, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi", "Identified NOT electron/Muon/Pion;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisPhi, axisTPCdEdx});
    histos.add("Tracks/GoodTrack/PID/Others/hTOFsignalVsP", "Identified NOT electron/Muon/Pion;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});

    if (doTwoTracks) {
      histos.add("EventTwoTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassTOF", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassPtCutITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/hDaughtersPvsITSclusterSize", ";Average ITS cluster size;Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
      histos.add("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos", ";Average ITS cluster size x cos(#lambda);Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
      histos.add("EventTwoTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});

      histos.add("EventTwoTracks/TwoElectrons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/TwoElectrons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoElectrons/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoElectrons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/TwoElectrons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/TwoElectrons/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingP", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPwide", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPt", ";Leading #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPhi", ";Leading #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingRapidity", ";Leading #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPvsOtherP", ";Leading #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPwideVsOtherPwide", ";Leading #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPtVsOtherPt", ";Leading #it{p_{T} (GeV/c); Other #it{p_{T} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPhiVsOtherPhi", ";Leading #phi (rad); Other #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingRapVsOtherRap", ";Leading #it{y} (-); Other #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP", ";Leading #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP", ";Other #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP", ";Leading #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});

      histos.add("EventTwoTracks/TwoMuons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/TwoMuons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoMuons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/TwoMuons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/TwoMuons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/TwoMuons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/TwoMuons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/TwoMuons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/TwoMuons/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventTwoTracks/TwoPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/TwoPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoPions/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/TwoPions/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/TwoPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/TwoPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/TwoPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/TwoPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/TwoPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/TwoPions/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/TwoPions/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/TwoPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventTwoTracks/ElectronMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/ElectronMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuon/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/ElectronMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/ElectronMuon/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventTwoTracks/ElectronPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/ElectronPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/ElectronPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/ElectronPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/ElectronPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/ElectronPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/ElectronPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventTwoTracks/MuonPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/MuonPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/MuonPion/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/MuonPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/MuonPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/MuonPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/MuonPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/MuonPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/MuonPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/MuonPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/MuonPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/MuonPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/MuonPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventTwoTracks/ElectronMuPi/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
      histos.add("EventTwoTracks/ElectronMuPi/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPt", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPvsOtherP", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPtVsOtherPt", ";Electron #it{p_{T} (GeV/c); #mu/#pi #it{p_{T} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPhiVsOtherPhi", ";Electron #phi (rad); #mu/#pi #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronRapVsOtherRap", ";Electron #it{y} (-); #mu/#pi #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsEP", ";Electron #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsOP", ";#mu/#pi #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEP", ";Electron #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});

      histos.add("EventTwoTracks/ElectronOther/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
      histos.add("EventTwoTracks/ElectronOther/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventTwoTracks/ElectronOther/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventTwoTracks/ElectronOther/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
      histos.add("EventTwoTracks/ElectronOther/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventTwoTracks/ElectronOther/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronOther/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventTwoTracks/ElectronOther/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventTwoTracks/ElectronOther/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPt", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronOther/hElectronPvsOtherP", ";Electron #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
      histos.add("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt", ";Electron #it{p_{T} (GeV/c); Other #it{p_{T} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi", ";Electron #phi (rad); Other #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
      histos.add("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap", ";Electron #it{y} (-); Other #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsEP", ";Electron #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsOP", ";#it{e}/#mu/#pi #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP", ";Electron #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});

      if (doPionStudy) {
        histos.add("EventTwoTracks/PionsSelection/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWideITS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUSmuMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUSITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLSITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutUS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutLS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutUSmuMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/PionsSelection/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
        histos.add("EventTwoTracks/PionsSelection/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
        histos.add("EventTwoTracks/PionsSelection/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
        histos.add("EventTwoTracks/PionsSelection/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventTwoTracks/PionsSelection/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
        histos.add("EventTwoTracks/PionsSelection/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
        histos.add("EventTwoTracks/PionsSelection/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXY", ";Daughter #it{p_{T}} (GeV/c);Daughter DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXYPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSize", ";Average ITS cluster size;Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
        histos.add("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSizeXcos", ";Average ITS cluster size x cos(#lambda);Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
      }

      if (doMuonStudy) {
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWideITS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutLS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        if (doJpsiMuMuTests) {
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSnegEta", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSposEta", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSnegRap", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSposRap", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap08", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap05", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap03", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNcls70", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNcls100", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNxRws70", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNxRws100", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        }
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutLS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        if (doJpsiMuMuTests) {
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSnegEta", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSposEta", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSnegRap", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSposRap", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap08", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap05", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap03", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNcls70", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNcls100", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNxRws70", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNxRws100", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi6", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi7", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi8", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi9", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi11", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi13", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi14", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi15", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi16", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi6", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi7", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi8", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi9", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi11", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi13", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi14", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi15", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi16", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC6", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC7", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC8", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC9", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC11", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC13", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC14", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC15", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC16", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC1", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC2", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC3", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC4", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC5", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC6", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC7", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC8", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC9", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC10", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC11", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC12", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC13", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC14", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC15", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
          histos.add("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC16", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        }
        histos.add("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutLSITScut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {axisAcoplanarity});
        histos.add("EventTwoTracks/MuonsSelection/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
        histos.add("EventTwoTracks/MuonsSelection/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
        histos.add("EventTwoTracks/MuonsSelection/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventTwoTracks/MuonsSelection/hJpsiPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventTwoTracks/MuonsSelection/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
        histos.add("EventTwoTracks/MuonsSelection/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
        histos.add("EventTwoTracks/MuonsSelection/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMom, axisMom});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {axisMomWide, axisMomWide});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisPt, axisPt});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {axisPhi, axisPhi});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        if (doJpsiMuMuTests) {
          histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut1", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut2", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut3", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut4", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut5", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        }
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        if (doJpsiMuMuTests) {
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut1", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut2", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut3", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut4", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
          histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut5", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        }
        histos.add("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {axisPt, axisModPhi});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {axisRap, axisRap});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXY", ";Daughter #it{p_{T}} (GeV/c);Daughter DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXYPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter DCA_{XY} (cm)", HistType::kTH2D, {axisPt, axisDCA});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSize", ";Average ITS cluster size;Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
        histos.add("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSizeXcos", ";Average ITS cluster size x cos(#lambda);Daughter #it{p} (GeV/c)", HistType::kTH2D, {axisAvgITSclsSizes, axisMomSigned});
        histos.add("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWidePtFitPlot", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWideCS", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsLP", ";Leading #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsOP", ";Other #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsLP", ";Leading #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {axisMom, axisTOFsignal});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
        histos.add("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {axisMom, axisNsigma});
      }
    }

    if (doFourTracks) {
      histos.add("EventFourTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventFourTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventFourTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventFourTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventFourTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventFourTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventFourTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventFourTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventFourTracks/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventFourTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventFourTracks/WithElectron/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventFourTracks/WithElectron/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventFourTracks/WithElectron/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventFourTracks/WithElectron/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventFourTracks/WithElectron/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventFourTracks/WithElectron/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventFourTracks/WithElectron/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventFourTracks/WithElectron/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventFourTracks/WithElectron/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventFourTracks/WithMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventFourTracks/WithMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventFourTracks/WithMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventFourTracks/WithMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventFourTracks/WithMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventFourTracks/WithMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventFourTracks/WithMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventFourTracks/WithMuon/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventFourTracks/WithMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventFourTracks/WithPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventFourTracks/WithPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventFourTracks/WithPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventFourTracks/WithPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventFourTracks/WithPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventFourTracks/WithPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventFourTracks/WithPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventFourTracks/WithPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventFourTracks/WithPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
      if (doFourTrackPsi2S) {
        histos.add("EventFourTracks/MuonsPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventFourTracks/MuonsPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventFourTracks/MuonsPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
        histos.add("EventFourTracks/MuonsPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
        histos.add("EventFourTracks/MuonsPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventFourTracks/MuonsPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
        histos.add("EventFourTracks/MuonsPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
        histos.add("EventFourTracks/MuonsPions/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});

        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
        histos.add("EventFourTracks/Psi2StoMuMuPiPi/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});

        histos.add("EventFourTracks/Psi2StoElElPiPi/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
        histos.add("EventFourTracks/Psi2StoElElPiPi/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      }
    }

    if (doSixTracks) {
      histos.add("EventSixTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventSixTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventSixTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventSixTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventSixTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventSixTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventSixTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventSixTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventSixTracks//hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventSixTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});

      histos.add("EventSixTracks/SixPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMass});
      histos.add("EventSixTracks/SixPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {axisInvMassWide});
      histos.add("EventSixTracks/SixPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMom});
      histos.add("EventSixTracks/SixPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {axisMomWide});
      histos.add("EventSixTracks/SixPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {axisPt});
      histos.add("EventSixTracks/SixPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {axisPhi});
      histos.add("EventSixTracks/SixPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {axisRap});
      histos.add("EventSixTracks/SixPions/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {axisInvMassWide, axisPt});
      histos.add("EventSixTracks/SixPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {axisMom, axisTPCdEdx});
    }

  } // end init

  // run (always called before process :( )
  void run(ProcessingContext&)
  {

    if (verboseInfo)
      printLargeMessage("RUN METHOD");
    printDebugMessage(Form("countCollisions = %d", countCollisions));

  } // end run

  std::vector<std::pair<int8_t, std::set<uint8_t>>> cutMyRequiredITSHits{};

  void mySetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
  {
    // layer 0 corresponds to the the innermost ITS layer
    cutMyRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  }

  void mySetITShitsRule(int matching)
  {
    switch (matching) {
      case 0: // Run3ITSibAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2});
        break;
      case 1: // Run3ITSibTwo
        mySetRequireHitsInITSLayers(2, {0, 1, 2});
        break;
      case 2: // Run3ITSallAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2, 3, 4, 5, 6});
        break;
      case 3: // Run3ITSall7Layers
        mySetRequireHitsInITSLayers(7, {0, 1, 2, 3, 4, 5, 6});
        break;
      default:
        LOG(fatal) << "You chose wrong ITS matching";
        break;
    }
  }

  bool isFulfillsITSHitRequirementsReinstatement(uint8_t itsClusterMap) const
  {
    constexpr uint8_t bit = 1;
    for (auto& itsRequirement : cutMyRequiredITSHits) {
      auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (bit << requiredLayer); });
      if ((itsRequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < itsRequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  template <typename T>
  bool isGlobalTrackReinstatement(T const& track)
  {
    // kInAcceptance copy
    if (track.pt() < cutMyGTptMin || track.pt() > cutMyGTptMax)
      return false;
    if (eta(track.px(), track.py(), track.pz()) < cutMyGTetaMin || eta(track.px(), track.py(), track.pz()) > cutMyGTetaMax)
      return false;
    // kPrimaryTracks
    // GoldenChi2 cut is only for Run 2
    if (abs(track.dcaZ()) > cutMyGTdcaZmax)
      return false;
    if (cutMyGTdcaXYusePt) {
      float maxDCA = 0.0105f + 0.0350f / pow(track.pt(), 1.1f); // ? not sure yet if will be used
      if (abs(track.dcaXY()) > maxDCA)
        return false;
    } else {
      if (abs(track.dcaXY()) > cutMyGTdcaXYmax)
        return false;
    }
    // kQualityTrack
    // TrackType is always 1 as per definition of processed Run3 AO2Ds
    // ITS
    if (cutMyHasITS && !track.hasITS())
      return false; // ITS refit
    if (track.itsNCls() < cutMyGTitsNClsMin)
      return false;
    if (track.itsChi2NCl() > cutMyGTitsChi2NclMax)
      return false;
    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
      return false;
    //  TPC
    if (cutMyHasTPC && !track.hasTPC())
      return false; // TPC refit
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutMyGTtpcNClsMin)
      return false; // tpcNClsFound()
    if (track.tpcNClsCrossedRows() < cutMyGTtpcNClsCrossedRowsMin)
      return false;
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutMyGTtpcNClsCrossedRowsOverNClsMin)
      return false;
    if (track.tpcChi2NCl() > cutMyGTtpcChi2NclMax)
      return false; // TPC chi2

    return true;
  }

  template <typename T>
  bool reinstallRun2JpsiTrackSelection(T const& track)
  {
    // kInAcceptance copy
    if (eta(track.px(), track.py(), track.pz()) < -0.8 || eta(track.px(), track.py(), track.pz()) > 0.8)
      return false;
    // kPrimaryTracks
    if (abs(track.dcaZ()) > 2.0)
      return false;
    float maxDCA = 0.0105f + 0.0350f / pow(track.pt(), 1.1f);
    if (abs(track.dcaXY()) > maxDCA)
      return false;
    // kQualityTrack
    // ITS
    if (!track.hasITS())
      return false; // ITS refit
    if (track.itsChi2NCl() > 36.)
      return false;
    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
      return false;
    // TPC
    if (!track.hasTPC())
      return false; // TPC refit
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < 70)
      return false; // tpcNClsFound()
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < 0.8)
      return false;
    if (track.tpcChi2NCl() > 4.)
      return false;
    // TOF
    if (!track.hasTOF())
      return false;

    return true;
  }

  template <typename C, typename T>
  bool reinstallRun2JpsiEventSelection(C const& collision, T const& trk1, T const& trk2, float rapMother, float aco)
  {
    // tracks
    if (!reinstallRun2JpsiTrackSelection(trk1))
      return false;
    if (!reinstallRun2JpsiTrackSelection(trk2))
      return false;
    if (trk1.sign() * trk2.sign() > 0)
      return false; // opposite sign
    if ((trk1.tpcNSigmaMu() * trk1.tpcNSigmaMu() + trk2.tpcNSigmaMu() * trk2.tpcNSigmaMu()) >
        (trk1.tpcNSigmaEl() * trk1.tpcNSigmaEl() + trk2.tpcNSigmaEl() * trk2.tpcNSigmaEl()))
      return false; // definitely muon than electron
    // event
    if (collision.posZ() > 15.)
      return false;
    if (rapMother < -0.8 || rapMother > 0.8)
      return false;
    if (aco > 4 * o2::constants::math::PI / 5) // max opening angle 144 degrees (I hope, check)
      return false;

    return true;
  }

  float getPhiModN(float phimodn, int sign, int fieldpolarity)
  {
    if (fieldpolarity < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (sign < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    return std::fmod(phimodn, o2::constants::math::PI / 9.0);
  }

  bool isNotCloseToTPCBorder(float phimodn, float trackpt, float cutWidth)
  {

    funcPhiCutL = new TF1("funcPhiCutL", Form("0.06/x+pi/18.0-%.f", cutWidth), 0, 100);
    funcPhiCutH = new TF1("funcPhiCutH", Form("0.1/x+pi/18.0+%.f", cutWidth), 0, 100);

    if (phimodn < funcPhiCutH->Eval(trackpt) && phimodn > funcPhiCutL->Eval(trackpt))
      return false; // reject track

    return true;
  }

  template <typename T>
  int whatPsi2Schannel(T const& trkDaug1, T const& trkDaug2, T const& trkDaug3, T const& trkDaug4, std::vector<int>& vecPIDidx)
  {
    TLorentzVector jpsi, daug[4];
    daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
    daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
    daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
    daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
    // Find index of the two largest values
    std::vector<std::pair<double, double>> vecPts;
    for (int i = 0; i < 4; i++) {
      vecPts.push_back(std::make_pair(static_cast<double>(daug[i].Pt()), i));
    }
    sort(vecPts.begin(), vecPts.end());
    int idx1L = vecPts[vecPts.size() - 1].second;
    int idx2L = vecPts[vecPts.size() - 2].second;
    int idx3L = vecPts[vecPts.size() - 3].second;
    int idx4L = vecPts[vecPts.size() - 4].second;
    // Create the jpsi
    jpsi = daug[idx1L] + daug[idx2L];
    // The two smallest-pT tracks should be pions
    if ((vecPIDidx[idx3L] == P_MUON || vecPIDidx[idx3L] == P_PION) && (vecPIDidx[idx4L] == P_MUON || vecPIDidx[idx4L] == P_PION)) {
      // Branch into Jpsi to mumu and Jpsi to elel
      if ((vecPIDidx[idx1L] == P_MUON || vecPIDidx[idx1L] == P_PION) && (vecPIDidx[idx2L] == P_MUON || vecPIDidx[idx2L] == P_PION)) {
        // Is jpsi mass?
        if (jpsi.M() < 2.9 || jpsi.M() > 3.3)
          return 0; // Not Psi2S
        return 1;
      } else if (vecPIDidx[idx1L] == P_ELECTRON && vecPIDidx[idx2L] == P_ELECTRON) {
        // Is jpsi mass?
        if (jpsi.M() < 2.75 || jpsi.M() > 3.3)
          return 0; // Not Psi2S
        return 2;
      } else {
        return 0; // Not Psi2S
      }
    } else {
      return 0; // Not Psi2S
    }
  }

  template <typename C, typename Ts>
  void fillHistograms(C reconstructedCollision, Ts reconstructedBarrelTracks)
  {

    if (isFirstReconstructedCollisions) {
      isFirstReconstructedCollisions = false;
      if (verboseInfo)
        printLargeMessage("START LOOPING OVER RECONSTRUCTED COLLISIONS");
    }

    histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(1);
    histos.get<TH1>(HIST("Events/hNreconstructedTracks"))->Fill(reconstructedBarrelTracks.size());

    // Loop over tracks without selections
    for (auto& track : reconstructedBarrelTracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/raw/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/raw/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/raw/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/raw/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH1>(HIST("Tracks/raw/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsModPhi"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      if (track.hasTOF())
        histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsModPhiTOF"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      histos.get<TH1>(HIST("Tracks/raw/hTrackDcaXY"))->Fill(track.dcaXY());
      histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsDcaXY"))->Fill(track.pt(), track.dcaXY());
      histos.get<TH1>(HIST("Tracks/raw/hTrackDcaZ"))->Fill(track.dcaZ());
      histos.get<TH1>(HIST("Tracks/raw/ITS/itsNCls"))->Fill(track.itsNCls());
      histos.get<TH1>(HIST("Tracks/raw/ITS/itsChi2NCl"))->Fill(track.itsChi2NCl());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcNClsFindable"))->Fill(track.tpcNClsFindable());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcNClsFound"))->Fill(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcCrossedRows"))->Fill(track.tpcNClsCrossedRows());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcCrossedRowsOverFindableCls"))->Fill((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcChi2NCl"))->Fill(track.tpcChi2NCl());
    } // Loop over tracks without selections

    int countPVGT = 0;
    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTothers = 0;
    int countPVGTpionsSelection = 0;
    int countPVGTmuonsSelection = 0;
    int countTOFtracks = 0;
    int countTPCcls70 = 0;
    int countTPCcls100 = 0;
    int countTPCxRws70 = 0;
    int countTPCxRws100 = 0;
    std::vector<int> vecPVidx;
    std::vector<int> vecPIDidx;
    // Loop over tracks with selections
    for (auto& track : reconstructedBarrelTracks) {
      if (track.isPVContributor() != 1)
        continue;
      if (cutMyGlobalTracksOnly) {
        if (isGlobalTrackReinstatement(track) != 1)
          continue;
      }
      countPVGT++;
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsModPhi"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      if (track.hasTOF())
        histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsModPhiTOF"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackDcaXY"))->Fill(track.dcaXY());
      histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsDcaXY"))->Fill(track.pt(), track.dcaXY());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackDcaZ"))->Fill(track.dcaZ());
      histos.get<TH1>(HIST("Tracks/GoodTrack/ITS/itsNCls"))->Fill(track.itsNCls());
      histos.get<TH1>(HIST("Tracks/GoodTrack/ITS/itsChi2NCl"))->Fill(track.itsChi2NCl());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcNClsFindable"))->Fill(track.tpcNClsFindable());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcNClsFound"))->Fill(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcCrossedRows"))->Fill(track.tpcNClsCrossedRows());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcCrossedRowsOverFindableCls"))->Fill((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcChi2NCl"))->Fill(track.tpcChi2NCl());
      if (track.hasTOF())
        countTOFtracks++;
      if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) > 70)
        countTPCcls70++;
      if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) > 100)
        countTPCcls100++;
      if (track.tpcNClsCrossedRows() > 70)
        countTPCxRws70++;
      if (track.tpcNClsCrossedRows() > 100)
        countTPCxRws100++;
      int hypothesisID = testPIDhypothesis(track, cutMyNsigmaTPCPIDselector, cutMyNsigmaTOFPIDselector, usePIDwithTOF, usePIDwithTOFsigmaAfterTPC);
      vecPIDidx.push_back(hypothesisID);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
        } else {
          countPVGTpions++;
        }
      } else {
        countPVGTothers++;
      }
      if (abs(track.tpcNSigmaPi()) < cutMyTPCnSigmaPi)
        countPVGTpionsSelection++;
      if (abs(track.tpcNSigmaMu()) < cutMyTPCnSigmaMu)
        countPVGTmuonsSelection++;

    } // Loop over tracks with selections

    histos.get<TH1>(HIST("Events/hNreconstructedPVGT"))->Fill(countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedNotPVGT"))->Fill(reconstructedBarrelTracks.size() - countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTelectrons"))->Fill(countPVGTelectrons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTmuons"))->Fill(countPVGTmuons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTpions"))->Fill(countPVGTpions);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTothers"))->Fill(countPVGTothers);

    if (countPVGTselected == 2 && doTwoTracks) {
      TLorentzVector mother, daug[2], motherOfPions, pion[2], motherOfMuons, muon[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      mother = daug[0] + daug[1];
      pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(211), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(211), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      motherOfPions = pion[0] + pion[1];
      muon[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(13), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      muon[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(13), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      motherOfMuons = muon[0] + muon[1];
      auto acoplanarity = calculateAcoplanarity(daug[0].Phi(), daug[1].Phi());
      auto sign = trkDaug1.sign() * trkDaug2.sign();
      bool passAvgITSclsSizesCut = passITSAvgClsSizesLowMomCut(trkDaug1, cutAvgITSclusterSize, cutPtAvgITSclusterSize) && passITSAvgClsSizesLowMomCut(trkDaug2, cutAvgITSclusterSize, cutPtAvgITSclusterSize);
      if (applyTauEventSelection) {
        if (sign > 0)
          return;
        if (acoplanarity > 4 * o2::constants::math::PI / 5)
          return; // max opening angle 144 degrees (I hope, check)
                  //        if (daug[0].Pt() < 0.2 || daug[1].Pt() < 0.2) return;
                  //        if (motherOfPions.M() > 0.55 || motherOfPions.M() < 1.05) return;
                  //        if (!trkDaug1.hasTOF() || !trkDaug2.hasTOF()) return;
      }

      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMass"))->Fill(motherOfPions.M());
      histos.get<TH1>(HIST("EventTwoTracks/hAcoplanarity"))->Fill(acoplanarity);
      histos.get<TH1>(HIST("EventTwoTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      histos.get<TH2>(HIST("EventTwoTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      if (motherOfPions.Pt() < 0.2) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassPtCut"))->Fill(motherOfPions.M());
        if (passAvgITSclsSizesCut) {
          histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassPtCutITScut"))->Fill(motherOfPions.M());
        }
      }
      if (countTOFtracks == 2) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassTOF"))->Fill(motherOfPions.M());
      }
      if (passAvgITSclsSizesCut) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassITScut"))->Fill(motherOfPions.M());
      }
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug1), trkDaug1.sign() * daug[0].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug2), trkDaug2.sign() * daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug1) * getCosLambda(trkDaug1), trkDaug1.sign() * daug[0].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug2) * getCosLambda(trkDaug2), trkDaug2.sign() * daug[1].P());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTelectrons == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(0);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingRapidity"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPvsOtherP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPwideVsOtherPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPtVsOtherPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()), ((daug[0].P() > daug[1].P()) ? daug[1].Pt() : daug[0].Pt()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPhiVsOtherPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()), ((daug[0].P() > daug[1].P()) ? daug[1].Phi() : daug[0].Phi()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingRapVsOtherRap"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()), ((daug[0].P() > daug[1].P()) ? daug[1].Rapidity() : daug[0].Rapidity()));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTmuons == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(1);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTelectrons == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(2);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 2) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(3);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(4);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(5);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWidePtCut"))->Fill(mother.M());
        }
      }
      if ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
        double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hElectronPt"))->Fill(electronPt);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPvsOtherP"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPwideVsOtherPwide"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPtVsOtherPt"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Pt() : daug[0].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPhiVsOtherPhi"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Phi() : daug[1].Phi(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Phi() : daug[0].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronRapVsOtherRap"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Rapidity() : daug[1].Rapidity(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Rapidity() : daug[0].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(12);
      }
      if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
        double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        if (countPVGTelectrons == 2)
          electronPt = (daug[0].Pt() > daug[1].Pt()) ? daug[0].Pt() : daug[1].Pt();
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hElectronPt"))->Fill(electronPt);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPvsOtherP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()), ((daug[0].P() > daug[1].P()) ? daug[1].Pt() : daug[0].Pt()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()), ((daug[0].P() > daug[1].P()) ? daug[1].Phi() : daug[0].Phi()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()), ((daug[0].P() > daug[1].P()) ? daug[1].Rapidity() : daug[0].Rapidity()));
        } else {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPvsOtherP"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Pt() : daug[0].Pt());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Phi() : daug[1].Phi(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Phi() : daug[0].Phi());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap"))->Fill((enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Rapidity() : daug[1].Rapidity(), (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[1].Rapidity() : daug[0].Rapidity());
        }
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(12);
      }
      if (countPVGTpionsSelection == 2 && doPionStudy) {
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMass"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWide"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherP"))->Fill(motherOfPions.P());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPwide"))->Fill(motherOfPions.P());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPt"))->Fill(motherOfPions.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherPhi"))->Fill(motherOfPions.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hMotherRapidity"))->Fill(motherOfPions.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hMotherMassVsPt"))->Fill(motherOfPions.M(), motherOfPions.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersP"))->Fill(pion[0].P(), pion[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPwide"))->Fill(pion[0].P(), pion[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPt"))->Fill(pion[0].Pt(), pion[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPhi"))->Fill(pion[0].Phi(), pion[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersRapidity"))->Fill(pion[0].Rapidity(), pion[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXY"))->Fill(trkDaug1.pt(), trkDaug1.dcaXY());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXY"))->Fill(trkDaug2.pt(), trkDaug2.dcaXY());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhi"))->Fill(pion[0].Pt(), getPhiModN(pion[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhi"))->Fill(pion[1].Pt(), getPhiModN(pion[1].Phi(), trkDaug2.sign(), 1));
        if (motherOfPions.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassPtCut"))->Fill(motherOfPions.M());
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCut"))->Fill(motherOfPions.M());
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhiPtCut"))->Fill(pion[0].Pt(), getPhiModN(pion[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsModPhiPtCut"))->Fill(pion[1].Pt(), getPhiModN(pion[1].Phi(), trkDaug2.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXYPtCut"))->Fill(trkDaug1.pt(), trkDaug1.dcaXY());
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPtvsDcaXYPtCut"))->Fill(trkDaug2.pt(), trkDaug2.dcaXY());
          if (sign < 0) {
            histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUS"))->Fill(motherOfPions.M());
            histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUSmuMass"))->Fill(motherOfMuons.M());
          }
          if (sign > 0)
            histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLS"))->Fill(motherOfPions.M());
          if (countTOFtracks == 2) {
            if (sign < 0) {
              histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutUS"))->Fill(motherOfPions.M());
              histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutUSmuMass"))->Fill(motherOfMuons.M());
            }
            if (sign > 0)
              histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWidePtCutLS"))->Fill(motherOfPions.M());
          }
          if (passAvgITSclsSizesCut) {
            if (sign < 0)
              histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutUSITScut"))->Fill(motherOfPions.M());
            if (sign > 0)
              histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWidePtCutLSITScut"))->Fill(motherOfPions.M());
          }
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhiPtCut"))->Fill(pion[0].Pt(), getPhiModN(pion[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhiPtCut"))->Fill(pion[1].Pt(), getPhiModN(pion[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hasTOF/hInvariantMassWide"))->Fill(motherOfPions.M());
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhi"))->Fill(pion[0].Pt(), getPhiModN(pion[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hasTOF/hDaughtersPtvsModPhi"))->Fill(pion[1].Pt(), getPhiModN(pion[1].Phi(), trkDaug2.sign(), 1));
        }
        if (passAvgITSclsSizesCut) {
          histos.get<TH1>(HIST("EventTwoTracks/PionsSelection/hInvariantMassWideITS"))->Fill(motherOfPions.M());
        }
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug1), trkDaug1.sign() * daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug2), trkDaug2.sign() * daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug1) * getCosLambda(trkDaug1), trkDaug1.sign() * daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/PionsSelection/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug2) * getCosLambda(trkDaug2), trkDaug2.sign() * daug[1].P());
      }
      if (countPVGTmuonsSelection == 2 && doMuonStudy) {
        float phiModNtrk1 = getPhiModN(muon[0].Phi(), trkDaug1.sign(), 1);
        float phiModNtrk2 = getPhiModN(muon[1].Phi(), trkDaug2.sign(), 1);
        float cutPhiModN1 = 0.01;
        float cutPhiModN2 = 0.03;
        float cutPhiModN3 = 0.06;
        float cutPhiModN4 = 0.1;
        float cutPhiModN5 = 0.2;
        float phiPos = 0.;
        float phiNeg = 0.;
        if (trkDaug1.sign() > 0) {
          phiPos = muon[0].Phi();
          phiNeg = muon[1].Phi();
        } else {
          phiPos = muon[1].Phi();
          phiNeg = muon[0].Phi();
        }
        float phiPosTPC = phiPos - o2::math_utils::angle2Alpha(phiPos);
        float phiNegTPC = phiNeg - o2::math_utils::angle2Alpha(phiNeg);
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMass"))->Fill(motherOfMuons.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWide"))->Fill(motherOfMuons.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hMotherP"))->Fill(motherOfMuons.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hMotherPwide"))->Fill(motherOfMuons.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hMotherPt"))->Fill(motherOfMuons.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hMotherPhi"))->Fill(motherOfMuons.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hMotherRapidity"))->Fill(motherOfMuons.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersP"))->Fill(muon[0].P(), muon[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPwide"))->Fill(muon[0].P(), muon[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPt"))->Fill(muon[0].Pt(), muon[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPhi"))->Fill(muon[0].Phi(), muon[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersRapidity"))->Fill(muon[0].Rapidity(), muon[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXY"))->Fill(trkDaug1.pt(), trkDaug1.dcaXY());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXY"))->Fill(trkDaug2.pt(), trkDaug2.dcaXY());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhi"))->Fill(muon[0].Pt(), phiModNtrk1);
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhi"))->Fill(muon[1].Pt(), phiModNtrk2);
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hMotherMassVsPt"))->Fill(motherOfMuons.M(), motherOfMuons.Pt());
        if (motherOfMuons.M() > 2.9 && motherOfMuons.M() < 3.2)
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hJpsiPt"))->Fill(motherOfMuons.Pt());
        if (doJpsiMuMuTests) {
          if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN1))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut1"))->Fill(muon[0].Pt(), phiModNtrk1);
          if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN1))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut1"))->Fill(muon[1].Pt(), phiModNtrk2);
          if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN2))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut2"))->Fill(muon[0].Pt(), phiModNtrk1);
          if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN2))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut2"))->Fill(muon[1].Pt(), phiModNtrk2);
          if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN3))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut3"))->Fill(muon[0].Pt(), phiModNtrk1);
          if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN3))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut3"))->Fill(muon[1].Pt(), phiModNtrk2);
          if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN4))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut4"))->Fill(muon[0].Pt(), phiModNtrk1);
          if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN4))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut4"))->Fill(muon[1].Pt(), phiModNtrk2);
          if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN5))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut5"))->Fill(muon[0].Pt(), phiModNtrk1);
          if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN5))
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiTPCbordersCut5"))->Fill(muon[1].Pt(), phiModNtrk2);
        }
        if (motherOfMuons.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassPtCut"))->Fill(motherOfMuons.M());
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCut"))->Fill(motherOfMuons.M());
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiPtCut"))->Fill(muon[0].Pt(), phiModNtrk1);
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsModPhiPtCut"))->Fill(muon[1].Pt(), phiModNtrk2);
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXYPtCut"))->Fill(trkDaug1.pt(), trkDaug1.dcaXY());
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPtvsDcaXYPtCut"))->Fill(trkDaug2.pt(), trkDaug2.dcaXY());
          if (sign < 0) {
            histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUS"))->Fill(motherOfMuons.M());
            if (doJpsiMuMuTests) {
              if (muon[0].Eta() < 0.0 && muon[1].Eta() < 0.0)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSnegEta"))->Fill(motherOfMuons.M());
              if (muon[0].Eta() > 0.0 && muon[1].Eta() > 0.0)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSposEta"))->Fill(motherOfMuons.M());
              if (motherOfMuons.Rapidity() < 0.0)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSnegRap"))->Fill(motherOfMuons.M());
              if (motherOfMuons.Rapidity() > 0.0)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSposRap"))->Fill(motherOfMuons.M());
              if (std::abs(motherOfMuons.Rapidity()) < 1.2)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap12"))->Fill(motherOfMuons.M());
              if (std::abs(motherOfMuons.Rapidity()) < 1.0)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap10"))->Fill(motherOfMuons.M());
              if (std::abs(motherOfMuons.Rapidity()) < 0.8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap08"))->Fill(motherOfMuons.M());
              if (std::abs(motherOfMuons.Rapidity()) < 0.5)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap05"))->Fill(motherOfMuons.M());
              if (std::abs(motherOfMuons.Rapidity()) < 0.3)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSrap03"))->Fill(motherOfMuons.M());
              if (countTPCcls70 == 2)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNcls70"))->Fill(motherOfMuons.M());
              if (countTPCcls100 == 2)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNcls100"))->Fill(motherOfMuons.M());
              if (countTPCxRws70 == 2)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNxRws70"))->Fill(motherOfMuons.M());
              if (countTPCxRws100 == 2)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcNxRws100"))->Fill(motherOfMuons.M());
              if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN1) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN1))
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut1"))->Fill(motherOfMuons.M());
              if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN2) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN2))
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut2"))->Fill(motherOfMuons.M());
              if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN3) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN3))
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut3"))->Fill(motherOfMuons.M());
              if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN4) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN4))
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut4"))->Fill(motherOfMuons.M());
              if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN5) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN5))
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUStpcBordersCut5"))->Fill(motherOfMuons.M());
              if (-8 * o2::constants::math::PI / 8 <= phiPos && phiPos <= -7 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi1"))->Fill(motherOfMuons.M());
              if (-7 * o2::constants::math::PI / 8 < phiPos && phiPos <= -6 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi2"))->Fill(motherOfMuons.M());
              if (-6 * o2::constants::math::PI / 8 < phiPos && phiPos <= -5 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi3"))->Fill(motherOfMuons.M());
              if (-5 * o2::constants::math::PI / 8 < phiPos && phiPos <= -4 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi4"))->Fill(motherOfMuons.M());
              if (-4 * o2::constants::math::PI / 8 < phiPos && phiPos <= -3 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi5"))->Fill(motherOfMuons.M());
              if (-3 * o2::constants::math::PI / 8 < phiPos && phiPos <= -2 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi6"))->Fill(motherOfMuons.M());
              if (-2 * o2::constants::math::PI / 8 < phiPos && phiPos <= -1 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi7"))->Fill(motherOfMuons.M());
              if (-1 * o2::constants::math::PI / 8 < phiPos && phiPos <= -0 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi8"))->Fill(motherOfMuons.M());
              if (0 * o2::constants::math::PI / 8 < phiPos && phiPos <= 1 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi9"))->Fill(motherOfMuons.M());
              if (1 * o2::constants::math::PI / 8 < phiPos && phiPos <= 2 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi10"))->Fill(motherOfMuons.M());
              if (2 * o2::constants::math::PI / 8 < phiPos && phiPos <= 3 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi11"))->Fill(motherOfMuons.M());
              if (3 * o2::constants::math::PI / 8 < phiPos && phiPos <= 4 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi12"))->Fill(motherOfMuons.M());
              if (4 * o2::constants::math::PI / 8 < phiPos && phiPos <= 5 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi13"))->Fill(motherOfMuons.M());
              if (5 * o2::constants::math::PI / 8 < phiPos && phiPos <= 6 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi14"))->Fill(motherOfMuons.M());
              if (6 * o2::constants::math::PI / 8 < phiPos && phiPos <= 7 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi15"))->Fill(motherOfMuons.M());
              if (7 * o2::constants::math::PI / 8 < phiPos && phiPos <= 8 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhi16"))->Fill(motherOfMuons.M());
              if (-8 * o2::constants::math::PI / 8 <= phiNeg && phiNeg <= -7 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi1"))->Fill(motherOfMuons.M());
              if (-7 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -6 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi2"))->Fill(motherOfMuons.M());
              if (-6 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -5 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi3"))->Fill(motherOfMuons.M());
              if (-5 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -4 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi4"))->Fill(motherOfMuons.M());
              if (-4 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -3 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi5"))->Fill(motherOfMuons.M());
              if (-3 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -2 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi6"))->Fill(motherOfMuons.M());
              if (-2 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -1 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi7"))->Fill(motherOfMuons.M());
              if (-1 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= -0 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi8"))->Fill(motherOfMuons.M());
              if (0 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 1 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi9"))->Fill(motherOfMuons.M());
              if (1 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 2 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi10"))->Fill(motherOfMuons.M());
              if (2 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 3 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi11"))->Fill(motherOfMuons.M());
              if (3 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 4 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi12"))->Fill(motherOfMuons.M());
              if (4 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 5 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi13"))->Fill(motherOfMuons.M());
              if (5 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 6 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi14"))->Fill(motherOfMuons.M());
              if (6 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 7 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi15"))->Fill(motherOfMuons.M());
              if (7 * o2::constants::math::PI / 8 < phiNeg && phiNeg <= 8 * o2::constants::math::PI / 8)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhi16"))->Fill(motherOfMuons.M());

              if (-8 * o2::constants::math::PI / 256 <= phiPosTPC && phiPosTPC <= -7 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC1"))->Fill(motherOfMuons.M());
              if (-7 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -6 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC2"))->Fill(motherOfMuons.M());
              if (-6 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -5 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC3"))->Fill(motherOfMuons.M());
              if (-5 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -4 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC4"))->Fill(motherOfMuons.M());
              if (-4 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -3 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC5"))->Fill(motherOfMuons.M());
              if (-3 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -2 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC6"))->Fill(motherOfMuons.M());
              if (-2 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -1 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC7"))->Fill(motherOfMuons.M());
              if (-1 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= -0 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC8"))->Fill(motherOfMuons.M());
              if (0 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 1 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC9"))->Fill(motherOfMuons.M());
              if (1 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 2 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC10"))->Fill(motherOfMuons.M());
              if (2 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 3 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC11"))->Fill(motherOfMuons.M());
              if (3 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 4 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC12"))->Fill(motherOfMuons.M());
              if (4 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 5 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC13"))->Fill(motherOfMuons.M());
              if (5 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 6 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC14"))->Fill(motherOfMuons.M());
              if (6 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 7 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC15"))->Fill(motherOfMuons.M());
              if (7 * o2::constants::math::PI / 256 < phiPosTPC && phiPosTPC <= 8 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmupPhiTPC16"))->Fill(motherOfMuons.M());
              if (-8 * o2::constants::math::PI / 256 <= phiNegTPC && phiNegTPC <= -7 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC1"))->Fill(motherOfMuons.M());
              if (-7 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -6 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC2"))->Fill(motherOfMuons.M());
              if (-6 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -5 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC3"))->Fill(motherOfMuons.M());
              if (-5 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -4 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC4"))->Fill(motherOfMuons.M());
              if (-4 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -3 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC5"))->Fill(motherOfMuons.M());
              if (-3 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -2 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC6"))->Fill(motherOfMuons.M());
              if (-2 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -1 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC7"))->Fill(motherOfMuons.M());
              if (-1 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= -0 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC8"))->Fill(motherOfMuons.M());
              if (0 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 1 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC9"))->Fill(motherOfMuons.M());
              if (1 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 2 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC10"))->Fill(motherOfMuons.M());
              if (2 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 3 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC11"))->Fill(motherOfMuons.M());
              if (3 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 4 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC12"))->Fill(motherOfMuons.M());
              if (4 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 5 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC13"))->Fill(motherOfMuons.M());
              if (5 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 6 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC14"))->Fill(motherOfMuons.M());
              if (6 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 7 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC15"))->Fill(motherOfMuons.M());
              if (7 * o2::constants::math::PI / 256 < phiNegTPC && phiNegTPC <= 8 * o2::constants::math::PI / 256)
                histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Ruben/hInvariantMassWidePtCutUSmumPhiTPC16"))->Fill(motherOfMuons.M());
            }
          }
          if (sign > 0)
            histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutLS"))->Fill(motherOfMuons.M());
          if (countTOFtracks == 2) {
            if (sign < 0) {
              histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUS"))->Fill(motherOfMuons.M());
              if (doJpsiMuMuTests) {
                if (muon[0].Eta() < 0.0 && muon[1].Eta() < 0.0)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSnegEta"))->Fill(motherOfMuons.M());
                if (muon[0].Eta() > 0.0 && muon[1].Eta() > 0.0)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSposEta"))->Fill(motherOfMuons.M());
                if (motherOfMuons.Rapidity() < 0.0)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSnegRap"))->Fill(motherOfMuons.M());
                if (motherOfMuons.Rapidity() > 0.0)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSposRap"))->Fill(motherOfMuons.M());
                if (std::abs(motherOfMuons.Rapidity()) < 1.2)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap12"))->Fill(motherOfMuons.M());
                if (std::abs(motherOfMuons.Rapidity()) < 1.0)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap10"))->Fill(motherOfMuons.M());
                if (std::abs(motherOfMuons.Rapidity()) < 0.8)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap08"))->Fill(motherOfMuons.M());
                if (std::abs(motherOfMuons.Rapidity()) < 0.5)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap05"))->Fill(motherOfMuons.M());
                if (std::abs(motherOfMuons.Rapidity()) < 0.3)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUSrap03"))->Fill(motherOfMuons.M());
                if (countTPCcls70 == 2)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNcls70"))->Fill(motherOfMuons.M());
                if (countTPCcls100 == 2)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNcls100"))->Fill(motherOfMuons.M());
                if (countTPCxRws70 == 2)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNxRws70"))->Fill(motherOfMuons.M());
                if (countTPCxRws100 == 2)
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcNxRws100"))->Fill(motherOfMuons.M());
                if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN1) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN1))
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut1"))->Fill(motherOfMuons.M());
                if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN2) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN2))
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut2"))->Fill(motherOfMuons.M());
                if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN3) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN3))
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut3"))->Fill(motherOfMuons.M());
                if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN4) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN4))
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut4"))->Fill(motherOfMuons.M());
                if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN5) && isNotCloseToTPCBorder(phiModNtrk2, muon[0].Pt(), cutPhiModN5))
                  histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutUStpcBordersCut5"))->Fill(motherOfMuons.M());
              }
            }
            if (sign > 0)
              histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWidePtCutLS"))->Fill(motherOfMuons.M());
          }
          if (passAvgITSclsSizesCut) {
            if (sign < 0)
              histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutUSITScut"))->Fill(motherOfMuons.M());
            if (sign > 0)
              histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWidePtCutLSITScut"))->Fill(motherOfMuons.M());
          }
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiPtCut"))->Fill(muon[0].Pt(), phiModNtrk1);
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiPtCut"))->Fill(muon[1].Pt(), phiModNtrk2);
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hInvariantMassWide"))->Fill(motherOfMuons.M());
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhi"))->Fill(muon[0].Pt(), phiModNtrk1);
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhi"))->Fill(muon[1].Pt(), phiModNtrk2);
          if (doJpsiMuMuTests) {
            if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN1))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut1"))->Fill(muon[0].Pt(), phiModNtrk1);
            if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN1))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut1"))->Fill(muon[1].Pt(), phiModNtrk2);
            if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN2))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut2"))->Fill(muon[0].Pt(), phiModNtrk1);
            if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN2))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut2"))->Fill(muon[1].Pt(), phiModNtrk2);
            if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN3))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut3"))->Fill(muon[0].Pt(), phiModNtrk1);
            if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN3))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut3"))->Fill(muon[1].Pt(), phiModNtrk2);
            if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN4))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut4"))->Fill(muon[0].Pt(), phiModNtrk1);
            if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN4))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut4"))->Fill(muon[1].Pt(), phiModNtrk2);
            if (isNotCloseToTPCBorder(phiModNtrk1, muon[0].Pt(), cutPhiModN5))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut5"))->Fill(muon[0].Pt(), phiModNtrk1);
            if (isNotCloseToTPCBorder(phiModNtrk2, muon[1].Pt(), cutPhiModN5))
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hasTOF/hDaughtersPtvsModPhiTPCbordersCut5"))->Fill(muon[1].Pt(), phiModNtrk2);
          }
        }
        if (reinstallRun2JpsiEventSelection(reconstructedCollision, trkDaug1, trkDaug2, motherOfMuons.Rapidity(), acoplanarity)) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWide"))->Fill(motherOfMuons.M());
          if (motherOfMuons.Pt() < 2.0) {
            histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWidePtFitPlot"))->Fill(motherOfMuons.M());
          }
          if (motherOfMuons.Pt() < 0.11) {
            histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/Run2Cuts/hInvariantMassWideCS"))->Fill(motherOfMuons.M());
          }
        }
        if (passAvgITSclsSizesCut) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonsSelection/hInvariantMassWideITS"))->Fill(motherOfMuons.M());
        }
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug1), trkDaug1.sign() * daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug2), trkDaug2.sign() * daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug1) * getCosLambda(trkDaug1), trkDaug1.sign() * daug[0].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug2) * getCosLambda(trkDaug2), trkDaug2.sign() * daug[1].P());
      }
    } else if (countPVGTselected == 4 && doFourTracks) {

      TLorentzVector mother, daug[4];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3];

      histos.get<TH1>(HIST("EventFourTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventFourTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventFourTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 4) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(6);
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 3) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(7);
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithElectron/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
      if (countPVGTpions == 3 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(8);
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithMuon/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
      // Hunting down psi2s: broad acceptance
      if (doFourTrackPsi2S) {
        if (countPVGTpions == 4 ||
            (countPVGTpions == 3 && countPVGTmuons == 1) ||
            (countPVGTpions == 2 && countPVGTmuons == 2) ||
            (countPVGTpions == 1 && countPVGTmuons == 3) ||
            countPVGTmuons == 4) {
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hInvariantMass"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hInvariantMassWide"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hMotherP"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hMotherPwide"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hMotherPt"))->Fill(mother.Pt());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hMotherPhi"))->Fill(mother.Phi());
          histos.get<TH1>(HIST("EventFourTracks/MuonsPions/hMotherRapidity"))->Fill(mother.Rapidity());
          histos.get<TH2>(HIST("EventFourTracks/MuonsPions/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        }
        // Hunting down psi2s: ideal case
        if (whatPsi2Schannel(trkDaug1, trkDaug2, trkDaug3, trkDaug4, vecPIDidx) == 1) {
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hInvariantMass"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hInvariantMassWide"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherP"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherPwide"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherPt"))->Fill(mother.Pt());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherPhi"))->Fill(mother.Phi());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherRapidity"))->Fill(mother.Rapidity());
          histos.get<TH2>(HIST("EventFourTracks/Psi2StoMuMuPiPi/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        }
        if (whatPsi2Schannel(trkDaug1, trkDaug2, trkDaug3, trkDaug4, vecPIDidx) == 2) {
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hInvariantMass"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hInvariantMassWide"))->Fill(mother.M());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherP"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherPwide"))->Fill(mother.P());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherPt"))->Fill(mother.Pt());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherPhi"))->Fill(mother.Phi());
          histos.get<TH1>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherRapidity"))->Fill(mother.Rapidity());
          histos.get<TH2>(HIST("EventFourTracks/Psi2StoElElPiPi/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        }
      }
    } else if (countPVGTselected == 6 && doSixTracks) {
      TLorentzVector mother, daug[6];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      const auto& trkDaug5 = reconstructedBarrelTracks.iteratorAt(vecPVidx[4]);
      const auto& trkDaug6 = reconstructedBarrelTracks.iteratorAt(vecPVidx[5]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      daug[4].SetPxPyPzE(trkDaug5.px(), trkDaug5.py(), trkDaug5.pz(), energy(pdg->Mass(trackPDG(trkDaug5)), trkDaug5.px(), trkDaug5.py(), trkDaug5.pz()));
      daug[5].SetPxPyPzE(trkDaug6.px(), trkDaug6.py(), trkDaug6.pz(), energy(pdg->Mass(trackPDG(trkDaug6)), trkDaug6.px(), trkDaug6.py(), trkDaug6.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3] + daug[4] + daug[5];

      histos.get<TH1>(HIST("EventSixTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventSixTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventSixTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 6) {
        histos.get<TH1>(HIST("Events/hChannelsRatio"))->Fill(9);
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventSixTracks/SixPions/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
    } else {
      printDebugMessage("Other particles");
    }

  } // end fillHistograms

  template <typename C, typename Ts>
  void fillPIDhistograms(C /*reconstructedCollision*/, Ts reconstructedBarrelTracks)
  {

    if (isFirstReconstructedCollisions) {
      isFirstReconstructedCollisions = false;
      if (verboseInfo)
        printLargeMessage("START LOOPING OVER RECONSTRUCTED COLLISIONS");
    }

    // Loop over tracks without selections
    for (auto& track : reconstructedBarrelTracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      if (track.hasTPC()) {
        //          histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        if (track.hasTOF())
          histos.get<TH2>(HIST("Tracks/raw/PID/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        if (track.sign() == 1) {
          //                  histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else if (track.sign() == -1) {
          //                  histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
    } // Loop over tracks without selections

    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTmuonsSelection = 0;
    std::vector<int> vecPVidx;
    // Loop over tracks with selections
    for (auto& track : reconstructedBarrelTracks) {
      if (track.isPVContributor() != 1)
        continue;
      if (cutMyGlobalTracksOnly) {
        if (isGlobalTrackReinstatement(track) != 1)
          continue;
      }
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      if (track.hasTPC()) {
        //          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        if (track.hasTOF())
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        if (track.sign() == 1) {
          //                  histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else if (track.sign() == -1) {
          //                  histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
      int hypothesisID = testPIDhypothesis(track, cutMyNsigmaTPCPIDselector, cutMyNsigmaTOFPIDselector, usePIDwithTOF, usePIDwithTOFsigmaAfterTPC);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
          //                  histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsMu"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPi"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsKa"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPr"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsMu"))->Fill(track.tofNSigmaEl(), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPi"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsKa"))->Fill(track.tofNSigmaEl(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPr"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPr());
          }
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
          //                  histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsEl"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPi"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsKa"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPr"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsEl"))->Fill(track.tofNSigmaMu(), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPi"))->Fill(track.tofNSigmaMu(), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsKa"))->Fill(track.tofNSigmaMu(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPr"))->Fill(track.tofNSigmaMu(), track.tofNSigmaPr());
          }
        } else {
          countPVGTpions++;
          //                  histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsEl"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsMu"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsKa"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsPr"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsEl"))->Fill(track.tofNSigmaPi(), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsMu"))->Fill(track.tofNSigmaPi(), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsKa"))->Fill(track.tofNSigmaPi(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsPr"))->Fill(track.tofNSigmaPi(), track.tofNSigmaPr());
          }
        }
      } else {
        //              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        if (track.hasTOF()) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        }
      }
      if (abs(track.tpcNSigmaMu()) < cutMyTPCnSigmaMu)
        countPVGTmuonsSelection++;

    } // Loop over tracks with selections

    if (countPVGTselected == 2 && doTwoTracks) {
      TLorentzVector daug[2], pion[2], muon[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(211), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(211), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      muon[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(13), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      muon[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(13), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (trkDaug1.hasTOF()) {
          histos.get<TH2>(HIST("EventTwoTracks/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
        }
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
        if (countPVGTmuonsSelection == 2 && doMuonStudy) {
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
        if (countPVGTmuons == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (trkDaug2.hasTOF()) {
          histos.get<TH2>(HIST("EventTwoTracks/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
        }
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
        if (countPVGTmuonsSelection == 2 && doMuonStudy) {
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
        if (countPVGTmuons == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
      }
      if (trkDaug1.hasTPC() && trkDaug2.hasTPC()) {
        if (countPVGTelectrons == 2) {
          if (daug[0].P() > daug[1].P()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          } else {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          }
        }
        if (countPVGTmuonsSelection == 2 && doMuonStudy) {
          if (daug[0].P() > daug[1].P()) {
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          } else {
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTPCnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/MuonsSelection/PID/hTOFnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          }
        }
        if ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double electronPID = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double electronNsigma = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double otherPt = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double otherPID = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double otherNsigmaMu = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcNSigmaMu() : trkDaug2.tpcNSigmaMu();
          double otherNsigmaPi = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsEP"))->Fill(electronPt, electronPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsOP"))->Fill(otherPt, otherPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEP"))->Fill(electronPt, electronNsigma);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsMP"))->Fill(otherPt, otherNsigmaMu);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPP"))->Fill(otherPt, otherNsigmaPi);
          if (trkDaug1.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug1.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug1.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug1.tofNSigmaPi());
            }
          }
          if (trkDaug2.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug2.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug2.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug2.tofNSigmaPi());
            }
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          double electronPt = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double electronPID = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double electronNsigma = (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double otherPt = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double otherPID = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double otherNsigmaMu = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcNSigmaMu() : trkDaug2.tpcNSigmaMu();
          double otherNsigmaPi = (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          if (countPVGTelectrons == 2) {
            electronPt = (daug[0].Pt() > daug[1].Pt()) ? daug[0].Pt() : daug[1].Pt();
            electronPID = (daug[0].Pt() > daug[1].Pt()) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
            electronNsigma = (daug[0].Pt() > daug[1].Pt()) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
            otherPt = (daug[0].Pt() > daug[1].Pt()) ? daug[1].Pt() : daug[0].Pt();
            otherPID = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcSignal() : trkDaug1.tpcSignal();
            otherNsigmaMu = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcNSigmaMu() : trkDaug1.tpcNSigmaMu();
            otherNsigmaPi = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcNSigmaPi() : trkDaug1.tpcNSigmaPi();
          }
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsEP"))->Fill(electronPt, electronPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsOP"))->Fill(otherPt, otherPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsEP"))->Fill(electronPt, electronNsigma);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsMP"))->Fill(otherPt, otherNsigmaMu);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsPP"))->Fill(otherPt, otherNsigmaPi);
          if (trkDaug1.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug1)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug1.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug1.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug1.tofNSigmaPi());
            }
          }
          if (trkDaug2.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug2)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug2.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug2.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug2.tofNSigmaPi());
            }
          }
        }
      }
    } else if (countPVGTselected == 4 && doFourTracks) {

      TLorentzVector daug[4];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }
    } else if (countPVGTselected == 6 && doSixTracks) {
      TLorentzVector daug[6];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      const auto& trkDaug5 = reconstructedBarrelTracks.iteratorAt(vecPVidx[4]);
      const auto& trkDaug6 = reconstructedBarrelTracks.iteratorAt(vecPVidx[5]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      daug[4].SetPxPyPzE(trkDaug5.px(), trkDaug5.py(), trkDaug5.pz(), energy(pdg->Mass(trackPDG(trkDaug5)), trkDaug5.px(), trkDaug5.py(), trkDaug5.pz()));
      daug[5].SetPxPyPzE(trkDaug6.px(), trkDaug6.py(), trkDaug6.pz(), energy(pdg->Mass(trackPDG(trkDaug6)), trkDaug6.px(), trkDaug6.py(), trkDaug6.pz()));

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }
      if (trkDaug5.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
      }
      if (trkDaug6.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
      }
    } else {
      printDebugMessage("Other particles");
    }

  } // end fillPIDhistograms

  void processDGrecoLevel(FullUDCollision const& reconstructedCollision,
                          FullUDTracks const& reconstructedBarrelTracks)
  {
    countCollisions++;

    if (doMainHistos)
      fillHistograms(reconstructedCollision, reconstructedBarrelTracks);
    if (doPIDhistos)
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDGrecoLevel

  void processSGrecoLevel(FullSGUDCollision const& reconstructedCollision,
                          FullUDTracks const& reconstructedBarrelTracks)
  {
    countCollisions++;

    int gapSide = reconstructedCollision.gapSide();
    int trueGapSide = sgSelector.trueGap(reconstructedCollision, cutMyGapSideFV0, cutMyGapSideFT0A, cutMyGapSideFT0C, cutMyGapSideZDC);
    histos.fill(HIST("Events/UDtableGapSide"), gapSide);
    histos.fill(HIST("Events/TrueGapSideDiffToTableValue"), gapSide - trueGapSide);
    if (useTrueGap)
      gapSide = trueGapSide;

    if (gapSide != whichGapSide)
      return;

    if (doMainHistos)
      fillHistograms(reconstructedCollision, reconstructedBarrelTracks);
    if (doPIDhistos)
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDGrecoLevel

  void processAnalysisFinished(FullUDCollision const&)
  {

    if (verboseInfo)
      printLargeMessage("END OF THIS DATAFRAME");
    printDebugMessage(Form("countCollisions = %d", countCollisions));
    isFirstReconstructedCollisions = true;

  } // end processAnalysisFinished

  PROCESS_SWITCH(UpcTauCentralBarrelRL, processDGrecoLevel, "Iterate UD tables with reconstructed data created by DG-Candidate-Producer", false);
  PROCESS_SWITCH(UpcTauCentralBarrelRL, processSGrecoLevel, "Iterate UD tables with reconstructed data created by SG-Candidate-Producer", false);
  PROCESS_SWITCH(UpcTauCentralBarrelRL, processAnalysisFinished, "Simply runs in the end of the dataframe", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTauCentralBarrelRL>(cfgc, TaskName{"upc-tau-rl"})};
}
