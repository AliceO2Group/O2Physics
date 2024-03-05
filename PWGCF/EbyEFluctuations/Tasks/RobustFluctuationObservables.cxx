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
/// \brief This task is an QA task to accumulate basic event- and track-level plots.
/// \author Igor Altsybeev, Igor.Altsybeev@cern.ch

#include <iostream>
#include <string>

#include "TF1.h"
#include "TGraphErrors.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "Common/DataModel/FT0Corrected.h"
#include "DataFormatsFT0/Digit.h"

#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct RobustFluctuationObservables {
  // for vertex vs time:
  bool flagShowInfo = false;
  int lastRunNumber = -1;
  int nBCsPerOrbit = 3564;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};

  int64_t bcSOR = -1;    // global bc of the start of the first orbit
  int64_t bcSORbis = -1; // global bc of the start of the first orbit - try alternative
  int64_t nBCsPerTF = 1; // 128*3564; // duration of TF in bcs

  //
  TF1* fPhiCutExpPosHigh;
  TF1* fPhiCutExpPosLow;
  TF1* fPhiCutExpNegHigh;
  TF1* fPhiCutExpNegLow;
  double constPhiShift = 0.175;
  double ptTPCsectorCut = 0.4;

  HistogramRegistry histosEvent{"histosEventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventCounters{"histosEventCounters", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventBcInTF{"histosEventSelectionBcInTF", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosFT0{"histosFT0", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosTracks{"histosTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosK0S{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ##### configurables
  Configurable<int> nBinsPt{"nBinsPt", 800, "N bins in pT histo"};
  Configurable<int> nBinsEta{"nBinsEta", 400, "N bins in eta histo"};
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "N bins in phi histo"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Accepted z-vertex range (cm)"};
  Configurable<int> flagPbPb{"flagPbPb", 0, "0 - pp, 1 - PbPb"};

  // orbit QA
  uint32_t orbitAtCollIndexZero = 0;

  // DF QA
  map<uint64_t, uint64_t> mapDF;

  OutputObj<TH1D> h2D_Orbit_vs_CollIndex_0{TH1D("h2D_Orbit_vs_CollIndex_DF0", "h2D_Orbit_vs_CollIndex_DF0;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_Orbit_vs_CollIndex_1{TH1D("h2D_Orbit_vs_CollIndex_DF1", "h2D_Orbit_vs_CollIndex_DF1;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_Orbit_vs_CollIndex_2{TH1D("h2D_Orbit_vs_CollIndex_DF2", "h2D_Orbit_vs_CollIndex_DF2;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_Orbit_vs_CollIndex_3{TH1D("h2D_Orbit_vs_CollIndex_DF3", "h2D_Orbit_vs_CollIndex_DF3;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_Orbit_vs_CollIndex_4{TH1D("h2D_Orbit_vs_CollIndex_DF4", "h2D_Orbit_vs_CollIndex_DF4;collision index;orbit", 10001, -0.5, 10000.5)};

  OutputObj<TH1D> h2D_BC_vs_CollIndex_0{TH1D("h2D_BC_vs_CollIndex_DF0", "h2D_BC_vs_CollIndex_DF0;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_BC_vs_CollIndex_1{TH1D("h2D_BC_vs_CollIndex_DF1", "h2D_BC_vs_CollIndex_DF1;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_BC_vs_CollIndex_2{TH1D("h2D_BC_vs_CollIndex_DF2", "h2D_BC_vs_CollIndex_DF2;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_BC_vs_CollIndex_3{TH1D("h2D_BC_vs_CollIndex_DF3", "h2D_BC_vs_CollIndex_DF3;collision index;orbit", 10001, -0.5, 10000.5)};
  OutputObj<TH1D> h2D_BC_vs_CollIndex_4{TH1D("h2D_BC_vs_CollIndex_DF4", "h2D_BC_vs_CollIndex_DF4;collision index;orbit", 10001, -0.5, 10000.5)};

  void init(InitContext const&)
  {
    // cuts on phi vs pt (to avoid TPC boundaries)
    fPhiCutExpPosHigh = new TF1("fPhiCutExpPosHigh", "[0]*exp([1]*x)+[2]", 0, 20);
    fPhiCutExpPosLow = new TF1("fPhiCutExpPosLow", "[0]*exp([1]*x)+[2]", 0, 20);

    fPhiCutExpNegHigh = new TF1("fPhiCutExpNegHigh", "[0]*exp([1]*x)+[2]", 0, 20);
    fPhiCutExpNegLow = new TF1("fPhiCutExpNegLow", "[0]*exp([1]*x)+[2]", 0, 20);

    fPhiCutExpPosHigh->SetParameters(0.3, -1, 0.055 + constPhiShift);
    fPhiCutExpPosLow->SetParameters(0.15, -1, -0.02 + constPhiShift);

    fPhiCutExpNegHigh->SetParameters(-0.15, -1, +0.02 + constPhiShift);
    fPhiCutExpNegLow->SetParameters(-0.3, -1, -0.055 + constPhiShift);

    // ### event-wise:
    // AxisSpec axisNcontrib{flagPbPb ? 8001 : 501, -0.5, flagPbPb ? 8000.5 : 500.5, "n vertex contributors"};
    // AxisSpec axisNtracks{flagPbPb ? 8001 : 501, -0.5, flagPbPb ? 8000.5 : 500.5, "n tracks"};
    AxisSpec axisNcontrib{500, -0.5, flagPbPb ? 7999.5 : 499.5, "n vertex contributors"};
    AxisSpec axisNtracks{500, -0.5, flagPbPb ? 7999.5 : 499.5, "n tracks"};
    AxisSpec axisCollIndex{flagPbPb ? 501 : 10001, -0.5, flagPbPb ? 500.5 : 10000.5, "CollIndex"};
    AxisSpec axisCollTime{1000, -50, 50, "CollTime"};
    AxisSpec axisCollTimeRes{2000, -20, 20, "CollTimeRes"};
    AxisSpec axisBC{3601, -0.5, 3600.5, "bc"};

    int myMaxOrbitsPerTF = 4 * 32; // 128 - in 2022, 32 - in 2023
    AxisSpec axisBCinTF{myMaxOrbitsPerTF * 4000 + 1, -0.5, myMaxOrbitsPerTF * 4000 + 0.5, "bc"};
    // AxisSpec axisOrbitInTF{myMaxOrbitsPerTF, -2.5, 130.5, "bc"};

    // AxisSpec axisOrbit{1025, -0.5, 1024.5, "orbit"};
    AxisSpec axisOrbit{100, -0.5, 9999.5, "orbit"};

    histosEvent.add("hRunNumber", "hRunNumber", kTH1D, {{6000, 534000.5, 540000.5, "hRunNumber"}});
    histosEvent.add("hMF", "hMF", kTH1D, {{100, -1, 1, "hMF"}});

    histosEvent.add("hCollIndexBef", "hCollIndexBef", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTimeBef", "hCollTimeBef", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeResBef", "hCollTimeResBef", kTH1D, {axisCollTimeRes});
    histosEvent.add("hNumContribBef", "hNumContribBef", kTH1D, {axisNcontrib});
    histosEvent.add("hNtracksBef", "hNtracksBef", kTH1D, {axisNtracks});
    histosEvent.add("hBC_Bef", "hBC_Bef", kTH1D, {axisBC});

    // after ev sel
    histosEvent.add("hCollIndexAft", "hCollIndexAft", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTimeAft", "hCollTimeAft", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeResAft", "hCollTimeResAft", kTH1D, {axisCollTimeRes});
    histosEvent.add("hNumContribAft", "hNumContribAft", kTH1D, {axisNcontrib});
    histosEvent.add("hNtracksAft", "hNtracksAft", kTH1D, {axisNtracks});
    histosEvent.add("hBC_Aft", "hBC_Aft", kTH1D, {axisBC});
    histosEvent.add("hBCFound_Aft", "hBCFound_Aft", kTH1D, {axisBC});
    // histosEvent.add("h2D_numContrib_vs_collIndex", "h2D_numContrib_vs_collIndex", kTH2D, {axisCollIndex, axisNcontrib});
    histosEvent.add("h2D_numContrib_vs_BC", "h2D_numContrib_vs_BC", kTH2D, {axisBC, axisNcontrib});
    histosEvent.add("h2D_diffFoundBC_vs_BC", "h2D_diffFoundBC_vs_BC", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});

    histosEvent.add("hOrbitStartFromCollIndexZeroAft", "hOrbitStartFromCollIndexZeroAft", kTH1D, {axisOrbit});
    histosEvent.add("h2D_Orbit_vs_CollIndex_Aft", "h2D_Orbit_vs_CollIndex_Aft", kTH2D, {axisCollIndex, axisOrbit});

    histosEvent.add("hNtrackshGlobalAft", "hNtrackshGlobalAft", kTH1D, {axisNtracks});
    histosEvent.add("hNtrackshGlobalAft_AfterTimeFrameCut", "hNtrackshGlobalAft_AfterTimeFrameCut", kTH1D, {axisNtracks});

    // hist vs bcInTF, Feb 2, 2024
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz", "hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz_ReallyAllContrib", "hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz_ReallyAllContrib;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_BEFORE_Vz", "hNumContrib_vs_bcInTF_BEFORE_Vz;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_ReallyAllContrib", "hNumContrib_vs_bcInTF_ReallyAllContrib;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF", "hNumContrib_vs_bcInTF;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_ReallyAllContrib_AfterTimeFrameCut", "hNumContrib_vs_bcInTF_ReallyAllContrib_AfterTimeFrameCut;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_After_ITS_ROF_cut", "hNumContrib_vs_bcInTF_After_ITS_ROF_cut;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hNumContrib_vs_bcInTF_AfterTimeFrameCut", "hNumContrib_vs_bcInTF_AfterTimeFrameCut;bc in TF; n vertex contributors", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hGlobalTracks_vs_bcInTF", "hGlobalTracks_vs_bcInTF;bc in TF; n tracks", kTH1D, {axisBCinTF});

    histosEventBcInTF.add("hIsTriggerTVX_vs_bcInTF_BEFORE_SEL8_AND_Vz", "hIsTriggerTVX_vs_bcInTF_BEFORE_SEL8_AND_Vz;bc in TF;IsTriggerTVX", kTH1D, {axisBCinTF});
    histosEventBcInTF.add("hIsTriggerTVX_vs_bcInTF", "hIsTriggerTVX_vs_bcInTF;bc in TF;IsTriggerTVX", kTH1D, {axisBCinTF});

    // nTracks vs BC
    histosEvent.add("h2D_nTracksBeforeCuts_vs_BC", "h2D_nTracksBeforeCuts_vs_BC", kTH2D, {axisBC, axisNcontrib});
    histosEvent.add("h2D_nTracksAfterEtaTPCcuts_vs_BC", "h2D_nTracksAfterEtaTPCcuts_vs_BC", kTH2D, {axisBC, axisNcontrib});
    histosEvent.add("h2D_nTracksITSonly_vs_BC", "h2D_nTracksITSonly_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithITS_vs_BC", "h2D_nTracksWithITS_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithITS7hits_vs_BC", "h2D_nTracksWithITS7hits_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithITSandTPC_vs_BC", "h2D_nTracksWithITSandTPC_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithTRD_vs_BC", "h2D_nTracksWithTRD_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithTOF_vs_BC", "h2D_nTracksWithTOF_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksWithTRDorTOF_vs_BC", "h2D_nTracksWithTRDorTOF_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksGlobal_vs_BC", "h2D_nTracksGlobal_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksGlobalWithITS7hits_vs_BC", "h2D_nTracksGlobalWithITS7hits_vs_BC", kTH2D, {axisBC, axisNtracks});
    histosEvent.add("h2D_nTracksGlobalWithTRDorTOF_vs_BC", "h2D_nTracksGlobalWithTRDorTOF_vs_BC", kTH2D, {axisBC, axisNtracks});

    // AxisSpec axisGlobalTracks{flagPbPb ? 2501 : 81, -0.5, flagPbPb ? 2500.5 : 80.5, "n global tracks"};
    histosEvent.add("h1D_EventCounter_vs_BC", "h1D_EventCounter_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks0Global_vs_BC", "h1D_nTracks0Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks1Global_vs_BC", "h1D_nTracks1Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks2Global_vs_BC", "h1D_nTracks2Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks3Global_vs_BC", "h1D_nTracks3Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks4Global_vs_BC", "h1D_nTracks4Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks5Global_vs_BC", "h1D_nTracks5Global_vs_BC;bc;entries", kTH1D, {axisBC});

    // mean pT vs BC
    AxisSpec axisMeanPt{10, 0., 4.0, "<p_{T}>"};
    histosEvent.add("h1D_vContributors_1_meanPt_vs_BC", "h1D_vContributors_1_meanPt_vs_BC;bc;entries", kTH2D, {axisBC, axisMeanPt});
    histosEvent.add("h1D_vContributors_2_meanPt_vs_BC", "h1D_vContributors_2_meanPt_vs_BC;bc;entries", kTH2D, {axisBC, axisMeanPt});
    histosEvent.add("h1D_vContributors_3_meanPt_vs_BC", "h1D_vContributors_3_meanPt_vs_BC;bc;entries", kTH2D, {axisBC, axisMeanPt});
    histosEvent.add("h1D_vContributors_4_meanPt_vs_BC", "h1D_vContributors_4_meanPt_vs_BC;bc;entries", kTH2D, {axisBC, axisMeanPt});
    histosEvent.add("h1D_vContributors_5_meanPt_vs_BC", "h1D_vContributors_5_meanPt_vs_BC;bc;entries", kTH2D, {axisBC, axisMeanPt});

    //
    histosEvent.add("hNumContribAfterTPCcuts", "hNumContribAfterTPCcuts", kTH1D, {axisNcontrib});
    histosEvent.add("hNumContribITS7hits", "hNumContribITS7hits", kTH1D, {axisNcontrib});

    // only verteces with nContr>=3 with 7 ITS clusters
    histosEvent.add("hCollIndex_vertNcontr3_withITS7hits", "hCollIndex_vertNcontr3_withITS7hits", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTime_vertNcontr3_withITS7hits", "hCollTime_vertNcontr3_withITS7hits", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeRes_vertNcontr3_withITS7hits", "hCollTimeRes_vertNcontr3_withITS7hits", kTH1D, {axisCollTimeRes});
    histosEvent.add("hNumContrib_vertNcontr3_withITS7hits", "hNumContrib_vertNcontr3_withITS7hits", kTH1D, {axisNcontrib});
    histosEvent.add("hBC_vertNcontr3_withITS7hits", "hBC_vertNcontr3_withITS7hits", kTH1D, {axisBC});

    // only verteces with nContr>=3 that have TOF or TRD track
    histosEvent.add("hCollIndex_vertNcontr3_TRDorTOF", "hCollIndex_vertNcontr3_TRDorTOF", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTime_vertNcontr3_TRDorTOF", "hCollTime_vertNcontr3_TRDorTOF", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeRes_vertNcontr3_TRDorTOF", "hCollTimeRes_vertNcontr3_TRDorTOF", kTH1D, {axisCollTimeRes});
    histosEvent.add("hNumContrib_vertNcontr3_TRDorTOF", "hNumContrib_vertNcontr3_TRDorTOF", kTH1D, {axisNcontrib});
    histosEvent.add("hBC_vertNcontr3_TRDorTOF", "hBC_vertNcontr3_TRDorTOF", kTH1D, {axisBC});

    // only verteces with nContr>=3 with 7 ITS clusters + has TOF or TRD track
    histosEvent.add("hCollIndex_vertNcontr3_withITS7hits_and_TRDorTOF", "hCollIndex_vertNcontr3_withITS7hits_and_TRDorTOF", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTime_vertNcontr3_withITS7hits_and_TRDorTOF", "hCollTime_vertNcontr3_withITS7hits_and_TRDorTOF", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeRes_vertNcontr3_withITS7hits_and_TRDorTOF", "hCollTimeRes_vertNcontr3_withITS7hits_and_TRDorTOF", kTH1D, {axisCollTimeRes});
    histosEvent.add("hNumContrib_vertNcontr3_withITS7hits_and_TRDorTOF", "hNumContrib_vertNcontr3_withITS7hits_and_TRDorTOF", kTH1D, {axisNcontrib});
    histosEvent.add("hBC_vertNcontr3_withITS7hits_and_TRDorTOF", "hBC_vertNcontr3_withITS7hits_and_TRDorTOF", kTH1D, {axisBC});

    // only events with cut on ITS RO frame
    histosEvent.add("hCollIndex_ITS_ROF_cut", "hCollIndex_ITS_ROF_cut", kTH1D, {axisCollIndex});
    histosEvent.add("hCollTime_ITS_ROF_cut", "hCollTime_ITS_ROF_cut", kTH1D, {axisCollTime});
    histosEvent.add("hCollTimeRes_ITS_ROF_cut", "hCollTimeRes_ITS_ROF_cut", kTH1D, {axisCollTimeRes});
    // histosEvent.add("hNumContrib_ITS_ROF_cut", "hNumContrib_ITS_ROF_cut", kTH1D, {axisNcontrib});
    // histosEvent.add("h2D_numContrib_vs_collIndex_ITS_ROF_cut", "h2D_numContrib_vs_collIndex_ITS_ROF_cut", kTH2D, {axisCollIndex, axisNcontrib});
    histosEvent.add("hBC_ITS_ROF_cut", "hBC_ITS_ROF_cut", kTH1D, {axisBC});

    // only BC with FT0
    histosEvent.add("hBC_vertNcontr3_with_FT0", "hBC_vertNcontr3_with_FT0", kTH1D, {axisBC});
    // only BC with FT0 and diff FT0-tracks vertex < 1 cm
    histosEvent.add("hBC_vertNcontr3_with_FT0_diffPV_1cm", "hBC_vertNcontr3_with_FT0", kTH1D, {axisBC});

    // histos from Alex D.
    AxisSpec axisCent{100, 0.f, 100.f, "centrality"};
    AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};
    AxisSpec axisMultFw{200, 0, 200000, "mult Fwd"};           //{1000, 0, 200000, "mult"};
    AxisSpec axisMult{200, 0.f, 5000.f, "multiplicity"};       //{1000, 0.f, 5000.f, "multiplicity"};
    AxisSpec axisMultAllTr{200, 0.f, 20000.f, "multiplicity"}; //{1000, 0.f, 5000.f, "multiplicity"};

    histosEvent.add("vtxCutsBef", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});

    histosEvent.add("multAllTr_vs_CentBef", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CBef", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0ABef", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0ABef", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVBef", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentBef", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CBef", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0ABef", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0ABef", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVBef", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentBef", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CBef", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0ABef", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0ABef", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentBef", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0ABef", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0ABef", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CBef", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentBef", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    // histosEvent.add("nTracksITSonlyvsnTracksWithITSandTPCBef", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
    // histosEvent.add("nTracksITSonlyvsnTracksWithITSandTPCNBef", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});

    histosEvent.add("vtxCutsAft", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
    histosEvent.add("multAllTr_vs_CentAft", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    histosEvent.add("multAllTr_vs_CentAft_TFcut", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_TFcut", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_TFcut", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_TFcut", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_TFcut", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_TFcut", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_TFcut", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_TFcut", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_TFcut", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_TFcut", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_TFcut", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_TFcut", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_TFcut", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_TFcut", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_TFcut", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_TFcut", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_TFcut", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_TFcut", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_TFcut", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    histosEvent.add("multAllTr_vs_CentAft_ITSROFcut", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_ITSROFcut", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_ITSROFcut", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_ITSROFcut", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_ITSROFcut", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_ITSROFcut", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_ITSROFcut", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_ITSROFcut", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_ITSROFcut", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_ITSROFcut", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_ITSROFcut", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_ITSROFcut", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_ITSROFcut", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_ITSROFcut", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_ITSROFcut", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_ITSROFcut", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_ITSROFcut", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_ITSROFcut", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_ITSROFcut", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    histosEvent.add("multAllTr_vs_CentAft_ITSROF_TF_cuts", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_ITSROF_TF_cuts", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_ITSROF_TF_cuts", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_ITSROF_TF_cuts", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_ITSROF_TF_cuts", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_ITSROF_TF_cuts", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_ITSROF_TF_cuts", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_ITSROF_TF_cuts", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_ITSROF_TF_cuts", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_ITSROF_TF_cuts", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_ITSROF_TF_cuts", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_ITSROF_TF_cuts", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_ITSROF_TF_cuts", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_ITSROF_TF_cuts", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_ITSROF_TF_cuts", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_ITSROF_TF_cuts", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_ITSROF_TF_cuts", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_ITSROF_TF_cuts", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_ITSROF_TF_cuts", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // nTracksGlobalPVAccepted >= 2
    histosEvent.add("multAllTr_vs_CentAft_2globalPVcontrib", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_2globalPVcontrib", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_2globalPVcontrib", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_2globalPVcontrib", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_2globalPVcontrib", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_2globalPVcontrib", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_2globalPVcontrib", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_2globalPVcontrib", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_2globalPVcontrib", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_2globalPVcontrib", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_2globalPVcontrib", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_2globalPVcontrib", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_2globalPVcontrib", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_2globalPVcontrib", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_2globalPVcontrib", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_2globalPVcontrib", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_2globalPVcontrib", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_2globalPVcontrib", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // nTracksGlobalPVAccepted with 7 ITS hits >= 2
    histosEvent.add("multAllTr_vs_CentAft_2globalPVcontrib_ITS7hits", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_2globalPVcontrib_ITS7hits", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib_ITS7hits", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_2globalPVcontrib_ITS7hits", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_2globalPVcontrib_ITS7hits", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_2globalPVcontrib_ITS7hits", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_2globalPVcontrib_ITS7hits", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_2globalPVcontrib_ITS7hits", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_2globalPVcontrib_ITS7hits", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_2globalPVcontrib_ITS7hits", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_2globalPVcontrib_ITS7hits", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_2globalPVcontrib_ITS7hits", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // nTracksGlobalPVAccepted with TRD or TOF >= 2
    histosEvent.add("multAllTr_vs_CentAft_2globalPVcontrib_TRDorTOF", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_2globalPVcontrib_TRDorTOF", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib_TRDorTOF", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_2globalPVcontrib_TRDorTOF", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_2globalPVcontrib_TRDorTOF", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_2globalPVcontrib_TRDorTOF", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_2globalPVcontrib_TRDorTOF", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_2globalPVcontrib_TRDorTOF", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_2globalPVcontrib_TRDorTOF", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_2globalPVcontrib_TRDorTOF", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_2globalPVcontrib_TRDorTOF", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_2globalPVcontrib_TRDorTOF", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // diffFoundBC_vs_BC cut
    histosEvent.add("multAllTr_vs_CentAft_diffFoundBC_vs_BC_0", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_diffFoundBC_vs_BC_0", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_diffFoundBC_vs_BC_0", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_diffFoundBC_vs_BC_0", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_diffFoundBC_vs_BC_0", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_diffFoundBC_vs_BC_0", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_diffFoundBC_vs_BC_0", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_diffFoundBC_vs_BC_0", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_diffFoundBC_vs_BC_0", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_diffFoundBC_vs_BC_0", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_diffFoundBC_vs_BC_0", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_diffFoundBC_vs_BC_0", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // FT0 present
    histosEvent.add("multAllTr_vs_CentAft_hasFT0_CorrectedValid", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_hasFT0_CorrectedValid", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_hasFT0_CorrectedValid", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_hasFT0_CorrectedValid", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_hasFT0_CorrectedValid", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_hasFT0_CorrectedValid", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_hasFT0_CorrectedValid", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_hasFT0_CorrectedValid", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_hasFT0_CorrectedValid", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_hasFT0_CorrectedValid", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_hasFT0_CorrectedValid", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_hasFT0_CorrectedValid", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // cut on diff FT0-tracks vZ < 1cm
    histosEvent.add("multAllTr_vs_CentAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_1cm", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_after_PV_FT0_diff_cut_1cm", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_after_PV_FT0_diff_cut_1cm", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_after_PV_FT0_diff_cut_1cm", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_after_PV_FT0_diff_cut_1cm", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_after_PV_FT0_diff_cut_1cm", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_after_PV_FT0_diff_cut_1cm", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // cut on diff FT0-tracks vZ TIGHT
    histosEvent.add("multAllTr_vs_CentAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_after_PV_FT0_diff_cut_TIGHT", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // ALL CUTS SIMULT
    histosEvent.add("multAllTr_vs_CentAft_ALL_CUTS", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_ALL_CUTS", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_ALL_CUTS", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_ALL_CUTS", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_ALL_CUTS", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_ALL_CUTS", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_ALL_CUTS", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_ALL_CUTS", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_ALL_CUTS", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_ALL_CUTS", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_ALL_CUTS", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_ALL_CUTS", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_ALL_CUTS", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_ALL_CUTS", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_ALL_CUTS", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_ALL_CUTS", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_ALL_CUTS", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_ALL_CUTS", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_ALL_CUTS", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // ALL CUTS SIMULT TIGHTER
    histosEvent.add("multAllTr_vs_CentAft_ALL_CUTS_Tighter", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0CAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multV0AAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multT0AAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMultAllTr});
    histosEvent.add("multAllTr_vs_multTrkPVAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMultAllTr});
    histosEvent.add("multGlobalTr_vs_CentAft_ALL_CUTS_Tighter", " multiplicity vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0CAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multV0AAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multT0AAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multGlobalTr_vs_multTrkPVAft_ALL_CUTS_Tighter", " multiplicity vs multiplicity PV", kTH2F, {axisMult, axisMult});
    histosEvent.add("multTrkPV_vs_CentAft_ALL_CUTS_Tighter", " multiplicity PV vs centrality T0C", kTH2F, {axisCent, axisMult});
    histosEvent.add("multTrkPV_vs_multT0CAft_ALL_CUTS_Tighter", " multiplicity PV vs multiplicity T0C", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multV0AAft_ALL_CUTS_Tighter", " multiplicity PV vs multiplicity V0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multTrkPV_vs_multT0AAft_ALL_CUTS_Tighter", " multiplicity PV vs multiplicity T0A", kTH2F, {axisMultFw, axisMult});
    histosEvent.add("multV0A_vs_CentAft_ALL_CUTS_Tighter", " multiplicity V0A vs centrality T0C", kTH2F, {axisCent, axisMultFw});
    histosEvent.add("multT0C_vs_multT0AAft_ALL_CUTS_Tighter", " multiplicity T0C vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0AAft_ALL_CUTS_Tighter", " multiplicity V0A vs multiplicity T0A", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multV0A_vs_multT0CAft_ALL_CUTS_Tighter", " multiplicity V0A vs multiplicity T0C", kTH2F, {axisMultFw, axisMultFw});
    histosEvent.add("multT0C_vs_CentAft_ALL_CUTS_Tighter", " multiplicity T0C vs centrality T0C", kTH2F, {axisCent, axisMultFw});

    // histosEvent.add("nTracksITSonlyvsnTracksWithITSandTPCAft", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});
    // histosEvent.add("nTracksITSonlyvsnTracksWithITSandTPCNAft", " multiplicity ITS vs multiplicity ITS+TPC", kTH2F, {axisMult, axisMult});

    AxisSpec axisCounterAllVsTF{10, -0.5, 9.5, "cut"};
    histosEventCounters.add("hNtracksAll_vs_variousCuts", "hNtracksAll_vs_variousCuts", kTH1D, {axisCounterAllVsTF});
    histosEventCounters.add("hNtracksGlobal_vs_variousCuts", "hNtracksGlobal_vs_variousCuts", kTH1D, {axisCounterAllVsTF});
    histosEventCounters.add("hNtotalCollisions_vs_variousCuts", "hNtotalCollisions_vs_variousCuts", kTH1D, {axisCounterAllVsTF});

    AxisSpec axisVtxChi2{500, 0., 100.f, "chi2"};
    histosEvent.add("vtxChi2Aft", "vtxChi2Aft; chi2; Counts", kTH1F, {axisVtxChi2});
    histosEvent.add("vtxChi2_vertNcontr3_TRDorTOF", "vtxChi2_vertNcontr3_TRDorTOF; chi2; Counts", kTH1F, {axisVtxChi2});

    // #### FT0 time
    const AxisSpec axisTimeFT0{500, -5., 5., "collision time (ns)"};
    const AxisSpec axisColTimeResFT0{!flagPbPb ? 500 : 2000, -0.5, 0.5, "(T0A - T0C)/2 (ns)"};
    const AxisSpec axisVertexFT0{300, -30., 30.};
    const AxisSpec axisVertexFT0diff{1200, -30., 30.};

    histosFT0.add("hT0A", "T0A;T0A time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0C", "T0C;T0C time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0AC", "T0AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0res", "FT0 resolution", kTH1F, {axisColTimeResFT0});
    // histos.add("hColTime", "", kTH1F, {axisTimeFT0});
    // #### FT0 vertex
    histosFT0.add("hT0vertex", "FT0 vertex;FT0 vertex (cm);counts", kTH1F, {axisVertexFT0});
    histosFT0.add("hPV", "PV;primary vertex (cm);counts", kTH1F, {axisVertexFT0});

    histosFT0.add("hT0vertexDiff", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_after_EvSelAndVz", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_after_ITSROFcut", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_after_TFcut", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_after_ITSROF_and_TFcut", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_2globalPVcontrib", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_2goodPVcontribTPC", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_4globalPVcontrib", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_2globalPVcontrib_ITS7hits", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});
    histosFT0.add("hT0vertexDiff_2globalPVcontrib_TOForTRD", "FT0V - PV;FT0 vertex -  PV (cm);counts", kTH1F, {axisVertexFT0diff});

    histosFT0.add("hVertex_T0_PV", "PV vs. FT0V;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_EvSelAndVz", "PV vs. FT0V after EvSelAndVz;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_ITSROFcut", "PV vs. FT0V after ITSROF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_TFcut", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_ITSROF_and_TFcut", "PV vs. FT0V after ITSROF and TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2goodPVcontribTPC", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_4globalPVcontrib", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib_ITS7hits", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib_TOForTRD", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});

    // IA:
    histosFT0.add("hT0_sum_AC", "hT0_sum_AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});

    // ### track-wise:
    AxisSpec axisPt{nBinsPt, 0, 20, "p_{T}"};
    AxisSpec axisEta{nBinsEta, -1.5, +1.5, "#eta"};
    AxisSpec axisPhi{nBinsPhi, 0, TMath::TwoPi(), "#varphi"};
    AxisSpec axisPhiSpecMod9{nBinsPhi, -2 * TMath::TwoPi() / 9, 2 * TMath::TwoPi() / 9, "#varphi"};
    histosTracks.add("etaHistogram", "etaHistogram", kTH1D, {axisEta});
    histosTracks.add("etaHistogramAfter08cut", "etaHistogramAfter08cut", kTH1D, {axisEta});
    histosTracks.add("ptHistogram", "ptHistogram", kTH1D, {axisPt});
    histosTracks.add("phiHistogram", "phiHistogram", kTH1D, {axisPhi});

    histosTracks.add("pidCombSigma", "pidCombSigma", kTH1D, {{200, 0, 10, "pidCombSigma"}});

    // AxisSpec axisNevents{10, 0, 10, "n events"};
    // histosTracks.add("hEventCounter", "hEventCounter", kTH1D, {axisNevents});

    histosTracks.add("hTpcNClsCrossedRows", "hTpcNClsCrossedRows", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRows"}});
    histosTracks.add("hTpcNClsCrossedRowsITS7hits", "hTpcNClsCrossedRowsITS7hits", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRowsITS7hits"}});
    histosTracks.add("hNumITSclusters", "hNumITSclusters", kTH1D, {{10, -0.5, 9.5, "NumITSclusters"}});
    histosTracks.add("hDcaXY", "hDcaXY", kTH1D, {{800, -4, 4, "DcaXY"}});
    histosTracks.add("hTrackLength", "hTrackLength", kTH1D, {{400, 0, 2000, "TrackLength"}});

    AxisSpec axisTrackChi2{200, 0., 10.f, "chi2"};
    // AxisSpec axisPtBinsForPhiGaps{{0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}, "p_{T}"};
    AxisSpec axisPtBinsForPhiGaps{{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0}, "p_{T}"};
    histosTracks.add("hChi2TPCperDOF", "hChi2TPCperDOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hChi2ITSperDOF", "hChi2ITSperDOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hTrackTime", "hTrackTime", kTH1D, {axisCollTime});
    histosTracks.add("hTrackTimeRes", "hTrackTimeRes", kTH1D, {axisCollTime});

    histosTracks.add("etaHistogram_AftCuts", "etaHistogram_AftCuts", kTH1D, {axisEta});
    histosTracks.add("etaHistogramITS7hits_AftCuts", "etaHistogramITS7hits_AftCuts", kTH1D, {axisEta});
    histosTracks.add("ptHistogram_AftCuts", "ptHistogram_AftCuts", kTH1D, {axisPt});
    histosTracks.add("phiHistogram_AftCuts", "phiHistogram_AftCuts", kTH1D, {axisPhi});

    histosTracks.add("hTpcNClsCrossedRows_AftCuts", "hTpcNClsCrossedRows_AftCuts", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRows_AftCuts"}});
    histosTracks.add("hTpcNClsCrossedRowsITS7hits_AftCuts", "hTpcNClsCrossedRowsITS7hits_AftCuts", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRowsITS7hits_AftCuts"}});
    histosTracks.add("hNumITSclusters_AftCuts", "hNumITSclusters_AftCuts", kTH1D, {{10, -0.5, 9.5, "NumITSclusters_AftCuts"}});
    histosTracks.add("hDcaXY_AftCuts", "hDcaXY_AftCuts", kTH1D, {{800, -4, 4, "DcaXY_AftCuts"}});
    histosTracks.add("hTrackLength_AftCuts", "hTrackLength_AftCuts", kTH1D, {{400, 0, 2000, "TrackLength_AftCuts"}});
    histosTracks.add("hChi2perDOF_AftCuts", "hChi2perDOF_AftCuts", kTH1D, {{200, 0, 20, "Chi2perDOF_AftCuts"}});

    histosTracks.add("hTpcNClsFound_AftCuts", "hTpcNClsFound_AftCuts", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound_AftCuts"}});
    histosTracks.add("hTpcNClsFound_ITS7hits_AftCuts", "hTpcNClsFound_ITS7hits_AftCuts", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound_ITS7hits_AftCuts"}});

    // ### chi2 after cuts:
    histosTracks.add("hChi2TPCperDOF_ITS7hits", "hChi2TPCperDOF_ITS7hits", kTH1D, {axisTrackChi2});
    histosTracks.add("hChi2ITSperDOF_ITS7hits", "hChi2ITSperDOF_ITS7hits", kTH1D, {axisTrackChi2});
    histosTracks.add("hTrackTime_ITS7hits", "hTrackTime", kTH1D, {axisCollTime});
    histosTracks.add("hTrackTimeRes_ITS7hits", "hTrackTimeRes_ITS7hits", kTH1D, {axisCollTime});

    histosTracks.add("hChi2TPCperDOF_TRDorTOF", "hChi2TPCperDOF_TRDorTOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hChi2ITSperDOF_TRDorTOF", "hChi2ITSperDOF_TRDorTOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hTrackTime_TRDorTOF", "hTrackTimeTRDorTOF", kTH1D, {axisCollTime});
    histosTracks.add("hTrackTimeResTRDorTOF", "hTrackTimeResTRDorTOF", kTH1D, {axisCollTime});

    // QA TPC sector boundaries
    histosTracks.add("posSelTrack_pt", "posSelTrack_pt", kTH1D, {axisPt});
    histosTracks.add("posSelTrack_phi", "posSelTrack_phi", kTH1D, {axisPhi});
    histosTracks.add("posSelTrack_eta", "posSelTrack_eta", kTH1D, {axisEta});
    histosTracks.add("posSelTrack_phi_pT_05_10", "posSelTrack_phi_pT_05_10", kTH1D, {axisPhi});
    histosTracks.add("posSelTrack_phi_pT_10_20", "posSelTrack_phi_pT_10_20", kTH1D, {axisPhi});
    histosTracks.add("posSelTrack_phi_pT_20_100", "posSelTrack_phi_pT_20_100", kTH1D, {axisPhi});
    histosTracks.add("posSelTrack_phi_vs_pt", "posSelTrack_phi_vs_pt", kTH2D, {axisPtBinsForPhiGaps, axisPhi});
    histosTracks.add("posSelTrack_phi_vs_pt_modPiOver9", "posSelTrack_phi_vs_pt_modPiOver9", kTH2D, {axisPtBinsForPhiGaps, axisPhiSpecMod9});
    histosTracks.add("posSelTrack_phi_vs_pt_afterCut", "posSelTrack_phi_vs_pt_afterCut", kTH2D, {axisPtBinsForPhiGaps, axisPhi});

    histosTracks.add("negSelTrack_pt", "negSelTrack_pt", kTH1D, {axisPt});
    histosTracks.add("negSelTrack_phi", "negSelTrack_phi", kTH1D, {axisPhi});
    histosTracks.add("negSelTrack_eta", "negSelTrack_eta", kTH1D, {axisEta});
    histosTracks.add("negSelTrack_phi_pT_05_10", "negSelTrack_phi_pT_05_10", kTH1D, {axisPhi});
    histosTracks.add("negSelTrack_phi_pT_10_20", "negSelTrack_phi_pT_10_20", kTH1D, {axisPhi});
    histosTracks.add("negSelTrack_phi_pT_20_100", "negSelTrack_phi_pT_20_100", kTH1D, {axisPhi});
    histosTracks.add("negSelTrack_phi_vs_pt", "negSelTrack_phi_vs_pt", kTH2D, {axisPtBinsForPhiGaps, axisPhi});
    histosTracks.add("negSelTrack_phi_vs_pt_modPiOver9", "negSelTrack_phi_vs_pt_modPiOver9", kTH2D, {axisPtBinsForPhiGaps, axisPhiSpecMod9});
    histosTracks.add("negSelTrack_phi_vs_pt_afterCut", "negSelTrack_phi_vs_pt_afterCut", kTH2D, {axisPtBinsForPhiGaps, axisPhi});

    // Configurable<float> vtxZ{"vtxZ", 10.f, ""};
    // Filter posZfilter = nabs(aod::collision::posZ) < vtxZ;

    // ##### v0s
    // K0s reconstruction
    // AxisSpec vertexZAxis = {400, -15., 15., "vrtx_{Z} [cm]"};
    AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    histosK0S.add("hK0Sradius", "hK0Sradius", kTH1D, {{800, 0, 100, "hK0Sradius"}});
    histosK0S.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
    histosK0S.add("hMassK0ShortAfterSelection", "hMassK0ShortAfterSelection", {HistType::kTH1F, {K0ShortMassAxis}});
    // for(int iPt=0; iPt<nPtRangesK0S; iPt++)
    // {
    // string hName = "hMassK0Short_pT"+to_string(iPt);
    // arrHistPtK0Snames[iPt] = hName;
    AxisSpec axisK0SptBins{{0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}, "p_{T}"};
    histosK0S.add("hMassK0ShortAfterSelectionVsPt", "hMassK0ShortAfterSelectionVsPt", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtEta01_08", "hMassK0ShortAfterSelectionVsPtEta01_08", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtEta08_01", "hMassK0ShortAfterSelectionVsPtEta08_01", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigma", "hMassK0ShortAfterSelectionVsPtPIDnSigma", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4_andROFcut", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4_andROFcut", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNITS7hitsContribAbove2", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNITS7hitsContribAbove2", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_nVertexContributorsWithTRDorTOF_Above2", "hMassK0S_nVertexContributorsWithTRDorTOF_Above2", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_nVertexContributorsWithTRDorTOF_ITS7hits_Above2", "hMassK0S_nVertexContributorsWithTRDorTOF_ITS7hits_Above2", kTH2D, {axisK0SptBins, K0ShortMassAxis});

    histosK0S.add("hMassK0S_cutsOnDaughters", "hMassK0S_cutsOnDaughters", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_bothTOF", "hMassK0S_cutsOnDaughters_bothTOF", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_bothTRD", "hMassK0S_cutsOnDaughters_bothTRD", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_TOForTRD", "hMassK0S_cutsOnDaughters_TOForTRD", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_norTOFnorTRD", "hMassK0S_cutsOnDaughters_norTOFnorTRD", kTH2D, {axisK0SptBins, K0ShortMassAxis});

    histosK0S.add("hMassK0S_cutsOnDaughters_ITSmin4hits", "hMassK0S_cutsOnDaughters_ITSmin4hits", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_ITSmin5hits", "hMassK0S_cutsOnDaughters_ITSmin5hits", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_ITSmin6hits", "hMassK0S_cutsOnDaughters_ITSmin6hits", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_chi2tpc2", "hMassK0S_cutsOnDaughters_chi2tpc2", kTH2D, {axisK0SptBins, K0ShortMassAxis});

    histosK0S.add("hMassK0S_cutsOnDaughters_v0Rad_0_2", "hMassK0S_cutsOnDaughters_v0Rad_0_2", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_v0Rad_2_4", "hMassK0S_cutsOnDaughters_v0Rad_2_4", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_v0Rad_4_6", "hMassK0S_cutsOnDaughters_v0Rad_4_6", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_v0Rad_6_10", "hMassK0S_cutsOnDaughters_v0Rad_6_10", kTH2D, {axisK0SptBins, K0ShortMassAxis});
    histosK0S.add("hMassK0S_cutsOnDaughters_v0Rad_10_20", "hMassK0S_cutsOnDaughters_v0Rad_10_20", kTH2D, {axisK0SptBins, K0ShortMassAxis});

    histosK0S.add("hMassK0S_cut_TPC_boundaries_for_Dau", "hMassK0S_cut_TPC_boundaries_for_Dau", kTH2D, {axisK0SptBins, K0ShortMassAxis});

    histosK0S.add("posDau_pt", "posDau_pt", kTH1D, {axisPt});
    histosK0S.add("posDau_phi", "posDau_phi", kTH1D, {axisPhi});
    histosK0S.add("posDau_eta", "posDau_eta", kTH1D, {axisEta});
    histosK0S.add("posDau_phi_pT_05_10", "posDau_phi_pT_05_10", kTH1D, {axisPhi});
    histosK0S.add("posDau_phi_pT_10_20", "posDau_phi_pT_10_20", kTH1D, {axisPhi});
    histosK0S.add("posDau_phi_pT_20_100", "posDau_phi_pT_20_100", kTH1D, {axisPhi});
    histosK0S.add("posDau_phi_vs_pt", "posDau_phi_vs_pt", kTH2D, {axisPtBinsForPhiGaps, axisPhi});
    histosK0S.add("posDau_phi_vs_pt_afterCut", "posDau_phi_vs_pt_afterCut", kTH2D, {axisPtBinsForPhiGaps, axisPhi});

    histosK0S.add("negDau_pt", "negDau_pt", kTH1D, {axisPt});
    histosK0S.add("negDau_phi", "negDau_phi", kTH1D, {axisPhi});
    histosK0S.add("negDau_eta", "negDau_eta", kTH1D, {axisEta});
    histosK0S.add("negDau_phi_pT_05_10", "negDau_phi_pT_05_10", kTH1D, {axisPhi});
    histosK0S.add("negDau_phi_pT_10_20", "negDau_phi_pT_10_20", kTH1D, {axisPhi});
    histosK0S.add("negDau_phi_pT_20_100", "negDau_phi_pT_20_100", kTH1D, {axisPhi});
    histosK0S.add("negDau_phi_vs_pt", "negDau_phi_vs_pt", kTH2D, {axisPtBinsForPhiGaps, axisPhi});
    histosK0S.add("negDau_phi_vs_pt_afterCut", "negDau_phi_vs_pt_afterCut", kTH2D, {axisPtBinsForPhiGaps, axisPhi});
  }

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  // Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  // Filter<aod::Tracks> ptFilter = aod::track::pt > 1;

  // void process(aod::TracksIU const& tracks)
  // void process(aod::Collision const& collision, aod::Tracks &tracks) {
  // void process(aod::Collision const& collision, soa::Filtered<aod::Tracks> &tracks) {
  // void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks)
  // void process(aod::Collision const& collision,   aod::fullTrack const& fullTrack)  // one can use the abbreviation FullTrack, which is a predefined join of all three track tables
  // A list of predefined joins is available in The Data Model section of these documentation pages).

  // The grouping works with any number of children. In the below example the process function is given three arguments. In this case process is run for each collision with the tracks and V0s belonging to the actual collision.
  // void process(aod::Collision const& collision, aod::Tracks const& tracks, aod::V0s const& v0s)

  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi>;

  // using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Cs>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::FT0sCorrected>;
  // using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::Origins>;

  // BCsWithTimestamps - is soa::Join<aod::BCs, aod::Timestamps>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>; //, aod::Run3MatchedToBCSparse>;

  void process(
    // aod::Collision const& collision
    // soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision
    Colls::iterator const& collision,
    aod::FT0s const& ft0s,
    // aod::BCs const& bcs,
    BCsRun3 const& bcs,
    aod::Origins const& origins,
    // Colls::iterator const& collision, aod::BCsWithTimestamps const& bcs,
    // , soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA> const& tracks
    // , soa::Join<aod::Tracks, aod::pidTOFEl, aod::pidTPCEl>::iterator const& track
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
              aod::TrackSelectionExtension,
              aod::TracksDCA, aod::pidTOFEl, aod::pidTPCEl> const& tracks,
    aod::V0Datas const& v0s, DaughterTracks const&)
  // aod::V0s_001 const& v0s)
  {
    auto bc = collision.bc_as<BCsRun3>();

    uint32_t orbit = bc.globalBC() / o2::constants::lhc::LHCMaxBunches;
    if (collision.index() == 0)
      orbitAtCollIndexZero = orbit;

    // cout << "orbit = " << orbit << endl;

    histosEvent.fill(HIST("hCollIndexBef"), collision.index());
    histosEvent.fill(HIST("hCollTimeBef"), collision.collisionTime());
    histosEvent.fill(HIST("hCollTimeResBef"), collision.collisionTimeRes());
    histosEvent.fill(HIST("hNumContribBef"), collision.numContrib());
    // histosEvent.fill(HIST("hVertexZ"), collision.posZ() );
    histosEvent.fill(HIST("hNtracksBef"), tracks.size());

    auto collBC = bc.globalBC() % 3564;
    uint64_t globalFoundBC = -1;
    // int64_t globalFoundBC = -1;
    if (collision.has_foundBC()) {
      auto bcFound = collision.foundBC_as<BCsRun3>(); // collision.foundBC();
      globalFoundBC = bcFound.globalBC() % 3564;
    }
    histosEvent.fill(HIST("hBC_Bef"), collBC);

    // #### code from Alex Dobrin:
    auto t0cCentr = 0; // flagPbPb ? collision.centFT0C() : 0;

    auto multV0A = collision.multFV0A();
    auto multT0A = collision.multFT0A();
    auto multT0C = collision.multFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    Int_t multTrk = tracks.size();
    double vZ = collision.posZ();

    // Find BC associated with collision
    // uint64_t mostProbableBC = -1000;
    // if (collision.has_foundBC()) {
    // foundBCId is stored in EvSels
    // auto bc = collision.template foundBC_as<TBCs>();
    // Obtain slice of compatible BCs
    // mostProbableBC = bc.globalBC();
    // }

    // ######## DF QA
    int64_t myDF_ID = -1;
    uint64_t DF_ID_raw = -1;
    if (mapDF.size() < 10) {
      for (auto const& origin : origins) {
        uint64_t DF_ID = origin.dataframeID();
        DF_ID_raw = DF_ID;
        // look only at the id of the first subDF in this DF
        if (origin.globalIndex() == 0) {
          if (mapDF.find(DF_ID) == mapDF.end()) {
            // not found
            mapDF.insert({DF_ID, mapDF.size()});
          } else {
            // found
          }
          myDF_ID = mapDF[DF_ID];
        }
        // cout << "DF globalIndex = " << origin.globalIndex() << ", ID = " << origin.dataframeID() << ", myDF_ID = " << myDF_ID << ", mapDF.size() = " << mapDF.size() << endl;
      }
    }

    if (myDF_ID >= 0 && myDF_ID < 5) {
      int diffOrbits = (int32_t)orbit - (int32_t)orbitAtCollIndexZero;
      TString strDF = Form("DF_%d", static_cast<int>(DF_ID_raw));
      if (myDF_ID == 0) {
        h2D_Orbit_vs_CollIndex_0->Fill(collision.index(), diffOrbits);
        h2D_Orbit_vs_CollIndex_0->SetTitle(strDF);
        h2D_BC_vs_CollIndex_0->Fill(collision.index(), collBC);
        h2D_BC_vs_CollIndex_0->SetTitle(strDF);
      } else if (myDF_ID == 1) {
        h2D_Orbit_vs_CollIndex_1->Fill(collision.index(), diffOrbits);
        h2D_Orbit_vs_CollIndex_1->SetTitle(strDF);
        h2D_BC_vs_CollIndex_1->Fill(collision.index(), collBC);
        h2D_BC_vs_CollIndex_1->SetTitle(strDF);
      } else if (myDF_ID == 2) {
        h2D_Orbit_vs_CollIndex_2->Fill(collision.index(), diffOrbits);
        h2D_Orbit_vs_CollIndex_2->SetTitle(strDF);
        h2D_BC_vs_CollIndex_2->Fill(collision.index(), collBC);
        h2D_BC_vs_CollIndex_2->SetTitle(strDF);
      } else if (myDF_ID == 3) {
        h2D_Orbit_vs_CollIndex_3->Fill(collision.index(), diffOrbits);
        h2D_Orbit_vs_CollIndex_3->SetTitle(strDF);
        h2D_BC_vs_CollIndex_3->Fill(collision.index(), collBC);
        h2D_BC_vs_CollIndex_3->SetTitle(strDF);
      } else if (myDF_ID == 4) {
        h2D_Orbit_vs_CollIndex_4->Fill(collision.index(), diffOrbits);
        h2D_Orbit_vs_CollIndex_4->SetTitle(strDF);
        h2D_BC_vs_CollIndex_4->Fill(collision.index(), collBC);
        h2D_BC_vs_CollIndex_4->SetTitle(strDF);
      }
    }

    if (0)
      LOGF(info, "collision.globalIndex() = %d, index() = %d, has_foundBC() = %d, tracks.size() = %d", //, globalBC=%.1f",
           collision.globalIndex(), collision.index(),
           //  DF_ID, myDF_ID,
           (int)collision.has_foundBC(), // mostProbableBC,
           tracks.size());               //, collision.globalBC() );

    int nTracksAll = 0;
    int nTracksAfterEtaTPCCuts = 0;
    int nTracksITSonly = 0, nTracksWithITS = 0, nTracksWithITSandTPC = 0;
    int nTracksWithTRD = 0, nTracksWithTOF = 0, nTracksWithTRDorTOF = 0;
    int nTracksWithITS7hits = 0;
    int nTracksGlobalAccepted = 0;
    int nTracksGlobalWithITS7hits = 0;
    int nTracksGlobalWithTRDorTOF = 0;

    int nTracksGlobalPVAccepted = 0;
    int nTracksGlobalPVwithITS7hits = 0;
    int nTracksGlobalPVwithTRDorTOF = 0;

    int counterPVcontributorsAfterTPCcuts = 0;
    int counterPVcontributorsITS7hits = 0;

    int counterVertexContributorsWithTRDorTOF = 0;
    int counterVertexContributorsWithTRDorTOF_ITS7hits = 0;

    int counterPVcontributorsNoTOFandTRD = 0;
    double meanPtForPVContributorsNoTOFandTRD = 0;

    // ### track pre-loop
    for (auto& track : tracks) {
      nTracksAll++;

      if (fabs(track.eta()) > 0.8)
        continue;

      if (track.tpcNClsFound() < 80)
        continue;
      if (track.tpcNClsCrossedRows() < 100)
        continue;

      nTracksAfterEtaTPCCuts++;

      if (track.isPVContributor()) {
        counterPVcontributorsAfterTPCcuts++;

        if (track.itsNCls() == 7)
          counterPVcontributorsITS7hits++;

        if (track.hasTRD() || track.hasTOF()) {
          counterVertexContributorsWithTRDorTOF++;
          if (track.itsNCls() == 7)
            counterVertexContributorsWithTRDorTOF_ITS7hits++;
        }

        if (!track.hasTRD() && !track.hasTOF()) {
          counterPVcontributorsNoTOFandTRD++;
          meanPtForPVContributorsNoTOFandTRD += track.pt();
        }

        if (track.isGlobalTrack()) {
          nTracksGlobalPVAccepted++;
          if (track.itsNCls() == 7)
            nTracksGlobalPVwithITS7hits++;
          if (track.hasTRD() || track.hasTOF())
            nTracksGlobalPVwithTRDorTOF++;
        }
      }

      if (track.hasITS() && !track.hasTPC()) // Flag to check if track has a TPC match
        nTracksITSonly++;

      if (track.hasITS())
        nTracksWithITS++;
      if (track.itsNCls() == 7)
        nTracksWithITS7hits++;

      if (track.hasITS() && track.hasTPC())
        nTracksWithITSandTPC++;

      if (track.hasTRD())
        nTracksWithTRD++;
      if (track.hasTOF())
        nTracksWithTOF++;
      if (track.hasTRD() || track.hasTOF())
        nTracksWithTRDorTOF++;

      if (track.isGlobalTrack()) {
        nTracksGlobalAccepted++;
        if (track.itsNCls() == 7)
          nTracksGlobalWithITS7hits++;
        if (track.hasTRD() || track.hasTOF())
          nTracksGlobalWithTRDorTOF++;
      }
    }

    histosEvent.fill(HIST("vtxCutsBef"), vZ);
    histosEvent.fill(HIST("multAllTr_vs_multT0CBef"), multT0C, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multV0ABef"), multV0A, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multT0ABef"), multT0A, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multTrkPVBef"), multNTracksPV, multTrk);
    histosEvent.fill(HIST("multGlobalTr_vs_multT0CBef"), multT0C, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multV0ABef"), multV0A, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multT0ABef"), multT0A, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVBef"), multNTracksPV, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multTrkPV_vs_multT0CBef"), multT0C, multNTracksPV);
    histosEvent.fill(HIST("multTrkPV_vs_multV0ABef"), multV0A, multNTracksPV);
    histosEvent.fill(HIST("multTrkPV_vs_multT0ABef"), multT0A, multNTracksPV);
    histosEvent.fill(HIST("multT0C_vs_multT0ABef"), multT0A, multT0C);
    histosEvent.fill(HIST("multV0A_vs_multT0ABef"), multT0A, multV0A);
    histosEvent.fill(HIST("multV0A_vs_multT0CBef"), multT0C, multV0A);
    if (flagPbPb) {
      histosEvent.fill(HIST("multAllTr_vs_CentBef"), t0cCentr, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_CentBef"), t0cCentr, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_CentBef"), t0cCentr, multNTracksPV);
      histosEvent.fill(HIST("multV0A_vs_CentBef"), t0cCentr, multV0A);
      histosEvent.fill(HIST("multT0C_vs_CentBef"), t0cCentr, multT0C);
    }

    //
    int runNumber = bc.runNumber();
    histosEvent.fill(HIST("hRunNumber"), runNumber);

    // const auto timestamp = collision.bc_as<aod::BCsWithTimestamps>().timestamp(); /// NB: in ms

    double magneticField = 0;

    // #### begin of the code from qaPrimVtxVsTime.cxx

    if (lastRunNumber != runNumber) {
      /// execute the code in this scope only once, i.e. when the current run is considered for the first time in this DF
      lastRunNumber = runNumber;
      int64_t tsSOR = 0;
      // int64_t tsEOR = 0;

      /// reject AO2Ds for which no CCDB access is possible
      if (runNumber < 500000) {
        LOG(warning) << ">>> run number " << runNumber << " < 500000. access to CCDB not possible. Exiting";
        return;
      }

      /// If we are here, the current run was never considered before.
      /// Let's add the TH2 that we need for the monitoring
      /// Let's define the x-axis according to the start-of-run (SOR) and end-of-run (EOR) times
      // o2::ccdb::CcdbApi ccdb_api;
      // ccdb_api.init(ccdburl);
      // std::map<string, string> metadataRCT, headers;
      // headers = ccdb_api.retrieveHeaders(Form("RCT/Info/RunInformation/%i", runNumber), metadataRCT, -1);
      // tsSOR = atol(headers["SOR"].c_str());
      // tsEOR = atol(headers["EOR"].c_str());
      // double minSec = floor(tsSOR / 1000.); /// round tsSOR to the highest integer lower than tsSOR
      // double maxSec = ceil(tsEOR / 1000.);  /// round tsEOR to the lowest integer higher than tsEOR
      // const AxisSpec axisSeconds{static_cast<int>(maxSec - minSec), minSec, maxSec, "seconds (from January 1st, 1970 at UTC)"};
      // histosEvent.add("hPosZvsTime", "", kTH2F, {axisSeconds, axisZvert});

      // ##### code from Evgeny, Feb 1, 2024
      int64_t ts = bcs.iteratorAt(0).timestamp();
      // access orbit-reset timestamp
      auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
      int64_t tsOrbitReset = (*ctpx)[0]; // us

      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", runNumber);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
      uint32_t nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF
      tsSOR = grpecs->getTimeStart();                 // ms
      // tsEOR = grpecs->getTimeEnd();                   // ms

      // calculate SOR and EOR orbits
      // int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      // int64_t orbitEOR = (tsEOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;

      // adjust to the nearest TF edge
      // orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF - 1;
      // orbitEOR = orbitEOR / nOrbitsPerTF * nOrbitsPerTF - 1;

      // first bc of the first orbit (should coincide with TF start)
      // bcSOR = orbitSOR * nBCsPerOrbit;

      // duration of TF in bcs
      nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;

      // IA - try without "-1"
      int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF; // - 1;
      bcSOR = orbitSOR * nBCsPerOrbit;

      magneticField = 1;
      static o2::parameters::GRPMagField* grpo = nullptr;
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", ts);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", ts);
        return;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with L3 ", ts, grpo->getL3Current());
      // taken from GRP onject definition of getNominalL3Field; update later to something smarter (mNominalL3Field = std::lround(5.f * mL3Current / 30000.f);)
      auto NominalL3Field = std::lround(5.f * grpo->getL3Current() / 30000.f);
      magneticField = 0.1 * (NominalL3Field);
      LOGF(info, "magneticField =  %f", magneticField);

      // o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, ts); // collision.bc().timestamp());
      // if (grpo != nullptr)
      // magneticField = grpo->getNominalL3Field();
      // magneticField = collision.magField();

      // const char* ccdbpath_grp = "GLO/GRP/GRP";
      // o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
      // if (grpo != nullptr) {
      //   magneticField = grpo->getNominalL3Field();
      //   LOGF(info, "Setting magnetic field to %f kG for run %d", magneticField, runNumber);
      // } else {
      //   LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", runNumber, bc.timestamp());
      // }
    }

    // bc in Time Frame:
    // int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
    // int64_t orbitInTF = bcInTF / nBCsPerOrbit;
    int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;

    histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz"), bcInTF, counterPVcontributorsAfterTPCcuts);
    histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz_ReallyAllContrib"), bcInTF, collision.numContrib());
    histosEventBcInTF.fill(HIST("hIsTriggerTVX_vs_bcInTF_BEFORE_SEL8_AND_Vz"), bcInTF, collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ? 1 : 0);

    // is FT0
    bool isFT0 = false;
    double ft0_posZ = 0;
    double diff_PV_ft0_tracks = 0;
    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      if (collision.t0CCorrectedValid() && collision.t0ACorrectedValid())
        isFT0 = true;

      if (isFT0) {
        ft0_posZ = ft0.posZ();
        diff_PV_ft0_tracks = ft0_posZ - collision.posZ();
        histosFT0.fill(HIST("hT0AC"), collision.t0AC());
        histosFT0.fill(HIST("hT0vertex"), ft0.posZ());
        histosFT0.fill(HIST("hVertex_T0_PV"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hPV"), collision.posZ());
        histosFT0.fill(HIST("hT0res"), collision.t0resolution());
        histosFT0.fill(HIST("hT0vertexDiff"), diff_PV_ft0_tracks);

        if (collision.t0ACorrectedValid() && collision.t0CCorrectedValid()) {
          histosFT0.fill(HIST("hT0A"), collision.t0ACorrected());
          histosFT0.fill(HIST("hT0C"), collision.t0CCorrected());
          histosFT0.fill(HIST("hT0_sum_AC"), collision.t0ACorrected() + collision.t0CCorrected());
        }
      }
    } // end of if (collision.has_foundFT0())

    histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 0, tracks.size());            // ALL tracks
    histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 0, nTracksGlobalAccepted); // global tracks
    histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 0, 1);                  // collisions counter

    // do counting also with TF border cut BEFORE ev. sel.
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 1, tracks.size());            // ALL tracks
      histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 1, nTracksGlobalAccepted); // global tracks
      histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 1, 1);                  // collisions counter
    }

    // ### event selection cuts
    if (!collision.sel8())
      return;

    histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_Vz"), bcInTF, counterPVcontributorsAfterTPCcuts);

    if (TMath::Abs(vZ) > vtxCut)
      return;

    if (isFT0) {
      histosFT0.fill(HIST("hVertex_T0_PV_after_EvSelAndVz"), ft0_posZ, collision.posZ());
      histosFT0.fill(HIST("hT0vertexDiff_after_EvSelAndVz"), diff_PV_ft0_tracks);

      if (counterPVcontributorsAfterTPCcuts >= 2) {
        histosFT0.fill(HIST("hVertex_T0_PV_2goodPVcontribTPC"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_2goodPVcontribTPC"), diff_PV_ft0_tracks);
      }
      if (nTracksGlobalPVAccepted >= 2) {
        histosFT0.fill(HIST("hVertex_T0_PV_2globalPVcontrib"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_2globalPVcontrib"), diff_PV_ft0_tracks);
      }

      if (nTracksGlobalPVAccepted >= 4) {
        histosFT0.fill(HIST("hVertex_T0_PV_4globalPVcontrib"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_4globalPVcontrib"), diff_PV_ft0_tracks);
      }
      if (nTracksGlobalPVwithITS7hits >= 2) {
        histosFT0.fill(HIST("hVertex_T0_PV_2globalPVcontrib_ITS7hits"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_2globalPVcontrib_ITS7hits"), diff_PV_ft0_tracks);
      }
      if (nTracksGlobalPVwithTRDorTOF >= 2) {
        histosFT0.fill(HIST("hVertex_T0_PV_2globalPVcontrib_TOForTRD"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_2globalPVcontrib_TOForTRD"), diff_PV_ft0_tracks);
      }
    }

    histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 2, tracks.size());            // ALL tracks
    histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 2, nTracksGlobalAccepted); // global tracks
    histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 2, 1);                  // collisions counter

    histosEvent.fill(HIST("hNumContribAfterTPCcuts"), counterPVcontributorsAfterTPCcuts);
    histosEvent.fill(HIST("hNumContribITS7hits"), counterPVcontributorsITS7hits);

    histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF"), bcInTF, counterPVcontributorsAfterTPCcuts);
    histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_ReallyAllContrib"), bcInTF, collision.numContrib());
    histosEventBcInTF.fill(HIST("hIsTriggerTVX_vs_bcInTF"), bcInTF, collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ? 1 : 0);

    histosEventBcInTF.fill(HIST("hGlobalTracks_vs_bcInTF"), bcInTF, nTracksGlobalAccepted);
    histosEvent.fill(HIST("hOrbitStartFromCollIndexZeroAft"), (int32_t)orbit - (int32_t)orbitAtCollIndexZero);
    histosEvent.fill(HIST("h2D_Orbit_vs_CollIndex_Aft"), collision.index(), (int32_t)orbit - (int32_t)orbitAtCollIndexZero);

    histosEvent.fill(HIST("hMF"), magneticField);
    int MFsign = magneticField > 0 ? +1 : -1;

    histosEvent.fill(HIST("vtxChi2Aft"), collision.chi2());

    // fill 2D histograms
    histosEvent.fill(HIST("vtxCutsAft"), vZ);
    histosEvent.fill(HIST("multAllTr_vs_multT0CAft"), multT0C, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multV0AAft"), multV0A, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multT0AAft"), multT0A, multTrk);
    histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft"), multNTracksPV, multTrk);
    histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft"), multT0C, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft"), multV0A, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft"), multT0A, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft"), multNTracksPV, nTracksGlobalAccepted);
    histosEvent.fill(HIST("multTrkPV_vs_multT0CAft"), multT0C, multNTracksPV);
    histosEvent.fill(HIST("multTrkPV_vs_multV0AAft"), multV0A, multNTracksPV);
    histosEvent.fill(HIST("multTrkPV_vs_multT0AAft"), multT0A, multNTracksPV);
    histosEvent.fill(HIST("multT0C_vs_multT0AAft"), multT0A, multT0C);
    histosEvent.fill(HIST("multV0A_vs_multT0AAft"), multT0A, multV0A);
    histosEvent.fill(HIST("multV0A_vs_multT0CAft"), multT0C, multV0A);
    if (flagPbPb) {
      histosEvent.fill(HIST("multAllTr_vs_CentAft"), t0cCentr, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_CentAft"), t0cCentr, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_CentAft"), t0cCentr, multNTracksPV);
      histosEvent.fill(HIST("multV0A_vs_CentAft"), t0cCentr, multV0A);
      histosEvent.fill(HIST("multT0C_vs_CentAft"), t0cCentr, multT0C);
    }

    histosEvent.fill(HIST("hNtrackshGlobalAft"), nTracksGlobalAccepted);

    // cut TF borders
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {

      histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 3, tracks.size());            // ALL tracks with TF cut
      histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 3, nTracksGlobalAccepted); // global tracks
      histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 3, 1);                  // collisions counter

      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_ReallyAllContrib_AfterTimeFrameCut"), bcInTF, collision.numContrib());
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_AfterTimeFrameCut"), bcInTF, counterPVcontributorsAfterTPCcuts);

      histosEvent.fill(HIST("hNtrackshGlobalAft_AfterTimeFrameCut"), nTracksGlobalAccepted);

      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_TFcut"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_TFcut"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_TFcut"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_TFcut"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_TFcut"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_TFcut"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_TFcut"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_TFcut"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_TFcut"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_TFcut"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_TFcut"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_TFcut"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_TFcut"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_TFcut"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_TFcut"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_TFcut"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_TFcut"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_TFcut"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_TFcut"), t0cCentr, multT0C);
      }

      if (isFT0) {
        histosFT0.fill(HIST("hVertex_T0_PV_after_TFcut"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_after_TFcut"), diff_PV_ft0_tracks);
      }
    }

    // histosEvent.fill(HIST("nTracksITSonlyvsnTracksWithITSandTPCAft"), nTracksWithITSandTPC, nTracksITSonly);
    // histosEvent.fill(HIST("nTracksITSonlyvsnTracksWithITSandTPCNAft"), nTracksWithITSandTPC, nTracksWithITS);

    histosEvent.fill(HIST("hCollIndexAft"), collision.index());
    histosEvent.fill(HIST("hCollTimeAft"), collision.collisionTime());
    histosEvent.fill(HIST("hCollTimeResAft"), collision.collisionTimeRes());
    histosEvent.fill(HIST("hNumContribAft"), collision.numContrib());
    // histosEvent.fill(HIST("hVertexZ"), collision.posZ() );
    histosEvent.fill(HIST("hNtracksAft"), tracks.size());
    // histosEvent.fill(HIST("h2D_numContrib_vs_collIndex"), collision.index(), collision.numContrib());
    histosEvent.fill(HIST("hBC_Aft"), collBC);
    histosEvent.fill(HIST("hBCFound_Aft"), globalFoundBC);
    histosEvent.fill(HIST("h2D_numContrib_vs_BC"), collBC, collision.numContrib());

    int64_t diffFoundBC_vs_BC = (int64_t)globalFoundBC - (int64_t)collBC;
    histosEvent.fill(HIST("h2D_diffFoundBC_vs_BC"), collBC, (int64_t)globalFoundBC - (int64_t)collBC);

    // with FT0 conditions
    if (isFT0) {
      histosEvent.fill(HIST("hBC_vertNcontr3_with_FT0"), collBC);
      if (fabs(diff_PV_ft0_tracks) < 1)
        histosEvent.fill(HIST("hBC_vertNcontr3_with_FT0_diffPV_1cm"), collBC);
    }

    // nTracks vs BC
    histosEvent.fill(HIST("h2D_nTracksBeforeCuts_vs_BC"), collBC, nTracksAll);
    histosEvent.fill(HIST("h2D_nTracksAfterEtaTPCcuts_vs_BC"), collBC, nTracksAfterEtaTPCCuts);
    histosEvent.fill(HIST("h2D_nTracksITSonly_vs_BC"), collBC, nTracksITSonly);
    histosEvent.fill(HIST("h2D_nTracksWithITS_vs_BC"), collBC, nTracksWithITS);
    histosEvent.fill(HIST("h2D_nTracksWithITS7hits_vs_BC"), collBC, nTracksWithITS7hits);
    histosEvent.fill(HIST("h2D_nTracksWithITSandTPC_vs_BC"), collBC, nTracksWithITSandTPC);
    histosEvent.fill(HIST("h2D_nTracksWithTRD_vs_BC"), collBC, nTracksWithTRD);
    histosEvent.fill(HIST("h2D_nTracksWithTOF_vs_BC"), collBC, nTracksWithTOF);
    histosEvent.fill(HIST("h2D_nTracksWithTRDorTOF_vs_BC"), collBC, nTracksWithTRDorTOF);
    histosEvent.fill(HIST("h2D_nTracksGlobal_vs_BC"), collBC, nTracksGlobalAccepted);

    histosEvent.fill(HIST("h2D_nTracksGlobalWithITS7hits_vs_BC"), collBC, nTracksGlobalWithITS7hits);
    histosEvent.fill(HIST("h2D_nTracksGlobalWithTRDorTOF_vs_BC"), collBC, nTracksGlobalWithTRDorTOF);

    // QA by hand
    histosEvent.fill(HIST("h1D_EventCounter_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 0)
      histosEvent.fill(HIST("h1D_nTracks0Global_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 1)
      histosEvent.fill(HIST("h1D_nTracks1Global_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 2)
      histosEvent.fill(HIST("h1D_nTracks2Global_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 3)
      histosEvent.fill(HIST("h1D_nTracks3Global_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 4)
      histosEvent.fill(HIST("h1D_nTracks4Global_vs_BC"), collBC);
    if (nTracksGlobalAccepted == 5)
      histosEvent.fill(HIST("h1D_nTracks5Global_vs_BC"), collBC);

    // QA mean pT for a given number of v contributors without TOF && TRD
    if ((counterPVcontributorsNoTOFandTRD == counterPVcontributorsAfterTPCcuts) && counterPVcontributorsNoTOFandTRD > 0) {
      double meanPt = meanPtForPVContributorsNoTOFandTRD / counterPVcontributorsNoTOFandTRD;
      if (counterPVcontributorsNoTOFandTRD == 1)
        histosEvent.fill(HIST("h1D_vContributors_1_meanPt_vs_BC"), collBC, meanPt);
      else if (counterPVcontributorsNoTOFandTRD == 2)
        histosEvent.fill(HIST("h1D_vContributors_2_meanPt_vs_BC"), collBC, meanPt);
      else if (counterPVcontributorsNoTOFandTRD == 3)
        histosEvent.fill(HIST("h1D_vContributors_3_meanPt_vs_BC"), collBC, meanPt);
      else if (counterPVcontributorsNoTOFandTRD == 4)
        histosEvent.fill(HIST("h1D_vContributors_4_meanPt_vs_BC"), collBC, meanPt);
      else if (counterPVcontributorsNoTOFandTRD == 5)
        histosEvent.fill(HIST("h1D_vContributors_5_meanPt_vs_BC"), collBC, meanPt);
    }

    // nTracks vs coll index
    // histosEvent.fill(HIST("h2D_nTracksBeforeCuts_vs_collIndex"), collision.index(), nTracksAll);
    // histosEvent.fill(HIST("h2D_nTracksAfterEtaTPCcuts_vs_collIndex"), collision.index(), nTracksAfterEtaTPCCuts);
    // histosEvent.fill(HIST("h2D_nTracksITSonly_vs_collIndex"), collision.index(), nTracksITSonly);
    // histosEvent.fill(HIST("h2D_nTracksWithITS_vs_collIndex"), collision.index(), nTracksWithITS);
    // histosEvent.fill(HIST("h2D_nTracksWithITS7hits_vs_collIndex"), collision.index(), nTracksWithITS7hits);
    // histosEvent.fill(HIST("h2D_nTracksWithITSandTPC_vs_collIndex"), collision.index(), nTracksWithITSandTPC);
    // histosEvent.fill(HIST("h2D_nTracksWithTRD_vs_collIndex"), collision.index(), nTracksWithTRD);
    // histosEvent.fill(HIST("h2D_nTracksWithTOF_vs_collIndex"), collision.index(), nTracksWithTOF);
    // histosEvent.fill(HIST("h2D_nTracksWithTRDorTOF_vs_collIndex"), collision.index(), nTracksWithTRDorTOF);
    // histosEvent.fill(HIST("h2D_nTracksGlobal_vs_collIndex"), collision.index(), nTracksGlobalAccepted);
    // histosEvent.fill(HIST("h2D_nTracksGlobalWithITS7hits_vs_collIndex"), collision.index(), nTracksGlobalWithITS7hits);
    // histosEvent.fill(HIST("h2D_nTracksGlobalWithTRDorTOF_vs_collIndex"), collision.index(), nTracksGlobalWithTRDorTOF);

    // test RO frame effect
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_After_ITS_ROF_cut"), bcInTF, counterPVcontributorsAfterTPCcuts);

      histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 4, tracks.size());            // ALL tracks outside ITS ROF cut
      histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 4, nTracksGlobalAccepted); // global tracks
      histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 4, 1);                  // collisions counter

      histosEvent.fill(HIST("hCollIndex_ITS_ROF_cut"), collision.index());
      histosEvent.fill(HIST("hCollTime_ITS_ROF_cut"), collision.collisionTime());
      histosEvent.fill(HIST("hCollTimeRes_ITS_ROF_cut"), collision.collisionTimeRes());
      // histosEvent.fill(HIST("h2D_numContrib_vs_collIndex_ITS_ROF_cut"), collision.index(), collision.numContrib());
      histosEvent.fill(HIST("hBC_ITS_ROF_cut"), collBC);

      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_ITSROFcut"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_ITSROFcut"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_ITSROFcut"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_ITSROFcut"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_ITSROFcut"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_ITSROFcut"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_ITSROFcut"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_ITSROFcut"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_ITSROFcut"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_ITSROFcut"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_ITSROFcut"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_ITSROFcut"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_ITSROFcut"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_ITSROFcut"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_ITSROFcut"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_ITSROFcut"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_ITSROFcut"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_ITSROFcut"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_ITSROFcut"), t0cCentr, multT0C);
      }

      if (isFT0) {
        histosFT0.fill(HIST("hVertex_T0_PV_after_ITSROFcut"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_after_ITSROFcut"), diff_PV_ft0_tracks);
      }

      // in addition: TF border cuts:
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        histosEvent.fill(HIST("multAllTr_vs_multT0CAft_ITSROF_TF_cuts"), multT0C, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multV0AAft_ITSROF_TF_cuts"), multV0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multT0AAft_ITSROF_TF_cuts"), multT0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_ITSROF_TF_cuts"), multNTracksPV, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_ITSROF_TF_cuts"), multT0C, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_ITSROF_TF_cuts"), multV0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_ITSROF_TF_cuts"), multT0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_ITSROF_TF_cuts"), multNTracksPV, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_ITSROF_TF_cuts"), multT0C, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_ITSROF_TF_cuts"), multV0A, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_ITSROF_TF_cuts"), multT0A, multNTracksPV);
        histosEvent.fill(HIST("multT0C_vs_multT0AAft_ITSROF_TF_cuts"), multT0A, multT0C);
        histosEvent.fill(HIST("multV0A_vs_multT0AAft_ITSROF_TF_cuts"), multT0A, multV0A);
        histosEvent.fill(HIST("multV0A_vs_multT0CAft_ITSROF_TF_cuts"), multT0C, multV0A);
        if (flagPbPb) {
          histosEvent.fill(HIST("multAllTr_vs_CentAft_ITSROF_TF_cuts"), t0cCentr, multTrk);
          histosEvent.fill(HIST("multGlobalTr_vs_CentAft_ITSROF_TF_cuts"), t0cCentr, nTracksGlobalAccepted);
          histosEvent.fill(HIST("multTrkPV_vs_CentAft_ITSROF_TF_cuts"), t0cCentr, multNTracksPV);
          histosEvent.fill(HIST("multV0A_vs_CentAft_ITSROF_TF_cuts"), t0cCentr, multV0A);
          histosEvent.fill(HIST("multT0C_vs_CentAft_ITSROF_TF_cuts"), t0cCentr, multT0C);
        }

        if (isFT0) {
          histosFT0.fill(HIST("hVertex_T0_PV_after_ITSROF_and_TFcut"), ft0_posZ, collision.posZ());
          histosFT0.fill(HIST("hT0vertexDiff_after_ITSROF_and_TFcut"), diff_PV_ft0_tracks);
        }
      }
    }

    // global PV contributors
    if (nTracksGlobalPVAccepted >= 2) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_2globalPVcontrib"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_2globalPVcontrib"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_2globalPVcontrib"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_2globalPVcontrib"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_2globalPVcontrib"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_2globalPVcontrib"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_2globalPVcontrib"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_2globalPVcontrib"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_2globalPVcontrib"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_2globalPVcontrib"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_2globalPVcontrib"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_2globalPVcontrib"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_2globalPVcontrib"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_2globalPVcontrib"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_2globalPVcontrib"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_2globalPVcontrib"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_2globalPVcontrib"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_2globalPVcontrib"), t0cCentr, multT0C);
      }
    }
    // global PV contributors with 7 ITS hits
    if (nTracksGlobalPVwithITS7hits >= 2) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_2globalPVcontrib_ITS7hits"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_2globalPVcontrib_ITS7hits"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_2globalPVcontrib_ITS7hits"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_2globalPVcontrib_ITS7hits"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_2globalPVcontrib_ITS7hits"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_2globalPVcontrib_ITS7hits"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_2globalPVcontrib_ITS7hits"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib_ITS7hits"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_2globalPVcontrib_ITS7hits"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_2globalPVcontrib_ITS7hits"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_2globalPVcontrib_ITS7hits"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_2globalPVcontrib_ITS7hits"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_2globalPVcontrib_ITS7hits"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_2globalPVcontrib_ITS7hits"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_2globalPVcontrib_ITS7hits"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_2globalPVcontrib_ITS7hits"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_2globalPVcontrib_ITS7hits"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_2globalPVcontrib_ITS7hits"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_2globalPVcontrib_ITS7hits"), t0cCentr, multT0C);
      }
    }

    // global PV contributors with TRD or TOF signal
    if (nTracksGlobalPVwithTRDorTOF >= 2) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_2globalPVcontrib_TRDorTOF"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_2globalPVcontrib_TRDorTOF"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_2globalPVcontrib_TRDorTOF"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_2globalPVcontrib_TRDorTOF"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_2globalPVcontrib_TRDorTOF"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_2globalPVcontrib_TRDorTOF"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_2globalPVcontrib_TRDorTOF"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_2globalPVcontrib_TRDorTOF"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_2globalPVcontrib_TRDorTOF"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_2globalPVcontrib_TRDorTOF"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_2globalPVcontrib_TRDorTOF"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_2globalPVcontrib_TRDorTOF"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_2globalPVcontrib_TRDorTOF"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_2globalPVcontrib_TRDorTOF"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_2globalPVcontrib_TRDorTOF"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_2globalPVcontrib_TRDorTOF"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_2globalPVcontrib_TRDorTOF"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_2globalPVcontrib_TRDorTOF"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_2globalPVcontrib_TRDorTOF"), t0cCentr, multT0C);
      }
    }

    // diffFoundBC_vs_BC cut
    if (diffFoundBC_vs_BC == 0) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_diffFoundBC_vs_BC_0"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_diffFoundBC_vs_BC_0"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_diffFoundBC_vs_BC_0"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_diffFoundBC_vs_BC_0"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_diffFoundBC_vs_BC_0"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_diffFoundBC_vs_BC_0"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_diffFoundBC_vs_BC_0"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_diffFoundBC_vs_BC_0"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_diffFoundBC_vs_BC_0"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_diffFoundBC_vs_BC_0"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_diffFoundBC_vs_BC_0"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_diffFoundBC_vs_BC_0"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_diffFoundBC_vs_BC_0"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_diffFoundBC_vs_BC_0"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_diffFoundBC_vs_BC_0"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_diffFoundBC_vs_BC_0"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_diffFoundBC_vs_BC_0"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_diffFoundBC_vs_BC_0"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_diffFoundBC_vs_BC_0"), t0cCentr, multT0C);
      }
    }

    // FT0 present
    if (isFT0) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_hasFT0_CorrectedValid"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_hasFT0_CorrectedValid"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_hasFT0_CorrectedValid"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_hasFT0_CorrectedValid"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_hasFT0_CorrectedValid"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_hasFT0_CorrectedValid"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_hasFT0_CorrectedValid"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_hasFT0_CorrectedValid"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_hasFT0_CorrectedValid"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_hasFT0_CorrectedValid"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_hasFT0_CorrectedValid"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_hasFT0_CorrectedValid"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_hasFT0_CorrectedValid"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_hasFT0_CorrectedValid"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_hasFT0_CorrectedValid"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_hasFT0_CorrectedValid"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_hasFT0_CorrectedValid"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_hasFT0_CorrectedValid"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_hasFT0_CorrectedValid"), t0cCentr, multT0C);
      }
    }

    // cut on diff b/n vertex from FT0 and track-based
    if (isFT0 && fabs(diff_PV_ft0_tracks) < 1) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_after_PV_FT0_diff_cut_1cm"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_after_PV_FT0_diff_cut_1cm"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_after_PV_FT0_diff_cut_1cm"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_1cm"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_after_PV_FT0_diff_cut_1cm"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_after_PV_FT0_diff_cut_1cm"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_after_PV_FT0_diff_cut_1cm"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_1cm"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_after_PV_FT0_diff_cut_1cm"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_after_PV_FT0_diff_cut_1cm"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_after_PV_FT0_diff_cut_1cm"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_after_PV_FT0_diff_cut_1cm"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_after_PV_FT0_diff_cut_1cm"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_after_PV_FT0_diff_cut_1cm"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_after_PV_FT0_diff_cut_1cm"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_after_PV_FT0_diff_cut_1cm"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_after_PV_FT0_diff_cut_1cm"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_after_PV_FT0_diff_cut_1cm"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_after_PV_FT0_diff_cut_1cm"), t0cCentr, multT0C);
      }

      // ALL OTHER CUTS:
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        histosEvent.fill(HIST("multAllTr_vs_multT0CAft_ALL_CUTS"), multT0C, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multV0AAft_ALL_CUTS"), multV0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multT0AAft_ALL_CUTS"), multT0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_ALL_CUTS"), multNTracksPV, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_ALL_CUTS"), multT0C, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_ALL_CUTS"), multV0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_ALL_CUTS"), multT0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_ALL_CUTS"), multNTracksPV, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_ALL_CUTS"), multT0C, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_ALL_CUTS"), multV0A, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_ALL_CUTS"), multT0A, multNTracksPV);
        histosEvent.fill(HIST("multT0C_vs_multT0AAft_ALL_CUTS"), multT0A, multT0C);
        histosEvent.fill(HIST("multV0A_vs_multT0AAft_ALL_CUTS"), multT0A, multV0A);
        histosEvent.fill(HIST("multV0A_vs_multT0CAft_ALL_CUTS"), multT0C, multV0A);
        if (flagPbPb) {
          histosEvent.fill(HIST("multAllTr_vs_CentAft_ALL_CUTS"), t0cCentr, multTrk);
          histosEvent.fill(HIST("multGlobalTr_vs_CentAft_ALL_CUTS"), t0cCentr, nTracksGlobalAccepted);
          histosEvent.fill(HIST("multTrkPV_vs_CentAft_ALL_CUTS"), t0cCentr, multNTracksPV);
          histosEvent.fill(HIST("multV0A_vs_CentAft_ALL_CUTS"), t0cCentr, multV0A);
          histosEvent.fill(HIST("multT0C_vs_CentAft_ALL_CUTS"), t0cCentr, multT0C);
        }
      }
    }

    // cut on diff b/n vertex from FT0 and track-based - TIGHTER
    if (isFT0 && diff_PV_ft0_tracks > -0.8 && diff_PV_ft0_tracks < 0.) {
      histosEvent.fill(HIST("multAllTr_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT"), multT0C, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT"), multV0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT"), multT0A, multTrk);
      histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_TIGHT"), multNTracksPV, multTrk);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT"), multT0C, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT"), multV0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT"), multT0A, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_after_PV_FT0_diff_cut_TIGHT"), multNTracksPV, nTracksGlobalAccepted);
      histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT"), multT0C, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_after_PV_FT0_diff_cut_TIGHT"), multV0A, multNTracksPV);
      histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT"), multT0A, multNTracksPV);
      histosEvent.fill(HIST("multT0C_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT"), multT0A, multT0C);
      histosEvent.fill(HIST("multV0A_vs_multT0AAft_after_PV_FT0_diff_cut_TIGHT"), multT0A, multV0A);
      histosEvent.fill(HIST("multV0A_vs_multT0CAft_after_PV_FT0_diff_cut_TIGHT"), multT0C, multV0A);
      if (flagPbPb) {
        histosEvent.fill(HIST("multAllTr_vs_CentAft_after_PV_FT0_diff_cut_TIGHT"), t0cCentr, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_CentAft_after_PV_FT0_diff_cut_TIGHT"), t0cCentr, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_CentAft_after_PV_FT0_diff_cut_TIGHT"), t0cCentr, multNTracksPV);
        histosEvent.fill(HIST("multV0A_vs_CentAft_after_PV_FT0_diff_cut_TIGHT"), t0cCentr, multV0A);
        histosEvent.fill(HIST("multT0C_vs_CentAft_after_PV_FT0_diff_cut_TIGHT"), t0cCentr, multT0C);
      }

      // ALL OTHER CUTS TIGHTER:
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        histosEvent.fill(HIST("multAllTr_vs_multT0CAft_ALL_CUTS_Tighter"), multT0C, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multV0AAft_ALL_CUTS_Tighter"), multV0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multT0AAft_ALL_CUTS_Tighter"), multT0A, multTrk);
        histosEvent.fill(HIST("multAllTr_vs_multTrkPVAft_ALL_CUTS_Tighter"), multNTracksPV, multTrk);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0CAft_ALL_CUTS_Tighter"), multT0C, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multV0AAft_ALL_CUTS_Tighter"), multV0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multT0AAft_ALL_CUTS_Tighter"), multT0A, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multGlobalTr_vs_multTrkPVAft_ALL_CUTS_Tighter"), multNTracksPV, nTracksGlobalAccepted);
        histosEvent.fill(HIST("multTrkPV_vs_multT0CAft_ALL_CUTS_Tighter"), multT0C, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multV0AAft_ALL_CUTS_Tighter"), multV0A, multNTracksPV);
        histosEvent.fill(HIST("multTrkPV_vs_multT0AAft_ALL_CUTS_Tighter"), multT0A, multNTracksPV);
        histosEvent.fill(HIST("multT0C_vs_multT0AAft_ALL_CUTS_Tighter"), multT0A, multT0C);
        histosEvent.fill(HIST("multV0A_vs_multT0AAft_ALL_CUTS_Tighter"), multT0A, multV0A);
        histosEvent.fill(HIST("multV0A_vs_multT0CAft_ALL_CUTS_Tighter"), multT0C, multV0A);
        if (flagPbPb) {
          histosEvent.fill(HIST("multAllTr_vs_CentAft_ALL_CUTS_Tighter"), t0cCentr, multTrk);
          histosEvent.fill(HIST("multGlobalTr_vs_CentAft_ALL_CUTS_Tighter"), t0cCentr, nTracksGlobalAccepted);
          histosEvent.fill(HIST("multTrkPV_vs_CentAft_ALL_CUTS_Tighter"), t0cCentr, multNTracksPV);
          histosEvent.fill(HIST("multV0A_vs_CentAft_ALL_CUTS_Tighter"), t0cCentr, multV0A);
          histosEvent.fill(HIST("multT0C_vs_CentAft_ALL_CUTS_Tighter"), t0cCentr, multT0C);
        }
      }
    }

    // only verteces with nContr>=3 with 7 ITS clusters
    if (counterPVcontributorsITS7hits >= 3) {
      histosEvent.fill(HIST("hCollIndex_vertNcontr3_withITS7hits"), collision.index());
      histosEvent.fill(HIST("hCollTime_vertNcontr3_withITS7hits"), collision.collisionTime());
      histosEvent.fill(HIST("hCollTimeRes_vertNcontr3_withITS7hits"), collision.collisionTimeRes());
      histosEvent.fill(HIST("hBC_vertNcontr3_withITS7hits"), collBC);
    }
    histosEvent.fill(HIST("hNumContrib_vertNcontr3_withITS7hits"), counterPVcontributorsITS7hits);

    // only verteces with nContr>=3 that have TOF or TRD
    if (counterVertexContributorsWithTRDorTOF >= 3) {
      histosEvent.fill(HIST("hCollIndex_vertNcontr3_TRDorTOF"), collision.index());
      histosEvent.fill(HIST("hCollTime_vertNcontr3_TRDorTOF"), collision.collisionTime());
      histosEvent.fill(HIST("hCollTimeRes_vertNcontr3_TRDorTOF"), collision.collisionTimeRes());
      histosEvent.fill(HIST("vtxChi2_vertNcontr3_TRDorTOF"), collision.chi2());
      histosEvent.fill(HIST("hBC_vertNcontr3_TRDorTOF"), collBC);
    }
    histosEvent.fill(HIST("hNumContrib_vertNcontr3_TRDorTOF"), counterVertexContributorsWithTRDorTOF);

    // only verteces with nContr>=3 with 7 ITS clusters + has TOF or TRD
    if (counterVertexContributorsWithTRDorTOF_ITS7hits >= 3) {
      histosEvent.fill(HIST("hCollIndex_vertNcontr3_withITS7hits_and_TRDorTOF"), collision.index());
      histosEvent.fill(HIST("hCollTime_vertNcontr3_withITS7hits_and_TRDorTOF"), collision.collisionTime());
      histosEvent.fill(HIST("hCollTimeRes_vertNcontr3_withITS7hits_and_TRDorTOF"), collision.collisionTimeRes());
      histosEvent.fill(HIST("hBC_vertNcontr3_withITS7hits_and_TRDorTOF"), collBC);
    }
    histosEvent.fill(HIST("hNumContrib_vertNcontr3_withITS7hits_and_TRDorTOF"), counterVertexContributorsWithTRDorTOF_ITS7hits);

    // `tracks` contains tracks belonging to `collision`
    // `v0s`    contains v0s    belonging to `collision`
    if (flagShowInfo)
      LOGF(info, "#########");
    if (flagShowInfo)
      LOGF(info, "Collision index : %d", collision.index());
    if (flagShowInfo)
      LOGF(info, "Number of tracks: %d", tracks.size());
    if (flagShowInfo)
      LOGF(info, "Number of v0s   : %d", v0s.size());

    // float vtxZ = 10.;
    if (flagShowInfo)
      LOGF(info, "Collision: %d [N = %d out of %d], vZ pos = %.1f ",                      //, globalBC=%.1f",
           collision.globalIndex(), tracks.size(), tracks.tableSize(), collision.posZ()); //, collision.globalBC() );

    // cout << //"collision.bc() = " << collision.bc().globalBC()  <<
    // "runNumber = " << collision.bc().runNumber() << endl;

    if (flagShowInfo)
      LOGF(info, "The collision time:   %f", collision.collisionTime());
    if (flagShowInfo)
      LOGF(info, "Tracks for this collision: %d", tracks.size());

    bool flagShowInfoFirstTrack = false;

    // #### loop over tracks
    for (auto const& track : tracks) {
      // eta
      histosTracks.fill(HIST("etaHistogram"), track.eta());
      if (fabs(track.eta()) > 0.8)
        continue;
      histosTracks.fill(HIST("etaHistogramAfter08cut"), track.eta());

      histosTracks.fill(HIST("hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());

      histosTracks.fill(HIST("hNumITSclusters"), track.itsNCls());
      histosTracks.fill(HIST("hDcaXY"), track.dcaXY());

      histosTracks.fill(HIST("ptHistogram"), track.pt());
      histosTracks.fill(HIST("phiHistogram"), track.phi());

      histosTracks.fill(HIST("hTrackLength"), track.length()); // from TrackExtras

      if (track.itsNCls() == 7)
        histosTracks.fill(HIST("hTpcNClsCrossedRowsITS7hits"), track.tpcNClsCrossedRows());

      if (flagShowInfoFirstTrack) {
        cout << "track.length() = " << track.length() << endl;
        flagShowInfoFirstTrack = false;
      }
      const float combNSigmaEl = std::sqrt(track.tofNSigmaEl() * track.tofNSigmaEl() + track.tpcNSigmaEl() * track.tpcNSigmaEl());
      // cout << "combNSigmaEl = " << combNSigmaEl << endl;
      histosTracks.fill(HIST("pidCombSigma"), combNSigmaEl);

      // ### my track cuts
      // if (fabs(track.eta()) > 0.8)
      // continue;
      if (track.tpcNClsCrossedRows() < 110)
        continue;

      histosTracks.fill(HIST("hChi2TPCperDOF"), track.tpcChi2NCl());
      histosTracks.fill(HIST("hChi2ITSperDOF"), track.itsChi2NCl());
      histosTracks.fill(HIST("hTrackTime"), track.trackTime());
      histosTracks.fill(HIST("hTrackTimeRes"), track.trackTimeRes());

      histosTracks.fill(HIST("hTpcNClsCrossedRows_AftCuts"), track.tpcNClsCrossedRows());
      histosTracks.fill(HIST("hTpcNClsFound_AftCuts"), track.tpcNClsFound());

      // QA TPC boundaries
      if (track.isGlobalTrack()) {
        double pt = track.pt();
        double eta = track.eta();
        double phi = track.phi();

        if (track.sign() > 0) {
          histosTracks.fill(HIST("posSelTrack_pt"), pt);
          histosTracks.fill(HIST("posSelTrack_phi"), phi);
          histosTracks.fill(HIST("posSelTrack_eta"), eta);

          histosTracks.fill(HIST("posSelTrack_phi_vs_pt"), pt, phi);
          histosTracks.fill(HIST("posSelTrack_phi_vs_pt_modPiOver9"), pt, fmod(phi + constPhiShift, TMath::Pi() / 9));

          // MFsign
          if ((pt > ptTPCsectorCut && (fmod(phi + constPhiShift, TMath::Pi() / 9) > (MFsign < 0 ? fPhiCutExpPosHigh : fPhiCutExpNegHigh)->Eval(pt) || fmod(phi + constPhiShift, TMath::Pi() / 9) < (MFsign < 0 ? fPhiCutExpPosLow : fPhiCutExpNegLow)->Eval(pt))) || pt <= ptTPCsectorCut)
            histosTracks.fill(HIST("posSelTrack_phi_vs_pt_afterCut"), pt, phi);

          if (pt > 0.5 && pt < 1.0)
            histosTracks.fill(HIST("posSelTrack_phi_pT_05_10"), phi);
          else if (pt > 1.0 && pt < 2.0)
            histosTracks.fill(HIST("posSelTrack_phi_pT_10_20"), phi);
          else if (pt > 2.0 && pt < 10.0)
            histosTracks.fill(HIST("posSelTrack_phi_pT_20_100"), phi);
        } else if (track.sign() < 0) {
          histosTracks.fill(HIST("negSelTrack_pt"), pt);
          histosTracks.fill(HIST("negSelTrack_phi"), phi);
          histosTracks.fill(HIST("negSelTrack_eta"), eta);

          histosTracks.fill(HIST("negSelTrack_phi_vs_pt"), pt, phi);
          histosTracks.fill(HIST("negSelTrack_phi_vs_pt_modPiOver9"), pt, fmod(phi + constPhiShift, TMath::Pi() / 9));
          if ((pt > ptTPCsectorCut && (fmod(phi + constPhiShift, TMath::Pi() / 9) > (MFsign < 0 ? fPhiCutExpNegHigh : fPhiCutExpPosHigh)->Eval(pt) || fmod(phi + constPhiShift, TMath::Pi() / 9) < (MFsign < 0 ? fPhiCutExpNegLow : fPhiCutExpPosLow)->Eval(pt))) || pt <= ptTPCsectorCut)
            histosTracks.fill(HIST("negSelTrack_phi_vs_pt_afterCut"), pt, phi);

          if (pt > 0.5 && pt < 1.0)
            histosTracks.fill(HIST("negSelTrack_phi_pT_05_10"), phi);
          else if (pt > 1.0 && pt < 2.0)
            histosTracks.fill(HIST("negSelTrack_phi_pT_10_20"), phi);
          else if (pt > 2.0 && pt < 10.0)
            histosTracks.fill(HIST("negSelTrack_phi_pT_20_100"), phi);
        }
      }

      if (track.itsNCls() == 7) {
        histosTracks.fill(HIST("hTpcNClsCrossedRowsITS7hits_AftCuts"), track.tpcNClsCrossedRows());
        histosTracks.fill(HIST("hTpcNClsFound_ITS7hits_AftCuts"), track.tpcNClsFound());
        histosTracks.fill(HIST("etaHistogramITS7hits_AftCuts"), track.eta());

        histosTracks.fill(HIST("hChi2TPCperDOF_ITS7hits"), track.tpcChi2NCl());
        histosTracks.fill(HIST("hChi2ITSperDOF_ITS7hits"), track.itsChi2NCl());
        histosTracks.fill(HIST("hTrackTime_ITS7hits"), track.trackTime());
        histosTracks.fill(HIST("hTrackTimeRes_ITS7hits"), track.trackTimeRes());
      }
      if (track.hasTRD() || track.hasTOF()) {
        histosTracks.fill(HIST("hChi2TPCperDOF_TRDorTOF"), track.tpcChi2NCl());
        histosTracks.fill(HIST("hChi2ITSperDOF_TRDorTOF"), track.itsChi2NCl());
        histosTracks.fill(HIST("hTrackTime_TRDorTOF"), track.trackTime());
        histosTracks.fill(HIST("hTrackTimeResTRDorTOF"), track.trackTimeRes());
      }

      histosTracks.fill(HIST("hNumITSclusters_AftCuts"), track.itsNCls());
      histosTracks.fill(HIST("hDcaXY_AftCuts"), track.dcaXY());

      histosTracks.fill(HIST("etaHistogram_AftCuts"), track.eta());

      histosTracks.fill(HIST("ptHistogram_AftCuts"), track.pt());
      histosTracks.fill(HIST("phiHistogram_AftCuts"), track.phi());

      histosTracks.fill(HIST("hTrackLength_AftCuts"), track.length()); // from TrackExtras
      // histosTracks.fill(HIST("hChi2perDOF_AftCuts"), track.chi2perNDF());

    } // end of track loop

    // ##### v0 analysis
    for (const auto& v0 : v0s) {
      double v0rad = v0.v0radius();
      histosK0S.fill(HIST("hK0Sradius"), v0rad);

      // Cut on dynamic columns
      if (v0.v0cosPA() < 0.99) // v0setting_cospa)
        continue;

      histosK0S.fill(HIST("hMassK0Short"), v0.mK0Short());

      if (v0rad < 0.8 || v0rad > 20) // v0setting_radius)
        continue;

      // if (v0rad < v0RadiusMin || v0rad > v0RadiusMax || v0.eta() > assocEtaMax || v0.eta() < assocEtaMin || v0.v0cosPA() < v0Cospa) {
      if (fabs(v0.eta()) > 0.75)
        continue;

      const auto& posDaughterTrack = v0.posTrack_as<DaughterTracks>();
      const auto& negDaughterTrack = v0.negTrack_as<DaughterTracks>();

      double etaDauPos = posDaughterTrack.eta();
      double etaDauNeg = negDaughterTrack.eta();

      if (fabs(etaDauPos) > 0.8 || fabs(etaDauNeg) > 0.8)
        continue;

      histosK0S.fill(HIST("hMassK0ShortAfterSelection"), v0.mK0Short());
      histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPt"), v0.pt(), v0.mK0Short());

      if (TMath::Abs(posDaughterTrack.tpcNSigmaPi()) > 3) // NSigmaTPCPion)
        continue;
      if (TMath::Abs(negDaughterTrack.tpcNSigmaPi()) > 3) // NSigmaTPCPion)
        continue;

      double ptDauPos = posDaughterTrack.pt();
      double phiDauPos = posDaughterTrack.phi();
      double ptDauNeg = negDaughterTrack.pt();
      double phiDauNeg = negDaughterTrack.phi();

      // QA daughter track params
      if (posDaughterTrack.tpcNClsFound() > 80 && negDaughterTrack.tpcNClsFound() > 80) {
        // pos
        histosK0S.fill(HIST("posDau_pt"), ptDauPos);
        histosK0S.fill(HIST("posDau_phi"), phiDauPos);
        histosK0S.fill(HIST("posDau_eta"), etaDauPos);
        histosK0S.fill(HIST("posDau_phi_vs_pt"), ptDauPos, phiDauPos);

        double dauPosTPCcutOk = false;
        if ((ptDauPos > ptTPCsectorCut && (fmod(phiDauPos + constPhiShift, TMath::Pi() / 9) > (MFsign < 0 ? fPhiCutExpPosHigh : fPhiCutExpNegHigh)->Eval(ptDauPos) || fmod(phiDauPos + constPhiShift, TMath::Pi() / 9) < (MFsign < 0 ? fPhiCutExpPosLow : fPhiCutExpNegLow)->Eval(ptDauPos))) || ptDauPos <= ptTPCsectorCut)
          dauPosTPCcutOk = true;

        if (dauPosTPCcutOk)
          histosK0S.fill(HIST("posDau_phi_vs_pt_afterCut"), ptDauPos, phiDauPos);

        if (ptDauPos > 0.5 && ptDauPos < 1.0)
          histosK0S.fill(HIST("posDau_phi_pT_05_10"), phiDauPos);
        else if (ptDauPos > 1.0 && ptDauPos < 2.0)
          histosK0S.fill(HIST("posDau_phi_pT_10_20"), phiDauPos);
        else if (ptDauPos > 2.0 && ptDauPos < 10.0)
          histosK0S.fill(HIST("posDau_phi_pT_20_100"), phiDauPos);

        // neg
        histosK0S.fill(HIST("negDau_pt"), ptDauNeg);
        histosK0S.fill(HIST("negDau_phi"), phiDauNeg);
        histosK0S.fill(HIST("negDau_eta"), etaDauNeg);
        histosK0S.fill(HIST("negDau_phi_vs_pt"), ptDauNeg, phiDauNeg);

        double dauNegTPCcutOk = false;
        if ((ptDauNeg > ptTPCsectorCut && (fmod(phiDauPos + constPhiShift, TMath::Pi() / 9) > (MFsign < 0 ? fPhiCutExpNegHigh : fPhiCutExpPosHigh)->Eval(ptDauNeg) || fmod(phiDauPos + constPhiShift, TMath::Pi() / 9) < (MFsign < 0 ? fPhiCutExpNegLow : fPhiCutExpPosLow)->Eval(ptDauNeg))) || ptDauNeg <= ptTPCsectorCut)
          dauNegTPCcutOk = true;

        if (dauNegTPCcutOk)
          histosK0S.fill(HIST("negDau_phi_vs_pt_afterCut"), ptDauNeg, phiDauPos);

        if (ptDauPos > 0.5 && ptDauPos < 1.0)
          histosK0S.fill(HIST("negDau_phi_pT_05_10"), phiDauNeg);
        else if (ptDauPos > 1.0 && ptDauPos < 2.0)
          histosK0S.fill(HIST("negDau_phi_pT_10_20"), phiDauNeg);
        else if (ptDauPos > 2.0 && ptDauPos < 10.0)
          histosK0S.fill(HIST("negDau_phi_pT_20_100"), phiDauNeg);

        // accept if both daughters are not at the TPC boundaries
        if (dauPosTPCcutOk && dauNegTPCcutOk)
          histosK0S.fill(HIST("hMassK0S_cut_TPC_boundaries_for_Dau"), v0.pt(), v0.mK0Short());
      }

      if (etaDauPos > 0.1 && etaDauPos<0.8 & etaDauNeg> 0.1 && etaDauNeg < 0.8)
        histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtEta01_08"), v0.pt(), v0.mK0Short());
      if (etaDauPos > -0.8 && etaDauPos<-0.1 & etaDauNeg> - 0.8 && etaDauNeg < -0.1)
        histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtEta08_01"), v0.pt(), v0.mK0Short());

      histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtPIDnSigma"), v0.pt(), v0.mK0Short());

      if (collision.numContrib() > 4) {
        histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4"), v0.pt(), v0.mK0Short());
        // if (rofCut){
        if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
          histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4_andROFcut"), v0.pt(), v0.mK0Short());
      }

      if (counterPVcontributorsITS7hits >= 3)
        histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtPIDnSigmaNITS7hitsContribAbove2"), v0.pt(), v0.mK0Short());
      if (counterVertexContributorsWithTRDorTOF >= 3)
        histosK0S.fill(HIST("hMassK0S_nVertexContributorsWithTRDorTOF_Above2"), v0.pt(), v0.mK0Short());
      if (counterVertexContributorsWithTRDorTOF_ITS7hits >= 3)
        histosK0S.fill(HIST("hMassK0S_nVertexContributorsWithTRDorTOF_ITS7hits_Above2"), v0.pt(), v0.mK0Short());

      // more cuts on daughters TPC clusters and eta:
      if (counterVertexContributorsWithTRDorTOF >= 3) {
        if (posDaughterTrack.tpcNClsFound() > 80 && negDaughterTrack.tpcNClsFound() > 80) {
          if (fabs(etaDauPos) < 0.8 && fabs(etaDauNeg) < 0.8 && fabs(etaDauPos) > 0.1 && fabs(etaDauNeg) > 0.1) {
            histosK0S.fill(HIST("hMassK0S_cutsOnDaughters"), v0.pt(), v0.mK0Short());

            // dependance on K0S radius:
            if (v0rad < 2)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_v0Rad_0_2"), v0.pt(), v0.mK0Short());
            else if (v0rad > 2 && v0rad < 4)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_v0Rad_2_4"), v0.pt(), v0.mK0Short());
            else if (v0rad > 4 && v0rad < 6)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_v0Rad_4_6"), v0.pt(), v0.mK0Short());
            else if (v0rad > 6 && v0rad < 10)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_v0Rad_6_10"), v0.pt(), v0.mK0Short());
            else if (v0rad > 10 && v0rad < 20)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_v0Rad_10_20"), v0.pt(), v0.mK0Short());

            if (posDaughterTrack.hasTOF() && negDaughterTrack.hasTOF())
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_bothTOF"), v0.pt(), v0.mK0Short());
            if (posDaughterTrack.hasTRD() && negDaughterTrack.hasTRD())
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_bothTRD"), v0.pt(), v0.mK0Short());

            if ((posDaughterTrack.hasTOF() || posDaughterTrack.hasTRD()) && (negDaughterTrack.hasTOF() || negDaughterTrack.hasTRD()))
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_TOForTRD"), v0.pt(), v0.mK0Short());
            else
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_norTOFnorTRD"), v0.pt(), v0.mK0Short());

            // check ITS clusters
            if (posDaughterTrack.itsNCls() >= 4 && negDaughterTrack.itsNCls() >= 4)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_ITSmin4hits"), v0.pt(), v0.mK0Short());
            if (posDaughterTrack.itsNCls() >= 5 && negDaughterTrack.itsNCls() >= 4)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_ITSmin5hits"), v0.pt(), v0.mK0Short());
            if (posDaughterTrack.itsNCls() >= 6 && negDaughterTrack.itsNCls() >= 6)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_ITSmin6hits"), v0.pt(), v0.mK0Short());
            if (posDaughterTrack.tpcChi2NCl() < 2 && negDaughterTrack.tpcChi2NCl() < 2)
              histosK0S.fill(HIST("hMassK0S_cutsOnDaughters_chi2tpc2"), v0.pt(), v0.mK0Short());
          }
        }
      }

    } // end of v0 loop
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<RobustFluctuationObservables>(cfgc)};
}
