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
/// \brief This task is a QA task to accumulate basic event- and track-level plots.
/// \author Igor Altsybeev, Igor.Altsybeev@cern.ch

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "DataFormatsParameters/GRPECSObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "TF1.h"
#include "TGraphErrors.h"

#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// declaratives to orginise initialization and filling of histograms for various sets of cuts
#define FILL_QA_HIST_1D(cutId, hName, x...)                  \
  if (!mHist1D.count(hName))                                 \
    LOGF(fatal, "AHTUNG! no key %s in mHist1D map!", hName); \
  pHist1D[cutId][mHist1D[hName]]->Fill(x);

#define FILL_QA_HIST_2D(cutId, hName, x, y...)               \
  if (!mHist2D.count(hName))                                 \
    LOGF(fatal, "AHTUNG! no key %s in mHist2D map!", hName); \
  pHist2D[cutId][mHist2D[hName]]->Fill(x, y);

#define ADD_QA_HIST_1D(hName, title, axisX...) \
  mHist1D[hName] = v1D.size();                 \
  v1D.push_back(histosEventNew.add<TH1>(TString::Format("%s/%s", hName, cutName.c_str()).Data(), title, kTH1D, {axisX}));

#define ADD_QA_HIST_2D(hName, title, axisX, axisY...) \
  mHist2D[hName] = v2D.size();                        \
  v2D.push_back(histosEventNew.add<TH2>(TString::Format("%s/%s", hName, cutName.c_str()).Data(), title, kTH2F, {axisX, axisY}));

// pt study:
#define ADD_PT_HIST_1D(hName, cutName, hTitle, axisX...)                                                                         \
  {                                                                                                                              \
    string str = hName;                                                                                                          \
    string fullName = str + "/" + cutName;                                                                                       \
    mPtStudyCuts[fullName] = pHistPtStudy1D.size();                                                                              \
    pHistPtStudy1D.push_back(histosDetailedPt.add<TH1>(TString::Format("%s", fullName.c_str()).Data(), hTitle, kTH1D, {axisX})); \
  }

#define FILL_PT_HIST_1D(cutName, x, w...)                           \
  if (!mPtStudyCuts.count(cutName))                                 \
    LOGF(fatal, "AHTUNG! no key %s in mPtStudyCuts map!", cutName); \
  pHistPtStudy1D[mPtStudyCuts[cutName]]->Fill(x, w);

// main class
struct RobustFluctuationObservables {
  // for vertex vs time:
  bool flagShowInfo = false;
  int lastRunNumber = -1;
  uint64_t nBCsPerOrbit = 3564;

  // bc position correlations
  int64_t prevOrbit = -1;
  int64_t prevBcInTF = -1;
  int64_t prevFoundBcInTF = -1;
  uint64_t prevBC = 9999;            //-1;
  uint64_t prevTF = 9999;            //-1;
  uint64_t prevGlobalFoundBC = 9999; //-1;
  float prev_vZ = -1;
  float prev_mult = -1;
  float prev_globTrkContrib = -1;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  int64_t bcSOR = -1; // global bc of the start of the first orbit
  int64_t orbitSOR = -1;
  // int64_t bcSORbis = -1; // global bc of the start of the first orbit - try alternative
  int64_t nBCsPerTF = 1; // 128*3564; // duration of TF in bcs
  uint64_t TFid = 0;     // count time frames in a given run
  bool flagWaitForNewTF = false;
  uint32_t nOrbitsPerTF = 0;

  // functions for parametrized cuts on pT vs phi plots
  TF1* fPhiCutExpPosHigh;
  TF1* fPhiCutExpPosLow;
  TF1* fPhiCutExpNegHigh;
  TF1* fPhiCutExpNegLow;
  double constPhiShift = 0.175;
  double ptTPCsectorCut = 0.4;

  // ##### hist registries
  HistogramRegistry histosDF{"histosDF", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEvent{"histosEventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventNew{"histosEventSelectionNew", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventTracksFoundBC{"histosEventTracksFoundBC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventCounters{"histosEventCounters", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosEventBcInTF{"histosEventSelectionBcInTF", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosFT0{"histosFT0", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosTracks{"histosTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosDetailedPt{"histosDetailedPt", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histosK0S{"kzeroShort", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ##### configurables and axes
  Configurable<int> nBinsPt{"nBinsPt", 800, "N bins in pT histo"};
  Configurable<int> nBinsEta{"nBinsEta", 400, "N bins in eta histo"};
  Configurable<int> nBinsPhi{"nBinsPhi", 360, "N bins in phi histo"};
  Configurable<float> vtxCut{"vtxCut", 10.0, "Accepted z-vertex range (cm)"};
  Configurable<int> flagPbPb{"flagPbPb", 0, "0 - pp, 1 - PbPb"};
  Configurable<int> flagIncludeQAHistBcInTF{"flagIncludeQAHistBcInTF", 1, "flagIncludeQAHistBcInTF: 0 - no, 1 - yes"};
  Configurable<int> flagIncludeQAHistBcVsNtracks{"flagIncludeQAHistBcVsNtracks", 1, "flagIncludeQAHistBcVsNtracks: 0 - no, 1 - yes"};
  Configurable<int> flagIncludeQAHistK0S{"flagIncludeQAHistK0S", 1, "flagIncludeQAHistK0S: 0 - no, 1 - yes"};
  Configurable<int> flagIncludeQATracksInFoundBC{"flagIncludeQATracksInFoundBC", 1, "flagIncludeQATracksInFoundBC: 0 - no, 1 - yes"};
  Configurable<int> flagIncludeQAHistPartsOfOrbit{"flagIncludeQAHistPartsOfOrbit", 1, "flagIncludeQAHistPartsOfOrbit: 0 - no, 1 - yes"};
  Configurable<int> flagIncludeQAROFinPtEtaBins{"flagIncludeQAROFinPtEtaBins", 1, "flagIncludeQAROFinPtEtaBins: 0 - no, 1 - yes"};

  Configurable<int> nBinsContrib{"nBinsContrib", 200, "N bins in n vertex contrib histo"};
  Configurable<int> nMaxContrib{"nMaxContrib", 200, "N max in nContrib histo"};

  Configurable<int> nBinsTracks{"nBinsTracks", 200, "N bins in n tracks histo"};
  Configurable<int> nMaxTracks{"nMaxTracks", 200, "N max in n tracks histo"};
  // AxisSpec axisMult{200, 0.f, 5000.f, "multiplicity"};       //{1000, 0.f, 5000.f, "multiplicity"};

  Configurable<int> nBinsCollIndex{"nBinsCollIndex", 200, "nBinsCollIndex"};
  Configurable<int> nMaxCollIndex{"nMaxCollIndex", 10000, "nMaxCollIndex"};

  Configurable<int> nBinsMultFwd{"nBinsMultFwd", 200, "N bins in mult fwd histo"};
  Configurable<float> nMaxMultFwd{"nMaxMultFwd", 20000, "N max in mult fwd histo"};

  // hand-made ITS ROF cut
  Configurable<int> nITSROF{"nITSROF", 6, "nITSROF"};
  Configurable<int> nITSROF_BC_offset{"nITSROF_BC_offset", 65, "nITSROF_BC_offset"};
  Configurable<uint64_t> nITSROF_BC_cutWidth{"nITSROF_BC_cutWidth", 40, "nITSROF_BC_cutWidth"};
  // Configurable<int> nITSROF_middle_cut_forITSonlyVert{"nITSROF_middle_cut_forITSonlyVert", 198/2 /*ROF=198 in pp*/, "nITSROF_middle_cut_forITSonlyVert"};
  // Configurable<int> nNoITSonlyVertices{"nNoITSonlyVertices", false, "nITSROF_middle_cut_forITSonlyVert"};

  // cuts on difference in vZ position by tracks vs T0
  Configurable<float> cutVzTrackT0diffLower{"cutVzTrackT0diffLower", -1., "cutVzTrackT0diffLower, cm"};
  Configurable<float> cutVzTrackT0diffUpper{"cutVzTrackT0diffUpper", 1., "cutVzTrackT0diffUpper, cm"};

  // splitting of the orbit into several BC ranges
  Configurable<std::vector<uint64_t>> vSplitBCpointsOfTheOrbit{"SplitBCpointsOfTheOrbit", {1200, 2000, 3000}, "BC split points of the orbit"};

  // orbit QA
  uint32_t orbitAtCollIndexZero = 0;

  // QA per-DF histograms
  map<uint64_t, uint64_t> mapDF;
  Configurable<int> nHistQAplotsDF{"nDFforQA", 5, "number of per-DF QA histograms"};
  vector<shared_ptr<TH1>> fV_h1D_Orbit_vs_CollIndex; //{(const int)nHistQAplotsDF, nullptr};

  // QA histograms vs various cuts
  map<string, int> mMyCuts;
  map<string, int> mHist1D;
  map<string, int> mHist2D;
  // std::shared_ptr<TH2> pHist[20][20];
  vector<vector<std::shared_ptr<TH1>>> pHist1D;
  vector<vector<std::shared_ptr<TH2>>> pHist2D;

  // QA for pt vs various cuts
  map<string, int> mPtStudyCuts;
  vector<std::shared_ptr<TH1>> pHistPtStudy1D;

  //
  TF1* funcCutEventsByMultPVvsV0A;
  TF1* funcCutEventsByMultPVvsT0C;

  TF1* funcCutEventsByMultPVvsV0A_anotherTrend; // Pb-Pb

  void init(InitContext const&)
  {
    // naming of the 2D correlation plots
    string strCutNames[] =
      {
        "Bef",
        "Aft",
        "TFcut",
        "ITSROFcut",
        // "ITSROF_TF_cuts",
        // "1globalPVcontrib",
        "1globalPVcontrib_ITS7hits",
        "1globalPVcontrib_TRDorTOF",
        "diffFoundBC_vs_BC_0",
        // "hasFT0_CorrectedValid",
        "PV_FT0_diff_cut",
        // "PV_FT0_diff_cut_TIGHT",
        "ALL_CUTS",
        // "ALL_CUTS_Tighter",
        "NoTFborder_FoundBCwithDiff0",
        "NoTFborder_NoFoundBCwithDiff0",
        "NoTFborder_DiffBnBcAtLeast10",
        // "NoTFborder_DiffBnBcAtLeast20",
        "NoTFborder_DiffBnBcAtLeast40",
        "antiITSROFcut",
        "CutEventsByMultPVvsV0A",
        "antiCutEventsByMultPVvsV0A",
        "ALL_CUTS_CutEventsByMultPVvsV0A",
        "Handmade_ITSROFcut",
        "Handmade_ITSROFcut_OnFoundBC",
        "antiHandmade_ITSROFcut",
        "ALL_CUTS_Handmade_ITSROFcut",
        "ALL_CUTS_Handmade_ITSROFcut_OnFoundBC",
        "ALL_CUTS_Handmade_ITSROFcut_anotherCutOnMultV0A",
        "ALL_CUTS_Handmade_ITSROFcut_kTVXinTRD",
        "Handmade_ITSROF_and_TF_cuts",

        "isITSonlyVertex",
        "antiIsITSonlyVertex",

        "kTVXinTRD",

        "ALL_CUTS_Handmade_ITSROFcut_Orbit_1_part",
        "ALL_CUTS_Handmade_ITSROFcut_Orbit_2_part",
        "ALL_CUTS_Handmade_ITSROFcut_Orbit_3_part",
        "ALL_CUTS_Handmade_ITSROFcut_Orbit_4_part",
      };

    // naming of the 1D pT plots
    string strPtSelFolderNames[] =
      {
        "allBC",
        "TFborder",
        "noTFborder",
        "ROFborder",
        "noROFborder",
        "noTFandROFborder",
      };

    string strPtSelNames[] =
      {
        "pt",
        "pt_w_nTPCcls",
        "pos_pt",
        "pos_pt_w_nTPCcls",
        "neg_pt",
        "neg_pt_w_nTPCcls",
      };

    funcCutEventsByMultPVvsV0A = new TF1("funcCutEventsByMultPVvsV0A", "[0]*x+[1]", 0, 200000);
    funcCutEventsByMultPVvsV0A->SetParameters(3200. / 160000, -240);

    funcCutEventsByMultPVvsT0C = new TF1("funcCutEventsByMultPVvsT0C", "[0]*x+[1]", 0, 10000);
    funcCutEventsByMultPVvsT0C->SetParameters(8. / 1600, -0.2);

    funcCutEventsByMultPVvsV0A_anotherTrend = new TF1("funcCutEventsByMultPVvsV0A_anotherTrend", "[0]*x+[1]", 0, 200000);
    funcCutEventsByMultPVvsV0A_anotherTrend->SetParameters(3400. / 160000, -410);

    AxisSpec axisNcontrib{nBinsContrib, -0.5, nMaxContrib - 0.5, "n vertex contributors"};
    AxisSpec axisNtracks{nBinsTracks, -0.5, nMaxTracks - 0.5, "n tracks"};
    AxisSpec axisMultAllTr{200, -0.5, 5 * nMaxTracks - 0.5, "multiplicity"};

    AxisSpec axisCollIndex{nBinsCollIndex, -0.5, nMaxCollIndex - 0.5, "collision index"};
    AxisSpec axisMultFw{nBinsMultFwd, 0., static_cast<float>(nMaxMultFwd), "mult Fwd"}; //{1000, 0, 200000, "mult"};

    AxisSpec axisCent{100, 0.f, 100.f, "centrality"};
    AxisSpec axisCentBins{{0, 5., 10., 20., 30., 40., 50., 60., 70., 80.}, "centrality percentile"};

    AxisSpec axisCollTime{1000, -50, 50, "CollTime"};
    AxisSpec axisCollTimeRes{2000, -20, 20, "CollTimeRes"};
    AxisSpec axisBC{3601, -0.5, 3600.5, "bc"};

    AxisSpec axisZvert{120, -30.f, 30.f, "Vtx z (cm)"};
    AxisSpec axisDFQA{10001, -0.5, 10000.5, "collision index"};

    // ### per-DF QA:
    for (int i = 0; i < nHistQAplotsDF; i++)
      fV_h1D_Orbit_vs_CollIndex.push_back(histosDF.add<TH1>(TString::Format("h1D_Orbit_vs_CollIndex_DF%d", i).Data(), TString::Format("h1D_Orbit_vs_CollIndex_DF%d;collision index;orbit", i).Data(), kTH1F, {axisDFQA}));

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
    int myMaxOrbitsPerTF = 4 * 32; // 128 - in 2022, 32 - in 2023
    AxisSpec axisBCinTF{myMaxOrbitsPerTF * 4000 + 1, -0.5, myMaxOrbitsPerTF * 4000 + 0.5, "bc"};
    AxisSpec axisOrbit{100, -0.5, 9999.5, "orbit"};

    histosEvent.add("hRunNumber", "hRunNumber", kTH1D, {{6000, 534000.5, 540000.5, "hRunNumber"}});
    histosEvent.add("hMF", "hMF", kTH1D, {{100, -1, 1, "hMF"}});

    histosEvent.add("vtxCutsBef", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});
    histosEvent.add("vtxCutsAft", "Vtx distribution; Vtx z [cm]; Counts", kTH1F, {axisZvert});

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
    histosEvent.add("h2D_diffFoundBC_vs_BC_inTF", "h2D_diffFoundBC_vs_BC_inTF", kTH2D, {axisBC, {201, -100.5, 100.5, "foundBC-BC"}});

    histosEvent.add("hOrbitStartFromCollIndexZeroAft", "hOrbitStartFromCollIndexZeroAft", kTH1D, {axisOrbit});
    histosEvent.add("h2D_Orbit_vs_CollIndex_Aft", "h2D_Orbit_vs_CollIndex_Aft", kTH2D, {axisCollIndex, axisOrbit});

    histosEvent.add("hNtrackshGlobalAft", "hNtrackshGlobalAft", kTH1D, {axisNtracks});
    histosEvent.add("hNtrackshGlobalAft_AfterTimeFrameCut", "hNtrackshGlobalAft_AfterTimeFrameCut", kTH1D, {axisNtracks});

    // bc "correlations"
    histosEvent.add("hBC_DIFF_to_previous", "hBC_DIFF_to_previous", kTH1D, {{3564 * 2 + 1, -3564.5, 3564.5, "diff BC"}});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC", "hBC_DIFF_to_previous_FOUND_BC", kTH1D, {{3564 * 2 + 1, -3564.5, 3564.5, "diff BC"}});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC", "hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC", kTH1D, {{3564 * 2 + 1, -3564.5, 3564.5, "diff BC"}});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut", "hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut", kTH1D, {{3564 * 2 + 1, -3564.5, 3564.5, "diff BC"}});

    histosEvent.add("hBC_DIFF_to_previous_vZvZ_2D", "hBC_DIFF_to_previous_vZvZ_2D", kTH2F, {axisZvert, axisZvert});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_vZvZ_2D_diff0", "hBC_DIFF_to_previous_FOUND_BC_vZvZ_2D_diff0", kTH2F, {axisZvert, axisZvert});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC_vZvZ_2D_diff0", "hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC_vZvZ_2D_diff0", kTH2F, {axisZvert, axisZvert});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_diff0", "hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_diff0", kTH2F, {axisZvert, axisZvert});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_mult_mult_2D_diff0", "hBC_DIFF_to_previous_FOUND_BC_mult_mult_2D_diff0", kTH2F, {axisNtracks, axisNtracks});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_globTrkContrib_globTrkContrib_2D_diff0", "hBC_DIFF_to_previous_FOUND_BC_globTrkContrib_globTrkContrib_2D_diff0", kTH2F, {axisNtracks, axisNtracks});
    histosEvent.add("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_NoBCdiff0", "hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_NoBCdiff0", kTH2F, {axisZvert, axisZvert});

    histosEvent.add("hBCinTF_DIFF_to_previous", "hBCinTF_DIFF_to_previous", kTH1D, {{3564 * 4 + 1, -3564.5, 2 * 3564 + 3564.5, "diff BC"}});
    histosEvent.add("hBCinTF_DIFF_to_previous_FOUND_BC", "hBCinTF_DIFF_to_previous_FOUND_BC", kTH1D, {{3564 * 4 + 1, -3564.5, 2 * 3564 + 3564.5, "diff BC"}});

    // hist vs bcInTF, Feb 2, 2024
    if (flagIncludeQAHistBcInTF) {
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
    }

    // nTracks vs BC
    if (flagIncludeQAHistBcVsNtracks) {
      histosEvent.add("h2D_nTracksBeforeCuts_vs_BC", "h2D_nTracksBeforeCuts_vs_BC", kTH2F, {axisBC, axisNcontrib});
      histosEvent.add("h2D_nTracksAfterEtaTPCcuts_vs_BC", "h2D_nTracksAfterEtaTPCcuts_vs_BC", kTH2F, {axisBC, axisNcontrib});
      histosEvent.add("h2D_nTracksITSonly_vs_BC", "h2D_nTracksITSonly_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithITS_vs_BC", "h2D_nTracksWithITS_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithITS7hits_vs_BC", "h2D_nTracksWithITS7hits_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithITSandTPC_vs_BC", "h2D_nTracksWithITSandTPC_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithTRD_vs_BC", "h2D_nTracksWithTRD_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithTOF_vs_BC", "h2D_nTracksWithTOF_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksWithTRDorTOF_vs_BC", "h2D_nTracksWithTRDorTOF_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksGlobal_vs_BC", "h2D_nTracksGlobal_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksGlobalWithITS7hits_vs_BC", "h2D_nTracksGlobalWithITS7hits_vs_BC", kTH2F, {axisBC, axisNtracks});
      histosEvent.add("h2D_nTracksGlobalWithTRDorTOF_vs_BC", "h2D_nTracksGlobalWithTRDorTOF_vs_BC", kTH2F, {axisBC, axisNtracks});

      // mean pT vs BC
      AxisSpec axisMeanPt{10, 0., 2.0, "<p_{T}>"};
      histosEvent.add("h1D_vContributors_1_meanPt_vs_BC", "h1D_vContributors_1_meanPt_vs_BC;bc;entries", kTH2F, {axisBC, axisMeanPt});
      histosEvent.add("h1D_vContributors_2_meanPt_vs_BC", "h1D_vContributors_2_meanPt_vs_BC;bc;entries", kTH2F, {axisBC, axisMeanPt});
      histosEvent.add("h1D_vContributors_3_meanPt_vs_BC", "h1D_vContributors_3_meanPt_vs_BC;bc;entries", kTH2F, {axisBC, axisMeanPt});
      histosEvent.add("h1D_vContributors_4_meanPt_vs_BC", "h1D_vContributors_4_meanPt_vs_BC;bc;entries", kTH2F, {axisBC, axisMeanPt});
      histosEvent.add("h1D_vContributors_5_meanPt_vs_BC", "h1D_vContributors_5_meanPt_vs_BC;bc;entries", kTH2F, {axisBC, axisMeanPt});
    }

    // AxisSpec axisGlobalTracks{flagPbPb ? 2501 : 81, -0.5, flagPbPb ? 2500.5 : 80.5, "n global tracks"};
    histosEvent.add("h1D_EventCounter_vs_BC", "h1D_EventCounter_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks0Global_vs_BC", "h1D_nTracks0Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks1Global_vs_BC", "h1D_nTracks1Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks2Global_vs_BC", "h1D_nTracks2Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks3Global_vs_BC", "h1D_nTracks3Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks4Global_vs_BC", "h1D_nTracks4Global_vs_BC;bc;entries", kTH1D, {axisBC});
    histosEvent.add("h1D_nTracks5Global_vs_BC", "h1D_nTracks5Global_vs_BC;bc;entries", kTH1D, {axisBC});

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
    histosEvent.add("hBC_Aft_rejectedByCutOnMultPVvsV0A", "hBC_Aft_rejectedByCutOnMultPVvsV0A", kTH1D, {axisBC});
    histosEvent.add("hBC_Aft_Handmade_ITSROFcut", "hBC_Aft_Handmade_ITSROFcut", kTH1D, {axisBC});

    histosEvent.add("hBC_Aft_rejectedByCutOnMultPVvsT0C", "hBC_Aft_rejectedByCutOnMultPVvsT0C", kTH1D, {axisBC});
    histosEvent.add("hBC_Aft_rejectedByCutOnMultPVvsT0C_after_Handmade_ITSROF_cut", "hBC_Aft_rejectedByCutOnMultPVvsT0C_after_Handmade_ITSROF_cut", kTH1D, {axisBC});
    histosEvent.add("hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV", "hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV", kTH1D, {axisBC});
    histosEvent.add("hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV_ANTI", "hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV_ANTI", kTH1D, {axisBC});

    histosEvent.add("hOccupancyInTimeRange", "hOccupancyInTimeRange", kTH1D, {{400, 0, 20000, "occupancy (nITStracks nCls>=5)"}});
    histosEvent.add("hOccupancyInTimeRangeIsAvailable", "hOccupancyInTimeRangeIsAvailable", kTH1D, {{203, -2.5, 200.5, "occupancy zoomed (nITStracks nCls>=5)"}});

    AxisSpec axisCentrForOccupQA{20, 0, 100, "centrality, %"};
    histosEvent.add("hCentralityVsOccupancy", "hCentralityVsOccupancy", kTH2D, {{400, 0, 20000, "occupancy zoomed (nITStracks nCls>=5)"}, axisCentrForOccupQA});

    // ##### tracks for All collisions
    if (flagIncludeQATracksInFoundBC) {
      histosEventTracksFoundBC.add("hFoundBC_nAllTracks", "hFoundBC_nAllTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTracksPV", "hFoundBC_nTracksPV", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nGlobalTracks", "hFoundBC_nGlobalTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nITStracks", "hFoundBC_nITStracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTPCtracks", "hFoundBC_nTPCtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTOFtracks", "hFoundBC_nTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTRDtracks", "hFoundBC_nTRDtracks", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_nTRDTOFtracks", "hFoundBC_nTRDTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nITSTPCTRDtracks", "hFoundBC_nITSTPCTRDtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nITSTPCTOFtracks", "hFoundBC_nITSTPCTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nITSTPCTRDTOFtracks", "hFoundBC_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_nTPCtracks_HandmadeITSROFcut", "hFoundBC_nTPCtracks_HandmadeITSROFcut", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTOFtracks_HandmadeITSROFcut", "hFoundBC_nTOFtracks_HandmadeITSROFcut", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_nTRDTOFtracks_HandmadeITSROFcut", "hFoundBC_nTRDTOFtracks_HandmadeITSROFcut", kTH1D, {axisBC});

      // ##### tracks for TVX collisions
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nAllTracks", "hFoundBC_kTVX_nAllTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTracksPV", "hFoundBC_kTVX_nTracksPV", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nGlobalTracks", "hFoundBC_kTVX_nGlobalTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITStracks", "hFoundBC_kTVX_nITStracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTPCtracks", "hFoundBC_kTVX_nTPCtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCtracks", "hFoundBC_kTVX_nITSTPCtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTOFtracks", "hFoundBC_kTVX_nTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTRDtracks", "hFoundBC_kTVX_nTRDtracks", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTRDTOFtracks", "hFoundBC_kTVX_nTRDTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks", "hFoundBC_kTVX_nITSTPCTRDtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTOFtracks", "hFoundBC_kTVX_nITSTPCTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDTOFtracks", "hFoundBC_kTVX_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTPCtracks_HandmadeITSROFcut", "hFoundBC_kTVX_nTPCtracks_HandmadeITSROFcut", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTOFtracks_HandmadeITSROFcut", "hFoundBC_kTVX_nTOFtracks_HandmadeITSROFcut", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nTRDTOFtracks_HandmadeITSROFcut", "hFoundBC_kTVX_nTRDTOFtracks_HandmadeITSROFcut", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSlayers_for_GlobalTracks", "hFoundBC_kTVX_nITSlayers_for_GlobalTracks", kTH1D, {axisBC});
      AxisSpec axisPtBinsForROFstudy{{0.1, 0.2, 0.3, 0.4, 0.5, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0}, "p_{T}"};
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt", "hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt", kTH2D, {axisBC, axisPtBinsForROFstudy});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCtracks_vs_pt", "hFoundBC_kTVX_nITSTPCtracks_vs_pt", kTH2D, {axisBC, axisPtBinsForROFstudy});

      // ##### tracks for TVXinTRD collisions
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nAllTracks", "hFoundBC_kTVXinTRD_nAllTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nTracksPV", "hFoundBC_kTVXinTRD_nTracksPV", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nGlobalTracks", "hFoundBC_kTVXinTRD_nGlobalTracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nITStracks", "hFoundBC_kTVXinTRD_nITStracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nTPCtracks", "hFoundBC_kTVXinTRD_nTPCtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nTOFtracks", "hFoundBC_kTVXinTRD_nTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nTRDtracks", "hFoundBC_kTVXinTRD_nTRDtracks", kTH1D, {axisBC});

      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nTRDTOFtracks", "hFoundBC_kTVXinTRD_nTRDTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nITSTPCTRDtracks", "hFoundBC_kTVXinTRD_nITSTPCTRDtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nITSTPCTOFtracks", "hFoundBC_kTVXinTRD_nITSTPCTOFtracks", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks", "hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

      // to study ITS ROF dips vs multiplicity
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000", kTH1D, {axisBC});
      histosEventTracksFoundBC.add("hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000", "hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000", kTH1D, {axisBC});

      // to study ITS ROF dips vs pt
      if (flagIncludeQAROFinPtEtaBins) {
        AxisSpec axisITSlayers{7, -0.5, 7 - 0.5, "layer id"};
        histosEventTracksFoundBC.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});

        // pions
        histosEventTracksFoundBC.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pi_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pi_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_pi_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});

        // protons
        histosEventTracksFoundBC.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02_04", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02_04", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta04_06", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta04_06", kTH2D, {axisBC, axisPtBinsForROFstudy});
        histosEventTracksFoundBC.add("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta06_08", "hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta06_08", kTH2D, {axisBC, axisPtBinsForROFstudy});

        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_p_pt02_05", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_p_pt05_10", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
        histosEventTracksFoundBC.add("hITSlayerCounts_vs_BC_p_pt10_50", "hITSlayerCounts_vs_BC", kTH2D, {axisBC, axisITSlayers});
      }

      // bcInTF
      AxisSpec axisBCinTF32orbits{32 * 3564, -0.5, 32 * 3564, "bcInTF"};
      histosEventTracksFoundBC.add("hBCinTF_kTVX_nTracksPV", "hBCinTF_kTVX_nTracksPV;bc in TF; n tracks", kTH1D, {axisBCinTF32orbits});
      histosEventTracksFoundBC.add("hBCinTF_kTVX_nITSTPCTRDtracks", "hBCinTF_kTVX_nITSTPCTRDtracks;bc in TF; n tracks", kTH1D, {axisBCinTF32orbits});

    } // end of if flagIncludeQATracksInFoundBC

    // #### add QA histos via loop
    int nCutsQA = sizeof(strCutNames) / sizeof(*strCutNames);
    for (int i = 0; i < nCutsQA; i++) {
      // add this cut name and id to map
      string cutName = strCutNames[i].c_str();
      mMyCuts[cutName] = i;

      // now add 1D histograms:
      vector<std::shared_ptr<TH1>> v1D;

      // ADD_QA_HIST_2D("multAllTr_vs_Cent", "multiplicity vs centrality T0C", axisCent, axisMultAllTr);
      ADD_QA_HIST_1D("nGlobalAcceptedTracks", "nGlobalAcceptedTracks", axisNtracks);
      ADD_QA_HIST_1D("bcWithCuts", "bcWithCuts", axisBC);

      // add vector with pointers to histos for a given cut to an "external" vector
      pHist1D.push_back(v1D);

      // now add 2D histograms:
      vector<std::shared_ptr<TH2>> v2D;

      // ADD_QA_HIST_2D("multAllTr_vs_Cent", "multiplicity vs centrality T0C", axisCent, axisMultAllTr);
      ADD_QA_HIST_2D("multAllTr_vs_multT0C", "multiplicity vs multiplicity T0C", axisMultFw, axisMultAllTr);
      ADD_QA_HIST_2D("multAllTr_vs_multV0A", "multiplicity vs multiplicity V0A", axisMultFw, axisMultAllTr);
      ADD_QA_HIST_2D("multAllTr_vs_multT0A", "multiplicity vs multiplicity T0A", axisMultFw, axisMultAllTr);
      ADD_QA_HIST_2D("multAllTr_vs_multTrkPV", "multiplicity vs multiplicity PV", axisNtracks, axisMultAllTr);
      // ADD_QA_HIST_2D("multGlobalTr_vs_Cent", "multiplicity vs centrality T0C", axisCent, axisMult);
      ADD_QA_HIST_2D("multGlobalTr_vs_multT0C", "multiplicity vs multiplicity T0C", axisMultFw, axisNtracks);
      ADD_QA_HIST_2D("multGlobalTr_vs_multV0A", "multiplicity vs multiplicity V0A", axisMultFw, axisNtracks);
      ADD_QA_HIST_2D("multGlobalTr_vs_multT0A", "multiplicity vs multiplicity T0A", axisMultFw, axisNtracks);
      ADD_QA_HIST_2D("multGlobalTr_vs_multTrkPV", "multiplicity vs multiplicity PV", axisNtracks, axisNtracks);
      // ADD_QA_HIST_2D("multTrkPV_vs_Cent", "multiplicity PV vs centrality T0C", axisCent, axisMult);
      ADD_QA_HIST_2D("multTrkPV_vs_multT0C", "multiplicity PV vs multiplicity T0C", axisMultFw, axisNtracks);
      ADD_QA_HIST_2D("multTrkPV_vs_multV0A", "multiplicity PV vs multiplicity V0A", axisMultFw, axisNtracks);
      ADD_QA_HIST_2D("multTrkPV_vs_multT0A", "multiplicity PV vs multiplicity T0A", axisMultFw, axisNtracks);
      // ADD_QA_HIST_2D("multV0A_vs_Cent", "multiplicity V0A vs centrality T0C", axisCent, axisMultFw);
      ADD_QA_HIST_2D("multT0C_vs_multT0A", "multiplicity T0C vs multiplicity T0A", axisMultFw, axisMultFw);
      ADD_QA_HIST_2D("multV0A_vs_multT0A", "multiplicity V0A vs multiplicity T0A", axisMultFw, axisMultFw);
      ADD_QA_HIST_2D("multV0A_vs_multT0C", "multiplicity V0A vs multiplicity T0C", axisMultFw, axisMultFw);
      // ADD_QA_HIST_2D("multT0C_vs_Cent", "multiplicity T0C vs centrality T0C", axisCent, axisMultFw);

      // add vector with pointers to histos for a given cut to an "external" vector
      pHist2D.push_back(v2D);
    }

    AxisSpec axisCounterAllVsTF{10, -0.5, 9.5, "cut"};
    histosEventCounters.add("hNtracksAll_vs_variousCuts", "hNtracksAll_vs_variousCuts", kTH1D, {axisCounterAllVsTF});
    histosEventCounters.add("hNtracksGlobal_vs_variousCuts", "hNtracksGlobal_vs_variousCuts", kTH1D, {axisCounterAllVsTF});
    histosEventCounters.add("hNtotalCollisions_vs_variousCuts", "hNtotalCollisions_vs_variousCuts", kTH1D, {axisCounterAllVsTF});

    AxisSpec axisVtxChi2{500, 0., 100.f, "chi2"};
    histosEvent.add("vtxChi2Aft", "vtxChi2Aft; chi2; Counts", kTH1F, {axisVtxChi2});
    histosEvent.add("vtxChi2_vertNcontr3_TRDorTOF", "vtxChi2_vertNcontr3_TRDorTOF; chi2; Counts", kTH1F, {axisVtxChi2});

    // #### FT0 time
    const AxisSpec axisTimeFT0{400, -2., 2., "collision time (ns)"};
    const AxisSpec axis2DTimeFT0{100, -2., 2., "collision time (ns)"};
    const AxisSpec axisColTimeResFT0{!flagPbPb ? 500 : 2000, -0.5, 0.5, "(T0A - T0C)/2 (ns)"};
    const AxisSpec axisVertexFT0{300, -30., 30.};
    const AxisSpec axisVertexFT0diff{1200, -30., 30.};

    histosFT0.add("hT0A", "T0A;T0A time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0C", "T0C;T0C time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0AC", "T0AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0res", "FT0 resolution", kTH1F, {axisColTimeResFT0});
    // histos.add("hColTime", "", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_sum_AC", "hT0_sum_AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_diff_AC", "hT0_diff_AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_diff_vs_sum_AC", ";T0A-T0C;T0A+T0C;counts", kTH2F, {axis2DTimeFT0, axis2DTimeFT0});

    histosFT0.add("hT0timeA_uncorr", "T0A;T0A time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0timeC_uncorr", "T0C;T0C time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_sum_AC_uncorr", "hT0_sum_AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_diff_AC_uncorr", "hT0_diff_AC;T0AC time (ns);counts", kTH1F, {axisTimeFT0});
    histosFT0.add("hT0_diff_vs_sum_AC_uncorr", ";T0A-T0C;T0A+T0C;counts", kTH2F, {axis2DTimeFT0, axis2DTimeFT0});

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
    histosFT0.add("hVertex_T0_PV_after_HandmadeITSROFcut", "PV vs. FT0V after handmade ITSROF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_HandmadeITSROFcut_ANTI", "PV vs. FT0V after handmade ITSROF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_TFcut", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_ITSROF_and_TFcut", "PV vs. FT0V after ITSROF and TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2goodPVcontribTPC", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_4globalPVcontrib", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib_ITS7hits", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_2globalPVcontrib_TOForTRD", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});

    histosFT0.add("hVertex_T0_PV_after_NoBCdiff0", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_NoBCdiff0_ANTI", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_ITSonlyVertex", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});
    histosFT0.add("hVertex_T0_PV_after_ITSonlyVertex_ANTI", "PV vs. FT0V after TF cut;FT0 vertex (cm);primary vertex (cm)", kTH2F, {axisVertexFT0, axisVertexFT0});

    AxisSpec axisVertexFT0diffFor2D{300, -15., 15., "FT0 vertex -  PV (cm)"};
    AxisSpec axisNtracksForV02D{200, -0.5, 200 - 0.5, "multPV"};
    histosFT0.add("hT0_2D_multPV_vs_vertexDiff", "multPVFT0V - PV diff;counts", kTH2F, {axisVertexFT0diffFor2D, axisNtracksForV02D});

    // ### track-wise:
    AxisSpec axisPt{nBinsPt, 0, 20, "p_{T}"};
    AxisSpec axisEta{nBinsEta, -1.5, +1.5, "#eta"};
    AxisSpec axisPhi{nBinsPhi, 0, TMath::TwoPi(), "#varphi"};
    AxisSpec axisPhiSpecMod9{nBinsPhi, -2 * TMath::TwoPi() / 9, 2 * TMath::TwoPi() / 9, "#varphi"};
    histosTracks.add("eta", "eta", kTH1D, {axisEta});
    histosTracks.add("etaAfter08cut", "etaAfter08cut", kTH1D, {axisEta});

    // ### detailed pt study
    AxisSpec axisLogPt{100, 0.05, 40, "p_{T}"};
    axisLogPt.makeLogarithmic();

    // #### add pt QA histos via loop
    int nFolderPt = sizeof(strPtSelFolderNames) / sizeof(*strPtSelFolderNames);
    int nCutsPtQA = sizeof(strPtSelNames) / sizeof(*strPtSelNames);
    // int counterPtCuts = 0;
    for (int i = 0; i < nFolderPt; i++) {
      for (int j = 0; j < nCutsPtQA; j++) {
        // add this cut name and id to map
        string cutName = (strPtSelFolderNames[i] + "/" + strPtSelNames[j]).c_str();

        // now add 1D histograms:
        ADD_PT_HIST_1D("All", cutName, "", axisLogPt);
        ADD_PT_HIST_1D("Global", cutName, "", axisLogPt);

        ADD_PT_HIST_1D("ITS7hits", cutName, "", axisLogPt);
        ADD_PT_HIST_1D("ITS7hits_TPC80cl", cutName, "", axisLogPt);

        ADD_PT_HIST_1D("ITS4567", cutName, "", axisLogPt);
        ADD_PT_HIST_1D("ITS4567_TPC80cl", cutName, "", axisLogPt);

        ADD_PT_HIST_1D("ITS567", cutName, "", axisLogPt);
        ADD_PT_HIST_1D("ITS567_TPC80cl", cutName, "", axisLogPt);

        ADD_PT_HIST_1D("ITS67", cutName, "", axisLogPt);
        ADD_PT_HIST_1D("ITS67_TPC80cl", cutName, "", axisLogPt);
      }
    }
    // ### end of detailed pt study

    histosTracks.add("phiHistogram", "phiHistogram", kTH1D, {axisPhi});

    histosTracks.add("pidCombSigma", "pidCombSigma", kTH1D, {{200, 0, 10, "pidCombSigma"}});

    histosTracks.add("hTpcNClsCrossedRows", "hTpcNClsCrossedRows", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRows"}});
    histosTracks.add("hTpcNClsCrossedRows_Global_beforeCrossedRowsCut", "hTpcNClsCrossedRows_Global_beforeCrossedRowsCut", kTH1D, {{170, -0.5, 169.5, "hTpcNClsCrossedRows_Global_beforeCrossedRowsCut"}});
    histosTracks.add("hTpcNClsCrossedRowsITS7hits", "hTpcNClsCrossedRowsITS7hits", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRowsITS7hits"}});
    histosTracks.add("hNumITSclusters", "hNumITSclusters", kTH1D, {{10, -0.5, 9.5, "NumITSclusters"}});
    histosTracks.add("hNumITSclusters_Global_beforeCrossedRowsCut", "hNumITSclusters_Global_beforeCrossedRowsCut", kTH1D, {{10, -0.5, 9.5, "hNumITSclusters_Global_beforeCrossedRowsCut"}});
    histosTracks.add("hNumITSclusters_ITSonlyVert", "hNumITSclusters_ITSonlyVert", kTH1D, {{10, -0.5, 9.5, "NumITSclusters"}});
    histosTracks.add("hDcaXY", "hDcaXY", kTH1D, {{800, -4, 4, "DcaXY"}});
    histosTracks.add("hDcaXY_Global_beforeCrossedRowsCut", "hDcaXY_Global_beforeCrossedRowsCut", kTH1D, {{800, -4, 4, "DcaXY"}});
    histosTracks.add("hDcaZ", "hDcaZ", kTH1D, {{800, -4, 4, "DcaZ"}});
    histosTracks.add("hDcaZ_Global_beforeCrossedRowsCut", "hDcaZ_Global_beforeCrossedRowsCut", kTH1D, {{800, -4, 4, "DcaZ"}});
    histosTracks.add("hTrackLength", "hTrackLength", kTH1D, {{400, 0, 2000, "TrackLength"}});

    AxisSpec axisTrackChi2{200, 0., 10.f, "chi2"};
    // AxisSpec axisPtBinsForPhiGaps{{0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0}, "p_{T}"};
    AxisSpec axisPtBinsForPhiGaps{{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0}, "p_{T}"};
    histosTracks.add("hChi2TPCperDOF", "hChi2TPCperDOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hChi2ITSperDOF", "hChi2ITSperDOF", kTH1D, {axisTrackChi2});
    histosTracks.add("hTrackTime", "hTrackTime", kTH1D, {axisCollTime});
    histosTracks.add("hTrackTimeRes", "hTrackTimeRes", kTH1D, {axisCollTime});

    histosTracks.add("etaAftCuts", "etaAftCuts", kTH1D, {axisEta});
    histosTracks.add("etaAftCutsGlobal", "etaAftCutsGlobal", kTH1D, {axisEta});
    histosTracks.add("etaITS7hits_AftCuts", "etaITS7hits_AftCuts", kTH1D, {axisEta});
    histosTracks.add("ptHistogram_AftCuts", "ptHistogram_AftCuts", kTH1D, {axisPt});
    histosTracks.add("phiHistogram_AftCuts", "phiHistogram_AftCuts", kTH1D, {axisPhi});

    histosTracks.add("hTpcNClsCrossedRows_AftCuts", "hTpcNClsCrossedRows_AftCuts", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRows"}});
    histosTracks.add("hTpcNClsCrossedRowsITS7hits_AftCuts", "hTpcNClsCrossedRowsITS7hits_AftCuts", kTH1D, {{170, -0.5, 169.5, "TpcNClsCrossedRowsITS7hits_AftCuts"}});
    histosTracks.add("hNumITSclusters_AftCuts", "hNumITSclusters_AftCuts", kTH1D, {{10, -0.5, 9.5, "NumITSclusters_AftCuts"}});
    histosTracks.add("hNumITSclusters_AftCuts_Global", "hNumITSclusters_AftCuts_Global", kTH1D, {{10, -0.5, 9.5, "NumITSclusters_AftCuts"}});
    histosTracks.add("hDcaXY_AftCuts", "hDcaXY_AftCuts", kTH1D, {{800, -4, 4, "DcaXY_AftCuts"}});
    histosTracks.add("hDcaZ_AftCuts", "hDcaZ_AftCuts", kTH1D, {{800, -4, 4, "DcaZ_AftCuts"}});
    histosTracks.add("hTrackLength_AftCuts", "hTrackLength_AftCuts", kTH1D, {{400, 0, 2000, "TrackLength_AftCuts"}});
    histosTracks.add("hChi2perDOF_AftCuts", "hChi2perDOF_AftCuts", kTH1D, {{200, 0, 20, "Chi2perDOF_AftCuts"}});

    histosTracks.add("hTpcNClsFound", "hTpcNClsFound", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound"}});
    histosTracks.add("hTpcNClsFound_Global_beforeCrossedRowsCut", "hTpcNClsFound_Global_beforeCrossedRowsCut", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound_Global_beforeCrossedRowsCut"}});
    histosTracks.add("hTpcNClsFound_AftCuts", "hTpcNClsFound_AftCuts", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound_AftCuts"}});
    histosTracks.add("hTpcNClsFound_ITS7hits_AftCuts", "hTpcNClsFound_ITS7hits_AftCuts", kTH1D, {{170, -0.5, 169.5, "hTpcNClsFound_ITS7hits_AftCuts"}});

    histosTracks.add("hNumITScls_vs_TPCcls_AftCuts_Global", "hNumITScls_vs_TPCcls_AftCuts_Global", kTH2D, {{8, -0.5, 7.5, "NumITSclusters"}, {170, -0.5, 169.5, "TpcNCls"}});
    histosTracks.add("hNumITScls_vs_TPCcrossedRows_AftCuts_Global", "hNumITScls_vs_TPCcrossedRows_AftCuts_Global", kTH2D, {{8, -0.5, 7.5, "NumITSclusters"}, {170, -0.5, 169.5, "TpcNClsCrossedRows"}});
    histosTracks.add("hNumITScls_vs_TPCcrossedRows_Global_beforeCrossedRowsCut", "hNumITScls_vs_TPCcrossedRows_Global_beforeCrossedRowsCut", kTH2D, {{8, -0.5, 7.5, "NumITSclusters"}, {170, -0.5, 169.5, "TpcNClsCrossedRows"}});

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

    // ##### v0s
    if (flagIncludeQAHistK0S) {
      AxisSpec K0ShortMassAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
      histosK0S.add("hK0Sradius", "hK0Sradius", kTH1D, {{800, 0, 100, "hK0Sradius"}});
      histosK0S.add("hMassK0Short", "hMassK0Short", {HistType::kTH1F, {K0ShortMassAxis}});
      histosK0S.add("hMassK0ShortAfterSelection", "hMassK0ShortAfterSelection", {HistType::kTH1F, {K0ShortMassAxis}});

      AxisSpec axisK0SptBins{{0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0}, "p_{T}"};
      histosK0S.add("hMassK0ShortAfterSelectionVsPt", "hMassK0ShortAfterSelectionVsPt", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtEta01_08", "hMassK0ShortAfterSelectionVsPtEta01_08", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtEta08_01", "hMassK0ShortAfterSelectionVsPtEta08_01", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigma", "hMassK0ShortAfterSelectionVsPtPIDnSigma", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4_andROFcut", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNcontribAbove4_andROFcut", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtPIDnSigmaNITS7hitsContribAbove2", "hMassK0ShortAfterSelectionVsPtPIDnSigmaNITS7hitsContribAbove2", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0S_nVertexContributorsWithTRDorTOF_Above2", "hMassK0S_nVertexContributorsWithTRDorTOF_Above2", kTH2D, {axisK0SptBins, K0ShortMassAxis});
      histosK0S.add("hMassK0S_nVertexContributorsWithTRDorTOF_ITS7hits_Above2", "hMassK0S_nVertexContributorsWithTRDorTOF_ITS7hits_Above2", kTH2D, {axisK0SptBins, K0ShortMassAxis});

      AxisSpec axisCentrality{5, 0, 100, "centrality, %"};

      histosK0S.add("hMassK0ShortAfterSelectionVsPt_AllOccupancies", "hMassK0ShortAfterSelectionVsPt_AllOccupancies", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_0_250", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_0_250", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_250_500", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_250_500", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_500_750", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_500_750", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_750_1500", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_750_1500", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_1500_3000", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_1500_3000", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_3000_4500", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_3000_4500", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_4500_6000", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_4500_6000", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});
      histosK0S.add("hMassK0ShortAfterSelectionVsPtAfterOccupCut_above_6000", "hMassK0ShortAfterSelectionVsPtAfterOccupCut_above_6000", kTH3D, {axisK0SptBins, K0ShortMassAxis, axisCentrality});

      histosK0S.add("nEventsVsCentrAllOccupancies", "nEventsVsCentrAllOccupancies", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_0_250", "nEventsVsCentrAfterOccupCut_0_250", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_250_500", "nEventsVsCentrAfterOccupCut_250_500", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_500_750", "nEventsVsCentrAfterOccupCut_500_750", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_750_1500", "nEventsVsCentrAfterOccupCut_750_1500", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_1500_3000", "nEventsVsCentrAfterOccupCut_1500_3000", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_3000_4500", "nEventsVsCentrAfterOccupCut_3000_4500", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_4500_6000", "nEventsVsCentrAfterOccupCut_4500_6000", kTH1D, {axisCentrality});
      histosK0S.add("nEventsVsCentrAfterOccupCut_above_6000", "nEventsVsCentrAfterOccupCut_above_6000", kTH1D, {axisCentrality});

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
    } // end of if flagIncludeQAHistK0S
  }

  using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::pidTPCPi>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::Mults, aod::FT0sCorrected>;
  // using Colls = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::FT0sCorrected>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>; //, aod::Run3MatchedToBCSparse>;

  void processRobustFluctuationObservables(
    Colls::iterator const& collision,
    aod::FT0s const&,
    BCsRun3 const& bcs,
    aod::Origins const& origins,
    soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection,
              aod::TrackSelectionExtension,
              aod::TracksDCA, aod::pidTPCEl, aod::pidTPCMu, aod::pidTPCPi,
              aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCDe, aod::pidTOFEl> const& tracks,
    aod::V0Datas const& v0s, DaughterTracks const&)
  {
    auto bc = collision.bc_as<BCsRun3>();

    uint32_t orbit = bc.globalBC() / o2::constants::lhc::LHCMaxBunches;
    if (collision.index() == 0)
      orbitAtCollIndexZero = orbit;

    histosEvent.fill(HIST("hCollIndexBef"), collision.index());
    histosEvent.fill(HIST("hCollTimeBef"), collision.collisionTime());
    histosEvent.fill(HIST("hCollTimeResBef"), collision.collisionTimeRes());
    histosEvent.fill(HIST("hNumContribBef"), collision.numContrib());
    // histosEvent.fill(HIST("hVertexZ"), collision.posZ() );
    histosEvent.fill(HIST("hNtracksBef"), tracks.size());

    auto collBC = bc.globalBC() % 3564;
    uint64_t globalFoundBC = 9999;
    // int64_t globalFoundBC = -1;
    if (collision.has_foundBC()) {
      auto bcFound = collision.foundBC_as<BCsRun3>(); // collision.foundBC();
      globalFoundBC = bcFound.globalBC() % 3564;
    }
    histosEvent.fill(HIST("hBC_Bef"), collBC);

    // cut ITS ROF by hand, with configurables
    bool flagBCisNotInHandmadeBoundariesITSROF = true;
    if (collBC > nITSROF_BC_cutWidth / 2 && collBC < (nBCsPerOrbit - nITSROF_BC_cutWidth / 2)) {
      for (int iROF = 0; iROF < nITSROF; iROF++) {
        int ROF_BC_center = nITSROF_BC_offset + iROF * nBCsPerOrbit / nITSROF;
        if (fabs(static_cast<int>(collBC) - ROF_BC_center) <= nITSROF_BC_cutWidth / 2) {
          flagBCisNotInHandmadeBoundariesITSROF = false;
          break;
        }
      }
    }

    // similar cut ITS ROF by hand - for FOUND BC, with configurables
    bool flagFoundBCisNotInHandmadeBoundariesITSROF = (globalFoundBC == 9999 ? false : true);
    if (globalFoundBC > nITSROF_BC_cutWidth / 2 && globalFoundBC < (nBCsPerOrbit - nITSROF_BC_cutWidth / 2)) {
      for (int iROF = 0; iROF < nITSROF; iROF++) {
        int ROF_BC_center = nITSROF_BC_offset + iROF * nBCsPerOrbit / nITSROF;
        if (fabs(static_cast<int>(globalFoundBC) - ROF_BC_center) <= nITSROF_BC_cutWidth / 2) {
          flagFoundBCisNotInHandmadeBoundariesITSROF = false;
          break;
        }
      }
    }

    double vZ = collision.posZ();

    // #### values for 2D histogram axes
    auto t0cCentr = collision.centFT0C();

    auto multV0A = collision.multFV0A();
    auto multT0A = collision.multFT0A();
    auto multT0C = collision.multFT0C();
    auto multNTracksPV = collision.multNTracksPV();

    Int_t multTrk = tracks.size();

    // ######## Data Frames QA
    int64_t myDF_ID = -1;
    uint64_t DF_ID_raw = 0;
    if (mapDF.size() < 10) {
      for (auto const& origin : origins) {
        uint64_t DF_ID = origin.dataframeID();
        DF_ID_raw = DF_ID;

        if (origin.globalIndex() == 0) // look only at the id of the first subDF in this DF
        {
          if (mapDF.find(DF_ID) == mapDF.end()) {
            // not found
            mapDF.insert({DF_ID, mapDF.size()});
          }
          myDF_ID = mapDF[DF_ID];
        }
        // cout << "DF globalIndex = " << origin.globalIndex() << ", ID = " << origin.dataframeID() << ", myDF_ID = " << myDF_ID << ", mapDF.size() = " << mapDF.size() << endl;
      }
    }

    if (myDF_ID >= 0 && myDF_ID < nHistQAplotsDF) {
      int diffOrbits = static_cast<int32_t>(orbit) - static_cast<int32_t>(orbitAtCollIndexZero);
      TString strDF = Form("DF_%d", static_cast<int>(DF_ID_raw));
      fV_h1D_Orbit_vs_CollIndex[myDF_ID]->Fill(collision.index(), diffOrbits);
      fV_h1D_Orbit_vs_CollIndex[myDF_ID]->SetTitle(strDF);
    }

    if (0)
      LOGF(info, "collision.globalIndex() = %d, index() = %d, has_foundBC() = %d, tracks.size() = %d", //, globalBC=%.1f",
           collision.globalIndex(), collision.index(),
           //  DF_ID, myDF_ID,
           (int)collision.has_foundBC(),
           tracks.size());

    // int nTracksBeforeCutsITSonly = 0;
    int nTracksBeforeCutsITSonlyPV = 0;

    int nTracksAll = 0;
    // int nTracksAfterEtaCuts = 0;
    // int nTracksAfterEtaCutsPV = 0;
    int nTracksAfterEtaCutsITSonly = 0;
    // int nTracksAfterEtaCutsITSonlyPV = 0;

    int nTracksAfterEtaTPCCuts = 0;
    int nTracksWithITS = 0, nTracksWithITSandTPC = 0;
    int nTracksWithTRD = 0, nTracksWithTOF = 0, nTracksWithTRDorTOF = 0;
    int nTracksWithITS7hits = 0;
    int nTracksGlobalAccepted = 0;
    int nTracksGlobalWithITS7hits = 0;
    int nTracksGlobalWithTRDorTOF = 0;

    int nTracksGlobalPVAccepted = 0;
    int nTracksGlobalPVwithITS7hits = 0;
    int nTracksGlobalPVwithTRDorTOF = 0;

    int counterPVcontributorsBeforeCuts = 0;
    int counterPVcontributorsAfterTPCcuts = 0;
    int counterPVcontributorsITS7hits = 0;

    int counterVertexContributorsWithTRDorTOF = 0;
    int counterVertexContributorsWithTRDorTOF_ITS7hits = 0;

    int counterPVcontributorsNoTOFandTRD = 0;
    double meanPtForPVContributorsNoTOFandTRD = 0;

    // ### track pre-loop
    for (auto& track : tracks) {
      nTracksAll++;
      if (track.isPVContributor())
        counterPVcontributorsBeforeCuts++;

      if (track.hasITS() && !track.hasTPC()) // Flag to check if track has a TPC match
      {
        // nTracksBeforeCutsITSonly++;
        if (track.isPVContributor())
          nTracksBeforeCutsITSonlyPV++;
      }

      if (fabs(track.eta()) > 0.8)
        continue;

      if (track.hasITS() && !track.hasTPC()) // Flag to check if track has a TPC match
      {
        nTracksAfterEtaCutsITSonly++;
        // if (track.isPVContributor())
        // nTracksAfterEtaCutsITSonlyPV++;
      }

      // nTracksAfterEtaCuts++;
      // if (track.isPVContributor())
      // nTracksAfterEtaCutsPV++;

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

    // flag for collisions with ITS-only vertex
    bool isITSonlyVertex = (nTracksBeforeCutsITSonlyPV == counterPVcontributorsBeforeCuts);

    histosEvent.fill(HIST("vtxCutsBef"), vZ);
    fillHistForThisCut("Bef", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    // routine for each new run
    int runNumber = bc.runNumber();
    histosEvent.fill(HIST("hRunNumber"), runNumber);
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

      // ##### code to find TF borders, Feb 1, 2024
      int64_t ts = bcs.iteratorAt(0).timestamp();
      // access orbit-reset timestamp
      auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", ts);
      int64_t tsOrbitReset = (*ctpx)[0]; // us

      std::map<std::string, std::string> metadata;
      metadata["runNumber"] = Form("%d", runNumber);
      auto grpecs = ccdb->getSpecific<o2::parameters::GRPECSObject>("GLO/Config/GRPECS", ts, metadata);
      nOrbitsPerTF = grpecs->getNHBFPerTF(); // assuming 1 orbit = 1 HBF
      tsSOR = grpecs->getTimeStart();        // ms
      // tsEOR = grpecs->getTimeEnd();                   // ms

      // duration of TF in bcs
      nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;

      // IA - try without "-1" shift
      orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
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

      TFid = 0;
      flagWaitForNewTF = false;
    }

    // bc in Time Frame:
    int64_t bcInTF = (bc.globalBC() - bcSOR) % nBCsPerTF;
    int orbitInTF = bcInTF / nBCsPerOrbit; //(bc.globalBC() - bcSOR) / nBCsPerOrbit; // - orbitSOR;
    // int TFid = (bc.globalBC() - bcSOR) / nBCsPerTF;// ( nOrbitsPerTF * nBCsPerTF);
    if (flagWaitForNewTF && prevBcInTF >= 0 && prevBcInTF > bcInTF) // orbitInTF == 0 && !firstTF )
    {
      TFid++;
      flagWaitForNewTF = false;
    }

    // tick flag if we are about to have a new TF
    if (orbitInTF > 0.5 * nOrbitsPerTF)
      flagWaitForNewTF = true;

    // if (TFid < 10 && (bcInTF < 50 || bcInTF > 3564 * 31 - 200))
    // cout << "QA: bcInTF = " << bcInTF << ", orbitSOR = " << orbitSOR << ", orbitInTF = " << orbitInTF << ", TFid = " << TFid << endl;

    if (flagIncludeQAHistBcInTF) {
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz"), bcInTF, counterPVcontributorsAfterTPCcuts);
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_SEL8_AND_Vz_ReallyAllContrib"), bcInTF, collision.numContrib());
      histosEventBcInTF.fill(HIST("hIsTriggerTVX_vs_bcInTF_BEFORE_SEL8_AND_Vz"), bcInTF, collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ? 1 : 0);
    }

    int64_t foundBcInTF = 0;
    if (collision.has_foundBC()) {
      auto bcFound = collision.foundBC_as<BCsRun3>();
      foundBcInTF = (bcFound.globalBC() - bcSOR) % nBCsPerTF;
    }

    // ### special sub-loop over tracks for numerator of tracks / T0ampl ratios
    if (collision.has_foundBC() && flagIncludeQATracksInFoundBC) {
      int nAllTracks = 0;
      int nTracksPV = 0;

      int nITStracks = 0;
      int nTPCtracks = 0;
      int nITSTPCtracks = 0;
      int nTPCtracks_HandmadeITSROFcut = 0;
      int nTOFtracks = 0;
      int nTOFtracks_HandmadeITSROFcut = 0;
      int nTRDtracks = 0;

      int nTRDTOFtracks_HandmadeITSROFcut = 0;

      int nTRDTOFtracks = 0;
      int nITSTPCTRDtracks = 0;
      int nITSTPCTOFtracks = 0;
      int nITSTPCTRDTOFtracks = 0;

      int nGlobalTracks = 0;

      double timeFromTOFtracks = 0;
      double timeFromTRDtracks = 0;

      for (auto& track : tracks) {

        nAllTracks++;
        if (!track.isPVContributor()) {
          continue;
        }
        nTracksPV++;

        if (track.isGlobalTrack())
          nGlobalTracks++;

        nITStracks += track.hasITS() && !track.hasTPC();
        nTPCtracks += track.hasTPC();
        nITSTPCtracks += track.hasITS() && track.hasTPC();
        nTOFtracks += track.hasTOF();
        nTRDtracks += track.hasTRD() && !track.hasTOF();

        nTRDTOFtracks += track.hasTRD() && track.hasTOF();
        nITSTPCTRDtracks += track.hasITS() && track.hasTPC() && track.hasTRD();
        nITSTPCTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF();
        nITSTPCTRDTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF() && track.hasTRD();

        float eta = track.eta();
        float pt = track.pt();

        if (track.hasITS() && track.hasTPC()) {
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks"), globalFoundBC, track.itsNCls());

          if (fabs(eta) < 0.8 && fabs(collision.posZ()) < 10) {
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCtracks_vs_pt"), globalFoundBC, pt);
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks_vs_pt"), globalFoundBC, pt, track.itsNCls());
          }
        }

        if (track.isGlobalTrack())
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSlayers_for_GlobalTracks"), globalFoundBC, track.itsNCls());

        // ROF study: eta and pt bins
        if (flagIncludeQAROFinPtEtaBins && track.itsNCls() >= 2 && track.tpcNClsFound() > 80) {

          bool isPion = (TMath::Abs(track.tpcNSigmaPi()) < 3 && TMath::Abs(track.tpcNSigmaKa()) > 3 && TMath::Abs(track.tpcNSigmaPr()) > 3) ? true : false;
          bool isProton = (TMath::Abs(track.tpcNSigmaPr()) < 3 && TMath::Abs(track.tpcNSigmaKa()) > 3 && TMath::Abs(track.tpcNSigmaPi()) > 3) ? true : false;
          if (fabs(eta) < 0.2) {
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02"), globalFoundBC, pt);
            if (isPion)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02"), globalFoundBC, pt);
            if (isProton)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02"), globalFoundBC, pt);

            // for this eta: study ITS cluster pattern
            for (int layer = 0; layer < 7; layer++) {
              if (track.itsClusterMap() & (1 << layer)) {
                if (pt > 0.2 && pt < 0.5) {
                  histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pt02_05"), globalFoundBC, layer);
                  if (isPion)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pi_pt02_05"), globalFoundBC, layer);
                  if (isProton)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_p_pt02_05"), globalFoundBC, layer);
                } else if (pt > 0.5 && pt < 1.0) {
                  histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pt05_10"), globalFoundBC, layer);
                  if (isPion)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pi_pt05_10"), globalFoundBC, layer);
                  if (isProton)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_p_pt05_10"), globalFoundBC, layer);
                } else if (pt > 1.0 && pt < 5.0) {
                  histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pt10_50"), globalFoundBC, layer);
                  if (isPion)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_pi_pt10_50"), globalFoundBC, layer);
                  if (isProton)
                    histosEventTracksFoundBC.fill(HIST("hITSlayerCounts_vs_BC_p_pt10_50"), globalFoundBC, layer);
                }
              }
            }
          } else if (fabs(eta) < 0.4) {
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta02_04"), globalFoundBC, pt);
            if (isPion)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta02_04"), globalFoundBC, pt);
            if (isProton)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta02_04"), globalFoundBC, pt);
          } else if (fabs(eta) < 0.6) {
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta04_06"), globalFoundBC, pt);
            if (isPion)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta04_06"), globalFoundBC, pt);
            if (isProton)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta04_06"), globalFoundBC, pt);
          } else if (fabs(eta) < 0.8) {
            histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITStracks_vs_BC_in_ptBins_eta06_08"), globalFoundBC, pt);
            if (isPion)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_pi_nITStracks_vs_BC_in_ptBins_eta06_08"), globalFoundBC, pt);
            if (isProton)
              histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_p_nITStracks_vs_BC_in_ptBins_eta06_08"), globalFoundBC, pt);
          }
        }

        // if (flagBCisNotInHandmadeBoundariesITSROF)
        if (flagFoundBCisNotInHandmadeBoundariesITSROF) {
          nTPCtracks_HandmadeITSROFcut += track.hasTPC();
          nTOFtracks_HandmadeITSROFcut += track.hasTOF();
          nTRDTOFtracks_HandmadeITSROFcut += track.hasTRD() && track.hasTOF();
        }

        // calculate average time using TOF and TRD tracks
        if (track.hasTOF()) {
          timeFromTOFtracks += track.trackTime();
        } else if (track.hasTRD()) {
          timeFromTRDtracks += track.trackTime();
        }
      }

      if (TFid < 2 && collBC < 200) {
        if (nTOFtracks > 0)
          LOGF(info, "QA: average timeFromTOFtracks = %.f", timeFromTOFtracks / nTOFtracks);
        if (nTRDtracks > 0)
          LOGF(info, "QA: average timeFromTRDtracks = %.f", timeFromTRDtracks / nTRDtracks);
      }

      // all collisions
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nAllTracks"), globalFoundBC, nAllTracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTracksPV"), globalFoundBC, nTracksPV);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nGlobalTracks"), globalFoundBC, nGlobalTracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nITStracks"), globalFoundBC, nITStracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTPCtracks"), globalFoundBC, nTPCtracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTOFtracks"), globalFoundBC, nTOFtracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTRDtracks"), globalFoundBC, nTRDtracks);

      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTPCtracks_HandmadeITSROFcut"), globalFoundBC, nTPCtracks_HandmadeITSROFcut);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTOFtracks_HandmadeITSROFcut"), globalFoundBC, nTOFtracks_HandmadeITSROFcut);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTRDTOFtracks_HandmadeITSROFcut"), globalFoundBC, nTRDTOFtracks_HandmadeITSROFcut);

      histosEventTracksFoundBC.fill(HIST("hFoundBC_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
      histosEventTracksFoundBC.fill(HIST("hFoundBC_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);

      // TVX collisions
      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nAllTracks"), globalFoundBC, nAllTracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTracksPV"), globalFoundBC, nTracksPV);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nGlobalTracks"), globalFoundBC, nGlobalTracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITStracks"), globalFoundBC, nITStracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTPCtracks"), globalFoundBC, nTPCtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCtracks"), globalFoundBC, nITSTPCtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTOFtracks"), globalFoundBC, nTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTRDtracks"), globalFoundBC, nTRDtracks);

        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTPCtracks_HandmadeITSROFcut"), globalFoundBC, nTPCtracks_HandmadeITSROFcut);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTOFtracks_HandmadeITSROFcut"), globalFoundBC, nTOFtracks_HandmadeITSROFcut);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTRDTOFtracks_HandmadeITSROFcut"), globalFoundBC, nTRDTOFtracks_HandmadeITSROFcut);

        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);

        // to study ITS ROF dips vs multiplicity
        if (nTracksPV >= 1 && nTracksPV < 20)
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_1_20"), globalFoundBC, nITSTPCTRDtracks);
        else if (nTracksPV >= 20 && nTracksPV < 100)
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_20_100"), globalFoundBC, nITSTPCTRDtracks);
        else if (nTracksPV >= 100 && nTracksPV < 400)
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_100_400"), globalFoundBC, nITSTPCTRDtracks);
        else if (nTracksPV >= 400 && nTracksPV < 1000)
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_400_1000"), globalFoundBC, nITSTPCTRDtracks);
        else if (nTracksPV >= 1000)
          histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVX_nITSTPCTRDtracks_mult_above_1000"), globalFoundBC, nITSTPCTRDtracks);

        // vs bcInTF:
        histosEventTracksFoundBC.fill(HIST("hBCinTF_kTVX_nTracksPV"), bcInTF, nTracksPV);
        histosEventTracksFoundBC.fill(HIST("hBCinTF_kTVX_nITSTPCTRDtracks"), bcInTF, nITSTPCTRDtracks);
      }

      // TVXinTRD collisions
      if (collision.alias_bit(kTVXinTRD)) {
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nAllTracks"), globalFoundBC, nAllTracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nTracksPV"), globalFoundBC, nTracksPV);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nGlobalTracks"), globalFoundBC, nGlobalTracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nITStracks"), globalFoundBC, nITStracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nTPCtracks"), globalFoundBC, nTPCtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nTOFtracks"), globalFoundBC, nTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nTRDtracks"), globalFoundBC, nTRDtracks);

        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nTRDTOFtracks"), globalFoundBC, nTRDTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTRDtracks"), globalFoundBC, nITSTPCTRDtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTOFtracks"), globalFoundBC, nITSTPCTOFtracks);
        histosEventTracksFoundBC.fill(HIST("hFoundBC_kTVXinTRD_nITSTPCTRDTOFtracks"), globalFoundBC, nITSTPCTRDTOFtracks);
      }
    } // end of if has_foundBC and && flagIncludeQATracksInFoundBC

    // is FT0
    bool isFT0 = false;
    double ft0_posZ = 0;
    double diff_PV_ft0_tracks = 0;
    // if (collision.has_foundFT0()) {
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
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

        histosFT0.fill(HIST("hT0_2D_multPV_vs_vertexDiff"), diff_PV_ft0_tracks, multNTracksPV);

        if (collision.t0ACorrectedValid() && collision.t0CCorrectedValid()) {
          double tT0A = collision.t0ACorrected();
          double tT0C = collision.t0CCorrected();
          histosFT0.fill(HIST("hT0A"), tT0A);
          histosFT0.fill(HIST("hT0C"), tT0C);
          histosFT0.fill(HIST("hT0_sum_AC"), tT0A + tT0C);
          histosFT0.fill(HIST("hT0_diff_AC"), tT0A - tT0C);
          histosFT0.fill(HIST("hT0_diff_vs_sum_AC"), tT0A - tT0C, tT0A + tT0C);

          // uncorrected?
          double tT0A_uncorr = ft0.timeA();
          double tT0C_uncorr = ft0.timeC();
          histosFT0.fill(HIST("hT0timeA_uncorr"), tT0A_uncorr);
          histosFT0.fill(HIST("hT0timeC_uncorr"), tT0C_uncorr);
          histosFT0.fill(HIST("hT0_sum_AC_uncorr"), tT0A_uncorr + tT0C_uncorr);
          histosFT0.fill(HIST("hT0_diff_AC_uncorr"), tT0A_uncorr - tT0C_uncorr);
          histosFT0.fill(HIST("hT0_diff_vs_sum_AC_uncorr"), tT0A_uncorr - tT0C_uncorr, tT0A_uncorr + tT0C_uncorr);
        }
      }
    } // end of if collision with FT0

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
    if (!collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      prevBC = collBC;
      prevBcInTF = bcInTF;
      prevOrbit = orbitInTF;
      prevTF = TFid;
      prevGlobalFoundBC = globalFoundBC;
      prevFoundBcInTF = foundBcInTF;
      // double prev_vZ_keepForBelow = prev_vZ;
      prev_vZ = vZ;
      prev_mult = nTracksAll;
      prev_globTrkContrib = nTracksGlobalPVAccepted;

      return;
    }

    // ##### check how often we analyze collision in the same BC (and also the vZ difference)
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      if (prevBC != 9999) {
        int32_t diff = static_cast<int32_t>(collBC) - static_cast<int32_t>(prevBC);
        histosEvent.fill(HIST("hBC_DIFF_to_previous"), diff);
        if (diff == 0)
          histosEvent.fill(HIST("hBC_DIFF_to_previous_vZvZ_2D"), vZ, prev_vZ);
      }
      if (prevBcInTF >= 0) {
        int32_t diffBcInTF = static_cast<int32_t>(bcInTF) - static_cast<int32_t>(prevBcInTF);
        histosEvent.fill(HIST("hBCinTF_DIFF_to_previous"), diffBcInTF);
      }
      if (prevFoundBcInTF >= 0) {
        int32_t diffGlobalBcInTF = static_cast<int32_t>(foundBcInTF) - static_cast<int32_t>(prevFoundBcInTF);
        histosEvent.fill(HIST("hBCinTF_DIFF_to_previous_FOUND_BC"), diffGlobalBcInTF);
      }

      // global found BC:
      if (prevGlobalFoundBC != 9999) {
        int32_t diff = static_cast<int32_t>(globalFoundBC) - static_cast<int32_t>(prevGlobalFoundBC);
        histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC"), diff);

        if (counterPVcontributorsAfterTPCcuts > 0) {
          histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC"), diff);
          if (isFT0 && cutVzTrackT0diffLower < diff_PV_ft0_tracks && diff_PV_ft0_tracks < cutVzTrackT0diffUpper) // cm
            histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut"), diff);
        }

        if (diff == 0) // QA vZ with prev_vZ correlations
        {
          histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_vZvZ_2D_diff0"), vZ, prev_vZ);
          if (counterPVcontributorsAfterTPCcuts > 0) {
            histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_1goodVertContribTPC_vZvZ_2D_diff0"), vZ, prev_vZ);
            if (isFT0 && cutVzTrackT0diffLower < diff_PV_ft0_tracks && diff_PV_ft0_tracks < cutVzTrackT0diffUpper) // cm
              histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_diff0"), vZ, prev_vZ);

            histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_mult_mult_2D_diff0"), nTracksAll, prev_mult);
            histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_globTrkContrib_globTrkContrib_2D_diff0"), nTracksGlobalPVAccepted, prev_globTrkContrib);
          }

          fillHistForThisCut("NoTFborder_FoundBCwithDiff0", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
        }
      }
    }

    // ##### take a look into the next collision (March 2024)
    bool flagNotTheSameFoundBC = false;
    bool flagDiffBnBcIs10 = false;
    // bool flagDiffBnBcIs20 = false;
    bool flagDiffBnBcIs40 = false;
    // int64_t nextOrbitInTF = -1;

    if ((collision.index() + 1 < collision.size())) {
      auto nextColl = collision + 1;

      if (TFid < 10 && collision.index() == 0) {
        LOGF(info, "QA: collision.size = %d, colIndexThis = %d / vZ = %.2f, colIndexNext = %d / vZ = %.2f",
             collision.size(), collision.index(), collision.posZ(), nextColl.index(), nextColl.posZ());
      }

      auto bcNext = nextColl.bc_as<BCsRun3>();
      // auto nextCollBC = bcNext.globalBC() % 3564;
      int64_t nextBcInTF = (bcNext.globalBC() - bcSOR) % nBCsPerTF;

      // check if Found BC is NOT the same for prev, current and next events:
      if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
        flagNotTheSameFoundBC = true;

        // check also the distance (in BC) b/n collisions
        if (prevTF == TFid) {
          int diffWrtPrev = bcInTF - prevBcInTF;
          int diffWrtNext = nextBcInTF - bcInTF;

          if (diffWrtPrev >= 10 && diffWrtNext >= 10)
            flagDiffBnBcIs10 = true;
          // if (diffWrtPrev >= 20 && diffWrtNext >= 20)
          // flagDiffBnBcIs20 = true;
          if (diffWrtPrev >= 40 && diffWrtNext >= 40)
            flagDiffBnBcIs40 = true;
        }
      }
    }

    // save some information from the current collision for comparision with the next one
    prevBC = collBC;
    prevBcInTF = bcInTF;
    prevOrbit = orbitInTF;
    prevTF = TFid;
    prevGlobalFoundBC = globalFoundBC;
    prevFoundBcInTF = foundBcInTF;
    double prev_vZ_keepForBelow = prev_vZ;
    prev_vZ = vZ;
    prev_mult = nTracksAll;
    prev_globTrkContrib = nTracksGlobalPVAccepted;

    if (flagIncludeQAHistBcInTF)
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_BEFORE_Vz"), bcInTF, counterPVcontributorsAfterTPCcuts);

    // ##### now the vZ cut
    if (TMath::Abs(vZ) > vtxCut)
      return;

    // fill after-vZ cut histograms
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
      if (isITSonlyVertex) {
        histosFT0.fill(HIST("hVertex_T0_PV_after_ITSonlyVertex"), ft0_posZ, collision.posZ());
      } else {
        histosFT0.fill(HIST("hVertex_T0_PV_after_ITSonlyVertex_ANTI"), ft0_posZ, collision.posZ());
      }
    }

    histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 2, tracks.size());            // ALL tracks
    histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 2, nTracksGlobalAccepted); // global tracks
    histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 2, 1);                  // collisions counter

    histosEvent.fill(HIST("hNumContribAfterTPCcuts"), counterPVcontributorsAfterTPCcuts);
    histosEvent.fill(HIST("hNumContribITS7hits"), counterPVcontributorsITS7hits);

    if (flagIncludeQAHistBcInTF) {
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF"), bcInTF, counterPVcontributorsAfterTPCcuts);
      histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_ReallyAllContrib"), bcInTF, collision.numContrib());
      histosEventBcInTF.fill(HIST("hIsTriggerTVX_vs_bcInTF"), bcInTF, collision.selection_bit(o2::aod::evsel::kIsTriggerTVX) ? 1 : 0);

      histosEventBcInTF.fill(HIST("hGlobalTracks_vs_bcInTF"), bcInTF, nTracksGlobalAccepted);
    }
    histosEvent.fill(HIST("hOrbitStartFromCollIndexZeroAft"), static_cast<int32_t>(orbit) - static_cast<int32_t>(orbitAtCollIndexZero));
    histosEvent.fill(HIST("h2D_Orbit_vs_CollIndex_Aft"), collision.index(), static_cast<int32_t>(orbit) - static_cast<int32_t>(orbitAtCollIndexZero));

    histosEvent.fill(HIST("hMF"), magneticField);
    int MFsign = magneticField > 0 ? +1 : -1;

    histosEvent.fill(HIST("vtxChi2Aft"), collision.chi2());
    histosEvent.fill(HIST("vtxCutsAft"), vZ);
    fillHistForThisCut("Aft", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
    histosEvent.fill(HIST("hNtrackshGlobalAft"), nTracksGlobalAccepted);

    // cut TF borders
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      fillHistForThisCut("TFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

      histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 3, tracks.size());            // ALL tracks with TF cut
      histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 3, nTracksGlobalAccepted); // global tracks
      histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 3, 1);                  // collisions counter

      if (flagIncludeQAHistBcInTF) {
        histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_ReallyAllContrib_AfterTimeFrameCut"), bcInTF, collision.numContrib());
        histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_AfterTimeFrameCut"), bcInTF, counterPVcontributorsAfterTPCcuts);
      }

      histosEvent.fill(HIST("hNtrackshGlobalAft_AfterTimeFrameCut"), nTracksGlobalAccepted);

      if (isFT0) {
        histosFT0.fill(HIST("hVertex_T0_PV_after_TFcut"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_after_TFcut"), diff_PV_ft0_tracks);
      }

      if (flagNotTheSameFoundBC) {
        histosEvent.fill(HIST("hBC_DIFF_to_previous_FOUND_BC_1goodCont_FTvertexCut_vZvZ_2D_NoBCdiff0"), vZ, prev_vZ_keepForBelow);
        fillHistForThisCut("NoTFborder_NoFoundBCwithDiff0", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

        if (isFT0)
          histosFT0.fill(HIST("hVertex_T0_PV_after_NoBCdiff0"), ft0_posZ, collision.posZ());
      } else {
        if (isFT0)
          histosFT0.fill(HIST("hVertex_T0_PV_after_NoBCdiff0_ANTI"), ft0_posZ, collision.posZ());
      }

      if (flagDiffBnBcIs10)
        fillHistForThisCut("NoTFborder_DiffBnBcAtLeast10", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      // if (flagDiffBnBcIs20)
      // fillHistForThisCut("NoTFborder_DiffBnBcAtLeast20", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      if (flagDiffBnBcIs40)
        fillHistForThisCut("NoTFborder_DiffBnBcAtLeast40", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
    }

    histosEvent.fill(HIST("hCollIndexAft"), collision.index());
    histosEvent.fill(HIST("hCollTimeAft"), collision.collisionTime());
    histosEvent.fill(HIST("hCollTimeResAft"), collision.collisionTimeRes());
    histosEvent.fill(HIST("hNumContribAft"), collision.numContrib());
    histosEvent.fill(HIST("hNtracksAft"), tracks.size());
    // histosEvent.fill(HIST("h2D_numContrib_vs_collIndex"), collision.index(), collision.numContrib());
    histosEvent.fill(HIST("hBC_Aft"), collBC);
    histosEvent.fill(HIST("hBCFound_Aft"), globalFoundBC);
    histosEvent.fill(HIST("h2D_numContrib_vs_BC"), collBC, collision.numContrib());

    int64_t diffFoundBC_vs_BC = static_cast<int64_t>(globalFoundBC) - static_cast<int64_t>(collBC);
    histosEvent.fill(HIST("h2D_diffFoundBC_vs_BC"), collBC, diffFoundBC_vs_BC);

    if (collision.has_foundBC())
      histosEvent.fill(HIST("h2D_diffFoundBC_vs_BC_inTF"), collBC, static_cast<int64_t>(foundBcInTF) - static_cast<int64_t>(bcInTF));

    // with FT0 conditions
    if (isFT0) {
      histosEvent.fill(HIST("hBC_vertNcontr3_with_FT0"), collBC);
      if (cutVzTrackT0diffLower < diff_PV_ft0_tracks && diff_PV_ft0_tracks < cutVzTrackT0diffUpper) //  < 1)
        histosEvent.fill(HIST("hBC_vertNcontr3_with_FT0_diffPV_1cm"), collBC);
    }

    // nTracks vs BC
    if (flagIncludeQAHistBcVsNtracks) {
      histosEvent.fill(HIST("h2D_nTracksBeforeCuts_vs_BC"), collBC, nTracksAll);
      histosEvent.fill(HIST("h2D_nTracksAfterEtaTPCcuts_vs_BC"), collBC, nTracksAfterEtaTPCCuts);
      histosEvent.fill(HIST("h2D_nTracksITSonly_vs_BC"), collBC, nTracksAfterEtaCutsITSonly);
      histosEvent.fill(HIST("h2D_nTracksWithITS_vs_BC"), collBC, nTracksWithITS);
      histosEvent.fill(HIST("h2D_nTracksWithITS7hits_vs_BC"), collBC, nTracksWithITS7hits);
      histosEvent.fill(HIST("h2D_nTracksWithITSandTPC_vs_BC"), collBC, nTracksWithITSandTPC);
      histosEvent.fill(HIST("h2D_nTracksWithTRD_vs_BC"), collBC, nTracksWithTRD);
      histosEvent.fill(HIST("h2D_nTracksWithTOF_vs_BC"), collBC, nTracksWithTOF);
      histosEvent.fill(HIST("h2D_nTracksWithTRDorTOF_vs_BC"), collBC, nTracksWithTRDorTOF);
      histosEvent.fill(HIST("h2D_nTracksGlobal_vs_BC"), collBC, nTracksGlobalAccepted);

      histosEvent.fill(HIST("h2D_nTracksGlobalWithITS7hits_vs_BC"), collBC, nTracksGlobalWithITS7hits);
      histosEvent.fill(HIST("h2D_nTracksGlobalWithTRDorTOF_vs_BC"), collBC, nTracksGlobalWithTRDorTOF);
    }

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
    if (flagIncludeQAHistBcVsNtracks && (counterPVcontributorsNoTOFandTRD == counterPVcontributorsAfterTPCcuts) && counterPVcontributorsNoTOFandTRD > 0) {
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

    // test RO frame effect
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      if (flagIncludeQAHistBcInTF)
        histosEventBcInTF.fill(HIST("hNumContrib_vs_bcInTF_After_ITS_ROF_cut"), bcInTF, counterPVcontributorsAfterTPCcuts);

      histosEventCounters.fill(HIST("hNtracksAll_vs_variousCuts"), 4, tracks.size());            // ALL tracks outside ITS ROF cut
      histosEventCounters.fill(HIST("hNtracksGlobal_vs_variousCuts"), 4, nTracksGlobalAccepted); // global tracks
      histosEventCounters.fill(HIST("hNtotalCollisions_vs_variousCuts"), 4, 1);                  // collisions counter

      histosEvent.fill(HIST("hCollIndex_ITS_ROF_cut"), collision.index());
      histosEvent.fill(HIST("hCollTime_ITS_ROF_cut"), collision.collisionTime());
      histosEvent.fill(HIST("hCollTimeRes_ITS_ROF_cut"), collision.collisionTimeRes());
      // histosEvent.fill(HIST("h2D_numContrib_vs_collIndex_ITS_ROF_cut"), collision.index(), collision.numContrib());
      histosEvent.fill(HIST("hBC_ITS_ROF_cut"), collBC);

      fillHistForThisCut("ITSROFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

      if (isFT0) {
        histosFT0.fill(HIST("hVertex_T0_PV_after_ITSROFcut"), ft0_posZ, collision.posZ());
        histosFT0.fill(HIST("hT0vertexDiff_after_ITSROFcut"), diff_PV_ft0_tracks);
      }

      // in addition: TF border cuts:
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        if (isFT0) {
          histosFT0.fill(HIST("hVertex_T0_PV_after_ITSROF_and_TFcut"), ft0_posZ, collision.posZ());
          histosFT0.fill(HIST("hT0vertexDiff_after_ITSROF_and_TFcut"), diff_PV_ft0_tracks);
        }
      }
    } else {
      // look into what is within the ITSROF cut borders
      fillHistForThisCut("antiITSROFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
    }

    // global PV contributors with 7 ITS hits
    if (nTracksGlobalPVwithITS7hits >= 1) // 2)
      fillHistForThisCut("1globalPVcontrib_ITS7hits", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    // global PV contributors with TRD or TOF signal
    if (nTracksGlobalPVwithTRDorTOF >= 1) // 2)
      fillHistForThisCut("1globalPVcontrib_TRDorTOF", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    // diffFoundBC_vs_BC cut
    if (diffFoundBC_vs_BC == 0)
      fillHistForThisCut("diffFoundBC_vs_BC_0", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    // cut on diff b/n vertex from FT0 and track-based
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      fillHistForThisCut("PV_FT0_diff_cut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

      // ALL OTHER CUTS:
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && flagNotTheSameFoundBC && !isITSonlyVertex) //&& nTracksGlobalPVAccepted >= 1)
      {
        fillHistForThisCut("ALL_CUTS", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

        // cut "strange" events below the main 2D diagonal trend
        if (multNTracksPV > funcCutEventsByMultPVvsV0A->Eval(multV0A))
          fillHistForThisCut("ALL_CUTS_CutEventsByMultPVvsV0A", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      }
    }

    // cut "strange" events below the main 2D diagonal trend
    if (multNTracksPV > funcCutEventsByMultPVvsV0A->Eval(multV0A))
      fillHistForThisCut("CutEventsByMultPVvsV0A", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
    else // QA for events we actually cut out
    {
      fillHistForThisCut("antiCutEventsByMultPVvsV0A", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      histosEvent.fill(HIST("hBC_Aft_rejectedByCutOnMultPVvsV0A"), collBC);
    }

    // the cut used for pp 2023 to capture the ITS ROF boundaries:
    if (multNTracksPV < funcCutEventsByMultPVvsT0C->Eval(multT0C)) {
      histosEvent.fill(HIST("hBC_Aft_rejectedByCutOnMultPVvsT0C"), collBC);
      // check if this vertex is ITS-only
      if (isITSonlyVertex)
        histosEvent.fill(HIST("hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV"), collBC);
      else
        histosEvent.fill(HIST("hBC_Aft_rejectedByCutOnMultPVvsT0C_AndIfITSonlyPV_ANTI"), collBC);
    }

    // cut the remaining "strange" events just below the main 2D diagonal trend on multPV vs V0A plot
    bool flagBelowRemainingTrendMultPVvsV0A = false;
    if (multNTracksPV > funcCutEventsByMultPVvsV0A_anotherTrend->Eval(multV0A))
      flagBelowRemainingTrendMultPVvsV0A = true;

    if (isITSonlyVertex)
      fillHistForThisCut("isITSonlyVertex", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
    else
      fillHistForThisCut("antiIsITSonlyVertex", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    if (collision.alias_bit(kTVXinTRD))
      fillHistForThisCut("kTVXinTRD", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

    if (flagBCisNotInHandmadeBoundariesITSROF)
      if (multNTracksPV < funcCutEventsByMultPVvsT0C->Eval(multT0C))
        histosEvent.fill(HIST("hBC_Aft_rejectedByCutOnMultPVvsT0C_after_Handmade_ITSROF_cut"), collBC);

    // QA after hand-made ITSROF cut
    if (flagBCisNotInHandmadeBoundariesITSROF) {
      histosEvent.fill(HIST("hBC_Aft_Handmade_ITSROFcut"), collBC);
      fillHistForThisCut("Handmade_ITSROFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      if (isFT0)
        histosFT0.fill(HIST("hVertex_T0_PV_after_HandmadeITSROFcut"), ft0_posZ, collision.posZ());

      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
        fillHistForThisCut("Handmade_ITSROF_and_TF_cuts", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

      if (isFT0 && (cutVzTrackT0diffLower < diff_PV_ft0_tracks) && (diff_PV_ft0_tracks < cutVzTrackT0diffUpper) && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && flagNotTheSameFoundBC && !isITSonlyVertex) {
        fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
        if (flagBelowRemainingTrendMultPVvsV0A)
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_anotherCutOnMultV0A", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

        if (collision.alias_bit(kTVXinTRD))
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_kTVXinTRD", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

        // study different BC ranges of the orbit:
        if (collBC < vSplitBCpointsOfTheOrbit->at(0))
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_Orbit_1_part", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
        else if (collBC >= vSplitBCpointsOfTheOrbit->at(0) && collBC < vSplitBCpointsOfTheOrbit->at(1))
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_Orbit_2_part", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
        else if (collBC >= vSplitBCpointsOfTheOrbit->at(1) && collBC < vSplitBCpointsOfTheOrbit->at(2))
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_Orbit_3_part", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
        else if (collBC >= vSplitBCpointsOfTheOrbit->at(2))
          fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_Orbit_4_part", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      }
    } else {
      fillHistForThisCut("antiHandmade_ITSROFcut", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
      if (isFT0)
        histosFT0.fill(HIST("hVertex_T0_PV_after_HandmadeITSROFcut_ANTI"), ft0_posZ, collision.posZ());
    }

    // QA after hand-made ITSROF cut ON FOUND BC
    if (flagFoundBCisNotInHandmadeBoundariesITSROF) {
      fillHistForThisCut("Handmade_ITSROFcut_OnFoundBC", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);

      if (isFT0 && (cutVzTrackT0diffLower < diff_PV_ft0_tracks) && (diff_PV_ft0_tracks < cutVzTrackT0diffUpper) && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && flagNotTheSameFoundBC && !isITSonlyVertex) {
        fillHistForThisCut("ALL_CUTS_Handmade_ITSROFcut_OnFoundBC", multNTracksPV, multTrk, nTracksGlobalAccepted, multT0A, multT0C, multV0A, t0cCentr, collBC);
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
    if (flagShowInfo) {
      LOGF(info, "#########");
      LOGF(info, "Collision index : %d", collision.index());
      LOGF(info, "Number of tracks: %d", tracks.size());
      LOGF(info, "Number of v0s   : %d", v0s.size());
      LOGF(info, "Collision: %d [N = %d out of %d], vZ pos = %.1f ",                      //, globalBC=%.1f",
           collision.globalIndex(), tracks.size(), tracks.tableSize(), collision.posZ()); //, collision.globalBC() );
      LOGF(info, "The collision time:   %f", collision.collisionTime());
      LOGF(info, "Tracks for this collision: %d", tracks.size());
    }

    bool flagShowInfoFirstTrack = false;

    // #### loop over tracks
    for (auto const& track : tracks) {
      // eta
      histosTracks.fill(HIST("eta"), track.eta());
      if (fabs(track.eta()) > 0.8)
        continue;
      histosTracks.fill(HIST("etaAfter08cut"), track.eta());

      histosTracks.fill(HIST("hTpcNClsCrossedRows"), track.tpcNClsCrossedRows());
      histosTracks.fill(HIST("hTpcNClsFound"), track.tpcNClsFound());
      histosTracks.fill(HIST("hNumITSclusters"), track.itsNCls());

      if (track.isGlobalTrack()) {
        histosTracks.fill(HIST("hTpcNClsFound_Global_beforeCrossedRowsCut"), track.tpcNClsFound());
        histosTracks.fill(HIST("hTpcNClsCrossedRows_Global_beforeCrossedRowsCut"), track.tpcNClsCrossedRows());
        histosTracks.fill(HIST("hNumITSclusters_Global_beforeCrossedRowsCut"), track.itsNCls());

        histosTracks.fill(HIST("hNumITScls_vs_TPCcrossedRows_Global_beforeCrossedRowsCut"), track.itsNCls(), track.tpcNClsCrossedRows());
        histosTracks.fill(HIST("hDcaXY_Global_beforeCrossedRowsCut"), track.dcaXY());
        histosTracks.fill(HIST("hDcaZ_Global_beforeCrossedRowsCut"), track.dcaZ());
      }

      if (isITSonlyVertex)
        histosTracks.fill(HIST("hNumITSclusters_ITSonlyVert"), track.itsNCls());
      histosTracks.fill(HIST("hDcaXY"), track.dcaXY());
      histosTracks.fill(HIST("hDcaZ"), track.dcaZ());

      // #### detailed pt study:
      if (1) {
        bool noTF = collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder);
        bool noROF = collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder);

        fillPtHistos("All", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);
        if (track.isGlobalTrack())
          fillPtHistos("Global", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);

        if (track.itsNCls() == 7)
          fillPtHistos("ITS7hits", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);
        if (track.itsNCls() == 7 && track.tpcNClsFound() >= 80)
          fillPtHistos("ITS7hits_TPC80cl", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);

        if (track.itsNCls() >= 4)
          fillPtHistos("ITS4567", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);
        if (track.itsNCls() >= 4 && track.tpcNClsFound() >= 80)
          fillPtHistos("ITS4567_TPC80cl", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);

        if (track.itsNCls() >= 5)
          fillPtHistos("ITS567", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);
        if (track.itsNCls() >= 5 && track.tpcNClsFound() >= 80)
          fillPtHistos("ITS567_TPC80cl", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);

        if (track.itsNCls() >= 6)
          fillPtHistos("ITS67", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);
        if (track.itsNCls() >= 6 && track.tpcNClsFound() >= 80)
          fillPtHistos("ITS67_TPC80cl", track.pt(), track.sign(), track.tpcNClsFound(), noTF, noROF);

      } // #### end of detailed pt study:

      histosTracks.fill(HIST("phiHistogram"), track.phi());

      histosTracks.fill(HIST("hTrackLength"), track.length()); // from TrackExtras

      if (track.itsNCls() == 7)
        histosTracks.fill(HIST("hTpcNClsCrossedRowsITS7hits"), track.tpcNClsCrossedRows());

      if (flagShowInfoFirstTrack) {
        LOGF(info, "track.length() = %.2f", track.length());
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

        histosTracks.fill(HIST("etaAftCutsGlobal"), eta);

        histosTracks.fill(HIST("hNumITSclusters_AftCuts_Global"), track.itsNCls());

        histosTracks.fill(HIST("hNumITScls_vs_TPCcls_AftCuts_Global"), track.itsNCls(), track.tpcNClsFound());
        histosTracks.fill(HIST("hNumITScls_vs_TPCcrossedRows_AftCuts_Global"), track.itsNCls(), track.tpcNClsCrossedRows());

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
      } // end of isGlobalTrack

      if (track.itsNCls() == 7) {
        histosTracks.fill(HIST("hTpcNClsCrossedRowsITS7hits_AftCuts"), track.tpcNClsCrossedRows());
        histosTracks.fill(HIST("hTpcNClsFound_ITS7hits_AftCuts"), track.tpcNClsFound());
        histosTracks.fill(HIST("etaITS7hits_AftCuts"), track.eta());

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
      histosTracks.fill(HIST("hDcaZ_AftCuts"), track.dcaZ());

      histosTracks.fill(HIST("etaAftCuts"), track.eta());

      histosTracks.fill(HIST("ptHistogram_AftCuts"), track.pt());
      histosTracks.fill(HIST("phiHistogram_AftCuts"), track.phi());

      histosTracks.fill(HIST("hTrackLength_AftCuts"), track.length()); // from TrackExtras
      // histosTracks.fill(HIST("hChi2perDOF_AftCuts"), track.chi2perNDF());

    } // end of track loop

    int occupancy = collision.trackOccupancyInTimeRange();
    histosEvent.fill(HIST("hOccupancyInTimeRange"), occupancy);
    histosEvent.fill(HIST("hOccupancyInTimeRangeIsAvailable"), occupancy);
    histosEvent.fill(HIST("hCentralityVsOccupancy"), occupancy, t0cCentr);

    // n events in occupancy classes
    if (occupancy >= 0 && t0cCentr < 95 && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      histosK0S.fill(HIST("nEventsVsCentrAllOccupancies"), t0cCentr);

      if (occupancy < 250)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_0_250"), t0cCentr);
      else if (occupancy >= 250 && occupancy < 500)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_250_500"), t0cCentr);
      else if (occupancy >= 500 && occupancy < 750)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_500_750"), t0cCentr);
      else if (occupancy >= 750 && occupancy < 1500)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_750_1500"), t0cCentr);
      else if (occupancy >= 1500 && occupancy < 3000)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_1500_3000"), t0cCentr);
      else if (occupancy >= 3000 && occupancy < 4500)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_3000_4500"), t0cCentr);
      else if (occupancy >= 4500 && occupancy < 6000)
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_4500_6000"), t0cCentr);
      else
        histosK0S.fill(HIST("nEventsVsCentrAfterOccupCut_above_6000"), t0cCentr);
    }

    // ##### v0 analysis
    if (flagIncludeQAHistK0S && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
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

        if (occupancy >= 0 && t0cCentr < 95 && collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder) && posDaughterTrack.tpcNClsFound() > 80 && negDaughterTrack.tpcNClsFound() > 80 && posDaughterTrack.itsNCls() >= 4 && negDaughterTrack.itsNCls() >= 4) {
          histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPt_AllOccupancies"), v0.pt(), v0.mK0Short(), t0cCentr);

          if (occupancy < 250)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_0_250"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 250 && occupancy < 500)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_250_500"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 500 && occupancy < 750)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_500_750"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 750 && occupancy < 1500)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_750_1500"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 1500 && occupancy < 3000)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_1500_3000"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 3000 && occupancy < 4500)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_3000_4500"), v0.pt(), v0.mK0Short(), t0cCentr);
          else if (occupancy >= 4500 && occupancy < 6000)
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_4500_6000"), v0.pt(), v0.mK0Short(), t0cCentr);
          else
            histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtAfterOccupCut_above_6000"), v0.pt(), v0.mK0Short(), t0cCentr);
        }

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

        if (etaDauPos > 0.1 && etaDauPos < 0.8 && etaDauNeg > 0.1 && etaDauNeg < 0.8)
          histosK0S.fill(HIST("hMassK0ShortAfterSelectionVsPtEta01_08"), v0.pt(), v0.mK0Short());
        if (etaDauPos > -0.8 && etaDauPos < -0.1 && etaDauNeg > -0.8 && etaDauNeg < -0.1)
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
    } // end of if (flagIncludeQAHistK0S)
  } // end of processRobustFluctuationObservables()

  // shortcut function to fill 2D histograms
  void fillHistForThisCut(string cutName, int multNTracksPV, int multTrk, int nTracksGlobalAccepted, double multT0A, double multT0C, double multV0A, double /*t0cCentr*/, int bc)
  {
    // registry.get<TH1>(HIST("eta"))->Fill(track.eta());
    // arrPointers[histId][cutId]->Fill(xval, yval, weight);
    if (!mMyCuts.count(cutName)) {
      LOGF(fatal, "AHTUNG! no key %s in mMyCuts map!", cutName);
      return;
    }
    // --> key exists
    int cutId = mMyCuts[cutName];
    // cout << "cutName = " << cutName << ", cutId = " << cutId << endl;
    // registry.add<thetype>(TString::Format("%s/%s", thedirectory, thename).Data(), thetitle, thekind, thebinning)

    FILL_QA_HIST_1D(cutId, "nGlobalAcceptedTracks", nTracksGlobalAccepted);
    FILL_QA_HIST_1D(cutId, "bcWithCuts", bc);

    FILL_QA_HIST_2D(cutId, "multAllTr_vs_multT0C", multT0C, multTrk);
    FILL_QA_HIST_2D(cutId, "multAllTr_vs_multV0A", multV0A, multTrk);
    FILL_QA_HIST_2D(cutId, "multAllTr_vs_multT0A", multT0A, multTrk);
    FILL_QA_HIST_2D(cutId, "multAllTr_vs_multTrkPV", multNTracksPV, multTrk);
    FILL_QA_HIST_2D(cutId, "multGlobalTr_vs_multT0C", multT0C, nTracksGlobalAccepted);
    FILL_QA_HIST_2D(cutId, "multGlobalTr_vs_multV0A", multV0A, nTracksGlobalAccepted);
    FILL_QA_HIST_2D(cutId, "multGlobalTr_vs_multT0A", multT0A, nTracksGlobalAccepted);
    FILL_QA_HIST_2D(cutId, "multGlobalTr_vs_multTrkPV", multNTracksPV, nTracksGlobalAccepted);
    FILL_QA_HIST_2D(cutId, "multTrkPV_vs_multT0C", multT0C, multNTracksPV);
    FILL_QA_HIST_2D(cutId, "multTrkPV_vs_multV0A", multV0A, multNTracksPV);
    FILL_QA_HIST_2D(cutId, "multTrkPV_vs_multT0A", multT0A, multNTracksPV);
    FILL_QA_HIST_2D(cutId, "multT0C_vs_multT0A", multT0A, multT0C);
    FILL_QA_HIST_2D(cutId, "multV0A_vs_multT0A", multT0A, multV0A);
    FILL_QA_HIST_2D(cutId, "multV0A_vs_multT0C", multT0C, multV0A);
    // if (flagPbPb) {
    //   FILL_QA_HIST_2D(cutId, "multAllTr_vs_Cent", t0cCentr, multTrk);
    //   FILL_QA_HIST_2D(cutId, "multGlobalTr_vs_Cent", t0cCentr, nTracksGlobalAccepted);
    //   FILL_QA_HIST_2D(cutId, "multTrkPV_vs_Cent", t0cCentr, multNTracksPV);
    //   FILL_QA_HIST_2D(cutId, "multV0A_vs_Cent", t0cCentr, multV0A);
    //   FILL_QA_HIST_2D(cutId, "multT0C_vs_Cent", t0cCentr, multT0C);
    // }
  }

  void fillPtHistos(string strTrackType, float pt, int charge, float w, bool noTF, bool noROF)
  {
    fillPtHistosThisCut(strTrackType, "allBC", pt, charge, w);

    // TF
    if (noTF)
      fillPtHistosThisCut(strTrackType, "noTFborder", pt, charge, w);
    else
      fillPtHistosThisCut(strTrackType, "TFborder", pt, charge, w);

    // ROF
    if (noROF)
      fillPtHistosThisCut(strTrackType, "noROFborder", pt, charge, w);
    else
      fillPtHistosThisCut(strTrackType, "ROFborder", pt, charge, w);

    // TF, ROF
    if (noTF && noROF)
      fillPtHistosThisCut(strTrackType, "noTFandROFborder", pt, charge, w);
  }

  void fillPtHistosThisCut(string strTrackType, string strEvSelType, float pt, int charge, float w)
  {
    string strPre = strTrackType + "/" + strEvSelType;

    FILL_PT_HIST_1D(strPre + "/pt", pt, 1);
    FILL_PT_HIST_1D(strPre + "/pt_w_nTPCcls", pt, w);

    if (charge > 0) {
      FILL_PT_HIST_1D(strPre + "/pos_pt", pt, 1);
      FILL_PT_HIST_1D(strPre + "/pos_pt_w_nTPCcls", pt, w);
    } else {
      FILL_PT_HIST_1D(strPre + "/neg_pt", pt, 1);
      FILL_PT_HIST_1D(strPre + "/neg_pt_w_nTPCcls", pt, w);
    }
  }

  PROCESS_SWITCH(RobustFluctuationObservables, processRobustFluctuationObservables, "Process RobustFluctuationObservables", true);
};

// #### Supplementary task to iterate over bunch crossings

struct RobustFluctuationBunchCrossingQA {
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  void init(InitContext&)
  {
    AxisSpec axisBC{3601, -0.5, 3600.5, "bc"};
    histos.add("hBC_allBC", "hBC_allBC", kTH1D, {axisBC});
    histos.add("hBC_hasFT0", "hBC_hasFT0", kTH1D, {axisBC});
    histos.add("hBC_kTVX", "hBC_kTVX", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD", "hBC_kTVXinTRD", kTH1D, {axisBC});

    histos.add("hBC_FT0_ampl", "hBC_FT0_ampl", kTH1D, {axisBC});
    histos.add("hBC_FT0_ampl_kTVX", "hBC_FT0_ampl_kTVX", kTH1D, {axisBC});
    histos.add("hBC_FT0_ampl_kTVXinTRD", "hBC_FT0_ampl_kTVXinTRD", kTH1D, {axisBC});

    // ##### tracks for All collisions
    histos.add("hBC_nAllTracks", "hBC_nAllTracks", kTH1D, {axisBC});
    histos.add("hBC_nTracksPV", "hBC_nTracksPV", kTH1D, {axisBC});
    // histos.add("hBC_nGlobalTracks", "hBC_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hBC_nITStracks", "hBC_nITStracks", kTH1D, {axisBC});
    histos.add("hBC_nTPCtracks", "hBC_nTPCtracks", kTH1D, {axisBC});
    histos.add("hBC_nTOFtracks", "hBC_nTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_nTRDtracks", "hBC_nTRDtracks", kTH1D, {axisBC});

    histos.add("hBC_nTRDTOFtracks", "hBC_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_nITSTPCTRDtracks", "hBC_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hBC_nITSTPCTOFtracks", "hBC_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_nITSTPCTRDTOFtracks", "hBC_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

    // ##### tracks for TVX collisions
    histos.add("hBC_kTVX_nAllTracks", "hBC_kTVX_nAllTracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nTracksPV", "hBC_kTVX_nTracksPV", kTH1D, {axisBC});
    // histos.add("hBC_kTVX_nGlobalTracks", "hBC_kTVX_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nITStracks", "hBC_kTVX_nITStracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nTPCtracks", "hBC_kTVX_nTPCtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nTOFtracks", "hBC_kTVX_nTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nTRDtracks", "hBC_kTVX_nTRDtracks", kTH1D, {axisBC});

    histos.add("hBC_kTVX_nTRDTOFtracks", "hBC_kTVX_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nITSTPCTRDtracks", "hBC_kTVX_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nITSTPCTOFtracks", "hBC_kTVX_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVX_nITSTPCTRDTOFtracks", "hBC_kTVX_nITSTPCTRDTOFtracks", kTH1D, {axisBC});

    // ##### tracks for TVXinTRD collisions
    histos.add("hBC_kTVXinTRD_nAllTracks", "hBC_kTVXinTRD_nAllTracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nTracksPV", "hBC_kTVXinTRD_nTracksPV", kTH1D, {axisBC});
    // histos.add("hBC_kTVXinTRD_nGlobalTracks", "hBC_kTVXinTRD_nGlobalTracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nITStracks", "hBC_kTVXinTRD_nITStracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nTPCtracks", "hBC_kTVXinTRD_nTPCtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nTOFtracks", "hBC_kTVXinTRD_nTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nTRDtracks", "hBC_kTVXinTRD_nTRDtracks", kTH1D, {axisBC});

    histos.add("hBC_kTVXinTRD_nTRDTOFtracks", "hBC_kTVXinTRD_nTRDTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nITSTPCTRDtracks", "hBC_kTVXinTRD_nITSTPCTRDtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nITSTPCTOFtracks", "hBC_kTVXinTRD_nITSTPCTOFtracks", kTH1D, {axisBC});
    histos.add("hBC_kTVXinTRD_nITSTPCTRDTOFtracks", "hBC_kTVXinTRD_nITSTPCTRDTOFtracks", kTH1D, {axisBC});
  }

  using BCsWithBcSelsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
  using Colls = soa::Join<aod::Collisions, aod::EvSels>;

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;
  void processRobustFluctuationBunchCrossingQA(
    // aod::Collisions const& cols
    Colls const& cols,
    FullTracksIU const& tracks,
    BCsWithBcSelsRun3 const& bcs,
    aod::FT0s const&)
  {
    if (bcs.size() == 0)
      return;

    for (auto& bc : bcs) {
      auto collBC = bc.globalBC() % 3564;

      histos.fill(HIST("hBC_allBC"), collBC);

      // float multFV0A = 0.f;
      float multFT0A = 0.f;
      float multFT0C = 0.f;
      // float multFDDA = 0.f;
      // float multFDDC = 0.f;

      if (bc.has_ft0()) {
        histos.fill(HIST("hBC_hasFT0"), collBC);

        auto ft0 = bc.ft0();
        for (auto amplitude : ft0.amplitudeA()) {
          multFT0A += amplitude;
        }
        for (auto amplitude : ft0.amplitudeC()) {
          multFT0C += amplitude;
        }

        histos.fill(HIST("hBC_FT0_ampl"), collBC, multFT0A + multFT0C);
      }

      // check TVX trigger
      if (bc.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        histos.fill(HIST("hBC_kTVX"), collBC);
        if (bc.has_ft0())
          histos.fill(HIST("hBC_FT0_ampl_kTVX"), collBC, multFT0A + multFT0C);
      }

      // check kTVXinTRD
      if (bc.alias_bit(kTVXinTRD)) {
        histos.fill(HIST("hBC_kTVXinTRD"), collBC);
        if (bc.has_ft0())
          histos.fill(HIST("hBC_FT0_ampl_kTVXinTRD"), collBC, multFT0A + multFT0C);
      }

      // if (bc.has_fdd()) { // to be used for event selection in the future?
      // }
    }

    // collisions & tracks loop
    for (auto& col : cols) {
      auto bc = col.bc_as<BCsWithBcSelsRun3>();
      auto collBC = bc.globalBC() % 3564;

      // count tracks of different types
      int nAllTracks = 0;
      int nTracksPV = 0;
      // int nGlobalTracks = 0;

      int nITStracks = 0;
      int nTPCtracks = 0;
      int nTOFtracks = 0;
      int nTRDtracks = 0;

      int nTRDTOFtracks = 0;
      int nITSTPCTRDtracks = 0;
      int nITSTPCTOFtracks = 0;
      int nITSTPCTRDTOFtracks = 0;

      auto tracksGrouped = tracks.sliceBy(perCollision, col.globalIndex());
      for (auto& track : tracksGrouped) {
        nAllTracks++;
        if (!track.isPVContributor()) {
          continue;
        }
        nTracksPV++;

        // if (track.isGlobalTrack())
        //   nGlobalTracks++;

        nITStracks += track.hasITS() && !track.hasTPC();
        nTPCtracks += track.hasTPC();
        nTOFtracks += track.hasTOF();
        nTRDtracks += track.hasTRD() && !track.hasTOF();

        nTRDTOFtracks += track.hasTRD() && track.hasTOF();
        nITSTPCTRDtracks += track.hasITS() && track.hasTPC() && track.hasTRD();
        nITSTPCTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF();
        nITSTPCTRDTOFtracks += track.hasITS() && track.hasTPC() && track.hasTOF() && track.hasTRD();
      }

      // all collisions
      histos.fill(HIST("hBC_nAllTracks"), collBC, nAllTracks);
      histos.fill(HIST("hBC_nTracksPV"), collBC, nTracksPV);
      // histos.fill(HIST("hBC_nGlobalTracks"), collBC, nGlobalTracks);
      histos.fill(HIST("hBC_nITStracks"), collBC, nITStracks);
      histos.fill(HIST("hBC_nTPCtracks"), collBC, nTPCtracks);
      histos.fill(HIST("hBC_nTOFtracks"), collBC, nTOFtracks);
      histos.fill(HIST("hBC_nTRDtracks"), collBC, nTRDtracks);

      histos.fill(HIST("hBC_nTRDTOFtracks"), collBC, nTRDTOFtracks);
      histos.fill(HIST("hBC_nITSTPCTRDtracks"), collBC, nITSTPCTRDtracks);
      histos.fill(HIST("hBC_nITSTPCTOFtracks"), collBC, nITSTPCTOFtracks);
      histos.fill(HIST("hBC_nITSTPCTRDTOFtracks"), collBC, nITSTPCTRDTOFtracks);

      // TVX collisions
      if (col.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        histos.fill(HIST("hBC_kTVX_nAllTracks"), collBC, nAllTracks);
        histos.fill(HIST("hBC_kTVX_nTracksPV"), collBC, nTracksPV);
        // histos.fill(HIST("hBC_kTVX_nGlobalTracks"), collBC, nGlobalTracks);
        histos.fill(HIST("hBC_kTVX_nITStracks"), collBC, nITStracks);
        histos.fill(HIST("hBC_kTVX_nTPCtracks"), collBC, nTPCtracks);
        histos.fill(HIST("hBC_kTVX_nTOFtracks"), collBC, nTOFtracks);
        histos.fill(HIST("hBC_kTVX_nTRDtracks"), collBC, nTRDtracks);

        histos.fill(HIST("hBC_kTVX_nTRDTOFtracks"), collBC, nTRDTOFtracks);
        histos.fill(HIST("hBC_kTVX_nITSTPCTRDtracks"), collBC, nITSTPCTRDtracks);
        histos.fill(HIST("hBC_kTVX_nITSTPCTOFtracks"), collBC, nITSTPCTOFtracks);
        histos.fill(HIST("hBC_kTVX_nITSTPCTRDTOFtracks"), collBC, nITSTPCTRDTOFtracks);
      }

      // TVXinTRD collisions
      if (col.alias_bit(kTVXinTRD)) {
        histos.fill(HIST("hBC_kTVXinTRD_nAllTracks"), collBC, nAllTracks);
        histos.fill(HIST("hBC_kTVXinTRD_nTracksPV"), collBC, nTracksPV);
        // histos.fill(HIST("hBC_kTVXinTRD_nGlobalTracks"), collBC, nGlobalTracks);
        histos.fill(HIST("hBC_kTVXinTRD_nITStracks"), collBC, nITStracks);
        histos.fill(HIST("hBC_kTVXinTRD_nTPCtracks"), collBC, nTPCtracks);
        histos.fill(HIST("hBC_kTVXinTRD_nTOFtracks"), collBC, nTOFtracks);
        histos.fill(HIST("hBC_kTVXinTRD_nTRDtracks"), collBC, nTRDtracks);

        histos.fill(HIST("hBC_kTVXinTRD_nTRDTOFtracks"), collBC, nTRDTOFtracks);
        histos.fill(HIST("hBC_kTVXinTRD_nITSTPCTRDtracks"), collBC, nITSTPCTRDtracks);
        histos.fill(HIST("hBC_kTVXinTRD_nITSTPCTOFtracks"), collBC, nITSTPCTOFtracks);
        histos.fill(HIST("hBC_kTVXinTRD_nITSTPCTRDTOFtracks"), collBC, nITSTPCTRDTOFtracks);
      }

    } // end of collisions loop
  }
  PROCESS_SWITCH(RobustFluctuationBunchCrossingQA, processRobustFluctuationBunchCrossingQA, "Process Run3 event selection", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<RobustFluctuationObservables>(cfgc),
    adaptAnalysisTask<RobustFluctuationBunchCrossingQA>(cfgc)};
}
