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

/// \file flowEventPlane.cxx
/// \brief Flow calculation using event plane.
/// \author Yash Patley <yash.patley@cern.ch>

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <map>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;

enum CorrectionType {
  kFineCorr = 0,
  kCoarseCorr,
  kNCorr
};

enum CollisionParameterType {
  kCent = 0,
  kVx,
  kVy,
  kVz
};

enum ZDCXYType {
  kXa = 0,
  kYa,
  kXc,
  kYc,
  kXYAC
};

struct FlowEventPlane {
  // Configurables
  // Collisions
  Configurable<float> cMinZVtx{"cMinZVtx", -10.0, "Min VtxZ cut"};
  Configurable<float> cMaxZVtx{"cMaxZVtx", 10.0, "Max VtxZ cut"};
  Configurable<float> cMinCent{"cMinCent", 0., "Minumum Centrality"};
  Configurable<float> cMaxCent{"cMaxCent", 100.0, "Maximum Centrality"};
  Configurable<bool> cSel8Trig{"cSel8Trig", true, "Sel8 (T0A + T0C) Selection Run3"};
  Configurable<bool> cTriggerTvxSel{"cTriggerTvxSel", false, "Trigger Time and Vertex Selection"};
  Configurable<bool> cTFBorder{"cTFBorder", false, "Timeframe Border Selection"};
  Configurable<bool> cNoItsROBorder{"cNoItsROBorder", false, "No ITSRO Border Cut"};
  Configurable<bool> cItsTpcVtx{"cItsTpcVtx", false, "ITS+TPC Vertex Selection"};
  Configurable<bool> cPileupReject{"cPileupReject", true, "Pileup rejection"};
  Configurable<bool> cZVtxTimeDiff{"cZVtxTimeDiff", false, "z-vtx time diff selection"};
  Configurable<bool> cIsGoodITSLayers{"cIsGoodITSLayers", true, "Good ITS Layers All"};
  Configurable<float> cMinOccupancy{"cMinOccupancy", 0, "Minimum FT0C Occupancy"};
  Configurable<float> cMaxOccupancy{"cMaxOccupancy", 1e6, "Maximum FT0C Occupancy"};

  // Tracks
  Configurable<float> cTrackMinPt{"cTrackMinPt", 0.15, "p_{T} minimum"};
  Configurable<float> cTrackMaxPt{"cTrackMaxPt", 2.0, "p_{T} maximum"};
  Configurable<float> cTrackEtaCut{"cTrackEtaCut", 0.8, "Pseudorapidity cut"};
  Configurable<bool> cTrackGlobal{"cTrackGlobal", true, "Global Track"};
  Configurable<float> cTrackDcaXYCut{"cTrackDcaXYCut", 0.1, "DcaXY Cut"};
  Configurable<float> cTrackDcaZCut{"cTrackDcaZCut", 1., "DcaXY Cut"};

  // Coarse binning factor
  Configurable<int> cAxisCBF{"cAxisCBF", 1, "Coarse Bin Factor"};

  // Cent Vx Vy Vz Bins
  Configurable<int> cAxisCentBins{"cAxisCentBins", 20, "NBins Centrality"};
  Configurable<int> cAxisVxyBins{"cAxisVxyBins", 20, "NBins Vx Vy"};
  Configurable<int> cAxisVzBins{"cAxisVzBins", 20, "NBins Vz"};
  Configurable<float> cAxisVxMin{"cAxisVxMin", -0.06, "Vx Min"};
  Configurable<float> cAxisVxMax{"cAxisVxMax", -0.02, "Vx Max"};
  Configurable<float> cAxisVyMin{"cAxisVyMin", -0.01, "Vy Min"};
  Configurable<float> cAxisVyMax{"cAxisVyMax", 0.006, "Vy Max"};

  // Corrections
  Configurable<std::vector<int>> cCorrFlagVector{"cCorrFlagVector", {0, 0, 0, 0, 0, 0}, "Correction Flag"};

  // CCDB
  Configurable<std::string> cCcdbUrl{"cCcdbUrl", "http://ccdb-test.cern.ch:8080", "url of ccdb"};
  Configurable<std::string> cCcdbPath{"cCcdbPath", "Users/y/ypatley/DFOO", "Path for ccdb-object"};

  // Initialize CCDB Service
  Service<o2::ccdb::BasicCCDBManager> ccdbService;

  // Histogram registry: an object to hold your histograms
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Global objects
  float cent = 0., mult = 0.;
  float posX = 0., posY = 0., posZ = 0.;

  std::array<float, 4> znaXWeigthEnergy = {1., 1., 1., 1.};
  std::array<float, 4> znaYWeigthEnergy = {1., 1., 1., 1.};
  std::array<float, 4> zncXWeigthEnergy = {1., 1., 1., 1.};
  std::array<float, 4> zncYWeigthEnergy = {1., 1., 1., 1.};

  std::vector<std::vector<std::string>> vCoarseCorrHistNames = {
    {"hXZNAVsCentVxVyVz"},
    {"hYZNAVsCentVxVyVz"},
    {"hXZNCVsCentVxVyVz"},
    {"hYZNCVsCentVxVyVz"}};
  std::vector<std::vector<std::string>> vFineCorrHistNames = {
    {"hXZNAVsCent", "hXZNAVsVx", "hXZNAVsVy", "hXZNAVsVz"},
    {"hYZNAVsCent", "hYZNAVsVx", "hYZNAVsVy", "hYZNAVsVz"},
    {"hXZNCVsCent", "hXZNCVsVx", "hXZNCVsVy", "hXZNCVsVz"},
    {"hYZNCVsCent", "hYZNCVsVx", "hYZNCVsVy", "hYZNCVsVz"}};

  // Map for Correction Type and Histogram Names
  std::map<CorrectionType, std::vector<std::vector<std::string>>> corrTypeHistNameMap = {{kFineCorr, vFineCorrHistNames}, {kCoarseCorr, vCoarseCorrHistNames}};

  // Structure to hold CCDB objects
  struct CcdbObjects {
    TList* ccdbList;
    TObject* obj;
    TProfile* hp;
    THnSparseF* hn;
  } ccdbObjects;

  void init(InitContext const&)
  {
    // Set CCDB url
    ccdbService->setURL(cCcdbUrl.value);
    ccdbService->setCaching(true);

    // Define axes
    const AxisSpec axisZDCEnergy{1000, 0, 5000, "ZD[AC] Signal"};

    const AxisSpec axisCent{100, 0., 100, "FT0C%"};
    const AxisSpec axisVx{cAxisVxyBins, cAxisVxMin, cAxisVxMax, "V_{X}(cm)"};
    const AxisSpec axisVy{cAxisVxyBins, cAxisVyMin, cAxisVyMax, "V_{Y}(cm)"};
    const AxisSpec axisVz{cAxisVzBins, cMinZVtx, cMaxZVtx, "V_{Z}(cm)"};

    const AxisSpec axisCoarseCent{cAxisCentBins / cAxisCBF, cMinCent, cMaxCent, "FT0C%"};
    const AxisSpec axisCoarseVx{cAxisVxyBins / cAxisCBF, cAxisVxMin, cAxisVxMax, "V_{x}"};
    const AxisSpec axisCoarseVy{cAxisVxyBins / cAxisCBF, cAxisVyMin, cAxisVyMax, "V_{y}"};
    const AxisSpec axisCoarseVz{cAxisVzBins / cAxisCBF, cMinZVtx, cMaxZVtx, "V_{z}"};

    const AxisSpec axisFineCent{cAxisCentBins, cMinCent, cMaxCent, "FT0C%"};
    const AxisSpec axisFineVx{cAxisVxyBins, cAxisVxMin, cAxisVxMax, "V_{x}"};
    const AxisSpec axisFineVy{cAxisVxyBins, cAxisVyMin, cAxisVyMax, "V_{x}"};
    const AxisSpec axisFineVz{cAxisVzBins, cMinZVtx, cMaxZVtx, "V_{z}"};

    const AxisSpec axisXa{40, -1, 1, "X^{ZNA}_{1}"};
    const AxisSpec axisYa{40, -1, 1, "Y^{ZNA}_{1}"};
    const AxisSpec axisXc{40, -1, 1, "X^{ZNC}_{1}"};
    const AxisSpec axisYc{40, -1, 1, "Y^{ZNC}_{1}"};

    const AxisSpec axisPsi{18, -PIHalf, PIHalf, "#Psi_{SP}"};

    const AxisSpec axisXYac{600, -6, 6, "Q^{t}Q^{p}"};
    const AxisSpec axisV1{400, -4, 4, "v_{1}"};

    const AxisSpec axisTrackPt{100, 0., 10., "p_{T} (GeV/#it{c})"};
    const AxisSpec axisTrackEta{16, -0.8, 0.8, "#eta"};

    const AxisSpec axisTrackDcaXY{60, -0.15, 0.15, "DCA_{XY}"};
    const AxisSpec axisTrackDcaZ{230, -1.15, 1.15, "DCA_{XY}"};

    // Create histograms
    histos.add("Event/hCent", "FT0C%", kTH1F, {axisCent});
    histos.add("Event/hVx", "V_{x}", kTH1F, {axisVx});
    histos.add("Event/hVy", "V_{y}", kTH1F, {axisVy});
    histos.add("Event/hVz", "V_{z}", kTH1F, {axisVz});
    histos.add("QA/hZNASignalSector1", "ZNA Signal Sector 1", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNASignalSector2", "ZNA Signal Sector 2", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNASignalSector3", "ZNA Signal Sector 3", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNASignalSector4", "ZNA Signal Sector 4", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNCSignalSector1", "ZNC Signal Sector 1", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNCSignalSector2", "ZNC Signal Sector 2", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNCSignalSector3", "ZNC Signal Sector 3", kTH1F, {axisZDCEnergy});
    histos.add("QA/hZNCSignalSector4", "ZNC Signal Sector 4", kTH1F, {axisZDCEnergy});
    histos.add("CorrHist/hWtXZNA", "X^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtYZNA", "Y^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtXZNC", "X^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hWtYZNC", "Y^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtXZNA", "X^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtYZNA", "Y^{ZNA}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtXZNC", "X^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hUWtYZNC", "Y^{ZNC}_{1}", kTHnSparseF, {axisCoarseCent, axisCoarseVx, axisCoarseVy, axisCoarseVz});
    histos.add("CorrHist/hXZNAVsCent", "X^{ZNA}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hXZNAVsVx", "X^{ZNA}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hXZNAVsVy", "X^{ZNA}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hXZNAVsVz", "X^{ZNA}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hYZNAVsCent", "Y^{ZNA}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hYZNAVsVx", "Y^{ZNA}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hYZNAVsVy", "Y^{ZNA}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hYZNAVsVz", "Y^{ZNA}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hXZNCVsCent", "X^{ZNC}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hXZNCVsVx", "X^{ZNC}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hXZNCVsVy", "X^{ZNC}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hXZNCVsVz", "X^{ZNC}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("CorrHist/hYZNCVsCent", "Y^{ZNC}_{1} Vs Cent", kTProfile, {axisFineCent});
    histos.add("CorrHist/hYZNCVsVx", "Y^{ZNC}_{1} Vs V_{x}", kTProfile, {axisFineVx});
    histos.add("CorrHist/hYZNCVsVy", "Y^{ZNC}_{1} Vs V_{y}", kTProfile, {axisFineVy});
    histos.add("CorrHist/hYZNCVsVz", "Y^{ZNC}_{1} Vs V_{z}", kTProfile, {axisFineVz});
    histos.add("Checks/hPsiSPA", "#Psi_{SP}^{A} distribution", kTH2F, {axisCent, axisPsi});
    histos.add("Checks/hPsiSPC", "#Psi_{SP}^{C} distribution", kTH2F, {axisCent, axisPsi});
    histos.add("Checks/hCosPsiSPAC", "Cos(#Psi_{SP}^{A} #minus #Psi_{SP}^{C}) distribution", kTProfile, {axisCent});
    histos.add("Checks/hSinPsiSPAC", "Sin(#Psi_{SP}^{A} #minus #Psi_{SP}^{C}) distribution", kTProfile, {axisCent});
    histos.add("Checks/hXaXc", "X^{ZNC}_{1} Vs X^{ZNA}_{1}", kTProfile, {axisCent});
    histos.add("Checks/hYaYc", "Y^{ZNC}_{1} Vs Y^{ZNA}_{1}", kTProfile, {axisCent});
    histos.add("TrackQA/hPtDcaXY", "DCA_{XY} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaXY});
    histos.add("TrackQA/hPtDcaZ", "DCA_{Z} vs p_{T}", kTH2F, {axisTrackPt, axisTrackDcaZ});
    histos.add("DF/hQaQc", "X^{A}_{1}X^{C}_{1} + Y^{A}_{1}Y^{C}_{1}", kTProfile, {axisCent});
    histos.add("DF/hAQu", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
    histos.add("DF/hCQu", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
    histos.add("DF/hAQuPos", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
    histos.add("DF/hCQuPos", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
    histos.add("DF/hAQuNeg", "u_{x}X^{A}_{1} + u_{y}Y^{A}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
    histos.add("DF/hCQuNeg", "u_{x}X^{C}_{1} + u_{y}Y^{C}_{1}", kTH3F, {axisCent, axisV1, axisTrackEta});
  }

  template <typename C>
  bool selCollision(C const& col)
  {
    if (col.posZ() <= cMinZVtx || col.posZ() >= cMaxZVtx) { // VtxZ selection
      return false;
    }

    if (cSel8Trig && !col.sel8()) { // Sel8 selection
      return false;
    }

    cent = col.centFT0C();
    if (cent <= cMinCent || cent >= cMaxCent) { // Centrality selection
      return false;
    }

    if (col.ft0cOccupancyInTimeRange() < cMinOccupancy || col.ft0cOccupancyInTimeRange() > cMaxOccupancy) { // Occupancy cut
      return false;
    }

    if (cTriggerTvxSel && !col.selection_bit(aod::evsel::kIsTriggerTVX)) { // Time and Vertex trigger
      return false;
    }

    if (cTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder)) { // Time frame border
      return false;
    }

    if (cNoItsROBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder)) { // ITS Readout frame border
      return false;
    }

    if (cItsTpcVtx && !col.selection_bit(aod::evsel::kIsVertexITSTPC)) { // ITS+TPC Vertex
      return false;
    }

    if (cPileupReject && !col.selection_bit(aod::evsel::kNoSameBunchPileup)) { // Pile-up rejection
      return false;
    }

    if (cZVtxTimeDiff && !col.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) { // ZvtxFT0 vs PV
      return false;
    }

    if (cIsGoodITSLayers && !col.selection_bit(aod::evsel::kIsGoodITSLayersAll)) { // All ITS layer active
      return false;
    }

    // Set Multiplicity
    mult = col.multNTracksHasTPC();

    return true;
  }

  // Track Selection
  template <typename T>
  bool selectTrack(T const& track)
  {
    if (track.pt() <= cTrackMinPt || track.pt() >= cTrackMaxPt || std::abs(track.eta()) >= cTrackEtaCut) {
      return false;
    }

    if (cTrackGlobal && !track.isGlobalTrackWoDCA()) {
      return false;
    }

    if (std::abs(track.dcaXY()) > cTrackDcaXYCut || std::abs(track.dcaZ()) > cTrackDcaZCut) {
      return false;
    }

    return true;
  }

  template <typename C, typename F>
  std::vector<float> getAvgCorrFactors(const C& ccdbObject, const F& corrType, const std::vector<float>& vCollParam)
  {
    std::vector<float> vAvgOutput = {0., 0., 0., 0.};
    std::vector<std::vector<std::string>> vHistNames = corrTypeHistNameMap.at(corrType);
    int binarray[4];
    int cntrx = 0;
    for (auto const& x : vHistNames) {
      int cntry = 0;
      for (auto const& y : x) {
        ccdbObjects.obj = reinterpret_cast<TObject*>(ccdbObject->FindObject(y.c_str()));
        if (corrType == kFineCorr) {
          ccdbObjects.hp = reinterpret_cast<TProfile*>(ccdbObjects.obj->Clone());
          vAvgOutput[cntrx] += ccdbObjects.hp->GetBinContent(ccdbObjects.hp->GetXaxis()->FindBin(vCollParam[cntry]));
          delete ccdbObjects.hp;
        } else {
          ccdbObjects.hn = reinterpret_cast<THnSparseF*>(ccdbObjects.obj->Clone());
          for (int i = 0; i < static_cast<int>(vHistNames.size()); ++i) {
            binarray[i] = ccdbObjects.hn->GetAxis(i)->FindBin(vCollParam[i]);
          }
          vAvgOutput[cntrx] += ccdbObjects.hn->GetBinContent(ccdbObjects.hn->GetBin(binarray));
          delete ccdbObjects.hn;
        }
        ++cntry;
      }
      ++cntrx;
    }
    return vAvgOutput;
  }

  void applyCorrection(const std::vector<float> inputParam, const int& runNumber, std::vector<float>& outputParam)
  {
    std::vector<int> vCorrFlags = static_cast<std::vector<int>>(cCorrFlagVector);
    int nitr = vCorrFlags.size();
    CorrectionType corrType = kFineCorr;
    std::string ccdbPath;

    for (int i = 0; i < nitr; ++i) {
      // Don't correct if corrFlag != 1
      if (vCorrFlags[i] != 1) {
        continue;
      }

      // Set correction type
      if (i % kNCorr == 0) {
        corrType = kCoarseCorr;
      } else {
        corrType = kFineCorr;
      }

      // Set ccdb path
      ccdbPath = static_cast<std::string>(cCcdbPath) + "/CorrItr_" + std::to_string(i + 1) + "/Run" + std::to_string(runNumber);

      // Get object from CCDB
      ccdbObjects.ccdbList = ccdbService->getForTimeStamp<TList>(ccdbPath, -1);

      // Check CCDB Object
      if (!ccdbObjects.ccdbList) {
        LOGF(warning, "CCDB OBJECT NOT FOUND");
        return;
      }

      // Get averages
      std::vector<float> vAvg = getAvgCorrFactors(ccdbObjects.ccdbList, corrType, inputParam);

      // Apply correction
      outputParam[kXa] -= vAvg[kXa];
      outputParam[kYa] -= vAvg[kYa];
      outputParam[kXc] -= vAvg[kXc];
      outputParam[kYc] -= vAvg[kYc];
    }
  }

  void fillCorrHist(const std::vector<float>& vCollParam, const std::vector<float>& vSP)
  {
    histos.fill(HIST("CorrHist/hWtXZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kXa]);
    histos.fill(HIST("CorrHist/hWtYZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kYa]);
    histos.fill(HIST("CorrHist/hWtXZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kXc]);
    histos.fill(HIST("CorrHist/hWtYZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz], vSP[kYc]);
    histos.fill(HIST("CorrHist/hUWtXZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtYZNA"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtXZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hUWtYZNC"), vCollParam[kCent], vCollParam[kVx], vCollParam[kVy], vCollParam[kVz]);
    histos.fill(HIST("CorrHist/hXZNAVsCent"), vCollParam[kCent], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVx"), vCollParam[kVx], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVy"), vCollParam[kVy], vSP[kXa]);
    histos.fill(HIST("CorrHist/hXZNAVsVz"), vCollParam[kVz], vSP[kXa]);
    histos.fill(HIST("CorrHist/hYZNAVsCent"), vCollParam[kCent], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVx"), vCollParam[kVx], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVy"), vCollParam[kVy], vSP[kYa]);
    histos.fill(HIST("CorrHist/hYZNAVsVz"), vCollParam[kVz], vSP[kYa]);
    histos.fill(HIST("CorrHist/hXZNCVsCent"), vCollParam[kCent], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVx"), vCollParam[kVx], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVy"), vCollParam[kVy], vSP[kXc]);
    histos.fill(HIST("CorrHist/hXZNCVsVz"), vCollParam[kVz], vSP[kXc]);
    histos.fill(HIST("CorrHist/hYZNCVsCent"), vCollParam[kCent], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVx"), vCollParam[kVx], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVy"), vCollParam[kVy], vSP[kYc]);
    histos.fill(HIST("CorrHist/hYZNCVsVz"), vCollParam[kVz], vSP[kYc]);
  }

  template <typename T>
  void fillTrackHist(const T& track)
  {
    histos.fill(HIST("TrackQA/hPtDcaZ"), track.pt(), track.dcaZ());
    histos.fill(HIST("TrackQA/hPtDcaXY"), track.pt(), track.dcaXY());
  }

  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;
  using CollisionsRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::MultsExtra>;
  using Tracks = soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr, aod::TrackCompColls>;

  void process(CollisionsRun3::iterator const& collision, BCsRun3 const& /* bcs*/, aod::Zdcs const&, Tracks const& tracks)
  {
    // Event selection
    if (!selCollision(collision)) {
      return;
    }
    posX = collision.posX();
    posY = collision.posY();
    posZ = collision.posZ();
    std::vector<float> vCollParam = {cent, posX, posY, posZ};

    histos.fill(HIST("Event/hCent"), cent);
    histos.fill(HIST("Event/hVx"), posX);
    histos.fill(HIST("Event/hVy"), posY);
    histos.fill(HIST("Event/hVz"), posZ);

    // Get bunch crossing
    auto bc = collision.foundBC_as<BCsRun3>();

    // check zdc
    if (!bc.has_zdc()) {
      return;
    }

    auto zdc = bc.zdc();
    auto znaEnergy = zdc.energySectorZNA();
    auto zncEnergy = zdc.energySectorZNC();
    auto znaEnergyCommon = zdc.energyCommonZNA();
    auto zncEnergyCommon = zdc.energyCommonZNC();

    // check energy deposits
    if (znaEnergyCommon <= 0 || zncEnergyCommon <= 0 || znaEnergy[0] <= 0 || znaEnergy[1] <= 0 || znaEnergy[2] <= 0 || znaEnergy[3] <= 0 || zncEnergy[0] <= 0 || zncEnergy[1] <= 0 || zncEnergy[2] <= 0 || zncEnergy[3] <= 0) {
      return;
    }

    // Fill QA histograms
    histos.fill(HIST("QA/hZNASignalSector1"), znaEnergy[0]);
    histos.fill(HIST("QA/hZNASignalSector2"), znaEnergy[1]);
    histos.fill(HIST("QA/hZNASignalSector3"), znaEnergy[2]);
    histos.fill(HIST("QA/hZNASignalSector4"), znaEnergy[3]);
    histos.fill(HIST("QA/hZNCSignalSector1"), zncEnergy[0]);
    histos.fill(HIST("QA/hZNCSignalSector2"), zncEnergy[1]);
    histos.fill(HIST("QA/hZNCSignalSector3"), zncEnergy[2]);
    histos.fill(HIST("QA/hZNCSignalSector4"), zncEnergy[3]);

    /*auto alphaZDC = 0.395;*/
    const double x[4] = {-1.75, 1.75, -1.75, 1.75};
    const double y[4] = {-1.75, -1.75, 1.75, 1.75};

    // Calculate X and Y
    float znaXSumNum = 0., znaXSumDnm = 0.;
    float znaYSumNum = 0., znaYSumDnm = 0.;
    float zncXSumNum = 0., zncXSumDnm = 0.;
    float zncYSumNum = 0., zncYSumDnm = 0.;

    // Loop over zdc sectors
    for (int i = 0; i < kXYAC; ++i) {
      znaXSumNum += znaXWeigthEnergy[i] * znaEnergy[i] * x[i];
      znaYSumNum += znaYWeigthEnergy[i] * znaEnergy[i] * y[i];
      znaXSumDnm += znaXWeigthEnergy[i] * znaEnergy[i];
      znaYSumDnm += znaYWeigthEnergy[i] * znaEnergy[i];
      zncXSumNum += zncXWeigthEnergy[i] * zncEnergy[i] * x[i];
      zncYSumNum += zncYWeigthEnergy[i] * zncEnergy[i] * y[i];
      zncXSumDnm += zncXWeigthEnergy[i] * zncEnergy[i];
      zncYSumDnm += zncYWeigthEnergy[i] * zncEnergy[i];
    }

    // Get X and Y for A and C side ZNA
    std::vector<float> vSP = {0, 0, 0, 0};
    vSP[kXa] = znaXSumNum / znaXSumDnm;
    vSP[kYa] = znaYSumNum / znaYSumDnm;
    vSP[kXc] = zncXSumNum / zncXSumDnm;
    vSP[kYc] = zncYSumNum / zncYSumDnm;

    // Do corrections
    int runNumber = collision.foundBC_as<BCsRun3>().runNumber();
    applyCorrection(vCollParam, runNumber, vSP);

    // Fill X and Y histograms
    fillCorrHist(vCollParam, vSP);
    float psiA = std::atan2(vSP[kYa], vSP[kXa]);
    float psiC = std::atan2(vSP[kYa], vSP[kXa]);
    histos.fill(HIST("Checks/hXaXc"), cent, (vSP[kXa] * vSP[kXc]));
    histos.fill(HIST("Checks/hYaYc"), cent, (vSP[kYa] * vSP[kYc]));
    histos.fill(HIST("Checks/hPsiSPA"), cent, psiA);
    histos.fill(HIST("Checks/hPsiSPC"), cent, psiC);
    histos.fill(HIST("Checks/hCosPsiSPAC"), cent, std::cos(psiA - psiC));
    histos.fill(HIST("Checks/hSinPsiSPAC"), cent, std::sin(psiA - psiC));

    // Directed flow
    float qac = (vSP[kXa] * vSP[kXc]) + (vSP[kYa] * vSP[kYc]);
    histos.fill(HIST("DF/hQaQc"), cent, qac);

    // Loop over tracks
    float ux = 0., uy = 0., v1a = 0., v1c = 0.;
    for (auto const& track : tracks) {
      // Select track
      if (!selectTrack(track)) {
        continue;
      }

      // Fill track QA
      fillTrackHist(track);

      // Get directed flow
      ux = std::cos(track.phi());
      uy = std::sin(track.phi());
      v1a = ux * vSP[kXa] + uy * vSP[kYa];
      v1c = ux * vSP[kXc] + uy * vSP[kYc];

      // Fill histogram
      histos.fill(HIST("DF/hAQu"), cent, v1a, track.eta());
      histos.fill(HIST("DF/hCQu"), cent, v1c, track.eta());
      if (track.sign() > 0) {
        histos.fill(HIST("DF/hAQuPos"), cent, v1a, track.eta());
        histos.fill(HIST("DF/hCQuPos"), cent, v1c, track.eta());
      } else {
        histos.fill(HIST("DF/hAQuNeg"), cent, v1a, track.eta());
        histos.fill(HIST("DF/hCQuNeg"), cent, v1c, track.eta());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEventPlane>(cfgc)};
}
