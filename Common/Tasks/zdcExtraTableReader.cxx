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

/// \file zdcExtraTableReader.cxx
/// \brief Task reading AOD/ZDCEXTRA table
/// \author Uliana Dmitrieva <uliana.dmitrieva@cern.ch>, INFN Torino

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/ZDCExtra.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <TList.h>
#include <TProfile3D.h>

#include <Rtypes.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

namespace
{

std::unordered_map<int, TH1*> gEventCounter;
std::unordered_map<int, TH2*> gCentroidZNA;
std::unordered_map<int, TH2*> gCentroidZNC;
std::unordered_map<int, TH1*> gPmcZNA;
std::unordered_map<int, TH1*> gPm1ZNA;
std::unordered_map<int, TH1*> gPm2ZNA;
std::unordered_map<int, TH1*> gPm3ZNA;
std::unordered_map<int, TH1*> gPm4ZNA;
std::unordered_map<int, TH1*> gSumZNA;
std::unordered_map<int, TH1*> gPmcZNC;
std::unordered_map<int, TH1*> gPm1ZNC;
std::unordered_map<int, TH1*> gPm2ZNC;
std::unordered_map<int, TH1*> gPm3ZNC;
std::unordered_map<int, TH1*> gPm4ZNC;
std::unordered_map<int, TH1*> gSumZNC;
std::unordered_map<int, TH2*> gQxVsCentZNA;
std::unordered_map<int, TH2*> gQyVsCentZNA;
std::unordered_map<int, TH2*> gQxVsVxZNA;
std::unordered_map<int, TH2*> gQyVsVxZNA;
std::unordered_map<int, TH2*> gQxVsVyZNA;
std::unordered_map<int, TH2*> gQyVsVyZNA;
std::unordered_map<int, TH2*> gQxVsVzZNA;
std::unordered_map<int, TH2*> gQyVsVzZNA;
std::unordered_map<int, TH2*> gQxVsCentZNC;
std::unordered_map<int, TH2*> gQyVsCentZNC;
std::unordered_map<int, TH2*> gQxVsVxZNC;
std::unordered_map<int, TH2*> gQyVsVxZNC;
std::unordered_map<int, TH2*> gQxVsVyZNC;
std::unordered_map<int, TH2*> gQyVsVyZNC;
std::unordered_map<int, TH2*> gQxVsVzZNC;
std::unordered_map<int, TH2*> gQyVsVzZNC;
std::unordered_map<int, TH2*> gQxQyVsCent;
std::unordered_map<int, TH2*> gQyQxVsCent;
std::unordered_map<int, TH2*> gQxQxVsCent;
std::unordered_map<int, TH2*> gQyQyVsCent;

std::unordered_map<int, THn*> gQx5DZNA;
std::unordered_map<int, THn*> gQy5DZNA;
std::unordered_map<int, THn*> gQx5DZNC;
std::unordered_map<int, THn*> gQy5DZNC;

std::unordered_map<int, TH1*> gPsiZNA;
std::unordered_map<int, TH1*> gPsiZNC;

std::unordered_map<int, TH1*> gVx;
std::unordered_map<int, TH1*> gVy;

// centroid stability vs timestamp
std::unordered_map<int, TH2*> gQxVsTimeZNA;
std::unordered_map<int, TH2*> gQyVsTimeZNA;
std::unordered_map<int, TH2*> gQxVsTimeZNC;
std::unordered_map<int, TH2*> gQyVsTimeZNC;

std::unordered_map<int, TProfile3D*> gShiftProfileZNA;
std::unordered_map<int, TProfile3D*> gShiftProfileZNC;

TH1* gCurrentEventCounter{nullptr};
TH2* gCurrentCentroidZNA{nullptr};
TH2* gCurrentCentroidZNC{nullptr};
TH1* gCurrentPmcZNA{nullptr};
TH1* gCurrentPm1ZNA{nullptr};
TH1* gCurrentPm2ZNA{nullptr};
TH1* gCurrentPm3ZNA{nullptr};
TH1* gCurrentPm4ZNA{nullptr};
TH1* gCurrentSumZNA{nullptr};
TH1* gCurrentPmcZNC{nullptr};
TH1* gCurrentPm1ZNC{nullptr};
TH1* gCurrentPm2ZNC{nullptr};
TH1* gCurrentPm3ZNC{nullptr};
TH1* gCurrentPm4ZNC{nullptr};
TH1* gCurrentSumZNC{nullptr};
TH2* gCurrentQxVsCentZNA{nullptr};
TH2* gCurrentQyVsCentZNA{nullptr};
TH2* gCurrentQxVsVxZNA{nullptr};
TH2* gCurrentQyVsVxZNA{nullptr};
TH2* gCurrentQxVsVyZNA{nullptr};
TH2* gCurrentQyVsVyZNA{nullptr};
TH2* gCurrentQxVsVzZNA{nullptr};
TH2* gCurrentQyVsVzZNA{nullptr};
TH2* gCurrentQxVsCentZNC{nullptr};
TH2* gCurrentQyVsCentZNC{nullptr};
TH2* gCurrentQxVsVxZNC{nullptr};
TH2* gCurrentQyVsVxZNC{nullptr};
TH2* gCurrentQxVsVyZNC{nullptr};
TH2* gCurrentQyVsVyZNC{nullptr};
TH2* gCurrentQxVsVzZNC{nullptr};
TH2* gCurrentQyVsVzZNC{nullptr};
TH2* gCurrentQxQyVsCent{nullptr};
TH2* gCurrentQyQxVsCent{nullptr};
TH2* gCurrentQxQxVsCent{nullptr};
TH2* gCurrentQyQyVsCent{nullptr};

THn* gCurrentQxZNA{nullptr};
THn* gCurrentQyZNA{nullptr};
THn* gCurrentQxZNC{nullptr};
THn* gCurrentQyZNC{nullptr};

TH1* gCurrentPsiZNA{nullptr} ;
TH1* gCurrentPsiZNC{nullptr};

TH1* gCurrentVx{nullptr};
TH1* gCurrentVy{nullptr};

TH2* gCurrentQxVsTimeZNA{nullptr};
TH2* gCurrentQyVsTimeZNA{nullptr};
TH2* gCurrentQxVsTimeZNC{nullptr};
TH2* gCurrentQyVsTimeZNC{nullptr};

TProfile3D* gCurrentShiftProfileZNA{nullptr};
TProfile3D* gCurrentShiftProfileZNC{nullptr};

} // namespace

// Helper for 4D recentering maps
double getMeanQFromMap(THn* h, double cent, double vx, double vy, double vz)
{
  if (!h) {
    LOGF(fatal, "[MeanQ] Null THn pointer");
  }

  TAxis* axCent = h->GetAxis(0);
  TAxis* axVx = h->GetAxis(1);
  TAxis* axVy = h->GetAxis(2);
  TAxis* axVz = h->GetAxis(3);

  if (!axCent || !axVx || !axVy || !axVz) {
    LOGF(fatal, "[MeanQ] One of THn axes is null");
  }

  int binCent = axCent->FindFixBin(cent);
  int binVx = axVx->FindFixBin(vx);
  int binVy = axVy->FindFixBin(vy);
  int binVz = axVz->FindFixBin(vz);

  int idx[4] = {binCent, binVx, binVy, binVz};
  return h->GetBinContent(idx);
}

// Helper for 1D recentering maps: returns mean Q for coordinate x
// If bin out of range, returns 0.0
double getMeanQ1D(TH1* h, double x)
{
  if (!h) {
    LOGF(fatal, "[MeanQ1D] Null TH1 pointer");
  }
  int bin = h->FindFixBin(x);
  if (bin < 1 || bin > h->GetNbinsX()) {
    return 0.0;
  }
  return h->GetBinContent(bin);
}

struct ZdcExtraTableReader {

  Configurable<int> nBinsZN{"nBinsZN", 2000, "n bins for ZN histograms"};
  Configurable<float> maxZN{"maxZN", 399.5, "Max ZN signal"};
  Configurable<bool> applyTdcCut{"applyTdcCut", true, "Flag for TDC cut"};
  Configurable<float> tdcZnMin{"tdcZnMin", -2.5, "Min ZN TDC cut"};
  Configurable<float> tdcZnMax{"tdcZnMax", 2.5, "Max ZN TDC cut"};
  Configurable<bool> plotPMs{"plotPMs", false, "Flag to plot individual PMs"};

  ConfigurableAxis qxyAxis{"qxyAxis", {100, -2.0f, 2.0f}, ""};

  Configurable<int> vxNbins{"vxNbins", 50, "Bins in Vx"};
  Configurable<float> vxMin{"vxMin", -0.1f, "Vx lower edge"};
  Configurable<float> vxMax{"vxMax", 0.1f, "Vx upper edge"};

  Configurable<int> vyNbins{"vyNbins", 50, "Bins in Vy"};
  Configurable<float> vyMin{"vyMin", -0.1f, "Vy lower edge"};
  Configurable<float> vyMax{"vyMax", 0.1f, "Vy upper edge"};

  Configurable<int> vzNbins{"vzNbins", 50, "Bins in Vz"};
  Configurable<float> vzMin{"vzMin", -10.0f, "Vz lower edge"};
  Configurable<float> vzMax{"vzMax", 10.0f, "Vz upper edge"};

  Configurable<int> centNbins{"centNbins", 16, "Bins in centrality"};
  Configurable<float> centMin{"centMin", 0.0f, "Centrality lower edge"};
  Configurable<float> centMax{"centMax", 80.0f, "Centrality upper edge"};

  Configurable<int> phiNbins{"phiNbins", 60, "Bins in phi"};

  Configurable<int> minNTowersFired{"minNTowersFired", 2, "Minimum number of towers fired for Q-vector determination"};

  Configurable<int> qNbins5D{"qNbins5D", 4, "Bins in each dimension for 5D histograms"};
  Configurable<bool> plot5D{"plot5D", false, "Flag to plot 5D histograms"};

  Configurable<int> calibrationStep{"calibrationStep", 1, "Calibration step"};
  Configurable<bool> isFineCalibrationStep{"isFineCalibrationStep", false, "Calibration: base  or refine"};

  Configurable<bool> applyBeamSpotCorrection{"applyBeamSpotCorrection", true, "Beam spot correction"};

  Configurable<bool> applySel8{"applySel8", true, "Event selection: Sel8"};
  Configurable<bool> applyZVtxCut{"applyZVtxCut", true, "Event selection: zVtx cut set in producer (tipically < 10 cm)"};
  Configurable<bool> applyOccupancyCut{"applyOccupancyCut", false, "Event selection: occupancy cut set in producer"};
  Configurable<bool> selectNoSameBunchPileupEvents{"selectNoSameBunchPileupEvents", false, "Event selection: no same bunch pileup"};
  Configurable<bool> selectGoodZvtxFT0vsPV{"selectGoodZvtxFT0vsPV", false, "Event selection: good Zvtx FT0 vs PV"};
  Configurable<bool> applyNoCollInTimeRangeStandard{"applyNoCollInTimeRangeStandard", false, "Event selection: no collision in time range standard"};
  Configurable<bool> selectVertexITSTPC{"selectVertexITSTPC", false, "Event selection: vertex ITS TPC"};
  Configurable<bool> selectGoodITSLayersAll{"selectGoodITSLayersAll", false, "Event selection: good ITS layers all"};

  Configurable<bool> applyShiftCorrection{"applyShiftCorrection", false, "Apply shift correction (Read from CCDB)"};
  Configurable<bool> fillShiftHistos{"fillShiftHistos", false, "Fill shift profiles (Write to output)"};
  Configurable<int> nHarmonics{"nHarmonics", 10, "Number of harmonics"};

  Configurable<std::string> qRecenteringCcdb{"qRecenteringCcdb", "Users/u/udmitrie/ZDC/LHC24ar_apass2", "Recentering maps containing step folder"};

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Struct to hold calibration data for a single step
  struct CalibStepData {
    // 5D maps (Base)
    THn* hMeanQxZNA{nullptr};
    THn* hMeanQyZNA{nullptr};
    THn* hMeanQxZNC{nullptr};
    THn* hMeanQyZNC{nullptr};

    // 1D maps (Refine)
    TH1* hMeanQxCentZNA{nullptr};
    TH1* hMeanQyCentZNA{nullptr};
    TH1* hMeanQxCentZNC{nullptr};
    TH1* hMeanQyCentZNC{nullptr};

    TH1* hMeanQxVzZNA{nullptr};
    TH1* hMeanQyVzZNA{nullptr};
    TH1* hMeanQxVzZNC{nullptr};
    TH1* hMeanQyVzZNC{nullptr};

    TH1* hMeanQxVxZNA{nullptr};
    TH1* hMeanQyVxZNA{nullptr};
    TH1* hMeanQxVxZNC{nullptr};
    TH1* hMeanQyVxZNC{nullptr};

    TH1* hMeanQxVyZNA{nullptr};
    TH1* hMeanQyVyZNA{nullptr};
    TH1* hMeanQxVyZNC{nullptr};
    TH1* hMeanQyVyZNC{nullptr};

    // Destructor to handle cleanup automatically
    ~CalibStepData()
    {
      delete hMeanQxZNA;
      delete hMeanQyZNA;
      delete hMeanQxZNC;
      delete hMeanQyZNC;

      delete hMeanQxCentZNA;
      delete hMeanQyCentZNA;
      delete hMeanQxCentZNC;
      delete hMeanQyCentZNC;

      delete hMeanQxVzZNA;
      delete hMeanQyVzZNA;
      delete hMeanQxVzZNC;
      delete hMeanQyVzZNC;

      delete hMeanQxVxZNA;
      delete hMeanQyVxZNA;
      delete hMeanQxVxZNC;
      delete hMeanQyVxZNC;

      delete hMeanQxVyZNA;
      delete hMeanQyVyZNA;
      delete hMeanQxVyZNC;
      delete hMeanQyVyZNC;
    }
  };

  // Cache container: Vector index = Step index (0-based, so step 1 is at index 0)
  std::vector<CalibStepData> calibCache;

  // Vertex correction cache
  TH1* hMeanVx{nullptr};
  TH1* hMeanVy{nullptr};

  // Phase shift correction cache
  TProfile3D* shiftProfileZNA{nullptr};
  TProfile3D* shiftProfileZNC{nullptr};

  HistogramRegistry histos{"histos"};

  enum EvSelBits { // same bits as in zdcExtraTableProducer.cxx
    ZVtxCut,
    Sel8,
    OccupancyCut,
    NoSameBunchPileup,
    IsGoodZvtxFT0vsPV,
    NoCollInTimeRangeStandard,
    IsVertexITSTPC,
    IsGoodITSLayersAll,
    AllEvents,
    NEventSelections
  };

  int currentRunNumber{-1};

  // Helper to safely clone a histogram and detach from file
  template <typename T>
  T* safeClone(TObject* obj)
  {
    if (!obj)
      return nullptr;
    T* cloned = dynamic_cast<T*>(obj->Clone());
    if (cloned) {

      if (dynamic_cast<TH1*>(cloned)) {
        dynamic_cast<TH1*>(cloned)->SetDirectory(nullptr);
      }
    }
    return cloned;
  }

  void clearCache()
  {
    delete hMeanVx;
    hMeanVx = nullptr;

    delete hMeanVy;
    hMeanVy = nullptr;

    delete shiftProfileZNA;
    shiftProfileZNA = nullptr;

    delete shiftProfileZNC;
    shiftProfileZNC = nullptr;

    calibCache.clear();
  }

  void initHistos(const int runNumber)
  {
    if (runNumber == currentRunNumber) {
      return;
    }
    currentRunNumber = runNumber;

    if (!gEventCounter.contains(runNumber)) {
      // if new run, initialize histograms

      const AxisSpec axisCounter{1, 0, +1, ""};
      const AxisSpec axisZN{nBinsZN, -0.5, maxZN, "(a.u.)"};
      const AxisSpec axisCent = {centNbins, centMin, centMax, "Centrality (%)"};
      const AxisSpec axisVx = {vxNbins, vxMin, vxMax, "V_{x} (cm)"};
      const AxisSpec axisVy = {vyNbins, vyMin, vyMax, "V_{y} (cm)"};
      const AxisSpec axisVz = {vzNbins, vzMin, vzMax, "V_{z} (cm)"};

      const AxisSpec axisCent5D = {qNbins5D, centMin, centMax, "Centrality (%)"};
      const AxisSpec axisVx5D = {qNbins5D, vxMin, vxMax, "V_{x} (cm)"};
      const AxisSpec axisVy5D = {qNbins5D, vyMin, vyMax, "V_{y} (cm)"};
      const AxisSpec axisVz5D = {qNbins5D, vzMin, vzMax, "V_{z} (cm)"};

      const AxisSpec axisQx{qxyAxis, "Q_{x}"};
      const AxisSpec axisQy{qxyAxis, "Q_{y}"};
      const AxisSpec axisQxQy = {qxyAxis, ""};

      const AxisSpec axisPhi = {phiNbins, -1.0f * o2::constants::math::PI, 1.0f * o2::constants::math::PI, "#phi"};

      const AxisSpec axisTime = {90, 0, 90, "Time (minutes)"}; // 90 minutes

      gEventCounter[runNumber] = histos.add<TH1>(Form("%i/eventCounter", runNumber), "Number of Event; ; #Events Passed Cut", kTH1D, {{NEventSelections, 0, NEventSelections}}).get();
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(AllEvents + 1, "allEvents");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(ZVtxCut + 1, "zVtxCut");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(Sel8 + 1, "Sel8");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(OccupancyCut + 1, "occupancyCut");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(NoSameBunchPileup + 1, "NoSameBunchPileup");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(IsGoodZvtxFT0vsPV + 1, "isGoodZvtxFT0vsPV");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(NoCollInTimeRangeStandard + 1, "noCollInTimeRangeStandard");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(IsVertexITSTPC + 1, "isVertexITSTPC");
      gEventCounter[runNumber]->GetXaxis()->SetBinLabel(IsGoodITSLayersAll + 1, "isGoodITSLayersAll");

      gCentroidZNA[runNumber] = histos.add<TH2>(Form("%i/CentroidZNA", runNumber), "ZNA Centroid; Q_{X}; Q_{Y}", kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}).get();
      gCentroidZNC[runNumber] = histos.add<TH2>(Form("%i/CentroidZNC", runNumber), "ZNC Centroid; Q_{X}; Q_{Y}", kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}).get();
      gPmcZNA[runNumber] = histos.add<TH1>(Form("%i/pmcZNA", runNumber), "; E_{PMC}^{ZNA} (TeV);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm1ZNA[runNumber] = histos.add<TH1>(Form("%i/pm1ZNA", runNumber), "; E_{PM1}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm2ZNA[runNumber] = histos.add<TH1>(Form("%i/pm2ZNA", runNumber), "; E_{PM2}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm3ZNA[runNumber] = histos.add<TH1>(Form("%i/pm3ZNA", runNumber), "; E_{PM3}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm4ZNA[runNumber] = histos.add<TH1>(Form("%i/pm4ZNA", runNumber), "; E_{PM4}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gSumZNA[runNumber] = histos.add<TH1>(Form("%i/sumZNA", runNumber), "; E_{sum PMs}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPmcZNC[runNumber] = histos.add<TH1>(Form("%i/pmcZNC", runNumber), "; E_{PMC}^{ZNC} (TeV);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm1ZNC[runNumber] = histos.add<TH1>(Form("%i/pm1ZNC", runNumber), "; E_{PM1}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm2ZNC[runNumber] = histos.add<TH1>(Form("%i/pm2ZNC", runNumber), "; E_{PM2}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm3ZNC[runNumber] = histos.add<TH1>(Form("%i/pm3ZNC", runNumber), "; E_{PM3}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm4ZNC[runNumber] = histos.add<TH1>(Form("%i/pm4ZNC", runNumber), "; E_{PM4}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gSumZNC[runNumber] = histos.add<TH1>(Form("%i/sumZNC", runNumber), "; E_{sum PMs}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gQxVsCentZNA[runNumber] = histos.add<TH2>(Form("%i/QxVsCentZNA", runNumber), "Q_{x}^{ZNA} vs Centrality", kTH2F, {axisCent, axisQx}).get();
      gQyVsCentZNA[runNumber] = histos.add<TH2>(Form("%i/QyVsCentZNA", runNumber), "Q_{y}^{ZNA} vs Centrality", kTH2F, {axisCent, axisQy}).get();
      gQxVsVxZNA[runNumber] = histos.add<TH2>(Form("%i/QxVsVxZNA", runNumber), "Q_{x}^{ZNA} vs V_{x}; V_{x} (cm); Q_{x}", kTH2F, {axisVx, axisQx}).get();
      gQyVsVxZNA[runNumber] = histos.add<TH2>(Form("%i/QyVsVxZNA", runNumber), "Q_{y}^{ZNA} vs V_{x}; V_{x} (cm); Q_{y}", kTH2F, {axisVx, axisQy}).get();
      gQxVsVyZNA[runNumber] = histos.add<TH2>(Form("%i/QxVsVyZNA", runNumber), "Q_{x}^{ZNA} vs V_{y}; V_{y} (cm); Q_{x}", kTH2F, {axisVy, axisQx}).get();
      gQyVsVyZNA[runNumber] = histos.add<TH2>(Form("%i/QyVsVyZNA", runNumber), "Q_{y}^{ZNA} vs V_{y}; V_{y} (cm); Q_{y}", kTH2F, {axisVy, axisQy}).get();
      gQxVsVzZNA[runNumber] = histos.add<TH2>(Form("%i/QxVsVzZNA", runNumber), "Q_{x}^{ZNA} vs V_{z}; V_{z} (cm); Q_{x}", kTH2F, {axisVz, axisQx}).get();
      gQyVsVzZNA[runNumber] = histos.add<TH2>(Form("%i/QyVsVzZNA", runNumber), "Q_{y}^{ZNA} vs V_{z}; V_{z} (cm); Q_{y}", kTH2F, {axisVz, axisQy}).get();
      gQxVsCentZNC[runNumber] = histos.add<TH2>(Form("%i/QxVsCentZNC", runNumber), "Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{x}", kTH2F, {axisCent, axisQx}).get();
      gQyVsCentZNC[runNumber] = histos.add<TH2>(Form("%i/QyVsCentZNC", runNumber), "Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{y}", kTH2F, {axisCent, axisQy}).get();
      gQxVsVxZNC[runNumber] = histos.add<TH2>(Form("%i/QxVsVxZNC", runNumber), "Q_{x}^{ZNC} vs V_{x}; V_{x} (cm); Q_{x}", kTH2F, {axisVx, axisQx}).get();
      gQyVsVxZNC[runNumber] = histos.add<TH2>(Form("%i/QyVsVxZNC", runNumber), "Q_{y}^{ZNC} vs V_{x}; V_{x} (cm); Q_{y}", kTH2F, {axisVx, axisQy}).get();
      gQxVsVyZNC[runNumber] = histos.add<TH2>(Form("%i/QxVsVyZNC", runNumber), "Q_{x}^{ZNC} vs V_{y}; V_{y} (cm); Q_{x}", kTH2F, {axisVy, axisQx}).get();
      gQyVsVyZNC[runNumber] = histos.add<TH2>(Form("%i/QyVsVyZNC", runNumber), "Q_{y}^{ZNC} vs V_{y}; V_{y} (cm); Q_{y}", kTH2F, {axisVy, axisQy}).get();
      gQxVsVzZNC[runNumber] = histos.add<TH2>(Form("%i/QxVsVzZNC", runNumber), "Q_{x}^{ZNC} vs V_{z}; V_{z} (cm); Q_{x}", kTH2F, {axisVz, axisQx}).get();
      gQyVsVzZNC[runNumber] = histos.add<TH2>(Form("%i/QyVsVzZNC", runNumber), "Q_{y}^{ZNC} vs V_{z}; V_{z} (cm); Q_{y}", kTH2F, {axisVz, axisQy}).get();
      gQxQyVsCent[runNumber] = histos.add<TH2>(Form("%i/QxQyVsCent", runNumber), "Q_{x}^{ZNC}Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{x}^{ZNA}Q_{y}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQyQxVsCent[runNumber] = histos.add<TH2>(Form("%i/QyQxVsCent", runNumber), "Q_{y}^{ZNC}Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{y}^{ZNA}Q_{x}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQxQxVsCent[runNumber] = histos.add<TH2>(Form("%i/QxQxVsCent", runNumber), "Q_{x}^{ZNC}Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{x}^{ZNA}Q_{x}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQyQyVsCent[runNumber] = histos.add<TH2>(Form("%i/QyQyVsCent", runNumber), "Q_{y}^{ZNC}Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{y}^{ZNA}Q_{y}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();

      gQx5DZNA[runNumber] = histos.add<THn>(Form("%i/Qx5DZNA", runNumber), "Qx recenter map ZNA", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQx}, true).get();
      gQy5DZNA[runNumber] = histos.add<THn>(Form("%i/Qy5DZNA", runNumber), "Qy recenter map ZNA", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQy}, true).get();
      gQx5DZNC[runNumber] = histos.add<THn>(Form("%i/Qx5DZNC", runNumber), "Qx recenter map ZNC", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQx}, true).get();
      gQy5DZNC[runNumber] = histos.add<THn>(Form("%i/Qy5DZNC", runNumber), "Qy recenter map ZNC", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQy}, true).get();

      gPsiZNA[runNumber] = histos.add<TH1>(Form("%i/PsiZNA", runNumber), ";#Phi_{ZNA} (rad)", kTH1F, {axisPhi}).get();
      gPsiZNC[runNumber] = histos.add<TH1>(Form("%i/PsiZNC", runNumber), ";#Phi_{ZNC} (rad)", kTH1F, {axisPhi}).get();

      gVx[runNumber] = histos.add<TH1>(Form("%i/Vx", runNumber), "V_{x} distribution; V_{x} (cm); Entries", kTH1F, {axisVx}).get();
      gVy[runNumber] = histos.add<TH1>(Form("%i/Vy", runNumber), "V_{y} distribution; V_{y} (cm); Entries", kTH1F, {axisVy}).get();

      gQxVsTimeZNA[runNumber] = histos.add<TH2>(Form("%i/QxVsTimeZNA", runNumber), "Q_{x}^{ZNA} vs Time; Time (minutes); Q_{x}", kTH2F, {axisTime, axisQx}).get();
      gQyVsTimeZNA[runNumber] = histos.add<TH2>(Form("%i/QyVsTimeZNA", runNumber), "Q_{y}^{ZNA} vs Time; Time (minutes); Q_{y}", kTH2F, {axisTime, axisQy}).get();
      gQxVsTimeZNC[runNumber] = histos.add<TH2>(Form("%i/QxVsTimeZNC", runNumber), "Q_{x}^{ZNC} vs Time; Time (minutes); Q_{x}", kTH2F, {axisTime, axisQx}).get();
      gQyVsTimeZNC[runNumber] = histos.add<TH2>(Form("%i/QyVsTimeZNC", runNumber), "Q_{y}^{ZNC} vs Time; Time (minutes); Q_{y}", kTH2F, {axisTime, axisQy}).get();

      gShiftProfileZNA[runNumber] = histos.add<TProfile3D>(Form("%i/ShiftProfileZNA", runNumber), "ZNA Shift Coeffs;Cent;Type;Harmonic", kTProfile3D, {axisCent, {2, 0, 2}, {nHarmonics, 0, static_cast<double>(nHarmonics)}}).get();
      gShiftProfileZNC[runNumber] = histos.add<TProfile3D>(Form("%i/ShiftProfileZNC", runNumber), "ZNC Shift Coeffs;Cent;Type;Harmonic", kTProfile3D, {axisCent, {2, 0, 2}, {nHarmonics, 0, static_cast<double>(nHarmonics)}}).get();
    }

    gCurrentEventCounter = gEventCounter[currentRunNumber];
    gCurrentCentroidZNA = gCentroidZNA[currentRunNumber];
    gCurrentCentroidZNC = gCentroidZNC[currentRunNumber];
    gCurrentPmcZNA = gPmcZNA[currentRunNumber];
    gCurrentPm1ZNA = gPm1ZNA[currentRunNumber];
    gCurrentPm2ZNA = gPm2ZNA[currentRunNumber];
    gCurrentPm3ZNA = gPm3ZNA[currentRunNumber];
    gCurrentPm4ZNA = gPm4ZNA[currentRunNumber];
    gCurrentSumZNA = gSumZNA[currentRunNumber];
    gCurrentPmcZNC = gPmcZNC[currentRunNumber];
    gCurrentPm1ZNC = gPm1ZNC[currentRunNumber];
    gCurrentPm2ZNC = gPm2ZNC[currentRunNumber];
    gCurrentPm3ZNC = gPm3ZNC[currentRunNumber];
    gCurrentPm4ZNC = gPm4ZNC[currentRunNumber];
    gCurrentSumZNC = gSumZNC[currentRunNumber];
    gCurrentQxVsCentZNA = gQxVsCentZNA[currentRunNumber];
    gCurrentQyVsCentZNA = gQyVsCentZNA[currentRunNumber];
    gCurrentQxVsVxZNA = gQxVsVxZNA[currentRunNumber];
    gCurrentQyVsVxZNA = gQyVsVxZNA[currentRunNumber];
    gCurrentQxVsVyZNA = gQxVsVyZNA[currentRunNumber];
    gCurrentQyVsVyZNA = gQyVsVyZNA[currentRunNumber];
    gCurrentQxVsVzZNA = gQxVsVzZNA[currentRunNumber];
    gCurrentQyVsVzZNA = gQyVsVzZNA[currentRunNumber];
    gCurrentQxVsCentZNC = gQxVsCentZNC[currentRunNumber];
    gCurrentQyVsCentZNC = gQyVsCentZNC[currentRunNumber];
    gCurrentQxVsVxZNC = gQxVsVxZNC[currentRunNumber];
    gCurrentQyVsVxZNC = gQyVsVxZNC[currentRunNumber];
    gCurrentQxVsVyZNC = gQxVsVyZNC[currentRunNumber];
    gCurrentQyVsVyZNC = gQyVsVyZNC[currentRunNumber];
    gCurrentQxVsVzZNC = gQxVsVzZNC[currentRunNumber];
    gCurrentQyVsVzZNC = gQyVsVzZNC[currentRunNumber];
    gCurrentQxQyVsCent = gQxQyVsCent[currentRunNumber];
    gCurrentQyQxVsCent = gQyQxVsCent[currentRunNumber];
    gCurrentQxQxVsCent = gQxQxVsCent[currentRunNumber];
    gCurrentQyQyVsCent = gQyQyVsCent[currentRunNumber];

    gCurrentQxZNA = gQx5DZNA[currentRunNumber];
    gCurrentQyZNA = gQy5DZNA[currentRunNumber];
    gCurrentQxZNC = gQx5DZNC[currentRunNumber];
    gCurrentQyZNC = gQy5DZNC[currentRunNumber];

    gCurrentPsiZNA = gPsiZNA[currentRunNumber];
    gCurrentPsiZNC = gPsiZNC[currentRunNumber];

    gCurrentVx = gVx[currentRunNumber];
    gCurrentVy = gVy[currentRunNumber];

    gCurrentQxVsTimeZNA = gQxVsTimeZNA[currentRunNumber];
    gCurrentQyVsTimeZNA = gQyVsTimeZNA[currentRunNumber];
    gCurrentQxVsTimeZNC = gQxVsTimeZNC[currentRunNumber];
    gCurrentQyVsTimeZNC = gQyVsTimeZNC[currentRunNumber];

    gCurrentShiftProfileZNA = gShiftProfileZNA[currentRunNumber];
    gCurrentShiftProfileZNC = gShiftProfileZNC[currentRunNumber];
  }

  // Optimized method to load ALL calibrations for the new run at once
  void loadCalibrations(int runNumber)
  {
    clearCache();

    // Vertex Calibration
    if (applyBeamSpotCorrection) {
      std::string folder = Form("%s/step0", qRecenteringCcdb.value.c_str());
      auto* lst = ccdb->getForRun<TList>(folder, runNumber);
      if (lst) {
        hMeanVx = safeClone<TH1>(lst->FindObject("hMeanVx"));
        hMeanVy = safeClone<TH1>(lst->FindObject("hMeanVy"));
      } else {
        LOGF(error, "  >> CCDB TList is NULL for path: %s. Check object type (TList vs TFile).", folder.c_str());
      }
    }

    // Step Calibrations
    std::size_t targetSteps = (calibrationStep > 0) ? static_cast<std::size_t>(calibrationStep.value) : 0;
    calibCache.resize(targetSteps);

    for (std::size_t stepIdx = 0; stepIdx < targetSteps; ++stepIdx) {
      int step = static_cast<int>(stepIdx + 1);

      // Load 5D (Base)
      std::string folderBase = Form("%s/step%d_base", qRecenteringCcdb.value.c_str(), step);
      auto* lstBase = ccdb->getForRun<TList>(folderBase, runNumber);

      if (!lstBase) {
        LOGF(error, "  >> CCDB TList is NULL for path: %s. Check object type (TList vs TFile).", folderBase.c_str());
        continue;
      }

      calibCache[stepIdx].hMeanQxZNA = safeClone<THn>(lstBase->FindObject("hMeanQxZNA"));
      calibCache[stepIdx].hMeanQyZNA = safeClone<THn>(lstBase->FindObject("hMeanQyZNA"));
      calibCache[stepIdx].hMeanQxZNC = safeClone<THn>(lstBase->FindObject("hMeanQxZNC"));
      calibCache[stepIdx].hMeanQyZNC = safeClone<THn>(lstBase->FindObject("hMeanQyZNC"));

      // Load 1D (Refine)
      if ((step != calibrationStep) || isFineCalibrationStep) {
        std::string folderRefine = Form("%s/step%d_refine", qRecenteringCcdb.value.c_str(), step);
        auto* lstRefine = ccdb->getForRun<TList>(folderRefine, runNumber);

        if (!lstRefine) {
          LOGF(error, "  >> CCDB TList is NULL for path: %s. Check object type (TList vs TFile).", folderRefine.c_str());
          continue;
        }

        calibCache[stepIdx].hMeanQxCentZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxCentZNA"));
        calibCache[stepIdx].hMeanQyCentZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyCentZNA"));
        calibCache[stepIdx].hMeanQxCentZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxCentZNC"));
        calibCache[stepIdx].hMeanQyCentZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyCentZNC"));

        calibCache[stepIdx].hMeanQxVzZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVzZNA"));
        calibCache[stepIdx].hMeanQyVzZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVzZNA"));
        calibCache[stepIdx].hMeanQxVzZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVzZNC"));
        calibCache[stepIdx].hMeanQyVzZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVzZNC"));

        calibCache[stepIdx].hMeanQxVxZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVxZNA"));
        calibCache[stepIdx].hMeanQyVxZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVxZNA"));
        calibCache[stepIdx].hMeanQxVxZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVxZNC"));
        calibCache[stepIdx].hMeanQyVxZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVxZNC"));

        calibCache[stepIdx].hMeanQxVyZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVyZNA"));
        calibCache[stepIdx].hMeanQyVyZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVyZNA"));
        calibCache[stepIdx].hMeanQxVyZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVyZNC"));
        calibCache[stepIdx].hMeanQyVyZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVyZNC"));
      }
    } // end of step loop

    if (applyShiftCorrection) {
      std::string folder = Form("%s/psiHarm", qRecenteringCcdb.value.c_str());

      //LOGF(info, "Loading Shift Correction from %s for runNumber %d", folder.c_str(), runNumber);

      // Attempt to fetch TList from CCDB
      auto* lst = ccdb->getForRun<TList>(folder, runNumber);

      if (!lst) {
        LOGF(error, "  >> CCDB TList is NULL for path: %s. Check object type (TList vs TFile).", folder.c_str());
        return;
      }

      // Important: Object names must match exactly what was saved
      shiftProfileZNA = safeClone<TProfile3D>(lst->FindObject("ShiftProfileZNA"));
      shiftProfileZNC = safeClone<TProfile3D>(lst->FindObject("ShiftProfileZNC"));

      if (shiftProfileZNA) {
        shiftProfileZNA->SetDirectory(nullptr); // Detach from file
        //  LOGF(info, "  >> ShiftProfileZNA found! Entries: %.0f, Mean: %f", shiftProfileZNA->GetEntries(), shiftProfileZNA->GetMean());
      } else {
        LOGF(error, "  >> ShiftProfileZNA NOT found in TList! Content follows:");
        lst->Print();
      }

      if (shiftProfileZNC) {
        shiftProfileZNC->SetDirectory(nullptr);
        // LOGF(info, "  >> ShiftProfileZNC found! Entries: %.0f", shiftProfileZNC->GetEntries());
      } else {
        LOGF(error, "  >> ShiftProfileZNC NOT found in TList!");
      }
    }
  } // end of loadCalibrations()

  ~ZdcExtraTableReader()
  {
    clearCache();
  }

  /// Initializes histograms and other resources before event processing.
  void init(InitContext&)
  {
    const AxisSpec axisCounter{1, 0, +1, ""};
    histos.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    // CCDB
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setFatalWhenNull(false);
  } // end of init()

  void process(aod::ZdcExtras::iterator const& zdc)
  {

    // Apply event selection

    if (applySel8 && !TESTBIT(zdc.selectionBits(), Sel8)) {
      return;
    }
    if (applyZVtxCut && !TESTBIT(zdc.selectionBits(), ZVtxCut)) {
      return;
    }
    if (applyOccupancyCut && !TESTBIT(zdc.selectionBits(), OccupancyCut)) {
      return;
    }
    if (selectNoSameBunchPileupEvents && !TESTBIT(zdc.selectionBits(), NoSameBunchPileup)) {
      return;
    }
    if (selectGoodZvtxFT0vsPV && !TESTBIT(zdc.selectionBits(), IsGoodZvtxFT0vsPV)) {
      return;
    }
    if (applyNoCollInTimeRangeStandard && !TESTBIT(zdc.selectionBits(), NoCollInTimeRangeStandard)) {
      return;
    }
    if (selectVertexITSTPC && !TESTBIT(zdc.selectionBits(), IsVertexITSTPC)) {
      return;
    }
    if (selectGoodITSLayersAll && !TESTBIT(zdc.selectionBits(), IsGoodITSLayersAll)) {
      return;
    }

    const int runNumber = zdc.runNumber();

    const uint64_t timestamp = zdc.timestamp(); // in milliseconds

    // Convert timestamp to hours from run start (approximate)
    // Store first timestamp of run to calculate relative time
    static std::unordered_map<int, uint64_t> runStartTime;
    if (!runStartTime.contains(runNumber)) {
      runStartTime[runNumber] = timestamp;
    }

    double timeInMinutes = (timestamp - runStartTime[runNumber]) / 60000.0; // ms -> minutes

    // Initialization block if Run Number changes
    if (runNumber != currentRunNumber) {
      initHistos(runNumber);       // Init output histograms
      loadCalibrations(runNumber); // Load all steps from CCDB once
      currentRunNumber = runNumber;
    }

    histos.fill(HIST("eventCounter"), 0.5);

    gCurrentEventCounter->Fill(AllEvents);

    // Fill histogram for all bits that passed
    for (int bit = 0; bit < NEventSelections + 1; bit++) {
      if (TESTBIT(zdc.selectionBits(), bit)) {
        gCurrentEventCounter->Fill(bit);
      }
    }

    bool isZNChit = false, isZNAhit = false;
    
    double tdcZNC = zdc.zncTdc();
    double tdcZNA = zdc.znaTdc();

    if (applyTdcCut) { // TDC window is set
      if ((tdcZNC >= tdcZnMin) && (tdcZNC <= tdcZnMax)) {
        isZNChit = true;
      }
      if ((tdcZNA >= tdcZnMin) && (tdcZNA <= tdcZnMax)) {
        isZNAhit = true;
      }
    } else { // if no window on TDC is set
      if (zdc.zncTowC() > -1.) {
        isZNChit = true;
      }
      if (zdc.znaTowC() > -1.) {
        isZNAhit = true;
      }
    }
    
    bool isZNASpDeterminable = false;
    bool isZNCSpDeterminable = false;

    constexpr float QvectorMaxValue = 990.0;

    if (isZNAhit) {
      int activeTowersZNA = (zdc.znaTow1() > 0.) + (zdc.znaTow2() > 0.) + (zdc.znaTow3() > 0.) + (zdc.znaTow4() > 0.);
      float znaSum = zdc.znaTow1() + zdc.znaTow2() + zdc.znaTow3() + zdc.znaTow4();
      if (activeTowersZNA >= minNTowersFired && znaSum > 0 && zdc.znaQx() < QvectorMaxValue)
        isZNASpDeterminable = true;
    }

    if (isZNChit) {
      int activeTowersZNC = (zdc.zncTow1() > 0.) + (zdc.zncTow2() > 0.) + (zdc.zncTow3() > 0.) + (zdc.zncTow4() > 0.);
      float zncSum = zdc.zncTow1() + zdc.zncTow2() + zdc.zncTow3() + zdc.zncTow4();
      if (activeTowersZNC >= minNTowersFired && zncSum > 0 && zdc.zncQx() < QvectorMaxValue)
        isZNCSpDeterminable = true;
    }

    if (plotPMs) {
      if (isZNAhit) {
        gCurrentPmcZNA->Fill(zdc.znaTowC());
        gCurrentPm1ZNA->Fill(zdc.znaTow1());
        gCurrentPm2ZNA->Fill(zdc.znaTow2());
        gCurrentPm3ZNA->Fill(zdc.znaTow3());
        gCurrentPm4ZNA->Fill(zdc.znaTow4());

        float znaSum = zdc.znaTow1() + zdc.znaTow2() + zdc.znaTow3() + zdc.znaTow4();
        gCurrentSumZNA->Fill(znaSum);
      }

      if (isZNChit) {
        gCurrentPmcZNC->Fill(zdc.zncTowC());
        gCurrentPm1ZNC->Fill(zdc.zncTow1());
        gCurrentPm2ZNC->Fill(zdc.zncTow2());
        gCurrentPm3ZNC->Fill(zdc.zncTow3());
        gCurrentPm4ZNC->Fill(zdc.zncTow4());

        float zncSum = zdc.zncTow1() + zdc.zncTow2() + zdc.zncTow3() + zdc.zncTow4();
        gCurrentSumZNC->Fill(zncSum);
      }
    } // end of if (plotPMs)

    double qxZNArec = 0., qyZNArec = 0.;
    double qxZNCrec = 0., qyZNCrec = 0.;
    
    double cent = zdc.centrality();
    double vx = zdc.vx();
    double vy = zdc.vy();
    double vz = zdc.vz();


    if (applyBeamSpotCorrection) {
      // Use cached vertex pointers
      if (hMeanVx && hMeanVy) {
        vx -= hMeanVx->GetBinContent(1);
        vy -= hMeanVy->GetBinContent(1);
      }
    }

    // -------- ZNA --------
    if (isZNASpDeterminable) {
      double qx = zdc.znaQx();
      double qy = zdc.znaQy();

      qxZNArec = qx;
      qyZNArec = qy;
      
      for (int step = 1; step <= calibrationStep; step++) {

        int cacheIdx = step - 1;
        // Check if index is valid within cached vector
        if (cacheIdx >= static_cast<int>(calibCache.size()))
          continue;

        const auto& calib = calibCache[cacheIdx];

        // Apply 5D Base calibration
        if (calib.hMeanQxZNA && calib.hMeanQyZNA) {
          qxZNArec -= getMeanQFromMap(calib.hMeanQxZNA, cent, vx, vy, vz);
          qyZNArec -= getMeanQFromMap(calib.hMeanQyZNA, cent, vx, vy, vz);
        }

        if ((step != calibrationStep) || isFineCalibrationStep) {
          // Apply 1D Refine calibration
          qxZNArec -= getMeanQ1D(calib.hMeanQxCentZNA, cent);
          qyZNArec -= getMeanQ1D(calib.hMeanQyCentZNA, cent);

          qxZNArec -= getMeanQ1D(calib.hMeanQxVzZNA, vz);
          qyZNArec -= getMeanQ1D(calib.hMeanQyVzZNA, vz);

          qxZNArec -= getMeanQ1D(calib.hMeanQxVxZNA, vx);
          qyZNArec -= getMeanQ1D(calib.hMeanQyVxZNA, vx);

          qxZNArec -= getMeanQ1D(calib.hMeanQxVyZNA, vy);
          qyZNArec -= getMeanQ1D(calib.hMeanQyVyZNA, vy);
        }
      }

      double valuesQxZNA[5] = {cent, vx, vy, vz, qxZNArec};
      double valuesQyZNA[5] = {cent, vx, vy, vz, qyZNArec};

      gCurrentCentroidZNA->Fill(qxZNArec, qyZNArec);

      gCurrentQxVsCentZNA->Fill(cent, qxZNArec);
      gCurrentQyVsCentZNA->Fill(cent, qyZNArec);
      gCurrentQxVsVxZNA->Fill(vx, qxZNArec);
      gCurrentQyVsVxZNA->Fill(vx, qyZNArec);
      gCurrentQxVsVyZNA->Fill(vy, qxZNArec);
      gCurrentQyVsVyZNA->Fill(vy, qyZNArec);
      gCurrentQxVsVzZNA->Fill(vz, qxZNArec);
      gCurrentQyVsVzZNA->Fill(vz, qyZNArec);

      // Fill time-dependent plots
      gCurrentQxVsTimeZNA->Fill(timeInMinutes, qxZNArec);
      gCurrentQyVsTimeZNA->Fill(timeInMinutes, qyZNArec);

      if (plot5D) {
        gCurrentQxZNA->Fill(valuesQxZNA);
        gCurrentQyZNA->Fill(valuesQyZNA);
      }

      //  Calculate raw/recentered angle
      double psiZNA = std::atan2(qyZNArec, qxZNArec);

      // Apply Correction (Read Mode)
      // Checks if correction is enabled AND if the map from CCDB was loaded successfully
      if (applyShiftCorrection && shiftProfileZNA) {
        double deltaPsi = 0.0;

        // Loop over harmonics (usually 1 to 10)
        for (int iHarm = 1; iHarm <= nHarmonics; iHarm++) {
          // Retrieve coefficients from TProfile3D
          // Axis mapping:
          // X: Centrality
          // Y: Type (0.5 for Sin, 1.5 for Cos)
          // Z: Harmonic index (iHarm - 0.5 maps to bin 1, 2, etc.)

          int binSin = shiftProfileZNA->FindFixBin(cent, 0.5, static_cast<double>(iHarm) - 0.5);
          int binCos = shiftProfileZNA->FindFixBin(cent, 1.5, static_cast<double>(iHarm) - 0.5);

          double coeffSin = shiftProfileZNA->GetBinContent(binSin);
          double coeffCos = shiftProfileZNA->GetBinContent(binCos);

          // Fourier flattening formula:
          // DeltaPsi = sum( (2/k) * ( <cos>*sin(k*psi) - <sin>*cos(k*psi) ) )
          // Note: signs depend on definition, this matches the standard correction logic
          deltaPsi += (2.0 / iHarm) * (-coeffSin * std::cos(iHarm * psiZNA) + coeffCos * std::sin(iHarm * psiZNA));
        }

        // DEBUG: Print only if shift is actually happening for first few events
        static int debugPrintCount = 0;
        constexpr int MaxDebugPrints = 10;
        constexpr double PsiTolerance = 1e-6;

        if (debugPrintCount < MaxDebugPrints && std::abs(deltaPsi) > PsiTolerance) {
          LOGF(info, "ZNA Shift: Cent %.1f, Raw %.3f (Delta %.4f)", cent, psiZNA, deltaPsi);
          debugPrintCount++;
        }

        // Apply the calculated shift
        psiZNA += deltaPsi;

        // Wrap angle to [-pi, pi] range
        psiZNA = RecoDecay::constrainAngle(psiZNA, -o2::constants::math::PI);
      }

      // Fill Shift Profiles (Write Mode)
      // Used to generate calibration for the next step or to verify correction (QA)
      if (fillShiftHistos && gCurrentShiftProfileZNA) {
        for (int iHarm = 1; iHarm <= nHarmonics; iHarm++) {
          // Fill Sin component (Y = 0.5)
          gCurrentShiftProfileZNA->Fill(cent, 0.5, static_cast<double>(iHarm) - 0.5, std::sin(iHarm * psiZNA));
          // Fill Cos component (Y = 1.5)
          gCurrentShiftProfileZNA->Fill(cent, 1.5, static_cast<double>(iHarm) - 0.5, std::cos(iHarm * psiZNA));
        }
      }

      // Fill final analysis histogram with the best available Psi (Raw or Corrected)
      gCurrentPsiZNA->Fill(psiZNA);
    }

    // -------- ZNC --------
    if (isZNCSpDeterminable) {
      double qx = zdc.zncQx();
      double qy = zdc.zncQy();

      qxZNCrec = qx;
      qyZNCrec = qy;
      
      // Iterate through steps using cached vector
      for (int step = 1; step <= calibrationStep; step++) {

        int cacheIdx = step - 1;
        if (cacheIdx >= static_cast<int>(calibCache.size()))
          continue;
        
        const auto& calib = calibCache[cacheIdx];

        // Apply 5D Base calibration
        if (calib.hMeanQxZNC && calib.hMeanQyZNC) {
          qxZNCrec -= getMeanQFromMap(calib.hMeanQxZNC, cent, vx, vy, vz);
          qyZNCrec -= getMeanQFromMap(calib.hMeanQyZNC, cent, vx, vy, vz);
        }

        if ((step != calibrationStep) || isFineCalibrationStep) {

          // Apply 1D Refine calibration
          qxZNCrec -= getMeanQ1D(calib.hMeanQxCentZNC, cent);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyCentZNC, cent);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVzZNC, vz);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVzZNC, vz);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVxZNC, vx);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVxZNC, vx);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVyZNC, vy);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVyZNC, vy);
        }
      }

      double valuesQxZNC[5] = {cent, vx, vy, vz, qxZNCrec};
      double valuesQyZNC[5] = {cent, vx, vy, vz, qyZNCrec};
      
      gCurrentCentroidZNC->Fill(qxZNCrec, qyZNCrec);

      gCurrentQxVsCentZNC->Fill(cent, qxZNCrec);
      gCurrentQyVsCentZNC->Fill(cent, qyZNCrec);
      gCurrentQxVsVxZNC->Fill(vx, qxZNCrec);
      gCurrentQyVsVxZNC->Fill(vx, qyZNCrec);
      gCurrentQxVsVyZNC->Fill(vy, qxZNCrec);
      gCurrentQyVsVyZNC->Fill(vy, qyZNCrec);
      gCurrentQxVsVzZNC->Fill(vz, qxZNCrec);
      gCurrentQyVsVzZNC->Fill(vz, qyZNCrec);

      // Fill time-dependent plots
      gCurrentQxVsTimeZNC->Fill(timeInMinutes, qxZNCrec);
      gCurrentQyVsTimeZNC->Fill(timeInMinutes, qyZNCrec);

      if (plot5D) {
        gCurrentQxZNC->Fill(valuesQxZNC);
        gCurrentQyZNC->Fill(valuesQyZNC);
      }

      //  Calculate raw/recentered angle
      double psiZNC = std::atan2(qyZNCrec, qxZNCrec);

      // Apply Correction (Read Mode)
      // Checks if correction is enabled AND if the map from CCDB was loaded successfully
      if (applyShiftCorrection && shiftProfileZNC) {
        double deltaPsi = 0.0;

        // Loop over harmonics (usually 1 to 10)
        for (int iHarm = 1; iHarm <= nHarmonics; iHarm++) {
          // Retrieve coefficients from TProfile3D
          // Axis mapping:
          // X: Centrality
          // Y: Type (0.5 for Sin, 1.5 for Cos)
          // Z: Harmonic index (iHarm - 0.5 maps to bin 1, 2, etc.)

          int binSin = shiftProfileZNC->FindFixBin(cent, 0.5, static_cast<double>(iHarm) - 0.5);
          int binCos = shiftProfileZNC->FindFixBin(cent, 1.5, static_cast<double>(iHarm) - 0.5);

          double coeffSin = shiftProfileZNC->GetBinContent(binSin);
          double coeffCos = shiftProfileZNC->GetBinContent(binCos);

          // Fourier flattening formula:
          // DeltaPsi = sum( (2/k) * ( <cos>*sin(k*psi) - <sin>*cos(k*psi) ) )
          // Note: signs depend on definition, this matches the standard correction logic
          deltaPsi += (2.0 / iHarm) * (-coeffSin * std::cos(iHarm * psiZNC) + coeffCos * std::sin(iHarm * psiZNC));
        }

        // Apply the calculated shift
        psiZNC += deltaPsi;

        // Wrap angle to [-pi, pi] range
        psiZNC = RecoDecay::constrainAngle(psiZNC, -o2::constants::math::PI);
      }

      // Fill Shift Profiles (Write Mode)
      // Used to generate calibration for the next step or to verify correction (QA)
      if (fillShiftHistos && gCurrentShiftProfileZNC) {
        for (int iHarm = 1; iHarm <= nHarmonics; iHarm++) {
          // Fill Sin component (Y = 0.5)
          gCurrentShiftProfileZNC->Fill(cent, 0.5, static_cast<double>(iHarm) - 0.5, std::sin(iHarm * psiZNC));
          // Fill Cos component (Y = 1.5)
          gCurrentShiftProfileZNC->Fill(cent, 1.5, static_cast<double>(iHarm) - 0.5, std::cos(iHarm * psiZNC));
        }
      }

      // Fill final analysis histogram with the best available Psi (Raw or Corrected)
      gCurrentPsiZNC->Fill(psiZNC);
    }

    if (isZNASpDeterminable && isZNCSpDeterminable) {
      gCurrentQxQyVsCent->Fill(cent, qxZNArec * qyZNCrec);
      gCurrentQyQxVsCent->Fill(cent, qyZNArec * qxZNCrec);
      gCurrentQxQxVsCent->Fill(cent, qxZNArec * qxZNCrec);
      gCurrentQyQyVsCent->Fill(cent, qyZNArec * qyZNCrec);
    }

    gCurrentVx->Fill(vx);
    gCurrentVy->Fill(vy);

  } // end of process() (= end of loop through collisions with ZDC)
}; // end of struct ZdcExtraTableReader

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ZdcExtraTableReader>(cfgc)};
}
