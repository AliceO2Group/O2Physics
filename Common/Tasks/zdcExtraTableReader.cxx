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

#include <CommonConstants/MathConstants.h>      

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <THn.h>
#include <TList.h>
#include <TMath.h>
#include <TProfile3D.h>

#include <memory>
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

TH1* gCurrentEventCounter;
TH2* gCurrentCentroidZNA;
TH2* gCurrentCentroidZNC;
TH1* gCurrentPmcZNA;
TH1* gCurrentPm1ZNA;
TH1* gCurrentPm2ZNA;
TH1* gCurrentPm3ZNA;
TH1* gCurrentPm4ZNA;
TH1* gCurrentSumZNA;
TH1* gCurrentPmcZNC;
TH1* gCurrentPm1ZNC;
TH1* gCurrentPm2ZNC;
TH1* gCurrentPm3ZNC;
TH1* gCurrentPm4ZNC;
TH1* gCurrentSumZNC;
TH2* gCurrentQxVsCentZNA;
TH2* gCurrentQyVsCentZNA;
TH2* gCurrentQxVsVxZNA;
TH2* gCurrentQyVsVxZNA;
TH2* gCurrentQxVsVyZNA;
TH2* gCurrentQyVsVyZNA;
TH2* gCurrentQxVsVzZNA;
TH2* gCurrentQyVsVzZNA;
TH2* gCurrentQxVsCentZNC;
TH2* gCurrentQyVsCentZNC;
TH2* gCurrentQxVsVxZNC;
TH2* gCurrentQyVsVxZNC;
TH2* gCurrentQxVsVyZNC;
TH2* gCurrentQyVsVyZNC;
TH2* gCurrentQxVsVzZNC;
TH2* gCurrentQyVsVzZNC;
TH2* gCurrentQxQyVsCent;
TH2* gCurrentQyQxVsCent;
TH2* gCurrentQxQxVsCent;
TH2* gCurrentQyQyVsCent;

THn* gCurrentQxZNA;
THn* gCurrentQyZNA;
THn* gCurrentQxZNC;
THn* gCurrentQyZNC;

TH1* gCurrentPsiZNA;
TH1* gCurrentPsiZNC;

TH1* gCurrentVx;
TH1* gCurrentVy;

TH2* gCurrentQxVsTimeZNA;
TH2* gCurrentQyVsTimeZNA;
TH2* gCurrentQxVsTimeZNC;
TH2* gCurrentQyVsTimeZNC;

TProfile3D* gCurrentShiftProfileZNA;
TProfile3D* gCurrentShiftProfileZNC;

} // namespace

double getMeanQFromMap(THn* h, double cent, double vx, double vy, double vz)
{
  if (!h) {
    // Commented out to reduce log spam in case of missing maps, enable for debugging
    // std::cerr << "[MeanQ] Null THn pointer (cent=" << cent << ", vx=" << vx << ", vy=" << vy << ", vz=" << vz << ")" << std::endl;
    return 0.0;
  }

  TAxis* axCent = h->GetAxis(0);
  TAxis* axVx = h->GetAxis(1);
  TAxis* axVy = h->GetAxis(2);
  TAxis* axVz = h->GetAxis(3);

  if (!axCent || !axVx || !axVy || !axVz) {
    std::cerr << "[MeanQ] One of THn axes is null" << std::endl;
    return 0.0;
  }

  int binCent = axCent->FindBin(cent);
  int binVx = axVx->FindBin(vx);
  int binVy = axVy->FindBin(vy);
  int binVz = axVz->FindBin(vz);

  int idx[4] = {binCent, binVx, binVy, binVz};
  double meanQ = h->GetBinContent(idx);

  return meanQ;
}

// Helper for 1D recentering maps: returns mean Q for coordinate x.
// If histogram is missing or bin out of range, returns 0.0.
double getMeanQ1D(TH1* h, double x)
{
  if (!h) {
    return 0.0;
  }
  int bin = h->FindBin(x);
  if (bin < 1 || bin > h->GetNbinsX()) {
    return 0.0;
  }
  return h->GetBinContent(bin);
}

struct ZdcExtraTableReader {

  Configurable<int> nBinsZN{"nBinsZN", 2000, "n bins for ZN histograms"};
  Configurable<float> maxZN{"maxZN", 399.5, "Max ZN signal"};
  Configurable<bool> tdcCut{"tdcCut", true, "Flag for TDC cut"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -2.5, "Min ZN TDC cut"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", 2.5, "Max ZN TDC cut"};
  Configurable<bool> plotPMs{"plotPMs", false, "Flag to plot individual PMs"};

  Configurable<int> qxyNbins{"qxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<float> qxyMin{"qxyMin", -2.0f, "Lower edge for QxQy histograms"};
  Configurable<float> qxyMax{"qxyMax", 2.0f, "Upper edge for QxQy histograms"};

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

  Configurable<int> qNbins5D{"qNbins5D", 4, "Bins in each dimension for 5D histograms"};
  Configurable<bool> plot5D{"plot5D", false, "Flag to plot 5D histograms"};

  Configurable<int> calibrationStep{"calibrationStep", 1, "Calibration step"};
  Configurable<bool> ifFineCalibration{"ifFineCalibration", false, "Calibration: base  or refine"};

  Configurable<bool> ifBeamSpotCorrection{"ifBeamSpotCorrection", true, "Beam spot correction"};

  Configurable<bool> ifSel8{"ifSel8", true, "Event selection: sel8"};
  Configurable<bool> ifVtxZle10{"ifVtxZle10", true, "Event selection: zVtx < 10 cm"};
  Configurable<bool> ifOccupancyCut{"ifOccupancyCut", false, "Event selection: occupancy cut"};
  Configurable<bool> ifNoSameBunchPileup{"ifNoSameBunchPileup", false, "Event selection: no same bunch pileup"};
  Configurable<bool> ifIsGoodZvtxFT0vsPV{"ifIsGoodZvtxFT0vsPV", false, "Event selection: good Zvtx FT0 vs PV"};
  Configurable<bool> ifNoCollInTimeRangeStandard{"ifNoCollInTimeRangeStandard", false, "Event selection: no collision in time range standard"};
  Configurable<bool> ifIsVertexITSTPC{"ifIsVertexITSTPC", false, "Event selection: vertex ITS TPC"};
  Configurable<bool> ifIsGoodITSLayersAll{"ifIsGoodITSLayersAll", false, "Event selection: good ITS layers all"};

  Configurable<bool> ifShiftCorrection{"ifShiftCorrection", false, "Apply shift correction (Read from CCDB)"};
  Configurable<bool> fillShiftHistos{"fillShiftHistos", true, "Fill shift profiles (Write to output)"};
  Configurable<int> nShift{"nShift", 10, "Number of harmonics"};

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
  };

  // Cache container: Vector index = Step index (0-based, so step 1 is at index 0)
  std::vector<CalibStepData> mCalibCache;

  // Vertex correction cache
  TH1* hMeanVx{nullptr};
  TH1* hMeanVy{nullptr};

  // Phase shift correction cache
  TProfile3D* hShiftZNA{nullptr};
  TProfile3D* hShiftZNC{nullptr};

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  enum EvSelBits { // TO DO: move to a common header file
    evSel_zvtx,
    evSel_sel8,
    evSel_occupancy,
    evSel_kNoSameBunchPileup,
    evSel_kIsGoodZvtxFT0vsPV,
    evSel_kNoCollInTimeRangeStandard,
    evSel_kIsVertexITSTPC,
    evSel_kIsGoodITSLayersAll,
    evSel_allEvents,
    nEventSelections
  };

  int mCurrentRunNumber{-1};

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
    if (hMeanVx) {
      delete hMeanVx;
      hMeanVx = nullptr;
    }
    if (hMeanVy) {
      delete hMeanVy;
      hMeanVy = nullptr;
    }

    for (const auto& step : mCalibCache) {
      delete step.hMeanQxZNA;
      delete step.hMeanQyZNA;
      delete step.hMeanQxZNC;
      delete step.hMeanQyZNC;

      delete step.hMeanQxCentZNA;
      delete step.hMeanQyCentZNA;
      delete step.hMeanQxCentZNC;
      delete step.hMeanQyCentZNC;

      delete step.hMeanQxVzZNA;
      delete step.hMeanQyVzZNA;
      delete step.hMeanQxVzZNC;
      delete step.hMeanQyVzZNC;

      delete step.hMeanQxVxZNA;
      delete step.hMeanQyVxZNA;
      delete step.hMeanQxVxZNC;
      delete step.hMeanQyVxZNC;

      delete step.hMeanQxVyZNA;
      delete step.hMeanQyVyZNA;
      delete step.hMeanQxVyZNC;
      delete step.hMeanQyVyZNC;

      if (hShiftZNA) {
        delete hShiftZNA;
        hShiftZNA = nullptr;
      }
      if (hShiftZNC) {
        delete hShiftZNC;
        hShiftZNC = nullptr;
      }
    }
    mCalibCache.clear();
  }

  void initHistos(const int& mRunNumber)
  {
    if (mRunNumber == mCurrentRunNumber) {
      return;
    }
    mCurrentRunNumber = mRunNumber;

    if (gEventCounter.find(mRunNumber) == gEventCounter.end()) {
      // if new run, initialize histograms

      const AxisSpec axisCounter{1, 0, +1, ""};
      const AxisSpec axisZN{nBinsZN, -0.5, maxZN, "(a.u.)"};
      const AxisSpec axisCent = {centNbins, centMin, centMax, "Centrality (\%)"};
      const AxisSpec axisVx = {vxNbins, vxMin, vxMax, "V_{x} (cm)"};
      const AxisSpec axisVy = {vyNbins, vyMin, vyMax, "V_{y} (cm)"};
      const AxisSpec axisVz = {vzNbins, vzMin, vzMax, "V_{z} (cm)"};

      const AxisSpec axisCent5D = {qNbins5D, centMin, centMax, "Centrality (\%)"};
      const AxisSpec axisVx5D = {qNbins5D, vxMin, vxMax, "V_{x} (cm)"};
      const AxisSpec axisVy5D = {qNbins5D, vyMin, vyMax, "V_{y} (cm)"};
      const AxisSpec axisVz5D = {qNbins5D, vzMin, vzMax, "V_{z} (cm)"};

      const AxisSpec axisQx = {qxyNbins, qxyMin, qxyMax, "Q_{x}"};
      const AxisSpec axisQy = {qxyNbins, qxyMin, qxyMax, "Q_{y}"};

      const AxisSpec axisQxQy = {qxyNbins, qxyMin, qxyMax, ""};
      const AxisSpec axisPhi = {phiNbins, -1.0f * o2::constants::math::PI, 1.0f * o2::constants::math::PI, "#phi"};

      const AxisSpec axisTime = {90, 0, 90, "Time (minutes)"}; // 90 minutes

      gEventCounter[mRunNumber] = histos.add<TH1>(Form("%i/eventCounter", mRunNumber), "Number of Event; ; #Events Passed Cut", kTH1D, {{nEventSelections, 0, nEventSelections}}).get();
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_allEvents + 1, "All events");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_zvtx + 1, "vtxZ");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_sel8 + 1, "Sel8");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_occupancy + 1, "kOccupancy");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_kNoSameBunchPileup + 1, "kNoSameBunchPileup");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_kIsGoodZvtxFT0vsPV + 1, "kIsGoodZvtxFT0vsPV");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_kNoCollInTimeRangeStandard + 1, "kNoCollInTimeRangeStandard");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_kIsVertexITSTPC + 1, "kIsVertexITSTPC");
      gEventCounter[mRunNumber]->GetXaxis()->SetBinLabel(evSel_kIsGoodITSLayersAll + 1, "kIsGoodITSLayersAll");

      gCentroidZNA[mRunNumber] = histos.add<TH2>(Form("%i/CentroidZNA", mRunNumber), "ZNA Centroid; Q_{X}; Q_{Y}", kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}).get();
      gCentroidZNC[mRunNumber] = histos.add<TH2>(Form("%i/CentroidZNC", mRunNumber), "ZNC Centroid; Q_{X}; Q_{Y}", kTH2F, {{50, -1.5, 1.5}, {50, -1.5, 1.5}}).get();
      gPmcZNA[mRunNumber] = histos.add<TH1>(Form("%i/pmcZNA", mRunNumber), "; E_{PMC}^{ZNA} (TeV);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm1ZNA[mRunNumber] = histos.add<TH1>(Form("%i/pm1ZNA", mRunNumber), "; E_{PM1}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm2ZNA[mRunNumber] = histos.add<TH1>(Form("%i/pm2ZNA", mRunNumber), "; E_{PM2}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm3ZNA[mRunNumber] = histos.add<TH1>(Form("%i/pm3ZNA", mRunNumber), "; E_{PM3}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm4ZNA[mRunNumber] = histos.add<TH1>(Form("%i/pm4ZNA", mRunNumber), "; E_{PM4}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gSumZNA[mRunNumber] = histos.add<TH1>(Form("%i/sumZNA", mRunNumber), "; E_{sum PMs}^{ZNA} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPmcZNC[mRunNumber] = histos.add<TH1>(Form("%i/pmcZNC", ж, ждать его на ленте, подстраиваться под заселение в отель и доплачивать авиакомпании.Поэтому сегодня все больше людей осознанно выбирают формат luggage free travel — путешествия только с ручной кладью.mRunNumber), "; E_{PMC}^{ZNC} (TeV);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm1ZNC[mRunNumber] = histos.add<TH1>(Form("%i/pm1ZNC", mRunNumber), "; E_{PM1}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm2ZNC[mRunNumber] = histos.add<TH1>(Form("%i/pm2ZNC", mRunNumber), "; E_{PM2}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm3ZNC[mRunNumber] = histos.add<TH1>(Form("%i/pm3ZNC", mRunNumber), "; E_{PM3}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gPm4ZNC[mRunNumber] = histos.add<TH1>(Form("%i/pm4ZNC", mRunNumber), "; E_{PM4}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gSumZNC[mRunNumber] = histos.add<TH1>(Form("%i/sumZNC", mRunNumber), "; E_{sum PMs}^{ZNC} (a.u.);", kTH1F, {{nBinsZN, -0.5, maxZN}}).get();
      gQxVsCentZNA[mRunNumber] = histos.add<TH2>(Form("%i/QxVsCentZNA", mRunNumber), "Q_{x}^{ZNA} vs Centrality", kTH2F, {axisCent, axisQx}).get();
      gQyVsCentZNA[mRunNumber] = histos.add<TH2>(Form("%i/QyVsCentZNA", mRunNumber), "Q_{y}^{ZNA} vs Centrality", kTH2F, {axisCent, axisQy}).get();
      gQxVsVxZNA[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVxZNA", mRunNumber), "Q_{x}^{ZNA} vs V_{x}; V_{x} (cm); Q_{x}", kTH2F, {axisVx, axisQx}).get();
      gQyVsVxZNA[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVxZNA", mRunNumber), "Q_{y}^{ZNA} vs V_{x}; V_{x} (cm); Q_{y}", kTH2F, {axisVx, axisQy}).get();
      gQxVsVyZNA[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVyZNA", mRunNumber), "Q_{x}^{ZNA} vs V_{y}; V_{y} (cm); Q_{x}", kTH2F, {axisVy, axisQx}).get();
      gQyVsVyZNA[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVyZNA", mRunNumber), "Q_{y}^{ZNA} vs V_{y}; V_{y} (cm); Q_{y}", kTH2F, {axisVy, axisQy}).get();
      gQxVsVzZNA[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVzZNA", mRunNumber), "Q_{x}^{ZNA} vs V_{z}; V_{z} (cm); Q_{x}", kTH2F, {axisVz, axisQx}).get();
      gQyVsVzZNA[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVzZNA", mRunNumber), "Q_{y}^{ZNA} vs V_{z}; V_{z} (cm); Q_{y}", kTH2F, {axisVz, axisQy}).get();
      gQxVsCentZNC[mRunNumber] = histos.add<TH2>(Form("%i/QxVsCentZNC", mRunNumber), "Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{x}", kTH2F, {axisCent, axisQx}).get();
      gQyVsCentZNC[mRunNumber] = histos.add<TH2>(Form("%i/QyVsCentZNC", mRunNumber), "Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{y}", kTH2F, {axisCent, axisQy}).get();
      gQxVsVxZNC[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVxZNC", mRunNumber), "Q_{x}^{ZNC} vs V_{x}; V_{x} (cm); Q_{x}", kTH2F, {axisVx, axisQx}).get();
      gQyVsVxZNC[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVxZNC", mRunNumber), "Q_{y}^{ZNC} vs V_{x}; V_{x} (cm); Q_{y}", kTH2F, {axisVx, axisQy}).get();
      gQxVsVyZNC[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVyZNC", mRunNumber), "Q_{x}^{ZNC} vs V_{y}; V_{y} (cm); Q_{x}", kTH2F, {axisVy, axisQx}).get();
      gQyVsVyZNC[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVyZNC", mRunNumber), "Q_{y}^{ZNC} vs V_{y}; V_{y} (cm); Q_{y}", kTH2F, {axisVy, axisQy}).get();
      gQxVsVzZNC[mRunNumber] = histos.add<TH2>(Form("%i/QxVsVzZNC", mRunNumber), "Q_{x}^{ZNC} vs V_{z}; V_{z} (cm); Q_{x}", kTH2F, {axisVz, axisQx}).get();
      gQyVsVzZNC[mRunNumber] = histos.add<TH2>(Form("%i/QyVsVzZNC", mRunNumber), "Q_{y}^{ZNC} vs V_{z}; V_{z} (cm); Q_{y}", kTH2F, {axisVz, axisQy}).get();
      gQxQyVsCent[mRunNumber] = histos.add<TH2>(Form("%i/QxQyVsCent", mRunNumber), "Q_{x}^{ZNC}Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{x}^{ZNA}Q_{y}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQyQxVsCent[mRunNumber] = histos.add<TH2>(Form("%i/QyQxVsCent", mRunNumber), "Q_{y}^{ZNC}Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{y}^{ZNA}Q_{x}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQxQxVsCent[mRunNumber] = histos.add<TH2>(Form("%i/QxQxVsCent", mRunNumber), "Q_{x}^{ZNC}Q_{x}^{ZNC} vs Centrality; Centrality (%); Q_{x}^{ZNA}Q_{x}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();
      gQyQyVsCent[mRunNumber] = histos.add<TH2>(Form("%i/QyQyVsCent", mRunNumber), "Q_{y}^{ZNC}Q_{y}^{ZNC} vs Centrality; Centrality (%); Q_{y}^{ZNA}Q_{y}^{ZNC}", kTH2F, {axisCent, {50, -1.5, 1.5}}).get();

      gQx5DZNA[mRunNumber] = histos.add<THn>(Form("%i/Qx5DZNA", mRunNumber), "Qx recenter map ZNA", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQx}, true).get();
      gQy5DZNA[mRunNumber] = histos.add<THn>(Form("%i/Qy5DZNA", mRunNumber), "Qy recenter map ZNA", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQy}, true).get();
      gQx5DZNC[mRunNumber] = histos.add<THn>(Form("%i/Qx5DZNC", mRunNumber), "Qx recenter map ZNC", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQx}, true).get();
      gQy5DZNC[mRunNumber] = histos.add<THn>(Form("%i/Qy5DZNC", mRunNumber), "Qy recenter map ZNC", kTHnF, {axisCent5D, axisVx5D, axisVy5D, axisVz5D, axisQy}, true).get();

      gPsiZNA[mRunNumber] = histos.add<TH1>(Form("%i/PsiZNA", mRunNumber), ";#Phi_{ZNA} (rad)", kTH1F, {axisPhi}).get();
      gPsiZNC[mRunNumber] = histos.add<TH1>(Form("%i/PsiZNC", mRunNumber), ";#Phi_{ZNC} (rad)", kTH1F, {axisPhi}).get();

      gVx[mRunNumber] = histos.add<TH1>(Form("%i/Vx", mRunNumber), "V_{x} distribution; V_{x} (cm); Entries", kTH1F, {axisVx}).get();
      gVy[mRunNumber] = histos.add<TH1>(Form("%i/Vy", mRunNumber), "V_{y} distribution; V_{y} (cm); Entries", kTH1F, {axisVy}).get();

      gQxVsTimeZNA[mRunNumber] = histos.add<TH2>(Form("%i/QxVsTimeZNA", mRunNumber), "Q_{x}^{ZNA} vs Time; Time (minutes); Q_{x}", kTH2F, {axisTime, axisQx}).get();
      gQyVsTimeZNA[mRunNumber] = histos.add<TH2>(Form("%i/QyVsTimeZNA", mRunNumber), "Q_{y}^{ZNA} vs Time; Time (minutes); Q_{y}", kTH2F, {axisTime, axisQy}).get();
      gQxVsTimeZNC[mRunNumber] = histos.add<TH2>(Form("%i/QxVsTimeZNC", mRunNumber), "Q_{x}^{ZNC} vs Time; Time (minutes); Q_{x}", kTH2F, {axisTime, axisQx}).get();
      gQyVsTimeZNC[mRunNumber] = histos.add<TH2>(Form("%i/QyVsTimeZNC", mRunNumber), "Q_{y}^{ZNC} vs Time; Time (minutes); Q_{y}", kTH2F, {axisTime, axisQy}).get();

      gShiftProfileZNA[mRunNumber] = histos.add<TProfile3D>(Form("%i/ShiftProfileZNA", mRunNumber), "ZNA Shift Coeffs;Cent;Type;Harmonic", kTProfile3D, {axisCent, {2, 0, 2}, {nShift, 0, (double)nShift}}).get();
      gShiftProfileZNC[mRunNumber] = histos.add<TProfile3D>(Form("%i/ShiftProfileZNC", mRunNumber), "ZNC Shift Coeffs;Cent;Type;Harmonic", kTProfile3D, {axisCent, {2, 0, 2}, {nShift, 0, (double)nShift}}).get();
    }

    gCurrentEventCounter = gEventCounter[mCurrentRunNumber];
    gCurrentCentroidZNA = gCentroidZNA[mCurrentRunNumber];
    gCurrentCentroidZNC = gCentroidZNC[mCurrentRunNumber];
    gCurrentPmcZNA = gPmcZNA[mCurrentRunNumber];
    gCurrentPm1ZNA = gPm1ZNA[mCurrentRunNumber];
    gCurrentPm2ZNA = gPm2ZNA[mCurrentRunNumber];
    gCurrentPm3ZNA = gPm3ZNA[mCurrentRunNumber];
    gCurrentPm4ZNA = gPm4ZNA[mCurrentRunNumber];
    gCurrentSumZNA = gSumZNA[mCurrentRunNumber];
    gCurrentPmcZNC = gPmcZNC[mCurrentRunNumber];
    gCurrentPm1ZNC = gPm1ZNC[mCurrentRunNumber];
    gCurrentPm2ZNC = gPm2ZNC[mCurrentRunNumber];
    gCurrentPm3ZNC = gPm3ZNC[mCurrentRunNumber];
    gCurrentPm4ZNC = gPm4ZNC[mCurrentRunNumber];
    gCurrentSumZNC = gSumZNC[mCurrentRunNumber];
    gCurrentQxVsCentZNA = gQxVsCentZNA[mCurrentRunNumber];
    gCurrentQyVsCentZNA = gQyVsCentZNA[mCurrentRunNumber];
    gCurrentQxVsVxZNA = gQxVsVxZNA[mCurrentRunNumber];
    gCurrentQyVsVxZNA = gQyVsVxZNA[mCurrentRunNumber];
    gCurrentQxVsVyZNA = gQxVsVyZNA[mCurrentRunNumber];
    gCurrentQyVsVyZNA = gQyVsVyZNA[mCurrentRunNumber];
    gCurrentQxVsVzZNA = gQxVsVzZNA[mCurrentRunNumber];
    gCurrentQyVsVzZNA = gQyVsVzZNA[mCurrentRunNumber];
    gCurrentQxVsCentZNC = gQxVsCentZNC[mCurrentRunNumber];
    gCurrentQyVsCentZNC = gQyVsCentZNC[mCurrentRunNumber];
    gCurrentQxVsVxZNC = gQxVsVxZNC[mCurrentRunNumber];
    gCurrentQyVsVxZNC = gQyVsVxZNC[mCurrentRunNumber];
    gCurrentQxVsVyZNC = gQxVsVyZNC[mCurrentRunNumber];
    gCurrentQyVsVyZNC = gQyVsVyZNC[mCurrentRunNumber];
    gCurrentQxVsVzZNC = gQxVsVzZNC[mCurrentRunNumber];
    gCurrentQyVsVzZNC = gQyVsVzZNC[mCurrentRunNumber];
    gCurrentQxQyVsCent = gQxQyVsCent[mCurrentRunNumber];
    gCurrentQyQxVsCent = gQyQxVsCent[mCurrentRunNumber];
    gCurrentQxQxVsCent = gQxQxVsCent[mCurrentRunNumber];
    gCurrentQyQyVsCent = gQyQyVsCent[mCurrentRunNumber];

    gCurrentQxZNA = gQx5DZNA[mCurrentRunNumber];
    gCurrentQyZNA = gQy5DZNA[mCurrentRunNumber];
    gCurrentQxZNC = gQx5DZNC[mCurrentRunNumber];
    gCurrentQyZNC = gQy5DZNC[mCurrentRunNumber];

    gCurrentPsiZNA = gPsiZNA[mCurrentRunNumber];
    gCurrentPsiZNC = gPsiZNC[mCurrentRunNumber];

    gCurrentVx = gVx[mCurrentRunNumber];
    gCurrentVy = gVy[mCurrentRunNumber];

    gCurrentQxVsTimeZNA = gQxVsTimeZNA[mCurrentRunNumber];
    gCurrentQyVsTimeZNA = gQyVsTimeZNA[mCurrentRunNumber];
    gCurrentQxVsTimeZNC = gQxVsTimeZNC[mCurrentRunNumber];
    gCurrentQyVsTimeZNC = gQyVsTimeZNC[mCurrentRunNumber];

    gCurrentShiftProfileZNA = gShiftProfileZNA[mCurrentRunNumber];
    gCurrentShiftProfileZNC = gShiftProfileZNC[mCurrentRunNumber];
  }

  // Optimized method to load ALL calibrations for the new run at once
  void loadCalibrations(int run)
  {
    clearCache();

    // Vertex Calibration
    if (ifBeamSpotCorrection) {
      std::string folder = Form("%s/step0", qRecenteringCcdb.value.c_str());
      TList* lst = ccdb->getForRun<TList>(folder, run);
      if (lst) {
        hMeanVx = safeClone<TH1>(lst->FindObject("hMeanVx"));
        hMeanVy = safeClone<TH1>(lst->FindObject("hMeanVy"));
      }
    }

    // Step Calibrations
    std::size_t targetSteps = (calibrationStep > 0) ? static_cast<std::size_t>(calibrationStep.value) : 0;
    mCalibCache.resize(targetSteps);

    for (std::size_t stepIdx = 0; stepIdx < targetSteps; ++stepIdx) {
      int step = static_cast<int>(stepIdx + 1);

      // Load 5D (Base)
      std::string folderBase = Form("%s/step%d_base", qRecenteringCcdb.value.c_str(), step);
      TList* lstBase = ccdb->getForRun<TList>(folderBase, run);
      if (lstBase) {
        mCalibCache[stepIdx].hMeanQxZNA = safeClone<THn>(lstBase->FindObject("hMeanQxZNA"));
        mCalibCache[stepIdx].hMeanQyZNA = safeClone<THn>(lstBase->FindObject("hMeanQyZNA"));
        mCalibCache[stepIdx].hMeanQxZNC = safeClone<THn>(lstBase->FindObject("hMeanQxZNC"));
        mCalibCache[stepIdx].hMeanQyZNC = safeClone<THn>(lstBase->FindObject("hMeanQyZNC"));
      }

      // Load 1D (Refine)
      if ((step != calibrationStep) || ifFineCalibration) {
        std::string folderRefine = Form("%s/step%d_refine", qRecenteringCcdb.value.c_str(), step);
        TList* lstRefine = ccdb->getForRun<TList>(folderRefine, run);
        if (lstRefine) {
          mCalibCache[stepIdx].hMeanQxCentZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxCentZNA"));
          mCalibCache[stepIdx].hMeanQyCentZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyCentZNA"));
          mCalibCache[stepIdx].hMeanQxCentZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxCentZNC"));
          mCalibCache[stepIdx].hMeanQyCentZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyCentZNC"));

          mCalibCache[stepIdx].hMeanQxVzZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVzZNA"));
          mCalibCache[stepIdx].hMeanQyVzZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVzZNA"));
          mCalibCache[stepIdx].hMeanQxVzZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVzZNC"));
          mCalibCache[stepIdx].hMeanQyVzZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVzZNC"));

          mCalibCache[stepIdx].hMeanQxVxZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVxZNA"));
          mCalibCache[stepIdx].hMeanQyVxZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVxZNA"));
          mCalibCache[stepIdx].hMeanQxVxZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVxZNC"));
          mCalibCache[stepIdx].hMeanQyVxZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVxZNC"));

          mCalibCache[stepIdx].hMeanQxVyZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQxVyZNA"));
          mCalibCache[stepIdx].hMeanQyVyZNA = safeClone<TH1>(lstRefine->FindObject("hMeanQyVyZNA"));
          mCalibCache[stepIdx].hMeanQxVyZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQxVyZNC"));
          mCalibCache[stepIdx].hMeanQyVyZNC = safeClone<TH1>(lstRefine->FindObject("hMeanQyVyZNC"));
        }
      }
    } // end of step loop

    if (ifShiftCorrection) {
      std::string folder = Form("%s/psiShift", qRecenteringCcdb.value.c_str());

      LOGF(info, "ZDC Analysis: Loading Shift Correction from %s for run %d", folder.c_str(), run);

      // Attempt to fetch TList from CCDB
      auto* lst = ccdb->getForRun<TList>(folder, run);

      if (lst) {
        // Important: Object names must match exactly what was saved
        hShiftZNA = safeClone<TProfile3D>(lst->FindObject("hShiftZNA"));
        hShiftZNC = safeClone<TProfile3D>(lst->FindObject("hShiftZNC"));

        if (hShiftZNA) {
          hShiftZNA->SetDirectory(nullptr); // Detach from file
          LOGF(info, "  >> ShiftProfileZNA found! Entries: %.0f, Mean: %f",
               hShiftZNA->GetEntries(), hShiftZNA->GetMean());
        } else {
          LOGF(error, "  >> ShiftProfileZNA NOT found in TList! Content follows:");
          lst->Print();
        }

        if (hShiftZNC) {
          hShiftZNC->SetDirectory(nullptr);
          LOGF(info, "  >> ShiftProfileZNC found! Entries: %.0f", hShiftZNC->GetEntries());
        } else {
          LOGF(error, "  >> ShiftProfileZNC NOT found in TList!");
        }

      } else {
        LOGF(error, "  >> CCDB TList is NULL for path: %s. Check object type (TList vs TFile).", folder.c_str());
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

  void process(aod::ZdcExtras::iterator const& zdc /*collision*/)
  {

    // Apply event selection

    if (ifkSel8 && !(zdc.selectionBits() & (1 << evSel_sel8))) {
      return;
    }
    if (ifkVtxZle10 && !(zdc.selectionBits() & (1 << evSel_zvtx))) {
      return;
    }
    if (ifkOccupancyCut && !(zdc.selectionBits() & (1 << evSel_occupancy))) {
      return;
    }
    if (ifkNoSameBunchPileup && !(zdc.selectionBits() & (1 << evSel_kNoSameBunchPileup))) {
      return;
    }
    if (ifkIsGoodZvtxFT0vsPV && !(zdc.selectionBits() & (1 << evSel_kIsGoodZvtxFT0vsPV))) {
      return;
    }
    if (ifkNoCollInTimeRangeStandard && !(zdc.selectionBits() & (1 << evSel_kNoCollInTimeRangeStandard))) {
      return;
    }
    if (ifkIsVertexITSTPC && !(zdc.selectionBits() & (1 << evSel_kIsVertexITSTPC))) {
      return;
    }
    if (ifkIsGoodITSLayersAll && !(zdc.selectionBits() & (1 << evSel_kIsGoodITSLayersAll))) {
      return;
    }

    //  LOGF(info, "zvtx       = %d", (zdc.selectionBits() & (1 << evSel_zvtx)) > 0);
    //  LOGF(info, "sel8       = %d", (zdc.selectionBits() & (1 << evSel_sel8)) > 0);
    //  LOGF(info, "occupancy  = %d", (zdc.selectionBits() & (1 << evSel_occupancy)) > 0);
    //  LOGF(info, "noSameBC   = %d", (zdc.selectionBits() & (1 << evSel_kNoSameBunchPileup)) > 0);
    //  LOGF(info, "FT0vsPV    = %d", (zdc.selectionBits() & (1 << evSel_kIsGoodZvtxFT0vsPV)) > 0);
    //  LOGF(info, "noCollTR   = %d", (zdc.selectionBits() & (1 << evSel_kNoCollInTimeRangeStandard)) > 0);
    //  LOGF(info, "vtxITSTPC  = %d", (zdc.selectionBits() & (1 << evSel_kIsVertexITSTPC)) > 0);
    //  LOGF(info, "ITS layers = %d", (zdc.selectionBits() & (1 << evSel_kIsGoodITSLayersAll)) > 0);

    const int mRunNumber = zdc.runNumber();

    const uint64_t timestamp = zdc.timestamp(); // in milliseconds

    // Convert timestamp to hours from run start (approximate)
    // Store first timestamp of run to calculate relative time
    static std::unordered_map<int, uint64_t> runStartTime;
    if (runStartTime.find(mRunNumber) == runStartTime.end()) {
      runStartTime[mRunNumber] = timestamp;
    }

    double timeInMinutes = (timestamp - runStartTime[mRunNumber]) / 60000.0; // ms -> minutes

    // Initialization block if Run Number changes
    if (mRunNumber != mCurrentRunNumber) {
      initHistos(mRunNumber);       // Init output histograms
      loadCalibrations(mRunNumber); // Load all steps from CCDB once
      mCurrentRunNumber = mRunNumber;
    }

    histos.fill(HIST("eventCounter"), 0.5);

    gCurrentEventCounter->Fill(evSel_allEvents);

    // Fill histogram for all bits that passed
    for (int bit = 0; bit < nEventSelections + 1; bit++) {
      if (zdc.selectionBits() & (1 << bit)) {
        gCurrentEventCounter->Fill(bit);
      }
    }

    bool isZNChit = false, isZNAhit = false;
    //
    double tdcZNC = zdc.zncTdc();
    double tdcZNA = zdc.znaTdc();

    if (tdcCut) { // TDC window is set
      if ((tdcZNC >= tdcZNmincut) && (tdcZNC <= tdcZNmaxcut)) {
        isZNChit = true;
      }
      if ((tdcZNA >= tdcZNmincut) && (tdcZNA <= tdcZNmaxcut)) {
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
    //

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

    //
    double cent = zdc.centrality();
    double vx = zdc.vx();
    double vy = zdc.vy();
    double vz = zdc.vz();

    //   double vx_corrected = vx;
    // double vy_corrected = vy;

    if (ifBeamSpotCorrection) {
      // Use cached vertex pointers
      if (hMeanVx && hMeanVy) {
        vx -= hMeanVx->GetBinContent(1);
        vy -= hMeanVy->GetBinContent(1);
      }
    }

    // -------- ZNA --------
    if (isZNAhit) {
      double qx = zdc.znaQx();
      double qy = zdc.znaQy();

      qxZNArec = qx;
      qyZNArec = qy;

      //
      for (int step = 1; step <= calibrationStep; step++) {

        int cacheIdx = step - 1;
        // Check if index is valid within cached vector
        if (cacheIdx >= static_cast<int>(mCalibCache.size()))
          continue;

        const auto& calib = mCalibCache[cacheIdx];

        // Apply 5D Base calibration
        if (calib.hMeanQxZNA && calib.hMeanQyZNA) {
          qxZNArec -= getMeanQFromMap(calib.hMeanQxZNA, cent, vx, vy, vz);
          qyZNArec -= getMeanQFromMap(calib.hMeanQyZNA, cent, vx, vy, vz);
        }

        if ((step != calibrationStep) || ifFineCalibration) {
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
      double psiZNA = TMath::ATan2(qyZNArec, qxZNArec);

      // Apply Correction (Read Mode)
      // Checks if correction is enabled AND if the map from CCDB was loaded successfully
      if (ifShiftCorrection && hShiftZNA) {
        double deltaPsi = 0.0;

        // Loop over harmonics (usually 1 to 10)
        for (int ishift = 1; ishift <= nShift; ishift++) {
          // Retrieve coefficients from TProfile3D
          // Axis mapping:
          // X: Centrality
          // Y: Type (0.5 for Sin, 1.5 for Cos)
          // Z: Harmonic index (ishift - 0.5 maps to bin 1, 2, etc.)

          int binSin = hShiftZNA->FindBin(cent, 0.5, (double)ishift - 0.5);
          int binCos = hShiftZNA->FindBin(cent, 1.5, (double)ishift - 0.5);

          double coeffSin = hShiftZNA->GetBinContent(binSin);
          double coeffCos = hShiftZNA->GetBinContent(binCos);

          // Fourier flattening formula:
          // DeltaPsi = sum( (2/k) * ( <cos>*sin(k*psi) - <sin>*cos(k*psi) ) )
          // Note: signs depend on definition, this matches the standard correction logic
          deltaPsi += (2.0 / ishift) * (-coeffSin * TMath::Cos(ishift * psiZNA) + coeffCos * TMath::Sin(ishift * psiZNA));
        }

        // DEBUG: Print only if shift is actually happening for first few events
        static int debugPrintCount = 0;
        constexpr int maxDebugPrints = 10;
        constexpr double psiTolerance = 1e-6;

        if (debugPrintCount < maxDebugPrints && std::abs(deltaPsi) > psiTolerance) {
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
        for (int ishift = 1; ishift <= nShift; ishift++) {
          // Fill Sin component (Y = 0.5)
          gCurrentShiftProfileZNA->Fill(cent, 0.5, (double)ishift - 0.5, TMath::Sin(ishift * psiZNA));
          // Fill Cos component (Y = 1.5)
          gCurrentShiftProfileZNA->Fill(cent, 1.5, (double)ishift - 0.5, TMath::Cos(ishift * psiZNA));
        }
      }

      // Fill final analysis histogram with the best available Psi (Raw or Corrected)
      gCurrentPsiZNA->Fill(psiZNA);
    }

    // -------- ZNC --------
    if (isZNChit) {
      double qx = zdc.zncQx();
      double qy = zdc.zncQy();

      qxZNCrec = qx;
      qyZNCrec = qy;

      //   LOGF(info, "Qx init = %f", qxZNCrec);

      // Iterate through steps using cached vector
      for (int step = 1; step <= calibrationStep; step++) {

        int cacheIdx = step - 1;
        if (cacheIdx >= static_cast<int>(mCalibCache.size()))
          continue;

        //     LOGF(info, "Go to step = %d", step);

        const auto& calib = mCalibCache[cacheIdx];

        // Apply 5D Base calibration
        if (calib.hMeanQxZNC && calib.hMeanQyZNC) {
          qxZNCrec -= getMeanQFromMap(calib.hMeanQxZNC, cent, vx, vy, vz);
          qyZNCrec -= getMeanQFromMap(calib.hMeanQyZNC, cent, vx, vy, vz);

          //  LOGF(info, "Go to base calibration step = %d", step);
          //  LOGF(info, "Qx after base calibration step %d = %f", step, qxZNCrec);
        }

        if ((step != calibrationStep) || ifFineCalibration) {

          // Apply 1D Refine calibration
          qxZNCrec -= getMeanQ1D(calib.hMeanQxCentZNC, cent);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyCentZNC, cent);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVzZNC, vz);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVzZNC, vz);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVxZNC, vx);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVxZNC, vx);

          qxZNCrec -= getMeanQ1D(calib.hMeanQxVyZNC, vy);
          qyZNCrec -= getMeanQ1D(calib.hMeanQyVyZNC, vy);
          // LOGF(info, "Go to refine calibration step = %d", step);
          //   LOGF(info, "Qx after fine calibration step %d = %f", step, qxZNCrec);
        }
      }

      double valuesQxZNC[5] = {cent, vx, vy, vz, qxZNCrec};
      double valuesQyZNC[5] = {cent, vx, vy, vz, qyZNCrec};

      //    LOGF(info, "Qx cal = %f", qxZNCrec);

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
      double psiZNC = TMath::ATan2(qyZNCrec, qxZNCrec);

      // Apply Correction (Read Mode)
      // Checks if correction is enabled AND if the map from CCDB was loaded successfully
      if (ifShiftCorrection && hShiftZNC) {
        double deltaPsi = 0.0;

        // Loop over harmonics (usually 1 to 10)
        for (int ishift = 1; ishift <= nShift; ishift++) {
          // Retrieve coefficients from TProfile3D
          // Axis mapping:
          // X: Centrality
          // Y: Type (0.5 for Sin, 1.5 for Cos)
          // Z: Harmonic index (ishift - 0.5 maps to bin 1, 2, etc.)

          int binSin = hShiftZNC->FindBin(cent, 0.5, (double)ishift - 0.5);
          int binCos = hShiftZNC->FindBin(cent, 1.5, (double)ishift - 0.5);

          double coeffSin = hShiftZNC->GetBinContent(binSin);
          double coeffCos = hShiftZNC->GetBinContent(binCos);

          // Fourier flattening formula:
          // DeltaPsi = sum( (2/k) * ( <cos>*sin(k*psi) - <sin>*cos(k*psi) ) )
          // Note: signs depend on definition, this matches the standard correction logic
          deltaPsi += (2.0 / ishift) * (-coeffSin * TMath::Cos(ishift * psiZNC) + coeffCos * TMath::Sin(ishift * psiZNC));
        }

        // Apply the calculated shift
        psiZNC += deltaPsi;

        // Wrap angle to [-pi, pi] range
        psiZNC = RecoDecay::constrainAngle(psiZNC, -o2::constants::math::PI);
      }

      // Fill Shift Profiles (Write Mode)
      // Used to generate calibration for the next step or to verify correction (QA)
      if (fillShiftHistos && gCurrentShiftProfileZNC) {
        for (int ishift = 1; ishift <= nShift; ishift++) {
          // Fill Sin component (Y = 0.5)
          gCurrentShiftProfileZNC->Fill(cent, 0.5, (double)ishift - 0.5, TMath::Sin(ishift * psiZNC));
          // Fill Cos component (Y = 1.5)
          gCurrentShiftProfileZNC->Fill(cent, 1.5, (double)ishift - 0.5, TMath::Cos(ishift * psiZNC));
        }
      }

      // Fill final analysis histogram with the best available Psi (Raw or Corrected)
      gCurrentPsiZNC->Fill(psiZNC);
    }

    if (isZNAhit && isZNChit) {
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
