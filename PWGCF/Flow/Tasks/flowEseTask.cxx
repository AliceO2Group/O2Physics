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

/// \author Junlee Kim (jikim1290@gmail.com)
/// \file flowEseTask.cxx
/// \brief Task for flow and event shape engineering correlation with other observation.
/// \since 2023-05-15
/// \version 1.0

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TF1.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TObject.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include <TString.h>
#include <TVector2.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;

struct FlowEseTask {
  //  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults, aod::Qvectors, aod::QvectorFT0CVecs, aod::QvectorTPCposVecs, aod::QvectorTPCnegVecs, aod::QvectorTPCallVecs>;
  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::TrackSelectionExtension>;
  using V0TrackCandidate = aod::V0Datas;

  HistogramRegistry histos{
    "histos",
    {},
    OutputObjHandlingPolicy::AnalysisObject};

  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL",
                                     "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> ccdbNoLaterThan{"ccdbNoLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::ccdb::CcdbApi ccdbApi;

  Configurable<float> cfgCentSel{"cfgCentSel", 80., "Centrality selection"};
  Configurable<int> cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 1: FT0C, 2: FT0M"};

  Configurable<bool> cfgPVSel{"cfgPVSel", false, "Additional PV selection flag for syst"};
  Configurable<float> cfgPV{"cfgPV", 8.0, "Additional PV selection range for syst"};
  Configurable<bool> cfgAddEvtSelPileup{"cfgAddEvtSelPileup", false, "flag for additional pileup selection"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};

  Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
  Configurable<float> cfgDCAPrToPVMin{"cfgDCAPrToPVMin", 0.05, "minimum DCA to PV for proton track"};
  Configurable<float> cfgDCAPiToPVMin{"cfgDCAPiToPVMin", 0.1, "minimum DCA to PV for pion track"};
  Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
  Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};

  Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
  Configurable<float> cfgV0EtaMin{"cfgV0EtaMin", -0.5, "maximum rapidity"};
  Configurable<float> cfgV0EtaMax{"cfgV0EtaMax", 0.5, "maximum rapidity"};
  Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};

  Configurable<float> cfgMinPt{"cfgMinPt", 0.15, "Minimum transverse momentum for track"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8, "Maximum pseudorapidiy for charged track"};

  Configurable<bool> cfgQAv0{"cfgQAv0", false, "QA plot"};

  Configurable<int> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
  Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
  Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
  Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
  Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};

  Configurable<int> cfgnMods{"cfgnMods", 1, "The number of modulations of interest starting from 2"};
  Configurable<int> cfgNQvec{"cfgNQvec", 7, "The number of total Qvectors for looping over the task"};
  Configurable<int> cfgEseHarmonic{"cfgEseHarmonic", 2, "Harmonic used only for ESE analysis (2 or 3)"};

  Configurable<std::string> cfgQvecDetName{"cfgQvecDetName", "FT0C", "The name of detector to be analyzed"};
  Configurable<std::string> cfgQvecRefAName{"cfgQvecRefAName", "TPCpos", "The name of detector for reference A"};
  Configurable<std::string> cfgQvecRefBName{"cfgQvecRefBName", "TPCneg", "The name of detector for reference B"};

  Configurable<bool> cfgPhiDepStudy{"cfgPhiDepStudy", false, "cfg for phi dependent study"};
  Configurable<bool> cfgUSESP{"cfgUSESP", false, "cfg for sp"};
  Configurable<float> cfgPhiDepSig{"cfgPhiDepSig", 0.2, "cfg for significance on phi dependent study"};

  Configurable<bool> cfgShiftCorr{"cfgShiftCorr", false, "additional shift correction"};
  Configurable<bool> cfgShiftCorrDef{"cfgShiftCorrDef", false, "additional shift correction definition"};
  Configurable<std::string> cfgShiftPath{"cfgShiftPath", "Users/j/junlee/Qvector/QvecCalib/Shift", "Path for Shift"};

  Configurable<bool> cfgEffCor{"cfgEffCor", false, "flag to apply efficiency correction"};
  Configurable<std::string> cfgEffCorPath{"cfgEffCorPath", "", "path for pseudo efficiency correction"};

  Configurable<bool> cfgAccCor{"cfgAccCor", false, "flag to apply acceptance correction"};
  Configurable<std::string> cfgAccCorPath{"cfgAccCorPath", "", "path for pseudo acceptance correction"};

  Configurable<bool> cfgCalcCum{"cfgCalcCum", false, "flag to calculate cumulants of cossin"};
  Configurable<bool> cfgCalcCum1{"cfgCalcCum1", false, "flag to calculate cumulants of coscos"};

  Configurable<bool> cfgRapidityDep{"cfgRapidityDep", false, "flag for rapidity dependent study"};
  Configurable<bool> cfgAccAzimuth{"cfgAccAzimuth", false, "flag for azimuth closure study"};

  Configurable<bool> cfgFullCheck{"cfgFullCheck", true, "flag for full hist"};
  Configurable<bool> cfgMultCor{"cfgMultCor", false, "flag for different Mult choice"};
  Configurable<std::vector<float>> cfgQ2PercentileCuts{"cfgQ2PercentileCuts",
                                                       {63.633132f, 92.417875f, 116.584440f, 139.217742f, 161.810257f, 185.609391f, 212.240958f, 244.765093f, 292.045106f,
                                                        57.593461f, 82.544764f, 102.874824f, 121.420485f, 139.519402f, 158.204340f, 178.728019f, 203.334468f, 238.370537f,
                                                        49.132498f, 70.010553f, 86.825107f, 102.036480f, 116.776740f, 131.910217f, 148.447115f, 168.145964f, 195.993568f,
                                                        38.966476f, 55.866796f, 69.630500f, 82.179489f, 94.412155f, 107.016988f, 120.821718f, 137.304969f, 160.612468f,
                                                        29.476997f, 42.628632f, 53.544991f, 63.651863f, 73.628370f, 84.009092f, 95.488749f, 109.295584f, 128.965135f,
                                                        21.658770f, 31.478901f, 39.740751f, 47.485032f, 55.216668f, 63.362554f, 72.466583f, 83.551694f, 99.548536f,
                                                        15.370223f, 22.386745f, 28.322567f, 33.920535f, 39.547603f, 45.520061f, 52.251446f, 60.535625f, 72.660540f,
                                                        10.351031f, 15.094899f, 19.114684f, 22.915848f, 26.751668f, 30.837985f, 35.471803f, 41.208444f, 49.696178f},
                                                       "p10-p90 q2 cuts in each 10% centrality interval from 0 to 80%"};
  Configurable<std::vector<float>> cfgQ3PercentileCuts{"cfgQ3PercentileCuts", {}, "p10-p90 q3 cuts in each 10% centrality interval from 0 to 80%"};

  ConfigurableAxis massAxis{"massAxis", {30, 1.1, 1.13}, "Invariant mass axis"};
  ConfigurableAxis ptAxis{"ptAxis", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "Transverse momentum bins"};
  ConfigurableAxis centAxis{"centAxis", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality interval"};
  ConfigurableAxis cosAxis{"cosAxis", {110, -1.05, 1.05}, "Cosine axis"};
  ConfigurableAxis rapAxis{"rapAxis", {10, -0.5, 0.5}, "Rapidity axis"};
  ConfigurableAxis qqAxis{"qqAxis", {100, -0.1, 0.1}, "qq axis"};
  ConfigurableAxis multAxis{"multAxis", {300, 0, 2700}, "multiplicity"};
  ConfigurableAxis q2QaAxis{"q2QaAxis", {8000, 0.0, 800.0}, "q2 axis for QA"};

  static constexpr float MinAmplitudeThreshold = 1e-5f;
  static constexpr int ShiftLevel = 10;
  static constexpr int LambdaId = 3122;
  static constexpr int SecondHarmonic = 2;
  static constexpr int ThirdHarmonic = 3;
  static constexpr std::array<int, 4> CorrLevel = {2, 3, 4, 1};
  static constexpr std::array<float, 10> CentBoundaries = {0.0f, 3.49f, 4.93f, 6.98f, 8.55f, 9.87f, 11.0f, 12.1f, 13.1f, 14.0f};
  static constexpr std::array<float, 9> CentValues = {2.5f, 7.5f, 15.0f, 25.0f, 35.0f, 45.0f, 55.0f, 65.0f, 75.0f};
  static constexpr float EtaAcceptance = 0.8f;
  static constexpr float CentUpperLimit = 80.0f;
  static constexpr int NEseCentBins = 8;
  static constexpr int NEseGroups = 10;
  static constexpr int NEseCutsPerCentBin = NEseGroups - 1;
  static constexpr float EseCentMin = 0.0f;
  static constexpr float EseCentMax = 80.0f;
  static constexpr float EseCentBinWidth = 10.0f;

  TH2* histEventCountEseGroup = nullptr;
  std::unordered_map<std::string, TObject*> histEse;

  EventPlaneHelper helperEP;

  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;

  int detId = -1;
  int refAId = -1;
  int refBId = -1;

  int qvecDetInd = -1;
  int qvecRefAInd = -1;
  int qvecRefBInd = -1;

  float centrality = -1.0f;

  double angle = 0.0;
  double psi = 0.0;
  double relphi = 0.0;
  double productPhi = 0.0;

  int currentRunNumber = -999;
  int lastRunNumber = -999;
  std::vector<TProfile3D*> shiftprofile;
  TProfile2D* effMap = nullptr;
  TProfile2D* accMap = nullptr;

  std::string fullCCDBShiftCorrPath;

  template <typename T>
  int getDetId(const T& name)
  {
    if (name.value == "FT0A") {
      return 1;
    }
    if (name.value == "FT0M") {
      return 2;
    }
    if (name.value == "FV0A") {
      return 3;
    }
    if (name.value == "TPCpos") {
      return 4;
    }
    if (name.value == "TPCneg") {
      return 5;
    }
    if (name.value == "TPCall") {
      return 6;
    }
    return 0;
  }

  int eseCentBin(float cent) const
  {
    if (cent < EseCentMin || cent >= EseCentMax) {
      return -1;
    }
    return static_cast<int>(cent / EseCentBinWidth);
  }

  std::string eseGroupSuffix(int group) const
  {
    if (group < 0 || group >= NEseGroups) {
      return {};
    }
    return Form("q%dp%02d_%02d", cfgEseHarmonic.value, group * 10, (group + 1) * 10);
  }

  const std::vector<float>& esePercentileCuts() const
  {
    return cfgEseHarmonic.value == SecondHarmonic ? cfgQ2PercentileCuts.value : cfgQ3PercentileCuts.value;
  }

  template <typename T>
  void addEseHistogram(const std::string& name, HistType type, const std::vector<AxisSpec>& axes)
  {
    histEse[name] = histos.add<T>(name, "", type, axes).get();
  }

  template <typename T>
  T* getEseHistogram(const std::string& name)
  {
    const auto histogram = histEse.find(name);
    if (histogram == histEse.end()) {
      LOGF(fatal, "Could not find ESE histogram %s", name.c_str());
      return nullptr;
    }
    return static_cast<T*>(histogram->second);
  }

  int eseGroup(float cent, double qn) const
  {
    const int centBin = eseCentBin(cent);
    if (centBin < 0 || !std::isfinite(qn)) {
      return -1;
    }
    const auto& cuts = esePercentileCuts();
    const int offset = centBin * NEseCutsPerCentBin;
    for (int iCut = 0; iCut < NEseCutsPerCentBin; ++iCut) {
      if (qn < cuts.at(offset + iCut)) {
        return iCut;
      }
    }
    return NEseGroups - 1;
  }

  template <typename TCollision>
  double getEseQ(TCollision const& collision)
  {
    const int harmonicIndex = cfgEseHarmonic.value - SecondHarmonic;
    if (collision.qvecFT0CReVec().size() <= static_cast<std::size_t>(harmonicIndex) || collision.qvecFT0CImVec().size() <= static_cast<std::size_t>(harmonicIndex)) {
      LOGF(fatal, "FT0C Q-vector table does not contain harmonic %d", cfgEseHarmonic.value);
    }
    const double qx = collision.qvecFT0CReVec()[harmonicIndex];
    const double qy = collision.qvecFT0CImVec()[harmonicIndex];
    if (cfgMultCor) {
      return std::sqrt(qx * qx + qy * qy) * collision.sumAmplFT0C() / std::sqrt(collision.multFT0C());
    }
    return std::sqrt(qx * qx + qy * qy) * std::sqrt(collision.sumAmplFT0C());
  }

  void init(o2::framework::InitContext&)
  {
    if (cfgEseHarmonic.value != SecondHarmonic && cfgEseHarmonic.value != ThirdHarmonic) {
      LOGF(fatal, "cfgEseHarmonic must be 2 or 3, got %d", cfgEseHarmonic.value);
    }
    const auto& eseCuts = esePercentileCuts();
    if (eseCuts.size() != NEseCentBins * NEseCutsPerCentBin) {
      LOGF(fatal, "cfgQ%dPercentileCuts must contain %d values, got %d", cfgEseHarmonic.value, NEseCentBins * NEseCutsPerCentBin, static_cast<int>(eseCuts.size()));
    }
    for (int iCent = 0; iCent < NEseCentBins; ++iCent) {
      const int offset = iCent * NEseCutsPerCentBin;
      for (int iCut = 0; iCut < NEseCutsPerCentBin; ++iCut) {
        const float cut = eseCuts.at(offset + iCut);
        if (!std::isfinite(cut) || (iCut > 0 && cut <= eseCuts.at(offset + iCut - 1))) {
          LOGF(fatal, "Invalid q%d percentile cut at centrality bin %d, cut %d", cfgEseHarmonic.value, iCent, iCut);
        }
      }
    }

    AxisSpec centQaAxis = {80, 0.0, 80.0};
    AxisSpec pVzQaAxis = {300, -15.0, 15.0};
    AxisSpec epAxis = {6, 0.0, o2::constants::math::TwoPI};
    AxisSpec epQaAxis = {100, -1.0 * o2::constants::math::PI, o2::constants::math::PI};

    AxisSpec pidAxis = {100, -10, 10};
    AxisSpec vertexAxis = {100, -20, 20};

    AxisSpec shiftAxis = {10, 0, 10, "shift"};
    AxisSpec basisAxis = {20, 0, 20, "basis"};
    AxisSpec eseGroupAxis = {10, 0.0, 10.0, Form("q_{%d} percentile group", cfgEseHarmonic.value)};

    histos.add("histQvecCent", "", {HistType::kTH2F, {q2QaAxis, centQaAxis}});
    histEventCountEseGroup = histos.add<TH2>(Form("histEventCountQ%dGroup", cfgEseHarmonic.value), "", HistType::kTH2F, {centQaAxis, eseGroupAxis}).get();
    histos.add(Form("histVertex"), "", {HistType::kTHnSparseF, {vertexAxis, vertexAxis, vertexAxis, centAxis}});
    for (int iGroup = 0; iGroup < NEseGroups; ++iGroup) {
      const auto suffix = eseGroupSuffix(iGroup);
      addEseHistogram<THnSparse>(Form("histV%d_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {centAxis, ptAxis, cosAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_lambda_cos_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_alambda_cos_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_lambda_cos2_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_alambda_cos2_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_lambda_cossin_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
      addEseHistogram<THnSparse>(Form("psi%d/h_alambda_cossin_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis});
    }
    histos.add("QA/CentDist", "", {HistType::kTH1F, {centQaAxis}});
    histos.add("QA/PVzDist", "", {HistType::kTH1F, {pVzQaAxis}});

    for (auto i = 2; i < cfgnMods + 2; i++) {
      histos.add(Form("psi%d/h_lambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_lambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      histos.add(Form("psi%d/h_alambda_cos2", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, epAxis}});
      if (cfgRapidityDep) {
        histos.add(Form("psi%d/h_lambda_cos2_rap", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, rapAxis}});
        histos.add(Form("psi%d/h_alambda_cos2_rap", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis, rapAxis}});
      }

      histos.add(Form("psi%d/h_lambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_cossin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_lambda_cossin_SP", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_cossin_SP", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      if (cfgAccAzimuth) {
        histos.add(Form("psi%d/h_lambda_coscos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
        histos.add(Form("psi%d/h_alambda_coscos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      }

      histos.add(Form("psi%d/h_lambda_vncos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_lambda_vnsin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_vncos", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add(Form("psi%d/h_alambda_vnsin", i), "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }
    histos.add("QA/ptspec_l", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspec_al", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspecCor_l", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});
    histos.add("QA/ptspecCor_al", "", {HistType::kTH3F, {massAxis, ptAxis, centAxis}});

    if (cfgCalcCum) {
      histos.add("psi2/QA/cosTheta_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }

    if (cfgCalcCum1) {
      histos.add("psi2/QA/cosTheta_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPhi_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_l", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosTheta_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/cosPhi_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_cosPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPhi_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/cosTheta_sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});

      histos.add("psi2/QA/sinPsi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
      histos.add("psi2/QA/sinPhi_al", "", {HistType::kTHnSparseF, {massAxis, ptAxis, cosAxis, centAxis}});
    }

    if (cfgQAv0) {

      histos.add("QA/nsigma_tpc_pt_ppr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_ppi", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpr", "", {HistType::kTH2F, {ptAxis, pidAxis}});
      histos.add("QA/nsigma_tpc_pt_mpi", "", {HistType::kTH2F, {ptAxis, pidAxis}});

      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/QA/EP_Det", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefA", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_RefB", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/qqAxis_Det_RefA_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_Det_RefB_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_RefA_RefB_xx", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});

        histos.add(Form("psi%d/QA/qqAxis_Det_RefA_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_Det_RefB_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});
        histos.add(Form("psi%d/QA/qqAxis_RefA_RefB_yy", i), "", {HistType::kTH2F, {centQaAxis, qqAxis}});

        histos.add(Form("psi%d/QA/EPRes_Det_RefA", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_Det_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_RefA_RefB", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EP_FT0C_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});
        histos.add(Form("psi%d/QA/EP_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, epQaAxis}});

        histos.add(Form("psi%d/QA/EPRes_FT0C_FT0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0C_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
        histos.add(Form("psi%d/QA/EPRes_FT0A_FV0A_shifted", i), "", {HistType::kTH2F, {centQaAxis, cosAxis}});
      }
      for (int iGroup = 0; iGroup < NEseGroups; ++iGroup) {
        const auto suffix = eseGroupSuffix(iGroup);
        addEseHistogram<TH2>(Form("psi%d/QA/EPRes_Det_RefA_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTH2F, {centQaAxis, cosAxis});
        addEseHistogram<TH2>(Form("psi%d/QA/EPRes_Det_RefB_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTH2F, {centQaAxis, cosAxis});
        addEseHistogram<TH2>(Form("psi%d/QA/EPRes_RefA_RefB_%s", cfgEseHarmonic.value, suffix.c_str()), HistType::kTH2F, {centQaAxis, cosAxis});
      }
    }

    if (doprocessMcItsTpc) {
      histos.add("hImpactParameter", "Impact parameter", kTH1F, {{200, 0.0f, 20.0f}});
      histos.add("hEventPlaneAngle", "hEventPlaneAngle", kTH1F, {{200, -1.0 * o2::constants::math::TwoPI, 1.0 * o2::constants::math::TwoPI}});
      histos.add("hEventPlaneAngleRec", "hEventPlaneAngleRec", kTH1F, {{200, -1.0 * o2::constants::math::TwoPI, 1.0 * o2::constants::math::TwoPI}});
      histos.add("hNchVsImpactParameter", "hNchVsImpactParameter", kTH2F, {{200, 0.0f, 20.0f}, {500, -0.5f, 5000.5f}});
      histos.add("hSparseMCGenWeight", "hSparseMCGenWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0f, o2::constants::math::PI}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecWeight", "hSparseMCRecWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0f, o2::constants::math::PI}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
      histos.add("hSparseMCRecAllTrackWeight", "hSparseMCRecAllTrackWeight", HistType::kTHnSparseF, {centAxis, {36, 0.0, o2::constants::math::PI}, {50, 0.0f, 1}, ptAxis, {8, -0.8, 0.8}});
    }

    if (cfgShiftCorrDef) {
      for (auto i = 2; i < cfgnMods + 2; i++) {
        histos.add(Form("psi%d/ShiftFIT", i), "", kTProfile3D, {centQaAxis, basisAxis, shiftAxis});
      }
      if (cfgEseHarmonic.value >= cfgnMods.value + 2) {
        histos.add(Form("psi%d/ShiftFIT", cfgEseHarmonic.value), "", kTProfile3D, {centQaAxis, basisAxis, shiftAxis});
      }
    }

    detId = getDetId(cfgQvecDetName);
    refAId = getDetId(cfgQvecRefAName);
    refBId = getDetId(cfgQvecRefBName);

    if (detId == refAId || detId == refBId || refAId == refBId) {
      LOGF(info, "Wrong detector configuration \n The FT0C will be used to get Q-Vector \n The TPCpos and TPCneg will be used as reference systems");
      detId = 0;
      refAId = 4;
      refBId = 5;
    }

    fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
    fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
    fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);

    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
  }

  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  ROOT::Math::PxPyPzMVector protonVec, pionVec, LambdaVec, protonBoostedVec, pionBoostedVec;

  template <typename TCollision>
  bool eventSelected(TCollision const& collision)
  {
    if (!collision.sel8()) {
      return 0;
    }

    if (cfgCentSel < centrality) {
      return 0;
    }
    /*
        auto multNTracksPV = collision.multNTracksPV();
        if (multNTracksPV < fMultPVCutLow->Eval(centrality)) {
          return 0;
        }
        if (multNTracksPV > fMultPVCutHigh->Eval(centrality)) {
          return 0;
        }
    */
    if (!collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return 0;
    }
    if (!collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return 0;
    }
    if (cfgPVSel && std::abs(collision.posZ()) > cfgPV) {
      return 0;
    }
    if (cfgAddEvtSelPileup && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return 0;
    }
    if (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy) {
      return 0;
    }

    return 1;
  } // event selection

  template <typename TCollision, typename V0>
  bool selectionV0(TCollision const& collision, V0 const& candidate, int lambdaTag)
  {
    if (candidate.v0radius() < cfgv0radiusMin) {
      return false;
    }
    if (lambdaTag) {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPrToPVMin) {
        return false;
      }
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPiToPVMin) {
        return false;
      }
    } else {
      if (std::abs(candidate.dcapostopv()) < cfgDCAPiToPVMin) {
        return false;
      }
      if (std::abs(candidate.dcanegtopv()) < cfgDCAPrToPVMin) {
        return false;
      }
    }
    if (candidate.v0cosPA() < cfgv0CosPA) {
      return false;
    }
    if (std::abs(candidate.dcaV0daughters()) > cfgDCAV0Dau) {
      return false;
    }
    if (candidate.pt() < cfgV0PtMin) {
      return false;
    }
    if (candidate.yLambda() < cfgV0EtaMin) {
      return false;
    }
    if (candidate.yLambda() > cfgV0EtaMax) {
      return false;
    }
    if (candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda > cfgV0LifeTime) {
      return false;
    }

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (track.tpcNClsFound() < cfgDaughTPCnclsMin) {
      return false;
    }
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgDaughPIDCutsTPCPr) {
      return false;
    }
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgDaughPIDCutsTPCPi) {
      return false;
    }
    if (track.eta() > cfgDaughEtaMax) {
      return false;
    }
    if (track.eta() < cfgDaughEtaMin) {
      return false;
    }
    if (pid == 0 && track.pt() < cfgDaughPrPt) {
      return false;
    }
    if (pid == 1 && track.pt() < cfgDaughPiPt) {
      return false;
    }

    return true;
  }

  double safeATan2(double y, double x)
  {
    if (x != 0) {
      return std::atan2(y, x);
    }
    if (y == 0) {
      return 0;
    }
    if (y > 0) {
      return o2::constants::math::PIHalf;
    }
    return -o2::constants::math::PIHalf;
  }

  template <typename TrackType>
  bool selectionTrack(TrackType const& track)
  {
    if (track.pt() < cfgMinPt) {
      return false;
    }
    if (std::abs(track.eta()) > cfgMaxEta) {
      return false;
    }
    if (!track.passedITSNCls()) {
      return false;
    }
    if (!track.passedITSChi2NDF()) {
      return false;
    }
    if (!track.passedITSHits()) {
      return false;
    }
    if (!track.passedTPCCrossedRowsOverNCls()) {
      return false;
    }
    if (!track.passedTPCChi2NDF()) {
      return false;
    }
    if (!track.passedDCAxy()) {
      return false;
    }
    if (!track.passedDCAz()) {
      return false;
    }

    return true;
  }

  template <typename TCollision>
  void fillShiftCorrection(TCollision const& collision, int nmode)
  {
    qvecDetInd = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefAInd = refAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefBInd = refBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    for (int ishift = 1; ishift <= ShiftLevel; ishift++) {
      if (nmode == CorrLevel[0]) {
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 0.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 1.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 2.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 3.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi2/ShiftFIT"), centrality, 4.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi2/ShiftFIT"), centrality, 5.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == CorrLevel[1]) {
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 0.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 1.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 2.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 3.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi3/ShiftFIT"), centrality, 4.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi3/ShiftFIT"), centrality, 5.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
      } else if (nmode == CorrLevel[2]) {
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 0.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 1.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 2.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 3.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode)));

        histos.fill(HIST("psi4/ShiftFIT"), centrality, 4.5, ishift - 0.5, std::sin(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
        histos.fill(HIST("psi4/ShiftFIT"), centrality, 5.5, ishift - 0.5, std::cos(ishift * static_cast<float>(nmode) * std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode)));
      }
    }
  }

  template <typename TCollision>
  void fillEPQA(TCollision const& collision, int nmode)
  {
    qvecDetInd = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefAInd = refAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefBInd = refBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;

    if (collision.qvecAmp()[detId] < MinAmplitudeThreshold || collision.qvecAmp()[refAId] < MinAmplitudeThreshold || collision.qvecAmp()[refBId] < MinAmplitudeThreshold) {
      return;
    }

    if (nmode == CorrLevel[0]) {
      histos.fill(HIST("psi2/QA/EP_Det"), centrality, std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefA"), centrality, std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi2/QA/EP_RefB"), centrality, std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi2/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefAInd]);
      histos.fill(HIST("psi2/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefBInd]);
      histos.fill(HIST("psi2/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[qvecRefAInd] * collision.qvecRe()[qvecRefBInd]);

      histos.fill(HIST("psi2/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefAInd]);
      histos.fill(HIST("psi2/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefBInd]);
      histos.fill(HIST("psi2/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[qvecRefAInd] * collision.qvecIm()[qvecRefBInd]);

      histos.fill(HIST("psi2/QA/EPRes_Det_RefA"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd])));
      histos.fill(HIST("psi2/QA/EPRes_Det_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
      histos.fill(HIST("psi2/QA/EPRes_RefA_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
    } else if (nmode == CorrLevel[1]) {
      histos.fill(HIST("psi3/QA/EP_Det"), centrality, std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefA"), centrality, std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi3/QA/EP_RefB"), centrality, std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi3/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefAInd]);
      histos.fill(HIST("psi3/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefBInd]);
      histos.fill(HIST("psi3/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[qvecRefAInd] * collision.qvecRe()[qvecRefBInd]);

      histos.fill(HIST("psi3/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefAInd]);
      histos.fill(HIST("psi3/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefBInd]);
      histos.fill(HIST("psi3/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[qvecRefAInd] * collision.qvecIm()[qvecRefBInd]);

      histos.fill(HIST("psi3/QA/EPRes_Det_RefA"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd])));
      histos.fill(HIST("psi3/QA/EPRes_Det_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
      histos.fill(HIST("psi3/QA/EPRes_RefA_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
    } else if (nmode == CorrLevel[2]) {
      histos.fill(HIST("psi4/QA/EP_Det"), centrality, std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefA"), centrality, std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode));
      histos.fill(HIST("psi4/QA/EP_RefB"), centrality, std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode));

      histos.fill(HIST("psi4/QA/qqAxis_Det_RefA_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefAInd]);
      histos.fill(HIST("psi4/QA/qqAxis_Det_RefB_xx"), centrality, collision.qvecRe()[qvecDetInd] * collision.qvecRe()[qvecRefBInd]);
      histos.fill(HIST("psi4/QA/qqAxis_RefA_RefB_xx"), centrality, collision.qvecRe()[qvecRefAInd] * collision.qvecRe()[qvecRefBInd]);

      histos.fill(HIST("psi4/QA/qqAxis_Det_RefA_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefAInd]);
      histos.fill(HIST("psi4/QA/qqAxis_Det_RefB_yy"), centrality, collision.qvecIm()[qvecDetInd] * collision.qvecIm()[qvecRefBInd]);
      histos.fill(HIST("psi4/QA/qqAxis_RefA_RefB_yy"), centrality, collision.qvecIm()[qvecRefAInd] * collision.qvecIm()[qvecRefBInd]);

      histos.fill(HIST("psi4/QA/EPRes_Det_RefA"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd])));
      histos.fill(HIST("psi4/QA/EPRes_Det_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
      histos.fill(HIST("psi4/QA/EPRes_RefA_RefB"), centrality, std::cos(std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) - std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd])));
    }

    if (cfgShiftCorr) {
      auto deltapsiFT0C = 0.0;
      auto deltapsiFT0A = 0.0;
      auto deltapsiFV0A = 0.0;

      auto psidefFT0C = std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode);
      auto psidefFT0A = std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode);
      auto psidefFV0A = std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode);
      for (int ishift = 1; ishift <= ShiftLevel; ishift++) {
        auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
        auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
        auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
        auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
        auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
        auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

        deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * std::cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * std::sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
        deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * std::cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * std::sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
        deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * std::cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * std::sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
      }
      if (nmode == CorrLevel[0]) {
        histos.fill(HIST("psi2/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi2/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi2/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi2/QA/EPRes_FT0C_FT0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0C_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi2/QA/EPRes_FT0A_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == CorrLevel[1]) {
        histos.fill(HIST("psi3/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi3/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi3/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi3/QA/EPRes_FT0C_FT0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0C_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi3/QA/EPRes_FT0A_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      } else if (nmode == CorrLevel[2]) {
        histos.fill(HIST("psi4/QA/EP_FT0C_shifted"), centrality, psidefFT0C + deltapsiFT0C);
        histos.fill(HIST("psi4/QA/EP_FT0A_shifted"), centrality, psidefFT0A + deltapsiFT0A);
        histos.fill(HIST("psi4/QA/EP_FV0A_shifted"), centrality, psidefFV0A + deltapsiFV0A);

        histos.fill(HIST("psi4/QA/EPRes_FT0C_FT0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFT0A - deltapsiFT0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0C_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0C + deltapsiFT0C - psidefFV0A - deltapsiFV0A)));
        histos.fill(HIST("psi4/QA/EPRes_FT0A_FV0A_shifted"), centrality, std::cos(static_cast<float>(nmode) * (psidefFT0A + deltapsiFT0A - psidefFV0A - deltapsiFV0A)));
      }
    }
  }

  template <typename TCollision>
  void fillEseEPQA(TCollision const& collision, int eseGroupIndex)
  {
    if (eseGroupIndex < 0 || !cfgQAv0) {
      return;
    }
    const int nmode = cfgEseHarmonic.value;
    const int detIndex = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    const int refAIndex = refAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    const int refBIndex = refBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    if (collision.qvecAmp()[detId] < MinAmplitudeThreshold || collision.qvecAmp()[refAId] < MinAmplitudeThreshold || collision.qvecAmp()[refBId] < MinAmplitudeThreshold) {
      return;
    }

    const double detPhase = std::atan2(collision.qvecIm()[detIndex], collision.qvecRe()[detIndex]);
    const double refAPhase = std::atan2(collision.qvecIm()[refAIndex], collision.qvecRe()[refAIndex]);
    const double refBPhase = std::atan2(collision.qvecIm()[refBIndex], collision.qvecRe()[refBIndex]);
    const auto suffix = eseGroupSuffix(eseGroupIndex);
    auto* histEseEPResDetRefA = getEseHistogram<TH2>(Form("psi%d/QA/EPRes_Det_RefA_%s", cfgEseHarmonic.value, suffix.c_str()));
    auto* histEseEPResDetRefB = getEseHistogram<TH2>(Form("psi%d/QA/EPRes_Det_RefB_%s", cfgEseHarmonic.value, suffix.c_str()));
    auto* histEseEPResRefARefB = getEseHistogram<TH2>(Form("psi%d/QA/EPRes_RefA_RefB_%s", cfgEseHarmonic.value, suffix.c_str()));
    histEseEPResDetRefA->Fill(centrality, std::cos(detPhase - refAPhase));
    histEseEPResDetRefB->Fill(centrality, std::cos(detPhase - refBPhase));
    histEseEPResRefARefB->Fill(centrality, std::cos(refAPhase - refBPhase));
  }

  template <typename TCollision, typename V0, typename TrackType>
  void fillHistograms(TCollision const& collision, V0 const& V0s, TrackType const& track, int nmode, int eseGroupIndex = -1, bool fillRegular = true)
  {
    qvecDetInd = detId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefAInd = refAId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    qvecRefBInd = refBId * 4 + 3 + (nmode - 2) * cfgNQvec * 4;
    const bool fillEse = eseGroupIndex >= 0 && nmode == cfgEseHarmonic.value;

    if (fillEse) {
      const auto suffix = eseGroupSuffix(eseGroupIndex);
      auto* histEseVn = getEseHistogram<THnSparse>(Form("histV%d_%s", cfgEseHarmonic.value, suffix.c_str()));
      const int harmonicIndex = nmode - 2;
      const double esePlane = helperEP.GetEventPlane(collision.qvecFT0CReVec()[harmonicIndex], collision.qvecFT0CImVec()[harmonicIndex], nmode);
      for (const auto& trk : track) {
        if (!selectionTrack(trk)) {
          continue;
        }
        const std::array<double, 3> values = {centrality, trk.pt(), std::cos(static_cast<float>(nmode) * (trk.phi() - esePlane))};
        histEseVn->Fill(values.data());
      }
    }

    if (!fillRegular && !cfgFullCheck) {
      return;
    }
    if (fillRegular) {
      histos.fill(HIST("histVertex"), collision.posX(), collision.posY(), collision.posZ(), collision.centFT0C());
    }

    THnSparse* histEseLambdaCos = nullptr;
    THnSparse* histEseAntiLambdaCos = nullptr;
    THnSparse* histEseLambdaCos2 = nullptr;
    THnSparse* histEseAntiLambdaCos2 = nullptr;
    THnSparse* histEseLambdaCosSin = nullptr;
    THnSparse* histEseAntiLambdaCosSin = nullptr;
    if (fillEse && cfgFullCheck) {
      const auto suffix = eseGroupSuffix(eseGroupIndex);
      histEseLambdaCos = getEseHistogram<THnSparse>(Form("psi%d/h_lambda_cos_%s", cfgEseHarmonic.value, suffix.c_str()));
      histEseAntiLambdaCos = getEseHistogram<THnSparse>(Form("psi%d/h_alambda_cos_%s", cfgEseHarmonic.value, suffix.c_str()));
      histEseLambdaCos2 = getEseHistogram<THnSparse>(Form("psi%d/h_lambda_cos2_%s", cfgEseHarmonic.value, suffix.c_str()));
      histEseAntiLambdaCos2 = getEseHistogram<THnSparse>(Form("psi%d/h_alambda_cos2_%s", cfgEseHarmonic.value, suffix.c_str()));
      histEseLambdaCosSin = getEseHistogram<THnSparse>(Form("psi%d/h_lambda_cossin_%s", cfgEseHarmonic.value, suffix.c_str()));
      histEseAntiLambdaCosSin = getEseHistogram<THnSparse>(Form("psi%d/h_alambda_cossin_%s", cfgEseHarmonic.value, suffix.c_str()));
    }

    for (const auto& v0 : V0s) {
      auto postrack = v0.template posTrack_as<TrackCandidates>();
      auto negtrack = v0.template negTrack_as<TrackCandidates>();

      double nTPCSigmaPosPr = postrack.tpcNSigmaPr();
      double nTPCSigmaNegPi = negtrack.tpcNSigmaPi();

      double nTPCSigmaNegPr = negtrack.tpcNSigmaPr();
      double nTPCSigmaPosPi = postrack.tpcNSigmaPi();

      if (fillRegular && cfgQAv0 && nmode == CorrLevel[0]) {
        histos.fill(HIST("QA/nsigma_tpc_pt_ppr"), postrack.pt(), nTPCSigmaPosPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_ppi"), postrack.pt(), nTPCSigmaPosPi);

        histos.fill(HIST("QA/nsigma_tpc_pt_mpr"), negtrack.pt(), nTPCSigmaNegPr);
        histos.fill(HIST("QA/nsigma_tpc_pt_mpi"), negtrack.pt(), nTPCSigmaNegPi);
      }

      int lambdaTag = 0;
      int aLambdaTag = 0;

      if (isSelectedV0Daughter(postrack, 0) && isSelectedV0Daughter(negtrack, 1)) {
        lambdaTag = 1;
      }
      if (isSelectedV0Daughter(negtrack, 0) && isSelectedV0Daughter(postrack, 1)) {
        aLambdaTag = 1;
      }

      if (lambdaTag == aLambdaTag) {
        continue;
      }

      if (!selectionV0(collision, v0, lambdaTag)) {
        continue;
      }

      if (lambdaTag) {
        protonVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        pionVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
      }
      if (aLambdaTag) {
        protonVec = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        pionVec = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
      }
      LambdaVec = protonVec + pionVec;
      LambdaVec.SetM(massLambda);

      ROOT::Math::Boost boost{LambdaVec.BoostToCM()};
      protonBoostedVec = boost(protonVec);

      angle = protonBoostedVec.Pz() / protonBoostedVec.P();
      psi = safeATan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode);
      relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - psi));
      productPhi = std::sin(static_cast<float>(nmode) * LambdaVec.Phi()) * collision.qvecRe()[qvecDetInd] -
                   std::cos(static_cast<float>(nmode) * LambdaVec.Phi()) * collision.qvecIm()[qvecDetInd];

      if (cfgShiftCorr) {
        auto deltapsiFT0C = 0.0;
        auto deltapsiFT0A = 0.0;
        auto deltapsiFV0A = 0.0;

        auto psidefFT0C = std::atan2(collision.qvecIm()[qvecDetInd], collision.qvecRe()[qvecDetInd]) / static_cast<float>(nmode);
        auto psidefFT0A = std::atan2(collision.qvecIm()[qvecRefAInd], collision.qvecRe()[qvecRefAInd]) / static_cast<float>(nmode);
        auto psidefFV0A = std::atan2(collision.qvecIm()[qvecRefBInd], collision.qvecRe()[qvecRefBInd]) / static_cast<float>(nmode);
        for (int ishift = 1; ishift <= ShiftLevel; ishift++) {
          auto coeffshiftxFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 0.5, ishift - 0.5));
          auto coeffshiftyFT0C = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 1.5, ishift - 0.5));
          auto coeffshiftxFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 2.5, ishift - 0.5));
          auto coeffshiftyFT0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 3.5, ishift - 0.5));
          auto coeffshiftxFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 4.5, ishift - 0.5));
          auto coeffshiftyFV0A = shiftprofile.at(nmode - 2)->GetBinContent(shiftprofile.at(nmode - 2)->FindBin(centrality, 5.5, ishift - 0.5));

          deltapsiFT0C += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0C * std::cos(ishift * static_cast<float>(nmode) * psidefFT0C) + coeffshiftyFT0C * std::sin(ishift * static_cast<float>(nmode) * psidefFT0C)));
          deltapsiFT0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFT0A * std::cos(ishift * static_cast<float>(nmode) * psidefFT0A) + coeffshiftyFT0A * std::sin(ishift * static_cast<float>(nmode) * psidefFT0A)));
          deltapsiFV0A += ((1 / (1.0 * ishift)) * (-coeffshiftxFV0A * std::cos(ishift * static_cast<float>(nmode) * psidefFV0A) + coeffshiftyFV0A * std::sin(ishift * static_cast<float>(nmode) * psidefFV0A)));
        }
        psi += deltapsiFT0C;
        relphi = TVector2::Phi_0_2pi(static_cast<float>(nmode) * (LambdaVec.Phi() - psidefFT0C - deltapsiFT0C));
      }

      if (cfgPhiDepStudy && cfgPhiDepSig * std::abs(std::sin(relphi)) > gRandom->Uniform(0, 1)) {
        continue;
      }

      if (fillRegular) {
        if (lambdaTag) {
          histos.fill(HIST("QA/ptspec_l"), v0.mLambda(), v0.pt(), centrality);
          if (cfgEffCor) {
            histos.fill(HIST("QA/ptspecCor_l"), v0.mLambda(), v0.pt(), centrality,
                        1.0 / effMap->GetBinContent(effMap->GetXaxis()->FindBin(v0.pt()), effMap->GetYaxis()->FindBin(centrality)));
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("QA/ptspec_al"), v0.mAntiLambda(), v0.pt(), centrality);
          if (cfgEffCor) {
            histos.fill(HIST("QA/ptspecCor_al"), v0.mAntiLambda(), v0.pt(), centrality,
                        1.0 / effMap->GetBinContent(effMap->GetXaxis()->FindBin(v0.pt()), effMap->GetYaxis()->FindBin(centrality)));
          }
        }
      }
      double weight = 1.0;
      weight *= cfgEffCor ? 1.0 / effMap->GetBinContent(effMap->GetXaxis()->FindBin(v0.pt()), effMap->GetYaxis()->FindBin(centrality)) : 1.;
      weight *= cfgAccCor ? 1.0 / accMap->GetBinContent(accMap->GetXaxis()->FindBin(v0.pt()), accMap->GetYaxis()->FindBin(v0.yLambda())) : 1.;

      if (fillEse && cfgFullCheck) {
        const double mass = lambdaTag ? v0.mLambda() : v0.mAntiLambda();
        const std::array<double, 4> cosValues = {mass, v0.pt(), angle * weight, centrality};
        const std::array<double, 4> cos2Values = {mass, v0.pt(), angle * angle, centrality};
        const std::array<double, 4> cosSinValues = {mass, v0.pt(), angle * std::sin(relphi) * weight, centrality};
        if (lambdaTag) {
          histEseLambdaCos->Fill(cosValues.data());
          histEseLambdaCos2->Fill(cos2Values.data());
          histEseLambdaCosSin->Fill(cosSinValues.data());
        } else {
          histEseAntiLambdaCos->Fill(cosValues.data());
          histEseAntiLambdaCos2->Fill(cos2Values.data());
          histEseAntiLambdaCosSin->Fill(cosSinValues.data());
        }
      }
      if (!fillRegular) {
        continue;
      }

      double qvecMag = 1.0;
      if (cfgUSESP) {
        qvecMag *= std::sqrt(std::pow(collision.qvecIm()[3 + (nmode - 2) * 28], 2) + std::pow(collision.qvecRe()[3 + (nmode - 2) * 28], 2));
      }

      if (nmode == CorrLevel[0] && cfgFullCheck) { ////////////
        if (lambdaTag) {
          histos.fill(HIST("psi2/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_lambda_cossin_SP"), v0.mLambda(), v0.pt(), angle * productPhi * weight, centrality);
          histos.fill(HIST("psi2/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_lambda_vnsin"), v0.mLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi2/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi2/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }

          if (cfgCalcCum) {
            histos.fill(HIST("psi2/QA/cosTheta_l"), v0.mLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_l"), v0.mLambda(), v0.pt(), std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_l"), v0.mLambda(), v0.pt(), std::cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_l"), v0.mLambda(), v0.pt(), std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_l"), v0.mLambda(), v0.pt(), std::sin(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_l"), v0.mLambda(), v0.pt(), angle * std::cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_l"), v0.mLambda(), v0.pt(), angle * std::cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_l"), v0.mLambda(), v0.pt(), angle * std::sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_l"), v0.mLambda(), v0.pt(), angle * std::sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_sinPsi_l"), v0.mLambda(), v0.pt(), std::cos(v0.phi() * 2.0) * std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_cosPsi_l"), v0.mLambda(), v0.pt(), std::sin(v0.phi() * 2.0) * std::cos(psi * 2.0), centrality);
          }
          if (cfgCalcCum1) {
            histos.fill(HIST("psi2/QA/cosTheta_l"), v0.mLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_l"), v0.mLambda(), v0.pt(), std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_l"), v0.mLambda(), v0.pt(), std::cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_cosPsi_l"), v0.mLambda(), v0.pt(), std::cos(v0.phi() * 2.0) * std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_l"), v0.mLambda(), v0.pt(), angle * std::cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_l"), v0.mLambda(), v0.pt(), angle * std::cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPhi_sinPsi_l"), v0.mLambda(), v0.pt(), std::sin(v0.phi() * 2.0) * std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_l"), v0.mLambda(), v0.pt(), angle * std::sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_l"), v0.mLambda(), v0.pt(), angle * std::sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_l"), v0.mLambda(), v0.pt(), std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_l"), v0.mLambda(), v0.pt(), std::sin(v0.phi() * 2.0), centrality);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi2/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi2/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_alambda_cossin_SP"), v0.mAntiLambda(), v0.pt(), angle * productPhi * weight, centrality);
          histos.fill(HIST("psi2/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi2/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi2/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi2/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }

          if (cfgCalcCum) {
            histos.fill(HIST("psi2/QA/cosTheta_al"), v0.mAntiLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_al"), v0.mAntiLambda(), v0.pt(), std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_al"), v0.mAntiLambda(), v0.pt(), std::cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_al"), v0.mAntiLambda(), v0.pt(), std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_al"), v0.mAntiLambda(), v0.pt(), std::sin(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_al"), v0.mAntiLambda(), v0.pt(), angle * std::cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_al"), v0.mAntiLambda(), v0.pt(), angle * std::cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_al"), v0.mAntiLambda(), v0.pt(), angle * std::sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_al"), v0.mAntiLambda(), v0.pt(), angle * std::sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_sinPsi_al"), v0.mAntiLambda(), v0.pt(), std::cos(v0.phi() * 2.0) * std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_cosPsi_al"), v0.mAntiLambda(), v0.pt(), std::sin(v0.phi() * 2.0) * std::cos(psi * 2.0), centrality);
          }
          if (cfgCalcCum1) {
            histos.fill(HIST("psi2/QA/cosTheta_al"), v0.mAntiLambda(), v0.pt(), angle, centrality);
            histos.fill(HIST("psi2/QA/cosPsi_al"), v0.mAntiLambda(), v0.pt(), std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosPhi_al"), v0.mAntiLambda(), v0.pt(), std::cos(v0.phi() * 2.0), centrality);

            histos.fill(HIST("psi2/QA/cosPhi_cosPsi_al"), v0.mAntiLambda(), v0.pt(), std::cos(v0.phi() * 2.0) * std::cos(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPhi_al"), v0.mAntiLambda(), v0.pt(), angle * std::cos(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_cosPsi_al"), v0.mAntiLambda(), v0.pt(), angle * std::cos(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPhi_sinPsi_al"), v0.mAntiLambda(), v0.pt(), std::sin(v0.phi() * 2.0) * std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPhi_al"), v0.mAntiLambda(), v0.pt(), angle * std::sin(v0.phi() * 2.0), centrality);
            histos.fill(HIST("psi2/QA/cosTheta_sinPsi_al"), v0.mAntiLambda(), v0.pt(), angle * std::sin(psi * 2.0), centrality);

            histos.fill(HIST("psi2/QA/sinPsi_al"), v0.mAntiLambda(), v0.pt(), std::sin(psi * 2.0), centrality);
            histos.fill(HIST("psi2/QA/sinPhi_al"), v0.mAntiLambda(), v0.pt(), std::sin(v0.phi() * 2.0), centrality);
          }
        }
      } else if (nmode == CorrLevel[1]) {
        if (lambdaTag) {
          histos.fill(HIST("psi3/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi3/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_lambda_vnsin"), v0.mLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi3/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi3/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi3/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi3/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi, weight);
          histos.fill(HIST("psi3/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi3/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi3/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi3/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }
        }
      } else if (nmode == CorrLevel[2]) {
        if (lambdaTag) {
          histos.fill(HIST("psi4/h_lambda_cos"), v0.mLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cos2"), v0.mLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_lambda_cossin"), v0.mLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_lambda_vncos"), v0.mLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_lambda_vnsin"), v0.mLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi4/h_lambda_cos2_rap"), v0.mLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi4/h_lambda_coscos"), v0.mLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }
        }
        if (aLambdaTag) {
          histos.fill(HIST("psi4/h_alambda_cos"), v0.mAntiLambda(), v0.pt(), angle * weight, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cos2"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, relphi);
          histos.fill(HIST("psi4/h_alambda_cossin"), v0.mAntiLambda(), v0.pt(), angle * std::sin(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_alambda_vncos"), v0.mAntiLambda(), v0.pt(), qvecMag * std::cos(relphi) * weight, centrality);
          histos.fill(HIST("psi4/h_alambda_vnsin"), v0.mAntiLambda(), v0.pt(), std::sin(relphi), centrality);

          if (cfgRapidityDep) {
            histos.fill(HIST("psi4/h_alambda_cos2_rap"), v0.mAntiLambda(), v0.pt(), angle * angle, centrality, v0.yLambda(), weight);
          }

          if (cfgAccAzimuth) {
            histos.fill(HIST("psi4/h_alambda_coscos"), v0.mAntiLambda(), v0.pt(), angle * std::cos(relphi), centrality, weight);
          }
        }
      } ////////// FIXME: not possible to get histograms using nmode
    }
  }

  void processData(EventCandidates::iterator const& collision,
                   TrackCandidates const& tracks, aod::V0Datas const& V0s,
                   aod::BCsWithTimestamps const&)
  {
    if (cfgCentEst == CorrLevel[3]) {
      centrality = collision.centFT0C();
    } else if (cfgCentEst == CorrLevel[0]) {
      centrality = collision.centFT0M();
    }
    if (!eventSelected(collision)) {
      return;
    }
    histos.fill(HIST("QA/CentDist"), centrality, 1.0);
    histos.fill(HIST("QA/PVzDist"), collision.posZ(), 1.0);
    const double eseQ = getEseQ(collision);
    histos.fill(HIST("histQvecCent"), eseQ, centrality);
    const int eseGroupIndex = eseGroup(centrality, eseQ);
    if (eseGroupIndex >= 0) {
      histEventCountEseGroup->Fill(centrality, eseGroupIndex + 0.5);
    }

    if (cfgShiftCorr) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      currentRunNumber = bc.runNumber();
      if (currentRunNumber != lastRunNumber) {
        shiftprofile.clear();
        const int maxShiftHarmonic = std::max(cfgnMods.value + 1, cfgEseHarmonic.value);
        for (int i = 2; i <= maxShiftHarmonic; i++) {
          fullCCDBShiftCorrPath = cfgShiftPath;
          fullCCDBShiftCorrPath += "/v";
          fullCCDBShiftCorrPath += std::to_string(i);
          auto objshift = ccdb->getForTimeStamp<TProfile3D>(fullCCDBShiftCorrPath, bc.timestamp());
          shiftprofile.push_back(objshift);
        }
        lastRunNumber = currentRunNumber;
      }
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (cfgEffCor) {
      effMap = ccdb->getForTimeStamp<TProfile2D>(cfgEffCorPath.value, bc.timestamp());
    }
    if (cfgAccCor) {
      accMap = ccdb->getForTimeStamp<TProfile2D>(cfgAccCorPath.value, bc.timestamp());
    }
    fillEseEPQA(collision, eseGroupIndex);
    if (cfgShiftCorrDef && cfgEseHarmonic.value >= cfgnMods.value + 2) {
      fillShiftCorrection(collision, cfgEseHarmonic.value);
    }
    if (eseGroupIndex >= 0 && cfgEseHarmonic.value >= cfgnMods.value + 2) {
      fillHistograms(collision, V0s, tracks, cfgEseHarmonic.value, eseGroupIndex, false);
    }
    for (int i = 2; i < cfgnMods + 2; i++) {
      if (cfgShiftCorrDef) {
        fillShiftCorrection(collision, i);
      }
      if (cfgQAv0) {
        fillEPQA(collision, i);
      }
      fillHistograms(collision, V0s, tracks, i, i == cfgEseHarmonic.value ? eseGroupIndex : -1);
    } // FIXME: need to fill different histograms for different harmonic
  }
  PROCESS_SWITCH(FlowEseTask, processData, "Process Event for data", true);

  using RecoTracks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  void processMcItsTpc(aod::McCollision const& mcCollision, soa::Join<aod::McParticles, aod::ParticlesToTracks> const& mcParticles, RecoTracks const&)
  {
    float imp = mcCollision.impactParameter();
    float evPhi = mcCollision.eventPlaneAngle() / 2.0;
    float centclass = -999;
    if (imp >= CentBoundaries[0] && imp < CentBoundaries[1]) {
      centclass = CentValues[0];
    }
    if (imp >= CentBoundaries[1] && imp < CentBoundaries[2]) {
      centclass = CentValues[1];
    }
    if (imp >= CentBoundaries[2] && imp < CentBoundaries[3]) {
      centclass = CentValues[2];
    }
    if (imp >= CentBoundaries[3] && imp < CentBoundaries[4]) {
      centclass = CentValues[3];
    }
    if (imp >= CentBoundaries[4] && imp < CentBoundaries[5]) {
      centclass = CentValues[4];
    }
    if (imp >= CentBoundaries[5] && imp < CentBoundaries[6]) {
      centclass = CentValues[5];
    }
    if (imp >= CentBoundaries[6] && imp < CentBoundaries[7]) {
      centclass = CentValues[6];
    }
    if (imp >= CentBoundaries[7] && imp < CentBoundaries[8]) {
      centclass = CentValues[7];
    }
    if (imp >= CentBoundaries[8] && imp < CentBoundaries[9]) {
      centclass = CentValues[8];
    }

    int nCh = 0;

    if (centclass > 0 && centclass < CentUpperLimit) {
      // event within range
      histos.fill(HIST("hImpactParameter"), imp);
      histos.fill(HIST("hEventPlaneAngle"), evPhi);
      for (auto const& mcParticle : mcParticles) {
        float deltaPhi = mcParticle.phi() - mcCollision.eventPlaneAngle();
        // focus on bulk: e, mu, pi, k, p
        int pdgCode = std::abs(mcParticle.pdgCode());
        if (pdgCode != LambdaId) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (std::abs(mcParticle.eta()) > EtaAcceptance) { // main acceptance
          continue;
        }
        histos.fill(HIST("hSparseMCGenWeight"), centclass, RecoDecay::constrainAngle(deltaPhi, 0, 2), std::pow(std::cos(2.0 * RecoDecay::constrainAngle(deltaPhi, 0, 2)), 2.0), mcParticle.pt(), mcParticle.eta());
        nCh++;
        bool validGlobal = false;
        bool validAny = false;
        if (mcParticle.has_tracks()) {
          auto const& tracks = mcParticle.tracks_as<RecoTracks>();
          for (auto const& track : tracks) {
            if (track.hasTPC() && track.hasITS()) {
              validGlobal = true;
            }
            if (track.hasTPC() || track.hasITS()) {
              validAny = true;
            }
          }
        }
        // if valid global, fill
        if (validGlobal) {
          histos.fill(HIST("hSparseMCRecWeight"), centclass, RecoDecay::constrainAngle(deltaPhi, 0, 2), std::pow(std::cos(2.0 * RecoDecay::constrainAngle(deltaPhi, 0, 2)), 2.0), mcParticle.pt(), mcParticle.eta());
        }
        if (validAny) {
          histos.fill(HIST("hSparseMCRecAllTrackWeight"), centclass, RecoDecay::constrainAngle(deltaPhi, 0, 2), std::pow(std::cos(2.0 * RecoDecay::constrainAngle(deltaPhi, 0, 2)), 2.0), mcParticle.pt(), mcParticle.eta());
          histos.fill(HIST("hEventPlaneAngleRec"), RecoDecay::constrainAngle(deltaPhi, 0, 2));
        }
        // if any track present, fill
      }
    }
    histos.fill(HIST("hNchVsImpactParameter"), imp, nCh);
  }
  PROCESS_SWITCH(FlowEseTask, processMcItsTpc, "Process MC for ITSTPC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<FlowEseTask>(cfgc)};
}
