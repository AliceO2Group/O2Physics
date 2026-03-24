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

/// \file radialFlowDecorr.cxx
/// \brief Analysis task for event-by-event radial-flow decorrelation measurement.
/// \author Somadutta Bhatta

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "FT0Base/Geometry.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Utils.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/PID.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"

#include "TDirectory.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct RadialFlowDecorr {

  static constexpr int KPidPionOne = 1;
  static constexpr int KPidKaonTwo = 2;
  static constexpr int KPidProtonThree = 3;
  static constexpr int KConstTen = 10;

  static constexpr int KnFt0cCell = 96;
  static constexpr int KIntM = 3;
  static constexpr int KIntK = 3;
  static constexpr int KNEta = 17;
  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr int KPiPlus = 211;
  static constexpr int KKPlus = 321;
  static constexpr int KProton = 2212;
  static constexpr float KBinOffset = 0.5f;
  static constexpr float KPhiMin = 0.f;
  static constexpr int KNbinsZvtx = 240;
  static constexpr float KZvtxMin = -12.f;
  static constexpr float KZvtxMax = 12.f;
  static constexpr float KPMin = 0.f;
  static constexpr float KPMax = 10.f;
  static constexpr int KNbinsPt = 200;
  static constexpr float KPtMin = 0.15f;
  static constexpr float KPtMax = 10.f;
  static constexpr float KEtaMin = -1.2f;
  static constexpr float KEtaMax = 1.2f;
  static constexpr int KNbinsPhi = 64;
  static constexpr float KEtaAxisMin = -0.8f;
  static constexpr float KEtaAxisMax = 0.8f;
  static constexpr int KNbinsPtRes = 50;
  static constexpr float KPtResMax = 1.f;
  static constexpr int KNbinsEtaRes = 100;
  static constexpr float KEtaResMax = 0.5f;
  static constexpr int KNbinsVz = 80;
  static constexpr float KVzMin = -40.f;
  static constexpr float KVzMax = 40.f;
  static constexpr float KVzResMax = 20.f;
  static constexpr int KNbinsEtaFine = 20;
  static constexpr float KEtaFineMax = 1.f;
  static constexpr int KNbinsDca = 400;
  static constexpr float KDcaMax = 0.2f;
  static constexpr float KCentMax = 90;

  enum PIDIdx {
    kInclusiveIdx = 0,
    kPiMinusIdx,
    kPiPlusIdx,
    kPiAllIdx,
    kKaMinusIdx,
    kKaPlusIdx,
    kKaAllIdx,
    kAntiPrIdx,
    kPrIdx,
    kPrAllIdx,
    KNsp
  };

  inline static const std::vector<std::string> pidSuffix = {"", "_PiMinus", "_PiPlus", "_PiAll", "_KaMinus", "_KaPlus", "_KaAll", "_AntiPr", "_Pr", "_PrAll"};

  struct PIDMeanSigmaMap {
    static constexpr int MaxCentBins = 100;
    double meanTOF[KNsp][MaxCentBins] = {{0.0}};
    double sigmaTOF[KNsp][MaxCentBins] = {{1.0}}; // Default sigma = 1
    double meanTPC[KNsp][MaxCentBins] = {{0.0}};
    double sigmaTPC[KNsp][MaxCentBins] = {{1.0}}; // Default sigma = 1
  };

  PIDMeanSigmaMap* pidMeanSigmaMap = nullptr;

  enum ECentralityEstimator {
    kCentFT0C = 1,
    kCentFT0M = 2,
    kCentFDDM = 3,
    kCentFV0A = 4
  };
  enum SystemType {
    kPbPb = 1,
    kNeNe = 2,
    kOO = 3,
    kpp = 4
  };
  static constexpr float KinvalidCentrality = -1.0f;
  inline static const std::vector<float> etaLw = {
    -0.8,
    -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
  inline static const std::vector<float> etaUp = {
    0.8,
    -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  Configurable<float> cfgVtxZCut{"cfgVtxZCut", 10.f, "z-vertex range"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "min pT"};
  Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "max pT"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "|η| cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutTracKDcaMaxZ{"cfgCutTracKDcaMaxZ", 2.0f, "Maximum DcaZ"};
  Configurable<float> cfgCutTracKDcaMaxXY{"cfgCutTracKDcaMaxXY", 0.2f, "Maximum DcaZ"};

  Configurable<bool> cfgPtDepDCAxy{"cfgPtDepDCAxy", false, "Use pt-dependent DCAxy cut"};
  Configurable<float> cfgDcaXyP0{"cfgDcaXyP0", 0.0026f, "p0 for DCAxy"};
  Configurable<float> cfgDcaXyP1{"cfgDcaXyP1", 0.005f, "p1 for DCAxy"};
  Configurable<float> cfgDcaXyP2{"cfgDcaXyP2", 1.01f, "p2 for DCAxy"};

  Configurable<bool> cfgPtDepDCAz{"cfgPtDepDCAz", false, "Use pt-dependent DCAz cut"};
  Configurable<float> cfgDcaZP0{"cfgDcaZP0", 0.0026f, "p0 for DCAz"};
  Configurable<float> cfgDcaZP1{"cfgDcaZP1", 0.005f, "p1 for DCAz"};
  Configurable<float> cfgDcaZP2{"cfgDcaZP2", 1.01f, "p2 for DCAz"};

  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
  Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
  Configurable<float> cfgnSigmaOtherParticles{"cfgnSigmaOtherParticles", 3.0f, "PID nSigma cut to remove other particles (default:3)"};
  Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
  Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
  Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};

  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut for inclusive hadron analysis"};
  Configurable<float> cfgCutPtUpperPID{"cfgCutPtUpperPID", 6.0f, "Higher pT cut for identified particle analysis"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<int> cfgNsubsample{"cfgNsubsample", 10, "Number of subsamples"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 1-->FT0C, 2-->FT0M, 3-->FDDM, 4-->FV0A"};
  Configurable<bool> cfgEvSelNoSameBunchPileup{"cfgEvSelNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgIsGoodZvtxFT0VsPV{"cfgIsGoodZvtxFT0VsPV", true, "Good Vertexing cut"};

  Configurable<float> cfgPupnSig{"cfgPupnSig", 6.0f, "Additional Pileup Cut"};
  Configurable<bool> cfgApplySigPupCut{"cfgApplySigPupCut", 0, "nSig Pileup Cut"};
  Configurable<bool> cfgApplyLinPupCut{"cfgApplyLinPupCut", 0, "Lin Pileup Cut"};
  Configurable<float> cfgLinPupParam0{"cfgLinPupParam0", 3.0f, "(Upper) Linear Pileup Cut Const"};
  Configurable<float> cfgLinPupParam1{"cfgLinPupParam1", 3.0f, "(Upper) Linear Pileup Slope"};
  Configurable<float> cfgLinPupParam2{"cfgLinPupParam2", 3.0f, "(Lower) Linear Pileup Cut Const"};
  Configurable<float> cfgLinPupParam3{"cfgLinPupParam3", 3.0f, "(Lower) Linear Pileup Slope"};

  Configurable<int> cfgNchPbMax{"cfgNchPbMax", 4000, "Max Nch range for PbPb collisions"};
  Configurable<int> cfgNchOMax{"cfgNchOMax", 800, "Max Nch range for OO collisions"};

  Configurable<int> cfgSys{"cfgSys", 1, "Efficiency to be used for which system? 1-->PbPb, 2-->NeNe, 3-->OO, 4-->pp"};
  Configurable<bool> cfgFlat{"cfgFlat", false, "Whether to use flattening weights"};
  Configurable<bool> cfgEff{"cfgEff", false, "Whether to use Efficiency weights"};
  Configurable<bool> cfgZDC{"cfgZDC", false, "Whether to use ZDC for pileup histograms"};

  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 5.0, 10, 20, 40, 60, 80, 100}, "centrality axis (percentile)"};

  const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
  const AxisSpec centAxis1Per{100, 0.0, 100.0, "Centrality (%)"};
  AxisSpec nChAxis{1, 0., 1., "Nch", "Nch"};
  AxisSpec nChAxis2{1, 0., 1., "Nch", "Nch"};

  const AxisSpec vzAxis{5, -12.5, 12.5, "Vz"};
  const AxisSpec chgAxis{3, -1.5, 1.5};
  const AxisSpec pTAxis{{0.0, 0.2, 0.4, 0.6, 0.8, 1, 3, 5, 7, 10}, "pT Axis"};
  const AxisSpec etaAxis{{-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, "Eta"};
  const AxisSpec phiAxis{KNbinsPhi, KPhiMin, TwoPI, "#phi"};
  const AxisSpec etaBinAxis{KNEta + 1, -0.5, KNEta + 0.5, "#eta bin Number"};
  const AxisSpec spBinAxis{KNsp + 1, -KBinOffset, static_cast<float>(KNsp) + KBinOffset, "species index Number"};

  const AxisSpec gapAxis{{-1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.85,
                          -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05,
                          0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75,
                          0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55},
                         "Gap"};

  const AxisSpec sumAxis{{-1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.85,
                          -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05,
                          0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75,
                          0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55},
                         "Sum"};

  Configurable<bool> cfgRunMCGetNSig{"cfgRunMCGetNSig", false, "Run MC pass to get mean of Nsig Plots"};
  Configurable<bool> cfgRunGetEff{"cfgRunGetEff", false, "Run MC pass to build efficiency/fake maps"};
  Configurable<bool> cfgRunGetMCFlat{"cfgRunGetMCFlat", false, "Run MC to Get Flattening Weights"};
  Configurable<bool> cfgRunMCMean{"cfgRunMCMean", false, "Run MC mean(pT) & mean(Et)"};
  Configurable<bool> cfgRunMCFluc{"cfgRunMCFluc", false, "Run MC fluctuations (C2, subevent)"};

  Configurable<bool> cfgRunDataGetNSig{"cfgRunDataGetNSig", false, "Run MC pass to get mean of Nsig Plots"};
  Configurable<bool> cfgRunGetDataFlat{"cfgRunGetDataFlat", false, "Run Data Get Flattening Weights"};
  Configurable<bool> cfgRunDataMean{"cfgRunDataMean", false, "Run DATA mean(pT) & mean(Et)"};
  Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};

  Service<ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  struct InternalState {
    std::array<TH3F*, KNsp> hEff{};
    std::array<TH3F*, KNsp> hFake{};
    std::array<THnSparseF*, KNsp> hFlatWeight{};

    std::vector<std::pair<float, float>> mLimitsNchCent;
    float mMinXNchCent = 0, mMaxXNchCent = 0;

    TProfile3D* pmeanTruNchEtabinSpbinStep2 = nullptr;
    TProfile3D* pmeanRecoNchEtabinSpbinStep2 = nullptr;
    TProfile3D* pmeanRecoEffcorrNchEtabinSpbinStep2 = nullptr;

    TProfile3D* pmeanMultTruNchEtabinSpbinStep2 = nullptr;
    TProfile3D* pmeanMultRecoNchEtabinSpbinStep2 = nullptr;
    TProfile3D* pmeanMultRecoEffcorrNchEtabinSpbinStep2 = nullptr;

    TProfile3D* pmeanNchEtabinSpbinStep2 = nullptr;
    TProfile3D* pmeanMultNchEtabinSpbinStep2 = nullptr;

    TProfile* pmeanFT0AmultpvStep2 = nullptr;
    TProfile* pmeanFT0CmultpvStep2 = nullptr;
  } state;
  o2::ft0::Geometry ft0Det;

  template <typename T>
  void fillNSigmaBefCut(const T& track, float cent)
  {
    float pt = track.pt();
    auto sign = track.sign();

    if (sign > 0) {
      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_PiPlus"), cent, pt, track.tpcNSigmaPi());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_PiPlus"), cent, pt, track.tofNSigmaPi());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_PiPlus"), cent, track.tofNSigmaPi(), track.tpcNSigmaPi());

      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_KaPlus"), cent, pt, track.tpcNSigmaKa());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_KaPlus"), cent, pt, track.tofNSigmaKa());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_KaPlus"), cent, track.tofNSigmaKa(), track.tpcNSigmaKa());

      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_Pr"), cent, pt, track.tpcNSigmaPr());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_Pr"), cent, pt, track.tofNSigmaPr());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_Pr"), cent, track.tofNSigmaPr(), track.tpcNSigmaPr());
    } else if (sign < 0) {
      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_PiMinus"), cent, pt, track.tpcNSigmaPi());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_PiMinus"), cent, pt, track.tofNSigmaPi());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_PiMinus"), cent, track.tofNSigmaPi(), track.tpcNSigmaPi());

      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_KaMinus"), cent, pt, track.tpcNSigmaKa());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_KaMinus"), cent, pt, track.tofNSigmaKa());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_KaMinus"), cent, track.tofNSigmaKa(), track.tpcNSigmaKa());

      histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_AntiPr"), cent, pt, track.tpcNSigmaPr());
      histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_AntiPr"), cent, pt, track.tofNSigmaPr());
      histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_AntiPr"), cent, track.tofNSigmaPr(), track.tpcNSigmaPr());
    }
    histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_PiAll"), cent, pt, track.tpcNSigmaPi());
    histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_PiAll"), cent, pt, track.tofNSigmaPi());
    histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_PiAll"), cent, track.tofNSigmaPi(), track.tpcNSigmaPi());

    histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_KaAll"), cent, pt, track.tpcNSigmaKa());
    histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_KaAll"), cent, pt, track.tofNSigmaKa());
    histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_KaAll"), cent, track.tofNSigmaKa(), track.tpcNSigmaKa());

    histos.fill(HIST("h3DnsigmaTpcVsPtBefCut_Cent_PrAll"), cent, pt, track.tpcNSigmaPr());
    histos.fill(HIST("h3DnsigmaTofVsPtBefCut_Cent_PrAll"), cent, pt, track.tofNSigmaPr());
    histos.fill(HIST("h3DnsigmaTpcVsTofBefCut_Cent_PrAll"), cent, track.tofNSigmaPr(), track.tpcNSigmaPr());
  }

  template <typename T>
  void fillNSigmaAftCut(const T& track, float cent, bool isSpecies[])
  {
    float pt = track.pt();
    float tpcPi = track.tpcNSigmaPi();
    float tofPi = track.tofNSigmaPi();

    float tpcKa = track.tpcNSigmaKa();
    float tofKa = track.tofNSigmaKa();

    float tpcPr = track.tpcNSigmaPr();
    float tofPr = track.tofNSigmaPr();

    if (isSpecies[kPiPlusIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_PiPlus"), cent, pt, tpcPi);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_PiPlus"), cent, pt, tofPi);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_PiPlus"), cent, tofPi, tpcPi);
    }
    if (isSpecies[kPiMinusIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_PiMinus"), cent, pt, tpcPi);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_PiMinus"), cent, pt, tofPi);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_PiMinus"), cent, tofPi, tpcPi);
    }
    if (isSpecies[kPiAllIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_PiAll"), cent, pt, tpcPi);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_PiAll"), cent, pt, tofPi);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_PiAll"), cent, tofPi, tpcPi);
    }
    if (isSpecies[kKaPlusIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_KaPlus"), cent, pt, tpcKa);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_KaPlus"), cent, pt, tofKa);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_KaPlus"), cent, tofKa, tpcKa);
    }
    if (isSpecies[kKaMinusIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_KaMinus"), cent, pt, tpcKa);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_KaMinus"), cent, pt, tofKa);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_KaMinus"), cent, tofKa, tpcKa);
    }
    if (isSpecies[kKaAllIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_KaAll"), cent, pt, tpcKa);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_KaAll"), cent, pt, tofKa);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_KaAll"), cent, tofKa, tpcKa);
    }
    if (isSpecies[kPrIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_Pr"), cent, pt, tpcPr);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_Pr"), cent, pt, tofPr);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_Pr"), cent, tofPr, tpcPr);
    }
    if (isSpecies[kAntiPrIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_AntiPr"), cent, pt, tpcPr);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_AntiPr"), cent, pt, tofPr);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_AntiPr"), cent, tofPr, tpcPr);
    }
    if (isSpecies[kPrAllIdx]) {
      histos.fill(HIST("h3DnsigmaTpcVsPtAftCut_Cent_PrAll"), cent, pt, tpcPr);
      histos.fill(HIST("h3DnsigmaTofVsPtAftCut_Cent_PrAll"), cent, pt, tofPr);
      histos.fill(HIST("h3DnsigmaTpcVsTofAftCut_Cent_PrAll"), cent, tofPr, tpcPr);
    } else {
      return;
    }
  }

  // Returns: 0 = Unknown/Reject, 1 = Pion, 2 = Kaon, 3 = Proton
  template <typename T>
  int identifyTrack(const T& candidate, int cent)
  {
    if (!candidate.hasTPC())
      return 0;

    float pt = candidate.pt();
    if (pt <= cfgCutPtLower || pt >= cfgCutPtUpperPID)
      return 0; // Out of bounds

    if (!pidMeanSigmaMap)
      return 0;

    int centBin = cent + 0.5;
    auto charge = candidate.sign();
    int piIdx = (charge > 0) ? kPiPlusIdx : kPiMinusIdx;
    int kaIdx = (charge > 0) ? kKaPlusIdx : kKaMinusIdx;
    int prIdx = (charge > 0) ? kPrIdx : kAntiPrIdx;

    // TPC
    float mPiTpc = pidMeanSigmaMap->meanTPC[piIdx][centBin];
    float sPiTpc = pidMeanSigmaMap->sigmaTPC[piIdx][centBin];

    float mKaTpc = pidMeanSigmaMap->meanTPC[kaIdx][centBin];
    float sKaTpc = pidMeanSigmaMap->sigmaTPC[kaIdx][centBin];

    float mPrTpc = pidMeanSigmaMap->meanTPC[prIdx][centBin];
    float sPrTpc = pidMeanSigmaMap->sigmaTPC[prIdx][centBin];

    // TOF
    float mPiTof = pidMeanSigmaMap->meanTOF[piIdx][centBin];
    float sPiTof = pidMeanSigmaMap->sigmaTOF[piIdx][centBin];

    float mKaTof = pidMeanSigmaMap->meanTOF[kaIdx][centBin];
    float sKaTof = pidMeanSigmaMap->sigmaTOF[kaIdx][centBin];

    float mPrTof = pidMeanSigmaMap->meanTOF[prIdx][centBin];
    float sPrTof = pidMeanSigmaMap->sigmaTOF[prIdx][centBin];

    static int debugLogCounter = 0;
    if (debugLogCounter < KConstTen) {
      LOGF(info, "[PID DEBUG] CentBin: %d, Charge: %d", centBin, charge);
      LOGF(info, "   -> TPC USED | Pi (\u03bc=%.3f, \u03c3=%.3f) | Ka (\u03bc=%.3f, \u03c3=%.3f) | Pr (\u03bc=%.3f, \u03c3=%.3f)",
           mPiTpc, sPiTpc, mKaTpc, sKaTpc, mPrTpc, sPrTpc);

      if (candidate.hasTOF()) {
        LOGF(info, "   -> TOF USED | Pi (\u03bc=%.3f, \u03c3=%.3f) | Ka (\u03bc=%.3f, \u03c3=%.3f) | Pr (\u03bc=%.3f, \u03c3=%.3f)",
             mPiTof, sPiTof, mKaTof, sKaTof, mPrTof, sPrTof);
      } else {
        LOGF(info, "   -> TOF USED | Track has no TOF signal.");
      }
      debugLogCounter++;
    }
    // Fetch Raw nSigma Values
    float rawTpcPi = candidate.tpcNSigmaPi();
    float rawTpcKa = candidate.tpcNSigmaKa();
    float rawTpcPr = candidate.tpcNSigmaPr();

    float rawTofPi = 0.f, rawTofKa = 0.f, rawTofPr = 0.f;
    if (candidate.hasTOF()) {
      rawTofPi = candidate.tofNSigmaPi();
      rawTofKa = candidate.tofNSigmaKa();
      rawTofPr = candidate.tofNSigmaPr();
    }

    // --- Low PT Regime ---
    if (pt <= cfgCutPtUpperTPC) {
      // Basic TPC passing check: |Raw - Mean| < (Cut * Sigma)
      bool inTpcPi = std::abs(rawTpcPi - mPiTpc) < (cfgnSigmaCutTPC * sPiTpc);
      bool inTpcKa = std::abs(rawTpcKa - mKaTpc) < (cfgnSigmaCutTPC * sKaTpc);
      bool inTpcPr = std::abs(rawTpcPr - mPrTpc) < (cfgnSigmaCutTPC * sPrTpc);

      // Combined passing check (adds TOF if available)
      bool passPi = inTpcPi && (!candidate.hasTOF() || std::abs(rawTofPi - mPiTof) < (cfgnSigmaCutTOF * sPiTof));
      bool passKa = inTpcKa && (!candidate.hasTOF() || std::abs(rawTofKa - mKaTof) < (cfgnSigmaCutTOF * sKaTof));
      bool passPr = inTpcPr && (!candidate.hasTOF() || std::abs(rawTofPr - mPrTof) < (cfgnSigmaCutTOF * sPrTof));

      // Uniqueness check: Must pass target cut, and NOT fall into the TPC range of the others
      if (passPi && !passKa && !passPr)
        return 1;
      if (passKa && !passPi && !passPr)
        return 2;
      if (passPr && !passPi && !passKa)
        return 3;

      return 0; // Ambiguous or failed all cuts
    }

    // --- High PT Regime---
    if (candidate.hasTOF() && pt > cfgCutPtUpperTPC) {
      // Calculate 2D Normalized Distance (Elliptical distance normalized by sigma)
      float dPi = std::hypot((rawTpcPi - mPiTpc) / sPiTpc, (rawTofPi - mPiTof) / sPiTof);
      float dKa = std::hypot((rawTpcKa - mKaTpc) / sKaTpc, (rawTofKa - mKaTof) / sKaTof);
      float dPr = std::hypot((rawTpcPr - mPrTpc) / sPrTpc, (rawTofPr - mPrTof) / sPrTof);

      // Count how many particles are within the ambiguity radius
      int competitors = (dPi < cfgnSigmaOtherParticles) +
                        (dKa < cfgnSigmaOtherParticles) +
                        (dPr < cfgnSigmaOtherParticles);

      // If 1 or fewer are in the ambiguity region, pick the absolute best match
      if (competitors <= 1) {
        if (dPi <= dKa && dPi <= dPr && dPi < cfgnSigmaCutCombTPCTOF)
          return 1;
        if (dKa <= dPi && dKa <= dPr && dKa < cfgnSigmaCutCombTPCTOF)
          return 2;
        if (dPr <= dPi && dPr <= dKa && dPr < cfgnSigmaCutCombTPCTOF)
          return 3;
      }
    }
    return 0; // Unknown/Reject
  }

  template <typename T>
  bool isEventSelected(const T& col)
  {
    histos.fill(HIST("hEvtCount"), 0.5);

    if (!col.sel8())
      return false;
    histos.fill(HIST("hEvtCount"), 1.5);

    if (std::abs(col.posZ()) > cfgCutVertex)
      return false;
    histos.fill(HIST("hEvtCount"), 2.5);

    if (cfgEvSelNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return false;
    histos.fill(HIST("hEvtCount"), 3.5);

    if (cfgUseGoodITSLayerAllCut && !col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    histos.fill(HIST("hEvtCount"), 4.5);

    if (cfgIsGoodZvtxFT0VsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    histos.fill(HIST("hEvtCount"), 5.5);
    return true;
  }

  bool isPassAddPileup(int multPV, int trksize, float cent)
  {
    auto checkLimits = [](float x, float y, const std::vector<std::pair<float, float>>& limits, float xM, float xMx) {
      if (limits.empty())
        return true;
      int bin = 1 + static_cast<int>((x - xM) / (xMx - xM) * (limits.size() - 2));
      if (bin < 1 || bin >= static_cast<int>(limits.size() - 1))
        return false;
      return (y >= limits[bin].first && y <= limits[bin].second);
    };
    if (cfgApplySigPupCut) {
      if (!checkLimits(cent, trksize, state.mLimitsNchCent, state.mMinXNchCent, state.mMaxXNchCent))
        return false;
      histos.fill(HIST("hEvtCount"), 6.5);
    }

    if (cfgApplyLinPupCut) {
      if (trksize > (cfgLinPupParam0 + cfgLinPupParam1 * multPV))
        return false;
      histos.fill(HIST("hEvtCount"), 7.5);
      if (trksize < (cfgLinPupParam2 + cfgLinPupParam3 * multPV))
        return false;
      histos.fill(HIST("hEvtCount"), 8.5);
    }
    return true;
  }

  template <typename T>
  bool isTrackSelected(const T& trk)
  {
    histos.fill(HIST("hTrkCount"), 0.5);

    if (trk.sign() == 0)
      return false;
    histos.fill(HIST("hTrkCount"), 1.5);

    if (!trk.has_collision())
      return false;
    histos.fill(HIST("hTrkCount"), 2.5);

    if (!trk.isPVContributor())
      return false;
    histos.fill(HIST("hTrkCount"), 3.5);

    if (!(trk.itsNCls() > cfgITScluster))
      return false;
    histos.fill(HIST("hTrkCount"), 4.5);

    if (!(trk.tpcNClsFound() >= cfgTPCcluster))
      return false;
    histos.fill(HIST("hTrkCount"), 5.5);

    if (!(trk.tpcNClsCrossedRows() >= cfgTPCnCrossedRows))
      return false;
    histos.fill(HIST("hTrkCount"), 6.5);

    if (trk.pt() < cfgCutPtLower || trk.pt() > cfgCutPtUpper || std::abs(trk.eta()) > cfgCutEta)
      return false;
    histos.fill(HIST("hTrkCount"), 7.5);

    if (!trk.isGlobalTrack())
      return false;
    histos.fill(HIST("hTrkCount"), 8.5);

    if (cfgPtDepDCAxy) {
      // Evaluates: P0 + P1 / (pt^P2)
      float maxDcaXY = cfgDcaXyP0 + cfgDcaXyP1 / std::pow(trk.pt(), cfgDcaXyP2);
      if (std::abs(trk.dcaXY()) > maxDcaXY) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 9.5);
    } else {
      if (std::abs(trk.dcaXY()) > cfgCutTracKDcaMaxXY) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 9.5);
    }
    if (cfgPtDepDCAz) {
      // Evaluates: P0 + P1 / (pt^P2)
      float maxDcaZ = cfgDcaZP0 + cfgDcaZP1 / std::pow(trk.pt(), cfgDcaZP2);
      if (std::abs(trk.dcaZ()) > maxDcaZ) {
        return false; // Reject track if DCA is too large
      }
      histos.fill(HIST("hTrkCount"), 10.5);
    } else {
      if (std::abs(trk.dcaZ()) > cfgCutTracKDcaMaxZ) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 10.5);
    }
    return true;
  }

  template <typename T>
  bool isParticleSelected(const T& particle)
  {
    auto* pd = pdg->GetParticle(particle.pdgCode());
    if (!pd)
      return false;
    if (std::abs(pd->Charge()) == 0)
      return false;
    if (particle.pt() < cfgCutPtLower || particle.pt() > cfgCutPtUpper || std::abs(particle.eta()) > cfgCutEta)
      return false;
    if (std::abs(particle.vz()) > cfgCutVertex)
      return false;
    return true;
  }

  float getCentrality(const auto& col) const
  {
    if (cfgCentralityChoice.value == kCentFT0C)
      return col.centFT0C();
    if (cfgCentralityChoice.value == kCentFT0M)
      return col.centFT0M();
    if (cfgCentralityChoice.value == kCentFDDM)
      return col.centFDDM();
    if (cfgCentralityChoice.value == kCentFV0A)
      return col.centFV0A();
    return KinvalidCentrality;
  }

  float getEfficiency(float mult, float pt, float eta, PIDIdx pidType, int effidx, bool cfgEff) const
  {
    if (!cfgEff) {
      return (effidx == 0) ? 1.0f : 0.0f;
    }
    TH3F* h = (effidx == 0) ? state.hEff[pidType] : state.hFake[pidType];

    if (!h)
      return -1;

    int ibx = h->GetXaxis()->FindBin(mult);
    int iby = h->GetYaxis()->FindBin(pt);
    int ibz = h->GetZaxis()->FindBin(eta);
    float val = h->GetBinContent(ibx, iby, ibz);

    if (effidx == 0)
      return (val > 0.f) ? val : 1.0f;
    return val;
  }

  float getFlatteningWeight(float vz, float chg, float pt, float eta, float phi, PIDIdx pidType, bool cfgflat) const
  {
    if (!cfgflat)
      return 1.0;
    THnSparseF* h = state.hFlatWeight[pidType];

    if (!h)
      return 0.0;
    int bins[5];
    bins[0] = h->GetAxis(0)->FindBin(vz);
    bins[1] = h->GetAxis(1)->FindBin(chg);
    bins[2] = h->GetAxis(2)->FindBin(pt);
    bins[3] = h->GetAxis(3)->FindBin(eta);
    bins[4] = h->GetAxis(4)->FindBin(phi);
    float val = h->GetBinContent(bins);

    return val;
  }

  std::vector<o2::detectors::AlignParam>* offsetFT0 = nullptr;
  uint64_t mLastTimestamp = 0;
  double getEtaFT0(uint64_t globalChno, int i)
  {
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }
    auto chPos = ft0Det.getChannelCenter(globalChno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    // i == 0 is FT0A (A-side, positive Z)
    // i == 1 is FT0C (C-side, negative Z)
    if (i == 1) {
      z = -std::abs(z);
    } else if (i == 0) {
      z = std::abs(z);
    }
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  void loadAlignParam(uint64_t timestamp)
  {
    if (timestamp == mLastTimestamp && offsetFT0 != nullptr) {
      return;
    }
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    if (!offsetFT0) {
      LOGF(fatal, "Failed to load valid FT0 alignment from CCDB!");
      return;
    }
    mLastTimestamp = timestamp;
    LOGF(info, "Successfully loaded new alignment parameters for timestamp %llu", timestamp);
    LOGF(info, "Offset for FT0A: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[0].getX(), (*offsetFT0)[0].getY(), (*offsetFT0)[0].getZ());
    LOGF(info, "Offset for FT0C: x = %.3f y = %.3f z = %.3f\n", (*offsetFT0)[1].getX(), (*offsetFT0)[1].getY(), (*offsetFT0)[1].getZ());
  }

  template <int KIntM, int KIntK>
  std::pair<float, float> calculateMeanAndC2FromSums(const double sumpmwk[KIntM][KIntK], const double sumwk[KIntK], float referenceMeanPt) const
  {
    if (sumwk[1] == 0.) {
      return {0.f, 0.f};
    }

    double tau1 = sumwk[2] / (sumwk[1] * sumwk[1]);
    double denom2 = 1. - tau1;

    if (std::abs(denom2) < KFloatEpsilon) {
      double pmk11safe = sumpmwk[1][1] / sumwk[1];
      return {static_cast<float>(pmk11safe), 0.f};
    }

    double pmk11 = sumpmwk[1][1] / sumwk[1];

    double pmk12 = 0.f;
    if (sumwk[2] != 0.f) {
      pmk12 = sumpmwk[1][2] / sumwk[2];
    }

    double pmk22 = 0.f;
    if (sumwk[2] != 0.f) {
      pmk22 = sumpmwk[2][2] / sumwk[2];
    }

    float calculatedMeanPt = pmk11;

    double p1kBar1 = pmk11 - referenceMeanPt;
    double p2kBar2 = pmk22 - 2.0f * pmk12 * referenceMeanPt + referenceMeanPt * referenceMeanPt;

    double p1kBar1sq = p1kBar1 * p1kBar1;
    double numerator2 = p1kBar1sq - (tau1 * p2kBar2);

    float twopcorr = numerator2 / denom2;
    return {calculatedMeanPt, twopcorr};
  }

  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                      aod::FT0sCorrected,
                                      aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
                                      aod::CentNTPVs>;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZCut;
  using AodCollisionsSel = soa::Filtered<GeneralCollisions>;

  using UnfilteredTracks = soa::Join<
    aod::Tracks,
    aod::TracksExtra,
    aod::TrackSelection,
    aod::TracksDCA,
    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl,
    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl>;
  Filter trackFilter = aod::track::pt > KPtMin&&
                                          aod::track::pt < KPtMax&&
                       requireGlobalTrackInFilter();
  using AodTracksSel = soa::Filtered<UnfilteredTracks>;
  using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
  using FilteredTCs = soa::Filtered<TCs>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  using MyRun3MCCollisions = soa::Join<
    aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
    aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
    aod::CentNGlobals, aod::McCollisionLabels>;

  PresliceUnsorted<MyRun3MCCollisions> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  void declareCommonQA()
  {
    histos.add("hVtxZ_after_sel", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hVtxZ", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis1Per}});
    histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis}, {nChAxis}});
    histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis}, {centAxis1Per}});

    histos.add("Hist2D_globalTracks_cent", "cent (%);N_{global}", kTH2F, {{centAxis1Per}, {nChAxis}});
    histos.add("Hist2D_PVTracks_cent", "cent (%);N_{PV}", kTH2F, {{centAxis1Per}, {nChAxis}});

    histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsPt, KPMin, KPMax}});
    histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("hEta", ";#eta", kTH1F, {{KNbinsEtaFine, KEtaMin, KEtaMax}});
    histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});

    histos.add("hEvtCount", "Number of Event;; Count", kTH1F, {{9, 0, 9}});
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(1, "all Events");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(3, "after VertexZ Cut");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(4, "after kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(5, "after kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(6, "after kIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(7, "after PVTracksCent Pileup Cut");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(8, "after Linear Pileup Cut (Up)");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(9, "after Linear Pileup Cut (Lw)");

    histos.add("hTrkCount", "Number of Tracks;; Count", kTH1F, {{11, 0, 11}});
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(1, "all Tracks");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(2, "after sign!=0");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(3, "after has_collision");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(4, "after isPVContributor");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(5, "after itsNCls");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(6, "after tpcNClsFound");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(7, "after tpcNClsCrossedRows");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(8, "after pT,#eta selections");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(9, "after isGlobalTrack");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(10, "after dcaXY");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(11, "after dcaZ");
  }

  void declareMCCommonHists()
  {
    for (const auto& suf : pidSuffix) {
      histos.add("h3_AllPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoMatchedToPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_AllReco" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoUnMatchedToPrimary_Secondary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoUnMatchedToPrimary_Fake" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoMatchedToPrimary_MisID" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    }

    histos.add("ptResolution", ";p_{T}^{MC};(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}", kTH2F, {{KNbinsPtRes, KPtMin, KPtMax}, {100, -0.2, 0.2}});
    histos.add("etaResolution", ";#eta^{MC};#eta^{reco}-#eta^{MC}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {100, -0.02, 0.02}});
    histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}});
    histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
    histos.add("vzResolution", ";Vz^{MC};(Vz^{reco}-Vz^{MC})/Vz^{MC}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {100, -0.1, 0.1}});
  }

  void declarenSigHists()
  {
    for (const auto& suf : pidSuffix) {
      histos.add("h3DnsigmaTpcVsPtBefCut_Cent" + suf, "TPC nSigma vs pT Before Cut;cent [%]; p_{T} (GeV/c);n#sigma_{TPC}", kTH3F, {{centAxis1Per}, {KNbinsPtRes, KPtMin, KPtMax}, {200, -10.f, 10.f}});
      histos.add("h3DnsigmaTofVsPtBefCut_Cent" + suf, "TOF nSigma vs pT Before Cut;cent [%]; p_{T} (GeV/c);n#sigma_{TOF}", kTH3F, {{centAxis1Per}, {KNbinsPtRes, KPtMin, KPtMax}, {200, -10.f, 10.f}});
      histos.add("h3DnsigmaTpcVsTofBefCut_Cent" + suf, "TPC vs TOF nSigma Before Cut;cent [%]; n#sigma_{TOF};n#sigma_{TPC}", kTH3F, {{centAxis1Per}, {200, -10.f, 10.f}, {200, -10.f, 10.f}});
      histos.add("h3DnsigmaTpcVsPtAftCut_Cent" + suf, "TPC nSigma vs pT After Cut;cent [%],; p_{T} (GeV/c);n#sigma_{TPC}", kTH3F, {{centAxis1Per}, {KNbinsPtRes, KPtMin, KPtMax}, {200, -10.f, 10.f}});
      histos.add("h3DnsigmaTofVsPtAftCut_Cent" + suf, "TOF nSigma vs pT After Cut;cent [%],; p_{T} (GeV/c);n#sigma_{TOF}", kTH3F, {{centAxis1Per}, {KNbinsPtRes, KPtMin, KPtMax}, {200, -10.f, 10.f}});
      histos.add("h3DnsigmaTpcVsTofAftCut_Cent" + suf, "TPC vs TOF nSigma After Cut;cent [%],; n#sigma_{TOF};n#sigma_{TPC}", kTH3F, {{centAxis1Per}, {200, -10.f, 10.f}, {200, -10.f, 10.f}});
    }
  }

  void declareMCGetFlatHists()
  {
    for (const auto& suf : pidSuffix) {
      histos.add("MCGen/hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("MCGen/hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("MCGen/hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
    }
  }

  void declareMCMeanHists()
  {
    histos.add("Eff_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("Eff_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("Eff_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Eff_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("Fake_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("Fake_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("Fake_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Fake_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("wgt_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("wgt_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("wgt_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("wgt_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("pmeanFT0Amultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
    histos.add("pmeanFT0Cmultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0C_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;id;#eta", kTProfile3D, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;id;#eta", kTH3F, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});

    histos.add<TProfile2D>("MCGen/Prof_Cent_Nsp_Nchrec", ";cent;isp", kTProfile2D, {{centAxis1Per}, {spBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Mult_Nsp_Nchrec", ";mult;isp", kTProfile2D, {{nChAxis}, {spBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Cent_Nsp_MeanpT", ";cent;isp", kTProfile2D, {{centAxis1Per}, {spBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Mult_Nsp_MeanpT", ";mult;isp", kTProfile2D, {{nChAxis}, {spBinAxis}});

    histos.add<TProfile3D>("pmeanTru_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanReco_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanRecoEffcorr_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});

    histos.add<TProfile3D>("pmeanMultTru_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanMultReco_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanMultRecoEffcorr_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});

    for (const auto& suf : pidSuffix) {
      histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});

      histos.add<TProfile3D>("Prof2D_MeanpTSub_Tru" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
      histos.add<TProfile3D>("Prof2D_MeanpTSub_Reco" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
      histos.add<TProfile3D>("Prof2D_MeanpTSub_RecoEffCorr" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
    }
  }

  void declareMCFlucHists()
  {
    histos.add<TProfile3D>("MCGen/Prof_Cent_NEta_Nsp_Nchrec", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Mult_NEta_Nsp_Nchrec", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Cent_NEta_Nsp_MeanpT", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Mult_NEta_Nsp_MeanpT", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});

    histos.add<TProfile3D>("MCGen/Prof_MeanpT_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_C2_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_C2Sub_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Cov_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0A_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0C_Cent_etabin_spbin", ";cent;eta;isp", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {spBinAxis}});

    histos.add<TProfile3D>("MCGen/Prof_MeanpT_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_C2_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_C2Sub_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Cov_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0A_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0C_Mult_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {etaBinAxis}, {spBinAxis}});

    for (const auto& suf : pidSuffix) {
      histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});

      histos.add<TProfile3D>(Form("MCGen/Prof_C2Sub2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>(Form("MCGen/Prof_GapSum2D%s", suf.c_str()), ";cent;gap;sum", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
      histos.add<TProfile3D>(Form("MCGen/Prof_Cov2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>(Form("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>(Form("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    }
  }

  void declareDataGetFlatHists()
  {
    for (const auto& suf : pidSuffix) {
      histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
    }
    histos.add("hnTrkPVZDC", ";N_{PV};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 3000}});
    histos.add("hNchZDC", ";N_{trk};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 30000}});
  }

  void declareDataMeanHists()
  {
    histos.add("pmeanFT0Amultpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0Cmultpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0C_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});

    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;channel id; #eta;amplitude", kTProfile3D, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;channel id; #eta", kTH3F, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});

    histos.add<TProfile2D>("Prof_Cent_Nsp_Nchrec", ";cent;Species;#LT N_{PV}#GT", kTProfile2D, {{centAxis1Per}, {spBinAxis}});
    histos.add<TProfile2D>("Prof_Mult_Nsp_Nchrec", ";N_{PV};Species;#LT N_{PV}#GT", kTProfile2D, {{nChAxis}, {spBinAxis}});
    histos.add<TProfile2D>("Prof_Cent_Nsp_MeanpT", ";cent;Species;#LT p_{T}#GT", kTProfile2D, {{centAxis1Per}, {spBinAxis}});
    histos.add<TProfile2D>("Prof_Mult_Nsp_MeanpT", ";N_{PV};Species;#LT p_{T}#GT", kTProfile2D, {{nChAxis}, {spBinAxis}});

    histos.add<TProfile3D>("pmean_nch_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanMult_nch_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("pmean_cent_etabin_spbin", ";Centrality (%) ;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("pmeanMult_cent_etabin_spbin", ";Centrality (%) ;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});

    for (const auto& suf : pidSuffix) {
      histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});

      histos.add<TProfile3D>("Prof2D_MeanpTSub" + suf, ";cent;#eta_{A} bin;#eta_{C} bin", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {{etaBinAxis}}});
    }
  }

  void declareDataFlucHists()
  {
    histos.add<TProfile3D>("Prof_MeanpT_Cent_etabin_spbin", ";cent;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_MeanpT_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_C2_Cent_etabin_spbin", ";cent;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_C2_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});

    histos.add<TProfile3D>("Prof_C2Sub_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_C2Sub_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_Cov_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_Cov_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});

    histos.add<TProfile3D>("Prof_CovFT0A_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_CovFT0A_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_CovFT0C_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {{etaBinAxis}}, {spBinAxis}});
    histos.add<TProfile3D>("Prof_CovFT0C_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {{etaBinAxis}}, {spBinAxis}});

    for (const auto& suf : pidSuffix) {
      histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      histos.add<TProfile3D>("Prof_C2Sub2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>("Prof_GapSum2D" + suf, ";cent;#Delta#eta (Gap);#Sigma#eta/2 (Sum)", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
      histos.add<TProfile3D>("Prof_Cov2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A} bin;#eta_{C} bin", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>("Prof_CovFT0A2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      histos.add<TProfile3D>("Prof_CovFT0C2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    }
  }

  THnSparseF* buildWeightMapFromRaw(THnSparseF* hRaw, const char* mapName)
  {
    if (!hRaw) {
      LOGF(error, "Raw eta-phi map for '%s' is null; no flattening will be applied.", mapName);
      return nullptr;
    }
    auto hWMap = reinterpret_cast<THnSparseF*>(hRaw->Clone(mapName));
    hWMap->SetTitle(Form("Flattening Weight Map %s (w_{#phi} = <N_{#phi}> / N_{#phi})", mapName));
    hWMap->Reset();
    auto axV = hRaw->GetAxis(0);   // Vz
    auto axChg = hRaw->GetAxis(1); // Charge
    auto axPt = hRaw->GetAxis(2);  // Pt
    auto axE = hRaw->GetAxis(3);   // Eta
    auto axP = hRaw->GetAxis(4);   // Phi

    int bins[5];
    for (int iv = 1; iv <= axV->GetNbins(); ++iv) {
      bins[0] = iv;
      for (int ichg = 1; ichg <= axChg->GetNbins(); ++ichg) {
        bins[1] = ichg;
        for (int ipt = 1; ipt <= axPt->GetNbins(); ++ipt) {
          bins[2] = ipt;
          for (int ie = 1; ie <= axE->GetNbins(); ++ie) {
            bins[3] = ie;
            double sum = 0.0;
            int nphi = axP->GetNbins();
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = ip;
              sum += hRaw->GetBinContent(bins);
            }
            const double avg = (nphi > 0 ? sum / nphi : 0.0);
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = ip;
              const double raw = hRaw->GetBinContent(bins);
              const double w = (avg > 0.0 && raw > 0.0) ? (avg / raw) : 1.0;
              hWMap->SetBinContent(bins, w);
            }
          }
        }
      }
    }

    LOGF(info, "Flattening weight map '%s' built.", mapName);
    return hWMap;
  }

  inline void loadTProfile3D(TDirectory* dir, const char* name, TProfile3D*& target)
  {
    if (!dir) {
      LOGF(error, "loadTProfile3D: directory is null for object %s", name);
      return;
    }
    auto* obj = dir->Get(name);
    if (!obj) {
      LOGF(error, "loadTProfile3D: object '%s' not found in directory %s", name, dir->GetName());
      return;
    }
    auto* prof = dynamic_cast<TProfile3D*>(obj);
    if (!prof) {
      LOGF(error, "loadTProfile3D: object '%s' is not a TProfile3D (it is %s)", name, obj->ClassName());
      return;
    }
    target = reinterpret_cast<TProfile3D*>(prof->Clone(Form("%s_clone", name)));
    target->SetDirectory(nullptr);
    LOGF(info, "Loaded TProfile3D '%s' with entries = %.0f", name, target->GetEntries());
  }

  void init(InitContext&)
  {
    if (cfgSys == kPbPb) {
      nChAxis = {cfgNchPbMax / 2, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchPbMax / 4, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else {
      nChAxis = {cfgNchOMax, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchOMax, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    }

    ccdb->setURL(cfgCCDBurl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    loadAlignParam(now);
    ft0Det.calculateChannelCenter();

    std::string sysDir = "";
    switch (cfgSys) {
      case kPbPb:
        sysDir = "PbPbTest";
        break;
      case kNeNe:
        sysDir = "NeNeTest";
        break;
      case kOO:
        sysDir = "OOTest";
        break;
      case kpp:
        sysDir = "ppTest";
        break;
      default:
        LOGF(fatal, "Invalid cfgSys value: %d", cfgSys.value);
    }
    std::string pathNsig = cfgCCDBUserPath.value + "/" + sysDir + "/Job0_nSigMaps";
    std::string pathEff = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_EffMaps";
    std::string pathMCFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_MCFlatMaps";
    std::string pathMCMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_MCMean";

    std::string pathDataNsig = cfgCCDBUserPath.value + "/" + sysDir + "/Job0_DatanSigMaps";
    std::string pathDataFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_DataFlatMaps";
    std::string pathDataMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_DataMean";

    declareCommonQA();
    if (cfgRunMCGetNSig || cfgRunGetEff || cfgRunDataGetNSig || cfgRunGetDataFlat) {
      declarenSigHists();
    }

    if (cfgRunMCMean || cfgRunMCFluc || cfgRunGetEff) {
      declareMCCommonHists();
    }
    if (cfgRunGetMCFlat) {
      declareMCGetFlatHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunMCMean) {
      declareMCMeanHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunMCFluc) {
      declareMCFlucHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunDataGetNSig) {
      declareDataGetFlatHists();
    }
    if (cfgRunGetDataFlat) {
      declareDataGetFlatHists();
    }
    if (cfgRunDataMean) {
      declareDataMeanHists();
    }
    if (cfgRunDataFluc) {
      declareDataFlucHists();
    }

    if (!cfgRunGetEff && !cfgRunMCGetNSig && (cfgEff)) {
      TList* lst = ccdb->getForTimeStamp<TList>(pathEff, now);

      if (!lst) {
        LOGF(fatal, "Efficiency maps required but CCDB list is null at %s!", pathEff.c_str());
      }

      LOGF(info, "Loading Eff/Fake maps from TList for all species...");

      auto loadEffFakeForPID = [&](PIDIdx pidType) {
        std::string suffix = pidSuffix[pidType];
        std::string hEffNumName = "h3_RecoMatchedToPrimary" + suffix;
        std::string hEffDenName = "h3_AllPrimary" + suffix;
        std::string hFakeNumSecName = "h3_RecoUnMatchedToPrimary_Secondary" + suffix;
        std::string hFakeNumFakName = "h3_RecoUnMatchedToPrimary_Fake" + suffix;
        std::string hFakeNumFakName2 = "h3_RecoMatchedToPrimary_MisID" + suffix;
        std::string hFakeDenName = "h3_AllReco" + suffix;

        auto* hNum = reinterpret_cast<TH3F*>(lst->FindObject(hEffNumName.c_str()));
        auto* hDen = reinterpret_cast<TH3F*>(lst->FindObject(hEffDenName.c_str()));

        if (hNum && hDen) {
          state.hEff[pidType] = reinterpret_cast<TH3F*>(hNum->Clone(Form("hEff%s", suffix.c_str())));
          state.hEff[pidType]->SetDirectory(nullptr);
          state.hEff[pidType]->Divide(hDen);
        } else {
          LOGF(error, "Missing CCDB objects for efficiency. Checked: %s, %s", hEffNumName.c_str(), hEffDenName.c_str());
        }

        auto* hNumS = reinterpret_cast<TH3F*>(lst->FindObject(hFakeNumSecName.c_str()));
        auto* hNumF = reinterpret_cast<TH3F*>(lst->FindObject(hFakeNumFakName.c_str()));
        auto* hNumF2 = reinterpret_cast<TH3F*>(lst->FindObject(hFakeNumFakName2.c_str()));
        auto* hDenF = reinterpret_cast<TH3F*>(lst->FindObject(hFakeDenName.c_str()));

        if (hNumS && hNumF && hDenF) {
          state.hFake[pidType] = reinterpret_cast<TH3F*>(hNumS->Clone(Form("hFake%s", suffix.c_str())));
          state.hFake[pidType]->Add(hNumF);
          if (pidType != kInclusiveIdx && hNumF2) {
            state.hFake[pidType]->Add(hNumF2);
          }
          state.hFake[pidType]->SetDirectory(nullptr);
          state.hFake[pidType]->Divide(hDenF);
        } else {
          LOGF(error, "Missing CCDB object(s) for fakes for %s in list.", suffix.c_str());
        }
      };

      for (int i = 0; i < KNsp; ++i) {
        loadEffFakeForPID(static_cast<PIDIdx>(i));
      }
    }

    bool requiresMCMap = (cfgRunGetEff || cfgRunGetMCFlat || cfgRunMCMean || cfgRunMCFluc);
    bool requiresDataMap = (cfgRunGetDataFlat || cfgRunDataMean || cfgRunDataFluc);

    if (requiresMCMap || requiresDataMap) {
      std::string currentPath = requiresMCMap ? pathNsig : pathDataNsig;
      TList* pidList = ccdb->getForTimeStamp<TList>(currentPath, now);

      if (!pidList) {
        LOGF(warn, "nSigma maps required but CCDB list is null at %s! Using raw values.", currentPath.c_str());
      } else {
        if (!pidMeanSigmaMap) {
          pidMeanSigmaMap = new PIDMeanSigmaMap();
        }
        LOGF(info, "Performing 2D Gaussian fits on PID maps from CCDB...");
        auto loadPIDMeans = [&](PIDIdx pidType) {
          std::string suffix = pidSuffix[pidType];
          std::string hName = "h3DnsigmaTpcVsTofBefCut_Cent" + suffix;
          auto* h3 = reinterpret_cast<TH3F*>(pidList->FindObject(hName.c_str()));
          if (!h3) {
            LOGF(warn, "  [!] PID Hist %s not found in CCDB list.", hName.c_str());
            return;
          }
          int nCentBins = std::min(h3->GetXaxis()->GetNbins(), PIDMeanSigmaMap::MaxCentBins - 1);
          LOGF(info, "  -> Species: %s (Bins: %d)", hName.c_str(), nCentBins);
          for (int iCent = 1; iCent <= nCentBins; ++iCent) {
            h3->GetXaxis()->SetRange(iCent, iCent);
            // Projecting: Z(TPC) vs Y(TOF). Result: X_axis=TOF, Y_axis=TPC
            std::unique_ptr<TH2D> h2(reinterpret_cast<TH2D*>(h3->Project3D("zy")));
            if (h2) {
              int binX, binY, binZ;
              h2->GetMaximumBin(binX, binY, binZ);
              double guessMeanTOF = h2->GetXaxis()->GetBinCenter(binX);
              double guessMeanTPC = h2->GetYaxis()->GetBinCenter(binY);
              TF2 f2("f2", "[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])", -3, 3, -3, 3);
              f2.SetParameters(h2->GetMaximum(), guessMeanTOF, 1.0, guessMeanTPC, 1.0);
              h2->Fit(&f2, "QRN"); // Q=Quiet, R=Range, N=NoDraw
              pidMeanSigmaMap->meanTOF[pidType][iCent - 1] = f2.GetParameter(1);
              pidMeanSigmaMap->meanTPC[pidType][iCent - 1] = f2.GetParameter(3);
              pidMeanSigmaMap->sigmaTOF[pidType][iCent - 1] = std::abs(f2.GetParameter(2));
              pidMeanSigmaMap->sigmaTPC[pidType][iCent - 1] = std::abs(f2.GetParameter(4));
              if (iCent % KConstTen == 0) {
                LOGF(info, "  Derived: For Species: %s (Bins: %d), Mean TOF = %.3f, Mean TPC = %.3f, Sigma TOF = %.3f, Sigma TPC = %.3f", hName.c_str(), iCent - 1, f2.GetParameter(1), f2.GetParameter(3), f2.GetParameter(2), f2.GetParameter(4));
              }
            }
          }
        };
        for (int i = 1; i < KNsp; ++i) {
          loadPIDMeans(static_cast<PIDIdx>(i));
        }

        auto loadLimits = [&](const char* name, std::vector<std::pair<float, float>>& limits, float& xMin, float& xMax) {
          auto* h2 = reinterpret_cast<TH2*>(pidList->FindObject(name));
          if (!h2)
            return;

          std::unique_ptr<TProfile> prof(h2->ProfileX("ptmp", 1, -1, "S"));

          int nBins = prof->GetNbinsX();
          xMin = prof->GetXaxis()->GetXmin();
          xMax = prof->GetXaxis()->GetXmax();

          limits.assign(nBins + 2, {-99999.f, 999999.f});

          for (int i = 1; i <= nBins; ++i) {
            float mean = prof->GetBinContent(i);
            float rms = prof->GetBinError(i);

            limits[i] = {mean - cfgPupnSig * rms, mean + cfgPupnSig * rms};
          }
        };
        loadLimits("Hist2D_globalTracks_cent", state.mLimitsNchCent, state.mMinXNchCent, state.mMaxXNchCent);
      }
    }
    if (!cfgRunGetEff && (cfgFlat)) {
      if (cfgRunDataMean || cfgRunDataFluc) {
        LOGF(info, "Data Run: Loading flattening maps from %s", pathDataFlat.c_str());
        TList* lstDataFlat = ccdb->getForTimeStamp<TList>(pathDataFlat, now);

        if (lstDataFlat) {
          for (int i = 0; i < KNsp; ++i) {
            std::string suffix = pidSuffix[i];
            std::string hName;

            if (cfgEff && cfgFlat) {
              hName = "hEtaPhiRecoWtd" + suffix;
            } else if (cfgEff) {
              hName = "hEtaPhiRecoEffWtd" + suffix;
            } else {
              hName = "hEtaPhiReco" + suffix;
            }
            auto* hRaw = reinterpret_cast<THnSparseF*>(lstDataFlat->FindObject(hName.c_str()));

            if (hRaw) {
              state.hFlatWeight[i] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
            } else {
              LOGF(error, "Data flattening map '%s' not found.", hName.c_str());
            }
          }
        } else {
          LOGF(error, "Could not retrieve Data Flattening TList from: %s", pathDataFlat.c_str());
        }
      }

      if (cfgRunMCMean || cfgRunMCFluc) {
        LOGF(info, "MC Run: Loading flattening maps from %s", pathMCFlat.c_str());
        TList* lstMCFlat = ccdb->getForTimeStamp<TList>(pathMCFlat, now);

        if (lstMCFlat) {
          auto loadFlatForPID = [&](PIDIdx pidType) {
            std::string suffix = pidSuffix[pidType];
            std::string hFlatSrcName;
            if (cfgEff && cfgFlat) {
              hFlatSrcName = "MCReco/hEtaPhiRecoWtd" + suffix;
            } else if (cfgEff) {
              hFlatSrcName = "MCReco/hEtaPhiRecoEffWtd" + suffix;
            } else {
              hFlatSrcName = "MCReco/hEtaPhiReco" + suffix;
            }

            auto* hRaw = reinterpret_cast<THnSparseF*>(lstMCFlat->FindObject(hFlatSrcName.c_str()));

            if (hRaw) {
              state.hFlatWeight[pidType] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
            } else {
              LOGF(warning, "MC flattening source '%s' not found in list.", hFlatSrcName.c_str());
            }
          };

          for (int i = 0; i < KNsp; ++i) {
            loadFlatForPID(static_cast<PIDIdx>(i));
          }
        } else {
          LOGF(error, "Could not retrieve MC Flattening TList from: %s", pathMCFlat.c_str());
        }
      }
    }

    auto loadTProfile3DFromList = [&](TList* sourceList, const char* objName, TProfile3D*& target) {
      if (!sourceList)
        return;
      auto* tp = reinterpret_cast<TProfile3D*>(sourceList->FindObject(objName));
      if (tp) {
        target = reinterpret_cast<TProfile3D*>(tp->Clone());
        target->SetDirectory(nullptr);
        LOGF(info, "Loaded %s from list", objName);
      } else {
        LOGF(error, "Histogram %s missing in CCDB TList", objName);
      }
    };

    auto loadTProfileFromList = [&](TList* sourceList, const char* objName, TProfile*& target) {
      if (!sourceList)
        return;
      auto* tp = reinterpret_cast<TProfile*>(sourceList->FindObject(objName));
      if (tp) {
        target = reinterpret_cast<TProfile*>(tp->Clone());
        target->SetDirectory(nullptr);
        LOGF(info, "Loaded %s from list", objName);
      } else {
        LOGF(error, "Histogram %s missing in CCDB TList", objName);
      }
    };

    if (cfgRunMCFluc) {
      LOGF(info, "Loading MC Mean profiles from CCDB path: %s", pathMCMean.c_str());
      TList* lstMCMean = ccdb->getForTimeStamp<TList>(pathMCMean, now);

      if (lstMCMean) {
        loadTProfileFromList(lstMCMean, "pmeanFT0Amultpv", state.pmeanFT0AmultpvStep2);
        loadTProfileFromList(lstMCMean, "pmeanFT0Cmultpv", state.pmeanFT0CmultpvStep2);

        loadTProfile3DFromList(lstMCMean, "pmeanTru_nch_etabin_spbin", state.pmeanTruNchEtabinSpbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanReco_nch_etabin_spbin", state.pmeanRecoNchEtabinSpbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanRecoEffcorr_nch_etabin_spbin", state.pmeanRecoEffcorrNchEtabinSpbinStep2);

        loadTProfile3DFromList(lstMCMean, "pmeanMultTru_nch_etabin_spbin", state.pmeanMultTruNchEtabinSpbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanMultReco_nch_etabin_spbin", state.pmeanMultRecoNchEtabinSpbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanMultRecoEffcorr_nch_etabin_spbin", state.pmeanMultRecoEffcorrNchEtabinSpbinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for MC Mean from: %s", pathMCMean.c_str());
      }
    }

    if (cfgRunDataFluc) {
      LOGF(info, "Loading Data Mean profiles from CCDB path: %s", pathDataMean.c_str());
      TList* lstDataMean = ccdb->getForTimeStamp<TList>(pathDataMean, now);

      if (lstDataMean) {
        loadTProfileFromList(lstDataMean, "pmeanFT0Amultpv", state.pmeanFT0AmultpvStep2);
        loadTProfileFromList(lstDataMean, "pmeanFT0Cmultpv", state.pmeanFT0CmultpvStep2);

        loadTProfile3DFromList(lstDataMean, "pmean_nch_etabin_spbin", state.pmeanNchEtabinSpbinStep2);
        loadTProfile3DFromList(lstDataMean, "pmeanMult_nch_etabin_spbin", state.pmeanMultNchEtabinSpbinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for Data Mean from: %s", pathDataMean.c_str());
      }
    }
    LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
  }

  void processMCGetMeanNsig(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks)
  {
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision))
      return;
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax)
      return;
    float multPV = mcCollision.multNTracksPV();

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index())
        continue;

      if (!isTrackSelected(track))
        continue;
      fillNSigmaBefCut(track, cent);
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCGetMeanNsig, "process MC to calculate Mean values of nSig Plots", cfgRunMCGetNSig);

  void processGetEffHists(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision))
      return;
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax)
      return;
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent))
      return;
    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId())
        continue;

      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
        continue;

      const int pdg = particle.pdgCode();
      const int absPdg = std::abs(pdg);
      float pt = particle.pt(), eta = particle.eta();

      bool isSpecies[KNsp] = {
        true,              // kInclusiveIdx
        pdg == -KPiPlus,   // kPiMinusIdx
        pdg == KPiPlus,    // kPiPlusIdx
        absPdg == KPiPlus, // kPiAllIdx
        pdg == -KKPlus,    // kKaMinusIdx
        pdg == KKPlus,     // kKaPlusIdx
        absPdg == KKPlus,  // kKaAllIdx
        pdg == -KProton,   // kAntiPrIdx
        pdg == KProton,    // kPrIdx
        absPdg == KProton  // kPrAllIdx
      };

      histos.fill(HIST("h3_AllPrimary"), multPV, pt, eta);
      if (isSpecies[kPiMinusIdx])
        histos.fill(HIST("h3_AllPrimary_PiMinus"), multPV, pt, eta);
      else if (isSpecies[kPiPlusIdx])
        histos.fill(HIST("h3_AllPrimary_PiPlus"), multPV, pt, eta);
      if (isSpecies[kPiAllIdx])
        histos.fill(HIST("h3_AllPrimary_PiAll"), multPV, pt, eta);

      if (isSpecies[kKaMinusIdx])
        histos.fill(HIST("h3_AllPrimary_KaMinus"), multPV, pt, eta);
      else if (isSpecies[kKaPlusIdx])
        histos.fill(HIST("h3_AllPrimary_KaPlus"), multPV, pt, eta);
      if (isSpecies[kKaAllIdx])
        histos.fill(HIST("h3_AllPrimary_KaAll"), multPV, pt, eta);

      if (isSpecies[kAntiPrIdx])
        histos.fill(HIST("h3_AllPrimary_AntiPr"), multPV, pt, eta);
      else if (isSpecies[kPrIdx])
        histos.fill(HIST("h3_AllPrimary_Pr"), multPV, pt, eta);
      if (isSpecies[kPrAllIdx])
        histos.fill(HIST("h3_AllPrimary_PrAll"), multPV, pt, eta);
    }

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index())
        continue;

      if (!isTrackSelected(track))
        continue;

      float pt = track.pt(), eta = track.eta();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      auto sign = track.sign();
      fillNSigmaBefCut(track, cent);

      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      fillNSigmaAftCut(track, cent, isSpecies);

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;

        if (isp == kInclusiveIdx) {
          histos.fill(HIST("h3_AllReco"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              histos.fill(HIST("ptResolution"), mcP.pt(), (pt - mcP.pt()) / mcP.pt());
              histos.fill(HIST("etaResolution"), mcP.eta(), eta - mcP.eta());
              histos.fill(HIST("etaTruthReco"), mcP.eta(), eta);
              histos.fill(HIST("vzResolution"), mcP.vz(), (vz - mcP.vz()) / mcP.vz());
              histos.fill(HIST("TruthTracKVz"), mcP.vz(), vz);
              histos.fill(HIST("h3_RecoMatchedToPrimary"), multPV, mcP.pt(), mcP.eta());
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), multPV, pt, eta);
          }
        } else if (isp == kPiMinusIdx) {
          histos.fill(HIST("h3_AllReco_PiMinus"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == -KPiPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PiMinus"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_PiMinus"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiMinus"), multPV, pt, eta);
            }
          } else { // No MC
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiMinus"), multPV, pt, eta);
          }
        } else if (isp == kPiPlusIdx) {
          histos.fill(HIST("h3_AllReco_PiPlus"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == KPiPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PiPlus"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_PiPlus"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiPlus"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiPlus"), multPV, pt, eta);
          }
        } else if (isp == kPiAllIdx) {
          histos.fill(HIST("h3_AllReco_PiAll"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (std::abs(mcP.pdgCode()) == KPiPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PiAll"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_PiAll"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiAll"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiAll"), multPV, pt, eta);
          }
        } else if (isp == kKaMinusIdx) {
          histos.fill(HIST("h3_AllReco_KaMinus"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == -KKPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_KaMinus"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_KaMinus"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaMinus"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaMinus"), multPV, pt, eta);
          }
        } else if (isp == kKaPlusIdx) {
          histos.fill(HIST("h3_AllReco_KaPlus"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == KKPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_KaPlus"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_KaPlus"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaPlus"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaPlus"), multPV, pt, eta);
          }
        } else if (isp == kKaAllIdx) {
          histos.fill(HIST("h3_AllReco_KaAll"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (std::abs(mcP.pdgCode()) == KKPlus) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_KaAll"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_KaAll"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaAll"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaAll"), multPV, pt, eta);
          }
        } else if (isp == kAntiPrIdx) {
          histos.fill(HIST("h3_AllReco_AntiPr"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == -KProton) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_AntiPr"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_AntiPr"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_AntiPr"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_AntiPr"), multPV, pt, eta);
          }
        } else if (isp == kPrIdx) {
          histos.fill(HIST("h3_AllReco_Pr"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (mcP.pdgCode() == KProton) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_Pr"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_Pr"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_Pr"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_Pr"), multPV, pt, eta);
          }
        } else if (isp == kPrAllIdx) {
          histos.fill(HIST("h3_AllReco_PrAll"), multPV, pt, eta);
          if (track.has_mcParticle()) {
            auto mcP = track.mcParticle();
            if (mcP.isPhysicalPrimary()) {
              if (std::abs(mcP.pdgCode()) == KProton) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PrAll"), multPV, mcP.pt(), mcP.eta());
              } else { // Misidentified
                histos.fill(HIST("h3_RecoMatchedToPrimary_MisID_PrAll"), multPV, pt, eta);
              }
            } else {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PrAll"), multPV, pt, eta);
            }
          } else {
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PrAll"), multPV, pt, eta);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate EffWeights", cfgRunGetEff);

  void processMCFlat(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks)
  {
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision))
      return;

    float cent = getCentrality(mcCollision);
    if (cent > KCentMax)
      return;

    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();

    if (!isPassAddPileup(multPV, mcTracks.size(), cent))
      return;
    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), mcCollision.multNTracksPV(), mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);
    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index())
        continue;

      if (!isTrackSelected(track))
        continue;

      float pt = track.pt(), eta = track.eta(), phi = track.phi();
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);
      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;

        float eff = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        float fake = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float w = (eff > KFloatEpsilon) ? (1.0f - fake) / eff : 0.0f;

        if (std::isfinite(w) && w > 0.f) {
          if (isp == kInclusiveIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
          } else if (isp == kPiMinusIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
          } else if (isp == kPiPlusIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
          } else if (isp == kPiAllIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
          } else if (isp == kKaMinusIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
          } else if (isp == kKaPlusIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
          } else if (isp == kKaAllIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
          } else if (isp == kAntiPrIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
          } else if (isp == kPrIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_Pr"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
          } else if (isp == kPrAllIdx) {
            histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd_PrAll"), vz, sign, pt, eta, phi, w);
            histos.fill(HIST("MCReco/hEtaPhiReco_PrAll"), vz, sign, pt, eta, phi, 1.0);
            histos.fill(HIST("MCReco/hEtaPhiRecoWtd_PrAll"), vz, sign, pt, eta, phi, w);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

  void processMCMean(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    double sumWiTruth[KNsp][KNEta]{}, sumWiptiTruth[KNsp][KNEta]{};
    double sumWiReco[KNsp][KNEta]{}, sumWiptiReco[KNsp][KNEta]{};
    double sumWiRecoEffCorr[KNsp][KNEta]{}, sumWiptiRecoEffCorr[KNsp][KNEta]{};
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision))
      return;
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax)
      return;
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent))
      return;

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), mcCollision.multNTracksPV(), mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    memset(sumWiTruth, 0, sizeof(sumWiTruth));
    memset(sumWiptiTruth, 0, sizeof(sumWiptiTruth));
    memset(sumWiReco, 0, sizeof(sumWiReco));
    memset(sumWiptiReco, 0, sizeof(sumWiptiReco));
    memset(sumWiRecoEffCorr, 0, sizeof(sumWiRecoEffCorr));
    memset(sumWiptiRecoEffCorr, 0, sizeof(sumWiptiRecoEffCorr));

    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId())
        continue;

      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
        continue;
      float pt = particle.pt(), eta = particle.eta();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      int pdgCode = particle.pdgCode();
      int absPdg = std::abs(pdgCode);

      bool isSpecies[KNsp] = {
        true,                // kInclusiveIdx
        pdgCode == -KPiPlus, // kPiMinusIdx
        pdgCode == KPiPlus,  // kPiPlusIdx
        absPdg == KPiPlus,   // kPiAllIdx
        pdgCode == -KKPlus,  // kKaMinusIdx
        pdgCode == KKPlus,   // kKaPlusIdx
        absPdg == KKPlus,    // kKaAllIdx
        pdgCode == -KProton, // kAntiPrIdx
        pdgCode == KProton,  // kPrIdx
        absPdg == KProton    // kPrAllIdx
      };

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta])
          continue;

        for (int isp = 0; isp < KNsp; ++isp) {
          if (isSpecies[isp]) {
            sumWiTruth[isp][ieta]++;
            sumWiptiTruth[isp][ieta] += pt;
          }
        }
      }
    }

    for (int isp = 0; isp < KNsp; ++isp) {
      histos.fill(HIST("MCGen/Prof_Cent_Nsp_Nchrec"), cent, isp, sumWiTruth[isp][0]);
      histos.fill(HIST("MCGen/Prof_Mult_Nsp_Nchrec"), multPV, isp, sumWiTruth[isp][0]);
      if (sumWiTruth[isp][0] > 1.0f) {
        histos.fill(HIST("MCGen/Prof_Cent_Nsp_MeanpT"), cent, isp, sumWiptiTruth[isp][0] / sumWiTruth[isp][0]);
        histos.fill(HIST("MCGen/Prof_Mult_Nsp_MeanpT"), multPV, isp, sumWiptiTruth[isp][0] / sumWiTruth[isp][0]);
      }
    }

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index())
        continue;

      if (!isTrackSelected(track))
        continue;
      float pt = track.pt(), eta = track.eta(), phi = track.phi();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);
      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;
        float eff = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        float fake = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PIDIdx>(isp), cfgFlat);
        float w = flatW * (1.0 - fake) / eff;
        if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon)
          continue;

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          sumWiReco[isp][ieta]++;
          sumWiptiReco[isp][ieta] += pt;
          sumWiRecoEffCorr[isp][ieta] += w;
          sumWiptiRecoEffCorr[isp][ieta] += w * pt;
        }

        if (isp == kInclusiveIdx) {
          histos.fill(HIST("Eff_cent"), cent, eff);
          histos.fill(HIST("Fake_cent"), cent, fake);
          histos.fill(HIST("wgt_cent"), cent, w);

          histos.fill(HIST("Eff_Ntrk"), multPV, eff);
          histos.fill(HIST("Fake_Ntrk"), multPV, fake);
          histos.fill(HIST("wgt_Ntrk"), multPV, w);

          histos.fill(HIST("Eff_pT"), pt, eff);
          histos.fill(HIST("Fake_pT"), pt, fake);
          histos.fill(HIST("wgt_pT"), pt, w);

          histos.fill(HIST("Eff_eta"), eta, eff);
          histos.fill(HIST("Fake_eta"), eta, fake);
          histos.fill(HIST("wgt_eta"), eta, w);
        }
        if (isp == kInclusiveIdx) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaAllIdx) {
          histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPrIdx) {
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kAntiPrIdx) {
          histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPrAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PrAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PrAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PrAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        }
      }
    }

    for (int isp = 0; isp < KNsp; ++isp) {
      histos.fill(HIST("MCReco/Prof_Cent_Nsp_Nchrec"), cent, isp, sumWiReco[isp][0]);
      histos.fill(HIST("MCReco/Prof_Mult_Nsp_Nchrec"), multPV, isp, sumWiReco[isp][0]);

      histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nsp_Nchrec"), cent, isp, sumWiRecoEffCorr[isp][0]);
      histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nsp_Nchrec"), multPV, isp, sumWiRecoEffCorr[isp][0]);

      if (sumWiReco[isp][0] > 1.0f) {
        histos.fill(HIST("MCReco/Prof_Cent_Nsp_MeanpT"), cent, isp, sumWiptiReco[isp][0] / sumWiReco[isp][0]);
        histos.fill(HIST("MCReco/Prof_Mult_Nsp_MeanpT"), multPV, isp, sumWiptiReco[isp][0] / sumWiReco[isp][0]);
      }
      if (sumWiRecoEffCorr[isp][0] > 1.0f) {
        histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nsp_MeanpT"), cent, isp, sumWiptiRecoEffCorr[isp][0] / sumWiRecoEffCorr[isp][0]);
        histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nsp_MeanpT"), multPV, isp, sumWiptiRecoEffCorr[isp][0] / sumWiRecoEffCorr[isp][0]);
      }
    }

    for (int ietaA = 0; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 0; ietaC < KNEta; ++ietaC) {
        for (int isp = 0; isp < KNsp; ++isp) {
          float nTruAB = sumWiTruth[isp][ietaA] + sumWiTruth[isp][ietaC];
          float nRecoAB = sumWiReco[isp][ietaA] + sumWiReco[isp][ietaC];
          float nCorrAB = sumWiRecoEffCorr[isp][ietaA] + sumWiRecoEffCorr[isp][ietaC];

          float mptsubTru = (sumWiptiTruth[isp][ietaA] + sumWiptiTruth[isp][ietaC]) / nTruAB;
          float mptsubReco = (sumWiptiReco[isp][ietaA] + sumWiptiReco[isp][ietaC]) / nRecoAB;
          float mptsubRecoEffCorr = (sumWiptiRecoEffCorr[isp][ietaA] + sumWiptiRecoEffCorr[isp][ietaC]) / nCorrAB;

          if (nTruAB > 0) {
            if (isp == kInclusiveIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kPiMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiMinus"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kPiPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiPlus"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kPiAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiAll"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kKaMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaMinus"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kKaPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaPlus"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kKaAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaAll"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_Pr"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kAntiPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_AntiPr"), cent, ietaA, ietaC, mptsubTru);
            else if (isp == kPrAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Tru_PrAll"), cent, ietaA, ietaC, mptsubTru);
          }

          if (nRecoAB > 0) {
            if (isp == kInclusiveIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kPiMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiMinus"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kPiPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiPlus"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kPiAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiAll"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kKaMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaMinus"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kKaPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaPlus"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kKaAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaAll"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_Pr"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kAntiPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_AntiPr"), cent, ietaA, ietaC, mptsubReco);
            else if (isp == kPrAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Reco_PrAll"), cent, ietaA, ietaC, mptsubReco);
          }

          if (nCorrAB > 0) {
            if (isp == kInclusiveIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kPiMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiMinus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kPiPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiPlus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kPiAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiAll"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kKaMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaMinus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kKaPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaPlus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kKaAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaAll"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_Pr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kAntiPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_AntiPr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
            else if (isp == kPrAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PrAll"), cent, ietaA, ietaC, mptsubRecoEffCorr);
          }
        }
      }

      for (int isp = 0; isp < KNsp; ++isp) {
        if (sumWiTruth[isp][ietaA] > 0) {
          float val = sumWiptiTruth[isp][ietaA] / sumWiTruth[isp][ietaA];
          histos.fill(HIST("pmeanTru_nch_etabin_spbin"), multPV, ietaA, isp, val);
          histos.fill(HIST("pmeanMultTru_nch_etabin_spbin"), multPV, ietaA, isp, sumWiTruth[isp][ietaA]);
        }
        if (sumWiReco[isp][ietaA] > 0) {
          float val = sumWiptiReco[isp][ietaA] / sumWiReco[isp][ietaA];
          histos.fill(HIST("pmeanReco_nch_etabin_spbin"), multPV, ietaA, isp, val);
          histos.fill(HIST("pmeanMultReco_nch_etabin_spbin"), multPV, ietaA, isp, sumWiReco[isp][ietaA]);
        }
        if (sumWiRecoEffCorr[isp][ietaA] > 0) {
          float val = sumWiptiRecoEffCorr[isp][ietaA] / sumWiRecoEffCorr[isp][ietaA];
          histos.fill(HIST("pmeanRecoEffcorr_nch_etabin_spbin"), multPV, ietaA, isp, val);
          histos.fill(HIST("pmeanMultRecoEffcorr_nch_etabin_spbin"), multPV, ietaA, isp, sumWiRecoEffCorr[isp][ietaA]);
        }
      }
    } // end ietaA

    double amplFT0A = 0, amplFT0C = 0;
    if (mcCollision.has_foundFT0()) {
      const auto& ft0 = mcCollision.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
        auto eta = getEtaFT0(chanelid, 0);
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        auto eta = getEtaFT0(globalId, 1);
        amplFT0C += ampl;
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, globalId, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, globalId, eta, ampl);
      }
    }

    histos.fill(HIST("pmeanFT0Amultpv"), multPV, amplFT0A);
    histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
    histos.fill(HIST("pmeanFT0Cmultpv"), multPV, amplFT0C);
    histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt and Eff Hists", cfgRunMCMean);

  void processMCFluc(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    if (!state.pmeanTruNchEtabinSpbinStep2 || !state.pmeanRecoNchEtabinSpbinStep2 || !state.pmeanRecoEffcorrNchEtabinSpbinStep2 ||
        !state.pmeanMultTruNchEtabinSpbinStep2 || !state.pmeanMultRecoNchEtabinSpbinStep2 || !state.pmeanMultRecoEffcorrNchEtabinSpbinStep2) {
      LOGF(warning, "MC fluc: Unified Mean pT or Mult map missing");
      return;
    }
    double sumPmwkTru[KNsp][KNEta][KIntM][KIntK]{};
    double sumWkTru[KNsp][KNEta][KIntK]{};
    double sumPmwkReco[KNsp][KNEta][KIntM][KIntK]{};
    double sumWkReco[KNsp][KNEta][KIntK]{};
    double sumPmwkRecoEffCor[KNsp][KNEta][KIntM][KIntK]{};
    double sumWkRecoEffCor[KNsp][KNEta][KIntK]{};

    double meanTru[KNsp][KNEta]{}, c2Tru[KNsp][KNEta]{};
    double meanReco[KNsp][KNEta]{}, c2Reco[KNsp][KNEta]{};
    double meanRecoEffCor[KNsp][KNEta]{}, c2RecoEffCor[KNsp][KNEta]{};

    double meanTruMult[KNsp][KNEta]{};
    double meanRecoMult[KNsp][KNEta]{};
    double meanRecoEffCorMult[KNsp][KNEta]{};

    double p1kBarTru[KNsp][KNEta]{}, p1kBarReco[KNsp][KNEta]{}, p1kBarRecoEffCor[KNsp][KNEta]{};
    double p1kBarTruMult[KNsp][KNEta]{}, p1kBarRecoMult[KNsp][KNEta]{}, p1kBarRecoEffCorMult[KNsp][KNEta]{};

    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision))
      return;
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax)
      return;
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent))
      return;

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    memset(sumPmwkTru, 0, sizeof(sumPmwkTru));
    memset(sumWkTru, 0, sizeof(sumWkTru));
    memset(sumPmwkReco, 0, sizeof(sumPmwkReco));
    memset(sumWkReco, 0, sizeof(sumWkReco));
    memset(sumPmwkRecoEffCor, 0, sizeof(sumPmwkRecoEffCor));
    memset(sumWkRecoEffCor, 0, sizeof(sumWkRecoEffCor));

    memset(meanTru, 0, sizeof(meanTru));
    memset(c2Tru, 0, sizeof(c2Tru));
    memset(meanReco, 0, sizeof(meanReco));
    memset(c2Reco, 0, sizeof(c2Reco));
    memset(meanRecoEffCor, 0, sizeof(meanRecoEffCor));
    memset(c2RecoEffCor, 0, sizeof(c2RecoEffCor));

    memset(meanTruMult, 0, sizeof(meanTruMult));
    memset(meanRecoMult, 0, sizeof(meanRecoMult));
    memset(meanRecoEffCorMult, 0, sizeof(meanRecoEffCorMult));

    memset(p1kBarTru, 0, sizeof(p1kBarTru));
    memset(p1kBarReco, 0, sizeof(p1kBarReco));
    memset(p1kBarRecoEffCor, 0, sizeof(p1kBarRecoEffCor));

    memset(p1kBarTruMult, 0, sizeof(p1kBarTruMult));
    memset(p1kBarRecoMult, 0, sizeof(p1kBarRecoMult));
    memset(p1kBarRecoEffCorMult, 0, sizeof(p1kBarRecoEffCorMult));

    double p1kBarFt0A = 0.0, p1kBarFt0C = 0.0;

    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId())
        continue;

      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
        continue;

      float pt = particle.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      float eta = particle.eta();
      int pdgCode = particle.pdgCode();
      int absPdg = std::abs(pdgCode);

      bool isSpecies[KNsp] = {
        true,                // kInclusiveIdx
        pdgCode == -KPiPlus, // kPiMinusIdx
        pdgCode == KPiPlus,  // kPiPlusIdx
        absPdg == KPiPlus,   // kPiAllIdx
        pdgCode == -KKPlus,  // kKaMinusIdx
        pdgCode == KKPlus,   // kKaPlusIdx
        absPdg == KKPlus,    // kKaAllIdx
        pdgCode == -KProton, // kAntiPrIdx
        pdgCode == KProton,  // kPrIdx
        absPdg == KProton    // kPrAllIdx
      };

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta])
          continue;
        for (int isp = 0; isp < KNsp; ++isp) {
          if (isSpecies[isp]) {
            for (int k = 0; k < KIntK; ++k) {
              for (int m = 0; m < KIntM; ++m) {
                sumPmwkTru[isp][ieta][m][k] += std::pow(pt, m);
              }
              sumWkTru[isp][ieta][k]++;
            }
          }
        }
      }
    } // end truth loop

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index())
        continue;

      if (!isTrackSelected(track))
        continue;

      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      float eta = track.eta();
      float phi = track.phi();
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;
        float eff = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        float fake = getEfficiency(multPV, pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PIDIdx>(isp), cfgFlat);
        float w = flatW * (1.0 - fake) / eff;

        if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon)
          continue;

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int k = 0; k < KIntK; ++k) {
            for (int m = 0; m < KIntM; ++m) {
              sumPmwkReco[isp][ieta][m][k] += std::pow(1.0, k) * std::pow(pt, m);
              sumPmwkRecoEffCor[isp][ieta][m][k] += std::pow(w, k) * std::pow(pt, m);
            }
            sumWkReco[isp][ieta][k] += std::pow(1.0, k);
            sumWkRecoEffCor[isp][ieta][k] += std::pow(w, k);
          }
        }

        if (isp == kInclusiveIdx) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPiAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kKaAllIdx) {
          histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPrIdx) {
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kAntiPrIdx) {
          histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        } else if (isp == kPrAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PrAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoWtd_PrAll"), vz, sign, pt, eta, phi, w);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PrAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
        }
      }
    } // trkslice

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      const int ibx = state.pmeanTruNchEtabinSpbinStep2->GetXaxis()->FindBin(mcCollision.multNTracksPV());
      const int iby = ieta + 1;

      for (int isp = 0; isp < KNsp; ++isp) {
        const int ibz = isp + 1;

        meanTruMult[isp][ieta] = sumWkTru[isp][ieta][1];
        meanRecoMult[isp][ieta] = sumWkReco[isp][ieta][1];
        meanRecoEffCorMult[isp][ieta] = sumWkRecoEffCor[isp][ieta][1];

        float mmptTru = state.pmeanTruNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
        float mmptReco = state.pmeanRecoNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
        float mmptRecoEffCor = state.pmeanRecoEffcorrNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

        float mmMultTru = state.pmeanMultTruNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
        float mmMultReco = state.pmeanMultRecoNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
        float mmMultRecoEffCor = state.pmeanMultRecoEffcorrNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

        if (std::isfinite(mmptTru))
          std::tie(meanTru[isp][ieta], c2Tru[isp][ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTru[isp][ieta], sumWkTru[isp][ieta], mmptTru);
        if (std::isfinite(mmptReco))
          std::tie(meanReco[isp][ieta], c2Reco[isp][ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkReco[isp][ieta], sumWkReco[isp][ieta], mmptReco);
        if (std::isfinite(mmptRecoEffCor))
          std::tie(meanRecoEffCor[isp][ieta], c2RecoEffCor[isp][ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCor[isp][ieta], sumWkRecoEffCor[isp][ieta], mmptRecoEffCor);

        if (mmptTru != 0.0f)
          p1kBarTru[isp][ieta] = meanTru[isp][ieta] - mmptTru;
        if (mmptReco != 0.0f)
          p1kBarReco[isp][ieta] = meanReco[isp][ieta] - mmptReco;
        if (mmptRecoEffCor != 0.0f)
          p1kBarRecoEffCor[isp][ieta] = meanRecoEffCor[isp][ieta] - mmptRecoEffCor;

        if (mmMultTru != 0.0f)
          p1kBarTruMult[isp][ieta] = meanTruMult[isp][ieta] - mmMultTru;
        if (mmMultReco != 0.0f)
          p1kBarRecoMult[isp][ieta] = meanRecoMult[isp][ieta] - mmMultReco;
        if (mmMultRecoEffCor != 0.0f)
          p1kBarRecoEffCorMult[isp][ieta] = meanRecoEffCorMult[isp][ieta] - mmMultRecoEffCor;
      }
    }

    double amplFT0A = 0, amplFT0C = 0;
    if (mcCollision.has_foundFT0()) {
      const auto& ft0 = mcCollision.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        float ampl = ft0.amplitudeC()[iCh];
        amplFT0C += ampl;
      }
    }

    for (int isp = 0; isp < KNsp; ++isp) {
      for (int ieta = 0; ieta < KNEta; ++ieta) {
        histos.fill(HIST("MCGen/Prof_Cent_NEta_Nsp_Nchrec"), cent, ieta, isp, sumWkTru[isp][ieta][1]);
        histos.fill(HIST("MCGen/Prof_Mult_NEta_Nsp_Nchrec"), multPV, ieta, isp, sumWkTru[isp][ieta][1]);

        histos.fill(HIST("MCReco/Prof_Cent_NEta_Nsp_Nchrec"), cent, ieta, isp, sumWkReco[isp][ieta][1]);
        histos.fill(HIST("MCReco/Prof_Mult_NEta_Nsp_Nchrec"), multPV, ieta, isp, sumWkReco[isp][ieta][1]);

        histos.fill(HIST("MCRecoEffCorr/Prof_Cent_NEta_Nsp_Nchrec"), cent, ieta, isp, sumWkRecoEffCor[isp][ieta][1]);
        histos.fill(HIST("MCRecoEffCorr/Prof_Mult_NEta_Nsp_Nchrec"), multPV, ieta, isp, sumWkRecoEffCor[isp][ieta][1]);

        if (sumWkTru[isp][ieta][1] > 1.0f) {
          histos.fill(HIST("MCGen/Prof_Cent_NEta_Nsp_MeanpT"), cent, ieta, isp, meanTru[isp][ieta]);
          histos.fill(HIST("MCGen/Prof_Mult_NEta_Nsp_MeanpT"), multPV, ieta, isp, meanTru[isp][ieta]);
        }
        if (sumWkReco[isp][ieta][1] > 1.0f) {
          histos.fill(HIST("MCReco/Prof_Cent_NEta_Nsp_MeanpT"), cent, ieta, isp, meanReco[isp][ieta]);
          histos.fill(HIST("MCReco/Prof_Mult_NEta_Nsp_MeanpT"), multPV, ieta, isp, meanReco[isp][ieta]);
        }
        if (sumWkRecoEffCor[isp][ieta][1] > 1.0f) {
          histos.fill(HIST("MCRecoEffCorr/Prof_Cent_NEta_Nsp_MeanpT"), cent, ieta, isp, meanRecoEffCor[isp][ieta]);
          histos.fill(HIST("MCRecoEffCorr/Prof_Mult_NEta_Nsp_MeanpT"), multPV, ieta, isp, meanRecoEffCor[isp][ieta]);
        }
      }
    }

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int isp = 0; isp < KNsp; ++isp) {
        if (std::isfinite(meanTru[isp][ieta])) {
          histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanTru[isp][ieta]);
          histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_spbin"), multPV, ieta, isp, meanTru[isp][ieta]);
        }
        if (std::isfinite(c2Tru[isp][ieta])) {
          histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2Tru[isp][ieta]);
          histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_spbin"), multPV, ieta, isp, c2Tru[isp][ieta]);
        }
        if (std::isfinite(meanReco[isp][ieta])) {
          histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanReco[isp][ieta]);
          histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_spbin"), multPV, ieta, isp, meanReco[isp][ieta]);
        }
        if (std::isfinite(c2Reco[isp][ieta])) {
          histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2Reco[isp][ieta]);
          histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_spbin"), multPV, ieta, isp, c2Reco[isp][ieta]);
        }
        if (std::isfinite(meanRecoEffCor[isp][ieta])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanRecoEffCor[isp][ieta]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_spbin"), multPV, ieta, isp, meanRecoEffCor[isp][ieta]);
        }
        if (std::isfinite(c2RecoEffCor[isp][ieta])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2RecoEffCor[isp][ieta]);
          histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_spbin"), multPV, ieta, isp, c2RecoEffCor[isp][ieta]);
        }
      }
    }

    p1kBarFt0A = amplFT0A - state.pmeanFT0AmultpvStep2->GetBinContent(state.pmeanFT0AmultpvStep2->GetXaxis()->FindBin(multPV));
    p1kBarFt0C = amplFT0C - state.pmeanFT0CmultpvStep2->GetBinContent(state.pmeanFT0CmultpvStep2->GetXaxis()->FindBin(multPV));

    for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
      int ietaC = KNEta - ietaA;
      for (int isp = 0; isp < KNsp; ++isp) {
        float c2SubTru = p1kBarTru[isp][ietaA] * p1kBarTru[isp][ietaC];
        float c2SubReco = p1kBarReco[isp][ietaA] * p1kBarReco[isp][ietaC];
        float c2SubRecoEffCor = p1kBarRecoEffCor[isp][ietaA] * p1kBarRecoEffCor[isp][ietaC];

        float covTru = p1kBarTruMult[isp][ietaA] * p1kBarTru[isp][ietaC];
        float covReco = p1kBarRecoMult[isp][ietaA] * p1kBarReco[isp][ietaC];
        float covRecoEffCor = p1kBarRecoEffCorMult[isp][ietaA] * p1kBarRecoEffCor[isp][ietaC];

        float covFT0ATru = p1kBarFt0A * p1kBarTru[isp][ietaC];
        float covFT0AReco = p1kBarFt0A * p1kBarReco[isp][ietaC];
        float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[isp][ietaC];

        float covFT0CTru = p1kBarFt0C * p1kBarTru[isp][ietaA];
        float covFT0CReco = p1kBarFt0C * p1kBarReco[isp][ietaA];
        float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[isp][ietaA];

        if (std::isfinite(c2SubTru)) {
          histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2SubTru);
          histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_spbin"), multPV, ietaA, isp, c2SubTru);
        }
        if (std::isfinite(c2SubReco)) {
          histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2SubReco);
          histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_spbin"), multPV, ietaA, isp, c2SubReco);
        }
        if (std::isfinite(c2SubRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2SubRecoEffCor);
          histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_spbin"), multPV, ietaA, isp, c2SubRecoEffCor);
        }
        if (std::isfinite(covTru)) {
          histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covTru);
          histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_spbin"), multPV, ietaA, isp, covTru);
        }
        if (std::isfinite(covReco)) {
          histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covReco);
          histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_spbin"), multPV, ietaA, isp, covReco);
        }
        if (std::isfinite(covRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covRecoEffCor);
          histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_spbin"), multPV, ietaA, isp, covRecoEffCor);
        }

        if (std::isfinite(covFT0ATru)) {
          histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0ATru);
          histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0ATru);
        }
        if (std::isfinite(covFT0AReco)) {
          histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0AReco);
          histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0AReco);
        }
        if (std::isfinite(covFT0ARecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0ARecoEffCor);
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0ARecoEffCor);
        }

        if (std::isfinite(covFT0CTru)) {
          histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CTru);
          histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0CTru);
        }
        if (std::isfinite(covFT0CReco)) {
          histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CReco);
          histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0CReco);
        }
        if (std::isfinite(covFT0CRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CRecoEffCor);
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_spbin"), multPV, ietaA, isp, covFT0CRecoEffCor);
        }
      }
    }

    for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 1; ietaC < KNEta; ++ietaC) {

        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        float gap = etaValA - etaValB;
        float sum = (etaValA + etaValB);
        for (int isp = 0; isp < KNsp; ++isp) {

          float c2SubTru = p1kBarTru[isp][ietaA] * p1kBarTru[isp][ietaC];
          float c2SubReco = p1kBarReco[isp][ietaA] * p1kBarReco[isp][ietaC];
          float c2SubRecoEffCor = p1kBarRecoEffCor[isp][ietaA] * p1kBarRecoEffCor[isp][ietaC];

          float covTru = p1kBarTruMult[isp][ietaA] * p1kBarTru[isp][ietaC];
          float covReco = p1kBarRecoMult[isp][ietaA] * p1kBarReco[isp][ietaC];
          float covRecoEffCor = p1kBarRecoEffCorMult[isp][ietaA] * p1kBarRecoEffCor[isp][ietaC];

          float covFT0ATru = p1kBarFt0A * p1kBarTru[isp][ietaC];
          float covFT0AReco = p1kBarFt0A * p1kBarReco[isp][ietaC];
          float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[isp][ietaC];

          float covFT0CTru = p1kBarFt0C * p1kBarTru[isp][ietaA];
          float covFT0CReco = p1kBarFt0C * p1kBarReco[isp][ietaA];
          float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[isp][ietaA];

          if (isp == kInclusiveIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kPiMinusIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_PiMinus"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_PiMinus"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_PiMinus"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kPiPlusIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_PiPlus"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_PiPlus"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_PiPlus"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kPiAllIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_PiAll"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_PiAll"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_PiAll"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kKaMinusIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_KaMinus"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_KaMinus"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_KaMinus"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kKaPlusIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_KaPlus"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_KaPlus"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_KaPlus"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kKaAllIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_KaAll"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_KaAll"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_KaAll"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kPrIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_Pr"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_Pr"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_Pr"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kAntiPrIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_AntiPr"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_AntiPr"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_AntiPr"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0CRecoEffCor);

          } else if (isp == kPrAllIdx) {
            if (std::isfinite(c2SubTru)) {
              histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, c2SubTru);
              histos.fill(HIST("MCGen/Prof_GapSum2D_PrAll"), cent, gap, sum, c2SubTru);
            }
            if (std::isfinite(c2SubReco)) {
              histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, c2SubReco);
              histos.fill(HIST("MCReco/Prof_GapSum2D_PrAll"), cent, gap, sum, c2SubReco);
            }
            if (std::isfinite(c2SubRecoEffCor)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, c2SubRecoEffCor);
              histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_PrAll"), cent, gap, sum, c2SubRecoEffCor);
            }

            if (std::isfinite(covTru))
              histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covTru);
            if (std::isfinite(covReco))
              histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covReco);
            if (std::isfinite(covRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covRecoEffCor);

            if (std::isfinite(covFT0ATru))
              histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0ATru);
            if (std::isfinite(covFT0AReco))
              histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0AReco);
            if (std::isfinite(covFT0ARecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0ARecoEffCor);

            if (std::isfinite(covFT0CTru))
              histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0CTru);
            if (std::isfinite(covFT0CReco))
              histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0CReco);
            if (std::isfinite(covFT0CRecoEffCor))
              histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0CRecoEffCor);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt fluc", cfgRunMCFluc);

  void processDataGetNSig(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, AodTracksSel const& tracks)
  {
    histos.fill(HIST("hVtxZ"), coll.posZ());
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;
    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    int ntrk = 0;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;
      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      float eta = track.eta();
      if (eta > etaLw[0] && eta < etaUp[0])
        ntrk++;
      fillNSigmaBefCut(track, cent);
    }

    if (cfgZDC) {
      const auto& foundBC = coll.foundBC_as<BCsRun3>();
      if (!foundBC.has_zdc()) {
        return;
      }
      auto zdc = foundBC.zdc();
      auto zdcAmp = zdc.energyCommonZNA() + zdc.energyCommonZNC();
      histos.fill(HIST("hnTrkPVZDC"), coll.multNTracksPV(), zdcAmp);
      histos.fill(HIST("hNchZDC"), ntrk, zdcAmp);
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataGetNSig, "process data to Get Nsigma cuts", cfgRunDataGetNSig);

  void processGetDataFlat(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, AodTracksSel const& tracks)
  {
    histos.fill(HIST("hVtxZ"), coll.posZ());
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent))
      return;

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    int ntrk = 0;
    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;
      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      float eta = track.eta();
      float phi = track.phi();
      auto sign = track.sign();

      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      if (eta > etaLw[0] && eta < etaUp[0])
        ntrk++;
      fillNSigmaBefCut(track, cent);
      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      fillNSigmaAftCut(track, cent, isSpecies);
      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;
        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        if (eff <= KFloatEpsilon)
          continue;

        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float w = (1.0f - fake) / eff;

        if (!std::isfinite(w) || w <= 0.f)
          continue;

        if (isp == kInclusiveIdx) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaAllIdx) {
          histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPrIdx) {
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
        } else if (isp == kAntiPrIdx) {
          histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPrAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PrAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PrAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PrAll"), vz, sign, pt, eta, phi, w);
        }
      }
    }

    if (cfgZDC) {
      const auto& foundBC = coll.foundBC_as<BCsRun3>();
      if (!foundBC.has_zdc()) {
        return;
      }
      auto zdc = foundBC.zdc();
      auto zdcAmp = zdc.energyCommonZNA() + zdc.energyCommonZNC();
      histos.fill(HIST("hnTrkPVZDC"), coll.multNTracksPV(), zdcAmp);
      histos.fill(HIST("hNchZDC"), ntrk, zdcAmp);
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetDataFlat, "process data to calculate Flattening maps", cfgRunGetDataFlat);

  void processDataMean(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FT0s const&, AodTracksSel const& tracks)
  {
    double sumWi[KNsp][KNEta]{}, sumWipti[KNsp][KNEta]{};

    if (!isEventSelected(coll))
      return;

    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent))
      return;

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;

      float p = track.p();
      float pt = track.pt();
      float eta = track.eta();
      float phi = track.phi();
      auto sign = track.sign();

      if (p < KFloatEpsilon)
        continue;

      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;

      histos.fill(HIST("hP"), p);
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);
      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;
        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        if (eff <= KFloatEpsilon)
          continue;

        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PIDIdx>(isp), cfgFlat);
        float w = flatWeight * (1.0f - fake) / eff;

        if (!std::isfinite(w) || w <= 0.f)
          continue;

        if (isp == kInclusiveIdx) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPiAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaMinusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaPlusIdx) {
          histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
        } else if (isp == kKaAllIdx) {
          histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPrIdx) {
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
        } else if (isp == kAntiPrIdx) {
          histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
        } else if (isp == kPrAllIdx) {
          histos.fill(HIST("hEtaPhiReco_PrAll"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_PrAll"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_PrAll"), vz, sign, pt, eta, phi, w);
        }

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          sumWi[isp][ieta] += w;
          sumWipti[isp][ieta] += w * pt;
        }
      }
    }

    for (int isp = 0; isp < KNsp; ++isp) {
      if (sumWi[isp][0] < 1.0f)
        continue;
      histos.fill(HIST("Prof_Cent_Nsp_Nchrec"), cent, isp, sumWi[isp][0]);
      histos.fill(HIST("Prof_Mult_Nsp_Nchrec"), coll.multNTracksPV(), isp, sumWi[isp][0]);
      histos.fill(HIST("Prof_Cent_Nsp_MeanpT"), cent, isp, sumWipti[isp][0] / sumWi[isp][0]);
      histos.fill(HIST("Prof_Mult_Nsp_MeanpT"), coll.multNTracksPV(), isp, sumWipti[isp][0] / sumWi[isp][0]);
    }

    for (int ietaA = 0; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 0; ietaC < KNEta; ++ietaC) {
        for (int isp = 0; isp < KNsp; ++isp) {
          if ((sumWi[isp][ietaA] < 1.0f) || (sumWi[isp][ietaC] < 1.0f))
            continue;

          double wCorrAB = sumWi[isp][ietaA] + sumWi[isp][ietaC];
          if (wCorrAB > 0) {
            float mptsub = (sumWipti[isp][ietaA] + sumWipti[isp][ietaC]) / wCorrAB;
            if (isp == kInclusiveIdx)
              histos.fill(HIST("Prof2D_MeanpTSub"), cent, ietaA, ietaC, mptsub);
            else if (isp == kPiMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_PiMinus"), cent, ietaA, ietaC, mptsub);
            else if (isp == kPiPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_PiPlus"), cent, ietaA, ietaC, mptsub);
            else if (isp == kPiAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_PiAll"), cent, ietaA, ietaC, mptsub);
            else if (isp == kKaMinusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_KaMinus"), cent, ietaA, ietaC, mptsub);
            else if (isp == kKaPlusIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_KaPlus"), cent, ietaA, ietaC, mptsub);
            else if (isp == kKaAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_KaAll"), cent, ietaA, ietaC, mptsub);
            else if (isp == kPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_Pr"), cent, ietaA, ietaC, mptsub);
            else if (isp == kAntiPrIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_AntiPr"), cent, ietaA, ietaC, mptsub);
            else if (isp == kPrAllIdx)
              histos.fill(HIST("Prof2D_MeanpTSub_PrAll"), cent, ietaA, ietaC, mptsub);
          }
          if (ietaA == ietaC) {
            double mpt = sumWipti[isp][ietaA] / sumWi[isp][ietaA];
            if (sumWi[isp][ietaA] >= 1.0f && std::isfinite(mpt)) {
              histos.fill(HIST("pmean_nch_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, mpt);
              histos.fill(HIST("pmeanMult_nch_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, sumWi[isp][ietaA]);
              histos.fill(HIST("pmean_cent_etabin_spbin"), cent, ietaA, isp, mpt);
              histos.fill(HIST("pmeanMult_cent_etabin_spbin"), cent, ietaA, isp, sumWi[isp][ietaA]);
            }
          }
        }
      }
    }
    double amplFT0A = 0, amplFT0C = 0;
    if (coll.has_foundFT0()) {
      const auto& ft0 = coll.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
        auto eta = getEtaFT0(chanelid, 0);
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        auto eta = getEtaFT0(globalId, 1);
        amplFT0C += ampl;
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, globalId, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, globalId, eta, ampl);
      }
    }

    histos.fill(HIST("pmeanFT0Amultpv"), coll.multNTracksPV(), amplFT0A);
    histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
    histos.fill(HIST("pmeanFT0Cmultpv"), coll.multNTracksPV(), amplFT0C);
    histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataMean, "process data to calculate mean pT", cfgRunDataMean);

  void processDataFluc(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FT0s const&, AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent))
      return;

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    if (!state.pmeanNchEtabinSpbinStep2 || !state.pmeanMultNchEtabinSpbinStep2) {
      LOGF(warning, "Data fluc: Mean pT or Mult map missing");
      return;
    }

    for (int isp = 0; isp < KNsp; ++isp) {
      auto pid = static_cast<PIDIdx>(isp);
      if (!state.hEff[pid] || !state.hFake[pid] || !state.hFlatWeight[pid]) {
        LOGF(warning, "Data fluc: Correction maps (Eff, Fake, or Flat) are null for species index %d", isp);
        return;
      }
    }

    double sumpmwk[KNsp][KNEta][KIntM][KIntK]{};
    double sumwk[KNsp][KNEta][KIntK]{};

    double mean[KNsp][KNEta]{}, c2[KNsp][KNEta]{};
    double p1kBar[KNsp][KNEta]{};
    double meanMult[KNsp][KNEta]{}, p1kBarMult[KNsp][KNEta]{};

    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;

      float p = track.p();
      float pt = track.pt();
      float eta = track.eta();
      float phi = track.phi();
      auto sign = track.sign();

      if (p < KFloatEpsilon)
        continue;

      if (pt <= cfgPtMin || pt > cfgPtMax)
        continue;
      int id = identifyTrack(track, cent);
      bool isPi = (id == KPidPionOne);
      bool isKa = (id == KPidKaonTwo);
      bool isPr = (id == KPidProtonThree);
      bool isSpecies[KNsp] = {
        true,
        isPi && sign < 0, isPi && sign > 0, isPi,
        isKa && sign < 0, isKa && sign > 0, isKa,
        isPr && sign < 0, isPr && sign > 0, isPr};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;
        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 0, cfgEff);
        if (eff <= KFloatEpsilon)
          continue;

        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PIDIdx>(isp), 1, cfgEff);
        float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PIDIdx>(isp), cfgFlat);
        float w = flatWeight * (1.0f - fake) / eff;

        if (!std::isfinite(w) || w <= 0.f)
          continue;

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int k = 0; k < KIntK; ++k) {
            for (int m = 0; m < KIntM; ++m) {
              sumpmwk[isp][ieta][m][k] += std::pow(w, k) * std::pow(pt, m);
            }
            sumwk[isp][ieta][k] += std::pow(w, k);
          }
        }
      }
    }

    double amplFT0A = 0, amplFT0C = 0;
    if (coll.has_foundFT0()) {
      const auto& ft0 = coll.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        float ampl = ft0.amplitudeC()[iCh];
        amplFT0C += ampl;
      }
    }
    double p1kBarFt0A = amplFT0A - state.pmeanFT0AmultpvStep2->GetBinContent(state.pmeanFT0AmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));
    double p1kBarFt0C = amplFT0C - state.pmeanFT0CmultpvStep2->GetBinContent(state.pmeanFT0CmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      const int ibx = state.pmeanNchEtabinSpbinStep2->GetXaxis()->FindBin(coll.multNTracksPV());
      const int iby = ieta + 1;

      for (int isp = 0; isp < KNsp; ++isp) {
        const int ibz = isp + 1;

        float mmpt = state.pmeanNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
        float mmMult = state.pmeanMultNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

        mean[isp][ieta] = sumpmwk[isp][ieta][1][1] / sumwk[isp][ieta][1];
        meanMult[isp][ieta] = sumwk[isp][ieta][1];

        if (std::isfinite(mmpt)) {
          std::tie(mean[isp][ieta], c2[isp][ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwk[isp][ieta], sumwk[isp][ieta], mmpt);
          p1kBar[isp][ieta] = mean[isp][ieta] - mmpt;
        }
        p1kBarMult[isp][ieta] = meanMult[isp][ieta] - mmMult;
      }
    }
    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int isp = 0; isp < KNsp; ++isp) {
        if (std::isfinite(mean[isp][ieta])) {
          histos.fill(HIST("Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, mean[isp][ieta]);
          histos.fill(HIST("Prof_MeanpT_Mult_etabin_spbin"), coll.multNTracksPV(), ieta, isp, mean[isp][ieta]);
        }
        if (std::isfinite(c2[isp][ieta])) {
          histos.fill(HIST("Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2[isp][ieta]);
          histos.fill(HIST("Prof_C2_Mult_etabin_spbin"), coll.multNTracksPV(), ieta, isp, c2[isp][ieta]);
        }
      }
    }

    for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
      int ietaC = KNEta - ietaA;
      for (int isp = 0; isp < KNsp; ++isp) {
        float c2Sub = p1kBar[isp][ietaA] * p1kBar[isp][ietaC];
        float covAC = p1kBarMult[isp][ietaA] * p1kBar[isp][ietaC];
        float covCA = p1kBar[isp][ietaA] * p1kBarMult[isp][ietaC];

        float covFT0A = p1kBarFt0A * p1kBar[isp][ietaC];
        float covFT0C = p1kBarFt0C * p1kBar[isp][ietaA];

        if (std::isfinite(c2Sub)) {
          histos.fill(HIST("Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2Sub);
          histos.fill(HIST("Prof_C2Sub_Mult_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, c2Sub);
        }
        if (std::isfinite(covAC)) {
          histos.fill(HIST("Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covAC);
          histos.fill(HIST("Prof_Cov_Mult_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, covAC);
        }
        if (std::isfinite(covCA)) {
          histos.fill(HIST("Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covCA);
          histos.fill(HIST("Prof_Cov_Mult_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, covCA);
        }
        if (std::isfinite(covFT0A)) {
          histos.fill(HIST("Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0A);
          histos.fill(HIST("Prof_CovFT0A_Mult_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, covFT0A);
        }
        if (std::isfinite(covFT0C)) {
          histos.fill(HIST("Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0C);
          histos.fill(HIST("Prof_CovFT0C_Mult_etabin_spbin"), coll.multNTracksPV(), ietaA, isp, covFT0C);
        }
      }
    }

    for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 1; ietaC < KNEta; ++ietaC) {

        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        float gap = etaValA - etaValB;
        float sum = (etaValA + etaValB);

        for (int isp = 0; isp < KNsp; ++isp) {

          float c2Sub = p1kBar[isp][ietaA] * p1kBar[isp][ietaC];
          float cov = p1kBarMult[isp][ietaA] * p1kBar[isp][ietaC];
          float covFT0A = p1kBarFt0A * p1kBar[isp][ietaC];
          float covFT0C = p1kBarFt0C * p1kBar[isp][ietaA];

          if (isp == kInclusiveIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kPiMinusIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_PiMinus"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_PiMinus"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kPiPlusIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_PiPlus"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_PiPlus"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kPiAllIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_PiAll"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_PiAll"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kKaMinusIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_KaMinus"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_KaMinus"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kKaPlusIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_KaPlus"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_KaPlus"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kKaAllIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_KaAll"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_KaAll"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kPrIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_Pr"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kAntiPrIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_AntiPr"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_AntiPr"), cent, etaValA, etaValB, covFT0C);
          } else if (isp == kPrAllIdx) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, c2Sub);
              histos.fill(HIST("Prof_GapSum2D_PrAll"), cent, gap, sum, c2Sub);
            }
            if (std::isfinite(cov))
              histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, cov);
            if (std::isfinite(covFT0A))
              histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0A);
            if (std::isfinite(covFT0C))
              histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_PrAll"), cent, etaValA, etaValB, covFT0C);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataFluc, "process data to calculate fluc pT", cfgRunDataFluc);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<RadialFlowDecorr>(cfgc)};
  return workflow;
}
