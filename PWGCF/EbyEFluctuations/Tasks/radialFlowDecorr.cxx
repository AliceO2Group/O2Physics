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

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
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

  static constexpr int KnFt0cCell = 96;
  static constexpr int KIntM = 3;
  static constexpr int KIntK = 3;

  static constexpr int KNEta = 17;
  static constexpr int KNpT = 3;

  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr int KPiPlus = 211;
  static constexpr int KKPlus = 321;
  static constexpr int KProton = 2212;
  static constexpr int KNsp = 4;

  static constexpr float KCentTestMin = 10.f;
  static constexpr float KCentTestMaxLo = 60.f;
  static constexpr float KCentTestMaxHi = 70.f;
  static constexpr float KCentCovCut = 1.0f;
  static constexpr float KBinOffset = 0.5f;

  static constexpr float KHalf = 0.5f;
  static constexpr float KPhiMin = 0.f;

  static constexpr int KNbinsZvtx = 240;
  static constexpr float KZvtxMin = -12.f;
  static constexpr float KZvtxMax = 12.f;
  static constexpr int KNbinsP = 100;
  static constexpr float KPMin = 0.f;
  static constexpr float KPMax = 10.f;
  static constexpr int KNbinsPt = 200;
  static constexpr float KPtMin = 0.f;
  static constexpr float KPtMax = 10.f;
  static constexpr int KNbinsEta = 120;
  static constexpr float KEtaMin = -1.2f;
  static constexpr float KEtaMax = 1.2f;
  static constexpr int KNbinsPhi = 64;
  static constexpr float KEtaAxisMin = -0.8f;
  static constexpr float KEtaAxisMax = 0.8f;
  static constexpr int KNbinsPhiFine = 16;

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
  static constexpr int KNbinsPtCoarse = 50;
  static constexpr float KPtMinDefault = 0.2f;
  static constexpr float KPtMidMax = 3.0f;
  static constexpr float KPtHighMax = 5.0f;
  static constexpr float KPtFullMax = 10.0f;
  static constexpr float KCentMax = 90;
  enum PID {
    numKInclusive = 0, // Suffix ""
    numKPion,          // Suffix "_Pi"
    numKKaon,          // Suffix "_Ka"
    numKProton,        // Suffix "_Pr"
    numKNumPID         // Total: 4
  };

  const std::vector<std::string> pidSuffix = {"", "_Pi", "_Ka", "_Pr"};

  enum ECentralityEstimator {
    kCentFT0C = 1,
    kCentFT0A = 2,
    kCentFT0M = 3,
    kCentFV0A = 4
  };
  enum SystemType {
    kPbPb = 1,
    kOO = 2,
    kpPb = 3,
    kpp = 4
  };
  static constexpr float KinvalidCentrality = -1.0f;
  const std::vector<float> etaLw = {
    -0.8,
    -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
  const std::vector<float> etaUp = {
    0.8,
    -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

  const std::vector<float> pTLw = {KPtMinDefault, KPtMinDefault, KPtMinDefault};
  const std::vector<float> pTUp = {KPtMidMax, KPtHighMax, KPtFullMax};

  Configurable<float> cfgVtxZCut{"cfgVtxZCut", 10.f, "z-vertex range"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "min pT"};
  Configurable<float> cfgPtMax{"cfgPtMax", 10.0f, "max pT"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "|η| cut"};
  Configurable<float> cfgDCAXY{"cfgDCAXY", 2.4f, "DCAxy cut"};
  Configurable<float> cfgDCAZ{"cfgDCAZ", 3.2f, "DCAz cut"};
  Configurable<float> cfgTPCClsMin{"cfgTPCClsMin", 70.f, "min TPC clusters"};
  Configurable<float> cfgChi2TPCMax{"cfgChi2TPCMax", 4.0f, "max TPC χ²"};
  Configurable<float> cfgPIDnSigmaCut{"cfgPIDnSigmaCut", 3.f, "TPC PID |nσ| cut"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
  Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};
  Configurable<float> cfgCutTracKDcaMaxZ{"cfgCutTracKDcaMaxZ", 2.0f, "Maximum DcaZ"};
  Configurable<float> cfgCutTracKDcaMaxXY{"cfgCutTracKDcaMaxXY", 0.2f, "Maximum DcaZ"};

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
  Configurable<float> cfgCutEtaLeft{"cfgCutEtaLeft", 0.8f, "Left end of eta gap"};
  Configurable<float> cfgCutEtaRight{"cfgCutEtaRight", 0.8f, "Right end of eta gap"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of subsamples"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 1-->FT0C, 2-->FT0A, 3-->FT0M, 4-->FV0A"};
  Configurable<bool> cfgEvSelkNoSameBunchPileup{"cfgEvSelkNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgEvSelkNoITSROFrameBorder{"cfgEvSelkNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
  Configurable<bool> cfgEvSelkNoTimeFrameBorder{"cfgEvSelkNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
  Configurable<bool> cfgIsGoodZvtxFT0VsPV{"cfgIsGoodZvtxFT0VsPV", true, "Good Vertexing cut"};

  Configurable<int> cfgNchPbMax{"cfgNchPbMax", 4000, "Max Nch range for PbPb collisions"};
  Configurable<int> cfgNchOMax{"cfgNchOMax", 600, "Max Nch range for OO collisions"};

  Configurable<int> cfgSys{"cfgSys", 1, "Efficiency to be used for which system? 1-->PbPb, 2-->OO, 3-->pPb, 4-->pp"};
  Configurable<bool> cfgFlat{"cfgFlat", false, "Whether to use flattening weights"};
  Configurable<bool> cfgEff{"cfgEff", false, "Whether to use Efficiency weights"};
  Configurable<bool> cfgZDC{"cfgZDC", false, "Whether to use ZDC for pileup histograms"};

  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 3.0, 5.0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "centrality axis (percentile)"};

  const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
  const AxisSpec centAxis1Per{101, -0.5, 100.5,
                              "Centrality (%)"
                              "Centrality (%)"};
  AxisSpec nChAxis{1, 0., 1., "Nch", "Nch"};
  AxisSpec nChAxis2{1, 0., 1., "Nch", "Nch"};

  const AxisSpec vzAxis{5, -12.5, 12.5, "Vz"};
  const AxisSpec chgAxis{3, -1.5, 1.5};
  const AxisSpec pTAxis{{0.0, 0.2, 0.5, 1, 3, 5, 7.5, 10}, "pT Axis"};
  const AxisSpec etaAxis{{-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, "Eta"};
  const AxisSpec gapAxis{{-1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55}, "Gaps"};
  const AxisSpec sumAxis{{-0.775, -0.725, -0.675, -0.625, -0.575, -0.525, -0.475, -0.425, -0.375, -0.325, -0.275, -0.225, -0.175, -0.125, -0.075, -0.025, 0.025, 0.075, 0.125, 0.175, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775}, "Sums"};

  Configurable<bool> cfgRunGetEff{"cfgRunGetEff", false, "Run MC pass to build efficiency/fake maps"};
  Configurable<bool> cfgRunGetMCFlat{"cfgRunGetMCFlat", false, "Run MC to Get Flattening Weights"};
  Configurable<bool> cfgRunMCMean{"cfgRunMCMean", false, "Run MC mean(pT) & mean(Et)"};
  Configurable<bool> cfgRunMCFluc{"cfgRunMCFluc", false, "Run MC fluctuations (C2, subevent)"};
  Configurable<bool> cfgRunGetDataFlat{"cfgRunGetDataFlat", false, "Run Data Get Flattening Weights"};
  Configurable<bool> cfgRunDataMean{"cfgRunDataMean", false, "Run DATA mean(pT) & mean(Et)"};
  Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};

  Service<ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::array<TH3F*, numKNumPID> hEff{};
  std::array<TH3F*, numKNumPID> hFake{};
  std::array<THnSparseF*, numKNumPID> hFlatWeight{};

  std::array<TProfile3D*, KNsp> pmeanTruNchEtabinPtbinStep2{};
  std::array<TProfile3D*, KNsp> pmeanRecoNchEtabinPtbinStep2{};
  std::array<TProfile3D*, KNsp> pmeanRecoEffcorrNchEtabinPtbinStep2{};

  std::array<TProfile3D*, KNsp> pmeanMultTruNchEtabinPtbinStep2{};
  std::array<TProfile3D*, KNsp> pmeanMultRecoNchEtabinPtbinStep2{};
  std::array<TProfile3D*, KNsp> pmeanMultRecoEffcorrNchEtabinPtbinStep2{};

  std::array<TProfile3D*, KNsp> pmeanNchEtabinPtbinStep2{};
  std::array<TProfile3D*, KNsp> pmeanMultNchEtabinPtbinStep2{};

  TProfile* pmeanFT0A_multpvStep2 = nullptr;
  TProfile* pmeanFT0C_multpvStep2 = nullptr;
  o2::ft0::Geometry ft0Det;

  template <typename T>
  static std::tuple<float, float, float> getAllCombinedNSigmas(const T& candidate)
  {
    return {
      std::hypot(candidate.tpcNSigmaPr(), candidate.tofNSigmaPr()), // Proton
      std::hypot(candidate.tpcNSigmaPi(), candidate.tofNSigmaPi()), // Pion
      std::hypot(candidate.tpcNSigmaKa(), candidate.tofNSigmaKa())  // Kaon
    };
  }

  template <typename T>
  bool isEventSelected(const T& col)
  {
    if (!col.sel8())
      return false;
    if (std::abs(col.posZ()) > cfgCutVertex)
      return false;
    if (cfgEvSelkNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return false;
    if (cfgUseGoodITSLayerAllCut && !col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return false;
    if (cfgIsGoodZvtxFT0VsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;

    return true;
  }

  template <typename T>
  bool isTrackSelected(const T& trk)
  {
    if (trk.sign() == 0)
      return false;
    if (!trk.has_collision())
      return false;
    if (!trk.isPVContributor())
      return false;
    if (!(trk.itsNCls() > cfgITScluster))
      return false;
    if (!(trk.tpcNClsFound() >= cfgTPCcluster))
      return false;
    if (!(trk.tpcNClsCrossedRows() >= cfgTPCnCrossedRows))
      return false;

    if (trk.pt() < cfgCutPtLower || trk.pt() > cfgCutPtUpper || std::abs(trk.eta()) > cfgCutEta)
      return false;
    if (std::abs(trk.dcaXY()) > cfgCutTracKDcaMaxXY || std::abs(trk.dcaZ()) > cfgCutTracKDcaMaxZ)
      return false;
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

  template <typename T>
  bool selectionProton(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0;

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      auto [combNSigmaPr, combNSigmaPi, combNSigmaKa] = getAllCombinedNSigmas(candidate);
      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaPr > combNSigmaPi) && !(combNSigmaPr > combNSigmaKa)) {
        if (combNSigmaPr < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionPion(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0;

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      auto [combNSigmaPr, combNSigmaPi, combNSigmaKa] = getAllCombinedNSigmas(candidate);
      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaPi > combNSigmaPr) && !(combNSigmaPi > combNSigmaKa)) {
        if (combNSigmaPi < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  template <typename T>
  bool selectionKaon(const T& candidate)
  {
    if (!candidate.hasTPC())
      return false;
    int flag = 0;

    if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
      if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC) {
        flag = 1;
      }
      if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOF) {
        flag = 1;
      }
    }
    if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
      auto [combNSigmaPr, combNSigmaPi, combNSigmaKa] = getAllCombinedNSigmas(candidate);
      int flag2 = 0;
      if (combNSigmaPr < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaPi < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (combNSigmaKa < cfgnSigmaOtherParticles)
        flag2 += 1;
      if (!(flag2 > 1) && !(combNSigmaKa > combNSigmaPi) && !(combNSigmaKa > combNSigmaPr)) {
        if (combNSigmaKa < cfgnSigmaCutCombTPCTOF) {
          flag = 1;
        }
      }
    }
    if (flag == 1)
      return true;
    else
      return false;
  }

  float getCentrality(const auto& col) const
  {
    if (cfgCentralityChoice.value == kCentFT0C)
      return col.centFT0C();
    if (cfgCentralityChoice.value == kCentFT0A)
      return col.centFT0A();
    if (cfgCentralityChoice.value == kCentFT0M)
      return col.centFT0M();
    if (cfgCentralityChoice.value == kCentFV0A)
      return col.centFV0A();
    return KinvalidCentrality;
  }

  float getEfficiency(float mult, float pt, float eta, PID pidType, int effidx, bool cfgEff) const
  {
    if (!cfgEff) {
      if (effidx == 0)
        return 1.0;
      if (effidx == 1)
        return 0.0;
    }
    TH3F* h = nullptr;
    if (effidx == 0)
      h = hEff[pidType];
    if (effidx == 1)
      h = hFake[pidType];

    if (!h)
      return -1;
    const int ibx = h->GetXaxis()->FindBin(mult);
    const int iby = h->GetYaxis()->FindBin(pt);
    const int ibz = h->GetZaxis()->FindBin(eta);
    float val = h->GetBinContent(ibx, iby, ibz);
    return val;
  }

  float getFlatteningWeight(float vz, float chg, float pt, float eta, float phi, PID pidType, bool cfgflat) const
  {
    if (!cfgflat)
      return 1.0;
    THnSparseF* h = hFlatWeight[pidType];

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

    // Fetches from the pre-calculated array, very fast
    auto chPos = ft0Det.getChannelCenter(globalChno);

    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();

    // Force the correct physical Z-direction:
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
    // If the timestamp hasn't changed and pointers are valid, skip the CCDB query
    if (timestamp == mLastTimestamp && offsetFT0 != nullptr) {
      return;
    }

    // Fetch the alignment parameters for the new timestamp
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    // offsetFV0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FV0/Calib/Align", timestamp);

    // Safety checks: Ensure the objects exist and contain A-side and C-side data
    if (offsetFT0 == nullptr || offsetFT0->size() < 2) {
      LOGF(fatal, "Could not load FT0/Calib/Align for timestamp %llu, or vector size < 2", timestamp);
    }
    // if (offsetFV0 == nullptr || offsetFV0->size() < 2) {
    //   LOGF(fatal, "Could not load FV0/Calib/Align for timestamp %llu, or vector size < 2", timestamp);
    // }

    // Update the cached timestamp
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

  using GeneralCollisions = soa::Join<
    aod::Collisions,
    aod::EvSels,
    aod::Mults,
    aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As,
    aod::CentNGlobals>;
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZCut;
  using AodCollisionsSel = soa::Filtered<GeneralCollisions>;

  using UnfilteredTracks = soa::Join<
    aod::Tracks,
    aod::TracksExtra,
    aod::TrackSelection,
    aod::TracksDCA,
    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  Filter trackFilter = nabs(aod::track::eta) < cfgEtaCut &&
                       aod::track::pt > cfgPtMin&&
                                          aod::track::pt < cfgPtMax&&
                                                             nabs(aod::track::dcaXY) < cfgDCAXY&& nabs(aod::track::dcaZ) < cfgDCAZ;
  using AodTracksSel = soa::Filtered<UnfilteredTracks>;
  using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
  using FilteredTCs = soa::Filtered<TCs>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  using MyRun3MCCollisions = soa::Join<
    aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
    aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As,
    aod::CentNGlobals, aod::McCollisionLabels>;

  using MyMCTracks = soa::Join<
    aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
    aod::McTrackLabels,
    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
    aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  PresliceUnsorted<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<MyRun3MCCollisions> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<TCs> trackPerMcParticle = aod::mctracklabel::mcParticleId;
  Preslice<MyMCTracks> perCollision = aod::track::collisionId;
  Preslice<FilteredTCs> trackPerCollision = aod::track::collisionId;

  void declareCommonQA()
  {
    histos.add("hZvtx_after_sel", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hVtxZ", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis1Per}});
    histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis2}, {nChAxis2}});
    histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis2}, {centAxis1Per}});
    histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsP, KPMin, KPMax}});
    histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("hEta", ";#eta", kTH1F, {{KNbinsEta, KEtaMin, KEtaMax}});
    histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});
  }
  void declareMCCommonHists()
  {
    for (const auto& suf : pidSuffix) {
      histos.add("h3_AllPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoMatchedToPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_AllReco" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoUnMatchedToPrimary_Secondary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("h3_RecoUnMatchedToPrimary_Fake" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
      histos.add("hTruth_ParticleWeight" + suf, ";cent;p_{T};#eta", kTH3F, {{centAxis1Per}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    }

    histos.add("ptResolution", ";p_{T}^{MC};p_{T}^{MC}-p_{T}^{reco}", kTH2F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, -KPtResMax, KPtResMax}});
    histos.add("ptTruthReco", ";p_{T}^{MC};p_{T}^{reco}", kTH2F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("etaResolution", ";#eta^{MC};#eta^{MC}-#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsPtRes, -KEtaResMax, KEtaResMax}});
    histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsPtRes, -KEtaFineMax, KEtaFineMax}, {KNbinsPtRes, -KEtaFineMax, KEtaFineMax}});
    histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
    histos.add("vzResolution", ";Vz^{MC};Vz^{MC}-Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, -KVzResMax, KVzResMax}});

    histos.add("h_AllPrimary", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("h_RecoMatchedToPrimary", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("h_RecoUnMatchedToPrimary", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("h_AllReco", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("h_AllRecoEffCorr", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});

    histos.add("hDCAxy_Unmatched", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_Unmatched", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_NotPrimary", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_NotPrimary", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_RecoMatched", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_RecoMatched", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_Reco", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_Reco", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
  }

  void declareMCGetFlatHists()
  {
    for (const auto& suf : pidSuffix) {
      std::string nameEff = "hEtaPhiReco" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;
      std::string nameEffWtd = "hEtaPhiRecoEffWtd" + suf;

      histos.add(nameEffWtd, nameEffWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEff, nameEff.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    }
  }

  void declareMCMeanHists()
  {
    histos.add("Eff_cent", ";cent;#epsilon", kTProfile, {centAxis1Per});
    histos.add("Fake_cent", ";cent;f_{fake}", kTProfile, {centAxis1Per});
    histos.add("wgt_cent", ";cent;w", kTProfile, {centAxis1Per});
    histos.add("Eff_Ntrk", ";N_{PV};#epsilon", kTProfile, {nChAxis2});
    histos.add("Fake_Ntrk", ";N_{PV};f_{fake}", kTProfile, {nChAxis2});
    histos.add("wgt_Ntrk", ";N_{PV};w", kTProfile, {nChAxis2});
    histos.add("Eff_pT", ";p_{T};#epsilon", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("Fake_pT", ";p_{T};f_{fake}", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("wgt_pT", ";p_{T};w", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Eff_eta", ";#eta;#epsilon", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("Fake_eta", ";#eta;f_{fake}", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("wgt_eta", ";#eta;w", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("pmeanFT0A_multpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0C_multpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0C_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});

    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;channel id; #eta;amplitude", kTProfile3D, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;channel id; #eta", kTH3F, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});

    for (const auto& suf : pidSuffix) {
      // Basic Profiles
      histos.add("MCGen/Prof_Cent_Nchrec" + suf, ";cent;#LT N_{PV}#GT", kTProfile, {centAxis1Per});
      histos.add("MCGen/Prof_Mult_Nchrec" + suf, ";N_{PV};#LT N_{PV}#GT", kTProfile, {nChAxis});

      histos.add("MCGen/Prof_Cent_MeanpT" + suf, ";cent;#LT p_{T}#GT", kTProfile, {centAxis1Per});
      histos.add("MCGen/Prof_Mult_MeanpT" + suf, ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});

      histos.add<TProfile3D>("pmeanTruNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanRecoNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanRecoEffcorrNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("pmeanMultTruNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanMultRecoNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanMultRecoEffcorrNchEtabinPtbin" + suf, ";N_{PV};#eta bin;p_{T} bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      for (const int& i : {0, 1, 2}) {
        std::string ptTag = "_ipt" + std::to_string(i);
        histos.add<TProfile3D>("Prof2D_MeanpT_Sub" + ptTag + "_Tru" + suf, ";cent;etaA;etaB", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}});
        histos.add<TProfile3D>("Prof2D_MeanpT_Sub" + ptTag + "_Reco" + suf, ";cent;etaA;etaB", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}});
        histos.add<TProfile3D>("Prof2D_MeanpT_Sub" + ptTag + "_RecoEffCorr" + suf, ";cent;etaA;etaB", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}});
      }
    }

    for (const auto& suf : pidSuffix) {
      std::string nameEff = "hEtaPhiReco" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;
      std::string nameEffWtd = "hEtaPhiRecoEffWtd" + suf;

      histos.add(nameEffWtd, nameEffWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEff, nameEff.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    }
  }

  void declareMCFlucHists()
  {

    for (const auto& suf : pidSuffix) {
      // --- 1D Full Event Calc Profiles ---
      histos.add<TProfile3D>("MCGen/Prof_MeanpT_Cent_etabin_ptbin" + suf, ";cent;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_MeanpT_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("MCGen/Prof_C2_Cent_etabin_ptbin" + suf, ";cent;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_C2_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      // --- 1D Sub-Event Covariances ---
      histos.add<TProfile3D>("MCGen/Prof_C2Sub_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_C2Sub_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("MCGen/Prof_Cov_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_Cov_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("MCGen/Prof_CovFT0A_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_CovFT0A_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_CovFT0C_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("MCGen/Prof_CovFT0C_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      for (const int& i : {0, 1, 2}) {
        std::string ptTag = "_ipt" + std::to_string(i);
        histos.add<TProfile3D>("MCGen/Prof" + ptTag + "_C2Sub2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("MCGen/Prof" + ptTag + "_Cov2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("MCGen/Prof" + ptTag + "_CovFT0A2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("MCGen/Prof" + ptTag + "_CovFT0C2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("MCGen/Prof" + ptTag + "_GapSum2D" + suf, ";cent;#Delta#eta (Gap);#Sigma#eta/2 (Sum)", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
      }

      std::string nameEff = "hEtaPhiReco" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;
      std::string nameEffWtd = "hEtaPhiRecoEffWtd" + suf;
      histos.add(nameEffWtd, nameEffWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEff, nameEff.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

      histos.add("MCGen/Prof_Cent_Nchrec" + suf, ";cent;#LT N_{PV}#GT", kTProfile, {centAxis1Per});
      histos.add("MCGen/Prof_Mult_Nchrec" + suf, ";N_{PV};#LT N_{PV}#GT", kTProfile, {nChAxis});
      histos.add("MCGen/Prof_Cent_MeanpT" + suf, ";cent;#LT p_{T}#GT", kTProfile, {centAxis1Per});
      histos.add("MCGen/Prof_Mult_MeanpT" + suf, ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});
    }
  }

  void declareDataGetFlatHists()
  {
    // 1. Species-dependent Sparse Histograms
    for (const auto& suf : pidSuffix) {
      std::string nameEff = "hEtaPhiReco" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;
      std::string nameEffWtd = "hEtaPhiRecoEffWtd" + suf;

      histos.add(nameEffWtd, nameEffWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEff, nameEff.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    }

    histos.add("hnTrkPVZDC", ";N_{PV};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 3000}});
    histos.add("hNchZDC", ";N_{trk};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 30000}});

    histos.add("hCentnTrk", ";Centrality (%);N_{trk}", kTH2F, {{centAxis1Per}, {nChAxis2}});
    histos.add("hCentnTrkPV", ";Centrality (%);N_{trk, PV}", kTH2F, {{centAxis1Per}, {nChAxis2}});
  }

  void declareDataMeanHists()
  {

    histos.add("pmeanFT0A_multpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0C_multpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0C_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});

    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;channel id; #eta;amplitude", kTProfile3D, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;channel id; #eta", kTH3F, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});

    for (const auto& suf : pidSuffix) {
      std::string nameReco = "hEtaPhiReco" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;
      std::string nameEffWtd = "hEtaPhiRecoEffWtd" + suf;

      histos.add(nameReco, nameReco.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEffWtd, nameEffWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

      histos.add("Prof_Cent_Nchrec" + suf, ";cent;#LT N_{PV}#GT", kTProfile, {centAxis1Per});
      histos.add("Prof_Mult_Nchrec" + suf, ";N_{PV};#LT N_{PV}#GT", kTProfile, {nChAxis});
      histos.add("Prof_Cent_MeanpT" + suf, ";cent;#LT p_{T}#GT", kTProfile, {centAxis1Per});

      histos.add<TProfile3D>("pmean_nch_etabin_ptbin" + suf, ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanMult_nch_etabin_ptbin" + suf, ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("pmean_cent_etabin_ptbin" + suf, ";Centrality (%) ;#eta-bin;p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("pmeanMult_cent_etabin_ptbin" + suf, ";Centrality (%) ;#eta-bin;p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      for (const int& i : {0, 1, 2}) {
        std::string ptTag = "_ipt" + std::to_string(i);
        std::string histName = "Prof2D_MeanpT_Sub" + ptTag + suf;
        histos.add<TProfile3D>(histName, ";cent;#eta_{A} bin;#eta_{B} bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}});
      }
    }
  }

  void declareDataFlucHists()
  {
    for (const auto& suf : pidSuffix) {

      // --- THnSparse QA Histograms ---
      std::string nameReco = "hEtaPhiReco" + suf;
      std::string nameEff = "hEtaPhiRecoEffWtd" + suf;
      std::string nameWtd = "hEtaPhiRecoWtd" + suf;

      histos.add(nameReco, nameReco.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameEff, nameEff.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
      histos.add(nameWtd, nameWtd.c_str(), kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

      // --- 1D Full Event Calc Profiles ---
      histos.add<TProfile3D>("Prof_MeanpT_Cent_etabin_ptbin" + suf, ";cent;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_MeanpT_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("Prof_C2_Cent_etabin_ptbin" + suf, ";cent;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_C2_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      // --- 1D Sub-Event Covariances ---
      histos.add<TProfile3D>("Prof_C2Sub_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_C2Sub_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      histos.add<TProfile3D>("Prof_Cov_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_Cov_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_CovFT0A_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_CovFT0A_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_CovFT0C_Cent_etabin_ptbin" + suf, ";Centrality;#eta-bin; p_{T}-bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
      histos.add<TProfile3D>("Prof_CovFT0C_Mult_etabin_ptbin" + suf, ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

      for (const int& i : {0, 1, 2}) {
        std::string ptTag = "_ipt" + std::to_string(i);
        histos.add<TProfile3D>("Prof" + ptTag + "_C2Sub2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("Prof" + ptTag + "_GapSum2D" + suf, ";cent;#Delta#eta (Gap);#Sigma#eta/2 (Sum)", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
        histos.add<TProfile3D>("Prof" + ptTag + "_Cov2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A} bin;#eta_{C} bin", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("Prof" + ptTag + "_CovFT0A2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
        histos.add<TProfile3D>("Prof" + ptTag + "_CovFT0C2D_Cent_etaA_etaC" + suf, ";cent;#eta_{A};#eta_{B}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      }
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
    auto axPt = hRaw->GetAxis(2);  // Charge
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
      nChAxis = {cfgNchPbMax / 4, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchPbMax / 20, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else if (cfgSys == kOO || cfgSys == kpPb) {
      nChAxis = {cfgNchOMax / 2, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchOMax / 5, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else {
      nChAxis = {cfgNchOMax / 2, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchOMax / 5, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    }

    ccdb->setURL(cfgCCDBurl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    // uint64_t currentTimestamp = bc.timestamp();
    loadAlignParam(now);
    ft0Det.calculateChannelCenter();

    std::string sysDir = "";
    switch (cfgSys) {
      case kPbPb:
        sysDir = "PbPbTest";
        break;
      case kOO:
        sysDir = "OOTest";
        break;
      case kpPb:
        sysDir = "pPbTest";
        break;
      case kpp:
        sysDir = "ppTest";
        break;
      default:
        LOGF(fatal, "Invalid cfgSys value: %d", cfgSys.value);
    }

    std::string pathEff = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_EffMaps";
    std::string pathMCFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_MCFlatMaps";
    std::string pathMCMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_MCMean";
    std::string pathDataFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_DataFlatMaps";
    std::string pathDataMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_DataMean";

    declareCommonQA();
    std::string userCcdbPath;
    if (cfgSys == kPbPb) {
      userCcdbPath = "/Users/s/somadutt/PbPbTest/";
    }
    if (cfgSys == kOO) {
      userCcdbPath = "/Users/s/somadutt/OOTest/";
    }
    if (cfgSys == kpPb) {
      userCcdbPath = "/Users/s/somadutt/pPbTest/";
    }
    if (cfgSys == kpp) {
      userCcdbPath = "/Users/s/somadutt/ppTest/";
    }

    if (cfgRunMCMean || cfgRunMCFluc || cfgRunGetEff) {
      declareMCCommonHists();
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
    if (cfgRunGetDataFlat) {
      declareDataGetFlatHists();
    }
    if (cfgRunGetMCFlat) {
      declareMCGetFlatHists();
    }
    if (cfgRunDataMean) {
      declareDataMeanHists();
    }
    if (cfgRunDataFluc) {
      declareDataFlucHists();
    }

    if (!cfgRunGetEff && (cfgEff)) {
      TList* lst = ccdb->getForTimeStamp<TList>(pathEff, now);

      if (!lst) {
        LOGF(fatal, "Efficiency maps required but CCDB list is null at %s!", pathEff.c_str());
      }

      LOGF(info, "Loading Eff/Fake maps from TList for all species...");

      auto loadEffFakeForPID = [&](PID pidType) {
        std::string suffix = pidSuffix[pidType];
        std::string hEffNumName = "h3_RecoMatchedToPrimary" + suffix;
        std::string hEffDenName = "h3_AllPrimary" + suffix;
        std::string hFakeNumSecName = "h3_RecoUnMatchedToPrimary_Secondary" + suffix;
        std::string hFakeNumFakName = "h3_RecoUnMatchedToPrimary_Fake" + suffix;
        std::string hFakeDenName = "h3_AllReco" + suffix;

        auto* hNum = reinterpret_cast<TH3F*>(lst->FindObject(hEffNumName.c_str()));
        auto* hDen = reinterpret_cast<TH3F*>(lst->FindObject(hEffDenName.c_str()));

        if (hNum && hDen) {
          hEff[pidType] = reinterpret_cast<TH3F*>(hNum->Clone(Form("hEff%s", suffix.c_str())));
          hEff[pidType]->SetDirectory(nullptr);
          hEff[pidType]->Divide(hDen);
        } else {
          LOGF(error, "Missing CCDB objects for efficiency. Checked: %s, %s", hEffNumName.c_str(), hEffDenName.c_str());
        }

        auto* hNumS = reinterpret_cast<TH3F*>(lst->FindObject(hFakeNumSecName.c_str()));
        auto* hNumF = reinterpret_cast<TH3F*>(lst->FindObject(hFakeNumFakName.c_str()));
        auto* hDenF = reinterpret_cast<TH3F*>(lst->FindObject(hFakeDenName.c_str()));

        if (hNumS && hNumF && hDenF) {
          hFake[pidType] = reinterpret_cast<TH3F*>(hNumS->Clone(Form("hFake%s", suffix.c_str())));
          hFake[pidType]->Add(hNumF);
          hFake[pidType]->SetDirectory(nullptr);
          hFake[pidType]->Divide(hDenF);
        } else {
          LOGF(error, "Missing CCDB object(s) for fakes for %s in list.", suffix.c_str());
        }
      };

      // Loop through all PID types: kInclusive, kPion, kKaon, KProton
      for (int i = 0; i < PID::numKNumPID; ++i) {
        loadEffFakeForPID(static_cast<PID>(i));
      }
    }

    if (!cfgRunGetEff && (cfgFlat)) {
      // --- 1. Load Data Flattening Maps ---
      if (cfgRunDataMean || cfgRunDataFluc) {
        LOGF(info, "Data Run: Loading flattening maps from %s", pathDataFlat.c_str());
        TList* lstDataFlat = ccdb->getForTimeStamp<TList>(pathDataFlat, now);

        if (lstDataFlat) {
          // Use a loop to load species-specific flattening weights if they exist in data
          for (int i = 0; i < PID::numKNumPID; ++i) {
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
              hFlatWeight[i] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
            } else {
              LOGF(error, "Data flattening map '%s' not found.", hName.c_str());
            }
          }
        } else {
          LOGF(error, "Could not retrieve Data Flattening TList from: %s", pathDataFlat.c_str());
        }
      }

      // --- 2. Load MC Flattening Maps ---
      if (cfgRunMCMean || cfgRunMCFluc) {
        LOGF(info, "MC Run: Loading flattening maps from %s", pathMCFlat.c_str());
        TList* lstMCFlat = ccdb->getForTimeStamp<TList>(pathMCFlat, now);

        if (lstMCFlat) {
          auto loadFlatForPID = [&](PID pidType) {
            std::string suffix = pidSuffix[pidType];
            std::string hFlatSrcName;
            if (cfgEff && cfgFlat) {
              hFlatSrcName = "hEtaPhiRecoWtd" + suffix;
            } else if (cfgEff) {
              hFlatSrcName = "hEtaPhiRecoEffWtd" + suffix;
            } else {
              hFlatSrcName = "hEtaPhiReco" + suffix;
            }

            auto* hRaw = reinterpret_cast<THnSparseF*>(lstMCFlat->FindObject(hFlatSrcName.c_str()));

            if (hRaw) {
              hFlatWeight[pidType] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
            } else {
              LOGF(warning, "MC flattening source '%s' not found in list.", hFlatSrcName.c_str());
            }
          };

          for (int i = 0; i < PID::numKNumPID; ++i) {
            loadFlatForPID(static_cast<PID>(i));
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
        loadTProfileFromList(lstMCMean, "pmeanFT0A_multpv", pmeanFT0A_multpvStep2);
        loadTProfileFromList(lstMCMean, "pmeanFT0C_multpv", pmeanFT0C_multpvStep2);

        for (int isp = 0; isp < KNsp; ++isp) {
          std::string suf = pidSuffix[isp];
          loadTProfile3DFromList(lstMCMean, ("pmeanTruNchEtabinPtbin" + suf).c_str(), pmeanTruNchEtabinPtbinStep2[isp]);
          loadTProfile3DFromList(lstMCMean, ("pmeanRecoNchEtabinPtbin" + suf).c_str(), pmeanRecoNchEtabinPtbinStep2[isp]);
          loadTProfile3DFromList(lstMCMean, ("pmeanRecoEffcorrNchEtabinPtbin" + suf).c_str(), pmeanRecoEffcorrNchEtabinPtbinStep2[isp]);

          loadTProfile3DFromList(lstMCMean, ("pmeanMultTruNchEtabinPtbin" + suf).c_str(), pmeanMultTruNchEtabinPtbinStep2[isp]);
          loadTProfile3DFromList(lstMCMean, ("pmeanMultRecoNchEtabinPtbin" + suf).c_str(), pmeanMultRecoNchEtabinPtbinStep2[isp]);
          loadTProfile3DFromList(lstMCMean, ("pmeanMultRecoEffcorrNchEtabinPtbin" + suf).c_str(), pmeanMultRecoEffcorrNchEtabinPtbinStep2[isp]);
        }
      } else {
        LOGF(error, "Could not retrieve TList for MC Mean from: %s", pathMCMean.c_str());
      }
    }

    if (cfgRunDataFluc) {
      LOGF(info, "Loading Data Mean profiles from CCDB path: %s", pathDataMean.c_str());
      TList* lstDataMean = ccdb->getForTimeStamp<TList>(pathDataMean, now);

      loadTProfileFromList(lstDataMean, "pmeanFT0A_multpv", pmeanFT0A_multpvStep2);
      loadTProfileFromList(lstDataMean, "pmeanFT0C_multpv", pmeanFT0C_multpvStep2);
      if (lstDataMean) {
        for (int isp = 0; isp < KNsp; ++isp) {
          std::string suf = pidSuffix[isp];
          loadTProfile3DFromList(lstDataMean, ("pmean_nch_etabin_ptbin" + suf).c_str(), pmeanNchEtabinPtbinStep2[isp]);
          loadTProfile3DFromList(lstDataMean, ("pmeanMult_nch_etabin_ptbin" + suf).c_str(), pmeanMultNchEtabinPtbinStep2[isp]);
        }
      } else {
        LOGF(error, "Could not retrieve TList for Data Mean from: %s", pathDataMean.c_str());
      }
    }
    LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
  }

  void processGetEffHists(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, /*soa::SmallGroups<MyRun3MCCollisions> const& collisions,*/ TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision())
          continue;
        histos.fill(HIST("hVtxZ"), col.posZ());
        if (!isEventSelected(col))
          continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (trackSlice.size() < 1 || partSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;
        float multPV = col.multNTracksPV();

        histos.fill(HIST("hZvtx_after_sel"), col.posZ());
        histos.fill(HIST("hCentrality"), cent);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, tracks.size());
        histos.fill(HIST("Hist2D_cent_nch"), multPV, cent);

        // --- Denominator: Truth Particles ---
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
            continue;

          const int absPdgId = std::abs(particle.pdgCode());
          float pt = particle.pt();
          float eta = particle.eta();
          float w = particle.weight();

          // Inclusive (Denominator)
          histos.fill(HIST("hTruth_ParticleWeight"), cent, pt, eta, w);
          histos.fill(HIST("h3_AllPrimary"), multPV, pt, eta);
          histos.fill(HIST("h_AllPrimary"), pt);

          // Species Specific Denominators
          if (absPdgId == KPiPlus) {
            histos.fill(HIST("hTruth_ParticleWeight_Pi"), cent, pt, eta, w);
            histos.fill(HIST("h3_AllPrimary_Pi"), multPV, pt, eta);
          } else if (absPdgId == KKPlus) {
            histos.fill(HIST("hTruth_ParticleWeight_Ka"), cent, pt, eta, w);
            histos.fill(HIST("h3_AllPrimary_Ka"), multPV, pt, eta);
          } else if (absPdgId == KProton) {
            histos.fill(HIST("hTruth_ParticleWeight_Pr"), cent, pt, eta, w);
            histos.fill(HIST("h3_AllPrimary_Pr"), multPV, pt, eta);
          }
        }

        // --- Numerator and Fakes: Reconstructed Tracks ---
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          float pt = track.pt();
          float eta = track.eta();
          float phi = track.phi();

          bool isPi = selectionPion(track);
          bool isKa = selectionKaon(track);
          bool isPr = selectionProton(track);

          // Inclusive QA
          histos.fill(HIST("h_AllReco"), pt);
          histos.fill(HIST("h3_AllReco"), multPV, pt, eta);
          histos.fill(HIST("hEtaPhiReco"), col.posZ(), track.sign(), pt, eta, phi);

          // Species QA (Fills based on your PID selection)
          if (isPi) {
            histos.fill(HIST("h3_AllReco_Pi"), multPV, pt, eta);
          }
          if (isKa) {
            histos.fill(HIST("h3_AllReco_Ka"), multPV, pt, eta);
          }
          if (isPr) {
            histos.fill(HIST("h3_AllReco_Pr"), multPV, pt, eta);
          }

          if (track.has_mcParticle()) {
            auto mcPart2 = track.mcParticle();
            const int absPdgId = std::abs(mcPart2.pdgCode());

            if (mcPart2.isPhysicalPrimary()) {

              histos.fill(HIST("ptResolution"), mcPart2.pt(), (pt - mcPart2.pt()) / mcPart2.pt());
              histos.fill(HIST("etaResolution"), mcPart2.eta(), eta - mcPart2.eta());
              histos.fill(HIST("etaTruthReco"), mcPart2.eta(), eta);
              histos.fill(HIST("vzResolution"), mcPart2.vz(), (col.posZ() - mcPart2.vz()) / mcPart2.vz());
              histos.fill(HIST("TruthTracKVz"), mcPart2.vz(), col.posZ());

              // Reconstructed Numerator (Inclusive)
              histos.fill(HIST("h3_RecoMatchedToPrimary"), multPV, mcPart2.pt(), mcPart2.eta());
              histos.fill(HIST("h_RecoMatchedToPrimary"), pt);

              // Species Matching (Efficiency Numerator)
              // We fill ONLY if the reconstructed PID matches the Truth PDG
              if (isPi && absPdgId == KPiPlus)
                histos.fill(HIST("h3_RecoMatchedToPrimary_Pi"), multPV, mcPart2.pt(), mcPart2.eta());
              if (isKa && absPdgId == KKPlus)
                histos.fill(HIST("h3_RecoMatchedToPrimary_Ka"), multPV, mcPart2.pt(), mcPart2.eta());
              if (isPr && absPdgId == KProton)
                histos.fill(HIST("h3_RecoMatchedToPrimary_Pr"), multPV, mcPart2.pt(), mcPart2.eta());

            } else {
              // Secondary (Contamination)
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), multPV, pt, eta);
              histos.fill(HIST("h_RecoUnMatchedToPrimary"), pt);
              if (isPi)
                histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_Pi"), multPV, pt, eta);
              if (isKa)
                histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_Ka"), multPV, pt, eta);
              if (isPr)
                histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_Pr"), multPV, pt, eta);
            }
          } else {
            // Fake Tracks (No MC matching)
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), multPV, pt, eta);
            if (isPi)
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_Pi"), multPV, pt, eta);
            if (isKa)
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_Ka"), multPV, pt, eta);
            if (isPr)
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_Pr"), multPV, pt, eta);
          }
        }
      }
    }
    LOGF(info, "FINISHED RUNNING processGetEffHists");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate Eff and Fakes", cfgRunGetEff);

  void processMCFlat(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, /*soa::SmallGroups<MyRun3MCCollisions> const& collisions,*/ TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision() || !isEventSelected(col))
          continue;
        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (trackSlice.size() < 1 || partSlice.size() < 1)
          continue;

        float multPV = col.multNTracksPV();
        float vz = col.posZ();

        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          float pt = track.pt();
          float eta = track.eta();
          float phi = track.phi();
          auto sign = track.sign();

          // 1. Inclusive Weighting (Always filled for selected tracks)
          float effIncl = getEfficiency(multPV, pt, eta, numKInclusive, 0, cfgEff);
          float fakeIncl = getEfficiency(multPV, pt, eta, numKInclusive, 1, cfgEff);
          float wIncl = (1.0 - fakeIncl) / effIncl;

          if (std::isfinite(wIncl) && wIncl > 0.f && effIncl > KFloatEpsilon) {
            histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, wIncl);
            histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi, 1.0);
          }

          // 2. Pion Weighting
          if (selectionPion(track)) {
            float effPi = getEfficiency(multPV, pt, eta, numKPion, 0, cfgEff);
            float fakePi = getEfficiency(multPV, pt, eta, numKPion, 1, cfgEff);
            float wPi = (1.0 - fakePi) / effPi;
            if (std::isfinite(wPi) && wPi > 0.f && effPi > KFloatEpsilon) {
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, wPi);
              histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi, 1.0);
            }
          }

          // 3. Kaon Weighting
          if (selectionKaon(track)) {
            float effKa = getEfficiency(multPV, pt, eta, numKKaon, 0, cfgEff);
            float fakeKa = getEfficiency(multPV, pt, eta, numKKaon, 1, cfgEff);
            float wKa = (1.0 - fakeKa) / effKa;
            if (std::isfinite(wKa) && wKa > 0.f && effKa > KFloatEpsilon) {
              histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, wKa);
              histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi, 1.0);
            }
          }

          // 4. Proton Weighting
          if (selectionProton(track)) {
            float effPr = getEfficiency(multPV, pt, eta, numKProton, 0, cfgEff);
            float fakePr = getEfficiency(multPV, pt, eta, numKProton, 1, cfgEff);
            float wPr = (1.0 - fakePr) / effPr;
            if (std::isfinite(wPr) && wPr > 0.f && effPr > KFloatEpsilon) {
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, wPr);
              histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi, 1.0);
            }
          }
        } // end track loop
      } // end col loop
    } // end mcColl loop
    LOGF(info, "FINISHED RUNNING processMCFlat");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

  void processMCMean(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    // Track-sum arrays using KNsp index (isp=0: Incl, 1: Pi, 2: Ka, 3: Pr)
    double sumWiTruth[KNsp][KNEta][KNpT]{}, sumWiptiTruth[KNsp][KNEta][KNpT]{};
    double sumWiReco[KNsp][KNEta][KNpT]{}, sumWiptiReco[KNsp][KNEta][KNpT]{};
    double sumWiRecoEffCorr[KNsp][KNEta][KNpT]{}, sumWiptiRecoEffCorr[KNsp][KNEta][KNpT]{};

    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision() || !isEventSelected(col))
          continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (trackSlice.size() < 1 || partSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;
        float multPV = col.multNTracksPV();
        float vz = col.posZ();

        // Reset local event sum
        memset(sumWiTruth, 0, sizeof(sumWiTruth));
        memset(sumWiptiTruth, 0, sizeof(sumWiptiTruth));
        memset(sumWiReco, 0, sizeof(sumWiReco));
        memset(sumWiptiReco, 0, sizeof(sumWiptiReco));
        memset(sumWiRecoEffCorr, 0, sizeof(sumWiRecoEffCorr));
        memset(sumWiptiRecoEffCorr, 0, sizeof(sumWiptiRecoEffCorr));

        // --- 1. Truth Loop ---
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
            continue;
          float pt = particle.pt(), eta = particle.eta();
          const int absPdgId = std::abs(particle.pdgCode());
          bool isSpecies[KNsp] = {true, (absPdgId == KPiPlus), (absPdgId == KKPlus), (absPdgId == KProton)};
          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              for (int isp = 0; isp < KNsp; ++isp) {
                if (isSpecies[isp]) {
                  sumWiTruth[isp][ieta][ipt]++;
                  sumWiptiTruth[isp][ieta][ipt] += pt;
                }
              }
            }
          }
        }

        for (int isp = 0; isp < KNsp; ++isp) {
          if (isp == numKInclusive) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec"), cent, sumWiTruth[0][0][0]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec"), multPV, sumWiTruth[0][0][0]);
            if (sumWiTruth[0][0][0] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT"), cent, sumWiptiTruth[0][0][0] / sumWiTruth[0][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT"), multPV, sumWiptiTruth[0][0][0] / sumWiTruth[0][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Pi"), cent, sumWiTruth[1][0][0]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Pi"), multPV, sumWiTruth[1][0][0]);

            if (sumWiTruth[1][0][0] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Pi"), cent, sumWiptiTruth[1][0][0] / sumWiTruth[1][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Pi"), multPV, sumWiptiTruth[1][0][0] / sumWiTruth[1][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Ka"), cent, sumWiTruth[2][0][0]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Ka"), multPV, sumWiTruth[2][0][0]);

            if (sumWiTruth[2][0][0] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Ka"), cent, sumWiptiTruth[2][0][0] / sumWiTruth[2][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Ka"), multPV, sumWiptiTruth[2][0][0] / sumWiTruth[2][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Pr"), cent, sumWiTruth[3][0][0]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Pr"), multPV, sumWiTruth[3][0][0]);

            if (sumWiTruth[3][0][0] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Pr"), cent, sumWiptiTruth[3][0][0] / sumWiTruth[3][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Pr"), multPV, sumWiptiTruth[3][0][0] / sumWiTruth[3][0][0]);
            }
          }
        }

        // --- 2. Reco Loop ---
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;
          float pt = track.pt(), eta = track.eta(), phi = track.phi();
          auto sign = track.sign();
          bool isSpecies[KNsp] = {true, selectionPion(track), selectionKaon(track), selectionProton(track)};
          for (int isp = 0; isp < KNsp; ++isp) {
            if (!isSpecies[isp])
              continue;
            float eff = getEfficiency(multPV, pt, eta, static_cast<PID>(isp), 0, cfgEff);
            float fake = getEfficiency(multPV, pt, eta, static_cast<PID>(isp), 1, cfgEff);
            float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
            float w = flatW * (1.0 - fake) / eff;
            if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon)
              continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                sumWiReco[isp][ieta][ipt]++;
                sumWiptiReco[isp][ieta][ipt] += pt;
                sumWiRecoEffCorr[isp][ieta][ipt] += w;
                sumWiptiRecoEffCorr[isp][ieta][ipt] += w * pt;

                if (ipt == 0) {
                  // Fill profiles vs. Centrality
                  histos.fill(HIST("Eff_cent"), cent, eff);
                  histos.fill(HIST("Fake_cent"), cent, fake);
                  histos.fill(HIST("wgt_cent"), cent, w);

                  // Fill profiles vs. Multiplicity (Ntrk)
                  histos.fill(HIST("Eff_Ntrk"), multPV, eff);
                  histos.fill(HIST("Fake_Ntrk"), multPV, fake);
                  histos.fill(HIST("wgt_Ntrk"), multPV, w);

                  // Fill profiles vs. pT
                  histos.fill(HIST("Eff_pT"), pt, eff);
                  histos.fill(HIST("Fake_pT"), pt, fake);
                  histos.fill(HIST("wgt_pT"), pt, w);

                  // Fill profiles vs. Eta
                  histos.fill(HIST("Eff_eta"), eta, eff);
                  histos.fill(HIST("Fake_eta"), eta, fake);
                  histos.fill(HIST("wgt_eta"), eta, w);
                }
              }
            }

            if (isp == numKInclusive) {

              histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKPion) { // Pion
              histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Pi"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKKaon) { // Kaon
              histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Ka"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKProton) { // Proton
              histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
            }
          }
        }

        for (int isp = 0; isp < KNsp; ++isp) {
          if (isp == numKInclusive) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec"), cent, sumWiReco[0][0][0]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec"), multPV, sumWiReco[0][0][0]);
            if (sumWiReco[0][0][0] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT"), cent, sumWiptiReco[0][0][0] / sumWiReco[0][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT"), multPV, sumWiptiReco[0][0][0] / sumWiReco[0][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Pi"), cent, sumWiReco[1][0][0]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Pi"), multPV, sumWiReco[1][0][0]);

            if (sumWiReco[1][0][0] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Pi"), cent, sumWiptiReco[1][0][0] / sumWiReco[1][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Pi"), multPV, sumWiptiReco[1][0][0] / sumWiReco[1][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Ka"), cent, sumWiReco[2][0][0]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Ka"), multPV, sumWiReco[2][0][0]);

            if (sumWiReco[2][0][0] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Ka"), cent, sumWiptiReco[2][0][0] / sumWiReco[2][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Ka"), multPV, sumWiptiReco[2][0][0] / sumWiReco[2][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Pr"), cent, sumWiReco[3][0][0]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Pr"), multPV, sumWiReco[3][0][0]);
            if (sumWiReco[3][0][0] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Pr"), cent, sumWiptiReco[3][0][0] / sumWiReco[3][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Pr"), multPV, sumWiptiReco[3][0][0] / sumWiReco[3][0][0]);
            }
          }

          if (isp == numKInclusive) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec"), cent, sumWiRecoEffCorr[0][0][0]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec"), multPV, sumWiRecoEffCorr[0][0][0]);
            if (sumWiRecoEffCorr[0][0][0] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT"), cent, sumWiptiRecoEffCorr[0][0][0] / sumWiRecoEffCorr[0][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT"), multPV, sumWiptiRecoEffCorr[0][0][0] / sumWiRecoEffCorr[0][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Pi"), cent, sumWiRecoEffCorr[1][0][0]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Pi"), multPV, sumWiRecoEffCorr[1][0][0]);
            if (sumWiRecoEffCorr[1][0][0] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Pi"), cent, sumWiptiRecoEffCorr[1][0][0] / sumWiRecoEffCorr[1][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Pi"), multPV, sumWiptiRecoEffCorr[1][0][0] / sumWiRecoEffCorr[1][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Ka"), cent, sumWiRecoEffCorr[2][0][0]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Ka"), multPV, sumWiRecoEffCorr[2][0][0]);
            if (sumWiRecoEffCorr[2][0][0] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Ka"), cent, sumWiptiRecoEffCorr[2][0][0] / sumWiRecoEffCorr[2][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Ka"), multPV, sumWiptiRecoEffCorr[2][0][0] / sumWiRecoEffCorr[2][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Pr"), cent, sumWiRecoEffCorr[3][0][0]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Pr"), multPV, sumWiRecoEffCorr[3][0][0]);
            if (sumWiRecoEffCorr[3][0][0] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Pr"), cent, sumWiptiRecoEffCorr[3][0][0] / sumWiRecoEffCorr[3][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Pr"), multPV, sumWiptiRecoEffCorr[3][0][0] / sumWiRecoEffCorr[3][0][0]);
            }
          }
        }

        for (int ietaA = 0; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 0; ietaC < KNEta; ++ietaC) {
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              for (int isp = 0; isp < KNsp; ++isp) {

                // 1. Truth Sub-event Mean
                double nTruAB = sumWiTruth[isp][ietaA][ipt] + sumWiTruth[isp][ietaC][ipt];
                if (nTruAB > 0) {
                  float mptsubTru = (sumWiptiTruth[isp][ietaA][ipt] + sumWiptiTruth[isp][ietaC][ipt]) / nTruAB;
                  if (isp == numKInclusive) { // Inclusive
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Tru"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Tru"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Tru"), cent, ietaA, ietaC, mptsubTru);
                  } else if (isp == numKPion) { // Pion
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Tru_Pi"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Tru_Pi"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Tru_Pi"), cent, ietaA, ietaC, mptsubTru);
                  } else if (isp == numKKaon) { // Kaon
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Tru_Ka"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Tru_Ka"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Tru_Ka"), cent, ietaA, ietaC, mptsubTru);
                  } else if (isp == numKProton) { // Proton
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Tru_Pr"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Tru_Pr"), cent, ietaA, ietaC, mptsubTru);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Tru_Pr"), cent, ietaA, ietaC, mptsubTru);
                  }
                }

                // 2. Reco Raw Sub-event Mean
                double nRecAB = sumWiReco[isp][ietaA][ipt] + sumWiReco[isp][ietaC][ipt];
                if (nRecAB > 0) {
                  float mptsubReco = (sumWiptiReco[isp][ietaA][ipt] + sumWiptiReco[isp][ietaC][ipt]) / nRecAB;
                  if (isp == numKInclusive) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Reco"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Reco"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Reco"), cent, ietaA, ietaC, mptsubReco);
                  } else if (isp == numKPion) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Reco_Pi"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Reco_Pi"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Reco_Pi"), cent, ietaA, ietaC, mptsubReco);
                  } else if (isp == numKKaon) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Reco_Ka"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Reco_Ka"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Reco_Ka"), cent, ietaA, ietaC, mptsubReco);
                  } else if (isp == numKProton) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Reco_Pr"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Reco_Pr"), cent, ietaA, ietaC, mptsubReco);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Reco_Pr"), cent, ietaA, ietaC, mptsubReco);
                  }
                }

                // 3. Reco Efficiency Corrected Sub-event Mean
                double wCorrAB = sumWiRecoEffCorr[isp][ietaA][ipt] + sumWiRecoEffCorr[isp][ietaC][ipt];
                if (wCorrAB > 0) {
                  float mptsubRecoEffCorr = (sumWiptiRecoEffCorr[isp][ietaA][ipt] + sumWiptiRecoEffCorr[isp][ietaC][ipt]) / wCorrAB;
                  if (isp == numKInclusive) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_RecoEffCorr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_RecoEffCorr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_RecoEffCorr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                  } else if (isp == numKPion) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_RecoEffCorr_Pi"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_RecoEffCorr_Pi"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_RecoEffCorr_Pi"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                  } else if (isp == numKKaon) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_RecoEffCorr_Ka"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_RecoEffCorr_Ka"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_RecoEffCorr_Ka"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                  } else if (isp == numKProton) {
                    if (ipt == 0)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_RecoEffCorr_Pr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 2)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_RecoEffCorr_Pr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    if (ipt == KNpT - 1)
                      histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_RecoEffCorr_Pr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                  }
                }

                // 4. pmean Profiles (Individual Bins)
                if (ietaA == ietaC) { // only fill once per eta bin
                  if (sumWiTruth[isp][ietaA][ipt] > 0) {
                    float val = sumWiptiTruth[isp][ietaA][ipt] / sumWiTruth[isp][ietaA][ipt];
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanTruNchEtabinPtbin"), multPV, ietaA, ipt, val);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanTruNchEtabinPtbin_Pi"), multPV, ietaA, ipt, val);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanTruNchEtabinPtbin_Ka"), multPV, ietaA, ipt, val);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanTruNchEtabinPtbin_Pr"), multPV, ietaA, ipt, val);
                  }

                  if (sumWiReco[isp][ietaA][ipt] > 0) {
                    float val = sumWiptiReco[isp][ietaA][ipt] / sumWiReco[isp][ietaA][ipt];
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanRecoNchEtabinPtbin"), multPV, ietaA, ipt, val);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanRecoNchEtabinPtbin_Pi"), multPV, ietaA, ipt, val);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanRecoNchEtabinPtbin_Ka"), multPV, ietaA, ipt, val);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanRecoNchEtabinPtbin_Pr"), multPV, ietaA, ipt, val);
                  }
                  if (sumWiRecoEffCorr[isp][ietaA][ipt] > 0) {
                    float val = sumWiptiRecoEffCorr[isp][ietaA][ipt] / sumWiRecoEffCorr[isp][ietaA][ipt];
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin"), multPV, ietaA, ipt, val);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin_Pi"), multPV, ietaA, ipt, val);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin_Ka"), multPV, ietaA, ipt, val);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin_Pr"), multPV, ietaA, ipt, val);
                  }

                  if (sumWiTruth[isp][ietaA][ipt] > 0) {
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanMultTruNchEtabinPtbin"), multPV, ietaA, ipt, sumWiTruth[isp][ietaA][ipt]);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanMultTruNchEtabinPtbin_Pi"), multPV, ietaA, ipt, sumWiTruth[isp][ietaA][ipt]);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanMultTruNchEtabinPtbin_Ka"), multPV, ietaA, ipt, sumWiTruth[isp][ietaA][ipt]);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanMultTruNchEtabinPtbin_Pr"), multPV, ietaA, ipt, sumWiTruth[isp][ietaA][ipt]);
                  }

                  if (sumWiReco[isp][ietaA][ipt] > 0) {
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanMultRecoNchEtabinPtbin"), multPV, ietaA, ipt, sumWiReco[isp][ietaA][ipt]);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanMultRecoNchEtabinPtbin_Pi"), multPV, ietaA, ipt, sumWiReco[isp][ietaA][ipt]);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanMultRecoNchEtabinPtbin_Ka"), multPV, ietaA, ipt, sumWiReco[isp][ietaA][ipt]);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanMultRecoNchEtabinPtbin_Pr"), multPV, ietaA, ipt, sumWiReco[isp][ietaA][ipt]);
                  }
                  if (sumWiRecoEffCorr[isp][ietaA][ipt] > 0) {
                    if (isp == numKInclusive)
                      histos.fill(HIST("pmeanMultRecoEffcorrNchEtabinPtbin"), multPV, ietaA, ipt, sumWiRecoEffCorr[isp][ietaA][ipt]);
                    else if (isp == numKPion)
                      histos.fill(HIST("pmeanMultRecoEffcorrNchEtabinPtbin_Pi"), multPV, ietaA, ipt, sumWiRecoEffCorr[isp][ietaA][ipt]);
                    else if (isp == numKKaon)
                      histos.fill(HIST("pmeanMultRecoEffcorrNchEtabinPtbin_Ka"), multPV, ietaA, ipt, sumWiRecoEffCorr[isp][ietaA][ipt]);
                    else if (isp == numKProton)
                      histos.fill(HIST("pmeanMultRecoEffcorrNchEtabinPtbin_Pr"), multPV, ietaA, ipt, sumWiRecoEffCorr[isp][ietaA][ipt]);
                  }
                }
              } // end isp
            } // end ipt
          } // end ietaC
        } // end ietaA

        double amplFT0A = 0, amplFT0C = 0;
        if (col.has_foundFT0()) {
          const auto& ft0 = col.foundFT0();
          for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
            auto chanelid = ft0.channelA()[iCh];
            float ampl = ft0.amplitudeA()[iCh];
            amplFT0A += ampl;
            auto eta = getEtaFT0(chanelid, 0);
            histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, iCh, eta, ampl);
            histos.fill(HIST("h3_cent_id_eta_FT0"), cent, iCh, eta, ampl);
          }
          for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
            auto chanelid = ft0.channelC()[iCh];
            auto globalId = chanelid + KnFt0cCell;
            float ampl = ft0.amplitudeC()[iCh];
            auto eta = getEtaFT0(globalId, 1);
            amplFT0C += ampl;
            histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, iCh, eta, ampl);
            histos.fill(HIST("h3_cent_id_eta_FT0"), cent, iCh, eta, ampl);
          }
        }

        histos.fill(HIST("pmeanFT0A_multpv"), multPV, amplFT0A);
        histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
        histos.fill(HIST("pmeanFT0C_multpv"), multPV, amplFT0C);
        histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt/Et and Eff Hists", cfgRunMCMean);

  void processMCFluc(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    // 1. Safety Check: Step 2 Mean Maps
    for (int isp = 0; isp < KNsp; ++isp) {
      if (!pmeanTruNchEtabinPtbinStep2[isp] || !pmeanRecoNchEtabinPtbinStep2[isp] || !pmeanRecoEffcorrNchEtabinPtbinStep2[isp] ||
          !pmeanMultTruNchEtabinPtbinStep2[isp] || !pmeanMultRecoNchEtabinPtbinStep2[isp] || !pmeanMultRecoEffcorrNchEtabinPtbinStep2[isp]) {
        LOGF(warning, "MC fluc: Mean pT or Mult map missing for species index %d", isp);
        return;
      }
    }

    // Expanded with KNsp index (isp=0: Incl, 1: Pi, 2: Ka, 3: Pr)
    double sumPmwkTru[KNsp][KNEta][KNpT][KIntM][KIntK]{};
    double sumWkTru[KNsp][KNEta][KNpT][KIntK]{};
    double sumPmwkReco[KNsp][KNEta][KNpT][KIntM][KIntK]{};
    double sumWkReco[KNsp][KNEta][KNpT][KIntK]{};
    double sumPmwkRecoEffCor[KNsp][KNEta][KNpT][KIntM][KIntK]{};
    double sumWkRecoEffCor[KNsp][KNEta][KNpT][KIntK]{};

    double meanTru[KNsp][KNEta][KNpT]{}, c2Tru[KNsp][KNEta][KNpT]{};
    double meanReco[KNsp][KNEta][KNpT]{}, c2Reco[KNsp][KNEta][KNpT]{};
    double meanRecoEffCor[KNsp][KNEta][KNpT]{}, c2RecoEffCor[KNsp][KNEta][KNpT]{};

    double meanTruMult[KNsp][KNEta][KNpT]{};
    double meanRecoMult[KNsp][KNEta][KNpT]{};
    double meanRecoEffCorMult[KNsp][KNEta][KNpT]{};

    double p1kBarTru[KNsp][KNEta][KNpT]{}, p1kBarReco[KNsp][KNEta][KNpT]{}, p1kBarRecoEffCor[KNsp][KNEta][KNpT]{};
    double p1kBarTruMult[KNsp][KNEta][KNpT]{}, p1kBarRecoMult[KNsp][KNEta][KNpT]{}, p1kBarRecoEffCorMult[KNsp][KNEta][KNpT]{};

    for (const auto& mcCollision : mcColl) {
      auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision() || !isEventSelected(col))
          continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;
        float multPV = col.multNTracksPV();
        // Reset local arrays
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

        // --- 1. Truth Loop ---
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
            continue;

          float pt = particle.pt();
          float eta = particle.eta();
          const int absPdgId = std::abs(particle.pdgCode());
          bool isSpecies[KNsp] = {true, (absPdgId == KPiPlus), (absPdgId == KKPlus), (absPdgId == KProton)};

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;

              for (int isp = 0; isp < KNsp; ++isp) {
                if (isSpecies[isp]) {
                  for (int k = 0; k < KIntK; ++k) {
                    for (int m = 0; m < KIntM; ++m) {
                      sumPmwkTru[isp][ieta][ipt][m][k] += std::pow(pt, m);
                    }
                    sumWkTru[isp][ieta][ipt][k]++;
                  }
                }
              }
            }
          }
        } // end truth loop

        // --- 2. Reco Loop ---
        float vz = col.posZ();
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          float pt = track.pt();
          float eta = track.eta();
          float phi = track.phi();
          float sign = track.sign();
          bool isSpecies[KNsp] = {true, selectionPion(track), selectionKaon(track), selectionProton(track)};

          for (int isp = 0; isp < KNsp; ++isp) {
            if (!isSpecies[isp])
              continue;

            float eff = getEfficiency(col.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
            float fake = getEfficiency(col.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
            float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
            float w = flatW * (1.0 - fake) / eff;

            if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon)
              continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                for (int k = 0; k < KIntK; ++k) {
                  for (int m = 0; m < KIntM; ++m) {
                    sumPmwkReco[isp][ieta][ipt][m][k] += std::pow(1.0, k) * std::pow(pt, m);
                    sumPmwkRecoEffCor[isp][ieta][ipt][m][k] += std::pow(w, k) * std::pow(pt, m);
                  }
                  sumWkReco[isp][ieta][ipt][k] += std::pow(1.0, k);
                  sumWkRecoEffCor[isp][ieta][ipt][k] += std::pow(w, k);
                }
              }
            }

            if (isp == numKInclusive) {
              histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKPion) { // Pion
              histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Pi"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKKaon) { // Kaon
              histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Ka"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);

            } else if (isp == numKProton) { // Proton
              histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
              histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
              histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
            }
          }
        } // trkslice

        // --- 3. FullEvent calculation & Covariances ---
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          for (int ipt = 0; ipt < KNpT; ++ipt) {

            // Safely get the X-axis bin using the Inclusive map [0]
            const int ibx = pmeanTruNchEtabinPtbinStep2[0]->GetXaxis()->FindBin(col.multNTracksPV());
            const int iby = ieta + 1;
            const int ibz = ipt + 1;

            for (int isp = 0; isp < KNsp; ++isp) {
              meanTruMult[isp][ieta][ipt] = sumWkTru[isp][ieta][ipt][1];
              meanRecoMult[isp][ieta][ipt] = sumWkReco[isp][ieta][ipt][1];
              meanRecoEffCorMult[isp][ieta][ipt] = sumWkRecoEffCor[isp][ieta][ipt][1];

              // Dynamically fetch from the arrays using the 'isp' index!
              float mmptTru = pmeanTruNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);
              float mmptReco = pmeanRecoNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);
              float mmptRecoEffCor = pmeanRecoEffcorrNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);

              float mmMultTru = pmeanMultTruNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);
              float mmMultReco = pmeanMultRecoNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);
              float mmMultRecoEffCor = pmeanMultRecoEffcorrNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);

              if (std::isfinite(mmptTru))
                std::tie(meanTru[isp][ieta][ipt], c2Tru[isp][ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTru[isp][ieta][ipt], sumWkTru[isp][ieta][ipt], mmptTru);
              if (std::isfinite(mmptReco))
                std::tie(meanReco[isp][ieta][ipt], c2Reco[isp][ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkReco[isp][ieta][ipt], sumWkReco[isp][ieta][ipt], mmptReco);
              if (std::isfinite(mmptRecoEffCor))
                std::tie(meanRecoEffCor[isp][ieta][ipt], c2RecoEffCor[isp][ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCor[isp][ieta][ipt], sumWkRecoEffCor[isp][ieta][ipt], mmptRecoEffCor);

              if (mmptTru != 0.0f)
                p1kBarTru[isp][ieta][ipt] = meanTru[isp][ieta][ipt] - mmptTru;
              if (mmptReco != 0.0f)
                p1kBarReco[isp][ieta][ipt] = meanReco[isp][ieta][ipt] - mmptReco;
              if (mmptRecoEffCor != 0.0f)
                p1kBarRecoEffCor[isp][ieta][ipt] = meanRecoEffCor[isp][ieta][ipt] - mmptRecoEffCor;

              if (mmMultTru != 0.0f)
                p1kBarTruMult[isp][ieta][ipt] = meanTruMult[isp][ieta][ipt] - mmMultTru;
              if (mmMultReco != 0.0f)
                p1kBarRecoMult[isp][ieta][ipt] = meanRecoMult[isp][ieta][ipt] - mmMultReco;
              if (mmMultRecoEffCor != 0.0f)
                p1kBarRecoEffCorMult[isp][ieta][ipt] = meanRecoEffCorMult[isp][ieta][ipt] - mmMultRecoEffCor;
            }
          }
        }

        double amplFT0A = 0, amplFT0C = 0;
        if (col.has_foundFT0()) {
          const auto& ft0 = col.foundFT0();
          for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
            // auto chanelid = ft0.channelA()[iCh];
            float ampl = ft0.amplitudeA()[iCh];
            amplFT0A += ampl;
          }
          for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
            // auto chanelid = ft0.channelC()[iCh];
            // auto globalId = chanelid + KnFt0cCell;
            float ampl = ft0.amplitudeC()[iCh];
            // auto eta = getEtaFT0(globalId, 1);
            amplFT0C += ampl;
          }
        }

        for (int isp = 0; isp < KNsp; ++isp) {
          if (isp == numKInclusive) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec"), cent, sumWkTru[isp][0][0][1]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec"), multPV, sumWkTru[isp][0][0][1]);
            if (sumWkTru[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT"), cent, meanTru[isp][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT"), multPV, meanTru[isp][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Pi"), cent, sumWkTru[isp][0][0][1]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Pi"), multPV, sumWkTru[isp][0][0][1]);

            if (sumWkTru[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Pi"), cent, meanTru[isp][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Pi"), multPV, meanTru[isp][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Ka"), cent, sumWkTru[isp][0][0][1]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Ka"), multPV, sumWkTru[isp][0][0][1]);

            if (sumWkTru[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Ka"), cent, meanTru[isp][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Ka"), multPV, meanTru[isp][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCGen/Prof_Cent_Nchrec_Pr"), cent, sumWkTru[isp][0][0][1]);
            histos.fill(HIST("MCGen/Prof_Mult_Nchrec_Pr"), multPV, sumWkTru[isp][0][0][1]);

            if (sumWkTru[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCGen/Prof_Cent_MeanpT_Pr"), cent, meanTru[isp][0][0]);
              histos.fill(HIST("MCGen/Prof_Mult_MeanpT_Pr"), multPV, meanTru[isp][0][0]);
            }
          }
        }

        for (int isp = 0; isp < KNsp; ++isp) {
          if (isp == numKInclusive) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec"), cent, sumWkReco[isp][0][0][1]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec"), multPV, sumWkReco[isp][0][0][1]);
            if (sumWkReco[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT"), cent, meanReco[isp][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT"), multPV, meanReco[isp][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Pi"), cent, sumWkReco[isp][0][0][1]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Pi"), multPV, sumWkReco[isp][0][0][1]);

            if (sumWkReco[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Pi"), cent, meanReco[isp][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Pi"), multPV, meanReco[isp][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Ka"), cent, sumWkReco[isp][0][0][1]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Ka"), multPV, sumWkReco[isp][0][0][1]);

            if (sumWkReco[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Ka"), cent, meanReco[isp][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Ka"), multPV, meanReco[isp][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCReco/Prof_Cent_Nchrec_Pr"), cent, sumWkReco[isp][0][0][1]);
            histos.fill(HIST("MCReco/Prof_Mult_Nchrec_Pr"), multPV, sumWkReco[isp][0][0][1]);
            if (sumWkReco[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCReco/Prof_Cent_MeanpT_Pr"), cent, meanReco[isp][0][0]);
              histos.fill(HIST("MCReco/Prof_Mult_MeanpT_Pr"), multPV, meanReco[isp][0][0]);
            }
          }

          if (isp == numKInclusive) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec"), cent, sumWkRecoEffCor[isp][0][0][1]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec"), multPV, sumWkRecoEffCor[isp][0][0][1]);
            if (sumWkRecoEffCor[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT"), cent, meanRecoEffCor[isp][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT"), multPV, meanRecoEffCor[isp][0][0]);
            }

          } else if (isp == numKPion) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Pi"), cent, sumWkRecoEffCor[isp][0][0][1]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Pi"), multPV, sumWkRecoEffCor[isp][0][0][1]);
            if (sumWkRecoEffCor[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Pi"), cent, meanRecoEffCor[isp][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Pi"), multPV, meanRecoEffCor[isp][0][0]);
            }

          } else if (isp == numKKaon) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Ka"), cent, sumWkRecoEffCor[isp][0][0][1]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Ka"), multPV, sumWkRecoEffCor[isp][0][0][1]);
            if (sumWkRecoEffCor[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Ka"), cent, meanRecoEffCor[isp][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Ka"), multPV, meanRecoEffCor[isp][0][0]);
            }
          } else if (isp == numKProton) {
            histos.fill(HIST("MCRecoEffCorr/Prof_Cent_Nchrec_Pr"), cent, sumWkRecoEffCor[isp][0][0][1]);
            histos.fill(HIST("MCRecoEffCorr/Prof_Mult_Nchrec_Pr"), multPV, sumWkRecoEffCor[isp][0][0][1]);
            if (sumWkRecoEffCor[isp][0][0][1] > 1.0f) {
              histos.fill(HIST("MCRecoEffCorr/Prof_Cent_MeanpT_Pr"), cent, meanRecoEffCor[isp][0][0]);
              histos.fill(HIST("MCRecoEffCorr/Prof_Mult_MeanpT_Pr"), multPV, meanRecoEffCor[isp][0][0]);
            }
          }
        }

        // --- 3. Fill 1D Profiles: Gen, Reco, and EffCorr Levels ---
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            for (int isp = 0; isp < KNsp; ++isp) {

              if (isp == numKInclusive) { // Inclusive (No suffix)
                // --- MCGen (Truth) ---
                if (std::isfinite(meanTru[0][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_ptbin"), cent, ieta, ipt, meanTru[0][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, meanTru[0][ieta][ipt]);
                }
                if (std::isfinite(c2Tru[0][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_ptbin"), cent, ieta, ipt, c2Tru[0][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, c2Tru[0][ieta][ipt]);
                }
                // --- MCReco ---
                if (std::isfinite(meanReco[0][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_ptbin"), cent, ieta, ipt, meanReco[0][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, meanReco[0][ieta][ipt]);
                }
                if (std::isfinite(c2Reco[0][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_ptbin"), cent, ieta, ipt, c2Reco[0][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, c2Reco[0][ieta][ipt]);
                }
                // --- MCRecoEffCorr ---
                if (std::isfinite(meanRecoEffCor[0][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_ptbin"), cent, ieta, ipt, meanRecoEffCor[0][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, meanRecoEffCor[0][ieta][ipt]);
                }
                if (std::isfinite(c2RecoEffCor[0][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_ptbin"), cent, ieta, ipt, c2RecoEffCor[0][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_ptbin"), col.multNTracksPV(), ieta, ipt, c2RecoEffCor[0][ieta][ipt]);
                }

              } else if (isp == numKPion) { // Pions (_Pi)
                // --- MCGen (Truth) ---
                if (std::isfinite(meanTru[1][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, meanTru[1][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, meanTru[1][ieta][ipt]);
                }
                if (std::isfinite(c2Tru[1][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, c2Tru[1][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, c2Tru[1][ieta][ipt]);
                }
                // --- MCReco ---
                if (std::isfinite(meanReco[1][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, meanReco[1][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, meanReco[1][ieta][ipt]);
                }
                if (std::isfinite(c2Reco[1][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, c2Reco[1][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, c2Reco[1][ieta][ipt]);
                }
                // --- MCRecoEffCorr ---
                if (std::isfinite(meanRecoEffCor[1][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, meanRecoEffCor[1][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, meanRecoEffCor[1][ieta][ipt]);
                }
                if (std::isfinite(c2RecoEffCor[1][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, c2RecoEffCor[1][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ieta, ipt, c2RecoEffCor[1][ieta][ipt]);
                }

              } else if (isp == numKKaon) { // Kaons (_Ka)
                // --- MCGen (Truth) ---
                if (std::isfinite(meanTru[2][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, meanTru[2][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, meanTru[2][ieta][ipt]);
                }
                if (std::isfinite(c2Tru[2][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, c2Tru[2][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, c2Tru[2][ieta][ipt]);
                }
                // --- MCReco ---
                if (std::isfinite(meanReco[2][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, meanReco[2][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, meanReco[2][ieta][ipt]);
                }
                if (std::isfinite(c2Reco[2][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, c2Reco[2][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, c2Reco[2][ieta][ipt]);
                }
                // --- MCRecoEffCorr ---
                if (std::isfinite(meanRecoEffCor[2][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, meanRecoEffCor[2][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, meanRecoEffCor[2][ieta][ipt]);
                }
                if (std::isfinite(c2RecoEffCor[2][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, c2RecoEffCor[2][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ieta, ipt, c2RecoEffCor[2][ieta][ipt]);
                }

              } else if (isp == numKProton) { // Protons (_Pr)
                // --- MCGen (Truth) ---
                if (std::isfinite(meanTru[3][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, meanTru[3][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, meanTru[3][ieta][ipt]);
                }
                if (std::isfinite(c2Tru[3][ieta][ipt])) {
                  histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, c2Tru[3][ieta][ipt]);
                  histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, c2Tru[3][ieta][ipt]);
                }
                // --- MCReco ---
                if (std::isfinite(meanReco[3][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, meanReco[3][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, meanReco[3][ieta][ipt]);
                }
                if (std::isfinite(c2Reco[3][ieta][ipt])) {
                  histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, c2Reco[3][ieta][ipt]);
                  histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, c2Reco[3][ieta][ipt]);
                }
                // --- MCRecoEffCorr ---
                if (std::isfinite(meanRecoEffCor[3][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, meanRecoEffCor[3][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, meanRecoEffCor[3][ieta][ipt]);
                }
                if (std::isfinite(c2RecoEffCor[3][ieta][ipt])) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, c2RecoEffCor[3][ieta][ipt]);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ieta, ipt, c2RecoEffCor[3][ieta][ipt]);
                }
              }
            }
          }
        }

        p1kBarFt0A = amplFT0A - pmeanFT0A_multpvStep2->GetBinContent(pmeanFT0A_multpvStep2->GetXaxis()->FindBin(col.multNTracksPV()));
        p1kBarFt0C = amplFT0C - pmeanFT0C_multpvStep2->GetBinContent(pmeanFT0C_multpvStep2->GetXaxis()->FindBin(col.multNTracksPV()));

        // --- 4. Symmetric Sub-Event (1D) Covariances ---
        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            for (int isp = 0; isp < KNsp; ++isp) {
              float c2SubTru = p1kBarTru[isp][ietaA][ipt] * p1kBarTru[isp][ietaC][ipt];
              float c2SubReco = p1kBarReco[isp][ietaA][ipt] * p1kBarReco[isp][ietaC][ipt];
              float c2SubRecoEffCor = p1kBarRecoEffCor[isp][ietaA][ipt] * p1kBarRecoEffCor[isp][ietaC][ipt];

              float covTru = p1kBarTruMult[isp][ietaA][ipt] * p1kBarTru[isp][ietaC][ipt];
              float covReco = p1kBarRecoMult[isp][ietaA][ipt] * p1kBarReco[isp][ietaC][ipt];
              float covRecoEffCor = p1kBarRecoEffCorMult[isp][ietaA][ipt] * p1kBarRecoEffCor[isp][ietaC][ipt];

              float covFT0ATru = p1kBarFt0A * p1kBarTru[isp][ietaC][ipt];
              float covFT0AReco = p1kBarFt0A * p1kBarReco[isp][ietaC][ipt];
              float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[isp][ietaC][ipt];

              float covFT0CTru = p1kBarFt0C * p1kBarTru[isp][ietaA][ipt];
              float covFT0CReco = p1kBarFt0C * p1kBarReco[isp][ietaA][ipt];
              float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[isp][ietaA][ipt];

              if (isp == numKInclusive) {
                if (std::isfinite(c2SubTru)) {
                  histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin_ptbin"), cent, ietaA, ipt, c2SubTru);
                  histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubTru);
                }
                if (std::isfinite(c2SubReco)) {
                  histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_ptbin"), cent, ietaA, ipt, c2SubReco);
                  histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubReco);
                }
                if (std::isfinite(c2SubRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_ptbin"), cent, ietaA, ipt, c2SubRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubRecoEffCor);
                }
                if (std::isfinite(covTru)) {
                  histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_ptbin"), cent, ietaA, ipt, covTru);
                  histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covTru);
                }
                if (std::isfinite(covReco)) {
                  histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_ptbin"), cent, ietaA, ipt, covReco);
                  histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covReco);
                }
                if (std::isfinite(covRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_ptbin"), cent, ietaA, ipt, covRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covRecoEffCor);
                }

                if (std::isfinite(covFT0ATru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0ATru);
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0ATru);
                }
                if (std::isfinite(covFT0AReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0AReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0AReco);
                }
                if (std::isfinite(covFT0ARecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0ARecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0ARecoEffCor);
                }

                if (std::isfinite(covFT0CTru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0CTru);
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0CTru);
                }
                if (std::isfinite(covFT0CReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0CReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0CReco);
                }
                if (std::isfinite(covFT0CRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0CRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, covFT0CRecoEffCor);
                }

              } else if (isp == numKPion) { // Pion
                if (std::isfinite(c2SubTru)) {
                  histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, c2SubTru);
                  histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, c2SubTru);
                }
                if (std::isfinite(c2SubReco)) {
                  histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, c2SubReco);
                  histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, c2SubReco);
                }
                if (std::isfinite(c2SubRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, c2SubRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, c2SubRecoEffCor);
                }
                if (std::isfinite(covTru)) {
                  histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covTru);
                  histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covTru);
                }
                if (std::isfinite(covReco)) {
                  histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covReco);
                  histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covReco);
                }
                if (std::isfinite(covRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covRecoEffCor);
                }

                if (std::isfinite(covFT0ATru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0ATru);
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0ATru);
                }
                if (std::isfinite(covFT0AReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0AReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0AReco);
                }
                if (std::isfinite(covFT0ARecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0ARecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0ARecoEffCor);
                }

                if (std::isfinite(covFT0CTru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0CTru);
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0CTru);
                }
                if (std::isfinite(covFT0CReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0CReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0CReco);
                }
                if (std::isfinite(covFT0CRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0CRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_ptbin_Pi"), col.multNTracksPV(), ietaA, ipt, covFT0CRecoEffCor);
                }

              } else if (isp == numKKaon) { // Kaon
                if (std::isfinite(c2SubTru)) {
                  histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, c2SubTru);
                  histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, c2SubTru);
                }
                if (std::isfinite(c2SubReco)) {
                  histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, c2SubReco);
                  histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, c2SubReco);
                }
                if (std::isfinite(c2SubRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, c2SubRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, c2SubRecoEffCor);
                }
                if (std::isfinite(covTru)) {
                  histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covTru);
                  histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covTru);
                }
                if (std::isfinite(covReco)) {
                  histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covReco);
                  histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covReco);
                }
                if (std::isfinite(covRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covRecoEffCor);
                }

                if (std::isfinite(covFT0ATru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0ATru);
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0ATru);
                }
                if (std::isfinite(covFT0AReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0AReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0AReco);
                }
                if (std::isfinite(covFT0ARecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0ARecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0ARecoEffCor);
                }

                if (std::isfinite(covFT0CTru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0CTru);
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0CTru);
                }
                if (std::isfinite(covFT0CReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0CReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0CReco);
                }
                if (std::isfinite(covFT0CRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0CRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_ptbin_Ka"), col.multNTracksPV(), ietaA, ipt, covFT0CRecoEffCor);
                }

              } else if (isp == numKProton) { // Proton
                if (std::isfinite(c2SubTru)) {
                  histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, c2SubTru);
                  histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, c2SubTru);
                }
                if (std::isfinite(c2SubReco)) {
                  histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, c2SubReco);
                  histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, c2SubReco);
                }
                if (std::isfinite(c2SubRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, c2SubRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, c2SubRecoEffCor);
                }
                if (std::isfinite(covTru)) {
                  histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covTru);
                  histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covTru);
                }
                if (std::isfinite(covReco)) {
                  histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covReco);
                  histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covReco);
                }
                if (std::isfinite(covRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covRecoEffCor);
                }

                if (std::isfinite(covFT0ATru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0ATru);
                  histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0ATru);
                }
                if (std::isfinite(covFT0AReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0AReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0AReco);
                }
                if (std::isfinite(covFT0ARecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0ARecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0ARecoEffCor);
                }

                if (std::isfinite(covFT0CTru)) {
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0CTru);
                  histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0CTru);
                }
                if (std::isfinite(covFT0CReco)) {
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0CReco);
                  histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0CReco);
                }
                if (std::isfinite(covFT0CRecoEffCor)) {
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0CRecoEffCor);
                  histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_ptbin_Pr"), col.multNTracksPV(), ietaA, ipt, covFT0CRecoEffCor);
                }
              }
            }
          }
        }

        // --- 5. Full 2D Covariances & GapSum2D Profiles ---
        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {

            // Gap and Sum calculations
            float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
            float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
            float gap = etaValA - etaValB;
            float sum = (etaValA + etaValB) / 2.0f;

            for (int ipt = 0; ipt < KNpT; ++ipt) {
              for (int isp = 0; isp < KNsp; ++isp) {

                float c2SubTru = p1kBarTru[isp][ietaA][ipt] * p1kBarTru[isp][ietaC][ipt];
                float c2SubReco = p1kBarReco[isp][ietaA][ipt] * p1kBarReco[isp][ietaC][ipt];
                float c2SubRecoEffCor = p1kBarRecoEffCor[isp][ietaA][ipt] * p1kBarRecoEffCor[isp][ietaC][ipt];

                float covTru = p1kBarTruMult[isp][ietaA][ipt] * p1kBarTru[isp][ietaC][ipt];
                float covReco = p1kBarRecoMult[isp][ietaA][ipt] * p1kBarReco[isp][ietaC][ipt];
                float covRecoEffCor = p1kBarRecoEffCorMult[isp][ietaA][ipt] * p1kBarRecoEffCor[isp][ietaC][ipt];

                float covFT0ATru = p1kBarFt0A * p1kBarTru[isp][ietaC][ipt];
                float covFT0AReco = p1kBarFt0A * p1kBarReco[isp][ietaC][ipt];
                float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[isp][ietaC][ipt];

                float covFT0CTru = p1kBarFt0C * p1kBarTru[isp][ietaA][ipt];
                float covFT0CReco = p1kBarFt0C * p1kBarReco[isp][ietaA][ipt];
                float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[isp][ietaA][ipt];

                if (isp == numKInclusive) { // Inclusive
                  if (ipt == 0) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt0_GapSum2D"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt0_GapSum2D"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_GapSum2D"), cent, gap, sum, c2SubRecoEffCor);
                    }

                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 2) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt1_GapSum2D"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt1_GapSum2D"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_GapSum2D"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 1) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt2_GapSum2D"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt2_GapSum2D"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_GapSum2D"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CRecoEffCor);
                  }

                } else if (isp == numKPion) { // Pion
                  if (ipt == 0) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt0_GapSum2D_Pi"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt0_GapSum2D_Pi"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_GapSum2D_Pi"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 2) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt1_GapSum2D_Pi"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt1_GapSum2D_Pi"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_GapSum2D_Pi"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 1) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt2_GapSum2D_Pi"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt2_GapSum2D_Pi"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_GapSum2D_Pi"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0CRecoEffCor);
                  }
                } else if (isp == numKKaon) { // Kaon
                  if (ipt == 0) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt0_GapSum2D_Ka"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt0_GapSum2D_Ka"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_GapSum2D_Ka"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 2) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt1_GapSum2D_Ka"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt1_GapSum2D_Ka"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_GapSum2D_Ka"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 1) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt2_GapSum2D_Ka"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt2_GapSum2D_Ka"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_GapSum2D_Ka"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0CRecoEffCor);
                  }
                } else if (isp == numKProton) { // Proton
                  if (ipt == 0) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt0_GapSum2D_Pr"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt0_GapSum2D_Pr"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_GapSum2D_Pr"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 2) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt1_GapSum2D_Pr"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt1_GapSum2D_Pr"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_GapSum2D_Pr"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CRecoEffCor);

                  } else if (ipt == KNpT - 1) {
                    if (std::isfinite(c2SubTru)) {
                      histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubTru);
                      histos.fill(HIST("MCGen/Prof_ipt2_GapSum2D_Pr"), cent, gap, sum, c2SubTru);
                    }
                    if (std::isfinite(c2SubReco)) {
                      histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubReco);
                      histos.fill(HIST("MCReco/Prof_ipt2_GapSum2D_Pr"), cent, gap, sum, c2SubReco);
                    }
                    if (std::isfinite(c2SubRecoEffCor)) {
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2SubRecoEffCor);
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_GapSum2D_Pr"), cent, gap, sum, c2SubRecoEffCor);
                    }
                    if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covTru);
                    if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covReco);
                    if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covRecoEffCor);

                    if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ATru);
                    if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0AReco);
                    if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                    if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CTru);
                    if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CReco);
                    if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0CRecoEffCor);
                  }
                }
              }
            }
          }
        }
      } // colSlice
    } // mcColl
    LOGF(info, "FINISHED RUNNING processMCFluc");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt fluc", cfgRunMCFluc);

  void processGetDataFlat(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, AodTracksSel const& tracks)
  {
    histos.fill(HIST("hVtxZ"), coll.posZ());
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);

    int ntrk = 0;
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

      // Count tracks in the primary eta acceptance
      if (eta > etaLw[0] && eta < etaUp[0])
        ntrk++;

      // Define species array (0: Inclusive, 1: Pion, 2: Kaon, 3: Proton)
      bool isSpecies[KNsp] = {true, selectionPion(track), selectionKaon(track), selectionProton(track)};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;

        // Fetch efficiency specifically for this particle species
        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
        float w = (1.0 - fake) / eff;

        if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
          continue;

        // Unrolled THnSparse / QA Fills
        if (isp == numKInclusive) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKPion) { // Pion
          histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pi"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKKaon) { // Kaon
          histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Ka"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKProton) { // Proton
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
        }
      }
    }

    histos.fill(HIST("hCentnTrk"), cent, ntrk);
    histos.fill(HIST("hCentnTrkPV"), cent, coll.multNTracksPV());

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
    // Expanded to 4 species (isp = 0: Incl, 1: Pi, 2: Ka, 3: Pr)
    double sumWi[KNsp][KNEta][KNpT]{}, sumWipti[KNsp][KNEta][KNpT]{};

    if (!isEventSelected(coll))
      return;

    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);

    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;

      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      auto sign = track.sign();

      if (p < KFloatEpsilon)
        continue;

      histos.fill(HIST("hP"), p);
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      // Define species array
      bool isSpecies[KNsp] = {true, selectionPion(track), selectionKaon(track), selectionProton(track)};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;

        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
        float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
        float w = flatWeight * (1.0 - fake) / eff;

        if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
          continue;

        // Unrolled THnSparse / QA Fills
        if (isp == numKInclusive) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKPion) { // Pion
          histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pi"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKKaon) { // Kaon
          histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Ka"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKProton) { // Proton
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
        }

        // Accumulate sum
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt])
              continue;
            sumWi[isp][ieta][ipt] += w;
            sumWipti[isp][ieta][ipt] += w * pt;
          }
        }
      }
    }

    // Full Event Means
    for (int isp = 0; isp < KNsp; ++isp) {
      if (isp == numKInclusive) {
        histos.fill(HIST("Prof_Cent_Nchrec"), cent, sumWi[0][0][0]);
        histos.fill(HIST("Prof_Mult_Nchrec"), coll.multNTracksPV(), sumWi[0][0][0]);
        if (sumWi[0][0][0] > 1.0f)
          histos.fill(HIST("Prof_Cent_MeanpT"), cent, sumWipti[0][0][0] / sumWi[0][0][0]);
      } else if (isp == numKPion) {
        histos.fill(HIST("Prof_Cent_Nchrec_Pi"), cent, sumWi[1][0][0]);
        histos.fill(HIST("Prof_Mult_Nchrec_Pi"), coll.multNTracksPV(), sumWi[1][0][0]);

        if (sumWi[1][0][0] > 1.0f)
          histos.fill(HIST("Prof_Cent_MeanpT_Pi"), cent, sumWipti[1][0][0] / sumWi[1][0][0]);
      } else if (isp == numKKaon) {
        histos.fill(HIST("Prof_Cent_Nchrec_Ka"), cent, sumWi[2][0][0]);
        histos.fill(HIST("Prof_Mult_Nchrec_Ka"), coll.multNTracksPV(), sumWi[2][0][0]);

        if (sumWi[2][0][0] > 1.0f)
          histos.fill(HIST("Prof_Cent_MeanpT_Ka"), cent, sumWipti[2][0][0] / sumWi[2][0][0]);
      } else if (isp == numKProton) {
        histos.fill(HIST("Prof_Cent_Nchrec_Pr"), cent, sumWi[3][0][0]);
        histos.fill(HIST("Prof_Mult_Nchrec_Pr"), coll.multNTracksPV(), sumWi[3][0][0]);

        if (sumWi[3][0][0] > 1.0f)
          histos.fill(HIST("Prof_Cent_MeanpT_Pr"), cent, sumWipti[3][0][0] / sumWi[3][0][0]);
      }
    }

    // Kinematic Bin Means (1D and 2D Sub-event)
    for (int ietaA = 0; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 0; ietaC < KNEta; ++ietaC) {
        for (int ipt = 0; ipt < KNpT; ++ipt) {
          for (int isp = 0; isp < KNsp; ++isp) {

            // --- 2D Sub-Event Calculations ---
            double wCorrAB = sumWi[isp][ietaA][ipt] + sumWi[isp][ietaC][ipt];
            if (wCorrAB > 0) {
              float mptsub = (sumWipti[isp][ietaA][ipt] + sumWipti[isp][ietaC][ipt]) / wCorrAB;
              if (isp == numKInclusive) {
                if (ipt == 0)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 2)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 1)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2"), cent, ietaA, ietaC, mptsub);
              } else if (isp == numKPion) {
                if (ipt == 0)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Pi"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 2)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Pi"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 1)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Pi"), cent, ietaA, ietaC, mptsub);
              } else if (isp == numKKaon) {
                if (ipt == 0)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Ka"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 2)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Ka"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 1)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Ka"), cent, ietaA, ietaC, mptsub);
              } else if (isp == numKProton) {
                if (ipt == 0)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt0_Pr"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 2)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt1_Pr"), cent, ietaA, ietaC, mptsub);
                if (ipt == KNpT - 1)
                  histos.fill(HIST("Prof2D_MeanpT_Sub_ipt2_Pr"), cent, ietaA, ietaC, mptsub);
              }
            }

            // --- 1D Individual Bin Calculations (Only do when A == B to avoid overfilling) ---
            if (ietaA == ietaC) {
              double mpt = sumWipti[isp][ietaA][ipt] / sumWi[isp][ietaA][ipt];
              if (sumWi[isp][ietaA][ipt] >= 1.0f && std::isfinite(mpt)) {
                if (isp == numKInclusive) {
                  histos.fill(HIST("pmean_nch_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_nch_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, sumWi[0][ietaA][ipt]);
                  histos.fill(HIST("pmean_cent_etabin_ptbin"), cent, ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_cent_etabin_ptbin"), cent, ietaA, ipt, sumWi[0][ietaA][ipt]);
                } else if (isp == numKPion) {
                  histos.fill(HIST("pmean_nch_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_nch_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, sumWi[1][ietaA][ipt]);
                  histos.fill(HIST("pmean_cent_etabin_ptbin_Pi"), cent, ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_cent_etabin_ptbin_Pi"), cent, ietaA, ipt, sumWi[1][ietaA][ipt]);
                } else if (isp == numKKaon) {
                  histos.fill(HIST("pmean_nch_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_nch_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, sumWi[2][ietaA][ipt]);
                  histos.fill(HIST("pmean_cent_etabin_ptbin_Ka"), cent, ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_cent_etabin_ptbin_Ka"), cent, ietaA, ipt, sumWi[2][ietaA][ipt]);
                } else if (isp == numKProton) {
                  histos.fill(HIST("pmean_nch_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_nch_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, sumWi[3][ietaA][ipt]);
                  histos.fill(HIST("pmean_cent_etabin_ptbin_Pr"), cent, ietaA, ipt, mpt);
                  histos.fill(HIST("pmeanMult_cent_etabin_ptbin_Pr"), cent, ietaA, ipt, sumWi[3][ietaA][ipt]);
                }
              }
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
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, iCh, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, iCh, eta, ampl);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        auto eta = getEtaFT0(globalId, 1);
        amplFT0C += ampl;
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, iCh, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, iCh, eta, ampl);
      }
    }

    histos.fill(HIST("pmeanFT0A_multpv"), coll.multNTracksPV(), amplFT0A);
    histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
    histos.fill(HIST("pmeanFT0C_multpv"), coll.multNTracksPV(), amplFT0C);
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

    // 1. Safety Check: Step 2 Mean Maps
    for (int isp = 0; isp < KNsp; ++isp) {
      if (!pmeanNchEtabinPtbinStep2[isp] || !pmeanMultNchEtabinPtbinStep2[isp]) {
        LOGF(warning, "Data fluc: Mean pT or Mult map missing for species index %d", isp);
        return;
      }
    }

    // 2. Safety Check: Correction Maps (Looping over Inclusive, Pi, Ka, Pr)
    for (int isp = 0; isp < KNsp; ++isp) {
      auto pid = static_cast<PID>(isp);
      if (!hEff[pid] || !hFake[pid] || !hFlatWeight[pid]) {
        LOGF(warning, "Data fluc: Correction maps (Eff, Fake, or Flat) are null for species index %d", isp);
        return;
      }
    }

    // Expanded arrays to handle KNsp species (0: Incl, 1: Pi, 2: Ka, 3: Pr)
    double sumpmwk[KNsp][KNEta][KNpT][KIntM][KIntK]{};
    double sumwk[KNsp][KNEta][KNpT][KIntK]{};

    double mean[KNsp][KNEta][KNpT]{}, c2[KNsp][KNEta][KNpT]{};
    double p1kBar[KNsp][KNEta][KNpT]{};
    double meanMult[KNsp][KNEta][KNpT]{}, p1kBarMult[KNsp][KNEta][KNpT]{};

    float vz = coll.posZ();

    // --- 1. Track Loop: Accumulate sum ---
    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;

      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      auto sign = track.sign();

      if (p < KFloatEpsilon)
        continue;

      bool isSpecies[KNsp] = {true, selectionPion(track), selectionKaon(track), selectionProton(track)};

      for (int isp = 0; isp < KNsp; ++isp) {
        if (!isSpecies[isp])
          continue;

        float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
        float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
        float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
        float w = flatWeight * (1.0 - fake) / eff;

        if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
          continue;

        // QA Fills
        if (isp == numKInclusive) {
          histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKPion) { // Pion
          histos.fill(HIST("hEtaPhiReco_Pi"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pi"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pi"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKKaon) { // Kaon
          histos.fill(HIST("hEtaPhiReco_Ka"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Ka"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Ka"), vz, sign, pt, eta, phi, w);
        } else if (isp == numKProton) { // Proton
          histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
          histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
          histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
        }

        // Kinematic Bin sum
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt])
              continue;
            for (int k = 0; k < KIntK; ++k) {
              for (int m = 0; m < KIntM; ++m) {
                sumpmwk[isp][ieta][ipt][m][k] += std::pow(w, k) * std::pow(pt, m);
              }
              sumwk[isp][ieta][ipt][k] += std::pow(w, k);
            }
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
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        amplFT0C += ampl;
      }
    }
    double p1kBarFt0A = amplFT0A - pmeanFT0A_multpvStep2->GetBinContent(pmeanFT0A_multpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));
    double p1kBarFt0C = amplFT0C - pmeanFT0C_multpvStep2->GetBinContent(pmeanFT0C_multpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));

    // --- 2. Step 2 Means and 1D Fluc Variables ---
    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int ipt = 0; ipt < KNpT; ++ipt) {

        // Use [0] to safely grab the X-axis from the array!
        const int ibx = pmeanNchEtabinPtbinStep2[0]->GetXaxis()->FindBin(coll.multNTracksPV());
        const int iby = ieta + 1;
        const int ibz = ipt + 1;

        for (int isp = 0; isp < KNsp; ++isp) {
          // Dynamically fetch from the array
          float mmpt = pmeanNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);
          float mmMult = pmeanMultNchEtabinPtbinStep2[isp]->GetBinContent(ibx, iby, ibz);

          mean[isp][ieta][ipt] = sumpmwk[isp][ieta][ipt][1][1] / sumwk[isp][ieta][ipt][1];
          meanMult[isp][ieta][ipt] = sumwk[isp][ieta][ipt][1];

          if (std::isfinite(mmpt)) {
            std::tie(mean[isp][ieta][ipt], c2[isp][ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwk[isp][ieta][ipt], sumwk[isp][ieta][ipt], mmpt);
            p1kBar[isp][ieta][ipt] = mean[isp][ieta][ipt] - mmpt;
          }
          p1kBarMult[isp][ieta][ipt] = meanMult[isp][ieta][ipt] - mmMult;
        }
      }
    }

    // --- 3. Fill 1D Profiles ---
    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        for (int isp = 0; isp < KNsp; ++isp) {
          if (isp == numKInclusive) {
            if (std::isfinite(mean[0][ieta][ipt])) {
              histos.fill(HIST("Prof_MeanpT_Cent_etabin_ptbin"), cent, ieta, ipt, mean[0][ieta][ipt]);
              histos.fill(HIST("Prof_MeanpT_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, mean[0][ieta][ipt]);
            }
            if (std::isfinite(c2[0][ieta][ipt])) {
              histos.fill(HIST("Prof_C2_Cent_etabin_ptbin"), cent, ieta, ipt, c2[0][ieta][ipt]);
              histos.fill(HIST("Prof_C2_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, c2[0][ieta][ipt]);
            }
          } else if (isp == numKPion) { // Pi
            if (std::isfinite(mean[1][ieta][ipt])) {
              histos.fill(HIST("Prof_MeanpT_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, mean[1][ieta][ipt]);
              histos.fill(HIST("Prof_MeanpT_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ieta, ipt, mean[1][ieta][ipt]);
            }
            if (std::isfinite(c2[1][ieta][ipt])) {
              histos.fill(HIST("Prof_C2_Cent_etabin_ptbin_Pi"), cent, ieta, ipt, c2[1][ieta][ipt]);
              histos.fill(HIST("Prof_C2_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ieta, ipt, c2[1][ieta][ipt]);
            }
          } else if (isp == numKKaon) { // Ka
            if (std::isfinite(mean[2][ieta][ipt])) {
              histos.fill(HIST("Prof_MeanpT_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, mean[2][ieta][ipt]);
              histos.fill(HIST("Prof_MeanpT_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ieta, ipt, mean[2][ieta][ipt]);
            }
            if (std::isfinite(c2[2][ieta][ipt])) {
              histos.fill(HIST("Prof_C2_Cent_etabin_ptbin_Ka"), cent, ieta, ipt, c2[2][ieta][ipt]);
              histos.fill(HIST("Prof_C2_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ieta, ipt, c2[2][ieta][ipt]);
            }
          } else if (isp == numKProton) { // Pr
            if (std::isfinite(mean[3][ieta][ipt])) {
              histos.fill(HIST("Prof_MeanpT_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, mean[3][ieta][ipt]);
              histos.fill(HIST("Prof_MeanpT_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ieta, ipt, mean[3][ieta][ipt]);
            }
            if (std::isfinite(c2[3][ieta][ipt])) {
              histos.fill(HIST("Prof_C2_Cent_etabin_ptbin_Pr"), cent, ieta, ipt, c2[3][ieta][ipt]);
              histos.fill(HIST("Prof_C2_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ieta, ipt, c2[3][ieta][ipt]);
            }
          }
        }
      }
    }

    for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
      int ietaC = KNEta - ietaA;
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        for (int isp = 0; isp < KNsp; ++isp) {
          float c2Sub = p1kBar[isp][ietaA][ipt] * p1kBar[isp][ietaC][ipt];
          float covAC = p1kBarMult[isp][ietaA][ipt] * p1kBar[isp][ietaC][ipt];
          float covCA = p1kBar[isp][ietaA][ipt] * p1kBarMult[isp][ietaC][ipt];

          float covFT0A = p1kBarFt0A * p1kBar[isp][ietaC][ipt];
          float covFT0C = p1kBarFt0C * p1kBar[isp][ietaA][ipt];

          if (isp == numKInclusive) {
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub_Cent_etabin_ptbin"), cent, ietaA, ipt, c2Sub);
              histos.fill(HIST("Prof_C2Sub_Mult_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, c2Sub);
            }
            if (std::isfinite(covAC)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin"), cent, ietaA, ipt, covAC);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, covAC);
            }
            if (std::isfinite(covCA)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin"), cent, ietaA, ipt, covCA);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, covCA);
            }

            if (std::isfinite(covFT0A)) {
              histos.fill(HIST("Prof_CovFT0A_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0A);
              histos.fill(HIST("Prof_CovFT0A_Mult_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, covFT0A);
            }
            if (std::isfinite(covFT0C)) {
              histos.fill(HIST("Prof_CovFT0C_Cent_etabin_ptbin"), cent, ietaA, ipt, covFT0C);
              histos.fill(HIST("Prof_CovFT0C_Mult_etabin_ptbin"), coll.multNTracksPV(), ietaA, ipt, covFT0C);
            }

          } else if (isp == numKPion) { // Pi
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, c2Sub);
              histos.fill(HIST("Prof_C2Sub_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, c2Sub);
            }
            if (std::isfinite(covAC)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covAC);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, covAC);
            }
            if (std::isfinite(covCA)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covCA);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, covCA);
            }
            if (std::isfinite(covFT0A)) {
              histos.fill(HIST("Prof_CovFT0A_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0A);
              histos.fill(HIST("Prof_CovFT0A_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, covFT0A);
            }
            if (std::isfinite(covFT0C)) {
              histos.fill(HIST("Prof_CovFT0C_Cent_etabin_ptbin_Pi"), cent, ietaA, ipt, covFT0C);
              histos.fill(HIST("Prof_CovFT0C_Mult_etabin_ptbin_Pi"), coll.multNTracksPV(), ietaA, ipt, covFT0C);
            }

          } else if (isp == numKKaon) { // Ka
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, c2Sub);
              histos.fill(HIST("Prof_C2Sub_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, c2Sub);
            }
            if (std::isfinite(covAC)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covAC);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, covAC);
            }
            if (std::isfinite(covCA)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covCA);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, covCA);
            }
            if (std::isfinite(covFT0A)) {
              histos.fill(HIST("Prof_CovFT0A_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0A);
              histos.fill(HIST("Prof_CovFT0A_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, covFT0A);
            }
            if (std::isfinite(covFT0C)) {
              histos.fill(HIST("Prof_CovFT0C_Cent_etabin_ptbin_Ka"), cent, ietaA, ipt, covFT0C);
              histos.fill(HIST("Prof_CovFT0C_Mult_etabin_ptbin_Ka"), coll.multNTracksPV(), ietaA, ipt, covFT0C);
            }

          } else if (isp == numKProton) { // Pr
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("Prof_C2Sub_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, c2Sub);
              histos.fill(HIST("Prof_C2Sub_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, c2Sub);
            }
            if (std::isfinite(covAC)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covAC);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, covAC);
            }
            if (std::isfinite(covCA)) {
              histos.fill(HIST("Prof_Cov_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covCA);
              histos.fill(HIST("Prof_Cov_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, covCA);
            }
            if (std::isfinite(covFT0A)) {
              histos.fill(HIST("Prof_CovFT0A_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0A);
              histos.fill(HIST("Prof_CovFT0A_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, covFT0A);
            }
            if (std::isfinite(covFT0C)) {
              histos.fill(HIST("Prof_CovFT0C_Cent_etabin_ptbin_Pr"), cent, ietaA, ipt, covFT0C);
              histos.fill(HIST("Prof_CovFT0C_Mult_etabin_ptbin_Pr"), coll.multNTracksPV(), ietaA, ipt, covFT0C);
            }
          }
        }
      }
    }

    // --- 5. Full 2D Covariances & GapSum Profiles ---
    for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 1; ietaC < KNEta; ++ietaC) {

        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        float gap = etaValA - etaValB;
        float sum = (etaValA + etaValB) / 2.0f;

        for (int ipt = 0; ipt < KNpT; ++ipt) {
          for (int isp = 0; isp < KNsp; ++isp) {

            float c2Sub = p1kBar[isp][ietaA][ipt] * p1kBar[isp][ietaC][ipt];
            float cov = p1kBarMult[isp][ietaA][ipt] * p1kBar[isp][ietaC][ipt];
            float covFT0A = p1kBarFt0A * p1kBar[isp][ietaC][ipt];
            float covFT0C = p1kBarFt0C * p1kBar[isp][ietaA][ipt];

            if (isp == numKInclusive) { // Inclusive
              if (ipt == 0) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt0_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt0_GapSum2D"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt0_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt0_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt0_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 2) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt1_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt1_GapSum2D"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt1_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt1_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt1_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 1) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt2_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt2_GapSum2D"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt2_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt2_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt2_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0C);
              }
            } else if (isp == numKPion) { // Pi
              if (ipt == 0) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt0_GapSum2D_Pi"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt0_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 2) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt1_GapSum2D_Pi"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt1_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 1) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt2_GapSum2D_Pi"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt2_Cov2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pi"), cent, etaValA, etaValB, covFT0C);
              }
            } else if (isp == numKKaon) { // Ka
              if (ipt == 0) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt0_GapSum2D_Ka"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt0_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 2) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt1_GapSum2D_Ka"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt1_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 1) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt2_GapSum2D_Ka"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt2_Cov2D_Cent_etaA_etaC_Ka"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0C);
              }
            } else if (isp == numKProton) { // Pr
              if (ipt == 0) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt0_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt0_GapSum2D_Pr"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt0_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt0_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt0_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 2) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt1_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt1_GapSum2D_Pr"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt1_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt1_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt1_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0C);
              } else if (ipt == KNpT - 1) {
                if (std::isfinite(c2Sub)) {
                  histos.fill(HIST("Prof_ipt2_C2Sub2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, c2Sub);
                  histos.fill(HIST("Prof_ipt2_GapSum2D_Pr"), cent, gap, sum, c2Sub);
                }
                if (std::isfinite(cov))
                  histos.fill(HIST("Prof_ipt2_Cov2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, cov);
                if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_ipt2_CovFT0A2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0A);
                if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_ipt2_CovFT0C2D_Cent_etaA_etaC_Pr"), cent, etaValA, etaValB, covFT0C);
              }
            }
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
