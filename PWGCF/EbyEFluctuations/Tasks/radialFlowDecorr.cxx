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
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
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

  static constexpr int KIntM = 3;
  static constexpr int KIntK = 3;

  static constexpr int KNEta = 17;
  static constexpr int KNpT = 3;

  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr int KPiPlus = 211;
  static constexpr int KKPlus = 321;
  static constexpr int KProton = 2212;

  static constexpr float KCentTestMin = 10.f;
  static constexpr float KCentTestMaxLo = 60.f;
  static constexpr float KCentTestMaxHi = 70.f;
  static constexpr float KCentCovCut = 1.0f;
  static constexpr float KBinOffset = 0.5f;
  static constexpr int KNchPbMax = 3600;
  static constexpr int KNchOMax = 600;

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
  enum PID { kInclusive = 0,
             kCombinedPID,
             kNumPID };
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
  const std::vector<std::string> pidSuffix = {"", "_PID"};

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
  Configurable<int> cfgSys{"cfgSys", 2, "Efficiency to be used for which system? 1-->PbPb, 2-->OO, 3-->pPb, 4-->pp"};
  Configurable<bool> cfgFlat{"cfgFlat", true, "Whether to use flattening weights or not"};

  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 3.0, 5.0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "centrality axis (percentile)"};
  const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};

  AxisSpec nChAxis{1, 0., 1., "Nch", "Nch"};
  AxisSpec nChAxis2{1, 0., 1., "Nch", "Nch"};

  const AxisSpec vzAxis{5, -12.5, 12.5,
                        "Vz"
                        "Vz"};
  const AxisSpec chgAxis{3, -1.5, 1.5};
  ConfigurableAxis cfgpTAxis{"cfgpTAxis", {0.0, 0.2, 0.5, 1, 3, 5, 7.5, 10}, "pT axis for flattening"};
  const AxisSpec pTAxis{cfgpTAxis, "pT"};

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

  std::array<TH3F*, kNumPID> hEff{};
  std::array<TH3F*, kNumPID> hFake{};
  std::array<THnSparseF*, kNumPID> hFlatWeight{};

  TProfile3D* pmeanTruNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanRecoNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanRecoMatchedNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanRecoEffcorrNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanRecoMatchedEffcorrNchEtabinPtbinStep2 = nullptr;

  TProfile3D* pmeanEtTruNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanEtRecoNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanEtRecoMatchedNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanEtRecoEffcorrNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanEtRecoMatchedEffcorrNchEtabinPtbinStep2 = nullptr;

  TProfile3D* pmeanNchEtabinPtbinStep2 = nullptr;
  TProfile3D* pmeanEtNchEtabinPtbinStep2 = nullptr;

  // Helper to calculate all three combined PID sigmas at once
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
    if (cfgEvSelkNoITSROFrameBorder && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder))
      return false;
    if (cfgEvSelkNoTimeFrameBorder && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))
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
    // if (dpt::isStrangeBaryonPDG(particle.pdgCode())) return false;
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

  float getEfficiency(float mult, float pt, float eta, PID pidType, int effidx) const
  {
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
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis}});
    histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis2}, {nChAxis2}});
    histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis2}, {centAxis}});
    histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsP, KPMin, KPMax}});
    histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("hEta", ";#eta", kTH1F, {{KNbinsEta, KEtaMin, KEtaMax}});
    histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});

    histos.add("hEtaPhiReco", "hEtaPhiReco", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1) / 2, KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiReco_PID", "hEtaPhiReco_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1) / 2, KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
  }
  void declareMCCommonHists()
  {

    histos.add("ptResolution", ";p_{T}^{MC};p_{T}^{MC}-p_{T}^{reco}", kTH2F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, -KPtResMax, KPtResMax}});
    histos.add("ptTruthReco", ";p_{T}^{MC};p_{T}^{reco}", kTH2F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("etaResolution", ";#eta^{MC};#eta^{MC}-#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsPtRes, -KEtaResMax, KEtaResMax}});
    histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsPtRes, -KEtaFineMax, KEtaFineMax}, {KNbinsPtRes, -KEtaFineMax, KEtaFineMax}});

    histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
    histos.add("vzResolution", ";Vz^{MC};Vz^{MC}-Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, -KVzResMax, KVzResMax}});

    histos.add("h3_AllPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoMatchedToPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Secondary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Fake", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_AllReco", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("h3_AllPrimary_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoMatchedToPrimary_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Secondary_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Fake_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_AllReco_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("h_AllPrimary", ";p_{T}", kTH1F, {{KNbinsP, cfgPtMin, cfgPtMax}});
    histos.add("h_RecoMatchedToPrimary", ";p_{T}", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("h_RecoUnMatchedToPrimary", ";p_{T}", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("h_AllReco", ";p_{T}", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("h_AllRecoEffCorr", ";p_{T}", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});

    histos.add("hReco_ParticleWeight", ";cent;p_{T};#eta", kTH3F, {{centAxis}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, -KEtaFineMax, KEtaFineMax}});
    histos.add("hTruth_ParticleWeight", ";cent;p_{T};#eta", kTH3F, {{centAxis}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsPtRes, -KEtaFineMax, KEtaFineMax}});

    histos.add("hDCAxy_Unmatched", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_Unmatched", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_NotPrimary", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_NotPrimary", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_RecoMatched", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_RecoMatched", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAxy_Reco", ";DCA_{xy} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
    histos.add("hDCAz_Reco", ";DCA_{z} (cm)", kTH1F, {{KNbinsDca, -KDcaMax, KDcaMax}});
  }

  void declareMCMeanHists()
  {
    histos.add("Eff_cent", ";cent;#epsilon", kTProfile, {centAxis});
    histos.add("Fake_cent", ";cent;f_{fake}", kTProfile, {centAxis});
    histos.add("wgt_cent", ";cent;w", kTProfile, {centAxis});
    histos.add("Eff_Ntrk", ";N_{PV};#epsilon", kTProfile, {nChAxis2});
    histos.add("Fake_Ntrk", ";N_{PV};f_{fake}", kTProfile, {nChAxis2});
    histos.add("wgt_Ntrk", ";N_{PV};w", kTProfile, {nChAxis2});
    histos.add("Eff_pT", ";p_{T};#epsilon", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("Fake_pT", ";p_{T};f_{fake}", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
    histos.add("wgt_pT", ";p_{T};w", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Eff_eta", ";#eta;#epsilon", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("Fake_eta", ";#eta;f_{fake}", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("wgt_eta", ";#eta;w", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    // MC mean profiles (pT & Et) for various selections
    histos.add("MCGen/Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanpT_Mult", ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});
    histos.add<TProfile3D>("pmeanTruNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoMatchedNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoMatchedEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add("MCGen/Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanEt_Mult", ";N_{PV};#LT E_{T}#GT", kTProfile, {nChAxis});
    histos.add<TProfile3D>("pmeanEtTruNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoMatchedNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoMatchedEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

    histos.addClone("MCGen/", "MCReco/");
    histos.addClone("MCGen/", "MCRecoMatched/");
    histos.addClone("MCGen/", "MCRecoEffCorr/");
    histos.addClone("MCGen/", "MCRecoMatchedEffCorr/");
  }
  void declareMCFlucHists()
  {
    // pT cumulants
    histos.add("MCGen/Prof_C2_Cent", ";cent;C_{2}", kTProfile, {centAxis});
    histos.add("MCGen/Prof_C2_Mult", ";N_{PV};C_{2}", kTProfile, {nChAxis});
    histos.add<TProfile3D>("MCGen/Prof_C2Sub_Mult_etabin_ptbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis2}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("MCGen/Prof_ipt0_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("MCGen/Prof_ipt1_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("MCGen/Prof_ipt2_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt0_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt0_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt1_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt1_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt2_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt2_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add("MCGen/Prof_C2Et_Cent", ";cent;C_{2}^{E_{T}}", kTProfile, {centAxis});
    histos.add("MCGen/Prof_C2Et_Mult", ";N_{PV};C_{2}^{E_{T}}", kTProfile, {nChAxis});
    histos.add<TProfile3D>("MCGen/Prof_C2EtSub_Mult_etabin_ptbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis2}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("MCGen/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("MCGen/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("MCGen/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt0_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt0_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt1_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt1_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("MCGen/Prof_ipt2_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("MCGen/Prof_ipt2_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});

    histos.add("MCGen/Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanpT_Mult", ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});

    histos.add<TProfile3D>("pmeanTruNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoMatchedNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanRecoMatchedEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add("MCGen/Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
    histos.add("MCGen/Prof_MeanEt_Mult", ";N_{PV};#LT E_{T}#GT", kTProfile, {nChAxis});
    histos.add<TProfile3D>("pmeanEtTruNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoMatchedNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add<TProfile3D>("pmeanEtRecoMatchedEffcorrNchEtabinPtbin", ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
  }
  void declareDataMeanHists()
  {
    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

    histos.add("Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
    histos.add("Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
    histos.add<TProfile3D>("pmean_nch_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    // Et
    histos.add("Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
    histos.add<TProfile3D>("pmeanEt_nch_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
  }

  void declareDataGetFlatHists()
  {
    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
  }
  void declareMCGetFlatHists()
  {
    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
  }

  void declareDataFlucHists()
  {
    histos.add("hEtaPhiRecoWtd", "hEtaPhiRecoWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd", "hEtaPhiRecoEffWtd", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoWtd_PID", "hEtaPhiRecoWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});
    histos.add("hEtaPhiRecoEffWtd_PID", "hEtaPhiRecoEffWtd_PID", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {(KNEta - 1), KEtaAxisMin, KEtaAxisMax}, {KNbinsPhiFine, KPhiMin, TwoPI}});

    histos.add<TProfile3D>("pmean_nch_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});
    histos.add("Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
    histos.add<TProfile3D>("pmeanEt_nch_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add("Prof_C2_Cent", ";cent;C_{2}", kTProfile, {centAxis});
    histos.add("Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
    histos.add("Prof_MeanpT_Mult", ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});
    histos.add<TProfile3D>("Prof_C2_Mult_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add<TProfile2D>("Prof_ipt0_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt0_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("Prof_ipt1_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt1_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("Prof_ipt2_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt2_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});

    histos.add("Prof_C2Et_Cent", ";cent;C_{2}^{E_{T}}", kTProfile, {centAxis});
    histos.add("Prof_MeanEt_Mult", ";N_{PV};#LT E_{T}#GT", kTProfile, {nChAxis});
    histos.add<TProfile3D>("Prof_C2Et_Mult_etabin_ptbin", ";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNpT + 1, -KBinOffset, KNpT + KBinOffset}});

    histos.add<TProfile2D>("Prof_ipt0_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt0_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("Prof_ipt1_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt1_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});
    histos.add<TProfile2D>("Prof_ipt2_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta - 1) / 2, 0, KEtaAxisMax}});
    histos.add("Prof_ipt2_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta - 1) / 2, 0., KEtaAxisMax}});

    histos.add<TProfile3D>("Prof_ipt0_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("Prof_ipt1_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("Prof_ipt2_C2Sub2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("Prof_ipt0_C2SubEt2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("Prof_ipt1_C2SubEt2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
    histos.add<TProfile3D>("Prof_ipt2_C2SubEt2D_Mult_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}, {KNEta - 1, KEtaAxisMin, KEtaAxisMax}});
  }

  THnSparseF* buildWeightMapFromRaw(THnSparseF* hRaw, const char* mapName)
  {
    if (!hRaw) {
      LOGF(error, "Raw eta-phi map for '%s' is null; no flattening will be applied.", mapName);
      return nullptr;
    }
    auto hWMap = reinterpret_cast<THnSparseF*>(hRaw->Clone(mapName));
    hWMap->SetTitle(Form("Flattening Weight Map %s (w_{#phi} = <N_{#phi}> / N_{#phi})", mapName));
    // hWMap->SetDirectory(nullptr);
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
    if (cfgSys == 1) {
      nChAxis = {KNchPbMax / 3, KBinOffset, KNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {KNchPbMax / 10, KBinOffset, KNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else {
      nChAxis = {KNchOMax, KBinOffset, KNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {KNchOMax / 6, KBinOffset, KNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    }

    ccdb->setURL(cfgCCDBurl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

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
    }
    if (cfgRunMCFluc) {
      declareMCFlucHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoMatched/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
      histos.addClone("MCGen/", "MCRecoMatchedEffCorr/");
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

    if (!cfgRunGetEff) {
      TList* lst = ccdb->getForTimeStamp<TList>(pathEff, now);

      if (!lst) {
        LOGF(fatal, "Efficiency maps required but CCDB list is null at %s!", pathEff.c_str());
      }

      LOGF(info, "Loading Eff/Fake maps from TList...");

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
          LOGF(error, "Missing CCDB objects for efficiency. Checked in list: %s, %s", hEffNumName.c_str(), hEffDenName.c_str());
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

      loadEffFakeForPID(kInclusive);
      loadEffFakeForPID(kCombinedPID);

      // --- 1. Load Data Flattening Maps (if DataMean or DataFluc) ---
      if (cfgRunDataMean || cfgRunDataFluc) {
        LOGF(info, "Data Run: Loading flattening maps from CCDB path: %s", pathDataFlat.c_str());

        TList* lstDataFlat = ccdb->getForTimeStamp<TList>(pathDataFlat, now);

        if (lstDataFlat) {
          auto* hRawIncl = reinterpret_cast<THnSparseF*>(lstDataFlat->FindObject("hEtaPhiRecoEffWtd"));
          if (hRawIncl) {
            hFlatWeight[kInclusive] = buildWeightMapFromRaw(hRawIncl, "hFlatWeight");
          } else {
            LOGF(error, "Data flattening 'hEtaPhiRecoEffWtd' not found in list from %s", pathDataFlat.c_str());
          }

          auto* hRawPID = reinterpret_cast<THnSparseF*>(lstDataFlat->FindObject("hEtaPhiRecoEffWtd_PID"));
          if (hRawPID) {
            hFlatWeight[kCombinedPID] = buildWeightMapFromRaw(hRawPID, "hFlatWeight_PID");
          } else {
            LOGF(error, "Data flattening 'hEtaPhiRecoEffWtd_PID' not found in list from %s", pathDataFlat.c_str());
          }
        } else {
          LOGF(error, "Could not retrieve TList for Data Flattening from: %s", pathDataFlat.c_str());
        }
      }

      // --- 2. Load MC Flattening Maps (if MCMean or MCFluc) ---
      if (cfgRunMCMean || cfgRunMCFluc) {
        LOGF(info, "MC Run: Loading flattening maps from MC Flat list (%s)...", pathMCFlat.c_str());

        TList* lstMCFlat = ccdb->getForTimeStamp<TList>(pathMCFlat, now);
        if (lstMCFlat) {
          auto loadFlatForPID = [&](PID pidType) {
            std::string suffix = pidSuffix[pidType];
            std::string hFlatSrcName = "hEtaPhiRecoEffWtd" + suffix;

            auto* hRaw = reinterpret_cast<THnSparseF*>(lstMCFlat->FindObject(hFlatSrcName.c_str()));

            if (hRaw) {
              hFlatWeight[pidType] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
            } else {
              LOGF(warning, "MC flattening source '%s' not found in list; skipping this PID.", hFlatSrcName.c_str());
            }
          };

          loadFlatForPID(kInclusive);
          loadFlatForPID(kCombinedPID);
        } else {
          LOGF(error, "Could not retrieve TList for MC Flattening from: %s", pathMCFlat.c_str());
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

    if (cfgRunMCFluc) {
      LOGF(info, "Loading MC Mean profiles from CCDB path: %s", pathMCMean.c_str());
      TList* lstMCMean = ccdb->getForTimeStamp<TList>(pathMCMean, now);

      if (lstMCMean) {
        loadTProfile3DFromList(lstMCMean, "pmeanTruNchEtabinPtbin", pmeanTruNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanRecoNchEtabinPtbin", pmeanRecoNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanRecoMatchedNchEtabinPtbin", pmeanRecoMatchedNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanRecoEffcorrNchEtabinPtbin", pmeanRecoEffcorrNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanRecoMatchedEffcorrNchEtabinPtbin", pmeanRecoMatchedEffcorrNchEtabinPtbinStep2);

        loadTProfile3DFromList(lstMCMean, "pmeanEtTruNchEtabinPtbin", pmeanEtTruNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanEtRecoNchEtabinPtbin", pmeanEtRecoNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanEtRecoMatchedNchEtabinPtbin", pmeanEtRecoMatchedNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanEtRecoEffcorrNchEtabinPtbin", pmeanEtRecoEffcorrNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstMCMean, "pmeanEtRecoMatchedEffcorrNchEtabinPtbin", pmeanEtRecoMatchedEffcorrNchEtabinPtbinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for MC Mean from: %s", pathMCMean.c_str());
      }
    }

    if (cfgRunDataFluc) {
      LOGF(info, "Loading Data Mean profiles from CCDB path: %s", pathDataMean.c_str());
      TList* lstDataMean = ccdb->getForTimeStamp<TList>(pathDataMean, now);

      if (lstDataMean) {
        loadTProfile3DFromList(lstDataMean, "pmean_nch_etabin_ptbin", pmeanNchEtabinPtbinStep2);
        loadTProfile3DFromList(lstDataMean, "pmeanEt_nch_etabin_ptbin", pmeanEtNchEtabinPtbinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for Data Mean from: %s", pathDataMean.c_str());
      }
    }
    LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
  }

  void processGetEffHists(aod::McCollisions const& mcColl, soa::SmallGroups<MyRun3MCCollisions> const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
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
        if (trackSlice.size() < 1)
          continue;

        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (partSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;

        histos.fill(HIST("hZvtx_after_sel"), col.posZ());

        histos.fill(HIST("hCentrality"), cent);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), col.multNTracksPV(), tracks.size());
        histos.fill(HIST("Hist2D_cent_nch"), col.multNTracksPV(), cent);

        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle))
            continue;
          if (!particle.isPhysicalPrimary())
            continue;

          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == KPiPlus);
          const bool isKaon = (absPdgId == KKPlus);
          const bool isProton = (absPdgId == KProton);
          const bool isPid = (isPion || isKaon || isProton);

          histos.fill(HIST("hTruth_ParticleWeight"), cent, particle.pt(), particle.eta(), particle.weight());
          histos.fill(HIST("h3_AllPrimary"), col.multNTracksPV(), particle.pt(), particle.eta());
          histos.fill(HIST("h_AllPrimary"), particle.pt());

          if (isPid) {
            histos.fill(HIST("h3_AllPrimary_PID"), col.multNTracksPV(), particle.pt(), particle.eta());
          }
        }
        histos.fill(HIST("TruthTracKVz"), mcCollision.posZ(), col.posZ());
        histos.fill(HIST("vzResolution"), mcCollision.posZ(), mcCollision.posZ() - col.posZ());

        // Reconstructed
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);
          const bool isPid = (isPion || isKaon || isProton);
          histos.fill(HIST("hP"), track.p());
          histos.fill(HIST("hPt"), track.pt());
          histos.fill(HIST("hEta"), track.eta());
          histos.fill(HIST("hPhi"), track.phi());
          histos.fill(HIST("h_AllReco"), track.pt());

          histos.fill(HIST("h3_AllReco"), col.multNTracksPV(), track.pt(), track.eta());
          histos.fill(HIST("hEtaPhiReco"), col.posZ(), track.sign(), track.pt(), track.eta(), track.phi());

          histos.fill(HIST("hDCAxy_Reco"), track.dcaXY());
          histos.fill(HIST("hDCAz_Reco"), track.dcaZ());
          if (isPid) {
            histos.fill(HIST("h3_AllReco_PID"), col.multNTracksPV(), track.pt(), track.eta());
            histos.fill(HIST("hEtaPhiReco_PID"), col.posZ(), track.sign(), track.pt(), track.eta(), track.phi());
          }

          if (track.has_mcParticle()) {
            auto mcPart2 = track.mcParticle();
            if (mcPart2.isPhysicalPrimary()) {
              const int absPdgId = std::abs(mcPart2.pdgCode());
              const bool isPionTrue = (absPdgId == kPiPlus);
              const bool isKaonTrue = (absPdgId == kKPlus);
              const bool isProtonTrue = (absPdgId == kProton);
              const bool isPidTrue = (isPionTrue || isKaonTrue || isProtonTrue);

              histos.fill(HIST("hReco_ParticleWeight"), cent, mcPart2.pt(), mcPart2.eta(), mcPart2.weight());
              histos.fill(HIST("ptResolution"), mcPart2.pt(), mcPart2.pt() - track.pt());
              histos.fill(HIST("ptTruthReco"), mcPart2.pt(), track.pt());
              histos.fill(HIST("etaResolution"), mcPart2.eta(), mcPart2.eta() - track.eta());
              histos.fill(HIST("etaTruthReco"), mcPart2.eta(), track.eta());
              histos.fill(HIST("h3_RecoMatchedToPrimary"), col.multNTracksPV(), mcPart2.pt(), mcPart2.eta());
              histos.fill(HIST("h_RecoMatchedToPrimary"), mcPart2.pt());

              histos.fill(HIST("hDCAxy_RecoMatched"), track.dcaXY());
              histos.fill(HIST("hDCAz_RecoMatched"), track.dcaZ());

              if (isPid && isPidTrue) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PID"), col.multNTracksPV(), mcPart2.pt(), mcPart2.eta());
              }

            } else {
              // Matched to secondary
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), col.multNTracksPV(), track.pt(), track.eta());
              histos.fill(HIST("h_RecoUnMatchedToPrimary"), track.pt());
              histos.fill(HIST("hDCAxy_Unmatched"), track.dcaXY());
              histos.fill(HIST("hDCAz_Unmatched"), track.dcaZ());
              if (isPid) {
                histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PID"), col.multNTracksPV(), track.pt(), track.eta());
              }
            }
          } else {
            // Fake track
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), col.multNTracksPV(), track.pt(), track.eta());
            histos.fill(HIST("h_RecoUnMatchedToPrimary"), track.pt());
            histos.fill(HIST("hDCAxy_NotPrimary"), track.dcaXY());
            histos.fill(HIST("hDCAz_NotPrimary"), track.dcaZ());
            if (isPid) {
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PID"), col.multNTracksPV(), track.pt(), track.eta());
            }
          }
        } // tracks
      } // cols
    } // mcColl
    LOGF(info, "FINISHED RUNNING processGetEffHists");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate Eff and Fakes", cfgRunGetEff);

  void processMCFlat(aod::McCollisions const& mcColl, soa::SmallGroups<MyRun3MCCollisions> const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision())
          continue;
        if (!isEventSelected(col))
          continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1)
          continue;

        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (partSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;

        // Reconstructed
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          float effIncl = getEfficiency(col.multNTracksPV(), track.pt(), track.eta(), kInclusive, 0);
          float fakeIncl = getEfficiency(col.multNTracksPV(), track.pt(), track.eta(), kInclusive, 1);
          float wIncl = (1.0 - fakeIncl) / effIncl;
          if (!std::isfinite(wIncl) || wIncl <= 0.f)
            continue;
          if (effIncl <= 0 || !std::isfinite(effIncl) || !std::isfinite(fakeIncl))
            continue;
          histos.fill(HIST("hEtaPhiRecoEffWtd"), col.posZ(), track.sign(), track.pt(), track.eta(), track.phi(), wIncl);
          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);
          const bool isPid = (isPion || isKaon || isProton);
          float effPid = getEfficiency(col.multNTracksPV(), track.pt(), track.eta(), kCombinedPID, 0);
          float fakePid = getEfficiency(col.multNTracksPV(), track.pt(), track.eta(), kCombinedPID, 1);
          float wPid = (1.0 - fakePid) / effPid;
          if (effPid >= 1 || fakePid >= 1 || !std::isfinite(effPid) || effPid <= KFloatEpsilon || !std::isfinite(fakePid))
            continue;

          if (isPid) {
            histos.fill(HIST("hEtaPhiRecoEffWtd_PID"), col.posZ(), track.sign(), track.pt(), track.eta(), track.phi(), wPid);
          }
        } // tracks
      } // cols
    } // mcColl
    LOGF(info, "FINISHED RUNNING processMCFlat");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

  void processMCMean(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    float sumWiTruth[KNEta][KNpT], sumWiptiTruth[KNEta][KNpT];
    float sumWiReco[KNEta][KNpT], sumWiptiReco[KNEta][KNpT];
    float sumWiRecoEffCorr[KNEta][KNpT], sumWiptiRecoEffCorr[KNEta][KNpT];
    float sumWiTruthEt[KNEta][KNpT], sumWiptiTruthEt[KNEta][KNpT];
    float sumWiRecoEt[KNEta][KNpT], sumWiptiRecoEt[KNEta][KNpT];
    float sumWiRecoEffCorrEt[KNEta][KNpT], sumWiptiRecoEffCorrEt[KNEta][KNpT];

    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());

      if (colSlice.size() != 1)
        continue;
      for (const auto& col : colSlice) {
        if (!col.has_mcCollision())
          continue;
        if (!isEventSelected(col))
          continue;
        histos.fill(HIST("hVtxZ"), col.posZ());
        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1)
          continue;

        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (partSlice.size() < 1)
          continue;

        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;

        histos.fill(HIST("hZvtx_after_sel"), col.posZ());
        histos.fill(HIST("hCentrality"), cent);
        histos.fill(HIST("Hist2D_globalTracks_PVTracks"), col.multNTracksPV(), tracks.size());
        histos.fill(HIST("Hist2D_cent_nch"), col.multNTracksPV(), cent);
        histos.fill(HIST("TruthTracKVz"), mcCollision.posZ(), col.posZ());
        histos.fill(HIST("vzResolution"), mcCollision.posZ(), mcCollision.posZ() - col.posZ());

        // LOGF(info, "Event Check: cent = %.1f, nTracks = %d", cent, (int)trackSlice.size());
        memset(sumWiTruth, 0, sizeof(sumWiTruth));
        memset(sumWiptiTruth, 0, sizeof(sumWiptiTruth));
        memset(sumWiReco, 0, sizeof(sumWiReco));
        memset(sumWiptiReco, 0, sizeof(sumWiptiReco));
        memset(sumWiRecoEffCorr, 0, sizeof(sumWiRecoEffCorr));
        memset(sumWiptiRecoEffCorr, 0, sizeof(sumWiptiRecoEffCorr));
        memset(sumWiTruthEt, 0, sizeof(sumWiTruthEt));
        memset(sumWiptiTruthEt, 0, sizeof(sumWiptiTruthEt));
        memset(sumWiRecoEt, 0, sizeof(sumWiRecoEt));
        memset(sumWiptiRecoEt, 0, sizeof(sumWiptiRecoEt));
        memset(sumWiRecoEffCorrEt, 0, sizeof(sumWiRecoEffCorrEt));
        memset(sumWiptiRecoEffCorrEt, 0, sizeof(sumWiptiRecoEffCorrEt));

        // Truth
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle))
            continue;
          if (!particle.isPhysicalPrimary())
            continue;

          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == kPiPlus);
          const bool isKaon = (absPdgId == kKPlus);
          const bool isProton = (absPdgId == kProton);

          float pt = particle.pt();
          float eta = particle.eta();
          float p = particle.p();

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              sumWiTruth[ieta][ipt]++;
              sumWiptiTruth[ieta][ipt] += pt;
              if (isPion || isKaon || isProton) {
                float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                               : o2::constants::physics::MassProton;
                float energy = std::sqrt(p * p + m * m);
                float et = energy * (pt / p); // E_T = E * sin(theta) = E * (pT / p)
                sumWiTruthEt[ieta][ipt]++;
                sumWiptiTruthEt[ieta][ipt] += et;
              }
            }
          }
        }

        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;

          float pt = track.pt();
          float eta = track.eta();
          float p = track.p();
          float phi = track.phi();

          histos.fill(HIST("hEtaPhiReco"), col.posZ(), track.sign(), track.pt(), eta, phi);

          float effIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 0);
          float fakeIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 1);
          float flatWeightIncl = getFlatteningWeight(col.posZ(), track.sign(), pt, eta, phi, kInclusive, cfgFlat);
          float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
          if (!std::isfinite(wIncl) || wIncl <= 0.f)
            continue;
          if (effIncl <= 0 || !std::isfinite(effIncl) || !std::isfinite(fakeIncl) || !std::isfinite(flatWeightIncl))
            continue;
          histos.fill(HIST("hEtaPhiRecoWtd"), col.posZ(), track.sign(), pt, eta, phi, flatWeightIncl);
          histos.fill(HIST("hEtaPhiRecoEffWtd"), col.posZ(), track.sign(), pt, eta, phi, wIncl);

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              sumWiReco[ieta][ipt] += 1.0;
              sumWiptiReco[ieta][ipt] += pt;
            }
          }

          if (effIncl <= 0 || !std::isfinite(wIncl) || !std::isfinite(fakeIncl) || !std::isfinite(flatWeightIncl))
            continue;

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              sumWiRecoEffCorr[ieta][ipt] += wIncl;
              sumWiptiRecoEffCorr[ieta][ipt] += wIncl * pt;
            }
          }

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);
          if (isPion || isKaon || isProton) {
            histos.fill(HIST("hEtaPhiReco_PID"), col.posZ(), track.sign(), track.pt(), eta, phi);
            float effPid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 0);
            float fakePid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 1);
            float flatWeightPid = getFlatteningWeight(col.posZ(), track.sign(), pt, eta, phi, kCombinedPID, cfgFlat);
            float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
            if (effPid >= 1 || fakePid >= 1 || !std::isfinite(effPid) || effPid <= KFloatEpsilon || !std::isfinite(fakePid) || !std::isfinite(flatWeightPid))
              continue;
            histos.fill(HIST("hEtaPhiRecoWtd_PID"), col.posZ(), track.sign(), track.pt(), eta, track.phi(), flatWeightPid);
            histos.fill(HIST("hEtaPhiRecoEffWtd_PID"), col.posZ(), track.sign(), track.pt(), eta, track.phi(), wPid);

            float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                           : o2::constants::physics::MassProton;
            float energy = std::sqrt(p * p + m * m);
            float et = energy * (pt / p); // E_T = E * sin(theta)
            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                sumWiRecoEt[ieta][ipt] += 1.0;
                sumWiptiRecoEt[ieta][ipt] += et;
              }
            }

            if (effPid <= KFloatEpsilon || !std::isfinite(wPid) || !std::isfinite(fakePid) || !std::isfinite(flatWeightPid))
              continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                sumWiRecoEffCorrEt[ieta][ipt] += wPid;
                sumWiptiRecoEffCorrEt[ieta][ipt] += wPid * et;
              }
            }
          }

          if (std::isfinite(wIncl)) {
            if (cent < KCentTestMin) {
              histos.fill(HIST("wgt_pT"), pt, wIncl);
              histos.fill(HIST("Eff_pT"), pt, effIncl);
              histos.fill(HIST("Fake_pT"), pt, fakeIncl);
              histos.fill(HIST("Eff_eta"), eta, effIncl);
              histos.fill(HIST("Fake_eta"), eta, fakeIncl);
              histos.fill(HIST("wgt_eta"), eta, wIncl);
            }
            histos.fill(HIST("Eff_cent"), cent, effIncl);
            histos.fill(HIST("Eff_Ntrk"), col.multNTracksPV(), effIncl);
            histos.fill(HIST("Fake_cent"), cent, fakeIncl);
            histos.fill(HIST("Fake_Ntrk"), col.multNTracksPV(), fakeIncl);
            histos.fill(HIST("wgt_cent"), cent, wIncl);
            histos.fill(HIST("wgt_Ntrk"), col.multNTracksPV(), wIncl);
          }

        } // end track loop

        if (std::isfinite(sumWiTruth[0][0])) {
          float meanPtTruth = sumWiptiTruth[0][0] / sumWiTruth[0][0];
          if (!std::isfinite(meanPtTruth))
            LOGF(info, "meanPtTruth = %.3f, num = %.3f, den =%.3f", meanPtTruth, sumWiptiTruth[0][0], sumWiTruth[0][0]);
          if (!std::isfinite(meanPtTruth))
            continue;
          histos.fill(HIST("MCGen/Prof_cent_Nchrec"), cent, sumWiTruth[0][0]);
          histos.fill(HIST("MCGen/Prof_MeanpT_Cent"), cent, meanPtTruth);
          histos.fill(HIST("MCGen/Prof_MeanpT_Mult"), col.multNTracksPV(), meanPtTruth);
        }
        if (std::isfinite(sumWiReco[0][0])) {
          float meanPtReco = sumWiptiReco[0][0] / sumWiReco[0][0];
          if (!std::isfinite(meanPtReco))
            LOGF(info, "meanPtReco = %.3f, num = %.3f, den =%.3f", meanPtReco, sumWiptiReco[0][0], sumWiReco[0][0]);
          if (!std::isfinite(meanPtReco))
            continue;
          histos.fill(HIST("MCReco/Prof_cent_Nchrec"), cent, sumWiReco[0][0]);
          histos.fill(HIST("MCReco/Prof_MeanpT_Cent"), cent, meanPtReco);
          histos.fill(HIST("MCReco/Prof_MeanpT_Mult"), col.multNTracksPV(), meanPtReco);
        }
        if (std::isfinite(sumWiRecoEffCorr[0][0])) {
          float meanpTeffcorr = sumWiptiRecoEffCorr[0][0] / sumWiRecoEffCorr[0][0];
          if (!std::isfinite(meanpTeffcorr))
            LOGF(info, "meanPtRecoEffcorr = %.3f, num = %.3f, den =%.3f", meanpTeffcorr, sumWiptiRecoEffCorr[0][0], sumWiRecoEffCorr[0][0]);
          if (!std::isfinite(meanpTeffcorr))
            continue;
          histos.fill(HIST("MCRecoEffCorr/Prof_cent_Nchrec"), cent, sumWiRecoEffCorr[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent"), cent, meanpTeffcorr);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult"), col.multNTracksPV(), meanpTeffcorr);
        }

        if (std::isfinite(sumWiTruthEt[0][0])) {
          float meanEt = sumWiptiTruthEt[0][0] / sumWiTruthEt[0][0];
          if (!std::isfinite(meanEt))
            LOGF(info, "meanEtTruthEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiTruthEt[0][0], sumWiTruthEt[0][0]);
          if (!std::isfinite(meanEt))
            continue;
          histos.fill(HIST("MCGen/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCGen/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }
        // "MCReco"
        if (std::isfinite(sumWiRecoEt[0][0])) {
          float meanEt = sumWiptiRecoEt[0][0] / sumWiRecoEt[0][0];
          if (!std::isfinite(meanEt))
            LOGF(info, "meanEtRecoEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiRecoEt[0][0], sumWiRecoEt[0][0]);
          if (!std::isfinite(meanEt))
            continue;
          histos.fill(HIST("MCReco/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCReco/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }
        // "MCRecoEffCorr"
        if (std::isfinite(sumWiRecoEffCorrEt[0][0])) {
          float meanEt = sumWiptiRecoEffCorrEt[0][0] / sumWiRecoEffCorrEt[0][0];
          if (!std::isfinite(meanEt))
            LOGF(info, "meanEtRecoEffcorrEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiRecoEffCorrEt[0][0], sumWiRecoEffCorrEt[0][0]);
          if (!std::isfinite(meanEt))
            continue;
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (std::isfinite(sumWiTruth[ieta][ipt]))
              histos.fill(HIST("pmeanTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiTruth[ieta][ipt] / sumWiTruth[ieta][ipt]);
            if (std::isfinite(sumWiReco[ieta][ipt]))
              histos.fill(HIST("pmeanRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiReco[ieta][ipt] / sumWiReco[ieta][ipt]);
            if (std::isfinite(sumWiRecoEffCorr[ieta][ipt]))
              histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiRecoEffCorr[ieta][ipt] / sumWiRecoEffCorr[ieta][ipt]);

            if (std::isfinite(sumWiTruthEt[ieta][ipt]))
              histos.fill(HIST("pmeanEtTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiTruthEt[ieta][ipt] / sumWiTruthEt[ieta][ipt]);
            if (std::isfinite(sumWiRecoEt[ieta][ipt]))
              histos.fill(HIST("pmeanEtRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiRecoEt[ieta][ipt] / sumWiRecoEt[ieta][ipt]);
            if (std::isfinite(sumWiRecoEffCorrEt[ieta][ipt]))
              histos.fill(HIST("pmeanEtRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiRecoEffCorrEt[ieta][ipt] / sumWiRecoEffCorrEt[ieta][ipt]);
          }
        }
      } // end col loop
    }
    LOGF(info, "FINISHED RUNNING processMCMean (pT + Et)");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt/Et and Eff Hists", cfgRunMCMean);

  void processMCFluc(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    double sumPmwkTru[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkTru[KNEta][KNpT][KIntK]{};
    double sumPmwkReco[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkReco[KNEta][KNpT][KIntK]{};
    double sumPmwkRecoEffCor[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkRecoEffCor[KNEta][KNpT][KIntK]{};
    double sumPmwkTruEt[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkTruEt[KNEta][KNpT][KIntK]{};
    double sumPmwkRecoEt[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkRecoEt[KNEta][KNpT][KIntK]{};
    double sumPmwkRecoEffCorEt[KNEta][KNpT][KIntM][KIntK]{};
    double sumWkRecoEffCorEt[KNEta][KNpT][KIntK]{};
    double meanTru[KNEta][KNpT]{}, c2Tru[KNEta][KNpT]{};
    double meanReco[KNEta][KNpT]{}, c2Reco[KNEta][KNpT]{};
    double meanRecoEffCor[KNEta][KNpT]{}, c2RecoEffCor[KNEta][KNpT]{};
    double meanTruEt[KNEta][KNpT]{}, c2TruEt[KNEta][KNpT]{};
    double meanRecoEt[KNEta][KNpT]{}, c2RecoEt[KNEta][KNpT]{};
    double meanRecoEffCorEt[KNEta][KNpT]{}, c2RecoEffCorEt[KNEta][KNpT]{};

    for (const auto& mcCollision : mcColl) {
      auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1)
        continue;
      // histos.fill(HIST("MCGen/hVtxZ"), mcCollision.posZ());
      for (const auto& col : colSlice) {

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1)
          continue;

        memset(sumPmwkTru, 0, sizeof(sumPmwkTru));
        memset(sumWkTru, 0, sizeof(sumWkTru));
        memset(sumPmwkReco, 0, sizeof(sumPmwkReco));
        memset(sumWkReco, 0, sizeof(sumWkReco));
        memset(sumPmwkRecoEffCor, 0, sizeof(sumPmwkRecoEffCor));
        memset(sumWkRecoEffCor, 0, sizeof(sumWkRecoEffCor));

        memset(sumPmwkTruEt, 0, sizeof(sumPmwkTruEt));
        memset(sumWkTruEt, 0, sizeof(sumWkTruEt));
        memset(sumPmwkRecoEt, 0, sizeof(sumPmwkRecoEt));
        memset(sumWkRecoEt, 0, sizeof(sumWkRecoEt));
        memset(sumPmwkRecoEffCorEt, 0, sizeof(sumPmwkRecoEffCorEt));
        memset(sumWkRecoEffCorEt, 0, sizeof(sumWkRecoEffCorEt));

        memset(meanTru, 0, sizeof(meanTru));
        memset(c2Tru, 0, sizeof(c2Tru));
        memset(meanReco, 0, sizeof(meanReco));
        memset(c2Reco, 0, sizeof(c2Reco));
        memset(meanRecoEffCor, 0, sizeof(meanRecoEffCor));
        memset(c2RecoEffCor, 0, sizeof(c2RecoEffCor));

        memset(meanTruEt, 0, sizeof(meanTruEt));
        memset(c2TruEt, 0, sizeof(c2TruEt));
        memset(meanRecoEt, 0, sizeof(meanRecoEt));
        memset(c2RecoEt, 0, sizeof(c2RecoEt));
        memset(meanRecoEffCorEt, 0, sizeof(meanRecoEffCorEt));
        memset(c2RecoEffCorEt, 0, sizeof(c2RecoEffCorEt));

        if (!col.has_mcCollision() || !isEventSelected(col))
          continue;
        float cent = getCentrality(col);
        if (cent > KCentMax)
          continue;

        // truth
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle))
            continue;
          if (!particle.isPhysicalPrimary())
            continue;
          float pt = particle.pt();
          float eta = particle.eta();
          float p = particle.p();

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              for (int k = 0; k < KIntK; ++k) {
                for (int m = 0; m < KIntM; ++m) {
                  sumPmwkTru[ieta][ipt][m][k] += std::pow(pt, m);
                }
                sumWkTru[ieta][ipt][k]++;
              }
            }
          }
          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == kPiPlus);
          const bool isKaon = (absPdgId == kKPlus);
          const bool isProton = (absPdgId == kProton);
          if (isPion || isKaon || isProton) {

            float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                           : o2::constants::physics::MassProton;
            float energy = std::sqrt(p * p + m * m);
            float et = energy * (pt / p);
            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                for (int k = 0; k < KIntK; ++k) {
                  for (int m = 0; m < KIntM; ++m) {
                    sumPmwkTruEt[ieta][ipt][m][k] += std::pow(et, m);
                  }
                  sumWkTruEt[ieta][ipt][k]++;
                }
              }
            }
          }
        } // end truth loop

        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track))
            continue;
          float pt = track.pt();
          float eta = track.eta();
          float p = track.p();
          float phi = track.phi();

          float effIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 0);
          float fakeIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 1);
          float flatWeightIncl = getFlatteningWeight(col.posZ(), track.sign(), pt, eta, phi, kInclusive, cfgFlat);
          float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
          if (!std::isfinite(wIncl) || wIncl <= 0.f)
            continue;
          if (effIncl <= 0 || !std::isfinite(effIncl) || !std::isfinite(fakeIncl) || !std::isfinite(flatWeightIncl))
            continue;

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta])
              continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                continue;
              for (int k = 0; k < KIntK; ++k) {
                for (int m = 0; m < KIntM; ++m) {
                  sumPmwkReco[ieta][ipt][m][k] += std::pow(1.0, k) * std::pow(pt, m);
                  sumPmwkRecoEffCor[ieta][ipt][m][k] += std::pow(wIncl, k) * std::pow(pt, m);
                }
                sumWkReco[ieta][ipt][k] += std::pow(1.0, k);
                sumWkRecoEffCor[ieta][ipt][k] += std::pow(wIncl, k);
              }
            }
          }

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);

          if (isPion || isKaon || isProton) {
            float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                           : o2::constants::physics::MassProton;
            float energy = std::sqrt(p * p + m * m);
            float et = energy * (pt / p); // E_T = E * sin(theta)
            float effPid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 0);
            float fakePid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 1);
            float flatWeightPid = getFlatteningWeight(col.posZ(), track.sign(), pt, eta, phi, kCombinedPID, cfgFlat);
            float wPid = flatWeightPid * (1 - fakePid) / effPid;
            if (effPid >= 1 || fakePid >= 1 || !std::isfinite(effPid) || effPid <= KFloatEpsilon || !std::isfinite(fakePid) || !std::isfinite(flatWeightPid))
              continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt])
                  continue;
                for (int k = 0; k < KIntK; ++k) {
                  for (int m = 0; m < KIntM; ++m) {
                    sumPmwkRecoEt[ieta][ipt][m][k] += std::pow(1.0, k) * std::pow(et, m);
                    sumPmwkRecoEffCorEt[ieta][ipt][m][k] += std::pow(wPid, k) * std::pow(et, m);
                  }
                  sumWkRecoEt[ieta][ipt][k] += std::pow(1.0, k);
                  sumWkRecoEffCorEt[ieta][ipt][k] += std::pow(wPid, k);
                }
              }
            }
          }

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              const int ibx = pmeanTruNchEtabinPtbinStep2->GetXaxis()->FindBin(col.multNTracksPV());
              const int iby = ieta + 1;
              const int ibz = ipt + 1;

              float mmptTru = pmeanTruNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
              float mmptReco = pmeanRecoNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
              float mmptRecoEffCor = pmeanRecoEffcorrNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
              float mmetTru = pmeanEtTruNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
              float mmetReco = pmeanEtRecoNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
              float mmetRecoEffCor = pmeanEtRecoEffcorrNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);

              if (std::isfinite(mmptTru))
                std::tie(meanTru[ieta][ipt], c2Tru[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTru[ieta][ipt], sumWkTru[ieta][ipt], mmptTru);
              if (std::isfinite(mmptReco))
                std::tie(meanReco[ieta][ipt], c2Reco[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkReco[ieta][ipt], sumWkReco[ieta][ipt], mmptReco);
              if (std::isfinite(mmptRecoEffCor))
                std::tie(meanRecoEffCor[ieta][ipt], c2RecoEffCor[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCor[ieta][ipt], sumWkRecoEffCor[ieta][ipt], mmptRecoEffCor);

              if (std::isfinite(mmetTru))
                std::tie(meanTruEt[ieta][ipt], c2TruEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTruEt[ieta][ipt], sumWkTruEt[ieta][ipt], mmetTru);
              if (std::isfinite(mmetReco))
                std::tie(meanRecoEt[ieta][ipt], c2RecoEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEt[ieta][ipt], sumWkRecoEt[ieta][ipt], mmetReco);
              if (std::isfinite(mmetRecoEffCor))
                std::tie(meanRecoEffCorEt[ieta][ipt], c2RecoEffCorEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCorEt[ieta][ipt], sumWkRecoEffCorEt[ieta][ipt], mmetRecoEffCor);
            }
          }
        }
        if (std::isfinite(c2Tru[0][0])) {
          histos.fill(HIST("MCGen/Prof_C2_Cent"), cent, c2Tru[0][0]);
          histos.fill(HIST("MCGen/Prof_C2_Mult"), col.multNTracksPV(), c2Tru[0][0]);
        }
        if (std::isfinite(c2TruEt[0][0])) {
          histos.fill(HIST("MCGen/Prof_C2Et_Cent"), cent, c2TruEt[0][0]);
          histos.fill(HIST("MCGen/Prof_C2Et_Mult"), col.multNTracksPV(), c2TruEt[0][0]);
        }
        // "MCReco"
        if (std::isfinite(c2Reco[0][0])) {
          histos.fill(HIST("MCReco/Prof_C2_Cent"), cent, c2Reco[0][0]);
          histos.fill(HIST("MCReco/Prof_C2_Mult"), col.multNTracksPV(), c2Reco[0][0]);
        }
        if (std::isfinite(c2RecoEt[0][0])) {
          histos.fill(HIST("MCReco/Prof_C2Et_Cent"), cent, c2RecoEt[0][0]);
          histos.fill(HIST("MCReco/Prof_C2Et_Mult"), col.multNTracksPV(), c2RecoEt[0][0]);
        }

        if (std::isfinite(c2RecoEffCor[0][0])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent"), cent, c2RecoEffCor[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult"), col.multNTracksPV(), c2RecoEffCor[0][0]);
        }
        if (std::isfinite(c2RecoEffCorEt[0][0])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_C2Et_Cent"), cent, c2RecoEffCorEt[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_C2Et_Mult"), col.multNTracksPV(), c2RecoEffCorEt[0][0]);
        }

        if (std::isfinite(sumWkTru[0][0][1])) {
          histos.fill(HIST("MCGen/Prof_cent_Nchrec"), cent, sumWkTru[0][0][1]);
          histos.fill(HIST("MCGen/Prof_MeanpT_Cent"), cent, meanTru[0][0]);
          histos.fill(HIST("MCGen/Prof_MeanpT_Mult"), col.multNTracksPV(), meanTru[0][0]);
        }
        if (std::isfinite(sumWkTruEt[0][0][1])) {
          histos.fill(HIST("MCGen/Prof_MeanEt_Cent"), cent, meanTruEt[0][0]);
          histos.fill(HIST("MCGen/Prof_MeanEt_Mult"), col.multNTracksPV(), meanTruEt[0][0]);
        }
        // "MCReco"
        if (std::isfinite(sumWkReco[0][0][1])) {
          histos.fill(HIST("MCReco/Prof_cent_Nchrec"), cent, sumWkReco[0][0][1]);
          histos.fill(HIST("MCReco/Prof_MeanpT_Cent"), cent, meanReco[0][0]);
          histos.fill(HIST("MCReco/Prof_MeanpT_Mult"), col.multNTracksPV(), meanReco[0][0]);
        }
        if (std::isfinite(sumWkRecoEt[0][0][1])) {
          histos.fill(HIST("MCReco/Prof_MeanEt_Cent"), cent, meanRecoEt[0][0]);
          histos.fill(HIST("MCReco/Prof_MeanEt_Mult"), col.multNTracksPV(), meanRecoEt[0][0]);
        }
        // "MCRecoEffCorr"
        if (std::isfinite(sumWkRecoEffCor[0][0][1])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_cent_Nchrec"), cent, sumWkRecoEffCor[0][0][1]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent"), cent, meanRecoEffCor[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult"), col.multNTracksPV(), meanRecoEffCor[0][0]);
        }
        if (std::isfinite(sumWkRecoEffCorEt[0][0][1])) {
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Cent"), cent, meanRecoEffCorEt[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Mult"), col.multNTracksPV(), meanRecoEffCorEt[0][0]);
        }

        for (int ieta = 0; ieta < KNEta; ++ieta) {
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (std::isfinite(sumWkTru[ieta][ipt][1]))
              histos.fill(HIST("pmeanTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanTru[ieta][ipt]);
            if (std::isfinite(sumWkReco[ieta][ipt][1]))
              histos.fill(HIST("pmeanRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanReco[ieta][ipt]);
            if (std::isfinite(sumWkRecoEffCor[ieta][ipt][1]))
              histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanRecoEffCor[ieta][ipt]);
            if (std::isfinite(sumWkTruEt[ieta][ipt][1]))
              histos.fill(HIST("pmeanEtTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanTruEt[ieta][ipt]);
            if (std::isfinite(sumWkRecoEt[ieta][ipt][1]))
              histos.fill(HIST("pmeanEtRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanRecoEt[ieta][ipt]);
            if (std::isfinite(sumWkRecoEffCorEt[ieta][ipt][1]))
              histos.fill(HIST("pmeanEtRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, meanRecoEffCorEt[ieta][ipt]);
          }
        }

        float p1kBarTru[KNEta][KNpT]{}, p1kBarReco[KNEta][KNpT]{}, p1kBarRecoEffCor[KNEta][KNpT]{};
        float p1kBarTruEt[KNEta][KNpT]{}, p1kBarRecoEt[KNEta][KNpT]{}, p1kBarRecoEffCorEt[KNEta][KNpT]{};
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            const int ibx = pmeanTruNchEtabinPtbinStep2->GetXaxis()->FindBin(col.multNTracksPV());
            const int iby = ieta + 1;
            const int ibz = ipt + 1;

            float mmptTru = pmeanTruNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
            float mmptReco = pmeanRecoNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
            float mmptRecoEffCor = pmeanRecoEffcorrNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
            float mmetTru = pmeanEtTruNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
            float mmetReco = pmeanEtRecoNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
            float mmetRecoEffCor = pmeanEtRecoEffcorrNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);

            if (mmptTru != 0.0f)
              p1kBarTru[ieta][ipt] = meanTru[ieta][ipt] - mmptTru;
            if (mmptReco != 0.0f)
              p1kBarReco[ieta][ipt] = meanReco[ieta][ipt] - mmptReco;
            if (mmptRecoEffCor != 0.0f)
              p1kBarRecoEffCor[ieta][ipt] = meanRecoEffCor[ieta][ipt] - mmptRecoEffCor;

            if (mmetTru != 0.0f)
              p1kBarTruEt[ieta][ipt] = meanTruEt[ieta][ipt] - mmetTru;
            if (mmetReco != 0.0f)
              p1kBarRecoEt[ieta][ipt] = meanRecoEt[ieta][ipt] - mmetReco;
            if (mmetRecoEffCor != 0.0f)
              p1kBarRecoEffCorEt[ieta][ipt] = meanRecoEffCorEt[ieta][ipt] - mmetRecoEffCor;
          }
        }

        // --- 1D Covariance ---
        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

          for (int ipt = 0; ipt < KNpT; ++ipt) {
            float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCGen/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 1:
                  histos.fill(HIST("MCGen/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 2:
                  histos.fill(HIST("MCGen/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCGen/Prof_ipt0_Cov_Eta"), valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCGen/Prof_ipt1_Cov_Eta"), valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCGen/Prof_ipt2_Cov_Eta"), valy, c2Sub);
                    break;
                }
              }
            }

            float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCGen/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCGen/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 1:
                  histos.fill(HIST("MCGen/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 2:
                  histos.fill(HIST("MCGen/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCGen/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCGen/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCGen/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }

        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

          for (int ipt = 0; ipt < KNpT; ++ipt) {
            float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCReco/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 1:
                  histos.fill(HIST("MCReco/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 2:
                  histos.fill(HIST("MCReco/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCReco/Prof_ipt0_Cov_Eta"), valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCReco/Prof_ipt1_Cov_Eta"), valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCReco/Prof_ipt2_Cov_Eta"), valy, c2Sub);
                    break;
                }
              }
            }

            float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCReco/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCReco/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 1:
                  histos.fill(HIST("MCReco/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 2:
                  histos.fill(HIST("MCReco/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCReco/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCReco/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCReco/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }

        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

          for (int ipt = 0; ipt < KNpT; ++ipt) {
            float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 1:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
                case 2:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov_Eta"), valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov_Eta"), valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov_Eta"), valy, c2Sub);
                    break;
                }
              }
            }

            float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              switch (ipt) {
                case 0:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 1:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
                case 2:
                  histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
                  break;
              }
              if (cent < KCentCovCut) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }

        // --- 2D Correlations ---
        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = KHalf * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

            for (int ipt = 0; ipt < KNpT; ++ipt) {
              float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
              if (std::isfinite(c2Sub)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                }
              }
              float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCGen/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCGen/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCGen/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }

        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = KHalf * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

            for (int ipt = 0; ipt < KNpT; ++ipt) {
              float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
              if (std::isfinite(c2Sub)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                }
              }
              float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCReco/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCReco/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCReco/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }

        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = KHalf * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

            for (int ipt = 0; ipt < KNpT; ++ipt) {
              float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
              if (std::isfinite(c2Sub)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 1:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                  case 2:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
                    break;
                }
              }
              float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) {
                switch (ipt) {
                  case 0:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 1:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                  case 2:
                    histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
                    break;
                }
              }
            }
          }
        }
      }
    }
    LOGF(info, "FINISHED RUNNING processMCFluc (pT + Et)");
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt/Et fluc", cfgRunMCFluc);

  void processGetDataFlat(AodCollisionsSel::iterator const& coll, aod::BCsWithTimestamps const&, AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;
      float p = track.p();
      float pt = track.pt();
      float eta = track.eta();
      float phi = track.phi();
      if (p < KFloatEpsilon)
        continue;
      float effIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 0);
      float fakeIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 1);
      float wIncl = (1.0 - fakeIncl) / effIncl;
      if (!std::isfinite(wIncl) || wIncl <= KFloatEpsilon || effIncl <= KFloatEpsilon)
        continue;
      histos.fill(HIST("hEtaPhiReco"), coll.posZ(), track.sign(), pt, eta, phi);
      histos.fill(HIST("hEtaPhiRecoEffWtd"), coll.posZ(), track.sign(), pt, eta, phi, wIncl);

      const bool isPion = selectionPion(track);
      const bool isKaon = selectionKaon(track);
      const bool isProton = selectionProton(track);
      if (isPion || isKaon || isProton) {
        float effPid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 0);
        float fakePid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 1);
        float wPid = (1.0 - fakePid) / effPid;
        if (!std::isfinite(wPid) || wPid <= KFloatEpsilon || effPid <= KFloatEpsilon)
          continue;
        histos.fill(HIST("hEtaPhiReco_PID"), coll.posZ(), track.sign(), pt, eta, phi);
        histos.fill(HIST("hEtaPhiRecoEffWtd_PID"), coll.posZ(), track.sign(), pt, eta, phi, wPid);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetDataFlat, "process data to calculate Flattening maps", cfgRunGetDataFlat);

  void processDataMean(AodCollisionsSel::iterator const& coll, aod::BCsWithTimestamps const&, AodTracksSel const& tracks)
  {
    float sumWi[KNEta][KNpT]{}, sumWipti[KNEta][KNpT]{};
    float sumWiEt[KNEta][KNpT]{}, sumWiEtVal[KNEta][KNpT]{};
    if (!isEventSelected(coll))
      return;

    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;

    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;
      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      if (p < KFloatEpsilon)
        continue;
      histos.fill(HIST("hP"), p);
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), track.phi());

      float effIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 0);
      float fakeIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 1);
      float flatWeightIncl = getFlatteningWeight(coll.posZ(), track.sign(), pt, eta, phi, kInclusive, cfgFlat);
      float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
      if (!std::isfinite(wIncl) || wIncl <= KFloatEpsilon || effIncl <= KFloatEpsilon)
        continue;

      histos.fill(HIST("hEtaPhiReco"), coll.posZ(), track.sign(), pt, eta, track.phi());
      histos.fill(HIST("hEtaPhiRecoEffWtd"), coll.posZ(), track.sign(), eta, pt, track.phi(), (1.0 - fakeIncl) / effIncl);
      histos.fill(HIST("hEtaPhiRecoWtd"), coll.posZ(), track.sign(), eta, pt, track.phi(), wIncl);

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta])
          continue;
        for (int ipt = 0; ipt < KNpT; ++ipt) {
          if (pt <= pTLw[ipt] || pt > pTUp[ipt])
            continue;
          sumWi[ieta][ipt] += wIncl;
          sumWipti[ieta][ipt] += wIncl * pt;
        }
      }

      const bool isPion = selectionPion(track);
      const bool isKaon = selectionKaon(track);
      const bool isProton = selectionProton(track);
      if (isPion || isKaon || isProton) {
        float effPid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 0);
        float fakePid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 1);
        float flatWeightPid = getFlatteningWeight(coll.posZ(), track.sign(), pt, eta, phi, kCombinedPID, cfgFlat);
        float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
        if (!std::isfinite(wPid) || wPid <= KFloatEpsilon || effPid <= KFloatEpsilon)
          continue;

        histos.fill(HIST("hEtaPhiReco_PID"), coll.posZ(), track.sign(), pt, eta, track.phi());
        histos.fill(HIST("hEtaPhiRecoEffWtd_PID"), coll.posZ(), track.sign(), eta, pt, track.phi(), (1.0 - fakePid) / effPid);
        histos.fill(HIST("hEtaPhiRecoWtd_PID"), coll.posZ(), track.sign(), eta, pt, track.phi(), wPid);

        float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                       : o2::constants::physics::MassProton;
        float energy = std::sqrt(p * p + m * m);
        float et = energy * (pt / p); // E_T = E * sin(theta)
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt])
              continue;
            sumWiEt[ieta][ipt] += wPid;
            sumWiEtVal[ieta][ipt] += wPid * et;
          }
        }
      }
    }
    histos.fill(HIST("Prof_cent_Nchrec"), cent, sumWi[0][0]);
    if (std::isfinite(sumWi[0][0]))
      histos.fill(HIST("Prof_MeanpT_Cent"), cent, sumWipti[0][0] / sumWi[0][0]);
    if (std::isfinite(sumWiEt[0][0]))
      histos.fill(HIST("Prof_MeanEt_Cent"), cent, sumWiEtVal[0][0] / sumWiEt[0][0]);

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        if (std::isfinite(sumWi[ieta][ipt]))
          histos.fill(HIST("pmean_nch_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, sumWipti[ieta][ipt] / sumWi[ieta][ipt]);
        if (std::isfinite(sumWiEt[ieta][ipt]))
          histos.fill(HIST("pmeanEt_nch_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, sumWiEtVal[ieta][ipt] / sumWiEt[ieta][ipt]);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataMean, "process data to calculate mean pT and Et", cfgRunDataMean);

  void processDataFluc(AodCollisionsSel::iterator const& coll, aod::BCsWithTimestamps const&, AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll))
      return;
    float cent = getCentrality(coll);
    if (cent > KCentMax)
      return;
    if (!pmeanNchEtabinPtbinStep2 || !pmeanEtNchEtabinPtbinStep2) {
      LOGF(warning, "Data fluc: Mean pT or Et map missing");
      return;
    }

    if (!hEff[kInclusive] || !hFake[kInclusive] || !hFlatWeight[kInclusive] || !hEff[kCombinedPID] || !hFake[kCombinedPID] || !hFlatWeight[kCombinedPID]) {
      LOGF(warning, "Data fluc: Inclusive or PID correction maps are null");
      return;
    }
    double sumpmwk[KNEta][KNpT][KIntM][KIntK]{};
    double sumwk[KNEta][KNpT][KIntK]{};
    double sumpmwkEt[KNEta][KNpT][KIntM][KIntK]{};
    double sumwkEt[KNEta][KNpT][KIntK]{};
    double mean[KNEta][KNpT]{}, c2[KNEta][KNpT]{};
    double p1kBar[KNEta][KNpT]{};
    double meanEt[KNEta][KNpT]{}, c2Et[KNEta][KNpT]{};
    double p1kBarEt[KNEta][KNpT]{};

    for (const auto& track : tracks) {
      if (!isTrackSelected(track))
        continue;
      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      if (p < KFloatEpsilon)
        continue;

      float effIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 0);
      float fakeIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 1);
      float flatWeightIncl = getFlatteningWeight(coll.posZ(), track.sign(), pt, eta, phi, kInclusive, cfgFlat);
      float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
      if (!std::isfinite(wIncl) || wIncl <= KFloatEpsilon || effIncl <= KFloatEpsilon)
        continue;
      histos.fill(HIST("hEtaPhiReco"), coll.posZ(), track.sign(), pt, eta, track.phi());
      histos.fill(HIST("hEtaPhiRecoEffWtd"), coll.posZ(), track.sign(), eta, pt, track.phi(), (1.0 - fakeIncl) / effIncl);
      histos.fill(HIST("hEtaPhiRecoWtd"), coll.posZ(), track.sign(), eta, pt, track.phi(), wIncl);

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta])
          continue;
        for (int ipt = 0; ipt < KNpT; ++ipt) {
          if (pt <= pTLw[ipt] || pt > pTUp[ipt])
            continue;
          for (int k = 0; k < KIntK; ++k) {
            for (int m = 0; m < KIntM; ++m)
              sumpmwk[ieta][ipt][m][k] += std::pow(wIncl, k) * std::pow(pt, m);
            sumwk[ieta][ipt][k] += std::pow(wIncl, k);
          }
        }
      }

      const bool isPion = selectionPion(track);
      const bool isKaon = selectionKaon(track);
      const bool isProton = selectionProton(track);
      if (isPion || isKaon || isProton) {
        float effPid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 0);
        float fakePid = getEfficiency(coll.multNTracksPV(), pt, eta, kCombinedPID, 1);
        float flatWeightPid = getFlatteningWeight(coll.posZ(), track.sign(), pt, eta, phi, kCombinedPID, cfgFlat);
        float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
        if (!std::isfinite(wPid) || wPid <= KFloatEpsilon || effPid <= KFloatEpsilon)
          continue;
        histos.fill(HIST("hEtaPhiReco_PID"), coll.posZ(), track.sign(), pt, eta, track.phi());
        histos.fill(HIST("hEtaPhiRecoEffWtd_PID"), coll.posZ(), track.sign(), eta, pt, track.phi(), (1.0 - fakePid) / effPid);
        histos.fill(HIST("hEtaPhiRecoWtd_PID"), coll.posZ(), track.sign(), eta, pt, track.phi(), wPid);

        float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus
                                                                       : o2::constants::physics::MassProton;

        float energy = std::sqrt(p * p + m * m);
        float et = energy * (pt / p); // E_T = E * sin(theta)
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta])
            continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt])
              continue;
            for (int k = 0; k < KIntK; ++k) {
              for (int m = 0; m < KIntM; ++m)
                sumpmwkEt[ieta][ipt][m][k] += std::pow(wPid, k) * std::pow(et, m);
              sumwkEt[ieta][ipt][k] += std::pow(wPid, k);
            }
          }
        }
      }
    }

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        const int ibx = pmeanNchEtabinPtbinStep2->GetXaxis()->FindBin(coll.multNTracksPV());
        const int iby = ieta + 1;
        const int ibz = ipt + 1;
        float mmpt = pmeanNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
        float mmet = pmeanEtNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);

        mean[ieta][ipt] = sumpmwk[ieta][ipt][1][1] / sumwk[ieta][ipt][1];
        meanEt[ieta][ipt] = sumpmwkEt[ieta][ipt][1][1] / sumwkEt[ieta][ipt][1];

        if (std::isfinite(mmpt)) {
          std::tie(mean[ieta][ipt], c2[ieta][ipt]) =
            calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwk[ieta][ipt], sumwk[ieta][ipt], mmpt);
          p1kBar[ieta][ipt] = mean[ieta][ipt] - mmpt;
        }
        if (std::isfinite(mmet)) {
          std::tie(meanEt[ieta][ipt], c2Et[ieta][ipt]) =
            calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwkEt[ieta][ipt], sumwkEt[ieta][ipt], mmet);
          p1kBarEt[ieta][ipt] = meanEt[ieta][ipt] - mmet;
        }
      }
    }

    if (std::isfinite(c2[0][0]))
      histos.fill(HIST("Prof_C2_Cent"), cent, c2[0][0]);
    if (std::isfinite(c2Et[0][0]))
      histos.fill(HIST("Prof_C2Et_Cent"), cent, c2Et[0][0]);
    if (std::isfinite(sumwk[0][0][1])) {
      histos.fill(HIST("Prof_MeanpT_Cent"), cent, mean[0][0]);
      histos.fill(HIST("Prof_MeanpT_Mult"), coll.multNTracksPV(), mean[0][0]);
    }
    if (std::isfinite(sumwkEt[0][0][1])) {
      histos.fill(HIST("Prof_MeanEt_Cent"), cent, meanEt[0][0]);
      histos.fill(HIST("Prof_MeanEt_Mult"), coll.multNTracksPV(), meanEt[0][0]);
    }

    for (int ieta = 0; ieta < KNEta; ++ieta) {
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        if (std::isfinite(c2[ieta][ipt]))
          histos.fill(HIST("Prof_C2_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, c2[ieta][ipt]);
        if (std::isfinite(c2Et[ieta][ipt]))
          histos.fill(HIST("Prof_C2Et_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, c2Et[ieta][ipt]);
      }
    }

    for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
      int ietaC = KNEta - ietaA;
      float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        float c2Sub = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
        float c2SubEt = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
        if (std::isfinite(c2Sub)) {
          switch (ipt) {
            case 0:
              histos.fill(HIST("Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
              break;
            case 1:
              histos.fill(HIST("Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
              break;
            case 2:
              histos.fill(HIST("Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
              break;
          }
          if (cent < KCentCovCut) {
            switch (ipt) {
              case 0:
                histos.fill(HIST("Prof_ipt0_Cov_Eta"), valy, c2Sub);
                break;
              case 1:
                histos.fill(HIST("Prof_ipt1_Cov_Eta"), valy, c2Sub);
                break;
              case 2:
                histos.fill(HIST("Prof_ipt2_Cov_Eta"), valy, c2Sub);
                break;
            }
          }
        }

        if (std::isfinite(c2SubEt)) {
          switch (ipt) {
            case 0:
              histos.fill(HIST("Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
              break;
            case 1:
              histos.fill(HIST("Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
              break;
            case 2:
              histos.fill(HIST("Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
              break;
          }
          if (cent < KCentCovCut) {
            switch (ipt) {
              case 0:
                histos.fill(HIST("Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
                break;
              case 1:
                histos.fill(HIST("Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
                break;
              case 2:
                histos.fill(HIST("Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
                break;
            }
          }
        }
      }
    }

    for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
      for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
        float valx = KHalf * (etaLw[ietaA] + etaUp[ietaA]);
        float valy = KHalf * (etaLw[ietaC] + etaUp[ietaC]);

        for (int ipt = 0; ipt < KNpT; ++ipt) {
          float covpt = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
          if (std::isfinite(covpt)) {
            switch (ipt) {
              case 0:
                histos.fill(HIST("Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);
                break;
              case 1:
                histos.fill(HIST("Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);
                break;
              case 2:
                histos.fill(HIST("Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);
                break;
            }
          }

          float covet = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
          if (std::isfinite(covet)) {
            switch (ipt) {
              case 0:
                histos.fill(HIST("Prof_ipt0_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
                break;
              case 1:
                histos.fill(HIST("Prof_ipt1_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
                break;
              case 2:
                histos.fill(HIST("Prof_ipt2_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
                break;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataFluc, "process data to calculate fluc pT and Et", cfgRunDataFluc);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<RadialFlowDecorr>(cfgc)};
  return workflow;
}
