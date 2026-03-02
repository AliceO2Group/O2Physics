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
  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr int KPiPlus = 211;
  static constexpr int KKPlus = 321;
  static constexpr int KProton = 2212;
  static constexpr int KNsp = 10;
  static constexpr float KBinOffset = 0.5f;
  static constexpr float KPhiMin = 0.f;
  static constexpr int KNbinsZvtx = 240;
  static constexpr float KZvtxMin = -12.f;
  static constexpr float KZvtxMax = 12.f;
  static constexpr float KPMin = 0.f;
  static constexpr float KPMax = 10.f;
  static constexpr int KNbinsPt = 200;
  static constexpr float KPtMin = 0.f;
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
  enum PID {
    kInclusive = 0, // Suffix ""
    kPiMinus,       // Suffix "_PiMinus"
    kPiPlus,        // Suffix "_PiPlus"
    kPiAll,         // Suffix "_PiAll"
    kKaMinus,       // Suffix "_KaMinus"
    kKaPlus,        // Suffix "_KaPlus"
    kKaAll,         // Suffix "_KaAll"
    kAntiProton,    // Suffix "_AntiPr"
    kProton,        // Suffix "_Pr"
    kAllProton,     // Suffix "_AllPr"
    kNumPID         // Total: 10
  };

  const std::vector<std::string> pidSuffix = {"", "_PiMinus", "_PiPlus", "_PiAll", "_KaMinus", "_KaPlus", "_KaAll", "_AntiPr", "_Pr", "_AllPr"};

  enum ECentralityEstimator {
    kCentFT0C = 1,
    kCentFV0A = 2
  };
  enum SystemType {
    kPbPb = 1,
    kNeNe = 2,
    kOO = 3,
    kpp = 4
  };
  static constexpr float KinvalidCentrality = -1.0f;
  const std::vector<float> etaLw = {
    -0.8,
    -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    const std::vector<float> etaUp = {
      0.8,
      -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

      Configurable<float> cfgVtxZCut{"cfgVtxZCut", 10.f, "z-vertex range"};
      Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "min pT"};
      Configurable<float> cfgPtMax{"cfgPtMax", 10.0f, "max pT"};
      Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "|η| cut"};
      Configurable<float> cfgDCAXY{"cfgDCAXY", 2.4f, "DCAxy cut"};
      Configurable<float> cfgDCAZ{"cfgDCAZ", 3.2f, "DCAz cut"};
      Configurable<float> cfgTPCClsMin{"cfgTPCClsMin", 70.f, "min TPC clusters"};
      Configurable<float> cfgChi2TPCMax{"cfgChi2TPCMax", 4.0f, "max TPC χ²"};
      Configurable<float> cfgCutTpcChi2NCl{"cfgCutTpcChi2NCl", 2.5f, "Maximum TPCchi2NCl"};
      Configurable<float> cfgCutItsChi2NCl{"cfgCutItsChi2NCl", 36.0f, "Maximum ITSchi2NCl"};

      Configurable<float> cfgPIDnSigmaCut{"cfgPIDnSigmaCut", 3.f, "TPC PID |nσ| cut"};

      Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
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
      Configurable<int> cfgNsubsample{"cfgNsubsample", 10, "Number of subsamples"};
      Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 1-->FT0C, 2-->FV0A"};
      Configurable<bool> cfgEvSelNoSameBunchPileup{"cfgEvSelNoSameBunchPileup", true, "Pileup removal"};
      Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
      Configurable<bool> cfgEvSelkNoITSROFrameBorder{"cfgEvSelkNoITSROFrameBorder", true, "ITSROFrame border event selection cut"};
      Configurable<bool> cfgEvSelkNoTimeFrameBorder{"cfgEvSelkNoTimeFrameBorder", true, "TimeFrame border event selection cut"};
      Configurable<bool> cfgIsGoodZvtxFT0VsPV{"cfgIsGoodZvtxFT0VsPV", true, "Good Vertexing cut"};

      Configurable<int> cfgNchPbMax{"cfgNchPbMax", 4000, "Max Nch range for PbPb collisions"};
      Configurable<int> cfgNchOMax{"cfgNchOMax", 600, "Max Nch range for OO collisions"};

      Configurable<int> cfgSys{"cfgSys", 1, "Efficiency to be used for which system? 1-->PbPb, 2-->NeNe, 3-->OO, 4-->pp"};
      Configurable<bool> cfgFlat{"cfgFlat", false, "Whether to use flattening weights"};
      Configurable<bool> cfgEff{"cfgEff", false, "Whether to use Efficiency weights"};
      Configurable<bool> cfgZDC{"cfgZDC", false, "Whether to use ZDC for pileup histograms"};

      Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
      Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

      ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 3.0, 5.0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "centrality axis (percentile)"};

      const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
      const AxisSpec centAxis1Per{101, -0.5, 100.5,"Centrality (%)"};
      AxisSpec nChAxis{1, 0., 1., "Nch", "Nch"};
      AxisSpec nChAxis2{1, 0., 1., "Nch", "Nch"};

      const AxisSpec vzAxis{5, -12.5, 12.5, "Vz"};
      const AxisSpec chgAxis{3, -1.5, 1.5};
      const AxisSpec pTAxis{{0.0, 0.2, 0.5, 1, 3, 5, 7.5, 10}, "pT Axis"};
      const AxisSpec etaAxis{{-0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9}, "Eta"};
      const AxisSpec phiAxis{KNbinsPhi, KPhiMin, TwoPI, "#phi"};
      const AxisSpec etaBinAxis{KNEta + 1, -KBinOffset, KNEta + KBinOffset, "#eta bin Number"};
      const AxisSpec gapAxis{{-1.55, -1.45, -1.35, -1.25, -1.15, -1.05, -0.95, -0.85, -0.75, -0.65, -0.55, -0.45, -0.35, -0.25, -0.15, -0.05, 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55}, "Gaps"};
      const AxisSpec sumAxis{{-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}, "Sums"};
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
        if (cfgEvSelNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
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
      // Filter trackFilter = nabs(aod::track::eta) < cfgEtaCut &&
      // aod::track::pt > cfgPtMin&&
      // aod::track::pt < cfgPtMax&&
      // nabs(aod::track::dcaXY) < cfgDCAXY&& nabs(aod::track::dcaZ) < cfgDCAZ;
      using AodTracksSel = soa::Filtered<UnfilteredTracks>;
      using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
      using FilteredTCs = soa::Filtered<TCs>;
      using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

      using MyRun3MCCollisions = soa::Join<
      aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
      aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As,
      aod::CentNGlobals, aod::McCollisionLabels>;

      // using MyMCTracks = soa::Join<
      // aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
      // aod::McTrackLabels,
      // aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
      // aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

      PresliceUnsorted<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
      PresliceUnsorted<MyRun3MCCollisions> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
      Preslice<FilteredTCs> trackPerCollision = aod::track::collisionId;

      void declareCommonQA()
      {
        histos.add("hZvtx_after_sel", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
        histos.add("hVtxZ", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
        histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis1Per}});
        histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis2}, {nChAxis2}});
        histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis2}, {centAxis1Per}});
        histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsPt, KPMin, KPMax}});
        histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
        histos.add("hEta", ";#eta", kTH1F, {{KNbinsEtaFine, KEtaMin, KEtaMax}});
        histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});
        histos.add("hEtaPhiReco", ";z_{vtx};sign;p_{T};#eta;#phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
      }

      void declareMCCommonHists()
      {
        for (const auto& suf : pidSuffix) {
          histos.add("h3_AllPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
          histos.add("h3_RecoMatchedToPrimary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
          histos.add("h3_AllReco" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
          histos.add("h3_RecoUnMatchedToPrimary_Secondary" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
          histos.add("h3_RecoUnMatchedToPrimary_Fake" + suf, ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, cfgPtMin, cfgPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
        }

        histos.add("ptResolution", ";p_{T}^{MC};(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}", kTH2F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}, {100, -0.2, 0.2}});
        histos.add("etaResolution", ";#eta^{MC};#eta^{reco}-#eta^{MC}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {100, -0.02, 0.02}});
        histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}});
        histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
        histos.add("vzResolution", ";Vz^{MC};(Vz^{reco}-Vz^{MC})/Vz^{MC}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {100, -0.1, 0.1}});

        histos.add("h_RecoUnMatchedToPrimary", ";p_{T}", kTH1F, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
      }
/*
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
        histos.add("Eff_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
        histos.add("Eff_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

        histos.add("Fake_cent", ";cent", kTProfile, {centAxis1Per});
        histos.add("Fake_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
        histos.add("Fake_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
        histos.add("Fake_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

        histos.add("wgt_cent", ";cent", kTProfile, {centAxis1Per});
        histos.add("wgt_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
        histos.add("wgt_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, cfgPtMin, cfgPtMax}});
        histos.add("wgt_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

        histos.add("pmeanFT0Amultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
        histos.add("pmeanFT0Cmultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
        histos.add("pmeanFT0A_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
        histos.add("pmeanFT0C_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
        histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;id;#eta", kTProfile3D, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});
        histos.add("h3_cent_id_eta_FT0", ";cent;id;#eta", kTH3F, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});

        // for (const char* lvl : {"MCGen/", "MCReco/", "MCRecoEffCorr/"}) {
        histos.add<TProfile2D>("MCGen/Prof_Cent_Nsp_Nchrec", ";cent;isp", kTProfile2D, {{centAxis1Per}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("MCGen/Prof_Mult_Nsp_Nchrec", ";mult;isp", kTProfile2D, {{nChAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("MCGen/Prof_Cent_Nsp_MeanpT", ";cent;isp", kTProfile2D, {{centAxis1Per}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("MCGen/Prof_Mult_Nsp_MeanpT", ";mult;isp", kTProfile2D, {{nChAxis}, {KNsp, -0.5, KNsp - 0.5}});
        // }

        histos.add<TProfile3D>("pmeanTru_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanReco_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanRecoEffcorr_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});

        histos.add<TProfile3D>("pmeanMultTru_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanMultReco_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanMultRecoEffcorr_nch_etabin_spbin", ";mult;eta;isp", kTProfile3D, {{nChAxis}, {KNEta, 0.5, KNEta + 0.5}, {KNsp, -0.5, KNsp - 0.5}});

        for (const auto& suf : pidSuffix) {
          histos.add<TProfile3D>("Prof2D_MeanpTSub_Tru" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
          histos.add<TProfile3D>("Prof2D_MeanpTSub_Reco" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
          histos.add<TProfile3D>("Prof2D_MeanpTSub_RecoEffCorr" + suf, ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
        }
      }

      void declareMCFlucHists()
      {
        histos.add<TProfile3D>("MCGen/Prof_MeanpT_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("MCGen/Prof_C2_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("MCGen/Prof_C2Sub_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("MCGen/Prof_Cov_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("MCGen/Prof_CovFT0A_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("MCGen/Prof_CovFT0C_Cent_etabin_spbin",";cent;eta;isp",kTProfile3D,{{centAxis1Per}, {etaBinAxis}, {KNsp, -0.5, KNsp - 0.5}});


        for (const auto& suf : pidSuffix) {
          histos.add<TProfile3D>(Form("MCGen/Prof_C2Sub2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
          histos.add<TProfile3D>(Form("MCGen/Prof_GapSum2D%s", suf.c_str()), ";cent;gap;sum", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
          histos.add<TProfile3D>(Form("MCGen/Prof_Cov2D_Cent_etaA_etaC%s", suf.c_str()), ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
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
        histos.add("hCentnTrk", ";Centrality (%);N_{trk}", kTH2F, {{centAxis1Per}, {nChAxis2}});
        histos.add("hCentnTrkPV", ";Centrality (%);N_{trk, PV}", kTH2F, {{centAxis1Per}, {nChAxis2}});
      }


      void declareDataMeanHists()
      {
        histos.add("pmeanFT0Amultpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
        histos.add("pmeanFT0A_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});
        histos.add("pmeanFT0Cmultpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
        histos.add("pmeanFT0C_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});

        histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;channel id; #eta;amplitude", kTProfile3D, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});
        histos.add("h3_cent_id_eta_FT0", ";cent;channel id; #eta", kTH3F, {{centAxis1Per}, {100, -0.5, 99.5}, {100, -5.0, 5.0}});

        histos.add<TProfile2D>("Prof_Cent_Nsp_Nchrec", ";cent;Species;#LT N_{PV}#GT", kTProfile2D, {{centAxis1Per}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("Prof_Mult_Nsp_Nchrec", ";N_{PV};Species;#LT N_{PV}#GT", kTProfile2D, {{nChAxis}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("Prof_Cent_Nsp_MeanpT", ";cent;Species;#LT p_{T}#GT", kTProfile2D, {{centAxis1Per}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile2D>("Prof_Mult_Nsp_MeanpT", ";N_{PV};Species;#LT p_{T}#GT", kTProfile2D, {{nChAxis}, {KNsp, -0.5, KNsp - 0.5}});

        histos.add<TProfile3D>("pmean_nch_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanMult_nch_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmean_cent_etabin_spbin", ";Centrality (%) ;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("pmeanMult_cent_etabin_spbin", ";Centrality (%) ;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});

        for (const auto& suf : pidSuffix) {
          histos.add("hEtaPhiReco" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
          histos.add("hEtaPhiRecoEffWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});
          histos.add("hEtaPhiRecoWtd" + suf, ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaAxis}, {phiAxis}});

          histos.add<TProfile3D>("Prof2D_MeanpTSub" + suf, ";cent;#eta_{A} bin;#eta_{C} bin", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}});

        }
      }

      void declareDataFlucHists()
      {
        histos.add<TProfile3D>("Prof_MeanpT_Cent_etabin_spbin", ";cent;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_MeanpT_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_C2_Cent_etabin_spbin", ";cent;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_C2_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});

        histos.add<TProfile3D>("Prof_C2Sub_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_C2Sub_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_Cov_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_Cov_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});

        histos.add<TProfile3D>("Prof_CovFT0A_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_CovFT0A_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_CovFT0C_Cent_etabin_spbin", ";Centrality;#eta-bin;Species", kTProfile3D, {{centAxis1Per}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});
        histos.add<TProfile3D>("Prof_CovFT0C_Mult_etabin_spbin", ";N_{PV};#eta-bin;Species", kTProfile3D, {{nChAxis}, {KNEta + 1, -KBinOffset, KNEta + KBinOffset}, {KNsp, -0.5, KNsp - 0.5}});

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
*/
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
          nChAxis = {cfgNchPbMax / 4, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
          nChAxis2 = {cfgNchPbMax / 20, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
        } else if (cfgSys == kNeNe || cfgSys == kOO) {
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

        std::string pathEff = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_EffMaps";
        std::string pathMCFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_MCFlatMaps";
        std::string pathMCMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_MCMean";
        std::string pathDataFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_DataFlatMaps";
        std::string pathDataMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_DataMean";

        declareCommonQA();

        if (cfgRunMCMean || cfgRunMCFluc || cfgRunGetEff) {
          declareMCCommonHists();
        }
/*
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
          histos.addClone("MCGen/", "MCReco/");
          histos.addClone("MCGen/", "MCRecoEffCorr/");
        }
        if (cfgRunDataMean) {
          declareDataMeanHists();
        }
        if (cfgRunDataFluc) {
          declareDataFlucHists();
        }
*/
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

          for (int i = 0; i < KNsp; ++i) {
            loadEffFakeForPID(static_cast<PID>(i));
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
                  hFlatWeight[i] = buildWeightMapFromRaw(hRaw, Form("hFlatWeight%s", suffix.c_str()));
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

              for (int i = 0; i < KNsp; ++i) {
                loadFlatForPID(static_cast<PID>(i));
              }
            } else {
              LOGF(error, "Could not retrieve MC Flattening TList from: %s", pathMCFlat.c_str());
            }
          }
        }

        auto loadTProfile3DFromList = [&](TList* sourceList, const char* objName, TProfile3D*& target) {
          if (!sourceList) return;
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
          if (!sourceList) return;
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
            loadTProfileFromList(lstMCMean, "pmeanFT0Amultpv", pmeanFT0AmultpvStep2);
            loadTProfileFromList(lstMCMean, "pmeanFT0Cmultpv", pmeanFT0CmultpvStep2);

            loadTProfile3DFromList(lstMCMean, "pmeanTru_nch_etabin_spbin", pmeanTruNchEtabinSpbinStep2);
            loadTProfile3DFromList(lstMCMean, "pmeanReco_nch_etabin_spbin", pmeanRecoNchEtabinSpbinStep2);
            loadTProfile3DFromList(lstMCMean, "pmeanRecoEffcorr_nch_etabin_spbin", pmeanRecoEffcorrNchEtabinSpbinStep2);

            loadTProfile3DFromList(lstMCMean, "pmeanMultTru_nch_etabin_spbin", pmeanMultTruNchEtabinSpbinStep2);
            loadTProfile3DFromList(lstMCMean, "pmeanMultReco_nch_etabin_spbin", pmeanMultRecoNchEtabinSpbinStep2);
            loadTProfile3DFromList(lstMCMean, "pmeanMultRecoEffcorr_nch_etabin_spbin", pmeanMultRecoEffcorrNchEtabinSpbinStep2);
          } else {
            LOGF(error, "Could not retrieve TList for MC Mean from: %s", pathMCMean.c_str());
          }
        }

        if (cfgRunDataFluc) {
          LOGF(info, "Loading Data Mean profiles from CCDB path: %s", pathDataMean.c_str());
          TList* lstDataMean = ccdb->getForTimeStamp<TList>(pathDataMean, now);

          if (lstDataMean) {
            loadTProfileFromList(lstDataMean, "pmeanFT0Amultpv", pmeanFT0AmultpvStep2);
            loadTProfileFromList(lstDataMean, "pmeanFT0Cmultpv", pmeanFT0CmultpvStep2);

            loadTProfile3DFromList(lstDataMean, "pmean_nch_etabin_spbin", pmeanNchEtabinSpbinStep2);
            loadTProfile3DFromList(lstDataMean, "pmeanMult_nch_etabin_spbin", pmeanMultNchEtabinSpbinStep2);
          } else {
            LOGF(error, "Could not retrieve TList for Data Mean from: %s", pathDataMean.c_str());
          }
        }
        LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
      }

      void processGetEffHists(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
      {
        /*
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
            histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, (float)tracks.size());
            histos.fill(HIST("Hist2D_cent_nch"), multPV, cent);

            for (const auto& particle : partSlice) {
              if (!isParticleSelected(particle) || !particle.isPhysicalPrimary())
              continue;

              const int pdgCode = particle.pdgCode();
              const int absPdg = std::abs(pdgCode);
              float pt = particle.pt();
              float eta = particle.eta();
              bool isSpecies[KNsp] = {
                true,                             // kInclusive
                pdgCode == kPiMinus,        // kPiMinus
                pdgCode == kPiPlus,        // kPiPlus
                absPdg == kPiPlus,        // kPiAll
                pdgCode == kKMinus,        // kKaMinus
                pdgCode == kKPlus,        // kKaPlus
                absPdg == kKPlus,        // kKaAll
                pdgCode == kProtonBar,          // kAntiProton
                pdgCode ==  kProton,          // kProton
                absPdg == kProton           // kAllProton
              };


              histos.fill(HIST("h3_AllPrimary"), multPV, pt, eta);

              if (isSpecies[kPiMinus])      histos.fill(HIST("h3_AllPrimary_PiMinus"), multPV, pt, eta);
              else if (isSpecies[kPiPlus])  histos.fill(HIST("h3_AllPrimary_PiPlus"),  multPV, pt, eta);
              if (isSpecies[kPiAll])        histos.fill(HIST("h3_AllPrimary_PiAll"),   multPV, pt, eta);

              if (isSpecies[kKaMinus])      histos.fill(HIST("h3_AllPrimary_KaMinus"), multPV, pt, eta);
              else if (isSpecies[kKaPlus])  histos.fill(HIST("h3_AllPrimary_KaPlus"),  multPV, pt, eta);
              if (isSpecies[kKaAll])        histos.fill(HIST("h3_AllPrimary_KaAll"),   multPV, pt, eta);

              if (isSpecies[kAntiProton])   histos.fill(HIST("h3_AllPrimary_AntiPr"),  multPV, pt, eta);
              else if (isSpecies[kProton])  histos.fill(HIST("h3_AllPrimary_Pr"),      multPV, pt, eta);
              if (isSpecies[kAllProton])    histos.fill(HIST("h3_AllPrimary_AllPr"),   multPV, pt, eta);
            }

            for (const auto& track : trackSlice) {
              if (!isTrackSelected(track))
              continue;

              float pt = track.pt();
              float eta = track.eta();
              float phi = track.phi();
              auto sign = track.sign();

              bool isPi = selectionPion(track);
              bool isKa = selectionKaon(track);
              bool isPr = selectionProton(track);

              bool isSpecies[KNsp] = {
                true,               // kInclusive
                isPi && sign < 0,   // kPiMinus
                isPi && sign > 0,   // kPiPlus
                isPi,               // kPiAll
                isKa && sign < 0,   // kKaMinus
                isKa && sign > 0,   // kKaPlus
                isKa,               // kKaAll
                isPr && sign < 0,   // kAntiProton
                isPr && sign > 0,   // kProton
                isPr                // kAllProton
              };

              for (int isp = 0; isp < KNsp; ++isp) {
                if (!isSpecies[isp])
                continue;

                if (isp == kInclusive) {
                  histos.fill(HIST("h3_AllReco"), multPV, pt, eta);
                  histos.fill(HIST("hEtaPhiReco"), col.posZ(), sign, pt, eta, phi);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();

                    if (mcP.isPhysicalPrimary()) {
                      histos.fill(HIST("ptResolution"), mcP.pt(), (pt - mcP.pt()) / mcP.pt());
                      histos.fill(HIST("etaResolution"), mcP.eta(), eta - mcP.eta());
                      histos.fill(HIST("etaTruthReco"), mcP.eta(), eta);
                      histos.fill(HIST("vzResolution"), mcP.vz(), (col.posZ() - mcP.vz()) / mcP.vz());
                      histos.fill(HIST("TruthTracKVz"), mcP.vz(), col.posZ());

                      histos.fill(HIST("h3_RecoMatchedToPrimary"), multPV, mcP.pt(), mcP.eta());
                    } else {
                      histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), multPV, pt, eta);
                      histos.fill(HIST("h_RecoUnMatchedToPrimary"), pt);
                    }
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), multPV, pt, eta);
                  }
                }
                else if (isp == kPiMinus) {
                  histos.fill(HIST("h3_AllReco_PiMinus"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kPiMinus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_PiMinus"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiMinus"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiMinus"), multPV, pt, eta);
                  }
                }
                else if (isp == kPiPlus) {
                  histos.fill(HIST("h3_AllReco_PiPlus"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kPiPlus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_PiPlus"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiPlus"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiPlus"), multPV, pt, eta);
                  }
                }
                else if (isp == kPiAll) {
                  histos.fill(HIST("h3_AllReco_PiAll"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && std::abs(mcP.pdgCode()) == kPiPlus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_PiAll"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PiAll"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_PiAll"), multPV, pt, eta);
                  }
                }
                else if (isp == kKaMinus) {
                  histos.fill(HIST("h3_AllReco_KaMinus"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kKMinus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_KaMinus"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaMinus"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaMinus"), multPV, pt, eta);
                  }
                }
                else if (isp == kKaPlus) {
                  histos.fill(HIST("h3_AllReco_KaPlus"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kKPlus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_KaPlus"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaPlus"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaPlus"), multPV, pt, eta);
                  }
                }
                else if (isp == kKaAll) {
                  histos.fill(HIST("h3_AllReco_KaAll"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && std::abs(mcP.pdgCode()) == kKPlus)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_KaAll"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_KaAll"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_KaAll"), multPV, pt, eta);
                  }
                }
                else if (isp == kAntiProton) {
                  histos.fill(HIST("h3_AllReco_AntiPr"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kProtonBar)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_AntiPr"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_AntiPr"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_AntiPr"), multPV, pt, eta);
                  }
                }
                else if (isp == kProton) {
                  histos.fill(HIST("h3_AllReco_Pr"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && mcP.pdgCode() == kProton)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_Pr"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_Pr"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_Pr"), multPV, pt, eta);
                  }
                }
                else if (isp == kAllProton) {
                  histos.fill(HIST("h3_AllReco_AllPr"), multPV, pt, eta);
                  if (track.has_mcParticle()) {
                    auto mcP = track.mcParticle();
                    if (mcP.isPhysicalPrimary() && std::abs(mcP.pdgCode()) == kProton)
                    histos.fill(HIST("h3_RecoMatchedToPrimary_AllPr"), multPV, mcP.pt(), mcP.eta());
                    else if (!mcP.isPhysicalPrimary())
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_AllPr"), multPV, pt, eta);
                  } else {
                    histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake_AllPr"), multPV, pt, eta);
                  }
                }
              } // end isp loop
            } // end track loop
          }
        }
        */
        LOGF(info, "FINISHED RUNNING processGetEffHists");
      }
      PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate EffWeights", cfgRunGetEff);
/*
      void processMCFlat(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& , aod::McParticles const& mcParticles)
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

              bool isPi = selectionPion(track);
              bool isKa = selectionKaon(track);
              bool isPr = selectionProton(track);

              bool isSpecies[KNsp] = {
                true,                      // kInclusive
                isPi && sign < 0,          // kPiMinus
                isPi && sign > 0,          // kPiPlus
                isPi,                      // kPiAll
                isKa && sign < 0,          // kKaMinus
                isKa && sign > 0,          // kKaPlus
                isKa,                      // kKaAll
                isPr && sign < 0,          // kAntiProton (Negative)
                isPr && sign > 0,          // kProton  (Positive)
                isPr                       // kAllProton
              };

              for (int isp = 0; isp < KNsp; ++isp) {
                if (!isSpecies[isp])
                continue;
                std::string suffix = pidSuffix[isp];

                float eff = getEfficiency(multPV, pt, eta, static_cast<PID>(isp), 0, cfgEff);
                float fake = getEfficiency(multPV, pt, eta, static_cast<PID>(isp), 1, cfgEff);
                float w = (1.0 - fake) / eff;

                if (std::isfinite(w) && w > 0.f && eff > KFloatEpsilon) {
                  if (isp == kInclusive) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kPiMinus) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kPiPlus) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kPiAll) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kKaMinus) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kKaPlus) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kKaAll) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kAntiProton) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kProton) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi, 1.0);
                  } else if (isp == kAllProton) {
                    histos.fill(HIST("hEtaPhiRecoEffWtd_AllPr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiReco_AllPr"), vz, sign, pt, eta, phi, 1.0);
                  }
                }
              }
            } // end track loop
          } // end col loop
        } // end mcColl loop
        LOGF(info, "FINISHED RUNNING processMCFlat");
      }
      PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

      void processMCMean(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& , aod::FT0s const&, aod::McParticles const& mcParticles)
      {
        // Track-sum arrays using KNsp index (isp=0: Incl, 1: Pi, 2: Ka, 3: Pr)
        double sumWiTruth[KNsp][KNEta]{}, sumWiptiTruth[KNsp][KNEta]{};
        double sumWiReco[KNsp][KNEta]{}, sumWiptiReco[KNsp][KNEta]{};
        double sumWiRecoEffCorr[KNsp][KNEta]{}, sumWiptiRecoEffCorr[KNsp][KNEta]{};

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
              if (pt <= KPtMin || pt > KPtMax)
              continue;
              int pdgCode = particle.pdgCode();
              int absPdg = std::abs(pdgCode);

              bool isSpecies[KNsp] = {
                true,                             // kInclusive
                pdgCode == kPiMinus,        // kPiMinus
                pdgCode == kPiPlus,        // kPiPlus
                absPdg == kPiPlus,        // kPiAll
                pdgCode == kKMinus,        // kKaMinus
                pdgCode == kKPlus,        // kKaPlus
                absPdg == kKPlus,        // kKaAll
                pdgCode == kProtonBar,          // kAntiProton
                pdgCode ==  kProton,          // kProton
                absPdg == kProton           // kAllProton
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

            // Full Event MCGen
            for (int isp = 0; isp < KNsp; ++isp) {
              histos.fill(HIST("MCGen/Prof_Cent_Nsp_Nchrec"), cent, isp, sumWiTruth[isp][0]);
              histos.fill(HIST("MCGen/Prof_Mult_Nsp_Nchrec"), multPV, isp, sumWiTruth[isp][0]);
              if (sumWiTruth[isp][0] > 1.0f) {
                histos.fill(HIST("MCGen/Prof_Cent_Nsp_MeanpT"), cent, isp, sumWiptiTruth[isp][0] / sumWiTruth[isp][0]);
                histos.fill(HIST("MCGen/Prof_Mult_Nsp_MeanpT"), multPV, isp, sumWiptiTruth[isp][0] / sumWiTruth[isp][0]);
              }
            }

            // --- 2. Reco Loop ---
            for (const auto& track : trackSlice) {
              if (!isTrackSelected(track))
              continue;
              float pt = track.pt(), eta = track.eta(), phi = track.phi();
              if (pt <= KPtMin || pt > KPtMax)
              continue;
              auto sign = track.sign();
              bool isPi = selectionPion(track);
              bool isKa = selectionKaon(track);
              bool isPr = selectionProton(track);

              bool isSpecies[KNsp] = {
                true,                      // kInclusive
                isPi && sign < 0,          // kPiMinus
                isPi && sign > 0,          // kPiPlus
                isPi,                      // kPiAll
                isKa && sign < 0,          // kKaMinus
                isKa && sign > 0,          // kKaPlus
                isKa,                      // kKaAll
                isPr && sign < 0,          // kAntiProton (Negative)
                isPr && sign > 0,          // kProton (Positive)
                isPr                       // kAllProton
              };

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
                  sumWiReco[isp][ieta]++;
                  sumWiptiReco[isp][ieta] += pt;
                  sumWiRecoEffCorr[isp][ieta] += w;
                  sumWiptiRecoEffCorr[isp][ieta] += w * pt;
                }

                if (isp == kInclusive) {
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

                if (isp == kInclusive) {
                  histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kPiMinus) {
                  histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kPiPlus) {
                  histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kPiAll) {
                  histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kKaMinus) {
                  histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kKaPlus) {
                  histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kKaAll) {
                  histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kProton) {
                  histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kAntiProton) {
                  histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                } else if (isp == kAllProton) {
                  histos.fill(HIST("hEtaPhiReco_AllPr"), vz, sign, pt, eta, phi);
                  histos.fill(HIST("hEtaPhiRecoWtd_AllPr"), vz, sign, pt, eta, phi, w);
                  histos.fill(HIST("hEtaPhiRecoEffWtd_AllPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
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

                  // 1. Truth Sub-event Mean
                  float nTruAB = sumWiTruth[isp][ietaA] + sumWiTruth[isp][ietaC];
                  float nRecoAB = sumWiReco[isp][ietaA] + sumWiReco[isp][ietaC];
                  float nCorrAB = sumWiRecoEffCorr[isp][ietaA] + sumWiRecoEffCorr[isp][ietaC];

                  float mptsubTru = (sumWiptiTruth[isp][ietaA] + sumWiptiTruth[isp][ietaC]) / nTruAB;
                  float mptsubReco = (sumWiptiReco[isp][ietaA] + sumWiptiReco[isp][ietaC]) / nRecoAB;
                  float mptsubRecoEffCorr = (sumWiptiRecoEffCorr[isp][ietaA] + sumWiptiRecoEffCorr[isp][ietaC]) / nCorrAB;

                  if (nTruAB > 0) {
                    if (isp == kInclusive)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kPiMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiMinus"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kPiPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiPlus"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kPiAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_PiAll"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kKaMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaMinus"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kKaPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaPlus"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kKaAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_KaAll"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_Pr"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kAntiProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_AntiPr"), cent, ietaA, ietaC, mptsubTru);
                    else if (isp == kAllProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Tru_AllPr"), cent, ietaA, ietaC, mptsubTru);
                  }

                  if (nRecoAB > 0) {
                    if (isp == kInclusive)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kPiMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiMinus"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kPiPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiPlus"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kPiAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_PiAll"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kKaMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaMinus"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kKaPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaPlus"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kKaAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_KaAll"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_Pr"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kAntiProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_AntiPr"), cent, ietaA, ietaC, mptsubReco);
                    else if (isp == kAllProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_Reco_AllPr"), cent, ietaA, ietaC, mptsubReco);
                  }

                  if (nCorrAB > 0) {
                    if (isp == kInclusive)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kPiMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiMinus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kPiPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiPlus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kPiAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_PiAll"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kKaMinus)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaMinus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kKaPlus)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaPlus"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kKaAll)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_KaAll"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_Pr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kAntiProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_AntiPr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                    else if (isp == kAllProton)
                    histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr_AllPr"), cent, ietaA, ietaC, mptsubRecoEffCorr);
                  }

                  if (ietaA == ietaC) { // only fill once per eta bin
                    if (sumWiTruth[isp][ietaA] > 0) {
                      float val = sumWiptiTruth[isp][ietaA] / sumWiTruth[isp][ietaA];
                      histos.fill(HIST("pmeanTru_nch_etabin_spbin"), multPV, ietaA, isp, val);
                    }
                    if (sumWiReco[isp][ietaA] > 0) {
                      float val = sumWiptiReco[isp][ietaA] / sumWiReco[isp][ietaA];
                      histos.fill(HIST("pmeanReco_nch_etabin_spbin"), multPV, ietaA, isp, val);
                    }
                    if (sumWiRecoEffCorr[isp][ietaA] > 0) {
                      float val = sumWiptiRecoEffCorr[isp][ietaA] / sumWiRecoEffCorr[isp][ietaA];
                      histos.fill(HIST("pmeanRecoEffcorr_nch_etabin_spbin"), multPV, ietaA, isp, val);
                      histos.fill(HIST("pmeanMultRecoEffcorr_nch_etabin_spbin"), multPV, ietaA, isp, sumWiRecoEffCorr[isp][ietaA]);
                    }
                  }
                } // end isp
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

            histos.fill(HIST("pmeanFT0Amultpv"), multPV, amplFT0A);
            histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
            histos.fill(HIST("pmeanFT0Cmultpv"), multPV, amplFT0C);
            histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
          }
        }
      }
      PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt and Eff Hists", cfgRunMCMean);

      void processMCFluc(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& , aod::FT0s const&, aod::McParticles const& mcParticles)
      {
        if (!pmeanTruNchEtabinSpbinStep2 || !pmeanRecoNchEtabinSpbinStep2 || !pmeanRecoEffcorrNchEtabinSpbinStep2 ||
          !pmeanMultTruNchEtabinSpbinStep2 || !pmeanMultRecoNchEtabinSpbinStep2 || !pmeanMultRecoEffcorrNchEtabinSpbinStep2) {
            LOGF(warning, "MC fluc: Unified Mean pT or Mult map missing");
            return;
          }

          // Expanded with KNsp index (isp=0: Incl, 1: Pi, 2: Ka, 3: Pr)
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
                if (pt <= KPtMin || pt > KPtMax)
                continue;
                float eta = particle.eta();
                int pdgCode = particle.pdgCode();
                int absPdg = std::abs(pdgCode);

                bool isSpecies[KNsp] = {
                  true,                             // kInclusive
                  pdgCode == kPiMinus,        // kPiMinus
                  pdgCode == kPiPlus,        // kPiPlus
                  absPdg == kPiPlus,        // kPiAll
                  pdgCode == kKMinus,        // kKaMinus
                  pdgCode == kKPlus,        // kKaPlus
                  absPdg == kKPlus,        // kKaAll
                  pdgCode == kProtonBar,          // kAntiProton
                  pdgCode ==  kProton,          // kProton
                  absPdg == kProton           // kAllProton
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

              // --- 2. Reco Loop ---
              float vz = col.posZ();
              for (const auto& track : trackSlice) {
                if (!isTrackSelected(track))
                continue;

                float pt = track.pt();
                if (pt <= KPtMin || pt > KPtMax)
                continue;
                float eta = track.eta();
                float phi = track.phi();
                auto sign = track.sign();
                bool isPi = selectionPion(track);
                bool isKa = selectionKaon(track);
                bool isPr = selectionProton(track);

                bool isSpecies[KNsp] = {
                  true,                      // kInclusive
                  isPi && sign < 0,          // kPiMinus
                  isPi && sign > 0,          // kPiPlus
                  isPi,                      // kPiAll
                  isKa && sign < 0,          // kKaMinus
                  isKa && sign > 0,          // kKaPlus
                  isKa,                      // kKaAll
                  isPr && sign < 0,          // kAntiProton (Negative)
                  isPr && sign > 0,          // kProton (Positive)
                  isPr                       // kAllProton
                };

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
                    for (int k = 0; k < KIntK; ++k) {
                      for (int m = 0; m < KIntM; ++m) {
                        sumPmwkReco[isp][ieta][m][k] += std::pow(1.0, k) * std::pow(pt, m);
                        sumPmwkRecoEffCor[isp][ieta][m][k] += std::pow(w, k) * std::pow(pt, m);
                      }
                      sumWkReco[isp][ieta][k] += std::pow(1.0, k);
                      sumWkRecoEffCor[isp][ieta][k] += std::pow(w, k);
                    }
                  }

                  if (isp == kInclusive) {
                    histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kPiMinus) {
                    histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kPiPlus) {
                    histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kPiAll) {
                    histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kKaMinus) {
                    histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kKaPlus) {
                    histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kKaAll) {
                    histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kProton) {
                    histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kAntiProton) {
                    histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  } else if (isp == kAllProton) {
                    histos.fill(HIST("hEtaPhiReco_AllPr"), vz, sign, pt, eta, phi);
                    histos.fill(HIST("hEtaPhiRecoWtd_AllPr"), vz, sign, pt, eta, phi, w);
                    histos.fill(HIST("hEtaPhiRecoEffWtd_AllPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                  }
                }
              } // trkslice

              for (int ieta = 0; ieta < KNEta; ++ieta) {
                const int ibx = pmeanTruNchEtabinSpbinStep2->GetXaxis()->FindBin(col.multNTracksPV());
                const int iby = ieta + 1;

                for (int isp = 0; isp < KNsp; ++isp) {
                  const int ibz = isp + 1;

                  meanTruMult[isp][ieta] = sumWkTru[isp][ieta][1];
                  meanRecoMult[isp][ieta] = sumWkReco[isp][ieta][1];
                  meanRecoEffCorMult[isp][ieta] = sumWkRecoEffCor[isp][ieta][1];

                  float mmptTru = pmeanTruNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
                  float mmptReco = pmeanRecoNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
                  float mmptRecoEffCor = pmeanRecoEffcorrNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

                  float mmMultTru = pmeanMultTruNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
                  float mmMultReco = pmeanMultRecoNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
                  float mmMultRecoEffCor = pmeanMultRecoEffcorrNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

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
              if (col.has_foundFT0()) {
                const auto& ft0 = col.foundFT0();
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
                    histos.fill(HIST("MCReco/Prof_Cent_NEta_Nsp_MeanpT"), cent, ieta, isp, meanRecoEffCor[isp][ieta]);
                    histos.fill(HIST("MCReco/Prof_Mult_NEta_Nsp_MeanpT"), multPV, ieta, isp, meanRecoEffCor[isp][ieta]);
                  }
                }
              }

              // --- 3. Fill 1D Profiles: Gen, Reco, and EffCorr Levels ---
              for (int ieta = 0; ieta < KNEta; ++ieta) {
                for (int isp = 0; isp < KNsp; ++isp) {
                  // --- MCGen (Truth) ---
                  if (std::isfinite(meanTru[isp][ieta])) {
                    histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanTru[isp][ieta]);
                    histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, meanTru[isp][ieta]);
                  }
                  if (std::isfinite(c2Tru[isp][ieta])) {
                    histos.fill(HIST("MCGen/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2Tru[isp][ieta]);
                    histos.fill(HIST("MCGen/Prof_C2_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, c2Tru[isp][ieta]);
                  }
                  // --- MCReco ---
                  if (std::isfinite(meanReco[isp][ieta])) {
                    histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanReco[isp][ieta]);
                    histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, meanReco[isp][ieta]);
                  }
                  if (std::isfinite(c2Reco[isp][ieta])) {
                    histos.fill(HIST("MCReco/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2Reco[isp][ieta]);
                    histos.fill(HIST("MCReco/Prof_C2_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, c2Reco[isp][ieta]);
                  }
                  // --- MCRecoEffCorr ---
                  if (std::isfinite(meanRecoEffCor[isp][ieta])) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin_spbin"), cent, ieta, isp, meanRecoEffCor[isp][ieta]);
                    histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, meanRecoEffCor[isp][ieta]);
                  }
                  if (std::isfinite(c2RecoEffCor[isp][ieta])) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin_spbin"), cent, ieta, isp, c2RecoEffCor[isp][ieta]);
                    histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin_spbin"), col.multNTracksPV(), ieta, isp, c2RecoEffCor[isp][ieta]);
                  }
                }
              }

              p1kBarFt0A = amplFT0A - pmeanFT0AmultpvStep2->GetBinContent(pmeanFT0AmultpvStep2->GetXaxis()->FindBin(col.multNTracksPV()));
              p1kBarFt0C = amplFT0C - pmeanFT0CmultpvStep2->GetBinContent(pmeanFT0CmultpvStep2->GetXaxis()->FindBin(col.multNTracksPV()));

              // --- 4. Symmetric Sub-Event (1D) Covariances ---
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
                    histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, c2SubTru);
                  }
                  if (std::isfinite(c2SubReco)) {
                    histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2SubReco);
                    histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, c2SubReco);
                  }
                  if (std::isfinite(c2SubRecoEffCor)) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin_spbin"), cent, ietaA, isp, c2SubRecoEffCor);
                    histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, c2SubRecoEffCor);
                  }
                  if (std::isfinite(covTru)) {
                    histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covTru);
                    histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covTru);
                  }
                  if (std::isfinite(covReco)) {
                    histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covReco);
                    histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covReco);
                  }
                  if (std::isfinite(covRecoEffCor)) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin_spbin"), cent, ietaA, isp, covRecoEffCor);
                    histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covRecoEffCor);
                  }

                  if (std::isfinite(covFT0ATru)) {
                    histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0ATru);
                    histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0ATru);
                  }
                  if (std::isfinite(covFT0AReco)) {
                    histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0AReco);
                    histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0AReco);
                  }
                  if (std::isfinite(covFT0ARecoEffCor)) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin_spbin"), cent, ietaA, isp, covFT0ARecoEffCor);
                    histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0ARecoEffCor);
                  }

                  if (std::isfinite(covFT0CTru)) {
                    histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CTru);
                    histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0CTru);
                  }
                  if (std::isfinite(covFT0CReco)) {
                    histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CReco);
                    histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0CReco);
                  }
                  if (std::isfinite(covFT0CRecoEffCor)) {
                    histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin_spbin"), cent, ietaA, isp, covFT0CRecoEffCor);
                    histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin_spbin"), col.multNTracksPV(), ietaA, isp, covFT0CRecoEffCor);
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

                    if (isp == kInclusive) {
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

                    } else if (isp == kPiMinus) {
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

                    } else if (isp == kPiPlus) {
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

                    } else if (isp == kPiAll) {
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

                    } else if (isp == kKaMinus) {
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

                    } else if (isp == kKaPlus) {
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

                    } else if (isp == kKaAll) {
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

                    } else if (isp == kProton) {
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

                    } else if (isp == kAntiProton) {
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

                    } else if (isp == kAllProton) {
                      if (std::isfinite(c2SubTru)) {
                        histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, c2SubTru);
                        histos.fill(HIST("MCGen/Prof_GapSum2D_AllPr"), cent, gap, sum, c2SubTru);
                      }
                      if (std::isfinite(c2SubReco)) {
                        histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, c2SubReco);
                        histos.fill(HIST("MCReco/Prof_GapSum2D_AllPr"), cent, gap, sum, c2SubReco);
                      }
                      if (std::isfinite(c2SubRecoEffCor)) {
                        histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, c2SubRecoEffCor);
                        histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D_AllPr"), cent, gap, sum, c2SubRecoEffCor);
                      }

                      if (std::isfinite(covTru))
                      histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covTru);
                      if (std::isfinite(covReco))
                      histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covReco);
                      if (std::isfinite(covRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covRecoEffCor);

                      if (std::isfinite(covFT0ATru))
                      histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0ATru);
                      if (std::isfinite(covFT0AReco))
                      histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0AReco);
                      if (std::isfinite(covFT0ARecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0ARecoEffCor);

                      if (std::isfinite(covFT0CTru))
                      histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0CTru);
                      if (std::isfinite(covFT0CReco))
                      histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0CReco);
                      if (std::isfinite(covFT0CRecoEffCor))
                      histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0CRecoEffCor);
                    }
                  }
                }
              }
            } // colSlice
          } // mcColl
          LOGF(info, "FINISHED RUNNING processMCFluc");
        }
        PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt fluc", cfgRunMCFluc);

        void processGetDataFlat(AodCollisionsSel::iterator const& coll, BCsRun3 const& , aod::Zdcs const& , AodTracksSel const& tracks)
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
            float pt = track.pt();
            if (pt <= KPtMin || pt > KPtMax)
            continue;
            float eta = track.eta();
            float phi = track.phi();
            auto sign = track.sign();
            if (eta > etaLw[0] && eta < etaUp[0])
            ntrk++;
            bool isPi = selectionPion(track);
            bool isKa = selectionKaon(track);
            bool isPr = selectionProton(track);

            bool isSpecies[KNsp] = {
              true,                      // kInclusive
              isPi && sign < 0,          // kPiMinus
              isPi && sign > 0,          // kPiPlus
              isPi,                      // kPiAll
              isKa && sign < 0,          // kKaMinus
              isKa && sign > 0,          // kKaPlus
              isKa,                      // kKaAll
              isPr && sign < 0,          // kAntiProton (Negative)
              isPr && sign > 0,          // kProton (Positive)
              isPr                       // kAllProton
            };

            for (int isp = 0; isp < KNsp; ++isp) {
              if (!isSpecies[isp])
              continue;

              // Fetch efficiency specifically for this particle species
              float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
              float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
              float w = (1.0 - fake) / eff;

              if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
              continue;
              if (isp == kInclusive) {
                histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiMinus) {
                histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiPlus) {
                histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiAll) {
                histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaMinus) {
                histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaPlus) {
                histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaAll) {
                histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
              } else if (isp == kProton) {
                histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
              } else if (isp == kAntiProton) {
                histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
              } else if (isp == kAllProton) {
                histos.fill(HIST("hEtaPhiReco_AllPr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_AllPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_AllPr"), vz, sign, pt, eta, phi, w);
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

        void processDataMean(AodCollisionsSel::iterator const& coll, BCsRun3 const&, aod::Zdcs const& , aod::FT0s const&, AodTracksSel const& tracks)
        {
          double sumWi[KNsp][KNEta]{}, sumWipti[KNsp][KNEta]{};

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

            float p = track.p();
            float pt = track.pt();
            float eta = track.eta();
            float phi = track.phi();
            auto sign = track.sign();

            if (p < KFloatEpsilon)
            continue;

            if (pt <= KPtMin || pt > KPtMax)
            continue;

            histos.fill(HIST("hP"), p);
            histos.fill(HIST("hPt"), pt);
            histos.fill(HIST("hEta"), eta);
            histos.fill(HIST("hPhi"), phi);

            bool isPi = selectionPion(track);
            bool isKa = selectionKaon(track);
            bool isPr = selectionProton(track);

            bool isSpecies[KNsp] = {
              true,                      // kInclusive
              isPi && sign < 0,          // kPiMinus
              isPi && sign > 0,          // kPiPlus
              isPi,                      // kPiAll
              isKa && sign < 0,          // kKaMinus
              isKa && sign > 0,          // kKaPlus
              isKa,                      // kKaAll
              isPr && sign < 0,          // kAntiProton (Negative)
              isPr && sign > 0,          // kProton (Positive)
              isPr                       // kAllProton
            };

            for (int isp = 0; isp < KNsp; ++isp) {
              if (!isSpecies[isp])
              continue;

              float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
              float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
              float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
              float w = flatWeight * (1.0 - fake) / eff;

              if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
              continue;

              if (isp == kInclusive) {
                histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiMinus) {
                histos.fill(HIST("hEtaPhiReco_PiMinus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiMinus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiPlus) {
                histos.fill(HIST("hEtaPhiReco_PiPlus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiPlus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kPiAll) {
                histos.fill(HIST("hEtaPhiReco_PiAll"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_PiAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_PiAll"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaMinus) {
                histos.fill(HIST("hEtaPhiReco_KaMinus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaMinus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaMinus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaPlus) {
                histos.fill(HIST("hEtaPhiReco_KaPlus"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaPlus"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaPlus"), vz, sign, pt, eta, phi, w);
              } else if (isp == kKaAll) {
                histos.fill(HIST("hEtaPhiReco_KaAll"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_KaAll"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_KaAll"), vz, sign, pt, eta, phi, w);
              } else if (isp == kProton) {
                histos.fill(HIST("hEtaPhiReco_Pr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_Pr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_Pr"), vz, sign, pt, eta, phi, w);
              } else if (isp == kAntiProton) {
                histos.fill(HIST("hEtaPhiReco_AntiPr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_AntiPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_AntiPr"), vz, sign, pt, eta, phi, w);
              } else if (isp == kAllProton) {
                histos.fill(HIST("hEtaPhiReco_AllPr"), vz, sign, pt, eta, phi);
                histos.fill(HIST("hEtaPhiRecoEffWtd_AllPr"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
                histos.fill(HIST("hEtaPhiRecoWtd_AllPr"), vz, sign, pt, eta, phi, w);
              }

              // Accumulate sum
              for (int ieta = 0; ieta < KNEta; ++ieta) {
                if (eta <= etaLw[ieta] || eta > etaUp[ieta])
                continue;
                sumWi[isp][ieta] += w;
                sumWipti[isp][ieta] += w * pt;
              }
            }
          }

          // Full Event Means
          for (int isp = 0; isp < KNsp; ++isp) {
            histos.fill(HIST("Prof_Cent_Nsp_Nchrec"), cent, isp, sumWi[isp][0]);
            histos.fill(HIST("Prof_Mult_Nsp_Nchrec"), coll.multNTracksPV(), isp, sumWi[isp][0]);
            if (sumWi[isp][0] > 1.0f)
            histos.fill(HIST("Prof_Cent_Nsp_MeanpT"), cent, isp, sumWipti[isp][0] / sumWi[isp][0]);
          }

          // Kinematic Bin Means (1D and 2D Sub-event)
          for (int ietaA = 0; ietaA < KNEta; ++ietaA) {
            for (int ietaC = 0; ietaC < KNEta; ++ietaC) {
              for (int isp = 0; isp < KNsp; ++isp) {
                double wCorrAB = sumWi[isp][ietaA] + sumWi[isp][ietaC];
                if (wCorrAB > 0) {
                  float mptsub = (sumWipti[isp][ietaA] + sumWipti[isp][ietaC]) / wCorrAB;
                  if (isp == kInclusive)
                  histos.fill(HIST("Prof2D_MeanpTSub"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kPiMinus)
                  histos.fill(HIST("Prof2D_MeanpTSub_PiMinus"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kPiPlus)
                  histos.fill(HIST("Prof2D_MeanpTSub_PiPlus"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kPiAll)
                  histos.fill(HIST("Prof2D_MeanpTSub_PiAll"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kKaMinus)
                  histos.fill(HIST("Prof2D_MeanpTSub_KaMinus"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kKaPlus)
                  histos.fill(HIST("Prof2D_MeanpTSub_KaPlus"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kKaAll)
                  histos.fill(HIST("Prof2D_MeanpTSub_KaAll"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kProton)
                  histos.fill(HIST("Prof2D_MeanpTSub_Pr"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kAntiProton)
                  histos.fill(HIST("Prof2D_MeanpTSub_AntiPr"), cent, ietaA, ietaC, mptsub);
                  else if (isp == kAllProton)
                  histos.fill(HIST("Prof2D_MeanpTSub_AllPr"), cent, ietaA, ietaC, mptsub);
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

          histos.fill(HIST("pmeanFT0Amultpv"), coll.multNTracksPV(), amplFT0A);
          histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
          histos.fill(HIST("pmeanFT0Cmultpv"), coll.multNTracksPV(), amplFT0C);
          histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
        }
        PROCESS_SWITCH(RadialFlowDecorr, processDataMean, "process data to calculate mean pT", cfgRunDataMean);

        void processDataFluc(AodCollisionsSel::iterator const& coll, BCsRun3 const& , aod::Zdcs const& , aod::FT0s const&, AodTracksSel const& tracks)
        {
          if (!isEventSelected(coll))
          return;
          float cent = getCentrality(coll);
          if (cent > KCentMax)
          return;

          if (!pmeanNchEtabinSpbinStep2 || !pmeanMultNchEtabinSpbinStep2) {
            LOGF(warning, "Data fluc: Mean pT or Mult map missing");
            return;
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

            if (pt <= KPtMin || pt > KPtMax)
            continue;

            bool isPi = selectionPion(track);
            bool isKa = selectionKaon(track);
            bool isPr = selectionProton(track);

            bool isSpecies[KNsp] = {
              true,                      // kInclusive
              isPi && sign < 0,          // kPiMinus
              isPi && sign > 0,          // kPiPlus
              isPi,                      // kPiAll
              isKa && sign < 0,          // kKaMinus
              isKa && sign > 0,          // kKaPlus
              isKa,                      // kKaAll
              isPr && sign < 0,          // kAntiProton (Negative)
              isPr && sign > 0,          // kProton (Positive)
              isPr                       // kAllProton
            };

            for (int isp = 0; isp < KNsp; ++isp) {
              if (!isSpecies[isp])
              continue;

              float eff = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 0, cfgEff);
              float fake = getEfficiency(coll.multNTracksPV(), pt, eta, static_cast<PID>(isp), 1, cfgEff);
              float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, static_cast<PID>(isp), cfgFlat);
              float w = flatWeight * (1.0 - fake) / eff;

              if (!std::isfinite(w) || w <= KFloatEpsilon || eff <= KFloatEpsilon)
              continue;

              // Kinematic Bin sum
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
          double p1kBarFt0A = amplFT0A - pmeanFT0AmultpvStep2->GetBinContent(pmeanFT0AmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));
          double p1kBarFt0C = amplFT0C - pmeanFT0CmultpvStep2->GetBinContent(pmeanFT0CmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));

          // --- 2. Step 2 Means and 1D Fluc Variables ---
          for (int ieta = 0; ieta < KNEta; ++ieta) {
            const int ibx = pmeanNchEtabinSpbinStep2->GetXaxis()->FindBin(coll.multNTracksPV());
            const int iby = ieta + 1;

            for (int isp = 0; isp < KNsp; ++isp) {
              const int ibz = isp + 1;

              float mmpt = pmeanNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);
              float mmMult = pmeanMultNchEtabinSpbinStep2->GetBinContent(ibx, iby, ibz);

              mean[isp][ieta] = sumpmwk[isp][ieta][1][1] / sumwk[isp][ieta][1];
              meanMult[isp][ieta] = sumwk[isp][ieta][1];

              if (std::isfinite(mmpt)) {
                std::tie(mean[isp][ieta], c2[isp][ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwk[isp][ieta], sumwk[isp][ieta], mmpt);
                p1kBar[isp][ieta] = mean[isp][ieta] - mmpt;
              }
              p1kBarMult[isp][ieta] = meanMult[isp][ieta] - mmMult;
            }
          }

          // --- 3. Fill 1D Profiles ---
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

          // --- 5. Full 2D Covariances & GapSum Profiles ---
          for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
            for (int ietaC = 1; ietaC < KNEta; ++ietaC) {

              float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
              float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
              float gap = etaValA - etaValB;
              float sum = (etaValA + etaValB) / 2.0f;
              for (int isp = 0; isp < KNsp; ++isp) {

                float c2Sub = p1kBar[isp][ietaA] * p1kBar[isp][ietaC];
                float cov = p1kBarMult[isp][ietaA] * p1kBar[isp][ietaC];
                float covFT0A = p1kBarFt0A * p1kBar[isp][ietaC];
                float covFT0C = p1kBarFt0C * p1kBar[isp][ietaA];

                if (isp == kInclusive) {
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
                } else if (isp == kPiMinus) {
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
                } else if (isp == kPiPlus) {
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
                } else if (isp == kPiAll) {
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
                } else if (isp == kKaMinus) {
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
                } else if (isp == kKaPlus) {
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
                } else if (isp == kKaAll) {
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
                } else if (isp == kProton) {
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
                } else if (isp == kAntiProton) {
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
                } else if (isp == kAllProton) {
                  if (std::isfinite(c2Sub)) {
                    histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, c2Sub);
                    histos.fill(HIST("Prof_GapSum2D_AllPr"), cent, gap, sum, c2Sub);
                  }
                  if (std::isfinite(cov))
                  histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, cov);
                  if (std::isfinite(covFT0A))
                  histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0A);
                  if (std::isfinite(covFT0C))
                  histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC_AllPr"), cent, etaValA, etaValB, covFT0C);
                }
              }
            }
          }
        }
        PROCESS_SWITCH(RadialFlowDecorr, processDataFluc, "process data to calculate fluc pT", cfgRunDataFluc);
*/
      };

      WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
      {
        WorkflowSpec workflow{adaptAnalysisTask<RadialFlowDecorr>(cfgc)};
        return workflow;
      }
