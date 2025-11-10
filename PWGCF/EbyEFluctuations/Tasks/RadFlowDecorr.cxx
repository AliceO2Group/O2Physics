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
/// \file RadFlowDecorr.cxx
/// \brief Task for Radial Flow Decorrelation Measurement
/// \author Somadutta Bhatta
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

#include "TDirectory.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/Configurable.h"
#include "Framework/Logger.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "DataFormatsParameters/GRPObject.h"
#include "CCDB/BasicCCDBManager.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackTPCITS.h"
#include "MathUtils/Utils.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/Centrality.h"
#include <cmath>
#include "Framework/Logger.h"
#include "Framework/HistogramRegistry.h"
#include "CommonConstants/MathConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct RadFlowDecorr {

  static constexpr int KIntM = 3;
  static constexpr int KIntK = 3;
  // static constexpr int KNEta = 33;
  static constexpr int KNEta = 17;
  static constexpr int KNpT = 3;

  // --- Linter Fixes: Magic Numbers ---
  static constexpr float kFloatEpsilon = 1e-6f;
//   static constexpr int kPiPlus = o2::constants::physics::Pdg::kPiPlus;
//   static constexpr int kKPlus = o2::constants::physics::Pdg::kKaPlus;
//   static constexpr int kProton = o2::constants::physics::Pdg::kProton;
  static constexpr float kCentTestMin = 10.f;
  static constexpr float kCentTestMaxLo = 60.f;
  static constexpr float kCentTestMaxHi = 70.f;
  static constexpr float kCentCovCut = 1.0f;
  static constexpr float kBinOffset = 0.5f;
  static constexpr float kHalf = 0.5f;
  // Histogram axis definitions
  static constexpr int kNbinsZvtx = 240;
  static constexpr float kZvtxMin = -12.f;
  static constexpr float kZvtxMax = 12.f;
  static constexpr int kNbinsP = 100;
  static constexpr float kPMin = 0.f;
  static constexpr float kPMax = 10.f;
  static constexpr int kNbinsPt = 200;
  static constexpr float kPtMin = 0.f;
  static constexpr float kPtMax = 10.f;
  static constexpr int kNbinsEta = 120;
  static constexpr float kEtaMin = -1.2f;
  static constexpr float kEtaMax = 1.2f;
  static constexpr int kNbinsPhi = 64;
  static constexpr float kEtaAxisMin = -0.8f;
  static constexpr float kEtaAxisMax = 0.8f;
  static constexpr int kNbinsPhiFine = 30;
  static constexpr int kNbinsPtRes = 50;
  static constexpr float kPtResMax = 1.f;
  static constexpr int kNbinsEtaRes = 100;
  static constexpr float kEtaResMax = 0.5f;
  static constexpr int kNbinsVz = 80;
  static constexpr float kVzMin = -40.f;
  static constexpr float kVzMax = 40.f;
  static constexpr float kVzResMax = 20.f;
  static constexpr int kNbinsEtaFine = 20;
  static constexpr float kEtaFineMax = 1.f;
  static constexpr int kNbinsDca = 400;
  static constexpr float kDcaMax = 0.2f;
  static constexpr int kNbinsPtCoarse = 50;

  enum PID { kInclusive = 0, kCombinedPID, kNumPID };
  const std::vector<std::string> pidSuffix = {"", "_PID"};

  const std::vector<float> etaLw = {
    -0.8,
    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0.0, 0.1,0.2, 0.3, 0.4, 0.5, 0.6,0.7};
    const std::vector<float> etaUp = {
      0.8,
      -0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

      const std::vector<float> pTLw = {0.2, 0.2, 0.2};
      const std::vector<float> pTUp = {3,   5.0, 10.0};
      //==============================
      // Configurables
      //==============================
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
      Configurable<float> cfgCutTrackDcaZ{"cfgCutTrackDcaZ", 2.0f, "Maximum DcaZ"};
      Configurable<float> cfgCutTrackDcaXY{"cfgCutTrackDcaXY", 0.2f, "Maximum DcaZ"};

      Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
      Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
      Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};
      Configurable<float> cfgCutPtUpperTPC{"cfgCutPtUpperTPC", 0.6f, "Upper pT cut for PID using TPC only"};
      Configurable<float> cfgnSigmaOtherParticles{"cfgnSigmaOtherParticles", 3.0f, "PID nSigma cut to remove other particles (default:3)"};
      Configurable<float> cfgnSigmaCutTPC{"cfgnSigmaCutTPC", 2.0f, "PID nSigma cut for TPC"};
      Configurable<float> cfgnSigmaCutTOF{"cfgnSigmaCutTOF", 2.0f, "PID nSigma cut for TOF"};
      Configurable<float> cfgnSigmaCutCombTPCTOF{"cfgnSigmaCutCombTPCTOF", 2.0f, "PID nSigma combined cut for TPC and TOF"};
      Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut"};
      Configurable<float> cfgCutPtLowerProt{"cfgCutPtLowerProt", 0.2f, "Lower pT cut"};
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


      //==============================
      // Services and registry
      //==============================
      Service<ccdb::BasicCCDBManager> ccdb;
      Service<o2::framework::O2DatabasePDG> pdg;
      HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
      //==============================
      // External objects (shared across processes)
      //==============================
      // ### MODIFICATION ###
      // Store maps in arrays, one for each PID type (kInclusive, kCombinedPID)
      std::array<TH3F*, kNumPID> hEff{};
      std::array<TH3F*, kNumPID> hFake{};
      std::array<TH3F*, kNumPID> hWeightMap3D{}; // (cent, eta, phi)
      // ####################
      // mean pT / Et profiles loaded in init()
      // --- pT Profiles ---
      TProfile3D* pmeanTruNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanRecoNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanRecoMatchedNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanRecoEffcorrNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanRecoMatchedEffcorrNchEtabinPtbinStep2 = nullptr;
      // --- Et Profiles ---
      TProfile3D* pmeanEtTruNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanEtRecoNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanEtRecoMatchedNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanEtRecoEffcorrNchEtabinPtbinStep2 = nullptr; // This one was already here
      TProfile3D* pmeanEtRecoMatchedEffcorrNchEtabinPtbinStep2 = nullptr;
      // --- Data Profiles ---
      TProfile3D* pmeanNchEtabinPtbinStep2 = nullptr;
      TProfile3D* pmeanEtNchEtabinPtbinStep2 = nullptr;
      //==============================
      // Basic helpers
      //==============================
      // Line 131
      template <typename T> bool isEventSelected(const T& col) {
        if (!col.sel8()) return false;
        if (std::abs(col.posZ()) > cfgCutVertex) return false;
        if (cfgEvSelkNoSameBunchPileup   && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) return false;
        if (cfgEvSelkNoITSROFrameBorder  && !col.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) return false;
        if (cfgEvSelkNoTimeFrameBorder   && !col.selection_bit(o2::aod::evsel::kNoTimeFrameBorder))  return false;
        return true;
      }


      // Line 137
      template <typename T> bool isTrackSelected(const T& trk) {
        if (trk.sign() == 0) return false;
        if (!trk.has_collision()) return false;
        if (!trk.isPVContributor()) return false;
        if (!(trk.itsNCls() > cfgITScluster)) return false;
        if (!(trk.tpcNClsFound() >= cfgTPCcluster)) return false;
        if (!(trk.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) return false;

        if (trk.pt() < cfgCutPtLower || trk.pt() > cfgCutPtUpper || std::abs(trk.eta()) > cfgCutEta) return false;
        if (std::abs(trk.dcaXY()) > cfgCutTrackDcaXY || std::abs(trk.dcaZ()) > cfgCutTrackDcaZ) return false;
        return true;
      }

      template <typename T>
      bool isParticleSelected(const T& particle)
      {
        auto* pd = pdg->GetParticle(particle.pdgCode());
        if (!pd) return false;
        // if (dpt::isStrangeBaryonPDG(particle.pdgCode())) return false;
        if (std::abs(pd->Charge()) ==0) return false;
        if (particle.pt() < cfgCutPtLower || particle.pt() > cfgCutPtUpper || std::abs(particle.eta()) > cfgCutEta) return false;
        if (std::abs(particle.vz()) > cfgCutVertex) return false;
        return true;
      }

      //
      // // Lines 153-155
      // template <typename T> bool selectionPion (const T& trk) { return std::abs(trk.tpcNSigmaPi()) < cfgPIDnSigmaCut; }
      // template <typename T> bool selectionKaon (const T& trk) { return std::abs(trk.tpcNSigmaKa()) < cfgPIDnSigmaCut; }
      // template <typename T> bool selectionProton(const T& trk) { return std::abs(trk.tpcNSigmaPr()) < cfgPIDnSigmaCut; }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      template <typename T>
      bool selectionProton(const T& candidate)
      {
        if (!candidate.hasTPC())
        return false;
        int flag = 0; //! pid check main flag

        if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
          if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC) {
            flag = 1;
          }
          if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < cfgnSigmaCutTOF) {
            flag = 1;
          }
        }
        if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
          float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
          float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
          float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

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
        int flag = 0; //! pid check main flag

        if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
          if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC) {
            flag = 1;
          }
          if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < cfgnSigmaCutTOF) {
            flag = 1;
          }
        }
        if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
          float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
          float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
          float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

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
        int flag = 0; //! pid check main flag

        if (candidate.pt() > cfgCutPtLower && candidate.pt() <= cfgCutPtUpperTPC) {
          if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC) {
            flag = 1;
          }
          if (candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgnSigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < cfgnSigmaCutTOF) {
            flag = 1;
          }
        }
        if (candidate.hasTOF() && candidate.pt() > cfgCutPtUpperTPC && candidate.pt() < cfgCutPtUpperPID) {
          float combNSigmaPr = std::sqrt(std::pow(candidate.tpcNSigmaPr(), 2.0) + std::pow(candidate.tofNSigmaPr(), 2.0));
          float combNSigmaPi = std::sqrt(std::pow(candidate.tpcNSigmaPi(), 2.0) + std::pow(candidate.tofNSigmaPi(), 2.0));
          float combNSigmaKa = std::sqrt(std::pow(candidate.tpcNSigmaKa(), 2.0) + std::pow(candidate.tofNSigmaKa(), 2.0));

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


      // Utility helpers
      //==============================
      float getCentrality(const auto& col) const {
        if (cfgCentralityChoice.value == 1) return col.centFT0C();
        if (cfgCentralityChoice.value == 2) return col.centFT0A();
        if (cfgCentralityChoice.value == 3) return col.centFT0M();
        if (cfgCentralityChoice.value == 4) return col.centFV0A();
        return -1.0f;
      }
      //
      float getEfficiency(float mult, float pt, float eta, PID pidType, int effidx) const {
        TH3F* h;
        if(effidx==0) h = hEff[pidType];
        if(effidx==1) h = hFake[pidType];

        if (!h) return -1;
        const int ibx = h->GetXaxis()->FindBin(mult);
        const int iby = h->GetYaxis()->FindBin(pt);
        const int ibz = h->GetZaxis()->FindBin(eta);
        float val = h->GetBinContent(ibx, iby, ibz);
        return val;
      }

      /*
      float getEfficiency(float mult, float pt, float eta, PID pidType, int effidx) const {
      TH3F* h;
      if(effidx==0) h = hEff[pidType];  // Efficiency map
      if(effidx==1) h = hFake[pidType]; // Fake map

      if (!h) {
      LOGF(error, "getEfficiency: Histogram is null for pidType %d, effidx %d", pidType, effidx);
      return -1.0;
    }

    const int ibx = h->GetXaxis()->FindBin(mult);
    const int iby = h->GetYaxis()->FindBin(pt);
    const int ibz = h->GetZaxis()->FindBin(eta);
    float val = h->GetBinContent(ibx, iby, ibz);

    if (effidx == 1) { // This is a FAKE map
    // For fakes, "it's ok if it is zero"
    // But it's not ok if it's < 0 or >= 1. Set those to 0.
    if (val < 0.f || val >= 1.0) {
    return -1;
  }
  return val; // Return the valid fake value (which can be 0)
}

// If we are here, this is an EFFICIENCY map (effidx == 0)
// For efficiency, "it's ok if it is 1"
// "but not vice versa" (i.e., val <= 0 is the problem)

if (val > 0.f || val<1.1) {
// val > 0 is valid (including 1.0).
return val;
}

// --- PROBLEM CASE: val <= 0. We must interpolate. ---
// We will scan along the pT axis (iby) to find neighbors

// 1. Find bin "before" (lower pT)
int yBinBefore = iby - 1;
float valBefore = 0.f;
while (yBinBefore >= 1) {
valBefore = h->GetBinContent(ibx, yBinBefore, ibz);
if (valBefore > 0.f) break; // Found a valid bin
yBinBefore--;
}

// 2. Find bin "after" (higher pT)
int yBinAfter = iby + 1;
float valAfter = 0.f;
int nBinsY = h->GetNbinsY();
while (yBinAfter <= nBinsY) {
valAfter = h->GetBinContent(ibx, yBinAfter, ibz);
if (valAfter > 0.f) break; // Found a valid bin
yBinAfter++;
}

// 3. Interpolate/Extrapolate
if (valBefore > 0.f && valAfter > 0.f) {
// We found valid bins on both sides. Do linear interpolation.
float m = (valAfter - valBefore) / (yBinAfter - yBinBefore);
float interpolatedVal = valBefore + m * (iby - yBinBefore);
return interpolatedVal;
} else if (valBefore > 0.f) {
// Only found a bin before. Extrapolate using that value.
return valBefore;
} else if (valAfter > 0.f) {
// Only found a bin after. Extrapolate using that value.
return valAfter;
} else {
return -1;
}
}
*/
// Getter for (Cent, eta, phi) maps (Flattening)
float getFlatteningWeight(float cent, float eta, float phi, PID pidType) const {
  TH3F* h = hWeightMap3D[pidType];
  if (!h) return -1;
  const int ibx = h->GetXaxis()->FindBin(cent);
  const int iby = h->GetYaxis()->FindBin(eta);
  const int ibz = h->GetZaxis()->FindBin(phi);
  float val = h->GetBinContent(ibx, iby, ibz);
  return val;
}


template<int KIntM, int KIntK>
std::pair<float, float> calculateMeanAndC2FromSums(const double sumpmwk[KIntM][KIntK], const double sumwk[KIntK], float referenceMeanPt) const
{
  // --- Safety Checks ---
  if (sumwk[1] == 0.) {
    return {0.f, 0.f}; // No tracks, return 0
  }

  double tau1 = sumwk[2] / (sumwk[1] * sumwk[1]);
  double denom2 = 1. - tau1;

  if (std::abs(denom2) < kFloatEpsilon) { // Protect against divide-by-zero (e.g., N=1)
    // Return mean, but C2 is undefined
    double pmk11_safe = sumpmwk[1][1] / sumwk[1];
    return {static_cast<float>(pmk11_safe), 0.f};
  }

  // --- Calculations ---
  double pmk11 = sumpmwk[1][1] / sumwk[1];

  double pmk12 = 0.;
  if (sumwk[2] != 0.) {
    pmk12 = sumpmwk[1][2] / sumwk[2];
  }

  double pmk22 = 0.;
  if (sumwk[2] != 0.) {
    pmk22 = sumpmwk[2][2] / sumwk[2];
  }

  float calculatedMeanPt = pmk11; // Mean pT is <p^1 w^1> / <w^1>

  double p1kBar1 = pmk11 - referenceMeanPt;
  double p2kBar2 = pmk22 - 2. * pmk12 * referenceMeanPt + referenceMeanPt * referenceMeanPt;

  double p1kBar1_sq = p1kBar1 * p1kBar1;
  double numerator2 = p1kBar1_sq - (tau1 * p2kBar2);

  float C2 = numerator2 / denom2;
  return {calculatedMeanPt, C2};
}

ConfigurableAxis cfgAxisCent{"cfgAxisCent",{0.0, 1.0, 3.0, 5.0, 10, 20, 30, 40, 50, 60, 70, 80, 100},"centrality axis (percentile)"}; // FT0*/FV0A style
const AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
static constexpr int kNbinsNch = 5000;
static constexpr float kNchMin = 0.5f;
static constexpr float kNchMax = 5000.5f;
static constexpr int kNbinsNchCoarse = 500;
static constexpr float kNchCoarseMax = 5000.5f;
ConfigurableAxis nChAxis{"nChAxis", {kNbinsNch, kNchMin, kNchMax}, "PV-contributor track multiplicity axis"};
ConfigurableAxis nChAxis2{"nChAxis2", {kNbinsNchCoarse, kNchMin, kNchCoarseMax}, "PV-contributor track multiplicity axis"};

// run switches
Configurable<bool> cfgRunGetEff{"cfgRunGetEff", false, "Run MC pass to build efficiency/fake maps"};
Configurable<bool> cfgRunMCMean{"cfgRunMCMean", false, "Run MC mean(pT) & mean(Et)"};
Configurable<bool> cfgRunMCFluc{"cfgRunMCFluc", false, "Run MC fluctuations (C2, subevent)"};
Configurable<bool> cfgRunGetFlat{"cfgRunGetFlat", false, "Run Data Get Flattening Weights"};
Configurable<bool> cfgRunDataMean{"cfgRunDataMean", false, "Run DATA mean(pT) & mean(Et)"};
Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};
// --- DATA: collisions with centralities and selections ---
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
aod::track::pt > cfgPtMin &&
aod::track::pt < cfgPtMax &&
nabs(aod::track::dcaXY) < cfgDCAXY &&
nabs(aod::track::dcaZ) < cfgDCAZ;
using AodTracksSel = soa::Filtered<UnfilteredTracks>;
// --- MC: reconstructed tracks with truth labels ---
using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
using FilteredTCs = soa::Filtered<TCs>;
// --- MC: collisions with labels ---
using MyRun3MCCollisions = soa::Join<
aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As,
aod::CentNGlobals, aod::McCollisionLabels>;
// --- MC: truth-matched tracks ---
using MyMCTracks = soa::Join<
aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA,
aod::McTrackLabels,
aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
// --- Preslices (consistent policies) ---
PresliceUnsorted<aod::McParticles> partPerMcCollision = aod::mcparticle::mcCollisionId;
PresliceUnsorted<MyRun3MCCollisions> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
PresliceUnsorted<TCs> trackPerMcParticle = aod::mctracklabel::mcParticleId;
Preslice<MyMCTracks> perCollision = aod::track::collisionId;
Preslice<FilteredTCs> trackPerCollision = aod::track::collisionId;

void declareCommonQA()
{
  histos.add("hZvtx_after_sel", ";z_{vtx} (cm)", kTH1F, {{kNbinsZvtx, kZvtxMin, kZvtxMax}});
  histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis}});
  histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis2}, {nChAxis2}});
  histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis2}, {centAxis}});
  histos.add("hP", ";p (GeV/c)", kTH1F, {{kNbinsP, kPMin, kPMax}});
  histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{kNbinsPt, kPtMin, kPtMax}});
  histos.add("hEta", ";#eta", kTH1F, {{kNbinsEta, kEtaMin, kEtaMax}});
  histos.add("hPhi", ";#phi", kTH1F, {{kNbinsPhi, 0., TwoPI}});

  histos.add("hCentEtaPhi", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent1EtaPhi", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent7EtaPhi", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCentEtaPhiWtd", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent1EtaPhiWtd", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent7EtaPhiWtd", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});

  histos.add("hCentEtaPhiWtd_PID", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent1EtaPhiWtd_PID", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCent7EtaPhiWtd_PID", ";#eta;#phi", kTH2F, {{(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});

  histos.add("hCentEtaPhiTrue", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCentEtaPhiReco", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCentEtaPhiRecoMatched", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});

  histos.add("hCentEtaPhiTrue_PID", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCentEtaPhiReco_PID", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
  histos.add("hCentEtaPhiRecoMatched_PID", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});

}
void declareMCCommonHists()
{

  histos.add("ptResolution", ";p_{T}^{MC};p_{T}^{MC}-p_{T}^{reco}", kTH2F, {{kNbinsPtRes, 0., cfgPtMax}, {kNbinsPtRes, -kPtResMax, kPtResMax}});
  histos.add("ptTruthReco", ";p_{T}^{MC};p_{T}^{reco}", kTH2F, {{kNbinsPtRes, 0., cfgPtMax}, {kNbinsPtRes, 0., cfgPtMax}});
  histos.add("etaResolution", ";#eta^{MC};#eta^{MC}-#eta^{reco}", kTH2F, {{kNbinsEtaRes, -kEtaFineMax, kEtaFineMax}, {kNbinsPtRes, -kEtaResMax, kEtaResMax}});
  histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{kNbinsPtRes, -kEtaFineMax, kEtaFineMax}, {kNbinsPtRes, -kEtaFineMax, kEtaFineMax}});

  histos.add("TruthTrackVZ", ";Vz^{MC};Vz^{Reco}", kTH2F, {{kNbinsVz, kVzMin, kVzMax}, {kNbinsVz, kVzMin, kVzMax}});
  histos.add("vzResolution", ";Vz^{MC};Vz^{MC}-Vz^{Reco}", kTH2F, {{kNbinsVz, kVzMin, kVzMax}, {kNbinsVz, -kVzResMax, kVzResMax}});

  histos.add("h3_AllPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoMatchedToPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoUnMatchedToPrimary_Secondary",";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoUnMatchedToPrimary_Fake", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_AllReco", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});

  histos.add("h3_AllPrimary_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoMatchedToPrimary_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoUnMatchedToPrimary_Secondary_PID",";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_RecoUnMatchedToPrimary_Fake_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("h3_AllReco_PID", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {kNbinsPtRes, 0., cfgPtMax}, {kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});

  histos.add("h_AllPrimary", ";p_{T}", kTH1F, {{kNbinsP, 0., cfgPtMax}});
  histos.add("h_RecoMatchedToPrimary", ";p_{T}", kTH1F, {{kNbinsPt, 0., kPtMax}});
  histos.add("h_RecoUnMatchedToPrimary",";p_{T}", kTH1F, {{kNbinsPt, 0., kPtMax}});
  histos.add("h_AllReco", ";p_{T}", kTH1F, {{kNbinsPt, 0., kPtMax}});

  histos.add("hReco_ParticleWeight", ";cent;p_{T};#eta", kTH3F, {{centAxis},{kNbinsPtRes,0.,cfgPtMax},{kNbinsPtRes,-kEtaFineMax,kEtaFineMax}});
  histos.add("hTruth_ParticleWeight", ";cent;p_{T};#eta", kTH3F, {{centAxis},{kNbinsPtRes,0.,cfgPtMax},{kNbinsPtRes,-kEtaFineMax,kEtaFineMax}});

  histos.add("hDCAxy_Unmatched", ";DCA_{xy} (cm)", kTH1F, {{kNbinsDca, -kDcaMax, kDcaMax}});
  histos.add("hDCAz_Unmatched", ";DCA_{z} (cm)", kTH1F, {{kNbinsDca, -kDcaMax, kDcaMax}});
  histos.add("hDCAxy_NotPrimary", ";DCA_{xy} (cm)", kTH1F, {{kNbinsDca, -kDcaMax, kDcaMax}});
  histos.add("hDCAz_NotPrimary", ";DCA_{z} (cm)", kTH1F, {{kNbinsDca, -kDcaMax, kDcaMax}});
}

void declareMCMeanHists()
{
  // Eff/Fake and weights vs observables
  histos.add("Eff_cent", ";cent;#epsilon", kTProfile, {centAxis});
  histos.add("Fake_cent", ";cent;f_{fake}", kTProfile, {centAxis});
  histos.add("wgt_cent", ";cent;w", kTProfile, {centAxis});
  histos.add("Eff_Ntrk", ";N_{PV};#epsilon", kTProfile, {nChAxis2});
  histos.add("Fake_Ntrk", ";N_{PV};f_{fake}", kTProfile, {nChAxis2});
  histos.add("wgt_Ntrk", ";N_{PV};w", kTProfile, {nChAxis2});
  histos.add("Eff_pT", ";p_{T};#epsilon", kTProfile, {{kNbinsPtRes, 0., cfgPtMax}});
  histos.add("Fake_pT", ";p_{T};f_{fake}", kTProfile, {{kNbinsPtRes, 0., cfgPtMax}});
  histos.add("wgt_pT", ";p_{T};w", kTProfile, {{kNbinsPtRes, 0., kPtMax}});
  histos.add("Eff_eta", ";#eta;#epsilon", kTProfile, {{kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("Fake_eta", ";#eta;f_{fake}", kTProfile, {{kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  histos.add("wgt_eta", ";#eta;w", kTProfile, {{kNbinsEtaFine, -kEtaFineMax, kEtaFineMax}});
  // MC mean profiles (pT & Et) for various selections
  histos.add("MCGen/Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanpT_Mult", ";N_{PV};#LT p_{T}#GT",kTProfile, {nChAxis});
  histos.add<TProfile3D>("pmeanTruNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanRecoNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanRecoMatchedNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanRecoEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanRecoMatchedEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  // Et versions
  histos.add("MCGen/Prof_MeanEt_Cent",";cent;#LT E_{T}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanEt_Mult",";N_{PV};#LT E_{T}#GT",kTProfile, {nChAxis});
  histos.add<TProfile3D>("pmeanEtTruNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanEtRecoNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanEtRecoMatchedNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanEtRecoEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  histos.add<TProfile3D>("pmeanEtRecoMatchedEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -0.5, KNEta + 0.5}, {KNpT + 1, -0.5, KNpT + 0.5}});
  // Clone sets for Reco/Matched/EffCorr variants


  histos.addClone("MCGen/", "MCReco/");
  histos.addClone("MCGen/", "MCRecoMatched/");
  histos.addClone("MCGen/", "MCRecoEffCorr/");
  histos.addClone("MCGen/", "MCRecoMatchedEffCorr/");
}
void declareMCFlucHists()
{
  static constexpr int kNbinsNchFluc = 1000;
  // pT cumulants
  histos.add("MCGen/Prof_C2_Cent", ";cent;C_{2}", kTProfile, {centAxis});
  histos.add("MCGen/Prof_C2_Mult", ";N_{PV};C_{2}", kTProfile, {nChAxis});
  histos.add<TProfile3D>("MCGen/Prof_C2Sub_Mult_etabin_ptbin",
  ";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,
  {{kNbinsNchFluc, kNchMin, kNchCoarseMax},{KNEta+1,-kBinOffset,KNEta+kBinOffset},{KNpT+1,-kBinOffset,KNpT+kBinOffset}});
  histos.add<TProfile3D>("MCGen/Prof_ipt0_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("MCGen/Prof_ipt1_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("MCGen/Prof_ipt2_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt0_Cov_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt0_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt1_Cov_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt1_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt2_Cov_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt2_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});
  histos.add("MCGen/Prof_C2Et_Cent", ";cent;C_{2}^{E_{T}}", kTProfile, {centAxis});
  histos.add("MCGen/Prof_C2Et_Mult", ";N_{PV};C_{2}^{E_{T}}", kTProfile, {nChAxis});
  histos.add<TProfile3D>("MCGen/Prof_C2EtSub_Mult_etabin_ptbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,
  {{kNbinsNchFluc, kNchMin, kNchCoarseMax},{KNEta+1,-kBinOffset,KNEta+kBinOffset},{KNpT+1,-kBinOffset,KNpT+kBinOffset}});
  histos.add<TProfile3D>("MCGen/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("MCGen/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("MCGen/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt0_CovEt_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt0_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt1_CovEt_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt1_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});
  histos.add<TProfile2D>("MCGen/Prof_ipt2_CovEt_Cent_eta",";cent;#eta", kTProfile2D, {{centAxis},{(KNEta-1)/2,0,kEtaAxisMax}});
  histos.add("MCGen/Prof_ipt2_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2,0.,kEtaAxisMax}});

  histos.add("MCGen/Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanpT_Mult", ";N_{PV};#LT p_{T}#GT",kTProfile, {nChAxis});

  histos.add<TProfile3D>("pmeanTruNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanRecoNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanRecoMatchedNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanRecoEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanRecoMatchedEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});

  histos.add("MCGen/Prof_MeanEt_Cent",";cent;#LT E_{T}#GT", kTProfile, {centAxis});
  histos.add("MCGen/Prof_MeanEt_Mult",";N_{PV};#LT E_{T}#GT",kTProfile, {nChAxis});
  histos.add<TProfile3D>("pmeanEtTruNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanEtRecoNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanEtRecoMatchedNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanEtRecoEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add<TProfile3D>("pmeanEtRecoMatchedEffcorrNchEtabinPtbin",";N_{PV};#eta-bin; p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
}
void declareDataMeanHists()
{
  histos.add("Prof_cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis});
  histos.add("Prof_MeanpT_Cent", ";cent;#LT p_{T}#GT", kTProfile, {centAxis});
  histos.add<TProfile3D>("pmean_nch_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  // Et
  histos.add("Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
  histos.add<TProfile3D>("pmeanEt_nch_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
}

void declareDataGetFlatHists()
{
  // histos.add("hCentEtaPhi", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), -0.8, 0.8}, {30, 0, TMath::TwoPi()}});
  histos.add("hCentEtaPhi_PID", ";cent;#eta;#phi", kTH3F, {{centAxis}, {(KNEta - 1), kEtaAxisMin, kEtaAxisMax}, {kNbinsPhiFine, 0, TwoPI}});
}

void declareDataFlucHists()
{
  histos.add<TProfile3D>("pmean_nch_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});
  histos.add("Prof_MeanEt_Cent", ";cent;#LT E_{T}#GT", kTProfile, {centAxis});
  histos.add<TProfile3D>("pmeanEt_nch_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});

  histos.add("Prof_C2_Cent", ";cent;C_{2}", kTProfile, {centAxis});
  histos.add("Prof_MeanpT_Cent",";cent;#LT p_{T}#GT", kTProfile, {centAxis});
  histos.add("Prof_MeanpT_Mult",";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});
  histos.add<TProfile3D>("Prof_C2_Mult_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});

  histos.add<TProfile2D>("Prof_ipt0_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt0_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});
  histos.add<TProfile2D>("Prof_ipt1_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt1_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});
  histos.add<TProfile2D>("Prof_ipt2_Cov_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt2_Cov_Eta", ";#eta;cov", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});

  histos.add("Prof_C2Et_Cent", ";cent;C_{2}^{E_{T}}", kTProfile, {centAxis});
  histos.add("Prof_MeanEt_Mult", ";N_{PV};#LT E_{T}#GT", kTProfile, {nChAxis});
  histos.add<TProfile3D>("Prof_C2Et_Mult_etabin_ptbin",";N_{PV};#eta-bin;p_{T}-bin", kTProfile3D,{{nChAxis}, {KNEta + 1, -kBinOffset, KNEta + kBinOffset}, {KNpT + 1, -kBinOffset, KNpT + kBinOffset}});

  histos.add<TProfile2D>("Prof_ipt0_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt0_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});
  histos.add<TProfile2D>("Prof_ipt1_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt1_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});
  histos.add<TProfile2D>("Prof_ipt2_CovEt_Cent_eta", ";cent;#eta", kTProfile2D, {{centAxis}, {(KNEta-1)/2, 0, kEtaAxisMax}});
  histos.add("Prof_ipt2_CovEt_Eta", ";#eta;cov^{E_{T}}", kTProfile, {{(KNEta-1)/2, 0., kEtaAxisMax}});

  histos.add<TProfile3D>("Prof_ipt0_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("Prof_ipt1_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("Prof_ipt2_C2Sub2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("Prof_ipt0_C2SubEt2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("Prof_ipt1_C2SubEt2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
  histos.add<TProfile3D>("Prof_ipt2_C2SubEt2D_Mult_etaA_etaC",";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis},{KNEta-1,kEtaAxisMin,kEtaAxisMax},{KNEta-1,kEtaAxisMin,kEtaAxisMax}});
}

TH3F* buildWeightMapFromRaw(TH3F* hRaw, const char* mapName)
{
  if (!hRaw) {
    LOGF(error, "Raw eta-phi map for '%s' is null; no flattening will be applied.", mapName);
    return nullptr;
  }
  auto hWeightMap = (TH3F*)hRaw->Clone(mapName);
  hWeightMap->SetTitle(Form("Flattening Weight Map %s (w_{#phi} = <N_{#phi}> / N_{#phi})", mapName));
  hWeightMap->SetDirectory(nullptr);
  hWeightMap->Reset();
  auto axC = hRaw->GetXaxis();
  auto axE = hRaw->GetYaxis();
  auto axP = hRaw->GetZaxis();
  for (int ic = 1; ic <= axC->GetNbins(); ++ic) {
    for (int ie = 1; ie <= axE->GetNbins(); ++ie) {
      // average over phi at fixed (cent,eta)
      double sum = 0.0;
      int nphi = axP->GetNbins();
      for (int ip = 1; ip <= nphi; ++ip) sum += hRaw->GetBinContent(ic, ie, ip);
      const double avg = (nphi > 0 ? sum / nphi : 0.0);
      for (int ip = 1; ip <= nphi; ++ip) {
        const double raw = hRaw->GetBinContent(ic, ie, ip);
        const double w = (avg > 0.0 && raw > 0.0) ? (avg / raw) : 1.0;
        hWeightMap->SetBinContent(ic, ie, ip, w);
      }
    }
  }
  LOGF(info, "Flattening weight map '%s' built.", mapName);
  return hWeightMap;
}

void init(InitContext&)
{
  // CCDB
  ccdb->setURL("http://alice-ccdb.cern.ch");
  ccdb->setCaching(true);
  ccdb->setLocalObjectValidityChecking();
  int64_t nowMs = std::chrono::duration_cast<std::chrono::milliseconds>(
    std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(nowMs);
    // Always-on common QA
    declareCommonQA();
    // Conditionally declare sets
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
    if (cfgRunGetFlat) {
      declareDataGetFlatHists();
    }
    if (cfgRunDataMean) {
      declareDataMeanHists();
    }
    if (cfgRunDataFluc) {
      declareDataFlucHists();
    }

    //========================
    // Load correction maps
    //========================
    const bool needEffMaps = cfgRunMCMean || cfgRunMCFluc || cfgRunDataMean || cfgRunDataFluc;
    if (needEffMaps && !cfgRunGetEff) {

      // --- 1. Load Efficiency and Fake maps (always from MC file) ---
      if (auto* f = TFile::Open("Job1_EffMaps.root", "READ")) {
        if (!f->IsZombie()) {
          if (auto* dir = dynamic_cast<TDirectory*>(f->Get("rad-flow-decorr"))) {
            // Helper lambda to load Eff/Fake maps for a given PID type
            auto loadEffFakeForPID = [&](PID pidType) {
              std::string suffix = pidSuffix[pidType];
              std::string hEff_NumName = "h3_RecoMatchedToPrimary" + suffix;
              std::string hEff_DenName = "h3_AllPrimary" + suffix;
              std::string hFake_NumSecName = "h3_RecoUnMatchedToPrimary_Secondary" + suffix;
              std::string hFake_NumFakName = "h3_RecoUnMatchedToPrimary_Fake" + suffix;
              std::string hFake_DenName = "h3_AllReco" + suffix;

              // --- Efficiency ---
              if (auto* hNum = dynamic_cast<TH3F*>(dir->Get(hEff_NumName.c_str()))) {
                hEff[pidType] = (TH3F*)hNum->Clone(Form("hEff%s", suffix.c_str()));
                hEff[pidType]->SetDirectory(nullptr);
                if (auto* hDen = dynamic_cast<TH3F*>(dir->Get(hEff_DenName.c_str()))) {
                  hDen->SetDirectory(nullptr);
                  hEff[pidType]->Divide(hDen);
                  delete hDen;
                } else {
                  LOGF(error, "Missing denominator %s for efficiency.", hEff_DenName.c_str());
                }
              } else {
                LOGF(error, "Missing numerator %s for efficiency.", hEff_NumName.c_str());
              }
              // --- Fakes ---
              if (auto* hNumS = dynamic_cast<TH3F*>(dir->Get(hFake_NumSecName.c_str()))) {
                auto* hNumF = dynamic_cast<TH3F*>(dir->Get(hFake_NumFakName.c_str()));
                if (hNumS && hNumF) {
                  hFake[pidType] = (TH3F*)hNumS->Clone(Form("hFake%s", suffix.c_str()));
                  hFake[pidType]->Add(hNumF);
                  hFake[pidType]->SetDirectory(nullptr);
                  if (auto* hDenF = dynamic_cast<TH3F*>(dir->Get(hFake_DenName.c_str()))) {
                    hDenF->SetDirectory(nullptr);
                    hFake[pidType]->Divide(hDenF);
                    delete hDenF;
                  } else {
                    LOGF(error, "Missing denominator %s for fakes.", hFake_DenName.c_str());
                  }
                } else {
                  LOGF(error, "Missing fake numerator(s) for %s in EffMaps file.", suffix.c_str());
                }
              }
            };
            // Load Eff/Fake maps for both Inclusive and CombinedPID
            loadEffFakeForPID(kInclusive);
            loadEffFakeForPID(kCombinedPID);
          } else {
            LOGF(error, "Directory 'rad-flow-decorr' not found in Job1_EffMaps.root");
          }
        } else {
          LOGF(error, "Job1_EffMaps.root is a zombie.");
        }
        f->Close();
        delete f;
      } else {
        LOGF(error, "Cannot open Job1_EffMaps.root for Eff/Fake maps");
      }

      // --- 2. Load Flattening maps (source depends on run mode) ---
      const bool isDataRun = cfgRunDataMean || cfgRunDataFluc;
      if (isDataRun) {
        // Load from Data flattening file (produced by processGetFlat)
        LOGF(info, "Data Run: Loading flattening maps from Job1_DataFlat.root");
        if (auto* fFlat = TFile::Open("Job1_DataFlat.root", "READ")) {
          if (!fFlat->IsZombie()) {
            if (auto* dirFlat = dynamic_cast<TDirectory*>(fFlat->Get("rad-flow-decorr"))) {
              // Inclusive
              if (auto* hRawIncl = dynamic_cast<TH3F*>(dirFlat->Get("hCentEtaPhi"))) {
                hWeightMap3D[kInclusive] = buildWeightMapFromRaw(hRawIncl, "hWeightMap3D");
              } else {
                LOGF(error, "Data flattening source 'hCentEtaPhi' not found in Job1_DataFlat.root");
              }
              // PID
              if (auto* hRawPID = dynamic_cast<TH3F*>(dirFlat->Get("hCentEtaPhi_PID"))) {
                hWeightMap3D[kCombinedPID] = buildWeightMapFromRaw(hRawPID, "hWeightMap3D_PID");
              } else {
                LOGF(error, "Data flattening source 'hCentEtaPhi_PID' not found in Job1_DataFlat.root");
              }
            } else {
              LOGF(error, "Directory 'rad-flow-decorr' not found in Job1_DataFlat.root");
            }
          } else {
            LOGF(error, "Job1_DataFlat.root is a zombie.");
          }
          fFlat->Close();
          delete fFlat;
        } else {
          LOGF(error, "Cannot open Job1_DataFlat.root");
        }
      } else {
        // Load from MC efficiency file (as before, but ONLY flattening)
        LOGF(info, "MC Run: Loading flattening maps from Job1_EffMaps.root");
        if (auto* f = TFile::Open("Job1_EffMaps.root", "READ")) {
          if (!f->IsZombie()) {
            if (auto* dir = dynamic_cast<TDirectory*>(f->Get("rad-flow-decorr"))) {
              auto loadFlatForPID = [&](PID pidType) {
                std::string suffix = pidSuffix[pidType];
                // std::string hFlat_SrcName = "hCentEtaPhiRecoMatched" + suffix;
                std::string hFlat_SrcName = "hCentEtaPhiReco" + suffix; // As in original code
                if (auto* hRaw = dynamic_cast<TH3F*>(dir->Get(hFlat_SrcName.c_str()))) {
                  hWeightMap3D[pidType] = buildWeightMapFromRaw(hRaw, Form("hWeightMap3D%s", suffix.c_str()));
                } else {
                  LOGF(warning, "MC flattening source '%s' not found; proceeding without flattening for this PID.", hFlat_SrcName.c_str());
                }
              };
              loadFlatForPID(kInclusive);
              loadFlatForPID(kCombinedPID);
            } else {
              LOGF(error, "Directory 'rad-flow-decorr' not found in Job1_EffMaps.root (for MC flat)");
            }
          } else {
            LOGF(error, "Job1_EffMaps.root is a zombie (for MC flat).");
          }
          f->Close();
          delete f;
        } else {
          LOGF(error, "Cannot open Job1_EffMaps.root (for MC flat)");
        }
      }
    }

    //========================
    // Load Step-2 mean profiles (pT and Et)
    //========================
    auto loadTProfile3D = [&](const char* filename, const char* histname, TProfile3D*& target) {
      if (auto* f = TFile::Open(filename, "READ")) {
        if (!f->IsZombie()) {
          if (auto* dir = dynamic_cast<TDirectory*>(f->Get("rad-flow-decorr"))) {
            if (auto* tp = dynamic_cast<TProfile3D*>(dir->Get(histname))) {
              target = (TProfile3D*)tp->Clone();
              target->SetDirectory(nullptr);
            } else {
              LOGF(error, "Histogram %s missing in %s", histname, filename);
            }
          } else {
            LOGF(error, "Directory 'rad-flow-decorr' missing in %s", filename);
          }
        } else {
          LOGF(error, "%s is a zombie", filename);
        }
        f->Close();
        delete f;
      } else {
        LOGF(error, "Cannot open %s", filename);
      }
    };
    // === THIS BLOCK IS UPDATED ===
    if (cfgRunMCFluc) {
      // pT profiles
      loadTProfile3D("Job2_MCMean.root", "pmeanTruNchEtabinPtbin", pmeanTruNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanRecoNchEtabinPtbin", pmeanRecoNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanRecoMatchedNchEtabinPtbin", pmeanRecoMatchedNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanRecoEffcorrNchEtabinPtbin", pmeanRecoEffcorrNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanRecoMatchedEffcorrNchEtabinPtbin", pmeanRecoMatchedEffcorrNchEtabinPtbinStep2);
      // Et profiles
      loadTProfile3D("Job2_MCMean.root", "pmeanEtTruNchEtabinPtbin", pmeanEtTruNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanEtRecoNchEtabinPtbin", pmeanEtRecoNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanEtRecoMatchedNchEtabinPtbin", pmeanEtRecoMatchedNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanEtRecoEffcorrNchEtabinPtbin", pmeanEtRecoEffcorrNchEtabinPtbinStep2);
      loadTProfile3D("Job2_MCMean.root", "pmeanEtRecoMatchedEffcorrNchEtabinPtbin", pmeanEtRecoMatchedEffcorrNchEtabinPtbinStep2);
    }
    // =============================
    if (cfgRunDataFluc) {
      loadTProfile3D("Job1_DataMean.root", "pmean_nch_etabin_ptbin", pmeanNchEtabinPtbinStep2);
      loadTProfile3D("Job1_DataMean.root", "pmeanEt_nch_etabin_ptbin", pmeanEtNchEtabinPtbinStep2);
    }}

  void processGetEffHists(aod::McCollisions const& mcColl, soa::SmallGroups<MyRun3MCCollisions> const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      if (colSlice.size() != 1) continue;

      for (const auto& col : colSlice) {
        if (!col.has_mcCollision()) continue;
        if (!isEventSelected(col)) continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1) continue;


        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (partSlice.size() < 1) continue;
        if (col.globalIndex() >= trackSlice.size()) {
          LOGF(warning, "Skipping invalid globalIndex=%d for tracks (tracks.size=%d)",col.globalIndex(), tracks.size());
          continue;
        }

        float cent = getCentrality(col);
        if (cent < 0 || cent > 80) continue;

        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle)) continue;
          if (!particle.isPhysicalPrimary()) continue;

          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == kPiPlus);
          const bool isKaon = (absPdgId == kKPlus);
          const bool isProton = (absPdgId == kProton);
          const bool isPid = (isPion || isKaon || isProton);

          histos.fill(HIST("hTruth_ParticleWeight"), cent, particle.pt(), particle.eta(), particle.weight());
          histos.fill(HIST("hCentEtaPhiTrue"), cent, particle.eta(), particle.phi());
          histos.fill(HIST("h3_AllPrimary"), col.multNTracksPV(), particle.pt(), particle.eta());

          if (cent < kCentTestMin) histos.fill(HIST("hCent1EtaPhi"), particle.eta(), particle.phi());
          if (cent > kCentTestMaxLo && cent < kCentTestMaxHi) histos.fill(HIST("hCent7EtaPhi"), particle.eta(), particle.phi());

          if (isPid) {
            histos.fill(HIST("hCentEtaPhiTrue_PID"), cent, particle.eta(), particle.phi());
            histos.fill(HIST("h3_AllPrimary_PID"), col.multNTracksPV(), particle.pt(), particle.eta());
          }

        }
        histos.fill(HIST("TruthTrackVZ"), mcCollision.posZ(), col.posZ());
        histos.fill(HIST("vzResolution"), mcCollision.posZ(), mcCollision.posZ() - col.posZ());


        // Reconstructed
        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track)) continue;

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);
          const bool isPid = (isPion || isKaon || isProton);

          histos.fill(HIST("h3_AllReco"), col.multNTracksPV(), track.pt(), track.eta());
          histos.fill(HIST("hCentEtaPhiReco"), cent, track.eta(), track.phi());

          if (isPid) {
            histos.fill(HIST("h3_AllReco_PID"), col.multNTracksPV(), track.pt(), track.eta());
            histos.fill(HIST("hCentEtaPhiReco_PID"), cent, track.eta(), track.phi());
          }

          if (track.has_mcParticle()) {
            auto mcPart2 = track.mcParticle();
            if (mcPart2.isPhysicalPrimary()) {
              const int absPdgId = std::abs(mcPart2.pdgCode());
              const bool isPionTrue = (absPdgId == kPiPlus);
              const bool isKaonTrue = (absPdgId == kKPlus);
              const bool isProtonTrue = (absPdgId == kProton);
              const bool isPidTrue = (isPionTrue || isKaonTrue || isProtonTrue);
              // Fill inclusive hists
              histos.fill(HIST("hReco_ParticleWeight"), cent, mcPart2.pt(), mcPart2.eta(), mcPart2.weight());
              histos.fill(HIST("ptResolution"), mcPart2.pt(), mcPart2.pt() - track.pt());
              histos.fill(HIST("ptTruthReco"), mcPart2.pt(), track.pt());
              histos.fill(HIST("etaResolution"), mcPart2.eta(), mcPart2.eta() - track.eta());
              histos.fill(HIST("etaTruthReco"), mcPart2.eta(), track.eta());
              histos.fill(HIST("h3_RecoMatchedToPrimary"), col.multNTracksPV(), mcPart2.pt(), mcPart2.eta());
              histos.fill(HIST("hCentEtaPhiRecoMatched"), cent, mcPart2.eta(), mcPart2.phi());

              if (isPid && isPidTrue) {
                histos.fill(HIST("h3_RecoMatchedToPrimary_PID"), col.multNTracksPV(), mcPart2.pt(), mcPart2.eta());
                histos.fill(HIST("hCentEtaPhiRecoMatched_PID"), cent, mcPart2.eta(), mcPart2.phi());
              }

            }
            else {
              // Matched to secondary
              histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), col.multNTracksPV(), track.pt(), track.eta());
              histos.fill(HIST("hDCAxy_Unmatched"), track.dcaXY());
              histos.fill(HIST("hDCAz_Unmatched"), track.dcaZ());
              if (isPid) {
                histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary_PID"), col.multNTracksPV(), track.pt(), track.eta());
              }
            }
          }
          else {
            // Fake track
            histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), col.multNTracksPV(), track.pt(), track.eta());
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
  PROCESS_SWITCH(RadFlowDecorr, processGetEffHists, "process MC to calculate Eff and Fakes", cfgRunGetEff);


  void processMCMean(aod::McCollisions const& mcColl, MyRun3MCCollisions const& collisions, TCs const& tracks, FilteredTCs const& /*filteredTracks*/, aod::McParticles const& mcParticles)
  {
    float sumWiTruth[KNEta][KNpT], sumWiptiTruth[KNEta][KNpT];
    float sumWiReco[KNEta][KNpT], sumWiptiReco[KNEta][KNpT];
    float sumWiRecoMatched[KNEta][KNpT], sumWiptiRecoMatched[KNEta][KNpT];
    float sumWiRecoEffCorr[KNEta][KNpT], sumWiptiRecoEffCorr[KNEta][KNpT];
    float sumWiRecoMatchedEffCorr[KNEta][KNpT], sumWiptiRecoMatchedEffCorr[KNEta][KNpT];
    float sumWiTruthEt[KNEta][KNpT], sumWiptiTruthEt[KNEta][KNpT];
    float sumWiRecoEt[KNEta][KNpT], sumWiptiRecoEt[KNEta][KNpT];
    float sumWiRecoMatchedEt[KNEta][KNpT], sumWiptiRecoMatchedEt[KNEta][KNpT];
    float sumWiRecoEffCorrEt[KNEta][KNpT], sumWiptiRecoEffCorrEt[KNEta][KNpT];
    float sumWiRecoMatchedEffCorrEt[KNEta][KNpT], sumWiptiRecoMatchedEffCorrEt[KNEta][KNpT];


    for (const auto& mcCollision : mcColl) {
      auto colSlice = collisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());

      if (colSlice.size() != 1) continue;
      for (const auto& col : colSlice) {
        if (!col.has_mcCollision()) continue;
        if (!isEventSelected(col)) continue;

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1) continue;


        auto partSlice = mcParticles.sliceBy(partPerMcCollision, mcCollision.globalIndex());
        if (partSlice.size() < 1) continue;
        if (col.globalIndex() >= trackSlice.size()) {
          LOGF(warning, "Skipping invalid globalIndex=%d for tracks (tracks.size=%d)",col.globalIndex(), tracks.size());
          continue;
        }

        float cent = getCentrality(col);
        if (cent < 0 || cent > 80) continue;


        if (col.globalIndex() >= trackSlice.size()) {
          LOGF(warning, "Skipping invalid globalIndex=%d for tracks (tracks.size=%d)",col.globalIndex(), tracks.size());
          continue;
        }

        LOGF(info, "Event Check: cent = %.1f, nTracks = %d", cent, (int)trackSlice.size());
        memset(sumWiTruth, 0, sizeof(sumWiTruth));
        memset(sumWiptiTruth, 0, sizeof(sumWiptiTruth));
        memset(sumWiReco, 0, sizeof(sumWiReco));
        memset(sumWiptiReco, 0, sizeof(sumWiptiReco));
        memset(sumWiRecoMatched, 0, sizeof(sumWiRecoMatched));
        memset(sumWiptiRecoMatched, 0, sizeof(sumWiptiRecoMatched));
        memset(sumWiRecoEffCorr, 0, sizeof(sumWiRecoEffCorr));
        memset(sumWiptiRecoEffCorr, 0, sizeof(sumWiptiRecoEffCorr));
        memset(sumWiRecoMatchedEffCorr, 0, sizeof(sumWiRecoMatchedEffCorr));
        memset(sumWiptiRecoMatchedEffCorr, 0, sizeof(sumWiptiRecoMatchedEffCorr));
        memset(sumWiTruthEt, 0, sizeof(sumWiTruthEt));
        memset(sumWiptiTruthEt, 0, sizeof(sumWiptiTruthEt));
        memset(sumWiRecoEt, 0, sizeof(sumWiRecoEt));
        memset(sumWiptiRecoEt, 0, sizeof(sumWiptiRecoEt));
        memset(sumWiRecoMatchedEt, 0, sizeof(sumWiRecoMatchedEt));
        memset(sumWiptiRecoMatchedEt, 0, sizeof(sumWiptiRecoMatchedEt));
        memset(sumWiRecoEffCorrEt, 0, sizeof(sumWiRecoEffCorrEt));
        memset(sumWiptiRecoEffCorrEt, 0, sizeof(sumWiptiRecoEffCorrEt));
        memset(sumWiRecoMatchedEffCorrEt, 0, sizeof(sumWiRecoMatchedEffCorrEt));
        memset(sumWiptiRecoMatchedEffCorrEt, 0, sizeof(sumWiptiRecoMatchedEffCorrEt));

        // Truth
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle)) continue;
          if (!particle.isPhysicalPrimary()) continue;

          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == kPiPlus);
          const bool isKaon = (absPdgId == kKPlus);
          const bool isProton = (absPdgId == kProton);

          float pt = particle.pt();
          float eta = particle.eta();
          float p = particle.p();
          const double wi = 1.0;

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
              sumWiTruth[ieta][ipt] += 1.0;
              sumWiptiTruth[ieta][ipt] += 1.0 * pt;
              if (isPion || isKaon || isProton) {
                float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus : o2::constants::physics::MassProton;
                float E = std::sqrt(p * p + m * m);
                float Et = E * (pt / p); // E_T = E * sin(theta) = E * (pT / p)
                sumWiTruthEt[ieta][ipt] += 1.0;
                sumWiptiTruthEt[ieta][ipt] += Et;
              }
            }
          }
        }

        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track)) continue;

          float pt = track.pt();
          float eta = track.eta();
          float p = track.p();
          float phi = track.phi();


          histos.fill(HIST("hCentEtaPhi"), cent, eta, phi);
          if (cent < kCentTestMin) histos.fill(HIST("hCent1EtaPhi"), eta, phi);
          if (cent > kCentTestMaxLo && cent < kCentTestMaxHi) histos.fill(HIST("hCent7EtaPhi"), eta, phi);

          float effIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 0);
          float fakeIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 1);
          float flatWeightIncl = getFlatteningWeight(cent, eta, phi, kInclusive);
          float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;

          histos.fill(HIST("hCentEtaPhiWtd"), cent, eta, track.phi(), flatWeightIncl);
          if (cent < kCentTestMin) histos.fill(HIST("hCent1EtaPhiWtd"), track.eta(), track.phi(), flatWeightIncl);
          if (cent > kCentTestMaxLo && cent < kCentTestMaxHi) histos.fill(HIST("hCent7EtaPhiWtd"), track.eta(), track.phi(), flatWeightIncl);

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
              sumWiReco[ieta][ipt] += 1.0;
              sumWiptiReco[ieta][ipt] += pt;
            }
          }

          if (effIncl <= 0 || !isfinite(wIncl) || !isfinite(fakeIncl) || !isfinite(flatWeightIncl))  continue;

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
              sumWiRecoEffCorr[ieta][ipt] += wIncl;
              sumWiptiRecoEffCorr[ieta][ipt] += wIncl * pt;
            }
          }

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);
          if (isPion || isKaon || isProton) {
            float effPid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 0);
            float fakePid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 1);
            float flatWeightPid = getFlatteningWeight(cent, eta, phi, kCombinedPID);
            float wPid =  flatWeightPid * (1.0 - fakePid) / effPid;

            histos.fill(HIST("hCentEtaPhiWtd_PID"), cent, eta, track.phi(), flatWeightPid);
            if (cent < kCentTestMin) histos.fill(HIST("hCent1EtaPhiWtd_PID"), track.eta(), track.phi(), flatWeightPid);
            if (cent > kCentTestMaxLo && cent < kCentTestMaxHi) histos.fill(HIST("hCent7EtaPhiWtd_PID"), track.eta(), track.phi(), flatWeightPid);

            float m = isPion ? o2::constants::physics::MassPiPlus :isKaon ? o2::constants::physics::MassKPlus :o2::constants::physics::MassProton;
            float E = std::sqrt(p * p + m * m);
            float Et = E * (pt / p); // E_T = E * sin(theta)
            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
                sumWiRecoEt[ieta][ipt] += 1.0;
                sumWiptiRecoEt[ieta][ipt] += Et;
              }
            }

            if (effPid <= kFloatEpsilon || !isfinite(wPid) || !isfinite(fakePid) || !isfinite(flatWeightPid)) continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
                sumWiRecoEffCorrEt[ieta][ipt] += wPid;
                sumWiptiRecoEffCorrEt[ieta][ipt] += wPid * Et;
              }
            }
          }

          if (isfinite(wIncl)){
            if(cent < kCentTestMin){
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

        if (isfinite(sumWiTruth[0][0])) {
          float meanpT_truth = sumWiptiTruth[0][0] / sumWiTruth[0][0];
          if(!isfinite(meanpT_truth)) LOGF(info, "meanpT_truth = %.3f, num = %.3f, den =%.3f", meanpT_truth, sumWiptiTruth[0][0], sumWiTruth[0][0]);
          if(!isfinite(meanpT_truth)) continue;
          histos.fill(HIST("MCGen/Prof_cent_Nchrec"), cent, sumWiTruth[0][0]);
          histos.fill(HIST("MCGen/Prof_MeanpT_Cent"), cent, meanpT_truth);
          histos.fill(HIST("MCGen/Prof_MeanpT_Mult"), col.multNTracksPV(), meanpT_truth);
        }
        if (isfinite(sumWiReco[0][0])) {
          float meanpT_reco = sumWiptiReco[0][0] / sumWiReco[0][0];
          if(!isfinite(meanpT_reco)) LOGF(info, "meanpT_reco = %.3f, num = %.3f, den =%.3f", meanpT_reco, sumWiptiReco[0][0], sumWiReco[0][0]);
          if(!isfinite(meanpT_reco)) continue;
          histos.fill(HIST("MCReco/Prof_cent_Nchrec"), cent, sumWiReco[0][0]);
          histos.fill(HIST("MCReco/Prof_MeanpT_Cent"), cent, meanpT_reco);
          histos.fill(HIST("MCReco/Prof_MeanpT_Mult"), col.multNTracksPV(), meanpT_reco);
        }
        if (isfinite(sumWiRecoEffCorr[0][0])) {
          float meanpT_effcorr = sumWiptiRecoEffCorr[0][0] / sumWiRecoEffCorr[0][0];
          if(!isfinite(meanpT_effcorr)) LOGF(info, "meanpT_recoEffcorr = %.3f, num = %.3f, den =%.3f", meanpT_effcorr, sumWiptiRecoEffCorr[0][0], sumWiRecoEffCorr[0][0]);
          if(!isfinite(meanpT_effcorr)) continue;
          histos.fill(HIST("MCRecoEffCorr/Prof_cent_Nchrec"), cent, sumWiRecoEffCorr[0][0]);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent"), cent, meanpT_effcorr);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult"), col.multNTracksPV(), meanpT_effcorr);
        }

        if (isfinite(sumWiTruthEt[0][0])) {
          float meanEt = sumWiptiTruthEt[0][0] / sumWiTruthEt[0][0];
          if(!isfinite(meanEt)) LOGF(info, "meanEtTruthEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiTruthEt[0][0], sumWiTruthEt[0][0]);
          if(!isfinite(meanEt)) continue;
          histos.fill(HIST("MCGen/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCGen/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }
        // "MCReco"
        if (isfinite(sumWiRecoEt[0][0])) {
          float meanEt = sumWiptiRecoEt[0][0] / sumWiRecoEt[0][0];
          if(!isfinite(meanEt)) LOGF(info, "meanEtRecoEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiRecoEt[0][0], sumWiRecoEt[0][0]);
          if(!isfinite(meanEt)) continue;
          histos.fill(HIST("MCReco/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCReco/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }
        // "MCRecoEffCorr"
        if (isfinite(sumWiRecoEffCorrEt[0][0])) {
          float meanEt = sumWiptiRecoEffCorrEt[0][0] / sumWiRecoEffCorrEt[0][0];
          if(!isfinite(meanEt)) LOGF(info, "meanEtRecoEffcorrEt = %.3f, num = %.3f, den =%.3f", meanEt, sumWiptiRecoEffCorrEt[0][0], sumWiRecoEffCorrEt[0][0]);
          if(!isfinite(meanEt)) continue;
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Cent"), cent, meanEt);
          histos.fill(HIST("MCRecoEffCorr/Prof_MeanEt_Mult"), col.multNTracksPV(), meanEt);
        }

        // --- Fill 3D pT Profiles ---
        for (int ieta = 0; ieta < KNEta; ++ieta){
          for (int ipt = 0; ipt < KNpT; ++ipt){
            if (isfinite(sumWiTruth[ieta][ipt])) histos.fill(HIST("pmeanTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiTruth[ieta][ipt] / sumWiTruth[ieta][ipt]);
            if (isfinite(sumWiReco[ieta][ipt])) histos.fill(HIST("pmeanRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,sumWiptiReco[ieta][ipt] / sumWiReco[ieta][ipt]);
            if (isfinite(sumWiRecoEffCorr[ieta][ipt])) histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,sumWiptiRecoEffCorr[ieta][ipt] / sumWiRecoEffCorr[ieta][ipt]);
            // --- Fill 3D Et Profiles ---
            if (isfinite(sumWiTruthEt[ieta][ipt])) histos.fill(HIST("pmeanEtTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiTruthEt[ieta][ipt] / sumWiTruthEt[ieta][ipt]);
            if (isfinite(sumWiRecoEt[ieta][ipt])) histos.fill(HIST("pmeanEtRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiRecoEt[ieta][ipt] / sumWiRecoEt[ieta][ipt]);
            if (isfinite(sumWiRecoEffCorrEt[ieta][ipt])) histos.fill(HIST("pmeanEtRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt, sumWiptiRecoEffCorrEt[ieta][ipt] / sumWiRecoEffCorrEt[ieta][ipt]);
          }
        }
      } // end col loop
    }

    LOGF(info, "FINISHED RUNNING processMCMean (pT + Et)");
  }
  PROCESS_SWITCH(RadFlowDecorr, processMCMean, "process MC to calculate mean pt/Et and Eff Hists", cfgRunMCMean);

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
      if (colSlice.size() != 1) continue;
      histos.fill(HIST("MCGen/hVtxZ"), mcCollision.posZ()); // This histogram was not declared, commenting out.
      for (const auto& col : colSlice) {

        auto trackSlice = tracks.sliceBy(trackPerCollision, col.globalIndex());
        if (trackSlice.size() < 1) continue;

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

        memset(meanTru, 0, sizeof(meanTru)); memset(c2Tru, 0, sizeof(c2Tru));
        memset(meanReco, 0, sizeof(meanReco)); memset(c2Reco, 0, sizeof(c2Reco));
        memset(meanRecoEffCor, 0, sizeof(meanRecoEffCor)); memset(c2RecoEffCor, 0, sizeof(c2RecoEffCor));

        memset(meanTruEt, 0, sizeof(meanTruEt)); memset(c2TruEt, 0, sizeof(c2TruEt));
        memset(meanRecoEt, 0, sizeof(meanRecoEt)); memset(c2RecoEt, 0, sizeof(c2RecoEt));
        memset(meanRecoEffCorEt, 0, sizeof(meanRecoEffCorEt)); memset(c2RecoEffCorEt, 0, sizeof(c2RecoEffCorEt));

        if (!col.has_mcCollision() || !isEventSelected(col)) continue;
        float cent = getCentrality(col);
        if (cent < 0 || cent > 80) continue;

        // truth
        for (const auto& particle : partSlice) {
          if (!isParticleSelected(particle)) continue;
          if (!particle.isPhysicalPrimary()) continue;
          float pt = particle.pt();
          float eta = particle.eta();
          float p = particle.p();
          const double wi = 1.0;

          // pT (Inclusive)
          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
              for (int k = 0; k < KIntK; ++k) {
                for (int m = 0; m < KIntM; ++m) {
                  sumPmwkTru[ieta][ipt][m][k] += std::pow(wi, k) * std::pow(pt, m);
                }
                sumWkTru[ieta][ipt][k] += std::pow(wi, k);
              }
            }
          }
          const int absPdgId = std::abs(particle.pdgCode());
          const bool isPion = (absPdgId == kPiPlus);
          const bool isKaon = (absPdgId == kKPlus);
          const bool isProton = (absPdgId == kProton);
          if (isPion || isKaon || isProton){

            float m = isPion ? o2::constants::physics::MassPiPlus : isKaon ? o2::constants::physics::MassKPlus : o2::constants::physics::MassProton;
            float E = std::sqrt(p * p + m * m);
            float Et = E * (pt / p);
            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
                for (int k = 0; k < KIntK; ++k) {
                  for (int m = 0; m < KIntM; ++m) {
                    sumPmwkTruEt[ieta][ipt][m][k] += std::pow(wi, k) * std::pow(Et, m);
                  }
                  sumWkTruEt[ieta][ipt][k] += std::pow(wi, k);
                }
              }
            }
          }

        } // end truth loop

        for (const auto& track : trackSlice) {
          if (!isTrackSelected(track)) continue;
          float pt = track.pt();
          float eta = track.eta();
          float p = track.p();
          float phi = track.phi();

          float effIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 0);
          float fakeIncl = getEfficiency(col.multNTracksPV(), pt, eta, kInclusive, 1);
          float flatWeightIncl = getFlatteningWeight(cent, eta, phi, kInclusive);
          float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
          if (!std::isfinite(wIncl) || wIncl <= 0.f) continue;
          if (effIncl <= 0 || !isfinite(effIncl) || !isfinite(fakeIncl) || !isfinite(flatWeightIncl))  continue;

          for (int ieta = 0; ieta < KNEta; ++ieta) {
            if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
            for (int ipt = 0; ipt < KNpT; ++ipt) {
              if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
              for (int k = 0; k < KIntK; ++k) {
                for (int m = 0; m < KIntM; ++m) {
                  sumPmwkReco[ieta][ipt][m][k] += std::pow(1.0f, k) * std::pow(pt, m);
                  sumPmwkRecoEffCor[ieta][ipt][m][k] += std::pow(wIncl, k) * std::pow(pt, m);
                }
                sumWkReco[ieta][ipt][k] += std::pow(1.0f, k);
                sumWkRecoEffCor[ieta][ipt][k] += std::pow(wIncl, k);
              }
            }
          }

          const bool isPion = selectionPion(track);
          const bool isKaon = selectionKaon(track);
          const bool isProton = selectionProton(track);

          if (isPion || isKaon || isProton){
            float m = isPion ? o2::constants::physics::MassPiPlus :isKaon ? o2::constants::physics::MassKPlus :o2::constants::physics::MassProton;
            float E = std::sqrt(p * p + m * m);
            float Et = E * (pt / p); // E_T = E * sin(theta)
            float effPid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 0);
            float fakePid = getEfficiency(col.multNTracksPV(), pt, eta, kCombinedPID, 1);
            float flatWeightPid = getFlatteningWeight(cent, eta, phi, kCombinedPID);
            float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
            if(effPid >= 1.f || fakePid >= 1.f || !isfinite(effPid) || effPid <= kFloatEpsilon || !isfinite(fakePid) || !isfinite(flatWeightPid)) continue;

            for (int ieta = 0; ieta < KNEta; ++ieta) {
              if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
              for (int ipt = 0; ipt < KNpT; ++ipt) {
                if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
                for (int k = 0; k < KIntK; ++k) {
                  for (int m = 0; m < KIntM; ++m) {
                    sumPmwkRecoEt[ieta][ipt][m][k] += std::pow(1.0f, k) * std::pow(Et, m);
                    sumPmwkRecoEffCorEt[ieta][ipt][m][k] += std::pow(wPid, k) * std::pow(Et, m);
                  }
                  sumWkRecoEt[ieta][ipt][k] += std::pow(1.0f, k);
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

              if (isfinite(mmptTru)) std::tie(meanTru[ieta][ipt], c2Tru[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTru[ieta][ipt], sumWkTru[ieta][ipt], mmptTru);
              if (isfinite(mmptReco)) std::tie(meanReco[ieta][ipt], c2Reco[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkReco[ieta][ipt], sumWkReco[ieta][ipt], mmptReco);
              if (isfinite(mmptRecoEffCor)) std::tie(meanRecoEffCor[ieta][ipt], c2RecoEffCor[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCor[ieta][ipt], sumWkRecoEffCor[ieta][ipt], mmptRecoEffCor);

              if (isfinite(mmetTru)) std::tie(meanTruEt[ieta][ipt], c2TruEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTruEt[ieta][ipt], sumWkTruEt[ieta][ipt], mmetTru);
              if (isfinite(mmetReco)) std::tie(meanRecoEt[ieta][ipt], c2RecoEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEt[ieta][ipt], sumWkRecoEt[ieta][ipt], mmetReco);
              if (isfinite(mmetRecoEffCor)) std::tie(meanRecoEffCorEt[ieta][ipt], c2RecoEffCorEt[ieta][ipt]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCorEt[ieta][ipt], sumWkRecoEffCorEt[ieta][ipt], mmetRecoEffCor);
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

        for (int ieta = 0; ieta < KNEta; ++ieta){
          for (int ipt = 0; ipt < KNpT; ++ipt){
            if (std::isfinite(sumWkTru[ieta][ipt][1])) histos.fill(HIST("pmeanTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanTru[ieta][ipt]);
            if (std::isfinite(sumWkReco[ieta][ipt][1])) histos.fill(HIST("pmeanRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanReco[ieta][ipt]);
            if (std::isfinite(sumWkRecoEffCor[ieta][ipt][1])) histos.fill(HIST("pmeanRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanRecoEffCor[ieta][ipt]);
            if (std::isfinite(sumWkTruEt[ieta][ipt][1])) histos.fill(HIST("pmeanEtTruNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanTruEt[ieta][ipt]);
            if (std::isfinite(sumWkRecoEt[ieta][ipt][1])) histos.fill(HIST("pmeanEtRecoNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanRecoEt[ieta][ipt]);
            if (std::isfinite(sumWkRecoEffCorEt[ieta][ipt][1])) histos.fill(HIST("pmeanEtRecoEffcorrNchEtabinPtbin"), col.multNTracksPV(), ieta, ipt,meanRecoEffCorEt[ieta][ipt]);
          }
        }

        float p1kBarTru[KNEta][KNpT]{}, p1kBarReco[KNEta][KNpT]{},p1kBarRecoEffCor[KNEta][KNpT]{};
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

            // pT deviations
            if (mmptTru != 0.0) p1kBarTru[ieta][ipt] = meanTru[ieta][ipt] - mmptTru;
            if (mmptReco != 0.0) p1kBarReco[ieta][ipt] = meanReco[ieta][ipt] - mmptReco;
            if (mmptRecoEffCor != 0.0) p1kBarRecoEffCor[ieta][ipt] = meanRecoEffCor[ieta][ipt] - mmptRecoEffCor;
            // Et deviations
            if (mmetTru != 0.0) p1kBarTruEt[ieta][ipt] = meanTruEt[ieta][ipt] - mmetTru;
            if (mmetReco != 0.0) p1kBarRecoEt[ieta][ipt] = meanRecoEt[ieta][ipt] - mmetReco;
            if (mmetRecoEffCor != 0.0)p1kBarRecoEffCorEt[ieta][ipt] = meanRecoEffCorEt[ieta][ipt] - mmetRecoEffCor;
          }
        }


        // 1D Covariance (vs eta)
        for (int ietaA = 1; ietaA <=(KNEta-1)/2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);
          {
            const int ipt = 0;
            float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCGen/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt0_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCGen/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCGen/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
            }
          }

          // ipt = 1
          {
            const int ipt = 1;
            float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCGen/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt0_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCGen/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCGen/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
            }
          }
          // ipt = 2
          {
            const int ipt = 2;
            float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCGen/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt2_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCGen/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCGen/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCGen/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
            }
          }
        }

        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = 0.5f * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);
            // ipt = 0
            {
              const int ipt = 0;
              float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
              if(std::isfinite(c2Sub)) histos.fill(HIST("MCGen/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
              if(std::isfinite(c2SubEt)) histos.fill(HIST("MCGen/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }
            // ipt = 1
            {
              const int ipt = 1;
              float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
              if(std::isfinite(c2Sub)) histos.fill(HIST("MCGen/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
              if(std::isfinite(c2SubEt)) histos.fill(HIST("MCGen/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }
            // ipt = 2
            {
              const int ipt = 2;
              float c2Sub = p1kBarTru[ietaA][ipt] * p1kBarTru[ietaC][ipt];
              if(std::isfinite(c2Sub)) histos.fill(HIST("MCGen/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarTruEt[ietaA][ipt] * p1kBarTruEt[ietaC][ipt];
              if(std::isfinite(c2SubEt)) histos.fill(HIST("MCGen/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }
          }
        }


        // --- MCReco Covariance (1D vs eta) ---
        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]); // 0.5 for averaging

          // ipt = 0
          {
            const int ipt = 0;
            float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCReco/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt0_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCReco/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCReco/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
            }
          }

          // ipt = 1
          {
            const int ipt = 1;
            float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCReco/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt1_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCReco/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCReco/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
            }
          }

          // ipt = 2
          {
            const int ipt = 2;
            float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCReco/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt2_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCReco/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCReco/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCReco/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
            }
          }
        }

        // --- MCReco Covariance (2D etaA vs etaC) ---
        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = 0.5f * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);

            // ipt = 0
            {
              const int ipt = 0;
              float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCReco/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCReco/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }

            // ipt = 1
            {
              const int ipt = 1;
              float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCReco/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCReco/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }

            // ipt = 2
            {
              const int ipt = 2;
              float c2Sub = p1kBarReco[ietaA][ipt] * p1kBarReco[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCReco/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEt[ietaA][ipt] * p1kBarRecoEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCReco/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }
          }
        }

        // --- MCRecoEffCorr Covariance (1D vs eta) ---
        for (int ietaA = 1; ietaA <= (KNEta - 1) / 2; ++ietaA) {
          int ietaC = KNEta - ietaA;
          float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);

          // ipt = 0
          {
            const int ipt = 0;
            float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
            }
          }

          // ipt = 1
          {
            const int ipt = 1;
            float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
            }
          }

          // ipt = 2
          {
            const int ipt = 2;
            float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
            if (std::isfinite(c2Sub)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2Sub);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_Cov_Eta"), valy, c2Sub);
            }
            float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
            if (std::isfinite(c2SubEt)) {
              histos.fill(HIST("MCRecoEffCorr/Prof_C2EtSub_Mult_etabin_ptbin"), col.multNTracksPV(), ietaA, ipt, c2SubEt);
              histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
              if (cent < kCentCovCut) histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
            }
          }
        }

        // --- MCRecoEffCorr Covariance (2D etaA vs etaC) ---
        for (int ietaA = 1; ietaA < KNEta; ++ietaA) {
          for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
            float valx = 0.5f * (etaLw[ietaA] + etaUp[ietaA]);
            float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);

            // ipt = 0
            {
              const int ipt = 0;
              float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt0_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }

            // ipt = 1
            {
              const int ipt = 1;
              float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt1_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }

            // ipt = 2
            {
              const int ipt = 2;
              float c2Sub = p1kBarRecoEffCor[ietaA][ipt] * p1kBarRecoEffCor[ietaC][ipt];
              if (std::isfinite(c2Sub)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, c2Sub);
              float c2SubEt = p1kBarRecoEffCorEt[ietaA][ipt] * p1kBarRecoEffCorEt[ietaC][ipt];
              if (std::isfinite(c2SubEt)) histos.fill(HIST("MCRecoEffCorr/Prof_ipt2_C2EtSub2D_Mult_etaA_etaC"), cent, valx, valy, c2SubEt);
            }
          }
        }

      }
    }
    LOGF(info, "FINISHED RUNNING processMCFluc (pT + Et)");
  }
  PROCESS_SWITCH(RadFlowDecorr, processMCFluc, "process MC to calculate pt/Et fluc", cfgRunMCFluc);

  void processGetFlat(AodCollisionsSel::iterator const& coll,aod::BCsWithTimestamps const&,AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll)) return;
    float cent = getCentrality(coll);
    if (cent < 0 || cent > 80) return;
    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) continue;
      float p = track.p(); // ADDED
      float eta = track.eta();
      float phi = track.phi();
      if (p < kFloatEpsilon) continue;
      histos.fill(HIST("hCentEtaPhi"), cent, eta, phi);
      const bool isPion = selectionPion(track);
      const bool isKaon = selectionKaon(track);
      const bool isProton = selectionProton(track);
      if (isPion || isKaon || isProton) {
        histos.fill(HIST("hCentEtaPhi_PID"), cent, eta, phi);
      }
    }
  }
  PROCESS_SWITCH(RadFlowDecorr, processGetFlat,"process real data to calculate mean pT and Et", cfgRunGetFlat);

  void processDataMean(AodCollisionsSel::iterator const& coll,aod::BCsWithTimestamps const&,AodTracksSel const& tracks)
  {
    float sumWi[KNEta][KNpT]{}, sumWipti[KNEta][KNpT]{};
    float sumWiEt[KNEta][KNpT]{}, sumWiEtVal[KNEta][KNpT]{};
    if (!isEventSelected(coll)) return;
    float cent = getCentrality(coll);
    if (cent < 0 || cent > 80) return;
    histos.fill(HIST("hZvtx_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);

    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) continue;
      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      if (p < kFloatEpsilon) continue;
      histos.fill(HIST("hP"), p);
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), track.phi());

      float effIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 0);
      float fakeIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 1);
      float flatWeightIncl = getFlatteningWeight(cent, eta, phi, kInclusive);
      float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
      if (!std::isfinite(wIncl) || wIncl <= kFloatEpsilon || effIncl <= kFloatEpsilon) continue;

      histos.fill(HIST("hCentEtaPhi"), cent, eta, track.phi());
      histos.fill(HIST("hCentEtaPhiWtd"), cent, eta, track.phi(), flatWeightIncl);

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
        for (int ipt = 0; ipt < KNpT; ++ipt) {
          if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
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
        float flatWeightPid = getFlatteningWeight(cent, eta, phi, kCombinedPID);
        float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
        if (!std::isfinite(wPid) || wPid <= kFloatEpsilon || effPid <= kFloatEpsilon) continue;

        histos.fill(HIST("hCentEtaPhiWtd_PID"), cent, eta, track.phi(), flatWeightPid);
        float m = isPion ? o2::constants::physics::MassPiPlus :
        isKaon ? o2::constants::physics::MassKPlus :
        o2::constants::physics::MassProton;
        float E = std::sqrt(p * p + m * m);
        float Et = E * (pt / p); // E_T = E * sin(theta)
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
            sumWiEt[ieta][ipt] += wPid;
            sumWiEtVal[ieta][ipt] += wPid * Et;
          }
        }
      }
    }
    histos.fill(HIST("Prof_cent_Nchrec"), cent, sumWi[0][0]);
    if (isfinite(sumWi[0][0])) histos.fill(HIST("Prof_MeanpT_Cent"), cent, sumWipti[0][0] / sumWi[0][0]);
    if (isfinite(sumWiEt[0][0])) histos.fill(HIST("Prof_MeanEt_Cent"), cent, sumWiEtVal[0][0] / sumWiEt[0][0]);

    for (int ieta = 0; ieta < KNEta; ++ieta){
      for (int ipt = 0; ipt < KNpT; ++ipt) {
        if (isfinite(sumWi[ieta][ipt]))
        histos.fill(HIST("pmean_nch_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt,
        sumWipti[ieta][ipt] / sumWi[ieta][ipt]);
        if (isfinite(sumWiEt[ieta][ipt]))
        histos.fill(HIST("pmeanEt_nch_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt,
        sumWiEtVal[ieta][ipt] / sumWiEt[ieta][ipt]);
      }
    }
  }
  PROCESS_SWITCH(RadFlowDecorr, processDataMean,"process real data to calculate mean pT and Et", cfgRunDataMean);


  void processDataFluc(AodCollisionsSel::iterator const& coll,aod::BCsWithTimestamps const&,AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll)) return;
    float cent = getCentrality(coll);
    if (cent < 0 || cent > 80) return;
    if (!pmeanNchEtabinPtbinStep2 || !pmeanEtNchEtabinPtbinStep2) {
      LOGF(warning, "Data fluc: Mean pT or Et map missing");
      return;
    }

    if (!hEff[kInclusive] || !hFake[kInclusive] || !hWeightMap3D[kInclusive] ||!hEff[kCombinedPID] || !hFake[kCombinedPID] || !hWeightMap3D[kCombinedPID]) {
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
      if (!isTrackSelected(track)) continue;
      float pt = track.pt();
      float eta = track.eta();
      float p = track.p();
      float phi = track.phi();
      if (p < kFloatEpsilon) continue;

      float effIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 0);
      float fakeIncl = getEfficiency(coll.multNTracksPV(), pt, eta, kInclusive, 1);
      float flatWeightIncl = getFlatteningWeight(cent, eta, phi, kInclusive);

      float wIncl = flatWeightIncl * (1.0 - fakeIncl) / effIncl;
      if (!std::isfinite(wIncl) || wIncl <= kFloatEpsilon || effIncl <= kFloatEpsilon) continue;

      for (int ieta = 0; ieta < KNEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
        for (int ipt = 0; ipt < KNpT; ++ipt) {
          if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
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
        float flatWeightPid = getFlatteningWeight(cent, eta, phi, kCombinedPID);

        float wPid = flatWeightPid * (1.0 - fakePid) / effPid;
        if (!std::isfinite(wPid) || wPid <= kFloatEpsilon || effPid <= kFloatEpsilon) continue;

        float m = isPion ? o2::constants::physics::MassPiPlus :
        isKaon ? o2::constants::physics::MassKPlus :
        o2::constants::physics::MassProton;

        float E = std::sqrt(p * p + m * m);
        float Et = E * (pt / p); // E_T = E * sin(theta)
        for (int ieta = 0; ieta < KNEta; ++ieta) {
          if (eta <= etaLw[ieta] || eta > etaUp[ieta]) continue;
          for (int ipt = 0; ipt < KNpT; ++ipt) {
            if (pt <= pTLw[ipt] || pt > pTUp[ipt]) continue;
            for (int k = 0; k < KIntK; ++k) {
              for (int m = 0; m < KIntM; ++m)
              sumpmwkEt[ieta][ipt][m][k] += std::pow(wPid, k) * std::pow(Et, m);
              sumwkEt[ieta][ipt][k] += std::pow(wPid, k);
            }
          }
        }
      }
    }
    for (int ieta = 0; ieta < KNEta; ++ieta)
    for (int ipt = 0; ipt < KNpT; ++ipt) {
      const int ibx = pmeanNchEtabinPtbinStep2->GetXaxis()->FindBin(coll.multNTracksPV());
      const int iby = ieta + 1;
      const int ibz = ipt + 1;
      float mmpt = pmeanNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);
      float mmet = pmeanEtNchEtabinPtbinStep2->GetBinContent(ibx, iby, ibz);

      mean[ieta][ipt] = sumpmwk[ieta][ipt][1][1]/sumwk[ieta][ipt][1];
      meanEt[ieta][ipt] = sumpmwkEt[ieta][ipt][1][1]/sumwkEt[ieta][ipt][1];

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

    if (std::isfinite(c2[0][0])) histos.fill(HIST("Prof_C2_Cent"), cent, c2[0][0]);
    if (std::isfinite(c2Et[0][0])) histos.fill(HIST("Prof_C2Et_Cent"), cent, c2Et[0][0]);
    if (std::isfinite(sumwk[0][0][1])) {
      histos.fill(HIST("Prof_MeanpT_Cent"), cent, mean[0][0]);
      histos.fill(HIST("Prof_MeanpT_Mult"), coll.multNTracksPV(), mean[0][0]);
    }
    if (std::isfinite(sumwkEt[0][0][1])) {
      histos.fill(HIST("Prof_MeanEt_Cent"), cent, meanEt[0][0]);
      histos.fill(HIST("Prof_MeanEt_Mult"), coll.multNTracksPV(), meanEt[0][0]);
    }

    for (int ieta = 0; ieta < KNEta; ++ieta)
    for (int ipt = 0; ipt < KNpT; ++ipt) {
      if (std::isfinite(c2[ieta][ipt])) histos.fill(HIST("Prof_C2_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, c2[ieta][ipt]);
      if (std::isfinite(c2Et[ieta][ipt])) histos.fill(HIST("Prof_C2Et_Mult_etabin_ptbin"), coll.multNTracksPV(), ieta, ipt, c2Et[ieta][ipt]);
    }

    for (int ietaA = 1; ietaA <= (KNEta-1)/2; ++ietaA) {
      int ietaC = KNEta - ietaA;
      float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);

      {
        const int ipt = 0;
        float c2Sub = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
        if (std::isfinite(c2Sub)) {
          histos.fill(HIST("Prof_ipt0_Cov_Cent_eta"), cent, valy, c2Sub);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt0_Cov_Eta"), valy, c2Sub);
        }
        float c2SubEt = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
        if (std::isfinite(c2SubEt)) {
          histos.fill(HIST("Prof_ipt0_CovEt_Cent_eta"), cent, valy, c2SubEt);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt0_CovEt_Eta"), valy, c2SubEt);
        }
      }

      {
        const int ipt = 1;
        float c2Sub = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
        if (std::isfinite(c2Sub)) {
          histos.fill(HIST("Prof_ipt1_Cov_Cent_eta"), cent, valy, c2Sub);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt1_Cov_Eta"), valy, c2Sub);
        }
        float c2SubEt = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
        if (std::isfinite(c2SubEt)) {
          histos.fill(HIST("Prof_ipt1_CovEt_Cent_eta"), cent, valy, c2SubEt);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt1_CovEt_Eta"), valy, c2SubEt);
        }
      }

      {
        const int ipt = 2;
        float c2Sub = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
        if (std::isfinite(c2Sub)) {
          histos.fill(HIST("Prof_ipt2_Cov_Cent_eta"), cent, valy, c2Sub);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt2_Cov_Eta"), valy, c2Sub);
        }
        float c2SubEt = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
        if (std::isfinite(c2SubEt)) {
          histos.fill(HIST("Prof_ipt2_CovEt_Cent_eta"), cent, valy, c2SubEt);
          if (cent < kCentCovCut) histos.fill(HIST("Prof_ipt2_CovEt_Eta"), valy, c2SubEt);
        }
      }
    }

    for (int ietaA = 1; ietaA < KNEta; ++ietaA){
      for (int ietaC = 1; ietaC < KNEta; ++ietaC) {
        float valx = 0.5f * (etaLw[ietaA] + etaUp[ietaA]);
        float valy = kHalf * (etaLw[ietaC] + etaUp[ietaC]);
        {
          const int ipt = 0;
          float covpt = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
          if (std::isfinite(covpt)) histos.fill(HIST("Prof_ipt0_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);

          float covet = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
          if (std::isfinite(covet)) histos.fill(HIST("Prof_ipt0_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
        }


        {
          const int ipt = 1;
          float covpt = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
          if (std::isfinite(covpt)) histos.fill(HIST("Prof_ipt1_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);

          float covet = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
          if (std::isfinite(covet)) histos.fill(HIST("Prof_ipt1_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
        }
        {
          const int ipt = 2;

          float covpt = p1kBar[ietaA][ipt] * p1kBar[ietaC][ipt];
          if (std::isfinite(covpt)) histos.fill(HIST("Prof_ipt2_C2Sub2D_Mult_etaA_etaC"), cent, valx, valy, covpt);

          float covet = p1kBarEt[ietaA][ipt] * p1kBarEt[ietaC][ipt];
          if (std::isfinite(covet)) histos.fill(HIST("Prof_ipt2_C2SubEt2D_Mult_etaA_etaC"), cent, valx, valy, covet);
        }
      }
    }
  }
  PROCESS_SWITCH(RadFlowDecorr, processDataFluc,"process real data to calculate fluc pT and Et", cfgRunDataFluc);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<RadFlowDecorr>(cfgc)};
  return workflow;
}
