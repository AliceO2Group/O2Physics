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
//
// Example analysis task to select clean V0 sample
// ========================
//
// This code loops over a V0Data table and produces some standard analysis output.
//
//    Comments, questions, complaints, suggestions?
//    Please write to: daiki.sekihata@cern.ch
//
#include <array>
#include <map>
#include <string>
#include <memory>
#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "DCAFitter/DCAFitterN.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

/*using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi,
                                aod::pidTPCFullKa, aod::pidTPCFullPr>;*/
using FullTracksExt = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi,
                                aod::pidTPCFullKa, aod::pidTPCFullPr>;

constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackTPCPID;

struct v0selector {

  // Configurables for curved QT cut
  //  Gamma cuts
  Configurable<float> cutAlphaG{"cutAlphaG", 0.4, "cutAlphaG"};
  Configurable<float> cutQTG{"cutQTG", 0.03, "cutQTG"};
  Configurable<float> cutAlphaGLow{"cutAlphaGLow", 0.4, "cutAlphaGLow"};
  Configurable<float> cutAlphaGHigh{"cutAlphaGHigh", 0.8, "cutAlphaGHigh"};
  Configurable<float> cutQTG2{"cutQTG2", 0.02, "cutQTG2"};
  // K0S cuts
  Configurable<float> cutQTK0SLow{"cutQTK0SLow", 0.1075, "cutQTK0SLow"};
  Configurable<float> cutQTK0SHigh{"cutQTK0SHigh", 0.215, "cutQTK0SHigh"};
  Configurable<float> cutAPK0SLow{"cutAPK0SLow", 0.199, "cutAPK0SLow"};
  Configurable<float> cutAPK0SHigh{"cutAPK0SHigh", 0.8, "cutAPK0SHigh"};
  // Lambda & A-Lambda cuts
  Configurable<float> cutQTL{"cutQTL", 0.03, "cutQTL"};
  Configurable<float> cutAlphaLLow{"cutAlphaLLow", 0.35, "cutAlphaLLow"};
  Configurable<float> cutAlphaLHigh{"cutAlphaLHigh", 0.7, "cutAlphaLHigh"};
  Configurable<float> cutAlphaALLow{"cutAlphaALow", -0.7, "cutAlphaALow"};
  Configurable<float> cutAlphaALHigh{"cutAlphaAHigh", -0.35, "cutAlphaAHigh"};
  Configurable<float> cutAPL1{"cutAPL1", 0.107, "cutAPL1"};
  Configurable<float> cutAPL2{"cutAPL2", -0.69, "cutAPL2"};
  Configurable<float> cutAPL3{"cutAPL3", 0.5, "cutAPL3"};

  enum { // Reconstructed V0
    kUndef = -1,
    kGamma = 0,
    kK0S = 1,
    kLambda = 2,
    kAntiLambda = 3,
    kOmega = 4
  };

  Produces<o2::aod::V0Bits> v0bits;

  // int checkV0(const array<float, 3>& ppos, const array<float, 3>& pneg)
  int checkV0(const float alpha, const float qt)
  {
    // float alpha = alphav0(ppos, pneg);
    // float qt = qtarmv0(ppos, pneg);
    // // Gamma cuts
    // const float cutAlphaG = 0.4;
    // const float cutQTG = 0.03;
    // const float cutAlphaG2[2] = {0.4, 0.8};
    // const float cutQTG2 = 0.02;

    // // K0S cuts
    // const float cutQTK0S[2] = {0.1075, 0.215};
    // const float cutAPK0S[2] = {0.199, 0.8}; // parameters for curved QT cut

    // // Lambda & A-Lambda cuts
    // const float cutQTL = 0.03;
    // const float cutAlphaL[2] = {0.35, 0.7};
    // const float cutAlphaAL[2] = {-0.7, -0.35};
    // const float cutAPL[3] = {0.107, -0.69, 0.5}; // parameters for curved QT cut

    if (qt < cutQTG) {
      if ((TMath::Abs(alpha) < cutAlphaG)) {
        return kGamma;
      }
    }
    if (qt < cutQTG2) {
      // additional region - should help high pT gammas
      if ((TMath::Abs(alpha) > cutAlphaGLow) && (TMath::Abs(alpha) < cutAlphaGHigh)) {
        return kGamma;
      }
    }

    // Check for K0S candidates
    float q = cutAPK0SLow * TMath::Sqrt(TMath::Abs(1 - alpha * alpha / (cutAPK0SHigh * cutAPK0SHigh)));
    if ((qt > cutQTK0SLow) && (qt < cutQTK0SHigh) && (qt > q)) {
      return kK0S;
    }

    // Check for Lambda candidates
    q = cutAPL1 * TMath::Sqrt(TMath::Abs(1 - ((alpha + cutAPL2) * (alpha + cutAPL2)) / (cutAPL3 * cutAPL3)));
    if ((alpha > cutAlphaLLow) && (alpha < cutAlphaLHigh) && (qt > cutQTL) && (qt < q)) {
      return kLambda;
    }

    // Check for AntiLambda candidates
    q = cutAPL1 * TMath::Sqrt(TMath::Abs(1 - ((alpha - cutAPL2) * (alpha - cutAPL2)) / (cutAPL3 * cutAPL3)));
    if ((alpha > cutAlphaALLow) && (alpha < cutAlphaALHigh) && (qt > cutQTL) && (qt < q)) {
      return kAntiLambda;
    }

    return kUndef;
  }

  // Configurables
  Configurable<float> v0max_mee{"v0max_mee", 0.1, "max mee for photon"};
  Configurable<float> maxpsipair{"maxpsipair", 1.6, "max psi_pair for photon"};
  Configurable<float> v0cospa{"v0cospa", 0.998, "V0 CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 0.3, "DCA V0 Daughters"};
  Configurable<float> v0Rmin{"v0Rmin", 0.0, "v0Rmin"};
  Configurable<float> v0Rmax{"v0Rmax", 90.0, "v0Rmax"};
  Configurable<float> dcamin{"dcamin", 0.0, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<bool> fillhisto{"fillhisto", false, "flag to fill histograms"};
  // cutNsigmaElTPC, cutNsigmaPiTPC, cutNsigmaPrTPC
  Configurable<float> cutNsigmaElTPC{"cutNsigmaElTPC", 5.0, "cutNsigmaElTPC"};
  Configurable<float> cutNsigmaPiTPC{"cutNsigmaPiTPC", 5.0, "cutNsigmaPiTPC"};
  Configurable<float> cutNsigmaPrTPC{"cutNsigmaPrTPC", 5.0, "cutNsigmaPrTPC"};

  HistogramRegistry registry{"registry"};
  void init(o2::framework::InitContext&)
  {
    if (fillhisto) {
      registry.add("hV0Candidate", "hV0Candidate", HistType::kTH1F, {{2, 0.5f, 2.5f}});
      registry.add("hMassGamma", "hMassGamma", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 0.0f, 0.1f}});
      registry.add("hGammaRxy", "hGammaRxy", HistType::kTH2F, {{1800, -90.0f, 90.0f}, {1800, -90.0f, 90.0f}});
      registry.add("hMassK0S", "hMassK0S", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 0.45, 0.55}});
      registry.add("hMassK0SPt", "hMassK0SPt", HistType::kTH2F, {{200, 0.0f, 20.0f}, {100, 0.45, 0.55}});
      registry.add("hMassK0SEta", "hMassK0SEta", HistType::kTH2F, {{20, -1, 1}, {100, 0.45, 0.55}});
      registry.add("hMassK0SPhi", "hMassK0SPhi", HistType::kTH2F, {{63, 0, 6.3}, {100, 0.45, 0.55}});
      registry.add("hMassLambda", "hMassLambda", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 1.05, 1.15f}});
      registry.add("hMassAntiLambda", "hAntiMassLambda", HistType::kTH2F, {{900, 0.0f, 90.0f}, {100, 1.05, 1.15f}});
      registry.add("hV0Pt", "pT", HistType::kTH1F, {{100, 0.0f, 10}});
      registry.add("hV0EtaPhi", "#eta vs. #varphi", HistType::kTH2F, {{63, 0, 6.3}, {20, -1.0f, 1.0f}});
      registry.add("hV0Radius", "hV0Radius", HistType::kTH1F, {{1000, 0.0f, 100.0f}});
      registry.add("hV0CosPA", "hV0CosPA", HistType::kTH1F, {{50, 0.95f, 1.0f}});
      registry.add("hDCAxyPosToPV", "hDCAxyPosToPV", HistType::kTH1F, {{1000, -5.0f, 5.0f}});
      registry.add("hDCAxyNegToPV", "hDCAxyNegToPV", HistType::kTH1F, {{1000, -5.0f, 5.0f}});
      registry.add("hDCAzPosToPV", "hDCAzPosToPV", HistType::kTH1F, {{1000, -5.0f, 5.0f}});
      registry.add("hDCAzNegToPV", "hDCAzNegToPV", HistType::kTH1F, {{1000, -5.0f, 5.0f}});
      registry.add("hDCAV0Dau", "hDCAV0Dau", HistType::kTH1F, {{1000, 0.0f, 10.0f}});
      registry.add("hV0APplot", "hV0APplot", HistType::kTH2F, {{200, -1.0f, +1.0f}, {250, 0.0f, 0.25f}});
      registry.add("hV0Psi", "hV0Psi", HistType::kTH2F, {{100, 0, TMath::PiOver2()}, {100, 0, 0.1}});
    }
  }

  void process(aod::V0Datas const& V0s, FullTracksExt const& tracks, aod::Collisions const&)
  {
    std::map<int, uint8_t> pidmap;

    for (auto& V0 : V0s) {
      // if (!(V0.posTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
      //   continue;
      // }
      // if (!(V0.negTrack_as<FullTracksExt>().trackType() & o2::aod::track::TPCrefit)) {
      //   continue;
      // }

      // printf("V0.collisionId = %d , collision.globalIndex = %d\n",V0.collisionId(),collision.globalIndex());
      if (fillhisto) {
        registry.fill(HIST("hV0Candidate"), 1);
      }
      if (fabs(V0.posTrack_as<FullTracksExt>().eta()) > 0.9) {
        continue;
      }
      if (fabs(V0.negTrack_as<FullTracksExt>().eta()) > 0.9) {
        continue;
      }

      if (V0.posTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      if (V0.negTrack_as<FullTracksExt>().tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }

      if (V0.posTrack_as<FullTracksExt>().tpcChi2NCl() > maxchi2tpc) {
        continue;
      }
      if (V0.negTrack_as<FullTracksExt>().tpcChi2NCl() > maxchi2tpc) {
        continue;
      }

      if (fabs(V0.posTrack_as<FullTracksExt>().dcaXY()) < dcamin) {
        continue;
      }
      if (fabs(V0.negTrack_as<FullTracksExt>().dcaXY()) < dcamin) {
        continue;
      }

      if (fabs(V0.posTrack_as<FullTracksExt>().dcaXY()) > dcamax) {
        continue;
      }
      if (fabs(V0.negTrack_as<FullTracksExt>().dcaXY()) > dcamax) {
        continue;
      }

      if (V0.posTrack_as<FullTracksExt>().sign() * V0.negTrack_as<FullTracksExt>().sign() > 0) { // reject same sign pair
        continue;
      }

      // if (V0.posTrack_as<FullTracksExt>().collisionId() != V0.negTrack_as<FullTracksExt>().collisionId()) {
      //   continue;
      // }

      // if (!V0.posTrack_as<FullTracksExt>().has_collision() || !V0.negTrack_as<FullTracksExt>().has_collision()) {
      //   continue;
      // }

      // auto const& collision = V0.collision_as<aod::Collisions>();

      //      if (V0.collisionId() != collision.globalIndex()) {
      //        continue;
      //      }

      float V0dca = V0.dcaV0daughters();
      float V0CosinePA = V0.v0cosPA();
      float V0radius = V0.v0radius();

      if (V0dca > dcav0dau) {
        continue;
      }

      if (V0CosinePA < v0cospa) {
        continue;
      }

      if (V0radius < v0Rmin || v0Rmax < V0radius) {
        continue;
      }
      if (fillhisto) {
        registry.fill(HIST("hV0Candidate"), 2);
      }
      float mGamma = V0.mGamma();
      float mK0S = V0.mK0Short();
      float mLambda = V0.mLambda();
      float mAntiLambda = V0.mAntiLambda();
      float psipair = V0.psipair();

      if (fillhisto) {
        registry.fill(HIST("hV0Pt"), V0.pt());
        registry.fill(HIST("hV0EtaPhi"), V0.phi(), V0.eta());
        registry.fill(HIST("hDCAxyPosToPV"), V0.posTrack_as<FullTracksExt>().dcaXY());
        registry.fill(HIST("hDCAxyNegToPV"), V0.negTrack_as<FullTracksExt>().dcaXY());
        registry.fill(HIST("hDCAzPosToPV"), V0.posTrack_as<FullTracksExt>().dcaZ());
        registry.fill(HIST("hDCAzNegToPV"), V0.negTrack_as<FullTracksExt>().dcaZ());
        registry.fill(HIST("hV0APplot"), V0.alpha(), V0.qtarm());
        registry.fill(HIST("hV0Radius"), V0radius);
        registry.fill(HIST("hV0CosPA"), V0CosinePA);
        registry.fill(HIST("hDCAV0Dau"), V0dca);
      }

      int v0id = checkV0(V0.alpha(), V0.qtarm());
      if (v0id < 0) {
        // printf("This is not [Gamma/K0S/Lambda/AntiLambda] candidate.\n");
        continue;
      }

      if (v0id == kGamma) { // photon conversion
        if (fillhisto) {
          registry.fill(HIST("hMassGamma"), V0radius, mGamma);
          registry.fill(HIST("hV0Psi"), psipair, mGamma);
        }
        if (mGamma < v0max_mee && TMath::Abs(V0.posTrack_as<FullTracksExt>().tpcNSigmaEl()) < cutNsigmaElTPC && TMath::Abs(V0.negTrack_as<FullTracksExt>().tpcNSigmaEl()) < cutNsigmaElTPC && TMath::Abs(psipair) < maxpsipair) {
          pidmap[V0.posTrackId()] |= (uint8_t(1) << kGamma);
          pidmap[V0.negTrackId()] |= (uint8_t(1) << kGamma);
          if (fillhisto) {
            registry.fill(HIST("hGammaRxy"), V0.x(), V0.y());
          }
        }
      } else if (v0id == kK0S) { // K0S-> pi pi
        if (fillhisto) {
          registry.fill(HIST("hMassK0S"), V0radius, mK0S);
          registry.fill(HIST("hMassK0SPt"), V0.pt(), mK0S);
          registry.fill(HIST("hMassK0SEta"), V0.eta(), mK0S);
          registry.fill(HIST("hMassK0SPhi"), V0.phi(), mK0S);
        }
        if ((0.48 < mK0S && mK0S < 0.51) && TMath::Abs(V0.posTrack_as<FullTracksExt>().tpcNSigmaPi()) < cutNsigmaPiTPC && TMath::Abs(V0.negTrack_as<FullTracksExt>().tpcNSigmaPi()) < cutNsigmaPiTPC) {
          pidmap[V0.posTrackId()] |= (uint8_t(1) << kK0S);
          pidmap[V0.negTrackId()] |= (uint8_t(1) << kK0S);
        }
      } else if (v0id == kLambda) { // L->p + pi-
        if (fillhisto) {
          registry.fill(HIST("hMassLambda"), V0radius, mLambda);
        }
        if (v0id == kLambda && (1.110 < mLambda && mLambda < 1.120) && TMath::Abs(V0.posTrack_as<FullTracksExt>().tpcNSigmaPr()) < cutNsigmaPrTPC && TMath::Abs(V0.negTrack_as<FullTracksExt>().tpcNSigmaPi()) < cutNsigmaPiTPC) {
          pidmap[V0.posTrackId()] |= (uint8_t(1) << kLambda);
          pidmap[V0.negTrackId()] |= (uint8_t(1) << kLambda);
        }
      } else if (v0id == kAntiLambda) { // Lbar -> pbar + pi+
        if (fillhisto) {
          registry.fill(HIST("hMassAntiLambda"), V0radius, mAntiLambda);
        }
        if ((1.110 < mAntiLambda && mAntiLambda < 1.120) && TMath::Abs(V0.posTrack_as<FullTracksExt>().tpcNSigmaPi()) < cutNsigmaPiTPC && TMath::Abs(V0.negTrack_as<FullTracksExt>().tpcNSigmaPr()) < cutNsigmaPrTPC) {
          pidmap[V0.posTrackId()] |= (uint8_t(1) << kAntiLambda);
          pidmap[V0.negTrackId()] |= (uint8_t(1) << kAntiLambda);
        }
      }

      // printf("posTrackId = %d\n",V0.posTrackId());
      // printf("negTrackId = %d\n",V0.negTrackId());

    } // end of V0 loop

    for (auto& track : tracks) {
      // printf("setting pidmap[%lld] = %d\n",track.globalIndex(),pidmap[track.globalIndex()]);
      v0bits(pidmap[track.globalIndex()]);
    } // end of track loop

  } // end of process
};

struct trackPIDQA {
  Configurable<float> dcamin{"dcamin", 0.0, "dcamin"};
  Configurable<float> dcamax{"dcamax", 1e+10, "dcamax"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 4.0, "max chi2/NclsTPC"};
  Configurable<bool> fillDQHisto{"fillDQHisto", false, "flag to fill dq histograms"};
  Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of dq histograms"};

  HistogramRegistry registry{"registry"};
  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  HistogramManager* fHistMan;
  void init(o2::framework::InitContext& context)
  {
    bool enableBarrelHistos = context.mOptions.get<bool>("processQA");
    if (enableBarrelHistos) {
      registry.add("hEventCounter", "hEventCounter", HistType::kTH1F, {{5, 0.5f, 5.5f}});
      registry.add("hTrackPt_all", "pT", HistType::kTH1F, {{100, 0.0, 10}});
      registry.add("hTrackEtaPhi_all", "#eta vs. #varphi", HistType::kTH2F, {{63, 0, 6.3}, {20, -1.0f, 1.0f}});
      registry.add("h2TPCdEdx_Pin_all", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});
      registry.add("h2TOFbeta_Pin_all", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});

      registry.add("hTrackPt", "pT", HistType::kTH1F, {{100, 0.0, 10}});
      registry.add("hTrackEtaPhi", "#eta vs. #varphi", HistType::kTH2F, {{63, 0, 6.3}, {20, -1.0f, 1.0f}});

      registry.add("h2TPCdEdx_Pin", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});
      registry.add("h2TPCdEdx_Pin_El", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});
      registry.add("h2TPCdEdx_Pin_Pi", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});
      registry.add("h2TPCdEdx_Pin_Ka", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});
      registry.add("h2TPCdEdx_Pin_Pr", "TPC dEdx vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, 0.0, 200.}});

      registry.add("h2TPCnSigma_Pin_El", "TPC n#sigma_{e} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TPCnSigma_Pin_Pi", "TPC n#sigma_{#pi} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TPCnSigma_Pin_Ka", "TPC n#sigma_{K} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TPCnSigma_Pin_Pr", "TPC n#sigma_{p} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});

      registry.add("h2TOFbeta_Pin", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});
      registry.add("h2TOFbeta_Pin_El", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});
      registry.add("h2TOFbeta_Pin_Pi", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});
      registry.add("h2TOFbeta_Pin_Ka", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});
      registry.add("h2TOFbeta_Pin_Pr", "TOF #beta vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {120, 0.0, 1.2}});

      registry.add("h2TOFnSigma_Pin_El", "TOF n#sigma_{e} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TOFnSigma_Pin_Pi", "TOF n#sigma_{#pi} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TOFnSigma_Pin_Ka", "TOF n#sigma_{K} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
      registry.add("h2TOFnSigma_Pin_Pr", "TOF n#sigma_{p} vs. p_{in}", HistType::kTH2F, {{1000, 0.0, 10}, {200, -10, +10}});
    }

    if (fillDQHisto) {
      TString histClasses = "";
      VarManager::SetDefaultVarNames();
      fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
      fHistMan->SetUseDefaultVariableNames(kTRUE);
      fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);
      histClasses += "V0Track_all;";
      histClasses += "V0Track_electron;";
      histClasses += "V0Track_pion;";
      histClasses += "V0Track_proton;";
      DefineHistograms(histClasses);
      VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
      fOutputList.setObject(fHistMan->GetMainHistogramList());
    }
  }

  void processQA(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
                 soa::Join<FullTracksExt, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta, aod::V0Bits> const& tracks)
  {
    registry.fill(HIST("hEventCounter"), 1.0); // all

    if (!collision.sel8()) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 2.0); // FT0VX i.e. FT0and

    if (collision.numContrib() < 0.5) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 3.0); // Ncontrib > 0

    if (abs(collision.posZ()) > 10.0) {
      return;
    }
    registry.fill(HIST("hEventCounter"), 4.0); //|Zvtx| < 10 cm

    for (auto& track : tracks) {
      if (!track.has_collision()) {
        continue;
      }

      if (track.tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }

      if (track.tpcChi2NCl() > maxchi2tpc) {
        continue;
      }

      if (fabs(track.dcaXY()) < dcamin) {
        continue;
      }

      if (fabs(track.dcaXY()) > dcamax) {
        continue;
      }

      if (fabs(track.eta()) > 0.9) {
        continue;
      }

      registry.fill(HIST("hTrackPt_all"), track.pt());
      registry.fill(HIST("hTrackEtaPhi_all"), track.phi(), track.eta());
      registry.fill(HIST("h2TPCdEdx_Pin_all"), track.tpcInnerParam(), track.tpcSignal());
      registry.fill(HIST("h2TOFbeta_Pin_all"), track.tpcInnerParam(), track.beta());
      if (track.pidbit() > 0) {
        registry.fill(HIST("hTrackPt"), track.pt());
        registry.fill(HIST("hTrackEtaPhi"), track.phi(), track.eta());
        registry.fill(HIST("h2TPCdEdx_Pin"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("h2TOFbeta_Pin"), track.tpcInnerParam(), track.beta());
      }

      if (static_cast<bool>(track.pidbit() & (1 << v0selector::kGamma))) {
        registry.fill(HIST("h2TPCdEdx_Pin_El"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("h2TOFbeta_Pin_El"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("h2TPCnSigma_Pin_El"), track.tpcInnerParam(), track.tpcNSigmaEl());
        registry.fill(HIST("h2TOFnSigma_Pin_El"), track.tpcInnerParam(), track.tofNSigmaEl());
      }
      if (static_cast<bool>(track.pidbit() & (1 << v0selector::kK0S))) {
        registry.fill(HIST("h2TPCdEdx_Pin_Pi"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("h2TOFbeta_Pin_Pi"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("h2TPCnSigma_Pin_Pi"), track.tpcInnerParam(), track.tpcNSigmaPi());
        registry.fill(HIST("h2TOFnSigma_Pin_Pi"), track.tpcInnerParam(), track.tofNSigmaPi());
      }
      if (static_cast<bool>(track.pidbit() & (1 << v0selector::kLambda))) {
        if (track.sign() > 0) {
          registry.fill(HIST("h2TPCdEdx_Pin_Pr"), track.tpcInnerParam(), track.tpcSignal());
          registry.fill(HIST("h2TOFbeta_Pin_Pr"), track.tpcInnerParam(), track.beta());
          registry.fill(HIST("h2TPCnSigma_Pin_Pr"), track.tpcInnerParam(), track.tpcNSigmaPr());
          registry.fill(HIST("h2TOFnSigma_Pin_Pr"), track.tpcInnerParam(), track.tofNSigmaPr());
        } else {
          registry.fill(HIST("h2TPCdEdx_Pin_Pi"), track.tpcInnerParam(), track.tpcSignal());
          registry.fill(HIST("h2TOFbeta_Pin_Pi"), track.tpcInnerParam(), track.beta());
          registry.fill(HIST("h2TPCnSigma_Pin_Pi"), track.tpcInnerParam(), track.tpcNSigmaPi());
          registry.fill(HIST("h2TOFnSigma_Pin_Pi"), track.tpcInnerParam(), track.tofNSigmaPi());
        }
      }
      if (static_cast<bool>(track.pidbit() & (1 << v0selector::kAntiLambda))) {
        if (track.sign() > 0) {
          registry.fill(HIST("h2TPCdEdx_Pin_Pi"), track.tpcInnerParam(), track.tpcSignal());
          registry.fill(HIST("h2TOFbeta_Pin_Pi"), track.tpcInnerParam(), track.beta());
          registry.fill(HIST("h2TPCnSigma_Pin_Pi"), track.tpcInnerParam(), track.tpcNSigmaPi());
          registry.fill(HIST("h2TOFnSigma_Pin_Pi"), track.tpcInnerParam(), track.tofNSigmaPi());
        } else {
          registry.fill(HIST("h2TPCdEdx_Pin_Pr"), track.tpcInnerParam(), track.tpcSignal());
          registry.fill(HIST("h2TOFbeta_Pin_Pr"), track.tpcInnerParam(), track.beta());
          registry.fill(HIST("h2TPCnSigma_Pin_Pr"), track.tpcInnerParam(), track.tpcNSigmaPr());
          registry.fill(HIST("h2TOFnSigma_Pin_Pr"), track.tpcInnerParam(), track.tofNSigmaPr());
        }
      }
      if (static_cast<bool>(track.pidbit() & (1 << v0selector::kOmega))) {
        registry.fill(HIST("h2TPCdEdx_Pin_Ka"), track.tpcInnerParam(), track.tpcSignal());
        registry.fill(HIST("h2TOFbeta_Pin_Ka"), track.tpcInnerParam(), track.beta());
        registry.fill(HIST("h2TPCnSigma_Pin_Ka"), track.tpcInnerParam(), track.tpcNSigmaKa());
        registry.fill(HIST("h2TOFnSigma_Pin_Ka"), track.tpcInnerParam(), track.tofNSigmaKa());
      }

      if (fillDQHisto) {
        VarManager::FillTrack<gkTrackFillMap>(track);
        fHistMan->FillHistClass("V0Track_all", VarManager::fgValues);
        if (static_cast<bool>(track.pidbit() & (1 << v0selector::kGamma))) {
          fHistMan->FillHistClass("V0Track_electron", VarManager::fgValues);
        }
        if (static_cast<bool>(track.pidbit() & (1 << v0selector::kK0S))) {
          fHistMan->FillHistClass("V0Track_pion", VarManager::fgValues);
        }
        if (static_cast<bool>(track.pidbit() & (1 << v0selector::kLambda))) {
          if (track.sign() > 0) {
            fHistMan->FillHistClass("V0Track_proton", VarManager::fgValues);
          } else {
            fHistMan->FillHistClass("V0Track_pion", VarManager::fgValues);
          }
        }
        if (static_cast<bool>(track.pidbit() & (1 << v0selector::kAntiLambda))) {
          if (track.sign() > 0) {
            fHistMan->FillHistClass("V0Track_pion", VarManager::fgValues);
          } else {
            fHistMan->FillHistClass("V0Track_proton", VarManager::fgValues);
          }
        }
      }

    } // end of track loop
  }   // end of process

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      fHistMan->AddHistClass(classStr.Data());

      // fill the THn histograms
      if (classStr.Contains("V0Track_electron")) {
        o2::aod::dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_electron");
      }
      if (classStr.Contains("V0Track_pion")) {
        o2::aod::dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_pion");
      }
      if (classStr.Contains("V0Track_proton")) {
        o2::aod::dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", "postcalib_proton");
      }

      TString histTrackName = fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        o2::aod::dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
      }
    }
  }

  void processDummy(soa::Join<aod::Collisions, aod::EvSels>::iterator const&)
  {
    // do nothing
  }

  PROCESS_SWITCH(trackPIDQA, processQA, "Run PID QA for barrel tracks", true);
  PROCESS_SWITCH(trackPIDQA, processDummy, "Dummy function", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0selector>(cfgc, TaskName{"v0-selector"}),
    adaptAnalysisTask<trackPIDQA>(cfgc, TaskName{"track-pid-qa"})};
}
