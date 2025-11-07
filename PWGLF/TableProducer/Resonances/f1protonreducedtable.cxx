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

/// \file f1protonreducedtable.cxx
/// \brief Selection of events with triplets and pairs for femtoscopic studies
///
/// \author Sourav Kundu, sourav.kundu@cern.ch

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/ReducedF1ProtonTables.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include <Framework/Configurable.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h>
#include <TDatabasePDG.h> // FIXME
#include <TMath.h>
#include <TPDGCode.h> // FIXME

#include <fairlogger/Logger.h>

#include <iostream>
#include <iterator>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct f1protonreducedtable {

  // Produce derived tables
  Produces<aod::RedF1PEvents> redf1pevents;
  Produces<aod::F1Tracks> f1track;
  Produces<aod::ProtonTracks> protontrack;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Configs for events
  Configurable<bool> ConfEvtSelectZvtx{"ConfEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};

  // event spherocity calculation
  Configurable<int> trackSphDef{"trackSphDef", 0, "Spherocity Definition: |pT| = 1 -> 0, otherwise -> 1"};
  Configurable<int> trackSphMin{"trackSphMin", 10, "Number of tracks for Spherocity Calculation"};

  // Configs for track PID
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", true, "Analysed skimmed events"};
  Configurable<bool> ConfUseManualPIDproton{"ConfUseManualPIDproton", true, "True: use home-made PID solution for proton "};
  Configurable<bool> ConfUseManualPIDkaon{"ConfUseManualPIDkaon", true, "True: use home-made PID solution for kaon "};
  Configurable<bool> ConfUseManualPIDpion{"ConfUseManualPIDpion", true, "True: use home-made PID solution for pion "};
  Configurable<bool> ConfUseManualPIDdaughterPion{"ConfUseManualPIDdaughterPion", true, "True: use home-made PID solution for pion from V0"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "ccdb-url"};
  Configurable<std::string> ConfPIDBBProton{"ConfPIDBBProton", "Users/a/ariedel/FinalTrigger/PIDProton", "Path to the CCDB ocject for proton BB param"};
  Configurable<std::string> ConfPIDBBAntiProton{"ConfPIDBBAntiProton", "Analysis/PWGHF/ML/HFTrigger/TPC/AntiProton", "Path to the CCDB ocject for anti-proton BB param"};
  Configurable<std::string> ConfPIDBBKaon{"ConfPIDBBKaon", "Analysis/PWGHF/ML/HFTrigger/TPC/Kaon", "Path to the CCDB ocject for kaon BB param"};
  Configurable<std::string> ConfPIDBBAntiKaon{"ConfPIDBBAntiKaon", "Analysis/PWGHF/ML/HFTrigger/TPC/AntiKaon", "Path to the CCDB ocject for anti-kaon BB param"};
  Configurable<std::string> ConfPIDBBPion{"ConfPIDBBPion", "Analysis/PWGHF/ML/HFTrigger/TPC/Pion", "Path to the CCDB ocject for pion BB param"};
  Configurable<std::string> ConfPIDBBAntiPion{"ConfPIDBBAntiPion", "Analysis/PWGHF/ML/HFTrigger/TPC/AntiPion", "Path to the CCDB ocject for anti-pion BB param"};
  Configurable<bool> ConfRejectNotPropagatedTracks{"ConfRejectNotPropagatedTracks", false, "True: reject not propagated tracks"};
  Configurable<float> ConfPIDCutsTPCF1Proton{"ConfPIDCutsTPCF1Proton", 2, "Particle PID selections using TPC"};
  Configurable<float> ConfPIDCutsTOFF1Proton{"ConfPIDCutsTOFF1Proton", 2, "Particle PID selections using TOF"};
  Configurable<int> strategyPIDPion{"strategyPIDPion", 0, "PID strategy Pion"};
  Configurable<int> strategyPIDKaon{"strategyPIDKaon", 0, "PID strategy Kaon"};
  Configurable<int> strategyPIDProton{"strategyPIDProton", 1, "PID strategy Proton"};
  Configurable<double> pionMomentumPID{"pionMomentumPID", 0.5, "pi momentum range for TPC PID selection"};
  Configurable<double> kaonMomentumPID{"kaonMomentumPID", 0.45, "ka momentum range for TPC PID selection"};
  Configurable<double> protonMomentumPID{"protonMomentumPID", 0.75, "pr momentum range for TPC PID selection"};

  // Configs for track cut
  Configurable<float> ConfPtCutsF1Proton{"ConfPtCutsF1Proton", 0.1, "Particle Momentum selections"};
  Configurable<float> ConfTrkEtaF1Proton{"ConfTrkEtaF1Proton", 0.85, "Eta"};
  Configurable<float> ConfTPCNClustersMinF1Proton{"ConfTPCNClustersMinF1Proton", 80, " Minimum number of TPC cluster"};
  Configurable<float> ConfTrkTPCcRowsMinF1Proton{"ConfTrkTPCcRowsMinF1Proton", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> ConfTrkTPCfClsF1Proton{"ConfTrkTPCfClsF1Proton", 0.83, "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> ConfTrkTPCsClsMaxF1Proton{"ConfTrkTPCsClsMaxF1Proton", 160, "Maximum number of shared TPC clusters"};
  Configurable<float> ConfTrkITSnclsMinF1Proton{"ConfTrkITSnclsMinF1Proton", 0, "Minimum number of ITS clusters"};
  Configurable<float> ConfTrkITSnclsIbMinF1Proton{"ConfTrkITSnclsIbMinF1Proton", 0, "Minimum number of ITS clusters in the inner barrel"};
  Configurable<float> ConfTrkDCAxyMaxF1Proton{"ConfTrkDCAxyMaxF1Proton", 0.15, "Maximum DCA_xy"};
  Configurable<float> ConfTrkDCAzMaxF1Proton{"ConfTrkDCAzMaxF1Proton", 0.3, "Maximum DCA_z"};

  // Checks taken from global track definition
  Configurable<bool> ConfTrkRequireChi2MaxTPC{"ConfTrkRequireChi2MaxTPC", false, "True: require max chi2 per TPC cluster"};
  Configurable<bool> ConfTrkRequireChi2MaxITS{"ConfTrkRequireChi2MaxITS", false, "True: require max chi2 per ITS cluster"};
  Configurable<float> ConfTrkMaxChi2PerClusterTPC{"ConfTrkMaxChi2PerClusterTPC", 4.0f, "Minimal track selection: max allowed chi2 per TPC cluster"};  // 4.0 is default
  Configurable<float> ConfTrkMaxChi2PerClusterITS{"ConfTrkMaxChi2PerClusterITS", 36.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default
  Configurable<bool> ConfTrkTPCRefit{"ConfTrkTPCRefit", false, "True: require TPC refit"};
  Configurable<bool> ConfTrkITSRefit{"ConfTrkITSRefit", false, "True: require ITS refit"};

  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0DCADaughMax{"ConfV0DCADaughMax", 1.8f, "Maximum DCA between the V0 daughters"};
  Configurable<float> ConfV0CPAMin{"ConfV0CPAMin", 0.985f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 0.2f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<float> ConfV0DecVtxMax{"ConfV0DecVtxMax", 100.f, "Maximum distance from primary vertex"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "Minimum V0 CosPA to PV"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 40, "Maximum V0 life time"};
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 2, "Sigma cut on KS0 mass"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.85f, "V0 Daugh sel: max eta"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 60.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughDCAMin{"ConfDaughDCAMin", 0.04f, "V0 Daugh sel:  Max. DCA Daugh to PV (cm)"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for KS0 daughters"};

  // Configs for F1 candidate
  Configurable<double> cMaxMassKKs0{"cMaxMassKKs0", 1.04, "Mass cut on K-KS0 pair"};
  Configurable<double> cMaxMassF1{"cMaxMassF1", 1.80001, "Mass cut on F1 resonance"};
  Configurable<double> cMinF1Pt{"cMinF1Pt", 1.0, "Minimum pT cut on F1"};
  Configurable<double> cMinKaonPt{"cMinKaonPt", 0.3, "Minimum pT cut on Kaon daughter"};
  Configurable<double> cMaxProtonPt{"cMaxProtonPt", 2.0, "Maximum pT cut on Proton"};

  // config Femto relative momentum
  Configurable<double> cMaxRelMom{"cMaxRelMom", 0.5, "Relative momentum cut"};

  // Histogram
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", ";; Number of events", 2, 0.0f, 2.0f)};
  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hEventstat", "hEventstat", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
                                             {"hInvMassf1", "hInvMassf1", {HistType::kTH2F, {{400, 1.1f, 1.9f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassf1Like", "hInvMassf1Like", {HistType::kTH2F, {{400, 1.1f, 1.9f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassf1kstar", "hInvMassf1kstar", {HistType::kTH3F, {{400, 1.1f, 1.9f}, {100, 0.0f, 10.0f}, {8, 0.0f, 0.8f}}}},
                                             {"hkstarDist", "hkstarDist", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
                                             {"hDCAxy", "hDCAxy", {HistType::kTH1F, {{100, -5.0f, 5.0f}}}},
                                             {"hDCAz", "hDCAz", {HistType::kTH1F, {{100, -5.0f, 5.0f}}}},
                                             {"hPhi", "hPhi", {HistType::kTH1F, {{1400, -7.0f, 7.0f}}}},
                                             {"hPhiSphero", "hPhiSphero", {HistType::kTH1F, {{1400, -7.0f, 7.0f}}}},
                                             {"hEta", "hEta", {HistType::kTH1F, {{20, -1.0f, 1.0f}}}},
                                             {"hNsigmaPtpionTPC", "hNsigmaPtpionTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtpionTOF", "hNsigmaPtpionTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH3F, {{200, -10.0f, 10.0f}, {200, -20.0f, 20.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtprotonTPC", "hNsigmaPtprotonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtprotonTOF", "hNsigmaPtprotonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassk0", "hInvMassk0", {HistType::kTH2F, {{200, 0.4f, 0.6f}, {100, 0.0f, 10.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(url.value);
    ccdbApi.init(url);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    zorroSummary.setObject(zorro.getZorroSummary());
    hProcessedEvents->GetXaxis()->SetBinLabel(1, "All Trigger events");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events with F1 Trigger");
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (ConfEvtSelectZvtx && std::abs(col.posZ()) > ConfEvtZvtx) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedTrack(T const& track)
  {
    const auto pT = track.pt();
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto tpcRClsC = track.tpcCrossedRowsOverFindableCls();
    const auto tpcNClsC = track.tpcNClsCrossedRows();
    const auto tpcNClsS = track.tpcNClsShared();
    const auto itsNCls = track.itsNCls();
    const auto itsNClsIB = track.itsNClsInnerBarrel();
    const auto dcaXY = track.dcaXY();
    const auto dcaZ = track.dcaZ();
    if (pT < ConfPtCutsF1Proton) {
      return false;
    }
    if (std::abs(eta) > ConfTrkEtaF1Proton) {
      return false;
    }
    if (tpcNClsF < ConfTPCNClustersMinF1Proton) {
      return false;
    }
    if (tpcRClsC < ConfTrkTPCfClsF1Proton) {
      return false;
    }
    if (tpcNClsC < ConfTrkTPCcRowsMinF1Proton) {
      return false;
    }
    if (tpcNClsS > ConfTrkTPCsClsMaxF1Proton) {
      return false;
    }
    if (itsNCls < ConfTrkITSnclsMinF1Proton) {
      return false;
    }
    if (itsNClsIB < ConfTrkITSnclsIbMinF1Proton) {
      return false;
    }
    if (std::abs(dcaXY) > ConfTrkDCAxyMaxF1Proton) {
      return false;
    }
    if (std::abs(dcaZ) > ConfTrkDCAzMaxF1Proton) {
      return false;
    }
    // TODO: which dca, put dcaxy for now
    if (ConfRejectNotPropagatedTracks && std::abs(dcaXY) > 1e3) {
      return false;
    }
    if (ConfTrkRequireChi2MaxTPC && track.tpcChi2NCl() >= ConfTrkMaxChi2PerClusterTPC) {
      return false;
    }
    if (ConfTrkRequireChi2MaxITS && track.itsChi2NCl() >= ConfTrkMaxChi2PerClusterITS) {
      return false;
    }
    if (ConfTrkTPCRefit && !track.hasTPC()) {
      return false;
    }
    if (ConfTrkITSRefit && !track.hasITS()) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionGlobalTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor())) {
      return false;
    }
    return true;
  }

  template <typename T>
  double updatePID(T const& track, double bgScaling, std::vector<double> BB)
  {
    double expBethe = tpc::BetheBlochAleph(static_cast<double>(track.tpcInnerParam() * bgScaling), BB[0], BB[1], BB[2], BB[3], BB[4]);
    double expSigma = expBethe * BB[5];
    return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track, float charge, double nsigmaV0Daughter)
  {
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto dcaXY = track.dcaXY();
    const auto sign = track.sign();

    if (charge < 0 && sign > 0) {
      return false;
    }
    if (charge > 0 && sign < 0) {
      return false;
    }
    if (std::abs(eta) > ConfDaughEta) {
      return false;
    }
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (std::abs(dcaXY) < ConfDaughDCAMin) {
      return false;
    }

    if (std::abs(nsigmaV0Daughter) > ConfDaughPIDCuts) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool SelectionPID(const T& candidate, int PIDstrategy, int particle, double updatensigma)
  {
    if (PIDstrategy == 1) {
      if (particle == 0) {
        if (std::abs(candidate.p()) < pionMomentumPID && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton) {
          return true;
        } else if (std::abs(candidate.p()) >= pionMomentumPID && candidate.hasTOF() && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaPi()) < ConfPIDCutsTOFF1Proton) {
          return true;
        }
      } else if (particle == 1) {
        if (std::abs(candidate.p()) < kaonMomentumPID && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton) {
          return true;
        } else if (std::abs(candidate.p()) >= kaonMomentumPID && candidate.hasTOF() && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaKa()) < ConfPIDCutsTOFF1Proton) {
          return true;
        }
      } else if (particle == 2) {
        if (std::abs(candidate.p()) < protonMomentumPID && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton) {
          return true;
        } else if (std::abs(candidate.p()) >= protonMomentumPID && candidate.hasTOF() && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaPr()) < ConfPIDCutsTOFF1Proton) {
          return true;
        }
      }
    } else if (PIDstrategy == 0) {
      if (candidate.hasTOF()) {
        if (particle == 0 && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaPi()) < ConfPIDCutsTOFF1Proton) {
          return true;
        } else if (particle == 1 && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaKa()) < ConfPIDCutsTOFF1Proton) {
          return true;
        } else if (particle == 2 && std::abs(updatensigma) < ConfPIDCutsTPCF1Proton && std::abs(candidate.tofNSigmaPr()) < ConfPIDCutsTOFF1Proton) {
          return true;
        }
      } else if (std::abs(updatensigma) < ConfPIDCutsTPCF1Proton) {
        return true;
      }
    }
    return false;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (fabs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    const float pT = candidate.pt();
    const std::vector<float> decVtx = {candidate.x(), candidate.y(), candidate.z()};
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = candidate.dcaV0daughters();
    const float cpav0 = candidate.v0cosPA();

    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header
    float lowmasscutks0 = 0.497 - 2.0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + 2.0 * cSigmaMassKs0;

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    for (size_t i = 0; i < decVtx.size(); i++) {
      if (decVtx.at(i) > ConfV0DecVtxMax) {
        return false;
      }
    }
    if (fabs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    return true;
  }

  float getkstar(const ROOT::Math::PtEtaPhiMVector part1,
                 const ROOT::Math::PtEtaPhiMVector part2)
  {
    const ROOT::Math::PtEtaPhiMVector trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax =
      beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay =
      beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    ROOT::Math::PxPyPzMVector PartOneCMS(part1);
    ROOT::Math::PxPyPzMVector PartTwoCMS(part2);
    const ROOT::Math::Boost boostPRF =
      ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    const ROOT::Math::PxPyPzMVector trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  std::vector<double> setValuesBB(o2::ccdb::CcdbApi& ccdbApi, aod::BCsWithTimestamps::iterator const& bunchCrossing, const std::string ccdbPath)
  {
    std::map<std::string, std::string> metadata;
    auto h = ccdbApi.retrieveFromTFileAny<TH1F>(ccdbPath, metadata, bunchCrossing.timestamp());
    // auto h = ccdb->getForTimeStamp<TH1F>(ccdbPath, bunchCrossing.timestamp()); // check if possible to use this without getting fatal
    if (!h) {
      std::vector<double> dummy;
      LOG(info) << "File from CCDB in path " << ccdbPath << " was not found for run " << bunchCrossing.runNumber() << "and timestamp" << bunchCrossing.timestamp() << ". Will use default PID task values!";
      return dummy;
    }
    LOG(info) << "File from CCDB in path " << ccdbPath << " was found for run " << bunchCrossing.runNumber() << "!";

    TAxis* axis = h->GetXaxis();
    std::vector<double> v{static_cast<double>(h->GetBinContent(axis->FindBin("bb1"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb2"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb3"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb4"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("bb5"))),
                          static_cast<double>(h->GetBinContent(axis->FindBin("Resolution")))};
    return v;
  }

  /// Compute the spherocity of an event
  /// Important here is that the filter on tracks does not interfere here!
  /// In Run 2 we used here global tracks within |eta| < 0.8
  /// \tparam T type of the tracks
  /// \param tracks All tracks
  /// \return value of the spherocity of the event
  template <typename T>
  float ComputeSpherocity(T const& tracks, int nTracksMin, int spdef)
  {
    // if number of tracks is not enough for spherocity estimation.
    int ntrks = tracks.size();
    if (ntrks < nTracksMin)
      return -99.;

    // start computing spherocity

    float ptSum = 0.;
    for (auto const& track : tracks) {

      qaRegistry.fill(HIST("hPhiSphero"), track.phi());

      if (spdef == 0) {
        ptSum += 1.;
      } else {
        ptSum += track.pt();
      }
    }

    float tempSph = 1.;
    for (int i = 0; i < 360 / 0.1; ++i) {
      float sum = 0., pt = 0.;
      float phiparm = (TMath::Pi() * i * 0.1) / 180.;
      float nx = TMath::Cos(phiparm);
      float ny = TMath::Sin(phiparm);
      for (auto const& trk : tracks) {
        pt = trk.pt();
        if (spdef == 0) {
          pt = 1.;
        }
        float phi = trk.phi();
        float px = pt * TMath::Cos(phi);
        float py = pt * TMath::Sin(phi);
        // sum += pt * abs(sin(phiparm - phi));
        sum += TMath::Abs(px * ny - py * nx);
      }
      float sph = TMath::Power((sum / ptSum), 2);
      if (sph < tempSph)
        tempSph = sph;
    }

    return TMath::Power(TMath::Pi() / 2., 2) * tempSph;
  }

  std::vector<double> BBProton, BBAntiproton, BBPion, BBAntipion, BBKaon, BBAntikaon;
  ROOT::Math::PtEtaPhiMVector F1Vector, F1VectorDummy, F1d1dummy, F1d2dummy, F1d3dummy, KKs0Vector, ProtonVectorDummy, ProtonVectorDummy2;
  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass();   // FIXME: Get from the common header
  double massKa = TDatabasePDG::Instance()->GetParticle(kKPlus)->Mass();    // FIXME: Get from the common header
  double massPr = TDatabasePDG::Instance()->GetParticle(kProton)->Mass();   // FIXME: Get from the common header
  double massK0s = TDatabasePDG::Instance()->GetParticle(kK0Short)->Mass(); // FIXME: Get from the common header

  double massF1{0.};
  double masskKs0{0.};
  double pT{0.};
  int currentRunNumber = -999;
  int lastRunNumber = -999;
  double betaX = 0;
  double betaY = 0;
  double betaZ = 0;
  double relativeMomentum = 999;
  // Pre-filters for primary track
  Filter acceptanceFilter = nabs(aod::track::eta) < ConfTrkEtaF1Proton && nabs(aod::track::pt) > ConfPtCutsF1Proton;
  Filter dcaFilter = nabs(aod::track::dcaXY) < ConfTrkDCAxyMaxF1Proton && nabs(aod::track::dcaZ) < ConfTrkDCAzMaxF1Proton;

  // using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using EventCandidates = aod::Collisions;
  using ResoV0s = aod::V0Datas;
  using PrimaryTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                         aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                         aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                         aod::pidTPCFullPr, aod::pidTOFFullPr>>;

  void processF1ProtonReducedTable(EventCandidates::iterator const& collision, aod::BCsWithTimestamps const&, PrimaryTrackCandidates const& tracks, ResoV0s const& V0s)
  {
    o2::aod::ITSResponse itsResponse;
    bool keepEventF1Proton = false;
    int numberF1 = 0;
    // keep track of indices
    std::vector<int64_t> PionIndex = {};
    std::vector<int64_t> KaonIndex = {};
    std::vector<int64_t> ProtonIndex = {};
    std::vector<int64_t> KshortPosDaughIndex = {};
    std::vector<int64_t> KshortNegDaughIndex = {};

    ////// track indices to check share track in table
    std::vector<int64_t> F1PionIndex = {};
    std::vector<int64_t> F1KaonIndex = {};
    std::vector<int64_t> F1ProtonIndex = {};
    std::vector<int64_t> F1KshortDaughterPositiveIndex = {};
    std::vector<int64_t> F1KshortDaughterNegativeIndex = {};

    // keep charge of track
    std::vector<float> PionCharge = {};
    std::vector<float> KaonCharge = {};
    std::vector<float> ProtonCharge = {};
    // std::vector<float> ProtonChargeFinal = {};

    // keep TPC PID of proton
    std::vector<float> ProtonTPCNsigma = {};
    // std::vector<float> ProtonTPCNsigmaFinal = {};

    // keep TOF PID of proton
    std::vector<float> ProtonTOFNsigma = {};
    // std::vector<float> ProtonTOFNsigmaFinal = {};

    // keep TOF Hit of proton
    std::vector<int> ProtonTOFHit = {};
    // std::vector<int> ProtonTOFHitFinal = {};

    // keep TOF Hit of pion
    std::vector<int> PionTOFHit = {};
    std::vector<int> PionTOFHitFinal = {};
    std::vector<float> PionTPC = {};
    std::vector<float> PionTPCFinal = {};

    // keep TOF Hit of kaon
    std::vector<int> KaonTOFHit = {};
    std::vector<int> KaonTOFHitFinal = {};
    std::vector<float> KaonTPC = {};
    std::vector<float> KaonTPCFinal = {};
    std::vector<float> KaonTPCPionHypo = {};
    std::vector<float> KaonTPCPionHypoFinal = {};

    // keep kaon-kshort mass of f1resonance
    std::vector<float> f1kaonkshortmass = {};

    // keep status of F1 signal unlike or wrong sign
    std::vector<float> f1signal = {};

    // Prepare vectors for different species
    std::vector<ROOT::Math::PtEtaPhiMVector> protons, kaons, pions, kshorts, f1resonance, f1resonanced1, f1resonanced2, f1resonanced3;
    // , protonsfinal;
    float kstar = 999.f;

    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    hProcessedEvents->Fill(0.5);
    bool zorroSelected = false;
    if (isSelectedEvent(collision)) {
      if (currentRunNumber != lastRunNumber) {
        if (cfgSkimmedProcessing) {
          zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), "fTriggerEventF1Proton");
          zorro.populateHistRegistry(qaRegistry, bc.runNumber());
        }
        if (ConfUseManualPIDproton || ConfUseManualPIDkaon || ConfUseManualPIDpion) {
          if (ConfUseManualPIDproton) {
            BBProton = setValuesBB(ccdbApi, bc, ConfPIDBBProton);
            BBAntiproton = setValuesBB(ccdbApi, bc, ConfPIDBBAntiProton);
          }
          if (ConfUseManualPIDpion) {
            BBPion = setValuesBB(ccdbApi, bc, ConfPIDBBPion);
            BBAntipion = setValuesBB(ccdbApi, bc, ConfPIDBBAntiPion);
          }
          if (ConfUseManualPIDkaon) {
            BBKaon = setValuesBB(ccdbApi, bc, ConfPIDBBKaon);
            BBAntikaon = setValuesBB(ccdbApi, bc, ConfPIDBBAntiKaon);
          }
        }
        lastRunNumber = currentRunNumber;
      }
      if (cfgSkimmedProcessing) {
        zorroSelected = zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC());
      } else {
        zorroSelected = true;
      }
      if (zorroSelected) {
        hProcessedEvents->Fill(1.5);
        for (auto& track : tracks) {
          if (!isSelectedTrack(track)) {
            continue;
          }
          if (!selectionGlobalTrack(track)) {
            continue;
          }
          qaRegistry.fill(HIST("hDCAxy"), track.dcaXY());
          qaRegistry.fill(HIST("hDCAz"), track.dcaZ());
          qaRegistry.fill(HIST("hEta"), track.eta());
          qaRegistry.fill(HIST("hPhi"), track.phi());
          double nTPCSigmaP[3]{track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
          double nTPCSigmaN[3]{track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr()};
          if (ConfUseManualPIDproton) {
            auto bgScalingProton = 1 / massPr; // momentum scaling?
            if (BBProton.size() == 6)
              nTPCSigmaP[2] = updatePID(track, bgScalingProton, BBProton);
            if (BBAntiproton.size() == 6)
              nTPCSigmaN[2] = updatePID(track, bgScalingProton, BBAntiproton);
          }
          if (ConfUseManualPIDkaon) {
            auto bgScalingKaon = 1 / massKa; // momentum scaling?
            if (BBKaon.size() == 6)
              nTPCSigmaP[1] = updatePID(track, bgScalingKaon, BBKaon);
            if (BBAntikaon.size() == 6)
              nTPCSigmaN[1] = updatePID(track, bgScalingKaon, BBAntikaon);
          }
          if (ConfUseManualPIDpion) {
            auto bgScalingPion = 1 / massPi; // momentum scaling?
            if (BBPion.size() == 6)
              nTPCSigmaP[0] = updatePID(track, bgScalingPion, BBPion);
            if (BBAntipion.size() == 6)
              nTPCSigmaN[0] = updatePID(track, bgScalingPion, BBAntipion);
          }

          if ((track.sign() > 0 && SelectionPID(track, strategyPIDPion, 0, nTPCSigmaP[0]) && (itsResponse.nSigmaITS<o2::track::PID::Pion>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Pion>(track) < 3.0)) || (track.sign() < 0 && SelectionPID(track, strategyPIDPion, 0, nTPCSigmaN[0]) && (itsResponse.nSigmaITS<o2::track::PID::Pion>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Pion>(track) < 3.0))) {
            ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), massPi);
            pions.push_back(temp);
            PionIndex.push_back(track.globalIndex());
            PionCharge.push_back(track.sign());
            auto PionTOF = 0;
            if (track.sign() > 0) {
              qaRegistry.fill(HIST("hNsigmaPtpionTPC"), nTPCSigmaP[0], track.pt());
              PionTPC.push_back(nTPCSigmaP[0]);
            }
            if (track.sign() < 0) {
              qaRegistry.fill(HIST("hNsigmaPtpionTPC"), nTPCSigmaN[0], track.pt());
              PionTPC.push_back(nTPCSigmaN[0]);
            }
            if (track.hasTOF()) {
              qaRegistry.fill(HIST("hNsigmaPtpionTOF"), track.tofNSigmaPi(), track.pt());
              PionTOF = 1;
            }
            PionTOFHit.push_back(PionTOF);
          }

          if ((track.pt() > cMinKaonPt && track.sign() > 0 && SelectionPID(track, strategyPIDKaon, 1, nTPCSigmaP[1]) && (itsResponse.nSigmaITS<o2::track::PID::Kaon>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track) < 3.0)) || (track.pt() > cMinKaonPt && track.sign() < 0 && SelectionPID(track, strategyPIDKaon, 1, nTPCSigmaN[1]) && (itsResponse.nSigmaITS<o2::track::PID::Kaon>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Kaon>(track) < 3.0))) {
            ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), massKa);
            kaons.push_back(temp);
            KaonIndex.push_back(track.globalIndex());
            KaonCharge.push_back(track.sign());
            auto KaonTOF = 0;
            if (track.sign() > 0) {
              qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), nTPCSigmaP[1], nTPCSigmaP[0], track.pt());
              KaonTPC.push_back(nTPCSigmaP[1]);
              KaonTPCPionHypo.push_back(nTPCSigmaP[0]);
            }
            if (track.sign() < 0) {
              qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), nTPCSigmaN[1], nTPCSigmaN[0], track.pt());
              KaonTPC.push_back(nTPCSigmaN[1]);
              KaonTPCPionHypo.push_back(nTPCSigmaN[0]);
            }
            if (track.hasTOF()) {
              qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track.tofNSigmaKa(), track.pt());
              KaonTOF = 1;
            }
            KaonTOFHit.push_back(KaonTOF);
          }

          if ((track.pt() < cMaxProtonPt && track.sign() > 0 && SelectionPID(track, strategyPIDProton, 2, nTPCSigmaP[2]) && (itsResponse.nSigmaITS<o2::track::PID::Proton>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Proton>(track) < 3.0)) || (track.pt() < cMaxProtonPt && track.sign() < 0 && SelectionPID(track, strategyPIDProton, 2, nTPCSigmaN[2]) && (itsResponse.nSigmaITS<o2::track::PID::Proton>(track) > -3.0 && itsResponse.nSigmaITS<o2::track::PID::Proton>(track) < 3.0))) {
            ROOT::Math::PtEtaPhiMVector temp(track.pt(), track.eta(), track.phi(), massPr);
            protons.push_back(temp);
            ProtonIndex.push_back(track.globalIndex());
            ProtonCharge.push_back(track.sign());
            if (track.sign() > 0) {
              qaRegistry.fill(HIST("hNsigmaPtprotonTPC"), nTPCSigmaP[2], track.pt());
              ProtonTPCNsigma.push_back(nTPCSigmaP[2]);
            }
            if (track.sign() < 0) {
              qaRegistry.fill(HIST("hNsigmaPtprotonTPC"), nTPCSigmaN[2], track.pt());
              ProtonTPCNsigma.push_back(nTPCSigmaN[2]);
            }
            if (track.hasTOF()) {
              qaRegistry.fill(HIST("hNsigmaPtprotonTOF"), track.tofNSigmaPr(), track.pt());
              ProtonTOFNsigma.push_back(track.tofNSigmaPr());
              ProtonTOFHit.push_back(1);
            }
            if (!track.hasTOF()) {
              ProtonTOFNsigma.push_back(999.0);
              ProtonTOFHit.push_back(0);
            }
          }
        } // track loop end
        for (auto& v0 : V0s) {

          if (!SelectionV0(collision, v0)) {
            continue;
          }
          auto postrack = v0.template posTrack_as<PrimaryTrackCandidates>();
          auto negtrack = v0.template negTrack_as<PrimaryTrackCandidates>();
          double nTPCSigmaPos[1]{postrack.tpcNSigmaPi()};
          double nTPCSigmaNeg[1]{negtrack.tpcNSigmaPi()};
          if (ConfUseManualPIDdaughterPion) {
            auto bgScalingPion = 1 / massPi; // momentum scaling?
            if (BBPion.size() == 6)
              nTPCSigmaPos[0] = updatePID(postrack, bgScalingPion, BBPion);
            if (BBAntipion.size() == 6)
              nTPCSigmaNeg[0] = updatePID(negtrack, bgScalingPion, BBAntipion);
          }
          if (!isSelectedV0Daughter(postrack, 1, nTPCSigmaPos[0])) {
            continue;
          }
          if (!isSelectedV0Daughter(negtrack, -1, nTPCSigmaNeg[0])) {
            continue;
          }
          qaRegistry.fill(HIST("hInvMassk0"), v0.mK0Short(), v0.pt());
          ROOT::Math::PtEtaPhiMVector temp(v0.pt(), v0.eta(), v0.phi(), massK0s);
          kshorts.push_back(temp);
          KshortPosDaughIndex.push_back(postrack.globalIndex());
          KshortNegDaughIndex.push_back(negtrack.globalIndex());
        }

        if (pions.size() != 0 && kaons.size() != 0 && kshorts.size() != 0) {
          for (auto ipion = pions.begin(); ipion != pions.end(); ++ipion) {
            for (auto ikaon = kaons.begin(); ikaon != kaons.end(); ++ikaon) {
              auto i1 = std::distance(pions.begin(), ipion);
              auto i2 = std::distance(kaons.begin(), ikaon);
              if (PionIndex.at(i1) == KaonIndex.at(i2))
                continue;
              for (auto ikshort = kshorts.begin(); ikshort != kshorts.end(); ++ikshort) {
                auto i3 = std::distance(kshorts.begin(), ikshort);
                if (PionIndex.at(i1) == KshortPosDaughIndex.at(i3))
                  continue;
                if (PionIndex.at(i1) == KshortNegDaughIndex.at(i3))
                  continue;
                KKs0Vector = kaons.at(i2) + kshorts.at(i3);
                if (KKs0Vector.M() > cMaxMassKKs0)
                  continue;
                F1Vector = KKs0Vector + pions.at(i1);
                if (F1Vector.M() > cMaxMassF1)
                  continue;
                if (F1Vector.Pt() < cMinF1Pt)
                  continue;

                // check if the pair is unlike or wrongsign
                auto pairsign = 100;
                if (PionCharge.at(i1) * KaonCharge.at(i2) > 0) {
                  qaRegistry.fill(HIST("hInvMassf1Like"), F1Vector.M(), F1Vector.Pt());
                  pairsign = -1;
                } else if (PionCharge.at(i1) * KaonCharge.at(i2) < 0) {
                  if (KaonCharge.at(i2) > 0) {
                    pairsign = 1;
                  }
                  if (KaonCharge.at(i2) < 0) {
                    pairsign = 2;
                  }
                }
                ROOT::Math::PtEtaPhiMVector temp(F1Vector.Pt(), F1Vector.Eta(), F1Vector.Phi(), F1Vector.M());
                f1resonance.push_back(temp);
                f1resonanced1.push_back(pions.at(i1));
                f1resonanced2.push_back(kaons.at(i2));
                f1resonanced3.push_back(kshorts.at(i3));
                f1signal.push_back(pairsign);
                f1kaonkshortmass.push_back(KKs0Vector.M());
                F1PionIndex.push_back(PionIndex.at(i1));
                F1KaonIndex.push_back(KaonIndex.at(i2));
                F1KshortDaughterPositiveIndex.push_back(KshortPosDaughIndex.at(i3));
                F1KshortDaughterNegativeIndex.push_back(KshortNegDaughIndex.at(i3));
                PionTOFHitFinal.push_back(PionTOFHit.at(i1));           // Pion TOF Hit
                KaonTOFHitFinal.push_back(KaonTOFHit.at(i2));           // Kaon TOF Hit
                PionTPCFinal.push_back(PionTPC.at(i1));                 // Pion TPC
                KaonTPCFinal.push_back(KaonTPC.at(i2));                 // Kaon TPC
                KaonTPCPionHypoFinal.push_back(KaonTPCPionHypo.at(i2)); // Kaon TPC
                if (pairsign > 0) {
                  qaRegistry.fill(HIST("hInvMassf1"), F1Vector.M(), F1Vector.Pt());
                  numberF1 = numberF1 + 1;
                  for (auto iproton = protons.begin(); iproton != protons.end(); ++iproton) {
                    auto i4 = std::distance(protons.begin(), iproton);
                    if (ProtonIndex.at(i4) == PionIndex.at(i1))
                      continue;
                    if (ProtonIndex.at(i4) == KaonIndex.at(i2))
                      continue;
                    if (ProtonIndex.at(i4) == KshortPosDaughIndex.at(i3))
                      continue;
                    if (ProtonIndex.at(i4) == KshortNegDaughIndex.at(i3))
                      continue;
                    kstar = getkstar(F1Vector, *iproton);
                    qaRegistry.fill(HIST("hkstarDist"), kstar);
                    if (kstar > cMaxRelMom) {
                      continue;
                    }
                    qaRegistry.fill(HIST("hInvMassf1kstar"), F1Vector.M(), F1Vector.Pt(), kstar);
                    keepEventF1Proton = true;
                    /*ProtonVectorDummy = protons.at(i4);
                    if (numberF1 == 1 && keepEventF1Proton) {
                                  //////////// Fill final proton information after pairing//////////
                                  ROOT::Math::PtEtaPhiMVector temp(ProtonVectorDummy.Pt(), ProtonVectorDummy.Eta(), ProtonVectorDummy.Phi(), massPr);
                                  protonsfinal.push_back(temp);                           // 4 vector
                                  ProtonChargeFinal.push_back(ProtonCharge.at(i4));       // Charge
                                  ProtonTOFHitFinal.push_back(ProtonTOFHit.at(i4));       // TOF Hit
                                  ProtonTOFNsigmaFinal.push_back(ProtonTOFNsigma.at(i4)); // Nsigma TOF
                                  ProtonTPCNsigmaFinal.push_back(ProtonTPCNsigma.at(i4)); // Nsigma TPC
                                  F1ProtonIndex.push_back(ProtonIndex.at(i4));            // proton index for share track
                      }*/
                  }
                } // pair sign
              }
            }
          }
        }
      }
    }
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    if (numberF1 > 0 && (f1resonance.size() == f1signal.size()) && (f1resonance.size() == f1kaonkshortmass.size()) && (f1resonance.size() == f1resonanced1.size()) && (f1resonance.size() == f1resonanced2.size()) && (f1resonance.size() == f1resonanced3.size())) {
      qaRegistry.fill(HIST("hEventstat"), 1.5);
      if (keepEventF1Proton) {
        qaRegistry.fill(HIST("hEventstat"), 2.5);
        auto eventspherocity = 0.0;
        // ComputeSpherocity(tracks, trackSphMin, trackSphDef);
        /////////// Fill collision table///////////////
        redf1pevents(bc.globalBC(), currentRunNumber, bc.timestamp(), collision.posZ(), collision.numContrib(), eventspherocity);
        auto indexEvent = redf1pevents.lastIndex();
        //// Fill track table for F1//////////////////
        for (auto if1 = f1resonance.begin(); if1 != f1resonance.end(); ++if1) {
          auto i5 = std::distance(f1resonance.begin(), if1);
          F1VectorDummy = f1resonance.at(i5);
          F1d1dummy = f1resonanced1.at(i5);
          F1d2dummy = f1resonanced2.at(i5);
          F1d3dummy = f1resonanced3.at(i5);
          f1track(indexEvent, f1signal.at(i5), F1VectorDummy.Px(), F1VectorDummy.Py(), F1VectorDummy.Pz(), F1d1dummy.Px(), F1d1dummy.Py(), F1d1dummy.Pz(), F1d2dummy.Px(), F1d2dummy.Py(), F1d2dummy.Pz(), F1d3dummy.Px(), F1d3dummy.Py(), F1d3dummy.Pz(), PionTOFHitFinal.at(i5), KaonTOFHitFinal.at(i5), PionTPCFinal.at(i5), KaonTPCFinal.at(i5), KaonTPCPionHypoFinal.at(i5), F1VectorDummy.M(), f1kaonkshortmass.at(i5), F1PionIndex.at(i5), F1KaonIndex.at(i5), F1KshortDaughterPositiveIndex.at(i5), F1KshortDaughterNegativeIndex.at(i5));
        }
        //// Fill track table for proton//////////////////
        for (auto iproton = protons.begin(); iproton != protons.end(); ++iproton) {
          auto i6 = std::distance(protons.begin(), iproton);
          ProtonVectorDummy2 = protons.at(i6);
          protontrack(indexEvent, ProtonCharge.at(i6), ProtonVectorDummy2.Px(), ProtonVectorDummy2.Py(), ProtonVectorDummy2.Pz(), ProtonTPCNsigma.at(i6), ProtonTOFHit.at(i6), ProtonTOFNsigma.at(i6), ProtonIndex.at(i6));
        }
      }
    }
  }
  PROCESS_SWITCH(f1protonreducedtable, processF1ProtonReducedTable, "Process for create reduced table for f1-p analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<f1protonreducedtable>(cfg)};
}
