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

/// \file f1protonInitializer.cxx
/// check if the event have f1-p candidate
/// \author Sourav Kundu <sourav.kundu@cern.ch>
#include "PWGLF/DataModel/LFF1Tables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/collisionCuts.h"

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers
struct f1protoninitializer {

  Produces<aod::F1Collisions> f1Collisions;

  // Pre-selection cuts
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cMaxDCAXY{"cMaxDCAXY", 1.0, "Maximum DCA XY cut of primary track"};
  Configurable<float> cMaxDCAZ{"cMaxDCAZ", 1.0, "Maximum DCA Z cut of primary track"};
  Configurable<float> cMinTPCncls{"cMinTPCncls", 100.0, "Minimum number of TPC cluster"};
  Configurable<float> cMinTPCncr{"cMinTPCncr", 100.0, "Minimum number of TPC crossed rows"};
  Configurable<float> cMinTPCncrOverncls{"cMinTPCncrOverncls", 0.8, "Minimum number of TPC crossed rows over cluster"};
  Configurable<float> cMaxFracTPCncls{"cMaxFracTPCcls", 1.0, "Maximum fraction of shared TPC cluster"};

  /// V0 selcection
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 0.3, "Minimum V0 CosPA to PV"};
  Configurable<double> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};

  /// F1 mass cut///////////
  Configurable<double> cSigmaMassKs0{"cSigmaMassKs0", 2, "Sigma cut on KS0 mass"};
  Configurable<double> cMaxMassKKs0{"cMaxMassKKs0", 1.04, "Mass cut on K-KS0 pair"};
  Configurable<double> cMaxMassF1{"cMaxMassF1", 1.80001, "Mass cut on F1 resonance"};
  Configurable<double> cMaxRelMom{"cMaxRelMom", 0.5, "Relative momentum cut"};
  Configurable<double> cMinF1Pt{"cMinF1Pt", 1.0, "Minimum pT cut on F1"};

  /// PID///////////
  Configurable<double> pionMomentumPID{"pionMomentumPID", 0.5, "pi momentum range for TPC PID selection"};
  Configurable<double> kaonMomentumPID{"kaonMomentumPID", 0.45, "ka momentum range for TPC PID selection"};
  Configurable<double> protonMomentumPID{"protonMomentumPID", 0.75, "pr momentum range for TPC PID selection"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3, "nsigma cut TPC"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3, "nsigma cut TOF"};
  Configurable<int> strategyPIDPion{"strategyPIDPion", 0, "PID strategy Pion"};
  Configurable<int> strategyPIDKaon{"strategyPIDKaon", 0, "PID strategy Kaon"};
  Configurable<int> strategyPIDProton{"strategyPIDProton", 0, "PID strategy Proton"};

  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hEventstat", "hEventstat", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
                                             {"hInvMassf1", "hInvMassf1", {HistType::kTH2F, {{400, 1.1f, 1.9f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassf1Like", "hInvMassf1Like", {HistType::kTH2F, {{400, 1.1f, 1.9f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassKKs0", "hInvMassKKs0", {HistType::kTH1F, {{200, 0.9f, 1.1f}}}},
                                             {"hkstarDist", "hkstarDist", {HistType::kTH1F, {{300, 0.0f, 3.0f}}}},
                                             {"hDCAxy", "hDCAxy", {HistType::kTH1F, {{100, -5.0f, 5.0f}}}},
                                             {"hDCAz", "hDCAz", {HistType::kTH1F, {{100, -5.0f, 5.0f}}}},
                                             {"hPhi", "hPhi", {HistType::kTH1F, {{70, 0.0f, 7.0f}}}},
                                             {"hEta", "hEta", {HistType::kTH1F, {{20, -1.0f, 1.0f}}}},
                                             {"hNsigmaPtpionTPC", "hNsigmaPtpionTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtpionTOF", "hNsigmaPtpionTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTPC", "hNsigmaPtkaonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtkaonTOF", "hNsigmaPtkaonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtprotonTPC", "hNsigmaPtprotonTPC", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hNsigmaPtprotonTOF", "hNsigmaPtprotonTOF", {HistType::kTH2F, {{200, -10.0f, 10.0f}, {100, 0.0f, 10.0f}}}},
                                             {"hInvMassk0", "hInvMassk0", {HistType::kTH2F, {{200, 0.4f, 0.6f}, {100, 0.0f, 10.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  //////////////// Primary track selection and Ks0 selection//////////////////////////
  template <typename T>
  bool SelectionTrack(const T& candidate)
  {
    if (!candidate.isGlobalTrack()) {
      return false;
    }
    if (candidate.tpcNClsCrossedRows() < cMinTPCncr || candidate.tpcNClsCrossedRows() / candidate.tpcNClsFound() < cMinTPCncrOverncls || candidate.tpcNClsFound() < cMinTPCncls || candidate.tpcFractionSharedCls() > cMaxFracTPCncls) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool SelectionPi(const T& candidate)
  {
    if (std::abs(candidate.p()) < pionMomentumPID && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
      return true;
    } else if (std::abs(candidate.p()) >= pionMomentumPID && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool SelectionKa(const T& candidate)
  {
    if (std::abs(candidate.p()) < kaonMomentumPID && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
      return true;
    } else if (std::abs(candidate.p()) >= kaonMomentumPID && candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool SelectionPr(const T& candidate)
  {
    if (std::abs(candidate.p()) < protonMomentumPID && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
      return true;
    } else if (std::abs(candidate.p()) >= protonMomentumPID && candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool SelectionPID(const T& candidate, int PID)
  {
    if (candidate.hasTOF()) {
      if (PID == 0 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF) {
        return true;
      } else if (PID == 1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF) {
        return true;
      } else if (PID == 2 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC && std::abs(candidate.tofNSigmaPr()) < nsigmaCutTOF) {
        return true;
      }
    } else {
      if (PID == 0 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
        return true;
      } else if (PID == 1 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
        return true;
      } else if (PID == 2 && std::abs(candidate.tpcNSigmaPr()) < nsigmaCutTPC) {
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
    float CtauK0s = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;
    float lowmasscutks0 = 0.497 - 2.0 * cSigmaMassKs0;
    float highmasscutks0 = 0.497 + 2.0 * cSigmaMassKs0;
    if (fabs(CtauK0s) > cMaxV0LifeTime || candidate.mK0Short() < lowmasscutks0 || candidate.mK0Short() > highmasscutks0) {
      return false;
    }
    return true;
  }

  /////////////////////////////////////////////////////////////
  double massPi = o2::constants::physics::MassPiPlus;
  double massKa = o2::constants::physics::MassKPlus;
  double massK0s = o2::constants::physics::MassK0Short;
  double massF1{0.};
  double masskKs0{0.};
  double pT{0.};
  TLorentzVector F1Vector, ProtonVector, F1ProtonVector, RelativeMomentumVector;
  double betaX = 0;
  double betaY = 0;
  double betaZ = 0;
  double relativeMomentum = 999;
  // Pre-filters for primary track
  Filter acceptanceFilter = nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT;
  Filter dcaFilter = nabs(aod::track::dcaXY) < cMaxDCAXY && nabs(aod::track::dcaZ) < cMaxDCAZ;

  using EventCandidates = aod::Collisions;
  using ResoV0s = aod::V0Datas;
  using PrimaryTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                         aod::pidTPCFullPi, aod::pidTOFFullPi,
                                                         aod::pidTPCFullKa, aod::pidTOFFullKa,
                                                         aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  /*  Partition<PrimaryTrackCandidates> pionTracks = nabs(aod::pidtpc::tpcNSigmaPi) < 5.0f;
      Partition<PrimaryTrackCandidates> kaonTracks = nabs(aod::pidtpc::tpcNSigmaKa) < 5.0f;
      Partition<PrimaryTrackCandidates> protonTracks = nabs(aod::pidtpc::tpcNSigmaPr) < 5.0f;*/
  /*auto pionTracksGrouped = pionTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto kaonTracksGrouped = kaonTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());
    auto protonTracksGrouped = protonTracks->sliceByCached(aod::track::collisionId, collision.globalIndex());*/
  // for (auto& track1 : pionTracksGrouped) {
  void processF1Proton(EventCandidates::iterator const& collision, PrimaryTrackCandidates const& tracks, ResoV0s const& V0s)
  // void processF1Proton(aod::Collision const& collision, PrimaryTrackCandidates const& tracks, ResoV0s const& V0s)
  {
    qaRegistry.fill(HIST("hEventstat"), 0.5);
    int numberPion = 0;
    int numberPiKpair = 0;
    int numberF1 = 0;
    int numberF1Proton = 0;
    int numberF1ProtonFemto = 0;
    bool triggerF1 = false;
    bool triggerF1Proton = false;
    bool triggerF1ProtonFemto = false;
    for (auto track1 : tracks) {
      if (!SelectionTrack(track1)) {
        continue;
      }
      qaRegistry.fill(HIST("hDCAxy"), track1.dcaXY());
      qaRegistry.fill(HIST("hDCAz"), track1.dcaZ());
      qaRegistry.fill(HIST("hEta"), track1.eta());
      qaRegistry.fill(HIST("hPhi"), track1.phi());
      if (!SelectionPi(track1) && strategyPIDPion == 0) {
        continue;
      }
      if (!SelectionPID(track1, 0) && strategyPIDPion == 1) {
        continue;
      }
      numberPion = numberPion + 1;
      qaRegistry.fill(HIST("hNsigmaPtpionTPC"), track1.tpcNSigmaPi(), track1.pt());
      if (track1.hasTOF()) {
        qaRegistry.fill(HIST("hNsigmaPtpionTOF"), track1.tofNSigmaPi(), track1.pt());
      }
      auto track1ID = track1.globalIndex();
      for (auto track2 : tracks) {
        if (!SelectionTrack(track2)) {
          continue;
        }
        if (!SelectionKa(track2) && strategyPIDKaon == 0) {
          continue;
        }
        if (!SelectionPID(track2, 1) && strategyPIDKaon == 1) {
          continue;
        }
        if (numberPion == 1) {
          qaRegistry.fill(HIST("hNsigmaPtkaonTPC"), track2.tpcNSigmaKa(), track2.pt());
          if (track2.hasTOF()) {
            qaRegistry.fill(HIST("hNsigmaPtkaonTOF"), track2.tofNSigmaKa(), track2.pt());
          }
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }

        auto collisionId1 = track1.collisionId();
        auto collisionId2 = track2.collisionId();

        if (collisionId1 != collisionId2) {
          continue;
        }
        int track1Sign = track1.sign();
        int track2Sign = track2.sign();
        numberPiKpair = numberPiKpair + 1;

        for (auto track3 : V0s) {
          if (!SelectionV0(collision, track3)) {
            continue;
          }
          auto collisionId3 = track3.collisionId();
          if (collisionId1 != collisionId3) {
            continue;
          }
          // LOG(info) << "Track3" << "\n";
          if (numberPiKpair == 1) {
            qaRegistry.fill(HIST("hInvMassk0"), track3.mK0Short(), track3.pt());
          }
          pT = RecoDecay::pt(std::array{track1.px() + track2.px() + track3.px(), track1.py() + track2.py() + track3.py()});
          auto arrMomF1 = std::array{
            std::array{track1.px(), track1.py(), track1.pz()},
            std::array{track2.px(), track2.py(), track2.pz()},
            std::array{track3.px(), track3.py(), track3.pz()}};
          auto arrMom23 = std::array{
            std::array{track2.px(), track2.py(), track2.pz()},
            std::array{track3.px(), track3.py(), track3.pz()}};
          masskKs0 = RecoDecay::m(arrMom23, std::array{massKa, massK0s});
          massF1 = RecoDecay::m(arrMomF1, std::array{massPi, massKa, massK0s});
          qaRegistry.fill(HIST("hInvMassKKs0"), masskKs0);
          if ((masskKs0 > cMaxMassKKs0) || (massF1 > cMaxMassF1) || (pT < cMinF1Pt)) {
            continue;
          }
          if (track1Sign * track2Sign > 0) {
            qaRegistry.fill(HIST("hInvMassf1Like"), massF1, pT);
            continue;
          }
          qaRegistry.fill(HIST("hInvMassf1"), massF1, pT);
          numberF1 = numberF1 + 1;
          F1Vector.SetXYZM(track1.px() + track2.px() + track3.px(), track1.py() + track2.py() + track3.py(), track1.pz() + track2.pz() + track3.pz(), massF1);

          ////////////// proton loop for F1-proton trigger/////////////////
          for (auto track4 : tracks) {
            auto collisionId4 = track4.collisionId();
            if (collisionId1 != collisionId4) {
              continue;
            }
            if (!SelectionTrack(track4)) {
              continue;
            }
            if (!SelectionPr(track4) && strategyPIDProton == 0) {
              continue;
            }
            if (!SelectionPID(track4, 2) && strategyPIDProton == 1) {
              continue;
            }
            if (numberF1 == 1) {
              qaRegistry.fill(HIST("hNsigmaPtprotonTPC"), track4.tpcNSigmaPr(), track4.pt());
              if (track4.hasTOF()) {
                qaRegistry.fill(HIST("hNsigmaPtprotonTOF"), track4.tofNSigmaPr(), track4.pt());
              }
            }
            numberF1Proton = numberF1Proton + 1;
            ProtonVector.SetXYZM(track4.px(), track4.py(), track4.pz(), 0.938);
            F1ProtonVector = F1Vector + ProtonVector;
            betaX = -F1ProtonVector.X() / F1ProtonVector.E();
            betaY = -F1ProtonVector.Y() / F1ProtonVector.E();
            betaZ = -F1ProtonVector.Z() / F1ProtonVector.E();
            F1Vector.Boost(betaX, betaY, betaZ);
            ProtonVector.Boost(betaX, betaY, betaZ);
            RelativeMomentumVector = ProtonVector - F1Vector;
            relativeMomentum = 0.5 * RelativeMomentumVector.P();
            qaRegistry.fill(HIST("hkstarDist"), relativeMomentum);
            if (relativeMomentum > cMaxRelMom) {
              continue;
            }
            numberF1ProtonFemto = numberF1ProtonFemto + 1;
          }
        }
      }
    }
    if (numberF1 > 0) {
      triggerF1 = true;
      qaRegistry.fill(HIST("hEventstat"), 1.5);
    }
    if (numberF1Proton > 0) {
      triggerF1Proton = true;
      qaRegistry.fill(HIST("hEventstat"), 2.5);
    }
    if (numberF1ProtonFemto > 0) {
      triggerF1ProtonFemto = true;
      qaRegistry.fill(HIST("hEventstat"), 3.5);
    }
    f1Collisions(triggerF1, triggerF1Proton, triggerF1ProtonFemto);
  }
  PROCESS_SWITCH(f1protoninitializer, processF1Proton, "Process for trigger", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<f1protoninitializer>(cfgc, TaskName{"lf-f1protoninitializer"}),
  };
}
