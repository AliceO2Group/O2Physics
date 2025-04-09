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
///
/// \file HeavyNeutralMeson.cxx
///
/// \brief This code loops over collisions to reconstruct heavy mesons (omega or eta') using EMCal clusters and V0s (PCM)
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include <vector>
#include <iostream>
#include <iterator>
#include <string>

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TRandom3.h"

#include "PWGEM/PhotonMeson/Utils/HNMUtilities.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "fairlogger/Logger.h"
#include "Framework/Configurable.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "CommonConstants/MathConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pwgem::photonmeson;

namespace o2::aod
{
using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
using MyCollision = MyCollisions::iterator;
using SelectedTracks = soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
} // namespace o2::aod

namespace HNMPID
{

enum TracksPID {
  kPion,
  kNTracksPID
};

enum PIDLimits { kTPCMin,
                 kTPCMax,
                 kTPCTOF,
                 kITSmin,
                 kITSmax,
                 kNPIDLimits
};

const std::vector<std::string> SpeciesName{"pion"}; // ToDo include charged pions
const std::vector<std::string> PtCutsName{"Pt min", "Pt max", "P TOF thres"};
const std::vector<std::string> PidCutsName{"TPC min", "TPC max", "TPCTOF max", "ITS min", "ITS max"};
const float pidcutsTable[kNTracksPID][kNPIDLimits]{{-4.f, 4.f, 4.f, -99.f, 99.f}};
const float ptcutsTable[kNTracksPID][3]{{0.35f, 6.f, 0.75f}};
const float TPCNClustersMin[1][kNTracksPID]{{80.0f}};
const float ITSNClustersMin[1][kNTracksPID]{{4}};

} // namespace HNMPID

struct HeavyNeutralMeson {

  // PID selections
  Configurable<LabeledArray<float>> ConfPIDCuts{"ConfPIDCuts", {HNMPID::pidcutsTable[0], HNMPID::kNTracksPID, HNMPID::kNPIDLimits, HNMPID::SpeciesName, HNMPID::PidCutsName}, "Heavy Neutral Meson PID nsigma selections"};
  Configurable<LabeledArray<float>> ConfPtCuts{"ConfPtCuts", {HNMPID::ptcutsTable[0], HNMPID::kNTracksPID, 3, HNMPID::SpeciesName, HNMPID::PtCutsName}, "Heavy Neutral Meson pT selections"};
  Configurable<float> ConfTrkEta{"ConfTrkEta", 0.9, "Eta"};
  Configurable<LabeledArray<float>> ConfTPCNClustersMin{"ConfTPCNClustersMin", {HNMPID::TPCNClustersMin[0], 1, HNMPID::kNTracksPID, std::vector<std::string>{"TPCNClusMin"}, HNMPID::SpeciesName}, "Mininum of TPC Clusters"};
  Configurable<float> ConfTrkTPCfCls{"ConfTrkTPCfCls", 0.83, "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> ConfTrkTPCcRowsMin{"ConfTrkTPCcRowsMin", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> ConfTrkTPCsClsSharedFrac{"ConfTrkTPCsClsSharedFrac", 1.f, "Fraction of shared TPC clusters"};
  Configurable<LabeledArray<float>> ConfTrkITSnclsMin{"ConfTrkITSnclsMin", {HNMPID::ITSNClustersMin[0], 1, HNMPID::kNTracksPID, std::vector<std::string>{"Cut"}, HNMPID::SpeciesName}, "Minimum number of ITS clusters"};
  Configurable<float> ConfTrkDCAxyMax{"ConfTrkDCAxyMax", 0.15, "Maximum DCA_xy"};
  Configurable<float> ConfTrkDCAzMax{"ConfTrkDCAzMax", 0.3, "Maximum DCA_z"};
  Configurable<float> ConfTrkMaxChi2PerClusterTPC{"ConfTrkMaxChi2PerClusterTPC", 4.0f, "Minimal track selection: max allowed chi2 per TPC cluster"};  // 4.0 is default of global tracks on 20.01.2023
  Configurable<float> ConfTrkMaxChi2PerClusterITS{"ConfTrkMaxChi2PerClusterITS", 36.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of global tracks on 20.01.2023
  Configurable<bool> ConfEvtSelectZvtx{"ConfEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};
  Configurable<int> cfgHNMMassCorrection{"cfgHNMMassCorrection", 1, "Use GG PDG mass to correct HNM mass (0 = off, 1 = subDeltaPi0, 2 = subLambda)"};
  static constexpr float defaultMassWindows[4] = {0.1, 0.17, 0.5, 6.};
  Configurable<LabeledArray<float>> massWindowLightNeutralMesons{"massWindowLightNeutralMesons", {defaultMassWindows, 4, {"pi0_min", "pi0_max", "eta_min", "eta_max"}}, "Mass window for selected light neutral meson decay daughters"};
  Configurable<bool> ConfDoEMCShift{"ConfDoEMCShift", false, "Apply SM-wise shift in eta and phi to EMCal clusters to align with TPC tracks"};
  Configurable<std::vector<float>> ConfEMCEtaShift{"ConfEMCEtaShift", {0.f}, "values for SM-wise shift in eta to be added to EMCal clusters to align with TPC tracks"};
  Configurable<std::vector<float>> ConfEMCPhiShift{"ConfEMCPhiShift", {0.f}, "values for SM-wise shift in phi to be added to EMCal clusters to align with TPC tracks"};
  std::array<float, 20> EMCEtaShift = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::array<float, 20> EMCPhiShift = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  template <typename T>
  bool isSelectedTrack(T const& track, HNMPID::TracksPID partSpecies)
  {
    const auto pT = track.pt();
    const auto eta = track.eta();
    const auto tpcNClsF = track.tpcNClsFound();
    const auto tpcRClsC = track.tpcCrossedRowsOverFindableCls();
    const auto tpcNClsC = track.tpcNClsCrossedRows();
    const auto tpcNClsSFrac = track.tpcFractionSharedCls();
    const auto itsNCls = track.itsNCls();
    const auto dcaXY = track.dcaXY();
    const auto dcaZ = track.dcaZ();

    if (pT < ConfPtCuts->get(partSpecies, "Pt min")) {
      return false;
    }
    if (pT > ConfPtCuts->get(partSpecies, "Pt max")) {
      return false;
    }
    if (std::abs(eta) > ConfTrkEta) {
      return false;
    }
    if (tpcNClsF < ConfTPCNClustersMin->get("TPCNClusMin", partSpecies)) {
      return false;
    }
    if (tpcRClsC < ConfTrkTPCfCls) {
      return false;
    }
    if (tpcNClsC < ConfTrkTPCcRowsMin) {
      return false;
    }
    if (tpcNClsSFrac > ConfTrkTPCsClsSharedFrac) {
      return false;
    }
    if (itsNCls < ConfTrkITSnclsMin->get(static_cast<uint>(0), partSpecies)) {
      return false;
    }
    if (std::abs(dcaXY) > ConfTrkDCAxyMax) {
      return false;
    }
    if (std::abs(dcaZ) > ConfTrkDCAzMax) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedTrackPID(T const& track, HNMPID::TracksPID partSpecies)
  {
    bool isSelected = false;

    float nSigmaTrackTPC = track.tpcNSigmaPi();
    float nSigmaTrackTOF = track.tofNSigmaPi();
    float nSigmaTrackITS = track.itsNSigmaPi();

    float nSigmaTrackTPCTOF = std::sqrt(std::pow(nSigmaTrackTPC, 2) + std::pow(nSigmaTrackTOF, 2));

    // check if track is selected
    auto TPCmin = ConfPIDCuts->get(partSpecies, HNMPID::kTPCMin);
    auto TPCmax = ConfPIDCuts->get(partSpecies, HNMPID::kTPCMax);
    auto TPCTOFmax = ConfPIDCuts->get(partSpecies, HNMPID::kTPCTOF);
    auto ITSmin = ConfPIDCuts->get(partSpecies, HNMPID::kITSmin);
    auto ITSmax = ConfPIDCuts->get(partSpecies, HNMPID::kITSmax);

    if (track.p() <= ConfPtCuts->get(partSpecies, "P TOF thres")) {
      if (nSigmaTrackTPC > TPCmin &&
          nSigmaTrackTPC < TPCmax &&
          nSigmaTrackITS > ITSmin &&
          nSigmaTrackITS < ITSmax) {
        isSelected = true;
      }
    } else {
      if (nSigmaTrackTPCTOF < TPCTOFmax) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (ConfEvtSelectZvtx && std::abs(col.posZ()) > ConfEvtZvtx) {
      return false;
    }
    if (ConfEvtOfflineCheck && !col.sel8()) {
      return false;
    }
    return true;
  }

  // Circumvent missing of different phi mappings, enforce [0, 2 * M_PI]
  // Tracks have domain [0, 2 * M_PI]
  // TLorentVectors have domain [-M_PI, M_PI]
  double translatePhi(double phi)
  {
    if (phi < 0) {
      phi += 2 * M_PI; // Add 2 pi to make it positive
    }
    return phi;
  }

  HistogramRegistry mHistManager{"HeavyNeutralMesonHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;

  emcal::Geometry* emcalGeom;

  // Prepare vectors for different species
  std::vector<ROOT::Math::PtEtaPhiMVector> pion, antipion;

  void init(InitContext const&)
  {
    emcalGeom = emcal::Geometry::GetInstanceFromRunNumber(300000);
    auto hCollisionCounter = mHistManager.add<TH1>("Event/hCollisionCounter", "Number of collisions;;#bf{#it{N}_{Coll}}", HistType::kTH1F, {{6, -0.5, 5.5}});
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "kTVXinEMC");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "PCM #omega");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "EMC #omega");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "PCM #eta'");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "EMC #eta'");

    mHistManager.add("Event/nGGs", "Number of (selected) #gamma#gamma paris;#bf{#it{N}_{#gamma#gamma}};#bf{#it{N}_{#gamma#gamma}^{selected}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nTracks", "Number of tracks;#bf{N_{tracks}};#bf{#it{N}_{Coll}}", HistType::kTH1F, {{51, -0.5, 50.5}});
    mHistManager.add("Event/nHeavyNeutralMesons", "Number of (selected) HNM candidates;#bf{#it{N}_{HNM}};#bf{#it{N}_{HNM}^{selected}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nClustersVsV0s", "Number of clusters and V0s in the collision;#bf{#it{N}_{clusters}};#bf{#it{N}_{V0s}}", HistType::kTH2F, {{26, -0.5, 25.5}, {26, -0.5, 25.5}});

    mHistManager.add("GG/invMassVsPt_PCM", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_PCMEMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_EMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});

    mHistManager.add("HeavyNeutralMeson/invMassVsPt_PCM", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{pT}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
    mHistManager.add("HeavyNeutralMeson/invMassVsPt_PCMEMC", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{pT}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
    mHistManager.add("HeavyNeutralMeson/invMassVsPt_EMC", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{pT}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});

    // event cuts
    mHistManager.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("EventCuts/fZvtxBefore", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});
    mHistManager.add("EventCuts/fZvtxAfter", "Zvtx after event cuts;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    // mom correlations p vs pTPC
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationPos", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationNeg", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});

    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsPion", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiPion", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{1000, 0.0f, 20.0f}, {1000, 0.0f, 20.0f}}});

    // all tracks
    mHistManager.add("TrackCuts/TracksBefore/fPtTrackBefore", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/TracksBefore/fEtaTrackBefore", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/TracksBefore/fPhiTrackBefore", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // TPC signal
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignal", "TPCSignal;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalP", "TPCSignalP;p (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    // TPC signal anti
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAnti", "TPCSignal;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiP", "TPCSignalP;p (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});

    // TPC signal particles
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalPion", "fTPCSignalPion;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiPion", "fTPCSignalAntiPion;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    // pions
    mHistManager.add("TrackCuts/Pion/fPPion", "Momentum of Pions at PV;p (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Pion/fPTPCPion", "Momentum of Pions at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Pion/fPtPion", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Pion/fMomCorPionDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/Pion/fMomCorPionRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/Pion/fEtaPion", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/Pion/fPhiPion", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/Pion/fNsigmaTPCvsPPion", "NSigmaTPC Pion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Pion/fNsigmaTOFvsPPion", "NSigmaTOF Pion;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Pion/fNsigmaTPCTOFvsPPion", "NSigmaTPCTOF Pion;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/Pion/fNsigmaTPCvsPPionP", "NSigmaTPC Pion P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Pion/fNsigmaTOFvsPPionP", "NSigmaTOF Pion P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Pion/fNsigmaTPCTOFvsPPionP", "NSigmaTPCTOF Pion P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/Pion/fDCAxyPion", "fDCAxy Pion;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Pion/fDCAzPion", "fDCAz Pion;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Pion/fTPCsClsPion", "fTPCsCls Pion;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Pion/fTPCcRowsPion", "fTPCcRows Pion;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Pion/fTrkTPCfClsPion", "fTrkTPCfCls Pion;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/Pion/fTPCnclsPion", "fTPCncls Pion;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // anti-pions
    mHistManager.add("TrackCuts/AntiPion/fPtAntiPion", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/AntiPion/fMomCorAntiPionDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/AntiPion/fMomCorAntiPionRatio", "Momentum correlation;p_{reco} (GeV/c); |p_{TPC} - p_{reco}| (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/AntiPion/fEtaAntiPion", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/AntiPion/fPhiAntiPion", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/AntiPion/fNsigmaTPCvsPAntiPion", "NSigmaTPC AntiPion;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiPion/fNsigmaTOFvsPAntiPion", "NSigmaTOF AntiPion;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiPion/fNsigmaTPCTOFvsPAntiPion", "NSigmaTPCTOF AntiPion;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiPion/fNsigmaTPCvsPAntiPionP", "NSigmaTPC AntiPion P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiPion/fNsigmaTOFvsPAntiPionP", "NSigmaTOF AntiPion P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiPion/fNsigmaTPCTOFvsPAntiPionP", "NSigmaTPCTOF AntiPion P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiPion/fDCAxyAntiPion", "fDCAxy AntiPion;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiPion/fDCAzAntiPion", "fDCAz AntiPion;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiPion/fTPCsClsAntiPion", "fTPCsCls AntiPion;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiPion/fTPCcRowsAntiPion", "fTPCcRows AntiPion;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiPion/fTrkTPCfClsAntiPion", "fTrkTPCfCls AntiPion;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/AntiPion/fTPCnclsAntiPion", "fTPCncls AntiPion;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // HNM
    // omega QA
    // daughter pos before
    mHistManager.add("TrackCuts/HMN/Before/PosDaughter/fInvMass", "Invariant mass HMN Pos Daugh;M_{#pi};Entries", HistType::kTH1F, {{500, 0, 1}});
    mHistManager.add("TrackCuts/HMN/Before/PosDaughter/fPt", "Transverse momentum HMN Pos Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/PosDaughter/fEta", "HMN Pos Daugh Eta;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/PosDaughter/fPhi", "Azimuthal angle of HMN Pos Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // daughter neg before
    mHistManager.add("TrackCuts/HMN/Before/NegDaughter/fInvMass", "Invariant mass HMN Neg Daugh;M_{#pi};Entries", HistType::kTH1F, {{500, 0, 1}});
    mHistManager.add("TrackCuts/HMN/Before/NegDaughter/fPt", "Transverse momentum HMN Neg Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/NegDaughter/fEta", "HMN Neg Daugh Eta;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/NegDaughter/fPhi", "Azimuthal angle of HMN Neg Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // HMNCand tracks before
    mHistManager.add("TrackCuts/HMN/Before/fInvMass_tracks", "Invariant mass HMNCand;M_{#pi#pi#gammg#gamma};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/Before/fPt_tracks", "Transverse momentum HMNCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/fEta_tracks", "Pseudorapidity of HMNCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/fPhi_tracks", "Azimuthal angle of HMNCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/HMN/Before/PCM/fInvMass", "Invariant mass HMNCand;M_{#pi#pi#gammg#gamma};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/Before/PCM/fPt", "Transverse momentum HMNCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/PCM/fEta", "Pseudorapidity of HMNCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/PCM/fPhi", "Azimuthal angle of HMNCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/HMN/Before/EMC/fInvMass", "Invariant mass HMNCand;M_{#pi#pi#gammg#gamma};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/Before/EMC/fPt", "Transverse momentum HMNCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/EMC/fEta", "Pseudorapidity of HMNCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/EMC/fPhi", "Azimuthal angle of HMNCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/HMN/Before/PCMEMC/fInvMass", "Invariant mass HMNCand;M_{#pi#pi#gammg#gamma};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/Before/PCMEMC/fPt", "Transverse momentum HMNCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/Before/PCMEMC/fEta", "Pseudorapidity of HMNCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/Before/PCMEMC/fPhi", "Azimuthal angle of HMNCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    if (ConfDoEMCShift.value) {
      for (unsigned short iSM = 0; iSM < 20; iSM++) {
        EMCEtaShift[iSM] = ConfEMCEtaShift.value[iSM];
        EMCPhiShift[iSM] = ConfEMCPhiShift.value[iSM];
        LOG(info) << "SM-wise shift in eta/phi for SM " << iSM << ": " << EMCEtaShift[iSM] << " / " << EMCPhiShift[iSM];
      }
    }
  }
  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void process(aod::MyCollision const& collision, aod::MyBCs const&, aod::SkimEMCClusters const& clusters, aod::V0PhotonsKF const& v0s, aod::SelectedTracks const& tracks)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<aod::SelectedTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

    mHistManager.fill(HIST("Event/hCollisionCounter"), 0.);

    // QA all evts
    mHistManager.fill(HIST("EventCuts/fMultiplicityBefore"), collision.multNTracksPV());
    mHistManager.fill(HIST("EventCuts/fZvtxBefore"), collision.posZ());

    // Ensure evts are consistent with Sel8 and Vtx-z selection
    if (!isSelectedEvent(collision)) {
      return;
    }

    // QA accepted evts
    mHistManager.fill(HIST("EventCuts/fMultiplicityAfter"), collision.multNTracksPV());
    mHistManager.fill(HIST("EventCuts/fZvtxAfter"), collision.posZ());

    // clean vecs
    // Pions for HNM
    pion.clear();
    antipion.clear();
    vHNMs.clear();
    // vGGs vector is cleared in reconstructGGs.

    if (collision.foundBC_as<aod::MyBCs>().alias_bit(kTVXinEMC)) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 1.);
    }

    auto v0sInThisCollision = v0s.sliceBy(perCollision_pcm, collision.globalIndex());
    auto clustersInThisCollision = clusters.sliceBy(perCollision_emc, collision.globalIndex());

    mHistManager.fill(HIST("Event/nClustersVsV0s"), clustersInThisCollision.size(), v0sInThisCollision.size());
    mHistManager.fill(HIST("Event/nTracks"), tracksWithItsPid.size());

    std::vector<hnmutilities::Photon> vGammas;
    hnmutilities::storeGammasInVector(clustersInThisCollision, v0sInThisCollision, vGammas, EMCEtaShift, EMCPhiShift);
    hnmutilities::reconstructGGs(vGammas, vGGs);
    vGammas.clear();
    processGGs(vGGs);

    bool isPion = false;

    for (const auto& track : tracksWithItsPid) {
      // General QA
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fPtTrackBefore"), track.pt());
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fEtaTrackBefore"), track.eta());
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fPhiTrackBefore"), track.phi());
      // Fill PID info
      if (track.sign() > 0) {
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignal"), track.tpcInnerParam(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalP"), track.p(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationPos"), track.p(), track.tpcInnerParam());
      }
      if (track.sign() < 0) {
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAnti"), track.tpcInnerParam(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiP"), track.p(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationNeg"), track.p(), track.tpcInnerParam());
      }

      isPion = (isSelectedTrackPID(track, HNMPID::kPion) && isSelectedTrack(track, HNMPID::kPion));

      if (isPion) {
        if (track.sign() > 0) { // part
          pion.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassPionCharged);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsPion"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalPion"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/Pion/fPPion"), track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fPTPCPion"), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/Pion/fPtPion"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Pion/fMomCorPionDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fMomCorPionRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fEtaPion"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Pion/fPhiPion"), track.phi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCvsPPion"), track.tpcInnerParam(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTOFvsPPion"), track.tpcInnerParam(), track.tofNSigmaPi());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2));
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCTOFvsPPion"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCvsPPionP"), track.p(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTOFvsPPionP"), track.p(), track.tofNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCTOFvsPPionP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Pion/fDCAxyPion"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Pion/fDCAzPion"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCsClsPion"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCcRowsPion"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Pion/fTrkTPCfClsPion"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCnclsPion"), track.tpcNClsFound());
        } else { // antipart
          antipion.emplace_back(track.pt(), track.eta(), track.phi(), o2::constants::physics::MassPionCharged);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiPion"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiPion"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fPtAntiPion"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fMomCorAntiPionDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fMomCorAntiPionRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fEtaAntiPion"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fPhiAntiPion"), track.phi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCvsPAntiPion"), track.tpcInnerParam(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTOFvsPAntiPion"), track.tpcInnerParam(), track.tofNSigmaPi());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCTOFvsPAntiPion"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCvsPAntiPionP"), track.p(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTOFvsPAntiPionP"), track.p(), track.tofNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCTOFvsPAntiPionP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAxyAntiPion"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAzAntiPion"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCsClsAntiPion"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCcRowsAntiPion"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTrkTPCfClsAntiPion"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCnclsAntiPion"), track.tpcNClsFound());
        }
      }
    }

    // reconstruct HMN candidates
    for (const auto& posPion : pion) {
      for (const auto& negPion : antipion) {
        hnmutilities::reconstructHeavyNeutralMesons(posPion, negPion, vGGs, vHNMs);

        ROOT::Math::PtEtaPhiMVector temp = posPion + negPion;

        mHistManager.fill(HIST("TrackCuts/HMN/Before/fInvMass_tracks"), temp.M());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/fPt_tracks"), temp.pt());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/fEta_tracks"), temp.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/fPhi_tracks"), translatePhi(temp.phi()));

        mHistManager.fill(HIST("TrackCuts/HMN/Before/PosDaughter/fInvMass"), posPion.M());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PosDaughter/fPt"), posPion.pt());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PosDaughter/fEta"), posPion.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PosDaughter/fPhi"), translatePhi(posPion.phi()));

        mHistManager.fill(HIST("TrackCuts/HMN/Before/NegDaughter/fInvMass"), negPion.M());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/NegDaughter/fPt"), negPion.pt());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/NegDaughter/fEta"), negPion.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/NegDaughter/fPhi"), translatePhi(negPion.phi()));
      }
    }

    processHNMs(vHNMs); // Contains QA of HMN properties
  }

  /// \brief Loop over the GG candidates, fill the mass/pt histograms and set the isPi0/isEta flags based on the reconstructed mass
  void processGGs(std::vector<hnmutilities::GammaGammaPair>& vGGs)
  {
    int nGGsBeforeMassCuts = vGGs.size();
    for (unsigned int iGG = 0; iGG < vGGs.size(); iGG++) {
      auto lightMeson = &vGGs.at(iGG);

      if (lightMeson->reconstructionType == photonpair::kPCMPCM) {
        mHistManager.fill(HIST("GG/invMassVsPt_PCM"), lightMeson->m(), lightMeson->pT());
      } else if (lightMeson->reconstructionType == photonpair::kEMCEMC) {
        mHistManager.fill(HIST("GG/invMassVsPt_EMC"), lightMeson->m(), lightMeson->pT());
      } else {
        mHistManager.fill(HIST("GG/invMassVsPt_PCMEMC"), lightMeson->m(), lightMeson->pT());
      }

      if (lightMeson->m() > massWindowLightNeutralMesons->get("pi0_min") && lightMeson->m() < massWindowLightNeutralMesons->get("pi0_max")) {
        lightMeson->isPi0 = true;
      } else if (lightMeson->m() > massWindowLightNeutralMesons->get("eta_min") && lightMeson->m() < massWindowLightNeutralMesons->get("eta_max")) {
        lightMeson->isEta = true;
      } else {
        vGGs.erase(vGGs.begin() + iGG);
        iGG--;
      }
    }
    mHistManager.fill(HIST("Event/nGGs"), nGGsBeforeMassCuts, vGGs.size());
  }

  /// \brief Loop over the heavy neutral meson candidates, fill the mass/pt histograms and set the trigger flags based on the reconstructed mass
  void processHNMs(std::vector<hnmutilities::HeavyNeutralMeson>& vHNMs)
  {
    int nHNMsBeforeMassCuts = vHNMs.size();
    for (unsigned int iHNM = 0; iHNM < vHNMs.size(); iHNM++) {
      auto heavyNeutralMeson = vHNMs.at(iHNM);

      float massHNM = heavyNeutralMeson.m(cfgHNMMassCorrection);
      if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
        mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt_PCM"), massHNM, heavyNeutralMeson.pT());
        // QA
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCM/fInvMass"), massHNM);
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCM/fPt"), heavyNeutralMeson.pT());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCM/fEta"), heavyNeutralMeson.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCM/fPhi"), translatePhi(heavyNeutralMeson.phi()));
      } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
        mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt_EMC"), massHNM, heavyNeutralMeson.pT());
        // QA
        mHistManager.fill(HIST("TrackCuts/HMN/Before/EMC/fInvMass"), massHNM);
        mHistManager.fill(HIST("TrackCuts/HMN/Before/EMC/fPt"), heavyNeutralMeson.pT());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/EMC/fEta"), heavyNeutralMeson.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/EMC/fPhi"), translatePhi(heavyNeutralMeson.phi()));
      } else {
        mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt_PCMEMC"), massHNM, heavyNeutralMeson.pT());
        // QA
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCMEMC/fInvMass"), massHNM);
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCMEMC/fPt"), heavyNeutralMeson.pT());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCMEMC/fEta"), heavyNeutralMeson.eta());
        mHistManager.fill(HIST("TrackCuts/HMN/Before/PCMEMC/fPhi"), translatePhi(heavyNeutralMeson.phi()));
      }
    }
    mHistManager.fill(HIST("Event/nHeavyNeutralMesons"), nHNMsBeforeMassCuts, vHNMs.size());
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMeson>(cfgc)};
}
