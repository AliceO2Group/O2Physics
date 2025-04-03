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
/// \file HeavyNeutralMesonFilter.cxx
///
/// \brief This code loops over collisions to filter events contaning heavy mesons (omega or eta') using EMCal clusters and V0s (PCM)
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt;  Maximilian Korwieser (maximilian.korwieser@cern.ch) - Technical University Munich
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

namespace CFTrigger
{
enum CFFemtoTriggers {
  kPPOmega,
  kPPEtaPrime,
  kOmegaD,
  kEtaPrimeD,
  kOmegaP,
  kEtaPrimeP,
  kNFemtoTriggers
};

enum FemtoPartners {
  kProton,
  kDeuteron,
  kPion,
  kNFemtoPartners
};

enum PIDLimits { kTPCMin,
                 kTPCMax,
                 kTPCTOF,
                 kITSmin,
                 kITSmax,
                 kNPIDLimits
};
const std::vector<std::string> SpeciesName{"proton", "Deuteron", "pion"}; // ToDo include charged pions

const std::vector<std::string> PtCutsName{"Pt min", "Pt max", "P TOF thres"};

const std::vector<std::string> PidCutsName{"TPC min", "TPC max", "TPCTOF max", "ITS min", "ITS max"};

const std::vector<std::string> FemtoFilterNames{"PPOmega", "PPEtaPrime", "Omegad", "EtaPrimed", "OmegaP", "EtaPrimeP"};

// configs for tracks
const float pidcutsTable[kNFemtoPartners][kNPIDLimits]{
  {-4.f, 4.f, 4.f, -99.f, 99.f},
  {-4.f, 4.f, 4.f, -6.f, 6.f},
  {-4.f, 4.f, 4.f, -99.f, 99.f}};

const float ptcutsTable[kNFemtoPartners][3]{
  {0.35f, 6.f, 0.75f},
  {0.55f, 2.f, 1.2f},
  {0.35f, 6.f, 0.75f}};

const float TPCNClustersMin[1][kNFemtoPartners]{
  {80.0f, 80.0f, 80.0f}};
const float ITSNClustersMin[1][kNFemtoPartners]{
  {4, 4, 4}};

static const float triggerSwitches[1][kNFemtoTriggers]{
  {1, 1, 1, 1, 1, 1}};
const float TriggerLimits[1][kNFemtoTriggers]{
  {1.f, 1.f, 1.f, 1.f, 1.f, 1.f}};
} // namespace CFTrigger

struct HeavyNeutralMesonFilter {

  Configurable<LabeledArray<float>> ConfTriggerSwitches{
    "ConfTriggerSwitches",
    {CFTrigger::triggerSwitches[0], 1, CFTrigger::kNFemtoTriggers, std::vector<std::string>{"Switch"}, CFTrigger::FemtoFilterNames},
    "Turn on specific trigger"};

  Configurable<bool> ConfKeepTwoBody{
    "ConfKeepTwoBody",
    true,
    "Turn on specific trigger selection"};

  // PID selections
  Configurable<LabeledArray<float>>
    ConfPIDCuts{
      "ConfPIDCuts",
      {CFTrigger::pidcutsTable[0], CFTrigger::kNFemtoPartners, CFTrigger::kNPIDLimits, CFTrigger::SpeciesName, CFTrigger::PidCutsName},
      "Femtopartner PID nsigma selections"};

  Configurable<LabeledArray<float>> ConfPtCuts{
    "ConfPtCuts",
    {CFTrigger::ptcutsTable[0], CFTrigger::kNFemtoPartners, 3, CFTrigger::SpeciesName, CFTrigger::PtCutsName},
    "Femtopartner pT selections"};

  Configurable<float> ConfTrkEta{
    "ConfTrkEta",
    0.9,
    "Eta"};

  Configurable<LabeledArray<float>> ConfTPCNClustersMin{
    "ConfTPCNClustersMin",
    {CFTrigger::TPCNClustersMin[0], 1, CFTrigger::kNFemtoPartners, std::vector<std::string>{"TPCNClusMin"}, CFTrigger::SpeciesName},
    "Mininum of TPC Clusters"};

  Configurable<float> ConfTrkTPCfCls{
    "ConfTrkTPCfCls",
    0.83,
    "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> ConfTrkTPCcRowsMin{
    "ConfTrkTPCcRowsMin",
    70,
    "Minimum number of crossed TPC rows"};
  Configurable<float> ConfTrkTPCsClsSharedFrac{
    "ConfTrkTPCsClsSharedFrac",
    1.f,
    "Fraction of shared TPC clusters"};

  Configurable<LabeledArray<float>> ConfTrkITSnclsMin{
    "ConfTrkITSnclsMin",
    {CFTrigger::ITSNClustersMin[0], 1, CFTrigger::kNFemtoPartners, std::vector<std::string>{"Cut"}, CFTrigger::SpeciesName},
    "Minimum number of ITS clusters"};

  Configurable<float> ConfTrkDCAxyMax{
    "ConfTrkDCAxyMax",
    0.15,
    "Maximum DCA_xy"};
  Configurable<float> ConfTrkDCAzMax{
    "ConfTrkDCAzMax",
    0.3,
    "Maximum DCA_z"};

  Configurable<float>
    ConfTrkMaxChi2PerClusterTPC{
      "ConfTrkMaxChi2PerClusterTPC",
      4.0f,
      "Minimal track selection: max allowed chi2 per TPC cluster"}; // 4.0 is default of
                                                                    // global tracks
                                                                    // on 20.01.2023
  Configurable<float>
    ConfTrkMaxChi2PerClusterITS{
      "ConfTrkMaxChi2PerClusterITS",
      36.0f,
      "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of
                                                                    // global tracks
                                                                    // on 20.01.2023

  Configurable<LabeledArray<float>> ConfKinematicLimits{
    "ConfKstarLimits",
    {CFTrigger::TriggerLimits[0], 1, CFTrigger::kNFemtoTriggers, std::vector<std::string>{"Limit"}, CFTrigger::FemtoFilterNames},
    "hypermomentum limit for two body trigger"};

  // Configs for events
  Configurable<bool> ConfEvtSelectZvtx{
    "ConfEvtSelectZvtx",
    true,
    "Event selection includes max. z-Vertex"};
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx",
                                  10.f,
                                  "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtOfflineCheck{
    "ConfEvtOfflineCheck",
    false,
    "Evt sel: check for offline selection"};

  template <typename T>
  bool isSelectedTrack(T const& track, CFTrigger::FemtoPartners partSpecies)
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
  bool isSelectedTrackPID(T const& track, CFTrigger::FemtoPartners partSpecies)
  {
    // nSigma should have entries [proton, deuteron, pion]
    bool isSelected = false;

    float nSigmaTrackTPC = -999.f;
    float nSigmaTrackTOF = -999.f;
    float nSigmaTrackITS = -999.f;

    float nSigmaTrackTPCTOF = -999.f;

    switch (partSpecies) {
      case CFTrigger::kProton:
        nSigmaTrackTPC = track.tpcNSigmaPr();
        nSigmaTrackTOF = track.tofNSigmaPr();
        nSigmaTrackITS = track.itsNSigmaPr();
        break;
      case CFTrigger::kDeuteron:
        nSigmaTrackTPC = track.tpcNSigmaDe();
        nSigmaTrackTOF = track.tofNSigmaDe();
        nSigmaTrackITS = track.itsNSigmaDe();
        break;
      case CFTrigger::kPion:
        nSigmaTrackTPC = track.tpcNSigmaPi();
        nSigmaTrackTOF = track.tofNSigmaPi();
        nSigmaTrackITS = track.itsNSigmaPi();
        break;
      default:
        LOG(fatal) << "Particle species not known";
    }

    nSigmaTrackTPCTOF = std::sqrt(std::pow(nSigmaTrackTPC, 2) + std::pow(nSigmaTrackTOF, 2));

    // check if track is selected
    auto TPCmin = ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMin);
    auto TPCmax = ConfPIDCuts->get(partSpecies, CFTrigger::kTPCMax);
    auto TPCTOFmax = ConfPIDCuts->get(partSpecies, CFTrigger::kTPCTOF);
    auto ITSmin = ConfPIDCuts->get(partSpecies, CFTrigger::kITSmin);
    auto ITSmax = ConfPIDCuts->get(partSpecies, CFTrigger::kITSmax);

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

  ROOT::Math::PxPyPzEVector getqij(const ROOT::Math::PtEtaPhiMVector parti,
                                   const ROOT::Math::PtEtaPhiMVector partj)
  {
    ROOT::Math::PxPyPzEVector vecparti(parti);
    ROOT::Math::PxPyPzEVector vecpartj(partj);
    ROOT::Math::PxPyPzEVector trackSum = vecparti + vecpartj;
    ROOT::Math::PxPyPzEVector trackDifference = vecparti - vecpartj;
    float scaling = trackDifference.Dot(trackSum) / trackSum.Dot(trackSum);
    return trackDifference - scaling * trackSum;
  }
  float getQ3(const ROOT::Math::PtEtaPhiMVector part1,
              const ROOT::Math::PtEtaPhiMVector part2,
              const ROOT::Math::PtEtaPhiMVector part3)
  {
    ROOT::Math::PxPyPzEVector q12 = getqij(part1, part2);
    ROOT::Math::PxPyPzEVector q23 = getqij(part2, part3);
    ROOT::Math::PxPyPzEVector q31 = getqij(part3, part1);
    float Q32 = q12.M2() + q23.M2() + q31.M2();
    return sqrt(-Q32);
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

  Produces<aod::HeavyNeutralMesonFilters> tags;

  HistogramRegistry mHistManager{"HeavyNeutralMesonFilterHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> cfgHNMMassCorrection{"cfgHNMMassCorrection", 1, "Use GG PDG mass to correct HNM mass (0 = off, 1 = subDeltaPi0, 2 = subLambda)"};
  static constexpr float defaultMassWindows[2][4] = {{0., 0.4, 0.6, 1.}, {0.4, 0.8, 0.8, 1.2}};
  Configurable<LabeledArray<float>> massWindowOmega{"massWindowOmega", {defaultMassWindows[0], 4, {"pi0_min", "pi0_max", "omega_min", "omega_max"}}, "Mass window for selected omegas and their decay pi0"};
  Configurable<LabeledArray<float>> massWindowEtaPrime{"massWindowEtaPrime", {defaultMassWindows[1], 4, {"eta_min", "eta_max", "etaprime_min", "etaprime_max"}}, "Mass window for selected eta' and their decay eta"};

  static constexpr float defaultMinPts[4] = {1.8, 1.8, 2.6, 2.6};
  static constexpr float defaultFemtoMinPts[4] = {1.8, 1.8, 2.6, 2.6};

  Configurable<LabeledArray<float>> minHNMPts{"minHNMPts", {defaultMinPts, 4, {"PCM_omega", "PCM_etaprime", "EMC_omega", "EMC_etaprime"}}, "Minimum pT values for the trigger decisions (GeV/c)"};

  Configurable<LabeledArray<float>> minFemtoHNMPts{"minFemtoHNMPts", {defaultFemtoMinPts, 4, {"PCM_omega", "PCM_etaprime", "EMC_omega", "EMC_etaprime"}}, "Minimum pT values for the femto trigger decisions (GeV/c)"};

  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;

  bool colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime = false;

  emcal::Geometry* emcalGeom;

  // Femto
  // Prepare vectors for different species
  std::vector<ROOT::Math::PtEtaPhiMVector> etaPrimeEMC, etaPrimePCM, omegaEMC, omegaPCM, proton, antiproton, deuteron, antideuteron, pion, antipion;
  float mMassProton = o2::constants::physics::MassProton;
  float mMassDeuteron = o2::constants::physics::MassDeuteron;
  float mMassOmega = 0.782;
  float mMassEtaPrime = 0.957;
  float mMassPionCharged = o2::constants::physics::MassPionCharged;

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

    // include all femto histograms
    mHistManager.add("fProcessedEvents", "CF - event filtered;;Events", HistType::kTH1F, {{15, -0.5, 14.5}});
    std::vector<std::string> eventTitles = {"all", "rejected", "PPOmega", "PPEtaPrime", "Omegad", "EtaPrimed", "OmegaP", "EtaPrimeP", "kTVXinEMC", "PCM #omega", "EMC #omega", "PCM #eta'", "EMC #eta'"};
    for (size_t iBin = 0; iBin < eventTitles.size(); iBin++) {
      mHistManager.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
    }

    // event cuts
    mHistManager.add("EventCuts/fMultiplicityBefore", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("EventCuts/fMultiplicityAfter", "Multiplicity after event cuts;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("EventCuts/fZvtxBefore", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});
    mHistManager.add("EventCuts/fZvtxAfter", "Zvtx after event cuts;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    // mom correlations p vs pTPC
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationPos", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationNeg", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});

    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsProton", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiProton", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsDeuteron", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiDeuteron", "fMomCorrelation;p (GeV/c);p_{TPC} (GeV/c)", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});

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
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalProton", "fTPCSignalProton;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiProton", "fTPCSignalAntiProton;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalDeuteron", "fTPCSignalDeuteron;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiDeuteron", "fTPCSignalAntiDeuteron;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalPion", "fTPCSignalPion;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiPion", "fTPCSignalAntiPion;p_{TPC} (GeV/c);dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});

    // proton
    mHistManager.add("TrackCuts/Proton/fPProton", "Momentum of protons at PV;p (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Proton/fPTPCProton", "Momentum of protons at TPC inner wall;p_{TPC} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Proton/fPtProton", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Proton/fMomCorProtonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/Proton/fMomCorProtonRatio", "Momentum correlation;p_{reco} (GeV/c); p_{TPC} - p_{reco} / p_{reco}", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/Proton/fEtaProton", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/Proton/fPhiProton", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/Proton/fNsigmaTPCvsPProton", "NSigmaTPC Proton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Proton/fNsigmaTOFvsPProton", "NSigmaTOF Proton;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Proton/fNsigmaTPCTOFvsPProton", "NSigmaTPCTOF Proton;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    mHistManager.add("TrackCuts/Proton/fNsigmaITSvsPProton", "NSigmaITS Proton;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    mHistManager.add("TrackCuts/Proton/fNsigmaTPCvsPProtonP", "NSigmaTPC Proton P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Proton/fNsigmaTOFvsPProtonP", "NSigmaTOF Proton P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Proton/fNsigmaTPCTOFvsPProtonP", "NSigmaTPCTOF Proton P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/Proton/fDCAxyProton", "fDCAxy Proton;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Proton/fDCAzProton", "fDCAz Proton;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Proton/fTPCsClsProton", "fTPCsCls Proton;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Proton/fTPCcRowsProton", "fTPCcRows Proton;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Proton/fTrkTPCfClsProton", "fTrkTPCfCls Proton;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/Proton/fTPCnclsProton", "fTPCncls Proton;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // antiproton
    mHistManager.add("TrackCuts/AntiProton/fPtAntiProton", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/AntiProton/fMomCorAntiProtonDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/AntiProton/fMomCorAntiProtonRatio", "Momentum correlation;p_{reco} (GeV/c); |p_{TPC} - p_{reco}| (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/AntiProton/fEtaAntiProton", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/AntiProton/fPhiAntiProton", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton", "NSigmaTPC AntiProton;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton", "NSigmaTOF AntiProton;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton", "NSigmaTPCTOF AntiProton;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaITSvsPAntiProton", "NSigmaITS AntiProton;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProtonP", "NSigmaTPC AntiProton P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProtonP", "NSigmaTOF AntiProton P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProtonP", "NSigmaTPCTOF AntiProton P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiProton/fDCAxyAntiProton", "fDCAxy AntiProton;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiProton/fDCAzAntiProton", "fDCAz AntiProton;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiProton/fTPCsClsAntiProton", "fTPCsCls AntiProton;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiProton/fTPCcRowsAntiProton", "fTPCcRows AntiProton;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiProton/fTrkTPCfClsAntiProton", "fTrkTPCfCls AntiProton;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/AntiProton/fTPCnclsAntiProton", "fTPCncls AntiProton;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // deuteron
    mHistManager.add("TrackCuts/Deuteron/fPtDeuteron", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/Deuteron/fMomCorDeuteronDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/Deuteron/fMomCorDeuteronRatio", "Momentum correlation;p_{reco} (GeV/c); |p_{TPC} - p_{reco}| (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/Deuteron/fEtaDeuteron", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/Deuteron/fPhiDeuteron", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteron", "NSigmaTPC Deuteron;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaTOFvsPDeuteron", "NSigmaTOF Deuteron;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaTPCTOFvsPDeuteron", "NSigmaTPCTOF Deuteron;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaITSvsPDeuteron", "NSigmaITS Deuteron;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    mHistManager.add("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteronP", "NSigmaTPC Deuteron vd P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaTOFvsPDeuteronP", "NSigmaTOF Deuteron P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/Deuteron/fNsigmaTPCTOFvsPDeuteronP", "NSigmaTPCTOF Deuteron P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/Deuteron/fDCAxyDeuteron", "fDCAxy Deuteron;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Deuteron/fDCAzDeuteron", "fDCAz Deuteron;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/Deuteron/fTPCsClsDeuteron", "fTPCsCls Deuteron;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Deuteron/fTPCcRowsDeuteron", "fTPCcRows Deuteron;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/Deuteron/fTrkTPCfClsDeuteron", "fTrkTPCfCls Deuteron;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/Deuteron/fTPCnclsDeuteron", "fTPCncls Deuteron;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

    // antideuteron
    mHistManager.add("TrackCuts/AntiDeuteron/fPtAntiDeuteron", "Transverse momentum of all processed tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/AntiDeuteron/fMomCorAntiDeuteronDif", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fMomCorAntiDeuteronRatio", "Momentum correlation;p_{reco} (GeV/c); (p_{TPC} - p_{reco}) (GeV/c)", {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fEtaAntiDeuteron", "Pseudorapidity of all processed tracks;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/AntiDeuteron/fPhiAntiDeuteron", "Azimuthal angle of all processed tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteron", "NSigmaTPC AntiDeuteron;p_{TPC} (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTOFvsPAntiDeuteron", "NSigmaTOF AntiDeuteron;p_{TPC} (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsPAntiDeuteron", "NSigmaTPCTOF AntiDeuteron;p_{TPC} (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaITSvsPAntiDeuteron", "NSigmaITS AntiDeuteron;p (GeV/c);n#sigma_{ITS}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteronP", "NSigmaTPC AntiDeuteron P;p (GeV/c);n#sigma_{TPC}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTOFvsPAntiDeuteronP", "NSigmaTOF AntiDeuteron P;p (GeV/c);n#sigma_{TOF}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
    mHistManager.add("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsPAntiDeuteronP", "NSigmaTPCTOF AntiDeuteron P;p (GeV/c);n#sigma_{comb}", {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

    mHistManager.add("TrackCuts/AntiDeuteron/fDCAxyAntiDeuteron", "fDCAxy AntiDeuteron;DCA_{XY};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiDeuteron/fDCAzAntiDeuteron", "fDCAz AntiDeuteron;DCA_{Z};Entries", HistType::kTH1F, {{500, -0.5f, 0.5f}});
    mHistManager.add("TrackCuts/AntiDeuteron/fTPCsClsAntiDeuteron", "fTPCsCls AntiDeuteron;TPC Shared Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiDeuteron/fTPCcRowsAntiDeuteron", "fTPCcRows AntiDeuteron;TPC Crossed Rows;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});
    mHistManager.add("TrackCuts/AntiDeuteron/fTrkTPCfClsAntiDeuteron", "fTrkTPCfCls AntiDeuteron;TPC Findable/CrossedRows;Entries", HistType::kTH1F, {{500, 0.0f, 3.0f}});
    mHistManager.add("TrackCuts/AntiDeuteron/fTPCnclsAntiDeuteron", "fTPCncls AntiDeuteron;TPC Clusters;Entries", HistType::kTH1F, {{163, -1.0f, 162.0f}});

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
    // Right now we loose the information about the used tracks and gg properties after creating the HeavyNeutralMeson struct: Maybe it is better to keep track of these?
    /*// daughter pos after
    mHistManager.add("TrackCuts/omega/After/PosDaughter/fPt", "Transverse momentum omegaCand Pos Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/omega/After/PosDaughter/fEta", "omegaCandPos Daugh Eta;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/omega/After/PosDaughter/fPhi", "Azimuthal angle of omegaCand Pos Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // daughter neg after
    mHistManager.add("TrackCuts/omega/After/NegDaughter/fPt", "Transverse momentum omegaCand Neg Daugh tracks;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/omega/After/NegDaughter/fEta", "omegaCand Neg Daugh Eta;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/omega/After/NegDaughter/fPhi", "Azimuthal angle of omegaCand Neg Daugh tracks;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});*/
    // omegaCand after
    // HMNCand full information
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
    // OmegaCand
    mHistManager.add("TrackCuts/HMN/After/Omega/PCM/fInvMass", "Invariant mass omegaCand;M_{#pi#pi};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/After/Omega/PCM/fPt", "Transverse momentum omegaCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/After/Omega/PCM/fEta", "Pseudorapidity of omegaCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/After/Omega/PCM/fPhi", "Azimuthal angle of omegaCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/HMN/After/Omega/EMC/fInvMass", "Invariant mass omegaCand;M_{#pi#pi};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/After/Omega/EMC/fPt", "Transverse momentum omegaCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/After/Omega/EMC/fEta", "Pseudorapidity of omegaCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/After/Omega/EMC/fPhi", "Azimuthal angle of omegaCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    // EtaPrimeCand
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/PCM/fInvMass", "Invariant mass EtaPrimeCand;M_{#pi#pi};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/PCM/fPt", "Transverse momentum EtaPrimeCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/PCM/fEta", "Pseudorapidity of EtaPrimeCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/PCM/fPhi", "Azimuthal angle of EtaPrimeCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/EMC/fInvMass", "Invariant mass EtaPrimeCand;M_{#pi#pi};Entries", HistType::kTH1F, {{5000, 0, 5}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/EMC/fPt", "Transverse momentum EtaPrimeCand;p_{T} (GeV/c);Entries", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/EMC/fEta", "Pseudorapidity of EtaPrimeCand;#eta;Entries", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/HMN/After/EtaPrime/EMC/fPhi", "Azimuthal angle of EtaPrimeCand;#phi;Entries", HistType::kTH1F, {{720, 0, TMath::TwoPi()}});

    // Trigger combs
    //  for ppomega
    mHistManager.add("ppomega/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("ppomega/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("ppomega/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppomega/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppomega/fProtonPtVsQ3_EMC", "pT (proton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppomega/fOmegaCandPtVsQ3_EMC", "pT (omega) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppomega/fAntiProtonPtVsQ3_EMC", "pT (antiproton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});

    mHistManager.add("ppomega/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppomega/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppomega/fProtonPtVsQ3_PCM", "pT (proton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppomega/fOmegaCandPtVsQ3_PCM", "pT (omega) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppomega/fAntiProtonPtVsQ3_PCM", "pT (antiproton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    //  for ppetaprime
    mHistManager.add("ppetaprime/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("ppetaprime/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("ppetaprime/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppetaprime/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppetaprime/fProtonPtVsQ3_EMC", "pT (proton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppetaprime/fEtaPrimeCandPtVsQ3_EMC", "pT (etaprime) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppetaprime/fAntiProtonPtVsQ3_EMC", "pT (antiproton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});

    mHistManager.add("ppetaprime/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppetaprime/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("ppetaprime/fProtonPtVsQ3_PCM", "pT (proton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppetaprime/fEtaPrimeCandPtVsQ3_PCM", "pT (etaprime) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});
    mHistManager.add("ppetaprime/fAntiProtonPtVsQ3_PCM", "pT (antiproton) vs Q_{3};Q_{3} (GeV/c);p_{T} (GeV/c)", {HistType::kTH2F, {{150, 0, 1.5}, {500, 0, 10}}});

    // two body
    // omegad
    mHistManager.add("omegad/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("omegad/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("omegad/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fomegaPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fDeuteronPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fAntiDeuteronPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fomegaPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fDeuteronPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegad/fAntiDeuteronPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // etaprimed
    mHistManager.add("etaprimed/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("etaprimed/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("etaprimed/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fetaprimePtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fDeuteronPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fAntiDeuteronPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fetaprimePtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fDeuteronPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimed/fAntiDeuteronPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // omegap
    mHistManager.add("omegap/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("omegap/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("omegap/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fomegaPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fProtonPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fAntiProtonPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fomegaPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fProtonPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("omegap/fAntiProtonPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});

    // etaprimep
    mHistManager.add("etaprimep/fMultiplicity", "Multiplicity of all processed events;Mult;Entries", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("etaprimep/fZvtx", "Zvtx of all processed events;Z_{vtx};Entries", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("etaprimep/fSE_particle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fSE_Antiparticle_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fetaprimePtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fProtonPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fAntiProtonPtVskstar_PCM", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fSE_particle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fSE_Antiparticle_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fetaprimePtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fProtonPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
    mHistManager.add("etaprimep/fAntiProtonPtVskstar_EMC", "Same Event distribution", HistType::kTH1F, {{8000, 0, 8}});
  }
  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void process(aod::MyCollision const& collision, aod::MyBCs const&, aod::SkimEMCClusters const& clusters, aod::V0PhotonsKF const& v0s, aod::SelectedTracks const& tracks)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<aod::SelectedTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

    mHistManager.fill(HIST("Event/hCollisionCounter"), 0.);

    // QA all evts
    mHistManager.fill(HIST("fProcessedEvents"), 0);
    mHistManager.fill(HIST("EventCuts/fMultiplicityBefore"), collision.multNTracksPV());
    mHistManager.fill(HIST("EventCuts/fZvtxBefore"), collision.posZ());

    // Ensure evts are consistent with Sel8 and Vtx-z selection
    if (!isSelectedEvent(collision)) {
      return;
    }

    // QA accepted evts
    mHistManager.fill(HIST("EventCuts/fMultiplicityAfter"), collision.multNTracksPV());
    mHistManager.fill(HIST("EventCuts/fZvtxAfter"), collision.posZ());

    colContainsPCMOmega = colContainsEMCOmega = colContainsPCMEtaPrime = colContainsEMCEtaPrime = false;
    bool keepFemtoEvent[CFTrigger::kNFemtoTriggers] = {false, false, false, false, false, false};
    int lowMomentumMultiplets[CFTrigger::kNFemtoTriggers] = {0, 0, 0, 0, 0, 0};

    // clean vecs
    // HNM candidates
    etaPrimeEMC.clear();
    etaPrimePCM.clear();
    omegaEMC.clear();
    omegaPCM.clear();
    // Femto partners
    proton.clear();
    antiproton.clear();
    deuteron.clear();
    antideuteron.clear();
    // Pions for HNM
    pion.clear();
    antipion.clear();
    vHNMs.clear();
    // vGGs vector is cleared in reconstructGGs.

    if (collision.foundBC_as<aod::MyBCs>().alias_bit(kTVXinEMC)) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 1.);
      mHistManager.fill(HIST("fProcessedEvents"), 8);
    }

    auto v0sInThisCollision = v0s.sliceBy(perCollision_pcm, collision.globalIndex());
    auto clustersInThisCollision = clusters.sliceBy(perCollision_emc, collision.globalIndex());

    mHistManager.fill(HIST("Event/nClustersVsV0s"), clustersInThisCollision.size(), v0sInThisCollision.size());
    mHistManager.fill(HIST("Event/nTracks"), tracksWithItsPid.size());

    hnmutilities::reconstructGGs(clustersInThisCollision, v0sInThisCollision, vGGs);
    processGGs(vGGs);
    // hnmutilities::reconstructHeavyNeutralMesons(tracks, vGGs, vHNMs);
    // processHNMs(vHNMs);

    bool isProton = false;
    bool isDeuteron = false;
    bool isPion = false;

    // #femtoPart
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

      // #fill protons and deuterons
      isProton = (isSelectedTrackPID(track, CFTrigger::kProton) && isSelectedTrack(track, CFTrigger::kProton));
      isDeuteron = (isSelectedTrackPID(track, CFTrigger::kDeuteron) && isSelectedTrack(track, CFTrigger::kDeuteron));
      isPion = (isSelectedTrackPID(track, CFTrigger::kPion) && isSelectedTrack(track, CFTrigger::kPion));

      if (track.sign() > 0) { // part
        if (isProton) {
          proton.emplace_back(track.pt(), track.eta(), track.phi(), mMassProton);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsProton"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalProton"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/Proton/fPProton"), track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fPTPCProton"), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/Proton/fPtProton"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Proton/fMomCorProtonDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fMomCorProtonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fEtaProton"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Proton/fPhiProton"), track.phi());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsPProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsPProton"), track.tpcInnerParam(), track.tofNSigmaPr());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2));
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsPProton"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaITSvsPProton"), track.p(), track.itsNSigmaPr());

          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsPProtonP"), track.p(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsPProtonP"), track.p(), track.tofNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsPProtonP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Proton/fDCAxyProton"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Proton/fDCAzProton"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCsClsProton"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCcRowsProton"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Proton/fTrkTPCfClsProton"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCnclsProton"), track.tpcNClsFound());
        }
        if (isDeuteron) {
          deuteron.emplace_back(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsDeuteron"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalDeuteron"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fPtDeuteron"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fMomCorDeuteronDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fMomCorDeuteronRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fEtaDeuteron"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fPhiDeuteron"), track.phi());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteron"), track.tpcInnerParam(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTOFvsPDeuteron"), track.tpcInnerParam(), track.tofNSigmaDe());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaDe(), 2) + std::pow(track.tofNSigmaDe(), 2));
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCTOFvsPDeuteron"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaITSvsPDeuteron"), track.p(), track.itsNSigmaDe());

          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsPDeuteronP"), track.p(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTOFvsPDeuteronP"), track.p(), track.tofNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCTOFvsPDeuteronP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Deuteron/fDCAxyDeuteron"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fDCAzDeuteron"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCsClsDeuteron"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCcRowsDeuteron"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTrkTPCfClsDeuteron"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCnclsDeuteron"), track.tpcNClsFound());
        }
        if (isPion) {
          pion.emplace_back(track.pt(), track.eta(), track.phi(), mMassPionCharged);
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
        }

      } else { // antipart
        if (isProton) {
          antiproton.emplace_back(track.pt(), track.eta(), track.phi(), mMassProton);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiProton"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiProton"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fPtAntiProton"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fMomCorAntiProtonDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fMomCorAntiProtonRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fEtaAntiProton"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fPhiAntiProton"), track.phi());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProton"), track.tpcInnerParam(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProton"), track.tpcInnerParam(), track.tofNSigmaPr());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProton"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaITSvsPAntiProton"), track.p(), track.itsNSigmaPr());

          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsPAntiProtonP"), track.p(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsPAntiProtonP"), track.p(), track.tofNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsPAntiProtonP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiProton/fDCAxyAntiProton"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fDCAzAntiProton"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCsClsAntiProton"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCcRowsAntiProton"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTrkTPCfClsAntiProton"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCnclsAntiProton"), track.tpcNClsFound());
        }
        if (isDeuteron) {
          antideuteron.emplace_back(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiDeuteron"), track.p(), track.tpcInnerParam());

          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiDeuteron"), track.tpcInnerParam(), track.tpcSignal());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fPtAntiDeuteron"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fMomCorAntiDeuteronDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fMomCorAntiDeuteronRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fEtaAntiDeuteron"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fPhiAntiDeuteron"), track.phi());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteron"), track.tpcInnerParam(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTOFvsPAntiDeuteron"), track.tpcInnerParam(), track.tofNSigmaDe());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaDe(), 2) + std::pow(track.tofNSigmaDe(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsPAntiDeuteron"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaITSvsPAntiDeuteron"), track.p(), track.itsNSigmaDe());

          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsPAntiDeuteronP"), track.p(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTOFvsPAntiDeuteronP"), track.p(), track.tofNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsPAntiDeuteronP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fDCAxyAntiDeuteron"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fDCAzAntiDeuteron"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCsClsAntiDeuteron"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCcRowsAntiDeuteron"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTrkTPCfClsAntiDeuteron"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCnclsAntiDeuteron"), track.tpcNClsFound());
        }
        if (isPion) {
          antipion.emplace_back(track.pt(), track.eta(), track.phi(), mMassPionCharged);
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

    // build triplets
    float Q3 = 999.f, kstar = 999.f;
    // omega
    if (ConfTriggerSwitches->get("Switch", "PPOmega") > 0.) {
      // ppomega trigger

      for (size_t i = 0; i < proton.size(); ++i) {
        for (size_t j = i + 1; j < proton.size(); ++j) {
          const auto& Proton1 = proton[i];
          const auto& Proton2 = proton[j];
          // PCM
          for (const auto& omegaParticles : omegaPCM) {

            Q3 = getQ3(Proton1, Proton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_particle_PCM"), Q3);
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_PCM"), Q3, Proton1.Pt());
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_PCM"), Q3, Proton2.Pt());
            mHistManager.fill(HIST("ppomega/fOmegaCandPtVsQ3_PCM"), Q3, omegaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPOmega)) {
              lowMomentumMultiplets[CFTrigger::kPPOmega] += 1;
            }
          }
          // EMC
          for (const auto& omegaParticles : omegaEMC) {

            Q3 = getQ3(Proton1, Proton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_particle_EMC"), Q3);
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_EMC"), Q3, Proton1.Pt());
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_EMC"), Q3, Proton2.Pt());
            mHistManager.fill(HIST("ppomega/fOmegaCandPtVsQ3_EMC"), Q3, omegaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPOmega)) {
              lowMomentumMultiplets[CFTrigger::kPPOmega] += 1;
            }
          }
        }
      }
      // apapomega trigger
      // PCM
      for (size_t i = 0; i < antiproton.size(); ++i) {
        for (size_t j = i + 1; j < antiproton.size(); ++j) {
          const auto& antiProton1 = antiproton[i];
          const auto& antiProton2 = antiproton[j];
          // PCM
          for (const auto& omegaParticles : omegaPCM) {

            Q3 = getQ3(antiProton1, antiProton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_Antiparticle_PCM"), Q3);
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_PCM"), Q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_PCM"), Q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppomega/fOmegaCandPtVsQ3_PCM"), Q3, omegaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPOmega)) {
              lowMomentumMultiplets[CFTrigger::kPPOmega] += 1;
            }
          }
          // EMC
          for (const auto& omegaParticles : omegaEMC) {

            Q3 = getQ3(antiProton1, antiProton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_Antiparticle_EMC"), Q3);
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_EMC"), Q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_EMC"), Q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppomega/fOmegaCandPtVsQ3_EMC"), Q3, omegaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPOmega)) {
              lowMomentumMultiplets[CFTrigger::kPPOmega] += 1;
            }
          }
        }
      }
    }

    // etaprime
    if (ConfTriggerSwitches->get("Switch", "PPEtaPrime") > 0.) {
      // ppetaprimetrigger
      for (size_t i = 0; i < proton.size(); ++i) {
        for (size_t j = i + 1; j < proton.size(); ++j) {
          const auto& Proton1 = proton[i];
          const auto& Proton2 = proton[j];
          // PCM
          for (const auto& etaParticles : etaPrimePCM) {

            Q3 = getQ3(Proton1, Proton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_particle_PCM"), Q3);
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_PCM"), Q3, Proton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_PCM"), Q3, Proton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fEtaPrimeCandPtVsQ3_PCM"), Q3, etaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPEtaPrime)) {
              lowMomentumMultiplets[CFTrigger::kPPEtaPrime] += 1;
            }
          }
          // EMC
          for (const auto& etaParticles : etaPrimeEMC) {

            Q3 = getQ3(Proton1, Proton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_particle_EMC"), Q3);
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_EMC"), Q3, Proton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_EMC"), Q3, Proton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fEtaPrimeCandPtVsQ3_EMC"), Q3, etaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPEtaPrime)) {
              lowMomentumMultiplets[CFTrigger::kPPEtaPrime] += 1;
            }
          }
        }
      }
      // apapetaprime trigger
      for (size_t i = 0; i < antiproton.size(); ++i) {
        for (size_t j = i + 1; j < antiproton.size(); ++j) {
          const auto& antiProton1 = antiproton[i];
          const auto& antiProton2 = antiproton[j];
          // PCM
          for (const auto& etaParticles : etaPrimePCM) {

            Q3 = getQ3(antiProton1, antiProton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_Antiparticle_PCM"), Q3);
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_PCM"), Q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_PCM"), Q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fEtaPrimeCandPtVsQ3_PCM"), Q3, etaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPEtaPrime)) {
              lowMomentumMultiplets[CFTrigger::kPPEtaPrime] += 1;
            }
          }
          // EMC
          for (const auto& etaParticles : etaPrimeEMC) {

            Q3 = getQ3(antiProton1, antiProton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_Antiparticle_EMC"), Q3);
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_EMC"), Q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_EMC"), Q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fEtaPrimeCandPtVsQ3_EMC"), Q3, etaParticles.Pt());

            if (Q3 < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kPPEtaPrime)) {
              lowMomentumMultiplets[CFTrigger::kPPEtaPrime] += 1;
            }
          }
        }
      }
    }

    // build pairs
    if (ConfTriggerSwitches->get("Switch", "Omegad") > 0.) {
      // PCM
      //  omegad trigger
      for (auto iomega = omegaPCM.begin(); iomega != omegaPCM.end(); ++iomega) {
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {
          kstar = getkstar(*iomega, *iDeuteron);
          mHistManager.fill(HIST("omegad/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fDeuteronPtVskstar_PCM"), kstar, (*iDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaD)) {
            lowMomentumMultiplets[CFTrigger::kOmegaD] += 1;
          }
        }
        // omegaAd trigger
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) {
          kstar = getkstar(*iomega, *iAntiDeuteron);
          mHistManager.fill(HIST("omegad/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fAntiDeuteronPtVskstar_PCM"), kstar, (*iAntiDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaD)) {
            lowMomentumMultiplets[CFTrigger::kOmegaD] += 1;
          }
        }
      }
      // EMC
      // omegad trigger
      for (auto iomega = omegaEMC.begin(); iomega != omegaEMC.end(); ++iomega) {
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {
          kstar = getkstar(*iomega, *iDeuteron);
          mHistManager.fill(HIST("omegad/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fDeuteronPtVskstar_EMC"), kstar, (*iDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaD)) {
            lowMomentumMultiplets[CFTrigger::kOmegaD] += 1;
          }
        }
        // omegaAd trigger
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) {
          kstar = getkstar(*iomega, *iAntiDeuteron);
          mHistManager.fill(HIST("omegad/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fAntiDeuteronPtVskstar_EMC"), kstar, (*iAntiDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaD)) {
            lowMomentumMultiplets[CFTrigger::kOmegaD] += 1;
          }
        }
      }
    }
    if (ConfTriggerSwitches->get("Switch", "EtaPrimed") > 0.) {
      // PCM
      // etaPrimed trigger
      for (auto ietaprime = etaPrimePCM.begin(); ietaprime != etaPrimePCM.end(); ++ietaprime) {
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {
          kstar = getkstar(*ietaprime, *iDeuteron);
          mHistManager.fill(HIST("etaprimed/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fDeuteronPtVskstar_PCM"), kstar, (*iDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeD)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeD] += 1;
          }
        }
        // etaPrimeAd trigger
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) {
          kstar = getkstar(*ietaprime, *iAntiDeuteron);
          mHistManager.fill(HIST("etaprimed/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fAntiDeuteronPtVskstar_PCM"), kstar, (*iAntiDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeD)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeD] += 1;
          }
        }
      }
      // EMC
      // etaPrimed trigger
      for (auto ietaprime = etaPrimeEMC.begin(); ietaprime != etaPrimeEMC.end(); ++ietaprime) {
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {
          kstar = getkstar(*ietaprime, *iDeuteron);
          mHistManager.fill(HIST("etaprimed/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fDeuteronPtVskstar_EMC"), kstar, (*iDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeD)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeD] += 1;
          }
        }
        // etaPrimeAd trigger
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) {
          kstar = getkstar(*ietaprime, *iAntiDeuteron);
          mHistManager.fill(HIST("etaprimed/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fAntiDeuteronPtVskstar_EMC"), kstar, (*iAntiDeuteron).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeD)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeD] += 1;
          }
        }
      }
    }
    if (ConfTriggerSwitches->get("Switch", "OmegaP") > 0.) {
      // PCM
      //  omegap trigger
      for (auto iomega = omegaPCM.begin(); iomega != omegaPCM.end(); ++iomega) {
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {
          kstar = getkstar(*iomega, *iProton);
          mHistManager.fill(HIST("omegap/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fProtonPtVskstar_PCM"), kstar, (*iProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaP)) {
            lowMomentumMultiplets[CFTrigger::kOmegaP] += 1;
          }
        }
        // omegaAp trigger
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) {
          kstar = getkstar(*iomega, *iAntiProton);
          mHistManager.fill(HIST("omegap/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fAntiProtonPtVskstar_PCM"), kstar, (*iAntiProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaP)) {
            lowMomentumMultiplets[CFTrigger::kOmegaP] += 1;
          }
        }
      }
      // EMC
      // omegap trigger
      for (auto iomega = omegaEMC.begin(); iomega != omegaEMC.end(); ++iomega) {
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {
          kstar = getkstar(*iomega, *iProton);
          mHistManager.fill(HIST("omegap/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fProtonPtVskstar_EMC"), kstar, (*iProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaP)) {
            lowMomentumMultiplets[CFTrigger::kOmegaP] += 1;
          }
        }
        // omegaAp trigger
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) {
          kstar = getkstar(*iomega, *iAntiProton);
          mHistManager.fill(HIST("omegap/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fAntiProtonPtVskstar_EMC"), kstar, (*iAntiProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kOmegaP)) {
            lowMomentumMultiplets[CFTrigger::kOmegaP] += 1;
          }
        }
      }
    }
    if (ConfTriggerSwitches->get("Switch", "EtaPrimeP") > 0.) {
      // PCM
      // etaPrimep trigger
      for (auto ietaprime = etaPrimePCM.begin(); ietaprime != etaPrimePCM.end(); ++ietaprime) {
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {
          kstar = getkstar(*ietaprime, *iProton);
          mHistManager.fill(HIST("etaprimep/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fProtonPtVskstar_PCM"), kstar, (*iProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeP)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeP] += 1;
          }
        }
        // etaPrimeAp trigger
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) {
          kstar = getkstar(*ietaprime, *iAntiProton);
          mHistManager.fill(HIST("etaprimep/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fAntiProtonPtVskstar_PCM"), kstar, (*iAntiProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeP)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeP] += 1;
          }
        }
      }
      // EMC
      // etaPrimep trigger
      for (auto ietaprime = etaPrimeEMC.begin(); ietaprime != etaPrimeEMC.end(); ++ietaprime) {
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {
          kstar = getkstar(*ietaprime, *iProton);
          mHistManager.fill(HIST("etaprimep/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fProtonPtVskstar_EMC"), kstar, (*iProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeP)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeP] += 1;
          }
        }
        // etaPrimeAp trigger
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) {
          kstar = getkstar(*ietaprime, *iAntiProton);
          mHistManager.fill(HIST("etaprimep/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fAntiProtonPtVskstar_EMC"), kstar, (*iAntiProton).Pt());
          if (kstar < ConfKinematicLimits->get(static_cast<uint>(0), CFTrigger::kEtaPrimeP)) {
            lowMomentumMultiplets[CFTrigger::kEtaPrimeP] += 1;
          }
        }
      }
    }

    // create tags for three body triggers
    if (lowMomentumMultiplets[CFTrigger::kPPOmega] > 0) {
      keepFemtoEvent[CFTrigger::kPPOmega] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 2);
      mHistManager.fill(HIST("ppomega/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("ppomega/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[CFTrigger::kPPEtaPrime] > 0) {
      keepFemtoEvent[CFTrigger::kPPEtaPrime] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 3);
      mHistManager.fill(HIST("ppetaprime/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("ppetaprime/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[CFTrigger::kOmegaD] > 0) {
      keepFemtoEvent[CFTrigger::kOmegaD] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 4);
      mHistManager.fill(HIST("omegad/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("omegad/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[CFTrigger::kEtaPrimeD] > 0) {
      keepFemtoEvent[CFTrigger::kEtaPrimeD] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 5);
      mHistManager.fill(HIST("etaprimed/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("etaprimed/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[CFTrigger::kOmegaP] > 0) {
      keepFemtoEvent[CFTrigger::kOmegaP] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 6);
      mHistManager.fill(HIST("omegap/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("omegap/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[CFTrigger::kEtaPrimeP] > 0) {
      keepFemtoEvent[CFTrigger::kEtaPrimeP] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 7);
      mHistManager.fill(HIST("etaprimep/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("etaprimep/fZvtx"), collision.posZ());
    }

    // #set flag for tag
    if (ConfKeepTwoBody.value) {
      tags(colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime,
           keepFemtoEvent[CFTrigger::kPPOmega] || keepFemtoEvent[CFTrigger::kOmegaP], keepFemtoEvent[CFTrigger::kPPEtaPrime] || keepFemtoEvent[CFTrigger::kEtaPrimeP],
           keepFemtoEvent[CFTrigger::kOmegaD], keepFemtoEvent[CFTrigger::kEtaPrimeD]);
    } else {
      tags(colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime,
           keepFemtoEvent[CFTrigger::kPPOmega], keepFemtoEvent[CFTrigger::kPPEtaPrime],
           keepFemtoEvent[CFTrigger::kOmegaD], keepFemtoEvent[CFTrigger::kEtaPrimeD]);
    }

    if (!keepFemtoEvent[CFTrigger::kPPOmega] && !keepFemtoEvent[CFTrigger::kOmegaP] && !keepFemtoEvent[CFTrigger::kPPEtaPrime] && !keepFemtoEvent[CFTrigger::kEtaPrimeP] &&
        !keepFemtoEvent[CFTrigger::kOmegaD] && !keepFemtoEvent[CFTrigger::kEtaPrimeD]) {
      mHistManager.fill(HIST("fProcessedEvents"), 1);
    }
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

      if (lightMeson->m() > massWindowOmega->get("pi0_min") && lightMeson->m() < massWindowOmega->get("pi0_max")) {
        lightMeson->isPi0 = true;
      } else if (lightMeson->m() > massWindowEtaPrime->get("eta_min") && lightMeson->m() < massWindowEtaPrime->get("eta_max")) {
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

      if (heavyNeutralMeson.gg->isPi0 && massHNM > massWindowOmega->get("omega_min") && massHNM < massWindowOmega->get("omega_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          if (heavyNeutralMeson.pT() > minFemtoHNMPts->get("PCM_omega")) {
            omegaPCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), translatePhi(heavyNeutralMeson.phi()), mMassOmega);
            // QA
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/PCM/fInvMass"), massHNM);
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/PCM/fPt"), heavyNeutralMeson.pT());
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/PCM/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/PCM/fPhi"), translatePhi(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > minHNMPts->get("PCM_omega")) {
            colContainsPCMOmega = true;
          }
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          if (heavyNeutralMeson.pT() > minFemtoHNMPts->get("EMC_omega")) {
            omegaEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), translatePhi(heavyNeutralMeson.phi()), mMassOmega);
            // QA
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/EMC/fInvMass"), massHNM);
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/EMC/fPt"), heavyNeutralMeson.pT());
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/EMC/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("TrackCuts/HMN/After/Omega/EMC/fPhi"), translatePhi(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > minHNMPts->get("EMC_omega")) {
            colContainsEMCOmega = true;
          }
        }
      } else if (heavyNeutralMeson.gg->isEta && massHNM > massWindowEtaPrime->get("etaprime_min") && massHNM < massWindowEtaPrime->get("etaprime_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          if (heavyNeutralMeson.pT() > minFemtoHNMPts->get("PCM_etaprime")) {
            etaPrimePCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), translatePhi(heavyNeutralMeson.phi()), mMassEtaPrime);
            // QA
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/PCM/fInvMass"), massHNM);
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/PCM/fPt"), heavyNeutralMeson.pT());
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/PCM/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/PCM/fPhi"), translatePhi(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > minHNMPts->get("PCM_etaprime")) {
            colContainsPCMEtaPrime = true;
          }
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          if (heavyNeutralMeson.pT() > minFemtoHNMPts->get("EMC_etaprime")) {
            etaPrimeEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), translatePhi(heavyNeutralMeson.phi()), mMassEtaPrime);
            // QA
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/EMC/fInvMass"), massHNM);
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/EMC/fPt"), heavyNeutralMeson.pT());
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/EMC/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("TrackCuts/HMN/After/EtaPrime/EMC/fPhi"), translatePhi(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > minHNMPts->get("EMC_etaprime")) {
            colContainsEMCEtaPrime = true;
          }
        }
      } else {
        vHNMs.erase(vHNMs.begin() + iHNM);
        iHNM--;
      }
    }
    mHistManager.fill(HIST("Event/nHeavyNeutralMesons"), nHNMsBeforeMassCuts, vHNMs.size());

    if (colContainsPCMOmega) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 2.);
      mHistManager.fill(HIST("fProcessedEvents"), 9);
    }
    if (colContainsEMCOmega) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 3.);
      mHistManager.fill(HIST("fProcessedEvents"), 10);
    }
    if (colContainsPCMEtaPrime) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 4.);
      mHistManager.fill(HIST("fProcessedEvents"), 11);
    }
    if (colContainsEMCEtaPrime) {
      mHistManager.fill(HIST("Event/hCollisionCounter"), 5.);
      mHistManager.fill(HIST("fProcessedEvents"), 12);
    }
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMesonFilter>(cfgc)};
}
