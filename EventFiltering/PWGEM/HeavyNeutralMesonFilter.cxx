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
/// \brief This code loops over collisions to filter events contaning heavy neutral mesons (omega or eta') using EMCal clusters and V0s (PCM)
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt; Maximilian Korwieser (maximilian.korwieser@cern.ch) - Technical University Munich
///

#include "EventFiltering/filterTables.h"
//
#include "PWGEM/PhotonMeson/Utils/HNMUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/Boost.h>
#include <Math/Vector4D.h> // IWYU pragma: keep
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TString.h>

#include <sys/types.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pwgem::photonmeson;

namespace o2::aod
{
using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::EMCALMatchedCollisions>;
using MyCollision = MyCollisions::iterator;
using SelectedTracks = soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA,
                                 aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                 aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe>;
} // namespace o2::aod

namespace hnmtrigger
{
enum FemtoTriggers {
  kPPOmega,
  kPPEtaPrime,
  kOmegaD,
  kEtaPrimeD,
  kOmegaP,
  kEtaPrimeP,
  kNFemtoTriggers
};

enum TracksPID {
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
const std::vector<std::string> speciesName{"proton", "Deuteron", "pion"};
const std::vector<std::string> pTCutsName{"Pt min", "Pt max", "P TOF thres"};
const std::vector<std::string> pidCutsName{"TPC min", "TPC max", "TPCTOF max", "ITS min", "ITS max"};
const std::vector<std::string> femtoFilterNames{"PPOmega", "PPEtaPrime", "Omegad", "EtaPrimed", "OmegaP", "EtaPrimeP"};

// configs for tracks
// these are need [[maybe_unused]] to silence a warning from clangd, since the compiler will inline them directly to the configs down below and then say: Variable 'X' is not needed and will not be emitted
[[maybe_unused]] const float pidcutsTable[kNFemtoPartners][kNPIDLimits]{
  {-4.f, 4.f, 4.f, -99.f, 99.f},
  {-4.f, 4.f, 4.f, -6.f, 6.f},
  {-4.f, 4.f, 4.f, -99.f, 99.f}};

[[maybe_unused]] const float ptcutsTable[kNFemtoPartners][3]{
  {0.35f, 6.f, 0.75f},
  {0.55f, 2.f, 1.2f},
  {0.35f, 6.f, 0.75f}};

[[maybe_unused]] const float nClusterMinTPC[1][kNFemtoPartners]{{80.0f, 80.0f, 80.0f}};
[[maybe_unused]] const float nClusterMinITS[1][kNFemtoPartners]{{4, 4, 4}};

[[maybe_unused]] static const float triggerSwitches[1][kNFemtoTriggers]{{1, 1, 1, 1, 1, 1}};
[[maybe_unused]] const float triggerLimits[1][kNFemtoTriggers]{{1.f, 1.f, 1.f, 1.f, 1.f, 1.f}};
} // namespace hnmtrigger

struct HeavyNeutralMesonFilter {
  Produces<aod::HeavyNeutralMesonFilters> tags;

  // --------------------------------> Configurables <------------------------------------
  // - Event selection cuts
  // - Track selection cuts
  // - Cluster shifts
  // - HNM mass selection windows
  // - HNM min pTs / k*'s
  // -------------------------------------------------------------------------------------
  // ---> Event selection
  Configurable<bool> confEvtSelectZvtx{"confEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtRequireSel8{"confEvtRequireSel8", false, "Evt sel: check for offline selection (sel8)"};
  Configurable<bool> confEvtRequireTVX{"confEvtRequireTVX", false, "Evt sel: require the TVX trigger"};
  Configurable<bool> confEvtRequireNoTFBorder{"confEvtRequireNoTFBorder", false, "Evt sel: exclude time frame border edges"};
  Configurable<bool> confEvtRequireNoITSROFBorder{"confEvtRequireNoITSROFBorder", false, "Evt sel: exclude ITS readout frame border edges"};

  // ---> Track selection
  Configurable<LabeledArray<float>> cfgPtCuts{"cfgPtCuts", {hnmtrigger::ptcutsTable[0], hnmtrigger::kNFemtoPartners, 3, hnmtrigger::speciesName, hnmtrigger::pTCutsName}, "Track pT selections"};
  Configurable<float> cfgTrkEta{"cfgTrkEta", 0.9, "Eta"};
  Configurable<LabeledArray<float>> cfgTPCNClustersMin{"cfgTPCNClustersMin", {hnmtrigger::nClusterMinTPC[0], 1, hnmtrigger::kNFemtoPartners, std::vector<std::string>{"TPCNClusMin"}, hnmtrigger::speciesName}, "Mininum of TPC Clusters"};
  Configurable<float> cfgTrkTPCfCls{"cfgTrkTPCfCls", 0.83, "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> cfgTrkTPCcRowsMin{"cfgTrkTPCcRowsMin", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgTrkTPCsClsSharedFrac{"cfgTrkTPCsClsSharedFrac", 1.f, "Fraction of shared TPC clusters"};
  Configurable<LabeledArray<float>> cfgTrkITSnclsMin{"cfgTrkITSnclsMin", {hnmtrigger::nClusterMinITS[0], 1, hnmtrigger::kNFemtoPartners, std::vector<std::string>{"Cut"}, hnmtrigger::speciesName}, "Minimum number of ITS clusters"};
  Configurable<float> cfgTrkDCAxyMax{"cfgTrkDCAxyMax", 0.15, "Maximum DCA_xy"};
  Configurable<float> cfgTrkDCAzMax{"cfgTrkDCAzMax", 0.3, "Maximum DCA_z"};
  Configurable<float> cfgTrkMaxChi2PerClusterTPC{"cfgTrkMaxChi2PerClusterTPC", 4.0f, "Minimal track selection: max allowed chi2 per TPC cluster"};  // 4.0 is default of global tracks on 20.01.2023
  Configurable<float> cfgTrkMaxChi2PerClusterITS{"cfgTrkMaxChi2PerClusterITS", 36.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of global tracks on 20.01.2023

  Configurable<LabeledArray<float>> cfgPIDCuts{"cfgPIDCuts", {hnmtrigger::pidcutsTable[0], hnmtrigger::kNFemtoPartners, hnmtrigger::kNPIDLimits, hnmtrigger::speciesName, hnmtrigger::pidCutsName}, "Femtopartner PID nsigma selections"}; // PID selections

  // ---> Configurables to allow for a shift in eta/phi of EMCal clusters to better align with extrapolated TPC tracks
  Configurable<bool> cfgDoEMCShift{"cfgDoEMCShift", false, "Apply SM-wise shift in eta and phi to EMCal clusters to align with TPC tracks"};
  Configurable<std::vector<float>> cfgEMCEtaShift{"cfgEMCEtaShift", {0.f}, "values for SM-wise shift in eta to be added to EMCal clusters to align with TPC tracks"};
  Configurable<std::vector<float>> cfgEMCPhiShift{"cfgEMCPhiShift", {0.f}, "values for SM-wise shift in phi to be added to EMCal clusters to align with TPC tracks"};
  static const int nSMs = 20;
  std::array<float, nSMs> emcEtaShift = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  std::array<float, nSMs> emcPhiShift = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  // ---> Shift the omega/eta' mass based on the difference of the reconstructed mass of the pi0/eta to its PDG mass to reduce smearing caused by EMCal/PCM in photon measurement
  Configurable<int> cfgHNMMassCorrection{"cfgHNMMassCorrection", 1, "Use GG PDG mass to correct HNM mass (0 = off, 1 = subDeltaPi0, 2 = subLambda)"};

  // ---> Mass windows for the selection of heavy neutral mesons (also based on mass of their light neutral meson decay daughter)
  static constexpr float DefaultMassWindows[2][4] = {{0., 0.4, 0.6, 1.}, {0.4, 0.8, 0.8, 1.2}};
  Configurable<LabeledArray<float>> cfgMassWindowOmega{"cfgMassWindowOmega", {DefaultMassWindows[0], 4, {"pi0_min", "pi0_max", "omega_min", "omega_max"}}, "Mass window for selected omegas and their decay pi0"};
  Configurable<LabeledArray<float>> cfgMassWindowEtaPrime{"cfgMassWindowEtaPrime", {DefaultMassWindows[1], 4, {"eta_min", "eta_max", "etaprime_min", "etaprime_max"}}, "Mass window for selected eta' and their decay eta"};

  // ---> Minimum pT values for the trigger decisions of the spectra and femto trigger. The femto triggers additionally require a given k*/Q3
  static constexpr float DefaultSpectraMinPts[4] = {1.8, 1.8, 2.6, 2.6};
  static constexpr float DefaultFemtoMinPts[4] = {1.8, 1.8, 2.6, 2.6};
  Configurable<LabeledArray<float>> cfgMinHNMPtsSpectrumTrigger{"cfgMinHNMPtsSpectrumTrigger", {DefaultSpectraMinPts, 4, {"PCM_omega", "PCM_etaprime", "EMC_omega", "EMC_etaprime"}}, "Minimum pT values for the spetra trigger decisions (GeV/c)"};
  Configurable<LabeledArray<float>> cfgMinHNMPtsFemtoTrigger{"cfgMinHNMPtsFemtoTrigger", {DefaultFemtoMinPts, 4, {"PCM_omega", "PCM_etaprime", "EMC_omega", "EMC_etaprime"}}, "Minimum pT values for the femto trigger decisions (GeV/c)"};
  Configurable<LabeledArray<float>> cfgKinematicLimits{"cfgKinematicLimits", {hnmtrigger::triggerLimits[0], 1, hnmtrigger::kNFemtoTriggers, std::vector<std::string>{"Limit"}, hnmtrigger::femtoFilterNames}, "Maximum K* (Q_3) for two (three) body femto trigger"};

  Configurable<LabeledArray<float>> cfgTriggerSwitches{"cfgTriggerSwitches", {hnmtrigger::triggerSwitches[0], 1, hnmtrigger::kNFemtoTriggers, std::vector<std::string>{"Switch"}, hnmtrigger::femtoFilterNames}, "Turn on specific trigger"};

  HistogramRegistry mHistManager{"HeavyNeutralMesonFilterHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Prepare vectors for different species
  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;
  std::vector<ROOT::Math::PtEtaPhiMVector> etaPrimeEMC, etaPrimePCM, omegaEMC, omegaPCM, proton, antiproton, deuteron, antideuteron, pion, antipion;
  float mMassProton = constants::physics::MassProton;
  float mMassDeuteron = constants::physics::MassDeuteron;
  float mMassPionCharged = constants::physics::MassPionCharged;

  Preslice<aod::V0PhotonsKF> perCollisionPCM = aod::v0photonkf::collisionId;
  Preslice<aod::SkimEMCClusters> perCollisionEMC = aod::skimmedcluster::collisionId;

  bool colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime = false;

  template <typename T>
  bool isSelectedTrack(T const& track, hnmtrigger::TracksPID partSpecies)
  {
    if (track.pt() < cfgPtCuts->get(partSpecies, "Pt min"))
      return false;
    if (track.pt() > cfgPtCuts->get(partSpecies, "Pt max"))
      return false;
    if (std::abs(track.eta()) > cfgTrkEta)
      return false;
    if (track.tpcNClsFound() < cfgTPCNClustersMin->get("TPCNClusMin", partSpecies))
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgTrkTPCfCls)
      return false;
    if (track.tpcNClsCrossedRows() < cfgTrkTPCcRowsMin)
      return false;
    if (track.tpcFractionSharedCls() > cfgTrkTPCsClsSharedFrac)
      return false;
    if (track.itsNCls() < cfgTrkITSnclsMin->get(static_cast<uint>(0), partSpecies))
      return false;
    if (std::abs(track.dcaXY()) > cfgTrkDCAxyMax)
      return false;
    if (std::abs(track.dcaZ()) > cfgTrkDCAzMax)
      return false;
    if (track.tpcChi2NCl() > cfgTrkMaxChi2PerClusterTPC)
      return false;
    if (track.itsChi2NCl() > cfgTrkMaxChi2PerClusterITS)
      return false;
    return true;
  }

  template <typename T>
  bool isSelectedTrackPID(T const& track, hnmtrigger::TracksPID partSpecies)
  {
    // nSigma should have entries [proton, deuteron, pion]
    bool isSelected = false;

    float nSigmaTrackTPC = -999.f;
    float nSigmaTrackTOF = -999.f;
    float nSigmaTrackITS = -999.f;

    switch (partSpecies) {
      case hnmtrigger::kProton:
        nSigmaTrackTPC = track.tpcNSigmaPr();
        nSigmaTrackTOF = track.tofNSigmaPr();
        nSigmaTrackITS = track.itsNSigmaPr();
        break;
      case hnmtrigger::kDeuteron:
        nSigmaTrackTPC = track.tpcNSigmaDe();
        nSigmaTrackTOF = track.tofNSigmaDe();
        nSigmaTrackITS = track.itsNSigmaDe();
        break;
      case hnmtrigger::kPion:
        nSigmaTrackTPC = track.tpcNSigmaPi();
        nSigmaTrackTOF = track.tofNSigmaPi();
        nSigmaTrackITS = track.itsNSigmaPi();
        break;
      default:
        LOG(fatal) << "Particle species not known";
    }

    float nSigmaTrackTPCTOF = std::sqrt(std::pow(nSigmaTrackTPC, 2) + std::pow(nSigmaTrackTOF, 2));

    if (track.p() <= cfgPtCuts->get(partSpecies, "P TOF thres")) {
      if (nSigmaTrackTPC > cfgPIDCuts->get(partSpecies, hnmtrigger::kTPCMin) &&
          nSigmaTrackTPC < cfgPIDCuts->get(partSpecies, hnmtrigger::kTPCMax) &&
          nSigmaTrackITS > cfgPIDCuts->get(partSpecies, hnmtrigger::kITSmin) &&
          nSigmaTrackITS < cfgPIDCuts->get(partSpecies, hnmtrigger::kITSmax)) {
        isSelected = true;
      }
    } else {
      if (nSigmaTrackTPCTOF < cfgPIDCuts->get(partSpecies, hnmtrigger::kTPCTOF)) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (confEvtSelectZvtx && std::abs(col.posZ()) > confEvtZvtx)
      return false;
    if (confEvtRequireSel8 && !col.sel8())
      return false;
    if (confEvtRequireTVX && !col.selection_bit(aod::evsel::kIsTriggerTVX))
      return false;
    if (confEvtRequireNoTFBorder && !col.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;
    if (confEvtRequireNoITSROFBorder && !col.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    return true;
  }

  float getkstar(const ROOT::Math::PtEtaPhiMVector part1,
                 const ROOT::Math::PtEtaPhiMVector part2)
  {
    const ROOT::Math::PtEtaPhiMVector trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    ROOT::Math::PxPyPzMVector partOneCMS(part1);
    ROOT::Math::PxPyPzMVector partTwoCMS(part2);
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    partOneCMS = boostPRF(partOneCMS);
    partTwoCMS = boostPRF(partTwoCMS);
    const ROOT::Math::PxPyPzMVector trackRelK = partOneCMS - partTwoCMS;
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
    float q32 = q12.M2() + q23.M2() + q31.M2();
    return std::sqrt(-q32);
  }

  void init(InitContext const&)
  {
    mHistManager.add("Event/nGGs", "Number of (selected) #gamma#gamma paris;#bf{#it{N}^{#gamma#gamma}};#bf{#it{N}_{selected}^{#gamma#gamma}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nHeavyNeutralMesons", "Number of (selected) HNM candidates;#bf{#it{N}^{HNM}};#bf{#it{N}_{selected}^{HNM}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nClustersVsV0s", "Number of clusters and V0s in the collision;#bf{#it{N}^{clusters}};#bf{#it{N}^{V0s}}", HistType::kTH2F, {{26, -0.5, 25.5}, {26, -0.5, 25.5}});
    mHistManager.add("Event/nEMCalEvents", "Number of collisions with a certain combination of EMCal triggers;;#bf{#it{N}_{collisions}}", HistType::kTH1F, {{5, -0.5, 4.5}});
    std::vector<std::string> nEventTitles = {"Cells & kTVXinEMC", "Cells & L0", "Cells & !kTVXinEMC & !L0", "!Cells & kTVXinEMC", "!Cells & L0"};
    for (size_t iBin = 0; iBin < nEventTitles.size(); iBin++)
      mHistManager.get<TH1>(HIST("Event/nEMCalEvents"))->GetXaxis()->SetBinLabel(iBin + 1, nEventTitles[iBin].data());
    mHistManager.add("Event/fMultiplicityBefore", "Multiplicity of all processed events;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("Event/fMultiplicityAfter", "Multiplicity after event cuts;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("Event/fZvtxBefore", "Zvtx of all processed events;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});
    mHistManager.add("Event/fZvtxAfter", "Zvtx after event cuts;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});
    mHistManager.add("fProcessedEvents", "CF - event filtered;;Events", HistType::kTH1F, {{12, -0.5, 11.5}});
    std::vector<std::string> pEventTitles = {"all", "rejected", "PCM #omega", "EMC #omega", "PCM #eta'", "EMC #eta'", "PPOmega", "PPEtaPrime", "Omegad", "EtaPrimed", "OmegaP", "EtaPrimeP"};
    for (size_t iBin = 0; iBin < pEventTitles.size(); iBin++)
      mHistManager.get<TH1>(HIST("fProcessedEvents"))->GetXaxis()->SetBinLabel(iBin + 1, pEventTitles[iBin].data());

    mHistManager.add("GG/invMassVsPt_PCM", "Invariant mass and pT of gg candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_PCMEMC", "Invariant mass and pT of gg candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_EMC", "Invariant mass and pT of gg candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});

    // Momentum correlations p vs p_TPC
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationPos", "fMomCorrelation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationNeg", "fMomCorrelation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});

    // All tracks
    mHistManager.add("TrackCuts/TracksBefore/fPtTrackBefore", "Transverse momentum of all processed tracks;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}_{tracks}}", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/TracksBefore/fEtaTrackBefore", "Pseudorapidity of all processed tracks;#eta;#bf{#it{N}_{tracks}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/TracksBefore/fPhiTrackBefore", "Azimuthal angle of all processed tracks;#phi;#bf{#it{N}_{tracks}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

    // TPC signal
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalTPCP", "TPCSignal;#bf{#it{p}_{TPC} (GeV/#it{c})};#bf{TPC d#it{E}/d#it{x}}", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignal", "TPCSignalP;#bf{#it{p} (GeV/#it{c})};#bf{TPC d#it{E}/d#it{x}}", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    // TPC signal antiparticles (negative charge)
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiTPCP", "TPCSignal;#bf{#it{p}_{TPC} (GeV/#it{c})};#bf{TPC d#it{E}/d#it{x}}", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAnti", "TPCSignalP;#bf{#it{p} (GeV/#it{c})};#bf{TPC d#it{E}/d#it{x}}", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});

    const int nTrackSpecies = 2 * hnmtrigger::kNFemtoPartners; // x2 because of anti particles
    const char* particleSpecies[nTrackSpecies] = {"Proton", "AntiProton", "Deuteron", "AntiDeuteron", "Pion", "AntiPion"};
    const char* particleSpeciesLatex[nTrackSpecies] = {"p", "#bar{p}", "d", "#bar{d}", "#pi^{+}", "#pi^{-}"};

    for (int iParticle = 0; iParticle < nTrackSpecies; iParticle++) {
      mHistManager.add(Form("TrackCuts/TracksBefore/fMomCorrelationAfterCuts%s", particleSpecies[iParticle]), Form("%s momentum correlation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", particleSpecies[iParticle]), {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
      mHistManager.add(Form("TrackCuts/TPCSignal/fTPCSignal%s", particleSpecies[iParticle]), Form("%s TPC energy loss;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};#bf{TPC d#it{E}/d#it{x}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});

      mHistManager.add(Form("TrackCuts/%s/fP", particleSpecies[iParticle]), Form("%s momentum at PV;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0, 10}});
      mHistManager.add(Form("TrackCuts/%s/fPt", particleSpecies[iParticle]), Form("%s transverse momentum;#bf{#it{p}_{T}^{%s} (GeV/#it{c})};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0, 10}});
      mHistManager.add(Form("TrackCuts/%s/fMomCorDif", particleSpecies[iParticle]), Form("Momentum correlation;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{p}_{TPC}^{%s} - #it{p}^{%s} (GeV/#it{c})}", particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
      mHistManager.add(Form("TrackCuts/%s/fMomCorRatio", particleSpecies[iParticle]), Form("Relative momentum correlation;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{p}_{TPC}^{%s} - #it{p}^{%s} / #it{p}^{%s}}", particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
      mHistManager.add(Form("TrackCuts/%s/fEta", particleSpecies[iParticle]), Form("%s pseudorapidity distribution;#eta;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -2, 2}});
      mHistManager.add(Form("TrackCuts/%s/fPhi", particleSpecies[iParticle]), Form("%s azimuthal angle distribution;#phi;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCvsTPCP", particleSpecies[iParticle]), Form("NSigmaTPC %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};#bf{n#sigma_{TPC}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTOFvsTPCP", particleSpecies[iParticle]), Form("NSigmaTOF %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};#bf{n#sigma_{TOF}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCTOFvsTPCP", particleSpecies[iParticle]), Form("NSigmaTPCTOF %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};n#sigma_{comb}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaITSvsP", particleSpecies[iParticle]), Form("NSigmaITS %s;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{n#sigma_{ITS}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCvsP", particleSpecies[iParticle]), Form("NSigmaTPC %s P;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{n#sigma_{TPC}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTOFvsP", particleSpecies[iParticle]), Form("NSigmaTOF %s P;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{n#sigma_{TOF}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCTOFvsP", particleSpecies[iParticle]), Form("NSigmaTPCTOF %s P;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{n#sigma_{comb}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

      mHistManager.add(Form("TrackCuts/%s/fDCAxy", particleSpecies[iParticle]), Form("fDCAxy %s;#bf{DCA_{xy}};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -0.5f, 0.5f}});
      mHistManager.add(Form("TrackCuts/%s/fDCAz", particleSpecies[iParticle]), Form("fDCAz %s;#bf{DCA_{z}};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -0.5f, 0.5f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCsCls", particleSpecies[iParticle]), Form("fTPCsCls %s;#bf{TPC Shared Clusters};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCcRows", particleSpecies[iParticle]), Form("fTPCcRows %s;#bf{TPC Crossed Rows};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTrkTPCfCls", particleSpecies[iParticle]), Form("fTrkTPCfCls %s;#bf{TPC Findable/CrossedRows};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0.0f, 3.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCncls", particleSpecies[iParticle]), Form("fTPCncls %s;#bf{TPC Clusters};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
    }

    // --> HNM QA
    // pi+ daughter
    mHistManager.add("HNM/Before/PosDaughter/fInvMass", "Invariant mass HMN Pos Daugh;#bf{#it{M}^{#pi^{+}} (GeV/#it{c}^{2})};#bf{#it{N}^{#pi^{+}}}", HistType::kTH1F, {{200, 0, 0.2}});
    mHistManager.add("HNM/Before/PosDaughter/fPt", "Transverse momentum HMN Pos Daugh tracks;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}^{#pi^{+}}}", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("HNM/Before/PosDaughter/fEta", "HMN Pos Daugh Eta;#eta;#bf{#it{N}^{#pi^{+}}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("HNM/Before/PosDaughter/fPhi", "Azimuthal angle of HMN Pos Daugh tracks;#phi;#bf{#it{N}^{#pi^{+}}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});
    // pi- daughter
    mHistManager.add("HNM/Before/NegDaughter/fInvMass", "Invariant mass HMN Neg Daugh;#bf{#it{M}^{#pi^{-}} (GeV/#it{c}^{2})};#bf{#it{N}^{#pi^{-}}}", HistType::kTH1F, {{200, 0, 0.2}});
    mHistManager.add("HNM/Before/NegDaughter/fPt", "Transverse momentum HMN Neg Daugh tracks;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}^{#pi^{-}}}", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("HNM/Before/NegDaughter/fEta", "HMN Neg Daugh Eta;#eta;#bf{#it{N}^{#pi^{-}}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("HNM/Before/NegDaughter/fPhi", "Azimuthal angle of HMN Neg Daugh tracks;#phi;#bf{#it{N}^{#pi^{-}}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});
    // Properties of the pi+pi- pair
    mHistManager.add("HNM/Before/PiPlPiMi/fInvMassVsPt", "Invariant mass and pT of #pi^+pi^- pairs;#bf{#it{M}^{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}} (GeV/#it{c})}", HistType::kTH2F, {{400, 0.2, 1.}, {250, 0., 25.}});
    mHistManager.add("HNM/Before/PiPlPiMi/fEta", "Pseudorapidity of HMNCand;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("HNM/Before/PiPlPiMi/fPhi", "Azimuthal angle of HMNCand;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

    for (const auto& BeforeAfterString : {"Before", "After"}) {
      for (const auto& iHNM : {"Omega", "EtaPrime"}) {
        for (const auto& MethodString : {"PCM", "EMC"}) {
          mHistManager.add(Form("HNM/%s/%s/%s/fInvMassVsPt", BeforeAfterString, iHNM, MethodString), "Invariant mass and pT of heavy neutral meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
          mHistManager.add(Form("HNM/%s/%s/%s/fEta", BeforeAfterString, iHNM, MethodString), "Pseudorapidity of HNM candidate;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{500, -2, 2}});
          mHistManager.add(Form("HNM/%s/%s/%s/fPhi", BeforeAfterString, iHNM, MethodString), "Azimuthal angle of HNM candidate;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});
        }
      }
    }
    mHistManager.add("HNM/Before/Omega/PCMEMC/fInvMassVsPt", "Invariant mass and pT of omega meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
    mHistManager.add("HNM/Before/Omega/PCMEMC/fEta", "Pseudorapidity of HMNCand;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("HNM/Before/Omega/PCMEMC/fPhi", "Azimuthal angle of HMNCand;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});
    mHistManager.add("HNM/Before/EtaPrime/PCMEMC/fInvMassVsPt", "Invariant mass and pT of eta' meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{600, 0.8, 1.2}, {250, 0., 25.}});
    mHistManager.add("HNM/Before/EtaPrime/PCMEMC/fEta", "Pseudorapidity of HMNCand;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("HNM/Before/EtaPrime/PCMEMC/fPhi", "Azimuthal angle of HMNCand;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

    // --> Two body femto histograms
    for (const auto& iFemtoPartner : {"p", "d"}) {
      for (const auto& iHNM : {"omega", "etaprime"}) {
        mHistManager.add(Form("%s%s/fMultiplicity", iHNM, iFemtoPartner), "Multiplicity of all processed events;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
        mHistManager.add(Form("%s%s/fZvtx", iHNM, iFemtoPartner), "Zvtx of all processed events;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});
        for (const auto& iEMCPCM : {"PCM", "EMC"}) {
          mHistManager.add(Form("%s%s/fSE_particle_%s", iHNM, iFemtoPartner, iEMCPCM), Form("Same Event distribution;#bf{#it{K}^{*} (GeV/#it{c})};#bf{#it{N}^{%s}}", iFemtoPartner), HistType::kTH1F, {{8000, 0, 8}});
          mHistManager.add(Form("%s%s/fSE_Antiparticle_%s", iHNM, iFemtoPartner, iEMCPCM), Form("Same Event distribution;#bf{#it{K}^{*} (GeV/#it{c})};#bf{#it{N}^{#bar{%s}}}", iFemtoPartner), HistType::kTH1F, {{8000, 0, 8}});
          mHistManager.add(Form("%s%s/f%sPtVskstar_%s", iHNM, iFemtoPartner, iHNM, iEMCPCM), Form("K* vs %s pt;#bf{#it{K}^{*} (GeV/#it{c})};#bf{#it{p}_{T}^{%s} (GeV/#it{c})}", iHNM, iHNM), HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
          mHistManager.add(Form("%s%s/f%sPtVskstar_%s", iHNM, iFemtoPartner, iFemtoPartner, iEMCPCM), Form("K* vs %s pt;#bf{#it{K}^{*} (GeV/#it{c})};#bf{#it{p}_{T}^{%s} (GeV/#it{c})}", iFemtoPartner, iFemtoPartner), HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
          mHistManager.add(Form("%s%s/fAnti%sPtVskstar_%s", iHNM, iFemtoPartner, iFemtoPartner, iEMCPCM), Form("K* vs #bar{%s} pt;#bf{#it{K}^{*} (GeV/#it{c})};#bf{#it{p}_{T}^{#bar{%s}} (GeV/#it{c})}", iFemtoPartner, iFemtoPartner), HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
          mHistManager.add(Form("%s%s/fInvMassVsKStar_%s", iHNM, iFemtoPartner, iEMCPCM), "Invariant mass and K* of heavy neutral meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{K}^{*} (GeV/#it{c})}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 1.}});
        }
      }
    }

    // --> Three body femto histograms
    for (const auto& iHNM : {"omega", "etaprime"}) {
      mHistManager.add(Form("pp%s/fMultiplicity", iHNM), "Multiplicity of all processed events;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
      mHistManager.add(Form("pp%s/fZvtx", iHNM), "Zvtx of all processed events;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});
      for (const auto& iEMCPCM : {"PCM", "EMC"}) {
        mHistManager.add(Form("pp%s/fSE_particle_%s", iHNM, iEMCPCM), "Same Event distribution;#bf{#it{Q}_{3} (GeV/#it{c})};#bf{#it{N}^{pp}}", HistType::kTH1F, {{8000, 0, 8}});
        mHistManager.add(Form("pp%s/fSE_Antiparticle_%s", iHNM, iEMCPCM), "Same Event distribution;#bf{#it{Q}_{3} (GeV/#it{c})};#bf{#it{N}^{#bar{p}#bar{p}}}", HistType::kTH1F, {{8000, 0, 8}});
        mHistManager.add(Form("pp%s/fProtonPtVsQ3_%s", iHNM, iEMCPCM), "pT (proton) vs Q_{3};#bf{#it{Q}_{3} (GeV/#it{c})};#bf{#it{p}_{T}^{p} (GeV/#it{c})}", HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
        mHistManager.add(Form("pp%s/f%sCandPtVsQ3_%s", iHNM, iHNM, iEMCPCM), Form("pT (%s) vs Q_{3};#bf{#it{Q}_{3} (GeV/#it{c})};#bf{#it{p}_{T}^{%s} (GeV/#it{c})}", iHNM, iHNM), HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
        mHistManager.add(Form("pp%s/fAntiProtonPtVsQ3_%s", iHNM, iEMCPCM), "pT (antiproton) vs Q_{3};#bf{#it{Q}_{3} (GeV/#it{c})};#bf{#it{p}_{T}^{#bar{p}} (GeV/#it{c})}", HistType::kTH2F, {{{150, 0, 1.5}, {500, 0, 10}}});
        mHistManager.add(Form("pp%s/fInvMassVsQ3_%s", iHNM, iEMCPCM), "Invariant mass and Q3 of heavy neutral meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{Q}_{3} (GeV/#it{c})}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 1.}});
      }
    }

    if (cfgDoEMCShift.value) {
      for (int iSM = 0; iSM < nSMs; iSM++) {
        emcEtaShift[iSM] = cfgEMCEtaShift.value[iSM];
        emcPhiShift[iSM] = cfgEMCPhiShift.value[iSM];
        LOG(info) << "SM-wise shift in eta/phi for SM " << iSM << ": " << emcEtaShift[iSM] << " / " << emcPhiShift[iSM];
      }
    }
  }

  void process(aod::MyCollision const& collision, aod::MyBCs const&, aod::SkimEMCClusters const& clusters, aod::V0PhotonsKF const& v0s, aod::SelectedTracks const& tracks)
  {
    // inlcude ITS PID information
    auto tracksWithItsPid = soa::Attach<aod::SelectedTracks, aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe>(tracks);

    // QA all evts
    mHistManager.fill(HIST("fProcessedEvents"), 0);
    mHistManager.fill(HIST("Event/fMultiplicityBefore"), collision.multNTracksPV());
    mHistManager.fill(HIST("Event/fZvtxBefore"), collision.posZ());

    // Ensure evts are consistent with Sel8 and Vtx-z selection
    bool keepFemtoEvent[hnmtrigger::kNFemtoTriggers] = {false, false, false, false, false, false}; // Set based on number of found pairs (see above) - used to flag femto events
    if (!isSelectedEvent(collision)) {
      tags(keepFemtoEvent[hnmtrigger::kOmegaP], keepFemtoEvent[hnmtrigger::kPPOmega], keepFemtoEvent[hnmtrigger::kOmegaD], keepFemtoEvent[hnmtrigger::kEtaPrimeP], keepFemtoEvent[hnmtrigger::kPPEtaPrime], keepFemtoEvent[hnmtrigger::kEtaPrimeD]);
      return;
    }
    // QA accepted evts
    mHistManager.fill(HIST("Event/fMultiplicityAfter"), collision.multNTracksPV());
    mHistManager.fill(HIST("Event/fZvtxAfter"), collision.posZ());

    colContainsPCMOmega = colContainsEMCOmega = colContainsPCMEtaPrime = colContainsEMCEtaPrime = false; // Used by spectrum trigger to flag events with high-pT omega/eta' candidates
    int lowMomentumMultiplets[hnmtrigger::kNFemtoTriggers] = {0, 0, 0, 0, 0, 0};                         // Number of found femto pairs/triplets for each femto trigger

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

    // ---------------------------------> EMCal event QA <----------------------------------
    // - Fill Event/nEMCalEvents histogram for EMCal event QA
    // -------------------------------------------------------------------------------------
    bool bcHasEMCCells = collision.isemcreadout();
    bool iskTVXinEMC = collision.foundBC_as<aod::MyBCs>().alias_bit(kTVXinEMC);
    bool isL0Triggered = collision.foundBC_as<aod::MyBCs>().alias_bit(kEMC7) || collision.foundBC_as<aod::MyBCs>().alias_bit(kEG1) || collision.foundBC_as<aod::MyBCs>().alias_bit(kEG2);

    if (bcHasEMCCells && iskTVXinEMC)
      mHistManager.fill(HIST("Event/nEMCalEvents"), 0);
    if (bcHasEMCCells && isL0Triggered)
      mHistManager.fill(HIST("Event/nEMCalEvents"), 1);
    if (bcHasEMCCells && !iskTVXinEMC && !isL0Triggered)
      mHistManager.fill(HIST("Event/nEMCalEvents"), 2);
    if (!bcHasEMCCells && iskTVXinEMC)
      mHistManager.fill(HIST("Event/nEMCalEvents"), 3);
    if (!bcHasEMCCells && isL0Triggered)
      mHistManager.fill(HIST("Event/nEMCalEvents"), 4);

    // --------------------------------> Process Photons <----------------------------------
    // - Slice clusters and V0s by collision ID to get the ones in this collision
    // - Store the clusters and V0s in the vGammas vector
    // - Reconstruct gamma-gamma pairs
    // -------------------------------------------------------------------------------------
    auto v0sInThisCollision = v0s.sliceBy(perCollisionPCM, collision.globalIndex());
    auto clustersInThisCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());
    mHistManager.fill(HIST("Event/nClustersVsV0s"), clustersInThisCollision.size(), v0sInThisCollision.size());

    std::vector<hnmutilities::Photon> vGammas;
    hnmutilities::storeGammasInVector(clustersInThisCollision, v0sInThisCollision, vGammas, emcEtaShift, emcPhiShift);
    hnmutilities::reconstructGGs(vGammas, vGGs);
    vGammas.clear();
    processGGs(vGGs);

    // ------------------------------> Loop over all tracks <-------------------------------
    // - Sort them into vectors based on PID ((anti)protons, (anti)deuterons, (anti)pions)
    // - Fill QA histograms for all tracks and per particle species
    // -------------------------------------------------------------------------------------
    for (const auto& track : tracksWithItsPid) {
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fPtTrackBefore"), track.pt());
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fEtaTrackBefore"), track.eta());
      mHistManager.fill(HIST("TrackCuts/TracksBefore/fPhiTrackBefore"), track.phi());
      if (track.sign() > 0) { // All particles (positive electric charge)
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalTPCP"), track.tpcInnerParam(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignal"), track.p(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationPos"), track.p(), track.tpcInnerParam());
      }
      if (track.sign() < 0) { // All anti-particles (negative electric charge)
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiTPCP"), track.tpcInnerParam(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAnti"), track.p(), track.tpcSignal());
        mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationNeg"), track.p(), track.tpcInnerParam());
      }

      // For each track, check if it fulfills track and PID criteria to be identified as a proton, deuteron or pion
      bool isProton = (isSelectedTrackPID(track, hnmtrigger::kProton) && isSelectedTrack(track, hnmtrigger::kProton));
      bool isDeuteron = (isSelectedTrackPID(track, hnmtrigger::kDeuteron) && isSelectedTrack(track, hnmtrigger::kDeuteron));
      bool isPion = (isSelectedTrackPID(track, hnmtrigger::kPion) && isSelectedTrack(track, hnmtrigger::kPion));

      if (track.sign() > 0) { // Positive charge -> Particles
        if (isProton) {
          proton.emplace_back(track.pt(), track.eta(), track.phi(), mMassProton);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsProton"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalProton"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/Proton/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Proton/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Proton/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Proton/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaPr());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2));
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaITSvsP"), track.p(), track.itsNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTOFvsP"), track.p(), track.tofNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/Proton/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Proton/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Proton/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Proton/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Proton/fTPCncls"), track.tpcNClsFound());
        }
        if (isDeuteron) {
          deuteron.emplace_back(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsDeuteron"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalDeuteron"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/Deuteron/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaDe());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaDe(), 2) + std::pow(track.tofNSigmaDe(), 2));
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaITSvsP"), track.p(), track.itsNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTOFvsP"), track.p(), track.tofNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Deuteron/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Deuteron/fTPCncls"), track.tpcNClsFound());
        }
        if (isPion) {
          pion.emplace_back(track.pt(), track.eta(), track.phi(), mMassPionCharged);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsPion"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalPion"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/Pion/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/Pion/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/Pion/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/Pion/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaPi());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2));
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaITSvsP"), track.p(), track.itsNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTOFvsP"), track.p(), track.tofNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/Pion/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/Pion/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/Pion/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/Pion/fTPCncls"), track.tpcNClsFound());
        }
      } else { // Negative charge -> Anti-particles
        if (isProton) {
          antiproton.emplace_back(track.pt(), track.eta(), track.phi(), mMassProton);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiProton"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiProton"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/AntiProton/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaPr());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPr(), 2) + std::pow(track.tofNSigmaPr(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaITSvsP"), track.p(), track.itsNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTOFvsP"), track.p(), track.tofNSigmaPr());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPr() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPr() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiProton/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiProton/fTPCncls"), track.tpcNClsFound());
        }
        if (isDeuteron) {
          antideuteron.emplace_back(track.pt(), track.eta(), track.phi(), mMassDeuteron);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiDeuteron"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiDeuteron"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaDe());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaDe(), 2) + std::pow(track.tofNSigmaDe(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaITSvsP"), track.p(), track.itsNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTOFvsP"), track.p(), track.tofNSigmaDe());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaDe() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaDe() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiDeuteron/fTPCncls"), track.tpcNClsFound());
        }
        if (isPion) {
          antipion.emplace_back(track.pt(), track.eta(), track.phi(), mMassPionCharged);
          mHistManager.fill(HIST("TrackCuts/TracksBefore/fMomCorrelationAfterCutsAntiPion"), track.p(), track.tpcInnerParam());
          mHistManager.fill(HIST("TrackCuts/TPCSignal/fTPCSignalAntiPion"), track.tpcInnerParam(), track.tpcSignal());

          mHistManager.fill(HIST("TrackCuts/AntiPion/fP"), track.p());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fPt"), track.pt());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fMomCorDif"), track.p(), track.tpcInnerParam() - track.p());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fMomCorRatio"), track.p(), (track.tpcInnerParam() - track.p()) / track.p());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fEta"), track.eta());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fPhi"), track.phi());

          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCvsTPCP"), track.tpcInnerParam(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTOFvsTPCP"), track.tpcInnerParam(), track.tofNSigmaPi());
          auto nSigmaTrackTPCTOF = std::sqrt(std::pow(track.tpcNSigmaPi(), 2) + std::pow(track.tofNSigmaPi(), 2));
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCTOFvsTPCP"), track.tpcInnerParam(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaITSvsP"), track.p(), track.itsNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTOFvsP"), track.p(), track.tofNSigmaPi());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCTOFvsP"), track.p(), std::sqrt(std::pow(track.tpcNSigmaPi() - nSigmaTrackTPCTOF, 2) + std::pow(track.tofNSigmaPi() - nSigmaTrackTPCTOF, 2)));

          mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAxy"), track.dcaXY());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAz"), track.dcaZ());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCsCls"), track.tpcNClsShared());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCcRows"), track.tpcNClsCrossedRows());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
          mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCncls"), track.tpcNClsFound());
        }
      }
    }

    // -------------------------> Reconstruct HNM candidates <------------------------------
    // - Based on the previously filled (anti)pion vectors
    // - Fill QA histograms for kinematics of the pions and their combinations
    // -------------------------------------------------------------------------------------
    for (const auto& posPion : pion) {
      for (const auto& negPion : antipion) {
        ROOT::Math::PtEtaPhiMVector vecPiPlPiMi = posPion + negPion;
        hnmutilities::reconstructHeavyNeutralMesons(vecPiPlPiMi, vGGs, vHNMs);

        mHistManager.fill(HIST("HNM/Before/PiPlPiMi/fInvMassVsPt"), vecPiPlPiMi.M(), vecPiPlPiMi.pt());
        mHistManager.fill(HIST("HNM/Before/PiPlPiMi/fEta"), vecPiPlPiMi.eta());
        mHistManager.fill(HIST("HNM/Before/PiPlPiMi/fPhi"), RecoDecay::constrainAngle(vecPiPlPiMi.phi()));

        mHistManager.fill(HIST("HNM/Before/PosDaughter/fInvMass"), posPion.M());
        mHistManager.fill(HIST("HNM/Before/PosDaughter/fPt"), posPion.pt());
        mHistManager.fill(HIST("HNM/Before/PosDaughter/fEta"), posPion.eta());
        mHistManager.fill(HIST("HNM/Before/PosDaughter/fPhi"), RecoDecay::constrainAngle(posPion.phi()));

        mHistManager.fill(HIST("HNM/Before/NegDaughter/fInvMass"), negPion.M());
        mHistManager.fill(HIST("HNM/Before/NegDaughter/fPt"), negPion.pt());
        mHistManager.fill(HIST("HNM/Before/NegDaughter/fEta"), negPion.eta());
        mHistManager.fill(HIST("HNM/Before/NegDaughter/fPhi"), RecoDecay::constrainAngle(negPion.phi()));
      }
    }

    // ---------------------------> Process HNM candidates <--------------------------------
    // - Fill invMassVsPt histograms separated into HNM types (based on GG mass) and gamma reco method
    // - Set colContains* flags for each HNM type to be used in the high-pt spectrum trigger
    // - Fill femto HNM vectors (omegaPCM, etaPrimePCM, omegaEMC, etaPrimeEMC)
    // -------------------------------------------------------------------------------------
    processHNMs(vHNMs);

    // ------------------------------> Build triplets <-------------------------------------
    // - Calculate Q3 for each triplet (p-p-omega, p-p-eta', anti-p-anti-p-omega, anti-p-anti-p-eta')
    // - Fill QA histograms for Q3 and pT of the triplet and its daughters
    // - Increment lowMomentumMultiplets for each triplet with Q3 < kinematic limit (used in femto trigger)
    // -------------------------------------------------------------------------------------
    if (cfgTriggerSwitches->get("Switch", "PPOmega") > 0.) { // -----> p-p-omega femtoscopy
      for (size_t i = 0; i < proton.size(); ++i) {
        for (size_t j = i + 1; j < proton.size(); ++j) {
          const auto& proton1 = proton[i];
          const auto& proton2 = proton[j];
          for (const auto& omegaParticles : omegaPCM) { // ---> PCM

            float q3 = getQ3(proton1, proton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_particle_PCM"), q3);
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_PCM"), q3, proton1.Pt());
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_PCM"), q3, proton2.Pt());
            mHistManager.fill(HIST("ppomega/fomegaCandPtVsQ3_PCM"), q3, omegaParticles.Pt());
            mHistManager.fill(HIST("ppomega/fInvMassVsQ3_PCM"), omegaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPOmega))
              lowMomentumMultiplets[hnmtrigger::kPPOmega] += 1;
          }
          for (const auto& omegaParticles : omegaEMC) { // ---> EMC

            float q3 = getQ3(proton1, proton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_particle_EMC"), q3);
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_EMC"), q3, proton1.Pt());
            mHistManager.fill(HIST("ppomega/fProtonPtVsQ3_EMC"), q3, proton2.Pt());
            mHistManager.fill(HIST("ppomega/fomegaCandPtVsQ3_EMC"), q3, omegaParticles.Pt());
            mHistManager.fill(HIST("ppomega/fInvMassVsQ3_EMC"), omegaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPOmega))
              lowMomentumMultiplets[hnmtrigger::kPPOmega] += 1;
          }
        }
      }
      for (size_t i = 0; i < antiproton.size(); ++i) { // -----> antip-antip-omega femtoscopy
        for (size_t j = i + 1; j < antiproton.size(); ++j) {
          const auto& antiProton1 = antiproton[i];
          const auto& antiProton2 = antiproton[j];
          for (const auto& omegaParticles : omegaPCM) { // ---> PCM

            float q3 = getQ3(antiProton1, antiProton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_Antiparticle_PCM"), q3);
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_PCM"), q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_PCM"), q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppomega/fomegaCandPtVsQ3_PCM"), q3, omegaParticles.Pt());
            mHistManager.fill(HIST("ppomega/fInvMassVsQ3_PCM"), omegaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPOmega))
              lowMomentumMultiplets[hnmtrigger::kPPOmega] += 1;
          }
          for (const auto& omegaParticles : omegaEMC) { // ---> EMC

            float q3 = getQ3(antiProton1, antiProton2, omegaParticles);

            mHistManager.fill(HIST("ppomega/fSE_Antiparticle_EMC"), q3);
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_EMC"), q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppomega/fAntiProtonPtVsQ3_EMC"), q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppomega/fomegaCandPtVsQ3_EMC"), q3, omegaParticles.Pt());
            mHistManager.fill(HIST("ppomega/fInvMassVsQ3_EMC"), omegaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPOmega))
              lowMomentumMultiplets[hnmtrigger::kPPOmega] += 1;
          }
        }
      }
    }
    if (cfgTriggerSwitches->get("Switch", "PPEtaPrime") > 0.) { // -----> p-p-eta' femtoscopy
      for (size_t i = 0; i < proton.size(); ++i) {
        for (size_t j = i + 1; j < proton.size(); ++j) {
          const auto& proton1 = proton[i];
          const auto& proton2 = proton[j];
          for (const auto& etaParticles : etaPrimePCM) { // ---> PCM

            float q3 = getQ3(proton1, proton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_particle_PCM"), q3);
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_PCM"), q3, proton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_PCM"), q3, proton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fetaprimeCandPtVsQ3_PCM"), q3, etaParticles.Pt());
            mHistManager.fill(HIST("ppetaprime/fInvMassVsQ3_PCM"), etaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPEtaPrime))
              lowMomentumMultiplets[hnmtrigger::kPPEtaPrime] += 1;
          }
          for (const auto& etaParticles : etaPrimeEMC) { // ---> EMC

            float q3 = getQ3(proton1, proton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_particle_EMC"), q3);
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_EMC"), q3, proton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fProtonPtVsQ3_EMC"), q3, proton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fetaprimeCandPtVsQ3_EMC"), q3, etaParticles.Pt());
            mHistManager.fill(HIST("ppetaprime/fInvMassVsQ3_EMC"), etaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPEtaPrime))
              lowMomentumMultiplets[hnmtrigger::kPPEtaPrime] += 1;
          }
        }
      }
      for (size_t i = 0; i < antiproton.size(); ++i) { // -----> antip-antip-eta' femtoscopy
        for (size_t j = i + 1; j < antiproton.size(); ++j) {
          const auto& antiProton1 = antiproton[i];
          const auto& antiProton2 = antiproton[j];
          for (const auto& etaParticles : etaPrimePCM) { // ---> PCM

            float q3 = getQ3(antiProton1, antiProton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_Antiparticle_PCM"), q3);
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_PCM"), q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_PCM"), q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fetaprimeCandPtVsQ3_PCM"), q3, etaParticles.Pt());
            mHistManager.fill(HIST("ppetaprime/fInvMassVsQ3_PCM"), etaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPEtaPrime))
              lowMomentumMultiplets[hnmtrigger::kPPEtaPrime] += 1;
          }
          for (const auto& etaParticles : etaPrimeEMC) { // ---> EMC

            float q3 = getQ3(antiProton1, antiProton2, etaParticles);

            mHistManager.fill(HIST("ppetaprime/fSE_Antiparticle_EMC"), q3);
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_EMC"), q3, antiProton1.Pt());
            mHistManager.fill(HIST("ppetaprime/fAntiProtonPtVsQ3_EMC"), q3, antiProton2.Pt());
            mHistManager.fill(HIST("ppetaprime/fetaprimeCandPtVsQ3_EMC"), q3, etaParticles.Pt());
            mHistManager.fill(HIST("ppetaprime/fInvMassVsQ3_EMC"), etaParticles.M(), q3);

            if (q3 < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kPPEtaPrime))
              lowMomentumMultiplets[hnmtrigger::kPPEtaPrime] += 1;
          }
        }
      }
    }

    // --------------------------------> Build Pairs <--------------------------------------
    // - Calculate k* for each pair ((anti)d-omega, (anti)d-eta', (anti)p-omega, (anti)p-eta')
    // - Fill QA histograms for k* and pT of the pairs
    // - Increment lowMomentumMultiplets for each triplet with k* < kinematic limit (used in femto trigger)
    // -------------------------------------------------------------------------------------
    if (cfgTriggerSwitches->get("Switch", "Omegad") > 0.) {
      for (auto iomega = omegaPCM.begin(); iomega != omegaPCM.end(); ++iomega) {            // -----> PCM
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) { // ---> d-omega femtoscopy

          float kstar = getkstar(*iomega, *iDeuteron);

          mHistManager.fill(HIST("omegad/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fdPtVskstar_PCM"), kstar, (*iDeuteron).Pt());
          mHistManager.fill(HIST("omegad/fInvMassVsKStar_PCM"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaD))
            lowMomentumMultiplets[hnmtrigger::kOmegaD] += 1;
        }
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) { // ---> antid-omega femtoscopy

          float kstar = getkstar(*iomega, *iAntiDeuteron);

          mHistManager.fill(HIST("omegad/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fAntidPtVskstar_PCM"), kstar, (*iAntiDeuteron).Pt());
          mHistManager.fill(HIST("omegad/fInvMassVsKStar_PCM"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaD))
            lowMomentumMultiplets[hnmtrigger::kOmegaD] += 1;
        }
      }
      for (auto iomega = omegaEMC.begin(); iomega != omegaEMC.end(); ++iomega) {            // -----> EMC
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) { // ---> d-omega femtoscopy

          float kstar = getkstar(*iomega, *iDeuteron);

          mHistManager.fill(HIST("omegad/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fdPtVskstar_EMC"), kstar, (*iDeuteron).Pt());
          mHistManager.fill(HIST("omegad/fInvMassVsKStar_EMC"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaD))
            lowMomentumMultiplets[hnmtrigger::kOmegaD] += 1;
        }
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) { // ---> antid-omega femtoscopy

          float kstar = getkstar(*iomega, *iAntiDeuteron);

          mHistManager.fill(HIST("omegad/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("omegad/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegad/fAntidPtVskstar_EMC"), kstar, (*iAntiDeuteron).Pt());
          mHistManager.fill(HIST("omegad/fInvMassVsKStar_EMC"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaD))
            lowMomentumMultiplets[hnmtrigger::kOmegaD] += 1;
        }
      }
    }
    if (cfgTriggerSwitches->get("Switch", "EtaPrimed") > 0.) {
      for (auto ietaprime = etaPrimePCM.begin(); ietaprime != etaPrimePCM.end(); ++ietaprime) { // -----> PCM
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {     // ---> d-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iDeuteron);

          mHistManager.fill(HIST("etaprimed/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fdPtVskstar_PCM"), kstar, (*iDeuteron).Pt());
          mHistManager.fill(HIST("etaprimed/fInvMassVsKStar_PCM"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeD))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeD] += 1;
        }
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) { // ---> antid-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iAntiDeuteron);

          mHistManager.fill(HIST("etaprimed/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fAntidPtVskstar_PCM"), kstar, (*iAntiDeuteron).Pt());
          mHistManager.fill(HIST("etaprimed/fInvMassVsKStar_PCM"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeD))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeD] += 1;
        }
      }
      for (auto ietaprime = etaPrimeEMC.begin(); ietaprime != etaPrimeEMC.end(); ++ietaprime) { // -----> EMC
        for (auto iDeuteron = deuteron.begin(); iDeuteron != deuteron.end(); ++iDeuteron) {     // ---> d-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iDeuteron);

          mHistManager.fill(HIST("etaprimed/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fdPtVskstar_EMC"), kstar, (*iDeuteron).Pt());
          mHistManager.fill(HIST("etaprimed/fInvMassVsKStar_EMC"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeD))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeD] += 1;
        }
        for (auto iAntiDeuteron = antideuteron.begin(); iAntiDeuteron != antideuteron.end(); ++iAntiDeuteron) { // ---> antid-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iAntiDeuteron);

          mHistManager.fill(HIST("etaprimed/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimed/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimed/fAntidPtVskstar_EMC"), kstar, (*iAntiDeuteron).Pt());
          mHistManager.fill(HIST("etaprimed/fInvMassVsKStar_EMC"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeD))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeD] += 1;
        }
      }
    }
    if (cfgTriggerSwitches->get("Switch", "OmegaP") > 0.) {
      for (auto iomega = omegaPCM.begin(); iomega != omegaPCM.end(); ++iomega) {  // -----> PCM
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) { // ---> p-omega femtoscopy

          float kstar = getkstar(*iomega, *iProton);

          mHistManager.fill(HIST("omegap/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fpPtVskstar_PCM"), kstar, (*iProton).Pt());
          mHistManager.fill(HIST("omegap/fInvMassVsKStar_PCM"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaP))
            lowMomentumMultiplets[hnmtrigger::kOmegaP] += 1;
        }
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) { // ---> antip-omega femtoscopy

          float kstar = getkstar(*iomega, *iAntiProton);

          mHistManager.fill(HIST("omegap/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_PCM"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fAntipPtVskstar_PCM"), kstar, (*iAntiProton).Pt());
          mHistManager.fill(HIST("omegap/fInvMassVsKStar_PCM"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaP))
            lowMomentumMultiplets[hnmtrigger::kOmegaP] += 1;
        }
      }
      for (auto iomega = omegaEMC.begin(); iomega != omegaEMC.end(); ++iomega) {  // -----> EMC
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) { // ---> p-omega femtoscopy

          float kstar = getkstar(*iomega, *iProton);

          mHistManager.fill(HIST("omegap/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fpPtVskstar_EMC"), kstar, (*iProton).Pt());
          mHistManager.fill(HIST("omegap/fInvMassVsKStar_EMC"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaP))
            lowMomentumMultiplets[hnmtrigger::kOmegaP] += 1;
        }
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) { // ---> antip-omega femtoscopy

          float kstar = getkstar(*iomega, *iAntiProton);

          mHistManager.fill(HIST("omegap/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("omegap/fomegaPtVskstar_EMC"), kstar, (*iomega).Pt());
          mHistManager.fill(HIST("omegap/fAntipPtVskstar_EMC"), kstar, (*iAntiProton).Pt());
          mHistManager.fill(HIST("omegap/fInvMassVsKStar_EMC"), (*iomega).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kOmegaP))
            lowMomentumMultiplets[hnmtrigger::kOmegaP] += 1;
        }
      }
    }
    if (cfgTriggerSwitches->get("Switch", "EtaPrimeP") > 0.) {
      for (auto ietaprime = etaPrimePCM.begin(); ietaprime != etaPrimePCM.end(); ++ietaprime) { // -----> PCM
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {               // ---> p-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iProton);

          mHistManager.fill(HIST("etaprimep/fSE_particle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fpPtVskstar_PCM"), kstar, (*iProton).Pt());
          mHistManager.fill(HIST("etaprimep/fInvMassVsKStar_PCM"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeP))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeP] += 1;
        }
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) { // ---> antip-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iAntiProton);

          mHistManager.fill(HIST("etaprimep/fSE_Antiparticle_PCM"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_PCM"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fAntipPtVskstar_PCM"), kstar, (*iAntiProton).Pt());
          mHistManager.fill(HIST("etaprimep/fInvMassVsKStar_PCM"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeP))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeP] += 1;
        }
      }
      for (auto ietaprime = etaPrimeEMC.begin(); ietaprime != etaPrimeEMC.end(); ++ietaprime) { // -----> EMC
        for (auto iProton = proton.begin(); iProton != proton.end(); ++iProton) {               // ---> p-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iProton);

          mHistManager.fill(HIST("etaprimep/fSE_particle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fpPtVskstar_EMC"), kstar, (*iProton).Pt());
          mHistManager.fill(HIST("etaprimep/fInvMassVsKStar_EMC"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeP))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeP] += 1;
        }
        for (auto iAntiProton = antiproton.begin(); iAntiProton != antiproton.end(); ++iAntiProton) { // ---> antip-eta' femtoscopy

          float kstar = getkstar(*ietaprime, *iAntiProton);

          mHistManager.fill(HIST("etaprimep/fSE_Antiparticle_EMC"), kstar);
          mHistManager.fill(HIST("etaprimep/fetaprimePtVskstar_EMC"), kstar, (*ietaprime).Pt());
          mHistManager.fill(HIST("etaprimep/fAntipPtVskstar_EMC"), kstar, (*iAntiProton).Pt());
          mHistManager.fill(HIST("etaprimep/fInvMassVsKStar_EMC"), (*ietaprime).M(), kstar);

          if (kstar < cfgKinematicLimits->get(static_cast<uint>(0), hnmtrigger::kEtaPrimeP))
            lowMomentumMultiplets[hnmtrigger::kEtaPrimeP] += 1;
        }
      }
    }

    // -----------------------------> Create femto tags <-----------------------------------
    // - Set keepFemtoEvent flags for each HNM type based on the lowMomentumMultiplets
    // - Fill histograms for the multiplicity and z-vertex of femto-accepted events
    // -------------------------------------------------------------------------------------
    if (lowMomentumMultiplets[hnmtrigger::kPPOmega] > 0) {
      keepFemtoEvent[hnmtrigger::kPPOmega] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 6);
      mHistManager.fill(HIST("ppomega/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("ppomega/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[hnmtrigger::kPPEtaPrime] > 0) {
      keepFemtoEvent[hnmtrigger::kPPEtaPrime] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 7);
      mHistManager.fill(HIST("ppetaprime/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("ppetaprime/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[hnmtrigger::kOmegaD] > 0) {
      keepFemtoEvent[hnmtrigger::kOmegaD] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 8);
      mHistManager.fill(HIST("omegad/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("omegad/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[hnmtrigger::kEtaPrimeD] > 0) {
      keepFemtoEvent[hnmtrigger::kEtaPrimeD] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 9);
      mHistManager.fill(HIST("etaprimed/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("etaprimed/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[hnmtrigger::kOmegaP] > 0) {
      keepFemtoEvent[hnmtrigger::kOmegaP] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 10);
      mHistManager.fill(HIST("omegap/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("omegap/fZvtx"), collision.posZ());
    }
    if (lowMomentumMultiplets[hnmtrigger::kEtaPrimeP] > 0) {
      keepFemtoEvent[hnmtrigger::kEtaPrimeP] = true;
      mHistManager.fill(HIST("fProcessedEvents"), 11);
      mHistManager.fill(HIST("etaprimep/fMultiplicity"), collision.multNTracksPV());
      mHistManager.fill(HIST("etaprimep/fZvtx"), collision.posZ());
    }

    // -----------------------------> Set trigger flags <-----------------------------------
    // - 4 high pT spectrum trigger flags (PCM & EMC * omega & eta')
    // - 4 femto trigger flags (p-omega, p-eta', d-omega || pp-omega, d-eta' || pp-eta')
    // -------------------------------------------------------------------------------------
    tags(keepFemtoEvent[hnmtrigger::kOmegaP], keepFemtoEvent[hnmtrigger::kPPOmega], keepFemtoEvent[hnmtrigger::kOmegaD], keepFemtoEvent[hnmtrigger::kEtaPrimeP], keepFemtoEvent[hnmtrigger::kPPEtaPrime], keepFemtoEvent[hnmtrigger::kEtaPrimeD]);
    // tags(colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime, keepFemtoEvent[hnmtrigger::kOmegaP], keepFemtoEvent[hnmtrigger::kEtaPrimeP],
    //      keepFemtoEvent[hnmtrigger::kPPOmega] || keepFemtoEvent[hnmtrigger::kOmegaD], keepFemtoEvent[hnmtrigger::kPPEtaPrime] || keepFemtoEvent[hnmtrigger::kEtaPrimeD]);

    if (!keepFemtoEvent[hnmtrigger::kPPOmega] && !keepFemtoEvent[hnmtrigger::kOmegaP] && !keepFemtoEvent[hnmtrigger::kPPEtaPrime] && !keepFemtoEvent[hnmtrigger::kEtaPrimeP] && !keepFemtoEvent[hnmtrigger::kOmegaD] && !keepFemtoEvent[hnmtrigger::kEtaPrimeD])
      mHistManager.fill(HIST("fProcessedEvents"), 1); // Fill "rejected", if no trigger selected the event
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

      if (lightMeson->m() > cfgMassWindowOmega->get("pi0_min") && lightMeson->m() < cfgMassWindowOmega->get("pi0_max")) {
        lightMeson->isPi0 = true;
      } else if (lightMeson->m() > cfgMassWindowEtaPrime->get("eta_min") && lightMeson->m() < cfgMassWindowEtaPrime->get("eta_max")) {
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
        if (heavyNeutralMeson.gg->isPi0) {
          mHistManager.fill(HIST("HNM/Before/Omega/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/Omega/PCM/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/Omega/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        } else if (heavyNeutralMeson.gg->isEta) {
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCM/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        }
      } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
        if (heavyNeutralMeson.gg->isPi0) {
          mHistManager.fill(HIST("HNM/Before/Omega/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/Omega/EMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/Omega/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        } else if (heavyNeutralMeson.gg->isEta) {
          mHistManager.fill(HIST("HNM/Before/EtaPrime/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/EMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        }
      } else {
        if (heavyNeutralMeson.gg->isPi0) {
          mHistManager.fill(HIST("HNM/Before/Omega/PCMEMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/Omega/PCMEMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/Omega/PCMEMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        } else if (heavyNeutralMeson.gg->isEta) {
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCMEMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCMEMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/Before/EtaPrime/PCMEMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        }
      }

      if (heavyNeutralMeson.gg->isPi0 && massHNM > cfgMassWindowOmega->get("omega_min") && massHNM < cfgMassWindowOmega->get("omega_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsFemtoTrigger->get("PCM_omega")) {
            omegaPCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), massHNM);
            mHistManager.fill(HIST("HNM/After/Omega/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
            mHistManager.fill(HIST("HNM/After/Omega/PCM/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("HNM/After/Omega/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsSpectrumTrigger->get("PCM_omega"))
            colContainsPCMOmega = true;
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsFemtoTrigger->get("EMC_omega")) {
            omegaEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), massHNM);
            mHistManager.fill(HIST("HNM/After/Omega/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
            mHistManager.fill(HIST("HNM/After/Omega/EMC/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("HNM/After/Omega/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsSpectrumTrigger->get("EMC_omega"))
            colContainsEMCOmega = true;
        }
      } else if (heavyNeutralMeson.gg->isEta && massHNM > cfgMassWindowEtaPrime->get("etaprime_min") && massHNM < cfgMassWindowEtaPrime->get("etaprime_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsFemtoTrigger->get("PCM_etaprime")) {
            etaPrimePCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), massHNM);
            mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
            mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsSpectrumTrigger->get("PCM_etaprime"))
            colContainsPCMEtaPrime = true;
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsFemtoTrigger->get("EMC_etaprime")) {
            etaPrimeEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), massHNM);
            mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
            mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fEta"), heavyNeutralMeson.eta());
            mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
          }
          if (heavyNeutralMeson.pT() > cfgMinHNMPtsSpectrumTrigger->get("EMC_etaprime"))
            colContainsEMCEtaPrime = true;
        }
      } else {
        vHNMs.erase(vHNMs.begin() + iHNM);
        iHNM--;
      }
    }
    mHistManager.fill(HIST("Event/nHeavyNeutralMesons"), nHNMsBeforeMassCuts, vHNMs.size());

    if (colContainsPCMOmega)
      mHistManager.fill(HIST("fProcessedEvents"), 2);
    if (colContainsEMCOmega)
      mHistManager.fill(HIST("fProcessedEvents"), 3);
    if (colContainsPCMEtaPrime)
      mHistManager.fill(HIST("fProcessedEvents"), 4);
    if (colContainsEMCEtaPrime)
      mHistManager.fill(HIST("fProcessedEvents"), 5);
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMesonFilter>(cfgc)}; }
