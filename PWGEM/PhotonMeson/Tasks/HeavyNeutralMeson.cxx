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

#include "PWGEM/PhotonMeson/Utils/HNMUtilities.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector4D.h"
#include "TMath.h"
#include "TRandom3.h"

#include "fairlogger/Logger.h"

#include <iterator>
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

namespace hnm
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

const std::vector<std::string> speciesName{"pion"}; // ToDo include charged pions
const std::vector<std::string> pTCutsName{"Pt min", "Pt max", "P TOF thres"};
const std::vector<std::string> pidCutsName{"TPC min", "TPC max", "TPCTOF max", "ITS min", "ITS max"};
const float pidcutsTable[kNTracksPID][kNPIDLimits]{{-4.f, 4.f, 4.f, -99.f, 99.f}};
const float ptcutsTable[kNTracksPID][3]{{0.35f, 6.f, 0.75f}};
const float nClusterMinTPC[1][kNTracksPID]{{80.0f}};
const float nClusterMinITS[1][kNTracksPID]{{4}};

} // namespace hnm

struct HeavyNeutralMeson {

  // --------------------------------> Configurables <------------------------------------
  // - Event selection cuts
  // - Track selection cuts
  // - Cluster shifts
  // - HNM mass selection windows
  // -------------------------------------------------------------------------------------
  // ---> Event selection
  Configurable<bool> confEvtSelectZvtx{"confEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtRequireSel8{"confEvtRequireSel8", false, "Evt sel: check for offline selection (sel8)"};

  // ---> Track selection
  Configurable<LabeledArray<float>> cfgPtCuts{"cfgPtCuts", {hnm::ptcutsTable[0], 1, 3, hnm::speciesName, hnm::pTCutsName}, "Track pT selections"};
  Configurable<float> cfgTrkEta{"cfgTrkEta", 0.9, "Eta"};
  Configurable<LabeledArray<float>> cfgTPCNClustersMin{"cfgTPCNClustersMin", {hnm::nClusterMinTPC[0], 1, 1, std::vector<std::string>{"TPCNClusMin"}, hnm::speciesName}, "Mininum of TPC Clusters"};
  Configurable<float> cfgTrkTPCfCls{"cfgTrkTPCfCls", 0.83, "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> cfgTrkTPCcRowsMin{"cfgTrkTPCcRowsMin", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgTrkTPCsClsSharedFrac{"cfgTrkTPCsClsSharedFrac", 1.f, "Fraction of shared TPC clusters"};
  Configurable<LabeledArray<float>> cfgTrkITSnclsMin{"cfgTrkITSnclsMin", {hnm::nClusterMinITS[0], 1, 1, std::vector<std::string>{"Cut"}, hnm::speciesName}, "Minimum number of ITS clusters"};
  Configurable<float> cfgTrkDCAxyMax{"cfgTrkDCAxyMax", 0.15, "Maximum DCA_xy"};
  Configurable<float> cfgTrkDCAzMax{"cfgTrkDCAzMax", 0.3, "Maximum DCA_z"};
  Configurable<float> cfgTrkMaxChi2PerClusterTPC{"cfgTrkMaxChi2PerClusterTPC", 4.0f, "Minimal track selection: max allowed chi2 per TPC cluster"};  // 4.0 is default of global tracks on 20.01.2023
  Configurable<float> cfgTrkMaxChi2PerClusterITS{"cfgTrkMaxChi2PerClusterITS", 36.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of global tracks on 20.01.2023

  Configurable<LabeledArray<float>> cfgPIDCuts{"cfgPIDCuts", {hnm::pidcutsTable[0], 1, hnm::kNPIDLimits, hnm::speciesName, hnm::pidCutsName}, "Femtopartner PID nsigma selections"}; // PID selections

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

  Configurable<float> cfgMaxMultiplicity{"cfgMaxMultiplicity", 5000, "Maximum number of tracks in a collision (can be used to increase the S/B -> Very experimental)"};
  Configurable<float> cfgMinGGPtOverHNMPt{"cfgMinGGPtOverHNMPt", 0., "Minimum ratio of the pT of the gamma gamma pair over the pT of the HNM (can be used to increase the S/B)"};

  HistogramRegistry mHistManager{"HeavyNeutralMesonHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Prepare vectors for different species
  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;
  std::vector<ROOT::Math::PtEtaPhiMVector> etaPrimeEMC, etaPrimePCM, omegaEMC, omegaPCM, proton, antiproton, deuteron, antideuteron, pion, antipion;
  float mMassProton = constants::physics::MassProton;
  float mMassDeuteron = constants::physics::MassDeuteron;
  float mMassOmega = 0.782;
  float mMassEtaPrime = 0.957;
  float mMassPionCharged = constants::physics::MassPionCharged;

  Preslice<aod::V0PhotonsKF> perCollisionPCM = aod::v0photonkf::collisionId;
  Preslice<aod::SkimEMCClusters> perCollisionEMC = aod::skimmedcluster::collisionId;

  template <typename T>
  bool isSelectedTrack(T const& track, hnm::TracksPID partSpecies)
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
  bool isSelectedTrackPID(T const& track, hnm::TracksPID partSpecies)
  {
    bool isSelected = false;

    float nSigmaTrackTPC = track.tpcNSigmaPi();
    float nSigmaTrackTOF = track.tofNSigmaPi();
    float nSigmaTrackITS = track.itsNSigmaPi();

    float nSigmaTrackTPCTOF = std::sqrt(std::pow(nSigmaTrackTPC, 2) + std::pow(nSigmaTrackTOF, 2));

    if (track.p() <= cfgPtCuts->get(partSpecies, "P TOF thres")) {
      if (nSigmaTrackTPC > cfgPIDCuts->get(partSpecies, hnm::kTPCMin) &&
          nSigmaTrackTPC < cfgPIDCuts->get(partSpecies, hnm::kTPCMax) &&
          nSigmaTrackITS > cfgPIDCuts->get(partSpecies, hnm::kITSmin) &&
          nSigmaTrackITS < cfgPIDCuts->get(partSpecies, hnm::kITSmax)) {
        isSelected = true;
      }
    } else {
      if (nSigmaTrackTPCTOF < cfgPIDCuts->get(partSpecies, hnm::kTPCTOF)) {
        isSelected = true;
      }
    }
    return isSelected;
  }

  template <typename T>
  bool isSelectedEvent(T const& col)
  {
    if (confEvtSelectZvtx && std::abs(col.posZ()) > confEvtZvtx) {
      return false;
    }
    if (confEvtRequireSel8 && !col.sel8()) {
      return false;
    }
    if (col.multNTracksPV() > cfgMaxMultiplicity) {
      return false;
    }
    return true;
  }

  void init(InitContext const&)
  {
    mHistManager.add("Event/nGGs", "Number of (selected) #gamma#gamma paris;#bf{#it{N}_{#gamma#gamma}};#bf{#it{N}_{#gamma#gamma}^{selected}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nHeavyNeutralMesons", "Number of (selected) HNM candidates;#bf{#it{N}_{HNM}};#bf{#it{N}_{HNM}^{selected}}", HistType::kTH2F, {{51, -0.5, 50.5}, {51, -0.5, 50.5}});
    mHistManager.add("Event/nClustersVsV0s", "Number of clusters and V0s in the collision;#bf{#it{N}_{clusters}};#bf{#it{N}_{V0s}}", HistType::kTH2F, {{26, -0.5, 25.5}, {26, -0.5, 25.5}});
    mHistManager.add("Event/nEMCalEvents", "Number of collisions with a certain combination of EMCal triggers;;#bf{#it{N}_{collisions}}", HistType::kTH1F, {{5, -0.5, 4.5}});
    std::vector<std::string> nEventTitles = {"Cells & kTVXinEMC", "Cells & L0", "Cells & !kTVXinEMC & !L0", "!Cells & kTVXinEMC", "!Cells & L0"};
    for (size_t iBin = 0; iBin < nEventTitles.size(); iBin++)
      mHistManager.get<TH1>(HIST("Event/nEMCalEvents"))->GetXaxis()->SetBinLabel(iBin + 1, nEventTitles[iBin].data());
    mHistManager.add("Event/fMultiplicityBefore", "Multiplicity of all processed events;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("Event/fMultiplicityAfter", "Multiplicity after event cuts;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("Event/fZvtxBefore", "Zvtx of all processed events;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});
    mHistManager.add("Event/fZvtxAfter", "Zvtx after event cuts;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, -15, 15}});

    mHistManager.add("GG/invMassVsPt_PCM", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_PCMEMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_EMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{pT}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});

    // Momentum correlations p vs p_TPC
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationPos", "fMomCorrelation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
    mHistManager.add("TrackCuts/TracksBefore/fMomCorrelationNeg", "fMomCorrelation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});

    // All tracks
    mHistManager.add("TrackCuts/TracksBefore/fPtTrackBefore", "Transverse momentum of all processed tracks;#bf{#it{p}_{T} (GeV/#it{c})};#bf{#it{N}_{tracks}}", HistType::kTH1F, {{500, 0, 10}});
    mHistManager.add("TrackCuts/TracksBefore/fEtaTrackBefore", "Pseudorapidity of all processed tracks;#eta;#bf{#it{N}_{tracks}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("TrackCuts/TracksBefore/fPhiTrackBefore", "Azimuthal angle of all processed tracks;#phi;#bf{#it{N}_{tracks}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

    // TPC signal
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalTPCP", "TPCSignal;#bf{#it{p}_{TPC} (GeV/#it{c})};dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignal", "TPCSignalP;#bf{#it{p} (GeV/#it{c})};dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    // TPC signal antiparticles (negative charge)
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAntiTPCP", "TPCSignal;#bf{#it{p}_{TPC} (GeV/#it{c})};dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});
    mHistManager.add("TrackCuts/TPCSignal/fTPCSignalAnti", "TPCSignalP;#bf{#it{p} (GeV/#it{c})};dE/dx", {HistType::kTH2F, {{500, 0.0f, 6.0f}, {2000, -100.f, 500.f}}});

    const int nTrackSpecies = 2; // x2 because of anti particles
    const char* particleSpecies[nTrackSpecies] = {"Pion", "AntiPion"};
    const char* particleSpeciesLatex[nTrackSpecies] = {"#pi^{+}", "#pi^{-}"};

    for (int iParticle = 0; iParticle < nTrackSpecies; iParticle++) {
      mHistManager.add(Form("TrackCuts/TracksBefore/fMomCorrelationAfterCuts%s", particleSpecies[iParticle]), "fMomCorrelation;#bf{#it{p} (GeV/#it{c})};#bf{#it{p}_{TPC} (GeV/#it{c})}", {HistType::kTH2F, {{500, 0.0f, 20.0f}, {500, 0.0f, 20.0f}}});
      mHistManager.add(Form("TrackCuts/TPCSignal/fTPCSignal%s", particleSpecies[iParticle]), Form("%s TPC energy loss;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};dE/dx", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0.0f, 6.0f}, {10000, -100.f, 500.f}}});

      mHistManager.add(Form("TrackCuts/%s/fP", particleSpecies[iParticle]), Form("%s momentum at PV;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0, 10}});
      mHistManager.add(Form("TrackCuts/%s/fPt", particleSpecies[iParticle]), Form("%s transverse momentum;#bf{#it{p}_{T}^{%s} (GeV/#it{c})};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0, 10}});
      mHistManager.add(Form("TrackCuts/%s/fMomCorDif", particleSpecies[iParticle]), Form("Momentum correlation;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{p}_{TPC}^{%s} - #it{p}^{%s} (GeV/#it{c})}", particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0, 10}, {600, -3, 3}}});
      mHistManager.add(Form("TrackCuts/%s/fMomCorRatio", particleSpecies[iParticle]), Form("Relative momentum correlation;#bf{#it{p}^{%s} (GeV/#it{c})};#bf{#it{p}_{TPC}^{%s} - #it{p}^{%s} / #it{p}^{%s}}", particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{500, 0, 10}, {200, -1, 1}}});
      mHistManager.add(Form("TrackCuts/%s/fEta", particleSpecies[iParticle]), Form("%s pseudorapidity distribution;#eta;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -2, 2}});
      mHistManager.add(Form("TrackCuts/%s/fPhi", particleSpecies[iParticle]), Form("%s azimuthal angle distribution;#phi;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCvsTPCP", particleSpecies[iParticle]), Form("NSigmaTPC %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};n#sigma_{TPC}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTOFvsTPCP", particleSpecies[iParticle]), Form("NSigmaTOF %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};n#sigma_{TOF}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCTOFvsTPCP", particleSpecies[iParticle]), Form("NSigmaTPCTOF %s;#bf{#it{p}_{TPC}^{%s} (GeV/#it{c})};n#sigma_{comb}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaITSvsP", particleSpecies[iParticle]), Form("NSigmaITS %s;#bf{#it{p}^{%s} (GeV/#it{c})};n#sigma_{ITS}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCvsP", particleSpecies[iParticle]), Form("NSigmaTPC %s P;#bf{#it{p}^{%s} (GeV/#it{c})};n#sigma_{TPC}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTOFvsP", particleSpecies[iParticle]), Form("NSigmaTOF %s P;#bf{#it{p}^{%s} (GeV/#it{c})};n#sigma_{TOF}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});
      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCTOFvsP", particleSpecies[iParticle]), Form("NSigmaTPCTOF %s P;#bf{#it{p}^{%s} (GeV/#it{c})};n#sigma_{comb}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, 0.f, 10.f}}});

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
    mHistManager.fill(HIST("Event/fMultiplicityBefore"), collision.multNTracksPV());
    mHistManager.fill(HIST("Event/fZvtxBefore"), collision.posZ());

    // Ensure evts are consistent with Sel8 and Vtx-z selection
    if (!isSelectedEvent(collision))
      return;

    // QA accepted evts
    mHistManager.fill(HIST("Event/fMultiplicityAfter"), collision.multNTracksPV());
    mHistManager.fill(HIST("Event/fZvtxAfter"), collision.posZ());

    // clean vecs
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
      bool isPion = (isSelectedTrackPID(track, hnm::kPion) && isSelectedTrack(track, hnm::kPion));

      if (track.sign() > 0) { // Positive charge -> Particles
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
    // -------------------------------------------------------------------------------------
    processHNMs(vHNMs);
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

      if (heavyNeutralMeson.gg->isPi0 && massHNM > cfgMassWindowOmega->get("omega_min") && massHNM < cfgMassWindowOmega->get("omega_max") && heavyNeutralMeson.gg->pT() / heavyNeutralMeson.pT() > cfgMinGGPtOverHNMPt) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          omegaPCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), mMassOmega);
          mHistManager.fill(HIST("HNM/After/Omega/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/After/Omega/PCM/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/After/Omega/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          omegaEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), mMassOmega);
          mHistManager.fill(HIST("HNM/After/Omega/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/After/Omega/EMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/After/Omega/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        }
      } else if (heavyNeutralMeson.gg->isEta && massHNM > cfgMassWindowEtaPrime->get("etaprime_min") && massHNM < cfgMassWindowEtaPrime->get("etaprime_max") && heavyNeutralMeson.gg->pT() / heavyNeutralMeson.pT() > cfgMinGGPtOverHNMPt) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM) {
          etaPrimePCM.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), mMassEtaPrime);
          mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/After/EtaPrime/PCM/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
          etaPrimeEMC.emplace_back(heavyNeutralMeson.pT(), heavyNeutralMeson.eta(), RecoDecay::constrainAngle(heavyNeutralMeson.phi()), mMassEtaPrime);
          mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
          mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fEta"), heavyNeutralMeson.eta());
          mHistManager.fill(HIST("HNM/After/EtaPrime/EMC/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
        }
      } else {
        vHNMs.erase(vHNMs.begin() + iHNM);
        iHNM--;
      }
    }
    mHistManager.fill(HIST("Event/nHeavyNeutralMesons"), nHNMsBeforeMassCuts, vHNMs.size());
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMeson>(cfgc)}; }
