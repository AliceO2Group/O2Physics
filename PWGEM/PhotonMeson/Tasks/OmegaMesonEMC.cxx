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
/// \file OmegaMesonEMC.cxx
///
/// \brief This code loops over collisions to reconstruct heavy mesons (omega or eta') using EMCal clusters
///
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include "PWGEM/PhotonMeson/Utils/HNMUtilities.h"
#include "PWGJE/DataModel/EMCALMatchedCollisions.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
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
using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::EMCALMatchedCollisions>>;
using MyCollision = MyCollisions::iterator;
using SelectedTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TrackSelection, aod::TracksDCA, aod::pidTPCFullPi>>;
} // namespace o2::aod

namespace hnm
{

enum TrackCuts { kpT,
                 kEta,
                 kTPCSigma,
                 kTrackCuts
};

const std::vector<std::string> chargedPionMinMaxName{"Min", "Max"};
const std::vector<std::string> chargedPionCutsName{"pT", "eta", "TPC sigma"};
const float chargedPionCutsTable[kTrackCuts][2]{{0.35f, 20.f}, {-.8f, .8f}, {-4.f, 4.f}};

} // namespace hnm

struct OmegaMesonEMC {

  // --------------------------------> Configurables <------------------------------------
  // - Event selection cuts
  // - Track selection cuts
  // - Cluster shifts
  // - HNM mass selection windows
  // -------------------------------------------------------------------------------------
  // ---> Event selection
  Configurable<bool> confEvtSelectZvtx{"confEvtSelectZvtx", true, "Event selection includes max. z-Vertex"};
  Configurable<float> confEvtZvtx{"confEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> confEvtRequireSel8{"confEvtRequireSel8", true, "Evt sel: check for sel8 trigger bit"};
  Configurable<bool> confEvtRequirekTVXinEMC{"confEvtRequirekTVXinEMC", false, "Evt sel: check for EMCal MB trigger kTVXinEMC"};
  Configurable<bool> confEvtRequireL0{"confEvtRequireL0", false, "Evt sel: check for EMCal L0 trigger"};
  Configurable<bool> confEvtRequireGoodZVertex{"confEvtRequireGoodZVertex", false, "Evt sel: check for EMCal good z-vertex"};
  Configurable<bool> confEvtRequireNoSameBunchPileUp{"confEvtRequireNoSameBunchPileUp", false, "Evt sel: check for no same bunch pile-up"};

  // ---> Track selection
  Configurable<LabeledArray<float>> cfgChargedPionCuts{"cfgChargedPionCuts", {hnm::chargedPionCutsTable[0], 3, 2, hnm::chargedPionCutsName, hnm::chargedPionMinMaxName}, "Charged pion track cuts"};
  Configurable<float> cfgTPCNClustersMin{"cfgTPCNClustersMin", 80, "Mininum of TPC Clusters"};
  Configurable<float> cfgTrkTPCfCls{"cfgTrkTPCfCls", 0.83, "Minimum fraction of crossed rows over findable clusters"};
  Configurable<float> cfgTrkTPCcRowsMin{"cfgTrkTPCcRowsMin", 70, "Minimum number of crossed TPC rows"};
  Configurable<float> cfgTrkTPCsClsSharedFrac{"cfgTrkTPCsClsSharedFrac", 1.f, "Fraction of shared TPC clusters"};
  Configurable<float> cfgTrkITSnclsMin{"cfgTrkITSnclsMin", 4, "Minimum number of ITS clusters"};
  Configurable<float> cfgTrkDCAxyMax{"cfgTrkDCAxyMax", 0.15, "Maximum DCA_xy"};
  Configurable<float> cfgTrkDCAzMax{"cfgTrkDCAzMax", 0.3, "Maximum DCA_z"};
  Configurable<float> cfgTrkMaxChi2PerClusterTPC{"cfgTrkMaxChi2PerClusterTPC", 4.0f, "Minimal track selection: max allowed chi2 per TPC cluster"};  // 4.0 is default of global tracks on 20.01.2023
  Configurable<float> cfgTrkMaxChi2PerClusterITS{"cfgTrkMaxChi2PerClusterITS", 36.0f, "Minimal track selection: max allowed chi2 per ITS cluster"}; // 36.0 is default of global tracks on 20.01.2023

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
  static constexpr float DefaultMassWindows[2] = {0.11, 0.16};
  Configurable<LabeledArray<float>> cfgMassWindowPi0{"cfgMassWindowPi0", {DefaultMassWindows, 2, {"min", "max"}}, "Mass window for selected decay pi0"};

  Configurable<int32_t> cfgMaxMultiplicity{"cfgMaxMultiplicity", 5000, "Maximum number of tracks in a collision (can be used to increase the S/B -> Very experimental)"};
  Configurable<float> cfgMinGGPtOverHNMPt{"cfgMinGGPtOverHNMPt", 0., "Minimum ratio of the pT of the gamma gamma pair over the pT of the HNM (can be used to increase the S/B)"};

  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < confEvtZvtx;
  Filter collisionMultFilter = (o2::aod::mult::multNTracksPV <= cfgMaxMultiplicity);

  Filter trackPtFilter = (o2::aod::track::pt > cfgChargedPionCuts->get(hnm::kpT, "Min")) && (o2::aod::track::pt < cfgChargedPionCuts->get(hnm::kpT, "Max"));
  Filter trackEtaFilter = (o2::aod::track::eta > cfgChargedPionCuts->get(hnm::kEta, "Min")) && (o2::aod::track::eta < cfgChargedPionCuts->get(hnm::kEta, "Max"));
  Filter trackDCAXYFilter = nabs(o2::aod::track::dcaXY) < cfgTrkDCAxyMax;
  Filter trackDCAZFilter = nabs(o2::aod::track::dcaZ) < cfgTrkDCAzMax;

  Filter trackTPCChi2Filter = o2::aod::track::tpcChi2NCl < cfgTrkMaxChi2PerClusterTPC;
  Filter trackITSChi2Filter = o2::aod::track::itsChi2NCl < cfgTrkMaxChi2PerClusterITS;

  Filter trackTPCSigmaFilterTPC = (o2::aod::pidtpc::tpcNSigmaPi > cfgChargedPionCuts->get(hnm::kTPCSigma, "Min")) && (o2::aod::pidtpc::tpcNSigmaPi < cfgChargedPionCuts->get(hnm::kTPCSigma, "Max"));

  template <typename T>
  bool isSelectedTrack(T const& track)
  {
    if (track.tpcNClsFound() < cfgTPCNClustersMin)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgTrkTPCfCls)
      return false;
    if (track.tpcNClsCrossedRows() < cfgTrkTPCcRowsMin)
      return false;
    if (track.tpcFractionSharedCls() > cfgTrkTPCsClsSharedFrac)
      return false;
    if (track.itsNCls() < cfgTrkITSnclsMin)
      return false;
    return true;
  }

  HistogramRegistry mHistManager{"HeavyNeutralMesonHistograms", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Prepare vectors for different species
  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;
  std::vector<ROOT::Math::PtEtaPhiMVector> pion, antipion;

  Preslice<aod::SkimEMCClusters> perCollisionEMC = aod::skimmedcluster::collisionId;

  void init(InitContext const&)
  {
    mHistManager.add("Event/nEMCalEvents", "Number of collisions with a certain combination of EMCal triggers;;#bf{#it{N}_{collisions}}", HistType::kTH1F, {{5, -0.5, 4.5}});
    std::vector<std::string> nEventTitles = {"Cells & kTVXinEMC", "Cells & L0", "Cells & !kTVXinEMC & !L0", "!Cells & kTVXinEMC", "!Cells & L0"};
    for (size_t iBin = 0; iBin < nEventTitles.size(); iBin++)
      mHistManager.get<TH1>(HIST("Event/nEMCalEvents"))->GetXaxis()->SetBinLabel(iBin + 1, nEventTitles[iBin].data());
    mHistManager.add("Event/fMultiplicity", "Multiplicity after event cuts;#bf{#it{N}_{tracks}};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{500, 0, 500}});
    mHistManager.add("Event/fZvtx", "Zvtx after event cuts;#bf{z_{vtx} (cm)};#bf{#it{N}_{collisions}}", HistType::kTH1F, {{300, -15, 15}});

    mHistManager.add("GG/invMassVsPt", "Invariant mass and pT of gg candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_selected", "Invariant mass and pT of gg candidates;#bf{#it{M}^{#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});

    const int nTrackSpecies = 2; // x2 because of anti particles
    const char* particleSpecies[nTrackSpecies] = {"Pion", "AntiPion"};
    const char* particleSpeciesLatex[nTrackSpecies] = {"#pi^{+}", "#pi^{-}"};

    for (int iParticle = 0; iParticle < nTrackSpecies; iParticle++) {
      mHistManager.add(Form("TrackCuts/%s/fPt", particleSpecies[iParticle]), Form("%s transverse momentum;#bf{#it{p}_{T}^{%s} (GeV/#it{c})};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0, 10}});
      mHistManager.add(Form("TrackCuts/%s/fEta", particleSpecies[iParticle]), Form("%s pseudorapidity distribution;#eta;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -2, 2}});
      mHistManager.add(Form("TrackCuts/%s/fPhi", particleSpecies[iParticle]), Form("%s azimuthal angle distribution;#phi;#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

      mHistManager.add(Form("TrackCuts/%s/fNsigmaTPCvsP", particleSpecies[iParticle]), Form("NSigmaTPC %s P;#bf{#it{p}^{%s} (GeV/#it{c})};n#sigma_{TPC}^{%s}", particleSpecies[iParticle], particleSpeciesLatex[iParticle], particleSpeciesLatex[iParticle]), {HistType::kTH2F, {{100, 0.0f, 10.0f}, {100, -10.f, 10.f}}});

      mHistManager.add(Form("TrackCuts/%s/fDCAxy", particleSpecies[iParticle]), Form("fDCAxy %s;#bf{DCA_{xy}};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -0.5f, 0.5f}});
      mHistManager.add(Form("TrackCuts/%s/fDCAz", particleSpecies[iParticle]), Form("fDCAz %s;#bf{DCA_{z}};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, -0.5f, 0.5f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCsCls", particleSpecies[iParticle]), Form("fTPCsCls %s;#bf{TPC Shared Clusters};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCcRows", particleSpecies[iParticle]), Form("fTPCcRows %s;#bf{TPC Crossed Rows};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTrkTPCfCls", particleSpecies[iParticle]), Form("fTrkTPCfCls %s;#bf{TPC Findable/CrossedRows};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{500, 0.0f, 3.0f}});
      mHistManager.add(Form("TrackCuts/%s/fTPCncls", particleSpecies[iParticle]), Form("fTPCncls %s;#bf{TPC Clusters};#bf{#it{N}^{%s}}", particleSpecies[iParticle], particleSpeciesLatex[iParticle]), HistType::kTH1F, {{163, -1.0f, 162.0f}});
    }

    // --> HNM QA
    // Properties of the pi+pi- pair
    mHistManager.add("Omega/Before/PiPlPiMi/fInvMassVsPt", "Invariant mass and pT of #pi^+pi^- pairs;#bf{#it{M}^{#pi^{+}#pi^{-}} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}} (GeV/#it{c})}", HistType::kTH2F, {{400, 0.2, 1.}, {250, 0., 25.}});
    mHistManager.add("Omega/Before/PiPlPiMi/fEta", "Pseudorapidity of HMNCand;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}}}", HistType::kTH1F, {{500, -2, 2}});
    mHistManager.add("Omega/Before/PiPlPiMi/fPhi", "Azimuthal angle of HMNCand;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});

    for (const auto& BeforeAfterString : {"Before", "After"}) {
      mHistManager.add(Form("Omega/%s/fInvMassVsPt", BeforeAfterString), "Invariant mass and pT of heavy neutral meson candidates;#bf{#it{M}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c}^{2})};#bf{#it{p}_{T}^{#pi^{+}#pi^{-}#gamma#gamma} (GeV/#it{c})}", HistType::kTH2F, {{800, 0.4, 1.2}, {250, 0., 25.}});
      mHistManager.add(Form("Omega/%s/fEta", BeforeAfterString), "Pseudorapidity of HNM candidate;#eta;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{500, -2, 2}});
      mHistManager.add(Form("Omega/%s/fPhi", BeforeAfterString), "Azimuthal angle of HNM candidate;#phi;#bf{#it{N}^{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH1F, {{720, 0, constants::math::TwoPI}});
    }
    if (cfgDoEMCShift.value) {
      for (int iSM = 0; iSM < nSMs; iSM++) {
        emcEtaShift[iSM] = cfgEMCEtaShift.value[iSM];
        emcPhiShift[iSM] = cfgEMCPhiShift.value[iSM];
        LOG(info) << "SM-wise shift in eta/phi for SM " << iSM << ": " << emcEtaShift[iSM] << " / " << emcPhiShift[iSM];
      }
    }
  }

  void process(aod::MyCollision const& collision, aod::MyBCs const&, aod::SkimEMCClusters const& clusters, aod::SelectedTracks const& tracks)
  {
    // clean vecs
    pion.clear();
    antipion.clear();
    vHNMs.clear();

    // ---------------------------------> EMCal event QA <----------------------------------
    // - Fill Event/nEMCalEvents histogram for EMCal event QA
    // -------------------------------------------------------------------------------------
    bool bcHasEMCCells = collision.isemcreadout();
    auto foundBC = collision.foundBC_as<aod::MyBCs>();
    bool iskTVXinEMC = foundBC.alias_bit(kTVXinEMC);
    bool isL0Triggered = foundBC.alias_bit(kEMC7) || foundBC.alias_bit(kDMC7) || foundBC.alias_bit(kEG1) || foundBC.alias_bit(kEG2);

    if (confEvtRequireSel8 && !collision.sel8())
      return; // Skip this collision if sel8 trigger bit is not set
    if (confEvtRequirekTVXinEMC && !iskTVXinEMC)
      return; // Skip this collision if kTVXinEMC trigger bit is not set
    if (confEvtRequireL0 && !isL0Triggered)
      return; // Skip this collision if L0 trigger bit is not set
    if (confEvtRequireGoodZVertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return; // Skip this collision if good z-vertex condition is not met
    if (confEvtRequireNoSameBunchPileUp && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return; // Skip this collision if no same bunch pileup condition is not met

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

    mHistManager.fill(HIST("Event/fMultiplicity"), collision.multNTracksPV());
    mHistManager.fill(HIST("Event/fZvtx"), collision.posZ());

    // --------------------------------> Process Photons <----------------------------------
    // - Slice clusters and V0s by collision ID to get the ones in this collision
    // - Store the clusters and V0s in the vGammas vector
    // - Reconstruct gamma-gamma pairs
    // -------------------------------------------------------------------------------------
    auto clustersInThisCollision = clusters.sliceBy(perCollisionEMC, collision.globalIndex());

    std::vector<hnmutilities::Photon> vGammas;
    hnmutilities::storeGammasInVector(clustersInThisCollision, vGammas, emcEtaShift, emcPhiShift);
    hnmutilities::reconstructGGs(vGammas, vGGs);
    vGammas.clear();

    for (unsigned int iGG = 0; iGG < vGGs.size(); iGG++) {
      auto lightMeson = &vGGs.at(iGG);

      mHistManager.fill(HIST("GG/invMassVsPt"), lightMeson->m(), lightMeson->pT());

      if (lightMeson->m() > cfgMassWindowPi0->get("min") && lightMeson->m() < cfgMassWindowPi0->get("max")) {
        lightMeson->isPi0 = true;
        mHistManager.fill(HIST("GG/invMassVsPt_selected"), lightMeson->m(), lightMeson->pT());
      } else {
        vGGs.erase(vGGs.begin() + iGG);
        iGG--;
      }
    }

    // ------------------------------> Loop over all tracks <-------------------------------
    // - Fill QA histograms for all tracks and per particle species
    // -------------------------------------------------------------------------------------
    for (const auto& track : tracks) {
      if (!isSelectedTrack(track))
        continue;             // Skip tracks that do not pass the selection criteria
      if (track.sign() > 0) { // Positive charge -> Particles
        pion.emplace_back(track.pt(), track.eta(), track.phi(), constants::physics::MassPionCharged);

        mHistManager.fill(HIST("TrackCuts/Pion/fPt"), track.pt());
        mHistManager.fill(HIST("TrackCuts/Pion/fEta"), track.eta());
        mHistManager.fill(HIST("TrackCuts/Pion/fPhi"), track.phi());

        mHistManager.fill(HIST("TrackCuts/Pion/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPi());

        mHistManager.fill(HIST("TrackCuts/Pion/fDCAxy"), track.dcaXY());
        mHistManager.fill(HIST("TrackCuts/Pion/fDCAz"), track.dcaZ());
        mHistManager.fill(HIST("TrackCuts/Pion/fTPCsCls"), track.tpcNClsShared());
        mHistManager.fill(HIST("TrackCuts/Pion/fTPCcRows"), track.tpcNClsCrossedRows());
        mHistManager.fill(HIST("TrackCuts/Pion/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
        mHistManager.fill(HIST("TrackCuts/Pion/fTPCncls"), track.tpcNClsFound());
      } else { // Negative charge -> Anti-particles
        antipion.emplace_back(track.pt(), track.eta(), track.phi(), constants::physics::MassPionCharged);

        mHistManager.fill(HIST("TrackCuts/AntiPion/fPt"), track.pt());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fEta"), track.eta());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fPhi"), track.phi());

        mHistManager.fill(HIST("TrackCuts/AntiPion/fNsigmaTPCvsP"), track.p(), track.tpcNSigmaPi());

        mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAxy"), track.dcaXY());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fDCAz"), track.dcaZ());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCsCls"), track.tpcNClsShared());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCcRows"), track.tpcNClsCrossedRows());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fTrkTPCfCls"), track.tpcCrossedRowsOverFindableCls());
        mHistManager.fill(HIST("TrackCuts/AntiPion/fTPCncls"), track.tpcNClsFound());
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

        mHistManager.fill(HIST("Omega/Before/PiPlPiMi/fInvMassVsPt"), vecPiPlPiMi.M(), vecPiPlPiMi.pt());
        mHistManager.fill(HIST("Omega/Before/PiPlPiMi/fEta"), vecPiPlPiMi.eta());
        mHistManager.fill(HIST("Omega/Before/PiPlPiMi/fPhi"), RecoDecay::constrainAngle(vecPiPlPiMi.phi()));
      }
    }

    // ---------------------------> Process HNM candidates <--------------------------------
    // - Fill invMassVsPt histograms separated into HNM types (based on GG mass) and gamma reco method
    // -------------------------------------------------------------------------------------
    for (unsigned int iHNM = 0; iHNM < vHNMs.size(); iHNM++) {
      auto heavyNeutralMeson = vHNMs.at(iHNM);
      float massHNM = heavyNeutralMeson.m(cfgHNMMassCorrection);

      mHistManager.fill(HIST("Omega/Before/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
      mHistManager.fill(HIST("Omega/Before/fEta"), heavyNeutralMeson.eta());
      mHistManager.fill(HIST("Omega/Before/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));

      if (heavyNeutralMeson.gg->pT() / heavyNeutralMeson.pT() > cfgMinGGPtOverHNMPt) {
        mHistManager.fill(HIST("Omega/After/fInvMassVsPt"), massHNM, heavyNeutralMeson.pT());
        mHistManager.fill(HIST("Omega/After/fEta"), heavyNeutralMeson.eta());
        mHistManager.fill(HIST("Omega/After/fPhi"), RecoDecay::constrainAngle(heavyNeutralMeson.phi()));
      }
    }
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OmegaMesonEMC>(cfgc, TaskName{"omega-meson-emc"})}; }
