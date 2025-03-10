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
/// \author Nicolas Strangmann (nicolas.strangmann@cern.ch) - Goethe University Frankfurt
///

#include <vector>

#include "PWGEM/PhotonMeson/Utils/HNMUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::pwgem::photonmeson;

using MyBCs = soa::Join<aod::BCs, aod::BcSels, aod::Timestamps>;
using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using SelectedTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TrackSelection>>;

struct HeavyNeutralMesonFilter {
  Produces<aod::HeavyNeutralMesonFilters> tags;

  HistogramRegistry mHistManager{"HeavyNeutralMesonFilterHistograms", {}, OutputObjHandlingPolicy::QAObject};

  Configurable<float> cfgTrackMinPt{"cfgTrackMinPt", 0.1, "Minimum momentum of tracks (GeV/c)"};
  Configurable<int> cfgHNMMassCorrection{"cfgHNMMassCorrection", 1, "Use GG PDG mass to correct HNM mass (0 = off, 1 = subDeltaPi0, 2 = subLambda)"};
  static constexpr float defaultMassWindows[2][4] = {{0., 0.4, 0.6, 1.}, {0.4, 0.8, 0.8, 1.2}};
  Configurable<LabeledArray<float>> massWindowOmega{"massWindowOmega", {defaultMassWindows[0], 4, {"pi0_min", "pi0_max", "omega_min", "omega_max"}}, "Mass window for selected omegas and their decay pi0"};
  Configurable<LabeledArray<float>> massWindowEtaPrime{"massWindowEtaPrime", {defaultMassWindows[1], 4, {"eta_min", "eta_max", "etaprime_min", "etaprime_max"}}, "Mass window for selected eta' and their decay eta"};

  static constexpr float defaultMinPts[4] = {1.8, 1.8, 2.6, 2.6};
  Configurable<LabeledArray<float>> minHNMPts{"minHNMPts", {defaultMinPts, 4, {"PCM_omega", "PCM_etaprime", "EMC_omega", "EMC_etaprime"}}, "Minimum pT values for the trigger decisions (GeV/c)"};

  Filter trackPtFilter = aod::track::pt > cfgTrackMinPt;

  std::vector<hnmutilities::GammaGammaPair> vGGs;
  std::vector<hnmutilities::HeavyNeutralMeson> vHNMs;

  bool colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime = false;

  emcal::Geometry* emcalGeom;

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

    mHistManager.add("GG/invMassVsPt_PCM", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{N}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_PCMEMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{N}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});
    mHistManager.add("GG/invMassVsPt_EMC", "Invariant mass and pT of gg candidates;#bf{#it{M}_{#gamma#gamma}};#bf{#it{N}_{#gamma#gamma}}", HistType::kTH2F, {{400, 0., 0.8}, {250, 0., 25.}});

    mHistManager.add("HeavyNeutralMeson/invMassVsPt_PCM", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{N}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
    mHistManager.add("HeavyNeutralMeson/invMassVsPt_PCMEMC", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{N}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
    mHistManager.add("HeavyNeutralMeson/invMassVsPt_EMC", "Invariant mass and pT of HNM candidates;#bf{#it{M}_{#pi^{+}#pi^{-}#gamma#gamma}};#bf{#it{N}_{#pi^{+}#pi^{-}#gamma#gamma}}", HistType::kTH2F, {{600, 0.6, 1.2}, {250, 0., 25.}});
  }

  Preslice<aod::V0PhotonsKF> perCollision_pcm = aod::v0photonkf::collisionId;
  Preslice<aod::SkimEMCClusters> perCollision_emc = aod::skimmedcluster::collisionId;

  void process(MyCollisions::iterator const& collision, MyBCs const&, aod::SkimEMCClusters const& clusters, aod::V0PhotonsKF const& v0s, SelectedTracks const& tracks)
  {
    mHistManager.fill(HIST("Event/hCollisionCounter"), 0.);

    if (collision.foundBC_as<MyBCs>().alias_bit(kTVXinEMC))
      mHistManager.fill(HIST("Event/hCollisionCounter"), 1.);

    auto v0sInThisCollision = v0s.sliceBy(perCollision_pcm, collision.globalIndex());
    auto clustersInThisCollision = clusters.sliceBy(perCollision_emc, collision.globalIndex());

    mHistManager.fill(HIST("Event/nClustersVsV0s"), clustersInThisCollision.size(), v0sInThisCollision.size());
    mHistManager.fill(HIST("Event/nTracks"), tracks.size());

    colContainsPCMOmega = colContainsEMCOmega = colContainsPCMEtaPrime = colContainsEMCEtaPrime = false;

    hnmutilities::reconstructGGs(clustersInThisCollision, v0sInThisCollision, vGGs);
    processGGs(vGGs);
    hnmutilities::reconstructHeavyNeutralMesons(tracks, vGGs, vHNMs);
    processHNMs(vHNMs);

    tags(colContainsPCMOmega, colContainsEMCOmega, colContainsPCMEtaPrime, colContainsEMCEtaPrime);
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
      } else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC) {
        mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt_EMC"), massHNM, heavyNeutralMeson.pT());
      } else {
        mHistManager.fill(HIST("HeavyNeutralMeson/invMassVsPt_PCMEMC"), massHNM, heavyNeutralMeson.pT());
      }

      if (heavyNeutralMeson.gg->isPi0 && massHNM > massWindowOmega->get("omega_min") && massHNM < massWindowOmega->get("omega_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM && heavyNeutralMeson.pT() > minHNMPts->get("PCM_omega"))
          colContainsPCMOmega = true;
        else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC && heavyNeutralMeson.pT() > minHNMPts->get("EMC_omega"))
          colContainsEMCOmega = true;
      } else if (heavyNeutralMeson.gg->isEta && massHNM > massWindowEtaPrime->get("etaprime_min") && massHNM < massWindowEtaPrime->get("etaprime_max")) {
        if (heavyNeutralMeson.gg->reconstructionType == photonpair::kPCMPCM && heavyNeutralMeson.pT() > minHNMPts->get("PCM_etaprime"))
          colContainsPCMEtaPrime = true;
        else if (heavyNeutralMeson.gg->reconstructionType == photonpair::kEMCEMC && heavyNeutralMeson.pT() > minHNMPts->get("EMC_etaprime"))
          colContainsEMCEtaPrime = true;
      } else {
        vHNMs.erase(vHNMs.begin() + iHNM);
        iHNM--;
      }
    }
    mHistManager.fill(HIST("Event/nHeavyNeutralMesons"), nHNMsBeforeMassCuts, vHNMs.size());

    if (colContainsPCMOmega)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 2.);
    if (colContainsEMCOmega)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 3.);
    if (colContainsPCMEtaPrime)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 4.);
    if (colContainsEMCEtaPrime)
      mHistManager.fill(HIST("Event/hCollisionCounter"), 5.);
  }
};

WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HeavyNeutralMesonFilter>(cfgc)};
}
