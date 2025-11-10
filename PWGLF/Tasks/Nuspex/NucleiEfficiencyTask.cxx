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
// Authors: Rafael Manhart,
// Date: 06.05.2024

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/DataModel/spectraTOF.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"
#include <TF1.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct NucleiEfficiencyTask {

  HistogramRegistry MC_gen_reg{"MC_particles_gen", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg{"MC_particles_reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> histPDG_gen{TH1F("PDG_gen", "PDG;PDG code", 18, 0.0, 18)};
  OutputObj<TH1F> histPDG_gen_reco{TH1F("PDG_gen_reco", "PDG;PDG code", 18, 0.0, 18)};
  OutputObj<TH1F> histPDG_reco{TH1F("PDG_reco", "PDG;PDG code", 18, 0.0, 18)};

  void init(o2::framework::InitContext&)
  {

    std::vector<double> ptBinning = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4., 5., 6., 8., 10., 12., 14.};
    std::vector<double> etaBinning = {-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    std::vector<double> PDGBinning = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0};

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pAxis = {ptBinning, "#it{p} (GeV/#it{c})"};
    AxisSpec centralityAxis = {100, 0.0, 105.0, "VT0C (%)"};
    AxisSpec etaAxis = {etaBinning, "#eta"};
    AxisSpec ImPaAxis = {100, 0.0, 105.0, "Impact parameter"};
    AxisSpec PDGBINNING = {PDGBinning, "PDG code"};

    // *********************** Generated **********************
    MC_gen_reg.add("histGenVtxMC", "MC generated vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_gen_reg.add("histCentrality", "Impact parameter", HistType::kTH1F, {centralityAxis});
    MC_gen_reg.add("hist_gen_p", "generated p distribution", HistType::kTH2F, {pAxis, PDGBINNING});
    MC_gen_reg.add("hist_gen_pT", "generated p_{T} distribution", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_gen_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_gen_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_gen_reg.add("histRapid", "#gamma", HistType::kTH2F, {{1000, -5.0, 5.0}, PDGBINNING});

    // *********************** Generated reco **********************

    MC_gen_reg.add("histGenVtxMC_reco", "MC generated (reco) vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_gen_reg.add("histCentrality_reco", "Centrality", HistType::kTH1F, {centralityAxis});
    MC_gen_reg.add("histEta_reco", "generated (reco) #eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_gen_reg.add("hist_gen_reco_p", "generated (reco) p distribution", HistType::kTH2F, {pAxis, PDGBINNING});
    MC_gen_reg.add("hist_gen_reco_pT", "generated (reco) p_{T} distribution", HistType::kTH2F, {ptAxis, PDGBINNING});

    // ********************** Reconstructed *********************
    MC_recon_reg.add("histRecVtxMC", "MC reconstructed vertex z position", HistType::kTH1F, {{400, -40., +40., "z position (cm)"}});
    MC_recon_reg.add("histCentrality", "Centrality", HistType::kTH1F, {centralityAxis});
    MC_recon_reg.add("histPhi", "#phi", HistType::kTH2F, {{100, 0., 2. * TMath::Pi()}, PDGBINNING});
    MC_recon_reg.add("histEta", "#eta", HistType::kTH2F, {{102, -2.01, 2.01}, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_vs_p", "ITS reconstructed p distribution", HistType::kTH2F, {pAxis, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_TPC_vs_p", "ITS_TPC reconstructed p distribution", HistType::kTH2F, {pAxis, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_TPC_TOF_vs_p", "ITS_TPC_TOF reconstructed p distribution", HistType::kTH2F, {pAxis, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_vs_pT", "ITS reconstructed p_{T} distribution", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_TPC_vs_pT", "ITS_TPC reconstructed p_{T} distribution", HistType::kTH2F, {ptAxis, PDGBINNING});
    MC_recon_reg.add("hist_rec_ITS_TPC_TOF_vs_pT", "ITS_TPC_TOF reconstructed p_{T} distribution", HistType::kTH2F, {ptAxis, PDGBINNING});
  }

  // ************************ Configurables ***********************
  Configurable<bool> event_selection_MC_sel8{"event_selection_MC_sel8", true, "Enable sel8 event selection in MC processing"};
  Configurable<bool> applyPvZCutGenColl{"applyPvZCutGenColl", true, "applyPvZCutGenColl"};
  Configurable<float> yMin{"yMin", -0.5, "Minimum rapidity"};
  Configurable<float> yMax{"yMax", 0.5, "Maximum rapidity"};
  Configurable<float> p_min{"p_min", 0.1f, "min track.pt()"};
  Configurable<float> p_max{"p_max", 1e+10f, "max track.pt()"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> minCentrality{"minCentrality", 0.0, "min Centrality used"};
  Configurable<float> maxCentrality{"maxCentrality", 80.0, "max Centrality used"};
  Configurable<bool> enable_Centrality_cut{"enable_Centrality_cut", true, "enable Centrality cut"};

  // Track filter
  Configurable<bool> requireITS{"requireITS", true, "Additional cut on the ITS requirement"};
  Configurable<bool> requireTPC{"requireTPC", true, "Additional cut on the TPC requirement"};
  Configurable<bool> passedITSRefit{"passedITSRefit", true, "Additional cut on the ITS refit requirement"};
  Configurable<bool> passedTPCRefit{"passedTPCRefit", true, "Additional cut on the TPC refit requirement"};
  Configurable<float> minReqClusterITS{"minReqClusterITS", 1.0, "min number of clusters required in ITS"};
  Configurable<float> minReqClusterITSib{"minReqClusterITSib", 1.0, "min number of clusters required in ITS inner barrel"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 0.0f, "min number of crossed rows TPC"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 70.0f, "min number of crossed rows TPC"};
  Configurable<float> minRatioCrossedRowsTPC{"minRatioCrossedRowsTPC", 0.8f, "min ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxRatioCrossedRowsTPC{"maxRatioCrossedRowsTPC", 2.0f, "max ratio of crossed rows over findable clusters TPC"};
  Configurable<float> maxChi2PerClusterTPC{"maxChi2PerClusterTPC", 4.f, "Cut on the maximum value of the chi2 per cluster in the TPC"};
  Configurable<float> minChi2PerClusterTPC{"minChi2PerClusterTPC", 0.5f, "Cut on the minimum value of the chi2 per cluster in the TPC"};
  Configurable<float> maxChi2PerClusterITS{"maxChi2PerClusterITS", 36.f, "Cut on the maximum value of the chi2 per cluster in the ITS"};
  Configurable<float> maxDcaXYFactor{"maxDcaXYFactor", 0.5f, "DCA xy factor"};
  Configurable<float> maxDCA_Z{"maxDCA_Z", 2.0f, "max DCA to vertex z"};
  Configurable<int> lastRequiredTrdCluster{"lastRequiredTrdCluster", -1, "Last cluster to required in TRD for track selection. -1 does not require any TRD cluster"};
  Configurable<bool> requireGoldenChi2{"requireGoldenChi2", false, "Enable the requirement of GoldenChi2"};

  Configurable<bool> calc_cent{"calc_cent", false, "Enable centrality processing"};

  Configurable<bool> eta_cut_MC_gen{"eta_cut_MC_gen", true, "Enable eta cut for generated MC"};
  Configurable<bool> use_pT_cut{"use_pT_cut", true, "0: p is used | 1: pT is used"};

  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove TF border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove TF border"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Remove TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Remove TF border"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove TF border"};

  //***********************************************************************************

  template <typename particleType>
  bool isInAcceptance(const particleType& particle)
  {
    if (particle.pt() < p_min || particle.pt() > p_max)
      return false;
    if (particle.eta() < -cfgCutEta || particle.eta() > cfgCutEta)
      return false;
    // if (particle.phi() < phiMin || particle.phi() > phiMax) return false;
    if (particle.y() < yMin || particle.y() > yMax)
      return false;
    if (!particle.isPhysicalPrimary())
      return false;

    return true;
  }

  //***********************************************************************************

  template <typename CollisionType>
  bool isEventSelected(CollisionType const& collision)
  {
    if (removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;

    return true;
  }

  //***********************************************************************************

  template <typename CollType>
  bool isCollisionSelected(const CollType& collision)
  {
    if (event_selection_MC_sel8 && !collision.sel8())
      return false;
    if (collision.posZ() < -cfgCutVertex || collision.posZ() > cfgCutVertex)
      return false;
    if (removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
      return false;
    if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
      return false;
    if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
      return false;
    if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
      return false;
    if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
      return false;

    return true;
  }

  //***********************************************************************************

  template <typename trackType>
  bool isTrackSelected(trackType& track)
  {
    if (!track.has_mcParticle())
      return false;

    const auto mcParticle = track.mcParticle();
    if (!isInAcceptance(mcParticle))
      return false; // pt eta phi y
    if (!track.has_collision())
      return false;

    float TPCnumberClsFound = track.tpcNClsFound();
    float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
    float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
    float Chi2perClusterTPC = track.tpcChi2NCl();
    float Chi2perClusterITS = track.itsChi2NCl();

    bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f))));

    if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z || TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS)
      return false;
    if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
      return false;
    if (requireGoldenChi2 && !(track.passedGoldenChi2()))
      return false;

    return true;
  }

  //***********************************************************************************
  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, aod::CentFT0Cs>;
  using CollisionCandidatesMC = o2::soa::Join<CollisionCandidates, o2::aod::McCollisionLabels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TrackSelection, o2::aod::TrackSelectionExtension, o2::aod::TracksDCA>;
  using TrackCandidatesMC = o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels>;

  SliceCache cache;
  Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  Preslice<o2::aod::McParticles> perCollisionMc = o2::aod::mcparticle::mcCollisionId;
  PresliceUnsorted<CollisionCandidatesMC> collPerCollMc = o2::aod::mccollisionlabel::mcCollisionId;

  void processMC(o2::aod::McCollisions const& mcCollisions,
                 // o2::soa::SmallGroups<CollisionCandidatesMC> const& collisions,
                 CollisionCandidatesMC const& collisions,
                 TrackCandidatesMC const& tracks,
                 o2::aod::McParticles const& mcParticles)
  {
    /// loop over generated collisions
    for (const auto& mcCollision : mcCollisions) {

      const auto groupedCollisions = collisions.sliceBy(collPerCollMc, mcCollision.globalIndex());
      const auto groupedMcParticles = mcParticles.sliceBy(perCollisionMc, mcCollision.globalIndex());

      if (groupedCollisions.size() < 1)
        continue;
      float centrality = -1.;

      /// loop over reconstructed collisions
      for (const auto& collision : groupedCollisions) {
        if (!isCollisionSelected(collision))
          continue;

        centrality = collision.centFT0C();
        if (centrality < minCentrality || centrality > maxCentrality)
          continue;

        MC_recon_reg.fill(HIST("histCentrality"), centrality);
        MC_recon_reg.fill(HIST("histRecVtxMC"), collision.posZ());

        const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

        // Track loop
        for (const auto& track : groupedTracks) {
          if (!isTrackSelected(track))
            continue;

          const auto& particle = track.mcParticle();
          // TLorentzVector lorentzVector_particle_MC{};

          int pdgbin = -10;
          switch (particle.pdgCode()) {
            case +211:
              histPDG_reco->AddBinContent(1);
              pdgbin = 0;
              break;
            case -211:
              histPDG_reco->AddBinContent(2);
              pdgbin = 1;
              break;
            case +321:
              histPDG_reco->AddBinContent(3);
              pdgbin = 2;
              break;
            case -321:
              histPDG_reco->AddBinContent(4);
              pdgbin = 3;
              break;
            case +2212:
              histPDG_reco->AddBinContent(5);
              pdgbin = 4;
              break;
            case -2212:
              histPDG_reco->AddBinContent(6);
              pdgbin = 5;
              break;
            case +1000010020:
              histPDG_reco->AddBinContent(7);
              pdgbin = 6;
              break;
            case -1000010020:
              histPDG_reco->AddBinContent(8);
              pdgbin = 7;
              break;
            case +1000010030:
              histPDG_reco->AddBinContent(9);
              pdgbin = 8;
              break;
            case -1000010030:
              histPDG_reco->AddBinContent(10);
              pdgbin = 9;
              break;
            case +1000020030:
              histPDG_reco->AddBinContent(11);
              pdgbin = 10;
              break;
            case -1000020030:
              histPDG_reco->AddBinContent(12);
              pdgbin = 11;
              break;
            case +1000020040:
              histPDG_reco->AddBinContent(13);
              pdgbin = 12;
              break;
            case -1000020040:
              histPDG_reco->AddBinContent(14);
              pdgbin = 13;
              break;
            default:
              pdgbin = -10;
              continue;
              break;
          }

          MC_recon_reg.fill(HIST("histPhi"), track.phi(), pdgbin);
          MC_recon_reg.fill(HIST("histEta"), track.eta(), pdgbin);

          if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
            if (track.hasITS()) {
              MC_recon_reg.fill(HIST("hist_rec_ITS_vs_p"), track.p() * 2, pdgbin);
              MC_recon_reg.fill(HIST("hist_rec_ITS_vs_pT"), track.pt() * 2, pdgbin);
              if (track.hasTPC()) {
                MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_vs_p"), track.p() * 2, pdgbin);
                MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_vs_pT"), track.pt() * 2, pdgbin);
                if (track.hasTOF()) {
                  MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_TOF_vs_p"), track.p() * 2, pdgbin);
                  MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_TOF_vs_pT"), track.pt() * 2, pdgbin);
                }
              }
            }
          } else {
            if (track.hasITS()) {
              MC_recon_reg.fill(HIST("hist_rec_ITS_vs_p"), track.p(), pdgbin);
              MC_recon_reg.fill(HIST("hist_rec_ITS_vs_pT"), track.pt(), pdgbin);
              if (track.hasTPC()) {
                MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_vs_p"), track.p(), pdgbin);
                MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_vs_pT"), track.pt(), pdgbin);
                if (track.hasTOF()) {
                  MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_TOF_vs_p"), track.p(), pdgbin);
                  MC_recon_reg.fill(HIST("hist_rec_ITS_TPC_TOF_vs_pT"), track.pt(), pdgbin);
                }
              }
            }
          }
        }

        // Skipping collisions without the generated collisions
        // Actually this should never happen, since we group per MC collision
        if (!collision.has_mcCollision()) {
          continue;
        } else {
          // skip generated collisions outside the allowed vtx-z range
          // putting this condition here avoids the particle loop a few lines below
          if (applyPvZCutGenColl) {
            const float genPvZ = mcCollision.posZ();
            if (genPvZ < -cfgCutVertex || genPvZ > cfgCutVertex)
              continue;
          }
        }

        MC_gen_reg.fill(HIST("histGenVtxMC_reco"), mcCollision.posZ());
        MC_gen_reg.fill(HIST("histCentrality_reco"), centrality);

        /// only to fill denominator of ITS-TPC matched primary tracks only in MC events with at least 1 reco. vtx
        for (const auto& particle : groupedMcParticles) { // Particle loop

          /// require generated particle in acceptance
          if (!isInAcceptance(particle))
            continue;

          int pdgbin = -10;
          switch (particle.pdgCode()) {
            case +211:
              histPDG_gen_reco->AddBinContent(1);
              pdgbin = 0;
              break;
            case -211:
              histPDG_gen_reco->AddBinContent(2);
              pdgbin = 1;
              break;
            case +321:
              histPDG_gen_reco->AddBinContent(3);
              pdgbin = 2;
              break;
            case -321:
              histPDG_gen_reco->AddBinContent(4);
              pdgbin = 3;
              break;
            case +2212:
              histPDG_gen_reco->AddBinContent(5);
              pdgbin = 4;
              break;
            case -2212:
              histPDG_gen_reco->AddBinContent(6);
              pdgbin = 5;
              break;
            case +1000010020:
              histPDG_gen_reco->AddBinContent(7);
              pdgbin = 6;
              break;
            case -1000010020:
              histPDG_gen_reco->AddBinContent(8);
              pdgbin = 7;
              break;
            case +1000010030:
              histPDG_gen_reco->AddBinContent(9);
              pdgbin = 8;
              break;
            case -1000010030:
              histPDG_gen_reco->AddBinContent(10);
              pdgbin = 9;
              break;
            case +1000020030:
              histPDG_gen_reco->AddBinContent(11);
              pdgbin = 10;
              break;
            case -1000020030:
              histPDG_gen_reco->AddBinContent(12);
              pdgbin = 11;
              break;
            case +1000020040:
              histPDG_gen_reco->AddBinContent(13);
              pdgbin = 12;
              break;
            case -1000020040:
              histPDG_gen_reco->AddBinContent(14);
              pdgbin = 13;
              break;
            default:
              pdgbin = -10;
              continue;
              break;
          }
          MC_gen_reg.fill(HIST("histEta_reco"), particle.eta(), pdgbin);
          MC_gen_reg.fill(HIST("hist_gen_reco_p"), particle.p(), pdgbin);
          MC_gen_reg.fill(HIST("hist_gen_reco_pT"), particle.pt(), pdgbin);
        }
      } /// end loop over reconstructed collisions

      // skip generated collisions outside the allowed vtx-z range
      // putting this condition here avoids the particle loop a few lines below
      if (applyPvZCutGenColl) {
        const float genPvZ = mcCollision.posZ();
        if (genPvZ < -cfgCutVertex || genPvZ > cfgCutVertex) {
          continue;
        }
      }

      // Loop on particles to fill the denominator
      for (const auto& mcParticle : groupedMcParticles) {
        if (!isInAcceptance(mcParticle))
          continue;

        MC_gen_reg.fill(HIST("histGenVtxMC"), mcCollision.posZ());
        // MC_gen_reg.fill(HIST("histCentrality"), mcParticle.impactParameter());

        int pdgbin = -10;
        switch (mcParticle.pdgCode()) {
          case +211:
            histPDG_gen->AddBinContent(1);
            pdgbin = 0;
            break;
          case -211:
            histPDG_gen->AddBinContent(2);
            pdgbin = 1;
            break;
          case +321:
            histPDG_gen->AddBinContent(3);
            pdgbin = 2;
            break;
          case -321:
            histPDG_gen->AddBinContent(4);
            pdgbin = 3;
            break;
          case +2212:
            histPDG_gen->AddBinContent(5);
            pdgbin = 4;
            break;
          case -2212:
            histPDG_gen->AddBinContent(6);
            pdgbin = 5;
            break;
          case +1000010020:
            histPDG_gen->AddBinContent(7);
            pdgbin = 6;
            break;
          case -1000010020:
            histPDG_gen->AddBinContent(8);
            pdgbin = 7;
            break;
          case +1000010030:
            histPDG_gen->AddBinContent(9);
            pdgbin = 8;
            break;
          case -1000010030:
            histPDG_gen->AddBinContent(10);
            pdgbin = 9;
            break;
          case +1000020030:
            histPDG_gen->AddBinContent(11);
            pdgbin = 10;
            break;
          case -1000020030:
            histPDG_gen->AddBinContent(12);
            pdgbin = 11;
            break;
          case +1000020040:
            histPDG_gen->AddBinContent(13);
            pdgbin = 12;
            break;
          case -1000020040:
            histPDG_gen->AddBinContent(14);
            pdgbin = 13;
            break;
          default:
            pdgbin = -10;
            continue;
            break;
        }

        MC_gen_reg.fill(HIST("histPhi"), mcParticle.phi(), pdgbin);
        MC_gen_reg.fill(HIST("histEta"), mcParticle.eta(), pdgbin);
        MC_gen_reg.fill(HIST("histRapid"), mcParticle.y(), pdgbin);
        MC_gen_reg.fill(HIST("hist_gen_p"), mcParticle.p(), pdgbin);
        MC_gen_reg.fill(HIST("hist_gen_pT"), mcParticle.pt(), pdgbin);
      }
    } /// end loop over generated collisions
  }
  PROCESS_SWITCH(NucleiEfficiencyTask, processMC, "process generated MC", true);
};

//***********************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiEfficiencyTask>(cfgc, TaskName{"nuclei-efficiency-hist"})};
}
