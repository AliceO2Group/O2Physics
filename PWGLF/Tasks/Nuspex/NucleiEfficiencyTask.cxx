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

#include <cmath>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <vector>
#include <TF1.h>
#include <string>

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/HistogramRegistry.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "TPDGCode.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/StaticFor.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "PWGLF/DataModel/spectraTOF.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/inelGt.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "Common/Core/RecoDecay.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;
// using CollisionCandidates = o2::soa::Join<>;

struct NucleiEfficiencyTask {

  HistogramRegistry MC_gen_reg{"MC_particles_gen", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_gen_reg_cent{"MC_particles_gen_cent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg{"MC_particles_reco", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry MC_recon_reg_cent{"MC_particles_reco_cent", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> histPDG_gen{TH1F("PDG_gen", "PDG;PDG code", 18, 0.0, 18)};
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
    MC_gen_reg_cent.add("hist_gen_p_cent", "generated p distribution vs impact param", HistType::kTH3F, {pAxis, PDGBINNING, ImPaAxis});
    MC_gen_reg_cent.add("hist_gen_pT_cent", "generated p_{T} distribution vs impact param", HistType::kTH3F, {ptAxis, PDGBINNING, ImPaAxis});

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
    MC_recon_reg_cent.add("hist_rec_ITS_vs_p_cent", "ITS reconstructed p distribution vs centrality", HistType::kTH3F, {pAxis, PDGBINNING, centralityAxis});
    MC_recon_reg_cent.add("hist_rec_ITS_TPC_vs_p_cent", "ITS_TPC reconstructed p distribution vs centrality", HistType::kTH3F, {pAxis, PDGBINNING, centralityAxis});
    MC_recon_reg_cent.add("hist_rec_ITS_TPC_TOF_vs_p_cent", "ITS_TPC_TOF reconstructed p distribution vs centrality", HistType::kTH3F, {pAxis, PDGBINNING, centralityAxis});
    MC_recon_reg_cent.add("hist_rec_ITS_vs_pT_cent", "ITS reconstructed p_{T} distribution vs centrality", HistType::kTH3F, {ptAxis, PDGBINNING, centralityAxis});
    MC_recon_reg_cent.add("hist_rec_ITS_TPC_vs_pT_cent", "ITS_TPC reconstructed p_{T} distribution vs centrality", HistType::kTH3F, {ptAxis, PDGBINNING, centralityAxis});
    MC_recon_reg_cent.add("hist_rec_ITS_TPC_TOF_vs_pT_cent", "ITS_TPC_TOF reconstructed p_{T} distribution vs centrality", HistType::kTH3F, {ptAxis, PDGBINNING, centralityAxis});
  }

  // ************************ Configurables ***********************
  Configurable<bool> event_selection_MC_sel8{"event_selection_MC_sel8", true, "Enable sel8 event selection in MC processing"};
  Configurable<bool> y_cut_MC_gen{"y_cut_MC_gen", true, "Enable rapidity cut for generated MC"};
  Configurable<float> yMin_gen{"yMin_gen", -0.5, "Maximum rapidity (generated)"};
  Configurable<float> yMax_gen{"yMax_gen", 0.5, "Minimum rapidity (generated)"};
  Configurable<float> yMin_reco{"yMin_reco", -0.5, "Maximum rapidity (reconstructed)"};
  Configurable<float> yMax_reco{"yMax_reco", 0.5, "Minimum rapidity (reconstructed)"};
  Configurable<float> p_min{"p_min", 0.1f, "min track.pt()"};
  Configurable<float> p_max{"p_max", 1e+10f, "max track.pt()"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
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

  template <typename CollisionType>
  int getCentralityMC(CollisionType const& collision)
  {
    float multiplicity{0.f};
    int centrality = 0;
    float collMult{0.f};
    collMult = collision.numContrib();

    if (collMult > multiplicity) {
      centrality = collision.centFT0C();
      multiplicity = collMult;
    }

    return centrality;
  }

  //***********************************************************************************

  template <typename McCollisionType, typename McParticlesType>
  void process_MC_gen(const McCollisionType& mcCollision, const McParticlesType& mcParticles)
  {
    MC_gen_reg.fill(HIST("histGenVtxMC"), mcCollision.posZ());
    MC_gen_reg.fill(HIST("histCentrality"), mcCollision.impactParameter());

    for (const auto& MCparticle : mcParticles) {
      if (!MCparticle.isPhysicalPrimary())
        continue;
      if ((MCparticle.y() > yMax_gen || MCparticle.y() < yMin_gen) && y_cut_MC_gen)
        continue;
      if ((TMath::Abs(MCparticle.eta()) > cfgCutEta) && eta_cut_MC_gen)
        continue;

      int pdgbin = -10;
      switch (MCparticle.pdgCode()) {
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

      MC_gen_reg.fill(HIST("histPhi"), MCparticle.phi(), pdgbin);
      MC_gen_reg.fill(HIST("histEta"), MCparticle.eta(), pdgbin);
      MC_gen_reg.fill(HIST("histRapid"), MCparticle.y(), pdgbin);
      MC_gen_reg.fill(HIST("hist_gen_p"), MCparticle.p(), pdgbin);
      MC_gen_reg.fill(HIST("hist_gen_pT"), MCparticle.pt(), pdgbin);

      if (calc_cent) {
        MC_gen_reg_cent.fill(HIST("hist_gen_p_cent"), MCparticle.p(), pdgbin, mcCollision.impactParameter());
        MC_gen_reg_cent.fill(HIST("hist_gen_pT_cent"), MCparticle.pt(), pdgbin, mcCollision.impactParameter());
      }
    }
  }

  //***********************************************************************************

  template <typename CollisionType, typename TracksType, typename mcParticlesType>
  void process_MC_reco(const CollisionType& collision, const TracksType& tracks, const mcParticlesType& /*mcParticles*/)
  {

    int centrality = getCentralityMC(collision);
    if (event_selection_MC_sel8 && !collision.sel8())
      return;
    if (collision.posZ() > cfgCutVertex)
      return;
    MC_recon_reg.fill(HIST("histRecVtxMC"), collision.posZ());
    if (!isEventSelected(collision))
      return;
    if (centrality < minCentrality || centrality > maxCentrality)
      return;
    MC_recon_reg.fill(HIST("histCentrality"), centrality);

    for (auto& track : tracks) {
      const auto particle = track.mcParticle();
      if (!particle.isPhysicalPrimary())
        continue;
      TLorentzVector lorentzVector_particle_MC{};

      int pdgbin = -10;
      switch (particle.pdgCode()) {
        case +211:
          histPDG_reco->AddBinContent(1);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          pdgbin = 0;
          break;
        case -211:
          histPDG_reco->AddBinContent(2);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassPiPlus);
          pdgbin = 1;
          break;
        case +321:
          histPDG_reco->AddBinContent(3);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          pdgbin = 2;
          break;
        case -321:
          histPDG_reco->AddBinContent(4);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassKPlus);
          pdgbin = 3;
          break;
        case +2212:
          histPDG_reco->AddBinContent(5);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          pdgbin = 4;
          break;
        case -2212:
          histPDG_reco->AddBinContent(6);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassProton);
          pdgbin = 5;
          break;
        case +1000010020:
          histPDG_reco->AddBinContent(7);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          pdgbin = 6;
          break;
        case -1000010020:
          histPDG_reco->AddBinContent(8);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassDeuteron);
          pdgbin = 7;
          break;
        case +1000010030:
          histPDG_reco->AddBinContent(9);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          pdgbin = 8;
          break;
        case -1000010030:
          histPDG_reco->AddBinContent(10);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), constants::physics::MassTriton);
          pdgbin = 9;
          break;
        case +1000020030:
          histPDG_reco->AddBinContent(11);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          pdgbin = 10;
          break;
        case -1000020030:
          histPDG_reco->AddBinContent(12);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassHelium3);
          pdgbin = 11;
          break;
        case +1000020040:
          histPDG_reco->AddBinContent(13);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          pdgbin = 12;
          break;
        case -1000020040:
          histPDG_reco->AddBinContent(14);
          lorentzVector_particle_MC.SetPtEtaPhiM(track.pt() * 2.0, track.eta(), track.phi(), constants::physics::MassAlpha);
          pdgbin = 13;
          break;
        default:
          pdgbin = -10;
          continue;
          break;
      }

      if (lorentzVector_particle_MC.Rapidity() < yMin_reco || lorentzVector_particle_MC.Rapidity() > yMax_reco)
        continue;

      float TPCnumberClsFound = track.tpcNClsFound();
      float TPC_nCls_Crossed_Rows = track.tpcNClsCrossedRows();
      float RatioCrossedRowsOverFindableTPC = track.tpcCrossedRowsOverFindableCls();
      float Chi2perClusterTPC = track.tpcChi2NCl();
      float Chi2perClusterITS = track.itsChi2NCl();

      bool insideDCAxy = (std::abs(track.dcaXY()) <= (maxDcaXYFactor.value * (0.0105f + 0.0350f / pow(track.pt(), 1.1f))));

      if (!(insideDCAxy) || TMath::Abs(track.dcaZ()) > maxDCA_Z || TPCnumberClsFound < minTPCnClsFound || TPC_nCls_Crossed_Rows < minNCrossedRowsTPC || RatioCrossedRowsOverFindableTPC < minRatioCrossedRowsTPC || RatioCrossedRowsOverFindableTPC > maxRatioCrossedRowsTPC || Chi2perClusterTPC > maxChi2PerClusterTPC || Chi2perClusterTPC < minChi2PerClusterTPC || Chi2perClusterITS > maxChi2PerClusterITS || !(track.passedTPCRefit()) || !(track.passedITSRefit()) || (track.itsNClsInnerBarrel()) < minReqClusterITSib || (track.itsNCls()) < minReqClusterITS || track.pt() < p_min || track.pt() > p_max)
        continue;
      if ((requireITS && !(track.hasITS())) || (requireTPC && !(track.hasTPC())))
        continue;
      if (requireGoldenChi2 && !(track.passedGoldenChi2()))
        continue;

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
      if (calc_cent) {
        if ((particle.pdgCode() == 1000020030) || (particle.pdgCode() == -1000020030) || (particle.pdgCode() == 1000020040) || (particle.pdgCode() == -1000020040)) {
          if (track.hasITS()) {
            MC_recon_reg_cent.fill(HIST("hist_rec_ITS_vs_p_cent"), track.p() * 2, pdgbin, collision.centFT0C());
            MC_recon_reg_cent.fill(HIST("hist_rec_ITS_vs_pT_cent"), track.pt() * 2, pdgbin, collision.centFT0C());
            if (track.hasTPC()) {
              MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_vs_p_cent"), track.p() * 2, pdgbin, collision.centFT0C());
              MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_vs_pT_cent"), track.pt() * 2, pdgbin, collision.centFT0C());
              if (track.hasTOF()) {
                MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_TOF_vs_p_cent"), track.p() * 2, pdgbin, collision.centFT0C());
                MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_TOF_vs_pT_cent"), track.pt() * 2, pdgbin, collision.centFT0C());
              }
            }
          }
        } else {
          if (track.hasITS()) {
            MC_recon_reg_cent.fill(HIST("hist_rec_ITS_vs_p_cent"), track.p(), pdgbin, collision.centFT0C());
            MC_recon_reg_cent.fill(HIST("hist_rec_ITS_vs_pT_cent"), track.pt(), pdgbin, collision.centFT0C());
            if (track.hasTPC()) {
              MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_vs_p_cent"), track.p(), pdgbin, collision.centFT0C());
              MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_vs_pT_cent"), track.pt(), pdgbin, collision.centFT0C());
              if (track.hasTOF()) {
                MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_TOF_vs_p_cent"), track.p(), pdgbin, collision.centFT0C());
                MC_recon_reg_cent.fill(HIST("hist_rec_ITS_TPC_TOF_vs_pT_cent"), track.pt(), pdgbin, collision.centFT0C());
              }
            }
          }
        }
      }
    }
  }

  //***********************************************************************************

  void processMCgen(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    process_MC_gen(mcCollision, mcParticles);
  }
  PROCESS_SWITCH(NucleiEfficiencyTask, processMCgen, "process generated MC", true);

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgCutVertex);
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta && requireGlobalTrackWoDCAInFilter());

  void processMCreco(soa::Filtered<soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0Cs>>::iterator const& collision,
                     soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::TrackSelectionExtension>> const& tracks,
                     aod::McParticles const& mcParticles)
  {
    process_MC_reco(collision, tracks, mcParticles);
  }
  PROCESS_SWITCH(NucleiEfficiencyTask, processMCreco, "process reconstructed MC", false);
};

//***********************************************************************************

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleiEfficiencyTask>(cfgc, TaskName{"nuclei-efficiency-hist"})};
}
