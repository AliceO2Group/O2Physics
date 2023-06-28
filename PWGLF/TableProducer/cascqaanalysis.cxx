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
/// \brief QA task for V0 analysis using derived data
///
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "PWGLF/DataModel/cascqaanalysis.h"
#include "TRandom.h"
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TParticle.h>
#include <TParticlePDG.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

struct cascqaanalysis {

  // Produces
  Produces<aod::MyCascades> mycascades;

  HistogramRegistry registry{"registry"};

  AxisSpec ptAxis = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
  AxisSpec rapidityAxis = {200, -2.0f, 2.0f, "y"};
  ConfigurableAxis centAxis{"FT0M",
                            {VARIABLE_WIDTH, 0., 0.01, 0.05, 0.1, 0.5, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100},
                            "FT0M (%)"};

  void init(InitContext const&)
  {
    TString hCandidateCounterLabels[5] = {"All candidates", "v0data exists", "passed topo cuts", "has associated MC particle", "associated with Xi(Omega)"};
    TString hNEventsMCLabels[4] = {"All", "z vrtx", "INEL>0", "Associated with rec. collision"};
    TString hNEventsLabels[4] = {"All", "sel8", "z vrtx", "INEL>0"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1F, {{4, 0.f, 4.f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hNEvents"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEvents"))->GetXaxis()->SetBinLabel(n, hNEventsLabels[n - 1]);
    }
    registry.add("hNAssocCollisions", "hNAssocCollisions", {HistType::kTH1F, {{5, -0.5f, 4.5f}}});
    registry.add("hNContributorsCorrelation", "hNContributorsCorrelation", {HistType::kTH2F, {{250, -0.5f, 249.5f, "Secondary Contributor"}, {250, -0.5f, 249.5f, "Main Contributor"}}});
    registry.add("hZCollision", "hZCollision", {HistType::kTH1F, {{200, -20.f, 20.f}}});
    registry.add("hZCollisionGen", "hZCollisionGen", {HistType::kTH1F, {{200, -20.f, 20.f}}});
    registry.add("hCentFT0M", "hCentFT0M", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hCentFV0A", "hCentFV0A", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hPtXiPlusTrue", "hPtXiPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtXiMinusTrue", "hPtXiMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtOmegaPlusTrue", "hPtOmegaPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtOmegaMinusTrue", "hPtOmegaMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtXiPlusTrueAssoiciatedWithSelColl", "hPtXiPlusTrueAssoiciatedWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtXiMinusTrueAssoiciatedWithSelColl", "hPtXiMinusTrueAssoiciatedWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtOmegaPlusTrueAssoiciatedWithSelColl", "hPtOmegaPlusTrueAssoiciatedWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});
    registry.add("hPtOmegaMinusTrueAssoiciatedWithSelColl", "hPtOmegaMinusTrueAssoiciatedWithSelColl", {HistType::kTH3F, {ptAxis, rapidityAxis, centAxis}});

    registry.add("hNEventsMC", "hNEventsMC", {HistType::kTH1F, {{4, 0.0f, 4.0f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hNEventsMC"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEventsMC"))->GetXaxis()->SetBinLabel(n, hNEventsMCLabels[n - 1]);
    }

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hCandidateCounter"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hCandidateCounter"))->GetXaxis()->SetBinLabel(n, hCandidateCounterLabels[n - 1]);
    }

    AxisSpec allTracks = {2000, 0, 2000, "N_{all tracks}"};
    AxisSpec secondaryTracks = {2000, 0, 2000, "N_{secondary tracks}"};
    registry.add("hINELgt0PrimariesSelection", "hINELgt0PrimariesSelection", {HistType::kTH2F, {allTracks, secondaryTracks}});
    registry.add("hDCAz_BefCut", "hDCAz_BefCut", HistType::kTH2F, {{400, -0.2, 0.2, "DCAz"}, {150, 0.0, 15.0, "p_{T} (GeV/c)"}});
    registry.add("hDCAz_AfterCut", "hDCAz_AfterCut", HistType::kTH2F, {{400, -0.2, 0.2, "DCAz"}, {150, 0.0, 15.0, "p_{T} (GeV/c)"}});
    registry.add("hDCAxy_BefCut", "hDCAxy_BefCut", HistType::kTH2F, {{400, -0.2, 0.2, "DCAxy"}, {150, 0.0, 15.0, "p_{T} (GeV/c)"}});
    registry.add("hDCAxy_AfterCut", "hDCAxy_AfterCut", HistType::kTH2F, {{400, -0.2, 0.2, "DCAxy"}, {150, 0.0, 15.0, "p_{T} (GeV/c)"}});
    registry.add("hNchMultFT0M", "hNchMultFT0M", HistType::kTH2F, {{300, 0.0f, 300.0f, "N_{ch}"}, {10000, 0.f, 10000.f, "FT0M signal"}});
    registry.add("hNchMultFV0A", "hNchMultFV0A", HistType::kTH2F, {{300, 0.0f, 300.0f, "N_{ch}"}, {15000, 0.f, 15000.f, "FV0A signal"}});
  }

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 20.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};
  Configurable<bool> INELgt0{"INELgt0", 1, "Apply INEL>0 selection"};

  // Selection criteria
  Configurable<float> scalefactor{"scalefactor", 1.0, "Scaling factor"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcacascdau{"dcacascdau", 2.0, "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.0, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.0, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.0, "DCA Bach To PV"};
  Configurable<float> v0radius{"v0radius", 0.0, "V0 Radius"};
  Configurable<float> cascradius{"cascradius", 0.0, "Casc Radius"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};

  Configurable<float> maxDCANsigmaScaling{"maxDCANsigmaScaling", 1.0f, "N of 7*sigma scaling factor for DCA to select primaries"};
  Configurable<float> DCASigma{"DCASigma", 0.004f, "7*sigma for DCA"};
  Configurable<float> DCAPtScaling{"DCAPtScaling", 0.013f, "pt scaling for DCA"};
  Configurable<float> maxDCAz{"maxDCAz", 0.5f, "DCA z cut to select primaries"};

  TRandom* fRand = new TRandom();

  Filter preFilter =
    (nabs(aod::cascdata::dcapostopv) > dcapostopv &&
     nabs(aod::cascdata::dcanegtopv) > dcanegtopv &&
     nabs(aod::cascdata::dcabachtopv) > dcabachtopv &&
     aod::cascdata::dcaV0daughters < dcav0dau &&
     aod::cascdata::dcacascdaughters < dcacascdau);

  template <class TCascTracksTo, typename TCascade>
  bool AcceptCascCandidate(TCascade const& cascCand, float const& pvx, float const& pvy, float const& pvz)
  {
    // Access daughter tracks
    auto v0index = cascCand.template v0_as<o2::aod::V0sLinked>();
    auto v0 = v0index.v0Data();
    auto posdau = v0.template posTrack_as<TCascTracksTo>();
    auto negdau = v0.template negTrack_as<TCascTracksTo>();
    auto bachelor = cascCand.template bachelor_as<TCascTracksTo>();

    // Basic set of selections
    if (cascCand.cascradius() > cascradius &&
        v0.v0radius() > v0radius &&
        cascCand.casccosPA(pvx, pvy, pvz) > casccospa &&
        cascCand.v0cosPA(pvx, pvy, pvz) > v0cospa &&
        TMath::Abs(posdau.eta()) < etadau &&
        TMath::Abs(negdau.eta()) < etadau &&
        TMath::Abs(bachelor.eta()) < etadau) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TTracks>
  bool AcceptEvent(TCollision const& collision, TTracks const& tracks, bool isFillEventSelectionQA)
  {
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 0.5);
    }
    // Event selection if required
    if (sel8 && !collision.sel8()) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 1.5);
    }

    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 2.5);
    }

    if (INELgt0 && !isINELgt0(tracks, isFillEventSelectionQA)) {
      return false;
    }
    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 3.5);
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hZCollision"), collision.posZ());
      registry.fill(HIST("hCentFT0M"), collision.centFT0M());
      registry.fill(HIST("hCentFV0A"), collision.centFV0A());
    }
    return true;
  }

  template <typename TTrack>
  bool isPrimaryTrack(TTrack track)
  {
    return (TMath::Abs(track.dcaXY()) < (maxDCANsigmaScaling * (DCASigma + DCAPtScaling / track.pt()))) && (TMath::Abs(track.dcaZ()) < maxDCAz);
  }

  template <typename TTracks>
  bool isINELgt0(TTracks tracks, bool isFillEventSelectionQA)
  {
    // INEL > 0 (at least 1 charged track in |eta| < 1.0)
    std::vector<float> TracksEta(tracks.size());
    int nTracks = 0;
    int nRejTracks = 0;

    for (const auto& track : tracks) {
      registry.fill(HIST("hDCAxy_BefCut"), track.dcaXY(), track.pt());
      registry.fill(HIST("hDCAz_BefCut"), track.dcaZ(), track.pt());

      if (!isPrimaryTrack(track)) {
        nRejTracks++;
        continue; // consider only primaries
      }
      TracksEta[nTracks++] = track.eta();

      registry.fill(HIST("hDCAxy_AfterCut"), track.dcaXY(), track.pt());
      registry.fill(HIST("hDCAz_AfterCut"), track.dcaZ(), track.pt());
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hINELgt0PrimariesSelection"), tracks.size(), nRejTracks);
    }

    auto etaConditionFunc = [](float elem) {
      return TMath::Abs(elem) < 1.0;
    };

    if (std::any_of(TracksEta.begin(), TracksEta.end(), etaConditionFunc)) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TMcParticles>
  bool isINELgt0mc(TMcParticles particles)
  {
    // INEL > 0 (at least 1 charged particle in |eta| < 1.0)
    typedef struct EtaCharge {
      double eta;
      int charge;
    } EtaCharge;
    EtaCharge etaCharge;
    std::vector<EtaCharge> ParticlesEtaAndCharge(particles.size());
    unsigned int nParticles = 0;
    for (const auto& particle : particles) {
      if (particle.isPhysicalPrimary() == 0)
        continue;           // consider only primaries
      etaCharge = {999, 0}; // refresh init. for safety
      TParticlePDG* p = TDatabasePDG::Instance()->GetParticle(particle.pdgCode());
      if (!p) {
        switch (std::to_string(particle.pdgCode()).length()) {
          case 10: // nuclei
          {
            etaCharge = {particle.eta(), static_cast<int>(particle.pdgCode() / 10000 % 1000)};
            ParticlesEtaAndCharge[nParticles++] = etaCharge;
            break;
          }
          default:
            break;
        }
      } else {
        etaCharge = {particle.eta(), static_cast<int>(p->Charge())};
        ParticlesEtaAndCharge[nParticles++] = etaCharge;
      }
    }

    ParticlesEtaAndCharge.resize(nParticles);

    auto etaChargeConditionFunc = [](EtaCharge elem) {
      return ((TMath::Abs(elem.eta) < 1.0) && (TMath::Abs(elem.charge) > 0.001));
    };

    if (std::any_of(ParticlesEtaAndCharge.begin(), ParticlesEtaAndCharge.end(), etaChargeConditionFunc)) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TCollision, typename TTracks>
  void fillMultHisto(TCollision const& collision, TTracks const& tracks)
  {
    double Nch = 0;
    for (const auto& track : tracks) {
      if (TMath::Abs(track.eta()) > 0.5)
        continue;
      if (!isPrimaryTrack(track))
        continue;
      Nch++;
    }
    registry.fill(HIST("hNchMultFT0M"), Nch, collision.multFT0A() + collision.multFT0C());
    registry.fill(HIST("hNchMultFV0A"), Nch, collision.multFV0A());
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                   soa::Filtered<aod::CascDataExt> const& Cascades,
                   aod::V0sLinked const&,
                   aod::V0Datas const&,
                   DauTracks const& Tracks)
  {
    if (!AcceptEvent(collision, Tracks, 1)) {
      return;
    }

    float lEventScale = scalefactor;

    for (const auto& casc : Cascades) {              // loop over Cascades
      registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates

      // Access daughter tracks
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        return; // skip those cascades for which V0 doesn't exist
      }
      registry.fill(HIST("hCandidateCounter"), 1.5); // v0data exists

      auto v0 = v0index.v0Data();
      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();
      auto bachelor = casc.bachelor_as<DauTracks>();

      // ITS N hits
      int posITSNhits = 0, negITSNhits = 0, bachITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (posdau.itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (negdau.itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
        if (bachelor.itsClusterMap() & (1 << i)) {
          bachITSNhits++;
        }
      }

      // c x tau
      float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      //
      float ctauXi = RecoDecay::getMassPDG(3312) * cascpos / (cascptotmom + 1e-13);
      float ctauOmega = RecoDecay::getMassPDG(3334) * cascpos / (cascptotmom + 1e-13);

      if (AcceptCascCandidate<DauTracks>(casc, collision.posX(), collision.posY(), collision.posZ())) {
        registry.fill(HIST("hCandidateCounter"), 2.5); // passed topo cuts
        // Fill table
        if (fRand->Rndm() < lEventScale) {
          mycascades(casc.globalIndex(), collision.posZ(), collision.centFT0M(), collision.centFV0A(), casc.sign(), casc.pt(), casc.yXi(), casc.yOmega(), casc.eta(),
                     casc.mXi(), casc.mOmega(), casc.mLambda(), casc.cascradius(), casc.v0radius(),
                     casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(), casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     posdau.eta(), negdau.eta(), bachelor.eta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negdau.tpcNSigmaPr(), posdau.tpcNSigmaPr(), negdau.tpcNSigmaPi(), posdau.tpcNSigmaPi(), bachelor.tpcNSigmaPi(), bachelor.tpcNSigmaKa(),
                     negdau.tofNSigmaPr(), posdau.tofNSigmaPr(), negdau.tofNSigmaPi(), posdau.tofNSigmaPi(), bachelor.tofNSigmaPi(), bachelor.tofNSigmaKa(),
                     posdau.tpcNClsFound(), negdau.tpcNClsFound(), bachelor.tpcNClsFound(),
                     posdau.hasTOF(), negdau.hasTOF(), bachelor.hasTOF(),
                     posdau.pt(), negdau.pt(), bachelor.pt(), -1, -1, casc.bachBaryonCosPA(), casc.bachBaryonDCAxyToPV());
        }
      }
    }
  }

  PROCESS_SWITCH(cascqaanalysis, processData, "Process Run 3 data", true);

  void processMCrec(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                    soa::Filtered<LabeledCascades> const& Cascades,
                    aod::V0sLinked const&,
                    aod::V0Datas const&,
                    DauTracks const& Tracks,
                    aod::McParticles const&)
  {
    if (!AcceptEvent(collision, Tracks, 1)) {
      return;
    }

    fillMultHisto(collision, Tracks);

    float lEventScale = scalefactor;

    for (const auto& casc : Cascades) {              // loop over Cascades
      registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates
      // Access daughter tracks
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        return; // skip those cascades for which V0 doesn't exist
      }

      registry.fill(HIST("hCandidateCounter"), 1.5); // v0data exists

      auto v0 = v0index.v0Data();
      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();
      auto bachelor = casc.bachelor_as<DauTracks>();

      // ITS N hits
      int posITSNhits = 0, negITSNhits = 0, bachITSNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (posdau.itsClusterMap() & (1 << i)) {
          posITSNhits++;
        }
        if (negdau.itsClusterMap() & (1 << i)) {
          negITSNhits++;
        }
        if (bachelor.itsClusterMap() & (1 << i)) {
          bachITSNhits++;
        }
      }

      // c x tau
      float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      //
      float ctauXi = RecoDecay::getMassPDG(3312) * cascpos / (cascptotmom + 1e-13);
      float ctauOmega = RecoDecay::getMassPDG(3334) * cascpos / (cascptotmom + 1e-13);

      if (AcceptCascCandidate<DauTracks>(casc, collision.posX(), collision.posY(), collision.posZ())) {
        registry.fill(HIST("hCandidateCounter"), 2.5); // passed topo cuts
        // Check mc association
        float lPDG = -1;
        float isPrimary = -1;
        if (casc.has_mcParticle()) {
          registry.fill(HIST("hCandidateCounter"), 3.5); // has associated MC particle
          auto cascmc = casc.mcParticle();
          if (TMath::Abs(cascmc.pdgCode()) == 3312 || TMath::Abs(cascmc.pdgCode()) == 3334) {
            registry.fill(HIST("hCandidateCounter"), 4.5); // associated with Xi or Omega
            lPDG = cascmc.pdgCode();
            isPrimary = cascmc.isPhysicalPrimary() ? 1 : 0;
          }
        }
        // Fill table
        if (fRand->Rndm() < lEventScale) {
          mycascades(casc.globalIndex(), collision.posZ(), collision.multFT0A() + collision.multFT0C(), collision.multFV0A(), casc.sign(), casc.pt(), casc.yXi(), casc.yOmega(), casc.eta(),
                     casc.mXi(), casc.mOmega(), casc.mLambda(), casc.cascradius(), casc.v0radius(),
                     casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(), casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     posdau.eta(), negdau.eta(), bachelor.eta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negdau.tpcNSigmaPr(), posdau.tpcNSigmaPr(), negdau.tpcNSigmaPi(), posdau.tpcNSigmaPi(), bachelor.tpcNSigmaPi(), bachelor.tpcNSigmaKa(),
                     negdau.tofNSigmaPr(), posdau.tofNSigmaPr(), negdau.tofNSigmaPi(), posdau.tofNSigmaPi(), bachelor.tofNSigmaPi(), bachelor.tofNSigmaKa(),
                     posdau.tpcNClsFound(), negdau.tpcNClsFound(), bachelor.tpcNClsFound(),
                     posdau.hasTOF(), negdau.hasTOF(), bachelor.hasTOF(),
                     posdau.pt(), negdau.pt(), bachelor.pt(), lPDG, isPrimary, casc.bachBaryonCosPA(), casc.bachBaryonDCAxyToPV());
        }
      }
    }
  }

  PROCESS_SWITCH(cascqaanalysis, processMCrec, "Process Run 3 mc, reconstructed", false);

  void processMCgen(aod::McCollision const& mcCollision,
                    aod::McParticles const& mcParticles,
                    const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>>& collisions,
                    DauTracks const& Tracks)
  {
    // All generated collisions
    registry.fill(HIST("hNEventsMC"), 0.5);

    // Generated with accepted z vertex
    if (TMath::Abs(mcCollision.posZ()) > cutzvertex) {
      return;
    }
    registry.fill(HIST("hNEventsMC"), 1.5);

    // Generated collision is INEL>=0
    if (INELgt0 && !isINELgt0mc(mcParticles)) {
      return;
    }
    registry.fill(HIST("hNEventsMC"), 2.5);

    registry.fill(HIST("hZCollisionGen"), mcCollision.posZ());

    // Histos of generated cascades from generated events with accepted z vrtx + INEL>0 (for signal loss correction)
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() == 0)
        continue; // Consider only primaries
      if (mcParticle.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlusTrue"), mcParticle.pt(), mcParticle.y(), 0); // MB will be used correction
      }
      if (mcParticle.pdgCode() == 3312) {
        registry.fill(HIST("hPtXiMinusTrue"), mcParticle.pt(), mcParticle.y(), 0);
      }
      if (mcParticle.pdgCode() == -3334) {
        registry.fill(HIST("hPtOmegaPlusTrue"), mcParticle.pt(), mcParticle.y(), 0);
      }
      if (mcParticle.pdgCode() == 3334) {
        registry.fill(HIST("hPtOmegaMinusTrue"), mcParticle.pt(), mcParticle.y(), 0);
      }
    }

    std::vector<int64_t> SelectedEvents(collisions.size());
    std::vector<int64_t> NumberOfContributors;
    int nevts = 0;
    int nAssocColl = 0;
    for (const auto& collision : collisions) {
      if (!AcceptEvent(collision, Tracks, 0)) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
      if (collision.mcCollision_as<aod::McCollisions>().globalIndex() == mcCollision.globalIndex()) {
        nAssocColl++;
        NumberOfContributors.push_back(collision.numContrib());
      }
    }
    SelectedEvents.resize(nevts);
    registry.fill(HIST("hNAssocCollisions"), nAssocColl);
    if (NumberOfContributors.size() == 2) {
      std::sort(NumberOfContributors.begin(), NumberOfContributors.end());
      registry.fill(HIST("hNContributorsCorrelation"), NumberOfContributors[0], NumberOfContributors[1]);
    }

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end(); // at least 1 selected reconstructed event has the same global index as mcCollision

    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed event passes the selection
      return;
    }

    registry.fill(HIST("hNEventsMC"), 3.5);

    // Histos of generated cascades from generated events with good z vrtx + INEL>0 + associated to the accepted reconstructed event (for signal loss + efficiency x acceptance correction)
    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() == 0)
        continue; // Consider only primaries
      if (mcParticle.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlusTrueAssoiciatedWithSelColl"), mcParticle.pt(), mcParticle.y(), 0); // MB will be used correction
      }
      if (mcParticle.pdgCode() == 3312) {
        registry.fill(HIST("hPtXiMinusTrueAssoiciatedWithSelColl"), mcParticle.pt(), mcParticle.y(), 0);
      }
      if (mcParticle.pdgCode() == -3334) {
        registry.fill(HIST("hPtOmegaPlusTrueAssoiciatedWithSelColl"), mcParticle.pt(), mcParticle.y(), 0);
      }
      if (mcParticle.pdgCode() == 3334) {
        registry.fill(HIST("hPtOmegaMinusTrueAssoiciatedWithSelColl"), mcParticle.pt(), mcParticle.y(), 0);
      }
    }
  }
  PROCESS_SWITCH(cascqaanalysis, processMCgen, "Process Run 3 mc, genereated", false);
};

struct myCascades {

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    TString PGDlabels[3] = {"Unknown", "3312", "3334"};
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hMassXi", "hMassXi", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}});
    registry.add("hMassOmega", "hMassOmega", {HistType::kTH1F, {{3000, 0.0f, 3.0f}}});
    registry.add("hCascRadius", "hCascRadius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
    registry.add("hV0Radius", "hV0Radius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
    registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hCtauXi", "hCtauXi", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauOmega", "hCtauOmega", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hTPCNSigmaPosPi", "hTPCNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTPCNSigmaNegPi", "hTPCNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTPCNSigmaPosPr", "hTPCNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTPCNSigmaNegPr", "hTPCNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTPCNSigmaBachPi", "hTPCNSigmaBachPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTOFNSigmaPosPi", "hTOFNSigmaPosPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTOFNSigmaNegPi", "hTOFNSigmaNegPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTOFNSigmaPosPr", "hTOFNSigmaPosPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTOFNSigmaNegPr", "hTOFNSigmaNegPr", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hTOFNSigmaBachPi", "hTOFNSigmaBachPi", {HistType::kTH1F, {{100, -10.0f, 10.0f}}});
    registry.add("hPosITSHits", "hPosITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("hNegITSHits", "hNegITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("hBachITSHits", "hBachITSHits", {HistType::kTH1F, {{8, -0.5f, 7.5f}}});
    registry.add("hIsPrimary", "hIsPrimary", {HistType::kTH1F, {{3, -1.5f, 1.5f}}});
    registry.add("hPDGcode", "hPDGcode", {HistType::kTH1F, {{3, -1.5f, 1.5f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hPDGcode"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hPDGcode"))->GetXaxis()->SetBinLabel(n, PGDlabels[n - 1]);
    }
    registry.add("hBachBaryonCosPA", "hBachBaryonCosPA", {HistType::kTH1F, {{100, 0.0f, 1.0f}}});
    registry.add("hBachBaryonDCAxyToPV", "hBachBaryonDCAxyToPV", {HistType::kTH1F, {{300, -3.0f, 3.0f}}});
  }

  void process(aod::MyCascades const& mycascades)
  {
    for (auto& candidate : mycascades) {

      registry.fill(HIST("hMassXi"), candidate.massxi());
      registry.fill(HIST("hMassOmega"), candidate.massomega());
      registry.fill(HIST("hPt"), candidate.pt());
      registry.fill(HIST("hCascRadius"), candidate.cascradius());
      registry.fill(HIST("hV0Radius"), candidate.v0radius());
      registry.fill(HIST("hCascCosPA"), candidate.casccospa());
      registry.fill(HIST("hV0CosPA"), candidate.v0cospa());
      registry.fill(HIST("hDCANegToPV"), candidate.dcanegtopv());
      registry.fill(HIST("hDCAPosToPV"), candidate.dcapostopv());
      registry.fill(HIST("hDCABachToPV"), candidate.dcabachtopv());
      registry.fill(HIST("hDCACascDaughters"), candidate.dcacascdaughters());
      registry.fill(HIST("hDCAV0Daughters"), candidate.dcav0daughters());
      registry.fill(HIST("hCtauXi"), candidate.ctauxi());
      registry.fill(HIST("hCtauOmega"), candidate.ctauomega());
      registry.fill(HIST("hTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
      registry.fill(HIST("hTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
      registry.fill(HIST("hTPCNSigmaPosPr"), candidate.ntpcsigmapospr());
      registry.fill(HIST("hTPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
      registry.fill(HIST("hTPCNSigmaBachPi"), candidate.ntpcsigmabachpi());
      registry.fill(HIST("hTOFNSigmaPosPi"), candidate.ntofsigmapospi());
      registry.fill(HIST("hTOFNSigmaNegPi"), candidate.ntofsigmanegpi());
      registry.fill(HIST("hTOFNSigmaPosPr"), candidate.ntofsigmapospr());
      registry.fill(HIST("hTOFNSigmaNegPr"), candidate.ntofsigmanegpr());
      registry.fill(HIST("hTOFNSigmaBachPi"), candidate.ntofsigmabachpi());
      registry.fill(HIST("hPosITSHits"), candidate.positshits());
      registry.fill(HIST("hNegITSHits"), candidate.negitshits());
      registry.fill(HIST("hBachITSHits"), candidate.bachitshits());
      registry.fill(HIST("hIsPrimary"), candidate.isPrimary());
      registry.fill(HIST("hBachBaryonCosPA"), candidate.bachBaryonCosPA());
      registry.fill(HIST("hBachBaryonDCAxyToPV"), candidate.bachBaryonDCAxyToPV());

      if (TMath::Abs(candidate.mcPdgCode()) == 3312 || TMath::Abs(candidate.mcPdgCode()) == 3334) {
        registry.fill(HIST("hPDGcode"), TMath::Abs(candidate.mcPdgCode()) == 3312 ? 0 : 1); // 0 if Xi, 1 if Omega
      } else {
        registry.fill(HIST("hPDGcode"), -1); // -1 if unknown
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascqaanalysis>(cfgc, TaskName{"lf-cascqaanalysis"}),
    adaptAnalysisTask<myCascades>(cfgc, TaskName{"lf-mycascades"})};
}
