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
#include "TRandom.h"
#include <TPDGCode.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;
using LabeledCascades = soa::Join<aod::CascDataExt, aod::McCascLabels>;

namespace o2::aod
{

namespace mycascades
{

DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(CollisionZ, zcoll, float);
DECLARE_SOA_COLUMN(MultFT0M, multFT0M, float);
DECLARE_SOA_COLUMN(MultFV0A, multFV0A, float);
DECLARE_SOA_COLUMN(Sign, sign, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(RapXi, rapxi, float);
DECLARE_SOA_COLUMN(RapOmega, rapomega, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(MassXi, massxi, float);
DECLARE_SOA_COLUMN(MassOmega, massomega, float);
DECLARE_SOA_COLUMN(MassLambdaDau, masslambdadau, float);
DECLARE_SOA_COLUMN(CascRadius, cascradius, float);
DECLARE_SOA_COLUMN(V0Radius, v0radius, float);
DECLARE_SOA_COLUMN(CascCosPA, casccospa, float);
DECLARE_SOA_COLUMN(V0CosPA, v0cospa, float);
DECLARE_SOA_COLUMN(DCAPosToPV, dcapostopv, float);
DECLARE_SOA_COLUMN(DCANegToPV, dcanegtopv, float);
DECLARE_SOA_COLUMN(DCABachToPV, dcabachtopv, float);
DECLARE_SOA_COLUMN(DCACascDaughters, dcacascdaughters, float);
DECLARE_SOA_COLUMN(DCAV0Daughters, dcav0daughters, float);
DECLARE_SOA_COLUMN(DCAV0ToPV, dcav0topv, float);
DECLARE_SOA_COLUMN(PosEta, poseta, float);
DECLARE_SOA_COLUMN(NegEta, negeta, float);
DECLARE_SOA_COLUMN(BachEta, bacheta, float);
DECLARE_SOA_COLUMN(PosITSHits, positshits, float);
DECLARE_SOA_COLUMN(NegITSHits, negitshits, float);
DECLARE_SOA_COLUMN(BachITSHits, bachitshits, float);
DECLARE_SOA_COLUMN(CtauXi, ctauxi, float);
DECLARE_SOA_COLUMN(CtauOmega, ctauomega, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachPi, ntpcsigmabachpi, float);
DECLARE_SOA_COLUMN(NTPCSigmaBachKa, ntpcsigmabachka, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachPi, ntofsigmabachpi, float);
DECLARE_SOA_COLUMN(NTOFSigmaBachKa, ntofsigmabachka, float);
DECLARE_SOA_COLUMN(PosNTPCClusters, posntpcscls, float);
DECLARE_SOA_COLUMN(NegNTPCClusters, negntpcscls, float);
DECLARE_SOA_COLUMN(BachNTPCClusters, bachntpcscls, float);
DECLARE_SOA_COLUMN(PosHasTOF, poshastof, float);
DECLARE_SOA_COLUMN(NegHasTOF, neghastof, float);
DECLARE_SOA_COLUMN(BachHasTOF, bachhastof, float);
DECLARE_SOA_COLUMN(PosPt, pospt, float);
DECLARE_SOA_COLUMN(NegPt, negpt, float);
DECLARE_SOA_COLUMN(BachPt, bachpt, float);
DECLARE_SOA_COLUMN(McPdgCode, mcPdgCode, float); // -1 unknown
DECLARE_SOA_COLUMN(IsPrimary, isPrimary, float); // -1 unknown, 0 not primary, 1 primary

} // namespace mycascades

DECLARE_SOA_TABLE(MyCascades, "AOD", "MYCASCADES", o2::soa::Index<>,
                  mycascades::CollisionId, mycascades::CollisionZ, mycascades::MultFT0M, mycascades::MultFV0A, mycascades::Sign, mycascades::Pt, mycascades::RapXi, mycascades::RapOmega, mycascades::Eta, mycascades::MassXi, mycascades::MassOmega, mycascades::MassLambdaDau, mycascades::CascRadius, mycascades::V0Radius,
                  mycascades::CascCosPA, mycascades::V0CosPA, mycascades::DCAPosToPV, mycascades::DCANegToPV,
                  mycascades::DCABachToPV, mycascades::DCACascDaughters, mycascades::DCAV0Daughters, mycascades::DCAV0ToPV, mycascades::PosEta, mycascades::NegEta,
                  mycascades::BachEta, mycascades::PosITSHits, mycascades::NegITSHits, mycascades::BachITSHits,
                  mycascades::CtauXi, mycascades::CtauOmega,
                  mycascades::NTPCSigmaNegPr, mycascades::NTPCSigmaPosPr, mycascades::NTPCSigmaNegPi, mycascades::NTPCSigmaPosPi, mycascades::NTPCSigmaBachPi, mycascades::NTPCSigmaBachKa,
                  mycascades::NTOFSigmaNegPr, mycascades::NTOFSigmaPosPr, mycascades::NTOFSigmaNegPi,
                  mycascades::NTOFSigmaPosPi, mycascades::NTOFSigmaBachPi, mycascades::NTOFSigmaBachKa,
                  mycascades::PosNTPCClusters, mycascades::NegNTPCClusters, mycascades::BachNTPCClusters,
                  mycascades::PosHasTOF, mycascades::NegHasTOF, mycascades::BachHasTOF,
                  mycascades::PosPt, mycascades::NegPt, mycascades::BachPt, mycascades::McPdgCode, mycascades::IsPrimary);

} // namespace o2::aod

struct cascqaanalysis {

  // Produces
  Produces<aod::MyCascades> mycascades;

  HistogramRegistry registry{"registry"};

  AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
  AxisSpec rapidityAxis = {200, -2.0f, 2.0f, "y"};
  AxisSpec centFT0MAxis = {100, 0.0f, 100.0f, "FT0M (%)"};

  void init(InitContext const&)
  {
    TString CandidateCounterLabels[5] = {"All candidates", "v0data exists", "passed topo cuts", "has associated MC particle", "associated with Xi(Omega)"};
    TString hNEventsMCLabels[2] = {"All", "Selected"};

    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{1, 0.f, 1.f}}});
    registry.add("hZCollision", "hZCollision", {HistType::kTH1F, {{200, -20.f, 20.f}}});
    registry.add("hCentFT0M", "hCentFT0M", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hCentFV0A", "hCentFV0A", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hPtXiPlusTrue", "hPtXiPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
    registry.add("hPtXiMinusTrue", "hPtXiMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
    registry.add("hPtOmegaPlusTrue", "hPtOmegaPlusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});
    registry.add("hPtOmegaMinusTrue", "hPtOmegaMinusTrue", {HistType::kTH3F, {ptAxis, rapidityAxis, centFT0MAxis}});

    registry.add("hNEventsMC", "hNEventsMC", {HistType::kTH1F, {{2, 0.0f, 2.0f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hNEventsMC"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hNEventsMC"))->GetXaxis()->SetBinLabel(n, hNEventsMCLabels[n - 1]);
    }

    registry.add("hCandidateCounter", "hCandidateCounter", {HistType::kTH1F, {{5, 0.0f, 5.0f}}});
    for (Int_t n = 1; n <= registry.get<TH1>(HIST("hCandidateCounter"))->GetNbinsX(); n++) {
      registry.get<TH1>(HIST("hCandidateCounter"))->GetXaxis()->SetBinLabel(n, CandidateCounterLabels[n - 1]);
    }
  }

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};

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
    registry.fill(HIST("hCandidateCounter"), 0.5); // all candidates
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

  template <typename TCollision>
  bool AcceptEvent(TCollision const& collision, bool isFillEventSelectionQA)
  {
    // Event selection if required
    if (sel8 && !collision.sel8()) {
      return false;
    }
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return false;
    }

    if (isFillEventSelectionQA) {
      registry.fill(HIST("hNEvents"), 0.5);
      registry.fill(HIST("hZCollision"), collision.posZ());
      registry.fill(HIST("hCentFT0M"), collision.centFT0M());
      registry.fill(HIST("hCentFV0A"), collision.centFV0A());
    }

    return true;
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                   soa::Filtered<aod::CascDataExt> const& Cascades,
                   aod::V0sLinked const&,
                   aod::V0Datas const&,
                   DauTracks const& tracks)
  {
    if (!AcceptEvent(collision, 1)) {
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
                     posdau.pt(), negdau.pt(), bachelor.pt(), -1, -1);
        }
      }
    }
  }

  PROCESS_SWITCH(cascqaanalysis, processData, "Process Run 3 data", true);

  void processMCrec(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision,
                    soa::Filtered<LabeledCascades> const& Cascades,
                    aod::V0sLinked const&,
                    aod::V0Datas const&,
                    DauTracks const& tracks,
                    aod::McParticles const&) // centrality-table is not ready for MC yet.
  {
    if (!AcceptEvent(collision, 1)) {
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
          mycascades(casc.globalIndex(), collision.posZ(), collision.centFT0M(), collision.centFV0A(), casc.sign(), casc.pt(), casc.yXi(), casc.yOmega(), casc.eta(),
                     casc.mXi(), casc.mOmega(), casc.mLambda(), casc.cascradius(), casc.v0radius(),
                     casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                     casc.dcapostopv(), casc.dcanegtopv(), casc.dcabachtopv(), casc.dcacascdaughters(), casc.dcaV0daughters(), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                     posdau.eta(), negdau.eta(), bachelor.eta(), posITSNhits, negITSNhits, bachITSNhits,
                     ctauXi, ctauOmega, negdau.tpcNSigmaPr(), posdau.tpcNSigmaPr(), negdau.tpcNSigmaPi(), posdau.tpcNSigmaPi(), bachelor.tpcNSigmaPi(), bachelor.tpcNSigmaKa(),
                     negdau.tofNSigmaPr(), posdau.tofNSigmaPr(), negdau.tofNSigmaPi(), posdau.tofNSigmaPi(), bachelor.tofNSigmaPi(), bachelor.tofNSigmaKa(),
                     posdau.tpcNClsFound(), negdau.tpcNClsFound(), bachelor.tpcNClsFound(),
                     posdau.hasTOF(), negdau.hasTOF(), bachelor.hasTOF(),
                     posdau.pt(), negdau.pt(), bachelor.pt(), lPDG, isPrimary);
        }
      }
    }
  }

  PROCESS_SWITCH(cascqaanalysis, processMCrec, "Process Run 3 mc, reconstructed", false);

  void processMCgen(aod::McCollision const& mcCollision,
                    aod::McParticles const& mcParticles,
                    const soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>>& collisions)
  {
    registry.fill(HIST("hNEventsMC"), 0.5);

    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      if (!AcceptEvent(collision, 0)) {
        continue;
      }
      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);

    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end(); // at least 1 selected reconstructed event has the same global index as mcCollision

    if (!evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }

    registry.fill(HIST("hNEventsMC"), 1.5);

    for (const auto& mcParticle : mcParticles) {
      if (mcParticle.isPhysicalPrimary() == 0)
        continue; // Consider only primaries
      if (mcParticle.pdgCode() == -3312) {
        registry.fill(HIST("hPtXiPlusTrue"), mcParticle.pt(), mcParticle.y(), 0); // MB will be used for the efficiency
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
