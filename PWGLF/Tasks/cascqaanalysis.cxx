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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;
using DauTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa, aod::pidTOFPi, aod::pidTOFPr, aod::pidTOFKa>;

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
                  mycascades::PosPt, mycascades::NegPt, mycascades::BachPt);

} // namespace o2::aod

struct cascqaanalysis {

  // Produces
  Produces<aod::MyCascades> mycascades;

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
    registry.add("hNEvents", "hNEvents", {HistType::kTH1I, {{1, 0.f, 1.f}}});
    registry.add("hZCollision", "hZCollision", {HistType::kTH1F, {{200, -20.f, 20.f}}});
    registry.add("hCentFT0M", "hCentFT0M", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
    registry.add("hCentFV0A", "hCentFV0A", {HistType::kTH1F, {{1000, 0.f, 100.f}}});
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
    nabs(aod::cascdata::dcapostopv) > dcapostopv&& nabs(aod::cascdata::dcanegtopv) > dcanegtopv&& nabs(aod::cascdata::dcabachtopv) > dcabachtopv&& aod::cascdata::dcaV0daughters < dcav0dau&& aod::cascdata::dcacascdaughters < dcacascdau;

  void process(soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms, aod::CentFV0As>::iterator const& collision, soa::Filtered<aod::CascDataExt> const& Cascades, aod::V0sLinked const&, aod::V0Datas const&, DauTracks const& tracks)
  {
    // Event selection
    if (sel8 && !collision.sel8()) {
      return;
    }
    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      return;
    }

    registry.fill(HIST("hNEvents"), 0.5);
    registry.fill(HIST("hZCollision"), collision.posZ());
    registry.fill(HIST("hCentFT0M"), collision.centFT0M());
    registry.fill(HIST("hCentFV0A"), collision.centFV0A());

    float lEventScale = scalefactor;

    for (auto& casc : Cascades) { // loop over Cascades

      // Access daughter tracks
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        return; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data();
      auto posdau = v0.posTrack_as<DauTracks>();
      auto negdau = v0.negTrack_as<DauTracks>();
      auto bachelor = casc.bachelor_as<DauTracks>();

      // c x tau
      float cascpos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      float cascptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      //
      float ctauXi = RecoDecay::getMassPDG(3312) * cascpos / (cascptotmom + 1e-13);
      float ctauOmega = RecoDecay::getMassPDG(3334) * cascpos / (cascptotmom + 1e-13);

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

      // Basic set of selections
      if (casc.cascradius() > cascradius && v0.v0radius() > v0radius &&
          casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospa && casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) > v0cospa &&
          TMath::Abs(posdau.eta()) < etadau && TMath::Abs(negdau.eta()) < etadau && TMath::Abs(bachelor.eta()) < etadau) {

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
                     posdau.pt(), negdau.pt(), bachelor.pt());
        }
      }
    }
  }
};

struct myCascades {

  HistogramRegistry registry{"registry"};

  void init(InitContext const&)
  {
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
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascqaanalysis>(cfgc, TaskName{"lf-cascqaanalysis"}),
    adaptAnalysisTask<myCascades>(cfgc, TaskName{"lf-mycascades"})};
}
