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
/// \brief this is a starting point for the third session of the tutorial
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

namespace o2::aod {

  namespace mycascades {

  DECLARE_SOA_INDEX_COLUMN(Collision, collision);
  DECLARE_SOA_COLUMN(Sign, sign, float);
  DECLARE_SOA_COLUMN(Pt, pt, float);
  DECLARE_SOA_COLUMN(RapXi, rapxi, float);
  DECLARE_SOA_COLUMN(RapOmega, rapomega, float);
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
  DECLARE_SOA_COLUMN(CtauXiMinus, ctauximinus, float);
  DECLARE_SOA_COLUMN(CtauXiPlus, ctauxiplus, float);
  DECLARE_SOA_COLUMN(CtauOmegaMinus, ctauomegaminus, float);
  DECLARE_SOA_COLUMN(CtauOmegaPlus, ctauomegaplus, float);
  DECLARE_SOA_COLUMN(NTPCSigmaNegPr, ntpcsigmanegpr, float);
  DECLARE_SOA_COLUMN(NTPCSigmaPosPr, ntpcsigmapospr, float);
  DECLARE_SOA_COLUMN(NTPCSigmaNegPi, ntpcsigmanegpi, float);
  DECLARE_SOA_COLUMN(NTPCSigmaPosPi, ntpcsigmapospi, float);
  DECLARE_SOA_COLUMN(NTPCSigmaBachPi, ntpcsigmabachpi, float);
  DECLARE_SOA_COLUMN(NTOFSigmaNegPr, ntofsigmanegpr, float);
  DECLARE_SOA_COLUMN(NTOFSigmaPosPr, ntofsigmapospr, float);
  DECLARE_SOA_COLUMN(NTOFSigmaNegPi, ntofsigmanegpi, float);
  DECLARE_SOA_COLUMN(NTOFSigmaPosPi, ntofsigmapospi, float);
  DECLARE_SOA_COLUMN(NTOFSigmaBachPi, ntofsigmabachpi, float);
  DECLARE_SOA_COLUMN(PosNTPCClusters, posntpcscls, float);
  DECLARE_SOA_COLUMN(NegNTPCClusters, negntpcscls, float);
  DECLARE_SOA_COLUMN(BachNTPCClusters, bachntpcscls, float);
  DECLARE_SOA_COLUMN(PosHasTOF, poshastof, float);
  DECLARE_SOA_COLUMN(NegHasTOF, neghastof, float);
  DECLARE_SOA_COLUMN(BachHasTOF, bachhastof, float);

  } // namespace myv0candidates

  DECLARE_SOA_TABLE(MyCascades, "AOD", "MYCASCADES", o2::soa::Index<>,
                    mycascades::CollisionId, mycascades::Sign, mycascades::Pt, mycascades::RapXi, mycascades::RapOmega,
                    mycascades::MassXi, mycascades::MassOmega, mycascades::MassLambdaDau, mycascades::CascRadius, mycascades::V0Radius,
                    mycascades::CascCosPA, mycascades::V0CosPA, mycascades::DCAPosToPV, mycascades::DCANegToPV,
                    mycascades::DCABachToPV, mycascades::DCACascDaughters, mycascades::DCAV0Daughters, mycascades::DCAV0ToPV, mycascades::PosEta, mycascades::NegEta,
                    mycascades::BachEta, mycascades::PosITSHits, mycascades::NegITSHits, mycascades::BachITSHits,
                    mycascades::CtauXiMinus, mycascades::CtauXiPlus, mycascades::CtauOmegaMinus, mycascades::CtauOmegaPlus,
                    mycascades::NTPCSigmaNegPr, mycascades::NTPCSigmaPosPr, mycascades::NTPCSigmaNegPi, mycascades::NTPCSigmaPosPi,
                    mycascades::NTPCSigmaBachPi, mycascades::NTOFSigmaNegPr, mycascades::NTOFSigmaPosPr, mycascades::NTOFSigmaNegPi,
                    mycascades::NTOFSigmaPosPi, mycascades::NTOFSigmaBachPi,
                    mycascades::PosNTPCClusters, mycascades::NegNTPCClusters, mycascades::BachNTPCClusters,
                    mycascades::PosHasTOF, mycascades::NegHasTOF, mycascades::BachHasTOF);

} // namespace o2::aod

struct cascpostprocessing {

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};

  // Selection criteria
  Configurable<float> minpt{"minpt", 0.0, "Min p_{T} (GeV/c)"};
  Configurable<float> rap{"rap", 0.5, "Rapidity"};
  Configurable<float> masswin{"masswin", 0.075, "Mass window limit"};
  Configurable<float> rejcomp{"rejcomp", 0.008, "Competing Cascade rejection"};
  Configurable<float> lambdamasswin{"lambdamasswin", 0.008, "V0 Mass window limit"};
  Configurable<float> v0radius{"v0radius", 1.2, "V0 Radius"};
  Configurable<float> cascradius{"cascradius", 0.6, "Casc Radius"};
  Configurable<double> casccospa{"casccospa", 0.97, "Casc CosPA"};
  Configurable<double> v0cospa{"v0cospa", 0.97, "V0 CosPA"};
  Configurable<float> dcanegtopv{"dcanegtopv", 0.03, "DCA Neg To PV"};
  Configurable<float> dcapostopv{"dcapostopv", 0.03, "DCA Pos To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.03, "DCA Bach To PV"};
  Configurable<float> dcacascdau{"dcacascdau", 1.3, "DCA Casc Daughters"};
  Configurable<float> dcav0dau{"dcav0dau", 1.5, "DCA V0 Daughters"};
  Configurable<float> dcav0topv{"dcav0topv", 0.06, "DCA V0 To PV"};
  Configurable<float> etadau{"etadau", 0.8, "Eta Daughters"};
  Configurable<float> proplifetimexi{"proplifetimexi", 14.73, "Xi Propagation Lifetime"};
  Configurable<float> proplifetimeomega{"proplifetimeomega", 7.38, "Omega Propagation Lifetime"};
  Configurable<float> nsigmatpc{"nsigmatpc", 6, "N sigma TPC"};
  Configurable<float> nsigmatof{"nsigmatof", 6, "N sigma TOF"};
  Configurable<float> hasTOF{"hasTOF", 0, "Has TOF"};

  HistogramRegistry registry{"registryts"};

  void init(InitContext const&)
  {
    registry.add("hPt", "hPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}});
    registry.add("hXiMinusInvMass", "hXiMinusInvMass", {HistType::kTH1F, {{1000, 1.0f, 2.0f}}});
    registry.add("hOmegaMinusInvMass", "hOmegaMinusInvMass", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hXiPlusInvMass", "hXiPlusInvMass", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hOmegaPlusInvMass", "hOmegaPlusInvMass", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hXiMinusInvMass_BefSels", "hXiMinusInvMass_BefSels", {HistType::kTH1F, {{1000, 1.0f, 2.0f}}});
    registry.add("hOmegaMinusInvMass_BefSels", "hOmegaMinusInvMass_BefSels", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hXiPlusInvMass_BefSels", "hXiPlusInvMass_BefSels", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hOmegaPlusInvMass_BefSels", "hOmegaPlusInvMass_BefSels", {HistType::kTH1F, {{3000, 1.0f, 2.0f}}});
    registry.add("hCascRadius", "hCascRadius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
    registry.add("hV0Radius", "hV0Radius", {HistType::kTH1D, {{100, 0.0f, 40.0f}}});
    registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hCtauXiMinus", "hCtauXiMinus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauXiPlus", "hCtauXiPlus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauOmegaMinus", "hCtauOmegaMinus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauOmegaPlus", "hCtauOmegaPlus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
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
  }

  void process(aod::MyCascades const& mycascades)
  {
    for (auto& candidate : mycascades) {
      
      // To have trace of how it was before selections
      if(candidate.sign()<0){
        registry.fill(HIST("hXiMinusInvMass_BefSels"), candidate.massxi());
        registry.fill(HIST("hOmegaMinusInvMass_BefSels"), candidate.massomega());
      }
      if(candidate.sign()>0){
        registry.fill(HIST("hXiPlusInvMass_BefSels"), candidate.massxi());
        registry.fill(HIST("hOmegaPlusInvMass_BefSels"), candidate.massomega());
      }
      
      // Apply selections
      if (candidate.cascradius() > cascradius && candidate.v0radius() > v0radius &&
          candidate.casccospa() > casccospa && candidate.v0cospa() > v0cospa &&
          TMath::Abs(candidate.poseta()) < etadau && TMath::Abs(candidate.negeta()) < etadau && TMath::Abs(candidate.bacheta()) < etadau &&
          TMath::Abs(candidate.dcanegtopv()) > dcanegtopv && TMath::Abs(candidate.dcapostopv()) > dcapostopv && TMath::Abs(candidate.dcabachtopv()) > dcabachtopv &&
          candidate.dcacascdaughters() < dcacascdau && candidate.dcav0daughters() < dcav0dau &&
          TMath::Abs(candidate.dcav0topv()) > dcav0topv &&
          TMath::Abs(candidate.masslambdadau() - RecoDecay::getMassPDG(3122)) < lambdamasswin) {


        // massa inv vs pT vs raggio
        // se MC mass true reco vs pT vs raggio
        // distrubuzione di tutte le variabili
        // e tutti anche prima

         

        /*registry.fill(HIST("hPt"), candidate.pt());
        registry.fill(HIST("hCascRadius"), candidate.cascradius());
        registry.fill(HIST("hV0Radius"), candidate.v0radius());
        registry.fill(HIST("hCascCosPA"), candidate.casccospa());
        registry.fill(HIST("hV0CosPA"), candidate.v0cospa());
        registry.fill(HIST("hDCANegToPV"), candidate.dcanegtopv());
        registry.fill(HIST("hDCAPosToPV"), candidate.dcapostopv());
        registry.fill(HIST("hDCABachToPV"), candidate.dcabachtopv());
        registry.fill(HIST("hDCACascDaughters"), candidate.dcacascdaughters());
        registry.fill(HIST("hDCAV0Daughters"), candidate.dcav0daughters());
        registry.fill(HIST("hCtauXiMinus"), candidate.ctauximinus());
        registry.fill(HIST("hCtauXiPlus"), candidate.ctauxiplus());
        registry.fill(HIST("hCtauOmegaMinus"), candidate.ctauomegaminus());
        registry.fill(HIST("hCtauOmegaPlus"), candidate.ctauomegaplus());
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
        registry.fill(HIST("hBachITSHits"), candidate.bachitshits());*/
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascpostprocessing>(cfgc, TaskName{"lf-cascpostprocessing"})
    };
}
