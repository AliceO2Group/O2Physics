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
/// \brief post processing
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \modified by Chiara De Martin (chiara.de.martin@cern.ch)
/// \since March 20, 2023

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

// constants
const float ctauxi = 4.91;     // from PDG
const float ctauomega = 2.461; // from PDG

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using DauTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTOFPi, aod::pidTOFPr>;

namespace o2::aod
{

namespace mycascades
{

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

} // namespace mycascades

DECLARE_SOA_TABLE(MyCascades, "AOD", "MYCASCADES", o2::soa::Index<>,
                  mycascades::CollisionId, mycascades::Sign, mycascades::Pt, mycascades::RapXi, mycascades::RapOmega,
                  mycascades::MassXi, mycascades::MassOmega, mycascades::MassLambdaDau, mycascades::CascRadius, mycascades::V0Radius,
                  mycascades::CascCosPA, mycascades::V0CosPA, mycascades::DCAPosToPV, mycascades::DCANegToPV,
                  mycascades::DCABachToPV, mycascades::DCACascDaughters, mycascades::DCAV0Daughters, mycascades::DCAV0ToPV, mycascades::PosEta, mycascades::NegEta,
                  mycascades::BachEta, mycascades::PosITSHits, mycascades::NegITSHits, mycascades::BachITSHits,
                  mycascades::CtauXiMinus, mycascades::CtauXiPlus, mycascades::CtauOmegaMinus, mycascades::CtauOmegaPlus,
                  mycascades::NTPCSigmaNegPr, mycascades::NTPCSigmaPosPr, mycascades::NTPCSigmaNegPi, mycascades::NTPCSigmaPosPi,
                  mycascades::NTPCSigmaBachPi, mycascades::NTPCSigmaBachKa, mycascades::NTOFSigmaNegPr, mycascades::NTOFSigmaPosPr, mycascades::NTOFSigmaNegPi,
                  mycascades::NTOFSigmaPosPi, mycascades::NTOFSigmaBachPi, mycascades::NTOFSigmaBachKa,
                  mycascades::PosNTPCClusters, mycascades::NegNTPCClusters, mycascades::BachNTPCClusters,
                  mycascades::PosHasTOF, mycascades::NegHasTOF, mycascades::BachHasTOF);

} // namespace o2::aod

struct cascpostprocessing {
  // Xi or Omega
  Configurable<bool> isXi{"isXi", 1, "Apply cuts for Xi identification"};

  // Event selection criteria
  Configurable<float> cutzvertex{"cutzvertex", 15.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> sel8{"sel8", 1, "Apply sel8 event selection"};

  // Selection criteria
  Configurable<float> minpt{"minpt", 0.0, "Min p_{T} (GeV/c)"};
  Configurable<float> rap{"rap", 0.5, "Rapidity"};
  Configurable<float> masswin{"masswin", 0.075, "Mass window limit"};
  Configurable<float> masswintpc{"masswintpc", 0.075, "Mass window limit for Nsigma TPC daughters"};
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
  Configurable<float> proplifetime{"proplifetime", 6, "ctau/<ctau>"};
  Configurable<float> nsigmatpcPi{"nsigmatpcPi", 6, "N sigma TPC Pion"};
  Configurable<float> nsigmatpcPr{"nsigmatpcPr", 6, "N sigma TPC Proton"};
  Configurable<float> nsigmatpcKa{"nsigmatpcKa", 6, "N sigma TPC Kaon"};
  Configurable<float> nsigmatofPi{"nsigmatofPi", 6, "N sigma TOF Pion"};
  Configurable<float> nsigmatofPr{"nsigmatofPr", 6, "N sigma TOF Proton"};
  Configurable<float> nsigmatofKa{"nsigmatofKa", 6, "N sigma TOF Kaon"};
  Configurable<float> hasTOF{"hasTOF", 0, "Has TOF"};

  HistogramRegistry registry{"registryts"};

  void init(InitContext const&)
  {

    AxisSpec ximassAxis = {200, 1.28f, 1.36f};
    AxisSpec omegamassAxis = {200, 1.59f, 1.75f};
    AxisSpec massAxis = ximassAxis;
    if (!isXi)
      massAxis = omegamassAxis;
    AxisSpec ptAxis = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPCAxis = {100, 0.0f, 10.0f, "#it{p} TPC (GeV/#it{c})"};
    AxisSpec etaAxis = {200, -2.0f, 2.0f, "#eta"};
    AxisSpec rapidityAxis = {200, -2.0f, 2.0f, "y"};
    AxisSpec phiAxis = {100, -TMath::Pi() / 2, 3. * TMath::Pi() / 2, "#varphi"};

    registry.add("hPt", "hPt", {HistType::kTH1F, {ptAxis}});
    registry.add("hCascMinusInvMassvsPt", "hCascMinusInvMassvsPt", HistType::kTH2F, {ptAxis, massAxis});
    registry.add("hCascPlusInvMassvsPt", "hCascPlusInvMassvsPt", HistType::kTH2F, {ptAxis, massAxis});
    registry.add("hXiMinusInvMassvsPt_BefSels", "hXiMinusInvMassvsPt_BefSels", HistType::kTH2F, {ptAxis, ximassAxis});
    registry.add("hOmegaMinusInvMassvsPt_BefSels", "hOmegaMinusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, omegamassAxis}});
    registry.add("hXiPlusInvMassvsPt_BefSels", "hXiPlusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, ximassAxis}});
    registry.add("hOmegaPlusInvMassvsPt_BefSels", "hOmegaPlusInvMassvsPt_BefSels", {HistType::kTH2F, {ptAxis, omegamassAxis}});

    // topo
    registry.add("hDCANegToPV", "hDCANegToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hDCAPosToPV", "hDCAPosToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hDCABachToPV", "hDCABachToPV", {HistType::kTH1F, {{200, -1.0f, 1.0f}}});
    registry.add("hCascCosPA", "hCascCosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hV0CosPA", "hV0CosPA", {HistType::kTH1F, {{100, 0.9f, 1.0f}}});
    registry.add("hCascRadius", "hCascRadius", {HistType::kTH1D, {{500, 0.0f, 50.0f}}});
    registry.add("hV0Radius", "hV0Radius", {HistType::kTH1D, {{500, 0.0f, 50.0f}}});
    registry.add("hDCACascDaughters", "hDCACascDaughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hDCAV0Daughters", "hDCAV0Daughters", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hDCAV0ToPV", "hDCAV0ToPV", {HistType::kTH1F, {{55, 0.0f, 2.20f}}});
    registry.add("hMassLambdaDau", "hMassLambdaDau", {HistType::kTH1F, {{60, 1.1f, 1.13f}}});

    // kine
    registry.add("hEtaMinus", "hEtaMinus", {HistType::kTH2F, {ptAxis, etaAxis}});
    registry.add("hEtaPlus", "hEtaPlus", {HistType::kTH2F, {ptAxis, etaAxis}});
    registry.add("hRapMinus", "hRapMinus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("hRapMinus1D", "hRapMinus1D", {HistType::kTH1F, {rapidityAxis}});
    registry.add("hRapPlus", "hRapPlus", {HistType::kTH2F, {ptAxis, rapidityAxis}});
    registry.add("hRapPlus1D", "hRapPlus1D", {HistType::kTH1F, {rapidityAxis}});
    registry.add("hPhiMinus", "hPhiMinus", {HistType::kTH2F, {ptAxis, phiAxis}});
    registry.add("hPhiPlus", "hPhiPlus", {HistType::kTH2F, {ptAxis, phiAxis}});
    registry.add("hCtauMinus", "hCtauMinus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});
    registry.add("hCtauPlus", "hCtauPlus", {HistType::kTH1F, {{100, 0.0f, 40.0f}}});

    // daughter tracks
    registry.add("hTPCNSigmaPosPi", "hTPCNSigmaPosPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaNegPi", "hTPCNSigmaNegPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaPosPr", "hTPCNSigmaPosPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaNegPr", "hTPCNSigmaNegPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaBachPi", "hTPCNSigmaBachPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTPCNSigmaBachKa", "hTPCNSigmaBachKa", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaPosPi", "hTOFNSigmaPosPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaNegPi", "hTOFNSigmaNegPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaPosPr", "hTOFNSigmaPosPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaNegPr", "hTOFNSigmaNegPr", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaBachPi", "hTOFNSigmaBachPi", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hTOFNSigmaBachKa", "hTOFNSigmaBachKa", {HistType::kTH1F, {{120, -6.0f, 6.0f}}});
    registry.add("hCascMinusEtaPos", "hCascMinusEtaPos", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hCascMinusEtaNeg", "hCascMinusEtaNeg", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
    registry.add("hCascMinusEtaBach", "hCascMinusEtaBach", {HistType::kTH1F, {{100, -1.0f, 1.0f}}});
  }

  void process(aod::MyCascades const& mycascades)
  {
    float invmass = 0;
    float ctau = 0;
    float rapidity = 0;
    bool isCandidate = 0;

    for (auto& candidate : mycascades) {

      // To have trace of how it was before selections
      if (candidate.sign() < 0) {
        registry.fill(HIST("hXiMinusInvMassvsPt_BefSels"), candidate.pt(), candidate.massxi());
        registry.fill(HIST("hOmegaMinusInvMassvsPt_BefSels"), candidate.pt(), candidate.massomega());
      }
      if (candidate.sign() > 0) {
        registry.fill(HIST("hXiPlusInvMassvsPt_BefSels"), candidate.pt(), candidate.massxi());
        registry.fill(HIST("hOmegaPlusInvMassvsPt_BefSels"), candidate.pt(), candidate.massomega());
      }

      // Apply selections
      if (TMath::Abs(candidate.poseta()) > etadau || TMath::Abs(candidate.negeta()) > etadau || TMath::Abs(candidate.bacheta()) > etadau)
        continue;
      if (TMath::Abs(candidate.dcanegtopv()) < dcanegtopv || TMath::Abs(candidate.dcapostopv()) < dcapostopv || TMath::Abs(candidate.dcabachtopv()) < dcabachtopv)
        continue;
      if (candidate.casccospa() < casccospa)
        continue;
      if (candidate.v0cospa() < v0cospa)
        continue;
      if (candidate.dcacascdaughters() > dcacascdau)
        continue;
      if (candidate.dcav0daughters() > dcav0dau)
        continue;
      if (candidate.cascradius() < cascradius)
        continue;
      if (candidate.v0radius() < v0radius)
        continue;
      if (TMath::Abs(candidate.dcav0topv()) < dcav0topv)
        continue;
      if (TMath::Abs(candidate.masslambdadau() - RecoDecay::getMassPDG(3122)) > lambdamasswin)
        continue;

      if (candidate.sign() < 0) {
        if (TMath::Abs(candidate.ntpcsigmapospr()) > nsigmatpcPr)
          continue;
        if (TMath::Abs(candidate.ntpcsigmanegpi()) > nsigmatpcPi)
          continue;
      } else if (candidate.sign() > 0) {
        if (TMath::Abs(candidate.ntpcsigmanegpr()) > nsigmatpcPr)
          continue;
        if (TMath::Abs(candidate.ntpcsigmapospi()) > nsigmatpcPi)
          continue;
      }
      // TOF PID only for daughters with pT > threshold (to be implemented)
      /*
  if (candidate.sign() < 0){
  if (TMath::Abs(candidate.ntofsigmapospr()) > nsigmatofPr) continue;
  if (TMath::Abs(candidate.ntofsigmanegpi()) > nsigmatofPi) continue;
  }
  else if (candidate.sign() > 0){
  if (TMath::Abs(candidate.ntofsigmanegpr()) > nsigmatofPr) continue;
  if (TMath::Abs(candidate.ntofsigmapospi()) > nsigmatofPi) continue;
  }
      */

      if (isXi) {
        if (TMath::Abs(candidate.ntpcsigmabachpi()) > nsigmatpcPi)
          continue;
        if (TMath::Abs(candidate.rapxi()) > rap)
          continue;
        if (candidate.ctauximinus() > proplifetime * ctauxi)
          continue;
        if (TMath::Abs(candidate.massxi() - RecoDecay::getMassPDG(3312)) > masswin)
          continue;
        if (TMath::Abs(candidate.massomega() - RecoDecay::getMassPDG(3334)) < rejcomp)
          continue;
        rapidity = candidate.rapxi();
        ctau = candidate.ctauximinus();
        invmass = candidate.massxi();
      } else {
        if (TMath::Abs(candidate.ntpcsigmabachka()) > nsigmatpcKa)
          continue;
        if (TMath::Abs(candidate.rapxi()) > rap)
          continue;
        if (candidate.ctauomegaminus() > proplifetime * ctauomega)
          continue;
        if (TMath::Abs(candidate.massomega() - RecoDecay::getMassPDG(3334)) > masswin)
          continue;
        if (TMath::Abs(candidate.massxi() - RecoDecay::getMassPDG(3312)) < rejcomp)
          continue;
        rapidity = candidate.rapomega();
        ctau = candidate.ctauomegaminus();
        invmass = candidate.massomega();
      }

      registry.fill(HIST("hPt"), candidate.pt());
      registry.fill(HIST("hDCANegToPV"), candidate.dcanegtopv());
      registry.fill(HIST("hDCAPosToPV"), candidate.dcapostopv());
      registry.fill(HIST("hDCABachToPV"), candidate.dcabachtopv());
      registry.fill(HIST("hCascCosPA"), candidate.casccospa());
      registry.fill(HIST("hV0CosPA"), candidate.v0cospa());
      registry.fill(HIST("hCascRadius"), candidate.cascradius());
      registry.fill(HIST("hV0Radius"), candidate.v0radius());
      registry.fill(HIST("hDCACascDaughters"), candidate.dcacascdaughters());
      registry.fill(HIST("hDCAV0Daughters"), candidate.dcav0daughters());
      registry.fill(HIST("hDCAV0ToPV"), candidate.dcav0topv());
      registry.fill(HIST("hMassLambdaDau"), candidate.masslambdadau());

      if (candidate.sign() > 0) {
        registry.fill(HIST("hCtauPlus"), ctau);
        // registry.fill(HIST("hEtaPlus"), candidate.pt(), candidate.eta());
        registry.fill(HIST("hRapPlus"), candidate.pt(), rapidity);
        registry.fill(HIST("hRapPlus1D"), rapidity);
        // registry.fill(HIST("hPhiPlus"), candidate.pt(), candidate.phi());
      } else {
        registry.fill(HIST("hCtauMinus"), ctau);
        // registry.fill(HIST("hEtaMinus"), candidate.pt(), candidate.eta());
        registry.fill(HIST("hRapMinus"), candidate.pt(), rapidity);
        registry.fill(HIST("hRapMinus1D"), rapidity);
        // registry.fill(HIST("hPhiMinus"), candidate.pt(), candidate.phi());
      }

      if (isXi) {
        if (TMath::Abs(candidate.massxi() - RecoDecay::getMassPDG(3312)) < masswintpc)
          isCandidate = 1;
      } else if (!isXi) {
        if (TMath::Abs(candidate.massomega() - RecoDecay::getMassPDG(3334)) < masswintpc)
          isCandidate = 1;
      }
      if (isCandidate) {
        if (candidate.sign() < 0) {
          registry.fill(HIST("hTPCNSigmaPosPr"), candidate.ntpcsigmapospr());
          registry.fill(HIST("hTPCNSigmaNegPi"), candidate.ntpcsigmanegpi());
          registry.fill(HIST("hTOFNSigmaPosPr"), candidate.ntofsigmapospr());
          registry.fill(HIST("hTOFNSigmaNegPi"), candidate.ntofsigmanegpi());
          registry.fill(HIST("hCascMinusEtaPos"), candidate.poseta());
          registry.fill(HIST("hCascMinusEtaNeg"), candidate.negeta());
          registry.fill(HIST("hCascMinusEtaBach"), candidate.bacheta());
        } else {
          registry.fill(HIST("hTPCNSigmaPosPi"), candidate.ntpcsigmapospi());
          registry.fill(HIST("hTPCNSigmaNegPr"), candidate.ntpcsigmanegpr());
          registry.fill(HIST("hTOFNSigmaPosPi"), candidate.ntofsigmapospi());
          registry.fill(HIST("hTOFNSigmaNegPr"), candidate.ntofsigmanegpr());
        }
        if (isXi) {
          registry.fill(HIST("hTPCNSigmaBachPi"), candidate.ntpcsigmabachpi());
          registry.fill(HIST("hTOFNSigmaBachPi"), candidate.ntofsigmabachpi());
        } else {
          registry.fill(HIST("hTPCNSigmaBachKa"), candidate.ntpcsigmabachka());
          registry.fill(HIST("hTOFNSigmaBachKa"), candidate.ntofsigmabachka());
        }
      }
      // registry.fill(HIST("hPosITSHits"), candidate.positshits());
      // registry.fill(HIST("hNegITSHits"), candidate.negitshits());
      // registry.fill(HIST("hBachITSHits"), candidate.bachitshits());

      if (candidate.sign() < 0) {
        registry.fill(HIST("hCascMinusInvMassvsPt"), candidate.pt(), invmass);
      }
      if (candidate.sign() > 0) {
        registry.fill(HIST("hCascPlusInvMassvsPt"), candidate.pt(), invmass);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<cascpostprocessing>(cfgc, TaskName{"lf-cascpostprocessing"})};
}
