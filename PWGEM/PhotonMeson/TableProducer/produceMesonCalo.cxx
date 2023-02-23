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

/// \brief perform calo photon analysis on calo photons from skimmergammacalo task
/// dependencies: skimmergammacalo
/// \author marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/DataModel/mesonTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

// includes for the R recalculation
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct produceMesonCalo {

  Produces<aod::CaloMeson> tableCaloMeson;

  HistogramRegistry spectra = {
    "spectra",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Configurable for histograms
  Configurable<int> nBinsMinv{"nBinsMinv", 800, "N bins for minv axis"};
  Configurable<float> minMinv{"minMinv", 0.0, "Minimum value for minv axis"};
  Configurable<float> maxMinv{"maxMinv", 0.8, "Maximum value for minv axis"};
  Configurable<int> nBinsPt{"nBinsPt", 180, "N bins for pT axis"};
  Configurable<float> minPt{"minPt", 0., "Minimum value for pT axis"};
  Configurable<float> maxPt{"maxPt", 60., "Maximum value for pT axis"};

  void init(o2::framework::InitContext&)
  {
    std::vector<double> ptBinning(nBinsPt, 0);

    for (int i = 0; i < nBinsPt; i++) {
      if (i < 100) {
        ptBinning.at(i) = 0.10 * i;
      } else if (i < 140) {
        ptBinning.at(i) = 10. + 0.25 * (i - 100);
      } else if (i < 180) {
        ptBinning.at(i) = 20. + 1.00 * (i - 140);
      } else {
        ptBinning.at(i) = maxPt;
      }
    }

    AxisSpec ptAxis = {ptBinning, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec minvAxis = {nBinsMinv, minMinv, maxMinv,
                         "#it{m}_{inv} (GeV/#it{c}^{2})"};
    AxisSpec etaAxis = {100, -0.8, 0.8, "#eta"};
    AxisSpec phiAxis = {360, 0, 2 * M_PI, "#varphi (rad)"};
    AxisSpec alphaAxis = {200, -1, +1, "#alpha"};
    AxisSpec oaAxis = {180, 0, M_PI, "#vartheta_{#gamma#gamma} (rad)"};

    HistogramConfigSpec defaultPtMinvHist(
      {HistType::kTH2F, {minvAxis, ptAxis}});

    HistogramConfigSpec defaultEtaPhiHist(
      {HistType::kTH2F, {etaAxis, phiAxis}});

    HistogramConfigSpec defaultPtMotherPtGammaHist(
      {HistType::kTH2F, {ptAxis, ptAxis}});

    HistogramConfigSpec defaultPtAlpha(
      {HistType::kTH2F, {ptAxis, alphaAxis}});

    HistogramConfigSpec defaultPtOA(
      {HistType::kTH2F, {ptAxis, oaAxis}});

    spectra.add("SameEvent_Minv_Pt", "SameEvent_Minv_Pt", defaultPtMinvHist, true);
    spectra.add("SameEvent_Eta_Phi", "SameEvent_Eta_Phi", defaultEtaPhiHist, true);
    spectra.add("SameEvent_Pt_Alpha", "SameEvent_Pt_Alpha", defaultPtAlpha, true);
    spectra.add("SameEvent_Pt_OA", "SameEvent_Pt_OA", defaultPtOA, true);
    spectra.add("SameEvent_PtMother_PtGamma", "SameEvent_PtMother_PtGamma", defaultPtMotherPtGammaHist, true);

    spectra.add("Photon_Eta_Phi", "Photon_Eta_Phi", defaultEtaPhiHist, true);
  }

  void
    processRec(aod::Collision const&,
               aod::SkimGammas const& skimgammas)
  {
    for (auto& [gamma0, gamma1] : // EMC-EMC
         combinations(o2::soa::CombinationsStrictlyUpperIndexPolicy(skimgammas,
                                                                    skimgammas))) {
      float openingAngle = acos((cos(gamma0.phi() - gamma1.phi()) +
                                 sinh(gamma0.eta()) * sinh(gamma1.eta())) /
                                (cosh(gamma0.eta()) * cosh(gamma1.eta())));
      float E = gamma0.e() + gamma1.e();
      float pt0 = gamma0.e() / cosh(gamma0.eta());
      float pt1 = gamma1.e() / cosh(gamma1.eta());
      float px =
        pt0 * cos(gamma0.phi()) + pt1 * cos(gamma1.phi());
      float py =
        pt0 * sin(gamma0.phi()) + pt1 * sin(gamma1.phi());
      float pz =
        pt0 * sinh(gamma0.eta()) + pt1 * sinh(gamma1.eta());
      float alpha = (gamma0.e() - gamma1.e()) != 0.
                      ? (gamma0.e() - gamma1.e()) / (gamma0.e() + gamma1.e())
                      : 0.;
      float Pt = sqrt(pt0 * pt0 + pt1 * pt1 +
                      2. * pt0 * pt1 *
                        cos(gamma0.phi() - gamma1.phi()));
      float minv =
        sqrt(2. * gamma0.e() * gamma1.e() * (1. - cos(openingAngle)));
      float eta = asinh(pz / Pt);
      float phi = atan2(py, px);
      phi = (phi < 0) ? phi + 2. * M_PI : phi;
      tableCaloMeson(gamma0.collisionId(), gamma0.globalIndex(), gamma1.globalIndex(),
                     openingAngle, px, py, pz, E, alpha, minv, eta, phi,
                     Pt);
      spectra.get<TH2>(HIST("SameEvent_Minv_Pt"))->Fill(minv, Pt);
      spectra.get<TH2>(HIST("SameEvent_Eta_Phi"))->Fill(eta, phi);
      spectra.get<TH2>(HIST("SameEvent_Pt_Alpha"))->Fill(Pt, alpha);
      spectra.get<TH2>(HIST("SameEvent_Pt_OA"))->Fill(Pt, openingAngle);
      spectra.get<TH2>(HIST("SameEvent_PtMother_PtGamma"))->Fill(Pt, pt0);
      spectra.get<TH2>(HIST("SameEvent_PtMother_PtGamma"))->Fill(Pt, pt1);

      spectra.get<TH2>(HIST("Photon_Eta_Phi"))->Fill(gamma0.eta(), gamma0.phi());
      spectra.get<TH2>(HIST("Photon_Eta_Phi"))->Fill(gamma1.eta(), gamma1.phi());
    }
  }
  PROCESS_SWITCH(produceMesonCalo, processRec,
                 "process only reconstructed info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<produceMesonCalo>(cfgc)};
  return workflow;
}
