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

/// \header file for histograms
/// \author daiki.sekihata@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_

#include <vector>
#include "TF1.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

using namespace o2::framework;
using namespace o2::aod::pwgem::photonmeson::utils::mcutil;

namespace o2::aod::pwgem::photonmeson::utils::nmhistogram
{
void addNMHistograms(HistogramRegistry* fRegistry, bool isMC, const char* pairname = "#gamma#gamma")
{
  // !!Don't change pt,eta,y binning. These binnings have to be consistent with binned data at skimming.!!
  std::vector<double> ptbins;
  for (int i = 0; i < 2; i++) {
    ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.05 GeV/c, every 0.05 GeV/c
  }
  for (int i = 2; i < 51; i++) {
    ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 4.9 GeV/c, every 0.1 GeV/c
  }
  for (int i = 51; i < 61; i++) {
    ptbins.emplace_back(0.5 * (i - 51) + 5.0); // from 5 to 9.5 GeV/c, every 0.5 GeV/c
  }
  for (int i = 61; i < 72; i++) {
    ptbins.emplace_back(1.0 * (i - 61) + 10.0); // from 10 to 20 GeV/c, every 1 GeV/c
  }
  const AxisSpec axis_pt{ptbins, Form("p_{T,%s} (GeV/c)", pairname)};
  const AxisSpec axis_mass{400, 0, 0.8, Form("m_{%s} (GeV/c^{2})", pairname)};

  if (isMC) {
    fRegistry->add("Pair/Pi0/hs_Primary", "rec. true pi0", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Pi0/hs_FromWD", "rec. true pi0 from weak decay", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Pi0/hs_FromHS", "rec. true pi0 from hadronic shower in material", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Pi0/hs_FromSameGamma", "Two clusters from same gamma that is a pi0 daughter (conversion)", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Eta/hs_FromSameGamma", "Two clusters from same gamma that is a eta daughter (conversion)", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Eta/hs_Primary", "rec. true eta", kTHnSparseD, {axis_mass, axis_pt}, true);

    const AxisSpec axis_rapidity{{0.0, +0.8, +0.9}, "rapidity |y|"};
    fRegistry->add("Generated/Pi0/hPt", "pT;p_{T} (GeV/c)", kTH1F, {axis_pt}, true);
    fRegistry->add("Generated/Pi0/hPtY", "Generated info", kTH2F, {axis_pt, axis_rapidity}, true);
    fRegistry->addClone("Generated/Pi0/", "Generated/Eta/");

    fRegistry->get<TH1>(HIST("Generated/Pi0/hPt"))->SetXTitle("p_{T} (GeV/c)");
    fRegistry->get<TH1>(HIST("Generated/Eta/hPt"))->SetXTitle("p_{T} (GeV/c)");
    fRegistry->get<TH2>(HIST("Generated/Pi0/hPtY"))->SetXTitle("p_{T} (GeV/c)");
    fRegistry->get<TH2>(HIST("Generated/Pi0/hPtY"))->SetYTitle("rapidity |y|");
    fRegistry->get<TH2>(HIST("Generated/Eta/hPtY"))->SetXTitle("p_{T} (GeV/c)");
    fRegistry->get<TH2>(HIST("Generated/Eta/hPtY"))->SetYTitle("rapidity |y|");
  } else {
    fRegistry->add("Pair/same/hs", "diphoton", kTHnSparseD, {axis_mass, axis_pt}, true);
    fRegistry->addClone("Pair/same/", "Pair/mix/");
  }
}

template <typename TDiphoton, typename TMCParitlce, typename TMCParticles, typename TMCCollisions>
void fillTruePairInfo(HistogramRegistry* fRegistry, TDiphoton const& v12, TMCParitlce const& mcparticle, TMCParticles const& mcparticles, TMCCollisions const&, const TF1* f1fd_k0s_to_pi0 = nullptr, float eventWeight = 1.f)
{
  int pdg = abs(mcparticle.pdgCode());
  float weight = eventWeight;
  switch (pdg) {
    case 111: {
      int motherid_strhad = IsFromWD(mcparticle.template emmcevent_as<TMCCollisions>(), mcparticle, mcparticles);
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Pi0/hs_Primary"), v12.M(), v12.Pt(), weight);
      } else if (motherid_strhad > 0) {
        auto str_had = mcparticles.iteratorAt(motherid_strhad);
        if (abs(str_had.pdgCode()) == 310 && f1fd_k0s_to_pi0 != nullptr) {
          weight *= f1fd_k0s_to_pi0->Eval(str_had.pt());
        }
        fRegistry->fill(HIST("Pair/Pi0/hs_FromWD"), v12.M(), v12.Pt(), weight);
      } else {
        fRegistry->fill(HIST("Pair/Pi0/hs_FromHS"), v12.M(), v12.Pt(), weight);
      }
      break;
    }
    case 221: {
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Eta/hs_Primary"), v12.M(), v12.Pt(), weight);
      }
      break;
    }
    default:
      break;
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::nmhistogram

#endif // PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
