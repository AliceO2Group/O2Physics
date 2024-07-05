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

#include <vector>
#include "TF1.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#ifndef PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_

using namespace o2::framework;
using namespace o2::aod::photonpair;
using namespace o2::aod::pwgem::mcutil;

namespace o2::aod::pwgem::photonmeson::utils::nmhistogram
{
void addNMHistograms(HistogramRegistry* fRegistry, bool do_v2, bool do_v3, bool isMC, const char* pairname = "#gamma#gamma", const char* epdetname = "")
{
  // !!Don't change pt,eta,y binning. These binnings have to be consistent with binned data at skimming.!!
  std::vector<double> ptbins;
  for (int i = 0; i < 2; i++) {
    ptbins.emplace_back(0.05 * (i - 0) + 0.0); // from 0 to 0.1 GeV/c, every 0.05 GeV/c
  }
  for (int i = 2; i < 52; i++) {
    ptbins.emplace_back(0.1 * (i - 2) + 0.1); // from 0.1 to 5 GeV/c, every 0.1 GeV/c
  }
  for (int i = 52; i < 62; i++) {
    ptbins.emplace_back(0.5 * (i - 52) + 5.0); // from 5 to 10 GeV/c, evety 0.5 GeV/c
  }
  for (int i = 62; i < 73; i++) {
    ptbins.emplace_back(1.0 * (i - 62) + 10.0); // from 10 to 20 GeV/c, evety 1 GeV/c
  }
  const AxisSpec axis_pt{ptbins, Form("p_{T,%s} (GeV/c)", pairname)};

  const AxisSpec axis_mass{400, 0, 0.8, Form("m_{%s} (GeV/c^{2})", pairname)};

  int nbin_sp2 = 1;
  int nbin_sp3 = 1;
  if (!isMC) {
    if (do_v2) {
      nbin_sp2 = 100;
    }
    if (do_v3) {
      nbin_sp3 = 100;
    }
  }
  const AxisSpec axis_sp2{nbin_sp2, -5.f, 5.f, Form("u_{2}^{%s} #upoint Q_{2}^{%s}", pairname, epdetname)};
  const AxisSpec axis_sp3{nbin_sp3, -5.f, 5.f, Form("u_{3}^{%s} #upoint Q_{3}^{%s}", pairname, epdetname)};

  if (isMC) {
    fRegistry->add("Pair/Pi0/hs_Primary", "rec. true pi0", kTHnSparseD, {axis_mass, axis_pt, axis_sp2, axis_sp3}, true);
    fRegistry->add("Pair/Pi0/hs_FromWD", "rec. true pi0 from weak decay", kTHnSparseD, {axis_mass, axis_pt, axis_sp2, axis_sp3}, true);
    fRegistry->add("Pair/Eta/hs_Primary", "rec. true eta", kTHnSparseD, {axis_mass, axis_pt, axis_sp2, axis_sp3}, true);

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
    fRegistry->add("Pair/same/hs", "diphoton", kTHnSparseD, {axis_mass, axis_pt, axis_sp2, axis_sp3}, true);
    fRegistry->addClone("Pair/same/", "Pair/mix/");
  }
}

template <typename TDiphoton, typename TMCParitlce, typename TMCParticles, typename TMCCollisions>
void fillTruePairInfo(HistogramRegistry* fRegistry, TDiphoton const& v12, TMCParitlce const& mcparticle, TMCParticles const& mcparticles, TMCCollisions const&, const TF1* f1fd_k0s_to_pi0 = nullptr)
{
  int pdg = abs(mcparticle.pdgCode());
  switch (pdg) {
    case 111: {
      int motherid_strhad = IsFromWD(mcparticle.template emmcevent_as<TMCCollisions>(), mcparticle, mcparticles);
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Pi0/hs_Primary"), v12.M(), v12.Pt(), 0.0, 0.0);
      } else if (motherid_strhad > 0) {
        float weight = 1.f;
        auto str_had = mcparticles.iteratorAt(motherid_strhad);
        if (abs(str_had.pdgCode()) == 310 && f1fd_k0s_to_pi0 != nullptr) {
          weight = f1fd_k0s_to_pi0->Eval(str_had.pt());
        }
        fRegistry->fill(HIST("Pair/Pi0/hs_FromWD"), v12.M(), v12.Pt(), 0.0, 0.0, weight);
      }
      break;
    }
    case 221: {
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Eta/hs_Primary"), v12.M(), v12.Pt(), 0.0, 0.0);
      }
      break;
    }
    default:
      break;
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::nmhistogram

#endif // PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
