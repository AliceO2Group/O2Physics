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
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"

#ifndef PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
#define PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_

using namespace o2::framework;
using namespace o2::aod::photonpair;
using namespace o2::aod::pwgem::mcutil;

namespace o2::aod::pwgem::photonmeson::utils::nmhistogram
{
void addNMHistograms(HistogramRegistry* fRegistry, bool doFlow, bool isMC, const char* pairname = "#gamma#gamma")
{
  std::vector<double> ptbins;
  for (int i = 0; i < 50; i++) {
    ptbins.emplace_back(0.1 * (i - 0) + 0.0); // from 0 to 5 GeV/c, every 0.1 GeV/c
  }
  for (int i = 50; i < 60; i++) {
    ptbins.emplace_back(0.5 * (i - 50) + 5.0); // from 5 to 10 GeV/c, evety 0.5 GeV/c
  }
  for (int i = 60; i < 71; i++) {
    ptbins.emplace_back(1.0 * (i - 60) + 10.0); // from 10 to 20 GeV/c, evety 1 GeV/c
  }
  const AxisSpec axis_pt{ptbins, Form("p_{T,%s} (GeV/c)", pairname)};

  const AxisSpec axis_mass{400, 0, 0.8, Form("m_{%s} (GeV/c^{2})", pairname)};
  const AxisSpec sp2_mass{100, -5, 5, Form("u_{2}^{%s} #upoint Q_{2}", pairname)};

  if (isMC) {
    fRegistry->add("Generated/Pi0/hPt", "pT;p_{T} (GeV/c)", kTH1F, {{2000, 0.0f, 20}});
    fRegistry->add("Generated/Pi0/hY", "y;rapidity y", kTH1F, {{40, -2.0f, 2.0f}});
    fRegistry->add("Generated/Pi0/hPhi", "#varphi;#varphi (rad.)", kTH1F, {{180, 0, 2 * M_PI}});
    fRegistry->add("Generated/Pi0/hPt_Acc", "pT;p_{T} (GeV/c)", kTH1F, {{2000, 0.0f, 20}});          // in pair acceptance
    fRegistry->add("Generated/Pi0/hY_Acc", "y;rapidity y", kTH1F, {{40, -2.0f, 2.0f}});              // in pair acceptance
    fRegistry->add("Generated/Pi0/hPhi_Acc", "#varphi;#varphi (rad.)", kTH1F, {{180, 0, 2 * M_PI}}); // in pair acceptance
    fRegistry->addClone("Generated/Pi0/", "Generated/Eta/");

    fRegistry->add("Pair/Pi0/hMggPt_Primary", "rec. true pi0", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Pi0/hMggPt_FromWD", "rec. true pi0 from weak decay", kTH2F, {axis_mass, axis_pt}, true);
    fRegistry->add("Pair/Eta/hMggPt_Primary", "rec. true eta", kTH2F, {axis_mass, axis_pt}, true);
  } else {
    if (doFlow) {
      std::string_view sp_names[4] = {"FT0M", "FT0A", "FT0C", "FV0A"};
      for (int i = 0; i < 4; i++) {
        fRegistry->add(Form("Pair/same/hs_same_SPQ2%s", sp_names[i].data()), "2photon", kTHnSparseF, {axis_mass, axis_pt, sp2_mass}, true);
      }
      fRegistry->add("Pair/mix/hMggPt", "2photon", kTH2F, {axis_mass, axis_pt}, true);
    } else {
      fRegistry->add("Pair/same/hMggPt", "2photon", kTH2F, {axis_mass, axis_pt}, true);
      fRegistry->addClone("Pair/same/", "Pair/mix/");
    }
  }
}

template <int ev_id, PairType pairtype, typename TCollision, typename TDiphoton>
void fillPairInfo(HistogramRegistry* fRegistry, TCollision const& collision, TDiphoton const& diphoton, const bool doFlow, const float weight = 1.f)
{
  static constexpr std::string_view event_pair_types[2] = {"same/", "mix/"};
  if (doFlow) {
    if constexpr (ev_id == 0) { // same
      std::array<float, 2> q2ft0m = {collision.q2xft0m(), collision.q2yft0m()};
      std::array<float, 2> q2ft0a = {collision.q2xft0a(), collision.q2yft0a()};
      std::array<float, 2> q2ft0c = {collision.q2xft0c(), collision.q2yft0c()};
      std::array<float, 2> q2fv0a = {collision.q2xfv0a(), collision.q2yfv0a()};
      std::array<float, 2> u2_gg = {static_cast<float>(std::cos(2 * diphoton.Phi())), static_cast<float>(std::sin(2 * diphoton.Phi()))};

      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0M"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0m));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0A"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0a));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FT0C"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2ft0c));
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hs_same_SPQ2FV0A"), diphoton.M(), diphoton.Pt(), RecoDecay::dotProd(u2_gg, q2fv0a));
    } else { // mix
      fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMggPt"), diphoton.M(), diphoton.Pt());
    }
  } else {
    fRegistry->fill(HIST("Pair/") + HIST(event_pair_types[ev_id]) + HIST("hMggPt"), diphoton.M(), diphoton.Pt());
  }
}

template <typename TDiphoton, typename TMCParitlce, typename TMCParticles>
void fillTruePairInfo(HistogramRegistry* fRegistry, TDiphoton const& v12, TMCParitlce const& mcparticle, TMCParticles const& mcparticles, const float weight = 1.f)
{
  int pdg = abs(mcparticle.pdgCode());
  switch (pdg) {
    case 111: {
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Pi0/hMggPt_Primary"), v12.M(), v12.Pt());
      } else if (IsFromWD(mcparticle.emmcevent(), mcparticle, mcparticles)) {
        fRegistry->fill(HIST("Pair/Pi0/hMggPt_FromWD"), v12.M(), v12.Pt());
      }
      break;
    }
    case 221: {
      if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
        fRegistry->fill(HIST("Pair/Eta/hMggPt_Primary"), v12.M(), v12.Pt());
      }
      break;
    }
    default:
      break;
  }
}

} // namespace o2::aod::pwgem::photonmeson::utils::nmhistogram

#endif // PWGEM_PHOTONMESON_UTILS_NMHISTOGRAMS_H_
