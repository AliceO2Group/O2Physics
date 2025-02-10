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
// ========================
//
// This code is to study MC truth. e.g. S/B
//    Please write to: daiki.sekihata@cern.ch

#include <vector>
#include <unordered_map>
#include "Math/Vector4D.h"

#include "Framework/StaticFor.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TableHelper.h"
#include "PWGEM/Dilepton/Utils/MCUtilities.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;

struct studyMCTruth {
  Configurable<int> pdg_lepton{"pdg_lepton", 11, "pdg code for desired lepton"};
  Configurable<float> max_rxy_PC{"max_rxy_PC", 2.4, "max. rxy for electron from photon conversions"}; // beam pipe:1.8 cm, ITSLayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f}; in cm
  Configurable<float> min_imp_par{"min_imp_par", -1.f, "min. impact parameter in fm"};                // [0, 4] fm for centrality FT0C 0-10%,  [8, 10] fm for centrality FT0C 30-50%,
  Configurable<float> max_imp_par{"max_imp_par", 999.f, "max. impact parameter in fm"};

  Configurable<float> min_pt_gen{"min_pt_gen", 0.4, "min. pT of single lepton"};
  Configurable<float> max_pt_gen{"max_pt_gen", 1e+10f, "max. pT of single lepton"};
  Configurable<float> min_eta_gen{"min_eta_gen", -0.8, "min. eta of for single lepton"};
  Configurable<float> max_eta_gen{"max_eta_gen", +0.8, "max. eta of for single lepton"};

  Configurable<float> min_pt_gen_pf{"min_pt_gen_pf", 0.05, "min. pT of single lepton used for prefilter"};
  Configurable<float> max_pt_gen_pf{"max_pt_gen_pf", 1e+10f, "max. pT of single lepton used for prefilter"};
  Configurable<float> min_eta_gen_pf{"min_eta_gen_pf", -0.9, "min. eta of single lepton used for prefilter"};
  Configurable<float> max_eta_gen_pf{"max_eta_gen_pf", +0.9, "max. eta of single lepton used for prefilter"};
  Configurable<std::vector<float>> max_mll_vec{"max_mll_vec", std::vector<float>{0.08, 0.10, 0.12}, "vector fo max mll for prefilter in ULS. Please use exactly 3 values."}; // currently, 3 thoresholds are allowed.

  ConfigurableAxis ConfMllBins{"ConfMllBins", {400, 0.f, 4.f}, "mll bins for output histograms"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {100, 0.f, 10.f}, "pTll bins for output histograms"};

  HistogramRegistry fRegistry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  float leptonMass = o2::constants::physics::MassElectron;
  void init(o2::framework::InitContext&)
  {
    if (std::abs(pdg_lepton.value) == 11) {
      leptonMass = o2::constants::physics::MassElectron;
    } else if (std::abs(pdg_lepton.value) == 13) {
      leptonMass = o2::constants::physics::MassMuon;
    }
    addHistograms();
  }

  static constexpr std::string_view dileptonSigns[3] = {"uls/", "lspp/", "lsmm/"};
  static constexpr std::string_view pfNames[4] = {"default/", "mllPF0/", "mllPF1/", "mllPF2/"};

  void addHistograms()
  {
    const AxisSpec axis_mll{ConfMllBins, "m_{ll} (GeV/c^{2})"};
    const AxisSpec axis_ptll{ConfPtllBins, "p_{T,ll} (GeV/c)"};

    fRegistry.add("Event/hZvtx", "MC Zvtx;Z_{vtx} (cm)", kTH1D, {{100, -50, +50}}, false);
    fRegistry.add("Event/hImpactParameter", "impact parameter;impact parameter b (fm)", kTH1D, {{200, 0, 20}}, false);
    fRegistry.add("Event/hPtY_pion", "#pi^{#pm} yield;y;p_{T} (GeV/c)", kTH2D, {{20, -1, +1}, {100, 0, 10}}, false);
    fRegistry.add("Event/hPtY_kaon", "K^{#pm} yield;y;p_{T} (GeV/c)", kTH2D, {{20, -1, +1}, {100, 0, 10}}, false);
    fRegistry.add("Event/hPtY_proton", "p(#bar{p}) yield;y;p_{T} (GeV/c)", kTH2D, {{20, -1, +1}, {100, 0, 10}}, false);

    fRegistry.add("Track/default/hPt", "p_{T,e}", kTH1D, {{1000, 0, 10}}, true);
    fRegistry.add("Track/default/hEtaPhi", "#eta_{e} vs. #varphi_{e};#varphi_{e} (rad.);#eta_{e};", kTH2D, {{90, 0, 2 * M_PI}, {60, -5, +1}}, true);
    fRegistry.addClone("Track/default/", "Track/mllPF0/");
    fRegistry.addClone("Track/default/", "Track/mllPF1/");
    fRegistry.addClone("Track/default/", "Track/mllPF2/");

    fRegistry.add("Pair/default/Pi0/uls/hMvsPt", "m_{ll} vs. p_{T,ll}", kTH2D, {axis_mll, axis_ptll}, true);
    fRegistry.addClone("Pair/default/Pi0/uls/", "Pair/default/Pi0/lspp/");
    fRegistry.addClone("Pair/default/Pi0/uls/", "Pair/default/Pi0/lsmm/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Eta/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/EtaPrime/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Rho/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Omega/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Phi/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/JPsi/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Psi2S/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Upsilon1S/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Upsilon2S/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Upsilon3S/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/ccbar/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/bbbar/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/Photon/");
    fRegistry.addClone("Pair/default/Pi0/", "Pair/default/comb_bkg/");
    fRegistry.addClone("Pair/default/", "Pair/mllPF0/");
    fRegistry.addClone("Pair/default/", "Pair/mllPF1/");
    fRegistry.addClone("Pair/default/", "Pair/mllPF2/");
    fRegistry.addClone("Pair/default/Pi0/uls/", "Pair/PF/all/uls/");
  }

  template <typename TTrack, typename TMCParticles>
  int FindLF(TTrack const& posmc, TTrack const& negmc, TMCParticles const& mcparticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 22, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 111, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 221, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 331, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 113, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 223, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 333, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 100443, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 100553, mcparticles),
      FindCommonMotherFrom2Prongs(posmc, negmc, -pdg_lepton, pdg_lepton, 200553, mcparticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isEleFromPCOnITSib(TMCParticle const& mcparticle, TMCParticles const& mcParticles)
  {
    if (!mcparticle.has_mothers()) {
      return false;
    }
    if (std::abs(mcparticle.pdgCode()) != 11) {
      return false;
    }
    if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
      return false;
    }
    if (std::sqrt(std::pow(mcparticle.vx(), 2) + std::pow(mcparticle.vy(), 2)) > max_rxy_PC) {
      return false;
    }
    const auto& mp = mcParticles.iteratorAt(mcparticle.mothersIds()[0]);
    if (std::abs(mp.pdgCode()) == 22) {
      return true;
    } else {
      return false;
    }
  }

  template <typename TMCParticle, typename TMCParticles>
  bool isSelectedMCParticle(TMCParticle const& mcparticle, TMCParticles const& mcParticles)
  {
    if (std::abs(mcparticle.pdgCode()) != pdg_lepton) {
      return false;
    }
    if (!mcparticle.has_mothers()) {
      return false;
    }
    if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator() || isEleFromPCOnITSib(mcparticle, mcParticles)) {
      return true;
    } else {
      return false;
    }
  }

  template <int pftype, int signtype, typename TMCLepton, typename TMCParticles>
  void fillTrueInfo(TMCLepton const& t1, TMCLepton const& t2, TMCParticles const& mcParticles)
  {
    if (!isSelectedMCParticle(t1, mcParticles) || !isSelectedMCParticle(t2, mcParticles)) {
      return;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonMass);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonMass);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    int mother_id = FindLF(t1, t2, mcParticles);
    int hfll_type = IsHF(t1, t2, mcParticles);
    if (mother_id < 0 && hfll_type < 0) {
      fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("comb_bkg/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
      return;
    }

    if (mother_id > 0) { // same mother (photon, LF, Quarkonia)
      const auto& mp = mcParticles.iteratorAt(mother_id);
      if (std::abs(mp.pdgCode()) == 22 && isEleFromPCOnITSib(t1, mcParticles) && isEleFromPCOnITSib(t2, mcParticles)) {
        fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Photon/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        return;
      }
      if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
        return;
      }
      if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
        return;
      }
      if (!mp.isPhysicalPrimary() && !mp.producedByGenerator()) {
        return;
      }
      switch (std::abs(mp.pdgCode())) {
        case 111:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Pi0/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 221:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Eta/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 331:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("EtaPrime/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 113:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Rho/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 223:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Omega/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 333:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Phi/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 443:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("JPsi/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 100443:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Psi2S/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 553:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Upsilon1S/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 100553:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Upsilon2S/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case 200553:
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("Upsilon3S/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        default:
          break;
      }
    } else if (hfll_type > 0) { // HFll
      switch (hfll_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("ccbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB): // ULS
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB): // LS
          fRegistry.fill(HIST("Pair/") + HIST(pfNames[pftype]) + HIST("bbbar/") + HIST(dileptonSigns[signtype]) + HIST("hMvsPt"), v12.M(), v12.Pt());
        default:
          break;
      }
    }
  }

  template <int pftype, typename TMCParticle>
  void fillTrackInfo(TMCParticle const& mcparticle)
  {
    fRegistry.fill(HIST("Track/") + HIST(pfNames[pftype]) + HIST("hPt"), mcparticle.pt());
    fRegistry.fill(HIST("Track/") + HIST(pfNames[pftype]) + HIST("hEtaPhi"), mcparticle.phi(), mcparticle.eta());
  }

  template <typename TMCCollisions, typename TMCParticles>
  void runMC(TMCCollisions const& mcCollisions, TMCParticles const& mcParticles)
  {
    std::unordered_map<int, uint8_t> map_pfb; // global index of mc particle -> prefilter bit
    for (const auto& mcCollision : mcCollisions) {
      fRegistry.fill(HIST("Event/hZvtx"), mcCollision.posZ());
      fRegistry.fill(HIST("Event/hImpactParameter"), mcCollision.impactParameter());

      auto mcParticles_per_mccollision = mcParticles.sliceBy(perMcCollision, mcCollision.globalIndex());
      for (const auto& mcparticle : mcParticles_per_mccollision) {
        if (mcparticle.isPhysicalPrimary() || mcparticle.producedByGenerator()) {
          if (std::abs(mcparticle.pdgCode()) == 211) {
            fRegistry.fill(HIST("Event/hPtY_pion"), mcparticle.y(), mcparticle.pt());
          } else if (std::abs(mcparticle.pdgCode()) == 321) {
            fRegistry.fill(HIST("Event/hPtY_kaon"), mcparticle.y(), mcparticle.pt());
          } else if (std::abs(mcparticle.pdgCode()) == 2212) {
            fRegistry.fill(HIST("Event/hPtY_proton"), mcparticle.y(), mcparticle.pt());
          }
        }
      }

      // store MC true information
      auto posLeptons_per_mccollision = mcPosLeptons.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto negLeptons_per_mccollision = mcNegLeptons.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto posLeptonsPF_per_mccollision = mcPosLeptonsPF.sliceBy(perMcCollision, mcCollision.globalIndex());
      auto negLeptonsPF_per_mccollision = mcNegLeptonsPF.sliceBy(perMcCollision, mcCollision.globalIndex());

      for (const auto& pos : posLeptons_per_mccollision) {
        map_pfb[pos.globalIndex()] = 0;
        if (!isSelectedMCParticle(pos, mcParticles)) {
          continue;
        }
        fillTrackInfo<0>(pos);
        // fRegistry.fill(HIST("Track/default/hPt"), pos.pt());
        // fRegistry.fill(HIST("Track/default/hEtaPhi"), pos.phi(), pos.eta());
      }
      for (const auto& neg : negLeptons_per_mccollision) {
        map_pfb[neg.globalIndex()] = 0;
        if (!isSelectedMCParticle(neg, mcParticles)) {
          continue;
        }
        fillTrackInfo<0>(neg);
        // fRegistry.fill(HIST("Track/default/hPt"), neg.pt());
        // fRegistry.fill(HIST("Track/default/hEtaPhi"), neg.phi(), neg.eta());
      }

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posLeptonsPF_per_mccollision, negLeptonsPF_per_mccollision))) { // ULS to set prefilter bits
        if (!isSelectedMCParticle(pos, mcParticles) || !isSelectedMCParticle(neg, mcParticles)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(pos.pt(), pos.eta(), pos.phi(), leptonMass);
        ROOT::Math::PtEtaPhiMVector v2(neg.pt(), neg.eta(), neg.phi(), leptonMass);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/PF/all/uls/hMvsPt"), v12.M(), v12.Pt());
        for (int i = 0; i < static_cast<int>(max_mll_vec->size()); i++) {
          if (v12.M() < max_mll_vec->at(i)) {
            map_pfb[pos.globalIndex()] |= (uint8_t(1) << i);
            map_pfb[neg.globalIndex()] |= (uint8_t(1) << i);
          }
        }
      } // end of ULS pair loop to set prefilter bits

      for (const auto& pos : posLeptons_per_mccollision) {
        if (!isSelectedMCParticle(pos, mcParticles)) {
          continue;
        }
        static_for<0, 3 - 1>([&](auto i) {
          constexpr int index = i.value;
          if ((map_pfb[pos.globalIndex()] & (uint8_t(1) << index)) == 0) {
            fillTrackInfo<index + 1>(pos);
            // fRegistry.fill(HIST("Track/") + HIST(pfNames[index + 1]) + HIST("hPt"), pos.pt());
            // fRegistry.fill(HIST("Track/") + HIST(pfNames[index + 1]) + HIST("hEtaPhi"), pos.phi(), pos.eta());
          }
        });
      }
      for (const auto& neg : negLeptons_per_mccollision) {
        if (!isSelectedMCParticle(neg, mcParticles)) {
          continue;
        }
        static_for<0, 3 - 1>([&](auto i) {
          constexpr int index = i.value;
          if ((map_pfb[neg.globalIndex()] & (uint8_t(1) << index)) == 0) {
            fillTrackInfo<index + 1>(neg);
            // fRegistry.fill(HIST("Track/") + HIST(pfNames[index + 1]) + HIST("hPt"), neg.pt());
            // fRegistry.fill(HIST("Track/") + HIST(pfNames[index + 1]) + HIST("hEtaPhi"), neg.phi(), neg.eta());
          }
        });
      }

      for (const auto& [pos, neg] : combinations(CombinationsFullIndexPolicy(posLeptons_per_mccollision, negLeptons_per_mccollision))) { // ULS
        if (!isSelectedMCParticle(pos, mcParticles) || !isSelectedMCParticle(neg, mcParticles)) {
          continue;
        }
        // LOGF(info, "mcCollision.globalIndex() = %d,  pos.globalIndex() = %d, map_pfb[pos.globalIndex()] = %d, neg.globalIndex() = %d, map_pfb[neg.globalIndex()] = %d", mcCollision.globalIndex(), pos.globalIndex(), map_pfb[pos.globalIndex()], neg.globalIndex(), map_pfb[neg.globalIndex()]);

        fillTrueInfo<0, 0>(pos, neg, mcParticles); // default
        static_for<0, 3 - 1>([&](auto i) {
          constexpr int index = i.value;
          if ((map_pfb[pos.globalIndex()] & (uint8_t(1) << index)) == 0 && (map_pfb[neg.globalIndex()] & (uint8_t(1) << index)) == 0) {
            fillTrueInfo<index + 1, 0>(pos, neg, mcParticles); // mllPF0,1,2
          }
        });
      } // end of ULS pair loop

      for (auto& [pos1, pos2] : combinations(CombinationsStrictlyUpperIndexPolicy(posLeptons_per_mccollision, posLeptons_per_mccollision))) { // LS++
        if (!isSelectedMCParticle(pos1, mcParticles) || !isSelectedMCParticle(pos2, mcParticles)) {
          continue;
        }
        fillTrueInfo<0, 1>(pos1, pos2, mcParticles); // default
        static_for<0, 3 - 1>([&](auto i) {
          constexpr int index = i.value;
          if ((map_pfb[pos1.globalIndex()] & (uint8_t(1) << index)) == 0 && (map_pfb[pos2.globalIndex()] & (uint8_t(1) << index)) == 0) {
            fillTrueInfo<index + 1, 1>(pos1, pos2, mcParticles); // mllPF0,1,2
          }
        });
      } // end of LS++ pair loop

      for (auto& [neg1, neg2] : combinations(CombinationsStrictlyUpperIndexPolicy(negLeptons_per_mccollision, negLeptons_per_mccollision))) { // LS--
        if (!isSelectedMCParticle(neg1, mcParticles) || !isSelectedMCParticle(neg2, mcParticles)) {
          continue;
        }
        fillTrueInfo<0, 2>(neg1, neg2, mcParticles); // default
        static_for<0, 3 - 1>([&](auto i) {
          constexpr int index = i.value;
          if ((map_pfb[neg1.globalIndex()] & (uint8_t(1) << index)) == 0 && (map_pfb[neg2.globalIndex()] & (uint8_t(1) << index)) == 0) {
            fillTrueInfo<index + 1, 2>(neg1, neg2, mcParticles); // mllPF0,1,2
          }
        });
      } // end of LS-- pair loop

      map_pfb.clear();
    } // end of mc collision loop

  } //  end of skimmingMC

  SliceCache cache;
  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;

  Filter collisionFilter = min_imp_par < o2::aod::mccollision::impactParameter && o2::aod::mccollision::impactParameter < max_imp_par;
  using FilteredMcCollisions = soa::Filtered<aod::McCollisions>;

  Partition<aod::McParticles> mcPosLeptons = o2::aod::mcparticle::pdgCode == -pdg_lepton && (min_pt_gen < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < max_pt_gen) && (min_eta_gen < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen);
  Partition<aod::McParticles> mcNegLeptons = o2::aod::mcparticle::pdgCode == pdg_lepton && (min_pt_gen < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < max_pt_gen) && (min_eta_gen < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen);
  Partition<aod::McParticles> mcPosLeptonsPF = o2::aod::mcparticle::pdgCode == -pdg_lepton && (min_pt_gen_pf < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < max_pt_gen_pf) && (min_eta_gen_pf < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen_pf);
  Partition<aod::McParticles> mcNegLeptonsPF = o2::aod::mcparticle::pdgCode == pdg_lepton && (min_pt_gen_pf < o2::aod::mcparticle::pt && o2::aod::mcparticle::pt < max_pt_gen_pf) && (min_eta_gen_pf < o2::aod::mcparticle::eta && o2::aod::mcparticle::eta < max_eta_gen_pf);

  void processMC(FilteredMcCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    runMC(mcCollisions, mcParticles);
  }
  PROCESS_SWITCH(studyMCTruth, processMC, "process", true);

  void processDummy(FilteredMcCollisions const&) {}
  PROCESS_SWITCH(studyMCTruth, processDummy, "process Dummy", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<studyMCTruth>(cfgc, TaskName{"study-mc-truth"})};
}
