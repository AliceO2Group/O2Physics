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
/// \file    alice3-dilepton.cxx
/// \author  s.scheid@cern.ch, daiki.sekihata@cern.ch
///

#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>
#include <Framework/O2DatabasePDGPlugin.h>

#include <Math/Vector4D.h>


#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct Alice3Dilepton {
  enum HFllType {
    kUndef = -1,
    kCe_Ce = 0,        // ULS
    kBe_Be = 1,        // ULS
    kBCe_BCe = 2,      // ULS
    kBCe_Be_SameB = 3, // ULS
    kBCe_Be_DiffB = 4, // LS
  };
  enum PairType {
    kULS = 0,
    kLSpp = 1,
    kLSnn = 2,
  };

  SliceCache cache_mc;
  SliceCache cache_rec;

  Service<o2::framework::O2DatabasePDG> inspdg;

  Configurable<int> pdg{"pdg", 11, "pdg code for analysis. dielectron:11, dimuon:13"};
  Configurable<bool> requireHFEid{"requireHFEid", true, "Require HFE identification for both leptons"};
  Configurable<float> ptMin{"pt-min", 0.f, "Lower limit in pT"};
  Configurable<float> ptMax{"pt-max", 5.f, "Upper limit in pT"};
  Configurable<float> etaMin{"eta-min", -5.f, "Lower limit in eta"};
  Configurable<float> etaMax{"eta-max", 5.f, "Upper limit in eta"};
  Configurable<bool> useGen{"use-gen", false, "Use generated (true) or smeared/reconstructed (false) values for fiducial cuts"};
  Configurable<bool> selectReconstructed{"selectReconstructed", true, "Select only reconstructed tracks (true) or ghosts (false)"};
  Configurable<float> nSigmaEleCutOuterTOF{"nSigmaEleCutOuterTOF", 3., "Electron inclusion in outer TOF"};
  Configurable<float> nSigmaEleCutInnerTOF{"nSigmaEleCutInnerTOF", 3., "Electron inclusion in inner TOF"};
  Configurable<float> nSigmaPionCutOuterTOF{"nSigmaPionCutOuterTOF", 3., "Pion exclusion in outer TOF"};
  Configurable<float> nSigmaPionCutInnerTOF{"nSigmaPionCutInnerTOF", 3., "Pion exclusion in inner TOF"};
  Configurable<float> nSigmaElectronRich{"nSigmaElectronRich", 3., "Electron inclusion RICH"};
  Configurable<float> nSigmaPionRich{"nSigmaPionRich", 3., "Pion exclusion RICH"};
  Configurable<int> otfConfig{"otfConfig", 0, "OTF configuration flag"};

  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisVx{100, -1, 1, "Vtx_{x}"};
    const AxisSpec axisVy{100, -1, 1, "Vtx_{y}"};
    const AxisSpec axisVz{100, -20, 20, "Vtx_{z}"};

    const AxisSpec axisM{500, 0, 5, "#it{m}_{ll} (GeV/#it{c}^{2})"};
    const AxisSpec axisPt{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisSigmaEl{200, -10, 10, "n#sigma_{El}"};
    const AxisSpec axisTrackLengthOuterTOF{300, 0., 300., "Track length (cm)"};
    const AxisSpec axisEta{1000, -5, 5, "#it{#eta}"};
    const AxisSpec axisDCAxy{1000, 0, 20, "DCA_{xy,ll} (#sigma)"};
    const AxisSpec axisPhi{360, 0, TMath::TwoPi(), "#it{#varphi} (rad.)"};
    const AxisSpec axisProdx{2000, -100, 100, "Prod. Vertex X (cm)"};
    const AxisSpec axisPrody{2000, -100, 100, "Prod. Vertex Y (cm)"};
    const AxisSpec axisProdz{2000, -100, 100, "Prod. Vertex Z (cm)"};

    if (doprocessGen) {
      registry.add("Generated/Event/VtxX", "Vertex X", kTH1F, {axisVx});
      registry.add("Generated/Event/VtxY", "Vertex Y", kTH1F, {axisVy});
      registry.add("Generated/Event/VtxZ", "Vertex Z", kTH1F, {axisVz});
      registry.add("Generated/Particle/Pt", "Particle Pt", kTH1F, {axisPt});
      registry.add("Generated/Particle/Eta", "Particle Eta", kTH1F, {axisEta});
      registry.add("Generated/Particle/Phi", "Particle Phi", kTH1F, {axisPhi});
      registry.add("Generated/Particle/Eta_Pt", "Eta vs. Pt", kTH2F, {axisPt, axisEta}, true);
      registry.add("Generated/Particle/prodVx", "Particle Prod. Vertex X", kTH1F, {axisProdx});
      registry.add("Generated/Particle/prodVy", "Particle Prod. Vertex Y", kTH1F, {axisPrody});
      registry.add("Generated/Particle/prodVz", "Particle Prod. Vertex Z", kTH1F, {axisProdz});
      registry.add("Generated/Particle/ParticlesPerEvent", "Particles per event", kTH1F, {{100, 0, 100}});
      registry.add("Generated/Particle/ParticlesFit", "Charged Particles in Fit acceptance per event", kTH1F, {{15000,0,15000}});

      registry.add("Generated/Pair/ULS/Tried", "Pair tries", kTH1F, {{10, -0.5, 9.5}});
      registry.add("Generated/Pair/ULS/Mass", "Pair Mass", kTH1F, {axisM});
      registry.add("Generated/Pair/ULS/Pt", "Pair Pt", kTH1F, {axisPt});
      registry.add("Generated/Pair/ULS/Eta", "Pair Eta", kTH1F, {axisEta});
      registry.add("Generated/Pair/ULS/Phi", "Pair Phi", kTH1F, {axisPhi});
      registry.add("Generated/Pair/ULS/Mass_Pt", "Pair Mass vs. Pt", kTH2F, {axisM, axisPt}, true);

      registry.addClone("Generated/Pair/ULS/", "Generated/Pair/LSpp/");
      registry.addClone("Generated/Pair/ULS/", "Generated/Pair/LSnn/");
    }

    if (doprocessRec) {
      registry.add("Reconstructed/Event/VtxX", "Vertex X", kTH1F, {axisVx});
      registry.add("Reconstructed/Event/VtxY", "Vertex Y", kTH1F, {axisVy});
      registry.add("Reconstructed/Event/VtxZ", "Vertex Z", kTH1F, {axisVz});

      registry.add("Reconstructed/Track/Pt", "Track Pt", kTH1F, {axisPt});
      registry.add("Reconstructed/Track/Eta", "Track Eta", kTH1F, {axisEta});
      registry.add("Reconstructed/Track/Phi", "Track Phi", kTH1F, {axisPhi});
      registry.add("Reconstructed/Track/Eta_Pt", "Eta vs. Pt", kTH2F, {axisPt, axisEta}, true);
      registry.add("Reconstructed/Track/SigmaOTofvspt", "Track #sigma oTOF", kTH2F, {axisPt, axisSigmaEl});
      registry.add("Reconstructed/Track/SigmaITofvspt", "Track #sigma iTOF", kTH2F, {axisPt, axisSigmaEl});
      registry.add("Reconstructed/Track/SigmaRichvspt", "Track #sigma RICH", kTH2F, {axisPt, axisSigmaEl});
      registry.add("Reconstructed/Track/outerTOFTrackLength", "Track length outer TOF", kTH1F, {axisTrackLengthOuterTOF});

      registry.addClone("Reconstructed/Track/", "Reconstructed/TrackPID/");

      registry.add("Reconstructed/Pair/ULS/Mass", "Pair Mass", kTH1F, {axisM});
      registry.add("Reconstructed/Pair/ULS/Pt", "Pair Pt", kTH1F, {axisPt});
      registry.add("Reconstructed/Pair/ULS/Eta", "Pair Eta", kTH1F, {axisEta});
      registry.add("Reconstructed/Pair/ULS/Phi", "Pair Phi", kTH1F, {axisPhi});
      registry.add("Reconstructed/Pair/ULS/Mass_Pt", "Pair Mass vs. Pt", kTH2F, {axisM, axisPt}, true);

      registry.addClone("Reconstructed/Pair/ULS/", "Reconstructed/Pair/LSpp/");
      registry.addClone("Reconstructed/Pair/ULS/", "Reconstructed/Pair/LSnn/");

      registry.add("ReconstructedFiltered/Pair/ULS/Mass", "Pair Mass", kTH1F, {axisM});
      registry.add("ReconstructedFiltered/Pair/ULS/Pt", "Pair Pt", kTH1F, {axisPt});
      registry.add("ReconstructedFiltered/Pair/ULS/Eta", "Pair Eta", kTH1F, {axisEta});
      registry.add("ReconstructedFiltered/Pair/ULS/Phi", "Pair Phi", kTH1F, {axisPhi});
      registry.add("ReconstructedFiltered/Pair/ULS/Mass_Pt", "Pair Mass vs. Pt", kTH2F, {axisM, axisPt}, true);

      registry.addClone("ReconstructedFiltered/Pair/ULS/", "ReconstructedFiltered/Pair/LSpp/");
      registry.addClone("ReconstructedFiltered/Pair/ULS/", "ReconstructedFiltered/Pair/LSnn/");

      HistogramConfigSpec hs_rec{HistType::kTHnSparseF, {axisM, axisPt, axisDCAxy}, 3};
      registry.add("Reconstructed/Pair/ULS/hs_rec", "", hs_rec);
      registry.add("Reconstructed/Pair/LSpp/hs_rec", "", hs_rec);
      registry.add("Reconstructed/Pair/LSnn/hs_rec", "", hs_rec);
      registry.get<THnSparse>(HIST("Reconstructed/Pair/ULS/hs_rec"))->Sumw2();
      registry.get<THnSparse>(HIST("Reconstructed/Pair/LSpp/hs_rec"))->Sumw2();
      registry.get<THnSparse>(HIST("Reconstructed/Pair/LSnn/hs_rec"))->Sumw2();
    }
  }

  template <typename TTrack>
  bool IsInAcceptance(TTrack const& track)
  {
    if (track.pt() < ptMin || ptMax < track.pt()) {
      return false;
    }
    if (track.eta() < etaMin || etaMax < track.eta()) {
      return false;
    }
    return true;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int IsSameMother(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    if (!p1.has_mothers())
      return -1;
    if (!p2.has_mothers())
      return -1;

    // LOGF(info,"original motherid1 = %d , motherid2 = %d", p1.mothersIds()[0], p2.mothersIds()[0]);

    int motherid1 = p1.mothersIds()[0];
    auto mother1 = mcparticles.iteratorAt(motherid1);
    int mother1_pdg = mother1.pdgCode();

    int motherid2 = p2.mothersIds()[0];
    auto mother2 = mcparticles.iteratorAt(motherid2);
    int mother2_pdg = mother2.pdgCode();

    // LOGF(info,"motherid1 = %d , motherid2 = %d", motherid1, motherid2);

    if (motherid1 != motherid2)
      return -1;
    if (mother1_pdg != mother2_pdg)
      return -1;

    if (std::abs(mother1_pdg) != 22        // photon
        && std::abs(mother1_pdg) != 111    // pi0
        && std::abs(mother1_pdg) != 221    // eta
        && std::abs(mother1_pdg) != 331    // eta'
        && std::abs(mother1_pdg) != 113    // rho
        && std::abs(mother1_pdg) != 223    // omega
        && std::abs(mother1_pdg) != 333    // phi
        && std::abs(mother1_pdg) != 443    // Jpsi
        && std::abs(mother1_pdg) != 100443 // psi2S
    ) {
      return -1;
    }

    // LOGF(info, "%d is found.", mother1_pdg);
    return motherid1;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int IsHFULS(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    // in total, 4 cases for ULS pairs
    // 0. prompt c->e+ and cbar->e-
    // 1. b->e- and bbar->e+ (different b and bbar)
    // 2. b->c->e+ and bbar->cbar->e- (different b and bbar)
    // 3. b->c->e+ and b->e- (1 same b (or bbar))
    if (!p1.has_mothers())
      return HFllType::kUndef;
    if (!p2.has_mothers())
      return HFllType::kUndef;

    int motherid_p1 = p1.mothersIds()[0];
    int motherid_p2 = p2.mothersIds()[0];
    if (motherid_p1 == motherid_p2) { // different mother
      return HFllType::kUndef;        // this never happens in correlated HF->ee decays
    }

    auto mother_p1 = mcparticles.iteratorAt(motherid_p1);
    auto mother_p2 = mcparticles.iteratorAt(motherid_p2);
    int mother1_pdg = mother_p1.pdgCode();
    int mother2_pdg = mother_p2.pdgCode();

    if (((500 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 599) || (5000 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 5999)) && ((500 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 599) || (5000 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 5999))) {
      return HFllType::kBe_Be; // bb->ee, decay type = 2
    }

    if (mother_p1.has_mothers() && mother_p2.has_mothers()) { // search for decay type 1,3,4
      int grand_motherid_p1 = mother_p1.mothersIds()[0];
      int grand_motherid_p2 = mother_p2.mothersIds()[0];
      auto grand_mother_p1 = mcparticles.iteratorAt(grand_motherid_p1);
      auto grand_mother_p2 = mcparticles.iteratorAt(grand_motherid_p2);
      int grand_mother1_pdg = grand_mother_p1.pdgCode();
      int grand_mother2_pdg = grand_mother_p2.pdgCode();

      if (((400 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 499) || (4000 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 4999)) && ((400 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 499) || (4000 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 4999))) { // mother is charm

        if (((500 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 599) || (5000 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 5999)) && ((500 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 599) || (5000 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 5999))) { // grand mother is beauty
          return kBCe_BCe;                                                                                                                                                                                                                                                                                                          // b->c->e and b->c->e, decay type = 1
        } else {
          return kCe_Ce; // prompt cc->ee, decay type = 0
        }
      }

      if (motherid_p1 == grand_motherid_p2 || grand_motherid_p1 == motherid_p2) {
        if (
          (((500 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 599) || (5000 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 5999)) && ((500 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 599) || (5000 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 5999))) ||
          (((500 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 599) || (5000 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 5999)) && ((500 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 599) || (5000 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 5999)))) {
          return HFllType::kBCe_Be_SameB; // b->c->e and c->e, decay type = 3
        }
      }
    }
    return HFllType::kUndef;
  }

  template <typename TMCParticle1, typename TMCParticle2, typename TMCParticles>
  int IsHFLS(TMCParticle1 const& p1, TMCParticle2 const& p2, TMCParticles const& mcparticles)
  {
    // in total, 1 case for LS pairs
    // 4. b->c->e+ and bbar->e+
    if (!p1.has_mothers())
      return HFllType::kUndef;
    if (!p2.has_mothers())
      return HFllType::kUndef;

    int motherid_p1 = p1.mothersIds()[0];
    int motherid_p2 = p2.mothersIds()[0];
    if (motherid_p1 == motherid_p2) { // different mother
      return HFllType::kUndef;        // this never happens in correlated HF->ee decays
    }

    auto mother_p1 = mcparticles.iteratorAt(motherid_p1);
    auto mother_p2 = mcparticles.iteratorAt(motherid_p2);
    int mother1_pdg = mother_p1.pdgCode();
    int mother2_pdg = mother_p2.pdgCode();

    if (mother_p1.has_mothers() && mother_p2.has_mothers()) { // search for decay type 4
      int grand_motherid_p1 = mother_p1.mothersIds()[0];
      int grand_motherid_p2 = mother_p2.mothersIds()[0];
      auto grand_mother_p1 = mcparticles.iteratorAt(grand_motherid_p1);
      auto grand_mother_p2 = mcparticles.iteratorAt(grand_motherid_p2);
      int grand_mother1_pdg = grand_mother_p1.pdgCode();
      int grand_mother2_pdg = grand_mother_p2.pdgCode();

      if (motherid_p1 != grand_motherid_p2 && grand_motherid_p1 != motherid_p2) { // different mother and grand mother
        if (
          (((500 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 599) || (5000 < std::abs(mother1_pdg) && std::abs(mother1_pdg) < 5999)) && ((500 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 599) || (5000 < std::abs(grand_mother2_pdg) && std::abs(grand_mother2_pdg) < 5999))) ||
          (((500 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 599) || (5000 < std::abs(mother2_pdg) && std::abs(mother2_pdg) < 5999)) && ((500 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 599) || (5000 < std::abs(grand_mother1_pdg) && std::abs(grand_mother1_pdg) < 5999)))) {
          return HFllType::kBCe_Be_DiffB; // b->c->e and c->e, decay type = 4
        }
      }
    }
    return HFllType::kUndef;
  }

  template <typename T1, typename T2>
  ROOT::Math::PtEtaPhiMVector buildPairDCA(T1 const& t1, T2 const& t2, float& pair_dca_xy)
  {

    const float dcaXY_t1 = t1.dcaXY();
    const float dcaXY_t2 = t2.dcaXY();
    const float dcaXY_res_t1 = sqrt(t1.cYY());
    const float dcaXY_res_t2 = sqrt(t2.cYY());

    pair_dca_xy = sqrt((dcaXY_t2 * dcaXY_t2 / dcaXY_res_t2 + dcaXY_t1 * dcaXY_t1 / dcaXY_res_t1) / 2.);
    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon); // reconstructed pt/eta/phi
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon); // reconstructed pt/eta/phi
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
    return v12;
  }

  template <PairType pairtype, typename TTracks, typename TMCTracks>
  void FillPairRecWithPrefilter(TTracks const& tracks1, TTracks const& tracks2, TMCTracks const& /*mcParticles*/)
  {
    std::vector<uint64_t> prefilteredTracks;
    if constexpr (pairtype == PairType::kULS) {
      for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!t1.has_mcParticle() || !t2.has_mcParticle()) {
          continue;
        }

        float pair_dca_xy = 999.f;
        ROOT::Math::PtEtaPhiMVector v12 = buildPairDCA(t1, t2, pair_dca_xy);
        // prefilter for low-mass pairs
        if (v12.M() > 0.10) {
          continue;
        }
        // prefilter small opening angle pairs
        if (std::cos(t1.phi() - t2.phi()) < 0.99) {
          continue;
        }
        prefilteredTracks.push_back(t1.globalIndex());
        prefilteredTracks.push_back(t2.globalIndex());

      } // end of unlike-sign pair loop

      for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        // Skipping tracks that are in the prefiltered list
        if (std::find(prefilteredTracks.begin(), prefilteredTracks.end(), t1.globalIndex()) != prefilteredTracks.end()) {
          continue;
        }
        if (std::find(prefilteredTracks.begin(), prefilteredTracks.end(), t2.globalIndex()) != prefilteredTracks.end()) {
          continue;
        }

        float pair_dca_xy = 999.f;
        ROOT::Math::PtEtaPhiMVector v12 = buildPairDCA(t1, t2, pair_dca_xy);

        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Mass"), v12.M());
        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Pt"), v12.Pt());
        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Eta"), v12.Eta());
        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/Mass_Pt"), v12.M(), v12.Pt());
        registry.fill(HIST("ReconstructedFiltered/Pair/ULS/hs_rec"), v12.M(), v12.Pt(), pair_dca_xy);
      }
    }
  }

  template <PairType pairtype, typename TTracks, typename TMCTracks>
  void FillPairRec(TTracks const& tracks1, TTracks const& tracks2, TMCTracks const& mcParticles)
  {
    if constexpr (pairtype == PairType::kULS) {
      for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!t1.has_mcParticle() || !t2.has_mcParticle()) {
          continue;
        }

        auto mct1 = t1.template mcParticle_as<aod::McParticles>();
        if (std::abs(mct1.pdgCode()) != pdg || !mct1.isPhysicalPrimary()) {
          continue;
        }

        auto mct2 = t2.template mcParticle_as<aod::McParticles>();
        if (std::abs(mct2.pdgCode()) != pdg || !mct2.isPhysicalPrimary()) {
          continue;
        }

        int motherid = IsSameMother(mct1, mct2, mcParticles);
        int hfee_type = IsHFULS(mct1, mct2, mcParticles);

        if (motherid < 0 && hfee_type == HFllType::kUndef) {
          continue;
        }
        // auto mother = mcparticles.iteratorAt(motherid);

        float pair_dca_xy = 999.f;
        ROOT::Math::PtEtaPhiMVector v12 = buildPairDCA(t1, t2, pair_dca_xy);

        registry.fill(HIST("Reconstructed/Pair/ULS/Mass"), v12.M());
        registry.fill(HIST("Reconstructed/Pair/ULS/Pt"), v12.Pt());
        registry.fill(HIST("Reconstructed/Pair/ULS/Eta"), v12.Eta());
        registry.fill(HIST("Reconstructed/Pair/ULS/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
        registry.fill(HIST("Reconstructed/Pair/ULS/Mass_Pt"), v12.M(), v12.Pt());
        registry.fill(HIST("Reconstructed/Pair/ULS/hs_rec"), v12.M(), v12.Pt(), pair_dca_xy);
      } // end of unlike-sign pair loop

    } else if constexpr (pairtype == PairType::kLSpp || pairtype == PairType::kLSnn) {
      for (auto& [t1, t2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!t1.has_mcParticle() || !t2.has_mcParticle()) {
          continue;
        }

        auto mct1 = t1.template mcParticle_as<aod::McParticles>();
        if (std::abs(mct1.pdgCode()) != pdg || !mct1.isPhysicalPrimary()) {
          continue;
        }

        auto mct2 = t2.template mcParticle_as<aod::McParticles>();
        if (std::abs(mct2.pdgCode()) != pdg || !mct2.isPhysicalPrimary()) {
          continue;
        }

        int motherid = -1;
        int hfee_type = IsHFLS(mct1, mct2, mcParticles);

        if (motherid < 0 && hfee_type == HFllType::kUndef) {
          continue;
        }
        // auto mother = mcparticles.iteratorAt(motherid);

        float pair_dca_xy = 999.f;
        ROOT::Math::PtEtaPhiMVector v12 = buildPairDCA(t1, t2, pair_dca_xy);

        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Reconstructed/Pair/LSpp/Mass"), v12.M());
          registry.fill(HIST("Reconstructed/Pair/LSpp/Pt"), v12.Pt());
          registry.fill(HIST("Reconstructed/Pair/LSpp/Eta"), v12.Eta());
          registry.fill(HIST("Reconstructed/Pair/LSpp/Mass_Pt"), v12.M(), v12.Pt());
          registry.fill(HIST("Reconstructed/Pair/LSpp/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
          registry.fill(HIST("Reconstructed/Pair/LSpp/hs_rec"), v12.M(), v12.Pt(), pair_dca_xy);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Reconstructed/Pair/LSnn/Mass"), v12.M());
          registry.fill(HIST("Reconstructed/Pair/LSnn/Pt"), v12.Pt());
          registry.fill(HIST("Reconstructed/Pair/LSnn/Eta"), v12.Eta());
          registry.fill(HIST("Reconstructed/Pair/LSnn/Mass_Pt"), v12.M(), v12.Pt());
          registry.fill(HIST("Reconstructed/Pair/LSnn/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
          registry.fill(HIST("Reconstructed/Pair/LSnn/hs_rec"), v12.M(), v12.Pt(), pair_dca_xy);
        }
      } // end of like-sign pair loop
    }
  }

  template <PairType pairtype, typename TTracks, typename TMCTracks>
  void FillPairGen(TTracks const& tracks1, TTracks const& tracks2, TMCTracks const& mcParticles)
  {
    if constexpr (pairtype == PairType::kULS) {
      for (auto& [t1, t2] : combinations(soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        registry.fill(HIST("Generated/Pair/ULS/Tried"), 0);
        if (std::abs(t1.pdgCode()) != pdg || std::abs(t2.pdgCode()) != pdg) {
          continue;
        }
        registry.fill(HIST("Generated/Pair/ULS/Tried"), 1);

        if (!t1.isPhysicalPrimary() || !t2.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("Generated/Pair/ULS/Tried"), 2);

        if (!IsInAcceptance(t1) || !IsInAcceptance(t2)) {
          continue;
        }
        registry.fill(HIST("Generated/Pair/ULS/Tried"), 3);

        const int motherid = IsSameMother(t1, t2, mcParticles);
        const int hfee_type = IsHFULS(t1, t2, mcParticles);

        if (requireHFEid.value && motherid < 0 && hfee_type == HFllType::kUndef) {
          continue;
        }
        registry.fill(HIST("Generated/Pair/ULS/Tried"), 4);
        // auto mother = mcparticles.iteratorAt(motherid);

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        registry.fill(HIST("Generated/Pair/ULS/Mass"), v12.M());
        registry.fill(HIST("Generated/Pair/ULS/Pt"), v12.Pt());
        registry.fill(HIST("Generated/Pair/ULS/Eta"), v12.Eta());
        registry.fill(HIST("Generated/Pair/ULS/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
        registry.fill(HIST("Generated/Pair/ULS/Mass_Pt"), v12.M(), v12.Pt());
      } // end of unlike-sign pair loop

    } else if constexpr (pairtype == PairType::kLSpp || pairtype == PairType::kLSnn) {
      for (auto& [t1, t2] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Tried"), 0);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Tried"), 0);
        }
        if (std::abs(t1.pdgCode()) != pdg || std::abs(t2.pdgCode()) != pdg) {
          continue;
        }
        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Tried"), 1);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Tried"), 1);
        }

        if (!t1.isPhysicalPrimary() || !t2.isPhysicalPrimary()) {
          continue;
        }
        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Tried"), 2);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Tried"), 2);
        }

        if (!IsInAcceptance(t1) || !IsInAcceptance(t2)) {
          continue;
        }
        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Tried"), 3);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Tried"), 3);
        }

        const int motherid = -1;
        const int hfee_type = IsHFLS(t1, t2, mcParticles);

        if (requireHFEid.value && motherid < 0 && hfee_type == HFllType::kUndef) {
          continue;
        }
        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Tried"), 4);
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Tried"), 4);
        }
        // auto mother = mcparticles.iteratorAt(motherid);

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), std::abs(pdg) == 11 ? o2::constants::physics::MassElectron : o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        if constexpr (pairtype == PairType::kLSpp) {
          registry.fill(HIST("Generated/Pair/LSpp/Mass"), v12.M());
          registry.fill(HIST("Generated/Pair/LSpp/Pt"), v12.Pt());
          registry.fill(HIST("Generated/Pair/LSpp/Eta"), v12.Eta());
          registry.fill(HIST("Generated/Pair/LSpp/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
          registry.fill(HIST("Generated/Pair/LSpp/Mass_Pt"), v12.M(), v12.Pt());
        } else if constexpr (pairtype == PairType::kLSnn) {
          registry.fill(HIST("Generated/Pair/LSnn/Mass"), v12.M());
          registry.fill(HIST("Generated/Pair/LSnn/Pt"), v12.Pt());
          registry.fill(HIST("Generated/Pair/LSnn/Eta"), v12.Eta());
          registry.fill(HIST("Generated/Pair/LSnn/Phi"), v12.Phi() < 0.f ? v12.Phi() + TMath::TwoPi() : v12.Phi());
          registry.fill(HIST("Generated/Pair/LSnn/Mass_Pt"), v12.M(), v12.Pt());
        }

      } // end of like-sign pair loop
    }
  }

  // Functions for pid
  template <typename TTrack>
  bool electronIDTOF(TTrack const& track)
  {
    bool isElectron = false;
    bool isEleOuterTOF = std::abs(track.nSigmaElectronOuterTOF()) < nSigmaEleCutOuterTOF;
    bool isNotPionOuterTOF = std::abs(track.nSigmaPionOuterTOF()) > nSigmaPionCutOuterTOF;
    isEleOuterTOF = isEleOuterTOF && isNotPionOuterTOF;
    bool isEleInnerTOF = std::abs(track.nSigmaElectronInnerTOF()) < nSigmaEleCutInnerTOF;
    bool isNotPionInnerTOF = std::abs(track.nSigmaPionInnerTOF()) > nSigmaPionCutInnerTOF;
    isEleInnerTOF = isEleInnerTOF && isNotPionInnerTOF;
    isElectron = (isEleOuterTOF || isEleInnerTOF);
    return isElectron;
  }

  template <typename TTrack>
  bool electronIDRICH(TTrack const& track)
  {
    bool isElectron = false;
    bool isEleRICH = std::abs(track.nSigmaElectronRich()) < nSigmaElectronRich;
    bool isNotPionRICH = std::abs(track.nSigmaPionRich()) > nSigmaPionRich;
    isElectron = isEleRICH && isNotPionRICH;
    return isElectron;
  }

  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Partition<aod::McParticles> pos_mcParticles = o2::aod::mcparticle::pdgCode == -pdg; //-11 or -13
  Partition<aod::McParticles> neg_mcParticles = o2::aod::mcparticle::pdgCode == pdg;  // 11 or 13

  void processGen(o2::aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles)
  {
    for (const auto& mccollision : mccollisions) {
      registry.fill(HIST("Generated/Event/VtxX"), mccollision.posX());
      registry.fill(HIST("Generated/Event/VtxY"), mccollision.posY());
      registry.fill(HIST("Generated/Event/VtxZ"), mccollision.posZ());

      auto mcParticles_per_coll = mcParticles.sliceBy(perMCCollision, mccollision.globalIndex());
      int nParticlesInEvent = 0;
      int nParticlesFIT=0;
      for (const auto& mcParticle : mcParticles_per_coll) {
        if (mcParticle.isPhysicalPrimary()) {
          if ((2.2<mcParticle.eta()&&mcParticle.eta()< 5.0) || (-3.4<mcParticle.eta()&&mcParticle.eta()<-2.3)) {
            auto pdgParticle = inspdg->GetParticle(mcParticle.pdgCode());
            if (pdgParticle) {
              float charge = pdgParticle->Charge() / 3.f; // Charge in units of |e|
              if (std::abs(charge) >= 1.) {
                nParticlesFIT++;
              }
            }
          }
        }
        if (std::abs(mcParticle.pdgCode()) != pdg) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (!IsInAcceptance(mcParticle)) {
          continue;
        }
        nParticlesInEvent++;

        registry.fill(HIST("Generated/Particle/Pt"), mcParticle.pt());
        registry.fill(HIST("Generated/Particle/Eta"), mcParticle.eta());
        registry.fill(HIST("Generated/Particle/Phi"), mcParticle.phi());
        registry.fill(HIST("Generated/Particle/Eta_Pt"), mcParticle.pt(), mcParticle.eta());

        registry.fill(HIST("Generated/Particle/prodVx"), mcParticle.vx());
        registry.fill(HIST("Generated/Particle/prodVy"), mcParticle.vy());
        registry.fill(HIST("Generated/Particle/prodVz"), mcParticle.vz());

      } // end of mc particle loop
      registry.fill(HIST("Generated/Particle/ParticlesPerEvent"), nParticlesInEvent);
      registry.fill(HIST("Generated/Particle/ParticlesFIT"), nParticlesFIT);

      auto neg_mcParticles_coll = neg_mcParticles->sliceByCached(o2::aod::mcparticle::mcCollisionId, mccollision.globalIndex(), cache_mc);
      auto pos_mcParticles_coll = pos_mcParticles->sliceByCached(o2::aod::mcparticle::mcCollisionId, mccollision.globalIndex(), cache_mc);

      FillPairGen<PairType::kULS>(neg_mcParticles_coll, pos_mcParticles_coll, mcParticles);
      FillPairGen<PairType::kLSpp>(pos_mcParticles_coll, pos_mcParticles_coll, mcParticles);
      FillPairGen<PairType::kLSnn>(neg_mcParticles_coll, neg_mcParticles_coll, mcParticles);

    } // end of mc collision loop
  } // end of processGen

  using MyTracksMC = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::McTrackLabels, aod::UpgradeTofs, aod::UpgradeTofMCs, aod::UpgradeRichs, aod::TracksAlice3>;
  using Alice3Collision = soa::Join<aod::Collisions, aod::OTFLUTConfigId>;

  // Filter trackFilter = etaMin < o2::aod::track::eta &&
  //                      o2::aod::track::eta < etaMax &&
  //                      ptMin < o2::aod::track::pt &&
  //                      o2::aod::track::pt < ptMax &&
  //                      o2::aod::track_alice3::isReconstructed == selectReconstructed;
  Filter trackFilter = o2::aod::track_alice3::isReconstructed == selectReconstructed;
  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Filter configFilter = (aod::upgrade_collision::lutConfigId == otfConfig);
  using MyFilteredAlice3Collision = soa::Filtered<Alice3Collision>;
  Preslice<MyFilteredTracksMC> perCollision = aod::track::collisionId;
  Partition<MyFilteredTracksMC> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracks = o2::aod::track::signed1Pt < 0.f;

  void processRec(MyFilteredAlice3Collision const& collisions,
                  MyFilteredTracksMC const& tracks,
                  const o2::aod::McCollisions&,
                  const aod::McParticles& mcParticles)
  {
    for (const auto& collision : collisions) {
      registry.fill(HIST("Reconstructed/Event/VtxX"), collision.posX());
      registry.fill(HIST("Reconstructed/Event/VtxY"), collision.posY());
      registry.fill(HIST("Reconstructed/Event/VtxZ"), collision.posZ());

      auto tracks_coll = tracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& track : tracks_coll) {
        if (!track.has_mcParticle()) {
          continue;
        }
        const auto mcParticle = track.mcParticle_as<aod::McParticles>();
        if (std::abs(mcParticle.pdgCode()) != pdg) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        if (useGen) {
          if (!IsInAcceptance(mcParticle)) {
            continue;
          }
        } else {
          if (!IsInAcceptance(track)) {
            continue;
          }
        }
        if (std::abs(mcParticle.pdgCode()) != pdg) {
          continue;
        }
        if (!mcParticle.isPhysicalPrimary()) {
          continue;
        }
        registry.fill(HIST("Reconstructed/Track/SigmaOTofvspt"), track.pt(), track.nSigmaElectronOuterTOF());
        registry.fill(HIST("Reconstructed/Track/SigmaITofvspt"), track.pt(), track.nSigmaElectronInnerTOF());
        registry.fill(HIST("Reconstructed/Track/SigmaRichvspt"), track.pt(), track.nSigmaElectronRich());
        registry.fill(HIST("Reconstructed/Track/outerTOFTrackLength"), track.outerTOFTrackLength());
        registry.fill(HIST("Reconstructed/Track/Pt"), track.pt());
        registry.fill(HIST("Reconstructed/Track/Eta"), track.eta());
        registry.fill(HIST("Reconstructed/Track/Phi"), track.phi());
        registry.fill(HIST("Reconstructed/Track/Eta_Pt"), track.pt(), track.eta());
        // implement pid

        bool isElectronTOF = electronIDTOF(track);
        bool isElectronRICH = electronIDRICH(track);

        if (isElectronTOF || isElectronRICH) {
          registry.fill(HIST("Reconstructed/TrackPID/SigmaOTofvspt"), track.pt(), track.nSigmaElectronOuterTOF());
          registry.fill(HIST("Reconstructed/TrackPID/SigmaITofvspt"), track.pt(), track.nSigmaElectronInnerTOF());
          registry.fill(HIST("Reconstructed/TrackPID/SigmaRichvspt"), track.pt(), track.nSigmaElectronRich());
          registry.fill(HIST("Reconstructed/TrackPID/outerTOFTrackLength"), track.outerTOFTrackLength());
          registry.fill(HIST("Reconstructed/TrackPID/Pt"), track.pt());
          registry.fill(HIST("Reconstructed/TrackPID/Eta"), track.eta());
          registry.fill(HIST("Reconstructed/TrackPID/Phi"), track.phi());
          registry.fill(HIST("Reconstructed/TrackPID/Eta_Pt"), track.pt(), track.eta());
        }
      } // end of track loop
      // How is PID selection applied to tracks?
      auto negTracks_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache_rec);
      auto posTracks_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache_rec);

      FillPairRec<PairType::kULS>(negTracks_coll, posTracks_coll, mcParticles);
      FillPairRec<PairType::kLSpp>(posTracks_coll, posTracks_coll, mcParticles);
      FillPairRec<PairType::kLSnn>(negTracks_coll, negTracks_coll, mcParticles);

    } // end of collision loop
  } // end of processRec

  PROCESS_SWITCH(Alice3Dilepton, processGen, "Run for generated particle", true);
  PROCESS_SWITCH(Alice3Dilepton, processRec, "Run for reconstructed track", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3Dilepton>(cfgc, TaskName{"alice3-dilepton"})};
}
