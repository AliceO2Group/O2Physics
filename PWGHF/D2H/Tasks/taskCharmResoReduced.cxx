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

/// \file taskCharmResoReduced.cxx
/// \brief Charmed Resonances analysis task
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"

// #include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::analysis;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

enum DecayTypeMc : uint8_t {
  Ds1ToDStarK0ToD0PiK0s = 1,
  Ds2StarToDplusK0sToPiKaPiPiPi,
  Ds1ToDStarK0ToDPlusPi0K0s,
  Ds1ToDStarK0ToD0PiK0sPart,
  Ds1ToDStarK0ToD0NoPiK0sPart,
  Ds1ToDStarK0ToD0PiK0sOneMu,
  Ds2StarToDplusK0sOneMu
};

namespace o2::aod
{
namespace hf_cand_reso_lite
{
DECLARE_SOA_COLUMN(PtBach0, ptBach0, float);                                               //! Transverse momentum of bachelor 0 (GeV/c)
DECLARE_SOA_COLUMN(PtBach1, ptBach1, float);                                               //! Transverse momentum of bachelor 1 (GeV/c)
DECLARE_SOA_COLUMN(MBach0, mBach0, float);                                                 //! Invariant mass of bachelor 0 (GeV/c)
DECLARE_SOA_COLUMN(MBach1, mBach1, float);                                                 //! Invariant mass of bachelor 1 (GeV/c)
DECLARE_SOA_COLUMN(MBachD0, mBachD0, float);                                               //! Invariant mass of D0 bachelor (of bachelor 0) (GeV/c)
DECLARE_SOA_COLUMN(M, m, float);                                                           //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(Pt, pt, float);                                                         //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(P, p, float);                                                           //! Momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(Y, y, float);                                                           //! Rapidity of candidate
DECLARE_SOA_COLUMN(Eta, eta, float);                                                       //! Pseudorapidity of candidate
DECLARE_SOA_COLUMN(Phi, phi, float);                                                       //! Azimuth angle of candidate
DECLARE_SOA_COLUMN(E, e, float);                                                           //! Energy of candidate (GeV)
DECLARE_SOA_COLUMN(Sign, sign, int8_t);                                                    //! Sign of candidate
DECLARE_SOA_COLUMN(CosThetaStar, cosThetaStar, float);                                     //! VosThetaStar of candidate (GeV)
DECLARE_SOA_COLUMN(MlScoreBkgBach0, mlScoreBkgBach0, float);                               //! ML score for background class of charm daughter
DECLARE_SOA_COLUMN(MlScorePromptBach0, mlScorePromptBach0, float);                         //! ML score for prompt class of charm daughter
DECLARE_SOA_COLUMN(MlScoreNonPromptBach0, mlScoreNonPromptBach0, float);                   //! ML score for non-prompt class of charm daughter
DECLARE_SOA_COLUMN(ItsNClsProngMinBach0, itsNClsProngMinBach0, int);                       //! minimum value of number of ITS clusters for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsProngMinBach0, tpcNClsCrossedRowsProngMinBach0, int); //! minimum value of number of TPC crossed rows for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(TpcChi2NClProngMaxBach0, tpcChi2NClProngMaxBach0, float);               //! maximum value of TPC chi2 for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(ItsNClsProngMinBach1, itsNClsProngMinBach1, int);                       //! minimum value of number of ITS clusters for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsProngMinBach1, tpcNClsCrossedRowsProngMinBach1, int); //! minimum value of number of TPC crossed rows for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(TpcChi2NClProngMaxBach1, tpcChi2NClProngMaxBach1, float);               //! maximum value of TPC chi2 for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(CpaBach1, cpaBach1, float);                                             //! Cosine of Pointing Angle of bachelor 1
DECLARE_SOA_COLUMN(DcaBach1, dcaBach1, float);                                             //! DCA of bachelor 1
DECLARE_SOA_COLUMN(RadiusBach1, radiusBach1, float);                                       //! Radius of bachelor 1
DECLARE_SOA_COLUMN(FlagMcMatchRec, flagMcMatchRec, int8_t);                                //! flag for decay channel classification reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int8_t);                                        //! debug flag for mis-association at reconstruction level
DECLARE_SOA_COLUMN(Origin, origin, int8_t);                                                //! Flag for origin of MC particle 1=promt, 2=FD
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                                   //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(SignD0, signD0, float);                                                 //! Flag to distinguish D0 and D0Bar

} // namespace hf_cand_reso_lite

DECLARE_SOA_TABLE(HfCandResoLites, "AOD", "HFCANDRESOLITE", //! Table with some B0 properties
                                                            // Candidate Properties
                  hf_cand_reso_lite::M,
                  hf_cand_reso_lite::Pt,
                  hf_cand_reso_lite::P,
                  hf_cand_reso_lite::Y,
                  hf_cand_reso_lite::Eta,
                  hf_cand_reso_lite::Phi,
                  hf_cand_reso_lite::E,
                  hf_cand_reso_lite::CosThetaStar,
                  hf_cand_reso_lite::Sign,
                  // Bachelors Properties
                  hf_cand_reso_lite::MBach0,
                  hf_cand_reso_lite::PtBach0,
                  hf_cand_reso_lite::MlScoreBkgBach0,
                  hf_cand_reso_lite::MlScorePromptBach0,
                  hf_cand_reso_lite::MlScoreNonPromptBach0,
                  hf_cand_reso_lite::ItsNClsProngMinBach0,
                  hf_cand_reso_lite::TpcNClsCrossedRowsProngMinBach0,
                  hf_cand_reso_lite::TpcChi2NClProngMaxBach0,
                  hf_cand_reso_lite::MBach1,
                  hf_cand_reso_lite::PtBach1,
                  hf_cand_reso_lite::CpaBach1,
                  hf_cand_reso_lite::DcaBach1,
                  hf_cand_reso_lite::RadiusBach1,
                  hf_cand_reso_lite::ItsNClsProngMinBach1,
                  hf_cand_reso_lite::TpcNClsCrossedRowsProngMinBach1,
                  hf_cand_reso_lite::TpcChi2NClProngMaxBach1,
                  // MC
                  hf_cand_reso_lite::FlagMcMatchRec,
                  hf_cand_reso_lite::DebugMcRec,
                  hf_cand_reso_lite::Origin,
                  hf_cand_reso_lite::PtGen,
                  hf_cand_reso_lite::SignD0);

DECLARE_SOA_TABLE(HfGenResoLites, "AOD", "HFGENRESOLITE", //! Table with some B0 properties
                  hf_cand_reso_lite::Pt,
                  hf_cand_reso_lite::Y,
                  hf_cand_reso_lite::Origin);

} // namespace o2::aod

enum DecayChannel : uint8_t {
  Ds1ToDstarK0s = 0,
  Ds2StarToDplusK0s,
  XcToDplusLambda,
  LambdaDminus
};

struct HfTaskCharmResoReduced {
  Produces<aod::HfCandResoLites> hfCandResoLite;
  Produces<aod::HfGenResoLites> hfGenResoLite;
  Configurable<float> ptMinReso{"ptMinReso", -1, "Discard events with smaller pT"};
  Configurable<bool> fillTrees{"fillTrees", true, "Fill output Trees"};
  Configurable<bool> fillSparses{"fillSparses", false, "Fill output Sparses"};
  Configurable<bool> useDeltaMass{"useDeltaMass", true, "Use Delta Mass for resonance invariant Mass calculation"};
  Configurable<bool> fillOnlySignal{"fillOnlySignal", false, "Flag to Fill only signal candidates (MC only)"};
  Configurable<float> yCandGenMax{"yCandGenMax", 0.5, "max. gen particle rapidity"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", -1, "max. cand. rapidity"};
  Configurable<float> etaTrackMax{"etaTrackMax", 0.8, "max. track pseudo-rapidity for acceptance calculation"};
  Configurable<float> ptTrackMin{"ptTrackMin", 0.1, "min. track transverse momentum for acceptance calculation"};
  Configurable<float> massResoMin{"massResoMin", 0.49, "min. mass of resonance"};
  Configurable<float> massResoMax{"massResoMax", 1.29, "max. mass of resonance"};
  // Configurables axis for histos
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisPtProng0{"axisPtProng0", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "prong0 bach. #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisPtProng1{"axisPtProng1", {VARIABLE_WIDTH, 0., 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 24.f, 50.f}, "prong1 bach. #it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axisInvMassReso{"axisInvMassReso", {200, 2.34, 2.74}, "inv. mass (DV_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng0{"axisInvMassProng0", {175, 1.70, 2.05}, "inv. mass (D) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisInvMassProng1{"axisInvMassProng1", {80, 0.46, 0.54}, "inv. mass ({V}_{0}) (GeV/#it{c}^{2})"};
  ConfigurableAxis axisCosThetaStar{"axisCosThetaStar", {40, -1, 1}, "cos(#vartheta*)"};
  ConfigurableAxis axisBkgBdtScore{"axisBkgBdtScore", {100, 0, 1}, "bkg BDT Score"};
  ConfigurableAxis axisNonPromptBdtScore{"axisNonPromptBdtScore", {100, 0, 1}, "non-prompt BDT Score"};
  ConfigurableAxis axisEta{"axisEta", {30, -1.5, 1.5}, "pseudorapidity"};
  ConfigurableAxis axisOrigin{"axisOrigin", {3, -0.5, 2.5}, "origin"};
  ConfigurableAxis axisFlag{"axisFlag", {65, -32.5, 32.5}, "mc flag"};

  using ReducedReso = soa::Join<aod::HfCandCharmReso, aod::HfResoIndices>;
  using ReducedResoWithMl = soa::Join<aod::HfCandCharmReso, aod::HfCharmResoMLs, aod::HfResoIndices>;
  using ReducedResoMc = soa::Join<aod::HfCandCharmReso, aod::HfResoIndices, aod::HfMcRecRedResos>;
  using ReducedResoWithMlMc = soa::Join<aod::HfCandCharmReso, aod::HfCharmResoMLs, aod::HfResoIndices, aod::HfMcRecRedResos>;

  // Histogram Registry
  HistogramRegistry registry;

  // init
  void init(InitContext&)
  {
    registry.add("hMass", "Charm resonance candidates inv. mass", {HistType::kTH1F, {axisInvMassReso}});
    registry.add("hMassProng0", "D daughters inv. mass", {HistType::kTH1F, {axisInvMassProng0}});
    registry.add("hMassProng1", "V0 daughter inv. mass", {HistType::kTH1F, {axisInvMassProng1}});
    registry.add("hPt", "Charm resonance candidates pT", {HistType::kTH1F, {axisPt}});
    registry.add("hPtProng0", "D daughters pT", {HistType::kTH1F, {axisPtProng0}});
    registry.add("hPtProng1", "V0 daughter pT", {HistType::kTH1F, {axisPtProng1}});
    registry.add("hNPvCont", "Collision number of PV contributors ; N contrib ; entries", {HistType::kTH1F, {{125, -0.5, 249.5}}});
    registry.add("hZvert", "Collision Z Vtx ; z PV [cm] ; entries", {HistType::kTH1F, {{120, -12., 12.}}});
    registry.add("hBz", "Collision Bz ; Bz [T] ; entries", {HistType::kTH1F, {{20, -10., 10.}}});
    registry.add("hSparse", "THn for production studies with cosThStar and BDT scores", HistType::kTHnSparseF, {axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar, axisBkgBdtScore, axisNonPromptBdtScore});

    if (doprocessDs1Mc || doprocessDs2StarMc || doprocessDs1McWithMl || doprocessDs2StarMcWithMl) {
      // gen histos
      registry.add("hYRecPrompt", "Charm resonance candidates pT", {HistType::kTH2F, {axisPt, axisEta}});
      registry.add("hYRecNonPrompt", "Charm resonance candidates pT", {HistType::kTH2F, {axisPt, axisEta}});
      registry.add("hYGenPrompt", "Prompt {D_{S}}^j particles (generated);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2F, {axisPt, axisEta}});
      registry.add("hYGenPromptWithProngsInAcceptance", "Prompt {D_{S}}^j particles (generated-daughters in acceptance);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2F, {axisPt, axisEta}});
      registry.add("hYGenNonPrompt", "NonPrompt {D_{S}}^j particles (generated);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2F, {axisPt, axisEta}});
      registry.add("hYGenNonPromptWithProngsInAcceptance", "NonPrompt {D_{S}}^j particles (generated-daughters in acceptance);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2F, {axisPt, axisEta}});
      if (fillSparses) {
        registry.add("hPtYGenSig", "{D_{S}}^j particles (generated);#it{p}_{T}({D_{S}}^j) (GeV/#it{c});#it{y}({D_{S}}^j)", {HistType::kTHnSparseF, {axisPt, axisEta, axisOrigin, axisFlag}});
        registry.add("hPtYWithProngsInAccepanceGenSig", "{D_{S}}^j particles (generated-daughters in acceptance);#it{p}_{T}({D_{S}}^j) (GeV/#it{c});#it{y}({D_{S}}^j)", {HistType::kTHnSparseF, {axisPt, axisEta, axisOrigin, axisFlag}});
      }
    }
  }

  // Fill histograms
  /// \tparam channel is the decay channel of the Resonance
  /// \param candidate is a candidate
  /// \param coll is a reduced collision
  /// \param bach0 is a bachelor of the candidate
  /// \param bach1 is a bachelor of the candidate
  template <bool doMc, bool withMl, DecayChannel channel, typename Cand, typename Coll, typename CharmBach, typename V0Bach>
  void fillCand(const Cand& candidate, const Coll& collision, const CharmBach& bach0, const V0Bach& bach1)
  {
    // Compute quantities to be saved
    float invMassReso{0}, pdgMassReso, invMassBach0, invMassBach1, pdgMassBach0, pdgMassBach1, sign, invMassD0, cosThetaStar;
    if (channel == DecayChannel::Ds1ToDstarK0s) {
      pdgMassReso = MassDS1;
      pdgMassBach0 = MassDStar;
      pdgMassBach1 = MassK0;
      invMassBach1 = bach1.invMassK0s();
      cosThetaStar = candidate.cosThetaStarDs1();
      if (bach0.dType() > 0) {
        invMassBach0 = bach0.invMassDstar();
        invMassD0 = bach0.invMassD0();
        sign = 1;
        if (useDeltaMass) {
          invMassReso = RecoDecay::m(std::array{bach0.pVectorProng0(), bach0.pVectorProng1(), bach0.pVectorProng2(), bach1.pVector()}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
        }
      } else {
        invMassBach0 = bach0.invMassAntiDstar();
        invMassD0 = bach0.invMassD0Bar();
        sign = -1;
        if (useDeltaMass) {
          invMassReso = RecoDecay::m(std::array{bach0.pVectorProng1(), bach0.pVectorProng0(), bach0.pVectorProng2(), bach1.pVector()}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
        }
      }
    } else if (channel == DecayChannel::Ds2StarToDplusK0s) {
      pdgMassReso = MassDS2Star;
      pdgMassBach0 = MassDPlus;
      pdgMassBach1 = MassK0;
      invMassBach0 = bach0.invMassDplus();
      invMassD0 = 0;
      invMassBach1 = bach1.invMassK0s();
      cosThetaStar = candidate.cosThetaStarDs2Star();
      if (useDeltaMass) {
        invMassReso = RecoDecay::m(std::array{bach0.pVectorProng0(), bach0.pVectorProng1(), bach0.pVectorProng2(), bach1.pVector()}, std::array{MassPiPlus, MassKPlus, MassPiPlus, MassK0});
      }
      if (bach0.dType() > 0) {
        sign = 1;
      } else {
        sign = -1;
      }
    }
    float y = RecoDecay::y(std::array{candidate.px(), candidate.py(), candidate.pz()}, pdgMassReso);
    float eta = RecoDecay::eta(std::array{candidate.px(), candidate.py(), candidate.pz()});
    float phi = RecoDecay::phi(candidate.px(), candidate.py());
    float p = RecoDecay::p(std::array{candidate.px(), candidate.py(), candidate.pz()});
    float e = RecoDecay::e(std::array{candidate.px(), candidate.py(), candidate.pz()}, pdgMassReso);
    if (useDeltaMass) {
      invMassReso = invMassReso - invMassBach0;
    } else {
      invMassReso = RecoDecay::m(std::array{bach0.pVector(), bach1.pVector()}, std::array{pdgMassBach0, pdgMassBach1});
    }
    if (invMassReso < massResoMin || invMassReso > massResoMax) {
      return;
    }
    invMassBach0 = invMassBach0 - invMassD0;
    float ptGen{-1.};
    int8_t origin{-1}, flagMcMatchRec{-1}, debugMcRec{-1}, signD0{0};
    if constexpr (doMc) {
      ptGen = candidate.ptMother();
      origin = candidate.origin();
      flagMcMatchRec = candidate.flagMcMatchRec();
      debugMcRec = candidate.debugMcRec();
      if (fillOnlySignal) {
        if (channel == DecayChannel::Ds1ToDstarK0s && !(std::abs(flagMcMatchRec) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0s || std::abs(flagMcMatchRec) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0sPart || std::abs(flagMcMatchRec) == DecayTypeMc::Ds1ToDStarK0ToD0NoPiK0sPart || std::abs(flagMcMatchRec) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0sOneMu)) {
          return;
        }
        if (channel == DecayChannel::Ds2StarToDplusK0s && !(std::abs(flagMcMatchRec) == DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi || std::abs(flagMcMatchRec) == DecayTypeMc::Ds2StarToDplusK0sOneMu)) {
          return;
        }
      }
      if (origin == 1) {
        registry.fill(HIST("hYRecPrompt"), candidate.pt(), y);
      } else if (origin == 2) {
        registry.fill(HIST("hYRecNonPrompt"), candidate.pt(), y);
      }
    }
    float mlScoreBkg{-1.}, mlScorePrompt{-1.}, mlScoreNonPrompt{-1.};
    if constexpr (withMl) {
      mlScoreBkg = bach0.mlScoreBkgMassHypo0();
      mlScorePrompt = bach0.mlScorePromptMassHypo0();
      mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo0();
    }
    // Collision properties
    registry.fill(HIST("hNPvCont"), collision.numContrib());
    registry.fill(HIST("hZvert"), collision.posZ());
    registry.fill(HIST("hBz"), collision.bz());
    // Candidate properties
    registry.fill(HIST("hMass"), invMassReso);
    registry.fill(HIST("hMassProng0"), invMassBach0);
    registry.fill(HIST("hMassProng1"), invMassBach1);
    registry.fill(HIST("hPt"), candidate.pt());
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    if (fillSparses) {
      registry.fill(HIST("hSparse"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), invMassReso, invMassBach0, invMassBach1, cosThetaStar, mlScoreBkg, mlScoreNonPrompt);
    }

    if (fillTrees) {
      hfCandResoLite(
        invMassReso,
        candidate.pt(),
        p,
        y,
        eta,
        phi,
        e,
        cosThetaStar,
        sign,
        // Bachelors Properties
        invMassBach0,
        bach0.pt(),
        mlScoreBkg,
        mlScorePrompt,
        mlScoreNonPrompt,
        bach0.itsNClsProngMin(),
        bach0.tpcNClsCrossedRowsProngMin(),
        bach0.tpcChi2NClProngMax(),
        invMassBach1,
        bach1.pt(),
        bach1.cpa(),
        bach1.dca(),
        bach1.v0Radius(),
        bach1.itsNClsProngMin(),
        bach1.tpcNClsCrossedRowsProngMin(),
        bach1.tpcChi2NClProngMax(),
        // MC
        flagMcMatchRec,
        debugMcRec,
        origin,
        ptGen,
        signD0);
    }
  } // fillCand

  // Process data
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param CharmBach is the reduced 3 prong table
  /// \param V0Bach is the reduced v0 table
  /// \param Cand is the candidates table
  template <bool doMc, bool withMl, DecayChannel channel, typename Coll, typename Candidates, typename CharmBach>
  void processData(Coll const&, Candidates const& candidates, CharmBach const&, aod::HfRedVzeros const&)
  {
    for (const auto& cand : candidates) {
      if (ptMinReso >= 0 && cand.pt() < ptMinReso) {
        continue;
      }
      float pdgMassReso{0};
      if (channel == DecayChannel::Ds1ToDstarK0s) {
        pdgMassReso = MassDS1;
      } else if (channel == DecayChannel::Ds2StarToDplusK0s) {
        pdgMassReso = MassDS2Star;
      }
      if (yCandRecoMax >= 0. && std::abs(RecoDecay::y(std::array{cand.px(), cand.py(), cand.pz()}, pdgMassReso)) > yCandRecoMax) {
        continue;
      }
      auto coll = cand.template hfRedCollision_as<Coll>();
      auto bach0 = cand.template prong0_as<CharmBach>();
      auto bach1 = cand.template prong1_as<aod::HfRedVzeros>();
      fillCand<doMc, withMl, channel>(cand, coll, bach0, bach1);
    }
  }

  /// Selection of resonance daughters in geometrical acceptance
  /// \param etaProng is the pseudorapidity of Resonance prong
  /// \param ptProng is the pT of Resonance prong
  /// \return true if prong is in geometrical acceptance
  template <typename T = float>
  bool isProngInAcceptance(const T& etaProng, const T& ptProng)
  {
    return std::abs(etaProng) <= etaTrackMax && ptProng >= ptTrackMin;
  }

  /// Fill particle histograms (gen MC truth)
  template <DecayChannel channel>
  void fillCandMcGen(aod::HfMcGenRedResos const& mcParticles)
  {
    for (const auto& particle : mcParticles) {
      auto ptParticle = particle.ptTrack();
      auto yParticle = particle.yTrack();
      auto originParticle = particle.origin();
      auto flag = particle.flagMcMatchGen();
      if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
        continue;
      }
      std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
      std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};
      bool prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);
      if ((channel == DecayChannel::Ds1ToDstarK0s && std::abs(flag) == DecayTypeMc::Ds1ToDStarK0ToD0PiK0s) ||
          (channel == DecayChannel::Ds2StarToDplusK0s && std::abs(flag) == DecayTypeMc::Ds2StarToDplusK0sToPiKaPiPiPi)) {
        if (originParticle == 1) { // prompt particles
          registry.fill(HIST("hYGenPrompt"), ptParticle, yParticle);
          if (prongsInAcc) {
            registry.fill(HIST("hYGenPromptWithProngsInAcceptance"), ptParticle, yParticle);
          }
        } else if (originParticle == 2) {
          registry.fill(HIST("hYGenNonPrompt"), ptParticle, yParticle);
          if (prongsInAcc) {
            registry.fill(HIST("hYGenNonPromptWithProngsInAcceptance"), ptParticle, yParticle);
          }
        }
      }
      if (fillSparses) {
        registry.fill(HIST("hPtYGenSig"), ptParticle, yParticle, originParticle, flag);
        if (prongsInAcc) {
          registry.fill(HIST("hPtYWithProngsInAccepanceGenSig"), ptParticle, yParticle, originParticle, flag);
        }
      }
      if (fillTrees) {
        hfGenResoLite(ptParticle, yParticle, originParticle);
      }
    }
  } // fillCandMcGen

  // process functions

  void processDs1Data(aod::HfRedCollisions const& collisions, ReducedReso const& candidates, aod::HfRed3PrNoTrks const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::Ds1ToDstarK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1Data, "Process data for Ds1 analysis without Ml", true);

  void processDs1DataWithMl(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates, soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::Ds1ToDstarK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1DataWithMl, "Process data for Ds1 analysis with Ml", false);

  void processDs2StarData(aod::HfRedCollisions const& collisions, ReducedReso const& candidates, aod::HfRed3PrNoTrks const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::Ds2StarToDplusK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarData, "Process data Ds2Star analysis without Ml", false);

  void processDs2StarDataWithMl(aod::HfRedCollisions const& collisions, ReducedResoWithMl const& candidates, soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::Ds2StarToDplusK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarDataWithMl, "Process data Ds2Star analysis with Ml", false);

  void processDs1Mc(aod::HfRedCollisions const& collisions, ReducedResoMc const& candidates, aod::HfMcGenRedResos const& mcParticles, aod::HfRed3PrNoTrks const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<true, false, DecayChannel::Ds1ToDstarK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::Ds1ToDstarK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1Mc, "Process Mc for Ds1 analysis without Ml", false);

  void processDs1McWithMl(aod::HfRedCollisions const& collisions, ReducedResoWithMlMc const& candidates, aod::HfMcGenRedResos const& mcParticles, soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<true, true, DecayChannel::Ds1ToDstarK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::Ds1ToDstarK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs1McWithMl, "Process Mc for Ds1 analysis with Ml", false);

  void processDs2StarMc(aod::HfRedCollisions const& collisions, ReducedResoMc const& candidates, aod::HfMcGenRedResos const& mcParticles, aod::HfRed3PrNoTrks const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<true, false, DecayChannel::Ds2StarToDplusK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::Ds2StarToDplusK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarMc, "Process Mc for Ds2Star analysis without Ml", false);

  void processDs2StarMcWithMl(aod::HfRedCollisions const& collisions, ReducedResoWithMlMc const& candidates, aod::HfMcGenRedResos const& mcParticles, soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs, aod::HfRedVzeros const& v0Bachs)
  {
    processData<true, true, DecayChannel::Ds2StarToDplusK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::Ds2StarToDplusK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoReduced, processDs2StarMcWithMl, "Process Mc for Ds2Star analysis with Ml", false);

}; // struct HfTaskCharmResoReduced
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmResoReduced>(cfgc)};
}
