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

/// \file taskCharmResoToDV0Reduced.cxx
/// \brief Charmed Resonances decaying in a D meson and a V0 analysis task
///
/// \author Luca Aglietta <luca.aglietta@cern.ch>, University and INFN Torino

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/D2H/DataModel/ReducedDataModel.h"
#include "PWGHF/Utils/utilsMcMatching.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Rtypes.h>

#include <array>
#include <cstdint>
#include <cstdlib>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

enum BachelorType : uint8_t {
  K0s = 0,
  Lambda,
  AntiLambda
};

namespace o2::aod
{
namespace hf_cand_reso_to_v0_lite
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
DECLARE_SOA_COLUMN(ItsNClsSoftPi, itsNClsSoftPi, int);                                     //! minimum value of number of ITS clusters for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsSoftPi, tpcNClsCrossedRowsSoftPi, int);               //! minimum value of number of TPC crossed rows for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(TpcChi2NClSoftPi, tpcChi2NClSoftPi, float);                             //! maximum value of TPC chi2 for the decay daughter tracks of bachelor 0
DECLARE_SOA_COLUMN(ItsNClsProngMinBach1, itsNClsProngMinBach1, int);                       //! minimum value of number of ITS clusters for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(TpcNClsCrossedRowsProngMinBach1, tpcNClsCrossedRowsProngMinBach1, int); //! minimum value of number of TPC crossed rows for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(TpcChi2NClProngMaxBach1, tpcChi2NClProngMaxBach1, float);               //! maximum value of TPC chi2 for the decay daughter tracks of bachelor 1
DECLARE_SOA_COLUMN(CpaBach1, cpaBach1, float);                                             //! Cosine of Pointing Angle of bachelor 1
DECLARE_SOA_COLUMN(DcaBach1, dcaBach1, float);                                             //! DCA of bachelor 1
DECLARE_SOA_COLUMN(RadiusBach1, radiusBach1, float);                                       //! Radius of bachelor 1
DECLARE_SOA_COLUMN(FlagMcMatch, flagMcMatch, int8_t);                                      //! flag for decay channel classification reconstruction level
DECLARE_SOA_COLUMN(DebugMcRec, debugMcRec, int);                                           //! debug flag for mis-association at reconstruction level
DECLARE_SOA_COLUMN(Origin, origin, int8_t);                                                //! Flag for origin of MC particle 1=promt, 2=FD
DECLARE_SOA_COLUMN(PtGen, ptGen, float);                                                   //! Transverse momentum of candidate (GeV/c)
DECLARE_SOA_COLUMN(InvMassGen, invMassGen, float);                                         //! Invariant mass of candidate (GeV/c2)
DECLARE_SOA_COLUMN(FlagCharmBach, flagCharmBach, int8_t);                                  //! Flag for charm bachelor classification
DECLARE_SOA_COLUMN(FlagCharmBachInterm, flagCharmBachInterm, int8_t);                      //! Flag for charm bachelor classification intermediate
DECLARE_SOA_COLUMN(NKinkedTracks, nKinkedTracks, int8_t);                                  //! Number of kinked tracks found in MC matching
} // namespace hf_cand_reso_to_v0_lite

DECLARE_SOA_TABLE(HfCandDV0Lites, "AOD", "HFCANDDV0LITE", //! Table with some Resonances properties
                                                          // Candidate Properties
                  hf_cand_reso_to_v0_lite::M,
                  hf_cand_reso_to_v0_lite::Pt,
                  hf_cand_reso_to_v0_lite::P,
                  hf_cand_reso_to_v0_lite::Y,
                  hf_cand_reso_to_v0_lite::Eta,
                  hf_cand_reso_to_v0_lite::Phi,
                  hf_cand_reso_to_v0_lite::E,
                  hf_cand_reso_to_v0_lite::CosThetaStar,
                  hf_cand_reso_to_v0_lite::Sign,
                  // Bachelors Properties
                  hf_cand_reso_to_v0_lite::MBach0,
                  hf_cand_reso_to_v0_lite::PtBach0,
                  hf_cand_reso_to_v0_lite::MlScoreBkgBach0,
                  hf_cand_reso_to_v0_lite::MlScorePromptBach0,
                  hf_cand_reso_to_v0_lite::MlScoreNonPromptBach0,
                  hf_cand_reso_to_v0_lite::ItsNClsProngMinBach0,
                  hf_cand_reso_to_v0_lite::TpcNClsCrossedRowsProngMinBach0,
                  hf_cand_reso_to_v0_lite::TpcChi2NClProngMaxBach0,
                  hf_cand_reso_to_v0_lite::ItsNClsSoftPi,
                  hf_cand_reso_to_v0_lite::TpcNClsCrossedRowsSoftPi,
                  hf_cand_reso_to_v0_lite::TpcChi2NClSoftPi,
                  hf_cand_reso_to_v0_lite::MBach1,
                  hf_cand_reso_to_v0_lite::PtBach1,
                  hf_cand_reso_to_v0_lite::CpaBach1,
                  hf_cand_reso_to_v0_lite::DcaBach1,
                  hf_cand_reso_to_v0_lite::RadiusBach1,
                  hf_cand_reso_to_v0_lite::ItsNClsProngMinBach1,
                  hf_cand_reso_to_v0_lite::TpcNClsCrossedRowsProngMinBach1,
                  hf_cand_reso_to_v0_lite::TpcChi2NClProngMaxBach1,
                  // MC
                  hf_cand_reso_to_v0_lite::FlagMcMatch,
                  hf_cand_reso_to_v0_lite::DebugMcRec,
                  hf_cand_reso_to_v0_lite::Origin,
                  hf_cand_reso_to_v0_lite::PtGen,
                  hf_cand_reso_to_v0_lite::InvMassGen,
                  hf_cand_reso_to_v0_lite::FlagCharmBach,
                  hf_cand_reso_to_v0_lite::FlagCharmBachInterm,
                  hf_cand_reso_to_v0_lite::NKinkedTracks);

DECLARE_SOA_TABLE(HfGenResoLites, "AOD", "HFGENRESOLITE", //! Table with some B0 properties
                  hf_cand_reso_to_v0_lite::Pt,
                  hf_cand_reso_to_v0_lite::Y,
                  hf_cand_reso_to_v0_lite::Origin,
                  hf_cand_reso_to_v0_lite::FlagMcMatch);

} // namespace o2::aod

enum DecayChannel : uint8_t {
  DstarK0s = 0,
  DplusK0s,
  DplusLambda,
  D0Lambda
};

struct HfTaskCharmResoToDV0Reduced {
  Produces<aod::HfCandDV0Lites> hfCandResoLite;
  Produces<aod::HfGenResoLites> hfGenResoLite;

  Configurable<bool> doWrongSign{"doWrongSign", false, "Flag to enable wrong sign candidates"};
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

  using ReducedReso3PrV0 = soa::Join<aod::HfCandCharmReso, aod::Hf3PrV0Ids>;
  using ReducedResoDstarV0 = soa::Join<aod::HfCandCharmReso, aod::HfDstarV0Ids>;
  using ReducedReso2PrV0 = soa::Join<aod::HfCandCharmReso, aod::Hf2PrV0Ids>;
  using ReducedReso3PrV0MC = soa::Join<aod::HfCandCharmReso, aod::Hf3PrV0Ids, aod::HfMcRecRedResos>;
  using ReducedResoDstarV0MC = soa::Join<aod::HfCandCharmReso, aod::HfDstarV0Ids, aod::HfMcRecRedResos>;
  using ReducedReso2PrV0MC = soa::Join<aod::HfCandCharmReso, aod::Hf2PrV0Ids, aod::HfMcRecRedResos>;

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

  // Histogram Registry
  HistogramRegistry registry;

  // init
  void init(InitContext&)
  {
    registry.add("hMass", "Charm resonance candidates inv. mass", {HistType::kTH1D, {axisInvMassReso}});
    registry.add("hMassProng0", "D daughters inv. mass", {HistType::kTH1D, {axisInvMassProng0}});
    registry.add("hMassProng1", "V0 daughter inv. mass", {HistType::kTH1D, {axisInvMassProng1}});
    registry.add("hPt", "Charm resonance candidates pT", {HistType::kTH1D, {axisPt}});
    registry.add("hPtProng0", "D daughters pT", {HistType::kTH1D, {axisPtProng0}});
    registry.add("hPtProng1", "V0 daughter pT", {HistType::kTH1D, {axisPtProng1}});
    registry.add("hNPvCont", "Collision number of PV contributors ; N contrib ; entries", {HistType::kTH1D, {{125, -0.5, 249.5}}});
    registry.add("hZvert", "Collision Z Vtx ; z PV [cm] ; entries", {HistType::kTH1D, {{120, -12., 12.}}});
    registry.add("hBz", "Collision Bz ; Bz [T] ; entries", {HistType::kTH1D, {{20, -10., 10.}}});
    registry.add("hSparse", "THn for production studies with cosThStar and BDT scores", HistType::kTHnSparseF, {axisPt, axisPtProng0, axisPtProng1, axisInvMassReso, axisInvMassProng0, axisInvMassProng1, axisCosThetaStar, axisBkgBdtScore, axisNonPromptBdtScore});

    if (doprocessDstarK0sMC || doprocessDplusK0sMC || doprocessDstarK0sMCWithMl || doprocessDplusK0sMCWithMl ||
        doprocessDplusLambdaMC || doprocessD0LambdaMC || doprocessDplusLambdaMCWithMl || doprocessD0LambdaMCWithMl) {
      // gen histos
      registry.add("hYRecPrompt", "Charm resonance candidates pT", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYRecNonPrompt", "Charm resonance candidates pT", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYGenAll", "Prompt {D_{S}}^j particles (generated);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYGenPrompt", "Prompt {D_{S}}^j particles (generated);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYGenPromptWithProngsInAcceptance", "Prompt {D_{S}}^j particles (generated-daughters in acceptance);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYGenNonPrompt", "NonPrompt {D_{S}}^j particles (generated);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2D, {axisPt, axisEta}});
      registry.add("hYGenNonPromptWithProngsInAcceptance", "NonPrompt {D_{S}}^j particles (generated-daughters in acceptance);#it{p}_{T}^{gen}({D_{S}}^j) (GeV/#it{c});#it{y}^{gen}({D_{S}}^j);entries", {HistType::kTH2D, {axisPt, axisEta}});
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
  template <bool DoMc, bool WithMl, DecayChannel Channel, typename Cand, typename Coll, typename CharmBach, typename V0Bach>
  void fillCand(const Cand& candidate, const Coll& collision, const CharmBach& bach0, const V0Bach& bach1)
  {
    // Base
    float massReso{0}, cosThetaStar{0};
    int8_t sign{0};
    int itsNClsSoftPi{0}, tpcNClsCrossedRowsSoftPi{0};
    float tpcChi2NClSoftPi{0.};
    if constexpr (Channel == DecayChannel::DstarK0s) {
      sign = bach0.sign();
      massReso = useDeltaMass ? candidate.invMass() + MassDStar : candidate.invMass();
      cosThetaStar = RecoDecay::cosThetaStar(std::array{bach0.pVector(), bach1.pVector()}, std::array{MassDStar, MassK0}, massReso, 0);
      itsNClsSoftPi = bach0.itsNClsSoftPi();
      tpcNClsCrossedRowsSoftPi = bach0.tpcNClsCrossedRowsSoftPi();
      tpcChi2NClSoftPi = bach0.tpcChi2NClSoftPi();
    } else if constexpr (Channel == DecayChannel::DplusK0s) {
      sign = bach0.sign();
      massReso = useDeltaMass ? candidate.invMass() + MassDPlus : candidate.invMass();
      cosThetaStar = RecoDecay::cosThetaStar(std::array{bach0.pVector(), bach1.pVector()}, std::array{MassDPlus, MassK0}, massReso, 0);
    } else if constexpr (Channel == DecayChannel::DplusLambda) {
      sign = bach0.sign();
      massReso = useDeltaMass ? candidate.invMass() + MassDPlus : candidate.invMass();
      cosThetaStar = RecoDecay::cosThetaStar(std::array{bach0.pVector(), bach1.pVector()}, std::array{MassDPlus, MassLambda0}, massReso, 0);
    } else if constexpr (Channel == DecayChannel::D0Lambda) {
      massReso = useDeltaMass ? candidate.invMass() + MassD0 : candidate.invMass();
      cosThetaStar = RecoDecay::cosThetaStar(std::array{bach0.pVector(), bach1.pVector()}, std::array{MassD0, MassLambda0}, massReso, 0);
    }
    float y = RecoDecay::y(std::array{candidate.px(), candidate.py(), candidate.pz()}, massReso);
    float eta = RecoDecay::eta(std::array{candidate.px(), candidate.py(), candidate.pz()});
    float phi = RecoDecay::phi(candidate.px(), candidate.py());
    float p = RecoDecay::p(std::array{candidate.px(), candidate.py(), candidate.pz()});
    float e = RecoDecay::e(std::array{candidate.px(), candidate.py(), candidate.pz()}, massReso);

    // MC Rec
    float ptGen{-1.}, invMassGen{-1};
    int8_t origin{0}, flagMcMatchRec{0}, flagCharmBach{0}, flagCharmBachInterm{0}, nKinkedTracks{0};
    int debugMcRec{-1};
    if constexpr (DoMc) {
      ptGen = candidate.ptGen();
      origin = candidate.origin();
      flagMcMatchRec = candidate.flagMcMatchRec();
      debugMcRec = candidate.debugMcRec();
      invMassGen = candidate.invMassGen();
      flagCharmBach = candidate.flagMcMatchRecD();
      flagCharmBachInterm = candidate.flagMcMatchChanD();
      nKinkedTracks = candidate.nTracksDecayed();
      if (fillOnlySignal) {
        if (Channel == DecayChannel::DstarK0s &&
            !hf_decay::hf_cand_reso::particlesToDstarK0s.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flagMcMatchRec)))) {
          return;
        }
        if (Channel == DecayChannel::DplusK0s &&
            !hf_decay::hf_cand_reso::particlesToDplusK0s.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flagMcMatchRec)))) {
          return;
        }
        if (Channel == DecayChannel::DplusLambda &&
            !hf_decay::hf_cand_reso::particlesToDplusLambda.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flagMcMatchRec)))) {
          return;
        }
        if (Channel == DecayChannel::D0Lambda &&
            !hf_decay::hf_cand_reso::particlesToD0Lambda.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flagMcMatchRec)))) {
          return;
        }
      }
      if (origin == RecoDecay::OriginType::Prompt) {
        registry.fill(HIST("hYRecPrompt"), candidate.pt(), y);
      } else if (origin == RecoDecay::OriginType::NonPrompt) {
        registry.fill(HIST("hYRecNonPrompt"), candidate.pt(), y);
      }
    }

    // Ml
    float mlScoreBkg{-1.}, mlScorePrompt{-1.}, mlScoreNonPrompt{-1.};
    if constexpr (WithMl) {
      if constexpr (Channel == DecayChannel::D0Lambda) {
        if (TESTBIT(bach1.v0Type(), BachelorType::Lambda) && !doWrongSign) {
          mlScoreBkg = bach0.mlScoreBkgMassHypo0();
          mlScorePrompt = bach0.mlScorePromptMassHypo0();
          mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo0();
        } else if (TESTBIT(bach1.v0Type(), BachelorType::AntiLambda) && !doWrongSign) {
          mlScoreBkg = bach0.mlScoreBkgMassHypo1();
          mlScorePrompt = bach0.mlScorePromptMassHypo1();
          mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo1();
        } else if (TESTBIT(bach1.v0Type(), BachelorType::Lambda) && doWrongSign) {
          mlScoreBkg = bach0.mlScoreBkgMassHypo1();
          mlScorePrompt = bach0.mlScorePromptMassHypo1();
          mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo1();
        } else if (TESTBIT(bach1.v0Type(), BachelorType::AntiLambda) && doWrongSign) {
          mlScoreBkg = bach0.mlScoreBkgMassHypo0();
          mlScorePrompt = bach0.mlScorePromptMassHypo0();
          mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo0();
        }
      } else {
        mlScoreBkg = bach0.mlScoreBkgMassHypo0();
        mlScorePrompt = bach0.mlScorePromptMassHypo0();
        mlScoreNonPrompt = bach0.mlScoreNonpromptMassHypo0();
      }
    }
    // Collision properties
    registry.fill(HIST("hNPvCont"), collision.numContrib());
    registry.fill(HIST("hZvert"), collision.posZ());
    registry.fill(HIST("hBz"), collision.bz());
    // Candidate properties
    registry.fill(HIST("hMass"), candidate.invMass());
    registry.fill(HIST("hMassProng0"), candidate.invMassProng0());
    registry.fill(HIST("hMassProng1"), candidate.invMassProng1());
    registry.fill(HIST("hPt"), candidate.pt());
    registry.fill(HIST("hPtProng0"), candidate.ptProng0());
    registry.fill(HIST("hPtProng1"), candidate.ptProng1());
    if (fillSparses) {
      registry.fill(HIST("hSparse"), candidate.pt(), candidate.ptProng0(), candidate.ptProng1(), candidate.invMass(), candidate.invMassProng0(), candidate.invMassProng1(), cosThetaStar, mlScoreBkg, mlScoreNonPrompt);
    }
    if (fillTrees) {
      hfCandResoLite(
        candidate.invMass(),
        candidate.pt(),
        p,
        y,
        eta,
        phi,
        e,
        cosThetaStar,
        sign,
        // Bachelors Properties
        candidate.invMassProng0(),
        bach0.pt(),
        mlScoreBkg,
        mlScorePrompt,
        mlScoreNonPrompt,
        bach0.itsNClsProngMin(),
        bach0.tpcNClsCrossedRowsProngMin(),
        bach0.tpcChi2NClProngMax(),
        itsNClsSoftPi,
        tpcNClsCrossedRowsSoftPi,
        tpcChi2NClSoftPi,
        candidate.invMassProng1(),
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
        invMassGen,
        flagCharmBach,
        flagCharmBachInterm,
        nKinkedTracks);
    }
  } // fillCand

  // Process data
  /// \tparam channel is the decay channel of the Resonance
  /// \param Coll is the reduced collisions table
  /// \param CharmBach is the reduced 3 prong table
  /// \param V0Bach is the reduced v0 table
  /// \param Cand is the candidates table
  template <bool DoMc, bool WithMl, DecayChannel Channel, typename Coll, typename Candidates, typename CharmBach>
  void processData(Coll const&, Candidates const& candidates, CharmBach const&, aod::HfRedVzeros const&)
  {
    for (const auto& cand : candidates) {
      if (ptMinReso >= 0 && cand.pt() < ptMinReso) {
        continue;
      }
      if ((massResoMin >= 0 && cand.invMass() < massResoMin) ||
          (massResoMax >= 0 && cand.invMass() > massResoMax)) {
        continue;
      }
      if (doWrongSign && cand.isWrongSign() == 0) {
        continue;
      }
      if (!doWrongSign && cand.isWrongSign() != 0) {
        continue;
      }

      float massReso{0};
      if (useDeltaMass) {
        switch (Channel) {
          case DecayChannel::DstarK0s:
            massReso = cand.invMass() + MassDStar;
            break;
          case DecayChannel::DplusK0s:
            massReso = cand.invMass() + MassDPlus;
            break;
          case DecayChannel::DplusLambda:
            massReso = cand.invMass() + MassDPlus;
            break;
          case DecayChannel::D0Lambda:
            massReso = cand.invMass() + MassD0;
            break;
          default:
            break;
        }
      } else {
        massReso = cand.invMass();
      }
      if (yCandRecoMax >= 0. && std::abs(RecoDecay::y(std::array{cand.px(), cand.py(), cand.pz()}, massReso)) > yCandRecoMax) {
        continue;
      }
      auto coll = cand.template hfRedCollision_as<Coll>();
      auto bach0 = cand.template prong0_as<CharmBach>();
      auto bach1 = cand.template prong1_as<aod::HfRedVzeros>();
      fillCand<DoMc, WithMl, Channel>(cand, coll, bach0, bach1);
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
  template <DecayChannel Channel>
  void fillCandMcGen(aod::HfMcGenRedResos const& mcParticles)
  {
    for (const auto& particle : mcParticles) {
      auto ptParticle = particle.ptTrack();
      auto yParticle = particle.yTrack();
      auto originParticle = particle.origin();
      auto flag = particle.flagMcMatchGen();
      std::array<float, 2> ptProngs = {particle.ptProng0(), particle.ptProng1()};
      std::array<float, 2> etaProngs = {particle.etaProng0(), particle.etaProng1()};
      bool const prongsInAcc = isProngInAcceptance(etaProngs[0], ptProngs[0]) && isProngInAcceptance(etaProngs[1], ptProngs[1]);
      if (Channel == DecayChannel::DstarK0s &&
          !hf_decay::hf_cand_reso::particlesToDstarK0s.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flag)))) {
        continue;
      }
      if (Channel == DecayChannel::DplusK0s &&
          !hf_decay::hf_cand_reso::particlesToDplusK0s.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flag)))) {
        continue;
      }
      if (Channel == DecayChannel::DplusLambda &&
          !hf_decay::hf_cand_reso::particlesToDplusLambda.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flag)))) {
        continue;
      }
      if (Channel == DecayChannel::D0Lambda &&
          !hf_decay::hf_cand_reso::particlesToD0Lambda.contains(static_cast<hf_decay::hf_cand_reso::DecayChannelMain>(std::abs(flag)))) {
        continue;
      }
      registry.fill(HIST("hYGenAll"), ptParticle, yParticle);
      if (yCandGenMax >= 0. && std::abs(yParticle) > yCandGenMax) {
        continue;
      }
      if (originParticle == RecoDecay::OriginType::Prompt) { // prompt particles
        registry.fill(HIST("hYGenPrompt"), ptParticle, yParticle);
        if (prongsInAcc) {
          registry.fill(HIST("hYGenPromptWithProngsInAcceptance"), ptParticle, yParticle);
        }
      } else if (originParticle == RecoDecay::OriginType::NonPrompt) {
        registry.fill(HIST("hYGenNonPrompt"), ptParticle, yParticle);
        if (prongsInAcc) {
          registry.fill(HIST("hYGenNonPromptWithProngsInAcceptance"), ptParticle, yParticle);
        }
      }
      if (fillSparses) {
        registry.fill(HIST("hPtYGenSig"), ptParticle, yParticle, originParticle, flag);
        if (prongsInAcc) {
          registry.fill(HIST("hPtYWithProngsInAccepanceGenSig"), ptParticle, yParticle, originParticle, flag);
        }
      }
      if (fillTrees) {
        hfGenResoLite(ptParticle, yParticle, originParticle, flag);
      }
    }
  } // fillCandMcGen

  // process functions
  void processDstarK0sData(aod::HfRedCollisions const& collisions,
                           ReducedResoDstarV0 const& candidates,
                           aod::HfRedDstarNoTrks const& charmBachs,
                           aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::DstarK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDstarK0sData, "Process data for DstarK0s analysis", true);

  void processDplusK0sData(aod::HfRedCollisions const& collisions,
                           ReducedReso3PrV0 const& candidates,
                           aod::HfRed3PrNoTrks const& charmBachs,
                           aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::DplusK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusK0sData, "Process data for DplusK0s analysis", false);

  void processDplusLambdaData(aod::HfRedCollisions const& collisions,
                              ReducedReso3PrV0 const& candidates,
                              aod::HfRed3PrNoTrks const& charmBachs,
                              aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::DplusLambda>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusLambdaData, "Process data for DplusLambda analysis", false);

  void processD0LambdaData(aod::HfRedCollisions const& collisions,
                           ReducedReso2PrV0 const& candidates,
                           aod::HfRed2PrNoTrks const& charmBachs,
                           aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, false, DecayChannel::D0Lambda>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processD0LambdaData, "Process data for D0Lambda analysis", false);

  // Process data with ML
  void processDstarK0sDataWithMl(aod::HfRedCollisions const& collisions,
                                 ReducedResoDstarV0 const& candidates,
                                 soa::Join<aod::HfRedDstarNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                                 aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::DstarK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDstarK0sDataWithMl, "Process data for DstarK0s analysis with Ml", false);

  void processDplusK0sDataWithMl(aod::HfRedCollisions const& collisions,
                                 ReducedReso3PrV0 const& candidates,
                                 soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                                 aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::DplusK0s>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusK0sDataWithMl, "Process data for DplusK0s analysis with Ml", false);

  void processDplusLambdaDataWithMl(aod::HfRedCollisions const& collisions,
                                    ReducedReso3PrV0 const& candidates,
                                    soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                                    aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::DplusLambda>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusLambdaDataWithMl, "Process data for DplusLambda analysis with Ml", false);

  void processD0LambdaDataWithMl(aod::HfRedCollisions const& collisions,
                                 ReducedReso2PrV0 const& candidates,
                                 soa::Join<aod::HfRed2PrNoTrks, aod::HfRed2ProngsMl> const& charmBachs,
                                 aod::HfRedVzeros const& v0Bachs)
  {
    processData<false, true, DecayChannel::D0Lambda>(collisions, candidates, charmBachs, v0Bachs);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processD0LambdaDataWithMl, "Process data for D0Lambda analysis with Ml", false);

  // MC
  void processDstarK0sMC(aod::HfRedCollisions const& collisions,
                         ReducedResoDstarV0MC const& candidates,
                         aod::HfRedDstarNoTrks const& charmBachs,
                         aod::HfRedVzeros const& v0Bachs,
                         aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, false, DecayChannel::DstarK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DstarK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDstarK0sMC, "Process MC for DstarK0s analysis", false);

  void processDplusK0sMC(aod::HfRedCollisions const& collisions,
                         ReducedReso3PrV0MC const& candidates,
                         aod::HfRed3PrNoTrks const& charmBachs,
                         aod::HfRedVzeros const& v0Bachs,
                         aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, false, DecayChannel::DplusK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DplusK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusK0sMC, "Process MC for DplusK0s analysis", false);

  void processDplusLambdaMC(aod::HfRedCollisions const& collisions,
                            ReducedReso3PrV0MC const& candidates,
                            aod::HfRed3PrNoTrks const& charmBachs,
                            aod::HfRedVzeros const& v0Bachs,
                            aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, false, DecayChannel::DplusLambda>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DplusLambda>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusLambdaMC, "Process MC for DplusLambda analysis", false);

  void processD0LambdaMC(aod::HfRedCollisions const& collisions,
                         ReducedReso2PrV0MC const& candidates,
                         aod::HfRed2PrNoTrks const& charmBachs,
                         aod::HfRedVzeros const& v0Bachs,
                         aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, false, DecayChannel::D0Lambda>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::D0Lambda>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processD0LambdaMC, "Process MC for D0Lambda analysis", false);

  // MC with Ml
  void processDstarK0sMCWithMl(aod::HfRedCollisions const& collisions,
                               ReducedResoDstarV0MC const& candidates,
                               soa::Join<aod::HfRedDstarNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                               aod::HfRedVzeros const& v0Bachs,
                               aod::HfMcGenRedResos const& charmBachsMc)
  {
    processData<true, true, DecayChannel::DstarK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DstarK0s>(charmBachsMc);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDstarK0sMCWithMl, "Process MC for DstarK0s analysis with Ml", false);

  void processDplusK0sMCWithMl(aod::HfRedCollisions const& collisions,
                               ReducedReso3PrV0MC const& candidates,
                               soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                               aod::HfRedVzeros const& v0Bachs,
                               aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, true, DecayChannel::DplusK0s>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DplusK0s>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusK0sMCWithMl, "Process MC for DplusK0s analysis with Ml", false);

  void processDplusLambdaMCWithMl(aod::HfRedCollisions const& collisions,
                                  ReducedReso3PrV0MC const& candidates,
                                  soa::Join<aod::HfRed3PrNoTrks, aod::HfRed3ProngsMl> const& charmBachs,
                                  aod::HfRedVzeros const& v0Bachs,
                                  aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, true, DecayChannel::DplusLambda>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::DplusLambda>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processDplusLambdaMCWithMl, "Process MC for DplusLambda analysis with Ml", false);

  void processD0LambdaMCWithMl(aod::HfRedCollisions const& collisions,
                               ReducedReso2PrV0MC const& candidates,
                               soa::Join<aod::HfRed2PrNoTrks, aod::HfRed2ProngsMl> const& charmBachs,
                               aod::HfRedVzeros const& v0Bachs,
                               aod::HfMcGenRedResos const& mcParticles)
  {
    processData<true, true, DecayChannel::D0Lambda>(collisions, candidates, charmBachs, v0Bachs);
    fillCandMcGen<DecayChannel::D0Lambda>(mcParticles);
  }
  PROCESS_SWITCH(HfTaskCharmResoToDV0Reduced, processD0LambdaMCWithMl, "Process MC for D0Lambda analysis with Ml", false);

}; // struct HfTaskCharmResoToDV0Reduced
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskCharmResoToDV0Reduced>(cfgc)};
}
