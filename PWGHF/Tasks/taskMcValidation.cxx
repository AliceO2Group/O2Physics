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

/// \file taskMcValidation.cxx
/// \brief Monte Carlo validation task
/// \note generated and reconstructed level validation
///
/// \author Antonio Palasciano <antonio.palasciano@cern.ch>, Università degli Studi di Bari & INFN, Sezione di Bari
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
/// \author Fabio Catalano <fabio.catalano@cern.ch>, CERN

#include "PWGHF/Core/CentralityEstimation.h"
#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"
#include "PWGHF/Utils/utilsEvSelHf.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DeviceSpec.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TPDGCode.h>
#include <TString.h>

#include <sys/types.h>

#include <Rtypes.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <numeric>
#include <string_view>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::hf_evsel;
using namespace o2::hf_centrality;
using namespace o2::hf_occupancy;

namespace
{
enum DecayChannels { DzeroToKPi = 0,
                     DstarToDzeroPi,
                     DplusToPiKPi,
                     DplusToPhiPiToKKPi, // resonant channel with Phi, Dplus -> PhiPi -> KKPi
                     DsToPhiPiToKKPi,    // resonant channel with Phi, Ds -> PhiPi -> KKPi
                     DsToK0starKToKKPi,  // resonant channel with K0*, Ds -> K0*K -> KKPi
                     Ds1ToDStarK0s,
                     Ds2StarToDPlusK0s,
                     D10ToDStarPi,
                     D2Star0ToDPlusPi,
                     B0ToDminusPi,
                     BplusToD0Pi,
                     BsToDsPi,
                     LcToPKPi,
                     LcToPiK0s,
                     XiCplusToPKPi,
                     XiCplusToXiPiPi,
                     XiCzeroToXiPi,
                     OmegaCToOmegaPi,
                     OmegaCToXiPi,
                     NChannels
}; // always keep nChannels at the end

constexpr int NCharmMesonChannels = 10;                                                 // number of charm meson channels
constexpr int NBeautyChannels = 3;                                                      // number of beauty hadron channels
constexpr int NCharmBaryonChannels = NChannels - NCharmMesonChannels - NBeautyChannels; // number of charm baryon channels
constexpr int NOriginTypes = 2;                                                         // number of origin types (prompt, non-prompt; only for charm hadrons)
constexpr std::array<int, NChannels> PDGArrayParticle = {o2::constants::physics::Pdg::kD0, o2::constants::physics::Pdg::kDStar,
                                                         o2::constants::physics::Pdg::kDPlus, o2::constants::physics::Pdg::kDPlus,
                                                         o2::constants::physics::Pdg::kDS, o2::constants::physics::Pdg::kDS, o2::constants::physics::Pdg::kDS1, o2::constants::physics::Pdg::kDS2Star,
                                                         o2::constants::physics::Pdg::kD10, o2::constants::physics::Pdg::kD2Star0,
                                                         o2::constants::physics::Pdg::kB0, o2::constants::physics::Pdg::kBPlus, o2::constants::physics::Pdg::kBS,
                                                         o2::constants::physics::Pdg::kLambdaCPlus, o2::constants::physics::Pdg::kLambdaCPlus,
                                                         o2::constants::physics::Pdg::kXiCPlus, o2::constants::physics::Pdg::kXiCPlus, o2::constants::physics::Pdg::kXiC0,
                                                         o2::constants::physics::Pdg::kOmegaC0, o2::constants::physics::Pdg::kOmegaC0};
constexpr std::array<unsigned int, NChannels> NDaughters = {2, 3, 3, 3, 3, 3, 5, 5, 4, 4, 4, 3, 4, 3, 3, 3, 5, 4, 4, 4};
constexpr std::array<int, NChannels> MaxDepthForSearch = {1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 2, 3, 2, 3, 2, 4, 3, 3, 3};
// keep coherent indexing with PDGArrayParticle
// FIXME: look for a better solution
constexpr std::array<std::array<int, 2>, NChannels> ArrPdgFinal2Prong = {{{+kPiPlus, -kKPlus}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}};
constexpr std::array<std::array<int, 3>, NChannels> ArrPdgFinal3Prong = {{{}, {+kPiPlus, -kKPlus, +kPiPlus}, {+kPiPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {+kKPlus, -kKPlus, +kPiPlus}, {}, {}, {}, {}, {}, {-kPiPlus, +kKPlus, +kPiPlus}, {}, {+kProton, -kKPlus, +kPiPlus}, {+kProton, -kPiPlus, +kPiPlus}, {+kProton, -kKPlus, +kPiPlus}, {}, {}, {}, {}}};
constexpr std::array<std::array<int, 4>, NChannels> ArrPdgFinal4Prong = {{{}, {}, {}, {}, {}, {}, {}, {}, {+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus}, {+kPiPlus, -kKPlus, +kPiPlus, -kPiPlus}, {-kPiPlus, +kKPlus, -kPiPlus, +kPiPlus}, {}, {-kKPlus, +kKPlus, -kPiPlus, +kPiPlus}, {}, {}, {}, {}, {+kPiPlus, -kPiPlus, -kPiPlus, +kProton}, {+kPiPlus, -kKPlus, -kPiPlus, +kProton}, {+kPiPlus, -kPiPlus, -kPiPlus, +kProton}}};
constexpr std::array<std::array<int, 5>, NChannels> ArrPdgFinal5Prong = {{{}, {}, {}, {}, {}, {}, {+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, {+kPiPlus, -kKPlus, +kPiPlus, +kPiPlus, -kPiPlus}, {}, {}, {}, {}, {}, {}, {}, {}, {+kPiPlus, +kPiPlus, -kPiPlus, -kPiPlus, +kProton}, {}, {}, {}}};
constexpr std::string_view Labels[NChannels] = {"D^{0} #rightarrow K#pi", "D*^{+} #rightarrow D^{0}#pi", "D^{+} #rightarrow K#pi#pi", "D^{+} #rightarrow KK#pi", "D_{s}^{+} #rightarrow #Phi#pi #rightarrow KK#pi",
                                                "D_{s}^{+} #rightarrow #bar{K}^{*0}K #rightarrow KK#pi", "D_{s}1 #rightarrow D*^{+}K^{0}_{s}", "D_{s}2* #rightarrow D^{+}K^{0}_{s}", "D1^{0} #rightarrow D*^{+}#pi",
                                                "D2^{*} #rightarrow D^{+}#pi",
                                                "B^{0} #rightarrow D^{-}#pi", "B^{+} #rightarrow D^{0}#pi", "B_{s}^{+} #rightarrow D_{s}^{-}#pi",
                                                "#Lambda_{c}^{+} #rightarrow pK#pi",
                                                "#Lambda_{c}^{+} #rightarrow pK^{0}_{s}", "#Xi_{c}^{+} #rightarrow pK#pi",
                                                "#Xi_{c}^{+} #rightarrow #Xi#pi#pi", "#Xi_{c}^{0} #rightarrow #Xi#pi", "#Omega_{c}^{0} #rightarrow #Omega#pi", "#Omega_{c}^{0} #rightarrow #Xi#pi"};
constexpr std::string_view ParticleNames[NChannels] = {"DzeroToKPi", "DstarToDzeroPi", "DplusToPiKPi", "DplusToPhiPiToKKPi", "DsToPhiPiToKKPi", "DsToK0starKToKKPi", "Ds1ToDStarK0s", "Ds2StarToDPlusK0s", "D10ToDStarPi", "D2Star0ToDPlusPi", "B0ToDminusPi", "BplusToD0Pi", "BsToDsPi",
                                                       "LcToPKPi", "LcToPiK0s", "XiCplusToPKPi", "XiCplusToXiPiPi", "XiCzeroToXiPi", "OmegaCToOmegaPi", "OmegaCToXiPi"};
constexpr std::string_view OriginNames[NOriginTypes] = {"Prompt", "NonPrompt"};
} // namespace

/// Generated Level Validation
///
/// - Number of HF quarks produced per collision
/// - Number of candidates per collision
/// - Momentum Conservation for these particles
struct HfTaskMcValidationGen {

  using BCsInfo = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels>;
  using CollisionsNoCents = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CollisionsFT0Cs = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Cs>;
  using CollisionsFT0Ms = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, aod::CentFT0Ms>;
  using McCollisionsCentFT0Ms = soa::Join<aod::McCollisions, aod::McCentFT0Ms>;
  PresliceUnsorted<CollisionsNoCents> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsFT0Cs> colPerMcCollisionFT0C = aod::mccollisionlabel::mcCollisionId;
  PresliceUnsorted<CollisionsFT0Ms> colPerMcCollisionFT0M = aod::mccollisionlabel::mcCollisionId;
  Preslice<aod::McParticles> mcParticlesPerMcCollision = aod::mcparticle::mcCollisionId;

  Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};
  Configurable<bool> rejectParticlesFromBkgEvent{"rejectParticlesFromBkgEvent", true, "Reject particles"};
  Configurable<bool> storeOccupancy{"storeOccupancy", false, "Store collision occupancy for dedicated studies"};

  HfEventSelectionMc hfEvSelMc; // mc event selection and monitoring

  AxisSpec axisNhadrons{NChannels, -0.5, static_cast<float>(NChannels) - 0.5};
  AxisSpec axisNquarks{20, -0.5, 19.5};
  AxisSpec axisResiduals{100, -0.01, 0.01};
  AxisSpec axisPt{100, 0., 50.};
  AxisSpec axisY{100, -5., 5.};
  AxisSpec axisCent{110, 0., 110.};
  AxisSpec axisOcc{3000, 0., 15000.};
  AxisSpec axisCharmMesonSpecies{NCharmMesonChannels, -0.5, static_cast<float>(NCharmMesonChannels) - 0.5};
  AxisSpec axisBeautySpecies{NBeautyChannels, -0.5, static_cast<float>(NBeautyChannels) - 0.5};
  AxisSpec axisCharmBaryonSpecies{NCharmBaryonChannels, -0.5, static_cast<float>(NCharmBaryonChannels) - 0.5};
  AxisSpec axisDecLen{100, 0., 10000.};

  HistogramRegistry registry{
    "registry",
    {{"hNevGen", "Generated events counter; Gen. events; entries", {HistType::kTH1F, {{1, -0.5, +0.5}}}},
     {"hMomentumCheck", "Mom. Conservation (1 = true, 0 = false) (#it{#epsilon} = 1 MeV/#it{c}); Mom. Conservation result; entries", {HistType::kTH1F, {{2, -0.5, +1.5}}}},
     {"hPtDiffMotherDaughterGen", "Pt Difference Mother-Daughters; #Delta#it{p}_{T}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPxDiffMotherDaughterGen", "Px Difference Mother-Daughters; #Delta#it{p}_{x}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPyDiffMotherDaughterGen", "Py Difference Mother-Daughters; #Delta#it{p}_{y}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPzDiffMotherDaughterGen", "Pz Difference Mother-Daughters; #Delta#it{p}_{z}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"hPDiffMotherDaughterGen", "P  Difference Mother-Daughters; #Delta#it{p}^{gen} (GeV/#it{c}); entries", {HistType::kTH1F, {axisResiduals}}},
     {"Quarks/hCountC", "Event counter - Number of charm quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountB", "Event counter - Number of beauty quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountCbar", "Event counter - Number of anti-charm quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hCountBbar", "Event counter - Number of anti-beauty quarks; Events Per Collision; entries", {HistType::kTH1F, {axisNquarks}}},
     {"Quarks/hPtVsYCharmQuark", "Y vs. Pt - charm quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {axisPt, axisY}}},
     {"Quarks/hPtVsYBeautyQuark", "Y vs. Pt - beauty quarks ; #it{p}_{T}^{gen} (GeV/#it{c}); #it{y}^{gen}", {HistType::kTH2F, {axisPt, axisY}}},
     {"PromptCharmMesons/hPromptMesonsPtDistr", "Pt distribution vs prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisCharmMesonSpecies, axisPt}}},
     {"PromptCharmMesons/hPromptMesonsPtCentDistr", "Pt vs Cent distribution vs prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Centrality (%)", {HistType::kTH3F, {axisCharmMesonSpecies, axisPt, axisCent}}},
     {"PromptCharmMesons/hPromptMesonsPtOccDistr", "Pt vs Occ distribution vs prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Occupancy", {HistType::kTH3F, {axisCharmMesonSpecies, axisPt, axisOcc}}},
     {"PromptCharmMesons/hPromptMesonsYDistr", "Y distribution vs prompt charm meson; ; #it{y}^{gen}", {HistType::kTH2F, {axisCharmMesonSpecies, axisY}}},
     {"PromptCharmMesons/hPromptMesonsDecLenDistr", "Decay length distribution vs prompt charm meson; ; decay length (#mum)", {HistType::kTH2F, {axisCharmMesonSpecies, axisDecLen}}},
     {"NonPromptCharmMesons/hNonPromptMesonsPtDistr", "Pt distribution vs non-prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisCharmMesonSpecies, axisPt}}},
     {"NonPromptCharmMesons/hNonPromptMesonsPtCentDistr", "Pt vs Cent distribution vs non-prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Centrality (%)", {HistType::kTH3F, {axisCharmMesonSpecies, axisPt, axisCent}}},
     {"NonPromptCharmMesons/hNonPromptMesonsPtOccDistr", "Pt vs Occ distribution vs non-prompt charm meson in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Occupancy", {HistType::kTH3F, {axisCharmMesonSpecies, axisPt, axisOcc}}},
     {"NonPromptCharmMesons/hNonPromptMesonsYDistr", "Y distribution vs non-prompt charm meson; ; #it{y}^{gen}", {HistType::kTH2F, {axisCharmMesonSpecies, axisY}}},
     {"NonPromptCharmMesons/hNonPromptMesonsDecLenDistr", "Decay length distribution vs non-prompt charm meson; ; decay length (#mum)", {HistType::kTH2F, {axisCharmMesonSpecies, axisDecLen}}},
     {"Beauty/hPtDistr", "Pt distribution vs beauty hadron in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisBeautySpecies, axisPt}}},
     {"Beauty/hPtCentDistr", "Pt vs Cent distribution vs beauty hadron in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Centrality (%)", {HistType::kTH3F, {axisBeautySpecies, axisPt, axisCent}}},
     {"Beauty/hPtOccDistr", "Pt vs Occ distribution vs beauty hadron in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Occupancy", {HistType::kTH3F, {axisBeautySpecies, axisPt, axisOcc}}},
     {"Beauty/hYDistr", "Y distribution vs beauty hadron; ; #it{y}^{gen}", {HistType::kTH2F, {axisBeautySpecies, axisY}}},
     {"Beauty/hDecLenDistr", "Decay length distribution vs beauty hadron; ; decay length (#mum)", {HistType::kTH2F, {axisBeautySpecies, axisDecLen}}},
     {"PromptCharmBaryons/hPromptBaryonsPtDistr", "Pt distribution vs prompt charm baryon in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisCharmBaryonSpecies, axisPt}}},
     {"PromptCharmBaryons/hPromptBaryonsPtCentDistr", "Pt vs Cent distribution vs prompt charm baryons in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Centrality (%)", {HistType::kTH3F, {axisCharmBaryonSpecies, axisPt, axisCent}}},
     {"PromptCharmBaryons/hPromptBaryonsYDistr", "Y distribution vs prompt charm baryon; ; #it{y}^{gen}", {HistType::kTH2F, {axisCharmBaryonSpecies, axisY}}},
     {"PromptCharmBaryons/hPromptBaryonsPtOccDistr", "Pt vs Occ distribution vs prompt charm baryons in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Occupancy", {HistType::kTH3F, {axisCharmBaryonSpecies, axisPt, axisOcc}}},
     {"PromptCharmBaryons/hPromptBaryonsDecLenDistr", "Decay length distribution vs prompt charm baryon; ; decay length (#mum)", {HistType::kTH2F, {axisCharmBaryonSpecies, axisDecLen}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsPtDistr", "Pt distribution vs non-prompt charm baryon in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c})", {HistType::kTH2F, {axisCharmBaryonSpecies, axisPt}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsPtCentDistr", "Pt vs Cent distribution vs prompt charm baryons in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Centrality (%)", {HistType::kTH3F, {axisCharmBaryonSpecies, axisPt, axisCent}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsPtOccDistr", "Pt vs Occ distribution vs prompt charm baryons in |#it{y}^{gen}|<0.5; ; #it{p}_{T}^{gen} (GeV/#it{c}); Occupancy", {HistType::kTH3F, {axisCharmBaryonSpecies, axisPt, axisOcc}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsYDistr", "Y distribution vs non-prompt charm baryon; ; #it{y}^{gen}", {HistType::kTH2F, {axisCharmBaryonSpecies, axisY}}},
     {"NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr", "Decay length distribution vs non-prompt charm baryon; ; decay length (#mum)", {HistType::kTH2F, {axisCharmBaryonSpecies, axisDecLen}}}}};

  void init(InitContext& initContext)
  {

    std::array<bool, 3> processes = {doprocessNoCentSel, doprocessCentFT0C, doprocessCentFT0M};
    if (std::accumulate(processes.begin(), processes.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for generated particles can be enabled at a time.");
    }

    // add per species histograms
    for (auto originName : OriginNames) {
      for (int iChannel = 0; iChannel < NCharmMesonChannels; iChannel++) { // Charm mesons
        registry.add(Form("%sCharmMesons/hCount%s%s", originName.data(), originName.data(), ParticleNames[iChannel].data()),
                     Form("Event counter - %s %s; Events Per Collision; entries", originName.data(), Labels[iChannel].data()), {HistType::kTH1F, {axisNhadrons}});
      }
      for (int iChannel = NCharmMesonChannels + NBeautyChannels; iChannel < NChannels; iChannel++) { // Charm baryons
        registry.add(Form("%sCharmBaryons/hCount%s%s", originName.data(), originName.data(), ParticleNames[iChannel].data()),
                     Form("Event counter - %s %s; Events Per Collision; entries", originName.data(), Labels[iChannel].data()), {HistType::kTH1F, {axisNhadrons}});
      }
    }
    for (int iChannel = NCharmMesonChannels; iChannel < NCharmMesonChannels + NBeautyChannels; iChannel++) { // Beauty mesons
      registry.add(Form("Beauty/hCount%s", ParticleNames[iChannel].data()),
                   Form("Event counter - %s; Events Per Collision; entries", Labels[iChannel].data()), {HistType::kTH1F, {axisNhadrons}});
    }

    for (auto iBin = 1; iBin <= NCharmMesonChannels; ++iBin) {
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH3>(HIST("PromptCharmMesons/hPromptMesonsPtCentDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH3>(HIST("PromptCharmMesons/hPromptMesonsPtOccDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsYDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH2>(HIST("PromptCharmMesons/hPromptMesonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH3>(HIST("NonPromptCharmMesons/hNonPromptMesonsPtCentDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH3>(HIST("NonPromptCharmMesons/hNonPromptMesonsPtOccDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsYDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmMesons/hNonPromptMesonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin - 1].data());
    }
    for (auto iBin = 1; iBin <= NBeautyChannels; ++iBin) {
      registry.get<TH2>(HIST("Beauty/hPtDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels - 1].data());
      registry.get<TH3>(HIST("Beauty/hPtCentDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels - 1].data());
      registry.get<TH3>(HIST("Beauty/hPtOccDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels - 1].data());
      registry.get<TH2>(HIST("Beauty/hYDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels - 1].data());
      registry.get<TH2>(HIST("Beauty/hDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels - 1].data());
    }
    for (auto iBin = 1; iBin <= NCharmBaryonChannels; ++iBin) {
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH3>(HIST("PromptCharmBaryons/hPromptBaryonsPtCentDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH3>(HIST("PromptCharmBaryons/hPromptBaryonsPtOccDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsYDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH2>(HIST("PromptCharmBaryons/hPromptBaryonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH3>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtCentDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH3>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtOccDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsYDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
      registry.get<TH2>(HIST("NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr"))->GetXaxis()->SetBinLabel(iBin, Labels[iBin + NCharmMesonChannels + NBeautyChannels - 1].data());
    }

    // inspect for which particle species the candidates were created and which zPvPosMax cut was set for reconstructed
    const auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (const DeviceSpec& device : workflows.devices) {
      if (device.name == "hf-task-mc-validation-rec") {
        hfEvSelMc.configureFromDevice(device);
        break;
      }
    }

    hfEvSelMc.addHistograms(registry); // particles monitoring
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename GenColl, typename Particles, typename RecoColls>
  void runCheckGenParticles(GenColl const& mcCollision, Particles const& mcParticles, RecoColls const& recoCollisions, BCsInfo const&, std::array<int, NChannels>& counterPrompt, std::array<int, NChannels>& counterNonPrompt)
  {
    if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
      return;
    }
    registry.fill(HIST("hNevGen"), 1);

    // Slice the collisions table to get the collision info for the current MC collision
    float centrality{105.f};
    float occupancy{0.f};
    if (storeOccupancy) {
      occupancy = o2::hf_occupancy::getOccupancyGenColl(recoCollisions, OccupancyEstimator::Its);
    }
    o2::hf_evsel::HfCollisionRejectionMask rejectionMask{};
    if constexpr (CentEstimator == CentralityEstimator::FT0C || CentEstimator == CentralityEstimator::FT0M || CentEstimator == CentralityEstimator::None) {
      rejectionMask = hfEvSelMc.getHfMcCollisionRejectionMask<BCsInfo, CentEstimator>(mcCollision, recoCollisions, centrality);
    }
    hfEvSelMc.fillHistograms<CentEstimator>(mcCollision, rejectionMask);
    if (rejectionMask != 0) {
      return;
    }

    int cPerCollision = 0;
    int cBarPerCollision = 0;
    int bPerCollision = 0;
    int bBarPerCollision = 0;

    for (const auto& particle : mcParticles) {

      if (rejectParticlesFromBkgEvent && particle.fromBackgroundEvent()) {
        continue;
      }

      if (!particle.has_mothers()) {
        continue;
      }

      int const particlePdgCode = particle.pdgCode();
      bool isDiffFromMothers = true;
      for (const auto& mother : particle.template mothers_as<Particles>()) {
        if (particlePdgCode == mother.pdgCode()) {
          isDiffFromMothers = false;
          break;
        }
      }
      if (isDiffFromMothers) {
        switch (particlePdgCode) {
          case kCharm:
            cPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kCharmBar:
            cBarPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYCharmQuark"), particle.pt(), particle.y());
            break;
          case kBottom:
            bPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
          case kBottomBar:
            bBarPerCollision++;
            registry.fill(HIST("Quarks/hPtVsYBeautyQuark"), particle.pt(), particle.y());
            break;
        }
      }

      // Checking the decay of the particles and the momentum conservation
      for (std::size_t iD = 0; iD < PDGArrayParticle.size(); ++iD) {
        if (std::abs(particlePdgCode) != PDGArrayParticle[iD]) {
          continue;
        }

        std::vector<int> listDaughters{};

        // Check that the decay channel is correct and retrieve the daughters
        if (NDaughters[iD] == 2) {
          if (!RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], ArrPdgFinal2Prong[iD], true, nullptr, MaxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
        }

        if (NDaughters[iD] == 3) {
          if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], ArrPdgFinal3Prong[iD], true, nullptr, MaxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
          if (iD == DstarToDzeroPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kD0, +kPiPlus}, true)) {
            continue;
          }
          if ((iD == DplusToPhiPiToKKPi || iD == DsToPhiPiToKKPi) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kPhi, +kPiPlus}, false) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kPhi, -kPiPlus}, false)) {
            continue;
          }
          if (iD == DsToK0starKToKKPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{-313, +kKPlus}, true)) {
            continue;
          }
          if (iD == LcToPiK0s &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kK0Short, +kProton}, false, nullptr, 2) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+kK0Short, -kProton}, false, nullptr, 2)) {
            continue;
          }
          if (iD == BplusToD0Pi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{-o2::constants::physics::Pdg::kD0, +kPiPlus}, true)) {
            continue;
          }
        }

        if (NDaughters[iD] == 4) {
          if (iD != B0ToDminusPi) {
            if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], ArrPdgFinal4Prong[iD], true, nullptr, MaxDepthForSearch[iD], &listDaughters)) {
              continue;
            }
          } else { // For B0 we consider flavour oscillations
            if (!RecoDecay::isMatchedMCGen<true, false>(mcParticles, particle, PDGArrayParticle[iD], ArrPdgFinal4Prong[iD], true, nullptr, MaxDepthForSearch[iD], &listDaughters)) {
              continue;
            }
          }
          if (iD == D10ToDStarPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kDStar, -kPiPlus}, true)) {
            continue;
          }
          if (iD == D2Star0ToDPlusPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+o2::constants::physics::Pdg::kDPlus, -kPiPlus}, true)) {
            continue;
          }
          if ((iD == XiCzeroToXiPi || iD == OmegaCToXiPi) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kXiMinus, +kPiPlus}, true)) {
            continue;
          }
          if (iD == OmegaCToOmegaPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kOmegaMinus, +kPiPlus}, true)) {
            continue;
          }
          if (iD == B0ToDminusPi &&
              !RecoDecay::isMatchedMCGen<true, false>(mcParticles, particle, PDGArrayParticle[iD], std::array{-o2::constants::physics::Pdg::kDPlus, +kPiPlus}, true)) {
            continue;
          }
          if (iD == BsToDsPi &&
              !RecoDecay::isMatchedMCGen<true, false>(mcParticles, particle, PDGArrayParticle[iD], std::array{-o2::constants::physics::Pdg::kDS, +kPiPlus}, true)) {
            continue;
          }
        }

        if (NDaughters[iD] == 5) {
          if (!RecoDecay::isMatchedMCGen<false, true>(mcParticles, particle, PDGArrayParticle[iD], ArrPdgFinal5Prong[iD], true, nullptr, MaxDepthForSearch[iD], &listDaughters)) {
            continue;
          }
          if (iD == Ds1ToDStarK0s &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kK0Short, +o2::constants::physics::Pdg::kDStar}, false, nullptr, 2) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+kK0Short, -o2::constants::physics::Pdg::kDStar}, false, nullptr, 2)) {
            continue;
          }
          if (iD == Ds2StarToDPlusK0s &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kK0Short, +o2::constants::physics::Pdg::kDPlus}, false, nullptr, 2) &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, -PDGArrayParticle[iD], std::array{+kK0Short, -o2::constants::physics::Pdg::kDPlus}, false, nullptr, 2)) {
            continue;
          }
          if (iD == XiCplusToXiPiPi &&
              !RecoDecay::isMatchedMCGen(mcParticles, particle, PDGArrayParticle[iD], std::array{+kXiMinus, +kPiPlus, +kPiPlus}, true, nullptr, 2)) {
            continue;
          }
        }

        // Check momentum conservation
        double sumPxDau = 0.;
        double sumPyDau = 0.;
        double sumPzDau = 0.;
        bool momentumCheck = true;
        for (const int listDaughter : listDaughters) {
          auto daughter = mcParticles.rawIteratorAt(listDaughter - mcParticles.offset());
          sumPxDau += daughter.px();
          sumPyDau += daughter.py();
          sumPzDau += daughter.pz();
        }
        auto pxDiff = particle.px() - sumPxDau;
        auto pyDiff = particle.py() - sumPyDau;
        auto pzDiff = particle.pz() - sumPzDau;
        if (std::abs(pxDiff) > 0.001 || std::abs(pyDiff) > 0.001 || std::abs(pzDiff) > 0.001) {
          momentumCheck = false;
        }
        double const pDiff = RecoDecay::p(pxDiff, pyDiff, pzDiff);
        double const ptDiff = RecoDecay::pt(pxDiff, pyDiff);
        registry.fill(HIST("hMomentumCheck"), static_cast<float>(momentumCheck));
        registry.fill(HIST("hPxDiffMotherDaughterGen"), pxDiff);
        registry.fill(HIST("hPyDiffMotherDaughterGen"), pyDiff);
        registry.fill(HIST("hPzDiffMotherDaughterGen"), pzDiff);
        registry.fill(HIST("hPDiffMotherDaughterGen"), pDiff);
        registry.fill(HIST("hPtDiffMotherDaughterGen"), ptDiff);

        int origin{0};
        if (iD < NCharmMesonChannels || iD >= NCharmMesonChannels + NBeautyChannels) { // Charm hadrons
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle);
          if (origin == RecoDecay::OriginType::Prompt) {
            counterPrompt[iD]++;
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            counterNonPrompt[iD]++;
          }
        }

        auto daughter0 = particle.template daughters_as<Particles>().begin();
        double const vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
        double const vertexPrimary[3] = {mcCollision.posX(), mcCollision.posY(), mcCollision.posZ()};
        auto decayLength = RecoDecay::distance(vertexPrimary, vertexDau);
        if (iD < NCharmMesonChannels) {
          if (origin == RecoDecay::OriginType::Prompt) { // Prompt charm mesons
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("PromptCharmMesons/hPromptMesonsPtDistr"), iD, particle.pt());
              registry.fill(HIST("PromptCharmMesons/hPromptMesonsPtCentDistr"), iD, particle.pt(), centrality);
              if (storeOccupancy) {
                registry.fill(HIST("PromptCharmMesons/hPromptMesonsPtOccDistr"), iD, particle.pt(), occupancy);
              }
            }
            registry.fill(HIST("PromptCharmMesons/hPromptMesonsYDistr"), iD, particle.y());
            registry.fill(HIST("PromptCharmMesons/hPromptMesonsDecLenDistr"), iD, decayLength * 10000);
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsPtDistr"), iD, particle.pt());
              registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsPtCentDistr"), iD, particle.pt(), centrality);
              if (storeOccupancy) {
                registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsPtOccDistr"), iD, particle.pt(), occupancy);
              }
            }
            registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsYDistr"), iD, particle.y());
            registry.fill(HIST("NonPromptCharmMesons/hNonPromptMesonsDecLenDistr"), iD, decayLength * 10000);
          }
        } else if (iD < NCharmMesonChannels + NBeautyChannels) { // Beauty mesons
          if (std::abs(particle.y()) < 0.5) {
            registry.fill(HIST("Beauty/hPtDistr"), iD - NCharmMesonChannels, particle.pt());
            registry.fill(HIST("Beauty/hPtCentDistr"), iD - NCharmMesonChannels, particle.pt(), centrality);
          }
          registry.fill(HIST("Beauty/hYDistr"), iD - NCharmMesonChannels, particle.y());
          registry.fill(HIST("Beauty/hDecLenDistr"), iD - NCharmMesonChannels, decayLength * 10000);
        } else { // Charm baryons
          if (origin == RecoDecay::OriginType::Prompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsPtDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt());
              registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsPtCentDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt(), centrality);
              if (storeOccupancy) {
                registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsPtOccDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt(), occupancy);
              }
            }
            registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsYDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.y());
            registry.fill(HIST("PromptCharmBaryons/hPromptBaryonsDecLenDistr"), iD - NCharmMesonChannels - NBeautyChannels, decayLength * 10000);
          } else if (origin == RecoDecay::OriginType::NonPrompt) {
            if (std::abs(particle.y()) < 0.5) {
              registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt());
              registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtCentDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt(), centrality);
              if (storeOccupancy) {
                registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsPtOccDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.pt(), occupancy);
              }
            }
            registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsYDistr"), iD - NCharmMesonChannels - NBeautyChannels, particle.y());
            registry.fill(HIST("NonPromptCharmBaryons/hNonPromptBaryonsDecLenDistr"), iD - NCharmMesonChannels - NBeautyChannels, decayLength * 10000);
          }
        }
      }
    } // end particles

    registry.fill(HIST("Quarks/hCountC"), cPerCollision);
    registry.fill(HIST("Quarks/hCountB"), bPerCollision);
    registry.fill(HIST("Quarks/hCountCbar"), cBarPerCollision);
    registry.fill(HIST("Quarks/hCountBbar"), bBarPerCollision);
  }

  void processNoCentSel(aod::McCollisions const& mcCollisions,
                        aod::McParticles const& mcParticles,
                        CollisionsNoCents const& recoCollisions,
                        BCsInfo const& bcInfo)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollision, mcCollision.globalIndex());
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      std::array<int, NChannels> counterPrompt{0}, counterNonPrompt{0};
      runCheckGenParticles<o2::hf_centrality::CentralityEstimator::None>(mcCollision, mcParticlesPerMcColl, recoCollsPerMcColl, bcInfo, counterPrompt, counterNonPrompt);
      static_for<0, NCharmMesonChannels - 1>([&](auto i) { // Charm mesons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmMesons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmMesons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
      static_for<NCharmMesonChannels, NCharmMesonChannels + NBeautyChannels - 1>([&](auto i) { // Beauty hadrons
        constexpr int Index = i.value;
        registry.fill(HIST("Beauty/hCount") + HIST(ParticleNames[Index]), counterPrompt[Index]);
      });
      static_for<NCharmMesonChannels + NBeautyChannels, NChannels - 1>([&](auto i) { // Charm baryons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmBaryons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmBaryons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
    }
  } // end processNoCentSel
  PROCESS_SWITCH(HfTaskMcValidationGen, processNoCentSel, "Process generated collisions information without centrality selection", true);

  void processCentFT0C(aod::McCollisions const& mcCollisions,
                       aod::McParticles const& mcParticles,
                       CollisionsFT0Cs const& recoCollisions,
                       BCsInfo const& bcInfo)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollisionFT0C, mcCollision.globalIndex());
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      std::array<int, NChannels> counterPrompt{0}, counterNonPrompt{0};
      runCheckGenParticles<o2::hf_centrality::CentralityEstimator::FT0C>(mcCollision, mcParticlesPerMcColl, recoCollsPerMcColl, bcInfo, counterPrompt, counterNonPrompt);
      static_for<0, NCharmMesonChannels - 1>([&](auto i) { // Charm mesons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmMesons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmMesons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
      static_for<NCharmMesonChannels, NCharmMesonChannels + NBeautyChannels - 1>([&](auto i) { // Beauty
        constexpr int Index = i.value;
        registry.fill(HIST("Beauty/hCount") + HIST(ParticleNames[Index]), counterPrompt[Index]);
      });
      static_for<NCharmMesonChannels + NBeautyChannels, NChannels - 1>([&](auto i) { // Charm baryons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmBaryons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmBaryons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
    }
  } // end processCentFT0C
  PROCESS_SWITCH(HfTaskMcValidationGen, processCentFT0C, "Process generated collisions information with centrality selection using FT0C", false);

  void processCentFT0M(McCollisionsCentFT0Ms const& mcCollisions,
                       aod::McParticles const& mcParticles,
                       CollisionsFT0Ms const& recoCollisions,
                       BCsInfo const& bcInfo)
  {
    for (const auto& mcCollision : mcCollisions) {
      const auto recoCollsPerMcColl = recoCollisions.sliceBy(colPerMcCollisionFT0M, mcCollision.globalIndex());
      const auto mcParticlesPerMcColl = mcParticles.sliceBy(mcParticlesPerMcCollision, mcCollision.globalIndex());
      std::array<int, NChannels> counterPrompt{0}, counterNonPrompt{0};
      runCheckGenParticles<o2::hf_centrality::CentralityEstimator::FT0M>(mcCollision, mcParticlesPerMcColl, recoCollsPerMcColl, bcInfo, counterPrompt, counterNonPrompt);
      static_for<0, NCharmMesonChannels - 1>([&](auto i) { // Charm mesons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmMesons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmMesons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
      static_for<NCharmMesonChannels, NCharmMesonChannels + NBeautyChannels - 1>([&](auto i) { // Beauty mesons
        constexpr int Index = i.value;
        registry.fill(HIST("Beauty/hCount") + HIST(ParticleNames[Index]), counterPrompt[Index]);
      });
      static_for<NCharmMesonChannels + NBeautyChannels, NChannels - 1>([&](auto i) { // Charm baryons
        constexpr int Index = i.value;
        registry.fill(HIST("PromptCharmBaryons/hCountPrompt") + HIST(ParticleNames[Index]), counterPrompt[Index]);
        registry.fill(HIST("NonPromptCharmBaryons/hCountNonPrompt") + HIST(ParticleNames[Index]), counterNonPrompt[Index]);
      });
    }
  } // end processCentFT0M
  PROCESS_SWITCH(HfTaskMcValidationGen, processCentFT0M, "Process generated collisions information with centrality selection using FT0M", false);
};

/// Reconstruction Level Validation
///
///   - Gen-Rec Level Momentum Difference per component;
///   - Gen-Rec Level Difference for secondary-vertex coordinates and decay length;
struct HfTaskMcValidationRec {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;

  Configurable<int> eventGeneratorType{"eventGeneratorType", -1, "If positive, enable event selection using subGeneratorId information. The value indicates which events to keep (0 = MB, 4 = charm triggered, 5 = beauty triggered)"};
  Configurable<bool> storeOccupancy{"storeOccupancy", false, "Store collision occupancy for dedicated studies"};

  std::array<std::shared_ptr<TH1>, NChannels> histDeltaPt, histDeltaPx, histDeltaPy, histDeltaPz, histDeltaSecondaryVertexX, histDeltaSecondaryVertexY, histDeltaSecondaryVertexZ, histDeltaDecayLength;
  std::array<std::array<std::shared_ptr<TH2>, 2>, NChannels> histPtCentReco;
  std::array<std::array<std::shared_ptr<TH2>, 2>, NChannels> histPtOccReco;
  std::array<std::array<std::array<std::shared_ptr<TH1>, 5>, 2>, NChannels> histPtDau, histEtaDau, histImpactParameterDau;
  std::array<std::shared_ptr<THnSparse>, 4> histOriginTracks;
  std::shared_ptr<TH2> histAmbiguousTracks, histTracks;
  std::shared_ptr<TH1> histContributors;

  using HfCand2ProngWithMCRec = soa::Join<aod::HfCand2Prong, aod::HfCand2ProngMcRec>;
  using HfCand3ProngWithMCRec = soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec>;
  using CandMcGen = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen>;
  using CollisionsWithMCLabels = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using CollisionsWithMCLabelsAndCentFT0C = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, CentFT0Cs>;
  using CollisionsWithMCLabelsAndCentFT0M = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels, CentFT0Ms>;
  using TracksWithSel = soa::Join<aod::TracksWMc, aod::TracksExtra, aod::TrackSelection, aod::TrackCompColls>;

  Partition<TracksWithSel> tracksFilteredGlobalTrackWoDCA = requireGlobalTrackWoDCAInFilter();
  Partition<TracksWithSel> tracksInAcc = requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks);

  Preslice<HfCand2ProngWithMCRec> cand2ProngPerCollision = aod::hf_cand::collisionId;
  Preslice<HfCand3ProngWithMCRec> cand3ProngPerCollision = aod::hf_cand::collisionId;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HfEventSelection hfEvSel; // event selection and monitoring

  AxisSpec axisDeltaMom{2000, -1., 1.};
  AxisSpec axisOrigin{4, -1.5, 2.5};
  AxisSpec axisEta{40, -1., 1.};
  AxisSpec axisPt{50, 0., 10.};
  AxisSpec axisPtD{100, 0., 50.};
  AxisSpec axisCent{110, 0., 110.};
  AxisSpec axisOcc{3000, 0., 15000.};
  AxisSpec axisDeltaVtx{200, -1, 1.};
  AxisSpec axisDecision{2, -0.5, 1.5};
  AxisSpec axisITShits{8, -0.5, 7.5};
  AxisSpec axisMult{200, 0., 200.};
  AxisSpec axisR{100, 0., 0.5};
  AxisSpec axisSmallNum{20, -0.5, 19.5};

  HistogramRegistry registry{
    "registry",
    {{"histNtracks", "Number of global tracks w/o DCA requirement;#it{N}_{tracks};entries", {HistType::kTH1F, {axisMult}}},
     {"hNevReco", "Reconstructed events counter; Reco. events; entries", {HistType::kTH1F, {{1, -0.5, +0.5}}}},
     {"histXvtxReco", "Position of reco PV in #it{X};#it{X}^{reco} (cm);entries", {HistType::kTH1F, {axisDeltaVtx}}},
     {"histYvtxReco", "Position of reco PV in #it{Y};#it{Y}^{reco} (cm);entries", {HistType::kTH1F, {axisDeltaVtx}}},
     {"histZvtxReco", "Position of reco PV in #it{Z};#it{Z}^{reco} (cm);entries", {HistType::kTH1F, {{200, -20, 20.}}}},
     {"histDeltaZvtx", "Residual distribution of PV in #it{Z} as a function of number of contributors;number of contributors;#it{Z}^{reco} - #it{Z}^{gen} (cm);entries", {HistType::kTH2F, {{100, -0.5, 99.5}, {1000, -0.5, 0.5}}}},
     {"TrackToCollChecks/histAmbiguousTrackNumCollisions", "Number of collisions associated to an ambiguous track;number of collisions;entries", {HistType::kTH1F, {{30, -0.5, 29.5}}}},
     {"TrackToCollChecks/histAmbiguousTrackZvtxRMS", "RMS of #it{Z}^{reco} of collisions associated to a track;RMS(#it{Z}^{reco}) (cm);entries", {HistType::kTH1F, {{100, 0., 0.5}}}},
     {"TrackToCollChecks/histFracGoodContributors", "Fraction of PV contributors originating from the correct collision;fraction;entries", {HistType::kTH1F, {{101, 0., 1.01}}}},
     {"TrackToCollChecks/histCollisionsSameBC", "Collisions in same BC;number of contributors collision 1;number of contributors collision 2;#it{R}_{xy} collision 1 (cm);#it{R}_{xy} collision 2 (cm);number of contributors from beauty collision 1;number of contributors from beauty collision 2;", {HistType::kTHnSparseF, {axisMult, axisMult, axisR, axisR, axisSmallNum, axisSmallNum}}}}};
  HistogramRegistry registryMesons{"registryMesons"};
  HistogramRegistry registryBaryons{"registryBaryons"};

  /// RMS calculation
  /// \param vec  vector of values to compute RMS
  template <typename T>
  T computeRMS(std::vector<T>& vec)
  {
    T sum = std::accumulate(vec.begin(), vec.end(), 0.0);
    T mean = sum / vec.size();

    std::vector<T> diff(vec.size());
    std::transform(vec.begin(), vec.end(), diff.begin(), [mean](T x) { return x - mean; });
    T sqSum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    T stdev = std::sqrt(sqSum / vec.size());

    return stdev;
  }

  /// Fill MC histograms with decay properties at reconstruction level
  /// \param candidate is candidate
  /// \param mother is mother particle
  /// \param whichHad int indicating charm-hadron and decay channel, see enum DecayChannels
  /// \param whichOrigin int indicating origin: prompt or non-prompt
  /// \param centrality is collision centrality
  /// \param occupancy is collision occupancy
  template <typename T, typename U>
  void fillHisto(const T& candidate, const U& mother, int whichHad, int whichOrigin, float centrality, int occupancy)
  {
    histDeltaPt[whichHad]->Fill(candidate.pt() - mother.pt());
    histDeltaPx[whichHad]->Fill(candidate.px() - mother.px());
    histDeltaPy[whichHad]->Fill(candidate.py() - mother.py());
    histDeltaPz[whichHad]->Fill(candidate.pz() - mother.pz());
    // Compare Secondary vertex and decay length with MC
    auto daughter0 = mother.template daughters_as<aod::McParticles>().begin();
    double const vertexDau[3] = {daughter0.vx(), daughter0.vy(), daughter0.vz()};
    double const vertexMoth[3] = {mother.vx(), mother.vy(), mother.vz()};
    auto decayLength = RecoDecay::distance(vertexMoth, vertexDau);

    histDeltaSecondaryVertexX[whichHad]->Fill(candidate.xSecondaryVertex() - vertexDau[0]);
    histDeltaSecondaryVertexY[whichHad]->Fill(candidate.ySecondaryVertex() - vertexDau[1]);
    histDeltaSecondaryVertexZ[whichHad]->Fill(candidate.zSecondaryVertex() - vertexDau[2]);
    histDeltaDecayLength[whichHad]->Fill(candidate.decayLength() - decayLength);
    std::array<double, 3> const momDau0 = {candidate.pxProng0(),
                                           candidate.pyProng0(),
                                           candidate.pzProng0()};
    std::array<double, 3> const momDau1 = {candidate.pxProng1(),
                                           candidate.pyProng1(),
                                           candidate.pzProng1()};
    histPtCentReco[whichHad][whichOrigin]->Fill(candidate.pt(), centrality);
    if (storeOccupancy) {
      histPtOccReco[whichHad][whichOrigin]->Fill(candidate.pt(), occupancy);
    }
    histPtDau[whichHad][whichOrigin][0]->Fill(RecoDecay::pt(momDau0));
    histEtaDau[whichHad][whichOrigin][0]->Fill(RecoDecay::eta(momDau0));
    histImpactParameterDau[whichHad][whichOrigin][0]->Fill(candidate.impactParameter0());
    histPtDau[whichHad][whichOrigin][1]->Fill(RecoDecay::pt(momDau1));
    histEtaDau[whichHad][whichOrigin][1]->Fill(RecoDecay::eta(momDau1));
    histImpactParameterDau[whichHad][whichOrigin][1]->Fill(candidate.impactParameter1());
  }

  void init(InitContext&)
  {
    std::array<bool, 3> procCollisions = {doprocessColl, doprocessCollWithCentFTOC, doprocessCollWithCentFTOM};
    if (std::accumulate(procCollisions.begin(), procCollisions.end(), 0) > 1) {
      LOGP(fatal, "At most one process function for collision study can be enabled at a time.");
    }

    std::array<bool, 3> procCollAccoc = {doprocessCollAssoc, doprocessCollAssocWithCentFTOC, doprocessCollAssocWithCentFTOM};
    if (std::accumulate(procCollAccoc.begin(), procCollAccoc.end(), 0) > 1) {
      LOGP(fatal, "At most one process for collision association study function can be enabled at a time.");
    }

    histOriginTracks[0] = registry.add<THnSparse>("TrackToCollChecks/histOriginNonAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});           // tracks not associated to any collision
    histOriginTracks[1] = registry.add<THnSparse>("TrackToCollChecks/histOriginAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});              // tracks associasted to a collision
    histOriginTracks[2] = registry.add<THnSparse>("TrackToCollChecks/histOriginGoodAssociatedTracks", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits});          // tracks associated to the correct collision considering only first reco collision (based on the MC collision index)
    histOriginTracks[3] = registry.add<THnSparse>("TrackToCollChecks/histOriginGoodAssociatedTracksAmbiguous", ";origin;#it{p}_{T}^{reco} (GeV/#it{c});#it{#eta}^{reco};#it{Z}_{vtx}^{reco}#minus#it{Z}_{vtx}^{gen} (cm); is PV contributor; has TOF; number of ITS hits", HistType::kTHnSparseF, {axisOrigin, axisPt, axisEta, axisDeltaVtx, axisDecision, axisDecision, axisITShits}); // tracks associated to the correct collision considering all ambiguous reco collisions (based on the MC collision index)
    for (std::size_t iHist{0}; iHist < histOriginTracks.size(); ++iHist) {
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(1, "no MC particle");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(2, "no quark");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(3, "charm");
      histOriginTracks[iHist]->GetAxis(0)->SetBinLabel(4, "beauty");
    }
    histAmbiguousTracks = registry.add<TH2>("TrackToCollChecks/histAmbiguousTracks", "Tracks that are ambiguous vs. origin;#it{p}_{T}^{reco} (GeV/#it{c});entries", HistType::kTH2F, {axisOrigin, axisPt});
    histTracks = registry.add<TH2>("histTracks", "Tracks vs. origin;#it{p}_{T}^{reco} (GeV/#it{c});entries", HistType::kTH2F, {axisOrigin, axisPt});
    histTracks->GetXaxis()->SetBinLabel(1, "no MC particle");
    histTracks->GetXaxis()->SetBinLabel(2, "no quark");
    histTracks->GetXaxis()->SetBinLabel(3, "charm");
    histTracks->GetXaxis()->SetBinLabel(4, "beauty");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(1, "no MC particle");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(2, "no quark");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(3, "charm");
    histAmbiguousTracks->GetXaxis()->SetBinLabel(4, "beauty");
    for (auto iHad = 0; iHad < NChannels; ++iHad) {
      if (iHad < NCharmMesonChannels) {
        histDeltaPt[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaPt", ParticleNames[iHad].data()), Form("Pt difference reco - MC %s; #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPx[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaPx", ParticleNames[iHad].data()), Form("Px difference reco - MC %s; #it{p}_{x}^{reco} - #it{p}_{x}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPy[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaPy", ParticleNames[iHad].data()), Form("Py difference reco - MC %s; #it{p}_{y}^{reco} - #it{p}_{y}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPz[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaPz", ParticleNames[iHad].data()), Form("Pz difference reco - MC %s; #it{p}_{z}^{reco} - #it{p}_{z}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaSecondaryVertexX[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaSecondaryVertexX", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta x (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaSecondaryVertexY[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaSecondaryVertexY", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta y (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaSecondaryVertexZ[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaSecondaryVertexZ", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta z (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaDecayLength[iHad] = registryMesons.add<TH1>(Form("%s/histDeltaDecayLength", ParticleNames[iHad].data()), Form("Decay length difference reco - MC (%s); #Delta L (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        for (auto iOrigin = 0; iOrigin < 2; ++iOrigin) {
          histPtCentReco[iHad][iOrigin] = registryMesons.add<TH2>(Form("%s/histPtCentReco%s", ParticleNames[iHad].data(), OriginNames[iOrigin].data()), Form("Pt Cent reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); Centrality (%%); entries", OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH2F, {axisPtD, axisCent});
          if (storeOccupancy) {
            histPtOccReco[iHad][iOrigin] = registryMesons.add<TH2>(Form("%s/histPtOccReco%s", ParticleNames[iHad].data(), OriginNames[iOrigin].data()), Form("Pt Cent reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); Occupancy; entries", OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH2F, {axisPtD, axisOcc});
          }
          for (unsigned int iDau = 0; iDau < NDaughters[iHad]; ++iDau) {
            histPtDau[iHad][iOrigin][iDau] = registryMesons.add<TH1>(Form("%s/histPtDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d Pt reco - %s %s; #it{p}_{T}^{dau, reco} (GeV/#it{c}); entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {axisPt});
            histEtaDau[iHad][iOrigin][iDau] = registryMesons.add<TH1>(Form("%s/histEtaDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d Eta reco - %s %s; #it{#eta}^{dau, reco}; entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {{100, -1., 1.}});
            histImpactParameterDau[iHad][iOrigin][iDau] = registryMesons.add<TH1>(Form("%s/histImpactParameterDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d DCAxy reco - %s %s; DCAxy (cm); entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
          }
        }
      } else if (iHad >= NCharmMesonChannels + NBeautyChannels) {
        histDeltaPt[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaPt", ParticleNames[iHad].data()), Form("Pt difference reco - MC %s; #it{p}_{T}^{reco} - #it{p}_{T}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPx[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaPx", ParticleNames[iHad].data()), Form("Px difference reco - MC %s; #it{p}_{x}^{reco} - #it{p}_{x}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPy[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaPy", ParticleNames[iHad].data()), Form("Py difference reco - MC %s; #it{p}_{y}^{reco} - #it{p}_{y}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaPz[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaPz", ParticleNames[iHad].data()), Form("Pz difference reco - MC %s; #it{p}_{z}^{reco} - #it{p}_{z}^{gen} (GeV/#it{c}); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaMom});
        histDeltaSecondaryVertexX[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaSecondaryVertexX", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta x (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaSecondaryVertexY[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaSecondaryVertexY", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta y (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaSecondaryVertexZ[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaSecondaryVertexZ", ParticleNames[iHad].data()), Form("Sec. Vertex difference reco - MC (MC matched) - %s; #Delta z (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        histDeltaDecayLength[iHad] = registryBaryons.add<TH1>(Form("%s/histDeltaDecayLength", ParticleNames[iHad].data()), Form("Decay length difference reco - MC (%s); #Delta L (cm); entries", Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
        for (auto iOrigin = 0; iOrigin < 2; ++iOrigin) {
          histPtCentReco[iHad][iOrigin] = registryBaryons.add<TH2>(Form("%s/histPtCentReco%s", ParticleNames[iHad].data(), OriginNames[iOrigin].data()), Form("Pt Cent reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); Centrality (%%); entries", OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH2F, {axisPtD, axisCent});
          if (storeOccupancy) {
            histPtOccReco[iHad][iOrigin] = registryBaryons.add<TH2>(Form("%s/histPtOccReco%s", ParticleNames[iHad].data(), OriginNames[iOrigin].data()), Form("Pt Cent reco %s %s; #it{p}_{T}^{reco} (GeV/#it{c}); Occupancy; entries", OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH2F, {axisPtD, axisOcc});
          }
          for (unsigned int iDau = 0; iDau < NDaughters[iHad]; ++iDau) {
            histPtDau[iHad][iOrigin][iDau] = registryBaryons.add<TH1>(Form("%s/histPtDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d Pt reco - %s %s; #it{p}_{T}^{dau, reco} (GeV/#it{c}); entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {axisPt});
            histEtaDau[iHad][iOrigin][iDau] = registryBaryons.add<TH1>(Form("%s/histEtaDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d Eta reco - %s %s; #it{#eta}^{dau, reco}; entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {{100, -1., 1.}});
            histImpactParameterDau[iHad][iOrigin][iDau] = registryBaryons.add<TH1>(Form("%s/histImpactParameterDau%d%s", ParticleNames[iHad].data(), iDau, OriginNames[iOrigin].data()), Form("Daughter %d DCAxy reco - %s %s; DCAxy (cm); entries", iDau, OriginNames[iOrigin].data(), Labels[iHad].data()), HistType::kTH1F, {axisDeltaVtx});
          }
        }
      }
    }
    histContributors = registry.add<TH1>("TrackToCollChecks/histContributors", "PV contributors from correct/wrong MC collision;;entries", HistType::kTH1F, {axisDecision});
    histContributors->GetXaxis()->SetBinLabel(1, "correct MC collision");
    histContributors->GetXaxis()->SetBinLabel(2, "wrong MC collision");
    hfEvSel.addHistograms(registry); // collision monitoring

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename Coll>
  void checkCollisions(Coll const& collision,
                       aod::McCollisions const&,
                       aod::BCsWithTimestamps const&)
  {
    // apply event selection
    if (!collision.has_mcCollision()) {
      return;
    }

    float centrality{-1.f};
    const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
    if (rejectionMask != 0) {
      /// at least one event selection not satisfied --> reject the candidate
      return;
    }

    auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
    if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
      return;
    }

    registry.fill(HIST("hNevReco"), 1);
    registry.fill(HIST("histXvtxReco"), collision.posX());
    registry.fill(HIST("histYvtxReco"), collision.posY());
    registry.fill(HIST("histZvtxReco"), collision.posZ());
    registry.fill(HIST("histDeltaZvtx"), collision.numContrib(), collision.posZ() - mcCollision.posZ());
  }

  template <o2::hf_centrality::CentralityEstimator CentEstimator, typename Colls>
  void checkCollisionAssociation(Colls const& collisions,
                                 TracksWithSel const&,
                                 aod::McParticles const& mcParticles,
                                 aod::McCollisions const&,
                                 aod::BCsWithTimestamps const&)
  {
    // loop over collisions
    for (const auto& collision : collisions) {
      // check that collision is selected by hf-track-index-skim-creator-tag-sel-collisions

      float centrality{-1.f};
      const auto rejectionMask = hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry);
      if (rejectionMask != 0) {
        /// at least one event selection not satisfied --> reject the candidate
        continue;
      }

      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
        continue;
      }
      auto tracksGlobalWoDCAColl1 = tracksFilteredGlobalTrackWoDCA->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      registry.fill(HIST("histNtracks"), tracksGlobalWoDCAColl1.size());
      auto tracksColl1 = tracksInAcc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      int nContributors = 0, nGoodContributors = 0;
      for (const auto& track : tracksColl1) {
        if (!track.isPVContributor()) {
          continue;
        }
        if (!track.has_mcParticle()) {
          continue;
        }
        nContributors++;
        auto particle = track.mcParticle();
        if (collision.mcCollisionId() == particle.mcCollisionId()) {
          nGoodContributors++;
        }
      }
      float const frac = (nContributors > 0) ? static_cast<float>(nGoodContributors) / nContributors : 1.;
      registry.fill(HIST("TrackToCollChecks/histFracGoodContributors"), frac);
      uint64_t const mostProbableBC = collision.bc().globalBC();
      for (auto collision2 = collision + 1; collision2 != collisions.end(); ++collision2) {
        uint64_t const mostProbableBC2 = collision2.bc().globalBC();
        if (mostProbableBC2 == mostProbableBC) {
          float const radColl1 = std::sqrt(collision.posX() * collision.posX() + collision.posY() * collision.posY());
          float const radColl2 = std::sqrt(collision2.posX() * collision2.posX() + collision2.posY() * collision2.posY());
          int nFromBeautyColl1 = 0, nFromBeautyColl2 = 0;
          for (const auto& trackColl1 : tracksColl1) {
            if (trackColl1.has_mcParticle() && trackColl1.isPVContributor()) {
              auto particleColl1 = trackColl1.mcParticle();
              auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particleColl1, true);
              if (origin == RecoDecay::NonPrompt) {
                nFromBeautyColl1++;
              }
            }
          }
          auto tracksColl2 = tracksInAcc->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);
          for (const auto& trackColl2 : tracksColl2) {
            if (trackColl2.has_mcParticle() && trackColl2.isPVContributor()) {
              auto particleColl2 = trackColl2.mcParticle();
              auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particleColl2, true);
              if (origin == RecoDecay::NonPrompt) {
                nFromBeautyColl2++;
              }
            }
          }
          registry.fill(HIST("TrackToCollChecks/histCollisionsSameBC"), collision.numContrib(), collision2.numContrib(), radColl1, radColl2, nFromBeautyColl1, nFromBeautyColl2);
          break;
        }
      }
    }

    // loop over tracks
    for (const auto& track : tracksFilteredGlobalTrackWoDCA) {
      // check number of ITS hits
      int nITSlayers = 0;
      uint8_t const itsHitMap = track.itsClusterMap();
      for (int iLayer = 0; iLayer < 7; ++iLayer) {
        if (TESTBIT(itsHitMap, iLayer)) {
          nITSlayers++;
        }
      }
      uint const index = uint(track.collisionId() >= 0);
      if (track.has_mcParticle()) {
        const auto& particle = track.mcParticle(); // get corresponding MC particle to check origin
        const auto& mcCollision = particle.mcCollision_as<aod::McCollisions>();
        if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
          continue;
        }
        auto origin = RecoDecay::getCharmHadronOrigin(mcParticles, particle, true);
        histTracks->Fill(origin, track.pt());
        bool const isAmbiguous = (track.compatibleCollIds().size() != 1);
        if (isAmbiguous) {
          registry.fill(HIST("TrackToCollChecks/histAmbiguousTrackNumCollisions"), track.compatibleCollIds().size());
          histAmbiguousTracks->Fill(origin, track.pt());
          std::vector<double> ambCollPosZ{};
          for (const auto& collIdx : track.compatibleCollIds()) {
            const auto& ambCollision = collisions.rawIteratorAt(collIdx);
            ambCollPosZ.push_back(ambCollision.posZ());
          }
          // here we are only interested to tracks associated to multiple vertices
          if (!ambCollPosZ.empty()) {
            registry.fill(HIST("TrackToCollChecks/histAmbiguousTrackZvtxRMS"), computeRMS(ambCollPosZ));
          }
        }
        float deltaZ = -999.f;
        if (index) {
          const auto& collision = track.collision_as<CollisionsWithMCLabels>();
          deltaZ = collision.posZ() - mcCollision.posZ();
          if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
            histOriginTracks[index + 1]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
          } else { // if the default associated collision is not the good one, check if the tracks is ambiguous
            if (isAmbiguous) {
              for (const auto& collIdx : track.compatibleCollIds()) {
                auto ambCollision = collisions.rawIteratorAt(collIdx);

                if (ambCollision.has_mcCollision() && ambCollision.mcCollisionId() == particle.mcCollisionId()) {
                  histOriginTracks[index + 2]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
                  break;
                }
              }
            }
          }
          if (track.isPVContributor()) {
            if (collision.has_mcCollision() && collision.mcCollisionId() == particle.mcCollisionId()) {
              histContributors->Fill(0);
            } else {
              histContributors->Fill(1);
            }
          }
        }
        histOriginTracks[index]->Fill(origin, track.pt(), track.eta(), deltaZ, track.isPVContributor(), track.hasTOF(), nITSlayers);
      } else {
        histOriginTracks[index]->Fill(-1.f, track.pt(), track.eta(), -999.f, track.isPVContributor(), track.hasTOF(), nITSlayers);
      }
    }
  }

  void processColl(CollisionsWithMCLabels::iterator const& collision,
                   aod::McCollisions const& mcCollisions,
                   aod::BCsWithTimestamps const& bcs)
  {
    checkCollisions<o2::hf_centrality::CentralityEstimator::None>(collision, mcCollisions, bcs);
  } // end process
  PROCESS_SWITCH(HfTaskMcValidationRec, processColl, "Process collision information without centrality selection", true);

  void processCollWithCentFTOC(CollisionsWithMCLabelsAndCentFT0C::iterator const& collision,
                               aod::McCollisions const& mcCollisions,
                               aod::BCsWithTimestamps const& bcs)
  {
    checkCollisions<o2::hf_centrality::CentralityEstimator::FT0C>(collision, mcCollisions, bcs);
  } // end process
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollWithCentFTOC, "Process collision information with centrality selection with FT0C", false);

  void processCollWithCentFTOM(CollisionsWithMCLabelsAndCentFT0M::iterator const& collision,
                               aod::McCollisions const& mcCollisions,
                               aod::BCsWithTimestamps const& bcs)
  {
    checkCollisions<o2::hf_centrality::CentralityEstimator::FT0M>(collision, mcCollisions, bcs);
  } // end process
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollWithCentFTOM, "Process collision information with centrality selection with FT0M", false);

  void processCollAssoc(CollisionsWithMCLabels const& collisions,
                        TracksWithSel const& tracks,
                        aod::McParticles const& mcParticles,
                        aod::McCollisions const& mcCollisions,
                        aod::BCsWithTimestamps const& bcs)
  {
    checkCollisionAssociation<o2::hf_centrality::CentralityEstimator::None>(collisions, tracks, mcParticles, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollAssoc, "Process collision-association information, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  void processCollAssocWithCentFTOC(CollisionsWithMCLabelsAndCentFT0C const& collisions,
                                    TracksWithSel const& tracks,
                                    aod::McParticles const& mcParticles,
                                    aod::McCollisions const& mcCollisions,
                                    aod::BCsWithTimestamps const& bcs)
  {
    checkCollisionAssociation<o2::hf_centrality::CentralityEstimator::FT0C>(collisions, tracks, mcParticles, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollAssocWithCentFTOC, "Process collision-association information with centrality selection with FT0C, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  void processCollAssocWithCentFTOM(CollisionsWithMCLabelsAndCentFT0M const& collisions,
                                    TracksWithSel const& tracks,
                                    aod::McParticles const& mcParticles,
                                    aod::McCollisions const& mcCollisions,
                                    aod::BCsWithTimestamps const& bcs)
  {
    checkCollisionAssociation<o2::hf_centrality::CentralityEstimator::FT0M>(collisions, tracks, mcParticles, mcCollisions, bcs);
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processCollAssocWithCentFTOM, "Process collision-association information with centrality selection with FT0M, requires extra table from TrackToCollisionAssociation task (fillTableOfCollIdsPerTrack=true)", false);

  template <CentralityEstimator CentEstimator, typename Coll>
  void processEff(HfCand2ProngWithMCRec const& cand2Prongs,
                  HfCand3ProngWithMCRec const& cand3Prongs,
                  aod::TracksWMc const&,
                  aod::McParticles const& mcParticles,
                  aod::McCollisions const&,
                  aod::BCsWithTimestamps const&,
                  Coll const& collisions,
                  Preslice<HfCand2ProngWithMCRec> cand2ProngsPerCollision,
                  Preslice<HfCand3ProngWithMCRec> cand3ProngsPerCollision)
  {
    // loop over collisions
    for (const auto& collision : collisions) {
      // apply event selection
      float centrality{105.f};
      int const occupancy = collision.trackOccupancyInTimeRange();
      hfEvSel.getHfCollisionRejectionMask<true, CentEstimator, aod::BCsWithTimestamps>(collision, centrality, ccdb, registry); // only needed to update centrality, no bitmask selection applied
      if (!collision.has_mcCollision()) {
        return;
      }
      auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
      if (eventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != eventGeneratorType) {
        return;
      }

      // group 2- and 3-prongs for collision
      auto thisCollId = collision.globalIndex();
      auto grouped2ProngCandidates = cand2Prongs.sliceBy(cand2ProngsPerCollision, thisCollId);
      auto grouped3ProngCandidates = cand3Prongs.sliceBy(cand3ProngsPerCollision, thisCollId);

      // loop over 2-prong candidates
      for (const auto& cand2Prong : grouped2ProngCandidates) {

        // determine which kind of candidate it is
        bool const isD0Sel = TESTBIT(cand2Prong.hfflag(), o2::aod::hf_cand_2prong::DecayType::D0ToPiK);
        if (!isD0Sel) {
          continue;
        }
        int whichHad = -1;
        if (std::abs(cand2Prong.flagMcMatchRec()) == o2::hf_decay::hf_cand_2prong::DecayChannelMain::D0ToPiK) {
          whichHad = DzeroToKPi;
        }
        int whichOrigin;
        if (cand2Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
          whichOrigin = 0;
        } else {
          whichOrigin = 1;
        }

        if (whichHad >= 0) {
          int indexParticle = -1;
          indexParticle = RecoDecay::getMother(mcParticles, cand2Prong.template prong0_as<aod::TracksWMc>().template mcParticle_as<aod::McParticles>(), PDGArrayParticle[whichHad], true);
          if (indexParticle < 0) {
            continue;
          }
          auto mother = mcParticles.rawIteratorAt(indexParticle);
          fillHisto(cand2Prong, mother, whichHad, whichOrigin, centrality, occupancy);
        }
      } // end loop on 2-prong candidates

      // loop over 3-prong candidates
      for (const auto& cand3Prong : grouped3ProngCandidates) {

        // determine which kind of candidate it is
        // FIXME: add D* and decays with cascades
        bool const isDPlusSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::DplusToPiKPi);
        bool const isDsSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::DsToKKPi);
        bool const isLcSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::LcToPKPi);
        bool const isXicSel = TESTBIT(cand3Prong.hfflag(), hf_cand_3prong::DecayType::XicToPKPi);
        if (!isDPlusSel && !isDsSel && !isLcSel && !isXicSel) {
          continue;
        }
        int whichHad = -1;
        if (isDPlusSel && std::abs(cand3Prong.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DplusToPiKPi) {
          whichHad = DplusToPiKPi;
        } else if (isDsSel && std::abs(cand3Prong.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::DsToPiKK) {
          if (cand3Prong.flagMcDecayChanRec() == o2::hf_decay::hf_cand_3prong::DecayChannelResonant::DsToPhiPi) {
            whichHad = DsToPhiPiToKKPi;
          }
          if (cand3Prong.flagMcDecayChanRec() == o2::hf_decay::hf_cand_3prong::DecayChannelResonant::DsToKstar0K) {
            whichHad = DsToK0starKToKKPi;
          }
          if (cand3Prong.flagMcDecayChanRec() == o2::hf_decay::hf_cand_3prong::DecayChannelResonant::DplusToPhiPi) {
            whichHad = DplusToPhiPiToKKPi;
          }
        } else if (isLcSel && std::abs(cand3Prong.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi) {
          whichHad = LcToPKPi;
        } else if (isXicSel && std::abs(cand3Prong.flagMcMatchRec()) == o2::hf_decay::hf_cand_3prong::DecayChannelMain::XicToPKPi) {
          whichHad = XiCplusToPKPi;
        }
        int whichOrigin;
        if (cand3Prong.originMcRec() == RecoDecay::OriginType::Prompt) {
          whichOrigin = 0;
        } else {
          whichOrigin = 1;
        }

        if (whichHad >= 0) {
          int indexParticle = -1;
          if (cand3Prong.template prong0_as<aod::TracksWMc>().has_mcParticle()) {
            indexParticle = RecoDecay::getMother(mcParticles, cand3Prong.template prong0_as<aod::TracksWMc>().template mcParticle_as<aod::McParticles>(), PDGArrayParticle[whichHad], true);
          }
          if (indexParticle < 0) {
            continue;
          }
          auto mother = mcParticles.rawIteratorAt(indexParticle);
          fillHisto(cand3Prong, mother, whichHad, whichOrigin, centrality, occupancy);
          std::array<double, 3> const momDau2 = {cand3Prong.pxProng2(),
                                                 cand3Prong.pyProng2(),
                                                 cand3Prong.pzProng2()};
          histPtDau[whichHad][whichOrigin][2]->Fill(RecoDecay::pt(momDau2));
          histEtaDau[whichHad][whichOrigin][2]->Fill(RecoDecay::eta(momDau2));
          histImpactParameterDau[whichHad][whichOrigin][2]->Fill(cand3Prong.impactParameter2());
        }
      } // end loop on 3-prong candidates
    } // end loop on collisions
  }
  void processEffNoCent(HfCand2ProngWithMCRec const& cand2Prongs,
                        HfCand3ProngWithMCRec const& cand3Prongs,
                        aod::TracksWMc const& mcTracks,
                        aod::McParticles const& mcParticles,
                        aod::McCollisions const& mcCollisions,
                        aod::BCsWithTimestamps const& bcs,
                        CollisionsWithMCLabels const& collsWithLabels)
  {
    processEff<CentralityEstimator::None, CollisionsWithMCLabels>(cand2Prongs, cand3Prongs, mcTracks, mcParticles, mcCollisions, bcs, collsWithLabels, cand2ProngPerCollision, cand3ProngPerCollision);
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processEffNoCent, "Compute charm-hadron efficiencies (not all of them are implemented), requires HF candidate creators w/o information on centrality", false);

  void processEffCentFT0C(HfCand2ProngWithMCRec const& cand2Prongs,
                          HfCand3ProngWithMCRec const& cand3Prongs,
                          aod::TracksWMc const& mcTracks,
                          aod::McParticles const& mcParticles,
                          aod::McCollisions const& mcCollisions,
                          aod::BCsWithTimestamps const& bcs,
                          CollisionsWithMCLabelsAndCentFT0C const& collsWithLabels)
  {
    processEff<CentralityEstimator::FT0C, CollisionsWithMCLabelsAndCentFT0C>(cand2Prongs, cand3Prongs, mcTracks, mcParticles, mcCollisions, bcs, collsWithLabels, cand2ProngPerCollision, cand3ProngPerCollision);
  }
  PROCESS_SWITCH(HfTaskMcValidationRec, processEffCentFT0C, "Compute charm-hadron efficiencies (not all of them are implemented), requires HF candidate creators with information on centrality from FT0C", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskMcValidationGen>(cfgc),
    adaptAnalysisTask<HfTaskMcValidationRec>(cfgc)};
  return workflow;
}
