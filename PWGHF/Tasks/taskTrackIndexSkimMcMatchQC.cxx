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

/// \file taskTrackIndexSkimMcMatchQC.cxx
/// @brief process function to perform the candidate matching after the trackIndexSkimCreator
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN Padova
/// \author Fabrizio Grosa <fgroa@cern.ch>, CERN

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/RecoDecay.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

namespace
{
static const double massPi = RecoDecay::getMassPDG(kPiPlus);
static const double massK = RecoDecay::getMassPDG(kKPlus);
static const double massProton = RecoDecay::getMassPDG(kProton);
} // namespace

//____________________________________________________________________________________________________________________________________________
// Struct to match candidates from skimming to MC
struct HfTrackIndexSkimMcMatchQc {

  HistogramRegistry registry{"registry"};
  void init(InitContext const&)
  {
    /// 2-prong signal
    // D0(bar) → π± K∓
    registry.add("hMassPromptD0ToPiKvsPt", "prompt D^{0} signal;inv. mass (#pi K) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.4, 2.2}, {48, 0., 24.}}});
    registry.add("hMassNonPromptD0ToPiKvsPt", "non-prompt D^{0} signal;inv. mass (#pi K) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.4, 2.2}, {48, 0., 24.}}});
    /// 3-prong signal
    // D± → π± K∓ π±
    registry.add("hMassPromptDPlusToPiKPivsPt", "prompt D^{#plus} signal;inv. mass (#pi K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.4, 2.2}, {48, 0., 24.}}});
    registry.add("hMassNonPromptDPlusToPiKPivsPt", "non-prompt D^{#plus} signal;inv. mass (#pi K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.4, 2.2}, {48, 0., 24.}}});
    // Ds± → K± K∓ π±
    registry.add("hMassPromptDsToKKPivsPt", "prompt D_{s}^{#plus} signal;inv. mass (K K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.6, 2.4}, {48, 0., 24.}}});
    registry.add("hMassNonPromptDsToKKPivsPt", "nonprompt D_{s}^{#plus} signal;inv. mass (K K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.6, 2.4}, {48, 0., 24.}}});
    // Λc± → p± K∓ π±
    registry.add("hMassPromptLcToPKPivsPt", "prompt #Lambda_{c}^{#plus} signal;inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.7, 2.5}, {48, 0., 24.}}});
    registry.add("hMassNonPromptLcToPKPivsPt", "non-prompt #Lambda_{c}^{#plus} signal;inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 1.7, 2.5}, {48, 0., 24.}}});
    // Ξc± → p± K∓ π±
    registry.add("hMassPromptXicToPKPivsPt", "prompt #Xi_{c}^{#plus} signal;inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 2.0, 2.8}, {48, 0., 24.}}});
    registry.add("hMassNonPromptXicToPKPivsPt", "non-prompt #Xi_{c}^{#plus} signal;inv. mass (p K #pi) (GeV/#it{c}^{2});#it{p}_{T} (GeV/#it{c});", {HistType::kTH2F, {{400, 2.0, 2.8}, {48, 0., 24.}}});
  }

  using BigTracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels>;

  void processMcMatch(aod::Hf2Prongs const& cand2Prongs,
                      aod::Hf3Prongs const& cand3Prongs,
                      aod::McParticles const& particlesMC,
                      BigTracksMC const&)
  {
    /// 2-prong signal
    for (const auto& cand2Prong : cand2Prongs) { // start loop over 2 prongs

      auto trackPos = cand2Prong.prong0_as<BigTracksMC>(); // positive daughter
      auto trackNeg = cand2Prong.prong1_as<BigTracksMC>(); // negative daughter

      std::array<float, 3> pVecPos = {trackPos.px(), trackPos.py(), trackPos.pz()};
      std::array<float, 3> pVecNeg = {trackNeg.px(), trackNeg.py(), trackNeg.pz()};
      auto pVec2Prong = RecoDecay::pVec(pVecPos, pVecNeg);
      auto pt2Prong = RecoDecay::pt(pVec2Prong);

      auto invMassD0 = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massPi, massK});
      auto invMassD0bar = RecoDecay::m(std::array{pVecPos, pVecNeg}, std::array{massK, massPi});

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;

      // D0(bar) → π± K∓
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, std::array{trackPos, trackNeg}, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign);
      if (indexRec > -1) {
        auto particle = particlesMC.rawIteratorAt(indexRec);
        flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);

        // look for the correct mass hypothesis
        // either Kπ or πK
        auto invMass = invMassD0;
        if (trackPos.signed1Pt() > 0) {
          auto pdgCode = trackPos.mcParticle().pdgCode();
          if (pdgCode == kKPlus) {
            invMass = invMassD0bar;
          }
        } else if (trackNeg.signed1Pt() > 0) { // probably a not necessary check...
          auto pdgCode = trackNeg.mcParticle().pdgCode();
          if (pdgCode == kKPlus) {
            invMass = invMassD0bar;
          }
        }

        if (flag == RecoDecay::OriginType::Prompt) {
          registry.fill(HIST("hMassPromptD0ToPiKvsPt"), invMass, pt2Prong);
        } else if (flag == RecoDecay::OriginType::NonPrompt) {
          registry.fill(HIST("hMassNonPromptD0ToPiKvsPt"), invMass, pt2Prong);
        }
      }
    } // end loop over 2-prong candidates

    /// 3-prong signal
    for (const auto& cand3Prong : cand3Prongs) { // start loop over 3 prongs

      auto trackFirst = cand3Prong.prong0_as<BigTracksMC>();  // first daughter
      auto trackSecond = cand3Prong.prong1_as<BigTracksMC>(); // second daughter
      auto trackThird = cand3Prong.prong2_as<BigTracksMC>();  // third daughter
      auto arrayDaughters = std::array{trackFirst, trackSecond, trackThird};

      std::array<float, 3> pVecFirst = {trackFirst.px(), trackFirst.py(), trackFirst.pz()};
      std::array<float, 3> pVecSecond = {trackSecond.px(), trackSecond.py(), trackSecond.pz()};
      std::array<float, 3> pVecThird = {trackThird.px(), trackThird.py(), trackThird.pz()};

      auto pVec3Prong = RecoDecay::pVec(pVecFirst, pVecSecond, pVecThird);
      auto pt3Prong = RecoDecay::pt(pVec3Prong);

      auto invMassDplus = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massPi});

      auto invMassDsToKKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massK, massK, massPi});
      auto invMassDsToPiKK = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massK});

      auto invMassLcToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassLcToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      auto invMassXicToPKPi = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massProton, massK, massPi});
      auto invMassXicToPiKP = RecoDecay::m(std::array{pVecFirst, pVecSecond, pVecThird}, std::array{massPi, massK, massProton});

      int8_t sign = 0;
      int8_t flag = RecoDecay::OriginType::None;
      // int8_t channel = -1;

      // D± → π± K∓ π±
      auto indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kDPlus, array{+kPiPlus, -kKPlus, +kPiPlus}, true, &sign, 2);
      if (indexRec >= 0) {
        // channel = kDplus;
        if (indexRec > -1) {
          auto particle = particlesMC.rawIteratorAt(indexRec);
          flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);
          if (flag == RecoDecay::OriginType::Prompt) {
            registry.fill(HIST("hMassPromptDPlusToPiKPivsPt"), invMassDplus, pt3Prong);
          } else if (flag == RecoDecay::OriginType::NonPrompt) {
            registry.fill(HIST("hMassNonPromptDPlusToPiKPivsPt"), invMassDplus, pt3Prong);
          }
        }
      }
      if (indexRec < 0) {
        // Ds± → K± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, 431, array{+kKPlus, -kKPlus, +kPiPlus}, true, &sign, 2); // TODO: replace hard coded pdg code
        if (indexRec >= 0) {
          // channel = kDs;
          if (indexRec > -1) {
            auto particle = particlesMC.rawIteratorAt(indexRec);
            flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);

            // look for the correct mass hypothesis
            // either KKπ or πKK
            auto invMass = invMassDsToKKPi;
            if (std::abs(trackFirst.mcParticle().pdgCode()) == kPiPlus) {
              invMass = invMassDsToPiKK;
            }

            if (flag == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hMassPromptDsToKKPivsPt"), invMass, pt3Prong);
            } else if (flag == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hMassNonPromptDsToKKPivsPt"), invMass, pt3Prong);
            }
          }
        }
      }
      if (indexRec < 0) {
        // Λc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          // channel = kLc;
          if (indexRec > -1) {
            auto particle = particlesMC.rawIteratorAt(indexRec);
            flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);

            // look for the correct mass hypothesis
            // either pKπ or πKp
            auto invMass = invMassLcToPKPi;
            if (std::abs(trackFirst.mcParticle().pdgCode()) == kPiPlus) {
              invMass = invMassLcToPiKP;
            }

            if (flag == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hMassPromptLcToPKPivsPt"), invMass, pt3Prong);
            } else if (flag == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hMassNonPromptLcToPKPivsPt"), invMass, pt3Prong);
            }
          }
        }
      }
      if (indexRec < 0) {
        // Ξc± → p± K∓ π±
        indexRec = RecoDecay::getMatchedMCRec(particlesMC, arrayDaughters, pdg::Code::kXiCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
        if (indexRec >= 0) {
          // channel = kXic;
          if (indexRec > -1) {
            auto particle = particlesMC.rawIteratorAt(indexRec);
            flag = RecoDecay::getCharmHadronOrigin(particlesMC, particle);

            // look for the correct mass hypothesis
            // either pKπ or πKp
            auto invMass = invMassXicToPKPi;
            if (std::abs(trackFirst.mcParticle().pdgCode()) == kPiPlus) {
              invMass = invMassXicToPiKP;
            }

            if (flag == RecoDecay::OriginType::Prompt) {
              registry.fill(HIST("hMassPromptXicToPKPivsPt"), invMass, pt3Prong);
            } else if (flag == RecoDecay::OriginType::NonPrompt) {
              registry.fill(HIST("hMassNonPromptXicToPKPivsPt"), invMass, pt3Prong);
            }
          }
        }
      }
    } // end loop over 3-prong candidates
  }   // end processMcMatch
};

//________________________________________________________________________________________________________________________
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<HfTrackIndexSkimMcMatchQc>(cfgc));

  return workflow;
}
