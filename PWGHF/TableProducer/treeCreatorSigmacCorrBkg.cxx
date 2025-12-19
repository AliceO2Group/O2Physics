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

/// \file treeCreatorSigmacCorrBkg.cxx
/// \brief Code to reconstruct correlated-background candidates for Σc0,++ analysis
/// \note Λc± candidates selected from the HFLcCandidateSelector.cxx
/// \note Σc0,++ candidates selected from the candidateCreatorSigmac0plusplus.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, INFN PADOVA

#include "PWGHF/Core/DecayChannels.h"
#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/D2H/Utils/utilsSigmac.h"
#include "PWGHF/DataModel/AliasTables.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

#include "Common/Core/RecoDecay.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <Rtypes.h>

#include <cstdint>
#include <cstdlib>

using namespace o2;
using namespace o2::framework; // for Produces, Configuable

namespace o2::aod
{
namespace hf_sigmac_bkg
{
const int pdgCodeLambdac2595 = 14122; // o2-linter: disable=pdg/explicit-code (PDG code needed only for this study)
const int pdgCodeLambdac2625 = 4124;  // o2-linter: disable=pdg/explicit-code (PDG code needed only for this study)
enum Decays { Sigmac2455Pi = 0,
              LambdacPiPi };
enum DecaysLambdac { PKPi = 0,
                     PiKP };
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Mass, mass, float);
DECLARE_SOA_COLUMN(DeltaMass, deltaMass, float);
DECLARE_SOA_COLUMN(Charge, charge, int8_t);
DECLARE_SOA_COLUMN(MotherPdg, motherPdg, int);
DECLARE_SOA_COLUMN(Decay, decay, int8_t);
DECLARE_SOA_COLUMN(DecayLambdac, decayLambdac, int8_t);
DECLARE_SOA_COLUMN(MlScoreFirstClass, mlScoreFirstClass, float); /// background score Λc
DECLARE_SOA_COLUMN(MlScoreThirdClass, mlScoreThirdClass, float); /// non-prompt score Λc
} // namespace hf_sigmac_bkg
DECLARE_SOA_TABLE(HfCorrBkgSc, "AOD", "HFCORRBKGSC",
                  hf_sigmac_bkg::Y,
                  hf_sigmac_bkg::Pt,
                  hf_sigmac_bkg::Mass,
                  hf_sigmac_bkg::DeltaMass,
                  hf_sigmac_bkg::Charge,
                  hf_sigmac_bkg::MotherPdg,
                  hf_sigmac_bkg::Decay,
                  hf_sigmac_bkg::DecayLambdac,
                  hf_sigmac_bkg::MlScoreFirstClass,
                  hf_sigmac_bkg::MlScoreThirdClass);
} // namespace o2::aod

struct HfTreeCreatorSigmacCorrBkg {

  Produces<o2::aod::HfCorrBkgSc> rowCorrBkgSc;

  /// Selection of candidates Λc+
  Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
  Configurable<float> yCandRecoMax{"yCandRecoMax", -1, "Maximum Sc candidate rapidity"};

  using RecoLcMc = soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc, aod::HfMlLcToPKPi>;
  using RecoScMc = soa::Join<aod::HfCandSc, aod::HfCandScMcRec>;
  using ParticlesLcSigmac = soa::Join<aod::McParticles, aod::HfCand3ProngMcGen, aod::HfCandScMcGen>;

  /// @brief init function
  void init(InitContext&) {}

  ///
  void fillTable(RecoScMc::iterator candidateSc, RecoLcMc::iterator candLcDauSc, int motherPdg, int motherDecay = -1)
  {
    const int8_t chargeSc = candidateSc.charge();                                                            // either Σc0 or Σc++
    const float rapidity = chargeSc == 0 ? HfHelper::ySc0(candidateSc) : HfHelper::yScPlusPlus(candidateSc); // NB: since in data we cannot tag Sc(2455) and Sc(2520), then we use only Sc(2455) for y selection on reconstructed signal
    float massSc = -1.f;
    float massLc = -1.f;
    float deltaMass = -1.f;
    const int8_t isCandPKPiPiKP = hf_sigmac_utils::isDecayToPKPiToPiKP(candLcDauSc, candidateSc);
    std::array<float, 2> outputMl{-1., -1.};
    /// rapidity selection on Σc0,++
    if (yCandRecoMax >= 0. && std::abs(rapidity) > yCandRecoMax) {
      return;
    }

    /// BDT scores
    if (!candLcDauSc.mlProbLcToPiKP().empty()) {
      outputMl.at(0) = candLcDauSc.mlProbLcToPiKP()[0]; /// bkg score
      outputMl.at(1) = candLcDauSc.mlProbLcToPiKP()[2]; /// non-prompt score
    }

    if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PKPi)) && std::abs(candLcDauSc.template prong0_as<aod::TracksWMc>().template mcParticle_as<ParticlesLcSigmac>().pdgCode()) == kProton) {
      massSc = HfHelper::invMassScRecoLcToPKPi(candidateSc, candLcDauSc);
      massLc = HfHelper::invMassLcToPKPi(candLcDauSc);
      deltaMass = massSc - massLc;

      /// fill the tree
      rowCorrBkgSc(rapidity, candidateSc.pt(), massSc, deltaMass, chargeSc, motherPdg, motherDecay, aod::hf_sigmac_bkg::DecaysLambdac::PKPi, outputMl.at(0), outputMl.at(1));
    }
    if ((TESTBIT(isCandPKPiPiKP, o2::aod::hf_cand_sigmac::Decays::PiKP)) && std::abs(candLcDauSc.template prong0_as<aod::TracksWMc>().template mcParticle_as<ParticlesLcSigmac>().pdgCode()) == kPiPlus) {
      massSc = HfHelper::invMassScRecoLcToPiKP(candidateSc, candLcDauSc);
      massLc = HfHelper::invMassLcToPiKP(candLcDauSc);
      deltaMass = massSc - massLc;

      /// fill the tree
      rowCorrBkgSc(rapidity, candidateSc.pt(), massSc, deltaMass, chargeSc, motherPdg, motherDecay, aod::hf_sigmac_bkg::DecaysLambdac::PiKP, outputMl.at(0), outputMl.at(1));
    }
  }

  /// @brief process function to loop over the Σc reconstructed candidates and match them to corr. background sources in MC
  void process(RecoScMc const& candidatesSc,
               ParticlesLcSigmac const& particles,
               RecoLcMc const&,
               aod::TracksWMc const&)
  {
    /// loop over reconstructed Σc candidates
    for (auto const& candidateSc : candidatesSc) {

      auto candLcDauSc = candidateSc.template prongLc_as<RecoLcMc>();
      auto candSoftPiDauSc = candidateSc.template prong1_as<aod::TracksWMc>();

      /// tag immediately the Σc0,++(2455) and Σc0,++(2520) signal
      auto flagMcDecayChanScAbs = std::abs(candidateSc.flagMcMatchRec());
      bool const isTrueSigmac0 = (flagMcDecayChanScAbs == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::Sc0ToPKPiPi);
      bool const isTrueSigmacPlusPlus = (flagMcDecayChanScAbs == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScplusplusToPKPiPi);
      bool const isTrueSigmacStar0 = (flagMcDecayChanScAbs == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStar0ToPKPiPi);
      bool const isTrueSigmacStarPlusPlus = (flagMcDecayChanScAbs == o2::hf_decay::hf_cand_sigmac::DecayChannelMain::ScStarPlusPlusToPKPiPi);
      if (isTrueSigmac0) {
        /// fill the output for the signal
        fillTable(candidateSc, candLcDauSc, o2::constants::physics::Pdg::kSigmaC0);

        /// the candidate that we reconstructed is a real Sigmac(2455, 2520), but later we look for correlated background sources
        /// let's continue
        continue;
      }
      if (isTrueSigmacPlusPlus) {
        /// fill the output for the signal
        fillTable(candidateSc, candLcDauSc, o2::constants::physics::Pdg::kSigmaCPlusPlus);

        /// the candidate that we reconstructed is a real Sigmac(2455, 2520), but later we look for correlated background sources
        /// let's continue
        continue;
      }
      if (isTrueSigmacStar0) {
        /// fill the output for the signal
        fillTable(candidateSc, candLcDauSc, o2::constants::physics::Pdg::kSigmaCStar0);

        /// the candidate that we reconstructed is a real Sigmac(2455, 2520), but later we look for correlated background sources
        /// let's continue
        continue;
      }
      if (isTrueSigmacStarPlusPlus) {
        /// fill the output for the signal
        fillTable(candidateSc, candLcDauSc, o2::constants::physics::Pdg::kSigmaCStarPlusPlus);

        /// the candidate that we reconstructed is a real Sigmac(2455, 2520), but later we look for correlated background sources
        /// let's continue
        continue;
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///                                                                                                                     ///
      ///   Look for the Σc0,++ correlated background sources.                                                                ///
      ///                                                                                                                     ///
      ///   Two sources are possible:                                                                                         ///
      ///                                                                                                                     ///
      ///     1) Λc±(2595, 2625) → Σc0,++(2455) π+,-                                                                          ///
      ///        In this case, we need that the candidate Σc0,++(2455) is formed by the Λc± daughter of a real Σc0,++(2455)   ///
      ///        paired with the bachelor pion of the Λc±(2595, 2625)                                                         ///
      ///                                                                                                                     ///
      ///     2) Λc±(2595, 2625) → Λc± π+ π-                                                                                  ///
      ///        It means that the reconstructed Σc candidate it's actually the pair of Λc± π+ or Λc± π-                      ///
      ///        coming from the same Λc±(2595) or Λc±(2625)                                                                  ///
      ///                                                                                                                     ///
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      /// check if the candidate Lc and soft pion daugthers are not real Lc or pion
      bool const isLambdac = std::abs(candLcDauSc.flagMcMatchRec()) == hf_decay::hf_cand_3prong::DecayChannelMain::LcToPKPi;
      bool isPion = false;
      if (candSoftPiDauSc.has_mcParticle()) {
        isPion = std::abs(candSoftPiDauSc.template mcParticle_as<ParticlesLcSigmac>().pdgCode()) == kPiPlus;
      }
      if (!isLambdac || !isPion) {
        continue;
      }

      /// First of all, search for an incomplete decay of Λc±(2595, 2625) into Λc± (→ pK-π+) π+ π-,
      /// we will distinguish the two background sources from the first mother of the Λc±
      /// This means that we look for cases in which we use 4 tracks out of 5, missing one pion.
      /// The possible combinations are: (1) pK-π+ π+; (2) pK-π+ π-
      /// The full chain has a max depth of 4:
      ///    i.  Λc±(2595, 2625) → Σc0,++(2455) π+,- (+1)
      ///    ii. Σc0,++(2455) → Λc+ π-,+ (+1)
      ///    ii. Λc+ → pK-π+ which can go in the direct channel direct (+1) or through a resonance (+2)
      auto arrayDaughters = std::array{candLcDauSc.template prong0_as<aod::TracksWMc>(),
                                       candLcDauSc.template prong1_as<aod::TracksWMc>(),
                                       candLcDauSc.template prong2_as<aod::TracksWMc>(),
                                       candidateSc.template prong1_as<aod::TracksWMc>()};
      auto arrayPdgDaughters1 = std::array{+kProton, -kKPlus, +kPiPlus, +kPiPlus}; /// (1)
      auto arrayPdgDaughters2 = std::array{+kProton, -kKPlus, +kPiPlus, -kPiPlus}; /// (2)
      int8_t sign = 0;
      int indexMother = -1;
      int motherPdg = -1;
      int motherDecay = -1;
      // look for Λc±(2595) - first daughter pdg-code combination
      indexMother = RecoDecay::getMatchedMCRec<false, false, true /*acceptIncompleteReco*/, true /*acceptTrackDecay*/, true /*acceptTrackIntWithMaterial*/>(particles, arrayDaughters, aod::hf_sigmac_bkg::pdgCodeLambdac2595, arrayPdgDaughters1, true, &sign, 4 /*depthMainMax*/);
      if (indexMother >= 0) {
        /// mother found!
        motherPdg = aod::hf_sigmac_bkg::pdgCodeLambdac2595;
      } else {
        // look for Λc±(2595) - second daughter pdg-code combination
        indexMother = RecoDecay::getMatchedMCRec<false, false, true /*acceptIncompleteReco*/, true /*acceptTrackDecay*/, true /*acceptTrackIntWithMaterial*/>(particles, arrayDaughters, aod::hf_sigmac_bkg::pdgCodeLambdac2595, arrayPdgDaughters2, true, &sign, 4 /*depthMainMax*/);
        if (indexMother >= 0) {
          /// mother found!
          motherPdg = aod::hf_sigmac_bkg::pdgCodeLambdac2595;
        } else {
          // look for Λc±(2625) - first daughter pdg-code combination
          indexMother = RecoDecay::getMatchedMCRec<false, false, true /*acceptIncompleteReco*/, true /*acceptTrackDecay*/, true /*acceptTrackIntWithMaterial*/>(particles, arrayDaughters, aod::hf_sigmac_bkg::pdgCodeLambdac2625, arrayPdgDaughters1, true, &sign, 4 /*depthMainMax*/);
          if (indexMother >= 0) {
            /// mother found!
            motherPdg = aod::hf_sigmac_bkg::pdgCodeLambdac2625;
          } else {
            // look for Λc±(2625) - second daughter pdg-code combination
            indexMother = RecoDecay::getMatchedMCRec<false, false, true /*acceptIncompleteReco*/, true /*acceptTrackDecay*/, true /*acceptTrackIntWithMaterial*/>(particles, arrayDaughters, aod::hf_sigmac_bkg::pdgCodeLambdac2625, arrayPdgDaughters2, true, &sign, 4 /*depthMainMax*/);
            if (indexMother >= 0) {
              /// mother found!
              motherPdg = aod::hf_sigmac_bkg::pdgCodeLambdac2625;
            } else {
              /// no mother found, it means that this is not a corr. bkg candidate
              /// let's skip it
              continue;
            }
          }
        }
      }
      // LOG(info) << "motherPdg: " << motherPdg;

      /// now that we matched a Λc±(2595, 2625), let's determine the precise decay channel
      /// by checking the mother of the Λc±
      auto arrayDaughtersLambdac = std::array{candLcDauSc.template prong0_as<aod::TracksWMc>(),
                                              candLcDauSc.template prong1_as<aod::TracksWMc>(),
                                              candLcDauSc.template prong2_as<aod::TracksWMc>()};
      int8_t signLambdac = 0;
      int const indexRecLc = RecoDecay::getMatchedMCRec<false, false, false, true, true>(particles, arrayDaughtersLambdac, o2::constants::physics::Pdg::kLambdaCPlus, std::array{+kProton, -kKPlus, +kPiPlus}, true, &signLambdac, 2);
      if (indexRecLc < 0) {
        /// this should never happen, since we check above that the isLambdac==true
        LOG(fatal) << "Generated Lambdac not found. Not expected. Aborting.";
      }
      auto particleLc = particles.rawIteratorAt(indexRecLc);
      if (particleLc.has_mothers()) {
        /// we should always enter here, since the Λc± is coming from a Λc±(2595, 2625) decay
        for (auto iMother = particleLc.mothersIds().front(); iMother <= particleLc.mothersIds().back(); ++iMother) {
          auto mother = particles.rawIteratorAt(iMother);
          int const pdgCodeMotherAbs = std::abs(mother.pdgCode());
          if (pdgCodeMotherAbs == o2::constants::physics::Pdg::kSigmaC0 || pdgCodeMotherAbs == o2::constants::physics::Pdg::kSigmaCPlusPlus) {
            /// the Λc± comes from a Σc0,++(2455)
            /// ==> we found a Λc±(2595, 2625) → Σc0,++(2455) π+,- decay!
            motherDecay = aod::hf_sigmac_bkg::Decays::Sigmac2455Pi;

            /// This should be enough, i.e. not necessary to check that the pion is not daughter of the same Sigmac
            /// This case should be already excluded by searching for real Sigmac (see the beginning of the process function)

          } else {
            /// considering all the checks done before, the only other possibility is that this Λc± directly comes from a Λc±(2595, 2625)
            motherDecay = aod::hf_sigmac_bkg::Decays::LambdacPiPi;
          }
          break;
        }
      } else {
        /// we should never eneter here
        LOG(fatal) << "Lambdac particle without mothers. Not expected. Aborting.";
      }

      /// we found a corr. bkg. candidate
      /// let's fill our output
      fillTable(candidateSc, candLcDauSc, motherPdg, motherDecay);

    } /// end loop over reconstructed Σc candidates
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTreeCreatorSigmacCorrBkg>(cfgc)};
}
