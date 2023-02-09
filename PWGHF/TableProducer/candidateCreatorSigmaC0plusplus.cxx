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

/// \file candidateCreatorSigmaC0plusplus.cxx
/// \brief Σc0,++ → Λc+(→pK-π+) π-,+ candidate builder
/// \note Λc± candidates selected from the HFLcCandidateSelector.cxx
///
/// \author Mattia Faggin <mfaggin@cern.ch>, University and INFN PADOVA

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_3prong;

struct candidateCreatorSigmaC0plusplus {

    /// Table with Σc0,++ info
    Produces<aod::HfCandSigmaC> rowScCandidates;

    /// Selection of candidates Λc+
    Configurable<int> selectionFlagLc{"selectionFlagLc", 1, "Selection Flag for Lc"};
    Configurable<double> cutYCandLcMax{"cutYCandLcMax", -1., "max. candLc. Lc rapidity"};
    Configurable<double> cutMinvpKpiCandLcMax{"cutMinvpKpiCandLcMax", 0.03, "max. spread (abs. value9) between PDG(Lc) and Minv(pKpi)"};
    Configurable<double> cutMinvpiKpCandLcMax{"cutMinvpiKpCandLcMax", 0.03, "max. spread (abs. value9) between PDG(Lc) and Minv(pK0s)"};

    /// Selections on candidate soft π-,+
    Configurable<float> softPiEta{"softPiEta", 0.9f, "Soft pion max value for pseudorapidity (abs vale)"};
    Configurable<int> softPiITSHitMap{"softPiITSHitMap", 127, "Soft pion ITS hitmap"};
    Configurable<int> softPiITSHitsMin{"softPiITSHitsMin", 1, "Minimum number of ITS layers crossed by the soft pion among those in \"softPiITSHitMap\""};
    Configurable<float> softPidcaXYmax{"softPidcaXYmax", 0.065, "Soft pion max dcaXY (cm)"};
    Configurable<float> softPidcaZmax{"softPidcaZmax", 0.065, "Soft pion max dcaZ (cm)"};

    /// Filter the candidate Λc+ used for the Σc0,++ creation
    Filter filterSelectCandidateLc = (aod::hf_sel_candidate_lc::isSelLcToPKPi >= selectionFlagLc || aod::hf_sel_candidate_lc::isSelLcToPiKP >= selectionFlagLc);

    /// Cut selection object for soft π-,+
    TrackSelection softPiCuts;

    /// @brief init function, to define the soft pion selections and histograms
    /// @param  
    void init(InitContext&) {

        ////////////////////////////////////////
        /// set the selections for soft pion ///
        ////////////////////////////////////////
        softPiCuts.SetEtaRange(-softPiEta, softPiEta);  // eta
        softPiCuts.SetMaxDcaXY(softPidcaXYmax); // dcaXY
        softPiCuts.SetMaxDcaZ(softPidcaZmax);   // dcaZ
        // ITS hitmap
        std::set<uint8_t> set_softPiITSHitMap; // = {};
        for (int ITSlayerId = 0; ITSlayerId < 7; ITSlayerId++) {
            if ((softPiITSHitMap & (1 << ITSlayerId)) > 0) {
              set_softPiITSHitMap.insert(static_cast<uint8_t>(ITSlayerId));
            }
        }
        LOG(info) << "### ITS hitmap for soft pion";
        LOG(info) << "    >>> set_softPiITSHitMap.size(): " << set_softPiITSHitMap.size();
        LOG(info) << "    >>> Custom ITS hitmap checked: ";
        for (std::set<uint8_t>::iterator it = set_softPiITSHitMap.begin(); it != set_softPiITSHitMap.end(); it++) {
          LOG(info) << "        Layer " << (int)(*it) << " ";
        }
        LOG(info) << "############";
        softPiCuts.SetRequireITSRefit();
        softPiCuts.SetRequireHitsInITSLayers(softPiITSHitsMin, set_softPiITSHitMap);

    }

    /// @brief process function for Σc0,++ → Λc+(→pK-π+) π- candidate reconstruction
    /// @param collision is a o2::aod::Collision
    /// @param tracks are the tracks (with dcaXY, dcaZ information) in the collision → soft-pion candidate tracks
    /// @param candidates are 3-prong candidates satisfying the analysis selections for Λc+ → pK-π+ (and charge conj.)
    void process(const o2::aod::Collision& collision, const soa::Join<aod::Tracks, aod::TracksDCA>& tracks, soa::Filtered<soa::Join<aod::HfCand3Prong, aod::HfSelLc>> const& candidates) {

      /// loop over Λc+ → pK-π+ (and charge conj.) candidates
      for(auto& candLc : candidates) {

        /// keep only the candidates flagged as possible Λc+ (and charge conj.) decaying into a charged pion, kaon and proton
        /// if not selected, skip it and go to the next one
        if (!(candLc.hfflag() & 1 << DecayType::LcToPKPi)) {
          continue;
        }
        /// keep only the candidates Λc+ (and charge conj.) within the desired rapidity
        /// if not selected, skip it and go to the next one
        if (cutYCandLcMax >= 0. && std::abs(yLc(candLc)) > cutYCandLcMax) {
          continue;
        }

        /// selection on the Λc+ inv. mass window we want to consider for Σc0,++ candidate creation
        auto statusSpreadMinvpKpiFromPDG = 0;
        auto statusSpreadMinvpiKpFromPDG = 0;
        if(candLc.isSelLcToPKPi() >= 1 && std::abs( invMassLcToPKPi(candLc) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus) ) <= cutMinvpKpiCandLcMax) {
          statusSpreadMinvpKpiFromPDG = 1;
        }
        if(candLc.isSelLcToPiKP() >= 1 && std::abs( invMassLcToPiKP(candLc) - RecoDecay::getMassPDG(pdg::Code::kLambdaCPlus) ) <= cutMinvpiKpCandLcMax) {
          statusSpreadMinvpiKpFromPDG = 1;
        }
        if (statusSpreadMinvpKpiFromPDG == 0 && statusSpreadMinvpiKpFromPDG == 0){
          /// none of the two possibilities are satisfied, therefore this candidate Lc can be skipped
          continue;
        }

        /// loop over tracks
        for(auto& trackSoftPi : tracks) {

          /////////////////////////////////////////////////////////////////////////////////
          ///                       Σc0,++ candidate creation                           ///
          ///                                                                           ///
          /// For each candidate Λc, let's loop over all the candidate soft-pion tracks ///
          /////////////////////////////////////////////////////////////////////////////////

          /// keep only soft-pion candidate tracks
          /// if not selected, skip it and go to the next one
          if( !softPiCuts.IsSelected(trackSoftPi) ) {
            continue;
          }

          /// Exclude the current candidate soft pion if it corresponds already to a candidate Lc prong
          int indexProng0 = candLc.prong0_as<aod::Tracks>().globalIndex();
          int indexProng1 = candLc.prong1_as<aod::Tracks>().globalIndex();
          int indexProng2 = candLc.prong2_as<aod::Tracks>().globalIndex();
          int indexSoftPi = trackSoftPi.globalIndex();
          if (indexSoftPi == indexProng0 || indexSoftPi == indexProng1 || indexSoftPi == indexProng2)
            continue;

          /// determine the Σc candidate charge
          int chargeLc = candLc.prong0_as<aod::Tracks>().sign() + candLc.prong1_as<aod::Tracks>().sign() + candLc.prong2_as<aod::Tracks>().sign();
          int chargeSoftPi = trackSoftPi.sign();
          int chargeSc = chargeLc + chargeSoftPi;
          if(std::abs(chargeSc) != 0 && std::abs(chargeSc) != 2) {
            /// this shall never happen
            LOG(fatal) << ">>> SigmaC candidate with charge +1 built, not possible! Charge Lc: " << chargeLc << ", charge soft pion: " << chargeSoftPi;
            continue;
          }

          /// fill the Σc0,++ candidate table
          rowScCandidates(/* general columns */
                          candLc.collisionId(),
                          /* 2-prong specific columns */
                          candLc.px(), candLc.py(), candLc.pz(),
                          trackSoftPi.px(), trackSoftPi.py(), trackSoftPi.pz(),
                          candLc.globalIndex(), trackSoftPi.globalIndex(),
                          candLc.hfflag(),
                          /* Σc0,++ specific columns */
                          chargeSc,
                          statusSpreadMinvpKpiFromPDG, statusSpreadMinvpiKpFromPDG);

        } /// end loop over tracks

      } /// end loop over candidtes

    } /// end process
};

struct candidateSigmaC0plusplusMcMatch {

    /// Table with MC info for reconstructed and geerated Σc0,++
    /// TODO
    /// [...]  

    /// @brief init function
    void init(InitContext const&) {}

    using LambdacMC = soa::Join<aod::HfCand3Prong, aod::HfSelLc, aod::HfCand3ProngMcRec>;
    using TracksMC = soa::Join<aod::Tracks, aod::McTrackLabels>;

    /// @brief process function for MC matching of Σc0,++ → Λc+(→pK-π+) π- reconstructed candidates and counting of generated ones
    /// @param candidatesSigmaC reconstructed Σc0,++ candidates
    /// @param particlesMc table of generaed particles
    void processMC(const aod::HfCandSigmaC& candidatesSigmaC, aod::McParticles const& particlesMc,
    LambdacMC const&, const TracksMC&) {

      Produces<aod::HfCandSigmaCMcRec> rowMCMatchSigmaCRec;
      Produces<aod::HfCandSigmaCMcGen> rowMCMatchSigmaCGen;

      int indexRec = -1;
      int8_t sign = 0;
      int8_t flag = 0;
      int8_t origin = 0;
      int chargeSc = 999;
      //std::vector<int> arrDaughIndex; /// index of daughters of MC particle

      /// Match reconstructed Σc0,++ candidates
      for(auto const& candSigmaC : candidatesSigmaC) {
        indexRec = -1;
        sign = 0;
        flag = 0;
        origin = 0;
        //arrDaughIndex.clear();

        /// skip immediately the candidate Σc0,++ w/o a Λc+ matched to MC
        auto candLc = candSigmaC.prong0_as<LambdacMC>();
        if (!(std::abs(candLc.flagMcMatchRec()) == 1 << DecayType::LcToPKPi)) {
          rowMCMatchSigmaCRec(flag, origin);
          continue;
        }

        /// matching to MC
        auto arrayDaughters = array{candLc.prong0_as<TracksMC>(),
                                    candLc.prong1_as<TracksMC>(),
                                    candLc.prong2_as<TracksMC>(),
                                    candSigmaC.prong1_as<TracksMC>()};
        chargeSc = candSigmaC.charge();
        if(chargeSc == 0) {
          /// candidate Σc0
          indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughters, pdg::Code::kSigmaC0, array{(int) kProton, (int) kKMinus, (int) kPiPlus, (int) kPiMinus}, true, &sign, 2);
          if (indexRec > -1) {
            flag = sign * (1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi);
          }
        } else if(std::abs(chargeSc) == 2) {
          /// candidate Σc++
          indexRec = RecoDecay::getMatchedMCRec(particlesMc, arrayDaughters, pdg::Code::kSigmaCPlusPlus, array{(int) kProton, (int) kKMinus, (int) kPiPlus, (int) kPiPlus}, true, &sign, 2);
          if (indexRec > -1) {
            flag = sign * (1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi);
          }
        }

        /// check the origin (prompt vs. non-prompt) 
        if(flag != 0) {
          auto particle = particlesMc.rawIteratorAt(indexRec);
          origin = RecoDecay::getCharmHadronOrigin(particlesMc, particle);
        }

        /// fill the table with results of reconstruction level MC matching
        rowMCMatchSigmaCRec(flag, origin);
      } /// end loop over reconstructed Σc0,++ candidates

      /// Match generated Σc0,++ candidates
      for(auto& particle : particlesMc) {
        flag = 0;
        origin = 0;

        if (RecoDecay::isMatchedMCGen(particlesMc, particle, pdg::Code::kSigmaC0, array{(int) kProton, (int) kKMinus, (int) kPiPlus, (int) kPiMinus}, true, &sign, 2)) {
          // generated Σc0
          flag = sign * (1 << aod::hf_cand_sc::DecayType::SigmaC0ToPKPiPi);
        } else if (RecoDecay::isMatchedMCGen(particlesMc, particle, pdg::Code::kSigmaCPlusPlus, array{(int) kProton, (int) kKMinus, (int) kPiPlus, (int) kPiPlus}, true, &sign, 2)) {
          // generated Σc++
          flag = sign * (1 << aod::hf_cand_sc::DecayType::SigmaCplusplusToPKPiPi);
        }

        /// fill the table with results of generation level MC matching
        rowMCMatchSigmaCGen(flag, origin);

      }/// end loop over particlesMc  
    } /// end processMC
    PROCESS_SWITCH(candidateSigmaC0plusplusMcMatch, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<candidateCreatorSigmaC0plusplus>(cfgc, TaskName{"hf-sigmac-candidate-creator"}),
    adaptAnalysisTask<candidateSigmaC0plusplusMcMatch>(cfgc, TaskName{"hf-sigmac-mc-match"})};
}