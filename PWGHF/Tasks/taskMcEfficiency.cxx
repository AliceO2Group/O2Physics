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

/// \file taskMcEfficiency.cxx
/// \brief Task to analyse the MC efficiency ("found / trackable")
///
/// \author Jan Fiete Grosse-Oetringhaus, CERN

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
//using namespace o2::aod::hf_cand_2prong;
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

struct HfTaskMcEfficiency {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  ConfigurableAxis axisPt{"axisPt", {10, 0, 10}, "pT axis"};
  ConfigurableAxis axisMass{"axisMass", {120, 1.5848, 2.1848}, "m_inv axis"};
  ConfigurableAxis axisPdg{"axisPdg", {VARIABLE_WIDTH, -421.5, 0, 421.5}, "PDG code axis"};
  ConfigurableAxis axisCPA{"axisCPA", {102, -1.02, 1.02}, "Cosine of pointing angle axis"};
  Configurable<std::vector<int>> pdgCodes{"pdgCodes", {pdg::kD0Bar, pdg::kD0}, "PDG codes to process"};

  Configurable<float> mcAcceptancePt{"mcAcceptancePt", 0.1, "MC Acceptance: lower pt limit"};
  Configurable<float> mcAcceptanceEta{"mcAcceptanceEta", 0.8, "MC Acceptance: upper eta limit"};

  enum HFStep { kHFStepMC = 0,
                kHFStepMcInRapidity,        // MC mothers in rapidity |y| < 0.5
                kHFStepAcceptance,          // MC mothers where all prongs are in the acceptance
                kHFStepTrackable,           // MC mothers where all prongs have a reconstructed track
                kHFStepAcceptanceTrackable, // MC mothers where all prongs are in the acceptance and have a reconstructed track
                kHFStepTrackableCuts,       // MC mothers where all prongs have a reconstructed track which passes the track selection
                kHFStepTracked,             // signal candidates which have been found
                kHFStepTrackedCuts,         // signal candidates which have been found and all prongs pass the track selection (cross-check)
                kHFStepTrackedSelected,     // signal candidates which pass the selector
                kHFStepTrackedDuplicates,   // signal candidates which pass the selector and are duplicated (only the second, third etc. ones are filled)
                kHFNSteps };

  enum TrackableStep { kTrackableAll = 0, // all tracks
                       kTrackableITS,     // track contains ITS information
                       kTrackableTPC,     // track contains TPC information
                       kTrackableITSTPC,  // track contains ITS and TPC information
                       kNTrackableSteps };

  using TracksWithSelection = soa::Join<aod::Tracks, aod::TrackSelection>;
  using TracksWithSelectionMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TrackSelection>;

  HistogramRegistry registry{"registry"};

  bool isThereD0 = (std::find(pdgCodes.value.begin(), pdgCodes.value.end(), pdg::kD0) != pdgCodes.value.end());
  bool isThereD0Bar = (std::find(pdgCodes.value.begin(), pdgCodes.value.end(), pdg::kD0Bar) != pdgCodes.value.end());
  bool isThereLc = (std::find(pdgCodes.value.begin(), pdgCodes.value.end(), pdg::kLambdaCPlus) != pdgCodes.value.end());

  void init(o2::framework::InitContext&)
  {
    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPt, axisMass, axisPdg, axisCPA, {2, -0.5, 1.5, "collision matched"}}, kHFNSteps});
    hCandidates->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    hCandidates->GetAxis(1)->SetTitle("m_{inv}");
    hCandidates->GetAxis(2)->SetTitle("PDG code");
    hCandidates->GetAxis(3)->SetTitle("CPA");

    registry.add("hDuplicateCount", "Duplicate count;frequency;count", {HistType::kTH1F, {{10, 0.5, 10.5}}});

    auto hTrackablePtEta = registry.add<StepTHn>("hTrackablePtEta", "Prongs kinematics at different steps", {HistType::kStepTHnF, {{20, 0, 10}, {40, -2, 2}}, kNTrackableSteps});
    hTrackablePtEta->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    hTrackablePtEta->GetAxis(1)->SetTitle("#eta");
  }

  template <typename T>
  inline bool checkTrack(T track)
  {
    // TODO configurable?
    return track.isGlobalTrackWoDCA();
  }

  template <bool mc, typename T1, typename T2, typename T3>
  void candidate3ProngLoop(T1& candidates, T2& tracks, T3& mcParticles)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<std::size_t, int> duplicates;

    for (const auto& candidate : candidates) { /// loop over candidates

      for (const auto pdgCode : pdgCodes.value) { /// loop on pdg codes
        auto decayType = -1;
        std::array<int, 3> pdgDaughters;

        if (pdgCode == pdg::kLambdaCPlus) {
          decayType = 1 << o2::aod::hf_cand_3prong::DecayType::LcToPKPi;
          pdgDaughters[0] = +kProton;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else {
          LOGP(fatal, "Not implemented for PDG {}", pdgCode);
        }

        if (!(candidate.hfflag() & decayType)) {
          continue;
        }

        auto trackPos = candidate.template prong0_as<TracksType>();
        auto trackNeg = candidate.template prong1_as<TracksType>();
        auto trackThird = candidate.template prong2_as<TracksType>();

        bool isHypoMass1TrackStep = true;
        bool isHypoMass2TrackStep = true;
        bool isHypoMass1SelStep = false;
        bool isHypoMass2SelStep = false;
        if (pdgCode == pdg::kLambdaCPlus) {
          isHypoMass1SelStep = candidate.isSelLcToPKPi();  /// from candidate selector!
          isHypoMass2SelStep = candidate.isSelLcToPiKP();  /// from candidate selector!
        }
        bool collisionMatched = false;
        if constexpr (mc) { /// info MC used
          int8_t sign = 0;
          int indexRec = -999;
          if (pdgCode == pdg::kLambdaCPlus) {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg, trackThird}, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
          }
          
          if (indexRec < 0) {
            continue;
          }

          if (pdgCode == pdg::kLambdaCPlus) {
            auto daughter = trackPos.mcParticle();
            if (std::abs( daughter.pdgCode() ) == kProton ) {
              isHypoMass1TrackStep = true;
              isHypoMass1SelStep = true;
              isHypoMass2TrackStep = false;
              isHypoMass2SelStep = false;
            } else if (std::abs( daughter.pdgCode() ) == kPiPlus) {
              isHypoMass1TrackStep = false;
              isHypoMass1SelStep = false;
              isHypoMass2TrackStep = true;
              isHypoMass2SelStep = true;
            } else {
              continue;
            }
          }

          collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRec).mcCollisionId();
        } /// end info MC used

        float massHypo1 = -1;
        float massHypo2 = -1;
        float cpa = candidate.cpa();
        float pt = candidate.pt();
        //bool selected = false;

        /// all candidates
        if (isHypoMass1TrackStep) {
          if (pdgCode == pdg::kLambdaCPlus) {
            massHypo1 = invMassLcToPKPi(candidate);
          }
          hCandidates->Fill(kHFStepTracked, pt, massHypo1, pdgCode, cpa, collisionMatched);
        }
        if (isHypoMass2TrackStep) {
          if (pdgCode == pdg::kLambdaCPlus) {
            massHypo2 = invMassLcToPiKP(candidate);
          }
          hCandidates->Fill(kHFStepTracked, pt, massHypo2, pdgCode, cpa, collisionMatched);
        }

        // check if prongs have passed track cuts
        if (checkTrack(trackPos) && checkTrack(trackNeg) && checkTrack(trackThird)) {
          if (isHypoMass1TrackStep) {
            hCandidates->Fill(kHFStepTrackedCuts, pt, massHypo1, pdgCode, cpa, collisionMatched);
          }
          if (isHypoMass2TrackStep) {
            hCandidates->Fill(kHFStepTrackedCuts, pt, massHypo2, pdgCode, cpa, collisionMatched);
          }
        }

        if (!isHypoMass1SelStep && !isHypoMass2SelStep) {
          continue;
        }

        // selected candidates
        if (isHypoMass1SelStep) {
          hCandidates->Fill(kHFStepTrackedSelected, pt, massHypo1, pdgCode, cpa, collisionMatched);
        }
        if (isHypoMass2SelStep) {
          hCandidates->Fill(kHFStepTrackedSelected, pt, massHypo2, pdgCode, cpa, collisionMatched);
        }

        // duplicates
        std::array<int, 3> prongIds = {candidate.prong0Id(), candidate.prong1Id(), candidate.prong2Id()};
        std::sort(prongIds.begin(), prongIds.end());
        std::string concat = std::to_string(prongIds[0]) + std::to_string(prongIds[1]) + std::to_string(prongIds[2]);
        std::size_t hash = std::hash<std::string>{}(concat); /// unique value for the 'concat' string
        if (duplicates.find(hash) != duplicates.end()) {
          if (isHypoMass1TrackStep) {
            hCandidates->Fill(kHFStepTrackedDuplicates, pt, massHypo1, pdgCode, cpa, collisionMatched);
          }
          if (isHypoMass2TrackStep) {
            hCandidates->Fill(kHFStepTrackedDuplicates, pt, massHypo2, pdgCode, cpa, collisionMatched);
          }
        }
        duplicates[hash]++;

      } /// end loop on pdg codes
    } /// end loop over candidates

    auto hDuplicateCount = registry.get<TH1>(HIST("hDuplicateCount"));
    for (const auto& i : duplicates) {
      hDuplicateCount->Fill(i.second);
    }
  }

  template <bool mc, typename T1, typename T2, typename T3>
  void candidate2ProngLoop(T1& candidates, T2& tracks, T3& mcParticles)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<int64_t, int> duplicates;

    for (const auto pdgCode : pdgCodes.value) {
      auto decayType = -1;
      std::array<int, 2> pdgDaughters;

      if (pdgCode == pdg::kD0) {
        decayType = 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK;
        pdgDaughters[0] = +kPiPlus;
        pdgDaughters[1] = -kKPlus;
      } else if (pdgCode == pdg::kD0Bar) {
        decayType = 1 << o2::aod::hf_cand_2prong::DecayType::D0ToPiK;
        pdgDaughters[0] = -kPiPlus;
        pdgDaughters[1] = +kKPlus;
      } else {
        LOGP(fatal, "Not implemented for PDG {}", pdgCode);
      }

      int nTracked = 0;
      int nSelected = 0;
      for (const auto& candidate : candidates) {
        if (!(candidate.hfflag() & decayType)) {
          continue;
        }

        auto trackPos = candidate.template prong0_as<TracksType>();
        auto trackNeg = candidate.template prong1_as<TracksType>();

        bool collisionMatched = false;
        if constexpr (mc) {
          auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg}, pdgCode, pdgDaughters, false);
          if (indexRec < 0) {
            continue;
          }

          collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRec).mcCollisionId();
        }

        float mass = -1;
        float cpa = candidate.cpa();
        float pt = candidate.pt();
        bool selected = false;
        if (pdgCode == pdg::kD0) {
          mass = invMassD0ToPiK(candidate);
          selected = candidate.isSelD0() >= selectionFlagD0;
        } else if (pdgCode == pdg::kD0Bar) {
          mass = invMassD0barToKPi(candidate);
          selected = candidate.isSelD0bar() >= selectionFlagD0bar;
        }
        LOGP(debug, "Candidate {} has prong {} and prong {} and pT {} and mass {}", candidate.globalIndex(), candidate.prong0Id(), candidate.prong1Id(), candidate.pt(), mass);

        // all candidates
        hCandidates->Fill(kHFStepTracked, pt, mass, pdgCode, cpa, collisionMatched);
        ++nTracked;

        // check if prongs have passed track cuts
        if (checkTrack(trackPos) && checkTrack(trackNeg)) {
          hCandidates->Fill(kHFStepTrackedCuts, pt, mass, pdgCode, cpa, collisionMatched);
        }

        if (!selected) {
          continue;
        }

        // selected candidates
        hCandidates->Fill(kHFStepTrackedSelected, pt, mass, pdgCode, cpa, collisionMatched);
        ++nSelected;

        // duplicates
        /// put two 32-bit indices in a 64-bit integer
        int64_t hash = 0;
        if (candidate.prong0Id() < candidate.prong1Id()) {
          hash = ((int64_t)candidate.prong0Id() << 32) | candidate.prong1Id();
        } else {
          hash = ((int64_t)candidate.prong1Id() << 32) | candidate.prong0Id();
        }
        if (duplicates.find(hash) != duplicates.end()) {
          hCandidates->Fill(kHFStepTrackedDuplicates, pt, mass, pdgCode, cpa, collisionMatched);
        }
        duplicates[hash]++;
      }

      LOGP(debug, "PDG code: {}: Tracked: {} | Selected: {}", pdgCode, nTracked, nSelected);
    }

    auto hDuplicateCount = registry.get<TH1>(HIST("hDuplicateCount"));
    for (const auto& i : duplicates) {
      hDuplicateCount->Fill(i.second);
    }
  }

  void processDataD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelection& tracks)
  {
    candidate2ProngLoop<false>(candidates, tracks, tracks); // NOTE third argument has to be provided but is not used as template argument is <false>
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataD0, "Process D0 data (no MC information needed)", false);

  template <typename C>
  void candidate2ProngMcLoop(C& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls)
  {
    candidate2ProngLoop<true>(candidates, tracks, mcParticles);

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    auto hTrackablePtEta = registry.get<StepTHn>(HIST("hTrackablePtEta"));

    // lists for optimization
    std::vector<bool> tracked(mcParticles.size(), false);
    std::vector<bool> hasITS(mcParticles.size(), false);
    std::vector<bool> hasTPC(mcParticles.size(), false);
    std::vector<bool> selected(mcParticles.size(), false);
    for (const auto& track : tracks) {
      if (track.mcParticleId() >= 0) {
        tracked[track.mcParticleId()] = true;
        if (checkTrack(track)) {
          selected[track.mcParticleId()] = true;
        }
        if (track.hasITS()) {
          hasITS[track.mcParticleId()] = true;
        }
        if (track.hasTPC()) {
          hasTPC[track.mcParticleId()] = true;
        }
      }
    }

    for (const auto pdgCode : pdgCodes.value) {
      auto mass = RecoDecay::getMassPDG(pdgCode);

      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() != pdgCode) {
          continue;
        }
        /// check if we end-up with the correct final state using MC info
        int8_t sign = 0;
        if(std::abs(mcParticle.pdgCode()) == pdg::kD0 && (isThereD0 || isThereD0Bar) && !RecoDecay::isMatchedMCGen(mcParticles, mcParticle, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign) ) {
          /// check if we have D0(bar) → π± K∓
          continue;
        }
        hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, pdgCode, 1.0, true);

        if (std::abs(mcParticle.y()) < 0.5 ) {
          hCandidates->Fill(kHFStepMcInRapidity, mcParticle.pt(), mass, pdgCode, 1.0, true);
        }

        if (mcParticle.daughtersIds().size() != 2) {
          LOGP(fatal, "Invalid numbers of daughters for D0(bar) {}: {}", mcParticle.globalIndex(), mcParticle.daughtersIds().size());
        }

        auto prong0Id = mcParticle.daughtersIds()[0];
        auto prong1Id = mcParticle.daughtersIds()[1];

        if (prong0Id < 0 || prong1Id < 0) {
          continue;
        }

        bool inAcceptance = true;
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        for (const auto& daughter : daughters) {
          if (daughter.pt() < mcAcceptancePt || std::abs(daughter.eta()) > mcAcceptanceEta) {
            inAcceptance = false;
          }
        }

        if (inAcceptance) {
          hCandidates->Fill(kHFStepAcceptance, mcParticle.pt(), mass, pdgCode, 1.0, true);
        }

        if (tracked[prong0Id] && tracked[prong1Id]) {
          hCandidates->Fill(kHFStepTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true);
          if (inAcceptance) {
            hCandidates->Fill(kHFStepAcceptanceTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true);
          } else {
            LOGP(debug, "Candidate {} not in acceptance but tracked.", mcParticle.globalIndex());
            for (const auto& daughter : daughters) {
              LOGP(debug, "   MC: pt={} eta={}", daughter.pt(), daughter.eta());
            }
          }
          for (const auto& daughter : daughters) {
            hTrackablePtEta->Fill(kTrackableAll, daughter.pt(), daughter.eta());
            if (hasITS[daughter.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableITS, daughter.pt(), daughter.eta());
            }
            if (hasTPC[daughter.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableTPC, daughter.pt(), daughter.eta());
            }
            if (hasITS[daughter.globalIndex()] && hasTPC[daughter.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableITSTPC, daughter.pt(), daughter.eta());
            }
          }
        }

        // with track cuts (see checkTrack)
        if (selected[prong0Id] && selected[prong1Id]) {
          hCandidates->Fill(kHFStepTrackableCuts, mcParticle.pt(), mass, pdgCode, 1.0, true);
          if (!inAcceptance) {
            LOGP(info, "Candidate {} not in acceptance but tracked and selected.", mcParticle.globalIndex());
            for (const auto& daughter : daughters) {
              LOGP(info, "   MC: pt={} eta={}", daughter.pt(), daughter.eta());
            }
          }
        }
      }
    }
  }

  void processMcD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls) {
    candidate2ProngMcLoop(candidates, tracks, mcParticles, colls);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcD0, "Process MC for D0 signal", true);

  /// 3-prong analyses

  template <typename C>
  void candidate3ProngMcLoop(C& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls)
  {
    candidate3ProngLoop<true>(candidates, tracks, mcParticles);

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    auto hTrackablePtEta = registry.get<StepTHn>(HIST("hTrackablePtEta"));

    // lists for optimization
    std::vector<bool> tracked(mcParticles.size(), false);
    std::vector<bool> hasITS(mcParticles.size(), false);
    std::vector<bool> hasTPC(mcParticles.size(), false);
    std::vector<bool> selected(mcParticles.size(), false);
    for (const auto& track : tracks) {
      if (track.mcParticleId() >= 0) {
        tracked[track.mcParticleId()] = true;
        if (checkTrack(track)) {
          selected[track.mcParticleId()] = true;
        }
        if (track.hasITS()) {
          hasITS[track.mcParticleId()] = true;
        }
        if (track.hasTPC()) {
          hasTPC[track.mcParticleId()] = true;
        }
      }
    }

    for (const auto pdgCode : pdgCodes.value) { /// loop over PDG codes
      auto mass = RecoDecay::getMassPDG(pdgCode);

      for (const auto& mcParticle : mcParticles) { /// loop over MC particles

        //////////////////////////
        ///   Step kHFStepMC   ///
        //////////////////////////
        if (std::abs(mcParticle.pdgCode()) != pdgCode) { /// abs. value because only "kLambdaCPlus" is defined, not "kAntiLambdaCPlus"
          continue;
        }
        /// check if we end-up with the correct final state using MC info
        int8_t sign = 0;
        if (std::abs(mcParticle.pdgCode()) == pdg::kLambdaCPlus && isThereLc && !RecoDecay::isMatchedMCGen(mcParticles, mcParticle, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2)) {
          /// check if we have Λc± → p± K∓ π± (either direct or resonant)
          continue;
        }
        hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, pdgCode, 1.0, true);

        ////////////////////////////////////
        ///   Step kHFStepMcInRapidity   ///
        ////////////////////////////////////
        if (std::abs(mcParticle.y()) < 0.5 ) {
          hCandidates->Fill(kHFStepMcInRapidity, mcParticle.pt(), mass, pdgCode, 1.0, true);
        }

        auto nDaughters = mcParticle.daughtersIds().size();
        if (nDaughters != 2 && nDaughters != 3) {
          /// # daughters==3: direct decay
          /// # daughters==2: resonant channel
          LOGP(fatal, "Invalid numbers of daughters for 3-prong candidate {}: {}", mcParticle.globalIndex(), mcParticle.daughtersIds().size());
        }

        int32_t prong0Id = -1, prong1Id = -1, prong2Id = -1;
        prong0Id = mcParticle.daughtersIds()[0];
        prong1Id = mcParticle.daughtersIds()[1];
        if (prong0Id < 0 || prong1Id < 0) {
          /// this excludes both resonant and direct channels
          continue;
        }
        if(nDaughters == 3) {
          prong2Id = mcParticle.daughtersIds()[2];
          if (prong2Id < 0) {
            /// this excludes the direct channel
            continue;
          }
        }

        //////////////////////////////////
        ///   Step kHFStepAcceptance   ///
        //////////////////////////////////
        sign = 0;
        bool inAcceptance = true;
        auto daughters = mcParticle.daughters_as<aod::McParticles>();
        nDaughters = daughters.size();

        //
        // Λc± → p± K∓ π± case
        bool isDirect = false;
        bool isResonanceFound = false;
        int resoId = -1;
        //std::array<int, 2> resoProngsId = {-1, -1};
        int otherDaughterId = -1;
        if (isThereLc) {
          if(nDaughters == 3) {
            // direct channel: Λc± → p± K∓ π±
            isDirect = true;
            for (const auto& daughter : daughters) {
              if (daughter.pt() < mcAcceptancePt || std::abs(daughter.eta()) > mcAcceptanceEta) {
                inAcceptance = false;
              }
            }
          } else {
            // this means resonant channels (i.e. 2 prongs)
            
            //bool isResonanceFound = false;
            //int resoId = -1;
            //int otherDaughterId = -1;
            std::vector<int> arrDaughIndex;
            std::vector<int> arrPDGDaugh;
            RecoDecay::getDaughters(mcParticle, &arrDaughIndex, array{0}, 1); /// get indices of Lc daughters (resonant channels)
            if (arrDaughIndex.size() != 2) {
              continue;
            }
            
            // resonant channel Λc± → p± K* (kProton, 313)
            if (!isResonanceFound) {
              for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) { /// loop of Lc daughters (resonant channel)
                auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
                arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
                if(arrPDGDaugh[jProng] == 313) {
                  /// this is the K*
                  resoId = arrDaughIndex[jProng];
                } else if (arrPDGDaugh[jProng] == kProton) {
                  /// this is the proton
                  otherDaughterId = arrDaughIndex[jProng];
                } else {
                  /// this means that the current case is not a Λc± → p± K*
                  /// let's move directly to the next case
                  break;
                }
                if(resoId > 0 && otherDaughterId > 0) {
                  /// we found the Λc± → p± K* resonant channel!
                  isResonanceFound = true;
                }
              } /// end loop of Lc daughters (resonant channel)
            }

            // resonant channel Λc± → Δ(1232)±± K∓ (2224, kKPlus)
            if (!isResonanceFound) {
              for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) { /// loop of Lc daughters (resonant channel)
                auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
                arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
                if(arrPDGDaugh[jProng] == 2224) {
                  /// this is the Δ(1232)±±
                  resoId = arrDaughIndex[jProng];
                } else if (arrPDGDaugh[jProng] == kKPlus) {
                  /// this is the kaon
                  otherDaughterId = arrDaughIndex[jProng];
                } else {
                  /// this means that the current case is not a Λc± → Δ(1232)±± K∓
                  /// let's move directly to the next case
                  break;
                }
                if(resoId > 0 && otherDaughterId > 0) {
                  /// we found the Λc± → Δ(1232)±± K∓ resonant channel!
                  isResonanceFound = true;
                }
              } /// end loop of Lc daughters (resonant channel)
            }

            // resonant channel Λc± → Λ(1520) π± (3124, kPiPlus)
            if (!isResonanceFound) {
              for (auto jProng = 0u; jProng < arrDaughIndex.size(); ++jProng) { /// loop of Lc daughters (resonant channel)
                auto daughJ = mcParticles.rawIteratorAt(arrDaughIndex[jProng]);
                arrPDGDaugh[jProng] = std::abs(daughJ.pdgCode());
                if(arrPDGDaugh[jProng] == 3124) {
                  /// this is the Λ(1520)
                  resoId = arrDaughIndex[jProng];
                } else if (arrPDGDaugh[jProng] == kPiPlus) {
                  /// this is the pion
                  otherDaughterId = arrDaughIndex[jProng];
                } else {
                  /// this means that the current case is not a Λ(1520) π±
                  /// let's move directly to the next case
                  break;
                }
                if(resoId > 0 && otherDaughterId > 0) {
                  /// we found the Λ(1520) π± resonant channel!
                  isResonanceFound = true;
                }
              } /// end loop of Lc daughters (resonant channel)
            }
            if(!isResonanceFound) {
              continue;
            }

            /// check if the prong of resonance daughter are in acceptance
            auto resoProngs = mcParticles.rawIteratorAt(resoId).daughters_as<aod::McParticles>();
            if(resoProngs.size() != 2) {
              continue;
            }
            for (const auto& resoProng : resoProngs) {
              if (resoProng.pt() < mcAcceptancePt || std::abs(resoProng.eta()) > mcAcceptanceEta) {
                inAcceptance = false;
              }
            }
            prong0Id = mcParticles.rawIteratorAt(resoId).daughtersIds()[0];
            prong1Id = mcParticles.rawIteratorAt(resoId).daughtersIds()[1];
            /// check if the other daughter is in acceptance
            auto otherDaug = mcParticles.rawIteratorAt(otherDaughterId);
            if (otherDaug.pt() < mcAcceptancePt || std::abs(otherDaug.eta()) > mcAcceptanceEta) {
              inAcceptance = false;
            }
            prong2Id = otherDaughterId;

          } /// end study resonant channels Λc± → p± K∓ π± decays
        }
        // end Λc± → p± K∓ π± case
        //
        else {
          /// this should be changed if we add more 3-prong particles, eg D+, Ds+, ...
          LOGP(fatal, ">>> test inAcceptance step for PDG {} not implemented", pdgCode);
        }

        if (inAcceptance) {
          hCandidates->Fill(kHFStepAcceptance, mcParticle.pt(), mass, pdgCode, 1.0, true);
        }

        /////////////////////////////////
        ///   Step kHFStepTrackable   ///
        /////////////////////////////////
        if (tracked[prong0Id] && tracked[prong1Id] && tracked[prong2Id]) {
          hCandidates->Fill(kHFStepTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true);

          ///////////////////////////////////////////
          ///   Step kHFStepAcceptanceTrackable   ///
          ///////////////////////////////////////////
          if (inAcceptance) {
            hCandidates->Fill(kHFStepAcceptanceTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true);
          } else {
            LOGP(debug, "Candidate {} not in acceptance but tracked.", mcParticle.globalIndex());
            //for (const auto& daughter : daughters) {
            //  LOGP(debug, "   MC: pt={} eta={}", daughter.pt(), daughter.eta());
            //}
          }

          /// +++ direct channel +++
          if (isDirect) {
            
            for (const auto& daughter : daughters) {
              //////////////////////////////
              ///   Step kTrackableAll   ///
              //////////////////////////////
              hTrackablePtEta->Fill(kTrackableAll, daughter.pt(), daughter.eta());
              //////////////////////////////
              ///   Step kTrackableITS   ///
              //////////////////////////////
              if (hasITS[daughter.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableITS, daughter.pt(), daughter.eta());
              }
              //////////////////////////////
              ///   Step kTrackableTPC   ///
              //////////////////////////////
              if (hasTPC[daughter.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableTPC, daughter.pt(), daughter.eta());
              }
              /////////////////////////////////
              ///   Step kTrackableITSTPC   ///
              /////////////////////////////////
              if (hasITS[daughter.globalIndex()] && hasTPC[daughter.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableITSTPC, daughter.pt(), daughter.eta());
              }
            }

          }
          /// +++ resonant channel +++
          else if (isResonanceFound) {
            
            /// +++ resonance daughters +++
            auto resoProngs = mcParticles.rawIteratorAt(resoId).daughters_as<aod::McParticles>();
            for (const auto& resoProng : resoProngs) {
              //////////////////////////////
              ///   Step kTrackableAll   ///
              //////////////////////////////
              hTrackablePtEta->Fill(kTrackableAll, resoProng.pt(), resoProng.eta());
              //////////////////////////////
              ///   Step kTrackableITS   ///
              //////////////////////////////
              if (hasITS[resoProng.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableITS, resoProng.pt(), resoProng.eta());
              }
              //////////////////////////////
              ///   Step kTrackableTPC   ///
              //////////////////////////////
              if (hasTPC[resoProng.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableTPC, resoProng.pt(), resoProng.eta());
              }
              /////////////////////////////////
              ///   Step kTrackableITSTPC   ///
              /////////////////////////////////
              if (hasITS[resoProng.globalIndex()] && hasTPC[resoProng.globalIndex()]) {
                hTrackablePtEta->Fill(kTrackableITSTPC, resoProng.pt(), resoProng.eta());
              }
            }

            /// +++ other daughter +++
            auto otherDaug = mcParticles.rawIteratorAt(otherDaughterId);
            //////////////////////////////
            ///   Step kTrackableAll   ///
            //////////////////////////////
            hTrackablePtEta->Fill(kTrackableAll, otherDaug.pt(), otherDaug.eta());
            //////////////////////////////
            ///   Step kTrackableITS   ///
            //////////////////////////////
            if (hasITS[otherDaug.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableITS, otherDaug.pt(), otherDaug.eta());
            }
            //////////////////////////////
            ///   Step kTrackableTPC   ///
            //////////////////////////////
            if (hasTPC[otherDaug.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableTPC, otherDaug.pt(), otherDaug.eta());
            }
            /////////////////////////////////
            ///   Step kTrackableITSTPC   ///
            /////////////////////////////////
            if (hasITS[otherDaug.globalIndex()] && hasTPC[otherDaug.globalIndex()]) {
              hTrackablePtEta->Fill(kTrackableITSTPC, otherDaug.pt(), otherDaug.eta());
            }
          }
        } /// end "if(tracked[...])"

        /////////////////////////////////////
        ///   Step kHFStepTrackableCuts   ///
        /////////////////////////////////////
        if (selected[prong0Id] && selected[prong1Id] && selected[prong2Id]) {
          hCandidates->Fill(kHFStepTrackableCuts, mcParticle.pt(), mass, pdgCode, 1.0, true);
          if (!inAcceptance) {
            LOGP(info, "Candidate {} not in acceptance but tracked and selected.", mcParticle.globalIndex());
          }
        } /// end "if(selected[...])"

      } /// end loop over MC particles
    } /// end loop over PDG codes

  } /// end candidate3ProngMcLoop

  void processDataLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc>& candidates, TracksWithSelection& tracks) {
    candidate3ProngLoop<false>(candidates, tracks, tracks);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataLc, "Process Lc data (no MC information needed)", false);

  void processMcLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc>& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls) {
    candidate3ProngMcLoop(candidates, tracks, mcParticles, colls);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcLc, "Process MC for Lc signal", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcEfficiency>(cfgc)};
}
