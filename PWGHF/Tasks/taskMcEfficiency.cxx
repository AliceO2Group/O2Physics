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
using namespace o2::analysis::hf_cuts_d0_to_pi_k;

struct HfTaskMcEfficiency {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  ConfigurableAxis axisPt{"axisPt", {10, 0, 10}, "pT axis"};
  ConfigurableAxis axisMass{"axisMass", {120, 1.5848, 2.1848}, "m_inv axis"};
  ConfigurableAxis axisPdg{"axisPdg", {VARIABLE_WIDTH, -4122.5, -421.5, 0, 421.5, 4122.5}, "PDG code axis"};
  ConfigurableAxis axisCPA{"axisCPA", {102, -1.02, 1.02}, "Cosine of pointing angle axis"};

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

  void init(InitContext&)
  {

    std::array<bool, 2> doprocessData{doprocessDataD0, doprocessDataLc};
    std::array<bool, 2> doprocessMc{doprocessMcD0, doprocessMcLc};
    if (std::accumulate(doprocessData.begin(), doprocessData.end(), 0) > 0 && std::accumulate(doprocessMc.begin(), doprocessMc.end(), 0) > 0) {
      LOGP(fatal, "Data and MC process functions cannot run simultaneously!");
    }

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPt, axisMass, axisPdg, axisCPA, {2, -0.5, 1.5, "collision matched"}, {RecoDecay::OriginType::NonPrompt + 1, RecoDecay::OriginType::None - 0.5, RecoDecay::OriginType::NonPrompt + 0.5}}, kHFNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
    hCandidates->GetAxis(2)->SetTitle("PDG code");
    hCandidates->GetAxis(3)->SetTitle("CPA");
    hCandidates->GetAxis(5)->SetTitle("Charm hadron origin");

    registry.add("hDuplicateCount", "Duplicate count;frequency;count", {HistType::kTH1F, {{10, 0.5, 10.5}}});

    auto hTrackablePtEta = registry.add<StepTHn>("hTrackablePtEta", "Prongs kinematics at different steps", {HistType::kStepTHnF, {{20, 0, 10}, {40, -2, 2}}, kNTrackableSteps});
    hTrackablePtEta->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hTrackablePtEta->GetAxis(1)->SetTitle("#eta");
  }

  template <typename T>
  inline bool checkTrack(T track)
  {
    // TODO configurable?
    return track.isGlobalTrackWoDCA();
  }

  template <bool mc, typename T1, typename T2, typename T3>
  void candidate3ProngLoop(T1& candidates, T2& tracks, T3& mcParticles, std::vector<int> pdgCodes)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<std::size_t, int> duplicates;

    for (const auto& candidate : candidates) { /// loop over candidates

      for (const auto pdgCode : pdgCodes) { /// loop on pdg codes
        auto decayType = -1;
        std::array<int, 3> pdgDaughters;

        if (pdgCode == pdg::kLambdaCPlus) {
          decayType = 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;
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
          isHypoMass1SelStep = candidate.isSelLcToPKPi(); /// from candidate selector!
          isHypoMass2SelStep = candidate.isSelLcToPiKP(); /// from candidate selector!
        }
        bool collisionMatched = false;
        int origin = RecoDecay::OriginType::None;
        if constexpr (mc) { /// info MC used
          int8_t sign = 0;
          int indexRec = -999;
          if (pdgCode == pdg::kLambdaCPlus) {
            indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg, trackThird}, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2);
          }

          if (indexRec < 0) {
            continue;
          }

          origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticles.rawIteratorAt(indexRec));

          if (pdgCode == pdg::kLambdaCPlus) {
            auto daughter = trackPos.mcParticle();
            if (std::abs(daughter.pdgCode()) == kProton) {
              isHypoMass1TrackStep = true;
              isHypoMass1SelStep = true;
              isHypoMass2TrackStep = false;
              isHypoMass2SelStep = false;
            } else if (std::abs(daughter.pdgCode()) == kPiPlus) {
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
        // bool selected = false;

        /// all candidates
        if (isHypoMass1TrackStep) {
          if (pdgCode == pdg::kLambdaCPlus) {
            massHypo1 = invMassLcToPKPi(candidate);
          }
          hCandidates->Fill(kHFStepTracked, pt, massHypo1, pdgCode, cpa, collisionMatched, origin);
        }
        if (isHypoMass2TrackStep) {
          if (pdgCode == pdg::kLambdaCPlus) {
            massHypo2 = invMassLcToPiKP(candidate);
          }
          hCandidates->Fill(kHFStepTracked, pt, massHypo2, pdgCode, cpa, collisionMatched, origin);
        }

        // check if prongs have passed track cuts
        if (checkTrack(trackPos) && checkTrack(trackNeg) && checkTrack(trackThird)) {
          if (isHypoMass1TrackStep) {
            hCandidates->Fill(kHFStepTrackedCuts, pt, massHypo1, pdgCode, cpa, collisionMatched, origin);
          }
          if (isHypoMass2TrackStep) {
            hCandidates->Fill(kHFStepTrackedCuts, pt, massHypo2, pdgCode, cpa, collisionMatched, origin);
          }
        }

        if (!isHypoMass1SelStep && !isHypoMass2SelStep) {
          continue;
        }

        // selected candidates
        if (isHypoMass1SelStep) {
          hCandidates->Fill(kHFStepTrackedSelected, pt, massHypo1, pdgCode, cpa, collisionMatched, origin);
        }
        if (isHypoMass2SelStep) {
          hCandidates->Fill(kHFStepTrackedSelected, pt, massHypo2, pdgCode, cpa, collisionMatched, origin);
        }

        // duplicates
        std::array<int, 3> prongIds = {candidate.prong0Id(), candidate.prong1Id(), candidate.prong2Id()};
        std::sort(prongIds.begin(), prongIds.end());
        std::string concat = std::to_string(prongIds[0]) + std::to_string(prongIds[1]) + std::to_string(prongIds[2]);
        std::size_t hash = std::hash<std::string>{}(concat); /// unique value for the 'concat' string
        if (duplicates.find(hash) != duplicates.end()) {
          if (isHypoMass1TrackStep) {
            hCandidates->Fill(kHFStepTrackedDuplicates, pt, massHypo1, pdgCode, cpa, collisionMatched, origin);
          }
          if (isHypoMass2TrackStep) {
            hCandidates->Fill(kHFStepTrackedDuplicates, pt, massHypo2, pdgCode, cpa, collisionMatched, origin);
          }
        }
        duplicates[hash]++;

      } /// end loop on pdg codes
    }   /// end loop over candidates

    auto hDuplicateCount = registry.get<TH1>(HIST("hDuplicateCount"));
    for (const auto& i : duplicates) {
      hDuplicateCount->Fill(i.second);
    }
  }

  template <bool mc, typename T1, typename T2, typename T3>
  void candidate2ProngLoop(T1& candidates, T2& tracks, T3& mcParticles, std::vector<int> pdgCodes)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<int64_t, int> duplicates;

    for (const auto pdgCode : pdgCodes) {
      auto decayType = -1;
      std::array<int, 2> pdgDaughters;

      if (pdgCode == pdg::kD0) {
        decayType = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
        pdgDaughters[0] = +kPiPlus;
        pdgDaughters[1] = -kKPlus;
      } else if (pdgCode == pdg::kD0Bar) {
        decayType = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
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
        int origin = RecoDecay::OriginType::None;
        if constexpr (mc) {
          auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg}, pdgCode, pdgDaughters, false);
          if (indexRec < 0) {
            continue;
          }
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticles.rawIteratorAt(indexRec));

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
        hCandidates->Fill(kHFStepTracked, pt, mass, pdgCode, cpa, collisionMatched, origin);
        ++nTracked;

        // check if prongs have passed track cuts
        if (checkTrack(trackPos) && checkTrack(trackNeg)) {
          hCandidates->Fill(kHFStepTrackedCuts, pt, mass, pdgCode, cpa, collisionMatched, origin);
        }

        if (!selected) {
          continue;
        }

        // selected candidates
        hCandidates->Fill(kHFStepTrackedSelected, pt, mass, pdgCode, cpa, collisionMatched, origin);
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
          hCandidates->Fill(kHFStepTrackedDuplicates, pt, mass, pdgCode, cpa, collisionMatched, origin);
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

  template <typename C>
  void candidate2ProngMcLoop(C& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls, std::vector<int> pdgCodes)
  {
    candidate2ProngLoop<true>(candidates, tracks, mcParticles, pdgCodes);

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

    for (const auto pdgCode : pdgCodes) {
      auto mass = RecoDecay::getMassPDG(pdgCode);

      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() != pdgCode) {
          continue;
        }
        /// check if we end-up with the correct final state using MC info
        int8_t sign = 0;
        if (std::abs(mcParticle.pdgCode()) == pdg::kD0 && !RecoDecay::isMatchedMCGen(mcParticles, mcParticle, pdg::Code::kD0, array{+kPiPlus, -kKPlus}, true, &sign)) {
          /// check if we have D0(bar) → π± K∓
          continue;
        }

        int origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle);

        hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);

        if (std::abs(mcParticle.y()) < 0.5) {
          hCandidates->Fill(kHFStepMcInRapidity, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);
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
          hCandidates->Fill(kHFStepAcceptance, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);
        }

        if (tracked[prong0Id] && tracked[prong1Id]) {
          hCandidates->Fill(kHFStepTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);
          if (inAcceptance) {
            hCandidates->Fill(kHFStepAcceptanceTrackable, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);
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
          hCandidates->Fill(kHFStepTrackableCuts, mcParticle.pt(), mass, pdgCode, 1.0, true, origin);
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

  /// 3-prong analyses

  template <typename C>
  void candidate3ProngMcLoop(C& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls, std::vector<int> pdgCodes)
  {
    candidate3ProngLoop<true>(candidates, tracks, mcParticles, pdgCodes);

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

    for (const auto pdgCode : pdgCodes) { /// loop over PDG codes
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
        std::unique_ptr<std::vector<int>> listIndexDaughters(new std::vector<int>{});
        if (std::abs(mcParticle.pdgCode()) == pdg::kLambdaCPlus && !RecoDecay::isMatchedMCGen(mcParticles, mcParticle, pdg::Code::kLambdaCPlus, array{+kProton, -kKPlus, +kPiPlus}, true, &sign, 2, listIndexDaughters.get())) {
          /// check if we have Λc± → p± K∓ π± (either direct or resonant)
          continue;
        }

        int origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle);

        hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);

        ////////////////////////////////////
        ///   Step kHFStepMcInRapidity   ///
        ////////////////////////////////////
        if (std::abs(mcParticle.y()) < 0.5) {
          hCandidates->Fill(kHFStepMcInRapidity, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);
        }

        auto nDaughters = listIndexDaughters.get()->size();
        if (nDaughters != 2 && nDaughters != 3) {
          /// # daughters==3: direct decay
          /// # daughters==2: resonant channel
          LOGP(fatal, "Invalid numbers of daughters for 3-prong candidate {}: {}", mcParticle.globalIndex(), listIndexDaughters.get()->size());
        }

        bool hasBadDaughter = false;
        for (const auto& prongId : *listIndexDaughters.get()) {
          if (prongId < 0) {
            hasBadDaughter = true;
            break;
          }
        }
        if (hasBadDaughter) {
          continue;
        }

        //////////////////////////////////
        ///   Step kHFStepAcceptance   ///
        //////////////////////////////////
        bool inAcceptance = true;
        for (const auto& prongId : *listIndexDaughters.get()) {
          auto daughter = mcParticles.rawIteratorAt(prongId);
          if (daughter.pt() < mcAcceptancePt || std::abs(daughter.eta()) > mcAcceptanceEta) {
            inAcceptance = false;
          }
        }

        if (inAcceptance) {
          hCandidates->Fill(kHFStepAcceptance, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);
        }

        /////////////////////////////////
        ///   Step kHFStepTrackable   ///
        /////////////////////////////////
        if (tracked[listIndexDaughters.get()->at(0)] && tracked[listIndexDaughters.get()->at(1)] && tracked[listIndexDaughters.get()->at(2)]) {
          hCandidates->Fill(kHFStepTrackable, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);

          ///////////////////////////////////////////
          ///   Step kHFStepAcceptanceTrackable   ///
          ///////////////////////////////////////////
          if (inAcceptance) {
            hCandidates->Fill(kHFStepAcceptanceTrackable, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);
          } else {
            LOGP(debug, "Candidate {} not in acceptance but tracked.", mcParticle.globalIndex());
          }

          for (const auto& prongId : *listIndexDaughters.get()) {
            auto daughter = mcParticles.rawIteratorAt(prongId);
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
        } /// end "if(tracked[...])"

        /////////////////////////////////////
        ///   Step kHFStepTrackableCuts   ///
        /////////////////////////////////////
        if (selected[listIndexDaughters.get()->at(0)] && selected[listIndexDaughters.get()->at(1)] && selected[listIndexDaughters.get()->at(2)]) {
          hCandidates->Fill(kHFStepTrackableCuts, mcParticle.pt(), mass, pdgCode * sign, 1.0, true, origin);
          if (!inAcceptance) {
            LOGP(info, "Candidate {} not in acceptance but tracked and selected.", mcParticle.globalIndex());
          }
        } /// end "if(selected[...])"

      } /// end loop over MC particles
    }   /// end loop over PDG codes

  } /// end candidate3ProngMcLoop

  // process functions for data
  void processDataD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelection& tracks)
  {
    std::vector<int> pdgCodes{pdg::kD0Bar, pdg::kD0};
    candidate2ProngLoop<false>(candidates, tracks, tracks, pdgCodes); // NOTE third argument has to be provided but is not used as template argument is <false>
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataD0, "Process D0 data (no MC information needed)", false);

  void processDataLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc>& candidates, TracksWithSelection& tracks)
  {
    std::vector<int> pdgCodes{pdg::kLambdaCPlus};
    candidate3ProngLoop<false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataLc, "Process Lc data (no MC information needed)", false);

  // process functions for MC
  void processMcD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls)
  {
    std::vector<int> pdgCodes{pdg::kD0Bar, pdg::kD0};
    candidate2ProngMcLoop(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcD0, "Process MC for D0 signal", true);

  void processMcLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc>& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels& colls)
  {
    std::vector<int> pdgCodes{pdg::kLambdaCPlus};
    candidate3ProngMcLoop(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcLc, "Process MC for Lc signal", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcEfficiency>(cfgc)};
}
