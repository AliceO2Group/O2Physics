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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "PWGHF/DataModel/TrackIndexSkimmingTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/StepTHn.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskMcEfficiency {
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection Flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection Flag for D0bar"};

  ConfigurableAxis axisPt{"axisPt", {10, 0, 10}, "pT axis"};
  ConfigurableAxis axisMass{"axisMass", {120, 1.5848, 2.1848}, "m_inv axis"};
  ConfigurableAxis axisPdg{"axisPdg", {VARIABLE_WIDTH, -4232.5, -4122.5, -431.5, -421.5, -411.5, 0, 411.5, 421.5, 431.5, 4122.5, 4232.5}, "PDG code axis"};
  ConfigurableAxis axisCPA{"axisCPA", {102, -1.02, 1.02}, "Cosine of pointing angle axis"};

  Configurable<float> mcAcceptancePt{"mcAcceptancePt", 0.1, "MC Acceptance: lower pt limit"};
  Configurable<float> mcAcceptanceEta{"mcAcceptanceEta", 0.8, "MC Acceptance: upper eta limit"};

  Service<o2::framework::O2DatabasePDG> pdg;

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

    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPt, axisMass, axisPdg, axisCPA, {2, -0.5, 1.5, "collision matched"}, {RecoDecay::OriginType::NonPrompt + 1, +RecoDecay::OriginType::None - 0.5, +RecoDecay::OriginType::NonPrompt + 0.5}}, kHFNSteps});
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
  bool checkTrack(T track)
  {
    // TODO configurable?
    return track.isGlobalTrackWoDCA();
  }

  template <bool Mc, bool HasDplus, bool HasDs, bool HasLc, bool HasXicPlus, typename T1, typename T2, typename T3>
  void candidate3ProngLoop(T1 const& candidates, T2 const& tracks, T3 const& mcParticles, std::vector<int> const& pdgCodes)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<std::size_t, int> duplicates;

    for (const auto& candidate : candidates) { /// loop over candidates

      for (const auto pdgCode : pdgCodes) { /// loop on pdg codes
        auto decayType = -1;
        std::array<int, 3> pdgDaughters{};

        if (pdgCode == Pdg::kDPlus) {
          decayType = 1 << aod::hf_cand_3prong::DecayType::DplusToPiKPi;
          pdgDaughters[0] = +kPiPlus;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kDS) {
          decayType = 1 << aod::hf_cand_3prong::DecayType::DsToKKPi;
          pdgDaughters[0] = +kKPlus;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kLambdaCPlus) {
          decayType = 1 << aod::hf_cand_3prong::DecayType::LcToPKPi;
          pdgDaughters[0] = +kProton;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kXiCPlus) {
          decayType = 1 << aod::hf_cand_3prong::DecayType::XicToPKPi;
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
        /// selections from candidate selectors
        if constexpr (HasDplus) {
          if (pdgCode == Pdg::kDPlus) {
            isHypoMass1SelStep = candidate.isSelDplusToPiKPi(); // only one mass hypo for D+
          }
        }
        if constexpr (HasDs) {
          if (pdgCode == Pdg::kDS) {
            isHypoMass1SelStep = candidate.isSelDsToKKPi();
            isHypoMass2SelStep = candidate.isSelDsToPiKK();
          }
        }
        if constexpr (HasLc) {
          if (pdgCode == Pdg::kLambdaCPlus) {
            isHypoMass1SelStep = candidate.isSelLcToPKPi();
            isHypoMass2SelStep = candidate.isSelLcToPiKP();
          }
        }
        if constexpr (HasXicPlus) {
          if (pdgCode == Pdg::kXiCPlus) {
            isHypoMass1SelStep = candidate.isSelXicToPKPi();
            isHypoMass2SelStep = candidate.isSelXicToPiKP();
          }
        }

        bool collisionMatched = false;
        int origin = RecoDecay::OriginType::None;
        if constexpr (Mc) { /// info MC used
          int8_t sign = 0;
          int const indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg, trackThird}, pdgCode, pdgDaughters, true, &sign, 2);

          if (indexRec < 0) {
            continue;
          }

          origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticles.rawIteratorAt(indexRec));

          if (pdgCode == Pdg::kLambdaCPlus) {
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

          if (pdgCode == Pdg::kXiCPlus) {
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
        float const cpa = candidate.cpa();
        float const pt = candidate.pt();
        // bool selected = false;

        /// all candidates
        if (isHypoMass1TrackStep) {
          if (pdgCode == Pdg::kLambdaCPlus) {
            massHypo1 = HfHelper::invMassLcToPKPi(candidate);
          } else if (pdgCode == Pdg::kXiCPlus) {
            massHypo1 = HfHelper::invMassXicToPKPi(candidate);
          } else if (pdgCode == Pdg::kDPlus) {
            massHypo1 = HfHelper::invMassDplusToPiKPi(candidate);
          } else if (pdgCode == Pdg::kDS) {
            massHypo1 = HfHelper::invMassDsToKKPi(candidate);
          }
          hCandidates->Fill(kHFStepTracked, pt, massHypo1, pdgCode, cpa, collisionMatched, origin);
        }
        if (isHypoMass2TrackStep) {
          if (pdgCode == Pdg::kLambdaCPlus) {
            massHypo2 = HfHelper::invMassLcToPiKP(candidate);
          } else if (pdgCode == Pdg::kXiCPlus) {
            massHypo2 = HfHelper::invMassXicToPiKP(candidate);
          } else if (pdgCode == Pdg::kDS) {
            massHypo2 = HfHelper::invMassDsToPiKK(candidate);
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
        std::string const concat = std::to_string(prongIds[0]) + std::to_string(prongIds[1]) + std::to_string(prongIds[2]);
        std::size_t const hash = std::hash<std::string>{}(concat); /// unique value for the 'concat' string
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
    } /// end loop over candidates

    auto hDuplicateCount = registry.get<TH1>(HIST("hDuplicateCount"));
    for (const auto& i : duplicates) {
      hDuplicateCount->Fill(i.second);
    }
  }

  template <bool Mc, typename T1, typename T2, typename T3>
  void candidate2ProngLoop(T1 const& candidates, T2 const& tracks, T3 const& mcParticles, std::vector<int> const& pdgCodes)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<int64_t, int> duplicates;

    for (const auto pdgCode : pdgCodes) {
      auto decayType = -1;
      std::array<int, 2> pdgDaughters{};

      if (pdgCode == Pdg::kD0) {
        decayType = 1 << aod::hf_cand_2prong::DecayType::D0ToPiK;
        pdgDaughters[0] = +kPiPlus;
        pdgDaughters[1] = -kKPlus;
      } else if (pdgCode == Pdg::kD0Bar) {
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
        if constexpr (Mc) {
          auto indexRec = RecoDecay::getMatchedMCRec(mcParticles, std::array{trackPos, trackNeg}, pdgCode, pdgDaughters, false);
          if (indexRec < 0) {
            continue;
          }
          origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticles.rawIteratorAt(indexRec));

          collisionMatched = candidate.template collision_as<aod::McCollisionLabels>().mcCollisionId() == mcParticles.iteratorAt(indexRec).mcCollisionId();
        }

        float mass = -1;
        float const cpa = candidate.cpa();
        float const pt = candidate.pt();
        bool selected = false;
        if (pdgCode == Pdg::kD0) {
          mass = HfHelper::invMassD0ToPiK(candidate);
          selected = candidate.isSelD0() >= selectionFlagD0;
        } else if (pdgCode == Pdg::kD0Bar) {
          mass = HfHelper::invMassD0barToKPi(candidate);
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
          hash = (static_cast<int64_t>(candidate.prong0Id()) << 32) | candidate.prong1Id();
        } else {
          hash = (static_cast<int64_t>(candidate.prong1Id()) << 32) | candidate.prong0Id();
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
  void candidate2ProngMcLoop(C const& candidates, TracksWithSelectionMC const& tracks, aod::McParticles const& mcParticles, aod::McCollisionLabels const&, std::vector<int> pdgCodes)
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
      auto mass = pdg->Mass(pdgCode);

      for (const auto& mcParticle : mcParticles) {
        if (mcParticle.pdgCode() != pdgCode) {
          continue;
        }
        /// check if we end-up with the correct final state using MC info
        int8_t sign = 0;
        if (std::abs(mcParticle.pdgCode()) == Pdg::kD0 && !RecoDecay::isMatchedMCGen(mcParticles, mcParticle, Pdg::kD0, std::array{+kPiPlus, -kKPlus}, true, &sign)) {
          /// check if we have D0(bar) → π± K∓
          continue;
        }

        int const origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle);

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

  template <bool HasDplus, bool HasDs, bool HasLc, bool HasXicPlus, typename C>
  void candidate3ProngMcLoop(C const& candidates, TracksWithSelectionMC const& tracks, aod::McParticles const& mcParticles, aod::McCollisionLabels const&, std::vector<int> pdgCodes)
  {
    candidate3ProngLoop<true, HasDplus, HasDs, HasLc, HasXicPlus>(candidates, tracks, mcParticles, pdgCodes);

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
      auto mass = pdg->Mass(pdgCode);

      for (const auto& mcParticle : mcParticles) { /// loop over MC particles

        //////////////////////////
        ///   Step kHFStepMC   ///
        //////////////////////////
        auto absPdg = std::abs(mcParticle.pdgCode());
        if (absPdg != pdgCode) { /// abs. value because only "kLambdaCPlus" is defined, not "kAntiLambdaCPlus"
          continue;
        }

        std::array<int, 3> pdgDaughters{};
        if (pdgCode == Pdg::kDPlus) {
          pdgDaughters[0] = +kPiPlus;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kDS) {
          pdgDaughters[0] = +kKPlus;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kLambdaCPlus) {
          pdgDaughters[0] = +kProton;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else if (pdgCode == Pdg::kXiCPlus) {
          pdgDaughters[0] = +kProton;
          pdgDaughters[1] = -kKPlus;
          pdgDaughters[2] = +kPiPlus;
        } else {
          LOGP(fatal, "Not implemented for PDG {}", pdgCode);
        }

        /// check if we end-up with the correct final state using MC info
        int8_t sign = 0;
        std::unique_ptr<std::vector<int>> const listIndexDaughters(new std::vector<int>{});
        if (!RecoDecay::isMatchedMCGen(mcParticles, mcParticle, pdgCode, pdgDaughters, true, &sign, 2, listIndexDaughters.get())) {
          /// check if we have Λc± → p± K∓ π± (either direct or resonant), or D± → π± K∓ π± (either direct or resonant) or Ds± → K± K∓ π± (either direct or resonant)
          continue;
        }

        int const origin = RecoDecay::getCharmHadronOrigin(mcParticles, mcParticle);

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
    } /// end loop over PDG codes

  } /// end candidate3ProngMcLoop

  // process functions for data
  void processDataD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                     TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kD0Bar, Pdg::kD0};
    candidate2ProngLoop<false>(candidates, tracks, tracks, pdgCodes); // NOTE third argument has to be provided but is not used as template argument is <false>
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataD0, "Process D0 data (no MC information needed)", false);

  void processDataDplus(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const& candidates,
                        TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus};
    candidate3ProngLoop<false, true, false, false, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDplus, "Process D+ data (no MC information needed)", false);

  void processDataDs(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const& candidates,
                     TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDS};
    candidate3ProngLoop<false, false, true, false, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDs, "Process Ds+ data (no MC information needed)", false);

  void processDataLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                     TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kLambdaCPlus};
    candidate3ProngLoop<false, false, false, true, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataLc, "Process Lc data (no MC information needed)", false);

  void processDataXic(soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi> const& candidates,
                      TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kXiCPlus};
    candidate3ProngLoop<false, false, false, false, true>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataXic, "Process Xic data (no MC information needed)", false);

  void processDataDplusDs(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelDsToKKPi> const& candidates,
                          TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kDS};
    candidate3ProngLoop<false, true, true, false, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDplusDs, "Process D+ and Ds+ data (no MC information needed)", false);

  void processDataDplusDsLc(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelDsToKKPi, aod::HfSelLc> const& candidates,
                            TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kDS, Pdg::kLambdaCPlus};
    candidate3ProngLoop<false, true, true, true, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDplusDsLc, "Process D+, Ds+, and Lc data (no MC information needed)", false);

  void processDataDplusLc(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelLc> const& candidates,
                          TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kLambdaCPlus};
    candidate3ProngLoop<false, true, false, true, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDplusLc, "Process D+ and Lc data (no MC information needed)", false);

  void processDataDsLc(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfSelLc> const& candidates,
                       TracksWithSelection const& tracks)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kDS, Pdg::kLambdaCPlus};
    candidate3ProngLoop<false, false, true, true, false>(candidates, tracks, tracks, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processDataDsLc, "Process Ds+ and Lc data (no MC information needed)", false);

  // process functions for MC
  void processMcD0(soa::Join<aod::HfCand2Prong, aod::HfSelD0> const& candidates,
                   TracksWithSelectionMC const& tracks,
                   aod::McParticles const& mcParticles,
                   aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kD0Bar, Pdg::kD0};
    candidate2ProngMcLoop(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcD0, "Process MC for D0 signal", true);

  void processMcDplus(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi> const& candidates,
                      TracksWithSelectionMC const& tracks,
                      aod::McParticles const& mcParticles,
                      aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus};
    candidate3ProngMcLoop<true, false, false, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDplus, "Process MC for D+ signal", false);

  void processMcDs(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const& candidates,
                   TracksWithSelectionMC const& tracks,
                   aod::McParticles const& mcParticles,
                   aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDS};
    candidate3ProngMcLoop<false, true, false, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDs, "Process MC for Ds+ signal", false);

  void processMcLc(soa::Join<aod::HfCand3Prong, aod::HfSelLc> const& candidates,
                   TracksWithSelectionMC const& tracks,
                   aod::McParticles const& mcParticles,
                   aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kLambdaCPlus};
    candidate3ProngMcLoop<false, false, true, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcLc, "Process MC for Lc signal", false);

  void processMcXic(soa::Join<aod::HfCand3Prong, aod::HfSelXicToPKPi> const& candidates,
                    TracksWithSelectionMC const& tracks,
                    aod::McParticles const& mcParticles,
                    aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kXiCPlus};
    candidate3ProngMcLoop<false, false, false, true>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcXic, "Process MC for Xic signal", false);

  void processMcDplusDs(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelDsToKKPi> const& candidates,
                        TracksWithSelectionMC const& tracks,
                        aod::McParticles const& mcParticles,
                        aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kDS};
    candidate3ProngMcLoop<true, true, false, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDplusDs, "Process MC for D+ and Ds+ signals", false);

  void processMcDplusDsLc(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelDsToKKPi, aod::HfSelLc> const& candidates,
                          TracksWithSelectionMC const& tracks,
                          aod::McParticles const& mcParticles,
                          aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kDS, Pdg::kLambdaCPlus};
    candidate3ProngMcLoop<true, true, true, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDplusDsLc, "Process MC for D+, Ds+, and Lc signals", false);

  void processMcDplusLc(soa::Join<aod::HfCand3Prong, aod::HfSelDplusToPiKPi, aod::HfSelLc> const& candidates,
                        TracksWithSelectionMC const& tracks,
                        aod::McParticles const& mcParticles,
                        aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDPlus, Pdg::kLambdaCPlus};
    candidate3ProngMcLoop<true, false, true, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDplusLc, "Process MC for D+ and Lc signals", false);

  void processMcDsLc(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi, aod::HfSelLc> const& candidates,
                     TracksWithSelectionMC const& tracks,
                     aod::McParticles const& mcParticles,
                     aod::McCollisionLabels const& colls)
  {
    std::vector<int> const pdgCodes{Pdg::kDS, Pdg::kLambdaCPlus};
    candidate3ProngMcLoop<false, true, true, false>(candidates, tracks, mcParticles, colls, pdgCodes);
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processMcDsLc, "Process MC for Ds+ and Lc signals", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcEfficiency>(cfgc)};
}
