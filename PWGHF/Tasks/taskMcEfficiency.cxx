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
using namespace o2::aod::hf_cand_2prong;
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
  void candidateLoop(T1& candidates, T2& tracks, T3& mcParticles)
  {
    using TracksType = std::decay_t<decltype(tracks)>;

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    std::map<int64_t, int> duplicates;

    for (const auto pdgCode : pdgCodes.value) {
      auto decayType = -1;
      std::array<int, 2> pdgDaughters;

      if (pdgCode == pdg::kD0) {
        decayType = 1 << DecayType::D0ToPiK;
        pdgDaughters[0] = +kPiPlus;
        pdgDaughters[1] = -kKPlus;
      } else if (pdgCode == pdg::kD0Bar) {
        decayType = 1 << DecayType::D0ToPiK;
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

  void processData(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelection& tracks)
  {
    candidateLoop<false>(candidates, tracks, tracks); // NOTE third argument has to be provided but is not used as template argument is <false>
  }
  PROCESS_SWITCH(HfTaskMcEfficiency, processData, "Process data (no MC information needed)", false);

  void processMc(soa::Join<aod::HfCand2Prong, aod::HfSelD0>& candidates, TracksWithSelectionMC& tracks, aod::McParticles& mcParticles, aod::McCollisionLabels&)
  {
    candidateLoop<true>(candidates, tracks, mcParticles);

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
        hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, pdgCode, 1.0, true);

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
  PROCESS_SWITCH(HfTaskMcEfficiency, processMc, "Process MC", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcEfficiency>(cfgc)};
}
