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

/// \file taskMcEfficiencyToXiPi.cxx
/// \brief Task to analyse the MC efficiency for toXiPi charm baryon decay ("found / trackable")
///
/// \author Federica Zanone, Heidelberg University

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskMcEfficiencyToXiPi {
  Configurable<float> rapidityCharmBaryonMax{"rapidityCharmBaryonMax", 0.5, "Max absolute value of rapidity for charm baryon"};
  Configurable<float> acceptanceEtaLf{"acceptanceEtaLf", 1.0, "Max absolute value of eta for LF daughters"};
  Configurable<float> acceptancePtPionFromCascade{"acceptancePtPionFromCascade", 0.2, "Min value of pt for pion <-- cascade"};
  Configurable<float> acceptanceEtaPionFromCharm{"acceptanceEtaPionFromCharm", 0.8, "Max absolute value of eta for pion <-- charm baryon"};
  Configurable<float> acceptancePtPionFromCharm{"acceptancePtPionFromCharm", 0.5, "Min value of pt for pion <-- charm baryon"};

  Configurable<int> nClustersTpcMin{"nClustersTpcMin", 70, "Minimum number of TPC clusters requirement for pion <-- charm baryon"};
  Configurable<int> nClustersItsMin{"nClustersItsMin", 3, "Minimum number of ITS clusters requirement for pion <- charm baryon"};

  Configurable<bool> matchOmegac{"matchOmegac", false, "Do MC studies for Omegac0"};
  Configurable<bool> matchXic{"matchXic", true, "Do MC studies for Xic0"};

  ConfigurableAxis axisPt{"axisPt", {200, 0, 20}, "pT axis"};
  ConfigurableAxis axisMass{"axisMass", {900, 2.1, 3}, "m_inv axis"};

  enum HFStep { kHFStepMC = 0,              // all generated mothers
                kHFStepMcInRapidity,        // MC mother in rapidity |y| < rapidityCharmBaryonMax=0.5
                kHFStepAcceptance,          // MC mother where all final state candidates pass eta and pt selection
                kHFStepTrackable,           // MC mother where all final state candidates have a reconstructed track
                kHFStepAcceptanceTrackable, // MC mother where all final state candidates have a reconstructed track and pass eta and pt selection
                kHFStepItsTpcTrackableCuts, // MC mother where charm bachelor passes isGlobalTrkWoDca selection
                kHFStepItsTrackableCuts,    // MC mother where charm bachelor passes isITSQualityTrack selection
                kHFStepTracked,             // signal candidates which have been found by candidateCreator for the desired decay channel
                kHFStepTrackedCuts,         // signal candidates which have been found and pass the candidateSelector cuts but PID
                kHFStepTrackedSelected,     // signal candidates which pass the selector and pass the candidateSelector cuts including PID
                kHFNSteps };

  enum TrackableStep { kTrackableAll = 0, // all tracks
                       kTrackableITS,     // track contains ITS information
                       kTrackableTPC,     // track contains TPC information
                       kTrackableITSTPC,  // track contains ITS and TPC information
                       kNTrackableSteps };

  using TracksWithSelectionMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::McTrackLabels, aod::TrackSelection>;
  using CandidateInfo = soa::Join<aod::HfCandToXiPi, aod::HfToXiPiMCRec, aod::HfSelToXiPi>;
  using ParticleInfo = soa::Join<aod::McParticles, aod::HfToXiPiMCGen>;

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    auto hCandidates = registry.add<StepTHn>("hCandidates", "Candidate count at different steps", {HistType::kStepTHnF, {axisPt, axisMass, {2, -0.5, 1.5, "collision matched"}, {RecoDecay::OriginType::NonPrompt + 1, RecoDecay::OriginType::None - 0.5, RecoDecay::OriginType::NonPrompt + 0.5}}, kHFNSteps});
    hCandidates->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hCandidates->GetAxis(1)->SetTitle("#it{m}_{inv} (GeV/#it{c}^{2})");
    hCandidates->GetAxis(2)->SetTitle("Collision matched");
    hCandidates->GetAxis(3)->SetTitle("Charm hadron origin");

    auto hTrackablePtEta = registry.add<StepTHn>("hTrackablePtEta", "Prongs kinematics at different steps", {HistType::kStepTHnF, {{200, 0, 20}, {400, -2, 2}}, kNTrackableSteps});
    hTrackablePtEta->GetAxis(0)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    hTrackablePtEta->GetAxis(1)->SetTitle("#eta");
  }

  template <typename T>
  inline bool checkTrackGlbTrk(T const& track)
  {
    return (track.isGlobalTrackWoDCA() && track.tpcNClsFound() > nClustersTpcMin && track.itsNCls() > nClustersItsMin);
  }

  template <typename T>
  inline bool checkTrackItsTrk(T const& track)
  {
    return (track.isQualityTrackITS() && track.itsNCls() > nClustersItsMin);
  }

  // candidates -> join candidateCreator, candidateCreator McRec and candidateSelector tables
  template <typename T1>
  void candidateRecLoop(T1 const& candidates, int pdgCode)
  {
    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));

    bool selectedKine = false;
    bool selectedPid = false;

    double mass = 0.0;
    double pt = 0.0;

    int decayFlag = 0;
    if (pdgCode == Pdg::kXiC0) {
      decayFlag = 1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi;
    } else if (pdgCode == Pdg::kOmegaC0) {
      decayFlag = 1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi;
    } else {
      LOGP(fatal, "Not implemented for PDG code: ", pdgCode);
    }

    for (const auto& candidate : candidates) {

      // only take into account matched decays
      if (std::abs(candidate.flagMcMatchRec()) != decayFlag) {
        continue;
      }

      selectedKine = false;
      selectedPid = false;
      mass = candidate.invMassCharmBaryon();
      pt = RecoDecay::sqrtSumOfSquares(candidate.pxCharmBaryon(), candidate.pyCharmBaryon());

      // all candidates (candidateCreator)
      hCandidates->Fill(kHFStepTracked, pt, mass, candidate.collisionMatched(), candidate.originRec());

      // check xi-pi candidate passed candidate selector cuts (except PID)
      if (candidate.resultSelections() && candidate.statusInvMassLambda() && candidate.statusInvMassCascade() && candidate.statusInvMassCharmBaryon()) {
        hCandidates->Fill(kHFStepTrackedCuts, pt, mass, candidate.collisionMatched(), candidate.originRec());
        selectedKine = true;
      }
      if (!selectedKine) {
        continue;
      }

      // check xi-pi candidate passed candidate selector cuts (PID included)
      if (candidate.statusPidCharmBaryon()) {
        hCandidates->Fill(kHFStepTrackedSelected, pt, mass, candidate.collisionMatched(), candidate.originRec());
        selectedPid = true;
      }
      if (!selectedPid) {
        continue;
      }
    }
  }

  // candidates -> join candidateCreator, candidateCreator McRec and candidateSelector tables
  // genParticles -> join aod::McParticles and candidateCreator McGen tables
  template <typename T1, typename T2>
  void candidateFullLoop(T1 const& candidates, T2 const& genParticles, TracksWithSelectionMC const& tracks, aod::McCollisionLabels const& colls, int pdgCode)
  {
    // fill hCandidates histogram
    candidateRecLoop(candidates, pdgCode);

    auto hCandidates = registry.get<StepTHn>(HIST("hCandidates"));
    auto hTrackablePtEta = registry.get<StepTHn>(HIST("hTrackablePtEta"));

    // lists for optimization
    std::vector<bool> tracked(genParticles.size(), false);
    std::vector<bool> hasITS(genParticles.size(), false);
    std::vector<bool> hasTPC(genParticles.size(), false);
    std::vector<bool> selectedIts(genParticles.size(), false);
    std::vector<bool> selectedItsTpc(genParticles.size(), false);
    for (const auto& track : tracks) {
      if (track.mcParticleId() >= 0) {
        tracked[track.mcParticleId()] = true;
        if (track.hasITS()) {
          hasITS[track.mcParticleId()] = true;
        }
        if (track.hasTPC()) {
          hasTPC[track.mcParticleId()] = true;
        }
        if (checkTrackGlbTrk(track)) {
          selectedItsTpc[track.mcParticleId()] = true;
        }
        if (checkTrackItsTrk(track)) {
          selectedIts[track.mcParticleId()] = true;
        }
      }
    }

    double mass = 0.;
    int decayFlagGen = 0;

    if (pdgCode == Pdg::kXiC0) {
      decayFlagGen = 1 << aod::hf_cand_toxipi::DecayType::XiczeroToXiPi;
      mass = MassXiC0;
    } else if (pdgCode == Pdg::kOmegaC0) {
      decayFlagGen = 1 << aod::hf_cand_toxipi::DecayType::OmegaczeroToXiPi;
      mass = MassOmegaC0;
    } else {
      LOGP(fatal, "Not implemented for PDG code: ", pdgCode);
    }

    for (const auto& mcParticle : genParticles) {

      // check if I am treating the desired charm baryon
      if (std::abs(mcParticle.pdgCode()) != pdgCode) {
        continue;
      }
      /// check if the charm baryon decays to the desired channel
      if (std::abs(mcParticle.flagMcMatchGen()) != decayFlagGen) {
        continue;
      }

      // all candidates
      hCandidates->Fill(kHFStepMC, mcParticle.pt(), mass, true, mcParticle.originGen()); // set matchedCollision to true by default at gen level

      // candidates with charm baryon within eta range
      if (std::abs(mcParticle.y()) < rapidityCharmBaryonMax) {
        hCandidates->Fill(kHFStepMcInRapidity, mcParticle.pt(), mass, true, mcParticle.originGen());
      }

      // exclude cases with undesired decays
      if (mcParticle.daughtersIds().size() != 2) {
        LOGP(fatal, "Invalid numbers of daughters for charm baryon {}: {}", mcParticle.globalIndex(), mcParticle.daughtersIds().size());
      }
      auto cascId = mcParticle.daughtersIds()[0];
      auto pionId = mcParticle.daughtersIds()[1];
      if (cascId < 0 || pionId < 0) {
        continue;
      }
      auto cascade = genParticles.rawIteratorAt(mcParticle.daughtersIds().front());
      auto pion = genParticles.rawIteratorAt(mcParticle.daughtersIds().back());
      if (std::abs(cascade.pdgCode()) == kPiPlus) { // check if std::abs(cascade.pdgCode()) is different wrt kXiMinus and equal to kPiPlus (daughters ID assignment swapped)
        std::swap(cascade, pion);
        std::swap(cascId, pionId);
      } else if (std::abs(cascade.pdgCode()) == kXiMinus && std::abs(pion.pdgCode()) == kPiPlus) {
        LOGP(debug, "Correct assignment of charm baryon daughters IDs - gen level");
      } else {
        LOGP(fatal, "Invalid charm baryon daughters PDG codes");
      }

      // check pion <-- charm baryon pt and eta
      bool inAcceptance = true;
      if (std::abs(pion.eta()) > acceptanceEtaPionFromCharm || pion.pt() < acceptancePtPionFromCharm) {
        inAcceptance = false;
      }

      // check LF daughters pt (pion<--cascade) and eta
      // first create cascade daughters objects
      if (cascade.daughtersIds().size() != 2) {
        LOGP(fatal, "Invalid numbers of daughters for cascade {}: {}", cascade.globalIndex(), cascade.daughtersIds().size());
      }
      auto lambdaId = cascade.daughtersIds()[0];
      auto pionFromCascadeId = cascade.daughtersIds()[1];
      if (lambdaId < 0 || pionFromCascadeId < 0) {
        continue;
      }
      auto lambda = genParticles.rawIteratorAt(cascade.daughtersIds().front());
      auto pionFromCascade = genParticles.rawIteratorAt(cascade.daughtersIds().back());
      if (std::abs(lambda.pdgCode()) == kPiPlus) { // check if std::abs(lambda.pdgCode()) is different wrt kLambda0 and equal to kPiPlus (daughters ID assignment swapped)
        std::swap(lambda, pionFromCascade);
        std::swap(lambdaId, pionFromCascadeId);
      } else if (std::abs(lambda.pdgCode()) == kLambda0 && std::abs(pionFromCascade.pdgCode()) == kPiPlus) {
        LOGP(debug, "Correct assignment of cascade daughters IDs - gen level");
      } else {
        LOGP(fatal, "Invalid cascade daughters PDG codes");
      }
      // then create lambda daughters objects
      if (lambda.daughtersIds().size() != 2) {
        LOGP(fatal, "Invalid numbers of daughters for lambda {}: {}", lambda.globalIndex(), lambda.daughtersIds().size());
      }
      auto protonId = lambda.daughtersIds()[0];
      auto pionFromLambdaId = lambda.daughtersIds()[1];
      if (protonId < 0 || pionFromLambdaId < 0) {
        continue;
      }
      auto proton = genParticles.rawIteratorAt(lambda.daughtersIds().front());
      auto pionFromLambda = genParticles.rawIteratorAt(lambda.daughtersIds().back());
      if (std::abs(proton.pdgCode()) == kPiPlus) { // check if std::abs(proton.pdgCode()) is different wrt kProton and equal to kPiPlus (daughters ID assignment swapped)
        std::swap(proton, pionFromLambda);
        std::swap(protonId, pionFromLambdaId);
      } else if (std::abs(proton.pdgCode()) == kProton && std::abs(pionFromLambda.pdgCode()) == kPiPlus) {
        LOGP(debug, "Correct assignment of lambda daughters IDs - gen level");
      } else {
        LOGP(fatal, "Invalid lambda daughters PDG codes");
      }
      // check on pt and eta
      if (std::abs(pionFromCascade.eta()) > acceptanceEtaLf || std::abs(pionFromLambda.eta()) > acceptanceEtaLf || std::abs(proton.eta()) > acceptanceEtaLf) {
        inAcceptance = false;
      }
      if (pionFromCascade.pt() < acceptancePtPionFromCascade) {
        inAcceptance = false;
      }
      // final state candidates pass eta and pt selection
      if (inAcceptance) {
        hCandidates->Fill(kHFStepAcceptance, mcParticle.pt(), mass, true, mcParticle.originGen());
      }

      if (tracked[pionId] && tracked[pionFromCascadeId] && tracked[pionFromLambdaId] && tracked[protonId]) {
        // final state candidates have a mc particleID != 0
        hCandidates->Fill(kHFStepTrackable, mcParticle.pt(), mass, true, mcParticle.originGen());
        if (inAcceptance) {
          hCandidates->Fill(kHFStepAcceptanceTrackable, mcParticle.pt(), mass, true, mcParticle.originGen());
        } else {
          LOGP(debug, "Candidate {} not in acceptance but tracked.", mcParticle.globalIndex());
          LOGP(debug, "MC cascade: pt={} eta={}", cascade.pt(), cascade.eta());
          LOGP(debug, "MC pion<--charm baryon: pt={} eta={}", pion.pt(), pion.eta());
        }

        // final state daughters info
        // pion <- charm baryon
        hTrackablePtEta->Fill(kTrackableAll, pion.pt(), pion.eta());
        if (hasITS[pion.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITS, pion.pt(), pion.eta());
        }
        if (hasTPC[pion.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableTPC, pion.pt(), pion.eta());
        }
        if (hasITS[pion.globalIndex()] && hasTPC[pion.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITSTPC, pion.pt(), pion.eta());
        }
        // pion <- cascade
        hTrackablePtEta->Fill(kTrackableAll, pionFromCascade.pt(), pionFromCascade.eta());
        if (hasITS[pionFromCascade.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITS, pionFromCascade.pt(), pionFromCascade.eta());
        }
        if (hasTPC[pionFromCascade.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableTPC, pionFromCascade.pt(), pionFromCascade.eta());
        }
        if (hasITS[pionFromCascade.globalIndex()] && hasTPC[pionFromCascade.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITSTPC, pionFromCascade.pt(), pionFromCascade.eta());
        }
        // pion <- lambda
        hTrackablePtEta->Fill(kTrackableAll, pionFromLambda.pt(), pionFromLambda.eta());
        if (hasITS[pionFromLambda.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITS, pionFromLambda.pt(), pionFromLambda.eta());
        }
        if (hasTPC[pionFromLambda.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableTPC, pionFromLambda.pt(), pionFromLambda.eta());
        }
        if (hasITS[pionFromLambda.globalIndex()] && hasTPC[pionFromLambda.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITSTPC, pionFromLambda.pt(), pionFromLambda.eta());
        }
        // proton <- lambda
        hTrackablePtEta->Fill(kTrackableAll, proton.pt(), proton.eta());
        if (hasITS[proton.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITS, proton.pt(), proton.eta());
        }
        if (hasTPC[proton.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableTPC, proton.pt(), proton.eta());
        }
        if (hasITS[proton.globalIndex()] && hasTPC[proton.globalIndex()]) {
          hTrackablePtEta->Fill(kTrackableITSTPC, proton.pt(), proton.eta());
        }
      } // close tracked if

      // with pion track cuts (see checkTrackGlbTrk and checkTrkItsTrk)
      if (selectedIts[pionId]) {
        hCandidates->Fill(kHFStepItsTrackableCuts, mcParticle.pt(), mass, true, mcParticle.originGen());
        if (!inAcceptance) {
          LOGP(debug, "Candidate {} has daughters not in acceptance but pion <-- charm tracked and selected (its only)", mcParticle.globalIndex());
        }
      }
      if (selectedItsTpc[pionId]) {
        hCandidates->Fill(kHFStepItsTpcTrackableCuts, mcParticle.pt(), mass, true, mcParticle.originGen());
        if (!inAcceptance) {
          LOGP(debug, "Candidate {} has daughters not in acceptance but pion <-- charm tracked and selected (its & tpc)", mcParticle.globalIndex());
        }
      }

    } // close loop mcParticles
  }   // close candidateMcLoop

  // process functions
  void process(CandidateInfo const& candidates,
               ParticleInfo const& genParticles,
               TracksWithSelectionMC const& tracks,
               aod::McCollisionLabels const& colls)
  {
    if (matchXic && matchOmegac) {
      LOGP(fatal, "Can't match Omegac0 and Xic0 at the same time, please choose one");
    } else if (!matchXic && !matchOmegac) {
      LOGP(fatal, "Please match either Omegac0 or Xic0");
    } else if (matchXic) {
      candidateFullLoop(candidates, genParticles, tracks, colls, Pdg::kXiC0);
    } else if (matchOmegac) {
      candidateFullLoop(candidates, genParticles, tracks, colls, Pdg::kOmegaC0);
    }
  }

}; // close struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<HfTaskMcEfficiencyToXiPi>(cfgc)};
}
