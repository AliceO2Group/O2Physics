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

/// \file treeCreatorDssToKKPi.cxx
/// \brief Writer of the 3 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables.
/// \note Extended from treeCreatorBplusToD0Pi.cxx
///
/// \author Stefano Politanò <stefano.politano@polito.it>, Politecnico & INFN, Torino

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::aod::hf_cand_3prong;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
// Candidate properties
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(DeltaMassPhi, deltaMassPhi, float);
DECLARE_SOA_COLUMN(AbsCos3PiK, absCos3PiK, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);

DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t); // is prompt or non-prompt, reco level
DECLARE_SOA_COLUMN(OriginMcGen, originMcGen, int8_t); // is prompt or non-prompt, Gen level
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandDsFull, "AOD", "HFCANDDsFull",
                  full::RSecondaryVertex,
                  full::M,
                  full::Pt,
                  full::Y,
                  full::DecayLengthXY,
                  full::DecayLengthXYNormalised,
                  full::CPA,
                  full::CPAXY,
                  full::DeltaMassPhi,
                  full::AbsCos3PiK,
                  full::MCflag, // is signal or background
                  full::OriginMcRec); // is prompt or non-prompt, reco level

DECLARE_SOA_TABLE(HfCandDsFullParticles, "AOD", "HFCANDDsFullP",
                  collision::BCId,
                  full::Pt,
                  full::Y,
                  full::MCflag,
                  full::OriginMcGen); // is prompt or non-prompt, Gen level

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorDsToKKPi {
  Produces<o2::aod::HfCandDsFull> rowCandidateFull;
  Produces<o2::aod::HfCandDsFullParticles> rowCandidateFullParticles;

  Configurable<int> isSignal{"isSignal", 1, "save only MC matched candidates"}; // 0 = all, 1 = signal

  void init(InitContext const&)
  {
  }

  template <typename T>
  auto fillTable(const T& candidate, int candFlag, int selection,
                 double invMass, double yDs, int8_t flagMc, int8_t origin)
  {
    if (selection >= 1) {
      rowCandidateFull(
        candidate.rSecondaryVertex(),
        invMass,
        candidate.pt(),
        yDs,
        candidate.decayLengthXY(),
        candidate.decayLengthXYNormalised(),
        candidate.cpa(),
        candidate.cpaXY(),
        deltaMassPhiDsToKKPi(candidate),
        std::abs(cos3PiKDsToKKPi(candidate)),
        flagMc,
        origin);
    }
  }

  void processData(soa::Join<aod::HfCand3Prong, aod::HfSelDsToKKPi> const& candidates)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      auto yD = yDs(candidate);
      fillTable(candidate, 0, candidate.isSelDsToKKPi(), yD, invMassDsToKKPi(candidate), 0, 0);
      fillTable(candidate, 1, candidate.isSelDsToPiKK(), yD, invMassDsToPiKK(candidate), 0, 0);
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDsToKKPi, processData, "Process data", true);

  void processMc(aod::Collisions const& collisions,
                 aod::McCollisions const&,
                 soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelDsToKKPi> const& candidates,
                 soa::Join<aod::McParticles, aod::HfCand3ProngMcGen> const& particles,
                 aod::BigTracksPID const&)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto const& candidate : candidates) {
      if (std::abs(candidate.flagMcMatchRec()) >= isSignal) {
      auto yD = yDs(candidate);
      fillTable(candidate, 0, candidate.isSelDsToKKPi(), invMassDsToKKPi(candidate), yD, candidate.flagMcMatchRec(), candidate.originMcRec());
      fillTable(candidate, 1, candidate.isSelDsToPiKK(), invMassDsToPiKK(candidate), yD, candidate.flagMcMatchRec(), candidate.originMcRec());
      }
    }

    // Filling particle properties
    rowCandidateFullParticles.reserve(particles.size());
    for (auto const& particle : particles) {
      if (std::abs(particle.flagMcMatchGen()) == 1 << DecayType::DsToKKPi) {
        rowCandidateFullParticles(
          particle.mcCollision().bcId(),
          particle.pt(),
          RecoDecay::y(array{particle.px(), particle.py(), particle.pz()}, RecoDecay::getMassPDG(particle.pdgCode())),
          particle.flagMcMatchGen(),
          particle.originMcGen());
      }
    }
  }

  PROCESS_SWITCH(HfTreeCreatorDsToKKPi, processMc, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorDsToKKPi>(cfgc));
  return workflow;
}
