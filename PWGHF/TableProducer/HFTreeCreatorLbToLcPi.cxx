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

/// \file HFTreeCreatorLbToLcPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from HFTreeCreatorD0ToKPi, HFTreeCreatorLcToPKPi, HFTreeCreatorXToJpsiPiPi
///
/// \author Panos Christakoglou <Panos.Christakoglou@cern.ch>, Nikhef
/// \author Maurice Jongerhuis <m.v.jongerhuis@students.uu.nl>, University Utrecht

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "PWGHF/DataModel/HFSecondaryVertex.h"
#include "PWGHF/DataModel/HFCandidateSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/TrackSelectorPID.h"
#include "ALICE3/DataModel/RICH.h"
#include "Common/Core/PID/PIDResponse.h"
#include "ReconstructionDataFormats/PID.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod::hf_cand;
using namespace o2::aod::hf_cand_lb;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(RSecondaryVertex, rSecondaryVertex, float);
DECLARE_SOA_COLUMN(PtProng0, ptProng0, float);
DECLARE_SOA_COLUMN(PProng0, pProng0, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised0, impactParameterNormalised0, float);
DECLARE_SOA_COLUMN(PtProng1, ptProng1, float);
DECLARE_SOA_COLUMN(PProng1, pProng1, float);
DECLARE_SOA_COLUMN(ImpactParameterNormalised1, impactParameterNormalised1, float);
DECLARE_SOA_COLUMN(CandidateSelFlag, candidateSelFlag, int8_t);
DECLARE_SOA_COLUMN(M, m, float);
DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(P, p, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(DecayLength, decayLength, float);
DECLARE_SOA_COLUMN(DecayLengthXY, decayLengthXY, float);
DECLARE_SOA_COLUMN(DecayLengthNormalised, decayLengthNormalised, float);
DECLARE_SOA_COLUMN(DecayLengthXYNormalised, decayLengthXYNormalised, float);
DECLARE_SOA_COLUMN(CPA, cpa, float);
DECLARE_SOA_COLUMN(CPAXY, cpaXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(MCflag, mcflag, int8_t);
DECLARE_SOA_COLUMN(NSigRICHPi0, nsigRICHPi0, float);
DECLARE_SOA_COLUMN(NSigfRICHPi0, nsigfRICHPi0, float);
DECLARE_SOA_COLUMN(NSigTOFPi0, nsigTOFPi0, float);
// Lc selection parameters
DECLARE_SOA_COLUMN(LcM, lcM, float);
DECLARE_SOA_COLUMN(LcCt, lcCt, float);
DECLARE_SOA_COLUMN(LcY, lcY, float);
DECLARE_SOA_COLUMN(LcE, lcE, float);
DECLARE_SOA_COLUMN(LcEta, lcEta, float);
DECLARE_SOA_COLUMN(LcCPA, lcCPA, float);
DECLARE_SOA_COLUMN(LcCPAXY, lcCPAXY, float);
DECLARE_SOA_COLUMN(LcChi2PCA, lcChi2PCA, float);
DECLARE_SOA_COLUMN(LcDecayLength, lcDecayLength, float);
DECLARE_SOA_COLUMN(LcDecayLengthXY, lcDecayLengthXY, float);
DECLARE_SOA_COLUMN(LcDecayLengthNormalised, lcDecayLengthNormalised, float);
DECLARE_SOA_COLUMN(LcImpactParameter0, lcImpactParameter0, float);
DECLARE_SOA_COLUMN(LcImpactParameter1, lcImpactParameter1, float);
DECLARE_SOA_COLUMN(LcImpactParameter2, lcImpactParameter2, float);
DECLARE_SOA_COLUMN(NSigRICHTrk1Pi, nSigRICHTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigRICHTrk1Pr, nSigRICHTrk1Pr, float);
DECLARE_SOA_COLUMN(NSigRICHTrk2Ka, nSigRICHTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigRICHTrk3Pi, nSigRICHTrk3Pi, float);
DECLARE_SOA_COLUMN(NSigRICHTrk3Pr, nSigRICHTrk3Pr, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk1Pi, nSigfRICHTrk1Pi, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk1Pr, nSigfRICHTrk1Pr, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk2Ka, nSigfRICHTrk2Ka, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk3Pi, nSigfRICHTrk3Pi, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk3Pr, nSigfRICHTrk3Pr, float);
DECLARE_SOA_COLUMN(NSigTOFTrk1Pi, nSigTOFrk1Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk1Pr, nSigTOFrk1Pr, float);
DECLARE_SOA_COLUMN(NSigTOFTrk2Ka, nSigTOFrk2Ka, float);
DECLARE_SOA_COLUMN(NSigTOFTrk3Pi, nSigTOFrk3Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk3Pr, nSigTOFrk3Pr, float);
// Events
DECLARE_SOA_COLUMN(IsEventReject, isEventReject, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandLbFull, "AOD", "HFCANDLbFull",
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  hf_cand::Chi2PCA,
                  full::ImpactParameterNormalised0,
                  full::PtProng0,
                  full::PProng0,
                  full::ImpactParameterNormalised1,
                  full::PtProng1,
                  full::PProng1,
                  hf_cand::PxProng0,
                  hf_cand::PyProng0,
                  hf_cand::PzProng0,
                  hf_cand::PxProng1,
                  hf_cand::PyProng1,
                  hf_cand::PzProng1,
                  hf_cand::ImpactParameter0,
                  hf_cand::ImpactParameter1,
                  hf_cand::ErrorImpactParameter0,
                  hf_cand::ErrorImpactParameter1,
                  full::NSigTOFPi0,
                  full::NSigRICHPi0,
                  full::NSigRICHTrk1Pi,
                  full::NSigRICHTrk1Pr,
                  full::NSigRICHTrk2Ka,
                  full::NSigRICHTrk3Pi,
                  full::NSigRICHTrk3Pr,
                  full::NSigfRICHPi0,
                  full::NSigfRICHTrk1Pi,
                  full::NSigfRICHTrk1Pr,
                  full::NSigfRICHTrk2Ka,
                  full::NSigfRICHTrk3Pi,
                  full::NSigfRICHTrk3Pr,
                  full::NSigTOFTrk1Pi,
                  full::NSigTOFTrk1Pr,
                  full::NSigTOFTrk2Ka,
                  full::NSigTOFTrk3Pi,
                  full::NSigTOFTrk3Pr,
                  full::LcM,
                  full::LcCt,
                  full::LcY,
                  full::LcE,
                  full::LcEta,
                  full::LcCPA,
                  full::LcCPAXY,
                  full::LcChi2PCA,
                  full::LcDecayLength,
                  full::LcDecayLengthXY,
                  full::LcDecayLengthNormalised,
                  full::LcImpactParameter0,
                  full::LcImpactParameter1,
                  full::LcImpactParameter2,
                  full::CandidateSelFlag,
                  full::M,
                  full::Pt,
                  full::P,
                  full::CPA,
                  full::CPAXY,
                  full::Ct,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

DECLARE_SOA_TABLE(HfCandLbFullEvents, "AOD", "HFCANDLbFullE",
                  collision::BCId,
                  collision::NumContrib,
                  collision::PosX,
                  collision::PosY,
                  collision::PosZ,
                  full::IsEventReject,
                  full::RunNumber);

DECLARE_SOA_TABLE(HfCandLbFullParticles, "AOD", "HFCANDLbFullP",
                  collision::BCId,
                  full::Pt,
                  full::Eta,
                  full::Phi,
                  full::Y,
                  full::MCflag);

} // namespace o2::aod

namespace o2::aod
{
namespace hf_track_index_alice3_pid
{
DECLARE_SOA_INDEX_COLUMN(Track, track); //!
DECLARE_SOA_INDEX_COLUMN(RICH, rich);   //!
DECLARE_SOA_INDEX_COLUMN(FRICH, frich); //!
} // namespace hf_track_index_alice3_pid

DECLARE_SOA_INDEX_TABLE_USER(HfTrackIndexALICE3PID, Tracks, "HFTRKIDXA3PID", //!
                             hf_track_index_alice3_pid::TrackId,
                             hf_track_index_alice3_pid::RICHId,
                             hf_track_index_alice3_pid::FRICHId);
} // namespace o2::aod

struct Alice3PidIndexBuilder {
  Builds<o2::aod::HfTrackIndexALICE3PID> index;
  void init(o2::framework::InitContext&) {}
};

/// Writes the full information in an output TTree
struct HfTreeCreatorLbToLcPi {
  Produces<o2::aod::HfCandLbFull> rowCandidateFull;
  Produces<o2::aod::HfCandLbFullEvents> rowCandidateFullEvents;
  Produces<o2::aod::HfCandLbFullParticles> rowCandidateFullParticles;

  void init(InitContext const&)
  {
  }

  using TracksPID = soa::Join<aod::BigTracksPID, aod::HfTrackIndexALICE3PID>;
  using ExtendedTracksPID = soa::Join<TracksPID, aod::TracksExtended>;

  void process(aod::Collisions const& collisions,
               aod::McCollisions const& mccollisions,
               soa::Join<aod::HfCandLb, aod::HfCandLbMCRec, aod::HFSelLbToLcPiCandidate> const& candidates,
               soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec, aod::HFSelLcCandidate> const& Lccandidates,
               soa::Join<aod::McParticles, aod::HfCandLbMCGen> const& particles,
               aod::BigTracksMC const& bigtracksmc,
               ExtendedTracksPID const&,
               aod::FRICHs const&,
               aod::RICHs const&)
  {

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (auto& candidate : candidates) {
      auto fillTable = [&](int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) {
          auto LcCand = candidate.index0_as<soa::Join<aod::HfCandProng3, aod::HfCandProng3MCRec, aod::HFSelLcCandidate>>();
          auto track0 = candidate.index1_as<ExtendedTracksPID>(); // daughter pion track
          auto track1 = LcCand.index0_as<ExtendedTracksPID>();    // granddaughter tracks (lc decay particles)
          auto track2 = LcCand.index1_as<ExtendedTracksPID>();
          auto track3 = LcCand.index2_as<ExtendedTracksPID>();

          auto RICHPi0 = -5000.0;
          auto RICHTrk1Pi = -5000.0;
          auto RICHTrk1p = -5000.0;
          auto RICHTrk2K = -5000.0;
          auto RICHTrk3Pi = -5000.0;
          auto RICHTrk3p = -5000.0;

          auto fRICHPi0 = -5000.0;
          auto fRICHTrk1Pi = -5000.0;
          auto fRICHTrk1p = -5000.0;
          auto fRICHTrk2K = -5000.0;
          auto fRICHTrk3Pi = -5000.0;
          auto fRICHTrk3p = -5000.0;

          if (track0.has_rich())
            RICHPi0 = track0.rich().richNsigmaPi();
          if (track1.has_rich())
            RICHTrk1Pi = track1.rich().richNsigmaPi();
          if (track1.has_rich())
            RICHTrk1p = track1.rich().richNsigmaPr();
          if (track2.has_rich())
            RICHTrk2K = track2.rich().richNsigmaKa();
          if (track3.has_rich())
            RICHTrk3Pi = track3.rich().richNsigmaPi();
          if (track3.has_rich())
            RICHTrk3p = track3.rich().richNsigmaPr();

          if (track0.has_frich())
            fRICHPi0 = track0.frich().frichNsigmaPi();
          if (track1.has_frich())
            fRICHTrk1Pi = track1.frich().frichNsigmaPi();
          if (track1.has_frich())
            fRICHTrk1p = track1.frich().frichNsigmaPr();
          if (track2.has_frich())
            fRICHTrk2K = track2.frich().frichNsigmaKa();
          if (track3.has_frich())
            fRICHTrk3Pi = track3.frich().frichNsigmaPi();
          if (track3.has_frich())
            fRICHTrk3p = track3.frich().frichNsigmaPr();

          rowCandidateFull(
            candidate.rSecondaryVertex(),
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.decayLengthNormalised(),
            candidate.decayLengthXYNormalised(),
            candidate.chi2PCA(),
            candidate.impactParameterNormalised0(),
            candidate.ptProng0(),
            RecoDecay::P(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.impactParameterNormalised1(),
            candidate.ptProng1(),
            RecoDecay::P(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
            candidate.pxProng0(),
            candidate.pyProng0(),
            candidate.pzProng0(),
            candidate.pxProng1(),
            candidate.pyProng1(),
            candidate.pzProng1(),
            candidate.impactParameter0(),
            candidate.impactParameter1(),
            candidate.errorImpactParameter0(),
            candidate.errorImpactParameter1(),
            track0.tofNSigmaPi(),
            RICHPi0,
            RICHTrk1Pi,
            RICHTrk1p,
            RICHTrk2K,
            RICHTrk3Pi,
            RICHTrk3p,
            fRICHPi0,
            fRICHTrk1Pi,
            fRICHTrk1p,
            fRICHTrk2K,
            fRICHTrk3Pi,
            fRICHTrk3p,
            track1.tofNSigmaPi(),
            track1.tofNSigmaPr(),
            track2.tofNSigmaKa(),
            track3.tofNSigmaPi(),
            track3.tofNSigmaPr(),
            o2::aod::hf_cand_prong3::InvMassLcpKpi(LcCand),
            o2::aod::hf_cand_prong3::CtLc(LcCand),
            o2::aod::hf_cand_prong3::YLc(LcCand),
            o2::aod::hf_cand_prong3::ELc(LcCand),
            LcCand.eta(),
            LcCand.cpa(),
            LcCand.cpaXY(),
            LcCand.chi2PCA(),
            LcCand.decayLength(),
            LcCand.decayLengthXY(),
            LcCand.decayLengthXYNormalised(),
            LcCand.impactParameter0(),
            LcCand.impactParameter1(),
            LcCand.impactParameter2(),
            FunctionSelection,
            FunctionInvMass,
            candidate.pt(),
            candidate.p(),
            candidate.cpa(),
            candidate.cpaXY(),
            FunctionCt,
            candidate.eta(),
            candidate.phi(),
            FunctionY,
            candidate.flagMCMatchRec());
        }
      };
      fillTable(candidate.isSelLbToLcPi(), InvMassLbToLcPi(candidate), CtLb(candidate), YLb(candidate));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<Alice3PidIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLbToLcPi>(cfgc));
  return workflow;
}