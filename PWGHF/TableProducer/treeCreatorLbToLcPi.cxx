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

/// \file treeCreatorLbToLcPi.cxx
/// \brief Writer of the 2 prong candidates in the form of flat tables to be stored in TTrees.
///        Intended for debug or for the local optimization of analysis on small samples.
///        In this file are defined and filled the output tables
/// \note Extended from treeCreatorD0ToKPi.cxx, treeCreatorLcToPKPi.cxx, treeCreatorXToJpsiPiPi.cxx
///
/// \author Panos Christakoglou <Panos.Christakoglou@cern.ch>, Nikhef
/// \author Maurice Jongerhuis <m.v.jongerhuis@students.uu.nl>, University Utrecht

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "ALICE3/DataModel/RICH.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

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
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(NSigRICHTrk0Pi, nsigRICHTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigfRICHTrk0Pi, nsigfRICHTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigTOFTrk0Pi, nsigTOFTrk0Pi, float);
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
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandLbFulls, "AOD", "HFCANDLBFULL",
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
                  full::NSigTOFTrk0Pi,
                  full::NSigRICHTrk0Pi,
                  full::NSigRICHTrk1Pi,
                  full::NSigRICHTrk1Pr,
                  full::NSigRICHTrk2Ka,
                  full::NSigRICHTrk3Pi,
                  full::NSigRICHTrk3Pr,
                  full::NSigfRICHTrk0Pi,
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
                  full::MCflag,
                  full::OriginMcRec);

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

struct HfTreeCreatorLbToLcPiAlice3PidIndexBuilder {
  Builds<o2::aod::HfTrackIndexALICE3PID> index;

  void init(InitContext&) {}
};

/// Writes the full information in an output TTree
struct HfTreeCreatorLbToLcPi {
  Produces<o2::aod::HfCandLbFulls> rowCandidateFull;

  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::HfTrackIndexALICE3PID>;

  void process(soa::Join<aod::HfCandLb, aod::HfCandLbMcRec, aod::HfSelLbToLcPi> const& candidates,
               soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc> const&,
               TracksWPid const&,
               aod::FRICHs const&,
               aod::RICHs const&)
  {

    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      auto fillTable = [&](int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        if (FunctionSelection >= 1) {
          auto candLc = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfCand3ProngMcRec, aod::HfSelLc>>();
          auto track0 = candidate.prong1_as<TracksWPid>(); // daughter pion track
          auto track1 = candLc.prong0_as<TracksWPid>();    // granddaughter tracks (lc decay particles)
          auto track2 = candLc.prong1_as<TracksWPid>();
          auto track3 = candLc.prong2_as<TracksWPid>();

          auto RICHTrk0Pi = -5000.0;
          auto RICHTrk1Pi = -5000.0;
          auto RICHTrk1P = -5000.0;
          auto RICHTrk2K = -5000.0;
          auto RICHTrk3Pi = -5000.0;
          auto RICHTrk3P = -5000.0;

          auto fRICHTrk0Pi = -5000.0;
          auto fRICHTrk1Pi = -5000.0;
          auto fRICHTrk1P = -5000.0;
          auto fRICHTrk2K = -5000.0;
          auto fRICHTrk3Pi = -5000.0;
          auto fRICHTrk3P = -5000.0;

          if (track0.has_rich())
            RICHTrk0Pi = track0.rich().richNsigmaPi();
          if (track1.has_rich()) {
            RICHTrk1Pi = track1.rich().richNsigmaPi();
            RICHTrk1P = track1.rich().richNsigmaPr();
          }
          if (track2.has_rich())
            RICHTrk2K = track2.rich().richNsigmaKa();
          if (track3.has_rich()) {
            RICHTrk3Pi = track3.rich().richNsigmaPi();
            RICHTrk3P = track3.rich().richNsigmaPr();
          }

          if (track0.has_frich())
            fRICHTrk0Pi = track0.frich().frichNsigmaPi();
          if (track1.has_frich()) {
            fRICHTrk1Pi = track1.frich().frichNsigmaPi();
            fRICHTrk1P = track1.frich().frichNsigmaPr();
          }
          if (track2.has_frich())
            fRICHTrk2K = track2.frich().frichNsigmaKa();
          if (track3.has_frich()) {
            fRICHTrk3Pi = track3.frich().frichNsigmaPi();
            fRICHTrk3P = track3.frich().frichNsigmaPr();
          }

          rowCandidateFull(
            candidate.rSecondaryVertex(),
            candidate.decayLength(),
            candidate.decayLengthXY(),
            candidate.decayLengthNormalised(),
            candidate.decayLengthXYNormalised(),
            candidate.chi2PCA(),
            candidate.impactParameterNormalised0(),
            candidate.ptProng0(),
            RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
            candidate.impactParameterNormalised1(),
            candidate.ptProng1(),
            RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
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
            RICHTrk0Pi,
            RICHTrk1Pi,
            RICHTrk1P,
            RICHTrk2K,
            RICHTrk3Pi,
            RICHTrk3P,
            fRICHTrk0Pi,
            fRICHTrk1Pi,
            fRICHTrk1P,
            fRICHTrk2K,
            fRICHTrk3Pi,
            fRICHTrk3P,
            track1.tofNSigmaPi(),
            track1.tofNSigmaPr(),
            track2.tofNSigmaKa(),
            track3.tofNSigmaPi(),
            track3.tofNSigmaPr(),
            hfHelper.invMassLcToPKPi(candLc),
            hfHelper.ctLc(candLc),
            hfHelper.yLc(candLc),
            hfHelper.eLc(candLc),
            candLc.eta(),
            candLc.cpa(),
            candLc.cpaXY(),
            candLc.chi2PCA(),
            candLc.decayLength(),
            candLc.decayLengthXY(),
            candLc.decayLengthXYNormalised(),
            candLc.impactParameter0(),
            candLc.impactParameter1(),
            candLc.impactParameter2(),
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
            candidate.flagMcMatchRec(),
            candidate.originMcRec());
        }
      };
      fillTable(candidate.isSelLbToLcPi(), hfHelper.invMassLbToLcPi(candidate), hfHelper.ctLb(candidate), hfHelper.yLb(candidate));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLbToLcPiAlice3PidIndexBuilder>(cfgc));
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLbToLcPi>(cfgc));
  return workflow;
}
