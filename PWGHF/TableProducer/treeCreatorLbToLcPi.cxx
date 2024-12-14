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

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::framework;

namespace o2::aod
{
namespace full
{
DECLARE_SOA_COLUMN(ImpactParameterXY, impactParameterXY, float);
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
DECLARE_SOA_COLUMN(CPA, cPA, float);
DECLARE_SOA_COLUMN(CPAXY, cPAXY, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(McFlag, mcFlag, int8_t);
DECLARE_SOA_COLUMN(OriginMcRec, originMcRec, int8_t);
DECLARE_SOA_COLUMN(SignTrk0, signTrk0, int16_t);
DECLARE_SOA_COLUMN(NSigTOFTrk0Pi, nSigTOFTrk0Pi, float);
DECLARE_SOA_COLUMN(NSigTPCTrk0Pi, nSigTPCTrk0Pi, float);
// Lc selection parameters
DECLARE_SOA_COLUMN(LcM, lcM, float);
DECLARE_SOA_COLUMN(LcCt, lcCt, float);
DECLARE_SOA_COLUMN(LcY, lcY, float);
DECLARE_SOA_COLUMN(LcE, lcE, float);
DECLARE_SOA_COLUMN(LcEta, lcEta, float);
DECLARE_SOA_COLUMN(LcVertexX, lcVertexX, float);
DECLARE_SOA_COLUMN(LcVertexY, lcVertexY, float);
DECLARE_SOA_COLUMN(LcVertexZ, lcVertexZ, float);
DECLARE_SOA_COLUMN(LcCPA, lcCPA, float);
DECLARE_SOA_COLUMN(LcCPAXY, lcCPAXY, float);
DECLARE_SOA_COLUMN(LcChi2PCA, lcChi2PCA, float);
DECLARE_SOA_COLUMN(LcDecayLength, lcDecayLength, float);
DECLARE_SOA_COLUMN(LcDecayLengthXY, lcDecayLengthXY, float);
DECLARE_SOA_COLUMN(LcDecayLengthNormalised, lcDecayLengthNormalised, float);
DECLARE_SOA_COLUMN(LcImpactParameter0, lcImpactParameter0, float);
DECLARE_SOA_COLUMN(LcImpactParameterError0, lcImpactParameterError0, float);
DECLARE_SOA_COLUMN(LcImpactParameter1, lcImpactParameter1, float);
DECLARE_SOA_COLUMN(LcImpactParameterError1, lcImpactParameterError1, float);
DECLARE_SOA_COLUMN(LcImpactParameter2, lcImpactParameter2, float);
DECLARE_SOA_COLUMN(LcImpactParameterError2, lcImpactParameterError2, float);
DECLARE_SOA_COLUMN(LcPx0, lcPx0, float);
DECLARE_SOA_COLUMN(LcPy0, lcPy0, float);
DECLARE_SOA_COLUMN(LcPz0, lcPz0, float);
DECLARE_SOA_COLUMN(LcPx1, lcPx1, float);
DECLARE_SOA_COLUMN(LcPy1, lcPy1, float);
DECLARE_SOA_COLUMN(LcPz1, lcPz1, float);
DECLARE_SOA_COLUMN(LcPx2, lcPx2, float);
DECLARE_SOA_COLUMN(LcPy2, lcPy2, float);
DECLARE_SOA_COLUMN(LcPz2, lcPz2, float);
DECLARE_SOA_COLUMN(LcSignProng0, lcSignProng0, int16_t);
DECLARE_SOA_COLUMN(LcSignProng1, lcSignProng1, int16_t);
DECLARE_SOA_COLUMN(LcSignProng2, lcSignProng2, int16_t);
DECLARE_SOA_COLUMN(LcNSigTPCPi0, lcNSigTPCPi0, float);
DECLARE_SOA_COLUMN(LcNSigTPCK0, lcNSigTPCK0, float);
DECLARE_SOA_COLUMN(LcNSigTPCPr0, lcNSigTPCPr0, float);
DECLARE_SOA_COLUMN(LcNSigTPCPi1, lcNSigTPCPi1, float);
DECLARE_SOA_COLUMN(LcNSigTPCK1, lcNSigTPCK1, float);
DECLARE_SOA_COLUMN(LcNSigTPCPr1, lcNSigTPCPr1, float);
DECLARE_SOA_COLUMN(LcNSigTPCPi2, lcNSigTPCPi2, float);
DECLARE_SOA_COLUMN(LcNSigTPCK2, lcNSigTPCK2, float);
DECLARE_SOA_COLUMN(LcNSigTPCPr2, lcNSigTPCPr2, float);
DECLARE_SOA_COLUMN(LcNSigTOFPr0, lcNSigTOFPr0, float);
DECLARE_SOA_COLUMN(LcNSigTOFK1, lcNSigTOFK1, float);
DECLARE_SOA_COLUMN(LcNSigTOFPi2, lcNSigTOFPi2, float);
} // namespace full

// put the arguments into the table
DECLARE_SOA_TABLE(HfCandLbFulls, "AOD", "HFCANDLBFULL",
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand::XSecondaryVertex, hf_cand::YSecondaryVertex, hf_cand::ZSecondaryVertex,
                  hf_cand::ErrorDecayLength, hf_cand::ErrorDecayLengthXY,
                  hf_cand::Chi2PCA,
                  full::ImpactParameterXY,
                  full::RSecondaryVertex,
                  full::DecayLength,
                  full::DecayLengthXY,
                  full::DecayLengthNormalised,
                  full::DecayLengthXYNormalised,
                  full::ImpactParameterNormalised0,
                  full::PtProng0,
                  full::PProng0,
                  full::ImpactParameterNormalised1,
                  full::PtProng1,
                  full::PProng1,
                  full::SignTrk0,
                  full::NSigTOFTrk0Pi,
                  full::NSigTPCTrk0Pi,
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
                  full::LcNSigTPCPi0, full::LcNSigTPCK0, full::LcNSigTPCPr0,
                  full::LcNSigTPCPi1, full::LcNSigTPCK1, full::LcNSigTPCPr1,
                  full::LcNSigTPCPi2, full::LcNSigTPCK2, full::LcNSigTPCPr2,
                  full::LcNSigTOFPr0,
                  full::LcNSigTOFK1,
                  full::LcNSigTOFPi2,
                  full::LcM,
                  full::LcCt,
                  full::LcY,
                  full::LcE,
                  full::LcEta,
                  full::LcVertexX,
                  full::LcVertexY,
                  full::LcVertexZ,
                  full::LcCPA,
                  full::LcCPAXY,
                  full::LcChi2PCA,
                  full::LcDecayLength,
                  full::LcDecayLengthXY,
                  full::LcDecayLengthNormalised,
                  full::LcImpactParameter0,
                  full::LcImpactParameterError0,
                  full::LcImpactParameter1,
                  full::LcImpactParameterError1,
                  full::LcImpactParameter2,
                  full::LcImpactParameterError2,
                  full::LcPx0, full::LcPy0, full::LcPz0,
                  full::LcPx1, full::LcPy1, full::LcPz1,
                  full::LcPx2, full::LcPy2, full::LcPz2,
                  full::LcSignProng0,
                  full::LcSignProng1,
                  full::LcSignProng2,
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
                  full::McFlag,
                  full::OriginMcRec);

} // namespace o2::aod

/// Writes the full information in an output TTree
struct HfTreeCreatorLbToLcPi {
  Produces<o2::aod::HfCandLbFulls> rowCandidateFull;
  HfHelper hfHelper;

  using TracksWPid = soa::Join<aod::Tracks, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

  void process(soa::Join<aod::HfCandLb, aod::HfSelLbToLcPi> const& candidates,
               soa::Join<aod::HfCand3Prong, aod::HfSelLc> const&,
               TracksWPid const&)
  {
    // Filling candidate properties
    rowCandidateFull.reserve(candidates.size());
    for (const auto& candidate : candidates) {
      auto fillTable = [&](int FunctionSelection,
                           float FunctionInvMass,
                           float FunctionCt,
                           float FunctionY) {
        auto candLc = candidate.prong0_as<soa::Join<aod::HfCand3Prong, aod::HfSelLc>>();
        auto track0 = candidate.prong1_as<TracksWPid>(); // daughter pion track
        auto track1 = candLc.prong0_as<TracksWPid>();    // granddaughter tracks (lc decay particles)
        auto track2 = candLc.prong1_as<TracksWPid>();
        auto track3 = candLc.prong2_as<TracksWPid>();

        auto tempConst = -1; // For data

        rowCandidateFull(
          candidate.posX(),
          candidate.posY(),
          candidate.posZ(),
          candidate.xSecondaryVertex(),
          candidate.ySecondaryVertex(),
          candidate.zSecondaryVertex(),
          candidate.errorDecayLength(),
          candidate.errorDecayLengthXY(),
          candidate.chi2PCA(),
          candidate.impactParameterXY(),
          candidate.rSecondaryVertex(),
          candidate.decayLength(),
          candidate.decayLengthXY(),
          candidate.decayLengthNormalised(),
          candidate.decayLengthXYNormalised(),
          candidate.impactParameterNormalised0(),
          candidate.ptProng0(),
          RecoDecay::p(candidate.pxProng0(), candidate.pyProng0(), candidate.pzProng0()),
          candidate.impactParameterNormalised1(),
          candidate.ptProng1(),
          RecoDecay::p(candidate.pxProng1(), candidate.pyProng1(), candidate.pzProng1()),
          track0.sign(),
          track0.tofNSigmaPi(),
          track0.tpcNSigmaPi(),
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
          track1.tpcNSigmaPi(),
          track1.tpcNSigmaKa(),
          track1.tpcNSigmaPr(),
          track2.tpcNSigmaPi(),
          track2.tpcNSigmaKa(),
          track2.tpcNSigmaPr(),
          track3.tpcNSigmaPi(),
          track3.tpcNSigmaKa(),
          track3.tpcNSigmaPr(),
          track1.tofNSigmaPr(),
          track2.tofNSigmaKa(),
          track3.tofNSigmaPi(),
          hfHelper.invMassLcToPKPi(candLc),
          hfHelper.ctLc(candLc),
          hfHelper.yLc(candLc),
          hfHelper.eLc(candLc),
          candLc.eta(),
          candLc.xSecondaryVertex(),
          candLc.ySecondaryVertex(),
          candLc.zSecondaryVertex(),
          candLc.cpa(),
          candLc.cpaXY(),
          candLc.chi2PCA(),
          candLc.decayLength(),
          candLc.decayLengthXY(),
          candLc.decayLengthXYNormalised(),
          candLc.impactParameter0(),
          candLc.errorImpactParameter0(),
          candLc.impactParameter1(),
          candLc.errorImpactParameter1(),
          candLc.impactParameter2(),
          candLc.errorImpactParameter2(),
          track1.px(), track1.py(), track1.pz(),
          track2.px(), track2.py(), track2.pz(),
          track3.px(), track3.py(), track3.pz(),
          track1.sign(), track2.sign(), track3.sign(),
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
          tempConst,
          tempConst);
      };
      fillTable(candidate.isSelLbToLcPi(), hfHelper.invMassLbToLcPi(candidate), hfHelper.ctLb(candidate), hfHelper.yLb(candidate));
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow;
  workflow.push_back(adaptAnalysisTask<HfTreeCreatorLbToLcPi>(cfgc));
  return workflow;
}
